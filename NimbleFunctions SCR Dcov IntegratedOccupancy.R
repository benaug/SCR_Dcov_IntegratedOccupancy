dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0), InSS = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(InSS==1){
      logProb <- log(pi.cell)
    }else{
      logProb <- -Inf
    }
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0), InSS = double(0)) {
    returnType(double(0))
    return(0)
  }
)

Getpd2.j <- nimbleFunction(
  run = function(pd2 = double(2), z = double(1)){
    returnType(double(1))
    M <- nimDim(pd2)[1]
    J <- nimDim(pd2)[2]
    tmp <- rep(1,J)
    for(i in 1:M){
      if(z[i]==1){
        tmp <- tmp*(1-pd2[i,])
      }
    }
    pd2.j <- 1 - tmp
    return(pd2.j)
  }
)

GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dBinomialVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){#skip calculation if z=0
      return(0)
    }else{
      logProb <- sum(dbinom(x, size = K, p = pd, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0), pd = double(1), K = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#Required custom update for number of individuals
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(c("N","y1","y2","pd1","pd2","pd2.j")) #nodes to copy back to mvSaved
    inds.detected <- control$inds.detected
    z.ups <- control$z.ups
    J2 <- control$J2
    M <- control$M
    #nodes used for update, calcNodes + z nodes
    y1.nodes <- model$expandNodeNames("y1")
    y2.nodes <- model$expandNodeNames("y2")
    pd1.nodes <- model$expandNodeNames("pd1")
    pd2.nodes <- model$expandNodeNames("pd2")
    pd2.j.nodes <- model$expandNodeNames("pd2.j")
    N.node <- model$expandNodeNames("N")
    z.nodes <- model$expandNodeNames("z")
  },
  run = function() {
    
    #track these "manually" so computations faster than nimble will do them
    nondetect.probs.initial <- 1 - model$pd2.j #p(no detect)
    
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        if(any(pick==inds.detected)){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y1 <- model$getLogProb(y1.nodes[pick])
          lp.initial.y2 <- model$getLogProb(y2.nodes)

          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0

          #turn pd off
          #don't use calculate for pd2.j, compute and insert manually
          model$calculate(pd1.nodes[pick])
          nondetect.probs.proposed <- nondetect.probs.initial/(1-model$pd2[pick,]) #divide these out before calculate, which sets to 0
          model$calculate(pd2.nodes[pick])
          # model$calculate(pd2.j.nodes)
          model$pd2.j <<- 1 - nondetect.probs.proposed
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y1 <- model$calculate(y1.nodes[pick]) #will always be 0
          lp.proposed.y2 <- model$calculate(y2.nodes)

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y1 + lp.proposed.y2) - (lp.initial.N + lp.initial.y1 + lp.initial.y2)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd1",1][pick,] <<- model[["pd1"]][pick,]
            mvSaved["pd2",1][pick,] <<- model[["pd2"]][pick,]
            mvSaved["pd2.j",1][1:J2] <<- model[["pd2.j"]][1:J2]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            nondetect.probs.initial  <- nondetect.probs.proposed
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd1"]][pick,] <<- mvSaved["pd1",1][pick,]
            model[["pd2"]][pick,] <<- mvSaved["pd2",1][pick,]
            model[["pd2.j"]][1:J2] <<- mvSaved["pd2.j",1][1:J2]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y1.nodes[pick])
            model$calculate(y2.nodes)
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y1 <- model$getLogProb(y1.nodes[pick]) #will always be 0
          lp.initial.y2 <- model$getLogProb(y2.nodes)
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn pd on
          model$calculate(pd1.nodes[pick])
          model$calculate(pd2.nodes[pick])
          #don't use calculate, compute and insert manually
          # model$calculate(pd2.j.nodes)
          nondetect.probs.proposed <- nondetect.probs.initial*(1-model$pd2[pick,])
          model$pd2.j <<- 1 - nondetect.probs.proposed
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y1 <- model$calculate(y1.nodes[pick])
          lp.proposed.y2 <- model$calculate(y2.nodes)
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y1 + lp.proposed.y2) - (lp.initial.N + lp.initial.y1 + lp.initial.y2)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["pd1",1][pick,] <<- model[["pd1"]][pick,]
            mvSaved["pd2",1][pick,] <<- model[["pd2"]][pick,]
            mvSaved["pd2.j",1][1:J2] <<- model[["pd2.j"]][1:J2]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            nondetect.probs.initial  <- nondetect.probs.proposed
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["pd1"]][pick,] <<- mvSaved["pd1",1][pick,]
            model[["pd2"]][pick,] <<- mvSaved["pd2",1][pick,]
            model[["pd2.j"]][1:J2] <<- mvSaved["pd2.j",1][1:J2]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y1.nodes[pick])
            model$calculate(y2.nodes)
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    # copy(from = model, to = mvSaved, row = 1, nodes = z.nodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)