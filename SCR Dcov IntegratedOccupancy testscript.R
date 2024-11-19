library(nimble)
library(coda)
source("sim.SCR.Dcov.IntegratedOccupancy.R")
source("NimbleModel SCR Dcov IntegratedOccupancy.R")
source("NimbleFunctions SCR Dcov IntegratedOccupancy.R")
source("sSampler Dcov.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

#simulate some data
p0.SCR <- 0.25 #baseline detection probability - SCR observation model
p0.USCR <- 0.25 #baseline detection probability - USCR observation model
sigma <- 0.50 #detection spatial scale
K1 <- 5 #number of SCR occasions
K2 <- 5 #number of USCR occasions
buff <- 3 #state space buffer. Should be at least 3 sigma (generally).

#make an SCR trapping array surrounded by USCR traps
X <- as.matrix(expand.grid(1:10,1:10)) #all traps
center <- colMeans(X)
dists <- sqrt((center[1]-X[,1])^2+(center[2]-X[,2])^2)
dist.cutoff <- 3.75
X1 <- X[dists<dist.cutoff,] #SCR trapping array
X2 <- X[dists>dist.cutoff,] #USCR trapping array

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
X[,1] <- X[,1]-x.shift
X[,2] <- X[,2]-y.shift
X1[,1] <- X1[,1]-x.shift
X1[,2] <- X1[,2]-y.shift
X2[,1] <- X2[,1]-x.shift
X2[,2] <- X2[,2]-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#create a density covariate
# D.cov <- rep(NA,n.cells)
# for(c in 1:n.cells){
#   D.cov[c] <-  7*dSS[c,1] - 0.5*dSS[c,1]^2 + 7*dSS[c,2] - 0.5*dSS[c,2]^2
# }
# D.cov <- as.numeric(scale(D.cov))
# 
# image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="Covariate Value")
# points(X,pch=4,cex=0.75,col="lightblue")

#simulate a D.cov, higher cov.pars for large scale cov
#change seed to get new D.cov. trial and error to create one with good trapping array coverage
# set.seed(13210) #pretty good one
set.seed(13216)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(1000,1000),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]>12] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]>12] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- 0
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density")
points(X1,pch=4,cex=0.75) #SCR traps
points(X2,pch=1,cex=0.75) #USCR traps

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(13435341) #change seed for new data set
data <- sim.SCR.Dcov.IntegratedOccupancy(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
                  p0.SCR=p0.SCR,p0.USCR=p0.USCR,sigma=sigma,K1=K1,K2=K2,X1=X1,X2=X2,
                  xlim=xlim,ylim=ylim,res=res)

points(data$s,pch=16)

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
                       n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
                       x.vals=data$x.vals,y.vals=data$y.vals)

#Data augmentation level
M <- 500

#trap operation matrix
X1 <- data$X1
X2 <- data$X2
X <- rbind(data$X1,data$X2)
J1 <- nrow(X1)
K1D1 <- rep(data$K1,J1)
J2 <- nrow(X2)
K1D2 <- rep(data$K2,J2)

#Augment and initialize
n1 <- data$n1
y2D <- matrix(0,M,J1)
y2D[1:n1,] <- apply(data$y1,c(1,2),sum)
xlim <- data$xlim
ylim <- data$ylim
s.init<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
idx <- which(rowSums(y2D)>0) #switch for those actually caught
for(i in idx){
  trps<- matrix(X[y2D[i,]>0,1:2],ncol=2,byrow=FALSE)
  if(nrow(trps)>1){
    s.init[i,]<- c(mean(trps[,1]),mean(trps[,2]))
  }else{
    s.init[i,]<- trps
  }
}

#sum USCR counts over occasions
y2 <- rowSums(data$y2)

#If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
e2dist  <-  function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
getCell  <-  function(s,res,cells){
  cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
}
alldists <- e2dist(s.init,data$dSS)
alldists[,data$InSS==0] <- Inf
for(i in 1:M){
  this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
  if(data$InSS[this.cell]==0){
    cands <- alldists[i,]
    new.cell <- which(alldists[i,]==min(alldists[i,]))
    s.init[i,] <- data$dSS[new.cell,]
  }
}

#plot to make sure initialized activity centers are in habitat
image(data$x.vals,data$y.vals,matrix(data$InSS,data$n.cells.x,data$n.cells.y))
points(s.init,pch=16)


z.init <- 1*(rowSums(y2D)>0)
z.data <- rep(NA,M)
z.data[1:n1] <- 1

#inits for nimble - MUST use z init and N init for data augmentation scheme to work. should use s.init, too.
Niminits <- list(z=z.init,N=sum(z.init>0),s=s.init,
                 p0.SCR=0.5,p0.USCR=0.5,sigma=1,
                 D0=sum(z.init)/(sum(data$InSS)*data$res^2),D.beta1=0)

#constants for Nimble
#here, you probably want to center your D.cov. The one I simulated for this testscript is already centered.
# D.cov.use <- data$D.cov - mean(data$D.cov) #plug this into constants$D.cov if centering
constants <- list(M=M,J1=J1,J2=J2,K1D1=K1D1,K1D2=K1D2,xlim=xlim,ylim=ylim,
                D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                res=data$res)

#supply data to nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(y1=y2D,y2=y2,z=z.data,X1=X1,X2=X2,dummy.data=dummy.data,cells=data$cells,InSS=data$InSS)

# set parameters to monitor
parameters <- c('N','lambda.N','p0.SCR','p0.USCR','sigma','D0',"D.beta1")
parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting
nt <- 1 #thinning rate for parameters
nt2 <- 5 #thinning rate for parameters2

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
config.nodes <- c('p0.SCR','p0.USCR','sigma','D0',"D.beta1")
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2,thin2=nt2,nodes=config.nodes,
                      useConjugacy=FALSE) #configure very slow if checking for conjugacy


###*required* sampler replacement
z.ups <- round(M*0.25) # how many z proposals per iteration??? 25% of M generally seems good, but no idea what is optimal
conf$removeSampler("N")
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(inds.detected=1:n1,z.ups=z.ups,J2=J2,M=M),
                silent = TRUE)

conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  type = 'sSampler',control=list(i=i,J1=J1,J2=J2,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                 xlim=data$xlim,ylim=data$ylim),silent = TRUE)
}

#often better to block these 2 sets of parameters
#AF_slice mixes better, but runs more slowly, with reduced runtime a function of how
#costly it is to evaluate the likelihoods involved.
#AF_slice is fast for D covs
conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
#AF_slice causes larger slowdown in run time. May still be worth it. But RW_block is faster.
#Not sure which is most efficient in terms of effective sample size/unit time
conf$addSampler(target = c("p0.SCR","p0.USCR","sigma"),
                type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples <-as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[burnin:nrow(mvSamples),]))

#truth
data$N
data$lambda

mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10


#compare expected D plot to truth
#image will show 
#posterior means
lambda.cell.post <- cellArea*mvSamples2[burnin2:nrow(mvSamples2),D0.idx]*mvSamples2[burnin2:nrow(mvSamples2),lambda.cell.idx]
lambda.cell.ests <- colMeans(lambda.cell.post)
#remove non-habitat
lambda.cell.ests[InSS==0] <- NA
lambda.cell[InSS==0] <- NA

par(mfrow=c(1,1),ask=FALSE)
zlim <- range(c(lambda.cell,lambda.cell.ests),na.rm=TRUE) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)

#cell ests vs. truth. should plot CIs...
plot(lambda.cell.ests~lambda.cell,pch=16) #remove non-habitat
abline(0,1,col="darkred",lwd=2)

