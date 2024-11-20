NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  #detection priors
  p0.SCR ~ dunif(0,1)
  p0.USCR ~ dunif(0,1)
  sigma ~ dunif(0,100)
  #--------------------------------------------------------------
  #Density model
  #splitting D0 off here removes need to consider s.cell likelihoods in D0 update
  D.intercept <- D0*cellArea
  # D.intercept <- exp(D.beta0)*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N) #realized N in state space
  for(i in 1:M){
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~  dunif(xlim[1],xlim[2])
    s[i,2] ~  dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]],InSS=InSS[s.cell[i]])
    #SCR Observation model, skipping z_i=0 calculations
    #SCR i x j detection probabilities
    pd1[i,1:J1] <- GetDetectionProb(s = s[i,1:2], X = X1[1:J1,1:2], J=J1,sigma=sigma, p0=p0.SCR, z=z[i])
    y1[i,1:J1] ~ dBinomialVector(pd=pd1[i,1:J1],K1D1[1:J1],z=z[i]) #vectorized obs mod
    #USCR i x j detection probabilities, skipping z_i=0 calculations
    pd2[i,1:J2] <- GetDetectionProb(s = s[i,1:2], X = X2[1:J2,1:2], J=J2,sigma=sigma, p0=p0.USCR, z=z[i])
  }
  #USCR Observation model, marginalized over individuals
  for(j in 1:J2){
    pd2.j[j] <- 1-(prod(1-pd2[1:M,j]*z[1:M]))
    y2[j] ~ dbinom(pd2.j[j],K1D2[j]) #occupancy data
  }
})
#custom Metropolis-Hastings update for N/z
