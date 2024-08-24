# different packages used:
library(ltm)
library(vegan)
library(gllvm)
library(boral); 
library(gmf)
library(glmmTMB)
library(vegan)

eta.simulate <- function(n_sites, n_species, predictors=FALSE, n_lv=5) {
  X <- matrix(1,n_sites,1)
  B <- runif(n_species, -1, 1)
  if (predictors) {
    x1 <- rnorm(n_sites)
    x2 <- rexp(n_sites)
    X <- cbind(X,x1,x2)
    B <- rbind(B, matrix(runif(2*n_species,-1,1),2,n_species)) 
  }
  lam <- matrix(runif(n_lv*n_species,-2,2),n_species,n_lv)
  diag(lam) <- abs(diag(lam)) # diagonal of the loading matrix positive
  lam[upper.tri(lam)] <- 0 # upper triangle of the loading matrix set to 0
  eta0 <- X%*%B
  return(list(X=X, B=B, lambda=lam, eta0=eta0))
}

y.simulate <- function(eta, model="poisson") {
  n = dim(eta)[1]
  p = dim(eta)[2]
  y = matrix(0,n,p)
  z <- 1 
  while (z > 0 || z2>0) {
    
    if(model=="poisson"){
      for (i in 1:n) {
        y[i,] <- rpois(p, lambda=exp(eta[i,]))
      }}
    
    if(model=="probit"){
      for (i in 1:n) {
        y[i,] <- rbinom(p, prob=pnorm(eta[i,]),size=1)
      }}
    if(model=="logit"){
      for (i in 1:n) {
        y[i,] <- rbinom(p, prob=plogis(eta[i,]),size=1)
      }}
    
    z <- sum(rowSums(y>0)<2)  # data with rows full of zeros ignored  
    z2 <- sum(colSums(y>0)<2)  # data with columns full of zeros ignored  
  }
  y
}


gllvm.sim <- function(true.vals, j, n.sim=100, link="probit") {
  eta0.true = true.vals$eta0
  b.true = true.vals$B
  lam.true = true.vals$lam
  n = dim(eta0.true)[1]
  p = dim(eta0.true)[2]
  num.lv = dim(lam.true)[2]
  
  proc.U = proc.lam = beta = t.fit = etas = NULL
  i = j
  while(i <= (n.sim-1+j)) {
    set.seed(i)
    # draw the latent variable scores:
    U.true = matrix(rnorm(n*num.lv), n, num.lv)
    eta.true = eta0.true + U.true %*% t(lam.true)
    yi = y.simulate(eta.true, model="poisson")
    proc.lami = proc.Ui = t.fiti = rep(NA, 6)
    betai = matrix(NA, nrow=p, ncol=6)
    VAi = LAi = EVAi = MCMCi = PQLi = NULL
    # fit the models:
    # gllvm-VA
    t.fiti[1] = system.time(try({
      VAi = gllvm(yi, family="poisson", num.lv = num.lv)
      lami = t(t(VAi$params$theta)*VAi$params$sigma.lv)
      Ui = VAi$lvs
      proc.Ui[1] = vegan::procrustes(scale(U.true, scale=TRUE),scale(Ui, scale=TRUE))$ss
      proc.lami[1] = vegan::procrustes(scale(lam.true, scale=TRUE),scale(lami, scale=TRUE))$ss
      betai[,1] = VAi$params$beta0
      },silent = TRUE))[3] 
    # gllvm-LA
    t.fiti[2] = system.time(try({
      LAi = gllvm(yi, family="poisson", method="LA", num.lv = num.lv)
      lami = t(t(LAi$params$theta)*LAi$params$sigma.lv)
      Ui = LAi$lvs
      proc.Ui[2] = vegan::procrustes(scale(U.true, scale=TRUE),scale(Ui, scale=TRUE))$ss
      proc.lami[2] = vegan::procrustes(scale(lam.true, scale=TRUE),scale(lami, scale=TRUE))$ss
      betai[,2] = LAi$params$beta0
      },silent = TRUE))[3]
    # gllvm-EVA
    t.fiti[3] = system.time(try({
      EVAi = gllvm(yi, family="poisson", method="EVA", num.lv = num.lv)
      lami = t(t(EVAi$params$theta)*EVAi$params$sigma.lv)
      Ui = EVAi$lvs
      proc.Ui[3] = vegan::procrustes(scale(U.true, scale=TRUE),scale(Ui, scale=TRUE))$ss
      proc.lami[3] = vegan::procrustes(scale(lam.true, scale=TRUE),scale(lami, scale=TRUE))$ss
      betai[,3] = EVAi$params$beta0
      },silent = TRUE))[3]
    # boral
    t.fiti[4] = system.time(try({
      MCMCi = boral(yi, lv.control = list(num.lv=num.lv), family="poisson")
      lami = MCMCi$lv.coefs.mean[,2:(num.lv+1)]
      Ui = MCMCi$lv.mean
      proc.Ui[4] = vegan::procrustes(scale(U.true, scale=TRUE),scale(Ui, scale=TRUE))$ss
      proc.lami[4] = vegan::procrustes(scale(lam.true, scale=TRUE),scale(lami, scale=TRUE))$ss
      betai[,4] = MCMCi$lv.coefs.mean[,1]
      },silent = TRUE))[3]
    # gmf
    t.fiti[5] = system.time(try({
      PQLi = gmf(yi, family=poisson(), intercept=TRUE, verbose=FALSE, p=num.lv)
      lami = PQLi$v
      Ui = PQLi$u
      proc.Ui[5] = vegan::procrustes(scale(U.true, scale=TRUE),scale(Ui, scale=TRUE))$ss
      proc.lami[5] = vegan::procrustes(scale(lam.true, scale=TRUE),scale(lami, scale=TRUE))$ss
      betai[,5] = PQLi$beta
      },silent = FALSE))[3]
    # glmmTMB:
    t.fiti[6] = system.time(try({
      y0 = y00 = data.frame(yi)
      y0$id = 1:n
      # need to reshape the data to long format for glmmTMB
      y0l = reshape(y0, idvar="id", timevar="Species", times=colnames(y00),
                   varying=list(colnames(y00)), v.names="abund", direction="long")
      RRi = glmmTMB(abund ~ Species + 0 + rr(Species + 0|id, d=num.lv),
                                       data = y0l, family=poisson())
      Ui = ranef(RRi)$cond$id[,1:num.lv]
      proc.Ui[6] = vegan::procrustes(scale(U.true, scale=TRUE),scale(Ui, scale=TRUE))$ss
      ord0 = paste("SpeciesX", 1:p, sep="") # proper order of cols
      b0i = fixef(RRi)$cond
      ord1 = names(b0i)
      betai[,6] = b0i[ord0]
      lami = RRi$obj$env$report(RRi$fit$parfull)$fact_load[[1]]
      rownames(lami) = ord1
      lami = lami[ord0,]
      proc.lami[6] = vegan::procrustes(scale(lam.true, scale=TRUE),scale(lami, scale=TRUE))$ss
    },silent=FALSE))[3]
    
    t.fit = rbind(t.fit, t.fiti)
    proc.U = rbind(proc.U, proc.Ui)
    proc.lam = rbind(proc.lam, proc.lami)
    beta = rbind(beta, betai)
    etas = rbind(etas, eta.true)
    i = i+1
  }
  colnames(proc.U) = colnames(proc.lam) = colnames(t.fit) = 
    colnames(beta) = c("VA", "LA", "EVA", "MCMC", "PQL", "LA2")
  return(list(t.fit=t.fit, proc.U=proc.U, proc.lam=proc.lam, beta=beta, etas=etas))
}

set.seed(19082024)
# Do the simulations in parallel on a cluster
library(doMPI, quietly=T)
cl <- startMPIcluster()
registerDoMPI(cl)

ss<-seq(1,1000,5)

# a) n=50 obs/rows, m=50 species/cols, 5 LVs.
ptime1 <- system.time({
  tvals50x50 <- eta.simulate(n_sites = 50, n_species = 50, predictors = FALSE, n_lv = 5)
  rev.50.50 <- foreach(j = ss,.packages=c("vegan","MASS")) %dopar% {
    gllvm.sim(tvals50x50, j, n.sim=5)
  };
})[3]
save(rev.50.50,file = "./lvmrev_50x50_pois.Rdata")
save(tvals50x50,file="./tvals50x50_pois.Rdata")

# b) n=50 obs/rows, m=150 species/cols, 5 LVs.
ptime2 <- system.time({
  tvals50x150 <- eta.simulate(n_sites = 50, n_species = 150, predictors = FALSE, n_lv = 5)
  rev.50.150 <- foreach(j = ss,.packages=c("vegan","MASS")) %dopar% {
    gllvm.sim(tvals50x150, j, n.sim=5)
  };
})[3]
save(rev.50.150,file = "./lvmrev_50x150_pois.Rdata")
save(tvals50x150,file="./tvals50x150_pois.Rdata")

# c) n=200 obs/rows, m=50 species/cols, 5 LVs.
ptime3 <- system.time({
  tvals200x50 <- eta.simulate(n_sites = 200, n_species = 50, predictors = FALSE, n_lv = 5)
  rev.200.50 <- foreach(j = ss,.packages=c("vegan","MASS")) %dopar% {
    gllvm.sim(tvals200x50, j, n.sim=5)
  };
})[3]
save(rev.200.50,file = "./lvmrev_200x50_pois.Rdata")
save(tvals200x50,file="./tvals200x50_pois.Rdata")

# a) n=200 obs/rows, m=150 species/cols, 5 LVs.
ptime4 <- system.time({
  tvals200x150 <- eta.simulate(n_sites = 200, n_species = 150, predictors = FALSE, n_lv = 5)
  rev.200.150 <- foreach(j = ss,.packages=c("vegan","MASS")) %dopar% {
    gllvm.sim(tvals200x150, j, n.sim=5)
  };
})[3]
save(rev.200.150,file = "./lvmrev_200x150_pois.Rdata")
save(tvals200x150,file="./tvals200x150_pois.Rdata")