#Testing for Github.
# I'm goign to test it some more.  What is up with this?
# Okay, now I'm making change.

# Implementing feature A.  
#10/03/2012
# These constants depend on the global variables gv.  
W.init <- 4 # Initial Cake size
kAllw <- seq(0.01, W.init, 0.01)   # State Space
kAllc <- kAllw     # Control Space
kNumW <- length(kAllw)
kNumC <- length(kAllc)
kNumStates <- kNumW
kSinkIndex <- length(kAllw) + 1
kSinkUtil <- -1e15
theta <- 0.62
beta <- 0.8
B <- 1/(1-beta)
A <- B*(log(1-beta)+((beta/(1-beta))*log(beta)))
V.true <- A+B*log(0.5*kAllw)

# Discretize epsilons accoding to a log normal.
# N <- 50  # Number of grid size
# MAX <- 10
# MIN <- 0
# data.grid <- MIN:N*(MAX/N)
# kAlleps <- rep(NA, length(data.grid)-1)
# kNumEps <- length(kAlleps)
# for(i in 1:(length(data.grid)-1)) {
#   kAlleps[i] <- ((data.grid[i+1]-data.grid[i])/2)+data.grid[i]
# }
# kAlleps.prob <- dlnorm(kAlleps)
# library(actuar)
# kAlleps.prob2 <- discretize(plnorm, MIN, MAX, step=(MAX-MIN)/N, method="rounding")
# test <- sample(kAlleps, 1000, replace=TRUE, prob=kAlleps.prob)
# mean(test)
# sum(kAlleps*kAlleps.prob2)

### Begin functions
utilityFunc <- function(curr.w, new.w, eps, theta){
  #util <- ifelse(curr.w > new.w,
  #               eps*log(curr.w-(theta*new.w)),
  #               -34.54)   # log(1e-15)
  util <- eps*log(curr.w-(theta*new.w))
  
  return(util)
}

fpCalc <- function(theta) {
  V <- matrix(0, nrow=kNumStates+1)
  V.temp <- V + 1
  V.c.vec <- rep(NA, kNumC)
  V.eps.vec <- rep(NA, kNumEps)
  c.max.vec <- rep(NA, kNumEps)
  G <- matrix(NA, nrow=kNumStates+1, ncol=kNumEps, dimnames=list(c(kAllw, "Sink"), kAlleps))
  V[kSinkIndex] <- kSinkUtil
  V[1] <- log(1e-15)
  G[1,] <- rep("Sink", kNumEps)
  G["Sink", ] <- rep("Sink", kNumEps)
  
  iter <- 1
  while(max(abs(V[-c(1,kSinkIndex)]-V.temp[-c(1,kSinkIndex)])) > 1e-7) {
    cat("iteration: ", iter, "\n")
    V.temp <- V
    
    kW <- 2     
    for(wSize in kAllw[-1]) {
      kEps <- 1
      for(eps in kAlleps) {
        w.vec <- rep(wSize, kNumC)      
        eps.vec <- rep(eps, kNumC)
        V.c.vec <- ifelse(kAllc < w.vec,
                          utilityFunc(w.vec, kAllc, eps.vec, theta) + beta*V[-kSinkIndex],
                          rep(kSinkUtil, kNumC))
        V.eps.vec[kEps] <- max(V.c.vec)
        c.max.vec[kEps] <- kAllc[which.max(V.c.vec)]
        kEps <- kEps + 1
      }
      
      V[kW] <- sum(kAlleps.prob2*V.eps.vec)  
      G[kW, ] <- c.max.vec
      
      kW <- kW + 1
    }            
    iter = iter +1
  }   
      
  return(list("V"=V, "G"=G))
}
out <- fpCalc(theta)
V <- out$V
G <- out$G
plot(kAllw, V.true, type='l', ylim=c(-60, 0), main="Integrated Value Function", ylab="V", xlab="w")
lines(kAllw, V, col="red")
legend("bottomright", c("Analytical Solution", "VF Iteration"), col=c("black", "red"), pch=c("-", "-"))

generate.data <- function(W.init, n.periods) {
  w.t <- rep(W.init, n.periods)
  c.t <- rep(NA, n.periods)  # Choose next period state
  consume.t <- rep(NA, n.periods)
 
  out <- fpCalc(theta)  
  V <- out$V
  G <- out$G
  
  errors <- sample(kAlleps, size=n.periods, replace=TRUE, prob=kAlleps.prob2)
      
  for(iperiod in 1:n.periods) {
    #V.vec <- utilityFunc(w.t[iperiod], kAllc, errors[iperiod], theta) + beta*V
    c.t[iperiod] <- as.numeric(G[which(rownames(G)==w.t[iperiod]), which(colnames(G)==errors[iperiod])])
    consume.t[iperiod] <- w.t[iperiod]-c.t[iperiod]
              
    if(iperiod < n.periods){
      w.t[iperiod+1] <- c.t[iperiod]
    }
  }
  
  mydata <- data.frame("w"=w.t, "c"=c.t, "consume"=consume.t)
  return(mydata)
}

nfxp.lik <- function(theta, w, w.prime, out=NULL){ 
  if(is.null(out)) {
    out <- fpCalc(theta)       
  } 
  
  V <- out$V
  G <- out$G
  
  f.eps <- rep(NA, length(w))
  
  #eps = beta*(w[-length(w)]-theta*w.prime[-length(w)])/(w.prime[-length(w.prime)]-theta*w.prime[-1]) * mean(kAlleps)   
  for(i in 1:length(w)) {
    eps.vec <- G[which(rownames(G)==w[i]),]
    eps <- colnames(G)[which(eps.vec == w.prime[i])][1]    
    f.eps[i] <- ifelse(is.na(eps), 1e-15, kAlleps.prob2[which(kAlleps==eps)])
  }
  
  return(sum(log(f.eps)))  
}

########################################
kAlleps <- c(0.4, 1.5)
kAlleps.prob2 <- c(0.7, 0.3)
kNumEps <- length(kAlleps)
W.init <- 4 # Initial Cake size
kAllw <- seq(0.01, W.init, 0.01)   # State Space
kAllc <- kAllw     # Control Space
kNumW <- length(kAllw)
kNumC <- length(kAllc)
kNumStates <- kNumW
kSinkIndex <- length(kAllw) + 1
kSinkUtil <- -1e15
theta <- 0.7   # when theta is between 0.84 and 0.99
beta <- 0.8
my.data <- generate.data(W.init, 35)
sol <- optimize(f=nfxp.lik,
  interval=c(0.01,0.99),  
  w=my.data[, "w"], 
  w.prime=my.data[, "c"], 
  #out=out,
  maximum=TRUE)
sol
