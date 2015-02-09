#10/03/2012
# These constants depend on the global variables gv.  
W.init <- 2 # Initial Cake size
kAllw <- seq(0.01, W.init, 0.01)   # State Space
kAllc <- kAllw     # Control Space
kNumW <- length(kAllw)
kNumC <- length(kAllc)
kNumStates <- kNumW
kSinkUtil <- 1e-15
theta <- 1
beta <- 0.9
B <- 1/(1-beta)
A <- B*(log(1-beta)+((beta/(1-beta))*log(beta)))
V.true <- A+B*log(kAllw)

kAlleps <- 1
kNumEps <- length(kAlleps)
kAlleps.prob2 <- 1

### Begin functions
utilityFunc <- function(curr.w, new.w, eps, theta){
  util <- ifelse(curr.w > new.w,
                 eps*log(curr.w-(theta*new.w)),
                 log(kSinkUtil))
  return(util)
}

fpCalc <- function(theta) {
  V <- matrix(0, nrow=kNumStates)
  V.temp <- V + 1
  V.c.vec <- rep(NA, kNumC)
  V.eps.vec <- rep(NA, kNumEps)
  c.max.vec <- rep(NA, kNumEps)
  G <- matrix(NA, nrow=kNumStates, ncol=kNumEps, dimnames=list(kAllw, kAlleps))
 
  iter <- 1
  while(max(abs(V-V.temp)) > 1e-20) {
    cat("iteration: ", iter, "\n")
    V.temp <- V
    
    kW <-1     
    for(wSize in kAllw) {
      kEps <- 1
      for(eps in kAlleps) {
        w.vec <- rep(wSize, kNumC)      
        eps.vec <- rep(eps, kNumC)      
        V.c.vec = utilityFunc(w.vec, kAllc, eps.vec, theta) + beta*V
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

plot(kAllw, V.true, type='l', ylim=c(-80, -25), main="Integrated Value Function", ylab="V", xlab="w")
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
    c.t[iperiod] <- G[which(rownames(G)==w.t[iperiod]), which(colnames(G)==errors[iperiod])]
    consume.t[iperiod] <- w.t[iperiod]-c.t[iperiod]
              
    if(iperiod < n.periods){
      w.t[iperiod+1] <- c.t[iperiod]
    }
  }
  
  mydata <- data.frame("w"=w.t, "c"=c.t, "consume"=consume.t)
  return(mydata)
}

my.data <- generate.data(W.init, 20)

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
    f.eps[i] <- kAlleps.prob2[which(kAlleps==eps)]
  }
  
  return(sum(log(f.eps)))  
}

nfxp.lik(1.5, my.data[9:20, "w"], my.data[9:20, "c"], out)

sol <- optimize(f=nfxp.lik,
  interval=c(-1,10),  
  w=my.data[1:15, "w"], 
  w.prime=my.data[1:15, "c"], 
  out=out,
  maximum=TRUE)
