#10/17/2012
# These constants depend on the global variables gv.  
W.init <- 2 # Initial Cake size
kAllw <- seq(0.01, W.init, 0.01)   # State Space
kAllc <- kAllw     # Control Space
kNumW <- length(kAllw)
kEvenW <- kAllw[seq(2, kNumW, 2)]
kNumC <- length(kAllc)
kNumStates <- kNumW
#kSinkIndex <- length(kAllw) + 1
kSinkUtil <- -1e15
theta <- 1
beta <- 0.8
B <- 1/(1-beta)
A <- B*(log(1-beta)+((beta/(1-beta))*log(beta)))
V.true <- A+B*log(kAllw)

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
#
### Begin functions
utilityFunc <- function(curr.w, new.w, eps, theta){
  util <- ifelse(curr.w > new.w,
                 eps*log(curr.w-(theta*new.w)),
                 log(1e-15))
  return(util)
}

fpCalc <- function(theta) {
  V.sp <- matrix(0, nrow=kNumStates)
  V.sp.t <- V.sp + 1
  V <- matrix(0, nrow=length(kEvenW))
  V.temp <- V + 1
  
  V.c.vec <- rep(NA, kNumC)
  V.eps.vec <- rep(NA, kNumEps)
  c.max.vec <- rep(NA, kNumEps)
  G <- matrix(NA, nrow=kNumStates, ncol=kNumEps, dimnames=list(kAllw, kAlleps))
     
  iter <- 1
  while(max(abs(V.sp-V.sp.t)) > 1e-20) {
    cat("iteration: ", iter, "\n")
    V.sp.t <- V.sp
    
    kW <- 1     
    for(wSize in kEvenW) {
      kEps <- 1
      for(eps in kAlleps) {
        w.vec <- rep(wSize, kNumC)      
        eps.vec <- rep(eps, kNumC)
        #V.c.vec <- ifelse(kAllc < w.vec,
        #                  utilityFunc(w.vec, kAllc, eps.vec, theta) + beta*V.sp,
        #                  rep(kSinkUtil, kNumC))
        V.c.vec <- utilityFunc(w.vec, kAllc, eps.vec, theta) + beta*V.sp.t
        V.eps.vec[kEps] <- max(V.c.vec)
        c.max.vec[kEps] <- kAllc[which.max(V.c.vec)]
        kEps <- kEps + 1
      }
      
      V[kW] <- sum(kAlleps.prob2*V.eps.vec)  
      G[kW, ] <- c.max.vec
      
      kW <- kW + 1
    }

    V.spline <- splinefun(kEvenW, V)            
    V.sp <- V.spline(kAllw)
    
    iter = iter +1
  }   
      
  return(list("V"=V.sp, "G"=G))
}

out <- fpCalc(theta)
kEvenW <- kAllw[seq(2, kNumW, 10)]
length(kEvenW)
V.20 <- fpCalc(theta)$V
kEvenW <- kAllw[seq(2, kNumW, 5)]
length(kEvenW)
V.40 <- fpCalc(theta)$V
kEvenW <- kAllw[seq(2, kNumW, 3.25)]
length(kEvenW)
V.60 <- fpCalc(theta)$V
kEvenW <- kAllw[seq(2, kNumW, 2.5)]
length(kEvenW)
V.80 <- fpCalc(theta)$V
kEvenW <- kAllw[seq(2, kNumW, 2)]
length(kEvenW)
V.100 <- fpCalc(theta)$V

G <- out$G
plot(kAllw, V.true, type='l', ylim=c(-40, -8), main="Integrated Bellman Function Using VF Iteration", ylab="V", xlab="w", sub="200 Discrete Points")
lines(kAllw, V, col="cyan")
lines(kAllw, V.20, col="red")
lines(kAllw, V.40, col="green")
lines(kAllw, V.60, col="blue")
lines(kAllw, V.80, col="yellow")
lines(kAllw, V.100, col="purple")
legend("bottomright", c("Analytical Solution", 
                        "No Splines", 
                        "20 Spline Nodes",
                        "40 Spline Nodes",
                        "60 Spline Nodes",
                        "80 Spline Nodes",
                        "100 Spline Nodes"), 
                      col=c("black", 
                            "cyan",
                            "red",
                            "green",
                            "blue",
                            "yellow",
                            "purple"), pch=c("-", "-", "-", "-", "-", "-"))

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
kAlleps <- 1
kAlleps.prob2 <- 1
kNumEps <- length(kAlleps)
W.init <- 4 # Initial Cake size
kAllw <- seq(0.01, W.init, 0.01)   # State Space
kAllc <- kAllw     # Control Space
kNumW <- length(kAllw)
kNumC <- length(kAllc)
kNumStates <- kNumW
kSinkIndex <- length(kAllw) + 1
kSinkUtil <- -1e15
theta <- 0.75   # when theta is between 0.84 and 0.99
beta <- 0.8

my.data <- generate.data(W.init, 35)
sol <- optimize(f=nfxp.lik,
  interval=c(0.01,0.99),  
  w=my.data[, "w"], 
  w.prime=my.data[, "c"], 
  #out=out,
  maximum=TRUE)
sol
