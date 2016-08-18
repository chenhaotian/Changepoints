co.intervals
x <- co.intervals(1:20)

## How robust they can be?
x <- rnorm(10000)
sd(x)
mad(x)
x[2] <- 100
sd(x)
mad(x)


## Hazard functions
exphazard <- function(x,rate=1){
    dexp(x,rate=rate)/(1-pexp(x,rate = rate))
}
exphazard(1:10)                         #constant hazard 


## integrate a scalar function over a multidimensional rectangle
library(cubature)
adaptIntegrate()

## direct calcualtion of predictivedistribution:
## P(x|hyperparameter)=Int_{theta}P(x|theta)P(theta|hyperparameter)

k <- function(k0,n){
    k0+n
}

mu <- function(mu0,k0,n,x){
    (k0*mu0+n*mean(x))/
        k(k0,n)
}

x <- rnorm(2,sd = 100)

mu0 <- 5
k0 <- 0.1
alpha0 <- 2
beta0 <- 3

mu(mu0,k0,2,x)
k(k0,2)

mu0 <- mu(mu0,k0,1,x[1])
k0 <- k(k0,1)
mu(mu0,k0,1,x[2])
k(k0,1)

alpha <- function(alpha0,n){
    alpha0+n/2
}

beta <- function(beta0,n,x,mu0,k0){
    beta0+1/2*sum((x-mean(x))^2)+k0*n*(mean(x)-mu0)^2/2/(k0+n)
}

alpha(alpha0,2)
beta(beta0,2,x,mu0,k0)

alpha0 <- alpha(alpha0,1)
beta0 <- beta(beta0,1,x[1],mu0,k0)
mu0 <- mu(mu0,k0,1,x[1])
k0 <- k(k0,1)

alpha(alpha0,1)
beta(beta0,1,x[2],mu0,k0)
mu(mu0,k0,1,x[2])
k(k0,1)

## online bayesian----------------
library(cubature)
X <- c(rnorm(1000),rnorm(700,mean=1),rnorm(700,mean = 2),rnorm(700),rnorm(700,mean = 1))
Y <- X+seq(1:length(X))*0.01
par(mfcol = c(2,1))
plot(x,type = "l")
plot(y,type = "l")

## Normal likelihood and normal-gamma prior distribution
## mu0, k0, alpha0, beta0: normal-gamma parameters
## MAXlength: indicating the maximum run length.
## lower/upperLimit: integrete limit in calculating predictive probability
## lambda: parameter of exponential hazard function.
onlinechangepoint <- function(X,MAXlength=1000L,
                              mu0=0,k0=1,alpha0=1/2,beta0=1,
                              lowerLimit=c(-1e2,0),upperLimit=c(1e2,1e2),
                              lambda=1000){
    if(MAXlength<2){
        stop("MAXlength must be a integer greater than 2!")
    }
    ## density function of normal-gamma distribution
    dnormalgamma <- function(mu,tau){
        dnorm(mu,mean = tracker["mu0"],sd=sqrt(1/(tracker["k0"]*tau)))*
            dgamma(tau,shape = tracker["alpha0"],rate = tracker["beta0"])
    }
    prediction <- function(eta){
        dnorm(x,mean = eta[1],sd = sqrt(1/eta[2]))*
            dnormalgamma(eta[1],eta[2])
    }
    hazard <- 1/lambda                  #constant hazard function
    
    ## initialize
    r <- 0:MAXlength                    #run lengths
    D_r <- rep(0,MAXlength+1)          #density of r, P(r_t|x_{1:t})
    D_rX <- c(1,rep(0,MAXlength))      #density of r and x_{1:t}, P(r_t,x_{1:t})
    P_x <- 0                           #probability of observing x, P(x|hyper)
    P_X <- 1                           #evidence, P(x_{1:t})


    r <- 0                             #run lengths
    D_r <- 1                           #density of r, P(r_t|x_{1:t})
    D_rX <- 1                          #density of r and x_{1:t}, P(r_t,x_{1:t})
    P_X <- 1                           #evidence, P(x_{1:t})

    tracker <- matrix(c(0,1,1,1,mu0,k0,alpha0,beta0),nrow = 1,
                      dimnames = list(NULL,
                                      c("r","D_r","D_rX","P_X","mu0","k0","alpha0","beta0")))

    P_x <- 0                           #probability of observing x, P(x|hyper)

    res <- matrix(0,nrow = length(X),ncol = MAXlength+1) #result container
    pos <- 1

    pb <- txtProgressBar(min = 1,max = length(X),style = 3)

    debugmu0 <- rep(0,length(X))
    
    for(x in X){
        lapply(R,function(r){
            
        })
        ## predictive probability
        P_x <- adaptIntegrate(f=prediction,lowerLimit = lowerLimit,upperLimit = upperLimit,maxEval = 10000)$integral
        ## growth probabilities and changepoint probability
        idx <- which(D_rX[-MAXlength]>0)
        if(length(idx)>0){
            ## growth probabilities
            D_rX[idx+1] <- D_rX[idx]*P_x*(1-hazard)
            ## changepoint probability
            D_rX[1] <- sum(D_rX[idx]*P_x*hazard)
        }
        ## evidence
        P_X <- sum(D_rX)
        ## run length distribution
        D_r <- D_rX/P_X
        ## update hyperparameters
        mu0 <- (k0*mu0+x)/(k0+1)
        alpha0 <- alpha0+1/2
        beta0 <- beta0+k0*(x-mu0)^2/2/(k0+1)
        k0 <- k0+1

        debugmu0[pos] <- mu0

        
        ## attach result
        res[pos,] <- D_r
        pos <- pos+1
        ## progress bar
        setTxtProgressBar(pb,pos)


    }
    ## res
    debugmu0
}

Z <- 1:10
testfun <- function(){
    print(z)
}

for(z in Z){
    testfun()
}
