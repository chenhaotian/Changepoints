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
X <- c(rnorm(100),rnorm(70,mean=1),rnorm(70,mean = 2),rnorm(70),rnorm(70,mean = 1))
Y <- X+seq(1:length(X))*0.01
par(mfcol = c(2,1))
plot(X,type = "l")
plot(Y,type = "l")

## Normal likelihood and normal-gamma prior distribution
## mu0, k0, alpha0, beta0: normal-gamma parameters
## MAXlength: indicating the maximum run length.
## lower/upperLimit: integrete limit in calculating predictive probability
## lambda: parameter of exponential hazard function.
## FILTER: if P(r_t|x_{1:t})<FILTER, this r_t will be omitted in the next calculation.
onlinechangepoint <- function(X,
                              mu0=0,k0=1,alpha0=1/2,beta0=1,
                              lowerLimit=c(-1,0),upperLimit=c(3,1e2),
                              lambda=1000, #exponential hazard
                              FILTER=1e-4){
    hazard <- 1/lambda                  #constant hazard function

    ## initialize
    x_snapshot <- 1                     #P(x_{1:t})
    r_snapshot <- matrix(c(0,1,0,0,0,mu0,k0,alpha0,beta0,1),nrow = 1,
           dimnames=list(NULL,
                         c("r","p","ppredict","pgrow","pshrink","mu0","k0","alpha0","beta0","prx"))) #r: run length, p: P(r_t,x_{1:t}), prx: P(r_t|x_{1:t})

    res <- list()
    
    pb <- txtProgressBar(min = 1,max = length(X),style = 3)
    pos <- 1
    
    ## online update
    for(x in X){
        ## 1. general calculation
        ## P(x_{t+1}|hyper)
        r_snapshot[,"ppredict"] <- as.vector(apply(r_snapshot,1,function(l){
            ## density function of normal-gamma distribution
            dnormalgamma <- function(mu,tau){
                dnorm(mu,mean = l["mu0"],sd=sqrt(1/(l["k0"]*tau)))*
                    dgamma(tau,shape = l["alpha0"],rate = l["beta0"])
            }
            prediction <- function(eta){
                dnorm(x,mean = eta[1],sd = sqrt(1/eta[2]))*
                    dnormalgamma(eta[1],eta[2])
            }
            adaptIntegrate(f=prediction,lowerLimit = lowerLimit,upperLimit = upperLimit,maxEval = 10000)$integral
        }))
        ## P(r+1,x_{1:t}) and P(0,x_{1:t})
        tmp <- r_snapshot[,"p"]*r_snapshot[,"ppredict"]
        r_snapshot[,"pgrow"] <- tmp*(1-hazard)
        r_snapshot[,"pshrink"] <- tmp*hazard
        ## 2. grow
        ## move one step further
        r_snapshot[,"r"] <- r_snapshot[,"r"]+1
        r_snapshot[,"p"] <- r_snapshot[,"pgrow"]
        ## update hyperparameters
        r_snapshot[,"mu0"] <- (r_snapshot[,"k0"]*r_snapshot[,"mu0"]+x)/
            (r_snapshot[,"k0",drop=TRUE]+1)
        r_snapshot[,"alpha0"] <- r_snapshot[,"alpha0"]+1/2
        r_snapshot[,"beta0"] <- r_snapshot[,"beta0"]+
            r_snapshot[,"k0"]*(x-r_snapshot[,"mu0"])^2/2/(r_snapshot[,"k0"]+1)
        r_snapshot[,"k0"] <- r_snapshot[,"k0"]+1
        ## 3. shrink
        r_snapshot <- rbind(
            matrix(c(0,sum(r_snapshot[,"pshrink"]),0,0,0,mu0,k0,alpha0,beta0,1),nrow = 1,dimnames=list(NULL,c("r","p","ppredict","pgrow","pshrink","mu0","k0","alpha0","beta0","prx"))),
            r_snapshot
        )
        ## 4. evidence P(x_{1:t}) and conditional probabiity P(r_t|x_{1:t})
        x_snapshot <- sum(r_snapshot[,"p"])
        r_snapshot[,"prx"] <- r_snapshot[,"p"]/x_snapshot
        ## 5. filter low probability run lengths
        r_snapshot <- r_snapshot[r_snapshot[,"prx"]>FILTER,,drop=FALSE]
        
        res <- c(res,list(r_snapshot[,c("r","prx")]))

        pos <- pos+1
        setTxtProgressBar(pb,pos)
    }

    ## res <- matrix(0,nrow = length(X),ncol = MAXlength+1) #result container
    ## pos <- 1


    res
}

res <- onlinechangepoint(X,
                         mu0=0.7,k0=1,alpha0=1/2,beta0=1,
                         lowerLimit=c(-10,0),upperLimit=c(10,1e2),
                         lambda=1000, #exponential hazard
                         FILTER=1e-4)

b <- c(1,2)
m <- c(3,7)
plot(1:10,type="l")
for(i in 1:2){
    points(b[i],m[i],pch=20)
}

pd <- data.frame()
for(i in 1:length(res)){
    pd <- rbind(pd,data.frame(x=i,y=res[[i]][,"r"],alpha=res[[i]][,"prx"]))
}
library(ggplot2)
ggplot(pd)+geom_point(aes(x=x,y=y,alpha=alpha),fill="black")

onlinechangepoint <- function(X,
                              mu0=0,k0=1,alpha0=1/2,beta0=1,
                              lowerLimit=c(-1e2,0),upperLimit=c(1e2,1e2),
                              lambda=1000,
                              FILTER=1e-4){
    ## density function of normal-gamma distribution
    dnormalgamma <- function(mu,tau){
        dnorm(mu,mean = l["mu0"],sd=sqrt(1/(l["k0"]*tau)))*
            dgamma(tau,shape = l["alpha0"],rate = l["beta0"])
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

    ## snapshot at t0
    snapshot <- matrix(c(0,1,1,1,mu0,k0,alpha0,beta0),nrow = 1,
                      dimnames = list(NULL,
                                      c("r","D_r","D_rX","P_X","mu0","k0","alpha0","beta0")))


    snapshot <- matrix(c(0,1,1,1,mu0,k0,alpha0,beta0),nrow = 1,
                       dimnames = list(NULL,
                                       c("r","D_r","D_rX","P_X","mu0","k0","alpha0","beta0")))
    

    x_snapshot <- 1                     #P(x_{1:t})
    ## prx: P(r|x_{1:t})
    ## p(r=0|x_0)=1
    r_snapshot <- matrix(c(0,1,0,0,0,mu0,k0,alpha0,beta0,1),nrow = 1,
           dimnames=list(NULL,
                         c("r","p","ppredict","pgrow","pshrink","mu0","k0","alpha0","beta0","prx")))
    for(x in X){
        ## 1. general calculation
        ## P(x_{t+1}|hyper)
        r_snapshot[,"ppredict"] <- as.vector(apply(r_snapshot,1,function(l){
            adaptIntegrate(f=prediction,lowerLimit = lowerLimit,upperLimit = upperLimit,maxEval = 10000)$integral
        }))
        ## P(r+1,x_{1:t}) and P(0,x_{1:t})
        tmp <- r_snapshot[,"p"]*r_snapshot[,"ppredict"]
        r_snapshot[,"pgrow"] <- tmp*hazard
        r_snapshot[,"pshrink"] <- tmp*(1-hazard)
        ## 2. grow
        ## move one step further
        r_snapshot[,"r"] <- r_snapshot[,"r"]+1
        r_snapshot[,"p"] <- r_snapshot[,"pgrow"]
        ## update hyperparameters
        r_snapshot[,"mu0"] <- (r_snapshot[,"k0"]*r_snapshot[,"mu0"]+x)/
            (r_snapshot[,"k0",drop=TRUE]+1)
        r_snapshot[,"alpha0"] <- r_snapshot[,"alpha0"]+1/2
        r_snapshot[,"beta0"] <- r_snapshot[,"beta0"]+
            r_snapshot[,"k0"]*(x-r_snapshot[,"mu0"])^2/2/(r_snapshot[,"k0"]+1)
        r_snapshot[,"k0"] <- r_snapshot[,"k0"]+1
        ## 3. shrink
        r_snapshot <- rbind(
            matrix(c(0,sum(r_snapshot[,"pshrink"]),0,0,0,mu0,k0,alpha0,beta0,1),nrow = 1,dimnames=list(NULL,c("r","p","ppredict","pgrow","pshrink","mu0","k0","alpha0","beta0","prx"))),
            r_snapshot
        )
        ## 4. evidence P(x_{1:t}) and conditional probabiity P(r_t|x_{1:t})
        x_snapshot <- sum(r_snapshot[,"p"])
        r_snapshot[,"prx"] <- r_snapshot[,"p"]/x_snapshot
        ## 5. filter low probability run lengths
        x_snapshot <- x[x_snapshot[,"prx"]>FILTER,,drop=FALSE]
    }

    res <- matrix(0,nrow = length(X),ncol = MAXlength+1) #result container
    pos <- 1

    pb <- txtProgressBar(min = 1,max = length(X),style = 3)

    debugmu0 <- rep(0,length(X))
    
    for(x in X){
        apply(snapshot,function(r){
            r
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
