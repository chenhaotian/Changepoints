## online bayesian-------------------------
library(cubature)                       #integrate over a hyper rectangle
library(ggplot2)

X <- c(rnorm(100,sd=0.5),rnorm(70,mean=5,sd=0.5),rnorm(70,mean = 2,sd=0.5),rnorm(70,sd=0.5),rnorm(70,mean = 7,sd=0.5))
Y <- c(rnorm(100,sd=0.5),rnorm(70,sd=1),rnorm(70,sd=3),rnorm(70,sd=1),rnorm(70,sd=0.5))
par(mfcol = c(2,1))
plot(X,type = "l")
plot(Y,type = "l")

## Normal likelihood and normal-gamma prior distribution
## mu0, k0, alpha0, beta0: normal-gamma parameters
## bpmethod: bayesian prediction method. 'bruteforce': calculate the bayesian prediction probability directly by integrate over (6). 'mean': use the posterior mean of \eta to approximate the calculation.
## lower/upperLimit: only meaningful when bpmethod='bruteforce', specifying integrete limit of the parameters in the likelihood function. Use this limit to approximate the whole interval in calculating predictive probability.
## lambda: parameter of the exponential hazard function.
## FILTER: if P(r_t|x_{1:t})<FILTER, this r_t will be omitted in the next calculation.
onlinechangepoint <- function(X,
                              mu0=0,k0=1,alpha0=1/2,beta0=1,
                              bpmethod=c("mean","bruteforce"),
                              lowerLimit=c(-1,0),upperLimit=c(3,1e2),
                              lambda=1000, #exponential hazard
                              FILTER=1e-4){
    
    match.arg(bpmethod,c("bruteforce","mean"))
    hazard <- 1/lambda                  #constant hazard function

    ## initialize
    x_snapshot <- 1                     #P(x_{1:t})
    r_snapshot <- matrix(c(0,1,0,0,0,mu0,k0,alpha0,beta0,1),nrow = 1,
           dimnames=list(NULL,
                         c("r","p","ppredict","pgrow","pshrink","mu0","k0","alpha0","beta0","prx"))) #r: run length, p: un-normalized P(r_t,x_{1:t}), prx: normalized P(r_t|x_{1:t})

    res <- list()
    
    pb <- txtProgressBar(min = 1,max = length(X),style = 3)
    pos <- 1
    
    ## online update
    for(x in X){
        ## 1. general calculation
        ## P(x_{t+1}|hyper)
        if(bpmethod=="bruteforce"){
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
                ## direct calcualtion of predictivedistribution:
                ## P(x|hyperparameter)=Int_{theta}P(x|theta)P(theta|hyperparameter)
                ## integrate a scalar function over a multidimensional rectangle: library(cubature)
                adaptIntegrate(f=prediction,lowerLimit = lowerLimit,upperLimit = upperLimit,maxEval = 10000)$integral
            }))
        }
        else if(bpmethod=="mean"){
            ## if x~gamma(shape=alpha,rate=beta), then mean(x)=alpha/beta
            r_snapshot[,"ppredict"] <- dnorm(x,mean = r_snapshot[,"mu0"],sd = sqrt(r_snapshot[,"beta0"]/r_snapshot[,"alpha0"]))
        }
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

## online changepoint detection for series X
resX <- onlinechangepoint(X,
                         mu0=0.7,k0=1,alpha0=1/2,beta0=1,
                         bpmethod = "mean", #the brute-force method is too time consuming
                         lambda=50, #exponential hazard
                         FILTER=1e-3)

resX <- onlinechangepoint(X,
                         mu0=0,k0=1,alpha0=1/2,beta0=1,
                         bpmethod = "mean", #the brute-force method is too time consuming
                         lambda=50, #exponential hazard
                         FILTER=1e-3)


## online changepoint detection for series Y
resY <- onlinechangepoint(Y,
                          mu0=0,k0=0.5,alpha0=1/2,beta0=1,
                          bpmethod = "mean",
                          lambda=50, #exponential hazard
                          FILTER=1e-3)

## plot for X (same as Y)
pd <- data.frame()
for(i in 1:length(resX)){
    pd <- rbind(pd,data.frame(x=i,y=resX[[i]][,"r"],alpha=resX[[i]][,"prx"]))
}
p1 <- ggplot(pd)+geom_point(aes(x=x,y=y,alpha=alpha),fill="black")+theme(legend.position = "none")
p2 <- ggplot(data.frame(x=1:length(X),y=X))+geom_line(aes(x=x,y=y))

p3 <- ggplot(data.frame(x=1:length(tail(price,999)),y=tail(price,999)))+geom_line(aes(x=x,y=y))

## multiplot() from R-cookbook
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(p2,p1,layout = matrix(c(1,2),nrow=2))




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
