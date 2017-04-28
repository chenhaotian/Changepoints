## online bayesian-------------------------
library(cubature)                       #integrate over a hyper rectangle
library(ggplot2)

X <- c(rnorm(100,sd=0.5),rnorm(70,mean=5,sd=0.5),rnorm(70,mean = 2,sd=0.5),rnorm(70,sd=0.5),rnorm(70,mean = 7,sd=0.5))
Y <- c(rnorm(100,sd=0.5),rnorm(70,sd=1),rnorm(70,sd=3),rnorm(70,sd=1),rnorm(70,sd=0.5))
Z <- c(rpois(100,lambda=5),rpois(70,lambda=10),rpois(70,lambda=30),rpois(70,lambda=10),rpois(70,lambda=5))
par(mfcol = c(3,1))
plot(X,type = "l")
plot(Y,type = "l")
plot(Z,type = "l")

## model: specifying model
##    nng:normal evidence and normal-gamma prior
##    pg: poisson evidence and gamma prior
## mu0, k0, alpha0, beta0: normal-gamma parameter, when model="nng"
## alpha0, beta0: gamma parameters when model="pg"(these two names alpha0 and beta0 are shared with "nng")
## bpmethod: bayesian prediction method. 'bruteforce': calculate the bayesian prediction probability directly by integrate over (6), Note:relevant code has been commented. 'mean': use the posterior mean of \eta to approximate the calculation.
## lower/upperLimit: only meaningful when bpmethod='bruteforce', specifying integrete limit of the parameters in the likelihood function. Use this limit to approximate the whole interval in calculating predictive probability:
##     int_{lowerLimit}^{upperLimit} {P(x_t|theta_t)P(theta_t|hyper_t)}
## lambda: parameter of the exponential hazard function.(transition probability)
## FILTER: if P(r_t|x_{1:t})<FILTER, this r_t will be omitted in the next calculation.
onlinechangepoint <- function(X,
                              model=c("nng","pg"),
                              mu0=0,k0=1,alpha0=1/2,beta0=1,
                              bpmethod=c("mean","bruteforce"),
                              lowerLimit=c(-1,0),upperLimit=c(3,1e2),
                              lambda=1000, #exponential hazard
                              FILTER=1e-4){

    model <- match.arg(model)
    bpmethod <- match.arg(bpmethod)
    hazard <- 1/lambda                  #constant hazard function(transition probability)

    if(model=="pg"){
        if(any(X<0)) stop("X must be integer greater than zero")
    }
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
        ## if(bpmethod=="bruteforce"){
        ##     r_snapshot[,"ppredict"] <- as.vector(apply(r_snapshot,1,function(l){
        ##         ## density function of normal-gamma distribution
        ##         dnormalgamma <- function(mu,tau){
        ##             dnorm(mu,mean = l["mu0"],sd=sqrt(1/(l["k0"]*tau)))*
        ##                 dgamma(tau,shape = l["alpha0"],rate = l["beta0"])
        ##         }
        ##         prediction <- function(eta){
        ##             dnorm(x,mean = eta[1],sd = sqrt(1/eta[2]))*
        ##                 dnormalgamma(eta[1],eta[2])
        ##         }
        ##         ## direct calcualtion of predictivedistribution:
        ##         ## P(x|hyperparameter)=Int_{theta}P(x|theta)P(theta|hyperparameter)
        ##         ## integrate a scalar function over a multidimensional rectangle: library(cubature)
        ##         adaptIntegrate(f=prediction,lowerLimit = lowerLimit,upperLimit = upperLimit,maxEval = 10000)$integral
        ##     }))
        ## }
        ## else if(bpmethod=="mean"){
        
        ## if x~gamma(shape=alpha,rate=beta), then mean(x)=alpha/beta
        r_snapshot[,"ppredict"] <- if(model=="nng") dnorm(x,mean = r_snapshot[,"mu0"],sd = sqrt(r_snapshot[,"beta0"]/r_snapshot[,"alpha0"])) else if(model=="pg") dpois(x,lambda = r_snapshot[,"alpha0"]/r_snapshot[,"beta0"])        

        ## P(r+1,x_{1:t}) and P(0,x_{1:t})
        tmp <- r_snapshot[,"p"]*r_snapshot[,"ppredict"]
        r_snapshot[,"pgrow"] <- tmp*(1-hazard)
        r_snapshot[,"pshrink"] <- tmp*hazard
        ## 2. grow
        ## move one step further
        r_snapshot[,"r"] <- r_snapshot[,"r"]+1
        r_snapshot[,"p"] <- r_snapshot[,"pgrow"]
        ## update hyperparameters

        ## get posterior distributions
        if(model=="nng"){
            r_snapshot[,"mu0"] <- (r_snapshot[,"k0"]*r_snapshot[,"mu0"]+x)/(r_snapshot[,"k0",drop=TRUE]+1)
            r_snapshot[,"alpha0"] <- r_snapshot[,"alpha0"]+1/2
            r_snapshot[,"beta0"] <- (r_snapshot[,"beta0"]+r_snapshot[,"k0"]*(x-r_snapshot[,"mu0"])^2/2/(r_snapshot[,"k0"]+1))
            r_snapshot[,"k0"] <- r_snapshot[,"k0"]+1
        }else if(model=="pg"){
            r_snapshot[,"alpha0"] <- r_snapshot[,"alpha0"]+x
            r_snapshot[,"beta0"] <- r_snapshot[,"beta0"]+1
        }
        ## r_snapshot[,"mu0"] <- if(model=="nng") (r_snapshot[,"k0"]*r_snapshot[,"mu0"]+x)/(r_snapshot[,"k0",drop=TRUE]+1)
        ##                       else if(model=="pg") r_snapshot[,"mu0"] #doesn't change
        ## r_snapshot[,"alpha0"] <- if(model=="nng") r_snapshot[,"alpha0"]+1/2
        ##                          else if(model=="pg") r_snapshot[,"alpha0"]+x
        ## r_snapshot[,"beta0"] <- if(model=="nng") (r_snapshot[,"beta0"]+r_snapshot[,"k0"]*(x-r_snapshot[,"mu0"])^2/2/(r_snapshot[,"k0"]+1))
        ##                         else if(model=="pg") r_snapshot[,"beta0"]+1
        ## r_snapshot[,"k0"] <- if(model=="nng") r_snapshot[,"k0"]+1
        ##                      else if(model=="pg") r_snapshot[,"k0"] #doesn't change
        
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

        if(any(is.na(r_snapshot))){
            print(r_snapshot)
            print(pos)
            stop()
        }

        pos <- pos+1
        setTxtProgressBar(pb,pos)
    }

    ## res <- matrix(0,nrow = length(X),ncol = MAXlength+1) #result container
    ## pos <- 1
    cat("\n")
    res
}

## online changepoint detection for series X
resX <- onlinechangepoint(X,
                          model = "nng",
                          mu0=0.7,k0=1,alpha0=1/2,beta0=1,
                          bpmethod = "mean", #the brute-force method is too time consuming
                          lambda=50, #exponential hazard
                          FILTER=1e-3)

## online changepoint detection for series Y
resY <- onlinechangepoint(Y,
                          model = "nng",
                          mu0=0,k0=0.5,alpha0=1/2,beta0=1,
                          bpmethod = "mean",
                          lambda=50, #exponential hazard
                          FILTER=1e-3)

## online changepoint detection for series Z
resZ <- onlinechangepoint(Z,
                          model = "pg",
                          alpha0=10,beta0=1,
                          bpmethod = "mean",
                          lambda=50, #exponential hazard
                          FILTER=1e-3)


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

## multi series merged plot
## example:
## x <- matrix(runif(20),ncol = 2,dimnames = list(NULL,c("y1","y2")))
## tsplot(x)
tsplot <- function(x){
    try(dev.off(which = dev.list()),silent = TRUE)
    x <- as.matrix(x)
    nlines <- ncol(x)
    par(mar=c(5, 4*nlines, 4, 4) + 0.1)
    N <- nrow(x)
    NAMES <- colnames(x)
    TIME <- 1:N

    PRE <- 2

    YLIM <- c(min(x[,1],na.rm = TRUE),max(x[,1],na.rm = TRUE))
    plot(TIME, x[,1], axes=F, ylim=YLIM, xlab="", ylab="",type="l",col="black", main="",xlim=c(1,N))
    ## points(TIME,x[,1],pch=20,col="black")
    axis(2, ylim=YLIM,col="black",lwd=2)
    mtext(2,text=NAMES[1],line=PRE)

    if(nlines>=2)
        for(i in 2:nlines){
            YLIM <- c(min(x[,i],na.rm = TRUE),max(x[,i],na.rm = TRUE))
            par(new=T)
            PRE <- PRE+1.5
            plot(TIME, x[,i], axes=F, ylim=YLIM, xlab="", ylab="", 
                 type="l",lty=i, main="",xlim=c(1,N),lwd=2)
            axis(2, ylim=YLIM,lwd=2,line=PRE)
            PRE <- PRE+2
            mtext(2,text=NAMES[i],line=PRE)
        }
    axis(1,pretty(range(TIME),10))
    mtext("TIME",side=1,col="black",line=2)
    legend(x=round(N*0.7),y=max(x[,nlines]),legend=NAMES,lty=1:nlines)
}

## plot for Y (same as X)
pd <- data.frame()
for(i in 1:length(resY)){
    pd <- rbind(pd,data.frame(x=i,y=resY[[i]][,"r"],alpha=resY[[i]][,"prx"]))
}
p1 <- ggplot(pd)+geom_point(aes(x=x,y=y,alpha=alpha),fill="black")+theme(legend.position = "none")
p2 <- ggplot(data.frame(x=1:length(Y),y=Y))+geom_line(aes(x=x,y=y))
multiplot(p2,p1,layout = matrix(c(1,2),nrow=2))

pd <- data.frame()
for(i in 1:length(resX)){
    pd <- rbind(pd,data.frame(x=i,y=resX[[i]][,"r"],alpha=resX[[i]][,"prx"]))
}
p1 <- ggplot(pd)+geom_point(aes(x=x,y=y,alpha=alpha),fill="black")+theme(legend.position = "none")
p2 <- ggplot(data.frame(x=1:length(X),y=X))+geom_line(aes(x=x,y=y))
multiplot(p2,p1,layout = matrix(c(1,2),nrow=2))



## price and volume data
load("price_volume")
Z <- unname(pv[,"volume"])
Z <- Z[Z!=0]
poisson_lambda <- round(mean(pv[,"volume"]))
plot(Z,type = "h")

## online changepoint detection for series Y
resZ <- onlinechangepoint(Z,
                          model = "pg",
                          alpha0=poisson_lambda,beta0=1,
                          bpmethod = "mean",
                          lambda=500, #exponential hazard
                          FILTER=1e-3)

pd <- data.frame()
for(i in 1:length(resZ)){
    pd <- rbind(pd,data.frame(x=i,y=resZ[[i]][,"r"],alpha=resZ[[i]][,"prx"]))
}

p1 <- ggplot(pd)+geom_point(aes(x=x,y=y,alpha=alpha),fill="black")+theme(legend.position = "none")
p2 <- ggplot(data.frame(x=1:length(Z),y=Z))+geom_bar(aes(x=x,y=y),stat = "identity")






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

