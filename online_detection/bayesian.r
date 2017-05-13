## online bayesian-------------------------
## model: specifying model
##    nng: normal evidence and normal-gamma prior
##    pg: poisson evidence and gamma prior
##    bb: binomial evidence and beta prior
##         ref: [1] https://en.wikipedia.org/wiki/Checking_whether_a_coin_is_fair
##    g: gamma evidence
##         ref: [1] DeGroot,Optimal Statistical Decisions, chapter 9
##              [2] https://en.wikipedia.org/wiki/Gamma_distribution
##              [3] Fink, D. 1995 A Compendium of Conjugate Priors. In progress report: Extension and enhancement of methods for setting data quality objectives. (DOE contract 95‑831).
## mu0, k0, alpha0, beta0: normal-gamma parameter, when model="nng"
##          alpha0, beta0: gamma parameters when model="pg"(these two names alpha0 and beta0 are shared with "nng")
##          alpha0, beta0: beta parameters when model="bb"
## bpmethod: bayesian prediction method. 'bruteforce': calculate the bayesian prediction probability directly by integrate over (6), Note:relevant code has been commented. 'mean': use the posterior mean of \eta to approximate the calculation.
## lower/upperLimit: only meaningful when bpmethod='bruteforce', specifying integrete limit of the parameters in the likelihood function. Use this limit to approximate the whole interval in calculating predictive probability:
##     int_{lowerLimit}^{upperLimit} {P(x_t|theta_t)P(theta_t|hyper_t)}
## lambda: parameter of the exponential hazard function.(transition probability)
## FILTER: if P(r_t|x_{1:t})<FILTER, this r_t will be omitted in the next calculation.
##
onlinechangepoint <- function(X,
                              model=c("nng","pg","bb"),
                              mu0=0,k0=1,alpha0=1/2,beta0=1,
                              bpmethod=c("mean","bruteforce"),
                              lowerLimit=c(-1,0),upperLimit=c(3,1e2),
                              lambda=1000, #exponential hazard
                              FILTER=1e-4){
    ## require(cubature)                       #integrate over a hyper rectangle, no longer needed after removing 'integrate over (6)' part
    model <- match.arg(model)
    bpmethod <- match.arg(bpmethod)
    hazard <- 1/lambda                  #constant hazard function(transition probability)
    if(model=="pg"){
        if(any(X<0)) stop("X must be integer greater than or equal to zero")
    }else if(model=="pg"){
        if(any(X!=0&X!=1)) stop("X must be a sequence of 0s and 1s, 1 for head, 0 for tail")
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

        ## --------------------!!DELETED!!--------------------
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
        ## --------------------!!DELETED!!--------------------
        
        ## if x~gamma(shape=alpha,rate=beta), then mean(x)=alpha/beta
        ## if x~beta(shape1=alpha,shape2=beta), then mean(x)=alpha/(alpha+beta)
        r_snapshot[,"ppredict"] <- if(model=="nng") dnorm(x,mean = r_snapshot[,"mu0"],sd = sqrt(r_snapshot[,"beta0"]/r_snapshot[,"alpha0"])) else if(model=="pg") dpois(x,lambda = r_snapshot[,"alpha0"]/r_snapshot[,"beta0"]) else if(model=="bb") dbinom(x,size = 1,prob = r_snapshot[,"alpha0"]/(r_snapshot[,"beta0"]+r_snapshot[,"alpha0"]))

        ## P(r+1,x_{1:t}) and P(0,x_{1:t})
        tmp <- r_snapshot[,"prx"]*r_snapshot[,"ppredict"]
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
        }else if(model=="bb"){
            if(x==0) r_snapshot[,"beta0"] <- r_snapshot[,"beta0"]+1
            else r_snapshot[,"alpha0"] <- r_snapshot[,"alpha0"]+1 #x==1
        }
        
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


## Forward algorithm for HMM
##     ref: [1] MLAPP p609-610, Algorithm 17.1
## return only the state distributions(which is different from Algorithm 17.1)
## X: numeric or matrix, observations, each row of X is an entry of a multi-variable observation
## transition: matrix, transition matrix, dim = NxN
## observation.model: character, model of observation distribution(evidence distribution), "n" normal, "p" poisson, "b" binomial. "m" multinomial.
## observation.params: list, observation model parameters for each state, length = N, each distributions' parameters must also be contained in a list.
## pi: initial state distribution
HMMFwd <- function(X,
                   transition=matrix(), #transition matrix in finite state HMM
                   observation.model=c("n","p","b","m"), #observation distribution
                   observation.params=list(),        #observation parameters
                   pi=numeric()    #initial state distribution
                   ){
    if(is.vector(X)) X <- matrix(X,ncol = 1)
    ## 1. parameter check
    if(length(unique(c(dim(transition),length(observation.params),length(pi))))>1){
        stop("dim(transition), length(observation.params) and length(pi) must equal to the number of hidden states")
    }
    if(length(observation.params)==0) stop("missing observation.params")
    if(unique(sapply(observation.params,class))!="list") stop("parameters for each observation distribution should also be contained in a list!")
    observation.model <- match.arg(observation.model)
    ## 2. preparation
    ## subroutine
    normalize <- function(l){
        l/sum(l)
    }
    N <- length(pi)                     #number of hidden states
    a <- matrix(0,nrow = nrow(X),ncol = N) #output, state distributions
    ## get observation distribution
    observation <- switch(observation.model,
                       n=dnorm,
                       p=dpois,
                       b=dbinom,
                       m=dmultinom
                       )
    ## 3.initialize
    evidence <- sapply(observation.params,function(l){
        par <- c(list(x=X[1,]),l)
        do.call(observation,par)
    },simplify = TRUE,USE.NAMES = FALSE)
    a[1,] <- normalize(evidence*pi)
    ## 4. main loop
    pb <- txtProgressBar(min = 2,max = nrow(X),style = 3)
    for(i in 2:nrow(X)){
        evidence <- sapply(observation.params,function(l){
            par <- c(list(x=X[i,]),l)
            do.call(observation,par)
        },simplify = TRUE,USE.NAMES = FALSE)
        a[i,] <- normalize(evidence*(base::crossprod(transition,a[i-1,])))
        setTxtProgressBar(pb,i)
    }
    cat("\n")
    a
}

## Forward-Backword algorithm for HMM
HMMFwdBack <- function(X,
                   transition=matrix(), #transition matrix in finite state HMM
                   observation.model=c("n","p","b","m"), #observation distribution
                   observation.params=list(),        #observation parameters
                   pi=numeric()    #initial state distribution
                   ){
    if(is.vector(X)) X <- matrix(X,ncol = 1)
    ## 1. parameter check
    if(length(unique(c(dim(transition),length(observation.params),length(pi))))>1){
        stop("dim(transition), length(observation.params) and length(pi) must equal to the number of hidden states")
    }
    if(length(observation.params)==0) stop("missing observation.params")
    if(unique(sapply(observation.params,class))!="list") stop("parameters for each observation distribution should also be contained in a list!")
    observation.model <- match.arg(observation.model)
    ## 2. preparation
    ## subroutine
    normalize <- function(l){
        l/sum(l)
    }
    N <- length(pi)                     #number of hidden states
    T <- nrow(X)                        #number of observations
    evidence <- matrix(0,nrow = T,ncol = N) #evidence probabilities for each state
    a <- matrix(0,nrow = T,ncol = N) #output, forward state distributions
    b <- matrix(0,nrow = T,ncol = N) #output, backward state result
    g <- matrix(0,nrow = T,ncol = N) #output, smoothed state distribution, g \propto a*b
    ## get observation distribution
    observation <- switch(observation.model,
                       n=dnorm,
                       p=dpois,
                       b=dbinom,
                       m=dmultinom
                       )

    ## Forward begin---------------------------------------------------------
    ## this part is the same with HMMFwd()
    ## 3.initialize forward filter
    evidence[1,] <- sapply(observation.params,function(l){
        par <- c(list(x=X[1,]),l)
        do.call(observation,par)
    },simplify = TRUE,USE.NAMES = FALSE)
    a[1,] <- normalize(evidence[1,]*pi)
    ## 4. forward loop
    cat("forward filtering...\n")
    pb <- txtProgressBar(min = 2,max = T,style = 3)
    for(i in 2:T){
        evidence[i,] <- sapply(observation.params,function(l){
            par <- c(list(x=X[i,]),l)
            do.call(observation,par)
        },simplify = TRUE,USE.NAMES = FALSE)
        a[i,] <- normalize(evidence[i,]*(base::crossprod(transition,a[i-1,])))
        setTxtProgressBar(pb,i)
    }
    cat("\n")
    ## Forward end-----------------------------------------------------------

    ## Backward begin--------------------------------------------------------
    b[T,] <- 1
    ## 5. backward loop
    cat("backward smoothing...\n")
    pb <- txtProgressBar(min = 1,max = T-1,style = 3)
    for(i in (T-1):1){
        b[i,] <- unname(drop(
                     transition%*%drop(evidence[i+1,]*b[i+1,])
                 ))
        setTxtProgressBar(pb,T-i)
    }
    cat("\n done\n")
    ## Backward end----------------------------------------------------------

    ## Smoothing-------------------------------------------------------------
    g <- t(apply(a*b,MARGIN = 1,normalize))
    invisible(list(a=a,b=b,g=g))
}

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
    require(grid)
    
    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        ## Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        ## Make each plot, in the correct location
        for (i in 1:numPlots) {
            ## Get the i,j matrix positions of the regions that contain this subplot
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

## temp function for plotting original series and detected run-length distribution
tmpplot <- function(resY,Y){
    require(ggplot2)
    pd <- data.frame()
    for(i in 1:length(resY)){
        pd <- rbind(pd,data.frame(x=i,y=resY[[i]][,"r"],alpha=resY[[i]][,"prx"]))
    }
    p1 <- ggplot(pd)+geom_point(aes(x=x,y=y,alpha=alpha),fill="black")+theme(legend.position = "none")
    p2 <- ggplot(data.frame(x=1:length(Y),y=Y))+geom_line(aes(x=x,y=y))
    multiplot(p2,p1,layout = matrix(c(1,2),nrow=2))
}

## HMMFwd examples------------------------------------

## 1. 'occasionally dishonest casino' from MLAPP p607

## initialize
pi <- c(dice1=0.5,dice2=0.5)
transition <- matrix(c(0.95,0.1,0.05,0.9),nrow = 2,ncol = 2,
                     dimnames = list(c("dice1","dice2"),c("dice1","dice2")))
observation.params <- list(
    dice1=list(size=1,prob=c(1/6,1/6,1/6,1/6,1/6,1/6)),
    dice2=list(size=1,prob=c(1/10,1/10,1/10,1/10,1/10,5/10))
)
## simulate 300 samples
X <- matrix(0,nrow = 300,ncol = 6)
a <- character(300)
state <- "dice1"
for(i in 1:300){
    if(state=="dice1"&runif(1)<=0.95){
        ## do nothing
    }else if(state=="dice1"){
        state <- "dice2"
    }else if(state=="dice2"&runif(1)<=0.9){
        ## do nothing
    }else{
        state <- "dice1"
    }
    a[i] <- state
    X[i,] <- drop(do.call(rmultinom,c(observation.params[[state]],list(n=1))))
}
## forward filter
res <- HMMFwd(X,transition,observation.model = "m",observation.params,pi)
res2 <- HMMFwdBack(X,transition,observation.model = "m",observation.params,pi)$g
## plot
tmp <- res[,2]
tmp[tmp<=0.5] <- 0
plot(tmp,type = "h")
abline(v=which(a=="dice2"),col="gray")

tmp2 <- res2[,2]
tmp2[tmp2<=0.5] <- 0
plot(tmp2,type = "h")
abline(v=which(a=="dice2"),col="gray")

## online changepoint examples----------------------------
X <- c(rnorm(100,sd=0.5),rnorm(70,mean=5,sd=0.5),rnorm(70,mean = 2,sd=0.5),rnorm(70,sd=0.5),rnorm(70,mean = 7,sd=0.5))
Y <- c(rnorm(100,sd=0.5),rnorm(70,sd=1),rnorm(70,sd=3),rnorm(70,sd=1),rnorm(70,sd=0.5))
Z <- c(rpois(100,lambda=5),rpois(70,lambda=7),rpois(70,lambda=5),rpois(70,lambda=6),rpois(70,lambda=5))
par(mfcol = c(3,1))
plot(X,type = "l")
plot(Y,type = "l")
plot(Z,type = "l")

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
                          lambda=10, #exponential hazard
                          FILTER=1e-3)

tmpplot(resX,X)
tmpplot(resY,Y)
tmpplot(resZ,Z)

## price and volume data
load("price_volume")
directions <- c(0,unname(sign(diff(pv[,"close"]))))
directions[directions==0] <- NA
directions <- zoo::na.locf(directions,na.rm = FALSE)
pv <- cbind(pv,directions)
pv <- na.omit(pv[pv[,"volume"]!=0,])

pv[,"volume"] <- log(pv[,"volume"])
me <- median(pv[,"volume"])
O <- rep(unname(ifelse(pv[,"directions"]==1,1,0)),times = ceiling(pv[,"volume"]/me))

unname(pv[,"volume"])*directions

resO <- onlinechangepoint(O,
                          model = "bb",
                          alpha0=50,beta0=50,
                          bpmethod = "mean",
                          lambda=100, #exponential hazard
                          FILTER=1e-3)

tmpplot(resO,O)

## generate coin toss sequence

Z <- unname(pv[,"volume"])
Z <- Z[Z!=0]


poisson_lambda <- round(mean(pv[,"volume"]))
plot(Z,type = "h")

## use log to scale large numbers
ceiling(log(Z)/median(log(Z)))

## gamma related estimate
x <- rgamma(10000,shape=5,rate=5)
library(MASS)    # may be loaded by default
fitdistr(x, "gamma", start=list(shape=1, rate=1))$estimate


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

