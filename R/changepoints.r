## Bayesian online change point detection and offline learning

## Inference in joint Gaussian distribution
## Get P(X1|X2=x2) where P(X1,X2) is jointly Gaussian
## returns the conditional mean and coviance matrix of P(x_1|x_2)
## 
##
## x2; numeric, an sample of X2, satisfying length(x2)<D, D is the joint dimension.
## mu: numeric, length D mean vector. mu=c(mu_X1,mu_X2)
## Precision: DxD precision matrix, satisfying Precision = inverse(Sigma). At least one of Sigma and Precision should be non-NULL
## Sigma: DxD covariance matrix. At least one of Sigma and Precision should be non-NULL.
##    Sigma = | sigma1   sigma_12 |
##            | sigma_21 sigma2   |
##
## Examples:
## tmp <- matrix(runif(100),20,5)
## S <- crossprod(tmp)                 #some synthetic covariance matrix
## P <- solve(S)
## m <- runif(5)
## x2 <- runif(3)
## inferenceJointGaussian(x2 = x2,mu = m,Precision = P)
## 
inferenceJointGaussian <- function(x2,mu,Precision=NULL,Sigma=NULL){
    if(is.null(Sigma) & is.null(Precision)) stop("Error in inferenceJointGaussian(): At least one of Sigma and Precision should be non-NULL")
    if(is.null(Precision)) Precision <- solve(Sigma)
    ## get the dimensions of x_1 and x_2
    D <- nrow(Precision)
    if(D!=length(mu)) stop("Error in inferenceJointGaussian(): Dimension doesn't mathc! nrow(Precision) != length(mu)")
    D1 <- D-length(x2)                  #dimension of X1
    if(D1<=0) stop("Error in inferenceJointGaussian(): length(mu2) should be strictly smaller than D")
    
    ## S11 <- Sigma[1:D1,1:D1,drop=FALSE]
    ## S12 <- Sigma[1:D1,(D1+1):D,drop=FALSE]
    ## S22 <- Sigma[(D1+1):D,(D1+1):D,drop=FALSE]
    ## S21 <- Sigma[(D1+1):D,1:D1,drop=FALSE]
    ## S11-S12%*%solve(S22)%*%S21
    
    P11 <- Precision[1:D1,1:D1,drop=FALSE]
    P12 <- Precision[1:D1,(D1+1):D,drop=FALSE]

    Sigma12 <- solve(P11)
    mu12 <- as.vector(mu[1:D1] - Sigma12%*%P12%*%(x2-mu[(D1+1):D]))

    list(mu12=mu12,Sigma12=Sigma12)
}

## Gaussian posterior distribution (NIW prior)
## Get P(theta|X) where X is Gaussian observations and P(theta) is NIW (Normal-Inverse-Wishart) with prior parameter m0,k0,v0,S0
## returns the posterior parameter mN,kN,vN,SN
## 
## x: observation matrix, or the ones that can be converted to matrix
## w: sample weights, default NULL
## m0,k0,v0,S0: prior Inverse-Whishart parameters
## 
posteriorGaussian <- function(x,w=NULL,m0,k0,v0,S0){
    x <- as.matrix(x)
    if(is.null(w)) w <- rep(1,nrow(x))
    N <- sum(w)                         #effective number of samples

    ## sample sufficient statistics
    mu <- colSums(x*w)/N                #sample mean
    x <- sweep(x,2,mu,`-`)
    S <- t(w*x)%*%x                     #sample scatter matrix centered on sample mean

    ## posterior parameters
    kN <- k0+N        
    mN <- (k0*m0+N*mu)/kN
    vN <- v0+N
    SN <- S0+S+k0*N/(kN)*(mu-m0)%*%t(mu-m0)

    list(mN=mN,kN=kN,vN=vN,SN=SN,N=N)
}

## Gaussian MAP esitmate (NIW prior)
## Get theta_MAP from NIW(theta),where theta_MAP = mod(NIW(theta))
## return the MAP estimates of mu and Sigma
##
## 
## m,k,v,S: NIW parameters
mapGaussian <- function(m,k,v,S){
    D <- length(m)                      #dimension
    list(muMAP=m,
         sigmaMAP=S/(v+D+2))
}

## Marginal likelihood of Gaussian observations (NIW prior)
## return the log marginal likelihood log P(x|model)
##
## x: observation matrix, or the ones that can be converted to matrix. each row of x is an observation
##
## Examples:
## m <- c(1,2,3)
## k <- 3
## v <- 3
## S <- matrix(c(2,0.1,0.2,0.1,2,0.3,0.2,0.3,2),3,3)
## x <- mvtnorm::rmvnorm(200,mean = c(2,3,2))
## marginalLikelihoodGaussian(x = x,m0=m,k0=k,v0=v,S0=S)
## 
marginalLikelihoodGaussian <- function(x,m0,k0,v0,S0,LOG=TRUE){
    x <- as.matrix(x)
    ## get the posterior parameters and sample size
    post <- posteriorGaussian(x=x,m0=m0,k0=k0,v0=v0,S0=S0)
    ## put the posteriors into marginalLikelihoodGaussian_byPosterior()
    marginalLikelihoodGaussian_byPosterior(N=post$N,
                                           mN=post$mN,
                                           kN=post$kN,
                                           vN=post$vN,
                                           SN=post$SN,
                                           m0,k0,v0,S0,LOG=LOG)
}
## Marginal likelihood of Gaussian observations (NIW prior)
## same as marginalLikelihoodGaussian(), difference is the input x is replaced by it's posteriors
marginalLikelihoodGaussian_byPosterior <- function(N,mN,kN,vN,SN,m0,k0,v0,S0,LOG=TRUE){
    ## copied from myTools.r
    ## tested, stable
    lmvgamma <- function(a,p){
        sapply(a,function(ai){
            p*(p-1)/4*log(pi)+ sum(lgamma(ai+(1-(1:p))/2))
        },simplify = TRUE)
    }
    det <- function(m){                 #det doesn't support scalar values, so write a wrapper around it
        if(is.matrix(m)){
            base::det(m)
        }else if(is.vector(m)){
            base::prod(m)
        }
    }
    p <- length(mN)                     #dimension
    
    ## method 1 use the equation
    logp <- -N*p/2*log(pi) + p/2*log(k0/kN) + v0/2*log(det(S0)) - vN/2*log(det(SN)) + lmvgamma(vN/2,p) - lmvgamma(v0/2,p)
    if(!LOG) logp <- exp(logp)
    
    ## method 2 use the t distribution
    ## to be completed
    
    logp
}

## Gaussian posterior prediction (NIW prior)
## return log P(x_new|x,m,k,v,S), note that when x is empty, "prediction problem" is equivilant to "marginal likelihood".
## 
## x_new: observation matrix, or the ones that can be converted to matrix.
## x: observation matrix, or the ones that can be converted to matrix. Default NULL
## m,k,v,S: prior Inverse-Whishart parameters
##
## Examples:
##
## m <- c(1,2,3)
## k <- 3
## v <- 3
## S <- matrix(c(2,0.1,0.2,0.1,2,0.3,0.2,0.3,2),3,3)
## x <- mvtnorm::rmvnorm(200,mean = c(2,3,2))
## xnew <- mvtnorm::rmvnorm(20,mean = c(2,3,2))
## posteriorPredictionGaussian(xnew = xnew,x = x,m=m,k=k,v=v,S=S)
posteriorPredictionGaussian <- function(xnew,x=NULL,m,k,v,S,LOG=TRUE){
    ## There are two methods
    ## method 1 and method 2 will produce identical result, 2 is more efficient
    
    ## method 1:
    ## exploit the fact that p(xnew|x,model) = p(xnew,x|model)/p(x|model)
    ## marginalLikelihoodGaussian(x=rbind(xnew,x),m0=m,k0=k,v0=v,S0=S,LOG=LOG) - marginalLikelihoodGaussian(x=x,m0=m,k0=k,v0=v,S0=S,LOG=LOG)

    ## method 2:
    ## use x to update model, then use the updated model to predix xnew
    if(!is.null(x)){
        x <- as.matrix(x)
        post <- posteriorGaussian(x=x,m0=m,k0=k,v0=v,S0=S)
        m <- post$mN
        k <- post$kN
        v <- post$vN
        S <- post$SN
    }
    xnew <- as.matrix(xnew)
    marginalLikelihoodGaussian(x=xnew,m0=m,k0=k,v0=v,S0=S,LOG=LOG)
}

## Linear Gaussian systems
## knowing that
##    P(X1) ~ N(mu1,Sigma1)
##    P(X2|X1) ~ N(AX1+b,Sigma21)
## get conditional distribution P(X1|X2=x2) ~ N(mu12,Sigma12) and marginal distribution P(X2) ~ N(mu2,Sigma2)
## returns mu12, Sigma12, mu2 and Sigma2
## 
##
## x2: numberic vector, an observation of X2
## mu1,Sigma1,A,b,Sigma21: as shown in above definition
## Precision1: inverse of Sigma1, either Precision1 or Sigma1 should be non-NULL
## Precision21: inverse of Sigma21, either Precision21 or Sigma21 should be non-NULL
##
## Examples:
## 
linearGaussian <- function(x2,mu1,Sigma1=NULL,Precision1=NULL,A,b,Sigma21=NULL,Precision21=NULL){
    ## again, only precision is needed
    if(is.null(Sigma1) & is.null(Precision1)) stop("Error in linearGaussian(): At least one of Sigma1 and Precision1 should be non-NULL")
    if(is.null(Precision1)) Precision1 <- solve(Sigma1)
    if(is.null(Sigma1)) Sigma1 <- solve(Precision1)
    if(is.null(Sigma21) & is.null(Precision21)) stop("Error in linearGaussian(): At least one of Sigma21 and Precision21 should be non-NULL")
    if(is.null(Precision21)) Precision21 <- solve(Sigma21)

    Precision12 <- Precision1 + t(A)%*%Precision21%*%A
    Sigma12 <- solve(Precision12)
    mu12 <- Sigma12 %*% (t(A)%*%Precision21%*%(x2-b) + Precision1%*%mu1)

    ## mu2 and Sigma2 can be easily calculated, so omit the output here.
    list(
        mu12=mu12,
        Sigma12=Sigma12## ,
        ## mu2 = A%*%mu1 + b,
        ## Sigma2 = A%*%Sigma1%*%t(A)+Sigma21
    )
}

## calculate the weighted scatter matrix
## copied from myTools.r
scatterMatrix <- function(x,w=1,mu=NULL){
    x <- as.matrix(x)
    if(length(w)==1L) w <- rep(w,nrow(x))
    if(!is.null(mu)){
        if(length(mu)!=ncol(x)) stop("Error in scatterMatrix(): mu must be the same length of observation dimensions")
        x <- sweep(x,2,mu,`-`)
    }
    t(w*x)%*%x
}

## function wrapper for setting a weakly informed (empirical) prior
## copied from myTools.r
priorNormal <- function(x){
    x <- as.matrix(x)
    D <- ncol(x)
    list(
        m0=colMeans(x),
        k0=0.001,                       #when k=0 the MAP of mu is the MLE of mu, but k=0 will result in a Inf in dNIW(), so use a small positive number instead of 0
        v0=ncol(x)+2,
        S0=diag(diag(scatterMatrix(x=x,w=1,mu=colMeans(x)))/nrow(x),nrow = D,ncol = D))
}

## log(sum(exp(l)))
## X: a logged matrix, nrow(X) = number of samples, ncol(X) = data dimension.
## return the log-sum-exp of each row of X
logsumexp <- function(X){
    X <- as.matrix(X)
    a <- apply(X,1,max)
    log(rowSums(exp(X-a)))+a
}

## log Multivariate Gamma Function
## Reference: https://en.wikipedia.org/wiki/Multivariate_gamma_function
## 
## a: numeric vector
## p: integer, the dimension
## 
lmvgamma <- function(a,p){
    sapply(a,function(ai){
        p*(p-1)/4*log(pi)+ sum(lgamma(ai+(1-(1:p))/2))
    },simplify = TRUE)
}

## log P(x_t | r_t, x^{(r_t)})
## return 1. the posterior parameters for
##     univariate normal: m and S
##     multivariate normal: m and S
##     ...
## and 2. the prediction probability
##     log P(x_t | r_t, x^{(r_t)})
## 
## there are two ways of using this function:
##     1. when x is not NULL, will generate the posterior parameters for each accumulation of x, and the probability of observing each x
##     2. prePos and xnew are not NULL and x is NULL, will generate the posterior parameter for data accumulated til xnew, and the probability of observing xnew
## x: numeric matrix, or an object that can be coerced to a matrix. Observation sequence, each row of x is a TRANSPOSE of an observation
## xnew: numeric vector, or an object that can be coerced to a vector. a new observations (only one observation)
## prePos: list, posterior parameters in previous time point
## prior: prior parameters, if model = "univariate-normal" or "multivariate-gaussian", prior should be a named list: list(m0,k0,v0,S0)
##                          if model = "poisson" or "exponential", prior should be a be a named list: list(scale)
##                          if model = "gamma" or "weibull", prior should be a named list:  list(shape,scale)
##                          if model = ...
##
## Example
## x <- rnorm(1000,mean = 10,sd = 5)
## out <- bcpPrediction(x=x,model = "univariate-normal",prior = list(m0=0,k0=1.2,v0=1.1,S0=2))
## out2 <- bcpPrediction(xnew = x[98],prePos = list(m=out$m[[97]],S=out$S[[97]],logS=out$logS[[97]]),model = "univariate-normal",prior = list(m0=0,k0=1.2,v0=1.1,S0=2))
## identical(out$m[[98]],out2$m)
## identical(out$S[[98]],out2$S)
## identical(out$p[[98]],out2$p)
bcpPrediction <- function(x,model=c("univariate-gaussian","multivariate-gaussian"),prior=list(m0=NULL,k0=NULL,v0=NULL,S0=NULL)){
    model <- match.arg(model)
    ## generate the posterior parameters for all data
    x <- as.matrix(x)
    nPars <- 1L:nrow(x)             #number of cache parameters per time point
    D <- ncol(x)                    #dimension
    if(model=="univariate-gaussian"){
        m0 <- prior$m0
        k0 <- prior$k0
        v0 <- prior$v0
        S0 <- prior$S0
        out <- list(m=lapply(nPars,FUN = function(l){rep(0,l)}),
                    S=lapply(nPars,FUN = function(l){rep(0,l)}),
                    logS=lapply(nPars,FUN = function(l){rep(0,l)}),
                    p=lapply(nPars,FUN = function(l){rep(0,l)}))
        k <- 1:nrow(x)+k0             #place holder for k
        v <- 1:nrow(x)+v0             #place holder for v
        out$m[[1]] <- (m0*k0 + x[1])/k[1]
        out$S[[1]] <- S0+ x[1]^2+ k0*(m0^2) - k[1]*(out$m[[1]]^2)
        out$logS[[1]] <- log(out$S[[1]])                                  #temporary variable, log(S) is used in calculating the predicted probability
        logs0 <- log(S0)                                                  #constant number log(S0)
        d2logpi <- -D/2 * log(pi)                                         #constant number -D/2 * log(pi)
        logkii <- D/2 * log(c(k0,k[1L:(nrow(x)-1L)])/k)                   #constant vector D/2 * log(k_{i-1}/k_i)
        vi2 <- c(v0,v[1L:(nrow(x)-1L)])/2                                 #constant vector v_{i-1}/2
        v2 <- v/2                                                         #constant vector v_i/2
        loggammavii <- lgamma(v2) -lgamma(vi2)                            #constant vector lmvGamma(v_i/2) - lmvGamma(v_{i-1}/2)
        out$p[[1]] <- d2logpi + logkii[1] + vi2[1] * logs0 - v2[1]*out$logS[[1]] + loggammavii[1]
        if(nrow(x)>=2){
            for(i in nPars[-1L]){
                out$m[[i]] <- c((m0*k0 + x[i])/k[1],
                (out$m[[i-1]]*k[1L:(i-1L)] + x[i]) / k[2L:i])
                out$S[[i]] <- c(S0+ x[i]^2+ k0*(m0^2) - k[1]*(out$m[[i]][1]^2),
                                out$S[[i-1]]+ x[i]^2+ k[1L:(i-1L)]*(out$m[[i-1]]^2) - k[2L:i] * (out$m[[i]][-1]^2))
                out$logS[[i]] <- log(out$S[[i]])
                ## -D/2 log(pi) + D/2 log(k_{i-1}/k_i) + v_{i-1}/2 log(|S_{i-1}|) - v_i/2 log(|S_i|) + lmvGamma(v_i/2) - lmvGamma(v_{i-1}/2)
                out$p[[i]] <- d2logpi + logkii[1L:i] + vi2[1L:i] * c(logs0,out$logS[[i-1]]) - v2[1L:i]*out$logS[[i]] + loggammavii[1L:i]
            }
        }
    }else if(model=="multivariate-gaussian"){
        ## out <- list(m=lapply(nPars,FUN = function(l){matrix(0,l,D)}),
        ##             S=lapply(nPars,FUN = function(l){array(0,dim = c(l,D,D))}))
        stop("multivariate-gaussian observation not supported yet")
    }

    return(out)
}

## loglikelihood for Left Truncated and Right Censored data
## return:
##  log P(q+1|q) = log S(q+1) - log S(q)  and
##  log P(0|q) = log f(q+1) - log S(q)
##  as an nx2 matrix, n is the number of inputs, the first column is log P(q+1|q), second column is log P(0|q)
## q: numeric vector
## model.pars: parameters for the selected model
##             weibull: shape parameter and scale parameter()
##
## Example
## logLTRC(1:10,shape = 2,scale = 10)
logLTRC <- function(q,cpt_model=c("weibull","gamma","exponential"),shape=NULL,scale=NULL){
    cpt_model <- match.arg(cpt_model)
    if(cpt_model=="weibull"){
        ## lower.tail=FALSE is the survival function

        ## ## when assuming p(0|i) is dead at time t+1
        ## matrix(c(pweibull(q+1,shape=shape,scale=scale,lower.tail = FALSE,log.p = TRUE),dweibull(q+1,shape=shape,scale=scale,log = TRUE))-
        ##        pweibull(q,shape=shape,scale=scale,lower.tail = FALSE,log.p = TRUE),ncol = 2)

        ## ## when assuming p(0|i) is dead at time (t,t+1]
        out <- matrix(0,nrow = length(q),ncol = 2)
        out[,1] <- pweibull(q+1,shape=shape,scale=scale,lower.tail = FALSE,log.p = TRUE)-pweibull(q,shape=shape,scale=scale,lower.tail = FALSE,log.p = TRUE)
        out[,2] <- log(-expm1(out[,1])) #IMPORTANT if use log(1-exp(out[,1])) will result in log(0) in some of the parameter settings, to avoid underflow, use expm1(x) to replace exp(x)-1
        out
    }else if(cpt_model=="gamma"){
       stop("In LTRC(): ",cpt_model," not supported yet!")
    }else if(cpt_model=="exponential"){
        stop("In LTRC(): ",cpt_model," not supported yet!")
    }else{
        stop("In LTRC(): ",cpt_model," not supported yet!")
    }
}

## generate random shape and scale sample, for the change point distribution
## when cpt_model = "weibull"
##      gamma distribution with prior parameter shape and scale for inverse-scale parameter
##      uniform distribution with prior start and end for shape parameter
## other cpt_model are not yet supported
rcpt <- function(cpt_model=c("weibull","gamma","exponential"),cpt_prior=list()){
    cpt_model <- match.arg(cpt_model)
    if(cpt_model=="weibull"){
        shape=runif(1,cpt_prior$start,cpt_prior$end)
        scale = 1/rgamma(1,cpt_prior$shape,cpt_prior$scale) #since gamma is for inverse-scale, now I need to inverse it back
        c(shape,scale)
    }else{
        stop("In rcpt(): ",cpt_model," not supported yet!")
    }
}

## constructor of a "bcp" object, used in offline learning
## use S3 with environment to achieve "pass by reference"
##
## x: numeric matrix, or an object that can be coerced to a matrix. Observation sequence, each row of x is a TRANSPOSE of an observation
## breaks: the ending index of each segments. breaks = NULL is equivalent to breaks = nrow(x)
## cpt_model: name of the change point distribution
## obs_model: name of the observation distribution
## obs_prior: prior parameters for observation model
##                  if obs_model = "univariate-normal" or "multivariate-gaussian", prior should be a named list: list(m0,k0,v0,S0)
##                  if obs_model = "poisson" or "exponential", prior should be a be a named list: list(scale)
##                  if obs_model = "gamma" or "weibull", prior should be a named list:  list(shape,scale)
##                  if obs_model = ...
## shape: initial shape parameter for change point distribution, no need to specify
## scale: initial scale parameter for change point distribution, no need to specify
bcp <- function(x=NULL,breaks=NULL,cpt_model=c("weibull","gamma","exponential"),obs_model=c("univariate-gaussian","multivariate-gaussian","multinomial","poisson","exponential","gamma","linear"),obs_prior=list(),cpt_prior=list(),shape=NULL,scale=NULL){

    cpt_model <- match.arg(cpt_model)
    obs_model <- match.arg(obs_model)

    object <- new.env(parent=globalenv())

    x <- as.matrix(x)
    object$x <- x

    if(is.null(breaks)){
        object$breaks <- c(0,nrow(x))
        object$Nsegs <- 1               #number of segments
        object$segLengths <- nrow(x)    #length of each segment
    }else{
        if(any(breaks<=0 | breaks > nrow(x)) | is.unsorted(breaks)) stop("Error in bcp(): values in 'breaks' has to be ordered and range between [1,nrow(x))")
        object$breaks <- c(0,breaks)
        object$Nsegs <- length(object$breaks)-1
        object$segLengths <- diff(object$breaks)
    }

    object$nPars <- numeric(0)             #number of unique run-lengths per time point
    for(i in 1:object$Nsegs)
        object$nPars <- c(object$nPars,1:object$segLengths[i])
    object$D <- ncol(x)                    #dimension
    object$T <- nrow(x)

    object$maxRun <- max(diff(object$breaks)) #maximum run length

    if(is.null(shape) & is.null(scale)){
        ss <- rcpt(cpt_model = cpt_model,cpt_prior = cpt_prior)
        object$shape <- ss[1]
        object$scale <- ss[2]
    }else if(!is.null(shape) & !is.null(scale)){
        object$shape <- shape
        object$scale <- scale
    }else{
        stop("Error in bcp(): shape and scale should both be NULL of non-NULL")
    }
    object$cpt_model <- cpt_model
    object$obs_model <- obs_model
    object$obs_prior <- obs_prior
    object$cpt_prior <- cpt_prior

    ## generate lower/upper search bound for L-BFGS-B algorithm
    if(cpt_model=="weibull"){
        object$lower <- c(cpt_prior$start,1e-7)
        object$upper <- c(cpt_prior$end,qgamma(0.99999999,shape = cpt_prior$shape,scale=cpt_prior$scale))
        ## object$upper <- c(cpt_prior$end,500)
    }else{
        stop("In bcp(): ",cpt_model," not supported yet!")
    }

    class(object) <- 'bcp'
    
    return(object)

}

## offline forward filtering
## x: numeric matrix, or an object that can be coerced to a matrix. Observation sequence, each row of x is a TRANSPOSE of an observation
## cpt_model: name of the change point distribution
## shape: shape parameter for change point distribution
## scale: scale parameter for change point distribution
## obs_model: name of the observation distribution
## obs_prior: prior parameters for observation model
##                  if obs_model = "univariate-normal" or "multivariate-gaussian", prior should be a named list: list(m0,k0,v0,S0)
##                  if obs_model = "poisson" or "exponential", prior should be a be a named list: list(scale)
##                  if obs_model = "gamma" or "weibull", prior should be a named list:  list(shape,scale)
##                  if obs_model = ...
bcpFiltering <- function(bcpObj,...) UseMethod("bcpFiltering")
bcpFiltering.bcp <- function(bcpObj){
    ## if(nrow(bcpObj$x)<2) stop("Error in bcpFiltering(): there should be at least 2 observations")
    
    bcpObj$a <- lapply(bcpObj$nPars,FUN = function(l){rep(0,l)}) #place holder for filtered result
    ## generate the transition matrix
    ## first column: P(q+1|q), second column: P(0|q)
    bcpObj$transition <- logLTRC(0L:(bcpObj$maxRun-1L),cpt_model = bcpObj$cpt_model,shape = bcpObj$shape,scale=bcpObj$scale)
    bcpObj$transition <- exp(bcpObj$transition)
    ## plot(transition[,2],type = "l") #increasing hazard
    
    ## generate the prediction probability: log P(x_t | r_t, x^{(r_t)})
    ## ### if(!exists("pxr",envir = bcpObj,inherits = FALSE)){ #must guarantee inherites=FALSE, or the namespace of parent frame will also be searched!!!!!!!!!!!!!
    bcpObj$pxr <- list()
    for(s in 1:bcpObj$Nsegs){
        bcpObj$pxr <- c(bcpObj$pxr,bcpPrediction(x=bcpObj$x[(bcpObj$breaks[s]+1):bcpObj$breaks[s+1],,drop=FALSE],model = bcpObj$obs_model,prior = bcpObj$obs_prior)$p) #only log probability is needed, so append a "$p" at the end
    }
    bcpObj$pxr_exp <- lapply(bcpObj$pxr,exp) #not log version
    
    ## debug code
    ## if(any(sapply(bcpObj$pxr,function(l){
    ##     any(is.nan(l)|is.na(l)|is.infinite(l))
    ## }))) stop("xxxx")
    
    tmp <- 0                            #temporary variable
    LL <- 0
    start <- 0                          #starting location of each segment
    for(s in 1:bcpObj$Nsegs){
        start <- bcpObj$breaks[s]
        bcpObj$a[[start+1]] <- 1
        LL <- LL+ bcpObj$pxr[[start+1]]               #observed data log likelihood
        for(i in 2:bcpObj$segLengths[s]){
            bcpObj$a[[start+i]] <- log(c(sum(bcpObj$transition[1:(i-1),2]*bcpObj$a[[start+i-1]]),bcpObj$transition[1:(i-1),1]*bcpObj$a[[start+i-1]]))+
                bcpObj$pxr[[start+i]]
            tmp <- logsumexp(matrix(bcpObj$a[[start+i]],nrow = 1)) #normalizing constant of the filtering process, which happens to be log P(x_t|x_1:t-1)
            LL <- LL+tmp
            bcpObj$a[[start+i]] <- bcpObj$a[[start+i]]-tmp
            bcpObj$a[[start+i]] <- exp(bcpObj$a[[start+i]])
        }
    }

    ## add prior likelihood
    if(bcpObj$cpt_model=="weibull"){
        logprior <- dunif(bcpObj$shape,min = bcpObj$cpt_prior$start,max = bcpObj$cpt_prior$end,log = TRUE)+
            dgamma(bcpObj$scale,shape = bcpObj$cpt_prior$shape,scale = bcpObj$cpt_prior$scale,log = TRUE)
    }else{
        stop("In bcpEM.bcp().obj(): ",bcpObj$cpt_model," not supported yet!")
    }
    
    LL <- LL+logprior
    
    ## invisible(gc())
    
    ## append resuts to bcp object
    bcpObj$LL <- LL
}

## offline forward filtering and backward smoothing
bcpSmoothing <- function(bcpObj,...) UseMethod("bcpSmoothing")
bcpSmoothing.bcp <- function(bcpObj){
    ## ----------------first part same as bcpFiltering----------------
    ## ---forward filtering---
    bcpFiltering(bcpObj)
    
    ## ---backward smoother---
    bcpObj$b <- lapply(bcpObj$nPars,FUN = function(l){rep(0,l)}) #place holder for backward result
    bcpObj$g <- lapply(bcpObj$nPars,FUN = function(l){rep(0,l)}) #place holder for smoothed result
    ## invisible(gc)
    start <- 0                          #starting location of each segment
    for(s in 1:bcpObj$Nsegs){
        start <- bcpObj$breaks[s]
        ## end <- bcpObj$breaks[s+1]
        bcpObj$b[[bcpObj$breaks[s+1]]] <- rep(1,bcpObj$segLengths[s])
        bcpObj$g[[bcpObj$breaks[s+1]]] <- bcpObj$a[[bcpObj$breaks[s+1]]]
        for(i in (bcpObj$segLengths[s]-1):1){
            ## transition: first column: P(q+1|q), second column: P(0|q)
            bcpObj$b[[start+i]] <- bcpObj$transition[1:i,1] * bcpObj$pxr_exp[[start+i+1]][-1] * bcpObj$b[[start+i+1]][-1]+
                bcpObj$transition[1:i,2] * bcpObj$pxr_exp[[start+i+1]][1]  * bcpObj$b[[start+i+1]][1]
            bcpObj$b[[start+i]] <- bcpObj$b[[start+i]]/sum(bcpObj$b[[start+i]])    #normalize to make sure it won't underflow
            bcpObj$g[[start+i]] <- bcpObj$b[[start+i]]*bcpObj$a[[start+i]]
            bcpObj$g[[start+i]] <- bcpObj$g[[start+i]]/sum(bcpObj$g[[start+i]])
        }
    }
}

## Two slice likelihoods, and filtered and smoothed
bcpTwoslice <- function(bcpObj,...) UseMethod("bcpTwoslice")
bcpTwoslice.bcp <- function(bcpObj){
    bcpSmoothing(bcpObj)

    ## can be see as a aparse version of the twoslice matrix
    ## because t in 1:T-1, so run length are in 0:T-2, i.e. T-1 unique different run lengths
    bcpObj$pii <- matrix(0,nrow = bcpObj$maxRun-1,ncol=2) #place holder for transition counts
    ## invisible(gc)
    ## place holder
    ## tmp <- matrix(0,nrow = bcpObj$maxRun-1,ncol=2)
    ## for(i in 1:(bcpObj$T-1)){
    ##     ## tansition: first column: P(q+1|q), second column: P(0|q)
    ##     tmp[1:i,] <- bcpObj$transition[1:i,] * bcpObj$a[[i]]
    ##     tmp[1:i,1] <- tmp[1:i,1] * bcpObj$b[[i+1]][-1] * bcpObj$pxr_exp[[i+1]][-1]
    ##     tmp[1:i,2] <- tmp[1:i,2] * bcpObj$b[[i+1]][1] * bcpObj$pxr_exp[[i+1]][1]
    ##     tmp[1:i,] <- tmp[1:i,]/sum(tmp[1:i,]) #normalize
    ##     bcpObj$pii[1:i,] <- bcpObj$pii[1:i,] + tmp[1:i,]     #the transition counts
    ## }

    ## system.time(for(s in 1:bcpObj$Nsegs){
    ##     start <- bcpObj$breaks[s]
    ##     for(i in 1:(bcpObj$segLengths[s]-1)){
    ##         ## tansition: first column: P(q+1|q), second column: P(0|q)
    ##         tmp[1:i,] <- bcpObj$transition[1:i,] * bcpObj$a[[start+i]]
    ##         tmp[1:i,1] <- tmp[1:i,1] * bcpObj$b[[start+i+1]][-1] * bcpObj$pxr_exp[[start+i+1]][-1]
    ##         tmp[1:i,2] <- tmp[1:i,2] * bcpObj$b[[start+i+1]][1] * bcpObj$pxr_exp[[start+i+1]][1]
    ##         tmp[1:i,] <- tmp[1:i,]/sum(tmp[1:i,]) #normalize
    ##         bcpObj$pii[1:i,] <- bcpObj$pii[1:i,] + tmp[1:i,]     #the transition counts
    ##         tmp[] <- 0
    ##     }
    ##             })
    for(s in 1:bcpObj$Nsegs){
        start <- bcpObj$breaks[s]
        for(i in 1:(bcpObj$segLengths[s]-1)){
            ## tansition: first column: P(q+1|q), second column: P(0|q)
            tmp <- bcpObj$transition[1:i,,drop=FALSE] * bcpObj$a[[start+i]]
            tmp[,1] <- tmp[,1] * bcpObj$b[[start+i+1]][-1] * bcpObj$pxr_exp[[start+i+1]][-1]
            tmp[,2] <- tmp[,2] * bcpObj$b[[start+i+1]][1] * bcpObj$pxr_exp[[start+i+1]][1]
            tmp <- tmp/sum(tmp) #normalize
            bcpObj$pii[1:i,] <- bcpObj$pii[1:i,] + tmp     #the transition counts
        }
    }
    
}


## Offline learning with exact method
## use generalized EM algorithm
## 
## nstart: number of random restarts
## maxit: maximum number of EM iterations for each (re)start
## deps: converge if LL-previousLL<deps
## 
bcpEM <- function(bcpObj,maxit=100,deps=1e-4,nstart=10L) UseMethod("bcpEM")
bcpEM.bcp <- function(bcpObj,maxit=100,deps=1e-4,nstart=10L){
    obj <- function(ss){
        if(bcpObj$cpt_model=="weibull"){
            logprior <- dunif(ss[1],min = bcpObj$cpt_prior$start,max = bcpObj$cpt_prior$end,log = TRUE)+
                dgamma(ss[2],shape = bcpObj$cpt_prior$shape,scale = bcpObj$cpt_prior$scale,log = TRUE)
        }else{
            stop("In bcpEM.bcp().obj(): ",bcpObj$cpt_model," not supported yet!")
        }
        ## optim is to minimize, so add an minus sign
        -sum(bcpObj$pii*logLTRC(0L:(bcpObj$maxRun-2L),cpt_model = bcpObj$cpt_model,shape = ss[1],scale=ss[2]))-logprior
    }

    maxLL <- -Inf                       #overall maximum likelihood
    LLs <- rep(0,nstart)
    bcpObj$MAP <- rep(0,0)                     #place holder for MAP shape and scale
    for(start in 1L:nstart){
        ## reset starting values since the second round(first round aleady has a starting value)
        if(start>=2){
            ss <- rcpt(cpt_model = bcpObj$cpt_model,cpt_prior = bcpObj$cpt_prior)
            bcpObj$shape <- ss[1]
            bcpObj$scale <- ss[2]
        }

        ## main EM loop start --------------------------------
        it <- 0
        maxLLit <- -Inf                 #current maximum likelihood/postori
        while(it <= maxit){
            ## 1. E step-------------------------------
            ## filter, smoother and twoslice
            bcpTwoslice(bcpObj)
            ## if(bcpObj$LL-maxLLit< -deps){
            if(bcpObj$LL < maxLLit){
                cat("diff",bcpObj$LL-maxLLit,"\n")
                stop("decreasing likelihood, LL: ",bcpObj$LL) #some numerical error is irreducible
            }
            else if(bcpObj$LL-maxLLit<deps) break 
            else maxLLit <- bcpObj$LL
            ## 2. M step------------------------------
            ## objective
            newSS <- optim(par = c(bcpObj$shape,bcpObj$scale),fn = obj,method = "L-BFGS-B",lower = bcpObj$lower,upper = bcpObj$upper)$par
            cat("iteration:",it,"  LL:",bcpObj$LL,"\n")
            cat("shape:",bcpObj$shape,"  scale:",bcpObj$scale,"\n")
            bcpObj$shape <- newSS[1]
            bcpObj$scale <- newSS[2]
            it <- it+1
        }
        ## main EM loop end --------------------------------
        
        if(maxLLit > maxLL){
            maxLL <- maxLLit
            bcpObj$MAP <- c(bcpObj$shape,bcpObj$scale) #get the best MAP estimate
        }
    }
    
}

bcpMCMC <- function(bcpObj,burnin=100,nSample=5000) UseMethod("bcpMCMC")
bcpMCMC.bcp <- function(bcpObj,burnin=100,nSample=5000){
}

## constructor of a "bcpo" object, used in online filtering
## use S3 with environment to achieve "pass by reference"
##
## cpt_model: name of the change point distribution
## shape: shape parameter for change point distribution
## scale: scale parameter for change point distribution
## obs_model: name of the observation distribution
## obs_prior: prior parameters for observation model
##                  if obs_model = "univariate-normal" or "multivariate-gaussian", prior should be a named list: list(m0,k0,v0,S0)
##                  if obs_model = "poisson" or "exponential", prior should be a be a named list: list(scale)
##                  if obs_model = "gamma" or "weibull", prior should be a named list:  list(shape,scale)
##                  if obs_model = ...
## l: inference lag, bcpOline will perform online filtering when when l=0, smoothing when l = Inf, fixed lag smoothing when 0<l<Inf. Default 0, filtering.
## 
bcpo <- function(shape=NULL,scale=NULL,cpt_model=c("weibull","gamma","exponential"),obs_model=c("univariate-gaussian","multivariate-gaussian","multinomial","poisson","exponential","gamma","linear"),obs_prior=list(),cpt_prior=list(),l=0){

    cpt_model <- match.arg(cpt_model)
    obs_model <- match.arg(obs_model)

    object <- new.env(parent=globalenv())

    if(is.null(shape) & is.null(scale)){
        ss <- rcpt(cpt_model = cpt_model,cpt_prior = cpt_prior)
        object$shape <- ss[1]
        object$scale <- ss[2]
    }else if(!is.null(shape) & !is.null(scale)){
        object$shape <- shape
        object$scale <- scale
    }else{
        stop("Error in bcpo(): shape and scale should both be NULL of non-NULL")
    }
    object$cpt_model <- cpt_model
    object$obs_model <- obs_model
    object$obs_prior <- obs_prior
    object$cpt_prior <- cpt_prior

    object$x <- NULL
    object$T <- 0L

    ## place holder for filtered/smoothing result and evidence likelihood
    object$a <- list()
    object$pxr <- list()
    object$pxr_exp <- list()

    ## place holder for changepoint probability
    object$pCPT <- numeric()

    ## place holder for changepoint indicator
    object$isCPT <- numeric()


    ## most recent posteriors
    if(obs_model %in% c("univariate-gaussian","multivariate-gaussian")){
        object$prePost <- list(m=list(obs_prior$m0),
                               k=obs_prior$k0,
                               v=obs_prior$v0,
                               S=list(obs_prior$S0))
    }else{
        stop("In bcpo(): ",obs_model," not supported yet!")
    }

    ## generate lower/upper search bound for L-BFGS-B algorithm
    if(cpt_model=="weibull"){
        object$lower <- c(cpt_prior$start,1e-7)
        object$upper <- c(cpt_prior$end,qgamma(0.99999999,shape = cpt_prior$shape,scale=cpt_prior$scale))
        ## object$upper <- c(cpt_prior$end,500)
    }else{
        stop("In bcpo(): ",cpt_model," not supported yet!")
    }

    object$l <- as.integer(l)                       #inference lag

    ## tracker for most recent filter/smoother index
    ## bcpOnline will only calculate the filter status between filterTracker and T
    ## bcpOnline will only calculate the (fixed-lag)smoother status between smootherTracker and T-l
    object$filterTracker <- 0L
    object$smootherTracker <- 0L

    object$transition <- matrix(numeric(),ncol = 2)

    class(object) <- 'bcpo'
    
    return(object)

}

## attach new observation(s) to an bcpo object
## newObs: numeric matrix of at least 1 row, or an object that can be coerced to a matrix.
## bcpoObj: an bcpo object
bcpoAttach <- function(bcpoObj,newObs) UseMethod("bcpoAttach")
bcpoAttach.bcpo <- function(bcpoObj,newObs){
    newObs <- as.matrix(newObs)

    if(is.null(bcpoObj$x)){
        bcpoObj$x <- newObs
    }else{
        bcpoObj$x <- rbind(bcpoObj$x,newObs)
    }

    currentMaxRun <- bcpoObj$T+nrow(newObs)-1L
    ## transition matrix P(r_i | r_{i-1})
    bcpoObj$transition <- rbind(bcpoObj$transition,
                                exp(logLTRC(bcpoObj$T:currentMaxRun,cpt_model = bcpoObj$cpt_model,shape = bcpoObj$shape,scale=bcpoObj$scale)))
    bcpoObj$T <- currentMaxRun+1L

    ## prediction probability log P(x_t | r_t, x^{(r_t)})
    ## and update posteriors
    if(bcpoObj$obs_model %in% c("univariate-gaussian","multivariate-gaussian")){
        for(i in 1:nrow(newObs)){                   #for each obs
            m <- list()
            S <- list()
            px <- numeric(length(bcpoObj$prePost$k))
            for(j in 1:length(bcpoObj$prePost$k)){                          #for each prior combination
                tmp <- posteriorGaussian(x=newObs[i,,drop=FALSE],m0=bcpoObj$prePost$m[[j]],k0=bcpoObj$prePost$k[[j]],v0=bcpoObj$prePost$v[[j]],S0=bcpoObj$prePost$S[[j]])
                m[[j]] <- tmp$mN
                S[[j]] <- tmp$SN
                px[j] <- marginalLikelihoodGaussian_byPosterior(N=tmp$N,
                              mN = tmp$mN,kN = tmp$kN,vN = tmp$vN,SN = tmp$SN,
                              m0 = bcpoObj$prePost$m[[j]],k0 = bcpoObj$prePost$k[[j]],v0 = bcpoObj$prePost$v[[j]],S0=bcpoObj$prePost$S[[j]],LOG = TRUE)
            }
            bcpoObj$prePost$k <- c(bcpoObj$obs_prior$k0,bcpoObj$prePost$k+1)
            bcpoObj$prePost$v <- c(bcpoObj$obs_prior$v0,bcpoObj$prePost$v+1)
            bcpoObj$prePost$m <- c(list(bcpoObj$obs_prior$m0),m)
            bcpoObj$prePost$S <- c(list(bcpoObj$obs_prior$S0),S)
            bcpoObj$pxr <- c(bcpoObj$pxr,list(px))
            bcpoObj$pxr_exp <- c(bcpoObj$pxr_exp,list(exp(px)))
        }
    }
    
    ## invisible(gc())                                #have to run garbage colloctor, because there are lots of garbages
    


}

## Online change point detection
##
## bcpoObj: an object of class "bcpo"
## newObs: an new observation, must be a matrix of at least one row, or anything that be converted to a matrix.
##
## attach the probability of change point bcpoObj$pCPT and the change point indicator bcpoObj$isCPT of each time index to bcpoObj
bcpOnline <- function(bcpoObj,newObs) UseMethod("bcpOnline")
bcpOnline.bcpo <- function(bcpoObj,newObs){
    newObs <- as.matrix(newObs)
    bcpoAttach(bcpoObj,newObs)           #attach new observation
    ## online filter ----------------------------------------
    magneticHead <- 0
    while(bcpoObj$T > bcpoObj$filterTracker){
        magneticHead <- bcpoObj$filterTracker+1L
        if(magneticHead==1){
            bcpoObj$a[[1]] <- 1
            bcpoObj$filterTracker <- magneticHead
            bcpoObj$pCPT <- 1
            bcpoObj$isCPT <- TRUE
            next
        }
        bcpoObj$a[[magneticHead]] <- log(c(sum(bcpoObj$transition[1:(magneticHead-1),2]*bcpoObj$a[[magneticHead-1]]),bcpoObj$transition[1:(magneticHead-1),1]*bcpoObj$a[[magneticHead-1]])) + bcpoObj$pxr[[magneticHead]]
        bcpoObj$a[[magneticHead]] <- bcpoObj$a[[magneticHead]]-logsumexp(matrix(bcpoObj$a[[magneticHead]],nrow = 1)) #minus the logsumexp to normalize
        bcpoObj$a[[magneticHead]] <- exp(bcpoObj$a[[magneticHead]])
        bcpoObj$filterTracker <- magneticHead
        bcpoObj$pCPT <- c(bcpoObj$pCPT,bcpoObj$a[[magneticHead]][1]) #attach change point probability
        bcpoObj$isCPT <- c(bcpoObj$isCPT,which.max(bcpoObj$a[[magneticHead]])==1)
    }

    ## online fixed-lag smoother ----------------------------
    magneticHead <- 0
    magneticEnd <- 0
    while(bcpoObj$smootherTracker < (bcpoObj$T-bcpoObj$l) & bcpoObj$l > 0L){
        magneticHead <- bcpoObj$smootherTracker+1L
        magneticEnd <- bcpoObj$smootherTracker+1L+bcpoObj$l
        bPre <- rep(1,magneticEnd)       #the last backward info
        ## pass backward info
        for(i in (magneticEnd-1L):magneticHead){
            b <- bcpoObj$transition[1:i,1] * bcpoObj$pxr_exp[[i+1]][-1] * bPre[-1]+
                bcpoObj$transition[1:i,2] * bcpoObj$pxr_exp[[i+1]][1]  * bPre[1]
            b <- b/sum(b)    #normalize to make sure it won't underflow
            bPre <- b
        }
        ## merge with forward info
        bcpoObj$a[[magneticHead]] <- b*bcpoObj$a[[magneticHead]] #merge backward info into a
        bcpoObj$a[[magneticHead]] <- bcpoObj$a[[magneticHead]]/sum(bcpoObj$a[[magneticHead]])
        bcpoObj$smootherTracker <- magneticHead
        bcpoObj$pCPT[magneticHead] <- bcpoObj$a[[magneticHead]][1] #update the believe of the lagged change point probability
        bcpoObj$isCPT[magneticHead] <- which.max(bcpoObj$a[[magneticHead]])==1 #update the changepoint ideicator
    }

    ## no need to run gc(), temporary objects will be deleted automatically when function exit
    ## invisible(gc())                     #garbage clean, there are lots of temporary variables generated
    
}
