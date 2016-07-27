# state variable priors

# 1.priors for graph and block 
# prior.g and prior.m are uniform distribution, thus ommited from computation
logprior.g <- function(){0}
logprior.m <- logprior.g

# 2.priors for precision matrix(hierarchical case)
# Talih 3.3.2
logprior.theta <- function(ThetaB,psi){
    ThetaB <- ThetaB*pi/2
    sum(log(pi/2)-log(cos(ThetaB))*2)-sum(
        rollapply(ThetaB,width = 2,FUN = function(b,phi,alpha,eta2){
            (tan(b[2])-phi*tan(b[1])-(1-phi)*alpha)^2/
                (2*eta2)
        },phi=psi[1],alpha=psi[2],eta2=psi[3],align = "right")
        )
}
logprior.sigma2 <- function(Sigma2B,psi){
    sum(-log(Sigma2B))-sum(
        rollapply(Sigma2B,width = 2,FUN = function(b,phi,alpha,eta2){
            (log(b[2])-phi*log(b[1])-(1-phi)*alpha)^2/
                (2*eta2)
        },phi=psi[1],alpha=psi[2],eta2=psi[3],align = "right")
        )
}

# 3.priors for hyper-parameters
# Talih 3.3.3
logprior.psitheta <- function(psi,std=1,dfr=1){
    ## phi: auto-regressive
    ## alpha: trend
    ## eta2: persistence
    phi <- psi[1];alpha <- psi[2];eta2 <- psi[3]
    log(phi+1)+log(dnorm(alpha,sd=std))+log(dchisq(eta2,df=dfr))
}
logprior.psisigma2 <- logprior.psitheta



if(FALSE){
prior.phi.sample <- function(n.sample=1000){
    M <- 4
    phi <- NULL
    i <- 1
    while(i<=n.sample){
        y <- runif(1,-1,1);u <- runif(1)
        if(M/2*u<=y+1){
            phi <- c(phi,y)
            i <- i+1}
    }
    phi
}
prior.alpha.sample <- function(n.sample=1000,std=1){
    rnorm(n=n.sample,sd=std)
}
prior.eta2.sample <- function(n.sample=1000,dfr=1){
    rchisq(n=n.sample,df = dfr)
}

# Talih 4
## Psi <- matrix(c(prior.alpha.sample(10000),prior.phi.sample(10000),prior.eta.sample(10000)),10000,3)
logprior.theta <- function(ThetaB,Psi){
    ## Psi: alpha,phi,eta
    ThetaB <- ThetaB*pi/2
    mean(
        exp(
            sum(log(pi/2*cos(ThetaB)^(-2)))-
            colSums(rollapply(ThetaB,P=Psi,width = 2,FUN = function(Data,P){
                (tan(Data[2])-P[,2]*tan(Data[1])-(1-P[,2])*P[,1])^2/
                    (2*P[,3])
            },align="left"))
        )
    )
}
logprior.sigma2 <- function(Sigma2B,Psi){
    ## Psi: alpha,phi,eta
    mean(
        exp(
            sum(log(Sigma2B))-
            colSums(rollapply(Sigma2B,P=Psi,width = 2,FUN = function(Data,P){
                (log(Data[2])-P[,2]*log(Data[1])-(1-P[,2])*P[,1])^2/ 
                    (2*P[,3])
            },align="left"))
        )
    )
}
}
