# likelihood function

# mapping edges to off-diagnal elements of precision matrix K
parameter.index.matrix <- function(a,b){
    if(b-a>=1){
        return(c(a+((a+1):b-1)*b,parameter.index.matrix(a+1,b)))
    }
    else return()
}
# generate precision matrix K
precision.matrix <- function(d,Index,Degree,Theta,Sigma2){
    offdiagbal <- parameter.index.matrix(1,d)
    K <- rep(0,d*d)
    K[parameter.index.matrix(1,d)[Index]] <- (-Theta)/Sigma2
    K <- matrix(K,d,d)
    K <- K+t(K)
    diag(K) <- Degree/Sigma2
    K
}
# loglikelihood
loglikelihood <- function(g,m,x,ThetaB,Sigma2B,B){
    ## m: (m1,m2,...,mB) length:B
    ## g: list(glist,ilist)  length:B, d and index are derrived form g
    ## ThetaB: (theta1,theta2,...,thetaB)
    ## Sigma2B: (sigma21,sigma22,...,sigma2B)
    ## x: observations, data.frame
    d <- ncol(x)
    x$M <- rep(1:B,times=m) #grouping
    #loglikelihood
    ll <- ddply(x,.(M),function(xm,theta,sigma2,ug,dimension){
        i <- xm$M[1]
        K <- precision.matrix(d=dimension,Index = ug$ilist[[i]],
                              Degree = pmax(degree(ug$glist[[i]]),1),
                              Theta = theta[i],Sigma2 = sigma2[i])
        x <- as.matrix(xm[,1:dimension])
        data.frame(loglh= sum(
                       -diag(x%*%K%*%t(x))/2
                       -dimension/2*log(2*pi)+1/2*log(det(K)))
                   )
    },theta=ThetaB,sigma2=Sigma2B,ug=g,dimension=d)
    sum(ll$loglh)
}


