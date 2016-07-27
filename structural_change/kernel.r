# proposal distributions of parameters and hyper-parameters

# 1.kernels for grahp and block
# they are symmetirc, thus kernel(new,old)==kernel(old,new)
logkernel.g <- function(){0}
logkernel.m <- logkernel.g

# 2.kernels for precision matrix
# Talih 4.3.1
logkernel.theta <- function(Theta.new,Theta.old=NULL){
    ## symmetirc on tan scale
    -2*log(cos(Theta.new*pi/2))
}
logkernel.sigma2 <- function(Sigma2.new,Sigma2.old=NULL){
    ## symmetric on log scale
    -log(Sigma2.new)
}

# 3.kernels for hyper-parameters
# Talih 4.3.2
logkernel.psitheta <- function(psi.new,psi.old=NULL){
    #phi is symmetric on tan scale
    #eta2 is symmetric on log scale
    #alpha is symmetric
    logkernel.theta(Theta.new=psi.new[1])+
        logkernel.sigma2(Sigma2.new=psi.new[3])
}
logkernel.psisigma2 <- logkernel.psitheta

