#reference: Talih. Structural learning with time varing components
library(gtools) #combinations
library(igraph)
library(plyr)   
library(mvtnorm) #multi-variate normal sample generator
library(zoo) #rollapply

source("ini.r")                         #graph generating functions
source("likelihood.r")                  #likelihood functions
source("prior.r")                       #prior functions
source("kernel.r")                      #kernels
source("perturb.r")                     #perturbing functions

# state variables:
## g: graph list
## m: block lengths
## ThetaB: theta in each block
## Sigma2B: sigma^2 in each block
## psi: hyper parameters for ThetaB and Sigma2B, including psitheta and psisigma2

# other objects:
## x: data, data.frame
## K: precision matrix
## d: dimension
## B: number of blocks
## complete: complete graph
## zeta: perturbing parameter

# function used to generate simulated data
generate.data <- function(d,b1,ThetaB1,Sigma2B1,psitheta,psisigma2,m){
    ## b1: edges of first block
    ## ThetaB1,Sigma2B1: theta and sigma2 of first block
    complete <- generate.ug(d,Plot = FALSE)
    B <- length(m)
    g <- initial.ug(complete$ug,lout = complete$lyt,b1,B)
    ThetaB <- atan(hierarchy(initial=tan(ThetaB1*pi/2),
                             psitheta,B))/pi*2
    Sigma2B <- exp(hierarchy(initial=log(Sigma2B1),
                             psisigma2,B))
    x <- NULL
    for(i in 1:B){
        K <- precision.matrix(d,g$ilist[[i]],
                              pmax(degree(g$glist[[i]]),1),
                              ThetaB[i],Sigma2B[i])
        x <- rbind(x,data.frame(rmvnorm(m[i],sigma = solve(K))))
    }
    list(x=x,complete=complete,g=g,ThetaB=ThetaB,Sigma2B=Sigma2B)
}

# metropolis hastings sampling
MAIN <- function(x,complete,g,m,B,ThetaB,Sigma2B,psitheta,psisigma2,n.sample){
    ## n.sample: sample amount
    out <- list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B)
    pb <- txtProgressBar(2,n.sample,style = 3)
    ## main loop
    for(i in 2:n.sample){
        ## g
        g.new <- perturb.g(complete,g,B)
        a <- exp(loglikelihood(g.new,m,x,ThetaB,Sigma2B,B)-
                 loglikelihood(g,m,x,ThetaB,Sigma2B,B)) #min is omitted
        if(runif(1)<=a){
            g <- g.new
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        }
        else
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        ## m
        m.new <- perturb.m(m)
        a <- exp(loglikelihood(g,m.new,x,ThetaB,Sigma2B,B)-
                 loglikelihood(g,m,x,ThetaB,Sigma2B,B))
        if(runif(1)<=a){
            m <- m.new
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        }
        else
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        ## ThetaB
        ThetaB.new <- perturb.theta(ThetaB,B,psitheta)
        b <- length(which(ThetaB.new==ThetaB))+1
        a <- exp(loglikelihood(g,m,x,ThetaB.new,Sigma2B,B)+
                 logprior.theta(ThetaB.new,psitheta)+
                 logkernel.theta(ThetaB[b])-
                 loglikelihood(g,m,x,ThetaB,Sigma2B,B)+
                 logprior.theta(ThetaB,psitheta)-
                 logkernel.theta(ThetaB.new[b])
                 )
        if(runif(1)<=a){
            ThetaB <- ThetaB.new
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        }
        else
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        ## Sigma2B
        Sigma2B.new <- perturb.sigma2(Sigma2B,B,psisigma2)
        b <- length(which(Sigma2B.new==Sigma2B))+1
        a <- exp(loglikelihood(g,m,x,ThetaB,Sigma2B.new,B)+
                 logprior.sigma2(Sigma2B.new,psisigma2)+
                 logkernel.sigma2(Sigma2B[b])-
                 loglikelihood(g,m,x,ThetaB,Sigma2B,B)-
                 logprior.sigma2(Sigma2B,psisigma2)-
                 logkernel.sigma2(Sigma2B.new[b])
                 )
        if(runif(1)<=a){
            Sigma2B <- Sigma2B.new
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        }
        else
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        ## if(FALSE){
        ## psitheta
        psitheta.new <- perturb.psitheta(psitheta)
        a <- exp(logprior.theta(ThetaB,psitheta.new)+
                 logprior.psitheta(psitheta.new)+
                 logkernel.psitheta(psitheta)-
                 logprior.theta(ThetaB,psitheta)-
                 logprior.psitheta(psitheta)-
                 logkernel.psitheta(psitheta.new)
                 )
        if(runif(1)<=a){
            psitheta <- psitheta.new
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        }
        else
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        ## psisigma2
        psisigma2.new <- perturb.psisigma2(psisigma2)
        a <- exp(logprior.sigma2(Sigma2B,psisigma2.new)+
                 logprior.psisigma2(psisigma2.new)+
                 logkernel.psisigma2(psisigma2)-
                 logprior.sigma2(Sigma2B,psisigma2)-
                 logprior.psisigma2(psisigma2)-
                 logkernel.psisigma2(psisigma2.new)
                 )
        if(runif(1)<=a){
            psisigma2 <- psisigma2.new
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        }
        else
            out <- c(out,
                     list(graph=g,blocks=m,theta=ThetaB,sigma2=Sigma2B))
        ## }
        ## -----------------------------------
        i <- i+4
        setTxtProgressBar(pb,i)
    }
    out
}

# test on simulated data---------------------------------------------
## 1.simulate 600 multivariate-normal samples
Data <- generate.data(5,3,0.7,3,c(0.75,4,0.25),c(0.75,1,0.25),c(200,100,300))
plot.g(Data$complete,Data$g,3)
## 2.generate initial state of markov chain
initial.state <- generate.data(d=5,b1=2,ThetaB1 = 0.8,Sigma2B1 = 4,
                               psitheta = c(0.8,2,0.5),
                               psisigma2 = c(0.8,2,0.5),
                               m=rep(150,4))
plot.g(initial.state$complete,initial.state$g,4)
## 3.generate 120000 sample
out <- MAIN(x=Data$x,complete=Data$complete,
            g=initial.state$g,m=rep(150,4),B=4,
            ThetaB=initial.state$ThetaB,Sigma2B=initial.state$Sigma2B,
            psitheta = c(0.8,2,0.5),
            psisigma2 = c(0.8,2,0.5),
            n.sample = 20000)

tmp <- ldply(out[(1:5000)*4-2],function(d){matrix(d,nrow = 1)})

save(out,file = "loop12000")

plot.g(Data$complete,initial.state$g,4)
plot.g(Data$complete,out[[20000-3]],4)
plot.g(Data$complete,Data$g,3)


out2 <- MAIN(x=Data$x,complete=Data$complete,
            g=out[[20000-3]],m=out[[20000-2]],B=4,
            ThetaB=out[[20000-1]],Sigma2B=out[[20000]],
            psitheta = c(0.8,2,0.5),
            psisigma2 = c(0.8,2,0.5),
            n.sample = 20000)

tmp <- ldply(out2[(1:5000)*4-2],function(d){matrix(d,nrow = 1)})
plot.g(Data$complete,out2[[20000-3]],4)
plot.g(Data$complete,Data$g,3)


out3 <- MAIN(x=Data$x,complete=Data$complete,
            g=out2[[20000-3]],m=out2[[20000-2]],B=4,
            ThetaB=out2[[20000-1]],Sigma2B=out2[[20000]],
            psitheta = c(0.8,2,0.5),
            psisigma2 = c(0.8,2,0.5),
            n.sample = 40000)

tmp <- ldply(out3[(1:10000)*4-2],function(d){matrix(d,nrow = 1)})
plot.g(Data$complete,out3[[40000-3]],4)
plot.g(Data$complete,Data$g,3)

par(mfcol = c(5,1),mar=rep(0,4))
for(i in 1:5 ) {
    plot(x=1:600,y=Data$x[,i],type = "l")
    abline(v=c(200,300),lty=3,col="red")
}

par(mfcol = c(1,3))
for(i in 1:3) plot(Data$g$glist[[i]],layout=Data$complete$lyt,main=i)

par(mfcol = c(1,3))
for(i in 3:4) plot(out2[[40000-3]]$glist[[i]],layout=Data$complete$lyt,main=i)


tmp <- ldply(out2[(1:10000)*4-2],function(d){matrix(d,nrow = 1)})
