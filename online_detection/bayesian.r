
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
