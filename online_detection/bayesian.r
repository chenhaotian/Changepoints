
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
