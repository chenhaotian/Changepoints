library(ggplot2)
library(plyr)
library(lubridate)


CUMSUM <- function(x){
    n <- length(x)
    res <- NULL
    for(i in seq_along(x)){
        res[i] <- sum(x[1:i])-i/n*sum(x)
    }
    res/sqrt(n)
}

#level--------------------------------------------
Cstar <- function(x,h=0.1,s=1,e=length(x)){
    y <- x[s:e] #仅看当前序列
    sigma <-sum((y-mean(y))^2)/length(y)
    n <- e-s+1
    upper <- n-round(n*h)
    lower <- round(n*h)
    cum <- CUMSUM(y)^2
    cum <- cum[-c(1:lower,upper:n)]
    res <- NULL
    for(i in (lower+1):(upper-1)){
        res[i-lower] <- i/n*(1-i/n)
    }
    C <- cum/sigma/na.omit(res)
    if(max(cum/sigma/res)<12.827)
        return()
        ## return("no change point under 0.05 level")
    else{
        changepoint <- lower+which.max(cum/sigma/res)-1+s
        print(paste("change point:",changepoint))
        res <- NULL
        res <- c(res,changepoint,
                 Cstar(x,s=changepoint+1,e=e),
                 Cstar(x,s=s,e=changepoint))
        return(res)
    }
}
plot.level <- function(y,cp){
    n <- length(y)
    cp <- c(1,cp,n)
    cp <- cp[order(cp)]
    g <- NULL
    for(i in 2:length(cp))
        g <- c(g,rep(i-1,cp[i]-cp[i-1]))
    g <- c(1,g)
    DATA <- data.frame(y,g)
    DATA <- ddply(.data = DATA,.variables = "g",.fun = function(.data){
        n <- nrow(.data)
        mu <- mean(.data$y)
        data.frame(y=.data$y,g=.data$g,mu=rep(mu,n))
    })
    DATA$x <- 1:n
    ggplot(DATA)+geom_path(aes(x=x,y=y),alpha=I(1/2))+geom_path(aes(x=x,y=mu))
}

#trend-------------------------------------------
Dstar <- function(x,h=0.1,s=1,e=length(x)){
    y <- x[s:e] #仅看当前序列
    n <- e-s+1
    mod <- lm(y~t,data = data.frame(y,t=1:n))
    rsd <- mod$residuals
    sigma <-sum(rsd^2)/(n-2)
    upper <- n-round(n*h)
    lower <- round(n*h)
    Tc <- NULL
    for(i in (lower+1):(upper-1)){
        Tc[i-lower] <- -(n^(-1/2))*sum(rsd[1:i])/sqrt(sigma)/
            sqrt(i/n*(1-i/n)*(1-3*i/n*(1-i/n)))
    }
    D <- Tc^2
    if(max(D)<14.692)
        return()
        ## return("no change point under 0.05 level")
    else{
        changepoint <- lower+which.max(D)-1+s
        print(paste("change point:",changepoint))
        res <- NULL
        res <- c(res,changepoint,
                 Dstar(x,s=changepoint+1,e=e),
                 Dstar(x,s=s,e=changepoint))
        return(res)
    }
}
plot.trend <- function(y,cp){
    n <- length(y)
    cp <- c(1,cp,n)
    cp <- cp[order(cp)]
    g <- NULL
    for(i in 2:length(cp))
        g <- c(g,rep(i-1,cp[i]-cp[i-1]))
    g <- c(1,g)
    DATA <- data.frame(y,g)
    DATA <- ddply(.data = DATA,.variables = "g",.fun = function(.data){
        n <- nrow(.data)
        mod <- lm(y~t,data = data.frame(y=.data$y,t=1:n))
        py <- mod$coefficients[1]+mod$coefficients[2]*(1:n)
        data.frame(y=.data$y,g=.data$g,py)
    })
    DATA$x <- 1:n
    ggplot(DATA)+geom_path(aes(x=x,y=y),alpha=I(1/2))+geom_path(aes(x=x,y=py))
}



## x <- c(rnorm(1000),rnorm(700,mean=1),rnorm(700,mean = 2),rnorm(700),rnorm(700,mean = 1))
## y <- x+seq(1:length(x))*0.01

## par(mfcol = c(2,1))
## plot(x,type = "l",main = "Raw Level")
## plot(y,type = "l",main = "Raw Trend")



## cp <- Cstar(x)
## plot.level(x,cp)

## cp <- Dstar(y)
## plot.trend(y,cp)
