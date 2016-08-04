#FNN: Fast Nearest Neighbor Search Algorithms and Applications
#install.packages("FNN")
library(FNN)
library(zoo)
library(compiler)
library(plyr)

#local reachability distance
#lrd=1/(sum(reach.dist)/k) 其中 reach.dist=max(dist,knn.dist)
#kni是从get.knn传出的list: kni <- get.knn(x)
#kni中有两个矩阵nn.index nn.dist 分别为k nearest neighbours 和 distance to k nearest neighbours
lrd <- function(x=c(),k=integer(),kni=list()){
    reach.dist <- pmax(as.matrix(dist(x)),matrix(kni$nn.dist[,k])) #自动按列循环匹配kni$nn.dist[,k]
    unlist(                                                        
        llply(.data = seq_along(x),.fun = function(.data,reach.dist,nn.index,k){
            k/sum(reach.dist[.data,nn.index[.data,]])
        },reach.dist=reach.dist,nn.index=kni$nn.index,k=k)
    )
}

#local outlier factor
lof <- function(x=c(),k=integer(),kni=list()){
    x.lrd <- lrd(x,k,kni)
    unlist(
        llply(.data = seq_along(x.lrd),.fun = function(.data,x.lrd,nn.index,k){
            sum(x.lrd[nn.index[.data,]])/x.lrd[.data]/k
        },x.lrd=x.lrd,nn.index=kni$nn.index,k=k)
    )
}

#k influence space
## reverse k nearest neighbor(RkNN): Given a point q, a reverse k nearest neighbor(RkNN) query retrieves all the data points that have q as one of their k nearest neighbors.
## kis = knn ∪ rknn
kis <- function(kni=list()){
    llply(.data =1:nrow(kni$nn.index),.fun = function(.data,nn.index){
        union(nn.index[.data,],as.vector(nn.index[nn.index[.data,],]))
    },nn.index=kni$nn.index)
}

#influenced outlierness
inflo <- function(x=c(),k=integer(),kni=list()){
    dens <- 1/kni$nn.dist[,k]           #density
    x.kis <- kis(kni)
    unlist(
        llply(.data =seq_along(x),.fun = function(.data,dens,x.kis,nn.index){
            current.kis <- x.kis[[.data]]
            sum(dens[current.kis])/length(current.kis)/dens[.data]
        },dens=dens,x.kis=x.kis,nn.index=kni$nn.index)
    )
}

switch.detector<-function(detectorname){
  cmpfun(
    switch(detectorname,
           lof=lof,
           inflo=inflo
    ))
}

#main
outlier.detection <- function(x=c(),k=10,detectorname="inflo",thres=10,Diff=TRUE,PLOT=TRUE){
    x.origin <- x
    if(Diff)
        x <- diff(x)
    kni <- get.knn(x,k=k)
    detector <- switch.detector(detectorname)
    x.dected <- detector(x,k,kni)
    outliers <- which(x.dected>thres)

    if(PLOT){
    par(mfcol = c(2,2))
    
    #origin line
    plot(x.origin,type = "l",ylab = "value",main = "Original Data")
    points(outliers,x.origin[outliers],pch="o",col="red")
    #diff line
    plot(x,type = "l",ylab="diff",main = "First Order Difference")
    points(outliers,x[outliers],pch="o",col="red")
    #boxplot
    y.xis <- rollapply(hist(x,plot = FALSE)$breaks,width = 2,by = 1,FUN = mean)
    x.xis <- hist(x,plot=FALSE)$density
    plot(x.xis,y.xis,type = "l",
         xlim=c(-range(x.xis)[2],range(x.xis)[2]),
         ylim=c(range(y.xis)[1]*1.1,range(y.xis)[2]*1.1),xlab="density",ylab="diff",main = "Boxplot of First Order Difference")
    lines(-x.xis,y.xis)
    points(rep(0,length(x)),x,pch="-")
    points(rep(0,length(outliers)),x[outliers],pch="o",col="red")
    #histogram
    hist(x,freq = FALSE,xlab="diff",main = "Density of First Order Difference")
    points(x[outliers],rep(0,length(outliers)),pch="*",col="red")
    ## dev.off()
    }
    return(outliers)
}



#例子-------------------------------------------------------------------------------

#加载数据
ram <- readLines("/tmp/memory",encoding = "UTF-8")
ram <- strsplit(ram,split = " ")
ram <- llply(.data = ram,.fun = function(.data){
    as.numeric(.data)
})

#检测outlier并画图
outlier.detection(ram[[1]])
outlier.detection(ram[[1]],detectorname="lof")
outlier.detection(ram[[2]])
outlier.detection(ram[[2]],detectorname="lof")
outlier.detection(ram[[3]])
outlier.detection(ram[[4]])
outlier.detection(ram[[5]])
outlier.detection(ram[[5]],detectorname="lof")
outlier.detection(ram[[6]])
outlier.detection(ram[[7]],detectorname = "lof",thres = 3)

outlier.detection(tmp$V1,thres = 3,detectorname="lof",Diff = FALSE)

plot(tmp$V1)
