# Bayesian online change point detection and offline learning



+ [Introduction](#introduction)
+ [Usage](#usage)
  + [R version](#R version)
    + [R Installation](#R Installation)
    + [R Functions](#R Functions)
    + [R Examples](#R Examples)
  + [C++ version](#C++ version) (in development)
+ [Derivation](#derivation)
+ [Reference](#reference)
+ [License](#license)



## Introduction

A comprehensive extension and upgrade to the model proposed in `Bayesian Online Changepoint Detection(Ryan Prescott Adams and David J.C. MacKay 2007)`.

Adam and MacKay's  paper treat change point detection as an online filtering problem on an infinite state hidden Markov model, where the duration of the segments(a change point is defined as the beginning of a new segment) are assumed exponentially distributed(so that the hazard rate, or the probability of a change point occurring is constant w.r.t time). The model provides a simple and intuitive structure for change point problems but suffers from it's over simplicity and requires too many assumptions. This repository aims at extend the model's capability so as to make it more applicable to real life scenarios. The extensions are:

+ For online inference, add support for fixed-lag smoothing, so that users can trade off latency with accuracy.
+ For offline learning, add support for learning the segment duration distribution directly from historical data, so that users do not need to make up the duration distribution themselves. There are two learning algorithms applied in this repository, one based on exact inference, the other based on simulation. Details is explained in the following sections.
+ For segment duration distribution, add support for distributions with non-constant hazard rates, such as Weibull, Gamma, log-Normal ...
+ For evidence distribution, add support for distributions with various structures, such as (Bayesian)linear regression, multivariate Gaussian, multinomial ...

**<span style="color:green">Note:</span>**

For now only **univariate-Gaussian** observation and **Weibull** segment duration is supported. Others are in development.



## Usage

The code has 2 versions, **R** and **C++**. Both versions provide the same functionalities. R version is better for educational purpose, but runs relatively slower.

### R version

#### R Installation

Load the script directly from source:
```R
source("https://raw.githubusercontent.com/chenhaotian/Changepoints/master/R/BOCPDOL.r")
```
#### R Functions

This section lists the only the general descriptions. See the section [R Examples](#R Examples) for implementation notes.

| Name            | Description                                                  |
| --------------- | ------------------------------------------------------------ |
| **bcp**         | Create a Bayesian change point object for **offline learning**. <br/>**Usage**: <br/>    bcp(x = NULL, breaks = NULL,  cpt_model = c("weibull", "gamma", "exponential"), obs_model = c("univariate-gaussian", "multivariate-gaussian", "multinomial", "poisson", "exponential",  "gamma","linear"), obs_prior = list(), cpt_prior = list(), shape = NULL, scale = NULL,)<br/>**Arguments**:<br/>    x: numeric matrix, or an object that can be coerced to a matrix. Observation sequence, each row of x is a TRANSPOSE of an observation vector.<br/>    breaks: if 'x' is an `rbind` of multiple set of observations, 'breaks' will be the ending row number of each sets. For example if 'x' is an `rbind` of 3 days time series with 5,9 and 3 observations each day. then 'breaks' will be c(5,14,17). When 'breaks' = NULL, the model will assume all observations came from only dataset. Default NULL.<br/>    cpt_model: the segment duration distribution to be used, must be one of "weibull","gamma" and "exponential"<br/>    obs_model: the observation distribution to be used, must be one of "univariate-gaussian","multivariate-gaussian","multinomial","poisson","exponential","gamma" and "linear".<br/>    obs_prior: hyper parameters for the observation distribution. See Examples for details. <br/>    cpt_prior: hyper parameters for the segment duration distribution. See Examples for details.<br/>    shape, scale: the initial shape and scale parameter for the segment duration distribution. No need to specify, default NULL. <br/>**Value**:<br/>    returns an object of class 'bcp'. This object will be further used in the learning and offline smoothing processes. |
| **bcpEM**       | Get MAP estimates of the segment duration with (generalized) EM algorithm. <br/>**Usage**:<br/>    bcpEM(bcpObj, maxit = 100, deps = 1e-04, nstart = 10L)<br/>**Arguments**:<br/>    bcpObj, and object created by bcp().<br/>    maxit, number of maximum EM iterations.<br/>    deps, learning precision, the learning stops when the difference between the most recent two observed data likelihoods is smaller than 'deps'.<br/>    nstart, number of random restarts, set 'nstart' bigger to avoid local optimal, but will cost more time.<br/>**Value**:<br/>    The MAP estimates of shape and scale parameters will be stored in `bcpObj$MAP`. |
| **bcpMCMC**     | Sample from the posterior distribution of the segment durations with MCMC.<br/>**Usage**:<br/>    bcpMCMC(bcpObj, burnin = 100, nSample = 5000)<br/>**Arguments**:<br/>    bcpObj: an object created by bcp().<br/>    burnin: number of burn-in samples.<br>    nSamples: number of samples to draw after burn-in.<br/> **Value**:<br/>    The posterior samples of shape and scale will be stored in `bcpObj$postSamples`. |
| **bcpo**        | Create a Bayesian change point object for **online inference**. <br/>**Usage**: <br/>    bcpo(shape = NULL, scale = NULL, cpt_model = c("weibull", "gamma", "exponential"), obs_model = c("univariate-gaussian", "multivariate-gaussian", "multinomial", "poisson", "exponential",  "gamma", "linear"), obs_prior = list(), cpt_prior = list(), l = 0)<br/>**Arguments**:<br/>    shape, scale: the shape and scale of the segment duration distribution. shape and scale can be learned from bcpEM() or bcpMCMC(), or be specified manually.<br/>    cpt_model, obs_model, obs_prior, cpt_prior: same as the ones defined in bcp().<br/>    l: inference lag, bcpOnline will perform online filtering when l=0, fixed-lag smoothing when l>0.<br/>**Value**:<br/>    returns an object of class 'bcpo'. This object will be further used in the bcpOnline() function for online filtering and fixed-lag smoothing. |
| **bcpOnline**   | Bayesian change point online filtering and fixed-lag smoothing.<br/>**Usage**:<br/>    bcpOnline(bcpoObj, newObs)<br/>**Arguments**:<br/>    bcpoObj: an object created by bcpo().<br/>    newObs: new observations matrix, or an object that can be coerced to a matrix.<br/>**Value**:<br/>    The filtered/fixed-lag smoothed change point probability will be attached to `bcpoObj$pCPT`;<br>    The change point indicator of each time point will be attached to `bcpoObj$isCPT` |
| **priorNormal** | Pick an weekly informed Normal-Inverse-Wishart prior for Gaussian distributions empirically.<br/>**Usage**:<br/>    priorNormal(x)<br/>**Arguments**:<br/>    x: numeric matrix, or an object that can be coerced to a matrix.<br/>**Value**:<br/>    A list of length 4, representing the NIW parameters 'm', 'k', 'v' and 'S' respectively. |



#### R Examples

**Example 1:** Learn segment duration distribution with EM algorithm, then perform online inference.

Say we have observed a time series data in 10 days, there are 200, 300, 220, 250, 190, 290, 310, 320, 230, 220, 220, 240, 310, 230 and 242 observations in each day. We want to learn the segment duration distribution from this 10 series. 

After learning the segment duration distribution, we want to apply the learned model to some new data for online inference.

Generate some synthetic data:

```R
## generate historical data 'x', x is normally distributed, the segment duration is Weibull with shape=2 and scale = 50
set.seed(2)
segLengths <- c(200,300,220,250,190,290,310,320,230,220,220,240,310,230,242) #
BREAKS <- cumsum(segLengths)
x <- unlist(sapply(segLengths,function(n){
    head(unlist(sapply(round(rweibull(round(n/2),shape = 2,scale = 50)),function(l){rnorm(l,mean = rnorm(1,mean = 0,sd = 40),sd = rgamma(1,shape = 10,rate = 3))})),n)
}))
## generate new data 'xnew', with the same distributions as 'x'
set.seed(3)
newx <- unlist(sapply(round(rweibull(10,shape = 2,scale = 50)),function(l){
    rnorm(l,mean = rnorm(1,mean = 0,sd = 40),sd = rgamma(1,shape = 10,rate = 3))
}))
```

Offline learning from historical data:

```R
## for learning purpose. Create an 'bcp' object named "bcpObj"
## where we use univariate-Gaussian as the observation model, Weibull as the segment duration model.
## hyper parameters of univariate-Gaussian is calculated imperically using priorNormal(x), 
## hyper parameters of Weibull is assumed uniform(start=1,end=10) for shape, and gamma(shape=1,scale=100) for scale. In this model we assume the Weibull parameters are marginally independent of each other.
bcpObj <- bcp(x,breaks = BREAKS,cpt_model = "weibull",obs_model = "univariate-gaussian",obs_prior = priorNormal(x),cpt_prior = list(start=1,end=10,shape=1,scale=100))
## Learning Weibull(shape,scale) with generalized EM algorithm, random restart 2 times to avoid local optimal
bcpEM(bcpObj,nstart = 2)
```

Online filtering on new observations

```R
## perfom online filtering(by setting l=0)
bcpoObj1 <- bcpo(shape = bcpObj$MAP[1],scale = bcpObj$MAP[2],cpt_model = cpt_model,obs_model = obs_model,obs_prior = obs_prior,cpt_prior = list(start=1,end=10,shape=1,scale=100),l=0)
## use a for loop to mimic the arrival of data stream
for(i in seq_along(newx)) bcpOnline(bcpoObj=bcpoObj1,newObs=newx[i])
```

Online fixed-lag smoothing on new observations:

```R
## perfom online fixed-lag smoothing(by setting l>0, in this case let's use l=10 and l=40)
bcpoObj2 <- bcpo(shape = bcpObj$MAP[1],scale = bcpObj$MAP[2],cpt_model = cpt_model,obs_model = obs_model,obs_prior = obs_prior,cpt_prior = list(start=1,end=10,shape=1,scale=100),l=10L)
bcpoObj3 <- bcpo(shape = bcpObj$MAP[1],scale = bcpObj$MAP[2],cpt_model = cpt_model,obs_model = obs_model,obs_prior = obs_prior,cpt_prior = list(start=1,end=10,shape=1,scale=100),l=40L)
## use a for loop to mimic the arrival of data stream
for(i in seq_along(newx)){
    bcpOnline(bcpoObj=bcpoObj2,newObs=newx[i])
    bcpOnline(bcpoObj=bcpoObj3,newObs=newx[i])
}
```

Compare the filtered and fixed-lag smoothed results:

```R
par(mfcol = c(3,1))
plot(newx,type="l",main = "lag = 0 (equivalent to filtering)")
abline(v=which(bcpoObj1$isCPT),lty=2,col="red")
plot(newx,type="l",main = "lag = 10")
abline(v=which(bcpoObj2$isCPT),lty=2,col="red")
plot(newx,type="l",main = "lag = 40")
abline(v=which(bcpoObj3$isCPT),lty=2,col="red")
```

![](./notes_pictures/example1.png)



**Example 2:** Learn with MCMC, then perform online inference.

**in development ...**



### C++ version

**in development ...**



## Derivation

### Model

Denote $r_t$ denote the time since the last change point, the "run length". And use $x^{(r_t)}$ to indicate the set of observations associated with the run $r_t$. As $r$ may be zero, the set $x^{(r)}$ may be empty:
$$
\begin{equation*}
x^{(r_t)} = \begin{cases}
  \emptyset               & r_t = 0\text{, }t=1,...T \\
  x_{t-r_{t}} : x_{t-1}   & r_t > 0\text{, }t=1,..,T
\end{cases}
\end{equation*}
$$

The conditional probability distributions are:
$$
\begin{align*}
& P(r_1 = 0) = 1\\
& P(r_t | r_{t-1}) =
  \begin{cases}
    \frac{S_d(r_{t-1})-S_d(r_t)}{S_d(r_{t-1})\text{}}       & r_t = 0 \text{, } t \ge 2 \\
    \frac{S_d(r_t)}{S_d(r_{t-1})}                    & r_t = r_{t-1}+1 \text{, } t \ge 2
  \end{cases} \\
& \text{ where } S_d(r_t) = P(d_t \ge r_t | \theta_d)\\
& P(\theta_d | \eta_d) = \mathcal{H}_d(\eta_d) \\
& P(\theta_x | \eta_x) = \mathcal{H}_x(\eta_x) \\
& P(x_t | \theta_x) = \mathcal{E}(\theta_x) \\
\end{align*}
$$


**Adams&Mackay's paper interpret the change point detection problem with an infinite state hidden markov model. Hereby the definitions:**

+ **$r_t$ is the time length since last changepoint, $x_{t}^{(r)}$ denote the data sets associated with run $r_t$. As $r$ may be zero, $x^{(r)}$ may be empty.The hidden states are the run lengths $r$ of current observation, they evolve as the following way:**
**![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/online_detection/state_transition.png)**
+ **For each parition $\rho$, the data within it are i.i.d. from some probability distribution $P (x_t | \eta_\rho)$**
+ **The discrete a priori(Latin phrase, means 'from the earlier',contrast to a posteriori,'from the latter') probability distribution over the interval between changepoints is denoted as $P_{gap}(g)$, i.e. the probability of current run ending up at length $g$.(recall geometric distribution)**
+ **Predictive distribution $P(x_{t+1}|r_t,x_t^{(r)})$**
  **Posterior distribution $P(r_t|x_{1:t})$**
  **Joint distribution $P(r_t,x_{1:t})$**

**With the above definitions, our goal here is to <span style="color:darkred">recursively(use recusive methods for the purpose of online calculation)</span> compute the probability distribution of the length of the current “run”, or time since the last changepoint, denote this posterior distribution as $P(r_t|x_{1:t})$**
**According to bayesian therom,**
**$$P(r_t|x_{1:t})=\frac{P(r_t,x_{1:t})}{P(x_{1:t})}, \tag{1}$$**
**we can devide the calculation of the conditional distribution into two subproblems:**
**1. Calculate the joint distributioin $P(r_t,x_{1:t})$.**
**2. Calculate the normalizing constant $P(x_{1:t})$.**
**As for subproblem 2, $P(x_{1:t})$ can be easily calculated by integrate out $r_t$ in subproblem 1. i.e. $P(x_{1:t})=\sum_{r_t}P(r_t,x_{1:t})$. So basically we need only focus on subproblem 1, the derivation of $P(r_t,x_{1:t})$.**
**Keep in mind that to make the online calculation possible, the derived represntation of subproblem 1 must be recursive. i.e. when the previous information is known, one need to derive a relation between previous information and current estimation. Here in subproblem 1, the previous information is $P(r_{t-1}|x_{1:t-1})$, and the current estimation is $P(r_t|x_{1:t})$. In order to relate these two distributions, first we need to <span style="color:darkred">disintegrate</span>  $P(r_t|x_{1:t})$ by $r_{t-1}$, and then use bayes rule to extract $P(r_{t-1}|x_{1:t-1})$:**
$$
\begin{align}
P(r_t,x_{1:t}) & = \sum_{r_{t-1}} P(r_t,r_{t-1},x_{1:t})\\
= & \sum_{r_{t-1}} P(r_t,x_t|r_{t-1},x_{1:t-1})P(r_{t-1},x_{1:t-1})\\
= & \sum_{r_{t-1}} P(x_t|r_t,r_{t-1},x_{1:t-1})P(r_t|r_{t-1},x_{1:t-1})P(r_{t-1},x_{1:t-1})\\
= & \sum_{r_{t-1}} P(x_t|x^{(r)}_t)P(r_t|r_{t-1})P(r_{t-1},x_{1:t-1}) \tag{2}
\end{align}
$$
**Remark:**
+ **By definition, $x^{(r)}_t \subset x_{1:t-1}$. Or more precisely, $x^{(r)}_t$ is the last $r_t$(may be zero) elements of $x_{1:t-1}$.**
+ **$P(x_{t}|x_{1:t-1},r_t\ or\ r_{t-1})$ relies only on $x_t^{(r)}$, i.e. $P(x_{t}|x_{1:t-1},r_t\ or\ r_{t-1})=P(x_{t}|x_t^{(r)})$. Intuitively, if $x_{t}$ and $x_{t-1}$ belong to the same partition, then  the prediction probability of $x_{t}$ is $P(x_{t}|\eta_{t-1})$, where $\eta_{t-1}$ is determined by $x^{(r)}_t$, so we have $P(x_{t}|\eta_{t-1})=P(x_{t}|x^{(r)}_t)$. On the other hand, if $x_{t}$ belongs to a new partition, then $x_{t}$ is irrevelant to $\eta_{t-1}$, the equation still holds: $P(x_{t})=P(x_{t}|\eta_{t-1})=P(x_{t}|x^{(r)}_t)$**
+ **$r_t$ only depends on $r_{t-1}$.**

**According to $(2)$, subproblem 1 is further devided into 3 minus subproblems:**
**(a). Construct a boundary condition, or the initial distribution of $P(r_{t-1}, x_{1:t-1})$. When $t=1$, $P(r_{t-1}, x_{1:t-1})=P(r_0,x_{0})$, note that $x_0$ is an empty set, so the problem can be simplified to the construction of $P(r_0)$.**
**(b). Construct a transition distribution $P(r_t|r_{t-1})$**
**(c). Construct a prediction distribution $P(x_t|x^{(r)}_t)$**
**As to (a), the simplist way is to assume a changepoint occured before the first data. In such cases we place all of the probability mass for the initial run length at zero, i.e. $P(r_0=0) = 1$.**
**As to (b), $P(r_t|r_{t-1}) $has non-zero mass at only two outcomes: the run length either continues to grow and $r_t = r_{t−1} + 1$ or a changepoint occurs and $r_t = 0$. <span style="color:darkred">This reasoning process is the same as the one used in inferring motality rate, which can be easily interpretated by a hazard function($hazard=\frac{density}{survival}$). </span> Here we have:**
$$
\begin{align}
P(r_t|r_{t-1}) & = &\\
& H(r_{t-1}+1), & if\ r_t=0,\\
& {1-H(r_{t-1}+1)}, &if\ r_t=r_{t-1}+1\\
& 0, & other. \tag{3}
\end{align}
$$
**Where**
$$
H(g=\tau)=\frac{P_{gap}(g=\tau)}{\sum_{i=\tau}^{infinity}{P_{gap}(g=i)}}. \tag{4}
$$
**Any prior information about run length transitions can be easily incoporated by the form of $P_{gap}(g)$. For example, if $P_{gap}(g)$ is an exponential(or geometric, if discrete) distribution with timesacle $\lambda$, then the hazard function is constant(independent of previous run length), $H(g)=1/\lambda$. It is the same as making no prior assumptions on run length transitions.**
**As to (c). Recall that in bayesian point of view, if every data point $x$ subject to a distribution $P(x|\eta)$, and $\eta$ subject to a prior distribution $P(\eta|\theta)$, where $\theta$ is the hyperparameters. Then the posterior distribution of $\eta$ can be shown satisfying $P(\eta|x)=P(x|\eta)P(\eta|\theta)/P(x)$. Here in this problem. let $\eta_t$ be the parameter inferred from $x^{(r)}_t$, then we have:**
$$
P(\eta_t|x^{(r)}_t)=P(x^{(r)}_t|\eta_t)P(\eta_t|\theta)/P(x_{1:t-1}). \tag{5}
$$
**And the prediction distribution can be easily calculated when disintegrated by $\eta_t$(bayesian prediction problem):**
$$
\begin{align}
P(x_t|x^{(r)}_t)&=\int_{\eta_t}P(x_t|\eta_t)P(\eta_t|x^{(r)}_t)\\
&=\int_{\eta_t}P(x_t|\eta_t) \tag{6}
\end{align}
$$
**So far we have only one problem left, i.e. the derivation of $(5)$. There are three parts in $(5)$, they are the likelihood function $P(x^{(r)}_t|\eta_t)$(or the sample distribution $P(x|\eta_t)$, the two are technically implying the same thing), the prior distribution $P(\eta_t|\theta)$ and the normalizing constant $P(x_{1:t-1})$. Because the normalizing constant is already dealt in subproblem 2, here we only need to find a proper likelihood function and a prior distribution.**
**Generally, we need to calculate the posterior distribution of $\eta_t$ directly from $(5)$, but sometimes it would be vary time consuming. Particularly, if both likelihood and prior distribution are from exponential family, the posterior distribution will also from exponential family, and allow inference with a finite number of sufficient statistics which can be calculated incrementally as data arrives.**
**Here is a general case where there is a normal ikelihood function with unkonwn mean $\mu$ and precision $\tau$(or inverse of variance,$1/\sigma^2$). Usually, if $\tau$ is known, $\mu$ follows a noramal distribution with mean $\mu_0$ and precision $k_0$, and if $\mu$ is known, $\tau$ follows a gamma distribution with shape $\alpha_0$ and rate(1/sclale) $\beta_0$. When both $\mu$ and $\tau$ are unknown, it is usual to assme they follow a normal-gamma distribution with parameters $\mu_0,k_0,\alpha_0$ and $\beta_0$ :**
$$
NG(\mu,\tau) = N(\mu|(k_0\tau)^{-1})Ga(\tau|\alpha_0,\beta_0)
$$
**see [exponettial family](https://en.wikipedia.org/wiki/Exponential_family#Bayesian_estimation:_conjugate_distributions) and [Conjugate Bayesian analysis of the Gaussian distribution](http://www-devel.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf) for more information.**
**** Note****
1. **In the above example c($k_0$,$\alpha_0$) and c($\mu_0$,$\beta_0$) act the same as effective number of observations $v$ and the amount contribute to the sufficient statistic $\chi$ repestectively.**
2. **In the bayesian prediction problem $(6)$, calculate the integration $\int_{\eta} p(x_t|\eta)p(\eta|x_{1:t-1})$ directly will be too time consuming. Alternatively, we can use the expected value of $\eta$ instead of the whole distribution to approximate the probability. i.e.** 
$$
\int_{\eta} p(x_t|\eta)p(\eta|x_{1:t-1})\approx p(x_t|E(\eta|x_{1:t-1}))
$$



## Reference

**[1]: Adams, Ryan Prescott, and David JC MacKay. "Bayesian online changepoint detection." arXiv preprint arXiv:0710.3742 (2007).**

**[2]: Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.2σ2 (2007): 16.**

**[3]: McLachlan, Geoffrey, and Thriyambakam Krishnan. *The EM algorithm and extensions*. Vol. 382. John Wiley & Sons, 2007.**

**[4]: Murphy, Kevin P. *Machine learning: a probabilistic perspective*. MIT press, 2012.**

**[5]: Gelman, Andrew, et al. *Bayesian data analysis*. Chapman and Hall/CRC, 2013.**

**[6]: Soland, Richard M. "Bayesian analysis of the Weibull process with unknown scale and shape parameters." *IEEE Transactions on Reliability* 18.4 (1969): 181-184.**

**[7]: Walker, Stephen G. "Sampling the Dirichlet mixture model with slices." *Communications in Statistics—Simulation and Computation®* 36.1 (2007): 45-54.**

**[8]: Van Gael, Jurgen, et al. "Beam sampling for the infinite hidden Markov model." *Proceedings of the 25th international conference on Machine learning*. ACM, 2008.**



## License

Do whatever you want with the files in this repository, under the condition that you obey the license of the included packages/softwares/publications.