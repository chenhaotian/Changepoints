# Bayesian online change point detection
An R implementation of the paper: Bayesian Online Changepoint Detection(Ryan Prescott Adams and David J.C. MacKay 2007), **with derivation process explanations, additional conjugate distributions and examples,** the following chapters will cover:
+ [Basic idea and function usage](#usage)
+ [Derivation process explanations](#derivation)
+ Normal evidence distribution and normal-gamma prior + [examples](#examples)
+ Poisson evidence distribution and gamma prior + examples
+ Binomial evidence distribution and beta prior + examples

## Usage
Load the function `onlinechangepoint()` directly from source:
```R
source("https://raw.githubusercontent.com/chenhaotian/Changepoints/master/online_detection/bayesian.r")
```
Parameters:
```R
args(onlinechangepoint)
```
+ X: numeric, observations.
+ model: character, specifying model
	+ nng: normal evidence(observation distribution) and normal-gamma prior
	+ pg: poisson evidence and gamma prior
	+ bb: binomial evidence and beta prior
	+ g: gamma evidence
+ mu0, k0, alpha0, beta0: numeric, specifying hyper parameters.
	+ mu0, k0, alpha0, beta0: normal-gamma parameter, when model="nng"
	+ alpha0, beta0: gamma parameters when model="pg"(these two names alpha0 and beta0 are shared with "nng")
	+ alpha0, beta0: beta parameters when model="bb"
+ lambda: numeric, parameter of the exponential hazard function.(act as transition distribution)
+ FILTER: if `P(r_t|x_{1:t})<FILTER`, this `r_t` will be omitted in next round of online calculation.

See [Examples](#examples) below for detailed explanation in function usage.

## Derivation
Adams&Mackay's paper interpret the change point detection problem with an infinite state hidden markov model. Hereby the definitions:
+ $r_t$ is the time lenght since last changepoint, $x_{t}^{(r)}$ denote the data sets associated with run $r_t$. **As $r$ may be zero, $x^{(r)}$ may be empty**.The hidden states are the run lengths $r$ of current observation, they evolve as the following way:
![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/online_detection/state_transition.png)
+ For each parition $\rho$, the data within it are i.i.d. from some probability distribution $P (x_t | \eta_\rho)$
+ The discrete **a priori**(Latin phrase, means 'from the earlier',contrast to **a posteriori**,'from the latter') probability distribution over the interval between changepoints is denoted as $P_{gap}(g)$, i.e. the probability of current run ending up at length $g$.(recall geometric distribution)
+ Predictive distribution $P(x_{t+1}|r_t,x_t^{(r)})$
  Posterior distribution $P(r_t|x_{1:t})$
  Joint distribution $P(r_t,x_{1:t})$

With the above definitions, our goal here is to <span style="color:darkred">**recursively**(use recusive methods for the purpose of online calculation)</span> compute the probability distribution of **the length of the current “run”,** or time since the last changepoint, denote this **posterior distribution** as $P(r_t|x_{1:t})$
According to bayesian therom,
$$P(r_t|x_{1:t})=\frac{P(r_t,x_{1:t})}{P(x_{1:t})}, \tag{1}$$
we can devide the calculation of the conditional distribution into two subproblems:
**1.** Calculate the joint distributioin $P(r_t,x_{1:t})$.
**2.** Calculate the normalizing constant $P(x_{1:t})$.
As for subproblem 2, $P(x_{1:t})$ can be easily calculated by integrate out $r_t$ in subproblem 1. i.e. $P(x_{1:t})=\sum_{r_t}P(r_t,x_{1:t})$. So basically we need only focus on subproblem 1, the derivation of $P(r_t,x_{1:t})$.
**Keep in mind that to make the online calculation possible, the derived represntation of subproblem 1 must be recursive.** i.e. when the previous information is known, one need to derive a relation between previous information and current estimation. Here in subproblem 1, the previous information is $P(r_{t-1}|x_{1:t-1})$, and the current estimation is $P(r_t|x_{1:t})$. In order to relate these two distributions, first we need to <span style="color:darkred">**disintegrate**</span>  $P(r_t|x_{1:t})$ by $r_{t-1}$, and then use bayes rule to extract $P(r_{t-1}|x_{1:t-1})$:
$$
\begin{align}
P(r_t,x_{1:t}) & = \sum_{r_{t-1}} P(r_t,r_{t-1},x_{1:t})\\
= & \sum_{r_{t-1}} P(r_t,x_t|r_{t-1},x_{1:t-1})P(r_{t-1},x_{1:t-1})\\
= & \sum_{r_{t-1}} P(x_t|r_t,r_{t-1},x_{1:t-1})P(r_t|r_{t-1},x_{1:t-1})P(r_{t-1},x_{1:t-1})\\
= & \sum_{r_{t-1}} P(x_t|x^{(r)}_t)P(r_t|r_{t-1})P(r_{t-1},x_{1:t-1}) \tag{2}
\end{align}
$$
Remark:
+ By definition, $x^{(r)}_t \subset x_{1:t-1}$. Or more precisely, $x^{(r)}_t$ is the last $r_t$(may be zero) elements of $x_{1:t-1}$.
+ $P(x_{t}|x_{1:t-1},r_t\ or\ r_{t-1})$ relies only on $x_t^{(r)}$, i.e. $P(x_{t}|x_{1:t-1},r_t\ or\ r_{t-1})=P(x_{t}|x_t^{(r)})$. Intuitively, if $x_{t}$ and $x_{t-1}$ belong to the same partition, then  the prediction probability of $x_{t}$ is $P(x_{t}|\eta_{t-1})$, where $\eta_{t-1}$ is determined by $x^{(r)}_t$, so we have $P(x_{t}|\eta_{t-1})=P(x_{t}|x^{(r)}_t)$. On the other hand, if $x_{t}$ belongs to a new partition, then $x_{t}$ is irrevelant to $\eta_{t-1}$, the equation still holds: $P(x_{t})=P(x_{t}|\eta_{t-1})=P(x_{t}|x^{(r)}_t)$
+ $r_t$ only depends on $r_{t-1}$.

According to $(2)$, subproblem 1 is further devided into 3 minus subproblems:
**(a).** Construct a **boundary condition**, or the initial distribution of $P(r_{t-1}, x_{1:t-1})$. When $t=1$, $P(r_{t-1}, x_{1:t-1})=P(r_0,x_{0})$, note that $x_0$ is an empty set, so the problem can be simplified to the construction of $P(r_0)$.
**(b).** Construct a transition distribution $P(r_t|r_{t-1})$
**(c).** Construct a prediction distribution $P(x_t|x^{(r)}_t)$
As to (a), the simplist way is to assume a changepoint occured before the first data. In such cases we place all of the probability mass for the initial run length at zero, i.e. $P(r_0=0) = 1$.
As to (b), $P(r_t|r_{t-1}) $has non-zero mass at only two outcomes: the run length either continues to grow and $r_t = r_{t−1} + 1$ or a changepoint occurs and $r_t = 0$. <span style="color:darkred">This reasoning process is the same as the one used in inferring **motality rate**, which can be easily interpretated by a **hazard function**($hazard=\frac{density}{survival}$). </span> Here we have:
$$
\begin{align}
P(r_t|r_{t-1}) & = &\\
& H(r_{t-1}+1), & if\ r_t=0,\\
& {1-H(r_{t-1}+1)}, &if\ r_t=r_{t-1}+1\\
& 0, & other. \tag{3}
\end{align}
$$
Where
$$
H(g=\tau)=\frac{P_{gap}(g=\tau)}{\sum_{i=\tau}^{infinity}{P_{gap}(g=i)}}. \tag{4}
$$
Any prior information about run length transitions can be easily incoporated by the form of $P_{gap}(g)$. For example, if $P_{gap}(g)$ is an exponential(or geometric, if discrete) distribution with timesacle $\lambda$, then the hazard function is constant(independent of previous run length), $H(g)=1/\lambda$. It is the same as making no prior assumptions on run length transitions.
As to (c). Recall that in bayesian point of view, if every data point $x$ subject to a distribution $P(x|\eta)$, and $\eta$ subject to a prior distribution $P(\eta|\theta)$, where $\theta$ is the hyperparameters. Then the posterior distribution of $\eta$ can be shown satisfying $P(\eta|x)=P(x|\eta)P(\eta|\theta)/P(x)$. Here in this problem. let $\eta_t$ be the parameter inferred from $x^{(r)}_t$, then we have:
$$
P(\eta_t|x^{(r)}_t)=P(x^{(r)}_t|\eta_t)P(\eta_t|\theta)/P(x_{1:t-1}). \tag{5}
$$
And the prediction distribution can be easily calculated when disintegrated by $\eta_t$(**bayesian prediction problem**):
$$
\begin{align}
P(x_t|x^{(r)}_t)&=\int_{\eta_t}P(x_t|\eta_t)P(\eta_t|x^{(r)}_t)\\
&=\int_{\eta_t}P(x_t|\eta_t) \tag{6}
\end{align}
$$
So far we have only one problem left, i.e. the derivation of $(5)$. There are three parts in $(5)$, they are the **likelihood function** $P(x^{(r)}_t|\eta_t)$(or the sample distribution $P(x|\eta_t)$, the two are technically implying the same thing), the **prior distribution** $P(\eta_t|\theta)$ and the normalizing constant $P(x_{1:t-1})$. Because the normalizing constant is already dealt in subproblem 2, here we only need to find a proper likelihood function and a prior distribution.
Generally, we need to calculate the posterior distribution of $\eta_t$ directly from $(5)$, but sometimes it would be vary time consuming. Particularly, if both likelihood and prior distribution are from exponential family, the posterior distribution will also from exponential family, and allow inference with a finite number of sufficient statistics which can be calculated incrementally as data arrives.
Here is a general case where there is a normal ikelihood function with unkonwn mean $\mu$ and precision $\tau$(or inverse of variance,$1/\sigma^2$). Usually, if $\tau$ is known, $\mu$ follows a noramal distribution with mean $\mu_0$ and precision $k_0$, and if $\mu$ is known, $\tau$ follows a gamma distribution with shape $\alpha_0$ and rate(1/sclale) $\beta_0$. When both $\mu$ and $\tau$ are unknown, it is usual to assme they follow a normal-gamma distribution with parameters $\mu_0,k_0,\alpha_0$ and $\beta_0$ :
$$
NG(\mu,\tau) = N(\mu|(k_0\tau)^{-1})Ga(\tau|\alpha_0,\beta_0)
$$
see [exponettial family](https://en.wikipedia.org/wiki/Exponential_family#Bayesian_estimation:_conjugate_distributions) and [Conjugate Bayesian analysis of the Gaussian distribution](http://www-devel.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf) for more information.
** Note**
1. In the above example c($k_0$,$\alpha_0$) and c($\mu_0$,$\beta_0$) act the same as effective number of observations $v$ and the amount contribute to the sufficient statistic $\chi$ repestectively.
2. In the bayesian prediction problem $(6)$, calculate the integration $\int_{\eta} p(x_t|\eta)p(\eta|x_{1:t-1})$ directly will be too time consuming. Alternatively, we can use the expected value of $\eta$ instead of the whole distribution to approximate the probability. i.e. 
$$
\int_{\eta} p(x_t|\eta)p(\eta|x_{1:t-1})\approx p(x_t|E(\eta|x_{1:t-1}))
$$

## Examples
Here are several examples in applying the algorithm.
###### 1. Level change(normal observation):
Different mean and same volatility.
```R
## Generate random samples with real changepoints: 100,170,240,310,380
X <- c(rnorm(100,sd=0.5),rnorm(70,mean=5,sd=0.5),rnorm(70,mean = 2,sd=0.5),rnorm(70,sd=0.5),rnorm(70,mean = 7,sd=0.5))
## online changepoint detection for series X
resX <- onlinechangepoint(X,
                          model = "nng",
                          mu0=0.7,k0=1,alpha0=1/2,beta0=1, #initial parameters
                          bpmethod = "mean",
                          lambda=50, #exponential hazard
                          FILTER=1e-3)
tmpplot(resX,X) # visualize original data and run length distribution
```
The last statement `tmpplot(rexX,X)` produce a graph comparing the original data(upper part) and the run lengthes(lower part):
![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/online_detection/online_level.png)
###### 2. Volatility change(normal observation):
Different volatility and same mean.
```R
## Generate random samples with real changepoints: 100,170,240,310,380.
Y <- c(rnorm(100,sd=0.5),rnorm(70,sd=1),rnorm(70,sd=3),rnorm(70,sd=1),rnorm(70,sd=0.5))
## online changepoint detection for series Y
resY <- onlinechangepoint(Y,
                          model = "nng",
                          mu0=0,k0=0.5,alpha0=1/2,beta0=1,
                          bpmethod = "mean",
                          lambda=50, #exponential hazard
                          FILTER=1e-3)
tmpplot(resY,Y)
```
Original data and inferred run length distribution:
![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/online_detection/online_variance.png)
###### 3. Level change(poisson observation):
```R
Z <- c(rpois(100,lambda=5),rpois(70,lambda=7),rpois(70,lambda=5),rpois(70,lambda=6),rpois(70,lambda=5))
## online changepoint detection for series Z
resZ <- onlinechangepoint(Z,
                          model = "pg",
                          alpha0=10,beta0=1,
                          bpmethod = "mean",
                          lambda=10, #exponential hazard
                          FILTER=1e-3)
tmpplot(resZ,Z)
```
Original data and inferred run length distribution:
![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/online_detection/online_poisson.png)


## Reference
[1] [Bayesian Online Changepoint Detection](http://arxiv.org/abs/0710.3742)
[2] [Checking wether a coin is fair](https://en.wikipedia.org/wiki/Checking_whether_a_coin_is_fair
)
[3] [DeGroot,Optimal Statistical Decisions, chapter 9]()
[4] [Gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution)
[5] [Fink, D. 1995 A Compendium of Conjugate Priors. In progress report: Extension and enhancement of methods for setting data quality objectives. (DOE contract 95‑831).]()
