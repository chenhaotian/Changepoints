### 5. Online Changepoint Detection

Changepoints are abrupt variations in the generative parameters of a data sequence.

We , using a simple message-passing algorithm.

Observations can be devided into non-overlapping paritions, data within each parition are i.i.d.

Exponential family likelihoods allow inference with a finite number of sufficient statistics which can be calculated incrementally as data arrives

**Boundary Condition: initial distribution**

definition:
+ $r_t$ is the time lenght since last changepoint, $x_{t}^{(r)}$ denote the data sets associated with run $r_t$. **As $r$ may be zero, $x^{(r)}$ may be empty**.
+ For each parition $\rho$, the data within it are i.i.d.
+ The discrete **a priori**(Latin phrase, means 'from the earlier',contrast to **a posteriori**,'from the latter') probability distribution over the interval between changepoints is denoted as $P_{gap}(g)$, i.e. the probability of current run ending up at length $g$.(recall geometric distribution)
+ Predictive distribution $P(x_{t+1}|r_t,x_t^{(r)})$
  Posterior distribution $P(r_t|x_{1:t})$
  Joint distribution $P(r_t,x_{1:t})$

target:
+ The goal is to (**recursively**) compute the probability distribution of **the length of the current “run,”** or time since the last changepoint, denote this **posterior distribution** as $P(r_t|x_{1:t})$
+ Generate a recursive message-passing algorithm for the joint distribution over the current run length and the data

intuitions:
Recall from state space models, to get a posterior estimation of true state, one need to prepare a predictive distribution and a observation distirbution.



assume:
+ $P(x_{t+1}|r_t,x^{(r)}_t)$ is known(or computable)