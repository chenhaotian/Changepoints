### 5. Online Changepoint Detection

Changepoints are abrupt variations in the generative parameters of a data sequence.

We , using a simple message-passing algorithm.

Observations can be devided into non-overlapping paritions, data within each parition are i.i.d.

Exponential family likelihoods allow inference with a finite number of sufficient statistics which can be calculated incrementally as data arrives

**Boundary Condition: initial distribution**

###### Definition:
+ $r_t$ is the time lenght since last changepoint, $x_{t}^{(r)}$ denote the data sets associated with run $r_t$. **As $r$ may be zero, $x^{(r)}$ may be empty**.
+ For each parition $\rho$, the data within it are i.i.d. from some probability distribution $P (x_t | \eta_\rho)$
+ The discrete **a priori**(Latin phrase, means 'from the earlier',contrast to **a posteriori**,'from the latter') probability distribution over the interval between changepoints is denoted as $P_{gap}(g)$, i.e. the probability of current run ending up at length $g$.(recall geometric distribution)
+ Predictive distribution $P(x_{t+1}|r_t,x_t^{(r)})$
  Posterior distribution $P(r_t|x_{1:t})$
  Joint distribution $P(r_t,x_{1:t})$



Generate a recursive message-passing algorithm for the joint distribution over the current run length and the data

###### Basic Principle:
As from the above definition, our goal here is to <span style="color:darkred">**recursively**(use recusive methods for the purpose of online calculation)</span> compute the probability distribution of **the length of the current “run,”** or time since the last changepoint, denote this **posterior distribution** as $P(r_t|x_{1:t})$
According to bayesian therom,
$$P(r_t|x_{1:t})=\frac{P(r_t,x_{1:t})}{P(x_{1:t})}, \tag{1}$$
we can devide the calculation of the conditional distribution into two subproblems according to $(1)$:
1.Calculate the joint distributioin $P(r_t|x_{1:t})$
2.Calculate the normalizing constant $P(x_{1:t})$
As for subproblem 2, $P(x_{1:t})$ can be easily calculated by integrate out $r_t$ in subproblem 1. i.e. $P(x_{1:t})=\sum_{r_t}P(r_t,x_{1:t})$. So basically we need only focus on subproblem 1, the derivation of $P(r_t|x_{1:t})$.
**Keep in mind that to make the online calculation possible, the derived represntation of subproblem 1 must be recursive.** i.e. when the previous information is known, one need to derive a relation between previous information and current estimation. Here in subproblem 1, the previous information is $P(r_{t-1}|x_{1:t-1})$, and the current estimation is $P(r_t|x_{1:t})$. In order to relate these two distributions, first we need to <span style="color:darkred">**disintegrate**</span>  $P(r_t|x_{1:t})$ by $r_{t-1}$, and then use bayes rule to extract $P(r_{t-1}|x_{1:t-1})$:
$$
\begin{align}
P(r_t|x_{1:t}) & = \sum_{r_{t-1}} P(r_t,r_{t-1},x_{1:t})\\
= & \sum_{r_{t-1}} P(r_t,x_t|r_{t-1},x_{1:t-1})P(r_{t-1},x_{1:t-1})\\
= & \sum_{r_{t-1}} P(x_t|r_t,r_{t-1},x_{1:t-1})P(r_t|r_{t-1},x_{1:t-1})P(r_{t-1},x_{1:t-1})\\
= & \sum_{r_{t-1}} P(x_t|x^{(r)}_t)P(r_t|r_{t-1})P(r_{t-1},x_{1:t-1}) \tag{2}
\end{align}
$$
Remark:
+ By definition, $x^{(r)}_t \subset x_{1:t-1}$. Or more precisely, $x^{(r)}_t$ is the last $r_t$ elements of $x_{1:t-1}$.
+ $P(x_{t}|x_{1:t-1},r_t\ or\ r_{t-1})$ relies only on $x_t^{(r)}$, i.e. $P(x_{t}|x_{1:t-1},r_t\ or\ r_{t-1})=P(x_{t}|x_t^{(r)})$. Intuitively, if $x_{t}$ and $x_{t-1}$ belong to the same partition, then  the prediction probability of $x_{t}$ is $P(x_{t}|\eta_{t-1})$, where $\eta_{t-1}$ is determined by $x^{(r)}_t$, so we have $P(x_{t}|\eta_{t-1})=P(x_{t}|x^{(r)}_t)$. On the other hand, if $x_{t}$ belongs to a new partition, then $x_{t}$ is irrevelant to $\eta_{t-1}$, the equation still holds: $P(x_{t})=P(x_{t}|\eta_{t-1})=P(x_{t}|x^{(r)}_t)$
+ $r_t$ only depends on $r_{t-1}$.

According to $(2)$, subproblem 1 is further devided into 3 minus subproblems:
(a). prediction distribution $P(x_t|x^{(r)}_t)$
(b). 


Where 

intuitions:
Recall from state space models, to get a posterior estimation of true state, one need to prepare a predictive distribution and a observation distirbution.

assume:
+ $P(x_{t+1}|r_t,x^{(r)}_t)$ is known(or computable)


