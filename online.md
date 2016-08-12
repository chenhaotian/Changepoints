### 5. Online Changepoint Detection

Changepoints are abrupt variations in the generative parameters of a data sequence.

We compute the probability distribution of **the length of the current “run,”** or time since the last changepoint, using a simple message-passing algorithm.

Observations can be devided into non-overlapping paritions, data within each parition are i.i.d.

Exponential family likelihoods allow inference with a finite number of sufficient statistics which can be calculated incrementally as data arrives

**Boundary Condition: initial distribution**
