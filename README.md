# Bayesian online change point detection 
An R implementation of the paper: Bayesian Online Changepoint Detection(Ryan Prescott Adams and David J.C. MacKay 2007), **with derivation process explanations, additional conjugate distributions and examples,** the following chapters will cover: 
+ [Basic idea and function usage](#usage) 
+ [Derivation process explanations](#derivation) 
+ Normal evidence distribution and normal-gamma prior + [examples](#examples) 
+ Poisson evidence distribution and gamma prior + examples 
+ Binomial evidence distribution and beta prior + examples 
+ [Additional readings](#reference) 
 
## Usage 
Load the function `onlinechangepoint()` directly from source: 
```R 
source("https://raw.githubusercontent.com/chenhaotian/Changepoints/master/online_detection/bayesian.r") 
``` 
Parameters: 
```R 
args(onlinechangepoint) 
``` 
+ **X** : numeric, observations. 
+ **model** : character, specifying model to be used, can be one of c("nng","pg","bb","g") 
	+ nng: normal evidence(observation distribution) and normal-gamma prior 
	+ pg: poisson evidence and gamma prior 
	+ bb: binomial evidence and beta prior 
	+ g: gamma evidence 
+ **mu0, k0, alpha0, beta0** : numeric, specifying hyper parameters. 
	+ mu0, k0, alpha0, beta0: normal-gamma parameter, when model="nng" 
	+ alpha0, beta0: gamma parameters when model="pg"(these two names alpha0 and beta0 are shared with "nng") 
	+ alpha0, beta0: beta parameters when model="bb" 
+ **lambda** : numeric, parameter of the exponential hazard function.(act as transition distribution) 
+ **FILTER** : if `P(r_t|x_{1:t})<FILTER`, this `r_t` will be omitted in next round of online calculation. 
 
See [Examples](#examples) below for detailed explanation in function usage. 
 
## Derivation 
Adams&Mackay's paper interpret the change point detection problem with an infinite state hidden markov model. Hereby the definitions: 
+ ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/onaescgymdkflzwjbxihvruqpt1.gif) is the time lenght since last changepoint, ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/hucpgriknmvosjbweatzdfqlyx2.gif) denote the data sets associated with run ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/onaescgymdkflzwjbxihvruqpt1.gif). **As ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/lvdwciohbxgsruznkjymfepatq3.gif) may be zero, ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/qovteagyfbhnlmdwzikjpxucrs4.gif) may be empty**.The hidden states are the run lengths ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/lvdwciohbxgsruznkjymfepatq3.gif) of current observation, they evolve as the following way: 
![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/online_detection/state_transition.png) 
+ For each parition ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/vpcfmgdnjbxokyizwsaqlturhe5.gif), the data within it are i.i.d. from some probability distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/brxqfcglmyvozkpwuajhdisnte6.gif) 
+ The discrete **a priori**(Latin phrase, means 'from the earlier',contrast to **a posteriori**,'from the latter') probability distribution over the interval between changepoints is denoted as ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/iuapmjgqwnkdesvtbhrxylfzco7.gif), i.e. the probability of current run ending up at length ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/iwaubeyrvgsdplczxjnqothkmf8.gif).(recall geometric distribution) 
+ Predictive distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/figyubmxcvkdwshoajntzplqre9.gif) 
  Posterior distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/vifmsqkablhcxgjetdzynrowpu10.gif) 
  Joint distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/apgcuhblqreiyjmskwxnovftdz11.gif) 
 
With the above definitions, our goal here is to <span style="color:darkred">**recursively**(use recusive methods for the purpose of online calculation)</span> compute the probability distribution of **the length of the current “run”,** or time since the last changepoint, denote this **posterior distribution** as ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/vifmsqkablhcxgjetdzynrowpu10.gif) 
According to bayesian therom, 

![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/voqyrljntidgfumachekxbpszw1.gif)
 
we can devide the calculation of the conditional distribution into two subproblems: 
**1.** Calculate the joint distributioin ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/apgcuhblqreiyjmskwxnovftdz11.gif). 
**2.** Calculate the normalizing constant ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/yxjdgnsclbtfupvmiozwrkhqea12.gif). 
As for subproblem 2, ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/yxjdgnsclbtfupvmiozwrkhqea12.gif) can be easily calculated by integrate out ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/onaescgymdkflzwjbxihvruqpt1.gif) in subproblem 1. i.e. ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/wkjnmvpfahsdtlucgoqbexzryi13.gif). So basically we need only focus on subproblem 1, the derivation of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/apgcuhblqreiyjmskwxnovftdz11.gif). 
**Keep in mind that to make the online calculation possible, the derived represntation of subproblem 1 must be recursive.** i.e. when the previous information is known, one need to derive a relation between previous information and current estimation. Here in subproblem 1, the previous information is ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/qyvtdwxlzbnipgeukrsocahmfj14.gif), and the current estimation is ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/vifmsqkablhcxgjetdzynrowpu10.gif). In order to relate these two distributions, first we need to <span style="color:darkred">**disintegrate**</span>  ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/vifmsqkablhcxgjetdzynrowpu10.gif) by ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/njpqafxlioedbuwyrkcthgvszm15.gif), and then use bayes rule to extract ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/qyvtdwxlzbnipgeukrsocahmfj14.gif): 

![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/wzhkofguclnjvtadsieybmxqpr2.gif)
 
Remark: 
+ By definition, ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/ydaixbelosuvwfgmqnpjrztckh16.gif). Or more precisely, ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/fsvyixwctndzhjarkbeomuglpq17.gif) is the last ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/onaescgymdkflzwjbxihvruqpt1.gif)(may be zero) elements of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/owxerpuytigfkvdcmbqhnazljs18.gif). 
+ ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/rpvtnuxmoibsawcdzkyfehlqgj19.gif) relies only on ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/wknvqcyfaxmgblerpjuoshzitd20.gif), i.e. ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/fyrnjimzpudthkwecsxqoavbgl21.gif). Intuitively, if ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/hlzketirbnsjocyudxqpfvwmag22.gif) and ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/xwlfikzhmdsrbvotygqcaneupj23.gif) belong to the same partition, then  the prediction probability of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/hlzketirbnsjocyudxqpfvwmag22.gif) is ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/avhyrfdqstzcgkiebunoxpmjlw24.gif), where ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/omnfiuqdgaryhwkbplvzjtsexc25.gif) is determined by ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/fsvyixwctndzhjarkbeomuglpq17.gif), so we have ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/bdlhozyvsgitpanuqrjwcefmkx26.gif). On the other hand, if ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/hlzketirbnsjocyudxqpfvwmag22.gif) belongs to a new partition, then ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/hlzketirbnsjocyudxqpfvwmag22.gif) is irrevelant to ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/omnfiuqdgaryhwkbplvzjtsexc25.gif), the equation still holds: ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/gzwftvspucdxyhaqeinmljrkbo27.gif) 
+ ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/onaescgymdkflzwjbxihvruqpt1.gif) only depends on ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/njpqafxlioedbuwyrkcthgvszm15.gif). 
 
According to ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/bdfmcwjxoqikurzgtlaevhpysn28.gif), subproblem 1 is further devided into 3 minus subproblems: 
**(a).** Construct a **boundary condition**, or the initial distribution of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/wrdjoebpcxlmakgnhfvyiszutq29.gif). When ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/xlqpjcsgtaedwfyimhuvokrznb30.gif), ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/mxgnluhipcsdraejzbqkfoytwv31.gif), note that ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/nwtdqulmaybepjghkcovfizrxs32.gif) is an empty set, so the problem can be simplified to the construction of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/rbxgijvzmdhcsuplqetyfankow33.gif). 
**(b).** Construct a transition distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/mpwxqrjhvcygftaoiludnszbek34.gif) 
**(c).** Construct a prediction distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/uthoepncmsqjxkivrwblgdfayz35.gif) 
As to (a), the simplist way is to assume a changepoint occured before the first data. In such cases we place all of the probability mass for the initial run length at zero, i.e. ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/dzmhbwpfolxvngjeraykiuqstc36.gif). 
As to (b), ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/ivwetmahsxkypngfzudcrjolqb37.gif)has non-zero mass at only two outcomes: the run length either continues to grow and ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/ntlzuegsqfijmrvdpyckowxhab38.gif) or a changepoint occurs and ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/qnmgblvaepjrzdcuixtsfyowhk39.gif). <span style="color:darkred">This reasoning process is the same as the one used in inferring **motality rate**, which can be easily interpretated by a **hazard function**(![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/dkjrxuqpcwhlbtvyfigzsenamo40.gif)). </span> Here we have: 

![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/zhlciuojqdnatxkbrgfwsevypm3.gif)
 
Where 

![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/pauqyjosrmnbxdtzgkhliewfcv4.gif)
 
Any prior information about run length transitions can be easily incoporated by the form of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/iuapmjgqwnkdesvtbhrxylfzco7.gif). For example, if ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/iuapmjgqwnkdesvtbhrxylfzco7.gif) is an exponential(or geometric, if discrete) distribution with timesacle ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/njmkvgiucltazqerpwhbyfodsx41.gif), then the hazard function is constant(independent of previous run length), ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/ctfmnyriqzsxgvwpkejdlhaubo42.gif). It is the same as making no prior assumptions on run length transitions. 
As to (c). Recall that in bayesian point of view, if every data point ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/xnkcvlqsfmzaubtwohepgdjiyr43.gif) subject to a distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/jbletsroiwvugdakxyhpzcnmqf44.gif), and ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/klxtherzamjspyfgqbduwnivco45.gif) subject to a prior distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/hqijrfundxbzcyltpvesmkwaog46.gif), where ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/sdiceqpgyfbkruhxmatlonwzjv47.gif) is the hyperparameters. Then the posterior distribution of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/klxtherzamjspyfgqbduwnivco45.gif) can be shown satisfying ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/zgpkdhatwbsijmuvrloncqeyfx48.gif). Here in this problem. let ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/orlawtuenszyvmipcxqkbhgfdj49.gif) be the parameter inferred from ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/fsvyixwctndzhjarkbeomuglpq17.gif), then we have: 

![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/duiowmsjnqxlbpytkzgevfcrah5.gif)
 
And the prediction distribution can be easily calculated when disintegrated by ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/orlawtuenszyvmipcxqkbhgfdj49.gif)(**bayesian prediction problem**): 

![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/bxsuecyhlmpjnwfaivokdtqgrz6.gif)
 
So far we have only one problem left, i.e. the derivation of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/jycgitkadxlmzsvefonuwqphbr50.gif). There are three parts in ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/jycgitkadxlmzsvefonuwqphbr50.gif), they are the **likelihood function** ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/uwftygdxkmrnibvzsjaeqlphoc51.gif)(or the sample distribution ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/iudzehmygfoalkvcbqpwsrntjx52.gif), the two are technically implying the same thing), the **prior distribution** ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/wjmroftgkaizcebyldsxhuqnpv53.gif) and the normalizing constant ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/ymoxkcnlwrasfebdghptivqjuz54.gif). Because the normalizing constant is already dealt in subproblem 2, here we only need to find a proper likelihood function and a prior distribution. 
Generally, we need to calculate the posterior distribution of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/orlawtuenszyvmipcxqkbhgfdj49.gif) directly from ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/jycgitkadxlmzsvefonuwqphbr50.gif), but sometimes it would be vary time consuming. Particularly, if both likelihood and prior distribution are from exponential family, the posterior distribution will also from exponential family, and allow inference with a finite number of sufficient statistics which can be calculated incrementally as data arrives. 
Here is a general case where there is a normal ikelihood function with unkonwn mean ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/gjrwhdtekqiaczlmxbuoyfspnv55.gif) and precision ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/edwqxsrjulvbhzkfctmyopiang56.gif)(or inverse of variance,![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/gpxtwluabfrdhcnokviqsyemzj57.gif)). Usually, if ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/edwqxsrjulvbhzkfctmyopiang56.gif) is known, ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/gjrwhdtekqiaczlmxbuoyfspnv55.gif) follows a noramal distribution with mean ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/krhqtvowjuynpbmfgiadszcxle58.gif) and precision ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/msrdgcnkivewubjhqatpxlzofy59.gif), and if ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/gjrwhdtekqiaczlmxbuoyfspnv55.gif) is known, ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/edwqxsrjulvbhzkfctmyopiang56.gif) follows a gamma distribution with shape ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/qgbcsmzdvwkytreonaxfphuijl60.gif) and rate(1/sclale) ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/zqxhasyonumilwedvbcjtrgpfk61.gif). When both ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/gjrwhdtekqiaczlmxbuoyfspnv55.gif) and ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/edwqxsrjulvbhzkfctmyopiang56.gif) are unknown, it is usual to assme they follow a normal-gamma distribution with parameters ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/dwshfucaqyjginlmxkpzvterbo62.gif) and ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/zqxhasyonumilwedvbcjtrgpfk61.gif) : 

![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/zvdigreapumtjlyhkqnxcfwobs7.gif)
 
see [exponettial family](https://en.wikipedia.org/wiki/Exponential_family#Bayesian_estimation:_conjugate_distributions) and [Conjugate Bayesian analysis of the Gaussian distribution](http://www-devel.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf) for more information. 
** Note** 
1. In the above example c(![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/msrdgcnkivewubjhqatpxlzofy59.gif),![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/qgbcsmzdvwkytreonaxfphuijl60.gif)) and c(![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/krhqtvowjuynpbmfgiadszcxle58.gif),![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/zqxhasyonumilwedvbcjtrgpfk61.gif)) act the same as effective number of observations ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/nbzwsdrgpkxjuofmvthilqecya63.gif) and the amount contribute to the sufficient statistic ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/gsxmvjpyrqndioactfzlwkehub64.gif) repestectively. 
2. In the bayesian prediction problem ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/lcsbxwgnpqmvyuotekdihzarfj65.gif), calculate the integration ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/fquwdhkgisyxpatnrzlcjomvbe66.gif) directly will be too time consuming. Alternatively, we can use the expected value of ![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/klxtherzamjspyfgqbduwnivco45.gif) instead of the whole distribution to approximate the probability. i.e.  

![](https://raw.githubusercontent.com/chenhaotian/Changepoints/master/equations/equations/kugtfebvdnymqzixwlsjoapcrh8.gif)
 
 
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

