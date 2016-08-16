# Changepoints

[TOC "float:left"]

### 1. Online Changepoint Detection


### 2. Level&Trend Change(CUMSUM)
An implementation of Page's mean shift detection algorithm(CUMSUM).

**[R code](https://github.com/chenhaotian/Changepoints/tree/master/leveltrend_change)**

**reference:**
[Changepoint Detection in Climate Time Series with Long-Term Trends](http://journals.ametsoc.org/doi/full/10.1175/JCLI-D-12-00704.1)

##### Raw data
Real change point positions: 1000,1700,2400,3100
![](./leveltrend_change/raw.png)
##### Detected level and trend
Detected change point positions: 996,1695,2400,3104
![](./leveltrend_change/level.png)
![](./leveltrend_change/trend.png)

### 3. Burst Detection

### 4. Distrubution Change

### 5. Structural Change
An implementation of  Makram Talih's graph structural change detection algorithm.

**[R code](https://github.com/chenhaotian/Changepoints/tree/master/structural_change)**

**reference:**
[Structural learning with time‐varying components: tracking the cross‐section of financial time series](https://www.researchgate.net/publication/4914219_Structural_learning_with_time-varying_components_Tracking_the_crosssection_of_financial_time_series)

##### Raw data

![](./structural_change/Screenshot from 2015-04-30 16:35:28.png)

##### Original gaussian graph

![](./structural_change/Screenshot from 2015-04-30 16:35:48.png)

##### Detucted gaussian graph

![](./structural_change/Screenshot from 2015-04-30 16:36:38.png)

### 6. Outlier

reference:
