# Batch forecast modelling with auto-ETS and auto-ARIMA models on 130 time-series data with cross validations.
* Written by: [Ivy Wu S.Y.](https://www.linkedin.com/in/ivy-wusumyi)
* Technologies: R, batch forecast, tseries, ggplot2, cross-validation


## 1. Introduction
Manually fitting a forecast model is time-consuming, alternatively, the model can be estimated using several automatic functions available in the R software to facilitate the whole process. In this batch modeling, automatic ETS and automatic ARIMA functions are used to produce batch forecasting models for 130 time series. The functions will automatically estimate the model parameters and return the most appropriate fitted model by comparing AICc value as an information criterion. Corss validation is used for selecting the best forcest model whereas the forecast accuracy will be evaluated with benchmark methods by horizons.


## 2. About the data
The data used is provided by the International Institute of Forecasters and available in the [M3-Competition package](https://cran.r-project.org/web/packages/Mcomp/index.html) in R. In this batch forecasting, 130 monthly time series are extracted from the M3-competition dataset. They are time series data with IDs ending with “0” from ID 1501 to 2800. There are 6 categories in total describing each time series including Micro, Industry, Macro, Finance, Demographic and Other. The in-sample data length of each time series varies from 48 to 126, whereas all time series have a fixed 18 length of out-sample data.


## 3. Model Selection Strategy 
<img src="/images/strategy.png?raw=true" width="550"><br/>
Given that the statistical software will automatically return the best fitted ETS and ARIMA model for each time series, it is necessary to define a model selection strategy to choose the best model for each time series. As such, each time series is separated into training and test set to perform cross-validation on ETS and ARIMA models. The models are hence compared with forecast accuracy and the best method is chosen for out-sample forecasting.

### 3.1 Time Series Cross Validation with Rolling Origins
To prevent biased estimation of the true error, cross-validation is used rather than simple validation as it provides more and better information about how the model performs on unseen data with repeated trails. Considering the time series characteristics of the data, rolling origin evaluation will be used such that the test set data comes chronologically after the training data. The size of the validation window (test set) is equal to 18 for better simulation and the rolling origins of 6 or 12 is used depending on the size of available data in each time series dataset. Once the parameters have been decided, the cross-validation process is run with multiple rolling times. Both ETS and ARIMA models are fitted with the train dataset and output forecasts per iterations. Each iterations’ forecast is compared with the test dataset (validation window) to calculate the MAPE (mean absolute percentage error). <br/>
<img src="/images/cross_validation.png?raw=true" width="550">

### 3.2. Model Selection Result
By evaluating the forecast accuracy (MAPE) of ETS and ARIMA models for 130 time series, the result indicates that ETS has a better model performance on 71 time series, whereas ARIMA has a better model performance on another 59 time series. ETS generally has a lower average MAPE compared to that of ARIMA’s. Although both methods have similar average error percentages at 1 to 3 horizons, ETS performs better in medium and long horizons across 130 time series.<br/>
<img src="/images/model_selection.png?raw=true" width="550">


## 4. Out-of-Sample Forecast
After the best model has been chosen for individual time series, the selected model should produce 18 months out-of-sample forecasts and evaluate its forecast accuracy using genuine forecasts (the actual out-of-sample data available in the M3-competition). 3 error measures (MPE, MAPE, MASE) will be used to determine the model performance, meanwhile benchmarking with 5 other forecast methods including naïve, moving average, simple exponential smoothing, pure ETS and pure ARIMA methods.

### 4.1 Evaluating Forecast Accuracy with Benchmarks by Horizons
The average performance of model selection strategy (the best model chosen between ETS or ARIMA) compared with the 5 benchmarks across 130 time series is given below.<br/>
<img src="/images/accuracy_horizons.png?raw=true" width="550">

In general, the model selection strategy tends to have a lower error rate in terms of average MAPE across all horizons. In other words, it performs better as the accuracy of forecasts is higher compared to that of alternative methods. The auto-ETS method outperforms the others in terms of average MASE, indicating it is a better model when all forecasts of different methods are calculated against a one-step in-sample naive forecast. It is also worth noting that all methods have a negative MPE value, indicating the forecasts are overestimated since the forecast values are always higher than that of the actual values. Even so, the model selection strategy seems to be relatively unbiased among all alternatives as its average MPE value is closest to 0. 

In the aspect of different planning horizons, the proposed model selection strategy has the best forecasting performance in short and medium horizons in terms of MAPE, whereas the auto-ETS performs slightly better in the long horizon. As expected, all methods perform badly in the long horizon since more uncertainties are involved when predicting the far future. On the contrary, they tend to perform better in the medium horizon with regards to a lower value of MAPE. Conversely, the MASE reveals that all methods tend to perform better in short horizon when benchmarking with in-sample naive forecasts.<br/>

### 4.2 A Deeper Look of Forecast Performance by Horizons
To have a deeper look on the strategy’s forecast performance, below visulisations are generated to evaluate the forecast accuracy of each error measure by each horizon.<br/>
<img src="/images/MAPE.png?raw=true" width="650">
<img src="/images/MPE.png?raw=true" width="650">
<img src="/images/MASE.png?raw=true" width="650">

### 4.3 A Deeper Look of Forecast Performance by Categories
A forecasting method might perform differently in predicting various categories of data. In this sense, the below analysis will focus on the MAPE forecast accuracy of the proposed strategy and 5 benchmarks on different types of times series. <br/>
<img src="/images/categories.png?raw=true" width="800">

It is obvious that all methods have a lower error rate in forecasting demographic data which rates between 2.75% to 5.44%. Alternatively, they seem to have a poor forecast accuracy in micro data. It is obvious that demographic data tends to be more “stable” without many fluctuations over the years whereas the micro data contains more uncertainties. 

### 4.4 A Deeper Look of Forecast Performance by Time Series Characteristics 
It is usually difficult to make projections from raw data since a time series can be mixed up with trend and seasonal variables. Different forecasting methods are created with strengths in different aspects. As such, the below analysis is carried out to investigate the forecast accuracy of the proposed strategy and benchmarks with various combinations of characteristics.<br/>
<img src="/images/characteristics.png?raw=true" width="800">

Generally, all models have similar forecast performance in “no trend nor seasonality” data across 130 time series, in particular moving average has the best performance. It is believed that the moving average model is more favorable in predicting fluctuating time series with no trend and seasonality as it has a longer memory nature. In contrast, the proposed strategy outperforms all alternatives in “seasonality plus damped trend” and “trend” data, whereas the auto-ARIMA performs the best in “seasonality and trend” data. Although auto-ETS has a better performance in predicting “damped trend” and “seasonality”, the proposed strategy has very similar performance compared to auto-ETS which accounts for the second-best model among others.

### To learn more about Ivy Wu S.Y., visit her [LinkedIn profile](https://www.linkedin.com/in/ivy-wusumyi)

All rights reserved 2022. All codes are developed and owned by Ivy Wu S.Y..