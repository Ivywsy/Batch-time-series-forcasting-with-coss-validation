library(forecast)
library(smooth)
library(Mcomp)
library(tseries)
library(parallel)
library(foreach)
library(doSNOW)
library(dplyr)


#input data of 130 time series
raw_data <- lapply(M3, I)
sub_data<- sub("^N", "", names(raw_data))
data_range <- sub_data > 1500 & sub_data <= 2800
data <- raw_data[data_range][grep("0$", names(raw_data[data_range]))] 
rm(raw_data,data_range,sub_data)

#-------Cross-validation, model selection, error calculation for selected model------------------
horizon <- 18
no_models <- 2 #Auto-ETS or Auto-ARIMA
ets_MAPEs <- arima_MAPEs <- data.frame()
ets_APEs <- arima_APEs <- matrix(NA, 130,horizon)
strategy_MAPE <- strategy_MPE <- strategy_MASE <- data.frame()
model_selected<-c()
ts_characteristics <- matrix(NA, 130, 4)

#parallel computing set up
#Detect the number of core and define cluster
cores <- detectCores()
cl <- makeCluster(cores-1, type = "SOCK")
registerDoSNOW(cl)

for (i in 1:130){
    print(i)
    
    y <- data[[i]]$x
    yout <- data[[i]]$xx
    ytlen = data[[i]]$n - horizon
    yt<- head(y, ytlen)
    
    #get the best model recommended by auto.arima()/ets() for each ts
    fit_ets <- ets(yt)
    fit_arima <- auto.arima(yt, method = "CSS")
    
    #check the characteristics of data for further forecast evaluations
    ts_characteristics[i,1] <- data[[i]]$type
    #check the time series characteristics(trend/ seasonality) for further forecast evaluations
    ts_characteristics[i,2] <- fit_ets$component[3] != "N" #Seasonality
    ts_characteristics[i,3] <- fit_ets$component[2] != "N" #Trend
    ts_characteristics[i,4] <- as.logical(fit_ets$component[4])  #Damped
    
    #check if the length of each time series is enough for 12 origins
    #otherwise, apply 6 origins only
    if(ytlen-12 < 25) {
        ytstart<-ytlen-6+1
    }else{
        ytstart<-ytlen-12+1  
    }
    
    origins <- ytstart:ytlen
    #parallel computing for cross validation per each origin
    MAPEs <- foreach(origin=origins, .combine = 'cbind', .packages = 'forecast') %dopar% {
        yt <- head(y, origin)
        yv <- y[(origin+1):(origin+horizon)]
        
        cvfit_ets <- ets(yt, model = paste(fit_ets$components[1],
                                           fit_ets$components[2],
                                           fit_ets$components[3],sep = ""),
                         damped = as.logical(fit_ets$components[4]))
        fcs_ets <- forecast:::forecast(cvfit_ets, h=horizon)$mean
        ets_APE <- 100 *(abs(yv - fcs_ets)/abs(yv))
        ets_MAPE <- mean(ets_APE)
        
        cvfit_arima <- Arima(yt, order=c(fit_arima$arma[1],fit_arima$arma[6],fit_arima$arma[2]),
                             seasonal = c(fit_arima$arma[3],fit_arima$arma[7],fit_arima$arma[4]),
                             method = "CSS")
        fcs_arima <- forecast:::forecast(cvfit_arima, h=horizon)$mean
        arima_APE<- 100 *(abs(yv - fcs_arima)/abs(yv))
        arima_MAPE <- mean(arima_APE)
        
        
        list(ets_MAPE, arima_MAPE,ets_APE,arima_APE)
    }
    
    #Separate ETS and Arima output from each clusters
    ets_MAPEs <- bind_rows(ets_MAPEs,MAPEs[1,])
    arima_MAPEs <- bind_rows(arima_MAPEs, MAPEs[2,])
    
    ets_APEs[i, 1:horizon] <-  colMeans(foreach(j=1:length(MAPEs[3,]), .combine = 'rbind') %dopar% {
        MAPEs[3,][[j]]
    }, na.rm = TRUE)
    
    arima_APEs[i, 1:horizon] <-  colMeans(foreach(j=1:length(MAPEs[4,]), .combine = 'rbind') %dopar% {
        MAPEs[4,][[j]]
    }, na.rm = TRUE)
    
    rm(MAPEs)
    
    #Calculate Naive method's MAE for the scaling denominator of MASE
    #(will be used in the following IF statement)
    naive_fit <- naive(y, h = horizon)
    naive_fcs <- forecast:::forecast(naive_fit)$mean
    
    
    #compare the MAPE of ETS and ARIMA model
    if (rowMeans(ets_MAPEs[i,],na.rm = TRUE) > rowMeans(arima_MAPEs[i,],na.rm = TRUE)){
        fit <- auto.arima(y)
        #Select ARIMA
        model_selected = c(model_selected,"ARIMA")
        
        #create ARIMA forecast for out-sample-data and calculate the accuracy of each horizon
        fcs <- forecast:::forecast(fit, h=horizon)$mean
        strategy_MAPE[i, 1:horizon]<- 100 *(abs(yout - fcs)/abs(yout))
        strategy_MPE[i, 1:horizon]<- 100 *((yout - fcs)/(yout))
        strategy_MASE[i, 1:horizon]<- abs(yout - fcs)/accuracy(naive_fit)[3]
        
    }else{
        fit <- ets(y)
        #Select ETS
        model_selected = c(model_selected,"ETS")
        
        #create ETS forecast for out-sample-data and calculate the accuracy of each horizon
        fcs <- forecast:::forecast(fit, h=horizon)$mean
        strategy_MAPE[i,1:horizon]<- 100 *(abs(yout - fcs)/abs(yout))
        strategy_MPE[i,1:horizon]<- 100 *((yout - fcs)/(yout))
        strategy_MASE[i, 1:horizon]<- abs(yout - fcs)/accuracy(naive_fit)[3]
    }
    
}
stopCluster(cl)


#Remove temp data
rm(fit,fit_ets,fit_arima)



#-----------MAPE plot (cross Validation output)---------------------------------
#ETS generally has better performance forecast compared to ARIMA, resulting 71 ETS
plot(1:18,type = "n",
     main = "Average training set MAPE of 130 time series",
     xlab="Horizon", ylab="MAPE %",
     ylim=c(0,40))
legend("topleft",legend=c("ETS","ARIMA"),col=2:3,lty=1)
lines(1:18, colMeans(ets_APEs), type="l",col=2)
lines(1:18, colMeans(arima_APEs), type="l",col=3)



#----------Benchmark methods----------------------------------------------------

horizon = 18
naive_fcs_MAPE<-sma_fcs_MAPE<-ses_fcs_MAPE<-autoets_fcs_MAPE<-autoarima_fcs_MAPE<-data.frame()
naive_fcs_MPE<-sma_fcs_MPE<-ses_fcs_MPE<-autoets_fcs_MPE<-autoarima_fcs_MPE<-data.frame()
naive_fcs_MASE<-sma_fcs_MASE<- ses_fcs_MASE<-autoets_fcs_MASE<-autoarima_fcs_MASE<-data.frame()

#parallel computing set up
#Detect the number of core and define cluster
cores <- detectCores()
cl <- makeCluster(cores-1, type = "SOCK")
registerDoSNOW(cl)

for (i in 1:130){
    print(i)
    
    y <- data[[i]]$x
    yout <- data[[i]]$xx
    
    #Naive method
    m1 <- forecast:::naive(y, h = horizon)
    #Moving Average method with k = 20
    m2 <- smooth:::sma(y,order = 20, h = horizon)
    #Simple Exponential Smoothing
    m3 <- forecast:::ses(y, h = horizon)
    #auto-ETS
    m4 <- forecast:::ets(y)
    #auto-ARIMA
    m5 <- forecast:::auto.arima(y)
    
    #parallel forecasting
    forecasts <- foreach(j = 1, combine = 'rbind', .packages = 'forecast') %dopar% {
        #Naive method
        naive_fcs <- forecast:::forecast(m1)$mean
        #Moving Average method with k = 20
        sma_fcs <- m2$forecast
        #Simple Exponential Smoothing
        ses_fcs <- forecast:::forecast(m3)$mean
        #auto-ETS
        autoets_fcs <- forecast:::forecast(m4)$mean
        #auto-ARIMA
        autoarima_fcs <- forecast:::forecast(m5)$mean
        
        list(naive_fcs, ses_fcs, ses_fcs, autoets_fcs, autoarima_fcs)
    }
    
    #Naive method
    naive_fcs <- forecasts[[1]][[1]]
    #naive: MAPE
    naive_fcs_MAPE[i, 1:horizon] <- 100 *(abs(yout - naive_fcs)/abs(yout))
    #naive: MPE
    naive_fcs_MPE[i, 1:horizon] <- 100 *((yout - naive_fcs)/(yout))
    #Naive: MASE (scaled by in-sample naive MAE)
    naive_fcs_MASE[i,1:horizon] <- abs(yout - naive_fcs)/accuracy(m1)[3]
    
    #Moving Average method with k = 20
    sma_fcs <- forecasts[[1]][[2]]
    #Moving Average: MAPE
    sma_fcs_MAPE[i, 1:horizon] <- 100 *(abs(yout - sma_fcs)/abs(yout))
    #Moving Average: MPE
    sma_fcs_MPE[i, 1:horizon] <- 100 *((yout - sma_fcs)/(yout))
    #Moving Average: MASE (scaled by in-sample naive MAE)
    sma_fcs_MASE[i,1:horizon] <- abs(yout - sma_fcs)/accuracy(m1)[3]
    
    #Simple Exponential Smoothing
    ses_fcs <- forecasts[[1]][[3]]
    #SES: MAPE
    ses_fcs_MAPE[i, 1:horizon] <- 100 *(abs(yout - ses_fcs)/abs(yout))
    #SES: MPE
    ses_fcs_MPE[i, 1:horizon] <- 100 *((yout - ses_fcs)/(yout))
    #SES: MASE (scaled by in-sample naive MAE)
    ses_fcs_MASE[i, 1:horizon] <- abs(yout - ses_fcs)/accuracy(m1)[3]
    
    #auto-ETS
    autoets_fcs <- forecasts[[1]][[4]]
    #auto-ETS: MAPE
    autoets_fcs_MAPE[i, 1:horizon] <- 100 *(abs(yout - autoets_fcs)/abs(yout))
    #auto-ETS: MPE
    autoets_fcs_MPE[i, 1:horizon] <- 100 *((yout - autoets_fcs)/(yout))
    #auto-ETS: MASE (scaled by in-sample naive MAE)
    autoets_fcs_MASE[i, 1:horizon] <- abs(yout - autoets_fcs)/accuracy(m1)[3]
    
    #auto-ARIMA
    autoarima_fcs <- forecasts[[1]][[5]]
    #auto-ARIMA: MAPE
    autoarima_fcs_MAPE[i, 1:horizon] <- 100 *(abs(yout - autoarima_fcs)/abs(yout))
    #auto-ARIMA: MPE
    autoarima_fcs_MPE[i, 1:horizon] <- 100 *((yout - autoarima_fcs)/(yout))
    #auto-ARIMA: MASE (scaled by in-sample naive MAE)
    autoarima_fcs_MASE[i, 1:horizon] <- abs(yout - autoarima_fcs)/accuracy(m1)[3]
    
    rm(forecasts)
}

stopCluster(cl)

rm(m1,m2,m3,m4,m5)

#-----------------------------------------------------------------------------
#MAPE
plot(1:18,type = "n",
     main = "Average MAPE of 130 time series",
     xlab="Horizon", ylab="MAPE",
     ylim=c(0,40))
axis(side = 1, at = c(1:18))
legend("topleft",legend=c("Proposed strategy","Naive","Simple Exponential Smoothing","Moving Average","Auto ETS","Auto ARIMA"),col=c(1:4,7,8),lty=1,cex = 0.8)
lines(1:18, colMeans(strategy_MAPE), type="l",col=1)
lines(1:18, colMeans(naive_fcs_MAPE), type="l",col=2)
lines(1:18, colMeans(ses_fcs_MAPE), type="l",col=3)
lines(1:18, colMeans(sma_fcs_MAPE), type="l",col=4)
lines(1:18, colMeans(autoets_fcs_MAPE), type="l",col=7)
lines(1:18, colMeans(autoarima_fcs_MAPE), type="l",col=8)


#MASE
plot(1:18,type = "n",
     main = "Average MASE of 130 time series",
     xlab="Horizon", ylab="MASE",
     ylim=c(0,5))
axis(side = 1, at = c(1:18))
legend("topleft",legend=c("Proposed strategy","Naive","Simple Exponential Smoothing","Moving Average","Auto ETS","Auto ARIMA"),col=c(1:4,7,8),lty=1,cex = 0.8)
lines(1:18, colMeans(strategy_MASE), type="l",col=1)
lines(1:18, colMeans(naive_fcs_MASE), type="l",col=2)
lines(1:18, colMeans(ses_fcs_MASE), type="l",col=3)
lines(1:18, colMeans(sma_fcs_MASE), type="l",col=4)
lines(1:18, colMeans(autoets_fcs_MASE), type="l",col=7)
lines(1:18, colMeans(autoarima_fcs_MASE), type="l",col=8)

#MPE
plot(1:18,type = "n",
     main = "Average MPE of 130 time series",
     xlab="Horizon", ylab="MPE",
     ylim=c(-30,20))
axis(side = 1, at = c(1:18))
legend("topleft",legend=c("Proposed strategy","Naive","Simple Exponential Smoothing","Moving Average","Auto ETS","Auto ARIMA"),col=c(1:4,7,8),lty=1,cex = 0.8)
lines(1:18, colMeans(strategy_MPE), type="l",col=1)
lines(1:18, colMeans(naive_fcs_MPE), type="l",col=2)
lines(1:18, colMeans(ses_fcs_MPE), type="l",col=3)
lines(1:18, colMeans(sma_fcs_MPE), type="l",col=4)
lines(1:18, colMeans(autoets_fcs_MPE), type="l",col=7)
lines(1:18, colMeans(autoarima_fcs_MPE), type="l",col=8)
abline(h = 0, col = "grey")


#----------Evaluation on 6 categories-----------------------------------------
#Filter the 130 ts into 6 categories and evaluate the average MAPE for each forecasting methods

categories <- unique(ts_characteristics[,1]) 
categories_table <- matrix(NA,nrow = 6,ncol = 6)
rownames(categories_table) <- c("Proposed Strategy","Naive","Simple Exponential Smoothing","Moving Average", "Auto ETS", "Auto ARIMA")
colnames(categories_table) <- categories
i <- 1
for (category in categories){
    filtered_strategy_MAPE <- filter(strategy_MAPE, ts_characteristics[,1] == category)
    filtered_naive_MAPE <- filter(naive_fcs_MAPE, ts_characteristics[,1] == category) 
    filtered_ses_MAPE <- filter(ses_fcs_MAPE, ts_characteristics[,1] == category) 
    filtered_sma_MAPE <- filter(sma_fcs_MAPE, ts_characteristics[,1] == category)
    filtered_autoets_MAPE <- filter(autoets_fcs_MAPE, ts_characteristics[,1] == category)
    filtered_autoarima_MAPE <- filter(autoarima_fcs_MAPE, ts_characteristics[,1] == category)
    
    categories_table[1,i] <-mean(rowMeans(filtered_strategy_MAPE,na.rm = TRUE))
    categories_table[2,i] <-mean(rowMeans(filtered_naive_MAPE))
    categories_table[3,i] <-mean(rowMeans(filtered_ses_MAPE))
    categories_table[4,i] <-mean(rowMeans(filtered_sma_MAPE))
    categories_table[5,i] <-mean(rowMeans(filtered_autoets_MAPE))
    categories_table[6,i] <-mean(rowMeans(filtered_autoarima_MAPE))
    
    i <- i + 1
}

rm(filtered_strategy_MAPE, filtered_naive_MAPE,filtered_ses_MAPE,filtered_sma_MAPE,filtered_autoets_MAPE,filtered_autoarima_MAPE)

#----------Evaluation on 6 ts characteristics-----------------------------------------
#Filter the 130 ts into 6 combinations of trend/seasonality 
#evaluate the average MAPE for each forecasting methods
conditions <- list(c(TRUE, TRUE, FALSE),    #Seasonality+Trend
                   c(TRUE, TRUE, TRUE),     #Seasonality+DTrend
                   c(FALSE, TRUE, TRUE),    #DTrend
                   c(FALSE, TRUE, FALSE),   #Trend
                   c(TRUE, FALSE, FALSE),   #Seasonality
                   c(FALSE, FALSE, FALSE))  #No Trend/Seasonality

ts_characteristic_table <- matrix(NA,nrow = 6,ncol = 6)
rownames(ts_characteristic_table) <- c("Proposed Strategy","Naive",
                                       "Simple Exponential Smoothing",
                                       "Moving Average",
                                       "Auto ETS", "Auto ARIMA")
colnames(ts_characteristic_table) <- c("Seasonality+Trend",
                                       "Seasonality+DTrend",
                                       "DTrend",
                                       "Trend",
                                       "Seasonality",
                                       "No Trend/Seasonality")
i <- 1
for (condition in conditions){
    filter_index <-
        ts_characteristics[,2] == condition[1] & 
        ts_characteristics[,3] == condition[2] & 
        ts_characteristics[,4] == condition[3]
    filtered_strategy_MAPE <- filter(strategy_MAPE, filter_index)
    filtered_naive_MAPE <- filter(naive_fcs_MAPE, filter_index) 
    filtered_ses_MAPE <- filter(ses_fcs_MAPE, filter_index) 
    filtered_sma_MAPE <- filter(sma_fcs_MAPE, filter_index)
    filtered_autoets_MAPE <- filter(autoets_fcs_MAPE, filter_index)
    filtered_autoarima_MAPE <- filter(autoarima_fcs_MAPE, filter_index)
    
    ts_characteristic_table[1,i] <-mean(rowMeans(filtered_strategy_MAPE,na.rm = TRUE))
    ts_characteristic_table[2,i] <-mean(rowMeans(filtered_naive_MAPE))
    ts_characteristic_table[3,i] <-mean(rowMeans(filtered_ses_MAPE))
    ts_characteristic_table[4,i] <-mean(rowMeans(filtered_sma_MAPE))
    ts_characteristic_table[5,i] <-mean(rowMeans(filtered_autoets_MAPE))
    ts_characteristic_table[6,i] <-mean(rowMeans(filtered_autoarima_MAPE))
    i <- i + 1
}

rm(filtered_strategy_MAPE,filtered_naive_MAPE,filtered_ses_MAPE,filtered_sma_MAPE,filtered_autoets_MAPE,filtered_autoarima_MAPE)

