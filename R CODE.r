###PART 01: TIME SERIES FORECASTING####

#Packages
library(imputeTS)
library(forecast)
library(ggplot2)
library(tseries)
library(tidyverse)
library(urca)
library(TSstudio)
library(zoo)
library(remindR)
library(timetk)
library(anomalize)
library(tibbletime)
library(janitor)
library(fpp)
library(seasonalview)

#Data Calling
setwd("~/Documents/Farzana UoL/Semester 2/Forecasting and advanced business analytics/Coursework/Datasets")
data1<-read.csv("PCE.csv", encoding = 'UTF-8')

head(data1)
nrow(data1) 
summary(data1) 
hist(data1$PCE, main = "Distribution of US Personal Consumption Expenditure", 
     ylab = "Frequency of PCE", xlab = "Personal Consumption Expenditure", col = "blue", 
     ylim = c(0,300), xlim = c(0,20000))

#DataType transformation
class(data1$DATE)
data1$DATE<-as.Date(data1$DATE, format ="%d/%m/%Y") 
class(data1$DATE)

#Time series object 
dataTS<-ts(data1$PCE,start=c(1959,1,1),end = c(2021,12,1),frequency=12)
dataTS = window(dataTS, start = c(2005,1), end =c(2021,12))

summary(dataTS) 
start(dataTS)
hend(dataTS)
frequency(dataTS)

autoplot(dataTS)+ggtitle("US Personal Consumption Expenditure 1959-2021")+
  labs(x="Time", y="PCE")

#Generating descriptive statistics
summary(dataTS)

#Dealing with missing data and outliers 
sum(is.na(dataTS))  
tsoutliers(dataTS) 

tsclean<-tsclean(dataTS, replace.missing = TRUE, iterate= 2, lambda = NULL)

autoplot(tsclean)+ggtitle("US Personal Consumption Expenditure 2020-2021 (Clean Dataset)")+
  labs(x="Time", y="PCE")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##Decomposition: Additive and multiplicative

decA<- decompose(tsclean) 
autoplot(decA)+ ggtitle("Additive Decomposition")
checkresiduals(remainder(decA))

decM<- decompose(tsclean, type="multiplicative")
autoplot(decM)+ ggtitle("Multiplicative Decomposition")
checkresiduals(remainder(decM))+ ggtitle("Multiplicative Decomposition Residuals")


#STL Loess 
s<- stl(tsclean, s.window = 'period')
autoplot(s)+ ggtitle("STL Decomposition")
checkresiduals(remainder(s))+ ggtitle("STL Decomposition Residuals")

#X-11

x11<-seas(tsclean, x11="")
autoplot(x11)+ggtitle("X-11 Decomposition")
checkresiduals(remainder(x11))+ ggtitle("X-11 Decomposition Residuals")


#Seasonality check
seasonplot(tsclean) 
ggseasonplot(tsclean, col = rainbow(12), year.labels = TRUE)+ggtitle("Season Plot") 
ggsubseriesplot(tsclean)+ggtitle("Series Plot") 


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

###END OF PRE-PROCESSING######

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Splitting the dataset 

split_pce<-ts_split(tsclean, sample.out = 41) 

training<-split_pce$train 
testing<-split_pce$test 

length(training) 
length(testing) 

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#################################################
#Simple forecasting models 
#################################################
###Simple forecasting methods 

# average method
fit_s1 <- meanf(training, h=41) 
print(summary(fit_s1))
checkresiduals(fit_s1)
accuracy(fit_s1, testing) 
# naive method
fit_s2 <- naive(training, h=151) 
print(summary(fit_s2)) 
checkresiduals(fit_s2)
accuracy(fit_s2, testing) 
# seasonal naive method
fit_s3 <- snaive(training, h=53)	
print(summary(fit_s3)) 
checkresiduals(fit_s3)
accuracy(fit_s3, testing) 

# drift method
fit_s4 <- rwf(training,drift=TRUE,h=41) 
print(summary(fit_s4)) 
checkresiduals(fit_s4)
accuracy(fit_s4, testing) 

test_forecast(actual=tsclean, forecast.obj = fit_s4, test=testing)

#################################################
#Exponential Smoothing Models 
#################################################

#Simple exponential smoothing model

fit_ses1<-ses(training, alpha=0.3, h=41)
print(summary(fit_ses1))  
checkresiduals(fit_ses1)
accuracy(fit_ses1, testing)

fit_ses2<-ses(training, initial = "optimal", h=41) 
print(summary(fit_ses2)) 
checkresiduals(fit_ses2)
accuracy(fit_ses2, testing)

fit_ets<-ets(training) 
print(summary(fit_ets)) 
checkresiduals(fit_ets)
fit_etsfor<-forecast(fit_ets, h=41)
accuracy(fit_etsfor, testing)

fcholt<- holt(training,h=41) 
checkresiduals(fcholt)
fcholt_for<-forecast(fcholt, h=41)
accuracy(fcholt_for, testing)

fchw<-hw(training, h=41) 
print(summary(fchw))
checkresiduals(fchw)
fchw_for<-forecast(fchw, h=41)
accuracy(fchw_for, testing)

test_forecast(actual=tsclean, forecast.obj = fit_etsfor, test=testing)

#################################################
#ARIMA Model 
#################################################

#Differencing

#First differencing
ddeseason=diff(tsclean)
plot(ddeseason,main="First differencing of PCE TS") 
#Second Differencing 
dddeseason=diff(diff(tsclean))
plot(dddeseason, main="Second differencing of PCE TS") 

##Non-stationarity check 

#Augmented Dickey Fuller test 
adf.test(tsclean) 
adf.test(ddeseason)
adf.test(dddeseason) 

#Phillips Perron Test
pp.test(tsclean) 
pp.test(ddeseason)
pp.test(dddeseason) 

#KPSS test 
kpss.test(tsclean) 
kpss.test(ddeseason) 
kpss.test(dddeseason) 

##ACF and PACF 
ggAcf(tsclean)+ggtitle("ACF of PCE ")
ggPacf(tsclean)+ggtitle("PACF of PCE") 

ggAcf(ddeseason)+ggtitle("ACF of PCE (First Differenced)") 
ggPacf(ddeseason)+ggtitle("PACF of PCE (First Differenced)") 

ggAcf(ddeseason)+ggtitle("ACF of PCE (Second Differenced)") 
ggPacf(ddeseason)+ggtitle("PACF of PCE (Second Differenced)")

fit_arima1<-arima(training,order=c(1,1,0)) 
print(summary(fit_arima1)) 
autoplot(fit_arima1) 
checkresiduals(fit_arima1)
fit_arima1for<-forecast(fit_arima1, h=41)
accuracy(fit_arima1for, testing)

fit_arima2<-arima(training,order=c(1,1,1)) 
print(summary(fit_arima2)) 
autoplot(fit_arima2) 
checkresiduals(fit_arima2)
fit_arima2for<-forecast(fit_arima2, h=41)
accuracy(fit_arima2for, testing)

fit_arima7<-arima(training,order=c(1,2,0)) 
print(summary(fit_arima7)) 
autoplot(fit_arima7) 
checkresiduals(fit_arima7)
fit_arima7for<-forecast(fit_arima7, h=41)
accuracy(fit_arima7for, testing)

fit_arima4<-arima(training,order=c(1,2,1)) 
print(summary(fit_arima4)) 
autoplot(fit_arima4) 
checkresiduals(fit_arima4)
fit_arima4for<-forecast(fit_arima4, h=41)
accuracy(fit_arima4for, testing)

fit_arima5<-arima(training,order=c(1,2,2)) 
print(summary(fit_arima5)) 
autoplot(fit_arima5) 
checkresiduals(fit_arima5)
fit_arima5for<-forecast(fit_arima5, h=41)
accuracy(fit_arima5for, testing)

fit_arima6<-arima(training,order=c(2,2,1)) 
print(summary(fit_arima6)) 
autoplot(fit_arima6) 
checkresiduals(fit_arima6)
fit_arima6for<-forecast(fit_arima6, h=41)
accuracy(fit_arima6for, testing)

fit_arima3<-auto.arima(training,stepwise = FALSE, approximation = FALSE,
                       trace=TRUE)  
print(summary(fit_arima3))    
autoplot(fit_arima3) 
checkresiduals(fit_arima3)
fit_arima3for<-forecast(fit_arima3,h=41)
accuracy(fit_arima3for, testing)

fit_sarima<- arima(training,order=c(1,1,1), seasonal = list(order=c(1,0,0)))
print(summary(fit_sarima)) 
autoplot(fit_sarima) 
checkresiduals(fit_sarima)
fit_sarimafor<-forecast(fit_sarima, h=41)
print(fit_sarimafor)
accuracy(fit_sarimafor, testing)


test_forecast(actual=tsclean, forecast.obj = fit_arima3for, test=testing)

#Combined graph of the original data and the predicted data of the best of the 3 types of models 
autoplot(tsclean)+autolayer(fit_s4$mean) + autolayer(fit_etsfor$mean) + autolayer(fit_arima3for$mean)+
  ggtitle("Model Comparison against pre-processed time series") 
autoplot(dataTS)+autolayer(fit_s4$mean) + autolayer(fit_etsfor$mean) + autolayer(fit_arima3for$mean)+ 
  ggtitle("Model Comparison against raw time series")

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Out of sample forecasting for October, 2022

#Generate the optimal fit

finalfit <- rwf(tsclean,drift=TRUE,h=10) 
print(finalfit)
autoplot(finalfit)+
  ggtitle("Forecasted PCE for October 2022 with from Random Walk with Drift")


#One-step ahead rolling forecast 

train<-window(tsclean, end=2018)
test<-window(tsclean, start=2019)

#Drift 
rf_fit_drift<-rwf(train,drift=TRUE)
refit_drift<-snaive(tsclean, model=rf_fit_drfit)
fc_drift<-window(fitted(refit_drift),start=2019)
fc_drift
accuracy(fc_drift,test)

#ETS
rf_fit_ets<-ets(train)
refit_ets<-ets(tsclean, model=rf_fit_ets, use.initial.values=TRUE)
fc_ets<-window(fitted(refit_ets),start=2019)
fc_ets
accuracy(fc_ets,test)

#Auto ARIMA 
rf_fit_arima<-auto.arima(training,stepwise = FALSE, approximation = FALSE,
                         trace=TRUE)
refit_arima<-Arima(tsclean, model=rf_fit_arima, use.initial.values=TRUE)
fc_arima<-window(fitted(refit_arima),start=2019)
fc_arima
accuracy(fc_arima,test)

autoplot(tsclean)+autolayer(fc_drift)+autolayer(fc_ets)+autolayer(fc_arima)+ 
  ggtitle("Rolling forecasting model comparison against pre-processed ts")
autoplot(dataTS)+autolayer(fc_drift)+autolayer(fc_ets)+autolayer(fc_arima)+ 
  ggtitle("Rolling forecasting model comparison against original ts")


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
