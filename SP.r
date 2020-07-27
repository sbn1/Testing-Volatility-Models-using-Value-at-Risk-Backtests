rm(list = ls())
library(rugarch)
library(tseries)
library(timeSeries)
library(forecast)
library(quantmod)

#### Getting Data on S&P500 for 15 years & converting to log returns
getSymbols("SPY" , src="yahoo", from = "2001-01-01", to = "2016-01-01")
SPY = SPY[,"SPY.Adjusted"]
SPY_returns = diff(log(SPY))*100
SPY_returns = SPY_returns[-1,]
rm(SPY)


#Ratio of training vs testing data
sample_split = 0.8

##### Specifying Combinations of Model Orders for mean and volatility, Error Types
ARMA_order = list(c(0,0), c(1,0), c(0,1), c(1,1), c(1,2), c(2,1), c(2,2), c(2,3), c(3,2), c(3,3))
GARCH_order = list(c(1,1), c(1,2), c(2,1), c(2,2))
error_type = c("norm", "std", "sstd")
criterion = data.frame(NA)
criterion[,1:9] = NA
criterion[1:(length(ARMA_order)*length(GARCH_order)*length(error_type)),] = NA
names(criterion) = c("ARMA Order p","ARMA Order q", "GARCH Order p", "GARCH Order q", "Error Type","AIC", "BIC", "SHI", "HQ")

####GJR-GARCH Loop for best fit
count = 1 #For cirtierion table
for (arma in ARMA_order) {
  for (garch in GARCH_order) {
    for (error in error_type) {
      criterion[count, 1:5] = c(arma, garch, error)
      spec = ugarchspec(variance.model = list(model = "gjrGARCH",
                                             garchOrder = garch),
                       mean.model = list(armaOrder = arma),
                       distribution.model = error)
      fit = ugarchfit(spec = spec, 
                      data = SPY_returns, 
                      out.sample = floor((1-sample_split)*length(SPY_returns))
                      )
     criterion[count, 6:9] = as.numeric(infocriteria(fit))
     print(count)
     count = count + 1
    }
  }
}
rm(arma, count, error, garch, fit, spec)
rm(ARMA_order, GARCH_order, error_type)
###Based on three criterion, ARMA(1,1)-GARCH(2,2) with skewed student errors seems to be the best model
## Refit based on this
sp_gjr_spec = ugarchspec(variance.model = list(model = "gjrGARCH",
                                               garchOrder = c(2,2)),
                         mean.model = list(armaOrder = c(1,1)),
                         distribution.model = "sstd")
sp_gjr_fit = ugarchfit(spec = sp_gjr_spec, 
                       data = SPY_returns, 
                       out.sample = floor((1-sample_split)*length(SPY_returns))
)
#Extract Fitted Values, VaR
gjr_var1 = quantile(sp_gjr_fit, 0.01)
gjr_fitted = fitted(sp_gjr_fit)
gjr_fitted_sigma = sigma(sp_gjr_fit)
sort(index(gjr_fitted), decreasing = TRUE)[1]
####Forecasting
gjr_forecast = ugarchforecast(sp_gjr_fit, n.ahead = 1, n.roll = (floor((1-sample_split)*length(SPY_returns)-1)))
forecast_time = sort(index(SPY_returns), decreasing = TRUE)[(floor((1-sample_split)*length(SPY_returns))):1]
forecast_time = as.Date(forecast_time, format = "%Y-%m-%d")

gjr_forecast_var = xts(t(quantile(gjr_forecast, 0.01)), order.by = forecast_time)
gjr_forecast_mu = xts(t(fitted(gjr_forecast)), order.by = forecast_time)
gjr_forecast_sigma = xts(t(sigma(gjr_forecast)), order.by = forecast_time)

### 1% VaR forecast 95% Confidence Intervall
gjr_sim_CI = data.frame(NA)
gjr_sim_CI[1:length(forecast_time),] = NA
gjr_sim_CI[,1:1000] = NA
for (j in 1:1000) {
    gjr_simulation = ugarchsim(sp_gjr_fit, n.start = 0, n.sim = length(forecast_time), m.sim = 1)
    gjr_sim_CI[,j] = as.numeric(quantile(gjr_simulation, 0.01))
    print(j)
}
rm(j, gjr_simulation)
gjr_sm_CI95 = apply(gjr_sim_CI, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
gjr_sim_CI = xts(t(gjr_sm_CI95), order.by = forecast_time)
rm(gjr_sm_CI95)

##Plotting
plot(SPY_returns, col = "gray80", ylim = c(-20,20), main = "S&P500 ARMA(1,1)-GJR-GARCH(2,2)")
lines(gjr_var1, col = "orangered", type = "l")
points(SPY_returns[SPY_returns < gjr_var1], col = "red4", pch = 19)
lines(gjr_fitted, col = "black")
lines(gjr_fitted_sigma, col = "gray40")
lines(gjr_forecast_var, col = "red1")
points((SPY_returns[(length(SPY_returns) - length(forecast_time) + 1):length(SPY_returns)])[(SPY_returns[(length(SPY_returns) - length(forecast_time) + 1):length(SPY_returns)])< gjr_forecast_var], col = "red4", pch = 15)
lines(gjr_forecast_mu, col = "blue")
lines(gjr_forecast_sigma, col = "green")
lines(gjr_sim_CI, col = "orange")
addLegend("topleft", 
       legend.names = c("Actual Return", "Fitted Return", "Fitted Standard Deviation", "1% VaR", "1% VaR Violations",
                  "Forecast Return", "Forecast 1% VaR", "Forecast Standard Deviation", "Forecast 1% VaR Violation",
                  "Forecast 1% VaR 95% CI"),
       col = c("gray80", "black", "gray40", "orangered", "red4",
               "blue", "red1", "green", "red4", "orange"),
       lty = c(1,1,1,1,NA,1,1,1,NA,1),
       pch = c(NA,NA,NA,NA,19,NA,NA,NA,15,NA),
       bg = "white", 
       bty = "o",
       ncol =2)

