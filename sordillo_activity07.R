#Ryan Sordillo
#4/9/2026
#Activity 07
#ENVST325

########Start of Tutorial################

library(dplyr)
library(ggplot2)
library(olsrr)
library(PerformanceAnalytics)

# read in greenhouse gas data from reservoirs
ghg <- read.csv("/cloud/project/activity07/Deemer_GHG_Data.csv")

# log transform methane fluxes
ghg$log.ch4 <- log(ghg$ch4+1)

ghg$log.age <- log(ghg$age)
ghg$log.DIP <- log(ghg$DIP+1)
ghg$log.precip <- log(ghg$precipitation)

unique(ghg$Region)

# binary variable for boreal region
ghg$BorealV <- ifelse(ghg$Region == "Boreal",1,0)
# binary variable for tropical region
ghg$TropicalV <- ifelse(ghg$Region == "Tropical",1,0)

# binary variable for alpine region
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)

# binary variable for known hydropower
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)

# multiple regression
# creates a model object
mod.full <- lm(log.ch4 ~ airTemp+
                 log.age+mean.depth+
                 log.DIP+
                 log.precip+ BorealV, data=ghg) #uses the data argument to specify dataframe
summary(mod.full)


res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)

# qq plot
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)

# shapiro-wilks test
shapiro.test(res.full)

plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:

reg.data <- data.frame(ghg$airTemp,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP,
                       ghg$log.precip)

# make a correlation matrix 
chart.Correlation(reg.data, histogram=TRUE, pch=19)

# run stepwise
full.step <- ols_step_forward_aic(mod.full)
# view table
full.step 

# check full model
full.step$model

# plot AIC over time
plot(full.step )


# prediction with interval for predicting a point
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="prediction")


##Question 3

# Log-transform depth to reduce skew from very deep reservoirs
ghg$log.depth <- log(ghg$mean.depth + 1)

# - Adds HydroV (directly answers the policy question)
# - Adds TropicalV (tropical reservoirs systematically higher CH4 - important given dataset composition)
# - Uses log.depth instead of raw mean.depth

mod.improved <- lm(log.ch4 ~ airTemp +
                     log.age +
                     log.depth +
                     log.DIP +
                     log.precip +
                     BorealV +
                     TropicalV,
                   data = ghg)

summary(mod.improved)


###Assumptions Check for improved model#####

res.imp <- rstandard(mod.improved)
fit.imp <- fitted.values(mod.improved)


qqnorm(res.imp, pch = 19, col = "grey50", main = "Q-Q Plot: Improved Model")
qqline(res.imp)
shapiro.test(res.imp)

plot(fit.imp, res.imp, pch = 19, col = "grey50",
     xlab = "Fitted Values", ylab = "Standardized Residuals",
     main = "Residuals vs Fitted: Improved Model")
abline(h = 0)

ols_vif_tol(mod.improved)

imp.step <- ols_step_forward_aic(mod.improved)
imp.step
plot(imp.step)
summary(imp.step$model)

predict.lm(imp.step$model,
           data.frame(airTemp    = 25,
                      log.age    = log(2),
                      log.depth  = log(16),
                      log.DIP    = log(51),
                      log.precip = log(1500),
                      BorealV    = 0,
                      TropicalV  = 1
                      ),
           interval = "prediction")

######Start of Tutorial#######

ETdat <- read.csv("/cloud/project/activity07/ETdata.csv")

unique(ETdat$crop)

library(lubridate)
library(ggplot2)
library(forecast)
library(dplyr)

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit



# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)


almondTrend <- almond_dec$trend
almondSeason <- almond_dec$seasonal

acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(almond_ts))

almond_y <- na.omit(almond_ts)
model1 <- arima(almond_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1

model4 <- arima(almond_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4

# calculate fit
AR_fit1 <- almond_y - residuals(model1) 
AR_fit4 <- almond_y - residuals(model4)
#plot data
plot(almond_y)
# plot fit
points(AR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")


newAlmond <- forecast(model4)
newAlmond

#make dataframe for plotting
newAlmondF <- data.frame(newAlmond)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newAlmondF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = almond, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(almond$date[1]),newAlmondF$dateF[24])+  # Plotting original data
  geom_line(data = newAlmondF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newAlmondF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")


#####Start of Homework#####

#1.
'The authors of the reservoir greenhouse gas study recommend using the following transformation for
 data:
Use the transformation and design a multiple regression analysis to present to water managers about
the impact of reservoir characteristics on carbon dioxide fluxes. In designing your regression, you
should consider the environmental conditions that impact carbon dioxide fluxes, the availability of
𝐶𝑂2
1/(𝐶𝑂2 + 1000)
data, and the assumptions of ordinary least squares regression. Report results using a regression
table (format in a formal presentation), , and the sample size. Write a paragraph summary and
interpretation of the findings that can be presented to water managers.'
#Transformation to CO2 flux

ghg$transformed_co2 <- 1 / (ghg$co2 + 1000)

sum(!is.na(ghg$co2))
summary(ghg$co2)


#Model for CO2
names(ghg)

mod.co2 <- lm(transformed_co2 ~ log.precip +
                BorealV+
                log.DIP+
                log.depth+
                log.age,
              data=ghg)

summary(mod.co2)

##Checking OLS Assumptions for co2 model##

res.imp <- rstandard(mod.co2)
fit.imp <- fitted.values(mod.co2)


qqnorm(res.imp, pch = 19, col = "grey50", main = "Q-Q Plot: CO2 Model")
qqline(res.imp)
shapiro.test(res.imp)

plot(fit.imp, res.imp, pch = 19, col = "grey50",
     xlab = "Fitted Values", ylab = "Standardized Residuals",
     main = "Residuals vs Fitted: CO2")
abline(h = 0)

ols_vif_tol(mod.co2)

imp.step <- ols_step_forward_aic(mod.co2)
imp.step
plot(imp.step)
summary(imp.step$model)

#3.
'Decompose the evapotranspiration time series for almonds, pistachios, fallow/idle fields, corn, and
table grapes. Evaluate differences in the observations, trends, and seasonality of the data between
the different crops. Write a summary of your evaluation for a water manager that is interested in
examining how irrigation can affect evapotranspiration. The manager also wants to understand what
crops have the greatest water consumption, the timing of high water consumption, and if there are
changes over time. Include plots of your decomposition.'


# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# average fields for each month for pistachios
pistachio <- ETdat %>% # ET data
  filter(crop == "Pistachios") %>% 
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# average fields for each month for fallow
fallow <- ETdat %>% # ET data
  filter(crop == "Fallow/Idle Cropland") %>% 
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# average fields for each month for corn
corn <- ETdat %>% # ET data
  filter(crop == "Corn") %>% # only use corn fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# average fields for each month for grapes
grape <- ETdat %>% # ET data
  filter(crop == "Grapes (Table/Raisin)") %>% 
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# Time Series

almond_ts2 <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

pistachio_ts <- ts(pistachio$ET.in, start = c(2016,1), frequency=12)
fallow_ts <- ts(fallow$ET.in, start = c(2016,1), frequency=12)
corn_ts <- ts(corn$ET.in, start = c(2016,1), frequency=12)
grape_ts <- ts(grape$ET.in, start = c(2016,1), frequency=12)

#Decompose TS
almond_dec <- decompose(almond_ts2)
pistachio_dec <- decompose(pistachio_ts)
fallow_dec <- decompose(fallow_ts)
corn_dec <- decompose(corn_ts)
grape_dec <- decompose(grape_ts)


plot(almond_dec, xlab="Year")
title("Almonds ET Decomposition")

plot(pistachio_dec, xlab="Year")
title("Pistachios ET Decomposition")

plot(fallow_dec, xlab="Year")
title("Fallow/Idle ET Decomposition")

plot(corn_dec, xlab="Year")
title("Corn ET Decomposition")

plot(grape_dec, xlab="Year")
title("Table Grapes ET Decomposition")


#4.
'Design an autoregressive model for pistachios and fallow/idle fields. Forecast future
evapotranspiration for each field so that water managers can include estimates in their planning.
Make a plot that includes historical and forecasted evapotranspiration for the crops to present to the
water manager. Include a brief explanation of your autoregressive models.
'

pistachio_y <- na.omit(pistachio_ts)

p_model1 <- arima(pistachio_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format

p_model4 <- arima(pistachio_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format

AR_fitp1 <- pistachio_y - residuals(p_model1) 
AR_fitp4 <- pistachio_y - residuals(p_model4)

points(AR_fitp1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fitp4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

newPistachio <- forecast(p_model4)

newPistachioF <- data.frame(newPistachio)
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newPistachioF$dateF <- ymd(paste(years,"/",month,"/",1))

#Make plot for pistachios including predictions
ggplot() +
  geom_line(data = pistachio, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(pistachio$date[1]),newPistachioF$dateF[24])+  # Plotting original data
  geom_line(data = newPistachioF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newPistachioF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)", title = "Pistachio AutoCorrelation Model")

##Fallow
fallow_y <- na.omit(fallow_ts)

f_model1 <- arima(fallow_y , # data 
                  order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format

f_model4 <- arima(fallow_y , # data 
                  order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format

AR_fitf1 <- fallow_y - residuals(f_model1) 
AR_fitf4 <- fallow_y - residuals(f_model4)

points(AR_fitf1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fitf4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

newFallow <- forecast(f_model4)

newFallowF <- data.frame(newFallow)
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newFallowF$dateF <- ymd(paste(years,"/",month,"/",1))

#Make plot for pistachios including predictions
ggplot() +
  geom_line(data = fallow, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(fallow$date[1]),newFallowF$dateF[24])+  # Plotting original data
  geom_line(data = newFallowF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newFallowF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)", title = "Fallow AutoCorrelation Model")