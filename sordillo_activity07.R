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
                     TropicalV +
                     HydroV,
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
                      TropicalV  = 1,
                      HydroV     = 1),
           interval = "prediction")










