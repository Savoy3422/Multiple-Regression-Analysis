### Load libraries ###
library(car) # for outliertest
library(MASS) # for variable selection (stepAIC)
library(olsrr) # for diagnoising collinearity function

#################
### Load data ###
#################
basedir = "/Volumes/qac/prj/webb/20180507-correlation-analyses/analysis/"
sledai = read.csv(paste0(basedir, "SLEDAI_combined.csv"), header=TRUE, sep=",")
sledai_N = sledai[21:39,] # obtain neutrophils group
sledai_p = sledai[1:20,] # obtain pDCs group

#################################################
### Neutrophils Data exploration/ Diagnostics ###
#################################################
### Checking to make sure there is a linear relationship between predictor variables and outcome variable ###
plot(sledai_N$X.ARID3a., sledai_N$SLEDAI) # checking for non-linearity, constant variance, and outliers
plot(sledai_N$X.IFNA., sledai_N$SLEDAI) # checking for non-linearity, constant variance, and outliers

### Assessing multicollinearity with corrlation matrix and VIF ###
cor(sledai_N[,3:5]) # correlation between ARID3a and IFNa
N_fit <- lm(sledai_N$SLEDAI ~ sledai_N$X.ARID3a. + sledai_N$X.IFNA. + sledai_N$X.ARID3a.* sledai_N$X.IFNA. ) # proposed model
ols_coll_diag(N_fit) # diagnoise collinearity (VIF above 4, check, VIF above 10 correction needed)
N_fit <- lm(SLEDAI ~ X.ARID3a. + X.IFNA., data=sledai_N) # new model in attempts to remove multicollinearity
ols_coll_diag(N_fit) # diagnoise collinearity (VIF above 4, check, VIF above 10 correction needed)

### Checking for outliers ###
d.N_fit = cooks.distance(N_fit)
cutoff <- 4/((nrow(sledai_N)-length(N_fit$coefficients)-2)) 
plot(N_fit, which=4, cook.levels=cutoff)
influencePlot(N_fit,	id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )

### Checking Multivariate normality ###
plot(sledai_N$SLEDAI, resid(N_fit),  ylab="Residuals", xlab="SLEDAI") # Plot residuals versus outcome variable
abline(0, 0)  
plot(N_fit$fitted.values, N_fit$residuals,  ylab="Residuals", xlab="Fitted") # Plot residuals versus fitted values
abline(0, 0)  
plot(sledai_N$X.ARID3a., N_fit$residuals,  ylab="Residuals", xlab="ARID3a") # Plot residuals versus ARID3a
abline(0, 0)  
plot(sledai_N$X.IFNA., N_fit$residuals,  ylab="Residuals", xlab="IFNa") # Plot residuals versus IFNa
abline(0, 0)  
qqPlot(N_fit, main="QQ Plot") #qq plot for studentized resid **Check normality
ks.test(N_fit$residuals, "pnorm")

### Checking Homoscedasticity ###
plot(N_fit$fitted.values,rstandard(N_fit)) # equally distributed across all values of the independent variables
abline(0,0)

########################
### Model Selections ###
########################
N_fit <- lm(sledai_N$SLEDAI ~ sledai_N$X.ARID3a. + sledai_N$X.IFNA.)
summary(N_fit) # show results
stepAIC(N_fit, direction="both")
N2_fit = lm(SLEDAI ~ X.ARID3a., data=sledai_N)
summary(N2_fit)


##########################################
### pDCs Data exploration/ Diagnostics ###
##########################################
### Checking to make sure there is a linear relationship between predictor variables and outcome variable ###
plot(sledai_p$X.ARID3a., sledai_p$SLEDAI) # checking for non-linearity, constant variance, and outliers
plot(sledai_p$X.IFNA., sledai_p$SLEDAI) # checking for non-linearity, constant variance, and outliers

### Assessing multicollinearity with corrlation matrix and VIF ###
cor(sledai_p[,3:5]) # correlation between ARID3a and IFNa
p_fit <- lm(SLEDAI ~ X.ARID3a. + X.IFNA. + X.ARID3a.* X.IFNA., data=sledai_p ) # proposed model
ols_coll_diag(p_fit) # diagnoise collinearity (VIF above 4, check, VIF above 10 correction needed)
p_fit <- lm(SLEDAI ~ X.ARID3a. + X.IFNA., data=sledai_p) # new model in attempts to remove multicollinearity
ols_coll_diag(p_fit) # diagnoise collinearity (VIF above 4, check, VIF above 10 correction needed)

### Checking for outliers ###
d.p_fit = cooks.distance(p_fit)
cutoff <- 4/((nrow(sledai_p)-length(p_fit$coefficients)-2)) 
plot(p_fit, which=4, cook.levels=cutoff)
influencePlot(p_fit,	id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )

### Checking Multivariate normality ###
plot(p_fit$fitted.values, p_fit$residuals,  ylab="Residuals", xlab="Fitted") # Plot residuals versus fitted values
abline(0, 0)  
plot(sledai_p$X.ARID3a., p_fit$residuals,  ylab="Residuals", xlab="ARID3a") # Plot residuals versus ARID3a
abline(0, 0)  
plot(sledai_p$X.IFNA., p_fit$residuals,  ylab="Residuals", xlab="IFNa") # Plot residuals versus IFNa
abline(0, 0)  
qqPlot(p_fit, main="QQ Plot") #qq plot for studentized resid **Check normality
ks.test(p_fit$residuals, "pnorm")

### Checking Homoscedasticity ###
plot(p_fit$fitted.values,rstandard(p_fit)) # equally distributed across all values of the independent variables
abline(0,0)

########################
### Model Selections ###
########################
p_fit <- lm(SLEDAI ~ X.ARID3a. + X.IFNA., data=sledai_p)
summary(p_fit) # show results
stepAIC(p_fit, direction="both")
p2_fit = lm(SLEDAI ~ X.IFNA., data=sledai_p)
summary(p2_fit)


