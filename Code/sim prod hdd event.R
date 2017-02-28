##################################################
##                                              ##
##  Title:      sim prod hdd event.R            ##
##  Author:     Andrew Bernath, Cadmus Group    ##
##  Created:    7/11/2016                       ##
##  Description:                                ##
##    Simulation of facility consumption where  ##
##    true model has the form:                  ##
##      kWh ~ int + prod + hdd + event +        ##
##            year 1 + prod*year1 + hdd*year1   ##
##                                              ##
##################################################

##  Clear variables
rm(list=ls())
options(scipen=999)

##  Include packages
library(plyr)
library(dplyr)
library(xlsx)
library(ggplot2)


##  Set project folder
projPath     <- "C:/Users/andrew.bernath/Documents/JSM Abstract/Simulation"
stopifnot(file.exists(projPath))

simData <- data.frame(read.csv(file.path(projPath,"data","Data for Sim.csv"), 
                               header=T,
                               stringsAsFactors=F))
head(simData)

##  Fix numbers
for(ii in 1:length(simData)) {
  simData[,ii] <- as.numeric(simData[,ii])
}



###############  CASE 3: TWO VARIABLES WITH YEAR 1 EVENT  ###############

##  Monthly data with two baseline years
n.period   <- 36

##  Confidence level and Z-Stat for testing true savings capture
conf   <- 0.80
z.stat <- qnorm(1-(1-conf)/2,0,1)

##  (PRODUCTION, HDD, & EVENT)  ##
##  Set parameters for model
coef.intercept <- 2532237
coef.prod      <- 4.5511
coef.hdd       <- 1383
coef.event.Y1  <- 225866
coef.year1     <- -270392
coef.prod.Y1   <- 0.1523
coef.hdd.Y1    <- -459.0

##  Facility data input
beta.true <- c(coef.intercept, coef.prod, coef.hdd, coef.event.Y1, 
               coef.year1, coef.prod.Y1, coef.hdd.Y1)
X         <- data.frame(intercept = rep(1, n.period), 
                        prod      = simData$prod,
                        hdd       = simData$HDD50,
                        event_Y1  = simData$event_Y1,
                        year1     = simData$year1, 
                        prod_Y1   = simData$prod_Y1,
                        hdd_Y1    = simData$HDD_Y1
)
X.pre     <- X[ 1:24,]
X.post    <- X[25:36,]
    


###############################################################################
#########################  BEGIN SIMULATION FUNCTION  #########################

##  Input model error
##  Output vector containing:
    # Savings estimates
    # Bias of savings estimates
    # Standard errors of savings estimates
    # CV for savings estimates
    # Model selection criteria (MSE, AIC, BIC, Adj. R^2)


savingsSim.func <- function(modError.in) {
  
  ######  TRUE MODEL  ######
  
  ##  Create vector of model errors
  mod.epsilon   <- rnorm(n.period, 0, modError.in)
  
  ##  Create BL and SEM response vectors with true parameters
    kWh.bl   <- (as.matrix(X[,1:3]) %*% beta.true[1:3]) + mod.epsilon
    kWh.meas <- (as.matrix(X)       %*% beta.true)      + mod.epsilon 
    true.sav <- sum(kWh.bl[25:36,]-(kWh.meas[25:36,]-coef.event.Y1*X$event_Y1[25:36]))
    consump  <- sum(kWh.meas[25:36,])
    
    
    
######  FORECAST MODEL  ######
    
    ##  Estimate model
    FC.mod       <- lm(kWh.meas[1:24] ~ prod + hdd, data=X.pre)
    FC.mod.sum   <- summary(FC.mod)
    FC.mod.b.hat <- FC.mod$coeff
    
    ##  Model selection criteria
    FC.mod.MSE   <- mean(residuals(FC.mod)^2)
    FC.mod.RMSE  <- sqrt(FC.mod.MSE)
    FC.mod.adjR2 <- FC.mod.sum$adj.r.squared
    FC.mod.AIC   <- AIC(FC.mod)
    FC.mod.BIC   <- BIC(FC.mod)
    
    ##  Compute savings estimate
    FC.mod.preds <- as.matrix(X.post[,1:3]) %*% FC.mod.b.hat
    FC.mod.sav   <- sum(FC.mod.preds - kWh.meas[25:36])

    ##  Estimate SE(savings)
    ##  Identity matrix
    ident.mat   <- diag(nrow(X.post))
    XX.inv      <- solve(t(data.matrix(X.pre[,1:3])) %*% data.matrix(X.pre[,1:3]))
    post.mats   <- data.matrix(X.post[,1:3]) %*% XX.inv %*% t(data.matrix(X.post[,1:3]))
    FC.mod.seSav <- sqrt(FC.mod.MSE*sum(post.mats+ident.mat))
    FC.mod.CV    <- abs(FC.mod.seSav/FC.mod.sav)
    
    FC.mod.incl <- 1
    if( true.sav < (FC.mod.sav-(z.stat*FC.mod.seSav)) |
        true.sav > (FC.mod.sav+(z.stat*FC.mod.seSav)) ) FC.mod.incl <- 0
    
    
######  SIMPLE PRE-POST MODEL  ######
    
    ##  Estimate model
    SPP.mod     <- lm(kWh.meas ~ X$prod + X$hdd + X$event_Y1 + X$year1)
    SPP.mod.sum <- summary(SPP.mod)
    
    ##  Estimate savings & SE(savings)
    SPP.mod.sav   <- -1*SPP.mod$coeff[5]*sum(X.post$year1)
    SPP.mod.seSav <- SPP.mod.sum$coeff[5,2]*sum(X.post$year1)
    SPP.mod.CV    <- abs(SPP.mod.seSav/SPP.mod.sav)
    
    ##  Model selection criteria
    SPP.mod.MSE   <- mean(residuals(SPP.mod)^2)
    SPP.mod.RMSE  <- sqrt(SPP.mod.MSE)
    SPP.mod.adjR2 <- SPP.mod.sum$adj.r.squared
    SPP.mod.AIC   <- AIC(SPP.mod)
    SPP.mod.BIC   <- BIC(SPP.mod)
    
    SPP.mod.incl <- 1
    
    if( true.sav < (SPP.mod.sav-(z.stat*SPP.mod.seSav)) |
        true.sav > (SPP.mod.sav+(z.stat*SPP.mod.seSav)) ) SPP.mod.incl <- 0
    
    
    
######  FULLY SPECIFIED PRE-POST MODEL  ######
    
    ##  Estimate model
    FPP.mod     <- lm(kWh.meas ~ X$prod + X$hdd + X$event_Y1 +
                        X$year1 + X$prod_Y1 + X$hdd_Y1)
    FPP.mod.sum <- summary(FPP.mod)
    
    ##  Estimate savings & SE(savings)
    FPP.mod.sav   <- -1*(FPP.mod$coeff[5]*sum(X.post$year1) + 
                           FPP.mod$coeff[6]*sum(X.post$prod_Y1) +
                           FPP.mod$coeff[7]*sum(X.post$hdd_Y1)
                        )
    FPP.cov.out   <- vcov(FPP.mod)
    FPP.mod.seSav <- sqrt(FPP.cov.out[5,5]*sum(X.post$year1)^2 +
                            FPP.cov.out[6,6]*sum(X.post$prod_Y1)^2 +
                            FPP.cov.out[7,7]*sum(X.post$hdd_Y1)^2 +
                            2*FPP.cov.out[5,6]*sum(X.post$year1)*sum(X.post$prod_Y1) +
                            2*FPP.cov.out[5,7]*sum(X.post$year1)*sum(X.post$hdd_Y1) +
                            2*FPP.cov.out[6,7]*sum(X.post$prod_Y1)*sum(X.post$hdd_Y1)
                          )
    FPP.mod.CV    <- abs(FPP.mod.seSav/FPP.mod.sav)
    
    ##  Model selection criteria
    FPP.mod.MSE   <- mean(residuals(FPP.mod)^2)
    FPP.mod.RMSE  <- sqrt(FPP.mod.MSE)
    FPP.mod.adjR2 <- FPP.mod.sum$adj.r.squared
    FPP.mod.AIC   <- AIC(FPP.mod)
    FPP.mod.BIC   <- BIC(FPP.mod)
    
    FPP.mod.incl <- 1
    if( true.sav < (FPP.mod.sav-(z.stat*FPP.mod.seSav)) |
        true.sav > (FPP.mod.sav+(z.stat*FPP.mod.seSav)) ) FPP.mod.incl <- 0
    
    
    
    ######  Create output vector
    
    sim.outVect <- data.frame(consump, true.sav, 
                              FC.mod.sav, SPP.mod.sav, FPP.mod.sav,
                              FC.mod.seSav, SPP.mod.seSav, FPP.mod.seSav,
                              FC.mod.CV, SPP.mod.CV, FPP.mod.CV,
                              FC.mod.MSE, SPP.mod.MSE, FPP.mod.MSE,
                              FC.mod.adjR2, SPP.mod.adjR2, FPP.mod.adjR2,
                              FC.mod.AIC, SPP.mod.AIC, FPP.mod.AIC,
                              FC.mod.BIC, SPP.mod.BIC, FPP.mod.BIC,
                              FC.mod.incl, SPP.mod.incl, FPP.mod.incl,
                              row.names=ii)
    
    return(sim.outVect)
    
}


##  Repeat the simulation 10,000 times for each model error input
N.sim <- 10000
modError <- 200000

for(ii in 1:N.sim) {
  if(ii==1) overspec.sim <- savingsSim.func(modError)
  else      overspec.sim <- rbind(overspec.sim, savingsSim.func(modError))
}   

##  Summarize results of sim
simSummary <- data.frame(t(colMeans(overspec.sim)))

##  Compute bias of mean savings estimate
simSummary$FC.mod.savBias  <- simSummary$FC.mod.sav  - simSummary$true.sav
simSummary$SPP.mod.savBias <- simSummary$SPP.mod.sav - simSummary$true.sav
simSummary$FPP.mod.savBias <- simSummary$FPP.mod.sav - simSummary$true.sav

simSummary 

write.xlsx(simSummary, file.path(projPath,"data","simData case3.xlsx"), 
           sheetName="Sim Out", col.names=T, append=F)
write.xlsx(overspec.sim, file.path(projPath,"data","simData case3.xlsx"), 
           sheetName="Sim Output - prod-hdd-event", col.names=T, append=T)
