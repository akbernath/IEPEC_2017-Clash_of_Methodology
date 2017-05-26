##############################################################################
## 
##  Title:      Case 1 - Correctly Specified.R
##  Author:     Andrew Bernath & Maggie Buffum, Cadmus Group
##  Created:    2/28/2017
##  Description:
##    Simulation of facility consumption where
##    true model has the form: 
##      kWh ~ int + prod1 + prod2 + event_pre + HDD + HDD*prod1 +
##              + event_post + prog_ind + prog*prod1 + prog*prod2 + 
##              + prog*HDD + prog*HDD*prod1 + ERROR
## 
##############################################################################

##############################################################################
##### Prepare workspace and format sim data
##############################################################################

  ################################
  ##### Prep Workspace
  ##  Clear variables
  rm(list=ls())
  options(scipen=999)
  
  ##  Include packages
  library(plyr)
  library(dplyr)
  library(xlsx)
  library(forecast)
  library(ggplot2)
  
  ##  Set project folder)
  if("Margaret.Buffum" %in% dir(file.path("C:","Users"))) projPath <- file.path("C:","Users","Margaret.Buffum" ,"Documents","GitHub","IEPEC SEM")
  if("andrew.bernath"  %in% dir(file.path("C:","Users"))) projPath <- file.path("C:","Users","andrew.bernath"  ,"Documents","Git","projects","IEPEC_2017-Clash_of_Methodology")

  stopifnot(file.exists(projPath))
  
  ################################
  ##### Read/Format Data
  ## Read in SIM data
  simData0 <- data.frame(read.xlsx(file.path(projPath,"Data","simData_complex.xlsx")
                                   , sheetName="Final Facility Data"
                                   , header=T
                                   , stringsAsFactors=F))
  simData <- simData0[which(!is.na(simData0$Start_Date)),]
  head(simData)
  
  ##  Read in parameters for model
  simCoeff <- data.frame(read.xlsx(file.path(projPath,"Data","simData_complex.xlsx")
                                   , sheetName="Step 2 - Final Model Spec"
                                   , header=T
                                   , stringsAsFactors=F
                                   , rowIndex=c(16:29)
                                   , colIndex=c(1:2)))
  names(simCoeff) <- c("Coefficient", "Value")
  head(simCoeff)
  simCoeff$Coefficient

  ##  Format numbers
  for(ii in 2:length(simData)) {
    simData[,ii] <- as.numeric(simData[,ii])
  }

  
##############################################################################
##### Create data??
##############################################################################
  
###############  CASE 3: TWO VARIABLES WITH YEAR 1 EVENT  ###############

##  Daily data with two baseline years and two program years
# n.period   <- 365

##  Confidence level and Z-Stat for testing true savings capture
conf   <- 0.80
z.stat <- qnorm(1-(1-conf)/2,0,1)


##  Facility data input
beta.true <- simCoeff$Value

X.full    <- data.frame( rep(1, nrow(simData))
                        , simData$prod1
                        , simData$prod2
                        , simData$noProd_ind
                        , simData$event1_pre
                        , simData$HDD55
                        , simData$prod1_x_HDD55
                        , simData$event2_post
                        , simData$prog_ind
                        , simData$prog_x_prod1
                        , simData$prog_x_prod2
                        , simData$prog_x_noProd
                        , simData$prog_x_prod1_x_HDD55
                        , simData$CDD65
                        , simData$prog_x_CDD65
)
names(X.full) <- c(simCoeff$Coefficient, "CDD65", "prog_x_CDD65")



###############################################################################
#########################  BEGIN SIMULATION FUNCTION  #########################

##  Input model error
##  Output vector containing:
    # Savings estimates
    # Bias of savings estimates
    # Standard errors of savings estimates
    # CV for savings estimates
    # Model selection criteria (MSE, AIC, BIC, Adj. R^2)

modError.in <- 0.015*mean(simData$sim_kWh[which(simData$prog_ind == 0)])
df.in <- X.full

savingsSim.func <- function(df.in, modError.in) {
  
######  TRUE MODEL  ######
  
    ##  Create vector of model errors
    mod.epsilon   <- rnorm(nrow(df.in), 0, modError.in)
    
    ##  Create BL and SEM response vectors with true parameters
    kWh.bl       <- (as.matrix(df.in[,1:7]) %*% simCoeff$Value[1:7]) + mod.epsilon
    kWh.meas     <- (as.matrix(df.in[,1:13]) %*% simCoeff$Value)     + mod.epsilon
    true.sav     <- sum(kWh.bl[which(df.in$prog_ind == 1)]) - 
                        sum(kWh.meas[which(df.in$prog_ind == 1)]) + 
                        simCoeff$Value[8]*sum(df.in$event2_post)
    true.pct     <- true.sav / sum(kWh.bl[which(df.in$prog_ind == 1)])
    consump.post <- sum(kWh.bl[which(df.in$prog_ind == 1)])
    
    ##  Subset design matrix for case specific changes
    X <- df.in
    
    ## Separate data into pre/post periods
    X.pre     <- X[which(X$prog_ind == 0),] # two years of pre data
    X.post    <- X[which(X$prog_ind == 1),] # two years of post data
    
    
  ######  FORECAST MODEL  ######
    
    ##  Estimate model
    FC.mod       <- lm(kWh.meas[which(X$prog_ind == 0)] ~ 
                         prod1 
                       + prod2 
                       + noProd_ind 
                       + event1_pre 
                       + HDD55 
                       + prod1_x_HDD55
                       + CDD65
                       , data=X[which(X$prog_ind == 0),])
    FC.mod.sum   <- summary(FC.mod)
    FC.mod.b.hat <- FC.mod$coeff
    
    ##  Model selection criteria
    FC.mod.RMSE    <- sqrt(mean(residuals(FC.mod)^2))
    FC.mod.relRMSE <- FC.mod.RMSE/mean(kWh.meas[which(X$prog_ind == 0)])
    FC.mod.adjR2   <- FC.mod.sum$adj.r.squared
    FC.mod.AIC     <- AIC(FC.mod)
    FC.mod.BIC     <- BIC(FC.mod)
    
    ##  Compute savings estimate
    FC.mod.preds  <- as.matrix(X.post[,c(1:7,14)]) %*% as.matrix(FC.mod.b.hat)
    FC.mod.sav    <- sum(FC.mod.preds - kWh.meas[which(X$prog_ind == 1)]) + 
                          simCoeff$Value[8]*sum(X$event2_post)
    FC.mod.bias   <- true.sav - FC.mod.sav
    FC.mod.pctErr <- FC.mod.bias/true.sav
      
    ##  Estimate SE(savings)
    ##  Identity matrix
    ident.mat   <- diag(nrow(X.post))
    XX.inv      <- solve(t(data.matrix(X.pre[,c(1:7,14)])) %*% data.matrix(X.pre[,c(1:7,14)]))
    post.mats   <- data.matrix(X.post[,c(1:7,14)]) %*% XX.inv %*% t(data.matrix(X.post[,c(1:7,14)]))
    FC.mod.seSav <- sqrt((FC.mod.RMSE^2)*sum(post.mats+ident.mat))
    FC.sav.CV    <- abs(FC.mod.seSav/FC.mod.sav)
    
    FC.mod.incl <- 1
    if( true.sav < (FC.mod.sav-(z.stat*FC.mod.seSav)) |
        true.sav > (FC.mod.sav+(z.stat*FC.mod.seSav)) ) FC.mod.incl <- 0
    
    
######  SIMPLE PRE-POST MODEL  ######
    
    ##  Estimate model
    SPP.mod     <- lm(kWh.meas ~ 
                        prod1 
                      + prod2 
                      + noProd_ind 
                      + event1_pre 
                      + HDD55 
                      + prod1_x_HDD55  
                      + event2_post 
                      + prog_ind
                      + CDD65
                      , data=X)
    SPP.mod.sum <- summary(SPP.mod)
    

    ##  Estimate savings & SE(savings)
    SPP.mod.sav    <- -1*SPP.mod$coeff[which(names(SPP.mod$coefficients) == "prog_ind")]*sum(X.post$prog_ind)
    SPP.mod.bias   <- true.sav - SPP.mod.sav
    SPP.mod.pctErr <- SPP.mod.bias/true.sav
    SPP.mod.seSav  <- SPP.mod.sum$coeff[which(names(SPP.mod$coefficients) == "prog_ind"),2]*sum(X.post$prog_ind)
    SPP.sav.CV     <- abs(SPP.mod.seSav/SPP.mod.sav)
    
    ##  Model selection criteria
    SPP.mod.RMSE    <- sqrt(mean(residuals(SPP.mod)^2))
    SPP.mod.relRMSE <- SPP.mod.RMSE/mean(kWh.meas)
    SPP.mod.adjR2   <- SPP.mod.sum$adj.r.squared
    SPP.mod.AIC     <- AIC(SPP.mod)
    SPP.mod.BIC     <- BIC(SPP.mod)
    
    SPP.mod.incl <- 1
    
    if( true.sav < (SPP.mod.sav-(z.stat*SPP.mod.seSav)) |
        true.sav > (SPP.mod.sav+(z.stat*SPP.mod.seSav)) ) SPP.mod.incl <- 0
    
    
    
######  FULLY SPECIFIED PRE-POST MODEL  ######
    
    ##  Estimate model
    FPP.mod     <- lm(kWh.meas ~ 
                        prod1 
                      + prod2 
                      + noProd_ind 
                      + event1_pre 
                      + HDD55 
                      + prod1_x_HDD55
                      + event2_post 
                      + prog_ind 
                      + prog_x_prod1 
                      + prog_x_prod2
                      + prog_x_noProd 
                      + prog_x_prod1_x_HDD55
                      + CDD65
                      + prog_x_CDD65
                      , data=X)
    FPP.mod.sum <- summary(FPP.mod)
    
    ##  Estimate savings & SE(savings)
    FPP.mod.sav   <- -1*(FPP.mod$coeff[9]*sum(X.post$prog_ind) + 
                         FPP.mod$coeff[10]*sum(X.post$prog_x_prod1) +
                         FPP.mod$coeff[11]*sum(X.post$prog_x_prod2) +
                         FPP.mod$coeff[12]*sum(X.post$prog_x_noProd) +
                         FPP.mod$coeff[13]*sum(X.post$prog_x_prod1_x_HDD55) +
                         FPP.mod$coeff[15]*sum(X.post$prog_x_CDD65)
                         )
    FPP.mod.bias   <- true.sav - FPP.mod.sav
    FPP.mod.pctErr <- FPP.mod.bias/true.sav
    FPP.cov.out    <- vcov(FPP.mod)
    FPP.mod.seSav  <- sqrt(FPP.cov.out[9,9]*sum(X.post$prog_ind)^2 +
                           FPP.cov.out[10,10]*sum(X.post$prog_x_prod1)^2 +
                           FPP.cov.out[11,11]*sum(X.post$prog_x_prod2)^2 +
                           FPP.cov.out[12,12]*sum(X.post$prog_x_noProd)^2 +
                           FPP.cov.out[13,13]*sum(X.post$prog_x_prod1_x_HDD55)^2 +
                           FPP.cov.out[15,15]*sum(X.post$prog_x_CDD65)^2 +
                             
                           2*FPP.cov.out[9,10]*sum(X.post$prog_ind)*sum(X.post$prog_x_prod1) +
                           2*FPP.cov.out[9,11]*sum(X.post$prog_ind)*sum(X.post$prog_x_prod2) +
                           2*FPP.cov.out[9,12]*sum(X.post$prog_ind)*sum(X.post$prog_x_noProd) +
                           2*FPP.cov.out[9,13]*sum(X.post$prog_ind)*sum(X.post$prog_x_prod1_x_HDD55) +
                           2*FPP.cov.out[9,15]*sum(X.post$prog_ind)*sum(X.post$prog_x_CDD65) +
                             
                           2*FPP.cov.out[10,11]*sum(X.post$prog_x_prod1)*sum(X.post$prog_x_prod2) +
                           2*FPP.cov.out[10,12]*sum(X.post$prog_x_prod1)*sum(X.post$prog_x_noProd) +
                           2*FPP.cov.out[10,13]*sum(X.post$prog_x_prod1)*sum(X.post$prog_x_prod1_x_HDD55) +
                           2*FPP.cov.out[10,15]*sum(X.post$prog_x_prod1)*sum(X.post$prog_x_CDD65) +
                             
                           2*FPP.cov.out[11,12]*sum(X.post$prog_x_prod2)*sum(X.post$prog_x_noProd) +
                           2*FPP.cov.out[11,13]*sum(X.post$prog_x_prod2)*sum(X.post$prog_x_prod1_x_HDD55) +
                           2*FPP.cov.out[11,15]*sum(X.post$prog_x_prod2)*sum(X.post$prog_x_CDD65) +
                             
                           2*FPP.cov.out[12,13]*sum(X.post$prog_x_noProd)*sum(X.post$prog_x_prod1_x_HDD55) +
                           2*FPP.cov.out[12,15]*sum(X.post$prog_x_noProd)*sum(X.post$prog_x_CDD65) +

                           2*FPP.cov.out[13,15]*sum(X.post$prog_x_prod1_x_HDD55)*sum(X.post$prog_x_CDD65)
                           )
    FPP.sav.CV    <- abs(FPP.mod.seSav/FPP.mod.sav)
    
    ##  Model selection criteria
    FPP.mod.RMSE    <- sqrt(mean(residuals(FPP.mod)^2))
    FPP.mod.relRMSE <- FPP.mod.RMSE/mean(kWh.meas)
    FPP.mod.adjR2   <- FPP.mod.sum$adj.r.squared
    FPP.mod.AIC     <- AIC(FPP.mod)
    FPP.mod.BIC     <- BIC(FPP.mod)
    
    FPP.mod.incl <- 1
    if( true.sav < (FPP.mod.sav-(z.stat*FPP.mod.seSav)) |
        true.sav > (FPP.mod.sav+(z.stat*FPP.mod.seSav)) ) FPP.mod.incl <- 0
    
    
    
    ######  Create output vector
    
    sim.outVect <- data.frame(consump.post, true.sav, 
                              FC.mod.sav, SPP.mod.sav, FPP.mod.sav,
                              FC.mod.seSav, SPP.mod.seSav, FPP.mod.seSav,
                              FC.sav.CV, SPP.sav.CV, FPP.sav.CV,
                              FC.mod.RMSE, SPP.mod.RMSE, FPP.mod.RMSE,
                              FC.mod.relRMSE, SPP.mod.relRMSE, FPP.mod.relRMSE,
                              FC.mod.bias, SPP.mod.bias, FPP.mod.bias,
                              FC.mod.pctErr, SPP.mod.pctErr, FPP.mod.pctErr,
                              FC.mod.adjR2, SPP.mod.adjR2, FPP.mod.adjR2,
                              FC.mod.AIC, SPP.mod.AIC, FPP.mod.AIC,
                              FC.mod.BIC, SPP.mod.BIC, FPP.mod.BIC,
                              FC.mod.incl, SPP.mod.incl, FPP.mod.incl,
                              row.names=ii)
    
    return(sim.outVect)
    
}


##  Repeat the simulation 10,000 times for each model error input
N.sim <- 10000
modError <- 0.02*mean(simData$sim_kWh[which(simData$prog_ind == 0)])  ## 2% of average daily kWh

for(ii in 1:N.sim) {
  if(ii==1) overspec.sim <- savingsSim.func(X.full, modError)
  else      overspec.sim <- rbind(overspec.sim, savingsSim.func(X.full, modError))
}   

##  Summarize results of sim
simSummary <- data.frame(t(colMeans(overspec.sim)))

##  Compute bias of mean savings estimate
# simSummary$FC.mod.savBias  <- simSummary$FC.mod.sav  - simSummary$true.sav
# simSummary$SPP.mod.savBias <- simSummary$SPP.mod.sav - simSummary$true.sav
# simSummary$FPP.mod.savBias <- simSummary$FPP.mod.sav - simSummary$true.sav

simSummary 

write.xlsx(simSummary, file.path(projPath,"Output","simData_complex - Case 5.xlsx"), 
           sheetName="Sim Out", col.names=T, row.names=F, append=F)
write.xlsx(overspec.sim, file.path(projPath,"Output","simData_complex - Case 5.xlsx"), 
           sheetName="Sim Output - prod-hdd-event", col.names=T, row.names=F, append=T)









