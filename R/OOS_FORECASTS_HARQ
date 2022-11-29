##########################################
# R Script Details:
##########################################

# Script name: priorityLASSO_HARQ.R

# Program aim: This R program estimates out-of-sample betas and forecasts for the HARQ model 
# based on the priority Lasso shrinkage norm. 

# written by: 

# Christis G. Katsouris (December 2019)
# Department of Economics
# University of Southampton
# Southampton, United Kingdom

######################################################################################

############################################################################
### REQUIRED PACKAGES
############################################################################

library(zoo)
library(dplyr)
library(lubridate)

library(xts)
library(highfrequency)
library(Matrix)

library(glmnet)
library("caret")
library("elasticnet")

library("lattice")
library("ggplot2")
library("lars")

source( "aggRV.R" )
source( "aggY.R" )
source( "HARDataCreation.R" )
source( "dataset_creation.R" )
source( "RQs_sqrt_demeaned_function.R" )
source( "HAR_estimate.R" )
source( "HARJ_estimate.R" )
source( "HARQ_estimate.R" )

###############################################################################  
## Function: OOS FORECASTS
###############################################################################

dataset <- read.csv( "dataset.csv", header = TRUE )
time    <- as.matrix( dataset$Date )

RV  <- read.csv( "RV.csv",  header = TRUE )
RV  <- as.matrix(RV)
RV  <- dataset_creation( dataset = RV,  time = time )


RQ  <- read.csv( "RQ.csv",  header = TRUE )
RQ  <- as.matrix(RQ)
RQ  <- dataset_creation( dataset = RQ,  time = time )



OOS_FORECASTS_HARQ <- function( RM1 = RV, RM2 = RQ,  nRoll=nRoll, windowType = "rolling" )
{# BEGIN OF FUNCTION
  
  RV    <- RM1
  RQ    <- RM2
  nRoll <- nRoll
  windowType <- windowType
  
  ######Initialization section ######
  periods <- c(1,5,22)
  h  <- 1
  iT <- length(RV[,2])
  
  lFit = list()
  
  ######Forecasting: One-Period-Ahead #######
  forecast.matrix.windows <- matrix(0 , nrow = nRoll , ncol = 27 ) 
  
  
  for (j in 1:nRoll)
  {# begin for-loop for one-period ahead forecasting 
    
    HARData <- HARDataCreation( RM = RV[j:(iT-nRoll+j),  ], periods = periods, h = h )
    Y       <- HARData$ymatrix
    X1      <- HARData$Design_matrix
    X1.corrected <- X1[-nrow(X1),] 
    Y.corrected  <- Y[-nrow(Y),]
    
    
    RJData   <- HARDataCreation( RM = RJ[j:(iT-nRoll+j),  ], periods = periods, h = h )
    X2       <- RJData$Design_matrix
    X2.corrected <- X2[-nrow(X2),] 
    
    forecasts <- HARJ_estimate(  x1 = X1.corrected, y1 = Y.corrected, x2 = X2.corrected  )
    forecast.matrix.windows[j, ] <- forecasts
    
  }# end for-loop for one-period ahead forecasting
  
  # the actual RV of the rolling window
  actual.RV   <- RV[(iT-nRoll+1):iT,  ]
  forecast.RV <-  forecast.matrix.windows  
  MSE <- ( actual.RV - forecast.RV )^2
  
  return( MSE )
  
}# END OF FUNCTION

###############################################################################

MSE_HARQ <- OOS_FORECASTS_HARQ( RM1 = RV, RM2 = RQ, nRoll=1000, windowType = "rolling" )

mean.MSE.windows <- matrix( 0, nrow = 1000, ncol =1 )

for (i in 1:1000)
{
  mean.MSE.windows[i,1] <- mean( MSE_HARQ[i , ] )
}

mean(  mean.MSE.windows )


write.csv( MSE_CHARJ, "MSE_RW_HARQ.csv" , row.names = FALSE )
write.csv( mean.MSE.windows, "mean_MSE_RW_HARQ.csv" , row.names = FALSE  )

