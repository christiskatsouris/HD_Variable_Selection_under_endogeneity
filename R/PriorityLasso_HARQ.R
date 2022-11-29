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
library(prioritylasso)

source( "aggRV.R" )
source( "aggY.R" )

source( "HARDataCreation.R" )
source( "dataset_creation.R" )
source( "RQs_sqrt_demeaned_function.R" )

##############################################################################
### Step 1: Input Dataset
##############################################################################

dataset <- read.csv( "dataset.csv", header = TRUE )
time    <- as.matrix( dataset$Date )

RV  <- read.csv( "RV.csv",  header = TRUE )
RV  <- as.matrix(RV)

RQ  <- read.csv( "RQ.csv",  header = TRUE )
RQ  <- as.matrix(RQ)

###############################################################################
### Function 1: fast_LASSO_HAR_MODEL
###############################################################################

fast_LASSO_HARQ_MODEL <- function(  x1 = X.RV.corrected, x2 = X.RQ.corrected, y = Y1, which.firm =  which.firm  )
{# BEGIN OF FUNCTION
  
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  y  <- as.matrix(y)
  k <- which.firm
  
  # At this step, I re-arrange the design matrix so that it has always the 3 lagged RVs of the particular firm at the first 3 columns
  x1.firm   <- x1[ , ( (3*k-2):(3*k) ) ]
  x2.firm   <- x2[ , ( 3*k - 2 ) ]
  
  x.firm <- cbind( x1.firm, x2.firm ) 
  
  x1.others <- x1[ , - ( (3*k-2):(3*k) ) ]
  x2.others <- x2[ , - ( 3*k - 2 ) ]
  
  x.others <- cbind( x1.others, x2.others )
  x.new    <- cbind( x.firm , x.others  ) 
  
  blocks  <- list( bp1 = 1:4, bp2 = 5:108 )
  fits <- prioritylasso( X = x.new , Y = y, family = "gaussian", type.measure = "mse", blocks = blocks, standardize = FALSE, block1.penalization = FALSE )
  
  res <- as.matrix( as.vector( fits$coefficients ) )
  return( res )
  
}# END OF FUNCTION

###############################################################################
### Function 2: FAST_HARestimate
###############################################################################

FAST_HARQestimate = function( RV = RV, RQ = RQ, periods = periods, h = h, which.firm =  which.firm )
{# BEGIN OF FUNCTION
  
  RV <- RV
  RQ <- RQ
  periods <- periods
  h <- h
  k <- which.firm
  
  RVData   <- HARDataCreation( RM = RV , periods = periods, h = h  )
  RQData   <- HARDataCreation( RM = RQ , periods = periods, h = h  )
  X.RV     <- as.matrix( RVData$Design_matrix )
  X.RQ     <- as.matrix( RQData$Design_matrix )
  
  design.matrix.RQ <- matrix(0, nrow = nrow(X.RQ), ncol = 27 )
  
  #At this step only select the lag 1 RQ for the correction
  for (j in 1:27)
  { 
    design.matrix.RQ[ ,j ] <- as.matrix( X.RQ[  , (3*j - 2) ] )
  }
  
  X.RQ.corrected <- design.matrix.RQ[-nrow(design.matrix.RQ),]
  X.RV.corrected <- X.RV[-nrow(X.RV),]
  
  Y.RV        <- RVData$ymatrix
  Y.corrected <- Y.RV[-nrow(Y.RV),]
  Y1 <- Y.corrected[,k]
  
  HARQ_Lasso  <- fast_LASSO_HARQ_MODEL( x1 = X.RV.corrected, x2 = X.RQ.corrected, y = Y1, which.firm =  k )
  return(list("coefficients" = HARQ_Lasso , "X.RV" = X.RV, "X.RQ" = design.matrix.RQ ) )
  
}# END OF FUNCTION

###############################################################################
### Function 3: OOS_HAR_Forecasts
###############################################################################

OOS_HARQ_Forecasts <- function( RV = RV, RQ = RQ, periods = periods, window.size = window.size, h = 1, which.firm = which.firm )
{# begin of HARForecast
  
  RV1 <- RV
  RQ1 <- RQ
  periods     <- periods
  window.size <- window.size
  k  <- which.firm
  h <- 1
  
  ######Initialization section ######
  T      <- length(RV1[,2])
  Lags   <- length(periods)
  
  ######Initialization end #######
  lFit = list()
  
  ######Forecasting #######
  nRoll <- ( T - window.size )
  mForecast           <- matrix( 0 , nrow = nRoll, ncol = 1 ) 
  beta.priority.lasso <- matrix( 0 , nrow = nRoll, ncol = 109) # 82+27
  
  for (j in 1:nRoll)
  {# begin for-loop for nAhead = 1
   
    model.Fit   <- FAST_HARQestimate( RV = RV1[j:(T-nRoll+j),  ], RQ = RQ1[j:(T-nRoll+j),  ] , periods = periods, h = h, which.firm =  k )
    beta.priority.lasso[j, ] <- model.Fit$coefficients
    
    X.RV <- model.Fit$X.RV
    X.RQ <- model.Fit$X.RQ
    
    
    x1.firm   <- X.RV[ , ( (3*k-2):(3*k) ) ]
    x2.firm   <- X.RQ[ , ( 3*k - 2 ) ]
    
    x.firm <- cbind( x1.firm, x2.firm ) 
    
    x1.others <- X.RV[ , - ( (3*k-2):(3*k) ) ]
    x2.others <- X.RQ[ , - ( 3*k - 2 ) ]
    
    x.others <- cbind( x1.others, x2.others )
    x.new    <- cbind( x.firm , x.others  ) 
    
    vCoef          <- model.Fit$coefficients
    mForecast[j,1] <- vCoef[1] + sum( vCoef[-1]*x.new [ nrow(x.new), ] )
    
  }# end for-loop for nAhead = 1 
  
  # the actual RV of the rolling window
  actual.RV   <- RV1[(T-nRoll+1):T, which.firm ]
  forecast.RV <-  mForecast  
  
  forecast.error.matrix <- cbind( actual.RV , forecast.RV )
  forecast.error <- ( actual.RV - forecast.RV )
  
  forecast.error.matrix.new <- cbind( forecast.error.matrix, forecast.error  )  
  return( list( "forecast_error_HARQ" = forecast.error.matrix.new, "betas_HARQ_priortity_lasso" =  beta.priority.lasso ) )
  
}# end of HARForecast  

###############################################################################

RV  <- dataset_creation( dataset = RV,  time = time )
RQ_demeaned <- RQs_sqrt_demeaned_function( RQ = RQ )
RQ          <- dataset_creation( dataset = RQ_demeaned, time = time )

periods     <- c(1, 5, 22)
window.size <- 4000

for (i in 1:1)
{  
  i <- 1
  forecast   <- OOS_HARQ_Forecasts( RV = RV, RQ = RQ, periods = periods, window.size = window.size, h = 1, which.firm = i )
  forecasts  <- forecast$forecast_error_HARQ
  betas      <- forecast$betas_HARQ_priortity_lasso
  
  name1  <- paste("forecast_HARQ_priority_lasso_w4000_firm_", i, sep="")
  file   <- paste( as.character(name1), ".csv", sep="" )
  write.csv( forecasts, file = paste( as.character(name1),".csv", sep="" ), row.names=FALSE )
  
  name2  <- paste("betas_HARQ_priority_lasso_w4000_firm_", i, sep="")
  file   <- paste( as.character(name2), ".csv", sep="" )
  write.csv( betas, file = paste( as.character(name2),".csv", sep="" ), row.names=FALSE )
}

###############################################################################
