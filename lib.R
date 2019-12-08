library(matlib)
library(corpcor)
library(MASS)
library(glmnet)
library(dplyr)
library(plyr)
library(ROCR)
library(ROCit)
library(ggplot2)
library(caret)
library(MESS)
library(e1071)
library(glasso)
library(CVglasso)
library(ggpubr)
# This library contains functions used for data generation and node-wide lasso.

# make_theta(): a function to create a theta matrix.
#   input: p, the dimension of theta matrix.
#   output: theta matrix

# generate_rv(): a function that generates r.v. according to sigma matrix.
#  Sigma matrix is inverse of the theta matrix.
#  input: 1) n: number of samples to generate, 
#         2) sigma matrix
#  output: dataframe containing the n samples.

# piece_lasso(): a function to regress one variable to 
#   the rest of variables in input df.
# input: 1) rv: dataframe containing both x and y.
#        2) i: the index of variable to be the response variable, 1<=i<=p.
#        3) l: the tuning parameter, lambda.
# output: df of predicted coefficients.

# full_lasso(): a function to iterate piece_lasso() over all i in [1, p]
# input: 1) rv: dataframe containing both x and y.
#        2) l: the turning parameter, lambda.
# output: pred: true and false df of coefficient matrix.
#        diagonal part is NA.

# find_nw(): a function to obtain node-wise lasso estimator,
#  based on coefficient matrix obtained by full_lasso().
# input: 1) lasso_out: True and False df of coefficient matrix.
#         This expects the output of full_lasso()
#        2) method: method to use:
#           'lasso1' or 'lasso2' is available.
# output: upper triangle matrix of coefficient matrix.
#         True if the two edge is predicted to be connected.
# ex: find_nw(lasso_out, 'lasso1')

# get_confusion_matrix(): a function that returns confusion matrix of 
#  estimated E and predicted E. E is an upper triangle matrix of coefficient,
#  whose entry is True if the edge is connected and False otherwise.
# input: 1) pred_matrix: an upper triangle matrix of estimated coefficients.
#        2) real_matrix: an upper triangle matrix of real coefficients.

make_theta = function(p, q){
  I = diag(p)
  theta = I
  for (row in 1:nrow(I)){
    for (col in 2:ncol(I)){
      if (row < col){
        theta[row, col] <- sample(x=c(0.5,0),prob=c(q, 1-q),size=1,replace=T)
      }
    }
  }
  theta[lower.tri(theta)]=t(theta)[lower.tri(theta)]
  repeat{
    if(is.positive.definite(theta)==T){
      break
    } else {
      theta=theta+I
    }
  }
  theta = cov2cor(theta)
  return(theta)
}
generate_rv = function(n, sigma){ #generate random samples according to the 
  p = ncol(sigma)
  rv=mvrnorm(n=n,mu=rep(0,p),Sigma=sigma)
  index=paste("X",1:p,sep="")
  colnames(rv)=paste("X",1:p,sep="")
  rownames(rv)=c(1:n)
  rv=data.frame(rv)
  return(rv)
}
piece_lasso=function(rv, i, lambda){
  y = rv %>% select(i)
  x = rv %>% select(-i)
  x = as.matrix(x)
  y = as.matrix(y)
  if (length(lambda) > 1) {
    lambda = cv.glmnet(x,y, alpha = 1, family = 'gaussian', lambda = lambda, intercept = FALSE)$lambda.1se
  }
  coef = coef(glmnet(x,y,alpha = 1, family = 'gaussian', lambda = lambda, intercept=FALSE))
  coef = as.matrix(coef)
  return(coef)
}
node_wise_lasso=function(rv, lambda){
  p = ncol(rv)
  pred = as.data.frame(t(piece_lasso(rv, 1, lambda)))
  for (i in seq(2,p)){
    t_coef = as.data.frame(t(piece_lasso(rv, i, lambda)))
    t_coef[1,1]=i
    pred = suppressMessages(join(t_coef, pred, type = "full"))
  }
  pred = pred[, -1]
  pred = pred[seq(dim(pred)[1],1),]
  rownames(pred) = paste("Beta(",1:p,",",sep="")
  colnames(pred) = paste(",",1:p,")",sep="") 
  return(pred)
}
graphical_lasso <- function(rv, lambda){
  if (length(lambda) > 1) {
    theta_hat <- CVglasso(rv, lam = lambda, trace = "none")$Omega
  } else {
    sigma_hat <- cov(rv)
    glasso_out <- glasso(sigma_hat, lambda)
    theta_hat <- glasso_out$wi
  }
  return(theta_hat)
}
get_true_edge <- function(theta){
  E <- (theta!=0) 
  E[lower.tri(E, diag = TRUE)] <- NA
  return(E)  
}
get_estimated_edge = function(lasso_out, method){
  lasso_out <- (lasso_out!=0) #turn the matrix to T&F to make count convient
  p = ncol(lasso_out)
  row_fill = 1
  col_fill = 2
  E_hat <- matrix(NA, nrow = p, ncol = p)
  result = lasso_out
  repeat{
    if(row_fill==p)
      break
    repeat{
      if(col_fill==p+1)
        break
      E_hat[row_fill, col_fill] <- FALSE
      if (method == 'lasso1'){
        if(result[row_fill,col_fill]== TRUE && result[col_fill,row_fill]==TRUE)
          E_hat[row_fill, col_fill] <- TRUE
      } else if (method == 'lasso2' || method == 'glasso'){
        # glasso can be in the previous part since lasso_out is symmetric.
        if(result[row_fill,col_fill] == TRUE | result[col_fill,row_fill]==TRUE)
          E_hat[row_fill, col_fill] <- TRUE
      } else {
        stop("unknown method")
      }
      col_fill=col_fill+1
    }
    row_fill=row_fill+1
    col_fill=row_fill+1
  }
  return(E_hat)
}

get_confusion_matrix = function(pred_matrix, real_matrix){
  pred_vector <- factor(pred_matrix, levels = c('TRUE', 'FALSE'))
  pred_vector <- pred_vector[!is.na(pred_vector)]
  real_vector <- factor(real_matrix, levels = c('TRUE', 'FALSE'))
  real_vector <- real_vector[!is.na(real_vector)]
  out <- confusionMatrix(pred_vector, real_vector, positive = 'TRUE')
  result <- out[['table']]
  return(result)
}

make_lambda_conf_matrix = function(E, rv, lambda_series, method){
  TP_result <- c()
  FP_result <- c()
  TN_result <- c()
  FN_result <- c()
  for(lambda in lambda_series){
    if (method %in% c('lasso1', 'lasso2')){
      lasso_out <- node_wise_lasso(rv, lambda)
    } else if ( method == 'glasso' ){
      lasso_out <- graphical_lasso(rv, lambda)
    } else {
      stop("undefined method: should be either 'lasso1', 'lasso2' or 'glasso'")
    }
    pred_E <- get_estimated_edge(lasso_out, method)
    conf_E <- get_confusion_matrix(pred_E, E)
    TP = conf_E[1,1]
    FP = conf_E[1,2]
    FN = conf_E[2,1]
    TN = conf_E[2,2]
    TP_result = c(TP_result, TP)
    FP_result = c(FP_result, FP)
    TN_result = c(TN_result, TN)
    FN_result = c(FN_result, FN)
  }
  result_matrix <- t(rbind(lambda_series, TP_result, FP_result, TN_result, FN_result))
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- c("lambda", "TP", "FP", "TN", "FN")
  result_df['method'] <- method
  return(result_df)
}

get_cv_conf_matrix = function(E, rv, lambda_series, method){
  TP_result <- c()
  FP_result <- c()
  TN_result <- c()
  FN_result <- c()
  if (method %in% c('lasso1', 'lasso2')){
    lasso_out <- node_wise_lasso(rv, lambda_series)
  } else if ( method == 'glasso' ){
    lasso_out <- graphical_lasso(rv, lambda_series)
  } else {
    stop("undefined method: should be either 'lasso1', 'lasso2' or 'glasso'")
  }
  pred_E <- get_estimated_edge(lasso_out, method)
  conf_E <- get_confusion_matrix(pred_E, E)
  TP = conf_E[1,1]
  FP = conf_E[1,2]
  FN = conf_E[2,1]
  TN = conf_E[2,2]
  TP_result = c(TP_result, TP)
  FP_result = c(FP_result, FP)
  TN_result = c(TN_result, TN)
  FN_result = c(FN_result, FN)
  result_matrix <- t(rbind(TP_result, FP_result, TN_result, FN_result))
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- c("TP", "FP", "TN", "FN")
  result_df['method'] <- method
  return(result_df)
}

get_tpr_fpr = function(lambda_conf_matrix){
  lcm <- lambda_conf_matrix
  TPR_result <- lcm$TP / (lcm$TP + lcm$FN)
  FPR_result <- lcm$FP / (lcm$FP + lcm$TN)
  result_matrix <- t(rbind(lcm$lambda, TPR_result, FPR_result))
  result_df <- as.data.frame(result_matrix)
  #rownames(result_df) <- lambda
  colnames(result_df) <- c("lambda", "TPR","FPR")
  result_df['method'] <- method
  return(result_df)
}

get_recovery <- function(conf_matrix){
  cm <- conf_matrix
  FPR_result <- cm$FP / (cm$FP + cm$TN)
  FNR_result <- cm$FN / (cm$FN + cm$TP)
  result_matrix <- t(rbind(FPR_result, FNR_result))
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- c("FPR", "FNR")
  return(result_df)
}

get_accuracy = function(conf_matrix){
  cm <- conf_matrix
  accuracy <- (cm$TP + cm$TN) / (cm$TP + cm$FN + cm$TN +cm$FP)
  result_matrix <- accuracy
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- c("accuracy")
  return(result_df)
}

get_f1 = function(lambda_conf_matrix){
  lcm <- lambda_conf_matrix
  precision <- lcm$TP / (lcm$TP + lcm$FP)
  recall <- lcm$TP / (lcm$TP + lcm$FN)
  f1 <- 2 * (precision * recall) / (precision + recall) 
  result_matrix <- t(rbind(lcm$lambda, f1))
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- c("lambda","f1_score")
  return(result_df)
}

