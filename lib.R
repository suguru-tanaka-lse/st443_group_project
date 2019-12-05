library('matlib')
library('corpcor')
library('MASS')
library('glmnet')
library('dplyr')
library('plyr')
library('ROCR')
library('ROCit')
library('ggplot2')
library('caret')
library(MESS)
library(e1071)
library(glasso)
# This library contains functions used for data generation and node-wide lasso.

# make_theta(): a function to create a theta matrix.
#   input: p, the dimension of theta matrix.
#   output: theta matrix

# generate_rvdf(): a function that generates r.v. according to sigma matrix.
#  Sigma matrix is inverse of the theta matrix.
#  input: 1) n: number of samples to generate, 
#         2) sigma matrix
#  output: dataframe containing the n samples.

# piece_lasso(): a function to regress one variable to 
#   the rest of variables in input df.
# input: 1) rvdf: dataframe containing both x and y.
#        2) i: the index of variable to be the response variable, 1<=i<=p.
#        3) l: the tuning parameter, lambda.
# output: df of predicted coefficients.

# full_lasso(): a function to iterate piece_lasso() over all i in [1, p]
# input: 1) rvdf: dataframe containing both x and y.
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
generate_rvdf = function(n, sigma){ #generate random samples according to the 
  p = ncol(sigma)
  rv=mvrnorm(n=n,mu=rep(0,p),Sigma=sigma)
  index=paste("X",1:p,sep="")
  colnames(rv)=paste("X",1:p,sep="")
  rownames(rv)=c(1:n)
  rvdf=data.frame(rv)
  return(rvdf)
}
piece_lasso=function(rvdf, i,l){ #l means lambda
  y = rvdf %>% select(i)
  x = rvdf %>% select(-i)
  x = as.matrix(x)
  y = as.matrix(y)
  coef = coef(glmnet(x,y,lambda=l, intercept=FALSE))
  coef = as.matrix(coef)
  return(coef)
}
full_lasso=function(rvdf, l){
  p = ncol(rvdf)
  pred = as.data.frame(t(piece_lasso(rvdf, 1, l)))
  for (i in seq(2,p)){
    t_coef = as.data.frame(t(piece_lasso(rvdf, i, l)))
    t_coef[1,1]=i
    pred = suppressMessages(join(t_coef, pred, type = "full"))
  }
  pred = pred[, -1]
  pred = pred[seq(dim(pred)[1],1),]
  rownames(pred) = paste("Beta(",1:p,",",sep="")
  colnames(pred) = paste(",",1:p,")",sep="") 
  pred = pred!=0 #turn the matrix to T&F to make count convient
  return(pred)
}
find_nw = function(lasso_out, method){
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
      } else if (method == 'lasso2'){
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

find_glasso <- function(rvdf, l){
  sigma_hat <- cov(rvdf)
  glasso_out <- glasso(sigma_hat, l)
  theta_hat <- glasso_out$wi
  theta_hat_bool= (theta_hat !=0)
  theta_hat_bool[lower.tri(theta_hat_bool, diag = TRUE)] <- NA
  return(theta_hat_bool)
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

get_lambda_conf_matrix = function(E, rvdf, lambda, method){
  TP_result <- c()
  FP_result <- c()
  TN_result <- c()
  FN_result <- c()
  for(l in lambda){
    if (method %in% c('lasso1', 'lasso2')){
      lasso_out <- full_lasso(rvdf, l)
      pred_E <- find_nw(lasso_out, method)
    } else if ( method == 'glasso' ){
      pred_E <- find_glasso(rvdf, l)
    } else {
      stop("undefined method: should be either 'lasso1', 'lasso2' or 'glasso'")
    }
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
  result_matrix <- t(rbind(lambda, TP_result, FP_result, TN_result, FN_result))
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- c("lambda", "TP", "FP", "TN", "FN")
  return(result_df)
}

get_ROC = function(lambda_conf_matrix){
  lcm <- lambda_conf_matrix
  TPR_result <- lcm$TP / (lcm$TP + lcm$FN)
  FPR_result <- lcm$FP / (lcm$FP + lcm$FN)
  result_matrix <- t(rbind(lcm$lambda, TPR_result, FPR_result))
  result_df <- as.data.frame(result_matrix)
  #rownames(result_df) <- lambda
  colnames(result_df) <- c("lambda", "TPR","FPR")
  return(result_df)
}

get_recovery <- function(lambda_conf_matrix){
  lcm <- lambda_conf_matrix
  FPR_result <- lcm$FP / (lcm$FP + lcm$FN)
  FNR_result <- lcm$FN / (lcm$FN + lcm$TP)
  result_matrix <- t(rbind(FPR_result, FNR_result))
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- c("FPR", "FNR")
  return(result_df)
}

get_accuracy = function(lambda_conf_matrix){
  lcm <- lambda_conf_matrix
  accuracy <- (lcm$TP + lcm$TN) / (lcm$TP + lcm$FN + lcm$TN +lcm$FP)
  result_matrix <- t(rbind(lcm$lambda, accuracy))
  result_df <- as.data.frame(result_matrix)
  colnames(result_df) <- c("lambda","accuracy")
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
# optimal_lambda <- (result_df %>% top_n(1, f1_score))[['lambda']]
# highest_f1 <- (result_df %>% top_n(1, f1_score))[['f1_score']]
# return(list(result_df = result_df, optimal_lambda = optimal_lambda, highest_f1 = highest_f1))
