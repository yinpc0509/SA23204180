#' @import DAAG
#' @import boot
#' @import bootstrap
#' @import coda
#' @import dplyr
#' @import ggplot2
#' @import kableExtra
#' @import microbenchmark
NULL

#' @importFrom microbenchmark microbenchmark
#' @title An Illustration Dataset
#' @name illustration_data
#' @description A dataset used to illustrate the performance of \code{vaccR} and \code{vaccC}.
#' @examples
#' \dontrun{
#' data(illustration_data)
#' attach(illustration_data)
#' tm <- microbenchmark::microbenchmark(
#'   vR = predictResponse(age, female, ily),
#'   vC = predictResponseCpp(age, female, ily)
#' )
#' print(summary(tm)[, c(1, 3, 5, 6)])
#' }
NULL

#' @title Benchmark R and Rcpp Functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions (\code{gibbsR} and \code{predictResponse}) and Cpp functions (\code{gibbsC} and \code{predictResponseCpp}).
#' @examples
#' \dontrun{
#' data(illustration_data)
#' attach(illustration_data)
#' tm1 <- microbenchmark::microbenchmark(
#'   rnR = gibbsR(100, 10),
#'   rnC = gibbsC(100, 10)
#' )
#' print(summary(tm1)[, c(1, 3, 5, 6)])
#' 
#' tm2 <- microbenchmark::microbenchmark(
#'   vR = predictResponse(age, female, ily),
#'   vC = predictResponseCpp(age, female, ily)
#' )
#' print(summary(tm2)[, c(1, 3, 5, 6)])
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma
#' @useDynLib StatComp
NULL




#' @title Predict Response Using Three Inputs in R
#' @description This function predicts a response using three input variables.
#'   The prediction model is described in http://www.babelgraph.org/wp/?p=358.
#' @param age the first predictor (numeric)
#' @param is_female the second predictor (logical)
#' @param has_ily the third predictor (logical)
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' prediction <- predictResponse(age, is_female, has_ily)
#' }
#' @export
predictResponse <- function(age, is_female, has_ily) {
  probability <- 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * has_ily
  probability <- probability * ifelse(is_female, 1.25, 0.75)
  probability <- pmax(0, probability)
  probability <- pmin(1, probability)
  probability
}


#' @title A Custom Gibbs Sampler using R
#' @description A custom implementation of Gibbs sampler using R
#' @param num_samples the number of samples to generate
#' @param thinning_factor the number of between-sample random numbers
#' @return a random sample matrix of size \code{num_samples}
#' @examples
#' \dontrun{
#' custom_samples <- customGibbsR(100, 10)
#' par(mfrow=c(2,1))
#' plot(custom_samples[, 1], type='l', main='Sampled Values for X')
#' plot(custom_samples[, 2], type='l', main='Sampled Values for Y')
#' }
#' @export
customGibbsR <- function(num_samples, thinning_factor) {
  sample_matrix <- matrix(nrow = num_samples, ncol = 2)
  x <- y <- 0
  
  for (i in 1:num_samples) {
    for (j in 1:thinning_factor) {
      x <- rgamma(1, 3, y * y + 4)
      y <- rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)))
    }
    sample_matrix[i, ] <- c(x, y)
  }
  
  sample_matrix
}


#' Convex Relaxation Function
#'
#' This function implements the convex relaxation algorithm for policy learning.
#' It iteratively updates the parameter matrix theta based on the features and rewards matrices.
#'
#' @param features_mat A matrix representing the features of the dataset.
#' @param rewards_mat A matrix representing the rewards associated with different actions.
#' @param learning_rate The learning rate for parameter updates (default is 0.1).
#' @param max_iter The maximum number of iterations for updating parameters (default is 10000).
#' @return A matrix representing the learned parameters theta.
#' @examples
#' \dontrun{
#' set.seed(123)
#' features_mat <- matrix(rnorm(250), ncol = 5)
#' beta <- matrix(rnorm(10), ncol = 2)
#' rewards_mat <- features_mat %*% beta + matrix(runif(100), ncol = 2)
#' theta_result_convex <- convex_relaxation(features_mat, rewards_mat)
#' print("Theta Result (Convex Relaxation):")
#' print(theta_result_convex)
#' }
#'
convex_relaxation <- function(features_mat, rewards_mat, learning_rate = 0.1, max_iter = 10000) {
  # n为样本数量，p为样本特征维度，k为动作数目
  n <- nrow(features_mat)
  p <- ncol(features_mat)
  k <- ncol(rewards_mat)
  # 初始化参数矩阵
  theta <- matrix(0, nrow = p, ncol = k)
  # 初始化动作选择矩阵
  act_star_array <- integer(n)
  # 根据奖励矩阵拿出每一时刻的最优动作
  for (i in 1:n) {
    act_star_array[i] <- which.max(rewards_mat[i,])
  }
  # 循环训练参数
  for (i in 1:max_iter) {
    # 随机指标
    idx <- sample(1:n, 1)
    # 得到具体动作数值
    act_star <- as.integer(act_star_array[idx])
    # 估计动作为特征映射后与奖励差值最小的
    act_hat <- which.min(rewards_mat[idx,] - features_mat[idx,] %*% theta)
    # 参数更新
    theta[, act_hat] <- theta[, act_hat] - (learning_rate) * features_mat[idx,]
    theta[, act_star] <- theta[, act_star] + (learning_rate) * features_mat[idx,]
  }
  
  return(theta)
}

#' Logistic Policy Function
#'
#' This function implements a simple logistic policy learning algorithm.
#' It iteratively updates the parameter matrix theta based on the features and rewards matrices.
#'
#' @param features_mat A matrix representing the features of the dataset.
#' @param rewards_mat A matrix representing the rewards associated with different actions.
#' @param max_iter The maximum number of iterations for updating parameters (default is 1000).
#' @param learning_rate The learning rate for parameter updates (default is 1).
#' @param reg_param The regularization parameter (default is 0).
#' @param theta_init Initial values for the parameter matrix theta (default is NULL).
#' @param intercept Include intercept term in the logistic regression (default is TRUE).
#' @return A matrix representing the learned parameters theta.
#' @examples
#' \dontrun{
#' set.seed(123)
#' features_mat <- matrix(rnorm(250), ncol = 5)
#' beta <- matrix(rnorm(10), ncol = 2)
#' rewards_mat <- features_mat %*% beta + matrix(runif(100), ncol = 2)
#' theta_result_logistic <- logistic_policy(features_mat, rewards_mat)
#' print("Theta Result (Logistic Policy):")
#' print(theta_result_logistic)
#' }
#'
logistic_policy <- function(features_mat, rewards_mat, max_iter = 1000, 
                            learning_rate = 1, reg_param = 0, theta_init = NULL, intercept = TRUE) {
  
  # 获取数据集的行数、特征数和动作数
  n <- nrow(features_mat)
  p <- ncol(features_mat)
  k <- ncol(rewards_mat)
  
  # 初始化参数矩阵 theta
  if (is.null(theta_init)) {
    theta <- matrix(0, nrow = p, ncol = k)
  } else {
    theta <- theta_init
  }
  
  # 迭代更新参数
  for (i in 1:max_iter) {
    
    # 计算分数矩阵
    scores_mat <- features_mat %*% theta
    
    # 计算指数分数矩阵
    exp_scores_mat <- exp(scores_mat)
    
    # 计算每行的指数分数占总和的比例
    exp_scores_mat_div_sum_exp <- exp_scores_mat / rowSums(exp_scores_mat, na.rm = TRUE)
    
    # 计算临时矩阵
    temp_mat <- (exp_scores_mat_div_sum_exp - exp_scores_mat_div_sum_exp^2) * rewards_mat
    
    # 初始化临时数组 temp
    temp <- array(0, dim = c(n, p, k))
    
    # 利用循环计算 temp 数组
    for (i in 1:n) {
      for (j in 1:p) {
        for (m in 1:k) {
          temp[i, j, k] = array(temp_mat, dim = c(n, 1, k))[i, 1, m] * array(features_mat, dim = c(n, p, 1))[i, j, 1]
        }
      }
    }
    
    # 计算梯度矩阵 grad_theta_mat
    grad_theta_mat <- apply(temp, c(2, 3), mean)
    
    # 更新参数 theta
    theta <- theta + (learning_rate) * (grad_theta_mat)
    
    # 对参数 theta 进行操作
    theta = t(t(theta) - theta[, 1])
  }
  
  # 返回最终学到的参数矩阵 theta
  return(theta)
}

