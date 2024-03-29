## -----------------------------------------------------------------------------
# convex_relaxation 函数
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

## -----------------------------------------------------------------------------
# 定义 logistic_policy 函数，该函数实现了一个简单的逻辑回归策略学习算法
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


## -----------------------------------------------------------------------------
# Example Usage for convex_relaxation:
set.seed(123)
features_mat <- matrix(rnorm(250), ncol = 5)
beta <- matrix(rnorm(10), ncol = 2)
rewards_mat <- features_mat %*% beta+matrix(runif(100), ncol = 2)

theta_result_convex <- convex_relaxation(features_mat, rewards_mat)
print("Theta Result (Convex Relaxation):")
print(theta_result_convex)

theta_result_logistic <- logistic_policy(features_mat, rewards_mat)
print("Theta Result (theta_result_logistic):")
print(theta_result_logistic)

