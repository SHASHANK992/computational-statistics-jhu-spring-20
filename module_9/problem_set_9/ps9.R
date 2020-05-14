### Problem 1
# part (a)
stomach <- c(25, 42, 45, 46, 51, # Load data
  103, 124, 146, 340, 396, 412, 876, 1112)
breast <- c(24, 40, 719, 727, 791,
  1166, 1235, 1581, 1804, 3460, 3808)
log_stomach <- log(stomach)
log_breast <- log(breast)

par(mfrow=c(1,2))
boxplot(stomach)
boxplot(breast)
dev.off()

par(mfrow=c(1,2))
boxplot(log(stomach))
boxplot(log(breast))
dev.off()

mean(log_stomach) # 4.96792
mean(log_breast) # 6.558603

### Finds the bootstrap t C.I. on log scale
bootstrap_t <- function(df, B = 1e3, alpha = 0.05) {
  log_df <- log(df) # Take log scale
  n <- length(df) # Initialize variables
  T_F <- mean(log_df)
  data_matrix <- matrix(NA, nrow = B, ncol = length(log_df))
  theta_matrix <- matrix(NA, nrow = B)
  R_X_F_matrix <- matrix(NA, nrow = B)
  set.seed(999) # Set seed
  for (i in 1:B) { # Bootstrap sample
    data_matrix[i,] <- sample(log_df, size = length(log_df), replace = TRUE)
    theta_matrix[i] <- mean(data_matrix[i,])
    V_F_star <- (1 / (n - 1)) * sum((data_matrix[i,] - theta_matrix[i])^2)
    R_X_F_matrix[i] <- (theta_matrix[i] - T_F) / sqrt(V_F_star)
  }
  
  theta_bar <- mean(theta_matrix) # mean bootstrap theta
  # variance of bootstrap theta
  var_theta_hat_star <- sum((theta_matrix - theta_bar)^2) / (B - 1)
  
  # Find quantiles
  xi <- quantile(R_X_F_matrix, probs = c(alpha / 2, 1 - alpha / 2))
  names(xi) <- NULL
  bootstrap_t <- c(T_F - sqrt(var_theta_hat_star) * xi[2],
                   T_F - sqrt(var_theta_hat_star) * xi[1])
  return(bootstrap_t)
}
stomach_bootstrap_t_CI <- bootstrap_t(df = stomach, B = 1e4, alpha = 0.05) # 4.772891 5.195032
breast_bootstrap_t_CI <- bootstrap_t(df = breast, B = 1e4, alpha = 0.05) # 5.891237 6.803640

# part (b)
# H_0: no diff. in mean survival times (mu_log_breast = mu_log_stomach)
# H_1: there is a diff. in mean survival times (mu_log_breast =/= mu_log_stomach)
n1 <- length(log_stomach)
n2 <- length(log_breast)

t0 <- mean(log_stomach) - mean(log_breast)
m <- choose((n1 + n2), n2)
m_combinations <- t(combn(c(log_stomach, log_breast), n2)) # m combinations of both data
stomach_breast <- c(log_stomach, log_breast) # Combined data
t_matrix <- matrix(NA, nrow = m) # Empty matrix for t_i

for (i in 1:m) { # Calculate t_i for all m combinations of the data
  ys <- m_combinations[i,]
  # Reference: https://stackoverflow.com/questions/5812478/how-i-can-select-rows-from-a-dataframe-that-do-not-match
  xs <- subset(stomach_breast, !(stomach_breast %in% ys))
  xbar_new <- mean(xs); ybar_new <- mean(ys)
  t_matrix[i] <- xbar_new - ybar_new
}
boxplot(t_matrix)
k <- which(t0 == sort(t_matrix, decreasing = TRUE))
k / m # 0.9925185
# rej. H_0 with sig. level. 99.25%

# part (c)
# percentile method
percentile_method <- function(df = breast, B = 1e3, alpha = 0.05, take_log = TRUE) {
  if (take_log == TRUE) {
    log_df <- log(df) # Take log scale
  } else {
    log_df <- df
  }
  data_matrix <- matrix(NA, nrow = B, ncol = length(log_df))
  theta_matrix <- matrix(NA, nrow = B)
  set.seed(888) # Set seed
  for (i in 1:B) { # Bootstrap sample
    data_matrix[i,] <- sample(log_df, size = length(log_df), replace = TRUE)
    theta_matrix[i] <- mean(data_matrix[i,])
  }
  
  percentile_CI <- quantile(theta_matrix, probs = c(alpha / 2, 1 - alpha / 2))
  names(percentile_CI) <- NULL
  
  return(percentile_CI)
}
breast_percentile_CI <- percentile_method(df = breast) # 5.617261 7.375876
exp(breast_percentile_CI) # 275.1347 1596.9896
exp(breast_bootstrap_t_CI) # 361.8524 901.1214
mean(breast) # 1395.909
mean(log_breast) # 6.558603
breast_percentile_CI_nonlog <- percentile_method(df = breast, # 778.1273 2145.3841
                                                 take_log = FALSE)

stomach_percentile_CI <- percentile_method(df = stomach) # 4.352572 5.631810
exp(stomach_percentile_CI) # 77.67798 279.16689
exp(stomach_bootstrap_t_CI) # 118.2606 180.3738
mean(stomach) # 286
mean(log_stomach) # 4.96792
stomach_percentile_CI_nonlog <- percentile_method(df = stomach, # 133.30 492.25
                                                 take_log = FALSE)
### Problem 2
# Cauchy distribution
xs <- seq(-4, 4, length.out = 1e4)
plot(xs, dnorm(x = xs, mean = 0, sd = 1), type = 'l',
     main = 'Normal and Cauchy', xlab = 'x', ylab = 'f(x)')
lines(xs, dcauchy(x = xs, location = 0, scale = 2), type = 'l', col = 'red')
legend("topleft", legend = c('Normal', 'Cauchy'), col = c('black', 'red'), lty = c(1,1))

par(mfrow=c(2,2))
cauchy_bootstrap <- function(sample_size = 10, B = 1e4, alpha = 0, beta = 2) {
  set.seed(888)
  cauchy_sample <- rcauchy(n = sample_size, location = alpha, scale = beta) # Original sample
  data_matrix <- matrix(NA, nrow = B, ncol = length(cauchy_sample)) # Initialize variables
  theta_matrix <- matrix(NA, nrow = B)
  set.seed(888) # Set seed
  for (i in 1:B) { # Bootstrap sample
    data_matrix[i,] <- sample(cauchy_sample, size = length(cauchy_sample), replace = TRUE)
    theta_matrix[i] <- mean(data_matrix[i,])
  }
  hist(theta_matrix, main = paste('Histogram of the Sample Mean (n=',sample_size,')'),
       xlab = 'Bootstrap Sample Mean',
       sub = paste('Original Sample Mean:', round(mean(cauchy_sample), 4),
                   ' Bootstrap Sample Mean:', round(mean(theta_matrix), 4)))
}
cauchy_bootstrap(sample_size = 1e1)
cauchy_bootstrap(sample_size = 5e1)
cauchy_bootstrap(sample_size = 1e2)
cauchy_bootstrap(sample_size = 1e3)

# Uniform distribution
xs <- seq(0, 10, length.out = 1e4)
plot(xs, dunif(xs, 0, 1e4), type = 'l', ylim = c(0, 0.00014))
segments(x0 = c(0,10), y0 = 0, x1 = c(0,10), y1 = dunif(xs, 0, 1e4))
abline(h = 0)

theta <- 10
set.seed(999)
unif_sample <- runif(n = 1e1, min = 0, max = theta)
sample(unif_sample, size = length(unif_sample), replace = TRUE)
theta_hat <- max(unif_sample)

unif_bootstrap <- function(sample_size = 100, B = 1e4, theta_max = 10) {
  set.seed(999) # Set seed
  unif_sample <- runif(n = sample_size, min = 0, max = theta_max) # Original sample
  data_matrix <- matrix(NA, nrow = B, ncol = sample_size) # Initialize variables
  theta_matrix <- matrix(NA, nrow = B)
  set.seed(999) # Set seed
  for (i in 1:B) { # Bootstrap sample
    data_matrix[i,] <- sample(unif_sample, size = sample_size, replace = TRUE)
    theta_matrix[i] <- max(data_matrix[i,])
  }
  hist(theta_matrix, main = paste('Histogram of the Maximum Order Statistic (n=',sample_size,')'),
       xlab = 'Bootstrap Maximum Order Statistic',
       sub = paste('Original Max. Order Statistic:', round(max(unif_sample), 4),
                   ' Mean Bootstrap Max. Order Statistic:', round(mean(theta_matrix), 4)))
}
par(mfrow = c(2,2))
unif_bootstrap(sample_size = 1e1, B = 1e4)
unif_bootstrap(sample_size = 1e2, B = 1e4)
unif_bootstrap(sample_size = 1e3, B = 1e4)
unif_bootstrap(sample_size = 1e4, B = 1e4)
dev.off()

### Problem 4
# part (a)
set.seed(999)
sample_size <- 1e2
random_sample <- rnorm(n = sample_size, mean = 0, sd = 1)
sample_mean <- mean(random_sample)

# part (b)
B <- 1e1
set.seed(999) # Set seed
data_matrix <- matrix(NA, nrow = B, ncol = sample_size) # Initialize variables
theta_matrix <- matrix(NA, nrow = B)
set.seed(999) # Set seed
for (i in 1:B) { # Bootstrap sample
  data_matrix[i,] <- sample(random_sample, size = sample_size, replace = TRUE)
  theta_matrix[i] <- mean(data_matrix[i,])
}
hist(theta_matrix, main = paste('Histogram of the Sample Mean (n=',sample_size,')'),
     xlab = 'Bootstrap Sample Mean',
     sub = paste('Original Sample Mean:', round(sample_mean, 4),
                 ' Mean Bootstrap Sample Mean:', round(mean(theta_matrix), 4)))
mu_bar_star <- mean(theta_matrix)
mu_bar_star - sample_mean
2 * sample_mean - mu_bar_star

(1 / (B - 1)) * sum((theta_matrix - mu_bar_star)^2)

# part (c)
concat_B <- rep(random_sample, B) # create n*B vector
set.seed(999)
balanced_bootstrap <- sample(concat_B, # permute n*B vector
                             size = length(concat_B), replace = FALSE)
# Reference: https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
balanced_groups <- split(balanced_bootstrap, # split into B groups
                         rep_len(1:B, length(balanced_bootstrap)))
balanced_theta <- unlist(lapply(balanced_groups, mean)) # calculate theta for each
names(balanced_theta) <- NULL

hist(balanced_theta, main = paste('Histogram of the Balanced Sample Mean (n=',sample_size,')'),
     xlab = 'Balanced Bootstrap Sample Mean',
     sub = paste('Original Sample Mean:', round(sample_mean, 4),
                 ' Mean Bootstrap Sample Mean:', round(mean(balanced_theta), 4)))
mu_bar_star <- mean(balanced_theta)
mu_bar_star - sample_mean
2 * sample_mean - mu_bar_star
(1 / (B - 1)) * sum((theta_matrix - mu_bar_star)^2)
