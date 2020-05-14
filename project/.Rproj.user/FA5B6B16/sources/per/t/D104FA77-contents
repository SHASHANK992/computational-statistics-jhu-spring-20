### 10.1
f12 <- read.csv('F12.txt', header = FALSE) # load data
f12 <- log(f12[[1]]) # take log

# custom histogram function
histogram <- function(df = f12, round_digits = 3, m = 10, bandwidth_type = 'Terrell') {
  D <- c( # histogram range
    -round(abs(range(df)[1]), digits = round_digits), # lower range
    round(range(df)[2], digits = round_digits) # upper range
  )
  # Reference: https://stackoverflow.com/questions/22051345/breaking-a-mathematical-range-in-equal-parts
  bin_width <- (D[2] - D[1]) / m # width of each bin

  # create (m + 1) break points for the m intervals
  break_pts <- matrix(NA, nrow = m + 1)
  break_pts[1] <- D[1]; break_pts[m + 1] <- D[2]
  for (i in 1:(m - 1)) { 
    break_pts[i + 1] <- D[1] + i * bin_width
  }

  hist(x = df, breaks = break_pts, freq = FALSE, # plot histogram
       main = paste('Histogram of log X', '(m=', m,')'),
       xlab = 'log F12', sub = bandwidth_type)
}

# normal pdf
normal <- function(z) { dnorm(z) }

# f-hat function for KDE + vectorized
f_hat <- function(x, X_i, h, K) {
  n <- length(X_i)
  (1 / h) * mean(K((x - X_i) / h))
}
f_hat_vec <- Vectorize(f_hat, vectorize.args = 'x')

# plot KDE
kde <- function(var1 = f12, h = silverman_h, K = normal,
                histogram_round = 3, histogram_m = 10, bw_name = 'Silverman') {
  # n <- length(var1)
  D <- c( # histogram range
    -round(abs(range(var1)[1]), digits = histogram_round), # lower range
    round(range(var1)[2], digits = histogram_round) # upper range
  )
  x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
  f_hat_pts <- f_hat_vec(x = x_pts, X_i = var1, h = h, K = K) # KDE
  histogram(df = var1, round_digits = histogram_round,
            m = histogram_m, bandwidth_type = bw_name) # plot
  lines(x_pts, f_hat_vec(x = x_pts, X_i = var1, h = h, K = K)) # add KDE
}

# part (a)
# Silverman's rule of thumb
n <- length(f12)
sigma_hat <- sd(f12)
silverman_h <- ((4 / (3 * n))^(1 / 5)) * sigma_hat

m_vec <- seq(10, 60, length.out = 6)
par(mfrow = c(3,2))
sapply(m_vec, function(x)
  kde(var1 = f12, h = silverman_h, K = normal, histogram_round = 3,
      histogram_m = x, bw_name = 'Silverman')
)

# Sheather-Jones approach
# Monte Carlo Integration
# g function
g <- function(y, X = f12, h_0 = silverman_h, A = 1e5) {
  n <- length(X)
  # A <- sum(
  #   sapply(X, function(z) {
  #     ( exp(-0.5 * ((y - z) / h_0)^2) *
  #       (1 - ((y - z) / h_0))^2 )^2
  #   })
  # )

  outer_matrix <- outer(X = X, Y = X, FUN = function(u, v) {
    (exp(-0.5 * ((y - u) / h_0)^2) *
       (1 - ((y - u) / h_0))^2) *
      (exp(-0.5 * ((y - v) / h_0)^2) *
         (1 - ((y - v) / h_0))^2)
  })
  # print(paste('A == B?', A == sum(diag(B_part))))
  L_val <- sum(outer_matrix)
  
  
  g_output <- 2 * A * (1 / (n * (h_0^3)))^2 * L_val
  return(g_output)
}
g_vec <- Vectorize(g)

test_values <- seq(-10,10,length.out = 21)
test_values[which(g_vec(test_values) > 1e-3)]
set.seed(0)
u_sample <- runif(1e5, -2,  2)
g_eval <- g_vec(y = u_sample)
R_f <- mean(g_eval)

R_K <- 1 / (2 * sqrt(pi))
h_hat <- (R_K / (n * R_f))^(1 / 5)
sheather_jones_h <- h_hat

par(mfrow = c(3,2))
sapply(m_vec, function(x)
  kde(var1 = f12, h = sheather_jones_h, K = normal, histogram_round = 3,
      histogram_m = x, bw_name = 'Sheather-Jones')
)

locfit::sjpi(x = f12, a = silverman_h)

# Terrell's maximal smoothing principal
R_K <- 1 / (2 * sqrt(pi))
terrell_h <- 3 * ((R_K / (35 * n))^(1 / 5)) * sigma_hat

par(mfrow = c(3,2))
sapply(m_vec, function(x)
  kde(var1 = f12, h = terrell_h, K = normal, histogram_round = 3,
      histogram_m = x, bw_name = 'Terrell')
)

# Combined histogram + kde estimates
histogram_combined <- function(df = f12, round_digits = 3, m = 40) {
  D <- c( # histogram range
    -round(abs(range(df)[1]), digits = round_digits), # lower range
    round(range(df)[2], digits = round_digits) # upper range
  )
  # Reference: https://stackoverflow.com/questions/22051345/breaking-a-mathematical-range-in-equal-parts
  bin_width <- (D[2] - D[1]) / m # width of each bin
  
  # create (m + 1) break points for the m intervals
  break_pts <- matrix(NA, nrow = m + 1)
  break_pts[1] <- D[1]; break_pts[m + 1] <- D[2]
  for (i in 1:(m - 1)) { 
    break_pts[i + 1] <- D[1] + i * bin_width
  }
  
  hist(x = df, breaks = break_pts, freq = FALSE, # plot histogram
       main = paste('Histogram of X', '(m=', m,')'),
       xlab = 'log F12')
}

kde_combined <- function(var1 = f12, h_s = silverman_h,
                         h_sj = sheather_jones_h, h_t = terrell_h,
                         K = normal, histogram_round = 3, histogram_m = 40) {
  D <- c( # histogram range
    -round(abs(range(var1)[1]), digits = histogram_round), # lower range
    round(range(var1)[2], digits = histogram_round) # upper range
  )
  x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
  histogram_combined(df = var1, round_digits = histogram_round, m = histogram_m) # plot
  
  # Silverman
  lines(x_pts, f_hat_vec(x = x_pts, X_i = var1, h = h_s, K = K), col = 'red', lty = 2)
  # Sheather-Jones
  lines(x_pts, f_hat_vec(x = x_pts, X_i = var1, h = h_sj, K = K), col = 'black', lty = 1)
  # Terrell
  lines(x_pts, f_hat_vec(x = x_pts, X_i = var1, h = h_t, K = K), col = 'blue', lty = 3)
  legend("topleft", legend = c('Silverman', 'Sheather-Jones', 'Terrell'),
         col = c('red', 'black', 'blue'), lty = c(2,1,3))
}

# part (b)
# Uniform
uniform <- function(z) {
  ifelse(abs(z) < 1, 1 / 2, 0)
}

# Epanechnikov
epanechnikov <- function(z) {
  ifelse(abs(z) < 1, (3 / 4) * (1 - z^2), 0)
}

# Triweight
triweight <- function(z) {
  ifelse(abs(z) < 1, (35 / 32) * ((1 - z^2)^3), 0)
}

par(mfrow = c(3,2))
sapply(m_vec, function(x)
  kde(var1 = f12, h = sheather_jones_h, K = uniform, histogram_round = 3,
      histogram_m = x, bw_name = 'Uniform')
)
sapply(m_vec, function(x)
  kde(var1 = f12, h = sheather_jones_h, K = epanechnikov, histogram_round = 3,
      histogram_m = x, bw_name = 'Epanechnikov')
)
sapply(m_vec, function(x)
  kde(var1 = f12, h = sheather_jones_h, K = triweight, histogram_round = 3,
      histogram_m = x, bw_name = 'Triweight')
)

# histogram estimator
hist_est <- function(x, var1 = f12, h = silverman_h) {
  var1_range <- range(var1) # range of data

  # Check if x outside range of f-hat
  if ((x > var1_range[2]) | (x < var1_range[1])) {
    return(0)
  }

  D <- var1_range[2] - var1_range[1] # length D of support
  # v_k <- D / m # volume of bin (length of interval)
  v_k <- h # volume of bin (length of interval)
  n <- length(var1)
  m <- ceiling(D / v_k)
  
  break_pts <- matrix(NA, nrow = m + 1) # Find break points
  break_pts[1] <- var1_range[1]; break_pts[m + 1] <- var1_range[2]
  for (i in 1:(m - 1)) { 
    break_pts[i + 1] <- var1_range[1] + i * v_k # possible bias towards right side
  } # create (m + 1) intervals
  break_pts <- as.vector(break_pts)
  
  m_bins <- matrix(0, nrow = m)
  for (i in var1) { # calculate the number of obs. in each bin
    lower_interval <- tail(which(i > break_pts), 1)
    m_bins[lower_interval] <- m_bins[lower_interval] + 1
    if (i == min(var1)) { # when i == min value
      m_bins[1] <- m_bins[1] + 1
    }
  }
  
  p_k <- m_bins / n # proportion per bin
  kth_bin <- tail(which(x >= break_pts), 1)
  f_hat <- p_k / v_k

  if ((kth_bin >= 1) & (kth_bin <= m)) {
    return(f_hat[kth_bin]) # return f-hat
  } else if (kth_bin == (m + 1)) {
    return(f_hat[m]) # edge case for last bin, return m'th bin
  } else {
    return(0)
  }
}
hist_est_vec <- Vectorize(hist_est, vectorize.args = c('x'))
xs <- seq(-3, 1.5, length.out = 1e3)
plot(xs, hist_est_vec(x = xs, var1 = f12, h = silverman_h),
     type = 'l', main = 'Histogram Estimator',
     xlab = latex2exp::TeX('$x$'), ylab = latex2exp::TeX('$\\hat{f}(x)$'))

