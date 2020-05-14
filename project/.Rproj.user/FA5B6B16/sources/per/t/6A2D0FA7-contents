# Fig. 4
### hist_est plots a histogram after separating the
### data into m intervals.
hist_est <- function(var1 = f12, m = 30) {
  var1_range <- range(var1) # range of data
  D <- var1_range[2] - var1_range[1] # ~length D of histogram
  v_k <- D / m # volume of bin (length of interval)
  n <- length(var1)
  
  break_pts <- matrix(NA, nrow = m + 1)
  break_pts[1] <- var1_range[1]; break_pts[m + 1] <- var1_range[2]
  for (i in 1:(m - 1)) { 
    break_pts[i + 1] <- var1_range[1] + i * v_k # possible bias towards right side
  } # create (m + 1) intervals
  break_pts <- as.vector(break_pts)
  
  m_bins <- matrix(0, nrow = m)
  for (i in var1) {
    lower_interval <- tail(which(i > break_pts), 1)
      if (i == break_pts[1]) { # edge case for smallest value
        m_bins[1] <- m_bins[1] + 1
      }
    m_bins[lower_interval] <- m_bins[lower_interval] + 1
  } # calculate the number of obs. in each bin
  
  p_k <- m_bins / n # proportion per bin
  f_hat <- p_k / v_k # est. histogram output
  
  # y_limits <- c(-0.05, max(f_hat_pts) + 0.5 * sd(f_hat_pts))
  y_limits <- c(0, max(f_hat) + 1e-2)
  # Reference: https://stackoverflow.com/questions/4785657/how-to-draw-an-empty-plot
  plot(1, type="n", xlab="", ylab="", xlim=var1_range, ylim=y_limits,
       main = paste0('Histogram Estimator (m=', m,')'))
  abline(h = 0)
  for (i in 1:length(f_hat)) {
    lines(c(break_pts[i], break_pts[i+1]), c(f_hat[i], f_hat[i]))
  } # plot horizontal and vertical lines of histogram
  segments(x0 = break_pts[-1], y0 = 0, x1 = break_pts[-1], y1 = f_hat)
  segments(x0 = break_pts[-length(break_pts)],
           y0 = 0, x1 = break_pts[-length(break_pts)], y1 = f_hat)
}

par(mfrow = c(3,2))
m_vec <- seq(from = 10, to = 60, length.out = 6)
for (i in m_vec) {
  hist_est(var1 = f12, m = i)
  lines(density(f12))
}
sapply(m_vec, function(x) hist_est(var1 = f12, m = x))

hist_est(var1 = f12, m = 20)

# Fig. 5
### hist_est2 provides the estimated probability based off
### the histogram estimate of the data. This is different
### from hist_est which plots the histogram.
hist_est2 <- function(x, var1 = f12, m = 30) {
  var1_range <- range(var1) # range of data
  
  # Check if x outside range of f-hat
  if ((x > var1_range[2]) | (x < var1_range[1])) {
    return(0)
  }
  
  D <- var1_range[2] - var1_range[1] # ~length D of histogram
  v_k <- D / m # volume of bin (length of interval)
  n <- length(var1)
  
  break_pts <- matrix(NA, nrow = m + 1)
  break_pts[1] <- var1_range[1]; break_pts[m + 1] <- var1_range[2]
  for (i in 1:(m - 1)) { 
    break_pts[i + 1] <- var1_range[1] + i * v_k # possible bias towards right side
  } # create (m + 1) intervals
  break_pts <- as.vector(break_pts)
  
  m_bins <- matrix(0, nrow = m)
  for (i in var1) {
    lower_interval <- tail(which(i > break_pts), 1)
    m_bins[lower_interval] <- m_bins[lower_interval] + 1
  } # calculate the number of obs. in each bin
  
  p_k <- m_bins / n # proportion per bin
  f_hat <- p_k / v_k # est. histogram output
  kth_bin <- tail(which(x >= break_pts), 1) # x's bin

  if ((kth_bin >= 1) & (kth_bin <= m)) {
    return(f_hat[kth_bin]) # return f-hat
  } else if (kth_bin == (m + 1)) {
    return(f_hat[m]) # edge case for last bin, return m'th bin
  } else {
    return(0)
  }
}

hist_est2_vec <- Vectorize(hist_est2, vectorize.args = c("x"))
xs <- seq(-3, 1.5, length.out = 1e4)
plot(xs, hist_est2_vec(x = xs, var1 = f12, m = 40),
     type = 'l', main = 'Histogram Estimator and K.D.E.',
     xlab = latex2exp::TeX('$x$'), ylab = latex2exp::TeX('$\\hat{f}(x)$'))
lines(density(f12), lty = 2, col = 'red')
legend("topleft", legend = c('Histogram Estimator', 'K.D.E.'),
       lty = c(1,2), col = c('black', 'red'))

# Fig. 6
plot(density(f12))
set.seed(0)
f12_sample <- sample(f12, 5)
points(f12_sample, rep(0, length(f12_sample)), pch = 3)
# hist_est2_vec(x = f12_sample)
# draw a normal kernel for each of the points
h <- density(f12)$bw
f12_sample[1]

# normal pdf
normal <- function(z) { dnorm(z) }

sum(normal((f12_sample[1] - f12) / h)) / (n * h)

lines(x = c(1,2), y = c(1,2))

xs <- c(1, 1.5, 1.7, 2.5)
plot(density(xs))
points(xs, rep(0, length(xs)), pch = 3)
zs <- seq((xs[1] - 1), (xs[1] + 1), length.out = 5)
lines(zs, dnorm(zs))

# Fig. 7
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

histogram <- function(df = f12, round_digits = 3,
                      m = 10, bandwidth_type = 'Terrell',
                      kernel_function = 'Normal') {
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
       main = paste0('Histogram of X (m=', m,')'),
       xlab = 'x', sub = paste0('Bandwidth type: ', bandwidth_type,
                                '; Kernel function: ', kernel_function))
}


kde <- function(var1 = f12, h = silverman_h, K = normal, histogram_round = 3,
                histogram_m = 10, bw_name = 'Silverman', K_name = 'Normal') {
  # n <- length(var1)
  D <- c( # histogram range
    -round(abs(range(var1)[1]), digits = histogram_round), # lower range
    round(range(var1)[2], digits = histogram_round) # upper range
  )
  x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
  f_hat_pts <- f_hat_vec(x = x_pts, X_i = var1, h = h, K = K) # KDE
  histogram(df = var1, round_digits = histogram_round,
            m = histogram_m, bandwidth_type = bw_name, kernel_function = K_name) # plot
  lines(x_pts, f_hat_vec(x = x_pts, X_i = var1, h = h, K = K)) # add KDE
}

# left-column: uniform, epanechnikov, triweight (m=40)
# right-column: normal (m=20, 40, 60) Silverman, Terrell, SJ

# Silverman's rule of thumb
n <- length(f12)
sigma_hat <- sd(f12)
silverman_h <- ((4 / (3 * n))^(1 / 5)) * sigma_hat

g <- function(y, X = f12, h_0 = silverman_h, A = 1e3) {
  n <- length(X)
  outer_matrix <- outer(X = X, Y = X, FUN = function(u, v) {
    (exp(-0.5 * ((y - u) / h_0)^2) *
       (1 - ((y - u) / h_0))^2) *
      (exp(-0.5 * ((y - v) / h_0)^2) *
         (1 - ((y - v) / h_0))^2)
  })
  L_val <- sum(outer_matrix)
  
  g_output <- (2 * A) * (1 / (n * (h_0^3)))^2 * L_val
  return(g_output)
}
g_vec <- Vectorize(g)

test_values <- seq(-10,10,length.out = 21)
test_values[which(g_vec(test_values, A = 1e3) > 1e-3)]
set.seed(0)
u_sample <- runif(1e4, -1e3, 1e3)
g_eval <- g_vec(y = u_sample, A = 1e3)
R_f <- mean(g_eval)

R_K <- 1 / (2 * sqrt(pi))
h_hat <- (R_K / (n * R_f))^(1 / 5)
# sheather_jones_h <- h_hat
sheather_jones_h <- bw.SJ(f12)
sheather_jones_h <- locfit::sjpi(f12, silverman_h)[1]

# Terrell's maximal smoothing principle
R_K <- 1 / (2 * sqrt(pi))
terrell_h <- 3 * ((R_K / (35 * n))^(1 / 5)) * sigma_hat

par(mfrow = c(1,2))
hist_est(var1 = f12, m = 40)
D <- range(f12)
x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
f_hat_pts_uniform <- f_hat_vec(x = x_pts, X_i = f12, h = sheather_jones_h, K = uniform)
lines(x_pts, f_hat_pts_uniform, lty = 1, col = 'blue') # add KDE
f_hat_pts_epanechnikov <- f_hat_vec(x = x_pts, X_i = f12, h = sheather_jones_h, K = epanechnikov)
lines(x_pts, f_hat_pts_epanechnikov, lty = 2, col = 'red') # add KDE
f_hat_pts_triweight <- f_hat_vec(x = x_pts, X_i = f12, h = sheather_jones_h, K = triweight)
lines(x_pts, f_hat_pts_triweight, lty = 3, col = 'green') # add KDE
legend("topleft", legend = c('Uniform', 'Epanechnikov', 'Triweight'),
       lty = c(1,2,3), col = c('blue', 'red', 'green'))

hist_est(var1 = f12, m = 40)
f_hat_pts_sj <- f_hat_vec(x = x_pts, X_i = f12, h = sheather_jones_h, K = normal)
lines(x_pts, f_hat_pts_sj, lty = 1, col = 'blue') # add KDE
f_hat_pts_1.5 <- f_hat_vec(x = x_pts, X_i = f12, h = 1.5, K = normal)
lines(x_pts, f_hat_pts_1.5, lty = 2, col = 'red') # add KDE
f_hat_pts_0.05 <- f_hat_vec(x = x_pts, X_i = f12, h = 0.05, K = normal)
lines(x_pts, f_hat_pts_0.05, lty = 3, col = 'green') # add KDE
legend("topleft", legend = c('Sheather-Jones', '1.5', '0.05'),
       lty = c(1,2,3), col = c('blue', 'red', 'green'))

# Fig. 8
normal
uniform
epanechnikov
triangle <- function(z) {
  ifelse(abs(z) < 1, 1 - abs(z), 0)
}

biweight <- function(z) {
  ifelse(abs(z) < 1, (15 / 16) * ((1 - z^2)^2), 0)
}

xs <- seq(-3, 3, length.out = 1e5)
plot(xs, uniform(xs), type = 'l', ylim = c(0,1.1),
     main = 'Kernel Functions', col = 'black',
     xlab = latex2exp::TeX('z'), ylab = latex2exp::TeX('K(z)'))
lines(xs, normal(xs), lty = 2, col = 'red')
lines(xs, epanechnikov(xs), lty = 3, col = 'green')
lines(xs, triangle(xs), lty = 4, col = 'blue')
lines(xs, biweight(xs), lty = 5, col = 'black')
lines(xs, triweight(xs), lty = 6, col = 'blue')
legend("topleft", legend = c('Uniform', 'Normal', 'Epanechnikov',
                             'Triangle', 'Biweght', 'Triweight'),
       lty = c(1:6), col = c('black', 'red', 'green', 'blue', 'black', 'blue'))

# Fig. 9
scale_bandwidth <- function(h_K, delta_K, delta_L) {
  h_L <- (h_K / delta_K) * delta_L
  return(h_L)
}

delta_normal <- (1 / (2 * sqrt(pi)))^(1 / 5)
delta_uniform <- (9 / 2)^(1 / 5)
delta_epanechnikov <- 15^(1 / 5)
delta_triangle <- 24^(1 / 5)
delta_biweight <- 35^(1/ 5)
delta_triweight <- (9450 / 143)^(1 / 5)

delta_vec <- c(delta_uniform, delta_epanechnikov,
               delta_triangle, delta_biweight, delta_triweight)
scale_vec <- sapply(delta_vec, function(x)
  scale_bandwidth(h_K = sheather_jones_h,
    delta_K = delta_normal, delta_L = x))

par(mfrow = c(1,2))
hist_est(var1 = f12, m = 40)
D <- range(f12)
x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
f_hat_pts_uniform <- f_hat_vec(x = x_pts, X_i = f12, h = sheather_jones_h, K = uniform)
lines(x_pts, f_hat_pts_uniform, lty = 1, col = 'blue') # add KDE
f_hat_pts_epanechnikov <- f_hat_vec(x = x_pts, X_i = f12, h = sheather_jones_h, K = epanechnikov)
lines(x_pts, f_hat_pts_epanechnikov, lty = 2, col = 'red') # add KDE
f_hat_pts_triweight <- f_hat_vec(x = x_pts, X_i = f12, h = sheather_jones_h, K = triweight)
lines(x_pts, f_hat_pts_triweight, lty = 3, col = 'green') # add KDE
legend("topleft", legend = c('Uniform', 'Epanechnikov', 'Triweight'),
       lty = c(1,2,3), col = c('blue', 'red', 'green'))
###
hist_est(var1 = f12, m = 40)
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = f12, h = sheather_jones_h, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 1, col = 'black', lwd = 2) # add KDE
f_hat_pts_uniform <- f_hat_vec(x = x_pts, X_i = f12, h = scale_vec[1], K = uniform)
lines(x_pts, f_hat_pts_uniform, lty = 1, col = 'blue') # add KDE
f_hat_pts_epanechnikov <- f_hat_vec(x = x_pts, X_i = f12, h = scale_vec[2], K = epanechnikov)
lines(x_pts, f_hat_pts_epanechnikov, lty = 2, col = 'red') # add KDE
f_hat_pts_triweight <- f_hat_vec(x = x_pts, X_i = f12, h = scale_vec[5], K = triweight)
lines(x_pts, f_hat_pts_triweight, lty = 3, col = 'green') # add KDE
legend("topleft", legend = c('Normal','Uniform', 'Epanechnikov', 'Triweight'),
       lty = c(1,1,2,3), lwd = c(2,1,1,1), col = c('black', 'blue', 'red', 'green'))
