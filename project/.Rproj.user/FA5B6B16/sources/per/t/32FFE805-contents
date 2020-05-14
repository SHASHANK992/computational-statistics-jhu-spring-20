X <- f12
h_0 <- silverman_h
A <- function(x, h_0 = silverman_h) {
  # (sqrt(pi) / (2 * h_0)) * ((h_0^2) - (x^2))
  (1 / (2 * pi)) * (sqrt(pi) / (h_0^3)) * ((8 * x^4) + (2 * (h_0^2) * (x^3)) -
                            (2 * (h_0^2) * x) + ((h_0^4) / 4) + ((h_0^2) / 2))
}

A_part <- sum(A(X))
B <- function(u, v, h_0 = silverman_h) {
  # D <- h_0 * sqrt(pi) * (1 / (2 * pi))
  D <- h_0 * sqrt(pi)
  H <- exp(-((u - v) / (2 * h_0))^2)

  UV <- -3 * (((u + v) / (2 * h_0))^4) + (((u + v) / (2 * h_0))^2) +
    (3 / 4) + (1 / (2 * (h_0^2))) * ((2 * u) + (2 * v) + (4 * u * v) - (u^2) - (v^2)) +
    (1 / (h_0^4)) *
      (((u + v)^3) + (u * v * ((u + v)^2)) + (2 * u * v * (u + v)) + (u*v)^2)

  return(D * H * UV)
}
B <- function(u, v, h_0 = silverman_h) {
  D <- h_0 * sqrt(pi) * (1 / (2 * pi))
  H <- exp(-((u - v) / (2 * h_0))^2)

  UV <- (1/h_0^4) * ( ((u+v) / 2)^4 + 3 * ((u+v) / 2)^2 * (h_0^2/2) + h_0^4/2 + ((u+v) / 2)^2 * h_0^2 + ((u+v) / 2)^2 * (h_0^2/2) + h_0^4/4) -
    ((2/h_0^4) * (u+v)) * ( ((u+v) / 2)^3 + 2 * ((u+v) / 2) * (h_0^2/2) + ((u+v) / 2) * (h_0^2/2) ) +
    ((1/h_0^4) * (u^2 + 4 * u * v + v^2) - (2 / h_0^2)) * (((u+v) / 2)^2 + (h_0^2 / 2)) +
    ((2/h_0^4) * (u^2 * v + u * v^2) + (2 / h_0^2) * (u + v)) * ((u+v)/2) -
    (1 / h_0^2) * (u^2 + v^2) + (1/h_0^4) * (u^2 * v^2) + 1

  return(D * H * UV)
}
B_vec <- Vectorize(B)

B_part <- outer(X = x, Y = x, FUN = B_vec)
diag(B_part) <- 0
B_part2 <- sum(B_part)
R_F <- (B_part2 + A_part) / ((n^2) * (h_0^6))
R_F <- sum(outer(X = x, Y = x, FUN = B_vec)) / ((n^2) * (h_0^6))
xs <- head(x,2)
sum(B_vec(u = xs, v = xs)) / ((2^2) * (h_0^6))



### hist est
hist_est <- function(var1 = f12, h = silverman_h, m = 30) {
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
    m_bins[lower_interval] <- m_bins[lower_interval] + 1
  } # calculate the number of obs. in each bin
  
  p_k <- m_bins / n # proportion per bin
  
  # y_limits <- c(-0.05, max(f_hat_pts) + 0.5 * sd(f_hat_pts))
  y_limits <- c(0, max(p_k) + 1e-2)
  # Reference: https://stackoverflow.com/questions/4785657/how-to-draw-an-empty-plot
  plot(1, type="n", xlab="", ylab="", xlim=var1_range, ylim=y_limits,
       main = 'Histogram Estimator')
  abline(h = 0)
  for (i in 1:length(p_k)) {
    lines(c(break_pts[i], break_pts[i+1]), c(p_k[i], p_k[i]))
  } # plot horizontal and vertical lines of histogram
  segments(x0 = break_pts[-1], y0 = 0, x1 = break_pts[-1], y1 = p_k)
  segments(x0 = break_pts[-length(break_pts)],
           y0 = 0, x1 = break_pts[-length(break_pts)], y1 = p_k)
}

hist_est(var1 = f12, h = silverman_h, m = 30)
