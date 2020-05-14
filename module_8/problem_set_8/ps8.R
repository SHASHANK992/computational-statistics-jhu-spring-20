### Problem 1
cvdataX <- scan(file.choose()) # Load variables
cvdataY <- scan(file.choose())
n <- length(cvdataX)

# part (a)
# Plot the (X,Y)
plot(cvdataX, cvdataY, main = 'CVdata',
     xlab = 'X', ylab = 'Y') # plot data
mu_hat <- mean(cvdataX) # MLE mu
x_seq <- seq(from = range(cvdataX)[1], to = range(cvdataX)[2],
             length.out = 100) # Generate x range
lines(x_seq, # Plot fitted line
      dnorm(x = x_seq, mean = mu_hat, sd = sqrt(2)),
      col = 'blue')
legend("topleft", legend = 'Fitted Line', lty = 1, col = 'blue')

g <- function(x) { # g function
  dnorm(x = x, mean = mu_hat, sd = sqrt(2))
  #(1 / (2 * sqrt(pi))) * exp((-1 / 4) * (x - mu_hat)^2)
}
R <- function(y, g) { # R function
  abs(y - g)
}

sum(R(y = cvdataY, g = g(cvdataX))) / n # Apparent error

# part (b)
S1 <- cvdataX[1:50]; S2 <- cvdataX[51:100] # Initialize variables
S1_y <- cvdataY[1:50]; S2_y <- cvdataY[51:100]
xbar_1 <- mean(S1); xbar_2 <- mean(S2)

g1 <- function(x, xbar = xbar_1) { # Create g_1X,Y(x)
  (1 / (2 * sqrt(pi))) * exp((-1 / 4) * (x - xbar)^2)
}
g2 <- function(x, xbar = xbar_2) { # Create g_2X,Y(x)
  (1 / (2 * sqrt(pi))) * exp((-1 / 4) * (x - xbar)^2)
}

(1 / n) * sum(R(y = S2_y, g = g1(x = S2)) +
  R(y = S1_y, g = g2(x = S1))) # Partitioned apparent error

### Problem 3
# part (b)
jackknife <- scan(file.choose()) # Load data
b2 <- function(Y) {
  y_bar <- mean(Y)
  sum((Y - y_bar)^4) /
    ((sum((jackknife - y_bar)^2))^2)
}
b2(Y = jackknife) # 0.02669995

b2_minus_j <- function(R_list, j = 0) {
  if (j == 0) { # Remove jth group, if j=0 remove none of the groups
    # Reference: https://stackoverflow.com/questions/1335830/why-cant-rs-ifelse-statements-return-vectors
    Z <- R_list
  } else {
    # Reference: https://stackoverflow.com/questions/652136/how-can-i-remove-an-element-from-a-list
    Z <- R_list[-j]
  }

  # Reference: https://stackoverflow.com/questions/14924935/using-r-convert-data-frame-to-simple-vector
  # Calculate b2 with jth group removed
  Z <- as.vector(unlist(Z), mode = 'numeric')
  z_bar <- mean(Z)
  Z_b2 <- sum((Z - z_bar)^4) / ((sum((Z - z_bar)^2))^2)
  
  return(Z_b2)
}

J <- function(R_j) {
  r <- length(R_j)
  b2_bar <- (1 / r) *
    sum(sapply(1:r,
               function(x) { b2_minus_j(R_list = R_j, j = x) }))
  jackknifed_stat <- r * b2 - (r - 1) * b2_bar
  return(jackknifed_stat)
}

se_jack <- function(Y = jackknife, # Vector of values
                    k = 1) { # Size of each group
  n <- length(Y) # Initialize variables
  r <- n / k
  
  # Reference: https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  # Split Y into r groups R = (r_1, r_2, ..., r_r)^T
  r_groups <- split(Y, cut(seq_along(Y), r, labels = FALSE))
  
  jackknifed_T <- J(R_j = r_groups) # J(T)
  T_stat <- b2_minus_j(R_list = r_groups, j = 0) # T

  numer <- sum( # numerator
    sapply(1:r, function(x) {
      T_j_star <- r * T_stat - (r - 1) *
        b2_minus_j(R_list = r_groups, j = x)
      (T_j_star - jackknifed_T)^2
    })
  )
  denom <- (r * (r - 1)) # denominator
  se_JT <- sqrt(numer / denom)

  return(se_JT)
}

se_jack(Y = jackknife, k = 1) # 0.003714044
se_jack(Y = jackknife, k = 5) # 0.003692711

# part (c)
se_jack(Y = jackknife, k = 10) # 0.003362023
se_jack(Y = jackknife, k = 20) # 0.004758237

plot(density(jackknife))
additional_samples <- rnorm(n = 1e4, mean = 0, sd = 1)
jackknife_2 <- c(jackknife, additional_samples)
plot(density(jackknife_2))
b2(jackknife_2) # 2.517476

se_jack(Y = jackknife_2, k = 1) # 2.653616
se_jack(Y = jackknife_2, k = 5) # 1.186968
se_jack(Y = jackknife_2, k = 10) # 0.8395211
se_jack(Y = jackknife_2, k = 20) # 0.5939255

additional_samples2 <- rnorm(n = 1e3, mean = 0, sd = 1)
jackknife_3 <- c(jackknife, additional_samples2)
plot(density(jackknife_3))
b2(jackknife_3) # 0.2688279

se_jack(Y = jackknife_3, k = 1) # 0.7922754
se_jack(Y = jackknife_3, k = 5) # 0.3549629
se_jack(Y = jackknife_3, k = 10) # 0.2515717
se_jack(Y = jackknife_3, k = 20) # 0.1787098

