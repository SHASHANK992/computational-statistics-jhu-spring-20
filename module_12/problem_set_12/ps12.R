### Q1
xs <- scan(file.choose()) # Load data
ys <- scan(file.choose())
# plot(xs,ys)

df <- data.frame(x = xs, y = ys)
df_sort <- df[order(df$x),] # sort values by x_i's
rownames(df_sort) <- NULL # reset indices

# part (a)
### Plot the CVRSS_k(s_hat_k) versus k and an explanation
### of the chosen span.
odd_values <- function(k) { 2 * k + 1 } # calculate odd values
k_seq <- 1:11
k_candidates <- odd_values(k_seq) # test spans

### The s_hat function calculates the value of s_hat for a given
### span k. It utlizes the formula in eq. 11.5 of the textbook.
### It requires an argument of both the x-value and the i'th row
### index for that x-value.
s_hat <- function(obs_data = df, x, i = 3, k = 3) {
  n <- nrow(obs_data) # Initialize variables
  j_min <- max(c(i - (k - 1) / 2, 1)); j_max <- min(c(i + (k - 1) / 2, n))
  
  s_hat_k <- mean(obs_data[j_min:j_max,'y'])
  
  return(s_hat_k)
}
s_hat_vec <- Vectorize(s_hat, vectorize.args = c('x', 'i'))

test_row <- df[10,]
s_hat(obs_data = df, x = test_row[,1], i = integer(rownames(test_row)), k = 3)

trunc_col <- function(j, k, n) {
  # for j'th row, which columns should be non-zero?
  k_half <- (k - 1) / 2
  # until row k_half + 1
  if (j <= k_half) {
    
  } else if ((k_half < j) & (j <= (n - k_half))) {
    # when j between k_half + 1 and n - k_half
    # e.g., n=100, j must be between 3 and 98 for k = 5
    
  } else {
    # when j > n - k_half
    
  }
  j
}

truncate_S <- function(n, k) {
  NA_mat <- matrix(data = NA, nrow = n, ncol = n)
  
  for(i in 1:n) {
    NA_mat[i,]
  }
  
}

cvrss <- function(obs_data = df, k) {
  n <- nrow(obs_data) # Initialize variables
  x_vec <- obs_data[,'x']; y_vec <- obs_data[,'y']
  
  
  
  mean(
    ((y_vec - s_hat_vec(obs_data = obs_data, x = x_vec, i = 1:n, k = k)) /
       (1 - S_ii))^2
  )
}


















