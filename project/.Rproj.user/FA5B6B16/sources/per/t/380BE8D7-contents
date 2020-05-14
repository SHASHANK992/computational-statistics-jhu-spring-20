### dataset too small
diabetes <- read.csv('diabetes.csv')
par(mfrow = c(2,3))
# pregnancies
m_vec <- seq(10, 25, length.out = 6)
sapply(m_vec, function(x) {
  hist_est(var1 = diabetes$Pregnancies, m = x)
  lines(density(diabetes$Pregnancies))
})

# glucose
par(mfrow = c(2,2))
# m_vec <- round(seq(30, 80, length.out = 9))
glucose_m <- c(8, 42, 68, 80)
sapply(glucose_m, function(x) {
  hist_est(var1 = diabetes$Glucose, m = x)
  lines(density(diabetes$Glucose))
})

# 
par(mfrow = c(3,3))
m_vec <- round(seq(10, 50, length.out = 9))
m_vec <- round(seq(30, 60, length.out = 9))
# glucose_m <- c(8, 42, 68, 80)
sapply(m_vec, function(x) {
  hist_est(var1 = diabetes$BloodPressure, m = x)
  lines(density(diabetes$BloodPressure))
})

### steps
# identify a suitable number of bins (possibly more than one set)
# parameters: binwidth, kernel function
# cross-validation: pseudo-likelihood, UCV
# Plug-in methods: Sheather-Jones/Silverman; package, closed-form solution, Monte Carlo
# Maximal smoothing principle: Terrell's method
# Idea: Fit first with kernel/S-J, scale for all the alternate kernels

# Blood pressure is intuitively related to diabetes
# at least bimodal, no typical pdf's

# compare default hist with chosen hist
hist(diabetes$BloodPressure, main = 'Histogram of Blood Pressure',
     xlab = 'Diastolic blood pressure (mm Hg)', freq = FALSE,
     ylim = c(0, 0.035))
lines(density(diabetes$BloodPressure))

# Fig. 13
par(mfrow = c(3,3)) # Plot all the variables find ones that maybe worth investigating
var_description <- c(
  'Number of times pregnant',
  'Plasma glucose concentration a 2 hours in an oral glucose tolerance test',
  'Diastolic blood pressure (mm Hg)',
  'Triceps skin fold thickness (mm)',
  '2-Hour serum insulin (mu U/ml)',
  'Body mass index (weight in kg/(height in m)^2)',
  'Diabetes pedigree function',
  'Age (years)',
  'Class variable (0 or 1)'
)
for (i in 1:ncol(diabetes)) {
  hist(diabetes[,i], freq = FALSE,
       main = paste0('Histogram of ', colnames(diabetes)[i]),
       xlab = var_description[i])
  lines(density(diabetes[,i]))
}

# Insulin has a more interesting kde
hist(diabetes$Insulin, freq = FALSE, ylim = c(0,0.01))
lines(density(diabetes$Insulin))

par(mfrow = c(3,3))
m_vec <- round(seq(10, 50, length.out = 9))
# m_vec <- round(seq(30, 60, length.out = 9))
# glucose_m <- c(8, 42, 68, 80)
sapply(m_vec, function(x) {
  hist_est(var1 = diabetes$Insulin, m = x)
  lines(density(diabetes$Insulin))
})

# skin thickness
hist(diabetes$SkinThickness, freq = FALSE)
lines(density(diabetes$SkinThickness))

par(mfrow = c(3,3))
m_vec <- round(seq(10, 50, length.out = 9))
# m_vec <- round(seq(30, 60, length.out = 9))
# glucose_m <- c(8, 42, 68, 80)
sapply(m_vec, function(x) {
  hist_est(var1 = diabetes$SkinThickness, m = x)
  lines(density(diabetes$SkinThickness))
})

### Insulin
# can try to overlap KDE for the two groups
diab_pos <- diabetes[diabetes$Outcome == 1,]
diab_neg <- diabetes[diabetes$Outcome == 0,]

plot(density(diab_neg$Insulin), main = 'K.D.E. of Insulin')
lines(density(diab_pos$Insulin), col = 'red')
legend("topright", legend = c('Positive', 'Negative'),
       col = c('red', 'black'), lty = c(1,1))

hist(diabetes$Insulin, main = 'Default Histogram for Insulin',
     xlab = '2-Hour serum insulin (mu U/ml)', freq = FALSE,
     ylim = c(0,0.01))
lines(density(diabetes$Insulin))

par(mfrow = c(3,3))
m_vec <- round(seq(20, 60, length.out = 9))
# m_vec <- round(seq(30, 60, length.out = 9))
# glucose_m <- c(8, 42, 68, 80)
sapply(m_vec, function(x) {
  hist_est(var1 = diabetes$Insulin, m = x)
  lines(density(diabetes$Insulin))
})

insulin_sub <- diabetes[diabetes$Insulin < 400,]
sapply(m_vec, function(x) {
  hist_est(var1 = insulin_sub$Insulin, m = x)
  lines(density(insulin_sub$Insulin))
})

# uninteresting

### BMI
plot(density(diab_pos$BMI), main = 'K.D.E. of BMI', col = 'red')
lines(density(diab_neg$BMI), col = 'black')
legend("topright", legend = c('Positive', 'Negative'),
       col = c('red', 'black'), lty = c(1,1))

hist(diabetes$BMI, main = 'Default Histogram for BMI',
     xlab = 'Body mass index (weight in kg/(height in m)^2)', freq = FALSE)
lines(density(diabetes$BMI))

par(mfrow = c(3,3))
m_vec <- round(seq(50, 120, length.out = 9))
# m_vec <- round(seq(30, 60, length.out = 9))
# glucose_m <- c(8, 42, 68, 80)
sapply(m_vec, function(x) {
  hist_est(var1 = diabetes$BMI, m = x)
  lines(density(diabetes$BMI))
})

### Age
plot(density(diab_neg$Age), main = 'K.D.E. of Age', col = 'black')
lines(density(diab_pos$Age), col = 'red')
legend("topright", legend = c('Positive', 'Negative'),
       col = c('red', 'black'), lty = c(1,1))

par(mfrow = c(3,3))
m_vec <- round(seq(25, 40, length.out = 9))
sapply(m_vec, function(x) {
  hist_est(var1 = diabetes$Age, m = x)
  lines(density(diabetes$Age))
})

# Fig. 14
par(mfrow = c(1,2))
hist(diabetes$Age, main = 'Default Histogram for Age',
     xlab = 'Age (years)', freq = FALSE)
lines(density(diabetes$Age))
hist_est(var1 = diabetes$Age, m = 32)
lines(density(diabetes$Age))

### Cross-Validation
# pseudo-likelihood
# Use S-J as cheat bandwidth?
silverman <- function(x = diabetes$Age) {
  n <- length(x)
  sigma_hat <- sd(x)
  h <- ((4 / (3 * n))^(1 / 5)) * sigma_hat
  return(h)
}
silverman_age <- silverman(x = diabetes$Age)
sheather_jones_age <- locfit::sjpi(x = diabetes$Age, a = silverman_age)[1]
cross_validation <- function(var1 = diabetes$Age, h = 0.3, K = normal) {
  n <- length(var1)
  i_not_j_vec <- 1:n
  cv_matrix <- matrix(NA, nrow = length(var1))
  K_i_not_j <- matrix(NA, nrow = length(var1))
  
  for (i in 1:n) {
    for (j in i_not_j_vec[-i]) {
      K_i_not_j[j] <- K((var1[i] - var1[j]) / h)
    }
    cv_matrix[i] <- (1 / (h * (n - 1))) * sum(K_i_not_j[,1], na.rm = TRUE)
    K_i_not_j <- matrix(NA, nrow = length(var1)) # reset to NA's
  }
  return(cv_matrix)
}
h_seq_age <- seq(sheather_jones_age - 2, sheather_jones_age + 2, length.out = 1e2)
cv_results_age <- sapply(h_seq_age, function(x) {
  cv_output <- cross_validation(var1 = diabetes$Age, h = sheather_jones_age, K = normal)
  prod(cv_output) # take the product of all the C-V values
})

plot(h_seq_age, cv_results_age, type = 'l',
     main = 'Pseudo-likelihood vs Bandwidth',
     xlab = 'Bandwidth h', ylab = 'pseudo-likelihood')

UCV <- function(x, h = sheather_jones_age) {
  # First calculate R(f-hat)
  n <- length(x)
  LHS <- (sqrt(pi) / (2 * pi * n * h))
  
  upper_lower_triangle <- outer(X = x, Y = x, FUN = function(u, v) {
    exp((-1 / (4 * (h^2))) * ((u - v)^2))
  })
  diag(upper_lower_triangle) <- 0
  RHS <- (sqrt(pi) / (2 * pi * (n^2) * h)) * sum(upper_lower_triangle)
  
  R_f <- LHS + RHS # R(f-hat)

  # Calculate the 2*E[f-hat] term using cross-validation and combine with R_f
  ucv <- R_f - 2 * mean(cross_validation(var1 = x, h = h, K = normal))
  
  return(ucv)
}

# test values about S-J
ucv_seq_age <- seq(sheather_jones_age - 3, sheather_jones_age + 3, length.out = 20)
ucv_values_age <- sapply(ucv_seq_age, function(x) {
  UCV(x = diabetes$Age, h = x)
})

ucv_h_age <- ucv_seq_age[which.min(ucv_values_age)] # 0.2517184
ucv_values_age[which.min(ucv_values_age)] # -0.0787289

par(mfrow = c(1,2))
plot(ucv_seq_age, ucv_values_age, type = 'l', # plot UCV results
     main = 'UCV vs. Bandwidth', xlab = 'Bandwidth h', ylab = 'UCV')
hist_est(var1 = diabetes$Age, m = 32)
D <- range(diabetes$Age)
x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
# Normal
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                              h = ucv_h_age, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 1, col = 'red', lwd = 1) # add KDE
legend("topright", legend = c('UCV'), col = 'red', lty = 1)


# Sheather-Jones / Silverman
silverman <- function(x = diabetes$Age) {
  n <- length(x)
  sigma_hat <- sd(x)
  h <- ((4 / (3 * n))^(1 / 5)) * sigma_hat
  return(h)
}
silverman_age <- silverman(x = diabetes$Age) # 3.298613
locfit::sjpi(x = diabetes$Age, a = silverman_age)[1] # 1.672771
bw.SJ(x = diabetes$Age) # 1.007274
MASS::width.SJ(x = diabetes$Age) # 4.029095

# Monte Carlo Integration
### The MC_SJ function uses Monte Carlo Integration to estimate
### the roughness of the second derivative of the K.D.E. f-hat.
### The function assumes that the kernel function is the Normal
### kernel. It also uses the Silverman bandwidth for a pilot
### bandwidth. There are parameters to determine the approximation
### for the infinite bounds (A=1e4) and the number of MC iterations
### to try (n_MC). The function will output the S-J estimate.
MC_SJ <- function(x = diabetes$Age, h_0 = silverman_age, A = 1e4, n_MC = 1e3) {
  n <- length(x) # Initialize variables
  set.seed(0)
  u_sample <- runif(n = n_MC, min = -A, max = A)

  # Calculate g(x), where X~Unif[-A,A]
  MC_est <- sapply(u_sample[1], function(y) {
    squared_sum <- outer(X = x, Y = x, FUN = function(u, v) {
        (1 / (2 * pi)) *
          exp(-0.5 * ((y - u) / h_0)^2) * (1 - ((y - u) / h_0)^2) *
          exp(-0.5 * ((y - v) / h_0)^2) * (1 - ((y - v) / h_0)^2)
      })
    squared_sum <- sum(squared_sum)
    (2 * A) * squared_sum
    }
  )

  MC_avg <- mean(MC_est) # Take avg of MC results

  # Calculate S-J estimate with MC
  roughness_estimate <- (1 / ((n^2) * (h_0^6))) * MC_avg
  sheather_jones_MC <- (R_K / (n * roughness_estimate))^(1 / 5)
  
  return(sheather_jones_MC)
}

sheather_jones_MC_age <- MC_SJ(x = diabetes$Age, h_0 = silverman_age, A = 1e4, n_MC = 1e3)

# Maximal smoothing principle
maximal_smoothing_principle <- function(x) {
  n <- length(x)
  sigma_hat <- sd(x)
  R_K <- 1 / (2 * sqrt(pi))
  terrell_h <- 3 * ((R_K / (35 * n))^(1 / 5)) * sigma_hat
  return(terrell_h)
}

terrell_age <- maximal_smoothing_principle(x = diabetes$Age) # 3.562298

# Bandwidth scaling
### The scale_bandwidth function is used to scale a bandwidth from
### kernel function K to a bandwidth from kernel function L.
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
  scale_bandwidth(h_K = sheather_jones_MC_age,
                  delta_K = delta_normal, delta_L = x))
scale_vec # 2.636304 3.354064 3.684642 3.973440 4.512032

# Compare results
par(mfrow = c(1,2))
hist_est(var1 = diabetes$Age, m = 32)
D <- range(diabetes$Age)
x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
# Normal
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                               h = sheather_jones_MC_age, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 1, col = 'black', lwd = 2) # add KDE
# Uniform
f_hat_pts_uniform <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                               h = scale_vec[1], K = uniform)
lines(x_pts, f_hat_pts_uniform, lty = 1, col = 'blue') # add KDE
# Epanechnikov
f_hat_pts_epanechnikov <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                                    h = scale_vec[2], K = epanechnikov)
lines(x_pts, f_hat_pts_epanechnikov, lty = 2, col = 'red') # add KDE
legend("topright", legend = c('Normal', 'Uniform', 'Epanechnikov'),
       col = c('black', 'blue', 'red'), lty = c(1,1,2), lwd = c(2,1,1))
# Triangle
f_hat_pts_triangle <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                                 h = scale_vec[3], K = triangle)
lines(x_pts, f_hat_pts_triangle, lty = 1, col = 'red') # add KDE
# Biweight
f_hat_pts_biweight <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                                 h = scale_vec[4], K = biweight)
lines(x_pts, f_hat_pts_biweight, lty = 2, col = 'blue') # add KDE
# Triweight
f_hat_pts_triweight <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                                 h = scale_vec[5], K = triweight)
lines(x_pts, f_hat_pts_triweight, lty = 3, col = 'green') # add KDE
legend("topright", legend = c('Normal', 'Triangle', 'Biweight', 'Triweight'),
       col = c('black', 'red', 'blue', 'green'), lty = c(1,1,2,3), lwd = c(2,1,1,1))

# Plot Normal + bandwidths
hist_est(var1 = diabetes$Age, m = 32)
D <- range(diabetes$Age)
x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
# Normal
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                              h = sheather_jones_MC_age, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 1, col = 'black', lwd = 1) # add KDE
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                              h = silverman_age, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 2, col = 'red', lwd = 1) # add KDE
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = diabetes$Age,
                              h = terrell_age, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 3, col = 'blue', lwd = 1) # add KDE
legend("topright", legend = c('Sheather-Jones', 'Silverman', 'Terrell'),
       lty = c(1,2,3), col = c('black', 'red', 'blue'))

# Explore difference
# can try to overlap KDE for the two groups
diab_pos <- diabetes[diabetes$Outcome == 1,]
diab_neg <- diabetes[diabetes$Outcome == 0,]

plot(density(diab_neg$Age), main = 'K.D.E. of Age (Positive and Negative for Diabetes)',
     xlab = 'Age (years)')
lines(density(diab_pos$Age), col = 'red')
legend("topright", legend = c('Positive', 'Negative'),
       col = c('red', 'black'), lty = c(1,1))

### Dataset 2
wine <- read.csv('winequality-red.csv')
par(mfrow = c(4,3))
var_names <- sub("[.]", " ", sub("[.]", " ", colnames(wine)))
var_description <- c(
  'fixed acidity',
  'amount of acetic acid',
  'citric acid',
  'amount of sugar after fermentation',
  'amount of salt',
  'free SO2',
  'total SO2',
  'density of wine',
  'acidic level pH',
  'sulphates',
  'percent alcohol',
  'quality score'
)
for (i in 1:ncol(wine)) {
  hist(wine[,i], freq = FALSE,
       main = paste0('Histogram of ', var_names[i]),
       xlab = var_description[i])
  lines(density(wine[,i]))
}

# 1
par(mfrow = c(3,3))
m_vec <- round(seq(10, 30, length.out = 9))
sapply(m_vec, function(x) {
  hist_est(var1 = wine$volatile.acidity, m = x)
  lines(density(wine$volatile.acidity))
})

par(mfrow = c(1,2))
hist(wine$volatile.acidity, main = 'Default Histogram for Volatile Acidity',
     xlab = 'Amount of Acetic Acid', freq = FALSE)
lines(density(wine$volatile.acidity))
hist_est(var1 = wine$volatile.acidity, m = 30)
lines(density(wine$volatile.acidity))

# 2
silverman_volatile_acid <- silverman(x = wine$volatile.acidity)
sheather_jones_vol_acid <- locfit::sjpi(x = wine$volatile.acidity, a = silverman_volatile_acid)[1]
h_seq_vol_a <- seq(sheather_jones_vol_acid - 0.2, sheather_jones_vol_acid + 0.1, length.out = 50)
cv_results_vol_a <- sapply(h_seq_vol_a, function(x) {
  prod(
    cross_validation(var1 = wine$volatile.acidity, h = x, K = normal)
  )
})

plot(h_seq_vol_a, cv_results_vol_a, type = 'l',
     main = 'Pseudo-likelihood vs. Bandwidth',
     xlab = 'Bandwidth h', ylab = 'pseudo-likelihood')
h_seq_vol_a[which.min(cv_results_vol_a)]

# UCV
ucv_seq_vol_a <- seq(sheather_jones_vol_acid - 3, sheather_jones_vol_acid + 3, length.out = 20)
ucv_values_vol_a <- sapply(ucv_seq_vol_a, function(x) {
  UCV(x = wine$volatile.acidity, h = x)
})

par(mfrow = c(1,2))
plot(ucv_seq_vol_a, ucv_values_vol_a, type = 'l', main = 'UCV vs. Bandwidth',
     xlab = 'Bandwidth h', ylab = 'UCV')
ucv_h_vol_a <- ucv_seq_vol_a[which.min(ucv_values_vol_a)] # 0.1907263
ucv_values_vol_a[which.min(ucv_values_vol_a)] # -1.468854

hist_est(var1 = wine$volatile.acidity, m = 30)
D <- range(wine$volatile.acidity)
x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
# Normal
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                              h = ucv_h_vol_a, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 1, col = 'red', lwd = 1) # add KDE
legend("topright", legend = c('UCV'), col = 'red', lty = 1)

# 3
silverman_volatile_acid <- silverman(x = wine$volatile.acidity) # 0.04337265
sheather_jones_vol_acid <- locfit::sjpi(x = wine$volatile.acidity, # 0.03283153
                                        a = silverman_volatile_acid)[1]
# sheather_jones_MC_vol_acid <- MC_SJ(x = wine$volatile.acidity,
#                                     h_0 = silverman_volatile_acid, A = 1e4, n_MC = 1e3)
sheather_jones_MC_vol_acid <- 0.02150129 # forced result using ad-hoc trick
bw.SJ(x = wine$volatile.acidity) # 0.03162348
MASS::width.SJ(x = wine$volatile.acidity) # 0.1264939

# 4
terrell_volatile_acid <- maximal_smoothing_principle(x = wine$volatile.acidity) # 0.04683978

# 5
scale_vec <- sapply(delta_vec, function(x)
  scale_bandwidth(h_K = sheather_jones_MC_vol_acid,
                  delta_K = delta_normal, delta_L = x))
scale_vec # 0.03741347 0.04759965 0.05229110 0.05638961 0.06403312

# 6
# Compare results
par(mfrow = c(1,2))
hist_est(var1 = wine$volatile.acidity, m = 30)
D <- range(wine$volatile.acidity)
x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
# Normal
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                              h = sheather_jones_MC_vol_acid, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 1, col = 'black', lwd = 2) # add KDE
# Uniform
f_hat_pts_uniform <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                               h = scale_vec[1], K = uniform)
lines(x_pts, f_hat_pts_uniform, lty = 1, col = 'blue') # add KDE
# Epanechnikov
f_hat_pts_epanechnikov <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                                    h = scale_vec[2], K = epanechnikov)
lines(x_pts, f_hat_pts_epanechnikov, lty = 2, col = 'red') # add KDE
legend("topright", legend = c('Normal', 'Uniform', 'Epanechnikov'),
       col = c('black', 'blue', 'red'), lty = c(1,1,2), lwd = c(2,1,1))
# Triangle
f_hat_pts_triangle <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                                h = scale_vec[3], K = triangle)
lines(x_pts, f_hat_pts_triangle, lty = 1, col = 'red') # add KDE
# Biweight
f_hat_pts_biweight <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                                h = scale_vec[4], K = biweight)
lines(x_pts, f_hat_pts_biweight, lty = 2, col = 'blue') # add KDE
# Triweight
f_hat_pts_triweight <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                                 h = scale_vec[5], K = triweight)
lines(x_pts, f_hat_pts_triweight, lty = 3, col = 'green') # add KDE
legend("topright", legend = c('Normal', 'Triangle', 'Biweight', 'Triweight'),
       col = c('black', 'red', 'blue', 'green'), lty = c(1,1,2,3), lwd = c(2,1,1,1))
dev.off()

# 7
# Plot Normal + bandwidths
hist_est(var1 = wine$volatile.acidity, m = 30)
D <- range(wine$volatile.acidity)
x_pts <- seq(D[1], D[2], length.out = 1e2) # sequence of x values for KDE
# Normal
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                              h = sheather_jones_MC_vol_acid, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 1, col = 'black', lwd = 1) # add KDE
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                              h = silverman_volatile_acid, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 2, col = 'red', lwd = 1) # add KDE
f_hat_pts_normal <- f_hat_vec(x = x_pts, X_i = wine$volatile.acidity,
                              h = terrell_volatile_acid, K = normal)
lines(x_pts, f_hat_pts_normal, lty = 3, col = 'blue', lwd = 1) # add KDE
legend("topright", legend = c('Sheather-Jones', 'Silverman', 'Terrell'),
       lty = c(1,2,3), col = c('black', 'red', 'blue'))

# volatile acidity
good_wine <- wine[wine$quality >= 7,]
bad_wine <- wine[wine$quality <= 6,]
plot(density(good_wine$volatile.acidity),
     main = 'K.D.E. of Volatile Acidity', col = 'red',
     xlab = 'Amount of Acetic Acid')
lines(density(bad_wine$volatile.acidity), col = 'black')
legend("topright", legend = c('Good Wine (7+)', 'Bad Wine (<7)'),
       col = c('red', 'black'), lty = c(1,1))
