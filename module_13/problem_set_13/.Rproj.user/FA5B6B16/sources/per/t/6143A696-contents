library(lattice); library(TeachingDemos); library(graphics) # load libraries

### Q1
univariate <- read.table('Univariate.txt') # load data
univariate <- as.vector(univariate[[1]])
# part (a)
par(mfrow = c(2, 2))
# plot 1 x_i vs. i
n <- length(univariate)
plot(1:n, univariate, main = 'x Against Its Index',
     ylab = 'x', xlab = 'Index')

# plot 2 x_i+1 vs. x_i
plot(univariate[1:(n - 1)], univariate[2:n],
     main = 'x Lagged', ylab = 'x[1:(n-1)]', xlab = 'x[2:n]')

# plot 3 histogram
num_breaks <- 1 + log2(n)
hist(univariate, main = 'Histogram of x', xlab = 'x',
     breaks = )

# plot 4 Q-Q plot
# prob plot
plot(univariate, pnorm(univariate), main = 'Probability Plot',
     xlab = latex2exp::TeX('x_i'), ylab = latex2exp::TeX('P(X $\\leq$ x_i)'))

# q-q
sort_univariate <- sort(univariate)

plot(qnorm((1:n)/n), sort_univariate, main = 'Q-Q Plot',
     xlab = 'Normal(0,1)', ylab = 'x')

# part (b)
par(mfrow = c(1, 2))
plot(sort_univariate, (1:n)/n, pch = 20, # broken-line ECDF
     main = 'Broken-Line ECDF', xlab = 'x',
     ylab = 'ECDF(x)')
lines(sort_univariate, (1:n)/n)

which(sort_univariate > median(sort_univariate))
which(sort_univariate < median(sort_univariate))

left_side <- (1:(n / 2)) / n
right_side <- ((n / 2):1) / n
mountain <- c(left_side, right_side)
plot(sort_univariate, mountain, pch = 20, # mountain plot
     main = 'Mountain Plot', xlab = 'x',
     ylab = 'Folded Percentile')
lines(sort_univariate, mountain)

### Q2
sample_size <- 2e2
mixture <- 0.8
set.seed(0)
norm_0_1 <- rnorm(n = sample_size * mixture, mean = 0, sd = 1)
norm_3_1 <- rnorm(n = sample_size * (1 - mixture), mean = 3, sd = 1)
population_data <- c(norm_0_1, norm_3_1)

par(mfrow = c(2, 2))
plot(density(population_data), main = 'Density of Mixture',
     xlab = 'Values')

bins_rec <- round(1 + log2(sample_size)) # 9 bins
get_break_pts <- function(df = population_data, m = 9) {
  # Initialize variables
  min_value <- range(df)[1]; max_value <- range(df)[2]
  
  bin_width <- (max_value - min_value) / m # bin width
  break_pts <- seq(from = min_value, # calculate break points
                   to = max_value, by = bin_width)
  
  hist(df, breaks = break_pts,
       main = paste0('Histogram of Data (m=', m, ')'),
       xlab = 'Generated Data')
}

break_pts_9 <- get_break_pts(df = population_data, m = 9)

get_break_pts_vec <- Vectorize(get_break_pts, vectorize.args = c('m'))
round(seq(from = 15, to = 50, length.out = 9))
par(mfrow = c(3,3))
get_break_pts_vec(df = population_data, m = round(seq(from = 15, to = 50, length.out = 9)))

get_break_pts_vec(df = population_data, m = c(28, 37))

### Q3
q3_qq <- function(n = 1e2, test_seed = 0) {
  set.seed(test_seed)
  normal_sample <- rnorm(n = n, mean = 0, sd = 1)
  normal_sample <- sort(normal_sample)
  df1 <- data.frame(x = qnorm(p = (1:n)/n), y = normal_sample)
  fit <- lm(formula = y~x, data = df1[1:(n-1),])
  plot(df1$x, df1$y, main = paste0('Q-Q Plot (seed=', test_seed,')'),
       xlab = 'Normal(0,1)', ylab = 'x')
  abline(fit)
}
# light tails -> small values above line, large values below line
q3_qq_vec <- Vectorize(q3_qq, vectorize.args = c('test_seed'))
par(mfrow = c(3, 3))
q3_qq_vec(n = 1e3, test_seed = 9:17)

q3_qq(n = 1e3, test_seed = 4)

### Q4
x <- read.table('multi1.txt', header = TRUE)[[1]] # load data
y <- read.table('multi2.txt', header = TRUE)[[1]]
z <- read.table('multi3.txt', header = TRUE)[[1]]
multiforgraph <- data.frame(x = x, y = y, z = z)

four_plot <- function(variable, var_name) {
  # plot 1 x_i vs. i
  n <- length(variable)
  plot(1:n, variable, main = paste0(var_name, ' Against Its Index'),
       ylab = var_name, xlab = 'Index')
  
  # plot 2 x_i+1 vs. x_i
  plot(variable[1:(n - 1)], variable[2:n],
       main = paste0(var_name, ' Lagged'),
       ylab = paste0(var_name, '[1:(n-1)]'),
       xlab = paste0(var_name, '[2:n]'))
  
  # plot 3 histogram
  # num_breaks <- 1 + log2(n)
  hist(variable, main = paste0('Histogram of ', var_name), xlab = var_name)
  
  # plot 4 Q-Q plot
  # q-q
  sort_variable <- sort(variable)
  plot(qnorm((1:n)/n), sort_variable, main = 'Q-Q Plot',
       xlab = 'Normal(0,1)', ylab = var_name)
  df1 <- data.frame(x = qnorm(p = (1:n)/n), y = sort_variable)
  fit <- lm(formula = y~x, data = df1[1:(n-1),])
  abline(fit)
}

par(mfrow = c(2, 2))
four_plot(variable = x, var_name = 'x')
four_plot(variable = y, var_name = 'y')
four_plot(variable = z, var_name = 'z')

par(mfrow = c(2, 2))
# x-variable Uniform Q-Q
sort_variable <- sort(x)
var_name <- 'x'
plot(qunif((1:n)/n), sort_variable, main = 'Q-Q Plot',
     xlab = 'Uniform(0,1)', ylab = var_name)
df1 <- data.frame(x = qunif(p = (1:n)/n), y = sort_variable)
fit <- lm(formula = y~x, data = df1[1:(n-1),])
abline(fit)

# y-variable Uniform Q-Q
sort_variable <- sort(y)
var_name <- 'y'
y_quantile <- qnorm((1:n)/n, mean(y), sd = sd(y))
plot(y_quantile, sort_variable, main = 'Q-Q Plot',
     xlab = latex2exp::TeX('Normal(\\bar{y}, s)'), ylab = var_name)
df1 <- data.frame(x = y_quantile, y = sort_variable)
fit <- lm(formula = y~x, data = df1[1:(n-1),])
abline(fit)

# z-variable Uniform Q-Q
sort_variable <- sort(z)
var_name <- 'z'
z_quantile <- qnorm((1:n)/n, mean(z), sd = sd(z))
plot(z_quantile, sort_variable, main = 'Q-Q Plot',
     xlab = latex2exp::TeX('Normal(\\bar{z}, s)'), ylab = var_name)
df1 <- data.frame(x = z_quantile, y = sort_variable)
fit <- lm(formula = y~x, data = df1[1:(n-1),])
abline(fit)
dev.off()

splom(x = multiforgraph[1:3])

image(x = multiforgraph)

### Q5
bears <- read.table('Bears.txt', header = TRUE)
bears_sub <- bears[,3:7]

# Stars
stars(bears_sub)

# Chernoff faces
faces(bears_sub)

# Parallel Coordinates
parallel(bears_sub)
