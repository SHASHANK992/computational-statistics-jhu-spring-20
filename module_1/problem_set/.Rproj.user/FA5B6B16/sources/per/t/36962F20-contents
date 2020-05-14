### Jared Yu
### Computational Statistics
### Problem Set 1

# Problem 1
favorite <- read.table(file = 'favorite.data', header = FALSE)
# part a
mean(favorite[,1])

# part b
median(favorite[,1])

# part c
sd(favorite[,1])

# part d
min(favorite[,1])

# part e
max(favorite[,1])

# part f
hist(favorite[,1], main = 'Histogram of Favorite Data', xlab = 'Value of Favorite')

# Problem 2
set.seed(1)
random_values <- rnorm(n = 10000, mean = 0, sd = 1)

# part a
hist(random_values, main = 'Histogram of Random Standard Normal Values', xlab = 'Value of Random Data')

# part b
mean(random_values); median(random_values); sd(random_values)

# Problem 3
a <- seq(from = 5, to = 160, by = 5); b <- seq(from = 87, to = 56); d <- a * b

# part a
d[15:17]

# part b
d[d>2000]

# part c
length(d[d>6000])

# Problem 4
sum_perfect_squares <- function(x) {
  sum(which(sqrt(1:x) %% 1 == 0))
}
# part a
sum_perfect_squares(100)

# part b
sum_perfect_squares(100000)

# Problem 5
find_perfect_squares <- function(x) {
  which(sqrt(1:x) %% 1 == 0)
}

# part a
find_perfect_squares(500)

# part b
perfect_squares_matrix <- matrix(find_perfect_squares(100000), ncol = 4)
sink('perfect_squares_matrix.txt')
perfect_squares_matrix
sink()

# part c
perfect_squares_matrix[15, 3]

# Problem 6
x <- seq(from = -pi, to = pi, length.out = 50)
y1 <- sin(x)
y2 <- cos(x)

# part a
plot(y1, x, main = 'Part a')

# part b
plot(y2, x, type = 'l', main = 'Part b')

# part c
plot(y2, x, type = 'l', main = 'Part c')
abline(a = 1, b = -1/3)











