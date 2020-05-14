### Q4 part(b)
orthogonal <- read.csv('Orthogonal.txt', header = FALSE)
orthogonal <- orthogonal[[1]]
hist(orthogonal)
plot(1:length(orthogonal), orthogonal)

q0 <- function(x) {
  (pi / 2)^(-1/2)
}

q1 <- function(x) {
  x * ((pi / 8)^(-1 / 2))
}

q2 <- function(x) {
  ((x^2) - (1 / 4)) / (sqrt(pi / 32))
}

q3 <- function(x) {
  ((x^3) - (x / 2)) / sqrt(pi / 128)
}

c0 <- mean(q0(orthogonal))
c1 <- mean(q1(orthogonal))
c2 <- mean(q2(orthogonal))
c3 <- mean(q3(orthogonal))

f <- function(x) {
  # c0 * q0(x) + c1 * q1(x) + c2 * q2(x)
  c0 * q0(x) + c1 * q1(x) + c2 * q2(x) + c3 * q3(x)
}

xs <- seq(-1.5, 1.5, length.out = 1e2)
plot(xs, f(xs), type = 'l', ylim = c(0, 1.5),
     main = 'Graph of Estimated Function and True Function',
     xlab = 'x', ylab = 'f(x)')
lines(xs, dnorm(xs, mean = 0, sd = 0.3), col = 'red')
legend("topright", legend = c('Estimated Function', 'True Function'),
       col = c('black', 'red'), lty = c(1,1))

### Q5
# part (a)
f <- function(x) {
  if ((x >= -1) & (x <= 0)) {
    (x + 1) + (x + 1)^3
  } else if ((x > 0) & (x <= 1)) {
    4 + (x - 1) + (x - 1)^3
  } else {
    0
  }
}

xs <- seq(-1.1, 1.1, length.out = 1e4)
f_vec <- Vectorize(f)
plot(xs, f_vec(xs), type = 'l')

# part (b)
s0 <- function(x) {
  (1 / 12) * ((x - 1)^3) + (5 / 12) * (x - 1) + (2 - x)
}

s1 <- function(x) {
  (1 /12) * ((3 - x)^3) + (1 / 3) * (x - 2) + (5 / 12) * (3 - x)
}

s2 <- function(x) {
  (1 / 4) * (x - 3) + (1 / 3) * (4 - x)
}

s0(1)
s0(2) == s1(2)
s1(3) == s2(3)
s2(4)

b = 0.01
x = c(seq(1,2,b), seq(2,3,b), seq(3,4,b))
y = c(s0(seq(1,2,b)), s1(seq(2,3,b)), s2(seq(3,4,b)))
plot(x,y, type = 'l')
abline(v = 1, col='red')
abline(v = 2, col='red')
abline(v = 3, col='red')
abline(v = 4, col='red')
