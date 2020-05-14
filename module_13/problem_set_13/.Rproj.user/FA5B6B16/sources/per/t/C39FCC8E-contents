library(ggplot2)

mystery_data <- data.frame(
  x = c(2, 0, 1, 6, 1, 0, 1, 0, 2, 0, 4, 2, 1, 1, 3, 2, 3, 0, 5, 1, 4, 1, 3, 4, 0, 1, 0, 3, 0, 3)
) 
a <- ggplot(data = mystery_data, aes(sample = x))+
  stat_qq(distribution = stats::qpois,
          dparams = list(lambda = mean(mystery_data$x))) +
  stat_qq_line(distribution = stats::qpois,
               dparams = list(lambda = mean(mystery_data$x))) +
  ggtitle("Poisson (lambda=mean)")

print(a)

ggplot(data = df,
       mapping = aes(sample = Freq)) + 
  stat_qq(distribution = stats::qpois,
          dparams = list(lambda = mean(df$Freq))) + 
  geom_step(data = data.frame(x = 6:16,
                              Freq = 6:16),
            mapping = aes(x = x,
                          y = Freq),
            colour = "red",
            alpha = 0.5) 