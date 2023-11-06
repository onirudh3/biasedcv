
# Libraries ---------------------------------------------------------------

library(np)
library(KernSmooth)

# Stuff -------------------------------------------------------------------

# Data
df <- data.frame(x = rnorm(100, mean = 0, sd = 1))

# Bandwidth with normal scale rule with a Gaussian kernel
n <- length(df$x)
hNS <- 1.059 * sqrt((n - 1) * var(df$x) / n) * (length(df$x)) ^ {-1 / 5}

# Histogram
hist(df$x, probability = T)
rug(df$x)
lines(bkde(df$x, bandwidth = hNS), col = 'yellow', lwd = 4, lty = 2)

