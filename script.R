
# Libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(np)
library(KernSmooth) # for bkde()

# For library "kedd"
# devtools::install_github("cran/kedd")
library(kedd) # for h.bcv()


# Fortune 500 companies data ----------------------------------------------

# Data
df <- read.csv("https://query.data.world/s/fpexibocaqatb6x6wghopimatt2ur4?dws=00000",
               header = TRUE, stringsAsFactors = FALSE)

# Our variable of interest is profits as percentage of total assets tables.2.data.1.value
# Other interesting variables are profits highlights.2.value, and profits as percentage of revenue tables.2.data.0.value
df <- subset(df, select = c("title", "tables.2.data.1.value"))

# Rename columns
df <- df %>%
  rename("company" = "title", "profit_perc" = "tables.2.data.1.value")

# Drop NA rows
df <- na.omit(df)

# Histogram
hist(df$profit_perc, probability = T)
rug(df$profit_perc)


# Bandwidths --------------------------------------------------------------

# Normal scale rule bandwidth with a Gaussian kernel
n <- length(df$profit_perc)
NSR <- 1.059 * sqrt((n - 1) * var(df$profit_perc) / n) * n ^ {-1 / 5}

# BCV with Gaussian kernel
gaus_bcv <- h.bcv(x = df$profit_perc, kernel = "gaussian")$h

# BCV with Epanechnikov kernel
epan_bcv <- h.bcv(x = df$profit_perc, kernel = "epanechnikov")$h

# BCV with Biweight kernel
biw_bcv <- h.bcv(x = df$profit_perc, kernel = "biweight")$h

# BCV with Triweight kernel
triw_bcv <- h.bcv(x = df$profit_perc, kernel = "triweight")$h


# Plot with comparison of estimators --------------------------------------

hist(df$profit_perc, probability = T, ylim = c(0, 0.25),
     main = "Kernel estimators for various bandwidths",
     xlab = "Profits as percentage of total assets")
rug(df$profit_perc)
lines(bkde(df$profit_perc, bandwidth = NSR), col = "#5d8aa8", lwd = 1.5,
      lty = 1)
lines(bkde(df$profit_perc, bandwidth = gaus_bcv), col = "#e32636", lwd = 1.5,
      lty = 1)
lines(bkde(df$profit_perc, bandwidth = epan_bcv), col = "#ffbf00", lwd = 1.5,
      lty = 1)
lines(bkde(df$profit_perc, bandwidth = biw_bcv), col = "#00ffff", lwd = 1.5,
      lty = 1)
lines(bkde(df$profit_perc, bandwidth = triw_bcv), col = "#00ff1f", lwd = 1.5,
      lty = 1)
legend("topright", legend = c("Gaussian NSR", "Gaussian BCV", "Epanechnikov BCV",
                              "Biweight BCV", "Triweght BCV"),
       col = c("#5d8aa8", "#e32636", "#ffbf00", "#00ffff", "#00ff1f"), lty = 1,
       cex = 0.8, lwd = 3, x.intersp = 1, text.width = 5)

