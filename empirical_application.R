####################################################################################
# Cross Validation Methods for Kernel Density Estimation: Biased Cross Validation
### data source: https://data.world/chasewillden/fortune-500-companies-2017
####################################################################################


#########################
## R installations and imports
#########################

install.packages("bandwidth")
install_github("cran/kedd")

library(ggplot2)
library(dplyr)
library(np) # Nonparametric Kernel Smoothing Methods for Mixed Data Types
library(KernSmooth)
?KernSmooth
library(remotes) # enable installation of github packages
library(kedd)

# Function to calculate proportions of nulls in each column of a data frame
pct_nulls <- function(data) {
  # Calculate the total number of rows in the data frame
  total_rows <- nrow(data)

  # Calculate the proportion of nulls in each column
  pct_nulls <- sapply(data, function(x) sum(is.na(x)) / total_rows)

  # Create a data frame to store the results
  result_df <- data.frame(
    Column = names(data),
    pct_null = pct_nulls
  )

  return(result_df)
}

# Function to find the optimal parameter h using normal scale rule, taken from in class TD
h_NSR <- function(data){
  # From TD
  # the var() function in R uses Bessel correction by multiplying variance
  # with n/(n-1) so we multiply by the inverse
  n <- length(data)
  return (1.059*sqrt((length(data)-1)*var(data)/n)*(length(data))^{-1/5})
}

#########################
## Access data
#########################

# Set the working directory to a specific folder
setwd("/Users/sam/Documents/TSE/M2/TSE M2 S1/Nonparametric econometric methods/Project/R")

# Replace "path/to/your/directory" with the actual path to the directory you want to set as your working directory.

# Alternatively, download from the web
df <- read.csv("https://query.data.world/s/fpexibocaqatb6x6wghopimatt2ur4?dws=00000", header=TRUE, stringsAsFactors=FALSE)

# Load the CSV file into a data frame
data_folder <- "data"
file_name <- "fortune500.csv"
file_path <- file.path(data_folder, file_name)
df <- read.csv(file_path)

#########################
## Explore data
#########################
head(df)
names(df)
nrow(df)

# get info about nulls per column
pct_nulls <- pct_nulls(df)
pct_nulls

##
revenue_col_value <- "highlights.0.value"
company_id <- "title"
company_revenues <- df %>%
  dplyr::select(c(company_id, revenue_col_value)) %>%
  dplyr::rename(revenue = revenue_col_value, company = company_id)

head(company_revenues)
paste("Shape of df: (", nrow(df),",", length(df), ")")

sapply(company_revenues, class)

### Plot the density of revenues - right skewed
ggplot(company_revenues, aes(x = revenue)) +
  geom_density() +
  labs(x = "Revenue", y = "Density", title = "Density Plot of Revenue")



revenue <- company_revenues$revenue

#########################
## Find an optimal h with NSR method and estimate density function with KDE
##### NEED TO TRANSFORM DATA TO BE SYMMETRIC to use this NSR method
##### because NSR relies on assumption that the true density is close to a
##### Gausisian density. If this is not the case, then the obtained smoothing
##### and corresponding density estimates are not reliable.... so requires
##### the transformation
#########################
log_revenue <- log(company_revenues$revenue)

boxplot(revenue)
boxplot(log_revenue) ## still a bit skewed

NSR_h <- h_NSR(data=revenue)
NSR_h_logged <- h_NSR(data=log_revenue)
paste("Optimal h with NSR using original data: ", NSR_h)
paste("Optimal h with NSR using logged data: ", NSR_h_logged)

# bkde(x, kernel = "normal", canonical = FALSE, bandwidth, gridsize = 401L, range.x, truncate = TRUE)
dens_est_h_NSR <- KernSmooth::bkde(x=log_revenue
                                   ,kernel = "normal"
                                   , canonical = FALSE
                                   , bandwidth=NSR_h)

plot(dens_est_h_NSR, type="l", ylab="Density", xlab="Revenue")


#########################
## Find an optimal h with BCV method and estimate density function with KDE
#########################
?kedd


gaus_bcv <- kedd::h.bcv(x=revenue, whichbcv=1, deriv.order=0, kernel="gaussian")
epan_bcv <- kedd::h.bcv(x=revenue, whichbcv=1, deriv.order=0, kernel="epanechnikov")
triw_bcv <- kedd::h.bcv(x=revenue, whichbcv=1, deriv.order=0, kernel="triweight")
biw_bcv <- kedd::h.bcv(x=revenue, whichbcv=1, deriv.order=0, kernel="biweight")



dens_est_h_gaus_bcv <- KernSmooth::bkde(x=revenue
                                   ,kernel = "normal"
                                   , canonical = FALSE
                                   , bandwidth=gaus_bcv$h)

dens_est_h_gaus_epan_bcv <- KernSmooth::bkde(x=revenue
                                        ,kernel = "epanech"
                                        , canonical = FALSE
                                        , bandwidth=epan_bcv$h)

dens_est_h_gaus_triw_bcv <- KernSmooth::bkde(x=revenue
                                             ,kernel = "triweight"
                                             , canonical = FALSE
                                             , bandwidth=triw_bcv$h)

dens_est_h_gaus_biw_bcv <- KernSmooth::bkde(x=revenue
                                             ,kernel = "biweight"
                                             , canonical = FALSE
                                             , bandwidth=biw_bcv$h)

plot(dens_est_h_NSR, type="l", ylab="Density", xlab="Revenue")
plot(dens_est_h_gaus_bcv, type="l", ylab="Density", xlab="Revenue")
plot(dens_est_h_gaus_epan_bcv, type="l", ylab="Density", xlab="Revenue")
plot(dens_est_h_gaus_triw_bcv, type="l", ylab="Density", xlab="Revenue")
plot(dens_est_h_gaus_biw_bcv, type="l", ylab="Density", xlab="Revenue")



#########################
## Plot densities all together
#########################

# Combine the density estimates into a data frame
densities <- data.frame(
  x = c(dens_est_h_NSR$x, dens_est_h_gaus_bcv$x, dens_est_h_gaus_epan_bcv$x,
        dens_est_h_gaus_triw_bcv$x, dens_est_h_gaus_biw_bcv$x),
  density = c(dens_est_h_NSR$y, dens_est_h_gaus_bcv$y, dens_est_h_gaus_epan_bcv$y,
              dens_est_h_gaus_triw_bcv$y, dens_est_h_gaus_biw_bcv$y),
  method = rep(c("NSR", "Gaussian (bcv)", "Epanechnikov (bcv)", "Triweight (bcv)", "Biweight (bcv)"), each = length(dens_est_h_NSR$x))
)

# Create the ggplot
ggplot(densities, aes(x = x, y = density, color = method)) +
  geom_line() +
  labs(y = "Density", x = "Revenue") +
  theme_minimal()


