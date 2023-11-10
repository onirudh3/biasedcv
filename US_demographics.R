####################################################################################
# Cross Validation Methods for Kernel Density Estimation: Biased Cross Validation
### data source: https://public.opendatasoft.com/explore/dataset/us-cities-demographics/table/
####################################################################################


#########################
## R installations and imports
#########################

# install.packages("bandwidth")
# install_github("cran/kedd")
# install.packages("reshape2")

library(ggplot2)
library(dplyr)
library(tidyr)
library(np) # Nonparametric Kernel Smoothing Methods for Mixed Data Types
library(KernSmooth) 
library(remotes) # enable installation of github packages
library(kedd)
library(reshape2)


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

# access data
data_folder <- "data"
file_name <- "us-cities-demographics.csv"
file_path <- file.path(data_folder, file_name)
df <- read.csv(file_path, header = TRUE, sep = ";")


#########################
## Explore data
#########################
head(df)
names(df)
nrow(df)

# get info about nulls per column
pct_nulls(df)

# Puerto Rico is the only one with NULL values of potential cols of interest
df %>% dplyr::filter(is.na(Number.of.Veterans)) %>% select(State) %>% distinct()
df %>% dplyr::filter(is.na(Foreign.born)) %>% select(State) %>% distinct()

## aggregate to city,state level
df_agg <- df %>%
  filter(State != "Puerto Rico") %>%
  group_by(City, State) %>%
  summarize(
    max_total_population = max(Total.Population),
    sum_number_of_veterans = sum(Number.of.Veterans),
    sum_foreign_born = sum(Foreign.born),
    sum_count = sum(Count)
  ) %>% 
  ungroup()

# check nulls on aggregated df
pct_nulls(df_agg)

# confirm data types
sapply(df_agg, class)

# View cities per state
city_per_state <- df_agg %>% group_by(State) %>% summarize(Cities = n_distinct(City)) %>% arrange(desc(Cities))
city_per_state

# get potential x column names
x_cols <- setdiff(names(df_agg), c("State", "City"))

# Create a grid of density plots
df_agg %>%
  # select(x_cols) %>%
  dplyr::select(all_of(x_cols)) %>%
  tidyr::pivot_longer(everything(), names_to = "Column") %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~ Column, scales = "free") +
  labs(x = "Value", y = "Density") +
  ggtitle("Density Plots of Selected Columns")

# Focus on count of veterans - Notice the distribution is very skewed
density_veterans <-ggplot(df_agg, aes(x = sum_number_of_veterans)) +
  geom_density() +
  labs(x = "Veterans", y = "Density")
density_veterans ### WE DONT want to use this in our report though bc ggplot uses KDE estimation


## Create transformations on veterans variable to try and create a more normal dist
veterans <- df_agg %>%
  select(sum_number_of_veterans) %>% 
  mutate(log_sum_number_of_veterans = log(sum_number_of_veterans))

# veterans data is skewed
boxplot_veterans <- ggplot(veterans, aes(y = sum_number_of_veterans)) +
  geom_boxplot() +
  labs(y = "Total Veterans")
boxplot_veterans

# log of veterans appears symmetric now
boxplot_log_veterans <- ggplot(veterans, aes(y = log_sum_number_of_veterans)) +
  geom_boxplot() +
  labs(y = "Log Total Veterans")
boxplot_log_veterans

#########################
## Find an optimal h with BCV method and estimate density function with KDE on the original data
#########################

# 2 different lists required for the parameter name differences in packages used
kernels_kedd <- c("gaussian", "epanechnikov", "triweight", "biweight")
kernels_kernel_smooth <- c("normal", "epanech", "triweight", "biweight")



# Create structures to hold estimated h values and densities on original data
optimal_h_original <- list()
bcv_densities_original <- list()

# define variable of interest
veteran <- veterans$sum_number_of_veterans

for (i in 1:length(kernels_kedd)) {
  # bcv cross validation
  bcv <- kedd::h.bcv(x = veteran,
                     whichbcv = 1,
                     deriv.order = 0,
                     kernel = kernels_kedd[i])
  
  # store estimated h
  optimal_h_original[[paste0(kernels_kedd[i], "_h")]] <- bcv$h  
  
  # generate density estimate
  density_values <- KernSmooth::bkde(x=veteran
                                       ,kernel = kernels_kernel_smooth[i]
                                       , canonical = FALSE
                                       , bandwidth=optimal_h_original[[paste0(kernels_kedd[i], "_h")]])
  # store estimated density as a df
  bcv_densities_original[[paste0(kernels_kedd[i], "_density")]] <- data.frame(
    x = density_values$x,
    density = density_values$y,
    method = rep(paste0(kernels_kedd[i], " h = ", round(optimal_h_original[[paste0(kernels_kedd[i], "_h")]], 3)), length(density_values$x))
  )
}

# plot densities all together
densities <- bind_rows(bcv_densities_original)
names(densities)

density_estimates_original_plot <- ggplot(densities, aes(x = x, y = density, color = method)) +
  geom_line() +
  labs(x = "Veterans", y = "Density")
density_estimates_original_plot

#########################
## Find an optimal h with NSR and BCV methods and estimate density function 
## with KDE on normal looking data
#########################

# define variable of interest
veteran <- veterans$log_sum_number_of_veterans

# Create structures to hold estimated h values for transformed (noraml) data
optimal_h_normal <- list()
bcv_densities_normal <- list()

#### NSR method
optimal_h_normal[["NSR_h"]] <- h_NSR(data=veteran) 
dens_est_h_NSR <- KernSmooth::bkde(x=veteran
                                   ,kernel = "normal"
                                   , canonical = FALSE
                                   , bandwidth=optimal_h_normal[["NSR_h"]])

bcv_densities_normal[[paste0("NSR_density")]] <- data.frame(
  x = dens_est_h_NSR$x,
  density = dens_est_h_NSR$y,
  method = rep(paste0("NSR h = ", round(optimal_h_normal[["NSR_h"]],3)), length(dens_est_h_NSR$x))
)

# use same bcv method as above except on logged data
for (i in 1:length(kernels_kedd)) {
  # bcv cross validation
  bcv <- kedd::h.bcv(x = veteran,
                     whichbcv = 1,
                     deriv.order = 0,
                     kernel = kernels_kedd[i])
  
  # store estimated h
  optimal_h_normal[[paste0(kernels_kedd[i], "_h")]] <- bcv$h  
  
  # generate density estimate
  density_values <- KernSmooth::bkde(x=veteran
                                     ,kernel = kernels_kernel_smooth[i]
                                     , canonical = FALSE
                                     , bandwidth=optimal_h_normal[[paste0(kernels_kedd[i], "_h")]])
  # store estimated density as a df
  bcv_densities_normal[[paste0(kernels_kedd[i], "_density")]] <- data.frame(
    x = density_values$x,
    density = density_values$y,
    method = rep(paste0(kernels_kedd[i], " h = ", round(optimal_h_normal[[paste0(kernels_kedd[i], "_h")]], 3)), length(density_values$x))
  )
}


# plot densities all together
densities <- bind_rows(bcv_densities_normal)
names(densities)
head(densities)

density_estimates_normal_plot <- ggplot(densities, aes(x = x, y = density, color = method)) +
  geom_line() +
  labs(x = "Veterans", y = "Density")
density_estimates_normal_plot






