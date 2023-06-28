
######################################
#  Portfolio Introduction
######################################

# Set working directory
workdir <-readline(prompt="working directory?")
setwd(workdir)

# Read the csv with AAPL and K stocks
monthly_prices <- read.csv("data.csv", row.names = 1, 
                           stringsAsFactors = FALSE,
                           colClasses=c("character", "double", "double")) 
# Sort the data by date
monthly_prices <- monthly_prices[order(
  as.Date(row.names(monthly_prices), format="%d-%b-%y")),] 

# Compute Continuous Returns
monthly_returns <- data.frame(apply(monthly_prices
                                    ,2,function(x) diff(log(x)) ))
# Note: the "apply" function utilizes a function over the rows (1) or columns (2)

# Sneak peak to the first five lines in our returns table
head(monthly_returns)



### Computing Descriptive Statistics for Stocks
# Mean Returns
mean_returns <- apply(monthly_returns, 2,mean) # here we can also use the colMeans function
mean_returns

# Variance (Sample)
monthly_var <- apply(monthly_returns, 2,var)
monthly_var

# Standard Deviation (Sample)
monthly_sigma <- apply(monthly_returns, 2,sd)
monthly_sigma

# Annual Statistics
mean_returns * 12
monthly_var * 12
monthly_sigma * sqrt(12)

# Calculating Covariance and Correlation
cov_mat<-cov(monthly_returns) # Covariance Matrix
cov_ab <- cov_mat["AAPL","K"] 
cov_ab # Covariance of Apple and K

# Correlation Coefficient Method 1 (Correlation Matrix)
cor(monthly_returns) # Correlation Matrix
corr_ab <- cor(monthly_returns)["AAPL","K"]
corr_ab # Correlation between Apple and K

## A Different View of the Correlation Coefficient
linear_model <- lm(monthly_returns[,"AAPL"] ~ monthly_returns[,"K"], 
                   data = monthly_returns)

# Scatter Plot
plot.ts(monthly_returns[,"AAPL"],monthly_returns[,"K"],  xy.labels = FALSE, 
        main = "Monthly Returns",
        ylab = "AAPL",
        xlab = "K")

# Add Trendline
abline(linear_model)

# Add Equation
equation <- paste(as.character(colnames(monthly_returns)[1]),
                  "=", round(linear_model$coefficients[1], digits = 3),
                  "+", round(linear_model$coefficients[2], digits = 3),
                  "x", as.character(colnames(monthly_returns)[2]),
                  sep = " ", collapse = NULL)
text(-0.1,0.1, equation)

# Check that R-sqr is equal to rho squared 
summary(linear_model)$r.squared
corr_ab ^ 2



###  Calculating Portfolio Means and Variances
# Input: proportions
Xa <- 0.5 # change if you want a different portfolio
proportions <- c(Xa, 1 - Xa) 

# Portfolio mean returns
port_returns <- sum(mean_returns * proportions)
port_returns

# Portfolio Return Variance
port_var <- sum(proportions * proportions * monthly_var,2 * prod(proportions) * cov_ab)

# alternative and more elegant approach to calculate portfolio variance
port_var <- proportions %*% cov_mat %*% proportions

# Portfolio standard deviation
port_sigma <- sqrt(port_var)

## Calculating the  Minimum-Variance Portfolio
# Analytical Approach
# This returns the weight of the 1st asset in the minimum variance portfolio
GMVP_AAPL <- as.numeric((monthly_var["K"] - cov_ab)/
  (monthly_var["AAPL"] + monthly_var["K"] - 2 * cov_ab))
GMVP_AAPL # AAPL's weight in the GMVP


# Calculating GMVP variance:
GMVP_prop <- c(XGMVP_AAPL = GMVP_AAPL, XGMVP_K = 1 - GMVP_AAPL)
GMVP_var <- GMVP_prop %*% cov_mat %*% GMVP_prop
GMVP_sigma<-sqrt(GMVP_var)
GMVP_var
GMVP_sigma


## Solver Approach
# This formula calculates the variance of a 2 asset portfolio
port_var_2Stocks <- function(X1, cov_mat) {
  prop<-c(X1,1-X1)
  prop %*% cov_mat %*% prop
}

# Minimizing the portfolio's variance:
GMVP.Result <- optimize(port_var_2Stocks, # the function to be optimized
                        interval = c(0,10), # the end-points of the interval to be searched for the minimum
                        tol = 0.00001, # the desired accuracy
                        maximum = FALSE, # FALSE FOR maximum, TRUE for minimum
                        cov_mat = cov_mat) # Constant Variables

GMVP.Result # We get the same variance result as the analytic approach
X1 <- GMVP.Result$minimum
X1


###   Portfolio Mean and Variance - Case of N Assets
## 4 assets example
# Input: Variance Covariance Matrix
var_cov_mat <- matrix(c(0.10,	0.01,	0.03,	0.05,
                        0.01,	0.30,	0.06,	-0.04,
                        0.03,	0.06,	0.40,	0.02,
                        0.05,	-0.04,	0.02,	0.50), 
                      ncol = 4)

# Adding column names to a matrix:
row.names(var_cov_mat) <- c("Asset 1",	"Asset 2", "Asset 3",	"Asset 4")
colnames(var_cov_mat) <-  row.names(var_cov_mat)

# Insert mean returns
mean_returns <- c(0.06, 0.08, 0.1, 0.15)

#Insert Portfolio weights
port_x_wei <- c(0.2, 0.3, 0.4, 0.1)
port_y_wei <- c(0.2, 0.1, 0.1, 0.6)

# Portfolios Mean Returns
port_x_ret <- mean_returns %*% port_x_wei
port_y_ret <- mean_returns %*% port_y_wei

# Portfolios variance calculation of n-assets
port_x_var <- port_x_wei %*% var_cov_mat %*% port_x_wei
port_y_var <- port_y_wei %*% var_cov_mat %*% port_y_wei

# Portfolios Standard deviation
port_x_sigma <- sqrt(port_x_var)
port_y_sigma <- sqrt(port_y_var)

# Covariance(X,Y)
cov_xy <- port_x_wei %*% var_cov_mat %*% port_y_wei

# Correlation Coefficient (X,Y)
correl_xy <- cov_xy / (port_x_sigma * port_y_sigma)

###   Envelope Portfolios

# Calculating returns of combinations of Portfolio X and Portfolio Y
# Z is a portfolio constructed of portfolio X and portfolio Y. 
X <- 0.3 #Proportion of X
port_z_weights <- c(X, 1-X) #Proportions of X and Y

# Portfolio Returns
port_z_ret <- port_z_weights %*% c(port_x_ret, port_y_ret)

# Portfolio Variance and sigma+
port_z_var <- port_z_weights[1] ^ 2 * port_x_var +
            port_z_weights[2] ^ 2 * port_y_var + 2 * prod(port_z_weights) * cov_xy

port_z_sigma <- sqrt(port_z_var)

## Plot
# Set an Array of Portfolio Weights
port_x_weight <- c( 1:20 / 10 - 0.6)
port_y_weight <- 1 - port_x_weight

# Portfolios Mean Returns
port_z_ret <- port_x_weight %*% port_x_ret + port_y_weight %*% port_y_ret

# Portfolios Sigma
port_z_sigma <- sqrt(port_x_weight ^ 2 %*% port_x_var + port_y_weight ^ 2 %*% 
                    port_y_var + (port_x_weight * port_y_weight) %*% cov_xy %*% 2)

# Plot
plot(data.frame(port_z_sigma, port_z_ret), type="l", xlab = "Standard Deviations",
     ylab = "Mean Returns", main = "Portfolio Means and Returns")
