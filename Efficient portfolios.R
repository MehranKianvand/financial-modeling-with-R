

###################################
#  EFFICIENT PORTFOLIOS
###################################

### Calculating the Efficient Frontier
## Inputs
# Variance-covariance matrix
var_cov_mat <- matrix(c(0.10,	0.01,	0.03,	0.05,
                        0.01,	0.30,	0.06,	-0.04,
                        0.03,	0.06,	0.40,	0.02,
                        0.05,	-0.04,0.02,	0.50),
                        nrow=4)

# Mean Returns
assets_mean_ret <- c(0.06,0.08,0.1,0.15)




## Find Efficient Portfolios
# Function: Calculate Envelope Portfolio proportions in one step
envelope_portfolio <- function(var_cov_mat, assets_mean_ret, constant){
  Z <- solve(var_cov_mat) %*% (assets_mean_ret - constant) # Z
  Z_Sum <- sum(Z) # Sum of Z's
  return( as.vector(Z/Z_Sum) ) # Efficient Portfolio proportions
}

# Computing Portfolios X and Y Asset Allocation
port_x_prop <- envelope_portfolio(var_cov_mat, assets_mean_ret, 0)
port_y_prop <- envelope_portfolio(var_cov_mat, assets_mean_ret, 0.04)

# Mean Returns
port_x_ret <- assets_mean_ret %*% port_x_prop
port_y_ret <- assets_mean_ret %*% port_y_prop

# Portfolios Variance and Covariance
port_x_var <- port_x_prop %*% var_cov_mat %*% port_x_prop
port_y_var <- port_y_prop %*% var_cov_mat %*% port_y_prop
cov_xy <- port_x_prop %*% var_cov_mat %*% port_y_prop

# Portfolios Standard Deviations
port_x_sigma <- sqrt(port_x_var)
port_y_sigma <- sqrt(port_y_var)

# Correlation Coefficient (Rho)
port_correl<-cov_xy / (port_x_sigma * port_y_sigma)


## Plotting the Envelope
# Setting Different Proportions of X
port_x_props <- c(c(1:25)/2 - 5)
port_y_props <- 1 - port_x_props

# Calculating Mean Returns
port_z_ret <- port_x_props * as.numeric(port_x_ret) + 
              port_y_props * as.numeric(port_y_ret)
port_z_sigma <- sqrt(port_x_props^2 * as.numeric(port_x_var) + 
                       port_y_props^2 * as.numeric(port_y_var) + 
                       (port_x_props * port_y_props) * as.numeric(cov_xy) * 2)

# Plot Envelope
plot(port_z_sigma, port_z_ret, main="The Envelope in Our Example", 
     type="l", xlim=c(0, 0.9), ylim=c(0, 0.2))

# Plot X and Y
lines(port_x_sigma, port_x_ret, type = "b", col = "red")
text(port_x_ret ~ port_x_sigma, labels = c("X"),cex = 0.9, font = 2, pos = 4)
lines(port_y_sigma, port_y_ret, type = "b", col = "blue")
text(port_y_ret ~ port_y_sigma, labels = c("Y"),cex = 0.9, font = 2, pos = 4)

#Assets Sigma
assets_sd <- sqrt(diag(var_cov_mat))

# Plot Assets
points(assets_sd, assets_mean_ret, col = "Green")
text(assets_mean_ret ~ assets_sd,
     labels = c("Asset 1", "Asset 2", "Asset 3", "Asset 4"),
     cex=0.9, font=2, pos=4)


### Computing the Global Minimum Variance Portfolio (GMVP)
## Computing the Global Minimum Variance Portfolio (GMVP)
# FUNCTION: GMVP as row
GMVP_as_row <- function(var_cov_mat){
  colSums( solve(var_cov_mat) ) / 
    sum( solve(var_cov_mat) )
}

GMVP_prop <- GMVP_as_row(var_cov_mat)

# GMVP Mean Return
GMVP_mean_ret <- GMVP_prop %*% assets_mean_ret

# GMVP Sigma
GMVP_Sigma <- sqrt(GMVP_prop %*% var_cov_mat %*% GMVP_prop)

# Add GMVP to the envelope
lines(GMVP_Sigma, GMVP_mean_ret, type = "b", col = "Orange")
text(GMVP_mean_ret ~ GMVP_Sigma, labels = c("GMVP"),cex = 0.9, font = 2, pos = 2)

###Efficient Portfolios Without Short Sales
# Input: Assets Mean Returns
assets_mean_ret <- c(0.02, 0.02, 0.08, 0.1)

# In this part, we set the constraints for the matrix as Ax => b (linear constraints)
# the constraints are: X1 => 0; X2 => 0; X3 => 0; X4 => 0 

# Setting Coefficients Matrix
A_matrix <- diag(1,nrow = 4) # identity matrix to enforce all asset proportions are > 0
A_matrix <- rbind(A_matrix, rep(1, 4)) # sum of proportions => 1
A_matrix <- rbind(A_matrix, rep(-1, 4)) # sum of proportions <= 1

# Setting Result vector
B_vector <- rep(0, 4)
B_vector <- c(B_vector, 1, -1)

# Setting Target Function: Sharpe
sharpe <- function(proportions, var_cov_mat, assets_mean_ret, constant){
  ex_return = proportions %*% (assets_mean_ret - constant)
  sigma = sqrt(proportions %*% var_cov_mat %*% proportions)
  return(ex_return / sigma)}

# Run (Note: Initial values must be in interior of the feasible region) 
result <- constrOptim(theta = c(0.25, 0.25, 0.25, 0.25), #  initial value                      
                      f = sharpe, # function to maximize
                      grad = NULL, # gradient of f
                      ui = A_matrix, # constraint matrix
                      ci = B_vector- 1e-05, # constraint vector minus infitisimal
                      mu = 1e-05, # result accuracy
                      control = list(fnscale = -1), # maximization 
                      outer.iterations = 100000,
                        method = "Nelder-Mead", # use this whithout gradient function
                      # target function variables:
                      var_cov_mat = var_cov_mat,
                      assets_mean_ret = assets_mean_ret,
                      constant = 0.02)
# Result
round(result$par, 4)
sharpe(result$par, var_cov_mat, assets_mean_ret, 0.02)



## No short (solver 2nd approach)
# Adding a rP = 9% Constraint
A_matrix <- rbind(A_matrix, assets_mean_ret) # Add a return = 9% costraint
B_vector <- c(B_vector, 0.09)

# Setting Target Function: Variance
port_var <- function(proportions, var_cov_mat){
  return(proportions %*% var_cov_mat %*% proportions)
}

# Run
result <- constrOptim(theta = c(0, 0, 0.25, 0.75), #  starting value 
                      f = port_var, # function to maximize
                      grad = NULL, # gradient of f
                      ui = A_matrix, # constraint matrix
                      ci = B_vector - 1e-05, # constraint vector 
                      mu = 1e-03, # result accuracy
                      outer.iterations = 1000000,
                      method = "Nelder-Mead", # use this whithout gradient function
                      # target function variables:
                      var_cov_mat = var_cov_mat)

# Result
round(result$par, 4) # Proportions
assets_mean_ret %*% result$par # Mean Return
sqrt(port_var(result$par, var_cov_mat)) # Sigma




### SOVER ALTERNATIVES
## Solving constrained Minimum Variance problem with "quadprog" package
#install.packages("quadprog")
library(quadprog)

# Set the constraints for the matrix as Ax => b (linear constraints)
A_matrix <- cbind(rep(1,4),# Sum of weights = 1 (meq constraint - see below)
                  diag(4), # each asset weight is => 0
                  assets_mean_ret) # return is => 9%
b_vector <- c(1, rep(0, 4), 0.09)

min_var <- solve.QP(Dmat = var_cov_mat,
                    dvec = matrix(assets_mean_ret, ncol = 1), # returns should be in a matrix
                    Amat = A_matrix, 
                    bvec = b_vector, 
                    meq = 1) # the first meq constraints are treated as equality constraints, all further as inequality constraints (defaults to 0).
min_var$solution #weights
assets_mean_ret %*% min_var$solution # Mean Return
sqrt(port_var(min_var$solution, var_cov_mat)) # Sigma

## Solving Maximum Sharpe Ratio
#install.packages("nloptr")
library(nloptr)

# Constraint Function
constr_fun <- function(proportions, var_cov_mat, assets_mean_ret, constant){
  return(c(sum(proportions) - 1, 
           -sum(proportions) + 1))} # Enforce sum of proportions = 100%

res <- nloptr(x0 = c(0.25, 0.25, 0.25, 0.25), # Initial Guess
              eval_f = function(proportions, var_cov_mat, assets_mean_ret, constant){
                -sharpe(proportions, var_cov_mat, assets_mean_ret, constant)}, # Function to MINIMIZE (-sharpe)
              eval_grad_f = NULL, 
              lb = rep(0,4), 
              ub = rep(1,4),
              eval_g_ineq = constr_fun,
              eval_jac_g_ineq = NULL, 
              eval_g_eq = NULL,
              eval_jac_g_eq = NULL, 
              opts = list("algorithm" = "NLOPT_LN_COBYLA",
                          "xtol_rel"=1.0e-5,
                          "maxeval"  = 10000),
              
              # Sharpe inputs:
              var_cov_mat = var_cov_mat, 
              assets_mean_ret = assets_mean_ret, 
              constant = 0.02
)

res$solution
sum(res$solution)
sharpe(res$solution, var_cov_mat, assets_mean_ret, 0.02) 

