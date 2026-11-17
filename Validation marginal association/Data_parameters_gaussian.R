# number of simulations
M <- 1000

# parameters for the marginal association coefficient estimation
time <- 1:20
R <- 10^5
n_rep <- 1
a <- 9

# maximum number of repetitive observations for each patient
n <- 12

# probability to observed each variable
p1 <- 1 - 0.28 # probability to observed x 
p2 <- 1 - 0.45 # probability to observed y

# fixed effects coefficients
true_betas <- c(10.086, 
                -0.128, 
                53.240, 
                -1.974) # beta coefficients

# variances of the error terms
sd_x <- 0.32
sd_y <- 3.58

# variance-covariance matrices of the random effects
D_x <- matrix(c(0.14, -0.013, -0.013, 0.017), 2, byrow=TRUE) # var-cov matrix of the random effects of x
D_y <- matrix(c(0.12, -0.006, -0.006, 0.34), 2, byrow=TRUE)  # var-cov matrix of the random effects of y


# scaling factor of the JSM
gamma <- 1.30 


# marginal coefficient over time
coeff_time_jsm <- function(time, D_x, sd_x, gamma){
  var_m <- lapply(time, function(x) D_x[1, 1] + 2*x*D_x[1, 2] + (x^2)*D_x[2, 2])
  var_x <- lapply(time, function(x) D_x[1, 1] + 2*x*D_x[1, 2] + (x^2)*D_x[2, 2] + sd_x^2)
  coeff <- mapply(function(x, y) gamma*(x/y), var_m, var_x)
  
  return(coeff)
}
real_coeff <- coeff_time_jsm(time, D_x, sd_x, gamma)