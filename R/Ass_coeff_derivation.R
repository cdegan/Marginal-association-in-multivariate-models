##################### FUNCTIONS TO DERIVE THE ASSOCIATION COEFFICIENT IN JSM #######################################

#------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT IN JSM (first formulation) WITH THE CLOSED FORMULA (NORMAL CASE) ------------------- 

# INPUT:
# model: JSM model
# time: time points for which I want the estimation (number/vector)

coeff_jsm_old <- function(model, time) {
  # resampling
  temp <- inla.hyperpar.sample(10^5, model)
  
  var_x <- lapply(time, function(x) 1/temp[, 1] + 1/temp[, 5] + (x^2 * 1/temp[, 6]) + (2*x*temp[, 7]*sqrt(1/temp[, 5])*sqrt(1/temp[,6])))
  var_m <- lapply(time, function(x) 1/temp[, 5] + (x^2 * 1/temp[, 6]) + (2*x*temp[, 7]*sqrt(1/temp[, 5])*sqrt(1/temp[,6])))
  gamma <- temp[, 11]
  
  # closed formula
  coeff_sampled <- mapply(function(y, z) gamma*(y/z), var_m, var_x)
  coeff <- cbind(time, t(apply(coeff_sampled, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))))
  
  return(coeff)
}

# OUTPUT:
# coeff: matrix of association coefficients (and its credible intervals) estimated for each time point



#-------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT WITH THE INTEGRATION IN JSM (first formulation) ----------------------

# INPUT:
# model: JSM model
# fix_eff: fixed_effects of JSM (vector)
# R: number of Montecarlo samples (number)
# n: number of hyperparameters resamples (number)
# a: value of x that we are interesting on (number)
# t: for which time points I want the estimations (number/vector)
# cores: number of cores available (number)
# gamma: gamma parameter of JSM (number)
# link_function: name of the link function used in the model (character)

marginal_jsm_old <- function(model, fix_eff, R, n, a, t, cores, gamma, link_function){
  temp <- inla.posterior.sample(n, model, selection = list(Int.x = 1, t.x = 1, Int.y = 1, t.y = 1)
  
  # marginal of the random effects b
  mean_b <- rep(0, 4)
  var_b <- lapply(temp, function(x){
    cov1 <- matrix(c(1/x$hyperpar[5], x$hyperpar[7]*sqrt(1/x$hyperpar[5]*1/x$hyperpar[6]), x$hyperpar[7]*sqrt(1/x$hyperpar[5]*1/x$hyperpar[6]), 1/x$hyperpar[6]), 2, 2)
    cov2 <- matrix(c(1/x$hyperpar[8], x$hyperpar[10]*sqrt(1/x$hyperpar[8]*1/x$hyperpar[9]), x$hyperpar[10]*sqrt(1/x$hyperpar[8]*1/x$hyperpar[9]), 1/x$hyperpar[9]), 2, 2)
    as.matrix(bdiag(cov1, cov2))
  })

  # sample-specific gamma
  gamma_temp <- sapply(temp, function(x) x$hyperpar[11])
  
  # marginal of x
  mean_x <- lapply(temp, \(x) as.vector(c(1, t) %*% x$latent[1:2, ]))
  sd_x <- unlist(lapply(temp, \(x) sqrt(1/x$hyperpar[5] + t^2/x$hyperpar[6] + 2*t*x$hyperpar[7] * sqrt(1/x$hyperpar[5]) * sqrt(1/x$hyperpar[6]) + 1/x$hyperpar[1])))
  
  # posterior distribution of the random effects b|x
  vector_bx <- lapply(var_b, function(x) x[1, ] + t*x[, 2])
  mean_bx <- var_bx <- r_bx <- mean_bx2 <- r_bx2 <- list(NA)
  for(i in 1:n){
    mean_bx[[i]] <- mean_b + vector_bx[[i]]*sd_x[i]^(-2)*as.vector(a-mean_x[[i]])
    var_bx[[i]] <- as.matrix(round(var_b[[i]] - (vector_bx[[i]] * sd_x[i]^(-2)) %*% t(vector_bx[[i]]), 6))
    r_bx[[i]] <- mvtnorm::rmvnorm(R, mean_bx[[i]], var_bx[[i]], method = "svd")
    mean_bx2[[i]] <- mean_b + vector_bx[[i]]*sd_x[i]^(-2)*as.vector((a+1)-mean_x[[i]])
    r_bx2[[i]] <- mvtnorm::rmvnorm(R, mean_bx2[[i]], var_bx[[i]], method = "svd")
  }
  
  # conditional mean of y|b 
  beta_y <- lapply(temp, function(x) as.vector(c(1, t) %*% x$latent[3:4, ]))
  nu <- mclapply(seq_len(n), function(i) {apply(r_bx[[i]], 1, function(x) beta_y[[i]] + x[3] + t * x[4] + gamma_temp[i] * (mean_x[[i]] + x[1] + t * x[2]))}, mc.cores = cores)
  nu2 <- mclapply(seq_len(n), function(i) {apply(r_bx2[[i]], 1, function(x) beta_y[[i]] + x[3] + t * x[4] + gamma_temp[i] * (mean_x[[i]] + x[1] + t * x[2]))}, mc.cores = cores)

  mu_cond_y <- mclapply(nu, link_function, mc.cores = cores)
  mu_cond_y2 <- mclapply(nu2, link_function, mc.cores = cores)
  
  # averaging
  average <- unlist(mclapply(mu_cond_y, mean, na.rm = TRUE, mc.cores = cores))
  average2 <- unlist(mclapply(mu_cond_y2, mean, na.rm = TRUE, mc.cores = cores))
  mu <- mean(average2 - average)
  result <- mean(mu)
  var_tot <- var(average2 - average)
  
  return(c(result = result, variance = var_tot))
}

# OUTPUT:
# list with first element the association coefficient estimations, and second element the variance


#------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT IN JSM (second formulation) WITH THE CLOSED FORMULA (NORMAL CASE) ------------------- 

# INPUT:
# model: JSM model
# time: time points for which I want the estimation (number/vector)

coeff_jsm_new <- function(model, time) {
  # resampling
  temp <- inla.hyperpar.sample(10^5, model)
  
  var_x <- lapply(time, function(x) 1/temp[, 1] + 1/temp[, 4] + (x^2 * 1/temp[, 5]) + (2*x*temp[, 6]*sqrt(1/temp[, 4])*sqrt(1/temp[, 5])))
  var_m <- lapply(time, function(x) 1/temp[, 4] + (x^2 * 1/temp[, 5]) + (2*x*temp[, 6]*sqrt(1/temp[, 4])*sqrt(1/temp[, 5])))
  gamma <- temp[, 10]
  
  # closed formula
  coeff_sampled <- mapply(function(y, z) gamma*(y/z), var_m, var_x)
  coeff <- cbind(time, t(apply(coeff_sampled, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))))
  
  return(coeff)
}

# OUTPUT:
# coeff: matrix of association coefficients (and its credible intervals) estimated for each time point



#-------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT WITH THE INTEGRATION IN JSM (second formulation) ----------------------

# INPUT:
# model: JSM model
# R: number of Montecarlo samples (number)
# n: number of hyperparameters resamples (number)
# a: value of x that we are interesting on (number)
# time: for which time points I want the estimations (number/vector)
# cores: number of cores available (number)
# link_function: name of the link function used in the model (character)

marginal_jsm_new <- function(model, R, n, a, time, cores, link_function){
  # sampling from the joint posterior distribution of parameters and hyperparameters
  fixed <- inla.posterior.sample(n, model, selection = list(Int.x = 1, t.x = 1, Int.y = 1, t.y = 1))
  
  # marginal of the random effects b
  mean_b <- rep(0, 4)
  var_b <- lapply(fixed, function(x){
    cov1 <- matrix(c(1/x$hyperpar[4], x$hyperpar[6]*sqrt(1/x$hyperpar[4]*1/x$hyperpar[5]), x$hyperpar[6]*sqrt(1/x$hyperpar[4]*1/x$hyperpar[5]), 1/x$hyperpar[5]), 2, 2)
    cov2 <- matrix(c(1/x$hyperpar[7], x$hyperpar[9]*sqrt(1/x$hyperpar[7]*1/x$hyperpar[8]), x$hyperpar[9]*sqrt(1/x$hyperpar[7]*1/x$hyperpar[8]), 1/x$hyperpar[8]), 2, 2)
    as.matrix(bdiag(cov1, cov2))
  })
  
  # gamma resampling
  gamma_temp <- lapply(fixed, \(x) x$hyperpar[10])
  
  # calculation of the integral
  results <- lapply(time, function(t){
    # marginal of x
    mean_x <- lapply(fixed, \(x) as.vector(c(1, t) %*% x$latent[1:2, ]))
    sd_x <- unlist(lapply(fixed, \(x) sqrt(1/x$hyperpar[4] + t^2/x$hyperpar[5] + 2*t*x$hyperpar[6] * sqrt(1/x$hyperpar[4]) * sqrt(1/x$hyperpar[5]) + 1/x$hyperpar[1])))
    
    # posterior distribution of the random effects b|x
    vector_bx <- lapply(var_b, \(x) x[1, ] + t * x[, 2])
    mean_bx <- var_bx <- r_bx <- mean_bx2 <- r_bx2 <- vector("list", n)
    for(i in seq_len(n)){
      mean_bx[[i]] <- mean_b + vector_bx[[i]] * sd_x[i]^(-2) * as.vector(a - mean_x[[i]])
      var_bx[[i]] <- round(var_b[[i]] - (vector_bx[[i]] * sd_x[i]^(-2)) %*% t(vector_bx[[i]]), 6)
      r_bx[[i]] <- mvtnorm::rmvnorm(R, mean_bx[[i]], var_bx[[i]], method = "svd")
      mean_bx2[[i]] <- mean_b + vector_bx[[i]] * sd_x[i]^(-2) * as.vector((a + 1) - mean_x[[i]])
      r_bx2[[i]] <- mvtnorm::rmvnorm(R, mean_bx2[[i]], var_bx[[i]], method = "svd")
    }
    
    # conditional mean of y|b
    beta_y <- lapply(fixed, \(x) c(1, t) %*% x$latent[3:4, ])
    nu <- mclapply(seq_along(r_bx), function(i) apply(r_bx[[i]], 1, function(x) beta_y[[i]] + x[3] + t * x[4] + gamma_temp[[i]] * (mean_x[[i]] + x[1] + t * x[2])), mc.cores = cores)
    nu2 <- mclapply(seq_along(r_bx2), function(i) apply(r_bx2[[i]], 1, function(x) beta_y[[i]] + x[3] + t * x[4] + gamma_temp[[i]] * (mean_x[[i]] + x[1] + t * x[2])), mc.cores = cores)
    mu_cond_y  <- mclapply(nu,  link_function, mc.cores = cores)
    mu_cond_y2 <- mclapply(nu2, link_function, mc.cores = cores)
    
    # averaging
    average  <- unlist(mclapply(mu_cond_y, mean, mc.cores = cores))
    average2 <- unlist(mclapply(mu_cond_y2, mean, mc.cores = cores))
    delta <- average2 - average
    
    q <- quantile(delta, c(.025, .975), na.rm = TRUE)
    c(time = t, result = mean(delta), variance = var(delta), q025 = q[1], q975 = q[2])
  })
  
  do.call(rbind, results)
}

# OUTPUT:
# table with nrow = length(time) and as columns time, result, variance, q025, q975
