
Ca_flush_model <- "model {
  
  # PRIORS
  # flush rates
  lambda_Ca  ~ dbeta(2, 2)
  lambda_cond ~ dbeta(2, 2)

  # Initial Ca in porewater
  Ca_a ~ dnorm(0, 0.001) T(0,)
  Ca_b ~ dnorm(0, 0.001)
  
  cond_init ~ dnorm(0, 0.001) T(0,) 
  
  # Equation for conductivity aymptote
  cond_a ~ dnorm(0, 0.001) T(0,)
  cond_b ~ dnorm(0, 0.001)
 
  # Ca asymptote for ASW treatments
  Ca_inf_ASW ~ dnorm(0, 0.001) T(0,)  
  Ca_inf_NaCl <- 0.01
  
  # Observation noise
  sigma_Ca ~ dexp(5)
  sigma_cond ~ dexp(5)
  
  # Carbon model
  C_max ~ dnorm(0, 0.001) T(0,)
  k_m ~ dnorm(0, 0.001) T(0,)
  C_init ~ dnorm(0, 0.001) T(0,)
  
  C_c ~ dnorm(0, 1)
  C_b ~ dnorm(0, 1) T(0,)
  C_a ~ dnorm(0, 1) T(0,)
  
  for (c in 1:n_experiments) {
    
    # Calcium flushing
    Ca_init[c] <- exp(Ca_a * log(cond_flood[c]) +  Ca_b)
    
    Ca_inf[c] <- is_nacl[c]*Ca_inf_NaCl + (1-is_nacl[c]) *Ca_inf_ASW
    
    Ca_mobile[1, c] <- (Ca_flood[c] - Ca_inf[c]) + Ca_init[c]
    
    Ca_pw[1, c] <- Ca_mobile[1, c] + Ca_inf[c]
    
    # Log-normal likelihood
    Ca_obs[1, c] ~ dlnorm(log(Ca_pw[1, c]), 1 / sigma_Ca^2)
    
  
    # Conductivity flusing 
    cond_inf[c] <- exp(cond_a * log(cond_flood[c]) + cond_b)
    
    cond_mobile[1, c] <- (cond_flood[c] - cond_inf[c]) + cond_init
    
    cond_pw[1, c] <- cond_mobile[1, c] + cond_inf[c]
    
    # Log-normal likelihood
    cond_obs[1, c] ~ dlnorm(log(cond_pw[1, c]), 1 / sigma_cond^2)
  
    # Carbon module 
    # Carbon pool is initial carbon + carbon released from flooding treatment
    C_pool[c] <- C_init + (C_max * cond_flood[c] / (k_m + cond_flood[c]))
    
    lambda_C_logit[1, c] <-  C_c + C_b * Ca_pw[1, c] + C_a * cond_pw[1, c]  
    lambda_C[1, c] <- 1/(1+exp(-R_logit))
  
    for (i in 2:N_obs[c]) {

      # Ca pool updates
      Ca_mobile[i, c] <- Ca_mobile[i-1, c] * lambda_Ca

      # Porewater Ca — mobile + exchange
      Ca_pw[i, c] <- Ca_mobile[i, c] + Ca_inf[c]
  
      # Log-normal likelihood
      Ca_obs[i, c] ~ dlnorm(log(Ca_pw[i, c]), 1 / sigma_Ca^2)
      
      # Conductivity pool pudates 
      cond_mobile[i, c] <- cond_mobile[i-1, c] * lambda_cond 
    
      cond_pw[i, c] <- cond_mobile[i, c] + cond_inf[c]
    
      # Log-normal likelihood
      cond_obs[i, c] ~ dlnorm(log(cond_pw[i, c]), 1 / sigma_cond^2)
  
    } # end loop i, rows
    
  } # end loop c, concentrations
  
} # end model
"