# 
library(rjags)
library(tidyverse)
library(ggmcmc)

all_experiments <- read_csv("BHM/Data/TEMPEST_lab_experiments_all_L2.csv")

# Start with just calcium soluability in the ASW batch experiments
all_experiments_simplest <- all_experiments %>% 
  filter(Exp_Type %in% c("ASW", "NaCl")) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Control", "Freshwater", "ASW",
                                                     "0", "0.1", "1", "5", "25", "100"))) %>% 
 # group_by(Treatment) %>% 
  group_by(Exp_Type, Treatment) %>% 
  mutate(treatment_id = cur_group_id()) %>% 
  ungroup()

treatment_solutions <- all_experiments_simplest %>% 
  group_by(treatment_id) %>% 
  mutate(n = n()-1) %>% 
  filter(Timepoint == 0) %>% 
  select(Exp_Type, Treatment, treatment_id, 
         Ca_ppm, conductivity_uSpercm, n) %>% 
  mutate(is_nacl = ifelse(Exp_Type == "NaCl", 1, 0))

response_Ca <- all_experiments_simplest %>% 
  filter(Timepoint != 0) %>% 
  mutate(Ca_ppm = ifelse(Ca_ppm == 0, 0.01, Ca_ppm)) %>% 
  select(treatment_id, 
         Timepoint,
         Ca_ppm) %>% 
  pivot_wider(values_from = Ca_ppm, names_from = treatment_id) %>% 
  select(-Timepoint)

response_cond <- all_experiments_simplest %>% 
  filter(Timepoint != 0) %>% 
  select(treatment_id, 
         Timepoint,
         conductivity_uSpercm) %>% 
  pivot_wider(values_from = conductivity_uSpercm, names_from = treatment_id) %>% 
  select(-Timepoint)

response_doc <- all_experiments_simplest %>% 
  filter(Timepoint != 0) %>% 
  select(treatment_id, 
         Timepoint,
         doc_mgperL) %>% 
  pivot_wider(values_from = doc_mgperL, names_from = treatment_id) %>% 
  select(-Timepoint)

jags_data <- list(Ca_obs = as.matrix(response_Ca),
                  cond_obs = as.matrix(response_cond),
                  Ca_flood = treatment_solutions$Ca_ppm,
                  cond_flood = treatment_solutions$conductivity_uSpercm,
                  # treatment_id = treatment_solutions$treatment_id,
                  N_obs = treatment_solutions$n,
                  is_nacl = treatment_solutions$is_nacl,
                  
                  n_experiments = max(treatment_solutions$treatment_id)
                  
                  )

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

library(rjags)

inits_1 <- list(lambda_Ca = 0.5, lambda_cond = 0.4, Ca_inf_ASW = 0.5, 
                sigma_Ca = 0.2, sigma_cond = 0.1, Ca_b = -4, 
                Ca_a = 0.7, cond_a=0.5, cond_b = 1)
inits_2 <- list(lambda_Ca = 0.3, lambda_cond = 0.2, Ca_inf_ASW = 0.8, 
                sigma_Ca = 0.3, sigma_cond = 0.5, Ca_b = -3, 
                Ca_a = 0.5,cond_a=0.1, cond_b = 5)
inits_3 <- list(lambda_Ca = 0.7, lambda_cond = 0.1, Ca_inf_ASW = 0.3, 
                sigma_Ca = 0.1, sigma_cond = 0.2, Ca_b = -5, 
                Ca_a = 0.6,cond_a=1, cond_b = 1)
inits_4 <- list(lambda_Ca = 0.6, lambda_cond = 0.2, Ca_inf_ASW = 0.5, 
                sigma_Ca = 0.4, sigma_cond = 0.4, Ca_b = -4, 
                Ca_a = 0.8, cond_a=0.1, cond_b = 2)

jm <- jags.model(
  file     = textConnection(Ca_flush_model),
  data     = jags_data,
  inits    = list(inits_1, inits_2, inits_3, inits_4),
  n.chains = 4,
  n.adapt  = 2000
)

update(jm, n.iter = 5000)

params <- c("lambda_Ca", "lambda_cond", "Ca_inf_ASW", 
            "sigma_Ca", "Ca_a", "Ca_b", "cond_a", "cond_b", 
            "cond_init", "sigma_cond")

samples <- coda.samples(
  model          = jm,
  variable.names = params,
  n.iter         = 10000,
  thin           = 5
)

library(ggmcmc)

tidy_jags <- ggs(samples)

ggplot(tidy_jags, aes(x = value)) +
  geom_density(aes(fill = as.character(Chain)), alpha = 0.5) +
  facet_wrap(.~Parameter, scale = "free")

ggplot(tidy_jags, aes(x = Iteration, y = value)) +
  geom_line(aes(color = as.character(Chain))) +
  facet_wrap(.~Parameter, scale = "free")

subset_jags <- tidy_jags %>% 
  pivot_wider(names_from = Parameter, 
              values_from = value
  ) %>% 
  slice_sample(n = 100)

Ca_flood = treatment_solutions$Ca_ppm
cond_flood = treatment_solutions$conductivity_uSpercm
N_obs = 6
n_conc = max(treatment_solutions$treatment_id)
is_nacl = treatment_solutions$is_nacl
Ca_inf_NaCl <- 0.01

for (j in 1:nrow(subset_jags)) {
  
  lambda_Ca <- subset_jags$lambda_Ca[j]
  lambda_cond <- subset_jags$lambda_cond[j]
  
  # Initial Ca in porewater
  Ca_a <- subset_jags$Ca_a[j]
  Ca_b <- subset_jags$Ca_b[j]
  
  # asymptotes
  cond_a <- subset_jags$cond_a[j]
  cond_b <- subset_jags$cond_b[j]
  
  cond_init <- subset_jags$cond_init[j]
  
  # Ca asymptote for ASW treatments
  Ca_inf_ASW <- subset_jags$Ca_inf_ASW[j]
  
  Ca_mobile <- matrix(rep(NA, n_conc*N_obs), 
                      ncol = n_conc,
                      nrow = N_obs
                      )
  cond_mobile <- matrix(rep(NA, n_conc*N_obs), 
                      ncol = n_conc,
                      nrow = N_obs)
  
  Ca_init <- rep(NA, n_conc)
  Ca_inf <- rep(NA, n_conc)
  cond_inf <- rep(NA, n_conc)
  
  Ca_pw <- matrix(rep(NA, n_conc*N_obs), 
                      ncol = n_conc,
                      nrow = N_obs)
  cond_pw <- matrix(rep(NA, n_conc*N_obs), 
                    ncol = n_conc,
                    nrow = N_obs)
  
  for (c in 1:n_conc) {
    
    # Calcium initial condition
    Ca_init[c] <- exp(Ca_a * log(cond_flood[c]) +  Ca_b)
    
    Ca_inf[c] <- is_nacl[c]*Ca_inf_NaCl + (1-is_nacl[c]) * Ca_inf_ASW
    
    Ca_mobile[1, c] <- (Ca_flood[c] - Ca_inf[c]) + Ca_init[c]
    
    Ca_pw[1, c] <- Ca_mobile[1, c] + Ca_inf[c]
    
    # Cond initial cond
    # Conductivity flusing 
    cond_inf[c] <- exp(cond_a * log(cond_flood[c]) + cond_b)
    
    cond_mobile[1, c] <- (cond_flood[c] - cond_inf[c]) + cond_init
    
    cond_pw[1, c] <- cond_mobile[1, c] + cond_inf[c]
    
    for (i in 2:N_obs) {
      
      # Ca pool updates
      Ca_mobile[i, c] <- Ca_mobile[i-1, c] * lambda_Ca
      
      # Porewater Ca — mobile + exchange
      Ca_pw[i, c] <- Ca_mobile[i, c] + Ca_inf[c]
      
      # Conductivity pool pudates 
      cond_mobile[i, c] <- cond_mobile[i-1, c] * lambda_cond 
      
      cond_pw[i, c] <- cond_mobile[i, c] + cond_inf[c]

    } # end loop i, rows
    
  } # end loop c, concentrations
  
  colnames(Ca_pw) <- 1:n_conc
  colnames(cond_pw) <- 1:n_conc

  df_ca <- Ca_pw %>%
    as.data.frame() %>% 
    mutate(Timepoint = 1:6) %>% 
    pivot_longer(
      cols = -Timepoint,           # Exclude the 'Row' column from pivoting
      names_to = "treatment_id",   # Name the new column for column identifiers
      values_to = "value"    # Name the new column for values
    ) %>% 
    mutate(analyte = "Ca_ppm")
  
  df_cond <- cond_pw %>%
    as.data.frame() %>% 
    mutate(Timepoint = 1:6) %>% 
    pivot_longer(
      cols = -Timepoint,           # Exclude the 'Row' column from pivoting
      names_to = "treatment_id",   # Name the new column for column identifiers
      values_to = "value"    # Name the new column for values
    ) %>% 
    mutate(analyte = "conductivity_uSpercm")
  
  df_long <- df_ca %>% 
    bind_rows(df_cond) %>% 
    mutate(j = j)
  
  if (j == 1) {
    df_long_output <- df_long
  } else {
    df_long_output <- bind_rows(df_long_output, df_long)
  }
  
}

treatment_indexes <- treatment_solutions %>% 
  select(Exp_Type, Treatment, treatment_id)

df_long_plot <- df_long_output %>% 
  group_by(treatment_id, Timepoint, analyte) %>% 
  summarise(upper_ci = quantile(value, 0.975),
            lower_ci = quantile(value, 0.025),
            value = median(value)) %>% 
  mutate(treatment_id = as.integer(treatment_id)) %>% 
  left_join(treatment_indexes)

all_experiments_simplest_long <- all_experiments_simplest %>% 
  select(Exp_Type, Treatment, treatment_id, conductivity_uSpercm, Ca_ppm, Timepoint) %>% 
  pivot_longer(values_to = "value",
               names_to = "analyte",
               -c(Exp_Type, Treatment, treatment_id, Timepoint))

ggplot(df_long_plot, aes(x = Timepoint, y = value, color = Treatment)) +
  geom_point(data = all_experiments_simplest_long) +
  geom_line() +
  geom_ribbon(aes(ymin=lower_ci, ymax = upper_ci, fill = Treatment), alpha=0.3) +
  facet_grid(analyte~Exp_Type, scale = "free_y") +
  scale_y_log10()
  