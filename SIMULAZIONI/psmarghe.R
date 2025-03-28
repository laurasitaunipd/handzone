
rm(list = ls())

library(future.apply)
library(tidyverse)

# Set up parallel backend
plan(multisession, workers = parallel::detectCores() - 2)


# Run simulations
set.seed(1234)


simulate_100 = function(i, n, eff_fum_abi, eff_ses_abi, eff_ses_fum ) {
  
  results = tibble()  # Inizializza dataframe vuoto
  
  for (j in 1:100) {  # Ogni worker fa 100 simulazioni
    
    n = n
    ses = rnorm(n, 0, 1)
    
    # Define the effects
    eff_fum_abi = eff_fum_abi
    eff_ses_abi = eff_ses_abi
    eff_ses_fum = eff_ses_fum
    
    
    # Transform logit into probabilities
    logit_fumo = 0 + ses * eff_ses_fum + rnorm(n,0,1)
    
    # Calculate the probability of smoking
    p_fumo = 1 / (1 + exp(-logit_fumo)) 
    fumatore = rbinom(n, size = 1, prob = p_fumo)
    
    # Calculate the probability of ability
    logit_abi = 0 + fumatore * eff_fum_abi +  ses * eff_ses_abi + rnorm(n,0,1)
    p_abi = 1 / (1 + exp(-logit_abi))
    abilita = rbinom(n, size = 1, prob = p_abi)
    
    # Create the dataframe
    df = data.frame(
      ses = ses,
      fumatore = fumatore,
      abilita = abilita
    )
    
    # Fit the propensity score model
    mprop = glm(fumatore ~ ses, data = df, family = binomial())
    
    # Calculate the propensity scores
    df$ps = predict(mprop, type = "response")
    
    df$iptw=ifelse(df$fumatore == 1, 1/df$ps, 1/(1-df$ps))
    
    
    # Fit the model for 'abilita' without any adjust
    m1 = glm(abilita ~ fumatore, data = df, family = binomial())
    m1_pval = summary(m1)$coefficients[2,4]
    
    # Fit the model for 'abilita' with covariates 'ses'
    m2 = update(m1, . ~ . + ses, data = df)
    m2_pval = summary(m2)$coefficients[2,4]
    
    # Errore!!
    # Fit the model for 'abilita' with IPTW weights -  
    m3_binom = glm(abilita ~ fumatore, weights = iptw, data = df, family = binomial())
    m3_binom_pval = summary(m3_binom)$coefficients[2,4]
    
    # Fit the model for 'abilita' with IPTW weights and additional covariates
    m4_binom = update(m3_binom, . ~ . + ses, data = df)
    m4_binom_pval = summary(m4_binom)$coefficients[2,4]
    
    # Errore!!
    # Fit the model for 'abilita' with IPTW weights -  
    m3_quasibinom = glm(abilita ~ fumatore, weights = iptw, data = df, family = quasibinomial())
    m3_quasibinom_pval = summary(m3_quasibinom)$coefficients[2,4]
    
    # Fit the model for 'abilita' with IPTW weights and additional covariates
    m4_quasibinom = update(m3_quasibinom, . ~ . + ses, data = df)
    m4_quasibinom_pval = summary(m4_quasibinom)$coefficients[2,4]
    
    
    # Store results
    results <- bind_rows(results,tibble(m1_pval, m2_pval, 
                                        m3_binom_pval, m4_binom_pval,
                                        m3_quasibinom_pval, m4_quasibinom_pval))
  }
  return(results)
}


# Set values
n = 100
eff_fum_abi = 0.5
eff_ses_abi = 0
eff_ses_fum = 0

# Run simulations in parallel
n_cores = parallel::detectCores() - 2

star = Sys.time()
results = future_lapply(1:n_cores, function(i) simulate_100(i, n, 
                                                            eff_fum_abi = eff_fum_abi, 
                                                            eff_ses_abi = eff_ses_abi, 
                                                            eff_ses_fum = eff_ses_fum), 
                        future.seed = TRUE)
stp = Sys.time()
stp-star

# Stop parallel backend
plan(sequential)

# Results
results_df=bind_rows(results)

apply(results_df, MARGIN = 2, FUN = function(x) mean(x < 0.05))
# (quante volte il pvalue è significativo)/(numero tot di simulazioni)

mean(results_df$m1_pval < .05)


