rm(list =ls())
set.seed(123)

library(MASS)

n <- 300 # set n

# Desired correlation matrix

fu_ab <- -.5 # set to 0 for type I error, set to any other value for power analysis

cor_mat <- matrix(c(
  1.00,  0.2, -0.30,  0.20,  # eta
  0.2,  1.00,  0.30,  -0.30,  # ses
  -0.30,  0.30,  1.00,  fu_ab,  # abilità
  0.20,  -0.30,  fu_ab,  1.00   # logit fumo
), nrow = 4, byrow = TRUE)


# Variable names
colnames(cor_mat) <- rownames(cor_mat) <- c("eta", "ses", "abilita", "fumo")
cor_mat

# expected beta

df <- data.frame(mvrnorm(n, mu = rep(0, 4), Sigma = cor_mat, empirical = T))
# Trasform logit into probabilities
df$p_fumo   <- exp(df$fumo) / (1 + exp(df$fumo))
df$fumatore <- ifelse(df$p_fumo > .5, 1, 0)

cor_mat2 <- cor(df[, c('eta', 'ses', 'abilita', 'fumatore')])
cor_mat2 |> round(2)

P <- solve(cor_mat2)  # inversione della matrice
beta_yx_zw <- -P[4, 3] / P[4, 4]
beta_yx_zw |> round(2)

# Simulation
iter <- 1000

pvals_zero <- numeric(iter)
pvals_adj <- numeric(iter)
pvals_weights <- numeric(iter)
pvals_comb <- numeric(iter)
beta_zero <- numeric(iter)
beta_adj <- numeric(iter)
beta_weights <- numeric(iter)
beta_comb <- numeric(iter)

df_simul <- data.frame(pvals_zero, pvals_adj, pvals_weights, pvals_comb, beta_zero, beta_adj, beta_weights, beta_comb)

for(i in 1:iter){
  df <- data.frame(mvrnorm(n, mu = rep(0, 4), Sigma = cor_mat, empirical = T))
  # Trasform logit into probabilities
  df$p_fumo   <- exp(df$fumo) / (1 + exp(df$fumo))
  df$fumatore <- rbinom(n, size = 1, prob = df$p_fumo)
  
  # adjusted models
  m1 <- lm(abilita ~ fumatore, data = df) #zero model
  m2 <- update(m1, . ~ . + ses + eta, data = df) # adjusted only
  
  # compute propensity scores
  m3 <- glm(fumatore ~ ses + eta, data = df, family = binomial()) # logistic regression to compute ps
  ps <- exp(predict(m3))/(1 + exp(predict(m3))) # PSs
  iptw <- ifelse(df$fumatore == "sì", 1/ps, 1/(1-ps)) #IPTWs
  
  # PS wrghted models
  m4 <- lm(abilita ~ fumatore, weights =  iptw, data = df) # PS model
  m5 <- update(m4,.  ~ . + ses + eta, data = df) # combined (PS and adjust) model
  
  # store the results
  df_simul$pvals_zero[i] <- summary(m1)$coefficients['fumatore', 'Pr(>|t|)']
  df_simul$pvals_adj[i] <- summary(m2)$coefficients['fumatore', 'Pr(>|t|)']
  df_simul$pvals_weights[i] <- summary(m4)$coefficients['fumatore', 'Pr(>|t|)']
  df_simul$pvals_comb[i] <- summary(m5)$coefficients['fumatore', 'Pr(>|t|)']
  df_simul$beta_zero[i] <- summary(m1)$coefficients['fumatore', 'Estimate']
  df_simul$beta_adj[i] <- summary(m2)$coefficients['fumatore', 'Estimate']
  df_simul$beta_weights[i] <- summary(m4)$coefficients['fumatore', 'Estimate']
  df_simul$beta_comb[i] <- summary(m5)$coefficients['fumatore', 'Estimate']
  
}

sum(df_simul$pvals_zero < .05)/iter
sum(df_simul$pvals_adj < .05)/iter
sum(df_simul$pvals_weights < .05)/iter
sum(df_simul$pvals_comb < .05)/iter # più potente

summary(df_simul[, 5:ncol(df_simul)])

# confronta con valori attesi
cor_mat2['abilita','fumatore' ] |> round(2)
beta_yx_zw |> round(2)
