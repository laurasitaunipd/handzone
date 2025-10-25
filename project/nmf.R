# installare
#install.packages("remotes")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("Biobase")
#remotes::install_github("jennalandy/causalLFO")
#install.packages("lsa")
#install.packages("RcppHungarian")


# pacchetti
library(NMF)
library(causalLFO)
library(tidyverse)
library(ggridges)
library(ggplot2)
# in caso aggiungi anche le ultime due packages come library se dà errore in B


############################# APPLICAZIONE SU 1 DATASET SIMULATO CON N=100
# A) SIMULAZIONE DEL DATASET
rm(list=ls())
set.seed(321)
N = 100; D = 96; K = 4; ATE = c(1000, 0, 0, 0)

# Simulate treatment assignment
Tr = sample(c(0, 1), N, replace = TRUE)

# Simulate latent factors P
true_P = matrix(rexp(D*K, rate = 1), nrow = D)
# Normalize factors to sum to 1
true_P = sweep(true_P, 2, colSums(true_P), '/')

# Simulate untreated factor loadings C
true_C = matrix(nrow = K, ncol = N)
true_C[1,] <- rgamma(N, shape = 1, scale = 1000) # larger scale for factor 1
true_C[2,] <- rexp(N, rate = 0.01)
true_C[3,] <- rexp(N, rate = 0.01)
true_C[4,] <- rexp(N, rate = 0.01)
true_C[4,sample(1:N, 10)] <- rnorm(10, mean = 1500, sd = 1000) # outliers for factor 4
data.frame(t(true_C)) %>%
  setNames(paste0("k", 1:K)) %>%   # rinomina le colonne come k1, k2, k3, k4
  pivot_longer(1:K, names_to = 'k', values_to = 'C') %>%
  ggplot(aes(x = C, y = as.factor(k))) +
  geom_density_ridges() +
  theme_bw() +
  labs(x = "Untreated latent outcome distribution", y = "Latent dimension")

# Add ATE to loadings of treated samples
for (k in 1:K) {
  true_C[k, Tr == 1] <- true_C[k, Tr == 1] + ATE[k]
}

# Comparing and visualizing true_C treated vs untreated
df <- data.frame(t(true_C)) %>% 
  setNames(paste0("k", 1:K)) %>%   # rinomina le colonne come k1, k2, k3, k4
  mutate(Tr = factor(Tr, levels = c(0,1), labels = c("Untreated", "Treated"))) %>% 
  pivot_longer(cols = paste0("k", 1:K), names_to = "k", values_to = "C")

grafico_datasim <- ggplot(df, aes(x = C, y = as.factor(k), fill = Tr)) +
  geom_density_ridges(alpha = 0.6) +
  theme_bw() +
  scale_y_discrete(labels = paste0("L", 1:K)) + 
  labs(
    x = "Individual contributions", 
    y = "Latent factors",
    fill = "Group"
  )

grafico_datasim

# Simulate M ~ Poisson(PC)
M = matrix(nrow = D, ncol = N)
for (i in 1:N) {
  M[,i] <- rpois(D, lambda = true_P %*% true_C[,i])
}

# B) STIME DI TRUE_P(LAMBDA) E TRUE_C(L)

all_data_fit <- all_data(M, Tr, rank = 4, reference_P = true_P)
impute_and_stabilize_fit <- impute_and_stabilize(M, Tr, rank = 4, reference_P = true_P)

reorder_sim_mat <- function(sim_mat) {
  best_match <- apply(sim_mat, 1, which.max)
  sim_reordered <- sim_mat[, best_match, drop = FALSE]
  rownames(sim_reordered) <- paste0("True L", 1:nrow(sim_reordered))
  colnames(sim_reordered) <- paste0("Estimated L", 1:ncol(sim_reordered))
  sim_reordered
}

sim_all <- reorder_sim_mat(all_data_fit$sim_mat)
sim_is  <- reorder_sim_mat(impute_and_stabilize_fit$sim_mat)

# visualizzo Chat(Lhat)
impute_and_stabilize_fit$Chat
class(impute_and_stabilize_fit$Chat)
dim(impute_and_stabilize_fit$Chat) # 4 fattori x 100 sogg
#View(impute_and_stabilize_fit$Chat)

# Estrai i fattori stimati (le Chat)
Chat_all <- all_data_fit$Chat
Chat_IS  <- impute_and_stabilize_fit$Chat

# ---- Funzione per creare il plot per un dato fattore k ----
plot_latent_factor <- function(k, true_C, all_data_fit, impute_and_stabilize_fit) {
  
  df_plot <- data.frame(
    True = rep(true_C[k,], 2),
    Estimated = c(all_data_fit$Chat[k,], impute_and_stabilize_fit$Chat[k,]),
    Method = rep(c("All Data", "Impute & Stabilize"), each = length(true_C[k,]))
  )
  
  ggplot(df_plot, aes(x = True, y = Estimated, color = Method)) +
    geom_point(alpha = 0.6, size = 2) +  # punti semitrasparenti
    geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +  # linea 1:1
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 13),
      legend.position = "top",
      legend.title = element_blank()
    ) +
    scale_color_manual(
      values = c("All Data" = "blue", "Impute & Stabilize" = "darkorange")
    ) +
    labs(
      x = paste("True latent factor", k),
      y = paste("Estimated latent factor", k)
    )
}


# Fattore 1
plot_latent_factor(1, true_C, all_data_fit, impute_and_stabilize_fit)

# Fattore 2
plot_latent_factor(2, true_C, all_data_fit, impute_and_stabilize_fit)

# Fattore 3
plot_latent_factor(3, true_C, all_data_fit, impute_and_stabilize_fit)

# Fattore 4
plot_latent_factor(4, true_C, all_data_fit, impute_and_stabilize_fit)


# C) ATE E CONFRONTO DEI DUE ALGORITMI 
true_ate <- list(ATE = ATE) # aggiungo ATE simulato

res_list <- list(
  "All Data" = all_data_fit,
  "Impute and Stabilize" = impute_and_stabilize_fit,
  "True values" = true_ate
)

plot_causalLFO_results(res_list) + scale_color_manual(
  values = c(
    "All Data"                 = "#e34a33",  
    "Impute and Stabilize"     = "#0570b0",  
    "True values"              = "#77dd77"  
  )
)



####################################################### PROVA SENZA REFERENCE_P
impute_and_stabilize_fit_noP <- impute_and_stabilize(
  M, Tr, rank = 4)
summary(impute_and_stabilize_fit_noP)

all_data_fit_noP <- all_data(
  M, Tr, rank = 4)
summary(all_data_fit_noP)

res_list <- list(
  'All Data' = all_data_fit,
  'All Data noP' = all_data_fit_noP,
  'Impute and Stabilize' = impute_and_stabilize_fit,
  'Impute and Stabilize noP' = impute_and_stabilize_fit_noP
)

plot_causalLFO_results(res_list) + scale_color_manual(
  values = c(
    "All Data"                 = "#e34a33",  
    "All Data noP"             = "#fdbb84",  
    "Impute and Stabilize"     = "#0570b0",  
    "Impute and Stabilize noP" = "#74a9cf"  
  )
)


############################# APPLICAZIONE SU 100 DATASET SIMULATI (usando le funzioni in simulate.r)

source("simulate.R")




