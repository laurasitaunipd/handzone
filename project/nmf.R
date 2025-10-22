# installare
#install.packages("remotes")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("Biobase")
#remotes::install_github("jennalandy/causalLFO")


# pacchetti
library(NMF)
library(causalLFO)
library(tidyverse)
library(ggridges)
library(ggplot2)


############################# APPLICAZIONE SU 1 DATASET SIMULATO CON N=100
# A) SIMULAZIONE DEL DATASET
rm(list=ls())
set.seed(321)
N = 100; D = 96; K = 3; ATE = c(1000, 0, 0)

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
true_C[3,sample(1:N, 10)] <- rnorm(10, mean = 1500, sd = 1000) # outliers for factor 3
data.frame(t(true_C)) %>%
  setNames(paste0("k", 1:K)) %>%   # rinomina le colonne come k1, k2, k3
  pivot_longer(1:K, names_to = 'k', values_to = 'C') %>%
  ggplot(aes(x = C, y = as.factor(k))) +
  geom_density_ridges() +
  theme_bw() +
  labs(x = "Untreated latent outcome distribution", y = "Latent dimension")

# Add ATE to loadings of treated samples
for (k in 1:K) {
  true_C[k, Tr == 1] <- true_C[k, Tr == 1] + ATE[k]
}

# comparing and visualizing true_c treated vs untreated
df <- data.frame(t(true_C)) %>% 
  setNames(paste0("k", 1:K)) %>%   # rinomina le colonne come k1, k2, k3
  mutate(Tr = factor(Tr, levels = c(0,1), labels = c("Untreated", "Treated"))) %>% 
  pivot_longer(cols = paste0("k", 1:K), names_to = "k", values_to = "C")

plotxslide <- ggplot(df, aes(x = C, y = as.factor(k), fill = Tr)) +
  geom_density_ridges(alpha = 0.6) +
  theme_bw() +
  labs(x = "Latent outcome distribution", 
       y = "Latent dimension",
       fill = "Group")

plotxslide

# Simulate M ~ Poisson(PC)
M = matrix(nrow = D, ncol = N)
for (i in 1:N) {
  M[,i] <- rpois(D, lambda = true_P %*% true_C[,i])
}

#View(M)

# B) STIME DI TRUE_P(LAMBDA) E TRUE_C(L)







############################# APPLICAZIONE SU 100 DATASET SIMULATI (usando le funzioni in simulate.r)

source("/Users/laura/Desktop/GITHUB/handzone/project/simulate.R")


