
#######################################
rm(list = ls())
N = 1e5

#######################################

z = rnorm(N,0,1)
x = round(pnorm(rnorm(N,0,1) + 0.7*z))
y = rnorm(N,0,1) + 0.4*x + 0.4*z

#######################################

# uncorrected

fit = lm(y ~ x)
summary(fit)

#######################################

# corrected

fit = lm(y ~ x + z)
summary(fit)

#######################################

# propensity scores

fitPs = glm(x ~ z, family=binomial(link="logit"))
pred = predict(fitPs,type="response")
pred[x==1] = 1/pred[x==1]
pred[x==0] = 1/(1-pred[x==0])

fit1 = lm(y ~ x)
summary(fit1)
fit2 = lm(y ~ x + z)
summary(fit2)
fit3 = lm(y ~ x, weights = pred)
summary(fit3)
fit4 = lm(y ~ x + z, weights = pred)
summary(fit4)


############################################ AGGIUNTE
p1 = tibble(car::compareCoefs(fit0,fit1,fit2,fit3))
p$sim = paste(i,j, sep = "_")
res = bind_rows(p,p1)
str(p)

####################################### SIM CHAT
rm(list = ls())
library(car)

# Numero di simulazioni
n_sim <- 100  # Cambia il numero di simulazioni a piacere

# Liste per salvare i risultati
coefs_list <- vector("list", n_sim)
se_list <- vector("list", n_sim)

for (i in 1:n_sim) {
  # Generazione dati
  N = 1e5
  z = rnorm(N, 0, 1)
  x = round(pnorm(rnorm(N, 0, 1) + 0.7 * z))
  y = rnorm(N, 0, 1) + 0.4 * x + 0.4 * z
  
  # Modelli di regressione
  fit1 = lm(y ~ x)
  fit2 = lm(y ~ x + z)
  
  fitPs = glm(x ~ z, family = binomial(link = "logit"))
  pred = predict(fitPs, type = "response")
  pred[x == 1] = 1 / pred[x == 1]
  pred[x == 0] = 1 / (1 - pred[x == 0])
  
  fit3 = lm(y ~ x, weights = pred)
  fit4 = lm(y ~ x + z, weights = pred)
  
  # Estrarre coefficienti e standard error
  models <- list(fit1, fit2, fit3, fit4)
  coefs <- lapply(models, coef)  # Lista di coefficienti
  ses <- lapply(models, function(m) sqrt(diag(vcov(m))))  # Lista di standard error
  
  coefs_list[[i]] <- coefs
  se_list[[i]] <- ses
}

# Funzione per mediare le liste di coefficienti/SE
mean_coef_se <- function(list_of_lists) {
  keys <- unique(unlist(lapply(list_of_lists, names)))  # Trova tutti i nomi dei coefficienti
  mean_list <- lapply(1:length(list_of_lists[[1]]), function(j) {
    temp <- do.call(rbind, lapply(list_of_lists, function(l) {
      vec <- rep(NA, length(keys))  # Crea un vettore pieno di NA
      names(vec) <- keys
      vec[names(l[[j]])] <- l[[j]]  # Assegna solo i valori esistenti
      return(vec)
    }))
    colMeans(temp, na.rm = TRUE)  # Media ignorando gli NA
  })
  return(mean_list)
}

# Calcolare la media dei coefficienti e degli standard error
mean_coefs <- mean_coef_se(coefs_list)
mean_ses <- mean_coef_se(se_list)

# Creare modelli fittizi con i valori medi
fake_model <- function(original_fit, mean_coefs, mean_ses) {
  new_fit <- original_fit
  new_fit$coefficients <- mean_coefs
  attr(new_fit, "vcov") <- diag(mean_ses^2)  # Assegna direttamente la matrice di varianza-covarianza
  return(new_fit)
}

fit1_avg <- fake_model(fit1, mean_coefs[[1]], mean_ses[[1]])
fit2_avg <- fake_model(fit2, mean_coefs[[2]], mean_ses[[2]])
fit3_avg <- fake_model(fit3, mean_coefs[[3]], mean_ses[[3]])
fit4_avg <- fake_model(fit4, mean_coefs[[4]], mean_ses[[4]])

# Confronto finale con compareCoefs
compareCoefs(fit1_avg, fit2_avg, fit3_avg, fit4_avg)
