set.seed(123)

library(MASS)

n <- 100
eff_fum <- 0.5

# Matrice di correlazione desiderata
cor_mat <- matrix(c(
  1.00,  0.00, -0.30,  0.20,  # eta
  0.00,  1.00,  0.30,  -0.30,  # status socio-economico
  -0.30,  0.30,  1.00,  eff_fum,  # abilità
  0.20,  -0.30,  eff_fum,  1.00   # logit fumo
), nrow = 4, byrow = TRUE)

# Nomi delle variabili
colnames(cor_mat) <- rownames(cor_mat) <- c("eta", "ses", "abilita", "logit_fumo")

# Simulazione
df <- mvrnorm(n, mu = rep(0, 4), Sigma = cor_mat, empirical = T)

# Estrai variabili
eta           <- df[, "eta"]
ses           <- df[, "ses"]
abilita       <- df[, "abilita"]
logit_fumo    <- df[, "logit_fumo"]

# Trasforma logit in probabilità di fumare
p_fumo   <- 1 / (1 + exp(-logit_fumo))
fumatore <- rbinom(n, size = 1, prob = p_fumo)

# Costruisci il dataset finale
df <- data.frame(
  eta,
  ses,
  abilita,  # opzionale: scala in stile QI
  logit_fumo,
  fumatore = factor(fumatore, labels = c("no", "sì"))
)

# Controllo delle correlazioni
round(cor(df[, 1:4]), 2)
t.test(abilita ~ fumatore, df = df)

head(df)

# senza aggiustamento
m1 <- lm(abilita ~ fumatore, data = df)
summary(m1)
# aggiustamento additivo
m2 <- update(m1, . ~ . + ses + eta, data = df)
summary(m2)

m <- glm(fumatore ~ ses + eta, data = df, family = binomial())

ps <- exp(predict(m3))/(1+exp(predict(m3)))
iptw <- ifelse(df$fumatore == "sì", 1/ps, 1/(1-ps))

# solo ps
m3 <- lm(abilita ~ fumatore, weights =  iptw, data = df)
summary(m3)

# additivo e ps
m4 <- update(m3,.  ~ . + ses + eta, data = df)
summary(m4)

BIC(m2,m3,m4)
AIC(m2,m3,m4)

car::compareCoefs(m1,m2,m4,m5)

