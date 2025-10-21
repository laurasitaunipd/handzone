library(ggplot2)
library(dplyr)

set.seed(42)  # Per riproducibilità

# Numero di osservazioni ridotto
N <- 500

# Simuliamo coscienziosità e intelligenza come variabili indipendenti
coscienziosita <- rnorm(N, mean = 0, sd = 1)
intelligenza <- rnorm(N, mean = 0, sd = 1)

# L'ammissione dipende dalla somma di coscienziosità e intelligenza (forte selezione)
p_ammissione <- plogis(3 * (coscienziosita + intelligenza - 0.5))  # funzione logistica più estrema
ammissione <- rbinom(N, 1, p_ammissione)  # 1 = ammesso, 0 = non ammesso

# Creiamo un data frame
dati <- data.frame(coscienziosita, intelligenza, ammissione)

# Calcoliamo la correlazione nei due gruppi
cor_ammessi <- cor(dati$coscienziosita[dati$ammissione == 1], 
                   dati$intelligenza[dati$ammissione == 1])

cor_non_ammessi <- cor(dati$coscienziosita[dati$ammissione == 0], 
                       dati$intelligenza[dati$ammissione == 0])

summary(lm(data = dati, coscienziosita ~ intelligenza + ammissione))
summary(lm(data = dati, coscienziosita ~ intelligenza))

# Grafico con linee di regressione e correlazioni
ggplot(dati, aes(x = coscienziosita, y = intelligenza, color = as.factor(ammissione))) +
  geom_point(alpha = 0.6) +
  
  # Linea di regressione unica (su tutti i dati, senza condizionare)
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed") +
  
  # Linee di regressione separate per gli ammessi e non ammessi
  geom_smooth(method = "lm", se = FALSE) +
  
  # Testo con i coefficienti di correlazione
  annotate("text", x = 0, y = 4, label = paste0("Ammessi: ", round(cor_ammessi, 2)), color = "blue", size = 5) +
  annotate("text", x = 0, y = -4, label = paste0("Non Ammessi: ", round(cor_non_ammessi, 2)), color = "red", size = 5) +
  
  theme_minimal() +
  scale_color_manual(values = c("red", "blue"), labels = c("Non ammessi", "Ammessi")) +
  
  labs(
    title = "Bias del Collider: Coscienziosità vs Intelligenza",
    subtitle = "Le linee colorate sono condizionate all'ammissione, la linea nera è generale",
    color = "Ammissione",
    x = "Coscienziosità",
    y = "Intelligenza"
  )