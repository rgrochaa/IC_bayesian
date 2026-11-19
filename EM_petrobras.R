#Algoritmo EM (parece certo)

em_mistura <- function(y, K, p, mu, sigma) {
  
  max_iter = 1000 
  tol = 0.001
  n <- length(y)
  loglik <- numeric(0)
  
  for (t in 1:max_iter) {
    
    #Passo E
    # dens[i, k] = p_k * phi_k(y_i) é o numerador do gamma, a densidade condicional
    dens <- matrix(NA_real_, n, K)
    
    for (k in 1:K) {
      dens[, k] <- p[k] * dnorm(y, mu[k], sigma[k]) #a cada iteração, o k tá fixo e ele vai preenchendo a coluna
      #variando o y
    }
    
    # fy[i] = f_Y(y_i) = soma da linha i, já vimos na teoria
    fy <- rowSums(dens)
    
    # gamma[i, k] = p_k * phi_k(y_i) / f_Y(y_i) o R casa a dimensão bonitinho 
    gamma <- dens / fy
    
    # log-verossimilhança observada em theta^(t), somamos porque o log transforma produtos em somas
    loglik[t] <- sum(log(fy))
    
    ## ---------- Passo M ----------
    # n_j efetivo: soma da coluna j de gamma; é o termo que aparece sempre, quantas observações
    #esperamos que sejam do grupo j, que nem escrevi nas anotações
    nj <- colSums(gamma)
    
    # novos pesos
    p <- nj / n
    
    # novas médias (média ponderada por gamma)
    mu_novo <- numeric(K)
    for (k in 1:K) {
      mu_novo[k] <- sum(gamma[, k] * y) / nj[k]
    }
    mu <- mu_novo
    
    # novos desvios padrão (usam as médias novas)
    sigma_novo <- numeric(K)
    for (k in 1:K) {
      soma_quadrados <- sum(gamma[, k] * (y - mu[k])^2)
      sigma_novo[k] <- sqrt(soma_quadrados / nj[k])
    }
    sigma <- sigma_novo
    
    #Critério de parada
    if (t > 1) {
      if (abs(loglik[t] - loglik[t - 1]) < tol) {
        break
      }
    }
  }
  
  resultado <- list(p = p, mu = mu, sigma = sigma, iter = t, gamma = gamma)
  
  return(resultado)

}

#Exemplo do algoritmo
set.seed(4)
n <- 1000
z <- sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
y <- numeric(n)
for (i in 1:n) {
  if (z[i] == 0) {
    y[i] <- rnorm(1, mean = 0, sd = 1)
  } else {
    y[i] <- rnorm(1, mean = 4, sd = 1.5)
  }
}

fit <- em_mistura(y, K = 2,
                  p = c(0.5, 0.5),
                  mu = c(-1, 3),
                  sigma = c(1, 4))

fit$p
fit$mu
fit$sigma
fit$iter







#Análise da Petrobras
dados <- getSymbols("PETR4.SA", src = "yahoo", from = "2018-01-01",
                    to = Sys.Date(), auto.assign = FALSE)

preco_fechamento <- dados$PETR4.SA.Close
log_retornos <- diff(log(preco_fechamento))
log_retornos <- na.omit(log_retornos)



y <- as.numeric(log_retornos)     #tem que ser vetor, não matriz        
s <- sd(y)

#Vamos usar bear, crab e bull. A maior parte do tempo é crab. A média no bear (bull) tem que ser menor (maior) 
#e vou assumir que tem mais volatilidade na queda (faz sentido?), por isso a maior variância.
p     <- c(0.2, 0.6, 0.2)
mu    <- c(-s / 2, 0, s / 2)
sigma <- c(1.5 * s, 0.5 * s, 1.0 * s)

fit <- em_mistura(y, K = 3, p = p, mu = mu, sigma = sigma)

crab <- which.min(fit$sigma) #o crab vai ser o grupo da menor variância
outros <- setdiff(1:3, crab) # as posições que sobram
bull <- outros[which.max(fit$mu[outros])] #dentre as posições que sobraram, qual tem a maior média?
bear <- outros[which.min(fit$mu[outros])] # e qual tem a menor?


#colocando os nomes bonitinhos
#São 3 grupos, o grupo do bear recebe o nome bear, e o mesmo para as outras.
nomes <- character(3)
nomes[bear] <- "bear"
nomes[crab] <- "crab"
nomes[bull] <- "bull"


resumo <- data.frame(
  regime = nomes[c(bear, crab, bull)],
  p      = fit$p[c(bear, crab, bull)],
  mu     = fit$mu[c(bear, crab, bull)],
  sigma  = fit$sigma[c(bear, crab, bull)]) ##A tabela cai sempre na ordem bear, crab e bull