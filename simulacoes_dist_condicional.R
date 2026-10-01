# Função para gerar as distribuições condicionais
gera_amostras <- function(p, tam_amostral, var, mi_0, mi_1) {
  
  # Começar com a Bernoulli
  bernoulli <- rbinom(
    n = tam_amostral,
    size = 1,
    prob = p
  )
  
  # Gerar as médias
  medias <- ifelse(bernoulli == 0, mi_0, mi_1)
  
  # Com as médias, gerar as observações normais
  normal <- rnorm(
    n = tam_amostral,
    mean = medias,
    sd = sqrt(var)
  )
  
  # Retornar um data.frame com as observações
  observacoes <- data.frame(
    normal = normal,
    bernoulli = bernoulli
  )
  
  return(observacoes)
}


# Função para calcular as estimativas
gera_estimativas <- function(observacoes) {
  
  # Estimador de mi_0
  estimativa_mi_0 <- (
    sum(observacoes$normal * observacoes$bernoulli) -
      nrow(observacoes) * mean(observacoes$normal)
  ) / (
    nrow(observacoes) *
      (mean(observacoes$bernoulli) - 1)
  )
  
  # Estimador de mi_1
  estimativa_mi_1 <- (
    sum(observacoes$normal * observacoes$bernoulli)
  ) / (
    nrow(observacoes) * mean(observacoes$bernoulli)
  )
  
  # Estimador de sigma
  lambda <- sum(
    (1 - observacoes$bernoulli) *
      (observacoes$normal - estimativa_mi_0)^2 +
      observacoes$bernoulli *
      (observacoes$normal - estimativa_mi_1)^2
  )
  
  estimativa_sigma <- sqrt(
    lambda / nrow(observacoes)
  )
  
  # Retornar os estimadores
  estimativas <- data.frame(
    estimativa_mi_0 = estimativa_mi_0,
    estimativa_mi_1 = estimativa_mi_1,
    estimativa_sigma = estimativa_sigma
  )
  
  return(estimativas)
}


#agora, gerar as simulações

gera_simulacoes <- function(num_simulacoes, tam_amostral, var, mi_0, mi_1, p){
  
  valores_estimador_mi_0 <- numeric(num_simulacoes)
  valores_estimador_mi_1 <- numeric(num_simulacoes)
  valores_estimador_sigma <- numeric(num_simulacoes)
  
  
  for (i in 1:num_simulacoes){
    observacoes <- gera_amostras(p = p, tam_amostral = tam_amostral,
                                 var = var, mi_0 = mi_0, mi_1 = mi_1)
    estimativas <- gera_estimativas(observacoes)
    valores_estimador_mi_0[i] <- estimativas$estimativa_mi_0
    valores_estimador_mi_1[i] <- estimativas$estimativa_mi_1
    valores_estimador_sigma[i] <- estimativas$estimativa_sigma
  }
  
  estimativas <- data.frame(estimativa_mi_0 = valores_estimador_mi_0,
                            estimativa_mi_1 = valores_estimador_mi_1,
                            estimativa_sigma = valores_estimador_sigma
                            )
  
  estatisticas <- data.frame(media = sapply(estimativas, mean),
                             variancia = sapply(estimativas, var),
                             minimo = sapply(estimativas, min),
                             max = sapply(estimativas, max))
                             
  
  return(list(
    estimativas = estimativas,
    estatisticas = estatisticas
  ))
  

}