

lambda1 = 0.9
lambda2 = -3
gamma1 = 1.2
gamma2 = -0.52



generate_data <- function(n, lambda1, lambda2,
                          gamma1 , gamma2 ,seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  expit <- plogis
  
  W1 <- rbeta(n, 2, 3)
  W2 <- rbeta(n, 2, 3)
  W3 <- rbeta(n, 2, 3)
  
  A <- rbinom(n, 1, expit(0.5*W1 + 0.5*W2 - 1))
  Z <- rbinom(n, 3, expit(-1.7 + 1.5*A + 0.5*W3^2))
  M <- rbinom(n, 3, expit(-1.5 + lambda1*Z + lambda2*A + 0.4*W2 + 0.2*W3))
  Y <- rbinom(n, 1, expit(0.4*M + gamma1*Z + gamma2*A - 0.5*cos(W1) - 1.5))
  
  Data <- data.frame(W1, W2, W3, A, Z, M, Y)
  return(Data)
}


fit_models <- function(Data) {
  
  f <- glm(Y ~ A + M + Z + W1 + W2 + W3,
           data = Data, family = binomial())
  
  g <- glm(cbind(M, 3-M) ~ A + Z + W1 + W2 + W3,
           data = Data, family = binomial())
  
  h <- glm(cbind(Z, 3-Z) ~ A + W1 + W2 + W3,
           data = Data, family = binomial())
  
  return(list(f = f, g = g, h = h))
}
Data <- generate_data(1000,lambda1 = lambda1, lambda2 = lambda2,
                      gamma1 = gamma1, gamma2 = gamma2, seed=1)
models <- fit_models(Data)



g_formula_RT <- function(a, aprime, atilde, astar, data, models) {
  S <- numeric(nrow(data))
  triplets <- expand.grid(m = 0:3, z = 0:3, zstar = 0:3)
  
  for (i in seq_len(nrow(triplets))) {
    m <- triplets$m[i]; z <- triplets$z[i]; zstar <- triplets$zstar[i]
    
    df_f <- transform(data, A = a, M = m, Z = z)
    df_g <- transform(data, A = atilde, Z = zstar)
    df_z <- transform(data, A = aprime)
    df_zstar <- transform(data, A = astar)
    
    mu <- predict(models$f, newdata = df_f, type = "response")
    pm <- predict(models$g, newdata = df_g, type = "response")
    pz <- predict(models$h, newdata = df_z, type = "response")
    pzstar <- predict(models$h, newdata = df_zstar, type = "response")
    
    S <- S + mu *
      dbinom(m, size = 3, prob = pm) *
      dbinom(z, size = 3, prob = pz) *
      dbinom(zstar, size = 3, prob = pzstar)
  }
  
  mean(S)
}


g_formula_S1<-function(a, data, models){
  f <- glm(Y ~ A + W1 + W2 + W3,
           data = Data, family = binomial())
  df_f <- transform(data, A = a)
  mu <- predict(f, newdata = df_f, type = "response")
  mean(mu)
}


g_formula_natural<-function(a,aprime, atilde, data, models){
  duo <- expand.grid(m = 0:3, z = 0:3)
  f <- glm(Y ~ A+M + W1 + W2 + W3,
           data = Data, family = binomial())
  g <- glm(cbind(M, 3-M)  ~ A +  + W1 + W2 + W3,
           data = Data, family = binomial())
  h <- glm(cbind(Z, 3-Z) ~ A + W1 + W2 + W3,
           data = Data, family = binomial())
  
  for (i in seq_len(nrow(duo))) {
    m <- duo$m[i]
    z <- duo$z[i]
    
    df_f <- transform(data, A = a, M = m, Z = z)
    df_g <- transform(data, A = atilde, Z = z)
    df_z <- transform(data, A = aprime)
    
    mu <- predict(models$f, newdata = df_f, type = "response")
    pm <- predict(models$g, newdata = df_g, type = "response")
    pz <- predict(models$h, newdata = df_z, type = "response")
    
    S <- S + mu *
      dbinom(m, size = 3, prob = pm) *
      dbinom(z, size = 3, prob = pz) 
  }
  
  mean(S)
}
  





true_g_formula_RT <- function(a, aprime, atilde, astar, n_mc = 1e6,
                              lambda1, lambda2,
                              gamma1 , gamma2 , seed = 123) {
  
  set.seed(seed)
  W1 <- rbeta(n_mc, 2, 3)
  W2 <- rbeta(n_mc, 2, 3)
  W3 <- rbeta(n_mc, 2, 3)
  
  S <- numeric(n_mc)
  triplets <- expand.grid(m = 0:3, z = 0:3, zstar = 0:3)
  
  for (i in seq_len(nrow(triplets))) {
    m <- triplets$m[i]; z <- triplets$z[i]; zstar <- triplets$zstar[i]
    
    mu <- plogis(0.4*m + gamma1*z + gamma2*a - 0.5*cos(W1) - 1.5)
    pm <- plogis(-1.5 + lambda1*zstar + lambda2*atilde + 0.4*W2 + 0.2*W3)
    pz <- plogis(-1.7 + 1.5*aprime + 0.5*W3^2)
    pzstar <- plogis(-1.7 + 1.5*astar + 0.5*W3^2)
    
    S <- S + mu *
      dbinom(m, 3, pm) *
      dbinom(z, 3, pz) *
      dbinom(zstar, 3, pzstar)
  }
  
  mean(S)
}



true_g_formula_natural<- function(a, aprime, atilde, n_mc = 1e6,
                                  lambda1, lambda2,
                                  gamma1 , gamma2 ,seed = 123) {
  
  set.seed(seed)
  W1 <- rbeta(n_mc, 2, 3)
  W2 <- rbeta(n_mc, 2, 3)
  W3 <- rbeta(n_mc, 2, 3)
  
  S <- numeric(n_mc)
  triplets <- expand.grid(m = 0:3, z = 0:3)
  
  for (i in seq_len(nrow(triplets))) {
    m <- triplets$m[i]; z <- triplets$z[i]
    
    mu <- plogis(0.4*m + gamma1*z + gamma2*a - 0.5*cos(W1) - 1.5)
    pm <- plogis(-1.5 + lambda1*z + lambda2*atilde + 0.4*W2 + 0.2*W3)
    pz <- plogis(-1.7 + 1.5*aprime + 0.5*W3^2)
    
    S <- S + mu *
      dbinom(m, 3, pm) *
      dbinom(z, 3, pz) 
  }
  
  mean(S)
}

#psi_true <- true_g_formula(1, 0, 1, 0)

#psi_hat <- g_formula(1, 0, 1, 0, data = Data, models = models)

#c(True = psi_true, Estimated = psi_hat, Bias = psi_hat - psi_true)

##Test recanting twin : 

## S0 1,1,1,1
## S1 0,1,1,1
## S2 0,0,1,1
## S3 0,0,1,0
## S4 0,0,0,0


psi_2_3 <- function(lambda1 , lambda2,
                    gamma1, gamma2 ){
  true_g_formula_natural(1,1,1,lambda1 = lambda1, lambda2 = lambda2,
                         gamma1 = gamma1, gamma2 = gamma2, n_mc = 1e6)- true_g_formula_RT(1,1,1,1,lambda1 = lambda1, lambda2 = lambda2,
                                                                                          gamma1 = gamma1, gamma2 = gamma2, n_mc = 1e4) + true_g_formula_RT(0,0,1,0,lambda1 = lambda1, lambda2 = lambda2,
                                                                                                                                                            gamma1 = gamma1, gamma2 = gamma2, n_mc = 1e4) - true_g_formula_natural(0,0,1,lambda1 = lambda1, lambda2 = lambda2,
                                                                                                                                                                                                                                   gamma1 = gamma1, gamma2 = gamma2, n_mc = 1e4)
}



S1S1 <- function(lambda1 , lambda2,
                    gamma1, gamma2 ){
  true_g_formula_natural(1,1,1,lambda1 = lambda1, lambda2 = lambda2,
                         gamma1 = gamma1, gamma2 = gamma2, n_mc = 1e6)- true_g_formula_RT(1,1,1,1,lambda1 = lambda1, lambda2 = lambda2,
                                                                                          gamma1 = gamma1, gamma2 = gamma2, n_mc = 1e4) 
}

S3S3 <- function(lambda1 , lambda2,
                 gamma1, gamma2 ){
  true_g_formula_natural(0,0,1,lambda1 = lambda1, lambda2 = lambda2,
                         gamma1 = gamma1, gamma2 = gamma2, n_mc = 1e6)- true_g_formula_RT(0,0,1,0,lambda1 = lambda1, lambda2 = lambda2,
                                                                                          gamma1 = gamma1, gamma2 = gamma2, n_mc = 1e4) 
}




S3S3(lambda1 = lambda1, lambda2 = lambda2,
     gamma1 = gamma1, gamma2 = gamma2)


S1S1(lambda1 = lambda1, lambda2 = lambda2,
     gamma1 = gamma1, gamma2 = gamma2)
  

psi_2_3(lambda1 = lambda1, lambda2 = lambda2,
        gamma1 = gamma1, gamma2 = gamma2)


### Looking for a DGP that gives psi2v3 =0:

