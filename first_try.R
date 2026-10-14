
generate_data <- function(n, theta = 0.4, rho = 0.65, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  Data <- data.frame(
    W1 = rbeta(n,2,3), W2 = rbeta(n,2,3), W3 = rbeta(n,2,3)
  )
  
  Data$A <- with(Data, rbinom(n,1,plogis(0.5*W1+0.5*W2-1)))
  
  p0 <- with(Data, plogis(-1.7+0.5*W3^2))
  p1 <- with(Data, plogis(-0.2+0.5*W3^2))
  p11 <- p0*p1 + rho*(pmin(p0,p1)-p0*p1)
  p00 <- 1-p0-p1+p11
  p01 <- p1-p11
  
  U <- runif(n)
  Data$Z0 <- as.integer(U >= p00+p01)
  Data$Z1 <- as.integer((U >= p00 & U < p00+p01) | U >= 1-p11)
  Data$Z <- ifelse(Data$A == 0, Data$Z0, Data$Z1)
  
  Data$M <- with(Data, rbinom(n,3,plogis(
    -1.5+1.2*Z+1.5*A+0.4*W2+0.2*W3
  )))
  
  pY <- with(Data, 0.13+0.06*A+0.10*Z+0.11*(M/3)^2+
               0.07*sin(pi*W2)^2+0.05*W3^2+
               0.18*(1+0.3*A)*(W1-theta)*Z*(M/3)^2)
  
  Data$Y <- rbinom(n,1,pY)
  Data
}

Data <- generate_data(10000, theta = 0.4, seed = 123)

                          
                          
fit_models <- function(Data) {
  
  f <- glm(Y ~ A + M + Z + W1 + W2 + W3,
           data = Data, family = binomial())
  
  g <- glm(cbind(M, 3-M) ~ A + Z + W1 + W2 + W3,
           data = Data, family = binomial())
  
  h <- glm(Z ~ A + W1 + W2 + W3,
           data = Data, family = binomial())
  
  return(list(f = f, g = g, h = h))
}

models <- fit_models(Data)



g_formula_RT <- function(a, aprime, atilde, astar, data, models) {
  S <- numeric(nrow(data))
  triplets <- expand.grid(m = 0:3, z = 0:1, zstar = 0:1)
  
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
      dbinom(z, size = 1, prob = pz) *
      dbinom(zstar, size = 1, prob = pzstar)
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
  duo <- expand.grid(m = 0:3, z = 0:1)
  f <- glm(Y ~ A+M + W1 + W2 + W3,
           data = Data, family = binomial())
  g <- glm(cbind(M, 3-M)  ~ A +  + W1 + W2 + W3,
           data = Data, family = binomial())
  h <- glm(Z ~ A + W1 + W2 + W3,
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
      dbinom(z, size = 1, prob = pz) 
  }
  
  mean(S)
}
  





true_g_formula_RT <- function(a, aprime, atilde, astar,theta=0.4, n_mc = 1e6,
                                seed = 123) {
  
  set.seed(seed)
  W1 <- rbeta(n_mc, 2, 3)
  W2 <- rbeta(n_mc, 2, 3)
  W3 <- rbeta(n_mc, 2, 3)
  
  S <- numeric(n_mc)
  triplets <- expand.grid(m = 0:3, z = 0:1, zstar = 0:1)
  
  for (i in seq_len(nrow(triplets))) {
    m <- triplets$m[i]; z <- triplets$z[i]; zstar <- triplets$zstar[i]

    mu <- 0.13+0.06*a+0.10*z+0.11*(m/3)^2+
      0.07*sin(pi*W2)^2+0.05*W3^2+
      0.18*(1+0.3*a)*(W1-theta)*z*(m/3)^2
    pm <- plogis(
      -1.5+1.2*zstar+1.5*atilde+0.4*W2+0.2*W3
    )
    pz <- plogis(-1.7 + 1.5*aprime + 0.5*W3^2)
    pzstar <- plogis(-1.7 + 1.5*astar + 0.5*W3^2)
    S <- S + mu *
      dbinom(m, 3, pm) *
      dbinom(z, 1, pz) *
      dbinom(zstar, 1, pzstar)
  }
  
  mean(S)
}



true_g_formula_natural<- function(a, aprime, atilde,theta=0.4, n_mc = 1e6,
                                  seed = 123) {
  
  set.seed(seed)
  W1 <- rbeta(n_mc, 2, 3)
  W2 <- rbeta(n_mc, 2, 3)
  W3 <- rbeta(n_mc, 2, 3)
  
  S <- numeric(n_mc)
  triplets <- expand.grid(m = 0:3, z = 0:1)
  
  for (i in seq_len(nrow(triplets))) {
    m <- triplets$m[i]; z <- triplets$z[i]
    
    mu <- 0.13+0.06*a+0.10*z+0.11*(m/3)^2+
      0.07*sin(pi*W2)^2+0.05*W3^2+
      0.18*(1+0.3*a)*(W1-theta)*z*(m/3)^2
    pm <- plogis(
      -1.5+1.2*z+1.5*atilde+0.4*W2+0.2*W3
    )
    pz <- plogis(-1.7 + 1.5*aprime + 0.5*W3^2)
    
    S <- S + mu *
      dbinom(m, 3, pm) *
      dbinom(z, 1, pz) 
  }
  
  mean(S)
}





true_g_formula_natural_ind<- function(a, aprime, atilde,astar,theta=0.4,rho=0.65, n_mc = 1e6,
                                  seed = 123){
  
  set.seed(seed)
  W1 <- rbeta(n_mc, 2, 3)
  W2 <- rbeta(n_mc, 2, 3)
  W3 <- rbeta(n_mc, 2, 3)
  
  S <- numeric(n_mc)
  triplets <- expand.grid(m = 0:3, z = 0:1, zstar = 0:1)
  
  for (i in seq_len(nrow(triplets))) {
    m <- triplets$m[i]; z <- triplets$z[i]; zstar <- triplets$zstar[i]
    
    mu <- 0.13+0.06*a+0.10*z+0.11*(m/3)^2+
      0.07*sin(pi*W2)^2+0.05*W3^2+
      0.18*(1+0.3*a)*(W1-theta)*z*(m/3)^2
    pm <- plogis(
      -1.5+1.2*zstar+1.5*atilde+0.4*W2+0.2*W3
    )
    
    p0 <-  plogis(-1.7+0.5*W3^2)
    p1 <-  plogis(-0.2+0.5*W3^2)
    p11 <- p0*p1 + rho*(pmin(p0,p1)-p0*p1)
    p00 <- 1-p0-p1+p11
    p01 <- p1-p11
    p10 <- p0-p0*p1 + rho*(pmin(p0,p1)-p0*p1)
    
    
    if (aprime == 0 && astar == 1) {
      pzzstar <- list(p00, p01, p10, p11)[[2*z + zstar + 1]]
      
    } else if (aprime == 1 && astar == 0) {
      pzzstar <- list(p00, p10, p01, p11)[[2*z + zstar + 1]]
      
    } else {
      pzzstar <- if (z == zstar) dbinom(z, 1, if (aprime == 0) p0 else p1) else 0
    }
    
    
    S <- S + mu *
      dbinom(m, 3, pm) *
      pzzstar
  }
  
  mean(S)
}



psi_true <- true_g_formula_RT(1, 0, 1, 1)

psi_hat <- g_formula_RT(1, 0, 1, 1, data = Data, models = models)

c(True = psi_true, Estimated = psi_hat, Bias = psi_hat - psi_true)

##Test recanting twin : 

## S0 1,1,1,1
## S1 0,1,1,1
## S2 0,0,1,1
## S3 0,0,1,0
## S4 0,0,0,0


psi_2_3 <- function(){
  true_g_formula_natural(1,1,1,n_mc = 1e4)- true_g_formula_RT(1,1,1,1, n_mc = 1e4) + true_g_formula_RT(0,0,1,0,n_mc = 1e4) - true_g_formula_natural(0,0,1, n_mc = 1e4)
}



S1S1 <- function( ){
  true_g_formula_natural(1,1,1, n_mc = 1e4)- true_g_formula_RT(1,1,1,1, n_mc = 1e4) 
}

S3S3 <- function(){
  true_g_formula_natural(0,0,1,n_mc = 1e4)- true_g_formula_RT(0,0,1,0, n_mc = 1e4) 
}


S2S2 <- function( ){
  true_g_formula_natural_ind(0,0,1,1 ,n_mc = 1e4)- true_g_formula_RT(0,0,1,1, n_mc = 1e4) 
}






S3S3()


S1S1()
  
S2S2()
psi_2_3()


### Looking for a DGP that gives psi2v3 =0:

