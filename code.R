#### Simulation Study ####

N = 2000 # num of data samples
J = 5 # num of TF
K = 5 # num of conditions
enrichment = c(0.5,0.1,0.7,0.02,0.4) # ChIPSeq enrichment for each TF


#
set.seed(4)

coef = matrix(sign(runif(K*(J+1))-0.5)*runif(K*(J+1), min = 0.5, max = 2), nrow = K)
coef[,3] = 0 # 2nd TF is NOT a regulator
coef[,5] = 0 # 4th TF is NOT a regulator
coef[,6] = 0 # 5th TF is NOT a regulator

standev = runif(K)

# ratio of dataset from each of K condition
ratio = runif(K)
ratio = ratio/sum(ratio)

for (k in 1:K) {
  
  n = floor(ratio[k]*N)
  
  intercpt = matrix(1, n, 1)
  TF = matrix(rnorm(n*J,5,1), ncol=J)
  
  TF = cbind(intercpt, TF)
  gene = TF %*% coef[k,] + rnorm(n, mean = 0, sd = standev[k])
  
  if (k==1) {
    X = TF
    y = gene
  } else{
    X = rbind(X,TF)
    y = rbind(y, gene)
  }
}


#### Gibbs sampler #####

# data dimension
N = dim(X)[1]
P = dim(X)[2]


# Instantiate params
g1 = 100 # hyperparameters
g2 = 1
m1 = 1
m2 = 1
s2 = 10
a = rep(1.5, J)
b = 1.5*(sum(enrichment)-2*enrichment)/(2*enrichment)


lambda = rgamma(K, shape=1.5, rate=1)
lambda = lambda/sum(lambda) # lambda sampled from Dirichlet(1.5)
e = sample(1:K, N, replace=T, prob = lambda)
gamma = matrix(runif(K*J), nrow = K)
r = matrix(rbinom(K*J, 1, 0.5), nrow = K)
c = 1/rgamma(K, shape=g1, rate=g2)

beta = matrix(rnorm(K*(J+1)), nrow = K) # plus 1 intercept 
sigma2 = rgamma(K, shape=2, rate=2)





Niter = 2000

# list to store sampled param
lambda_l = vector("list", length=Niter)
e_l = vector("list", length=Niter)
gamma_l = vector("list", length=Niter)
r_l = vector("list", length=Niter)
c_l = vector("list", length=Niter)
beta_l = vector("list", length=Niter)
sigma2_l = vector("list", length=Niter)


lambda_l[[1]] = lambda
e_l[[1]] = e
gamma_l[[1]] = gamma
r_l[[1]] = r
c_l[[1]] = c
beta_l[[1]] = beta
sigma2_l[[1]] = sigma2


# proposed param storage
lambda_bar = numeric(K)
e_bar = integer(N)
gamma_bar = matrix(numeric(K*J), nrow = K)
r_bar = matrix(numeric(K*J), nrow = K)
c_bar = numeric(K)

beta_bar = matrix(numeric(K*(J+1)), nrow = K) # plus 1 intercept 
sigma2_bar = numeric(K)


for (ii in 2:Niter) {
  
  # print current iteration number
  cat(sprintf("Current iteration:%s/%s\n",ii,Niter))
  
  # sample lambda: sample from dirichlet
  for (k in 1:K) {
    lambda_bar[k] = rgamma(1, shape=1.5+sum(e==k), rate=1)
  }
  lambda_bar = lambda_bar/sum(lambda_bar)
  lambda_l[[ii]] = lambda_bar
  
  
  # sample e: sample from categorical
  for (n in 1:N) {
    proba = log(lambda_bar)-0.5*log(sigma2_l[[ii-1]])-
      0.5*(rep(y[n],K)-X[n,]%*%t(beta_l[[ii-1]]))^2/sigma2_l[[ii-1]]
    
    e_bar[n] = sample(1:K, 1, replace=T, prob = exp(proba))
  }
  e_l[[ii]] = e_bar
  
  
  # sample gamma
  for (k in 1:K) {
    for (j in 1:J) {
      r_kj = r_l[[ii-1]][k,j]
      gamma_bar[k,j] = rbeta(1, r_kj+a[j], b[j]-r_kj+1)
    }
  }
  gamma_l[[ii]] = gamma_bar
  
  
  # sample r: Metropolis algo using Bernoulli(0.5) as kernel
  for (k in 1:K) {
    for (j in 1:J) {
      r_kj = r_l[[ii-1]][k,j]
      r_kj_p = rbinom(1, 1, 0.5)
      
      gamma_kj = gamma_bar[k,j]
      c_k = c_l[[ii-1]][k]
      beta_kj = beta_l[[ii-1]][k,j]
      
      logA = (r_kj_p+a[j]-1)*log(gamma_kj)+(b[j]-r_kj_p)*log(1-gamma_kj)-
        (1-r_kj_p)/2*log(c_k)-beta_kj^2/(2*c_k^(1-r_kj_p)*s2) - 
        ((r_kj+a[j]-1)*log(gamma_kj)+(b[j]-r_kj)*log(1-gamma_kj)-
        (1-r_kj)/2*log(c_k)-beta_kj^2/(2*c_k^(1-r_kj)*s2))
      
      if (runif(1) < min(1, exp(logA))) {r_kj = r_kj_p}
      
      r_bar[k,j] = r_kj
    }
  }
  r_l[[ii]] = r_bar
  
  
  
  # sample c
  for (k in 1:K) {
    r_k = r_bar[k,]
    beta_k = beta_l[[ii-1]][k,2:(J+1)]
    
    c_bar[k] = 1/rgamma(1, g1+(J-sum(r_k))/2, g2+sum((beta_k[r_k==0])^2)/(2*s2))
  }
  c_l[[ii]] = c_bar
  
  
  
  # sample beta
  for (k in 1:K) {
    r_k = r_bar[k]
    c_k = c_bar[k]
    sigma2_k = sigma2_l[[ii-1]][k]
    sigma_inv = diag(1/(rep(c_k, J+1)^c(0,1-r_k)*s2))
    
    XtX = t(X[e_bar==k,]) %*% X[e_bar==k,]
    Xty = t(X[e_bar==k,]) %*% y[e_bar==k]
    
    mu = solve(1/sigma2_k*XtX + sigma_inv, Xty)
    mu = 1/sigma2_k * mu
    
    R = chol(sigma_inv+1/sigma2_k*XtX)
    
    beta_bar[k,] = mu + solve(R, rnorm(J+1))
  }
  
  beta_l[[ii]] = beta_bar
  
  
  
  # sample sigma2
  for (k in 1:K) {
    beta_k = beta_bar[k,]
    sigma2_bar[k] = 1/rgamma(1, m1+sum(e_bar==k)/2, m2+crossprod(y[e_bar==k]-X[e_bar==k,]%*%beta_k)/2)
  }
  sigma2_l[[ii]] = sigma2_bar
  
}


#### Plot ####
library(purrr)

# plot trace, histogram and autocorrelation of sigma2[1]
data = sigma2_l %>% map(~.x[1])
data = unlist(data)
sigma2 = data[1001:Niter]


par(mfrow = c(1,3))

plot(sigma2, type = 'l')
hist(sigma2, probability = T)
acf(sigma2, lag.max = 50, plot=T);



# plot distribution of beta
par(mfrow=c(K,J+1)) # row-wise subplot
for (k in 1:K) {
  for (j in 1:(J+1)) {
    data = beta_l %>% map(~.x[k,j]) # extract row k col j of each matrix in list
    
    data = unlist(data)
    data_o = data[1001:Niter]
    
    hist(data_o, probability=T, main=NULL, xlab = sprintf('condition %s beta %s',k,j-1))
    abline(v=mean(data_o), col='blue', lwd=2) #sample mean
    
    abline(v=coef[k,j], col='red', lwd=2) #ground truth
  }
}
