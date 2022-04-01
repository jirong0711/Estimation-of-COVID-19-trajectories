jags.model="
model{
  for(i in 1:N){
    for(j in 1:Ts){
      Y[i,j] ~ dnorm(theta1[i]*((1+xi[i]*exp(-theta2[i]*(j-theta3[i])))^(-1/xi[i])), global.sigma)
    }
    theta1[i] ~ dnorm(alpha1+X[i,]%*%beta1, 1/sig1.sq)
    theta2[i] ~ dnorm(alpha2+X[i,]%*%beta2, 1/sig2.sq)
    theta3[i] ~ dnorm(alpha3+X[i,]%*%beta3, 1/sig3.sq)
    
    xi[i] ~ dlnorm(0,1)
  }
  for(k in 1:p){
  beta1[k] ~ dnorm(0, 1/(sig1.sq*tau1.sq[k]*lambda1[k]))
  beta2[k] ~ dnorm(0, 1/(sig2.sq*tau2.sq[k]*lambda2[k]))
  beta3[k] ~ dnorm(0, 1/(sig3.sq*tau3.sq[k]*lambda3[k]))
  
  lamb1[k] ~ dscaled.gamma(1,1)
  lambda1[k] = 1/sqrt(lamb1[k])
  lamb2[k] ~ dscaled.gamma(1,1)
  lambda2[k] = 1/sqrt(lamb2[k])
  lamb3[k] ~ dscaled.gamma(1,1)
  lambda3[k] = 1/sqrt(lamb3[k])  

  gamma.tau1[k] ~ dscaled.gamma(1,1)
  tau1.sq[k] = 1/sqrt(gamma.tau1[k])
  gamma.tau2[k] ~ dscaled.gamma(1,1)
  tau2.sq[k] = 1/sqrt(gamma.tau2[k])
  gamma.tau3[k] ~ dscaled.gamma(1,1)
  tau3.sq[k] = 1/sqrt(gamma.tau3[k])  
  }
  
  global.sigma ~ dgamma(.005, .005)
  alpha1 ~ dnorm(0,0.01)
  alpha2 ~ dnorm(0,0.01)
  alpha3 ~ dnorm(0,0.01)
  sig1.sq ~ dgamma(.005, .005)
  sig2.sq ~ dgamma(.005, .005)
  sig3.sq ~ dgamma(.005, .005)
}"

library(rjags)
load.module("glm")
dat <- list(Y=Y, N=dim(Y)[1], Ts=dim(Y)[2], X=X, p=dim(X)[2])
model <- jags.model(textConnection(jags.model),
                    data=dat, n.chains=3)
update(model, 10000)
samp <- coda.samples(model, variable.names=c("beta","sigma"),
                     n.iter=20000, progress.bar="none")