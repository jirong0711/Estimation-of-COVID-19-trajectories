COVID19_curve <- function(Y,X,seed.no=1,burn=2000,nmc=2000,thin=10, varrho = 0, pro.var.theta.2 = 0.0002, pro.var.theta.3 = 0.05, mu = 0, rho.sq = 1){
  # varrho: hyper-parameter for the diffuse prior error variances ~ IG(varrho,varrho)
  # pro.var.theta.2 and pro.var.theta.2: proposal variances for the MH algorithms to sample from theta2.i and theta3.i  
  library(mvtnorm)
  library(nleqslv)
  {
    N <- nrow(Y); Ts <- ncol(Y)
    p <- ncol(X) #covariates
    y <- function(i){
      return(as.numeric(Y[i,]))
    }
  }
  
  # MCMC Simulation Times
  {
    burn = burn; nmc = nmc
    thin = thin; B = burn + nmc
  }
  
  # Room for parameters
  {
    #Model parameters
    theta1.vec <- matrix(rep(0, N*B), nrow=N)
    theta2.vec <- matrix(rep(0, N*B), nrow=N)
    theta3.vec <- matrix(rep(0, N*B), nrow=N)
    xi.vec <- matrix(rep(0, N*B), nrow=N)
    global.sigma <- rep(0, B)
    
    #Model Prior
    alpha1<- rep(0, B)
    alpha2 <- rep(0, B)
    alpha3 <- rep(0, B)
    local.sigma1 <- rep(0, B)
    local.sigma2 <- rep(0, B)
    local.sigma3 <- rep(0, B)
    
    #Horseshoe Prior
    beta1.vec <- matrix(rep(0, p*B), nrow=p)
    beta2.vec <- matrix(rep(0, p*B), nrow=p)
    beta3.vec <- matrix(rep(0, p*B), nrow=p)
    lambda1.vec <- matrix(rep(0, p*B), nrow=p)
    lambda2.vec <- matrix(rep(0, p*B), nrow=p)
    lambda3.vec <- matrix(rep(0, p*B), nrow=p)
    tau1 <- rep(0, B)
    tau2 <- rep(0, B)
    tau3 <- rep(0, B)
  }
  
  
  #Global function
  norm.f <- function(x, nul=rep(0, length(x))){
    return(sqrt(sum((x-nul)^2)))
  }
  
  fcal <- function(t, theta1, theta2, theta3, xi){
    return(theta1*((1+xi*exp(-theta2*(t-theta3)))^(-1/xi)))
  }
  
  f.i.vec = function(i,theta1.i,theta2.i,theta3.i, xi.i){
    res = matrix(data = fcal(t=c(1:Ts), theta1 = theta1.i, theta2 = theta2.i, theta3 = theta3.i, xi = xi.i), nrow = Ts, ncol = 1) 
    return(res)
  } 
  
  #Global Variable
  onevec.N <- rep(1,N); I.N <- diag(onevec.N)
  
  # Functions for Step 1
  h <- function(t, theta2, theta3, xi){
    return((1+xi*exp(-theta2*(t-theta3)))^(-1/xi))
  }
  h.i.vec <- function(i, theta2.i, theta3.i, xi.i){
    return(matrix(h(t=1:Ts, theta2=theta2.i, theta3=theta3.i, xi=xi.i), nrow=Ts))
  }
  r <- function(i, theta2.i, theta3.i, xi.i){
    return(t(y(i))%*%h.i.vec(i, theta2.i, theta3.i, xi.i))
  }
  r.vec <- function(theta2.vec, theta3.vec, xi.vec){
    temp <- sapply(1:N, function(x) r(x, theta2.vec[x], theta3.vec[x], xi.vec[x]))
    return(as.matrix(temp))
  }
  
  H <- function(theta2.vec, theta3.vec, xi.vec){
    temp = c()
    for (i in 1:N) {
      temp[i] = norm.f(x = h.i.vec(i = i, theta2.i = theta2.vec[i], theta3.i = theta3.vec[i], xi.i = xi.vec[i]))^2
    }
    res = diag(temp)
    return(res)
  }
  
  
  #Functions for step2-theta2
  f.2.i <- function(i,theta1.i, old.theta2.i, theta2.i,theta3.i, xi.i){
    newf <- f.i.vec(i = i, theta2.i = theta2.i, theta1.i = theta1.i, theta3.i = theta3.i, xi.i = xi.i)
    ff <- f.i.vec(i = i, theta2.i = old.theta2.i, theta1.i = theta1.i, theta3.i = theta3.i, xi.i = xi.i)
    return(newf-ff)
  }
  
  s.2.i <- function(i,theta1.i, old.theta2.i, theta2.i,theta3.i, xi.i){
    newf <- norm.f(x=f.i.vec(i = i, theta2.i = theta2.i, theta1.i = theta1.i, theta3.i = theta3.i, xi.i = xi.i))^2
    ff <- norm.f(x=f.i.vec(i = i, theta2.i = old.theta2.i, theta1.i = theta1.i, theta3.i = theta3.i, xi.i = xi.i))^2
    return(newf-ff)
  }
  
  loglike.2.i <- function(i,theta1.i, old.theta2.i, theta2.i,theta3.i, xi.i, global.sigma){
    ff <- sum(t(y(i))%*%f.2.i(i=i,theta1.i=theta1.i, old.theta2.i=old.theta2.i, theta2.i=theta2.i,theta3.i=theta3.i, xi.i=xi.i))
    ss <- s.2.i(i=i,theta1.i=theta1.i, old.theta2.i=old.theta2.i, theta2.i=theta2.i,theta3.i=theta3.i, xi.i=xi.i)
    return((1/global.sigma)*(ff-(1/2)*ss))
  }
  
  theta2.ratio <- function(i,theta1.i, old.theta2.i, theta2.i,theta3.i, xi.i, global.sigma, alpha2, local.sigma2, beta2.vec){
    ratio <- exp(loglike.2.i(i=i,theta1.i=theta1.i, old.theta2.i=old.theta2.i, theta2.i=theta2.i,theta3.i=theta3.i, xi.i=xi.i, global.sigma = global.sigma)+
                   (1/local.sigma2)*((alpha2+X[i,]%*%beta2.vec)*(theta2.i-old.theta2.i)- (1/2)*(theta2.i^2-old.theta2.i^2)))
    return(max(min(ratio,1),1e-10))
  }
  
  #Functions for step2-theta3
  f.3.i <- function(i,theta1.i, theta2.i, old.theta3.i, theta3.i, xi.i){
    newf <- f.i.vec(i = i, theta2.i = theta2.i, theta1.i = theta1.i, theta3.i = theta3.i, xi.i = xi.i)
    ff <- f.i.vec(i = i, theta2.i = theta2.i, theta1.i = theta1.i, theta3.i = old.theta3.i, xi.i = xi.i)
    return(newf-ff)
  }
  
  s.3.i <- function(i,theta1.i, theta2.i, old.theta3.i, theta3.i, xi.i){
    newf <- norm.f(x=f.i.vec(i = i, theta2.i = theta2.i, theta1.i = theta1.i, theta3.i = theta3.i, xi.i = xi.i))^2
    ff <- norm.f(x=f.i.vec(i = i, theta2.i = theta2.i, theta1.i = theta1.i, theta3.i = old.theta3.i, xi.i = xi.i))^2
    return(newf-ff)
  }
  
  loglike.3.i <- function(i,theta1.i, theta2.i, old.theta3.i, theta3.i, xi.i, global.sigma){
    f <- sum(t(y(i))%*%f.3.i(i=i,theta1.i=theta1.i, theta2.i=theta2.i, old.theta3.i=old.theta3.i,theta3.i=theta3.i, xi.i=xi.i))
    s <- s.3.i(i=i,theta1.i=theta1.i, theta2.i=theta2.i, old.theta3.i=old.theta3.i,theta3.i=theta3.i, xi.i=xi.i)
    return((1/global.sigma)*(f-(1/2)*s))
  }
  
  theta3.ratio <- function(i,theta1.i,theta2.i, old.theta3.i,theta3.i, xi.i, global.sigma, alpha3, local.sigma3, beta3.vec){
    ratio <- exp(loglike.3.i(i=i,theta1.i=theta1.i, theta2.i=theta2.i, old.theta3.i=old.theta3.i,theta3.i=theta3.i, xi.i=xi.i, global.sigma = global.sigma)+
                   (1/local.sigma3)*((alpha3+X[i,]%*%beta3.vec)*(theta3.i-old.theta3.i)- (1/2)*(theta3.i^2-old.theta3.i^2)))
    return(max(min(ratio,1),1e-10))
  }
  
  #Functions for step3
  f.eta.i <- function(i, eta.i, old.eta.i, theta1.i,theta2.i, theta3.i){
    xi.i <- exp(eta.i) 
    old.xi.i <- exp(old.eta.i) 
    temp1 <- f.i.vec(i = i, xi.i = xi.i, theta1.i = theta1.i, theta2.i = theta2.i, theta3.i = theta3.i)
    temp2 <- f.i.vec(i = i, xi.i = old.xi.i, theta1.i = theta1.i, theta2.i = theta2.i, theta3.i = theta3.i)        
    return(temp1-temp2)
  }
  
  s.eta.i <- function(i, eta.i, old.eta.i, theta1.i,theta2.i, theta3.i){
    xi.i <- exp(eta.i) 
    old.xi.i <- exp(old.eta.i) 
    ff <- norm.f(x=f.i.vec(i = i, xi.i = xi.i, theta1.i = theta1.i, theta2.i = theta2.i, theta3.i = theta3.i))^2
    ss <- norm.f(x=f.i.vec(i = i, xi.i = old.xi.i, theta1.i = theta1.i, theta2.i = theta2.i, theta3.i = theta3.i))^2
    return(ff-ss)
  }
  
  loglike.eta.i <- function(i, eta.i, old.eta.i, theta1.i,theta2.i, theta3.i, global.sigma){
    ff <- sum(t(y(i))%*%f.eta.i(i = i,  old.eta.i=old.eta.i, eta.i=eta.i, theta1.i = theta1.i, theta2.i = theta2.i, theta3.i = theta3.i))
    ss <- s.eta.i(i = i,  old.eta.i=old.eta.i, eta.i=eta.i, theta1.i = theta1.i, theta2.i = theta2.i, theta3.i = theta3.i)
    temp <- (1/global.sigma)*(ff-(1/2)*ss)
    return(temp)
  }
  
  eta.ratio <- function(i, eta.i, old.eta.i, theta1.i,theta2.i, theta3.i, global.sigma){
    ratio <- exp(loglike.eta.i(i = i,  old.eta.i=old.eta.i, eta.i=eta.i, theta1.i = theta1.i, theta2.i = theta2.i, theta3.i = theta3.i, global.sigma = global.sigma))
    return(max(min(ratio,1),1e-10))
  }
  
  
  #Function for step6
  beta.rmvt = function(Q, b){
    #Sampling from multivariate gaussian distribution with cholesky decomposition
    p = dim(Q)[1]
    L = t(chol(Q))
    z = rnorm(n = p ,mean = 0,sd = 1)
    y = solve(t(L), z)
    v = solve(L, b)
    theta = solve(t(L), v)
    beta = y + theta
    return(beta)
  }
  
  #Function for step 7 & 8
  ss.invgam = function(old.x, f,s){
    library(truncdist)
    u = runif(n = 1, min = 0, max = old.x/(1+old.x)) 
    x = 1/rtrunc(n=1, spec = "gamma", shape = f, rate = s, a = 0, b = (1-u)/u )
    return(x)
  }
  
  #Set initial value
  {
    for (i in 1:N){
      temp.d = data.frame(y=y(i), t=c(1:Ts))
      theta1.vec[i,1] = max(temp.d$y)
      theta3.vec[i,1] = which.max(diff(temp.d$y))
      xi.vec[i,1] = 0.1
      f = function(t, theta1,theta2, theta3, xi){
        res = theta1/( 1 + xi*exp(-theta2*(t - theta3)))^(1/xi)
        return(res)
      }
      eps = 0.00000001
      temp.ft = function(x){
        nomi = log(f(t = theta3.vec[i,1] + eps, theta1 = theta1.vec[i,1], theta2 = x, theta3 = theta3.vec[i,1], xi = xi.vec[i,1])/
                     f(t = theta3.vec[i,1], theta1 = theta1.vec[i,1], theta2 = x, theta3 = theta3.vec[i,1], xi = xi.vec[i,1]))
        denom = eps
        res = x/2 - nomi/denom
        return(res)
      }
      theta2.vec[i,1] = nleqslv(x = 0.1, fn = temp.ft)$x
    }
    global.sigma[1] = 1
    alpha1[1] = 1
    alpha2[1] = 1
    alpha3[1] = 1
    local.sigma1[1] = 1
    local.sigma2[1] = 1
    local.sigma3[1] = 1

    beta1.vec[,1] = rep(0,p)
    beta2.vec[,1] = rep(0,p)
    beta3.vec[,1] = rep(0,p)
    lambda1.vec[,1] = rep(1,p)
    lambda2.vec[,1] = rep(1,p)
    lambda3.vec[,1] = rep(1,p)
    tau1[1] = 1
    tau2[1] = 1
    tau3[1] = 1
  }
  
  #seed
  {
    set.seed(seed.no)  
  } 
  
  #Gibbs Sampling
  for (b in 1:(B-1)){
    print(b)
    # Step1: Sample theta1 
    {
      sigma.theta1.vec <- solve((1/global.sigma[b])*H(theta2.vec[,b], theta3.vec[,b], xi.vec[,b]) + 
                                  ((1/local.sigma1[b])*I.N))
      mu.theta1.vec <- sigma.theta1.vec%*%((1/global.sigma[b])*r.vec(theta2.vec[,b], theta3.vec[,b], xi.vec[,b])+
                                             (1/local.sigma1[b])*(onevec.N*alpha1[b] + X%*%beta1.vec[,b]))
      theta1.vec[,b+1] <- mvtnorm::rmvnorm(n=1, mean = mu.theta1.vec, 
                                           sigma = sigma.theta1.vec)
    }
    # Step2: Sample theta2 & theta3 using MH with Gaussian proposal densities
    {
      ##M-H algorithm for theta2
      for(i in 1:N){
        #Step1: generate proposals
        new.theta2.i <- rnorm(n = 1, mean = theta2.vec[i,b], sd = sqrt(pro.var.theta.2))
        
        #Step2: Calculate acceptance probability
        theta2.prob <- theta2.ratio(i = i, theta2.i = new.theta2.i, 
                                    old.theta2.i = theta2.vec[i,b],
                                    theta1.i = theta1.vec[i,(b+1)],
                                    theta3.i = theta3.vec[i,b],
                                    xi.i = xi.vec[i,b],
                                    global.sigma = global.sigma[b],
                                    alpha2 = alpha2[b],
                                    local.sigma2 = local.sigma2[b], 
                                    beta2.vec = beta2.vec[,b])
        #Step3: Reject or Accept
        u <- runif(1)
        theta2.vec[i,(b+1)] <- ifelse(u <= theta2.prob, new.theta2.i, theta2.vec[i,b])
      }

      ##M-H algorithm for theta3
      for(i in 1:N){
        #Step1: generate proposals
        new.theta3.i <- rnorm(n = 1, mean = theta3.vec[i,b], sd = sqrt(pro.var.theta.2))
        
        #Step2: Calculate acceptance probability
        theta3.prob <- theta3.ratio(i = i,  theta1.i = theta1.vec[i,(b+1)],
                                    theta2.i = theta2.vec[i,(b+1)],
                                    old.theta3.i = theta3.vec[i,b],
                                    theta3.i = new.theta3.i,
                                    xi.i = xi.vec[i,b],
                                    global.sigma = global.sigma[b],
                                    alpha3 = alpha3[b],
                                    local.sigma3 = local.sigma3[b], 
                                    beta3.vec = beta3.vec[,b])
        #Step3: Reject or Accept
        u <- runif(1)
        theta3.vec[i,(b+1)] <- ifelse(u <= theta3.prob, new.theta3.i, theta3.vec[i,b])
      }
      
    } 
    #Step3: Sample xi (Elliptical Slice Sampling)
    {
      ##ESS: delicate consideration and Gaussian prior assumed
      for(i in 1:N) {
        # Step 1:  change variable
        old.eta.i <- log(xi.vec[i,b])
        #(ESS) Step 2-a: cho-ose an ellipse
        nu <- rnorm(1, mean = mu, sd = sqrt(rho.sq))
        #(ESS) Step 2-b: define a criterion function
        #(ESS) Step 2-c: choose a threshold
        u <- runif(1)
        #(ESS) Step 2-d: draw an initial proposal
        phi <- runif(1, min= -pi, max=pi)
        eta.star.i <- (old.eta.i-mu)*cos(phi)+(nu-mu)*sin(phi)+mu #assume N(0,sigma)
        #(ESS) Step 2-e: accept procedure
        if(u < eta.ratio(i = i, eta.i = eta.star.i, old.eta.i = old.eta.i, theta1.i = theta1.vec[i,(b+1)], 
                         theta2.i = theta2.vec[i,(b+1)],theta3.i = theta3.vec[i,(b+1)], global.sigma = global.sigma[b])) {
          new.eta.i <- eta.star.i
        } else {
          phi.min = -pi ; phi.max = pi 
          while(u >= eta.ratio(i = i, eta.i = eta.star.i, old.eta.i = old.eta.i, theta1.i = theta1.vec[i,(b+1)], 
                               theta2.i = theta2.vec[i,(b+1)],theta3.i = theta3.vec[i,(b+1)], global.sigma = global.sigma[b])){
            if(phi > 0){
              phi.max = phi
            } else {
              phi.min = phi
            }
            phi <- runif(1, min=phi.min, max=phi.max)
            eta.star.i <- (old.eta.i-mu)*cos(phi)+(nu-mu)*sin(phi)+mu 
          }
          new.eta.i <- eta.star.i
        }
        #(ESS)  step3: variable change
        xi.vec[i, (b+1)] = exp(new.eta.i)
      }
    }
    #Step4: Sample global.sigma 
    {
      gs.shape <- (N*Ts)/2 + varrho
      temp.gs <- numeric(N)
      for(i in 1:N){
        #qq
        temp.gs[i] <-norm.f(x=y(i), nul=f.i.vec(i = i, theta1.i = theta1.vec[i,(b+1)], theta2.i = theta2.vec[i,(b+1)], theta3.i = theta3.vec[i,(b+1)], xi.i = xi.vec[i,(b+1)]))^2
      }
      rate <- (1/2)*sum(temp.gs) + varrho
      global.sigma[b+1] <- 1/rgamma(1, shape = gs.shape, rate=temp.gs)
    }
    #Step5: Sample alpha1,2,3
    {
      #alpha1
      mean.alpha1 <- (1/N)*t(onevec.N)%*%(theta1.vec[,(b+1)]-X%*%beta1.vec[,b])
      sd.alpha1 <- sqrt(local.sigma1[b]/N)
      alpha1[b+1] <- rnorm(1, mean=mean.alpha1, sd=sd.alpha1)
      
      #alpha2 
      mean.alpha2 <- (1/N)*t(onevec.N)%*%(theta2.vec[,(b+1)]-X%*%beta2.vec[,b])
      sd.alpha2 <- sqrt(local.sigma2[b]/N)
      alpha2[b+1] <- rnorm(1, mean=mean.alpha2, sd=sd.alpha2)
      
      #alpha3
      mean.alpha3 <- (1/N)*t(onevec.N)%*%(theta3.vec[,(b+1)]-X%*%beta3.vec[,b])
      sd.alpha3 <- sqrt(local.sigma3[b]/N)
      alpha3[b+1] <- rnorm(1, mean=mean.alpha3, sd=sd.alpha3)  
      
    }
    #Step6-1: Sample beta1,2,3
    {
      #beta1 
      ## with precision matrix
      Q1  <- (1/local.sigma1[b])*(t(X)%*%X + diag(1/(tau1[b]*lambda1.vec[,b])^2))
      b1 <- (1/local.sigma1[b])*(t(X)%*%(theta1.vec[,(b+1)]-onevec.N*alpha1[b+1]))
      beta1.vec[,(b+1)] = beta.rmvt(Q = Q1, b = b1)
      
      #beta2
      Q2  <- (1/local.sigma2[b])*(t(X)%*%X + diag(1/(tau2[b]*lambda2.vec[,b])^2))
      b2 <- (1/local.sigma2[b])*(t(X)%*%(theta2.vec[,(b+1)]-onevec.N*alpha2[b+1]))
      beta2.vec[,(b+1)] = beta.rmvt(Q = Q2, b = b2)
      
      #beta3
      Q3  <- (1/local.sigma3[b])*(t(X)%*%X + diag(1/(tau3[b]*lambda3.vec[,b])^2))
      b3 <- (1/local.sigma3[b])*(t(X)%*%(theta3.vec[,(b+1)]-onevec.N*alpha3[b+1]))
      beta3.vec[,(b+1)] = beta.rmvt(Q = Q3, b = b3)
    }

    #Step6-2: Sample local.sigma1,2,3
    {
      #local.sigma1
      ls1.shape <- (N+p)/2 + varrho
      ls1.f <- norm.f(x = theta1.vec[,(b+1)], nul = onevec.N*alpha1[b+1]+X%*%beta1.vec[,(b+1)])
      ls1.s <- t(beta1.vec[,(b+1)])%*%diag(1/(tau1[b]*lambda1.vec[,b])^2)%*%beta1.vec[,(b+1)]
      ls1.rate <- (1/2)*(ls1.f^2 + ls1.s) + varrho
      local.sigma1[b+1] = 1/rgamma(n = 1, shape = ls1.shape, rate = ls1.rate) 
      
      #local.sigma2
      ls2.shape <- (N+p)/2 + varrho
      ls2.f <- norm.f(x = theta2.vec[,(b+1)], nul = onevec.N*alpha2[b+1]+X%*%beta2.vec[,(b+1)])
      ls2.s <- t(beta2.vec[,(b+1)])%*%diag(1/(tau2[b]*lambda2.vec[,b])^2)%*%beta2.vec[,(b+1)]
      ls2.rate <- (1/2)*(ls2.f^2 + ls2.s) + varrho
      local.sigma2[b+1] = 1/rgamma(n = 1, shape = ls2.shape, rate = ls2.rate)
      
      #local.sigma3
      ls3.shape <- (N+p)/2 + varrho
      ls3.f <- norm.f(x = theta3.vec[,(b+1)], nul = onevec.N*alpha3[b+1]+X%*%beta3.vec[,(b+1)])
      ls3.s <- t(beta3.vec[,(b+1)])%*%diag(1/(tau3[b]*lambda3.vec[,b])^2)%*%beta3.vec[,(b+1)]
      ls3.rate <- (1/2)*(ls3.f^2 + ls3.s) + varrho 
      local.sigma3[b+1] = 1/rgamma(n = 1, shape = ls3.shape, rate = ls3.rate)
    }

    #Step7: Sample lambda1,2,3 with Slice sampling
    {
      #lambda1
      eta1.vec = lambda1.vec[,b]^2
      updated.eta1.vec = c()
      for (j in 1:p){
        updated.eta1.vec[j] = ss.invgam(old.x = eta1.vec[j],f = 1,
                                        s = beta1.vec[j,(b+1)]^2/(2*local.sigma1[b+1]*tau1[b]^2))
      }
      lambda1.vec[,(b+1)] = sqrt(updated.eta1.vec)
      
      
      #lambda2
      eta2.vec = lambda2.vec[,b]^2
      updated.eta2.vec = c()
      for (j in 1:p){
        updated.eta2.vec[j] = ss.invgam(old.x = eta2.vec[j],f= 1,
                                        s = beta2.vec[j,(b+1)]^2/(2*local.sigma2[b+1]*tau2[b]^2))
      }
      lambda2.vec[,(b+1)] = sqrt(updated.eta2.vec)
      
      #lambda3
      eta3.vec = lambda3.vec[,b]^2
      updated.eta3.vec = c()
      for (j in 1:p){
        updated.eta3.vec[j] = ss.invgam(old.x = eta3.vec[j],f = 1,
                                        s = beta3.vec[j,(b+1)]^2/(2*local.sigma3[b+1]*tau3[b]^2))
      }
      lambda3.vec[,(b+1)] = sqrt(updated.eta3.vec)
    }

    #Step8: Sample tau1,2,3
    {
      #tau1
      omega1 = tau1[b]^2
      f <-  (p+1)/2
      s <- (1/(2*local.sigma1[b+1]))*t(beta1.vec[,(b+1)])%*%diag(1/lambda1.vec[,(b+1)]^2)%*%beta1.vec[,(b+1)]
      updated.omega1 <- ss.invgam(old.x = omega1, f=f, s=s)
      tau1[b+1] = sqrt(updated.omega1)
      
      #tau2
      omega2 = tau2[b]^2
      f <-  (p+1)/2
      s <- (1/(2*local.sigma2[b+1]))*t(beta2.vec[,(b+1)])%*%diag(1/lambda2.vec[,(b+1)]^2)%*%beta2.vec[,(b+1)]
      updated.omega2 <- ss.invgam(old.x = omega2, f=f, s=s)
      tau2[b+1] = sqrt(updated.omega2)
      
      #tau3
      omega3 = tau3[b]^2
      f <-  (p+1)/2
      s <- (1/(2*local.sigma3[b+1]))*t(beta3.vec[,(b+1)])%*%diag(1/lambda3.vec[,(b+1)]^2)%*%beta3.vec[,(b+1)]
      updated.omega3 <- ss.invgam(old.x = omega3, f=f, s=s)
      tau3[b+1] = sqrt(updated.omega3)
      
    }

  }
  remainder = b%%1000 
  if(remainder == 0){
    mc.index = seq(from = burn + 1, to = burn + nmc, by = thin)
  }
  
  #Result
  {
    mc.index = seq(from = burn + 1, to = burn + nmc, by = thin)
    thinned.theta1.vec = theta1.vec[,mc.index]
    thinned.theta2.vec = theta2.vec[,mc.index]
    thinned.theta3.vec = theta3.vec[,mc.index]
    thinned.xi.vec = xi.vec[,mc.index]
    thinned.alpha1 = alpha1[mc.index]
    thinned.alpha2 = alpha2[mc.index]
    thinned.alpha3 = alpha3[mc.index]
    thinned.beta1.vec = beta1.vec[,mc.index]
    thinned.beta2.vec = beta2.vec[,mc.index]
    thinned.beta3.vec = beta3.vec[,mc.index]
    thinned.global.sigma = global.sigma[mc.index]
    thinned.local.sigma1 = local.sigma1[mc.index]
    thinned.local.sigma2 = local.sigma2[mc.index]
    thinned.local.sigma3 = local.sigma3[mc.index]
    
    val <- list(thinned.theta1.vec = thinned.theta1.vec,
               thinned.theta2.vec = thinned.theta2.vec,
               thinned.theta3.vec = thinned.theta3.vec,
               thinned.xi.vec = thinned.xi.vec,
               thinned.alpha1 = thinned.alpha1,
               thinned.alpha2 = thinned.alpha2,
               thinned.alpha3 = thinned.alpha3,
               thinned.beta1.vec = thinned.beta1.vec,
               thinned.beta2.vec = thinned.beta2.vec,
               thinned.beta3.vec = thinned.beta3.vec, 
               thinned.global.sigma = thinned.global.sigma,
               thinned.local.sigma1 = thinned.local.sigma1,
               thinned.local.sigma2 = thinned.local.sigma2,
               thinned.local.sigma3 = thinned.local.sigma3,
               mu = mu,
               rho.sq = rho.sq)
  } 
  return(val)
}
