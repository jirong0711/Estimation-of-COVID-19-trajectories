#Concept: Slice Sampling - Exponential


slice <- function(n, theta_init, target) {
  u = theta = numeric(n)
  theta[1] = theta_init
  u[1] = runif(1, 0, target(theta[1]))
  for (i in 2:n){
    u[i] = runif(1,0, target(theta[i-1]))
    theta[i] = runif(1, 0, -log(u[i]))
  }
  return(list(theta=theta,u=u))
}
res = slice(1000000, 0.1)
hist(res$theta, freq=F, breaks=50)



#Basic slice sampler
SliceSampler = function(x, pstar, w, n.sims) {
  
  in.slice = function(x, y) { y <= pstar(x) }
  sims = matrix(NA, n.sims, 2)
  
  for (s in 1:n.sims) {
    y = runif(1, 0, pstar(x)) 
    l = x - w * runif(1)
    u = l + w
    
    l.in = in.slice(l, y)
    if (l.in) {
      while (l.in) {
        l = l - w
        if (l < 0) { 
          l = 0
          break
        }
        l.in = in.slice(l, y)
      }
    }
    u.in = in.slice(u, y)
    if (u.in) {
      while (u.in) {
        u = u + w
        if (u > 1) { 
          u = 1
          break
        }
        u.in = in.slice(u, y)
      }
    }
    x.old = x 
    x.in = FALSE
    while (!x.in) {
      x = runif(1, l, u)
      x.in = in.slice(x, y)
      if (x > x.old) {
        u = x
      } else {
        l = x
      }
    }
    sims[s, 1] = x
    sims[s, 2] = y
  }
  return(sims)
}

#Implementation
n.sims = 30000
x.init = .5
mix.beta = function(x) {
  .45 * dbeta(x, 2, 10) + .45 * dbeta(x, 10, 2) + .1 * dbeta(x, 3, 3)
}
w = .2
slice.samps.mix = SliceSampler(x.init, mix.beta, w, n.sims)
