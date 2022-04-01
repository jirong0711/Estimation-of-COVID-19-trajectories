##data
setwd("/Users/jihyunlee/Covid/Covid-week1/dataset/")
time_series_data <- read.csv("new_time_series_data.csv")
Y = as.matrix(time_series_data[,-1])
design_matrix <- read.csv("new_design.csv")
X = as.matrix(design_matrix[,-c(1:2)])

##MODEL
library(rstan)
model_string <- "
  data{
    int<lower=1> N;
    int<lower=1> Ts;
    matrix[N,Ts] Y;
    int<lower=1> p;
    matrix[N, p] X;
  }
  parameters{
    matrix[3, N] theta;
    vector<lower=0>[N] xi;
    
    vector[3] alpha;
    matrix[3,p] beta;
    real<lower=0> global_sigma2;
    vector<lower=0>[3] local_sigma2;
    vector<lower=0>[3] tau;
    matrix<lower=0>[3,p] lambda;
  }
  model{
    for(n in 1:N){
      for(t in 1:Ts){
        Y[n,t] ~ normal(theta[1,n]*((1+xi[n]*exp(-theta[2,n]*(t-theta[3,n])))^(-1/xi[n])),
                        global_sigma2);
      }
      for(i in 1:3){
        theta[i,n] ~ normal(alpha[i]+dot_product(X[n,], beta[i,]), local_sigma2[i]);
      }
      xi[n] ~ lognormal(0,1);


    }
    for(i in 1:3){
      tau[i] ~ cauchy(0,1);
      for(j in 1:p){
        beta[i,j] ~ normal(0, (tau[i]^2)*(lambda[i,j]^2)*local_sigma2[i]);
        lambda[i,j] ~ cauchy(0,1);
      }
        local_sigma2[i] ~ gamma(7.5, 1);
        alpha[i] ~ normal(0,100);

    }
    global_sigma2 ~ gamma(7.5, 1);
  }
"
X = as.matrix(design_matrix[,-c(1:2)])
Y = as.matrix(time_series_data[,-1])
N = dim(Y)[1]; Ts = dim(Y)[2]
p = dim(X)[2]

model <- stan(model_code = model_string,
               data=list(Y=Y, X=X, N=N, Ts=Ts, p=p),
                init = function(){
                  list(xi = runif(N, 0.01, 5),
                       lambda = matrix(runif(3*p, 0.01, 5), nrow=3),
                       global_sigma = runif(1, 0.01, 1),
                       local_sigma = runif(3, 0.01, 1),
                       alpha = runif(3, -10, 10),
                       tau = runif(3, 0.01, 10)
                       )
                },
               chains = 3, iter = 5000*2, seed = 1)
  
