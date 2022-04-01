visualizations_individual <- function(i, res_cov, time_series_data, gamma = c(0.9999, 0.999, 0.99, 0.9)){
  # i: the index of the country to be displayed
  # res_cov: the output of BHRM_cov function
  # time_series data: N-by-T+2 dataframe storing the number of cumulative infected cases for each countries
  # gamma: a progression constant
  # Data import
  pred.T.days = 10
  par(mfrow = c(1,1))
  {
    
    # Import time series data
    Y = as.matrix(time_series_data[,-c(1,2)]) 
    
    # i = index for country ; j = index for covariate
    y = function(i){
      y = as.numeric(Y[i,]) 
      return(y)
    }
    
    # covariate model
    {
      
      thinned.theta1.vec_cov = res_cov$thinned.theta1.vec
      thinned.theta2.vec_cov = res_cov$thinned.theta2.vec
      thinned.theta3.vec_cov = res_cov$thinned.theta3.vec
      thinned.xi.vec = res_cov$thinned.xi.vec
    }
    
  }
  
  # Flat time
  {

    flat_time_points = rep(0, length(gamma))
    S = ncol(thinned.theta1.vec_cov)
    for (e in 1:length(gamma)){
      temp = c()
      for (s in 1:S){
        temp[s] = flat_time_point(theta1 = thinned.theta1.vec_cov[i,s], 
                                  theta2 = thinned.theta2.vec_cov[i,s],
                                  theta3 = thinned.theta3.vec_cov[i,s],
                                  xi = thinned.xi.vec[i,s],
                                  gamma = gamma[e])
      }
      flat_time_points[e] = mean(temp)
    }
    
  }
  
  # cov_model
  {
    mean_flat_time = max(flat_time_points)
    S = ncol(thinned.theta1.vec_cov)
    t.values = seq(from = 1, to = mean_flat_time + 10, length.out = 2000)
    y.temp = matrix(rep(0, S*2000), nrow = S, ncol = 2000)
    f = function(t, theta1,theta2, theta3, xi){
      # (Original) Richard model
      # theta.1: K : final epidemic size : real number
      # theta.2: r : intrinsic growth rate : real number
      # theta.3: tau : disease turning point : real number
      # xi: shape parameter (measuring deviation from symmetric, xi = 1 => SYMMETRIC)
      res = theta1/( 1 + xi*exp(-theta2*(t - theta3)))^(1/xi)
      return(res)
    }
    
    
    for (t in 1:2000){
      for (s in 1:S){
        y.temp[s,t] = f(t = t.values[t], 
                        theta1 = thinned.theta1.vec_cov[i,s], 
                        theta2 = thinned.theta2.vec_cov[i,s], 
                        theta3 = thinned.theta3.vec_cov[i,s],
                        xi = thinned.xi.vec[i,s])  
      }
    }
    
    cred.interval = function(x){
      q.05 = function(x){
        res = quantile(x = x, probs = 0.05)
        return(res)
      }
      q.95 = function(x){
        res = quantile(x = x, probs = 0.95)
        return(res)
      } 
      res.05 = apply(x, MARGIN = 2, FUN = q.05)
      res.95 = apply(x, MARGIN = 2, FUN = q.95)
      res = rbind(res.05,res.95)
      return(res)
    }
    cred_cov = cred.interval(y.temp)
    post_mean_cov = colMeans(y.temp) 
    
  }
  
  # Dataframe
  {
    
    post_mean = c(post_mean_cov) 
    lower_bound = c(cred_cov[1,])
    upper_bound = c(cred_cov[2,])
    models = c(rep("cov", length(post_mean_cov)))
    t.values = c(t.values)
    temp.df = data.frame(cbind(t.values, post_mean, lower_bound, upper_bound)) 
    temp.df$models = c(rep("cov", length(post_mean_cov)))
    
    # https://janhove.github.io/reporting/2017/05/12/visualising-models-2
    library(ggplot2)
    df_cov = temp.df[temp.df$models == "cov",] 
    df_cov_polygon = as.data.frame(cbind(c(df_cov$t.values,rev(df_cov$t.values)),
                                         c(df_cov$lower_bound,rev(df_cov$upper_bound)))  )  
    names(df_cov_polygon) = c("t.values","bounds")
    df_obs = as.data.frame(cbind(c(1:length(y(i))),as.numeric(y(i))))  
    names(df_obs) = c("t.values","obs")
    
    # ggplot
    {
      library("ggplot2")
      library("scales") # https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
      library("ggthemes")
      
      g.res = ggplot() +
        #geom_polygon(data = df_cov_polygon, mapping = aes(x = t.values, y = bounds),fill = "pink", alpha = 0.5) + 
        geom_line(data = df_cov, mapping = aes(x = t.values, y = post_mean), col = "red", size = 1.2) + 
        geom_point(data = df_obs, mapping = aes(x = t.values, y = obs),size = 2) +
        xlab("Days") + ylab("Cumulative number of infected cases") +
        scale_y_continuous(trans = 'sqrt',labels=comma)+
        theme(
          axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          plot.title = element_text(size = 14, face = "bold"))+
        scale_x_continuous(breaks=seq(0,floor(mean_flat_time) + 10,50)) +
        theme_hc()+ scale_colour_hc() +
        geom_vline(xintercept = c(floor(flat_time_points)), linetype="dashed", 
                   color = "violet", size=1) +
        geom_vline(xintercept = c(length(y(i))), linetype="dotted", 
                   color = "grey", size=1)+
        geom_hline(yintercept = c(mean(thinned.theta1.vec_cov[i,])), linetype="dashed", 
                   color = "blue", size=1)+
        ggtitle(paste("Extrapolated infection trajectory for ",  as.character(time_series_data[i,2])))
      
      floor(flat_time_points)
      mean(thinned.theta1.vec_cov[i,])
    }
  }
  flat_time = as.Date("2020-1-22") + rev(floor(flat_time_points))  
  
  list(figure = g.res, flat_time = flat_time, epidemic_size = round(mean(thinned.theta1.vec_cov[i,])))
}


