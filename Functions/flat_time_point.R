flat_time_point = function(theta1,theta2,theta3,xi,gamma = 0.99){
  nomi = log((1/xi) * ( (1/(gamma))^{xi} -1 ) )
  denom = theta2
  res = theta3 - nomi/denom
  return(res)
}

