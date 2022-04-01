ESS <- function(initial_theta, prior, lnpdf,
                angle_range = 0) {
  library(MASS)
  D <- length(initial_theta)
  cur_lnpdf <- lnpdf(D)
  
  # Set up the ellipse and the slice threshold
  if(is.null(dim(prior))) nu <- prior
  else{
    if((dim(prior)[1] == D) | (dim(prior)[2] == D)){
      return(print("Prior must be given by a D-element sample or DxD chol(Sigma)"))
    }
    nu <- mvrnorm(1, 0, prior)
  }
  hh <- log(runif(1)) + cur_lnpdf
  # Set up a bracket of angles and pick a first proposal.
  # "phi = (theta'-theta)" is a change in angle.
  if(angle_range == 0){
    # Bracket whole ellipse with both edges at first proposed point
    phi = runif(1, min=0, max=2*pi)
    phi_min = phi-2*pi
    phi_max = phi
  }
  else{
    phi_min = angle_range*runif(1)
    phi_max = phi_min + angle_range
    phi = runif(1, min = phi_min, max = phi_max)
  }
  
  # Slice sampling loop
  while(TRUE){
    # Compute xx for proposed angle difference and check if it's on the slice
    xx_prop = initial_theta*cos(phi)+nu*sin(phi)
    cur_lnpdf = lnpdf(xx_prop)
    # New point is on slice, ** EXIT LOOP **
    if(cur_lnpdf > hh) {break}
    # Shrink slice to rejected point
    if(phi > 0) {phi_max = phi}
    else if(phi < 0) {phi_min = phi}
    else{return(print("BUG DETECTED: Shrunk to current position and still not acceptable."))}
    # Propose new angle difference
    phi = runif(1, min=phi_min, max=phi_max) 
  }
  return(list(xx=xx_prop, ln=cur_lnpdf))
}

