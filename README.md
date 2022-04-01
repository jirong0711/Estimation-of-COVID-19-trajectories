<h1> Estimation of COVID-19 trajectories </h1>
<p align="left">
  <img src="https://img.shields.io/badge/R-276DC3?style=flat-square&logo=R&logoColor=white"/></a>&nbsp
  <img src="https://img.shields.io/badge/Jags-007396?style=flat-square"/></a>&nbsp 
  <img src="https://img.shields.io/badge/Stan-ffb13b?style=flat-square"/></a>&nbsp 
</p>


<b><i[Estimation of COVID-19 trajectories </i></b> is to estimate COVID-19 cumulative trajectories using data integration. This project follows a paper <b><i>'Estimation of COVID-19 spread curves integrating global data and borrowing information'</i></b>(https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0236860) written by Se Yoon Lee, Bowen Lei and Bani Mallick. This paper deals with efficiency and excellences of using data integration in estimating cumulatives of pandemic, that is COVID-19. 



<h2> Goals </h2>

- Understand overall structure of Bayesian model
- Comprehend property of shirnkage priors, especially **Horseshoe prior**
- Learn a number of sampling methods to estimate posterior distributions of BHRM parameters, especially **Slice Sampling and Elliptical Slice Sampling (ESS) **
- Implement Slice Sampling and ESS in **R and matlab**
- Fit model in **R, Stan and JAGs**
- Visualize fitted values of cumulative trajectories in R
- Predict expected numbers of infections in each country 


<h2> Files </h2>

### Data
- **design_matrix**: design matrix ('X') for BHRM model consisting of information of countries (e.g. health, age, population). 
- **time_series_data**: Cumulative trajectories, a response value for BHRM

### ESS
- **Slice Sampling Paper**: Original paper of Slice Sampling
- **Slice Sampling Slides and MSMC**: Referencs for learning base algorithm of SS
- **Slice Sampling and ESS.R**: Implementation of SS and ESS in R
- **ESS.matlab**: Implementation of ESS in matlab


### Functions
- **covid_curve**: a function of BHRM curves to estimate trajectories
- **flat_time_point**: a funciton of predicting trajectories in given time interval (gamma)
- **visualization**: a function of drawing BHRM curves

### Stan & JAGs
- **Stan**: Implemetation of covid_curve in terms of Stan
- **JAGS**: Implemetation of covid_curve in terms of JAGS

### Outputs
- **Elliptical Slice Sampling.pdf**: Organize concepts of SS and ESS and implementations in pdf format.
- **Estimation of covid curve.rmd**: Perform fitting models, predict and visualize outputs.
- **Estimation.pdf**: Organize the steps of constructing BHRM model.
- **Final.pdf**: Final Report of Project
- **Literature Review**: A literature review of main paper, Estimation of COVID-19 spread curves integrating global data and borrowing information
