# Future yield

The scripts in this repository apply multivariable fractional polynomial regression (MFP) to model the nonlinear relationship between yield and climate for 1173 counties in 16 states in the central United States.  Response curves calibrated on historical data are applied to climate projections to estimate yields for corn, soy, and winter wheat to the end of the 21st century.

The `FY_do.R` script contains the entire workflow needed to reproduce the analyses in [PAPER].  The `FY_func.R` script contains all functions used in this workflow.  We have also included a script that loads essential data (`FY_load.R`) and a script that documents any cleaning that was done to the data (`FY_clean.R`).  A summary of our data and results can be found in the `mfp_results.html` document.  

