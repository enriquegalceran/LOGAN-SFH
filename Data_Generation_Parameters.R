# Data_Generation_Parameters.R

"
    How to set parameters:
    Parameters are supposed to be one of the following:
      1) A vector of possible values
      2) A distribution of possible values
      3) A single value
    
    

"


Parameters <- list(
  name="snorm_burst",
  mfunc=massfunc_snorm_burst,
  mSFR=10,
  mpeak=seq(7,14.2,0.3),
  mperiod=seq(0.5,1.5, 0.2),
  mskew=seq(-0.5, 1, 0.1),
  mburstage=0.1,
  mburst=5,
  zfunc=Zfunc_massmap_box,
  Zstart=1e-4,
  yield=0.03,
  Zfinal=seq(0.02, 0.08, 0.02)
)