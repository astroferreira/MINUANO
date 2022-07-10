MINUANO - Pipeline 

NIRCAM

1. Detector reduction *uncal -> *rate, *rateint. In parallel: 
 
    find uncals/ -name '*uncal.fits' | xargs -n1 -P24 -I {} strun STAGE1_params.asdf {}

    PERFORMANCE for Pointing 5 (90 uncals):

        real    19m50.942s
        user    297m27.575s
        sys     67m7.145s

2. IMAGE Reduction 2

  time find calibrated/ -name '*rate.fits' | xargs -n1 -P24 -I {} strun STAGE2_params.asdf {}

  real    6m45.502s
  user    95m26.210s
  sys     41m47.545s


3. 


real    28m43.633s
user    33m52.800s
sys     21m16.182s
