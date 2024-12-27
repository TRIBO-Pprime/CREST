STA_SCRI
========

 **DEF_SIZE**
   1024           0512                      Number of points along $x$ and $y$ respectively
   200.e-6        100.e-6                   Width (m) and height (m)
   .False.                                  $x$, $y$ periodic surface

 **ACF_THEO**
   30.e-6         10.e-6         0.20       Correlation lengths $\tau1 \ge \tau2$ and cut plane $z$
   40.0                                     Roughness orientation $\alpha$
   .true.                                   *Acf* apodization

 **NB_PROCS**
   -1                                       0 (1 thread), -1 (max thread)

 **STA_THEO**
  -03.00          +15.00                    Skewness *Ssk* and kurtosis *Sku*

 **SAVE_ACF**
   "out/avv.sur"

 **CALC_FFH**                                   Digital filter $fhi = \sqrt{ \lvert FFT(imp\_acv) \rvert }$

 **CALC_Z_I**                                   Starting heights stored in *prof_x*

 **DIGI_FIL**                                   Apply digital filter $prof\_x = FFT^{-1}( FFT(prof\_x) * fhi )$

 **CALC_ORD**                                   Store order of *prof_x* in vector *order*

 **CALC_Z_F**                                   Heights with right stat moments and ordered

 **SMOOTH__**
   0.001                                    Cutoff of the Gaussian filter

 **STA_LOOP**
      050                                   Loop max number
      *DIGI_FIL*                              Digital filter
      *CALC_ORD*                              Calc height order
      *CALC_Z_F*                              Apply order to the heights
      *CALC_ACF*                              Calculate the height *Acf*
      *PLT__ACF*                              ACF actions
      +0.1                                  Criterion: stop loop if threshold reached
      "---"                                 Print or not *Acf* along principal axes
 **END_LOOP**

 **SUB_SURF**
   +400

 **CALC_ACF**                                   Calculate the height *Acf*

 **PLT__ACF**
   -03.                                     Negative if no action required
   "xy-"                                    Print or not *Acf* along principal axes

 **SAVE_PRF**
   "out/surf_final.sur"

END_SCRI
========
