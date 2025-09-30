STA_SCRI
========
*
 **READ_PRF**                               Read surface
   "sur/tooth.sur"
*
 **NB_PROCS**
   -1                                       0 (1 thread), -1 (max thread)
*
 **REPR_IMG**
   +0                                       Surface reproduction step
*
 **CALC_FFH**                               Digital filter $fhi = \sqrt{ \lvert FFT(imp\_acv) \rvert }$
   .true.
*
 **MAKE_TEX**                               Initialize with a smooth surface
   "circle"
   +4
*
 **STA_LOOP**
      030                                   Loop max number
      *DIGI_FIL*                               Digital filter
      *CALC_ORD*                               Calc height order
      *CALC_Z_F*                               Apply order to the heights
         .false.
      *CALC_ACF*                               Calculate the height *Acf*
         .false.
      *PLT__ACF*                               ACF actions
      +0.1                                     Criterion: stop loop if threshold reached
      "---"                                    Print or not *Acf* along principal axes
 **END_LOOP**
*
 **CALC_ACF**                               Calculate the height *Acf*
   .false.
*
 **PLT__ACF**
   -03.                                     Negative if no action required
   "xy-"                                    Print or not *Acf* along principal axes
*
 **REPR_IMG**
   +4                                       Surface reproduction step
*
 **SAVE_PRF**                               Save surface
   "out/tooth_reproduced.sur"
*
END_SCRI
========
