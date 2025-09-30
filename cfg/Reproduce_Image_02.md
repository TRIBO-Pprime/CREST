STA_SCRI
========
*
 **READ_PRF**
   "sur/tooth.sur"
*
 **NB_PROCS**
   -1                                       0 (1 thread), -1 (max thread)
*
================================================== step 1
*
**REPR_IMG**
   +1                                       Surface reproduction step
   -0.01                                    Cutoff of the Gaussian filter: lag/cutoff wavelength --- minus sign for top hat
*
 **CALC_FFH**                               Digital filter $fhi = \sqrt{ \lvert FFT(imp\_acv) \rvert }$
   .true.
*
 **CALC_Z_I**                               Starting heights stored in *prof_x*
*
 **STA_LOOP**
      030                                   Loop max number
      *SMOOTH__*
         -0.0100                               Cutoff of the Gaussian filter: lag/cutoff wavelength --- minus sign for top hat
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
 **SAVE_PRF**
   "out/tooth_LF.sur"
*
================================================== step 2
*
 **REPR_IMG**
   +2                                       Surface reproduction step
*
 **CALC_FFH**                               Digital filter $fhi = \sqrt{ \lvert FFT(imp\_acv) \rvert }$
   .true.
*
 **CALC_Z_I**                               Starting heights stored in *prof_x*
*
 **STA_LOOP**
      030                                   Loop max number
      *SMOOTH__*
         -0.0100                               Cutoff of the Gaussian filter: lag/cutoff wavelength --- minus sign for top hat
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
 **SAVE_PRF**
   "out/tooth_HF.sur"
*
================================================== step 3
*
 **REPR_IMG**
   +3                                       Surface reproduction step
*
 **CALC_FFH**                               Digital filter $fhi = \sqrt{ \lvert FFT(imp\_acv) \rvert }$
   .false.
*
 **STA_LOOP**
      010                                   Loop max number
      *SMOOTH__*
         -0.100                                Cutoff of the Gaussian filter: lag/cutoff wavelength --- minus sign for top hat
      *DIGI_FIL*                               Digital filter
      *CALC_ORD*                               Calc height order
      *CALC_Z_F*                               Apply order to the heights
         .false.
      *CALC_ACF*                               Calculate the height *Acf*
         .false.
      *PLT__ACF*                               ACF actions
         +0.1                                  Criterion: stop loop if threshold reached
         "---"                                 Print or not *Acf* along principal axes
 **END_LOOP**
*
 **CALC_ACF**                               Calculate the height *Acf*
   .false.
*
 **PLT__ACF**
   -03.                                     Negative if no action required
   "xy-"                                    Print or not *Acf* along principal axes
*
================================================== step 4
*
 **REPR_IMG**
   +4                                       Surface reproduction step
*
 **SAVE_PRF**
   "out/tooth_reproduced.sur"
*
END_SCRI
========
