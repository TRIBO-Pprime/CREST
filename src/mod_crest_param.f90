!< author: Arthur Francisco
!  version: 1.0.0
!  date: october, 23 2024
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.2em;">
!        **Global variables and types**
!  </span>

module crest_param
use data_arch,   only : I4, R8
use stat_mom,    only : moment_stat
use pikaia_oop,  only : pikaia_class
use surfile,     only : SCALE_SURF


implicit none

public

   integer(kind=I4) :: JOB                                     !! *JOB file: script*
   integer(kind=I4) :: SPY                                     !! *SPY file*
   integer(kind=I4) :: STA                                     !! *result file*

   integer(kind=I4), parameter :: TER =  6                     !! *terminal output*

   integer(kind=I4) :: line_read, save_line_read               !! *line number in the script"

   integer(kind=I4) :: i_iter, nb_iter                         !! *current and number of iterations*

   type(SCALE_SURF) :: scale_img                               !! *.SUR surface properties*

   type acf_param

      real(kind=R8) :: cl1                                     !! *correlation principal length at z=acf__z*
      real(kind=R8) :: cl2                                     !! *correlation secondary length at z=acf__z*
      real(kind=R8) :: cut                                     !! *acf cutting plane z, for correlation lengths determination*
      real(kind=R8) :: ang                                     !! *roughness orientation*

   endtype acf_param


   type param_crest

      real(kind=R8), allocatable, dimension(:)   :: vect_h     !! *vector used to store the heights that meet the stat moments*

      real(kind=R8), allocatable, dimension(:,:) :: surf       !! *surface array*
      real(kind=R8), allocatable, dimension(:,:) :: surf_copy  !! *surface array copy*
      real(kind=R8), allocatable, dimension(:,:) :: surf_LF    !! *surface low frequencies*
      real(kind=R8), allocatable, dimension(:,:) :: surf_HF    !! *surface high frequencies*

      real(kind=R8), allocatable, dimension(:,:) :: imp_acf    !! *imposed autocorrelation*
      real(kind=R8), allocatable, dimension(:,:) :: fhi        !! *digital filter*

      real(kind=R8), allocatable, dimension(:,:) :: acf_surf   !! *calculated autocorrelation*

      integer(kind=I4), allocatable, dimension(:) :: order     !! *vector that stores heights order*

      logical(kind=I4), allocatable, dimension(:,:) :: surf_msk!! *surface high frequencies*

      type(MOMENT_STAT) :: m_ini                               !! *surface to be reproduced stat moments*
      type(MOMENT_STAT) :: m__LF                               !! *surface to be reproduced stat moments, low  frequencies*
      type(MOMENT_STAT) :: m__HF                               !! *surface to be reproduced stat moments, high frequencies*
      type(MOMENT_STAT) :: m_end                               !! *final stat moments*
      type(MOMENT_STAT) :: m_inp                               !! *input stat moments for genetic algo optimizer*
      type(MOMENT_STAT) :: m_stt                               !! *starting stat moments*

      type(pikaia_class) :: pik_class                          !! **PIKAIA** *class instanciation*

      type(ACF_PARAM)  :: orig_surf                            !! *original surface ACF properties*
      type(ACF_PARAM)  :: curr_surf                            !! *current  surface ACF properties*

      integer(kind=I4) :: func_gen                             !! *mathematical function used to generate the heights*
      integer(kind=I4) :: nparam                               !! *number of parameters for the mathematical function*

      integer(kind=I4) :: nb_threads                           !! *number of concurrent threads*

      integer(kind=I4) :: width                                !! *surface nb points along x*
      integer(kind=I4) :: height                               !! *surface nb points along y*
      integer(kind=I4) :: npts                                 !! *surface nb points*
      integer(kind=I4) :: sub_width                            !! *subsurface nb points along x*
      integer(kind=I4) :: sub_height                           !! *subsurface nb points along y*
      integer(kind=I4) :: sub_npts                             !! *subsurface nb points*

      logical(kind=I4) :: reajust_skku                         !! *if Ssk**2 +1 > Sku, modify Sku*

      logical(kind=I4) :: periodic                             !! *is the surface periodic?*

      logical(kind=I4) :: apod                                 !! *apodize imposed acf?*

      logical(kind=I4) :: calc_mstt                            !! *calculate starting moments?*

      logical(kind=I4) :: calc_zf                              !! *calculate final heights?*

      real(kind=R8) :: cutoff                                  !! *Gaussian filter cutoff*

      real(kind=R8) :: surf_width                              !! *surface width (m)*
      real(kind=R8) :: surf_height                             !! *surface height (m)*

      real(kind=R8) :: sub_surf_width                          !! *subsurface width (m)*
      real(kind=R8) :: sub_surf_height                         !! *subsurface height (m)*

      real(kind=R8) :: surf_dx                                 !! *surface increment along x (m)*
      real(kind=R8) :: surf_dy                                 !! *surface increment along y (m)*

      real(kind=R8) :: crt_acf                                 !! *acf criterion: mean absolute difference between imposed and calculated acf allowed*

      real(kind=R8) :: res_acf                                 !! *store mean absolute difference between imposed and calculated acf*

   endtype param_crest

   type(param_crest) :: PARAM

   integer(kind=I4), parameter :: FCT_TANG = 1                 !! *tangent function for height generation*
   integer(kind=I4), parameter :: FCT_EXPO = 2                 !! *exponential function for height generation*


endmodule crest_param
