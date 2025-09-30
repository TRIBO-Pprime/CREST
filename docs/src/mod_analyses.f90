!< author: Arthur Francisco
!  version: 1.0.0
!  date: june, 24 2025
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.2em;">
!        **Routines for surface analyses**
!  </span>

module analyses

use data_arch,       only : I4, R8
use miscellaneous,   only : get_unit
use filter,          only : median_filter, fft_filter
use abbott,          only : abbott_param
use morpho,          only : topology, calcul_normales, surf_area
use grad_curv,       only : peaks_and_pits_curvatures
use stat_mom,        only : moment_stat, calc_moments, calc_median
use asfc,            only : calcul_asfc_hermite, indice_fractal
use anisotropy,      only : correlation_parameters, multiple_anisotropy, ellipse_acf
use fftw3,           only : init_fftw3, fftw_plan_with_nthreads, tab_init_fftw3, end_fftw3, tab_end_fftw3, NB_THREADS_FFT, PAD_FFT, FFTW_MEASURE
use files,           only : make_path, path2vec, vec2path, filename, dir_separator, mkdir, dirname
use func_acf,        only : acf_wiener
use crest_param,     only : PARAM, JOB, SPY, STA, TER, LINE_READ, SCALE_IMG

implicit none

private

public :: surface_analysis

contains

   subroutine surface_analysis(app)
   !================================================================================================
   !< @note
   !<
   !< The function *analyses* determinates ISO 25178 parameters of the current surface.
   !<
   !< The analysis is always performed on the whole surface.
   !< The results are written in the file which unit is *STA*.
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in), optional :: app   !! *append results to csv*

      integer(kind=I4) :: i, ib, ie, nn, pp, append

      real(kind=R8)     :: ra_t, md
      real(kind=R8)     :: dx, dy, si, fft_cutoff

      type(MOMENT_STAT) :: mx

      real(kind=R8), dimension(:,:), allocatable :: tab_tmp

      real(kind=R8), dimension(:), allocatable :: vec_heights, tab_results
      real(kind=R8), dimension(1:20)           :: ana_res

      character(len= :), allocatable :: result_str, head

      character(len=18) :: str

      if ( .not.present(app) ) then

         read(JOB,*) append ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "append ", append

      else

         append = app

      endif

      ! result file
      call get_unit( STA )

      if (append == 1) then

         open(unit = STA, file = "out/res.csv", status = 'unknown', position = 'append')

      else

         open(unit = STA, file = "out/res.csv", status = 'unknown')

      endif

      !---------------------------------------------------------

      ana_res = 0

      nn = PARAM%width
      pp = PARAM%height

      dx = PARAM%surf_dx
      dy = PARAM%surf_dy
      si = 1

      allocate( tab_results(1:100) )

      allocate( vec_heights(1:nn * pp) )

      vec_heights(1:nn * pp) = reshape( PARAM%surf(1:nn, 1:pp), [ nn * pp ] )

      !-----------------------------------

      head = 'smrk1,smrk2,spk,svk,Sk,pente,residus,coeffa_tan,coeffb_tan,'

      call abbott_param( tab     = vec_heights(1:nn * pp),           &  !
                         lg      = nn * pp,                          &  !
                         nom     = "out/abbott_res.txt",             &  !
                         curves  = [.false., .false., .false.],      &  !
                         results = ana_res(1:11),                    &  !
                         omp     = .true. )                             !

      ib = 1
      ie = ib + 9 - 1
      tab_results(ib:ie) = [ ana_res( 1), &  ! smrk1, iso 25178
                             ana_res( 2), &  ! smrk2, iso 25178
                             ana_res( 3), &  ! spk  , iso 25178
                             ana_res( 4), &  ! svk  , iso 25178
                                             ! 5 et 6 pour off1 et off2
                             ana_res( 7), &  ! sk   , iso 25178
                             ana_res( 8), &  ! core slope
                             ana_res( 9), &  ! adjustment factor (tangent fit)
                             ana_res(10), &  ! coeffa_tan        (tangent fit)
                             ana_res(11) ]   ! coeffb_tan        (tangent fit)

      deallocate( vec_heights )

      !-----------------------------------

      head = head//'Snb1,Smc1,Sk1,Snb2,Smc2,Sk2,Sdq,Scq,Sh3z,Sv3z,'

      call topology( tab  = PARAM%surf(1:nn, 1:pp),    &  !
                     long = nn,                        &  !
                     larg = pp,                        &  !
                     res  = ana_res(1:6) )                !

      fft_cutoff = dx / 5.e-6 ! 5.e-6 = 5 µm

      allocate( tab_tmp(1:nn, 1:pp) )

      call fft_filter(tab       = PARAM%surf(1:nn, 1:pp),      &  ! in
                      long      = nn,                          &  ! in
                      larg      = pp,                          &  ! in
                      cutoff    = fft_cutoff,                  &  ! in
                      bf_tab    = tab_tmp(1:nn, 1:pp),         &  ! out
                      multi_fft = .false.)                        ! in

      call peaks_and_pits_curvatures( heights      = tab_tmp(1:nn, 1:pp),  &  !
                                      nx           = nn,                   &  !
                                      ny           = pp,                   &  !
                                      dx           = dx,                   &  !
                                      dy           = dy,                   &  !
                                      S_param_grad = ana_res(07),          &  !
                                      S_param_curv = ana_res(08),          &  !
                                      peak_curv    = ana_res(09),          &  !
                                      pits_curv    = ana_res(10) )            !

      deallocate( tab_tmp )

      ib = ie + 1
      ie = ib + 10 - 1

      tab_results(ib:ie) = ana_res(1:10)

      !-----------------------------------

      head = head//'Sv,Sp,Smd,Sa,Sm,Sq,Ssk,Sku,Sks,'

      call calc_moments( tab    = reshape( PARAM%surf(1:nn, 1:pp), [nn * pp] ),   &  !
                         mx     = mx,                                             &  !
                         nb_mom = 4 )                                                !

      call calc_median( tab = reshape( PARAM%surf(1:nn, 1:pp), [nn * pp] ),       &  !
                        md  = md )                                                   !

      ra_t = sum( abs(PARAM%surf(1:nn, 1:pp) - mx%mu ) / (nn * pp) )

      ana_res(1:8) = [ minval( PARAM%surf(1:nn, 1:pp) ) - mx%mu,   &  !
                       maxval( PARAM%surf(1:nn, 1:pp) ) - mx%mu,   &  !
                                                     md - mx%mu,   &  !
                            ra_t, mx%mu, mx%si, mx%sk, mx%ku ]        !

      ib = ie + 1
      ie = ib + 8 - 1

      tab_results(ib:ie) = ana_res(1:8)

      ib = ie + 1
      ie = ib

      tab_results(ib) = ana_res(8)/( ana_res(7)**2 + 1 )     ! kind of kurtosis excess

      !-----------------------------------

      head = head//'Smbd,ordorig,R2adj,'

      call indice_fractal( tab_in = PARAM%surf(1:nn, 1:pp),  &  !
                           long   = nn,                      &  !
                           larg   = pp,                      &  !
                           indf   = ana_res(1:3) )              !

      ib = ie + 1
      ie = ib + 3 - 1

      tab_results(ib:ie) = ana_res(1:3)

      !-----------------------------------

      head = head//'Sh,Sdr,'

      call calcul_normales( tab_in     = PARAM%surf(1:nn, 1:pp),   &  !
                            long       = nn,                       &  !
                            larg       = pp,                       &  !
                            scale_xyz  = [ dx, dy, si ],           &  !
                            cone_angle = 5._R8,                    &  !
                            hori       = ana_res(1) )                 !

      call surf_area( tab_in     = PARAM%surf(1:nn, 1:pp),         &  !
                      long       = nn,                             &  !
                      larg       = pp,                             &  !
                      scale_xyz  = [ dx, dy, si ],                 &  !
                      aire       = ana_res(2) )                       !



      ib = ie + 1
      ie = ib + 2 - 1

      tab_results(ib:ie) = ana_res(1:2)

      !-----------------------------------

      head = head//'Sasfc,R2adj,'

      call calcul_asfc_hermite( tab_in   = PARAM%surf(1:nn, 1:pp),    &  !
                                scal     = SCALE_IMG,                 &  !
                                asfc_res = ana_res(1:2),              &  !
                                omp      = .true. )                      !

      ib = ie + 1
      ie = ib + 2 - 1

      tab_results(ib:ie) = ana_res(1:2)

      !-----------------------------------

      head = head//'Rmax,Sal,Stri,Std,d.sl,b.sl,r.sl,r.cv,bmp,smp,rmp,bml,sml,rml,bms,sms,rms'

      if ( sum(PARAM%acf_surf(1:nn, 1:pp)) == 0 ) then

         call acf_wiener(  tab_in = PARAM%surf(1:nn, 1:pp),          &  ! IN
                          tab_out = PARAM%acf_surf(1:nn, 1:pp),      &  ! OUT
                                w = nn,                              &  ! IN
                                h = pp  )                               ! IN

      endif

      call ellipse_acf( tabin = PARAM%acf_surf(1:nn, 1:pp),          &  ! IN
                         long = nn,                                  &  ! IN
                         larg = pp,                                  &  ! IN
                        p_acv = ana_res(1:8),                        &  ! OUT -> correlation lengths
                          cut = PARAM%curr_surf%cut,                 &  ! IN  -> z cut plane
                     scale_xy = [PARAM%surf_dx, PARAM%surf_dy],      &  ! IN  -> lags along x and y
                          omp = .true. )                                ! IN  -> use multithread?

      PARAM%curr_surf%cl1 = ana_res(1)
      PARAM%curr_surf%cl2 = ana_res(2)
      PARAM%curr_surf%ang = ana_res(4)

      call multiple_anisotropy( tabin     = PARAM%surf(1:nn, 1:pp),  &  ! IN
                                long      = nn,                      &  ! IN
                                larg      = pp,                      &  ! IN
                                scale_xy  = [ dx, dy ],              &  ! IN
                                multi_fft = .false.,                 &  ! IN
                                vec_ani   = ana_res(9:17) )             ! OUT

      ib = ie + 1
      ie = ib + 17 - 1

      tab_results(ib:ie) = ana_res(1:17)

      if (append /= 1) write(STA,'(a)') trim(head)


      write(str, '(E18.6)') tab_results(1)
      result_str = str
      do i = 2, ie

         write(str, '(E18.6)') tab_results(i)

         result_str = result_str//','//str

      enddo

      write(STA, *) result_str

      deallocate( head )
      deallocate( result_str )
      deallocate( tab_results )

   return
   endsubroutine surface_analysis

endmodule analyses
