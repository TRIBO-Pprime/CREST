!< author: Arthur Francisco
!  version: 1.0.0
!  date: october, 23 2024
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.2em;">
!        **Routines for acf calculations**
!  </span>

module func_acf
use data_arch,     only : I4, R8, HIG_E8, EPS_R8, UN, PI_R8
use fftw3,         only : calc_fftw3_real_bwd, calc_fftw3_real_fwd, tab_calc_fftw3_real_bwd, tab_calc_fftw3_real_fwd, fftw_plan_with_nthreads, init_fftw3_real, end_fftw3, PAD_FFT, extend, &  !
                          SINGL_FFTW_ALLOCATED, NB_THREADS_FFT, FFT_DIM, FFTW_ESTIMATE, FFTW_MEASURE, FFTW_EXHAUSTIVE
use stat_mom,      only : calc_moments
use crest_param,   only : PARAM, SPY, TER
use stat_mom,      only : moment_stat
use miscellaneous, only : trans_corner2center, trans_center2corner

implicit none

private

public :: calc_imp_acf, acf_wiener, apod2

   contains

   subroutine acf_wiener(tab_in, tab_out, w, h, multi_fft)
   !================================================================================================
   !<@note Function that returns the *acf* of an array.
   !< \[
   !< \begin{align*}
   !<    acf(i,j) &= (z \ast z)(i,j) = \sum_{k,l}^{n,n} z(k+1-i,l+1-j)z(k,l)  \\
   !<    TF(acf)  &= ACF = Z \cdot Z                                          \\
   !<    acf      &= TF^{-1}(ACF) = TF^{-1}(Z^2)
   !< \end{align*}
   !< \]
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                      :: w           !! *2D array length*
   integer(kind=I4), intent(in )                      :: h           !! *2D array width*
   real   (kind=R8), intent(in ), dimension(1:w, 1:h) :: tab_in      !! *input array*
   real   (kind=R8), intent(out), dimension(1:w, 1:h) :: tab_out     !! *acf of the input array*
   logical(kind=I4), intent(in ), optional            :: multi_fft   !! *run parallel acfs?*

      integer(kind=I4) :: lo2, la2

      real   (kind=R8) :: tmp

      logical(kind=I4) :: parallel_fft

      integer(kind=I4), dimension(1:2) :: loc_max

      complex(kind=R8), dimension(:,:), allocatable :: tab_cmpl
      real   (kind=R8), dimension(:,:), allocatable :: tab_real

      allocate( tab_cmpl(1:w, 1:h) )

      allocate( tab_real(1:w, 1:h) )

      ! check for simultaneous fftw calculations
      !.........................................
      parallel_fft = .false.

      if ( present(multi_fft) ) parallel_fft = multi_fft
      !.........................................

      ! DFFT real -> complex
      !.........................................
      if ( parallel_fft ) then

         call tab_calc_fftw3_real_fwd( tab_in = tab_in  (1:w, 1:h),     &  ! IN
                                       tab_ou = tab_cmpl(1:w, 1:h),     &  ! OUT
                                         long = w,                      &  ! IN
                                         larg = h )                        ! IN

      else

         call calc_fftw3_real_fwd(      tab_in = tab_in  (1:w, 1:h),    &  ! IN
                                        tab_ou = tab_cmpl(1:w, 1:h),    &  ! OUT
                                          long = w,                     &  ! IN
                                          larg = h,                     &  ! IN
                                  planner_flag = FFTW_MEASURE )            ! IN

      endif
      !.........................................

      tab_cmpl(1:w, 1:h) = cmplx( abs( tab_cmpl(1:w, 1:h) )**2, 0, kind = R8 )


      ! IFFT complex -> real
      !.........................................
     if ( parallel_fft ) then

         call tab_calc_fftw3_real_bwd( tab_in = tab_cmpl(1:w, 1:h),     &  ! IN
                                       tab_ou = tab_real(1:w, 1:h),     &  ! OUT
                                         long = w,                      &  ! IN
                                         larg = h )                        ! IN
      else

         call calc_fftw3_real_bwd(      tab_in = tab_cmpl(1:w, 1:h),    &  ! IN
                                        tab_ou = tab_real(1:w, 1:h),    &  ! OUT
                                          long = w,                     &  ! IN
                                          larg = h,                     &  ! IN
                                  planner_flag = FFTW_MEASURE )            ! IN

      endif
      !.........................................

      ! the maximum is placed in the array center
      !.........................................
      call trans_corner2center(  tab_in  = tab_real(1:w, 1:h),  &  ! IN
                                 tab_out = tab_out (1:w, 1:h),  &  ! OUT
                                 long    = w,                   &  ! IN
                                 larg    = h  )                    ! IN
      !.........................................


      ! the maximum is 1
      !.........................................
      loc_max(1:2) = maxloc( tab_out(1:w, 1:h) )
      lo2 = loc_max(1)
      la2 = loc_max(2)

      tmp = tab_out(lo2, la2)

      tab_out(1:w, 1:h) = tab_out(1:w, 1:h) / tmp
      !.........................................

      deallocate(tab_cmpl)
      deallocate(tab_real)

   return
   endsubroutine acf_wiener


   real(kind=R8) function autocov_impo(xi, xj, tau1, tau2, alpha, ang)
   !================================================================================================
   !<@note Function that returns \( \exp \left(\alpha \sqrt{\left(\frac{x}{\tau_1}\right)^2+
   !<                                                        \left(\frac{y}{\tau_2}\right)^2}
   !<                                      \right) \)
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real(kind=R8), intent(in) :: tau1   !! *correlation length along \(x\)*
   real(kind=R8), intent(in) :: tau2   !! *correlation length along \(y\)*
   real(kind=R8), intent(in) :: alpha  !! *log(z)* where *z* is often 0.2
   real(kind=R8), intent(in) :: xi     !! *\(x\) coordinate*
   real(kind=R8), intent(in) :: xj     !! *\(y\) coordinate*
   real(kind=R8), intent(in) :: ang    !! *angle* (rad)

      real(kind=R8) :: x, y, r

      x = +cos(ang) * xi + sin(ang) * xj
      y = -sin(ang) * xi + cos(ang) * xj

      r = sqrt( (x / tau1)**2 + (y / tau2)**2 )

      autocov_impo = exp( alpha * r )

   return
   endfunction autocov_impo


   subroutine calc_imp_acf(long, larg, tau1, tau2, alpha, ang, tab_acf, apod)
   !================================================================================================
   !<@note Function that returns the theoretical autocorrelation function in an array.<br/>
   !< The autocorrelation function is supposed to be obtained from a real surface which must be periodic
   !< or nearly periodic (because of the use of FFTs).
   !< In addition, the surface is supposed to be 0 mean and normalized (\(\sigma = 1 \)),
   !< therefore *acf* is zero-mean and normalized so that its max value is 1.<br/>
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in) :: long   !! *surface acf width*
   integer(kind=I4), intent(in) :: larg   !! *surface acf height*
   logical(kind=I4), intent(in) :: apod   !! *apodization?*
   real   (kind=R8), intent(in) :: tau1   !! *first correlation length*
   real   (kind=R8), intent(in) :: tau2   !! *surface second correlation length*
   real   (kind=R8), intent(in) :: alpha  !! *parameter that controls the expondential decrease*
   real   (kind=R8), intent(in) :: ang    !! *acf ellipsis angle*
   real   (kind=R8), dimension(1:long, 1:larg), intent(out) :: tab_acf  !! *resulting acf*

      integer(kind=I4) :: i, j, long2, larg2
      real   (kind=R8) :: xi, xj, s, c, coeff

      real(kind=R8), dimension(:,:), allocatable:: tab_tmp

      ! acf array, centered and normalized
      !.........................................
      c = cos(ang) ; s = sin(ang)

      long2 = long / 2 +1
      larg2 = larg / 2 +1

      do j = 1, larg
      do i = 1, long

         xi = real(i - long2, kind=R8) * PARAM%surf_dx      ! dimensioned coordinate x
         xj = real(j - larg2, kind=R8) * PARAM%surf_dy      ! dimensioned coordinate y

         tab_acf(i, j) = autocov_impo( xi    = xi,       &  ! IN
                                       xj    = xj,       &  ! IN
                                       tau1  = tau1,     &  ! IN
                                       tau2  = tau2,     &  ! IN
                                       alpha = alpha,    &  ! IN
                                       ang   = ang )        ! IN

      enddo
      enddo
      !.........................................

      ! For long correlation lengths and roughness orientation, the acf is far from periodic
      ! Furthermore, far from the center, respecting the acf becomes less important. A windowing
      ! can be determined so that at a given distance from the center, the acf is lessened.
      !.........................................
      if ( apod ) then

         allocate( tab_tmp(1:long, 1:larg) )

         coeff = 0.4 * PARAM%surf_width * c / tau1

         ! along the primary axis (longest correlation length) the acf is reduce beyond
         ! 0.4 * image width * cos(ang)
         ! (0.4 * image width is less than half width)

         call apod2(  tab_in = tab_acf(1:long, 1:larg),     &  ! IN
                     tab_out = tab_tmp(1:long, 1:larg),     &  ! OUT
                        long = long,                        &  ! IN
                        larg = larg,                        &  ! IN
                        tau1 = coeff * tau1 ,               &  ! IN
                        tau2 = coeff * tau2 ,               &  ! IN
                         ang = ang )                           ! IN

         tab_acf(1:long, 1:larg) = tab_tmp(1:long, 1:larg)

         deallocate( tab_tmp )

      endif
      !.........................................

      ! acf centered
      tab_acf(1:long, 1:larg) = tab_acf(1:long, 1:larg) - sum( tab_acf(1:long, 1:larg) ) / (long * larg)

      ! acf scaled (maximum = 1)
      tab_acf(1:long, 1:larg) = tab_acf(1:long, 1:larg) / tab_acf(long2, larg2)

   return
   endsubroutine calc_imp_acf


   subroutine apod2(tab_in, tab_out, long, larg, tau1, tau2, ang)
   !================================================================================================
   !<@note Function that returns an apodized array.<br/>
   !< To prevent gaps from appearing after FFT (because of non periodic waves), the surface must
   !< be transformed, but not too much. Here a modified Tukey window is determined. The starting
   !< surface is not modified below the "correlation lengths". Above the correlation lengths, a
   !< smooth decrease is forced with a cosine squared.
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in )                            :: long     !! *surface acf length*
   integer(kind=I4), intent(in )                            :: larg     !! *surface acf width*
   real   (kind=R8), intent(in )                            :: tau1     !! *surface first correlation length*
   real   (kind=R8), intent(in )                            :: tau2     !! *surface second correlation length*
   real   (kind=R8), intent(in )                            :: ang      !! *ellipsis angle*
   real   (kind=R8), intent(in ), dimension(1:long, 1:larg) :: tab_in   !! *input acf*
   real   (kind=R8), intent(out), dimension(1:long, 1:larg) :: tab_out  !! *apodized acf*

      real   (kind=R8) :: r2, c0, s0, rd, rr, theta, theta_diag, x, y, t, a_min, sum_inn, sum_tab, sum_int
      integer(kind=I4) :: i, j, k, long2, larg2, npt_out, n

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      ! first bissector angle (rad)
      theta_diag = atan2( PARAM%surf_height, PARAM%surf_width )

      ! sine and cosine of the ellipsis angle
      c0 = cos(ang) ; s0 = sin(ang)

      long2 = long / 2
      larg2 = larg / 2

      if ( long == 2 * (long/2) ) long2 = long/2 + 1
      if ( larg == 2 * (larg/2) ) larg2 = larg/2 + 1

      tab_out(1:long, 1:larg) = tab_in(1:long, 1:larg)

      rr = 1.e6_R8

      do j = 1, larg
      do i = 1, long

         x = ( i - long2 ) * PARAM%surf_dx
         y = ( j - larg2 ) * PARAM%surf_dy

         ! Distance as expressed in the theoretical acf. When r2 is 1., the acf is 0.2.
         r2 = ( (+ c0 * x + s0 * y) / tau1 ) ** 2 +  &  !
              ( (- s0 * x + c0 * y) / tau2 ) ** 2

         ! Below the correlation length, no transformation
         if ( r2 <= 1._R8 ) cycle

         ! The correlation length is exceeded.
         ! The angle corresponding to the position (i,j) is calculated
         theta = atan2( y, x ) ; if ( theta < 0. ) theta = theta + 2 * Pi_R8

         t = tan(theta)

         ! According the location of (i,j) (right, top, left or bottom) the line that begin at the surface center, passing by (i,j),
         ! ends on one of the four borders.
         if ( theta > 2 * Pi_R8 - theta_diag .or.  theta <=           + theta_diag ) then ; x = ( long - long2 ) * PARAM%surf_dx ; y = x * t ; endif

         if ( theta >           + theta_diag .and. theta <=     PI_R8 - theta_diag ) then ; y = ( larg - larg2 ) * PARAM%surf_dy ; x = y / t ; endif

         if ( theta >     PI_R8 - theta_diag .and. theta <=     PI_R8 + theta_diag ) then ; x = (    1 - long2 ) * PARAM%surf_dx ; y = x * t ; endif

         if ( theta >     PI_R8 + theta_diag .and. theta <= 2 * Pi_R8 - theta_diag ) then ; y = (    1 - larg2 ) * PARAM%surf_dy ; x = y / t ; endif

         ! The same distance as above is calculated, from the center to the surface edge, then ...
         rd = ( (+ c0 * x + s0 * y) / tau1 ) ** 2 +  &  !
              ( (- s0 * x + c0 * y) / tau2 ) ** 2

         rr = min( rr, sqrt(rd) )

      enddo
      enddo

      rr = 0.99 * rr

      sum_inn = 0
      sum_tab = 0
      sum_int = 0
      npt_out = 0

      n = 2

      do j = 1, larg
      do i = 1, long

         x = ( i - long2 ) * PARAM%surf_dx
         y = ( j - larg2 ) * PARAM%surf_dy

         ! Distance as expressed in the theoretical acf. When r2 is 1., the acf is 0.2.
         r2 = ( (+ c0 * x + s0 * y) / tau1 ) ** 2 +  &  !
              ( (- s0 * x + c0 * y) / tau2 ) ** 2

         ! Below the correlation length, no transformation
         if ( r2 <= 1._R8 ) then

            sum_inn = sum_inn + tab_in(i, j)

            cycle

         endif

         r2 = sqrt(r2)

         if ( r2 >= rr ) then

            npt_out = npt_out + 1

            cycle

         endif

         ! ... the modified Tuckey window can be determined.

         sum_tab = sum_tab + ( cos( 0.5_R8 * PI_R8 * ( r2 - 1. ) / ( rr - 1. ) ) ** n ) * tab_in(i, j)
         sum_int = sum_int + ( r2 - 1. ) / ( rr - 1. )

      enddo
      enddo

      a_min = -( sum_inn + sum_tab ) / ( npt_out + sum_int )

      do j = 1, larg
      do i = 1, long

         x = ( i - long2 ) * PARAM%surf_dx
         y = ( j - larg2 ) * PARAM%surf_dy

         ! Distance as expressed in the theoretical acf. When r2 is 1., the acf is 0.2.
         r2 = ( (+ c0 * x + s0 * y) / tau1 ) ** 2 +  &  !
              ( (- s0 * x + c0 * y) / tau2 ) ** 2

         ! Below the correlation length, no transformation
         if ( r2 <= 1._R8 ) cycle

         r2 = sqrt(r2)

         if ( r2 >= rr ) then

            tab_out(i, j) = a_min

            cycle

         endif

         ! ... the modified Tuckey window can be determined.
         tab_out(i, j) = ( cos( 0.5_R8 * PI_R8 * ( r2 - 1. ) / ( rr - 1. ) ) ** n ) * tab_in(i, j) + a_min * ( r2 - 1. ) / ( rr - 1. )

      enddo
      enddo

      allocate( tab_tmp(1:long, 1:larg) )

      do k = 1, 10

         tab_tmp(1:long, 1:larg) = tab_out(1:long, 1:larg)

         do j = 1 + 1, larg - 1
         do i = 1 + 1, long - 1

            x = ( i - long2 ) * PARAM%surf_dx
            y = ( j - larg2 ) * PARAM%surf_dy

            ! Distance as expressed in the theoretical acf. When r2 is 1., the acf is 0.2.
            r2 = ( (+ c0 * x + s0 * y) / tau1 ) ** 2 +  &  !
                 ( (- s0 * x + c0 * y) / tau2 ) ** 2

            r2 = sqrt(r2)

            if ( r2 <= 0.98_R8 .or. ( r2 >= 1.02_R8 .and. r2 <= 0.98_R8 * rr ) .or. r2 >= 1.02_R8 * rr ) cycle

            tab_out(i, j) = ( 2*tab_tmp(i, j) +tab_tmp(i +1, j   ) +tab_tmp(i -1, j   ) + &  !
                                               tab_tmp(i   , j +1) +tab_tmp(i   , j -1) +( tab_tmp(i +1, j -1) +tab_tmp(i -1, j -1) + &  !
                                                                                           tab_tmp(i -1, j +1) +tab_tmp(i +1, j +1) )/sqrt(2._R8) )/( 6. +4./sqrt(2._R8) )

         enddo
         enddo

      enddo

      deallocate( tab_tmp )

   return
   endsubroutine apod2

endmodule func_acf
