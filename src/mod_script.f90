!< author: Arthur Francisco
!  version: 1.0.0
!  date: october, 23 2024
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.2em;">
!        **Routines to decode the script**
!  </span>

module script
!$ use omp_lib
use data_arch,       only : I4, R4, R8, UN, PI_R8, EPS_R8
use skku_profiles,   only : calculs_skku_generique, build_heights
use crest_param,     only : PARAM, JOB, SPY, TER, FCT_TANG, FCT_EXPO
use stat_mom,        only : moment_stat, calc_moments, scramble
use sort_arrays,     only : sort_array2, init_order
use miscellaneous,   only : get_unit, trans_corner2center, trans_center2corner, progress_bar_terminal
use surfile,         only : init_scal, write_surf, SCALE_SURF
use func_acf,        only : calc_imp_acf, acf_wiener
use fftw3,           only : tab_init_fftw3_real, calc_fftw3_real_bwd, calc_fftw3_real_fwd, fftw_plan_with_nthreads, end_fftw3, tab_end_fftw3_real, extend, apod, &  !
                            SINGL_FFTW_ALLOCATED, NB_THREADS_FFT, FFT_DIM, FFTW_ESTIMATE, FFTW_MEASURE, FFTW_EXHAUSTIVE
use filter,          only : fft_filter, soften
use gnufor,          only : write_xyy_plots, run_gnuplot
use anisotropy,      only : ellipse_acf
use pikaia_oop,      only : pikaia_class
use files,           only : str_remove_chars
implicit none

integer(kind=I4) :: LINE_READ, SAVE_LINE_READ

integer(kind=I4) :: I_ITER, NB_ITER

private

public :: read_job

contains

   subroutine read_job(job_file)
   !================================================================================================
   !<@note Function that reads a script file. Keywords are identified and corresponding actions are
   ! triggered.
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   character(len=512), intent(in) :: job_file   !! job file with macros to execute

      integer  (kind=I4) :: vide
      character(len=512) :: keyword
      logical  (kind=I4) :: first

      ! script file
      call get_unit( JOB )
      open(unit = JOB, file = trim(job_file), status = 'old')

      ! witness file
      call get_unit( SPY )
      open(unit = SPY, file = "out/spy.txt", status = 'unknown')

      LINE_READ = 0
      first     = .true. ! determine if calc_z_f has been already used

      ! default function used: the tangent function
      PARAM%func_gen = FCT_TANG

      ! default is 2 bounds when generating profiles
      PARAM%nparam = 2

o:    do

         keyword = repeat( ' ', len(keyword) )

         read(JOB, *, iostat = vide) keyword ; LINE_READ = LINE_READ + 1

         ! remove unwanted characters from keyword
         keyword = str_remove_chars(string = trim(keyword), chars = '- # *')

         write(SPY, '(a6,I4.4,a)') "line: ", LINE_READ, ' ', trim(keyword)

         selectcase( keyword(1:8) )

            case('STA_SCRI')
               ! start script
               call sta_scri()

            case('DEF_SIZE')
               ! image size
               call def_size()

            case('NB_PROCS')
               ! number of threads
               call nb_procs()

            case('STA_THEO')
               ! desired stat moments
               call sta_theo()

            case('ACF_THEO')
               ! desired acf
               call acf_theo()

            case('CALC_FFH')
               ! deigital filter
               call calc_ffh( calc_m_stt = first )

            case('CALC_Z_I')
               ! starting heights
               call calc_z_i()

            case('DIGI_FIL')
               ! apply digital filter
               call digi_fil()

            case('CALC_ORD')
               ! determine heights order
               call calc_ord()

            case('CALC_Z_F')
               ! final heights
               call calc_z_f( to_be_made = first )
               ! now, calc_z_f is considered as already ran
               first = .false.

            case('STA_LOOP')
               ! loop start
               call sta_loop()

            case('END_LOOP')
               ! loop end
               call end_loop()

            case('SUB_SURF')
               ! extract the best surface
               call sub_surf()

            case('SMOOTH__')
               ! low-pass filter
               call smooth__()

            case('SAVE_PRF')
               ! save image
               call save_img( tab = PARAM%surf )     ! IN

            case('CALC_ACF')
               ! determine the surface acf
               call calc_acf()

            case('PLT__ACF')
               ! print the correlation graphs and/or determine
               ! if the stop criterion is reached.
               call plt__acf()

            case('SAVE_ACF')
               ! save the acf surface
               call save_img( tab = PARAM%imp_acf )    ! IN

            case('END_SCRI')
               ! close the script reading
               call end_scri()

               exit o

         endselect

      enddo o

      close(JOB)

   return
   endsubroutine read_job


   subroutine plt__acf()
   !================================================================================================
   !<@note Function that calculates the mean absolute difference between the desired Acf and
   ! the one obtained.
   ! However, the important zone where both should match is above acf__z - where the correlation
   ! lengths are determined.
   !
   ! If the mean absolute difference is below the criterion, the loops to improve the acf are
   ! stopped.
   !
   ! The function can also plot the acfs.
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: i, w, h, mw, mh, uplot, ll
      logical(kind=I4) :: is_x, is_y
      real   (kind=R8) :: lim_crit_acf, crit_acf, dxy, l1, l2

      character(len=512) :: plt_acf

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      real(kind=R8), dimension(1:8) :: res

      read(JOB,*) lim_crit_acf ; LINE_READ = LINE_READ +1 ; write(SPY,*) LINE_READ, "Acf criterion ", lim_crit_acf
      read(JOB,*) plt_acf      ; LINE_READ = LINE_READ +1 ; write(SPY,*) LINE_READ, trim(plt_acf)

      PARAM%crt_acf = lim_crit_acf ! mean absolute difference limit

      w = PARAM%width
      h = PARAM%height

      ! Mean absolute difference between the calculated acf and the theoretical acf
      !...................................................................................
      allocate( tab_tmp(1:w, 1:h) )

         tab_tmp(1:w, 1:h) = 0

         where( PARAM%imp_acf > PARAM%acf__z ) tab_tmp = abs( PARAM%acf_surf - PARAM%imp_acf )

         crit_acf = 100 * sum( tab_tmp(1:w, 1:h) ) / count( tab_tmp(1:w, 1:h) > 0 )

         PARAM%res_acf = crit_acf

      deallocate( tab_tmp )

      write(TER,*) "acf difference ", crit_acf
      write(SPY,*) 'acf difference ', crit_acf

      ! if the acf criterion is reached, the loops stop: a means is to modify the max number
      ! of loops, so that the main loop is exited.
      if (lim_crit_acf > 0. .and. lim_crit_acf > crit_acf) NB_ITER = 1
      !...................................................................................

      ! Graphs?
      !...................................................................................
      ! if 'x' is present, plot the graph along the principal axis
      ! if 'y' is present, plot the graph along the secondary axis
      is_x = ( index(trim(plt_acf), 'x') /= 0 )
      is_y = ( index(trim(plt_acf), 'y') /= 0 )

      if ( .not.( is_x .or. is_y ) ) return

      ! ensure that w and h are odd, or act accordingly
      mw = w / 2
      mh = h / 2

      if ( w == 2 * (w/2) ) mw = w/2 + 1
      if ( h == 2 * (h/2) ) mh = h/2 + 1

      call ellipse_acf( tabin = PARAM%acf_surf(1:w, 1:h),         &  ! IN
                         long = w,                                &  ! IN
                         larg = h,                                &  ! IN
                        p_acv = res(1:8),                         &  ! OUT -> correlation lengths
                          cut = PARAM%acf__z,                     &  ! IN  -> z cut plane
                     scale_xy = [PARAM%surf_dx, PARAM%surf_dy],   &  ! IN  -> lags along x and y
                          omp = .true. )                             ! IN  -> use multithread?

      write(TER,*) res(1:2), res(4)
      write(SPY,*) 'acf lengths and roughness orientation ', res(1:2), res(4)

      ! parameters for the plot
      ll = 2 * min(mw, mh) - 3

      dxy = sqrt( PARAM%surf_dx**2 + PARAM%surf_dy**2 )

      call get_unit( uplot )

      if ( is_x ) call graph(axis = 1)

      if ( is_y ) call graph(axis = 2)

   contains

      subroutine graph(axis)
      !================================================================================================
      !<@note Function that plots the graphs to compare the ACF along the primary and/or secondary axes.
      !<@endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      integer(kind=I4), intent(in) :: axis   !! *1 or 2 for primary or secondary axis*

         character(len=256) :: file_acf, file_gpl, title
         character(len=003) :: xlab

         real(kind=R8) :: angle

         real(kind=R8), allocatable, dimension(:) :: profile_acf_surf, profile_imp_acf

         allocate( profile_acf_surf(1:ll) )
         allocate( profile_imp_acf (1:ll) )

         if ( axis == 1 ) then

            angle = PARAM%a_acf

            file_acf = 'out/gpl/acfx.txt'
            file_gpl = 'out/gpl/acfx.gpl'

            title = '"ACF comparison along primary axis X"'
            xlab  = '"X"'

         else

            angle = PARAM%a_acf + 90.

            file_acf = 'out/gpl/acfy.txt'
            file_gpl = 'out/gpl/acfy.gpl'

            title = '"ACF comparison along secondary axis Y"'
            xlab  = '"Y"'

         endif

         ! extract the ACF profile along a particular direction
         call profile_at_angle( tab = PARAM%acf_surf(1:w, 1:h), profile = profile_acf_surf(1:ll), theta = angle )
         call profile_at_angle( tab = PARAM%imp_acf (1:w, 1:h), profile = profile_imp_acf (1:ll), theta = angle )


         open(uplot, file = trim(file_acf))

            write(uplot, *) 'X', '"calculated acf"', '"theoretical acf"'

            do i = 1, ll

               write(uplot, *) (i - ll/2) * dxy, real(profile_acf_surf(i), kind=R4),    &  !
                                                 real(profile_imp_acf (i), kind=R4)        !

               if ( i - ll/2 < 0 ) then

                  if ( profile_acf_surf(i) < PARAM%acf__z .and. profile_acf_surf(i + 1) > PARAM%acf__z ) l1 = (i - ll/2) * dxy

               endif

               if ( i - ll/2 > 0 .and. i < ll ) then

                  if ( profile_acf_surf(i) > PARAM%acf__z .and. profile_acf_surf(i + 1) < PARAM%acf__z ) l2 = (i - ll/2) * dxy

               endif

            enddo

         close(uplot)

         open(uplot, file = trim(file_gpl))

            write(uplot, '(a)') 'set title ' // trim(title)
            write(uplot, '(a)') 'set xlabel ' // trim(xlab)
            write(uplot, '(a)') 'set ylabel "ACF"'

            write(uplot, '(a,f4.2,a,f5.2,a)') "set arrow from graph 0, first ", PARAM%acf__z,                                          &  !
                                                         " to graph 1, first ", PARAM%acf__z,  ' nohead lc rgb "#000000" front'           !

            write(uplot, '(a,E8.2,a,E8.2,a,f5.2,a)') "set arrow from ", l1, ", graph 0 to ",                                           &  !
                                                                        l1, ",", PARAM%acf__z, ' nohead lc rgb "#000000" front'           !

            write(uplot, '(a,E8.2,a,E8.2,a,f5.2,a)') "set arrow from ", l2, ", graph 0 to ",                                           &  !
                                                                        l2, ",", PARAM%acf__z, ' nohead lc rgb "#000000" front'           !

            write(uplot, '(a,E8.2,a,E8.2,a,f5.2)') 'set label "L1 = ', res(axis), '" at ', l2, ',', PARAM%acf__z + 0.1

            write(uplot, '(a,i2,a)' ) 'plot "' // trim(file_acf) // '" using 1:2 with lines title "acf real surface", "' // trim(file_acf) // '" using 1:3 with lines title "theoretical acf"'

            write(uplot, '(a)' ) 'pause -1  "Hit return to continue"'
            write(uplot, '(a)' ) 'q'

         close(uplot)

         call system('gnuplot "' // trim(file_gpl) // '"')

         deallocate( profile_acf_surf )
         deallocate( profile_imp_acf  )

      return
      endsubroutine graph


      subroutine profile_at_angle(tab, profile, theta)
      !================================================================================================
      !<@note Function that extract the ACF profile along a particular direction
      !<@endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      real(kind=R8), intent(in), dimension(1:w, 1:h) :: tab
      real(kind=R8), intent(in) :: theta
      real(kind=R8), intent(out), dimension(1:ll) :: profile

         integer(kind=I4) :: p, nx, ny
         real   (kind=R8) :: r, x, y, xb, yb, xm, ym, xp, yp, h1, h2, h3, h4, hh

         do p = -ll / 2, ll / 2      !  identifying a point on the diameter

            r = p                    !  corresponding algebraic radius

            !  projection on x and y of the point marked by its radius and angle
            !  by taking the lower integer, we have the number of the bottom row and left-hand column of the rectangle
            !  the remainder (x-nx) represents the abscissa of the point in the rectangle with sides=1
            !  the 0.9999 coefficient is used to avoid falling right on an existing point

            x = mw + r * cos(theta * PI_R8 / 180) * 0.9999_R8 ; nx = floor(x) ; xb = x -nx
            y = mh + r * sin(theta * PI_R8 / 180) * 0.9999_R8 ; ny = floor(y) ; yb = y -ny

            xm = UN -xb ; xp = xb
            ym = UN -yb ; yp = yb

            if ( nx+1 <= w .and. ny+1 <= h .and.   &  !
                 nx   >= 1 .and. ny   >= 1) then

               ! attention r may be greater than lo2 or la2
               h1 = tab(nx   , ny   )
               h2 = tab(nx +1, ny   )
               h3 = tab(nx +1, ny +1)
               h4 = tab(nx   , ny +1)

               hh = h1 * xm * ym + &  !
                    h2 * xp * ym + &  !
                    h3 * xp * yp + &  !
                    h4 * xm * yp      !

               profile(p + ll / 2 + 1) = hh

            endif

         enddo

      return
      endsubroutine profile_at_angle

   endsubroutine plt__acf


   subroutine calc_acf()
   !================================================================================================
   !<@note Function that returns the autocorrelation function of a surface PARAM%acf_surf
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h

      w = PARAM%width
      h = PARAM%height

      call acf_wiener(  tab_in = PARAM%surf(1:w, 1:h),          &  ! IN
                       tab_out = PARAM%acf_surf(1:w, 1:h),      &  ! OUT
                             w = w,                             &  ! IN
                             h = h  )                              ! IN

   return
   endsubroutine calc_acf


   subroutine sta_loop()
   !! Starting the loop
   implicit none

      read(JOB,*) NB_ITER ; LINE_READ = LINE_READ +1

      write(SPY,*) LINE_READ, I_ITER, '/', NB_ITER

      I_ITER         = 1         ! the account begins
      SAVE_LINE_READ = LINE_READ ! remember where to go when rewinding

   return
   endsubroutine sta_loop


   subroutine end_loop()
   !! The loop ends here
   implicit none

      integer(kind=I4) :: i_ligne

      if ( I_ITER < NB_ITER ) then

         rewind(JOB)       ! the maximum number of loops is not reached,
                           ! go to the begining of the script

      else

         I_ITER = NB_ITER  ! the maximum number of loops is reached

         return

      endif

       ! return to the beginning of the loop
      LINE_READ = SAVE_LINE_READ

      do i_ligne = 1, SAVE_LINE_READ

         read(JOB,*)

      enddo

      I_ITER = I_ITER + 1

   return
   endsubroutine end_loop


   subroutine sub_surf()
   !================================================================================================
   !<@note Function that returns the best subsurface from the final surface.
   !
   ! We are here because a non periodic resulting surface is required. To do that, a wider
   ! periodic surface is created, and it matches the required moments and acf.
   !
   ! However, sub-sampling the surface into a smaller surface that matches the required size
   ! will result in degraded moments and acf. Hence, several locations are tested to find the
   ! best subsurface.
   !
   ! Note that the right moments can always be obtained by substitution, respecting the order of heights.
   ! However, the acf will be slightly impacted.
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: i, w, h, l, we, he, dw, dh, ndw, ndh, idw, jdh, ib, ie, jb, je, best_idw, best_jdh

      integer(kind=I4) :: n_seek, inc_dw, inc_dh, nn_res, res_ratio

      real(kind=R8)    :: best_acf, res_acf, size_ratio

      character(len=100) :: text

      type(MOMENT_STAT)  :: m_res

      integer(kind=I4), dimension(1:2) :: best_ind

      integer(kind=I4), allocatable, dimension(:)   :: order_tmp

      real(kind=R8),    allocatable, dimension(:,:) :: sav_surf, surf_tmp, acf_tmp, tab_res_acf

      real(kind=R8),    allocatable, dimension(:)   :: tab_tmp

      read(JOB,*)  n_seek ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "n_seek  ", n_seek
         ! n_seek is the number of subsurfaces explored

      if ( PARAM%periodic ) return

      ! reset the FFTW configuration because from now on, the acf will be calculated on smaller surfaces
      call end_fftw3()

      ! extended surface size
      we = PARAM%width
      he = PARAM%height

      allocate( sav_surf( 1:we, 1:he ) )

      ! save the extended surface
      sav_surf( 1:we, 1:he ) = PARAM%surf( 1:we, 1:he )

      ! reset previous arrays
      deallocate( PARAM%surf     )
      deallocate( PARAM%imp_acf  )
      deallocate( PARAM%acf_surf )

      call fftw_plan_with_nthreads( nthreads = 1 )

      NB_THREADS_FFT = omp_get_max_threads()

      ! the new image size is the one of the subsurface, because it is the size defined by the user
      w  = PARAM%sub_width
      h  = PARAM%sub_height
      l  = PARAM%sub_npts

      call tab_init_fftw3_real( long = w, larg = h, plan_flag = FFTW_MEASURE )

      PARAM%width  = w
      PARAM%height = h
      PARAM%npts   = l

      allocate( PARAM%surf( 1:w, 1:h ) )
      allocate( PARAM%acf_surf( 1:w, 1:h ) )

      ! build a smaller set of heights that match the required moments
      call build_heights(     vec_out = PARAM%vect_h(1:l),                                     &  !
                         use_fct_expo = ( PARAM%m_end%ku < 1.10 * PARAM%m_end%sk**2 + 1. ),    &  !
                             stats_in = PARAM%m_end,                                           &  !
                                   lg = l )                                                       !


      allocate( PARAM%imp_acf ( 1:w, 1:h ) )

      ! recalculate the theoretical acf, for the new size surface
      call calc_imp_acf( long = w,                                &  ! IN
                         larg = h,                                &  ! IN
                         apod = PARAM%apod,                       &  ! IN
                         tau1 = PARAM%l_acf1,                     &  ! IN
                         tau2 = PARAM%l_acf2,                     &  ! IN
                        alpha = log(PARAM%acf__z),                &  ! IN
                          ang = PARAM%a_acf * PI_R8/180,          &  ! IN
                      tab_acf = PARAM%imp_acf(1:w, 1:h))             ! OUT

      ! difference between old and new size: hence, there will be dw*dh possible locations
      ! for the subsurface
      dw = we - w
      dh = he - h

      ! initialize the best acf result (mean absolute difference) and subsurface locations
      best_acf = 1.e6_R8
      best_idw = 0
      best_jdh = 0

      size_ratio = dw / dh

      ! Number of locations along x and y. It respects the total number to be explored, as set
      ! by the user and the size ratio of the surface
      ndw = nint( sqrt( n_seek / size_ratio ) )
      ndh = nint( real(n_seek, kind=R8) / ndw )

      ! ndw and ndh are modified to be mumtiples of the number of threads, but the product
      ! should not be too far from n_seek
      if ( size_ratio > 1. ) then

         ndw = NB_THREADS_FFT * ( int(ndw / NB_THREADS_FFT) + 1 )
         ndh = NB_THREADS_FFT * ( int(ndh / NB_THREADS_FFT)     )

      else

         ndw = NB_THREADS_FFT * ( int(ndw / NB_THREADS_FFT)     )
         ndh = NB_THREADS_FFT * ( int(ndh / NB_THREADS_FFT) + 1 )

      endif

      ! don't exceed the maximums
      ndw = min( ndw, dw )
      ndh = min( ndh, dh )

      ! increments for looping
      inc_dw = dw / ndw
      inc_dh = dh / ndh

      ! result storage
      allocate( tab_res_acf(1:ndw, 1:ndh) )

      tab_res_acf(1:ndw, 1:ndh) = -1

      allocate( surf_tmp( 1:w, 1:h ) )
      allocate(  acf_tmp( 1:w, 1:h ) )

      allocate( order_tmp(1:l) )
      allocate(   tab_tmp(1:l) )

      call progress_bar_terminal(val = 0, max_val = ndw * ndh, init = .true.)

      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NB_THREADS_FFT)

      !$OMP DO SCHEDULE (STATIC, ndw/NB_THREADS_FFT) PRIVATE(ib, ie, jb, je, jdh, surf_tmp, tab_tmp, order_tmp, acf_tmp, i, m_res, res_acf, nn_res, res_ratio, text)

      do idw = 1, ndw

         ib = idw * inc_dw
         ie = ib + w - 1

         do jdh = 1, ndh

            jb = jdh * inc_dh
            je = jb + h - 1

            surf_tmp( 1:w, 1:h ) = sav_surf( ib:ie, jb:je )

            !.......................................................................................
            tab_tmp(1:l) = reshape( surf_tmp(1:w, 1:h), [l] )

            ! store the subsurface height order
            call init_order( order = order_tmp(1:l),        &  ! OUT
                                 n = l )                       ! IN

            call sort_array2( tab_inout =   tab_tmp(1:l),       &  ! INOUT
                                   tab0 = order_tmp(1:l),       &  ! INOUT
                                      n = l )                      ! IN

            ! replace old heights with new ones that matches required moments
            do i = 1, l

               tab_tmp( order_tmp(i) ) = PARAM%vect_h(i)

            enddo

            surf_tmp(1:w, 1:h) = reshape( tab_tmp(1:l), [w, h] )

            call calc_moments(   tab = tab_tmp(1:l),    &  ! IN
                                  mx = m_res,           &  ! OUT
                              nb_mom = 4 )                 ! IN

            surf_tmp(1:w, 1:h) = ( surf_tmp(1:w, 1:h) - m_res%mu ) / m_res%si
            !.......................................................................................

            ! calculate the acf
            call acf_wiener(  tab_in = surf_tmp(1:w, 1:h),     &  ! IN
                             tab_out = acf_tmp(1:w, 1:h),      &  ! OUT
                                   w = w,                      &  ! IN
                                   h = h,                      &  ! IN
                           multi_fft = .true.  )                  ! IN

            ! calculate the result (mean absolute difference between theoretical acf and calculated acf)
            call calc_res_acf(     acf_surf =       acf_tmp(1:w, 1:h),     &  ! IN
                                    imp_acf = PARAM%imp_acf(1:w, 1:h),     &  ! IN
                                     acf__z = PARAM%acf__z,                &  ! IN
                                   crit_acf = res_acf,                     &  ! OUT
                                          w = w,                           &  ! IN
                                          h = h )                             ! IN

            tab_res_acf( idw, jdh ) = res_acf

            ! update progressbar
            !$OMP CRITICAL

               nn_res = count( tab_res_acf > 0 )

               call progress_bar_terminal(val = nn_res, max_val = ndw * ndh, init = .false.)

            !$OMP END CRITICAL

         enddo

      enddo

      !$OMP END DO

      !$OMP END PARALLEL

      ! find the best subsurface
      best_ind(1:2) = minloc( tab_res_acf )

      idw = best_ind(1)
      jdh = best_ind(2)

      ib = idw * ( dw / ndw )
      ie = ib + w - 1

      jb = jdh * ( dh / ndh )
      je = jb + h - 1

      PARAM%surf( 1:w, 1:h ) = sav_surf( ib:ie, jb:je )

      ! redo calculations on the best subsurface
      !.......................................................................................
      tab_tmp(1:l) = reshape( PARAM%surf(1:w, 1:h), [l] )

      call init_order( order = PARAM%order(1:l),        &  !
                           n = l )                         !

      call sort_array2( tab_inout =     tab_tmp(1:l),       &  !
                             tab0 = PARAM%order(1:l),       &  !
                                n = l )                        !

      do i = 1, l

         tab_tmp( PARAM%order(i) ) = PARAM%vect_h(i)

      enddo

      PARAM%surf(1:w, 1:h) = reshape( tab_tmp(1:l), [w, h] )

      call calc_moments(   tab = tab_tmp(1:l),    &  !
                            mx = m_res,           &  !
                        nb_mom = 4 )                 !

      PARAM%surf(1:w, 1:h) = ( PARAM%surf(1:w, 1:h) - m_res%mu ) / m_res%si
      !.......................................................................................

      call tab_end_fftw3_real()

      deallocate( sav_surf )
      deallocate( tab_tmp )

      deallocate( surf_tmp )
      deallocate( acf_tmp )
      deallocate( order_tmp )
      deallocate( tab_res_acf )

   return
   endsubroutine sub_surf


   subroutine calc_res_acf(acf_surf, imp_acf, crit_acf, acf__z, w, h)
   !================================================================================================
   !<@note Function that returns *crit_acf* the mean absolute difference between theoretical
   ! and calculated acfs
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in ) :: w                                !! *surface acf width (points)*
   integer(kind=I4), intent(in ) :: h                                !! *surface acf height (points)*
   real(kind=R8),    intent(in ) :: acf__z                           !! *plane elevation z where correlation lengths are calculated*
   real(kind=R8),    intent(out) :: crit_acf                         !! *mean absolute difference between theoretical and calculated acfs*
   real(kind=R8),    intent(in ), dimension(1:w, 1:h) :: acf_surf    !! *calculated surface acf*
   real(kind=R8),    intent(in ), dimension(1:w, 1:h) :: imp_acf     !! *required surface acf*

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      allocate( tab_tmp(1:w, 1:h) )

      tab_tmp(1:w, 1:h) = 0

      where( imp_acf > acf__z ) tab_tmp = abs( acf_surf - imp_acf )

      crit_acf = 100 * sum( tab_tmp(1:w, 1:h) ) / count( tab_tmp(1:w, 1:h) > 0 )

      deallocate( tab_tmp )

   return
   endsubroutine calc_res_acf


   subroutine save_img(tab)
   !================================================================================================
   !<@note Function that save an array *tab* as a digital surf file.
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real(kind=R8), intent(in), dimension(:,:) :: tab !! *a surface to save as a .sur file*

      integer(kind=I4) :: w ! image width
      integer(kind=I4) :: h ! image height

      integer(kind=I4), dimension(1:2) :: shape_tab

      character(len=512) :: nom_surf

      type(SCALE_SURF) :: scale_img

      shape_tab = shape( tab )
      w = shape_tab(1)
      h = shape_tab(2)

      nom_surf = repeat (" ", len(nom_surf) )

      read(JOB,*) nom_surf ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, trim(nom_surf)

      call init_scal( scal = scale_img,          &           ! out; creates a surface type, containing ...
                        nx = w,                  &           !  in; ... the number of points along x ...
                        ny = h,                  &           !  in; ... the number of points along y ...
                        lx = PARAM%surf_width,   &           !  in; ... the length (default unit : m) ...
                        ly = PARAM%surf_height,  &           !  in; ... the width ...
                    unit_z = 'm'        )                    !  in; ... and the unit along z.

      call write_surf( nom_fic = trim(nom_surf),          &  !
                         tab_s = tab(1:w, 1:h),           &  !
                          scal = scale_img )                 !

   return
   endsubroutine save_img


   subroutine smooth__()
   !================================================================================================
   !<@note Function that applies a low-pass filter to the surface PARAM%surf
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h

      real(kind=R8), dimension(:,:), allocatable :: tab_tmp

      type(MOMENT_STAT) :: m_res

      real(kind=R8) :: cutoff

      read(JOB,*) cutoff ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "cutoff ", cutoff

      PARAM%cutoff = cutoff

      w = PARAM%width
      h = PARAM%height

      allocate( tab_tmp(1:w, 1:h) )

      call fft_filter(       tab = PARAM%surf(1:w, 1:h),    &  ! IN
                            long = w,                       &  ! IN
                            larg = h,                       &  ! IN
                          cutoff = PARAM%cutoff,            &  ! IN
                          bf_tab = tab_tmp(1:w, 1:h),       &  ! OUT
                       multi_fft = .false.,                 &  ! IN
                             pad = -1._R8,                  &  ! IN
                             ext = 'constant' )                ! IN

      call calc_moments(    tab = tab_tmp(1:w, 1:h),        &  ! IN
                             mx = m_res,                    &  ! OUT
                         nb_mom = 2 )                          ! IN

      ! the surface has been modified, so recenter and rescale
      PARAM%surf(1:w, 1:h) = ( tab_tmp(1:w, 1:h) - m_res%mu ) / m_res%si

      deallocate( tab_tmp )

   return
   endsubroutine smooth__


   subroutine calc_z_f(to_be_made)
   !================================================================================================
   !<@note Function that returns PARAM%surf, the surface made of heights with the required statistical
   ! moments, in the right order.
   !
   ! - The heights come from the vector PARAM%vect_h
   ! - the heights order is stored in the vector PARAM%order
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   logical(kind=I4),  intent(in) :: to_be_made !! *whether to determine the heights, or reuse them*

      integer(kind=I4) :: i
      integer(kind=I4) :: w, h, l

      real(kind=R8), allocatable, dimension(:) :: tab_tmp

      type(MOMENT_STAT) :: m_res

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      ! final set of heights are generated to meet the desired statistical moments
      ! It is done once.
      if ( to_be_made ) then

         write(SPY,*) 'calc_z_f -> final set of heights are generated to meet the desired statistical moments'

         call build_heights(     vec_out = PARAM%vect_h(1:l),                                     &  ! OUT
                            use_fct_expo = ( PARAM%m_end%ku < 1.10 * PARAM%m_end%sk**2 + 1. ),    &  ! IN
                                stats_in = PARAM%m_end,                                           &  ! IN
                                      lg = l )                                                       ! IN

      endif

      ! The heights stored in PARAM%vect_h(1:lg) are reinjected in PARAM%surf, with respect to the order PARAM%order

      write(SPY,*) 'calc_z_f -> substitution of PARAM%surf with PARAM%vect_h with respect to PARAM%order'

      allocate( tab_tmp(1:l) )

      do i = 1, l

         tab_tmp( PARAM%order(i) ) = PARAM%vect_h(i)

      enddo

      PARAM%surf(1:w, 1:h) = reshape( tab_tmp(1:l), [w, h] )

      call calc_moments(   tab = tab_tmp(1:l),    &  ! IN
                            mx = m_res,           &  ! OUT
                        nb_mom = 2 )                 ! IN

      PARAM%surf(1:w, 1:h) = ( PARAM%surf(1:w, 1:h) - m_res%mu ) / m_res%si

      deallocate( tab_tmp )

   return
   endsubroutine calc_z_f


   subroutine calc_ord()
   !================================================================================================
   !<@note Function that returns the vector PARAM%order that contains the heights order.
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: l, w, h

      real(kind=R8), allocatable, dimension(:) :: tab_tmp

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      allocate( tab_tmp(1:l) )

      tab_tmp(1:l) = reshape( PARAM%surf(1:w, 1:h), [l] )

      call init_order( order = PARAM%order(1:l),        &  ! OUT
                           n = l )                         ! IN

      call sort_array2( tab_inout =     tab_tmp(1:l),       &  ! IN
                             tab0 = PARAM%order(1:l),       &  ! OUT
                                n = l )                        ! IN

      deallocate( tab_tmp )

   return
   endsubroutine calc_ord


   subroutine digi_fil()
   !================================================================================================
   !<@note Function that applies the digital filter to the random heights
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h

      complex(kind=R8), dimension(:,:), allocatable :: cmple
      complex(kind=R8), dimension(:,:), allocatable :: ftab         ! FFT(tab_prf)

      type(MOMENT_STAT) :: m_res

      write(SPY,*) 'digi_fil -> 1 - extends then windows PARAM%surf to calculate its FFT ftab'
      write(SPY,*) 'digi_fil -> 2 - multiplies ftab and the digital filter PARAM%fhi, then FFT-1'
      write(SPY,*) 'digi_fil -> 3 - retrieves PARAM%surf by removing the padded extension'

      w = PARAM%width
      h = PARAM%height

      allocate( cmple(1:w, 1:h) )     !

      allocate(  ftab(1:w, 1:h) )     !

      cmple(1:w, 1:h) = cmplx( PARAM%surf(1:w, 1:h), 0, kind = R8 )

      call calc_fftw3_real_fwd( tab_in = PARAM%surf(1:w, 1:h),     &  ! IN
                                tab_ou = cmple(1:w, 1:h),          &  ! OUT
                                  long = w,                        &  ! IN
                                  larg = h )                          ! IN

      ftab(1:w, 1:h) = cmple(1:w, 1:h)

      where( abs(cmple(1:w, 1:h)) > 100 * EPS_R8 )

         ftab(1:w, 1:h) = cmple(1:w, 1:h) / abs( cmple(1:w, 1:h) )

      elsewhere

         ftab(1:w, 1:h) = cmplx( UN, 0._R8, kind=R8 )

      endwhere

      write(SPY,*) 'ftab normalized'

      !--------------------------------------------------------------------
      ! FFT of the filter * FFT of the random heights
      !--------------------------------------------------------------------

      cmple(1:w, 1:h) = PARAM%fhi(1:w, 1:h) * ftab(1:w, 1:h)

      call calc_fftw3_real_bwd( tab_in = cmple(1:w, 1:h),        &  ! IN
                                tab_ou = PARAM%surf(1:w, 1:h),   &  ! OUT
                                  long = w,                      &  ! IN
                                  larg = h  )                       ! IN

      ! signal centré et normé
      call calc_moments(    tab = PARAM%surf(1:w, 1:h),   &  ! IN
                             mx = m_res,                  &  ! OUT
                         nb_mom = 4 )                        ! IN

      PARAM%surf(1:w, 1:h) = ( PARAM%surf(1:w, 1:h) - m_res%mu ) / m_res%si

      write(SPY,*) 'sk ku fin ', m_res%sk, m_res%ku

      deallocate( cmple )
      deallocate( ftab )

   return
   endsubroutine digi_fil


   subroutine calc_z_i()
   !================================================================================================
   !<@note Function that returns the starting surface of random heights
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer (kind=I4) :: w, h, l

      real(kind=R8), allocatable, dimension(:) :: tab_tmp

      type(MOMENT_STAT) :: m_res

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      allocate( tab_tmp(1:l) )

      ! starting set of heights are generated to meet the statistical moments prescribed by calc_ffh

      write(SPY,*) 'calc_z_i -> starting set of heights are generated to meet the prescribed statistical moments by calc_ffh'

      call build_heights(     vec_out = tab_tmp(1:l),                                        &  ! OUT
                         use_fct_expo = ( PARAM%m_stt%ku < 1.34 * PARAM%m_stt%sk**2 + 1. ),  &  ! IN
                             stats_in = PARAM%m_stt,                                         &  ! IN
                                   lg = l )                                                     ! IN

      call calc_moments(   tab = tab_tmp(1:l),      &  ! IN
                            mx = m_res,             &  ! OUT
                        nb_mom = 4 )                   ! IN

      write(TER,*) 'starting statistical moments ', m_res%mu, m_res%va, m_res%sk, m_res%ku

      tab_tmp(1:l) = ( tab_tmp(1:l) - m_res%mu ) / m_res%si

      call scramble( tab = tab_tmp(1:l),   &  ! INOUT
                      lg = l )                ! IN


      PARAM%surf(1:w, 1:h) = reshape( tab_tmp(1:l), [w, h] )

      deallocate( tab_tmp )

   return
   endsubroutine calc_z_i


   subroutine calc_ffh(calc_m_stt)
   !================================================================================================
   !<@note Function that returns ...
   !
   ! - the digital filter to apply to the height distribution \( \text{PARAM%fhi} = \sqrt{ \left| FFT(\text{PARAM%imp_acf}) \right| } \)
   ! - the starting statistical moments PARAM%m_stt%sk, PARAM%m_stt%ku
   ! - whether the exponential function will be used, PARAM%reajust_skku
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   logical(kind=I4), intent(inout) :: calc_m_stt  !! *compute starting moments ?*

      integer(kind=I4) :: w, h, l

      complex(kind=R8), dimension(:,:), allocatable :: cmple
      real   (kind=R8), dimension(:,:), allocatable :: tab_tmp

      type(MOMENT_STAT) :: m_h

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      allocate( cmple(1:w, 1:h)  )

      write(SPY,*) 'calc_ffh -> PARAM%fhi = sqrt( abs( FFT(PARAM%imp_acf) ) )'

      call calc_fftw3_real_fwd( tab_in = PARAM%imp_acf(1:w, 1:h),  &  ! IN
                                tab_ou = cmple(1:w, 1:h),          &  ! OUT
                                  long = w,                        &  ! IN
                                  larg = h )                          ! IN

      PARAM%fhi(1:w, 1:h) = sqrt( abs( cmple(1:w, 1:h) ) )

      if ( calc_m_stt ) then

         ! determine starting statistical moments, if not already done
         !.......................................................................
         write(SPY,*) 'calc_ffh -> PARAM%m_stt calculated'

         cmple(1:w, 1:h) = cmplx( PARAM%fhi(1:w, 1:h), 0, kind = R8 )

         allocate( tab_tmp(1:w, 1:h) )

         call calc_fftw3_real_bwd( tab_in = cmple(1:w, 1:h),        &  ! IN
                                   tab_ou = tab_tmp(1:w, 1:h),      &  ! OUT
                                     long = w,                      &  ! IN
                                     larg = h)                         ! IN

         call calc_moments( tab     = tab_tmp(1:w, 1:h),      &  ! IN
                             mx     = m_h,                    &  ! OUT
                             nb_mom = 4 )                        ! IN

         PARAM%m_stt%mu = 0
         PARAM%m_stt%va = 1
         PARAM%m_end%mu = 0
         PARAM%m_end%va = m_h%va

         PARAM%m_stt%sk = sqrt( UN*l ) *  PARAM%m_end%sk      /  m_h%sk
         PARAM%m_stt%ku =          l   * (PARAM%m_end%ku -3.) / (m_h%ku - 3.) + 3.

         PARAM%reajust_skku = .false.

         if ( PARAM%m_stt%ku < PARAM%m_stt%sk**2 + 1.) then

            PARAM%m_stt%ku = PARAM%m_stt%sk**2 + 1.

            PARAM%reajust_skku = .true.

         endif

         deallocate( tab_tmp )
         !.......................................................................

      else

         write(SPY,*) 'calc_ffh -> PARAM%m_stt NOT calculated, set to PARAM%m_end'

         PARAM%m_stt%sk = PARAM%m_end%sk
         PARAM%m_stt%ku = PARAM%m_end%ku

      endif

      deallocate( cmple )

   return
   endsubroutine calc_ffh


   subroutine sta_scri()
      !! Start the script reading
   implicit none

      ! Initializes the state of the pseudorandom number generator used by RANDOM_NUMBER.
      call random_init(repeatable = .false., image_distinct = .true.)

   return
   endsubroutine sta_scri


   subroutine def_size()
      !! Geometrical characteristics of the numerical surface
   implicit none

      integer(kind=I4) :: w, h

      logical(kind=I4) :: period

      real(kind=R8)    :: lw, lh

      read(JOB,*)  w,  h ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "width, height  ", w, h
      read(JOB,*) lw, lh ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "img width and height (m) ", lw, lh
      read(JOB,*) period ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "periodic surface? ", period

      PARAM%periodic = period

      PARAM%width  = w
      PARAM%height = h
      PARAM%npts   = w * h

      PARAM%surf_width  = lw
      PARAM%surf_height = lh

      PARAM%surf_dx = lw / w
      PARAM%surf_dy = lh / h

   return
   endsubroutine def_size


   subroutine acf_theo()
   !================================================================================================
   !<@note Function that returns the theoretical acf PARAM%imp_acf.
   !
   ! If the surface to generate is non periodic, the starting surface is extended. The final surface
   ! will be a part of it. Indeed the extended surface will be periodic, because the use of FFTs.

   ! If a roughness orientation is chosen, in addition with long correlation lengths, a windowing
   ! should be applied to the acf to prevent from artifacts (vertical and horizontal lines)
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h
      logical(kind=I4) :: with_apod

      real(kind=R8) :: ratio, a, b, c, s, lx, ly

      real(kind=R8), dimension(1:8) :: res

      read(JOB,*) PARAM%l_acf1, PARAM%l_acf2, PARAM%acf__z ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "l_acf ", PARAM%l_acf1, PARAM%l_acf2, PARAM%acf__z

      if (PARAM%l_acf1 < PARAM%l_acf2) stop "inverser l_acf1, l_acf2"

      read(JOB,*) PARAM%a_acf  ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "a_acf ", PARAM%a_acf

      write(SPY,*) 'acf_theo -> determines the theoretical acf with a padded size, so correlation lengths are adjusted accordingly'

      read(JOB,*) with_apod    ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "apod ", with_apod

      PARAM%apod = with_apod

      w = PARAM%width
      h = PARAM%height

      if ( .not.PARAM%periodic ) then

         PARAM%sub_width  = w
         PARAM%sub_height = h
         PARAM%sub_npts   = w * h

         PARAM%sub_surf_width  = PARAM%surf_width
         PARAM%sub_surf_height = PARAM%surf_height

         ! the surface must be extended

         a = PARAM%l_acf1
         b = PARAM%l_acf2

         c = cos( PARAM%a_acf * PI_R8 / 180 )
         s = sin( PARAM%a_acf * PI_R8 / 180 )

         lx = sqrt( (a * c)**2 + (b * s)**2 )
         ly = sqrt( (a * s)**2 + (b * c)**2 )

         w = w + nint( lx / PARAM%surf_dx )
         h = h + nint( ly / PARAM%surf_dy )

         ! update sizes
         PARAM%width  = w
         PARAM%height = h
         PARAM%npts   = w * h

         PARAM%surf_width  = w * PARAM%surf_dx
         PARAM%surf_height = h * PARAM%surf_dy

      endif

      allocate( PARAM%surf(1:w, 1:h) )

      allocate( PARAM%fhi(1:w, 1:h) )

      allocate( PARAM%imp_acf(1:w, 1:h) )

      allocate( PARAM%order(1:w * h) )

      allocate( PARAM%vect_h(1:w * h) )

      allocate( PARAM%acf_surf(1:w, 1:h) )

      ratio = PARAM%l_acf2 / PARAM%l_acf1

      ! acf_theo is the theoretical acf for the normal surface
      call calc_imp_acf( long = w,                                &  ! IN
                         larg = h,                                &  ! IN
                         apod = with_apod,                        &  ! IN
                         tau1 = PARAM%l_acf1,                     &  ! IN
                         tau2 = PARAM%l_acf2,                     &  ! IN
                        alpha = log(PARAM%acf__z),                &  ! IN
                          ang = PARAM%a_acf * PI_R8/180,          &  ! IN
                      tab_acf = PARAM%imp_acf(1:w, 1:h))             ! OUT

      call ellipse_acf( tabin = PARAM%imp_acf(1:w, 1:h),          &  ! IN
                         long = w,                                &  ! IN
                         larg = h,                                &  ! IN
                        p_acv = res(1:8),                         &  ! OUT
                          cut = PARAM%acf__z,                     &  ! IN
                     scale_xy = [PARAM%surf_dx, PARAM%surf_dy],   &  ! IN
                          omp = .true. )                             ! IN

   return
   endsubroutine acf_theo


   subroutine nb_procs()
      !! Number of concurrent threads
   implicit none

      integer(kind=I4) :: nb_th

      read(JOB,*) nb_th ; LINE_READ = LINE_READ + 1 ; write(SPY, *) LINE_READ, 'nb_procs', nb_th

      select case( nb_th )

         case( 0) ! no multihreading
            PARAM%nb_threads = 1
            NB_THREADS_FFT   = 1

         case(-1) ! determined by system
            PARAM%nb_threads = omp_get_num_procs()
            NB_THREADS_FFT   = PARAM%nb_threads

         case default
            stop 'Bad choice "nb_procs" in "mod_script"'

      endselect

   return
   endsubroutine nb_procs


   subroutine sta_theo()
      !! Required statistical moments
   implicit none

      read(JOB,*) PARAM%m_end%sk, PARAM%m_end%ku ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "m_end ", PARAM%m_end%sk, PARAM%m_end%ku

   return
   endsubroutine sta_theo


   subroutine end_scri()
      !! End of script
   implicit none

      close(SPY)

      if ( allocated( PARAM%imp_acf  ) ) deallocate( PARAM%imp_acf  )
      if ( allocated( PARAM%acf_surf ) ) deallocate( PARAM%acf_surf )
      if ( allocated( PARAM%surf     ) ) deallocate( PARAM%surf     )
      if ( allocated( PARAM%vect_h   ) ) deallocate( PARAM%vect_h   )
      if ( allocated( PARAM%fhi      ) ) deallocate( PARAM%fhi      )
      if ( allocated( PARAM%order    ) ) deallocate( PARAM%order    )


      call end_fftw3()

   return
   endsubroutine end_scri


endmodule script

