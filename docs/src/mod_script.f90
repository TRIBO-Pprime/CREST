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
use skku_profiles,   only : build_heights
use crest_param,     only : PARAM, JOB, SPY, STA, TER, FCT_TANG, FCT_EXPO, LINE_READ, SAVE_LINE_READ, I_ITER, NB_ITER, SCALE_IMG
use stat_mom,        only : moment_stat, calc_moments, scramble, std_array
use sort_arrays,     only : sort_array2, init_order
use miscellaneous,   only : get_unit, trans_corner2center, trans_center2corner, progress_bar_terminal
use surfile,         only : init_scal, read_surf, write_surf, SCALE_SURF
use func_acf,        only : calc_imp_acf, acf_wiener, apod2
use fftw3,           only : tab_init_fftw3_real, calc_fftw3_real_bwd, calc_fftw3_real_fwd, fftw_plan_with_nthreads, end_fftw3, tab_end_fftw3_real, extend, apod, &  !
                            SINGL_FFTW_ALLOCATED, NB_THREADS_FFT, FFT_DIM, FFTW_ESTIMATE, FFTW_MEASURE, FFTW_EXHAUSTIVE, BACKWARD, FORWARD, calc_fftw3
use filter,          only : fft_filter, soften
use gnufor,          only : write_xyy_plots, run_gnuplot
use anisotropy,      only : ellipse_acf
use pikaia_oop,      only : pikaia_class
use files,           only : str_remove_chars, dir_separator, mkdir
use analyses,        only : surface_analysis

implicit none

private

public :: read_job

contains

   subroutine read_job(job_file)
   !================================================================================================
   !<@note Function that reads a script file. Keywords are identified and corresponding actions are
   !< triggered.
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   character(len=512), intent(in) :: job_file   !! *job file with macros to execute*

      integer  (kind=I4) :: vide
      character(len=512) :: keyword
      logical  (kind=I4) :: first

      ! script file
      call get_unit( JOB )
      open(unit = JOB, file = trim(job_file), status = 'old')

      call check_empty_lines()

      ! witness file
      call get_unit( SPY )
      open(unit = SPY, file = "out/spy.txt", status = 'unknown')

      LINE_READ = 0
      first     = .true. ! determine if calc_z_f has been already used

      ! default function used: the tangent function
      PARAM%func_gen = FCT_TANG

      ! default is 2 bounds when generating profiles
      PARAM%nparam = 2

      PARAM%calc_mstt = .true.

      PARAM%calc_zf = .true.

      ! default cutting plane for correlation lengths
      PARAM%orig_surf%cut = 0.5
      PARAM%curr_surf%cut = 0.5

o:    do

         keyword = repeat( ' ', len(keyword) )

         read(JOB, *, iostat = vide) keyword ; LINE_READ = LINE_READ + 1

         ! remove unwanted characters from keyword
         keyword = str_remove_chars(string = trim(keyword), chars = '- # *')

         write(SPY, '(a6,I4.4,a)') "line: ", LINE_READ, ' ', trim(keyword)

         selectcase( keyword(1:8) )

            case('ACF_THEO')
               ! desired acf
               call acf_theo()

            case('ADD_NOIS')
               ! add noise to current surface
               call add_nois()

            case('ALLOCATE')
               ! allocate tabs in PARAM type
               call alloc_tabs()

            case('ANALYSIS')
               call surface_analysis()

            case('APOD_ACF')
               ! apodize acf
               call apod_acf()

            case('APOD_SUR')
               ! apodize surface
               call apod_sur()

            case('CALC_ACF')
               ! determine the surface acf
               call calc_acf()

            case('CALC_FFH')
               ! deigital filter
               call calc_ffh()

            case('CALC_ORD')
               ! determine height order
               call calc_ord()

            case('CALC_Z_F')
               ! final heights
               call calc_z_f()

            case('CALC_Z_I')
               ! starting heights
               call calc_z_i()

            case('DEF_SIZE')
               ! image size
               call def_size()

            case('DIGI_FIL')
               ! apply digital filter
               call digi_fil()

            case('END_LOOP')
               ! loop end
               call end_loop()

            case('END_SCRI')
               ! close the script reading
               call end_scri()

               exit o

            case('MAKE_MSK')
               ! turn image into mask
               call make_msk()

            case('MAKE_SCR')
               ! make secratches
               call make_scratches()

            case('MAKE_TEX')
               ! build a texture
               call make_tex()

            case('NB_PROCS')
               ! number of threads
               call nb_procs()

            case('PLT__ACF')
               ! print the correlation graphs and/or determine
               ! if the stop criterion is reached.
               call plt__acf()

            case('READ_PRF')
               ! read image
               call read_img()

            case('REPR_IMG')
               ! reproduce image
               call repr_img()

            case('SAVE_ACF')
               ! save the acf surface
               if ( sum( PARAM%imp_acf ) == 0 ) then
                  call save_img( tab = PARAM%acf_surf )     ! IN: real acf
               else
                  call save_img( tab = PARAM%imp_acf )      ! IN: wanted acf
               endif

            case('SAVE_PRF')
               ! save image
               call save_img( tab = PARAM%surf )            ! IN: surface to save

            case('SMOOTH__')
               ! low-pass filter
               call smooth__()

            case('SPCT_SUR')
               ! surface spectrum frequencies
               call spct_sur()

            case('STA_LOOP')
               ! loop start
               call sta_loop()

            case('STA_SCRI')
               ! start script
               call sta_scri()

            case('STA_THEO')
               ! desired stat moments
               call sta_theo()

            case('STAT_SUR')
               ! surface moments as reference
               call stat_sur()

            case('SUB_SURF')
               ! extract the best surface
               call sub_surf()

         endselect

      enddo o

      close(JOB)

   contains

      subroutine check_empty_lines()
         !! Check for empty lines in the script
      implicit none

         do

            read(JOB, '(A)', iostat = vide) keyword

            if (vide < 0) then

               rewind(JOB)
               return

            endif

            if (len_trim(keyword) == 0) stop 'Empty line in script file'

         enddo

      return
      endsubroutine check_empty_lines

   endsubroutine read_job


   subroutine alloc_tabs()
   !================================================================================================
   !<@note Function that allocates arrays in global variable [[PARAM]]
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h

      w = PARAM%width
      h = PARAM%height

      allocate( PARAM%fhi(1:w, 1:h) )       ; PARAM%fhi(1:w, 1:h) = 0         ! *digital filter*

      allocate( PARAM%imp_acf(1:w, 1:h) )   ; PARAM%imp_acf(1:w, 1:h) = 0     ! *imposed autocorrelation*

      allocate( PARAM%order(1:w * h) )      ; PARAM%order(1:w * h) = 0        ! *vector that stores heights order*

      allocate( PARAM%vect_h(1:w * h) )     ; PARAM%vect_h(1:w * h) = 0       ! *vector used to store the heights that meet the stat moments*

      allocate( PARAM%acf_surf(1:w, 1:h) )  ; PARAM%acf_surf(1:w, 1:h) = 0    ! *calculated autocorrelation*

      allocate( PARAM%surf_copy(1:w, 1:h) ) ; PARAM%surf_copy(1:w, 1:h) = 0   ! *surface array copy*

      allocate( PARAM%surf_LF(1:w, 1:h) )   ; PARAM%surf_LF(1:w, 1:h) = 0     ! *surface low frequencies*

      allocate( PARAM%surf_HF(1:w, 1:h) )   ; PARAM%surf_HF(1:w, 1:h) = 0     ! *surface high frequencies*

   return
   endsubroutine alloc_tabs


   subroutine plt__acf()
   !================================================================================================
   !<@note Function that calculates the mean absolute difference between the desired Acf and
   !< the one obtained.
   !< However, the important zone where both should match is above the cut - where the correlation
   !< lengths are determined.
   !<
   !< If the mean absolute difference is below the criterion, the loops to improve the acf are
   !< stopped.
   !<
   !< The function can also plot the acfs.
   !<
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
      if ( PARAM%curr_surf%cut > 0 ) then

         allocate( tab_tmp(1:w, 1:h) )

            tab_tmp(1:w, 1:h) = 0

            where( PARAM%imp_acf > PARAM%curr_surf%cut ) tab_tmp = abs( PARAM%acf_surf - PARAM%imp_acf )

            crit_acf = 100 * sum( tab_tmp(1:w, 1:h) ) / count( tab_tmp(1:w, 1:h) > 0 )

            PARAM%res_acf = crit_acf

         deallocate( tab_tmp )

         write(TER,*) "acf difference ", crit_acf
         write(SPY,*) 'acf difference ', crit_acf

         ! if the acf criterion is reached, the loops stop: a means is to modify the max number
         ! of loops, so that the main loop is exited.
         if (lim_crit_acf > 0. .and. lim_crit_acf > crit_acf) NB_ITER = 1

      endif

      !...................................................................................

      ! Graphs?
      !...................................................................................
      ! if 'x' is present, plot the graph along the principal axis
      ! if 'y' is present, plot the graph along the secondary axis
      is_x = ( index(trim(plt_acf), 'x') /= 0 )
      is_y = ( index(trim(plt_acf), 'y') /= 0 )

      if ( .not.( is_x .or. is_y ) ) return

      mw = w / 2 + 1
      mh = h / 2 + 1

      call ellipse_acf( tabin = PARAM%acf_surf(1:w, 1:h),         &  ! IN
                         long = w,                                &  ! IN
                         larg = h,                                &  ! IN
                        p_acv = res(1:8),                         &  ! OUT -> correlation lengths
                          cut = PARAM%curr_surf%cut,              &  ! IN  -> z cut plane
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

         integer(kind=I4) :: exit_status

         logical(kind=I4) :: dir_exists

         character(len=256) :: file_acf, file_gpl, title, cwd
         character(len=003) :: xlab

         character(len=1) :: os_sep

         real(kind=R8) :: angle

         real(kind=R8), allocatable, dimension(:) :: profile_acf_surf, profile_imp_acf

         allocate( profile_acf_surf(1:ll) )
         allocate( profile_imp_acf (1:ll) )

         ! if "out/gpl" does not exist, create it
         inquire(file = "out/gpl", exist = dir_exists)
         if ( .not. dir_exists ) then

            os_sep = dir_separator()
            call getcwd( cwd )
            call mkdir(wkd = trim(cwd), directory = "out/gpl", sep = os_sep, exit_status = exit_status)

         endif

         if ( axis == 1 ) then

            angle = PARAM%curr_surf%ang

            file_acf = 'out/gpl/acfx.txt'
            file_gpl = 'out/gpl/acfx.gpl'

            title = '"ACF comparison along primary axis X"'
            xlab  = '"X"'

         else

            angle = PARAM%curr_surf%ang + 90.

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

                  if ( profile_acf_surf(i) < PARAM%curr_surf%cut .and. profile_acf_surf(i + 1) > PARAM%curr_surf%cut ) l1 = (i - ll/2) * dxy

               endif

               if ( i - ll/2 > 0 .and. i < ll ) then

                  if ( profile_acf_surf(i) > PARAM%curr_surf%cut .and. profile_acf_surf(i + 1) < PARAM%curr_surf%cut ) l2 = (i - ll/2) * dxy

               endif

            enddo

         close(uplot)

         open(uplot, file = trim(file_gpl))

            write(uplot, '(a)') 'set title ' // trim(title)
            write(uplot, '(a)') 'set xlabel ' // trim(xlab)
            write(uplot, '(a)') 'set ylabel "ACF"'

            write(uplot, '(a,f4.2,a,f5.2,a)') "set arrow from graph 0, first ", PARAM%curr_surf%cut,                                          &  !
                                                         " to graph 1, first ", PARAM%curr_surf%cut,  ' nohead lc rgb "#000000" front'           !

            write(uplot, '(a,E8.2,a,E8.2,a,f5.2,a)') "set arrow from ", l1, ", graph 0 to ",                                           &  !
                                                                        l1, ",", PARAM%curr_surf%cut, ' nohead lc rgb "#000000" front'           !

            write(uplot, '(a,E8.2,a,E8.2,a,f5.2,a)') "set arrow from ", l2, ", graph 0 to ",                                           &  !
                                                                        l2, ",", PARAM%curr_surf%cut, ' nohead lc rgb "#000000" front'           !

            write(uplot, '(a,E8.2,a,E8.2,a,f5.2)') 'set label "L1 = ', res(axis), '" at ', l2, ',', PARAM%curr_surf%cut + 0.1

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
   !<@note Function that returns the autocorrelation function of a surface in PARAM%acf_surf
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h
      logical(kind=I4) :: set_acf

      ! if set_acf is true, the acf becomes the prescribed one
      read(JOB,*) set_acf ; LINE_READ = LINE_READ +1 ; write(SPY,*) LINE_READ, "Set ACF ", set_acf

      w = PARAM%width
      h = PARAM%height

      call acf_wiener(  tab_in = PARAM%surf(1:w, 1:h),          &  ! IN
                       tab_out = PARAM%acf_surf(1:w, 1:h),      &  ! OUT
                             w = w,                             &  ! IN
                             h = h  )                              ! IN

      if ( set_acf ) PARAM%imp_acf(1:w, 1:h) = PARAM%acf_surf(1:w, 1:h)

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
   !<
   !< We are here because a non periodic resulting surface is required. To do that, a wider
   !< periodic surface is created, and it matches the required moments and acf.
   !<
   !< However, sub-sampling the surface into a smaller surface that matches the required size
   !< will result in degraded moments and acf. Hence, several locations are tested to find the
   !< best subsurface.
   !<
   !< Note that the right moments can always be obtained by substitution, respecting the order of heights.
   !< However, the acf will be slightly impacted.
   !<
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
                         use_fct_expo = ( PARAM%m_end%ku < 1.34 * PARAM%m_end%sk**2 + 1.8 ),   &  !
                             stats_in = PARAM%m_end,                                           &  !
                                   lg = l )                                                       !


      allocate( PARAM%imp_acf ( 1:w, 1:h ) )

      ! recalculate the theoretical acf, for the new size surface
      call calc_imp_acf( long = w,                                &  ! IN
                         larg = h,                                &  ! IN
                         apod = PARAM%apod,                       &  ! IN
                         tau1 = PARAM%curr_surf%cl1,              &  ! IN
                         tau2 = PARAM%curr_surf%cl2,              &  ! IN
                        alpha = log(PARAM%curr_surf%cut),         &  ! IN
                          ang = PARAM%curr_surf%ang * PI_R8/180,  &  ! IN
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
                                     acf__z = PARAM%curr_surf%cut,         &  ! IN
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
   !< and calculated acfs, above *z* (usually 0.2 as recommended by iso 25178)
   !<
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


   subroutine read_img()
   !================================================================================================
   !<@note Function that reads a digital surf file and returns the surface in PARAM%surf
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      character(len=512) :: nom_surf

      nom_surf = repeat (" ", len(nom_surf) )

      read(JOB,*) nom_surf ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, trim(nom_surf)

      ! read the surface, no scaling, no centering
      call read_surf( nom_fic = trim(nom_surf), &  ! IN
                        tab_s = PARAM%surf,     &  ! OUT
                         scal = SCALE_IMG )        ! OUT

      PARAM%width       = SCALE_IMG%xres
      PARAM%height      = SCALE_IMG%yres
      PARAM%npts        = SCALE_IMG%xres * SCALE_IMG%yres
      PARAM%surf_dx     = SCALE_IMG%lx / SCALE_IMG%xres
      PARAM%surf_dy     = SCALE_IMG%ly / SCALE_IMG%yres
      PARAM%surf_width  = SCALE_IMG%lx
      PARAM%surf_height = SCALE_IMG%ly

      write(*,*) "width ",   PARAM%width,      " height ", PARAM%height, " npts ", PARAM%npts
      write(*,*) "surf dx ", PARAM%surf_dx,    " dy ",     PARAM%surf_dy
      write(*,*) "width ",   PARAM%surf_width, " height ", PARAM%surf_height

   return
   endsubroutine read_img


   subroutine make_msk()
   !================================================================================================
   !<@note Function that reads a digital surf file and turns it into a mask
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h

      w = PARAM%width
      h = PARAM%height

      allocate(PARAM%surf_msk(1:w, 1:h))

      PARAM%surf_msk(1:w, 1:h) = nint( PARAM%surf(1:w, 1:h) )

   return
   endsubroutine make_msk


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

      type(SCALE_SURF) :: scale_img_tmp

      shape_tab = shape( tab )
      w = shape_tab(1)
      h = shape_tab(2)

      nom_surf = repeat (" ", len(nom_surf) )

      read(JOB,*) nom_surf ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, trim(nom_surf)

      call init_scal( scal = scale_img_tmp,      &           ! out; creates a surface type, containing ...
                        nx = w,                  &           !  in; ... the number of points along x ...
                        ny = h,                  &           !  in; ... the number of points along y ...
                        lx = PARAM%surf_width,   &           !  in; ... the length (default unit : m) ...
                        ly = PARAM%surf_height,  &           !  in; ... the width ...
                    unit_z = 'm'        )                    !  in; ... and the unit along z.

      call write_surf( nom_fic = trim(nom_surf),          &  !
                         tab_s = tab(1:w, 1:h),           &  !
                          scal = scale_img_tmp )             !

   return
   endsubroutine save_img


   subroutine spct_sur(file_spct, apod)
   !================================================================================================
   !<@note Returns the default surface spectrum
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   character(len=*), intent(in), optional :: file_spct   !! *txt file containing the surface FFT module*
   logical(kind=I4), intent(in), optional :: apod        !! *window applied to surface?*

      integer(kind=I4) :: w, h, ww, hh, i, j, np, txt_file
      logical(kind=I4) :: apod_surf

      character(len=512) :: str

      real(kind=R8),    dimension(:,:), allocatable :: tab_tmp, tab_freq1, tab_freq2, tab_spc
      complex(kind=R8), dimension(:,:), allocatable :: tab_cmpl, cmpl1

      w = PARAM%width
      h = PARAM%height

      if ( present(file_spct) ) then

         ! use in the program -> arguments passed to subroutine
         str = file_spct
         apod_surf = apod

      else

         ! use in script -> arguments passed to batch file
         read(JOB,*) str       ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, 'file_spct: ', trim(str)
         read(JOB,*) apod_surf ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, 'apod surf: ', apod_surf

      endif

      allocate( tab_tmp(1:w, 1:h) )
      allocate( tab_freq1(1:w, 1:h), tab_freq2(1:w, 1:h) )
      allocate( tab_cmpl(1:w, 1:h), cmpl1(1:w, 1:h) )

      if ( apod_surf ) then

         tab_tmp(1:w, 1:h) = PARAM%surf(1:w, 1:h)

         call apod_sur()

      endif

      cmpl1(1:w, 1:h) = cmplx( PARAM%surf(1:w, 1:h), 0, kind = R8 )

      call calc_fftw3(   sens = FORWARD,                    &  ! IN
                       tab_in = cmpl1(1:w, 1:h),            &  ! IN
                       tab_ou = tab_cmpl(1:w, 1:h),         &  ! OUT
                         long = w,                          &  ! IN
                         larg = h)                             ! IN


      tab_freq1(1:w, 1:h) = log10( abs( tab_cmpl(1:w, 1:h) ) + UN )

      call trans_corner2center(  tab_in  = tab_freq1(1:w, 1:h),  &  ! IN
                                 tab_out = tab_freq2(1:w, 1:h),  &  ! OUT
                                 long    = w,                    &  ! IN
                                 larg    = h  )                     ! IN

      np = 50     ! reduce image size
      ww = w / np
      hh = h / np

      allocate( tab_spc(1:ww, 1:hh) )

      do i = 1, ww
      do j = 1, hh

         tab_spc(i, j) = sum( tab_freq2( (i - 1) * np + 1 : i * np,  &  !
                                         (j - 1) * np + 1 : j * np  ) ) !

      enddo
      enddo

      call get_unit( txt_file )

      ! Ouvrir un fichier binaire
      open( newunit = txt_file, file = trim(str), action = "write" )

         ! Écrire les dimensions (pour Python)
         write(txt_file, *) ww, hh

         ! Écrire les données en mémoire contiguë
         write(txt_file, *) tab_spc(1:ww, 1:hh)

      close(txt_file)

      if ( apod_surf ) PARAM%surf(1:w, 1:h) = tab_tmp(1:w, 1:h)

      deallocate( tab_tmp, tab_freq1, tab_freq2, tab_spc )
      deallocate( tab_cmpl, cmpl1 )

   return
   endsubroutine spct_sur


   subroutine stat_sur()
   !================================================================================================
   !<@note Define surface statistical moments as reference
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h, l
      logical(kind=I4) :: apodize

      read(JOB,*) apodize ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, 'apod? ', apodize

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      PARAM%surf_copy(1:w, 1:h) = PARAM%surf(1:w, 1:h)

      if ( apodize ) then

         ! apodize surface for the acf
         call apod_sur()

         call acf_wiener(  tab_in = PARAM%surf(1:w, 1:h),          &  ! IN
                          tab_out = PARAM%acf_surf(1:w, 1:h),      &  ! OUT
                                w = w,                             &  ! IN
                                h = h  )                              ! IN

         ! retrieve original surface
         PARAM%surf(1:w, 1:h) = PARAM%surf_copy(1:w, 1:h)

      endif

      ! do not append to result file, make new one
      call surface_analysis( app = 0 )

      ! calculate statistical moments that become prescribed
      call std_array(tab = PARAM%surf(1:w, 1:h), mx = PARAM%m_end)

      ! determine final heights
      PARAM%vect_h(1:l) = reshape( PARAM%surf(1:w, 1:h), [l] )

      ! determine height order
      call sort_array2(tab_inout = PARAM%vect_h(1:l), n = l)

      ! prescribe moments
      PARAM%m_ini = PARAM%m_end

      ! prescribe acf
      PARAM%imp_acf(1:w, 1:h) = PARAM%acf_surf(1:w, 1:h)

   return
   endsubroutine stat_sur


   subroutine smooth__()
   !================================================================================================
   !<@note Function that applies a low-pass filter to the surface PARAM%surf
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h

      real(kind=R8), dimension(:,:), allocatable :: tab_tmp

      real(kind=R8) :: cutoff

      read(JOB,*) cutoff ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "cutoff ", cutoff

      PARAM%cutoff = cutoff

      w = PARAM%width
      h = PARAM%height

      ! forces FFT reinitialization
      FFT_DIM = 0

      allocate( tab_tmp(1:w, 1:h) )

      ! remark: if cutoff is negative, the filter is a top-hat instead of Gaussian
      call fft_filter(       tab = PARAM%surf(1:w, 1:h),    &  ! IN
                            long = w,                       &  ! IN
                            larg = h,                       &  ! IN
                          cutoff = PARAM%cutoff,            &  ! IN
                          bf_tab = tab_tmp(1:w, 1:h),       &  ! OUT
                       multi_fft = .false.,                 &  ! IN
                             pad = -1._R8,                  &  ! IN
                             ext = 'constant' )                ! IN

      ! standardize surface
      call std_array( tab = tab_tmp(1:w, 1:h) )

      PARAM%surf(1:w, 1:h) = tab_tmp(1:w, 1:h)

      deallocate( tab_tmp )

      FFT_DIM = 0

   return
   endsubroutine smooth__


   subroutine apod_sur()
   !================================================================================================
   !<@note Function that apodize the reference surface for acf calculus purposes
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      w = PARAM%width
      h = PARAM%height

      allocate(  tab_tmp(1:w, 1:h) )

      call apod(  tab_in = PARAM%surf(1:w, 1:h),   &  !  IN
                 tab_out = tab_tmp(1:w, 1:h),      &  !  OUT
                    long = w,                      &  !  IN
                    larg = h,                      &  !  IN
                type_apo = 'tuckey' )                 !  IN

      ! standardize surface
      call std_array( tab_tmp(1:w, 1:h) )

      PARAM%surf(1:w, 1:h) = tab_tmp(1:w, 1:h)

      deallocate( tab_tmp )

   return
   endsubroutine apod_sur


   subroutine apod_acf()
   !================================================================================================
   !<@note Function that apodize the acf to prevent spectral leakage
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h

      real(kind=R8) :: tau1, tau2, ang, coeff, c, s, R1, R2

      real(kind=R8), dimension(1:8) :: ana_res

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      logical(kind=I4) :: set_acf

      read(JOB,*) set_acf ; LINE_READ = LINE_READ +1 ; write(SPY,*) LINE_READ, "Set ACF ", set_acf

      w = PARAM%width
      h = PARAM%height

      call ellipse_acf( tabin = PARAM%acf_surf(1:w, 1:h),            &  ! IN
                         long = w,                                   &  ! IN
                         larg = h,                                   &  ! IN
                        p_acv = ana_res(1:8),                        &  ! OUT -> correlation lengths
                          cut = 0.5_R8,                              &  ! IN  -> z cut plane
                     scale_xy = [PARAM%surf_dx, PARAM%surf_dy],      &  ! IN  -> lags along x and y
                          omp = .true. )                                ! IN  -> use multithread?

      tau1 = ana_res(1)
      tau2 = ana_res(2)
      ang  = ana_res(4) * PI_R8 / 180

      c = cos(ang)
      s = sin(ang)

      ! ellipsis bounding box (half dimensions)
      R1 = sqrt( (c * tau1)**2 + (s * tau2)**2 )
      R2 = sqrt( (s * tau1)**2 + (c * tau2)**2 )

      ! scale factor for the apodization to take place within the surface
      ! 0.4 * image width is less than half width
      coeff = min( 0.4 * PARAM%surf_width / R1, 0.4 * PARAM%surf_height / R2 )

      allocate(  tab_tmp(1:w, 1:h) )

      ! modified Tuckey windowing
      call apod2(  tab_in = PARAM%acf_surf(1:w, 1:h),    &  ! IN
                  tab_out = tab_tmp(1:w, 1:h),           &  ! OUT
                     long = w,                           &  ! IN
                     larg = h,                           &  ! IN
                     tau1 = coeff * tau1 ,               &  ! IN
                     tau2 = coeff * tau2 ,               &  ! IN
                      ang = ang )                           ! IN

      if ( set_acf ) then

         ! acf calculated becomes prescribed acf
         PARAM%imp_acf(1:w, 1:h) = tab_tmp(1:w, 1:h)

      else

         PARAM%acf_surf(1:w, 1:h) = tab_tmp(1:w, 1:h)

      endif

      deallocate( tab_tmp )

   return
   endsubroutine apod_acf


   subroutine repr_img()
   !================================================================================================
   !<@note Function that set parameters for image reproduction
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h, l

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp, surf_tmp

      integer(kind=I4) :: reproduction_step

      real(kind=R8) :: cutoff, dx, dy

      type(SCALE_SURF) :: scale_img_tmp

      type(MOMENT_STAT) :: m_res

      read(JOB,*) reproduction_step ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "reproduction_step ", reproduction_step

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      dx = PARAM%surf_dx
      dy = PARAM%surf_dy

      allocate(  tab_tmp(1:w, 1:h) )
      allocate( surf_tmp(1:w, 1:h) )

      select case (reproduction_step)

         case(0:1)

            PARAM%orig_surf%cut = 0.5_R8
            PARAM%curr_surf%cut = 0.5_R8

            call alloc_tabs()

            PARAM%surf_copy(1:w, 1:h) = PARAM%surf(1:w, 1:h)

            call calc_moments(    tab = PARAM%surf(1:w, 1:h),      &  ! IN
                                   mx = PARAM%m_ini,               &  ! OUT
                               nb_mom = 4 )                           ! IN

            call apod(  tab_in = PARAM%surf(1:w, 1:h),   &  !
                       tab_out = tab_tmp(1:w, 1:h),      &  !
                          long = w,                      &  !
                          larg = h,                      &  !
                      type_apo = 'tuckey' )                 !

            call std_array( tab_tmp(1:w, 1:h) )

            call acf_wiener(  tab_in = tab_tmp(1:w, 1:h),             &  ! IN
                             tab_out = PARAM%acf_surf(1:w, 1:h),      &  ! OUT
                                   w = w,                             &  ! IN
                                   h = h  )                              ! IN

            call surface_analysis( app = 0 )

            !-------------- NORMALIZATION ------------------------

            call std_array( PARAM%surf(1:w, 1:h), PARAM%m_ini )

         case(2)

            !------------ SAVE LOW FREQ SURF ---------------------

            PARAM%surf_LF(1:w, 1:h) = PARAM%surf(1:w, 1:h)

            !------------ REPRODUCE HIGH FREQ ---------------------

            PARAM%surf(1:w, 1:h) = PARAM%surf_HF(1:w, 1:h)

            PARAM%m_end = PARAM%m__HF

            PARAM%vect_h(1:l) = reshape( PARAM%surf(1:w, 1:h), [l] )

            call sort_array2(tab_inout = PARAM%vect_h(1:l), n = l)

            surf_tmp(1:w, 1:h) = PARAM%surf(1:w, 1:h)

         case(3)

            !------------ SAVE HIGH FREQ SURF ---------------------

            PARAM%surf_HF(1:w, 1:h) = PARAM%surf(1:w, 1:h)

            PARAM%surf(1:w, 1:h) = PARAM%surf_LF(1:w, 1:h) + PARAM%surf_HF(1:w, 1:h) * PARAM%m__HF%si / PARAM%m__LF%si

            call std_array( PARAM%surf(1:w, 1:h) )

            PARAM%vect_h(1:l) = reshape( PARAM%surf_copy(1:w, 1:h), [l] )
            ! heights are sorted
            call sort_array2(tab_inout = PARAM%vect_h(1:l), n = l)

            call calc_moments(   tab = PARAM%vect_h(1:l),      &  ! IN
                                  mx = m_res,                  &  ! OUT
                              nb_mom = 2 )                        ! IN

            PARAM%vect_h(1:l) = ( PARAM%vect_h(1:l) - m_res%mu ) / m_res%si

            surf_tmp(1:w, 1:h) = PARAM%surf_copy(1:w, 1:h)

            call std_array( surf_tmp(1:w, 1:h) )

         case(4)

            call std_array( PARAM%surf(1:w, 1:h) )

            PARAM%surf(1:w, 1:h) = PARAM%surf(1:w, 1:h) * PARAM%m_ini%si + PARAM%m_ini%mu

            call surface_analysis( app = 1 )

      endselect

      select case (reproduction_step)

         case(0)

            PARAM%imp_acf(1:w, 1:h) = PARAM%acf_surf(1:w, 1:h)

            ! the calculated moments become the prescribed ones
            PARAM%m_end = PARAM%m_ini

            ! heights are stored because they are the prescribed ones
            PARAM%vect_h(1:l) = reshape( PARAM%surf(1:w, 1:h), [l] )

            ! shuffle the set, then ...
            call scramble( tab = PARAM%vect_h(1:l),   &  ! INOUT
                            lg = l )                     ! IN

            ! ... define an initial random surface ...
            PARAM%surf(1:w, 1:h) = reshape( PARAM%vect_h(1:l), [w, h] )

            ! and sort heights
            call sort_array2(tab_inout = PARAM%vect_h(1:l), n = l)

         case(1)

            !------------ SUBTRACT LOW FREQ ----------------------

            read(JOB,*) cutoff ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "cutoff ", cutoff

            call fft_filter(       tab = PARAM%surf(1:w, 1:h),    &  ! IN
                                  long = w,                       &  ! IN
                                  larg = h,                       &  ! IN
                                cutoff = cutoff,                  &  ! IN
                                bf_tab = PARAM%surf_LF(1:w, 1:h), &  ! OUT
                             multi_fft = .false.,                 &  ! IN
                                   ext = 'constant' )                ! IN

            PARAM%surf_HF(1:w, 1:h) = PARAM%surf(1:w, 1:h) - PARAM%surf_LF(1:w, 1:h)

            !------------ MOMENTS LF & HF  ----------------------
            call std_array( tab = PARAM%surf_HF(1:w, 1:h), mx = PARAM%m__HF )

            call std_array( tab = PARAM%surf_LF(1:w, 1:h), mx = PARAM%m__LF )

            !------------ PRINT LF & HF SURFS ----------------------
            call init_scal( scal = scale_img_tmp,      &           ! out; creates a surface type, containing ...
                              nx = w,                  &           !  in; ... the number of points along x ...
                              ny = h,                  &           !  in; ... the number of points along y ...
                              lx = PARAM%surf_width,   &           !  in; ... the length (default unit : m) ...
                              ly = PARAM%surf_height,  &           !  in; ... the width ...
                          unit_z = 'm'        )                    !  in; ... and the unit along z.

            call write_surf( nom_fic = "out/BF.sur",            &  !
                               tab_s = PARAM%surf_LF(1:w, 1:h), &  !
                                scal = scale_img_tmp )             !

            call write_surf( nom_fic = "out/HF.sur",            &  !
                               tab_s = PARAM%surf_HF(1:w, 1:h), &  !
                                scal = scale_img_tmp )             !

            !------------ REPRODUCE LOW FREQ ----------------------

            PARAM%surf(1:w, 1:h) = PARAM%surf_LF(1:w, 1:h)

            PARAM%m_end = PARAM%m__LF

            PARAM%vect_h(1:l) = reshape( PARAM%surf(1:w, 1:h), [l] )

            call sort_array2(tab_inout = PARAM%vect_h(1:l), n = l)

            surf_tmp(1:w, 1:h) = PARAM%surf(1:w, 1:h)

         case default

            continue

      endselect

      if ( reproduction_step > 0 .and. reproduction_step < 4 ) then

         !------------ DESIRED ACF: PARAM%imp_acf -------------

         call apod(  tab_in = surf_tmp(1:w, 1:h),     &  !
                    tab_out = tab_tmp(1:w, 1:h),      &  !
                       long = w,                      &  !
                       larg = h,                      &  !
                   type_apo = 'tuckey' )                 !

         call std_array( tab_tmp(1:w, 1:h) )

         call acf_wiener(  tab_in = tab_tmp(1:w, 1:h),             &  ! IN
                          tab_out = PARAM%imp_acf(1:w, 1:h),       &  ! OUT
                                w = w,                             &  ! IN
                                h = h  )                              ! IN

      endif

      deallocate(  tab_tmp )
      deallocate( surf_tmp )

   return
   endsubroutine repr_img


   subroutine calc_z_f()
   !================================================================================================
   !<@note Function that returns PARAM%surf, the surface made of heights with the required statistical
   !< moments, in the right order.
   !<
   !< - The heights come from the vector PARAM%vect_h
   !< - the heights order is stored in the vector PARAM%order
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: i
      integer(kind=I4) :: w, h, l

      real(kind=R8), allocatable, dimension(:) :: tab_tmp

      read(JOB,*) PARAM%calc_zf ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "PARAM%calc_zf ", PARAM%calc_zf

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      ! final set of heights are generated to meet the desired statistical moments
      ! It is done once.
      if ( PARAM%calc_zf ) then

         write(SPY,*) 'calc_z_f -> final set of heights are generated to meet the desired statistical moments'

         call build_heights(     vec_out = PARAM%vect_h(1:l),                                     &  ! OUT
                            use_fct_expo = ( PARAM%m_end%ku < 1.34 * PARAM%m_end%sk**2 + 1.8),    &  ! IN
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

      call std_array( tab = tab_tmp(1:l) )

      deallocate( tab_tmp )

   return
   endsubroutine calc_z_f


   subroutine calc_ord()
   !================================================================================================
   !<@note Function that returns the vector PARAM%order that contains the heights order.
   !<
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
   !<
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

      call std_array( tab = PARAM%surf(1:w, 1:h), mx = m_res )

      write(SPY,*) 'sk ku fin ', m_res%sk, m_res%ku

      deallocate( cmple )
      deallocate( ftab )

   return
   endsubroutine digi_fil


   subroutine add_nois()
   !================================================================================================
   !<@note Function that returns the starting surface of random heights
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer (kind=I4) :: w, h, l

      real(kind=R8) :: size_noise

      real(kind=R8), allocatable, dimension(:) :: tab_tmp

      type(MOMENT_STAT) :: mx

      read(JOB,*) mx%sk, mx%ku, size_noise ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "ssk, sku, size ", mx%sk, mx%ku, size_noise

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      allocate( tab_tmp(1:l) )

      call build_heights(     vec_out = tab_tmp(1:l),                      &  ! OUT
                         use_fct_expo = ( mx%ku < 1.34 * mx%sk**2 + 1.8 ), &  ! IN
                             stats_in = mx,                                &  ! IN
                                   lg = l )                                   ! IN

      call std_array( tab = tab_tmp(1:l) )

      call scramble( tab = tab_tmp(1:l),   &  ! INOUT
                      lg = l )                ! IN

      PARAM%surf(1:w, 1:h) = PARAM%surf(1:w, 1:h) + size_noise * reshape( tab_tmp(1:l), [w, h] )

      deallocate( tab_tmp )


   return
   endsubroutine add_nois


   subroutine calc_z_i()
   !================================================================================================
   !<@note Function that returns the starting surface of random heights
   !<
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
                         use_fct_expo = ( PARAM%m_stt%ku < 1.34 * PARAM%m_stt%sk**2 + 1.8 ), &  ! IN
                             stats_in = PARAM%m_stt,                                         &  ! IN
                                   lg = l )                                                     ! IN

      call std_array( tab = tab_tmp(1:l), mx = m_res )

      write(TER,*) 'starting statistical moments ', m_res%mu, m_res%va, m_res%sk, m_res%ku

      call scramble( tab = tab_tmp(1:l),   &  ! INOUT
                      lg = l )                ! IN


      PARAM%surf(1:w, 1:h) = reshape( tab_tmp(1:l), [w, h] )

      deallocate( tab_tmp )

   return
   endsubroutine calc_z_i


   subroutine calc_ffh()
   !================================================================================================
   !<@note Function that returns ...
   !<
   !< - the digital filter to apply to the height distribution \( \text{PARAM%fhi} = \sqrt{ \left| FFT(\text{PARAM%imp_acf}) \right| } \)
   !< - the starting statistical moments PARAM%m_stt%sk, PARAM%m_stt%ku
   !< - whether the exponential function will be used, PARAM%reajust_skku
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h, l

      complex(kind=R8), dimension(:,:), allocatable :: cmple
      real   (kind=R8), dimension(:,:), allocatable :: tab_tmp

      type(MOMENT_STAT) :: m_h

      read(JOB,*) PARAM%calc_mstt ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "PARAM%calc_mstt ", PARAM%calc_mstt

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

      if ( PARAM%calc_mstt ) then

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
   !<
   !< If the surface to generate is non periodic, the starting surface is extended. The final surface
   !< will be a part of it. Indeed the extended surface will be periodic, because the use of FFTs.
   !<
   !< If a roughness orientation is chosen, in addition with long correlation lengths, a windowing
   !< should be applied to the acf to prevent from artifacts (vertical and horizontal lines)
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: w, h
      logical(kind=I4) :: with_apod

      real(kind=R8) :: ratio, a, b, c, s, lx, ly

      real(kind=R8), dimension(1:8) :: res

      read(JOB,*) PARAM%curr_surf%cl1, PARAM%curr_surf%cl2, PARAM%curr_surf%cut ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "l_acf ", PARAM%curr_surf%cl1, PARAM%curr_surf%cl2, PARAM%curr_surf%cut

      if (PARAM%curr_surf%cl1 < PARAM%curr_surf%cl2) stop "inverser cl1, cl2"

      read(JOB,*) PARAM%curr_surf%ang  ; LINE_READ = LINE_READ +1 ; write(SPY,*) "line: ", LINE_READ, "ang ", PARAM%curr_surf%ang

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

         a = PARAM%curr_surf%cl1
         b = PARAM%curr_surf%cl2

         c = cos( PARAM%curr_surf%ang * PI_R8 / 180 )
         s = sin( PARAM%curr_surf%ang * PI_R8 / 180 )

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

      ratio = PARAM%curr_surf%cl2 / PARAM%curr_surf%cl1

      ! acf_theo is the theoretical acf for the normal surface
      call calc_imp_acf( long = w,                                &  ! IN
                         larg = h,                                &  ! IN
                         apod = with_apod,                        &  ! IN
                         tau1 = PARAM%curr_surf%cl1,              &  ! IN
                         tau2 = PARAM%curr_surf%cl2,              &  ! IN
                        alpha = log(PARAM%curr_surf%cut),         &  ! IN
                          ang = PARAM%curr_surf%ang * PI_R8/180,  &  ! IN
                      tab_acf = PARAM%imp_acf(1:w, 1:h))             ! OUT

      call ellipse_acf( tabin = PARAM%imp_acf(1:w, 1:h),          &  ! IN
                         long = w,                                &  ! IN
                         larg = h,                                &  ! IN
                        p_acv = res(1:8),                         &  ! OUT
                          cut = PARAM%curr_surf%cut,              &  ! IN
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

      close(STA)

      if ( allocated( PARAM%imp_acf   ) ) deallocate( PARAM%imp_acf   )
      if ( allocated( PARAM%acf_surf  ) ) deallocate( PARAM%acf_surf  )
      if ( allocated( PARAM%surf_LF   ) ) deallocate( PARAM%surf_LF   )
      if ( allocated( PARAM%surf_HF   ) ) deallocate( PARAM%surf_HF   )
      if ( allocated( PARAM%surf_msk  ) ) deallocate( PARAM%surf_msk  )
      if ( allocated( PARAM%surf      ) ) deallocate( PARAM%surf      )
      if ( allocated( PARAM%surf_copy ) ) deallocate( PARAM%surf_copy )
      if ( allocated( PARAM%vect_h    ) ) deallocate( PARAM%vect_h    )
      if ( allocated( PARAM%fhi       ) ) deallocate( PARAM%fhi       )
      if ( allocated( PARAM%order     ) ) deallocate( PARAM%order     )

      call end_fftw3()

   return
   endsubroutine end_scri


   subroutine make_tex()
   !================================================================================================
   !< @note Function that creates a periodic macro-texture: knowing the FFT of an analytical texture
   !        a surface is created and then transformed by a FFT. As a result it is periodic.
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      real(kind=R8), allocatable, dimension(:)   :: xloc, yloc, haut

      real   (kind=R8) :: ca, sa, a, b, c, x0, y0, h0
      integer(kind=I4) :: kk, w2, h2

      complex  (kind=R8) :: eix, eiy, ftmp
      real     (kind=R8) :: x1, y1, ax, ay, sx, sy, coeff
      integer  (kind=I4) :: ksi, eta

      complex(kind=R8), allocatable, dimension(:,:) :: cmpl_in, cmpl_ou, ft

      real(kind=R8) :: alpha                 ! texture orientation
      real(kind=R8) :: sigmax, sigmay        ! texture stdv

      integer(kind=I4) :: texture            ! number of texture elements

      character(len=016) :: dimples          ! texture elementar shapes: "circle" or "square"

      real(kind=R8) :: dx, dy

      integer(kind=I4) :: w, h, l

      type(MOMENT_STAT) :: m_res

      w = PARAM%width
      h = PARAM%height
      l = PARAM%npts

      dx = PARAM%surf_dx
      dy = PARAM%surf_dy

      read(JOB,*) dimples ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, "dimples ", dimples
      read(JOB,*) texture ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, "texture ", texture

      sigmax = PARAM%curr_surf%cl1 / dx / w * 2
      sigmay = PARAM%curr_surf%cl2 / dx / h * 2
      alpha  = PARAM%curr_surf%ang

      w2 = w / 2
      h2 = h / 2
      if ( w == 2 * (w / 2) + 1 ) w2 = w / 2 +1
      if ( h == 2 * (h / 2) + 1 ) h2 = h / 2 +1

      allocate( cmpl_in(1:w, 1:h),     &  !
                cmpl_ou(1:w, 1:h) )       !

      allocate( ft(-w/2:w2 -1,         &  !
                   -h/2:h2 -1) )          !

      allocate( xloc(1:TEXTURE),       &  !
                yloc(1:TEXTURE),       &  !
                haut(1:TEXTURE) )         !


      ca = cos(alpha * PI_R8 / 180)
      sa = sin(alpha * PI_R8 / 180)

      ! Gaussian coefficients
      a = 0.5*( (ca / sigmax) * (ca / sigmax) + (sa / sigmay) * (sa / sigmay) )
      c = 0.5*( (sa / sigmax) * (sa / sigmax) + (ca / sigmay) * (ca / sigmay) )
      b = 1.0*( (sa / sigmax) * (ca / sigmax) - (sa / sigmay) * (ca / sigmay) )

      call random_number( xloc(1:TEXTURE ) ) ! localisation in [0,1]x[0,1]
      call random_number( yloc(1:TEXTURE ) )

      call random_number( haut(1:TEXTURE ) )

      haut(1:TEXTURE ) = 2 * haut(1:TEXTURE ) - 1._R8

      coeff = 2 * pi_R8

      do ksi = -w / 2, w2 - 1
      do eta = -h / 2, h2 - 1

         x1 = (+ksi * ca + eta * sa) * coeff
         y1 = (-ksi * sa + eta * ca) * coeff

         ftmp = 0
         sx = sigmax
         sy = sigmay

         do kk = 1, TEXTURE ! nombre de pics sur le niveau niv

            x0 = xloc(kk)
            y0 = yloc(kk)
            h0 = haut(kk)   ! d'un niveau à l'autre, on atténue les pics d'un facteur "per_tex"

            ax = -x1 * x0
            ay = -y1 * y0

            eix = cmplx(cos(ax), sin(ax), kind=R8)
            eiy = cmplx(cos(ay), sin(ay), kind=R8)

            select case( DIMPLES(1:6) )

               case("square")
                  ftmp = ftmp + 2 * h0 * pi_R8 * sx * sy * sinc( sx * x1 ) * eix &  !
                                                         * sinc( sy * y1 ) * eiy
               case("circle")
                  ftmp = ftmp + 2 * h0 * pi_R8 * sx * sy * exp( -0.5 * (sx**2) * (x1**2) ) * eix &  !
                                                         * exp( -0.5 * (sy**2) * (y1**2) ) * eiy
               case default
                  stop "make_terra, no dimple type selected"

            endselect

         enddo

         ft(ksi, eta) = ftmp

      enddo
      enddo

      cmpl_in(1:w2,      1:h2) = ft(0:w2 - 1,      0:h2 - 1)   ! SO = NE
      cmpl_in(1:w2, h2 + 1:h ) = ft(0:w2 - 1, -h / 2:   -1 )   ! SE = NO

      cmpl_in(w2 + 1:w,     1:h2) = ft(-w/2:-1,    0:h2 -1)    ! NO = SE
      cmpl_in(w2 + 1:w, h2 +1:h ) = ft(-w/2:-1, -h/2:   -1)    ! NE = SO

      FFT_DIM = 0 ! forces fftw desalloc and reinit

      call calc_fftw3(sens = BACKWARD,                &  !
                    tab_in = cmpl_in(1:w, 1:h),       &  !
                    tab_ou = cmpl_ou(1:w, 1:h),       &  !
                      long = w,                       &  !
                      larg = h)                          !

      FFT_DIM = 0 ! forces fftw desalloc and reinit

      PARAM%surf(1:w, 1:h) = real(cmpl_ou(1:w, 1:h), kind=R8) * sqrt(real(w * h, kind=R8))

      call calc_moments(    tab = PARAM%surf(1:w, 1:h),      &  ! IN
                             mx = m_res,                     &  ! OUT
                         nb_mom = 3 )                           ! IN

      if ( m_res%sk > 0. ) PARAM%surf(1:w, 1:h) = - PARAM%surf(1:w, 1:h)

      deallocate( xloc, yloc, haut )

      deallocate( cmpl_in, cmpl_ou )

   contains
      !-----------------------------------------
      function sinc(x)
      implicit none
      real(kind=R8) :: sinc
      real(kind=R8), intent(in) :: x
         if ( abs(x)>1.0e-8) then
            sinc = sin(x)/x
         else
            sinc = UN
         endif
      return
      endfunction sinc
      !-----------------------------------------
   endsubroutine make_tex


   subroutine make_scratches()
   !================================================================================================
   !< @note This subroutine initializes a real matrix `tab` of dimensions `nx` by `ny` with ones
   !        and adds `sn` linear scratches with value 0. Each scratch has a random length up to `sl`,
   !        a constant width `sw`, and a random orientation. Scratches can touch the matrix boundaries.
   !
   !  @endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: nx   ! Number of rows
      integer(kind=I4) :: ny   ! Number of columns
      integer(kind=I4) :: sn   ! Number of scratches
      integer(kind=I4) :: sw   ! Scratch width
      integer(kind=I4) :: sl   ! Maximum scratch length

      integer(kind=I4) :: i, j, n
      real(kind=R8)    :: theta, length, depth, x0, y0, x1, y1
      real(kind=R8)    :: rand, r1, r2

      nx = PARAM%width
      ny = PARAM%height

      read(JOB,*)  sn, sw, sl ; LINE_READ = LINE_READ + 1 ; write(SPY,*) "line: ", LINE_READ, "sn, sw, sl  ", sn, sw, sl

      ! Initialize matrix with ones
      PARAM%surf(1:nx, 1:ny) = 1.0_R8

      ! Seed random number generator
      call random_seed()

      ! Generate each scratch
      do n = 1, sn

         ! Random starting point
         call random_number(rand)
         x0 = rand * nx
         call random_number(rand)
         y0 = rand * ny

         ! Random depth
         ! Génération de deux variables gamma
         call gamma_random(shape = 3._R8, output = r1)
         call gamma_random(shape = 1._R8, output = r2)
         ! Transformation pour obtenir la loi bêta
         depth = r1 / (r1 + r2)

         ! Random orientation (angle in radians)
         call random_number(rand)
         theta = rand * 2.0 * PI_R8

         ! Random length up to sl
         call random_number(rand)
         length = ( 2 * (rand - 0.5_R8) / 10 + 0.9_R8 ) * sl

         ! Calculate end point of the scratch
         x1 = x0 + length * cos(theta)
         y1 = y0 + length * sin(theta)

         ! Draw the scratch with width sw
         do i = max(1, floor(min(x0, x1)) - sw), min(nx, ceiling(max(x0, x1)) + sw)

            do j = max(1, floor(min(y0, y1)) - sw), min(ny, ceiling(max(y0, y1)) + sw)

               ! Check if point (i,j) is within sw/2 distance from the line segment
               if ( point_to_line_distance(real(i, kind=R8), real(j, kind=R8), x0, y0, x1, y1) <= sw/2.0 ) then

                  PARAM%surf(i, j) = depth

               endif

            enddo

         enddo

      enddo

   contains

      subroutine gamma_random(shape, output)
      implicit none
      real(kind=R8), intent( in) :: shape
      real(kind=R8), intent(out) :: output

         real(kind=R8) :: b, c
         real(kind=R8) :: u, v, w, y

         ! algorithme de marsaglia pour générer une variable gamma
         b = shape - 1.0/3.0
         c = 1.0 / sqrt(9.0 * b)

         do

            call random_number(u)
            call random_number(v)

            v = 1.0 + c * log(v) / sqrt(u)

            if (v <= 0.0) cycle

            w = v * v * v

            call random_number(y)

            if (y < 1.0 - 0.0331 * (log(v) ** 4)) exit

            if (log(y) < 0.5 * log(v) * log(v) + b * (1.0 - v + log(v))) exit

         enddo

         output = b * w

      return
      endsubroutine gamma_random
      !-----------------------------------------
      real(kind=R8) function point_to_line_distance(px, py, x0, y0, x1, y1)
         !! Calculate the shortest distance from a point to a line segment.
      implicit none
      real(kind=R8), intent(in) :: px  !! x-coordinate of the point
      real(kind=R8), intent(in) :: py  !! y-coordinate of the point
      real(kind=R8), intent(in) :: x0  !! x-coordinate of the line segment start
      real(kind=R8), intent(in) :: y0  !! y-coordinate of the line segment start
      real(kind=R8), intent(in) :: x1  !! x-coordinate of the line segment end
      real(kind=R8), intent(in) :: y1  !! y-coordinate of the line segment end

         real(kind=R8) :: dx, dy, t, proj_x, proj_y

         dx = x1 - x0
         dy = y1 - y0

         ! Handle case where start and end points are the same
         if ( (dx**2 + dy**2) < 1.e-12 ) then

            point_to_line_distance = sqrt((px - x0)**2 + (py - y0)**2)

            return

         endif

         ! Project point onto the line
         t = max(0.0_R8, min(1.0_R8, ((px - x0)*dx + (py - y0)*dy)/(dx*dx + dy*dy)))

         proj_x = x0 + t * dx
         proj_y = y0 + t * dy

         ! Calculate distance
         point_to_line_distance = sqrt( (px - proj_x)**2 + (py - proj_y)**2 )

      return
      endfunction point_to_line_distance
      !-----------------------------------------

   endsubroutine make_scratches


endmodule script

