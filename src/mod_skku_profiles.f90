!< author: Arthur Francisco
!  version: 1.0.0
!  date: may, 03 2019
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.2em;">
!        **Routines to generate a (Sk, Ku) series**
!  </span>

module skku_profiles
use data_arch,   only : I4, R8, UN, EPS_R8, PI_R8, HIG_E8
use crest_param, only : PARAM, JOB, SPY, TER, FCT_TANG, FCT_EXPO
use stat_mom,    only : moment_stat, calc_moments, scramble, std_array
use sort_arrays, only : sort_array2, init_order
use pikaia_oop,  only : pikaia_class
implicit none

!> {!src/inc_doc/profile_generation.md!}

private

public :: build_heights


contains

   subroutine build_heights(vec_out, use_fct_expo, stats_in, lg)
   !================================================================================================
   !<@note Function that returns a set of heights that matches desired statistical moments.
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer (kind=I4), intent(in )                    :: lg              !! *length of the height vector*
   type(MOMENT_STAT), intent(in )                    :: stats_in        !! *input statistical moments*
   logical (kind=I4), intent(in )                    :: use_fct_expo    !! *should exponential function rather than tangent function be used?*
   real    (kind=R8), intent(out), dimension(1:lg)   :: vec_out         !! *height vector*

      integer(kind=I4) :: istat
      real(kind=R8)    :: cost_val

      type(MOMENT_STAT) :: m_tmp

      real(kind=R8), dimension(:), allocatable :: xlower
      real(kind=R8), dimension(:), allocatable :: xupper
      real(kind=R8), dimension(:), allocatable :: xresul

      ! put input parameters in global variables, so that they can be used in the function "fitness_skku_anal"
      PARAM%m_inp%sk = stats_in%sk
      PARAM%m_inp%ku = stats_in%ku

      ! if the Pearson limit is to close to the point (Ssk, Sku), an exponential function is used
      if ( use_fct_expo ) then

         PARAM%func_gen = FCT_EXPO
         PARAM%nparam   = 3

      else

         PARAM%func_gen = FCT_TANG
         PARAM%nparam   = 2

      endif

      ! Genetic algorithm is used to determinate the tangent parameters \alpha and \beta so that, the set of lg heights
      ! will match the statistical moments.
      !..............................................................................

      ! initialization
      allocate( xresul(1:PARAM%nparam) ) ; xresul = 0.0_R8
      allocate( xlower(1:PARAM%nparam) ) ; xlower = 1.e-6_R8
      allocate( xupper(1:PARAM%nparam) ) ; xupper = 1.0_R8

      call pikaia_skku_solver( pik_class = PARAM%pik_class,          &  ! INOUT
                                    step = 'init',                   &  ! IN
                                      xl = xlower(1:PARAM%nparam),   &  ! IN
                                      xu = xupper(1:PARAM%nparam),   &  ! IN
                                      xx = xresul(1:PARAM%nparam),   &  ! IN
                                  nparam = PARAM%nparam,             &  ! IN
                                    cost = cost_func_skku,           &  ! IN
                                   istat = istat,                    &  ! OUT
                                       f = cost_val )                   ! IN

      call pikaia_skku_solver( pik_class = PARAM%pik_class,          &  ! INOUT
                                    step = 'solv',                   &  ! IN
                                      xl = xlower(1:PARAM%nparam),   &  ! IN
                                      xu = xupper(1:PARAM%nparam),   &  ! IN
                                      xx = xresul(1:PARAM%nparam),   &  ! OUT
                                  nparam = PARAM%nparam,             &  ! IN
                                    cost = cost_func_skku,           &  ! IN
                                   istat = istat,                    &  ! OUT
                                       f = cost_val )                   ! IN
      !..............................................................................

      ! the parameters habe been found, let generate lg heights
      !..............................................................................
      call profil_theo_trie_1D( tab = vec_out(1:lg),               &  ! OUT
                                 lg = lg,                          &  ! IN
                                  x = xresul(1:PARAM%nparam),      &  ! IN
                                 mx = m_tmp )                         ! OUT
      !..............................................................................

      ! PARAM%func_gen value is retrieved
      !PARAM%func_gen = fct_sav

      deallocate( xresul )
      deallocate( xlower )
      deallocate( xupper )

      call std_array( tab = vec_out(1:lg), mx = m_tmp)

      ! the parameter found can lead to inverted heights
      if (stats_in%sk * m_tmp%sk < 0.) then

         vec_out(1:lg) = -vec_out(1:lg)

      endif

      ! heights are sorted
      call sort_array2(tab_inout = vec_out(1:lg), n = lg)

   return
   endsubroutine build_heights


   subroutine cost_func_skku(me, x, f)
   !! Quantify de distance between desired moments and calculated moments
   implicit none
   class(pikaia_class), intent(inout)               :: me
   real(kind=R8)      , intent(in   ), dimension(:) :: x
   real(kind=R8)      , intent(  out)               :: f

      f = fitness_skku_anal(n = PARAM%nparam, x = x(1:PARAM%nparam) )

   return
   endsubroutine cost_func_skku


   real(kind=R8) function fitness_skku_anal(n, x)
   !================================================================================================
   !<@note Generic cost function: difference between the imposed statistical moments and those
   !< obtained. The optimization problem must be turned into a maximization problem (as often
   !< in the optimization routines).
   !<
   !< The closer cost to 100 the better series.
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer(kind=I4), intent(in)                 :: n !! *number of unknowns*
   real   (kind=R8), intent(in), dimension(1:n) :: x !! *vector of unknowns*

      real(kind=R8) :: sk, ku

      select case (PARAM%func_gen)
         case (FCT_TANG) ; call calculs_skku_tan (bounds = x(1:PARAM%nparam), lg = PARAM%npts, ssk = sk, sku = ku)
         case (FCT_EXPO) ; call calculs_skku_exp3(bounds = x(1:PARAM%nparam), lg = PARAM%npts, ssk = sk, sku = ku)
      endselect

      fitness_skku_anal = ( abs(sk**2 - PARAM%m_inp%sk**2) +                     &  !
                            abs(ku    - PARAM%m_inp%ku   ) )/PARAM%m_inp%ku         !

      fitness_skku_anal = 100._r8 / (1._r8 + fitness_skku_anal)

   return
   endfunction fitness_skku_anal


   subroutine pikaia_skku_solver(pik_class, step, xl, xu, nparam, cost, istat, f, xx)
   !================================================================================================
   !<@note This is a refactoring of the PIKAIA unconstrained optimization code from the High Altitude Observatory.
   !< The original code is public domain and was written by Paul Charbonneau & Barry Knapp.
   !<
   !< The present code is the awesome modern Fortran version written by Jabob Williams:
   !<
   !< [OOP Pikaia, Jacob Williams](https://github.com/jacobwilliams/pikaia)
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   type(pikaia_class), intent(inout) :: pik_class                       !! **PIKAIA** *class instanciation*
   character(len=4),   intent(in   ) :: step                            !! *init* or *solv*
   integer(kind=I4),   intent(in   ) :: nparam                          !! *number of parameters*
   real(kind=R8),      intent(in   ), dimension(1:nparam) :: xl         !! *lower bonds of xx*
   real(kind=R8),      intent(in   ), dimension(1:nparam) :: xu         !! *upper bonds of xx*
   real(kind=R8),      intent(  out), dimension(1:nparam) :: xx         !! *chromosom for* **PIKAIA**
   integer(kind=I4),   intent(  out) :: istat
   real(kind=R8),      intent(  out) :: f

   interface
      subroutine cost(me, x, f)
      use data_arch,   only : R8
      use pikaia_oop,  only : pikaia_class
      implicit none
      class(pikaia_class), intent(inout)               :: me
      real(kind=R8)      , intent(in   ), dimension(:) :: x
      real(kind=R8)      , intent(  out)               :: f
      endsubroutine cost
   endinterface

      select case(step)

         case('init')

            !initialize the class:
            call pik_class%init(   n = nparam,              &  ! IN           ; the parameter space dimension, i.e., the number of adjustable parameters (size of the x vector).
                                  xl = xl,                  &  ! IN, DIM(n)   ;  vector of lower bounds for x
                                  xu = xu,                  &  ! IN, DIM(n)   ;  vector of upper bounds for x
                                   f = cost,                &  !              ; user-supplied scalar function of n variables, which must have the pikaia_func procedure interface.
                              status = istat ,              &  ! OUT          ; status output flag (0 if there were no errors)
                             !iter_f = report_iteration,    &  !     OPT      ; user-supplied subroutine that will report the best solution for each generation. It must have the iter_func procedure interface.
                                  np = 100,                 &  ! IN, OPT      ; number of individuals in a population (default is 100)
                                ngen = 1000,                &  ! IN, OPT      ; maximum number of iterations
                                  nd = 9,                   &  ! IN           ; number of significant digits (i.e., number of genes) retained in chromosomal encoding
                              pcross = 0.85_R8,             &  ! IN, OPT      ; crossover probability; must be <= 1.0 (default is 0.85). If crossover takes place, either one or two splicing points are used, with equal probabilities
                              pmutmn = 0.0005_R8,           &  ! IN, OPT      ; minimum mutation rate; must be >= 0.0 (default is 0.0005)
                              pmutmx = 0.25_R8,             &  ! IN, OPT      ; maximum mutation rate; must be <= 1.0 (default is 0.25)
                                pmut = 0.005_R8,            &  ! IN, OPT      ; initial mutation rate; should be small (default is 0.005) (Note: the mutation rate is the probability that any one gene locus will mutate in any one generation.)
                                imut = 2,                   &  ! IN, OPT      ; mutation mode; 1/2/3/4/5 (default is 2).
                                                               !              1=one-point mutation, fixed rate.
                                                               !              2=one-point, adjustable rate based on fitness.
                                                               !              3=one-point, adjustable rate based on distance.
                                                               !              4=one-point+creep, fixed rate.
                                                               !              5=one-point+creep, adjustable rate based on fitness.
                                                               !              6=one-point+creep, adjustable rate based on distance.
                                fdif = 1._R8,               &  ! IN, OPT      ; relative fitness differential; range from 0 (none) to 1 (maximum). (default is 1.0)
                                irep = 3,                   &  ! IN, OPT      ; reproduction plan; 1/2/3=Full generational replacement/Steady-state-replace-random/Steady- state-replace-worst (default is 3)
                              ielite = 0,                   &  ! IN, OPT      ; elitism flag; 0/1=off/on (default is 0) (Applies only to reproduction plans 1 and 2)
                                ivrb = 0,                   &  ! IN, OPT      ; printed output 0/1/2=None/Minimal/Verbose
                     convergence_tol = 1.0e-6_R8,           &  ! IN, OPT      ; convergence tolerance; must be > 0.0 (default is 0.0001)
                  convergence_window = 200,                 &  ! IN, OPT      ; convergence window; must be >= 0 This is the number of consecutive solutions within the tolerance for convergence to be declared (default is 20)
                  initial_guess_frac = 0.1_R8,              &  ! IN, OPT      ; raction of the initial population to set equal to the initial guess. Range from 0 (none) to 1.0 (all). (default is 0.1 or 10%).
                               iseed = 999)                    ! IN, OPT      ; random seed value; must be > 0 (default is 999)

         case('solv')

            call pik_class%solve( x = xx(1:nparam),      &  ! INOUT, DIM(*) ;
                                  f = f,                 &  !   OUT         ;
                             status = istat,             &  !   OUT         ;
                                omp = .true. )              ! IN, OPTIONAL

         case default

            stop 'Wrong choice in "pikaia_skku_solver"'

      endselect

   return
   endsubroutine pikaia_skku_solver


   subroutine calculs_skku_tan(bounds, lg, ssk, sku)
   !================================================================================================
   !<@note Function to calculate the skewness and kurtosis of a **tangent** series
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real   (kind=R8), intent(in), dimension(1:2) :: bounds  !! *defines the function limits* [-pi/2.(1-bounds(1)), +pi/2.(1-bounds(2)]
   integer(kind=I4), intent(in)                 :: lg      !! *vec size*
   real   (kind=R8), intent(out)                :: ssk     !! *theoretical Ssk*
   real   (kind=R8), intent(out)                :: sku     !! *theoretical Sku*

      real(kind=R8)    :: xa, xb, mu, si, sk, ku, a, b
      real(kind=R8)    :: h, hh, b1, b2, alp, bet
      integer(kind=I4) :: i, ia, ib, deb, fin

      !---------------------------------------------
      ! WXMAXIMA file
      !---------------------------------------------
      !kill(all);
      !
      !f11(x):=tan(x)$
      !
      !assume(u<0.)$
      !assume(u>-%pi/2)$
      !assume(v>0.)$
      !assume(v<%pi/2)$
      !
      !I11:integrate(f11(x),x,u,v)$
      !I11:subst(-%pi/2+%pi/2*xa,u,I11)$
      !I11:subst(+%pi/2-%pi/2*xb,v,I11)$
      !I11:expand(trigsimp(I11));
      !
      !f21(x):=f11(x)-mu$
      !I21:integrate(expand(f21(x)^2),x,u,v)$
      !I21:subst(-%pi/2+%pi/2*xa,u,I21)$
      !I21:subst(+%pi/2-%pi/2*xb,v,I21)$
      !I21:expand(trigsimp(I21));
      !
      !f31(x):=f21(x)/si$
      !I31:integrate(f31(x)^3,x,u,v)$
      !I31:subst(-%pi/2+%pi/2*xa,u,I31)$
      !I31:subst(+%pi/2-%pi/2*xb,v,I31)$
      !I31:expand(trigsimp(I31));
      !
      !I41:integrate(f31(x)^4,x,u,v)$
      !I41:subst(-%pi/2+%pi/2*xa,u,I41)$
      !I41:subst(+%pi/2-%pi/2*xb,v,I41)$
      !I41:expand(trigsimp(I41));
      !---------------------------------------------

      ia = 256 ! ia and ib define the interval edges to be excluded ...
      ib = 256 ! ... because of high variations of the function.

      deb = 1  +ia ! start
      fin = lg -ib ! end

      a = bounds(1)
      b = bounds(2)

      hh = (2._R8-a-b)/(lg-1)
      h  = (PI_R8/2)*hh

      xa = a +ia*hh
      xb = b +ib*hh

      b1 = -PI_R8/2 *(UN-a)
      b2 = +PI_R8/2 *(UN-b)

      alp = -(b2-lg*b1)/(b2-b1)
      bet =      (lg-1)/(b2-b1)

      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      mu = log(sin((PI_R8*xa)/2.))-log(sin((PI_R8*xb)/2.))

      mu = (UN/h)*mu +add_tang(1, deb, fin, alp, bet, mu=0._R8, si=1._R8)
      do i = 1, ia-1
         mu = mu +tang(i*UN, 1, alp, bet, mu=0._R8, si=1._R8)
      enddo
      do i = lg, lg -(ib -2), -1
         mu = mu +tang(i*UN, 1, alp, bet, mu=0._R8, si=1._R8)
      enddo
      mu = mu/lg
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      si = 2.*mu*log(sin((PI_R8*xb)/2.))+cos((PI_R8*xb)/2.)/sin((PI_R8*xb)/2.)-(PI_R8*mu**2.*xb)/2.+(PI_R8*xb)/2.                      &  !
          -2.*mu*log(sin((PI_R8*xa)/2.))+cos((PI_R8*xa)/2.)/sin((PI_R8*xa)/2.)-(PI_R8*mu**2.*xa)/2.+(PI_R8*xa)/2.+PI_R8*mu**2.-PI_R8      !

      si = (UN/h)*si +add_tang(2, deb, fin, alp, bet, mu, si=1._R8)
      do i = 1, ia-1
         si = si+ tang(i*UN, 2, alp, bet, mu, si=1._R8)
      enddo
      do i = lg, lg -(ib -2), -1
         si = si +tang(i*UN, 2, alp, bet, mu, si=1._R8)
      enddo
      si = si/lg
      si = sqrt(si)
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      sk =  log(sin((PI_R8*xb)/2.)**2.)/(2.*si**3.)                                 &  !
            -(3.*mu**2.*log(sin((PI_R8*xb)/2.)))/si**3.                             &  !
            -(3.*mu*cos((PI_R8*xb)/2.))/(si**3.*sin((PI_R8*xb)/2.))                 &  !
            +1/(2.*si**3.*sin((PI_R8*xb)/2.)**2.)+(PI_R8*mu**3.*xb)/(2.*si**3.)     &  !
            -(3.*PI_R8*mu*xb)/(2.*si**3.)                                           &  !
            -log(sin((PI_R8*xa)/2.)**2.)/(2.*si**3.)                                &  !
            +(3.*mu**2.*log(sin((PI_R8*xa)/2.)))/si**3.                             &  !
            -(3.*mu*cos((PI_R8*xa)/2.))/(si**3.*sin((PI_R8*xa)/2.))                 &  !
            -1/(2.*si**3.*sin((PI_R8*xa)/2.)**2.)+(PI_R8*mu**3.*xa)/(2.*si**3.)     &  !
            -(3.*PI_R8*mu*xa)/(2.*si**3.)-(PI_R8*mu**3.)/si**3.+(3.*PI_R8*mu)/si**3.   !

      sk = (UN/h)*sk +add_tang(3, deb, fin, alp, bet, mu, si)
      do i = 1, ia-1
         sk = sk +tang(i*UN, 3, alp, bet, mu, si)
      enddo
      do i = lg, lg -(ib -2), -1
         sk = sk +tang(i*UN, 3, alp, bet, mu, si)
      enddo
      sk = sk/lg
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ku = -(2.*mu*log(sin((PI_R8*xb)/2.)**2.))/si**4.+(4.*mu**3.*log(sin((PI_R8*xb)/2.)))/si**4.                                   &  !
           +(6.*mu**2.*cos((PI_R8*xb)/2.))/(si**4.*sin((PI_R8*xb)/2.))-(4.*cos((PI_R8*xb)/2.))/(3.*si**4.*sin((PI_R8*xb)/2.))       &  !
           -(2.*mu)/(si**4.*sin((PI_R8*xb)/2.)**2.)+cos((PI_R8*xb)/2.)/(3.*si**4.*sin((PI_R8*xb)/2.)**3.)                           &  !
           -(PI_R8*mu**4.*xb)/(2.*si**4.)+(3.*PI_R8*mu**2.*xb)/si**4.-(PI_R8*xb)/(2.*si**4.)                                        &  !
           +(2.*mu*log(sin((PI_R8*xa)/2.)**2.))/si**4.-(4.*mu**3.*log(sin((PI_R8*xa)/2.)))/si**4.                                   &  !
           +(6.*mu**2.*cos((PI_R8*xa)/2.))/(si**4.*sin((PI_R8*xa)/2.))-(4.*cos((PI_R8*xa)/2.))/(3.*si**4.*sin((PI_R8*xa)/2.))       &  !
           +(2.*mu)/(si**4.*sin((PI_R8*xa)/2.)**2.)+cos((PI_R8*xa)/2.)/(3.*si**4.*sin((PI_R8*xa)/2.)**3.)                           &  !
           -(PI_R8*mu**4.*xa)/(2.*si**4.)+(3.*PI_R8*mu**2.*xa)/si**4.-(PI_R8*xa)/(2.*si**4.)                                        &  !
           +(PI_R8*mu**4.)/si**4.-(6.*PI_R8*mu**2.)/si**4.+PI_R8/si**4.                                                                !

      ku = (UN/h)*ku +add_tang(4, deb, fin, alp, bet, mu, si)
      do i = 1, ia-1
         ku = ku +tang(i*UN, 4, alp, bet, mu, si)
      enddo
      do i = lg, lg -(ib -2), -1
         ku = ku +tang(i*UN, 4, alp, bet, mu, si)
      enddo
      ku = ku/lg
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ssk = sk
      sku = ku

   return
   endsubroutine calculs_skku_tan


   real(kind=R8) function add_tang(n, deb, fin, alp, bet, mu, si)
   !================================================================================================
   !<@note Function that adds to the series mean the border integrals as explained in the docs
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real   (kind=R8), intent(in) :: alp    !! *offset so that points are in [b1,b2]*
   real   (kind=R8), intent(in) :: bet    !! *reduction so that points are in [b1,b2]*
   real   (kind=R8), intent(in) :: mu     !! *numerical mean*
   real   (kind=R8), intent(in) :: si     !! *numerical standard deviation*
   integer(kind=I4), intent(in) :: n      !! *statistical moment degree, n=3 for sk and n=4 for ku*
   integer(kind=I4), intent(in) :: fin    !! *last integration point*
   integer(kind=I4), intent(in) :: deb    !! *first integration point*

      real(kind=R8) :: xdeb, xfin

      xdeb = deb
      xfin = fin
      add_tang = (UN/12)*( +9*(tang(xdeb +0.0_R8, n, alp, bet, mu, si)+tang(xfin -0.0_R8, n, alp, bet, mu, si)) &
                           +1*(tang(xdeb +1.0_R8, n, alp, bet, mu, si)+tang(xfin -1.0_R8, n, alp, bet, mu, si)) &
                           -4*(tang(xdeb +0.5_R8, n, alp, bet, mu, si)+tang(xfin -0.5_R8, n, alp, bet, mu, si)) )
   return
   endfunction add_tang


   real(kind=R8) function tang(xi, n, alp, bet, mu, si)
   !================================================================================================
   !<@note Profile function based on the tangent function
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real   (kind=R8), intent(in) :: alp    !! *offset so that points are in [b1,b2]*
   real   (kind=R8), intent(in) :: bet    !! *reduction so that points are in [b1,b2]*
   real   (kind=R8), intent(in) :: mu     !! *numerical mean*
   real   (kind=R8), intent(in) :: si     !! *numerical standard deviation*
   real   (kind=R8), intent(in) :: xi     !! *abscissa*
   integer(kind=I4), intent(in) :: n      !! *statistical moment degree, n=3 for sk and n=4 for ku*

      real(kind=R8) :: tmp

      tmp = (xi +alp)/bet
      tang = ( (tan(tmp) -mu)/si )**n

   return
   endfunction tang


   subroutine calculs_skku_exp3(bounds, lg, ssk, sku)
   !================================================================================================
   !<@note Function to calculate the skewness and kurtosis of an **exponential** series.<br/>
   !< The principle is the same as [[calculs_skku_tan]], however it fits better some particular
   !< series quite binary (roughly two heights).
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real   (kind=R8), intent(in), dimension(1:3) :: bounds  !! *interval limits* [-(1/bounds(1)-1), +(1/bounds(2)-1)]
   integer(kind=I4), intent(in)                 :: lg      !! *vec size*
   real   (kind=R8), intent(out)                :: ssk     !! *theoretical Ssk*
   real   (kind=R8), intent(out)                :: sku     !! *theoretical Sku*

      real(kind=R8) :: xa, xb, kk, kk2, kk3, kk4, mu, si, sk, ku, a, b, c
      real(kind=R8) :: mu2, mu3, mu4, si3, si4
      real(kind=R8) :: h, hh, b1, b2, alp, bet, gam
      real(kind=R8) :: exp1b, exp1a, exp2b, exp2a, exp3b, exp3a, exp4b, exp4a, tmp1a, tmp1b, tmp2a, tmp2b, tmp3a, tmp3b, tmp4a, tmp4b

      real(kind=R8) :: tmp1a2, tmp1b2, tmp1a3, tmp1b3,tmp1a4, tmp1b4, tmp1a5 , tmp1b5 , tmp1a6 , tmp1b6
      real(kind=R8) :: tmp1a7, tmp1b7, tmp1a8, tmp1b8,tmp1a9, tmp1b9, tmp1a10, tmp1b10, tmp1a13, tmp1b13

      integer(kind=I4) :: deb, fin

      deb = 1
      fin = lg

      a = bounds(1)
      b = bounds(2)
      c = bounds(3)

      hh = (-2._R8 + UN/a + UN/b) / (lg - 1)
      h  = hh

      xa = a
      xb = b

      b1 = -(UN - a)/a
      b2 = +(UN - b)/b

      alp = -(b2 - lg*b1)/(b2 - b1)
      bet =      (lg - 1)/(b2 - b1)

      kk  = c / (b2 - b1)**3
      gam = kk

      kk2 = kk**2
      kk3 = kk**3
      kk4 = kk**4

      tmp1a = -1*(UN - UN/xa) ; tmp1a = min(+0.9*HIG_E8, tmp1a) ; tmp1a = max(-0.9*HIG_E8, tmp1a)
      tmp1b = -1*(UN - UN/xb) ; tmp1b = min(+0.9*HIG_E8, tmp1b) ; tmp1b = max(-0.9*HIG_E8, tmp1b)

      tmp2a = -2*(UN - UN/xa) ; tmp2a = min(+0.9*HIG_E8, tmp2a) ; tmp2a = max(-0.9*HIG_E8, tmp2a)
      tmp2b = -2*(UN - UN/xb) ; tmp2b = min(+0.9*HIG_E8, tmp2b) ; tmp2b = max(-0.9*HIG_E8, tmp2b)

      tmp3a = -3*(UN - UN/xa) ; tmp3a = min(+0.9*HIG_E8, tmp3a) ; tmp3a = max(-0.9*HIG_E8, tmp3a)
      tmp3b = -3*(UN - UN/xb) ; tmp3b = min(+0.9*HIG_E8, tmp3b) ; tmp3b = max(-0.9*HIG_E8, tmp3b)

      tmp4a = -4*(UN - UN/xa) ; tmp4a = min(+0.9*HIG_E8, tmp4a) ; tmp4a = max(-0.9*HIG_E8, tmp4a)
      tmp4b = -4*(UN - UN/xb) ; tmp4b = min(+0.9*HIG_E8, tmp4b) ; tmp4b = max(-0.9*HIG_E8, tmp4b)

      tmp1a2  = tmp1a**2  ; tmp1b2  = tmp1b**2
      tmp1a3  = tmp1a**3  ; tmp1b3  = tmp1b**3
      tmp1a4  = tmp1a**4  ; tmp1b4  = tmp1b**4
      tmp1a5  = tmp1a**5  ; tmp1b5  = tmp1b**5
      tmp1a6  = tmp1a**6  ; tmp1b6  = tmp1b**6
      tmp1a7  = tmp1a**7  ; tmp1b7  = tmp1b**7
      tmp1a8  = tmp1a**8  ; tmp1b8  = tmp1b**8
      tmp1a9  = tmp1a**9  ; tmp1b9  = tmp1b**9
      tmp1a10 = tmp1a**10 ; tmp1b10 = tmp1b**10
      tmp1a13 = tmp1a**13 ; tmp1b13 = tmp1b**13

      exp1a = exp(-tmp1a) ! exp( 1 * (1 - 1/xa) )
      exp1b = exp(-tmp1b) ! exp( 1 * (1 - 1/xb) )
      exp2a = exp(-tmp2a) ! exp( 2 * (1 - 1/xa) )
      exp2b = exp(-tmp2b) ! exp( 2 * (1 - 1/xb) )
      exp3a = exp(-tmp3a) ! exp( 3 * (1 - 1/xa) )
      exp3b = exp(-tmp3b) ! exp( 3 * (1 - 1/xb) )
      exp4a = exp(-tmp4a) ! exp( 4 * (1 - 1/xa) )
      exp4b = exp(-tmp4b) ! exp( 4 * (1 - 1/xb) )

      mu = -exp1a+exp1b-(kk*tmp1a4)/4.0D0+(kk*tmp1b4)/4.0D0-1.0D0/xa+1.0D0/xb

      mu = (UN/h)*mu + add_expo3(1, deb, fin, alp, bet, gam, mu=0._R8, si=1._R8)

      mu = mu/lg

      mu2 = mu**2
      mu3 = mu**3
      mu4 = mu**4

      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      si = kk*(-2.4D+1)-exp2a/2.0D0-exp2b/2.0D0+exp1a*(kk*1.2D+1+mu*2.0D0+2.0D0)+exp1b*(kk*1.2D+1-mu*2.0D0+2.0D0)+tmp1a*(mu*2.0D0+mu2+1.0D0)+tmp1b*(mu*(-2.0D0)+mu2+1.0D0)+(kk2*tmp1a7)/7.0D0+(kk2*tmp1b7)/7.0D0+(kk*tmp1a4*(mu+1.0D0))/2.0D0-(kk*tmp1b4*(mu-1.0D0))/2.0D0+kk*exp1a*tmp1a2*6.0D0+kk*exp1a*tmp1a3*2.0D0+kk*exp1b*tmp1b2*6.0D0+kk*exp1b*tmp1b3*2.0D0+kk*exp1a*tmp1a*1.2D+1+kk*exp1b*tmp1b*1.2D+1-3.0D0

      si = (UN/h)*si + add_expo3(2, deb, fin, alp, bet, gam, mu, si=1._R8)

      si = si/lg
      si = sqrt(si)

      si3 = si**3
      si4 = si**4

      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      sk = -1.0D0/si3*(kk*(-2.79D+2/8.0D0)-mu*(9.0D0/2.0D0)+exp3a/3.0D0-kk*mu*3.6D+1-exp2a*(kk*(9.0D0/8.0D0)+mu*(3.0D0/2.0D0)+3.0D0/2.0D0)+tmp1a*(mu*3.0D0+mu2*3.0D0+mu3+1.0D0)+(kk3*tmp1a10)/1.0D+1-kk2*2.16D+3-mu2*3.0D0+exp1a*(kk*3.6D+1+mu*6.0D0+kk*mu*3.6D+1+kk2*2.16D+3+mu2*3.0D0+3.0D0)-kk*exp2a*tmp1a2*(9.0D0/4.0D0)-kk*exp2a*tmp1a3*(3.0D0/2.0D0)+kk*tmp1a4*(mu+1.0D0)**2*(3.0D0/4.0D0)+kk2*tmp1a7*(mu+1.0D0)*(3.0D0/7.0D0)+kk2*exp1a*tmp1a4*9.0D+1+kk2*exp1a*tmp1a5*1.8D+1+kk2*exp1a*tmp1a6*3.0D0-kk*exp2a*tmp1a*(9.0D0/4.0D0)+kk*exp1a*tmp1a*(kk*6.0D+1+mu+1.0D0)*3.6D+1+kk*exp1a*tmp1a2*(kk*6.0D+1+mu+1.0D0)*1.8D+1+kk*exp1a*tmp1a3*(kk*6.0D+1+mu+1.0D0)*6.0D0-1.1D+1/6.0D0)+1.0D0/si3*(kk*(-2.79D+2/8.0D0)+mu*(9.0D0/2.0D0)+exp3b/3.0D0+kk*mu*3.6D+1-exp2b*(kk*(9.0D0/8.0D0)-mu*(3.0D0/2.0D0)+3.0D0/2.0D0)-tmp1b*(mu*3.0D0-mu2*3.0D0+mu3-1.0D0)+(kk3*tmp1b10)/1.0D+1-kk2*2.16D+3-mu2*3.0D0+exp1b*(kk*3.6D+1-mu*6.0D0-kk*mu*3.6D+1+kk2*2.16D+3+mu2*3.0D0+3.0D0)-kk*exp2b*tmp1b2*(9.0D0/4.0D0)-kk*exp2b*tmp1b3*(3.0D0/2.0D0)+kk*tmp1b4*(mu-1.0D0)**2*(3.0D0/4.0D0)-kk2*tmp1b7*(mu-1.0D0)*(3.0D0/7.0D0)+kk2*exp1b*tmp1b4*9.0D+1+kk2*exp1b*tmp1b5*1.8D+1+kk2*exp1b*tmp1b6*3.0D0-kk*exp2b*tmp1b*(9.0D0/4.0D0)+kk*exp1b*tmp1b2*(kk*6.0D+1-mu+1.0D0)*1.8D+1+kk*exp1b*tmp1b3*(kk*6.0D+1-mu+1.0D0)*6.0D0+kk*exp1b*tmp1b*(kk*6.0D+1-mu+1.0D0)*3.6D+1-1.1D+1/6.0D0)

      sk = (UN/h)*sk + add_expo3(3, deb, fin, alp, bet, gam, mu, si)

      sk = sk/lg
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      ku =  1.0D0/si4*(kk*(-6.77962962962963D+1)-mu*(2.2D+1/3.0D0)-exp4a/4.0D0-kk*mu*(2.79D+2/2.0D0)+exp3a*(kk*(8.0D0/2.7D+1)+mu*(4.0D0/3.0D0)+4.0D0/3.0D0)+exp1a*(kk*7.2D+1+mu*1.2D+1+kk*mu*1.44D+2+kk*mu2*7.2D+1+kk2*mu*8.64D+3+kk2*8.64D+3+kk3*1.45152D+6+mu2*1.2D+1+mu3*4.0D0+4.0D0)-kk*mu2*7.2D+1-kk2*mu*8.64D+3+(kk4*tmp1a13)/1.3D+1-kk2*8.60625D+3-kk3*1.45152D+6-mu2*9.0D0-mu3*4.0D0-exp2a*(kk*(9.0D0/2.0D0)+mu*6.0D0+kk*mu*(9.0D0/2.0D0)+kk2*(1.35D+2/4.0D0)+mu2*3.0D0+3.0D0)+tmp1a*(mu*4.0D0+mu2*6.0D0+mu3*4.0D0+mu4+1.0D0)+kk*exp3a*tmp1a2*(4.0D0/3.0D0)+kk*exp3a*tmp1a3*(4.0D0/3.0D0)+kk*tmp1a4*(mu+1.0D0)**3+kk3*tmp1a10*(mu+1.0D0)*(2.0D0/5.0D0)-kk2*exp2a*tmp1a4*(4.5D+1/2.0D0)-kk2*exp2a*tmp1a5*9.0D0-kk2*exp2a*tmp1a6*3.0D0+kk3*exp1a*tmp1a7*2.88D+2+kk3*exp1a*tmp1a8*3.6D+1+kk3*exp1a*tmp1a9*4.0D0+kk2*tmp1a7*(mu+1.0D0)**2*(6.0D0/7.0D0)+kk*exp3a*tmp1a*(8.0D0/9.0D0)+kk*exp1a*tmp1a*(kk*1.2D+2+mu*2.0D0+kk*mu*1.2D+2+kk2*2.016D+4+mu2+1.0D0)*7.2D+1-kk*exp2a*tmp1a2*(kk*1.5D+1+mu*2.0D0+2.0D0)*(9.0D0/2.0D0)-kk*exp2a*tmp1a3*(kk*1.5D+1+mu*2.0D0+2.0D0)*3.0D0+kk2*exp1a*tmp1a4*(kk*1.68D+2+mu+1.0D0)*3.6D+2+kk2*exp1a*tmp1a5*(kk*1.68D+2+mu+1.0D0)*7.2D+1+kk2*exp1a*tmp1a6*(kk*1.68D+2+mu+1.0D0)*1.2D+1+kk*exp1a*tmp1a2*(kk*1.2D+2+mu*2.0D0+kk*mu*1.2D+2+kk2*2.016D+4+mu2+1.0D0)*3.6D+1+kk*exp1a*tmp1a3*(kk*1.2D+2+mu*2.0D0+kk*mu*1.2D+2+kk2*2.016D+4+mu2+1.0D0)*1.2D+1-kk*exp2a*tmp1a*(kk*1.5D+1+mu*2.0D0+2.0D0)*(9.0D0/2.0D0)-2.5D+1/1.2D+1)+1.0D0/si4*(kk*(-6.77962962962963D+1)+mu*(2.2D+1/3.0D0)-exp4b/4.0D0+kk*mu*(2.79D+2/2.0D0)+exp3b*(kk*(8.0D0/2.7D+1)-mu*(4.0D0/3.0D0)+4.0D0/3.0D0)+exp1b*(kk*7.2D+1-mu*1.2D+1-kk*mu*1.44D+2+kk*mu2*7.2D+1-kk2*mu*8.64D+3+kk2*8.64D+3+kk3*1.45152D+6+mu2*1.2D+1-mu3*4.0D0+4.0D0)-kk*mu2*7.2D+1+kk2*mu*8.64D+3+(kk4*tmp1b13)/1.3D+1-kk2*8.60625D+3-kk3*1.45152D+6-mu2*9.0D0+mu3*4.0D0-exp2b*(kk*(9.0D0/2.0D0)-mu*6.0D0-kk*mu*(9.0D0/2.0D0)+kk2*(1.35D+2/4.0D0)+mu2*3.0D0+3.0D0)+tmp1b*(mu*(-4.0D0)+mu2*6.0D0-mu3*4.0D0+mu4+1.0D0)+kk*exp3b*tmp1b2*(4.0D0/3.0D0)+kk*exp3b*tmp1b3*(4.0D0/3.0D0)-kk*tmp1b4*(mu-1.0D0)**3-kk3*tmp1b10*(mu-1.0D0)*(2.0D0/5.0D0)-kk2*exp2b*tmp1b4*(4.5D+1/2.0D0)-kk2*exp2b*tmp1b5*9.0D0-kk2*exp2b*tmp1b6*3.0D0+kk3*exp1b*tmp1b7*2.88D+2+kk3*exp1b*tmp1b8*3.6D+1+kk3*exp1b*tmp1b9*4.0D0+kk2*tmp1b7*(mu-1.0D0)**2*(6.0D0/7.0D0)+kk*exp3b*tmp1b*(8.0D0/9.0D0)+kk*exp1b*tmp1b*(kk*1.2D+2-mu*2.0D0-kk*mu*1.2D+2+kk2*2.016D+4+mu2+1.0D0)*7.2D+1-kk*exp2b*tmp1b2*(kk*1.5D+1-mu*2.0D0+2.0D0)*(9.0D0/2.0D0)-kk*exp2b*tmp1b3*(kk*1.5D+1-mu*2.0D0+2.0D0)*3.0D0+kk*exp1b*tmp1b2*(kk*1.2D+2-mu*2.0D0-kk*mu*1.2D+2+kk2*2.016D+4+mu2+1.0D0)*3.6D+1+kk*exp1b*tmp1b3*(kk*1.2D+2-mu*2.0D0-kk*mu*1.2D+2+kk2*2.016D+4+mu2+1.0D0)*1.2D+1+kk2*exp1b*tmp1b4*(kk*1.68D+2-mu+1.0D0)*3.6D+2+kk2*exp1b*tmp1b5*(kk*1.68D+2-mu+1.0D0)*7.2D+1+kk2*exp1b*tmp1b6*(kk*1.68D+2-mu+1.0D0)*1.2D+1-kk*exp2b*tmp1b*(kk*1.5D+1-mu*2.0D0+2.0D0)*(9.0D0/2.0D0)-2.5D+1/1.2D+1)

      ku = (UN/h)*ku + add_expo3(4, deb, fin, alp, bet, gam, mu, si)

      ku = ku/lg

      ssk = sk
      sku = ku

   contains

      real(kind=R8) function expo3(xi, n, alp, bet, gam, mu, si)
      !================================================================================================
      !<@note Profile function based on the exponential function
      !<
      !<@endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      real   (kind=R8), intent(in) :: alp    !! *offset so that points are in [b1,b2]*
      real   (kind=R8), intent(in) :: bet    !! *reduction so that points are in [b1,b2]*
      real   (kind=R8), intent(in) :: gam    !! **
      real   (kind=R8), intent(in) :: mu     !! *numerical mean*
      real   (kind=R8), intent(in) :: si     !! *numerical standard deviation*
      real   (kind=R8), intent(in) :: xi     !! *abscissa*
      integer(kind=I4), intent(in) :: n      !! *statistical moment degree, n=3 for sk and n=4 for ku*

         real(kind=R8) :: tmp

         tmp  = (xi + alp)/bet

         expo3 = ( (sign(UN, tmp) * (UN - exp(-abs(tmp))) + gam * tmp**3 - mu)/si )**n

      return
      endfunction expo3

      real(kind=R8) function add_expo3(n, deb, fin, alp, bet, gam, mu, si)
      !================================================================================================
      !<@note Function that adds to the series mean the border integrals as explained in the modules
      !< presentation.
      !<
      !<@endnote
      !------------------------------------------------------------------------------------------------
      implicit none
      real   (kind=R8), intent(in) :: alp    !! *offset so that points are in [b1,b2]*
      real   (kind=R8), intent(in) :: bet    !! *reduction so that points are in [b1,b2]*
      real   (kind=R8), intent(in) :: gam    !! **
      real   (kind=R8), intent(in) :: mu     !! *numerical mean*
      real   (kind=R8), intent(in) :: si     !! *numerical standard deviation*
      integer(kind=I4), intent(in) :: n      !! *statistical moment degree, n=3 for sk and n=4 for ku*
      integer(kind=I4), intent(in) :: fin    !! *last integration point*
      integer(kind=I4), intent(in) :: deb    !! *first integration point*

         real(kind=R8) :: xdeb, xfin

         xdeb = deb
         xfin = fin
         add_expo3 = (UN/12)*( +9*(expo3(xdeb +0.0_R8, n, alp, bet, gam, mu, si)+expo3(xfin -0.0_R8, n, alp, bet, gam, mu, si)) &
                               +1*(expo3(xdeb +1.0_R8, n, alp, bet, gam, mu, si)+expo3(xfin -1.0_R8, n, alp, bet, gam, mu, si)) &
                               -4*(expo3(xdeb +0.5_R8, n, alp, bet, gam, mu, si)+expo3(xfin -0.5_R8, n, alp, bet, gam, mu, si)) )
      return
      endfunction add_expo3

   endsubroutine calculs_skku_exp3


   subroutine profil_theo_trie_1D(tab, lg, x, mx)
   !================================================================================================
   !<@note Function that generates the heights when the function limits have been determined.
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer (kind=I4), intent(in )                  :: lg    !! *height vector size*
   real    (kind=R8), intent(out), dimension(1:lg) :: tab   !! *height vector*
   real    (kind=R8), intent(in ), dimension( :  ) :: x     !! *unknowns: height function limits*
   type(MOMENT_STAT), intent(out)                  :: mx    !! *resulting statistical moments*

      real   (kind=R8) :: b1, b2, alp, bet, tmp
      integer(kind=I4) :: i

      select case(PARAM%func_gen)

         case(FCT_TANG)

            b1 = -PI_R8/2 *(UN-x(1))
            b2 = +PI_R8/2 *(UN-x(2))
            alp = -(b2-lg*b1)/(b2-b1)
            bet =     (lg- 1)/(b2-b1)

            do i = 1, lg
               tab(i) = tan( (i*UN+alp)/bet )
            enddo

         case(FCT_EXPO)

            b1 = -(UN-x(1))/x(1)
            b2 = +(UN-x(2))/x(2)
            alp = -(b2-lg*b1)/(b2-b1)
            bet =     (lg- 1)/(b2-b1)
            do i = 1, lg
               tmp = (i*UN+alp)/bet
               tmp = max(-0.9*HIG_E8, tmp)
               tmp = min(+0.9*HIG_E8, tmp)

               tab(i) = sign(UN, tmp) * (UN - exp(-abs(tmp))) + (x(3) / (b2-b1)**3) * tmp**3
            enddo

      endselect

      call std_array(tab = tab(1:lg), mx = mx)

      mx%mu = 0._R8
      mx%si = 1._R8

   return
   endsubroutine profil_theo_trie_1D

endmodule skku_profiles
