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
use stat_mom,    only : moment_stat, calc_moments, scramble
use sort_arrays, only : sort_array2, init_order
use pikaia_oop,  only : pikaia_class
implicit none

!> {!src/inc_doc/profile_generation.md!}

private

public :: calculs_skku_generique, build_heights


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

      integer(kind=I4) :: istat, fct_sav
      real(kind=R8)    :: cost_val

      type(MOMENT_STAT) :: m_tmp

      real(kind=R8), dimension(1:PARAM%nparam) :: xlower
      real(kind=R8), dimension(1:PARAM%nparam) :: xupper
      real(kind=R8), dimension(1:PARAM%nparam) :: xresul

      ! put input parameters in global variables, so that they can be used in the function "fitness_skku_anal"
      PARAM%m_inp%sk = stats_in%sk
      PARAM%m_inp%ku = stats_in%ku

      ! save PARAM%func_gen value
      fct_sav = PARAM%func_gen

      ! if the Pearson limit is to close to the point (Ssk, Sku), an exponential function is used
      if ( use_fct_expo ) PARAM%func_gen = FCT_EXPO

      ! Genetic algorithm is used to determinate the tangent parameters \alpha and \beta so that, the set of lg heights
      ! will match the statistical moments.
      !..............................................................................

      ! initialization
      xresul(1:PARAM%nparam) = 0.0_R8
      xlower(1:PARAM%nparam) = 0.0_R8
      xupper(1:PARAM%nparam) = 1.0_R8

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
      PARAM%func_gen = fct_sav

      ! height moments calculation
      call calc_moments(   tab = vec_out(1:lg),      &  ! IN
                            mx = m_tmp,              &  ! OUT
                        nb_mom = 4 )                    ! IN

      ! scale and center
      vec_out(1:lg) = ( vec_out(1:lg) - m_tmp%mu ) / m_tmp%si

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
         case (FCT_TANG) ; call calculs_skku_tan(bounds = x, lg = PARAM%npts, ssk = sk, sku = ku)
         case (FCT_EXPO) ; call calculs_skku_exp(bounds = x, lg = PARAM%npts, ssk = sk, sku = ku)
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


   subroutine calculs_skku_exp(bounds, lg, ssk, sku)
   !================================================================================================
   !<@note Function to calculate the skewness and kurtosis of an **exponential** series.<br/>
   !< The principle is the same as [[calculs_skku_tan]], however it fits better some particular
   !< series quite binary (roughly two heights).
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real   (kind=R8), intent(in), dimension(1:2) :: bounds  !! *interval limits* [-(1/bounds(1)-1), +(1/bounds(2)-1)]
   integer(kind=I4), intent(in)                 :: lg      !! *vec size*
   real   (kind=R8), intent(out)                :: ssk     !! *theoretical Ssk*
   real   (kind=R8), intent(out)                :: sku     !! *theoretical Sku*

      real(kind=R8) :: xa, xb, xk, mu, si, sk, ku, a, b, pente
      real(kind=R8) :: h, hh, b1, b2, alp, bet
      real(kind=R8) :: exp1b, exp1a, exp2b, exp2a, exp3b, exp3a, exp4b, exp4a, tmp1a, tmp1b, tmp2a, tmp2b, tmp3a, tmp3b, tmp4a, tmp4b
      integer(kind=I4) :: deb, fin

      pente = UN

      !---------------------------------------------
      ! WXMaxima file
      !---------------------------------------------
      !kill(all);
      !load (f90)$
      !
      !f11(x):=-1+exp(+xk*x)$
      !g11(x):=+1-exp(-xk*x)$
      !
      !f21(x):=f11(x)-mu$
      !g21(x):=g11(x)-mu$
      !
      !f31(x):=f21(x)/si$
      !g31(x):=g21(x)/si$
      !
      !assume(v>0.)$
      !assume(u>0.)$
      !/* [wxMaxima: input   end   ] */
      !
      !/* [wxMaxima: input   start ] */
      !I11:integrate(f11(x),x,-u,0)$
      !J11:integrate(g11(x),x,+0,v)$
      !I11:subst(1./xa-1.,u,I11)$
      !J11:subst(1./xb-1.,v,J11)$
      !I11:I11+J11$
      !I11:expand(trigsimp(I11))$
      !f90(I11);
      !/* [wxMaxima: input   end   ] */
      !
      !/* [wxMaxima: input   start ] */
      !I21:integrate(f21(x)^2,x,-u,0)$
      !J21:integrate(g21(x)^2,x,+0,v)$
      !I21:subst(1./xa-1.,u,I21)$
      !J21:subst(1./xb-1.,v,J21)$
      !I21:I21+J21$
      !I21:expand(trigsimp(I21))$
      !f90(I21);
      !/* [wxMaxima: input   end   ] */
      !
      !/* [wxMaxima: input   start ] */
      !I31:integrate(f31(x)^3,x,-u,0)$
      !J31:integrate(g31(x)^3,x,+0,v)$
      !I31:subst(1./xa-1.,u,I31)$
      !J31:subst(1./xb-1.,v,J31)$
      !I31:I31+J31$
      !I31:expand(trigsimp(I31))$
      !f90(I31);
      !/* [wxMaxima: input   end   ] */
      !
      !/* [wxMaxima: input   start ] */
      !I41:integrate(f31(x)^4,x,-u,0)$
      !J41:integrate(g31(x)^4,x,+0,v)$
      !I41:subst(1./xa-1.,u,I41)$
      !J41:subst(1./xb-1.,v,J41)$
      !I41:I41+J41$
      !I41:expand(trigsimp(I41))$
      !f90(I41);
      !---------------------------------------------

      a = max(bounds(1), 1.e-6_R8)
      b = max(bounds(2), 1.e-6_R8)

      hh = (-2._R8+UN/a+UN/b)/(lg-1)
      h  = hh

      xa = a
      xb = b
      xk = pente

      b1 = -(UN-a)/a
      b2 = +(UN-b)/b

      alp = -(b2-lg*b1)/(b2-b1)
      bet =      (lg-1)/(b2-b1)

      deb = 1
      fin = lg

      tmp1a = 1*(xk-  xk/xa) ; tmp1a = max(-0.9*HIG_E8, tmp1a) ; tmp1a = min(+0.9*HIG_E8, tmp1a)
      tmp1b = 1*(xk-  xk/xb) ; tmp1b = max(-0.9*HIG_E8, tmp1b) ; tmp1b = min(+0.9*HIG_E8, tmp1b)

      tmp2a = 2*(xk-  xk/xa) ; tmp2a = max(-0.9*HIG_E8, tmp2a) ; tmp2a = min(+0.9*HIG_E8, tmp2a)
      tmp2b = 2*(xk-  xk/xb) ; tmp2b = max(-0.9*HIG_E8, tmp2b) ; tmp2b = min(+0.9*HIG_E8, tmp2b)

      tmp3a = 3*(xk-  xk/xa) ; tmp3a = max(-0.9*HIG_E8, tmp3a) ; tmp3a = min(+0.9*HIG_E8, tmp3a)
      tmp3b = 3*(xk-  xk/xb) ; tmp3b = max(-0.9*HIG_E8, tmp3b) ; tmp3b = min(+0.9*HIG_E8, tmp3b)

      tmp4a = 4*(xk-  xk/xa) ; tmp4a = max(-0.9*HIG_E8, tmp4a) ; tmp4a = min(+0.9*HIG_E8, tmp4a)
      tmp4b = 4*(xk-  xk/xb) ; tmp4b = max(-0.9*HIG_E8, tmp4b) ; tmp4b = min(+0.9*HIG_E8, tmp4b)

      exp1a = exp(tmp1a)
      exp1b = exp(tmp1b)
      exp2a = exp(tmp2a)
      exp2b = exp(tmp2b)
      exp3a = exp(tmp3a)
      exp3b = exp(tmp3b)
      exp4a = exp(tmp4a)
      exp4b = exp(tmp4b)

      mu = exp1b/xk-exp1a/xk+1/xb-1/xa
      mu = (UN/h)*mu +add_expo(1, deb, fin, alp, bet, mu=0._R8, si=1._R8)
      mu = mu/lg
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      si = -2*mu*exp1b/xk+2*exp1b/xk-exp2b/xk/2.0+2*mu*exp1a/xk+2*exp1a/xk-exp2a/xk/2.0      &  !
           -3/xk+mu**2/xb-2*mu/xb+1/xb+mu**2/xa+2*mu/xa+1/xa-2*mu**2-2                          !
      si = (UN/h)*si +add_expo(2, deb, fin, alp, bet, mu, si=1._R8)
      si = si/lg
      si = sqrt(si)
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      sk = 3*mu**2*exp1b/(si**3*xk)-6*mu*exp1b/(si**3*xk)+3*exp1b/(si**3*xk)+3.0*mu*exp2b/(2.0*si**3*xk)          &  !
           +(-3.0)*exp2b/(2.0*si**3*xk)+exp3b/(si**3*xk)/3.0-3*mu**2*exp1a/(si**3*xk)-6*mu*exp1a/(si**3*xk)       &  !
           -3*exp1a/(si**3*xk)+3.0*mu*exp2a/(2.0*si**3*xk)+3.0*exp2a/(2.0*si**3*xk)-exp3a/(si**3*xk)/3.0          &  !
           +9*mu/(si**3*xk)-mu**3/(si**3*xb)+3*mu**2/(si**3*xb)-3*mu/(si**3*xb)+1/(si**3*xb)-mu**3/(si**3*xa)     &  !
           -3*mu**2/(si**3*xa)-3*mu/(si**3*xa)-1/(si**3*xa)+2*mu**3/si**3+6*mu/si**3                                 !
      sk = (UN/h)*sk +add_expo(3, deb, fin, alp, bet, mu, si)
      sk = sk/lg
      !. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ku = -4*mu**3*exp1b/(si**4*xk)+12*mu**2*exp1b/(si**4*xk)-12*mu*exp1b/(si**4*xk)+4*exp1b/(si**4*xk)          &  !
           -3*mu**2*exp2b/(si**4*xk)+6*mu*exp2b/(si**4*xk)-3*exp2b/(si**4*xk)+(-4.0)*mu*exp3b/(3.0*si**4*xk)      &  !
           +4.0*exp3b/(3.0*si**4*xk)-exp4b/(si**4*xk)/4.0+4*mu**3*exp1a/(si**4*xk)+12*mu**2*exp1a/(si**4*xk)      &  !
           +12*mu*exp1a/(si**4*xk)+4*exp1a/(si**4*xk)-3*mu**2*exp2a/(si**4*xk)-6*mu*exp2a/(si**4*xk)              &  !
           -3*exp2a/(si**4*xk)+4.0*mu*exp3a/(3.0*si**4*xk)+4.0*exp3a/(3.0*si**4*xk)-exp4a/(si**4*xk)/4.0          &  !
           -18*mu**2/(si**4*xk)+(-25.0)/(6.0*si**4*xk)+mu**4/(si**4*xb)-4*mu**3/(si**4*xb)+6*mu**2/(si**4*xb)     &  !
           -4*mu/(si**4*xb)+1/(si**4*xb)+mu**4/(si**4*xa)+4*mu**3/(si**4*xa)+6*mu**2/(si**4*xa)                   &  !
           +4*mu/(si**4*xa)+1/(si**4*xa)-2*mu**4/si**4-12*mu**2/si**4-2/si**4                                        !
      ku = (UN/h)*ku +add_expo(4, deb, fin, alp, bet, mu, si)
      ku = ku/lg

      ssk = sk
      sku = ku

   return
   endsubroutine calculs_skku_exp


   subroutine calculs_skku_generique(bounds, lg, ssk, sku)
   !================================================================================================
   !<@note Function that calls the right series generator.
   !<
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   real   (kind=R8), intent(in), dimension(:)   :: bounds  !! *interval limits*
   integer(kind=I4), intent(in)                 :: lg      !! *vec size*
   real   (kind=R8), intent(out)                :: ssk     !! *theoretical Ssk*
   real   (kind=R8), intent(out)                :: sku     !! *theoretical Sku*

      select case (PARAM%func_gen)
         case(FCT_TANG) ; call calculs_skku_tan(bounds, lg, ssk, sku)
         case(FCT_EXPO) ; call calculs_skku_exp(bounds, lg, ssk, sku)
      endselect
   return
   endsubroutine calculs_skku_generique


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


   real(kind=R8) function add_expo(n, deb, fin, alp, bet, mu, si)
   !================================================================================================
   !<@note Function that adds to the series mean the border integrals as explained in the modules
   !< presentation.
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
      add_expo = (UN/12)*( +9*(expo(xdeb +0.0_R8, n, alp, bet, mu, si)+expo(xfin -0.0_R8, n, alp, bet, mu, si)) &
                           +1*(expo(xdeb +1.0_R8, n, alp, bet, mu, si)+expo(xfin -1.0_R8, n, alp, bet, mu, si)) &
                           -4*(expo(xdeb +0.5_R8, n, alp, bet, mu, si)+expo(xfin -0.5_R8, n, alp, bet, mu, si)) )
   return
   endfunction add_expo


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


   real(kind=R8) function expo(xi, n, alp, bet, mu, si)
   !================================================================================================
   !<@note Profile function based on the exponential function
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

      real(kind=R8) :: tmp, pente

      pente = UN
      tmp  = (xi +alp)/bet
      tmp = min(+0.9*HIG_E8, tmp)
      tmp = max(-0.9*HIG_E8, tmp)
      expo = ( (sign(UN, tmp)*(UN -exp(-pente*abs(tmp))) -mu)/si )**n

   return
   endfunction expo


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

      real   (kind=R8) :: b1, b2, alp, bet, tmp, pente
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

            pente = UN
            b1 = -(UN-x(1))/x(1)
            b2 = +(UN-x(2))/x(2)
            alp = -(b2-lg*b1)/(b2-b1)
            bet =     (lg- 1)/(b2-b1)
            do i = 1, lg
               tmp = (i*UN+alp)/bet
               tmp = max(-0.9*HIG_E8, tmp)
               tmp = min(+0.9*HIG_E8, tmp)
               tab(i) = sign(UN, tmp)*(UN -exp(-pente*abs(tmp)))
            enddo

      endselect

      call calc_moments(tab = tab(1:lg), mx = mx, nb_mom = 4)

      tab(1:lg) = (tab(1:lg)  -mx%mu) / mx%si ! normalization

      mx%mu = 0._R8
      mx%si = 1._R8

   return
   endsubroutine profil_theo_trie_1D

endmodule skku_profiles
