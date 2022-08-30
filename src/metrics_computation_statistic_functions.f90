!----------------------------------------------------------------------
! File containing wrappers for f77 functions and data for statistic fucntions
! evaluation
!
! It is structured so that first come the cleaned up functions or the functions which
! are wraps and call either some library function or some hack. The hacks follow at the end
! of the file and at some point should either be cleaned up or eliminated.
!
! WARNING:
!         * instability for Phi =1 is removed by setting it to 1-epsilon
!         * f77 routine CUMNOR needs to be standardized to IEEE T float
!         * Need to set references right in the file
!         * dmdbnor is a hack and needs fixing (with its subprocedures)
!         * zdev comes from the original program, it must be changed at some point
!           but we need an elemental procedure, which the ACC does not provide yet
!
!
!
! DATE       PROGRAMMER        DESCRIPTION OF CHANGE
! ====       ==========        =====================
! 11/25/02   L. Pesce (UC)     Creation
! 03/13/03   L. Pesce (UC)     added one_minus_phi & phi1minusphi2 and so forth
!                              Eliminated the
!                              phi - epsilon term in phi, which was used
!                              when phi = 1. I have
!                              eliminated it because the cancellation error
!                              is solved by the new functions
! 03/20/03   L.Pesce (UC)      Added the function no_zero_or_one to prevent any of the
!                              functions from sending back 0 or 1, which would destabilize
!                              the proproc calculations
!                              Moved the INTERFACE statements for the subroutine CUMNOR into
!                              the procedures which call it, otherwise there is a calling mistake
!                              if the run time checks are enables -->> intel bug??
! 10/20/03  L. Pesce (UC)      Added elemental function zdev to compute normal deviates
! 10/22/03  L. Pesce (UC)      Added elemental function bivar_normal_distrib, which calls dmdbnor to
!                              compute bivariate normal density
! 07/08/04  L.Pesce (UC)       Added function norm_dist, and g
! 01/30/06  LP(U)              included SUBROUTINE CUMNOR (FROM ANDERSON CANCER CENTER) AND DEPENDENCIES.
! 12/10/08  LP(U)              Included subroutine to generate a random integer
!  --------------------------------------------------

module statistic_functions
 use data_types, only: double, operator(.speq.) , operator(.spne.)
 use computation_constants
 implicit none

 private
 ! All functions related to the normal distribution are defined as functions
 ! of the normal deviate
 public phi, one_minus_phi, phi1_minus_phi2, bootstrap_set
 public zdev, compute_zdev_plus, bivar_normal_distrib, norm_dist, g

 real(kind=double),parameter:: z_threshold = 38.0_double ! value for which
                              ! in double precision 1-p is equal to 0 with 15 digits
! Modified on June 2008 by LP (UC) to avoid it getting too close to machine infinities. The accuracy is anyway good
! enough because it would be essentially a 1 or  zero when applying phi().
 real(kind=double),parameter:: deviate_infinity = z_threshold !
                   ! value of the normal deviate when the distribution is equal to zero
                   ! or one (exactly would be -/+ infinity)
contains

!------------------------------------------------------------------
subroutine bootstrap_set(n1,n2, nboot, bset)
!------------------------------------------------------------------
! PURPOSE: generate a set, of size nboot, of indices, with replacement,
!          between  n1 and n2.
!          Specifically between min(n1,n2) and max(n1,n2).
!          Usually the number bset will be equal to
!          abs(n1-n2) + 1 and n2 would be the size of the array and
!          n1 will be equal to 1. However, arrays can have any types of
!          indices and the bootstrap can be run with different sizes of
!          the bootstrap set.

integer, intent(IN):: n1 ! lower bound for the random integer
integer, intent(IN):: n2 ! upper bound for the random integer
integer, intent(IN):: nboot ! number of integers to be selected with replacement
integer, dimension(nboot):: bset ! set of the indices for the bootstrap set

integer :: lower_int, upper_int! internally used upper and lower bounds
integer:: i

! Make sure the integers are in the right order
lower_int = min(n1,n2)
upper_int = max(n1,n2)

! set the array, we use a do loop to avoid potential issues with parallelization with
! random number generators, i.e., initiating the SAME sequence a number of times.
do i=1,nboot
    bset(i) = random_integer(n1,n2)
enddo

!------------------------------------------------------------------
end subroutine bootstrap_set
!------------------------------------------------------------------


!------------------------------------------------------------------
integer function random_integer(n1,n2)
!------------------------------------------------------------------
!PURPOSE:   generate a random integer between n1 and n2, extremes included
!           Specifically between min(n1,n2) and max(n1,n2). n1 and n2
!           can be positive or negative and don't need to have the correct ordering.
!ALGORITHM: use the gfortran random number generator random_number, which
!           is based upon a very reliable algorithm (check the manual of
!           gfortran). It assumes that r is [0,1). Int is used because
!           it conserves the sign of the difference, but constrains only the modulus
!           Note that the input is reordered to insure that sign permutations and
!           errors in ordering do not affect the outcome. Given that we first
!           reorder the input, we don't need to worry about the sign of the integers
!           and we can safely use the floor function (as opposed to int())

integer, intent(IN):: n1 ! first bound for the random integer
integer, intent(IN):: n2 ! second bound for the random integer

real(kind=double):: r ! real random number
integer :: lower_int, upper_int! internally used upper and lower bounds

! Make sure the integers are in the right order
lower_int = min(n1,n2)
upper_int = max(n1,n2)

! Generate a random number between [0,1), 1 excluded.
call random_number(r)

! create the integer. Note that we do not need to worry whether the
! numbers are positive or negative because we determine which is the
! smallest value, therefore the interval is above it and the difference
! between the two is by construct positive
random_integer = lower_int + floor( r*(upper_int - lower_int + 1) )

!------------------------------------------------------------------

end function random_integer
!------------------------------------------------------------------




!------------------------------------------------------------------
 subroutine compute_zdev_plus(p, zdev_plus, ierror)
!------------------------------------------------------------------
! Purpose: compute the deviate, but to a preset accuracy. Normally
! deviate functions are inaccurate, and carry a somewhat large error
! Here we use a root finding scheme to reduce the error. This
! under the assumption that the phi function is more accurate
! which is normally the case. Since it is slow, it should be used only
! if needed. It can most likely be made run faster, but for now I don't
! see why should we do that.
use gen_numerics, only: zbren ! to determine the root
implicit none

real(kind=double), intent(in):: p ! cumulative distribution value
real(kind=double), intent(out):: zdev_plus ! the value of the deviate
integer, intent(out):: ierror ! just in case the root-finding fails

real(kind=double) :: deviate2, deviate1 ! bracketing values for  the zero
real(kind=double), parameter:: max_err = 1.0e-10 ! error in the reconstruction of the phi
                      ! normally phi functions can't be computed more accuratly than
                      ! 10^-10

! First compute the deviate in a rough way, using some function
deviate1 = zdev(p)


if( phi(deviate1) < p ) then
! In this case we know that the maximum error is .00045
! so we take more than double that to bracket the exact root
! in case of need we iterate until it is bracketed
! Note that the iteration are done with a larger step to converge more
! quickly in case someting went wrong
      deviate2 = deviate1 + 0.001_double
      do
         if( phi(deviate2) > p ) exit
         deviate2 = deviate2 + 0.005_double
      enddo
elseif( phi(deviate1) > p ) then
     deviate2 = deviate1
     deviate1 = deviate1 - 0.001_double
     do
         if( phi(deviate1) < p ) exit
         deviate1 = deviate1 - 0.005_double
     enddo
else ! It is already the right value
 zdev_plus = deviate1
 ierror = 0
 return
endif

! Find the root with a superlinear method. Newton-Raphson would be a little faster, just in
! in case you need that
 call zbren(p, scalar_phi, deviate1, deviate2, max_err, zdev_plus, ierror)

!------------------------------------------------------------------
end subroutine compute_zdev_plus
!------------------------------------------------------------------

!------------------------------------------------------------------
real(kind=double) function scalar_phi(z)
!------------------------------------------------------------------
! PURPOSE: since the phi function is a pure function, it can't be used
!          as argument for procedure calls on some compilers (It doesn't
!          seem to be part of any standard). This function returns a
!          strictily scalar value, which seems to work.
!------------------------------------------------------------------
  implicit none
  real(kind=double):: z ! argument of the normal distribution

  scalar_phi = phi(z)

!------------------------------------------------------------------
end function scalar_phi
!------------------------------------------------------------------


!------------------------------------------------------------------
 elemental  real(kind=double) function g(x)
!------------------------------------------------------------------
 ! standard normal density
 implicit none
 real(kind=double), intent(IN):: x

 g =  1.0_double/SQRT(2.0_double*pi) * EXP( - x**2 / 2.0_double)


!------------------------------------------------------------------
 end function g
!------------------------------------------------------------------

!------------------------------------------------------------------
 pure subroutine bivar_normal_distrib(x, y, rho, p, ierror)
!------------------------------------------------------------------
! return the value of bivariate binormal distribution for x,y. with
! Error message.
! NOTE: The constant used to constraint the value of rho is chosen in order to render the
!       error in the estimate of  this term ~10^-6 when rho~1. It seems to work, but it should be
!       changed for a better function
! HISTORY: June 4th 2009 LP UC: added 4 9s to round_rho to improve the estimates of partial_aucs
 implicit none
 real(kind=double), intent(IN):: x
 real(kind=double), intent(IN):: y
 real(kind=double), intent(IN):: rho

 real(kind=double), intent(out) :: p
 integer,           intent(out) :: ierror
 real(kind=double):: round_rho ! rho bounded so that it does not create problems to dmdbnor

 round_rho = sign(  min( .9999999999_double, abs(rho) ) , rho    )

 call dmdbnor(x,y,round_rho,p,ierror)

 if (ierror /= 0) ierror = 1 ! set to failed to compute the function

!------------------------------------------------------------------
 end subroutine bivar_normal_distrib
!------------------------------------------------------------------


!----------------------------------------------------------
 elemental real(kind=double) function zdev(p)
!----------------------------------------------------------
!     PURPOSE
!       COMPUTES Z_DEVIATE = P**(-1) (Y), THE ARGUMENT Z_DEVIATE SUCH THAT Y=P(Z_DEVIATE)=
!       THE  PROBABILITY THAT THE RANDOM VARIABLE U, DISTRIBUTED NORMALLY
!       N(0,1), IS LESS THAN OR EQUAL TO Z_DEVIATE.  F(Z_DEVIATE).
!
!     DESCRIPTION OF PARAMETERS
!       P     - INPUT PROBABILITY.
!
!     REMARKS
!       MAXIMUM ERROR IS 0.00045.
!       IF p = 0, ZDEV IS SET TO -deviate_infinity.
!       IF p=1  , ZDEV IS SET TO +deviate_infinity.
!
!     NOTE1:  ORIGINAL PROGRAM SET Z_DEVIATE TO + OR -(10)**74
!     NOTE2:  ORIGINAL PROGRAM WAS SUBRoutINE AND RETURNED ALSO DENSITY AND ERROR CODE
!
!     SUBRoutINES AND SUBPROGRAMS REQUIRED.
!     NONE
!
!     METHOD
!       BASED ON APPROXIMATIONS IN C. HASTINGS, APPROXIMATIONS FOR
!       DIGITAL COMPUTERS, PRINCETON UNIV. PRESS, PRINCETON, N.J., 1955.
!       SEE EQUATION 26.2.23, HANDBOOK OF MATHEMATICAL FUNCTIONS,
!       ABRAMOWITZ AND STEGUN, DOVER PUBLICATIONS, INC., NEW YORK.
!
! WARNING: Need to have check on input and error returns
!     ------------------------------------------------------------------
!
 implicit none
 real(kind=double), intent(IN):: p

 real(kind=double):: wrk_p ! working value for the distribution
 real(kind=double):: t,t2
 real(kind=double):: epsilon_p ! number negligible compared to 1_double


! Check if the probability is a meaningful value [0-1]
! and return specific values if the value is an infinity
! like for p=0 or p=1 (+ infinity and - infinity respectively
 epsilon_p = epsilon(1.0_double)

 IF( p > epsilon_p .and. p < 1.0_double - epsilon_p) THEN ! usual values
       IF( p > 0.5_double) then
         wrk_p = 1.0_double - p
       ELSE
         wrk_p = p
       ENDIF
       t2 = LOG(   1.0_double / ( wrk_p * wrk_p )   )
       t  = SQRT(t2)
       zdev = t - &
               ( 2.515517_double  +  0.802853_double * t  +  0.010328_double * t2 ) &
               /&
               (   1.0_double  +  1.432788_double * t  +  0.189269_double * t2  + &
                   0.001308_double * t * t2 )
       IF(p <= 0.5_double)  zdev = -zdev

 ELSEIF(p > 1.0_double - epsilon_p) THEN ! Plus "infinity"
       zdev = deviate_infinity
 ELSEIF( p < epsilon_p) THEN ! Minus "infinity"
       zdev = - deviate_infinity
 ELSE ! meaningless input value for probability
       zdev = huge(p)
 ENDIF

 END FUNCTION zdev

!------------------------------------------------------------------------

! All the functions related to the normal distribution are considered as functions
! of the normal deviate, no variance or mean value are expected.

 real(kind=double) function one_minus_phi(z)
 ! Computes the difference between 1 and the cumulative distribution function
 ! for the normal distribution for value z
 implicit none
 ! DUMMY VARIABLES
 real(kind=double),intent(IN)::z ! the  argument the normal devite for
                                 ! which z -phi(z) is of interest
 ! work variables
 real(kind=double):: phi  ! used in the call with the canned routine
                                 ! check interface for details


 !  We use the asymptotic approximation from
 ! Abramovitz & Stegun if z is large, otherwise some other
 ! approximation. Note that it is done like this to eliminated
 ! cancellation errors for number too close to 1.
 ! For the goal of the code PROPROC the test for large negative numbers
 ! is not needed, but we included for general numerical stability

   if( abs(z) <= z_threshold) then
      call cumnor(z, phi, one_minus_phi) ! This uses a special procedure from Anderson
                                  ! cancer center, it takes 1 - phi for numerical
                                  ! accuracy
   elseif( z < - z_threshold) then ! z -> -Infinity
       one_minus_phi = 1.0_double - exp_asympthotic_term_for_phi(z )
   else ! z ->  +infinity
        one_minus_phi =  exp_asympthotic_term_for_phi(z)
   endif

 end function one_minus_phi
!----------------------------------------------------------------


!-----------------------------------------------------------------
!-----------------------------------------------------------------
 real(kind=double) function phi1_minus_phi2(z1,z2)
 ! Computes the difference between the cumulative distribution function for the normal
 ! distribution for 2 values of phi
 implicit none
 ! DUMMY VARIABLES
 real(kind=double),intent(IN)::z1,z2 ! the 2 values of the normal devite for
                                    ! which the difference is of interest

 ! INTERNAL VARIABLES


 !  We use the asymptotic approximation from
 ! Abramovitz & Stegun if z is large, otherwise some other
 ! approximation. Note that it is done like this to eliminated
 ! cancellation errors for number too close to 1.
 ! The algorithm is done to avoid:

 ! The only important cancellation is when both terms contain a 1.0
 ! The ifs are structured like this to be more efficient, having the most usual choice first
 ! even if this costs repetition

 if( abs(z1) <= z_threshold .and. abs(z2) <= z_threshold) then
       phi1_minus_phi2 = phi(z1) - phi(z2)
 elseif ( z1 >   z_threshold .and. z2 >  z_threshold)  then ! exponential term opposite sing at + or - INF
       phi1_minus_phi2 =  - exp_asympthotic_term_for_phi(z1) + exp_asympthotic_term_for_phi(z2)
 elseif ( z1 < - z_threshold .and. z2 < - z_threshold   ) then
       phi1_minus_phi2 =   exp_asympthotic_term_for_phi(z1) - exp_asympthotic_term_for_phi(z2)
 elseif ( z1 < - z_threshold .and. z2 >  z_threshold   ) then
       phi1_minus_phi2 =  exp_asympthotic_term_for_phi(z1) &
                          - 1.0_double  +  exp_asympthotic_term_for_phi(z2)
 elseif ( z1 >  z_threshold .and. z2 < - z_threshold   ) then
       phi1_minus_phi2 =  1.0_double - exp_asympthotic_term_for_phi(z1) - exp_asympthotic_term_for_phi(z2)
 else ! if one of the two is large, there is no avoidable cancellation
       phi1_minus_phi2 = phi(z1) - phi(z2)
 endif

   ! Check if the outcome of the calculation is numerically indistinguishable from 0 or 1
    phi1_minus_phi2 = no_zero_or_one(phi1_minus_phi2)


 end function phi1_minus_phi2

!---------------------------------------------------
 subroutine norm_dist(z, phi, one_minus_phi)
!---------------------------------------------------
! PURPOSE: Returns 2 IEEE T float values for the cumulative distribution
!     of a normal random variable  N[0:1] (phi) and its complement to 1 (one_minus_phi).
!     See any elementary statistics book about it.
! ALGORITHM: wrap written by Lorenzo Pesce 7/8/2004 to use canned  routine from
!     Anderson Cancer Center into the program PROPROC or the Ab & Steg asympthotic
!     approximation of phi(z). We use phi and one_minus_phi to minimize cancellation
!     of precision errors.

implicit none

! INPUT: value of the normal deviate z
  real(kind=double), intent (IN) :: z
! outPUT: cumulative distribution function and its complement to 1.
  real(kind=double), intent (out):: phi, one_minus_phi ! phi(z) and  1-phi(z)


  ! When the value of abs(z) becomes too large cumnor produces garbage
  if (z > z_threshold) then
   phi = 1.0_double
   one_minus_phi = 0.0_double
  elseif (z < -z_threshold) then
   phi = 0.0_double
   one_minus_phi = 1.0_double
  else
    call cumnor(z,phi,one_minus_phi) ! This uses a special procedure from Anderson
                                       ! cancer center
  endif

  end subroutine norm_dist
!---------------------------------------------------------

!---------------------------------------------------
 real(kind=double) elemental function phi(z)
!---------------------------------------------------
! Returns a IEEE T float (if the canned routine does so)
!  value for the cumulative distribution  of the normal random variable  N[0:1].
! written by Lorenzo Pesce 11/7/2002 to use canned  routine from
! Anderson Cancer Center into the program PROPROC or the AB & Steg asympthotic
! approximation of phi(z)

implicit none

! INPUT: value of the normal deviate

  real(kind=double), intent (IN) :: z
  real(kind=double):: one_minus_phi ! phi(z) and  1-phi(z)



  if( abs(z) <= z_threshold) then
      call cumnor(z, phi, one_minus_phi) ! This uses a special procedure from Anderson
                                ! cancer center
   elseif( z > z_threshold) then
      phi =  1.0_double - exp_asympthotic_term_for_phi(z)
   else
      phi =   exp_asympthotic_term_for_phi(z)
   endif

   ! Check if the outcome of the calculation is numerically indistinguishable from 0 or 1
!    phi = no_zero_or_one(phi)


  end function phi
!---------------------------------------------------------

 real (kind=double) elemental function Phi_A_Stegun(z)
!---------------------------------------------------
! Returns the cumulative distribution function as a IEEE T float value
! Uses algorithm from Abramovitz & Stegun  (to be checked)

  implicit none
  real(kind=double), intent (IN) :: z
  real(kind=double) :: az ! modulus of z, to remove sign
  real(kind=double) :: T ! work variable
  real(kind=double) :: D ! work variable
  real(kind=double) :: P ! work variable, very related to phi

  Az = ABS(z)

  T = 1.0_double / (1.0_double + 0.2316419_double*az)

  IF ( az < 18.0_double) THEN
         D = 0.3989423_double * EXP(- z*z / 2.0_double)
  ELSE
         D = 0.0_double
  ENDIF

  P = 1.0_double - D * T * &
     (    &
       ( ((1.330274_double * T - 1.821256_double) * T &
          + 1.781478_double) * T - 0.3565638_double &
       ) * T &
       + 0.3193815_double   &
     )

  IF ( z < 0.0_double ) then
      Phi_A_Stegun = 1.0_double - P
  ELSE
     Phi_A_Stegun = p
  ENDIF

  end function Phi_A_Stegun

!------------------------------------------------------------------------
!-----------------------------------------------------------------------
 real(kind=double) elemental function exp_asympthotic_term_for_phi(z)
! Use asymptotic form from Abramovitz and stegun to compute
! One minus phi to the lowest order.
! for z -> -oo   phy_asymptotic = phi(z)
! for z -> +oo   phy_asymptotic = 1 - phi(z)

 implicit none

 real(kind=double),intent(IN):: z
 real(kind=double):: abs_z

 abs_z = abs(z)

 exp_asympthotic_term_for_phi =   exp(-abs_z**2/2.0_double) / (  sqrt(2*pi) * abs_z  )

 end function  exp_asympthotic_term_for_phi

!------------------------------------------------------------------------
!-----------------------------------------------------------------------
  real(kind=double) elemental function no_zero_or_one(phi)
  ! this function checks if the value of its input is close enough to 0 or 1 to be indistinguishable
  ! from them, in which case it changes is of an infinitesimal amount, to prevent numerical errors
  !it mostly has to do with calculations of functions related to the cumulative distribution function

 implicit none

 real(kind=double),  intent(IN)::phi !this is the value than needs to be checked
 intrinsic:: tiny

  if (phi .speq. 1.0_double) then
    no_zero_or_one = phi - 2* TINY(phi)
  elseif (phi .speq. 0.0_double) then
     no_zero_or_one = 2*TINY(phi)
  else
     no_zero_or_one = phi
  endif

  end function no_zero_or_one


!-------------------------------------------------------------------
!-------------------------------------------------------------------
! HACKS, HERE FOLLOW THE FUNCTIONS THAT SHOULD EITHER BE POLISHED OR
! SHOULD BE TAKEN FROM A LIBRARY (BEST THE SECOND OPTION)
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!------------------------------------------------------------------
    elemental  SUBROUTINE dmdbnor (X,Y,RHO,P,IER)
!------------------------------------------------------------------
!
!   FUNCTION            - BIVARIATE NORMAL PROBABILITY DISTRIBUTION
!                           FUNCTION
!   USAGE               - CALL MDBNOR(X,Y,RHO,P,IER)
!   PARAMETERS   X      - INPUT  UPPER LIMIT OF INTEGRATION FOR THE
!                           FIRST VARIABLE
!                Y      - INPUT  UPPER LIMIT OF INTEGRATION FOR THE
!                           SECOND VARIABLE
!                RHO    - INPUT  CORRELATION COEFFICIENT
!                P      - outPUT PROBABILITY THAT THE FIRST VARIABLE
!                           IS LESS THAN OR EQUAL TO X AND THAT THE
!                           SECOND VARIABLE IS LESS THAN OR EQUAL TO Y
!                IER    - ERROR PARAMETER
!                         TERMINAL ERROR = 128+N
!                         N=1 INDICATES THE ABSOLUTE VALUE OF RHO IS
!                             GREATER THAN OR EQUAL TO ONE
!   PRECISION           - DOUBLE
!   REQD. IMSL RoutINES - DMDTNF,DUERTST
!
!                                  C0 = SQRT(.5)
  implicit none

  real(kind=double), intent(IN):: X
  real(kind=double), intent(IN):: Y
  real(kind=double), intent(IN):: RHO
  real(kind=double), intent(out):: P
  integer, intent(out):: ier

  integer:: IAX,IAY,IND,MMM,NNN
  real(kind=double):: EPS,F1,XY,AX,AY,TY,TX,XXX,YYY,QX,QY
  ! Lorenzo Pesce 3/11/04 to eliminate obsolete syntax
  !real(kind=double), parameter :: C0 = 0.707106781186547524_double ! 1/sqrt(2)
  real(kind=double), parameter :: C1 = .99999999999_double ! Added 5 9s, LP - UC June 2005
  real(kind=double), parameter :: C2 = 6.5_double
  real(kind=double), parameter :: C3 = 1.193e-7_double

 ! initialize variables
  NNN = 0
  EPS = 0.0_double
  IER = 0
  P = 0.0_double
! Check if  RHO is  out OF RANGE
  IF ( ABS(RHO) >=  C1) then
        IER = 129
        RETURN
  ENDIF

!  FOR LARGE POSITIVE  X OR Y VALUE - USE GAUSSIAN
!  APPROXIMATION

      IF( X > C2) THEN
             IF( Y > C2) THEN
                P = 1.0_double - abs(epsilon(min(x,y)))
                RETURN
             ELSE
                P = PHI(Y)
                GO TO 8000
             ENDIF
      ELSE
             IF( Y > C2) THEN
                P = PHI(X)
                GO TO 8000
             ELSE
             ENDIF
      ENDIF

!  FOR LARGE NEGATIVE  X OR Y VALUE IT IS BASICALLY 0
       IF(X < -C2 .OR. Y < -C2)THEN
         P = 0.0_double + abs(epsilon(min(x,y)))
         RETURN
       ENDIF

! Neither X nor Y is very large either positive or negative
! so a more precise approximation is needed

      F1 = 1.0_double/SQRT(1.0_double - RHO*RHO)
      XY = X*Y
      IAX = 0
      IAY = 0
      IND = 0
      IF (XY .spne. 0.0_double) THEN
        AX = F1*(Y/X - RHO)
        AY = F1*(X/Y - RHO)
        GO TO 25
      ENDIF

      IF (X .spne. 0.0_double) GO TO 15
      IF (Y .spne. 0.0_double) GO TO 20
!                                                                2 1/2
!                                  FOR X=Y=0 AX=AY=(1-RHO)/(1-RHO )
      AX = F1*(1.0_double - RHO)
      AY = AX
      GO TO 25
!                                  FOR Y=0,X LESS THAN 0     TY = -1/4
!                                  FOR Y=0,X GREATER THAN 0  TY =  1/4
   15 TY = 0.25_double
      IF (X < 0.0_double) TY = -TY
      AX = -F1*RHO
      IND = 1
      GO TO 25
!                                  FOR X=0,Y LESS THAN 0     TX = -1/4
!                                  FOR X=0,Y GREATER THAN 0  TX =  1/4
   20 TX = 0.25_double
      IF (Y < 0.0_double) TX = -TX
      AY = -F1 * RHO
      GO TO 35

   25 IF (AX >= 0.0_double) GO TO 30
      IAX = 1
      AX = -AX
   30 XXX=X
       MMM=0
   33       CALL DMDTNF(XXX,AX,EPS,TX,NNN)
       IF(NNN.GT.100)THEN
         IF(XXX > 0.0_double)THEN
           XXX = XXX - 0.000001_double
         ELSE
           XXX = XXX + 0.000001_double
         ENDIF
         MMM=MMM+1
         IF(MMM > 30)THEN
           IER=1000
           RETURN
         ELSE
           GO TO 33
         ENDIF
       ENDIF
      IF (IAX .NE. 0) TX = -TX
      IF (IND .NE. 0) GO TO 45
   35 IF (AY >= 0.0_double) GO TO 40
      IAY = 1
      AY = -AY
   40 YYY = Y
       MMM = 0
   43       CALL DMDTNF(YYY,AY,EPS,TY,NNN)
       IF(NNN > 100)THEN
         IF(YYY > 0.0_double)THEN
           YYY=YYY - 0.000001_double
         ELSE
           YYY=YYY + 0.000001_double
         ENDIF
         MMM=MMM+1
         IF(MMM > 30)THEN
           IER=1000
           RETURN
         ELSE
           GO TO 43
         ENDIF
       ENDIF
      IF (IAY .NE. 0) TY = -TY
   45 QX = PHI(X)
      QY = PHI(Y)
!                                  NOW EVALUATE P
      P = 0.5_double*(QX + QY) - TX - TY
      IF (XY <= 0.0_double .AND. ( (XY .spne. 0.0_double) .OR. X+Y<0.0_double))  P = P - 0.5_double
      P =  MIN(MAX(0.0_double,P),1.0_double)
 8000 CONTINUE
        IF(P < C3) P = 0.0_double
       IF((1.0_double - P) .LT. C3) P = 1.0_double

  END  SUBRoutINE  dmdbnor
!--------------------------------------------------------------
!--------------------------------------------------------------

!--------------------------------------------------------------------
  ELEMENTAL SUBRoutINE DMDTNF  (Y,Z,EPS,T,nnn)
!--------------------------------------------------------------------
!
!   FUNCTION            - INTEGRATE T(Y,Z) FOR NON-CENTRAL T USAGE.
!   USAGE               - CALL MDTNF(Y,Z,EPS,T)
!   PARAMETERS   Y      - INPUT PARAMETER.  SEE DOCUMENTATION FOR
!                           THE DEFINITION.
!                Z      - INPUT.  INTEGRATION IS FROM 0 TO Z.
!                EPS    - INPUT.  ACCURACY SHOULD NOT BE LESS THAN EPS.
!                           IF EPS=0 IS ENTERED, EPS=.000001 IS USED.
!                T      - RESULTANT VALUE OF THE INTEGRAL.
!   PRECISION           - DOUBLE
!

 implicit none

 real(kind=double), intent(IN):: Y
 real(kind=double), intent(IN):: Z
 real(kind=double), intent(IN):: EPS
 real(kind=double), intent(out):: T
 integer, intent(out):: nnn

  ! Lorenzo Pesce, 03/11/04 eliminated archaic syntax (data structures)
 real(kind=double),parameter:: C = .1591549_double
 real(kind=double),parameter:: EXPOV = 174.673_double

 real(kind=double):: EP1,B,A,TA,HSQB,BEXP,ASQ,A4,B4
 real(kind=double):: A4B4,AHSQB,AB4,F,SUM,G,G1,BER,TER,D1,D2,D,AEPS

      NNN=0
      EP1 = EPS
      IF(EPS .speq. 0.0_double) EP1 = .000001_double
      T = 0.0_double
      B = ABS(Y)
      A = ABS(Z)
      IF(A .speq. 0.0_double) GO TO 35
      TA = ATAN(A)
      IF (A*B <= 4.0_double) GO TO 10
      T = phi(B)
      T = C*(TA + ATAN(1.0_double/A)) - .5_double*(T-.5_double)
      GO TO 30
!                                  APPROXIMATION FOR SMALL Y*Z
   10 HSQB = .5_double*B*B
      IF (HSQB > EXPOV) GO TO 35
      BEXP = EXP(-HSQB)
      ASQ = A*A
      A4 = ASQ*ASQ
      B4 = HSQB * HSQB
      A4B4 = A4 * B4
      AHSQB = A * HSQB
      AB4 = A*B4*.5_double
      F = 1.0_double
      SUM = 0.0_double
      G = 3.0_double
!                                  BEGIN SERIES EXPANSION
   15 G1 = G
      BER = 0.0_double
      TER = AB4
   20 BER = BER+TER
      IF(TER .LE. BER*EP1) GO TO 25
!                                  DEVELOP COEFFICIENT SERIES
      TER = TER*HSQB/G1
      G1 = G1+1.0_double
      GO TO 20
   25 D1 = (BER+AHSQB)/F
      D2 = BER*ASQ/(F+2.0_double)
      D = D1-D2
      SUM = SUM+D
      T = TA-SUM*BEXP
      AEPS = EP1*T
      AHSQB = AHSQB*A4B4/((G-1.0_double)*G)
      AB4 = AB4*A4B4/((G +1.0_double)*G)
      F = F+4.0_double
      G = G+2.0_double
!                                  SHOULD SERIES EXPANSION BE TERMINATED
       IF (D2*BEXP .GE. AEPS)THEN
         NNN=NNN+1
         IF(NNN.GT.100)RETURN
         GO TO 15
       ENDIF
      T = T * C
   30 IF(Z .LT. 0.0_double) T = -T
   35 RETURN
!----------------------------------------------------------
      END SUBRoutINE  DMDTNF
!----------------------------------------------------------
!----------------------------------------------------------


! THIS IS THE FUNCTION CUMNOR AND ITS DEPENDENCIES.
! IT WAS HACKED TO F90 AT THE END OF JANUARY 2006 BY LORENZO PESCE AT THE UNIVERSITY
! OF CHICAGO. NOTE THAT A NUMBER OF "STRANGE FUNCTIONS" THAT APPEAR HERE SHOULD PROBABLY BE
! SCRAPPED IN FORTRAN 90 (LIKE THE IPMPAR) BOTH BECAUSE THEY ARE ARCHAIC AND
! BECAUSE THE MIGHT AFFECT PRECISION.
!
!

   pure subroutine cumnor(arg, result, ccum)
!**********************************************************************
!

!     SUBROUINE CUMNOR(X,RESULT,CCUM)
!
!
!                              Function
!
!
!     Computes the cumulative  of    the  normal   distribution,   i.e.,
!     the integral from -infinity to x of
!          (1/sqrt(2*pi)) exp(-u*u/2) du
!
!     X --> Upper limit of integration.
!                                        X is DOUBLE PRECISION
!
!     RESULT <-- Cumulative normal distribution.
!                                        RESULT is DOUBLE PRECISION
!
!     CCUM <-- Compliment of Cumulative normal distribution.
!                                        CCUM is DOUBLE PRECISION
!
!
!     Renaming of function ANORM from:
!
!     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
!     Package of Special Function Routines and Test Drivers"
!     acm Transactions on Mathematical Software. 19, 22-32.
!
!     with slight modifications to return ccum and to deal with
!     machine constants.
!
!**********************************************************************
!
!
! Original Comments:
!------------------------------------------------------------------
!
! This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!   The main computation evaluates near-minimax approximations
!   derived from those in "Rational Chebyshev approximations for
!   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
!   This transportable program uses rational functions that
!   theoretically approximate the normal distribution function to
!   at least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent
!   constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants.
!
!   MIN   = smallest machine representable number.
!
!   EPS   = argument below which anorm(x) may be represented by
!           0.5  and above which  x*x  will not underflow.
!           A conservative value is the largest machine number X
!           such that   1.0 + X = 1.0   to machine precision.
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ANORM = 0     for  ARG .LE. XLOW.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: March 15, 1992
!
!------------------------------------------------------------------
 use data_types
 implicit none

 real(kind = double), intent(in) ::  arg
 real(kind = double), intent(out) ::  result
 real(kind = double), intent(out) ::  ccum

! interal variables
 integer:: i
 real(kind = double) ::  del, eps, temp, x, xden, xnum, y, xsq, min
!------------------------------------------------------------------
!  Mathematical constants
!
!  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
!  THRSH is the argument for which anorm = 0.75.
 real(kind = double) , parameter:: one    = 1.0_double
 real(kind = double) , parameter:: half   = 0.5_double
 real(kind = double) , parameter:: zero   = 0.0_double
 real(kind = double) , parameter:: sixten = 1.6_double
 real(kind = double) , parameter:: sqrpi  = 3.9894228040143267794e-1_double
 real(kind = double) , parameter:: thrsh  = 0.66291_double
 real(kind = double) , parameter:: root32 = 5.656854248_double
!------------------------------------------------------------------
!  Coefficients for approximation in first interval
!------------------------------------------------------------------
 real(kind = double) , parameter, dimension(5) ::  a = (/ &
         2.2352520354606839287_double,    &
         1.6102823106855587881e2_double,  &
         1.0676894854603709582e3_double,  &
         1.8154981253343561249e4_double,  &
         6.5682337918207449113e-2_double  &
         /)

 real(kind = double) , parameter, dimension(4) ::  b = (/ &
         4.7202581904688241870e1_double,  &
         9.7609855173777669322e2_double,  &
         1.0260932208618978205e4_double,  &
         4.5507789335026729956e4_double   &
         /)
!------------------------------------------------------------------
!  Coefficients for approximation in second interval
!------------------------------------------------------------------
 real(kind = double) , parameter,  dimension(9) ::  c = (/ &
         3.9894151208813466764e-1_double,  &
         8.8831497943883759412e0_double,   &
         9.3506656132177855979e1_double,   &
         5.9727027639480026226e2_double,   &
         2.4945375852903726711e3_double,   &
         6.8481904505362823326e3_double,   &
         1.1602651437647350124e4_double,   &
         9.8427148383839780218e3_double,   &
         1.0765576773720192317e-8_double   &
        /)

 real(kind = double) , parameter,  dimension(8) ::  d = (/ &
         2.2266688044328115691e01_double, &
         2.3538790178262499861e02_double, &
         1.5193775994075548050e03_double, &
         6.4855582982667607550e03_double, &
         1.8615571640885098091e04_double, &
         3.4900952721145977266e04_double, &
         3.8912003286093271411e04_double, &
         1.9685429676859990727e04_double  &
         /)
!------------------------------------------------------------------
!  Coefficients for approximation in third interval
!------------------------------------------------------------------
 real(kind = double) , parameter, dimension(6) ::  p = (/ &
         2.1589853405795699e-1_double,          &
         1.274011611602473639e-1_double,        &
         2.2235277870649807e-2_double,          &
         1.421619193227893466e-3_double,        &
         2.9112874951168792e-5_double,          &
         2.307344176494017303e-2_double         &
         /)
 real(kind = double) , parameter,  dimension(5) ::  q = (/ &
         1.28426009614491121e00_double,         &
         4.68238212480865118e-1_double,         &
         6.59881378689285515e-2_double,         &
         3.78239633202758244e-3_double,         &
         7.29751555083966205e-5_double          &
         /)

!------------------------------------------------------------------

!------------------------------------------------------------------
!  Set Machine dependent constants
!-----------------------------------------------------------------
  eps =  spmpar(1)*0.5D0
  min =  spmpar(2)

!------------------------------------------------------------------
  x = arg
  y = abs(x)
  if  (y <= thrsh) THEN
!------------------------------------------------------------------
!  Evaluate  anorm  for  |X| <= 0.66291
!------------------------------------------------------------------
         xsq = zero
         if (y > eps) xsq = x*x
         xnum = a(5)*xsq
         xden = xsq
         do  i = 1,3
              xnum = (xnum+a(i))*xsq
              xden = (xden+b(i))*xsq
         enddo
         result = x* (xnum+a(4))/ (xden+b(4))
         temp = result
         result = half + temp
         ccum = half - temp
!------------------------------------------------------------------
!  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
!------------------------------------------------------------------
   elseif (y <= root32) THEN
          xnum = c(9)*y
          xden = y
          do  i = 1,7
              xnum = (xnum+c(i))*y
              xden = (xden+d(i))*y
          enddo
          result = (xnum+c(8))/ (xden+d(8))
          xsq = aint(y*sixten)/sixten
          del = (y-xsq)* (y+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          if (x > zero) then
              temp = result
              result = ccum
              ccum = temp
          endif
!------------------------------------------------------------------
!  Evaluate  anorm  for |X| > sqrt(32)
!------------------------------------------------------------------
   else
          result = zero
          xsq = one/ (x*x)
          xnum = p(6)*xsq
          xden = xsq
          do i = 1,4
              xnum = (xnum+p(i))*xsq
              xden = (xden+q(i))*xsq
          enddo
          result = xsq * (xnum+p(5))/ (xden+q(5))
          result = (sqrpi-result)/y
          xsq = aint(x*sixten)/sixten
          del = (x-xsq)* (x+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          if (x > zero) then
              temp = result
              result = ccum
              ccum = temp
          endif

   endif

   if ( result < min) result = 0.0D0
   if ( ccum   < min) ccum = 0.0D0
!------------------------------------------------------------------
  end subroutine cumnor
!------------------------------------------------------------------
!------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 pure real(kind=double) function spmpar(i)
!-----------------------------------------------------------------------
!
!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
!
!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
!
!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
!
!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!
!-----------------------------------------------------------------------
!     WRITTEN BY
!        ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN VIRGINIA
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
!     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
!     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
!-----------------------------------------------------------------------
!     .. Scalar Arguments ..
 use data_types
 implicit none

 integer, intent(in):: i
!     ..
!     .. Local Scalars ..
 real(kind = double) , parameter:: one    = 1.0_double
 real(kind=double):: b, binv, bm1, w, z
 integer:: emax, emin, ibeta, m

!     .. Executable Statements ..

 if ( i <= 1) then
      b = ipmpar(4)
      m = ipmpar(8)
      spmpar = b** (1-m)
 elseif (i <= 2) then
      b = ipmpar(4)
      emin = ipmpar(9)
      binv = one/b
      w = b** (emin+2)
      spmpar = ((w*binv)*binv)*binv
 else
      ibeta = ipmpar(4)
      m = ipmpar(8)
      emax = ipmpar(10)
!
      b = ibeta
      bm1 = ibeta - 1
      z = b** (m-1)
      w = ((z-one)*b+bm1)/ (b*z)
!
      z = b** (emax-2)
      spmpar = ((w*z)*b)*b
 endif

!-----------------------------------------------------------------------
  end function spmpar
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 pure integer function ipmpar(i)
!-----------------------------------------------------------------------
! NOTE: this function is probably useless in fortran 90, but for now I don't
!       care enough to kill it
!     IPMPAR PROVIDES THE INTEGER MA!HINE CONSTANTS FOR THE COMPUTER
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
!
!  INTEGERS.
!
!     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
!
!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
!
!     IPMPAR(1) = A, THE BASE.
!
!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
!
!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS.
!
!     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRE!ISION FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
!
!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
!
!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
!
!     IPMPAR(4) = B, THE BASE.
!
!  SINGLE-PRECISION
!
!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
!
!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
!
!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION
!
!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
!
!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
!
!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
!
!-----------------------------------------------------------------------
!
!     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE
!     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM
!     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN
!     COLUMN 1.)
!
!-----------------------------------------------------------------------
!
!     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
!     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
!     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
!     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
!
!-----------------------------------------------------------------------
!     .. Scalar Arguments ..
 use data_types
 implicit none
 integer, intent(in):: i

!     .. Local Arrays ..
! NOTE: the rubbish for the antique architectures was exterminated to save disk space.
! MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES,
 integer, parameter, dimension(10):: imach = (/ &
          2,          &
          31,         &
          2147483647, &
          2,          &
          24,         &
          -125,       &
          128,        &
          53,         &
          -1021,      &
          1024        &
          /)

!
      ipmpar = imach(i)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      end function ipmpar
!-----------------------------------------------------------------------

!-------------------------------------------------------------------
!-------------------------------------------------------------------
end module statistic_functions

