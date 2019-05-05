! MOdule containing generic numerical subroutine, like
! non linear root finding, bisection rules and so forth
module gen_numerics
 
 USE data_types, only: double, operator(.speq.) 

 IMPLICIT NONE

 PRIVATE
 PUBLIC zbren, brent, mnbrak


 CONTAINS
!------------------------------------------------------------
!------------------------------------------------------------
!-------------------------------------------------------
!----------------------------------------------------
SUBRoutINE MNBRAK(AX, BX, CX, FA, FB, FC, FUNC, IERROR)
!----------------------------------------------------
! This is part of the "numerical recipes package for line minimization.

! Given a function FUNC, and given the distinct initial points AX and BX
! this routine searches in the downhill direction (defined by the function as
! evaluated at the initial points) and returns the new points AX,BX,CX which
! bracket a minimum of the function. It also returns the function values
! (Fa,FB,FC) at the 3 points

implicit none

real(kind=double), parameter:: GOLD   = 1.618034_double
real(kind=double), parameter:: GLIMIT = 100.0_double
real(kind=double), parameter:: TINY   = 1.e-20_double

integer, intent(out):: IERROR

real(kind=double),intent(INout):: AX, FA 
real(kind=double),intent(INout):: BX, FB 
real(kind=double),intent(out)  :: CX, FC 

real(kind=double) :: r,q,u,ulim,fu

EXTERNAL func
real(kind=double) func

! Check if the initial points make sense
FA=FUNC(AX,IERROR)
if(IERROR.ne.0) return
FB=FUNC(BX,IERROR)
if(IERROR.ne.0) return

! Keep always B as the smallest function value
IF(FB > FA) CALL SWAP_A_B(AX,FA,BX,FB)

CX=BX+GOLD*(BX-AX)
FC=FUNC(CX,IERROR)
if(IERROR.ne.0) return
  


! LOOK FOR THE PROBABLE MINIMUM USING PARABOLAS
1     IF(FB.GE.FC)THEN
  R=(BX-AX)*(FB-FC)
  Q=(BX-CX)*(FB-FA)
  U = BX  - ((BX-CX)*Q-(BX-AX)*R) / (2.0_double*SIGN(MAX(ABS(Q-R),TINY),Q-R))
  ULIM=BX+GLIMIT*(CX-BX)

  IF((BX-U)*(U-CX).GT.0.0_double)THEN

    FU=FUNC(U,IERROR)

    if(IERROR.ne.0) return
    IF(FU.LT.FC)THEN
AX=BX
FA=FB
BX=U
FB=FU
GO TO 1
    ELSE IF(FU.GT.FB)THEN
CX=U
FC=FU
GO TO 1
    ENDIF
    U=CX+GOLD*(CX-BX)
    FU=FUNC(U,IERROR)
    if(IERROR.ne.0) return
  ELSE IF((CX-U)*(U-ULIM).GT.0.0_double)THEN
    FU=FUNC(U,IERROR)
    if(IERROR.ne.0) return
    IF(FU.LT.FC)THEN
BX=CX
CX=U
U=CX+GOLD*(CX-BX)
FB=FC
FC=FU
FU=FUNC(U,IERROR)
    if(IERROR.ne.0) return
    ENDIF
  ELSE IF((U-ULIM)*(ULIM-CX).GE.0.0_double)THEN
    U=ULIM
    FU=FUNC(U,IERROR)
    if(IERROR.ne.0) return
  ELSE
    U=CX+GOLD*(CX-BX)
    FU=FUNC(U,IERROR)
    if(IERROR.ne.0) return
  ENDIF

  AX=BX
  BX=CX
  CX=U
  FA=FB
  FB=FC
  FC=FU
  GO TO 1
ENDIF

CONTAINS

 subroutine swap_a_b(pta,fpta,ptb,fptb)
 implicit none
 ! Swap points a & b in order to have  
  real(kind=double):: pta, fpta
  real(kind=double):: ptb, fptb
  real(kind=double):: dummy
 
  dummy = pta
  pta = ptb
  ptb = dummy

  dummy = fpta
  fpta = fptb
  fptb = dummy

 end subroutine swap_a_b

END SUBRoutINE MNBRAK
!------------------------------------------------------------
!------------------------------------------------------------



!------------------------------------------------------------
!------------------------------------------------------------
real (kind=double) function BRENT(AX,BX,CX,FUNC,TOL,XMIN,IERROR)
!--------------------------------------------------------------
! Given a function FUNC and given a braketing triplet of abscissas AX,BX, and CX
! (BX is between AX and CX and FUNC(BX) is less that both FUNC(AX) and FUNC(CX)), this
! Routine isolates the minimum to a fractional precision about TOL using 
! Brent's method (see "Numerical recipes: the art of scientific computing" for
! details).  The abscissa of the minimum is returned as XMIN and the minimum
! function value is returned as BRENT, the returned function value 

! NOTE: THIS PROGRAM WAS ABUSED BY XP IN SOME UNKNOWN WAYS, AND THAT NEEDS TO
! BE CLEARED, SINCE HIS COMMENTS ARE ABSENT AND HIS VARIABLE NAMES ARE ABSURD
! WARNING: Program does not check if the input data are meaningful, it assumes
! that they were checked by the subroutine above

 implicit none

real(kind=double), intent(IN):: AX,BX,CX
real(kind=double), intent(out):: XMIN
real(kind=double), intent(IN):: TOL
integer, intent(out):: ierror


real(kind=double), parameter ::  CGOLD = .381966_double
integer, parameter :: ITMAX = 100
real(kind=double), parameter:: ZEPS = 1.0e-10_double

real(kind=double):: v,w,x,e,fx,fv,fw,xm,tol1,tol2,r,q,p,etemp,d,u,fu
real(kind=double):: a, b
integer:: iter

external func
real(kind=double) func



a = MIN(AX,CX)
b = MAX(AX,CX)


V = BX
W = V
X = V
E = 0.0_double
FX = FUNC(X,ierror)
if (ierror .ne. 0) return
FV = FX
FW = FX
do ITER = 1,ITMAX
  XM = 0.5_double * (A+B)
  TOL1 = TOL*ABS(X) + ZEPS
  TOL2 = 2.0_double * TOL1
  if (ABS(X-XM) .le. (TOL2-.5_double*(B-A))) goto 3
  if (ABS(E) .gt. TOL1) then
    R = (X-W) * (FX-FV)
    Q = (X-V) * (FX-FW)
    P = (X-V)*Q - (X-W)*R
    Q = 2.0_double * (Q-R)
    if (Q .gt. 0.0_double) P = -P
    Q = ABS(Q)
    ETEMP = E
    E = D
    if (ABS(P).ge.ABS(.5_double*Q*ETEMP) .or. P.le.Q*(A-X) .or.  P.ge.Q*(B-X)) goto 1
    D = P / Q
    U = X + D
    if (U-A.lt.TOL2 .or. B-U.lt.TOL2) D = SIGN(TOL1,XM-X)
    goto 2
   end if
 1 if (X .ge. XM) then
    E = A - X
   else
    E = B - X
   end if
   D = CGOLD * E
 2 if (ABS(D) .ge. TOL1) then
    U = X + D
   else
     U = X + SIGN(TOL1,D)
   end if
   FU = FUNC(U,ierror)
   if (ierror .ne. 0) return
   if (FU <= FX) then
    if (U >= X) then
         A = X
    else
       B = X
    end if
    V = W
    FV = FW
    W = X
    FW = FX
    X = U
    FX = FU
  else
    if (U .lt. X) then
         A = U
    else
         B = U
    end if
    if (FU <= FW .or. (W .speq. X)  ) then
              V = W
              FV = FW
              W = U
              FW = FU
    elseif (FU <= FV .or. (V .speq. X) .or. V <= W) then
              V = U
              FV = FU
    end if
  end if
 enddo

IERROR = 1
 3    XMIN = X
BRENT = FX

end function brent
!------------------------------------------------------
!-------------------------------------------------------
!----------------------------------------------------------------
 subroutine zbren(FUNC_VALUE,FUNC, X1,X2,TOL,zbrent,IERROR)
!----------------------------------------------------------------
! From numerical recipes, uses Van Wijngaarden-Dekker-Brent Method
! to find a root between 2 braketing values (probably).
! The root is supposed to be of the form FUNC(x) = FUNC_VALUE
! Assumes there is only one minimum and that the function is continuous.
! Warning: function does not check for meaningfulness of input parameters
!          I did not check the algorithm of this function (LP UC), but I checked
!          its results and they seemed correct

 implicit none

 real(kind=double), intent(in):: FUNC_VALUE ! The value in the fuction that we
         ! are trying to reproduce
 real(kind=double), intent(in) :: x1,x2
 real(kind=double), intent(in) :: tol


 real(kind=double), intent(out) :: zbrent
 integer, intent(out):: IERROR

 real(kind=double):: tol1,xm,s,p,q,r
 real(kind=double):: a,b,c,d,e,fa,fb,fc
 integer:: iter

 real(kind=double) func
 external func

      
 integer, parameter :: ITMAX = 100 
 real(kind=double),parameter:: EPS=1.e-12_double

 IERROR = 0

 a = x1
 b = x2


 fa = func(a) - FUNC_VALUE
 fb = func(b) - FUNC_VALUE


 if(fb*fa > 0.0_double) IERROR = 1 ! The root wasn't bracketed
       
 fc = fb

 ITMAX_LOOP: DO iter=1,itmax
    IF(fb*fc > 0.0_double) THEN ! check if the root is between b & c
        c  = a
        fc = fa
        d  = b - a
        e  = d
    ENDIF
    IF( ABS(FC) < ABS(FB)) THEN
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
    ENDIF

    tol1 = 2.0_double * eps * ABS(b) + 0.5_double*tol
    xm   = .5_double * (c - b)

    IF( ABS(xm) <= tol1 .OR. ( fb .speq. 0.0_double) )THEN
       zbrent = b
       return
    ENDIF
    IF( ABS(E) >= TOL1 .AND. ABS(FA) > ABS(FB)) THEN
        S=FB/FA
        IF(A .speq. C) THEN
           P = 2.0_double * XM * S
           Q = 1.0_double - S
        ELSE
           Q = FA / FC
           R = FB / FC
           P = S * (  2.0_double * XM * Q * (Q - R) - (B - A)*(R - 1.0_double)   )
           Q = (Q - 1.0_double) * (R - 1.0_double) * (S - 1.0_double)
        ENDIF
        IF(P > 0.0_double) Q = -Q
        P= ABS(P)
        IF(  2.0_double * P < MIN( 3.0_double * XM * Q - ABS(TOL1 * Q), ABS(E * Q) )  ) THEN
           E = D
           D = P / Q
        ELSE
           D = XM
           E = D
        ENDIF
    ELSE
        D = XM
        E = D
    ENDIF
    A  = B
    FA = FB
    IF(ABS(D) > TOL1) THEN
        B = B + D
    ELSE
        B = B + SIGN(TOL1,XM)
    ENDIF
    FB = FUNC(b) - FUNC_VALUE
 ENDDO ITMAX_LOOP

 ZBRENT = B
 IERROR    = 1

!----------------------------------------------------------------
      END SUBRoutINE ZBREN
!----------------------------------------------------------------
!----------------------------------------------------------------


END MODULE gen_numerics
