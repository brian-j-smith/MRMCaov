! DESCRIPTION: subroutines which perform linear algebra operations 
! For the ROC library (performs the functions that are not performed by
! lapack or at least that are still used by some programs and we did
! not move to lapack just yet.

 module l_algebra
 
 use data_types
 implicit none

 private
 
 public:: pseudoinverse ! compute the Moore-Penrose pseudoinverse, for numerically singular hessian matrices
 public::  M2TOM1 !create a vector (upper triangular) of lenght N*(N+1)/2 from a symmetric matrix N*N
 public::  M1TOM2 !create a symmetric matrix N*N from a vector of lenght N*(N+1)/2. i.e. it is 
                 !   the opposite operation of M2TOM1
 public::  SINV   !(Invert a given symmetric positive definite matrix), uses MFSD
 public:: qsortd ! subroutine that applies the quick sorting algorithm, which scales as n*log(n) -- APRIL 2007, LP Uchicago
 public:: median ! subroutine that computes the median -- MAY 2007, LP Uchicago -- taken from Alan Miller's

 contains 

! NOTE1:Pretty much all of these routines are either stolen from Lapack or BLAS
! at some point we should just link Lapack and BLAS and get rid of all 
! this garbage

!---------------------------------------------------------------
SUBROUTINE median(x_in, n, xmed)
!---------------------------------------------------------------
!---------------------------------------------------------------
! PURPOSE: Find the median of X(1), ... , X(N), using as much of the quicksort
!          algorithm as is needed to isolate it.
! NOTE:    If the array X is desired to be partially ordered the array x_in can be called x and the
!          definition of the array x that follows be removed.
! ORIGIN: from Alan Miller's median.f90 code http://users.bigpond.net.au/amiller/


!     Latest revision - 26 November 1996 ( May 2007 by LP - uchicago)
IMPLICIT NONE

INTEGER, INTENT(IN)                :: n
REAL(kind=double), INTENT(IN), DIMENSION(n) :: x_in
REAL(kind=double), INTENT(OUT)                  :: xmed

! Local variables

REAL(kind=double), DIMENSION(n) :: x
REAL(kind=double)    :: temp, xhi, xlo, xmax, xmin
LOGICAL              :: odd
INTEGER              :: hi, lo, nby2, nby2p1, mid, i, j, k

x = x_in ! load the data into the new array, to prevent it from affecting the input

nby2 = n / 2
nby2p1 = nby2 + 1
odd = .true.

!     HI & LO are position limits encompassing the median.

IF (n == 2 * nby2) odd = .false.
lo = 1
hi = n
IF (n < 3) THEN
  IF (n < 1) THEN
    xmed = 0.0
    RETURN
  END IF
  xmed = x(1)
  IF (n == 1) RETURN
  xmed = 0.5_double*(xmed + x(2))
  RETURN
END IF

!     Find median of 1st, middle & last values.

10 mid = (lo + hi)/2
xmed = x(mid)
xlo = x(lo)
xhi = x(hi)
IF (xhi < xlo) THEN          ! Swap xhi & xlo
  temp = xhi
  xhi = xlo
  xlo = temp
END IF
IF (xmed > xhi) THEN
  xmed = xhi
ELSE IF (xmed < xlo) THEN
  xmed = xlo
END IF

! The basic quicksort algorithm to move all values <= the sort key (XMED)
! to the left-hand end, and all higher values to the other end.

i = lo
j = hi
50 DO
  IF (x(i) >= xmed) EXIT
  i = i + 1
END DO
DO
  IF (x(j) <= xmed) EXIT
  j = j - 1
END DO
IF (i < j) THEN
  temp = x(i)
  x(i) = x(j)
  x(j) = temp
  i = i + 1
  j = j - 1

!     Decide which half the median is in.

  IF (i <= j) GO TO 50
END IF

IF (.NOT. odd) THEN
  IF (j == nby2 .AND. i == nby2p1) GO TO 130
  IF (j < nby2) lo = i
  IF (i > nby2p1) hi = j
  IF (i /= j) GO TO 100
  IF (i == nby2) lo = nby2
  IF (j == nby2p1) hi = nby2p1
ELSE
  IF (j < nby2p1) lo = i
  IF (i > nby2p1) hi = j
  IF (i /= j) GO TO 100

! Test whether median has been isolated.

  IF (i == nby2p1) RETURN
END IF
100 IF (lo < hi - 1) GO TO 10

IF (.NOT. odd) THEN
  xmed = 0.5_double*(x(nby2) + x(nby2p1))
  RETURN
END IF
temp = x(lo)
IF (temp > x(hi)) THEN
  x(lo) = x(hi)
  x(hi) = temp
END IF
xmed = x(nby2p1)
RETURN

! Special case, N even, J = N/2 & I = J + 1, so the median is
! between the two halves of the series.   Find max. of the first
! half & min. of the second half, then average.

130 xmax = x(1)
DO k = lo, j
  xmax = MAX(xmax, x(k))
END DO
xmin = x(n)
DO k = i, hi
  xmin = MIN(xmin, x(k))
END DO
xmed = 0.5_double*(xmin + xmax)


!---------------------------------------------------------------
END SUBROUTINE median
!---------------------------------------------------------------

!---------------------------------------------------------------
!---------------------------------------------------------------
pure subroutine qsortd(x, ind, n)
!---------------------------------------------------------------
! Code converted using TO_F90 by Alan Miller
! Date: 2002-12-18  Time: 11:55:47

IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL (kind=double), dimension(n), INTENT(IN)  :: x
INTEGER, dimension(n),  INTENT(OUT)   :: ind


!***************************************************************************

!                                                         ROBERT RENKA
!                                                 OAK RIDGE NATL. LAB.

! THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (kind=double)
! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

!                      X - VECTOR OF LENGTH N TO BE SORTED.

!                    IND - VECTOR OF LENGTH >= N.

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

!*********************************************************************

INTEGER   :: iu(21), il(21)
INTEGER   :: m, i, j, k, l, ij, it, itt, indx
REAL      :: r
REAL (kind=double) :: t

! LOCAL PARAMETERS -

! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X

IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

DO  i = 1, n
  ind(i) = i
END DO
m = 1
i = 1
j = n
r = .375

! TOP OF LOOP

20 IF (i >= j) GO TO 70
IF (r <= .5898437) THEN
  r = r + .0390625
ELSE
  r = r - .21875
END IF

! INITIALIZE K

30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

ij = i + int(r*(j-i)) 
it = ind(ij)
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) > t) THEN
  ind(ij) = indx
  ind(i) = it
  it = indx
  t = x(it)
END IF

! INITIALIZE L

l = j

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T

indx = ind(j)
IF (x(indx) >= t) GO TO 50
ind(ij) = indx
ind(j) = it
it = indx
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) <= t) GO TO 50
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
GO TO 50

! INTERCHANGE ELEMENTS K AND L

40 itt = ind(l)
ind(l) = ind(k)
ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T

50 l = l - 1
indx = ind(l)
IF (x(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60 k = k + 1
indx = ind(k)
IF (x(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED

IF (l-i > j-k) THEN
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  GO TO 80
END IF

il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70 m = m - 1
IF (m == 0) RETURN
i = il(m)
j = iu(m)

80 IF (j-i >= 11) GO TO 30
IF (i == 1) GO TO 20
i = i - 1

! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90 i = i + 1
IF (i == j) GO TO 70
indx = ind(i+1)
t = x(indx)
it = indx
indx = ind(i)
IF (x(indx) <= t) GO TO 90
k = i

100 ind(k+1) = ind(k)
k = k - 1
indx = ind(k)
IF (t < x(indx)) GO TO 100

ind(k+1) = it
GO TO 90

!---------------------------------------------------------------
END SUBROUTINE qsortd
!---------------------------------------------------------------



!---------------------------------------------------------------
!---------------------------------------------------------------
SUBRoutINE M2TOM1(N,XX,VXX)
!---------------------------------------------
!
!THIS SUBRoutINE STORES THE UPPER TRIANGULAR PART OF 
!A SYMMETRIC (N BY N) MATRIX XX COLUMNWISE AS A VECTOR
!VXX WITH LENGTH N*(N+1)/2.
!
implicit none
integer, intent(IN):: N
real(kind=double), DIMENSION(N,N), INTENT(IN):: XX
real(kind=double), DIMENSION(N*(N+1)/2),INTENT(out):: VXX

integer:: K,J,I

K=0
DO J=1,N
  DO I=1,J
    K=K+1
    VXX(K)=XX(I,J)
  enddo
enddo

!---------------------------------------------------------------
END subroutine M2TOM1
!------------------------------------------------------------------
!------------------------------------------------------------------


!---------------------------------------------------------------
!---------------------------------------
SUBRoutINE M1TOM2(N,VXX,XX)
!---------------------------------------
!
!THIS SUBRoutINE CREATES A SYMETRIC (N BY N) MATRIX FROM 
!THE VECTOR VXX WITH LENGTH N*(N+1)/2, WHICH REPRESENTS 
!THE UPPER TRIANGULAR PART OF XX BY COLUMNS.
!
implicit none

integer, intent(IN):: N
real(kind=double), DIMENSION(N,N), INTENT(out):: XX
real(kind=double), DIMENSION(N*(N+1)/2),INTENT(IN):: VXX

integer:: K,J,I

K=0
DO J=1,N
  DO I=1,J
    K=K+1
    XX(I,J)=VXX(K)
    IF(I /= J)XX(J,I)=XX(I,J)
  ENDDO
enddo      
!---------------------------------------------------------------
END SUBRoutINE M1TOM2
!---------------------------------------------------------------
!---------------------------------------------------------------


subroutine pseudoinverse(hessian, num_par, idebug, cov, ierror, err_msg)
! PURPOSE: Implement the Moore-Penrose pseudoinverse when the Hessian matrix
!          (more properly  the Fisher information since it is the expected one) is
!          singular or numerically singular. The matrix becomes singular for small values
!          of d_a. When a small change in d_a is made, d(d_a), it is enough to change
!          all the cutoffs by another quantity f( d(d)a) ), identical for all the cutoffs, 
!          to leave the likelihood functions unchanged. This means that this direction has no
!          effect on the likelihood function or the curve shape, and as such can't affect its
!          variance either, since the curve won't vary at all when these quantities are varied.
!          On the other hand, they will make the Hessian singular. Therefore, this direction is
!          removed in the computation of the variances.
! ALGORITHM: It uses the DSYEV subroutine from LAPACK to diagonalize the negative Hessian. It then eliminates
!            the smallest eigenvalue, after checking that the diagonalization was successful and that all the
!            eigenvalues that were kept were positive
! NOTE:     Currently it just eliminates the smallest eigenvalue. More than one eigenvalue might be eliminated 

use debugging, only: DisplayState
use io, only: line_length

implicit none

! negative of the hessian matrix in the case of MLE estimates
integer, intent(in):: num_par ! number of parameters in estimation
real(kind=double),dimension(num_par,num_par),intent(in):: hessian
integer, intent(in):: idebug ! number of parameters in estimation
! covariance matrix, a buld using the Moore-Penrose pseudoinverse
real(kind=double),dimension(num_par,num_par),intent(out):: cov
integer, intent(out):: ierror ! Error flag
                                    ! 0 => computation successful
                                    ! -1 => wrong input for diagonalization, either someone changed the code, or
                                    !       something went wrong so that the value of some input parameter for the 
                                    !       diagonalization routine assumed an illegal value, checking the value
                                    !       of the offending parameter is likely to lead to the originating problem
                                    ! + 1 => non converged, the matrix could not be diagonalized. Not likely to happen
                                    !        (it never did for us). Write out the hessian matrix and try to figure out
                                    !        why it can't be diagonalized, for example using other algorithms or 
                                    !        programs to diagonalize it, they might be able to tell you more about 
                                    !        the source of the problem. Note that we are using an algorithm for a
                                    !        symmetrical square matrix. If the matrix isn't symmetrical square,
                                    !        someone tampered with the algorithm.
                                    ! +2  => We don't care for the smallest eigenvalue, its value is "noise", but if 
                                    !        any other one is negative (- hessian), then the solution we have is not
                                    !        the maximum of the likelihood function and relevant variances might be
                                    !        negative and likely to be a poor estimates of the variance too, since
                                    !        we aren't at the maximum. 
character(len = line_length), intent(out):: err_msg    ! description of the error occurred, if any

! internal variables
real(kind=double),dimension(num_par,num_par):: upper_diag
! eigenvalues of hessian matrix, used to build the pseudoinverse
real(kind=DOUBLE), dimension(num_par):: eigv 
! transformation matrix
real(kind=double),dimension(num_par,num_par):: T
real(kind=double),dimension(num_par,num_par):: TEMP
 
! diagonalization routine variables
! workspace array, 
real(kind=DOUBLE), dimension(3*num_par-1):: work 
integer:: lwork ! size of workspace array

integer:: i ! loop counter
character(len = line_length):: msg    ! string for log file

lwork = 3*num_par - 1
upper_diag = hessian

! Subroutine from the Lapack library, diagonalizes the matrix
call dsyev('V','U', num_par, upper_diag, num_par,eigv, work, lwork, ierror)

! Prepare to return a warning flag with the error of the diagonalization routine
if(ierror < 0) then ! Diagonalization could not be performed, because the call was made 
    ! incorrectly, this shuld never happen, unless someone changed the code
    write(err_msg,*) "ERROR: pseudoinverse: DSYEV call has wrong value in ",-ierror,"-th parameter"
    ierror = -1
elseif(ierror > 0) then ! diagonalization did not converge
    write(err_msg,*) "ERROR: pseudoinverse: diagonalization not converge, DSYEV err code: ",ierror
    ierror = 1
elseif(  any( eigv(2:num_par) < 0.0_double )  ) then ! The larger iegenvalues have to be positive, otherwise
    ! it isn't a maximum of the log likelihood function and some of the estimated variances could be negative
    ! which don't make any sense.
    write(err_msg,*) "ERROR:: pseudoinverse: One or more large eigenvalues of Hessian are negative"
    ierror = 2
    if(idebug == 1) then ! Write some log information if asked to do so
        call DisplayState(" EIGENVALUES of hessian matrix")
        do i=1,num_par
           write(msg,*) i, eigv(i)
           call DisplayState(msg)
        enddo
    endif
else !ierror == 0 successful diagonalization
    if(idebug == 1) then ! Write some log information if asked to do so
       call DisplayState(" EIGENVALUES of hessian matrix")
       do i=1,num_par
          write(msg,*) i, eigv(i)
          call DisplayState(msg)
       enddo
    endif   
    ! build the matrix that transformes from diagonal form into
    ! our original basis (and transposed goes backward)
    T = upper_diag
    ! Eliminate the smallest eigenvalue. Theoretically one could set a cutoff in size and decide
    ! to eliminated one or more if they are too small. Given that we know that for the proper binorm 
    ! algorithm the singularity is caused by the flat coordinate (d_a, cutoff....), we know that
    ! we need to eliminate only one, and we know which one, the smallest (obviously, since the other
    ! ones are positive, since it is the negative of a the hessian at a maximum, and this one is zero or
    ! too close to it. Sometimes it is negative and very small, but it is never larger than numerical 
    ! noise). Since occasionally (for very flat likelihoods) the inverse seems to be unstable, but it isn't
    ! really (it is more due to a conservative check in the inversion algorithm), here we check whether the
    ! last eigenvalue is too small or not.
    ! LP modified in Sept 2005 to fix some inconsistencies in unbalanced datasets (90 positives and 10 negatives
    ! kind
    TEMP = 0.0_double
    ! Construct the Moore-Penrose inverse diagonal note that we leave the first diagonal element as zero
    ! Note that we don't need abs, because if it is negative, we are sure that we don't want it (when it is 
    ! negative it is always very small (10^-13 or so).
    if( eigv(1) > eigv(2) / 100000.0_double) then
          TEMP(1,1) = 1.0_double/eigv(1)
          if(idebug == 1) then ! Write some log information if asked to do so
                  write(msg,*) "Smallest  Eigv. ", eigv(1), "NOT REMOVED "
                  call DisplayState(msg)
          endif
    else
          if(idebug == 1) then ! Write some log information if asked to do so
                  write(msg,*) "Smallest (Removed) Eigv. ", eigv(1)
                  call DisplayState(msg)
                  write(msg,*) "Next to smallest Eigv. ", eigv(2)
                  call DisplayState(msg)
           endif   
    endif

    do i = 2, num_par
       TEMP(i,i) = 1.0_double/eigv(i)
    enddo
    TEMP = matmul(T,TEMP)
    cov =  matmul(TEMP,  transpose(T))
  endif
  
!---------------------------------------------------

end subroutine pseudoinverse
!---------------------------------------------------
!------------------------------------------------------ 

 
!---------------------------------------------------
!------------------------------------------------------
SUBRoutINE SINV(A,N,EPS,IER)   
!---------------------------------------------------
! INVERT A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX     
! 
!     USAGE   
! CALL SINV(A,N,EPS,IER)  
! 
!     DESCRIPTION OF PARAMETERS 
implicit none

real(kind=double),intent(INout),dimension(*):: A   
!  UPPER TRIANGULAR PART OF THE GIVEN SYMMETRIC POSITIVE     
!   DEFINITE N BY N COEFFICIENT MATRIX.  ON RETURN A    
!   CONTAINS THE RESULTANT UPPER TRIANGULAR MATRIX.     
integer, intent(IN):: N !THE NUMBER OF ROW (COLUMNS) IN MATRIX A   
real(kind=double),intent(IN):: EPS   
!   AN INPUT CONSTANT WHICH IS USED AS RELATIVE TOLERANCE     
!   FOR TEST ON LOSS OF SIGNIFICANCE. 
integer,intent(out):: IER 
!   RESULTING ERROR PARAMETER CODED AS FOLLOWS:   
!   IER=0  - NO ERROR     
!   IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAMETER N OR  
!BECAUSE SOME RADICAND IS NONPOSITIVE (MATRIX A   
!IS NOT POSITIVE DEFINITE, POSSIBLY DUE TO LOSS   
!OF SIGNIFICANCE)   
!   IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFICANCE.    
!THE RADICAND FORMED AT FACTORIZATION STEP K+1    
!WAS STILL POSITIVE BUT NO LONGER GREATER THAN    
!ABS(EPS*A(K+1,K+1)).     
! 
!     REMARKS 
! THE UPPER TRIANGULAR PART OF GIVEN MATRIX IS ASSUMED TO BE  
! STORED COLUMNWISE IN N*(N+1)/2 SUCCESSIVE STORAGE LOCATIONS.  IN  
! THE SAME STORAGE LOCATIONS THE RESULTING UPPER TRIANGULAR MATRIX  
! IS STORED COLUMNWISE TOO.     
! THE PROCEDURE GIVES RESULTS IF N IS GREATER THAN 0 AND ALL  
! CALCULATED RADICANDS ARE POSITIVE.  
! 
!     SUBRoutINES AND FUNCTION SUBPROGRAMS REQUIRED.    
! MFSD  
! 
!     METHOD  
! SOLUTION IS DONE USING THE FACTORIZATION BY SUBRoutINE MFSD.
!

!     FACTORIZE GIVEN MATRIX BY MEANS OF SUBRoutINE MFSD
! 
!     A=TRANSPOSE(T)*T    
! 

 integer:: ipiv, ind, i, kend, lanf, j,k, min, lhor,lver,l
 real(kind=double):: din, work
  
CALL MFSD(A,N,EPS,IER)   

IF(IER < 0) RETURN
 
! 
!     INVERT UPPER TRIANGULAR MATRIX T
!     PREPARE INVERSION-LOOP    
! 
IPIV=N*(N+1)/2 ! This is a division with intergers, handle with care
IND=IPIV  

INVERSION1: DO I=1,N   
  DIN=1.0_double/A(IPIV) 
  A(IPIV)=DIN  
  MIN=N   
  KEND=I-1
  LANF=N-KEND    
  IF( KEND > 0) then  
     J=IND 
     ROW_LOOP1: DO  K=1,KEND  
  WORK=0.0_double 
  MIN=MIN-1    
  LHOR=IPIV
  LVER=J 
  DO  L=LANF,MIN 
     LVER=LVER+1
     LHOR=LHOR+L    
     WORK=WORK+ A(LVER)*A(LHOR) 
  ENDDO
  A(J)=-WORK*DIN 
  J=J-MIN 
     ENDDO ROW_LOOP1
  ENDIF
  IPIV=IPIV-MIN 
  IND=IND-1  
ENDDO INVERSION1

!     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T) 
!     INVERSE(A)=INVERSE(T)*TRANSPOSE(INVERSE(T)) 
!     INITIALIZE MULTIPLICATION LOOP  
! 
INVERSION2:  DO I=1,N  
    IPIV=IPIV+I   
    J=IPIV    
!  INITIALIZE ROW-LOOP
    ROW_LOOP2:DO K=I,N   
WORK=0.0_double
LHOR=J     
DO L=K,N
  LVER=LHOR+K-I 
  WORK=WORK+A(LHOR)*A(LVER) 
  LHOR=LHOR+L
ENDDO
A(J)=WORK
J=J+K
    ENDDO ROW_LOOP2
 ENDDO INVERSION2  

!---------------------------------------------------------------
END  SUBRoutINE SINV   
!---------------------------------------------------------------
!---------------------------------------------------------------


!---------------------------------------------------------------
!---------------------------------------------------------------
SUBRoutINE MFSD(A,N,EPS,IER)  
!---------------------------------------------------------------
!  FACTOR A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX     
! ORIGIN: UKNOWN SO FAR

! NOTE1/WARNING: LP -> The description of this subroutine is not accurate
!  it has been hacked without reporting the hacking 
! 
implicit none
!     DESCRIPTION & DEWFINITION OF DUMMY ARGUMENTS
  integer,Intent(IN):: N !THE NUMBER OF ROW (COLUMNS) IN MATRIX A 
  real(kind=double),INTENT(IN):: EPS ! INPUT CONSTANT USED AS
!    RELATIVE TOLERANCE TEST ON LOSS OF SIGNIFICANCE.
  real(kind=double),intent(INout),dimension(*):: A 
!   AS INPUT AN UPPER TRIANGULAR MATRIX
!   DEFINITE N BY N COEFFICIENT MATRIX.  
!   ON RETURN A CONTAINS THE RESULTANT UPPER TRIANGULAR MATRIX.   
   integer,intent(out)::IER  ! ERROR PARAMETER CODED AS FOLLOWS:  
!   IER=0  - NO ERROR     
!   IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAMETER N OR  
!BECAUSE SOME RADICAND IS NONPOSITIVE (MATRIX A   
!IS NOT POSITIVE DEFINITE, POSSIBLY DUE TO LOSS   
!OF SIGNIFICANCE)   
!   IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFICANCE.    
!THE RADICAND FORMED AT FACTORIZATION STEP K+1    
!WAS STILL POSITIVE BUT NO LONGER GREATER THAN    
!ABS(EPS*A(K+1,K+1)).     
! 
!     REMARKS 
! THE UPPER TRIANGULAR PART OF GIVEN MATRIX IS ASSUMED TO BE  
! STORED COLUMNWISE IN N*(N+1)/2 SUCCESSIVE STORAGE LOCATIONS.  IN  
! THE SAME STORAGE LOCATIONS THE RESULTING UPPER TRIANGULAR MATRIX  
! IS STORED COLUMNWISE TOO.     
! THE PROCEDURE GIVES RESULTS IF N IS GREATER THAN 0 AND ALL  
! CALCULATED RADICANDS ARE POSITIVE.  
! THE PRODUCT OF RETURNED DIAGONAL TERMS IS EQUAL TO THE SQUARE     
! ROOT OF THE DETERMINANT OF THE GIVEN MATRIX.    
! 
!     SUBRoutINES AND FUNCTION SUBPROGRAMS REQUIRED     
! NONE  
! 
!     METHOD  
! SOLUTION IS DONE USING THE SQUARE-ROOT METHOD OF CHOLESKY.  
! THE GIVEN MATRIX IS REPRESENTED AS THE PRODUCT OF 2 TRIANGULAR    
! MATRICES, WHERE THE LEFT HAND FACTOR IS THE TRANSPOSE OF THE
! RETURNED RIGHT HAND FACTOR.   
! 
real (kind=double):: dpiv ! = 1.0e0_double
integer:: kpiv,k,ind,i,l,lanf,lind
real(kind=double):: tol,dsum


! 
!     TEST IF INPUT PARAMETER N IS meaningful 
!     
IF(N-1 < 0) THEN
  IER = -1
  RETURN
ENDIF

!  initialize the error message
IER=0
!     INITIALIZE DIAGONAL-LOOP  
KPIV=0    
DIAGONAL_LOOP: DO K=1,N 
    KPIV=KPIV+K   
    IND=KPIV     
!   CALCULATE TOLERANCE   
    TOL=ABS(EPS*A(KPIV)) 
!   START FACTORIZATION-LOOP OVER K-TH ROW  
    KROW_FACTORIZATION_LOOP: DO I=K,N    
DSUM=0.0_double
IF(k /= 1) then   
   DO  L=1,k-1
LANF = KPIV - L 
LIND = IND - L  
DSUM = DSUM + A(LANF) * A(LIND) 
   ENDDO
ENDIF
!     TRANSFORM ELEMENT A(IND)   
DSUM = A(IND) - DSUM
I_equal_K: IF(I == K) THEN  
!     TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE. 
   IF(  DSUM   <= TOL ) then
IF(DSUM <= 0.0_double ) THEN
   IER=-1; RETURN
ENDIF     
IF (IER <= 0) THEN
  IER = k-1
ENDIF
    ENDIF   
!  COMPUTE PIVOT ELEMENT   
   DPIV=SQRT(DSUM)
   A(KPIV)=DPIV
   DPIV=1.0_double/DPIV
!   CYCLE KROW_FACTORIZATION_LOOP   
ELSE I_equal_k ! Condition of not diagonal
!     CALCULATE TERMS IN ROW
   A(IND) = DSUM * DPIV
ENDIF I_equal_k
IND=IND+I 
  ENDDO KROW_FACTORIZATION_LOOP
ENDDO DIAGONAL_LOOP   
    
!---------------------------------------------------------------
END SUBRoutINE MFSD   
!---------------------------------------------------------------
!---------------------------------------------------------------



end module l_algebra

