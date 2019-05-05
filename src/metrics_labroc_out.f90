! Written by L Pesce - U Chicago, starting fall (autumn) 2009 
! Module that returns specific quantities for a proper ROC curve.
! This quantities are things like:
! 1) Sequence of points on the curve, for plotting purposes
! 2) FPF and TPF for a series of values (normally the cutoffs or category
!    boundaries in the MLE calculation


 module labroc_out

 use data_types, only: double, operator (.speq.), operator (.spne.)
 use computation_constants, only: pi 
 use labroc_functions, only: fpf_CVBM, tpf_CVBM
 use io
 use error_flags 

 implicit none

 private ! default is that objects of this module are not accessible
 public points_on_curve_CVBM ! procedure returns a fpf, tpf sequence for
                               ! plotting purposes. Points are evenly spread

 public points_at_cutoffs_CVBM ! procedure returns a sequence of fpf and tpf for a 
                          ! prespecified (input) sequence of cutoffs (category or run
                          ! boundaries))

 contains

 ! --------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------
  subroutine points_at_cutoffs_CVBM(a_par_in, b_par_in, cutoffs, num_cutoffs, CurvePoints, ierror)
 ! ----------------------------------------------------------------------
 !PURPOSE:  procedure returns a sequence of fpf and tpf for a prespecified (input) sequence of cutoffs (category
 !           boundaries).
 !ALGORITHM: Verify whether the cutoffs and parameters have acceptable values, then simply compute the FPF, TPF values
            
 
 implicit none
 
 real(kind=double), intent(IN):: a_par_in, b_par_in ! curve parameters
 integer, intent(IN) :: num_cutoffs ! number of cutoffs at which (FPF, TPF)  values are  desired
 real(kind=double), dimension(num_cutoffs), intent(in) :: cutoffs ! the actual cutoffs values
 real(kind=double), dimension(2, num_cutoffs + 1), intent(OUT) :: CurvePoints ! the actual
                  !  array with the fpf, tpf values
 integer, intent(OUT):: ierror  ! fit_OK -> OK; fit_failed -> Failed; bad_input -> wrong input
 
 real(kind=double):: fpf, tpf, junk
 real(kind=double):: a_par, b_par ! curve parameters
 real(kind=double):: x

 integer:: i

 ! We are producing points for the purpose of plotting them, so we certainly don't
 ! need an accuracy of more than six digits (which is numerically unattainable anyway
 ! because of the double precision arithmetic and the numerical precision of the 
 ! various functions used)


 if( b_par_in <  0.0_double) then ! unacceptable value
   ierror = bad_input
   return
 else
   b_par = b_par_in
 endif

 ! Take care of a_par. for now we accept any possible value even if negative values are generally looked upon with suspicion.
   a_par = a_par_in

! All cutoffs should be acceptable for this routine.


! Compute the fpf, tpf values for the cutoffs
! Note that the first point is the 1.0 1.0 point
  CurvePoints(1,1) = 1.0_double
  CurvePoints(2,1) = 1.0_double
  do i = 1, num_cutoffs
        x = cutoffs(i)
        call fpf_CVBM(a_par, b_par, x, fpf, junk)
        call tpf_CVBM(a_par, b_par, x, tpf, junk)
        CurvePoints(1,i+1) = fpf
        CurvePoints(2,i+1) = tpf
  enddo



 ! ----------------------------------------------------------------------
  end subroutine points_at_cutoffs_CVBM
 ! ----------------------------------------------------------------------
 ! ----------------------------------------------------------------------



 ! ----------------------------------------------------------------------
 ! ----------------------------------------------------------------------
 subroutine points_on_curve_CVBM(a_par_in, b_par_in, num_pts, CurvePoints,ierror)
 ! ----------------------------------------------------------------------
 !PURPOSE:   return a sequence of fpf,tpf values for the conventional binormal ROC curve with
 !           parameters  a_par, b_par.
 !ALGORITHM: the points will be put at the intersection of a sequence of equally
 !           spaced radii with origin in the point (1,0) of the ROC plot. This choice
 !           has been made to make plotting of skewed curve smooth.  
 !NOTES    : included 0,0 and 1,1.

 implicit none
 
 real(kind=double), intent(IN):: a_par_in, b_par_in ! curve parameters
 integer, intent(IN) :: num_pts ! number of curve points whose value is desired
 real(kind=double), dimension(2, num_pts), intent(OUT) :: CurvePoints ! the actual
                  !  array with the fpf, tpf values
 integer, intent(OUT):: ierror  ! fit_OK -> OK; fit_failed -> Failed; bad_input -> wrong input


 real(kind = double) :: pi2   = pi/2.0_double
 real(kind=double):: target_angle, current_angle ! used to determine the radii 
 real(kind=double):: fpf, tpf, junk, temp
 real(kind=double):: a_par, b_par ! curve parameters
real(kind=double):: x, x_min, x_max ! curve cutoffs
real(kind=double):: x_upper, x_lower ! curve cutoffs

 integer, parameter :: max_iter = 1000
 integer :: i, j_iter ! index for the points to be plotted, and iteration counter

 
 ! Take care of values which produce numerically unstable data points. In this case
 ! we are producing points for the purpose of plotting them, so we certainly don't
 ! need an accuracy of more than six digits (which is numerically unattainable anyway
 ! because of the double precision arithmetic and the numerical precision of the 
 ! various functions used)

 ! CHECK the value of c_par.  c_par == 1 creates numerical problems because the equations lose
 ! meaning.

 if( b_par_in <  0.0_double) then ! unacceptable value we are not going to plot bordeline absurd curves
   ierror = bad_input
   return
 else
   b_par = b_par_in
 endif

 ! Take care of a_par. for now we accept any possible value even if negative values are generally looked upon with suspicion.
   a_par = a_par_in

! All cutoffs should be acceptable for this routine.

! Determine the minimum value of x that it is worthy using.
! Since the curve does not need to be proper and can cross the 45 degree line, and be worse than random everywhere, it 
! follows that we need to check for both FPF and TPF to be about 1.

 call find_x_min(a_par, b_par, x_min, ierror)
 if (ierror == fit_OK) call find_x_max(a_par, b_par, x_max, ierror)
 if(ierror /= fit_OK) return ! failed to bracket solutions, return error message
  
! Verify that the lower boundary for v_c is not too negative, so that the value of
! tpf and fpf doesn't change while the cutoff is halved
 temp = x_min / 2.0_double
 do i = 1, max_iter
          call fpf_CVBM(a_par, b_par, temp, fpf, junk)
          call tpf_CVBM(a_par, b_par, temp, tpf, junk)
          if( tpf > .999999999999999_double .and.  fpf > .999999999999999_double ) then
             x_min = temp
             temp = x_min / 2.0_double
          else
             exit
          endif
 enddo

! If we could not make the value of fpf(min) change, it means that there was a 
! problem
if (i > max_iter) then
    ierror  = fit_failed
    return
endif

! Verify that the maximum is not too positive, for the same reason as done for the minimum
 temp = x_max / 2.0_double
 do i = 1, max_iter
          call fpf_CVBM(a_par, b_par, temp, fpf, junk)
          call tpf_CVBM(a_par, b_par, temp, tpf, junk)
          if( fpf < .000000000000001_double .and.  tpf < .000000000000001_double) then
             x_max = temp
             temp =  x_max / 2.0_double
          else
             exit
          endif
 enddo


! If we could not make the value of fpf(max) change, it means that there was a 
! problem
if (i > max_iter) then
    ierror  = fit_failed
    return
endif


 ! Compute the points on the curve
 
 x = x_max

 target_angle = 0.0_double  ! start sampling angles in ROC space starting from 0,
                            ! which is a point on the TPF = 0 axis


 plotpoints: do i = 1, num_pts
    ! Set value to bracket the cutoff value for this point we are seeking
    x_upper = x
    x_lower = x_min

     ! Start a bisection iteration to find the cutoff value, and thus the FPF/TPF value
     ! that intersects that specific radius from (1,0)
     find_cut_off: do j_iter = 1, max_iter
          x = (x_upper + x_lower ) / 2.0_double
          ! Find fpf and tpf at the new bisected point
          call fpf_CVBM(a_par, b_par, x, fpf, junk)
          call tpf_CVBM(a_par, b_par, x, tpf, junk)

          ! Compute the angle of the current point to see if it is the right
          ! point
          if( (1.0_double - fpf) > 1.0e-10_double) then
                current_angle = ATAN(tpf/(1.0_double - fpf))
          else ! if fpf == 1, then it is the upper right corner
                current_angle = pi2
          endif
          
          ! Determine if the current angle is too large or too small and 
          ! seek another bisection if needed
          if(   (current_angle > target_angle - pi2 / ( (num_pts - 1) * 10) ) &
              .and.  &
                (current_angle < target_angle + pi2 / ( (num_pts - 1) * 10) ) &
           ) then
                target_angle = target_angle + pi2/ (num_pts -1)
                CurvePoints(1,i) = fpf
                CurvePoints(2,i) = tpf
                exit find_cut_off
           elseif( current_angle > target_angle + pi2/ ( (num_pts - 1) * 10) ) THEN
                x_lower = x
                x = ( x + x_upper)/2.0_double
           else
                x_upper = x
                x = ( x + x_lower)/2.0_double
          endif
      enddo find_cut_off



      ! Check whether the algorithm converged for this point
      if(j_iter > max_iter) then
        ierror = fit_failed
        return
      endif

 enddo plotpoints

 ! set the first point to 0,0 and the last point to 1,1
 CurvePoints(1,1) = 0.0_double
 CurvePoints(2,1) = 0.0_double

 CurvePoints(1,num_pts) = 1.0_double
 CurvePoints(2,num_pts) = 1.0_double


 end subroutine points_on_curve_CVBM

!-------------------------------------------------------
  subroutine  find_x_max(a_par, b_par, x_max , ierror)
! PURPOSE: find a cutoffs value that produces a value of FPF and TPF very close
!          to 0. 
! ALGORITHM: Stewise increase (superlinear increase)

  implicit none

  real(kind=double), intent(in):: a_par, b_par ! curve parameters
  real(kind=double), intent(out):: x_max ! the actual minimum
  integer, intent(out):: ierror ! error flag

  real(kind=double):: tpf, one_minus_tpf
  real(kind=double):: step 
  integer, parameter:: num_iter = 10000 ! it never takes too long, but just in case
  integer:: i ! loop counter
 
 ! size of largest change for the FPF and TPF  function as a function of x
 step =  max (1.0_double / b_par , 1.0_double) 

 x_max = 4.0_double ! the maximum is always >= 4, because of the negative cases

 do i =  1, num_iter
       call tpf_CVBM(a_par, b_par, x_max, tpf, one_minus_tpf)
       if( tpf <  1.0e-8_double)  then
            exit
       else
           x_max = x_max + i * step ! make change more than linear
       endif
 enddo 

 if ( i > num_iter) then 
     ierror = fit_failed
 else
     ierror = fit_OK
 endif

!-------------------------------------------------------
 end subroutine find_x_max
!-------------------------------------------------------


!-------------------------------------------------------
 subroutine  find_x_min(a_par, b_par, x_min, ierror)
! PURPOSE: find a cutoffs value that produces a value of FPF and TPF very close
!          to 1. 
! ALGORITHM: Stewise increase (superlinear increase)

  implicit none

  real(kind=double), intent(in):: a_par, b_par ! curve parameters
  real(kind=double), intent(out):: x_min ! the actual minimum
  integer, intent(out):: ierror ! error flag

  real(kind=double):: tpf, one_minus_tpf
  real(kind=double):: step 
  integer, parameter:: num_iter = 10000 ! it never takes too long, but just in case
  integer:: i ! loop counter
 
 ! size of largest change for the FPF and TPF  function as a function of x
 step =  max (1.0_double / b_par , 1.0_double) 

 x_min = - 4.0_double ! the minimum is always <= 0

 do i =  1, num_iter
       call tpf_CVBM(a_par, b_par, x_min, tpf, one_minus_tpf)
       if( tpf > 1.0_double - 1.0e-8_double)  then
            exit
       else
           x_min = x_min - i * step ! make change more than linear
       endif
 enddo 

  if ( i > num_iter) then 
     ierror = fit_failed
  else
     ierror = fit_OK
  endif

!-------------------------------------------------------
 end subroutine find_x_min
!-------------------------------------------------------


 end module labroc_out
