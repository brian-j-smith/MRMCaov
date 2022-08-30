! This module contains functions which are specific to non parametric ROC analysis
!
!    1) The code is ANSI Fortran 95, it will not compile under
!        strict Fortran 90 or previous versions of the Fortran language
!    2) Deprecated forms of fortran 77 are not used
!    3) The module has the "private" status, so only the methods which are
!       declared as public can be used by external programs


module roc_nonparametric

use data_types ! where data types are defined
use l_algebra, only: qsortd ! quick sort algorithm
use error_flags ! module containing the values associated with specific errors, some are fitting routine specific

! GENERAL SETTING FOR THE WHOLE MODULE
implicit none

private ! All the variables and functions that are not set as public
        ! are not acessible
! Procedures for sensitivity- and specificity-like measures
public  empirical_TPF ! Compute the TPF (sensitivity) associated with a specific threshold value
public  empirical_FPF ! Compute the FPF (1 - specificity) associated with a specific threshold value

public exact_CI_empirical_FPF ! compute exact confidence intervals for the non-parametric estimate of FPF
public exact_CI_empirical_TPF ! compute exact confidence intervals for the non-parametric estimate of TPF

! Global measures, mostly empirical AUC ~ Mann-Whitney form of the Wilcoxon statistics
public  wilcoxon, fast_wilcoxon ! compute wilcoxon from two arrays of data, one for the positive and one for the
                 ! negatives.
public  wilcoxon_cat ! compute wilcoxon from operating points or cumulative fractions
public  ZhouAndGatsonis ! provide Wilcoxon estimates and variance-covariance matrix for a set of
                        ! partially to fully paired (by case) treatment datasets (i.e., the treatments are applied to the same
                        ! cases. It is slightly different from DeLong and DeLong.
public  m_mod_one_shot  ! Dr galls moment based method restricted to straight modality comparison
public  gandpboot        ! bootstrap based method after Dr G and P paper.
public  DeLongAndDeLong ! provide Wilcoxon estimates and variance-covariance matrix for a set of
                        ! fully paired (by case) treatment datasets (i.e., the treatments are applied to the same
                        ! cases

! Curve plotting and calculation of empirical operating points
public  empirical_operating_points_cat ! compute operating points from categorical data
public  empirical_operating_points_list ! compute operating points from list data

contains

!------------------------------------------------------------------
 real(kind=double) pure function wilcoxon_cat(num_categories, op_pts)

! NO SE OUTPUT IS ONLY THE MANN-WHITNEY WILCOXON, used FPF and TPF

! Compute the trapezoidal area using the positive and negative fractions
! This is the same as the wilcoxon statistics if the data are categorical

 implicit none

 integer, intent(IN):: num_categories
 real(kind=double),dimension(2, num_categories+1),intent(IN) ::  op_pts
                  ! Positive and negative fractions, assumes that point 1 is FPF = TPF = 0
                  ! and point n+1 is FPF = TPF = 1 or the opposite.
                  ! array is : ( (FPF_1,TPF_1), (FPF_2,TPF_2), ...)


 integer:: i  ! Loop counter

 wilcoxon_cat = 0.0_double

 ! The larger the category index, the larger the positivity, the smaller the cumulative fraction
 ! because they are selected first. The second term is a sum because we need to add the 1/2 term
 ! for equality to the full term for the smaller (smaller category index means th
 do i = 2, num_categories + 1
       wilcoxon_cat =  wilcoxon_cat + &
                            ( op_pts(1,i) - op_pts(1,i-1) ) * &
                             ( op_pts(2,i) + op_pts(2,i-1) )
 end do

! The division is because one should use the average of the positive fractions not their sum,
! so first we sum them then we divide them by 2
 wilcoxon_cat = wilcoxon_cat/2.0_double

 end function wilcoxon_cat


!----------------------------------------------------------------------
!------------------------------------------------------------------------------
 subroutine empirical_operating_points_list( mn, ms, neg_cases, pos_cases, positiveislarge, &
                                           num_pts, operatingpts,ierror)

!------------------------------------------------------------------------------
! PURPOSE: produce a number of values that can help create a relationship between the test result
!          values and the latent variable space.
! NOTE:    The subroutine prints directly because currently there seem to be no purpose in returning the
!          data to the calling program.
! NOTE1:   The algorithm tests for large AUC values because when the AUC is too large the estimate of beta
!          tends to become unstable.
! WARNING: This subroutine produces an output file and uses it. Handle with care

 use categorization, only: catgrz ! for sorting purposes

 integer, intent(IN):: mn      !number of actually negative cases
 integer, intent(IN):: ms      !number of actually positive cases
 real (kind=double),dimension(mn), intent(IN) :: neg_cases
 real (kind=double),dimension(ms), intent(IN) :: pos_cases
 integer, intent(IN):: positiveislarge ! whether positivity is for more positive or more negative values

 integer, intent(out) :: num_pts ! number of empirical operating points found
 real(kind=double), dimension(2, mn+ms), intent(OUT) :: operatingpts ! the actual
                  !  array with the fpf, tpf values of the empirical operating points
 integer, intent(OUT):: ierror  ! 0 -> OK; 1 -> Failed; 2 -> wrong input

 integer, dimension(2,mn+ms):: cat_data
 integer:: num_cat ! total number of truth runs
 real(kind=double):: delta_fpf, delta_tpf
 integer:: i

! set the coefficient used to rank cases correctly, up if positivity if for more positive values, otherwise down.

if(mn <= 0 .or. ms <= 0) then
    ierror = 2
else
    ierror = 0 ! I can't find a way to make this routine fail when the input is correct
endif

! construct ordinal categories from the real data calling catgrz.
! - no debugging information to external file required ("0" at the 4th variable in the call)
! - no restriction in the number of truth runs (setting 9th variable to mn+ms)

call catgrz(positiveislarge, mn, ms , 0 , cat_data, neg_cases, pos_cases, num_cat, mn+ms )

! define the increase in FPF and TPF obtained by adding one case to the "positive group" for each of them
delta_fpf = 1.0_double/ real(mn)
delta_tpf = 1.0_double/ real(ms)

! Initialize first point
operatingpts(1,1) = 0.0_double
operatingpts(2,1) = 0.0_double

! proceed with the positive cases, downwards
do i = 1, num_cat
    operatingpts(1,i+1) = operatingpts(1,i)   +   cat_data(1,num_cat-i+1) * delta_fpf
    operatingpts(2,i+1) = operatingpts(2,i)   +   cat_data(2,num_cat-i+1) * delta_tpf
enddo

num_pts = num_cat+1 ! it includes the (0,0) and (1,1)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
end subroutine  empirical_operating_points_list
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


!---------------------------------------------------------------
 pure subroutine empirical_operating_points_cat(num_categories, num_cases , categories, op_pts)
! ----------------------------------------------
! computes the cumultive truth fraction for the categories from input
! it assumes that the categories are swept from larger values to smaller
! values as. Note that in this way the ordering of the operating points will be
! the opposite as it is for the categories themselves (i.e., category i will match to
! operating point num_categories + 2 - i -- the first is (0,0)  and the last one is (1,1))

 implicit none

 integer, intent(in):: num_categories
 integer, intent(in):: num_cases
 integer, intent(in), dimension(num_categories) :: categories

 real(kind = double), dimension(num_categories+1), intent(out):: op_pts !(empirical operating
        ! points from (0,0) to (1,1)

 integer:: i ! loop counter

 op_pts  = 0.0_double

 do i = 2, num_categories + 1
      op_pts(i) = op_pts(i-1) + real( categories(num_categories-i+2), kind = double ) / num_cases
 enddo

!------------------------------------------------
 end subroutine empirical_operating_points_cat
!------------------------------------------------
!-------------------------------------------------

!------------------------------------------------
 pure subroutine exact_CI_empirical_FPF(mn, k , confidence_level, type_of_CI,  lb, ub, ierror)
!------------------------------------------------
!PURPOSE:   Compute the variance of the  value of FPF (1 - specificity) for a specific threshold
!           positivity is assumed to be for larger values. We use the number of cases called positive
!           to avoid issues with rounding errors. The calling program has to take care of rounding.
!ALGORITHM: calls empirical_k_CI that uses exact CIs from binomial distribution.

 implicit none

! Declaration of the interface
 integer, intent(IN):: mn ! number of actually-negative cases
 integer, intent(IN):: k  ! number of actually-negative cases cases that were called positive at a specific
                          ! threshold or decision scheme (e.g., a combination of one or more thresholds and a random
                          ! number). We use the integer to avoid issues that could be created by rounding errors
                          ! by forcing the calling program to take care of it.
 real(kind=double), intent(IN):: confidence_level ! e.g., 95% confidence interval a number from 0 to 1.
 integer, intent(IN):: type_of_CI ! -1 -> lower bound is 0 find upper (find the "inferiory" CI) ,
                                  !  0 -> find upper and lower (find the "non-equality" CI), and
                                  ! +1 -> upper bound is 1, find lower (find the superiority CI)

 real(kind=double), intent(OUT)::lb, ub ! estimated lower and upper bound of the CI.
 integer, intent(OUT) :: ierror ! error value for the CI calculation :
                                ! fit_OK -> Procedure did not detect any computation issues
                                ! bad_input -> input values are not acceptable

call empirical_k_CI(mn, k , confidence_level, type_of_CI,  lb, ub, ierror )

!------------------------------------------------
 end subroutine exact_CI_empirical_FPF
!------------------------------------------------
!-------------------------------------------------


!------------------------------------------------
 pure subroutine exact_CI_empirical_TPF(ms, k , confidence_level, type_of_CI,  lb, ub, ierror)
!------------------------------------------------
!PURPOSE:   Compute the variance of the  value of TPF (sensitivity) for a specific threshold
!           positivity is assumed to be for larger values. We use the number of cases called positive
!           to avoid issues with rounding errors.
!ALGORITHM: calls empirical_k_CI that uses exact CIs from binomial distribution.

 implicit none

! Declaration of the interface
 integer, intent(IN):: ms ! number of actually-positive cases
 integer, intent(IN):: k  ! number of actually-positive cases cases that were called positive at a specific
                          ! threshold or decision scheme (e.g., a combination of one or more thresholds and a random
                          ! number). We use the integer to avoid issues that could be created by rounding errors
                          ! by forcing the calling program to take care of it.
 real(kind=double), intent(IN):: confidence_level ! e.g., 95% confidence interval a number from 0 to 1.
 integer, intent(IN):: type_of_CI ! -1 -> lower bound is 0 find upper (find the "inferiory" CI) ,
                                  !  0 -> find upper and lower (find the "non-equality" CI), and
                                  ! +1 -> upper bound is 1, find lower (find the superiority CI)

 real(kind=double), intent(OUT)::lb, ub ! estimated lower and upper bound of the CI.
 integer, intent(OUT) :: ierror ! error value for the CI calculation :
                                ! fit_OK -> Procedure did not detect any computation issues
                                ! bad_input -> input values are not acceptable

call empirical_k_CI(ms, k , confidence_level, type_of_CI,  lb, ub, ierror )

!------------------------------------------------
 end subroutine exact_CI_empirical_TPF
!------------------------------------------------
!-------------------------------------------------


!------------------------------------------------
pure subroutine empirical_k_CI(N, k , confidence_level, type_of_CI,  lb, ub, ierror )
!------------------------------------------------
! PURPOSE: Compute the confidence intervals for a proportion using exact confidence intervals
!          starting from the number of positives observed (k). Here it is meant to be called by
!          specific wrappers for FPF and TPF because it would do exactly the same thing, so we have
!          just one function here. We decided to create to different functions to add some flexibility
!          avoid potential portability issues
!ALGORITHM:from Fleiss "statistical methods for rates and proportions", third edition, Wiley, page 22.
!          the search of largest (smallest) value of p (the probability of observing a success) that
!          has at at least a 5% chance of generating at least as few (at most as many) successes as k
!          is done by first simply bounding the value and applying bisection. More refined approaches can
!          be used, but it did not seem necessary at this point.We use the log of the probability of each
!          observation as basis of the calculations to avoid near constant over- ad under-flows.
!WARNING : The procedure does not have access to the data values therefore it assumes that the k successes
!          were produced by a legitimate algorithm (usually a combination of two cutoffs and a random number
!          generator can help select any point on an empirical ROC curve).

 implicit none

! Declaration of the interface
 integer, intent(IN):: N  ! Total number of actually-negative cases
 integer, intent(IN):: k  ! number of actually-negative cases cases that were called positive at a specific
                          ! threshold or decision scheme (e.g., a combination of one or more thresholds and a random
                          ! number). We use the integer to avoid issues that could be created by rounding errors
                          ! by forcing the calling program to take care of it.
 real(kind=double), intent(IN):: confidence_level ! e.g., 95% confidence interval a number from 0 to 1.
 integer, intent(IN):: type_of_CI ! -1 -> lower bound is 0 find upper (find the "inferiory" CI) ,
                                  !  0 -> find upper and lower (find the "non-equality" CI), and
                                  ! +1 -> upper bound is 1, find lower (find the superiority CI)

 real(kind=double), intent(OUT)::lb, ub ! estimated lower and upper bound of the CI.
 integer, intent(OUT) :: ierror ! error value for the CI calculation :
                                ! fit_OK -> Procedure did not detect any computation issues
                                ! bad_input -> input values are not acceptable

 ! work variables
 integer i ! Loop counter
 real(kind=double) :: tail ! size of the tail of the binomial distribution at fixed number of successes, for which we need to find
                           ! the probability of success -- more about this in the following variables.
 real(kind=double) ::est_p! estimated value of p.
 real(kind=double), dimension(0:N) :: lcoeff    ! array of factorial related coefficients, we compute them first to avoid doing it a
                                                ! million times over. coff(i) = log( N!/(k! (N-k)!)) (binomial coeffs). We use a real
                                                ! array both to avoid multiple conversions and to reduce the chance of overflow.
 real(kind=double) :: lb_kn, lb1, lb1_kn, lb2, lb2_kn ! Bounding values for the lower bound in the success probability with their associated
                                                      ! binomial tails for generating at least k successes (from k to N)
 real(kind=double) :: ub_0k, ub1, ub1_0k, ub2, ub2_0k ! Bounding values for the upper bound in the success probability with their associated
                                                      ! binomial tails for generating at most k successes ( from 0 to k)
 real(kind=double), parameter:: tol = .001_double ! accepted error in the calculation of the tail, doesn't need to be very small. We tried a few
                                                  ! values and this one seems fine
 integer, parameter :: max_iter = 1000 ! max number of iterations

 ! check whether values are acceptable, N needs to be larger than 1 or no variance can be estimated
 if(   k > N .or. k < 0 .or. N <= 1 & ! unacceptable integer values
       .or. confidence_level > 1.0 .or. confidence_level < 0.0 & ! unacceptale confidence values
       .or.(  (type_of_CI.ne.1) .and. (type_of_CI.ne.0) .and. (type_of_CI.ne.1)) & ! wrong input value for the type of CI
    ) then
         ! Initialize the return values for some absurd values -- this 666 comes from Kevin Schatz, UoI :-)
         ! values are in wrong order to make sure user has to realize there is a problem
         lb = +666.0_double
         ub = -666.0_double
         ierror = bad_input
         return
 elseif( confidence_level .speq. 1.0_double ) then
        !Only way to be sure is to cover every possible answer. Test done for completeness. It is absurd but acceptable
         lb = 0.0_double
         ub = 1.0_double
         ierror = 0
         return
 elseif( confidence_level .speq. 0.0_double ) then
         !This one is a little obscure to me, but this seems that as the confidence level descreases
         !it should collapse around the estimate. Test done for completeness. it is absurd, but acceptable.
         lb = real (k/N, kind=double)
         ub = lb
         ierror = 0
         return
 else
         ! It looks like we have a "reasonable" input.
         ierror = 0
 endif

 ! compute the log binomial coefficients, they are rarely all needed in the calcs, but I can't find an elegant to avoid computing all of them
 ! anyway and since the  memory is more or less irrelevant, it cleaner to just keep them all. Double precision overflows for N > ~ 1000, so
 ! we used logarithms.
 lcoeff(0) = 0.0_double
 do i = 1, N
       lcoeff(i) = lcoeff(i-1) + log(real(N-i+1, kind=double))  -log(real(i,kind=double))
 enddo


! first estimate the expected value of p, the proportion.
 est_p = real(k, kind=double) /real(N, kind=double)



! We need to know if there is one tail or two tails for that specific confidence level.
 if(type_of_CI == 0) then
        tail = (1.0 - confidence_level)/2.0_double
 else
        tail = (1.0 - confidence_level)
 endif

! We have already corrected for the tail size, now we have to see if the lower bound tail is
! present:
! - if we look only for an upper bound (type_of_CI = -1), it means we are not interested in a
!   lower bound, i.e., the bound can be set to zero (the lowest possible value)
! - if we observed zero successes, then the lower bound can only be zero (because the cumulative
!   (probability of measures at least as large as 0 becomes 1 and there are no smaller
!    values than 0).
 search_lb: if(type_of_CI == -1 .or. k == 0) then
        lb = 0.0_double
 else
        lb1 = 0.0_double ! lowest possible value. We could be more precise, but that would save us at most
                         ! a couple of iterations.
        lb1_kn = 0.0_double ! if  the probability of success is 0, it impossible to generate
                            ! k > 0 successes
        lb2 = est_p ! the lower bound has to be smaller than the estimated value.
        lb2_kn = bin_dist_kn(lb2)

        ! iterate between the two values to find a value that it is close enough
        do i = 1, max_iter
             ! bisect and compute the tail probability
             lb = (lb1 + lb2)/2.0_double
             lb_kn = bin_dist_kn(lb)
             ! check whether the tail probability is close to the desired value or the bounds are converging to
             ! each other.
             if( abs(lb_kn - tail)/tail < tol .or. abs(lb1 - lb2 ) < min(.1_double/N, tol*(lb1+lb2)/2) ) then
                  exit
             elseif( lb_kn < tail ) then ! lb is too small, so we set it as new lower limit
                  lb1 = lb
                  lb1_kn = lb_kn
             else ! lb is too large, so we set it as upper bound
                  lb2 = lb
                  lb2_kn = lb_kn
             endif
        enddo
 endif search_lb

! We have already corrected for the tail size, now we have to see if the upperer bound tail is
! present:
! - if we look only for an lower bound (type_of_CI = +1), it means we are not interested in a
!   upper bound, i.e., the bound can be set to 1 (the largest possible value)
! - if we observed N successes, then the upper bound can only be 1 (because the cumulative
!   (probability of measures at most as extreme as N becomes 1 and there are no larger
!    values than 1).

 search_ub: if(type_of_CI == +1 .or. k == N) then
        ub = 1.0_double
 else
        ub1 = est_p ! the upper bound has to be larger than the estimated value
        ub1_0k = bin_dist_0k(ub1)
        ub2 = 1.0_double ! set the current bound to its highest possible value, 1.
        ub2_0k = 0.0_double ! since k is different from N, a probability of less than N successes is zero for ub = 1

        ! iterate between the two values to find a value that it is close enough
        do i = 1, max_iter
             ! bisect and compute the tail probability
             ub = (ub1 + ub2)/2.0_double
             ub_0k = bin_dist_0k(ub)
             ! check whether the tail probability is close to the desired value or the bounds are converging to
             ! each other.
             if( abs(ub_0k - tail)/tail < tol .or. abs(ub1 - ub2 ) < min(.1_double/N, tol*(ub1+ub2)/2) ) then
                  exit
             elseif( ub_0k < tail ) then ! ub is too large, so it is the new upper bound
                  ub2 = ub
                  ub2_0k = ub_0k
             else ! ub is too small, so it is the new lower bound
                  ub1 = ub
                  ub1_0k = ub_0k
             endif
        enddo
 endif search_ub


contains

 real(kind=double) pure function bin_dist_0K(p)
 ! PURPOSE: compute the binomial distribution from 0 to k. The number of successes and
 !          samples and factorial are taken directly from the functions that contains
 !          this function
 real(kind=double), intent(IN):: p ! input values, the probability of a success.

 real(kind=double)::  p_at, q_at ! probability at the power of
 integer:: j ! loop counter

 ! initialize with the value probabide by zero successes and then add up to k

 p_at = 0.0_double
 q_at = N*log(1.0_double - p)
 bin_dist_0k = exp(p_at + q_at + lcoeff(0))


 do j = 1, k
     p_at = p_at+log(p)
     q_at = q_at - log(1.0_double - p)
     bin_dist_0k = bin_dist_0k + exp(  p_at + q_at + lcoeff(j)  )
 enddo


end function bin_dist_0k


 real(kind=double) pure function bin_dist_kn(p)
! PURPOSE: compute the binomial distribution from k to N. The number of successes and
 !          samples and factorial are taken directly from the functions that contains
 !          this function
 real(kind=double), intent(IN):: p ! input values, the probability of a success.

 real(kind=double)::  p_at, q_at ! probability at the power of
 integer:: j ! loop counter

! initialize with the value probabide by  k successes
 p_at = k*log(p)
 q_at = (N-k)*log(1.0_double - p)
 bin_dist_kn = exp(  p_at + q_at + lcoeff(k)  )


 do j = k+1, N
     p_at = p_at+log(p)
     q_at = q_at - log(1.0_double - p)
     bin_dist_kn = bin_dist_kn   +  exp(p_at + q_at + lcoeff(j))
 enddo


end function bin_dist_kn


!------------------------------------------------
 end subroutine empirical_k_CI
!------------------------------------------------
!------------------------------------------------------------------------------------

!------------------------------------------------
 pure subroutine empirical_TPF(ms, act_pos, threshold , emp_TPF, ierror)
!------------------------------------------------
!PURPOSE:   Compute the the value of TPF (sensitivity) for a specific threshold
!           positivity is assumed to be for larger values
!ALGORITHM: straight comparison of values, considers ties as  .5

 implicit none

! Declaration of the interface
 integer, intent(IN):: ms ! number of actually-negative cases
 real(kind=double), dimension (ms),intent(IN):: act_pos     ! actually-negative cases
 real(kind=double), intent(IN):: threshold ! value under which the FPF value is desired

 real(kind=double),intent(out):: emp_TPF    ! estimated TPF or sensitivity value
 integer, intent(OUT) :: ierror ! bad_input  <- unacceptable input values

 integer i ! Loop counter

 if(ms <= 0) then
     ierror = bad_input
     emp_TPF = 666.0_double
     return
 else
     ierror = 0
 endif

 emp_TPF  = 0.0_double
 do i=1,ms
           if(act_pos(I) > threshold) then
                  emp_TPF = emp_TPF + 1.0_double
           elseif(act_pos(I) .speq. threshold) then
                  emp_TPF = emp_TPF + 0.5_double
           else
                  ! One should add 0.0_double, i.e., do nothing
           endif
 enddo

  emp_TPF  =  emp_TPF/ real( ms, kind=double )


!------------------------------------------------
 end subroutine empirical_TPF
!------------------------------------------------
!------------------------------------------------------------------------------------

!------------------------------------------------
 pure subroutine empirical_FPF(mn, act_neg, threshold , emp_FPF, ierror)
!------------------------------------------------
!PURPOSE:   Compute the the value of FPF (1 - specificity) for a specific threshold
!           positivity is assumed to be for larger values
!ALGORITHM: straight comparison of values, considers ties as .5

 implicit none

! Declaration of the interface
 integer, intent(IN):: mn ! number of actually-negative cases
 real(kind=double), dimension (mn),intent(IN):: act_neg     ! actually-negative cases
 real(kind=double), intent(IN):: threshold ! value under which the FPF value is desired

 real(kind=double),intent(out):: emp_FPF    ! estimated FPF or 1 - specificity
 integer, intent(OUT) :: ierror ! bad_input  <- unacceptable input values

 integer i ! Loop counter

 if(mn <= 0) then
     emp_FPF = 666.0_double
     ierror = bad_input
     return
 else
     ierror = 0
 endif

 emp_FPF  = 0.0_double
 do i=1,mn
           if(act_neg(I) > threshold) then
                  emp_FPF = emp_FPF + 1.0_double
           elseif(act_neg(I) .speq. threshold) then
                  emp_FPF = emp_FPF + 0.5_double
           else
                  ! One should add 0.0_double, i.e., do nothing
           endif
 enddo
  emp_FPF  =  emp_FPF/ real( mn, kind=double )


!------------------------------------------------
 end subroutine empirical_FPF
!------------------------------------------------
!------------------------------------------------------------------------------------


!------------------------------------------------
 pure subroutine wilcoxon(mn, ms, act_neg, act_pos, wilcoxon_stat, var_wilcoxon)
!------------------------------------------------
! PURPOSE: Compute the  Mann-Whitney form of the  Wilcoxon statistic,
!          which is known to be equivalent to the Area Under the trapezoid
!          determined by the Emprical ROC points. Also called C-statistic.
!ALGORITHM:Two versions implemented:
!          straight comparison (scaling as MN*MS) currently inactive
!          and
!          qsort algorithm scaling as Mn*Log Mn + Ms*Log Ms + Mn+Ms
!
!          Errors are available as Hanley's exponetial-based (active)
!          or as DeLong and DeLong sort-of jackknife, which can also provide
!          analysis for correlated  modalities (in another procedure).


 implicit none

! Declaration of the interface
 integer, intent(IN):: mn, ms ! number of actually-negative and actually-positive
 real(kind=double), dimension (mn),intent(IN):: act_neg     ! actually-negative cases
 real(kind=double), dimension (ms),intent(IN):: act_pos ! actually-positive cases

 real(kind=double),intent(out):: wilcoxon_stat    ! Wilcoxon statistics
 real(kind=double),intent(out):: var_wilcoxon ! standard error of the wilcoxon statistics

! Internal variables. Needed only if placement values are used
! real(kind=double), dimension(mn):: V01 ! Array with the act-negative cases placement values
! real(kind=double), dimension(ms):: V10 ! Array with the act-positive cases placement values

 integer I,J ! Loop counters, needed only for straight comparison


!  call compute_V10_V01(mn, ms, act_neg, act_pos, V01, V10)

!  wilcoxon_stat =  sum(V01)/mn ! sum(V10)/ms should produce the same answer

! The centering is better done before the calculations for stability reasons. This error can
! be shown to be equivalent to some sort of Jackknife.

! se_wilcoxon =      dot_product(V01 - wilcoxon_stat, V01 - wilcoxon_stat)/ ( (mn-1)*ms ) + &
!                    dot_product(V10 - wilcoxon_stat, V10 - wilcoxon_stat)/ ( (ms-1)*mn )

! Straight comparison and Hanley et al. error estimation.

 wilcoxon_stat  = 0.0_double
 DO I=1,mn
   DO J=1,ms
           IF(act_neg(I) < act_pos(J)) THEN
               wilcoxon_stat = wilcoxon_stat + 1.0_double
           ELSEIF(act_neg(I) .speq. act_pos(J)) THEN
               wilcoxon_stat = wilcoxon_stat + 0.5_double
           ELSE
               ! One should add 0.0_double, i.e., do nothing
           ENDIF
    ENDDO
  ENDDO
 wilcoxon_stat = wilcoxon_stat /( real( mn, kind=double )* ms)

 var_wilcoxon = Var_hanley_exponential(wilcoxon_stat, mn, ms)


!------------------------------------------------
 end subroutine wilcoxon
!------------------------------------------------
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
pure real(kind=double) function Var_hanley_exponential(W, mn, ms) result(Var_W)
!------------------------------------------------------------------------------------
!PURPOSE: COMPUTATION OF THE STANDARD ERROR THE WILCONXON STAT. IT FOLLOWS THE FORMULA
!         FOR AN EXPONENTIAL DISTRIBUTION FROM:
!         HANLEY, J.A. AND MCNEIL, B.J. - RADIOLOGY, 143:29-36, 1982.

 integer,intent(IN):: mn, ms     ! number of actually-negative and actually-positive cases
 real(kind=double), intent(IN):: W ! value of the Wilcoxon statistic


! Internal variables
 real (kind=double):: q1,q2 ! quantities related to the mathematical definition of the
                            ! standard error term for the wilcoxon in the original hanley & McNeil paper

 q1 = W / (2.0_double - W)
 q2 = (2.0_double * W**2)  /  (1.0_double + W)

 Var_W =  ( W*(1.0_double-W) + (ms-1)*(q1-W**2) + (mn-1)*(q2-W**2) )  /  (real(mn, kind=double)*ms)

!------------------------------------------------------------------------------------
 end function Var_hanley_exponential
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------
 subroutine gandpboot(mn, ms, num_mod, act_neg, act_pos, des_neg, des_pos, U_vec , U_vec_cov, n_boot)
!--------------------------------------------------------------------------------------------------------
!PURPOSE:  Compute array with the Mann-Whitney form of the Wilcoxon for a set
!          of partially-paired modalities (aka treatments or fixed factors). It also provides the
!          variance covariance matrix for that set (check WARNING below!). Notice that this implementation is for num_mod
!          modalities, i.e., any number.
!          Wilcoxon is known to be equivalent to the Area Under the trapezoid
!          determined by the Emprical ROC points. Also related to the so-called C-statistic.
!          This measure has a number of issues and advantages, be mindful of it.
!ALGORITHM:uses the a set of n_boot bootstrap sets to compute the variance covariance matrix
!          Note that we use sparse matrices that in principle can create some memory problems if the
!          datasets are really huge. In case of need we can adopt a more complex, but compact notation.
!          See below for some possible solutions of memory issues (We do not believe that this will ever be of any
!          relevance).
!          The algorithm could be written as a single double loop over r and s, but I find it more straightforward in this
!          form and I can't see many efficiency gains from making a single loop.
! SOURCE:  Academic Radiology 15,1234
! NOTE:    It works also for fully paired datasets. We kept the general notation
!          of DeLong and DeLong because the notation!
!          modalities are indexed using r and s and cases are indexed using k. We also consider either cases that are in both
!          the modalities or all the cases that are in each modalities.
! NOTE1:   Assumes that for the input data larger values have larger "positivity", i.e., they are more likely to be
!          actually-positive according to the algorithm being tested. If your dataset has not this characterist, just find the
!          elementary transformation necessary to do so. (if you can't find it, it is a good sign that you should not be meddling
!          with this code).
! NOTE2:   r is assumed to be the first index (therefore the only one for runs over only one index)
! WARNING: the variance covariance matrix has in the i>=j elements the variance covariance matrix and in the i< j elements
!          it will have the Var{U_i - U_j). This is done to have a more stable estimate of that variance as opposed to
!          Var{i} + var{j} - 2*cov{i,j}
use statistic_functions
implicit none

integer, intent(IN):: num_mod ! number of modalities or treatments analized.

integer, intent(IN):: mn ! total number of distinct negative cases (i.e., every case that has a value for at least one of
                         ! the num_mod modalities)
integer, intent(IN):: ms ! total number of distinct positive cases (i.e., every case that has a value for at least one of
                         ! the num_mod modalities)

real(kind=double), dimension(mn,num_mod), intent(IN):: act_neg ! actually-negative input data for the two modalities to be analyzed
real(kind=double), dimension(ms,num_mod), intent(IN):: act_pos ! actually-positive input data for the two modalities to be analyzed

! Design matrices. Here we assume that if there are values different from 0 or 1, there is an input error (in general there are
! algorithms that allow the use of different flags for the design matrix, for example to indicate clustering, however, ROCKIT
! cannot make use of them and therefore will not accept them.
integer, dimension(mn,num_mod), intent(IN):: des_neg ! actually-negative design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed
integer, dimension(ms,num_mod), intent(IN):: des_pos ! actually-positive design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed?
Real(kind=double),dimension(num_mod), intent(out):: U_vec    ! Array of Wilcoxon statistics, one per treatment
real(kind=double),dimension(num_mod,num_mod), intent(out):: U_vec_cov ! Variance-covariance matrix of the U_vec in the i>=j
                                                                 ! elements the variance covariance matrix and in the i< j elements
                                                                 ! it will have the Var{U_i - U_j). This is done to
                                                                 ! have a more stable estimate of that variance as opposed to
                                                                 ! Var{i} + var{j} - 2*cov{i,j}

integer, intent(IN):: n_boot ! number of bootstrap sets

! Internal variables:

!complete set description for r and s
integer:: mn_r, mn_rs, mn_s ! number of actually-negative cases in modality r, in both modalities (rs) and  in modality s
integer:: ms_r, ms_rs, ms_s ! number of actually-positive cases in modality r, in both modalities (rs) and  in modality s
real(kind=double), dimension(mn):: act_neg_r ! actually-negative input data for modality r only
real(kind=double), dimension(mn,2):: act_neg_rs! actually-negative input data for modality r AND s
real(kind=double), dimension(mn):: act_neg_s ! actually-negative input data for modality r ( or s the current index)
real(kind=double), dimension(ms):: act_pos_r ! actually-negative input data for modality r only
real(kind=double), dimension(ms,2):: act_pos_rs! actually-negative input data for modality r AND s
real(kind=double), dimension(ms):: act_pos_s ! actually-negative input data for modality r ( or s the current index)

!Bootstrap set description for r and s
integer, dimension(max(mn,ms)):: boot_set
real(kind=double), dimension(mn,2):: act_neg_boot! actually-negative input data for modality r AND s
real(kind=double), dimension(ms,2):: act_pos_boot! actually-negative input data for modality r AND s

real(kind=double), dimension(n_boot):: wr ! Array with the values for W for modality r, for each bootstrap iteration
real(kind=double), dimension(n_boot):: ws ! Array with the values for W for modality s, for each bootstrap iteration
real(kind=double), dimension(n_boot):: w_r_s ! Array with the values for W for r minus W for s, for each bootstrap iteration

real(kind=double):: wr_avg, ws_avg, w_r_s_avg, wr_delta, ws_delta   ! bootstrap averages and deltas
real(kind=double):: s_r, s_s, s_rs, s_r_s   ! var-cov elements and variance of the difference

integer:: r_mod, s_mod, j_case ! Loop counters for modality r, s and case i positive case j negative

integer:: ib

! Compute the values of the wilcoxon statistics for all the modalities
r_loop: do r_mod = 1, num_mod ! Loop over modalities
         s_loop: do s_mod = r_mod, num_mod ! Loop over modalities
              ! initialize for each combination the cases in modality r (or s) will always
              ! be the same, but the ones only in r or in both might change, therefore it is
              ! easier to reset it every time. Cases are separated in the three groups (r,s,rs)
              mn_r = 0; mn_rs = 0; mn_s = 0
              act_neg_r = 0.0_double; act_neg_rs = 0.0_double; act_neg_s = 0.0_double
              ! note that the the act_neg_{r.s} array is compactly stored as opposed to act_neg
              loop_negatives: do j_case = 1, mn
                         ! Check whether the case in in both modalities or else
                         if(des_neg(j_case,r_mod) == 1 .and. des_neg(j_case,s_mod) == 1) then
                                mn_rs = mn_rs + 1
                                act_neg_rs(mn_rs,1) = act_neg(j_case,r_mod)
                                act_neg_rs(mn_rs,2) = act_neg(j_case,s_mod)
                          elseif(des_neg(j_case,r_mod) == 1 .and. des_neg(j_case,s_mod) == 0) then
                                mn_r = mn_r + 1
                                act_neg_r(mn_r) = act_neg(j_case,r_mod)
                         elseif(des_neg(j_case,r_mod) == 0 .and. des_neg(j_case,s_mod) == 1) then
                                mn_s = mn_s + 1
                                act_neg_s(mn_s) = act_neg(j_case,s_mod)
                         endif
              enddo loop_negatives
              ms_r = 0; ms_rs = 0; ms_s = 0
              act_pos_r = 0.0_double; act_pos_rs = 0.0_double; act_pos_s = 0.0_double
              ! note that the the act_pos_{r,s} array is compactly stored as opposed to act_pos
              loop_positives: do j_case = 1, ms
                         ! Check whether the case in in both modalities or else
                         if(des_pos(j_case,r_mod) == 1 .and. des_pos(j_case,s_mod) == 1) then
                                ms_rs = ms_rs + 1
                                act_pos_rs(ms_rs,1) = act_pos(j_case,r_mod)
                                act_pos_rs(ms_rs,2) = act_pos(j_case,s_mod)
                         elseif(des_pos(j_case,r_mod) == 1 .and. des_pos(j_case,s_mod) == 0) then
                                ms_r = ms_r + 1
                                act_pos_r(ms_r) = act_pos(j_case,r_mod)
                         elseif(des_pos(j_case,r_mod) == 0 .and. des_pos(j_case,s_mod) == 1) then
                                ms_s = ms_s + 1
                                act_pos_s(ms_s) = act_pos(j_case,s_mod)
                         endif
              enddo loop_positives
              ! COMPUTE THE VALUES FOR THE DATASET WITHOUT BOOTSTRAP
              if(r_mod == s_mod) then
                    ! Load the remaining data to compute the wilcoxon, since usually the paired part is the largest
                    ! to make it run faster we load everything in the paired arrays
                    act_neg_rs(mn_rs+1:mn_rs+1+mn_s,2) = act_neg_s(1:mn_s)
                    act_pos_rs(ms_rs+1:ms_rs+1+ms_s,2) = act_pos_s(1:ms_s)
                    call fast_wilcoxon(mn_s+mn_rs, ms_s+ms_rs, act_neg_rs(1:mn_rs+mn_s,2), &
                                       act_pos_rs(1:ms_rs+ms_s,2), U_vec(s_mod))
              endif

              ! BOOTSTRAP LOOP
              do ib=1,n_boot
                     ! NEGATIVES
                     ! build the only r part
                     call bootstrap_set(1, mn_r, mn_r, boot_set(1:mn_r))
                     forall(j_case=1:mn_r) act_neg_boot(j_case,1) = act_neg_r(boot_set(j_case))
                     ! build the r,s part
                     call bootstrap_set(1, mn_rs, mn_rs, boot_set(1:mn_rs))
                     forall(j_case=1:mn_rs) act_neg_boot(mn_r+j_case,1) = act_neg_rs(boot_set(j_case),1)
                     forall(j_case=1:mn_rs) act_neg_boot(j_case,2) = act_neg_rs(boot_set(j_case),2)
                     ! build the only s part
                     call bootstrap_set(1, mn_s, mn_s, boot_set(1:mn_s))
                     forall(j_case=1:mn_s) act_neg_boot(mn_rs+j_case,2) = act_neg_s(boot_set(j_case))

                     ! POSITIVES
                     ! build the only r part
                     call bootstrap_set(1, ms_r, ms_r, boot_set(1:ms_r))
                     forall(j_case=1:ms_r) act_pos_boot(j_case,1) = act_pos_r(boot_set(j_case))
                     ! build the r,s part
                     call bootstrap_set(1, ms_rs, ms_rs, boot_set(1:ms_rs))
                     forall(j_case=1:ms_rs) act_pos_boot(ms_r+j_case,1) = act_pos_rs(boot_set(j_case),1)
                     forall(j_case=1:ms_rs) act_pos_boot(j_case,2) = act_pos_rs(boot_set(j_case),2)
                     ! build the only s part
                     call bootstrap_set(1, ms_s, ms_s, boot_set(1:ms_s))
                     forall(j_case=1:ms_s) act_pos_boot(ms_rs+j_case,2) = act_pos_s(boot_set(j_case))

                     call fast_wilcoxon(mn_r+mn_rs, ms_r+ms_rs, act_neg_boot(1:mn_r+mn_rs,1), &
                                       act_pos_boot(1:ms_r+ms_rs,1), wr(ib))

                     call fast_wilcoxon(mn_rs+mn_s, ms_rs+ms_s, act_neg_boot(1:mn_rs+mn_s,2), &
                                       act_pos_boot(1:ms_rs+ms_s,2), ws(ib))
                     w_r_s(ib) = wr(ib) - ws(ib)
              enddo

              wr_avg = sum(wr)/n_boot
              ws_avg = sum(ws)/n_boot
              w_r_s_avg = sum(w_r_s)/n_boot

              s_r   = 0.0_double
              s_s   = 0.0_double
              s_rs  = 0.0_double
              s_r_s = 0.0_double
              ! this could be rendered more stable by randking or partitioning the differences
              ! thus far this did not seem to be necessary.
              do ib=1,n_boot
                     wr_delta = wr(ib) - wr_avg
                     ws_delta = ws(ib) - ws_avg
                     s_r_s = s_r_s + (w_r_s(ib) - w_r_s_avg)**2
                     s_r = s_r + wr_delta**2
                     s_s = s_s + ws_delta**2
                     s_rs = s_rs + wr_delta*ws_delta
              enddo
              if(r_mod == s_mod) U_vec_cov(r_mod,r_mod) = s_r/(n_boot-1)
              U_vec_cov(s_mod,r_mod) = s_rs/(n_boot-1)
              if( r_mod /= s_mod) U_vec_cov(r_mod,s_mod) =  s_r_s/(n_boot-1)

         enddo s_loop
enddo r_loop

!---------------------------------------------------------------------------
end subroutine gandpboot
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

!------------------------------------------------------------------------------------
 subroutine fast_wilcoxon(mn, ms, Y, X, wilc)
!-------------------------------------------------------------------------------------------
! PURPOSE  : Compute the wilcoxon statistics very quickly, returns not SE, just the wilcoxon.
! ALGORITHM:Uses quick sort algorithm to create two ranked sequences, one for X (the positives) and one for
!           Y (the negatives) -- Notation follows approximately the original article. Then the two sequences are
!           scanned (together) from the smallest value to the largest value.
! NOTE:     Positivity is assumed to be for larger values, one can submit -X and -Y if this is not the case.
! NOTE1:    Here we follow the The Long and the Long notation and give 10 to the positives and 01 to the negatives.
!           I usually prefer to do the opposite so be mindful when looking at other functions.
! SOURCE:   Some pages of notes that I have already shredded.
! AUTHOR:   LP U CHICAGO 2008 CE
!-------------------------------------------------------------------------------------------

! Interface definition
integer, intent(IN):: mn, ms ! Number of act-negative and of act-positive cases, respectively
real(kind=double), intent(IN), dimension(mn):: Y ! Array with the act-negative cases
real(kind=double), intent(IN), dimension(ms):: X ! Array with the act-positive cases
real(kind=double), intent(OUT):: wilc

! Internal variables
integer, dimension(mn):: index_Y ! Array with the sorting indexes of the act-neg cases, i.e., the first element of the
                                 ! array contains the location of the smallest value in the array X.
integer, dimension(ms):: index_X ! Array with the sorting indexes of the act-pos cases
integer:: i_neg, i_pos ! value of the current count for the act-negative and act-positive cases arrays
integer:: index_i_neg, index_i_pos ! location of the current counter index for the negative and positive cases
integer:: i_neg_old, i_pos_old ! counters for the negative and positive cases used to locate the beginning of a tie between act-pos and
                               ! act-neg cases (see body of procedure to see why we do this).
logical:: prev_pos ! previous increment was a positive

! Sort act-negative and act-positive cases and build the sorting index arrays.
call qsortd(Y, index_Y, mn)
call qsortd(X, index_X, ms)

! Start from the smallest values. We will count how many act-neg are smaller or equal then the current act-pos and
! by reciprocity how many act-pos are larger or equal than the current act-neg.
! We  work the cases from smallest to largest using the indexes i_neg and i_pos to know which rank are we at and the
! convert to the  indices index_i_neg and index_i_pos tell us where the corresponding values are and which V01 and V10
! functions to update. In this way the resulting functions will have the same ordering as the input arrays.

! Initialize variables used in the loop to the first two cases.
i_neg = 1 ; index_i_neg = index_Y(i_neg)
i_pos = 1;  index_i_pos = index_X(i_pos)
i_neg_old = 1
i_pos_old = 1
wilc = 0.0_double
prev_pos = .false.



! Endless loop, exit is produced by reaching the end of both arrays
do_wilc: do
! First case: the act-positive is larger than the act-negative
      if( X(index_i_pos) > Y(index_i_neg)) then
               ! find ALL the remaining negatives that are smaller than this positive
               i_neg_old = i_neg
               i_neg = i_neg + 1
               find_smaller_negs: do
                  if( i_neg > mn ) then
                        i_neg = i_neg - 1 ! go back to the previous which has to be the last one in the tie-stretch
                        exit find_smaller_negs
                  elseif( X(index_i_pos) <= Y(index_Y(i_neg)) ) then
                        i_neg = i_neg - 1 ! go back to the previous which has to be the last one in the tie-stretch
                        exit find_smaller_negs
                  else
                       i_neg = i_neg + 1
                  endif
               enddo find_smaller_negs

               if(i_neg == mn) then
                 wilc = wilc + (ms-i_pos+1)*real(mn,kind=double)
                exit do_wilc
             else
                 wilc = wilc + real(i_neg,kind=double)
                 ! increase the counter for the act-negative and update its value.
                 i_neg = i_neg + 1  ! it can't be a tie from the previous stretch because tie-stretches are loaded completely.
                 index_i_neg = index_Y(i_neg) ! update the "current negative index"
             endif
             prev_pos = .true.
      ! Start of a tie-stretch, act-pos and act-neg cases have the same value
      elseif( X(index_i_pos) == Y(index_i_neg)) then
             ! If they are equal we will keep incrementing until all the cases with values tied to these ones have been located.
             ! First load the current values as starting the starting point of the tie-stretch.
             i_pos_old = i_pos
             i_pos = i_pos + 1
             ! Find the end of the tie-stretch for both the act-negative and the act-positive
             find_equals_pos: do
                 if( i_pos > ms ) then
                        i_pos = i_pos - 1 ! go back to the previous which has to be the last one in the tie-stretch
                        exit find_equals_pos
                 elseif( X(index_X(i_pos)) > Y(index_i_neg)) then
                        i_pos = i_pos - 1 ! go back to the previous which has to be the last one in the tie-stretch
                        exit find_equals_pos
                 else ! move one more up
                       i_pos = i_pos + 1
                 endif
             enddo find_equals_pos
             i_neg_old = i_neg
             i_neg = i_neg + 1
             find_equals_neg: do
                 if( i_neg > mn ) then
                       i_neg = i_neg - 1 ! go back to the previous which has to be the last one in the tie-stretch
                       exit find_equals_neg
                 elseif( X(index_X(i_pos)) < Y(index_Y(i_neg))) then
                       i_neg = i_neg - 1 ! go back to the previous which has to be the last one in the tie-stretch
                       exit find_equals_neg
                  else
                       i_neg = i_neg + 1
                  endif
             enddo find_equals_neg
             ! Add the the ties .5*(number of positive tied)*(number of negative tied), plus all the 1 comparisons,
             ! which is equal too (number of positive tied)*(number of negative before the tie == i_neg_old - 1)
             ! However, if before the tie we had a positive, we need to remove some double counting
             if(prev_pos) then
                 wilc = wilc + 0.5_double*real(i_pos + 1 - i_pos_old, kind=double)* &
                        real(i_neg - 1 + i_neg_old, kind=double) - i_neg_old + 1
             else
                 wilc = wilc + 0.5_double*(i_pos + 1 - i_pos_old)*(i_neg - 1 + i_neg_old)
             endif

             ! All the ties between act-neg and act-pos are exhausted, we move to a larger positive, if possible
             if(i_pos < ms) then
                  ! Move to the next act-positive, update the current act-positive locator.
                  i_pos = i_pos + 1
                  index_i_pos = index_X(i_pos) ! update the index
                  !Locate the top negative in the tie-stretch
                  index_i_neg = index_Y(i_neg) !
             else ! We swept all the positives and all the negatives remaining negatives are larger thus inv-placement is zero
                 exit do_wilc
             endif
             prev_pos = .false.
      else
             i_pos_old = i_pos
             i_pos = i_pos + 1
             ! Find the end of the tie-stretch for both the act-negative and the act-positive
             find_smaller_pos: do
                 if( i_pos > ms ) then
                        exit do_wilc ! we are done here
                 elseif( X(index_X(i_pos)) >= Y(index_i_neg)) then
                         exit find_smaller_pos
                 else ! move one more up
                       wilc = wilc + real(i_neg-1,kind=double)
                       i_pos = i_pos + 1
                 endif
             enddo find_smaller_pos
             index_i_pos = index_X(i_pos) ! update the index
             prev_pos = .false.
      endif
end do do_wilc

wilc = wilc / ( real(mn, kind=double)*ms)


 end subroutine fast_wilcoxon
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
 subroutine m_mod_one_shot(mn, ms, num_mod, act_neg, act_pos, des_neg, des_pos, U_vec , U_vec_cov)
!---------------------------------------------------------------------------
!PURPOSE:  Compute array with the Mann-Whitney form of the Wilcoxon for a set
!          of partially-paired modalities (aka treatments or fixed factors). It also provides the
!          variance covariance matrix for that set. Notice that this implementation is for num_mod
!          modalities, i.e., any number. It is based on the moment approach proposed by Gallas.
!          for the multi-reader-multi-case scenario.
!          Wilcoxon is known to be equivalent to the Area Under the trapezoid
!          determined by the Emprical ROC points. Also related to the so-called C-statistic.
!          This measure has a number of issues and advantages, be mindful of it.
!ALGORITHM:uses the function compute_V10_V01(Mn, Ms, act_neg, act_pos, V01, V10) to
!          to compute the placement values for the act-positive and 1- placement
!          for the act-negative cases for each modality. It uses the results from
!          the Gallas and Brown paper to estimate the variance-covariance matrix.
!          Note that we use sparse matrices that in principle can create some memory problems if the
!          datasets are really huge. In case of need we can adopt a more complex, but compact notation.
!          See below for some possible solutions of memory issues (We do not believe that this will ever be of any
!          relevance).
!          The algorithm could be written as a single double loop over r and s, but I find it more straightforward in this
!          form and I can't see many efficiency gains from making a single loop.
! SOURCE:  Brandon D. Gallas, David G. Brown Neural Networks 21 (2008) 387-397 and some pages of notes that I have
!          already shredded and incinerated.
! NOTE:    It works also for fully paired datasets. We kept the general notation
!          of DeLong and DeLong to have all the partially-paired, non-parametric models based on the same notation. Moreover, while
!          the notation used by Gallas et al. is more elegant, it is less numerically efficient. Therefore
!          modalities are indexed using r and s and cases are indexed using k. We rewrote the equations so that we consider only cases
!          that are in both the modalities or all the cases that are in each modalities.
! NOTE1:   Assumes that for the input data larger values have larger "positivity", i.e., they are more likely to be
!          actually-positive according to the algorithm being tested. If your dataset has not this characterist, just find the
!          elementary transformation necessary to do so. (if you can't find it, it is a good sign that you should not be meddling
!          with this code).
! NOTE2:   r is assumed to be the first index (therefore the only one for runs over only one index)

integer, intent(IN):: num_mod ! number of modalities or treatments analized.

integer, intent(IN):: mn ! total number of distinct negative cases (i.e., every case that has a value for at least one of
                         ! the num_mod modalities)
integer, intent(IN):: ms ! total number of distinct positive cases (i.e., every case that has a value for at least one of
                         ! the num_mod modalities)

real(kind=double), dimension(mn,num_mod), intent(IN):: act_neg ! actually-negative input data for the two modalities to be analyzed
real(kind=double), dimension(ms,num_mod), intent(IN):: act_pos ! actually-positive input data for the two modalities to be analyzed

! Design matrices. Here we assume that if there are values different from 0 or 1, there is an input error (in general there are
! algorithms that allow the use of different flags for the design matrix, for example to indicate clustering, however, ROCKIT
! cannot make use of them and therefore will not accept them.
integer, dimension(mn,num_mod), intent(IN):: des_neg ! actually-negative design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed
integer, dimension(ms,num_mod), intent(IN):: des_pos ! actually-positive design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed
real(kind=double),dimension(num_mod), intent(out):: U_vec    ! Array of Wilcoxon statistics, one per treatment
real(kind=double),dimension(num_mod,num_mod), intent(out):: U_vec_cov    ! Variance-covariance matrix of the U_vec

! Internal variables.
integer:: mn_rs, ms_rs ! number of actually-negative and actually-positive in both modalities (rs) and  in modality s
integer, dimension(num_mod):: mn_mod! number of actually-negative cases in each modality
integer, dimension(num_mod):: ms_mod ! number of actually-positive cases in each modality

real(kind=double), dimension(mn):: act_neg_r ! actually-negative input data for modality r ( or s the current index)
real(kind=double), dimension(ms):: act_pos_r ! actually-positive input data for modality r (the current index)


real(kind=double), dimension(mn,num_mod):: V01 ! Array with the act-negative cases placement values, by modality, it is a partially
                 ! empty arrays, because not all modalities will have all cases, but unless the number of modalities and cases is very large
                 ! this should not be an issue. It is stored densely within each modality
real(kind=double), dimension(ms,num_mod):: V10 ! Array with the act-positive cases placement values, by treatment. Look at V01 for
                 ! more information.

real(kind=double), dimension(mn):: RP_V01 ! Array with the scalar products of the V01s, before ranking -- done to stabilize variance calcs
real(kind=double), dimension(ms):: RP_V10 ! Array with the scalar products of the V01s, before ranking -- done to stabilize variance calcs
real(kind=double):: SrSs_rs ! Scalar product of matrix Sr times Ss over the values they have in common
integer, dimension(mn):: Ind_RP_V01 ! Ranking indices used in the stabilization loop
integer, dimension(ms):: Ind_RP_V10 ! Ranking indices used in the stabilization loop
real(kind=double):: V10r_V10s, V01r_V01s ! the act-negative and act-positive components of the variance.


integer:: r_mod, s_mod, i_case, j_case ! Loop counters for modality r, s and case i positive case j negative
integer:: i_neg_r, i_pos_r, i_neg_s, i_pos_s ! counters for cases in modalities
real(kind=double):: mrs, mrs_rs, cmn, cms, c4m ! constants for computation of variances, from Gallas and Brown and some calculations

! Compute the placement values for all the modments/modalities
! Note that for each modality, the placement values are stored in compact form, i.e., the
! all the cases in that modality follow each other in the array, independently from their
! ordering in the possibly sparse act_neg and act_pos arrays. These arrays could be made
! smaller either by constructing a type which is made of arrays of different sizes (for each
! modality) or find the largest one. Since this might affect portability, not to mention
! understandability I preferred to leave a potentially sparse array. If you have
! memory issues, you can take either of these two routes or a different one.
V01 = 0.0_double
V10 = 0.0_double
mn_mod = 0
ms_mod = 0

do r_mod = 1, num_mod ! Loop over modalities
         ! Count how many actually-negative cases are in modality r_mod and load them
         act_neg_r = 0.0_double
         act_pos_r = 0.0_double
         do j_case  = 1, mn ! Loop over all the different cases present in the dataset
         ! note that the the act_neg_r array is compactly stored as opposed to act_neg
              if(des_neg(j_case,r_mod) == 1) then
                   mn_mod(r_mod) = mn_mod(r_mod) + 1
                   act_neg_r(mn_mod(r_mod)) = act_neg(j_case,r_mod)
              endif
         enddo
         ! Count how many actually-positive cases are in modality r_mod and load them
         do i_case  = 1, ms
              if(des_pos(i_case,r_mod) == 1) then
                   ms_mod(r_mod) = ms_mod(r_mod) + 1
                   act_pos_r(ms_mod(r_mod)) = act_pos(i_case,r_mod)
              endif
          enddo
          ! Compute the placement values
          call compute_V10_V01(mn_mod(r_mod), ms_mod(r_mod), act_neg_r(1:mn_mod(r_mod)), act_pos_r(1:ms_mod(r_mod)),&
                               V01(1:mn_mod(r_mod),r_mod), V10(1:ms_mod(r_mod),r_mod))
          ! Compute the U statistics.
          U_vec(r_mod) =  sum( V10(1:ms_mod(r_mod),r_mod) )/ms_mod(r_mod) ! sum(V10)/ms_r should produce the same answer
enddo

! Compute the variance covariance matrix according to Zhou and Gatsonis.
! We added a centering step for stability reasons.
! This variance estimation can be shown to be equivalent to some sort of Jackknife.

! Notice that some summations over the cases in modality r could be moved to the external loop. Given that
! the important loop for the variance is done over the cases that are in common between the two modalities
! r and s, it does not seem worthy the additional complexity and might be easily counterproductive.
r_loop: do r_mod = 1, num_mod
          s_loop: do s_mod = r_mod, num_mod
                  SrSs_rs = 0.0_double ! set the scalar product of the success matrix
                  ! build the cross products for the actually-negative, taking into account, only for the cases that
                  ! are in common between the two modalities
                  mn_rs = 0 ; i_neg_r = 0; i_neg_s = 0
                  loop_negatives: do j_case = 1, mn
                         ! Check whether the case in in both modalities, which is the only instance in which the
                         ! cases enter into the covariance calculation -- see reference
                         if(des_neg(j_case,r_mod) == 1 .and. des_neg(j_case,s_mod) == 1) then
                                mn_rs = mn_rs + 1
                                ! The number of cases in each modality (independently from the other one) needs to be updated
                                ! to locate the correct placement value in the V01 array.
                                i_neg_r  = i_neg_r + 1
                                i_neg_s  = i_neg_s + 1
                                ! values are first stored in this array, which is later ranked for stability reasons, it has as
                                ! many elements as the number of cases in common to both modalities.
                                RP_V01(mn_rs)= V01(i_neg_r,r_mod) * V01(i_neg_s,s_mod)


                                ! compute the scalar product of the success matrices More care should be paid to the equal signs
                                ! because conversion and I/O can play tricks. However, there is no simple
                                ! solution to this problem so for now we hope people do I/O sensibly.

                                loop_positives_SrSs: do i_case = 1, ms
                                       if(des_pos(i_case,r_mod) == 1 .and. des_pos(i_case,s_mod) == 1) then

                                           ! For both modalities the positive case is larger than the negative case,
                                           ! on average the most likely events (most modalities are better than random)
                                           if( act_pos(i_case,r_mod)  > act_neg(j_case,r_mod)               &
                                               .and.                                                        &
                                               act_pos(i_case,s_mod)  > act_neg(j_case,s_mod)               &
                                              ) then

                                                       SrSs_rs = SrSs_rs + 1.0_double
                                           ! For one modalities the positive case is larger than the negative case, for
                                           ! the other they are the same.
                                           elseif(        ( act_pos(i_case,r_mod)  > act_neg(j_case,r_mod)   &
                                                    .and.                                                    &
                                                        act_pos(i_case,s_mod)  == act_neg(j_case,s_mod)    ) &
                                               .or.                                                          &
                                                      ( act_pos(i_case,r_mod)  == act_neg(j_case,r_mod)      &
                                                     .and.                                                   &
                                                        act_pos(i_case,s_mod)  > act_neg(j_case,s_mod)    )  &
                                               ) then

                                                       SrSs_rs = SrSs_rs + 0.5_double
                                           ! both both modalies the positive and negative cases have the same value
                                           elseif(     act_pos(i_case,r_mod) == act_neg(j_case,r_mod)         &
                                                   .and.                                                      &
                                                      act_pos(i_case,s_mod)  == act_neg(j_case,s_mod)         &
                                                  ) then

                                                       SrSs_rs = SrSs_rs + 0.25_double
                                           endif
                                    endif
                                enddo loop_positives_SrSs

                         ! the value does not enter in the summation to compute the variance, however, we need to update the
                         ! counter for the V01 array
                         elseif(des_neg(j_case,r_mod) == 1 .and. des_neg(j_case,s_mod) == 0) then
                                i_neg_r = i_neg_r + 1
                         elseif(des_neg(j_case,r_mod) == 0 .and. des_neg(j_case,s_mod) == 1) then
                                i_neg_s = i_neg_s + 1
                         endif
                  enddo loop_negatives
                  ! build the cross products for the actually-positive, taking into account, only for the cases that
                  ! are in common between the two modalities
                  i_pos_r  = 0 ; ms_rs = 0 ; i_pos_s  = 0
                  loop_positives: do i_case = 1, ms
                         ! Check whether the case in in both modalities, which is the only instance in which the
                         ! cases enter into the covariance calculation -- see reference
                         if(des_pos(i_case,r_mod) == 1 .and. des_pos(i_case,s_mod) == 1) then
                                ms_rs = ms_rs + 1
                                ! The number of cases in each modality (independently from the other one) needs to be updated
                                ! to locate the correct placement value in the V01 array.
                                i_pos_r  = i_pos_r + 1
                                i_pos_s  = i_pos_s + 1
                                ! values are first stored in this array, which is later ranked for stability reasons, it has as
                                ! many elements as the number of cases in common to both modalities.
                                RP_V10(ms_rs)= V10(i_pos_r,r_mod) * V10(i_pos_s,s_mod)
                         ! the value does not enter in the summation to compute the variance, however, we need to update the
                         ! counter for the V01 array
                         elseif(des_pos(i_case,r_mod) == 1 .and. des_pos(i_case,s_mod) == 0) then
                                i_pos_r = i_pos_r + 1
                         elseif(des_pos(i_case,r_mod) == 0 .and. des_pos(i_case,s_mod) == 1) then
                                i_pos_s = i_pos_s + 1
                         endif
                  enddo loop_positives
                  ! Rank the cross products, this is a stabilization step added for numerical reasons. It is usually unnecessary
                  ! but better safe than sorry.
                  call qsortd(RP_V01(1:mn_rs), Ind_RP_V01(1:mn_rs), mn_rs)
                  call qsortd(RP_V10(1:ms_rs), Ind_RP_V10(1:ms_rs), ms_rs)
                  ! Sum variance components from the smallest to the largest, so we use a goto to keep ordering
                  V01r_V01s = 0.0_double
                  do j_case = 1, mn_rs
                        V01r_V01s = V01r_V01s + RP_V01( Ind_RP_V01(j_case) )
                  enddo
                  V10r_V10s = 0.0_double
                  do i_case = 1, ms_rs
                        V10r_V10s = V10r_V10s + RP_V10( Ind_RP_V10(i_case))
                  enddo


                  ! Combine terms, first compute one of the
                  mrs    = ( mn_mod(r_mod) * real(mn_mod(s_mod), kind=double) ) * &
                           ( ms_mod(r_mod) * real(ms_mod(s_mod), kind=double) )
                  mrs_rs = mn_rs*real(ms_rs,kind=double)
                  cmn    = mn_mod(r_mod) * real(mn_mod(s_mod),kind=double)
                  cms    = ms_mod(r_mod) * real(ms_mod(s_mod),kind=double)
                  c4m    = -(mn_rs*cms + ms_rs* cmn  - mrs_rs)/                   &
                            (mrs * (mrs - mn_rs*cms - ms_rs* cmn + mrs_rs) )
                  U_vec_cov(r_mod, s_mod) =                                                      &
                                            + V01r_V01s * ( 1.0_double / cmn - cms  * c4m )      &
                                            + V10r_V10s * ( 1.0_double / cms - cmn  * c4m )      &
                                            + U_vec(r_mod)* U_vec(s_mod) * mrs *c4m              &
                                            + SrSs_rs * ( c4m - 1.0_double/mrs)
                  ! load the value to the symmetrical location.
                  U_vec_cov(s_mod, r_mod) =  U_vec_cov(r_mod, s_mod)
        enddo s_loop

 enddo r_loop

!---------------------------------------------------------------------------
end subroutine m_mod_one_shot
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
 subroutine ZhouAndGatsonis(mn, ms, num_mod, act_neg, act_pos, des_neg, des_pos, U_vec , U_vec_cov)
!---------------------------------------------------------------------------
!PURPOSE:  Compute array with the Mann-Whitney form of the Wilcoxon for a set
!          of partially-paired modalities (aka treatments or fixed factors). It also provides the
!          variance covariance matrix for that set. Notice that this implementation is for num_mod
!          modalities, i.e., any number, as opposed to the original 2.
!          Wilcoxon is known to be equivalent to the Area Under the trapezoid
!          determined by the Emprical ROC points. Also related to the so-called C-statistic.
!          This measure has a number of issues and advantages, be mindful of it.
!ALGORITHM:uses the function compute_V10_V01(Mn, Ms, act_neg, act_pos, V01, V10) to
!          to compute the placement values for the act-positive and 1- placement
!          for the act-negative cases for each modality. It uses the results from
!          the Zhou and Gatsonis paper to estimate the variance-covariance matrix.
!          Note that we use sparse matrices that in principle can create some memory problems if the
!          datasets are really huge. In case of need we can adopt a more complex, but compact notation.
!          See below for some possible solutions of memory issues (We do not believe that this will ever be of any
!          relevance).
!          The algorithm could be written as a single double loop over r and s, but I find it more straightforward in this
!          form and I can't see many efficiency gains from making a single loop.
! SOURCE:  Zhou and Gatsonis Statistics in Medicine, Vol 15, 1687-1693 (1996) and some pages of notes that I have
!          already shredded and incinerated.
! NOTE:    It works also for fully paired datasets, where it is identical to De Long and DeLong. We kept the general notation
!          of DeLong and DeLong because the notation of Zhou and Gatsonis is both more confusing and less general. Therefore
!          modalities are indexed using r and s and cases are indexed using k. We also consider either cases that are in both
!          the modalities or all the cases that are in each modalities. Zhou and Gastonis consider the cases in both modalities
!          or only in each of them. We find our notation simpler and much more straightforward to implement (not to mention it
!          makes the algorithm look a lot more similar to DeLong and DeLong).
! NOTE1:   Assumes that for the input data larger values have larger "positivity", i.e., they are more likely to be
!          actually-positive according to the algorithm being tested. If your dataset has not this characterist, just find the
!          elementary transformation necessary to do so. (if you can't find it, it is a good sign that you should not be meddling
!          with this code).
! NOTE2:   r is assumed to be the first index (therefore the only one for runs over only one index)

integer, intent(IN):: num_mod ! number of modalities or treatments analized.

integer, intent(IN):: mn ! total number of distinct negative cases (i.e., every case that has a value for at least one of
                         ! the num_mod modalities)
integer, intent(IN):: ms ! total number of distinct positive cases (i.e., every case that has a value for at least one of
                         ! the num_mod modalities)

real(kind=double), dimension(mn,num_mod), intent(IN):: act_neg ! actually-negative input data for the two modalities to be analyzed
real(kind=double), dimension(ms,num_mod), intent(IN):: act_pos ! actually-positive input data for the two modalities to be analyzed

! Design matrices. Here we assume that if there are values different from 0 or 1, there is an input error (in general there are
! algorithms that allow the use of different flags for the design matrix, for example to indicate clustering, however, ROCKIT
! cannot make use of them and therefore will not accept them.
integer, dimension(mn,num_mod), intent(IN):: des_neg ! actually-negative design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed
integer, dimension(ms,num_mod), intent(IN):: des_pos ! actually-positive design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed
real(kind=double),dimension(num_mod), intent(out):: U_vec    ! Array of Wilcoxon statistics, one per treatment
real(kind=double),dimension(num_mod,num_mod), intent(out):: U_vec_cov    ! Variance-covariance matrix of the U_vec

! Internal variables.
integer:: mn_r, mn_rs, mn_s ! number of actually-negative cases in modality r, in both modalities (rs) and  in modality s
integer:: ms_r, ms_rs, ms_s ! number of actually-positive cases in modality r, in both modalities (rs) and  in modality s

real(kind=double), dimension(mn):: act_neg_r ! actually-negative input data for modality r ( or s the current index)
real(kind=double), dimension(ms):: act_pos_r ! actually-positive input data for modality r (the current index)


real(kind=double), dimension(mn,num_mod):: V01 ! Array with the act-negative cases placement values, by modality, it is a partially
                 ! empty arrays, because not all modalities will have all cases, but unless the number of modalities and cases is very large
                 ! this should not be an issue. It is stored densely within each modality
real(kind=double), dimension(ms,num_mod):: V10 ! Array with the act-positive cases placement values, by treatment. Look at V01 for
                 ! more information.

real(kind=double), dimension(mn):: RP_V01 ! Array with the scalar products of the V01s, before ranking -- done to stabilize variance calcs
real(kind=double), dimension(ms):: RP_V10 ! Array with the scalar products of the V01s, before ranking -- done to stabilize variance calcs
integer, dimension(mn):: Ind_RP_V01 ! Ranking indices used in the stabilization loop
integer, dimension(ms):: Ind_RP_V10 ! Ranking indices used in the stabilization loop
real(kind=double):: eps_10, eps_01 ! the act-negative and act-positive components of the variance.

integer:: r_mod, s_mod, i_case, j_case ! Loop counters for modality r, s and case i positive case j negative

! Compute the placement values for all the modments/modalities
! Note that for each modality, the placement values are stored in compact form, i.e., the
! all the cases in that modality follow each other in the array, independently from their
! ordering in the possibly sparse act_neg and act_pos arrays. These arrays could be made
! smaller either by constructing a type which is made of arrays of different sizes (for each
! modality) or find the largest one. Since this might affect portability, not to mention
! understandability I preferred to leave a potentially sparse array. If you have
! memory issues, you can take either of these two routes or a different one.
V01 = 0.0_double
V10 = 0.0_double
do r_mod = 1, num_mod ! Loop over modalities
         ! Count how many actually-negative cases are in modality r_mod and load them
         mn_r = 0 ! initialize
         act_neg_r = 0.0_double
         act_pos_r = 0.0_double
         do j_case  = 1, mn ! Loop over all the different cases present in the dataset
         ! note that the the act_neg_r array is compactly stored as opposed to act_neg
              if(des_neg(j_case,r_mod) == 1) then
                   mn_r = mn_r + 1
                   act_neg_r(mn_r) = act_neg(j_case,r_mod)
              endif
         enddo
         ! Count how many actually-positive cases are in modality r_mod and load them
         ms_r = 0
         do i_case  = 1, ms
              if(des_pos(i_case,r_mod) == 1) then
                   ms_r = ms_r + 1
                   act_pos_r(ms_r) = act_pos(i_case,r_mod)
              endif
          enddo
          ! Compute the placement values
          call compute_V10_V01(mn_r, ms_r, act_neg_r(1:mn_r), act_pos_r(1:ms_r), V01(1:mn_r,r_mod), V10(1:ms_r,r_mod))
          ! Compute the U statistics.
          U_vec(r_mod) =  sum( V10(1:ms_r,r_mod) )/ms_r ! sum(V10)/ms_r should produce the same answer
enddo

! Compute the variance covariance matrix according to Zhou and Gatsonis.
! We added a centering step for stability reasons.
! This variance estimation can be shown to be equivalent to some sort of Jackknife.

! Notice that some summations over the cases in modality r could be moved to the external loop. Given that
! the important loop for the variance is done over the cases that are in common between the two modalities
! r and s, it does not seem worthy the additional complexity and might be easily counterproductive.
r_loop: do r_mod = 1, num_mod
          s_loop: do s_mod = r_mod, num_mod
                  ! build the cross products for the actually-negative, taking into account, only for the cases that
                  ! are in common between the two modalities
                  mn_r  = 0 ; mn_rs = 0 ; mn_s  = 0
                  loop_negatives: do j_case = 1, mn
                         ! Check whether the case in in both modalities, which is the only instance in which the
                         ! cases enter into the covariance calculation -- see reference
                         if(des_neg(j_case,r_mod) == 1 .and. des_neg(j_case,s_mod) == 1) then
                                mn_rs = mn_rs + 1
                                ! The number of cases in each modality (independently from the other one) needs to be updated
                                ! to locate the correct placement value in the V01 array.
                                mn_r  = mn_r + 1
                                mn_s  = mn_s + 1
                                ! values are first stored in this array, which is later ranked for stability reasons, it has as
                                ! many elements as the number of cases in common to both modalities.
                                RP_V01(mn_rs)= ( V01(mn_r,r_mod) - U_vec(r_mod) ) * ( V01(mn_s,s_mod) - U_vec(s_mod) )
                         ! the value does not enter in the summation to compute the variance, however, we need to update the
                         ! counter for the V01 array
                         elseif(des_neg(j_case,r_mod) == 1 .and. des_neg(j_case,s_mod) == 0) then
                                mn_r = mn_r + 1
                         elseif(des_neg(j_case,r_mod) == 0 .and. des_neg(j_case,s_mod) == 1) then
                                mn_s = mn_s + 1
                         endif
                  enddo loop_negatives
                  ! build the cross products for the actually-positive, taking into account, only for the cases that
                  ! are in common between the two modalities
                  ms_r  = 0 ; ms_rs = 0 ; ms_s  = 0
                  loop_positives: do i_case = 1, ms
                         ! Check whether the case in in both modalities, which is the only instance in which the
                         ! cases enter into the covariance calculation -- see reference
                         if(des_pos(i_case,r_mod) == 1 .and. des_pos(i_case,s_mod) == 1) then
                                ms_rs = ms_rs + 1
                                ! The number of cases in each modality (independently from the other one) needs to be updated
                                ! to locate the correct placement value in the V01 array.
                                ms_r  = ms_r + 1
                                ms_s  = ms_s + 1
                                ! values are first stored in this array, which is later ranked for stability reasons, it has as
                                ! many elements as the number of cases in common to both modalities.
                                RP_V10(ms_rs)= ( V10(ms_r,r_mod) - U_vec(r_mod) ) * ( V10(ms_s,s_mod) - U_vec(s_mod) )
                         ! the value does not enter in the summation to compute the variance, however, we need to update the
                         ! counter for the V01 array
                         elseif(des_pos(i_case,r_mod) == 1 .and. des_pos(i_case,s_mod) == 0) then
                                ms_r = ms_r + 1
                         elseif(des_pos(i_case,r_mod) == 0 .and. des_pos(i_case,s_mod) == 1) then
                                ms_s = ms_s + 1
                         endif
                  enddo loop_positives
                  ! Rank the cross products, this is a stabilization step added for numerical reasons. It is usually unnecessary
                  ! but better safe than sorry.
                  call qsortd(RP_V01(1:mn_rs), Ind_RP_V01(1:mn_rs), mn_rs)
                  call qsortd(RP_V10(1:ms_rs), Ind_RP_V10(1:ms_rs), ms_rs)
                  ! Sum variance components from the smallest to the largest
                  eps_01 = 0.0_double
                  do j_case = 1, mn_rs
                        eps_01 = eps_01 + RP_V01( Ind_RP_V01(j_case) )
                  enddo
                  eps_01 = eps_01 / (mn_rs - 1)

                  eps_10 = 0.0_double
                  do i_case = 1, ms_rs
                        eps_10 = eps_10 + RP_V10( Ind_RP_V10(i_case))
                  enddo
                  eps_10 = eps_10 / (ms_rs - 1)

                  ! Combine eps10 and eps01 according to eqn. (3) page 1689 of Zhou and Gatsonis (following my notation which has the
                  ! advantage of being understandable.
                  U_vec_cov(r_mod, s_mod) =   eps_01 *  mn_rs / (real(mn_r,kind=double) * mn_s)     &
                                           +  eps_10 *  ms_rs / (real(ms_r,kind=double) * ms_s)
                  ! load the value to the symmetrical location.
                  U_vec_cov(s_mod, r_mod) =  U_vec_cov(r_mod, s_mod)
        enddo s_loop

 enddo r_loop


!---------------------------------------------------------------------------
end subroutine ZhouAndGatsonis
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------


!------------------------------------------------
 subroutine DeLongAndDeLong(mn, ms, num_mod, act_neg ,act_pos, U_vec , U_vec_cov)
!------------------------------------------------
!PURPOSE:  Compute array with the Mann-Whitney form of the Wilcoxon for a set
!          of fully-paired treatments (fixed factors). It also provides the
!          variance covariance matrix for that set.
!          Wilcoxon is known to be equivalent to the Area Under the trapezoid
!          determined by the Emprical ROC points. Also called C-statistic.
!ALGORITHM:uses the function compute_V10_V01(Mn, Ms, act_neg, act_pos, V01, V10) to
!          to compute the placement values for the act-positive and 1- placement
!          for the act-negative cases for each treatment, and the uses the results from
!          the De Long and De Long paper to estimate the variance-covariance matrix.
! SOURCE:  DeLong, DeLong and Clarke-Pearson, Biometrics 44(3),837 (1988) and some pages of notes that I have
!          already shredded.
! NOTE:    It works only for fully paired datasets.
!          !!! changed some integers to real to avoid overflows.
 implicit none

! Declaration of the interface
 integer, intent(IN):: mn, ms ! number of actually-negative and actually-positive
 integer, intent(IN) :: num_mod ! number of treatments
 real(kind=double), dimension (mn,num_mod),intent(IN):: act_neg     ! matrix with the values for the actually-negative cases
                                      ! for each treatment
 real(kind=double), dimension (ms,num_mod),intent(IN):: act_pos     ! matrix with the values for the actually-pos cases
                                      ! for each treatment

 real(kind=double),dimension(num_mod), intent(out):: U_vec    ! Array of Wilcoxon statistics, one per treatment
 real(kind=double),dimension(num_mod,num_mod), intent(out):: U_vec_cov    ! Variance-covariance matrix of the U_vec

 ! Internal variables.
 real(kind=double), dimension(mn,num_mod):: V01 ! Array with the act-negative cases placement values, by treatment
 real(kind=double), dimension(ms,num_mod):: V10 ! Array with the act-positive cases placement values, by treatment

 real(kind=double), dimension(mn):: RP_V01 ! Array with the scalar products of the V01s, ranked -- to stabilize variance calcs
 real(kind=double), dimension(ms):: RP_V10 ! Array with the scalar products of the V01s, ranked -- to stabilize variance calcs
 integer, dimension(mn):: Ind_RP_V01 ! Ranking indices
 integer, dimension(ms):: Ind_RP_V10 ! Ranking indices
 real(kind=double):: VarS, VarN ! the act-negative and act-positive components of the variance.

 integer:: i_mod, j_mod, k ! Loop counters


! Compute the placement values for all the treatments/modalities
 do i_mod = 1, num_mod
          call compute_V10_V01(mn, ms, act_neg(:,i_mod), act_pos(:,i_mod), V01(:,i_mod), V10(:,i_mod))
          U_vec(i_mod) =  sum( V01(:,i_mod) )/mn ! sum(V10)/ms should produce the same answer
  enddo

! the variance covariance matrix according to DeLong and DeLong and Clarke-Pearson.
! The centering is better done before the calculations for stability reasons. This error can
! be shown to be equivalent to some sort of Jackknife.


 do i_mod = 1, num_mod
    do j_mod = i_mod, num_mod
        ! build by negative and by positive cross products of V's
        forall(k= 1: mn)
                RP_V01(k)= ( V01(k,i_mod) - U_vec(i_mod) ) * ( V01(k,j_mod) - U_vec(j_mod) )
        end forall
        forall(k= 1: ms)
                RP_V10(k)= ( V10(k,i_mod) - U_vec(i_mod) ) * ( V10(k,j_mod) - U_vec(j_mod) )
        end forall
        ! Rank the cross products
        call qsortd(RP_V01, Ind_RP_V01, mn)
        call qsortd(RP_V10, Ind_RP_V10, ms)
        ! Sum variance components from the smallest to the largest
        VarN = 0.0_double
        do k = 1, mn
          VarN = VarN + RP_V01( Ind_RP_V01(k))
        enddo
        VarS = 0.0_double
        do k = 1, ms
          VarS = VarS + RP_V10( Ind_RP_V10(k))
        enddo
        ! LP Dec 2008, introduced the real(**) to deal with very large numbers that make the multiplication
        ! go out of range
        U_vec_cov(i_mod,j_mod) = VarN/( real(mn-1,kind=double)*mn) +  VarS/( real(ms-1,kind=double)*ms)


!         U_vec_cov(i_mod, j_mod) = &
!                    dot_product( (V01(:,i_mod) - U_vec(i_mod)) , (V01(:,j_mod) - U_vec(j_mod)) ) / ( (mn-1)*ms ) + &
!                    dot_product( (V10(:,i_mod) - U_vec(i_mod)) , (V10(:,j_mod) - U_vec(j_mod)) ) / ( (ms-1)*mn )

         if(i_mod /= j_mod) U_vec_cov(j_mod, i_mod) =  U_vec_cov(i_mod, j_mod)

    enddo
 enddo


!------------------------------------------------
 end subroutine DeLongAndDeLong
!------------------------------------------------


!------------------------------------------------------------------------------------
 subroutine compute_V10_V01(mn, ms, Y, X, V01, V10)
!-------------------------------------------------------------------------------------------
! PURPOSE  :Compute the functions V10 and V01, which are the placement value and the 1 - placement value
!           of each positive and negative case in relation to the other distribution, respectively.
! ALGORITHM:Uses quick sort algorithm to create two ranked sequences, one for X (the positives) and one for
!           Y (the negatives) -- Notation follows approximately the original article. Then the two sequences are
!           scanned (together) from the smallest value to the largest value.
! NOTE:     Positivity is assumed to be for larger values, one can submit -X and -Y if this is not the case.
! NOTE1:    Here we follow the The Long and the Long notation and give 10 to the positives and 01 to the negatives.
!           I usually prefer to do the opposite so be mindful when looking at other functions.
! SOURCE:   DeLong, DeLong and Clarke-Pearson, Biometrics 44(3),837 (1988) and some pages of notes that I have
!           already shredded.
! AUTHOR:   LP U CHICAGO 2008 CE
!-------------------------------------------------------------------------------------------

! Interface definition
integer, intent(IN):: mn, ms ! Number of act-negative and of act-positive cases, respectively
real(kind=double), intent(IN), dimension(mn):: Y ! Array with the act-negative cases
real(kind=double), intent(IN), dimension(ms):: X ! Array with the act-positive cases

real(kind=double), intent(OUT), dimension(mn):: V01 ! Array with the act-negative cases 1 - placement values
real(kind=double), intent(OUT), dimension(ms):: V10 ! Array with the act-positive cases placement values

! Internal variables
integer, dimension(mn):: index_Y ! Array with the sorting indexes of the act-neg cases, i.e., the first element of the
                                 ! array contains the location of the smallest value in the array X.
integer, dimension(ms):: index_X ! Array with the sorting indexes of the act-pos cases
integer:: i ! loop counter
integer:: i_neg, i_pos ! value of the current count for the act-negative and act-positive cases arrays
integer:: index_i_neg, index_i_pos ! location of the current counter index for the negative and positive cases
integer:: i_neg_old, i_pos_old ! counters for the negative and positive cases used to locate the beginning of a tie between act-pos and
                               ! act-neg cases (see body of procedure to see why we do this).
logical:: previous_was_tie ! Logical variable that defines whether the previous run was a tie,
real(kind=double):: incr_pos, incr_neg ! increments for positive and negatives after ties are not the same and the some other uses I don't
                                       ! real want to talk about right now.

! Sort act-negative and act-positive cases and build the sorting index arrays.
call qsortd(Y, index_Y, mn)
call qsortd(X, index_X, ms)

! Start from the smallest values. We will count how many act-neg are smaller or equal then the current act-pos and
! by reciprocity how many act-pos are larger or equal than the current act-neg.
! We  work the cases from smallest to largest using the indexes i_neg and i_pos to know which rank are we at and the
! convert to the  indices index_i_neg and index_i_pos tell us where the corresponding values are and which V01 and V10
! functions to update. In this way the resulting functions will have the same ordering as the input arrays.

! Initialize variables used in the loop to the first two cases.
i_neg = 1
i_pos = 1
i_neg_old = 1
i_pos_old = 1
index_i_neg = index_Y(i_neg)
index_i_pos = index_X(i_pos)
V01(index_i_neg) = real(ms, kind=double) !It is the smallest, we set to the max and take out all the act-pos cases that are smaller than this one
V10(index_i_pos) = 0.0_double ! we will add all the act-neg cases smaller than this one, which is the smallest.
previous_was_tie = .false. ! There was no previous, the logic of the algorithm requires stating false, complain to /dev/null if you don't like it

! Endless loop, exit is produced by reaching the end of both arrays

build_Vs: do
! First case: the act-positive is larger than the act-negative
      if( X(index_i_pos) > Y(index_i_neg)) then
             ! Notice that the algorithm exist from a tie sequence only by increasing the act-pos, if this is not possible, it exits.
             if(previous_was_tie) then ! it means the previous comparison was a tie between act-pos and act-neg
                 ! Add all the act-neg cases that were in the previous tie stretch to the current act-pos V10 -- i_neg was the last of the stretch
                 incr_pos = (i_neg - i_neg_old + 1)*0.5_double
                 ! Subtract all the act-pos cases that were in the tie stretch to the next act-neg, which will be larger then the tie because the
                 ! tie sequence runs through all the elements of each tie stretch. i_pos was updated previously, so the current i_pos is outside
                 ! the tie-stretch
                 incr_neg = (i_pos - i_pos_old)*0.5_double
                 previous_was_tie = .false. ! Current stretch is not a tie, reset the flag.
             else
                 incr_pos = 1.0_double
                 incr_neg = 0.0_double ! We are not testing the negative, but only updating if coming out of a tie-stretch
             endif

             ! Check if it is the last act-negative, otherwise the remaining act-positives will all be larger than all the act-negatives
             if(i_neg == mn) then
                 forall(i=i_pos:ms) V10(index_X(i)) = real(mn,kind=double)
                 exit build_Vs ! the loop is finished, both act-neg and act-pos sequences have been swept.
             else
                 ! Update the placement value for the act-pos
                 V10(index_i_pos) =  V10(index_i_pos) + incr_pos
                 ! increase the counter for the act-negative and update its value.
                 i_neg = i_neg + 1  ! it can't be a tie from the previous stretch because tie-stretches are loaded completely.
                 V01(index_Y(i_neg))= V01(index_i_neg) - incr_neg ! We need to subtract all the cases that were in the tie stretch
                                                                  ! --if it was a tie-stretch
                 index_i_neg = index_Y(i_neg) ! update the "current negative index"
             endif
! Start of a tie-stretch, act-pos and act-neg cases have the same value
      elseif( X(index_i_pos) == Y(index_i_neg)) then
             ! If they are equal we will keep incrementing until all the cases with values tied to these ones have been located.
             ! First load the current values as starting the starting point of the tie-stretch.
             i_pos_old = i_pos
             i_pos = i_pos + 1
             ! Find the end of the tie-stretch for both the act-negative and the act-positive
             find_equals_pos: do
                 if( i_pos > ms ) then
                        i_pos = i_pos - 1 ! go back to the previous which has to be the last one in the tie-stretch
                        exit find_equals_pos
                 elseif( X(index_X(i_pos)) > Y(index_Y(i_neg))) then
                        i_pos = i_pos - 1 ! go back to the previous which has to be the last one in the tie-stretch
                        exit find_equals_pos
                 else ! move one more up
                       i_pos = i_pos + 1
                 endif
             enddo find_equals_pos
             i_neg_old = i_neg
             i_neg = i_neg + 1
             find_equals_neg: do
                 if( i_neg > mn ) then
                       i_neg = i_neg - 1 ! go back to the previous which has to be the last one in the tie-stretch
                       exit find_equals_neg
                 elseif( X(index_X(i_pos)) < Y(index_Y(i_neg))) then
                       i_neg = i_neg - 1 ! go back to the previous which has to be the last one in the tie-stretch
                       exit find_equals_neg
                  else
                       i_neg = i_neg + 1
                  endif
             enddo find_equals_neg

             ! Update the placement values, first set the stating value, then copy it to the rest of the act-pos tie-stretch
             V10(index_X(i_pos_old)) =  V10(index_X(i_pos_old)) + 0.5_double*(i_neg + 1 - i_neg_old) ! For each act-neg add  1/2 case more
             do i = i_pos_old+1, i_pos
                  V10(index_X(i)) =  V10(index_X(i_pos_old)) ! They are ties, so they have the same value
             enddo

             ! Update the placement values, first set the stating value, then copy it to the rest of the act-neg tie-stretch
             V01(index_Y(i_neg_old)) =  V01(index_Y(i_neg_old)) - 0.5_double*(i_pos + 1 - i_pos_old) ! There is 1/2 case more per act-pos
             do i = i_neg_old+1, i_neg
                  V01(index_Y(i)) =  V01(index_Y(i_neg_old))  ! They are ties, so they have the same value
             enddo
             ! All the ties between act-neg and act-pos are exhausted, we move to a larger positive, if possible
             if(i_pos < ms) then
                  ! Move to the next act-positive, update the current act-positive locator.
                  i_pos = i_pos + 1
                  V10(index_X(i_pos)) =  V10(index_X(i_pos_old)) ! set the current act-positive placement value to the previous
                  index_i_pos = index_X(i_pos) ! update the index
                  !Locate the top negative in the tie-stretch
                  index_i_neg = index_Y(i_neg) !
             else ! We swept all the positives and all the negatives remaining negatives are larger thus inv-placement is zero
                  forall(i=i_neg+1:mn) V01(index_Y(i)) = 0.0_double
                  ! All cases are swept
                  exit build_Vs
             endif
             previous_was_tie = .true.! it is the end of this update, so for the next cycle the previous one was a tie.
! The current act-negative is larger then the current act-positive
      else
             ! Note that the previous step can't be a tie because ties either exit by updating the act-positives or the loop was concluded
             V01(index_i_neg) =  V01(index_i_neg) - 1.0_double  !take current act-pos case out of the placement value because <current act-neg
             if(i_pos == ms) then ! It was the last act-pos, all the following act-neg > all act-pos
                  forall(i=i_neg+1:mn) V01(index_Y(i)) = 0.0_double
                  exit build_Vs ! the loop is finished
             else ! increase the positives by one
                  i_pos = i_pos + 1
                  V10(index_X(i_pos)) =  V10(index_i_pos) ! set the next positive placement value to the current
                  index_i_pos = index_X(i_pos) ! update the index
             endif
      endif
end do build_Vs

 V10 = V10 / real(mn, kind=double)
 V01 = V01 / real(ms, kind=double)

 end subroutine compute_V10_V01
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------


 end module roc_nonparametric


