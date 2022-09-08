! This module contains functions which are specific to PROPER ROC analysis

! PROGRAMMING NOTES:
!    1) The code is ANSI Fortran 95, it will not compile unders
!        strict Fortran 90 or previous versions of the fortran language
!    2) Deprecated forms of fortran 77 are not used
!    3) The module has the "private" status, so only the methods which are
!       declared as public can be used by external programs
!    4) In ROC, the ROC curve values (FPF and  TPF) increase when the
!       operating points TEST RESULTS VALUES are
!       decresed. We decided to keep the indexing so that TEST RESULT VALUES
!       are in ascending order, which  makes FPF and TPF and related quantities
!        decrease with it.

! IT ASSUMES THAT THE CATEGORIES ARE ORDERED FROM THE LOWEST PROBABILITY OF BEING
! POSITIVE (CAT 1) TO THE HIGHEST (CAT NUM_CAT)
! NOTES:
! 1-18-06 - LP - Introduced the partial area computation
! 2-2-06  -LP - introduces the partial area for single point fits (one_point_fit_pa)
! 10-09-08 LP  - changed the names of functions to become something like d_fpf_d_d_a_PBM, i.e., leaving the proproc at the end as a
!                qualifier, the definition at the beginning (d for derivative, var for variance) and ordering the rest from general
!                tro specific. The word "proper" has been replaced with "proproc" to make a clear distinction from other proper models.
! April 30 2009 - LP -> changed ierror = 2 to ierror = bad_input. Added also checkes on input values for variances for variance estimation
!                       routines
 module proproc_functions
! NOTE that only the constants for actually negative and actually positive
! cases are loaded from the module array_dimensions. The only is used to that
! a strick control of data and methods is enforced
use data_types, only: double, operator (.speq.), operator (.spne.), check_number
use error_flags
use computation_constants, only:pi, p_q_closeto_zero
use proproc_computation_constants
use array_dimensions, only: act_neg, act_pos ! where array related quantities are
                                ! variable to be eliminated that has to do with the number of
                                ! parameters (XP is such an awful programmer that nothing makes
                                ! real sense in his junk)
use io

! GENERAL SETTING FOR THE WHOLE MODULE
implicit none

private ! All the variables and functions that are not set as public
        ! are not acessible

! CLEAN THIS ONE UP ONCE PROPROC IS FINISHED
! PUBLIC PROCEDURES / VALUES
public double ! Note that we are passing this value, which is inherited from data_types
public line_length ! Note that we are passing this value, which is inherited from io
! PUBLIC PROCEDURES
public pbmroc_mle ! main proper model procedure
public one_point_fit ! compute fit when only one pt is available (undetermined case)
public fpf_find_vc ! find vc, knowing d_a, c, fpf. Previously called one_pt_fpf_find_vc
public tpf_PBM ! find tpf, knowing the curve parameters and the cutoff value
public tpf_find_vc ! find vc, knowing d_a, c, tpf. Previously called  one_pt_tpf_find_vc
public fpf_find_tpf_PBM ! find TPF, knowing d_a, c, fpf
public var_fpf_find_tpf_PBM ! variance of TPF, knowing d_a, c, fpf.
public tpf_find_fpf_PBM ! find fpf, knowing d_a, c, tpf.
public fpf_PBM ! find fpf, knowing the curve parameters and the cutoff value
public density_vc_PBM  ! the density of vc, for positive or negative cases
public d_tpf_d_d_a_PBM,d_fpf_d_d_a_PBM,d_tpf_d_c_PBM,d_fpf_d_c_PBM ! derivatives of tpf and fpf
public partialauc_PBM, var_partialauc_PBM, one_point_fit_pa
public auc_PBM, var_auc_PBM ! compute the area under the curve
! Export the meaning of the different values that the fit error flag can have
public bad_input          ! categorical input has unacceptable values
public fit_OK             ! fit was successful
public fit_OK_pseudo      ! fit was successful, variances computed using pseudoinverse
public fit_failed         ! fit failed in MLE estimation
public fit_undetermined   ! Note enough data: there is only one point, and it is inside the ROC space
public fit_perfect        ! positives and negatives are perfectly separated, it is more of a warning
public fit_fail_init_est  ! initial estimates did not converge
public fit_fail_variances! estimates of variances did not converge
public fit_fail_var_auc_small! estimates of variances did not converge because area is too small (see below)

!public  hessian_matrixAUC: compute the hessisian in the gradient AUC direction, currently not used


contains

!---------------------------------------------------------------------------------------------------
subroutine one_point_fit_pa(mn, ms, num_categ, catn_in, cats_in, idebug, &
                         fraction_1, fraction_2, fpf_flag, partial_auc, ierror, err_msg)
! PURPOSE: estimate the partial area of a curve for which only one point is available. Clearly this estimation
!          carries a large error. We currently don't provide the variance, because
!          it is likely to be unreliable. In the end the only reasonable use of this procedure is
!          for jackknifing in ANOVA, for the few very rare datasets where this happens. In one curve fits
!          It is very irresponsible to make use of this fitting routine. Note that it could return any roc index
!          using this principle, for example it could mediate partial areas or TPF for fpf.
!          A possible alternative for the variance is  [max(partial_auc) - mean(partial_auc)]^2 + rule_of_thumb_var, so adding
!          the variance caused by the infinite number of possible curves that fit the point exactly to the
!          "natural" variances for a curve with this area. It is difficult to concoct a simulation to verify this,
!          since the whole calculation is doomed to be extremely unreliable because it takes very bad data
!          collection schemes to obtain points like these.
!ALGORITHM: first find the fits with d_a = 0. These are the largest and smallest values of c. Then it makes a
!           grid in the c space (it uses c because curves with different c are different, which is not necessarily
!           the case for d_a). It computes the exact fit for all of this curves. It finds the largest and smallest
!           partial_auc and uses the mean as its MLE estimate. This basically assumes that all areas are equally probably.
!           It should be noted that the likelihood for these one point exact fits is always the same, so the
!           density is in fact the same for each curve.

! Note:    It assumes that in fact there are 2 categories, without checking it
! Note2:   It doesn't return the curve parameters because we can't really  make claims about the curve shape
!          when we are fitting a single point
! Note3:   the algorithm is not meant to be fast, only to be easy and reliable. Speed is totally irrelevant

use debugging ! module containing the procedures to log the debugging information on the run
use problem_data ! where the ROC basic categorical information is stored.
use statistic_functions, ONLY : compute_zdev_plus
implicit none

! Subroutine arguments
integer, intent(in):: mn ! number of actually negative cases
integer, intent(in):: ms ! number of actually positive cases
integer, intent(in) :: num_categ ! Number of categories as created by catgrz
integer, dimension(num_categ), intent(in):: catn_in, cats_in ! arrays containing categorical data
integer, intent(in) :: idebug               ! 0 = no debug; 1  = debug
  real(kind=double), intent(in):: fraction_1, fraction_2 ! These are called fractions
       ! because the can be FPF or TPF depending upon which area are we computing
 integer, intent(in):: fpf_flag ! If it is 1, it means that the area computed will
       ! be between FPF_1 and FPF_2, otherwise it means that it will be between TPF_1 and
       ! TPF_2
  real(kind=double), intent(out):: partial_auc
integer, intent(out)                                   :: ierror ! Error flag about the MLE fit
! Note that the error values are set in this routine, or initialize_d_a_c so if their numbers have
! to be changed, they have to be changed here, the rest of the subroutines use their own numbering
! specific per routine. Look above for the different meaning. Not only 0 is successful fit
! ierror = 2, failed initial fits for d_a = 0
character(len = line_length), intent(out):: err_msg    ! description of the error occurred. If the fit worked, it
                                                       ! contains information about how were the variances computed

! Internal variables
real(kind = double):: d_a_par, c_par ! binormal proper model parameters
integer:: i, j ! loop counters
!character(len = line_length):: msg    ! character string used to display state information
real(kind = double) :: p, q ! location of the single point in the ROC plot, FPF, TPF
real(kind = double) :: vc ! cutoff value, used to solve the proper equations
real(kind = double) :: tpf, one_minus! location of the single point in the ROC plot
real(kind = double) :: dev_p, dev_q ! deviates, used if c =~ 0, since we can solve it analytically

integer, parameter:: n_c_par = 100 ! number of values of c to sample to find the maximum and minumum partial_auc
real(kind = double) :: delta_c ! step in the sampling of c
real(kind=double), parameter:: tol = 1.0e-10_double ! tolerance in the computation of tpf and fpf
real(kind=double):: d_a_small, d_a_large ! values of d_a that bracket the current root
integer, parameter:: max_iter = 1000 ! maximum number of iterations used to find the root
real(kind=double):: max_partial_auc, min_partial_auc ! Maximum and minimum, used to compute the mean
real(kind=double), parameter:: c_small = 1.0e-6_double ! values below this one are considered 0, note that we are not
! being as accurate as we could be (the normal proproc is much more accurate that this),
! the reason is that this kind of dataset carries such a large error that it is useless to
! worry too much for accuracy, which costs time and programming care. In fact we are being too careful for the
! needs of this algorithm.

! Debug is used to pass the debugging flag around optimizers that have no business knowing its value
if(idebug == 1) then
 debug = .true.
else
 debug = .false.
endif

! Prepare the variables that will be stored in the moduli from file "proproc_other_modules.f90"
! Since those are in a module while these are explicit, we seem to need to double the
! names at this level. The moduli are used to work around the various optimizers (initial estimates
! of cutoffs and full MLE). NOTE: all the calls that don't have to "pass-through" an optimizer are
! made explicit, even if this might seem confusing, I prefer explicit calls when I can do them.
num_cat = num_categ
allocate(catn(num_cat), cats(num_cat))
num_normal_cases   = mn
num_abnormal_cases = ms
catn(1:num_cat) = catn_in(1:num_cat)
cats(1:num_cat) = cats_in(1:num_cat)
pass_idebug  = idebug

! Compute the operating point, the single operating point
p = real(catn(2)) / num_normal_cases
q = real(cats(2)) / num_abnormal_cases


! First do the negative c part
call one_pt_d_a_zero_fit(p/2.0_double, q/2.0_double, d_a_par, c_par, vc, ierror, err_msg)

if(ierror /=0 ) then ! change to the appropriate error flag and return
 ierror = 2
 deallocate(catn , cats)
 return
endif

! Compute the value of partial_auc for this parameters
call partialauc_PBM(d_a_par, c_par, fraction_1, fraction_2, fpf_flag, partial_auc, ierror)
! We don't care for this error for the time being
ierror = 0
err_msg = ""
! initialize the variables with the largest and smallest partial_auc
max_partial_auc = partial_auc
min_partial_auc = partial_auc

! compute the step in the c grid
delta_c = c_par / (n_c_par - 1)

do i = 1, n_c_par - 1 ! look for the largest and smallest partial_auc on the grid
    c_par = c_par - delta_c ! increment the parameter c, to sample the grid
   ! set some rought bracketing for d_a
    d_a_small = 1.0e-10_double
    d_a_large = 4.0_double
    ! Make sure that small roundoffs didn't change the sign of this value of c
    if(c_par > c_small)  c_par = - c_par
!   if c is too close to zero, we better assume that it is in fact zero
!   and solve the equation analytically, otherwise, solve it iteratively
!   Note again that in this situation there is no need to be very accurate since
!   there is a very large uncertainty associates with these estimates. If very
!   accurate data exist which have this form, clearly ROC analysis and proper models
!   are not appropriate (exceptions shall be dealt with by their creators)
    if( abs (c_par) <= c_small) then
       call compute_zdev_plus(p, dev_p, ierror)
       call compute_zdev_plus(q, dev_q, ierror)
       d_a_par = dev_q - dev_p
       vc = d_a_par/2.0_double - dev_q
    else
       do j = 1, max_iter ! iterate to find the values of d_a and vc that solves FPF = p AND TPF = q
          ! Please note that this is sort of a generalization of a bisection algorithm
          call fpf_find_vc(d_a_par, c_par,p, tol, vc) ! find vc for the current values of d_a and c using FPF
          call tpf_PBM(d_a_par, c_par, vc, tpf, one_minus) ! compute the TPF using this last value of vc
          if( abs(tpf - q) < tol) then ! Check whether also the TPF is correct
              exit   ! calculation converged, we have found d_a
          elseif ( tpf < q)  then ! FPF is too small, it means that the d_a root is larger than the current value
             d_a_small = d_a_par
             d_a_par = (d_a_par + d_a_large)/2.0_double
          else! FPF is too large, it means that the d_a root is smallerer than the current value
             d_a_large = d_a_par
             d_a_par = (d_a_par + d_a_small)/2.0_double
         endif
      enddo
    endif
    ! determine the maximum and minimum values of partial_auc for the different exact fits
    call partialauc_PBM(d_a_par, c_par, fraction_1, fraction_2, fpf_flag, partial_auc, ierror)
    ierror = 0
    ! Check if the current value is a the largest or the smallest so far
    max_partial_auc = max(max_partial_auc, partial_auc)
    min_partial_auc = min(min_partial_auc, partial_auc)
enddo

! then do the positive c part (one could do the 2 in one loop).
! Note that the fit is computed in the same way as for c negative, apart that instead of p and q one has p+1 and q+1
call one_pt_d_a_zero_fit( (p + 1.0_double) / 2.0_double, (q + 1.0_double) / 2.0_double, d_a_par, c_par, vc, ierror, err_msg)
if(ierror /=0 ) then ! change to the appropriate error flag and return
 ierror = 2
 deallocate(catn , cats)
 return
endif

! determine the maximum and minimum values of partial_auc for the different exact fits
call partialauc_PBM(d_a_par, c_par, fraction_1, fraction_2, fpf_flag, partial_auc, ierror)
ierror = 0
! initialize  the maximum and minimum values of partial_auc
max_partial_auc = max(max_partial_auc, partial_auc)
min_partial_auc = min(min_partial_auc, partial_auc)

! compute the step in the c grid
delta_c = c_par / (n_c_par - 1)

! See comments for the negative c part is some steps aren't clear, the solution is simmetrical
do i = 1, n_c_par - 1
    c_par = c_par - delta_c ! increment the parameter c, to sample the grid
   ! set some rought bracketing for d_a
    d_a_small = 1.0e-10_double
    d_a_large = 4.0_double
   if(c_par < c_small)  c_par = - c_par
!   if c is too close to zero, we better assume that it is in fact zero
!   and solve the equation analytically, otherwise, solve it iteratively
    if( abs (c_par) <= c_small) then
       call compute_zdev_plus(p, dev_p, ierror)
       call compute_zdev_plus(q, dev_q, ierror)
       d_a_par = dev_q - dev_p
       vc = d_a_par/2.0_double - dev_q
    else
        do j = 1, max_iter
            call fpf_find_vc(d_a_par, c_par,p, tol, vc)
            call tpf_PBM(d_a_par, c_par, vc, tpf, one_minus)
            if( abs(tpf - q) < tol) then
                exit
            elseif ( tpf < q)  then
                d_a_small = d_a_par
                d_a_par = (d_a_par + d_a_large)/2.0_double
            else
                d_a_large = d_a_par
                d_a_par = (d_a_par + d_a_small)/2.0_double
            endif
        enddo
    endif
    ! determine the maximum and minimum values of partial_auc for the different exact fits
    call partialauc_PBM(d_a_par, c_par, fraction_1, fraction_2, fpf_flag, partial_auc, ierror)
    ierror = 0
    max_partial_auc = max(max_partial_auc, partial_auc)
    min_partial_auc = min(min_partial_auc, partial_auc)
enddo

partial_auc = ( min_partial_auc + max_partial_auc) / 2.0_double
deallocate(catn , cats)


!---------------------------------------------------------------------------------------------------
end subroutine one_point_fit_pa
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine one_point_fit(mn, ms, num_categ, catn_in, cats_in, idebug, &
                         auc , var_auc, ierror, err_msg)
! PURPOSE: estimate the area of a curve for which only one point is available. Clearly this estimation
!          carries a large error. We currently use the rule of thumb to determine the variance, but
!          it is likely to be an underestimate. In the end the only reasonable use of this procedure is
!          for jackknifing in ANOVA, for the few very rare datasets where this happens. In one curve fits
!          It is very irresponsible to make use of this fitting routine. Note that it could return any roc index
!          using this principle, for example it could mediate partial areas or TPF for FPF.
!          A possible alternative for the variance is  [max(auc) - mean(auc)]^2 + rule_of_thumb_var, so adding
!          the variance caused by the infinite number of possible curves that fit the point exactly to the
!          "natural" variances for a curve with this area. It is difficult to concoct a simulation to verify this,
!          since the whole calculation is doomed to be extremely unreliable because it takes very bad data
!          collection schemes to obtain points like these.
!ALGORITHM: first find the fits with d_a = 0. These are the largest and smallest values of c. Then it makes a
!           grid in the c space (it uses c because curves with different c are different, which is not necessarily
!           the case for d_a). It computes the exact fit for all of this curves. It finds the largest and smallest
!           auc and uses the mean as its MLE estimate. This basically assumes that all areas are equally probably.
!           It should be noted that the likelihood for these one point exact fits is always the same, so the
!           density is in fact the same for each curve.

! Note:    It assumes that in fact there are 2 categories, without checking it
! Note2:   It doesn't return the curve parameters because we can't really  make claims about the curve shape
!          when we are fitting a single point
! Note3:   the algorithm is not meant to be fast, only to be easy and reliable. Speed is totally irrelevant

use debugging ! module containing the procedures to log the debugging information on the run
use problem_data ! where the ROC basic categorical information is stored.
use statistic_functions, ONLY : compute_zdev_plus
implicit none

! Subroutine arguments
integer, intent(in):: mn ! number of actually negative cases
integer, intent(in):: ms ! number of actually positive cases
integer, intent(in) :: num_categ ! Number of categories as created by catgrz
integer, dimension(num_categ), intent(in):: catn_in, cats_in ! arrays containing categorical data
integer, intent(in) :: idebug               ! 0 = no debug; 1  = debug
real(kind = double),intent(out)                        :: auc ! AUC, area under the curve
real(kind = double),intent(out)                        :: var_auc ! var of  AUC, area under the curve

integer, intent(out)                                   :: ierror ! Error flag about the MLE fit
! Note that the error values are set in this routine, or initialize_d_a_c so if their numbers have
! to be changed, they have to be changed here, the rest of the subroutines use their own numbering
! specific per routine. Look above for the different meaning. Not only 0 is successful fit
! ierror = 2, failed initial fits for d_a = 0
character(len = line_length), intent(out):: err_msg    ! description of the error occurred. If the fit worked, it
                                                       ! contains information about how were the variances computed

! Internal variables
real(kind = double):: d_a_par, c_par ! binormal proper model parameters
integer:: i, j ! loop counters
real(kind = double) :: p, q ! location of the single point in the ROC plot, FPF, TPF
real(kind = double) :: vc ! cutoff value, used to solve the proper equations
real(kind = double) :: tpf, one_minus! location of the single point in the ROC plot
real(kind = double) :: dev_p, dev_q ! deviates, used if c =~ 0, since we can solve it analytically

integer, parameter:: n_c_par = 100 ! number of values of c to sample to find the maximum and minumum auc
real(kind = double) :: delta_c ! step in the sampling of c
real(kind=double), parameter:: tol = 1.0e-10_double ! tolerance in the computation of tpf and fpf
real(kind=double):: d_a_small, d_a_large ! values of d_a that bracket the current root
integer, parameter:: max_iter = 1000 ! maximum number of iterations used to find the root
real(kind=double):: max_auc, min_auc ! Maximum and minimum, used to compute the mean
real(kind=double), parameter:: c_small = 1.0e-6_double ! values below this one are considered 0, note that we are not
! being as accurate as we could be (the normal proproc is much more accurate that this),
! the reason is that this kind of dataset carries such a large error that it is useless to
! worry too much for accuracy, which costs time and programming care. In fact we are being too careful for the
! needs of this algorithm.

! Debug is used to pass the debugging flag around optimizers that have no business knowing its value
if(idebug == 1) then
 debug = .true.
else
 debug = .false.
endif

! Prepare the variables that will be stored in the moduli from file "proproc_other_modules.f90"
! Since those are in a module while these are explicit, we seem to need to double the
! names at this level. The moduli are used to work around the various optimizers (initial estimates
! of cutoffs and full MLE). NOTE: all the calls that don't have to "pass-through" an optimizer are
! made explicit, even if this might seem confusing, I prefer explicit calls when I can do them.
num_cat = num_categ
allocate(catn(num_cat), cats(num_cat))
num_normal_cases   = mn
num_abnormal_cases = ms
catn(1:num_cat) = catn_in(1:num_cat)
cats(1:num_cat) = cats_in(1:num_cat)
pass_idebug  = idebug

! Compute the operating point, the single operating point
p = real(catn(2)) / num_normal_cases
q = real(cats(2)) / num_abnormal_cases

! Check whether the single point is on the 45 degrees (no better than chance) line
! in that case, there is only one possible fit :-) (and the data kind of .....)
! Please notice that if the point is below the 45 degrees line we can't fit any
! proper curve to it, so we assume it comes from the chance line.
! We used this form of comparison because even if p is larger than q, it isn't numerically meaningful
! Please note that if you have very accurate readings and only one point and even more close to the
! -45 degree line, it means that we have a very serious data collection problem or the problem is
! not an ROC problem, but a categorical analysis problem.
if(  q < p*(1.0_double + 1.0e-10_double) ) then
    auc = .5_double
    deallocate(catn , cats)
    var_auc =  approx_var_auc(auc, num_normal_cases, num_abnormal_cases)
    return
endif


! First do the negative c part
call one_pt_d_a_zero_fit(p/2.0_double, q/2.0_double, d_a_par, c_par, vc, ierror, err_msg)

if(ierror /=0 ) then ! change to the appropriate error flag and return
 ierror = 2
 deallocate(catn , cats)
 return
endif

! Compute the value of AUC for this parameters
call auc_PBM(d_a_par, c_par, auc, ierror)
! We don't care for this error
ierror = 0
! initialize the variables with the largest and smallest auc
max_auc = auc
min_auc = auc

! compute the step in the c grid
delta_c = c_par / (n_c_par - 1)

do i = 1, n_c_par - 1 ! look for the largest and smallest AUC on the grid
    c_par = c_par - delta_c ! increment the parameter c, to sample the grid
   ! set some rought bracketing for d_a
    d_a_small = 1.0e-10_double
    d_a_large = 4.0_double
    ! Make sure that small roundoffs didn't change the sign of this value of c
    if(c_par > c_small)  c_par = - c_par
!   if c is too close to zero, we better assume that it is in fact zero
!   and solve the equation analytically, otherwise, solve it iteratively
!   Note again that in this situation there is no need to be very accurate since
!   there is a very large uncertainty associates with these estimates. If very
!   accurate data exist which have this form, clearly ROC analysis and proper models
!   are not appropriate (exceptions shall be dealt with by their creators)
    if( abs (c_par) <= c_small) then
       call compute_zdev_plus(p, dev_p, ierror)
       call compute_zdev_plus(q, dev_q, ierror)
       d_a_par = dev_q - dev_p
       vc = d_a_par/2.0_double - dev_q
    else
       do j = 1, max_iter ! iterate to find the values of d_a and vc that solves FPF = p AND TPF = q
          ! Please note that this is sort of a generalization of a bisection algorithm
          call fpf_find_vc(d_a_par, c_par,p, tol, vc) ! find vc for the current values of d_a and c using FPF
          call tpf_PBM(d_a_par, c_par, vc, tpf, one_minus) ! compute the TPF using this last value of vc
          if( abs(tpf - q) < tol) then ! Check whether also the TPF is correct
              exit   ! calculation converged, we have found d_a
          elseif ( tpf < q)  then ! FPF is too small, it means that the d_a root is larger than the current value
             d_a_small = d_a_par
             d_a_par = (d_a_par + d_a_large)/2.0_double
          else! FPF is too large, it means that the d_a root is smallerer than the current value
             d_a_large = d_a_par
             d_a_par = (d_a_par + d_a_small)/2.0_double
         endif
      enddo
    endif
    ! determine the maximum and minimum values of auc for the different exact fits
    call auc_PBM(d_a_par, c_par, auc, ierror)
    ierror = 0
    ! Check if the current value is a the largest or the smallest so far
    max_auc = max(max_auc, auc)
    min_auc = min(min_auc, auc)
enddo

! then do the positive c part (one could do the 2 in one loop).
! Note that the fit is computed in the same way as for c negative, apart that instead of p and q one has p+1 and q+1
call one_pt_d_a_zero_fit( (p + 1.0_double) / 2.0_double, (q + 1.0_double) / 2.0_double, d_a_par, c_par, vc, ierror, err_msg)
if(ierror /=0 ) then ! change to the appropriate error flag and return
 ierror = 2
 deallocate(catn , cats)
 return
endif

! determine the maximum and minimum values of auc for the different exact fits
call auc_PBM(d_a_par, c_par, auc, ierror)
ierror = 0
! initialize  the maximum and minimum values of auc
max_auc = max(max_auc, auc)
min_auc = min(min_auc, auc)

! compute the step in the c grid
delta_c = c_par / (n_c_par - 1)

! See comments for the negative c part is some steps aren't clear, the solution is simmetrical
do i = 1, n_c_par - 1
    c_par = c_par - delta_c ! increment the parameter c, to sample the grid
   ! set some rought bracketing for d_a
    d_a_small = 1.0e-10_double
    d_a_large = 4.0_double
   if(c_par < c_small)  c_par = - c_par
!   if c is too close to zero, we better assume that it is in fact zero
!   and solve the equation analytically, otherwise, solve it iteratively
    if( abs (c_par) <= c_small) then
       call compute_zdev_plus(p, dev_p, ierror)
       call compute_zdev_plus(q, dev_q, ierror)
       d_a_par = dev_q - dev_p
       vc = d_a_par/2.0_double - dev_q
    else
        do j = 1, max_iter
            call fpf_find_vc(d_a_par, c_par,p, tol, vc)
            call tpf_PBM(d_a_par, c_par, vc, tpf, one_minus)
            if( abs(tpf - q) < tol) then
                exit
            elseif ( tpf < q)  then
                d_a_small = d_a_par
                d_a_par = (d_a_par + d_a_large)/2.0_double
            else
                d_a_large = d_a_par
                d_a_par = (d_a_par + d_a_small)/2.0_double
            endif
        enddo
    endif
    ! determine the maximum and minimum values of auc for the different exact fits
    call auc_PBM(d_a_par, c_par, auc, ierror)
    ierror = 0
    max_auc = max(max_auc, auc)
    min_auc = min(min_auc, auc)
enddo

auc = ( min_auc + max_auc) / 2.0_double
var_auc =  approx_var_auc(auc, num_normal_cases, num_abnormal_cases)
deallocate(catn , cats)


!---------------------------------------------------------------------------------------------------
end subroutine one_point_fit
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine fpf_find_tpf_PBM(d_a_par, c_par, fpf , tpf, ierror)
!---------------------------------------------------------------------------------------------------
! PURPOSE: given a value for the FPF, find the corresponding TPF, given the parameters d_a and c
! Warning the error currently does not return failed. Not clear yet how would it.

real(kind=double), intent(in):: d_a_par, c_par ! parameters of the current fit, used as input
real(kind=double), intent(in):: fpf ! value of the fpf for the single point available for the fit
real(kind=double), intent(OUT):: tpf ! value of TPF for that fpf
integer, intent(OUT)::  ierror ! error flag. If ierror = bad_input, the value of FPF is out of bounds

real(kind=double), parameter:: tol = 10e-6_double ! Accepted tolerance in the estimation of FPF
real(kind=double):: vc ! value of vc corresponding to FPF
real(kind=double):: junk

! Check for impossible values of the input parameters
if(fpf > 1.0_double .or. fpf < .0_double .or. d_a_par < .0_double .or. abs(c_par) > 1.0_double) then
  ierror = bad_input
  return
elseif(fpf > 0.999_double) then
  tpf = 1
  return
elseif(fpf < .000001_double) then
  tpf = 0.000001_double
else
  ! Initialize vc, using a proxy for the infinite bound for the proper binormal model
  vc = - sign(10.0_double, c_par)
  ! find the cutoff value that corresponds to FPF
  call fpf_find_vc(d_a_par, c_par, fpf, tol,  vc)
  call tpf_PBM(d_a_par, c_par, vc ,tpf, junk)
endif

!---------------------------------------------------------------------------------------------------
end subroutine fpf_find_tpf_PBM
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine var_fpf_find_tpf_PBM(d_a_par, c_par, var_d_a, var_c, cov_d_a_c, fpf , var_tpf, ierror)
!---------------------------------------------------------------------------------------------------
! PURPOSE: given a value for the FPF, find the variance in the estimate of the corresponding TPF, given the parameters d_a and c
! Warning the error currently does not return failed. Not clear yet how would it.
! CHANGES: LP April 30 2009 -> changed ierror = 2 to ierror = bad_input. Added also checkes on input values for variances

real(kind=double), intent(IN):: d_a_par, c_par ! parameters of the current fit, used as input
real(kind=double), intent(IN):: var_d_a, var_c, cov_d_a_c ! parameters of the current fit, used as input
real(kind=double), intent(IN):: fpf ! fixed, known value of the fpf
real(kind=double), intent(OUT):: var_tpf ! value of variance of  TPF for that fpf
integer, intent(OUT)::  ierror ! error flag. If ierror = bad_input, there is an input error, ierror = 3, the variances are unreliable

real(kind=double), parameter:: tol = 10e-6_double ! Accepted tolerance in the estimation of FPF
real(kind=double):: vc ! value of vc corresponding to FPF

real(kind=double):: d_tpf_d_d_a, d_tpf_d_c ! derivaties of TPF as a function of a fixed FPF in the two parameters

! Check for impossible values of the input parameters
if(fpf > 1.0_double .or. fpf < .0_double .or. d_a_par < .0_double .or. abs(c_par) > 1.0_double .or. &
   var_d_a < 0.0_double .or. var_c < 0.0_double .or. cov_d_a_c**2 > var_d_a * var_c ) then
  ierror = bad_input
  return
elseif(fpf > 0.999_double) then
  ierror = 3
elseif(fpf < .000001_double) then
  ierror = 3
else
  ierror = 1
endif

! Initialize vc, using a proxy for the infinite bound for the proper binormal model
vc = - sign(10.0_double, c_par)
! find the cutoff value that corresponds to FPF
call fpf_find_vc(d_a_par, c_par, fpf, tol,  vc)

! Compute the derivatives, needed for the delta method. note that density(vc) = - dTPF/dvc or -dFPF/dvc
! but since there are two densities the signs cancel, therefore the minus sign obtained from the derivation
! of a implicit function ( FPF_0 - FPF(da,c,vc) == 0) remains
d_tpf_d_d_a = d_tpf_d_d_a_PBM(d_a_par, c_par, vc)  - density_vc_PBM(d_a_par, c_par, vc, 1) * &
              d_fpf_d_d_a_PBM(d_a_par, c_par, vc) / density_vc_PBM(d_a_par, c_par, vc, 0)

d_tpf_d_c = d_tpf_d_c_PBM(d_a_par, c_par, vc)  - density_vc_PBM(d_a_par, c_par, vc,1) * &
              d_fpf_d_c_PBM(d_a_par, c_par, vc) / density_vc_PBM(d_a_par, c_par, vc,0)

! Apply the delta method
var_tpf = d_tpf_d_d_a**2 *var_d_a + d_tpf_d_c**2 *var_c + 2*d_tpf_d_d_a*d_tpf_d_c *cov_d_a_c

!---------------------------------------------------------------------------------------------------
end subroutine var_fpf_find_tpf_PBM
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine tpf_find_fpf_PBM(d_a_par, c_par, tpf , fpf, ierror)
!---------------------------------------------------------------------------------------------------
! PURPOSE: given a value for the TPF, find the corresponding FPF, given the parameters d_a and c
! Warning the error currently does not return failed. Not clear yet how would it
! April 30 2009 - LP -> changed ierror = 2 to ierror = bad_input.

real(kind=double), intent(in):: d_a_par, c_par ! parameters of the current fit, used as input
real(kind=double), intent(in):: tpf ! value of the tpf for the single point available for the fit
real(kind=double), intent(OUT):: fpf ! on return the value of the FPF for that TPF.
integer, intent(OUT)::  ierror ! error flag. If ierror = bad_input, the value of FPF is out of bounds

real(kind=double), parameter:: tol = 10e-6_double ! Accepted tolerance in the estimation of FPF
real(kind=double):: vc ! value of vc corresponding to FPF
real(kind=double):: junk

! Check for impossible values of the input parameters
if(tpf > 1.0_double .or. tpf < .0_double .or. d_a_par < .0_double .or. abs(c_par) > 1.0_double) then
  ierror = bad_input
  return
elseif(tpf > 0.99999_double) then
  fpf = 0.99999
  return
elseif(tpf < .00001_double) then
  fpf = 0.0_double
else
  ! Initialize vc using a proxy for the infinite bound of the proper binormal model
  vc = - sign(10.0_double, c_par)
  ! find the cutoff value that corresponds to FPF
  call tpf_find_vc(d_a_par, c_par, tpf, tol,  vc)
  call fpf_PBM(d_a_par, c_par, vc ,fpf, junk)
endif



!---------------------------------------------------------------------------------------------------
end subroutine tpf_find_fpf_PBM
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine var_tpf_find_fpf_PBM(d_a_par, c_par, var_d_a, var_c, cov_d_a_c, tpf , var_fpf, ierror)
!---------------------------------------------------------------------------------------------------
! PURPOSE: given a value for the TPF, find the variance in the estimate of the corresponding FPF, given the parameters d_a and c
! Warning the error currently does not return failed. Not clear yet how would it.
! CHANGES: LP April 30 2009 -> changed ierror = 2 to ierror = bad_input. Added also checkes on input values for variances

real(kind=double), intent(IN):: d_a_par, c_par ! parameters of the current fit, used as input
real(kind=double), intent(IN):: var_d_a, var_c, cov_d_a_c ! parameters of the current fit, used as input
real(kind=double), intent(IN):: tpf ! value of the tpf for the single point available for the fit
real(kind=double), intent(OUT):: var_fpf ! value of variance of the FPF for that tpf
integer, intent(OUT)::  ierror ! error flag. If ierror = bad_input, the value of FPF is out of bounds,
!                                          ierror = 3, the variances are unreliable

real(kind=double), parameter:: tol = 10e-6_double ! Accepted tolerance in the estimation of FPF
real(kind=double):: vc ! value of vc corresponding to TPF

real(kind=double):: d_fpf_d_d_a, d_fpf_d_c ! derivaties of FPF as a function of a fixed TPF in the two parameters

! Check for impossible values of the input parameters
if(tpf > 1.0_double .or. tpf < .0_double .or. d_a_par < .0_double .or. abs(c_par) > 1.0_double .or. &
   var_d_a < 0.0_double .or. var_c < 0.0_double .or. cov_d_a_c**2 > var_d_a * var_c ) then
       ierror = bad_input
       return
elseif(tpf > 0.999_double) then
       ierror = 3
elseif(tpf < .000001_double) then
        ierror = 3
else ! initialize error message
         ierror = 0
endif

! Initialize vc, using a proxy for the infinite bound for the proper binormal model
vc = - sign(10.0_double, c_par)
! find the cutoff value that corresponds to TPF
call tpf_find_vc(d_a_par, c_par, tpf, tol,  vc)

! Compute the derivatives, needed for the delta method. note that density(vc) = - dTPF/dvc or -dFPF/dvc
! but since there are two densities the signs cancel, therefore the minus sign obtained from the derivation
! of a implicit function ( TPF_0 - TPF(da,c,vc) == 0) remains



d_fpf_d_d_a = d_fpf_d_d_a_PBM(d_a_par, c_par, vc)  - density_vc_PBM(d_a_par, c_par, vc, 0) * &
              d_tpf_d_d_a_PBM(d_a_par, c_par, vc) / density_vc_PBM(d_a_par, c_par, vc, 1)


d_fpf_d_c = d_fpf_d_c_PBM(d_a_par, c_par, vc)  - density_vc_PBM(d_a_par, c_par, vc,0) * &
              d_tpf_d_c_PBM(d_a_par, c_par, vc) / density_vc_PBM(d_a_par, c_par, vc,1)


! Apply the delta method
var_fpf = d_fpf_d_d_a**2 *var_d_a + d_fpf_d_c**2 *var_c + 2*d_fpf_d_d_a*d_fpf_d_c *cov_d_a_c

!---------------------------------------------------------------------------------------------------
end subroutine var_tpf_find_fpf_PBM
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine fpf_find_vc(d_a_par, c_par, p, tol,  vc)
! Purpose: find the value of vc that produces a specific fpf, given d_a and c
! NOTE: error messaging should be introduced, the routine could use a more subtle non linear method
!       but i don't care enough for this problem to do that.
! NOTE2: the routine assumes that the value of c that is sent to it is not zero. this will be
!       eventually corrected (error messages could do the trick)
! WARNING: routine doesn't deal with wrong input or input right at tbe boundary (FPF = 0 or 1)
! HISTORY: 1-22-06 introduced if statements to deal with c =~ 0.
!---------------------------------------------------------------------------------------------------
use proproc_one_pt ! module that contains the parameter values for the function called by the root finder
use gen_numerics, only: zbren !  root finder
use statistic_functions, only : phi
implicit none

real(kind=double), intent(in):: d_a_par, c_par ! parameters of the current fit, used as input
real(kind=double), intent(in):: p ! value of the fpf for the single point available for the fit
real(kind=double), intent(in):: tol ! absolute error tolerance
real(kind=double), intent(inout)::  vc ! cutoff value that satisfies FPF(d_a, c, vc ) == p

!internal variables
real(kind = double) :: fpf,  one_minus! value of the fpf and its complement, for a value of vc
real(kind=double):: vc_up, vc_down ! 2 values of vc that bracket the root FPF(d_a, c, vc ) == p
real(kind=double)::  delta_vc ! The step to be used to seek the bracketing values
real(kind=double)::  vc_bound ! maximum or minimum value of vc allowed by the proper model (see proproc paper)
real(kind=double)::  p_old ! p_old is used to check if fpf changes
integer :: ierror ! Error flag, might need to be returned to the calling program if this computation becomes
         ! problematic

! the finite bound of the cutoffs values, following mathematica definition, unless c is very close to 0. In which
! case we use the properties of normal density adding the precision we want to use (so we don't need to
! check further than 10 sigmas (in fact 8 would be enough too).
 if( abs(c_par) > 1.0e-2_double) then
     vc_bound = d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * sign( max(1.0e-8_double, abs(c_par)),c_par))
 else
     vc_bound = sign(10.0_double, c_par)
 endif

! compute the FPF for the initial value of vc, as defined by the calling program
 call fpf_PBM(d_a_par, c_par, vc ,fpf, one_minus)

! bracket the zero in vc for the equation:  FPF(d_a, c, vc ) == p
 if (fpf  < p) then
      vc_up = vc ! if the value computed now is too small it means that this vc is too large
      ! so we have the upper bound (remember that FPF is a monotonically decreasing function in vc

      ! Check whether we are moving towards a finite bound and set the incremental step accordingly, in
      ! order not to step out of bounds and not to take forever. Note that the step is regulated by the
      ! value of c
      if(c_par >= 1.0e-2_double ) then
             delta_vc = - .1_double/ abs( 1.0_double - sign( min(.999999999_double, abs(c_par)),c_par))
      else ! -abs is to make sure that rounding does not make the sign wrong
             delta_vc = - abs(vc_bound - vc_up)/100
      endif
!     since the increment is negative we sum it and look for the lower bound
      vc_down = vc_up + delta_vc

      do ! dangerous infinte loop, but so far we didn't get stuck
            p_old = fpf
            call fpf_PBM(d_a_par, c_par, vc_down ,fpf, one_minus)
            ! It can happen that we either hit the value of that the value sought after is
            ! a boundary value, like for fpf=0. In that case it can't be bracketed, but it
            ! can be found
            if (fpf  .speq. p)  then
              vc  = vc_down
              return
            endif
            if (fpf  > p)  exit ! We found the lower bound in vc, so we exit the loop
            vc_down = vc_down + delta_vc ! Increase vc_down to bracket
            ! Check whether the values have been changing any, otherwise accelerate the change, but make sure not to be
            ! too close to the mathematical limits for p, vc, and so on. The step is increasing exponentially.
            if ( abs(fpf - p_old)/max(p_old, 1.0e-6_double) < .1_double .and. (p < .99_double .or. p > .1_double) ) &
                  delta_vc = delta_vc*10
      enddo
 else
      vc_down = vc  ! if the value computed now is too large it means that this value of vc is too small
      ! so we have the lower bound (remember that FPF is a monotonically decreasing function in vc

      ! Check whether we are moving towards a finite bound and set the incremental step accordingly, in
      ! order not to step out of bounds and not to take forever. Note that the step is regulated by the
      ! value of c
      if(c_par < 1.0e-2_double ) then
             delta_vc =    .1/ abs( 1.0_double - sign( min(.999999999_double, abs(c_par)),c_par))
      else
             delta_vc =  abs(vc_bound - vc_down)/100
      endif
             vc_up = vc_down + delta_vc

      do ! infinite loop!
            p_old = fpf
            call fpf_PBM(d_a_par, c_par, vc_up ,fpf, one_minus)
            if (fpf  < p)  exit
            ! It can happen that we either hit the value of that the value sought after is
            ! a boundary value, like for fpf=0. In that case it can't be bracketed, but it
            ! can be found
            if (fpf  .speq. p)  then
              vc  = vc_up
              return
            endif
            vc_up = vc_up + delta_vc ! Increase vc_up to bracket
            ! Check whether the values have been changing any, otherwise accelerate the change
            if ( abs(fpf - p_old)/max(p_old,1.0e-6_double) < .1_double .and. (p < .99_double .or. p > .1_double) ) &
                 delta_vc = delta_vc*10
      enddo

 endif

! Load the values in the module that sends the data to the function through
! the optimizer
  d_a = d_a_par
  c = c_par
! find the root
  call  zbren( p , fpf_vc, vc_up,  vc_down,  tol, vc ,ierror)


!---------------------------------------------------------------------------------------------------
end subroutine fpf_find_vc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine tpf_find_vc(d_a_par, c_par, p, tol,  vc)
! PURPOSE: find the value of vc that produces a specific tpf, given d_a and c
! NOTE: error messaging should be introduced, the routine could use a more subtle non linear method
!       but i don't care enough for this problem to do that.
! NOTE2: the routine assumes that the value of c that is sent to it is not zero. this will be
!       eventually corrected (error messages could do the trick)
! WARNING: routine doesn't deal with wrong input or input right at tbe boundary (TPF = 0 or 1)
! HISTORY: 1-22-06 introduced if statements to deal with c =~ 0.
!---------------------------------------------------------------------------------------------------
use proproc_one_pt ! module that contains the parameter values for the function called by the root finder
use gen_numerics, only: zbren !  root finder
implicit none

real(kind=double), intent(in):: d_a_par, c_par ! parameters of the current fit, used as input
real(kind=double), intent(in):: p ! value of the tpf for the single point available for the fit
real(kind=double), intent(in):: tol ! absolute error tolerance
real(kind=double), intent(inout)::  vc ! cutoff value that satisfies tpf(d_a, c, vc ) == p

!internal variables
real(kind = double) :: tpf,  one_minus! value of the tpf and its complement, for a value of vc
real(kind=double):: vc_up, vc_down ! 2 values of vc that bracket the root tpf(d_a, c, vc ) == p
real(kind=double)::  delta_vc ! The step to be used to seek the bracketing values
real(kind=double)::  vc_bound ! maximum or minimum value of vc allowed by the proper model (see proproc paper)
real(kind=double)::  p_old ! p_old is used to check if tpf changes
integer :: ierror ! Error flag, might need to be returned to the calling program if this computation becomes
 ! problematic

! Find the finite bound of the cutoffs values, unless c is very close to 0, in which
! case we use the properties of normal density: we don't really need to
! check further than 10 sigmas (in fact 8 would be enough too).
 if( abs(c_par) > 1.0e-2_double) then
     vc_bound = d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * sign( max(1.0e-8_double, abs(c_par)),c_par))
 else
     vc_bound = sign(10.0_double, c_par)
 endif


! compute the tpf for the initial value of vc, as defined by the calling program
 call tpf_PBM(d_a_par, c_par, vc ,tpf, one_minus)

! bracket the zero in vc for the equation:  tpf(d_a, c, vc ) == p
 if (tpf  < p) then
      vc_up = vc ! if the value computed now is too small it means that this vc is too large
      ! so we have the upper bound (remember that tpf is a monotonically decreasing function in vc

      ! Check whether we are moving towards a finite bound and set the incremental step accordingly, in
      ! order not to step out of bounds and not to take forever. Note that the step is regulated by the
      ! value of c
      if(c_par >= 0 ) then
             delta_vc = -  .1/ abs( 1.0_double + sign( min(.999999999_double, abs(c_par)),c_par))
      else
             delta_vc = (vc_bound - vc_up)/100
      endif
!     since the increment is negative we sum it and look for the lower bound
      vc_down = vc_up + delta_vc

      do ! dangerous infinte loop, but so far we didn't get stuck
            p_old = tpf
            call tpf_PBM(d_a_par, c_par, vc_down ,tpf, one_minus)
            ! It can happen that we either hit the value of that the value sought after is
            ! a boundary value, like for fpf=0. In that case it can't be bracketed, but it
            ! can be found
            if (tpf  .speq. p)  then
              vc  = vc_down
              return
            endif
            if (tpf  > p)  exit ! We found the lower bound in vc, so we exit the loop
            vc_down = vc_down + delta_vc ! Increase vc_down to bracket
            ! Check whether the values have been changing any, otherwise accelerate the change, but make sure not to be
            ! too close to the mathematical limits for p, vc, and so on. The step is increasing exponentially.
            if ( abs(tpf - p_old)/p_old < .1_double .and. (p < .99_double .or. p > .1_double) )  delta_vc = delta_vc*10
      enddo
 else
      vc_down = vc  ! if the value computed now is too large it means that this value of vc is too small
      ! so we have the lower bound (remember that tpf is a monotonically decreasing function in vc

      ! Check whether we are moving towards a finite bound and set the incremental step accordingly, in
      ! order not to step out of bounds and not to take forever. Note that the step is regulated by the
      ! value of c
      if(c_par < 0 ) then
             delta_vc =    .1/ abs( 1.0_double + sign( min(.999999999_double, abs(c_par)),c_par))
      else
             delta_vc =  (vc_bound - vc_down)/100
      endif
             vc_up = vc_down + delta_vc

      do ! infinite loop!
            p_old = tpf
            call tpf_PBM(d_a_par, c_par, vc_up ,tpf, one_minus)
            ! It can happen that we either hit the value of that the value sought after is
            ! a boundary value, like for fpf=0. In that case it can't be bracketed, but it
            ! can be found
            if (tpf  .speq. p)  then
              vc  = vc_down
              return
            endif
            if (tpf  < p)  exit
            vc_up = vc_up + delta_vc ! Increase vc_up to bracket
            ! Check whether the values have been changing any, otherwise accelerate the change
            if ( abs(tpf - p_old)/p_old < .1_double .and. (p < .99_double .or. p > .1_double) )  delta_vc = delta_vc*10
      enddo

 endif

! Load the values in the module that sends the data to the function through
! the optimizer
  d_a = d_a_par
  c = c_par
! find the root
  call  zbren( p , tpf_vc, vc_up,  vc_down,  tol, vc ,ierror)


!---------------------------------------------------------------------------------------------------
end subroutine tpf_find_vc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
real(kind=double) function tpf_vc(vc)
!---------------------------------------------------------------------------------------------------
! PURPOSE: compute the tpf only as a function of a cutoff, d_a and c are passed through a module
!          it is meant to be used with an optimizer that calls a single variable functions and can't
!          (and probably shouldn't) see the values of the other parameters
   use proproc_one_pt
   implicit none

   real(kind=double), intent(in):: vc ! cutoff value
   real(kind=double):: one_minus ! one minus tpf, not used, but computed. In fact it is optional, but
                       ! some compiler choke on it (compiler is wrong) so we wrote it this way

   call tpf_PBM(d_a, c, vc ,tpf_vc, one_minus)

!---------------------------------------------------------------------------------------------------
 end function tpf_vc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
real(kind=double) function fpf_vc(vc)
!---------------------------------------------------------------------------------------------------
! PURPOSE: compute the fpf only as a function of a cutoff, d_a and c are passed through a module
!          it is meant to be used with an optimizer that calls a single variable functions and can't
!          (and probably shouldn't) see the values of the other parameters
   use proproc_one_pt
   implicit none

   real(kind=double), intent(in):: vc ! cutoff value
   real(kind=double):: one_minus ! one minus fpf, not used, but computed. In fact it is optional, but
                       ! some compiler choke on it (compiler is wrong) so we wrote it this way

   call fpf_PBM(d_a, c, vc ,fpf_vc, one_minus)

!---------------------------------------------------------------------------------------------------
 end function fpf_vc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine one_pt_d_a_zero_fit(p, q, d_a_par, c_par, vc, ierror, err_msg)
! Purpose: compute the exact fit through one point using the proper model. Note that the one pt
! fit form is formally identical for c negative and c positive, if  p -> p + 1.
! ALGORITHM: take the equations for the proper model, for FPF and TPF. Set d_a == 0, the resulting equations
!            can be inverted so that one can create a function f(TPF) = H ( f(FPF),  c), which allows to solve for
!            c
!---------------------------------------------------------------------------------------------------
use statistic_functions, ONLY : compute_zdev_plus

implicit none

real(kind = double), intent(in) :: p ! orizontal data coordinate (FPF)
real(kind = double), intent(in) :: q ! vertical data coordinate (TPF)
real(kind=double), intent(out):: d_a_par, c_par, vc ! parameters of the fit
integer, intent(out)                                   :: ierror ! Error flag
character(len = line_length), intent(out):: err_msg    ! description of the error occurred


real(kind=double):: first_arg_dev, second_arg_dev ! normal deviates for the 2 arguments in the proper model


! Determine accuratly the first deviate necessary to compute the exact solution
call compute_zdev_plus(p, first_arg_dev, ierror)
if(ierror /= 0) then
  ierror = 1
  err_msg = " Failed to compute the exact fit for d_a = 0 "
  return
endif

! Determine the second deviate accurately
call compute_zdev_plus(q, second_arg_dev, ierror)
if(ierror /= 0) then
  ierror = 1
  err_msg = " Failed to compute the exact fit for d_a = 0"
  return
endif

! analitical form of the root in vc
vc  =  - ( first_arg_dev +  second_arg_dev ) / 2.0_double

! analitical form of the root in c
c_par = first_arg_dev / vc + 1.0_double

! as assumed the value of d_a is zero
d_a_par = 0.0_double


end subroutine one_pt_d_a_zero_fit
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
subroutine pbmroc_mle(mn, ms, num_categ, catn_in, cats_in, idebug, &
                       d_a_par, c_par, auc, variance_auc, vc_cutoffs_out, log_like, ierror, &
                       cov_out, hessian_out)
! PURPOSE: compute the MLE estimate using the proper binormal model starting from categorial data
!          the program assumes that positivity (higher likelihood of having a positive reading) is
!          for categories with larger indices. If this is not the case, the categorical data need to
!          be reverse before feeding them to this subroutine.
!ALGORITHM:
!          call compute_MLE that returns the MLE estimate with the hessian matrix, that is used to
!          produce the standard errors using the Kramer-Rao bound.
!          Since there can be multiple maxima (mostly for very unbablanced datasets), the algorithms first
!          attempts using 3 points spread in the d_a, c space, then if 2 maxima are found, it looks for a third
!          in between them.
! NOTE:    The multiple maxima algorithm is important for most dataset situations since there is a probability
!          that there will be multiple maxima (about 5-20% for datasets smaller than 50+50 cases, then
!          a little smaller than  .5% for larger ones). On the other hand, the search for a triple maximum is
!          a safety add on. Triple maxima are somewhat more rare, so in case large speed increases were
!          needed, this part could be dropped.
! HISTORY: 03-14/06 LP, UC, introduced checking of input, specifically if the categorical input data contains
!                   only empty cells or contains  negative values.
!         03-21-06 LP, UC, introduce a collapser to take care of non reduced to corners data. Note that we return
!                   the reduced hessian and variance, that is the hessian and variance referred to the collapsed data
!                   one could conceivable build the hessian for the original system, but since we cannot see much
!                   use for it (in the sense that people who study the hessian presumable will send the collapsed
!                   data) we are not going to spend time to figure out how to do it and programming it.
!         06-04-08 LP, UC removed err_msg from the return call to make calling across languages less problematic
!         28-12-09 LP, UC added output to log file in case the input is unacceptable, to harmonize with cvbm
use debugging ! module containing the procedures to log the debugging information on the run
use problem_data ! where the ROC basic categorical information is stored.
use statistic_functions, ONLY : compute_zdev_plus
use categorization, ONLY:  LABROC4_Collapser

implicit none

! Subroutine arguments
integer, intent(in):: mn ! number of actually negative cases
integer, intent(in):: ms ! number of actually positive cases
integer, intent(in) :: num_categ ! Number of categories as created by catgrz
integer, dimension(num_categ), intent(in):: catn_in, cats_in ! arrays containing categorical data
integer, intent(in) :: idebug               ! 0 = no debug; 1  = debug
real(kind = double),intent(out)                        :: d_a_par, c_par ! MLE of the parameters
real(kind = double),intent(out)                        :: auc ! AUC, area under the curve
real(kind = double),intent(out)                        :: variance_auc ! estimated variance of AUC

real(kind = double), dimension(num_categ-1), intent(out) :: vc_cutoffs_out  ! cutoff parameter values at the maximum found
real(kind = double), intent(out)                       :: log_like ! value of the log likelihood function at the final point
integer, intent(out)                                   :: ierror ! Error flag about the MLE fit
! Note that the error values are set in this routine, or initialize_d_a_c so if their numbers have
! to be changed, they have to be changed here, the rest of the subroutines use their own numbering
! specific per routine. Look above for the different meaning. Not only 0 is successful fit
! June 2008-LP removed err_msg from the return call to prevent problems with calls from R and other programs
! character(len = line_length), intent(out):: err_msg    ! description of the error occurred. If the fit worked, it
                                                       ! contains information about how were the variances computed
real(kind = double), dimension(num_categ+1,num_categ+1), intent(out)  :: cov_out ! these are used because the number
real(kind = double), dimension(num_categ+1,num_categ+1), intent(out), optional  :: hessian_out ! of categories can be
     ! different inside proproc because of collapsing. Note that only the reduced hessian will be returned, the
     ! rest will be set to garbage

! Internal variables
integer:: i,j ! loop counter for logging warnings about cutoffs
character(len = line_length):: msg    ! character string used to display state information
character(len = line_length):: err_msg    ! character string used to display state information
! The algorithms first choses 3 points, to span the configuration space for multiple maxima, then
! if more than one maximum is found, it looks for additional maxima in between the 2 most distant
! maxima. Multiple Maxima problems aren't really solvable, but this fix seems to work
integer, parameter:: num_init_pts = 3 ! Number of initial estimates
integer, parameter:: num_pts_bt = 5 ! number of points to put in between the most different final estimates,
                                    ! used only if needed
integer, dimension(num_init_pts + num_pts_bt )  :: pts_error ! whether a specific fit has failed
real(kind = double), dimension(num_init_pts + num_pts_bt )  :: d_a_par_array, c_par_array ! MLE of the parameters
real(kind=double), dimension(num_init_pts + num_pts_bt ) ::  log_like_array ! the log like for the 2 possible initial estimates


! Some variables used if multiple maxima are found
real(kind = double) :: d_a_min, d_a_max, c_min, c_max ! if final estimates vary a lot, take the most different and
           ! look in between to make sure that there aren't multiple maxima.
real(kind=double):: auc_min, auc_max ! areas for the 2 extremal maxima
real(kind = double):: da1,da2,c1,c2 ! Coefficients for the curve to find points in between very different maxima
real(kind=double):: d_a_phase ! for the ellipsis, whether it is pi or 2 pi, depending upon the value of d_a and c
                          ! for the 2 extremal fits
real(kind=double):: temp
integer, dimension(num_categ):: new_cat_index ! if the data needs to be collapsed, the list of where to
                                              ! put the old data
integer:: icat ! loop counter
! These are the internal variables, defined using num_cat, the internal number of categories after
! collapsing
real(kind = double), allocatable, dimension(:,:)  :: cov ! internal values, can be smaller
real(kind = double), allocatable, dimension(:,:)  :: hessian ! if collapsed
real(kind = double), allocatable, dimension(:)    :: vc_cutoffs  ! cutoff parameter values at the maximum found
real(kind = double), allocatable, dimension(:, : ) :: vc_cutoffs_array  ! MLE cutoff parameter values

! Debug is used to pass the debugging flag around optimizers that have no business knowing its value
if(idebug == 1) then
 debug = .true.
else
 debug = .false.
endif

! First check whether the input is acceptable, that is it can't contain negative number
! of cases in any category for any truth state and the number of actually positive and
! actually negative cases in each category can't be  equal to zero for all categories
if(  ( any(catn_in < 0) .or. any(cats_in < 0) ) .or.  & ! negative # of cases
     all( (catn_in + cats_in) == 0 ) .or.  & ! categorical data arrays empty
     (sum(catn_in) /= mn .or. sum(cats_in) /= ms ) & ! mismatch in # of cases
     )then
    ierror = bad_input
    if(idebug == 1) then
        continue! write(msg,*) "Categorical data given to proproc is unacceptable"
        call DisplayState(msg)
        continue! write (msg,*) "CATEGORICAL DATA USED IN THE MLE, total cat = ", num_categ
        call DisplayState(msg)
        do i = 1, num_categ
           continue! write (msg,*) catn_in(i), cats_in(i)
           call DisplayState(msg)
        enddo
        continue! write (msg,*) "END OF CATEGORICAL DATA"
        call DisplayState(msg)
    endif
    ! Mathematically it should be undetermined, but we put an impossible value
    ! to prevent interfacing problems from arising
    if (present(hessian_out))  hessian_out = -666.0_double
    d_a_par = -666.0_double
    c_par = -666.0_double
    auc = -666.0_double
    variance_auc    = -666.0_double
    vc_cutoffs_out     = -666.0_double
    log_like        = -666.0_double
    cov_out         = -666.0_double
    return
endif
! COLLAPSE (IF NEEDED) AND LOAD THE COLLAPSED DATA INTO NEW ARRAYS
! Prepare the variables that will be stored in the moduli from file "proproc_other_modules.f90"
! Since those are in a module while these are explicit, we seem to need to double the
! names at this level. The moduli are used to work around the various optimizers (initial estimates
! of cutoffs and full MLE). NOTE: all the calls that don't have to "pass-through" an optimizer are
! made explicit, even if this might seem confusing, I prefer explicit calls when I can do them.

call LABROC4_Collapser(catn_in, cats_in, num_categ, idebug, new_cat_index)

num_cat = new_cat_index(num_categ) ! the total number of collapsed categories is the category assigned to the
       ! last category
allocate(catn(num_cat),  cats(num_cat))
allocate(    cov(num_cat + 1, num_cat + 1) )
allocate(hessian(num_cat + 1, num_cat + 1) )
allocate( vc_cutoffs(num_cat - 1) )
allocate( vc_cutoffs_array(num_cat - 1,num_init_pts + num_pts_bt ) )

num_normal_cases   = mn
num_abnormal_cases = ms

catn = 0
cats = 0
pass_idebug  = idebug

! Collapse the data
do icat = 1, num_categ
  catn(new_cat_index(icat)) =   catn(new_cat_index(icat)) +  catn_in(icat)
  cats(new_cat_index(icat)) =   cats(new_cat_index(icat)) +  cats_in(icat)
enddo


pts_error = fit_OK ! initialize array of error messages for point fits

! get 3 initial estimates for d_a and c. Notice that initialize_d_a_c always returns 3, so num_init_pts isn't a
! variable that can be changed here unless we change that procedure as well
call initialize_d_a_c(mn, ms, num_cat, catn, cats, idebug, d_a_par_array(1:num_init_pts), c_par_array(1:num_init_pts), &
                      ierror, err_msg)


initial_est_OK: if(ierror == fit_OK) then
 initial_estimates_loop: do j = 1, num_init_pts ! Loop over the num_init_pts initial estimates
    ! Load the jth initial estimate
    d_a_par = d_a_par_array(j)
    c_par = c_par_array(j)

    call initialize_cutoffs(d_a_par,c_par, mn, ms, num_cat, catn, cats, vc_cutoffs, idebug, ierror, err_msg)
    ! Note that here we might convert local error flags with proproc_<vers> error flags
    initial_est_cutoffs: if(ierror /= 0) then ! Initial estimates of cutoffs failed
             ! report a failure in the initialization of the cutoffs
             ierror = fit_fail_init_est ! failed initial estimates
             err_msg = " initialize_cutoffs, " // err_msg(:line_length - 21)
             pts_error(j) = 1 ! this specific value failed
    else ! Initial estimates of cutoff successful
             ! Check if there are any 2 cutoffs which have almost identical values
             ! Normally the algorithm deals with it. It never happens for the initial estimates
             ! since we changed the algorithm - LP summer 2005 - UC
             if(idebug == 1 ) then
                 if( any(vc_cutoffs(2:num_cat-1) .speq. vc_cutoffs(1:num_cat-2))  ) then
                         call DisplayState( ' WARNING:: initialize_cutoffs: Some cutoffs almost identical')
                         continue! write(msg,"('          for d_a & c equal to::',2(1x,d14.6))") d_a_par, c_par
                         call DisplayState(msg)
                         do i = 1, num_cat-1
                                        continue! write(msg,"(i3,1(1x,e22.15))") i, vc_cutoffs(i)
                                        call DisplayState(msg)
                        enddo
                  endif
              endif
             ! Since we checked both the initial estimates of the parameters and of the cutoffs, we can do MLE
             call compute_MLE(mn, ms, num_cat, catn, cats,idebug, d_a_par, c_par,  &
                              vc_cutoffs, log_like, hessian, ierror, err_msg)
             if(ierror == 0) then
                 ierror = fit_OK
             else
                 pts_error(j) = 1 ! this specific value failed
                 if(idebug == 1) then
                    continue! write(msg,*) "WARNING: Optimization for initial estimate ", num_init_pts, " failed "
                    call DisplayState(msg)
                    call DisplayState(err_msg)
                 endif
                 ierror = fit_failed
             endif
    endif initial_est_cutoffs
    ! Log the current estimate
    if(idebug == 1) then
       continue! write(msg, *) "For initial estimates:", d_a_par_array(j), c_par_array(j)
       call DisplayState(msg)
    endif

    ! replace the initial estimates with the current MLE estimates
    d_a_par_array(j) =  d_a_par
    c_par_array(j)   =  c_par
    log_like_array(j) = log_like
    vc_cutoffs_array(1:num_cat-1,j) = vc_cutoffs(1:num_cat-1)

    if(idebug == 1) then
       if( pts_error(j) == fit_OK) then ! check if the fit worked for  this point
         continue! write(msg, *) "Values of D_a, c, and Log Likelihood at Local Maximum"
         call DisplayState(msg)
         continue! write(msg, *)  d_a_par, c_par, log_like
         call DisplayState(msg)
       else
         continue! write(msg, *) "  NO FIT WAS PRODUCED because"
         call DisplayState(msg)
         continue! write(msg, *) " ", err_msg(1:76)
         call DisplayState(msg)
       endif
    endif

 enddo initial_estimates_loop

 ! Check whether there was at least one initial
 ! estimate that produced an acceptable fit
 ! (normally at most one fails, and only for very
 !  bad datasets)
 any_fit: if( .not. any(pts_error == fit_OK) ) then
     ierror = fit_failed
     err_msg = "No initial estimates converged to a maximum of Log Likelihood"
 else
   ! NOTE: We are leaving the value of the failed fit so that the algorithm can search
   !       around it, in case we have missed a true maximum because one of the initial
   !       estimates wasn't too good (usually this means too skewed)

   ! verify that the final estimates are quite similar, if they aren't, it might be that
   ! there are multiple maxima (not just 2), so we need to redo some calculations between those final
   ! estimates, to see if the final results is between them or not
   if (   &
      maxval(d_a_par_array(1:num_init_pts)) - minval(d_a_par_array(1:num_init_pts))  > .1_double   .or.   &
      abs( maxval(c_par_array(1:num_init_pts)) - minval(c_par_array(1:num_init_pts)) )  > .1_double &
   ) then



      ! We use c to detemine the 2 most extreme value, normally d_c and c vary
      ! inversely (so it would be the same), but c is a more reliable parameter
      d_a_min  = d_a_par_array ( minloc( c_par_array(1:num_init_pts), dim = 1) )
      c_min    = c_par_array   ( minloc( c_par_array(1:num_init_pts), dim = 1) )

      d_a_max  = d_a_par_array ( maxloc( c_par_array(1:num_init_pts), dim = 1) )
      c_max    = c_par_array   ( maxloc( c_par_array(1:num_init_pts), dim = 1) )

      if(idebug == 1) then
         continue! write(msg,*) "ATTENTION: Log Likelihood might display multiple Maxima "
         call DisplayState(msg)
         continue! write(msg,*) " Maximal final estimate ", d_a_max, c_max
         call DisplayState(msg)
         continue! write(msg,*) " Minimal final estimate ", d_a_min, c_min
         call DisplayState(msg)
      endif

      ! Create max_pts_bt points between the 2 farthest initial estimates, they sit on
      ! an ellipse. Note that if both maxima are on the d_a axis, we need to draw a different
      ! kind of ellipse (and there are 2 options, whether the points have the same value
      ! of c or not)
      if (d_a_max < .5_double * abs(c_max) .and.  d_a_min < .5_double * abs(c_min) .and. &
          c_max*c_min < 0.0_double ) then ! c's have the opposite signs
          call auc_PBM(d_a_min, c_min, auc_min, ierror)
          call auc_PBM(d_a_max, c_max, auc_max, ierror)
          call compute_zdev_plus( (auc_min + auc_max)*.5_double, temp, ierror)
          da1 = d_a_min
          da2 = temp * sqrt(2.0_double) - d_a_min
          c1  = c_max
          c2 = c_min - c_max
          d_a_phase = 2.0_double
      elseif (d_a_max > d_a_min) then
          da1 = d_a_min
          da2 = d_a_max - d_a_min
          c1  = c_max
          c2 = c_min - c_max
          d_a_phase = 1.0_double
      else
          da1 = d_a_max
          da2 = d_a_min - d_a_max
          c1  = c_min
          c2 = c_max - c_min
          d_a_phase = 1.0_double
      endif

      if(idebug == 1) then
          continue! write(msg,*) "Start with new set of initial points between extremal final estimates "
          call DisplayState(msg)
      endif
      do j = 1, num_pts_bt ! Compute points on an ellipse between the most different final estimates
         d_a_par = da1 + da2 * cos ( j* pi * d_a_phase/ ( 2 * (1 + num_pts_bt )) - (d_a_phase - 1) * pi / 2 )
         c_par  = c1 + c2 * sin ( j* pi/ ( 2 * (1 + num_pts_bt )) )

         if(idebug == 1) then
             continue! write(msg,*) "  ", d_a_par, c_par
             call DisplayState(msg)
         endif

         call initialize_cutoffs(d_a_par,c_par, mn, ms, num_cat, catn, cats, vc_cutoffs, idebug, ierror, err_msg)
         ! Note that here we might convert local error flags with proproc_<vers> error flags
         if(ierror /= 0) then ! Initial estimates of cutoffs failed
             ! report a failure in the initialization of the cutoffs
             ierror = fit_fail_init_est ! failed initial estimates
             err_msg = " initialize_cutoffs, " // err_msg(:line_length - 21)
         else ! Initial estimates of cutoff successful
             ! Check if there are any 2 cutoffs which have almost identical values
             ! Normally the algorithm deals with it. It never happens for the initial estimates
             ! since we changed the algorithm - LP summer 2005 - UC
             if(idebug == 1 ) then
                 if( any(vc_cutoffs(2:num_cat-1) .speq. vc_cutoffs(1:num_cat-2))  ) then
                         call DisplayState( 'WARNING:: initialize_cutoffs: Some cutoffs almost identical')
                         continue! write(msg,"('          for d_a & c equal to::',2(1x,d14.6))") d_a_par, c_par
                         call DisplayState(msg)
                         do i = 1, num_cat-1
                                        continue! write(msg,"(i2,3(1x,e22.15))") i, vc_cutoffs(i)
                                        call DisplayState(msg)
                        enddo
                  endif
             endif
             ! Since we checked both the initial estimates of the parameters and of the cutoffs, we can do MLE
             call compute_MLE(mn, ms, num_cat, catn, cats,idebug, d_a_par, c_par,  &
                              vc_cutoffs, log_like, hessian, ierror, err_msg)
             if(ierror == 0) then
                 ierror = fit_OK
             else
                 ierror = fit_failed
                 if(idebug == 1) then
                    continue! write(msg,*) "WARNING: Optimization for initial estimate ", num_init_pts + j, "failed "
                    call DisplayState(msg)
                    call DisplayState(err_msg)
                 endif
             endif
         endif

         d_a_par_array(num_init_pts + j)                =  d_a_par
         c_par_array(num_init_pts + j)                  =  c_par
         log_like_array(num_init_pts + j)               = log_like
         vc_cutoffs_array(1:num_cat-1,num_init_pts + j) = vc_cutoffs(1:num_cat-1)

         if(idebug == 1) then
             continue! write(msg, *) "Values of D_a, c, and Log Likelihood at Local Maximum"
             call DisplayState(msg)
             continue! write(msg, *)  d_a_par, c_par, log_like
             call DisplayState(msg)
         endif


     enddo

     ! Now select the best between all the computed points, dim = 1 is to force output to be scalar (to avoid rank mismatches)
     j = maxloc(log_like_array, dim = 1)  ! determine which one  has the largest Log Likelihood
     d_a_par    = d_a_par_array(j)
     c_par      = c_par_array(j)
     vc_cutoffs(1 : num_cat - 1) = vc_cutoffs_array(1:num_cat-1,j)

   else

     ! Now select the best, note that we use only the first 3, since we didn't notice any multiple maxima
     ! dim = 1 is to force output to be scalar (to avoid rank mismatches)
     j = maxloc(log_like_array(1:num_init_pts), dim = 1)  ! determine which one of the 2 has the largest Log Likelihood
     d_a_par    = d_a_par_array(j)
     c_par      = c_par_array(j)
     vc_cutoffs(1 : num_cat - 1) = vc_cutoffs_array(1:num_cat-1,j)

   endif

 ! rerun from the maximum to recompute the hessian and make sure that the system is truly consistent
 ! This hardly costs any time at all, since all the parameters are already optimal, it runs in a single
 ! step
 call compute_MLE(mn, ms, num_cat, catn, cats,idebug, d_a_par, c_par,  &
                   vc_cutoffs, log_like, hessian, ierror, err_msg)

 ! Log the chosen value
 if(idebug == 1) then
       continue! write(msg, *) "Final value of D_a, c, and Log Likelihood"
       call DisplayState(msg)
       continue! write(msg, *)  d_a_par, c_par, log_like
       call DisplayState(msg)
       continue! write(msg, *) "Final value of cutoffs:"
       call DisplayState(msg)
       do j = 1, num_cat -1
            continue! write(msg, *)  j, vc_cutoffs(j)
            call DisplayState(msg)
       enddo
  endif

endif any_fit ! sequence done only if at least one of the initial points worked out

endif initial_est_OK ! end of sequence done if we could produce initial estimates



! Estimate the area under the curve (AUC) and other parameters, depending upon the error flag
!Initialize arrays to non-sense in case quantities can't be computed, applies to all values
! of error flags
if (present(hessian_out))  hessian_out = -666.0_double
cov_out = -666.0_double
vc_cutoffs_out = -666.0_double


select case (ierror)
case(fit_OK) ! fit successful, compute the area
     call auc_PBM(d_a_par, c_par, auc, ierror)
     ! If ierror = 0, it is OK, if it is == 2, it jsut mean that c is too small, so it is also OK
     ! only think worthy reporting is when it is 1, which means that F is truly unstable. It never happens.
     if(ierror == 1 ) then
           continue! write(err_msg, *) "Computation of auc with the proper model failed"
           call print_warning_line (err_msg)
     endif
     call compute_covariance_matrix(d_a_par, c_par, num_cat+1, idebug, hessian, cov, ierror, err_msg)
    if(auc <= .50005_double) then ! One can't compute variances using the hessian when auc is too small
                 ! No matter what is the outcome of the diagnostics.
                 err_msg = " Variances cannot be computed because AUC is too small "
                 ierror = fit_fail_var_auc_small
     elseif(ierror == 1) then ! Variances could not be computed, also convert to error numbering for proproc to return
          if(idebug == 1) call DisplayState(err_msg)
                 ierror = fit_fail_variances
     else ! Variances can be computed
          if( ierror == 2) then
                ierror = fit_OK_pseudo
          else
                ierror = fit_OK
          end if
          ! Ignore error message because in this case we are unlikely to call it with the wrong input
          ! format, so if the covariances are wrong, the problem should have been caught before.
          call  var_auc_PBM(d_a_par, c_par, cov(1,1), cov(2,2), cov(1,2), variance_auc, ierror)
          ! Load reduced variance.
          cov_out(1:num_cat+1,1:num_cat+1) = cov(1:num_cat+1,1:num_cat+1)
     endif
     ! Load the reduced hessian into the hessian for output. Reduced refers to
     ! the fact that only the labroc4 corners are included. If there are extra places
     ! the values they were initialized with remains
     if (present(hessian_out)) hessian_out(1:num_cat+1,1:num_cat+1) = hessian(1:num_cat+1,1:num_cat+1)

     ! Added by LP at U of C August 2009 to prevent calling the reconstruction of cutoffs when it is not
     ! necessary. This steps recreates the cutoffs for the full non-collapsed data (i.e., it is possible
     ! that the input data was not reduced to the truth runs)
     if(num_cat /= num_categ) then
               call compute_orig_cutoffs(d_a_par, c_par, num_categ, num_cat, new_cat_index, vc_cutoffs, &
                                catn_in, cats_in, catn, cats, idebug, vc_cutoffs_out)
     else
               vc_cutoffs_out = vc_cutoffs
     endif

case(fit_perfect)
      ! There are infinite possible solutions which produce the exact same curve
      ! we return one (here we returned a right skewed version because it takes
      ! and infinite value for d_a to have an exact area of 1, while the
      ! largest value of c is exactly 1)
      d_a_par      = 0.0_double
      c_par        = 1.0_double !
      auc          = 1.0_double
      variance_auc = 0.0_double ! this is the result of the MLE approach, not that any rational
                                 ! person should believe it.
      ! the other quantities are meaningless or should be considered as such -- reset by LP October 2009
      vc_cutoffs_out  = -666.0_double
      log_like        = -666.0_double
      cov_out         = -666.0_double
      if (present(hessian_out))  hessian_out = -666.0_double
case(fit_failed, fit_undetermined, fit_fail_init_est)
! Mathematically it should be undetermined, but we put an impossible/absurd/crazy value
! to prevent interfacing problems from arising
      d_a_par         = -666.0_double
      c_par           = -666.0_double
      auc             = -666.0_double
      variance_auc    = -666.0_double
      vc_cutoffs_out  = -666.0_double
      log_like        = -666.0_double
      cov_out         = -666.0_double
      if (present(hessian_out))  hessian_out = -666.0_double
case default ! This happens if there is a programming error, obviously...
! Mathematically it should be undetermined, but we put an impossible/absurd/crazy value
! to prevent interfacing problems from arising
      d_a_par         = -666.0_double
      c_par           = -666.0_double
      auc             = -666.0_double
      variance_auc    = -666.0_double
      vc_cutoffs_out  = -666.0_double
      log_like        = -666.0_double
      cov_out         = -666.0_double
      if (present(hessian_out))  hessian_out = -666.0_double
end select


deallocate(catn, cats) ! Deallocate the arrays in the local data modules
deallocate(cov , hessian, vc_cutoffs, vc_cutoffs_array)

!---------------------------------------------------------------------------------------------------
end subroutine pbmroc_mle
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine compute_orig_cutoffs(d_a_par, c_par, num_cat_orig, num_cat, new_cat_index, vc_cutoffs, &
                                catn_orig, cats_orig, catn, cats, idebug, vc_cutoffs_orig)
!---------------------------------------------------------------------------------------------------
! PURPOSE: reconstruct cutoffs for an array that needed to be collapsed because was sent to the MLE routine
!          without being  collapsed.
!HISTORY: November 17 2009, LP added debug output, changed the calculation of intermediate FPF (was wrong) and
!                              added checking for many empty categories at the boundary,
!                              changed internal variables following scheme "" -> orig
use debugging
implicit none

real(kind = double),intent(in)                :: d_a_par, c_par ! MLE of the parameters
integer, intent(in):: num_cat_orig ! number of categories in non-collapsed matrix
integer, intent(in):: num_cat ! number of categories in collapsed matrix
integer, intent(in):: idebug ! whether to log information or not
integer, dimension(num_cat_orig), intent(in):: catn_orig, cats_orig ! arrays containing categorical data non-collapsed
integer, dimension(num_cat), intent(in):: catn, cats ! arrays containing categorical data collapsed
real(kind = double), dimension(num_cat-1), intent(in) :: vc_cutoffs  ! cutoff parameter values at the maximum found
real(kind = double), dimension(num_cat_orig-1), intent(out) :: vc_cutoffs_orig  ! cutoff parameter values at the maximum found
integer, dimension(num_cat_orig), intent(in):: new_cat_index ! if the data needs to be collapsed, the list of where to
                                                          ! put the old data

real(kind = double) :: cutoff, tpf , fpf, one_minus, prev_fpf, prev_tpf
real(kind = double), parameter:: tol =  1.0e-10_double ! tolerance in the computation of tpf and fpf


integer:: icat ! category (cutoff) index, for the non-collapsed
integer:: num_cases ! number of cases before this category in this truth run, used if collapsing is present
integer:: old_cat_coll ! the last value of the collapsed category index
character(len = line_length):: msg    ! character string used to display state information

! Write to file debug information -- Nov 2009
 if(idebug==1) then
      continue! write(msg,"(' ')")
      call DisplayState(msg)
      continue! write(msg,"('CUTOFFS ARE RECONSTRUCTED (EXACTLY) FROM THE TRUTH-RUNS')")
      call DisplayState(msg)
      continue! write(msg,"('Orig  c.psed   #orig+ /    #orig- /      FPF             TPF  ')")
      call DisplayState(msg)
      continue! write(msg,"('cat   cat      #c.psed+    #c.psed-  ')")
      call DisplayState(msg)
 endif

! set the number of cases inside of each of the
num_cases = 0
old_cat_coll = new_cat_index(1) ! added Nov 2009, to check for changing collapsed category
! Rebuild the cutoffs if needed.
non_collapsed_categories: do icat = 1, num_cat_orig - 1
       if(old_cat_coll /= new_cat_index(icat)) then ! Added Nov 2009
           old_cat_coll = new_cat_index(icat)
           num_cases = 0
       endif
      ! If the first category is empty, its cutoff corresponds to 1,1, that is the smallest value of vc
      if( icat == 1 .and.  catn_orig(icat) == 0 .and. cats_orig(icat) == 0) then
          ! if c < 0, then the lower cutoff is bounded and known, we just use a couple of max and min to
          ! prevent it from overflowing.
          if (c_par <= 0) then
              cutoff = &
                  d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * min(c_par,-.00001_double) )
          else ! by semplifying the FPF equation when FPF == 1, one can come up with this value
             cutoff =  - 10.0_double / ( 1 - sign( min(.999999999_double, abs(c_par)),c_par))
          endif
       ! If the last categories are empty, this cutoff corresponds to (0,0) changed to all in Nov 2009
      elseif( all(catn_orig(icat+1:num_cat_orig) == 0) .and. all(cats_orig(icat+1:num_cat_orig) == 0)  ) then
          ! by semplifying the TPF equation when TPF == 0, one can come up with this value
          if (c_par <= 0) then
             cutoff = 10.0_double / ( 1 + sign( min(.999999999_double, abs(c_par)),c_par)) + d_a_par*sqrt(1.0_double + c_par**2) / &
                     ( 2.0_double * max(-c_par,.00001_double) * (1.0_double + c_par**2) )
          else ! if c > 0, then the upper cutoff is bounded and known, we just use a couple of max and min to
          ! prevent it from overflowing.
             cutoff =  d_a_par *  SQRT ( 1.0_double + c_par**2) / ( 4.0_double * max(c_par,.00001_double) )
          endif
      else
             ! If both the positive cases and the negative cases are different from 0 for this category, then
             ! it can't be collapsed with any other category, so it has the cutoff of the category where it
             ! was collapsed to.
             if(catn_orig(icat) /= 0 .and. cats_orig(icat) /= 0) then
                  cutoff = vc_cutoffs(new_cat_index(icat) )
                  num_cases = 0 ! This is a new category, so we reset this counter to zero
             ! For cutoffs in the middle, if the category is empty, they get the same cutoff as the previous
             ! category, the extremes are dealth with at the beginning
             elseif(catn_orig(icat) == 0 .and. cats_orig(icat) == 0) then
                  cutoff = cutoff
                  num_cases = num_cases ! This is an empty cat, this counter is not affected
             elseif(cats_orig(icat) == 0) then  ! There are only negative cases
                  ! If collapsed and non collapsed have the same number of cases, then this category
                  ! wasn't collapsed to categories with number of cases > 0, and it takes the same cutoff
                  ! if other empty categories follow, they will also have the same cutoff
                  if  (catn_orig(icat) ==  catn(new_cat_index(icat) ) ) then
                       cutoff = vc_cutoffs(new_cat_index(icat) )
                       num_cases = 0 ! This is a new category, so we reset this counter to zero
                  else ! We are sure that this category is collapsed together with some other one
                      if( new_cat_index(icat) == 1) then ! this is the first cutoff and it is not empty. We know that the
                          ! point before this point is (1.0, 1.0) because this is the first category
                          prev_fpf = 1.0_double
                      else ! We need to compute the previous value, the array can't get messed up because
                          call fpf_PBM(d_a_par, c_par, vc_cutoffs(new_cat_index(icat)  - 1) , prev_fpf, one_minus)
                      endif
                          ! Fixed on August 2009 by LP at U of C. To prevent the search for a cutoff above vc_cutoff(numcat
                          ! First compute the theoretical operating point for the category to which this one was collapsed
                          ! need to make sure that the collapsed categories are not the last ones because otherwise there
                          ! no cutoff after this one and the next operating point is (0,0)
                          if( new_cat_index(icat) == num_cat) then
                                 fpf = 0.0_double
                          else
                              call fpf_PBM(d_a_par, c_par,  vc_cutoffs(new_cat_index(icat) ) , fpf, one_minus)
                          endif
                          num_cases = num_cases + catn_orig(icat) ! num_cases must start here, because it is the first
                          ! Compute the theoretical fpf for this point, by partitioning the probability between
                          ! this cutoff and the next one
                          fpf = prev_fpf - (prev_fpf - fpf) * num_cases / catn(new_cat_index(icat) )
                          call fpf_find_vc(d_a_par, c_par, fpf, tol, cutoff) ! find the cutoff
                  endif
             elseif(catn_orig(icat) == 0) then ! There are only positive cases
                  ! If collapsed and non collapsed have the same number of cases, then this category
                  ! wasn't collapsed to categories with number of cases > 0, and it takes the same cutoff
                  if  (cats_orig(icat) ==  cats(new_cat_index(icat) ) ) then
                       cutoff = vc_cutoffs(new_cat_index(icat) )
                       num_cases = 0 ! This is a new category, so we reset this counter to zero
                  else ! We are sure that this category is collapsed together with some other one
                      if( new_cat_index(icat) == 1) then ! this is the first cutoff and it is not empty. We know that the
                          ! point before this point is 1.0 1.0 because this is the first category
                           prev_tpf = 1.0_double
                      else ! We need to compute the previous value, the array can't get messed up because
                          call tpf_PBM(d_a_par, c_par,  vc_cutoffs(new_cat_index(icat)  - 1) , prev_tpf, one_minus)
                      endif
                          ! Fixed on August 2009 by LP at U of C. To prevent the search for a cutoff above vc_cutoff(numcat
                          ! First compute the theoretical operating point for the category to which this one was collapsed
                          ! need to make sure that the collapsed categories are not the last ones because otherwise there
                          ! no cutoff after this one and the next operating point is (0,0)
                          if( new_cat_index(icat) == num_cat) then
                                 tpf = 0.0_double
                          else
                                 call tpf_PBM(d_a_par, c_par,  vc_cutoffs(new_cat_index(icat)) , tpf, one_minus)
                          endif
                          num_cases = num_cases + cats_orig(icat) ! num_cases must start here, because it is the first
                          ! Compute the theoretical fpf for this point, by partitioning the probability between
                          ! this cutoff and the next one
                          tpf = prev_tpf - (prev_tpf - tpf) * num_cases / cats(new_cat_index(icat) )
                          call tpf_find_vc(d_a_par, c_par, tpf, tol, cutoff) ! find the cutoff
                  endif
             endif
       endif
       vc_cutoffs_orig(icat) = cutoff

 ! if needed report log data -- November 2009
 if(idebug==1) then
      call fpf_PBM(d_a_par, c_par,  vc_cutoffs_orig(icat) , fpf, one_minus)
      call tpf_PBM(d_a_par, c_par,  vc_cutoffs_orig(icat) , tpf, one_minus)
      continue! write(msg,"(1x,I3,4x,I3,4x,I3,'/',I3,4x,I3,'/',I3,4x,2(f14.10,2x) )") icat, new_cat_index(icat), &
      continue!           catn_orig(icat), catn(new_cat_index(icat)), cats_orig(icat), cats(new_cat_index(icat)), &
      continue!           tpf, fpf
      call DisplayState(msg)
 endif



enddo non_collapsed_categories



!---------------------------------------------------------------------------------------------------
 end subroutine compute_orig_cutoffs
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine compute_covariance_matrix(d_a_par, c_par, num_param,idebug, hessian,cov,ierror, err_msg)
!---------------------------------------------------------------------------------------------------
! PURPOSE:  COMPUTATION OF VARIANCES by  Computing the inverse of the hessian to obtain the variance covariance
!            matrix
! ALGORITHM: First the program attempts a straigth inversion, but sometimes the proper model hessian is
!            singular, because the parameter d_a becomes almost redundant (technically changing d_a and
!            translating all the cutoffs leaves the Log Likelihood numerically unchanged), since there is
!            a redundant direction, the hessian becomes singular and can't be inverted. In this case we remove the
!            redundant direction, but diagonalizing the matrix and eliminating the smallest eigenvalue, which is
!            the one associated with the almost flat direction in the likelihood space (so the second derivative
!            is zero and accordingly its eigenvalues is also zero or close to it).

use l_algebra, only: pseudoinverse, m1tom2, m2tom1, sinv
use debugging
implicit none

real(kind = double),intent(in)                        :: d_a_par, c_par ! MLE of the parameters
integer, intent(in):: num_param ! number of rows or columns in the hessian and covariance matrix
integer, intent(in):: idebug ! whether to log information or not

real(kind = double), dimension(num_param,num_param), intent(in)  :: hessian
real(kind = double), dimension(num_param,num_param), intent(out)  :: cov
integer, intent(out)                                   :: ierror ! Error flag about variance-covariance
                                                       ! matrix calculation, it works as follows
! ERROR FLAGS to return
 integer, parameter :: hessian_variance    = 0 ! variance computed using the inverse of the negative hessian
 integer, parameter :: pseudo_inv_variance = 2 ! variance computed using the Moore-Penrose pseudoinverse, specifically
          ! dropping the smallest eigenvalue - see header description of the procedure
 integer, parameter :: failed_variance     = 1 ! the variance could not be computed

character(len = line_length), intent(out):: err_msg    ! description of the error occurred

! internal variables
real(kind = double), dimension( (num_param)*(num_param) ):: xxdum
real(kind=double), parameter :: inversion_tolerance = 1.0e-3_double
integer:: inversion_error ! whether the inversion worked or not, index of the matrix being singular (see "ALGORITHM"
                         ! above
real(kind = double), dimension(num_param,num_param) :: identity ! matrix to check the quality of the inversion
integer:: i,j ! loop counters
character(len = line_length):: msg    ! messages to be reported is idebug is on


! Copy the negative of the upper triangular hessian in a packed  array
call m2tom1(num_param,- hessian,xxdum)
call sinv(xxdum,num_param,inversion_tolerance,inversion_error) ! invert the matrix (if not singular)
call m1tom2(num_param,xxdum,cov) ! Copy from packed upper triangular to full matrix

! If asked, report on the inversion of the full hessian
if(idebug == 1 .and. inversion_error /= 0) then
  if(inversion_error < 0) then
      continue! write(msg,*) "WARNING:: FAILED to invert Full Hessian to obtain covariance matrix"
  elseif(inversion_error > 0) then
      continue! write(msg,*)"WARNING:: unstable inversion of Hessian  SINV error flag:",inversion_error
  endif
  call DisplayState(msg)
endif

! If either inversion was unstable, or we would expect it to be unstable, compute and return the pseudoinverse
! to cut out the direction that does not affect the LL, and write that results on top of cov
! Note that .5 is sort of an empirical value, based upon the fact that there is hardly ever a fit in this area
! (say 1 time in 10,000,000). Note that pseudoinverse will check whether in fact it is a good idea to
! remove the last eigenvalue (using a threshold). Sometimes the hessian is kind of flat, but the inverse
! isn't really too bad, and it is better than the alternative of removing one piece.
if(inversion_error /=0 .or. d_a_par < .5_double* abs(c_par) ) then
   if(idebug == 1) call DisplayState(" VARIANCES computed using the pseudoinverse                   ")
   call  pseudoinverse(-hessian, num_param, idebug, cov, inversion_error, err_msg)
   ! Note that the value of the inversion error flag is not inquired for different source
   ! of failure. we leave this to the caller. The truth is that we never had this problems,
   ! currently they are only theoretically possible. If we will find situations where it
   ! does happen frequently enough, we might want to do something about it in this routine
   if(inversion_error /= 0) then
      ! return the same error message as returned here by pseudoinverse
      ierror  = failed_variance ! We could not compute variances
      if(idebug == 1) call DisplayState(err_msg)
   else
      continue! write(err_msg,*) "Variances computed with Pseudo-inverse of Hessian matrix"
      ierror  = pseudo_inv_variance ! we used a pseudoinverse, we should flag it
   endif
else ! The message is returned by proproc containing info about the computation of the variances
   continue! write(err_msg,*) "Variances computed inverting the negative Hessian"
   if(idebug == 1) call DisplayState(err_msg)
   ierror = hessian_variance ! Simply inverting the - Hessian produced numerically stable variances
endif

! Now use whatever cov we decided to use to
! Compute Hessian * Hessian^-1 to see whether we get the identity matrix
! Of course, we can't do that if the computation of the variance failed, so we have to check that
if(idebug == 1 .and. (ierror == hessian_variance .or. ierror == pseudo_inv_variance) ) then
    identity = - matmul(cov,hessian)
    call DisplayState("   COVARIANCE MATRIX TIMES HESSIAN                   ")
    if(ierror == hessian_variance) then
       call DisplayState("     It should be identity")
    else
       call DisplayState("     Deviation from identity is a measure of the size of the coordinate removed")
       call DisplayState("     Elements linked to c (parameter nr 2) should be close to identity matrix")
    endif
    do i=1, num_param
      do j=1, num_param
         if( abs(identity(i,j)) > 1.0d-10)  then
             continue! write(msg,"('    ', 2(1x,i6),1x, g22.15)")i,j, identity(i,j)
             call DisplayState(msg)
         endif
      enddo
    enddo
endif

!---------------------------------------------------------------------------------------------------
end subroutine compute_covariance_matrix
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine initialize_d_a_c(mn, ms, num_cat, catn, cats, idebug, &
                       d_a_par_array, c_par_array, ierror, err_msg)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
! PURPOSE: initialize the  values of d_a and c, the 2 parameters of the proper binormal model.
!          Moreover checks whether the data received is compatible with the model or if it
!          has analytic solutions (in d_a, c)
!ALGORITHM:  call least squares functions to determine the value of d_a and c, while it
!           uses some simple rules to determine whether the datasets are suitable for the MLE
!           It provides 2 initial estimates, one obtained with transforming a conventional fit,
!           and the other using a skewed fit, skewed the same way as the conventional fit was
!           This is done to prevent the algorithm to get stuck in the wrong maximum if multiple
!           maxima are available
! CHANGES: Feb 2007 - removed a bug in the definition of the sign of c (see in the body)
!          LP UC 28/12/2009 changed the algorithm for degenerate data, to harmonize with labroc4 and to prevent
!                           perverse fits from being considered perfect.
use debugging ! module containing the procedures to log the debugging information on the run
! LP - Feb 2008 introduce the module roc_nonparametric
use roc_nonparametric, only: empirical_operating_points_cat
use labroc_functions, only: least_squares_a_b
use statistic_functions, only: phi
implicit none

! Subroutine arguments
integer, intent(in):: mn ! number of actually negative cases
integer, intent(in):: ms ! number of actually positive cases
integer, intent(in) :: num_cat ! Number of categories as created by catgrz
integer, dimension(num_cat), intent(in):: catn, cats ! arrays containing categorical data
integer, intent(in) :: idebug               ! 0 = no debug; 1  = debug


integer, intent(out)                                   :: ierror ! Error flag about the MLE fit
character(len = line_length), intent(out):: err_msg    ! description of the error occurred, follows the numeration
          ! for the return values from proproc_<vers>, lower routines follows their own error numbering

character(len = line_length) :: msg    ! text for output to out files
real(kind = double), dimension(3), intent(out)  :: d_a_par_array, c_par_array ! initial estimates of the params
                 ! 2 values to diminish risk of multiple maxima


real(kind = double)  :: d_a_par, c_par ! Initial estimates, internal
real(kind = double)  :: d_a_par_0, c_par_0 ! Initial estimates from the unconstrained conventional fit
real(kind=double) :: cumul_fraction(act_neg:act_pos,num_cat+1)! NOTE that the
           ! cumulative fraction that corresponds to category i has index num_cat - i+1
           ! because categories start at + infinity
real(kind = double)   :: a_par, b_par, a_z ! contentional binormal model parameters, used in initial estimation
real(kind = double)   :: auc ! value of the area under the curve for the proper model
integer::i ! loop counter
integer:: skew_pos ! whether the initial estimates have to be taken skewed positive or skewed negative
                   ! note that it is only used if we don't keep the conventional fit
! Compute the positive and negative fractions corresponding to the categories
! found and print them to the log file, if asked to do so

call empirical_operating_points_cat( num_cat, mn , catn, cumul_fraction(act_neg,1:num_cat+1))
call empirical_operating_points_cat( num_cat, ms , cats, cumul_fraction(act_pos,1:num_cat+1))

if(idebug == 1) then
      continue! write (msg,*) "NUMBER OF CATEGORIES ", num_cat
      call DisplayState(msg)
      do i = 1, num_cat + 1
         continue! write (msg,*) cumul_fraction(act_neg,i), cumul_fraction(act_pos,i)
         call DisplayState(msg)
      enddo

      continue! write (msg,*) "CATEGORICAL DATA USED IN THE MLE, total cat = ", num_cat
      call DisplayState(msg)
      do i = 1, num_cat
         continue! write (msg,*) catn(i), cats(i)
         call DisplayState(msg)
      enddo
      continue! write (msg,*) "END OF CATEGORICAL DATA"
      call DisplayState(msg)

endif

!   CHECK WHETHER THERE ISN'T ANY POINT INSIDE THE ROC PLOT, I.E., THAT HAS BOTH
!   0 << FPF AND TPF << 1,  THE CURVE IS DEGENERATE WITH PERFECT SEPARATION for the
!   proper binormal model
call check_degen_PBM(num_cat, cumul_fraction,ierror,err_msg)


! Capure standard exceptions
select case (ierror)
! These are all fits that render looking unecessary/useless to attempt
! to fit the data. In this case the value of the parameters provided by
! the function that checkes degeneracy is the fit we plan to return
case(fit_perfect, fit_undetermined, bad_input)
     return
case default

end select

! Compute approximate a & b parameters of the standard bivariate binormal model using a least
! squares regression on the TPF and FPF - to compute initial estimates

! START WITH AN UNCONSTRAINED CONVENTIONAL FIT

call least_squares_a_b(cumul_fraction(act_neg,:), cumul_fraction(act_pos,:), &
                        num_cat+1, mn, ms, a_par,b_par, 0,ierror)

if (ierror /= 0 ) then
     ierror = fit_fail_init_est ! failed fit message because it couldn't make it
     err_msg = "failed the initial estimates of a and b to obtain d_a and c "
     return
endif

! "Convert" a and b in d_a and c (the conversion is not entirely possible because the 2 systems are
! not equivalent. Make sure that if d_a and c are both small the transformation will not end up
! in a very skewed curve with area of 1
d_a_par = sqrt(2.0_double) * a_par/ sqrt(1.0_double + b_par**2)

! We need the value of a_z to decide how to deal with small c's
a_z =  phi(a_par/sqrt(1.0_double + b_par**2))

if( abs(a_z - .5) > .0001_double) then
     ! correct c from getting too close to 0 when the area is small (for numerical stability)
     c_par   = sign( min(abs(( b_par - 1.0_double )  /  ( b_par + 1.0_double )),.99_double), b_par - 1.0_double)
else
     ! Feb 2007 - changed from sign(...,c_par) because c_par in the new formulation is undefined at this stage
     c_par = sign(.0001_double, 1.0_double - b_par)
endif

! Write out details of fit is required to do so
if(idebug == 1) then
   call DisplayStateMsg('CONVENTIONAL UNCONSTR. L-SQUARES PARAM. INITIAL ESTIMATE')
   continue! write (msg,"( 7x,'A',13x,'B',13x,'AZ')")
   call DisplayStateMsg(msg)
   continue! write(msg,"( 3(2x,f12.6))") a_par,b_par, a_z
   call DisplayStateMsg(msg)

   call auc_PBM(d_a_par, c_par, auc, ierror)

   ierror = fit_OK ! We don't care if the computation of the area had some precision issues

   continue! write (msg,"( 7x,'D_a',13x,'C',13x,'AUC')")
   call DisplayStateMsg(msg)
   continue! write(msg,"( 3(2x,f12.6))") d_a_par, c_par, auc
   call DisplayStateMsg(msg)

 endif

 ! Make sure that the parameters aren't in a numerically possibly unstable region
 d_a_par_0 = d_a_par
 c_par_0 = c_par
 d_a_par_array(1) =   max(d_a_par_0,.01_double)
 c_par_array(1)   = sign(    max( abs(c_par_0), c_almost_zero*100000 ),   c_par_0 )



! D == 0 Least squares fits (one + c and one - c). Note that occasionally the sign of c
! will be wrong from the initial estimates, so to be sure , we shoot 2 estimates
 do i = 1, 2
     if( i  ==  1) then
          skew_pos = 1
     else
          skew_pos = 0
     endif

     ! Make the fit for the desired symmetry
     call least_squares_d_a_zero(skew_pos,cumul_fraction(act_neg,:), cumul_fraction(act_pos,:), &
                                 num_cat+1, mn, ms, c_par)
     ! correct c from getting too close to 1, and d_a from getting too close to 0
     ! d_a is corrected like this to prevent messing things up when the area is too close to .5
     c_par   = sign(  min( abs(c_par), .99_double)   , c_par)
     c_par   =  sign(     max(   abs(c_par), c_almost_zero*100000  ),   c_par    )
     d_a_par  = .3_double * abs(c_par) ! This is a simple solution to prevent d_a to both become too
                ! small and to become too large compared to d_a when the curve is too close to
                ! the chance diagonal (AUC = .5 ; d_a = 0, c = 0). Note that when d_a is smaller than
                ! about .5 c, it hardly affects curve shapes and likelihoods at all, so it doesn't really
                ! matter what value it is used.

     if(idebug ==1 ) then
        call DisplayStateMsg('ASYMMETRIC LEAST SQUARES INITIAL ESTIMATES')

        call auc_PBM(d_a_par, c_par, auc, ierror)
        ierror = fit_OK ! We don't care if the computation of the area had some precision issues

        continue! write (msg,"( 7x,'D_a',13x,'C',13x,'AUC')")
        call DisplayStateMsg(msg)
        continue! write(msg,"( 3(2x,f12.6))") d_a_par,c_par, auc
        call DisplayStateMsg(msg)
     endif

     ! load the initial estimate
     d_a_par_array(i+1) = d_a_par
     c_par_array(i+1) = c_par

enddo

!---------------------------------------------------------------------------------------------------
end subroutine initialize_d_a_c
!---------------------------------------------------------------------------------------------------

!------------------------------------------------------------------
!------------------------------------------------------------------
 Subroutine least_squares_d_a_zero(c_positive,FPF,TPF, num_pts, mn,ms, c_par)
!------------------------------------------------------------------
! Use least squares on the TPF and FPF to make an estimate of the
! parameters a, b for a binormal ROC model, using linear least
! squares over the deviates.
! NOTE: this is the same as doing a least squares fitting a & b with
! the constraint a = 0
! Warning:: NEEDS TO BE CLEANED and get some error messaging!
! CHANGES: LP June 2007 made correction smaller and symmetrized for mn and ms to
!          reduce the differences between symmetrical datasets (+c, -c) when they
!          are very small.

use statistic_functions, ONLY : zdev
implicit none

integer, intent(in):: c_positive !whether c is consider close to +1 (==1) or
                     ! to -1 ( == 0)
integer, intent(IN):: num_pts ! number of points, including 0,0 and 1,1
real(kind=double), dimension(num_pts+1), intent(IN):: FPF ! because it contains also 0,0 and 1,1
real(kind=double), dimension(num_pts+1), intent(IN):: TPF ! because it contains also 0,0 and 1,1
integer, intent(IN):: mn,ms ! # normal and abnormal cases
real(kind=double), intent(out):: c_par ! parameter c of the proper binormal model


real(kind=double), dimension(num_pts-2):: x, xs
real(kind=double):: sumxy, sumx2, slope ! Terms used by the least squares estimate
real(kind=double):: start, unit ! the min &  max values that the x variable can have in the model
                    ! which is different depending upon the sign of c (2 for positive
                    ! (one for negative)
integer i

! Fix the fit accordingly to the sign of c, the equations are
! different, do not copy (1,1) which goes to infinity, but keep
! (0,0) as a constraint (not sure about this one)
! The two estremal points (0,0) and (1,1) have to be eliminated
! from the fit

   if(c_positive==1) then
      x(1:num_pts-2) = fpf(2:num_pts-1) + 1.0_double
      xs(1:num_pts-2) = tpf(2:num_pts-1) + 1.0_double
      start = 1.0_double
      unit = 2.0_double
   else
      x(1:num_pts-2) = fpf(2:num_pts-1)
      xs(1:num_pts-2) = tpf(2:num_pts-1)
      start = 0.0_double
      unit = 1.0_double
   endif


!
!     TEST FOR 1.'S AND CORRECT BY SUBTRACTING 1/2N. TEST FOR 0'S
!     AND CORRECT BY CHANGING TO 1/2N.  THEN CALL SUBRoutINE ZDEV TO
!     TRANSFORM THE CUMULATIVE PROBABILITIES IN THE ARRAYS TO STANDARD
!     NORMAL DEVIATES.
  do i = 1 , num_pts - 2
         IF(x(i) .speq. start)    x(i) = start + 1.0_double/(10.0_double*real(mn+ms))
         IF(x(i) .speq. unit)    x(i)  = unit - 1.0_double/(10.0_double*real(mn+ms))

         IF(xs(i) .speq. start)  xs(i) = start + 1.0_double/(10.0_double*real(mn+ms))
         IF(xs(i) .speq. unit ) xs(i)  = unit - 1.0_double/(10.0_double*real(mn+ms))
  enddo


!     ADJUST THE SEQUENCE OF CUTOFFS TO AVOID IDENTICAL VALUES
if(c_positive==1) then
  do i = 2, num_pts - 2
    IF( x(i)  <=  x(i-1)  )    x(i)  = max( x(i), x(i-1)  ) + 1.0_double/(100.0_double* real(mn+ms))
    IF( xs(i) <= xs(i-1) )     xs(i) = max( xs(i), xs(i-1) ) + 1.0_double/(100.0_double* real(ms+ms))
  enddo
else
  do i = num_pts - 2, 2
    IF( x(i)  <=  x(i-1)  )    x(i-1)  = min( x(i), x(i-1)  )  - 1.0_double/(100.0_double* real(mn+ms))
    IF( xs(i) <= xs(i-1) )     xs(i-1) = min( xs(i), xs(i-1) ) - 1.0_double/(100.0_double* real(ms+ms))
  enddo
endif

! Compute the deviates for the proper fit for d_a equal zero
   x(1: num_pts - 2) =  x(1:num_pts - 2) / 2.0_double
  xs(1: num_pts - 2) = xs(1:num_pts - 2) / 2.0_double

   x(1: num_pts - 2) = zdev(  x(1:num_pts - 2) )
  xs(1: num_pts - 2) = zdev( xs(1:num_pts - 2) )

!
! CALCULATE LEAST SQUARES SOLUTIONS
!

  sumxy = dot_product(  x(1: num_pts - 2), xs(1: num_pts - 2)  )
  sumx2 = dot_product(  x(1: num_pts - 2),  x(1: num_pts - 2)  )

  slope  = sumxy/sumx2

! Determine the parameter c (since d_a is set to ~ 0 for this fit as an assumption)
! Note that to correct for pathological datasets one has to enforce the correct sign
! to the parameter c (sometimes too many points are below the guessing line when modeling
! c positive, for example)

  if (c_positive==1) then
        c_par  =  max(c_almost_zero, min(c_almost_one,(slope - 1.0_double) / ( slope + 1.0_double)) )
  else
        c_par  =  min(-c_almost_zero,max(-c_almost_one,(slope - 1.0_double) / ( slope + 1.0_double)) )
  endif

!------------------------------------------------------------------
 end subroutine least_squares_d_a_zero
!------------------------------------------------------------------
!------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine check_degen_PBM(num_cat,cumul_fraction,ierror,err_msg)
! PURPOSE: determine if data are perfectly  separated, i.e. all positives are above all negatives
!          or whether the proper fit of the data will produce a perfect fit because this last would
!          be the analytical solution, so it is useless to run an MLE. This is the only situation, barred
!          a sequence of points on the FPF = TPF axis that would produce a perfect fit. However, only
!          the former might be problematic and it is considered "degenerate" and asymptotic.
! NOTE:    We do not return ierror as degenerate as there is only one type and perfect fit seems a more
!          appropriate definition (LP October 2009, U of C)
!
!         LP UC 12/29/2009 changed the selection making it dependend upon number of cats, to harmonize with cvbm
implicit none

integer, intent(in)                      :: num_cat ! number of categories
real(kind=double), intent(in)            :: cumul_fraction(act_neg:act_pos,num_cat+1)! operating points
integer, intent(out)                     :: ierror ! Error flag with information about the MLE fit
character(len = line_length), intent(out):: err_msg    ! description of the error occurred

real(kind=double), parameter :: zero_frac = 1.0e-8_double !The cutoff value under which we consider the difference between

intrinsic all

select case (num_cat)
case(4:)
     ierror = fit_OK
case(3) ! three categories => two operating points.
     if( cumul_fraction(act_neg,2)<zero_frac .and. cumul_fraction(act_pos,num_cat)>1.0_double-zero_frac) then
           ierror = fit_perfect
           err_msg = "The data is degenerate and can be fit exactly by a perfect curve (AUC = 1)"
     else !The dataset is "normal" and therefore the variable values are already set
           ierror = fit_OK
     endif
case (2)!Number of categories is insufficient for the fit and extremely large number of CvBM curves can be fit to this type of data
     !the return values are undefined as b could have values from 0 to +oo and a from  zdev(cumul_fraction(act_pos, 2) ) to +oo
     ! (or something like that, I am not going to find out what pathological entities are compatible with this).
     ! of course the data could be perfectly separated or perverse
     !Check whether a perfect fit is possible: first point has FPF = 0 and second point has TPF = 1.
     if( cumul_fraction(act_neg,2)<zero_frac .and. cumul_fraction(act_pos,num_cat)>1.0_double-zero_frac) then
           ierror = fit_perfect
           err_msg = "The data is degenerate and can be fit exactly by a perfect curve (AUC = 1)"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
     !Check whether the  fit is undetermined, which requires a single point and in the upper triangle of the ROC plot
     elseif( cumul_fraction(act_pos,2) > cumul_fraction(act_neg,2)) then
           ierror = fit_undetermined
           err_msg = "The data is degenerate and can be fit exactly by a perverse curve (AUC = 0)"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
     ! the point is below the chance line, since the fit is proper, it will produce the chance line.
     else
         ierror = fit_OK
         err_msg = " fit should be chance line"
     endif
case(:1)
        ! Less than two categories cannot represent a detection problem and therefore the data is not just insufficient, it is wrong.
       err_msg = "ROC data need to have al least two categories"
       ierror = bad_input
end select

!-----------------------------------------------------
end subroutine check_degen_PBM
!-----------------------------------------------------


!----------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine ddlike_PBM(num_param, param_vec, nf, grad, h, uip, urp, lk)
!--------------------------------------------------------------------
! PURPOSE: compute the derivatives  (gradient and hessian) of minus log likelihood function for the
!          proper binormal model

use ll_derivatives! data module with the values and first  derivatives of p_i's and q_i's as a function
                   ! of the parameters
use problem_data ! data module were problem data constants are stored
implicit none

integer, intent(IN):: num_param ! number of parameters
integer, intent(INout):: nf ! flag for the optmizer, it also contains error (if ==0)
integer, dimension(1), intent(INout) :: uip ! unused

real(kind=double), dimension(num_param), intent(IN):: param_vec ! array with the parameter values (d_a, c, cutoffs)
real(kind=double), dimension(num_param), optional, intent(out):: urp ! not used, it is vector used by the optimizer
real(kind=double), dimension(num_param), intent(out):: grad ! gradient
real(kind=double), dimension(num_param*(num_param+1)/2), intent(out):: h ! hessian,triangular

! internal variables
real(kind=double), dimension(num_param,num_param) :: full_h ! hessian, full
integer:: i,j, k ! loop counters

! lk, the last parameter is taken for compatibility with the optimizer used, but it isn't used
! in the internal calls, so it is optional (the same applies to the following routines
optional lk
external lk

! Call the function that constructs the gradient, this function will also make sure
! that the matrix of the drivatives of dp in d{parameter) is updated to the current
! iteration
call dlike_PBM(num_param, param_vec, nf, grad, uip)

! Note that the forall statement does not use the fact that a lot of the matrix elements are zero
! this moderately affects speed (since the really expensive part is to compute these elements, not
! to multiply them, but heavily affects the readibility of the code. We might change it in the future
 do i = 1 , num_param
    do j = 1, i
         full_h(i,j) = &
                 - sum( dp_dpar(1:num_param-1,i) * dp_dpar(1:num_param-1,j) / p_val(1:num_param-1) ) * num_normal_cases &
                 - sum( dq_dpar(1:num_param-1,i) * dq_dpar(1:num_param-1,j) / q_val(1:num_param-1) ) * num_abnormal_cases
         full_h(j,i) = full_h(i,j)
  enddo
enddo

! Load full matrix in diagonal form, as used by the TOMS611 optimizer
h(1) = -full_h(1,1)
k = 2
do i = 2, num_param
      h(k: k + i - 1) = -full_h(i,1:i)
      k = k + i
enddo

!--------------------------------------------------------------------------------------------
end subroutine ddlike_PBM
!---------------------------------------------------------------------------------------------

!-------------------------------------------------
subroutine dlike_PBM(num_param,param_vec,nf,grad, uip,urp,lk)
!------------------------------------------------
! PURPOSE:  Gradient of minus log likelihood function for the proper binormal model
!           variable and array shapes are done to match the interface with the trust zone Newton-Raphson
!           algorithm TOMS611

use problem_data ! data module were problem data constants are stored
use statistic_functions, only: phi, g ! normal distribution and density
use ll_derivatives ! data module with the values and first  derivatives of p_i's and q_i's as a function
                   ! of the parameters

implicit none

integer, intent(IN):: num_param ! number of parameters used in optimization  cutoffs + curve parameters
integer, intent(inout):: nf  ! variable used by optimizer, if it is set to zero it means wrong step
integer, dimension(1), intent(inout) :: uip ! non used vector

real(kind=double), dimension(num_param), intent(IN):: param_vec ! vector with the parameters for this step
                                                                ! first d_a, then c, then the cutoffs
real(kind=double), dimension(num_param), optional, intent(inout):: urp ! non used
real(kind=double), dimension(num_param), intent(out):: grad ! gradient of the log likelihood function

! internal variables.
real(kind=double):: log_like_val ! Log likelihood value

! variables related to the equations for the FPF, the TPF and their derivatives
real(kind=double):: af1,at1,af2,at2
real(kind=double):: df1,dt1,df2,dt2

! model parameters
real(kind=double):: d_a_par
real(kind=double):: c_par

real(kind=double):: d1_d_a, d2_d_a ! Derivatives of the arguments of the phi functions in d_a
real(kind=double):: d1_c, d2_c     ! Derivatives of the arguments of the phi functions in c
real(kind=double):: temp

real(kind=double) :: dfpf_dd_a, dfpf_dc, dfpf_dvc ! derivatives of the fpf as a function of the parameters
real(kind=double) :: dtpf_dd_a, dtpf_dc, dtpf_dvc ! derivatives of the tpf as a function of the parameters

real(kind=double), parameter:: min_p = 1.0e-10_double ! minimum cateogry probability, to be used only for
                                          ! the gradient calculation to prevent overflow

integer:: icat ! loop counter over categories

optional lk
external lk

character(len=line_length) ::  msg !  used to log warnings if required to do so


! If the function value  was just computed for this step, we can use the values needed by the gradient
! that were stored by the call to like_PBM, if not, we need to recompute it to get a new gradient
! - Not a large speed impact.
! Note that this makes sense only if these functions are called by the algorithm TOMS611, otherwise this
! function needs to call like_PBM all the time to make sure that the p_val and q_val (with the FPF and TPF values
! for the current cutoffs and parameter values) arrays are properly computed
if (uip(1) /=  nf) call like_PBM(num_param,param_vec,nf,log_like_val,uip)

! Load  the curve parameter values from the vector
d_a_par = param_vec(1)
c_par   = param_vec(2)

! Make sure that there is no potential overflow in the p values, this would make the
! gradient explode and kill the steps in the optimizer
p_val = max(min_p, p_val)
q_val = max(min_p, q_val)



! Build the derivatives of the likelihood function

! derivative of the 1st phi arguments (both for FPF and for TPF) in d_a_par (apart from sign)
d1_d_a = sqrt(1.0_double + c_par**2) / 2.0_double
! derivative of the second phi arguments (both for FPF and for TPF) in d_a_par (apart from sign)
if( abs(c_par)  < c_almost_zero) then
       d2_d_a = 0.0_double ! to prevent overflow, all these terms would cancel anyway
else
       d2_d_a = d1_d_a / c_par
endif

! Part of the derivatives of the fpf/tpf phi arguments in d(c) that doesn't depend upon vc
d1_c   =  ( d_a_par * c_par ) / ( 2.0_double * sqrt( 1.0_double + c_par**2) )
if( ABS(c_par)  < c_almost_zero) then
     d2_c = 0.0_double ! to prevent overflow, all these terms would cancel anyway
else
     d2_c   =  - ( d_a_par * SQRT( 1.0_double + c_par**2) ) / ( 2.0_double * c_par**2) &
               + d_a_par / ( 2.0_double * sqrt( 1.0_double + c_par**2) )
endif

! initialize the array with the gradient, and the derivatives of the p_i's and q_i's, that will be used to
! compute the hessian (if needed), and for the gradient (after this).
grad    = 0.0_double
dp_dpar = 0.0_double
dq_dpar = 0.0_double

!Build the terms that depend upon the cutoffs
do icat = 1, num_cat - 1
    ! Derivatives with respect to the  argument of the phi functions. They will be used to compute
    ! the derivatives as a function of the parameters =>
    !    d(p)/d(par) = d(phi1)/d(arg1) * d(arg1) / d(par) +  d(phi2)/d(arg2) * d(arg2) / d(par)

    ! compute the argument values for the first phi functions of the FPF and the TPF
    af1 =  z1_fpf(d_a_par, c_par, param_vec(2+icat) )
    at1 =  z1_tpf(d_a_par, c_par, param_vec(2+icat) )
    ! Since they are phi functions, the derivative is the gaussian of their argument
    df1 = g( af1)
    dt1 = g( at1)

    ! Compute the derivatives of the second phi functions, remember that in the proper model these might be
    ! zero if the value of c is zero (or numerically too close to zero)
    if( ABS(c_par)  < c_almost_zero) then
        df2 = 0.0_double
        dt2 = 0.0_double
    else ! same comments as for the first 2 phi functions
       af2 =  z2_fpf(d_a_par, c_par, param_vec(2+icat) )
       at2 =  z2_tpf(d_a_par, c_par, param_vec(2+icat) )
       df2 = g( af2)
       dt2 = g( at2)
    endif

   ! Add the first derivative term that depends upon the positive cases, note that the values of p_val
   ! were computed by the function like_PBM
   temp =  - catn(icat) /  p_val(icat)  +  catn(icat+1) / p_val(icat+1)

   dfpf_dd_a =    df1 * (-d1_d_a) + df2 * d2_d_a
   dfpf_dc   =    df1 * ( param_vec(2+icat) - d1_c ) + df2 * (param_vec(2+icat) + d2_c)
   dfpf_dvc  =  - (1.0_double - c_par) * (df1 + df2)
   ! Note that the derivatives in the curve parameters sum over all the cutoffs, while the cutoff derivatives,
   ! of course, regard only the specific cutoff whose derivative is being computed.
   grad(1)      =  grad(1) + temp * dfpf_dd_a
   grad(2)      =  grad(2) + temp * dfpf_dc
   grad(icat+2) =          + temp * dfpf_dvc  ! each cutoff is separated from the other ones

   ! Save the derivatives of dp_i in d{theta_j}
!   if(icat < 1) then
       dp_dpar(icat+1, 1)      =  dp_dpar(icat+1,1)      + dfpf_dd_a
       dp_dpar(icat+1, 2)      =  dp_dpar(icat+1,2)      + dfpf_dc
       dp_dpar(icat+1, 2+icat) =  dp_dpar(icat+1,2+icat) + dfpf_dvc
!   endif
   dp_dpar(icat,1)      = dp_dpar(icat,1)  -  dfpf_dd_a
   dp_dpar(icat,2)      = dp_dpar(icat,2) -  dfpf_dc
   dp_dpar(icat,2+icat) = dp_dpar(icat,2+icat) - dfpf_dvc

   ! Add the first derivative term that depends upon the negative cases
   ! the minus of the term from p cancels with the minus of the derivative of
   ! the argument in d{param_vec(2+i)} = - ( 1 + c)
   temp = - cats(icat)   / q_val(icat)  + cats(icat+1) / q_val(icat+1)

   dtpf_dd_a = + (  dt1 * d1_d_a + dt2 * d2_d_a )
   dtpf_dc   = + (  dt1 * ( - param_vec(2+icat) + d1_c ) + dt2 * ( - param_vec(2+icat) + d2_c) )
   dtpf_dvc  = - (1.0_double + c_par) * (dt1 + dt2)


   grad(1)      = grad(1) + temp * dtpf_dd_a
   grad(2)      = grad(2) + temp * dtpf_dc
   grad(icat+2) = grad(icat+2)  + temp * dtpf_dvc

   ! Load the derivatives for the P_i's and q_i's, that are needed to compute the expected hessian
   dq_dpar(icat+1,1)      =  dq_dpar(icat+1,1)      + dtpf_dd_a
   dq_dpar(icat+1,2)      =  dq_dpar(icat+1,2)      +  dtpf_dc
   dq_dpar(icat+1,2+icat) =  dq_dpar(icat+1,2+icat) + dtpf_dvc

   dq_dpar(icat,1)      = dq_dpar(icat,1) - dtpf_dd_a
   dq_dpar(icat,2)      = dq_dpar(icat,2) - dtpf_dc
   dq_dpar(icat,2+icat) = dq_dpar(icat,2+icat) - dtpf_dvc


! If the gradient is meaningless,  return an error message. In this case we will make TOMS611 reject the
! current step
! underflow part (if present in check_number) is necessary or useful at all, work in progress)
   if(.NOT. check_number( grad(icat+2)) .AND. (grad(icat+2) .spne. 0.0_double) ) then
     if(pass_idebug ==1) then
        call DisplayState("ERROR DLIKE_PBM:: OVERFLOW IN DERIVATIVES OF CUTOFFS ")
        continue! write(msg,*)  icat, param_vec(2+icat), grad(icat+2)
        call DisplayState(msg)
        continue! write(msg,*)   cats(icat), cats(icat+1), catn(icat),catn(icat+1)
        call DisplayState(msg)
     endif

     log_like_val   =  huge(param_vec(2+icat))
     nf = 0
     return
   endif


enddo

! take the minus the gradient, since this is what the optimizer is expecting
grad = - grad


!------------------------------------------------------------------------------
end subroutine dlike_PBM
!-------------------------------------------------------------------------------


!------------------------------------------------
!------------------------------------------------
subroutine like_PBM(num_param,param_vec,nf,log_like,uip,urp,unusedvariable)
!------------------------------------------------
! PURPOSE: compute  - LIKELIHOOD FUNCTION  to use to find the cutoffs that
!          maximize the log likelihood function given the curve parameters
!          d_a & c
!ALGORITHM: rewritten entirely by LP in summer 2004 to make use of the
!          symmetry in the proper function and in the cumulative distribution
!          function. This in turn allows to remove the setting of boundaries for
!          the maximum and minimum values of v_c (apart from the mathematical
!          ones)

use statistic_functions, only: phi,g
use problem_data
use ll_derivatives

implicit none

integer, intent(IN):: num_param ! number of parameters (num_cat + 1)
real(kind=double), dimension(num_param), intent(IN):: param_vec ! array with the parameters
     ! in this case (1) is d_a_par, (2) is c_par, the rest are the cutoffs
integer, intent(INOUT):: nf ! parameter for the thermonuclear optimizer, if set to zero means failure
real(kind=double), intent(out):: log_like ! function value
integer, dimension(1), intent(INout) :: uip ! parameters for the TMO, non used
real(kind=double), dimension(num_param), optional,intent(out):: urp
! This declaration has something to do with the Quasi Newton optimizer
! we have to deal with it in a more decent form, once we decided what to
! do with it
optional unusedvariable
external unusedvariable

! curve parameters
real(kind=double):: d_a_par
real(kind=double):: c_par

real(kind=double) :: vc_bound ! finite boundary for vc (see Metz & Pan)
                                ! sometimes an upper sometimes a lower
character(len=line_length) :: err_msg ! error message from subroutines, Not used by this program
integer:: icat ! loop counter over caregories

! values of FPFs and TPFs and their complements to 1, to be used when computing the log likelihood (to compute
! the p's and the q's
real(kind=double):: fpf, one_minus_fpf, last_fpf, last_one_minus_fpf
real(kind=double):: tpf, one_minus_tpf, last_tpf, last_one_minus_tpf
real(kind=double):: min_prob ! minimum accepted value for a probability

character(len=line_length) ::  msg !  used to log warnings

! load the values of the curve parameters
d_a_par = param_vec(1)
c_par   = param_vec(2)

min_prob = 10*tiny(1.0_double) ! better a parameter, but f95 does not allow it, minimal value of probability, values
    ! smaller that this one are rounded up to this value
uip(1) = nf ! variable for the TOMS611

! Check the parameter values, if they are in the mathematically meaningful range
! currently it is a very simple check, but is used to be complex, so we left it into a separate function, so that
! it can made complex again if needed.
if( .not. test_parameters(d_a_par, c_par) ) then
      if(pass_idebug  == 1) then
          continue! write(msg,"('d_a & c unacceptable:',2(2x,f14.8) )") d_a_par, c_par    ; call DisplayState(msg)
      endif
      log_like   =  huge( d_a_par )
      nf = 0
      return
endif

! the finite bound of the cutoffs values.
vc_bound = d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * sign( max(.00000000001_double, abs(c_par)),c_par))

! Check if the cutoffs passed by the calling program in the current step
! are mathematically meaningful
if(.not.check_cutoffs(err_msg)) then
      if(pass_idebug  == 1) call DisplayState(err_msg)
      log_like   =  huge(param_vec(1))
      nf = 0
      return
endif

! Initialize variables fpf and tpf values, since we start with
! category 1, its lower cutoff is the one that includes all values
! so the previous one is the smallest possible value
last_fpf =  1.0_double
last_tpf =  1.0_double
last_one_minus_fpf = 0.0_double
last_one_minus_tpf = 0.0_double

construct_p_i_and_q_i: do icat = 1, num_cat-1
       ! Evaluate the p (and 1-p) , q_j (and 1-q). The 1-p (or q) are needed
       ! for numerical stability
       call fpf_PBM(d_a_par, c_par, param_vec(2+icat),fpf, one_minus_fpf)
       call tpf_PBM(d_a_par, c_par, param_vec(2+icat),tpf, one_minus_tpf)
        ! construct the differences the p(i) = TPF(param_vec(2+icat-1)) - TPF(param_vec(2+icat))
        ! if numerical is lost and so the 2 p's or q's are identical, first
        ! try using the complement, if that does not work, set them to the
        ! numerical basis minimum value. Note that this does not take into
        !account situations where the 2 cutoffs are so close that their
        ! p values are almost identical and accordingly leave no numerical
        ! accuracy to the difference. Note however that when this happens
        ! the it means that the populations of 2 points differ by something smaller
        ! than 10^-8, which is much smaller than all the errors that we can
        ! possibly make (in my opinion, no proof)

        if(last_fpf .spne. fpf) then
           p_val(icat)  = last_fpf - fpf
        elseif(icat == 1) then ! the first term is exactly 1 - so...
           p_val(icat) = max(one_minus_fpf, min_prob)
        else ! since it is 1 - x, signs swap
           p_val(icat) = max(one_minus_fpf - last_one_minus_fpf, min_prob)
        endif

        if(last_tpf .spne. tpf) then
           q_val(icat)  = last_tpf - tpf
        elseif(icat == 1) then
            q_val(icat)  =  max(one_minus_tpf, min_prob)
        else
            q_val(icat)  =  max(one_minus_tpf - last_one_minus_tpf, min_prob)
        endif

        ! save the components that make the TPF and FPF of the current cutoff
        ! which will serve as lower cutoff (which has a larger TPF and FPF)
        ! for the following category
        last_tpf = tpf
        last_one_minus_tpf = one_minus_tpf
        last_fpf = fpf
        last_one_minus_fpf = one_minus_fpf
enddo construct_p_i_and_q_i

! Since the last category is between the n-1 cutoff and vcmax, and at vcmax
! FPF = TPF = 0, we use only the cutoff from the previous category
p_val(num_cat)  = max(last_fpf, min_prob)
q_val(num_cat) =  max(last_tpf, min_prob)

! build the log likelihood and take the minus for the minimization procedure
! (if we use a minimization to maximize the log likelihood). Note that this could be computed more
! quickly by not always using logs, it becomes harder to make it stable, and it doesn't really save much
! cpu time, so we didn't do it.


log_like = - ( dot_product(catn(1:num_cat),log(p_val(1:num_cat))) + dot_product(cats(1:num_cat),log(q_val(1:num_cat))) )

! verify that no overflows or underflows happened
! if it did, reject the step made by the calling program (here set for the TOMS611 optimizer)
if( .not. check_number(log_like)) then
     debugging_info:  if(pass_idebug  == 1)  then
       continue! write(msg,*) "value of d_a & c ", d_a_par, c_par    ; call DisplayState(msg)
       continue! write(msg,*) "bound for cutoffs: ", vc_bound        ; call DisplayState(msg)
       continue! write(msg,*) "WARNING:: OVERFLOW OR UNDERFLOW IN LIKE_PBM ", log_like ; call DisplayState(msg)
       ! write out the values of the cutoffs
       continue! write(msg,*) " CAT CUTOFF_VALUE "        ; call DisplayState(msg)
       do icat=1,num_cat-1
             continue! write(msg,*) icat, param_vec(2+icat) ; call DisplayState(msg)
       enddo
       ! Add the values of TPFs and FPFs
       continue! write(msg,*) " CUTOFF_VALUE FPF's TPF's and their complements to 1 " ; call DisplayState(msg)
       do icat=1,num_cat-1
          call fpf_PBM(d_a_par, c_par, param_vec(2+icat), fpf, one_minus_fpf)
          call tpf_PBM(d_a_par, c_par, param_vec(2+icat), tpf, one_minus_tpf)
          continue! write(msg,*) icat, fpf, tpf  ; call DisplayState(msg)
          continue! write(msg,*) icat, one_minus_fpf, one_minus_tpf ; call DisplayState(msg)
       enddo
       ! Print pieces of the likelihood function
       continue! write(msg,*) "CAT TRUTH       #CASES     P(OR Q)     CATN(I)*LOG(P) "  ; call DisplayState(msg)
       do icat=1,num_cat
            if(catn(icat) /= 0 ) then
                     continue! write(msg,"(2x,i2,' ACT-NEG ',i3, 2(1x,d14.6))")icat, catn(icat), p_val(icat),&
                     continue!                                            catn(icat) * LOG( p_val(icat))
                     call DisplayState(msg)
            endif
            if(cats(icat) /= 0) then
                     continue! write(msg,"( 2x,i2,' ACT-POS ',i3,2(1x,d14.6))")icat,cats(icat), q_val(icat), &
                     continue!                                          cats(icat) * LOG( q_val(icat))
                     call DisplayState(msg)
            endif
       enddo
     endif debugging_info
     ! Set the parameters to reject this step
     log_like   = huge(d_a_par)
     nf = 0
endif


contains
    !----------------------------------------------------------
    !-----------------------------------------------------------
    function check_cutoffs(err_msg) result (test_res)
    ! checks if the cutoffs are within the region where they should be and whether they
    ! are in the correct order. It is an internal function, this is why nothing is passed
    implicit none
    logical test_res
    character(len=line_length), intent(out):: err_msg

    ! Check whether the current move attempted by the optimizer is going out of the space where the
    ! cutoff parameters are mathematically meaningful
    ! nf == 0 is a flag for the optimzier TMO.f telling that the move was into a non
    ! meaningful place. Here we assume that the cutoffs are ranked, if they aren't they will
    ! fail the next test anyway

     if(  ( c_par < -c_almost_zero .and. param_vec(2+1)     < vc_bound ) .OR. &
        ( c_par > +c_almost_zero .and. param_vec(2+num_cat-1) > vc_bound )       ) then
        err_msg = "OPTIMIZER produced cutoffs outside of definition domain"
        test_res = .false.
        return
     else
        test_res = .true.
     endif

     ! Check if the optimizer produced a new set of betas
     ! which are in the wrong order. This would mean that the
     ! category boundaries swapped, which is nonsense from the
     ! ROC model point of view. so we return
     ! a huge value, telling that this configuration is in the
     ! wrong configuration space to the optimizer
     if ( any(param_vec(2+1:num_cat) >= param_vec(2+2:num_cat+1)) ) then
         err_msg = "Rejected step: CUTOFFS in wrong order"
         test_res = .false.
     else
         err_msg = "CUTOFFS ARE INSIDE THE DOMAIN AND IN CORRECT ORDER"
     end if

    end function check_cutoffs

!--------------------------------------------------------------------
end subroutine like_PBM
!--------------------------------------------------------------------

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
 subroutine var_partialauc_PBM(d_a_in,c_in, fraction_1, fraction_2, fpf_flag, &
                                          var_d_a,var_c,cov_d_a_c, var_p_auc, ierror)
!---------------------------------------------------------------------------------
! PURPOSE:  compute the variance of the  partial area either vertically (from FPF_1 to FPF_2 or
!           horizontally (from TPF_1 to TPF_2). Follows the equations from
!           Pan & Metz, academic radiology, 4,380, but uses numerical derivatives
! ALGORITHM: Derive the Variance of the partial area from the variances and
!            covariances of the parameters d_a & c  using the
!            first order approximation, known as delta method (see Papoulis).
!            Is starts from the equation:
!
!            p_auc = partial_AUC(...) using numerical first derivatives to compute
!
!            Then uses is first derivatives to compute:
!            Var(AREA) = (D{AREA}/D{d_a})^2  VAR(d_a) + 2  (D{AREA}/D{d_a})*(D{AREA}/D{c}) * COV(d_a,c) +
!            (D{AREA}/D{c})^2 * VAR(c) + ignored terms
!NOTE:       while the function accepts any value of da and c that is acceptable mathematically (i.e., da >=0
!            and -1 <= c <= 1), it makes little sense to use the delta method when the partial AUC is far from linear
!            over the range of da and c values covered by their standard errors.
! April 30 2009 - LP -> changed ierror = 2 to ierror = bad_input. Added also checkes on input values for variances for variance estimation
!                       routines
 implicit none
 ! parameters d_a, and c
 real(kind=double), intent(IN):: d_a_in,c_in
 real(kind=double), intent(in):: fraction_1, fraction_2 ! These are called fractions
       ! because the can be FPF or TPF depending upon which area are we computing
 integer, intent(in):: fpf_flag ! If it is 1, it means that the area computed will
       ! be between FPF_1 and FPF_2, otherwise it means that it will be between TPF_1 and
       ! TPF_2
 ! variance d_a, variance c, covariance d_a,c
 real(kind=double), intent(IN):: var_d_a,var_c,cov_d_a_c
 real(kind=double), intent(OUT):: var_p_auc
 integer,           intent(out):: ierror   ! = 0, computatio is OK
                                           ! = 1 failed calculation
                                           ! = bad_input, input has unacceptable values
                                           !   or fractions have wrong order.
                                           ! = 3 , fractions are almost identical or identical
                                           !
 real(kind=double) :: p_auc ! value of the partial area
 real(kind=double) :: p_auc_mh, p_auc_ph ! value of the partial area
 real(kind=double) :: h ! local parameter for derivative calculations

 real(kind=double) :: d_p_auc_d_d_a ! value of the derivative of the partial area in d_a
 real(kind=double) :: d_p_auc_d_d_a_o ! value of the derivative of the partial area in d_a
 real(kind=double) :: d_p_auc_d_c   ! value of the derivative of the partial area in c
 real(kind=double) :: d_p_auc_d_c_O   ! value of the derivative of the partial area in c

 real(kind=double) :: d_a_par, c_par ! value of the parameters of the proper binormal models as
                                     ! used internally by the function. This involves some bounding
                                     ! and rounding to avoid numerical instabilities.
 real(kind=double) :: d_a_par_mh, d_a_par_ph ! the values of d_a for the upper and lover bounds,as used
                                     ! in the numerical intergration. These values are checked for consistency
 real(kind=double) :: c_par_mh, c_par_ph ! same as was above for d_a

 integer:: i ! iteration index


! call partial auc to verify whether the input is correct, this is done to keep the two functions
! consistent

ierror = 0
call partialauc_PBM(d_a_par, c_par, fraction_1, fraction_2, fpf_flag, p_auc, ierror)
if(ierror /= 0) then
    var_p_auc = 666.0_double
    return
endif

! varify that input variances are meaningful
if( var_d_a < 0.0_double .or. var_c < 0.0_double .or. cov_d_a_c**2 > var_d_a * var_c ) then
    ierror = bad_input
    var_p_auc = 666.0_double
    return
endif

! Set some relatively strict boundaries for d_a and c, to avoid having messy boundaries
! for the numerical derivatives. These derivaties are prone to instabilities because of
! cancellation errors and the possibility of cusps. Of course the calculation will not be
! too good if we have a cusp, but using the delta-method around a cusp doesn't work in
! the firs place, so don't do it.
d_a_par = min( max(d_a_in, .0001_double), 5.0_double)
c_par   = min( max(c_in, -.9999_double), .9999_double)

! FIND DERIVATIVE IN d_d_A
! This is simply a cut in half the interval until the derivative doesn't change much kind
! of algorithm. 20 was chosen because it works for the situations we tested

h =  .01_double ! Suitably small value to start with

find_d_d_a: do i = 1, 20

    ! the increment is used to set two values of d_a, so that it can be checked whether they are
    ! withint the acceptable bounds. These are a little looser than the one set for the initial
    ! value
    d_a_par_mh = max(d_a_par - h, 0.00001_double)
    d_a_par_ph = min(d_a_par + h, 5.1_double)

    ! compute the two partial AUC values and compute the numerical derivative.
    call partialauc_PBM( d_a_par_mh, c_par, fraction_1, fraction_2, fpf_flag, p_auc_mh, ierror)
    call partialauc_PBM(d_a_par_ph, c_par, fraction_1, fraction_2, fpf_flag, p_auc_ph, ierror)
    d_p_auc_d_d_a =  (p_auc_ph - p_auc_mh)/ (d_a_par_ph - d_a_par_mh)


    ! Check whether the derivative is converging in relative or absolute terms, the two constants
    ! are somewhat arbitrary, but seem to work fine.
    if(i >= 2) then
         if( abs( d_p_auc_d_d_a_o - d_p_auc_d_d_a) < .000001_double  .or. &
             abs( (d_p_auc_d_d_a_o - d_p_auc_d_d_a)/d_p_auc_d_d_a ) < .0001_double &
            ) exit find_d_d_a
    endif


    ! load the current term to be checked against in the next iteration
    d_p_auc_d_d_a_o  = d_p_auc_d_d_a

    ! halv the interval size (boundaries permitting)
    h  = h / 2.0_double


enddo find_d_d_a


! FIND DERIVATIVE IN d_C

! This is simply a cut in half the interval until the derivative doesn't change much kind
! of algorithm all works more or less the same as the partial derivative in d_a, so look
! for those comments if anything is unclear in the following

h = .01_double ! Suitably small value to start with

find_d_c: do i = 1, 20

    c_par_mh = max(c_par - h, -.9999_double)
    c_par_ph = min(c_par + h, .9999_double)

    call partialauc_PBM(d_a_par, c_par_mh, fraction_1, fraction_2, fpf_flag, p_auc_mh, ierror)
    call partialauc_PBM(d_a_par, c_par_ph, fraction_1, fraction_2, fpf_flag, p_auc_ph, ierror)

    d_p_auc_d_c = (p_auc_ph - p_auc_mh)/ (c_par_ph - c_par_mh)

    ! the derivative in partial AUC is essentially zero if  da > 2 * c (and dead zero if c = 0).
    ! However, the derivative can appear to be zero also outside of this value because of numerical
    ! cancellation and the limited accuracy of the partial AUC function. Therefore it is important
    ! to check whether the two partial AUC values are too close and then use the previous step as
    ! the best estimate. This is works very well for the situations we tested.
    if( ( abs(p_auc_mh - p_auc_ph) < 1.0e-7_double) .and. &
        ( c_par > .99_double .or. d_a_par < 2.0_double * c_par)  ) then
           d_p_auc_d_c  = d_p_auc_d_c_o
           exit find_d_c
    endif

    ! Check whether the derivative is converging in relative or absolute terms
    if(i >= 2) then
         if( abs(d_p_auc_d_c_o - d_p_auc_d_c) < .000001_double .or. &
             abs( (d_p_auc_d_c_o - d_p_auc_d_c)/d_p_auc_d_c) < .0001_double &
            ) exit find_d_c
    endif


    d_p_auc_d_c_o  = d_p_auc_d_c
    h  = h / 2.0_double

enddo find_d_c

  ! apply the delta method
  var_p_auc =  d_p_auc_d_d_a**2 * var_d_a   +   d_p_auc_d_c**2 * var_c  + &
                2.0_double * d_p_auc_d_d_a * d_p_auc_d_c * cov_d_a_c


!--------------------------------------------------------------------
end subroutine var_partialauc_PBM
!--------------------------------------------------------------------
!--------------------------------------------------------------------



!--------------------------------------------------------------------
!--------------------------------------------------------------------
subroutine partialauc_PBM(d_a_par,c_in, fraction_1, fraction_2, fpf_flag, &
                                   partial_auc, ierror)
! PURPOSE: compute the partial area either vertically (from FPF_1 to FPF_2 or
!          horizontally (from TPF_1 to TPF_2). Follows the equations from
!          Pan & Metz, academic radiology, 4,380. For the proper binormal model
!ALGORITHM:It uses symmetry to compute the horizontal partial areas (the
!          partial areas defined by TPF boundaries). For this reason they
!          are computed as vertial partial areas too.
! NOTE:    Follow are much as possible the same variable naming as in the
!          paper mentioned in the "PURPOSE" section.
! WARNING: The equation for the horizontal partial AUC in the paper is incorrect
!          but it is not used here because of the symmetry.
!
! April 30 2009 - LP -> changed ierror = 2 to ierror = bad_input.
! June 3rd 2009 - LP -> changed the assignment of c_par to include a cutoff to prevent it from being truly
!                        0 or +- 1. This avoids getting tangled in trying to separate numerical instability from
!                        analytical divergence in the numerical calculation of derivatives, that uses
!                        this function and is used to computer the variance of the partial auc (by another function)
!                        Added also a check on the input parameters to make sure they are acceptable values
!                        Note that the numerical derivative could fail because of the rounding in the partial AUC
!                        done here. Currently the rounding is on the 8th decimal place, which is the expected accuracy
!                        of the bivariate normal distribution calculations (however, one should be mindful with rounding
!                        errors, that can be considerable for this specifici calculation).
!July 2009 - LP       -> Changed the return values of error (ierror) to make it consistent with the other functions
 use statistic_functions, only: phi, bivar_normal_distrib
 implicit none

  real(kind=double), intent(in):: d_a_par ! Curve parameters
  real(kind=double), intent(in):: c_in ! Curve parameters
  real(kind=double), intent(in):: fraction_1, fraction_2 ! These are called fractions
       ! because the can be FPF or TPF depending upon which area are we computing
  integer, intent(in):: fpf_flag ! If it is 1, it means that the area computed will
       ! be between FPF_1 and FPF_2, otherwise it means that it will be between TPF_1 and
       ! TPF_2
  real(kind=double), intent(out):: partial_auc
  integer,           intent(out):: ierror  ! = 0, computatio is OK
                                           ! = 1, failed calculation, not used here.
                                           ! = bad_input, one or more fractions are larger than 1 or
                                           !      smaller than 0 or are in wrong order or unacceptable parameter
                                           !      values
                                           ! = 3 , fractions are almost identical or identical


  ! INTERNAL VARIABLES
  real(kind=double), parameter:: tol = 1.0e-13_double ! tolerance in the computation of tpf and fpf
  real(kind = double):: rho
  real(kind = double):: c_par ! curve parameter for internal use
  real(kind = double):: fpf_1, fpf_2 ! The internal calculations are done using only FPF because the
         ! horizontal partial areas are computed using symmetry, thus again as vertical
  real(kind = double):: y_1, y_2 ! see paper
  real(kind = double):: vc_1, vc_2 ! see paper
  real(kind = double):: x_1, x_1p ! see paper
  real(kind = double):: x_2, x_2p ! see paper
  real(kind = double):: F1, F2, F3, F4, F5, F6, F7, F8 ! the 8 different kinds of bivariate normal
     ! distributions used in the calculation of the partial area

  ierror = 0
  partial_auc = 0.0_double

 ! Verify that input cutoff values for TPF or FPF are meaningful or if the parameters are out of bounds.
 if( fraction_1 > 1.0_double .or. fraction_1 < 0.0_double .or. &
     fraction_2 > 1.0_double .or. fraction_2 < 0.0_double .or. &
     d_a_par    < 0.0_double .or. abs(c_in)  > 1.0_double ) then
         ierror = bad_input
         return
 elseif(fraction_1 >= fraction_2) then
         ierror = bad_input
         return
 elseif(fraction_1 .speq. fraction_2) then
         ierror = 3
         return
 endif

 ! Note that the vertical (from TPF_1 to TPF_2)  partial area is identical to a partial area
 ! on the horizonal computed between  FPF_1 = 1 - TPF_2   FPF_2 = 1 - TPF_1, with c = - c
 if( fpf_flag ==1) then
      c_par = sign( min( max(.00000000001_double, abs(c_in)),.9999999999_double),c_in)
      fpf_1 = fraction_1
      fpf_2 = fraction_2
 else
      c_par = - sign( min( max(.00000000001_double, abs(c_in)),.9999999999_double),c_in)
      fpf_1 = 1.0_double - fraction_2
      fpf_2 = 1.0_double - fraction_1
 endif

! Compute the partial area following the paper by Pan and Metz (see header of subroutine)
  rho =  - 1.0_double / sqrt(2.0_double) *  ( ( 1.0_double + c_par ) / sqrt( 1.0_double + c_par**2 ) )
  y_1 =  d_a_par / sqrt(2.0_double)
  y_2 = ( d_a_par / sqrt(2.0_double) ) * ( ( 1.0_double + c_par**2) / (2.0_double * c_par) )

  ! Find the cutoff corresponding to those values of FPF_1 and FPF_2, then compute
  ! the arguments of the phi functions of the proper model
  if( abs(c_par) > 1.0e-4) then
     vc_1 = .9_double * d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * abs(c_par))
  else
     vc_1 = 20.1_double ! value is chosen to match the bound of fpf_find_vc
  endif
  ! Find the vc values corresponding to the cutoffs in the FPF/TPF
  call fpf_find_vc(d_a_par, c_par, fpf_1 , tol,  vc_1)
  vc_2 = vc_1 ! use privious value as initial value
  call fpf_find_vc(d_a_par, c_par, fpf_2 , tol,  vc_2)


  ! Compute the variables for the bi and univariate normal distributions
  ! used to compute the partial area( see aformentioned paper for details)
  x_1  =  z1_fpf(d_a_par, c_par, vc_1)
  x_1p =  z2_fpf(d_a_par, c_par, vc_1)
  x_2  =  z1_fpf(d_a_par, c_par, vc_2)
  x_2p =  z2_fpf(d_a_par, c_par, vc_2)


! For now we ignore the return error messages, they are never used anyway.
  call  bivar_normal_distrib( x_2,   y_1 , rho, F1, ierror)
  call  bivar_normal_distrib( x_1,   y_1 , rho, F2, ierror)
  call  bivar_normal_distrib( x_2,   y_2 , rho, F3, ierror)
  call  bivar_normal_distrib( x_1,   y_2 , rho, F4, ierror)
  call  bivar_normal_distrib( x_2p, -y_2 , rho, F5, ierror)
  call  bivar_normal_distrib( x_1p, -y_2 , rho, F6, ierror)
  call  bivar_normal_distrib( x_2p, -y_1 , rho, F7, ierror)
  call  bivar_normal_distrib( x_1p, -y_1 , rho, F8, ierror)

  partial_auc = (F1 + F3 + F5 + F7) - (F2 + F4 + F6 + F8)

  ! Add term depending upon the sign of c (created by the step function)
  if( c_par > 0.0_double ) partial_auc = partial_auc  - fpf_2 + fpf_1

 ! round to the 8th decimal place. This choice seems to be strike a reasonable
 ! balance between the use of this function in itself and to compute numerical
 ! derivatives
   partial_auc = real(nint(partial_auc * 1.0e8_double)) / 1.0e8_double

!--------------------------------------------------------------------
end subroutine partialauc_PBM
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!--------------------------------------------------------------------
pure subroutine auc_PBM(d_a_par,c_par, auc, ierror)
! PURPOSE: Computes the AUC (area under the curve) of the proper binormal model
!          as difference between phi (the normal cumulative distribution function)
!          and F (the bivariate, normal distribution function), see papers and the
!          function body for the exact form
! NOTE:    return error messages when the bivariate calculation fails or is unstable
!          Note that the normal distibution is never unstable for any real value (as far
!          as we have tested them), so its error is not used
!HISTORY: changed by LP May 2009 to add checks on parameter values

 use statistic_functions, only: phi, bivar_normal_distrib
 implicit none

  real(kind=double), intent(in):: d_a_par
  real(kind=double), intent(in):: c_par


  real(kind=double), intent(out):: auc

  ! ierror:: 0 -> OK, 1 == failed, 2 == c too small, used only phi part, bad_input wrong input
  integer,           intent(out):: ierror


  real(kind=double):: F
  real(kind=double):: rho


! Check for impossible values of the input parameters
if(d_a_par < .0_double .or. abs(c_par) > 1.0_double) then
    ierror = bad_input
    return
endif

  rho = - (1.0_double-c_par*c_par) / (1.0_double+c_par*c_par)

  call  bivar_normal_distrib( - d_a_par / SQRT(2.0_double), 0.0_double , rho, F, ierror)

  if(ierror == 0 ) then
      auc =  phi(d_a_par/SQRT(2.0_double)  )  +  2.0_double * F
  elseif( abs(rho) >= .999999) then
      auc =  phi(d_a_par/SQRT(2.0_double)  )
      ierror = 2
  else
      auc =  phi(d_a_par/SQRT(2.0_double)  )
      ierror = 1
  endif

!------------------------------------------------------------------------
end subroutine auc_PBM
!------------------------------------------------------------------------
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!------------------------------------------------------------------------

!------------------------------------------------------------------------
subroutine compute_MLE(mn, ms, num_cat, catn, cats, idebug, d_a_par, c_par, vc_cutoffs, &
                        log_like, hessian, ierror, err_msg)
!------------------------------------------------------------------------
! PURPOSE: compute the Maximum Likelihood estimate of the curve parameters using the initial estimates of
!          the cutoffs and parameter curves as sent by the calling program. It returns also the hessian matrix
!          to be used to compute variances.
!ALGORITHM: the optimization is done using TOMS611, a Newton-Raphson trust zone optimizer. Sometimes the initial
!          estimates are slightly changed and the calculation is repeated. this has to to with especially unstable
!          regions around the guessing line (AUC = .5, minimum value possible for a proper model). It happens very
!          rarely and when it does, it cycles very few times. The algorithm does check whether the optimization
!          produced meaningful numbers.
!          If the optimization wasn't successful, the algorithm returns Log_like = - Huge(1.0_double), i.e. the
!          smallest value possible.
! NOTE:   This procedure does not check if the parameters are in a specific numerical
!         space, but only if they are mathematically acceptable
! HISTORY: 1-24-06, LP, changed max_iter, made it dependent upon the number of parameters. When the number of
!          parameters is large, it does not help to iterate. In these situations if one initial estimate fails
!          it is better to drop it and try another one. We have never seen multiple maxima in this situation so
!          it is useless to obtain multiple solutions. It is ad hoc and might need improvement.
 use debugging, only: DisplayState
 use ll_derivatives
 use toms611, only: deflt, humsl
 use l_algebra, only: m1tom2
 implicit none

 integer, intent(in):: mn ! number of actually negative cases
 integer, intent(in):: ms ! number of actually positive cases
 integer, intent(in) :: num_cat ! Number of categories as created by catgrz
 integer, dimension(num_cat), intent(in):: catn, cats ! arrays containing categorical data
 integer, intent(in):: idebug ! whether to print log information or not. This applies only to log, Warnings and
                              ! errors are printed anyway

 real(kind = double),intent(inout)            :: d_a_par, c_par ! input: initial estimates, output: final
 real(kind = double), dimension(num_cat-1), intent(inout) :: vc_cutoffs  ! in: initial est. out: MLE cutoffs
 real(kind = double), intent(out)           :: log_like ! value of the log likelihood function at the final point
 real(kind = double), dimension(num_cat+1,num_cat+1), intent(out)  :: hessian
 integer, intent(out)                                   :: ierror ! Error flag about the MLE fit
 character(len = line_length), intent(out):: err_msg

 real(kind=double), dimension(num_cat+1,num_cat+1):: gradient
 real(kind=double), dimension( (num_cat+1)*(num_cat+2)/2) :: t_hessian ! hessian,triangular
 ! parameters defined for the transactions on mathematical software algorithm 611 optimizer
 integer :: lv  ! size of the work array of the optimizer
 real(kind=double), allocatable, dimension(:):: v
 integer, parameter :: liv = 60
 integer, dimension(1):: uip
 integer, dimension(liv):: iv
 integer:: num_param ! number of parameters in the parameter vector
 real(kind=double),dimension(num_cat+1):: urp ! It can be used to transfer information between
             ! the functions called by the TOMS611, currently non used
 real(kind=double),dimension(num_cat+1):: param_vec
 real(kind=double), dimension(num_cat+1):: d_scale
 real(kind=double):: d_a_par_in, c_par_in
 integer:: iter ! how many times the optimizer started on a specific d_a, c optmization
 integer:: max_iter ! Maximum number of iterations that will be attempted

 integer::i, nf ! nf is a flag for the optimization procedure

 real(kind=double) :: vc_bound ! finite boundary for vc (see Metz & Pan)
                                ! sometimes an upper sometimes a lower

 real(kind=double) :: fpf, tpf
 character(len = line_length):: msg ! text line to be sent use to log information about the run


 intrinsic random_number


 ! Set the return error to the default (OK)
 ierror = 0

 num_param = num_cat + 1 ! the optimization parameters are d_a, c, and num_cat - 1 categories

 ! This is a little ad hoc. Basically the idea is that if there are many categories
 ! it becomes really difficult to converge when the the initial values are not
 ! very good. Since in these situation there aren't any multiple maxima, it does not
 ! make much sense to try too hard to converge either, so the values are dropped.
 if (num_param > 100) then
    max_iter = 2
 else
    max_iter = 20
 endif

 ! Allocate the matrices that will contain the derivatives of the
 ! p_i's and q_i's in the parameters and the values of the p_i's and q_i's
 allocate(dp_dpar(num_cat, num_param),dq_dpar(num_cat, num_param))
 allocate(p_val(num_cat), q_val(num_cat))

 ! size of optimizer's work arrays - refer to the optimizer for information about it, and set it to zero
 lv = 78 + (num_param)*(num_param + 12) ! size of the
 allocate(v(lv))
 v = 0

 ! Store the initial value of d_a and c for use in case the first attempts of the MLE fail
 ! and the procedure need to be iterated with different initial estimates (sometimes needed
 ! when curves are very close to the chance line)
 d_a_par_in = d_a_par
 c_par_in = c_par


 ! Copy the initial estimates into the parameter vector for the optimizer
 param_vec(1) = d_a_par
 param_vec(2) = c_par
 param_vec(3:num_param) = vc_cutoffs(1: num_cat - 1) ! cutoffs are in between categories

! This loop is rarely used. Only for curve that are barely or not at all above the guessing line (AUC = .5)
! Since the maximum tends to be a cups in such a situation, surrounded by a couple of lines where the hessian
! is singular, it happens that the algorithm can use being restarted from a slightly different point. It seems
! like the problem has more to do with cutoffs than proper parameters. It never uses more than a couple of
! iterations even the rare times it iterates.
optimize: do  iter = 1 ,  max_iter
   ! Call the routine that initializes the TOMS611 optimizer.
   call deflt(2,iv,liv,lv,v)
   iv(21) = 0
   iv(17) = 5000 ! Maximum number of function evaluations to attempt, only very small curves use more than a few 10's
   iv(18) = 5000 ! Maximum number of overall optimization steps " "
   iv(1)= 12 ! tell sumsl that deflt (initialization) was already called
   iv(16) = 1 ! scaling factor scheme to use for the quasi-Newton, for the definition of the trust zone
   d_scale(1:num_param)  = 1.0_double ! scaling factors to define the trust zone, start simple
   ! Call the Trust zone Newton-Raphson (Fisher scoring) to find the Maximum of the Log Likelihood function
   call humsl(num_param, d_scale ,param_vec, like_PBM, ddlike_PBM,iv,liv, lv,v, uip,urp, like_PBM)
   ! Copy the final values from the optimization vector into their variables
   d_a_par = param_vec(1)
   c_par   = param_vec(2)
   vc_cutoffs = param_vec(3:num_param)
   ! in the proper model the cutoffs become bounded in one direction, either for large positive or large
   ! negative values. The reason has to do with the fact that log likelihood ratios normally are bounded
   vc_bound = d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * sign( max(.00000000001_double, abs(c_par)),c_par))
   ! Check whether the optimization was successful, iv(1) is the error flag of the
   ! optimization, details of the optimization run are in the "optimization.info" file
   if( iter < max_iter .AND. & ! If this is the last iteration, we aren't going to try anything more
             (iv(1) == 8 .OR. iv(1) == 9  .OR. iv(1) == 10 .OR. iv(1)== 63 .OR. & ! unstable opt
             ( c_par < -c_almost_zero .AND. vc_cutoffs(1) < vc_bound ) .OR. & ! cutoffs out of
             ( c_par > +c_almost_zero .AND. vc_cutoffs(num_cat-1) > vc_bound )  .OR. &! domain
             .not.check_cutoff_order() ) & ! cutoffs have wrong order
   ) then
       ! send a warning in the log file about a failed attempt in optimization
       call print_convergence(iv(1), ierror,err_msg) ! explicit what kind of convergence was obtained by the optimizer
       if (idebug == 1) then
          call DisplayState(err_msg)
          ! LP April 2007 fixed to prevent crashing while printing
          continue! write(err_msg,"(I2,' TOMS611 call:: da & c =', 2(1x,e16.8))") iter, d_a_par, c_par
          call DisplayState(err_msg)
       endif
       ! Restart from the original initial values, but change d_a a little. If d_a is too small, it can't be a relative
       ! change, so we add some absolute value - probably it would also work with the absolute only
       if(d_a_par_in  > .2_double) then ! changing d_a a little took care of the problem
                param_vec(1)  = d_a_par_in*(1.0_double + real(iter)/max_iter)
                       ! (Note that d_a is only increased, since all problems seem to happen with small d_a
                       ! for all the 5,000,000 datasets we tried. - 5-05-05)
        else
               param_vec(1)  = d_a_par_in + (.5_double + real(iter)/max_iter)
        endif

       param_vec(2)  = c_par_in
       ! Recompute the cutoffs using the new values of d_a
       call initialize_cutoffs(param_vec(1),param_vec(2), mn, ms, &
                               num_cat, catn, cats, vc_cutoffs, idebug, ierror, err_msg)
       ! Check whether the new cutoffs initialization succeeded
       if(ierror /= 0) then
          ! LP put a huge value to reject this step (here we are optimzing - log likelihood
          log_like= huge(d_a_par)
          err_msg = " initialize_cutoffs failed " // err_msg(:line_length - 27)
          if(idebug == 1) call DisplayState(err_msg)
          cycle ! try again unless iteration limit is exceeded
       else
          param_vec(3:num_param) = vc_cutoffs(1:num_cat-1)
       endif
   else
        exit optimize ! either it worked or we exceeded max_iter or the error can't be fixed
   endif
enddo optimize

call print_convergence(iv(1),ierror,err_msg) ! explicit what kind of convergence was obtained
                                             ! by the optimizer

! Update the finite bound value for the cutoffs
vc_bound = d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * sign( max(.00000000001_double, abs(c_par)),c_par))

! Note that the procedure returns a - "infinite" log likelihood if something went wrong
if(ierror /= 0) then
          ! MLE failed return the error code and the error message from the optimizer
          ! We might want to do something else in the future, so I left the line
          ierror = 1
          log_like = -huge(1.0_double)
! Check whether the parameter are in bound. It used to happen in the past, it never happened with the
! most recent versions, but better check
elseif(.not.test_parameters(d_a_par, c_par) ) then
        ierror = 1
        err_msg = "Error:: compute_MLE: parameters out of bounds "
        log_like = -huge(1.0_double)
! Check whether the optimizer returned values ouside the mathematical definition
! domain for the cutoffs (Cutoffs are supposed to be ranked, which is checked
! right after this instruction
elseif(  ( c_par < -c_almost_zero .and. vc_cutoffs(1)   < vc_bound ) .or. &
     ( c_par > +c_almost_zero .and. vc_cutoffs(num_cat-1) > vc_bound )       ) then
        err_msg = "ERROR:: compute_MLE -> OPTIMZER RETURNED PARAMETERS OUTSIDE DOMAIN"
        ierror = 1 ! return an error to calling program
        log_like = -huge(1.0_double)
! Check if the optimizer finished the optimization with a set of betas which are in the wrong order
! This would mean that the category boundaries swapped, which is nonsense from the
! ROC model point of view.
! If 2 points have identical TPF/FPF, clearly the flipping  in the values is nothing to care for, so we
! can switch them  to the right order and keep going. This never happened in the newer versions of proproc
elseif (.not.check_cutoff_order()) then
        if(idebug == 1) then
           call DisplayState("ERROR:: compute_MLE -> OPTIMZER RETURNED CUTOFFS WITH WRONG ORDER")
           continue! write(msg,*) " i     cutoff    fpf          tpf  "
           call DisplayState(msg)
           do i=1,num_cat-1
                  call fpf_PBM( d_a_par, c_par, vc_cutoffs(i), fpf)
                  call tpf_PBM( d_a_par, c_par, vc_cutoffs(i), tpf)
                  continue! write(err_msg,"(i2,3(1x,e22.15))") i, vc_cutoffs(i), fpf,tpf
                  call DisplayState(msg)
           enddo
        endif
        err_msg = "ERROR:: compute_MLE -> OPTIMZER RETURNED CUTOFFS WITH WRONG ORDER"
        ierror = 1 ! return an error to calling program
        log_like = -huge(1.0_double)
else ! The fit worked out
        nf = 1 ! This is a parameter which is used in optimization, here it does not serve any purpose
                  ! but needs to be set because it is intent(inout) in like_PBM. Call can be improved.

        ! compute the likelihood value at the maximum, to return to the calling program
        call like_PBM(num_param, param_vec, nf , log_like, uip)

        log_like = - log_like ! the optimizer is a mimimizer so it works with - log, so we need to correct that

        ! Compute the hessian to estimate the variances
        call ddlike_PBM(num_param, param_vec, nf, gradient, t_hessian, uip)

        ! Load full matrix from the upper triangular form - used by the optimizer
        call M1TOM2(num_param, t_hessian, hessian)
        hessian = - hessian ! Since we optimize the -likelihood, we need to change sign to return the
                   ! actual likelihood, note that the inversion will have to take the minus again, but
                   ! I preferred to return straight object as opposed to save millisenconds fiddling with
                   ! signs
endif

 ! Deallocate all the arrays in the module used to transfer data, to prevent memory
 ! problems (happen with some compilers, for example Lahey...)
 deallocate(dp_dpar, dq_dpar)
 deallocate(p_val,q_val)
 deallocate(v)

! --------
 contains

!------------------------------------------------------------
  subroutine print_convergence(flag,errflg,err_msg)
!------------------------------------------------------------
! Interpretes the flag that describes the quality of the optimization phase just concluded by the program
! sumsl (ot HUMSL - if hessian is used for scoring) and convert it into more understandable sentences
  implicit none
  integer, intent(IN):: flag
  integer, intent(out):: errflg ! if 0, fit is acceptable, if 1 it isn't
  character(len=line_length), intent(out):: err_msg

  select case(flag)
   case(3)
     err_msg = " Achieved x-convergence for parameter vector "
     errflg = 0
   case(4)
     err_msg = " Achieved rel. function convergence for the log likelihood "
     errflg = 0
   case(5)
     err_msg = " Achieved convergence both for parameter vector and log likelihood "
     errflg = 0
   case(6)
     err_msg = " Achieved absolute function convergence for the log likelihood "
     errflg = 0
   case(7)
     err_msg = " Hessian is locally singular and so steps are unlikely to produce a better LL"
     errflg = 0
   case(8)
     err_msg = " ERROR:: convergence at noncritical point, # "
     errflg = 1
   case(9:10)
     err_msg = " ERROR: convergence failed, iteration limit "
     errflg = 1
   case(18)
     err_msg = "ERROR:: scaling vector, d_scale,  has negative value!!! "
     errflg = 1
   case(63)
     err_msg = "ERROR:: cannot compute Log Likelihood at initial estimates"
     errflg = 1
   case(65)
     err_msg = "ERROR:: cannot compute GRAD{LL} at initial estimates "
     errflg = 1
   case default
     err_msg = "ERROR:: error in optimization process, please send input file to roc@bsd.uchicago.edu "
     errflg = 1
   end select

!------------------------------------------------------------
  end subroutine print_convergence
!------------------------------------------------------------

!------------------------------------------------------------
 logical function  check_cutoff_order() result (test_value)
 ! PURPOSE: internal function check whether the order of the cutoffs is ascending.
 !          if it isn't, it means that something went wrong in the optimization (likely an
 !          overstep). If the wrong order produces small differences, the cutoffs are
 !          reordered, otherwise and error is returned
 !NOTE:     Although it is a function it does change some values, note however that these are
 !          not calling arguments, but passed by the fact that the function is internal. Moreover
 !          it does not make any change which has a measurable effect on the Likelihood function

  implicit none

  integer:: i
  real(kind=double) :: swap
  real(kind=double) :: fpf1, tpf1
  real(kind=double) :: fpf2, tpf2

  test_value = .true.
  do i=2,num_cat-1
   if ( (vc_cutoffs(i-1) >= vc_cutoffs(i) ) ) then ! Check whether the cutoffs are in the right order
      ! if they are, make sure that it matters, i.e., that the FPF's and TPF's are in fact different
      call fpf_PBM( d_a_par, c_par, vc_cutoffs(i-1), fpf1)
      call tpf_PBM( d_a_par, c_par, vc_cutoffs(i-1), tpf1)
      call fpf_PBM( d_a_par, c_par, vc_cutoffs(i  ), fpf2)
      call tpf_PBM( d_a_par, c_par, vc_cutoffs(i  ), tpf2)

      if(  abs(fpf1 - fpf2) < p_q_closeto_zero .and. abs(tpf1 - tpf2) < p_q_closeto_zero  ) then
          swap      =  vc_cutoffs(i)
          vc_cutoffs(i)   =  vc_cutoffs(i-1)
          vc_cutoffs(i-1) = swap
      else
        test_value = .false.
        exit
     endif
   endif
 enddo

 end function check_cutoff_order

!------------------------------------------------------------
 end subroutine compute_MLE
!------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine initialize_cutoffs(d_a_par_in, c_par_in, num_normal_cases, num_abnormal_cases, num_cat, &
                               catn, cats, beta, idebug,  ierror, err_msg)
!---------------------------------------------------------------------------------------------------
! PURPOSE: compute initial values for the cutoffs in the
!         latent variable space (v_c) to be returned as
!         the array beta

! ALGORITHM: The subroutines determines the initial estimates by
!            assuming that each cutoff value is partitioning the
!            cases into 2 categories, the cases with classification
!            smaller than the cutoff  and the ones larger
!            First an upper bound for the largest cutoff is found,
!            then the algorithm proceeds using the previous cutoff
!            as upper bound for the following one. This insures also
!            that the cutoffs have the correct ranking.
!            The initial estimate is found as the value that maximizes
!            the log likelihood for the 2 category data. The optimization
!            is done using the Brent method directly on the log likelihood
!            function
! NOTE0:     IF cutoffs are almost identical, they are moved apart to prevent ovelapping
!            and consequent divergence of functions and derivative values. It should be reminded that
!            these are initial values, so they don't have to be precise and can't be precise anyway
! NOTE1:     the Log likelihood was used instead of the derivative to avoid the
!            overflow/underflow issues in computing the ratio of sums and differences
!            of normal  distribution functions and densities, for which no elegant
!            solution was found (I do not like if statements).
! NOTE2:     The Accuracy of the initial estimates affects the results a lot since
!            we change in a chaotic way the starting point of the cutoff
!            optimization for similar parameters , and since the Likelihood
!            function can be rather flat this effect is magnified.
! NOTE3:     the algorithm uses the symmetry between +c and -c (that is
!            if you change the positives with the negatives and change the sign
!            of the cases, one obtains a curve with same area but opposite value of
!            c. Or alternatively, changing z -> -z and phi with 1- phi for all the
!            terms involved in the TPF or FPF for each operating point, the equations for
!            + c and -c are the same.
! WARNING:  the last cutoff is pointless
! HISTORY: september 2006 - changed the search of the next cutoff by adding a small value (the approximation used
!           (in detetermining the value itself) to prevent cutoffs from ovelapping.

 use debugging, only: DisplayState
 use gen_numerics, only: brent
 use proproc_initial_cutoffs ! needed to exchange data around the canned linear minimization brent
 implicit none
 real(kind=double), intent(in):: d_a_par_in, c_par_in
 integer, intent(in):: num_normal_cases
 integer, intent(in):: num_abnormal_cases
 integer, intent(in):: num_cat ! Number of categories as found by catgrz
 integer, intent(in), dimension(1:num_cat):: catn, cats ! arrays containing categorical data

 real(kind=double), dimension(num_cat-1), intent(out):: beta ! initial estimates
                              ! of cutoffs
 integer, intent(in):: idebug ! whether to print log information of warnings
 integer, intent(out):: ierror ! Error flag, for failed initialization
                               ! of cutoffs 0 OK, 1 failed
 character(len=line_length), intent(out) :: err_msg


 real(kind=double),parameter:: max_err = 1.0e-3_double ! max uncertainty in est.
                            ! initial value for cutoff, the older algorithm needed a very
                            ! precise estimate (10-10), this isn't true anymore (although
                            ! is has no speed impact on the code, since it is done only once
 real(kind=double) :: previous_bound, middle, next_bound ! cutoff values that bound
                                             ! a maximum of the likelihood function
 real(kind=double):: log_like ! value of the log likelihood at the initial estimate of the
                             ! cutoff, currently not used
 real(kind=double):: inf_bound, model_bound ! value for the maximum or minimum value of the
                                                 ! cutoffs
 integer:: j_cut, j ! loop counter over the cutoffs
 integer:: frst_cat, lst_cat, incr


 ierror = 0 ! set to zero the error flag
 ! initialize the values in the module initil_cutoffs, that will be used around brent (minimization)
 d_a_par = d_a_par_in
 c_par = c_par_in
 beta = 0.0_double ! Initialize the cutoffs array

 ! the algorithm proceds from the finite upper bound (c is never really zero, it is not allowed by the code)
 ! in the direction of the other cutoffs until it has found all of them. This means that we start from the
 ! last category moving down. the opposite is true for c negative

 if (c_par >= 0 ) then
     frst_cat = num_cat
     lst_cat  = 2
     incr     = -1
     j_cut    = num_cat - 1
     cat_neg = 0
     cat_pos = 0
 else
     frst_cat = 1
     lst_cat  = num_cat - 1
     incr     = 1
     j_cut    = 1
     cat_neg = num_normal_cases
     cat_pos = num_abnormal_cases
 endif

 !This is the analytical extremal values for the cutoffs, often the proper functions are unstable in this
 ! neighborhood and this is why we look for a numerical one slightly more inside the domain than this value,
 ! and use it unless we can't find it
 model_bound = d_a_par * sqrt( 1.0_double + c_par**2) / (4.0_double * sign( max(.00000000001_double, abs(c_par)),c_par))
 ! Find, if possible, the value after which the 2 cat LL for the first cutoff starts to decrease
 call set_infinite_upper_bound(d_a_par,c_par, cat_neg - incr*catn(frst_cat) ,cat_pos -incr * cats(frst_cat)  , &
                      num_normal_cases, num_abnormal_cases, inf_bound,ierror,err_msg)

 if(ierror /= 0) then
         if(idebug == 1) call DisplayState(err_msg)
         ! If the routine failed, default to the analythical maximum value, and use half of it as backup (never used)
         inf_bound = model_bound
         previous_bound =  model_bound*.5
 else
 ! Note that the analytical upper (or lower) bound is stored, so in case the other value is not working it can be used
 ! as well.
         previous_bound =  inf_bound
         inf_bound = model_bound
 endif

! the loop is build in this way so that we can go from the last category to the first or from the first to
! the last, depending upon the sign of c. This makes the algorithm symmetrical, removing bias in the sign of c
 CATEGORY_BOUNDARIES: do j =  frst_cat, lst_cat, incr
    cat_neg = cat_neg - catn(j)*incr   ! compute the total number of cases before (or after depending upon sign of c)
    cat_pos = cat_pos - cats(j)*incr   ! the cutoff whose initial estimates we are trying to find
    ! bracket the minimum for the next cutoff
    call set_next_bound(d_a_par,c_par, cat_neg, cat_pos, num_normal_cases, num_abnormal_cases, &
                      previous_bound,  middle, next_bound, ierror, err_msg)

! Check if there was a bracketing problem with the maximum of the 2 category likelihood function for the first cutoff
    if(ierror /= 0 .AND. j == frst_cat) then ! In this case, we can try to see if we can use a the "backup" estimate of
                                             ! the first cutoff. Normally it works (in fact this in itself is used
                                             ! very rarely and typically for very small areas and noisy data)
          call set_next_bound(d_a_par,c_par, cat_neg, cat_pos, &
                         num_normal_cases, num_abnormal_cases, &
                         inf_bound,  middle, next_bound, ierror, err_msg)
    endif

! If the maximum could not be bracketed, we return a failure to the calling program
! with the same error code and error message
    if(ierror /= 0)  return

    ! Determine the optimal point using the function brent, which in turn uses function
    ! minus_two_cat_ll because it is a minimization algorithm
    log_like = brent(previous_bound,middle,next_bound, minus_two_cat_ll,max_err,beta(j_cut),ierror)
    if(ierror /= 0 .or.  .not.check_number(log_like) )  return ! Brent failed, need to return - introduced july 1st 2005
    previous_bound = beta(j_cut) + max_err*incr ! Use this optimal value as upper bound for the next cutoff, but add some to
    ! prevent possible ovelapping of cutoffs, that makes the starting of the algorithm very unstable - sept 9 2006
    j_cut = j_cut + incr
 end do CATEGORY_BOUNDARIES

!----------------------------------------------------
!----------------------------------------------------
end subroutine initialize_cutoffs
!----------------------------------------------------
!----------------------------------------------------

!-----------------------------------------------------
!-----------------------------------------------------
function minus_two_cat_ll (cut_off,ierror) result( ll)

 use problem_data
 use proproc_initial_cutoffs
 implicit none

! PURPOSE: generate the negative of the two category log likelihood to be used by optimization
!          procedures
! WARNING: need to change calling scheme, I do not like to change or return arguments for functions
 real(kind=double), intent(in):: cut_off
 integer, intent(out):: ierror
 real(kind = double):: ll

 call two_cat_ll(d_a_par, c_par,cat_neg,cat_pos, num_normal_cases,num_abnormal_cases, cut_off, ll, ierror)
 ll = -ll

end function minus_two_cat_ll

!-------------------------------------------------------
!-------------------------------------------------------
 subroutine set_next_bound(d_a_par, c_par, cut_neg, cut_pos, num_neg, num_pos,upper_bound, &
                            middle,lower_bound, ierror, err_msg)
! PURPOSE: Determine a lower bound for a 2 category log likelihood function
!          for an roc curve determined by d_a_ and c_par.
! ALGORITHM: Use the function value to determine whether the function is decreasing or increasing
!            uses the fact the the c determines how quickly will the function change

! NOTE: Need to check that the algorithm does not produce oscillations around values when it searches
!       for a lower bound. In particular this can happen when the function is flat around the upper bound,
!       because of the very low probability of some categorical data (normally caused by a very-very bad fit).
!       this might be related to a small asymmetry noticed in the d_a,c space for the log likelihood. RIght now,
!       If the algorithm overshoots a maximum and then comes back, if it finds a value which is the same as the
!       log likelihood at the upper bound it tries to go outside a little bit and if it is still flat,
!       it takes that value as middle point (point with highest likelihood). This should not create problems
!       To the optimizer.
 implicit none

 real(kind=double), intent(IN) :: d_a_par, c_par ! ROC curve parameters
 integer, intent(IN) :: cut_neg, cut_pos ! number of act. neg. and act. pos above this cutoff
 integer, intent(IN) :: num_neg, num_pos ! total number of act. neg. and act. pos
 real(kind=double), intent(IN) :: upper_bound

 real(kind=double), intent(out) :: lower_bound, middle
 integer, intent(out):: ierror ! Error flag
 character(len=line_length), intent(out):: err_msg

 real(kind=double):: scale ! scale factor to use while trying to bound the maximum
 real(kind=double):: old_lbound ! value used in the iterative scheme
 real(kind=double):: ll, old_ll, upper_ll, middle_ll  ! 2 category log likelihood
 integer:: iter ! iteration counter
 integer, parameter:: max_iter = 500 ! Max number of iterations used to determine the lower
                                 ! bound of this cutoff
 logical:: away ! we are moving away from upper bound or not

 ierror = 0 ! initialize error flag


! Take the slowest scale (either the TPF or FPF, depending upon the sign of c
 scale =  1.0 / (1.0 + min(abs(c_par),  c_almost_one) ) * sign( 1.0_double,  c_par)

 ! determine the value of the LL for the upper bound
 call two_cat_ll(d_a_par, c_par,cut_neg,cut_pos, num_neg,num_pos, upper_bound, upper_ll, ierror)
 ! if c is negative increment should increase the function, the opposite is it is positive
 old_lbound = upper_bound - scale
 away = .true.
 do iter = 1, max_iter ! make sure that the first point has a larger LL than the upper bound
      call two_cat_ll(d_a_par, c_par,cut_neg,cut_pos, num_neg,num_pos, old_lbound, old_ll, ierror)
! for desperate debugging, when initial estimates don't seem to work
!           write(*,"('m ', i3, (1x,d14.8),1x, d22.15)") iter, old_lbound, old_ll
      if (old_ll > upper_ll) then
          exit ! Theoretically the search could have overshoot the part where the function varies and ended up
               ! in a part of the cutoff space where the function doesn't change anymore where at the same time
               ! the 2 cat LL is larger. This never happened in practice so we ignored it
               ! if you determine problems here, you can set a step back and forward a little to
               ! check whether the function changes right around this point and then take the
               ! relevant countermeasures when it does not happen (namely go back) - 5-06-05
      elseif(  (old_ll .speq. upper_ll) .and. away) then
                  old_lbound = old_lbound - scale*iter ! geometric progression away
                  away = .true.
      else
        old_lbound = (old_lbound + upper_bound)/2
        if (upper_bound .speq. old_lbound) exit
        away = .false.
      endif
 enddo

 if(iter > max_iter) then
      continue! write(err_msg,"('set_next_bound:: middle not found, param:',2(1x,f6.3))") d_a_par, c_par
      ierror = 1
      return
 else
       middle = old_lbound
       middle_ll = old_ll
 endif

 call two_cat_ll(d_a_par, c_par,cut_neg,cut_pos, num_neg,num_pos, old_lbound, old_ll, ierror)

 do  iter=1, max_iter
       if(iter < 10) then
          lower_bound = old_lbound - scale*iter
       else
          lower_bound = old_lbound - scale*iter**2
       endif
       call two_cat_ll(d_a_par, c_par,cut_neg,cut_pos, num_neg,num_pos, lower_bound, ll, ierror)
       if(ierror /= 0) then
             continue! write(err_msg,"('set_next_bound:: D(LL) lost meaning, param:',2(1x,f6.3))") d_a_par, c_par
             ierror = 1
             return
       endif
       if( ll < old_ll ) then
           if(middle_ll < old_ll) then
              middle = old_lbound
              middle_ll = old_ll
           endif
!           write(*,"('lls    ', 3(1x,d16.8))")  upper_ll, middle_ll, ll
           exit ! if the derivative is negative, we passed the maximum, which means
           ! that we have found the other bound
       else
!          write(*,"('l ',i3,1x,3(f16.8))") iter, lower_bound, ll, old_ll
          old_ll = ll
          old_lbound = lower_bound
       endif
 enddo

! verify that convergence was obtained
 if(iter > max_iter) then
       continue! write(err_msg,"('set_next_bound:: not found, param:',2(1x,f6.3))") d_a_par, c_par
       ierror = 1
       return
 endif


 end subroutine set_next_bound

!--------------------------------------------------
!--------------------------------------------------

!--------------------------------------------------
!--------------------------------------------------
 subroutine set_infinite_upper_bound(d_a_par, c_par, cut_neg, cut_pos, num_neg, num_pos, &
                                     upper_bound, ierror, err_msg)
 ! PURPOSE: Determine an upper bound for a 2 category log likelihood function
 !          for the  operating point whose cutoff is closest to the finite boundary of the cutoff space,
 !          as determined by fpf and tpf, for a population roc curve  determined by d_a_ and c_par

 implicit none

 real(kind=double), intent(IN) :: d_a_par, c_par
 integer, intent(IN) :: cut_neg, cut_pos ! number of act. neg. and act. pos above this cutoff
 integer, intent(IN) :: num_neg, num_pos ! total number of act. neg. and act. pos
 real(kind=double), intent(out) :: upper_bound
 integer, intent(out):: ierror
 character(len=line_length), intent(out):: err_msg

 real(kind=double):: scale
 real(kind=double):: old_ubound
 real(kind=double):: ll, old_ll  ! derivative of the 2 category log likelihood
 integer:: iter !
 integer, parameter:: max_iter = 1000 ! Max number of iterations used to determine the upper
                                      ! bound of this cutoff. We use a lot of iterations (which are
                                      ! never used) and a short step to maximize the stability of the
                                      ! algorithm


 ierror = 0 ! initialize error flag

 old_ubound = 0.0_double
 upper_bound = old_ubound
! scale =  sign(.5_double,c_par) - before we used a fixed value, but this occasionally turned out to be unstable
! Take the slowest scale
 scale =  sign(1.0 / (1.0 + min(abs(c_par),  c_almost_one) ),  c_par)

 call two_cat_ll(d_a_par, c_par,cut_neg,cut_pos, num_neg,num_pos, old_ubound, old_ll, ierror)
 do  iter = 1, max_iter
      upper_bound =  upper_bound + scale
      call two_cat_ll(d_a_par, c_par,cut_neg,cut_pos, num_neg,num_pos, upper_bound, ll, ierror)
      if(ierror /= 0) then
          err_msg = "ERROR:: set_infinite_upper_bound: 2 category LL lost meaning "
          return
      endif
      if(ll < old_ll) then
            exit ! if the increment is negative, we passed the maximum, which means
                 ! that we have found the upper bound
      endif
     old_ubound = upper_bound
     old_ll = ll
  enddo
 ! verify that convergence was obtained
 if(iter > max_iter) then
          err_msg = "ERROR:: set_infinite_upper_bound: did not find upper_bound "
          ierror = 1
 endif

 end subroutine set_infinite_upper_bound

!--------------------------------------------------
!--------------------------------------------------
  subroutine two_cat_ll(d_a_par, c_par, cut_neg, cut_pos, num_neg,num_pos, vc_cutoff, &
                           log_like, ierror)
! PURPOSE : computes the likelihood function under the assumption that this is the only
!           cutoff and so that it is separating the data in only 2 categories

! NOTE:    the procedure simply reports to the caller whether the function was deep in the tails
!          of the probability distribution for some probability. This could be improved upon by
!          computing the logarithm of the probability directly given that asymptotic forms are
!          of the form exp(-Z^2)/Sqrt(2*Pi) abs(Z) [ or sum of 2 such terms),
!          so their logarithm can be computed analytically and it is stable for a larger range
! Warning:: assumes that its parameters are meaningful
! Warning:: error flag needs to be improved, by declaring at which stage is the failure
!           happening
 implicit none

 real(kind=double),intent(IN):: d_a_par,c_par
 integer, intent(IN) :: cut_neg, cut_pos ! number of act. neg. and act. pos above this cutoff
 integer, intent(IN) :: num_neg, num_pos ! total number of act. neg. and act. pos
 real(kind=double),intent(IN):: vc_cutoff! cutoff value in the vc space where the log_like with
                                         ! respect to the cutoff has to be computed

 real(kind=double),intent(out):: log_like
 integer,intent(out):: ierror ! error flag

 real(kind=double):: p, one_minus_p  ! probability, for the actually negative (for cat and complement)
 real(kind=double):: q, one_minus_q  ! probability, for the actually positive
 real(kind=double):: smallest_prob

 ! Evaluate the p (and  1-p ) ,  q_j (and  1-q ). The 1-p (or q) terms are needed for numerical
 ! stability

 call fpf_PBM(d_a_par, c_par, vc_cutoff, p, one_minus_p)
 call tpf_PBM(d_a_par, c_par, vc_cutoff, q, one_minus_q)

 smallest_prob = TINY(1.0_double)

 ! Make sure that there are no terms equal to zero, which would make the logarithm overflow
 if( p < smallest_prob) then
  p = smallest_prob
 elseif( one_minus_p < smallest_prob) then
  one_minus_p =  smallest_prob
 endif

if( q < smallest_prob) then
  q = smallest_prob
 elseif( one_minus_q < smallest_prob) then
  one_minus_q =  smallest_prob
 endif


 log_like =  cut_neg * log(p) + ( num_neg - cut_neg) * log(one_minus_p) + &
             cut_pos * log(q) + ( num_pos - cut_pos) * log(one_minus_q)


 if( check_number(log_like) ) then
       ierror = 0
 else
       ierror = 1
 endif

!-------------------------------------------------------
 end subroutine two_cat_ll
!-------------------------------------------------------

 logical function test_parameters(d_a_par, c_par)
!-------------------------------------------------------

! DESCRIPTION: CHECK IF PARAMETERS ARE MEANINGFUL and report back to the calling program.
! PURPOSE    : Prevent mathematically absurd or numerically unstable values for the parameters
!              from entering the cycle.
! ALGORITHM  : check the limits of d_a > 0 and abs(c) smaller than one.
! NOTE       : This algorithm used to be much more complex, but insted the complexity is left to
!              the log likelihood function (in dealing with numerical issues). this is based on the
!              assumption that the values that are numerically destabilizing can in fact be replaced
!              by other ones without affecting the result in a numerically meaningful way
 implicit none

 real(kind=double), intent(in)  :: d_a_par  ! d_a as passed by the call
 real(kind=double), intent(in)  :: c_par    ! same but for c

 if( d_a_par > 0.0_double .and. abs(c_par) < 1.0_double ) then
  test_parameters = .true.
 else
  test_parameters = .false.
 endif

!-------------------------------------------------------
 end function test_parameters
!-------------------------------------------------------


!-------------------------------------------------------
!-------------------------------------------------------
 subroutine fpf_PBM(d_a, c, v_c, fpf, one_minus_fpf)
!-------------------------------------------------------
!-------------------------------------------------------
! PURPOSE: compute the TPF value for a cutoff value (v_c), and
!          its complemt to one,  given the curve parameters (d_a, c)
!          according to the proper binormal model

!ALGORITHM: Follow the equations from Metz & Pan, but use the
!           independently computed phi and one_minus_phi to
!           ensure that the obtained values are both numerically
!           significant. When the cancellation of precisions is due
!           to using identical but opposite arguments in the phi
!           function, attempt a fix using the trapezoid rule
!           (the fix is only attempted because the cancellation
!            will appear also in the dz term
!NOTE:      the function does not attempt to change values when they are too
!           close to zero, this correction should be done by the calling program
!           according to needs
 use statistic_functions, only: norm_dist,g
 implicit none

 real(kind= double), intent(in):: d_a ! curve parameter d_a
 real(kind= double), intent(in):: c ! curve parameter c
 real(kind= double), intent(in):: v_c ! cutoff in the v space (se M&P paper)


 real(kind= double), intent(out):: fpf ! the value for the false positive fraction
 real(kind= double), optional, intent(out):: one_minus_fpf ! the value for its complement to 1
  ! it is optional because for io purposes and other situations we don't need the complement


 real(kind= double):: z1, z2 ! the 2 z's of the proper TPF phi's
                             ! see Metz & Pan paper

 real(kind= double):: phi1, one_minus_phi1
 real(kind= double):: phi2, one_minus_phi2

 real(kind=double):: p, one_minus_p

 ! compute z for the first z in the definition of the proper fpf
 z1 =   z1_fpf(d_a, c, v_c)

 ! If c is very close to 0 the proper model concide with the conventional
 ! model, so we can drop the second phi function
 if( abs(c)  < c_almost_zero ) then
     call norm_dist(z1,fpf, one_minus_fpf)
 else

     z2 = z2_fpf(d_a, c, v_c)


     ! Consider that  changing the sign of c is equivalent to change
     ! all the phi distributions with 1 - phi, i.e. z -> -z

     call norm_dist(-SIGN(1.0_double,c) * z1  ,phi1, one_minus_phi1)
     call norm_dist(-SIGN(1.0_double,c) * z2  ,phi2, one_minus_phi2)


     p =  phi1 + phi2
     ! if there are no important cancellation of significant figures
     if( one_minus_phi1 .spne. phi2) then
        one_minus_p = one_minus_phi1 - phi2
     ! Otherwise use the trapezoid rule to compute the integral.
     ! note that if c < 0 z1+z2 < 0 otherwise z1+z2 > 0 (as expected
     ! because of symmetry). Insure that z1 =! -z2 (limiting value
     ! but it can happen because of cancellation problems), which would
     ! underflow one_minus_p
    else
         one_minus_p =  .5_double * ( g(z1) + g(z2) ) * abs( z1 + z2)
    endif



     if ( c < 0) then
         fpf = p
         if(present(one_minus_fpf)) one_minus_fpf = one_minus_p
     else
         fpf = one_minus_p
         if(present(one_minus_fpf)) one_minus_fpf = p
     endif

 endif


!-------------------------------------------------------
  end subroutine fpf_PBM
!-------------------------------------------------------
!-------------------------------------------------------

!-------------------------------------------------------
!-------------------------------------------------------
 real(kind=double) pure function d_fpf_d_d_a_PBM(d_a, c, v_c)
!-------------------------------------------------------
!-------------------------------------------------------
! PURPOSE: compute the d(FPF)/d[d_a] value for a cutoff value (v_c)  according to the proper binormal model

!ALGORITHM:
 use statistic_functions, only: g
 implicit none

 real(kind= double), intent(in):: d_a ! curve parameter d_a
 real(kind= double), intent(in):: c ! curve parameter c
 real(kind= double), intent(in):: v_c ! cutoff in the v space (se M&P paper)


 real(kind= double):: z1, z2 ! the 2 z's of the proper TPF phi's
                             ! see Metz & Pan paper


 ! compute z for the first z in the definition of the proper fpf
 z1 =   z1_fpf(d_a, c, v_c)

 ! If c is very close to 0 the proper model concide with the conventional
 ! model, so we can drop the second phi function
 if( abs(c)  < c_almost_zero ) then
     d_fpf_d_d_a_PBM = - ( sqrt(1.0_double + c**2)/ 2.0_double ) * g(z1)
 else
     z2 = z2_fpf(d_a, c, v_c)
     d_fpf_d_d_a_PBM = - ( sqrt(1.0_double + c**2)/ 2.0_double ) * g(z1)  &
                          + ( sqrt(1.0_double + c**2)/ (2.0_double * c) ) * g(z2)
 endif

!-------------------------------------------------------
 end function d_fpf_d_d_a_PBM
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
 real(kind=double) pure function d_fpf_d_c_PBM(d_a, c, v_c)
!-------------------------------------------------------
!-------------------------------------------------------
! PURPOSE: compute the d(FPF)d[c] value for a cutoff value (v_c)  according to the proper binormal model

!ALGORITHM:
 use statistic_functions, only: g
 implicit none

 real(kind= double), intent(in):: d_a ! curve parameter d_a
 real(kind= double), intent(in):: c ! curve parameter c
 real(kind= double), intent(in):: v_c ! cutoff in the v space (se M&P paper)


 real(kind= double):: z1, z2 ! the 2 z's of the proper TPF phi's
                             ! see Metz & Pan paper


 ! compute z for the first z in the definition of the proper fpf
 z1 =   z1_fpf(d_a, c, v_c)

 ! If c is very close to 0 the proper model concide with the conventional
 ! model, so we can drop the second phi function
 if( abs(c)  < c_almost_zero ) then
     d_fpf_d_c_PBM =  ( v_c - d_a * c / (2.0_double * sqrt(1.0_double + c**2) ) ) * g(z1)

 else
     z2 = z2_fpf(d_a, c, v_c)
     d_fpf_d_c_PBM =  ( v_c - d_a * c / (2.0_double * sqrt(1.0_double + c**2) ) ) * g(z1)  + &
      ( v_c - d_a*sqrt(1.0_double + c**2)/(2.0_double * c**2) + d_a/(2.0_double * sqrt(1.0_double + c**2)) )*g(z2)
 endif

!-------------------------------------------------------
 end function d_fpf_d_c_PBM
!-------------------------------------------------------
!-------------------------------------------------------

!-------------------------------------------------------
!-------------------------------------------------------
 real(kind=double) pure function d_tpf_d_d_a_PBM(d_a, c, v_c)
!-------------------------------------------------------
!-------------------------------------------------------
! PURPOSE: compute the d(TPF)/d[d_a] value for a cutoff value (v_c)  according to the proper binormal model

!ALGORITHM:
 use statistic_functions, only: g
 implicit none

 real(kind= double), intent(in):: d_a ! curve parameter d_a
 real(kind= double), intent(in):: c ! curve parameter c
 real(kind= double), intent(in):: v_c ! cutoff in the v space (se M&P paper)


 real(kind= double):: z1, z2 ! the 2 z's of the proper TPF phi's
                             ! see Metz & Pan paper


 ! compute z for the first z in the definition of the proper fpf
 z1 =   z1_tpf(d_a, c, v_c)

 ! If c is very close to 0 the proper model concide with the conventional
 ! model, so we can drop the second phi function
 if( abs(c)  < c_almost_zero ) then
     d_tpf_d_d_a_PBM =  ( sqrt(1.0_double + c**2)/ 2.0_double ) * g(z1)
 else
     z2 = z2_tpf(d_a, c, v_c)
     d_tpf_d_d_a_PBM =  ( sqrt(1.0_double + c**2)/ 2.0_double ) * g(z1) + &
                           ( sqrt(1.0_double + c**2)/ (2.0_double * c) ) * g(z2)
 endif

!-------------------------------------------------------
 end function d_tpf_d_d_a_PBM
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
 real(kind=double) pure function d_tpf_d_c_PBM(d_a, c, v_c)
!-------------------------------------------------------
!-------------------------------------------------------
! PURPOSE: compute the d(TPF)d[c] value for a cutoff value (v_c)  according to the proper binormal model

!ALGORITHM:
 use statistic_functions, only: g
 implicit none

 real(kind= double), intent(in):: d_a ! curve parameter d_a
 real(kind= double), intent(in):: c ! curve parameter c
 real(kind= double), intent(in):: v_c ! cutoff in the v space (se M&P paper)


 real(kind= double):: z1, z2 ! the 2 z's of the proper TPF phi's
                             ! see Metz & Pan paper


 ! compute z for the first z in the definition of the proper fpf
 z1 =   z1_tpf(d_a, c, v_c)

 ! If c is very close to 0 the proper model concide with the conventional
 ! model, so we can drop the second phi function
 if( abs(c)  < c_almost_zero ) then
     d_tpf_d_c_PBM = ( -v_c + d_a * c / (2.0_double * sqrt(1.0_double + c**2) ) ) * g(z1)
 else
     z2 = z2_tpf(d_a, c, v_c)
     d_tpf_d_c_PBM = ( -v_c + d_a * c / (2.0_double * sqrt(1.0_double + c**2) ) ) * g(z1)  + &
      ( -v_c - d_a*sqrt(1.0_double + c**2)/(2.0_double * c**2) + d_a/(2.0_double * sqrt(1.0_double + c**2)) )*g(z2)
 endif

!-------------------------------------------------------
 end function d_tpf_d_c_PBM
!-------------------------------------------------------
!-------------------------------------------------------


!-------------------------------------------------------
 real (kind=double) function density_vc_PBM(d_a, c, v_c, pos)
!-------------------------------------------------------
! PURPOSE: compute the density conditional on the truth of cases being positive or
!          negative, equivalent to the derivative of TPF and FPF in vc, respectively
! CREATION: June 2007 LP UCHICAGO

 use statistic_functions, only: norm_dist, g

 implicit none

 real(kind=double), intent(IN):: d_a ! curve parameter d_a
 real(kind=double), intent(IN):: c ! curve parameter c
 real(kind=double), intent(IN):: v_c ! cutoff in the v space (se M&P paper)
 integer , intent(IN) :: pos ! whether the density is for positive (1) or negative cases (0)

 real(kind=double) :: k1, k2, k3


 if(pos==1)then
     k1 = -(1.0_double + c)
     k2 =  d_a * sqrt(1.0_double + c**2) / 2.0_double
     k3 =  k2 / sign( max(abs(c), 1.0e-8_double), c)
 else
     k1 = -(1.0_double - c)
     k2 = - d_a * sqrt(1.0_double + c**2) / 2.0_double
     k3 = - k2 / sign( max(abs(c), 1.0e-8_double), c)
 endif

 density_vc_PBM = abs(k1) * ( g(k1*v_c + k2) + g(k1*v_c + k3) )


!-------------------------------------------------------
end function density_vc_PBM
!-------------------------------------------------------
!-------------------------------------------------------


!-------------------------------------------------------
 subroutine tpf_PBM(d_a, c, v_c, tpf, one_minus_tpf)
!-------------------------------------------------------
!-------------------------------------------------------
! PURPOSE: compute the TPF value for a cutoff value (v_c), and
!          its complemt to one,  given the curve parameters (d_a, c)
!          according to the proper binormal model

!ALGORITHM: Follow the equations from Metz & Pan, but use the
!           independently computed phi and one_minus_phi to
!           ensure that the obtained values are both numerically
!           significant. When the cancellation of precisions is due
!           to using identical but opposite arguments in the phi
!           function, attempt a fix using the trapezoid rule
!           (the fix is only attempted because the cancellation
!            will appear also in the dz term
!NOTE:      the function does not attempt to change values when they are too
!           close to zero, this correction should be done by the calling program
!           according to needs
 use statistic_functions, only: norm_dist,g
 implicit none

 real(kind=double), intent(IN):: d_a ! curve parameter d_a
 real(kind=double), intent(IN):: c ! curve parameter c
 real(kind=double), intent(IN):: v_c ! cutoff in the v space (se M&P paper)


 real(kind=double), intent(out):: tpf ! the value for the false positive fraction
 real(kind=double), optional, intent(out):: one_minus_tpf ! the value for its complement to 1

 real(kind=double):: z1, z2 ! the 2 z's of the proper TPF phi's
                             ! see Metz & Pan paper

 real(kind=double):: phi1, one_minus_phi1
 real(kind=double):: phi2, one_minus_phi2

 real(kind=double):: p, one_minus_p


 ! compute z for the first z in the definition of the proper tpf
 z1 =  - (1.0_double + c) * v_c + d_a * SQRT(1 + c**2) / 2.0_double


 ! If c is very close to 0 the proper model concide with the conventional
 ! model, so we can drop the second phi function
 if( abs(c)  < c_almost_zero ) then
     call norm_dist(z1,tpf, one_minus_tpf)
 else

     z2 =   - (1.0_double + c) * v_c +  d_a* SQRT(1 + c**2) / (2.0_double * c)

    ! Consider that  changing the sign of c is equivalent to change
    ! all the phi distributions with 1 - phi
    call norm_dist(-SIGN(1.0_double,c) * z1, phi1, one_minus_phi1)
    call norm_dist(-SIGN(1.0_double,c) * z2, phi2, one_minus_phi2)

    p =  phi1 + phi2

    ! there are no important cancellation of significant figures
    ! if p is different enough from one
!    if( one_minus_phi1 .spne. phi2) then
     if( abs(one_minus_phi1 - phi2) > 1.0e-7_double) then
       one_minus_p = one_minus_phi1 - phi2
    ! Otherwise use the trapezoid rule to compute the integral.
    ! note that if c < 0 z1+z2 < 0 otherwise z1+z2 > 0 (as expected
    ! because of symmetry). Insure that z1 =! -z2 (limiting value
    ! but it can happen because of cancellation problems), which would
    ! underflow one_minus_p. Note that they sum because it computes 1 - Phi(z1)  - Phi(z2)
    else
      one_minus_p =  .5_double * ( g(z1) + g(z2) ) * abs( z1 + z2)
    endif

    if ( c < 0) then
        tpf = p
        if(present(one_minus_tpf)) one_minus_tpf = one_minus_p
     else
        tpf = one_minus_p
        if(present(one_minus_tpf)) one_minus_tpf = p
    endif

 endif


!-------------------------------------------------------
  end subroutine tpf_PBM
!-------------------------------------------------------
!-------------------------------------------------------

!---------------------------------------------------------------------
!---------------------------------------------------------------------
 subroutine hessian_matrixAUC(num_normal_cases,num_abnormal_cases,num_cat,hessian, &
                           d_a_par,c_par, v_c_par)

! PURPOSE: computes the hessian matrix of the -log likelihood function but changing the coordinates
!          system from d_a, c to the gradient of the AUC(d_a, c) and its normal
! NOTE:    It should be checked. It is not used because it did not provide the numerical advantages
!          we were hoping to get (better variances). It showed that the problem isn't with the hessian
!          inversion but more with the delta method and the non quadratic nature of the proper
!          log likelihood. We left it in in case someone wanted to try to use it. It provides exactly the
!          same answers as the straight hessian.
!---------------------------------------------------------------------
!---------------------------------------------------------------------
 use statistic_functions, only: phi,g
 implicit none

 integer, intent(IN):: num_cat
 integer, intent(IN) :: num_normal_cases, num_abnormal_cases
 real(kind=double), dimension(num_cat+1,num_cat+1), intent(out):: hessian
 real(kind=double), intent(IN):: d_a_par, c_par

 real(kind=double), dimension(num_cat):: v_c_par

 integer::  num_cutoffs ! To avoid the inherent confusion due to the fact that
                        ! there are num_cat -1 cutoffs, num_cat categories and num_cat + 1
                        ! parameters (and also num_cat + 1 cutoffs if one count the
                        ! the points at infinity), we rename num_cat - 1 num_cutoffs

 real(kind=double), dimension(num_cat-1):: fpf, tpf
 real(kind=double), dimension(num_cat-1):: d_fpf_d_d_a, d_tpf_d_d_a
 real(kind=double), dimension(num_cat-1):: d_fpf_d_c  , d_tpf_d_c
 real(kind=double), dimension(num_cat-1):: d_fpf_d_v_c, d_tpf_d_v_c


 real(kind=double), dimension(num_cat):: expected_p, expected_q
 real(kind=double), dimension(num_cat):: d_expected_p_d_d_a, d_expected_q_d_d_a, dummy
 real(kind=double), dimension(num_cat):: d_expected_p_d_c, d_expected_q_d_c

 integer:: i

 real(kind=double) :: arg_fpf, arg_tpf ! arguments of the phi functions used to
                                       ! compute the fpf and tpf

 real(kind=double):: d_fpf_d_arg, d_tpf_d_arg ! derivative of the phi functions
                                   ! with respect to their argument
 real(kind=double):: one_minus_fpf, one_minus_tpf
 ! Derivatives of the arguments
 real(kind=double):: d_arg_fpf_d_d_a, d_arg_fpf_d_c, d_arg_fpf_d_v_c
 real(kind=double):: d_arg_tpf_d_d_a, d_arg_tpf_d_c, d_arg_tpf_d_v_c

 real(kind=double):: d_auc_d_c, d_auc_d_d_a, r

 fpf = 0.0_double
 tpf = 0.0_double
 hessian = 0.0_double ! initialize the covariance matrix to zero, since most the terms will be zero
                  ! anyway
 num_cutoffs = num_cat -1

! if(ABS(c_par) < c_almost_zero) then
!     d_auc_d_d_a =  1.0_double / SQRT(2.0_double) *  g( d_a_par / SQRT(2.0_double) )
!     d_auc_d_c  = 0.0_double
! else
     r  = - min( ( 1.0_double -  c_par**2 ) / ( 1.0_double +  c_par**2 ), .999999999_double)
     d_auc_d_d_a =  1.0_double / SQRT(2.0_double) * g( d_a_par / SQRT(2.0_double) ) * &
                       ( 1.0_double  - &
                         2.0_double*phi( r * d_a_par / SQRT(2.0_double - 2.0_double * r**2)  ) &
                       )

     d_auc_d_c = EXP ( - d_a_par**2 / ( 4.0_double * ( 1.0_double - r**2) ) ) &
                  / ( 2.0_double * pi *SQRT( 1.0_double - r**2))  * &
                  4.0_double * c_par / ( 1.0_double + c_par**2)**2
! endif

 ! Do loop over the cutoff values to build the the TPFs, FPFs and their
 ! derivatives. which will be user to construct the derivatives of the
 ! p_i's and q_i's.

fpfs_AND_tpfs: do i = 1, num_cutoffs

       call fpf_PBM(d_a_par, c_par, v_c_par(i), fpf(i), one_minus_fpf)
       call tpf_PBM(d_a_par, c_par, v_c_par(i), tpf(i), one_minus_tpf)

     ! Arguments of the phi function for the tpf and fpf respectively
     ! Note that the terms are first initialized with the value from the
     ! first phi function and then the other terms are added as required
     arg_fpf = - ( 1.0_double - c_par) * v_c_par(i) - &
                (d_a_par / 2.0_double) * SQRT( 1.0_double + c_par**2)
     arg_tpf = - ( 1.0_double + c_par) * v_c_par(i)  + &
                (d_a_par / 2.0_double) * SQRT( 1.0_double + c_par**2)



     ! build the derivatives of the FPF's and tpf's: first derivative in
     !                        the argument of the phi functions
     d_fpf_d_arg = g(arg_fpf)
     d_tpf_d_arg = g(arg_tpf)

     !                   derivatives of argument of phi in the arguments
     d_arg_fpf_d_d_a = - SQRT( 1.0_double + c_par**2) / 2.0_double
     d_arg_fpf_d_c   = + v_c_par(i) - d_a_par * c_par / ( 2.0_double * SQRT( 1.0_double + c_par**2) )
     d_arg_fpf_d_v_c = - ( 1.0_double - c_par)


     d_arg_tpf_d_d_a = - d_arg_fpf_d_d_a
     d_arg_tpf_d_c   = - d_arg_fpf_d_c
     d_arg_tpf_d_v_c = - ( 1.0_double + c_par)

     d_fpf_d_d_a(i) = d_fpf_d_arg * d_arg_fpf_d_d_a
     d_fpf_d_c(i)   = d_fpf_d_arg * d_arg_fpf_d_c
     d_fpf_d_v_c(i) = d_fpf_d_arg * d_arg_fpf_d_v_c ! only one derivative is different from zero
                                                    ! since it depends upon only one cutoff

     d_tpf_d_d_a(i) = d_tpf_d_arg * d_arg_tpf_d_d_a
     d_tpf_d_c(i)   = d_tpf_d_arg * d_arg_tpf_d_c
     d_tpf_d_v_c(i) = d_tpf_d_arg * d_arg_tpf_d_v_c ! " "


! compute the terms from the second phi and the Heavyside
! function. Note that the heavyside has zero derivative, apart
! from for c=0, where it cancels out with the second phi
! function anyway
  if (ABS(c_par) > c_almost_zero) then

     ! Arguments of the phi function for the tpf and fpf respectively
     ! Note that the terms are first initialized with the value from the
     ! first phi function and then the other terms are added as required
     arg_fpf = -  ( 1.0_double - c_par)*v_c_par(i) +     &
                  (d_a_par / ( c_par * 2.0_double) ) * SQRT( 1.0_double + c_par**2)
     arg_tpf = -  ( 1.0_double + c_par)*v_c_par(i) +     &
                  (d_a_par/( 2.0_double * c_par)) * SQRT( 1.0_double + c_par**2)


     ! build the derivatives of the FPF's and TPF's: first derivative in
     !                        the argument of the phi functions
     d_fpf_d_arg = g(arg_fpf)
     d_tpf_d_arg = g(arg_tpf)

     !                   derivatives of argument of phi in the arguments
     d_arg_fpf_d_d_a = + SQRT( 1.0_double + c_par**2) / ( 2.0_double * c_par)
     d_arg_fpf_d_c   = + v_c_par(i) - d_a_par * SQRT( 1.0_double + c_par**2)  / ( 2.0_double * c_par**2) &
                       + d_a_par / ( 2.0_double * SQRT( 1.0_double + c_par**2) )
     d_arg_fpf_d_v_c = - ( 1.0_double - c_par)

     d_arg_tpf_d_d_a = + d_arg_fpf_d_d_a
     d_arg_tpf_d_c   = - v_c_par(i) - d_a_par * SQRT( 1.0_double + c_par**2)  / ( 2.0_double * c_par**2) &
                       + d_a_par / ( 2.0_double * SQRT( 1.0_double + c_par**2) )
     d_arg_tpf_d_v_c = - ( 1.0_double + c_par)


     d_fpf_d_d_a(i) =  d_fpf_d_d_a(i) + d_fpf_d_arg * d_arg_fpf_d_d_a
     d_fpf_d_c(i)   =  d_fpf_d_c(i)   + d_fpf_d_arg * d_arg_fpf_d_c
     d_fpf_d_v_c(i) =  d_fpf_d_v_c(i) + d_fpf_d_arg * d_arg_fpf_d_v_c ! only one
                          !derivative is different from zero since it depends upon only one cutoff

     d_tpf_d_d_a(i) =  d_tpf_d_d_a(i) + d_tpf_d_arg * d_arg_tpf_d_d_a
     d_tpf_d_c(i)   =  d_tpf_d_c(i)   + d_tpf_d_arg * d_arg_tpf_d_c
     d_tpf_d_v_c(i) =  d_tpf_d_v_c(i) + d_tpf_d_arg * d_arg_tpf_d_v_c ! " "

   end if

 end do fpfs_AND_tpfs

! Use the FPF's and TPF's to build the p_i & q_i (differences of FPF and TPF respectively
! which are the expected values for the population in each category.


 expected_p(1)           = 1.0_double         - fpf(1)
 expected_p(2:num_cat-1) = fpf(1:num_cat-2)   - fpf(2:num_cat-1)
 expected_p(num_cat)     = fpf(num_cat-1)  ! last term is zero

 expected_q(1)           = 1.0_double           - tpf(1)
 expected_q(2:num_cat-1) = tpf(1:num_cat-2)     - tpf(2:num_cat-1)
 expected_q(num_cat)     = tpf(num_cat-1)   ! last term is zero


! Use the d{FPF}'s and d{TPF}'s to build the d{p}_i & d{q}_i (differences of FPF and TPF respectively
! which are the expected values for the population in each category.
! NOTE: the derivatives for the cutoffs have Kronecker delta selection rules (d(i,i) & d(i+1) )
!       and as such needs to be addressed differently
! NOTE: The derivatives in d_a accept a scaling factor to condition the hessian matrix. see definition
!       of the variable.

 d_expected_p_d_d_a(1)           = - d_fpf_d_d_a(1)
 d_expected_p_d_d_a(2:num_cat-1) = + d_fpf_d_d_a(1:num_cat-2) - d_fpf_d_d_a(2:num_cat-1)
 d_expected_p_d_d_a(num_cat)     = + d_fpf_d_d_a(num_cat -1)

 d_expected_p_d_c(1)             = - d_fpf_d_c(1)
 d_expected_p_d_c(2:num_cat-1)   = + d_fpf_d_c(1:num_cat-2) - d_fpf_d_c(2:num_cat-1)
 d_expected_p_d_c(num_cat)       = + d_fpf_d_c(num_cat -1)

 d_expected_q_d_d_a(1)           = - d_tpf_d_d_a(1)
 d_expected_q_d_d_a(2:num_cat-1) = + d_tpf_d_d_a(1:num_cat-2) - d_tpf_d_d_a(2:num_cat-1)
 d_expected_q_d_d_a(num_cat)     = + d_tpf_d_d_a(num_cat -1)

 d_expected_q_d_c(1)             = - d_tpf_d_c(1)
 d_expected_q_d_c(2:num_cat-1)   = + d_tpf_d_c(1:num_cat-2) - d_tpf_d_c(2:num_cat-1)
 d_expected_q_d_c(num_cat)       = + d_tpf_d_c(num_cat -1)



! CHANGE THE DERIVATIVES IN D_A & C, into D_AUC, and D_normal_AUC, but KEEP THE NAMES
! (I KNOW, VERY CONFUSING)
! d_auc_d_d_a = d_auc_d_d_a
! d_auc_d_c = d_auc_d_c

 dummy =  d_expected_p_d_d_a / d_auc_d_d_a + d_expected_p_d_c / d_auc_d_c
 d_expected_p_d_c =  - d_expected_p_d_d_a * d_auc_d_c + d_expected_p_d_c * d_auc_d_d_a
 d_expected_p_d_d_a = dummy


 dummy =  d_expected_q_d_d_a / d_auc_d_d_a + d_expected_q_d_c / d_auc_d_c
 d_expected_q_d_c =  -d_expected_q_d_d_a * d_auc_d_c + d_expected_q_d_c * d_auc_d_d_a
 d_expected_q_d_d_a = dummy


 ! Check that terms aren't too small or negative T
 ! WARNING THIS NEEDS TO BE ADDRESSED BETTER, SINCE
 ! IT IMPLIES THAT BOUNDARIES are bad!
 ! Please notice that the way the output is done is incompatible with the rest of the code. We did not want to
 ! spend time upgrading a function that we are not using. If you want to play with it, fix it.
  if ( ANY(expected_p  <= 0.0_double )  .OR.  ANY(expected_q < 0.0_double ) ) then
 !    write(log_file,*)
 !    write(log_file,*) "WARNING:: Negative value for a p or q! "
 !    write(log_file,*) "          while computing the covariance matrix! "
 !    write(log_file,*) "          VARIANCE MIGHT BE UNRELIABLE! "
 !    write(log_file,*) "          d_a: " ,d_a_par,", c: ", c_par
 !    write(log_file,*) "          DOMAIN OF DEFINITION OF cutoffs is "
 !        if (ABS(c_par) < c_almost_zero ) then
 !              write(log_file,*) '          [ - INFINITY , + INFINITY ] '
 !        elseif ( c_par < - c_almost_zero ) then
 !              write(log_file,"('           [ ',d14.6, ' + INFINITY ] ')") &
 !                   d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * c_par)
 !        else
 !              write(log_file,"('           [ - INFINITY , ',d14.6, ' ]')") &
 !                   d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * c_par)
 !        endif
 !    write(log_file,*)
 !    write(log_file,*) "          i      cut_off         p_i            q_i"
 !    do i=1,num_cat
 !      write(log_file,"(10x, i2, 3(1x,d16.8))") i, v_c_par(i), expected_p(i), expected_q(i)
 !    enddo
     do i=1,num_cat
       expected_p(i) = MAX( expected_p(i), 5.0e-10_double )
       expected_q(i) = MAX( expected_q(i), 5.0e-10_double )
     enddo
 endif

! Construct the 2nd derivative matrix at the expected value. This renders the derivatives equal to:
!  d^2(LL)/d(theta_j)  d(theta_k) =
!                    - NUM_NEG_CASES * Sum_i [ 1 / p_i * d(p_i)/d(theta_j) * d(p_i)/d(theta_k)  ] +
!                    - NUM_POS_CASES * Sum_i [ 1 / q_i * d(q_i)/d(theta_j) * d(q_i)/d(theta_k)  ] +

!***** respect D_A_PAR and D_A_PAR **********

 hessian(1,1) = - num_normal_cases   * DOT_PRODUCT( d_expected_p_d_d_a/expected_p , d_expected_p_d_d_a )  &
            - num_abnormal_cases * DOT_PRODUCT( d_expected_q_d_d_a/expected_q , d_expected_q_d_d_a )

!
!***** respect C_PAR and C_PAR **********
!

 hessian(2,2) = - num_normal_cases   * DOT_PRODUCT( d_expected_p_d_c/expected_p , d_expected_p_d_c )  &
            - num_abnormal_cases * DOT_PRODUCT( d_expected_q_d_c/expected_q , d_expected_q_d_c )

!
!**** respect D_A_PAR and C_PAR **********
!

 hessian(1,2) = - num_normal_cases   * DOT_PRODUCT( d_expected_p_d_d_a/expected_p , d_expected_p_d_c )  &
            - num_abnormal_cases * DOT_PRODUCT( d_expected_q_d_d_a/expected_q , d_expected_q_d_c )

!
!**** respect D_A_PAR and v_c l **********
 forall(i = 1:num_cutoffs)
     hessian(1,2+i) = - num_normal_cases * (- d_expected_p_d_d_a(i)   * d_fpf_d_v_c(i) /expected_p(i)     &
                                        + d_expected_p_d_d_a(i+1) * d_fpf_d_v_c(i) /expected_p(i+1) ) &
                     -num_abnormal_cases *(- d_expected_q_d_d_a(i)   * d_tpf_d_v_c(i) /expected_q(i)     &
                                        + d_expected_q_d_d_a(i+1) * d_tpf_d_v_c(i) /expected_q(i+1)  )
 end forall


!**** respect D_C_PAR and v_c l **********
 forall(i = 1:num_cutoffs)
      hessian(2,2+i) = -num_normal_cases * ( -d_expected_p_d_c(i)   * d_fpf_d_v_c(i) /expected_p(i)     &
                                         + d_expected_p_d_c(i+1) * d_fpf_d_v_c(i) /expected_p(i+1) ) &
                - num_abnormal_cases * ( -d_expected_q_d_c(i)   * d_tpf_d_v_c(i) /expected_q(i)     &
                                         + d_expected_q_d_c(i+1) * d_tpf_d_v_c(i) /expected_q(i+1) )
 end forall

!
!**** respect v_c l and v_c l **********
 forall(i = 1:num_cutoffs)
        hessian(2+i,2+i) =  - num_normal_cases *  d_fpf_d_v_c(i) **2     *            &
                                ( + 1.0_double /expected_p(i) + 1.0_double /expected_p(i+1) )   &
                        - num_abnormal_cases *  d_tpf_d_v_c(i) **2     * &
                                ( + 1.0_double /expected_q(i) + 1.0_double /expected_q(i+1) )
 end forall

!
!**** respect v_c l and v_c m **********
 forall(i = 1:num_cutoffs-1)
       hessian(2+i,3+i) = - num_normal_cases *  d_fpf_d_v_c(i) * (- d_fpf_d_v_c(i+1) ) / expected_p(i+1) &
                      - num_abnormal_cases *  d_tpf_d_v_c(i) * (- d_tpf_d_v_c(i+1) ) / expected_q(i+1)
 end forall
!
!**** symmetrize the hessian matrix *******************

 do i=2,num_cat+1
      hessian(i,1:i-1) =  hessian(1:i-1,i)
 enddo

      hessian(1,:) =  hessian(1,:) * .5_double
      hessian(:,1) =  hessian(:,1) * .5_double
      hessian(2,:) =  hessian(2,:) !* 1.0e6_double
      hessian(:,2) =  hessian(:,2) !* 1.0e6_double

!------------------------------------------------------------------------------
  end subroutine hessian_matrixAUC
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
 real(kind=double) elemental function approx_var_auc(auc, mn, ms)
!------------------------------------------------------------------------------
! PURPOSE: compute the variance of the area using Charles Metz's rule of thumb rule.
implicit none
real(kind=double), intent(in) :: auc ! area  under the curve
integer, intent(in):: mn ! number of negative cases
integer, intent(in):: ms ! number of positive cases

  approx_var_auc =  auc * ( 1.0_double - auc) *  &
                   0.25_double * ( 1.0_double / mn +  1.0_double / ms )  !

!------------------------------------------------------------------------------
 end  function approx_var_AUC
!------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
 pure subroutine var_auc_PBM(d_a_par,c_par,var_d_a,var_c,cov_d_a_c, var_auc, ierror)
!---------------------------------------------------------------------------------
! Derive the Standard error of the area (Variance) from the variances and
! covariances of the parameters d_a & c  using the
! second order approximation (see Papoulis).
! Is starts from the equation:
!
! Area = Phi( d_a/ Sqrt[2]) + 2 F (- d_a/Sqrt[2], 0; - (1 - c^2)/(1 + c^2) )
!
! Then uses is first derivatives to compute:
! SE(AREA) = (D{AREA}/D{d_a})^2  VAR(d_a) + 2  (D{AREA}/D{d_a})*(D{AREA}/D{c}) * COV(d_a,c) +
!            (D{AREA}/D{c})^2 * VAR(c) + o
! CHANGES: LP April 30 2009 -> changed ierror = 2 to ierror = bad_input. Added also checkes on input values for variances

 use statistic_functions, only : phi,g
 implicit none
 ! parameters d_a, and c
 real(kind=double),  intent(IN):: d_a_par,c_par
 ! variance d_a, variance c, covariance d_a,c
 real(kind=double),  intent(IN):: var_d_a,var_c,cov_d_a_c
 real(kind=double), intent(OUT):: var_auc
 integer,           intent(OUT):: ierror

 Real(kind=double) d_auc_d_d_a, d_auc_d_c

 real(kind=double) r

! varify that input values are meaningful
if( d_a_par < 0.0_double .or. abs(c_par) > 1.0_double .or. &
    var_d_a < 0.0_double .or. var_c < 0.0_double .or. cov_d_a_c**2 > var_d_a * var_c ) then
        ierror = bad_input
        return
else
        ierror = 0
endif


 ! Clearly in this case the result is independent of c and the formula for d_a
 ! is also simplified since the term which depends upon the bivariate binormal
 ! distribution is zero
 if(ABS(c_par) < c_almost_zero) then
       d_auc_d_d_a =  1.0_double / SQRT(2.0_double) *  g( d_a_par / SQRT(2.0_double) )
       var_auc =  d_auc_d_d_a**2 * var_d_a
 else
       r  = - min( ( 1.0_double -  c_par**2 ) / ( 1.0_double +  c_par**2 ), .999999999_double)
       d_auc_d_d_a =  1.0_double / SQRT(2.0_double) * g( d_a_par / SQRT(2.0_double) ) * &
                       ( 1.0_double  - &
                         2.0_double * phi( r * d_a_par / SQRT(2.0_double - 2.0_double * r**2)  ) &
                       )

       d_auc_d_c = EXP ( - d_a_par**2 / ( 4.0_double * ( 1.0_double - r**2) ) ) &
                  / ( pi *SQRT( (1.0_double - r**2)  ))  * &
                  4.0_double * c_par / ( 1.0_double + c_par**2)**2
       var_auc =  d_auc_d_d_a**2 * var_d_a   +   d_auc_d_c**2 * var_c  + &
                2.0_double * d_auc_d_d_a * d_auc_d_c * cov_d_a_c
 endif

!---------------------------------------------------------------------------------
 end subroutine var_auc_PBM
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
 real(kind=double) pure function z1_fpf(d_a,c,v_c)
!---------------------------------------------------------------------------------
 ! compute the argument of the first phi function used in the
 ! proper model to compute the fpf
 implicit none
 real(kind= double), intent(IN):: d_a ! curve parameter d_a
 real(kind= double), intent(IN):: c ! curve parameter c
 real(kind= double), intent(IN):: v_c ! cutoff in the v space (se M&P paper)

 z1_fpf =  - (1.0_double - c) * v_c - d_a * SQRT(1.0_double + c**2) / 2.0_double

!---------------------------------------------------------------------------------
 end function z1_fpf
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
 real(kind=double) pure function z2_fpf(d_a,c,v_c)
!---------------------------------------------------------------------------------
 ! compute the argument of the second phi function used in the
 ! proper model to compute the fpf
 implicit none
 real(kind= double), intent(IN):: d_a ! curve parameter d_a
 real(kind= double), intent(IN):: c ! curve parameter c
 real(kind= double), intent(IN):: v_c ! cutoff in the v space (se M&P paper)

 z2_fpf =  - (1.0_double - c) * v_c + d_a * SQRT(1.0_double + c**2) / ( 2.0_double * sign( max(.00000000001_double, abs(c)),c))

!---------------------------------------------------------------------------------
 end function z2_fpf
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
 real(kind=double) pure function z1_tpf(d_a,c,v_c)
!---------------------------------------------------------------------------------
 ! compute the argument of the first phi function used in the
 ! proper model to compute the fpf
 implicit none
 real(kind= double), intent(IN):: d_a ! curve parameter d_a
 real(kind= double), intent(IN):: c ! curve parameter c
 real(kind= double), intent(IN):: v_c ! cutoff in the v space (se M&P paper)

 z1_tpf =  - (1.0_double + c) * v_c + d_a * SQRT(1.0_double + c**2) / 2.0_double

!---------------------------------------------------------------------------------
 end function z1_tpf
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
 real(kind=double) pure function z2_tpf(d_a,c,v_c)
!---------------------------------------------------------------------------------
 ! compute the argument of the second phi function used in the
 ! proper model to compute the fpf
 implicit none
 real(kind= double), intent(IN):: d_a ! curve parameter d_a
 real(kind= double), intent(IN):: c ! curve parameter c
 real(kind= double), intent(IN):: v_c ! cutoff in the v space (se M&P paper)

 z2_tpf =  - (1.0_double + c) * v_c + d_a * SQRT(1.0_double + c**2) / ( 2.0_double * sign( max(.00000000001_double, abs(c)),c))

!---------------------------------------------------------------------------------
 end function z2_tpf
!---------------------------------------------------------------------------------

 end module proproc_functions
