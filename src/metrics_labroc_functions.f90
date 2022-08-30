! contains functions which are specific to LABROC (LABROC4 for be precise implementation of the
! conventional binormal model (CvBM, previous implementations involve ROCFIT, RSCORE and
! some other less known programs. This model is based on two normal distributions
! One centered at 0 with stadard deviation 1 (for the normal cases) and one centered
! at a/b with standard deviation 1/b.
  !

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

! NOTES: the Conventional binormal model will occasionally be called CvBM. This is done
!        to prevent confusion with the CBM (contaminated binormal model)

! MODULE WAS CREATED IN FEBRUARY 2009 BY LORENZO L PESCE AT THE UNIVERSITY OF CHICAGO

 module labroc_functions
! NOTE that only the constants for actually negative and actually positive
! cases are loaded from the module array_dimensions. The only is used to that
! a strick control of data and methods is enforced
use data_types, only: double, operator (.speq.), operator (.spne.), check_number
use computation_constants
use array_dimensions, only: act_neg, act_pos ! where array related quantities are
                                ! variable to be eliminated that has to do with the number of
                                ! parameters (XP is such an awful programmer that nothing makes
                                ! real sense in his junk)
use io
use error_flags


! GENERAL SETTING FOR THE WHOLE MODULE
implicit none

private ! All the variables and functions that are not set as public
        ! are not acessible


! PUBLIC PROCEDURES / VALUES
public double ! Note that we are passing this value, which is inherited from data_types
public line_length ! Note that we are passing this value, which is inherited from io
! PUBLIC PROCEDURES
public cvbmroc_mle ! main proper model procedure
public tpf_CVBM ! find TPF, knowing a, b, fpf
public fpf_CVBM ! find FPF, knowing a, b, tpf
public fpf_find_tpf_CVBM ! find TPF, knowing a, b, fpf
public tpf_find_fpf_CVBM ! find FPF, knowing a, b, tpf
public var_tpf_find_fpf_CVBM ! error of FPF, knowing a, b, tpf
public var_fpf_find_tpf_CVBM ! error of TPF, knowing a, b, fpf
public partialauc_CVBM, var_partialauc_CVBM
public auc_CVBM, var_auc_CVBM ! compute the area under the curve
public  least_squares_a_b! least squares for initial estimates, for the Conventional binormal model

! Export the meaning of the different values that the fit error flag can have
public bad_input          ! categorical input has unacceptable values
public fit_OK             ! fit was successful
public fit_OK_pseudo      ! fit was successful, variances computed using pseudoinverse
public fit_failed         ! fit failed in MLE estimation
public fit_undetermined   ! Note enough data: there is only one point, and it is inside the ROC space
public fit_perfect        ! positives and negatives are perfectly separated, it is more of a warning
public fit_fail_init_est  ! initial estimates did not converge
public fit_fail_variances ! estimates of variances did not converge
public fit_fail_var_auc_small! estimates of variances did not converge because area is too small (see below)
public fit_degenerate     ! the fit is degenerate
public fit_perverse     ! the fit is degenerate

contains

!---------------------------------------------------------------------------------------------------
subroutine cvbmroc_mle(mn, ms, num_categ, catn_in, cats_in, idebug, &
                       a_par, b_par, auc, variance_auc, cutoffs_out, log_like, ierror, &
                       cov_out, hessian_out)
! PURPOSE: compute the MLE estimate using the conventional binormal model starting from categorial data
!          the program assumes that positivity (higher likelihood of having a positive reading) is
!          for categories with larger indices. If this is not the case, the categorical data need to
!          be reverse before feeding them to this subroutine.
!          It is essentially designed after the proper binormal model fitting algorithm.
!          Originally, this algorithm was called LABROC4, but we decided to split the discretization step
!          (which we call LABROC4) from the fitting of the conventional binormal model. Historically, the
!          two were developed intertwined and the joined algorithm was called LABROC4 (and was incomprehensible)
!ALGORITHM:
!          call compute_MLE that returns the MLE estimate with the hessian matrix, that is used to
!          produce the standard errors using the Kramer-Rao bound.
!          There could be multiple maxima (mostly for very unbablanced datasets -- from proproc simulations),
!          the algorithms first attempts using 3 points spread in the a, b space, then if more than  1 maxima are found,
!          it takes the one with the largest log-likelihood. It is possible that there are more than 2 maxima, but
!          from proproc simulations it does not appear to be a serious issue.

! HISTORY: ! Created August 5th 2009
use debugging ! module containing the procedures to log the debugging information on the run
use problem_data ! where the ROC basic categorical information is stored.
use statistic_functions, ONLY : compute_zdev_plus
use categorization, ONLY:  LABROC4_Collapser

! Subroutine arguments
integer, intent(in):: mn ! number of actually negative cases
integer, intent(in):: ms ! number of actually positive cases
integer, intent(in) :: num_categ ! Number of categories as created by catgrz
integer, dimension(num_categ), intent(in):: catn_in, cats_in ! arrays containing categorical data
integer, intent(in) :: idebug               ! 0 = no debug; 1  = debug
real(kind = double),intent(out)                        :: a_par, b_par ! MLE of the parameters
real(kind = double),intent(out)                        :: auc ! AUC, area under the curve also known as Az (only for this model)
real(kind = double),intent(out)                        :: variance_auc ! estimated variance of AUC

real(kind = double), dimension(num_categ-1), intent(out) :: cutoffs_out  ! cutoff parameter values at the maximum found
real(kind = double), intent(out)                       :: log_like ! value of the log likelihood function at the final point
integer, intent(out)                                   :: ierror ! Error flag about the MLE fit
! Note that the error values are set in this routine, or initialize_a_b so if their numbers have
! to be changed, they have to be changed here, the rest of the subroutines use their own numbering
! specific per routine. Look above for the different meaning. Not only 0 is successful fit
real(kind = double), dimension(num_categ+1,num_categ+1), intent(out)  :: cov_out ! these are used because the number
real(kind = double), dimension(num_categ+1,num_categ+1), intent(out), optional  :: hessian_out ! of categories can be
     ! different inside  because of collapsing. Note that only the reduced hessian will be returned, the
     ! rest will be set to garbage

! Internal variables
integer:: i,j ! loop counter for logging warnings about cutoffs
character(len = line_length):: msg    ! character string used to display state information
character(len = line_length):: err_msg    ! character string used to display state information
! The algorithms first choses 3 points, to span the configuration space for multiple maxima, then
! if more than one maximum is found, the largest is picked.
integer, parameter:: num_init_pts = 3 ! Number of initial estimates
integer, dimension(num_init_pts)  :: pts_error ! whether a specific fit has failed
logical, dimension(num_init_pts) :: pts_OK_mask  ! Whether the MLE succeeded or failed for this initial estimate
real(kind = double), dimension(num_init_pts)  :: a_par_array, b_par_array ! MLE of the parameters
real(kind=double), dimension(num_init_pts) ::  log_like_array ! the log like for the 2 possible initial estimates

! Some variables used if multiple maxima are possible
real(kind = double) :: a_min, a_max, b_min, b_max ! judging from simulations with proproc, we decided to simply
! flag this even, but we are not trying to determine where there are in fact multiple maxima and whether we have
! missed the largest. We just take the largest log-likelihood

integer, dimension(num_categ):: new_cat_index ! if the data needs to be collapsed, the list of where to
                                              ! put the old data
integer:: auc_error ! error return by AUC computation routines, as of 12/31/2009 ignored.
integer:: icat ! loop counter

! These are the internal variables, defined using num_cat, the internal number of categories after
! collapsing
real(kind = double), allocatable, dimension(:,:)  :: cov ! internal values, can be smaller
real(kind = double), allocatable, dimension(:,:)  :: hessian ! if collapsed
real(kind = double), allocatable, dimension(:)    :: cutoffs  ! cutoff parameter values at the maximum found
real(kind = double), allocatable, dimension(:, : ) :: cutoffs_array  ! MLE cutoff parameter values

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
        write(msg,*) "Categorical data given to cvbmroc is unacceptable"
        call DisplayState(msg)
        write (msg,*) "CATEGORICAL DATA USED IN THE MLE, total cat = ", num_categ
        call DisplayState(msg)
        do i = 1, num_categ
           write (msg,*) catn_in(i), cats_in(i)
           call DisplayState(msg)
        enddo
        write (msg,*) "END OF CATEGORICAL DATA"
        call DisplayState(msg)
    endif
    ! Mathematically it should be undetermined, but we put an impossible value
    ! to prevent interfacing problems from arising
    if (present(hessian_out))  hessian_out = -666.0_double
    a_par           = -666.0_double
    b_par           = -666.0_double
    auc             = -666.0_double
    variance_auc    = -666.0_double
    cutoffs_out     = -666.0_double
    log_like        = -666.0_double
    cov_out         = -666.0_double
    return
endif
! COLLAPSE (IF NEEDED) AND LOAD THE COLLAPSED DATA INTO NEW ARRAYS
! Prepare the variables that will be stored in the moduli from file "other_modules.f90"
! Since those are in a module while these are explicit, we seem to need to double the
! names at this level. The moduli are used to work around the various optimizers (initial estimates
! of cutoffs and full MLE). NOTE: all the calls that don't have to "pass-through" an optimizer are
! made explicit, even if this might seem confusing, I prefer explicit calls when I can do them.

call LABROC4_Collapser(catn_in, cats_in, num_categ, idebug, new_cat_index)

num_cat = new_cat_index(num_categ) ! the total number of collapsed categories is the category assigned to the
       ! last category
allocate(catn(num_cat),  cats(num_cat))
allocate(cov(num_cat + 1, num_cat + 1) )
allocate(hessian(num_cat + 1, num_cat + 1) )
allocate(cutoffs(num_cat - 1) )
allocate(cutoffs_array(num_cat - 1,num_init_pts) )

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
pts_OK_mask =.false. ! Initialize array of error messages to failed

! get 3 initial estimates for d_a and c. Notice that initialize_d_a_c always returns 3, so num_init_pts isn't a
! variable that can be changed here unless we change that procedure as well
call initialize_a_b(mn, ms, num_cat, catn, cats, idebug, a_par_array(1:num_init_pts), b_par_array(1:num_init_pts), &
                    ierror, err_msg)

! If the initial estimates were successful, then we compute the MLE, if not, we skip to the final wrap up section
initial_est_OK: if(ierror == fit_OK) then

 initial_estimates_loop: do j = 1, num_init_pts ! Loop over the num_init_pts initial estimates
    ! Load the jth initial estimate
    a_par = a_par_array(j)
    b_par = b_par_array(j)

    call initialize_cutoffs(a_par,b_par, mn, ms, num_cat, catn, cats, cutoffs, idebug, ierror, err_msg)

    initial_est_cutoffs: if(ierror /= 0) then ! Initial estimates of cutoffs failed
             ! report a failure in the initialization of the cutoffs
             ierror = fit_fail_init_est ! failed initial estimates
             err_msg = " initialize_cutoffs, " // err_msg(:line_length - 21)
             pts_error(j) = fit_fail_init_est ! this specific initial estimate failed, we can't do MLE
    else ! Initial estimates of cutoff successful
             ! Check if there are any 2 cutoffs which have almost identical values
             ! It should not happen because initialize_cutoffs should prevent it
             if(idebug == 1 ) then
                 if( any(cutoffs(2:num_cat-1) .speq. cutoffs(1:num_cat-2))  ) then
                         call DisplayState( ' WARNING:: initialize_cutoffs: Some cutoffs almost identical')
                         write(msg,"('          for a & b equal to::',2(1x,d14.6))") a_par, b_par
                         call DisplayState(msg)
                         do i = 1, num_cat-1
                                        write(msg,"(i3,1(1x,e22.15))") i, cutoffs(i)
                                        call DisplayState(msg)
                        enddo
                  endif
              endif
             ! Since we checked both the initial estimates of the parameters and of the cutoffs, we can do MLE
             call compute_MLE(mn, ms, num_cat, catn, cats, idebug, a_par, b_par,  &
                              cutoffs, log_like, hessian, ierror, err_msg)
             if(ierror == 0) then
                 ierror = fit_OK
                 pts_error(j) = fit_OK ! this specific value failed
                 pts_OK_mask(j) = .true.
              else
                 ierror = fit_failed
                 pts_error(j) = fit_failed ! this specific value failed
                 pts_OK_mask(j) = .false.
                 if(idebug == 1) then
                    write(msg,*) "WARNING: Optimization for initial estimate ", j , " failed "
                    call DisplayState(msg)
                    call DisplayState(err_msg)
                 endif
             endif
    endif initial_est_cutoffs
    ! Log the current estimate, starting with the value of the initial estimates
    if(idebug == 1) then
       write(msg, *) "For initial estimates:", a_par_array(j), b_par_array(j)
       call DisplayState(msg)
    endif

    ! replace the initial estimates with the current MLE estimates
    a_par_array(j)   =  a_par
    b_par_array(j)   =  b_par
    log_like_array(j) = log_like
    cutoffs_array(1:num_cat-1,j) = cutoffs(1:num_cat-1)

    if(idebug == 1) then
       if( pts_OK_mask(j) ) then ! check if the fit worked for  this point
           write(msg, *) "Values of a, b, and Log Likelihood at Local Maximum"
           call DisplayState(msg)
           write(msg, *)  a_par, b_par, log_like
           call DisplayState(msg)
       else
           write(msg, *) "  NO FIT WAS PRODUCED because"
           call DisplayState(msg)
           write(msg, *) " ", err_msg(1:76)
           call DisplayState(msg)
       endif
    endif

! Added on may 14th 2010, to allow (if active) the program to skip fits at the first good one
! this is done because we do not expect multiple maxima to be a serious issue for CvBM
! if(pts_OK_mask(j)) exit initial_estimates_loop

 enddo initial_estimates_loop

 ! Check whether there was at least one initial
 ! estimate that produced an acceptable fit
 ! (normally at most one fails, and only for very
 !  bad datasets)
 any_fit: if( .not. any(pts_OK_mask) ) then
     ierror = fit_failed
     err_msg = "No initial estimates converged to a maximum of Log Likelihood"
 else
     ! verify that the final estimates are very similar, if they aren't, it might be that
     ! there are multiple maxima (more than found so far), we are not looking for more maxima
     ! but the users should be warned. We seek only  the MLE calculations that convereged
     ! If the maxima are rather close, we do not warn the user.
     if (&
        maxval(a_par_array(1:num_init_pts),pts_OK_mask) - minval(a_par_array(1:num_init_pts),pts_OK_mask) > .1_double &
        .or.&
        maxval(b_par_array(1:num_init_pts),pts_OK_mask) - minval(b_par_array(1:num_init_pts),pts_OK_mask) > .1_double &
     ) then
           a_min    = minval(a_par_array(1:num_init_pts),pts_OK_mask)
           b_min    = minval(b_par_array(1:num_init_pts),pts_OK_mask)

           a_min    = maxval(a_par_array(1:num_init_pts),pts_OK_mask)
           b_min    = maxval(b_par_array(1:num_init_pts),pts_OK_mask)

           if(idebug == 1) then
                write(msg,*) "ATTENTION: Log Likelihood might display multiple Maxima "
                call DisplayState(msg)
                write(msg,*) " Maximal final estimate ", a_max, b_max
                call DisplayState(msg)
                write(msg,*) " Minimal final estimate ", a_min, b_min
                call DisplayState(msg)
                write(msg,*) "ATTENTION: largest value of LL will be selected"
                call DisplayState(msg)
           endif
      endif

      ! Now select the best between all the computed points, dim = 1 is to force output to be scalar (to avoid rank mismatches)
      j = maxloc(log_like_array, dim = 1, mask=pts_OK_mask)  ! determine which one  has the largest Log Likelihood
      a_par                = a_par_array(j)
      b_par                = b_par_array(j)
      cutoffs(1:num_cat-1) = cutoffs_array(1:num_cat-1,j)

      ! Rerun from the maximum to recompute the hessian and make sure that the system is truly consistent
      ! This hardly costs any time at all, since all the parameters are already optimal, it runs in a single
      ! step -- NOTE that IERROR is set by the MLE function at the maximum
      call compute_MLE(mn, ms, num_cat, catn, cats,idebug, a_par, b_par,cutoffs, log_like, hessian, ierror, err_msg)

      ! Log the chosen value
      if(idebug == 1) then
           write(msg, *) "Final value of a, b, and Log Likelihood"
           call DisplayState(msg)
           write(msg, *)  a_par, b_par, log_like
           call DisplayState(msg)
           write(msg, *) "Final value of cutoffs:"
           call DisplayState(msg)
           do j = 1, num_cat -1
                write(msg, *)  j, cutoffs(j)
                call DisplayState(msg)
           enddo
      endif

endif any_fit ! sequence done only if at least one of the initial points worked out


endif initial_est_OK ! end of sequence done if we could produce initial estimates


! Estimate the area under the curve (AUC) and other parameters, depending upon the error flag

! Initialize arrays to nonsense, in case quantities can't be computed, applies to all values of error flags
if (present(hessian_out))  hessian_out = -666.0_double
cov_out = -666.0_double
cutoffs_out = -666.0_double


! Prepare the output and return values depending upon the error messages. Note that many of the error messages were
! defined in the initialization phase.

if(idebug == 1) then
           write(msg, *) ierror, fit_undetermined
           call DisplayState(msg)
endif

select case (ierror)
case(fit_OK) ! fit successful, compute the area
    call auc_CVBM(a_par,b_par, auc, auc_error)
    call compute_covariance_matrix(a_par, b_par, num_cat+1, idebug, hessian, cov, ierror, err_msg)
    if(auc <= .0005_double) then ! One can't compute variances using the hessian when auc is too small
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
          call  var_auc_CVBM(a_par, b_par, cov(1,1), cov(2,2), cov(1,2), variance_auc, ierror)
          ! Load reduced variance.
          cov_out(1:num_cat+1,1:num_cat+1) = cov(1:num_cat+1,1:num_cat+1)
     endif
     ! Load the reduced hessian into the hessian for output. Reduced refers to
     ! the fact that only the labroc4 corners are included. If there are extra places
     ! the values they were initialized with remains
     if (present(hessian_out)) hessian_out(1:num_cat+1,1:num_cat+1) = hessian(1:num_cat+1,1:num_cat+1)

     ! Prevents calling the reconstruction of cutoffs when it is not
     ! necessary. This steps recreates the cutoffs for the full non-collapsed data (i.e., it is possible
     ! that the input data was not reduced to the truth runs)
     if(num_cat /= num_categ) then
               call compute_orig_cutoffs(a_par, b_par, num_categ, num_cat, new_cat_index, cutoffs, &
                                catn_in, cats_in, catn, cats, idebug, cutoffs_out)
     else
               cutoffs_out = cutoffs
     endif

case(fit_perfect)
       ! there are infinite possible fits that make a perfect CvBM fit, they all have in common
       ! and infinite value of a and any value of b. We return the values chosen by the initialization
       ! routine
       auc     = 1.0_double
       variance_auc = 0.0_double ! this is the result of the MLE approach, not that any rational
                                 ! person should believe it.
       ! the other quantities are meaningless or should be considered as such
       cutoffs_out  = -666.0_double
       log_like     = -666.0_double
       cov_out      = -666.0_double
       if (present(hessian_out))  hessian_out = -666.0_double
case(fit_perverse)
       ! there are infinite possible fits that make a perfect CvBM fit, they all have in common
       ! and infinite value of a and any value of b.We return the values chosen by the initialization
       ! routine
       auc     = 0.0_double
       variance_auc = 0.0_double ! this is the result of the MLE approach, not that any rational
                                 ! person should believe it.
       ! the other quantities are meaningless or should be considered as such
       cutoffs_out  = -666.0_double
       log_like     = -666.0_double
       cov_out      = -666.0_double
       if (present(hessian_out))  hessian_out = -666.0_double
case(fit_degenerate)
       ! Mathematically it should be undetermined, but we return whatever value the initialization
       ! routine selected, which is usually some sort of reasonable average. Of course it is not very
       ! rational to use this estimate for anything more than a rough guess (check called routines to
       ! make sure
       a_par = a_par_array(1)
       b_par = b_par_array(1)
       call auc_CVBM(a_par,b_par, auc, auc_error)
       ! The remaining quantities are not very useful and should not be trusted too much (in principle
       ! they all exist and are well-defined because the likelihood function is a well-behaved function
       ! around a perfect fit) so we don't provide them. Have fun deriving them if you wish to do so:-)
       variance_auc    = -666.0_double
       cutoffs_out     = -666.0_double
       log_like        = -666.0_double
       cov_out         = -666.0_double
       if (present(hessian_out))  hessian_out = -666.0_double
case(fit_undetermined)
       ! The remaining quantities are not very useful, we don't believe the variance should be provided
       ! here because it is too unreliable.
       a_par         = -666.0_double
       b_par           = -666.0_double
       auc             = -666.0_double
       variance_auc    = -666.0_double
       cutoffs_out     = -666.0_double
       log_like        = -666.0_double
       cov_out         = -666.0_double
       if (present(hessian_out))  hessian_out = -666.0_double
case(fit_failed, fit_fail_init_est)
       ! Mathematically it should be undetermined, but we put an impossible/absurd/crazy value
       ! to prevent interfacing problems from arising
       a_par         = -666.0_double
       b_par           = -666.0_double
       auc             = -666.0_double
       variance_auc    = -666.0_double
       cutoffs_out     = -666.0_double
       log_like        = -666.0_double
       cov_out         = -666.0_double
       if (present(hessian_out))  hessian_out = -666.0_double
case default ! This happens if there is a programming error, obviously...
       ! Mathematically it should be undetermined, but we put an impossible/absurd/crazy value
       ! to prevent interfacing problems from arising
       a_par         = -666.0_double
       b_par           = -666.0_double
       auc             = -666.0_double
       variance_auc    = -666.0_double
       cutoffs_out     = -666.0_double
       log_like        = -666.0_double
       cov_out         = -666.0_double
       if (present(hessian_out))  hessian_out = -666.0_double
end select

deallocate(catn, cats) ! Deallocate the arrays in the local data modules
deallocate(cov , hessian, cutoffs, cutoffs_array)

!---------------------------------------------------------------------------------
end subroutine cvbmroc_mle
!--------------------------------------------------------------------------------
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine compute_orig_cutoffs(a_par, b_par, num_cat_orig, num_cat, new_cat_index, cutoffs, &
                                catn_orig, cats_orig, catn, cats, idebug, cutoffs_orig)
!---------------------------------------------------------------------------------------------------
! PURPOSE: reconstruct cutoffs for an array that needed to be collapsed because was sent to the MLE routine
!          without being  collapsed.
!HISTORY:  Created on MARCH 2010

use debugging
use statistic_functions, only: phi, compute_zdev_plus
implicit none

real(kind = double),intent(in)                :: a_par, b_par ! MLE of the parameters
integer, intent(in):: num_cat_orig ! number of categories in non-collapsed matrix
integer, intent(in):: num_cat ! number of categories in collapsed matrix
integer, intent(in):: idebug ! whether to log information or not
integer, dimension(num_cat_orig), intent(in):: catn_orig, cats_orig ! arrays containing categorical data non-collapsed
integer, dimension(num_cat), intent(in):: catn, cats ! arrays containing categorical data collapsed
real(kind = double), dimension(num_cat-1), intent(in) :: cutoffs  ! cutoff parameter values at the maximum found
real(kind = double), dimension(num_cat_orig-1), intent(out) :: cutoffs_orig  ! cutoff parameter values at the maximum found
integer, dimension(num_cat_orig), intent(in):: new_cat_index ! if the data needs to be collapsed, the list of where to
                                                          ! put the old data

real(kind = double) :: cutoff, tpf , fpf, one_minus, prev_fpf, prev_tpf
real(kind = double), parameter:: tol =  1.0e-10_double ! tolerance in the computation of tpf and fpf


integer:: icat ! category (cutoff) index, for the non-collapsed
integer:: num_cases ! number of cases before this category in this truth run, used if collapsing is present
integer:: old_cat_coll ! the last value of the collapsed category index
character(len = line_length):: msg    ! character string used to display state information
integer:: ierror


! Write to file debug information -- Nov 2009
 if(idebug==1) then
      write(msg,"(' ')")
      call DisplayState(msg)
      write(msg,"('CUTOFFS ARE RECONSTRUCTED (EXACTLY) FROM THE TRUTH-RUNS')")
      call DisplayState(msg)
      write(msg,"('Orig  c.psed   #orig+ /    #orig- /      FPF             TPF  ')")
      call DisplayState(msg)
      write(msg,"('cat   cat      #c.psed+    #c.psed-  ')")
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
      ! If the first category is empty, its cutoff corresponds to 1,1, that is the smallest value of the cutoff
      if( icat == 1 .and.  catn_orig(icat) == 0 .and. cats_orig(icat) == 0) then
              cutoff =  min( -6.0_double, (a_par - 6.0_double)/max(b_par,1.0e-8_double))
       ! If the last categories are empty, this cutoff corresponds to (0,0)
      elseif( all(catn_orig(icat+1:num_cat_orig) == 0) .and. all(cats_orig(icat+1:num_cat_orig) == 0)  ) then
              cutoff =  max( 6.0_double, (a_par + 6.0_double)/max(b_par,1.0e-8_double))
      else
             ! If both the positive cases and the negative cases are different from 0 for this category, then
             ! it can't be collapsed with any other category, so it has the cutoff of the category where it
             ! was collapsed to.
             if(catn_orig(icat) /= 0 .and. cats_orig(icat) /= 0) then
                  cutoff    = cutoffs(new_cat_index(icat) )
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
                       cutoff = cutoffs(new_cat_index(icat) )
                       num_cases = 0 ! This is a new category, so we reset this counter to zero
                  else ! We are sure that this category is collapsed together with some other one
                      if( new_cat_index(icat) == 1) then ! this is the first cutoff and it is not empty. We know that the
                          ! point before this point is (1.0, 1.0) because this is the first category
                          prev_fpf = 1.0_double
                      else ! We need to compute the previous value, the array can't get messed up because
                          call fpf_CVBM(a_par, b_par, cutoffs(new_cat_index(icat)  - 1) , prev_fpf, one_minus)
                      endif
                          ! Fixed on August 2009 by LP at U of C. To prevent the search for a cutoff above vc_cutoff(numcat
                          ! First compute the theoretical operating point for the category to which this one was collapsed
                          ! need to make sure that the collapsed categories are not the last ones because otherwise there
                          ! no cutoff after this one and the next operating point is (0,0)
                          if( new_cat_index(icat) == num_cat) then
                                 fpf = 0.0_double
                          else
                              call fpf_CVBM(a_par, b_par,  cutoffs(new_cat_index(icat) ) , fpf, one_minus)
                          endif
                          num_cases = num_cases + catn_orig(icat) ! num_cases must start here, because it is the first
                          ! Compute the theoretical fpf for this point, by partitioning the probability between
                          ! this cutoff and the next one
                          fpf = prev_fpf - (prev_fpf - fpf) * num_cases / catn(new_cat_index(icat) )
                          call compute_zdev_plus(fpf, cutoff,ierror)
                          cutoff = - cutoff
                  endif
             elseif(catn_orig(icat) == 0) then ! There are only positive cases
                  ! If collapsed and non collapsed have the same number of cases, then this category
                  ! wasn't collapsed to categories with number of cases > 0, and it takes the same cutoff
                  if  (cats_orig(icat) ==  cats(new_cat_index(icat) ) ) then
                       cutoff = cutoffs(new_cat_index(icat) )
                       num_cases = 0 ! This is a new category, so we reset this counter to zero
                  else ! We are sure that this category is collapsed together with some other one
                      if( new_cat_index(icat) == 1) then ! this is the first cutoff and it is not empty. We know that the
                          ! point before this point is 1.0 1.0 because this is the first category
                           prev_tpf = 1.0_double
                      else ! We need to compute the previous value, the array can't get messed up because
                          call tpf_CVBM(a_par, b_par,  cutoffs(new_cat_index(icat)  - 1) , prev_tpf, one_minus)
                      endif
                          ! Fixed on August 2009 by LP at U of C. To prevent the search for a cutoff above vc_cutoff(numcat
                          ! First compute the theoretical operating point for the category to which this one was collapsed
                          ! need to make sure that the collapsed categories are not the last ones because otherwise there
                          ! no cutoff after this one and the next operating point is (0,0)
                          if( new_cat_index(icat) == num_cat) then
                                 tpf = 0.0_double
                          else
                                 call tpf_CVBM(a_par, b_par,  cutoffs(new_cat_index(icat)) , tpf, one_minus)
                          endif
                          num_cases = num_cases + cats_orig(icat) ! num_cases must start here, because it is the first
                          ! Compute the theoretical fpf for this point, by partitioning the probability between
                          ! this cutoff and the next one
                          tpf = prev_tpf - (prev_tpf - tpf) * num_cases / cats(new_cat_index(icat) )
                          call compute_zdev_plus(fpf, cutoff,ierror)
                          cutoff = (a_par - cutoff)/max(b_par,1.0e-8_double)
                  endif
             endif
       endif
       cutoffs_orig(icat) = cutoff

 ! if needed report log data
 if(idebug==1) then
      call fpf_CVBM(a_par, b_par,  cutoffs_orig(icat) , fpf, one_minus)
      call tpf_CVBM(a_par, b_par,  cutoffs_orig(icat) , tpf, one_minus)
      write(msg,"(1x,I3,4x,I3,4x,I3,'/',I3,4x,I3,'/',I3,4x,2(f14.10,2x) )") icat, new_cat_index(icat), &
                catn_orig(icat), catn(new_cat_index(icat)), cats_orig(icat), cats(new_cat_index(icat)), &
                tpf, fpf
      call DisplayState(msg)
 endif



enddo non_collapsed_categories

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
 end subroutine compute_orig_cutoffs
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine compute_covariance_matrix(a_par, b_par, num_param,idebug, hessian,cov,ierror, err_msg)
!---------------------------------------------------------------------------------------------------
! PURPOSE:  COMPUTATION OF VARIANCES by  Computing the inverse of the hessian to obtain the variance covariance
!            matrix
! ALGORITHM: First the program attempts a straigth inversion, but the hessian might become
!            singular, because of collinearity in the parameters.
!            In this case we remove the redundant direction, by diagonalizing the matrix and eliminating the smallest eigenvalue,
!            which is the one associated with the almost flat direction in the likelihood space (so the second derivative
!            is zero and accordingly its eigenvalues is also zero or close to it).
! NOTE:      LP March 2010: Currently the values of the parameters are no used, but I left them there to keep it similar to the
!            one by proproc and because in the future we might find situations where it would be useful to use such values to
!            stabilize the hessian.
use l_algebra, only: pseudoinverse, m1tom2, m2tom1, sinv
use debugging
implicit none

real(kind = double),intent(in)                        :: a_par, b_par ! MLE of the parameters
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
      write(msg,*) "WARNING:: FAILED to invert Full Hessian to obtain covariance matrix"
  elseif(inversion_error > 0) then
      write(msg,*)"WARNING:: unstable inversion of Hessian  SINV error flag:",inversion_error
  endif
  call DisplayState(msg)
endif

! If either inversion was unstable compute and return the pseudoinverse
! remove the last eigenvalue (using a threshold). Sometimes the hessian is kind of flat, but the inverse
! isn't really too bad, and it is better than the alternative of removing one piece.
if(inversion_error /=0 ) then
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
      write(err_msg,*) "Variances computed with Pseudo-inverse of Hessian matrix"
      ierror  = pseudo_inv_variance ! we used a pseudoinverse, we should flag it
   endif
else ! The message is returned by proproc containing info about the computation of the variances
   write(err_msg,*) "Variances computed inverting the negative Hessian"
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
             write(msg,"('    ', 2(1x,i6),1x, g22.15)")i,j, identity(i,j)
             call DisplayState(msg)
         endif
      enddo
    enddo
endif

!---------------------------------------------------------------------------------------------------
end subroutine compute_covariance_matrix
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------
subroutine compute_MLE(mn, ms, num_cat, catn, cats, idebug, a_par, b_par, cutoffs, &
                        log_like, hessian, ierror, err_msg)
!------------------------------------------------------------------------
! PURPOSE: compute the Maximum Likelihood estimate of the curve parameters using the initial estimates of
!          the cutoffs and parameter curves as sent by the calling program. It returns also the hessian matrix
!          to be used to compute variances.
!ALGORITHM: the optimization is done using TOMS611, a Newton-Raphson trust zone optimizer. Sometimes the initial
!          estimates are slightly changed and the calculation is repeated. this has to to with especially unstable
!          regions. It happens very rarely and when it does, it cycles very few times. The algorithm checks
!          whether the optimization produced meaningful numbers.
!          If the optimization wasn't successful, the algorithm returns Log_like = - Huge(1.0_double), i.e., the
!          smallest value possible so that it would be automatically excluded as a candidate solution (worst possible
!          outcome)
! NOTE:   This procedure does not check if the parameters are in a specific numerical space, but only if they are
!         mathematically acceptable, i.e., it checks if they make sense not if they are "good".
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

 real(kind = double),intent(inout)            :: a_par, b_par ! input: initial estimates, output: final
 real(kind = double), dimension(num_cat-1), intent(inout) :: cutoffs  ! in: initial est. out: MLE cutoffs
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
 real(kind=double):: a_par_in, b_par_in
 integer:: iter ! how many times the optimizer started on a specific a, c optmization
 integer:: max_iter ! Maximum number of iterations that will be attempted

 integer::i, nf ! nf is a flag for the optimization procedure

 real(kind=double) :: fpf, tpf
 character(len = line_length):: msg ! text line to be sent use to log information about the run


 intrinsic random_number


 ! Set the return error to the default (OK)
 ierror = 0

 num_param = num_cat + 1 ! the optimization parameters are a, b, and num_cat - 1 cutoffs

 ! This is a little ad hoc. Basically the idea is that if there are many categories
 ! it becomes really difficult to converge when the the initial values are not
 ! very good. Since in these situation there aren't any multiple maxima, it does not
 ! make much sense to try too hard to converge either, so the non-converged values are dropped
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

 ! Store the initial value of a and b for use in case the first attempts of the MLE fail
 ! and the procedure need to be iterated with different initial estimates (sometimes needed
 ! when curves are very close to the chance line)
 a_par_in = a_par
 b_par_in = b_par

 ! Copy the initial estimates into the parameter vector for the optimizer
 param_vec(1) = a_par
 param_vec(2) = b_par
 param_vec(3:num_param) = cutoffs(1: num_cat - 1) ! cutoffs are in between categories

! Maximum might be on a cups in such a situation, the hessian
! is singular, then we restart the algorithm from a slightly different point.
optimize: do  iter = 1 ,  max_iter
   ! Call the routine that initializes the TOMS611 optimizer.
   call deflt(2,iv,liv,lv,v)
   iv(21) = 0
   iv(17) = 500 ! Maximum number of function evaluations to attempt, only very small curves use more than a few 10's
   iv(18) = 500 ! Maximum number of overall optimization steps " "
   iv(1)= 12 ! tell sumsl that deflt (initialization) was already called
   iv(16) = 1 ! scaling factor scheme to use for the quasi-Newton, for the definition of the trust zone
   d_scale(1:num_param)  = 1.0_double ! scaling factors to define the trust zone, start simple
   ! Call the Trust zone Newton-Raphson (Fisher scoring) to find the Maximum of the Log Likelihood function
   call humsl(num_param, d_scale ,param_vec, like_CVBM, ddlike_CVBM,iv,liv, lv,v, uip,urp, like_CVBM)
   ! Copy the final values from the optimization vector into their variables
   a_par = param_vec(1)
   b_par   = param_vec(2)
   cutoffs = param_vec(3:num_param)
   ! Check whether the optimization was successful, iv(1) is the error flag of the
   ! optimization, details of the optimization run are in the "optimization.info" file
   if( iter < max_iter .AND. & ! If this is the last iteration, we aren't going to try anything more
             (iv(1) == 8 .OR. iv(1) == 9  .OR. iv(1) == 10 .OR. iv(1)== 63 .OR. & ! unstable opt
             .not.check_cutoff_order() ) & ! cutoffs have wrong order
   ) then
       ! send a warning in the log file about a failed attempt in optimization
       call print_convergence(iv(1), ierror,err_msg) ! explicit what kind of convergence was obtained by the optimizer
       if (idebug == 1) then
          call DisplayState(err_msg)
          ! LP April 2007 fixed to prevent crashing while printing
          write(err_msg,"(I2,' TOMS611 call:: a & b =', 2(1x,e16.8))") iter, a_par, b_par
          call DisplayState(err_msg)
       endif
       ! Restart from the original initial values, but change a a little, for even attempts increasing,
       ! for odd attempts decreasing, for a
       param_vec(1)  = a_par_in + (2*modulo(iter,2)-1) * real(iter) /max_iter
       ! b is not changed
       param_vec(2)  = b_par_in
       ! Recompute the cutoffs using the new values of a
       call initialize_cutoffs(param_vec(1),param_vec(2), mn, ms, &
                               num_cat, catn, cats, cutoffs, idebug, ierror, err_msg)
       ! Check whether the new cutoffs initialization succeeded
       if(ierror /= 0) then
          ! LP put a huge value to reject this step (here we are optimzing - log likelihood
          log_like= huge(a_par)
          err_msg = " initialize_cutoffs failed " // err_msg(:line_length - 27)
          if(idebug == 1) call DisplayState(err_msg)
          cycle ! try again unless iteration limit is exceeded
       else
          param_vec(3:num_param) = cutoffs(1:num_cat-1)
       endif
   else
        exit optimize ! either it worked or we exceeded max_iter or the error can't be fixed
   endif
enddo optimize

call print_convergence(iv(1),ierror,err_msg) ! explicit what kind of convergence was obtained
                                             ! by the optimizer

! Note that the procedure returns a - "infinite" log likelihood if something went wrong
if(ierror /= 0) then
          ! MLE failed return the error code and the error message from the optimizer
          ! We might want to do something else in the future, so I left the line
          ierror = 1
          log_like = -huge(1.0_double)
! Check whether the parameter are in bound. It used to happen in the past, it never happened with the
! most recent versions, but better check
elseif(.not.test_parameters(a_par, b_par) ) then
        ierror = 1
        err_msg = "Error:: compute_MLE: parameters out of bounds "
        log_like = -huge(1.0_double)
! Check if the optimizer finished the optimization with a set of cutoffs which are in the wrong order
! This would mean that the category boundaries swapped, which is nonsense from the
! ROC model point of view.
! If 2 points have identical TPF/FPF, clearly the flipping  in the values is nothing to care for, so we
! can switch them  to the right order and keep going. This never happened in the most recent versions
elseif (.not.check_cutoff_order()) then
        if(idebug == 1) then
           call DisplayState("ERROR:: compute_MLE -> OPTIMZER RETURNED CUTOFFS WITH WRONG ORDER")
           write(msg,*) " i     cutoff    fpf          tpf  "
           call DisplayState(msg)
           do i=1,num_cat-1
                  call fpf_CVBM( a_par, b_par, cutoffs(i), fpf)
                  call tpf_CVBM( a_par, b_par, cutoffs(i), tpf)
                  write(err_msg,"(i2,3(1x,e22.15))") i, cutoffs(i), fpf,tpf
                  call DisplayState(msg)
           enddo
        endif
        err_msg = "ERROR:: compute_MLE -> OPTIMZER RETURNED CUTOFFS WITH WRONG ORDER"
        ierror = 1 ! return an error to calling program
        log_like = -huge(1.0_double)
else ! The fit worked out
        nf = 1 ! This is a parameter which is used in optimization, here it does not serve any purpose
                  ! but needs to be set because it is intent(inout) in like_CVBM. Call can be improved.

        ! compute the likelihood value at the maximum, to return to the calling program
        call like_CVBM(num_param, param_vec, nf , log_like, uip)

        log_like = - log_like ! the optimizer is a mimimizer so it works with - log, so we need to correct that

        ! Compute the hessian to estimate the variances
        call ddlike_CVBM(num_param, param_vec, nf, gradient, t_hessian, uip)

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
  subroutine print_convergence(flag, errflg, err_msg)
!------------------------------------------------------------
! Interpretes the flag that describes the quality of the optimization phase just concluded by the program
! sumsl (ot HUMSL - if hessian is used for scoring) and convert it into more understandable sentences
! Should be wrapped together with TOMS611, but I don't want to do it right now.
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
     err_msg = "ERROR:: error in optimization process, please send input file to lpesce@uchicago.edu "
     errflg = 1
   end select

!------------------------------------------------------------
  end subroutine print_convergence
!------------------------------------------------------------


!----------------------------------------------------------------------------------
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
   if ( (cutoffs(i-1) >= cutoffs(i) ) ) then ! Check whether the cutoffs are in the right order
      ! if they are, make sure that it matters, i.e., that the FPF's and TPF's are in fact different
      call fpf_CVBM(a_par, b_par, cutoffs(i-1),fpf1)
      call tpf_CVBM(a_par, b_par, cutoffs(i-1),tpf1)
      call fpf_CVBM(a_par, b_par, cutoffs(i  ),fpf2)
      call tpf_CVBM(a_par, b_par, cutoffs(i  ),tpf2)

      if(  abs(fpf1 - fpf2) < p_q_closeto_zero .and. abs(tpf1 - tpf2) < p_q_closeto_zero  ) then
          swap      =  cutoffs(i)
          cutoffs(i)   =  cutoffs(i-1)
          cutoffs(i-1) = swap
      else
        test_value = .false.
        exit
     endif
   endif
 enddo

 end function check_cutoff_order

!--------------------------------------------------------------------------------------------------
 end subroutine compute_MLE
!--------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
 logical function test_parameters(a_par, b_par)
!-------------------------------------------------------

! DESCRIPTION: CHECK IF PARAMETERS ARE MEANINGFUL and report back to the calling program.
! PURPOSE    : Prevent mathematically absurd or numerically unstable values for the parameters
!              from entering the cycle.
! ALGORITHM  : check the limits of b > 0 (i.e., negative a is accepted because it is a possible if
!              hideous fit.
! Creted march 12th 2010 by LP at U of C
 implicit none

 real(kind=double), intent(in)  :: a_par  ! a as passed by the call
 real(kind=double), intent(in)  :: b_par    ! same but for b

 if( b_par > 0.0_double ) then
  test_parameters = .true.
 else
  test_parameters = .false.
 endif

!----------------------------------------------------------------------------------
 end function test_parameters
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine ddlike_CVBM(num_param, param_vec, nf, grad, h, uip, urp, lk)
!--------------------------------------------------------------------
! PURPOSE: compute the derivatives  (gradient and hessian) of minus log likelihood function for the
!          conventional binormal model

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
call dlike_CVBM(num_param, param_vec, nf, grad, uip)

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

!---------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
end subroutine ddlike_CVBM
!---------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
subroutine dlike_CVBM(num_param,param_vec,nf,grad, uip,urp,lk)
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
! PURPOSE:  Gradient of minus log likelihood function for the conventional binormal model
!           variable and array shapes are done to match the interface with the trust zone Newton-Raphson
!           algorithm TOMS611. Note that Likelihood = Product_i { N!/(k_i!l_i!)  * p_i^k_i * q_i^l_i}
!           Therefore the part of the log-likelihood relevant to the optimizaton is
!           LL = sum_i { log(p_i) * k_i + log(q_i) * l_i }

use problem_data ! data module were problem data constants are stored
use statistic_functions, only: g ! normal distribution and density
use ll_derivatives ! data module with the values and first  derivatives of p_i's and q_i's as a function
                   ! of the parameters

implicit none

integer, intent(IN):: num_param ! number of parameters used in optimization  cutoffs + curve parameters
integer, intent(INOUT):: nf  ! variable used by optimizer, if it is set to zero it means wrong step
integer, dimension(1), intent(INOUT) :: uip ! non used vector

real(kind=double), dimension(num_param), intent(IN):: param_vec ! vector with the parameters for this step
                                                                ! first d_a, then c, then the cutoffs
real(kind=double), dimension(num_param), optional, intent(inout):: urp ! non used
real(kind=double), dimension(num_param), intent(out):: grad ! gradient of the log likelihood function

! internal variables.
real(kind=double):: log_like_val ! Log likelihood value

! variables related to the equations for  the TPF and their derivatives -- FPF is trivial
real(kind=double):: dt1

! model parameters
real(kind=double):: a_par
real(kind=double):: b_par

real(kind=double):: temp

real(kind=double) :: dfpf_dcut ! derivatives of the fpf as a function of the parameters,NOTE: the d in a and b are zero
real(kind=double) :: dtpf_da, dtpf_db, dtpf_dcut ! derivatives of the tpf as a function of the parameters

real(kind=double), parameter:: min_p = 1.0e-10_double ! minimum cateogry probability, to be used only for
                                                      ! the gradient calculation to prevent overflow

integer:: icat ! loop counter over categories

optional lk
external lk

character(len=line_length) ::  msg !  used to log warnings if required to do so


! If the function value  was just computed for this step, we can use the values needed by the gradient
! that were stored by the call to like_CVBM, if not, we need to recompute it to get a new gradient
! - Not a large speed impact.
! Note that this makes sense only if these functions are called by the algorithm TOMS611, otherwise this
! function needs to call like_CVBM all the time to make sure that the p_val and q_val (with the FPF and TPF values
! for the current cutoffs and parameter values) arrays are properly computed
if (uip(1) /=  nf) call like_CVBM(num_param,param_vec,nf,log_like_val,uip)

! Load  the curve parameter values from the optimizer's parameter vector
a_par = param_vec(1)
b_par = param_vec(2)

! To reduce the potential for overflow in the p values, we bound them from below (from above they are bounded by
! 1 by construct. There is no evidence that this might affect negatively the properties of the estimator.
p_val = max(min_p, p_val)
q_val = max(min_p, q_val)

! Build the derivatives of the likelihood function

! Initialize the array with the gradient, and the derivatives of the p_i's and q_i's, that will be used to
! compute the hessian (if needed)
grad    = 0.0_double
dp_dpar = 0.0_double
dq_dpar = 0.0_double

! Compute the summations over the categories for the derivatives
do icat = 1, num_cat - 1

   ! Actually-negative cases term, i.e., - k_i/p_i + k_{i+1}/p_{i+1}
   temp =  - catn(icat) /  p_val(icat)  +  catn(icat+1) / p_val(icat+1)

   ! Derivatives with respect to the  argument of the phi functions. They will be used to compute
   ! the derivatives as a function of the parameters, we do it only for TPF =>
   !    d(p)/d(par) = d(phi)/d(arg) * d(arg) / d(par)

   ! Go over the FPF, only the cutoffs remains for this model, the other ones are zero
   dfpf_dcut  =  - g( -param_vec(2+icat) )

   grad(icat+2) =  temp * dfpf_dcut  ! each cutoff is separated from the other ones

   dp_dpar(icat+1, 2+icat) =  dp_dpar(icat+1,2+icat) + dfpf_dcut
   dp_dpar(icat  , 2+icat) =  dp_dpar(icat  ,2+icat) - dfpf_dcut

   ! Actually-positive cases terms
   temp = - cats(icat)   / q_val(icat)  + cats(icat+1) / q_val(icat+1)

   ! Since they are phi functions, the derivative is the gaussian of their argument
   dt1 = g( a_par - b_par*param_vec(2+icat) )

   dtpf_da   =  dt1
   dtpf_db   =  dt1 * (-param_vec(2+icat))
   dtpf_dcut =  dt1 * (-b_par)

   grad(1)      = grad(1) + temp * dtpf_da
   grad(2)      = grad(2) + temp * dtpf_db
   grad(icat+2) = grad(icat+2)  + temp * dtpf_dcut

   ! Load the derivatives for the P_i's and q_i's, that are needed to compute the expected hessian
   dq_dpar(icat+1,1)      = dq_dpar(icat+1,1)      + dtpf_da
   dq_dpar(icat+1,2)      = dq_dpar(icat+1,2)      + dtpf_db
   dq_dpar(icat+1,2+icat) = dq_dpar(icat+1,2+icat) + dtpf_dcut

   dq_dpar(icat,1)        = dq_dpar(icat,1)        - dtpf_da
   dq_dpar(icat,2)        = dq_dpar(icat,2)        - dtpf_db
   dq_dpar(icat,2+icat)   = dq_dpar(icat,2+icat)   - dtpf_dcut


! If the gradient is meaningless,  return an error message. In this case we will make TOMS611 reject the
! current step
! underflow part (if present in check_number) is necessary or useful at all, work in progress)
   if(.NOT. check_number( grad(icat+2)) .AND. (grad(icat+2) .spne. 0.0_double) ) then
     if(pass_idebug ==1) then
        call DisplayState("ERROR DLIKE_CVBM:: OVERFLOW IN DERIVATIVES OF CUTOFFS ")
        write(msg,*)  icat, param_vec(2+icat), grad(icat+2)
        call DisplayState(msg)
        write(msg,*)   cats(icat), cats(icat+1), catn(icat),catn(icat+1)
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
!------------------------------------------------------------------------------
end subroutine dlike_CVBM
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------



!------------------------------------------------
!------------------------------------------------
subroutine like_CVBM(num_param, param_vec, nf, log_like, uip, urp, unusedvariable)
!------------------------------------------------
! PURPOSE: compute  - LIKELIHOOD FUNCTION  to use to find the cutoffs that
!          maximize the log likelihood function given the curve parameters
!          a & b
!ALGORITHM: adapted from the function from proproc by LP in spring 2010

use statistic_functions, only: phi,g
use problem_data
use ll_derivatives

implicit none

integer, intent(IN):: num_param ! number of parameters (num_cat + 1)
real(kind=double), dimension(num_param), intent(IN):: param_vec ! array with the parameters
     ! in this case (1) is a_par, (2) is b_par, the rest are the cutoffs
integer, intent(INOUT):: nf ! parameter for the optimizer,TOMS611, if set to zero means failure
real(kind=double), intent(out):: log_like ! function value
integer, dimension(1), intent(INout) :: uip ! parameters for the TOMS611, non used
real(kind=double), dimension(num_param), optional,intent(out):: urp
! This declaration has something to do with the Quasi Newton optimizer
! we have to deal with it in a more decent form, once we decided what to
! do with it
optional unusedvariable
external unusedvariable

! curve parameters
real(kind=double):: a_par
real(kind=double):: b_par

character(len=line_length) :: err_msg ! error message from subroutines, Not used by this program
integer:: icat ! loop counter over caregories

! values of FPFs and TPFs and their complements to 1, to be used when computing the log likelihood (to compute
! the p's and the q's
real(kind=double):: fpf, one_minus_fpf, last_fpf, last_one_minus_fpf
real(kind=double):: tpf, one_minus_tpf, last_tpf, last_one_minus_tpf
real(kind=double):: min_prob ! minimum accepted value for a probability

character(len=line_length) ::  msg !  used to log warnings

! load the values of the curve parameters
a_par = param_vec(1)
b_par   = param_vec(2)


min_prob = 10*tiny(1.0_double) ! better a parameter, but f95 does not allow it, minimal value of probability, values
    ! smaller that this one are rounded up to this value
uip(1) = nf ! variable for the TOMS611

! Check the parameter values, if they are in the mathematically meaningful range
! currently it is a very simple check, but is used to be complex, so we left it into a separate function, so that
! it can made complex again if needed.
if( .not.test_parameters(a_par, b_par) ) then
      if(pass_idebug  == 1) then
          write(msg,"('a & b unacceptable:',2(2x,f14.8) )") a_par, b_par    ; call DisplayState(msg)
      endif
      log_like   =  huge( a_par )
      nf = 0
      return
endif



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
       call fpf_CVBM(a_par, b_par, param_vec(2+icat),fpf, one_minus_fpf)
       call tpf_CVBM(a_par, b_par, param_vec(2+icat),tpf, one_minus_tpf)
        ! construct the differences the p(i) = TPF(param_vec(2+icat-1)) - TPF(param_vec(2+icat))
        ! if numerical is lost and so the 2 p's or q's are identical, first
        ! try using the complement, if that does not work, set them to the
        ! numerical basis minimum value. Note that this does not take into
        ! account situations where the 2 cutoffs are so close that their
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
       write(msg,*) "value of a & b ", a_par, b_par    ; call DisplayState(msg)
       write(msg,*) "WARNING:: OVERFLOW OR UNDERFLOW IN LIKE_CVBM ", log_like ; call DisplayState(msg)
       ! write out the values of the cutoffs
       write(msg,*) " CAT CUTOFF_VALUE "        ; call DisplayState(msg)
       do icat=1,num_cat-1
             write(msg,*) icat, param_vec(2+icat) ; call DisplayState(msg)
       enddo
       ! Add the values of TPFs and FPFs
       write(msg,*) " CUTOFF_VALUE FPF's TPF's and their complements to 1 " ; call DisplayState(msg)
       do icat=1,num_cat-1
          call fpf_CVBM(a_par, b_par, param_vec(2+icat), fpf, one_minus_fpf)
          call tpf_CVBM(a_par, b_par, param_vec(2+icat), tpf, one_minus_tpf)
          write(msg,*) icat, fpf, tpf  ; call DisplayState(msg)
          write(msg,*) icat, one_minus_fpf, one_minus_tpf ; call DisplayState(msg)
       enddo
       ! Print pieces of the likelihood function
       write(msg,*) "CAT TRUTH       #CASES     P(OR Q)     CATN(I)*LOG(P) "  ; call DisplayState(msg)
       do icat=1,num_cat
            if(catn(icat) /= 0 ) then
                     write(msg,"(2x,i2,' ACT-NEG ',i3, 2(1x,d14.6))")icat, catn(icat), p_val(icat),&
                                                                catn(icat) * LOG( p_val(icat))
                     call DisplayState(msg)
            endif
            if(cats(icat) /= 0) then
                     write(msg,"( 2x,i2,' ACT-POS ',i3,2(1x,d14.6))")icat,cats(icat), q_val(icat), &
                                                              cats(icat) * LOG( q_val(icat))
                     call DisplayState(msg)
            endif
       enddo
     endif debugging_info
     ! Set the parameters to reject this step
     log_like   = huge(a_par)
     nf = 0
endif


contains
    !----------------------------------------------------------
    !-----------------------------------------------------------
    function check_cutoffs(err_msg) result (test_res)
    ! checks if the cutoffs are in the correct order. It is an internal function, this is why nothing is passed
    implicit none
    logical test_res
    character(len=line_length), intent(out):: err_msg

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
         test_res = .true.
         err_msg = "CUTOFFS ARE INSIDE THE DOMAIN AND IN CORRECT ORDER"
     end if

    end function check_cutoffs

!--------------------------------------------------------------------
end subroutine like_CVBM
!--------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine initialize_cutoffs(a_par_in, b_par_in, num_normal_cases, num_abnormal_cases, num_cat, &
                               catn, cats, cutoffs, idebug,  ierror, err_msg)
!---------------------------------------------------------------------------------------------------
! PURPOSE: compute initial values for the cutoffs in the latent variable space (x) to be returned as
!          the array cutoffs, the estimates do not need to be very accurate, but they need to be
!          stable and prevent overlaps or other situations that would destabilize the following
!          optimization

! ALGORITHM: The subroutines determines the initial estimates by assuming that each cutoff value is partitioning the
!            cases into 2 categories, the cases with classification  smaller than the cutoff  and the ones larger
!            First a lower bound for the smallest cutoff is found (i.e., more or less corresponding to (1,1)),
!            then the algorithm proceeds iteratively using the previous cutoff as lower bound for the current one.
!            This insures also that the cutoffs have the correct ranking and renders relatively easier to force them
!            apart; identical  cutoffs  was one of the main reasons old implementations became unstable.
!            The initial estimate is found as the value that maximizes the log likelihood for the 2 category data.
!            The optimization is done using the Brent method directly on the log likelihood
!            function after the maximum has been bracketed.
! NOTE0:     If cutoffs are almost identical, they are moved apart to prevent ovelapping
!            and consequent divergence of functions and derivative values. It should be reminded that
!            these are initial values, so they don't have to be precise, but only try to avoid selecting a bad initial
!            values
! WARNING:  the last cutoff is meaningless
! HISTORY:  January 2010 - LP, U of C, creation; March 2010 completion

 use debugging, only: DisplayState
 use gen_numerics, only: brent
 use labroc_initial_cutoffs ! needed to exchange data around the canned linear minimization brent
 implicit none

real(kind=double), intent(in):: a_par_in, b_par_in
 integer, intent(in):: num_normal_cases
 integer, intent(in):: num_abnormal_cases
 integer, intent(in):: num_cat ! Number of categories as found by catgrz
 integer, intent(in), dimension(1:num_cat):: catn, cats ! arrays containing categorical data

 real(kind=double), dimension(num_cat-1), intent(out):: cutoffs ! initial estimates
                              ! of cutoffs
 integer, intent(in) :: idebug ! whether to print log information of warnings
 integer, intent(out):: ierror ! Error flag, for failed initialization
                               ! of cutoffs 0 OK, 1 failed
 character(len=line_length), intent(out) :: err_msg


 real(kind=double):: max_err!
                            ! initial value for cutoff, the older algorithm needed a very
                            ! precise estimate (10-10), this isn't true anymore (although
                            ! is hasn't much of a speed impact on the code, since it is done only once
 real(kind=double) :: previous_bound, middle, next_bound ! cutoff values that bound
                                             ! a maximum of the likelihood function
 real(kind=double):: log_like  ! value of the log likelihood at the initial estimate of the
                               ! cutoff, currently not used
 real(kind=double):: inf_bound ! Minimum value of the cutoffs -- boundary for search algorithm.
 integer:: j_cut, j ! loop counter over the cutoffs
 integer:: frst_cat, lst_cat, incr
 character(len=line_length) :: msg


 ierror = 0 ! set to zero the error flag
 ! initialize the values in the module initil_cutoffs, that will be used around brent (minimization)
 a_par = a_par_in
 b_par = b_par_in
 cutoffs = 0.0_double ! Initialize the cutoffs array

 ! the algorithm proceds from up to down, in some extreme situations it might be necessary to do otherwise as
 ! it happens for PROPROC.

 max_err = min(1.0e-3_double, 1.0_double/max(num_normal_cases,num_abnormal_cases)) ! max uncertainty in est.
 frst_cat = 1
 lst_cat  = num_cat - 1
 incr     = 1
 j_cut    = 1
 cat_neg = num_normal_cases
 cat_pos = num_abnormal_cases

if(idebug==1) then
      write (msg, "('INITIAL EST. CUTOFFS for a =', d12.6,'b =', d12.6)")  a_par_in, b_par_in
      call DisplayState(msg)
endif

 ! Find the value after which all ensities are zero (we go 4 stds from the means, as we don't need to be precise and
 ! stability is more important than precision for initial estimates). Set a safety infinite bound for a more extreme value
 inf_bound = min( -5.0_double, a_par/b_par - 5.0_double/b_par)
 previous_bound = min( -4.99_double, a_par/b_par - 4.99_double/b_par)


! the loop is build in this way so that we can go from the last category to the first or from the first to
! the last, depending upon the sign of c. This makes the algorithm symmetrical, removing bias in the sign of c
 CATEGORY_BOUNDARIES: do j =  frst_cat, lst_cat, incr
    cat_neg = cat_neg - catn(j)*incr   ! compute the total number of cases after (i.e., cases that are more positive) than
    cat_pos = cat_pos - cats(j)*incr   ! the cutoff whose initial estimates we are trying to find

    ! bracket current cutoff
    call set_next_bound(a_par, b_par, cat_neg, cat_pos, num_normal_cases, num_abnormal_cases, &
                      previous_bound, middle, next_bound, ierror, err_msg)

                   !write (msg,"('bracketing',i2,3(1x,d14.6))")j, previous_bound, middle, next_bound
                   !call DisplayState(msg)

! Check if there was a bracketing problem with the maximum of the 2 category likelihood function for the first cutoff
    if(ierror /= 0 .AND. j == frst_cat) then ! In this case, we can try to see if we can use a the "backup" estimate of
                                             ! the first cutoff. Normally it works (in fact this in itself is used
                                             ! very rarely and typically for very small areas and noisy data)
          if(idebug==1) then
                   write (msg,*) "bracketing failed!"
                   call DisplayState(msg)
                   call DisplayState(err_msg) ! Display error message from boundary setting routine
          endif
          call set_next_bound(a_par,b_par, cat_neg, cat_pos, &
                         num_normal_cases, num_abnormal_cases, &
                         inf_bound,  middle, next_bound, ierror, err_msg)
    endif

! If the maximum could not be bracketed, we return a failure to the calling program
! with the same error code and error message
    if(ierror /= 0) then
          if(idebug==1) then
            write (msg,*) "Initialize_cuttoffs:: Cutoff", j, "could not be bracketed"
            call DisplayState(msg)
            call DisplayState(err_msg) ! Display error message from boundary setting routine
          endif
          return
    endif
    ! Determine the optimal point using the function brent, which in turn uses function
    ! minus_two_cat_ll because it is a minimization algorithm
    log_like = brent(previous_bound, middle, next_bound, minus_two_cat_ll,max_err,cutoffs(j_cut),ierror)
    if(ierror /= 0 .or.  .not.check_number(log_like) )  return ! Brent failed, need to return - introduced july 1st 2005
    previous_bound = cutoffs(j_cut) + 10*max_err*incr ! Use this optimal value as upper bound for the next cutoff, but add some to
    ! prevent possible ovelapping of cutoffs, which makes the following optimization algorithm very unstable - march 11th 2010
    if(idebug==1) then
            write (msg,*) "Cutoff", j, " = ", cutoffs(j_cut)
            call DisplayState(msg)
    endif
    j_cut = j_cut + incr
 end do CATEGORY_BOUNDARIES


!----------------------------------------------------
!----------------------------------------------------
end subroutine initialize_cutoffs
!----------------------------------------------------
!----------------------------------------------------

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
 subroutine set_next_bound(a_par, b_par, cut_neg, cut_pos, num_neg, num_pos,previous_bound, &
                            middle_bound,bracketing_bound, ierror, err_msg)
! PURPOSE:   Starting from a known value (the previous cutoff -- previous_bound) find another value that brackets
!            the maximum for a two category likelihood roc of a curve determined by a_par and b_par and the categorical data
! ALGORITHM: Use the function value to determine whether the function is decreasing or increasing
!            The algorithm first finds a cutoff whose likelihood is larger than the previous_bound and then looks for
!            a value futher away whose likelihood is smaller than this value, thus bracketing a maximum (function is
!            continuous and differentiable)
! NOTE:      Need to check that the algorithm does not produce oscillations around values when it searches
!            for a lower bound. In particular this can happen when the function is flat around the upper bound,
!            likely because of a numerical issue with the likelihood function or
!            because of the very low probability of some categorical data (normally caused by the fact that the data
!            is fit very poorly by those values of a and b or by a CvBM in general.
!History:    Modified by LP on March 1st 2010 to make it more suitable for CvBM (originally designed for PROPROC)

 implicit none

 real(kind=double), intent(IN) :: a_par, b_par ! ROC curve parameters
 integer, intent(IN) :: cut_neg, cut_pos ! number of act. neg. and act. pos above this cutoff
 integer, intent(IN) :: num_neg, num_pos ! total number of act. neg. and act. pos
 real(kind=double), intent(IN) :: previous_bound

 real(kind=double), intent(out) :: bracketing_bound, middle_bound
 integer, intent(out):: ierror ! Error flag
 character(len=line_length), intent(out):: err_msg

 real(kind=double):: scale ! scale factor to use while trying to bound the maximum
 real(kind=double):: old_lbound ! value used in the iterative scheme
 real(kind=double):: ll, old_ll, previous_ll, middle_ll  ! 2 category log likelihood
  integer:: iter ! iteration counter
 integer, parameter:: max_iter = 500 ! Max number of iterations used to determine the lower
                                 ! bound of this cutoff
 logical:: away ! we are moving away from upper bound or not

 ierror = 0 ! initialize error flag


! Take the slowest scale (either the TPF or FPF, depending upon the value of b
 if(b_par > 1.0_double) then
    scale = +1.0 / b_par
 else
    scale = +1.0
 endif

 !determine the value of the LL at the previous bound, i.e., where we start for the bracketing
 call two_cat_ll(a_par, b_par,cut_neg,cut_pos, num_neg,num_pos, previous_bound, previous_ll, ierror)

! First we seek a point that has higher likelihood than the starting point
 old_lbound = previous_bound + scale
! The algorithm sometimes might need to backtrack, for example if it overshot the maximum by a lot
 away = .true.
 do iter = 1, max_iter
      call two_cat_ll(a_par, b_par, cut_neg, cut_pos, num_neg, num_pos, old_lbound, old_ll, ierror)
      ! for desperate debugging, when initial estimates don't seem to work
      ! write(88,"('m ', i3, (1x,d14.8),1x, d22.15)") iter, old_lbound, old_ll
      if (old_ll > previous_ll) then
          exit ! Theoretically the search could have overshoot the part where the function varies and ended up
               ! in a part of the cutoff space where the function doesn't change anymore and would make the
               ! following part of the search fail because it can't decrease. This never happened in practice so we ignored it
               ! if you determine problems here, you can set a step back and forward a little to
               ! check whether the function changes right around this point and then take the
               ! relevant countermeasures when it does not happen, as of right now it does not seem to be necessary
      elseif(  (old_ll .speq. previous_ll) .and. away) then
                  old_lbound = old_lbound + scale*(1.0_double + log(real(iter))) ! super-linear progression away
                  ! to take care of situations with very large b values.
                  away = .true.
      else ! If it is smaller or non-defined, it is better to go backwards, as we have probably overshot the maximum
        old_lbound = (old_lbound + previous_bound)/2.0_double
        if (previous_bound .speq. old_lbound) exit ! middle point and upper bound set to the same value
        away = .false.
      endif
 enddo

! Check whether a point with larger Log likelihood was found.
 if(iter > max_iter .or. ierror /= 0) then
      write(err_msg,"('set_next_bound:: middle_bound not found, param:',2(1x,f6.3))") a_par, b_par
      ierror = 1
      return
 else
       middle_bound = old_lbound
       middle_ll = old_ll
 endif

! Insure that the correct value has been loaded
 call two_cat_ll(a_par, b_par, cut_neg, cut_pos, num_neg, num_pos, old_lbound, old_ll, ierror)

! Now we search a point with likelihood lower than the middle point, if we do the maximum is bracketed
 do  iter=1, max_iter
       bracketing_bound = old_lbound + scale*(1.0_double + log(real(iter))) ! super-linear progression away
                  ! to take care of situations with very large b values.
       call two_cat_ll(a_par, b_par, cut_neg, cut_pos, num_neg, num_pos, bracketing_bound, ll, ierror)
       if(ierror /= 0) then
             write(err_msg,"('set_next_bound:: D(LL) lost meaning, param:',2(1x,f6.3))") a_par, b_par
             ierror = 1
             return
       endif
       ! Check whether we have found a point with likelihood smaller than the middle point
       if( ll < middle_ll ) then
           if(middle_ll < old_ll) then ! If the last point has larger likelihood than the middle point
              ! update the middle point (it is "closer" to the maximum)
              middle_bound = old_lbound
              middle_ll = old_ll
           endif
           exit
       else! function is still increasing, we try another point
          old_ll = ll
          old_lbound = bracketing_bound
       endif
 enddo

! verify that convergence was obtained
 if(iter > max_iter) then
       write(err_msg,"('set_next_bound:: not found, param:',2(1x,f6.3))") a_par, b_par
       ierror = 1
       return
 endif


 end subroutine set_next_bound

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
function minus_two_cat_ll (cut_off,ierror) result( ll)

 use problem_data
 use labroc_initial_cutoffs
 implicit none

! PURPOSE: generate the negative of the two category log likelihood to be used by optimization
!          procedures
! WARNING: need to change calling scheme, I do not like to change or return arguments for functions
 real(kind=double), intent(in):: cut_off
 integer, intent(out):: ierror
 real(kind = double):: ll

 call two_cat_ll(a_par, b_par,cat_neg,cat_pos, num_normal_cases,num_abnormal_cases, cut_off, ll, ierror)
 ll = -ll

end function minus_two_cat_ll

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
  subroutine two_cat_ll(a_par, b_par, cut_neg, cut_pos, num_neg, num_pos, cutoff, &
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

 real(kind=double),intent(IN):: a_par,b_par
 integer, intent(IN) :: cut_neg, cut_pos ! number of act. neg. and act. pos above this cutoff
 integer, intent(IN) :: num_neg, num_pos ! total number of act. neg. and act. pos
 real(kind=double),intent(IN):: cutoff! cutoff value in the vc space where the log_like with
                                         ! respect to the cutoff has to be computed

 real(kind=double),intent(out):: log_like
 integer,intent(out):: ierror ! error flag

 real(kind=double):: p, one_minus_p  ! probability, for the actually negative (for cat and complement)
 real(kind=double):: q, one_minus_q  ! probability, for the actually positive
 real(kind=double):: smallest_prob

 ! Evaluate the p (and  1-p ) ,  q_j (and  1-q ). The 1-p (or q) terms are needed for numerical
 ! stability

 call fpf_CVBM(a_par, b_par, cutoff, p, one_minus_p)
 call tpf_CVBM(a_par, b_par, cutoff, q, one_minus_q)

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
!-------------------------------------------------------
 end subroutine two_cat_ll
!-------------------------------------------------------
!-------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine initialize_a_b(mn, ms, num_cat, catn, cats, idebug, &
                       a_par_array, b_par_array, ierror, err_msg)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
! PURPOSE: initialize the  values of a and b, the 2 parameters of the conventional binormal model.
!          Moreover checks whether the data received is compatible with the model or if it
!          has analytic solutions, which may be asymptotic and therefore hard to identify numerically
!          (in a, b). In case an analytic solution is possible the model might return the relevant values
!          of the parameters, which should not be taken too seriously.
!ALGORITHM: call least squares functions to determine the value of a and b, while it
!           uses some simple rules to determine whether the datasets are suitable for the MLE
! CHANGES:  Created August 5th 2009
use debugging ! module containing the procedures to log the debugging information on the run

use roc_nonparametric, only: empirical_operating_points_cat
use statistic_functions, only: phi, compute_zdev_plus
implicit none

! Subroutine arguments
integer, intent(in):: mn ! number of actually negative cases
integer, intent(in):: ms ! number of actually positive cases
integer, intent(in) :: num_cat ! Number of categories as created by catgrz
integer, dimension(num_cat), intent(in):: catn, cats ! arrays containing categorical data
integer, intent(in) :: idebug               ! 0 = no debug; 1  = debug


integer, intent(out)  :: ierror ! Error flag
                                !    <- 0  The initialization run smoothly
                                !    <- fit_perfect the data is consistent with a perfect fit, no need to run a MLE
                                !    <- fit_undetermined the data is insufficient to compute an MLE, an "average" solution might be provided
                                !    <- fit_perverse the data is consistent with a perverse fit, no need to run MLE
                                !    <- fit_degenerate the data is consistent with an exact degenerate fit, no need to run an MLE
                                !    <- bad_input the data is no acceptable from the ROC point of view.
                                !    <- fit_fail_init_est   failed fit message because it couldn't produce initial estimates.
character(len = line_length), intent(out):: err_msg    ! description of the error occurred, follows the numeration
          ! for the return values from error_flags, lower routines follows their own error numbering

character(len = line_length) :: msg    ! text for output to out files
real(kind = double), dimension(3), intent(out)  :: a_par_array, b_par_array ! initial estimates of the params
                 ! 3 values to diminish risk of multiple maxima


real(kind = double)  :: a_par, b_par, a_z ! Initial estimates, internal
real(kind=double) :: cumul_fraction(act_neg:act_pos,num_cat+1)! NOTE that the
           ! cumulative fraction that corresponds to category i has index num_cat - i+1
           ! because categories start at + infinity
integer::i ! loop counter
real(kind=double)  :: deviate

a_par_array = -666.0_double
b_par_array = -666.0_double

call empirical_operating_points_cat( num_cat, mn , catn, cumul_fraction(act_neg,1:num_cat+1))
call empirical_operating_points_cat( num_cat, ms , cats, cumul_fraction(act_pos,1:num_cat+1))

if(idebug == 1) then
      write (msg,*) "NUMBER OF CATEGORIES ", num_cat
      call DisplayState(msg)
      do i = 1, num_cat + 1
         write (msg,*) cumul_fraction(act_neg,i), cumul_fraction(act_pos,i)
         call DisplayState(msg)
      enddo

      write (msg,*) "CATEGORICAL DATA USED IN THE MLE, total cat = ", num_cat
      call DisplayState(msg)
      do i = 1, num_cat
         write (msg,*) catn(i), cats(i)
         call DisplayState(msg)
      enddo
      write (msg,*) "END OF CATEGORICAL DATA"
      call DisplayState(msg)

endif

!Check whether, after collapsing the empty/useless categories (categories that do not
! affect the MLE estimation) the fit is known analytically. For the proper binormal model
! this happens only if the points are on the boundaries and the FPF = 0 or TPF=1 axes
! (if more than one point is avaible). The CvBM is more "snaky" and more perfect fits
! are possible. Historically these were called degenerate datasets.

call check_degen_CvBM(num_cat, cumul_fraction, a_par, b_par, ierror, err_msg)


! Capure standard exceptions
select case (ierror)
! These are all fits that render looking unecessary/useless to attempt
! to fit the data. In this case the value of the parameters provided by
! the function that checkes degeneracy is the fit we plan to return
case(fit_perfect, fit_undetermined, fit_perverse, fit_degenerate, bad_input)
     a_par_array(1)  = a_par
     b_par_array (1) = b_par
     if (idebug==1) then
       write (msg,*) "Data cannot be fit, error = ", ierror
       call DisplayState(msg)
       write (msg,*) "a , b = ", a_par_array(1), b_par_array(1)
       call DisplayState(msg)
     endif
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
     err_msg = "failed the initial estimates of a and b"
     return
endif

! We need the value of a_z for the skewed curves and the output
call auc_CVBM(a_par,b_par, a_z, ierror)

! Write out details of fit is required to do so
if(idebug == 1) then
   write(msg,*) 'CONVENTIONAL BINORMAL MODEL L-SQUARES PARAM. INITIAL ESTIMATE'
   call DisplayState(msg)
   write (msg,"( 7x,'A',13x,'B',13x,'AZ')")
   call DisplayState(msg)
   write(msg,"( 3(2x,f12.6))") a_par,b_par, a_z
   call DisplayState(msg)
endif

! Make sure that the parameters aren't in a numerically possibly unstable region
! they are initial estimates anyway. Having exceedingly weird initial estimates can
! destabilize the remaining of the initialization algorithm for no particular reason
! since MLE is not necessarily close to the least squares, in particular if the fits
! are strange-looking.
 a_par_array(1) =   a_par
 b_par_array(1)   = min( max(b_par,.01_double), 100.0_double)


! LOOK FOR TWO SKEWED INITIAL ESTIMATES (WILL HAVE HOOKS BECAUSE MODEL ISN'T PROPER)
! Compute the deviate associated with the current value of Az
call compute_zdev_plus(a_z, deviate, ierror)
! We don't expect any errorsr here. Az is certatainly between 0 and 1
! Fits (one b = .4  and one b = 2.4 ), i.e. skewed and hooky.
 do i = 1, 2
     ! Set the values of a and b
     if( i  ==  1) then
         b_par = .5_double
     else
         b_par = 2.0_double
     endif

     a_par = deviate * sqrt(1.0_double + b_par**2)


     if(idebug ==1 ) then
           call DisplayState('ASYMMETRIC FIXED B INITIAL ESTIMATES')
            write (msg,"( 7x,'A',13x,'B',13x,'AZ')")
            call DisplayState(msg)
            call auc_CVBM(a_par,b_par,a_z,ierror)
            write(msg,"( 3(2x,f12.6))") a_par, b_par, a_z
            call DisplayState(msg)
            ierror = fit_OK ! We don't care if the computation of the area had some precision issues
     endif

     ! load the initial estimate -- make sure they are acceptable values. In principle we don't need to
     ! check b_par as it should have been set right above here, but it is better to make sure, just in
     ! case someone modifies the algorithm without enforcing stability.

     a_par_array(i+1) =   max(a_par,.05_double)
     b_par_array(i+1)   = min( max(b_par,.05_double), 5.0_double)

enddo

!---------------------------------------------------------------------------------------------------
end subroutine initialize_a_b
!---------------------------------------------------------------------------------------------------

!------------------------------------------------------------------
!------------------------------------------------------------------
 Subroutine least_squares_a_b(FPF,TPF, num_pts, mn,ms, a,b, b_equal_one, ierror)
!------------------------------------------------------------------
! Use least squares on the TPF and FPF to make an estimate of the
! parameters a, b for a binormal ROC model, using linear least
! squares over the deviates.

! Warning:: could be cleaned (arcaic syntax and arguable numerical decisions
!           but since it seems to perform its job as good as needed we will
!           leave it to the future generations (LP - winter 2005)

use statistic_functions, ONLY : zdev
implicit none

integer, intent(IN):: num_pts
real(kind=double), dimension(num_pts), intent(IN):: FPF ! Include 0 and 1
real(kind=double), dimension(num_pts), intent(IN):: TPF ! Include 0 and 1
integer, intent(IN):: mn,ms ! # normal and abnormal cases
real(kind=double), intent(out):: a,b
integer, intent(IN) :: b_equal_one ! check if the fit requires b = 1, which
                                   ! forces the algorithm to search for the best symmetrical
                                   ! least squares fit, i.e. b = 1
integer, intent(out):: ierror ! error flag, equal to 1 if the fit failed to produce meaningful values

real(kind=double), dimension(num_pts-2):: x, xs
real(kind=double):: sumx,sumy,sumxy,sumx2,xmean,ymean,xk

integer:: i
integer:: reg_min, reg_max ! Number of points used in the regression, which may be reduced if many points are available
                           ! and there are problematic points with (FPF or TPF)= (0 or 1)

! eliminate the (0,0) and (1,1) points, which are at infinity
! (in the deviate world) and contain no information
x(1:num_pts-2) = fpf(2:num_pts-1)
xs(1:num_pts-2) = tpf(2:num_pts-1)

! If many points are available, eliminate points at the boundary, which are liable to produce very "weird"
! fits.
reg_min = 1
reg_max = num_pts-2
over_pts: do i = 1 , num_pts - 2
         if (x(i) .speq. 0.0_double) then
            if( reg_max - reg_min > 5) then
                 reg_min = reg_min+1
                 cycle over_pts
            else
                 x(i) =  0.5_double/real(mn)
            endif
         endif

         if( abs(x(i) - 1.0_double) < 1.0e-05_double) then
            if( reg_max - reg_min > 5) then
                 reg_max = reg_max-1
                 cycle over_pts
            else
                 x(i) = 1.0_double - 0.5_double/real(mn)
            endif
         endif

         x(i)= zdev(x(i))


         if (xs(i) .speq. 0.0_double) then
            if( reg_max - reg_min > 5) then
                 reg_min = reg_min+1
                 cycle over_pts
            else
                 xs(i) =  0.5_double/real(ms)
            endif
         endif

         if( abs(xs(i) - 1.0_double) < 1.0e-05_double) then
            if( reg_max - reg_min > 5) then
                 reg_max = reg_max-1
                 cycle over_pts
            else
                 xs(i) = 1.0_double - 0.5_double/real(ms)
            endif
         endif

         xs(i)= zdev(xs(i))

enddo over_pts


!     ADJUST THE SEQUENCE OF CUTOFFS TO AVOID inverted  VALUES
!     where a higher cutoff, because of roundoff errors, is lower
!     or equal to a lower cutoff. It increases it by 1/100 of the
!     delta between the previous and the following value to prevent
!     stepping over the next one as well. 5-5-05
!         write(88, *) x(reg_min), xs(reg_min)
do i= reg_min+1, reg_max
     if( x(i)  <=  x(i-1)  )  then
         if ( i < num_pts-2 ) then ! if there is a point after this one, better
              ! make sure that we don't overshoot the next point
              x(i)  = x(i-1)  + 0.01_double*abs( x(i+1) - x(i-1) )
         else
              x(i)  = x(i-1)  + 0.01_double*abs( x(i-1) )
         endif
     endif
     if( xs(i)  <=  xs(i-1)  )  then
         if ( i < num_pts-2 ) then ! if there is a point after this one, better
               ! make sure that we don't overshoot the next point
               xs(i)  = xs(i-1)  + 0.01_double*abs( xs(i+1) - xs(i-1) )
         else
               xs(i)  = xs(i-1)  + 0.01_double*abs( xs(i-1) )
         endif
     endif
enddo

!
!     CALCULATE LEAST SQUARES SOLUTIONS
!
sumx  = 0.0_double
sumy  = 0.0_double
sumxy = 0.0_double
sumx2 = 0.0_double

! Number of points used in the regression
xk = real(reg_max-reg_min+1)


do i = reg_min,reg_max
   sumx  = sumx  + x(i)
   sumx2 = sumx2 + x(i)**2
   sumy  = sumy  + xs(i)
   sumxy = sumxy + xs(i)*x(i)
enddo

xmean = sumx / xk
ymean = sumy / xk

! if there is only one point inside the problem is underdetermined
! and we need to fix something to solve it
if ( b_equal_one == 1 .or. (xk .speq. 1.0_double) ) then
  b = 1.0_double
  a = ymean  -  xmean
else
    b = (xk * sumxy  -  sumx * sumy) / (xk * sumx2  -  sumx**2)
    a = ymean  -  b * xmean
endif

if(check_number(a) .and. check_number(b) )  then
  ierror = 0
else
  ierror = 1
endif



!------------------------------------------------------------------
 end subroutine least_squares_a_b
!------------------------------------------------------------------
!------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
subroutine check_degen_CvBM(num_cat, cumul_fraction, a_par, b_par, ierror, err_msg)
! PURPOSE:  determine if the data can be fit exactly by a pathological ROC curve. Here we define pathological
!           in the sense of being made of horizontal and vertical line segments. This routine is specific to
!           CvBM, other models may have different issues. These pathological curves not only seem to defeat
!           the purpose of fitting a smooth curve (as they aren't smooth) but are produced by asymptotic values
!           of the fitting parameters, i.e., a = +- oo, b = 0 or b = +oo. Therefore, they tend to produce numerical
!           issues as their true solution "cannot" be represented in finite precision. While well designed numerical
!           approaches can deal with this, we prefer to provide a reasonable solution (i.e., a solution that is for
!           all pratical purposes is the same as the exact one), while at the same time flagging the pathological
!           nature of the fit. In principle this routine could be skipped.
! ALGORITHM:The search for pathological behavior is based on the number of categories of the data, starting from
!           the largest that can produce degeneracy (5, for CvBM, for PBM would be 2), and then looking for the conditions
!           that produce such a fit, which can be reduced to: can three vertial or horizontal lines --two of one type
!           and one of the other-- fit exactly through the data (plus (0,0) and (1,1))?
!           The corresponding values of b and a will be set to be very large, but not infinite, and very small, but not zero
!           to prevent potential infinities and non-numerical quantities, that can be handled without specific systems
!           but might create problems across platforms.
! NOTE1:    For "normal" datasets some harmless values are returned, but no attempt is made to return "good" ones.
! NOTE2:    the algorithm assumes that the categories have been collapsed, i.e., only the corners of the
!           operating points are left following the LABROC4 algorithm. Therefore, no more than two points can exist on an
!           axis (or any horizontal or vertical line) -- note that for this purpose (0,0) and (1,1) count.
! NOTE3:    decisions whether a point is on axis or not or whether two points are aligned are based on a somewhat arbitrary
!           cutoff. This is done to prevent possible rounding errors in the determination of the operating points. For very
!           large datasets (in the hundreds of millions) this might create problems. Another approach could have been to use
!           the original categories, but we preferred to make a procedure that could work with operating points (to be honest
!           (I don't remember why and it doesn't look like a good decision to me right now...).
use statistic_functions, only: zdev ! To compute normal deviates, necessary to produce a reasonable value of a.
implicit none

integer, intent(in)                      :: num_cat ! number of categories
real(kind=double), intent(in)            :: cumul_fraction(act_neg:act_pos,num_cat+1)! operating points
real(kind=double), intent(out)           :: a_par, b_par ! the a and b parameters for the degenerate data (or harmless values
                                                         ! for normal data
integer, intent(out)                     :: ierror ! Error flag with information about the data. see program for details
character(len = line_length), intent(out):: err_msg    !Description of the error occurred. Note that this can produce errors
                                                       !if the subroutine is called, say, from R as strings aren't handled
                                                       !very intelligently
real(kind=double), parameter :: zero_frac = 1.0e-8_double !The cutoff value under which we consider the difference between
                                                           !operating points  negligible.
real(kind=double), parameter :: as_zero = 1.0e-4_double! smallest value of a parameter we are going to return
real(kind=double), parameter :: as_infinite = 1.0e4_double! largest value of a parameter we are going to return

intrinsic all

! initialize to non-degenerate, "normal" data
ierror = 0
err_msg = ""
! return some harmless values of a and b if the data is found to be normal.
a_par = 1.0_double
b_par = 1.0_double

! Notes for the following:
! - The first operating point, with index 1, is  (0,0) and the last, the num_cat+1, is (1,1). This is true for every set.
! - The data are assumed to be collapsed, or reduced to LABROC4 corners, at the time this was written there was a procedure
!   called LABROC4_Collapser that could perform that transformation (date is 10-21-09)
! - the return variables corresponding to normal data have been already set, therefore when normal data is found, not action is
!   taken

! The main search is based upon the category number as the typology of degeneracy, perfection and imperfection depends heavily upon
! the number of categories present.


select case (num_cat)

case  (6:) !If collapsed data have more than 5 categories, they can't be degenerate because at most two points can be on
!an axis and 3 points cannot be aligned horizontally or vertically (or they would be collapsed into two).

case(5) ! For degeracy to be possible with 5 categories, the first and the last need *both* to be made of only positive
! or only negative cases (i.e., the need to be on *opposite* axes). Moreover, the remaining two points need to be
! aligned as well. Note that they can't be all aligned or they would be collapsed in two 3 categories. No perfect, perverse
! or undeterminate fits are possible.
     !If the last two points are on the TPF = 0 and TPF = 1, to be degenerate, it is necessary that the two points inside
     !the plot need to be aligned vertically -- there has to be two points inside the plot if there are 5 categories.
     if( ( cumul_fraction(act_pos,2)<zero_frac .and. cumul_fraction(act_pos,num_cat)>1.0_double-zero_frac )           &
        .and.                                                                                                         &
         abs(cumul_fraction(act_neg,4)-cumul_fraction(act_neg, 3)) < zero_frac                                        &
     ) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (vertical) and can be fitted exactly by a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_infinite
            a_par = -zdev(cumul_fraction(act_neg,3))*b_par
    !If the last two points are on the FPF = 0 and FPF = 1, to be degenerate, it is necessary that the two points inside
    !the plot need to be aligned horizontally -- there has to be two points inside the plot if there are 5 categories.
     elseif( ( cumul_fraction(act_neg,2)<zero_frac .and. cumul_fraction(act_neg,num_cat)>1.0_double-zero_frac )         &
       .and.                                                                                                            &
             abs(cumul_fraction(act_pos,4)-cumul_fraction(act_pos,3)) < zero_frac                                       &
     ) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (horizontal) and can be fitted exactly by a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_zero
            a_par = zdev(cumul_fraction(act_pos,3))
     else !The dataset is "normal" and therefore the variable values are already set

     endif
case(4) !For degeracy to be possible with 4 categories, either the first and the last are *both* to made of only positive
! or only negative cases or at least one is on the axis and the other two are aligned. No perfect, perverse of undeterminate
! fits are possible.
    !If the last two points are on the TPF = 0 and TPF = 1,then it is degenerate for sure. (The other point has to be
    !in the middle of the ROC plot, but its location does not affect the set degenerate nature.)
     if( cumul_fraction(act_pos,2)<zero_frac .and. cumul_fraction(act_pos,num_cat)>1.0_double-zero_frac ) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (vertical) and can be fitted exactly by a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_infinite
            a_par = -zdev(cumul_fraction(act_neg,3))*b_par
    !If the last two points are on the FPF = 0 and FPF = 1,then it is degenerate for sure. The other point has to be
    !in the middle.
     elseif( cumul_fraction(act_neg,2)<zero_frac .and. cumul_fraction(act_neg,num_cat)>1.0_double - zero_frac ) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (horizontal) and can be fitted exactly by a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_zero
            a_par = zdev(cumul_fraction(act_pos, 3) )
    !If *EITHER*, but not both, of the last two points are on the TPF = 0 and TPF = 1, to be degenerate, it is necessary that
    !the other two points, which are inside, are aligned vertically. We don't need to check for not both to be on an axis
    !because it was already done in the previous steps.
     elseif( cumul_fraction(act_pos,2)<zero_frac .and.                                                                   &
             abs(cumul_fraction(act_neg,num_cat)-cumul_fraction(act_neg,3)) < zero_frac )then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (vertical) and can be fitted exactly by a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_infinite
            a_par = -zdev(cumul_fraction(act_neg, 3))*b_par
     elseif( cumul_fraction(act_pos,num_cat) > 1.0_double - zero_frac .and.                                              &
             abs(cumul_fraction(act_neg,3)-cumul_fraction(act_neg,2)) < zero_frac ) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (vertical) and can be fitted exactly by a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_infinite
            a_par = -zdev(cumul_fraction(act_neg, 3))*b_par
    !If *EITHER*, but not both, of the last two points are on the FPF = 0 and FPF = 1, to be degenerate, it is necessary that
    !the other two points, which are inside, are aligned horizontally. We don't need to check for not both to be on an axis
    !because it was already done.
     elseif( cumul_fraction(act_neg,2)<zero_frac .and.                                                                   &
         abs(cumul_fraction(act_pos,4)-cumul_fraction(act_pos,3)) < zero_frac ) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (horizontal) and can be fitted exactly by a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_zero
            a_par = zdev(cumul_fraction(act_pos, 3) )
     elseif( cumul_fraction(act_neg,num_cat)>1.0_double-zero_frac .and.                                                  &
          abs(cumul_fraction(act_pos,3)-cumul_fraction(act_pos,2)) < zero_frac ) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (horizontal) and can be fitted exactly by a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_zero
            a_par = zdev(cumul_fraction(act_pos, 3) )
     else !The dataset is "normal" and therefore the variable values are already set

     endif
case(3) ! three categories => two operating points. Fits can be normal, degenerate, indeterminate, perfect and perverse.
    !Check whether a perfect fit is possible: first point has FPF = 0 and second point has TPF = 1.
     if( cumul_fraction(act_neg,2)<zero_frac .and. cumul_fraction(act_pos,num_cat)>1.0_double-zero_frac) then
           ierror = fit_perfect
           err_msg = "The data is degenerate and can be fit exactly by a perfect curve (AUC = 1)"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
           b_par = 1.0_double ! As long as this value is not too close to 0 or oo any value will work as long as the value
                              ! of a that represents infinity is chosen accordingly.
           a_par = as_infinite
     !Check whether a perverse fit is possible: first point has TPF = 0 and second point has FPF = 1.
     elseif( cumul_fraction(act_pos,2)<zero_frac .and. cumul_fraction(act_neg,num_cat)>1.0_double-zero_frac) then
           ierror = fit_perverse
           err_msg = "The data is degenerate and can be fit exactly by a perverse curve (AUC = 0)"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
           b_par = 1.0_double
           a_par = -as_infinite
     ! With 3 categories, if the two points are aligned, the fit is degenerate, whether they are on opposing axes or not.
     elseif( abs(cumul_fraction(act_neg,num_cat)-cumul_fraction(act_neg,2)) < zero_frac )then
           ierror = fit_degenerate
           err_msg = "The data is degenerate (vertical) and can be fitted exactly by a degenerate curve"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
           b_par = as_infinite
           a_par = -zdev(cumul_fraction(act_neg,3))*b_par
     elseif( abs(cumul_fraction(act_pos,num_cat)-cumul_fraction(act_pos,2)) < zero_frac )then
           ierror = fit_degenerate
           err_msg = "The data is degenerate (horizontal) and can be fitted with a degenerate curve"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
           b_par = as_zero
           a_par = zdev(cumul_fraction(act_pos, num_cat))
     !If the last two points are on the TPF = 0 and TPF = 1,then it is undetermined, because the points cannot be aligned
     !as we have tested this already in previous step. We return an average value for a
     elseif( cumul_fraction(act_pos, 2)<zero_frac .and. cumul_fraction(act_pos, num_cat)>1.0_double-zero_frac) then
           ierror = fit_undetermined
           err_msg = "The data is insufficient to identify a specific degenerate curve, return average"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
           b_par = as_infinite
           a_par = -zdev(  (cumul_fraction(act_neg,2)+cumul_fraction(act_neg,num_cat))*.5_double  )*b_par
     !If the last two points are on the FPF = 0 and FPF = 1,then it is undetermined, because the points cannot be aligned
     !as we have tested this already in previous steps. We return an average value for a
     elseif( cumul_fraction(act_neg,2)<zero_frac .and. cumul_fraction(act_neg,num_cat)>1.0_double-zero_frac ) then
           ierror = fit_undetermined
           err_msg = "The data is insufficient to identify a specific degenerate curve, return average"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
           b_par = as_zero
           a_par = zdev( (cumul_fraction(act_pos,2)+cumul_fraction(act_pos,num_cat))*.5_double )
     !If *EITHER*, but not both, of the last two points are on the TPF = 0 and TPF = 1, it is degenerate (alignment doesn't matter)
     ! We don't need to check for not both to be on an axis because it was already done in previous steps
     elseif( cumul_fraction(act_pos,2)<zero_frac .or. cumul_fraction(act_pos,num_cat)>1.0_double-zero_frac) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (horizontal) and can be fitted with a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_infinite
            a_par = -zdev(cumul_fraction(act_neg,num_cat))*b_par
     !If *EITHER*, but not both, of the last two points are on the FPF = 0 and FPF = 1, it is degenerate. (alignment doesn't matter)
     ! We don't need to check for not both to be on an axis because it was already done in previous steps
     elseif(cumul_fraction(act_neg,2)<zero_frac .or. cumul_fraction(act_neg, num_cat)>1.0_double-zero_frac) then
            ierror = fit_degenerate
            err_msg = "The data is degenerate (vertical) and can be fitted with a degenerate curve"
            ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
            b_par = as_zero
            a_par = zdev(cumul_fraction(act_pos, 3) )
     else !The dataset is "normal" and therefore the variable values are already set

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
           b_par = 1.0_double ! As long as this value is not too close to 0 or oo any value will work as long as the value
                              ! of a that represents infinity is chosen accordingly.
           a_par = as_infinite
     !Check whether a perverse fit is possible: first point has TPF = 0 and second point has FPF = 1.
     elseif( cumul_fraction(act_pos,2)<zero_frac .and. cumul_fraction(act_neg,num_cat)>1.0_double-zero_frac) then
           ierror = fit_perverse
           err_msg = "The data is degenerate and can be fit exactly by a perverse curve (AUC = 0)"
           ! The values aren't "exact -- that would require infinities -- but they are close enough for pratical purposes
           b_par = 1.0_double
           a_par = -as_infinite
     ! With 2 categories, if the two points are aligned, the fit is degenerate, whether they are on opposing axes or not.
     else
         ierror = fit_undetermined
         err_msg = "Error: Undetermined problem, CvBM fit requires >= 3 categories "
         ! Nonsense values for a and b
         a_par = -666.0_double
         b_par = -666.0_double
     endif
case(:1)
        ! Less than two categories cannot represent a detection problem and therefore the data is not just insufficient, it is wrong.
        ierror = bad_input
         ! Nonsense values for a and b
         a_par = -666.0_double
         b_par = -666.0_double
end select



!-----------------------------------------------------
end subroutine check_degen_CvBM
!-----------------------------------------------------



!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
 subroutine var_partialauc_CVBM(a_in, b_in, fraction_1, fraction_2, fpf_flag, &
                                          var_a, var_b, cov_a_b, var_p_auc, ierror)
!---------------------------------------------------------------------------------
! PURPOSE:   compute the variance of the partial area either vertically (from FPF_1 to FPF_2 or
!            horizontally (from TPF_1 to TPF_2). Follows the equations from
!            Zhou Obuchowski and McClish (page 122) but uses the notation of
!            Pan & Metz, academic radiology, 4,380, to keep consistency with the other
!            functions in this module. The covariances of the parameters a & b are utilized via the
!            first order approximation known as delta method (see Papoulis).
! ALGORITHM: It starts from the equation:
!
!             p_auc = partial_AUC(...) then uses  analytical first derivatives to compute:
!
!             Var(AREA) = (D{AREA}/D{a})^2  VAR(a) + 2  (D{AREA}/D{a})*(D{AREA}/D{b}) * COV(a,b) +
!            (D{AREA}/D{b})^2 * VAR(b) + ignored terms.
!
!             which is a first order expansion as described in Papoulis. The symmetry of the ROC
!             curves is used, so only the vertical AUC algorithm is utilized, see partial AUC
!             procedure for more details.
! NOTE:

use statistic_functions, only: phi , compute_zdev_plus
implicit none
! parameters a, and b
real(kind=double), intent(IN):: a_in, b_in
real(kind=double), intent(in):: fraction_1, fraction_2 ! These are called fractions
       ! because the can be FPF or TPF depending upon which area are we computing
integer, intent(in):: fpf_flag ! If it is 1, it means that the area computed will
       ! be between FPF_1 and FPF_2, otherwise it means that it will be between TPF_1 and
       ! TPF_2
! variance a, variance b, covariance a,b as estimated by some other procedure. The procedure
! labroc4 provides such variances covariances.
 real(kind=double), intent(IN):: var_a
 real(kind=double), intent(IN):: var_b
 real(kind=double), intent(IN):: cov_a_b
 real(kind=double), intent(OUT):: var_p_auc ! estimated variance of the partial AUC.
 integer,           intent(OUT):: ierror  ! = 0, computatio is OK
                                          ! See partialauc_CVBM for remaining errors
                                          ! will be set to bad_input if the variances and covariances have
                                          ! unacceptable values.
! internal variables.
 real(kind=double):: fpf_1, fpf_2 ! local values, only for vertical p-AUC, see above for details.
 real(kind=double):: x_1, x_2     ! see references and function header
 real(kind=double):: a_par, b_par ! curve parameters for internal use
 real(kind=double):: k0, k1, k2 ! local constants for derivative calculations without any meaning worthy
                                ! mentioning at this time
 real(kind=double) :: h_1, h_2  ! local parameter for derivative calculations, see XOM book.
 real(kind=double) :: d_p_auc_d_a ! value of the derivative of the partial area in a
 real(kind=double) :: d_p_auc_d_b   ! value of the derivative of the partial area in b
 real(kind=double) :: p_auc ! partial AUC value

ierror = 0
! Varify that the input parameters were acceptable values using the partial AUC function
call partialauc_CVBM(a_in, b_in, fraction_1, fraction_2, fpf_flag, p_auc, ierror)
if(ierror /= 0) return

! Check whether the values of the variances and covariances are acceptable
if( var_a < 0.0_double .or. var_b < 0.0_double .or. cov_a_b**2 > var_a * var_b ) then
         ierror = bad_input
         return
else
         ierror = 0 ! initialize
end if


! Note that the horizontal (from TPF_1 to TPF_2)  partial area is identical to a partial area
! computed vertically between  FPF_1 = 1 - TPF_2   FPF_2 = 1 - TPF_1, with  a = a/b and  b = 1/b
 if( fpf_flag ==1) then
      a_par = a_in
      b_par = b_in
      fpf_1 = fraction_1
      fpf_2 = fraction_2
 else
      a_par = a_in/b_in
      b_par = 1.0_double/b_in
      fpf_1 = 1.0_double - fraction_2
      fpf_2 = 1.0_double - fraction_1
 endif

! Defined some numerical constants, which are useful to simplify the final equations
! Following naming from  Zhou, Obuchowski McClish
 k0 =  1.0_double + b_par**2
 k1 =  exp( -a_par**2 / (2.0_double* k0) )
 k2 =  2*pi*k0

! compute the deviate and verify that the calculation was successful for both FPF values
call compute_zdev_plus(fpf_1, x_1, ierror)
if(ierror /= 0) then
        ierror = 1
        return
endif
call compute_zdev_plus(fpf_2, x_2, ierror)
if(ierror /= 0) then
        ierror = 1
        return
endif

! Define a few more constants see reference above
h_1 = sqrt(k0) * ( x_1 + a_par*b_par/k0)
h_2 = sqrt(k0) * ( x_2 + a_par*b_par/k0)

! Compute the partial derivatives in a and b. Please notive that while we are following ZOM for the form, we are
! using the same notation as the other functions. This can be computes using some simple calculus as well, by
! doing a derivative under the sign of the integral of the bivariate normal distribution and then doing the
! resulting integration.
d_p_auc_d_a = (k1/sqrt(k2))*( phi(h_2) - phi(h_1) )
d_p_auc_d_b = (k1/k2) * ( exp(-h_1**2/2.0_double) -exp(-h_2**2/2.0_double) ) - &
              (a_par*b_par*k1/(sqrt(k2)*k0))* ( phi(h_2) - phi(h_1) )

! Apply the delta method (first order expansion)
  var_p_auc =  d_p_auc_d_a**2 * var_a   +   d_p_auc_d_b**2 * var_b  + &
                2.0_double * d_p_auc_d_a * d_p_auc_d_b * cov_a_b


!--------------------------------------------------------------------
end subroutine var_partialauc_CVBM
!--------------------------------------------------------------------
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!--------------------------------------------------------------------
subroutine partialauc_CVBM(a_in, b_in, fraction_1, fraction_2, fpf_flag, &
                                   partial_auc, ierror)
! PURPOSE: compute the partial area either vertically (from FPF_1 to FPF_2 or
!          horizontally (from TPF_1 to TPF_2). Follows the equations from
!          Pan & Metz, academic radiology, 4,380. For the conventional binormal model
!ALGORITHM:It uses symmetry to compute the horizontal partial areas (the
!          partial areas defined by TPF boundaries). For this reason they
!          are computed as vertial partial areas too.
! NOTE:    Follow are much as possible the same variable naming as in the
!          paper mentioned in the "PURPOSE" section.
!WARNING:  Both the horizontal and vertical partial AUC equations in the paper mentioned
!          above are incorrect. However, the equations implemented here have been checked very carefully and
!          compared with numerical simulations.
! HISTORY: LP June 16th 2009: allowed the value b = 0; introduced a minimum value for rho for numerical stability and
!          introduced rounding in the partial area value to prevent "strange" results when exact cancellation should
!          produce zeros or near zeros
!          LP July 2009: changed return errors to make them consistent with other functions (removed ierror =4)
 use statistic_functions, only: phi, compute_zdev_plus, bivar_normal_distrib
 implicit none

 real(kind=double), intent(in):: a_in ! Curve parameters on input (before using symmetry to compute values)
 real(kind=double), intent(in):: b_in ! Curve parameters
 real(kind=double), intent(in):: fraction_1, fraction_2 ! These are called fractions
       ! because the can be FPF or TPF depending upon which area are we computing
  integer, intent(in):: fpf_flag ! If it is 1, it means that the area computed will
       ! be between FPF_1 and FPF_2, otherwise it means that it will be between TPF_1 and
       ! TPF_2
  real(kind=double), intent(out):: partial_auc
  integer,           intent(out):: ierror  ! = 0, computatio is OK
                                           ! = 1 failed calculation: normal deviate could not be computed or other numerical error
                                           ! = bad_input, one or more fractions are larger than 1 or
                                           !      smaller than 0 or the values of the parameters
                                           !      (a and b) are unacceptable
                                           ! = 3 , fractions are almost identical or identical


 ! INTERNAL VARIABLES
 real(kind=double), parameter:: tol = 1.0e-13_double ! tolerance in the computation of tpf and fpf
 real(kind = double):: rho
 real(kind = double):: a_par, b_par ! curve parameters for internal use, corresponding to vertical partial AUC
 real(kind = double):: fpf_1, fpf_2 ! The internal calculations are done using only FPF because the
         ! horizontal partial areas are computed using symmetry, thus again as vertical
 real(kind = double):: y ! see paper
 real(kind = double):: x_1, x_2 ! see paper
 real(kind = double):: F1, F2 ! normal bivariate distributions used in the calculation of the partial area
                               ! see paper.

!Verify that input cutoff values for TPF or FPF are meaningful and also check whether the parameter values
!are acceptable.
 if( fraction_1 > 1.0_double .or. fraction_1 < 0.0_double .or. &
     fraction_2 > 1.0_double .or. fraction_2 < 0.0_double .or. &
     a_in < 0.0_double .or. b_in < 0.0_double) then
           ierror = bad_input
           return
 elseif(fraction_1 >= fraction_2) then
           ierror = bad_input
           return
 elseif(fraction_1 .speq. fraction_2) then
           ierror = 3
           return
 else ! initialize error flag
           ierror = 0
 endif


 ! Note that the horizontal (from TPF_1 to TPF_2)  partial area is identical to a partial area
 ! computed horizontally between  FPF_1 = 1 - TPF_2   FPF_2 = 1 - TPF_1, with  a = a/b and  b = 1/b
 if( fpf_flag ==1) then
      b_par = b_in
      a_par = a_in
      fpf_1 = fraction_1
      fpf_2 = fraction_2
 else
      b_par = 1.0_double/max(b_in,1.0e-8_double)
      a_par = min(a_in*b_par, 1.0e4_double) ! LP June 2009 to remove 0/0 uncertainty and divergence problems
                         ! Note the that we do 1/b first. Max a is at 10,000 for numerical reasons.
      fpf_1 = 1.0_double - fraction_2
      fpf_2 = 1.0_double - fraction_1
 endif

! Compute the partial area following the paper by Pan and Metz (see header of subroutine)
! after correcting the mistakes (leave parameters rho and y unchanged, change the sign
! of the x values)
 rho =  - min( max( b_par / sqrt( 1.0_double + b_par**2 ), 1.0e-8_double), .999999999_double )
 y   =  + a_par / sqrt( 1.0_double + b_par**2 )


! compute the deviate and verify that the calculation was successful for both FPF values
! AGAIN notice that the paper uses -x_1, but the paper is wrong.
call compute_zdev_plus(fpf_1, x_1, ierror)
if(ierror /= 0) then
        ierror = 1
        return
endif

call compute_zdev_plus(fpf_2, x_2, ierror)
if(ierror /= 0) then
        ierror = 1
        return
endif

! Compute the bivariate terms, the same as the paper.
! For now we ignore the return error messages from the bivatiate function, they do not seem to
! affect the outcome of the estimation.
  call  bivar_normal_distrib( x_2,   y , rho, F2, ierror)
  call  bivar_normal_distrib( x_1,   y , rho, F1, ierror)

  ierror = 0
  partial_auc = F2  - F1

  partial_auc = real(nint(partial_auc * 1.0e8_double)) / 1.0e8_double


!--------------------------------------------------------------------
end subroutine partialauc_CVBM
!--------------------------------------------------------------------
!--------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine var_fpf_find_tpf_CVBM(a_par, b_par, var_a, var_b, cov_a_b, fpf, var_tpf, ierror)
!---------------------------------------------------------------------------------------------------
! PURPOSE: given a value for the FPF, find the variance of the corresponding TPF, given the parameters a and b
!---------------------------------------------------------------------------------
! Derive the Standard error of the TPF from the variances and
! covariances of the parameters a & b  using the  first order approximation
! (see Papoulis) known as the delta method.
! Is starts from the equation:
!
! tpf = Phi( a + b* Phi-1(FPF))
!
! d_tpf_d_a = g( a + b* Phi-1(FPF))
! d_tpf_d_b = - g( a + b* Phi-1(FPF))* Phi-1(FPF)

! Then uses is first derivatives to compute:
! VAR(TPF) = (D{TPF}/D{a})^2  VAR(a) + 2  (D{TPF}/D{d_a})*(D{TPF}/D{b}) * COV(a,b) +
!            (D{TPF}/D{b})^2 * VAR(b) + o()
!
! Note that the inverse of FPF does not carry any dependence upon the parameters because of the way the conventional
! binormal model is defined (the distribution for the normals is always centered on 0 and has variance of 1).
! CHANGES: 10/25/2010, LP: changed the sign of deviate because the result of the inversion is minus the inverse of
!          the FPF function [FPF = Phi(-x) =>  -x = Inv_Phi(FPF)

 use statistic_functions, only : phi, compute_zdev_plus, g
 implicit none
 ! parameters a, and b
 real(kind=double), intent(IN):: a_par, b_par

 ! variance a, variance b, covariance a,b
 real(kind=double), intent(IN) :: var_a, var_b, cov_a_b
 real(kind=double), intent(IN) :: fpf ! value of the fpf for the single point available for the fit
 real(kind=double), intent(OUT):: var_tpf
 integer,           intent(out):: ierror  ! = 0, computatio is OK
                                          ! = bad_input, FPF  larger than 1 or
                                          !      smaller than 0 or the values of the parameters
                                          !      (a and b) are unacceptable
                                          ! = 4  normal deviate could not be computed or other numerical error

 real(kind=double):: d_tpf_d_a, d_tpf_d_b
 real(kind=double):: d_phi_d_z

real(kind=double) :: normal_deviate ! inverse of the cumulative normal distribution, Phi

! Check for impossible values of the input parameters
if(fpf > 1.0_double .or. fpf < 0.0_double .or. a_par < 0.0_double .or. b_par <= 0.0_double .or. &
   var_a < 0.0_double .or. var_b < 0.0_double .or. cov_a_b**2 > var_a * var_b ) then
   ierror = bad_input
   return
else
   ierror = 0 ! initialize return error value
end if

! compute the deviate and verify that the calculation was successful
call compute_zdev_plus(fpf, normal_deviate,ierror)
if(ierror /= 0) then
  ierror = 4
  return
endif


! Compute the derivatives of AUC as a function of the parameters a and b
! derivative of Phi(z) with respect of its argument, z. NOTE: the sign is
! the inverse of what usually found on equations because the result of the
! straight inversion, compute_zdev_plus above, is -normal_deviate
 d_phi_d_z =  g(a_par  + b_par*normal_deviate)

 d_tpf_d_a =  + d_phi_d_z
 d_tpf_d_b =  + d_phi_d_z * normal_deviate

! apply the delta method
 var_tpf =  d_tpf_d_a**2 * var_a + d_tpf_d_b**2 * var_b + 2.0_double * d_tpf_d_a * d_tpf_d_b * cov_a_b

!---------------------------------------------------------------------------------
 end subroutine var_fpf_find_tpf_CVBM
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine fpf_find_tpf_CVBM(a_par, b_par, fpf , tpf, ierror)
!---------------------------------------------------------------------------------------------------
! PURPOSE: given a value for the FPF, find the corresponding TPF, given the parameters a and b

use statistic_functions, ONLY :phi, compute_zdev_plus
implicit none

real(kind=double), intent(IN):: a_par, b_par ! parameters of the current fit, used as input
real(kind=double), intent(IN):: fpf ! value of the fpf for the single point available for the fit
real(kind=double), intent(OUT):: tpf ! value of TPF for that fpf
 integer,           intent(out):: ierror  ! = 0, computatio is OK
                                          ! = bad_input, FPF  larger than 1 or
                                          !      smaller than 0 or the values of the parameters
                                          !      (a and b) are unacceptable
                                          ! = 4  normal deviate could not be computed or other numerical error


real(kind=double) :: normal_deviate ! inverse of the cumulative normal distribution, Phi

! Check for impossible values of the input parameters
if(fpf > 1.0_double .or. fpf < 0.0_double .or. a_par < 0.0_double .or. b_par <= 0.0_double) then
   ierror = bad_input
   return
else
   ierror = 0 ! initialize
end if

! compute the deviate and verify that the calculation was successful
call compute_zdev_plus(fpf, normal_deviate, ierror)
if(ierror /= 0) then
  ierror = 4
  return
endif

tpf =  phi(a_par + b_par* normal_deviate)


!---------------------------------------------------------------------------------------------------
end subroutine fpf_find_tpf_CVBM
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine var_tpf_find_fpf_CVBM(a_par, b_par, var_a, var_b, cov_a_b, tpf, var_fpf, ierror)
!---------------------------------------------------------------------------------------------------
! PURPOSE: given a value for the TPF, find the variance of the corresponding FPF, given the parameters a and b
!---------------------------------------------------------------------------------
! Derive the Standard error of FPF from the variances and
! covariances of the parameters a & b  using the  first order approximation
! (see Papoulis) known as the delta method.
! Is starts from the equation:
!
! tpf = Phi( a - b* Phi-1(fpf)) =>
! fpf = Phi( (Phi-1(tpf) - a)/b  ) ! because of change in sign, it is Phi(-z)
!
! d_fpf_d_a = -g( ( Phi-1(tpf) - a)/b ) / b
! d_fpf_d_b = -g( ( Phi-1(tpf) - a)/b ) * ( Phi-1(tpf) - a)/ b^2
!
! Then uses is first derivatives to compute:
! VAR(FPF) = (D{FPF}/D{a})^2  VAR(a) + 2  (D{FPF}/D{a})*(D{FPF}/D{b}) * COV(a,b) +
!            (D{FPF}/D{b})^2 * VAR(b) + o()
! Note that the inverse on TPF does not carry any dependence upon the parameters

 use statistic_functions, only : phi, compute_zdev_plus, g
 implicit none
 ! parameters a, and b
 real(kind=double), intent(IN):: a_par, b_par

 ! variance a, variance b, covariance a,b
 real(kind=double), intent(IN):: var_a,var_b,cov_a_b
 real(kind=double), intent(IN) :: tpf ! value of the fpf for the single point available for the fit
 real(kind=double), intent(OUT):: var_fpf
 integer,           intent(out):: ierror  ! = 0, computatio is OK
                                          ! = bad_input, TPF  larger than 1 or
                                          !      smaller than 0 or the values of the parameters
                                          !      (a and b) are unacceptable
                                          ! = 4  normal deviate could not be computed or other numerical error

 real(kind=double):: d_fpf_d_a, d_fpf_d_b
 real(kind=double):: d_phi_d_z

real(kind=double) :: normal_deviate ! inverse of the cumulative normal distribution, Phi

! Check for impossible values of the input parameters
if(tpf > 1.0_double .or. tpf < .0_double .or. a_par < 0.0_double .or. b_par <= 0.0_double .or. &
   var_a < 0.0_double .or. var_b < 0.0_double .or. cov_a_b**2 > var_a * var_b ) then
         ierror = bad_input
         return
else
         ierror = 0 ! initialize
end if

! compute the deviate and verify that the calculation was successful
call compute_zdev_plus(tpf, normal_deviate,ierror)
if(ierror /= 0) then
  ierror = 4
  return
endif

! Compute the derivatives of AUC as a function of the parameters a and b
! derivative of Phi(z) with respect of its argument, z
 d_phi_d_z =  g( (normal_deviate - a_par) /  b_par )

 d_fpf_d_a = - d_phi_d_z / b_par
 d_fpf_d_b = - d_phi_d_z * ( normal_deviate - a_par) /  b_par**2

! apply the delta method
 var_fpf =  d_fpf_d_a**2 * var_a + d_fpf_d_b**2 * var_b + 2.0_double * d_fpf_d_a * d_fpf_d_b * cov_a_b

!---------------------------------------------------------------------------------
 end subroutine var_tpf_find_fpf_CVBM
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine tpf_find_fpf_CVBM(a_par, b_par, tpf , fpf, ierror)
!---------------------------------------------------------------------------------------------------
! PURPOSE: given a value for the FPF, find the corresponding TPF, given the parameters a and b

use statistic_functions, ONLY : phi, compute_zdev_plus
implicit none

real(kind=double), intent(IN):: a_par, b_par ! parameters of the current fit, used as input
real(kind=double), intent(IN):: tpf ! value of the fpf for the single point available for the fit
real(kind=double), intent(OUT):: fpf ! value of TPF for that fpf
integer,           intent(out):: ierror  ! = 0, computatio is OK
                                          ! = bad_input, TPF  larger than 1 or
                                          !      smaller than 0 or the values of the parameters
                                          !      (a and b) are unacceptable
                                          ! = 4  normal deviate could not be computed or other numerical error


real(kind=double) :: normal_deviate ! inverse of the cumulative normal distribution, Phi

! Check for impossible values of the input parameters
if(tpf > 1.0_double .or. tpf < .0_double .or. a_par < 0.0_double .or. b_par <= 0.0_double) then
   ierror = bad_input
   return
else
   ierror = 0
end if

!Compute the deviate and verify that the calculation was successful
call compute_zdev_plus(tpf, normal_deviate,ierror)
if(ierror /= 0) then
  ierror = 4
  return
end if

!Invert the argument of the deviate of the TPF to obtain the argument of the FPF phi function
normal_deviate = ( -normal_deviate + a_par ) /  b_par

!Compute the fpf value
fpf =  phi(-normal_deviate)


!---------------------------------------------------------------------------------------------------
end subroutine tpf_find_fpf_CVBM
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
 pure subroutine var_auc_CVBM(a_par, b_par, var_a, var_b, cov_a_b, var_auc, ierror)
!---------------------------------------------------------------------------------
! Derive the Standard error of the area (Variance) from the variances and
! covariances of the parameters a & b  using the  first order approximation
! (see Papoulis) known as the delta method.
! Is starts from the equation:
!
! Area = Phi( a/ Sqrt[1+b])
!
! Then uses is first derivatives to compute:
! VAR(AREA) = (D{AREA}/D{a})^2  VAR(a) + 2  (D{AREA}/D{d_a})*(D{AREA}/D{b}) * COV(a,b) +
!            (D{AREA}/D{b})^2 * VAR(b) + o()

 use statistic_functions, only : phi,g

 implicit none
 ! parameters a, and b
 real(kind=double), intent(IN) :: a_par
 real(kind=double), intent(IN) :: b_par
 ! variance a, variance b, covariance a,b
 real(kind=double), intent(IN) :: var_a
 real(kind=double), intent(IN) :: var_b
 real(kind=double), intent(IN) :: cov_a_b
 real(kind=double), intent(OUT):: var_auc
 integer,           intent(OUT):: ierror  ! = 0, computatio is OK
                                          ! = bad_input, TPF  larger than 1 or


 real(kind=double):: d_auc_d_a, d_auc_d_b
 real(kind=double):: d_phi_d_z

! Check for impossible values of the input parameters, including variances
if(a_par < 0.0_double .or. b_par <= 0.0_double .or. &
   var_a < 0.0_double .or. var_b < 0.0_double .or. cov_a_b**2 > var_a * var_b ) then
       ierror = bad_input
       return
else
        ierror = 0
end if



! Compute the derivatives of AUC as a function of the parameters a and b
! derivative of Phi(z) with respect of its argument, z
 d_phi_d_z =  g( a_par / SQRT(1.0_double+b_par**2) )

 d_auc_d_a =  d_phi_d_z / sqrt(1.0_double+b_par**2)
 d_auc_d_b =  d_phi_d_z * (- a_par * b_par)/ (1.0_double+b_par**2)**1.5_double

! apply the delta method
 var_auc =  d_auc_d_a**2 * var_a + d_auc_d_b**2 * var_b + 2.0_double * d_auc_d_a * d_auc_d_b * cov_a_b

!---------------------------------------------------------------------------------
 end subroutine var_auc_CVBM
!---------------------------------------------------------------------------------

!--------------------------------------------------------------------
!--------------------------------------------------------------------
pure subroutine auc_CVBM(a_par,b_par, auc, ierror)
! PURPOSE: Computes the AUC (area under the curve) of the convetional binormal model
!          CvBM
!          using the function phi (the normal cumulative distribution function)
! NOTE:    No reason to set up an error message was found, with the exception of flawed input values
!          on 12/31/2009 LP at UC decided to allow negative values of a because they represent hideous
!          curves, but possible ones, who am I to judge?

 use statistic_functions, only: phi
 implicit none

 real(kind=double), intent(in):: a_par ! the parameter a for the ...
 real(kind=double), intent(in):: b_par ! the parameter b for the...

 integer,           intent(out):: ierror  ! = 0, computatio is OK
 real(kind=double), intent(out):: auc                                           ! = bad_input, TPF  larger than 1 or

! Check for impossible values of the input parameters
if(b_par <= 0.0_double) then
   ierror = bad_input
   return
else
   ierror = 0
   auc =  phi(a_par/sqrt(1.0_double + b_par**2))
end if


!------------------------------------------------------------------------
end subroutine auc_CVBM
!------------------------------------------------------------------------
!------------------------------------------------------------------------

!-------------------------------------------------------
!-------------------------------------------------------
 subroutine fpf_CVBM(a, b, x, fpf, one_minus_fpf)
!-------------------------------------------------------
!-------------------------------------------------------
!PURPOSE: compute the FPF value for a cutoff value (x), and its complement to one,  given the curve parameters (a, b)
!          according to the conventional binormal model

!ALGORITHM: Follow the usual equations, but use the independently computed phi and one_minus_phi to
!           ensure that the obtained values are both numerically meaningful (to avoid overflows and underflows, mostly).
!NOTE:      the function does not attempt to change values when they are too
!           close to zero, this correction should be done by the calling program
!           according to needs
 use statistic_functions, only: norm_dist
 implicit none

 real(kind= double), intent(in):: a ! curve parameter a
 real(kind= double), intent(in):: b ! curve parameter b
 real(kind= double), intent(in):: x ! cutoff in the v space (se M&P paper)


 real(kind= double), intent(out):: fpf ! the value for the false positive fraction
 real(kind= double), optional, intent(out):: one_minus_fpf ! the value for its complement to 1
  ! it is optional because for io purposes and other situations we don't need the complement

  call norm_dist(-x, fpf, one_minus_fpf)


!-------------------------------------------------------
  end subroutine fpf_CVBM
!-------------------------------------------------------
!-------------------------------------------------------

!-------------------------------------------------------
!-------------------------------------------------------
 subroutine tpf_CVBM(a, b, x, tpf, one_minus_tpf)
!-------------------------------------------------------
!-------------------------------------------------------
!PURPOSE: compute the TPF value for a cutoff value (x), and its complement to one,  given the curve parameters (a, b)
!          according to the conventional binormal model

!ALGORITHM: Follow the usual equations, but use the independently computed phi and one_minus_phi to
!           ensure that the obtained values are both numerically meaningful (to avoid overflows and underflows, mostly).
!NOTE:      the function does not attempt to change values when they are too
!           close to zero, this correction should be done by the calling program
!           according to needs

 use statistic_functions, only: norm_dist
 implicit none

 real(kind=double), intent(in):: a ! curve parameter a
 real(kind=double), intent(in):: b ! curve parameter b
 real(kind=double), intent(in):: x ! cutoff in the v space (se M&P paper)


 real(kind=double), intent(out):: tpf ! the value for the false positive fraction
 real(kind=double), optional, intent(out):: one_minus_tpf ! the value for its complement to 1
  ! it is optional because for io purposes and other situations we don't need the complement

real(kind=double) x_arg ! argument of the comulative distribution function

  x_arg = a - b*x
  call norm_dist(x_arg, tpf, one_minus_tpf)


!-------------------------------------------------------
  end subroutine tpf_CVBM
!-------------------------------------------------------
!-------------------------------------------------------

end  module labroc_functions
