! This module contains the values of the error flags used by the programs
! CREATED by  Lorenzo L Pesce, The University of Chicago August 2009 CE

module error_flags

use data_types, only: double

implicit none

private 

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
public fit_degenerate     ! a perfect fit is possible and it is a snaky combination of straight segments.
public fit_perverse       ! a perfectly perverse fit is possible with a combination of two straight segments. AUC  = 0

integer, parameter :: bad_input          = -1! the categorical data send to the proproc subroutine is bad ROC data
integer, parameter :: fit_OK             = 0 ! fit was successful
integer, parameter :: fit_failed         = 1 ! fit failed in MLE estimation
integer, parameter :: fit_undetermined   = 2 ! Note enough data: there is only one point, and it is inside the ROC space
integer, parameter :: fit_perfect        = 3 ! positives and negatives are perfectly separated, it is more of a warning
integer, parameter :: fit_fail_init_est  = 4 ! initial estimates did not converge
integer, parameter :: fit_fail_variances = 5 ! estimates of variances did not converge
integer, parameter :: fit_OK_pseudo      = 6 ! fit was successful, variances are pseudovariances
integer, parameter :: fit_fail_var_auc_small = 7 ! estimates of var are bad if auc is too small because the 
           ! log likelihood becomes cusp like (at least it looks like). Not worthy worrying for bad datasets
           ! of bad method
integer, parameter :: fit_degenerate = 8     ! the data is such that a snaky fit made of straight segments, as produced by
                                             ! some asymptotic values of a and b for the conventional binormal model, is an
                                             ! exact fit to the data, as such it is also the MLE fit as the perfect fit has the
                                             ! highest possible likelihood.
integer, parameter :: fit_perverse = 9     ! the data is such that a fit made of two straight segments with AUC = 0 is possible

end module error_flags
