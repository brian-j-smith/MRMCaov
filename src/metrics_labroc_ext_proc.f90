 subroutine cvbmroc(num_normal_cases, num_abnormal_cases, NegInput, PosInput,&
                        a_par, b_par,auc, var_AUC)

!
! HISTORY: LP at U of C added the 666.0 values to the returns when error are present.
 use data_types
 use computation_constants, only: ncategory
 use array_dimensions
 use io ! contains parameters and procedures to handle the I/O, warning messaging and error messaging
 use categorization, only: catgrz
 use roc_nonparametric, only: wilcoxon
 use labroc_functions, only: cvbmroc_mle, &
           fit_OK, fit_OK_pseudo, & ! return flag values
           fit_failed, fit_undetermined, fit_perfect, fit_fail_init_est, fit_fail_variances, &
           fit_fail_var_auc_small, fit_degenerate ! return flag values
 use labroc_functions, only: partialauc_cvbm
 use statistic_functions, only: phi, bivar_normal_distrib
 use labroc_out, only: points_on_curve_cvbm,  points_at_cutoffs_cvbm


 implicit none

integer, intent(IN):: num_normal_cases, num_abnormal_cases
! array which contains the input data divided by actually positive
! and actually negative, the ncase constants is used for I/O only and can be
! easily changed.
real (kind=double),dimension(num_normal_cases),intent(IN):: NegInput
real (kind=double),dimension(num_abnormal_cases),intent(IN):: PosInput
real(kind = double),intent(OUT):: a_par
real(kind = double),intent(OUT):: b_par
real(kind=double),intent(OUT) ::  auc ! Area Under the Curve, the area of the fitted ROC curve
real(kind=double),intent(OUT) ::  var_auc ! var of  Area Under the Curve

! # of categories once collapsed
 integer,dimension(act_neg:act_pos,ncategory) :: cat_data
 integer, dimension(act_neg:act_pos, max(num_normal_cases, num_abnormal_cases)):: case_cat ! category by case data
 real(kind=double):: fvalue

 real(kind=double), allocatable, dimension(:) :: cutoffs

 real(kind=double), dimension(:,:), allocatable :: hessian, cov

! VARIABLES TO COMPUTE THE PARTIAL AUC, DEACTIVATED ON AUGUST 2009 BECAUSE THEY AREN'T CURRENTLY USED
! integer:: is_FPF ! logical flag for partial area, either on the FPF (1) or the TPF (0)
! real(kind=double):: lower_FP   ! the lower boundary from where to compute the partial area
! real(kind=double):: upper_FP   ! the upper boundary from where to compute the partial area
! real(kind=double) ::  partial_auc ! partial Area Under the Curve, the area of the fitted ROC curve


 INTEGER, parameter:: positiveislarge = 1 ! are larger values associated with
                           ! higher likelihood of positivity

 integer:: ierror ! error flag for computation errrors
 character(len=line_length) :: err_msg! Error message from called functions, text message

 integer :: num_cat ! Number of categories as found by catgrz

 integer,parameter :: idebug = 0 ! 1 => write log files, 0 or other => don't


!Construct ordinal categories from the real data calling catgrz
 call catgrz(positiveislarge, num_normal_cases,num_abnormal_cases,idebug, cat_data, &
              NegInput, PosInput, num_cat, ncategory, case_cat)

!Allocate the arrays necessary to compute the variance/covariance matrix
 allocate(cutoffs(num_cat-1))
 allocate(hessian(num_cat+1,num_cat+1))
 allocate(cov(num_cat+1,num_cat+1))


 call cvbmroc_mle(num_normal_cases, num_abnormal_cases, num_cat, &
                  cat_data(act_neg,1:num_cat), cat_data(act_pos,1:num_cat), idebug, &
                  a_par, b_par, auc, var_auc, cutoffs, fvalue, ierror, cov, hessian)


!Check the outcome of the MLE optimization
 select case (ierror)
 case(fit_OK , fit_OK_pseudo) ! Successful and variances computed either with normal inverse or with pseudoinverse
 case(fit_failed)! failed fit, deliver the fatal messages
 case (fit_undetermined) ! Not enough data for the fit
 case (fit_degenerate) ! Not enough data for the fit
 case (fit_perfect ) ! perfect fit, not really an error
 !I have decided to consider it a warning because it mostly points to a problem in the dataset
 !collection and I want it to be flagged clearly to the not necessarily careful user - LP UC summer 2005
 !call print_warning_line(err_msg)
       a_par = huge(0.0_double)
       b_par = 1.0_double
       auc   = 1.0_double
       var_auc = 0.0_double ! the MLE estimate of variance in this case is 0, which does not make a lot of
                            ! sense of course
 case (fit_fail_init_est) ! failed initial estimates
       a_par = 666.0_double
       b_par = 666.0_double
       auc   = 666.0_double
       var_auc = 666.0_double
 case (fit_fail_variances) ! failed variances
       var_auc = 666.0_double
 case (fit_fail_var_auc_small) ! failed variances
    var_auc =    auc * ( 1.0_double - AUC) * &
                         .25_double *  ( 1.0_double / num_normal_cases +  1.0_double / num_abnormal_cases )
 case default ! unknown error
 end select

!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
 end subroutine  cvbmroc
!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------


subroutine cvbmrocpartial(a_par, b_par, min, max, flag, auc, ierror)

  use data_types
  use labroc_functions, only: partialauc_CVBM

  real(kind=double), intent(in):: a_par, b_par
  real(kind=double), intent(in):: min, max
  integer, intent(in):: flag
  real(kind=double), intent(OUT):: auc
  integer, intent(OUT):: ierror

  call partialauc_CVBM(a_par, b_par, min, max, flag, auc, ierror)

end subroutine cvbmrocpartial


subroutine cvbmrocfpf2tpf(a_par, b_par, fpf, tpf, ierror)

  use data_types
  use labroc_functions, only: fpf_find_tpf_CVBM

  real(kind=double), intent(in):: a_par, b_par
  real(kind=double), intent(in):: fpf
  real(kind=double), intent(OUT):: tpf
  integer, intent(OUT):: ierror

  call fpf_find_tpf_CVBM(a_par, b_par, fpf, tpf, ierror)

end subroutine cvbmrocfpf2tpf


subroutine cvbmroctpf2fpf(a_par, b_par, tpf, fpf, ierror)

  use data_types
  use labroc_functions, only: tpf_find_fpf_CVBM

  real(kind=double), intent(in):: a_par, b_par
  real(kind=double), intent(in):: tpf
  real(kind=double), intent(OUT):: fpf
  integer, intent(OUT):: ierror

  call tpf_find_fpf_CVBM(a_par, b_par, tpf, fpf, ierror)

end subroutine cvbmroctpf2fpf
