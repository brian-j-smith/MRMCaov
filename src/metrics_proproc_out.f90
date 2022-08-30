! Written by L Pesce - U Chicago, starting fall (autumn) 2005
! Module that returns specific quantities for a proper ROC curve.
! This quantities are things like:
! 1) Sequence of points on the curve, for plotting purposes
! 2) FPF and TPF for a series of values (normally the cutoffs or category
!    boundaries in the MLE calculation


 module proproc_out

 use data_types, only: double, operator (.speq.), operator (.spne.)
 use computation_constants, only: pi
 use proproc_computation_constants, only: c_almost_zero
 use proproc_functions, only: fpf_PBM, tpf_PBM
 use io
 use error_flags

 implicit none

 private ! default is that objects of this module are not accessible
 public points_on_curve_PBM ! procedure returns a fpf, tpf sequence for
                               ! plotting purposes. Points are evenly spread

 public points_at_cutoffs_PBM ! procedure returns a sequence of fpf and tpf for a
                          ! prespecified (input) sequence of cutoffs (category
                          ! boundaries, sometimes called betas)

 public print_beta_vs_test_value

 contains

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
real (kind=double) function bnd_val(x,n)
! Check whether the values is bounded in modulus between 10^n and 10^-n to prevent Fortran from printing output
! of the form 10-345 numbers that other programs do not like
! ----------------------------------------------------------------------

real(kind=double), intent(IN):: x
integer, intent(IN):: n
real(kind=double):: abs_x

abs_x = max(abs(x), 10.0_double**(-n))
abs_x = min(abs_x, 10.0_double**(+n))

bnd_val = sign(abs_x, x)


! ----------------------------------------------------------------------
end function bnd_val
! ----------------------------------------------------------------------

!----------------------------------------------------------------------
!------------------------------------------------------------------------------
 subroutine print_beta_vs_test_value(num_cat, mn, ms, neg_cases, pos_cases, positiveislarge, catn, cats,  &
            max_case_per_cat_neg, max_case_per_cat_pos, &
            case_cat_neg, case_cat_pos,  d_a_par, c_par, vc_cutoffs, cov, ierror, err_msg)
!------------------------------------------------------------------------------
! PURPOSE: produce a number of values that can help create a relationship between the test result
!          values and the latent variable space.
! NOTE:    The subroutine prints directly because currently there seem to be no purpose in returning the
!          data to the calling program.
! NOTE1:   The algorithm tests for large AUC values because when the AUC is too large the estimate of beta
!          tends to become unstable.
! WARNING: This subroutine produces an output file and uses it. Handle with care

 use l_algebra, only: qsortd ! for sorting purposes
 use proproc_functions, only: auc_PBM ! to check whether the methods are likely to work

 integer, intent(IN):: num_cat !number of categories
 integer, intent(IN):: mn      !number of actually negative cases
 integer, intent(IN):: ms      !number of actually positive cases


 real (kind=double),dimension(mn), intent(IN) :: neg_cases
 real (kind=double),dimension(ms), intent(IN) :: pos_cases
 integer, intent(IN):: positiveislarge ! whether positivity is for more positive or more negative values

 integer,  dimension(num_cat), intent(IN) :: catn, cats ! categorical data: number of negative and positive cases in each category
 integer, intent(IN):: max_case_per_cat_neg ! the number of negative cases in the category with the most negative cases
 integer, intent(IN):: max_case_per_cat_pos ! the number of positive cases in the category with the most postive cases

 integer, dimension(mn), intent(IN)  :: case_cat_neg ! category by case data
 integer, dimension(ms), intent(IN):: case_cat_pos ! category by case data

 real(kind=double),  dimension(num_cat-1), intent(IN) :: vc_cutoffs
 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters

 real(kind=double),  dimension(num_cat+1,num_cat+1), intent(IN) :: cov ! Variance-covariance matrix as estimated by the MLE

 integer, intent(OUT):: ierror !number of categories
 character(len = line_length), intent(out):: err_msg    ! description of the error occurred. If the fit worked, it
                                                       ! contains information about how were the variances computed
 real(kind=double), dimension(num_cat, max_case_per_cat_neg):: neg_val_by_cat, neg_scratch ! Array where the values for the negative
                     !  cases are loaded in their categories
 real(kind=double), dimension(num_cat, max_case_per_cat_pos):: pos_val_by_cat, pos_scratch! Array where the values for the positive
                     !  cases are loaded in their categories
 integer,  dimension(num_cat) :: neg_cat_pointer, pos_cat_pointer ! used to find which case of this category we are working with

 integer, dimension(max_case_per_cat_neg+max_case_per_cat_pos):: rank_within_cat ! Array of the ranks of cases by truth within category,
                ! It is oversized because fortran does allow us to pick the maximum of the two in the declaration

 real(kind=double),  dimension(0: num_cat) :: cutoffs  ! array with the cutoff values, including the extremes

 real(kind=double):: ranking_coeff ! if positiviy is for more negative values, in each category data needs to be ranked from large to small.

 real(kind=double):: auc ! area under the curve, used to assess whether the estimate of the correlation is likely to be
                         ! reliable

 integer:: i, icat

 character(len = line_length):: msg    ! buffer for printing out results

  ! Initialize
  neg_val_by_cat = 0.0_double
  pos_val_by_cat = 0.0_double
  neg_scratch = 0.0_double
  pos_scratch = 0.0_double
  neg_cat_pointer = 0
  pos_cat_pointer = 0

  ! set the coefficient used to rank cases correctly, up if positivity if for more positive values, otherwise down.
  if(positiveislarge==1) then
         ranking_coeff = 1.0_double
  else
         ranking_coeff = -1.0_double
  endif

  !If the area is large, the variable space usually becomes unstable, moreover there are very few points so the estimation
  ! becomes unstable. in this situation we print a warning
  call auc_PBM(d_a_par, c_par, auc, ierror)
  if(auc > .99_double) then
     write(err_msg, *) "# Estimate of beta values is unreliable. AUC is too large"
     call print_score_to_latent_line(err_msg)
  endif


 ! Arrange negative cases by  category and rank them inside the category
 do i = 1, mn
       icat =  case_cat_neg(i) ! points to the category of current case
       neg_cat_pointer(icat)  =  neg_cat_pointer(icat)  + 1 ! increase the number of cases in this cat by one
       neg_scratch(icat, neg_cat_pointer(icat)) = neg_cases(i)
 enddo

! Sort negative cases within the categories
 do icat = 1, num_cat
    call qsortd(ranking_coeff * neg_scratch(icat,1:catn(icat)), rank_within_cat, catn(icat))
    do i = 1, catn(icat)
       neg_val_by_cat(icat, i) =  neg_scratch(icat, rank_within_cat(i))
    enddo
 enddo


 ! Arrange positive cases by  category and rank them inside the category
 do i = 1, ms
       icat =  case_cat_pos(i) ! points to the category of current case
       pos_cat_pointer(icat)  =  pos_cat_pointer(icat)  + 1 ! increase the number of cases in this cat by one
       pos_scratch(icat, pos_cat_pointer(icat)) = pos_cases(i)
 enddo

 ! Sort positive cases within the categories
 do icat = 1, num_cat
    call qsortd(ranking_coeff * pos_scratch(icat,1:cats(icat)), rank_within_cat, cats(icat))
    do i = 1, cats(icat)
       pos_val_by_cat(icat, i) =  pos_scratch(icat, rank_within_cat(i))
    enddo
 enddo

 ! PRINT OUT THE DATA POINTS, with test result value and categories, sorted
 write(msg, *) "# TEST RESULT VALUES DIVIDED BY CATEGORY/TRUTH RUN "
 call print_score_to_latent_line(msg)
 ! Actually negative cases
 write(msg, *) "#  Actually negative cases "
 call print_score_to_latent_line(msg)
 do icat = 1, num_cat
         do i = 1, catn(icat)
              write(msg, *) icat, neg_val_by_cat(icat, i)
              call print_score_to_latent_line(msg)
         enddo
 enddo
 write(msg, *) " "
 call print_score_to_latent_line(msg)
 write(msg, *) "#  Actually positive cases "
 call print_score_to_latent_line(msg)
 do icat = 1, num_cat
         do i = 1, cats(icat)
              write(msg, *) icat, pos_val_by_cat(icat, i)
              call print_score_to_latent_line(msg)
         enddo
 enddo
 write(msg, *) "# END OF CATEGORY BY TEST RESULT DATA"
 call print_score_to_latent_line(msg)


 ! Compute the largest and smallest value of the cutoff that will be used. They are chosen
 ! so that one has FPF and TPF = 1 and the other one FPF and TPF = 0.
 if (ABS(c_par) < c_almost_zero ) then
        call find_v_c_min(d_a_par, c_par, cutoffs(0), ierror)
        if (ierror == 0) call find_v_c_max(d_a_par, c_par, cutoffs(num_cat), ierror)
 elseif ( c_par <= - c_almost_zero ) then
        cutoffs(0) =  d_a_par * SQRT ( 1.0_double + c_par**2) / &
                       ( 4.0_double * sign( max(.00000000001_double, abs(c_par)),c_par))
        call find_v_c_max(d_a_par, c_par, cutoffs(num_cat), ierror)
 else
        call find_v_c_min(d_a_par, c_par, cutoffs(0), ierror)
        cutoffs(num_cat) =  d_a_par * SQRT ( 1.0_double + c_par**2) / &
                            ( 4.0_double * sign( max(.00000000001_double, abs(c_par)),c_par))
 endif

 if(ierror /= 0) then
      err_msg = "print_beta...: could not set vc boundaries"
      return
 endif

 ! Load the cutoffs, only the part we want
 cutoffs(1:num_cat-1) = vc_cutoffs(1:num_cat-1)


 ! Print header for file with test results value to beta relationship
 write(msg, *) "# EST. OF RELATIONSHIP BETWEEN TEST VALUE RESULTS AND LATENT SPACE"
 call print_score_to_latent_line(msg)
 write(msg, *) "Trth  MEAN          MEDIAN         E[v]          E[v^2]      E[Ln Beta]   ln[M-beta]   SE ln[M-beta]&
               &  SE E[ln-beta]"
 call print_score_to_latent_line(msg)

! And the FPF/TPF file
 write(msg, *) "# EST. OF RELATIONSHIP BETWEEN TEST VALUE RESULTS AND FPF and TPF values"
 call print_score_to_FPF_TPF_line(msg)

 ! Print warning message because the estimates, especially the CIs can be seriously affected for large or small values.
 write(msg, *) "WARNING: ESTIMATES FOR VERY LARGE AND VERY SMALL VALUES ARE UNRELIABLE, ESPECIALLY CONFIDENCE INTERVALS"
 call print_score_to_FPF_TPF_line(msg)

 write(msg, *) "Trth  MEAN          MEDIAN       FPF LBOUND    FPF          FPF UBOUND    TPF LBOUND     TPF    &
               &       TPF UBOUND "
 call print_score_to_FPF_TPF_line(msg)

 ! contruct the test result value to latent variable relationship
 do icat = 1, num_cat
        call compute_values_by_cat(icat, num_cat, d_a_par, c_par, cutoffs(icat-1), cutoffs(icat), catn(icat), cats(icat), &
                      neg_val_by_cat(icat, 1:catn(icat) ),  pos_val_by_cat(icat, 1:cats(icat) ), cov)
 enddo


 !------------------------------------------------------------------------------
  end subroutine print_beta_vs_test_value
 !------------------------------------------------------------------------------


 ! ----------------------------------------------------------------------
 ! ----------------------------------------------------------------------
 subroutine compute_values_by_cat(curr_cat, num_cat, d_a_par, c_par, cat_lower_cutoff, cat_upper_cutoff, neg_in_cat, pos_in_cat, &
                               & neg_val, pos_val, cov)
 ! ----------------------------------------------------------------------
 ! PURPOSE: create an array that connects the different values for the test result value with
 !          the value of beta from the estimate proper binormal distribution.
 !          The method needs to deal with situations where cases of multiple truth are present,
 !          and when labroc5 collapsing might had been implemented
 ! ALGORITHM: A)determine subcategories to labroc4 data, if present.
 !            B) Compute means and medians for the data in the category
 !            C) Compute, analytically, the means and medians for the estimated distributions, conditional
 !               to the values being within the category boundaries.
 ! WARNING: It assumes that test result values arrays are ranked
 ! NOTE: the algorithm can be rewritten more compactly by replacing the FPF and TPF functions with PF functions that make use of
 !       truth. However, I did not really want to to this.Also it could be made faster by doing the calculations as an array instead of
 !       by interval. I see no reason to seek speed and I prefer the interval version because it seems more logical to me.
 use proproc_functions, only: fpf_find_vc, tpf_find_vc
 use l_algebra, only: median ! to compute the median

 implicit none

 integer, intent(IN):: curr_cat ! the index of the current category
 integer, intent(IN):: num_cat !  total number of categories
 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN) :: cat_upper_cutoff, cat_lower_cutoff ! the cutoff estimate, right above and
                       ! below the current
                         ! category
 integer,  intent(IN) :: neg_in_cat, pos_in_cat !  number of negative and positive cases in this category
 real(kind=double), dimension(neg_in_cat), intent(IN) ::  neg_val
 real(kind=double), dimension(pos_in_cat), intent(IN) ::  pos_val
 real(kind=double),  dimension(num_cat+1,num_cat+1), intent(IN) :: cov ! Variance-covariance matrix as estimated by the MLE

 real(kind=double) :: cat_upper_fpf, cat_lower_fpf ! FPF values at the boundaries of the categories, for the
  ! estimated population -- as opposed to the empirical data
 real(kind=double) :: cat_upper_tpf, cat_lower_tpf! TPF values at the boundaries of the categories, for the
  ! estimated population -- as opposed to the empirical data

 real(kind=double) :: lfpf, ltpf ! fpf and tpf at the lower bound of this category
 real(kind=double) :: ufpf, utpf ! fpf and tpf at the upper bound of this category
 real(kind=double) :: delta_p ! difference between fractions (positive or negative) at the boundary of a category
 real(kind=double) :: frac ! fraction of cases between the lower bound and the current value within the category
 real(kind=double) :: old_frac ! fraction of cases between the lower bound and the previous subcategory boundary
 real(kind=double) :: e_val_vc, e_val_vc2 ! expected value of vc and vc^2
 real(kind=double) :: e_val_logbeta ! expected value of log beta
 real(kind=double) :: median_vc ! vc value for median
 real(kind=double) :: rho ! prevalence of positive cases in specific category
 real(kind=double), parameter :: tol = 1.0e-6_double ! tolerace of numerical calculations
 real(kind=double) :: current_lbound ! the current lower bound, for the subcategory
 real(kind=double) :: current_ubound! the current upper bound, for the subcategory
 real(kind=double) :: data_mean, data_median ! the mean and median for the test result values
 real(kind=double), allocatable, dimension(:):: median_array ! array to compute medians

 logical:: first_lbound, last_ubound  ! to take care of situations where the two cutoffs of the current subcategory
                                      ! are the (0,0) or (1,1) cutoffs because in that case their derivatives should not
                                      ! be computed

 integer:: i ! index
 integer:: sub_cat ! Categories that have different values within a single labroc4 category. The labroc4 theorem
   ! guarantees that the value of the maximum likelihood estimate is independent from how many categories are
   ! set within a truth run (a sequence of values of the same truth). Thus they are all put into one, for numerical
   ! reasons. Once the MLE is done, the theorem also provides a way to compute where the boundaries between the
   ! subcategories would have been, namely splitting the probability evenly over the cases.
   ! (i.e., P(case_i) = P_cat/Num_cases_in_cat

 character(len = line_length):: msg    ! buffer for printing out results

 ! Compute fpf and tpf values at the category boundaries.
 call fpf_PBM(d_a_par, c_par, cat_lower_cutoff, lfpf)
 call tpf_PBM(d_a_par, c_par, cat_lower_cutoff, ltpf)

 call fpf_PBM(d_a_par, c_par, cat_upper_cutoff, ufpf)
 call tpf_PBM(d_a_par, c_par, cat_upper_cutoff, utpf)

! if the cutoffs are the one related to (0,0) or (1,1) and they behave in unusual ways because they can be
! functions of the other parameters of being points at infinity
 if(curr_cat == 1) then
     first_lbound = .true.
 else
     first_lbound = .false.
 endif
 last_ubound = .false.

 ! start looking with categories where the points have the same truth
 ! these can't be produced by a labroc5 collapsing because cases of the
 ! same truth belong to the same category already, all is needed is to
 ! split them if they have different values
  current_lbound = cat_lower_cutoff
 if( pos_in_cat == 0) then ! This is a labroc4 category with all actually negative cases
     rho = 0.0_double ! there are only negative cases in this category
     old_frac = 0.0_double ! the previous subcategory of this category was an empty category because this is the first case
     call fpf_PBM(d_a_par, c_par, cat_upper_cutoff, cat_upper_fpf)
     call fpf_PBM(d_a_par, c_par, cat_lower_cutoff, cat_lower_fpf)
     delta_p = cat_lower_fpf - cat_upper_fpf ! The larger the cutoff the smaller the fraction
     sub_cat = 1 ! This is the first subcategory because this is the first "unique" value
     ! Cases have all the same value, so the mean and median are identical and equal to that value
     data_mean = neg_val(1)
     data_median = neg_val(1)
     lfpf = cat_lower_fpf ! The next subgroup of cases will be in a category whose lower bound is the
     ! category bound because it is the first subbroup
     do i = 2, neg_in_cat  ! loop over the negative cases in the category
     ! to determine the possible subcategories
           if ( neg_val(i) .spne. neg_val(i-1) ) then ! Check if a new subcategory was found
               ! when found, compute the values of the previous one (i-1)
               frac = real(i-1)/real(neg_in_cat) ! compute delta FPF from LB of category (max FPF in category)
               ufpf = cat_lower_fpf - frac*delta_p ! fpf associated with that fraction
               ! Find the vc value correspondending to the point just found
               call fpf_find_vc(d_a_par, c_par, max( 1.0e-6_double, ufpf), tol, current_ubound)
               ! prevent numerical approximations to render an inner category cutoff larger that the
               ! category cutoff. Notice that this should happen only for extremely large datasets
               ! because the approximation of the inversion routine is 10^-6. Ideally one could set
               ! the tolerance based upon the number of cases in problem, but I did not want to do
               ! that because it would have rendered the call of the function less focused. One could
               ! have also define a tolerance for the whole program based on the number of cases.
               current_ubound = min(current_ubound, cat_upper_cutoff)
               call compute_values(" - ")
               ! Compute the estimates values and print them
               call fpf_PBM(d_a_par, c_par, current_ubound, ufpf) ! value of fpf for the upper
               ! set new subcategory boundaries by shifting the ones of the old one
               first_lbound = .false. ! we moved away from the first boundary, if we were there.
               old_frac = frac ! the fraction of the current boundary is the previous fraction of the next one
               current_lbound = current_ubound
               lfpf = ufpf
               data_mean = neg_val(i)
               data_median = neg_val(i)
               sub_cat = sub_cat + 1
           endif
     enddo
     ! Take care of the last subcategory of this category, which has the upper bound of the whole category
     if(curr_cat == num_cat) last_ubound = .true. ! Upper cutoff is the largest possible value for the cutoff
     current_ubound = cat_upper_cutoff
     call compute_values(" - ")
     call fpf_PBM(d_a_par, c_par, current_ubound, ufpf) ! value of fpf for the upper
 elseif( neg_in_cat == 0) then! This is a labroc4 category with all actually positive cases
     old_frac = 0.0_double ! the previous subcategory of this category was an empty category because this is the first case
     rho = 1.0_double ! there are only positive cases in this category
     call tpf_PBM(d_a_par, c_par, cat_upper_cutoff, cat_upper_tpf)
     call tpf_PBM(d_a_par, c_par, cat_lower_cutoff, cat_lower_tpf)
     delta_p = cat_lower_tpf - cat_upper_tpf ! The larger the cutoff the smaller the fraction
     sub_cat = 1 ! This is the first subcategory because this is the first "unique" value
     ! Cases have all the same value, so the mean and median are...
     data_mean = pos_val(1)
     data_median = pos_val(1)
     ltpf = cat_lower_tpf ! The next subgroup of cases will be in a category whose lower bound is the
     ! category bound because it is the first subbroup
     do i = 2, pos_in_cat  ! loop over the negative cases in the category
     ! to determine the possible subcategories
           if ( pos_val(i) .spne. pos_val(i-1) ) then ! Check if a new subcategory was found
               ! when found, compute the values of the previous one (i-1)
               frac = real(i-1)/real(pos_in_cat) ! compute delta FPF from LB of category (max FPF in category)
               utpf = cat_lower_tpf - frac*delta_p ! fpf associated
               ! Find the vc value correspondending to the point just found
               call tpf_find_vc(d_a_par, c_par, max( 1.0e-6_double, utpf), tol, current_ubound)
               ! prevent numerical approximations to render an inner category cutoff larger that the
               ! category cutoff. Notice that this should happen only for extremely large datasets
               ! because the approximation of the inversion routine is 10^-6. Ideally one could set
               ! the tolerance based upon the number of cases in problem, but I did not want to do
               ! that because it would have rendered the call of the function less focused. One could
               ! have also define a tolerance for the whole program based on the number of cases.
               current_ubound = min(current_ubound, cat_upper_cutoff)
               call compute_values(" + ")
               call tpf_PBM(d_a_par, c_par, current_ubound, utpf) ! value of fpf for the upper
               ! Initialize new subcategory
               first_lbound = .false. ! we moved away from the first boundary, if we were there.
               old_frac = frac ! the fraction of the current boundary is the previous fraction of the next one
               current_lbound = current_ubound
               ltpf = utpf
               data_mean = pos_val(i)
               data_median = pos_val(i)
               sub_cat = sub_cat + 1
           endif
     enddo
     ! Take care of the last subcategory of this category
     if(curr_cat == num_cat) last_ubound = .true. ! Upper cutoff is maximum possible value for vc
     current_ubound = cat_upper_cutoff
     call compute_values(" + ")
     call tpf_PBM(d_a_par, c_par, current_ubound, utpf) ! value of fpf for the upper
 else ! general case, the category is mixed truth. This category contains both
     ! actually positive and actually negative cases, then it can be either a labroc4
     ! kind of category (all cases have the same test result value) or a labroc5
     ! category (the values are different, but the were collapsed to reduce the number of cases).
     ! LABROC5 categories can only be mixed truth. Since this is an approximation, the LABROC4
     ! theorem does not apply and the category cannot be split into subcategories simply using
     ! the probabilities in the cutoffs
     if(curr_cat == num_cat) last_ubound = .true. ! Upper cutoff is maximum possible value
     old_frac = 0.0_double ! the previous subcategory of this category was an empty category because this is the first case
     frac = 1.0_double ! the current subcategory contains all cases in the category
     data_mean =  real( sum(pos_val) + sum(neg_val) ) / real( neg_in_cat + pos_in_cat)
     ! Compute the data median
     if(allocated(median_array)) deallocate(median_array)
     allocate(median_array(neg_in_cat+pos_in_cat))
     median_array(1:neg_in_cat) = neg_val
     median_array(neg_in_cat+1: neg_in_cat+pos_in_cat) = pos_val
     call median(median_array, neg_in_cat+pos_in_cat, data_median)
     rho =  real(pos_in_cat, kind=double)/(pos_in_cat + neg_in_cat)
     current_ubound = cat_upper_cutoff
     call compute_values(" +-")
 endif

 contains

 subroutine compute_values(cat_type)
    ! PURPOSE: this routine computes and displays values that link the test result space to the different representations in latent space.
    !          It makes use of the change of variables to estimate the standard errors.
    !Note:     if the category is mixed, there is no splitting, so the jacobian is the identity matrix.

     character(len=3), intent(IN):: cat_type ! what kind of category is this one: pos, neg, mix
     real(kind=double), dimension(4,2) :: d_vc_d_th1 ! Derivative of the lower bound (current_lbound) and upper (current_ubound) around the current
                            ! test result value as a function of the MLE parameters, da, c, and the two cutoffs around the current
                            ! truth run (cat_lower_cutoff, cat_upper_cutoff). Variable names are the ones used by LP in July 2007.
     real(kind=double):: se_log_beta_median
     real(kind=double):: se_Ev_log_beta

     real(kind=double), dimension(2) :: PF ! Value of FPF (first) and TPF for the median
     real(kind=double), dimension(2) :: ubound_PF ! Lower bound for 95% CI for FPF (first) and TPF (second)
     real(kind=double), dimension(2) :: lbound_PF ! Upper bound for 95% CI for FPF (first) and TPF (second)

   ! The jacobian has a use only if the categories are made of only negative or only positive cases otherwise the truth runs cannot be split and
   ! thus the Jacobian is the identity matrix. The same is true if the boundaries of the subcategory are the same as the boundaries of the
   ! truth run and it is a waste of time to compute the jacobian -- and can run into numerical issues).
   if(     (  .not.(rho .speq. 0.0_double) .and. .not.(rho .speq. 1.0_double) )  & ! if rho is neither 0 nor 1, it is a mixed category
       .or.  & ! in either case
           (  (cat_lower_cutoff .speq. current_lbound) .and. (cat_upper_cutoff .speq. current_ubound) ) & ! if the bounds are the same, it can't be split
    )then
            d_vc_d_th1(:,1) = (/ 0.0_double, 0.0_double, 1.0_double, 0.0_double /)
            d_vc_d_th1(:,2) = (/ 0.0_double, 0.0_double, 0.0_double, 1.0_double /)
   else ! non-mixed category, with subcategories.
             call compute_th0_th1_jacobian(d_a_par, c_par, cat_lower_cutoff, cat_upper_cutoff, current_lbound, &
                     current_ubound, old_frac, frac, rho, d_vc_d_th1)
   endif


   ! Check whether the boundaries of the subcategories are the end of the category (truth run), in which case the
   ! subcategory is the same as the category one, thus the Jacobian is the identity between the two and zero
   ! otherwise.
   if(cat_lower_cutoff .speq. current_lbound) then
                 d_vc_d_th1(:,1) = (/ 0.0_double, 0.0_double, 1.0_double, 0.0_double /)
   elseif( cat_upper_cutoff .speq. current_ubound ) then
             d_vc_d_th1(:,2) = (/ 0.0_double, 0.0_double, 0.0_double, 1.0_double /)
   endif

   call compute_expected_val_logbeta(d_a_par, c_par, current_lbound, current_ubound, rho , e_val_vc, e_val_vc2, &
                                      e_val_logbeta)
   call compute_median_vc(d_a_par, c_par, current_lbound, current_ubound, rho, median_vc)

   call compute_se_logbeta_median(d_a_par, c_par, current_lbound, current_ubound, median_vc, curr_cat, num_cat,&
                                   cov, rho, d_vc_d_th1,se_log_beta_median)

   call compute_se_EV_logbeta(d_a_par, c_par, current_lbound, current_ubound, first_lbound, last_ubound, e_val_vc, e_val_vc2, &
                                 curr_cat, num_cat, cov, rho, d_vc_d_th1, se_Ev_log_beta)

   Write(msg, "(a3,8(2x,e12.6))") cat_type, bnd_val(data_mean,99), bnd_val(data_median,99),&
                                  bnd_val(e_val_vc,99), bnd_val(e_val_vc2,99), &
                                  bnd_val(e_val_logbeta,99), bnd_val(logbeta_vc(d_a_par, c_par, median_vc),99), &
                                  bnd_val(se_log_beta_median,99), bnd_val(se_Ev_log_beta,99)
   call print_score_to_latent_line(msg)


   call compute_FPF_TPF_median(d_a_par, c_par, current_lbound, current_ubound, median_vc, curr_cat, num_cat,&
                                   cov, rho, d_vc_d_th1, PF, lbound_PF, ubound_PF)


   Write(msg, "(a3,9(2x,e12.6))") cat_type, bnd_val(data_mean,99),bnd_val(data_median,99),&
                                  bnd_val(lbound_PF(1),10),bnd_val(PF(1),10), bnd_val(ubound_PF(1),10), &
                                  bnd_val(lbound_PF(2),10),bnd_val(PF(2),10), bnd_val(ubound_PF(2),10)
   call print_score_to_FPF_TPF_line(msg)


  end subroutine compute_values

! ----------------------------------------------------------------------
end subroutine compute_values_by_cat
! ----------------------------------------------------------------------

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------
 subroutine compute_FPF_TPF_median(d_a_par, c_par, vc1, vc2, median_vc, &
                                 curr_cat, num_cat, cov, rho, d_vc_d_th1, PF, lbound_PF, ubound_PF)
! PURPOSE:  Compute the value of FPF and TPF correspondent to the median, that is the median FPF and TPF.
!           The median is used because it was shown that the estimated done with the expected values and the median are
!           essentially identical and the expected values are lot more complicated to compute, especially for quantities like
!           FPF and TPF. It returns both the estimates values and the lower and upper boundaries for a 95% CI. The
!           CIs are computed using a logistic transformation from the FPF/TPF space to Logit(FPF) Logit (TPF). The SEs are
!           computed for the logits, and so are the CIs. The CIs are then transformed back into the FPF and TPF space.
!           This is a classical trick to compute reliable CIs when the variable of interested is bounded.
! ALGORITHM: The variances are computed using the delta method starting from the Information matrix (inverse of the expected Hessian)
!           of the likelihood function. This is defined in terms of d_a, c, and the cutoffs. The variables are transformed using their
!           first derivatives, so the variances of the logits are computed.
!           It should be notices that the logistic transformation is very non-linear in the tails (for FPF or TPF close to 1 or to 0) and
!           accordingly the inverse is very flat for large values the net effect of this is that for very small or large values the variance
!           estimation, being a linear approximation, fails, sometimes rather miserably. The CIs are within the bounds, but that is because
!           it is a mathematical necessity of the transformations that were applied. The lack of quality of the CIs is compounded by the
!           tendency of the variances and derivatives of FPF and TPF to be also a little unstable in this region (because usually it is the
!           fixed boundary so they don't smoothly go to zero, but rather tend to become unstable and so do the estimated variances of the
!           more extreme cutoffs).
!           The CIs and values can be much better estimated extrapolating from the neighboring values because usually only the last values
!           is unstable. It is fairly easy to find the instability simply by plotting the values.
! WARNING:   Check the ALGORITHM section for instability issues with CIs and how to correct for them.
! PROGRAMMING NOTE: we use PF meaning both FPF and TPF depending from the contenxt (can be both is the class is mixed)
 use proproc_functions, only: density_vc_PBM
 use proproc_functions, only: d_tpf_d_d_a_PBM,d_fpf_d_d_a_PBM ! derivatives of tpf and fpf
 use proproc_functions, only: d_tpf_d_c_PBM,d_fpf_d_c_PBM ! derivatives of tpf and fpf
 use statistic_functions, only: compute_zdev_plus, g, phi

 implicit none

 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: vc1, vc2 ! cutoffs that bound the current test result value
 real(kind=double), intent(IN):: median_vc ! median vc value
 integer, intent(IN):: curr_cat !  current category
 integer, intent(IN):: num_cat !  total number of categories
 real(kind=double),  dimension(num_cat+1,num_cat+1), intent(IN) :: cov ! Variance-covariance matrix as estimated by the MLE
 real(kind=double) :: rho ! fraction of actually positive cases
 real(kind=double), intent(IN), dimension(4,2) :: d_vc_d_th1 ! Derivative of the lower bound (current_lbound) and upper (current_ubound)
                                                              ! around the current
                            ! test result value as a function of the MLE parameters, da, c, and the two cutoffs around the current
                            ! truth run (cat_lower_cutoff, cat_upper_cutoff). Variable names are the ones used by LP in July 2007.
 real(kind=double), dimension(2), intent(OUT) :: PF ! Value of FPF (first) and TPF for the median
 real(kind=double), dimension(2), intent(OUT) :: lbound_PF ! Lower bound for 95% CI for FPF (first) and TPF (second)
 real(kind=double), dimension(2), intent(OUT) :: ubound_PF ! Upper bound for 95% CI for FPF (first) and TPF (second)

! Internal variables
 real(kind=double), dimension(2) :: t ! Value of logit of FPF (first) and TPF for the median
 real(kind=double), dimension(2) :: SE_t ! Value of the SE of the logit of FPF (first) and TPF for the median
 real(kind=double):: d_G_d_m ! deriviative of the implicit function, see under "ALGORITHM"
 real(kind=double):: d_G_d_vc1, d_G_d_vc2 ! deriviative of the implicit function, see under "ALGORITHM"
 real(kind=double):: d_G_d_d_a, d_G_d_c! deriviative of the implicit function, see under "ALGORITHM"

 real(kind=double):: d_m_d_vc1, d_m_d_vc2, d_m_d_d_a, d_m_d_c ! derivatives of m in the parameters of the interval

 real(kind=double), dimension(2) :: d_PF_d_m ! Derivative of FPF and TPF in the median
 real(kind=double):: a, b ! CvBM parameters

 real(kind=double), dimension(2):: d_PF_d_d_a_m_fix, d_PF_d_c_m_fix ! Derivative of FPF and TPF in the median, keeping the median fixed

! real(kind=double), dimension(2) :: d_PF_d_vc1, d_PF_d_vc2, d_PF_d_d_a, d_PF_d_c ! derivatives of log(beta(median)) in the parameters of the interval
 real(kind=double), dimension(2) :: d_t_d_cut1, d_t_d_cut2, d_t_d_d_a, d_t_d_c ! derivatives of log(beta(median)) in the parameters of the MLE estimation
                                           ! that affect this specific estimation.

 !real(kind=double), dimension(2) :: pf_corr ! When converting from FPF/TPF to their logit the net effect on the estimation of the error based on the delta method is
                                            ! dividing the SE by 1/[FPF*(1-FPF)]  (or equivalently for TPF). This term represents this correction.

 integer:: truth_pos ! Whether the current category is made of true positive or true negative cases.

 integer:: i ! Loop counter
 integer:: ierror ! error flag from subroutines


 ierror = 0 ! NEED TO SET PROPER ERROR CODES AND CHECK INPUT.

 ! As a first step we check whether the category is made of only positives or only negatives or mixed.
 ! Mixed categories cannot be split, so the jacobian has to be the indentity matrix and the part we compute is
 ! 1 between the cutoffs and zero otherwise.
 if ( rho .speq. 1.0_double) then
       truth_pos = 1
 elseif ( rho .speq. 0.0_double) then
       truth_pos = 0
 else ! Mixed category
       truth_pos = 2
 endif

! Now we use the implicit function G to compute  the derivatives of the median in the parameters
d_G_d_m    =  - mix_density2(median_vc)
d_G_d_vc1 =  0.5_double * mix_density2(vc1)
d_G_d_vc2 =  0.5_double * mix_density2(vc2)
d_G_d_d_a  =  d_pf_d_d_a2(median_vc) - 0.5_double*d_pf_d_d_a2(vc1) - 0.5_double*d_pf_d_d_a2(vc2)
d_G_d_c    =  d_pf_d_c2(median_vc)   - 0.5_double*d_pf_d_c2(vc1)   - 0.5_double*d_pf_d_c2(vc2)

d_m_d_vc1 = - d_G_d_vc1 / d_G_d_m
d_m_d_vc2 = - d_G_d_vc2 / d_G_d_m
d_m_d_d_a = - d_G_d_d_a / d_G_d_m
d_m_d_c   = - d_G_d_c / d_G_d_m


! Compute a and b
 b = (1.0_double + c_par) / (1.0_double - c_par)
 a = d_a_par * sqrt(1.0_double + b**2) / sqrt(2.0_double)

! Compute the FPF/FPF  function in the median, and compute their logit or deviate or other function
! if this transformation is required.

call fpf_PBM(d_a_par, c_par, median_vc, PF(1))
call tpf_PBM(d_a_par, c_par, median_vc, PF(2))

!If So desired here are the transformations. All these operations if explicitly used should
! be automatized in some form. As of June 2008 the Logit seem worthless for this kind of calc.
! Logit transformation
!t(:) = log( PF(:)/(1.0_double - PF(:)) )
!Deviate transformation
!call compute_zdev_plus(PF(1), t(1), ierror)
!call compute_zdev_plus(PF(2), t(2), ierror)


! NEED TO TRAP ERROR!!!! STILL .... 2001-2008 (Benny)

! Compute the derivative of the FPF/FPF  function in the median

d_PF_d_m(1) = -density_vc_PBM(d_a_par, c_par, median_vc,  0)
d_PF_d_m(2) = -density_vc_PBM(d_a_par, c_par, median_vc,  1)

! Compute the partial derivative of the the  function FPF/TPF in d_a and c

d_PF_d_d_a_m_fix(1) = d_fpf_d_d_a_PBM(d_a_par, c_par, median_vc)
d_PF_d_d_a_m_fix(2) = d_tpf_d_d_a_PBM(d_a_par, c_par, median_vc)

d_PF_d_c_m_fix(1)  = d_fpf_d_c_PBM(d_a_par, c_par, median_vc)
d_PF_d_c_m_fix(2)  = d_tpf_d_c_PBM(d_a_par, c_par, median_vc)


! Construct the derivatives of the logit  PF { FPF/TPF } in the MLE parameters
! First build the derivatives of PF, then transform them to the logit derivatives

 d_t_d_cut1(:) = d_PF_d_m(:) * (  d_m_d_vc1 * d_vc_d_th1(3,1)  + d_m_d_vc2 * d_vc_d_th1(3,2))
 d_t_d_cut2(:) = d_PF_d_m(:) * (  d_m_d_vc1 * d_vc_d_th1(4,1)  + d_m_d_vc2 * d_vc_d_th1(4,2))
 d_t_d_d_a(:) = d_PF_d_d_a_m_fix(:) + d_PF_d_m(:)*(d_m_d_vc1*d_vc_d_th1(1,1) + d_m_d_vc2*d_vc_d_th1(1,2) + d_m_d_d_a)
 d_t_d_c(:) = d_PF_d_c_m_fix(:) + d_PF_d_m(:)*(d_m_d_vc1*d_vc_d_th1(2,1) + d_m_d_vc2*d_vc_d_th1(2,2) + d_m_d_c)

! Apply the delta method to FPF and TPF, the conversion to t will be done when computing the CIs

do i=1, 2 ! Loop over FPF and TPF

   ! If it is the first category, then the cutoff before the current is a fixed point
   ! so it does not contribute to the derivation
   if(curr_cat == 1) then
       se_t(i) = sqrt( &
           d_t_d_d_a(i)**2 *  cov(1,1) + &
           d_t_d_c(i)**2   *  cov(2,2) + &
           d_t_d_cut2(i)**2 *  cov(curr_cat  + 2,curr_cat  + 2) + &
           2.0_double * d_t_d_d_a(i) * d_t_d_c(i)   * cov(1,2)               + &
           2.0_double * d_t_d_d_a(i) * d_t_d_cut2(i) *  cov(1, curr_cat  + 2) + &
           2.0_double * d_t_d_c(i)   * d_t_d_cut2(i) *  cov(2, curr_cat  + 2)   &
         )
   ! If it is the last category, then the cutoff after the current is a fixed point
   ! so it does not contribute to the derivation
   elseif(curr_cat == num_cat) then
       se_t(i) = sqrt( &
           d_t_d_d_a(i)**2 *  cov(1,1) + &
           d_t_d_c(i)**2   *  cov(2,2) + &
           d_t_d_cut1(i)**2 *  cov(curr_cat - 1 + 2, curr_cat - 1 + 2) + &
           2.0_double * d_t_d_d_a(i) * d_t_d_c(i)   * cov(1,2) + &
           2.0_double * d_t_d_d_a(i) * d_t_d_cut1(i) *  cov(1, curr_cat - 1 + 2) + &
           2.0_double * d_t_d_c(i)   * d_t_d_cut1(i) *  cov(2, curr_cat - 1 + 2) &
         )
   else ! Any category in the middle has all the cutoffs
       se_t(i) = sqrt( &
           d_t_d_d_a(i)**2 *  cov(1,1) + &
           d_t_d_c(i)**2   *  cov(2,2) + &
           d_t_d_cut1(i)**2 *  cov(curr_cat - 1 + 2, curr_cat - 1 + 2) + &
           d_t_d_cut2(i)**2 *  cov(curr_cat  + 2,curr_cat  + 2) + &
           2.0_double * d_t_d_d_a(i) * d_t_d_c(i)   *  cov(1,2)                 + &
           2.0_double * d_t_d_d_a(i) * d_t_d_cut1(i) *  cov(1, curr_cat - 1 + 2) + &
           2.0_double * d_t_d_d_a(i) * d_t_d_cut2(i) *  cov(1, curr_cat  + 2)    + &
           2.0_double * d_t_d_c(i)   * d_t_d_cut1(i) *  cov(2, curr_cat - 1 + 2) + &
           2.0_double * d_t_d_c(i)   * d_t_d_cut2(i) *  cov(2, curr_cat  + 2)    + &
           2.0_double * d_t_d_cut1(i) * d_t_d_cut2(i) *  cov(curr_cat - 1 + 2, curr_cat  + 2) &
         )
   endif

enddo


! Apply the transformation from SE(FPF) to SE(logit(PF) and compute the 95% CI based on the assumption that
! the variable is approximately normal.pf_corr is necessary if a stranformation like logit is performed
! Correction for the logit transformation or probit trnasformation, if used
! pf_corr(:) =  PF(:)*(1.0_double - PF(:))
! pf_corr(:) =   g(t(:))

lbound_PF(:) = t(:) - 1.96_double*se_t(:) !/pf_corr(:)
ubound_PF(:) = t(:) + 1.96_double*se_t(:) !/pf_corr(:)

lbound_PF(:) = PF(:) - 1.96_double*se_t(:) !/pf_corr(:)
ubound_PF(:) = PF(:) + 1.96_double*se_t(:) !/pf_corr(:)

!write(*,*) curr_cat, num_cat
!write(*,*) t(1), se_t(1), lbound_PF(1), ubound_PF(1)

! Apply the inverse transformation to obtain the 95% CI for the FPF and TPF
do i=1,2
!  lbound_PF(i) = inv_logit(lbound_PF(i))
!  ubound_PF(i) = inv_logit(ubound_PF(i))
!  lbound_PF(i) = phi(lbound_PF(i))
!  Ubound_pf(i) = phi(ubound_PF(i))
  lbound_PF(i) =  Max(lbound_PF(i), 0.0_double)
  Ubound_pf(i) =  Min(ubound_PF(i),1.0_double)
enddo


! INTERNAL FUNCTIONS
contains
 real(kind=double) function inv_logit(x)
   real(kind=double), intent(IN):: x

  inv_logit = exp(x)/(1.0_double + exp(x))

 end function inv_logit

 real(kind=double) function mix_density2(vc) result(mix_density)
 ! function tht computes the density as a function of vc, depending upon the truth of the category that is being
 ! analyzed
   real(kind=double), intent(IN):: vc

   if(truth_pos == 1) then
      mix_density = density_vc_PBM(d_a_par, c_par, vc,  1)
   elseif(truth_pos == 0) then
      mix_density = density_vc_PBM(d_a_par, c_par, vc,  0)
   else
      mix_density = rho*density_vc_PBM(d_a_par,c_par,vc,1) + (1.0_double-rho)*density_vc_PBM(d_a_par,c_par,vc,0)
   endif

 end function mix_density2

 real(kind=double) function d_pf_d_d_a2(vc) result(d_pf_d_d_a)
 ! function tht computes the derivative in d_a  of FPF or TPF, depending upon the truth of the category that is being
 ! analyzed
   real(kind=double), intent(IN):: vc

   if(truth_pos == 1) then
      d_pf_d_d_a =  d_tpf_d_d_a_PBM(d_a_par, c_par, vc)
   elseif(truth_pos == 0) then
      d_pf_d_d_a =  d_fpf_d_d_a_PBM(d_a_par, c_par, vc)
   else
      d_pf_d_d_a = rho*d_tpf_d_d_a_PBM(d_a_par,c_par,vc) + (1.0_double-rho)*d_fpf_d_d_a_PBM(d_a_par,c_par,vc)
   endif

 end function d_pf_d_d_a2

 real(kind=double) function d_pf_d_c2(vc) result(d_pf_d_c)
 ! function tht computes the derivative in c  of FPF or TPF, depending upon the truth of the category that is being
 ! analyzed
   real(kind=double), intent(IN):: vc


   if(truth_pos == 1) then
      d_pf_d_c =  d_tpf_d_c_PBM(d_a_par, c_par, vc)
   elseif(truth_pos == 0) then
      d_pf_d_c =  d_fpf_d_c_PBM(d_a_par, c_par, vc)
   else
      d_pf_d_c = rho*d_tpf_d_c_PBM(d_a_par,c_par,vc) + (1.0_double-rho)*d_fpf_d_c_PBM(d_a_par,c_par,vc)
   endif

 end function d_pf_d_c2

! --------------------------------------------------------------------------------------------
 end subroutine compute_FPF_TPF_median
! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------




! ---------------------------------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------
 subroutine compute_se_EV_logbeta(d_a_par, c_par, vc1, vc2, vc1_min, max_vc2, e_val_vc, e_val_vc2, &
                                 curr_cat, num_cat, cov, rho, d_vc_d_th1, se_log_beta)
! PURPOSE:  Compute the standard error of the expected value of log beta.
! ALGORITHM:It is based on the delta method and makes use of the fact that
!           the expected value of log(beta) can be written as E{log(beta)} = 2 a / (1 + b) E{v} - 2 c E{v^2} + log(b). Thus it is made of
!           3 terms , here called A, B, and C. The derivatives are built backwards by building those terms and summing them up.
!           The terms I and I2 are the integrals that make up the expected values and so DI and DI2 are their derivatives.
! NOTE:     The jacobian is 2*4 because we only need the derivatives of the sub-category cutoffs (vc1,vc2) in d_a, c, cat_ubound, cat_lbound.
!           the dependence upon d_a and c is voided because they are not computed as a function of d_a and c. However, the correlation that
!           is present in the inverse of the Expected Hessian (variance covariance matrix) is used when the delta method is applied.

 use proproc_functions, only: density_vc_PBM
 use proproc_functions, only: d_tpf_d_d_a_PBM,d_fpf_d_d_a_PBM ! derivatives of tpf and fpf
 use proproc_functions, only: d_tpf_d_c_PBM,d_fpf_d_c_PBM ! derivatives of tpf and fpf

 implicit none

 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: e_val_vc, e_val_vc2 ! expected values for vc and vc^2
 real(kind=double), intent(IN):: vc1, vc2 ! cutoffs that bound the current test result value
 logical, intent(IN):: vc1_min, max_vc2 ! whether one of the two cutoffs is the boundaries of the ROC space
 integer, intent(IN):: num_cat !  total number of categories
 integer, intent(IN):: curr_cat !  current category
 real(kind=double),  dimension(num_cat+1,num_cat+1), intent(IN) :: cov ! Variance-covariance matrix as estimated by the MLE
 real(kind=double), intent(IN) :: rho ! fraction of actually positive cases
 real(kind=double), intent(IN), dimension(4,2) :: d_vc_d_th1 ! Derivative of the lower bound (current_lbound) and upper (current_ubound)
                                                              ! around the current
                            ! test result value as a function of the MLE parameters, da, c, and the two cutoffs around the current
                            ! truth run (cat_lower_cutoff, cat_upper_cutoff). Variable names are the ones used by LP in July 2007.
 real(kind=double), intent(OUT)::  se_log_beta ! standard error of log beta of the mean

! Internal variables

 real(kind=double):: a, b ! CvBM parameters

 real(kind=double):: fpf1,fpf2, tpf1, tpf2 ! TPF and FPF for subinterval boundaries
 real(kind=double):: norm ! the integral of the density between the subinterval boundaries to compute the conditional expected values
 real(kind=double):: d_norm_d_d_a, d_norm_d_c, d_norm_d_vc1, d_norm_d_vc2 ! the derivativer of the integral of the density between
                          ! the subinterval boundaries to compute the conditional expected values
 real(kind=double):: k1, k1p, k1m, k2, k3 ! FPF/TPF coefficients

 real(kind=double):: coeff_d_a, coeff_c ! whether one of the two cutoffs is the boundaries of the ROC space
 real(kind=double):: d_k2_d_d_a, d_k2_d_c, d_k3_d_d_a, d_k3_d_c ! derivatives of the FPF/TPF coefficients in the parameters of the subinterval, note that
                     ! the derivatives of k1 aren't performed because they are trivial, and the ones in vc aren't necessary because the integrals as
                     ! I*(alpha1, alpha2, vc1, vc2) [ * = 1 or 2] are a function of the cutoffs and not of a function of them.
 real(kind=double):: d_Ev_d_d_a, d_Ev_d_c, d_Ev_d_vc1, d_Ev_d_vc2 ! derivatives of E(v) in the parameters of the subinterval
 real(kind=double):: d_Ev2_d_d_a, d_Ev2_d_c, d_Ev2_d_vc1, d_Ev2_d_vc2 ! derivatives of E(v2) in the parameters of the subinterval

 real(kind=double):: d_A_d_d_a, d_A_d_c, d_A_d_vc1, d_A_d_vc2 ! derivatives of 1st the three functions that sum up to E(log[beta])
 real(kind=double):: d_B_d_d_a, d_B_d_c, d_B_d_vc1, d_B_d_vc2 ! derivatives of 2nd the three functions that sum up to E(log[beta])
 real(kind=double):: d_C_d_d_a, d_C_d_c, d_C_d_vc1, d_C_d_vc2 ! derivatives of 3rd the three functions that sum up to E(log[beta])

 real(kind=double):: d_lb_d_vc1, d_lb_d_vc2, d_lb_d_d_a, d_lb_d_c ! derivatives of log(beta(median)) in the parameters of the interval
 real(kind=double):: d_logbeta_d_cut1, d_logbeta_d_cut2, d_logbeta_d_d_a, d_logbeta_d_c ! derivatives of log(beta(median)) in the parameters of the MLE estimation
                                           ! that affect this specific estimation.

 integer:: coeff ! Whether the current category is made of true positive or true negative cases or mixed


 !character(len = line_length):: msg    ! buffer for printing out results

 ! As a first step we check whether the category is made of only positives or only negatives or mixed.
 ! Mixed categories cannot be split, so the jacobian has to be the indentity matrix and the part we compute is
 ! 1 between the cutoffs and zero otherwise.
 if ( rho .speq. 1.0_double) then
       coeff = 1
 elseif ( rho .speq. 0.0_double) then
       coeff = -1
 else ! Mixed category
       coeff = 0
 endif

! Compute a and b
 b = (1.0_double + c_par) / (1.0_double - c_par)
 a = d_a_par * sqrt(1.0_double + b**2) / sqrt(2.0_double)

! Compute coefficients (which correspond to alpha1 and alpha2 in the integrals)
k2 = d_a_par * sqrt( 1.0_double + c_par**2) / 2.0_double
k3 = k2  / sign( max(abs(c_par), 1.0e-8_double), c_par) ! Max to avoid divisions by zero


! Compute the derivatives of the coefficients (which correspond to alpha1 and alpha2 in the integrals)
d_k2_d_d_a = sqrt(1.0_double + c_par**2)/ 2.0_double

d_k2_d_c = d_a_par * c_par/ (sqrt(1.0_double + c_par**2) * 2.0_double)

d_k3_d_d_a = d_k2_d_d_a / sign(  max( abs(c_par), 1.0e-8_double ),c_par) ! Max to avoid divisions by zero

d_k3_d_c = d_a_par / (sqrt(1.0_double + c_par**2) * 2.0_double) - &
           d_a_par * sqrt(1.0_double + c_par**2)/ ( 2.0_double * max( c_par**2, 1.0e-16_double ) ) ! Max to avoid divisions by zero


! Construct the derivaties of E(v) in the subinterval parameters. First take the derivative of the numerator of
! E(v) = Integral{vc1,vc2} f(vc|t) vc dvc / Integral{vc1,vc2} f(vc|t) dvc  that is
! Integral{vc1,vc2} f(vc|t) vc dvc . Considering that it is a ratio of two functions, use the rule for the
! derivatio of a ratio.

! We first compute the normalization constants, PF(vc1) - PF(vc2).
if(coeff == -1) then
      call fpf_PBM(d_a_par, c_par, vc1, fpf1)
      call fpf_PBM(d_a_par, c_par, vc2, fpf2)
      norm = fpf1 - fpf2
elseif(coeff == 1) then
      call tpf_PBM(d_a_par, c_par, vc1, tpf1)
      call tpf_PBM(d_a_par, c_par, vc2, tpf2)
      norm = tpf1 - tpf2
else
      call fpf_PBM(d_a_par, c_par, vc1, fpf1)
      call fpf_PBM(d_a_par, c_par, vc2, fpf2)
      call tpf_PBM(d_a_par, c_par, vc1, tpf1)
      call tpf_PBM(d_a_par, c_par, vc2, tpf2)
      norm = rho * ( tpf1 - tpf2) + (1.0_double - rho) * (fpf1 - fpf2)
endif



! in case vc1 is the first or last possible ROC point, then the vc value is a function of d_a and c
coeff_d_a =  sqrt(1.0_double + c_par**2) / (4.0_double*c_par)
coeff_c =  d_a_par/ ( 4.0_double*sqrt(1.0_double + c_par**2)) - &
                  d_a_par * sqrt(1.0_double + c_par**2) / (4.0_double*c_par**2)

! Then compute the derivative of the normalization constants, dPF(vc1)/dtheta - dPF(vc2)/dtheta, that is the denominator
! of the conditional expeced value Exp{v| vc1<vc<vc2} / Delta probability
if(coeff == -1 ) then ! Actually negative cases
       d_norm_d_d_a =  d_fpf_d_d_a_PBM(d_a_par, c_par, vc1 ) - d_fpf_d_d_a_PBM(d_a_par, c_par, vc2)
       d_norm_d_c   =  d_fpf_d_c_PBM(d_a_par, c_par, vc1) - d_fpf_d_c_PBM(d_a_par, c_par, vc2)
       d_norm_d_vc1 = -density_vc_PBM(d_a_par, c_par, vc1,  0)
       d_norm_d_vc2 = +density_vc_PBM(d_a_par, c_par, vc2,  0)

  elseif(coeff == 1) then
       d_norm_d_d_a =  d_tpf_d_d_a_PBM(d_a_par, c_par, vc1) - d_tpf_d_d_a_PBM(d_a_par, c_par, vc2)
       d_norm_d_c   =  d_tpf_d_c_PBM(d_a_par, c_par, vc1) - d_tpf_d_c_PBM(d_a_par, c_par, vc2)
       d_norm_d_vc1 = -density_vc_PBM(d_a_par, c_par, vc1,  1)
       d_norm_d_vc2 = +density_vc_PBM(d_a_par, c_par, vc2,  1)

else
       d_norm_d_d_a =  rho*  ( d_tpf_d_d_a_PBM(d_a_par,c_par,vc1)   - d_tpf_d_d_a_PBM(d_a_par,c_par,vc2)  ) + &
          (1.0_double - rho)*( d_fpf_d_d_a_PBM(d_a_par, c_par, vc1) - d_fpf_d_d_a_PBM(d_a_par, c_par, vc2))

       d_norm_d_c =  rho*(d_tpf_d_c_PBM(d_a_par,c_par,vc1) - d_tpf_d_c_PBM(d_a_par,c_par,vc2)  ) + &
          (1.0_double-rho)*(d_fpf_d_c_PBM(d_a_par,c_par,vc1) - d_fpf_d_c_PBM(d_a_par,c_par,vc2)  )

       d_norm_d_vc1 = -rho*density_vc_PBM(d_a_par,c_par,vc1,1)-(1.0_double-rho)*density_vc_PBM(d_a_par,c_par,vc1,0)
       d_norm_d_vc2 = +rho * density_vc_PBM(d_a_par,c_par,vc2,1) +(1.0_double - rho) * density_vc_PBM(d_a_par,c_par,vc2,0)

endif

!Compute the derivative of the numerator of E(v), d{Integral{vc1,vc2} f(vc|t) vc dvc}d{theta} .
if( coeff == 1 .or. coeff == -1 ) then
     ! the c_par appears with + for the positives and - for the negatives
     k1 = - (1.0_double + coeff * c_par)

     d_Ev_d_d_a =   d_Integral1_d_alpha2( k1 , coeff*k2, vc1, vc2) * coeff*d_k2_d_d_a  &
                  + d_Integral1_d_alpha2( k1 , k3      , vc1, vc2) * d_k3_d_d_a

     d_Ev_d_c =   d_Integral1_d_alpha1( k1 , coeff*k2 , vc1, vc2) * (-coeff) &
                + d_Integral1_d_alpha2( k1 , coeff*k2 , vc1, vc2) * coeff*d_k2_d_c &
                + d_Integral1_d_alpha1( k1 , k3       , vc1, vc2) * (-coeff) &
                + d_Integral1_d_alpha2( k1 , k3       , vc1, vc2) * d_k3_d_c

     d_Ev_d_vc1 = +(d_Integral1_d_vc( k1 , coeff*k2 , vc1) + d_Integral1_d_vc( k1 , k3 , vc1) )
     d_Ev_d_vc2 = -(d_Integral1_d_vc( k1 , coeff*k2 , vc2) + d_Integral1_d_vc( k1 , k3 , vc2) )


else
     k1p = - (1.0_double + c_par)
     k1m = - (1.0_double - c_par)

     d_Ev_d_d_a =   rho*(  d_Integral1_d_alpha2( k1p , k2 , vc1, vc2)*d_k2_d_d_a  &
                          +d_Integral1_d_alpha2( k1p , k3 , vc1, vc2)*d_k3_d_d_a ) + &
                  (1.0_double-rho)*( -d_Integral1_d_alpha2( k1m , -k2 , vc1, vc2)*d_k2_d_d_a  &
                                     +d_Integral1_d_alpha2( k1m ,  k3 , vc1, vc2)*d_k3_d_d_a )

     d_Ev_d_c =  rho*    ( - d_Integral1_d_alpha1( k1p ,  k2 , vc1, vc2) &
                           + d_Integral1_d_alpha2( k1p ,  k2 , vc1, vc2)*d_k2_d_c &
                           - d_Integral1_d_alpha1( k1p ,  k3 , vc1, vc2) &
                           + d_Integral1_d_alpha2( k1p ,  k3 , vc1, vc2)*d_k3_d_c ) + &
        (1.0_double-rho)*( + d_Integral1_d_alpha1( k1m , -k2 , vc1, vc2) &
                           - d_Integral1_d_alpha2( k1m , -k2 , vc1, vc2)*d_k2_d_c &
                           + d_Integral1_d_alpha1( k1m ,  k3 , vc1, vc2) &
                           + d_Integral1_d_alpha2( k1m ,  k3 , vc1, vc2)*d_k3_d_c )

     d_Ev_d_vc1 = rho*( d_Integral1_d_vc( k1p , k2 , vc1) + d_Integral1_d_vc( k1p , k3 , vc1) ) + &
                 (1.0_double-rho) * (d_Integral1_d_vc( k1m , -k2 , vc1) + d_Integral1_d_vc( k1m , k3 , vc1)  )

     d_Ev_d_vc2 = - rho*(d_Integral1_d_vc( k1p , k2 , vc2) + d_Integral1_d_vc( k1p , k3 , vc2))  &
                  - (1.0_double-rho) * (d_Integral1_d_vc( k1m , -k2 , vc2) + d_Integral1_d_vc( k1m , k3 , vc2)  )


endif


! Construct the full derivative of the ratio that produces the expected value. Notice that
! e_val_vc is the numerator divided by the denominator and this is why it needs to be multiplied by norm
! in order to obtain the numerator only

! write(*,"('DEvv0', 4(1x,d12.6))") d_Ev_d_d_a, d_Ev_d_c, d_Ev_d_vc1, d_Ev_d_vc2

 d_Ev_d_d_a =  (d_Ev_d_d_a  - e_val_vc * d_norm_d_d_a)/ norm
 d_Ev_d_c   =  (d_Ev_d_c    - e_val_vc * d_norm_d_c  )/ norm
 d_Ev_d_vc1 =  (d_Ev_d_vc1  - e_val_vc * d_norm_d_vc1)/ norm
 d_Ev_d_vc2 =  (d_Ev_d_vc2  - e_val_vc * d_norm_d_vc2)/ norm

 !write(*,"('DEvv1', 4(1x,d12.6))") d_Ev_d_d_a, d_Ev_d_c, d_Ev_d_vc1, d_Ev_d_vc2

! Construct the derivaties of E(v^2) in the subinterval parameters. First take the derivative of the numerator of
! E(v^2) = Integral{vc1,vc2} f(vc|t) vc^2 dvc / Integral{vc1,vc2} f(vc|t) dvc  that is
! Integral{vc1,vc2} f(vc|t) vc dvc . Again, considering that it is a ratio of two functions, use the rule for the
! derivatio of a ratio. The denominator and its derivatives are the same as for E(v).

!Compute the derivative of the numerator of E(v^2), d{Integral{vc1,vc2} f(vc|t) vc^2 dvc}d{theta} .
if( coeff == 1 .or. coeff == -1 ) then
     k1 = - (1.0_double + coeff * c_par)

     d_Ev2_d_d_a =   d_Integral2_d_alpha2( k1 , coeff*k2 , vc1, vc2)*coeff*d_k2_d_d_a  &
                   + d_Integral2_d_alpha2( k1 ,       k3 , vc1, vc2)*d_k3_d_d_a

     d_Ev2_d_c =     d_Integral2_d_alpha1( k1 , coeff*k2 , vc1, vc2)*(-coeff) &
                   + d_Integral2_d_alpha2( k1 , coeff*k2 , vc1, vc2)*coeff*d_k2_d_c &
                   + d_Integral2_d_alpha1( k1 , k3       , vc1, vc2)*(-coeff) &
                   + d_Integral2_d_alpha2( k1 , k3       , vc1, vc2)*d_k3_d_c

     d_Ev2_d_vc1 =  d_Integral2_d_vc( k1 , coeff*k2 , vc1) + d_Integral2_d_vc( k1 , k3 , vc1)
     d_Ev2_d_vc2 = -(d_Integral2_d_vc( k1 , coeff*k2 , vc2) + d_Integral2_d_vc( k1 , k3 , vc2))
else
     k1p = - (1.0_double + c_par)
     k1m = - (1.0_double - c_par)

     d_Ev2_d_d_a =     rho  *(   d_Integral2_d_alpha2( k1p ,  k2 , vc1, vc2)*d_k2_d_d_a  &
                                +d_Integral2_d_alpha2( k1p ,  k3 , vc1, vc2)*d_k3_d_d_a ) + &
            (1.0_double-rho)*( - d_Integral2_d_alpha2( k1m , -k2 , vc1, vc2)*d_k2_d_d_a  &
                               + d_Integral2_d_alpha2( k1m ,  k3 , vc1, vc2)*d_k3_d_d_a )

     d_Ev2_d_c =         rho*( - d_Integral2_d_alpha1( k1p , k2 , vc1, vc2) &
                               + d_Integral2_d_alpha2( k1p , k2 , vc1, vc2)*d_k2_d_c &
                               - d_Integral2_d_alpha1( k1p , k3 , vc1, vc2) &
                               + d_Integral2_d_alpha2( k1p , k3 , vc1, vc2)*d_k3_d_c ) + &
            (1.0_double-rho)*( + d_Integral2_d_alpha1( k1m , -k2 , vc1, vc2) &
                               - d_Integral2_d_alpha2( k1m , -k2 , vc1, vc2)*d_k2_d_c &
                               + d_Integral2_d_alpha1( k1m ,  k3 , vc1, vc2) &
                               + d_Integral2_d_alpha2( k1m ,  k3 , vc1, vc2)*d_k3_d_c )

     d_Ev2_d_vc1 =   rho        *(d_Integral2_d_vc( k1p , k2  , vc1) + d_Integral2_d_vc( k1p , k3 , vc1)) + &
                (1.0_double-rho)*(d_Integral2_d_vc( k1m , -k2 , vc1) + d_Integral2_d_vc( k1m , k3 , vc1)  )

     d_Ev2_d_vc2 = - rho          *(d_Integral2_d_vc( k1p ,  k2 , vc2) + d_Integral2_d_vc( k1p , k3 , vc2)   )  &
                - (1.0_double-rho)*(d_Integral2_d_vc( k1m , -k2 , vc2) + d_Integral2_d_vc( k1m , k3 , vc2)  )

endif


! Construct the full derivative of the ratio that produces the expected value of v^2. Notice that
! e_val_vc is the numerator divided by the denominator and this is why it needs to be multiplied by norm
! in order to obtain the numerator only


 d_Ev2_d_d_a =  (d_Ev2_d_d_a  - e_val_vc2 * d_norm_d_d_a)/ norm
 d_Ev2_d_c   =  (d_Ev2_d_c    - e_val_vc2 * d_norm_d_c  )/ norm
 d_Ev2_d_vc1 =  (d_Ev2_d_vc1  - e_val_vc2 * d_norm_d_vc1)/ norm
 d_Ev2_d_vc2 =  (d_Ev2_d_vc2  - e_val_vc2 * d_norm_d_vc2)/ norm

! construct the derivatives of the three terms that compose E{Log(beta()} in the interval parameters: d_A_d_da, d_B_d_d_a, d_c_d_d_a and so on....

d_A_d_d_a =  2.0_double/(1.0_double+b) * ( sqrt( (1.0_double + b**2)/ 2.0_double) * e_val_vc   +   a * d_Ev_d_d_a )

d_A_d_c = ( Sqrt(8.0_double)*d_a_par/(1.0_double-c_par)**2  ) * &
          (  b/( (1.0_double + b)*sqrt(1.0_double + b**2)) - sqrt(1.0_double + b**2)/(1.0_double + b)**2 )*e_val_vc + &
          (2.0_double * a / (1.0_double + b) ) * d_Ev_d_c

d_A_d_vc1 = ( 2.0_double * a / (1.0_double + b) ) * d_Ev_d_vc1

d_A_d_vc2 = ( 2.0_double * a / (1.0_double + b) ) * d_Ev_d_vc2

! write(*,"('EvA', 4(1x,d12.6))") d_A_d_d_a, d_A_d_c, d_A_d_vc1, d_A_d_vc2


d_B_d_d_a = - 2.0_double*c_par * d_Ev2_d_d_a
d_B_d_c   = - 2.0_double*c_par * d_Ev2_d_c   - 2.0_double*e_val_vc2
d_B_d_vc1 = - 2.0_double*c_par * d_Ev2_d_vc1
d_B_d_vc2 = - 2.0_double*c_par * d_Ev2_d_vc2

 !write(*,"('EvB', 4(1x,d12.6))") d_B_d_d_a, d_B_d_c, d_B_d_vc1, d_B_d_vc2


d_C_d_d_a = 0.0_double
d_C_d_c   = 2.0_double / (  b*(1.0_double - c_par)**2  )
d_C_d_vc1 = 0.0_double
d_C_d_vc2 = 0.0_double

! construct the derivatives of E{Log(beta()} in the interval parameters

 d_lb_d_vc1 = d_A_d_vc1 + d_B_d_vc1 + d_C_d_vc1
 d_lb_d_vc2 = d_A_d_vc2 + d_B_d_vc2 + d_C_d_vc2
 d_lb_d_d_a = d_A_d_d_a + d_B_d_d_a + d_C_d_d_a
 d_lb_d_c   = d_A_d_c   + d_B_d_c   + d_C_d_c


! NOTE THAT IT IS IMPOSSIBLE THAT THEY ARE BOTH THE FIRST AND LAST UNLESS
! THERE IS ONLY ONE CATEGORY, WHICH CAN'T BE ROC BECAUSE IT MAKES NO DISTINCTIONS
! Take care of when the first cutoff is the smallest possible value
if(vc1_min .and. (c_par < 0.0_double) ) then
     d_lb_d_d_a = d_lb_d_d_a + d_lb_d_vc1*coeff_d_a
     d_lb_d_c   = d_lb_d_c + d_lb_d_vc1*coeff_c
     d_lb_d_vc1 = 0.0_double
elseif(vc1_min .and. (c_par > -1.0e-8_double) ) then
     d_lb_d_vc1 = 0.0_double
endif

! Take care of when the last cutoff is the smallest possible value
if(max_vc2 .and. (c_par > 0.0_double) ) then
     d_lb_d_d_a = d_lb_d_d_a + d_lb_d_vc2*coeff_d_a
     d_lb_d_c   = d_lb_d_c + d_lb_d_vc2*coeff_c
     d_lb_d_vc2 = 0.0_double
elseif(max_vc2 .and. (c_par < +1.0e-8_double) ) then
     d_lb_d_vc2 = 0.0_double
endif

! write(*,"('Ev', 4(1x,d12.6))") d_lb_d_d_a, d_lb_d_c, d_lb_d_vc1, d_lb_d_vc2



! Construct the derivatives of Log(beta(median)) in the MLE parameters
! Note that if the subcategory cutoff is the boundary so is the category boundary and so it is zero too.
 d_logbeta_d_cut1 = (d_lb_d_vc1 *  d_vc_d_th1(3,1)  +  d_lb_d_vc2*d_vc_d_th1(3,2) )
 d_logbeta_d_cut2 = (d_lb_d_vc1 *  d_vc_d_th1(4,1)  +  d_lb_d_vc2*d_vc_d_th1(4,2) )

 d_logbeta_d_d_a  = d_lb_d_vc1 *  d_vc_d_th1(1,1)  +   d_lb_d_vc2*d_vc_d_th1(1,2) + d_lb_d_d_a
 d_logbeta_d_c    = d_lb_d_vc1 *  d_vc_d_th1(2,1)  +   d_lb_d_vc2*d_vc_d_th1(2,2) + d_lb_d_c

! Apply the delta method

! If it is the first category, then the cutoff before the current is a fixed point
! so it does not contribute to the derivation
if(curr_cat == 1) then
    se_log_beta = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_cut2**2 *  cov(curr_cat  + 2,curr_cat  + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   * cov(1,2)               + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_cut2 *  cov(1, curr_cat  + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_cut2 *  cov(2, curr_cat  + 2)   &
        )
! If it is the last category, then the cutoff after the current is a fixed point
! so it does not contribute to the derivation
elseif(curr_cat == num_cat) then
    se_log_beta = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_cut1**2 *  cov(curr_cat - 1 + 2, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   * cov(1,2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_cut1 *  cov(1, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_cut1 *  cov(2, curr_cat - 1 + 2) &
       )
else ! Any category in the middle has all the cutoffs
    se_log_beta = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_cut1**2 *  cov(curr_cat - 1 + 2, curr_cat - 1 + 2) + &
         d_logbeta_d_cut2**2 *  cov(curr_cat  + 2,curr_cat  + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   *  cov(1,2)                 + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_cut1 *  cov(1, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_cut2 *  cov(1, curr_cat  + 2)    + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_cut1 *  cov(2, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_cut2 *  cov(2, curr_cat  + 2)    + &
         2.0_double * d_logbeta_d_cut1 * d_logbeta_d_cut2 *  cov(curr_cat - 1 + 2, curr_cat  + 2) &
       )
 endif




! --------------------------------------------------------------------------------------------
 end subroutine compute_se_EV_logbeta
! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------
real(kind=double) function  d_Integral1_d_alpha1(alpha1, alpha2, vc1, vc2)
! --------------------------------------------------------------------------------------------
! PURPOSE: compute the derivative in alpha1 of the general integral type that is the basic ingredient to compute expected values
!          of v. Refer to the LLPesce formulas
! NOTE:    it deals roughly with preventing divisions by zero, the calling program should make sure that
!          alpha1 is not becoming zero.
 use statistic_functions, only : phi,g

 implicit none

 real(kind=double), intent(IN):: alpha1, alpha2 ! constants used to compute the integral
 real(kind=double), intent(IN):: vc1, vc2 ! integration limits
 real(kind=double):: gamma1, gamma2 ! arguments of special functions

 gamma1 = alpha1*vc1  + alpha2
 gamma2 = alpha1*vc2  + alpha2


 d_integral1_d_alpha1 = +g(gamma1)*vc1*(alpha2 - gamma1) - g(gamma2)*vc2*(alpha2 - gamma2) &
                       + Integral1(alpha1, alpha2, vc1, vc2)

 d_integral1_d_alpha1 = - d_integral1_d_alpha1 / min( alpha1, -1.0e-10_double ) ! Min to avoid divisions by zero, alpha1 is
                         ! always negative


! --------------------------------------------------------------------------------------------
   end function d_Integral1_d_alpha1
! --------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------
real(kind=double) function  d_Integral1_d_alpha2(alpha1, alpha2, vc1, vc2)
! --------------------------------------------------------------------------------------------
! PURPOSE: compute the derivative in alpha2 of the general integral type that is the basic ingredient to compute expected values
!          of v. Refer to the LLPesce formulas
! NOTE:    it deals roughly with preventing divisions by zero, the calling program should make sure that
!          alpha1 is not becoming zero.
 use statistic_functions, only : phi,g

 implicit none

 real(kind=double), intent(IN):: alpha1, alpha2 ! constants used to compute the integral
 real(kind=double), intent(IN):: vc1, vc2 ! integration limits
 real(kind=double):: gamma1, gamma2 ! arguments of special functions

 gamma1 = alpha1*vc1  + alpha2
 gamma2 = alpha1*vc2  + alpha2

 d_integral1_d_alpha2 = g(gamma1)*(alpha2 - gamma1) - g(gamma2)*(alpha2 - gamma2) + phi(gamma1) - phi(gamma2)

 d_integral1_d_alpha2 = - d_integral1_d_alpha2 / min( alpha1, -1.0e-10_double ) ! Max to avoid divisions by zero

! --------------------------------------------------------------------------------------------
   end function d_Integral1_d_alpha2
! --------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------
real(kind=double) function  d_Integral1_d_vc(alpha1, alpha2, vc)
! --------------------------------------------------------------------------------------------
! PURPOSE: compute the derivative in vc of the general integral type that is the basic ingredient to compute expected values
!          of v. Refer to the LLPesce formulas. If called for the first integration limit (usually called vc1) the
!          derivative should be used as it, if called for the second, it should be used with the minus sign.
! NOTE:    it deals roughly with preventing divisions by zero, the calling program should make sure that
!          alpha1 is not becoming zero.
! NOTE2:   I decided to use it with the plus and minus signs externally instead of making a logical switch because it seemed
!          too stupid to waste numerical efficiency like this here.

 use statistic_functions, only : phi,g

 implicit none

 real(kind=double), intent(IN):: alpha1, alpha2 ! constants used to compute the integral
 real(kind=double), intent(IN):: vc ! integration limit currently used (one of the two)
 real(kind=double):: gamma ! argument of special functions

 gamma = alpha1*vc  + alpha2

 d_integral1_d_vc = g(gamma)*(gamma - alpha2)


! --------------------------------------------------------------------------------------------
   end function d_Integral1_d_vc
! --------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------
real(kind=double) function  d_Integral2_d_alpha1(alpha1, alpha2, vc1, vc2)
! --------------------------------------------------------------------------------------------
! PURPOSE: compute the derivative in alpha1 of the general integral type that is the basic ingredient to compute expected values
!          of v^2. Refer to the LLPesce formulas
! NOTE:    it deals roughly with preventing divisions by zero, the calling program should make sure that
!          alpha1 is not becoming zero.
 use statistic_functions, only : phi,g

 implicit none

 real(kind=double), intent(IN):: alpha1, alpha2 ! constants used to compute the integral
 real(kind=double), intent(IN):: vc1, vc2 ! integration limits
 real(kind=double):: gamma1, gamma2 ! arguments of special functions

 gamma1 = alpha1*vc1  + alpha2
 gamma2 = alpha1*vc2  + alpha2


 d_integral2_d_alpha1 =   (alpha2**2 + 1.0_double) * ( g(gamma1)* vc1  -  g(gamma2)*vc2 ) &
              - g(gamma1)*vc1*(  1.0_double + gamma1*(alpha2 - alpha1* vc1)  ) &
              + g(gamma2)*vc2*(  1.0_double + gamma2*(alpha2 - alpha1* vc2)  )

 d_integral2_d_alpha1 = d_integral2_d_alpha1 / max( alpha1**2, 1.0e-20_double ) ! Max to avoid divisions by zero

 ! Add the derivative of the coefficient in dalpha1
 d_integral2_d_alpha1 =  d_integral2_d_alpha1 - ( 2.0_double / alpha1 ) *  Integral2(alpha1, alpha2, vc1, vc2)

! --------------------------------------------------------------------------------------------
   end function d_Integral2_d_alpha1
! --------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------
real(kind=double) function  d_Integral2_d_alpha2(alpha1, alpha2, vc1, vc2)
! --------------------------------------------------------------------------------------------
! PURPOSE: compute the derivative in alpha2 of the general integral type that is the basic ingredient to compute expected values
!          of v^2. Refer to the LLPesce formulas
! NOTE:    it deals roughly with preventing divisions by zero, the calling program should make sure that
!          alpha1 is not becoming zero.
 use statistic_functions, only : phi,g

 implicit none

 real(kind=double), intent(IN):: alpha1, alpha2 ! constants used to compute the integral
 real(kind=double), intent(IN):: vc1, vc2 ! integration limits
 real(kind=double):: gamma1, gamma2 ! arguments of special functions

 gamma1 = alpha1*vc1  + alpha2
 gamma2 = alpha1*vc2  + alpha2

 d_integral2_d_alpha2 =  2.0_double*alpha2 * (  phi(gamma1)  -  phi(gamma2) ) &
                        + (1.0_double + alpha2**2) * ( g(gamma1)  -  g(gamma2) ) &
                        + g(gamma1)*(  1.0_double - gamma1*(alpha2 - alpha1* vc1)  ) &
                        - g(gamma2)*(  1.0_double - gamma2*(alpha2 - alpha1* vc2)  )

 d_integral2_d_alpha2 = d_integral2_d_alpha2 / max( alpha1**2, 1.0e-20_double ) ! Max to avoid divisions by zero


! --------------------------------------------------------------------------------------------
   end function d_Integral2_d_alpha2
! --------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------
real(kind=double) function  d_Integral2_d_vc(alpha1, alpha2, vc)
! --------------------------------------------------------------------------------------------
! PURPOSE: compute the derivative in vc of the general integral type that is the basic ingredient to compute expected values
!          of v^2. Refer to the LLPesce formulas. If called for the first integration limit (usually called vc1) the
!          derivative should be used as it, if called for the second, it should be used with the minus sign.
! NOTE:    it deals roughly with preventing divisions by zero, the calling program should make sure that
!          alpha1 is not becoming zero.
! NOTE2:   I decided to use it with the plus and minus signs externally instead of making a logical switch because it seemed
!          too stupid to waste numerical efficiency like this here.

 use statistic_functions, only : phi,g

 implicit none

 real(kind=double), intent(IN):: alpha1, alpha2 ! constants used to compute the integral
 real(kind=double), intent(IN):: vc ! integration limit currently used (one of the two)
 real(kind=double):: gamma ! argument of special functions

 gamma = alpha1*vc  + alpha2

 d_integral2_d_vc =  (1.0_double + alpha2**2) *  g(gamma)*alpha1 &
                      - g(gamma)* alpha1*(  1.0_double + gamma*(alpha2 - alpha1* vc)  )

 d_integral2_d_vc = d_integral2_d_vc / max( alpha1**2, 1.0e-20_double ) ! Max to avoid divisions by zero


! --------------------------------------------------------------------------------------------
   end function d_Integral2_d_vc
! --------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------
 subroutine compute_se_logbeta_median(d_a_par, c_par, vc1, vc2, median_vc, &
                                 curr_cat, num_cat, cov, rho, d_vc_d_th1, se_log_beta)
! PURPOSE: This routine computes the Jacobian of the transformation from the variables used in the MLE, the truth runs cutoffs,
!          d_a and c, to the variables used in the computation of the latent variables quantities, which are the bounds around a
!          test value result, d_a and c. The second is a subset of the first.
!ALGORITHM: first one notices that only the derivatives of the bounds around the test result values are needed because d_a and c
!           appear in both representations and they don't depend upon the truth runs cutoffs. That part of the jacobian would be
!           the identity matrix and there is no need to compute it or store it.
!           Here we call G the function that implicitly defines the cutoffs around the truth run. G(vc*, cut1, cut2, d_a, c) =0
!           PF(vc*) - PF(cut1 = cat_ubound)*[1 - kij] - PF(cut2 = cat_lbound)*kij  ==0     ; where kij = the sum of all cases in this
!           truth run that relate to latent values smaller than vc*. We need this function to compute the derivatives implicitely because
!           we cannot invert the function analytically (or at least I did not find a way to do it).
! NOTE:     The jacobian is 2*4 because we only need the derivatives of the sub-category cutoffs (vc1,vc2) in d_a, c, cat_ubound, cat_lbound.
!           the dependence upon d_a and c is voided because they are not computed as a function of d_a and c. However, the correlation that
!           is present in the inverse of the Expected Hessian (variance covariance matrix) is used when the delta method is applied.

 use proproc_functions, only: density_vc_PBM
 use proproc_functions, only: d_tpf_d_d_a_PBM,d_fpf_d_d_a_PBM ! derivatives of tpf and fpf
 use proproc_functions, only: d_tpf_d_c_PBM,d_fpf_d_c_PBM ! derivatives of tpf and fpf


 implicit none

 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: median_vc ! median vc value
 real(kind=double), intent(IN):: vc1, vc2 ! cutoffs that bound the current test result value
 integer, intent(IN):: num_cat !  total number of categories
 integer, intent(IN):: curr_cat !  current category
 real(kind=double),  dimension(num_cat+1,num_cat+1), intent(IN) :: cov ! Variance-covariance matrix as estimated by the MLE
 real(kind=double) :: rho ! fraction of actually positive cases
 real(kind=double), intent(IN), dimension(4,2) :: d_vc_d_th1 ! Derivative of the lower bound (current_lbound) and upper (current_ubound)
                                                              ! around the current
                            ! test result value as a function of the MLE parameters, da, c, and the two cutoffs around the current
                            ! truth run (cat_lower_cutoff, cat_upper_cutoff). Variable names are the ones used by LP in July 2007.
 real(kind=double), intent(OUT)::  se_log_beta ! standard error of log beta of the mean
! Internal variables
 real(kind=double):: d_G_d_m ! deriviative of the implicit function, see under "ALGORITHM"
 real(kind=double):: d_G_d_vc1, d_G_d_vc2 ! deriviative of the implicit function, see under "ALGORITHM"
 real(kind=double):: d_G_d_d_a, d_G_d_c! deriviative of the implicit function, see under "ALGORITHM"

 real(kind=double):: d_m_d_vc1, d_m_d_vc2, d_m_d_d_a, d_m_d_c ! derivatives of m in the parameters of the interval

 real(kind=double):: d_lb_dm ! Derivative of logbeta in the median
 real(kind=double):: a, b ! CvBM parameters

 real(kind=double):: d_lb_d_d_a_m_fix, d_lb_d_c_m_fix ! Derivative of logbeta in the median, keeping the median fixed

 real(kind=double):: d_lb_d_vc1, d_lb_d_vc2, d_lb_d_d_a, d_lb_d_c ! derivatives of log(beta(median)) in the parameters of the interval
 real(kind=double):: d_logbeta_d_cut1, d_logbeta_d_cut2, d_logbeta_d_d_a, d_logbeta_d_c ! derivatives of log(beta(median)) in the parameters of the MLE estimation
                                           ! that affect this specific estimation.

 integer:: truth_pos ! Whether the current category is made of true positive or true negative cases.

 ! As a first step we check whether the category is made of only positives or only negatives or mixed.
 ! Mixed categories cannot be split, so the jacobian has to be the indentity matrix and the part we compute is
 ! 1 between the cutoffs and zero otherwise.
 if ( rho .speq. 1.0_double) then
       truth_pos = 1
 elseif ( rho .speq. 0.0_double) then
       truth_pos = 0
 else ! Mixed category
       truth_pos = 2
 endif

! Now we use the implicit function G to compute  the derivatives of the median in the parameters
d_G_d_m    =  - mix_density(median_vc)
d_G_d_vc1 =  0.5_double * mix_density(vc1)
d_G_d_vc2 =  0.5_double * mix_density(vc2)
d_G_d_d_a  =  d_pf_d_d_a(median_vc) - 0.5_double*d_pf_d_d_a(vc1) - 0.5_double*d_pf_d_d_a(vc2)
d_G_d_c    =  d_pf_d_c(median_vc)   - 0.5_double*d_pf_d_c(vc1)   - 0.5_double*d_pf_d_c(vc2)

d_m_d_vc1 = - d_G_d_vc1 / d_G_d_m
d_m_d_vc2 = - d_G_d_vc2 / d_G_d_m
d_m_d_d_a = - d_G_d_d_a / d_G_d_m
d_m_d_c   = - d_G_d_c / d_G_d_m


! Compute a and b
 b = (1.0_double + c_par) / (1.0_double - c_par)
 a = d_a_par * sqrt(1.0_double + b**2) / sqrt(2.0_double)

! Compute the derivative of the the log(beta) function ( 2 a m / (1+b) * median_vc - 2 c median_vc^2 ) in the median

d_lb_dm = 2.0_double * a / (1.0_double + b) - 4.0_double*c_par*median_vc

! Compute the derivative of the the log(beta) function ( 2 a m / (1+b) * median_vc - 2 c median_vc^2 ) id d_a and c, keeping
! the median constant

d_lb_d_d_a_m_fix = ( sqrt( 2.0_double * (1.0_double + b**2)) / ( 1.0_double + b)  ) * median_vc

d_lb_d_c_m_fix  = ( 2.0_double * sqrt(2.0_double)* d_a_par * median_vc / (1.0_double - c_par)**2   ) *    &
                    (     b / ( (1.0_double + b) * sqrt(1.0_double + b**2) )                              &
                         -sqrt(1.0_double + b**2) / (1.0_double + b)**2  )                                &
                  - 2.0_double * median_vc**2 + 2.0_double / (b*(1.0_double - c_par)**2 )

 ! construct the derivatives of Log(beta(median)) in the interval parameters

 d_lb_d_vc1 = d_lb_dm * d_m_d_vc1
 d_lb_d_vc2 = d_lb_dm * d_m_d_vc2
 d_lb_d_d_a = d_lb_dm * d_m_d_d_a + d_lb_d_d_a_m_fix
 d_lb_d_c   = d_lb_dm * d_m_d_c + d_lb_d_c_m_fix

! write(*,"('M', 4(1x,d12.6))") d_lb_d_d_a, d_lb_d_c, d_lb_d_vc1, d_lb_d_vc2

! Construct the derivatives of Log(beta(median)) in the MLE parameters

 d_logbeta_d_cut1 = d_lb_d_vc1 *  d_vc_d_th1(3,1)  +  d_lb_d_vc2*d_vc_d_th1(3,2)
 d_logbeta_d_cut2 = d_lb_d_vc1 *  d_vc_d_th1(4,1)  +  d_lb_d_vc2*d_vc_d_th1(4,2)

 d_logbeta_d_d_a  = d_lb_d_vc1 *  d_vc_d_th1(1,1)  +   d_lb_d_vc2*d_vc_d_th1(1,2) + d_lb_d_d_a
 d_logbeta_d_c    = d_lb_d_vc1 *  d_vc_d_th1(2,1)  +   d_lb_d_vc2*d_vc_d_th1(2,2) + d_lb_d_c



! Apply the delta method

! If it is the first category, then the cutoff before the current is a fixed point
! so it does not contribute to the derivation
if(curr_cat == 1) then
    se_log_beta = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_cut2**2 *  cov(curr_cat  + 2,curr_cat  + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   * cov(1,2)               + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_cut2 *  cov(1, curr_cat  + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_cut2 *  cov(2, curr_cat  + 2)   &
        )
! If it is the last category, then the cutoff after the current is a fixed point
! so it does not contribute to the derivation
elseif(curr_cat == num_cat) then
    se_log_beta = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_cut1**2 *  cov(curr_cat - 1 + 2, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   * cov(1,2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_cut1 *  cov(1, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_cut1 *  cov(2, curr_cat - 1 + 2) &
       )
else ! Any category in the middle has all the cutoffs
    se_log_beta = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_cut1**2 *  cov(curr_cat - 1 + 2, curr_cat - 1 + 2) + &
         d_logbeta_d_cut2**2 *  cov(curr_cat  + 2,curr_cat  + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   *  cov(1,2)                 + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_cut1 *  cov(1, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_cut2 *  cov(1, curr_cat  + 2)    + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_cut1 *  cov(2, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_cut2 *  cov(2, curr_cat  + 2)    + &
         2.0_double * d_logbeta_d_cut1 * d_logbeta_d_cut2 *  cov(curr_cat - 1 + 2, curr_cat  + 2) &
       )
 endif



 contains

 real(kind=double) function mix_density(vc)
 ! function tht computes the density as a function of vc, depending upon the truth of the category that is being
 ! analyzed
   real(kind=double), intent(IN):: vc

   if(truth_pos == 1) then
      mix_density = density_vc_PBM(d_a_par, c_par, vc,  1)
   elseif(truth_pos == 0) then
      mix_density = density_vc_PBM(d_a_par, c_par, vc,  0)
   else
      mix_density = rho*density_vc_PBM(d_a_par,c_par,vc, 1) + (1.0_double-rho)*density_vc_PBM(d_a_par,c_par,vc,0)
   endif

 end function mix_density

 real(kind=double) function d_pf_d_d_a(vc)
 ! function tht computes the derivative in d_a  of FPF or TPF, depending upon the truth of the category that is being
 ! analyzed
   real(kind=double), intent(IN):: vc

   if(truth_pos == 1) then
      d_pf_d_d_a =  d_tpf_d_d_a_PBM(d_a_par, c_par, vc)
   elseif(truth_pos == 0) then
      d_pf_d_d_a =  d_fpf_d_d_a_PBM(d_a_par, c_par, vc)
   else
      d_pf_d_d_a = rho*d_tpf_d_d_a_PBM(d_a_par,c_par,vc) + (1.0_double-rho)*d_fpf_d_d_a_PBM(d_a_par,c_par,vc)
   endif

 end function d_pf_d_d_a

 real(kind=double) function d_pf_d_c(vc)
 ! function tht computes the derivative in c  of FPF or TPF, depending upon the truth of the category that is being
 ! analyzed
   real(kind=double), intent(IN):: vc


   if(truth_pos == 1) then
      d_pf_d_c =  d_tpf_d_c_PBM(d_a_par, c_par, vc)
   elseif(truth_pos == 0) then
      d_pf_d_c =  d_fpf_d_c_PBM(d_a_par, c_par, vc)
   else
      d_pf_d_c = rho*d_tpf_d_c_PBM(d_a_par,c_par,vc) + (1.0_double-rho)*d_fpf_d_c_PBM(d_a_par,c_par,vc)
   endif

 end function d_pf_d_c




! --------------------------------------------------------------------------------------------
 end subroutine compute_se_logbeta_median
! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------
 subroutine compute_th0_th1_jacobian(d_a_par, c_par, cat_lbound, cat_ubound, vc1, vc2, kijm1, kij, rho, d_vc_d_th1)
! PURPOSE: This routine computes the Jacobian of the transformation from the variables used in the MLE, the truth runs cutoffs,
!          d_a and c, to the variables used in the computation of the latent variables quantities, which are the bounds around a
!          test value result, d_a and c. The second is a subset of the first.
!ALGORITHM: first one notices that only the derivatives of the bounds around the test result values are needed because d_a and c
!           appear in both representations and they don't depend upon the truth runs cutoffs. That part of the jacobian would be
!           the identity matrix and there is no need to compute it or store it.
!           Here we call G the function that implicitly defines the cutoffs around the truth run. G(vc*, cut1, cut2, d_a, c) =0
!           PF(vc*) - PF(cut1 = cat_ubound)*[1 - kij] - PF(cut2 = cat_lbound)*kij  ==0     ; where kij = the sum of all cases in this
!           truth run that relate to latent values smaller than vc*. We need this function to compute the derivatives implicitely because
!           we cannot invert the function analytically (or at least I did not find a way to do it).
! NOTE:     The jacobian is 2*4 because we only need the derivatives of the sub-category cutoffs (vc1,vc2) in d_a, c, cat_ubound, cat_lbound.
!           the dependence upon d_a and c is voided because they are not computed as a function of d_a and c. However, the correlation that
!           is present in the inverse of the Expected Hessian (variance covariance matrix) is used when the delta method is applied.

 use proproc_functions, only: density_vc_PBM
 use proproc_functions, only: d_tpf_d_d_a_PBM,d_fpf_d_d_a_PBM ! derivatives of tpf and fpf
 use proproc_functions, only: d_tpf_d_c_PBM,d_fpf_d_c_PBM ! derivatives of tpf and fpf


 implicit none

 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: cat_lbound, cat_ubound ! cutoffs that bound the truth run where the current test result value lies
 real(kind=double), intent(IN):: vc1, vc2 ! cutoffs that bound the current test result value
 real(kind=double), intent(IN):: kijm1, kij ! fraction of cases smaller thant the current current lower and upper bounds.
 real(kind=double) :: rho ! fraction of actually positive cases
 real(kind=double), intent(OUT), dimension(4,2) :: d_vc_d_th1 ! Derivative of the lower bound (current_lbound) and upper (current_ubound) around the current
                            ! test result value as a function of the MLE parameters, da, c, and the two cutoffs around the current
                            ! truth run (cat_lower_cutoff, cat_upper_cutoff). Variable names are the ones used by LP in July 2007.
 ! Internal variables
 real(kind=double):: d_G_d_vc1, d_G_d_vc2 ! deriviative of the implicit function, see under "ALGORITHM"
 real(kind=double):: d_G_d_cut1, d_G_d_cut2 ! deriviative of the implicit function, see under "ALGORITHM"
 real(kind=double):: d_G_d_d_a, d_G_d_c! deriviative of the implicit function, see under "ALGORITHM"

 integer:: truth_pos ! Whether the current category is made of true positive or true negative cases.

 ! As a first step we check whether the category is made of only positives or only negatives or mixed.
 ! Mixed categories cannot be split, so the jacobian has to be the indentity matrix and the part we compute is
 ! 1 between the cutoffs and zero otherwise.
 if ( rho .speq. 1.0_double) then
       truth_pos = 1
 elseif ( rho .speq. 0.0_double) then
       truth_pos = 0
 else
       d_vc_d_th1(:,1) = (/ 0.0_double, 0.0_double, 1.0_double, 0.0_double /)
       d_vc_d_th1(:,2) = (/ 0.0_double, 0.0_double, 0.0_double, 1.0_double /)
       return
 endif

! Now we use the implicit function G to compute  starting with vc1
d_G_d_vc1  =                        - density_vc_PBM(d_a_par, c_par, vc1,  truth_pos)
d_G_d_cut1 =  (1.0_double - kijm1) *  density_vc_PBM(d_a_par, c_par, cat_lbound,  truth_pos)
d_G_d_cut2 =   kijm1               *  density_vc_PBM(d_a_par, c_par, cat_ubound,  truth_pos)
d_G_d_d_a  =  d_pf_d_d_a(vc1) - (1.0_double - kijm1)*d_pf_d_d_a(cat_lbound) - kijm1*d_pf_d_d_a(cat_ubound)
d_G_d_c    =  d_pf_d_c(vc1)   - (1.0_double - kijm1)*d_pf_d_c(cat_lbound)   - kijm1*d_pf_d_c(cat_ubound)

d_vc_d_th1(1,1) = - d_G_d_d_a / d_G_d_vc1
d_vc_d_th1(2,1) = - d_G_d_c / d_G_d_vc1
d_vc_d_th1(3,1) = - d_G_d_cut1 / d_G_d_vc1
d_vc_d_th1(4,1) = - d_G_d_cut2 / d_G_d_vc1


! Now we use the implicit function G to compute vc2
d_G_d_vc2  =                      - density_vc_PBM(d_a_par, c_par, vc2,  truth_pos)
d_G_d_cut1 =  (1.0_double - kij) *  density_vc_PBM(d_a_par, c_par, cat_lbound,  truth_pos)
d_G_d_cut2 =   kij               *  density_vc_PBM(d_a_par, c_par, cat_ubound,  truth_pos)
d_G_d_d_a  =  d_pf_d_d_a(vc2) - (1.0_double - kij)*d_pf_d_d_a(cat_lbound) - kij*d_pf_d_d_a(cat_ubound)
d_G_d_c    =  d_pf_d_c(vc2)   - (1.0_double - kij)*d_pf_d_c(cat_lbound)   - kij*d_pf_d_c(cat_ubound)

d_vc_d_th1(1,2) = - d_G_d_d_a / d_G_d_vc2
d_vc_d_th1(2,2) = - d_G_d_c / d_G_d_vc2
d_vc_d_th1(3,2) = - d_G_d_cut1 / d_G_d_vc2
d_vc_d_th1(4,2) = - d_G_d_cut2 / d_G_d_vc2


 contains
 real(kind=double) function d_pf_d_d_a(vc)
 ! function tht computes the derivative in d_a  of FPF or TPF, depending upon the truth of the category that is being
 ! analyzed
   real(kind=double), intent(IN):: vc

   if(truth_pos==1) then
      d_pf_d_d_a =  d_tpf_d_d_a_PBM(d_a_par, c_par, vc)
   else
      d_pf_d_d_a =  d_fpf_d_d_a_PBM(d_a_par, c_par, vc)
   endif

 end function d_pf_d_d_a

 real(kind=double) function d_pf_d_c(vc)
 ! function tht computes the derivative in c  of FPF or TPF, depending upon the truth of the category that is being
 ! analyzed
   real(kind=double), intent(IN):: vc


   if(truth_pos==1) then
      d_pf_d_c =  d_tpf_d_c_PBM(d_a_par, c_par, vc)
   else
      d_pf_d_c =  d_fpf_d_c_PBM(d_a_par, c_par, vc)
   endif

 end function d_pf_d_c


 end subroutine compute_th0_th1_jacobian


 ! --------------------------------------------------------------------------------------------
 real(kind=double)  function se_logbeta_median(d_a_par, c_par, vc1, vc2, median_vc_j, cov, rho, curr_cat, num_cat, kij)
 ! --------------------------------------------------------------------------------------------
 ! PURPOSE: estimate the standard error of log beta of the median value within the interval chosen
 ! ALGORITHM: uses the delta method computing the derivatives of d(log[beta])/d(th) (th = d_a, c, vc1, vc2)
 !            using the fact that log(beta)(m) = f(m(ths), ths) = 2 a m / (1 + b) - 2 c m^2 ,
 !            where m is the median and thetas is the vector of the estimated parameters. then
 !            d f(m(ths),ths) / d(th_i) =  d f  |            * d m(ths)    + d f    |
 !                                         d m  | m = m(ths)   d th_i        d th_i | m = m(ths)
 !            With this equation all the terms are known, since we can obtain d Median / d th_i using the
 !            derivation rules for implicit functions on the definition of median.

 use proproc_functions, only: density_vc_PBM
 use proproc_functions, only: d_tpf_d_d_a_PBM,d_fpf_d_d_a_PBM ! derivatives of tpf and fpf
 use proproc_functions, only: d_tpf_d_c_PBM,d_fpf_d_c_PBM ! derivatives of tpf and fpf
 implicit none

 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: vc1, vc2 ! cutoffs that bound the CATEGORY where the value of this subcategory "falls"
 real(kind=double) :: median_vc_j !  vc value for median
 integer, intent(IN):: curr_cat ! the index of the current category
 integer, intent(IN):: num_cat !  total number of categories
 real(kind=double),  dimension(num_cat+1,num_cat+1), intent(IN) :: cov ! Variance-covariance matrix as estimated by the MLE
 real(kind=double) :: rho ! fraction of actually positive cases
 real(kind=double) :: kij ! this number represents the fraction of cases that are smaller the median of this subcategory. Some
                          ! subcategories are closer to one cutoff some closer to the other one so the value is not always
                          ! .5, but it always splits the cases that are in this category

 real(kind=double):: a, b ! the original parameters a and b for the binomal distribution

 real(kind=double):: d_g_d_m_j ! derivative of the median function ( PF1 + PF2 - 2PF_Median = 0)

 real(kind=double):: d_f_d_m_j ! derivative of log(beta(m) ) in d(m), where m is the median

 real(kind=double):: d_f_d_d_a ! derivative of log(beta(m) ) in d(d_a), keeping m constant.
 real(kind=double):: d_f_d_c ! derivative of log(beta(m) ) in d(c), keeping m constant.

 real(kind=double):: d_m_j_d_vc1 ! derivative of m  in d(vc1)
 real(kind=double):: d_m_j_d_vc2 ! derivative of m  in d(vc2)
 real(kind=double):: d_m_j_d_d_a, d_m_j_d_c ! derivative of m  in d_a and c

 real(kind=double):: d_logbeta_d_d_a, d_logbeta_d_c, d_logbeta_d_vc1, d_logbeta_d_vc2


 ! Compute a and b
 b = (1.0_double + c_par) / (1.0_double - c_par)
 a = d_a_par * sqrt(1.0_double + b**2) / sqrt(2.0_double)

 ! compute the derivative of log(beta) in terms of the median
 d_f_d_m_j = 2.0_double*a / (1.0_double + b) - 4 * c_par * median_vc_j

 ! Compute the derivatives of log(beta) as a function of the parameters, assuming median constant
 ! in this case d[log(beta)]/d(vc) = 0 for every vc
 d_f_d_d_a = sqrt(2.0_double * ( 1.0_double + b**2) ) * median_vc_j / ( 1.0_double + b)

 d_f_d_c =  ( 2.0_double * d_a_par * median_vc_j / ( 1.0_double - c_par)**2  ) * &
                 (&
                    b * sqrt(2.0_double) / ( Sqrt(1.0_double + b**2)* ( 1.0_double + b) ) - &
                    sqrt(2.0_double*(1.0_double+b**2) )/(1.0_double + b)**2                 &
                 ) - &
            2.0_double*median_vc_j**2 + 2.0_double / ( b*(1.0_double - c_par)**2 )


 ! Now compute the derivative of the implicit function g = PF(vc_i-1) + PF(vc_i) - 2 PF(Median) = 0
 ! to compute the derivatives of d(median)/d(theta) using the implicit derivation rules.

 d_g_d_m_j =  - (  rho              * density_vc_PBM(d_a_par, c_par, median_vc_j,  1) + &
                   (1.0_double-rho) * density_vc_PBM(d_a_par, c_par, median_vc_j, 0) &
                 )


 ! Compute d(median_vc_j) / d(parameter) = - dg/d(par) / dg/d(m)
 d_m_j_d_vc1 = - (1.0_double - kij) * ( &
                         rho                * density_vc_PBM(d_a_par, c_par, vc1, 1) + &
                         (1.0_double - rho) * density_vc_PBM(d_a_par, c_par, vc1, 0) &
                      ) / d_g_d_m_j


 d_m_j_d_vc2 = - kij * ( &
                         rho              * density_vc_PBM(d_a_par, c_par, vc2, 1) + &
                         (1.0_double-rho) * density_vc_PBM(d_a_par, c_par, vc2, 0)  &
                  ) / d_g_d_m_j


 d_m_j_d_d_a =  - ( &
              rho * ( - (1.0_double - kij) * d_tpf_d_d_a_PBM(d_a_par, c_par, vc1)                                &
                      -  kij               * d_tpf_d_d_a_PBM(d_a_par, c_par, vc2)                                &
                      +                      d_tpf_d_d_a_PBM(d_a_par, c_par, median_vc_j)                        &
                     ) +                                                                           &
              (1.0_double - rho) * (  - (1.0_double - kij) * d_fpf_d_d_a_PBM(d_a_par, c_par, vc1)                    &
                                      - kij                * d_fpf_d_d_a_PBM(d_a_par, c_par, vc2)                    &
                                      +                      d_fpf_d_d_a_PBM(d_a_par, c_par, median_vc_j)            &
                    )  &
                   ) /  d_g_d_m_j

 d_m_j_d_c =  - ( &
              rho * (  - (1.0_double - kij) * d_tpf_d_c_PBM(d_a_par, c_par, vc1)                                &
                       - kij                * d_tpf_d_c_PBM(d_a_par, c_par, vc2)                                &
                       +                      d_tpf_d_c_PBM(d_a_par, c_par, median_vc_j)                        &
                    ) +                                                                         &
              (1.0_double - rho) * ( - (1.0_double - kij) * d_fpf_d_c_PBM(d_a_par, c_par, vc1)                    &
                                     - kij                * d_fpf_d_c_PBM(d_a_par, c_par, vc2)                    &
                                     +                      d_fpf_d_c_PBM(d_a_par, c_par, median_vc_j)            &
                    )  &
                   ) /  d_g_d_m_j

! Apply the formula reported on the algorithm section to obtain the derivates of log(beta) in the
! parameters. Remember that d_f_d_vc = 0 for any vc because it does not appear explicitely in the
! formula to compute log(beta).

d_logbeta_d_d_a = d_f_d_m_j * d_m_j_d_d_a  + d_f_d_d_a
d_logbeta_d_c   = d_f_d_m_j * d_m_j_d_c    + d_f_d_c
d_logbeta_d_vc1 = d_f_d_m_j * d_m_j_d_vc1
d_logbeta_d_vc2 = d_f_d_m_j * d_m_j_d_vc2


! Apply the delta method

! If it is the first category, then the cutoff before the current is a fixed point
! so it does not contribute to the derivation
if(curr_cat == 1) then
    se_logbeta_median = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_vc2**2 *  cov(curr_cat  + 2,curr_cat  + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   * cov(1,2)               + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_vc2 *  cov(1, curr_cat  + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_vc2 *  cov(2, curr_cat  + 2)   &
        )
! If it is the last category, then the cutoff after the current is a fixed point
! so it does not contribute to the derivation
elseif(curr_cat == num_cat) then
    se_logbeta_median = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_vc1**2 *  cov(curr_cat - 1 + 2, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   * cov(1,2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_vc1 *  cov(1, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_vc1 *  cov(2, curr_cat - 1 + 2) &
       )
else ! Any category in the middle has all the cutoffs
    se_logbeta_median = sqrt( &
         d_logbeta_d_d_a**2 *  cov(1,1) + &
         d_logbeta_d_c**2   *  cov(2,2) + &
         d_logbeta_d_vc1**2 *  cov(curr_cat - 1 + 2, curr_cat - 1 + 2) + &
         d_logbeta_d_vc2**2 *  cov(curr_cat  + 2,curr_cat  + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_c   *  cov(1,2)                 + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_vc1 *  cov(1, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_d_a * d_logbeta_d_vc2 *  cov(1, curr_cat  + 2)    + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_vc1 *  cov(2, curr_cat - 1 + 2) + &
         2.0_double * d_logbeta_d_c   * d_logbeta_d_vc2 *  cov(2, curr_cat  + 2)    + &
         2.0_double * d_logbeta_d_vc1 * d_logbeta_d_vc2 *  cov(curr_cat - 1 + 2, curr_cat  + 2) &
       )
 endif

 ! --------------------------------------------------------------------------------------------
 end function se_logbeta_median
 ! --------------------------------------------------------------------------------------------




 ! --------------------------------------------------------------------------------------------
 real(kind=double) pure function logbeta_vc(d_a_par, c_par, v)
 ! --------------------------------------------------------------------------------------------
 ! PURPOSE: compute the values of log likelihood ratio given a value of v and values of d_a_par
 !          and c_par. Tramnsformations from Metz and Pan
 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: v ! cutoffs in the parametrization used for PROPROC2

 real(kind=double):: y ! variable y (see Metz and Pan)
 real(kind=double):: a, b ! the original parameters a and b for the binomal distribution


 ! Compute a and b
 b = (1.0_double + c_par) / (1.0_double - c_par)
 a = d_a_par * sqrt(1.0_double + b**2) / sqrt(2.0_double)
 ! compute y
 y = -2.0_double * ( -a*v + (b - 1.0_double)*v**2 ) / (1.0_double + b)
 ! compute log likelihood ration
 logbeta_vc = y + log(b)

 ! --------------------------------------------------------------------------------------------
 end function logbeta_vc
 ! --------------------------------------------------------------------------------------------


 ! --------------------------------------------------------------------------------------------
 subroutine compute_median_vc(d_a_par, c_par, vc_lower, vc_upper, rho, median_vc)
 ! --------------------------------------------------------------------------------------------
 !PURPOSE: compute the median of vc between an interval vc_lower and vc_upper

 use gen_numerics, only: zbren !  root finder
 use proproc_median_calc ! data that will be passed to the function that computes the medians
 implicit none

 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: vc_lower, vc_upper ! cutoff values between which the integration has to be
                                  ! performed to compute the expected value
 real(kind=double), intent(IN):: rho ! prevalence (fraction of positive cases)
 real(kind=double), intent(OUT):: median_vc ! your guess

 real(kind=double):: fpf1,fpf2, tpf1, tpf2 ! TPF and FPF at boundary
 real(kind=double):: target_val ! TPF and FPF at boundary

 real(kind=double), parameter :: tol = 1.0e-6_double ! tolerace of numerical calculation
 integer:: ierror ! computation error


 ! Load the parameters for the calculation of the median into the module proproc_median_calc
 d_a = d_a_par
 c   = c_par
 r   = rho

! Computethe target values, the errors created when rho = 0 or 1 should be negligible
 call fpf_PBM(d_a_par, c_par, vc_lower, fpf1)
 call fpf_PBM(d_a_par, c_par, vc_upper, fpf2)
 call tpf_PBM(d_a_par, c_par, vc_lower, tpf1)
 call tpf_PBM(d_a_par, c_par, vc_upper, tpf2)

 target_val = 0.5_double*( (1.0_double - rho) *(fpf1 + fpf2)  +   rho*(tpf1 + tpf2)     )

 call  zbren( target_val , median_func , vc_upper,  vc_lower,  tol, median_vc , ierror)

 ! --------------------------------------------------------------------------------------------
 end subroutine compute_median_vc
 ! --------------------------------------------------------------------------------------------

 ! --------------------------------------------------------------------------------------------
 real(kind=double) function  median_func(vc)
 ! --------------------------------------------------------------------------------------------
 ! PURPOSE: Compute the function needed to determine the median cutoff value, it is a weighted
 !          mean of positive and negative distributions.
 ! NOTE:    The if clause is done to avoid unstable values and useless computations, given that
 !          the computations of FPF and TPF is somewhat expensive and this on is nested into an
 !          iterative procedur.
 use proproc_median_calc ! data that will be passed to the function that computes the medians

 implicit none
 real(kind=double), intent(IN) :: vc
 real(kind=double):: fpf, tpf ! TPF and FPF at boundary


 if(r .speq. 0.0_double) then
        call fpf_PBM(d_a, c, vc, fpf)
        median_func =  fpf

 elseif(r .speq. 1.0_double) then
        call tpf_PBM(d_a, c, vc, tpf)
        median_func =   tpf
 else
        call fpf_PBM(d_a, c, vc, fpf)
        call tpf_PBM(d_a, c, vc, tpf)
        median_func = (1.0_double - r)* fpf + r * tpf
 endif


 ! --------------------------------------------------------------------------------------------
 end function  median_func
 ! --------------------------------------------------------------------------------------------


 ! --------------------------------------------------------------------------------------------
 subroutine compute_expected_val_logbeta(d_a_par, c_par, vc_lower, vc_upper, rho, e_val_vc, e_val_vc2,&
                                         expected_value_logbeta)
 ! --------------------------------------------------------------------------------------------
 !PURPOSE: compute the expected value of vc between an interval vc_lower and vc_upper
 ! NOTE1:   k11 is prevented from becoming 0 to avoid divisions by zero
 ! NOTE2:   c_par is prevented from being 0 to avoid divisions by zero
 ! The approximations made should have negligible effects on the calculations. unless an infinite range for vc is used
 ! while the curve is a perfectly separated degenerate curve.

 implicit none

 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: vc_lower, vc_upper ! cutoff values between which the integration has to be
                                  ! performed to compute the expected value
 real(kind=double), intent(IN):: rho ! prevalence (fraction of positive cases)
 real(kind=double), intent(OUT):: e_val_vc, e_val_vc2 ! Expected value of vc and vc squared
 real(kind=double), intent(OUT):: expected_value_logbeta ! your guess

 real(kind=double):: k11, k21, k22, k13, ky ! These are the constants used in the proper binormal model inside of
               ! the two phi functions that define FPF and TPF.
               ! i.e., FPF(vc) = Phi(k11*vc - k22) + Phi(k11*vc + k13) - H(c)
               !       TPF(vc) = Phi(k21*vc + k22) + Phi(k21*vc + k13) - H(c)
 real(kind=double) :: a, b !  parameters of the conventional binormal model
 real(kind=double):: fpf1,fpf2, tpf1, tpf2 ! nomralization for expected values



! k11 and k21 can be prevented from becoming zero because the function tends asymptotically to zero when they become very
! small because the curve becomes perfectly separated and thus expected values tend to be zero over any non infinite
! range -- This is not true if the range is the full range, but it is not relevant because this is not what this subroutine
! is supposed to compute.
 k11 = - max( ( 1.0_double - c_par ), 1.0e-8_double ) ! Prevent k11 from becoming 0, which would
               ! create problems with division by zero.
 k21 = - max( ( 1.0_double + c_par ), 1.0e-8_double ) ! Prevent k21 from becoming 0, which would
               ! create problems with division by zero.
 k22 = + ( d_a_par / 2.0_double ) * sqrt( 1.0_double + c_par**2 )

 b = (c_par + 1.0_double) / (1.0_double - c_par)

 a =  d_a_par * sqrt ( 1.0_double + b**2) / sqrt(2.0_double)

 ky = 2.0_double * a / (1.0_double + b)


if(rho .speq. 0.0_double) then
! The expected value is computed only using the negative cases
! Compute the expected value using the LLPesce's formula. The second part vanishes with c ->0, but the numerical approximation
! can become problematic, so we eliminate it (divisions by zero)

         ! determine normalization for expected value (probability over integration interval)
         call fpf_PBM(d_a_par, c_par, vc_lower, fpf1)
         call fpf_PBM(d_a_par, c_par, vc_upper, fpf2)



         if ( abs(c_par) > 1.0e-8_double)  then
             k13 = k22 / c_par
             e_val_vc  =  Integral1(k11, -k22, vc_lower, vc_upper) + Integral1(k11,  k13, vc_lower, vc_upper)
             e_val_vc2 =  Integral2(k11, -k22, vc_lower, vc_upper) + Integral2(k11,  k13, vc_lower, vc_upper)
         else ! if c == 0 the second order term in the relationship between vc and log beta vanishes with its coefficient
             e_val_vc  =  Integral1(k11, -k22, vc_lower, vc_upper)
             e_val_vc2 =  Integral2(k11, -k22, vc_lower, vc_upper)
         endif
        ! normalize to take the expected value over the range vc1 vc2
         e_val_vc  = e_val_vc / (fpf1 - fpf2)
         e_val_vc2  = e_val_vc2 / (fpf1 - fpf2)

elseif(rho .speq. 1.0_double) then
 ! The expected value is computed only using the positive cases
 ! Compute the expected value using the LLPesce's formula. The second part vanishes with c ->0, but the numerical approximation
 ! does not work.
         ! determine normalization for expected value (probability over integration interval)
         call tpf_PBM(d_a_par, c_par, vc_lower, tpf1)
         call tpf_PBM(d_a_par, c_par, vc_upper, tpf2)

         if ( abs(c_par) > 1.0e-8_double )  then
             k13 = k22 / c_par
             e_val_vc  = Integral1(k21, k22, vc_lower, vc_upper) + Integral1(k21,  k13, vc_lower, vc_upper)
             e_val_vc2 = Integral2(k21, k22, vc_lower, vc_upper) + Integral2(k21,  k13, vc_lower, vc_upper)
         else
             e_val_vc  = Integral1(k21, k22, vc_lower, vc_upper)
             e_val_vc2 = Integral2(k21, k22, vc_lower, vc_upper)
         endif
         ! normalize to take the expected value over the range vc1 vc2
         e_val_vc  = e_val_vc / (tpf1 - tpf2)
         e_val_vc2  = e_val_vc2 / (tpf1 - tpf2)
else
 ! The expected value is computed using a mix of positive and negative cases
 ! Compute the expected value using the LLPesce's formula. The second part vanishes with c ->0, but the numerical approximation
 ! does not work.
        ! determine normalization for expected value (probability over integration interval)
        ! positive and negative normalize independently because the expected values are taken independently)
        call fpf_PBM(d_a_par, c_par, vc_lower, fpf1)
        call fpf_PBM(d_a_par, c_par, vc_upper, fpf2)
        call tpf_PBM(d_a_par, c_par, vc_lower, tpf1)
        call tpf_PBM(d_a_par, c_par, vc_upper, tpf2)

        if ( abs(c_par) > 1.0e-8_double ) then
             k13 = k22 / c_par
                     e_val_vc2 = (1.0_double-rho) * ( &
                                   Integral2(k11,-k22,vc_lower,vc_upper) + Integral2(k11,k13,vc_lower,vc_upper) &
                                 )+&
                                 rho * ( &
                                  Integral2(k21, k22, vc_lower, vc_upper) + Integral2(k21,  k13, vc_lower, vc_upper) &
                                 )
                     e_val_vc  = (1.0_double-rho) * ( &
                                   Integral1(k11,-k22,vc_lower,vc_upper) + Integral1(k11,k13,vc_lower,vc_upper) &
                                 )+&
                                 rho * ( &
                                  Integral1(k21, k22, vc_lower, vc_upper) + Integral1(k21,  k13, vc_lower, vc_upper) &
                                 )
        else
                     e_val_vc2 = (1.0_double-rho) * ( Integral2(k11,-k22,vc_lower,vc_upper) )+ &
                                 rho              * ( Integral2(k21, k22, vc_lower, vc_upper) )
                     e_val_vc = (1.0_double-rho) * Integral1(k11,-k22, vc_lower, vc_upper) + &
                                 rho             * Integral1(k21, k22, vc_lower, vc_upper)
        endif
         ! normalize to take the expected value over the range vc1 vc2
        e_val_vc2 = e_val_vc2 / ( rho*(tpf1-tpf2) + (1.0_double - rho)*(fpf1-fpf2))
        e_val_vc = e_val_vc / ( rho*(tpf1-tpf2) + (1.0_double - rho)*(fpf1-fpf2))

endif

 ! transform from vc to log beta
 expected_value_logbeta =  ky * e_val_vc  - 2.0_double * c_par * e_val_vc2   + log(b)


 ! --------------------------------------------------------------------------------------------
 end subroutine compute_expected_val_logbeta
 ! --------------------------------------------------------------------------------------------


 ! --------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------
    real(kind=double) function  Integral2(alpha1, alpha2, vc1, vc2)
 ! --------------------------------------------------------------------------------------------
 ! PURPOSE: compute the general integral type that is the basic ingredient to compute expected values
 !          of v^2. Refer to the LLPesce formula
 ! NOTE:    it deals roughly with preventing divisions by zero, the calling program should make sure that
 !          alpha1 is not becoming zero.
 use statistic_functions, only : phi

 implicit none

 real(kind=double), intent(IN):: alpha1, alpha2 ! constants used to compute the integral
 real(kind=double), intent(IN):: vc1, vc2 ! integration limits
 real(kind=double):: gamma1, gamma2 ! arguments of special functions

 gamma1 = alpha1*vc1  + alpha2
 gamma2 = alpha1*vc2  + alpha2


 integral2 =   &
              (alpha2**2 + 1.0_double) * ( phi(gamma1)  -  phi(gamma2) ) + &
              (1.0_double/sqrt(2.0_double*pi))* ( &
                            exp(- gamma1**2 / 2.0_double) * (2.0_double*alpha2 - gamma1) - &
                            exp(- gamma2**2 / 2.0_double) * (2.0_double*alpha2 - gamma2)   &
              )


 integral2 = integral2 / max( alpha1**2, 1.0e-20_double ) ! Avoid divisions by zero


 ! --------------------------------------------------------------------------------------------
   end function Integral2
 ! --------------------------------------------------------------------------------------------


 ! --------------------------------------------------------------------------------------------
 subroutine compute_expected_val_vc(d_a_par, c_par, vc_lower, vc_upper, rho, expected_value_vc)
 ! --------------------------------------------------------------------------------------------
 !PURPOSE: compute the expected value of vc between an interval vc_lower and vc_upper
 ! NOTE1:   k11 is prevented from becoming 0 to avoid divisions by zero
 ! NOTE1:   c_par is prevented from being 0 to avoid divisions by zero
 ! The approximations made should have negligible effects on the calculations. unless an infinite range for vc is used
 ! while the curve is a perfectly separated degenerate curve.

 implicit none

 real(kind=double), intent(IN):: d_a_par, c_par ! curve parameters
 real(kind=double), intent(IN):: vc_lower, vc_upper ! cutoff values between which the integration has to be
                                  ! performed to compute the expected value
 real(kind=double), intent(IN):: rho ! prevalence (fraction of positive cases)
 real(kind=double), intent(OUT):: expected_value_vc ! your guess

 real(kind=double):: k11, k21, k22, k13 ! These are the constants used in the proper binormal model inside of
               ! the two phi functions that define FPF and TPF.
               ! i.e., FPF(vc) = Phi(k11*vc - k22) + Phi(k11*vc + k13) - H(c)
               !       TPF(vc) = Phi(k21*vc + k22) + Phi(k21*vc + k13) - H(c)
 real(kind=double):: fpf1,fpf2, tpf1, tpf2 ! nomralization for expected values

! k11 and k21 can be prevented from becoming zero because the function tends asymptotically to zero when they become very
! small because the curve becomes perfectly separated and thus expected values tend to be zero over any non infinite
! range -- This is not true if the range is the full range, but it is not relevant because this is not what this subroutine
! is supposed to compute.
 k11 = - max( ( 1.0_double - c_par ), 1.0e-8_double ) ! Prevent k11 from becoming 0, which would
               ! create problems with division by zero.
 k21 = - max( ( 1.0_double + c_par ), 1.0e-8_double ) ! Prevent k21 from becoming 0, which would
               ! create problems with division by zero.
 k22 = + ( d_a_par / 2.0_double ) * sqrt( 1.0_double + c_par**2 )


 if(rho .speq. 0.0_double) then
  ! The expected value is computed only using the negative cases
 ! Compute the expected value using the LLPesce's formula. The second part vanishes with c ->0, but the numerical approximation
 ! does not work.
         ! determine normalization for expected value (probability over integration interval)
         call fpf_PBM(d_a_par, c_par, vc_lower, fpf1)
         call fpf_PBM(d_a_par, c_par, vc_upper, fpf2)

         if ( abs(c_par) > 1.0e-8_double)  then
             k13 = k22 / c_par
             expected_value_vc = (Integral1(k11, -k22, vc_lower, vc_upper)+Integral1(k11,  k13, vc_lower, vc_upper))/(fpf1-fpf2)
         else
             expected_value_vc = (Integral1(k11, -k22, vc_lower, vc_upper)  ) /(fpf1 - fpf2)
         endif

 elseif(rho .speq. 1.0_double) then
 ! The expected value is computed only using the positive cases
 ! Compute the expected value using the LLPesce's formula. The second part vanishes with c ->0, but the numerical approximation
 ! does not work.
         ! determine normalization for expected value (probability over integration interval)
         call tpf_PBM(d_a_par, c_par, vc_lower, tpf1)
         call tpf_PBM(d_a_par, c_par, vc_upper, tpf2)

         if ( abs(c_par) > 1.0e-8_double )  then
             k13 = k22 / c_par
             expected_value_vc = (Integral1(k21, +k22,vc_lower,vc_upper)+Integral1(k21, k13,vc_lower,vc_upper))/(tpf1-tpf2)
         else
             expected_value_vc = (Integral1(k21, +k22, vc_lower, vc_upper)) / (tpf1 - tpf2)
         endif

 else
 ! The expected value is computed using a mix of positive and negative cases
 ! Compute the expected value using the LLPesce's formula. The second part vanishes with c ->0, but the numerical approximation
 ! does not work.
        ! determine normalization for expected value (probability over integration interval)
        ! positive and negative normalize independently because the expected values are taken independently)
        call fpf_PBM(d_a_par, c_par, vc_lower, fpf1)
        call fpf_PBM(d_a_par, c_par, vc_upper, fpf2)
        call tpf_PBM(d_a_par, c_par, vc_lower, tpf1)
        call tpf_PBM(d_a_par, c_par, vc_upper, tpf2)

         if ( abs(c_par) > 1.0e-8_double ) then
             k13 = k22 / c_par
             expected_value_vc = &
                (1.0_double - rho)*(Integral1(k11,-k22,vc_lower,vc_upper)+Integral1(k11,k13,vc_lower,vc_upper) )/(fpf1-fpf2)  + &
                rho               *(Integral1(k21,+k22,vc_lower,vc_upper)+Integral1(k21,k13,vc_lower,vc_upper) )/(tpf1-tpf2)
         else
             expected_value_vc = &
                (1.0_double - rho) * Integral1(k11, -k22, vc_lower, vc_upper)/(fpf1-fpf2)  + &
                rho                * Integral1(k21, +k22, vc_lower, vc_upper)/(tpf1-tpf2)
         endif
 endif


 ! --------------------------------------------------------------------------------------------
 end subroutine compute_expected_val_vc
 ! --------------------------------------------------------------------------------------------


 ! --------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------
    real(kind=double) pure function  Integral1(alpha1, alpha2, vc1, vc2)
 ! --------------------------------------------------------------------------------------------
 ! PURPOSE: compute the general integral type that is the basic ingredient to compute expected values
 !          of cutoffs. Refer to the LLPesce formula
 ! NOTE:    it deals roughly with preventing divisions by zero, the calling program should make sure that
 !          alpha1 is not becoming zero.
 use statistic_functions, only : phi

 implicit none

 real(kind=double), intent(IN):: alpha1, alpha2 ! constants used to compute the integral
 real(kind=double), intent(IN):: vc1, vc2 ! integration limits


 integral1 = (  &
                exp(- (alpha1 * vc1 + alpha2)**2 / 2.0_double) -  &
                exp(- (alpha1 * vc2 + alpha2)**2 / 2.0_double)    &
              ) /sqrt(2.0_double*pi) + &
              alpha2 * ( phi(alpha1*vc1 + alpha2)  -  phi(alpha1*vc2 + alpha2) )

 integral1 = integral1 / max( abs(alpha1), 1.0e-10_double ) ! Avoid divisions by zero


 ! --------------------------------------------------------------------------------------------
   end function Integral1
 ! --------------------------------------------------------------------------------------------

 ! --------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------
  subroutine points_at_cutoffs_PBM(d_a_par_in, c_par_in, cutoffs, num_cutoffs, CurvePoints, ierror)
 ! ----------------------------------------------------------------------
 ! PURPOSE:  procedure returns a sequence of fpf and tpf for a
 !           prespecified (input) sequence of cutoffs (category
 !           boundaries, sometimes called betas)
 !ALGORITHM: Verify whether the cutoffs and parameters have acceptable values,
 !           then simply compute the FPF, TPF values
 !

 implicit none

 real(kind=double), intent(IN):: d_a_par_in, c_par_in ! curve parameters
 integer, intent(IN) :: num_cutoffs ! number of cutoffs at which (FPF, TPF)  values are  desired
 real(kind=double), dimension(num_cutoffs), intent(in) :: cutoffs ! the actual cutoffs values
 real(kind=double), dimension(2, num_cutoffs + 1), intent(OUT) :: CurvePoints ! the actual
                  !  array with the fpf, tpf values
 integer, intent(OUT):: ierror  ! fit_OK -> OK; fit_failed -> Failed; bad_input -> wrong input

 real(kind=double):: fpf, tpf, junk
 real(kind=double):: d_a_par, c_par ! curve parameters
 real(kind=double):: v_c, v_c_max, v_c_min ! Maximum and minimum values possible for the cutoffs. Some models
                    ! like the proper binormal model, have a bounded cutoff space

 integer:: i

 ! Take care of values which produce numerically unstable data points. In this case
 ! we are producing points for the purpose of plotting them, so we certainly don't
 ! need an accuracy of more than six digits (which is numerically unattainable anyway
 ! because of the double precision arithmetic and the numerical precision of the
 ! various functions used)

 ierror = fit_OK ! initialize to OK

 ! CHECK the value of c_par.  c_par == 1 creates numerical problems because the equations lose
 ! meaning.
 if( abs(c_par_in) >  1.0_double) then ! unacceptable value
   ierror = bad_input
   return
 elseif( abs(c_par_in) .speq. 1.0_double) then
   c_par = sign(.999999_double, c_par)
 else
   c_par = c_par_in
 endif

 ! Take care of d_a_par. d_a == 0
 if( d_a_par_in <  0.0_double) then
   ierror = bad_input
   return
 elseif( d_a_par_in .speq. 0.0_double) then
   d_a_par = .0000001_double
 else
   d_a_par = d_a_par_in
 endif

! Verify that the cutoffs are compatible with the d_a, c values

 ! Find max and min values for cutoffs for specific parameter values
 if (ABS(c_par) < c_almost_zero ) then ! If c is close to 0, any value is acceptable
       v_c_min = -huge(1.0_double)
       v_c_max = huge(1.0_double)
 elseif ( c_par <= - c_almost_zero ) then ! if c is negative, v_c is bounded from below
       v_c_min =  d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * c_par)
       v_c_max = huge(1.0_double)
 else! if c is positive, v_c is bounded from above
       v_c_min = -huge(1.0_double)
       v_c_max =  d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * c_par)
 endif

 ! check whether any of the cutoffs is outside the valid domain
 if( any(cutoffs < v_c_min)  .or.  any(cutoffs > v_c_max)   ) then
   ierror = bad_input
   return
 endif


! Compute the fpf, tpf values for the cutoffs
! Note that the first point is the 1.0 1.0 point
  CurvePoints(1,1) = 1.0_double
  CurvePoints(2,1) = 1.0_double
  do i = 1, num_cutoffs
        v_c = cutoffs(i)
        call fpf_PBM(d_a_par, c_par, v_c, fpf, junk)
        call tpf_PBM(d_a_par, c_par, v_c, tpf, junk)
        CurvePoints(1,i+1) = fpf
        CurvePoints(2,i+1) = tpf
  enddo



 ! ----------------------------------------------------------------------
  end subroutine points_at_cutoffs_PBM
 ! ----------------------------------------------------------------------



 ! ----------------------------------------------------------------------
 ! ----------------------------------------------------------------------
 subroutine points_on_curve_PBM(d_a_par_in, c_par_in, num_pts, CurvePoints,ierror)
 ! ----------------------------------------------------------------------
 !PURPOSE: return a sequence of fpf,tpf values for the proper ROC curve with
 !         parameters  d_a_par, c_par.
 !ALGORITHM: the points will be put at the intersection of a sequence of equally
 !          spaced radii with origin in the point (1,0) of the ROC plot. This choice
 !          has been made to make plotting of skewed curve smooth.
 !NOTES     : modified by LP on September 2008 to include 0,0 and 1,1 exactly.

 implicit none

 real(kind=double), intent(IN):: d_a_par_in, c_par_in ! curve parameters
 integer, intent(IN) :: num_pts ! number of curve points whose value is desired
 real(kind=double), dimension(2, num_pts), intent(OUT) :: CurvePoints ! the actual
                  !  array with the fpf, tpf values
 integer, intent(OUT):: ierror  ! fit_OK -> OK; fit_failed -> Failed; bad_input -> wrong input


 real(kind = double) :: pi2   = pi/2.0_double
 real(kind=double):: v_c_min, v_c_max  ! grid values where to sample the cutoff
                         ! parameter to compute the (FPF, TPF) array
 real(kind=double):: upper_v_c, lower_v_c, v_c ! current value of v_c
 real(kind=double):: target_angle, current_angle ! used to determine the radii
 real(kind=double):: fpf, tpf, junk, temp
 real(kind=double):: d_a_par, c_par ! curve parameters

 integer, parameter :: max_iter = 1000
 integer :: i, j_iter ! index for the points to be plotted, and iteration counter


 ! Take care of values which produce numerically unstable data points. In this case
 ! we are producing points for the purpose of plotting them, so we certainly don't
 ! need an accuracy of more than six digits (which is numerically unattainable anyway
 ! because of the double precision arithmetic and the numerical precision of the
 ! various functions used)

 ierror = fit_OK ! initialize error message to OK.

 ! CHECK the value of c_par.  c_par == 1 creates numerical problems because the equations lose
 ! meaning.
 if( abs(c_par_in) >  1.0_double) then ! unacceptable value
   ierror = bad_input
   return
 elseif( abs(c_par_in) .speq. 1.0_double) then
   c_par = sign(.999999_double, c_par_in)
 else
   c_par = c_par_in
 endif

 ! Take care of d_a_par. d_a == 0
 if( d_a_par_in <  0.0_double) then
   ierror = bad_input
   return
 elseif( d_a_par_in .speq. 0.0_double) then
   d_a_par = .0000001_double
 else
   d_a_par = d_a_par_in
 endif


! Determine the minimum value of v_c that it is worthy using. Here we
! look for a value for which FPF =~ 1. Since TPF > FPF v v_c it follows
! that all v_c's smaller than this one will map to the same point on the curve
 if (ABS(c_par) < c_almost_zero ) then
        call find_v_c_min(d_a_par, c_par, v_c_min, ierror)
        if (ierror == 0) call find_v_c_max(d_a_par, c_par, v_c_max, ierror)
 elseif ( c_par <= - c_almost_zero ) then
        v_c_min =  d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * sign( max(.00000000001_double, abs(c_par)),c_par))
        call find_v_c_max(d_a_par, c_par, v_c_max, ierror)
 else
        call find_v_c_min(d_a_par, c_par, v_c_min, ierror)
        v_c_max =  d_a_par * SQRT ( 1.0_double + c_par**2) / ( 4.0_double * sign( max(.00000000001_double, abs(c_par)),c_par))
 endif

 if(ierror /= 0) then
     ierror = fit_failed
     return ! failed to bracket solutions, return error message
 endif

! Verify that the lower boundary for v_c is not too negative, so that the value of
! tpf and fpf doesn't change while the cutoff is halved
 temp = v_c_min / 2.0_double
 do i = 1, max_iter
          call fpf_PBM(d_a_par, c_par, temp, fpf, junk)
          call tpf_PBM(d_a_par, c_par, temp, tpf, junk)
          if( tpf > .999999999999999_double .and.  fpf > .999999999999999_double ) then
             v_c_min = temp
             temp =  v_c_min / 2.0_double
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
 temp = v_c_max / 2.0_double
 do i = 1, max_iter
          call fpf_PBM(d_a_par, c_par, temp, fpf, junk)
          call tpf_PBM(d_a_par, c_par, temp, tpf, junk)
          if( fpf < .000000000000001_double .and.  tpf < .000000000000001_double) then
             v_c_max = temp
             temp =  v_c_max / 2.0_double
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

 v_c = v_c_max

 target_angle = 0.0_double  ! start sampling angles in ROC space starting from 0,
                            ! which is a point on the TPF = 0 axis


 plotpoints: do i = 1, num_pts
    ! Set value to bracket the cutoff value for this point we are seeking
    upper_v_c = v_c
    lower_v_c = v_c_min

     ! Start a bisection iteration to find the cutoff value, and thus the FPF/TPF value
     ! that intersects that specific radius from (1,0)
     find_cut_off: do j_iter = 1, max_iter
          v_c = (upper_v_c + lower_v_c ) / 2.0_double
          ! Find fpf and tpf at the new bisected point
          call fpf_PBM(d_a_par, c_par, v_c, fpf, junk)
          call tpf_PBM(d_a_par, c_par, v_c, tpf, junk)

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
                lower_v_c = v_c
                v_c = ( v_c + upper_v_c)/2.0_double
           else
                upper_v_c = v_c
                v_c = ( v_c + lower_v_c)/2.0_double
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


 end subroutine points_on_curve_PBM

!-------------------------------------------------------
  subroutine  find_v_c_max(d_a_par, c_par, v_c_max, ierror)
! PURPOSE: find a cutoffs value that produces a value of TPF very close
!          to 0. Since TPF > FPF, it follows that also FPF = 0.
! ALGORITHM: Stewise increase (superlinear increase)

  implicit none

  real(kind=double), intent(in):: d_a_par, c_par ! curve parameters
  real(kind=double), intent(out):: v_c_max ! the actual minimum
  integer, intent(out):: ierror ! error flag

  real(kind=double):: tpf, one_minus_tpf
  real(kind=double):: step
  integer, parameter:: num_iter = 10000 ! it never takes too long, but just in case
  integer:: i ! loop counter

 ! size of change in the fpf function as a function of v_c
 step =  1.0_double / ( max(1.0_double + c_par, 1.0e-4_double) )

 v_c_max = 0.0_double ! the maximum is always >= 0

  do i =  1, num_iter
       call tpf_PBM(d_a_par, c_par, v_c_max, tpf, one_minus_tpf)
       if( tpf <  1.0e-8_double)  then
            exit
       else
           v_c_max = v_c_max + i * step ! make change more than linear
       endif
  enddo

  if ( i > num_iter) then
     ierror = 1
  else
    ierror = 0
  endif

!-------------------------------------------------------
 end subroutine find_v_c_max
!-------------------------------------------------------


!-------------------------------------------------------
  subroutine  find_v_c_min(d_a_par, c_par, v_c_min, ierror)
! PURPOSE: find a cutoffs value that produces a value of FPF very close
!          to 1. Since TPF > FPF, it follows that also TPF = 1.
! ALGORITHM: Stewise increase (superlinear increase)

  implicit none

  real(kind=double), intent(in):: d_a_par, c_par ! curve parameters
  real(kind=double), intent(out):: v_c_min ! the actual minimum
  integer, intent(out):: ierror ! error flag

  real(kind=double):: fpf, one_minus_fpf
  real(kind=double):: step
  integer, parameter:: num_iter = 10000 ! it never takes too long, but just in case
  integer:: i ! loop counter

 ! size of change in the fpf function as a function of v_c
 step =  1.0_double / ( max(1.0_double - c_par, 1.0e-4_double) )

  v_c_min = 0.0_double ! the minimum is always <= 0

  do i =  1, num_iter
       call fpf_PBM(d_a_par, c_par, v_c_min, fpf, one_minus_fpf)
       if( fpf > 1.0_double - 1.0e-8_double)  then
            exit
       else
           v_c_min = v_c_min - i * step ! make change more than linear
       endif
  enddo

  if ( i > num_iter) then
     ierror = 1
  else
    ierror = 0
  endif

!-------------------------------------------------------
 end subroutine find_v_c_min
!-------------------------------------------------------


 end module proproc_out
