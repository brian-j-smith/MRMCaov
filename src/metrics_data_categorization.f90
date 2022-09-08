! THIS MODULE contains fortran 95 ANSI procedures (subroutines and functions) to
! categorize  continuously distributed test results via the creation of truth-state
! runs of rank ordered data and subsequent reduction in the number of truth-state runs
! via the so called LABROC5 algorithm.

! DEFINITIONS:
!               ACTUALLY POSITIVE: a case, with a specific TEST RESULT VALUES, which is known by some other
!                   RELIABLE analysis to have the disease, or signal that we are looking for.
!                   Also called "abnormal cases" or signal cases
!               ACTUALLY NEGATIVE: a case, with a specific TEST RESULT VALUES, which is known by some other
!                   RELIABLE analysis to show *no* signs of the disease, or signal we are looking for.
!                    Also called "normal cases" or noise cases (all what is measure is noise with regards to this
!                    experiment)
!               CATEGORY: a category is normally defined as the set of cases with some caracteristics. The characteristic
!                    in the case of this categorization module is to have TEST RESULT VALUE between two boundaries, which
!                    include at least one TRUTH RUN, and all of the included TRUTH RUNS are included completely.
!                    When ratings are used the category is simply the set of cases which are classified with that rating.
!               CATEGORY BOUNDARY: the TEST RESULT VALUE that separates the TEST RESULTS VALUES of one category to the ones
!                    of another one. Typically we use the mean of the smallest test value results of one and the largest
!                    of the other one (chosen so that they are contiguous)
!               CORNER: when the raw data is plotted on an ROC plot, the LBAROC4 points (or CATEGORY BOUNDARIES between
!                        TRUTH RUNS) look like corners. For this reason they are sometimes calles LABROC4 corners
!               FPF or FALSE Positive Fraction: the fraction of actually negative cases (or normal cases) that have
!                   TEST RESULT VALUE larger than a threshold
!               OPERATING POINT: an operating point is a point on an ROC curve, and it is normally defined with it
!                   coordinates as (FPF,TPF). It corresponds to a Category Boundary in terms of TEST RESULT VALUE, since
!                   the FPF and TPF of the OPERATING POINT are the fraction of cases which have TEST RESULT VALUE larger
!                   than the last category boundary included in the calculation of the FPF & TPF.
!               TPF or TRUE Positive Fraction: the fraction of actually positive cases (or abnormal cases) that have
!                   TEST RESULT VALUE larger than a threshold. Typically considered with the corresponding FPF
!               TEST RESULT VALUE: The test result value is simply the measurement that is made, it could be a
!                   glucose concentration level or a category for a radiography. Every case will have a
!                   TEST RESULT VALUE otherwise it cannot be used to asses the methodology we are trying to
!                   asses
!              TRUTH STATE: if the case ware found to be actually positive (or abnormal, or with "signal" present)
!                    or negative (or "normal", or only with "noise" present). One needs some sort of truth to able
!                    to tell if a method is working or not.
!              TRUTH RUN OR TRUTH STATE RUN: A truth state run is a contiguous sequence of cases either of the
!                   same truth state (actualy positive or actually negative) or of different truth states but
!                   with the same test result value. Their are the atoms of a categorization algorithm for an
!                   ROC curve, because separation between truth states, and thus True Positive Fraction (TPF) vs
!                   False Positive Fraction tradeoff analysis, cannot be done within these truth state runs. I.e,
!                   we cannot set a threshold which would allow us in the detection to tell any of this states from
!                   the other using the test result value


! NOTE: The algorithm here implemented is slightly different from the one described in the
!       reference. Where the 2 algorithms differ there is a comment pointing it out.
! NOTE: If one desires to modify the categorization scheme, we recommend to reinstate the 2 consistency checks
!       which appear in the module and that are commented right now. Search for Write(* or stop.
! REFERENCE: Charles E. Metz, Benjamin A. Herman and Jong-Her Shen
!            "MAXIMUM LIKELIHOOD ESTIMATION OF RECEIVER OPERATING CHARACTERISTIC (ROC)
!             CURVES FROM CONTINUOUSLY-DISTRIBUTED DATA"
!                   - STATISTICS IN MEDICINE 17, 1033-1053 (1998)
!            NOTE: THIS PAPER WILL REFERRED FROM NOW ON AS the Continuous ROC paper

! PROGRAMMING NOTES: 1) The code is ANSI Fortran 95, it will not compile unders strict Fortran 90 or previous
!                        versions of the fortran language
!                    2) Deprecated forms of fortran 77 are not used
!                    3) The module has the "private" status, so only the methods which are declared as public can
!                       be used by external programs
!                    4) In ROC, the ROC curve values (FPF and  TPF) increase when the operating points TEST RESULTS VALUES are
!                       decresed. We decided to keep the indexing so that TEST RESULT VALUES are in ascending order, which
!                       makes FPF and TPF and related quantities decrease with it.

! WARNING:          1) the array which contains the final categories is dimensioned from above using a the integer
!                      MAX_NUM_CATEGORIES. This is done for a number of reasons, among which are:
!                      A)  because we do not know how many categories we are going to
!                          have until we analyze the data and it is complicated to have a subroutine dimension an array
!                          for the main program (We will deal with it at some point in the future)
!                      B)  The speed of the calculation scales at least with the cube of the number of categories, while the
!                          accuracy of the calculation does not change too much over 20 categories, so it can be a good idea
!                          to collapse them to 20.
!                     WE RECOMMEND TO CALL CATGRZ USING about 400 AS MAX_NUM_CATEGORIES (but it works also for more or less)
!                       AS FAR AS WE KNOW THE ERROR IN THE ROC AREA INDUCED BY COLLAPSING THE CATEGORIES TO about 20
!                       CAN BE A LARGE AS 0.01. SO IF YOU HAVE LARGE DATASETS WITH SMALL SE BE CAREFUL.
! CHANGES:
! 1/11/06 LP : Added the optional array CASE_CAT in the call to catgrz that contains the category in which each case is allocated.
!              To be used by MRMC.
!              Added val_index to compute CASE_CAT
! 3/21/06 LP: Added public procedure  LABROC4_Collapser. Subroutine analyzes a set of categorical data (array with number of
!             actually positive cases per category  and array with number of actually negative cases per category) to find
!             the labroc4 corners.
MODULE CATEGORIZATION


! NOTE that only the constants for actually negative and actually positive cases are loaded from the module
! array_dimensions. The only is used to that a strick control of data and methods is enforced

USE data_types ! where data types are defined
USE array_dimensions, ONLY: ACT_NEG, ACT_POS ! where array related quantities are
USE l_algebra, only: qsortd

! GENERAL SETTING FOR THE WHOLE MODULE
IMPLICIT NONE

PRIVATE ! All the variables and functions that are not set as public are not acessible
PUBLIC CATGRZ ! only the highest soubroutine is callable from outside
PUBLIC LABROC4_Collapser
public countm1m12m2 ! determine, from a dataset of mn negative cases and ms positive cases with missing data how many
                    ! cases have measurements only for modality one, how many only for modality two and how many for both
                    ! it is a key component of ROCKIT and other partially-paired ROC curve estimation procedures.
public find_truthruns_d ! determine how many truthruns are in a dataset where some of the values have to be ignored, as described
                    ! by the design matrix
public categorize_mod1_mod2 ! categorize partially-paired two modality data into data for mod1, mod2, and mod12

! Global variables within the module
 INTEGER, PARAMETER :: TST_RES_VAL = - 1 !index to locate the test value result in multidimensional arrays like
                        ! TRUTH_RUNS_CATEGORY. It is global in this module because these arrays are used
                        ! everywhere and it is a parameter, so there is no risk to mess it up. It is private to
                        ! this module because all of these arrays are internal to the module.

! Subroutines and functions code follows
CONTAINS

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine categorize_mod1_mod2 (mn, ms, data_neg, data_pos, des_neg, des_pos, num_cat1, num_cat2, &
                                 catn1, catn12, catn2, cats1, cats12, cats2, ierror)
!---------------------------------------------------------------------------
! PURPOSE: convert input of partially-paired two modality data into categorical data split into unpaired and
!          fully paired data
! NOTE:    the number of categories is taken externally, both to have no issues with dimensioning of arrays and
!          to leave the decisions about the maximum number of truth runs accepted to the calling program. This
!          number has to be meaninfgful and smaller than or equal to the number of categories actually found in
!          the data.

integer, intent(IN):: mn ! number of negative cases
integer, intent(IN):: ms ! number of positive cases

real(kind=double), dimension(2, mn), intent(IN):: data_neg ! actually-negative input data for the two modalities to be analyzed
real(kind=double), dimension(2, ms), intent(IN):: data_pos ! actually-positive input data for the two modalities to be analyzed

! Design matrices. Here we assume that if there are values different from 0 or 1, there is an input error (in general there are
! algorithms that allow the use of different flags for the design matrix, for example to indicate clustering, however, ROCKIT
! cannot make use of them and therefore will not accept them.
integer, dimension(2, mn), intent(IN):: des_neg ! actually-negative design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed
integer, dimension(2, ms), intent(IN):: des_pos ! actually-positive design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed

integer, intent(IN):: num_cat1, num_cat2 ! number of categories in modality one and two.

integer, dimension(num_cat1), intent(OUT):: catn1, cats1 ! categorical data by truth, independent part for mod1
integer, dimension(num_cat2), intent(OUT):: catn2, cats2 ! categorical data by truth, independent part for mod2
integer, dimension(num_cat1,num_cat2), intent(OUT) :: catn12, cats12   ! categorical data by modality and truth, fully-paired part

integer, intent(OUT):: ierror ! error flag: 0 -> fit worked properly ; 2 -> wrong design matrix values on input

integer:: mn1, mn12, mn2 ! number of actually-negative cases only in modality 1, in both modalities and only in modality 2
integer:: ms1, ms12, ms2 ! number of actually-positive cases only in modality 1, in both modalities and only in modality 2

real(kind=double), allocatable, dimension(:) :: data_neg1, data_neg2 ! actually-negative for the two modalities, separated
real(kind=double), allocatable, dimension(:) :: data_pos1, data_pos2 ! actually-negative for the two modalities, separated
integer, allocatable, dimension(:,:) :: case_cat1, case_cat2 ! category where each case belongs, from catgrz

integer, dimension(2,mn):: index_neg  ! has the location in the array data_neg1, data_neg2 of the case at location
                                                  ! (modality,i) in data_neg -- some locations might be empty because data is missing
integer, dimension(2,ms):: index_pos  ! has the location in the array data_pos1, data_pos2 of the case at location
                                                  ! (modality,i) in data_pos -- some locations might be empty because data is missing
integer:: index1, index2 ! location of a datavalue in arrays data_neg1 (data_pos1) and data_neg2 (data_pos2)

integer:: num_cat1_loc, num_cat2_loc ! maximum number of categories found in the data
integer, dimension(2,num_cat1):: cat1 ! categorical data for modality 1
integer, dimension(2,num_cat2):: cat2 ! categorical data for modality 2

integer:: i

! Count how many cases are in modality 1, 2 or both
call countm1m12m2(mn, ms, des_neg, des_pos,  mn1, mn12, mn2, ms1, ms12, ms2, ierror)
if(ierror /= 0) return


! Allocate arrays with the data by modality. These arrays are necessary to categorize the data by modality before creating the
! categorical data which includes the information about how the data is paired. However the categorization has to be done
! independently and with all the data because there is no reason to expect that the two scales are identical or even comparable.
allocate(data_neg1(mn1+mn12))
allocate(data_pos1(ms1+ms12))
allocate(data_neg2(mn2+mn12))
allocate(data_pos2(ms2+ms12))


! Start to load the data by modality with the actually-negative data, keep track of the ordering of the data
data_neg1 = 0.0_double
data_neg2 = 0.0_double
index_neg = 0
index1 = 0
index2 = 0

! Load the data for the two different modalities, storing the locations in the new arrays of the input data
do i = 1, mn
     if(des_neg(1,i) == 1 .and. des_neg(2,i) == 1) then
            index1 = index1+1
            data_neg1(index1) = data_neg(1,i)
            index_neg(1,i) = index1
            index2 = index2+1
            data_neg2(index2) = data_neg(2,i)
            index_neg(2,i) = index2
     elseif(des_neg(1,i) == 1 .and. des_neg(2,i) == 0) then
            index1 = index1+1
            data_neg1(index1) = data_neg(1,i)
            index_neg(1,i) = index1
     elseif(des_neg(1,i) == 0 .and. des_neg(2,i) == 1) then
            index2 = index2+1
            data_neg2(index2) = data_neg(2,i)
            index_neg(2,i) = index2
     endif
enddo

! Actually-positive data
data_pos1 = 0.0_double
data_pos2 = 0.0_double
index_pos = 0
index1 = 0
index2 = 0

! Load the data for the two different modalities, storing the locations in the new arrays of the input data
do i = 1, ms
     if(des_pos(1,i) == 1 .and. des_pos(2,i) == 1) then
            index1 = index1+1
            data_pos1(index1) = data_pos(1,i)
            index_pos(1,i) = index1
            index2 = index2+1
            data_pos2(index2) = data_pos(2,i)
            index_pos(2,i) = index2
     elseif(des_pos(1,i) == 1 .and. des_pos(2,i) == 0) then
            index1 = index1+1
            data_pos1(index1) = data_pos(1,i)
            index_pos(1,i) = index1
     elseif(des_pos(1,i) == 0 .and. des_pos(2,i) == 1) then
            index2 = index2+1
            data_pos2(index2) = data_pos(2,i)
            index_pos(2,i) = index2
     endif
enddo

! Categorize data, independently for the two categories because in general there is no relationship between the
! two

! modality 1
allocate(   case_cat1(2,max(ms1+ms12,mn1,mn12))   )
call catgrz(1, mn1+mn12, ms1+ms12, 0, cat1, data_neg1, data_pos1, num_cat1_loc, num_cat1, case_cat1)
if (num_cat1 /= num_cat1_loc) then
    ierror = 2 ! the input value is incorrect
    return
endif
! modality two
allocate(   case_cat2(2,max(ms2+ms12,mn2,mn12))   )
call catgrz(1, mn2+mn12, ms2+ms12, 0, cat2, data_neg2, data_pos2, num_cat2_loc, num_cat2, case_cat2)
if (num_cat2 /= num_cat2_loc) then
    ierror = 2 ! the input value is incorrect
    return
endif

catn1 = 0; cats1 = 0; catn2 = 0; cats2 = 0
catn12 = 0; cats12 = 0

! Load the CATEGORICAL data.
! Note that one index 1 to 2 refers to modality, while the other 1 to 2 refers to truth, for some arrays
! the dimension is truth (case_cat1, case_cat2)  for other is modality (index_neg, index_pos)
do i = 1, mn
     if(des_neg(1,i) == 1 .and. des_neg(2,i) == 1) then
            catn12( case_cat1(1,index_neg(1,i)), case_cat2(1,index_neg(2,i)) ) =  &
                                catn12( case_cat1(1,index_neg(1,i)),case_cat2(1,index_neg(2,i)) ) + 1
     elseif(des_neg(1,i) == 1 .and. des_neg(2,i) == 0) then
            catn1( case_cat1(1,index_neg(1,i)) ) = catn1( case_cat1(1,index_neg(1,i)) )  + 1
     elseif(des_neg(1,i) == 0 .and. des_neg(2,i) == 1) then
            catn2( case_cat2(1,index_neg(2,i)) ) = catn2( case_cat2(1,index_neg(2,i)) )  + 1
     endif
enddo

do i = 1, ms
     if(des_pos(1,i) == 1 .and. des_pos(2,i) == 1) then
            cats12( case_cat1(2,index_pos(1,i)), case_cat2(2,index_pos(2,i)) ) =  &
                                cats12(case_cat1(2,index_pos(1,i)), case_cat2(2,index_pos(2,i)) ) + 1
     elseif(des_pos(1,i) == 1 .and. des_pos(2,i) == 0) then
            cats1( case_cat1(2,index_pos(1,i)) ) = cats1(case_cat1(2,index_pos(1,i)) )  + 1
     elseif(des_pos(1,i) == 0 .and. des_pos(2,i) == 1) then
            cats2( case_cat2(2,index_pos(2,i)) ) = cats2( case_cat2(2,index_pos(2,i)) )  + 1
     endif
enddo




!---------------------------------------------------------------------------
end subroutine categorize_mod1_mod2
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------



!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine countm1m12m2(mn, ms, des_neg, des_pos, mn1, mn12, mn2, ms1, ms12, ms2, ierror)
!----------------------------------------------------------------------------------------
! PURPOSE: compute the number of positive and negative cases that are available for modality 1 only, modality 2 only or
!          both modalities at the same time. A lot of algorithms that deal with missing data for partially-paired modality
!          comparisons are based on using these three groups of cases. This is often used to dimension arrays or normalization
!          factors.
integer, intent(IN):: mn ! number of negative cases
integer, intent(IN):: ms ! number of positive cases

! Design matrices. Here we assume that if there are values different from 0 or 1, there is an input error (in general there are
! algorithms that allow the use of different flags for the design matrix, for example to indicate clustering, however, ROCKIT
! cannot make use of them and therefore will not accept them.
integer, dimension(2, Mn), intent(IN):: des_neg ! actually-negative design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed
integer, dimension(2, Ms), intent(IN):: des_pos ! actually-positive design matrix (whether a case is present (1) or absent (0)
                                                ! for each the two modalities to be analyzed
integer, intent(OUT):: mn1, mn12, mn2 ! number of actually-negative cases only in modality 1, in both modalities and only in modality 2
integer, intent(OUT):: ms1, ms12, ms2 ! number of actually-positive cases only in modality 1, in both modalities and only in modality 2

integer, intent(OUT):: ierror ! error flag: 0 -> fit worked properly ; 2 -> wrong design matrix values on input

integer:: i

! Count the number of cases in groups 1, 12 and 2 (see in "ALGORITHM") and
! check arrays to see whether cases are missing (it is not accepted if cases are missing from both modalities
! or if values are different from 0 or 1). We decided to opt for a straightforward implementation as opposed to some fancier
! array manipulation, e.g., any(MASK..) or all, to make future maintenance by non experts of fortran 90 simpler.

! initialize counts and count actually-negatives
mn1  = 0 ; mn12 = 0 ; mn2  = 0
do i = 1, mn
     if(des_neg(1,i) == 1 .and. des_neg(2,i) == 1) then
              mn12 = mn12 + 1
     elseif(des_neg(1,i) == 1 .and. des_neg(2,i) == 0) then
              mn1 = mn1 + 1
     elseif(des_neg(1,i) == 0 .and. des_neg(2,i) == 1) then
              mn2 = mn2 + 1
     else
     ! if some values in the array des_neg is unacceptable (diff from 0 or 1) or
     ! ifthere are cases for which cases are missing for both modalities. If is
     ! the case it is safer to assume that there is an input error as opposed to a deliberate
     ! design
            ierror = 2
            return
     endif
enddo

! initialize counts and count actually-positives
ms1  = 0 ; ms12 = 0 ; mn2  = 0
do i = 1, ms
     if(des_pos(1,i) == 1 .and. des_pos(2,i) == 1) then
              ms12 = ms12 + 1
     elseif(des_pos(1,i) == 1 .and. des_pos(2,i) == 0) then
              ms1 = ms1 + 1
     elseif(des_pos(1,i) == 0 .and. des_pos(2,i) == 1) then
              ms2 = ms2 + 1
     else
     ! if some values in the array des_neg is unacceptable (diff from 0 or 1) or
     ! ifthere are cases for which cases are missing for both modalities. If is
     ! the case it is safer to assume that there is an input error as opposed to a deliberate
     ! design
            ierror = 2
            return
     endif
enddo




!----------------------------------------------------------------------------------------
end subroutine countm1m12m2
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine find_truthruns_d(mn_tot, ms_tot, data_neg_tot, data_pos_tot, des_neg, des_pos, num_truth_runs,&
                            max_num_truth_runs, ierror)
!----------------------------------------------------------------------------------------
! PURPOSE: Find truth-runs from a dataset with a design matrix, i.e., some of the input
!          data array might no contain information relevant to this fit. The design matrix
!          contains the cases that should be used to determine the number of categories for
!          this problem
! NOTE:    it returns the total number of truth runs if it is smaller than max_num_truth_runs
! NOTE:    it assumes that positivity is for larger values.

integer, intent(IN):: mn_tot ! number of negative cases
integer, intent(IN):: ms_tot ! number of positive cases

real(kind=double), dimension(mn_tot), intent(IN):: data_neg_tot ! actually-negative input data, including cases not used here
real(kind=double), dimension(ms_tot), intent(IN):: data_pos_tot ! actually-positive input data  including cases not used here

! Design matrices. Here we assume that if there are values different from 0 or 1, there is an input error
integer, dimension(mn_tot), intent(IN):: des_neg ! actually-negative design matrix (whether a case is present (1) or absent (0)
integer, dimension(ms_tot), intent(IN):: des_pos ! actually-positive design matrix (whether a case is present (1) or absent (0)

integer, intent(OUT):: num_truth_runs ! number of truth runs found in the dataset described by the datamatrices data* and
                                      ! design matrices des*
integer, intent(IN):: max_num_truth_runs ! Maximum number of truth runs to be considered, it should be set to mn+ms to make sure
                              ! that no collapsing is done to reduce the number to a smaller value

integer, intent(OUT):: ierror ! error flag: 0 -> fit worked properly ; 2 -> wrong design matrix values on input

real(kind=double), allocatable, dimension(:) :: data_neg! actually-negative used here, full array
real(kind=double), allocatable, dimension(:) :: data_pos! actually-positive used here, full array
integer, dimension(max_num_truth_runs):: cat ! categorical data

integer:: mn, ms

integer:: i, index

! count the cases actually used, first the actually-negative
mn = 0
do i = 1, mn_tot
     if(des_neg(i) == 1) then
              mn = mn + 1
     elseif(des_neg(i) == 0) then
     else
     ! if some values in the array des_neg is unacceptable (diff from 0 or 1) or
            ierror = 2
            return
     endif
enddo
! actually-positive
ms = 0
do i = 1, ms_tot
     if(des_pos(i) == 1) then
              ms = ms + 1
     elseif(des_pos(i) == 0) then
     else
     ! if some values in the array des_neg is unacceptable (diff from 0 or 1) or
            ierror = 2
            return
     endif
enddo

! allocate full arrays
allocate(data_neg(mn))
allocate(data_pos(ms))


! Load the data for the two different modalities, storing the locations in the new arrays of the input data
index = 0
do i = 1, mn_tot
     if(des_neg(i) == 1) then
            index = index+1
            data_neg(index) = data_neg_tot(i)
     endif
enddo

! Load the data for the two different modalities, storing the locations in the new arrays of the input data
index = 0
do i = 1, ms_tot
     if(des_pos(i) == 1) then
            index = index+1
            data_pos(index) = data_pos_tot(i)
     endif
enddo

! categorize data to find how many categories are present (if smaller than max_num_truth_runs).
call catgrz(1, mn, ms, 1, cat, data_neg, data_pos, num_truth_runs, max_num_truth_runs)


!----------------------------------------------------------------------------------------
end subroutine find_truthruns_d
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
SUBROUTINE LABROC4_Collapser(catn, cats, num_cat, debug_flag, new_cat_index)
! PURPOSE: receive categorical data as input and eliminate empty categories or
!          categories that have the same truth as their neighboring categories.
!          This is done because the elimination of the category boundaries between
!          those categories does not affect the Maximum Likelihood estimation of
!          the curve parameters or of the other cutoffs.
!          In general we call categorize if we the method is applied to list data
!          while we call it collapse if it is applied to categorical data. the
!          operation is mathematically almost isomorph, but computationally they
!          work differently because they are used differently.
!          NOTE: No check on the number of categories is done here. Most fitting routines
!          have a maximum number of categories that they can deal with and those numbers
!          are different.
! Created march 21st 2006 by Lorenzo Pesce, starting from the code of Kevin Schartz.
! ---------------------------------------------------------------------
!
!
! Example:
!
!          Original categorical data matrix:
!
!              1     2     3     4     5
!           _______________________________
!           |     |     |     |     |     |
!           |  1  |  1  |  0  |  2  |  0  |
!           |_____|_____|_____|_____|_____|
!           |     |     |     |     |     |
!           |  0  |  0  |  0  |  0  |  1  |
!           |_____|_____|_____|_____|_____|
!
!
!           After the first pass: Deleting columns
!
!              1     2     3     5
!           _________________________
!           |     |     |     |     |
!           |  1  |  1  |  2  |  0  |
!           |_____|_____|_____|_____|
!           |     |     |     |     |
!           |  0  |  0  |  0  |  1  |
!           |_____|_____|_____|_____|
!
!
!           Second pass: Collapsing columns
!
!              1     3
!           _____________
!           |     |     |
!           |  4  |  0  |
!           |_____|_____|
!           |     |     |
!           |  0  |  1  |
!           |_____|_____|
!
!           This is the final form of the "matrix".
!           new_cat_index will be (1,1,1,1,2)
!           the other arrays won't be affected.
!
! =====================================================================

 USE Debugging

 implicit none

! dummy variables.
 integer, intent(in)::  num_cat                     ! Number of rating categories (number of cutpoints + 1)
 integer, dimension(num_cat), intent(in) ::  catn   ! Array of frequencies for normal cases
 integer, dimension(num_cat), intent(in) ::  cats   ! Array of frequencies for abnormal cases
 integer, intent(in):: debug_flag ! whether (1) or not (0) to write debugging information out
 integer, dimension(num_cat), intent(out)::  new_cat_index ! Keep track of the categories to which old
                                                    ! categories are collapsed.

! internal variables
 integer:: new_cat ! index of the current category in the collapsed array
 integer:: icat ! loop index, on the input categories
 integer:: prev_cat !  last non empty category

! Start writing the data on entry in the collapser, if requested to do so
 if (debug_flag==1) then
      continue! write (debugUnit, *) '--------------------------------------------------'
      continue! write (debugUnit, *) 'Start LABROC4_Collapser ...'
 endif

 ! The data in the first category goes to the new first category, independently from
 ! its nature (including being empty)
 new_cat_index (1) = 1
 new_cat = 1
 prev_cat = 1
! Determine the location of the corners.
 category_loop:  do icat =  2, num_cat
         ! If the current category is empty, it gets the same category number as the
         ! previous category
         if(  catn(icat) == 0 .and. cats(icat) == 0 ) then
                 prev_cat = prev_cat ! the last non zero category remains the same
                 new_cat = new_cat

         ! If the category has the same positive truth as the previous one, it goes together
         ! with the previous category as well
         elseif(   ( catn(icat) == 0 .and. catn(prev_cat) == 0 ) .and. &
                   ( cats(icat) /= 0 .and. cats(prev_cat) /= 0 )         ) then
                 prev_cat = icat
                 new_cat = new_cat

         ! If the category has the same negative truth as the previous one, it goes together
         ! with the previous category as well
         elseif(   ( cats(icat) == 0 .and. cats(prev_cat) == 0 ) .and. &
                   ( catn(icat) /= 0 .and. catn(prev_cat) /= 0 )         ) then
                 prev_cat = icat
                 new_cat = new_cat

         ! Check whether all the previous categories (one or more)  were empty categories, in which case the current
         ! category needs to be merged with the previous one.
         elseif( all( catn(1:icat-1) == 0) .and. all(cats(1:icat-1) == 0 ) ) then
                 prev_cat = prev_cat
                 new_cat = new_cat

         ! Otherwise, we have a new category index because the current category is different from
         ! from the previous one in labroc4 sense.
         else
                 prev_cat = icat
                 new_cat = new_cat + 1
         endif

         ! Load the collapsed category index for the current category
         new_cat_index (icat) = new_cat

  enddo category_loop


  if (debug_flag==1) then

      continue! write (debugUnit, *) '--------------------------------------------------'
      continue! write (debugUnit, *) 'Number of categories', num_cat
      continue! write (debugUnit, *) 'Indices of categories after collapsing,'
      continue! write (debugUnit, *) 'categories with the same index are collapsed'
      continue! write (debugUnit, *) '--------------------------------------------------'
      continue! write (debugUnit, *) ' Input   # neg.    # pos.   New Cat.     '
      continue! write (debugUnit, *) ' cat.      cases     cases    index '
      do icat = 1, num_cat
              continue! write (debugUnit,"(4(i6,3x))") icat,  catn(icat), cats(icat), new_cat_index(icat)
      enddo

      continue! write (debugUnit, *) '--------------------------------------------------'

  endif

!----------------------------------------------------------------------------------------
end subroutine LABROC4_Collapser
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------


!----------------------------------------------------
 SUBROUTINE CATGRZ(POSITIVEISLARGE, NUM_NORMAL_CASES,NUM_ABNORMAL_CASES,idebug, CAT0, &
                   NEG_INPUT, POS_INPUT,NUM_CATEGORIES,MAX_NUM_CATEGORIES,CASE_CAT)
!----------------------------------------------------
!
! CATEGORIZE POSSIBLY CONTINUOUS TEST RESULT DATA  INTO CATEGORICAL  DATA USING THE LABROC4/LABROC5 ALGORITHM
! SEE Continuous ROC paper for algorithmic and validation information
 !CHAGES: LP 2008 dec U of C : changed the size of the array CASE_CAT, to make unnecessary to use a predefined
!                          maximum number of cases
! CHANGE: LP U of C January 19th 2010: split the summations as the old one was wrongly only on the negative cases
 implicit none

! DUMMY VARIABLES (SUBRoutINE PARAMETERS)
 INTEGER, INTENT(IN):: POSITIVEISLARGE ! Likelihood of abnormal TEST RESULT VALUE associated with larger values
 INTEGER, INTENT(IN):: NUM_NORMAL_CASES
 INTEGER, INTENT(IN):: NUM_ABNORMAL_CASES
 INTEGER, INTENT(IN):: idebug ! whether to write a log file or not
 REAL(kind=DOUBLE),INTENT(IN),DIMENSION(NUM_NORMAL_CASES)::NEG_INPUT ! negative TEST RESULT VALUEs are stored
 REAL(kind=DOUBLE),INTENT(IN),DIMENSION(NUM_ABNORMAL_CASES)::POS_INPUT ! negative TEST RESULT VALUEs are stored

 INTEGER, INTENT(IN) :: MAX_NUM_CATEGORIES ! MAXIMUM NUMBER OF CATEGORIES ALLOWED BY THE DIMENSIONING IN THE
                                           ! MAIN PROGRAM, THE MODULE WILL SEEK TO PRODUCE MAX_NUM_CATEGORIES.
                                           ! IF LESS THAT THAT ARE AVAILABLE (SAY N_CAT) IT WILL RETURN N_CAT
                                           ! IF MORE ARE AVAILABLE, IT WILL RETURN MAX_NUM_CATEGORIES

 INTEGER, INTENT(out), DIMENSION(ACT_NEG:ACT_POS, MAX_NUM_CATEGORIES) :: CAT0 ! CONTAINS THE CATEGORIES
                    ! CREATED BY THIS CATEGORIZATION ALGORITHM ON EXIT
 INTEGER, INTENT(out):: NUM_CATEGORIES ! THE NUMBER OF CATEGORIES FOUND
 INTEGER, OPTIONAL, INTENT(out),DIMENSION(ACT_NEG:ACT_POS,max(NUM_NORMAL_CASES,NUM_ABNORMAL_CASES))::CASE_CAT
                   ! This array stores for each case the category where it is allocated. Mostly to be
                   ! used by MRMC schemes, and this is why it is optional

! INTERNAL VARIABLES
!   NOTE: 6/4/03 THESE ARRAYS ARE NOT DIMENSIONED EXACTLY TO FIT THEIR DATA, BUT RATHER TO FIT THE SIZE OF THE
!                INPUT CURRENTLY USED. THIS PREVENTS FROM HAVING TO CREATE HUGE ARAYS TO ACCOMODATE POSSIBLE LARGE
!                INPUTS, WHICH WOULD CREATE PROBLEMS FOR PEOPLE RUNNING SMALL DATASETS ON SMALL OR CROWDED
!                COMPUTERS, OR TO HAVE TO USE POYNTERS OR LINKED LISTS ALL OVER THE PLACE. THE SMALL DIMENSIONALITY
!                OF THESE ARRAYS (ESSENTIALLY 1-DIMENSIONAL) RENDER THE MEMORY WASTE OF THIS CONFIGUARATION
!                NON RRELEVANT. ARRAYS ARE NOT DEALLOCATED BECAUSE F90 DOES THAT AUTOMATICALLY.

 REAL(kind=DOUBLE),DIMENSION(2,NUM_NORMAL_CASES + NUM_ABNORMAL_CASES)     :: ALL_CASES
 REAL(kind=DOUBLE), DIMENSION(ACT_NEG:ACT_POS,2 , max(NUM_NORMAL_CASES,NUM_ABNORMAL_CASES)):: val_index
 REAL(kind=DOUBLE), DIMENSION(ACT_NEG:ACT_POS,2 , max(NUM_NORMAL_CASES,NUM_ABNORMAL_CASES)):: val_index2
 REAL(kind=DOUBLE),DIMENSION(TST_RES_VAL:act_pos,NUM_NORMAL_CASES + NUM_ABNORMAL_CASES) :: VALUE_CATEGORIES
 REAL(kind=DOUBLE),DIMENSION(TST_RES_VAL:act_pos,NUM_NORMAL_CASES + NUM_ABNORMAL_CASES) :: TRUTH_RUNS_CATEGORIES
              ! values of the
              ! test result value which separate two different truth runs, with the number of
              ! actually normal and abnormal cases in that truth run
 REAL(kind=DOUBLE), DIMENSION(2, NUM_NORMAL_CASES + NUM_ABNORMAL_CASES) :: CITY_BLOCK_DIST ! THIS ARRAY IS A DOUBLE ARRAY:
                                                  ! city_block_dist(1,:) = TPF+FPF FOR THIS CATEGORY BOUNDARY
                                                  ! city_block_dist(2,:) = IS A FLAG FOR THE PROCEDURE WHICH CREATES THE
                                                  !              FINAL CATEGORY. IF IT IS + 1 THE OPERATING
                                                  !              POINT IS TAKEN, IF NOT IT IS REJECTED
 REAL(kind=DOUBLE), DIMENSION(2, NUM_NORMAL_CASES + NUM_ABNORMAL_CASES) :: FPF_TPF
 REAL(kind=DOUBLE), ALLOCATABLE, DIMENSION(:)   :: TST_RES_VAL_CATEGORY_BOUNDARY
                                                  !This array contains the test result values that separate categories from
                                                  ! from each other


 INTEGER:: NUM_CASES ! Total number of cases, positive plus negative
 INTEGER:: NUM_DIFFERENT_VALUES ! number of different TEST RESULT VALUEs of this data set
 INTEGER:: NUM_TRUTH_RUNS ! Number of truth state runs found in the raw data.
 INTEGER:: NUM_BOUNDARIES ! number of boundaries between truth runs or categories (N - 1)
 INTEGER:: I, INDEX_POS, INDEX_NEG ! Innocent loop variable and some index variables

! EXECUTABLE INSTRUCTIONS START

! open output file which contains information about the steps the algorithm takes to
! create categorical data from the continuoud test results values


 NUM_CASES = NUM_NORMAL_CASES + NUM_ABNORMAL_CASES

! Create an array where the normal and abnormal cases follow
! each other but have a different index (-1 or +1, respectively). OPERATION DONE USING ARRAY SYNTAX
! IF SMALLER VALUES SHOW STRONGER EVIDENCE OF SIGNAL(abnormality), TAKE THE OPPOSITE OF EACH VALUE
! SO THAT NOW LARGER VALUES WILL SHOW STRONGER EVIDENCE OF SIGNAL
! Note that val_index has the same shape as in, but the same data as ALL_CASES. the
! second dimension is to allow for values and index

 IF (POSITIVEISLARGE==0) THEN
   ALL_CASES(1,1:NUM_NORMAL_CASES)              = -NEG_INPUT(1:NUM_NORMAL_CASES)
   ALL_CASES(1,NUM_NORMAL_CASES+1:NUM_CASES)    = -POS_INPUT(1:NUM_ABNORMAL_CASES)
 ELSE
   ALL_CASES(1,1:NUM_NORMAL_CASES)              =  NEG_INPUT(1:NUM_NORMAL_CASES)
   ALL_CASES(1,NUM_NORMAL_CASES+1:NUM_CASES)    =  POS_INPUT(1:NUM_ABNORMAL_CASES)
 ENDIF
 ! Load the truth
 ALL_CASES(2,1:NUM_NORMAL_CASES) =  -1.0_double
 ALL_CASES(2,NUM_NORMAL_CASES+1:NUM_CASES) = 1.0_double


! Sort the values independently from their truth state, sort in ascending value
! keep track of the original order in the array val_index
 CALL SORT(NUM_CASES,ALL_CASES)


! Collapse all the cases with the same test result value in the array value_categories
! The ascending order in test result values is maintained in the array VALUE_CATEGORIES
CALL COLLAPSE_IDENTICAL_VALUES( VALUE_CATEGORIES, NUM_DIFFERENT_VALUES, ALL_CASES, NUM_CASES )


! Take value_categories and find the truth state runs of this data set, the boundaries between these
! truth state runs, and the city block distance between these points and 0,0 on an ROC plot
! (which corresponds to infinity on the test results value axis.
! The ith boundary is between the current truth runs and the one with index i + 1, i.e, the cases have
! TEST RESULT VALUE larger equal to the value with their index and smaller the one before them
! (apart from the first one which is simply any measured value largeer that the value of the boundary).
! Ascending order is preserved.
CALL  CREATE_TRUTH_RUNS( NUM_TRUTH_RUNS,CITY_BLOCK_DIST(1,:),TRUTH_RUNS_CATEGORIES,  &
      NUM_NORMAL_CASES,NUM_ABNORMAL_CASES,NUM_DIFFERENT_VALUES,VALUE_CATEGORIES, FPF_TPF )


! NUM_TRUTH_RUNS IS THE NUMBER OF POINTS GENERATED BY THE LABROC4 ALGORITHM. NOW
! THE SUBRoutINE APPLIES THE LABROC5 ALGORITHM TO  COLLAPSE THE TOTAL NUMBER OF CATEGORIES IF
! THEY ARE TOO MANY. IF THE NUMBER OF TRUTH RUNS IS MORE THAN MAX_NUM_CATEGORIES, MAX_NUM_CATEGORIES POINTS
! WILL BE SELECTED. OTHERWISE ALL OF THEM WILL BE KEPT. 8/12/03 LP-UC


! INITIALIZE CITYBLOCK DISTANCE TO "REJECTED OPERATING POINT" STATUS
 CITY_BLOCK_DIST(2,:) = -1.0_double

IF (NUM_TRUTH_RUNS > MAX_NUM_CATEGORIES) THEN
     NUM_CATEGORIES = MAX_NUM_CATEGORIES
     CALL FIND_LABROC5_POINTS(NUM_TRUTH_RUNS, NUM_CATEGORIES, idebug, CITY_BLOCK_DIST)
ELSE ! In this case the number of categories is the same as the number of truth runs
     ! so all the points need to be selcted (set array value to 1, see initialization above
     NUM_CATEGORIES = NUM_TRUTH_RUNS
     CITY_BLOCK_DIST(2,2:NUM_TRUTH_RUNS) = 1.0_double
ENDIF


! The number of category boundaries is one less than the number of categories, in the do loop
! the boundaries are loaded, than 1 is added to the total number to get the number of categories
ALLOCATE(TST_RES_VAL_CATEGORY_BOUNDARY(NUM_CATEGORIES-1))



NUM_BOUNDARIES = 0
DO  I = 1, NUM_TRUTH_RUNS

  IF (CITY_BLOCK_DIST(2,I) > 0.0_double ) THEN
     NUM_BOUNDARIES = NUM_BOUNDARIES + 1
     TST_RES_VAL_CATEGORY_BOUNDARY(NUM_BOUNDARIES) = TRUTH_RUNS_CATEGORIES(TST_RES_VAL,I)
  ENDIF

ENDDO

! THIS CHECK SHOULD BE REINSTATED ANY TIMES ONE WANTS TO FIDDLE WITH THE CATEGORIZATION
! SCHEME
!    CHECK IF THE algorithm used is actually consistent, SINCE THE 2 NUMBERS COME FROM
!    DIFFFERENT CALCULATIONS.
!IF (NUM_CATEGORIES /= NUM_BOUNDARIES + 1) THEN
!  WRITE(*,*) "ERROR::SUBROUTINE CATGRZ"
!  WRITE(*,*) "ERROR:: MISMATCH BETWEEN NUMBER OF CATEGORIES AND NUMBER OF BOUNDARIES ", &
!              NUM_CATEGORIES, NUM_BOUNDARIES
!  STOP
!ENDIF

! CONVERT TEST RESULTS VALUES INTO CATEGORICAL DATA, I.E., ALL THE CASES WHICH LAY BETWEEN TWO BOUNDARIES
! ARE COUNTED AND ATTRIBUTED TO  THE CATEGORY WITH INDEX EQUAL TO THE SMALLEST BOUNDARY

CALL CREATE_ORDINAL_CATEGORIES(NUM_TRUTH_RUNS,TRUTH_RUNS_CATEGORIES,NUM_CATEGORIES, &
     TST_RES_VAL_CATEGORY_BOUNDARY,CAT0)

! Reverse the order of the matrix of the boundaries to
!the original value (see the beginning of this subroutine)
IF (POSITIVEISLARGE==0)  TST_RES_VAL_CATEGORY_BOUNDARY = -TST_RES_VAL_CATEGORY_BOUNDARY


! Verify whether the program requested the category by case data
category_by_case: IF( PRESENT(CASE_CAT) ) then

   IF (POSITIVEISLARGE==0) THEN
       val_index(act_neg, 1, 1:NUM_NORMAL_CASES)    = -NEG_INPUT(1:NUM_NORMAL_CASES)
       val_index(act_pos, 1, 1: NUM_ABNORMAL_CASES) = -POS_INPUT(1:NUM_ABNORMAL_CASES)
   ELSE
       val_index(act_neg, 1, 1:NUM_NORMAL_CASES)    =  NEG_INPUT(1:NUM_NORMAL_CASES)
       val_index(act_pos, 1, 1: NUM_ABNORMAL_CASES) =  POS_INPUT(1:NUM_ABNORMAL_CASES)
   ENDIF
   ! Load the original indexing
   val_index(act_neg, 2, 1:NUM_NORMAL_CASES)    = (/  (1.0_double*I, I = 1, NUM_NORMAL_CASES) /)
   val_index(act_pos, 2, 1: NUM_ABNORMAL_CASES) = (/  (1.0_double*I, I = 1, NUM_ABNORMAL_CASES) /)


   ! Sort the truths independtly, to be used later to build the case indexing
   CALL SORT(NUM_NORMAL_CASES,val_index(act_neg,:,:))
   CALL SORT(NUM_ABNORMAL_CASES,val_index(act_pos,:,:))


   ! Use a scratch val_index to load in the first index the case index, we don't need the value
   ! anymore (we can find the value with the index, so we don't need it)
   val_index2(act_neg, 1, 1: NUM_NORMAL_CASES)   =  val_index(act_neg, 2, 1: NUM_NORMAL_CASES)
   val_index2(act_pos, 1, 1: NUM_ABNORMAL_CASES) =  val_index(act_pos, 2, 1: NUM_ABNORMAL_CASES)

   ! Load the category number in the cases, since categories are in order and the values are
   ! also in increasing value, we can apply the categories straight to the arrays, without
   ! worrieg about how we made them
   INDEX_NEG = 1
   INDEX_POS = 1
   DO I = 1, NUM_CATEGORIES
      if ( CAT0(ACT_NEG,I) > 0)   val_index2(act_neg, 2, INDEX_NEG: INDEX_NEG + CAT0(ACT_NEG,I) - 1 ) = real(I)
      if ( CAT0(ACT_POS,I) > 0)   val_index2(act_pos, 2, INDEX_POS: INDEX_POS + CAT0(ACT_POS,I) - 1 ) = real(I)
      INDEX_NEG = INDEX_NEG +  CAT0(ACT_NEG,I)
      INDEX_POS = INDEX_POS +  CAT0(ACT_POS,I)
   ENDDO


   ! Sort the arrays by case index. Now we have the category by case data we need
   CALL SORT(NUM_NORMAL_CASES,val_index2(act_neg,:,:))
   CALL SORT(NUM_ABNORMAL_CASES,val_index2(act_pos,:,:))


! CHANGE: LP U of C January 19th 2010: split the summations as the old one was wrongly only on the negative cases
   forall(i=1:NUM_NORMAL_CASES)
     CASE_CAT(ACT_NEG, i)   = int( val_index2(ACT_NEG, 2, i) )
   end forall


   forall(i=1:NUM_ABNORMAL_CASES)
     CASE_CAT(ACT_POS, i)   = int( val_index2(ACT_POS, 2, i) )
   end forall

endif category_by_case


! Close the file with the categorization data

 DEALLOCATE(TST_RES_VAL_CATEGORY_BOUNDARY)


!----------------------------------------------------
END SUBRoutINE CATGRZ
!----------------------------------------------------
!---------------------------------------------------

!------------------------------------------------
!-------------------------------------------------
 SUBRoutINE CREATE_TRUTH_RUNS (NUM_TRUTH_RUNS,CITY_BLOCK_DIST,TRUTH_RUNS_CATEGORIES,&
      NUM_NORMAL_CASES,NUM_ABNORMAL_CASES,NUM_DIFFERENT_VALUES,VALUE_CATEGORIES, FPF_TPF)
!------------------------------------------------
! Takes the array value_categories,
! VALUE_CATEGORIES(TST_RES_VAL,I) =  TEST RESULT VALUE
! VALUE_CATEGORIES(act_neg,I) =  number of normal cases with this test result value
! VALUE_CATEGORIES(act_pos,I) =  number of negative cases with this test result value
! and uses it to identify the truth state runs present in these data. Se beginning of the module
! for explanation

 implicit none
! It assumes that the larger the test result value the most likely is that the case will be abnormal or positive


! INPUT DUMMY VARIABLES
 INTEGER,INTENT(IN):: NUM_NORMAL_CASES
 INTEGER,INTENT(IN):: NUM_ABNORMAL_CASES
 INTEGER,INTENT(IN):: NUM_DIFFERENT_VALUES ! number of different TEST RESULT VALUEs of this data set
 REAL(kind=DOUBLE),DIMENSION (TST_RES_VAL:act_pos,NUM_NORMAL_CASES+NUM_ABNORMAL_CASES), INTENT(IN):: VALUE_CATEGORIES


! outPUT DUMMY VARIABLES
 INTEGER,INTENT(out):: NUM_TRUTH_RUNS
 REAL(kind=DOUBLE),DIMENSION (NUM_NORMAL_CASES+NUM_ABNORMAL_CASES),INTENT(out)::CITY_BLOCK_DIST
                       ! City block distance
 REAL(kind=DOUBLE),DIMENSION (TST_RES_VAL:act_pos,NUM_NORMAL_CASES+NUM_ABNORMAL_CASES), INTENT(out)::TRUTH_RUNS_CATEGORIES
              ! values of the
              ! test result value which separate two different truth runs, with the number of
              ! actually normal and abnormal cases in that truth run
 REAL(kind=DOUBLE),DIMENSION (2,NUM_NORMAL_CASES+NUM_ABNORMAL_CASES),INTENT(out):: FPF_TPF ! This contains the FPF and TPF for
              ! EACH operating point, which would be a boundary between categories (and in terms of indices, contains the
              ! category with the same index


 ! Work variables
 INTEGER I_VALUE
 REAL(kind=DOUBLE)::NUM_NORMAL_IN_TRUTH_RUN  ! counter for the number of normal cases in a truth run
                                             ! it is real because of the array it is stored into
 REAL(kind=DOUBLE)::NUM_ABNORMAL_IN_TRUTH_RUN ! counter for the number of abnormal cases in a truth run
                                              ! it is real because of the array it is stored into
 REAL(kind=DOUBLE)::FPF_INC ! the value the FTF increments every time there is a normal case
 REAL(kind=DOUBLE)::TPF_INC ! the value the TPF increments every time there is a abnormal case
 REAL(kind=DOUBLE)::TPF
 REAL(kind=DOUBLE)::FPF

 ! This is how much are the FPF and TPF increased by adding one normal or abnormal case respectively.
 FPF_INC = 1.0_double/NUM_NORMAL_CASES
 TPF_INC = 1.0_double/NUM_ABNORMAL_CASES

! INITIALIZE VARIABLES IN WHICH ITERATIVE SUMS ARE MADE
 NUM_NORMAL_IN_TRUTH_RUN = 0.0_double
 NUM_ABNORMAL_IN_TRUTH_RUN = 0.0_double
 FPF = 0.0_double
 TPF = 0.0_double
 NUM_TRUTH_RUNS = 0


 ! LOOP OVER THE DIFFERENT VALUES TO FIND THE TRUTH RUNS. IN ROC ANALYSIS THER REFERENCE POINT IS
 ! THE ONES WHERE THE PROBABILITY OF ABNORMAILTY IS LARGEST, WHICH IN THIS CASE IS ASSUMED TO BE FOR LARGER VALUES
 ! SO, SINCE THE ARRAYS ARE ORDERED IN ASCENDING ORDER, THE COUNTER WILL MAKE A COUNT DOWN. THIS MEANS THAT THE TRUTH RUNS
 ! ARRAYS WILL BE IN OPPOSITVE ORDER AS THE OTHER ARRAYS. FOR THIS REASON AT THE END THE NEW ARRAYS ORDERING IS INVERTED

 ! THE ALGORITHM BASICALLY LOOKS FOR A CHANGE IN TRUTH RUN (SEE COMMENTS AT THE BEGINNING OF THE MODULE, THE REFERENCES THERE
 ! AND THE COMMENTS RIGHT BELOW HERE), WHILE SUMMING ALL THE VALUES FOUND BETWEEN THE PREVIOUS TRUTH RUN BOUNDARY THE ONE
 ! IT IS LOOKING FOR. ONCE THE BOUNDARY IS FOUND THE SUMMED VALUES ARE LOADED IN THE NEW TRUTH RUN DATA ARRAYS.
 OVER_DIFFERENT_VALUES: DO I_VALUE = NUM_DIFFERENT_VALUES, 2 , -1
   ! If one value has both normal and abnormal cases, than it is a truth run of its own
    FPF = FPF + VALUE_CATEGORIES(act_neg,I_VALUE) * FPF_INC
    TPF = TPF + VALUE_CATEGORIES(act_pos,I_VALUE) * TPF_INC
    NUM_NORMAL_IN_TRUTH_RUN = NUM_NORMAL_IN_TRUTH_RUN + VALUE_CATEGORIES(act_neg,I_VALUE)
    NUM_ABNORMAL_IN_TRUTH_RUN = NUM_ABNORMAL_IN_TRUTH_RUN + VALUE_CATEGORIES(act_pos,I_VALUE)

   ! THERE ARE TWO SITUATIONS WHERE A BOUNDAY BETWEEN TWO NEIGHBORING VALUE-CATEOGRIES COULD BE:
   ! A)  IF THE NUMBER OF BOTH THE NORMAL AND THE ABNORMAL CASES IS DIFFERENT FROM 0
   !     IN EITHER OF THE TWO VALUE-CATEGORIES.
   ! B)  IF THE TWO VALUE-CATEGORIES HAVE A DIFFERENT TRUTH STATE.
   !     ONCE CONDITION A) IS VERIFIED TO BE FALSE  ALL WE HAVE TO CHECK IF EITHER THE TWO VALUE-CATEGORIES,
   !     BUT NOT BOTH, DOES NOT HAVE CASES FOR ONE OF THE TRUTH STATES. THIS DECISION STATEMENT IN FACT PARTIALLY
   !     INCLUDES POINT A) (WHEN ONE IS MIXED, AND THE OTHER ONE IS NORMAL), BUT THIS HAS NO EFFECT ON THE CATEGORIZATION
   !     AS LONG AS THERE IS AN .OR. STATEMENT IN BETWEEN

   IF(   ( & ! FIRST CHECK IF EITHER THE VALUE_CATEGORY HAS *BOTH* TRUTH STATES
              ( (VALUE_CATEGORIES(act_neg,I_VALUE) .spne. 0.0_double)     &
               .AND.                                                    &
               (VALUE_CATEGORIES(act_pos,I_VALUE) .spne. 0.0_double) )    &
            .OR.                                                        &
               ( (VALUE_CATEGORIES(act_neg,I_VALUE-1) .spne. 0.0_double)  &
              .AND.                                                     &
                 (VALUE_CATEGORIES(act_pos,I_VALUE-1) .spne. 0.0_double) ) &
          )                                                             &
       .OR.                                                             &
          ( & ! THEN CHECK IF THE TWO HAVE DIFFERENT TRUTH STATE AS AN ALTERNATIVE
            (VALUE_CATEGORIES(act_neg,I_VALUE) .spne. 0.0_double)         &
          .NEQV.                                                        &
            (VALUE_CATEGORIES(act_neg,I_VALUE-1) .spne. 0.0_double)       &
          )                                                             &
    ) THEN
              NUM_TRUTH_RUNS = NUM_TRUTH_RUNS + 1
              FPF_TPF(1,NUM_TRUTH_RUNS)= FPF
              FPF_TPF(2,NUM_TRUTH_RUNS)= TPF
              CITY_BLOCK_DIST(NUM_TRUTH_RUNS) = FPF + TPF
              ! THE BOUNDARY IS BETWEEN THE SMALLEST TEST RESULT VALUE OF THE CURRENT TRUTH RUN AND
              ! THE NEXT SMALLER VALUE
              TRUTH_RUNS_CATEGORIES (TST_RES_VAL,NUM_TRUTH_RUNS) =              &
                    (  VALUE_CATEGORIES(TST_RES_VAL,I_VALUE) +                &
                    VALUE_CATEGORIES(TST_RES_VAL,I_VALUE-1)  ) / 2.0_double

              TRUTH_RUNS_CATEGORIES (act_neg,NUM_TRUTH_RUNS) = NUM_NORMAL_IN_TRUTH_RUN
              TRUTH_RUNS_CATEGORIES (act_pos,NUM_TRUTH_RUNS) = NUM_ABNORMAL_IN_TRUTH_RUN

              NUM_NORMAL_IN_TRUTH_RUN = 0.0_double
              NUM_ABNORMAL_IN_TRUTH_RUN = 0.0_double
    ENDIF

 ENDDO OVER_DIFFERENT_VALUES

! Note that in any case here we have to increase NUM_TRUTH_RUNS by one. Either there was a boundary right before
! Which means that also this one is a truth run, or there wasn't, which means that we need to create a truth tun
! for this value category and the ones which came before it

 FPF = FPF + VALUE_CATEGORIES(act_neg, 1) * FPF_INC
 TPF = TPF + VALUE_CATEGORIES(act_pos, 1) * TPF_INC
 NUM_NORMAL_IN_TRUTH_RUN = NUM_NORMAL_IN_TRUTH_RUN + VALUE_CATEGORIES(act_neg,1)
 NUM_ABNORMAL_IN_TRUTH_RUN = NUM_ABNORMAL_IN_TRUTH_RUN + VALUE_CATEGORIES(act_pos,1)

 NUM_TRUTH_RUNS = NUM_TRUTH_RUNS + 1

 FPF_TPF(1,NUM_TRUTH_RUNS)= FPF
 FPF_TPF(2,NUM_TRUTH_RUNS)= TPF
 CITY_BLOCK_DIST(NUM_TRUTH_RUNS) = FPF + TPF
 TRUTH_RUNS_CATEGORIES (TST_RES_VAL,NUM_TRUTH_RUNS) =  VALUE_CATEGORIES(TST_RES_VAL, 1 )
 TRUTH_RUNS_CATEGORIES (act_neg,NUM_TRUTH_RUNS) = NUM_NORMAL_IN_TRUTH_RUN
 TRUTH_RUNS_CATEGORIES (act_pos,NUM_TRUTH_RUNS) = NUM_ABNORMAL_IN_TRUTH_RUN

! Since at this point the truth runs are order in descending order of TEST RESULT VALUEs while
! the rest of the arrays are in ascending order, we turn them in ascending order too
 FORALL(I_VALUE = 1:NUM_TRUTH_RUNS)
    FPF_TPF(1,I_VALUE)  =   FPF_TPF(1,NUM_TRUTH_RUNS - I_VALUE + 1)
    FPF_TPF(2,I_VALUE)  =   FPF_TPF(2,NUM_TRUTH_RUNS - I_VALUE + 1)
    CITY_BLOCK_DIST(I_VALUE) =  CITY_BLOCK_DIST(NUM_TRUTH_RUNS - I_VALUE + 1)
    TRUTH_RUNS_CATEGORIES (TST_RES_VAL,I_VALUE) = TRUTH_RUNS_CATEGORIES (TST_RES_VAL,NUM_TRUTH_RUNS - I_VALUE + 1)
    TRUTH_RUNS_CATEGORIES (act_neg,I_VALUE) = TRUTH_RUNS_CATEGORIES (act_neg,NUM_TRUTH_RUNS - I_VALUE + 1)
    TRUTH_RUNS_CATEGORIES (act_pos,I_VALUE) = TRUTH_RUNS_CATEGORIES (act_pos,NUM_TRUTH_RUNS - I_VALUE + 1)
 END FORALL


!-------------------------------------------------
 END SUBRoutINE CREATE_TRUTH_RUNS
!------------------------------------------------
!-------------------------------------------------


!------------------------------------------------
!-------------------------------------------------
SUBRoutINE COLLAPSE_IDENTICAL_VALUES (VALUE_CATEGORIES,NUM_DIFFERENT_VALUES,ALL_CASES, NUM_CASES)
!------------------------------------------------
! Subroutine stores in the array VALUE_CATEGORIES the number of normal and abormal cases per each value
! of the test result value. In this way tied values will be in the same index
! so that:
! VALUE_CATEGORIES(TST_RES_VAL,I) =  TEST RESULT VALUE
! VALUE_CATEGORIES(act_neg,I) =  number of normal cases with this test result value
! VALUE_CATEGORIES(act_pos,I) =  number of negative cases with this test result value
! Note that this transformation preserves the information contained in the array ALL_CASES

 implicit none

! These are the input variables
 INTEGER,INTENT(IN):: NUM_CASES !total number of cases, inclusing positive and negative
 REAL(kind=DOUBLE),DIMENSION (2,NUM_CASES),INTENT(IN)::ALL_CASES ! Sorted array with all the cases
           ! in ascending order by TEST RESULT VALUE, second dimension in truth state (+1 or -1)

 REAL(kind=DOUBLE),DIMENSION (TST_RES_VAL:act_pos,NUM_CASES), INTENT(out):: VALUE_CATEGORIES ! computed array, see
            ! header of subroutine for its shape and meaning
 INTEGER,INTENT(out):: NUM_DIFFERENT_VALUES ! number of different TEST RESULT VALUEs of this data set

 ! Work variables
 INTEGER:: ICASE ! index pointing to the current case


 REAL(KIND=DOUBLE):: CURRENT_TEST_VALUE_RESULT

 ! initialize the array
 VALUE_CATEGORIES = 0.0_double
 ! inizialize the counter to one, since the different values are at least one
 NUM_DIFFERENT_VALUES = 1
 ! initialize the first value to the largest value, and load it into the relevant array
 CURRENT_TEST_VALUE_RESULT = ALL_CASES(1,1)
 VALUE_CATEGORIES (TST_RES_VAL,1) = CURRENT_TEST_VALUE_RESULT
 ! loop over all the cases to categorize them by their value
 COLLAPSING : DO ICASE = 1,NUM_CASES
   ! Check if the value of the currently analyzed cases is the same as the value of the current
   ! value-category, if not, create a new value-category
   IF(  ALL_CASES(1, ICASE) .spne. CURRENT_TEST_VALUE_RESULT ) THEN
      NUM_DIFFERENT_VALUES = NUM_DIFFERENT_VALUES + 1
      CURRENT_TEST_VALUE_RESULT = ALL_CASES(1,ICASE)
      VALUE_CATEGORIES (TST_RES_VAL,NUM_DIFFERENT_VALUES) = CURRENT_TEST_VALUE_RESULT
   ENDIF
   ! since we are pointing to the right category by now, increase the positive or negative count of
   ! this category by one, depending on the truth state of the current case
   IF ( ALL_CASES(2,ICASE) < 0.0_double) THEN
      VALUE_CATEGORIES (act_neg,NUM_DIFFERENT_VALUES) = &
                              VALUE_CATEGORIES (act_neg,NUM_DIFFERENT_VALUES) + 1
   ELSE
      VALUE_CATEGORIES (act_pos,NUM_DIFFERENT_VALUES) =  &
                              VALUE_CATEGORIES (act_pos,NUM_DIFFERENT_VALUES) + 1
   ENDIF

 ENDDO COLLAPSING


END SUBROUTINE COLLAPSE_IDENTICAL_VALUES
!------------------------------------------------
!------------------------------------------------


!--------------------------------------------
  SUBROUTINE SORT(NUM_VALUES,THE_ARRAY)
!--------------------------------------------
! SORT THE TWO VALUED REAL (KIND = *) ARRAY THE_ARRAY OF DIMENSIONS (2,NUM_VALUES)
! IN ASCENDING ORDER ACCORDING TO ITS FIRST VALUES:: THE_ARRAY(1,:)
! THE FUNCTIONING OF THIS SUBROUTINE DEPENDS UPON THE MODULES USED BY THE CATEGORIZATION MODULE

 implicit none

 INTEGER,INTENT(IN):: NUM_VALUES
 REAL(kind=DOUBLE),INTENT(INOUT),DIMENSION(2,NUM_VALUES)::THE_ARRAY ! ON entry the raw data with their qualification
                    ! at the second index, on exit the same values but in ascending order (according to the first
                    ! value

 ! INTERNAL VARIABLES
 INTEGER::I
! INTEGER,DIMENSION(1)::J ! This is a hack, because some compilers (INTEL at least) cannot tell a size one
!                         ! one dimensional array from a scalar.
 REAL(kind=DOUBLE), DIMENSION(2,NUM_VALUES)::TEMP1 ! Copy of raw data
 INTEGER, DIMENSION(NUM_VALUES)::VAL_IND ! indexes
! REAL(kind=DOUBLE),DIMENSION(2):: TEMP


! Loop over all the values and set the current value as the smallest
! No need to cycle over the last one, it will be ranked outmatically
! NOTE: a forall statement would not work, because of array operations precedence which would not
!       allow THE_ARARAY TO BE UPDATED on the right.
! DO I = 1, NUM_VALUES - 1
!      J =  I - 1 + MINLOC( THE_ARRAY(1,I:NUM_VALUES) ) !MINLOC gives the minimum of the current
!      TEMP =  THE_ARRAY(:,I)
!      THE_ARRAY(:,I) =  THE_ARRAY(:, J(1))
!      THE_ARRAY(:, J(1)) =  TEMP
! ENDDO

! Load the data in a temp array
TEMP1 = THE_ARRAY
! Find the indexing using quicksort
call qsortd(TEMP1(1,:), VAL_IND, NUM_VALUES)

forall(I=1:NUM_VALUES) THE_ARRAY(:,I) = TEMP1(:,VAL_IND(I))



!---------------------!--------------------------------------------
 END SUBRoutINE SORT
!---------------------!--------------------------------------------
!---------------------!--------------------------------------------

!-----------------------------------------

   SUBRoutINE  FIND_LABROC5_POINTS(NUM_TRUTH_RUNS, NUM_CATEGORIES, idebug, CITY_BLOCK_DIST)
!-----------------------------------------
!
!  Select a NUM_CATEGORIES - 1 big subset of the labroc4 corners (which are contained in the array CITY_BLOCK_DISTANCE)
!  (NOTE: this is because there are n-1 boundaries between n categories.)
!
!  The LABROC5 algorithm is based upon city block distances (the sum of the TPF and FPF for each corner) after they are
!  streched/compressed according to a normal distribution
! (see reference the Continuous ROC paper at the beginning of this module)

!  Note that the city block distance goes inversely with the boundary/threshold value, since we are
! assuming that the larger values have higher probability of abnormality or positivity.
! Thus, the lower is the threshold the more  cases will be included, and
!  accordingly the larger the TPF and FPF.  This is the reason why arrays are worked on "backwards", i.e., to produce
! FPF and TPF increasing in magnitude.

! Notation follow the reference " Continuous ROC paper" at the beginning of module whenever possible.

! ALGORITHM NOTES:
! The algorithm is a 2 step process which attempts to set the points as uniformly distributed as possible on a stretched
! scale (the scale is stretched so that distances between operating points that lie close to the -45 degrees diagonal
! on the ROC plot are stretched, while for points that are far from that diagoal they are compressed.)

! The two steps are:
! 1) set points close to some ideal uniform distances and at least as far from each other as some value (this is the "coverage
!    phase", because it is meant to cover the ROC space, possibly at the expense of leaving large gaps between chosen points)
! 2) If we do not get enough points after step a), some points are added in between the previous selected points,
!    trying to remove the largest gaps (see the subroutine "ADD_ONE_POINT_IN_BETWEEN for details about this.) "filling phase"

! The algorithm here described is slightly different from the original one in the sense that instead of starting from
! one corner ( (0,0) or (1,1) ) it starts from the point closest to the -45 degrees diagonal on the
! ROC plot). Moreover, the algorithm looks for the  point closest to the ideal, under the condition that the new point is more
! than some number of categories dependent distance from the previously selected one. This distance is half
! then the previous algorithm and measured from the previously selected point not from the previous ideal point.
! After the previous phase points are chosen accoding to the add_one_point_in_between subroutine, which is different from
! the one described in the paper too. (see the subroutine "ADD_ONE_POINT_IN_BETWEEN" for details)

 USE statistic_functions ! for the distribution functions

 implicit none

 INTEGER,INTENT(IN):: NUM_TRUTH_RUNS ! Total # of truth runs included the point which is in (1,1)
 INTEGER,INTENT(INout):: NUM_CATEGORIES
 INTEGER,INTENT(IN):: idebug
 REAL(kind=DOUBLE),DIMENSION (2,NUM_TRUTH_RUNS), INTENT(INout)::CITY_BLOCK_DIST !
                                                  ! city_block_dist(1,:) = TPF+FPF FOR THIS CATEGORY BOUNDARY
                                                  ! city_block_dist(2,:) = IS A FLAG FOR THE PROCEDURE WHICH CREATES THE
                                                  !              FINAL CATEGORY. IT IS SET TO + 1 IF THE OPERATING
                                                  !              POINT IS SELECTED.

! work variables
 REAL(kind=DOUBLE):: D_R_minus_1 ! max city block distance, (see reference the Continuous ROC paper at tge beginning of module)
 REAL(kind=DOUBLE):: D_1         ! min city block distance, (see reference the Continuous ROC paper at tge beginning of module)

 REAL(KIND=DOUBLE),PARAMETER :: RANGE_D = 1.96_double ! This is the value the city block distances
                  ! are normalized to before being mapped to a normal distribution. RIght now it is set to 2 sigmas (1.96)

 INTEGER:: new_points_search_start ! The subroutine uses this variable to decide from which truth state run to start
                                   ! selecting the LABROC5 points. The coverage algorithm might start from 1 or 2 depending
                                   ! upon where the first point is (it might be preselcted, before the stretching.)
 INTEGER:: new_points_search_end   ! This is the same as tnew_points_seach_start, but for large indices.
                                   ! this value is NUM_TRUTH_RUNS - 1 at most since the last point is
                                   ! the fixed point (1,1), and it should not be selected.

 REAL(kind=DOUBLE):: STEP          ! If we look for N operating points, STEP is the nearest neighbor distance for N uniformly
                                   ! distributed points
 REAL(kind=DOUBLE),DIMENSION(1:NUM_TRUTH_RUNS):: D_R_PRIME ! array containing the stretched city block distances
 REAL(KIND=DOUBLE):: D_R_PRIME_IDEAL_LARGE
 REAL(KIND=DOUBLE):: D_R_PRIME_IDEAL_SMALL
 INTEGER:: LAST_LARGE_TRUTH_RUN
 INTEGER:: LAST_SMALL_TRUTH_RUN
 INTEGER,DIMENSION(1)::  FORTRAN_BUG ! This array is used to correct for the incapability of the intel compiler to use
                       ! size one one dimensional arrays as scalars
 INTEGER::   POINT_CLOSEST_CENTER ! location of point closest to the -45 degrees diagonal

 LOGICAL:: FOUND_LARGE, FOUND_SMALL ! if one point was found for either larger or smaller stretched city block distances
 INTEGER:: SELECTED_BOUNDARIES ! Number of boundaries selected by the LABROC5 algorithm

 INTEGER:: I   ! LOOP COUNTER
 INTEGER:: POINTS_TO_BE_ADDED  ! THIS IS THE NUMBER OF OPERATING POINT WHICH NEEDS TO BE ADDED IN CASES
            ! THE FIRST PART OF THE SEARCH DOES NOT PROVIDE ENOUGH, SEE HEADER OF SUBRoutINE FOR CLARIFICATION

 D_R_PRIME = 0.0_double


! THIS CHECK SHOUD BE REINSTATED IF ONE WANTS TO MODIFY THE CATEGORIZATION SCHEME TO MAKE SURE
! THAT IT IS WORKING PROPERLY
! FIRST CHECK IF WE HAVE ENOUGH TRUTH RUNS TO FIND THE NUMBER OF CATEGORIES WE ARE LOOKING
! FOR (I.E., CHECK FOR CONSISTENCY IN INPUT)
! IF(NUM_CATEGORIES > NUM_TRUTH_RUNS) THEN
!  WRITE(*,*) "ERROR::SUBROUTINE FIND_LABROC5_POINTS"
!  WRITE(*,*) "ERROR:: THE NUMBER OF CATEGORIES REQUESTED > NUMBER OF CATEGORIES AVAILABLE"
!  STOP
! ENDIF



SELECTED_BOUNDARIES = 0 ! we haven't selected anything so far

D_1                               =  0.0_double
NEW_POINTS_SEARCH_END             =  NUM_TRUTH_RUNS


! Note that the largest meaningful value for a new point is one less that the total
! number of truth state runs, since the last one is *always* 1,1 by contruct so it has 2 as
! city block distance and it is of no insterest (it rovides no information)

D_R_minus_1                       =  2.0_double
D_R_PRIME(1)                      = 1.0_double
NEW_POINTS_SEARCH_START           = 2


! Create the array D'r (D_R_prime) as defined in the "continuoous ROC paper". This is equivalent to map the city block
! distances to a normal distribution, centered around 1, which is the city block distance of any point on the -45 degrees
! diagonal on an ROC plot. This mapping stretches distances as pointed out earlier

! NOTE: The denominators ( e.g., "( 0.5_double - PHI(- RANGE_D) )"  ) correspond  in eqn. 16 of Continuous ROC  paper
! to .475. Which was an approximation, we use the exact value instead
! Which allows us also more flexibility (if we change the value of D_range, it will authomatically
! be corrected

 WHERE( city_block_dist(1,NEW_POINTS_SEARCH_START:NEW_POINTS_SEARCH_END) < 1.0_double)
       D_R_PRIME(NEW_POINTS_SEARCH_START:NEW_POINTS_SEARCH_END) = &
           (    PHI(  - RANGE_D * ( 1.0_double - city_block_dist (1,NEW_POINTS_SEARCH_START:NEW_POINTS_SEARCH_END)  ) / &
                               ( 1.0_double -                                     D_1                             )   &
                    ) &
               - 0.5_double   )/ ( 0.5_double - PHI(- RANGE_D) )
 ELSEWHERE
       D_R_PRIME(NEW_POINTS_SEARCH_START:NEW_POINTS_SEARCH_END) = &
           (    PHI(  RANGE_D * ( city_block_dist (1,NEW_POINTS_SEARCH_START:NEW_POINTS_SEARCH_END)  - 1.0_double ) / &
                               ( D_R_minus_1                                                        - 1.0_double )   &
                    )  &
               - 0.5_double   )/   ( PHI(RANGE_D) - 0.5_double )
 END WHERE


! Select the point closest to the -45 degrees diagonal. NOTE that the awkward syntax is caused by the incapability of
! the INTEL compiler to consider size one 1-dimensional arrays as scalars. I kept 2 variables because I hope that the bug
! will be removed and I will be able to live with a scalar.

  FORTRAN_BUG =  MINLOC( ABS( D_R_PRIME(NEW_POINTS_SEARCH_START:NEW_POINTS_SEARCH_END) ))
  POINT_CLOSEST_CENTER =  FORTRAN_BUG(1)
  ! Check if the middle point isn't also an extreme (for very skewed data sets)
  IF( CITY_BLOCK_DIST(2,POINT_CLOSEST_CENTER) .speq. -1.0_double) THEN
    CITY_BLOCK_DIST(2,POINT_CLOSEST_CENTER) = 1.0_double ! setting this value to 1 impliest that it is selected
    SELECTED_BOUNDARIES =   SELECTED_BOUNDARIES + 1      ! One more point was se;ected
  ENDIF


! Set the step size, this is used to determine how far each point should be from the previous one. Note that
! the size of the step depends upon the number of points that we still need to select and not on the total
! number of points. If we could select uniformly distributed points, this would be their nearest-neighbor distance.

  step = 2.0_double/ ( NUM_CATEGORIES - SELECTED_BOUNDARIES)

! Now start the first part of the search, the "coverage phase". The algorithm looks for two points: one with city block
! distance larger than the point closest to the -45 degrees diagonal
! and one with distance smaller than the point closest to the -45 degrees diagonal. After that it keeps searching for 2 points
! one aobve and one below.
! First it sets ideal points,  then it looks for the labroc4 points closest to them, then it checks if they are far  enough from
! the previously selected point. Note that the search above and below are independent and one can end before the other one
! depending upon how many points do they have.

! Set the first 2  candidate ideal points
  D_R_PRIME_IDEAL_LARGE = D_R_PRIME(POINT_CLOSEST_CENTER) + step
  LAST_LARGE_TRUTH_RUN =  POINT_CLOSEST_CENTER - 1

  D_R_PRIME_IDEAL_SMALL = D_R_PRIME(POINT_CLOSEST_CENTER) - step
  LAST_SMALL_TRUTH_RUN = POINT_CLOSEST_CENTER + 1

! To understand this algorithm one has to keep in mind that when the category boundary goes up, the index goes up, but the
! city block distance goes down and the points becomes closer to (0,0) on the ROC plot. This is because the cumulative
! distributions of the ROC curves are taken from + infinity to the current threshold, or category boundary related to the
! operating point.

 TRUTH_RUNS: DO
     IF(LAST_LARGE_TRUTH_RUN >= 1) THEN
        CALL  FIND_ONE_BOUNDARY( LAST_LARGE_TRUTH_RUN,   NEW_POINTS_SEARCH_START, - 1, &
               D_R_PRIME_IDEAL_LARGE , FOUND_LARGE)
     ELSE
        FOUND_LARGE = .FALSE.
     ENDIF

     IF (FOUND_LARGE) THEN
            ! Increment the next ideal point, making sure that it is not too close to the 1,1 vertex
            D_R_PRIME_IDEAL_LARGE = MIN(D_R_PRIME_IDEAL_LARGE + step, 1.0_double - .5_double * step)
            CITY_BLOCK_DIST(2,LAST_LARGE_TRUTH_RUN) = 1.0_double ! Tag the cutoff J as a labroc5 cutoff
            SELECTED_BOUNDARIES = SELECTED_BOUNDARIES + 1
            LAST_LARGE_TRUTH_RUN = LAST_LARGE_TRUTH_RUN - 1
     ENDIF

    ! If we have found enough so far, no need to look for another small one
    IF( SELECTED_BOUNDARIES >= NUM_CATEGORIES - 1)   EXIT TRUTH_RUNS

    IF(LAST_SMALL_TRUTH_RUN <= NEW_POINTS_SEARCH_END) THEN
       CALL FIND_ONE_BOUNDARY( LAST_SMALL_TRUTH_RUN,   NEW_POINTS_SEARCH_END, +1, &
                D_R_PRIME_IDEAL_SMALL ,FOUND_SMALL)
    ELSE
       FOUND_SMALL = .FALSE.
    ENDIF

     IF (FOUND_SMALL) THEN
            ! Decrement the next ideal point, making sure that it is not too close to the 0,0 vertex
            D_R_PRIME_IDEAL_SMALL = MAX(D_R_PRIME_IDEAL_SMALL - step, - 1.0_double + .5_double * step)
            CITY_BLOCK_DIST(2,LAST_SMALL_TRUTH_RUN) = 1.0_double ! Tag the cutoff J as a labroc5 cutoff
            SELECTED_BOUNDARIES = SELECTED_BOUNDARIES + 1
            LAST_SMALL_TRUTH_RUN = LAST_SMALL_TRUTH_RUN + 1
     ENDIF

    ! Check if we have all the points we need or if we are not going to find any anymore with the "coverage" algorithm
    IF( SELECTED_BOUNDARIES >= NUM_CATEGORIES - 1 .OR. (.NOT.FOUND_SMALL .AND. .NOT.FOUND_LARGE) ) THEN
      EXIT TRUTH_RUNS
    ENDIF

 ENDDO TRUTH_RUNS

!  CHECK IF THERE ARE ENOUGH OPERATING POINTS. IF THERE AREN'T START THE "FILLING PHASE" UNTIL ENOUGH
!  POINTS ARE SELECTED.


 IF(SELECTED_BOUNDARIES < NUM_CATEGORIES - 1) THEN
    POINTS_TO_BE_ADDED = (NUM_CATEGORIES - 1) - SELECTED_BOUNDARIES
    DO I = 1, POINTS_TO_BE_ADDED
      CALL ADD_ONE_POINT_IN_BETWEEN(NEW_POINTS_SEARCH_START, NEW_POINTS_SEARCH_END, CITY_BLOCK_DIST,D_R_PRIME,idebug)
      SELECTED_BOUNDARIES = SELECTED_BOUNDARIES + 1
    ENDDO

 ENDIF

CONTAINS
     SUBRoutINE FIND_ONE_BOUNDARY(START_TRUTH_RUN,END_TRUTH_RUN,INCREMENT,D_R_PRIME_IDEAL,FOUND)
     ! THIS IS AN INTERNAL PROCEDURE TO THE SUBRoutINE FIND_LABROC5_POINTS.
     ! IT LOOKS FOR A POINT WHICH IS AS CLOSE AS POSSIBLE TO AN IDEAL POINT BUT FARTHER AWAY FROM
     ! THE PREVIOUSLY SELECTED POINT BY AT LEAST HALF THE QUANTITY CALLED STEP Defined in the calling subroutine.)
     ! IF IT FINDS A POINT, IT RETURNS THE ADDRESS OF THE POINT AND SETS THE VARIABLE "FOUND" TO .TRUE.
     ! NOTE THAT IT WORKS BOTH BY SEARCHING UP AND DOWN THE STRETCHED DISTANCE ARRAY, SINCE IT TAKES
     ! AN INCREMENT AS INPUT.
     ! NOTE: it is an internal subroutine because in this way it can share the subroutine's global variables and
     !       we can make the call leaner.

      implicit none

     INTEGER,INTENT(INout):: START_TRUTH_RUN
     INTEGER,INTENT(IN):: INCREMENT
     INTEGER,INTENT(IN):: END_TRUTH_RUN
     REAL(KIND=DOUBLE),INTENT(IN):: D_R_PRIME_IDEAL
     LOGICAL,INTENT(out):: FOUND
     INTEGER :: CURRENT_TRUTH_RUN

      FOUND = .FALSE.
      ! loops over the remaining LABROC4 points, moving up or down depending upon the call
      DO CURRENT_TRUTH_RUN = START_TRUTH_RUN, END_TRUTH_RUN - INCREMENT , INCREMENT
          ! The algorithm keeps skipping points until neither of the three following conditions is fulfilled:
           IF(                                                                    &
               ! 1) the following point is closer to the ideal point than this point
               (  ABS( D_R_PRIME( CURRENT_TRUTH_RUN + INCREMENT )  - D_R_PRIME_IDEAL )    &
               <                                                                   &
                 ABS( D_R_PRIME(CURRENT_TRUTH_RUN) - D_R_PRIME_IDEAL )                     )       &
              .OR. &
               ! 2) THE POINT IS CLOSER TO THE PREVIOUSLY SELECTED POINT THAN ...
               (  ABS(D_R_PRIME( CURRENT_TRUTH_RUN) - D_R_PRIME( START_TRUTH_RUN - INCREMENT)) < .5_double* step)     &
              .OR. &
               ! 3) IT IS TOO CLOSE TO THE ROC PLOT CORNERS.
               ( 1.0_double -  ABS(D_R_PRIME(CURRENT_TRUTH_RUN)) <  .5_double * step )     &
             ) CYCLE
             ! If this pooint wasn't skipped, it means it is good.
             FOUND = .TRUE.
             EXIT

      ENDDO

      IF (FOUND) THEN
         START_TRUTH_RUN = CURRENT_TRUTH_RUN
      ELSE
         IF( CURRENT_TRUTH_RUN  == END_TRUTH_RUN & ! The do loop was exited at the largest value, which
            .AND.   &                              ! then it is at the minimum distance (if acceptable). We take it.
            (ABS(D_R_PRIME( CURRENT_TRUTH_RUN) - D_R_PRIME( START_TRUTH_RUN - INCREMENT)) < .5_double* step)     &
            .AND. &
            ( 1.0_double -  ABS(D_R_PRIME(CURRENT_TRUTH_RUN)) <  .5_double * step )     &
           ) THEN
            START_TRUTH_RUN = END_TRUTH_RUN
            FOUND = .TRUE.
         ENDIF
      ENDIF

     END SUBRoutINE FIND_ONE_BOUNDARY

 END SUBRoutINE FIND_LABROC5_POINTS
!--------------------------------------------------------------
!---------------------------------------------------------------

!--------------------------------------------------
SUBRoutINE CREATE_ORDINAL_CATEGORIES(NUM_TRUTH_RUNS,TRUTH_RUNS_CATEGORIES,NUM_CATEGORIES, &
     TST_RES_VAL_CATEGORY_BOUNDARY,CAT0)
!--------------------------------------------------
!
!  CONVERT TEST RESULT VALUES INTO ORDINAL CATEGORICAL DATA.
!  It is based on the assumtion that the arrays are sorted in INCREASING test result value starting from
! the index 1. The array TST_RES_VAL_CATEGORY_BOUNDARY has element which separate category in the following way
! category(i) = x/x, TST_RES_VAL_CATEGORY_BOUNDARY[i] < x < TST_RES_VAL_CATEGORY_BOUNDARY[i-1]
! It sums all the actually postive together and all the actually negative together, for all the test result values which are
! between 2 boundary values which were selected before by previously called subroutines.

 implicit none

  INTEGER, INTENT(in):: NUM_TRUTH_RUNS
  REAL(kind=DOUBLE), INTENT(IN),DIMENSION(TST_RES_VAL:act_pos,NUM_TRUTH_RUNS):: TRUTH_RUNS_CATEGORIES
              ! values of the test result value which separate two different truth runs, with the number of
              ! actually normal and abnormal cases in that truth run
  INTEGER,INTENT(IN):: NUM_CATEGORIES ! total number of categories
  REAL(kind=DOUBLE), INTENT(IN),DIMENSION(NUM_CATEGORIES-1):: TST_RES_VAL_CATEGORY_BOUNDARY
                                                  !This array contains the test result values that separate categories from
                                                  ! from each other
  INTEGER, INTENT(out),DIMENSION(ACT_NEG:ACT_POS,NUM_CATEGORIES):: CAT0 ! ORDINAL CATEGORICAL DATA: CONTAINTS THE
                  ! NUMBER OF ACTUALLY POSITIVE (ABNORMAL) AND ACTUALLY NEGATIVE (NORMAL) CASES PER CATEGORY.

  INTEGER:: I_CATEGORY,I_TRUTH_RUN

! Initialize the category array
  CAT0 = 0

! Start the category counter, which is used to update the boundary/threshold value and the category array
  I_CATEGORY = NUM_CATEGORIES

 ! Since in ROC integration is started at + infinity, the arrays are scanned backwards
 DO I_TRUTH_RUN = NUM_TRUTH_RUNS, 1, -1
     ! when test result values are larger than the following category boundary, the category count is increased by 1
     ! all the cases betwee 2 boundaries are summmed together in the corresponding category.
     IF(TRUTH_RUNS_CATEGORIES(TST_RES_VAL,I_TRUTH_RUN) < TST_RES_VAL_CATEGORY_BOUNDARY (I_CATEGORY-1) ) I_CATEGORY = I_CATEGORY-1
     IF (I_CATEGORY == 1) EXIT
     CAT0(ACT_NEG,I_CATEGORY) =  CAT0(ACT_NEG,I_CATEGORY) +  INT( TRUTH_RUNS_CATEGORIES(ACT_NEG,I_TRUTH_RUN) )
     CAT0(ACT_POS,I_CATEGORY) =  CAT0(ACT_POS,I_CATEGORY) +  INT( TRUTH_RUNS_CATEGORIES(ACT_POS,I_TRUTH_RUN) )
 ENDDO
  ! sum all the remaining cases in the last category
  CAT0(ACT_NEG,1) =  INT( SUM(TRUTH_RUNS_CATEGORIES(ACT_NEG,1:I_TRUTH_RUN) ) )
  CAT0(ACT_POS,1) =  INT( SUM(TRUTH_RUNS_CATEGORIES(ACT_POS,1:I_TRUTH_RUN) ) )

 END SUBRoutINE CREATE_ORDINAL_CATEGORIES
!-------------------------------------------------------
!------------------------------------------------------

!-----------------------------------------------
  SUBRoutINE ADD_ONE_POINT_IN_BETWEEN(FIRST_TRUTH_RUN,LAST_TRUTH_RUN,CITY_BLOCK_DIST,D_R_PRIME,idebug)
!-----------------------------------------------
! THE SUBRoutINE SELECTS A NEW POINT FROM THE SET OF THE  LABROC4 POINTS THAT HAVE NOT BEEN SELECTED YET.
! IT LOOKS FOR THE LARGEST GAP BETWEEN 2 PREVIOUSLY SELECTED POINTS, AND TRIES TO PICK A POINT THERE, UNLESS THE MINIMAL
! DISTANCE OF THIS NEW POINT TO ALREADY SELECTED POINTS (I.E., THE DISTANCE BETWEEN THIS POINT AND THE CLOSEST OF THE ALREADY
! SELECTED POINTS) IS LESS THAN HALF OF THE LARGEST MINIMAL DISTANCE BETWEEN ANY UNSELECTED POINT AND AN ALREADY
! SELECTED ONE (I.E., LOOK AT ALL THE NON SELECTED POINTS, AND FIND THE ONE WHICH IS FATHER APART FROM ITS CLOSEST ALREADY
! SEECTED POINT).

! THIS ALGORITHM FUNCTIONS BY ASSUMING THAT THE POINT WITH A NEGATIVE CITY_BLOCK_DISTANCE(2, POINT_INDEX) ARE POINTS
! WHICH WERE NOT PREVIOUSLY SELECTED, WHILE POINTS WITH A POSITIVE VALUE WERE PREVIOUSLY SELECTED.

! NOTE THAT THE SUBRoutINE USES ASSUMED SIZE FOR THE ARRAY D_R_PRIME TO STREAMLINE THE CALL
!      HANDLE WITH CARE OR YOUR MIGHT OVERRUN THE BOUNDARIES.

 implicit none

 INTEGER, INTENT(IN) :: FIRST_TRUTH_RUN! SMALLEST INDEX OF A TRUTH RUN TO CONTROL FOR NON SELECTED POINTS
 INTEGER, INTENT(IN) :: LAST_TRUTH_RUN ! LARGEST INDEX OF A TRUTH RUN TO CONTROL FOR NON SELECTED POINTS
 REAL(kind=DOUBLE),INTENT(IN), DIMENSION(*) :: D_R_PRIME ! ARRAY OF STRETCHED DISTANCES
 REAL(kind=DOUBLE), INTENT(INout),DIMENSION (2,*) :: CITY_BLOCK_DIST ! ARRAY OF CITY BLOCK DISTANCES, SEE CALLING RoutINE FOR
                   ! DESCRIPTION
 INTEGER, intent(in):: idebug

! INTERNAL VARIABLES
 INTEGER:: I_TRUTH_RUN ! COUNTER FOR THE CURRENTLY CHECKED TRUTH RUN
 INTEGER:: MAX_DELTA_LOWER_BOUND
 INTEGER:: MAX_DELTA_UPPER_BOUND
 INTEGER:: LAST_DELTA_LOWER_BOUND ! THE INDEX OF THE LOWER BOUND OF THE CURRENT INTERVAL OF NON SELECTED POINTS
 INTEGER:: LAST_DELTA_UPPER_BOUND ! THE INDEX OF THE UPPER BOUND OF THE CURRENT INTERVAL OF NON SELECTED POINTS
 REAL(kind=DOUBLE):: MAX_DELTA    ! MAXIMUM GAP BETWEEN TWO PREVIOUSLY SELECTED POINTS
 REAL(kind=DOUBLE):: LAST_DELTA   ! SIZE OF THE GAP BETWEEN SELECTED POINTS CURRENTLY ANALIZED
 REAL(kind=DOUBLE):: MIN_DELTA    ! MINIMUM DISTANCE FROM ITS BOUNDARIES FOR A NON SELECTED POINT IN THE LARGEST GAP TO BE
                                  ! A SELECTION CANDIDATE. THE DISTANCE DEPENDS UPON THE POINT MOST DISTANT FROM ITS BOUNDARY
                                  ! OR FROM THE LAST SELECTION CANDIDATE (SEE CODE)
 LOGICAL:: FOUND_LOWER_BOUND ! lOGICAL TELLING IF A LOWER BOUND WAS FOUND BY THE DO LOOP THAT SEARCHES FOR A LOWER BOUND
 INTEGER:: THE_POINT ! THE ADDRESS OF THE POINT WHICH WILL BE SELECTED
 REAL(kind=DOUBLE):: MAX_MIN ! THE LARGEST MINIMUM DISTANCE (DISTANCE FROM THE CLOSEST SELECTED POINT)
 REAL(kind=DOUBLE):: THIS_MIN_DIST
 INTEGER:: I_THIS_MIN_DIST

 MAX_DELTA_LOWER_BOUND = FIRST_TRUTH_RUN
 MAX_DELTA_UPPER_BOUND = LAST_TRUTH_RUN
 MAX_DELTA = 0.0_double
 MAX_MIN = 0.0_double

 I_TRUTH_RUN = FIRST_TRUTH_RUN

 SCAN_UNUSED_POINTS: DO

 ! FIND THE NEXT LOWER BOUND OF AN INTERVAL CONTANING NON SELECTED POINTS, I.E., POINTS WITH NEGATIVE SECOND
 ! ELEMENT OF subarray ARRAY CITY_BLOCK_DIST(2,LAST_DELTA_LOWER_BOUND:LAST_TRUTH_RUN)
    FOUND_LOWER_BOUND = .FALSE.

    LOWER_BOUND: DO LAST_DELTA_LOWER_BOUND = I_TRUTH_RUN, LAST_TRUTH_RUN
      ! TAKE POINT AS LOWER BOUND IF:
      IF( &
           !1) THIS IS THE FIRST POINT AND THE CITY_BLOCK_DISTANCE(2,LAST_DELTA_LOWER_BOUND) IS  NEGATIVE
            ( (LAST_DELTA_LOWER_BOUND == 2 ) .AND.  CITY_BLOCK_DIST(2, LAST_DELTA_LOWER_BOUND) < 0) &
         .OR. &
           !2) CHECK IF THE PREVIOUS POINT WAS POSITIVE AND THIS NEGATIVE, SO WE KNOW IT IS A LOWER BOUND FOR
           !   A SEQUENCE OF NON SELECTED POINTS
            (  CITY_BLOCK_DIST(2, LAST_DELTA_LOWER_BOUND - 1) > 0  .AND. CITY_BLOCK_DIST(2, LAST_DELTA_LOWER_BOUND) < 0) &
        )  THEN
           FOUND_LOWER_BOUND = .TRUE.
           EXIT LOWER_BOUND
      ENDIF
    ENDDO LOWER_BOUND

    ! IF IT CAN'T FIND A LOWER BOUND, IT MEANS THAT THERE ARE NO MORE INTERVALS.
    IF (.NOT.FOUND_LOWER_BOUND)  EXIT SCAN_UNUSED_POINTS

 ! FIND THE UPPER BOUND OF THE PREVIOUS INTERVAL CONTANING NON SELECTED POINTS
 ! If there is a lower bound, there must be an upper one, since even if the last lower boundary was the last point
 ! this last point will also be the upper bound.
    UPPER_BOUND:  DO LAST_DELTA_UPPER_BOUND = LAST_DELTA_LOWER_BOUND, LAST_TRUTH_RUN
      ! we have found the upper bound if:
      IF( &
        !1) We are at the last point and it has a negative value
            ( (LAST_DELTA_UPPER_BOUND == LAST_TRUTH_RUN) .AND.  CITY_BLOCK_DIST(2, LAST_DELTA_UPPER_BOUND) < 0) &
         .OR. &
        !2) It has a negative value and the following has a positive one
            (  CITY_BLOCK_DIST(2, LAST_DELTA_UPPER_BOUND) < 0  .AND. CITY_BLOCK_DIST(2, LAST_DELTA_UPPER_BOUND + 1) > 0) &
        )  EXIT UPPER_BOUND
    ENDDO UPPER_BOUND

     I_TRUTH_RUN = LAST_DELTA_UPPER_BOUND + 1! The next lower bound can be the same as the previous upper + 1, but no less than
                                             ! that

! DETERMINE THE SIZE OF THE CURRENT INTERVAL. SINCE THE LOWER AND UPPER BOUNDS PREVIOUSLY FOUND ARE
! POINTS WHICH WERE NOT SELECTED PREVIOUSLY, AND WE WANT THE DISTANCE BETWEEN THE POINTS WHICH WERE SELECTED ALREADY,
! WE TAKE THE LAST_DELTA_LOWER_BOUND-1 AND LAST_DELTA_UPPER_BOUND+1. NOTE THAT WE COMPUTE THE LOWER MINUS THE UPPER BECAUSE
! THE DISTANCE ON THE ROC PLOT GOES IN OPPOSITVE ORDER THAN THE DISTANCE ON THE REST RESULT VALUE

     IF(LAST_DELTA_UPPER_BOUND == LAST_TRUTH_RUN) THEN
        LAST_DELTA =  D_R_PRIME(LAST_DELTA_LOWER_BOUND - 1)  + 1.0_DOUBLE
     ELSE
        LAST_DELTA =  D_R_PRIME(LAST_DELTA_LOWER_BOUND - 1) -   D_R_PRIME(LAST_DELTA_UPPER_BOUND + 1)
     ENDIF
     IF (LAST_DELTA >= MAX_DELTA) THEN ! IF THE CURRENT GAP IS THE LARGEST, RECORD ITS SIZE AND EXTREMES
           MAX_DELTA = LAST_DELTA
           MAX_DELTA_LOWER_BOUND = LAST_DELTA_LOWER_BOUND
           MAX_DELTA_UPPER_BOUND = LAST_DELTA_UPPER_BOUND
     ENDIF

! CHECK ALL THE POINTS INSIDE THIS INTERVAL TO SEE IF ANY OF THEM IS THE MOST DISTANT FROM ITS CLOSEST ALREADY
! SELECTED POINT
     MAXIMUM_MINIMUM_DISTANCE: DO I_THIS_MIN_DIST = LAST_DELTA_LOWER_BOUND, LAST_DELTA_UPPER_BOUND
        ! FIND THE DISTANCE OF THE CURRENT POINT FROM THE CLOSEST EXTREME OF THIS INTERVAL
        IF(LAST_DELTA_UPPER_BOUND + 1 < LAST_TRUTH_RUN) THEN
              THIS_MIN_DIST =  MIN( -D_R_PRIME(LAST_DELTA_UPPER_BOUND + 1) + D_R_PRIME(I_THIS_MIN_DIST), &
                     -D_R_PRIME(I_THIS_MIN_DIST) + D_R_PRIME(LAST_DELTA_LOWER_BOUND - 1)  )
        ELSE
              THIS_MIN_DIST =  MIN( 1.0_double + D_R_PRIME(I_THIS_MIN_DIST), &
                    -D_R_PRIME(I_THIS_MIN_DIST) + D_R_PRIME(LAST_DELTA_LOWER_BOUND - 1)  )
        ENDIF
        ! CHECK IF THIS DISTANCE IS THE LARGEST ONE SO FAR. IF IT IS RECORD THIS AS THE CURRENTLY CANDIDATE FOR
        ! THE NEXT POINT TO BE SELECTED
       IF(THIS_MIN_DIST >= MAX_MIN) THEN
        THE_POINT = I_THIS_MIN_DIST
        MAX_MIN = THIS_MIN_DIST
       ENDIF
     ENDDO MAXIMUM_MINIMUM_DISTANCE


 ENDDO SCAN_UNUSED_POINTS

 MIN_DELTA = MAX_MIN*.5_double ! SET REFERENCE DISTANCE FOR THE NEXT STEP

 ! LOOK INSIDE OF THE LARGEST GAP, AND FIND THE POINT WHICH IS MOST DISTANT FROM ITS CLOSEST EXTREME

  DO I_TRUTH_RUN = MAX_DELTA_LOWER_BOUND, MAX_DELTA_UPPER_BOUND
     LAST_DELTA = MIN(- D_R_PRIME(I_TRUTH_RUN) +  D_R_PRIME(MAX_DELTA_LOWER_BOUND-1), &
                       - D_R_PRIME(MAX_DELTA_UPPER_BOUND+1) + D_R_PRIME(I_TRUTH_RUN) )
     ! IF THIS MINIMUM DISTANCE IS LARGER THAN THE PRESELECTED VALUE, THE POINT BECOMES THE SELECTION CANDIDATE
     IF (LAST_DELTA >= MIN_DELTA) THEN
           MIN_DELTA = LAST_DELTA
           THE_POINT = I_TRUTH_RUN
     ENDIF

  ENDDO

  ! THE POINT WHICH IS RECORDED AS SELECTION CANDIDATE AT THIS POINT IS THE POINT WE WANT
  CITY_BLOCK_DIST(2, THE_POINT) = 1.0_double


!------------------------------------------------
  END SUBRoutINE ADD_ONE_POINT_IN_BETWEEN

!-------------------------------------------------------
!------------------------------------------------------
END MODULE CATEGORIZATION
