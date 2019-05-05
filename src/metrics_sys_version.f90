! File containing the procedures and modules used to report the current version
!
! Record of revisions:
! DATE       PERSON            DESCRIPTION OF CHANGE
! ====       ==========        =====================
! 09/02/09   L. Pesce (UC)     Creation

module libROC_version
  use data_types
  use io, only: line_length, line_format
  implicit none
  
  private
  
  public get_version_number

  ! THESE ARE THE CURRENT VERSION NUBERS, THE LIST OF THE CHANGES THAT MOVE FROM ONE VERSION TO THE NEXT IS AT THE BOTTOM OF THE MODULE, AFTER THE 
  ! "end module" STATEMENT. FIRST VERSION IS 0.9.0. 
  integer, parameter:: major = 1    ! major version number: significan change e.g., introduction of single curve, then ROCKIT, then MRMC
  integer, parameter:: minor = 0     ! minor version number: introduction of functionalities within tools, e.g. labroc added to single curve or making a change that is
                                     ! is not backward compatible.
  integer, parameter:: revision = 3 ! revision: bug removal => BACKWARD COMPATIBILITY IS TO BE EXPECTED

 
  contains

  !-----------------------------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------------------------
  subroutine get_version_number (maj, min, rev)
  !PURPOSE:  returns three numbers containing the sequence of major, minor and revision numbers as "major.minor.revision"
  !-----------------------------------------------------------------------------------------------------------------------------

   integer, intent(OUT) :: maj
   integer, intent(OUT) :: min
   integer, intent(OUT) :: rev

     maj = major
     min = minor
     rev = revision
  

  !-----------------------------------------------------------------------------------------------------------------------------
  end subroutine get_version_number

  !-----------------------------------------------------------------------------------------------------------------------------

end module libROC_version

! VERSION NUMBERS
! 0.9.0  -> Libraries at September 2nd 2009 start of the version "control"
! 0.10.0 -> October 2009 Introduced the non-parametric computation of FPF and TPF. Introduced their Confidence intervals using
!           "exact" methods -- see ../SRC_NONPARAM/roc_nonparametric.f90. The methods guarantee both the coverage and the correct
!           range of values (i.e., the do not overshoot 0 or 1)
! 0.10.1 -> added fit_degenerate/fit_perverse to the possible outcomes of a call to cvbmroc_0_0
! 0.10.2 -> added initialization of all return variables to 666.0 to failed or crazy fits to proproc_0_0/cvbmroc_0_0 
! 0.11.0 -> change the names of the libraries returning empirical operating points. There was one based on list data and one based on categorical
!           data, the first was part of SRC_PROPROC/proproc_out.f90 and the second was part of SRC_NONPARAM/roc_nonparametric.f90, now the are both 
!           in SRC_NONPARAM and are called empirical_operating_points_list and empirical_operating_points_cat. November 2009 
! 0.11.1 -> November 17th 2009, LP, changed the function compute_orig_cutoffs in ../SRC_PROPROC/proproc_functions.f90 to remove bugs caused by calling
!           proproc_0_0 with non-collapsed categorical data. The issues were not in the collapsing itself, therefore MLE curve estimates are unaffected
!           but in the cutoff reconstruction phase. A debug output was added from this function, allowing users to see whether a cutoff recontruction
!           was attempted and what its outcome was. The function was tested on various kinds of crummy input on OS X 10.4.11 (Xeon, rocuno).
!0.11.2 ->  added file closing in main drivers (proproc_0_0 and CvBMROC_0_0) to make sure that when the functions are called from any system they produce
!           a debug file. This is necessary at least on OS X 10.4.11 using gcc-4.3.0
!0.11.3 ->  Removed a bug from  the calling to the cutoff mapping for FPF, TPF and LOG_LIKE estimatio (had old input_data split by positives and negatives)
!           The bug was in the function catgrz in the module SRC_SYS/categorization.f90
!0.11.4 ->  Modyfied the function that computes the initial estimates for the conventional binormal model and eliminated the roc_functions.f90 module from
!           SRC_PROPROC and removed references to module computation_constants from module proproc_computation_constants.
!0.12.0 -> Introducting of CvBM ROC to replace LABROC4, with all its related functions.
!0.13.0 -> Stripped the executable binaries from the library and moved them into the SRC_BINS directory
!1.0.0 -> changed naming conventions to functions: proproc/proper -> PBM and labroc goes to CVBM. 
!  labroc/labroc4 refers now to the categorization that preserves the corners scheme
!  made the naming also more consistent, following function_model scheme (e.g., TPF_PBM,  FPF_find_TPF_PBM and so on)
!  NOTE: only the functions that are visible outside the modules follow the new naming conventions
!        changed the name of the likelihood related functions too like1 -> like_CVBM or like_PBM
!1.0.1 -> 4/14/2010 corrected a bug in the CvBM routine, it used to set the whole logical mask array after each fit instead
!         of the specific element.
!1.0.2 -> 5/02/2010 improved the  CvBM routine by changing the initial estimates
!         of the cutoffs (starting value changed and fixed-spacer now is 
!         dependent upon the size of the dataset
!1.0.3 -> 10/20/2010 fixed CVBM TPF function and its Var
