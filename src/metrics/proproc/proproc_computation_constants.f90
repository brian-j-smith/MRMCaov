! File containing numerical constants used in the proper algorithm.
! In most situations the use of the constants has been greatly reduced
! and their value is less arbirary than it used to be

! CONTAINS
! VARIABLE                  OPERATION
! =========                ==========
! c_almost_zero             used todetermine if b is too small to be
!                           considered different from zero
! p_q_closeto_one
! p_q_closeto_zero
! pi                        You guess
!
! DATE       PROGRAMMER        DESCRIPTION OF CHANGE
! ====       ==========        =====================
! 11/25/02   L. Pesce (UC)     Creation
! 12/02/02   L. Pesce (UC)     Added pi (using mathematica & memory)
! 02/02/08   L. Pesce (UC)     Eliminated some constants, moved the to computation
!                              constants from the ROC project
! 03/05/10   L. Pesce (UC)     Removed references to computation constants to avoid double referencing 
!                              and render dependencies more explicit and straightforward.
!--------------------------------------------------
module proproc_computation_constants
!-----------------------------------------------------
use data_types,only: double
! the files need some additional polishing...
implicit none
private
public c_almost_zero, c_almost_one
public d_a_almost_zero, d_a_almost_Infinity
public max_c_2_d_a
! The value of c_almost_zero (and of its complement c_almost_one) is
! defined so as the largest value after which the bivariate binormal 
! distribution subroutine used fails when rho = ( 1 - c**2) / ( 1 + c**2)
real(kind = double), parameter:: c_almost_zero = 1.0e-11_double
real(kind = double), parameter:: c_almost_one = 1.0_double - c_almost_zero
real(kind = double), parameter:: d_a_almost_zero = 1.0e-8_double
real(kind = double), parameter:: d_a_almost_Infinity = 5.0_double
! value is chosen so that the most unstable likelihoods (the one with 2
! category data) have a smooth profile before the cutoff
real(kind = double), parameter:: max_c_2_d_a = 4.0e8_double ! maximal allowed ratio : c/d_a 
      ! It is used to prevent the optimizer for looking for solutions in the
      ! small d_a large c space where almost identical curves have different
      ! parameters

end module proproc_computation_constants
