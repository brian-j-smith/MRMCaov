! File containing numerical constants used in various ROC algorithms.
! Either they are common constants or specific constants used to make the
! calculations more stable (e.g., preventing terms from becoming 0 or
! infinity)


! CONTAINS
! VARIABLE                  OPERATION
! =========                ==========
! par_almost_zero           used todetermine if |parameter| is too small to be
!                           considered different from zero
! p_q_closeto_one
! p_q_closeto_zero
! pi                        You guess
!
! DATE       PROGRAMMER        DESCRIPTION OF CHANGE
! ====       ==========        =====================
! 11/25/02   L. Pesce (UC)     Creation
! 12/02/02   L. Pesce (UC)     Added pi (using mathematica & memory)
! 10/24/07   L. Pesce (UC)     removed some terms and created the version for the ROC-project
!--------------------------------------------------
module computation_constants
!-----------------------------------------------------
use data_types,only: double
implicit none
private
public par_almost_zero
public p_q_closeto_zero,p_q_closeto_one, pi
public ncategory
integer, parameter :: ncategory = 400! maximum number of categories
                                     ! used by the labroc5 algorithm
real(kind = double), parameter:: par_almost_zero = 1.0e-11_double
! value is chosen so that the most unstable likelihoods (the one with 2
! category data) have a smooth profile before the cutoff
real(kind = double), parameter:: p_q_closeto_zero  = 1.0e-14_double
! Ideally P_q_close to one would be relate to p_q_closeto_zero, but because of
! compiler issues we had to define them independently...
real(kind = double), parameter:: p_q_closeto_one   = .99999999999999_double
real(kind = double), parameter:: pi   = 3.14159265358979323846_double

end module computation_constants
