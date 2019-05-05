! This file contains modules whose purpose is to share information between the procedures that 
! work around the optimization procedure used for the MLE. In this way we don't have to worry
! too much about what exactly does the optimizer allow the programs to exchange explicitly.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module labroc_initial_cutoffs
! this module is used to pass data around the optimizer brent
! to compute the initial estimates of the cutoffs. I did not want
! to change the canned procedure because it makes more sense to
! hide the data that should not be seen by the optimizer 
! LP - UC Janurary 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE data_types, only: double
  implicit none 

  real(kind=double):: a_par
  real(kind=double):: b_par
  real(kind=double):: pf,pt 
  integer:: cat_pos, cat_neg ! number of cases above the cutoff being initialized,
                           ! either positive or negative

end module labroc_initial_cutoffs

module labroc_one_pt
! module to work around optimizers for the one data fit
  USE data_types, only: double
  implicit none 
 real(kind = double):: a
 real(kind = double):: b

end module labroc_one_pt

module labroc_median_calc
! module to work around optimizers for the one data fit
  USE data_types, only: double
  implicit none 
 real(kind = double):: a
 real(kind = double):: b
 real(kind = double):: r

end module labroc_median_calc

