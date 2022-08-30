! The general parts were moved here by Lorenzo Pesce (U of C) on january 2010
! This file contains modules whose purpose is to share information between the procedures that
! work around the optimization procedure used for the MLE. In this way we don't have to worry
! too much about what exactly does the optimizer allow the programs to exchange explicitly.
module problem_data
! the data which are specific to the problem at hand
! and shared by all the procedures (and modified by none once they
! are set). The module is used to work around the various optimization
! routines to avoid exposing information that they don't need to see
  use debugging ! To import the logging procedures
  use io, only: line_length ! size of the line in the code

  use data_types, only : double

  implicit none

  integer:: num_normal_cases
  integer:: num_abnormal_cases
  integer:: num_cat ! Number of categories as found by catgrz
  integer, allocatable,  dimension(:):: catn, cats ! arrays containing categorical data
  integer:: pass_idebug ! whether debug data should be written to file or not

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module problem_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ll_derivatives
! derivatives of the p_i and q_i as a function of the
! parameters. These are necessary to compute both the
! gradient and the hessian at the expected values (or
! (expected information matrix, AKA Fisher information
! matrix). Module allows to work around the optimizer
! and save some CPU time by  not recomputing them
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use data_types, only: double
real(kind=double), allocatable, dimension(:):: p_val, q_val ! p & q values
real(kind=double), allocatable, dimension(:,:)  :: dp_dpar
real(kind=double), allocatable, dimension(:,:)  :: dq_dpar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ll_derivatives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
