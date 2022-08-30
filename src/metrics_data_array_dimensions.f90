! Module which contains the sizes of the arrays used in the
! code. Its destiny is to disappear once dynamic allocation is
! in place

 module array_dimensions
 implicit none

! integer, parameter :: ncase = 50000 ! Maximum number of cases accepted
 integer, parameter :: act_neg = 0    ! this is to locate in arrays the
                                      ! actually negative, or actually
                                      ! normal cases
 integer, parameter :: act_pos = 1    ! this is to locate in arrays the
                                      ! actually positive, or actually
                                      ! abnormal cases
 end module array_dimensions
