!----------------------------------------------------------------------
! File containing data the types standardization module
! WARNING: 
!         *There is no comparison between mixed types
!         *Check is needed to find out if IEEE double aren't quadruple
!          on a specific processor
!        * Need to add treatment for values like NAN INF and such 
!
! CONTAINS
! PROCEDURE                OPERATION
! ========                ==========
! eqsingle,nesingle        compare real single precision numbers (s_floating)
! eqdouble,nedouble        compare real double precision numbers (t_floating)

! DATE       PROGRAMMER        DESCRIPTION OF CHANGE
! ====       ==========        =====================
! 11/14/02   L. Pesce (UC)     Creation
! 11/19/02   L. Pesce (UC)     added comparison between reals
!                              for both double and single
! 05/18/04   L. Pesce (UC)     added overflow/underflow check (need to verify if it
!                              is stable across platforms and compilers  
!--------------------------------------------------
module data_types
!-----------------------------------------------------
implicit none
private
public :: single,double
public :: check_number ! identify overflows and underflows
public :: operator(.speq.),operator(.spne.)
! specify the data types for the programs
integer, parameter :: single = selected_real_kind(p=6) ! IEEE S_float +- 10^38, 6 digits
integer, parameter :: double = selected_real_kind(p=15) ! IEEE T_float +- < 10^308, 15 digits

!----------- Interfaces for comparison between reals
interface operator(.speq.)
 module procedure eqsingle, eqdouble ! remember to correct both if you wanna
                                     ! change anything
end interface operator(.speq.)

interface operator(.spne.)
 module procedure nesingle, nedouble ! remember to change both if you wanna
                                     ! change anything
end interface operator(.spne.)
!-------------------------------------------------


contains
!-------------------------------------------------
logical elemental function eqsingle(a,b)
!-------------------------------------------------
! Compares the difference between two numbers a,b, with the 
! possible numerical difference identifiable by the model
! used by the real a
implicit none

real(kind= single),intent(in):: a,b

 eqsingle = .not. nesingle(a,b)

!---------------------------------------------------
end function eqsingle
!---------------------------------------------------
!-------------------------------------------------
logical elemental function eqdouble(a,b)
!-------------------------------------------------
! Compares the difference between two numbers a,b, with the 
! possible numerical difference identifiable by the model
! used by the real a
implicit none

real(kind= double),intent(in):: a,b

eqdouble = .not. nedouble(a,b)

!---------------------------------------------------
end function eqdouble
!---------------------------------------------------

!-------------------------------------------------
logical elemental function nesingle(a,b)
!-------------------------------------------------
! Compares the difference between two numbers a,b, with the 
! possible numerical difference identifiable by the model
! used by the real a
implicit none
real(kind=single),intent(in):: a,b
intrinsic spacing 

 nesingle =  abs (a - b)  > 2.0_single * max(  spacing(a),spacing(b)  ) 

!---------------------------------------------------
end function nesingle
!---------------------------------------------------
!-------------------------------------------------
logical elemental function nedouble(a,b)
!-------------------------------------------------
! Compares the difference between two numbers a,b, with the 
! possible numerical difference identifiable by the model
! used by the real a
implicit none
real(kind=double),intent(in):: a,b
intrinsic spacing 

 nedouble =  abs (a - b)  > 2.0_double * max(  spacing(a),spacing(b)  )  

!---------------------------------------------------
end function nedouble
!---------------------------------------------------
!---------------------------------------------------
logical elemental function check_number(a)
! check in the numerical representation chosen (likely double)
! if the number is between huge and tiny which means that the
! number is neither an overflow nor an underflow. 
!WARNING: this might not work on all the architectures since 
!   the f95 standard does not tell what should the processor
!   do when such an exception occurs. If the underflow is reckoned as
!   zero, we can go ahead since 0 behaves properly. 5-5-05
implicit none
real(kind=double), intent(in):: a
!---------------------------------------------------

if( abs(a) < huge(1.0_double) .and. &
   ( abs(a)  > tiny(1.0_double) .or. (a .speq. 0.0_double) ) &
   ) then
 check_number = .true.
else
 check_number = .false.
endif


end function check_number
!--------------------------------------------------
!--------------------------------------------------
end module data_types
