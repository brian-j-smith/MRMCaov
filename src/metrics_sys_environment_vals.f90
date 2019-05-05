! Contains functions and variables that allow to inquire about the computational environment
! example would be inquiry about the time and the date when the program was started
! STARTED BY Lorenzo Pesce in october 2007 for the restauration of the ROCKIT project

module environment_vals

private !Default is that data and procedures are not available

public the_date ! character variable that contains date data 
public the_time ! charaster variable containing time data

character (len=9):: the_date ! stores the date at which a calculation was done
character (len=8):: the_time ! stores the time at which a calculation was done

contains ! procedures defined by the module

subroutine mydateandtime(current_date,current_time) 
! PURPOSE:  return 2 strings, one with the date in format yy--mm--dd
!                            and one with the time in format hh:mm:ss

!ALGORITHM: makes use of internal function two_digits to convert integers bigger 
!           than 9
! NOTE:     Imported nearly identically from LP version of ROCKIT, "fixed" in 2004

 character (len=9), intent(OUT):: current_date
 character (len=8), intent(OUT):: current_time

 integer, dimension(8):: values ! array that stores the values from the intrinsic function date_and_time
                                ! used to convert numeric months into strigs
 character  (len = 3), parameter,dimension(12)::&
      months =(/"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",&
               "Sep", "Oct", "Nov", "Dec"/)

 intrinsic date_and_time
 
 call date_and_time( values = values)

 current_date = two_digits(values(3)) // "-" // months(values(2)) // "-" //&
            two_digits(values(1) - 2000)
 current_time = two_digits(values(5))//":"//two_digits(values(6)) &
            &//":"//two_digits(values(7))

 contains
 elemental function two_digits(the_number)
 
   character(len = 2) two_digits
   integer,intent(IN):: the_number
   integer:: dec,unit ! split the numbers in first and second digit

   dec = the_number/10
   unit= the_number - dec*10

   two_digits = achar(48+dec) // achar(48+unit)
 end function two_digits

 end subroutine mydateandtime


end module environment_vals


