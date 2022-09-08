! File containing the procedures and modules used for the basic functions of reading input and printing output
!
! Record of revisions:
! DATE       PROGRAMMER        DESCRIPTION OF CHANGE
! ====       ==========        =====================
! 11/13/02   L. Pesce (UC)     Creation
! 01/05/04   L. Pesce (UC)     added new line character
! 06/06/2007 L. Pesce (UC)     changed line_lenght and line_format
! 10/24/07   L. Pesce (UC)     included file into the SYS_LIBS group and made it less proproc dependent
! 06/20/08   L. Pesce (UC)     included a unit for the input file
module io
  use data_types
  implicit none

  private

  public line_length ! length of character strings for the list input
  public line_format ! length of character strings formats for I/O for list input
  public print_fatal_error_line ! send to unit err_file a fatal error message
  public print_warning_line ! send to unit war_file a warning about calculation in progress
  public print_log_line ! send to unit log_file information about the evolving calculation
  public print_res_line ! send to unit res_file, the results of the calculation
  public print_plot_line ! send a line to the unit where plotting data are
  public print_MLE_par_line ! send a line to the unit where the MLE raw
  public print_AUC_SE_AUC_line ! send a line to the unit where the MLE raw
  public print_score_to_latent_line ! send a line to the unit where empirical variables are connected to latent ones
  public print_score_to_FPF_TPF_line ! send a line to the unit where the MLE raw parameters are reported
  public log_file
  public input_file_default

  integer, parameter :: line_length = 120 ! 80 ! length of character types
  character(len=6), parameter :: line_format = "(120a)" ! "(80a)"  Note that the value should be the same as in
                       ! the string line_length, this is the format for one line of characters

  ! NOTE, IMPORTANT: The unit numbers are defined with the array used_file_units, in this way it is much easier to
  !                  check, by the program, whether a unit is already allocated for some purpose in order to avoid
  !                  conflict. If a new unit is added, its unit number should be added to used_file_units, and then a
  !                  new mnemonic should be added.
  integer, parameter, dimension (9) :: used_file_units =(/ 8, 9, 6,12,13,15,16, 18, 19 /)
  integer, parameter:: log_file = used_file_units(7) ! send the log info straight to file standard IO
  integer, parameter:: err_file = used_file_units(7) ! send the errors straight to error IO, now standard I/O
  integer, parameter:: war_file = used_file_units(7) ! send the warning straight to error IO, now standard I/O
  integer, parameter:: res_file = used_file_units(9) ! results file is the 9th of the used_file units, or 19
  integer, parameter:: plot_file = used_file_units(1) ! data for plotting pusposes
  integer, parameter:: MLE_par_file  = used_file_units(8) ! File contains all the parameters estimate by the MLE
  integer, parameter:: AUC_SE_AUC_file  = used_file_units(6) ! contains only the estimated AUC and it SE
  integer, parameter:: score_to_latent_file  = used_file_units(4) ! contains information to to link test result value space to latent
  integer, parameter:: score_to_FPF_TPF_file  = used_file_units(5) ! contains information to to link test result value space to latent
  integer, parameter:: input_file_default = used_file_units(2) ! default Unit for input file

  contains

  !---------------------------------------
  subroutine print_score_to_FPF_TPF_line(err_msg)
  ! Prints a line of text to the file that contains info between score (test result values)
  ! and FPF TPF values
  implicit none
  character(len=line_length), intent(in):: err_msg

      continue! write(score_to_FPF_TPF_file,line_format ) err_msg

  end subroutine print_score_to_FPF_TPF_line
  !---------------------------------------

  !---------------------------------------
  subroutine print_score_to_latent_line(err_msg)
  ! Prints a line of text to the file that contains info between score (test result values)
  ! and latent spaces
  implicit none
  character(len=line_length), intent(in):: err_msg

      continue! write(score_to_latent_file,line_format ) err_msg

  end subroutine print_score_to_latent_line
  !---------------------------------------


  !---------------------------------------
  subroutine print_AUC_SE_AUC_line(err_msg)
  ! Prints a line of text to the file that contains AUC and its Se
  implicit none
  character(len=line_length), intent(in):: err_msg

      continue! write(AUC_SE_AUC_file,line_format ) err_msg

  end subroutine print_AUC_SE_AUC_line
  !---------------------------------------

  !---------------------------------------
  subroutine print_MLE_par_line(err_msg)
  ! Prints a line of text into the file containing the parameters estimated by the MLE
  implicit none
  character(len=line_length), intent(in):: err_msg

      continue! write(MLE_par_file,line_format ) err_msg

  end subroutine print_MLE_par_line
  !---------------------------------------

  !---------------------------------------
  subroutine print_fatal_error_line(err_msg)
  ! Subroutine which will send fatal error messages to whatever the error
  ! destination will be chosen to be, currently unit err_file
  implicit none
  character(len=line_length), intent(in):: err_msg

      continue! write(err_file,line_format ) err_msg

  end subroutine print_fatal_error_line
  !---------------------------------------

  !---------------------------------------
  subroutine print_warning_line(err_msg)
  ! Subroutine which will send non fatal warnings to whatever the warning
  ! destination will be chosen to be, currently unit war_file
  implicit none
  character(len=line_length), intent(in):: err_msg

      continue! write(war_file,line_format) err_msg

  end subroutine print_warning_line

  !-----------------------------------------------------

 !---------------------------------------
  subroutine print_log_line(msg)
  ! Subroutine which will send non informational lines to whatever the warning
  ! destination will be chosen to be, currently unit log_file
  implicit none
  character(len=line_length), intent(in):: msg

      continue! write(log_file,line_format) msg

  end subroutine print_log_line
  !-----------------------------------------------------

 !---------------------------------------
  subroutine print_plot_line(msg)
  ! Subroutine which will send plotting related data to whatever
  ! destination will be chosen to be, currently unit log_file
  implicit none
  character(len=line_length), intent(in):: msg

      continue! write(plot_file,line_format) msg

  end subroutine print_plot_line
  !-----------------------------------------------------

 !---------------------------------------
  subroutine print_res_line(msg)
  ! Subroutine which will send non informational lines to whatever the results
  ! destination will be chosen to be, currently file res_file
  implicit none
  character(len=line_length), intent(in):: msg

      continue! write(res_file,line_format) msg

  end subroutine print_res_line
  !-----------------------------------------------------



end module io
