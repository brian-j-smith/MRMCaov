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
  public connect_files ! function called to connect the units used for IO
  public connect_input_file ! function called to connect the input unit
  public print_fatal_error_line ! send to unit err_file a fatal error message
  public print_warning_line ! send to unit war_file a warning about calculation in progress
  public print_log_line ! send to unit log_file information about the evolving calculation
  public print_res_line ! send to unit res_file, the results of the calculation
  public print_plot_line ! send a line to the unit where plotting data are
  public print_MLE_par_line ! send a line to the unit where the MLE raw
  public print_AUC_SE_AUC_line ! send a line to the unit where the MLE raw
  public print_score_to_latent_line ! send a line to the unit where empirical variables are connected to latent ones
  public print_score_to_FPF_TPF_line ! send a line to the unit where the MLE raw parameters are reported
  public get_input_line !  read a line of input of size line_length
  public rewind_input_file ! reposition input file at the beginning
  public log_file
  public input_file_default

  integer, parameter :: line_length = 120 ! 80 ! length of character types
  character(len=6), parameter :: line_format = "(120a)" ! "(80a)"  Note that the value should be the same as in
                       ! the string line_length, this is the format for one line of characters 
  
  ! NOTE: log_file == err_file == war_file for connect_files() to work. If one desires these
  ! to be different, also the subroutine must be changed
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
  integer:: input_file ! Unit for input file DEFINED BY THE CONNECT_INPUT_FILE PROCEDURE
 
  contains

  !-----------------------------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------------------------
  subroutine rewind_input_file() 
  !PURPOSE:  reposition input file at the beginning
  !-----------------------------------------------------------------------------------------------------------------------------

    rewind(input_file)

  !-----------------------------------------------------------------------------------------------------------------------------
  end subroutine rewind_input_file
  !-----------------------------------------------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------------------------
  subroutine get_input_line(input_line,  stat) 
  !PURPOSE:  attempt to read a line from the previously connected unit input_file, return an error if this is not possible
  !-----------------------------------------------------------------------------------------------------------------------------

  character (len = line_length),intent(out) :: input_line ! the content of the input line read, if no error 
  integer, intent(OUT):: stat ! whether the reading was successful and if not what went wrong

  read(input_file,line_format,iostat = stat) input_line


  !-----------------------------------------------------------------------------------------------------------------------------
  end subroutine get_input_line
  !-----------------------------------------------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------------------------
  subroutine connect_input_file(input_unit_in, infile,stat) 
  !-----------------------------------------------------------------------------------------------------------------------------
  !PURPOSE:  To connect the input file that will be read by the input reading routines to the specified unit (input_unit_in)
  !          this unit will me mapped to input_file, which is the unit from which all the other procedures read. 
  !NOTE:     Often the unit will be the input_file_default unit, as defined by the module parameters above.
  !NOTE   :  It did not seem to make sense to allow for multiple input files, so I allowed only one. However, the routine can
  !          sequentially connect different files to different units. 

  implicit none

  integer, intent(IN), optional:: input_unit_in ! whether the opening of the file was successful, i.e. the file was there and it was
  character(len=line_length),optional, intent(INOUT) :: infile ! character string containing the input file name 
  integer, intent(OUT):: stat ! whether the opening of the file was successful, i.e. the file was there and it was
                              ! accessible to read 

  if (.not.present(input_unit_in) .and. .not.present(infile)) then ! The unit is not specified and also the filename is not specified, therefore
          ! we assume that the input will be done using Standard input (5). There is no need to connect it
          input_file = 5   
  elseif (present(input_unit_in) .and. .not.present(infile)) then ! The unit is specified but the filename is not specified, it must be a 
         ! progamming mistake
         stat = 1
  else ! We have at least an input file, so we need to connect it
          if (.not.present(input_unit_in)) then ! The unit is not specified, therefore we use the default input unit defined in the module
                   input_file = input_file_default ! define the global variable that will read the input files      
          else 
               ! First verify that the unit is not already in use or could be used creating a conflict
               if( any(used_file_units == input_unit_in) ) then
                   stat = 1
                   return
               else
                 input_file = input_unit_in ! define the global variable that will read the input files       
               endif
          endif   
          ! Conect the units
          open(unit=input_file,file=trim(infile),IOSTAT=stat,STATUS='old', form = 'formatted')
          if(stat == 0) then
            inquire(unit=input_file,NAME=infile) ! Verify that the connected name is consistent with the 
                                               ! file name -- sometimes there are issues with directories and folder
                                               ! names being attached to it
            infile=trim(infile) ! Return the name of the input file
          endif
  endif

  end subroutine connect_input_file 
  !---------------------------------------

  !----------------------------------------------------------------------------------
  subroutine connect_files(which_files)  
  !PURPOSE:  to connect the output files.
  !NOTE: This works only if log_file == err_file == war_file (see above)
  !      since in fortran one cannot connect multiple units to the same
  !      file.
  !NOTE: July 2008, UC, LP added the integer flag call to allow the opening of different files
  !               currently 1 is for all files, 0 is only for it.out  
  implicit none

  integer, intent(IN):: which_files ! which files to open, flag.

  open (unit = log_file, file = 'it.out', form = 'formatted')
  open (unit = res_file, file = 'MLE_fit.out', form = 'formatted')

   if(which_files == 1) then
       open (unit = plot_file, file = 'roc_plot.out', form = 'formatted')
       open (unit = MLE_par_file, file = 'MLE_par.out', form = 'formatted')
       open (unit = AUC_SE_AUC_file, file = 'AUC_SE_AUC.out', form = 'formatted')
       open (unit = score_to_latent_file, file = 'Score_to_latent.out', form = 'formatted')
       open (unit = score_to_FPF_TPF_file, file = 'Score_to_FPF_TPF.out', form = 'formatted')
   endif

  end subroutine connect_files 
  !------------------------------------------------------------------------------------------------------

  !---------------------------------------
  subroutine print_score_to_FPF_TPF_line(err_msg)
  ! Prints a line of text to the file that contains info between score (test result values)
  ! and FPF TPF values
  implicit none
  character(len=line_length), intent(in):: err_msg
  
      write(score_to_FPF_TPF_file,line_format ) err_msg
  
  end subroutine print_score_to_FPF_TPF_line 
  !---------------------------------------

  !---------------------------------------
  subroutine print_score_to_latent_line(err_msg)
  ! Prints a line of text to the file that contains info between score (test result values)
  ! and latent spaces
  implicit none
  character(len=line_length), intent(in):: err_msg
  
      write(score_to_latent_file,line_format ) err_msg
  
  end subroutine print_score_to_latent_line 
  !---------------------------------------


  !---------------------------------------
  subroutine print_AUC_SE_AUC_line(err_msg)
  ! Prints a line of text to the file that contains AUC and its Se 
  implicit none
  character(len=line_length), intent(in):: err_msg
  
      write(AUC_SE_AUC_file,line_format ) err_msg
  
  end subroutine print_AUC_SE_AUC_line 
  !---------------------------------------

  !---------------------------------------
  subroutine print_MLE_par_line(err_msg)
  ! Prints a line of text into the file containing the parameters estimated by the MLE
  implicit none
  character(len=line_length), intent(in):: err_msg
  
      write(MLE_par_file,line_format ) err_msg
  
  end subroutine print_MLE_par_line 
  !---------------------------------------
 
  !---------------------------------------
  subroutine print_fatal_error_line(err_msg)
  ! Subroutine which will send fatal error messages to whatever the error
  ! destination will be chosen to be, currently unit err_file
  implicit none
  character(len=line_length), intent(in):: err_msg
  
      write(err_file,line_format ) err_msg
  
  end subroutine print_fatal_error_line
  !---------------------------------------
  
  !---------------------------------------
  subroutine print_warning_line(err_msg)
  ! Subroutine which will send non fatal warnings to whatever the warning
  ! destination will be chosen to be, currently unit war_file
  implicit none
  character(len=line_length), intent(in):: err_msg
  
      write(war_file,line_format) err_msg
  
  end subroutine print_warning_line

  !-----------------------------------------------------
 
 !---------------------------------------
  subroutine print_log_line(msg)
  ! Subroutine which will send non informational lines to whatever the warning
  ! destination will be chosen to be, currently unit log_file
  implicit none
  character(len=line_length), intent(in):: msg
  
      write(log_file,line_format) msg
  
  end subroutine print_log_line
  !-----------------------------------------------------

 !---------------------------------------
  subroutine print_plot_line(msg)
  ! Subroutine which will send plotting related data to whatever
  ! destination will be chosen to be, currently unit log_file
  implicit none
  character(len=line_length), intent(in):: msg
  
      write(plot_file,line_format) msg
  
  end subroutine print_plot_line
  !-----------------------------------------------------

 !---------------------------------------
  subroutine print_res_line(msg)
  ! Subroutine which will send non informational lines to whatever the results
  ! destination will be chosen to be, currently file res_file
  implicit none
  character(len=line_length), intent(in):: msg
  
      write(res_file,line_format) msg
  
  end subroutine print_res_line
  !-----------------------------------------------------


  
end module io

