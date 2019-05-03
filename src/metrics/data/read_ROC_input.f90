!  module containing the procedures used for reading input that was written in list form
!  as per june 2005. the description of the list input form is right at the beginning of the contains section
!  of the module. basically it reads an input file where the cases are input one at a time, and it does not
!  matter whether they are integer or real numbers. first it reads the actually negative cases, then the
!  actually positive ones. note that is takes both the american "." and the european "," as separators 
!  between the integer part and the decimals. it reads both e format and d formats. 
! note:  this module still needed a lot of polishing, but it should work 
!        the error flag reporting is working fine, but it could use some polishing
! note2: the algorithm wasn't redesigned in f95. Accordingly is uglier than it need to be, especially the
!        handling of strings. A lot of the operations done on buffers aren't necessary in fortran 90, 
!        where we can seek directly the location of a specific character 
!note3:  all routines use the same unit to read the input. If multiple files need to be called, the
!        file unit needs to be closed in between

module read_ROC_input
 use data_types   
 use io, only: line_length, line_format, get_input_line! size of the input line, and input file unit
 use io, only: rewind_input_file 


 implicit none
  
  ! character constants used for comparison in the files
 character, parameter :: iblank = " " 
 character, parameter :: inewline = char(13) ! beginning of new line, technically this is an extension of the standard, but in this 
                    ! the standard is really stupid, so I don't mind.
 character, parameter :: istar = '*'
 character(len=*), parameter :: small = 'sS' ! Characters tell positivity is for small values
 character(len=*), parameter :: large = 'lL' ! Characters tell positivity is for large values
 character, dimension(line_length) :: lstring
 character, dimension(line_length) :: line
 integer :: length ! lenght of token currently being considered
 

! Return flag codes for the kind of file that the input file seems to be
 integer, parameter:: labmrmc_file    =  1 ! code for input file being labmrmc format
 integer, parameter:: rockit_file     =  2 ! code for input file being rockit format

 integer, parameter:: unknown    =  0 ! ERROR: the input file is of unknown format
 integer, parameter:: unreadable = -1 ! ERROR: The input file can't be read
 integer, parameter:: empty      = -2 ! ERROR: The input file seems to be empty
 integer, parameter:: incomplete = -3 ! ERROR: The input file seems to be empty
 
  private
! OLD and NEW mean that there are currently 2 versions of the same procedure, of which one is
!             legacy and will eventually be eliminated.
!  procedure    operation
! =========    ==========
 public:: read_labmrmc_file ! read the labmrmc file type
 public:: parse_labmrmc_file ! parses the labmrmc input file, checks what format it is and returns
                          ! data file characteristic
 public:: parse_rockit_file ! parses the ROCKIT input file and data file characteristic
 public:: read_rockit_file ! read the ROCKIT input file and data file characteristic
 public::  parse_labroc_list_data ! oversees labroc single file list input 
 public::  read_labroc_list_data ! oversees labroc single file list input 
! Variables and values
 public:: labmrmc_file, unknown, unreadable, empty, incomplete ! return codes for input file format types - see above 
 private:: get_title! read title line 
 private:: get_reader! read reader name
 private:: get_positivity ! determine what the positivity flags are
 private:: get_modality_names! read title line 
 private:: get_word! NEW: extracts a word (as a sequence of non blank characters separated by blanks or ends)
 ! Older labmrmc procedures that need some clean up
 private:: read_header! read title line and file description for labroc files 
 private:: readin! read numerical input data part 
 private:: tonum ! convert a strig made of digits into a real number, here are contained the details of 
                 ! the numerical formats that we are willing to accept
 
contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine read_rockit_file(num_modalities, num_normal_cases, num_abnormal_cases, &
                            case_data, design_matrix, dataset_title, modality_names, positiveislarge,&
                            ierror, err_msg)
!--------------------------------------------------------------------------------------------
!PURPOSE: read and store the data from  a ROCKIT input file.
! NOTE:   the procedur won't count missing data points or create a design matrix. This file simply 
!         conunts the total number of cases. The number of empty locations will be considered in the 
!         design matrix in read_rockit_file

integer, intent(in):: num_modalities ! number of modalities
integer, intent(in):: num_normal_cases 
integer, intent(in):: num_abnormal_cases 

character(len = line_length), intent(out) :: dataset_title
character(len = line_length), dimension (num_modalities), intent(out) :: modality_names 
              ! this is why it is outrageously large
integer, dimension(num_modalities), intent(OUT) :: positiveislarge ! where is positivity for the different modalities (1 for large, 0 for small)
integer, intent(OUT):: ierror
character(len=line_length), intent(out) :: err_msg ! error message from subroutines
real(kind=double), dimension(num_modalities,num_normal_cases + num_abnormal_cases) :: case_data ! data by case, by modality
integer, dimension(num_modalities,num_normal_cases + num_abnormal_cases) :: design_matrix ! data by case, by modality, whether
                  ! value is present

!! Work variables
integer:: num_mod


integer:: neg_count, pos_count ! the number of positive and negative cases found for this reader

! Set input result to default value, that is the file is correctly a ROCKIT file
ierror = 0

!!! READ AND CHECK TITLE LINE
call get_title(dataset_title,ierror, err_msg)
! Check if the title is as expected, the error flags have many values, meaning different 
! errors, currently all problems are considered irreparable, but it doesn't need to be that
! way
if ( ierror < 0 .or. ierror == 1) then ! if the reading of the first line failed, this file is not usable
   write(err_msg,*) "First line (title line) was empty or unreadable unknown kind of input file"
   return
endif


!! READ NAMES AND CHECK THE NUMBER OF MODALITIES AND READ THEIR NAMES
num_mod = num_modalities
call get_modality_names(modality_names, num_mod, ierror, err_msg)

! Check the content of the modality names slot is as expected,
if ( ierror < 0 .or. ierror == 1 .or. num_mod .ne. num_modalities) then 
  ierror = incomplete
  write(err_msg,*) " the modality names could not be read: check file"
  return
endif

! Determine where is positivity for the different modalities

call get_positivity(num_modalities, positiveislarge, ierror, err_msg)

if(ierror /= 0) then
   return
endif

!!! READ CASES !!!!!!!!!!!
! Count the number of negative cases and  check whether the input format is correct
neg_count = num_normal_cases

call get_cases_miss(num_modalities, neg_count, case_data, design_matrix, ierror, err_msg)

! There was a mistake in the cases input
if(ierror == 1) then
     err_msg = "for Neg cases " // err_msg(1:59) 
     return  
endif

! Count the number of positive cases and  check whether the input format is correct
pos_count = num_abnormal_cases

! Pass only the upper part of the array
call get_cases_miss(num_modalities, neg_count, case_data(:,num_normal_cases+1:), design_matrix(:,num_normal_cases+1:),&
                     ierror, err_msg)

!If there was a mistake in the cases input
if(ierror == 1) then
     err_msg = "for Pos cases " // err_msg(1:59) 
     return  
endif

err_msg = ""

!--------------------------------------------------------------------------------------------
end subroutine read_rockit_file
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine parse_rockit_file(file_type, num_modalities, num_normal_cases, num_abnormal_cases, err_msg)
!--------------------------------------------------------------------------------------------
!PURPOSE: read and parse a ROCKIT input file to return the TOTAL number of cases and modalities.
!         makes sure that the file in fact has the correct format.
! NOTE:   the procedur won't count missing data points or create a design matrix. This file simply 
!         conunts the total number of cases. The number of empty locations will be considered in the 
!         design matrix in read_rockit_file
integer, intent(out):: file_type !

integer, intent(out):: num_modalities ! number of modalities
integer, intent(out):: num_normal_cases 
integer, intent(out):: num_abnormal_cases 
character(len=line_length), intent(out) :: err_msg ! error message from subroutines

character(len = line_length) :: dataset_title
character(len = line_length), dimension (100) :: modality_names ! Scratch file where to put the names for now
              ! this is why it is outrageously large

integer:: ierror

integer, dimension(:), allocatable :: positiveislarge ! where is positivity for the different modalities
real(kind=double), dimension(:,:), allocatable :: case_data ! data by case, by modality
integer, dimension(:,:), allocatable :: design_matrix ! data by case, by modality: whether a case is present

integer:: neg_count, pos_count ! the number of positive and negative cases found for this reader


! Set input result to default value, that is the file is correctly a ROCKIT file
file_type =  rockit_file

!!! READ AND CHECK TITLE LINE


call get_title(dataset_title,ierror, err_msg)


! Check if the title is as expected, the error flags have many values, meaning different 
! errors, currently all problems are considered irreparable, but it doesn't need to be that
! way
if ( ierror < 0) then ! if the reading of the first line failed, this file is not usable
  file_type = ierror
  return
elseif(ierror == 1) then ! title line is empty, probably it isn't ROCKIT format
   file_type = unknown
   write(err_msg,*) "First line (title line) was empty, unknown kind of input file"
   return
endif

!! DETERMINE THE NUMBER OF MODALITIES AND READ THEIR NAMES
num_modalities = 100 ! Initialize the array size to A ludicrously large value for now


call get_modality_names(modality_names, num_modalities, ierror, err_msg)

! Check the content of the modality names slot is as expected,
if ( ierror < 0 .or. ierror == 1) then ! if the reading of the modality names line did not work out.
  file_type = incomplete
  write(err_msg,*) " the modality names could not be read: check file"
  return
endif

! Determine where is positivity for the different modalities
allocate(positiveislarge(num_modalities)) 

call get_positivity(num_modalities, positiveislarge, ierror, err_msg)

if(ierror /= 0) then
   file_type = incomplete
   return
endif


deallocate(positiveislarge ) 

!!! READ CASES !!!!!!!!!!!
! Count the number of negative cases and  check whether the input format is correct
num_normal_cases = 3
neg_count = num_normal_cases
allocate(case_data(num_modalities, neg_count ))
allocate(design_matrix(num_modalities, neg_count ))

call get_cases_miss(num_modalities, neg_count, case_data, design_matrix, ierror, err_msg)

deallocate(case_data, design_matrix)

! There was a mistake in the cases input
if(ierror == 1) then
     file_type = incomplete
     err_msg = "for Neg cases " // err_msg(1:59) 
     return  
! Now we know the number of cases
elseif(ierror == 2) then
    num_normal_cases = neg_count
endif


! Count the number of positive cases and  check whether the input format is correct
num_abnormal_cases = 3
pos_count = num_abnormal_cases
allocate(case_data(num_modalities, pos_count ))
allocate(design_matrix(num_modalities, pos_count ))

call get_cases_miss(num_modalities, pos_count, case_data, design_matrix, ierror, err_msg)


deallocate(case_data, design_matrix)
!If there was a mistake in the cases input
if(ierror == 1) then
     file_type = incomplete
     err_msg = "for Pos cases " // err_msg(1:59) 
     return  
! The number of cases wasn't 3 therefore we update the value
elseif(ierror == 2) then
     num_abnormal_cases = pos_count
endif


err_msg = ""

call rewind_input_file()

!--------------------------------------------------------------------------------------------
end subroutine parse_rockit_file
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------



!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine read_labmrmc_file(num_rdr, num_modalities, num_normal_cases, num_abnormal_cases, &
                         input_data, dataset_title, reader_names, modality_names, positiveislarge, &
                         ierror, err_msg)
!--------------------------------------------------------------------------------------------
!PURPOSE: read and parse a labmrmc input file to return the number of readers, cases and modalities.
!         makes sure that the file in fact has the correct format.

integer, intent(in):: num_rdr ! number of readers in the input file
integer, intent(in):: num_modalities ! number of modalities
integer, intent(in):: num_normal_cases 
integer, intent(in):: num_abnormal_cases 
real(kind=double), intent(out), dimension(num_modalities, num_rdr, num_normal_cases + num_abnormal_cases):: input_data
integer, intent(out):: ierror
character(len=line_length), intent(out) :: err_msg ! error message from subroutines

character(len = line_length), intent(out) :: dataset_title
character(len = line_length), intent(out), dimension(num_rdr) :: reader_names
character(len = line_length), intent(out),  dimension (num_modalities) :: modality_names  !the names for modalities
integer, intent(out), dimension(num_modalities) :: positiveislarge ! where is positivity for the different modalities

character(len=line_length) :: msg ! message
character(len=line_length) :: scratch !
integer:: neg_count, pos_count ! the number of positive and negative cases found for this reader
integer:: num_mod, irdr

!!! READ AND CHECK TITLE LINE
call get_title(dataset_title,ierror, err_msg)
! Check if the title is as expected, the error flags have many values, meaning different 
! errors, currently all problems are considered irreparable, but it doesn't need to be that
! way
if ( ierror < 0) then ! if the reading of the first line failed, this file is not usable
  return
elseif(ierror == 1) then ! title line is empty, probably it isn't LABMRMC format
   ierror = unknown
   write(err_msg,*) "First line (title line) was empty, unknown kind of input file"
   return
endif

!!! READ FIRST READERS NAME AND CHECK IF FILE IS OK
call get_reader(reader_names(1), ierror, err_msg)
! Check the content of the first reader name slot is as expected
if ( ierror < 0 .or. ierror > 1) then ! if the reading of the reader line did not work out.
  ierror = incomplete
  write(err_msg,*) " the first reader name could not be read: check file"
  return
endif

!! READ MODALITY NAMES
num_mod = num_modalities
call get_modality_names(modality_names, num_mod, ierror, err_msg)

! Check the content of the modality names slot is as expected,
if ( ierror < 0 .or. ierror == 1 .or. num_mod .ne. num_modalities) then 
  ierror = incomplete
  write(err_msg,*) " the modality names could not be read: check file"
  return
endif


call get_positivity(num_modalities, positiveislarge, ierror, err_msg)

if(ierror /= 0) then
   ierror = incomplete
   return
endif


!!! READ CASES FOR FIRST READER !!!!!!!!!!!
! Start with the first reader
! Count the number of negative cases and  check whether the input format is correct
neg_count = num_normal_cases

call get_cases(num_modalities, neg_count, input_data(:,1,1:num_normal_cases), ierror, err_msg)

! There was a mistake in the cases input
if(ierror == 1) then
     ierror = incomplete
     err_msg = "for Rdr 1, Neg cases " // err_msg(1:59) 
     return  
endif

! Count the number of positive cases and  check whether the input format is correct
pos_count = num_abnormal_cases

call get_cases(num_modalities, pos_count, input_data(:,1,num_normal_cases+1:), ierror, err_msg)

 ! There was a mistake in the cases input
if(ierror == 1) then
     ierror = incomplete
     err_msg = "for Rdr 1, Pos cases " // err_msg(1:59) 
     return  
endif

err_msg = ""

!!! LOOP OVER THE READERS DATA, CHECKING WHETHER THEY HAVE ALL THE SAME NUMBER
!!! OF CASES AND COUNTING THEM AS WELL.
!!! WHILE THE LOOP IS THEORETICALLY INFINITE, IT SHOULD BE NOTICED THAT AN
!!! INFINITE LOOP IS POSSIBLE ONLY IF THE INPUT FILE IS INFINITE
do irdr = 2, num_rdr 
    call get_reader(reader_names(irdr), ierror, err_msg)
    ! Check the content of the first reader name slot is as expected
     if ( ierror < 0 .or. ierror ==  1) then ! if the reading of the reader line did not work out.
           ierror = incomplete
           write(err_msg,*) "the reader ", num_rdr," name could not be read: check file"
           return
    endif
    ! READ AND PARSE THE CASES ASSOCIATED WITH THIS READER TO VERIFY WHETHER THEIR ARE
    ! CORRECT AND IN THE RIGHT NUMBER 
    !Check whether there was a mistake in the NEGATIVE cases input.
    neg_count = num_normal_cases
    call get_cases(num_modalities, neg_count, input_data(:,irdr,1:num_normal_cases), ierror, err_msg)
    if(ierror /= 0) then
        ierror = incomplete
        msg = err_msg
        write(err_msg,*) "for Rdr ", num_rdr, " Neg cases ", msg(1:55) 
        return  
    endif
    !Check whether there was a mistake in the POSITIVE cases input.
    pos_count = num_abnormal_cases
    call get_cases(num_modalities, pos_count, input_data(:, irdr,num_normal_cases+1:), ierror, err_msg)
    if(ierror /= 0) then
        ierror = incomplete
        msg = err_msg
        write(err_msg,*) "for Rdr ", num_rdr, " Pos cases ", msg(1:55) 
        return  
    endif
enddo

! Double check if the file is finished
call get_reader(scratch, ierror, err_msg)
if(ierror == 2) then ! the input file is finished
  ierror = 1
  err_msg = "missing end of file at second pass "
endif

!--------------------------------------------------------------------------------------------
end subroutine read_labmrmc_file
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine parse_labmrmc_file(file_type, num_rdr, num_modalities, num_normal_cases, num_abnormal_cases, err_msg)
!--------------------------------------------------------------------------------------------
!PURPOSE: read and parse a labmrmc input file to return the number of readers, cases and modalities.
!         makes sure that the file in fact has the correct format.

integer, intent(out):: file_type !

integer, intent(out):: num_rdr ! number of readers in the input file
integer, intent(out):: num_modalities ! number of modalities
integer, intent(out):: num_normal_cases 
integer, intent(out):: num_abnormal_cases 
character(len=line_length), intent(out) :: err_msg ! error message from subroutines

character(len = line_length) :: dataset_title
character(len = line_length) :: reader_name
character(len = line_length), dimension (100) :: modality_names ! Scratch file where to put the names for now
              ! this is why it is outrageously large
integer:: ierror

integer, dimension(:), allocatable :: positiveislarge ! where is positivity for the different modalities
real(kind=double), dimension(:,:), allocatable :: case_data ! data by case, by modality
character(len=line_length) :: msg ! message

integer:: neg_count, pos_count ! the number of positive and negative cases found for this reader

! Set input result to default value, that is the file is correctly a labmrmc file
file_type =  labmrmc_file

!!! READ AND CHECK TITLE LINE
call get_title(dataset_title,ierror, err_msg)
! Check if the title is as expected, the error flags have many values, meaning different 
! errors, currently all problems are considered irreparable, but it doesn't need to be that
! way
if ( ierror < 0) then ! if the reading of the first line failed, this file is not usable
  file_type = ierror
  return
elseif(ierror == 1) then ! title line is empty, probably it isn't LABMRMC format
   file_type = unknown
   write(err_msg,*) "First line (title line) was empty, unknown kind of input file"
   return
endif

!!! READ FIRST READERS NAME AND CHECK IF FILE IS OK
call get_reader(reader_name, ierror, err_msg)
! Check the content of the first reader name slot is as expected
if ( ierror < 0 .or. ierror > 1) then ! if the reading of the reader line did not work out.
  file_type = incomplete
  write(err_msg,*) " the first reader name could not be read: check file"
  return
else
   num_rdr = 1
endif

!! DETERMINE THE NUMBER OF MODALITIES AND READ THEIR NAMES
num_modalities = 100 ! Initialize the array size to A ludicrously large value for now

call get_modality_names(modality_names, num_modalities, ierror, err_msg)

! Check the content of the modality names slot is as expected,
if ( ierror < 0 .or. ierror == 1) then ! if the reading of the modality names line did not work out.
  file_type = incomplete
  write(err_msg,*) " the modality names could not be read: check file"
  return
endif

! Determine where is positivity for the different modalities
allocate(positiveislarge(num_modalities)) 

call get_positivity(num_modalities, positiveislarge, ierror, err_msg)

if(ierror /= 0) then
   file_type = incomplete
   return
endif

deallocate(positiveislarge ) 

!!! READ CASES FOR FIRST READER !!!!!!!!!!!
! Start with the first reader
! Count the number of negative cases and  check whether the input format is correct
num_normal_cases = 10
neg_count = num_normal_cases
allocate(case_data(num_modalities, neg_count ))

call get_cases(num_modalities, neg_count, case_data, ierror, err_msg)

deallocate(case_data)

! There was a mistake in the cases input
if(ierror == 1) then
     file_type = incomplete
     err_msg = "for Rdr 1, Neg cases " // err_msg(1:59) 
     return  
! Now we know the number of cases
elseif(ierror == 2) then
    num_normal_cases = neg_count
endif

! Count the number of positive cases and  check whether the input format is correct
num_abnormal_cases = 10
pos_count = num_abnormal_cases
allocate(case_data(num_modalities, pos_count ))

call get_cases(num_modalities, pos_count, case_data, ierror, err_msg)
deallocate(case_data)
 ! There was a mistake in the cases input
if(ierror == 1) then
     file_type = incomplete
     err_msg = "for Rdr 1, Pos cases " // err_msg(1:59) 
     return  
! The number of cases was wrong, let's retry
elseif(ierror == 2) then
     num_abnormal_cases = pos_count
endif

err_msg = ""

!!! LOOP OVER THE READERS DATA, CHECKING WHETHER THEY HAVE ALL THE SAME NUMBER
!!! OF CASES AND COUNTING THEM AS WELL.
!!! WHILE THE LOOP IS THEORETICALLY INFINITE, IT SHOULD BE NOTICED THAT AN
!!! INFINITE LOOP IS POSSIBLE ONLY IF THE INPUT FILE IS INFINITE
reader_count_loop: do 
    !!! READ  READER LINE TO FIND NAME NAME AND CHECK IF FILE IS FINISHED OR INCORRECT
    call get_reader(reader_name, ierror, err_msg)
    ! Check the content of the first reader name slot is as expected
     if ( ierror < 0 .or. ierror ==  1) then ! if the reading of the reader line did not work out.
           file_type = incomplete
           write(err_msg,*) "the reader ", num_rdr," name could not be read: check file"
           return
    elseif(ierror == 2) then ! the input file is finished
           ierror = 0
           exit reader_count_loop
    else
          num_rdr = num_rdr + 1  
    endif
    ! READ AND PARSE THE CASES ASSOCIATED WITH THIS READER TO VERIFY WHETHER THEIR ARE
    ! CORRECT AND IN THE RIGHT NUMBER 
    !Check whether there was a mistake in the NEGATIVE cases input.
    neg_count = num_normal_cases
    allocate(case_data(num_modalities, neg_count ))
    call get_cases(num_modalities, neg_count, case_data, ierror, err_msg)
    deallocate(case_data)
    if(ierror /= 0) then
        file_type = incomplete
        msg = err_msg
        write(err_msg,*) "for Rdr ", num_rdr, " Neg cases ", msg(1:55) 
        return  
    endif
    !Check whether there was a mistake in the POSITIVE cases input.
    pos_count = num_abnormal_cases
    allocate(case_data(num_modalities, pos_count ))
    call get_cases(num_modalities, pos_count, case_data, ierror, err_msg)
    deallocate(case_data)
    if(ierror /= 0) then
        file_type = incomplete
        msg = err_msg
        write(err_msg,*) "for Rdr ", num_rdr, " Pos cases ", msg(1:55) 
        return  
    endif
enddo reader_count_loop


!--------------------------------------------------------------------------------------------
end subroutine parse_labmrmc_file
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine get_cases_miss(num_modalities, num_cases, case_data, design_matrix, ierror, err_msg)
! PURPOSE: read the cases of one truth (positive or negative or other) for a group of
!          modalities. It assumes that the first line in the input file is where the
!          cases start. It also assumes that the number of cases will be finished when
!          the input will have one or more "*" on a line
!          What is after the num_modalities numbers is not read or used so it can be comment
!          Some values for cases might be missing and be indicated by a "#".
implicit none

integer, intent(in)  :: num_modalities 
integer, intent(inout)  :: num_cases ! On entry, the number of cases before this run, or 
                                     ! exit the number of cases actually found
real(kind=double), dimension(num_modalities, num_cases ):: case_data ! the cases data. Note
                           ! that it is incomplete if num_cases on entry is actually smaller
                           ! than the number found here.
integer, dimension(num_modalities, num_cases ):: design_matrix ! whether a value for a case
                                     ! is present or not for a modality
integer, intent(out) :: ierror        ! 0 ->  read successfully
                                      ! 1 ->  Case input missing/incorrect
                                      ! 2 -> num_cases different from entry value of num_cases
character(len=line_length), intent(out) :: err_msg ! error message from subroutines


! internal variables
integer:: istat
character(len=line_length)  :: input_line ! character buffer that contains an input line
character(len=1), dimension(line_length)  :: buffer ! character buffer that contains an input line
character(len=line_length)  :: word ! character buffer that contains an input line
integer:: icase, imod,i
real(kind=double):: case_value ! value for this case modality combination


ierror = 0
input_line = ""
icase = 0


cases_loop: do
    icase = icase + 1
    ! Read the input like as a character line

    call get_input_line(input_line,  istat) 
    ! eliminate all the spaces to the left of the flags
    input_line = trim(adjustl(input_line))
    ! check whether the reading worked correctly
    if (istat == - 1 .or. len_trim(input_line) == 0 ) then
         ierror = 1 ! the record is empty
         write(err_msg,*) "Case value line ", icase,  "is empty"  
         return
    elseif (istat >  0) then
         ierror = 1
         write(err_msg,*) "For case line ", icase,  "input can't be read or EOF"  
         return
    elseif ( verify(  input_line(1:1)  , '*') == 0  ) then
        icase = icase - 1 ! Last line wasn't a case
        exit cases_loop
    end if
    mod_loop: do imod = 1, num_modalities
          ! Extract the imod's word and check whether it worked
          call get_word(input_line, word, ierror)
          if (ierror == 1) then
               write(err_msg,*)  "Case line ", icase, " Modality ", imod, " has no input value "
               ierror = 1
               return
         endif
         ! Try to turn the extracted word into a number
         buffer = "" 
         do i = 1, len_trim(word) 
                buffer(i) = word(i:i)
         enddo
         ! Verify whether cases is missing data for this modality
         if ( verify(  word(1:1)  , '#') == 0  ) then
             if( icase < num_cases) then ! Load the case value only if we are not overshooting, if we are either the data
                 ! file is wrong or we are parsing the data
                 case_data(imod, icase) =  -666.0_double
                 design_matrix(imod, icase) =  0 ! case is missing
             endif
         else ! Either it is numerical data or it is an input error
             call tonum(len_trim(word), buffer, case_value, ierror)
             if (ierror == 1) then ! Could not read number, input file error
                write(err_msg,*)  "Case ", icase, " Modality ", imod, " value is wrong = ", trim(word)
                ierror = 1
                return
             endif
             if( icase <= num_cases) then ! Load the case value only if we are not overshooting, if we are either the data
                case_data(imod, icase) =  case_value !  case is missing
                design_matrix(imod, icase) =  1 !
             endif
         endif
    enddo   mod_loop

enddo cases_loop


if (icase /= num_cases ) then
  ierror  = 2
  write(err_msg,*) "Case_data has wrong size = ", num_cases, " ", icase
  num_cases = icase
endif


!--------------------------------------------------------------------------------------------
end subroutine get_cases_miss
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine get_cases(num_modalities, num_cases, case_data, ierror, err_msg)
! PURPOSE: read the cases of one truth (positive or negative or other) for a group of
!          modalities. It assumes that the first line in the input file is where the
!          cases start. It also assumes that the number of cases will be finished when
!          the input will have one or more "*" on a line
!          What is after the num_modalities numbers is not read or used so ti can be comment
implicit none

integer, intent(in)  :: num_modalities 
integer, intent(inout)  :: num_cases ! On entry, the number of cases before this run, or 
                                     ! exit the number of cases actually found
real(kind=double), dimension(num_modalities, num_cases ):: case_data ! the cases data. Note
                           ! that it is incomplete if num_cases on entry is actually smaller
                           ! than the number found here.
integer, intent(out) :: ierror        ! 0 ->  read successfully
                                      ! 1 ->  Case input missing/incorrect
                                      ! 2 -> num_cases different from entry value of num_cases
character(len=line_length), intent(out) :: err_msg ! error message from subroutines


! internal variables
integer:: istat
character(len=line_length)  :: input_line ! character buffer that contains an input line
character(len=1), dimension(line_length)  :: buffer ! character buffer that contains an input line
character(len=line_length)  :: word ! character buffer that contains an input line
integer:: icase, imod,i
real(kind=double):: case_value ! value for this case modality combination


ierror = 0
input_line = ""
icase = 0

cases_loop: do
    icase = icase + 1
    ! Read the input like as a character line

    call get_input_line(input_line,  istat) 
    ! eliminate all the spaces to the left of the flags
    input_line = trim(adjustl(input_line))
    ! check whether the reading worked correctly
    if (istat == - 1 .or. len_trim(input_line) == 0 ) then
         ierror = 1 ! the record is empty
         write(err_msg,*) "Case value line ", icase,  "is empty"  
         return
    elseif (istat >  0) then
         ierror = 1
         write(err_msg,*) "For case line ", icase,  "input can't be read or EOF"  
         return
    elseif ( verify(  input_line(1:1)  , '*') == 0  ) then
        icase = icase - 1 ! Last line wasn't a case
        exit cases_loop
    end if

    mod_loop: do imod = 1, num_modalities
          ! Extract the imod's word and check whether it worked
          call get_word(input_line, word, ierror)
          if (ierror == 1) then
               write(err_msg,*)  "Case line ", icase, " Modality ", imod, " has no input value "
               ierror = 1
               return
         endif
         ! Try to turn the extracted word into a number 
         do i = 1, line_length 
                buffer(i) = word(i:i)
         enddo

         call tonum(len_trim(word), buffer, case_value, ierror)

         if (ierror == 1) then
               write(err_msg,*)  "Case ", icase, " Modality ", imod, " value is wrong = ", trim(word)
               ierror = 1
               return
         endif

         if( icase < num_cases) case_data(imod, icase) =  case_value ! Load the case value only if we are not overshooting

    enddo   mod_loop

enddo cases_loop


if (icase /= num_cases ) then
  ierror  = 2
  write(err_msg,*) "Case_data has wrong size = ", num_cases, " ", icase
  num_cases = icase
endif


!--------------------------------------------------------------------------------------------
end subroutine get_cases
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------



!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine get_positivity(num_modalities, positiveislarge, ierror, err_msg)
! PURPOSE: Read the following line looking for an flags that signal positivity 
!          The flag starts with any of the flags for larger or small otherwise there 
!          is an input mistake. Flags are in  the header of the module. 
! NOTE:    Read only num_modalities tokens, assumes the rest is comments (to allow for comments
!          (there is no reason why there should be more info there if the format is in fact labmrmc

integer, intent(in):: num_modalities 
integer, dimension(num_modalities),intent(out):: positiveislarge ! where is positivity for the different modalities
integer, intent(out):: ierror        ! 0 ->  read successfully
                                     ! 1 ->  positivity information missing
                                     ! 2 ->  Line could not be read (corrupted file and so on)
character(len=line_length), intent(out) :: err_msg ! error message from subroutines
! internal variables
integer:: istat
character(len=line_length)  :: input_line ! character buffer that contains an input line
character(len=line_length)  :: word ! character buffer that contains an input line
integer:: i

ierror = 0

! Read the input like as a character line
 call get_input_line(input_line,  istat) 


! check whether the reading worked correctly
if (istat == - 1) then
    ierror = 1 ! the record is empty
    err_msg = "The positivity line is absent"
    return
elseif (istat >  0) then
    err_msg = "The positivity name can't be read " 
    ierror = 2
    return
end if

if(len_trim(input_line) == 0) then
  err_msg = " The positivity line is empty "
  ierror = 1
  return
else 
  ierror = 0 ! initialize the value to default
endif

! eliminate all the spaces to the left of the flags
input_line = trim(adjustl(input_line))

! Extract the first num_modalities words, and see if they are correct and contain info about
! positivity
do i = 1, num_modalities
    ! Get the first word from input file, and eliminate it form the string input_line
    call get_word(input_line, word, ierror)
    if (ierror == 1) then
       write(err_msg,*)  "Missing positivity flag for modality ", i
       ierror = 1
       return
    elseif ( verify(  word(1:1), small) == 0  ) then
       positiveislarge(i) = 0
    elseif ( verify( word(1:1), large) == 0) then
       positiveislarge(i) = 1
    else
        write(err_msg,*)  "wrong key code for positivity for modality ", i
        ierror = 1 ! wrong input
        return
    end if
enddo

end subroutine get_positivity

!---------------------------------------
subroutine get_modality_names(modality_names, num_modalities, ierror, err_msg)
!---------------------------------------
!PURPOSE: this subroutine reads in a free-text description of the 
!         different modality being analyzed. The modality names are expected to
!         be between the separator (see below) to allow a very flexible choice of names
!NOTE:    all terms that can be considered character are acceptable 
!         elements for the modality names (different languages, or habits might need them)
!         All the modality names must stay on a single line (for now)
!NOTE:    obviously names cannot contain the separator, but they can contain any other character
implicit none

integer, intent(inout):: num_modalities !
character (len = line_length), dimension(num_modalities), intent(out) :: modality_names
integer, intent(out):: ierror ! error flag, telling whether the file was incomplete or some error happened
        ! 0 ->  read successfully
        ! 1 ->  format is not labmrmc
        ! Title empty is is separated because it might be an uknown format or ???
character(len=line_length), intent(out) :: err_msg ! error message from subroutines
character (len = 1), parameter :: separator = '"' ! character that bounds modality names

! internal variables
integer:: istat
character(len=line_length)  :: input_line ! character buffer that contains an input line
character(len=line_length)  :: new_line ! character buffer that contains an input line
integer:: char_index

modality_names = ""


call get_input_line(input_line,  istat) 

! eliminate all the spaces to the left of the names
input_line = trim( adjustl(input_line) )

!     check whether the reading of the title worked correctly
if (istat == -1) then
    ierror = empty ! the file is empty
    err_msg = "The modality name line is absent"
    return
elseif (istat >  0) then
    err_msg = "The modality name can't be read " 
    ierror = unreadable
    return
end if

if(len_trim(input_line) == 0) then
  err_msg = " The modality name line is empty "
  ierror = 1
else 
  ierror = 0 ! initialize the value to default
endif

! extract the modality names, and their number

num_modalities = 0
do ! iterate until all the modalities are extracted
   ! work on the first separator (beginning of modality name)
   new_line = iblank ! set to "" the string
   char_index =  index(input_line, separator) ! find the first separator
   if(char_index == 0) exit ! there are no more modalities
   num_modalities = num_modalities + 1
   new_line = input_line(char_index+1:) ! eliminate all what is before the 1st separator
   input_line = new_line ! reload the line without the separator
   ! work on the second separator (end of modality name)
   new_line = iblank ! set to "" the string
   char_index =  index(input_line, separator) ! Find the second separator for this name
   if(char_index == 0 .or. len_trim(input_line(1:char_index-1) ) == 0 ) then ! second separator missing or name empty 
      ierror = 1
      err_msg = "Modality names format incompatible with labmrmc format"
      exit
   endif
   modality_names(num_modalities) =  input_line(1:char_index-1) ! extract the modality name
   new_line = input_line(char_index+1:) ! eliminate that modality nam
   input_line = new_line
enddo

if(num_modalities == 0) then
        ierror  = 1
        err_msg = "Modality names format incompatible with labmrmc format"
endif

!------------------------------------------------------------------------
end subroutine get_modality_names
!--------------------------------------------------------------------------
!-------------------------------------------------------------------------


!---------------------------------------
subroutine get_reader(reader_name,ierror, err_msg)
!---------------------------------------
!PURPOSE: this subroutine reads in a free-text description of the 
!         the reader being read. It also determines if the reader is the last
!         reader, that is if it starts with "#"

!NOTE:    all terms that can be considered character are acceptable 
!         elements for the reader name (different languages, or habits might need them)

implicit none

character (len = line_length),intent(out) :: reader_name
integer, intent(out):: ierror ! error flag, telling whether the file was incomplete or some error happened
        ! 0 -> reader was read successfully
        ! 1 -> reader is empty; negative is reading error, and follows the module general error flags (see above)
        ! 2 -> reader is last reader (first character is "#"
character(len=line_length), intent(out) :: err_msg ! error message from subroutines

character(len=line_length)  :: input_line ! character buffer that contains an input line

! internal variables
integer:: istat

!read input line, expected to be the reader

 call get_input_line(input_line,  istat) 

!     check whether the reading of the reader worked correctly
if (istat == -1) then
    ierror = 1 ! the file is empty
    err_msg = "The reader name line is absent"
    return
elseif (istat >  0) then
    err_msg = "The reader name can't be read " 
    ierror = 2
    return
end if

call  get_word( input_line, reader_name, ierror) 

if(ierror == 1) then
  err_msg = " The reader name is empty "
  ierror = 1
elseif ( verify(  reader_name(1:1)  , '#') == 0  ) then
  err_msg = " End of input file "
  ierror = 2
else
  ierror = 0 ! initialize the value to default
endif


 
!------------------------------------------------------------------------
end subroutine get_reader
!--------------------------------------------------------------------------
!-------------------------------------------------------------------------


!---------------------------------------
subroutine get_title(dataset_title,ierror, err_msg)
!---------------------------------------
!PURPOSE: this subroutine reads in a free-text description of the 
!         data, returns an error if the description line is not
!         there. 
!NOTE:    all terms that can be considered character are acceptable 
!         elements for the title (different languages, or habits might need them)

implicit none

character (len = line_length),intent(out) :: dataset_title
integer, intent(out):: ierror ! error flag, telling whether the file was incomplete or some error happened
        ! 0 -> title was read successfully
        ! 1 -> title is empty; negative is reading error, and follows the module general error flags (see above)
        ! Title empty is is separated because it might be an uknown format or ???
character(len=line_length), intent(out) :: err_msg ! error message from subroutines


! internal variables
integer:: istat


!read *the first* input line, that should be the title
 call get_input_line(dataset_title,  istat) 


!     check whether the reading of the title worked correctly
if (istat == -1) then
    ierror = empty ! the file is empty
    err_msg = "The input file seems to be empty: check file name "
    return
elseif (istat >  0) then
    err_msg = "The title of the input file can't be read " 
    ierror = unreadable
    return
end if

if(len_trim(dataset_title) == 0) then
  err_msg = " The title line is empty: check file format "
  ierror = 1
else 
  ierror = 0 ! initialize the value to default
endif

 
!------------------------------------------------------------------------
end subroutine get_title
!--------------------------------------------------------------------------
!-------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine get_word(input_line, word, ierror)
! PURPOSE: extract a word from an array of characters.
! NOTE:    words are defined as sequences of characters separated by other either by
!          the array boundaries, or by spaces of tabs
! NOTE2:   The word is removed from the array of characters, that is accordingly modified

 character(len=line_length), intent(inout)  :: input_line ! character buffer that contains an input line
 character(len = 2) :: blanks ! containing blank values
 character(len=line_length), intent(out)  :: word ! character buffer that contains an input line
 integer, intent(out):: ierror ! 0 -> word identified
                               ! 1 -> input character string is empty 


 integer:: istart, iend

 word = ""
 ierror = 0

 blanks = char(32)//char(9) ! Blanks are either blank or tab, one can here whatever suitable

 if(len_trim(input_line) == 0) then
    ierror = 1
    return
 endif

! Find the first non-blank character
  istart =  verify(  input_line, blanks)
! Find the end of the word, the result of scan is the first blank
  iend = istart - 2 + scan(  input_line(istart:), blanks)
! Load the word
  word = input_line(istart:iend)
! eliminate the word from input_line
 input_line = input_line(iend+1:)
 input_line(line_length -iend + 1:) = " "

end subroutine get_word
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
! BEGINS LABROC4 INPUT PROCEDURES
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine parse_labroc_list_data(num_normal_cases, num_abnormal_cases, ierror, err_msg)
! PURPOSE: check that the file has the right format and find what the number of cases is
! NOTE: Dec 2008 LP UC
! description of the list input as read by read_labroc_list_data
! 
! *******************************************************************************
!
!   required input and description of output:    
!
! ------------------------------------------------------------------------------
!
! For each data set, input the following on successive lines:  
!    1. up to line_length characters of free text to describe the data set.  
!    2. an alphabetic code word indicating whether small or large test-result values are 
!       associated with actually positive cases:   
!            "small" --> small values are associated with actually positive cases.
!            "large" --> large values are associated with actually positive cases.
!       (note: program reads only the first character of the input word. Caps or small is the same)
!    3 - . Two sequences of continuously-distributed test-result data. 
!          the data for actually negative cases should be entered first, followed by the data for actually positive
!          cases. the values in each sequence can be input in any order using 'free-format'
!          (integer, d format, f format, e format, comma separated format);
!          the only requirements are that the values must be separated by at least one blank space, and each
!          sequence must be terminated by an asterisk "*" (we do reccommend to use a new-line to separate
!          values, but it is not necessary) . each of the two sequences can contain up to any number of  entries. 
!          (note: this program can accept up to 6 digits to the left of the 
!          decimal point and up to 8 digits to the right of the decimal point.)
!
! INPUT EXAMPLE
! NOTE: the lines containing "BEGINNING" and "END" are not part of the
!       input file
!
! vvvvvvvvvvBEGINNING vvvvvvvvvvvvvvvvvvvvvvvvv
! sample input
! large  
! 2.5  
! -3.01   
! +4.2  
! -0.21
! 5.2
! 9.3
! -4.3
! -9.9
! +3.3
! 4.3
! *
! 0.221    
! 4.232    
! 0.453
! +3.21 
! -2.12
! -5.43 
! -3.2
! 8.3  
! 7.2  
! 9.1  
! 3.444 
! 8.334
! *
! ^END^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

implicit none

!     dummy arguments

integer, intent(out):: num_normal_cases
integer, intent(out):: num_abnormal_cases
integer, intent(out):: ierror ! error flag, ==1 wrong input, ==2 incomplete => needs improvement!
character(len=line_length), intent(out) :: err_msg ! error message from subroutines

character(len = line_length):: dataset_title
integer:: positiveislarge 

call read_header(positiveislarge, dataset_title ,ierror, err_msg)
!

if (ierror > 0) then
  ! error flag value for incomplete input =2, wrong input = 1
  return
endif

! For now it collapses all the error messages of read in in 2
! groups: incomplete input => 2, faulty input => 1
call parse_cases(num_normal_cases,ierror,err_msg)

if (ierror > 0) then
      err_msg = "Negative cases, " // err_msg(1:59)
      if(ierror /= 2) ierror = 1
      return
endif


call parse_cases(num_abnormal_cases,ierror,err_msg)
if (ierror > 0) then
      err_msg = "Positive cases, " // err_msg(1:59)
      if(ierror /= 2) ierror = 1
      return
endif

call rewind_input_file()

end subroutine parse_labroc_list_data
!------------------------------------------------
!------------------------------------------------

!----------------------------------
subroutine parse_cases(number_of_cases, ierror, err_msg)
!----------------------------------
! 
! PURPOSE  :this subroutine reads in a sequence of input data in free 
!         format. the only format requirements are that  (1)
!         any two input values must be separated by at least one  
!         iblank column and (2) the input data sequence must be terminated 
!         by an asterisk (*).
!         negative values are acceptable. it returns the number of cases found
! ALGORITHM: Read one input line, check whether the reading worked out, extract all values in 
!            that line, checking if they worked out, go to the next line.
! NOTE:  This code is complicated because it allows multiple inputs on the same line, which I think should
!        note be allowed, I kept it for now for legacy. LP - UC 2005
! NOTE: The algorithm in fact skips blank lines, as if they did not exist. This might be undesirable
! NOTE: Dec 2008 LP UC
implicit none

                                                                  ! each case
integer,intent(out):: number_of_cases ! number of cases of current truth read from the list input
integer, intent(out):: ierror ! error flag 
                       ! 0 => fine;  1 => can't read data;  2 => file incomplete;  3 => wrong format
character(len=line_length), intent(out) :: err_msg ! error message from subroutine

character(len=line_length) :: input_line ! character buffer that contains an input line
character(len=line_length) :: word ! character buffer that contains a word
character(len=1), dimension(line_length)  :: buffer ! character buffer that contains an input line
real(kind=double):: rvalue ! test result value currently being extracted
integer:: number_of_lines ! Number of lines read so far for this truth state
integer:: i ! iterator

integer:: istat ! failed reading flag
! 

number_of_cases = 0
number_of_lines = 0 ! this is used for error messaging

read_cases_of_this_truth: do
     number_of_lines = number_of_lines + 1
     call get_input_line(input_line,  istat)
     if (istat == -1) then ! Check whether the reading proceeded without errors
        ierror = 2 ! the input file is incomplete
        if (number_of_lines == 1) then
          write(err_msg, *) "Input file contains zero cases  "
        else
          write(err_msg, *) "Input file is truncated at line ", number_of_lines
        endif
        return
     elseif (istat >  0) then
        write(err_msg,*) "input line ", number_of_lines, "can't be read "
        ierror = 1  
        return
     end if
     read_this_line: do ! pick all the values in the line
            call get_word(input_line, word, ierror)
            if (ierror > 0) exit read_this_line ! We are done with this line, let's check the next
            if ( verify(  word(1:1)  , '*') == 0  ) then
                   ierror = 0
                   exit read_cases_of_this_truth ! stop reading and return the values of this truth state 
            endif
            ! Try to turn the extracted word into a number 
            do i = 1, line_length 
                buffer(i) = word(i:i)
            enddo

            call tonum(len_trim(word), buffer, rvalue, ierror)
            if (ierror > 0) then
                   ierror = 3 ! flag for wrong data format
                   ! This is a bit of a mess, used to avoid varialble formats. The truth is that
                   ! the lstring vectors should be eliminated, since f90 has better ways of dealing
                   ! with strings.
                   write(err_msg, *) number_of_lines,(lstring(i), i = 1,length) 
                   err_msg = adjustl(err_msg)
                   err_msg = 'wrong format at line ' // err_msg(1: line_length - 21)
                   return
            endif
            ! There was no error, so we load the new value
            number_of_cases = number_of_cases + 1
      enddo read_this_line
 enddo read_cases_of_this_truth

end subroutine parse_cases
!------------------------------------------------------------------
!------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine read_labroc_list_data(num_normal_cases, num_abnormal_cases, positiveislarge, dataset_title,&
                                input_data, ierror, err_msg)

! description of the list input as read by read_labroc_list_data
! 
! *******************************************************************************
!
!   required input and description of output:    
!
! ------------------------------------------------------------------------------
!
! For each data set, input the following on successive lines:  
!    1. up to line_length characters of free text to describe the data set.  
!    2. an alphabetic code word indicating whether small or large test-result values are 
!       associated with actually positive cases:   
!            "small" --> small values are associated with actually positive cases.
!            "large" --> large values are associated with actually positive cases.
!       (note: program reads only the first character of the input word. Caps or small is the same)
!    3 - . Two sequences of continuously-distributed test-result data. 
!          the data for actually negative cases should be entered first, followed by the data for actually positive
!          cases. the values in each sequence can be input in any order using 'free-format'
!          (integer, d format, f format, e format, comma separated format);
!          the only requirements are that the values must be separated by at least one blank space, and each
!          sequence must be terminated by an asterisk "*" (we do reccommend to use a new-line to separate
!          values, but it is not necessary) . each of the two sequences can contain up to any number of entries. 
!          (note: this program can accept up to 6 digits to the left of the 
!          decimal point and up to 8 digits to the right of the decimal point.)
!
! INPUT EXAMPLE
! NOTE: the lines containing "BEGINNING" and "END" are not part of the
!       input file
!
! vvvvvvvvvvBEGINNING vvvvvvvvvvvvvvvvvvvvvvvvv
! sample input
! large  
! 2.5  
! -3.01   
! +4.2  
! -0.21
! 5.2
! 9.3
! -4.3
! -9.9
! +3.3
! 4.3
! *
! 0.221    
! 4.232    
! 0.453
! +3.21 
! -2.12
! -5.43 
! -3.2
! 8.3  
! 7.2  
! 9.1  
! 3.444 
! 8.334
! *
! ^END^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

implicit none

!     dummy arguments
integer, intent(out):: positiveislarge 
integer, intent(in):: num_normal_cases
integer, intent(in):: num_abnormal_cases
real(kind=double),dimension(num_normal_cases+ num_abnormal_cases),intent(out):: input_data
character(len = line_length), intent(out):: dataset_title
integer, intent(out):: ierror ! error flag, ==1 wrong input, ==2 incomplete => needs improvement!
character(len=line_length), intent(out) :: err_msg ! error message from subroutines
integer:: mn, ms


call read_header(positiveislarge, dataset_title ,ierror, err_msg)
!

if (ierror > 0) then
  ! error flag value for incomplete input =2, wrong input = 1
  return
endif

! For now it collapses all the error messages of read in in 2
! groups: incomplete input => 2, faulty input => 1
call readin(mn,input_data(1:num_normal_cases),ierror,err_msg)

if(mn /= num_normal_cases) then
   ierror  = 1
   write(err_msg,*) "number of - cases found is different from expected", mn, num_normal_cases
endif

if (ierror > 0) then
      err_msg = "Negative cases, " // err_msg(1:59)
      if(ierror /= 2) ierror = 1
      return
endif


call readin(ms,input_data(num_normal_cases+1:num_normal_cases+num_abnormal_cases),ierror,err_msg)

if(ms /= num_abnormal_cases) then
   ierror  = 1
   write(err_msg,*) "number of + cases found is different from expected", ms, num_abnormal_cases
endif

if (ierror > 0) then
      err_msg = "Positive cases, " // err_msg(1:59)
      if(ierror /= 2) ierror = 1
      return
endif


end subroutine read_labroc_list_data
!------------------------------------------------
!------------------------------------------------


!---------------------------------------
subroutine read_header(positiveislarge,dataset_title,ierror, err_msg)
!---------------------------------------
!
!     this subroutine reads in a free-text description of the 
!     data and an alphabetic code word which indicates whether 
!     large values or small values are associated with actual
!     abnormality
!

implicit none

character (len = line_length),intent(out) :: dataset_title
integer,intent(out)::positiveislarge ! whether positivity is associated with larger or smaller values of
                                     ! the test value result, the number associated with each case in the
                                     ! input
integer, intent(out):: ierror ! error flag, telling whether the file was incomplete or some error happened
        ! needs improvement
character(len=line_length), intent(out) :: err_msg ! error message from subroutines

! internal variables

character(len=line_length)  :: input_line ! character buffer that contains an input line
integer:: istat

ierror = 0 ! initialize the value to default

!read *the first* input line, that should be the title
line = iblank
call get_input_line(dataset_title,  istat) 

!
!     check whether the reading of the title worked correctly
if (istat == -1) then
    ierror = 2 ! the file is empty
    err_msg = "The input file seems to be empty "
    return
elseif (istat >  0) then
    err_msg = "The title of the run can't be read " 
    ierror = 1  
    return
end if

if(len_trim(dataset_title) == 0) dataset_title = " THIS DATASET DID NOT HAVE A TITLE "

! read *the second* input line, the flag for positivity

input_line = ""

call get_input_line(input_line,  istat) 

 
!read (input_file,line_format,iostat = istat) (line(i), i = 1,line_length)
!
!     check how is the reading going

if (istat == -1) then
    ierror = 2 ! the file is empty
    err_msg = "The input file seems to be contain only the title of the dataset "
    return
elseif (istat >  0) then
    err_msg = "Only the title can be read for this dataset " 
    ierror = 1  
    return
end if


! check if the input is empty
if ( len_trim(input_line) == 0) then
  err_msg = "the line telling where is positive was empty, cannot proceed"
  ierror = 2 ! flag for incomplete input
  return
end if

input_line = adjustl(input_line)

! look if the second line starts with any of the flags for larger or small (see above) otherwise there 
! is an input mistake. 
if     ( verify(  input_line(1:1), small) == 0  ) then
  positiveislarge = 0
elseif ( verify( input_line(1:1), large) == 0) then
  positiveislarge = 1
else
  err_msg = "wrong key code for location of higher probability of positive test results"
  ierror = 1 ! wrong input
  return
end if
 
!------------------------------------------------------------------------
end subroutine read_header
!--------------------------------------------------------------------------
!-------------------------------------------------------------------------

!----------------------------------
subroutine readin(number_of_cases, test_result_values, ierror, err_msg)
!----------------------------------
! 
! PURPOSE  :this subroutine reads in a sequence of input data in free 
!         format. the only format requirements are that  (1)
!         any two input values must be separated by at least one  
!         iblank column and (2) the input data sequence must be terminated 
!         by an asterisk (*).
!         negative values are acceptable.
! ALGORITHM: Read one input line, check whether the reading worked out, extract all values in 
!            that line, checking if they worked out, go to the next line.
! NOTE:  This code is complicated because it allows multiple inputs on the same line, which I think should
!        note be allowed, I kept it for now for legacy. LP - UC 2005
! NOTE: The algorithm in fact skips blank lines, as if they did not exist. This might be undesirable

implicit none

real(kind=double), intent(out), dimension(*):: test_result_values ! array with the value associated with
                                                                  ! each case
integer,intent(out):: number_of_cases ! number of cases of current truth read from the list input
integer, intent(out):: ierror ! error flag 
                       ! 0 => fine;  1 => can't read data;  2 => file incomplete;  3 => wrong format
character(len=line_length), intent(out) :: err_msg ! error message from subroutine


character(len=line_length) :: input_line ! character buffer that contains an input line
character(len=line_length) :: word ! character buffer that contains a word
character(len=1), dimension(line_length)  :: buffer ! character buffer that contains an input line
real(kind=double):: rvalue ! test result value currently being extracted
integer:: number_of_lines ! Number of lines read so far for this truth state
integer:: i ! iterator

integer:: istat ! failed reading flag
! 

number_of_cases = 0
number_of_lines = 0 ! this is used for error messaging

read_cases_of_this_truth: do
     number_of_lines = number_of_lines + 1
     call get_input_line(input_line,  istat)
     if (istat == -1) then ! Check whether the reading proceeded without errors
        ierror = 2 ! the input file is incomplete
        if (number_of_lines == 1) then
          write(err_msg, *) "Input file contains zero cases  "
        else
          write(err_msg, *) "Input file is truncated at line ", number_of_lines
        endif
        return
     elseif (istat >  0) then
        write(err_msg,*) "input line ", number_of_lines, "can't be read "
        ierror = 1  
        return
     end if
     read_this_line: do ! pick all the values in the line
            call get_word(input_line, word, ierror)
            if (ierror > 0) exit read_this_line ! We are done with this line, let's check the next
            if ( verify(  word(1:1)  , '*') == 0  ) then
                   ierror = 0
                   exit read_cases_of_this_truth ! stop reading and return the values of this truth state 
            endif
            ! Try to turn the extracted word into a number 
            do i = 1, line_length 
                buffer(i) = word(i:i)
            enddo

            call tonum(len_trim(word), buffer, rvalue, ierror)
            if (ierror > 0) then
                   ierror = 3 ! flag for wrong data format
                   ! This is a bit of a mess, used to avoid varialble formats. The truth is that
                   ! the lstring vectors should be eliminated, since f90 has better ways of dealing
                   ! with strings.
                   write(err_msg, *) number_of_lines,(lstring(i), i = 1,length) 
                   err_msg = adjustl(err_msg)
                   err_msg = 'wrong format at line ' // err_msg(1: line_length - 21)
                   return
            endif
            ! There was no error, so we load the new value
            number_of_cases = number_of_cases + 1
            test_result_values(number_of_cases) = rvalue
      enddo read_this_line
 enddo read_cases_of_this_truth

end subroutine readin
!------------------------------------------------------------------
!------------------------------------------------------------------



!-----------------------------------------------------------
subroutine tonum(len,linput,rnum,ierr)
!-----------------------------------------------------------
!
! PURPOSE: Attempt to convert character string to real number
!          routine completely overhauled to allow it to read
!          real number formats with d and e (and , as well
!          european convention for decimals (",") 
!          returns error if conversion can't be done
! NOTE:    It will in fact work with mixed formats, which isn't necessarily a good practice...
! NOTE:    cate must be paid when dealing with comma delimited formats, since the european dot
!          is a comma.
implicit none

character, dimension(line_length), intent(in):: linput ! character string where the currently read line is
integer, intent(in):: len ! lenght of currently extracted word
integer, intent(out):: ierr ! error flag
real(kind=double), intent(out):: rnum  ! the actual number, once the word is converted into reals

character (len=line_length) buff,temp
integer:: ie, i_dot, in

rnum = 0.0_double
buff(:)=''

! Here we are assuming that since linput was read, it can also be written
write(buff,*) (linput(in),in=1,len)

 buff = adjustl(buff)
 ! Check if it is e or a d format
 ie = index(buff,'e')
 if (ie == 0) ie = index(buff,'E') 
 if (ie == 0) ie = index(buff,'d') ! in case the format out d format
 if (ie == 0) ie = index(buff,'D')

 ! check whether there is a dot, a comma, or neither
 i_dot = index(buff,'.')
 if (i_dot == 0) then
    i_dot=index(buff,',') ! check for the european format
    ! if it is a comma, convert it to a dot
    if(i_dot > 0)  buff = buff(1:i_dot-1)//'.'//buff(i_dot+1:)
 endif
 
! If there is no dot and no e, then it must be either an integer
! or it is an unknown format
 if (i_dot == 0 .and. ie == 0) then
        buff(len+1:)    ='.'
 else if (i_dot == 0) then
       temp             = buff(ie:len+1)
       buff(ie:ie)      = '.'
       buff(ie+1:len+2) = temp
 endif

 read(buff,*,iostat=ierr) rnum
 ! Check whether the resulting buffer is in fact a number and can be read as such
 if (ierr /= 0) ierr = 1
    
end subroutine tonum
!-----------------------------------------------------
!-----------------------------------------------------


!-----------------------------------------------------



  
end module read_ROC_input
