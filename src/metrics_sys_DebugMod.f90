! Last modified December 5 by Lorenzo Pesce and Kevin Schartz to replace
! the non-standard "\" format with "ADVANCE='NO'" and calls to
! AINT with calls to FLOOR.


! Last modified on August 7, 2006 to add the individual debugging flag variables for
! separate functions such as input, jackknifing, ANOVA, the various curve-fitters, etc.

! Last modified on April 14, 2006, to remove all array declarations of the form
! 0:N that began array indices at 0.  Calling programs have been modified so that
! the standard Fortran form of arrays beginning with an index of 1 is observed.
! Thus, it is not necessary for these routines to have arrays beginning at 0.


! This module contains variables used for debugging.
! Schartz, 4/3/2003
MODULE Debugging

    USE AIF_Format  ! Schartz, 6/17/2003

    ! Debugging variables:
    INTEGER, PARAMETER :: debugUnit = 3954    ! debug unit number
    LOGICAL :: Debug          ! Global debugging flag
    LOGICAL :: UIDebug        ! debug user interface components (e.g., dialog boxes)

    LOGICAL, PUBLIC :: DebugInput = .FALSE._1
    LOGICAL, PUBLIC :: DebugJackknife = .FALSE._1
    LOGICAL, PUBLIC :: DebugRSCORE = .FALSE._1
    LOGICAL, PUBLIC :: DebugSens = .FALSE._1
    LOGICAL, PUBLIC :: DebugPSV = .FALSE._1
    LOGICAL, PUBLIC :: DebugCollapse = .FALSE._1
    LOGICAL, PUBLIC :: DebugPROPROC = .FALSE._1
    LOGICAL, PUBLIC :: DebugBIGAMMA = .FALSE._1
    LOGICAL, PUBLIC :: DebugCBM = .FALSE._1
    LOGICAL, PUBLIC :: DebugANOVA = .FALSE._1

    ! Define an INTERFACE to allow for overloading of definitions so that a single calling
    ! statement can be used to refer to the INTEGER, REAL*8, and LOGICAL variants.  By
    ! putting the INTERFACE statement in the module, the generic name DisplayState is available
    ! to any procedure which accesses the module via a USES statement.
    INTERFACE DisplayState

        MODULE PROCEDURE DisplayStateMsg
        MODULE PROCEDURE DisplayStateInt
        MODULE PROCEDURE DisplayStateRL8
        MODULE PROCEDURE DisplayStateLog
        MODULE PROCEDURE DisplayStateText
        MODULE PROCEDURE DisplayStateIntArray1
        MODULE PROCEDURE DisplayStateIntArray01
        MODULE PROCEDURE DisplayStateRL8Array1
        MODULE PROCEDURE DisplayStateRL8Array01
        MODULE PROCEDURE DisplayStateRL8Array2
        MODULE PROCEDURE DisplayStateRL8Array2R
        MODULE PROCEDURE DisplayStateIntArray2
        MODULE PROCEDURE DisplayStateIntArray2R
        MODULE PROCEDURE DisplayStateRL8Array02

    END INTERFACE

CONTAINS
!********************************************************************************************
!********************************************************************************************
!																						API
!********************************************************************************************
!********************************************************************************************

!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateMsg ( tempText )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Modified March 9, 2006 to accommodate zero-length strings.
    IF (length > 0 ) THEN
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)
    ELSE
        continue! WRITE(debugUnit,*)
    ENDIF


END SUBROUTINE DisplayStateMsg



!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateInt ( tempText, intVal )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    INTEGER :: intVal
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
    continue! WRITE(debugUnit, '('//AFMT(length)//'," ",I10)') TRIM(text), intVal

END SUBROUTINE DisplayStateInt

!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateRL8 ( tempText, rl8Val )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    REAL(KIND=selected_real_kind(p=15)), INTENT(IN) :: rl8Val
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
    continue! WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), rl8Val


END SUBROUTINE DisplayStateRL8


!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateLog ( tempText, logVal )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    LOGICAL :: logVal
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length !, length1, length2

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
    !WRITE(debugUnit, '('//AFMT(length)//'," ",L1)') text, logVal

    IF (logVal) THEN
      continue! WRITE(debugUnit, '('//AFMT(length)//',A6)') TRIM(text), ": TRUE"
    ELSE
      continue! WRITE(debugUnit, '('//AFMT(length)//',A7)') TRIM(text), ": FALSE"
    ENDIF


END SUBROUTINE DisplayStateLog


!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateText ( tempText, tempTextVal )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    CHARACTER*80 :: text
    CHARACTER*255 :: textVal
    CHARACTER(LEN=*) :: tempText,tempTextVal
    INTEGER :: length1, length2

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)
    textVal = TRIM(tempTextVal)

    ! Get the length of the text argument:
    length1 = LEN_TRIM(text)
    length2 = LEN_TRIM(textVal)

    ! Write the text & value to the debugging window:
    continue! WRITE(debugUnit, '('//AFMT(length1)//'," ",'//AFMT(length2)//')') &
    continue! &  TRIM(text), TRIM(textVal)


END SUBROUTINE DisplayStateText

!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateIntArray1 ( tempText, IntVal, Index1 )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    INTEGER :: Index1
    INTEGER :: IntVal(Index1)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length, I
    INTEGER, PARAMETER :: MaxColumns = 10

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), IntVal
    IF (Index1 <= MaxColumns) THEN

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        DO I = 1, Index1
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") I
        END DO
        continue! WRITE(debugUnit,*)

        DO I = 1, Index1
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") IntVal(I)
        END DO
        continue! WRITE(debugUnit,*)

    ELSE ! write as rows rather than columns

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        DO I = 1, Index1
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") I
            continue! WRITE(debugUnit,'(I10, " ")') IntVal(I)
        END DO
        continue! WRITE(debugUnit,*)

    ENDIF

END SUBROUTINE DisplayStateIntArray1


!*******************************************************************************
!===============================================================================
SUBROUTINE DisplayStateIntArray01 ( tempText, IntVal, LowerIndex, UpperIndex )
!===============================================================================
!*******************************************************************************

    IMPLICIT NONE

    INTEGER :: LowerIndex,UpperIndex
    INTEGER :: IntVal(UpperIndex)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length, I
    INTEGER, PARAMETER :: MaxColumns = 10

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), IntVal
    IF (UpperIndex <= MaxColumns) THEN

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        DO I = LowerIndex, UpperIndex
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") I
        END DO
        continue! WRITE(debugUnit,*)

        DO I = LowerIndex, UpperIndex
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") IntVal(I)
        END DO
        continue! WRITE(debugUnit,*)

    ELSE ! write as rows rather than columns

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        DO I = LowerIndex, UpperIndex
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") I
            continue! WRITE(debugUnit,'(I10, " ")') IntVal(I)
        END DO
        continue! WRITE(debugUnit,*)

    ENDIF

END SUBROUTINE DisplayStateIntArray01


!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateRL8Array1 ( tempText, rl8Val, UpperIndex )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    INTEGER :: UpperIndex
    REAL(KIND=selected_real_kind(p=15)), INTENT(IN) :: rl8Val(UpperIndex)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length, I
    INTEGER, PARAMETER :: MaxColumns = 10

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), rl8Val
    IF (UpperIndex <= MaxColumns) THEN

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        DO I = 1, UpperIndex
            continue! WRITE(debugUnit,'(I12, " ")', ADVANCE = "NO") I
        END DO
        continue! WRITE(debugUnit,*)

        DO I = 1, UpperIndex
            continue! WRITE(debugUnit,'(F12.5, " ")', ADVANCE = "NO") rl8Val(I)
        END DO
        continue! WRITE(debugUnit,*)

    ELSE ! write as rows rather than columns

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        DO I = 1, UpperIndex
            continue! WRITE(debugUnit,'(I12, " ")', ADVANCE = "NO") I
            continue! WRITE(debugUnit,'(F12.5, " ")') rl8Val(I)
        END DO
        continue! WRITE(debugUnit,*)

    ENDIF

END SUBROUTINE DisplayStateRL8Array1

!*******************************************************************************
!===============================================================================
SUBROUTINE DisplayStateRL8Array01 ( tempText, rl8Val, LowerIndex, UpperIndex )
!===============================================================================
!*******************************************************************************

    IMPLICIT NONE

    INTEGER :: LowerIndex,UpperIndex
    REAL(KIND=selected_real_kind(p=15)), INTENT(IN) :: rl8Val(UpperIndex)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length, I
    INTEGER, PARAMETER :: MaxColumns = 10

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), rl8Val
    IF (UpperIndex <= MaxColumns) THEN

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        DO I = LowerIndex, UpperIndex
            continue! WRITE(debugUnit,'(I12, " ")', ADVANCE = "NO") I
        END DO
        continue! WRITE(debugUnit,*)

        DO I = LowerIndex, UpperIndex
            continue! WRITE(debugUnit,'(F12.5, " ")', ADVANCE = "NO") rl8Val(I)
        END DO
        continue! WRITE(debugUnit,*)

    ELSE ! write as rows rather than columns

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        DO I = LowerIndex, UpperIndex
            continue! WRITE(debugUnit,'(I12, " ")', ADVANCE = "NO") I
            continue! WRITE(debugUnit,'(F12.5, " ")') rl8Val(I)
        END DO
        continue! WRITE(debugUnit,*)

    ENDIF

END SUBROUTINE DisplayStateRL8Array01


!***************************************************************************************
!=======================================================================================
SUBROUTINE DisplayStateRL8Array2 ( tempText, rl8Val, UpperRowIndex, UpperColumnIndex )
!=======================================================================================
!***************************************************************************************

    IMPLICIT NONE

    INTEGER :: UpperRowIndex, UpperColumnIndex
    REAL(KIND=selected_real_kind(p=15)), INTENT(IN) :: rl8Val(UpperRowIndex,UpperColumnIndex)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length, I, J

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), rl8Val
    continue! WRITE(debugUnit,*)
    continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

    continue! WRITE(debugUnit,'("      ")', ADVANCE = "NO")
    DO J = 1, UpperColumnIndex
        continue! WRITE(debugUnit,'(I12, " ")', ADVANCE = "NO") J
    END DO
    continue! WRITE(debugUnit,*)

    DO I = 1, UpperRowIndex
        continue! WRITE(debugUnit,'(I4,"  ")', ADVANCE = "NO") I
        DO J = 1, UpperColumnIndex
            continue! WRITE(debugUnit,'(F12.5, " ")', ADVANCE = "NO") rl8Val(I,J)
        END DO
        continue! WRITE(debugUnit,*)
    END DO
    continue! WRITE(debugUnit,*)

END SUBROUTINE DisplayStateRL8Array2


!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateRL8Array2R ( tempText, rl8Val, UpperRowIndex, UpperColumnIndex, Reverse )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    INTEGER :: UpperRowIndex, UpperColumnIndex
    REAL(KIND=selected_real_kind(p=15)), INTENT(IN) :: rl8Val(UpperRowIndex,UpperColumnIndex)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    LOGICAL :: Reverse
    INTEGER :: length, I, J

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), rl8Val

    IF (Reverse) THEN

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        continue! WRITE(debugUnit,'("      ")', ADVANCE = "NO")
        DO I = 1, UpperRowIndex
            continue! WRITE(debugUnit,'(I12, " ")', ADVANCE = "NO") I
        END DO
        continue! WRITE(debugUnit,*)

        DO J = 1, UpperColumnIndex
            continue! WRITE(debugUnit,'(I4,"  ")', ADVANCE = "NO") J
            DO I = 1, UpperRowIndex
            continue! WRITE(debugUnit,'(F12.5, " ")', ADVANCE = "NO") rl8Val(I,J)
            END DO
            continue! WRITE(debugUnit,*)
        END DO
        continue! WRITE(debugUnit,*)

    ELSE

        ! Call other subroutine to print without reversing arrays:
        CALL DisplayStateRL8Array2(tempText,rl8Val,UpperRowIndex,UpperColumnIndex)

    ENDIF

END SUBROUTINE DisplayStateRL8Array2R

!***********************************************************************
!=======================================================================
SUBROUTINE DumpUnitInfo ( fileUnit )
!=======================================================================
!***********************************************************************
    INTEGER :: fileUnit
    LOGICAL :: lopened,lnamed
    CHARACTER(64):: fname

    continue! WRITE(*,*) ' --------------------------------------------------'
   ! WRITE(*,*) ' Displaying unit connection status for unit: ', fileUnit
    INQUIRE(UNIT=fileUnit,OPENED=lopened,NAMED=lnamed,NAME=fname)
    continue! WRITE(debugUnit,*) ' UNIT = ',fileUnit
    continue! WRITE(debugUnit,*) ' OPENED = ',lopened
    continue! WRITE(debugUnit,*) ' NAMED = ',lnamed
    continue! WRITE(debugUnit,*) ' NAME = ',fname
    continue! WRITE(debugUnit,*) ' --------------------------------------------------'

END SUBROUTINE DumpUnitInfo


!***************************************************************************************
!=======================================================================================
SUBROUTINE DisplayStateIntArray2 ( tempText, intVal, UpperRowIndex, UpperColumnIndex )
!=======================================================================================
!***************************************************************************************

    IMPLICIT NONE

    INTEGER :: UpperRowIndex, UpperColumnIndex
    INTEGER :: intVal(UpperRowIndex,UpperColumnIndex)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length, I, J

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), intVal
    continue! WRITE(debugUnit,*)
    continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

    continue! WRITE(debugUnit,'("      ")', ADVANCE = "NO")
    DO J = 1, UpperColumnIndex
        continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") J
    END DO
    continue! WRITE(debugUnit,*)

    DO I = 1, UpperRowIndex
        continue! WRITE(debugUnit,'(I4,"  ")', ADVANCE = "NO") I
        DO J = 1, UpperColumnIndex
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") intVal(I,J)
        END DO
        continue! WRITE(debugUnit,*)
    END DO
    continue! WRITE(debugUnit,*)

END SUBROUTINE DisplayStateIntArray2

!***********************************************************************
!=======================================================================
SUBROUTINE DisplayStateIntArray2R ( tempText, intVal, UpperRowIndex, UpperColumnIndex, Reverse )
!=======================================================================
!***********************************************************************

    IMPLICIT NONE

    INTEGER :: UpperRowIndex, UpperColumnIndex
    INTEGER :: intVal(UpperRowIndex,UpperColumnIndex)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    LOGICAL :: Reverse
    INTEGER :: length, I, J

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), intVal

    IF (Reverse) THEN

        continue! WRITE(debugUnit,*)
        continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

        continue! WRITE(debugUnit,'("      ")', ADVANCE = "NO")
        DO I = 1, UpperRowIndex
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") I
        END DO
        continue! WRITE(debugUnit,*)

        DO J = 1, UpperColumnIndex
            continue! WRITE(debugUnit,'(I4,"  ")', ADVANCE = "NO") J
            DO I = 1, UpperRowIndex
            continue! WRITE(debugUnit,'(I10, " ")', ADVANCE = "NO") intVal(I,J)
            END DO
            continue! WRITE(debugUnit,*)
        END DO
        continue! WRITE(debugUnit,*)

    ELSE

        ! Call other subroutine to print without reversing arrays:
        CALL DisplayStateIntArray2(tempText,intVal,UpperRowIndex,UpperColumnIndex)

    ENDIF

END SUBROUTINE DisplayStateIntArray2R

!***************************************************************************************
!=======================================================================================
SUBROUTINE DisplayStateRL8Array02 ( tempText, rl8Val, UpperRowIndex, LowerColumnIndex, UpperColumnIndex )
!=======================================================================================
!***************************************************************************************

    IMPLICIT NONE

    INTEGER :: UpperRowIndex, LowerColumnIndex,UpperColumnIndex
    REAL(KIND=selected_real_kind(p=15)), INTENT(IN) :: rl8Val(UpperRowIndex,UpperColumnIndex)
    CHARACTER*80 :: text
    CHARACTER(LEN=*) :: tempText
    INTEGER :: length, I, J

    ! Write the passed text string to a local variable:
    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ! Write the text & value to the debugging window:
!    WRITE(debugUnit, '('//AFMT(length)//'," ",F15.8)') TRIM(text), rl8Val
    continue! WRITE(debugUnit,*)
    continue! WRITE(debugUnit, '('//AFMT(length)//')') TRIM(text)

    continue! WRITE(debugUnit,'("      ")', ADVANCE = "NO")
    DO J = LowerColumnIndex, UpperColumnIndex
        continue! WRITE(debugUnit,'(I12, " ")', ADVANCE = "NO") J
    END DO
    continue! WRITE(debugUnit,*)

    DO I = 1, UpperRowIndex
        continue! WRITE(debugUnit,'(I4,"  ")', ADVANCE = "NO") I
        DO J = LowerColumnIndex, UpperColumnIndex
            continue! WRITE(debugUnit,'(F12.5, " ")', ADVANCE = "NO") rl8Val(I,J)
        END DO
        continue! WRITE(debugUnit,*)
    END DO
    continue! WRITE(debugUnit,*)

END SUBROUTINE DisplayStateRL8Array02


!********************************************************************************************
!********************************************************************************************
END MODULE Debugging


