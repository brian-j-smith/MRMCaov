! This module contains predefined arrays and functions that are used to facilitate the
! construction of dynamic FORMAT statements based on the length of the quantity to
! be written.
!
! Created by: Kevin Schartz
! Created on: June 17,2003
!
! Based on: A_Format_255.F90
!
! Last modified by: Kevin Schartz
! Last modified on: March 16, 2005

! Last modified December 1, 2006 by Lorenzo Pesce to replace
! calls to AINT with calls to FLOOR.  Also replaced
! "REAL*8" and "REAL(KIND=8" with "REAL(KIND=selected_real_kind(p=15))"

MODULE AIF_Format

  CHARACTER(4) AFMT(0:255)
  CHARACTER(2) IFMT(0:9)
  CHARACTER(5) FFMT0(0:9),FFMT1(0:9),FFMT2(0:9),FFMT3(0:9),FFMT4(0:9),FFMT5(0:9)

  DATA AFMT /&
             '  A0',& ! added for zero-length strings
             '  A1','  A2','  A3','  A4','  A5','  A6','  A7','  A8','  A9',' A10',&
             ' A11',' A12',' A13',' A14',' A15',' A16',' A17',' A18',' A19',' A20',&
             ' A21',' A22',' A23',' A24',' A25',' A26',' A27',' A28',' A29',' A30',&
             ' A31',' A32',' A33',' A34',' A35',' A36',' A37',' A38',' A39',' A40',&
             ' A41',' A42',' A43',' A44',' A45',' A46',' A47',' A48',' A49',' A50',&
             ' A51',' A52',' A53',' A54',' A55',' A56',' A57',' A58',' A59',' A60',&
             ' A61',' A62',' A63',' A64',' A65',' A66',' A67',' A68',' A69',' A70',&
             ' A71',' A72',' A73',' A74',' A75',' A76',' A77',' A78',' A79',' A80',&
             ' A81',' A82',' A83',' A84',' A85',' A86',' A87',' A88',' A89',' A90',&
             ' A91',' A92',' A93',' A94',' A95',' A96',' A97',' A98',' A99','A100',&
             'A101','A102','A103','A104','A105','A106','A107','A108','A109','A110',&
             'A111','A112','A113','A114','A115','A116','A117','A118','A119','A120',&
             'A121','A122','A123','A124','A125','A126','A127','A128','A129','A130',&
             'A131','A132','A133','A134','A135','A136','A137','A138','A139','A140',&
             'A141','A142','A143','A144','A145','A146','A147','A148','A149','A150',&
             'A151','A152','A153','A154','A155','A156','A157','A158','A159','A160',&
             'A161','A162','A163','A164','A165','A166','A167','A168','A169','A170',&
             'A171','A172','A173','A174','A175','A176','A177','A178','A179','A180',&
             'A181','A182','A183','A184','A185','A186','A187','A188','A189','A190',&
             'A191','A192','A193','A194','A195','A196','A197','A198','A199','A200',&
             'A201','A202','A203','A204','A205','A206','A207','A208','A209','A210',&
             'A211','A212','A213','A214','A215','A216','A217','A218','A219','A220',&
             'A221','A222','A223','A224','A225','A226','A227','A228','A229','A230',&
             'A231','A232','A233','A234','A235','A236','A237','A238','A239','A240',&
             'A241','A242','A243','A244','A245','A246','A247','A248','A249','A250',&
             'A251','A252','A253','A254','A255'/

  DATA IFMT /'I0','I1','I2','I3','I4','I5','I6','I7','I8','I9'/

  DATA FFMT0 /' F0.0',' F1.0',' F2.0',' F3.0',' F4.0',' F5.0',' F6.0',' F7.0',' F8.0',' F9.0'/
  DATA FFMT1 /' F1.1',' F2.1',' F3.1',' F4.1',' F5.1',' F6.1',' F7.1',' F8.1',' F9.1','F10.1'/
  DATA FFMT2 /' F2.2',' F3.2',' F4.2',' F5.2',' F6.2',' F7.2',' F8.2',' F9.2','F10.2','F11.2'/
  DATA FFMT3 /' F3.3',' F4.3',' F5.3',' F6.3',' F7.3',' F8.3',' F9.3','F10.3','F11.3','F12.3'/
  DATA FFMT4 /' F4.4',' F5.4',' F6.4',' F7.4',' F8.4',' F9.4','F10.4','F11.4','F12.4','F13.4'/
  DATA FFMT5 /' F5.5',' F6.5',' F7.5',' F8.5',' F9.5','F10.5','F11.5','F12.5','F13.5','F14.5'/

CONTAINS
!********************************************************************************************
!********************************************************************************************
!                                           API
!********************************************************************************************
!********************************************************************************************

!***********************************************************************
!=======================================================================
  CHARACTER*2 FUNCTION ICode (intVal)
!=======================================================================
!***********************************************************************
    INTEGER, INTENT(IN) :: intVal
    INTEGER :: length

    ! modified by Schartz, March 15, 2005
!    length = (intVal/10)+1
    length = FLOOR(   LOG10( DBLE(intVal) )   ) + 1

    ICode = IFMT(length)

  END FUNCTION ICode

!***********************************************************************
!=======================================================================
  CHARACTER*5 FUNCTION FCode (rl8Val,numDecimal)
!=======================================================================
! Change by LP on 12/1/06 -> replaced aint with floor
!***********************************************************************
    REAL(KIND=selected_real_kind(p=15)), INTENT(IN) :: rl8Val
    INTEGER, INTENT(IN) :: numDecimal ! number of decimal places

    INTEGER :: length

    ! modified by Schartz, March 16, 2005
!    length = (AINT(rl8Val)/10)+2
    length =  FLOOR(   LOG10(rl8Val)   ) + 2

    ! Modified to handle rl8Val arguments < 1.0
    ! Schartz, March 16, 2005
    IF (length < 1) THEN
      length = 1
    ENDIF

    SELECT CASE(numDecimal)

      CASE(0)
        FCode = FFMT0(length)
      CASE(1)
        FCode = FFMT1(length)
      CASE(2)
        FCode = FFMT2(length)
      CASE(3)
        FCode = FFMT3(length)
      CASE(4)
        FCode = FFMT4(length)
      CASE(5)
        FCode = FFMT5(length)

    END SELECT

  END FUNCTION FCode

!***********************************************************************
!=======================================================================
  CHARACTER*4 FUNCTION ACode (tempText)
!=======================================================================
!***********************************************************************
  CHARACTER*(*) :: tempText
  CHARACTER*80 :: text
  INTEGER :: length

    text = TRIM(tempText)

    ! Get the length of the text argument:
    length = LEN_TRIM(text)

    ACode = AFMT(length)

  END FUNCTION ACode

!********************************************************************************************
!********************************************************************************************
END MODULE
