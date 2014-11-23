MODULE IO

  USE sweeper

  IMPLICIT NONE

  PUBLIC :: processCmdLine
  PUBLIC :: closeFiles

  INTEGER,PARAMETER :: inpFileUnitNo=123

  CONTAINS
!===============================================================================
    SUBROUTINE processCmdLine()
      ! Local Variables
      CHARACTER(LEN=16) :: arg_in

      CALL GET_COMMAND_ARGUMENT(1,arg_in)
      OPEN(FILE=TRIM(ADJUSTL(arg_in)),UNIT=inpFileUnitNo)

    END SUBROUTINE processCmdLine
!===============================================================================
    SUBROUTINE closeFiles()

      CLOSE(UNIT=inpFileUnitNo)

    END SUBROUTINE closeFiles
END MODULE IO

