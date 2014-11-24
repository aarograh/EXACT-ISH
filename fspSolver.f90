MODULE fspSolver

  USE sweeper

  PUBLIC :: fspSolverType

  TYPE :: fspSolverType
    INTEGER :: igstt = 0
    INTEGER :: igstp = 0
    CLASS(SourceType),POINTER :: thisSource
    CLASS(sweeperType),POINTER :: thisSweeper => NULL()
    CONTAINS
      PROCEDURE,PASS :: initialize
  END TYPE fspSolverType

  CONTAINS
!===============================================================================
    SUBROUTINE initialize(mySolver,thisSweeper)
      CLASS(fspSolverType),INTENT(INOUT) :: mySolver
      CLASS(sweeperType),TARGET,INTENT(INOUT) :: thisSweeper



    END SUBROUTINE initialize
END MODULE fspSolver
