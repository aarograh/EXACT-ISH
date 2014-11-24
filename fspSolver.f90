MODULE fspSolver

  USE sweeper

  PUBLIC :: fspSolverType

  TYPE :: fspSolverType
    INTEGER :: igstt = 0
    INTEGER :: igstp = 0
    CLASS(SourceType),POINTER :: thisSource
    CLASS(sweeperType),POINTER :: thisSweeper => NULL()
    CONTAINS
      PROCEDURE,PASS :: initialize => initializeFspSolver
  END TYPE fspSolverType

  CONTAINS
!===============================================================================
    SUBROUTINE initializeFspSolver(mySolver,thisSweeper)
      CLASS(fspSolverType),INTENT(INOUT) :: mySolver
      CLASS(sweeperType),TARGET,INTENT(INOUT) :: thisSweeper

      mySolver%thisSweeper => thisSweeper
      mySolver%igstt = thisSweeper%igstt
      mySolver%igstp = thisSweeper%igstp

      CALL thisSweeper%initialize(mySolver%thisSource)

    END SUBROUTINE initializeFspSolver
END MODULE fspSolver
