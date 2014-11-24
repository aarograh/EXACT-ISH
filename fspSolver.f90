MODULE fspSolver

  USE sweeper

  IMPLICIT NONE

  PUBLIC :: fspSolverType

  TYPE :: fspSolverType
    INTEGER :: igstt = 0
    INTEGER :: igstp = 0
    DOUBLE PRECISION,POINTER :: psi(:) => NULL()
    CLASS(SourceType),POINTER :: thisSource
    CLASS(sweeperType),POINTER :: thisSweeper => NULL()
    CONTAINS
      PROCEDURE,PASS :: initialize => initializeFspSolver
      PROCEDURE,PASS :: step => stepFspSolver
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
!===============================================================================
    SUBROUTINE stepFspSolver(mySolver)
      CLASS(fspSolverType),INTENT(INOUT) :: mySolver
      ! Local Variables
      INTEGER :: ig

      DO ig=mySolver%igstt,mySolver%igstp
        ! Set up source
        CALL mySolver%thisSource%initExtSource(ig)
        CALL mySolver%thisSource%computeMGFS(ig,mySolver%psi)
        CALL mySolver%thisSource%updateInScatter(ig,mySolver%igstt,mySolver%igstp)
        CALL mySolver%thisSweeper%setExtSource(mySolver%thisSource)
        ! Perform sweep
        CALL mySolver%thisSweeper%sweep(ig,3,1.0D-03)
      ENDDO

    END SUBROUTINE stepFspSolver
END MODULE fspSolver
