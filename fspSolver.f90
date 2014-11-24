MODULE fspSolver

  USE sweeper
  USE IO

  IMPLICIT NONE

  PUBLIC :: fspSolverType

  TYPE :: fspSolverType
    DOUBLE PRECISION,POINTER :: psi(:) => NULL()
    CLASS(SourceType),POINTER :: thisSource
    CLASS(sweeperType),POINTER :: sweeper => NULL()
    CONTAINS
      PROCEDURE,PASS :: initialize => initializeFspSolver
      PROCEDURE,PASS :: solve => solveFspSolver
      PROCEDURE,PASS :: step => stepFspSolver
  END TYPE fspSolverType

  CONTAINS
!===============================================================================
    SUBROUTINE initializeFspSolver(solver)
      CLASS(fspSolverType),INTENT(INOUT) :: solver

      ALLOCATE(solver%sweeper)
      CALL populateData(solver%sweeper)

      CALL solver%sweeper%initialize(solver%thisSource)

    END SUBROUTINE initializeFspSolver
!===============================================================================
    SUBROUTINE solveFspSolver(solver)
      CLASS(fspSolverType),INTENT(INOUT) :: solver

      CALL solver%step()

    END SUBROUTINE solveFspSolver
!===============================================================================
    SUBROUTINE stepFspSolver(solver)
      CLASS(fspSolverType),INTENT(INOUT) :: solver
      ! Local Variables
      INTEGER :: ig

      DO ig=solver%sweeper%igstt,solver%sweeper%igstp
        ! Set up source
        CALL solver%thisSource%initExtSource(ig)
        CALL solver%thisSource%computeMGFS(ig,solver%psi)
        CALL solver%thisSource%updateInScatter( &
          ig,solver%sweeper%igstt,solver%sweeper%igstp)
        CALL solver%sweeper%setExtSource(solver%thisSource)
        ! Perform sweep
        CALL solver%sweeper%sweep(ig,3,1.0D-03)
      ENDDO

    END SUBROUTINE stepFspSolver
END MODULE fspSolver
