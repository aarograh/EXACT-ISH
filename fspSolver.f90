MODULE fspSolver

  USE sweeper
  USE openmp
  USE openacc
  USE IO

  IMPLICIT NONE

  PUBLIC :: fspSolverType

  TYPE :: fspSolverType
    DOUBLE PRECISION,POINTER :: psi(:) => NULL()
    CLASS(SourceType),POINTER :: source
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
      CALL populateData(solver%sweeper,solver%psi)

      CALL solver%sweeper%initialize(solver%source)

      !TODO: put branching statements here to associate %sweep
      ! with different kernels.
      ! Branching should be based on something passed into
      ! initializeFspSolver from main_aaron or main_jipu
      solver%sweeper%sweep2D_prodquad => sweep2D_prodquad_P0

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
        CALL solver%source%initExtSource(ig)
        CALL solver%source%computeMGFS(ig,solver%psi)
        CALL solver%source%updateInScatter( &
          ig,solver%sweeper%igstt,solver%sweeper%igstp)
        CALL solver%sweeper%setExtSource(solver%source)
        ! Perform sweep
        CALL solver%sweeper%sweep(ig,3,1.0D-03)
      ENDDO

    END SUBROUTINE stepFspSolver
END MODULE fspSolver
