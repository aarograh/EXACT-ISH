MODULE fspSolver

  USE sweeper
  USE sweeperUtils
  USE openmp
  USE openacc
  USE IO

  IMPLICIT NONE
  PRIVATE

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
    SUBROUTINE initializeFspSolver(solver,sweepType)
      CLASS(fspSolverType),INTENT(INOUT) :: solver
      INTEGER,INTENT(IN) :: sweepType

      ALLOCATE(solver%sweeper)
      CALL populateData(solver%sweeper,solver%psi)

      CALL solver%sweeper%initialize(solver%source)

      SELECTCASE(sweepType)
        CASE(BASESOLVER)
          solver%sweeper%sweep2D_prodquad => sweep2D_prodquad_P0
        CASE(VECTORIPOL)
          solver%sweeper%sweep2D_prodquad => sweep2D_prodquad_P0_vecoripol
        CASE DEFAULT
          WRITE(*,*) 'Something went wrong when selecting the solver type.'
          STOP 666
      END SELECT

    END SUBROUTINE initializeFspSolver
!===============================================================================
    SUBROUTINE solveFspSolver(solver)
      CLASS(fspSolverType),INTENT(INOUT) :: solver

      CALL solver%step()

    END SUBROUTINE solveFspSolver
!===============================================================================
    SUBROUTINE stepFspSolver(solver)
      CLASS(fspSolverType),INTENT(INOUT) :: solver

      ! Perform sweep
      CALL solver%sweeper%sweep(1,1.0D-03,solver%source,solver%psi)

    END SUBROUTINE stepFspSolver
END MODULE fspSolver
