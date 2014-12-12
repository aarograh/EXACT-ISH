MODULE fspSolver

  USE sweeper
  USE sweeperUtils
  USE openmp
  USE PGIutils
  USE openacc
  USE IO

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: fspSolverType

  TYPE :: fspSolverType
    INTEGER :: sweepType = 0
    DOUBLE PRECISION,POINTER :: psi(:) => NULL()
    CLASS(SourceType),POINTER :: source => NULL()
    CLASS(sweeperType),POINTER :: sweeper => NULL()
    CLASS(SourceType_PGI),POINTER :: PGIsource => NULL()
    CLASS(sweeperType_PGI),POINTER :: PGIsweeper => NULL()
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

      IF(sweepType < PGIBASESOLVER) THEN
        ALLOCATE(solver%sweeper)
        CALL populateData(solver%sweeper,solver%psi)
        CALL solver%sweeper%initialize(solver%source)
      ELSE
        ALLOCATE(solver%PGIsweeper)
        CALL populateData_PGI(solver%PGIsweeper,solver%psi)
        CALL solver%PGIsweeper%initialize(solver%PGIsource)
      ENDIF

      SELECTCASE(sweepType)
        CASE(BASESOLVER)
          solver%sweeper%sweep2D_prodquad => sweep2D_prodquad_P0
        CASE(VECTORIPOL)
          solver%sweeper%sweep2D_prodquad => sweep2D_prodquad_P0_vectoripol1_2
        CASE(PGIBASESOLVER)
!          solver%PGIsweeper%sweep2D_prodquad => sweep2D_prodquad_P0_PGI
          CONTINUE !PGI has no clue what procedure pointers are
        CASE DEFAULT
          WRITE(*,*) 'Something went wrong when selecting the solver type.'
          STOP 666
      END SELECT

      solver%sweepType = sweepType

    END SUBROUTINE initializeFspSolver
!===============================================================================
    SUBROUTINE solveFspSolver(solver)
      CLASS(fspSolverType),INTENT(INOUT) :: solver

      CALL solver%step()

    END SUBROUTINE solveFspSolver
!===============================================================================
    SUBROUTINE stepFspSolver(solver)
      CLASS(fspSolverType),INTENT(INOUT) :: solver
      CLASS(sweeperType_PGI),POINTER :: ptr

ptr => solver%PGISweeper
      ! Perform sweep
      IF(solver%sweepType < PGIBASESOLVER) THEN
        CALL solver%sweeper%sweep(1,1.0D-03,solver%source,solver%psi)
      ELSE
!        CALL solver%PGIsweeper%sweep(1,1.0D-03,solver%PGIsource,solver%psi)
! PGI can't handle procedure pointers without crying like a small child...
        CALL ptr%sweep(1,1.0D-03,solver%PGIsource,solver%psi)
      ENDIF

    END SUBROUTINE stepFspSolver
END MODULE fspSolver
