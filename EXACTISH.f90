PROGRAM EXACTISH

  USE IO
  USE fspSolver

  IMPLICIT NONE

  INTEGER :: solveType=0
  TYPE(fspSolverType) :: solver

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Processing command line arguments...'
  WRITE(*,*) '======================================'
  CALL processCmdLine(solveType)

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Initializing solvers...'
  WRITE(*,*) '======================================'
  CALL solver%initialize(solveType)

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Performing transport sweeps...'
  WRITE(*,*) '======================================'
  CALL solver%solve()

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Validating Solution...'
  WRITE(*,*) '======================================'
  IF(solveType < PGIBASESOLVER) THEN
    CALL validate(solver%sweeper%phis)
  ELSE
    CALL validate(solver%PGIsweeper%phis)
  ENDIF

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Closing files...'
  WRITE(*,*) '======================================'
  CALL closeFiles()

END PROGRAM EXACTISH
