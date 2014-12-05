PROGRAM jipu_main

  USE IO
  USE fspSolver

  IMPLICIT NONE

  TYPE(fspSolverType) :: solver

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Processing command line arguments...'
  WRITE(*,*) '======================================'
  CALL processCmdLine()

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Initializing solvers...'
  WRITE(*,*) '======================================'
  CALL solver%initialize(BASESOLVER)

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Performing transport sweeps...'
  WRITE(*,*) '======================================'
  CALL solver%solve()

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Validating Solution...'
  WRITE(*,*) '======================================'
  CALL validate(solver%sweeper)

  WRITE(*,*)
  WRITE(*,*) '======================================'
  WRITE(*,*) 'Closing files...'
  WRITE(*,*) '======================================'
  CALL closeFiles()

END PROGRAM jipu_main
