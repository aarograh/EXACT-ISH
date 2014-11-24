PROGRAM EXACTISH

  USE sweeper
  USE fspSolver
  USE IO

  IMPLICIT NONE

  TYPE(sweeperType),TARGET :: mySweeper
  TYPE(fspSolverType) :: mySolver

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Processing command line arguments...'
  WRITE(*,*) '======================================'
  CALL processCmdLine()

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Parsing data and setting up solvers...'
  WRITE(*,*) '======================================'
  CALL populateData(mySweeper)

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Initializing solvers...'
  WRITE(*,*) '======================================'
  CALL mySolver%initialize(mySweeper)

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Performing transport sweeps...'
  WRITE(*,*) '======================================'
  CALL mySolver%solve()

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Closing files...'
  WRITE(*,*) '======================================'
  CALL closeFiles()

  WRITE(*,*) 'Hello World'

END PROGRAM EXACTISH
