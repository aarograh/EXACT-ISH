PROGRAM EXACTISH

  USE sweeper
  USE IO

  TYPE(sweeperType) :: mySweeper

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Processing command line arguments...'
  WRITE(*,*) '======================================'
  CALL processCmdLine()

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Parsing data and setting up solvers...'
  WRITE(*,*) '======================================'
  CALL populateData(mySweeper)

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Performing transport sweeps...'
  WRITE(*,*) '======================================'

  WRITE(*,*) '======================================'
  WRITE(*,*) 'Closing files...'
  WRITE(*,*) '======================================'
  CALL closeFiles()

  WRITE(*,*) 'Hello World'

END PROGRAM EXACTISH
