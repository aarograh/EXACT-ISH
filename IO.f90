MODULE IO

  USE sweeper

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: processCmdLine
  PUBLIC :: closeFiles
  PUBLIC :: populateData
  PUBLIC :: validate
  PUBLIC :: BASESOLVER
  PUBLIC :: VECTORIPOL

  INTEGER,PARAMETER :: inpFileUnitNo=123
  INTEGER,PARAMETER :: solFileUnitNo=124
  INTEGER,PARAMETER :: outFileUnitNo=125

  INTEGER,PARAMETER :: BASESOLVER = 1
  INTEGER,PARAMETER :: VECTORIPOL = 2

  CONTAINS
!===============================================================================
    SUBROUTINE processCmdLine(sweepType)
      INTEGER,INTENT(INOUT) :: sweepType
      ! Local Variables
      CHARACTER(LEN=24) :: arg_in

      CALL GET_COMMAND_ARGUMENT(1,arg_in)
      OPEN(FILE=TRIM(ADJUSTL(arg_in))//'.dump',UNIT=inpFileUnitNo, &
        FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='OLD')
      OPEN(FILE=TRIM(ADJUSTL(arg_in))//'.sol',UNIT=solFileUnitNo, &
        FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='OLD')
      OPEN(FILE=TRIM(ADJUSTL(arg_in))//'.out',UNIT=outFileUnitNo)

      IF(COMMAND_ARGUMENT_COUNT() == 2) THEN
        arg_in = ''
        CALL GET_COMMAND_ARGUMENT(2,arg_in)
        READ(arg_in,*) sweepType
      ELSE
        sweepType = BASESOLVER
      ENDIF

    END SUBROUTINE processCmdLine
!===============================================================================
    SUBROUTINE closeFiles()

      CLOSE(UNIT=inpFileUnitNo)
      CLOSE(UNIT=solFileUnitNo)
      CLOSE(UNIT=outFileUnitNo)

    END SUBROUTINE closeFiles
!===============================================================================
    SUBROUTINE populateData(sweeper,psi)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      DOUBLE PRECISION,POINTER,INTENT(INOUT) :: psi(:)
      ! Local variables
      INTEGER :: i,ii,iii,n1,n2,n3,n4,n5

      ! Read number of groups and size of XS mesh
      READ(inpFileUnitNo) n1,n2,n3,n4
      sweeper%ng=n3-n2+1
      sweeper%igstt=n2
      sweeper%igstp=n3
      ALLOCATE(sweeper%myXSMesh(n1))

      ! Read XSMeshType data
      sweeper%nxsreg = n1
      DO i=1,n1
        READ(inpFileUnitNo) sweeper%myXSMesh(i)%nreg
        ALLOCATE(sweeper%myXSMesh(i)%ireg(sweeper%myXSMesh(i)%nreg))
        READ(inpFileUnitNo) sweeper%myXSMesh(i)%ireg
        ALLOCATE(sweeper%myXSMesh(i)%xsmactr(sweeper%ng))
        READ(inpFileUnitNo) sweeper%myXSMesh(i)%xsmactr
        READ(inpFileUnitNo) n2,n3
        IF(n2 > 0 .AND. n2 <= n3 .AND. n3 <= sweeper%ng) THEN
          ALLOCATE(sweeper%myXSMesh(i)%xsmacchi(sweeper%ng))
          READ(inpFileUnitNo) sweeper%myXSMesh(i)%xsmacchi
        ENDIF
        ALLOCATE(sweeper%myXSMesh(i)%xsmacsc(sweeper%ng,0:n4))
        DO ii=0,n4
          DO iii=1,sweeper%ng
            READ(inpFileUnitNo) sweeper%myXSMesh(i)%xsmacsc(iii,ii)%gmin
            READ(inpFileUnitNo) sweeper%myXSMesh(i)%xsmacsc(iii,ii)%gmax
            ALLOCATE(sweeper%myXSMesh(i)%xsmacsc(iii,ii)%from( &
              sweeper%myXSMesh(i)%xsmacsc(iii,ii)%gmin: &
              sweeper%myXSMesh(i)%xsmacsc(iii,ii)%gmax))
!IF(i == 1) WRITE(*,*) i,ii,iii,':',sweeper%myXSMesh(i)%xsmacsc(iii,ii)%gmin, &
!sweeper%myXSMesh(i)%xsmacsc(iii,ii)%gmax
            READ(inpFileUnitNo) sweeper%myXSMesh(i)%xsmacsc(iii,ii)%from
          ENDDO
        ENDDO
      ENDDO

      ! Read ModMeshType data
      ALLOCATE(sweeper%myModMesh)
      READ(inpFileUnitNo) sweeper%myModMesh%nmesh
      ALLOCATE(sweeper%myModMesh%ifrstfsreg(sweeper%myModMesh%nmesh))
      DO i=1,sweeper%myModMesh%nmesh
        READ(inpFileUnitNo) sweeper%myModMesh%ifrstfsreg(i)
      ENDDO
      READ(inpFileUnitNo) n1,n2
      ALLOCATE(sweeper%myModMesh%neigh(n1,n2))
      READ(inpFileUnitNo) sweeper%myModMesh%neigh

      ! Read psi
      READ(inpFileUnitNo) n1
      ALLOCATE(psi(n1))
      READ(inpFileUnitNo) psi

      ! Read ModMeshRayPtrArryType Data
      READ(inpFileUnitNo) n1
      ALLOCATE(sweeper%rtmesh(n1))
      DO i=1,n1
        READ(inpFileUnitNo) n2
        ALLOCATE(sweeper%rtmesh(i)%rtdat)
        ALLOCATE(sweeper%rtmesh(i)%rtdat%angles(n2))
        DO ii=1,n2
          READ(inpFileUnitNo) n3
          ALLOCATE(sweeper%rtmesh(i)%rtdat%angles(ii)%rays(n3))
          DO iii=1,n3
            READ(inpFileUnitNo) sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%nseg
            READ(inpFileUnitNo) n4
            ALLOCATE(sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%ireg(n4))
            READ(inpFileUnitNo) sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%ireg
            READ(inpFileUnitNo) n4
            ALLOCATE(sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%hseg(n4))
            READ(inpFileUnitNo) sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%hseg
          ENDDO
        ENDDO
      ENDDO

      ! Read CoreLongRays data
      READ(inpFileUnitNo) n1
      ALLOCATE(sweeper%longRayDat%nlongrays(n1))
      ALLOCATE(sweeper%longRayDat%angles(n1))
      DO i=1,n1
        READ(inpFileUnitNo) sweeper%longRayDat%nlongrays(i)
        ALLOCATE(sweeper%longRayDat%angles(i)%longrays( &
          sweeper%longRayDat%nlongrays(i)))
        DO ii=1,sweeper%longRayDat%nlongrays(i)
          READ(inpFileUnitNo) sweeper%longRayDat%angles(i)%longrays(ii)%nmods
          READ(inpFileUnitNo) sweeper%longRayDat%angles(i)%longrays(ii)%ifirstModMesh
          READ(inpFileUnitNo) sweeper%longRayDat%angles(i)%longrays(ii)%iside
          READ(inpFileUnitNo) sweeper%longRayDat%angles(i)%longrays(ii)%firstModRay
          READ(inpFileUnitNo) sweeper%longRayDat%angles(i)%longrays(ii)%BCIndex
        ENDDO
      ENDDO

      ! Read ModularRayType data
      ALLOCATE(sweeper%modRayDat)
      READ(inpFileUnitNo) sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
      n1=sweeper%modRayDat%iangstp-sweeper%modRayDat%iangstt+1
      ALLOCATE(sweeper%modRayDat%angles(n1))
      DO i=1,n1
        READ(inpFileUnitNo) sweeper%modRayDat%angles(i)%dlr
        READ(inpFileUnitNo) n2
        sweeper%modRayDat%angles(i)%nmodrays=n2
        ALLOCATE(sweeper%modRayDat%angles(i)%rays(n2))
        DO ii=1,n2
          READ(inpFileUnitNo) sweeper%modRayDat%angles(i)%rays(ii)%nextray
          READ(inpFileUnitNo) sweeper%modRayDat%angles(i)%rays(ii)%nextsurf
        ENDDO
      ENDDO

      ! Read Quadrature information
      READ(inpFileUnitNo) sweeper%modRayDat%angquad%npol
      READ(inpFileUnitNo) sweeper%modRayDat%angquad%nazi
      READ(inpFileUnitNo) n1
      ALLOCATE(sweeper%modRayDat%angquad%walpha(n1))
      READ(inpFileUnitNo) sweeper%modRayDat%angquad%walpha
      READ(inpFileUnitNo) n1
      ALLOCATE(sweeper%modRayDat%angquad%wtheta(n1))
      READ(inpFileUnitNo) sweeper%modRayDat%angquad%wtheta
      READ(inpFileUnitNo) n1
      ALLOCATE(sweeper%modRayDat%angquad%sinpolang(n1))
      READ(inpFileUnitNo) sweeper%modRayDat%angquad%sinpolang
      ALLOCATE(sweeper%modRayDat%angquad%rsinpolang(n1))
      READ(inpFileUnitNo) sweeper%modRayDat%angquad%rsinpolang

      ! Read pz data
      READ(inpFileUnitNo) sweeper%pz

      ! Read Exponential Tables data
      ALLOCATE(sweeper%expTableDat)
      READ(inpFileUnitNo) sweeper%expTableDat%minVal
      READ(inpFileUnitNo) sweeper%expTableDat%maxVal
      READ(inpFileUnitNo) n1
      READ(inpFileUnitNo) n2
      ALLOCATE(sweeper%expTableDat%table2D(1:2,n1:n2))
      READ(inpFileUnitNo) sweeper%expTableDat%table2D
      READ(inpFileUnitNo) sweeper%expTableDat%rdx

      ! Read AngFluxBC data
      READ(inpFileUnitNo) n1
      ALLOCATE(sweeper%phiang(n1))
      DO i=1,n1
        READ(inpFileUnitNo) n2
        ALLOCATE(sweeper%phiang(i)%angle(n2))
        IF(i == 1) ALLOCATE(sweeper%phiang1g_out%angle(n2))
        DO ii=1,n2
          READ(inpFileUnitNo) n3
          ALLOCATE(sweeper%phiang(i)%angle(ii)%face(n3))
          IF(i == 1) ALLOCATE(sweeper%phiang1g_out%angle(ii)%face(n3))
          DO iii=1,n3
            READ(inpFileUnitNo) n4,n5
            ALLOCATE(sweeper%phiang(i)%angle(ii)%face(iii)%angFlux(n4,0:n5-1))
            READ(inpFileUnitNo) sweeper%phiang(i)%angle(ii)%face(iii)%angFlux
            IF(i == 1) &
              ALLOCATE(sweeper%phiang1g_out%angle(ii)%face(iii)%angFlux(n4,0:n5-1))
          ENDDO
        ENDDO
      ENDDO

      ! Read UpdateBC data
      READ(inpFileUnitNo) sweeper%UpdateBC%offset
      READ(inpFileUnitNo) sweeper%UpdateBC%nfaces
      READ(inpFileUnitNo) sweeper%UpdateBC%bcType
      READ(inpFileUnitNo) sweeper%UpdateBC%nangles
      READ(inpFileUnitNo) sweeper%UpdateBC%iangstt
      READ(inpFileUnitNo) sweeper%UpdateBC%iangstp
      READ(inpFileUnitNo) n1,n2
      ALLOCATE(sweeper%UpdateBC%iang2irefl(n1,n2))
      READ(inpFileUnitNo) sweeper%UpdateBC%iang2irefl

      ! Read in initial phis
      READ(inpFileUnitNo) n1,n2
      sweeper%nreg = n1
      ALLOCATE(sweeper%phis(n1,n2))
      READ(inpFileUnitNo) sweeper%phis

      ! Miscellaneous
      READ(inpFileUnitNo) sweeper%maxsegray
      READ(inpFileUnitNo) sweeper%myModMesh%nmesh,sweeper%imeshstt
      ALLOCATE(sweeper%vol(sweeper%nreg))
      READ(inpFileUnitNo) sweeper%vol

    END SUBROUTINE populateData
!===============================================================================
    SUBROUTINE validate(sweeper)
      CLASS(sweeperType),INTENT(IN) :: sweeper
      ! Local Variables
      INTEGER :: nreg,ngroups,ireg,ig
      DOUBLE PRECISION :: compval,diff,maxdiff,rmsdiff

      READ(solFileUnitNo) nreg
      READ(solFileUnitNo) ngroups

      maxdiff = 0.0D0
      rmsdiff = 0.0D0
      DO ig=1,ngroups
        DO ireg=1,nreg
          READ(solFileUnitNo) compval
          diff = (sweeper%phis(ireg,ig) - compval)/compval
          maxdiff = MAX(maxdiff,ABS(diff))
          rmsdiff = rmsdiff + diff*diff
        ENDDO !ireg
      ENDDO !ig

      rmsdiff = SQRT(rmsdiff)

      WRITE(*,*)
      WRITE(*,*) 'RMS Difference = ',rmsdiff
      WRITE(*,*) 'Max Difference = ',maxdiff

    END SUBROUTINE validate
END MODULE IO

