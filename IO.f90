MODULE IO

  USE sweeper

  IMPLICIT NONE

  PUBLIC :: processCmdLine
  PUBLIC :: closeFiles
  PUBLIC :: populateData

  INTEGER,PARAMETER :: inpFileUnitNo=123
  INTEGER,PARAMETER :: solFileUnitNo=124

  CONTAINS
!===============================================================================
    SUBROUTINE processCmdLine()
      ! Local Variables
      CHARACTER(LEN=16) :: arg_in

      CALL GET_COMMAND_ARGUMENT(1,arg_in)
      OPEN(FILE=TRIM(ADJUSTL(arg_in))//'.dump',UNIT=inpFileUnitNo)
      OPEN(FILE=TRIM(ADJUSTL(arg_in))//'.sol',UNIT=solFileUnitNo)

    END SUBROUTINE processCmdLine
!===============================================================================
    SUBROUTINE closeFiles()

      CLOSE(UNIT=inpFileUnitNo)
      CLOSE(UNIT=solFileUnitNo)

    END SUBROUTINE closeFiles
!===============================================================================
    SUBROUTINE populateData(sweeper)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      ! Local variables
      INTEGER :: i,ii,iii,n1,n2,n3,n4,n5

      ! Read number of groups and size of XS mesh
      READ(123,*) n1,n2,n3
      sweeper%ng=n3-n2+1
      sweeper%igstt=n2
      sweeper%igstp=n3
      ALLOCATE(sweeper%myXSMesh(n1))

      ! Read XSMeshType data
      DO i=1,n1
        READ(123,*) sweeper%myXSMesh(i)%nreg
        ALLOCATE(sweeper%myXSMesh(i)%ireg(sweeper%myXSMesh(i)%nreg))
        READ(123,*) sweeper%myXSMesh(i)%ireg
        ALLOCATE(sweeper%myXSMesh(i)%xsmactr(sweeper%ng))
        READ(123,*) sweeper%myXSMesh(i)%xsmactr
        ALLOCATE(sweeper%myXSMesh(i)%xsmacchi(sweeper%ng))
        READ(123,*) sweeper%myXSMesh(i)%xsmacchi
        ALLOCATE(sweeper%myXSMesh(i)%xsmacsc(sweeper%ng,0:0))
        DO ii=1,sweeper%ng
          READ(123,*) sweeper%myXSMesh(i)%xsmacsc(ii,0)%gmin
          READ(123,*) sweeper%myXSMesh(i)%xsmacsc(ii,0)%gmax
          ALLOCATE(sweeper%myXSMesh(i)%xsmacsc(ii,0)%from( &
            sweeper%myXSMesh(i)%xsmacsc(ii,0)%gmin: &
            sweeper%myXSMesh(i)%xsmacsc(ii,0)%gmax))
          READ(123,*) sweeper%myXSMesh(i)%xsmacsc(ii,0)%from
        ENDDO
      ENDDO

      ! Read ModMeshType data
      ALLOCATE(sweeper%myModMesh)
      READ(123,*) sweeper%myModMesh%nmesh
      ALLOCATE(sweeper%myModMesh%ifrstfsreg(sweeper%myModMesh%nmesh))
      DO i=1,sweeper%myModMesh%nmesh
        READ(123,*) sweeper%myModMesh%ifrstfsreg(i)
      ENDDO
      READ(123,*) n1,n2
      ALLOCATE(sweeper%myModMesh%neigh(n1,n2))
      READ(123,*) sweeper%myModMesh%neigh

      ! Read ModMeshRayPtrArryType Data
      READ(123,*) n1
      ALLOCATE(sweeper%rtmesh(n1))
      DO i=1,n1
        READ(123,*) n2
        ALLOCATE(sweeper%rtmesh(i)%rtdat)
        ALLOCATE(sweeper%rtmesh(i)%rtdat%angles(n2))
        DO ii=1,n2
          READ(123,*) n3
          ALLOCATE(sweeper%rtmesh(i)%rtdat%angles(ii)%rays(n3))
          DO iii=1,n3
            READ(123,*) sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%nseg
            READ(123,*) n4
            ALLOCATE(sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%ireg(n4))
            READ(123,*) sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%ireg
            READ(123,*) n4
            ALLOCATE(sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%hseg(n4))
            READ(123,*) sweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%hseg
          ENDDO
        ENDDO
      ENDDO

      ! Read CoreLongRays data
      READ(123,*) n1
      ALLOCATE(sweeper%longRayDat%nlongrays(n1))
      ALLOCATE(sweeper%longRayDat%angles(n1))
      DO i=1,n1
        READ(123,*) sweeper%longRayDat%nlongrays(i)
        ALLOCATE(sweeper%longRayDat%angles(i)%longrays( &
          sweeper%longRayDat%nlongrays(i)))
        DO ii=1,sweeper%longRayDat%nlongrays(i)
          READ(123,*) sweeper%longRayDat%angles(i)%longrays(ii)%nmods
          READ(123,*) sweeper%longRayDat%angles(i)%longrays(ii)%ifirstModMesh
          READ(123,*) sweeper%longRayDat%angles(i)%longrays(ii)%iside
          READ(123,*) sweeper%longRayDat%angles(i)%longrays(ii)%firstModRay
          READ(123,*) sweeper%longRayDat%angles(i)%longrays(ii)%BCIndex
        ENDDO
      ENDDO

      ! Read ModularRayType data
      ALLOCATE(sweeper%modRayDat)
      READ(123,*) sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
      n1=sweeper%modRayDat%iangstp-sweeper%modRayDat%iangstt+1
      ALLOCATE(sweeper%modRayDat%angles(n1))
      DO i=1,n1
        READ(123,*) sweeper%modRayDat%angles(i)%dlr
        READ(123,*) n2
        sweeper%modRayDat%angles(i)%nmodrays=n2
        ALLOCATE(sweeper%modRayDat%angles(i)%rays(n2))
        DO ii=1,n2
          READ(123,*) sweeper%modRayDat%angles(i)%rays(ii)%nextray
          READ(123,*) sweeper%modRayDat%angles(i)%rays(ii)%nextsurf
        ENDDO
      ENDDO

      ! Read Quadrature information
      READ(123,*) sweeper%modRayDat%angquad%npol
      READ(123,*) sweeper%modRayDat%angquad%nazi
      READ(123,*) n1
      ALLOCATE(sweeper%modRayDat%angquad%walpha(n1))
      READ(123,*) sweeper%modRayDat%angquad%walpha        
      READ(123,*) n1
      ALLOCATE(sweeper%modRayDat%angquad%wtheta(n1))
      READ(123,*) sweeper%modRayDat%angquad%wtheta
      READ(123,*) n1
      ALLOCATE(sweeper%modRayDat%angquad%sinpolang(n1))
      READ(123,*) sweeper%modRayDat%angquad%sinpolang       
      ALLOCATE(sweeper%modRayDat%angquad%rsinpolang(n1))
      READ(123,*) sweeper%modRayDat%angquad%rsinpolang       

      ! Read pz data
      READ(123,*) sweeper%pz

      ! Read Exponential Tables data
      ALLOCATE(sweeper%expTableDat)
      READ(123,*) n1
      READ(123,*) n2
      ALLOCATE(sweeper%expTableDat%table2D(1:2,n1:n2))
      READ(123,*) sweeper%expTableDat%table2D
      READ(123,*) sweeper%expTableDat%rdx

      ! Read AngFluxBC data
      READ(123,*) n1
      ALLOCATE(sweeper%phiang(n1))
      DO i=1,n1
        READ(123,*) n2
        ALLOCATE(sweeper%phiang(i)%angle(n2))
        IF(i == 1) ALLOCATE(sweeper%phiang1g_out%angle(n2))
        DO ii=1,n2
          READ(123,*) n3
          ALLOCATE(sweeper%phiang(i)%angle(ii)%face(n3))
          IF(i == 1) ALLOCATE(sweeper%phiang1g_out%angle(ii)%face(n3))
          DO iii=1,n3
            READ(123,*) n4,n5
            ALLOCATE(sweeper%phiang(i)%angle(ii)%face(iii)%angFlux(n4,n5))
            READ(123,*) sweeper%phiang(i)%angle(ii)%face(iii)%angFlux
            IF(i == 1) &
              ALLOCATE(sweeper%phiang1g_out%angle(ii)%face(iii)%angFlux(n4,n5))
          ENDDO
        ENDDO
      ENDDO

      ! Read in initial phis
      READ(123,*) n1,n2
      sweeper%nreg = n1
      ALLOCATE(sweeper%phis(n1,n2))
      READ(123,*) sweeper%phis

      ! Miscellaneous
      READ(123,*) sweeper%maxsegray
      READ(123,*) sweeper%myModMesh%nmesh,sweeper%imeshstt
      ALLOCATE(sweeper%vol(sweeper%nreg))
      READ(123,*) sweeper%vol

    END SUBROUTINE populateData
END MODULE IO

