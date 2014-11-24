MODULE IO

  USE sweeper

  IMPLICIT NONE

  PUBLIC :: processCmdLine
  PUBLIC :: closeFiles
  PUBLIC :: populateData

  INTEGER,PARAMETER :: inpFileUnitNo=123

  CONTAINS
!===============================================================================
    SUBROUTINE processCmdLine()
      ! Local Variables
      CHARACTER(LEN=16) :: arg_in

      CALL GET_COMMAND_ARGUMENT(1,arg_in)
      OPEN(FILE=TRIM(ADJUSTL(arg_in)),UNIT=inpFileUnitNo)

    END SUBROUTINE processCmdLine
!===============================================================================
    SUBROUTINE closeFiles()

      CLOSE(UNIT=inpFileUnitNo)

    END SUBROUTINE closeFiles
!===============================================================================
    SUBROUTINE populateData(mySweeper)
      TYPE(sweeperType),INTENT(INOUT) :: mySweeper
      ! Local variables
      INTEGER :: i,ii,iii,n1,n2,n3,n4

      ! Read number of groups and size of XS mesh
      READ(123,*) n1,n2,n3
      mySweeper%ngroups=n3-n2+1
      ALLOCATE(mySweeper%myXSMesh(n1))

      ! Read XSMeshType data
      DO i=1,n1
        READ(123,*) mySweeper%myXSMesh(i)%nreg
        ALLOCATE(mySweeper%myXSMesh(i)%ireg(mySweeper%myXSMesh(i)%nreg))
        READ(123,*) mySweeper%myXSMesh(i)%ireg
        ALLOCATE(mySweeper%myXSMesh(i)%xsmactr(mySweeper%ngroups))
        READ(123,*) mySweeper%myXSMesh(i)%xsmactr
        ALLOCATE(mySweeper%myXSMesh(i)%xsmacchi(mySweeper%ngroups))
        READ(123,*) mySweeper%myXSMesh(i)%xsmacchi
        ALLOCATE(mySweeper%myXSMesh(i)%xsmacsc(mySweeper%ngroups,0:0))
        DO ii=1,mySweeper%ngroups
          READ(123,*) mySweeper%myXSMesh(i)%xsmacsc(ii,0)%gmin
          READ(123,*) mySweeper%myXSMesh(i)%xsmacsc(ii,0)%gmax
          ALLOCATE(mySweeper%myXSMesh(i)%xsmacsc(ii,0)%from( &
            mySweeper%myXSMesh(i)%xsmacsc(ii,0)%gmin: &
            mySweeper%myXSMesh(i)%xsmacsc(ii,0)%gmax))
          READ(123,*) mySweeper%myXSMesh(i)%xsmacsc(ii,0)%from
        ENDDO
      ENDDO

      ! Read ModMeshType data
      ALLOCATE(ModMeshType :: mySweeper%myModMesh)
      READ(123,*) mySweeper%myModMesh%nmesh
      ALLOCATE(mySweeper%myModMesh%ifrstfsreg(mySweeper%myModMesh%nmesh))
      DO i=1,mySweeper%myModMesh%nmesh
        READ(123,*) mySweeper%myModMesh%ifrstfsreg(i)
      ENDDO
      READ(123,*) n1,n2
      ALLOCATE(mySweeper%myModMesh%neigh(n1,n2))
      READ(123,*) mySweeper%myModMesh%neigh

      ! Read ModMeshRayPtrArryType Data
      READ(123,*) n1
      ALLOCATE(mySweeper%rtmesh(n1))
      DO i=1,n1
        READ(123,*) n2
        ALLOCATE(mySweeper%rtmesh(i)%rtdat)
        ALLOCATE(mySweeper%rtmesh(i)%rtdat%angles(n2))
        DO ii=1,n2
          READ(123,*) n3
          ALLOCATE(mySweeper%rtmesh(i)%rtdat%angles(ii)%rays(n3))
          DO iii=1,n3
            READ(123,*) n4
            ALLOCATE(mySweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%ireg(n4))
            READ(123,*) mySweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%ireg
            READ(123,*) n4
            ALLOCATE(mySweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%hseg(n4))
            READ(123,*) mySweeper%rtmesh(i)%rtdat%angles(ii)%rays(iii)%hseg
          ENDDO
        ENDDO
      ENDDO

      ! Read CoreLongRays data
      READ(123,*) n1
      ALLOCATE(mySweeper%longRayDat%nlongrays(n1))
      ALLOCATE(mySweeper%longRayDat%angles(n1))
      DO i=1,n1
        READ(123,*) mySweeper%longRayDat%nlongrays(i)
        ALLOCATE(mySweeper%longRayDat%angles(i)%longrays( &
          mySweeper%longRayDat%nlongrays(i)))
        DO ii=1,mySweeper%longRayDat%nlongrays(i)
          READ(123,*) mySweeper%longRayDat%angles(i)%longrays(ii)%ifirstModMesh
          READ(123,*) mySweeper%longRayDat%angles(i)%longrays(ii)%iside
          READ(123,*) mySweeper%longRayDat%angles(i)%longrays(ii)%firstModRay
          READ(123,*) mySweeper%longRayDat%angles(i)%longrays(ii)%BCIndex
        ENDDO
      ENDDO

      ! Read ModularRayType data
      ALLOCATE(mySweeper%modRayDat)
      READ(123,*) mySweeper%modRayDat%iangstt,mySweeper%modRayDat%iangstp
      n1=mySweeper%modRayDat%iangstp-mySweeper%modRayDat%iangstt+1
      ALLOCATE(mySweeper%modRayDat%angles(n1))
      DO i=1,n1
        READ(123,*) mySweeper%modRayDat%angles(i)%dlr
        READ(123,*) n2
        mySweeper%modRayDat%angles(i)%nmodrays=n2
        ALLOCATE(mySweeper%modRayDat%angles(i)%rays(n2))
        DO ii=1,n2
          READ(123,*) mySweeper%modRayDat%angles(i)%rays(ii)%nextray
          READ(123,*) mySweeper%modRayDat%angles(i)%rays(ii)%nextsurf
        ENDDO
      ENDDO

      ! Read Quadrature information
      READ(123,*) n1
      ALLOCATE(mySweeper%modRayDat%angquad%walpha(n1))
      READ(123,*) mySweeper%modRayDat%angquad%walpha        
      READ(123,*) n1
      ALLOCATE(mySweeper%modRayDat%angquad%wtheta(n1))
      READ(123,*) mySweeper%modRayDat%angquad%wtheta
      READ(123,*) n1
      ALLOCATE(mySweeper%modRayDat%angquad%sinpolang(n1))
      READ(123,*) mySweeper%modRayDat%angquad%sinpolang       
      ALLOCATE(mySweeper%modRayDat%angquad%rsinpolang(n1))
      READ(123,*) mySweeper%modRayDat%angquad%rsinpolang       

      ! Read pz data
      READ(123,*) mySweeper%pz

      ! Read Exponential Tables data
      ALLOCATE(mySweeper%expTableDat)
      READ(123,*) n1
      READ(123,*) n2
      ALLOCATE(mySweeper%expTableDat%table2D(1:2,n1:n2))
      READ(123,*) mySweeper%expTableDat%table2D

    END SUBROUTINE populateData
END MODULE IO

