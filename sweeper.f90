MODULE sweeper
!TODO: CALL sweep2D_prodquad
!TODO: mimic sweep loops
!TODO: Add timers for these loops

  USE sweeperUtils

  IMPLICIT NONE

  PUBLIC :: sweeperType

  TYPE :: sweeperType
    LOGICAL :: hasSource=.FALSE.
    INTEGER :: ninners=0
    INTEGER :: ngroups=0
    INTEGER :: igstt=0
    INTEGER :: igstp=0
    INTEGER :: activeg=0
    INTEGER :: nsweeps=0
    INTEGER :: nxsreg=0
    INTEGER :: nreg=0
    DOUBLE PRECISION :: pz=0.0D0
    DOUBLE PRECISION,ALLOCATABLE :: phis1g(:)
    DOUBLE PRECISION,ALLOCATABLE :: phis1gd(:)
    DOUBLE PRECISION,POINTER :: phis(:,:)
    DOUBLE PRECISION,POINTER :: xstr(:) => NULL()
    DOUBLE PRECISION,POINTER :: vol(:) => NULL()
    DOUBLE PRECISION,POINTER :: qbar(:) => NULL()
    DOUBLE PRECISION,POINTER :: qext(:) => NULL()
    TYPE(AngFluxBC),POINTER :: phiang1g_in => NULL() 
    TYPE(AngFluxBC) :: phiang1g_out
    TYPE(AngFluxBC),POINTER :: phiang(:) 
    TYPE(SourceType_P0),POINTER :: mySrc => NULL()
    TYPE(ModMeshType),POINTER :: myModMesh => NULL() 
    TYPE(ModularRayType),POINTER :: modRayDat
    TYPE(CoreLongRayType) :: longRayDat 
    TYPE(ModMeshRayPtrArryType),POINTER :: rtmesh(:) => NULL()
    TYPE(XSMeshType),POINTER :: myXSMesh(:)
    TYPE(ExpTableType),ALLOCATABLE :: expTableDat
!    TYPE(UpdateBCType_MOC) :: updateBC !maybe, UpdateBC_MOC.f90: define this, %Start() and %Finish() methods.  Might be an MPI thing that I don't need
    PROCEDURE(absintfc_sweep),POINTER :: sweep => NULL()
    PROCEDURE(absintfc_setExtSource),POINTER :: setExtSource => NULL()
    CONTAINS
      PROCEDURE,PASS :: initialize => initializeSweeper
!      PROCEDURE,PASS :: calcFissionSrc
  END TYPE sweeperType

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_sweep(sweeper,igroup,ninners,tol)
      IMPORT sweeperType
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: igroup
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
    END SUBROUTINE absintfc_sweep
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_setExtSource(thisTS,source)
      IMPORT sweeperType,SourceType
      CLASS(sweeperType),INTENT(INOUT) :: thisTS
      CLASS(SourceType),POINTER,INTENT(IN) :: source
    END SUBROUTINE absintfc_setExtSource
  END INTERFACE

  CONTAINS
!===============================================================================
    SUBROUTINE initializeSweeper(thisSweeper,source)
      CLASS(sweeperType),INTENT(INOUT) :: thisSweeper
      CLASS(sourceType),POINTER,INTENT(INOUT) :: source

      ! Set up source
      ALLOCATE(SourceType_P0 :: source)
      SELECTTYPE(source); TYPE IS(SourceType_P0)
        thisSweeper%mySrc => source
      END SELECT
      source%nreg = thisSweeper%nreg
      source%nxsreg = thisSweeper%nxsreg
      source%ng = thisSweeper%ngroups
      source%phis => thisSweeper%phis
      source%myXSMesh => thisSweeper%myXSMesh
      source%qext => thisSweeper%qext

    END SUBROUTINE initializeSweeper
!===============================================================================
!    SUBROUTINE calcFissionSrc(sweeper,psi,keff)
!      CLASS(sweeperType),INTENT(INOUT) :: sweeper
!      DOUBLE PRECISION,INTENT(INOUT) :: psi
!      DOUBLE PRECISION,OPTIONAL,INTENT(IN) :: keff
!      ! Local Variables
!      INTEGER :: ig,ix,ir,ireg
!      DOUBLE PRECISION :: xsnfg,rkeff
!
!      rkeff = 1.0D0
!      IF(PRESENT(keff)) rkeff = 1.0D0/keff
!
!      psi = 0.0D0
!      DO ig=1,sweeper%ngroups
!        DO ix=1,sweeper%nxsreg
!          IF(LBOUND(sweeper%myXSMesh(ix)%xsmacnf) <= ig .AND. &
!            ig <= UBOUND(sweeper%myXSMesh(ix)%xsmacnf)) THEN
!            xsnfg = sweeper%myXSMesh(ix)%xsmacnf(ig)*rkeff
!            DO ir=1,sweeper%myXSMesh(ix)%nreg
!              ireg = sweeper%myXSMesh(ix)%ireg(ir)
!              psi(ireg) = psi(ireg) + xsnfg*sweeper%phis(ireg,ig)
!            ENDDO
!          endif
!          
!
!    END SUBROUTINE calcFissionSrc
!===============================================================================
    SUBROUTINE setExtSource_MOCP0(thisTS,source)
      CLASS(sweeperType),INTENT(INOUT) :: thisTS
      CLASS(SourceType),POINTER,INTENT(IN) :: source

      NULLIFY(thisTS%qext)
      thisTS%hasSource=.FALSE.
      SELECTTYPE(source); TYPE IS(SourceType_P0)
        thisTS%mySrc => source
        thisTS%qext => source%qext
        thisTS%hasSource=.TRUE.
      ENDSELECT

    END SUBROUTINE setExtSource_MOCP0
!===============================================================================
    SUBROUTINE MOCSolver_Sweep1G(sweeper,igroup,ninners,tol)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: igroup
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
    END SUBROUTINE MOCSolver_Sweep1G
!===============================================================================
    SUBROUTINE MOCSOlver_Setup1GFSP(thisXSMesh,nxsreg,phis1g,nreg,xstr,qbar,ig)
      INTEGER,INTENT(IN) :: nxsreg
      TYPE(XSMeshType),INTENT(IN) :: thisXSMesh(nxsreg)
      INTEGER,INTENT(IN) :: nreg
      INTEGER,INTENT(IN) :: ig
      DOUBLE PRECISION,INTENT(IN) :: phis1g(nreg)
      DOUBLE PRECISION,INTENT(INOUT) :: xstr(nreg)
      DOUBLE PRECISION,INTENT(INOUT) :: qbar(nreg)
      ! Local Variables
      INTEGER :: ix,ir,ireg
      DOUBLE PRECISION :: xstrg

      DO ix=1,nxsreg
        xstrg=thisXSMesh(ix)%xsmactr(ig)
        DO ir=1,thisXSMesh(ix)%nreg
          ireg=thisXSMesh(ix)%ireg(ir)
          xstr(ireg)=xstrg
        ENDDO
      ENDDO

    END SUBROUTINE MOCSolver_Setup1GFSP
END MODULE sweeper
