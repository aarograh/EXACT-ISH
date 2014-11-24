MODULE sweeperUtils

  IMPLICIT NONE

  PUBLIC :: AngFluxBC
  PUBLIC :: SourceType
  PUBLIC :: SourceType_P0
  PUBLIC :: ModularRayType
  PUBLIC :: CoreLongRayType
  PUBLIC :: ModMeshRayPtrArryType
  PUBLIC :: XSMeshType
  PUBLIC :: ExpTableType

  TYPE :: AngFluxBCFace
    SEQUENCE
    DOUBLE PRECISION,ALLOCATABLE :: angflux(:,:)
  END TYPE AngFluxBCFace

  TYPE :: AngFluxBCAng
    SEQUENCE
    TYPE(AngFluxBCFace),ALLOCATABLE :: face(:)
  END TYPE AngFluxBCAng

  TYPE :: AngFluxBC
    SEQUENCE
    TYPE(AngFluxBCAng),ALLOCATABLE :: angle(:)
  END TYPE AngFluxBC

  TYPE,ABSTRACT :: SourceType
    INTEGER :: nreg=0
    INTEGER :: nxsreg=0
    INTEGER :: ng=0
    DOUBLE PRECISION,POINTER :: phis(:,:) => NULL()
    DOUBLE PRECISION,POINTER :: qext(:) => NULL()
    DOUBLE PRECISION,POINTER :: split(:) => NULL()
    DOUBLE PRECISION,POINTER :: qextmg(:,:) => NULL()
    TYPE(XSMeshType),POINTER :: myXSMesh(:) => NULL()
    CONTAINS
      PROCEDURE(absintfc_initExtSource),PASS,DEFERRED :: initExtSource
      PROCEDURE(absintfc_computeMGFS),PASS,DEFERRED :: computeMGFS
      PROCEDURE(absintfc_updateInScatter),PASS,DEFERRED :: updateInScatter
  END TYPE SourceType

  TYPE,EXTENDS(SourceType) :: SourceType_P0
    DOUBLE PRECISION,POINTER :: qi1g(:) => NULL()
    CONTAINS
      PROCEDURE,PASS :: updateSelfScatter_P0
      PROCEDURE,PASS :: initExtSource => initExtSource_P0
      PROCEDURE,PASS :: computeMGFS => computeMGFS_P0
      PROCEDURE,PASS :: updateInScatter => updateInScatter_P0
  END TYPE SourceType_P0

  TYPE :: ModMeshType
    INTEGER :: nmesh
    INTEGER,ALLOCATABLE :: ifrstfsreg(:)
    INTEGER,ALLOCATABLE :: neigh(:,:)
  END TYPE ModMeshType

  TYPE :: ModRayLineType
    INTEGER :: nextray(2)=0
    INTEGER :: nextsurf(2)=0
  END TYPE ModRayLineType

  TYPE :: ModAngRayType
    DOUBLE PRECISION :: dlr=0.0D0
    INTEGER :: nmodrays
    TYPE(ModRayLineType),ALLOCATABLE :: rays(:)
  END TYPE ModAngRayType

  TYPE :: AngQuadType
    DOUBLE PRECISION,ALLOCATABLE :: walpha(:)
    DOUBLE PRECISION,ALLOCATABLE :: wtheta(:)
    DOUBLE PRECISION,ALLOCATABLE :: sinpolang(:)
    DOUBLE PRECISION,ALLOCATABLE :: rsinpolang(:)
  END TYPE AngQuadType

  TYPE :: ModularRayType
    INTEGER :: iangstt
    INTEGER :: iangstp
    TYPE(ModAngRayType),ALLOCATABLE :: angles(:)
    TYPE(AngQuadType) :: angquad
  END TYPE ModularRayType

  TYPE :: LongRayType_Base
    SEQUENCE
    INTEGER :: ifirstModMesh=0
    INTEGER :: iside(2)=0
    INTEGER :: firstModRay=0
    INTEGER :: BCIndex(2)=0
  END TYPE LongRayType_Base

  TYPE :: AngLongRayType
    SEQUENCE
    TYPE(LongRayType_Base),ALLOCATABLE :: longrays(:)
  END TYPE AngLongRayType

  TYPE :: CoreLongRayType
    INTEGER,ALLOCATABLE :: nlongrays(:)
    TYPE(AngLongRayType),ALLOCATABLE :: angles(:)
  END TYPE CoreLongRayType

  TYPE :: RayType
    INTEGER,ALLOCATABLE :: ireg(:)
    DOUBLE PRECISION,ALLOCATABLE :: hseg(:)
  END TYPE RayType

  TYPE :: AngleRayType
    TYPE(RayType),ALLOCATABLE :: rays(:)
  END TYPE AngleRayType

  TYPE :: ModMeshRayType
    TYPE(AngleRayType),ALLOCATABLE :: angles(:)
  END TYPE ModMeshRayType

  TYPE :: ModMeshRayPtrArryType
    TYPE(ModMeshRayType),POINTER :: rtdat => NULL()
  END TYPE ModMeshRayPtrArryType

  TYPE :: ScatMatType
    SEQUENCE
    INTEGER :: gmin
    INTEGER :: gmax
    DOUBLE PRECISION,ALLOCATABLE :: from(:)
  END TYPE ScatMatType

  TYPE :: XSMeshType
    INTEGER :: nreg=0
    INTEGER,ALLOCATABLE :: ireg(:)
    DOUBLE PRECISION,ALLOCATABLE :: xsmactr(:)
    DOUBLE PRECISION,ALLOCATABLE :: xsmacchi(:)
    TYPE(ScatMatType),ALLOCATABLE :: xsmacsc(:,:)
  END TYPE XSMeshType

  TYPE :: ExpTableType
    DOUBLE PRECISION,ALLOCATABLE :: table2D(:,:)
  END TYPE ExpTableType

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_initExtSource(thisSrc,ig)
      IMPORT :: Sourcetype
      CLASS(SourceType),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig
    END SUBROUTINE absintfc_initExtSource 
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_computeMGFS(thisSrc,ig,psi)
      IMPORT SourceType
      CLASS(SourceType),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig
      DOUBLE PRECISION,INTENT(IN) :: psi(:)
    END SUBROUTINE absintfc_computeMGFS
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_updateInScatter(thisSrc,ig,igstt,igstp)
      IMPORT SourceType
      CLASS(SourceType),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig
      INTEGER,INTENT(IN) :: igstt
      INTEGER,INTENT(IN) :: igstp
    END SUBROUTINE absintfc_updateInScatter
  END INTERFACE

  CONTAINS
!===============================================================================
!TODO: MOCSweeper_P0.f90:47
    SUBROUTINE updateSelfScatter_P0(thisSrc,igroup,qbar,phis1g)
      CLASS(SourceType_P0),INTENT(IN) :: thisSrc
      INTEGER,INTENT(IN) :: igroup
      DOUBLE PRECISION,INTENT(INOUT) :: qbar(:)
      DOUBLE PRECISION,INTENT(IN) :: phis1g
    END SUBROUTINE updateSelfScatter_P0
!===============================================================================
    SUBROUTINE initExtSource_P0(thisSrc,ig)
      CLASS(SourceType_P0),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig

      thisSrc%qi1g=thisSrc%qextmg(:,ig)

    END SUBROUTINE initExtSource_P0
!===============================================================================
    SUBROUTINE computeMGFS_P0(thisSrc,ig,psi)
      CLASS(SourceType_P0),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig
      DOUBLE PRECISION,INTENT(IN) :: psi(:)
      ! Local variables
      INTEGER :: ix,ixo,ir,ireg
      DOUBLE PRECISION :: chireg

      DO ix=1,thisSrc%nxsreg
        IF(ALLOCATED(thisSrc%myXSMesh(ix)%xsmacchi)) THEN
          chireg=thisSrc%myXSMesh(ix)%xsmacchi(ig)
          DO ir=1,thisSrc%myXSMesh(ix)%nreg
            ireg=thisSrc%myXSMesh(ix)%ireg(ir)
            thisSrc%qi1g(ireg)=thisSrc%qi1g(ireg)+psi(ireg)*chireg
          ENDDO
        ENDIF
      ENDDO

    END SUBROUTINE computeMGFS_P0
!===============================================================================
    SUBROUTINE updateInScatter_P0(thisSrc,ig,igstt,igstp)
      CLASS(SourceType_P0),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig
      INTEGER,INTENT(IN) :: igstt
      INTEGER,INTENT(IN) :: igstp
      ! Local Variables
      INTEGER :: ix,ig2,ireg,ir
      DOUBLE PRECISION :: xss_ig2_to_ig

      DO ix=1,thisSrc%nxsreg
        DO ig2=1,thisSrc%ng
          IF(igstt <= ig2 .AND. ig2 <= igstp) THEN
            IF(thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%gmin <= ig2 .AND. &
              ig2 <= thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%gmax .AND. ig /= ig2) THEN
              xss_ig2_to_ig=thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%from(ig2)
              DO ir=1,thisSrc%myXSMesh(ix)%nreg
                ireg=thisSrc%myXSMesh(ix)%ireg(ir)
                thisSrc%qi1g(ireg)=thisSrc%qi1g(ireg)+ &
                  xss_ig2_to_ig*thisSrc%phis(ireg,ig2)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

    END SUBROUTINE updateInScatter_P0 
END MODULE sweeperUtils
