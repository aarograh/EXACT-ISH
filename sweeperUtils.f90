MODULE sweeperUtils

  PUBLIC :: AngFluxBC
  PUBLIC :: SourceType
  PUBLIC :: SourceType_P0
  PUBLIC :: ModularRayType
  PUBLIC :: CoreLongRayType
  PUBLIC :: ModMeshRayPtrArryType

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
!    TYPE(XSMeshType),POINTER :: myXSMesh(:) => NULL()
  END TYPE SourceType

  TYPE,EXTENDS(SourceType) :: SourceType_P0
    DOUBLE PRECISION,POINTER :: qi1g(:)
    CONTAINS
      PROCEDURE,PASS :: updateSelfScatter_P0
  END TYPE SourceType_P0

  TYPE :: ModMeshType
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

  CONTAINS
!===============================================================================
!TODO: MOCSweeper_P0.f90:47
    SUBROUTINE updateSelfScatter_P0(thisSrc,igroup,qbar,phis1g)
      CLASS(SourceType_P0),INTENT(IN) :: thisSrc
      INTEGER,INTENT(IN) :: igroup
      DOUBLE PRECISION,INTENT(INOUT) :: qbar(:)
      DOUBLE PRECISION,INTENT(IN) :: phis1g
    END SUBROUTINE updateSelfScatter_P0

END MODULE sweeperUtils
