MODULE sweeperUtils

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: PI
  PUBLIC :: AngFluxBC
  PUBLIC :: SourceType
  PUBLIC :: SourceType_P0
  PUBLIC :: ModularRayType
  PUBLIC :: CoreLongRayType
  PUBLIC :: ModMeshRayPtrArryType
  PUBLIC :: XSMeshType
  PUBLIC :: ExpTableType
  PUBLIC :: UpdateBCType_MOC
  PUBLIC :: LongRayType_Base
  PUBLIC :: ModMeshType

  DOUBLE PRECISION :: PI=3.141592653589793D0
  INTEGER,PARAMETER :: SINGLE_LEVEL_EXP_TABLE=1
  INTEGER,PARAMETER :: TWO_LEVEL_EXP_TABLE=2
  INTEGER,PARAMETER :: LINEAR_EXP_TABLE=3

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
    DOUBLE PRECISION,POINTER :: phisd(:,:) => NULL()
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
    DOUBLE PRECISION,POINTER :: qimg(:,:) => NULL()
    CONTAINS
      PROCEDURE,PASS :: updateSelfScatter => updateSelfScatter_P0
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
    INTEGER :: npol
    INTEGER :: nazi
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
    INTEGER :: nmods=0
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
    INTEGER :: nseg
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
    INTEGER :: tableType=LINEAR_EXP_TABLE
    DOUBLE PRECISION :: minVal=1.0D0
    DOUBLE PRECISION :: maxVal=0.0D0
    DOUBLE PRECISION :: rdx=0.0D0
    DOUBLE PRECISION,ALLOCATABLE :: table2D(:,:)
    DOUBLE PRECISION,ALLOCATABLE :: table3D(:,:,:)
    CONTAINS
      PROCEDURE,PASS :: EXPT
      PROCEDURE,PASS :: EXPT_vectoripol
  END TYPE ExpTableType

  TYPE :: UpdateBCType_MOC
    INTEGER :: offset=-1
    INTEGER :: nfaces=-1
    INTEGER :: bcType(6)=-1
    INTEGER :: nangles=-1
    INTEGER :: iangstt=-1
    INTEGER :: iangstp=-1
    INTEGER,ALLOCATABLE :: iang2irefl(:,:)
    CONTAINS
      PROCEDURE,PASS :: start => UpdateBC_Start
      PROCEDURE,PASS :: finish => UpdateBC_Finish
  END TYPE UpdateBCType_MOC

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
    SUBROUTINE updateSelfScatter_P0(thisSrc,ig,qbar,phis1g)
      CLASS(SourceType_P0),INTENT(IN) :: thisSrc
      INTEGER,INTENT(IN) :: ig
      DOUBLE PRECISION,INTENT(INOUT) :: qbar(:)
      DOUBLE PRECISION,INTENT(IN) :: phis1g(:)
      ! Local variables
      DOUBLE PRECISION,PARAMETER :: r4pi=0.25D0/3.141592653589793D0
      INTEGER :: ix,ir,ireg
      DOUBLE PRECISION :: xstrg,xssgg,rxstrg4pi

      qbar = thisSrc%qi1g
      ! Assumes no XS splitting (See SourceTypes.f90:470
      DO ix=1,thisSrc%nxsreg
!WRITE(*,*) ix,ig,':',SIZE(thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%from)
        !xssgg = thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%from(ig)
        xstrg = thisSrc%myXSMesh(ix)%xsmactr(ig)
        rxstrg4pi = r4pi/xstrg
        DO ir=1,thisSrc%myXSMesh(ix)%nreg
          ireg = thisSrc%myXSMesh(ix)%ireg(ir)
          qbar(ireg) = (qbar(ireg))*rxstrg4pi
        ENDDO !ir
      ENDDO !ix

    END SUBROUTINE updateSelfScatter_P0
!===============================================================================
    SUBROUTINE initExtSource_P0(thisSrc,ig)
      CLASS(SourceType_P0),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig

      thisSrc%qi1g = thisSrc%qextmg(:,ig)

    END SUBROUTINE initExtSource_P0
!===============================================================================
    SUBROUTINE computeMGFS_P0(thisSrc,ig,psi)
      CLASS(SourceType_P0),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig
      DOUBLE PRECISION,INTENT(IN) :: psi(:)
      ! Local variables
      INTEGER :: ix,ir,ireg
      DOUBLE PRECISION :: chireg

      DO ix=1,thisSrc%nxsreg
        IF(ALLOCATED(thisSrc%myXSMesh(ix)%xsmacchi)) THEN
          chireg = thisSrc%myXSMesh(ix)%xsmacchi(ig)
          DO ir=1,thisSrc%myXSMesh(ix)%nreg
            ireg = thisSrc%myXSMesh(ix)%ireg(ir)
            thisSrc%qi1g(ireg) = thisSrc%qi1g(ireg) + psi(ireg)*chireg
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
              ig2 <= thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%gmax) THEN
              xss_ig2_to_ig = thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%from(ig2)
              DO ir=1,thisSrc%myXSMesh(ix)%nreg
                ireg = thisSrc%myXSMesh(ix)%ireg(ir)
                thisSrc%qi1g(ireg) = thisSrc%qi1g(ireg) +  &
                  xss_ig2_to_ig*(thisSrc%phis(ireg,ig2)-thisSrc%phisd(ireg,ig2))
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

    END SUBROUTINE updateInScatter_P0
!===============================================================================
    FUNCTION EXPT(ET,x) RESULT(ans)
      CLASS(ExpTableType),INTENT(IN) :: ET
      DOUBLE PRECISION,INTENT(IN) :: x
      DOUBLE PRECISION :: ans
      ! Local Variables

      IF(ET%minVal <= x .AND. x <= ET%maxVal) THEN
        SELECTCASE(ET%tableType)
          CASE(SINGLE_LEVEL_EXP_TABLE)
          CASE(TWO_LEVEL_EXP_TABLE)
          CASE(LINEAR_EXP_TABLE)
            ans = EXPT_Linear(ET,x)
          CASE DEFAULT
            ans=1.0D0 - EXP(X)
        END SELECT
      ELSE
        IF(x < -700.0D0) THEN
          ans = 1.0D0
        ELSE
          ans = 1.0D0 - EXP(x)
        ENDIF
      ENDIF

    END FUNCTION EXPT
!===============================================================================
    FUNCTION EXPT_vectoripol(ET,x) RESULT(ans)
      CLASS(ExpTableType),INTENT(IN) :: ET
      DOUBLE PRECISION,INTENT(IN) :: x(:)
      DOUBLE PRECISION :: ans(SIZE(x))
      ! Local Variables
      INTEGER :: i,j
      DO i=1,SIZE(x)
        IF(ET%minVal <= x(i) .AND. x(i) <= ET%maxVal) THEN
          SELECTCASE(ET%tableType)
            CASE(SINGLE_LEVEL_EXP_TABLE)
            CASE(TWO_LEVEL_EXP_TABLE)
            CASE(LINEAR_EXP_TABLE)
              ans(i)=EXPT_Linear(ET,x(i))
!               j = FLOOR(x(i)*ET%rdx)
!               ans(i) = ET%table2D(1,j)*x(i) + ET%table2D(2,j)
            CASE DEFAULT
              ans(i)=1.0D0 - EXP(x(i))
          END SELECT
        ELSE
          IF(x(i) < -700.0D0) THEN
            ans(i) = 1.0D0
          ELSE
            ans(i) = 1.0D0 - EXP(x(i))
          ENDIF
        ENDIF
      ENDDO

    END FUNCTION EXPT_vectoripol
!===============================================================================
    FUNCTION EXPT_Linear(ET,x) RESULT(ans)
      CLASS(ExpTableType),INTENT(IN) :: ET
      DOUBLE PRECISION,INTENT(IN) :: x
      DOUBLE PRECISION :: ans
      INTEGER :: i

      i = FLOOR(x*ET%rdx)
      ans = ET%table2D(1,i)*x + ET%table2D(2,i)

    END FUNCTION EXPT_Linear
!===============================================================================
    SUBROUTINE UpdateBC_Start(thisBCU,iang,outgoing,incoming)
      CLASS(UpdateBCType_MOC),INTENT(IN) :: thisBCU
      INTEGER,INTENT(IN) :: iang
      TYPE(AngFluxBC),INTENT(INOUT) :: outgoing
      TYPE(AngFluxBC),INTENT(INOUT) :: incoming
      ! Local Variables
      INTEGER,PARAMETER :: VACUUMBC=0
      INTEGER,PARAMETER :: REFLECTIVEBC=1
      INTEGER,PARAMETER :: PERIODICBC=2
      INTEGER,PARAMETER :: PARALLELBC=3
      INTEGER,PARAMETER :: REENTRANTBC=4
      INTEGER,PARAMETER :: ANTISYMMBC=5
      INTEGER :: iface,irefl

      DO iface=1,thisBCU%nfaces
        irefl = thisBCU%iang2irefl(iface,iang)

        SELECTCASE(thisBCU%bcType(iface))
          CASE(VACUUMBC)
            CONTINUE
          CASE(REFLECTIVEBC)
            incoming%angle(irefl)%face(iface)%angflux = &
              outgoing%angle(iang)%face(iface)%angflux
          CASE(PERIODICBC)
            STOP 666
          CASE(PARALLELBC)
            incoming%angle(irefl)%face(iface)%angflux = &
              outgoing%angle(iang)%face(iface)%angflux
          CASE(ANTISYMMBC)
            STOP 667
        END SELECT

      ENDDO !iface

   END SUBROUTINE UpdateBC_Start
!===============================================================================
    SUBROUTINE UpdateBC_Finish(thisBCU)
      CLASS(UpdateBCType_MOC),INTENT(IN) :: thisBCU

    END SUBROUTINE UpdateBC_Finish

END MODULE sweeperUtils
