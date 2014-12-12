MODULE PGIutils

  USE sweeperUtils

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SourceType_PGI
  PUBLIC :: initExtSource_PGI
  PUBLIC :: computeMGFS_PGI
  PUBLIC :: updateInScatter_PGI

  TYPE :: SourceType_PGI
    INTEGER :: nreg=0
    INTEGER :: nxsreg=0
    INTEGER :: ng=0
    DOUBLE PRECISION,POINTER :: phis(:,:) => NULL()
    DOUBLE PRECISION,POINTER :: qi1g(:) => NULL()
    DOUBLE PRECISION,POINTER :: qext(:) => NULL()
    DOUBLE PRECISION,POINTER :: split(:) => NULL()
    DOUBLE PRECISION,POINTER :: qextmg(:,:) => NULL()
    TYPE(XSMeshType),POINTER :: myXSMesh(:) => NULL()
    CONTAINS
      PROCEDURE,PASS :: updateSelfScatter => updateSelfScatter_PGI
  END TYPE SourceType_PGI

  CONTAINS
!===============================================================================
    SUBROUTINE updateSelfScatter_PGI(thisSrc,ig,qbar,phis1g)
      CLASS(SourceType_PGI),INTENT(IN) :: thisSrc
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
        xssgg = thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%from(ig)
        xstrg = thisSrc%myXSMesh(ix)%xsmactr(ig)
        rxstrg4pi = r4pi/xstrg
        DO ir=1,thisSrc%myXSMesh(ix)%nreg
          ireg = thisSrc%myXSMesh(ix)%ireg(ir)
          qbar(ireg) = (qbar(ireg) + xssgg*phis1g(ireg))*rxstrg4pi
        ENDDO !ir
      ENDDO !ix

    END SUBROUTINE updateSelfScatter_PGI
!===============================================================================
    SUBROUTINE initExtSource_PGI(thisSrc,ig)
      CLASS(SourceType_PGI),INTENT(INOUT) :: thisSrc
      INTEGER,INTENT(IN) :: ig

      thisSrc%qi1g = thisSrc%qextmg(:,ig)

    END SUBROUTINE initExtSource_PGI
!===============================================================================
    SUBROUTINE computeMGFS_PGI(thisSrc,ig,psi)
      CLASS(SourceType_PGI),INTENT(INOUT) :: thisSrc
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

    END SUBROUTINE computeMGFS_PGI
!===============================================================================
    SUBROUTINE updateInScatter_PGI(thisSrc,ig,igstt,igstp)
      CLASS(SourceType_PGI),INTENT(INOUT) :: thisSrc
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
              xss_ig2_to_ig = thisSrc%myXSMesh(ix)%xsmacsc(ig,0)%from(ig2)
              DO ir=1,thisSrc%myXSMesh(ix)%nreg
                ireg = thisSrc%myXSMesh(ix)%ireg(ir)
                thisSrc%qi1g(ireg) = thisSrc%qi1g(ireg) +  &
                  xss_ig2_to_ig*thisSrc%phis(ireg,ig2)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

    END SUBROUTINE updateInScatter_PGI
END MODULE PGIutils
