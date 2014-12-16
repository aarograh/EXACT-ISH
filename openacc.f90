MODULE openacc
!TODO:
!just be sure to change the source everywhere at once, then xsT
!then segement info
!then long ray
  USE sweeperUtils
  USE sweeper
  USE PGIutils

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sweeperType_PGI
  PUBLIC :: sweep2D_prodquad_P0_PGI
  PUBLIC :: MOCSolver_Sweep1G_PGI
  PUBLIC :: MOCSolver_SweepMG_PGI
  PUBLIC :: MOCSolver_SweepACC_PGI

  TYPE :: sweeperType_PGI
    LOGICAL :: hasSource=.FALSE.
    INTEGER :: sweepType=0
    INTEGER :: ninners=0
    INTEGER :: ng=0
    INTEGER :: igstt=0
    INTEGER :: igstp=0
    INTEGER :: activeg=0
    INTEGER :: nsweeps=0
    INTEGER :: nxsreg=0
    INTEGER :: nreg=0
    INTEGER :: imeshstt=0
    INTEGER :: maxsegray=0
    DOUBLE PRECISION :: pz=0.0D0
    DOUBLE PRECISION,ALLOCATABLE :: phis1g(:)
    DOUBLE PRECISION,ALLOCATABLE :: phis1gd(:)
    DOUBLE PRECISION,POINTER :: phis(:,:)
    DOUBLE PRECISION,POINTER :: xstr(:) => NULL()
    DOUBLE PRECISION,POINTER :: xstrmg(:,:) => NULL()
    DOUBLE PRECISION,POINTER :: vol(:) => NULL()
    DOUBLE PRECISION,POINTER :: qbar(:) => NULL()
    DOUBLE PRECISION,POINTER :: qbarmg(:,:) => NULL()
    DOUBLE PRECISION,POINTER :: qext(:) => NULL()
    TYPE(AngFluxBC),POINTER :: phiang1g_in => NULL()
    TYPE(AngFluxBC) :: phiang1g_out
    TYPE(AngFluxBC),ALLOCATABLE :: phiangmg_out(:)
    TYPE(AngFluxBC),POINTER :: phiang(:)
    CLASS(SourceType_PGI),POINTER :: mySrc => NULL()
    TYPE(ModMeshType),POINTER :: myModMesh => NULL()
    TYPE(ModularRayType),POINTER :: modRayDat
    TYPE(CoreLongRayType) :: longRayDat
    TYPE(ModMeshRayPtrArryType),POINTER :: rtmesh(:) => NULL()
    TYPE(XSMeshType),POINTER :: myXSMesh(:)
    TYPE(ExpTableType),ALLOCATABLE :: expTableDat
    TYPE(UpdateBCType_MOC) :: updateBC !maybe, UpdateBC_MOC.f90: define this, %Start() and %Finish() methods.  Might be an MPI thing that I don't need
!    PROCEDURE(absintfc_sweep),POINTER :: sweep => NULL()
!    PROCEDURE(absintfc_setExtSource),POINTER :: setExtSource => NULL()
    PROCEDURE(absintfc_sweep2Dprodquad),POINTER :: sweep2D_prodquad => NULL()
    CONTAINS
      PROCEDURE,PASS :: initialize => initializeSweeper_PGI
      PROCEDURE,PASS :: sweep => MOCSolver_Sweep1G_PGI
      PROCEDURE,PASS :: setExtSource => setExtSource_MOCP0_PGI
! Not needed unless we decide to add a power iteration
!      PROCEDURE,PASS :: calcFissionSrc
  END TYPE sweeperType_PGI

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_sweep(sweeper,ninners,tol,source,psi)
      IMPORT sweeperType_PGI,SourceType_PGI
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
      CLASS(SourceType_PGI),INTENT(INOUT) :: source
      DOUBLE PRECISION,INTENT(INOUT) :: psi(:)
    END SUBROUTINE absintfc_sweep
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_setExtSource(sweeper,source)
      IMPORT sweeperType_PGI,SourceType_PGI
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      CLASS(SourceType_PGI),POINTER,INTENT(IN) :: source
    END SUBROUTINE absintfc_setExtSource
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_sweep2Dprodquad(sweeper,i)
      IMPORT sweeperType_PGI
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: i
    END SUBROUTINE absintfc_sweep2Dprodquad
  END INTERFACE

  CONTAINS
!===============================================================================
    SUBROUTINE initializeSweeper_PGI(sweeper,source,sweepType)
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      CLASS(SourceType_PGI),POINTER,INTENT(INOUT) :: source
      INTEGER,INTENT(IN) :: sweepType
      ! Local Variables
      INTEGER :: ig,iang,iface,i1,i2,i3,i4

      ! Set up source
      ALLOCATE(source)
      sweeper%mySrc => source
      source%nreg = sweeper%nreg
      source%nxsreg = sweeper%nxsreg
      source%ng = sweeper%ng
      source%phis => sweeper%phis
      source%myXSMesh => sweeper%myXSMesh
      ALLOCATE(source%qextmg(source%nreg,source%ng))
      source%qextmg = 0.0D0
      IF(sweepType < 12) THEN
        ALLOCATE(source%qi1g(sweeper%nreg))
        source%qi1g = 0.0D0
        source%qext => source%qi1g
      ELSE
        ALLOCATE(source%qimg(sweeper%nreg,sweeper%ng))
      ENDIF

      ! Allocate sweeper variables
      ALLOCATE(sweeper%phis1g(sweeper%nreg))
      sweeper%phis1g = 0.0D0
      ALLOCATE(sweeper%phis1gd(sweeper%nreg))
      sweeper%phis1gd = 0.0D0
      IF(sweepType < 12) THEN
        ALLOCATE(sweeper%xstr(sweeper%nreg))
        ALLOCATE(sweeper%qbar(sweeper%nreg))
        sweeper%qbar = 0.0D0
      ELSE
        ALLOCATE(sweeper%xstrmg(sweeper%nreg,sweeper%ng))
        sweeper%xstrmg = 0.0D0
        ALLOCATE(sweeper%qbarmg(sweeper%nreg,sweeper%ng))
        sweeper%qbarmg = 0.0D0
        ALLOCATE(sweeper%phiangmg_out(sweeper%ng))
        DO ig=1,sweeper%ng
          i1 = SIZE(sweeper%phiang1g_out%angle)
          ALLOCATE(sweeper%phiangmg_out(ig)%angle(i1))
          DO iang=1,i1
            i2 = SIZE(sweeper%phiang1g_out%angle(iang)%face)
            ALLOCATE(sweeper%phiangmg_out(ig)%angle(iang)%face(i2))
            DO iface=1,i2
              i3 = SIZE(sweeper%phiang1g_out%angle(iang)%face(iface)%angFlux,DIM=1)
              i4 = SIZE(sweeper%phiang1g_out%angle(iang)%face(iface)%angFlux,DIM=2)
              ALLOCATE(sweeper%phiangmg_out(ig)%angle(iang)%face(iface)%angFlux( &
                i3,0:i4-1))
              sweeper%phiangmg_out(ig)%angle(iang)%face(iface)%angFlux = 0.0D0
            ENDDO !iface
          ENDDO !iang
        ENDDO !ig
      ENDIF

    END SUBROUTINE initializeSweeper_PGI
!===============================================================================
    SUBROUTINE setExtSource_MOCP0_PGI(sweeper,source)
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      CLASS(SourceType_PGI),POINTER,INTENT(IN) :: source

      NULLIFY(sweeper%qext)
      sweeper%hasSource=.FALSE.
      sweeper%mySrc => source
      sweeper%qext => source%qext
      sweeper%hasSource=.TRUE.

    END SUBROUTINE setExtSource_MOCP0_PGI
!===============================================================================
    SUBROUTINE MOCSOlver_Setup1GFSP_PGI(thisXSMesh,nxsreg,phis1g,nreg,xstr,qbar,ig)
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

    END SUBROUTINE MOCSolver_Setup1GFSP_PGI
!===============================================================================
    SUBROUTINE MOCSolver_Sweep1G_PGI(sweeper,ninners,tol,source,psi)
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
      CLASS(SourceType_PGI),POINTER,INTENT(INOUT) :: source
      DOUBLE PRECISION,INTENT(INOUT) :: psi(:)
      ! Local Variables
      INTEGER :: i,ig
      DOUBLE PRECISION :: timeStt,timeStp,timeTotal

      timeTotal = 0.0D0

      ! Group loop.  This is actualy in FixedSrcSolver in MPACT
      DO ig=1,sweeper%ng
        ! Set up source
        CALL initExtSource_PGI(source,ig)
        CALL computeMGFS_PGI(source,ig,psi)
        CALL updateInScatter_PGI(source,ig,sweeper%igstt,sweeper%igstp)

        ! This is the real beginning of the sweep routines in MPACT
        IF(1 <= ig .AND. ig <= sweeper%ng) THEN
          sweeper%activeg = ig
          sweeper%phiang1g_in => sweeper%phiang(ig)
          sweeper%phis1g = sweeper%phis(:,ig)

          DO i=1,ninners
            sweeper%nsweeps = sweeper%nsweeps + 1
            !IF(i == ninners) sweeper%sweepCur=.TRUE.
            CALL sweeper%mySrc%updateSelfScatter(ig,sweeper%qbar,sweeper%phis1g)
            CALL MOCSolver_Setup1GFSP_PGI(sweeper%myXSMesh,sweeper%nxsreg, &
              sweeper%phis1g,sweeper%nreg,sweeper%xstr,sweeper%qbar,ig)
            sweeper%phis1gd = sweeper%phis1g

            sweeper%phis1g = 0.0D0

            ! There's an call to sweeper%UpdateBC%finishi here for i > 1

            ! Assumes sweeper%sweepCur == .FALSE.
            CALL CPU_TIME(timeStt)
            CALL sweep2D_prodquad_P0_PGI(sweeper,i)
            CALL CPU_TIME(timeStp)
            timeStp = timeStp - timeStt
            timeTotal = timeTotal + timeStp
            WRITE(*,*) 'Iteration Time: ',timeStp
            WRITE(*,*) 'Accumulated Sweep Time: ',timeTotal
            WRITE(*,*)
          ENDDO !i

          sweeper%phis(:,ig) = sweeper%phis1g

          ! Write to output file for comparison
          IF(ig == 1) WRITE(125,*) SHAPE(sweeper%phis)
          DO i=1,sweeper%nreg
            WRITE(125,*) sweeper%phis(i,ig)
          ENDDO

          ! Update boundary surface flux here, if sweep Cur and associated coarse mesh
          CALL sweeper%UpdateBC%Finish()

          ! hasSource = .FALSE.
        ENDIF
      ENDDO !ig

    END SUBROUTINE MOCSolver_Sweep1G_PGI
!===============================================================================
    SUBROUTINE sweep2D_prodquad_P0_PGI(sweeper,i)
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: i
      ! Local Variables
      INTEGER :: iang,ipol,ilray,imray,imod,im,iside,inextsurf,ifrstreg
      INTEGER :: imseg,iseg,iseg1,iseg2,ibc1,ibc2,nseglray,is1,is2,npol,ireg
      INTEGER :: irg_seg(0:sweeper%maxsegray),ithd
      INTEGER :: ifrstreg_proc,ireg1,ireg2
      DOUBLE PRECISION :: phid1,phid2,wsum,rpol
      DOUBLE PRECISION :: wtangazi,wtang(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION :: phio1(0:sweeper%maxsegray),phio2(1:sweeper%maxsegray+1)
      DOUBLE PRECISION :: tphi(sweeper%nreg,1)
      DOUBLE PRECISION :: tau_seg(sweeper%maxsegray)
      DOUBLE PRECISION :: &
        exparg(sweeper%maxsegray,SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION,ALLOCATABLE :: phibar(:)
      TYPE(LongRayType_Base) :: ilongRay

      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)
      wsum = 4.0D0*PI
      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

      ALLOCATE(phibar(sweeper%nreg))
      tphi(:,ithd) = 0.0D0

      DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
        wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
          sweeper%modRayDat%angquad%walpha(iang)*PI
        DO ipol=1,npol
          wtang(ipol) = wtangazi*sweeper%modRayDat%angquad%wtheta(ipol)* &
            sweeper%modRayDat%angquad%sinpolang(ipol)
        ENDDO !ipol

        phibar = 0.0D0

        DO ilray=1,sweeper%longRayDat%nlongrays(iang)
          ilongRay = sweeper%longRayDat%angles(iang)%longrays(ilray)
          im = ilongRay%ifirstModMesh
          iside = ilongRay%iside(1)
          imray = ilongRay%firstModRay
          ibc1 = ilongRay%BCIndex(1)
          ibc2 = ilongRay%BCIndex(2)
          is1 = ilongRay%iside(1)
          is2 = ilongRay%iside(2)
          iseg = 0

          DO imod=1,ilongRay%nmods
            ifrstreg = sweeper%myModMesh%ifrstfsreg(im)

            DO imseg=1,sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%nseg
              ireg = ifrstreg - ifrstreg_proc + &
                sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%ireg(imseg)
              iseg = iseg + 1
              tau_seg(iseg) = -sweeper%xstr(ireg)* &
                sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%hseg(imseg)
              irg_seg(iseg) = ireg
            ENDDO !imseg

            inextsurf = sweeper%modRayDat%angles(iang)%rays(imray)%nextsurf(1)
            imray = sweeper%modRayDat%angles(iang)%rays(imray)%nextray(1)
            im = sweeper%myModMesh%neigh(inextsurf,im)
          ENDDO !imod

          nseglray = iseg

          DO ipol=1,npol
            rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
            DO iseg=1,nseglray
              exparg(iseg,ipol) = sweeper%expTableDat%EXPT(tau_seg(iseg)*rpol)
            ENDDO !iseg
          ENDDO !ipol

          DO ipol=1,npol
            phio1(0) = &
              sweeper%phiang1g_in%angle(iang)%face(is1)%angflux(ipol,ibc1)
            phio2(nseglray+1) = &
              sweeper%phiang1g_in%angle(iang)%face(is2)%angflux(ipol,ibc2)
            iseg2 = nseglray + 1

            DO iseg1=1,nseglray
              iseg2 = iseg2 - 1

              ireg1 = irg_seg(iseg1)
              phid1 = phio1(iseg1-1) - sweeper%qbar(ireg1)
              phid1 = phid1*exparg(iseg1,ipol)
              !phio1 stores the outgoing angular flux to be used for the next
              !segment as incoming angular flux.
              phio1(iseg1) = phio1(iseg1-1) - phid1
              phibar(ireg1) = phibar(ireg1) + phid1*wtang(ipol)

              ireg2 = irg_seg(iseg2)
              phid2 = phio2(iseg2+1) - sweeper%qbar(ireg2)
              phid2 = phid2*exparg(iseg2,ipol)
              !phio1 stores the outgoing angular flux to be used for the next
              !segment as incoming angular flux.
              phio2(iseg2) = phio2(iseg2+1) - phid2
              phibar(ireg2) = phibar(ireg2) + phid2*wtang(ipol)
            ENDDO !iseg

            sweeper%phiang1g_out%angle(iang)%face(is1)%angflux(ipol,ibc1) = &
              phio2(1)
            sweeper%phiang1g_out%angle(iang)%face(is2)%angflux(ipol,ibc2) = &
              phio1(nseglray)
          ENDDO !ipol
        ENDDO !ilray

        !This is to sum over the angles owned by the current proc.
        !MPACT says it's for polar angles, which I think is not true.
        tphi(:,ithd) = tphi(:,ithd) + phibar
        CALL sweeper%UpdateBC%Start(iang,sweeper%phiang1g_out,sweeper%phiang1g_in)
      ENDDO !iang

      DEALLOCATE(phibar)

      sweeper%phis1g = sweeper%phis1g + tphi(:,ithd)
      sweeper%phis1g = sweeper%phis1g/(sweeper%xstr*sweeper%vol/sweeper%pz) + &
        sweeper%qbar*wsum

    END SUBROUTINE sweep2D_prodquad_P0_PGI
!===============================================================================
    SUBROUTINE MOCSolver_SweepMG_PGI(sweeper,ninners,tol,source,psi)
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
      CLASS(SourceType_PGI),POINTER,INTENT(INOUT) :: source
      DOUBLE PRECISION,INTENT(INOUT) :: psi(:)
      ! Local Variables
      INTEGER :: i,ig
      INTEGER :: iang,ipol,ilray,imray,imod,im,iside,inextsurf,ifrstreg
      INTEGER :: imseg,iseg,iseg1,iseg2,ibc1,ibc2,nseglray,is1,is2,npol,ireg
      INTEGER :: irg_seg(0:sweeper%maxsegray),ithd
      INTEGER :: ifrstreg_proc,ireg1,ireg2
      DOUBLE PRECISION :: timeStt,timeStp,timeTotal
      DOUBLE PRECISION :: phid1,phid2,wsum,rpol
      DOUBLE PRECISION :: wtangazi,wtang(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION :: phio1,phio1d,phio2,phio2d
      DOUBLE PRECISION :: tau_seg(sweeper%maxsegray)
      DOUBLE PRECISION :: exparg(sweeper%maxsegray)
      TYPE(LongRayType_Base) :: ilongRay

      timeTotal = 0.0D0

WRITE(*,*) ASSOCIATED(sweeper%qbarmg)
      DO ig=1,sweeper%ng
        ! Set up external and fission source
        source%qi1g => source%qimg(:,ig)
        CALL initExtSource_PGI(source,ig)
        CALL computeMGFS_PGI(source,ig,psi)
      ENDDO !ig

      DO i=1,ninners
        sweeper%nsweeps = sweeper%nsweeps + 1
        !IF(i == ninners) sweeper%sweepCur=.TRUE.
        DO ig=1,sweeper%ng
          ! Update in-scatter and self-scatter source
          CALL updateInScatter_PGI(source,ig,sweeper%igstt,sweeper%igstp)
          CALL sweeper%mySrc%updateSelfScatter(ig,sweeper%qbarmg(:,ig),sweeper%phis(:,ig))
          CALL MOCSolver_Setup1GFSP_PGI(sweeper%myXSMesh,sweeper%nxsreg, &
            sweeper%phis(:,ig),sweeper%nreg,sweeper%xstrmg(:,ig),sweeper%qbarmg(:,ig),ig)
          sweeper%phis1gd = sweeper%phis(:,ig)
        ENDDO !ig

        sweeper%phis = 0.0D0

        CALL CPU_TIME(timeStt)
!-----------------------------------------------------------------------------
        ithd = 1
        npol = SIZE(sweeper%modRayDat%angquad%wtheta)
        wsum = 4.0D0*PI
        ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)
          
        DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
          wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
            sweeper%modRayDat%angquad%walpha(iang)*PI
          DO ipol=1,npol
            wtang(ipol) = wtangazi*sweeper%modRayDat%angquad%wtheta(ipol)* &
              sweeper%modRayDat%angquad%sinpolang(ipol)
          ENDDO !ipol
        
          DO ilray=1,sweeper%longRayDat%nlongrays(iang)
            ilongRay = sweeper%longRayDat%angles(iang)%longrays(ilray)
            iside = ilongRay%iside(1)
            ibc1 = ilongRay%BCIndex(1)
            ibc2 = ilongRay%BCIndex(2)
            is1 = ilongRay%iside(1)
            is2 = ilongRay%iside(2)

            DO ig=1,sweeper%ng
              sweeper%xstr => sweeper%xstrmg(:,ig)
              sweeper%qbar => sweeper%qbarmg(:,ig)
              sweeper%activeg = ig
              imray = ilongRay%firstModRay
              im = ilongRay%ifirstModMesh
              iseg = 0
        
              DO imod=1,ilongRay%nmods
                ifrstreg = sweeper%myModMesh%ifrstfsreg(im)
        
                DO imseg=1,sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%nseg
                  ireg = ifrstreg - ifrstreg_proc + &
                    sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%ireg(imseg)
                  iseg = iseg + 1
                  tau_seg(iseg) = -sweeper%xstr(ireg)* &
                    sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%hseg(imseg)
                  irg_seg(iseg) = ireg
                ENDDO !imseg
        
                inextsurf = sweeper%modRayDat%angles(iang)%rays(imray)%nextsurf(1)
                imray = sweeper%modRayDat%angles(iang)%rays(imray)%nextray(1)
                im = sweeper%myModMesh%neigh(inextsurf,im)
              ENDDO !imod

              nseglray = iseg
              DO ipol=1,npol
                rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
                DO iseg=1,nseglray
                  exparg(iseg) = sweeper%expTableDat%EXPT(tau_seg(iseg)*rpol)
                ENDDO !iseg

                phio1 = sweeper%phiang(ig)%angle(iang)%face(is1)%angflux(ipol,ibc1)
                phio2 = sweeper%phiang(ig)%angle(iang)%face(is2)%angflux(ipol,ibc2)
                iseg2 = nseglray + 1
        
                DO iseg1=1,nseglray
                  phio1d = phio1
                  phio2d = phio2
                  iseg2 = iseg2 - 1
          
                  ireg1 = irg_seg(iseg1)
                  phid1 = (phio1d - sweeper%qbar(ireg1)) * exparg(iseg1)
                  !phio1 stores the outgoing angular flux to be used for the next
                  !segment as incoming angular flux.
                  phio1 = phio1d - phid1
                  sweeper%phis(ireg1,ig) = sweeper%phis(ireg1,ig) + phid1*wtang(ipol)
        
                  ireg2 = irg_seg(iseg2)
                  phid2 = (phio2d - sweeper%qbar(ireg2)) * exparg(iseg2)
                  !phio1 stores the outgoing angular flux to be used for the next
                  !segment as incoming angular flux.
                  phio2 = phio2d - phid2
                  sweeper%phis(ireg2,ig) = sweeper%phis(ireg2,ig) + phid2*wtang(ipol)
                ENDDO !iseg
        
                sweeper%phiangmg_out(ig)%angle(iang)%face(is1)%angflux(ipol,ibc1) = phio2
                sweeper%phiangmg_out(ig)%angle(iang)%face(is2)%angflux(ipol,ibc2) = phio1
              ENDDO !ipol
            ENDDO !ig
          ENDDO !ilray
        
          DO ig=1,sweeper%ng
            CALL sweeper%UpdateBC%Start(iang,sweeper%phiangmg_out(ig),sweeper%phiang(ig))
          ENDDO
        ENDDO !iang
        
        DO ig=1,sweeper%ng
          sweeper%phis(:,ig) = sweeper%phis(:,ig)/(sweeper%xstrmg(:,ig)* &
            sweeper%vol/sweeper%pz) + sweeper%qbarmg(:,ig)*wsum
        ENDDO !ig

        CALL CPU_TIME(timeStp)
!-------------------------------------------------------------------------------
        timeStp = timeStp - timeStt
        timeTotal = timeTotal + timeStp
        WRITE(*,*) 'Iteration Time: ',timeStp
        WRITE(*,*) 'Accumulated Sweep Time: ',timeTotal
        WRITE(*,*)
      ENDDO !i

      ! Write to output file for comparison
      WRITE(125,*) SHAPE(sweeper%phis)
      DO ig=1,sweeper%ng
        DO i=1,sweeper%nreg
          WRITE(125,*) sweeper%phis(i,ig)
        ENDDO
      ENDDO

      ! Update boundary surface flux here, if sweep Cur and associated coarse mesh
      CALL sweeper%UpdateBC%Finish()

    END SUBROUTINE MOCSolver_SweepMG_PGI
!===============================================================================
    SUBROUTINE MOCSolver_SweepACC_PGI(sweeper,ninners,tol,source,psi)
      CLASS(sweeperType_PGI),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
      CLASS(SourceType_PGI),POINTER,INTENT(INOUT) :: source
      DOUBLE PRECISION,INTENT(INOUT) :: psi(:)
      ! Local Variables
      LOGICAL :: UpdateBC
      INTEGER :: i,ig,iangstt,iangstp,maxneigh,maxsegray,nang,nsegs,imodray
      INTEGER :: nrays,ng,i1,irefl,iface,nfaces
      INTEGER :: iang,ipol,ilray,imray,imod,im,iside,inextsurf
      INTEGER :: imseg,iseg,iseg1,iseg2,ibc1,ibc2,nseglray,is1,is2,npol,ireg
      INTEGER :: irg_seg(0:sweeper%maxsegray),ithd,ifrstreg(SIZE(sweeper%myModMesh%ifrstfsreg))
      INTEGER :: ifrstreg_proc,ireg1,ireg2,nlrays,iray,nfacerays
      INTEGER :: nlongrays(sweeper%modRayDat%iangstt:sweeper%modRayDat%iangstp)
      INTEGER :: iang2irefl(sweeper%UpdateBC%nfaces,sweeper%UpdateBC%nangles)
      INTEGER,ALLOCATABLE :: lrayiside(:,:),BCIndex(:,:),firstModRay(:),ifirstModMesh(:)
      INTEGER,ALLOCATABLE :: ifirstfsreg(:),nseg(:),rtmeshireg(:),nextsurf(:)
      INTEGER,ALLOCATABLE :: neigh(:,:),nmods(:),nmodray(:),nmodrayseg(:),nextray(:)
      DOUBLE PRECISION :: timeStt,timeStp,timeTotal
      DOUBLE PRECISION :: phid1,phid2,wsum,rpol
      DOUBLE PRECISION :: wtangazi,wtang(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION :: phio1,phio1d,phio2,phio2d
      DOUBLE PRECISION :: tau_seg(sweeper%maxsegray)
      DOUBLE PRECISION :: exparg(sweeper%maxsegray)
      DOUBLE PRECISION :: dlr(sweeper%modRayDat%iangstt:sweeper%modRayDat%iangstp)
      DOUBLE PRECISION :: walpha(sweeper%modRayDat%iangstt:sweeper%modRayDat%iangstp)
      DOUBLE PRECISION :: wtheta(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION :: sinpolang(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION :: rsinpolang(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION :: xstrmg(sweeper%nreg,sweeper%ng),qbarmg(sweeper%nreg,sweeper%ng)
      DOUBLE PRECISION :: phis(sweeper%nreg,sweeper%ng)
      DOUBLE PRECISION,ALLOCATABLE :: angflux(:,:,:,:,:),hseg(:)
      TYPE(LongRayType_Base) :: ilongRay

      timeTotal = 0.0D0

      DO ig=1,sweeper%ng
        ! Set up external and fission source
        source%qi1g => source%qimg(:,ig)
        CALL initExtSource_PGI(source,ig)
        CALL computeMGFS_PGI(source,ig,psi)
      ENDDO !ig

      DO i=1,ninners
        sweeper%nsweeps = sweeper%nsweeps + 1
        !IF(i == ninners) sweeper%sweepCur=.TRUE.
        DO ig=1,sweeper%ng
          ! Update in-scatter and self-scatter source
          CALL updateInScatter_PGI(source,ig,sweeper%igstt,sweeper%igstp)
          CALL sweeper%mySrc%updateSelfScatter(ig,sweeper%qbarmg(:,ig),sweeper%phis(:,ig))
          CALL MOCSolver_Setup1GFSP_PGI(sweeper%myXSMesh,sweeper%nxsreg, &
            sweeper%phis(:,ig),sweeper%nreg,sweeper%xstrmg(:,ig),sweeper%qbarmg(:,ig),ig)
          sweeper%phis1gd = sweeper%phis(:,ig)
        ENDDO !ig

        sweeper%phis = 0.0D0

        ithd = 1
        wsum = 4.0D0*PI
        ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

        iangstt=sweeper%modRayDat%iangstt
        iangstp=sweeper%modRayDat%iangstp
        npol=SIZE(sweeper%modRayDat%angquad%wtheta)

        ! Fill pre-sized arrays
        phis = 0.0D0
        xstrmg = sweeper%xstrmg
        qbarmg = sweeper%qbarmg
        dlr = sweeper%modRayDat%angles(:)%dlr
        walpha = sweeper%modRayDat%angquad%walpha
        wtheta = sweeper%modRayDat%angquad%wtheta
        sinpolang = sweeper%modRayDat%angquad%sinpolang
        rsinpolang = sweeper%modRayDat%angquad%rsinpolang
        iang2irefl = sweeper%updateBC%iang2irefl
        ifrstreg = sweeper%myModMesh%ifrstfsreg
        ! Set some counters
        ng = sweeper%ng
        nlongrays = sweeper%longRayDat%nlongrays
        nlrays = SUM(nlongrays)
        nfaces = sweeper%updateBC%nfaces
        ALLOCATE(lrayiside(2,nlrays)); ALLOCATE(BCIndex(2,nlrays))
        ALLOCATE(firstModRay(nlrays))
        ALLOCATE(ifirstModMesh(nlrays))
        ALLOCATE(nmods(nlrays))
        ALLOCATE(nmodray(0:nlrays)); nmodray = 0
        maxsegray = 0
        nfacerays = 0
        DO iang=iangstt,iangstp
          DO ilray=1,nlongrays(iang)
            maxsegray = maxsegray + sweeper%longRayDat%angles(iang)%longrays(ilray)%nmods
          ENDDO !ilray
          DO ig=1,ng
            DO iface=1,nfaces
              nfacerays = MAX(nfacerays,SIZE(sweeper%phiang(ig)%angle(iang)%face(iface)%angflux,DIM=2))
            ENDDO !iface
          ENDDO !ig
        ENDDO !iang
        ALLOCATE(nmodrayseg(maxsegray))
        ALLOCATE(nextsurf(maxsegray))
        ALLOCATE(nextray(maxsegray))
        ALLOCATE(neigh(LBOUND(sweeper%myModMesh%neigh,DIM=1):UBOUND(sweeper%myModMesh%neigh,DIM=1), &
          LBOUND(sweeper%myModMesh%neigh,DIM=2):UBOUND(sweeper%myModMesh%neigh,DIM=2)))
        neigh = sweeper%myModMesh%neigh
        maxneigh = MAXVAL(neigh)
        maxsegray = sweeper%maxsegray
        nang = iangstp - iangstt + 1
        ALLOCATE(rtmeshireg(nlrays*maxsegray))
        ALLOCATE(hseg(nlrays*maxsegray))
        ALLOCATE(angflux(npol,0:nfacerays-1,nfaces,ng,nang)); angflux = 0.0D0
        
        iray = 0
        nsegs = 0
        imodray = 0
        DO iang=iangstt,iangstp
          DO ilray=1,nlongrays(iang)
            iray = iray + 1
            ilongRay = sweeper%longRayDat%angles(iang)%longrays(ilray)
            lrayiside(:,iray) = ilongRay%iside(1:2)
            BCIndex(:,iray) = ilongRay%BCIndex(1:2)
            firstModRay(iray) = iLongRay%firstModRay
            ifirstModMesh(iray) = ilongRay%ifirstModMesh
            nmods(iray) = iLongRay%nmods
            iseg = 0
            imray = ilongRay%firstModRay
            im = ilongRay%ifirstModMesh
            DO imod = 1,ilongRay%nmods
              imodray = imodray + 1
              nmodrayseg(imodray) = sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%nseg
              DO imseg=1,sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%nseg
                iseg = iseg + 1
                ireg = ifrstreg(im) - ifrstreg_proc + &
                  sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%ireg(imseg)
                rtmeshireg((iray-1)*maxsegray+iseg) = ireg
                hseg((iray-1)*maxsegray+iseg) = sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%hseg(imseg)
              ENDDO
              nextsurf(imodray) = sweeper%modRayDat%angles(iang)%rays(imray)%nextsurf(1)
              nextray(imodray) = sweeper%modRayDat%angles(iang)%rays(imray)%nextray(1)
              imray = nextray(imodray)
              im = neigh(nextsurf(imodray),im)
            ENDDO !imod
            nsegs = nsegs + maxsegray
            nmodray(iray) = imodray
          ENDDO !ilray
          DO ig=1,ng
            DO iface=1,nfaces
              ipol = SIZE(sweeper%phiang(ig)%angle(iang)%face(iface)%angflux,DIM=1)
              i1 = SIZE(sweeper%phiang(ig)%angle(iang)%face(iface)%angflux,DIM=2)
              angflux(1:ipol,0:i1-1,iface,ig,iang) = sweeper%phiang(ig)%angle(iang)%face(iface)%angflux
            ENDDO !iface
          ENDDO !ig
        ENDDO !iang

        WRITE(*,FMT='(a,i0,a,i0,a)') 'Solving ',nlrays,' rays and ',sweeper%nreg,' regions...'
        CALL CPU_TIME(timeStt)
!-----------------------------------------------------------------------------
!$acc parallel private(phio1,phio2,phio1d,phio2d,ig,iang,imod,imseg,ireg,wtang) 
!$acc loop
        DO iray = 1,nlrays
          nrays = 0
          DO iang=iangstt,iangstp
            nrays = nrays + nlongrays(iang)
            IF(nrays >= iray) THEN
              IF(nrays == iray) updateBC = .TRUE.
              EXIT
            ENDIF
          ENDDO !iang
          wtangazi = dlr(iang)*walpha(iang)*PI
          wtang = wtangazi*wtheta*sinpolang
!$acc loop private(tau_seg,irg_seg)
          DO ig=1,ng
            i1 = (iray-1)*maxsegray 
            imray = firstModRay(iray)
            im = ifirstModMesh(iray)
            iseg = 0
            imodray = nmodray(iray-1)
            DO imod=1,nmods(iray)
              imodray = imodray + 1
              DO imseg = 1,nmodrayseg(imodray)
                iseg = iseg + 1
                i1 = i1 + 1
                ireg = rtmeshireg(i1)
                tau_seg(iseg) = -xstrmg(ireg,ig)*hseg(i1)
                irg_seg(iseg) = ireg
              ENDDO !imseg
              imray = nextray(imodray) 
              im = neigh(nextsurf(imodray),im)
            ENDDO !imod

            ibc1 = BCIndex(1,iray)
            ibc2 = BCIndex(2,iray)
            is1 = lrayiside(1,iray)
            is2 = lrayiside(2,iray)

            nseglray = iseg
!$acc loop private(exparg)
            DO ipol=1,npol
              DO iseg = 1,nseglray
!                exparg(iseg) = sweeper%expTableDat%EXPT(tau_seg(iseg)*rsinpolang(ipol))
                exparg(iseg) = 1.0D0 - EXP(tau_seg(iseg)*rsinpolang(ipol))
              ENDDO !iseg

              phio1 = angflux(ipol,ibc1,is1,ig,iang)
              phio2 = angflux(ipol,ibc2,is2,ig,iang)
              iseg2 = nseglray + 1
        
              DO iseg1=1,nseglray
                phio1d = phio1
                phio2d = phio2
                iseg2 = iseg2 - 1
          
                ireg1 = irg_seg(iseg1)
                phid1 = (phio1d - qbarmg(ireg1,ig)) * exparg(iseg1)
                !phio1 stores the outgoing angular flux to be used for the next
                !segment as incoming angular flux.
                phio1 = phio1d - phid1
                phis(ireg1,ig) = phis(ireg1,ig) + phid1*wtang(ipol)
        
                ireg2 = irg_seg(iseg2)
                phid2 = (phio2d - qbarmg(ireg2,ig)) * exparg(iseg2)
                !phio1 stores the outgoing angular flux to be used for the next
                !segment as incoming angular flux.
                phio2 = phio2d - phid2
                phis(ireg2,ig) = phis(ireg2,ig) + phid2*wtang(ipol)
              ENDDO !iseg
        
               angflux(ipol,ibc1,is1,ig,iang) = phio2
               angflux(ipol,ibc2,is2,ig,iang) = phio1
            ENDDO !ipol
            IF(updateBC) THEN
              DO iface=1,nfaces
                irefl = iang2irefl(iface,iang)
                 angflux(:,:,iface,ig,irefl) = angflux(:,:,iface,ig,iang)
              ENDDO !iface
            ENDIF
          ENDDO !ig
        ENDDO !iray 
!$acc end parallel
              
        sweeper%phis = phis
        
        DO ig=1,sweeper%ng
          sweeper%phis(:,ig) = sweeper%phis(:,ig)/(sweeper%xstrmg(:,ig)* &
            sweeper%vol/sweeper%pz) + sweeper%qbarmg(:,ig)*wsum
        ENDDO !ig

!-------------------------------------------------------------------------------
        CALL CPU_TIME(timeStp)
        timeStp = timeStp - timeStt
        timeTotal = timeTotal + timeStp
        WRITE(*,*) 'Iteration Time: ',timeStp
        WRITE(*,*) 'Accumulated Sweep Time: ',timeTotal
        WRITE(*,*)
      ENDDO !i

      ! Write to output file for comparison
      WRITE(125,*) SHAPE(sweeper%phis)
      DO ig=1,sweeper%ng
        DO i=1,sweeper%nreg
          WRITE(125,*) sweeper%phis(i,ig)
        ENDDO
      ENDDO

      ! Update boundary surface flux here, if sweep Cur and associated coarse mesh
      CALL sweeper%UpdateBC%Finish()

    END SUBROUTINE MOCSolver_SweepACC_PGI
END MODULE openacc
