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
      IF(sweepType < 12) THEN
        ALLOCATE(source%qi1g(sweeper%nreg))
        source%qext => source%qi1g
      ELSE
        ALLOCATE(source%qimg(sweeper%nreg,sweeper%ng))
      ENDIF

      ! Allocate sweeper variables
      ALLOCATE(sweeper%phis1g(sweeper%nreg))
      ALLOCATE(sweeper%phis1gd(sweeper%nreg))
      ALLOCATE(sweeper%qbar(sweeper%nreg))
      IF(sweepType < 12) THEN
        ALLOCATE(sweeper%xstr(sweeper%nreg))
      ELSE
        ALLOCATE(sweeper%xstrmg(sweeper%nreg,sweeper%ng))
        ALLOCATE(sweeper%qbarmg(sweeper%nreg,sweeper%ng))
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
      DOUBLE PRECISION :: phio1(0:sweeper%maxsegray),phio2(1:sweeper%maxsegray+1)
      DOUBLE PRECISION :: tphi(sweeper%nreg,1)
      DOUBLE PRECISION :: tau_seg(sweeper%maxsegray,sweeper%ng)
      DOUBLE PRECISION :: &
        exparg(sweeper%maxsegray,SIZE(sweeper%modRayDat%angquad%wtheta),sweeper%ng)
      DOUBLE PRECISION,ALLOCATABLE :: phibar(:,:)
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
          CALL sweeper%mySrc%updateSelfScatter(ig,sweeper%qbar,sweeper%phis1g)
          CALL MOCSolver_Setup1GFSP_PGI(sweeper%myXSMesh,sweeper%nxsreg, &
            sweeper%phis1g,sweeper%nreg,sweeper%xstrmg(:,ig),sweeper%qbarmg(:,ig),ig)
          sweeper%phis1gd = sweeper%phis(:,ig)
        ENDDO !ig

        sweeper%phis = 0.0D0

        CALL CPU_TIME(timeStt)
!-----------------------------------------------------------------------------
        ithd = 1
        npol = SIZE(sweeper%modRayDat%angquad%wtheta)
        wsum = 4.0D0*PI
        ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)
          
        ALLOCATE(phibar(sweeper%nreg,sweeper%ng))
          
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
                DO ig=1,sweeper%ng
                  tau_seg(iseg,ig) = -sweeper%xstrmg(ireg,ig)* &
                    sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%hseg(imseg)
                ENDDO !ig
                irg_seg(iseg) = ireg
              ENDDO !imseg
        
              inextsurf = sweeper%modRayDat%angles(iang)%rays(imray)%nextsurf(1)
              imray = sweeper%modRayDat%angles(iang)%rays(imray)%nextray(1)
              im = sweeper%myModMesh%neigh(inextsurf,im)
            ENDDO !imod
        
            nseglray = iseg
        
            DO ipol=1,npol
              rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
              DO ig=1,sweeper%ng
                DO iseg=1,nseglray
                  exparg(iseg,ipol,ig) = sweeper%expTableDat%EXPT(tau_seg(iseg,ig)*rpol)
                ENDDO !iseg
              ENDDO !ig
            ENDDO !ipol
        
            DO ig=1,sweeper%ng
              sweeper%activeg = ig
              DO ipol=1,npol
                phio1(0) = &
                  sweeper%phiang(ig)%angle(iang)%face(is1)%angflux(ipol,ibc1)
                phio2(nseglray+1) = &
                  sweeper%phiang(ig)%angle(iang)%face(is2)%angflux(ipol,ibc2)
                iseg2 = nseglray + 1
        
                DO iseg1=1,nseglray
                  iseg2 = iseg2 - 1
          
                  ireg1 = irg_seg(iseg1)
                  phid1 = phio1(iseg1-1) - sweeper%qbarmg(ireg1,ig)
                  phid1 = phid1*exparg(iseg1,ipol,ig)
                  !phio1 stores the outgoing angular flux to be used for the next
                  !segment as incoming angular flux.
                  phio1(iseg1) = phio1(iseg1-1) - phid1
                  phibar(ireg1,ig) = phibar(ireg1,ig) + phid1*wtang(ipol)
        
                  ireg2 = irg_seg(iseg2)
                  phid2 = phio2(iseg2+1) - sweeper%qbarmg(ireg2,ig)
                  phid2 = phid2*exparg(iseg2,ipol,ig)
                  !phio1 stores the outgoing angular flux to be used for the next
                  !segment as incoming angular flux.
                  phio2(iseg2) = phio2(iseg2+1) - phid2
                  phibar(ireg2,ig) = phibar(ireg2,ig) + phid2*wtang(ipol)
                ENDDO !iseg
        
                sweeper%phiangmg_out(ig)%angle(iang)%face(is1)%angflux(ipol,ibc1) = &
                  phio2(1)
                sweeper%phiangmg_out(ig)%angle(iang)%face(is2)%angflux(ipol,ibc2) = &
                  phio1(nseglray)
              ENDDO !ipol
            ENDDO !ig
          ENDDO !ilray
        
          sweeper%phis = sweeper%phis + phibar
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
END MODULE openacc
