MODULE sweeper
!TODO: Add timers for these loops

  USE sweeperUtils

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sweeperType
  PUBLIC :: sweep2D_prodquad_P0
  PUBLIC :: expoa,expob
  PUBLIC :: expoab

  !For polar angle dependent exponential table inlining.
  DOUBLE PRECISION,ALLOCATABLE,SAVE :: expoa(:,:),expob(:,:)
  DOUBLE PRECISION,ALLOCATABLE,SAVE :: expoab(:)

  TYPE :: sweeperType
    LOGICAL :: hasSource=.FALSE.
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
    !for use in group inner
    DOUBLE PRECISION,POINTER :: xstrmg(:,:) => NULL()
    DOUBLE PRECISION,POINTER :: vol(:) => NULL()
    DOUBLE PRECISION,POINTER :: qbar(:) => NULL()
    !for use in group inenr
    DOUBLE PRECISION,POINTER :: qbarmg(:,:) => NULL()
    DOUBLE PRECISION,POINTER :: qext(:) => NULL()
    TYPE(AngFluxBC),POINTER :: phiang1g_in => NULL()
    TYPE(AngFluxBC) :: phiang1g_out
    TYPE(AngFluxBC),ALLOCATABLE :: phiangmg_out(:)
    TYPE(AngFluxBC),POINTER :: phiang(:)
    CLASS(SourceType_P0),POINTER :: mySrc => NULL()
    TYPE(ModMeshType),POINTER :: myModMesh => NULL()
    TYPE(ModularRayType),POINTER :: modRayDat
    TYPE(CoreLongRayType) :: longRayDat
    TYPE(ModMeshRayPtrArryType),POINTER :: rtmesh(:) => NULL()
    TYPE(XSMeshType),POINTER :: myXSMesh(:)
    TYPE(ExpTableType),ALLOCATABLE :: expTableDat
    TYPE(UpdateBCType_MOC) :: updateBC !maybe, UpdateBC_MOC.f90: define this, %Start() and %Finish() methods.  Might be an MPI thing that I don't need
    PROCEDURE(absintfc_sweep),POINTER :: sweep => NULL()
    PROCEDURE(absintfc_setExtSource),POINTER :: setExtSource => NULL()
    PROCEDURE(absintfc_sweep2Dprodquad),POINTER :: sweep2D_prodquad => NULL()
    CONTAINS
      PROCEDURE,PASS :: initialize => initializeSweeper
! Not needed unless we decide to add a power iteration
!      PROCEDURE,PASS :: calcFissionSrc
  END TYPE sweeperType

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_sweep(sweeper,ninners,tol,source,psi)
      IMPORT sweeperType,SourceType
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
      CLASS(SourceType),INTENT(INOUT) :: source
      DOUBLE PRECISION,INTENT(INOUT) :: psi(:)
    END SUBROUTINE absintfc_sweep
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_setExtSource(sweeper,source)
      IMPORT sweeperType,SourceType
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      CLASS(SourceType),POINTER,INTENT(IN) :: source
    END SUBROUTINE absintfc_setExtSource
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE absintfc_sweep2Dprodquad(sweeper,i)
      IMPORT sweeperType
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: i
    END SUBROUTINE absintfc_sweep2Dprodquad
  END INTERFACE

  CONTAINS
!===============================================================================
    SUBROUTINE initializeSweeper(sweeper,source)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      CLASS(sourceType),POINTER,INTENT(INOUT) :: source

      ! Local Variables
      INTEGER :: ig,iang,iface,i1,i2,i3,i4,i5

      ! Set up source
      ALLOCATE(SourceType_P0 :: source)
      SELECTTYPE(source); TYPE IS(SourceType_P0)
        sweeper%mySrc => source
        ALLOCATE(source%qi1g(sweeper%nreg))
        ALLOCATE(source%qimg(sweeper%nreg,sweeper%ng))
        source%qi1g = 0.0D0
        source%qext => source%qi1g
      END SELECT
      source%nreg = sweeper%nreg
      source%nxsreg = sweeper%nxsreg
      source%ng = sweeper%ng
      source%phis => sweeper%phis
      ALLOCATE(source%phisd(source%nreg,source%ng))
      !temporary set to zero for use in group inner
      source%phisd = 0.0D0
      source%myXSMesh => sweeper%myXSMesh
      ALLOCATE(source%qextmg(source%nreg,source%ng))
      source%qextmg = 0.0D0

      ! Allocate sweeper variables
      ALLOCATE(sweeper%phis1g(sweeper%nreg))
      sweeper%phis1g = 0.0D0
      ALLOCATE(sweeper%phis1gd(sweeper%nreg))
      sweeper%phis1gd = 0.0D0
      ALLOCATE(sweeper%qbar(sweeper%nreg))
      sweeper%qbar = 0.0D0
      ALLOCATE(sweeper%xstr(sweeper%nreg))
      sweeper%xstr = 0.0D0

      !for use in group inner
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

      ! Set up expoa and expob and expoab in case of table1D for exponential inlining
      IF(ALLOCATED(sweeper%expTableDat%table3D) .AND. .NOT.ALLOCATED(expoa) .AND. &
         .NOT.ALLOCATED(expob) .AND. .NOT.ALLOCATED(expoab)) THEN
        ALLOCATE(expoa(1:SIZE(sweeper%expTableDat%table3D,DIM=3), &
          LBOUND(sweeper%expTableDat%table3D,DIM=2):UBOUND(sweeper%expTableDat%table3D,DIM=2)))
        ALLOCATE(expob(1:SIZE(sweeper%expTableDat%table3D,DIM=3), &
          LBOUND(sweeper%expTableDat%table3D,DIM=2):UBOUND(sweeper%expTableDat%table3D,DIM=2)))

        i4=SIZE(expoa,DIM=2)
        i5=SIZE(expoa,DIM=1)
        i3=-i4*i5*2

        ALLOCATE(expoab(i3+1:0))

        expoa=TRANSPOSE(sweeper%expTableDat%table3D(1,:,:))*0.001D0
        expob=TRANSPOSE(sweeper%expTableDat%table3D(2,:,:)) !expoa(ipol,ix1)

        DO i1=-LBOUND(expoa,DIM=2),0
          DO i2=1,SIZE(expoa,DIM=1)
            i3=i3+1
            expoab(i3)=expoa(i2,i1)
            i3=i3+1
            expoab(i3)=expob(i2,i1)
          ENDDO
        ENDDO
      ENDIF

      sweeper%sweep => MOCSolver_SweepMG
      sweeper%setExtSource => setExtSource_MOCP0

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
!      DO ig=1,sweeper%ng
!        DO ix=1,sweeper%nxsreg
!          IF(LBOUND(sweeper%myXSMesh(ix)%xsmacnf) <= ig .AND. &
!            ig <= UBOUND(sweeper%myXSMesh(ix)%xsmacnf)) THEN
!            xsnfg = sweeper%myXSMesh(ix)%xsmacnf(ig)*rkeff
!            DO ir=1,sweeper%myXSMesh(ix)%nreg
!              ireg = sweeper%myXSMesh(ix)%ireg(ir)
!              psi(ireg) = psi(ireg) + xsnfg*sweeper%phis(ireg,ig)
!            ENDDO
!          endif
!TODO: calcMacroChi?
!
!
!    END SUBROUTINE calcFissionSrc
!===============================================================================
    SUBROUTINE setExtSource_MOCP0(sweeper,source)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      CLASS(SourceType),POINTER,INTENT(IN) :: source

      NULLIFY(sweeper%qext)
      sweeper%hasSource=.FALSE.
      SELECTTYPE(source); TYPE IS(SourceType_P0)
        sweeper%mySrc => source
        sweeper%qext => source%qext
        sweeper%hasSource=.TRUE.
      ENDSELECT

    END SUBROUTINE setExtSource_MOCP0
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
!===============================================================================
    SUBROUTINE MOCSolver_Sweep1G(sweeper,ninners,tol,source,psi)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
      CLASS(SourceType),POINTER,INTENT(INOUT) :: source
      DOUBLE PRECISION,INTENT(INOUT) :: psi(:)
      ! Local Variables
      INTEGER :: i,ig
      DOUBLE PRECISION :: timeStt,timeStp,timeTotal

      timeTotal = 0.0D0

      ! Group loop.  This is actualy in FixedSrcSolver in MPACT
      DO ig=1,sweeper%ng
        ! Set up source
!         CALL source%initExtSource(ig)
!         CALL source%computeMGFS(ig,psi)
        CALL source%updateInScatter( &
          ig,sweeper%igstt,sweeper%igstp)
        CALL sweeper%setExtSource(source)

        ! This is the real beginning of the sweep routines in MPACT
        IF(1 <= ig .AND. ig <= sweeper%ng) THEN
          sweeper%activeg = ig
          sweeper%phiang1g_in => sweeper%phiang(ig)
          sweeper%phis1g = sweeper%phis(:,ig)

          DO i=1,ninners
            sweeper%nsweeps = sweeper%nsweeps + 1
            !IF(i == ninners) sweeper%sweepCur=.TRUE.
            CALL sweeper%mySrc%updateSelfScatter(ig,sweeper%qbar,sweeper%phis1g)
            CALL MOCSolver_Setup1GFSP(sweeper%myXSMesh,sweeper%nxsreg, &
              sweeper%phis1g,sweeper%nreg,sweeper%xstr,sweeper%qbar,ig)
            sweeper%phis1gd = sweeper%phis1g

            sweeper%phis1g = 0.0D0
            ! Assumes sweeper%sweepCur == .FALSE.
            CALL CPU_TIME(timeStt)
            CALL sweeper%sweep2D_prodquad(i)
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

    END SUBROUTINE MOCSolver_Sweep1G
!===============================================================================
    SUBROUTINE MOCSolver_SweepMG(sweeper,ninners,tol,source,psi)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
      CLASS(SourceType),POINTER,INTENT(INOUT) :: source
      DOUBLE PRECISION,INTENT(INOUT) :: psi(:)
      ! Local Variables
      INTEGER :: i,ig
      DOUBLE PRECISION :: timeStt,timeStp,timeTotal

      timeTotal = 0.0D0

      ! Group loop.  This is actualy in FixedSrcSolver in MPACT
      ! Set up generic external source and fission source
      CALL source%initExtSource()
      CALL source%computeMGFS(psi)

      DO i=1,ninners
        sweeper%nsweeps = sweeper%nsweeps + 1
        DO ig=1,sweeper%ng
          CALL source%updateInScatter(ig,sweeper%igstt,sweeper%igstp)
          CALL sweeper%mySrc%updateSelfScatter(ig,sweeper%qbar,sweeper%phis(:,ig))
          CALL MOCSolver_Setup1GFSP(sweeper%myXSMesh,sweeper%nxsreg, &
            sweeper%phis(:,ig),sweeper%nreg,sweeper%xstrmg(:,ig),sweeper%qbarmg(:,ig),ig)
          source%phisd(:,ig)=sweeper%phis(:,ig)
        ENDDO !ig

        CALL sweeper%setExtSource(source)

        sweeper%phis = 0.0D0

        ! Assumes sweeper%sweepCur == .FALSE.
        CALL CPU_TIME(timeStt)
        CALL sweeper%sweep2D_prodquad(i)
        CALL CPU_TIME(timeStp)
        timeStp = timeStp - timeStt
        timeTotal = timeTotal + timeStp
        WRITE(*,*) 'Iteration Time: ',timeStp
        WRITE(*,*) 'Accumulated Sweep Time: ',timeTotal
        WRITE(*,*)
      ENDDO !i

      ! Write to output file for comparison
      WRITE(125) SIZE(sweeper%phis,1)
      WRITE(125) SIZE(sweeper%phis,2)
      DO ig=1,sweeper%ng
        DO i=1,sweeper%nreg
          WRITE(125) sweeper%phis(i,ig)
        ENDDO
      ENDDO

      ! Update boundary surface flux here, if sweep Cur and associated coarse mesh
      CALL sweeper%UpdateBC%Finish()

      ! hasSource = .FALSE.

    END SUBROUTINE MOCSolver_SweepMG
!===============================================================================
    SUBROUTINE sweep2D_prodquad_P0(sweeper,i)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
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

    END SUBROUTINE sweep2D_prodquad_P0
END MODULE sweeper
