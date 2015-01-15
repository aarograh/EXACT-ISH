MODULE openmp

  USE sweeper
  USE sweeperUtils

  IMPLICIT NONE
  PRIVATE

  !This is the same as the base sweeper
  PUBLIC :: sweep2D_prodquad_P0_vectoripol0
  !This version vectorize the first block
  PUBLIC :: sweep2D_prodquad_P0_vectoripol1
  !Inline the polar angle dependent table loopup
  PUBLIC :: sweep2D_prodquad_P0_vectoripol1_2
  !Vectorize the second block
  PUBLIC :: sweep2D_prodquad_P0_vectoripol2
  !Vectorize the third blcok
  PUBLIC :: sweep2D_prodquad_P0_vectoripol3
  !Vectorize all three blocks with inlined expT lookup - best performance
  PUBLIC :: sweep2D_prodquad_P0_vectoripol3_2
  !With OpenMP
  PUBLIC :: sweep2D_prodquad_P0_vectoripol3_3
  !group inner
  PUBLIC :: sweep2D_prodquad_P0_GI1
  !group innermost
  PUBLIC :: sweep2D_prodquad_P0_GI2


  CONTAINS
!===============================================================================
    !This version is the same as the base sweeper
    SUBROUTINE sweep2D_prodquad_P0_vectoripol0(sweeper,i)
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
      DOUBLE PRECISION :: timerStt,timerStp

      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)
      wsum = 4.0D0*PI
      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

      ALLOCATE(phibar(sweeper%nreg))
      tphi(:,ithd) = 0.0D0

      CALL CPU_TIME(timerStt)
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
      CALL CPU_TIME(timerStp)
      WRITE(*,*) "sweep time: ", timerStp-timerStt

      DEALLOCATE(phibar)

      sweeper%phis1g = sweeper%phis1g + tphi(:,ithd)
      sweeper%phis1g = sweeper%phis1g/(sweeper%xstr*sweeper%vol/sweeper%pz) + &
        sweeper%qbar*wsum

    ENDSUBROUTINE sweep2D_prodquad_P0_vectoripol0
!===============================================================================
    !This version vectorize the first block
    SUBROUTINE sweep2D_prodquad_P0_vectoripol1(sweeper,i)
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

      DOUBLE PRECISION :: timerStt,timerStp
      DOUBLE PRECISION,ALLOCATABLE :: arg(:)

      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)
      wsum = 4.0D0*PI
      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

      ALLOCATE(phibar(sweeper%nreg))
      tphi(:,ithd) = 0.0D0

      CALL CPU_TIME(timerStt)

      DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
        wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
          sweeper%modRayDat%angquad%walpha(iang)*PI
!         DO ipol=1,npol
!           wtang(ipol) = wtangazi*sweeper%modRayDat%angquad%wtheta(ipol)* &
!             sweeper%modRayDat%angquad%sinpolang(ipol)
!         ENDDO !ipol

        !vectorizing the first block here
        wtang = wtangazi*sweeper%modRayDat%angquad%wtheta* &
          sweeper%modRayDat%angquad%sinpolang

!         !variation 1: including the bounds
!         wtang(1:npol) = wtangazi*sweeper%modRayDat%angquad%wtheta(1:npol)* &
!           sweeper%modRayDat%angquad%sinpolang(1:npol)

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
      CALL CPU_TIME(timerStp)
      WRITE(*,*) "sweep time: ", timerStp-timerStt

      DEALLOCATE(phibar)

      sweeper%phis1g = sweeper%phis1g + tphi(:,ithd)
      sweeper%phis1g = sweeper%phis1g/(sweeper%xstr*sweeper%vol/sweeper%pz) + &
        sweeper%qbar*wsum

    ENDSUBROUTINE sweep2D_prodquad_P0_vectoripol1
!===============================================================================
    !This version vectorize the first block.
    !Also it integrates the polar angle dependent exponential table.
    !And functions are inlined to reduce the cost.
    SUBROUTINE sweep2D_prodquad_P0_vectoripol1_2(sweeper,i)
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
        exparg(SIZE(sweeper%modRayDat%angquad%wtheta),sweeper%maxsegray)
      DOUBLE PRECISION,ALLOCATABLE :: phibar(:)
      TYPE(LongRayType_Base) :: ilongRay

      DOUBLE PRECISION :: timerStt,timerStp
      DOUBLE PRECISION,ALLOCATABLE :: arg(:)
      DOUBLE PRECISION :: xval
      INTEGER :: ix

      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)
      wsum = 4.0D0*PI
      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

      ALLOCATE(phibar(sweeper%nreg))
      tphi(:,ithd) = 0.0D0

      CALL CPU_TIME(timerStt)

      DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
        wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
          sweeper%modRayDat%angquad%walpha(iang)*PI
!         DO ipol=1,npol
!           wtang(ipol) = wtangazi*sweeper%modRayDat%angquad%wtheta(ipol)* &
!             sweeper%modRayDat%angquad%sinpolang(ipol)
!         ENDDO !ipol

        !vectorizing the first block here
        wtang = wtangazi*sweeper%modRayDat%angquad%wtheta* &
          sweeper%modRayDat%angquad%sinpolang

!         !variation 1: including the bounds
!         wtang(1:npol) = wtangazi*sweeper%modRayDat%angquad%wtheta(1:npol)* &
!           sweeper%modRayDat%angquad%sinpolang(1:npol)

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

!           DO ipol=1,npol
!             rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
!             DO iseg=1,nseglray
!               exparg(iseg,ipol) = sweeper%expTableDat%EXPT(tau_seg(iseg)*rpol)
!             ENDDO !iseg
!           ENDDO !ipol

          DO iseg=1,nseglray
            xval=tau_seg(iseg)*1000.0D0
            ix=xval
            ix=MAX(ix,-16000)
            exparg(:,iseg)=expoa(:,ix)*xval+expob(:,ix)
          ENDDO

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
              phid1 = phid1*exparg(ipol,iseg1)
              !phio1 stores the outgoing angular flux to be used for the next
              !segment as incoming angular flux.
              phio1(iseg1) = phio1(iseg1-1) - phid1
              phibar(ireg1) = phibar(ireg1) + phid1*wtang(ipol)

              ireg2 = irg_seg(iseg2)
              phid2 = phio2(iseg2+1) - sweeper%qbar(ireg2)
              phid2 = phid2*exparg(ipol,iseg2)
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
      CALL CPU_TIME(timerStp)
      WRITE(*,*) "sweep time: ", timerStp-timerStt

      DEALLOCATE(phibar)

      sweeper%phis1g = sweeper%phis1g + tphi(:,ithd)
      sweeper%phis1g = sweeper%phis1g/(sweeper%xstr*sweeper%vol/sweeper%pz) + &
        sweeper%qbar*wsum

    ENDSUBROUTINE sweep2D_prodquad_P0_vectoripol1_2
!===============================================================================
    !This version vectorize the second block
    SUBROUTINE sweep2D_prodquad_P0_vectoripol2(sweeper,i)
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

      DOUBLE PRECISION :: timerStt,timerStp
      DOUBLE PRECISION,ALLOCATABLE :: arg(:)

      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)
      wsum = 4.0D0*PI
      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

      ALLOCATE(phibar(sweeper%nreg))
      tphi(:,ithd) = 0.0D0

      CALL CPU_TIME(timerStt)

      DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
        wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
          sweeper%modRayDat%angquad%walpha(iang)*PI
!         DO ipol=1,npol
!           wtang(ipol) = wtangazi*sweeper%modRayDat%angquad%wtheta(ipol)* &
!             sweeper%modRayDat%angquad%sinpolang(ipol)
!         ENDDO !ipol

        !vectorizing the first block here
        wtang = wtangazi*sweeper%modRayDat%angquad%wtheta* &
          sweeper%modRayDat%angquad%sinpolang

!         !variation 1: including the bounds
!         wtang(1:npol) = wtangazi*sweeper%modRayDat%angquad%wtheta(1:npol)* &
!           sweeper%modRayDat%angquad%sinpolang(1:npol)

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

! !the following blcok uses scalar lookup.
!           DO ipol=1,npol
!             rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
!             DO iseg=1,nseglray
!               exparg(iseg,ipol) = sweeper%expTableDat%EXPT(tau_seg(iseg)*rpol)
!             ENDDO !iseg
!           ENDDO !ipol

!the following block uses vector fetch.
          ALLOCATE(arg(npol))
          DO iseg=1,nseglray
            arg=tau_seg(iseg)*sweeper%modRayDat%angquad%rsinpolang
            exparg(iseg,:) = sweeper%expTableDat%EXPT_vectoripol(arg)
          ENDDO !iseg
          DEALLOCATE(arg)

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
      CALL CPU_TIME(timerStp)
      WRITE(*,*) "sweep time: ", timerStp-timerStt

      DEALLOCATE(phibar)

      sweeper%phis1g = sweeper%phis1g + tphi(:,ithd)
      sweeper%phis1g = sweeper%phis1g/(sweeper%xstr*sweeper%vol/sweeper%pz) + &
        sweeper%qbar*wsum
    ENDSUBROUTINE sweep2D_prodquad_P0_vectoripol2
!===============================================================================
    !This version vectorize the third block
    SUBROUTINE sweep2D_prodquad_P0_vectoripol3(sweeper,i)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: i
      ! Local Variables
      INTEGER :: iang,ipol,ilray,imray,imod,im,iside,inextsurf,ifrstreg
      INTEGER :: imseg,iseg,iseg1,iseg2,ibc1,ibc2,nseglray,is1,is2,npol,ireg
      INTEGER :: irg_seg(0:sweeper%maxsegray),ithd
      INTEGER :: ifrstreg_proc,ireg1,ireg2
      DOUBLE PRECISION,ALLOCATABLE :: phid1(:),phid2(:)
      DOUBLE PRECISION :: wsum,rpol
      DOUBLE PRECISION :: wtangazi,wtang(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION,ALLOCATABLE :: phio1(:,:),phio2(:,:)
      DOUBLE PRECISION :: tphi(sweeper%nreg,1)
      DOUBLE PRECISION :: tau_seg(sweeper%maxsegray)
      DOUBLE PRECISION :: &
        exparg(SIZE(sweeper%modRayDat%angquad%wtheta),sweeper%maxsegray)
      DOUBLE PRECISION,ALLOCATABLE :: phibar(:,:)
      TYPE(LongRayType_Base) :: ilongRay

      DOUBLE PRECISION :: timerStt,timerStp,temp
      DOUBLE PRECISION,ALLOCATABLE :: arg(:)
      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)
      ALLOCATE(phid1(npol))
      ALLOCATE(phid2(npol))
      ALLOCATE(phio1(npol,0:sweeper%maxsegray),phio2(npol,1:sweeper%maxsegray+1))
      ALLOCATE(phibar(npol,sweeper%nreg))

      wsum = 4.0D0*PI
      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

      tphi(:,ithd) = 0.0D0

      CALL CPU_TIME(timerStt)

      DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
        wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
          sweeper%modRayDat%angquad%walpha(iang)*PI

        !vectorizing the first block here
        wtang = wtangazi*sweeper%modRayDat%angquad%wtheta* &
          sweeper%modRayDat%angquad%sinpolang

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
!the following blcok uses scalar lookup.
          DO ipol=1,npol
            rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
            DO iseg=1,nseglray
              exparg(ipol,iseg) = sweeper%expTableDat%EXPT(tau_seg(iseg)*rpol)
            ENDDO !iseg
          ENDDO !ipol

! !the following block uses vector fetch.
!           ALLOCATE(arg(npol))
!           DO iseg=1,nseglray
!             arg=tau_seg(iseg)*sweeper%modRayDat%angquad%rsinpolang
!             exparg(:,iseg) = sweeper%expTableDat%EXPT_vectoripol(arg)
!           ENDDO !iseg
!           DEALLOCATE(arg)

            phio1(:,0) = &
              sweeper%phiang1g_in%angle(iang)%face(is1)%angflux(:,ibc1)
            phio2(:,nseglray+1) = &
              sweeper%phiang1g_in%angle(iang)%face(is2)%angflux(:,ibc2)
            iseg2 = nseglray + 1

            DO iseg1=1,nseglray
              iseg2 = iseg2 - 1

              ireg1 = irg_seg(iseg1)
              phid1 = phio1(:,iseg1-1) - sweeper%qbar(ireg1)
              phid1 = phid1*exparg(:,iseg1)
              !phio1 stores the outgoing angular flux to be used for the next
              !segment as incoming angular flux.
              phio1(:,iseg1) = phio1(:,iseg1-1) - phid1
              phibar(:,ireg1) = phibar(:,ireg1) + phid1*wtang

              ireg2 = irg_seg(iseg2)
              phid2 = phio2(:,iseg2+1) - sweeper%qbar(ireg2)
              phid2 = phid2*exparg(:,iseg2)
              !phio1 stores the outgoing angular flux to be used for the next
              !segment as incoming angular flux.
              phio2(:,iseg2) = phio2(:,iseg2+1) - phid2
              phibar(:,ireg2) = phibar(:,ireg2) + phid2*wtang
            ENDDO !iseg

            sweeper%phiang1g_out%angle(iang)%face(is1)%angflux(:,ibc1) = &
              phio2(:,1)
            sweeper%phiang1g_out%angle(iang)%face(is2)%angflux(:,ibc2) = &
              phio1(:,nseglray)
        ENDDO !ilray

        !This is to sum over the azimuthal angles owned by the current proc
        !MPACT says it's for polar angles, which I think is not true.
        tphi(:,ithd) = tphi(:,ithd) + SUM(phibar,DIM=1)

        CALL sweeper%UpdateBC%Start(iang,sweeper%phiang1g_out,sweeper%phiang1g_in)
      ENDDO !iang

      CALL CPU_TIME(timerStp)
      WRITE(*,*) "sweep time: ", timerStp-timerStt

      DEALLOCATE(phibar)
      DEALLOCATE(phio1)
      DEALLOCATE(phio2)
      DEALLOCATE(phid1)
      DEALLOCATE(phid2)

      sweeper%phis1g = sweeper%phis1g + tphi(:,ithd)
      sweeper%phis1g = sweeper%phis1g/(sweeper%xstr*sweeper%vol/sweeper%pz) + &
        sweeper%qbar*wsum

    ENDSUBROUTINE sweep2D_prodquad_P0_vectoripol3
!===============================================================================
    !This version vectorize all three blocks
    !inline the exponetial table lookup.
    !This has given the best performance ever since
    SUBROUTINE sweep2D_prodquad_P0_vectoripol3_2(sweeper,i)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: i
      ! Local Variables
      INTEGER :: iang,ipol,ilray,imray,imod,im,iside,inextsurf,ifrstreg
      INTEGER :: imseg,iseg,iseg1,iseg2,ibc1,ibc2,nseglray,is1,is2,npol,ireg
      INTEGER :: irg_seg(0:sweeper%maxsegray),ithd
      INTEGER :: ifrstreg_proc,ireg1,ireg2
      DOUBLE PRECISION,ALLOCATABLE :: phid1(:),phid2(:)
      DOUBLE PRECISION :: wsum,rpol
      DOUBLE PRECISION :: wtangazi,wtang(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION,ALLOCATABLE :: phio1(:,:),phio2(:,:)
      DOUBLE PRECISION :: tphi(sweeper%nreg,1)
      DOUBLE PRECISION :: tau_seg(sweeper%maxsegray)
      DOUBLE PRECISION :: &
        exparg(SIZE(sweeper%modRayDat%angquad%wtheta),sweeper%maxsegray)
      DOUBLE PRECISION,ALLOCATABLE :: phibar(:,:)
      TYPE(LongRayType_Base) :: ilongRay

      DOUBLE PRECISION :: timerStt,timerStp,temp
      DOUBLE PRECISION,ALLOCATABLE :: arg(:)

      DOUBLE PRECISION :: xval
      INTEGER :: ix

      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)
      ALLOCATE(phid1(npol))
      ALLOCATE(phid2(npol))
      ALLOCATE(phio1(npol,0:sweeper%maxsegray),phio2(npol,1:sweeper%maxsegray+1))
      ALLOCATE(phibar(npol,sweeper%nreg))

      wsum = 4.0D0*PI
      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

      tphi(:,ithd) = 0.0D0

      CALL CPU_TIME(timerStt)

      DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
        wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
          sweeper%modRayDat%angquad%walpha(iang)*PI

        !vectorizing the first block here
        wtang = wtangazi*sweeper%modRayDat%angquad%wtheta* &
          sweeper%modRayDat%angquad%sinpolang

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

          DO iseg=1,nseglray
            xval=tau_seg(iseg)*1000.0D0
            ix=xval
            ix=MAX(ix,-16000)
            exparg(:,iseg)=expoa(:,ix)*xval+expob(:,ix)
          ENDDO


            phio1(:,0) = &
              sweeper%phiang1g_in%angle(iang)%face(is1)%angflux(:,ibc1)
            phio2(:,nseglray+1) = &
              sweeper%phiang1g_in%angle(iang)%face(is2)%angflux(:,ibc2)
            iseg2 = nseglray + 1

            DO iseg1=1,nseglray
              iseg2 = iseg2 - 1

              ireg1 = irg_seg(iseg1)
              phid1 = phio1(:,iseg1-1) - sweeper%qbar(ireg1)
              phid1 = phid1*exparg(:,iseg1)
              !phio1 stores the outgoing angular flux to be used for the next
              !segment as incoming angular flux.
              phio1(:,iseg1) = phio1(:,iseg1-1) - phid1
              phibar(:,ireg1) = phibar(:,ireg1) + phid1*wtang

              ireg2 = irg_seg(iseg2)
              phid2 = phio2(:,iseg2+1) - sweeper%qbar(ireg2)
              phid2 = phid2*exparg(:,iseg2)
              !phio1 stores the outgoing angular flux to be used for the next
              !segment as incoming angular flux.
              phio2(:,iseg2) = phio2(:,iseg2+1) - phid2
              phibar(:,ireg2) = phibar(:,ireg2) + phid2*wtang
            ENDDO !iseg

            sweeper%phiang1g_out%angle(iang)%face(is1)%angflux(:,ibc1) = &
              phio2(:,1)
            sweeper%phiang1g_out%angle(iang)%face(is2)%angflux(:,ibc2) = &
              phio1(:,nseglray)
        ENDDO !ilray

        !This is to sum over the azimuthal angles owned by the current proc
        !MPACT says it's for polar angles, which I think is not true.
        tphi(:,ithd) = tphi(:,ithd) + SUM(phibar,DIM=1)

        CALL sweeper%UpdateBC%Start(iang,sweeper%phiang1g_out,sweeper%phiang1g_in)
      ENDDO !iang

      CALL CPU_TIME(timerStp)
      WRITE(*,*) "sweep time: ", timerStp-timerStt

      DEALLOCATE(phibar)
      DEALLOCATE(phio1)
      DEALLOCATE(phio2)
      DEALLOCATE(phid1)
      DEALLOCATE(phid2)

      sweeper%phis1g = sweeper%phis1g + tphi(:,ithd)
      sweeper%phis1g = sweeper%phis1g/(sweeper%xstr*sweeper%vol/sweeper%pz) + &
        sweeper%qbar*wsum

    ENDSUBROUTINE sweep2D_prodquad_P0_vectoripol3_2

!===============================================================================
    !Based on the previous performances.
    !This version tried to include omp
    SUBROUTINE sweep2D_prodquad_P0_vectoripol3_3(sweeper,i)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: i

      ! Local Variables
      INTEGER :: iang,ipol,ilray,imray,imod,im,iside,inextsurf,ifrstreg
      INTEGER :: imseg,iseg,iseg1,iseg2,ibc1,ibc2,nseglray,is1,is2,npol,ireg
      INTEGER :: irg_seg(0:sweeper%maxsegray),ithd
      INTEGER :: ifrstreg_proc,ireg1,ireg2
      DOUBLE PRECISION,ALLOCATABLE :: phid1(:),phid2(:)
      DOUBLE PRECISION :: wsum,rpol
      DOUBLE PRECISION :: wtangazi,wtang(SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION,ALLOCATABLE :: phio1(:,:),phio2(:,:)
      DOUBLE PRECISION,ALLOCATABLE :: tphi(:,:)
      DOUBLE PRECISION :: tau_seg(sweeper%maxsegray)
      DOUBLE PRECISION :: &
        exparg(SIZE(sweeper%modRayDat%angquad%wtheta),sweeper%maxsegray)
      DOUBLE PRECISION,ALLOCATABLE :: phibar(:,:)
      TYPE(LongRayType_Base) :: ilongRay

      DOUBLE PRECISION :: timerStt,timerStp,temp
      DOUBLE PRECISION,ALLOCATABLE :: arg(:)

      DOUBLE PRECISION :: xval
      INTEGER :: ix,OMP_GET_THREAD_NUM,nproc
      DOUBLE PRECISION :: OMP_GET_WTIME,oStt,oStp

      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)

      wsum = 4.0D0*PI
      nproc=8
! PRINT*,sweeper%longRayDat%nlongrays

      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)
      ALLOCATE(tphi(sweeper%nreg,nproc))
      CALL CPU_TIME(timerStt)
      oStt=OMP_GET_WTIME()

! OPEN(FILE='a.out',UNIT=456)

!$OMP PARALLEL NUM_THREADS(nproc) DEFAULT(PRIVATE) &
!$OMP SHARED(sweeper,npol,wsum,ifrstreg_proc,tphi,expoa,expob,PI)
!$    ithd=OMP_GET_THREAD_NUM()+1
      ALLOCATE(phid1(npol))
      ALLOCATE(phid2(npol))
      ALLOCATE(phio1(npol,0:sweeper%maxsegray))
      ALLOCATE(phio2(npol,1:sweeper%maxsegray+1))

      ALLOCATE(phibar(npol,sweeper%nreg))
      tphi(:,ithd) = 0.0D0

      DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
        wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
          sweeper%modRayDat%angquad%walpha(iang)*PI
        !vectorizing the first block here
        wtang = wtangazi*sweeper%modRayDat%angquad%wtheta* &
          sweeper%modRayDat%angquad%sinpolang

        phibar = 0.0D0

!$OMP DO SCHEDULE(DYNAMIC,4)
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
          !Evaluate exponential function
          DO iseg=1,nseglray
            xval=tau_seg(iseg)*1000.0D0
            ix=INT(xval)
            ix=MAX(ix,-16000)
            exparg(:,iseg)=expoa(:,ix)*xval+expob(:,ix)
          ENDDO


          phio1(:,0) = &
            sweeper%phiang1g_in%angle(iang)%face(is1)%angflux(:,ibc1)
          phio2(:,nseglray+1) = &
            sweeper%phiang1g_in%angle(iang)%face(is2)%angflux(:,ibc2)
          iseg2 = nseglray + 1

          DO iseg1=1,nseglray
            iseg2 = iseg2 - 1

            ireg1 = irg_seg(iseg1)
            phid1 = phio1(:,iseg1-1) - sweeper%qbar(ireg1)
            phid1 = phid1*exparg(:,iseg1)
            !phio1 stores the outgoing angular flux to be used for the next
            !segment as incoming angular flux.
            phio1(:,iseg1) = phio1(:,iseg1-1) - phid1
            phibar(:,ireg1) = phibar(:,ireg1) + phid1*wtang

            ireg2 = irg_seg(iseg2)
            phid2 = phio2(:,iseg2+1) - sweeper%qbar(ireg2)
            phid2 = phid2*exparg(:,iseg2)
            !phio1 stores the outgoing angular flux to be used for the next
            !segment as incoming angular flux.
            phio2(:,iseg2) = phio2(:,iseg2+1) - phid2
            phibar(:,ireg2) = phibar(:,ireg2) + phid2*wtang
          ENDDO !iseg
          sweeper%phiang1g_out%angle(iang)%face(is1)%angflux(:,ibc1) = &
            phio2(:,1)
          sweeper%phiang1g_out%angle(iang)%face(is2)%angflux(:,ibc2) = &
            phio1(:,nseglray)
        ENDDO !ilray
!$OMP END DO NOWAIT

        !This is to sum over the azimuthal angles owned by the current proc
        !MPACT says it's for polar angles, which I think is not true.
        tphi(:,ithd) = tphi(:,ithd) + SUM(phibar,DIM=1)

!$OMP BARRIER

!$OMP SINGLE
        CALL sweeper%UpdateBC%Start(iang,sweeper%phiang1g_out,sweeper%phiang1g_in)
!$OMP END SINGLE
      ENDDO !iang


      DEALLOCATE(phibar)
      DEALLOCATE(phio1)
      DEALLOCATE(phio2)
      DEALLOCATE(phid1)
      DEALLOCATE(phid2)
!$OMP END PARALLEL

      CALL CPU_TIME(timerStp)
      oStp=OMP_GET_WTIME()
      WRITE(*,*) "sweep time: ", timerStp-timerStt
      WRITE(*,*) "sweep wall time: ", oStp-oStt

      !Reduce over threads
      DO ithd=1,nproc
!          PRINT*,ithd,nproc
        sweeper%phis1g=sweeper%phis1g+tphi(:,ithd)
      ENDDO

      sweeper%phis1g = sweeper%phis1g/(sweeper%xstr*sweeper%vol/sweeper%pz) + &
        sweeper%qbar*wsum
      DEALLOCATE(tphi)
    ENDSUBROUTINE sweep2D_prodquad_P0_vectoripol3_3
!===============================================================================
    !This version for the group inner design
    SUBROUTINE sweep2D_prodquad_P0_GI1(sweeper,i)
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
      DOUBLE PRECISION :: hseg(sweeper%maxsegray)
      DOUBLE PRECISION :: &
        exparg(sweeper%maxsegray,SIZE(sweeper%modRayDat%angquad%wtheta))
      DOUBLE PRECISION,ALLOCATABLE :: phibar(:)
      TYPE(LongRayType_Base) :: ilongRay

      INTEGER :: ig

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
!               tau_seg(iseg) = -sweeper%xstr(ireg)* &
!                 sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%hseg(imseg)
              hseg(iseg)=sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%hseg(imseg)
              irg_seg(iseg) = ireg
            ENDDO !imseg

            inextsurf = sweeper%modRayDat%angles(iang)%rays(imray)%nextsurf(1)
            imray = sweeper%modRayDat%angles(iang)%rays(imray)%nextray(1)
            im = sweeper%myModMesh%neigh(inextsurf,im)
          ENDDO !imod

          nseglray = iseg

!           DO ipol=1,npol
!             rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
!             DO iseg=1,nseglray
!               exparg(iseg,ipol) = sweeper%expTableDat%EXPT(tau_seg(iseg)*rpol)
!             ENDDO !iseg
!           ENDDO !ipol

          DO ig=1,sweeper%ng

            sweeper%xstr => sweeper%xstrmg(:,ig)
            sweeper%qbar => sweeper%qbarmg(:,ig)
            sweeper%activeg=ig

            DO ipol=1,npol
              rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
              DO iseg=1,nseglray
                ireg=irg_seg(iseg)
                exparg(iseg,ipol) = sweeper%expTableDat%EXPT(-sweeper%xstr(ireg)*hseg(iseg)*rpol)
              ENDDO !iseg
            ENDDO !ipol

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
                phid1 = phid1*exparg(iseg1,ipol)
                !phio1 stores the outgoing angular flux to be used for the next
                !segment as incoming angular flux.
                phio1(iseg1) = phio1(iseg1-1) - phid1
                sweeper%phis(ireg1,ig) = sweeper%phis(ireg1,ig) + phid1*wtang(ipol)

                ireg2 = irg_seg(iseg2)
                phid2 = phio2(iseg2+1) - sweeper%qbarmg(ireg2,ig)
                phid2 = phid2*exparg(iseg2,ipol)
                !phio1 stores the outgoing angular flux to be used for the next
                !segment as incoming angular flux.
                phio2(iseg2) = phio2(iseg2+1) - phid2
                sweeper%phis(ireg2,ig) = sweeper%phis(ireg2,ig) + phid2*wtang(ipol)
              ENDDO !iseg

              sweeper%phiangmg_out(ig)%angle(iang)%face(is1)%angflux(ipol,ibc1) = &
                phio2(1)
              sweeper%phiangmg_out(ig)%angle(iang)%face(is2)%angflux(ipol,ibc2) = &
                phio1(nseglray)
            ENDDO !ipol
          ENDDO !ig
        ENDDO !ilray

        !This is to sum over the angles owned by the current proc.
        !MPACT says it's for polar angles, which I think is not true.
!         tphi(:,ithd) = tphi(:,ithd) + phibar
        DO ig=1,sweeper%ng
          CALL sweeper%UpdateBC%Start(iang,sweeper%phiangmg_out(ig),sweeper%phiang(ig))
        ENDDO !ig
      ENDDO !iang

      DO ig=1,sweeper%ng
        sweeper%phis(:,ig) = sweeper%phis(:,ig)/(sweeper%xstrmg(:,ig)*sweeper%vol/sweeper%pz) + &
          sweeper%qbarmg(:,ig)*wsum
      ENDDO !ig

    ENDSUBROUTINE sweep2D_prodquad_P0_GI1
!===============================================================================
    !This version is the same as the base sweeper
    SUBROUTINE sweep2D_prodquad_P0_GI2(sweeper,i)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: i
      ! Local Variables
      INTEGER :: iang,ipol,ilray,imray,imod,im,iside,inextsurf,ifrstreg
      INTEGER :: imseg,iseg,iseg1,iseg2,ibc1,ibc2,nseglray,is1,is2,npol,ireg
      INTEGER :: irg_seg(0:sweeper%maxsegray),ithd
      INTEGER :: ifrstreg_proc,ireg1,ireg2
      DOUBLE PRECISION :: phid1,phid2,wsum,rpol
      DOUBLE PRECISION :: wtangazi,wtang(SIZE(sweeper%modRayDat%angquad%wtheta))
!       DOUBLE PRECISION :: phio1(0:sweeper%maxsegray),phio2(1:sweeper%maxsegray+1)
      DOUBLE PRECISION :: tphi(sweeper%nreg,1)
      DOUBLE PRECISION :: tau_seg(sweeper%maxsegray)
      DOUBLE PRECISION :: hseg(sweeper%maxsegray)
!       DOUBLE PRECISION :: &
!         exparg(sweeper%maxsegray,SIZE(sweeper%modRayDat%angquad%wtheta))
!       DOUBLE PRECISION,ALLOCATABLE :: phibar(:)
      TYPE(LongRayType_Base) :: ilongRay
      DOUBLE PRECISION,ALLOCATABLE :: phis(:,:)
      DOUBLE PRECISION,ALLOCATABLE :: phio1(:,:),phio2(:,:)
      DOUBLE PRECISION :: exparg1,exparg2

      INTEGER :: ig

      ithd = 1
      npol = SIZE(sweeper%modRayDat%angquad%wtheta)
      wsum = 4.0D0*PI
      ifrstreg_proc = sweeper%myModMesh%ifrstfsreg(sweeper%imeshstt)

      ALLOCATE(phio1(1:npol,1:sweeper%ng))
      ALLOCATE(phio2(1:npol,1:sweeper%ng))
      ALLOCATE(phis(sweeper%nreg,sweeper%ng))

      DO iang=sweeper%modRayDat%iangstt,sweeper%modRayDat%iangstp
        wtangazi = sweeper%modRayDat%angles(iang)%dlr* &
          sweeper%modRayDat%angquad%walpha(iang)*PI
        DO ipol=1,npol
          wtang(ipol) = wtangazi*sweeper%modRayDat%angquad%wtheta(ipol)* &
            sweeper%modRayDat%angquad%sinpolang(ipol)
        ENDDO !ipol

        phis = 0.0D0

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
!               tau_seg(iseg) = -sweeper%xstr(ireg)* &
!                 sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%hseg(imseg)
              hseg(iseg)=sweeper%rtmesh(im)%rtdat%angles(iang)%rays(imray)%hseg(imseg)
              irg_seg(iseg) = ireg
            ENDDO !imseg

            inextsurf = sweeper%modRayDat%angles(iang)%rays(imray)%nextsurf(1)
            imray = sweeper%modRayDat%angles(iang)%rays(imray)%nextray(1)
            im = sweeper%myModMesh%neigh(inextsurf,im)
          ENDDO !imod

          nseglray = iseg

          iseg2=nseglray+1
          DO ig=1,sweeper%ng
            DO ipol=1,npol
              phio1(ipol,ig)=sweeper%phiang(ig)%angle(iang)%face(is1)%angflux(ipol,ibc1)
              phio2(ipol,ig)=sweeper%phiang(ig)%angle(iang)%face(is2)%angflux(ipol,ibc2)
            ENDDO !ipol
          ENDDO !ig

          DO iseg1=1,nseglray
            ireg1=irg_seg(iseg1)
            iseg2=iseg2-1
            ireg2=irg_seg(iseg2)
            DO ig=1,sweeper%ng
              DO ipol=1,npol
                rpol=sweeper%modRayDat%angquad%rsinpolang(ipol)
                exparg1=sweeper%expTableDat%EXPT(-sweeper%xstrmg(ireg1,ig)*hseg(iseg1)*rpol)
                exparg2=sweeper%expTableDat%EXPT(-sweeper%xstrmg(ireg2,ig)*hseg(iseg2)*rpol)
                !forward direction
                phid1=phio1(ipol,ig)-sweeper%qbarmg(ireg1,ig)
                phid1=phid1*exparg1
                phio1(ipol,ig)=phio1(ipol,ig)-phid1
                phis(ireg1,ig)=phis(ireg1,ig)+phid1*wtang(ipol)
                !backward direction
                phid2=phio2(ipol,ig)-sweeper%qbarmg(ireg2,ig)
                phid2=phid2*exparg2
                phio2(ipol,ig)=phio2(ipol,ig)-phid2
                phis(ireg2,ig)=phis(ireg2,ig)+phid2*wtang(ipol)
              ENDDO !ipol
            ENDDO !ig
          ENDDO !iseg

          !The following part could be avoided if we could change the
          !interface of sweeper%UpdateBC%Start(iang,
          !                         outgoing%angle(iang)%face(iface)%angflux,
          !                         incoming%angle(irefl)%face(iface)%angflux)
          DO ig=1,sweeper%ng
            DO ipol=1,npol
              sweeper%phiangmg_out(ig)%angle(iang)%face(is1)%angflux(ipol,ibc1) = &
                phio2(ipol,ig)
              sweeper%phiangmg_out(ig)%angle(iang)%face(is2)%angflux(ipol,ibc2) = &
                phio1(ipol,ig)
            ENDDO !ipol
          ENDDO !ig
        ENDDO !ilray

        DO ig=1,sweeper%ng
          CALL sweeper%UpdateBC%Start(iang,sweeper%phiangmg_out(ig),sweeper%phiang(ig))
        ENDDO !ig
      ENDDO !iang

      sweeper%phis=phis
      DEALLOCATE(phis)
      DEALLOCATE(phio1)
      DEALLOCATE(phio2)

      DO ig=1,sweeper%ng
        sweeper%phis(:,ig) = sweeper%phis(:,ig)/(sweeper%xstrmg(:,ig)*sweeper%vol/sweeper%pz) + &
          sweeper%qbarmg(:,ig)*wsum
      ENDDO !ig

    ENDSUBROUTINE sweep2D_prodquad_P0_GI2
ENDMODULE openmp
