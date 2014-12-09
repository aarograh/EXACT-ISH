MODULE openmp

  USE sweeper
  USE sweeperUtils

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sweep2D_prodquad_P0_vectoripol


  CONTAINS
!===============================================================================
    SUBROUTINE sweep2D_prodquad_P0_vectoripol(sweeper,i)
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
        !vector: wtang(:), wtheta(:), %sinpolang(:) of size npol
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

!the following blcok uses scalar fetch, but one sweep takes 0.08 sec.
          DO ipol=1,npol
            rpol = sweeper%modRayDat%angquad%rsinpolang(ipol)
            DO iseg=1,nseglray
              exparg(iseg,ipol) = sweeper%expTableDat%EXPT(tau_seg(iseg)*rpol)
            ENDDO !iseg
          ENDDO !ipol

! !the following block uses vector fetch, but one sweep now takes 0.11 sec.
! ALLOCATE(arg(npol))
!           DO iseg=1,nseglray
!             arg=tau_seg(iseg)*sweeper%modRayDat%angquad%rsinpolang
!             exparg(iseg,:) = sweeper%expTableDat%EXPT_vectoripol(arg)
!           ENDDO !iseg
! DEALLOCATE(arg)

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

      CALL CPU_TIME(timerStp)
      WRITE(*,*) "sweep time: ", timerStp-timerStt
    ENDSUBROUTINE sweep2D_prodquad_P0_vectoripol

END MODULE openmp
