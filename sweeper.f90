MODULE sweeper
!TODO: Check out initExtSource, computMGFS, and updateInScatter routines in FixedSrcSolver/f90:876-879
!TODO: Also check out setExtSource in FixedSrcSolver.f90:897
!TODO: call sweep
!TODO: check MOCSweeper_P0.f90:468-477
!TODO: CALL sweep2D_prodquad
!TODO: dump the following:
!  nreg
!  iangstt,iangstp
!  azimuthal and polar weights
!  longRayDat%nlongrays
!  longRayDat%anglaes%longrays%ifirstModMesh, %iside, %BCIndex, %iside
!  myModMesh%imesh%ifrstfsreg
!  trmesh%rtdat%angles%%rays%nseg
!  rtmesh%rtdat%angles%rays%ireg
!  xstr
!  rtmesh%rtdat%angles%rays%hseg
!  modRayDat%angles%rays%nextsurf
!  modRayDat%angles%rays%nextray
!  myModMesh%neigh
!TODO: look into expoa and expob
!TODO: mimic sweep loops
!TODO: Look into UpdateBC%Start
!TODO: Look into final line of sweep2D_prodquad_P0 routine


END MODULE sweeper
