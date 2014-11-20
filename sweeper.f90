MODULE sweeper
!TODO: Check out initExtSource, computMGFS, and updateInScatter routines in FixedSrcSolver/f90:876-879
!TODO: Also check out setExtSource in FixedSrcSolver.f90:897
!TODO: Need to mimic MOCSolver_Setup1GFSP in MOCSweeper_P0.f90:472
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
!TODO: Add timers for these loops

  USE sweeperUtils

  PUBLIC :: sweeperType
    INTEGER :: ninners=0
    INTEGER :: ngroups=0
    INTEGER :: activeg=0
    INTEGER :: nsweeps=0
    INTEGER :: nxsreg=0
    INTEGER :: nreg=0
    DOUBLE PRECISION,ALLOCATABLE :: phis1g(:)
    DOUBLE PRECISION,ALLOCATABLE :: phis1gd(:)
    DOUBLE PRECISION,ALLOCATABLE :: phis(:,:)
    DOUBLE PRECISION,POINTER :: xstr(:) => NULL()
    DOUBLE PRECISION,POINTER :: vol(:) => NULL()
    DOUBLE PRECISION,POINTER :: qbar(:) => NULL()
    TYPE(AngFluxBC),POINTER :: phiang1g_in => NULL() 
    TYPE(AngFluxBC) :: phiang1g_out
    TYPE(AngFluxBC),TARGET,ALLOCATABLE :: phiang(:) 
    TYPE(SourceType_P0),POINTER :: mySrc => NULL()
    TYPE(ModMeshType),POINTER :: myModMesh => NULL() 
    TYPE(ModularRayType),POINTER :: modRayDat
    TYPE(CoreLongRayType),SAVE :: longRayDat 
    TYPE(ModMeshRayPtrArryType),POINTER :: rtmesh(:) => NULL()
!    TYPE(UpdateBCType_MOC) :: updateBC !maybe, UpdateBC_MOC.f90: define this, %Start() and %Finish() methods.  Might be an MPI thing that I don't need
    PROCEDURE(absintf_sweep),POINTER :: sweep => NULL()
  TYPE :: sweeperType
  END TYPE sweeperType

  ABSTRACT INTERFACE
    SUBROUTINE absintf_sweep(sweeper,igroup,ninners,tol)
      IMPORT sweeperType
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: igroup
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
    END SUBROUTINE absintf_sweep
  END INTERFACE

  CONTAINS
!===============================================================================
    SUBROUTINE sweepP0_1G(sweeper,igroup,ninners,tol)
      CLASS(sweeperType),INTENT(INOUT) :: sweeper
      INTEGER,INTENT(IN) :: igroup
      INTEGER,INTENT(IN) :: ninners
      DOUBLE PRECISION,INTENT(IN) :: tol
    END SUBROUTINE sweepP0_1G
END MODULE sweeper
