MODULE  MF_Euler_dgDiscretizationModule

  ! --- AMReX Modules ---

  USE amrex_amrcore_module, ONLY: &
    amrex_get_finest_level
  USE amrex_amr_module, ONLY: &
    amrex_geom
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum, &
    amrex_parallel_communicator

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nDF
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit, &
    OffGridFlux_Euler_X1_Inner, &
    OffGridFlux_Euler_X1_Outer, &
    OffGridFlux_Euler_X2_Inner, &
    OffGridFlux_Euler_X2_Outer, &
    OffGridFlux_Euler_X3_Inner, &
    OffGridFlux_Euler_X3_Outer
  USE Euler_DiscontinuityDetectionModule, ONLY: &
    DetectShocks_Euler
  USE Euler_MeshRefinementModule, ONLY: &
    FaceRatio

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    thornado2amrex_X_F, &
    AllocateArray_X, &
    DeallocateArray_X
  USE MF_FieldsModule_Euler, ONLY: &
    FluxRegister_Euler, &
    OffGridFlux_Euler_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    swX, &
    UseFluxCorrection_Euler, &
    DEBUG
  USE FillPatchModule, ONLY: &
    FillPatch
  USE AverageDownModule, ONLY: &
    AverageDown

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_Euler_MF

  INTERFACE ComputeIncrement_Euler_MF
    MODULE PROCEDURE ComputeIncrement_Euler_MF_SingleLevel
    MODULE PROCEDURE ComputeIncrement_Euler_MF_MultipleLevels
  END INTERFACE ComputeIncrement_Euler_MF

CONTAINS


  SUBROUTINE ComputeIncrement_Euler_MF_MultipleLevels &
    ( MF_uGF, MF_uCF, MF_uDF, MF_duCF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF (0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF (0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF (0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCF(0:)

    INTEGER :: iLevel, iErr

    OffGridFlux_Euler_MF = Zero

    DO iLevel = 0, nLevels-1

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        WRITE(*,'(2x,A,I3.3)') &
          'CALL ComputeIncrement_Euler_MF_SingleLevel, iLevel: ', &
          iLevel

      END IF

      CALL ComputeIncrement_Euler_MF_SingleLevel &
             ( iLevel, MF_uGF, MF_uCF, MF_uDF, MF_duCF(iLevel) )

    END DO

    CALL AverageDown( MF_uGF, MF_duCF )

  END SUBROUTINE ComputeIncrement_Euler_MF_MultipleLevels


  SUBROUTINE ComputeIncrement_Euler_MF_SingleLevel &
    ( iLevel, MF_uGF, MF_uCF, MF_uDF, MF_duCF )

!    DO iLevel = 0, nLevels-1
!
!      ! --- Apply boundary conditions to interior domains ---
!
!      CALL FillPatch( iLevel, MF_uGF )
!      CALL FillPatch( iLevel, MF_uGF, MF_uDF )
!      CALL FillPatch &
!             ( iLevel, MF_uGF, MF_uCF, &
!               MF_uDF, ApplyPositivityLimiter_Option = .TRUE. )
!
!      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )
!
!      DO WHILE( MFI % next() )
!
!        uGF  => MF_uGF(iLevel) % DataPtr( MFI )
!        uCF  => MF_uCF(iLevel) % DataPtr( MFI )
!        uDF  => MF_uDF(iLevel) % DataPtr( MFI )
!
!        iLo_MF = LBOUND( uGF )
!
!        BX = MFI % tilebox()
!
!        iX_B0 = BX % lo
!        iX_E0 = BX % hi
!        iX_B1 = BX % lo - swX
!        iX_E1 = BX % hi + swX
!
!        CALL AllocateArray_X &
!               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
!                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
!                 G )
!
!        CALL AllocateArray_X &
!               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
!                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
!                 U )
!
!        CALL AllocateArray_X &
!               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
!                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDF ], &
!                 D )
!
!        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )
!
!        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )
!
!        CALL amrex2thornado_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )
!
!        ! --- Apply boundary conditions to physical boundaries ---
!
!        CALL ConstructEdgeMap( iLevel, BX, Edge_Map )
!
!        CALL ApplyBoundaryConditions_Euler_MF &
!               ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )
!
!        CALL DetectShocks_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )
!
!        CALL thornado2amrex_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )
!
!        CALL DeallocateArray_X &
!               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
!                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDF ], &
!                 D )
!
!        CALL DeallocateArray_X &
!               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
!                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
!                 U )
!
!        CALL DeallocateArray_X &
!               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
!                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
!                 G )
!
!      END DO
!
!      CALL amrex_mfiter_destroy( MFI )
!
!    END DO

    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCF

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: duCF    (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X1(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X2(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X3(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G             (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U             (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D             (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU            (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X1(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X2(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X3(:,:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    TYPE(amrex_multifab) :: SurfaceFluxes(1:nDimsX)
    INTEGER              :: iDimX, nGhost(nDimsX), nDOFX_X(3)
    LOGICAL              :: Nodal(nDimsX)

    TYPE(EdgeMap) :: Edge_Map

    ! --- Maybe don't need to apply boundary conditions since
    !     they're applied in the shock detector ---

    ! --- Apply boundary conditions to interior domains ---

    CALL FillPatch( iLevel, MF_uGF )
    CALL FillPatch( iLevel, MF_uGF, MF_uDF )
    CALL FillPatch &
           ( iLevel, MF_uGF, MF_uCF, &
             MF_uDF, ApplyPositivityLimiter_Option = .TRUE. )

    CALL MF_duCF % SetVal( Zero )

    CALL CreateMesh_MF( iLevel, MeshX )

    nDOFX_X(1) = nDOFX_X1
    nDOFX_X(2) = nDOFX_X2
    nDOFX_X(3) = nDOFX_X3

    nGhost = 0

    DO iDimX = 1, nDimsX

      Nodal        = .FALSE.
      Nodal(iDimX) = .TRUE.

      CALL amrex_multifab_build &
             ( SurfaceFluxes(iDimX), &
               MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX_X(iDimX) * nCF, nGhost, Nodal )

      CALL SurfaceFluxes(iDimX) % SetVal( Zero )

    END DO

! #if defined( THORNADO_OMP )
!     !$OMP PARALLEL &
!     !$OMP PRIVATE( MFI, BX, uGF, uCF, uDF, duCF, G, U, D, dU, &
!     !$OMP          uSurfaceFlux_X1, uSurfaceFlux_X2, uSurfaceFlux_X3, &
!     !$OMP           SurfaceFlux_X1,  SurfaceFlux_X2,  SurfaceFlux_X3, &
!     !$OMP           iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, Edge_Map )
! #endif

    CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF  => MF_uGF(iLevel) % DataPtr( MFI )
      uCF  => MF_uCF(iLevel) % DataPtr( MFI )
      uDF  => MF_uDF(iLevel) % DataPtr( MFI )
      duCF => MF_duCF        % DataPtr( MFI )

      uSurfaceFlux_X1 => SurfaceFluxes(1) % DataPtr( MFI )
      IF( nDimsX .GT. 1 ) uSurfaceFlux_X2 => SurfaceFluxes(2) % DataPtr( MFI )
      IF( nDimsX .GT. 2 ) uSurfaceFlux_X3 => SurfaceFluxes(3) % DataPtr( MFI )

      iLo_MF = LBOUND( uGF )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDF ], &
               D )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               dU )

      CALL AllocateArray_X &
             ( [ 1       , iX_B0(1)  , iX_B0(2), iX_B0(3), 1   ], &
               [ nDOFX_X1, iX_E0(1)+1, iX_E0(2), iX_E0(3), nCF ], &
               SurfaceFlux_X1 )

      CALL AllocateArray_X &
             ( [ 1       , iX_B0(1), iX_B0(2)  , iX_B0(3), 1   ], &
               [ nDOFX_X2, iX_E0(1), iX_E0(2)+1, iX_E0(3), nCF ], &
               SurfaceFlux_X2 )

      CALL AllocateArray_X &
             ( [ 1       , iX_B0(1), iX_B0(2), iX_B0(3)  , 1   ], &
               [ nDOFX_X3, iX_E0(1), iX_E0(2), iX_E0(3)+1, nCF ], &
               SurfaceFlux_X3 )

      CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

      CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

      CALL amrex2thornado_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

      ! --- Apply boundary conditions to physical boundaries ---

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL ComputeIncrement_Euler_DG_Explicit &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
               SuppressBC_Option = .TRUE., &
               SurfaceFlux_X1_Option = SurfaceFlux_X1, &
               SurfaceFlux_X2_Option = SurfaceFlux_X2, &
               SurfaceFlux_X3_Option = SurfaceFlux_X3 )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, duCF, dU )

      CALL IncrementOffGridTally_Euler( iLevel, iX_B0, iX_E0 )

      iLo_MF = LBOUND( uSurfaceFlux_X1 )

      CALL thornado2amrex_X_F &
             ( nDOFX_X1, nCF, [ iX_B0(1)  , iX_B0(2), iX_B0(3) ], &
                              [ iX_E0(1)+1, iX_E0(2), iX_E0(3) ], iLo_MF, &
                              [ iX_B0(1)  , iX_B0(2), iX_B0(3) ], &
                              [ iX_E0(1)+1, iX_E0(2), iX_E0(3) ], &
                              uSurfaceFlux_X1, SurfaceFlux_X1 )

      IF( nDimsX .GT. 1 )THEN

        iLo_MF = LBOUND( uSurfaceFlux_X2 )

        CALL thornado2amrex_X_F &
               ( nDOFX_X2, nCF, [ iX_B0(1), iX_B0(2)  , iX_B0(3) ], &
                                [ iX_E0(1), iX_E0(2)+1, iX_E0(3) ], iLo_MF, &
                                [ iX_B0(1), iX_B0(2)  , iX_B0(3) ], &
                                [ iX_E0(1), iX_E0(2)+1, iX_E0(3) ], &
                                uSurfaceFlux_X2, SurfaceFlux_X2 )

      END IF

      IF( nDimsX .GT. 2 )THEN

        iLo_MF = LBOUND( uSurfaceFlux_X3 )

        CALL thornado2amrex_X_F &
               ( nDOFX_X3, nCF, [ iX_B0(1), iX_B0(2), iX_B0(3)   ], &
                                [ iX_E0(1), iX_E0(2), iX_E0(3)+1 ], iLo_MF, &
                                [ iX_B0(1), iX_B0(2), iX_B0(3)   ], &
                                [ iX_E0(1), iX_E0(2), iX_E0(3)+1 ], &
                                uSurfaceFlux_X3, SurfaceFlux_X3 )

      END IF

      CALL DeallocateArray_X &
             ( [ 1       , iX_B0(1), iX_B0(2), iX_B0(3)  , 1   ], &
               [ nDOFX_X3, iX_E0(1), iX_E0(2), iX_E0(3)+1, nCF ], &
               SurfaceFlux_X3 )

      CALL DeallocateArray_X &
             ( [ 1       , iX_B0(1), iX_B0(2)  , iX_B0(3), 1   ], &
               [ nDOFX_X2, iX_E0(1), iX_E0(2)+1, iX_E0(3), nCF ], &
               SurfaceFlux_X2 )

      CALL DeallocateArray_X &
             ( [ 1       , iX_B0(1)  , iX_B0(2), iX_B0(3), 1   ], &
               [ nDOFX_X1, iX_E0(1)+1, iX_E0(2), iX_E0(3), nCF ], &
               SurfaceFlux_X1 )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               dU )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDF ], &
               D )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
               U )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
               G )

    END DO ! MFI

    CALL amrex_mfiter_destroy( MFI )

! #if defined( THORNADO_OMP )
!     !$OMP END PARALLEL
! #endif

    CALL amrex_parallel_reduce_sum( OffGridFlux_Euler_MF(:,iLevel), nCF )

#if defined( THORNADO_USE_MESHREFINEMENT )

    IF( UseFluxCorrection_Euler )THEN

      IF( iLevel .GT. 0 ) &
        CALL FluxRegister_Euler &
               ( iLevel   ) % FineAdd_DG ( SurfaceFluxes, nCF, FaceRatio )

      IF( iLevel .LT. amrex_get_finest_level() ) &
        CALL FluxRegister_Euler( iLevel+1 ) % CrseInit_DG( SurfaceFluxes, nCF )

    END IF ! UseFluxCorrection_Euler

#endif

    DO iDimX = 1, nDimsX

      CALL amrex_multifab_destroy( SurfaceFluxes(iDimX) )

    END DO

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ComputeIncrement_Euler_MF_SingleLevel


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE IncrementOffGridTally_Euler( iLevel, iX_B0, iX_E0 )

    INTEGER , INTENT(in) :: iLevel
    INTEGER , INTENT(in) :: iX_B0(3), iX_E0(3)

    ! --- dM = Minterior - Minitial + ( OffGrid_Outer - OffGrid_Inner ) ---

    IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo(1) ) &
      OffGridFlux_Euler_MF(:,iLevel) &
        = OffGridFlux_Euler_MF(:,iLevel) - OffGridFlux_Euler_X1_Inner
    IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi(1) ) &
      OffGridFlux_Euler_MF(:,iLevel) &
        = OffGridFlux_Euler_MF(:,iLevel) + OffGridFlux_Euler_X1_Outer

    IF( iX_B0(2) .EQ. amrex_geom(iLevel) % domain % lo(2) ) &
      OffGridFlux_Euler_MF(:,iLevel) &
        = OffGridFlux_Euler_MF(:,iLevel) - OffGridFlux_Euler_X2_Inner
    IF( iX_E0(2) .EQ. amrex_geom(iLevel) % domain % hi(2) ) &
      OffGridFlux_Euler_MF(:,iLevel) &
        = OffGridFlux_Euler_MF(:,iLevel) + OffGridFlux_Euler_X2_Outer

    IF( iX_B0(3) .EQ. amrex_geom(iLevel) % domain % lo(3) ) &
      OffGridFlux_Euler_MF(:,iLevel) &
        = OffGridFlux_Euler_MF(:,iLevel) - OffGridFlux_Euler_X3_Inner
    IF( iX_E0(3) .EQ. amrex_geom(iLevel) % domain % hi(3) ) &
      OffGridFlux_Euler_MF(:,iLevel) &
        = OffGridFlux_Euler_MF(:,iLevel) + OffGridFlux_Euler_X3_Outer

  END SUBROUTINE IncrementOffGridTally_Euler


END MODULE MF_Euler_dgDiscretizationModule
