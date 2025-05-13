MODULE  MF_MHD_dgDiscretizationModule

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
    nDimsX, &
    swX
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3
  USE MeshModule, ONLY: &
    MeshX
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    nDM
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE MHD_DiscretizationModule_Relativistic, ONLY: &
    ComputeIncrement_MHD_DG_Explicit, &
    OffGridFlux_MHD_X1_Inner, &
    OffGridFlux_MHD_X1_Outer, &
    OffGridFlux_MHD_X2_Inner, &
    OffGridFlux_MHD_X2_Outer, &
    OffGridFlux_MHD_X3_Inner, &
    OffGridFlux_MHD_X3_Outer
  USE MHD_DiscontinuityDetectionModule, ONLY: &
    DetectShocks_MHD
  USE MHD_MeshRefinementModule, ONLY: &
    FaceRatio, &
    WeightsX_X1c, &
    WeightsX_X2c, &
    WeightsX_X3c, &
    nFineX_X1, &
    nFineX_X2, &
    nFineX_X3, &
    vpLX_X1_Refined, &
    vpLX_X2_Refined, &
    vpLX_X3_Refined

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
  USE MF_FieldsModule_MHD, ONLY: &
    FluxRegister_MHD, &
    OffGridFlux_MHD_MF
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD_MF
  USE MF_MHD_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_MHD_MF
  USE MF_EdgeMapModule_MHD, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    UseFluxCorrection_MHD, &
    IsPeriodic, &
    DEBUG, &
    EvolveOnlyMagnetic, &
    UseDivergenceCleaning, &
    DampingParameter, &
    CleaningSpeed, &
    UsePowellSource
  USE FillPatchModule_MHD, ONLY: &
    FillPatch
  USE AverageDownModule_MHD, ONLY: &
    AverageDown

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_MHD_MF

  INTERFACE ComputeIncrement_MHD_MF
    MODULE PROCEDURE ComputeIncrement_MHD_MF_SingleLevel
    MODULE PROCEDURE ComputeIncrement_MHD_MF_MultipleLevels
  END INTERFACE ComputeIncrement_MHD_MF

CONTAINS


  SUBROUTINE ComputeIncrement_MHD_MF_MultipleLevels &
    ( t, CFL, MF_uGF, MF_uCM, MF_uDM, MF_duCM )

    REAL(DP),             INTENT(in   ) :: t      (0:)
    REAL(DP),             INTENT(in   ) :: CFL
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF (0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM (0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM (0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCM(0:)

    INTEGER :: iLevel, iErr

    OffGridFlux_MHD_MF = Zero

    DO iLevel = 0, nLevels-1

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        WRITE(*,'(2x,A,I3.3)') &
          'CALL ComputeIncrement_MHD_MF_SingleLevel, iLevel: ', &
          iLevel

      END IF

      CALL ComputeIncrement_MHD_MF_SingleLevel &
             ( t(iLevel), CFL, iLevel, MF_uGF, MF_uCM, MF_uDM, MF_duCM(iLevel) )

    END DO

    CALL AverageDown( MF_uGF, MF_duCM )

  END SUBROUTINE ComputeIncrement_MHD_MF_MultipleLevels


  SUBROUTINE ComputeIncrement_MHD_MF_SingleLevel &
    ( t, CFL, iLevel, MF_uGF, MF_uCM, MF_uDM, MF_duCM )

    REAL(DP)            , INTENT(in)    :: t
    REAL(DP)            , INTENT(in)    :: CFL
    INTEGER,              INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCM

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDM     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: duCM    (:,:,:,:)
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

    ! --- Apply boundary conditions to interior domains ---

    CALL FillPatch( iLevel, MF_uGF )
    CALL FillPatch( iLevel, MF_uGF, MF_uDM )
    CALL FillPatch( iLevel, MF_uGF, MF_uCM )
    CALL ApplyPositivityLimiter_MHD_MF &
           ( t, iLevel, MF_uGF(iLevel), MF_uCM(iLevel), MF_uDM(iLevel), &
             swX_Option = swX )

    CALL MF_duCM % SetVal( Zero )

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
               nDOFX_X(iDimX) * nCM, nGhost, Nodal )

      CALL SurfaceFluxes(iDimX) % SetVal( Zero )

    END DO

! #if defined( THORNADO_OMP )
!     !$OMP PARALLEL &
!     !$OMP PRIVATE( MFI, BX, uGF, uCM, uDM, duCM, G, U, D, dU, &
!     !$OMP          uSurfaceFlux_X1, uSurfaceFlux_X2, uSurfaceFlux_X3, &
!     !$OMP           SurfaceFlux_X1,  SurfaceFlux_X2,  SurfaceFlux_X3, &
!     !$OMP           iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, Edge_Map )
! #endif

    CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF  => MF_uGF(iLevel) % DataPtr( MFI )
      uCM  => MF_uCM(iLevel) % DataPtr( MFI )
      uDM  => MF_uDM(iLevel) % DataPtr( MFI )
      duCM => MF_duCM        % DataPtr( MFI )

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
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               U )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDM ], &
               D )

      CALL AllocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               dU )

      CALL AllocateArray_X &
             ( [ 1       , iX_B0(1)  , iX_B0(2), iX_B0(3), 1   ], &
               [ nDOFX_X1, iX_E0(1)+1, iX_E0(2), iX_E0(3), nCM ], &
               SurfaceFlux_X1 )

      CALL AllocateArray_X &
             ( [ 1       , iX_B0(1), iX_B0(2)  , iX_B0(3), 1   ], &
               [ nDOFX_X2, iX_E0(1), iX_E0(2)+1, iX_E0(3), nCM ], &
               SurfaceFlux_X2 )

      CALL AllocateArray_X &
             ( [ 1       , iX_B0(1), iX_B0(2), iX_B0(3)  , 1   ], &
               [ nDOFX_X3, iX_E0(1), iX_E0(2), iX_E0(3)+1, nCM ], &
               SurfaceFlux_X3 )

      CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

      CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCM, U )

      CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

      ! --- Apply boundary conditions to physical boundaries ---

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( t, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL DetectShocks_MHD( iX_B0, iX_E0, iX_B1, iX_E1, G, U, &
                             EvolveOnlyMagnetic, D )

      CALL ComputeIncrement_MHD_DG_Explicit &
             ( t, CFL, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
               SuppressBC_Option = .TRUE., &
               EvolveOnlyMagnetic_Option = EvolveOnlyMagnetic, &
               UseDivergenceCleaning_Option = UseDivergenceCleaning, &
               CleaningSpeed_Option = CleaningSpeed, &
               DampingParameter_Option = DampingParameter, &
               UsePowellSource_Option = UsePowellSource, &
               SurfaceFlux_X1_Option = SurfaceFlux_X1, &
               SurfaceFlux_X2_Option = SurfaceFlux_X2, &
               SurfaceFlux_X3_Option = SurfaceFlux_X3 )

      CALL thornado2amrex_X &
             ( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

      CALL thornado2amrex_X &
             ( nCM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, duCM, dU )

      CALL IncrementOffGridTally_MHD( iLevel, iX_B0, iX_E0 )

      iLo_MF = LBOUND( uSurfaceFlux_X1 )

      CALL thornado2amrex_X_F &
             ( nDOFX_X1, nCM, [ iX_B0(1)  , iX_B0(2), iX_B0(3) ], &
                              [ iX_E0(1)+1, iX_E0(2), iX_E0(3) ], iLo_MF, &
                              [ iX_B0(1)  , iX_B0(2), iX_B0(3) ], &
                              [ iX_E0(1)+1, iX_E0(2), iX_E0(3) ], &
                              uSurfaceFlux_X1, SurfaceFlux_X1 )

      IF( nDimsX .GT. 1 )THEN

        iLo_MF = LBOUND( uSurfaceFlux_X2 )

        CALL thornado2amrex_X_F &
               ( nDOFX_X2, nCM, [ iX_B0(1), iX_B0(2)  , iX_B0(3) ], &
                                [ iX_E0(1), iX_E0(2)+1, iX_E0(3) ], iLo_MF, &
                                [ iX_B0(1), iX_B0(2)  , iX_B0(3) ], &
                                [ iX_E0(1), iX_E0(2)+1, iX_E0(3) ], &
                                uSurfaceFlux_X2, SurfaceFlux_X2 )

      END IF

      IF( nDimsX .GT. 2 )THEN

        iLo_MF = LBOUND( uSurfaceFlux_X3 )

        CALL thornado2amrex_X_F &
               ( nDOFX_X3, nCM, [ iX_B0(1), iX_B0(2), iX_B0(3)   ], &
                                [ iX_E0(1), iX_E0(2), iX_E0(3)+1 ], iLo_MF, &
                                [ iX_B0(1), iX_B0(2), iX_B0(3)   ], &
                                [ iX_E0(1), iX_E0(2), iX_E0(3)+1 ], &
                                uSurfaceFlux_X3, SurfaceFlux_X3 )

      END IF

      CALL DeallocateArray_X &
             ( [ 1       , iX_B0(1), iX_B0(2), iX_B0(3)  , 1   ], &
               [ nDOFX_X3, iX_E0(1), iX_E0(2), iX_E0(3)+1, nCM ], &
               SurfaceFlux_X3 )

      CALL DeallocateArray_X &
             ( [ 1       , iX_B0(1), iX_B0(2)  , iX_B0(3), 1   ], &
               [ nDOFX_X2, iX_E0(1), iX_E0(2)+1, iX_E0(3), nCM ], &
               SurfaceFlux_X2 )

      CALL DeallocateArray_X &
             ( [ 1       , iX_B0(1)  , iX_B0(2), iX_B0(3), 1   ], &
               [ nDOFX_X1, iX_E0(1)+1, iX_E0(2), iX_E0(3), nCM ], &
               SurfaceFlux_X1 )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
               dU )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDM ], &
               D )

      CALL DeallocateArray_X &
             ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
               [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
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

    CALL amrex_parallel_reduce_sum( OffGridFlux_MHD_MF(:,iLevel), nCM )

#if defined( THORNADO_USE_MESHREFINEMENT )

    IF( UseFluxCorrection_MHD )THEN

      IF( iLevel .GT. 0 ) &
        CALL FluxRegister_MHD(iLevel  ) &
               % FineAdd_DG &
                  ( SurfaceFluxes, nCM, FaceRatio, &
                    nDOFX_X1, nDOFX_X2, nDOFX_X3, &
                    nFineX_X1, nFineX_X2, nFineX_X3, &
                    WeightsX_X1c,  WeightsX_X2c, WeightsX_X3c, &
                    vpLX_X1_Refined, vpLX_X2_Refined, vpLX_X3_Refined )

      IF( iLevel .LT. amrex_get_finest_level() ) &
        CALL FluxRegister_MHD(iLevel+1) &
               % CrseInit_DG &
                   ( SurfaceFluxes, nCM, &
                     nDOFX_X1, nDOFX_X2, nDOFX_X3, &
                     WeightsX_X1c,  WeightsX_X2c, WeightsX_X3c )

    END IF ! UseFluxCorrection_MHD

#endif

    DO iDimX = 1, nDimsX

      CALL amrex_multifab_destroy( SurfaceFluxes(iDimX) )

    END DO

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ComputeIncrement_MHD_MF_SingleLevel


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE IncrementOffGridTally_MHD( iLevel, iX_B0, iX_E0 )

    INTEGER , INTENT(in) :: iLevel
    INTEGER , INTENT(in) :: iX_B0(3), iX_E0(3)

    ! --- dM = Minterior - Minitial + ( OffGrid_Outer - OffGrid_Inner ) ---

    IF( .NOT. IsPeriodic(1) )THEN

      IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo(1) ) &
        OffGridFlux_MHD_MF(:,iLevel) &
          = OffGridFlux_MHD_MF(:,iLevel) - OffGridFlux_MHD_X1_Inner
      IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi(1) ) &
        OffGridFlux_MHD_MF(:,iLevel) &
          = OffGridFlux_MHD_MF(:,iLevel) + OffGridFlux_MHD_X1_Outer

    END IF

    IF( .NOT. IsPeriodic(2) )THEN

      IF( iX_B0(2) .EQ. amrex_geom(iLevel) % domain % lo(2) ) &
        OffGridFlux_MHD_MF(:,iLevel) &
          = OffGridFlux_MHD_MF(:,iLevel) - OffGridFlux_MHD_X2_Inner
      IF( iX_E0(2) .EQ. amrex_geom(iLevel) % domain % hi(2) ) &
        OffGridFlux_MHD_MF(:,iLevel) &
          = OffGridFlux_MHD_MF(:,iLevel) + OffGridFlux_MHD_X2_Outer

    END IF

    IF( .NOT. IsPeriodic(3) )THEN

      IF( iX_B0(3) .EQ. amrex_geom(iLevel) % domain % lo(3) ) &
        OffGridFlux_MHD_MF(:,iLevel) &
          = OffGridFlux_MHD_MF(:,iLevel) - OffGridFlux_MHD_X3_Inner
      IF( iX_E0(3) .EQ. amrex_geom(iLevel) % domain % hi(3) ) &
        OffGridFlux_MHD_MF(:,iLevel) &
          = OffGridFlux_MHD_MF(:,iLevel) + OffGridFlux_MHD_X3_Outer
    END IF

  END SUBROUTINE IncrementOffGridTally_MHD


END MODULE MF_MHD_dgDiscretizationModule
