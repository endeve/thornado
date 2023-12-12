MODULE MF_XCFC_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFZ, &
    iE_E0, &
    iE_B0, &
    nDOFE, &
    nE, &
    swX
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Psi, &
    nGF
  USE GeometryFieldsModuleE, ONLY: &
    iGE_Ep3, &
    uGE
  USE GeometryBoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting, &
    ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nCF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nPF
  USE RadiationFieldsModule, ONLY: &
    iCR_N, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nCR, &
    nSp => nSpecies
  USE XCFC_UtilitiesModule, ONLY: &
    nMF, &
    iGS_E, &
    iGS_S1, &
    iGS_S2, &
    iGS_S3, &
    iGS_S, &
    nGS, &
    swX_GS, &
    MultiplyWithPsi6, &
    ComputeGravitationalMass, &
    UpdateConformalFactorAndMetric_XCFC, &
    UpdateLapseShiftCurvature_XCFC
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE Euler_XCFC_UtilitiesModule, ONLY: &
    ComputeConformalFactorSourcesAndMg_XCFC_Euler, &
    ComputePressureTensorTrace_XCFC_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    FourPi
  USE MF_UtilitiesModule, ONLY: &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z, &
    amrex2thornado_X, &
    thornado2amrex_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MultiplyWithPsi6_MF
  PUBLIC :: ApplyBoundaryConditions_Geometry_XCFC_MF
  PUBLIC :: ComputeGravitationalMass_MF
  PUBLIC :: UpdateConformalFactorAndMetric_XCFC_MF
  PUBLIC :: UpdateLapseShiftCurvature_XCFC_MF
  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_MF
  PUBLIC :: ComputePressureTensorTrace_XCFC_MF

  INTERFACE MultiplyWithPsi6_MF
    MODULE PROCEDURE MultiplyWithPsi6_MF_X
    MODULE PROCEDURE MultiplyWithPsi6_MF_Z
  END INTERFACE MultiplyWithPsi6_MF

CONTAINS


  SUBROUTINE MultiplyWithPsi6_MF_X &
    ( MF_uGF, MF_uCF, Power, OnlyLeafElements_Option )

    INTEGER             , INTENT(in)    :: Power
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    LOGICAL             , INTENT(in), OPTIONAL :: OnlyLeafElements_Option

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    LOGICAL :: OnlyLeafElements

    OnlyLeafElements = .TRUE.
    IF( PRESENT( OnlyLeafElements_Option ) ) &
      OnlyLeafElements = OnlyLeafElements_Option

    DO iLevel = 0, nLevels-1

      IF( OnlyLeafElements ) &
        CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G, U, uGF, uCF, uFM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        IF( OnlyLeafElements ) &
          uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

        IF( OnlyLeafElements )THEN

          CALL MultiplyWithPsi6 &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Power, Mask_Option = uFM )

        ELSE

          CALL MultiplyWithPsi6 &
                 ( iX_B1, iX_E1, iX_B1, iX_E1, G, U, Power )

        END IF

        CALL thornado2amrex_X &
               ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE MultiplyWithPsi6_MF_X


  SUBROUTINE MultiplyWithPsi6_MF_Z &
    ( iE_B, iE_E, iE_B1, iE_E1, MF_uGF, MF_uCR, Power )

    INTEGER             , INTENT(in)    :: iE_B, iE_E, iE_B1, iE_E1, Power
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER  :: iLevel
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), nE
    INTEGER  :: iZ_B1(4), iZ_E1(4)

    nE = iE_E - iE_B + 1

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G, U, uGF, uCR, uFM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iZ_B1, iZ_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        iZ_B1(1)   = iE_B1
        iZ_B1(2:4) = iX_B1
        iZ_E1(1)   = iE_E1
        iZ_E1(2:4) = iX_E1

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_Z &
               ( [ 1    , iZ_B1(1), iZ_B1(2), iZ_B1(3), iZ_B1(4), 1  , 1   ], &
                 [ nDOFZ, iZ_E1(1), iZ_E1(2), iZ_E1(3), iZ_E1(4), nCR, nSp ], &
                 U )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_Z &
               ( nCR, nSp, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, LBOUND( uCR ), iZ_B1, iZ_E1, uCR, U )

        CALL MultiplyWithPsi6 &
               ( iE_B1, iE_E1, iX_B0, iX_E0, iX_B1, iX_E1, &
                 G, U, Power, Mask_Option = uFM )

        CALL thornado2amrex_Z &
               ( nCR, nSp, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, LBOUND( uCR ), iZ_B1, iZ_E1, uCR, U )

        CALL DeallocateArray_Z &
               ( [ 1    , iZ_B1(1), iZ_B1(2), iZ_B1(3), iZ_B1(4), 1  , 1   ], &
                 [ nDOFZ, iZ_E1(1), iZ_E1(2), iZ_E1(3), iZ_E1(4), nCR, nSp ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE MultiplyWithPsi6_MF_Z


  SUBROUTINE UpdateConformalFactorAndMetric_XCFC_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:) ! Metric Fields
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: M(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uMF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX_B(3), iX_E(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, M, G, uMF, uGF, uFM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iX_B, iX_E )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX
        iX_B  = iX_B0 - swX_GS
        iX_E  = iX_E0 + swX_GS

        CALL AllocateArray_X &
               ( [ 1    , iX_B(1), iX_B(2), iX_B(3), 1   ], &
                 [ nDOFX, iX_E(1), iX_E(2), iX_E(3), nMF ], &
                 M )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL amrex2thornado_X &
               ( nMF, iX_B, iX_E, LBOUND( uMF ), iX_B, iX_E, uMF, M )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL UpdateConformalFactorAndMetric_XCFC &
               ( iX_B0, iX_E0, iX_B1, iX_E1, M, G, Mask_Option = uFM )

        CALL thornado2amrex_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B, iX_E, uGF, G )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B(1), iX_B(2), iX_B(3), 1   ], &
                 [ nDOFX, iX_E(1), iX_E(2), iX_E(3), nMF ], &
                 M )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE UpdateConformalFactorAndMetric_XCFC_MF


  SUBROUTINE UpdateLapseShiftCurvature_XCFC_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: M(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uMF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX_B(3), iX_E(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uMF, uGF, uFM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iX_B, iX_E )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX
        iX_B  = iX_B0 - swX_GS
        iX_E  = iX_E0 + swX_GS

        CALL AllocateArray_X &
               ( [ 1    , iX_B(1), iX_B(2), iX_B(3), 1   ], &
                 [ nDOFX, iX_E(1), iX_E(2), iX_E(3), nMF ], &
                 M )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL amrex2thornado_X &
               ( nMF, iX_B, iX_E, LBOUND( uMF ), iX_B, iX_E, uMF, M )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL UpdateLapseShiftCurvature_XCFC &
               ( iX_B0, iX_E0, iX_B1, iX_E1, M, G, Mask_Option = uFM )

        CALL thornado2amrex_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B, iX_E, uGF, G )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B(1), iX_B(2), iX_B(3), 1   ], &
                 [ nDOFX, iX_E(1), iX_E(2), iX_E(3), nMF ], &
                 M )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE UpdateLapseShiftCurvature_XCFC_MF


  SUBROUTINE ComputeGravitationalMass_MF( MF_uGS, GravitationalMass )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:)
    REAL(DP)            , INTENT(out) :: GravitationalMass

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: GS(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLevel

    ! --- Assuming 1D spherical symmetry ---

    GravitationalMass = Zero

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGS % BA, MF_uGS % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, GS, uGS, uFM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1 ) &
      !$OMP REDUCTION( +:GravitationalMass )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGS(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGS => MF_uGS(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGS ], &
                 GS )

        CALL amrex2thornado_X &
               ( nGS, iX_B0, iX_E0, LBOUND( uGS ), iX_B0, iX_E0, uGS, GS )

        CALL ComputeGravitationalMass &
               ( iX_B0, iX_E0, iX_B1, iX_E1, GS, GravitationalMass, &
                 Mask_Option = uFM )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGS ], &
                 GS )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

    CALL amrex_parallel_reduce_sum( GravitationalMass )

  END SUBROUTINE ComputeGravitationalMass_MF


  SUBROUTINE ApplyBoundaryConditions_Geometry_XCFC_MF( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    CALL ApplyBoundaryConditions_X1_Inner( MF_uGF )
    CALL ApplyBoundaryConditions_X1_Outer( MF_uGF )

  END SUBROUTINE ApplyBoundaryConditions_Geometry_XCFC_MF


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

#else

  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_MF &
    ( MF_uGF, MF_uCF, MF_uGS )

#endif

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: F (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: GS(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
#ifndef THORNADO_NOTRANSPORT
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
#endif
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER  :: iLevel
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

! #if defined( THORNADO_OMP )
!       !$OMP PARALLEL &
!       !$OMP PRIVATE( BX, MFI, G, F, S, uGF, uCF, uGS, uFM, &
!       !$OMP          iX_B0, iX_E0, iX_B1, iX_E1 )
! #endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 F )

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGS ], &
                 GS )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, F )

        CALL amrex2thornado_X &
               ( nGS, iX_B0, iX_E0, LBOUND( uGS ), iX_B0, iX_E0, uGS, GS )

        CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, F, GS, Mask_Option = uFM )

        CALL thornado2amrex_X &
               ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, F )

        CALL thornado2amrex_X &
               ( nGS, iX_B0, iX_E0, LBOUND( uGS ), iX_B0, iX_E0, uGS, GS )

#ifndef THORNADO_NOTRANSPORT

        uCR => MF_uCR(iLevel) % DataPtr( MFI )

        CALL ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment &
               ( iLevel, iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uCR, uGS )

#endif

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGS ], &
                 GS )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 F )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

! #if defined( THORNADO_OMP )
!       !$OMP END PARALLEL
! #endif

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_MF


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE ComputePressureTensorTrace_XCFC_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

#else

  SUBROUTINE ComputePressureTensorTrace_XCFC_MF &
    ( MF_uGF, MF_uCF, MF_uGS )

#endif

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: GS(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
#ifndef THORNADO_NOTRANSPORT
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
#endif
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER  :: iLevel
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

! #if defined( THORNADO_OMP )
!       !$OMP PARALLEL &
!       !$OMP PRIVATE( BX, MFI, G, U, S, uGF, uCF, uGS, uFM, &
!       !$OMP          iX_B0, iX_E0, iX_B1, iX_E1 )
! #endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGS ], &
                 GS )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

        CALL amrex2thornado_X &
               ( nGS, iX_B0, iX_E0, LBOUND( uGS ), iX_B0, iX_E0, uGS, GS )

        CALL ComputePressureTensorTrace_XCFC_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, GS, Mask_Option = uFM )

        CALL thornado2amrex_X &
               ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

        CALL thornado2amrex_X &
               ( nGS, iX_B0, iX_E0, LBOUND( uGS ), iX_B0, iX_E0, uGS, GS )

#ifndef THORNADO_NOTRANSPORT

        uCR => MF_uCR(iLevel) % DataPtr( MFI )

        CALL ComputePressureTensorTrace_XCFC_TwoMoment &
               ( iLevel, iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uCR, uGS )

#endif

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGS ], &
                 GS )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

! #if defined( THORNADO_OMP )
!       !$OMP END PARALLEL
! #endif

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ComputePressureTensorTrace_XCFC_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE ApplyBoundaryConditions_X1_Inner( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    ! --- Inner Boundary: Reflecting ---

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

          CALL AllocateArray_X &
                 ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                   [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                   G )

          CALL amrex2thornado_X &
                 ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

          CALL ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G )

          CALL thornado2amrex_X &
                 ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

          CALL DeallocateArray_X &
                 ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                   [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                   G )

        END IF

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyBoundaryConditions_X1_Inner


  SUBROUTINE ApplyBoundaryConditions_X1_Outer( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    ! --- Outer Boundary: Extrapolate fields to face ---

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) )THEN

          CALL AllocateArray_X &
                 ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                   [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                   G )

          CALL amrex2thornado_X &
                 ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

          CALL ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G )

          CALL thornado2amrex_X &
                 ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

          CALL DeallocateArray_X &
                 ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                   [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                   G )

        END IF

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyBoundaryConditions_X1_Outer


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment &
    ( iLevel, iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uCR, uGS )

    INTEGER , INTENT(in)    :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: uCF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: uCR(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: uGS(iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iCF, ErrorExists
    INTEGER  :: iS, iE, iN_E, iN_Z
    INTEGER  :: iD_N, iD_G1, iD_G2, iD_G3

    REAL(DP) :: Psi6
    REAL(DP) :: uPF(nPF), LorentzFactor

    REAL(DP) :: E, S_i(3), E_int, S_i_int(3)
    REAL(DP) :: N, G_d_1, G_d_2, G_d_3, VdotG
    REAL(DP) :: V_u_1, V_u_2, V_u_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3

    INTEGER :: ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))

    ASSOCIATE( dZ1 => MeshE  % Width )

    ErrorExists = 0

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      ITERATION(iNX,iX1,iX2,iX3) = 0
      iErr     (iNX,iX1,iX2,iX3) = 0

      ! Assume Psi^(iStage) ~ Psi^(iStage+1) for Poseidon BCs

      Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

      DO iCF = 1, nCF

        uCF    (iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) &
          = uCF(iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) / Psi6

      END DO

      CALL ComputePrimitive_Euler &
             ( uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX), &
               uPF(iPF_D ), &
               uPF(iPF_V1), &
               uPF(iPF_V2), &
               uPF(iPF_V3), &
               uPF(iPF_E ), &
               uPF(iPF_Ne), &
               uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
               uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
               uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
               ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
               iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

      DO iCF = 1, nCF

        uCF    (iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) &
          = uCF(iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) * Psi6

      END DO

      ErrorExists = ErrorExists + iErr(iNX,iX1,iX2,iX3)

      V_u_1 = uPF(iPF_V1)
      V_u_2 = uPF(iPF_V2)
      V_u_3 = uPF(iPF_V3)

      V_d_1 = V_u_1 * uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX)
      V_d_2 = V_u_2 * uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX)
      V_d_3 = V_u_3 * uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX)

      LorentzFactor &
        = One / SQRT( One - ( V_u_1 * V_d_1 + V_u_2 * V_d_2 + V_u_3 * V_d_3 ) )

      E       = Zero
      E_int   = Zero
      S_i     = Zero
      S_i_int = Zero

      DO iS   = 1, nSp
      DO iE   = 1, nE
      DO iN_E = 1, nDOFE

        iN_Z = ( iNX - 1 ) * nDOFE + iN_E

        iD_N  = MapIndex( iS, iE, iN_Z, iCR_N  )
        iD_G1 = MapIndex( iS, iE, iN_Z, iCR_G1 )
        iD_G2 = MapIndex( iS, iE, iN_Z, iCR_G2 )
        iD_G3 = MapIndex( iS, iE, iN_Z, iCR_G3 )

        N     = uCR(iX1,iX2,iX3,iD_N )
        G_d_1 = uCR(iX1,iX2,iX3,iD_G1)
        G_d_2 = uCR(iX1,iX2,iX3,iD_G2)
        G_d_3 = uCR(iX1,iX2,iX3,iD_G3)

        VdotG = V_u_1 * G_d_1 + V_u_2 * G_d_2 + V_u_3 * G_d_3

        E_int      = LorentzFactor * N + VdotG
        S_i_int(1) = LorentzFactor * V_d_1 * N + G_d_1
        S_i_int(2) = LorentzFactor * V_d_2 * N + G_d_2
        S_i_int(3) = LorentzFactor * V_d_3 * N + G_d_3

        E &
          = E &
              + FourPi * dZ1(iE) * WeightsE(iN_E) * uGE(iN_E,iE,iGE_Ep3) &
              * E_int

        S_i &
          = S_i &
              + FourPi * dZ1(iE) * WeightsE(iN_E) * uGE(iN_E,iE,iGE_Ep3) &
              * S_i_int

      END DO
      END DO
      END DO

      uGS(iX1,iX2,iX3,nDOFX*(iGS_E -1)+iNX) &
        = uGS(iX1,iX2,iX3,nDOFX*(iGS_E -1)+iNX) + E

      uGS(iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) &
        = uGS(iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) + S_i(1)

      uGS(iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) &
        = uGS(iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) + S_i(2)

      uGS(iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) &
        = uGS(iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) + S_i(3)

    END DO
    END DO
    END DO
    END DO

    IF( ErrorExists .GT. 0 )THEN

      WRITE(*,*) 'ERROR'
      WRITE(*,*) '-----'
      WRITE(*,*) &
        '    MODULE: MF_TwoMoment_XCFC_UtilitiesModule'
      WRITE(*,*) &
        'SUBROUTINE: ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment'

      CALL CreateMesh_MF( iLevel, MeshX )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

        CALL DescribeError_Euler &
               ( iErr(iNX,iX1,iX2,iX3), &
                 Int_Option &
                   = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                       iX_B0(1), iX_B0(2), iX_B0(3), &
                       iX_E0(1), iX_E0(2), iX_E0(3), &
                       iNX, iX1, iX2, iX3 ], &
                 Real_Option &
                   = [ MeshX(1) % Center(iX1), &
                       MeshX(2) % Center(iX2), &
                       MeshX(3) % Center(iX3), &
                       MeshX(1) % Width (iX1), &
                       MeshX(2) % Width (iX2), &
                       MeshX(3) % Width (iX3), &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ], &
                 Char_Option = [ 'NA' ], &
                 Message_Option &
                   = 'Calling from ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment' )

      END DO
      END DO
      END DO
      END DO

      CALL DestroyMesh_MF( MeshX )

    END IF

    END ASSOCIATE

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment


  SUBROUTINE ComputePressureTensorTrace_XCFC_TwoMoment &
    ( iLevel, iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uCR, uGS )

    INTEGER , INTENT(in)    :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: uGF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: uCF(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: uCR(iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: uGS(iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iCF, ErrorExists
    INTEGER  :: iS, iE, iN_E, iN_Z
    INTEGER  :: iD_N, iD_G1, iD_G2, iD_G3

    REAL(DP) :: Psi6
    REAL(DP) :: uPF(nPF), LorentzFactor

    REAL(DP) :: S, S_int
    REAL(DP) :: N, G_d_1, G_d_2, G_d_3, VdotG
    REAL(DP) :: V_u_1, V_u_2, V_u_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3

    INTEGER :: ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))
    INTEGER :: iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                 iX_B0(2):iX_E0(2), &
                                 iX_B0(3):iX_E0(3))

    ASSOCIATE( dZ1 => MeshE  % Width )

    ErrorExists = 0

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      ITERATION(iNX,iX1,iX2,iX3) = 0
      iErr     (iNX,iX1,iX2,iX3) = 0

      Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

      ! --- Compute trace of stress tensor ---

      DO iCF = 1, nCF

        uCF    (iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) &
          = uCF(iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) / Psi6

      END DO

      CALL ComputePrimitive_Euler &
             ( uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX), &
               uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX), &
               uPF(iPF_D ), &
               uPF(iPF_V1), &
               uPF(iPF_V2), &
               uPF(iPF_V3), &
               uPF(iPF_E ), &
               uPF(iPF_Ne), &
               uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
               uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
               uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
               ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
               iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

      DO iCF = 1, nCF

        uCF    (iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) &
          = uCF(iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) * Psi6

      END DO

      ErrorExists = ErrorExists + iErr(iNX,iX1,iX2,iX3)

      V_u_1 = uPF(iPF_V1)
      V_u_2 = uPF(iPF_V2)
      V_u_3 = uPF(iPF_V3)

      V_d_1 = V_u_1 * uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX)
      V_d_2 = V_u_2 * uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX)
      V_d_3 = V_u_3 * uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX)

      LorentzFactor &
        = One / SQRT( One - ( V_u_1 * V_d_1 + V_u_2 * V_d_2 + V_u_3 * V_d_3 ) )

      S = Zero

      DO iS   = 1, nSp
      DO iE   = 1, nE
      DO iN_E = 1, nDOFE

        iN_Z = ( iNX - 1 ) * nDOFE + iN_E

        iD_N  = MapIndex( iS, iE, iN_Z, iCR_N  )
        iD_G1 = MapIndex( iS, iE, iN_Z, iCR_G1 )
        iD_G2 = MapIndex( iS, iE, iN_Z, iCR_G2 )
        iD_G3 = MapIndex( iS, iE, iN_Z, iCR_G3 )

        N     = uCR(iX1,iX2,iX3,iD_N )
        G_d_1 = uCR(iX1,iX2,iX3,iD_G1)
        G_d_2 = uCR(iX1,iX2,iX3,iD_G2)
        G_d_3 = uCR(iX1,iX2,iX3,iD_G3)

        VdotG = V_u_1 * G_d_1 + V_u_2 * G_d_2 + V_u_3 * G_d_3

        S_int = LorentzFactor * N + VdotG

        S &
          = S + FourPi * dZ1(iE) * WeightsE(iN_E) * uGE(iN_E,iE,iGE_Ep3) &
            * S_int

      END DO
      END DO
      END DO

      uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) &
        = uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) + S

    END DO
    END DO
    END DO
    END DO

    IF( ErrorExists .NE. 0 )THEN

      WRITE(*,*) &
        'ERROR'
      WRITE(*,*) &
        '-----'
      WRITE(*,*) &
        '    MODULE: Poseidon_UtilitiesModule'
      WRITE(*,*) &
        'SUBROUTINE: ComputePressureTensorTrace_XCFC_TwoMoment'

      CALL CreateMesh_MF( iLevel, MeshX )

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1       , nDOFX

        Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

        CALL DescribeError_Euler &
               ( iErr(iNX,iX1,iX2,iX3), &
                 Int_Option &
                   = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                       iX_B0(1), iX_B0(2), iX_B0(3), &
                       iX_E0(1), iX_E0(2), iX_E0(3), &
                       iNX, iX1, iX2, iX3 ], &
                 Real_Option &
                   = [ MeshX(1) % Center(iX1), &
                       MeshX(2) % Center(iX2), &
                       MeshX(3) % Center(iX3), &
                       MeshX(1) % Width (iX1), &
                       MeshX(2) % Width (iX2), &
                       MeshX(3) % Width (iX3), &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                       uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                       uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ], &
                 Char_Option = [ 'NA' ], &
                 Message_Option &
                   = 'Calling from ComputePressureTensorTrace_XCFC_TwoMoment' )

      END DO
      END DO
      END DO
      END DO

      CALL DestroyMesh_MF( MeshX )

    END IF

    END ASSOCIATE

  END SUBROUTINE ComputePressureTensorTrace_XCFC_TwoMoment


  INTEGER FUNCTION MapIndex( iS, iE, iN_Z, iCR ) RESULT( iD )

    INTEGER, INTENT(in) :: iS, iE, iN_Z, iCR

    iD = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
           + ( iCR - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
           + ( iE - 1 ) * nDOFZ + iN_Z

    RETURN
  END FUNCTION MapIndex


END MODULE MF_XCFC_UtilitiesModule
