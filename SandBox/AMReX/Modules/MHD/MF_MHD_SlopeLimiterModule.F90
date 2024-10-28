!> Module to apply slope-limiter to AMReX MultiFabs.
!> @todo Fix issue of multiple grids giving different results.
MODULE MF_MHD_SlopeLimiterModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    swX
  USE MeshModule, ONLY: &
    MeshX
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    nDM
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE MHD_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_MHD, &
    FinalizeSlopeLimiter_MHD, &
    ApplySlopeLimiter_MHD

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    AllocateArray_X, &
    DeallocateArray_X
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    DEBUG
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap
  USE AverageDownModule_MHD, ONLY: &
    AverageDown
  USE FillPatchModule_MHD, ONLY: &
    FillPatch

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_MHD_MF
  PUBLIC :: FinalizeSlopeLimiter_MHD_MF
  PUBLIC :: ApplySlopeLimiter_MHD_MF

  INTERFACE ApplySlopeLimiter_MHD_MF
    MODULE PROCEDURE ApplySlopeLimiter_MHD_MF_MultipleLevels
    MODULE PROCEDURE ApplySlopeLimiter_MHD_MF_SingleLevel
  END INTERFACE ApplySlopeLimiter_MHD_MF

  LOGICAL :: UseSlopeLimiter

CONTAINS


  SUBROUTINE InitializeSlopeLimiter_MHD_MF

    TYPE(amrex_parmparse) :: PP

    CHARACTER(:), ALLOCATABLE :: SlopeLimiterMethod
    REAL(DP)                  :: BetaTVD, BetaTVB
    REAL(DP)                  :: SlopeTolerance
    LOGICAL                   :: UseTroubledCellIndicator
    REAL(DP)                  :: LimiterThresholdParameter
    LOGICAL                   :: UseConservativeCorrection

    UseSlopeLimiter           = .TRUE.
    SlopeLimiterMethod        = 'TVD'
    BetaTVD                   = 1.75_DP
    BetaTVB                   = 0.00_DP
    SlopeTolerance            = 1.0e-6_DP
    UseTroubledCellIndicator  = .TRUE.
    LimiterThresholdParameter = 0.03_DP
    UseConservativeCorrection = .TRUE.
    CALL amrex_parmparse_build( PP, 'SL' )
      CALL PP % query( 'UseSlopeLimiter_MHD', &
                        UseSlopeLimiter )
      CALL PP % query( 'SlopeLimiterMethod_MHD', &
                        SlopeLimiterMethod )
      CALL PP % query( 'BetaTVD_MHD', &
                        BetaTVD )
      CALL PP % query( 'BetaTVB_MHD', &
                        BetaTVB )
      CALL PP % query( 'SlopeTolerance_MHD', &
                        SlopeTolerance )
      CALL PP % query( 'UseTroubledCellIndicator_MHD', &
                        UseTroubledCellIndicator )
      CALL PP % query( 'LimiterThresholdParameter_MHD', &
                        LimiterThresholdParameter )
      CALL PP % query( 'UseConservativeCorrection_MHD', &
                        UseConservativeCorrection )
    CALL amrex_parmparse_destroy( PP )

    CALL InitializeSlopeLimiter_MHD &
           ( BetaTVD_Option &
               = BetaTVD, &
             BetaTVB_Option &
               = BetaTVB, &
             SlopeTolerance_Option &
               = SlopeTolerance, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter, &
             SlopeLimiterMethod_Option &
               = SlopeLimiterMethod, &
             UseConservativeCorrection_Option &
               = UseConservativeCorrection, &
             Verbose_Option &
               = amrex_parallel_ioprocessor() )

  END SUBROUTINE InitializeSlopeLimiter_MHD_MF


  SUBROUTINE FinalizeSlopeLimiter_MHD_MF

    CALL FinalizeSlopeLimiter_MHD

  END SUBROUTINE FinalizeSlopeLimiter_MHD_MF


  SUBROUTINE ApplySlopeLimiter_MHD_MF_MultipleLevels &
    ( t, MF_uGF, MF_uCM, MF_uDM )

    REAL(DP)            , INTENT(in)    :: t(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM(0:)

    INTEGER :: iLevel, iErr

    DO iLevel = 0, nLevels-1

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(2x,A,I3.3)') &
            'CALL ApplySlopeLimiter_MHD_MF_SingleLevel, iLevel: ', &
            iLevel

        END IF

      END IF

      CALL ApplySlopeLimiter_MHD_MF_SingleLevel &
             ( t(iLevel), iLevel, MF_uGF, MF_uCM, MF_uDM )

    END DO

    CALL AverageDown( MF_uGF, MF_uCM )
    CALL AverageDown( MF_uGF, MF_uDM )

  END SUBROUTINE ApplySlopeLimiter_MHD_MF_MultipleLevels


  SUBROUTINE ApplySlopeLimiter_MHD_MF_SingleLevel &
    ( t, iLevel, MF_uGF, MF_uCM, MF_uDM )

    REAL(DP)            , INTENT(in)    :: t
    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM(0:)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDM (:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER       :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                     iLo_MF(4), iApplyBC(3)
    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    ! --- Apply boundary conditions to interior domains ---

    CALL FillPatch( iLevel, MF_uGF )
    CALL FillPatch( iLevel, MF_uGF, MF_uDM )
    CALL FillPatch( iLevel, MF_uGF, MF_uCM )

    CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( BX, MFI, uGF, uCF, uDM, G, U, D, &
    !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, iApplyBC, Edge_Map )
#endif

    CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF(iLevel) % DataPtr( MFI )
      uCM => MF_uCM(iLevel) % DataPtr( MFI )
      uDM => MF_uDM(iLevel) % DataPtr( MFI )

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

      CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

      CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCM, U )

      CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

      ! --- Apply boundary conditions to physical boundaries ---

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( t, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL Edge_Map % GetBC( iApplyBC )

      CALL ApplySlopeLimiter_MHD &
             ( t, iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
               SuppressBC_Option = .TRUE., iApplyBC_Option = iApplyBC )

      CALL thornado2amrex_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCM, U )

      CALL thornado2amrex_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

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

    END DO

    CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
    !$OMP END PARALLEL
#endif

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ApplySlopeLimiter_MHD_MF_SingleLevel


END MODULE MF_MHD_SlopeLimiterModule
