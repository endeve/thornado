MODULE MF_Euler_PositivityLimiterModule

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

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nDF
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D, &
    Max_D, &
    Min_T, &
    Max_T, &
    Min_Y, &
    Max_Y
  USE Euler_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_Euler, &
    FinalizePositivityLimiter_Euler, &
    ApplyPositivityLimiter_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    One
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    AllocateArray_X, &
    DeallocateArray_X
  USE InputParsingModule, ONLY: &
    UsePositivityLimiter_Euler, &
    Min_1_Euler, &
    Min_2_Euler, &
    D_Min_Euler_PL, &
    IntE_Min_Euler_PL, &
    EquationOfState, &
    nLevels, &
    UseTiling, &
    DEBUG
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap, &
    ApplyBoundaryConditions_Euler_MF
  USE FillPatchModule, ONLY: &
    FillPatch

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_Euler_MF
  PUBLIC :: FinalizePositivityLimiter_Euler_MF
  PUBLIC :: ApplyPositivityLimiter_Euler_MF

  INTERFACE ApplyPositivityLimiter_Euler_MF
    MODULE PROCEDURE ApplyPositivityLimiter_Euler_MF_MultipleLevels
    MODULE PROCEDURE ApplyPositivityLimiter_Euler_MF_SingleLevel
  END INTERFACE ApplyPositivityLimiter_Euler_MF

CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler_MF

    IF( TRIM( EquationOfState ) .EQ. 'TABLE' )THEN

      CALL InitializePositivityLimiter_Euler &
             ( UsePositivityLimiter_Option = UsePositivityLimiter_Euler, &
               Verbose_Option = amrex_parallel_ioprocessor(), &
               Min_1_Option = ( One + EPSILON(One) ) * Min_D, &
               Min_2_Option = ( One + EPSILON(One) ) * Min_T, &
               Min_3_Option = ( One + EPSILON(One) ) * Min_Y, &
               Max_1_Option = ( One - EPSILON(One) ) * Max_D, &
               Max_2_Option = ( One - EPSILON(One) ) * Max_T, &
               Max_3_Option = ( One - EPSILON(One) ) * Max_Y )

    ELSE

      CALL InitializePositivityLimiter_Euler &
             ( UsePositivityLimiter_Option = UsePositivityLimiter_Euler, &
               Verbose_Option = amrex_parallel_ioprocessor(), &
               Min_1_Option = Min_1_Euler, &
               Min_2_Option = Min_2_Euler, &
               D_Min_Euler_PL_Option    = D_Min_Euler_PL, &
               IntE_Min_Euler_PL_Option = IntE_Min_Euler_PL )

    END IF

  END SUBROUTINE InitializePositivityLimiter_Euler_MF


  SUBROUTINE FinalizePositivityLimiter_Euler_MF

    CALL FinalizePositivityLimiter_Euler

  END SUBROUTINE FinalizePositivityLimiter_Euler_MF


  SUBROUTINE ApplyPositivityLimiter_Euler_MF_MultipleLevels &
    ( MF_uGF, MF_uCF, MF_uDF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)

    INTEGER, PARAMETER :: nCycles = 1

    INTEGER :: iCycle, iLevel, iErr

    DO iCycle = 1, nCycles

      DO iLevel = 0, nLevels-1

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,'(2x,A,I3.3)') &
              'CALL ApplyPositivityLimiter_Euler_MF_SingleLevel, iLevel: ', &
              iLevel

          END IF

        END IF ! DEBUG

        CALL FillPatch( iLevel, 0.0_DP, MF_uGF, MF_uCF )

        CALL ApplyPositivityLimiter_Euler_MF_SingleLevel &
               ( iLevel, MF_uGF, MF_uCF, MF_uDF )

      END DO ! iLevel

      ! --- Ensure underlying coarse cells are consistent with
      !     cells on refined level ---

    END DO ! iCycle

  END SUBROUTINE ApplyPositivityLimiter_Euler_MF_MultipleLevels


  SUBROUTINE ApplyPositivityLimiter_Euler_MF_SingleLevel &
    ( iLevel, MF_uGF, MF_uCF, MF_uDF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER       :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                     iLo_MF(4)
    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter_Euler ) RETURN

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF(iLevel) % DataPtr( MFI )
      uCF => MF_uCF(iLevel) % DataPtr( MFI )
      uDF => MF_uDF(iLevel) % DataPtr( MFI )

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

      CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

      CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

      CALL amrex2thornado_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

      ! --- Apply boundary conditions to physical boundaries
      !     (needed for AMR) ---

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL ApplyPositivityLimiter_Euler( iX_B1, iX_E1, iX_B1, iX_E1, G, U, D )

      CALL thornado2amrex_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

      CALL thornado2amrex_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

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

    END DO

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ApplyPositivityLimiter_Euler_MF_SingleLevel


END MODULE MF_Euler_PositivityLimiterModule
