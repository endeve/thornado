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
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

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
    Zero, &
    One
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    AllocateArray_X, &
    DeallocateArray_X
  USE InputParsingModule, ONLY: &
    EquationOfState, &
    nLevels, &
    UseTiling, &
    DEBUG
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE MF_EdgeMapModule, ONLY: &
    ConstructEdgeMap, &
    EdgeMap

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_Euler_MF
  PUBLIC :: FinalizePositivityLimiter_Euler_MF
  PUBLIC :: ApplyPositivityLimiter_Euler_MF

  INTERFACE ApplyPositivityLimiter_Euler_MF
    MODULE PROCEDURE ApplyPositivityLimiter_Euler_MF_MultipleLevels
    MODULE PROCEDURE ApplyPositivityLimiter_Euler_MF_SingleLevel
  END INTERFACE ApplyPositivityLimiter_Euler_MF

  LOGICAL  :: UsePositivityLimiter

CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler_MF &
    ( D_Min_Euler_PL_Option, IntE_Min_Euler_PL_Option )

    REAL(DP), INTENT(in), OPTIONAL :: D_Min_Euler_PL_Option
    REAL(DP), INTENT(in), OPTIONAL :: IntE_Min_Euler_PL_Option

    TYPE(amrex_parmparse) :: PP

    REAL(DP) :: Min_1, Min_2, D_Min_Euler_PL, IntE_Min_Euler_PL

    UsePositivityLimiter = .TRUE.
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % query( 'UsePositivityLimiter_Euler', &
                        UsePositivityLimiter )
    CALL amrex_parmparse_destroy( PP )

    IF( TRIM( EquationOfState ) .EQ. 'TABLE' )THEN

      CALL InitializePositivityLimiter_Euler &
             ( UsePositivityLimiter_Option = UsePositivityLimiter, &
               Verbose_Option = amrex_parallel_ioprocessor(), &
               Min_1_Option = ( One + EPSILON(One) ) * Min_D, &
               Min_2_Option = ( One + EPSILON(One) ) * Min_T, &
               Min_3_Option = ( One + EPSILON(One) ) * Min_Y, &
               Max_1_Option = ( One - EPSILON(One) ) * Max_D, &
               Max_2_Option = ( One - EPSILON(One) ) * Max_T, &
               Max_3_Option = ( One - EPSILON(One) ) * Max_Y )

    ELSE

      Min_1 = 1.0e-12_DP
      Min_2 = 1.0e-12_DP
      CALL amrex_parmparse_build( PP, 'PL' )
        CALL PP % query( 'Min_1_Euler', &
                          Min_1 )
        CALL PP % query( 'Min_2_Euler', &
                          Min_2 )
      CALL amrex_parmparse_destroy( PP )

      D_Min_Euler_PL = Zero
      IF( PRESENT( D_Min_Euler_PL_Option ) ) &
        D_Min_Euler_PL = D_Min_Euler_PL_Option

      IntE_Min_Euler_PL = Zero
      IF( PRESENT( IntE_Min_Euler_PL_Option ) ) &
        IntE_Min_Euler_PL = IntE_Min_Euler_PL_Option

      CALL InitializePositivityLimiter_Euler &
             ( UsePositivityLimiter_Option = UsePositivityLimiter, &
               Verbose_Option              = amrex_parallel_ioprocessor(), &
               Min_1_Option                = Min_1, &
               Min_2_Option                = Min_2, &
               D_Min_Euler_PL_Option       = D_Min_Euler_PL, &
               IntE_Min_Euler_PL_Option    = IntE_Min_Euler_PL )

    END IF

  END SUBROUTINE InitializePositivityLimiter_Euler_MF


  SUBROUTINE FinalizePositivityLimiter_Euler_MF

    CALL FinalizePositivityLimiter_Euler

  END SUBROUTINE FinalizePositivityLimiter_Euler_MF


  SUBROUTINE ApplyPositivityLimiter_Euler_MF_MultipleLevels &
    ( MF_uGF, MF_uCF, MF_uDF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)

    INTEGER :: iLevel, iErr

    DO iLevel = 0, nLevels-1

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(2x,A,I3.3)') &
            'CALL ApplyPositivityLimiter_Euler_MF_SingleLevel, iLevel: ', &
            iLevel

        END IF

      END IF ! DEBUG

      CALL ApplyPositivityLimiter_Euler_MF_SingleLevel &
             ( iLevel, MF_uGF(iLevel), MF_uCF(iLevel), MF_uDF(iLevel) )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyPositivityLimiter_Euler_MF_MultipleLevels


  SUBROUTINE ApplyPositivityLimiter_Euler_MF_SingleLevel &
    ( iLevel, MF_uGF, MF_uCF, MF_uDF )

    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDF (:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER       :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                     iLo_MF(4)
    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( BX, MFI, uGF, uCF, uDF, G, U, D, &
    !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, Edge_Map )
#endif

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCF => MF_uCF % DataPtr( MFI )
      uDF => MF_uDF % DataPtr( MFI )

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

      CALL ApplyPositivityLimiter_Euler &
             ( iX_B1, iX_E1, iX_B1, iX_E1, G, U, D )

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

#if defined( THORNADO_OMP )
    !$OMP END PARALLEL
#endif

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ApplyPositivityLimiter_Euler_MF_SingleLevel


END MODULE MF_Euler_PositivityLimiterModule
