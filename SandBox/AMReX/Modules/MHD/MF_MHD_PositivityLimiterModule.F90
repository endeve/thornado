MODULE MF_MHD_PositivityLimiterModule

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
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    nDM
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE EquationOfStateModule, ONLY: &
    EquationOfState
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D, &
    Max_D, &
    Min_T, &
    Max_T, &
    Min_Y, &
    Max_Y
  USE MHD_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_MHD, &
    FinalizePositivityLimiter_MHD, &
    ApplyPositivityLimiter_MHD

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_UtilitiesModule_MHD, ONLY: &
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
  USE MF_EdgeMapModule_MHD, ONLY: &
    ConstructEdgeMap, &
    EdgeMap

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_MHD_MF
  PUBLIC :: FinalizePositivityLimiter_MHD_MF
  PUBLIC :: ApplyPositivityLimiter_MHD_MF

  INTERFACE ApplyPositivityLimiter_MHD_MF
    MODULE PROCEDURE ApplyPositivityLimiter_MHD_MF_MultipleLevels
    MODULE PROCEDURE ApplyPositivityLimiter_MHD_MF_SingleLevel
  END INTERFACE ApplyPositivityLimiter_MHD_MF

  LOGICAL  :: UsePositivityLimiter

CONTAINS


  SUBROUTINE InitializePositivityLimiter_MHD_MF &
    ( D_Min_MHD_PL_Option, IntE_Min_MHD_PL_Option )

    REAL(DP), INTENT(in), OPTIONAL :: D_Min_MHD_PL_Option
    REAL(DP), INTENT(in), OPTIONAL :: IntE_Min_MHD_PL_Option

    TYPE(amrex_parmparse) :: PP

    REAL(DP) :: Min_1, Min_2, D_Min_MHD_PL, IntE_Min_MHD_PL

    UsePositivityLimiter = .TRUE.
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % query( 'UsePositivityLimiter_MHD', &
                        UsePositivityLimiter )
    CALL amrex_parmparse_destroy( PP )

    IF( TRIM( EquationOfState ) .EQ. 'TABLE' )THEN

      CALL InitializePositivityLimiter_MHD &
             ( UsePositivityLimiter_Option = UsePositivityLimiter, &
               Verbose_Option = amrex_parallel_ioprocessor(), &
               Min_1_Option = ( One + EPSILON(One) ) * Min_D, &
               Min_2_Option = ( One + EPSILON(One) ) * Min_T, &
               Min_3_Option = ( One + EPSILON(One) ) * Min_Y, &
               Max_1_Option = ( One - EPSILON(One) ) * Max_D, &
               Max_2_Option = ( One - EPSILON(One) ) * Max_T, &
               Max_3_Option = ( One - EPSILON(One) ) * Max_Y )

    ELSE

      Min_1 = 1.0e-13_DP
      Min_2 = 1.0e-13_DP
      CALL amrex_parmparse_build( PP, 'PL' )
        CALL PP % query( 'Min_1_MHD', &
                          Min_1 )
        CALL PP % query( 'Min_2_MHD', &
                          Min_2 )
      CALL amrex_parmparse_destroy( PP )

      D_Min_MHD_PL = Zero
      IF( PRESENT( D_Min_MHD_PL_Option ) ) &
        D_Min_MHD_PL = D_Min_MHD_PL_Option

      IntE_Min_MHD_PL = Zero
      IF( PRESENT( IntE_Min_MHD_PL_Option ) ) &
        IntE_Min_MHD_PL = IntE_Min_MHD_PL_Option

      CALL InitializePositivityLimiter_MHD &
             ( UsePositivityLimiter_Option = UsePositivityLimiter, &
               Verbose_Option              = amrex_parallel_ioprocessor(), &
               Min_1_Option                = Min_1, &
               Min_2_Option                = Min_2, &
               D_Min_MHD_PL_Option       = D_Min_MHD_PL, &
               IntE_Min_MHD_PL_Option    = IntE_Min_MHD_PL )

    END IF

  END SUBROUTINE InitializePositivityLimiter_MHD_MF


  SUBROUTINE FinalizePositivityLimiter_MHD_MF

    CALL FinalizePositivityLimiter_MHD

  END SUBROUTINE FinalizePositivityLimiter_MHD_MF


  SUBROUTINE ApplyPositivityLimiter_MHD_MF_MultipleLevels &
    ( t, MF_uGF, MF_uCM, MF_uDM, swX_Option )

    REAL(DP)            , INTENT(in)    :: t(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swX_Option(3)

    INTEGER :: iLevel, iErr, swXX(3)

    swXX = 0
    IF( PRESENT( swX_Option ) ) &
      swXX = swX_Option

    DO iLevel = 0, nLevels-1

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(2x,A,I3.3)') &
            'CALL ApplyPositivityLimiter_MHD_MF_SingleLevel, iLevel: ', &
            iLevel

        END IF

      END IF ! DEBUG

      CALL ApplyPositivityLimiter_MHD_MF_SingleLevel &
             ( t(iLevel), iLevel, MF_uGF(iLevel), MF_uCM(iLevel), &
               MF_uDM(iLevel), &
               swX_Option = swXX )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyPositivityLimiter_MHD_MF_MultipleLevels


  SUBROUTINE ApplyPositivityLimiter_MHD_MF_SingleLevel &
    ( t, iLevel, MF_uGF, MF_uCM, MF_uDM, swX_Option )

    REAL(DP)            , INTENT(in)    :: t
    INTEGER             , INTENT(in)    :: iLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM
    INTEGER             , INTENT(in), OPTIONAL :: swX_Option(3)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDM (:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER       :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iX_B(3), iX_E(3), &
                     iLo_MF(4), swXX(3)
    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    swXX = 0
    IF( PRESENT( swX_Option ) ) &
      swXX = swX_Option

    CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( BX, MFI, uGF, uCM, uDM, G, U, D, &
    !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iX_B, iX_E, iLo_MF, Edge_Map )
#endif

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )
      uCM => MF_uCM % DataPtr( MFI )
      uDM => MF_uDM % DataPtr( MFI )

      iLo_MF = LBOUND( uGF )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX
      iX_B  = BX % lo - swXX
      iX_E  = BX % hi + swXX

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

      CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uGF, G )

      CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uCM, U )

      CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uDM, D )

      ! --- Apply boundary conditions to physical boundaries
      !     (needed for AMR) ---

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_MHD_MF &
             ( t, iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL ApplyPositivityLimiter_MHD &
             ( iX_B , iX_E , iX_B1, iX_E1, G, U, D )

      CALL thornado2amrex_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uCM, U )

      CALL thornado2amrex_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uDM, D )

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

  END SUBROUTINE ApplyPositivityLimiter_MHD_MF_SingleLevel


END MODULE MF_MHD_PositivityLimiterModule
