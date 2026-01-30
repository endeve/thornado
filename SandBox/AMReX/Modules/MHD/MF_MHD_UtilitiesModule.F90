MODULE MF_MHD_UtilitiesModule

  ! --- AMReX Modules ---

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
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    swX
  USE MeshModule, ONLY: &
    MeshX
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    Min_T
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    nGF
  USE MagnetofluidFieldsModule, ONLY: &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_Ne, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3, &
    iCM_Chi, &
    nCM, &
    iPM_D, &
    iPM_V1, &
    iPM_V2, &
    iPM_V3, &
    iPM_E, &
    iPM_Ne, &
    iPM_B1, &
    iPM_B2, &
    iPM_B3, &
    iPM_Chi, &
    nPM, &
    iAM_P, &
    iAM_Ye, &
    nAM, &
    iDM_MinE, &
    nDM
  USE MHD_UtilitiesModule, ONLY: &
    ComputeTimeStep_MHD, &
    ComputeFromConserved_MHD, &
    ComputeConserved_MHD

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    One
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    EvolveOnlyMagnetic, &
    UseDivergenceCleaning, &
    CleaningSpeed, &
    DampingTimeScaleFactor, &
    DEBUG
  USE MF_UtilitiesModule_MHD, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    AllocateArray_X, &
    DeallocateArray_X
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_TimersModule_MHD, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_ComputeTimeStep_MHD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeFromConserved_MHD_MF
  PUBLIC :: ComputeTimeStep_MHD_MF
  PUBLIC :: ComputeConserved_MHD_MF
  PUBLIC :: ComputeDiagnosticFields_MHD_MF

CONTAINS


  SUBROUTINE ComputeFromConserved_MHD_MF &
    ( MF_uGF, MF_uCM, MF_uPM, MF_uAM, &
      swXX_Option )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:), MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uPM(0:), MF_uAM(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAM(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: A(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iX_B(3), iX_E(3), swXX(3)

    swXX = 0
    IF( PRESENT( swXX_Option ) ) &
      swXX = swXX_Option

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( MFI, BX, uGF, uCM, uPM, uAM, G, U, P, A, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, iX_B, iX_E )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )
        uPM => MF_uPM(iLevel) % DataPtr( MFI )
        uAM => MF_uAM(iLevel) % DataPtr( MFI )

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
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPM ], &
                 P )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAM ], &
                 A )

        iX_B = iX_B0 - swXX
        iX_E = iX_E0 + swXX

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uGF, G )

        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uCM, U )

        CALL ComputeFromConserved_MHD &
               ( iX_B, iX_E, iX_B1, iX_E1, G, U, P, A, EvolveOnlyMagnetic )

        CALL thornado2amrex_X( nPM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uPM, P )

        CALL thornado2amrex_X( nAM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uAM, A )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAM ], &
                 A )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPM ], &
                 P )

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

    END DO

  END SUBROUTINE ComputeFromConserved_MHD_MF


  SUBROUTINE ComputeTimeStep_MHD_MF &
    ( MF_uGF, MF_uCM, CFL, TimeStepMin )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:), MF_uCM(0:)
    REAL(DP)            , INTENT(in)  :: CFL
    REAL(DP)            , INTENT(out) :: TimeStepMin(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4), iErr

    REAL(DP) :: TimeStep(0:nLevels-1)

    CALL TimersStart_AMReX( Timer_AMReX_ComputeTimeStep_MHD )

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() ) &
        WRITE(*,'(A)') 'CALL ComputeTimeStep_MHD_MF'

    END IF

    TimeStepMin = HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( MFI, BX, uGF, uCM, G, U, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, TimeStep ) &
      !$OMP REDUCTION( MIN:TimeStepMin )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

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
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
                 U )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCM, U )

        CALL ComputeTimeStep_MHD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep(iLevel), &
                 UseDivergenceCleaning, CleaningSpeed, DampingTimeScaleFactor, EvolveOnlyMagnetic )

        TimeStepMin(iLevel) = MIN( TimeStepMin(iLevel), TimeStep(iLevel) )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCM ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyMesh_MF( MeshX )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStepMin, SIZE( TimeStepMin ) )

    CALL TimersStop_AMReX( Timer_AMReX_ComputeTimeStep_MHD )

  END SUBROUTINE ComputeTimeStep_MHD_MF


  SUBROUTINE ComputeConserved_MHD_MF( MF_uGF, MF_uPM, MF_uAM, MF_uCM )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPM(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iNX, iX1, iX2, iX3

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( MFI, BX, uGF, uPM, uAM, uCM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uPM => MF_uPM(iLevel) % DataPtr( MFI )
        uAM => MF_uAM(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          CALL ComputeConserved_MHD &
                 ( uPM(iX1,iX2,iX3,nDOFX*(iPM_D       -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_V1      -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_V2      -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_V3      -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_E       -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_Ne      -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_B1      -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_B2      -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_B3      -1)+iNX), &
                   uPM(iX1,iX2,iX3,nDOFX*(iPM_Chi     -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_D       -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_S1      -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_S2      -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_S3      -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_E       -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_Ne      -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_B1      -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_B2      -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_B3      -1)+iNX), &
                   uCM(iX1,iX2,iX3,nDOFX*(iCM_Chi     -1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha   -1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1  -1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2  -1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3  -1)+iNX), &
                   uAM(iX1,iX2,iX3,nDOFX*(iAM_P       -1)+iNX), &
                   EvolveOnlyMagnetic )

        END DO
        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO

  END SUBROUTINE ComputeConserved_MHD_MF


  SUBROUTINE ComputeDiagnosticFields_MHD_MF &
    ( MF_uGF, MF_uCM, MF_uDM, swXX_Option )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:), MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    TYPE(amrex_multifab) :: MF_uPM(0:nLevels-1)
    TYPE(amrex_multifab) :: MF_uAM(0:nLevels-1)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDM(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: A(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iX_B(3), iX_E(3), swXX(3)

    swXX = 0
    IF( PRESENT( swXX_Option ) ) &
      swXX = swXX_Option

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uPM(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nPM, swX )

      CALL amrex_multifab_build &
             ( MF_uAM(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nAM, swX )

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( MFI, BX, uGF, uCM, uPM, uAM, uDM, G, U, P, A, D, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, iX_B, iX_E )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )
        uPM => MF_uPM(iLevel) % DataPtr( MFI )
        uAM => MF_uAM(iLevel) % DataPtr( MFI )
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
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPM ], &
                 P )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAM ], &
                 A )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDM ], &
                 D )

        iX_B = iX_B0 - swXX
        iX_E = iX_E0 + swXX

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uGF, G )

        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uCM, U )

        CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uDM, D )

        CALL ComputeFromConserved_MHD &
               ( iX_B, iX_E, iX_B1, iX_E1, G, U, P, A, EvolveOnlyMagnetic )

        CALL ComputeDiagnosticFields_MHD &
               ( iX_B, iX_E, iX_B1, iX_E1, P, A, D )

        CALL thornado2amrex_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uDM, D )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDM ], &
                 D )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAM ], &
                 A )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPM ], &
                 P )

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

      CALL amrex_multifab_destroy( MF_uAM(iLevel) )
      CALL amrex_multifab_destroy( MF_uPM(iLevel) )

    END DO

  END SUBROUTINE ComputeDiagnosticFields_MHD_MF


  SUBROUTINE ComputeDiagnosticFields_MHD &
    ( iX_B0, iX_E0, iX_B1, iX_E1, P, A, D )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: P(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                               A(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( P(iNX,iX1,iX2,iX3,iPM_D   ), &
               Min_T, &
               A(iNX,iX1,iX2,iX3,iAM_Ye  ), &
               D(iNX,iX1,iX2,iX3,iDM_MinE) )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeDiagnosticFields_MHD


END MODULE MF_MHD_UtilitiesModule
