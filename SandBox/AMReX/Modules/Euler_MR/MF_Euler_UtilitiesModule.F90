MODULE MF_Euler_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX
  USE MeshModule, ONLY: &
    MeshX
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    Min_T
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    nGF
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
    nPF, &
    iAF_P, &
    iAF_Ye, &
    nAF, &
    iDF_MinE, &
    nDF
  USE Euler_UtilitiesModule, ONLY: &
    ComputeTimeStep_Euler, &
    ComputeFromConserved_Euler, &
    ComputeConserved_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    One
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    swX, &
    nX, &
    xL, &
    xR, &
    DEBUG
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    AllocateArray_X, &
    DeallocateArray_X
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_ComputeTimeStep_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeFromConserved_Euler_MF
  PUBLIC :: ComputeTimeStep_Euler_MF
  PUBLIC :: ComputeConserved_Euler_MF
  PUBLIC :: ComputeDiagnosticFields_Euler_MF

CONTAINS


  SUBROUTINE ComputeFromConserved_Euler_MF &
    ( MF_uGF, MF_uCF, MF_uPF, MF_uAF, &
      swXX_Option )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:), MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uPF(0:), MF_uAF(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

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
      !$OMP PRIVATE( MFI, BX, uGF, uCF, uPF, uAF, G, U, P, A, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, iX_B, iX_E )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uAF => MF_uAF(iLevel) % DataPtr( MFI )

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
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 P )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAF ], &
                 A )

        iX_B = iX_B0 - swXX
        iX_E = iX_E0 + swXX

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uCF, U )

        CALL ComputeFromConserved_Euler &
               ( iX_B, iX_E, iX_B1, iX_E1, G, U, P, A )

        CALL thornado2amrex_X( nPF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uPF, P )

        CALL thornado2amrex_X( nAF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uAF, A )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAF ], &
                 A )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 P )

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

    END DO

  END SUBROUTINE ComputeFromConserved_Euler_MF


  SUBROUTINE ComputeTimeStep_Euler_MF &
    ( MF_uGF, MF_uCF, CFL, TimeStepMin )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:), MF_uCF(0:)
    REAL(DP)            , INTENT(in)  :: CFL
    REAL(DP)            , INTENT(out) :: TimeStepMin(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4), iErr

    REAL(DP) :: TimeStep(0:nLevels-1)

    CALL TimersStart_AMReX( Timer_AMReX_ComputeTimeStep_Euler )

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() ) &
        WRITE(*,'(A)') 'CALL ComputeTimeStep_Euler_MF'

    END IF

    TimeStepMin = HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( MFI, BX, uGF, uCF, G, U, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, TimeStep ) &
      !$OMP REDUCTION( MIN:TimeStepMin )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

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
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, U )

        CALL ComputeTimeStep_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep(iLevel) )

        TimeStepMin(iLevel) = MIN( TimeStepMin(iLevel), TimeStep(iLevel) )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
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

    CALL TimersStop_AMReX( Timer_AMReX_ComputeTimeStep_Euler )

  END SUBROUTINE ComputeTimeStep_Euler_MF


  SUBROUTINE ComputeConserved_Euler_MF( MF_uGF, MF_uPF, MF_uAF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: A(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iNX, iX1, iX2, iX3

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( MFI, BX, uGF, uPF, uAF, uCF, G, P, A, U, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uAF => MF_uAF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

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
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 P )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAF ], &
                 A )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nPF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uPF, P )

        CALL amrex2thornado_X( nAF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uAF, A )

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          CALL ComputeConserved_Euler &
                 ( P(iNX,iX1,iX2,iX3,iPF_D ), &
                   P(iNX,iX1,iX2,iX3,iPF_V1), &
                   P(iNX,iX1,iX2,iX3,iPF_V3), &
                   P(iNX,iX1,iX2,iX3,iPF_V3), &
                   P(iNX,iX1,iX2,iX3,iPF_E ), &
                   P(iNX,iX1,iX2,iX3,iPF_Ne), &
                   U(iNX,iX1,iX2,iX3,iCF_D ), &
                   U(iNX,iX1,iX2,iX3,iCF_S1), &
                   U(iNX,iX1,iX2,iX3,iCF_S2), &
                   U(iNX,iX1,iX2,iX3,iCF_S3), &
                   U(iNX,iX1,iX2,iX3,iCF_E ), &
                   U(iNX,iX1,iX2,iX3,iCF_Ne), &
                   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                   A(iNX,iX1,iX2,iX3,iAF_P) )

        END DO
        END DO
        END DO
        END DO

        CALL thornado2amrex_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAF ], &
                 A )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 P )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO

  END SUBROUTINE ComputeConserved_Euler_MF


  SUBROUTINE ComputeDiagnosticFields_Euler_MF &
    ( MF_uGF, MF_uCF, MF_uDF, swXX_Option )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:), MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)
    INTEGER             , INTENT(in), OPTIONAL :: swXX_Option(3)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    TYPE(amrex_multifab) :: MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab) :: MF_uAF(0:nLevels-1)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDF(:,:,:,:)

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
             ( MF_uPF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nPF, swX )

      CALL amrex_multifab_build &
             ( MF_uAF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nAF, swX )

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( MFI, BX, uGF, uCF, uPF, uAF, uDF, G, U, P, A, D, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iLo_MF, iX_B, iX_E )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uAF => MF_uAF(iLevel) % DataPtr( MFI )
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
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 P )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAF ], &
                 A )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDF ], &
                 D )

        iX_B = iX_B0 - swXX
        iX_E = iX_E0 + swXX

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uCF, U )

        CALL amrex2thornado_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uDF, D )

        CALL ComputeFromConserved_Euler &
               ( iX_B, iX_E, iX_B1, iX_E1, G, U, P, A )

        CALL ComputeDiagnosticFields_Euler &
               ( iX_B, iX_E, iX_B1, iX_E1, G, U, P, A, D )

        CALL thornado2amrex_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, uDF, D )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nDF ], &
                 D )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nAF ], &
                 A )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 P )

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

      CALL amrex_multifab_destroy( MF_uAF(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF(iLevel) )

    END DO

  END SUBROUTINE ComputeDiagnosticFields_Euler_MF


  SUBROUTINE ComputeDiagnosticFields_Euler &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A, D )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout) :: U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(out)   :: P(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:), &
                               A(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:), &
                               D(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1       , nDOFX

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( P(iNX,iX1,iX2,iX3,iPF_D   ), &
               Min_T, &
               A(iNX,iX1,iX2,iX3,iAF_Ye  ), &
               D(iNX,iX1,iX2,iX3,iDF_MinE) )

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeDiagnosticFields_Euler


END MODULE MF_Euler_UtilitiesModule
