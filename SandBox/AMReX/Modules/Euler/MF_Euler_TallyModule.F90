MODULE MF_Euler_TallyModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_imultifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    swX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE MeshModule, ONLY: &
    MeshType
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm, &
    nGF
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nCF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE InputParsingModule, ONLY: &
    nLevels, &
    ProgramName, &
    UseTiling
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    AllocateArray_X, &
    DeallocateArray_X

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE ADM_Mass_Module, ONLY: &
    Calc_ADM_Mass

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_Euler_MF
  PUBLIC :: ComputeTally_Euler_MF
  PUBLIC :: IncrementOffGridTally_Euler_MF
  PUBLIC :: FinalizeTally_Euler_MF

  LOGICAL :: SuppressTally_Euler

  INTEGER, PARAMETER :: SL = 256

  CHARACTER(SL)    :: BaryonicMass_FileName
  REAL(DP), PUBLIC :: BaryonicMass_Initial
  REAL(DP), PUBLIC :: BaryonicMass_OffGrid
  REAL(DP)         :: BaryonicMass_Interior
  REAL(DP)         :: BaryonicMass_Interior_OMP
  REAL(DP)         :: BaryonicMass_Change

  CHARACTER(SL)    :: EulerMomentumX1_FileName
  REAL(DP), PUBLIC :: EulerMomentumX1_Initial
  REAL(DP), PUBLIC :: EulerMomentumX1_OffGrid
  REAL(DP)         :: EulerMomentumX1_Interior
  REAL(DP)         :: EulerMomentumX1_Interior_OMP
  REAL(DP)         :: EulerMomentumX1_Change

  CHARACTER(SL)    :: EulerMomentumX2_FileName
  REAL(DP), PUBLIC :: EulerMomentumX2_Initial
  REAL(DP), PUBLIC :: EulerMomentumX2_OffGrid
  REAL(DP)         :: EulerMomentumX2_Interior
  REAL(DP)         :: EulerMomentumX2_Interior_OMP
  REAL(DP)         :: EulerMomentumX2_Change

  CHARACTER(SL)    :: EulerMomentumX3_FileName
  REAL(DP), PUBLIC :: EulerMomentumX3_Initial
  REAL(DP), PUBLIC :: EulerMomentumX3_OffGrid
  REAL(DP)         :: EulerMomentumX3_Interior
  REAL(DP)         :: EulerMomentumX3_Interior_OMP
  REAL(DP)         :: EulerMomentumX3_Change

  CHARACTER(SL)    :: EulerEnergy_FileName
  REAL(DP), PUBLIC :: EulerEnergy_Initial
  REAL(DP), PUBLIC :: EulerEnergy_OffGrid
  REAL(DP)         :: EulerEnergy_Interior
  REAL(DP)         :: EulerEnergy_Interior_OMP
  REAL(DP)         :: EulerEnergy_Change

  CHARACTER(SL)    :: ElectronNumber_FileName
  REAL(DP), PUBLIC :: ElectronNumber_Initial
  REAL(DP), PUBLIC :: ElectronNumber_OffGrid
  REAL(DP)         :: ElectronNumber_Interior
  REAL(DP)         :: ElectronNumber_Interior_OMP
  REAL(DP)         :: ElectronNumber_Change

  CHARACTER(SL)    :: ADMMass_FileName
  REAL(DP), PUBLIC :: ADMMass_Initial
  REAL(DP), PUBLIC :: ADMMass_OffGrid
  REAL(DP), PUBLIC :: ADMMass_Interior
  REAL(DP)         :: ADMMass_Change

CONTAINS


  SUBROUTINE InitializeTally_Euler_MF &
    ( InitializeFromCheckpoint_Option )

    LOGICAL, INTENT(in), OPTIONAL :: &
      InitializeFromCheckpoint_Option

    CHARACTER(:), ALLOCATABLE :: TallyFileNameRoot_Euler

    CHARACTER(SL) :: FileNameRoot

    LOGICAL :: InitializeFromCheckpoint

    TYPE(amrex_parmparse) :: PP

    CHARACTER(SL) :: TimeLabel

    InitializeFromCheckpoint = .FALSE.
    IF( PRESENT( InitializeFromCheckpoint_Option ) ) &
      InitializeFromCheckpoint = InitializeFromCheckpoint_Option

    TallyFileNameRoot_Euler = TRIM( ProgramName )
    SuppressTally_Euler     = .FALSE.
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL pp % query( 'TallyFileNameRoot_Euler', &
                        TallyFileNameRoot_Euler )
      CALL pp % query( 'SuppressTally_Euler', &
                        SuppressTally_Euler )
    CALL amrex_parmparse_destroy( PP )

    IF( SuppressTally_Euler ) RETURN

    IF( amrex_parallel_ioprocessor() )THEN

      FileNameRoot = TRIM( TallyFileNameRoot_Euler )

      TimeLabel     &
        = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'

      ! --- Baryonic Mass ---

      BaryonicMass_FileName &
        = TRIM( FileNameRoot ) // '_BaryonicMass.dat'

      CALL CreateFile &
             ( BaryonicMass_FileName, UnitsDisplay % MassLabel, TimeLabel )

      ! --- Euler Momentum (X1) ---

      EulerMomentumX1_FileName &
        = TRIM( FileNameRoot ) // '_EulerMomentumX1.dat'

      CALL CreateFile &
             ( EulerMomentumX1_FileName, &
               UnitsDisplay % MomentumX1Label, TimeLabel )

      ! --- Euler Momentum (X2) ---

      EulerMomentumX2_FileName &
        = TRIM( FileNameRoot ) // '_EulerMomentumX2.dat'

      CALL CreateFile &
             ( EulerMomentumX2_FileName, &
               UnitsDisplay % MomentumX2Label, TimeLabel )

      ! --- Euler Momentum (X3) ---

      EulerMomentumX3_FileName &
        = TRIM( FileNameRoot ) // '_EulerMomentumX3.dat'

      CALL CreateFile &
             ( EulerMomentumX3_FileName, &
               UnitsDisplay % MomentumX3Label, TimeLabel )

      ! --- Euler Energy ---

      EulerEnergy_FileName &
        = TRIM( FileNameRoot ) // '_EulerEnergy.dat'

      CALL CreateFile &
             ( EulerEnergy_FileName, &
               UnitsDisplay % EnergyGlobalLabel, TimeLabel )

      ! --- Electron Number ---

      ElectronNumber_FileName &
        = TRIM( FileNameRoot ) // '_ElectronNumber.dat'

#ifdef MICROPHYSICS_WEAKLIB

      CALL CreateFile &
             ( ElectronNumber_FileName, &
               UnitsDisplay % ParticleDensityLabel, TimeLabel )

#endif

      ! --- ADM Mass ---

      ADMMass_FileName &
        = TRIM( FileNameRoot ) // '_ADMMass.dat'

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      CALL CreateFile &
             ( ADMMass_FileName, UnitsDisplay % EnergyGlobalLabel, TimeLabel )

#endif

    END IF

    BaryonicMass_Interior = Zero
    BaryonicMass_Change   = Zero

    EulerMomentumX1_Interior = Zero
    EulerMomentumX1_Change   = Zero

    EulerMomentumX2_Interior = Zero
    EulerMomentumX2_Change   = Zero

    EulerMomentumX3_Interior = Zero
    EulerMomentumX3_Change   = Zero

    EulerEnergy_Interior = Zero
    EulerEnergy_Change   = Zero

    ElectronNumber_Interior = Zero
    ElectronNumber_Change   = Zero

    ADMMass_Change = Zero

    IF( .NOT. InitializeFromCheckpoint )THEN

      BaryonicMass_Initial = Zero
      BaryonicMass_OffGrid = Zero

      EulerMomentumX1_Initial = Zero
      EulerMomentumX1_OffGrid = Zero

      EulerMomentumX2_Initial = Zero
      EulerMomentumX2_OffGrid = Zero

      EulerMomentumX3_Initial = Zero
      EulerMomentumX3_OffGrid = Zero

      EulerEnergy_Initial = Zero
      EulerEnergy_OffGrid = Zero

      ElectronNumber_Initial = Zero
      ElectronNumber_OffGrid = Zero

      ADMMass_Initial  = Zero
      ADMMass_OffGrid  = Zero
      ADMMass_Interior = Zero

    END IF

  END SUBROUTINE InitializeTally_Euler_MF


  SUBROUTINE ComputeTally_Euler_MF &
    ( Time, MF_uGF, MF_uCF, SetInitialValues_Option, &
      WriteTally_Option, FixInteriorADMMass_Option, Verbose_Option )

    REAL(DP),             INTENT(in) :: Time  (0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:)
    LOGICAL,              INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL,              INTENT(in), OPTIONAL :: WriteTally_Option
    LOGICAL,              INTENT(in), OPTIONAL :: FixInteriorADMMass_Option
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: SetInitialValues, FixInteriorADMMass
    LOGICAL :: WriteTally
    LOGICAL :: Verbose

    INTEGER                       :: iX_B0(3), iX_E0(3)
    INTEGER                       :: iLevel, iLo_MF(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF     (:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G     (:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U     (:,:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_FineMask

    TYPE(MeshType) :: MeshX(3)
    INTEGER        :: iNX, iX1, iX2, iX3
    REAL(DP)       :: d3X

    IF( SuppressTally_Euler ) RETURN

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    WriteTally = .TRUE.
    IF( PRESENT( WriteTally_Option ) ) &
      WriteTally = WriteTally_Option

    FixInteriorADMMass = .FALSE.
    IF( PRESENT( FixInteriorADMMass_Option ) ) &
      FixInteriorADMMass = FixInteriorADMMass_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    BaryonicMass_Interior    = Zero
    EulerMomentumX1_Interior = Zero
    EulerMomentumX2_Interior = Zero
    EulerMomentumX3_Interior = Zero
    EulerEnergy_Interior     = Zero
    ElectronNumber_Interior  = Zero
    IF( .NOT. FixInteriorADMMass ) &
      ADMMass_Interior = Zero

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uCF % BA, MF_uCF % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

      BaryonicMass_Interior_OMP    = Zero
      EulerMomentumX1_Interior_OMP = Zero
      EulerMomentumX2_Interior_OMP = Zero
      EulerMomentumX3_Interior_OMP = Zero
      EulerEnergy_Interior_OMP     = Zero
      ElectronNumber_Interior_OMP  = Zero

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( iX_B0, iX_E0, iLo_MF, &
      !$OMP          BX, MFI, FineMask, uGF, uCF, G, U, d3X ) &
      !$OMP REDUCTION( +:BaryonicMass_Interior_OMP, &
      !$OMP              EulerMomentumX1_Interior_OMP, &
      !$OMP              EulerMomentumX2_Interior_OMP, &
      !$OMP              EulerMomentumX3_Interior_OMP, &
      !$OMP              EulerEnergy_Interior_OMP, &
      !$OMP              ElectronNumber_Interior_OMP )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask   % DataPtr( MFI )
        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uCF      => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCF ], &
                 U )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uGF, G )
        CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uCF, U )

        d3X =   MeshX(1) % Width(iX_B0(1)) &
              * MeshX(2) % Width(iX_B0(2)) &
              * MeshX(3) % Width(iX_B0(3))

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          BaryonicMass_Interior_OMP &
            = BaryonicMass_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_D)

          EulerMomentumX1_Interior_OMP &
            = EulerMomentumX1_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_S1)

          EulerMomentumX2_Interior_OMP &
            = EulerMomentumX2_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_S2)

          EulerMomentumX3_Interior_OMP &
            = EulerMomentumX3_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_S3)

          EulerEnergy_Interior_OMP &
            = EulerEnergy_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_E)

          ElectronNumber_Interior_OMP &
            = ElectronNumber_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCF_Ne)

        END DO
        END DO
        END DO
        END DO

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCF ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      BaryonicMass_Interior &
        = BaryonicMass_Interior    + BaryonicMass_Interior_OMP
      EulerMomentumX1_Interior &
        = EulerMomentumX1_Interior + EulerMomentumX1_Interior_OMP
      EulerMomentumX2_Interior &
        = EulerMomentumX2_Interior + EulerMomentumX2_Interior_OMP
      EulerMomentumX3_Interior &
        = EulerMomentumX3_Interior + EulerMomentumX3_Interior_OMP
      EulerEnergy_Interior &
        = EulerEnergy_Interior     + EulerEnergy_Interior_OMP
      ElectronNumber_Interior &
        = ElectronNumber_Interior  + ElectronNumber_Interior_OMP

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    IF( .NOT. FixInteriorADMMass ) &
      CALL Calc_ADM_Mass( ADMMass_Interior )

#else

    ADMMass_Interior = Zero

#endif

    CALL amrex_parallel_reduce_sum( BaryonicMass_Interior    )
    CALL amrex_parallel_reduce_sum( EulerMomentumX1_Interior )
    CALL amrex_parallel_reduce_sum( EulerMomentumX2_Interior )
    CALL amrex_parallel_reduce_sum( EulerMomentumX3_Interior )
    CALL amrex_parallel_reduce_sum( EulerEnergy_Interior     )
    CALL amrex_parallel_reduce_sum( ElectronNumber_Interior  )

    IF( SetInitialValues )THEN

      BaryonicMass_Initial    = BaryonicMass_Interior
      EulerMomentumX1_Initial = EulerMomentumX1_Interior
      EulerMomentumX2_Initial = EulerMomentumX2_Interior
      EulerMomentumX3_Initial = EulerMomentumX3_Interior
      EulerEnergy_Initial     = EulerEnergy_Interior
      ElectronNumber_Initial  = ElectronNumber_Interior
      ADMMass_Initial         = ADMMass_Interior

    END IF

    ! --- dM = Minterior - Minitial + ( OffGrid_Inner - OffGrid_Outer ) ---

    BaryonicMass_Change &
      = BaryonicMass_Interior &
          - ( BaryonicMass_Initial    + BaryonicMass_OffGrid )

    EulerMomentumX1_Change &
      = EulerMomentumX1_Interior &
          - ( EulerMomentumX1_Initial + EulerMomentumX1_OffGrid )

    EulerMomentumX2_Change &
      = EulerMomentumX2_Interior &
          - ( EulerMomentumX2_Initial + EulerMomentumX2_OffGrid )

    EulerMomentumX3_Change &
      = EulerMomentumX3_Interior &
          - ( EulerMomentumX3_Initial + EulerMomentumX3_OffGrid )

    EulerEnergy_Change &
      = EulerEnergy_Interior &
          - ( EulerEnergy_Initial     + EulerEnergy_OffGrid )

    ElectronNumber_Change &
      = ElectronNumber_Interior &
          - ( ElectronNumber_Initial  + ElectronNumber_OffGrid )

    ADMMass_Change &
      = ADMMass_Interior &
          - ( ADMMass_Initial         + ADMMass_OffGrid )

    IF( WriteTally ) &
      CALL WriteTally_Euler( Time(0) )

    IF( Verbose ) CALL DisplayTally( Time(0) )

  END SUBROUTINE ComputeTally_Euler_MF


  SUBROUTINE IncrementOffGridTally_Euler_MF( dM )

    REAL(DP), INTENT(in) :: dM(1:,0:)

    INTEGER :: iLevel

    IF( SuppressTally_Euler ) RETURN

    DO iLevel = 0, nLevels-1

      BaryonicMass_OffGrid &
        = BaryonicMass_OffGrid    + dM(iCF_D ,iLevel)

      EulerMomentumX1_OffGrid &
        = EulerMomentumX1_OffGrid + dM(iCF_S1,iLevel)

      EulerMomentumX2_OffGrid &
        = EulerMomentumX2_OffGrid + dM(iCF_S2,iLevel)

      EulerMomentumX3_OffGrid &
        = EulerMomentumX3_OffGrid + dM(iCF_S3,iLevel)

      EulerEnergy_OffGrid &
        = EulerEnergy_OffGrid     + dM(iCF_E ,iLevel)

      ElectronNumber_OffGrid &
        = ElectronNumber_OffGrid  + dM(iCF_Ne,iLevel)

      ADMMass_OffGrid &
        = Zero

    END DO

  END SUBROUTINE IncrementOffGridTally_Euler_MF


  SUBROUTINE FinalizeTally_Euler_MF

    IF( SuppressTally_Euler ) RETURN

  END SUBROUTINE FinalizeTally_Euler_MF


  ! --- PRIVATE Subroutines ---


  SUBROUTINE WriteTally_Euler( Time )

    REAL(DP), INTENT(in) :: Time

    IF( amrex_parallel_ioprocessor() )THEN

      ! --- Baryonic Mass ---

      CALL WriteTallyToFile &
             ( BaryonicMass_FileName, Time, UnitsDisplay % TimeUnit, &
               BaryonicMass_Interior, &
               BaryonicMass_Initial, &
               BaryonicMass_OffGrid, &
               BaryonicMass_Change, &
               UnitsDisplay % MassUnit )

      ! --- Euler Momentum (X1) ---

      CALL WriteTallyToFile &
             ( EulerMomentumX1_FileName, Time, UnitsDisplay % TimeUnit, &
               EulerMomentumX1_Interior, &
               EulerMomentumX1_Initial, &
               EulerMomentumX1_OffGrid, &
               EulerMomentumX1_Change, &
               UnitsDisplay % MomentumX1Unit )

      ! --- Euler Momentum (X2) ---

      CALL WriteTallyToFile &
             ( EulerMomentumX2_FileName, Time, UnitsDisplay % TimeUnit, &
               EulerMomentumX2_Interior, &
               EulerMomentumX2_Initial, &
               EulerMomentumX2_OffGrid, &
               EulerMomentumX2_Change, &
               UnitsDisplay % MomentumX2Unit )

      ! --- Euler Momentum (X3) ---

      CALL WriteTallyToFile &
             ( EulerMomentumX3_FileName, Time, UnitsDisplay % TimeUnit, &
               EulerMomentumX3_Interior, &
               EulerMomentumX3_Initial, &
               EulerMomentumX3_OffGrid, &
               EulerMomentumX3_Change, &
               UnitsDisplay % MomentumX3Unit )

      ! --- Euler Energy ---

      CALL WriteTallyToFile &
             ( EulerEnergy_FileName, Time, UnitsDisplay % TimeUnit, &
               EulerEnergy_Interior, &
               EulerEnergy_Initial, &
               EulerEnergy_OffGrid, &
               EulerEnergy_Change, &
               UnitsDisplay % EnergyGlobalUnit )

#ifdef MICROPHYSICS_WEAKLIB

      ! --- Electron Number ---

      CALL WriteTallyToFile &
             ( ElectronNumber_FileName, Time, UnitsDisplay % TimeUnit, &
               ElectronNumber_Interior, &
               ElectronNumber_Initial, &
               ElectronNumber_OffGrid, &
               ElectronNumber_Change, &
               UnitsDisplay % ParticleDensityUnit )

#endif

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      ! --- ADM Mass ---

      CALL WriteTallyToFile &
             ( ADMMass_FileName, Time, UnitsDisplay % TimeUnit, &
               ADMMass_Interior, &
               ADMMass_Initial, &
               ADMMass_OffGrid, &
               ADMMass_Change, &
               UnitsDisplay % EnergyGlobalUnit )

#endif

    END IF

  END SUBROUTINE WriteTally_Euler


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A8,A,ES8.2E2,x,A)') &
        '', 'Euler Tally. t = ', &
        Time / UnitsDisplay % TimeUnit, &
        UnitsDisplay % TimeLabel

      CALL WriteTallyToScreen( 'Baryonic Mass', &
                               BaryonicMass_Interior, &
                               BaryonicMass_Initial, &
                               BaryonicMass_OffGrid, &
                               BaryonicMass_Change, &
                               UnitsDisplay % MassUnit, &
                               UnitsDisplay % MassLabel )

      CALL WriteTallyToScreen( 'Euler Momentum (X1)', &
                               EulerMomentumX1_Interior, &
                               EulerMomentumX1_Initial, &
                               EulerMomentumX1_OffGrid, &
                               EulerMomentumX1_Change, &
                               UnitsDisplay % MomentumX1Unit, &
                               UnitsDisplay % MomentumX1Label )

      CALL WriteTallyToScreen( 'Euler Momentum (X2)', &
                               EulerMomentumX2_Interior, &
                               EulerMomentumX2_Initial, &
                               EulerMomentumX2_OffGrid, &
                               EulerMomentumX2_Change, &
                               UnitsDisplay % MomentumX2Unit, &
                               UnitsDisplay % MomentumX2Label )

      CALL WriteTallyToScreen( 'Euler Momentum (X3)', &
                               EulerMomentumX3_Interior, &
                               EulerMomentumX3_Initial, &
                               EulerMomentumX3_OffGrid, &
                               EulerMomentumX3_Change, &
                               UnitsDisplay % MomentumX3Unit, &
                               UnitsDisplay % MomentumX3Label )

      CALL WriteTallyToScreen( 'Euler Energy', &
                               EulerEnergy_Interior, &
                               EulerEnergy_Initial, &
                               EulerEnergy_OffGrid, &
                               EulerEnergy_Change, &
                               UnitsDisplay % EnergyGlobalUnit, &
                               UnitsDisplay % EnergyGlobalLabel )

#ifdef MICROPHYSICS_WEAKLIB

      CALL WriteTallyToScreen( 'Electron Number', &
                               ElectronNumber_Interior, &
                               ElectronNumber_Initial, &
                               ElectronNumber_OffGrid, &
                               ElectronNumber_Change, &
                               UnitsDisplay % ParticleDensityUnit, &
                               UnitsDisplay % ParticleDensityLabel )

#endif

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      CALL WriteTallyToScreen( 'ADM Mass', &
                               ADMMass_Interior, &
                               ADMMass_Initial, &
                               ADMMass_OffGrid, &
                               ADMMass_Change, &
                               UnitsDisplay % EnergyGlobalUnit, &
                               UnitsDisplay % EnergyGlobalLabel )

#endif

      WRITE(*,*)

    END IF

  END SUBROUTINE DisplayTally


  RECURSIVE SUBROUTINE CheckFileExistenceAndAppend( FileName, IntSuffix_Option )

    CHARACTER(LEN=SL), INTENT(inout) :: FileName
    INTEGER          , INTENT(inout), OPTIONAL :: IntSuffix_Option

    LOGICAL :: IsFile
    INTEGER :: IntSuffix
    INTEGER :: SL_T

    IntSuffix = 1
    IF( PRESENT( IntSuffix_Option ) ) &
      IntSuffix = IntSuffix_Option

    SL_T = LEN( TRIM( FileName ) )

    INQUIRE( FILE = TRIM( FileName ), EXIST = IsFile )

    IF( IsFile )THEN

      IF( FileName(SL_T-3:SL_T) .EQ. '.dat' )THEN

        WRITE(FileName,'(A,A,I2.2)') TRIM( FileName ), '_', IntSuffix

      ELSE

        WRITE(FileName(SL_T-1:SL_T),'(I2.2)') IntSuffix

      END IF

      IntSuffix = IntSuffix + 1

      CALL CheckFileExistenceAndAppend &
             ( FileName, IntSuffix_Option = IntSuffix )

    END IF

  END SUBROUTINE CheckFileExistenceAndAppend


  SUBROUTINE CreateFile( FileName, UnitsLabel, TimeLabel )

    CHARACTER(*), INTENT(inout) :: FileName
    CHARACTER(*), INTENT(in)    :: UnitsLabel, TimeLabel

    INTEGER       :: FileUnit
    CHARACTER(SL) :: InteriorLabel, InitialLabel, OffGridLabel, ChangeLabel

    InteriorLabel &
      = 'Interior [' // TRIM( UnitsLabel ) // ']'
    OffGridLabel  &
      = 'Off Grid [' // TRIM( UnitsLabel ) // ']'
    InitialLabel  &
      = 'Initial ['  // TRIM( UnitsLabel ) // ']'
    ChangeLabel   &
      = 'Change ['   // TRIM( UnitsLabel ) // ']'

    CALL CheckFileExistenceAndAppend( FileName )

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ) )

    WRITE(FileUnit,'(5(A25,x))') &
      TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
      TRIM( InitialLabel ), TRIM( ChangeLabel )

    CLOSE( FileUnit )

  END SUBROUTINE CreateFile


  SUBROUTINE WriteTallyToScreen &
    ( FieldName, Interior, Initial, OffGrid, Change, Units, Label )

    CHARACTER(*), INTENT(in) :: FieldName, Label
    REAL(DP)    , INTENT(in) :: Interior, Initial, OffGrid, Change, Units

    CHARACTER(32) :: FMT

    FMT = '(6x,A40,ES15.7E3,x,A)'

    WRITE(*,*)
    WRITE(*,TRIM(FMT)) &
      TRIM( FieldName ) // ' Interior.: ', Interior / Units, TRIM( Label )
    WRITE(*,TRIM(FMT)) &
      TRIM( FieldName ) // ' Initial..: ', Initial  / Units, TRIM( Label )
    WRITE(*,TRIM(FMT)) &
      TRIM( FieldName ) // ' Off Grid.: ', OffGrid  / Units, TRIM( Label )
    WRITE(*,TRIM(FMT)) &
      TRIM( FieldName ) // ' Change...: ', Change   / Units, TRIM( Label )

  END SUBROUTINE WriteTallyToScreen


  SUBROUTINE WriteTallyToFile &
    ( FileName, Time, TimeUnit, &
      Interior, Initial, OffGrid, Change, Units )

    CHARACTER(*), INTENT(in) :: FileName
    REAL(DP)    , INTENT(in) :: Time, TimeUnit, &
                                Interior, Initial, OffGrid, Change, Units

    INTEGER       :: FileUnit
    CHARACTER(32) :: FMT

    FMT = '(5(ES25.16E3,1x))'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, TRIM(FMT) ) &
      Time     / TimeUnit, &
      Interior / Units, &
      OffGrid  / Units, &
      Initial  / Units, &
      Change   / Units

    CLOSE( FileUnit )

  END SUBROUTINE WriteTallyToFile


END MODULE MF_Euler_TallyModule
