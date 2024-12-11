MODULE MF_MHD_TallyModule

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
    iAM_Tem33, &
    nAM

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

  PUBLIC :: InitializeTally_MHD_MF
  PUBLIC :: ComputeTally_MHD_MF
  PUBLIC :: IncrementOffGridTally_MHD_MF
  PUBLIC :: FinalizeTally_MHD_MF

  LOGICAL :: SuppressTally

  INTEGER, PARAMETER :: SL = 256

  CHARACTER(SL)    :: BaryonicMass_FileName
  REAL(DP), PUBLIC :: BaryonicMass_Initial
  REAL(DP), PUBLIC :: BaryonicMass_OffGrid
  REAL(DP)         :: BaryonicMass_Interior
  REAL(DP)         :: BaryonicMass_Interior_OMP
  REAL(DP)         :: BaryonicMass_Change

  CHARACTER(SL)    :: MHDMomentumX1_FileName
  REAL(DP), PUBLIC :: MHDMomentumX1_Initial
  REAL(DP), PUBLIC :: MHDMomentumX1_OffGrid
  REAL(DP)         :: MHDMomentumX1_Interior
  REAL(DP)         :: MHDMomentumX1_Interior_OMP
  REAL(DP)         :: MHDMomentumX1_Change

  CHARACTER(SL)    :: MHDMomentumX2_FileName
  REAL(DP), PUBLIC :: MHDMomentumX2_Initial
  REAL(DP), PUBLIC :: MHDMomentumX2_OffGrid
  REAL(DP)         :: MHDMomentumX2_Interior
  REAL(DP)         :: MHDMomentumX2_Interior_OMP
  REAL(DP)         :: MHDMomentumX2_Change

  CHARACTER(SL)    :: MHDMomentumX3_FileName
  REAL(DP), PUBLIC :: MHDMomentumX3_Initial
  REAL(DP), PUBLIC :: MHDMomentumX3_OffGrid
  REAL(DP)         :: MHDMomentumX3_Interior
  REAL(DP)         :: MHDMomentumX3_Interior_OMP
  REAL(DP)         :: MHDMomentumX3_Change

  CHARACTER(SL)    :: MHDEnergy_FileName
  REAL(DP), PUBLIC :: MHDEnergy_Initial
  REAL(DP), PUBLIC :: MHDEnergy_OffGrid
  REAL(DP)         :: MHDEnergy_Interior
  REAL(DP)         :: MHDEnergy_Interior_OMP
  REAL(DP)         :: MHDEnergy_Change

  CHARACTER(SL)    :: ElectronNumber_FileName
  REAL(DP), PUBLIC :: ElectronNumber_Initial
  REAL(DP), PUBLIC :: ElectronNumber_OffGrid
  REAL(DP)         :: ElectronNumber_Interior
  REAL(DP)         :: ElectronNumber_Interior_OMP
  REAL(DP)         :: ElectronNumber_Change

  CHARACTER(SL)    :: MHDMagFieldX1_FileName
  REAL(DP), PUBLIC :: MHDMagFieldX1_Initial
  REAL(DP), PUBLIC :: MHDMagFieldX1_OffGrid
  REAL(DP)         :: MHDMagFieldX1_Interior
  REAL(DP)         :: MHDMagFieldX1_Interior_OMP
  REAL(DP)         :: MHDMagFieldX1_Change

  CHARACTER(SL)    :: MHDMagFieldX2_FileName
  REAL(DP), PUBLIC :: MHDMagFieldX2_Initial
  REAL(DP), PUBLIC :: MHDMagFieldX2_OffGrid
  REAL(DP)         :: MHDMagFieldX2_Interior
  REAL(DP)         :: MHDMagFieldX2_Interior_OMP
  REAL(DP)         :: MHDMagFieldX2_Change

  CHARACTER(SL)    :: MHDMagFieldX3_FileName
  REAL(DP), PUBLIC :: MHDMagFieldX3_Initial
  REAL(DP), PUBLIC :: MHDMagFieldX3_OffGrid
  REAL(DP)         :: MHDMagFieldX3_Interior
  REAL(DP)         :: MHDMagFieldX3_Interior_OMP
  REAL(DP)         :: MHDMagFieldX3_Change

  CHARACTER(SL)    :: MHDCleaningField_FileName
  REAL(DP), PUBLIC :: MHDCleaningField_Initial
  REAL(DP), PUBLIC :: MHDCleaningField_OffGrid
  REAL(DP)         :: MHDCleaningField_Interior
  REAL(DP)         :: MHDCleaningField_Interior_OMP
  REAL(DP)         :: MHDCleaningField_Change

  CHARACTER(SL)    :: ADMMass_FileName
  REAL(DP), PUBLIC :: ADMMass_Initial
  REAL(DP), PUBLIC :: ADMMass_OffGrid
  REAL(DP), PUBLIC :: ADMMass_Interior
  REAL(DP)         :: ADMMass_Change

  CHARACTER(SL)    :: Tem33_FileName
  REAL(DP), PUBLIC :: Tem33_Initial
  REAL(DP), PUBLIC :: Tem33_OffGrid
  REAL(DP)         :: Tem33_Interior
  REAL(DP)         :: Tem33_Interior_OMP
  REAL(DP)         :: Tem33_Change

CONTAINS


  SUBROUTINE InitializeTally_MHD_MF &
    ( SuppressTally_Option, InitializeFromCheckpoint_Option )

    LOGICAL, INTENT(in), OPTIONAL :: &
      SuppressTally_Option, &
      InitializeFromCheckpoint_Option

    CHARACTER(:), ALLOCATABLE :: TallyFileNameRoot_MHD

    CHARACTER(SL) :: FileNameRoot

    LOGICAL :: InitializeFromCheckpoint

    TYPE(amrex_parmparse) :: PP

    CHARACTER(SL) :: TimeLabel

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

    InitializeFromCheckpoint = .FALSE.
    IF( PRESENT( InitializeFromCheckpoint_Option ) ) &
      InitializeFromCheckpoint = InitializeFromCheckpoint_Option

    TallyFileNameRoot_MHD = TRIM( ProgramName )
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL pp % query( 'TallyFileNameRoot_MHD', &
                        TallyFileNameRoot_MHD )
    CALL amrex_parmparse_destroy( PP )

    IF( SuppressTally ) RETURN

    IF( amrex_parallel_ioprocessor() )THEN

      FileNameRoot = TRIM( TallyFileNameRoot_MHD )

      TimeLabel     &
        = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'

      ! --- Baryonic Mass ---

      BaryonicMass_FileName &
        = TRIM( FileNameRoot ) // '_BaryonicMass.dat'

      CALL CreateFile &
             ( BaryonicMass_FileName, UnitsDisplay % MassLabel, TimeLabel )

      ! --- MHD Momentum (X1) ---

      MHDMomentumX1_FileName &
        = TRIM( FileNameRoot ) // '_MHDMomentumX1.dat'

      CALL CreateFile &
             ( MHDMomentumX1_FileName, &
               UnitsDisplay % MomentumX1Label, TimeLabel )

      ! --- MHD Momentum (X2) ---

      MHDMomentumX2_FileName &
        = TRIM( FileNameRoot ) // '_MHDMomentumX2.dat'

      CALL CreateFile &
             ( MHDMomentumX2_FileName, &
               UnitsDisplay % MomentumX2Label, TimeLabel )

      ! --- MHD Momentum (X3) ---

      MHDMomentumX3_FileName &
        = TRIM( FileNameRoot ) // '_MHDMomentumX3.dat'

      CALL CreateFile &
             ( MHDMomentumX3_FileName, &
               UnitsDisplay % MomentumX3Label, TimeLabel )

      ! --- MHD Energy ---

      MHDEnergy_FileName &
        = TRIM( FileNameRoot ) // '_MHDEnergy.dat'

      CALL CreateFile &
             ( MHDEnergy_FileName, &
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

      ! --- Tem33 ---

      Tem33_FileName &
        = TRIM( FileNameRoot ) // '_Tem33.dat'


    END IF

    BaryonicMass_Interior = Zero
    BaryonicMass_Change   = Zero

    MHDMomentumX1_Interior = Zero
    MHDMomentumX1_Change   = Zero

    MHDMomentumX2_Interior = Zero
    MHDMomentumX2_Change   = Zero

    MHDMomentumX3_Interior = Zero
    MHDMomentumX3_Change   = Zero

    MHDEnergy_Interior = Zero
    MHDEnergy_Change   = Zero

    ElectronNumber_Interior = Zero
    ElectronNumber_Change   = Zero

    ADMMass_Change = Zero

    Tem33_Interior = Zero
    Tem33_Change   = Zero

    IF( .NOT. InitializeFromCheckpoint )THEN

      BaryonicMass_Initial = Zero
      BaryonicMass_OffGrid = Zero

      MHDMomentumX1_Initial = Zero
      MHDMomentumX1_OffGrid = Zero

      MHDMomentumX2_Initial = Zero
      MHDMomentumX2_OffGrid = Zero

      MHDMomentumX3_Initial = Zero
      MHDMomentumX3_OffGrid = Zero

      MHDEnergy_Initial = Zero
      MHDEnergy_OffGrid = Zero

      ElectronNumber_Initial = Zero
      ElectronNumber_OffGrid = Zero

      ADMMass_Initial  = Zero
      ADMMass_OffGrid  = Zero
      ADMMass_Interior = Zero

      Tem33_Initial = Zero
      Tem33_OffGrid = Zero

    END IF

  END SUBROUTINE InitializeTally_MHD_MF


  SUBROUTINE ComputeTally_MHD_MF &
    ( Time, MF_uGF, MF_uCM, MF_uAM, SetInitialValues_Option, &
      WriteTally_Option, FixInteriorADMMass_Option, Verbose_Option )

    REAL(DP),             INTENT(in) :: Time  (0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCM(0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uAM(0:)
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
    REAL(DP), CONTIGUOUS, POINTER :: uCM     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAM     (:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G     (:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U     (:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: A     (:,:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_FineMask

    TYPE(MeshType) :: MeshX(3)
    INTEGER        :: iNX, iX1, iX2, iX3
    REAL(DP)       :: d3X

    IF( SuppressTally ) RETURN

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
    MHDMomentumX1_Interior = Zero
    MHDMomentumX2_Interior = Zero
    MHDMomentumX3_Interior = Zero
    MHDEnergy_Interior     = Zero
    ElectronNumber_Interior  = Zero
    Tem33_Interior           = Zero
    IF( .NOT. FixInteriorADMMass ) &
      ADMMass_Interior = Zero

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uCM % BA, MF_uCM % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

      BaryonicMass_Interior_OMP    = Zero
      MHDMomentumX1_Interior_OMP = Zero
      MHDMomentumX2_Interior_OMP = Zero
      MHDMomentumX3_Interior_OMP = Zero
      MHDEnergy_Interior_OMP     = Zero
      ElectronNumber_Interior_OMP  = Zero
      Tem33_Interior_OMP           = Zero

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( iX_B0, iX_E0, iLo_MF, &
      !$OMP          BX, MFI, FineMask, uGF, uCM, G, U, d3X ) &
      !$OMP REDUCTION( +:BaryonicMass_Interior_OMP, &
      !$OMP              MHDMomentumX1_Interior_OMP, &
      !$OMP              MHDMomentumX2_Interior_OMP, &
      !$OMP              MHDMomentumX3_Interior_OMP, &
      !$OMP              MHDEnergy_Interior_OMP, &
      !$OMP              ElectronNumber_Interior_OMP )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask   % DataPtr( MFI )
        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uCM      => MF_uCM(iLevel) % DataPtr( MFI )
        uAM      => MF_uAM(iLevel) % DataPtr( MFI )

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
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCM ], &
                 U )

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nAM ], &
                 A )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uGF, G )
        CALL amrex2thornado_X( nCM, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uCM, U )
        CALL amrex2thornado_X( nAM, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uAM, A )

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
                    * U(iNX,iX1,iX2,iX3,iCM_D)

          MHDMomentumX1_Interior_OMP &
            = MHDMomentumX1_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_S1)

          MHDMomentumX2_Interior_OMP &
            = MHDMomentumX2_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_S2)

          MHDMomentumX3_Interior_OMP &
            = MHDMomentumX3_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_S3)

          MHDEnergy_Interior_OMP &
            = MHDEnergy_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_E)

          ElectronNumber_Interior_OMP &
            = ElectronNumber_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * U(iNX,iX1,iX2,iX3,iCM_Ne)

          Tem33_Interior_OMP &
            = Tem33_Interior_OMP &
                + d3X &
                    * WeightsX_q(iNX) &
                    * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * A(iNX,iX1,iX2,iX3,iAM_Tem33)

        END DO
        END DO
        END DO
        END DO

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nAM ], &
                 A )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCM ], &
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
      MHDMomentumX1_Interior &
        = MHDMomentumX1_Interior + MHDMomentumX1_Interior_OMP
      MHDMomentumX2_Interior &
        = MHDMomentumX2_Interior + MHDMomentumX2_Interior_OMP
      MHDMomentumX3_Interior &
        = MHDMomentumX3_Interior + MHDMomentumX3_Interior_OMP
      MHDEnergy_Interior &
        = MHDEnergy_Interior     + MHDEnergy_Interior_OMP
      ElectronNumber_Interior &
        = ElectronNumber_Interior  + ElectronNumber_Interior_OMP
      Tem33_Interior &
        = Tem33_Interior + Tem33_Interior_OMP

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
    CALL amrex_parallel_reduce_sum( MHDMomentumX1_Interior )
    CALL amrex_parallel_reduce_sum( MHDMomentumX2_Interior )
    CALL amrex_parallel_reduce_sum( MHDMomentumX3_Interior )
    CALL amrex_parallel_reduce_sum( MHDEnergy_Interior     )
    CALL amrex_parallel_reduce_sum( ElectronNumber_Interior  )
    CALL amrex_parallel_reduce_sum( Tem33_Interior )

    IF( SetInitialValues )THEN

      BaryonicMass_Initial    = BaryonicMass_Interior
      MHDMomentumX1_Initial = MHDMomentumX1_Interior
      MHDMomentumX2_Initial = MHDMomentumX2_Interior
      MHDMomentumX3_Initial = MHDMomentumX3_Interior
      MHDEnergy_Initial     = MHDEnergy_Interior
      ElectronNumber_Initial  = ElectronNumber_Interior
      ADMMass_Initial         = ADMMass_Interior
      Tem33_Initial         = Tem33_Interior

    END IF

    ! --- dM = Minterior - Minitial + ( OffGrid_Outer - OffGrid_Inner ) ---

    BaryonicMass_Change &
      = BaryonicMass_Interior &
          - BaryonicMass_Initial    + BaryonicMass_OffGrid

    MHDMomentumX1_Change &
      = MHDMomentumX1_Interior &
          - MHDMomentumX1_Initial + MHDMomentumX1_OffGrid

    MHDMomentumX2_Change &
      = MHDMomentumX2_Interior &
          - MHDMomentumX2_Initial + MHDMomentumX2_OffGrid

    MHDMomentumX3_Change &
      = MHDMomentumX3_Interior &
          - MHDMomentumX3_Initial + MHDMomentumX3_OffGrid

    MHDEnergy_Change &
      = MHDEnergy_Interior &
          - MHDEnergy_Initial     + MHDEnergy_OffGrid

    ElectronNumber_Change &
      = ElectronNumber_Interior &
          - ElectronNumber_Initial  + ElectronNumber_OffGrid

    ADMMass_Change &
      = ADMMass_Interior &
          - ADMMass_Initial + ADMMass_OffGrid

    Tem33_Change &
      = Tem33_Interior &
          - Tem33_Initial + Tem33_OffGrid

    IF( WriteTally ) &
      CALL WriteTally_MHD( Time(0) )

    IF( Verbose ) CALL DisplayTally( Time(0) )

  END SUBROUTINE ComputeTally_MHD_MF


  SUBROUTINE IncrementOffGridTally_MHD_MF( dM )

    REAL(DP), INTENT(in) :: dM(1:,0:)

    INTEGER :: iLevel

    IF( SuppressTally ) RETURN

    DO iLevel = 0, nLevels-1

      BaryonicMass_OffGrid &
        = BaryonicMass_OffGrid    + dM(iCM_D ,iLevel)

      MHDMomentumX1_OffGrid &
        = MHDMomentumX1_OffGrid + dM(iCM_S1,iLevel)

      MHDMomentumX2_OffGrid &
        = MHDMomentumX2_OffGrid + dM(iCM_S2,iLevel)

      MHDMomentumX3_OffGrid &
        = MHDMomentumX3_OffGrid + dM(iCM_S3,iLevel)

      MHDEnergy_OffGrid &
        = MHDEnergy_OffGrid     + dM(iCM_E ,iLevel)

      ElectronNumber_OffGrid &
        = ElectronNumber_OffGrid  + dM(iCM_Ne,iLevel)

      ADMMass_OffGrid &
        = Zero

      Tem33_OffGrid &
        = Tem33_OffGrid + dM(iAM_Tem33,iLevel)

    END DO

  END SUBROUTINE IncrementOffGridTally_MHD_MF


  SUBROUTINE FinalizeTally_MHD_MF

    IF( SuppressTally ) RETURN

  END SUBROUTINE FinalizeTally_MHD_MF


  ! --- PRIVATE Subroutines ---


  SUBROUTINE WriteTally_MHD( Time )

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

      ! --- MHD Momentum (X1) ---

      CALL WriteTallyToFile &
             ( MHDMomentumX1_FileName, Time, UnitsDisplay % TimeUnit, &
               MHDMomentumX1_Interior, &
               MHDMomentumX1_Initial, &
               MHDMomentumX1_OffGrid, &
               MHDMomentumX1_Change, &
               UnitsDisplay % MomentumX1Unit )

      ! --- MHD Momentum (X2) ---

      CALL WriteTallyToFile &
             ( MHDMomentumX2_FileName, Time, UnitsDisplay % TimeUnit, &
               MHDMomentumX2_Interior, &
               MHDMomentumX2_Initial, &
               MHDMomentumX2_OffGrid, &
               MHDMomentumX2_Change, &
               UnitsDisplay % MomentumX2Unit )

      ! --- MHD Momentum (X3) ---

      CALL WriteTallyToFile &
             ( MHDMomentumX3_FileName, Time, UnitsDisplay % TimeUnit, &
               MHDMomentumX3_Interior, &
               MHDMomentumX3_Initial, &
               MHDMomentumX3_OffGrid, &
               MHDMomentumX3_Change, &
               UnitsDisplay % MomentumX3Unit )

      ! --- MHD Energy ---

      CALL WriteTallyToFile &
             ( MHDEnergy_FileName, Time, UnitsDisplay % TimeUnit, &
               MHDEnergy_Interior, &
               MHDEnergy_Initial, &
               MHDEnergy_OffGrid, &
               MHDEnergy_Change, &
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

      ! --- Tem33 ---

      CALL WriteTallyToFile &
             ( Tem33_FileName, Time, UnitsDisplay % TimeUnit, &
               Tem33_Interior, &
               Tem33_Initial, &
               Tem33_OffGrid, &
               Tem33_Change, &
               UnitsDisplay % EnergyGlobalUnit )

    END IF

  END SUBROUTINE WriteTally_MHD


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A8,A,ES8.2E2,x,A)') &
        '', 'MHD Tally. t = ', &
        Time / UnitsDisplay % TimeUnit, &
        UnitsDisplay % TimeLabel

      CALL WriteTallyToScreen( 'Baryonic Mass', &
                               BaryonicMass_Interior, &
                               BaryonicMass_Initial, &
                               BaryonicMass_OffGrid, &
                               BaryonicMass_Change, &
                               UnitsDisplay % MassUnit, &
                               UnitsDisplay % MassLabel )

      CALL WriteTallyToScreen( 'MHD Momentum (X1)', &
                               MHDMomentumX1_Interior, &
                               MHDMomentumX1_Initial, &
                               MHDMomentumX1_OffGrid, &
                               MHDMomentumX1_Change, &
                               UnitsDisplay % MomentumX1Unit, &
                               UnitsDisplay % MomentumX1Label )

      CALL WriteTallyToScreen( 'MHD Momentum (X2)', &
                               MHDMomentumX2_Interior, &
                               MHDMomentumX2_Initial, &
                               MHDMomentumX2_OffGrid, &
                               MHDMomentumX2_Change, &
                               UnitsDisplay % MomentumX2Unit, &
                               UnitsDisplay % MomentumX2Label )

      CALL WriteTallyToScreen( 'MHD Momentum (X3)', &
                               MHDMomentumX3_Interior, &
                               MHDMomentumX3_Initial, &
                               MHDMomentumX3_OffGrid, &
                               MHDMomentumX3_Change, &
                               UnitsDisplay % MomentumX3Unit, &
                               UnitsDisplay % MomentumX3Label )

      CALL WriteTallyToScreen( 'MHD Energy', &
                               MHDEnergy_Interior, &
                               MHDEnergy_Initial, &
                               MHDEnergy_OffGrid, &
                               MHDEnergy_Change, &
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

      CALL WriteTallyToScreen( 'Tem33', &
                               Tem33_Interior, &
                               Tem33_Initial, &
                               Tem33_OffGrid, &
                               Tem33_Change, &
                               UnitsDisplay % EnergyGlobalUnit, &
                               UnitsDisplay % EnergyGlobalLabel )

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


END MODULE MF_MHD_TallyModule
