MODULE MF_TwoMoment_TallyModule

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
    nDOFX, nDOFE, nDOFZ, &
    iE_B0, iE_E0
  USE ReferenceElementModule, ONLY: &
    Weights_q
  USE UnitsModule, ONLY: &
    UnitsActive, &
    SpeedOfLight, &
    PlanckConstant, &
    UnitsDisplay
  USE MeshModule, ONLY: &
    MeshType, &
    MeshE
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE GeometryFieldsModuleE, ONLY: &
    uGE, iGE_Ep2, iGE_Ep3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, LeptonNumber, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, Zero, One, FourPi
  USE InputParsingModule, ONLY: &
    nLevels, &
    ProgramName, &
    nE, &
    UseTiling, &
    iRestart
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_Z, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_TwoMoment_MF
  PUBLIC :: ComputeTally_TwoMoment_MF
  PUBLIC :: IncrementOffGridTally_TwoMoment_MF
  PUBLIC :: IncrementPositivityLimiterTally_TwoMoment_MF
  PUBLIC :: FinalizeTally_TwoMoment_MF
  PUBLIC :: WriteTallyCheckpoint_TwoMoment
  PUBLIC :: ReadTallyCheckpoint_TwoMoment

  LOGICAL :: SuppressTally_TwoMoment

  INTEGER, PARAMETER :: SL = 256

  REAL(DP) :: hc3

  CHARACTER(SL) :: TallyChkFileNameRoot

  CHARACTER(SL)    :: NeutrinoLeptonNumber_FileName
  REAL(DP), PUBLIC :: NeutrinoLeptonNumber_Initial
  REAL(DP), PUBLIC :: NeutrinoLeptonNumber_OffGrid
  REAL(DP), PUBLIC :: NeutrinoLeptonNumber_Interior
  REAL(DP)         :: NeutrinoLeptonNumber_Change

  CHARACTER(SL)    :: NeutrinoEnergy_FileName
  REAL(DP), PUBLIC :: NeutrinoEnergy_Initial
  REAL(DP), PUBLIC :: NeutrinoEnergy_OffGrid
  REAL(DP)         :: NeutrinoEnergy_Interior
  REAL(DP)         :: NeutrinoEnergy_Change

  CHARACTER(SL)    :: NeutrinoMomentumX1_FileName
  REAL(DP), PUBLIC :: NeutrinoMomentumX1_Initial
  REAL(DP), PUBLIC :: NeutrinoMomentumX1_OffGrid
  REAL(DP)         :: NeutrinoMomentumX1_Interior
  REAL(DP)         :: NeutrinoMomentumX1_Change

  CHARACTER(SL)    :: NeutrinoMomentumX2_FileName
  REAL(DP), PUBLIC :: NeutrinoMomentumX2_Initial
  REAL(DP), PUBLIC :: NeutrinoMomentumX2_OffGrid
  REAL(DP)         :: NeutrinoMomentumX2_Interior
  REAL(DP)         :: NeutrinoMomentumX2_Change

  CHARACTER(SL)    :: NeutrinoMomentumX3_FileName
  REAL(DP), PUBLIC :: NeutrinoMomentumX3_Initial
  REAL(DP), PUBLIC :: NeutrinoMomentumX3_OffGrid
  REAL(DP)         :: NeutrinoMomentumX3_Interior
  REAL(DP)         :: NeutrinoMomentumX3_Change

  REAL(DP), PUBLIC :: NeutrinoEnergy_PL
  REAL(DP), PUBLIC :: NeutrinoMomentumX1_PL
  REAL(DP), PUBLIC :: NeutrinoMomentumX2_PL
  REAL(DP), PUBLIC :: NeutrinoMomentumX3_PL

CONTAINS


  SUBROUTINE InitializeTally_TwoMoment_MF &
    ( InitializeFromCheckpoint_Option )

    LOGICAL, INTENT(in), OPTIONAL :: InitializeFromCheckpoint_Option

    CHARACTER(:), ALLOCATABLE :: TallyFileNameRoot_TwoMoment
    CHARACTER(SL)             :: FileNameRoot
    CHARACTER(SL)             :: TimeLabel
    LOGICAL                   :: InitializeFromCheckpoint
    TYPE(amrex_parmparse)     :: PP

    InitializeFromCheckpoint = .FALSE.
    IF( PRESENT( InitializeFromCheckpoint_Option ) ) &
      InitializeFromCheckpoint = InitializeFromCheckpoint_Option

    TallyFileNameRoot_TwoMoment = TRIM( ProgramName )
    SuppressTally_TwoMoment     = .FALSE.
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % query( 'TallyFileNameRoot_TwoMoment', &
                        TallyFileNameRoot_TwoMoment )
      CALL PP % query( 'SuppressTally_TwoMoment', &
                        SuppressTally_TwoMoment )
    CALL amrex_parmparse_destroy( PP )

    !IF( amrex_parallel_ioprocessor() ) &
    !  WRITE(*,'(A,L2,2x,A)') 'InitializeTally_TwoMoment_MF: Suppress =', &
    !    SuppressTally_TwoMoment, TRIM( TallyFileNameRoot_TwoMoment )

    IF( SuppressTally_TwoMoment ) RETURN

    IF( UnitsActive )THEN
      hc3 = ( PlanckConstant * SpeedOfLight )**3
    ELSE
      hc3 = One
    END IF

    NeutrinoLeptonNumber_Interior = Zero
    NeutrinoLeptonNumber_Change   = Zero

    NeutrinoEnergy_Interior = Zero
    NeutrinoEnergy_Change   = Zero

    NeutrinoMomentumX1_Interior = Zero
    NeutrinoMomentumX1_Change   = Zero

    NeutrinoMomentumX2_Interior = Zero
    NeutrinoMomentumX2_Change   = Zero

    NeutrinoMomentumX3_Interior = Zero
    NeutrinoMomentumX3_Change   = Zero

    IF( .NOT. InitializeFromCheckpoint )THEN

      NeutrinoLeptonNumber_Initial = Zero
      NeutrinoLeptonNumber_OffGrid = Zero

      NeutrinoEnergy_Initial = Zero
      NeutrinoEnergy_OffGrid = Zero

      NeutrinoMomentumX1_Initial = Zero
      NeutrinoMomentumX1_OffGrid = Zero

      NeutrinoMomentumX2_Initial = Zero
      NeutrinoMomentumX2_OffGrid = Zero

      NeutrinoMomentumX3_Initial = Zero
      NeutrinoMomentumX3_OffGrid = Zero

      NeutrinoEnergy_PL     = Zero
      NeutrinoMomentumX1_PL = Zero
      NeutrinoMomentumX2_PL = Zero
      NeutrinoMomentumX3_PL = Zero

    END IF

    FileNameRoot = TRIM( TallyFileNameRoot_TwoMoment )

    TallyChkFileNameRoot = TRIM( FileNameRoot )

    NeutrinoLeptonNumber_FileName &
      = TRIM( FileNameRoot ) // '_NeutrinoLeptonNumber.dat'
    NeutrinoEnergy_FileName &
      = TRIM( FileNameRoot ) // '_NeutrinoEnergy.dat'
    NeutrinoMomentumX1_FileName &
      = TRIM( FileNameRoot ) // '_NeutrinoMomentumX1.dat'
    NeutrinoMomentumX2_FileName &
      = TRIM( FileNameRoot ) // '_NeutrinoMomentumX2.dat'
    NeutrinoMomentumX3_FileName &
      = TRIM( FileNameRoot ) // '_NeutrinoMomentumX3.dat'

    !IF( InitializeFromCheckpoint ) RETURN

    IF( amrex_parallel_ioprocessor() )THEN

      TimeLabel = 'Time [' // TRIM( UnitsDisplay % TimeLabel ) // ']'

      CALL CreateFile( NeutrinoLeptonNumber_FileName, '', TimeLabel )

      CALL CreateFile( NeutrinoEnergy_FileName, &
                       UnitsDisplay % EnergyGlobalLabel, TimeLabel )

      CALL CreateFile( NeutrinoMomentumX1_FileName, '', TimeLabel )
      CALL CreateFile( NeutrinoMomentumX2_FileName, '', TimeLabel )
      CALL CreateFile( NeutrinoMomentumX3_FileName, '', TimeLabel )

    END IF

    IF( InitializeFromCheckpoint ) CALL ReadTallyCheckpoint_TwoMoment

  END SUBROUTINE InitializeTally_TwoMoment_MF


  SUBROUTINE FinalizeTally_TwoMoment_MF
  !Return nothing for now? Not sure, have to cvheck

  END SUBROUTINE FinalizeTally_TwoMoment_MF


  SUBROUTINE ComputeTally_TwoMoment_MF &
    ( Time, MF_uGF, MF_uCF, MF_uCR, SetInitialValues_Option, &
      WriteTally_Option, Verbose_Option )

    REAL(DP),             INTENT(in) :: Time  (0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCR(0:)
    LOGICAL,              INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL,              INTENT(in), OPTIONAL :: WriteTally_Option
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: SetInitialValues, WriteTally, Verbose

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: iZ_B0(4), iZ_E0(4)
    INTEGER :: iLo_GF(4), iLo_CF(4), iLo_CR(4)

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask
    TYPE(MeshType)        :: MeshX(3)

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    INTEGER,  CONTIGUOUS, POINTER :: FineMask(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: M(:,:,:,:,:,:,:)

    REAL(DP) :: dN, dE, dG1, dG2, dG3

    !IF( amrex_parallel_ioprocessor() ) &
    !  WRITE(*,'(A,L2)') 'ComputeTally_TwoMoment_MF: Suppress =', &
    !    SuppressTally_TwoMoment

    IF( SuppressTally_TwoMoment ) RETURN

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    WriteTally = .TRUE.
    IF( PRESENT( WriteTally_Option ) ) &
      WriteTally = WriteTally_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    NeutrinoLeptonNumber_Interior = Zero
    NeutrinoEnergy_Interior       = Zero
    NeutrinoMomentumX1_Interior   = Zero
    NeutrinoMomentumX2_Interior   = Zero
    NeutrinoMomentumX3_Interior   = Zero

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        FineMask => iMF_FineMask   % DataPtr( MFI )
        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uCF      => MF_uCF(iLevel) % DataPtr( MFI )
        uCR      => MF_uCR(iLevel) % DataPtr( MFI )

        iLo_GF = LBOUND( uGF )
        iLo_CF = LBOUND( uCF )
        iLo_CR = LBOUND( uCR )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        iZ_B0(1)   = iE_B0
        iZ_E0(1)   = iE_E0
        iZ_B0(2:4) = iX_B0
        iZ_E0(2:4) = iX_E0

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCF ], &
                 U )

        CALL AllocateArray_Z &
               ( [ 1       , iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), &
                   1       , 1        ], &
                 [ nDOFZ   , iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4), &
                   nCR     , nSpecies ], &
                 M )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_GF, iX_B0, iX_E0, uGF, G )
        CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo_CF, iX_B0, iX_E0, uCF, U )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B0, iZ_E0, iLo_CR, iZ_B0, iZ_E0, uCR, M )

        CALL ComputeTally_TwoMoment_Box &
               ( iZ_B0, iZ_E0, MeshX, G, U, M, &
                 FineMask(iX_B0(1):iX_E0(1), &
                          iX_B0(2):iX_E0(2), &
                          iX_B0(3):iX_E0(3), 1:1), &
                 dN, dE, dG1, dG2, dG3 )

        NeutrinoLeptonNumber_Interior = NeutrinoLeptonNumber_Interior + dN
        NeutrinoEnergy_Interior       = NeutrinoEnergy_Interior       + dE
        NeutrinoMomentumX1_Interior   = NeutrinoMomentumX1_Interior   + dG1
        NeutrinoMomentumX2_Interior   = NeutrinoMomentumX2_Interior   + dG2
        NeutrinoMomentumX3_Interior   = NeutrinoMomentumX3_Interior   + dG3

        CALL DeallocateArray_Z &
               ( [ 1       , iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), &
                   1       , 1        ], &
                 [ nDOFZ   , iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4), &
                   nCR     , nSpecies ], &
                 M )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCF ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

    NeutrinoLeptonNumber_Interior = FourPi * NeutrinoLeptonNumber_Interior / hc3
    NeutrinoEnergy_Interior       = FourPi * NeutrinoEnergy_Interior       / hc3
    NeutrinoMomentumX1_Interior   = FourPi * NeutrinoMomentumX1_Interior   / hc3
    NeutrinoMomentumX2_Interior   = FourPi * NeutrinoMomentumX2_Interior   / hc3
    NeutrinoMomentumX3_Interior   = FourPi * NeutrinoMomentumX3_Interior   / hc3

    CALL amrex_parallel_reduce_sum( NeutrinoLeptonNumber_Interior )
    CALL amrex_parallel_reduce_sum( NeutrinoEnergy_Interior       )
    CALL amrex_parallel_reduce_sum( NeutrinoMomentumX1_Interior   )
    CALL amrex_parallel_reduce_sum( NeutrinoMomentumX2_Interior   )
    CALL amrex_parallel_reduce_sum( NeutrinoMomentumX3_Interior   )

    IF( SetInitialValues )THEN

      NeutrinoLeptonNumber_Initial = NeutrinoLeptonNumber_Interior
      NeutrinoEnergy_Initial       = NeutrinoEnergy_Interior
      NeutrinoMomentumX1_Initial   = NeutrinoMomentumX1_Interior
      NeutrinoMomentumX2_Initial   = NeutrinoMomentumX2_Interior
      NeutrinoMomentumX3_Initial   = NeutrinoMomentumX3_Interior

    END IF

    NeutrinoLeptonNumber_Change &
      = NeutrinoLeptonNumber_Interior &
          - ( NeutrinoLeptonNumber_Initial + NeutrinoLeptonNumber_OffGrid )

    NeutrinoEnergy_Change &
      = NeutrinoEnergy_Interior &
          - ( NeutrinoEnergy_Initial       + NeutrinoEnergy_OffGrid )

    NeutrinoMomentumX1_Change &
      = NeutrinoMomentumX1_Interior &
          - ( NeutrinoMomentumX1_Initial   + NeutrinoMomentumX1_OffGrid )

    NeutrinoMomentumX2_Change &
      = NeutrinoMomentumX2_Interior &
          - ( NeutrinoMomentumX2_Initial   + NeutrinoMomentumX2_OffGrid )

    NeutrinoMomentumX3_Change &
      = NeutrinoMomentumX3_Interior &
          - ( NeutrinoMomentumX3_Initial   + NeutrinoMomentumX3_OffGrid )

    IF( WriteTally ) CALL WriteTally_TwoMoment( Time(0) )

    IF( Verbose ) CALL DisplayTally( Time(0) )

  END SUBROUTINE ComputeTally_TwoMoment_MF


  SUBROUTINE ComputeTally_TwoMoment_Box &
    ( iZ_B0, iZ_E0, MeshX, G, U, M, FineMask, N_Box, E_Box, G1_Box, G2_Box, G3_Box )

    INTEGER,        INTENT(in)  :: iZ_B0(4), iZ_E0(4)
    TYPE(MeshType), INTENT(in)  :: MeshX(3)
    REAL(DP),       INTENT(in)  :: G(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(DP),       INTENT(in)  :: U(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(DP),       INTENT(in)  :: M(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):, &
                                     iZ_B0(4):,1:,1:)
    INTEGER,        INTENT(in)  :: FineMask(iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(DP),       INTENT(out) :: N_Box, E_Box, G1_Box, G2_Box, G3_Box

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iNodeE, iNodeX, iNodeZ
    REAL(DP) :: d4Z, W

    REAL(DP) :: P(1:nDOFX, &
                  iZ_B0(2):iZ_E0(2), &
                  iZ_B0(3):iZ_E0(3), &
                  iZ_B0(4):iZ_E0(4), &
                  1:nPF)

    N_Box  = Zero
    E_Box  = Zero
    G1_Box = Zero
    G2_Box = Zero
    G3_Box = Zero

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      IF( IsNotLeafElement( FineMask(iZ2,iZ3,iZ4,1) ) ) CYCLE

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( U(iNodeX,iZ2,iZ3,iZ4,iCF_D ),       &
                 U(iNodeX,iZ2,iZ3,iZ4,iCF_S1),       &
                 U(iNodeX,iZ2,iZ3,iZ4,iCF_S2),       &
                 U(iNodeX,iZ2,iZ3,iZ4,iCF_S3),       &
                 U(iNodeX,iZ2,iZ3,iZ4,iCF_E ),       &
                 U(iNodeX,iZ2,iZ3,iZ4,iCF_Ne),       &
                 P(iNodeX,iZ2,iZ3,iZ4,iPF_D ),       &
                 P(iNodeX,iZ2,iZ3,iZ4,iPF_V1),       &
                 P(iNodeX,iZ2,iZ3,iZ4,iPF_V2),       &
                 P(iNodeX,iZ2,iZ3,iZ4,iPF_V3),       &
                 P(iNodeX,iZ2,iZ3,iZ4,iPF_E ),       &
                 P(iNodeX,iZ2,iZ3,iZ4,iPF_Ne),       &
                 G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      IF( IsNotLeafElement( FineMask(iZ2,iZ3,iZ4,1) ) ) CYCLE

      DO iZ1 = iZ_B0(1), iZ_E0(1)

        d4Z = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4)

        DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE

          iNodeZ = ( iNodeX - 1 ) * nDOFE + iNodeE

          W = d4Z * Weights_q(iNodeZ) * G(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

          N_Box                                            &
            = N_Box                                        &
                + W * uGE(iNodeE,iZ1,iGE_Ep2)              &
                    * LeptonNumber(iS)                     &
                    * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS)

          E_Box                                            &
            = E_Box                                        &
                + W * uGE(iNodeE,iZ1,iGE_Ep3)              &
                    * ( M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)    &
                        + P(iNodeX,iZ2,iZ3,iZ4,iPF_V1)         &
                            * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) &
                        + P(iNodeX,iZ2,iZ3,iZ4,iPF_V2)         &
                            * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) &
                        + P(iNodeX,iZ2,iZ3,iZ4,iPF_V3)         &
                            * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) )

          G1_Box                                           &
            = G1_Box                                       &
                + W * uGE(iNodeE,iZ1,iGE_Ep3)              &
                    * ( M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)   &
                        + G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)  &
                            * P(iNodeX,iZ2,iZ3,iZ4,iPF_V1)    &
                            * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )

          G2_Box                                           &
            = G2_Box                                       &
                + W * uGE(iNodeE,iZ1,iGE_Ep3)              &
                    * ( M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)   &
                        + G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)  &
                            * P(iNodeX,iZ2,iZ3,iZ4,iPF_V2)    &
                            * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )

          G3_Box                                           &
            = G3_Box                                       &
                + W * uGE(iNodeE,iZ1,iGE_Ep3)              &
                    * ( M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)   &
                        + G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)  &
                            * P(iNodeX,iZ2,iZ3,iZ4,iPF_V3)    &
                            * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )

        END DO
        END DO

      END DO

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dZ1, bla bla.

  END SUBROUTINE ComputeTally_TwoMoment_Box


  SUBROUTINE IncrementOffGridTally_TwoMoment_MF( dM )

    REAL(DP), INTENT(in) :: dM(1:,0:)

    INTEGER :: iLevel

    IF( SuppressTally_TwoMoment ) RETURN

    DO iLevel = 0, nLevels-1

      NeutrinoLeptonNumber_OffGrid &
        = NeutrinoLeptonNumber_OffGrid + FourPi * dM(iCR_N     ,iLevel) / hc3

      NeutrinoEnergy_OffGrid &
        = NeutrinoEnergy_OffGrid       + FourPi * dM(nCR+iCR_N ,iLevel) / hc3

      NeutrinoMomentumX1_OffGrid &
        = NeutrinoMomentumX1_OffGrid   + FourPi * dM(nCR+iCR_G1,iLevel) / hc3

      NeutrinoMomentumX2_OffGrid &
        = NeutrinoMomentumX2_OffGrid   + FourPi * dM(nCR+iCR_G2,iLevel) / hc3

      NeutrinoMomentumX3_OffGrid &
        = NeutrinoMomentumX3_OffGrid   + FourPi * dM(nCR+iCR_G3,iLevel) / hc3

    END DO

  END SUBROUTINE IncrementOffGridTally_TwoMoment_MF


  SUBROUTINE IncrementPositivityLimiterTally_TwoMoment_MF( dM )


    REAL(DP), INTENT(in) :: dM(1:,0:)

    INTEGER :: iLevel

    IF( SuppressTally_TwoMoment ) RETURN

    DO iLevel = 0, nLevels-1

      NeutrinoEnergy_PL     = NeutrinoEnergy_PL     + dM(iCR_N ,iLevel)
      NeutrinoMomentumX1_PL = NeutrinoMomentumX1_PL + dM(iCR_G1,iLevel)
      NeutrinoMomentumX2_PL = NeutrinoMomentumX2_PL + dM(iCR_G2,iLevel)
      NeutrinoMomentumX3_PL = NeutrinoMomentumX3_PL + dM(iCR_G3,iLevel)

    END DO

  END SUBROUTINE IncrementPositivityLimiterTally_TwoMoment_MF


  SUBROUTINE WriteTally_TwoMoment( Time )

    REAL(DP), INTENT(in) :: Time

    IF( .NOT. amrex_parallel_ioprocessor() ) RETURN

    CALL WriteTallyToFile &
           ( NeutrinoLeptonNumber_FileName, Time, UnitsDisplay % TimeUnit, &
             NeutrinoLeptonNumber_Interior, NeutrinoLeptonNumber_Initial,  &
             NeutrinoLeptonNumber_OffGrid , NeutrinoLeptonNumber_Change,   &
             One )

    CALL WriteTallyToFile &
           ( NeutrinoEnergy_FileName, Time, UnitsDisplay % TimeUnit, &
             NeutrinoEnergy_Interior, NeutrinoEnergy_Initial,        &
             NeutrinoEnergy_OffGrid , NeutrinoEnergy_Change,         &
             UnitsDisplay % EnergyGlobalUnit )

    CALL WriteTallyToFile &
           ( NeutrinoMomentumX1_FileName, Time, UnitsDisplay % TimeUnit, &
             NeutrinoMomentumX1_Interior, NeutrinoMomentumX1_Initial,    &
             NeutrinoMomentumX1_OffGrid , NeutrinoMomentumX1_Change,     &
             One )

    CALL WriteTallyToFile &
           ( NeutrinoMomentumX2_FileName, Time, UnitsDisplay % TimeUnit, &
             NeutrinoMomentumX2_Interior, NeutrinoMomentumX2_Initial,    &
             NeutrinoMomentumX2_OffGrid , NeutrinoMomentumX2_Change,     &
             One )

    CALL WriteTallyToFile &
           ( NeutrinoMomentumX3_FileName, Time, UnitsDisplay % TimeUnit, &
             NeutrinoMomentumX3_Interior, NeutrinoMomentumX3_Initial,    &
             NeutrinoMomentumX3_OffGrid , NeutrinoMomentumX3_Change,     &
             One )

  END SUBROUTINE WriteTally_TwoMoment


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    IF( .NOT. amrex_parallel_ioprocessor() ) RETURN

    IF( NeutrinoEnergy_Interior .NE. Zero ) &
      WRITE(*,'(6x,A40,2ES15.7E3)') 'Energy  Change | PL  (/Interior).: ', &
        NeutrinoEnergy_Change / NeutrinoEnergy_Interior, &
        NeutrinoEnergy_PL     / NeutrinoEnergy_Interior

    WRITE(*,*)
    WRITE(*,'(6x,A,ES13.6E3,x,A)') &
      'TwoMoment Tally. t = ', &
      Time / UnitsDisplay % TimeUnit, TRIM( UnitsDisplay % TimeLabel )

    CALL WriteTallyToScreen &
           ( 'Neutrino Lepton Number', &
             NeutrinoLeptonNumber_Interior, NeutrinoLeptonNumber_Initial, &
             NeutrinoLeptonNumber_OffGrid , NeutrinoLeptonNumber_Change,  &
             One, '' )

    CALL WriteTallyToScreen &
           ( 'Neutrino Energy', &
             NeutrinoEnergy_Interior, NeutrinoEnergy_Initial, &
             NeutrinoEnergy_OffGrid , NeutrinoEnergy_Change,  &
             UnitsDisplay % EnergyGlobalUnit, &
             TRIM( UnitsDisplay % EnergyGlobalLabel ) )

    WRITE(*,*)

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

    InteriorLabel = 'Interior [' // TRIM( UnitsLabel ) // ']'
    OffGridLabel  = 'Off Grid [' // TRIM( UnitsLabel ) // ']'
    InitialLabel  = 'Initial ['  // TRIM( UnitsLabel ) // ']'
    ChangeLabel   = 'Change ['   // TRIM( UnitsLabel ) // ']'

    CALL CheckFileExistenceAndAppend( FileName )

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ) )

    WRITE( FileUnit, '(5(A25,x))' ) &
      TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
      TRIM( InitialLabel ), TRIM( ChangeLabel )

    CLOSE( FileUnit )

  END SUBROUTINE CreateFile

  SUBROUTINE WriteTallyToScreen &
    ( FieldName, Interior, Initial, OffGrid, Change, Units, Label )

    CHARACTER(*), INTENT(in) :: FieldName, Label
    REAL(DP),     INTENT(in) :: Interior, Initial, OffGrid, Change, Units

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
    ( FileName, Time, TimeUnit, Interior, Initial, OffGrid, Change, Units )

    CHARACTER(*), INTENT(in) :: FileName
    REAL(DP),     INTENT(in) :: Time, TimeUnit, &
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

  SUBROUTINE WriteTallyCheckpoint_TwoMoment( StepNumber )

    INTEGER, INTENT(in) :: StepNumber

    INTEGER       :: FileUnit
    CHARACTER(SL) :: FileName

    IF( SuppressTally_TwoMoment ) RETURN
    IF( .NOT. amrex_parallel_ioprocessor() ) RETURN

    WRITE( FileName, '(A,A,I8.8,A)' ) &
      TRIM( TallyChkFileNameRoot ), '_TwoMomentTally_', StepNumber, '.dat'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ), ACTION = 'WRITE' )

    WRITE( FileUnit, '(ES24.16)' ) NeutrinoLeptonNumber_Initial
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoLeptonNumber_OffGrid
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoEnergy_Initial
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoEnergy_OffGrid
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX1_Initial
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX1_OffGrid
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX2_Initial
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX2_OffGrid
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX3_Initial
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX3_OffGrid
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoEnergy_PL
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX1_PL
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX2_PL
    WRITE( FileUnit, '(ES24.16)' ) NeutrinoMomentumX3_PL

    CLOSE( FileUnit )

  END SUBROUTINE WriteTallyCheckpoint_TwoMoment


  SUBROUTINE ReadTallyCheckpoint_TwoMoment

    INTEGER       :: FileUnit
    LOGICAL       :: IsFile
    CHARACTER(SL) :: FileName

    IF( SuppressTally_TwoMoment ) RETURN

    WRITE( FileName, '(A,A,I8.8,A)' ) &
      TRIM( TallyChkFileNameRoot ), '_TwoMomentTally_', iRestart, '.dat'

    INQUIRE( FILE = TRIM( FileName ), EXIST = IsFile )

    IF( .NOT. IsFile )THEN
      IF( amrex_parallel_ioprocessor() ) &
        WRITE(*,'(A)') &
          '  WARNING: TwoMoment tally checkpoint not found: ' &
          // TRIM( FileName ) // ' -- Initial/OffGrid start at zero'
      RETURN
    END IF

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ), ACTION = 'READ' )

    READ( FileUnit, * ) NeutrinoLeptonNumber_Initial
    READ( FileUnit, * ) NeutrinoLeptonNumber_OffGrid
    READ( FileUnit, * ) NeutrinoEnergy_Initial
    READ( FileUnit, * ) NeutrinoEnergy_OffGrid
    READ( FileUnit, * ) NeutrinoMomentumX1_Initial
    READ( FileUnit, * ) NeutrinoMomentumX1_OffGrid
    READ( FileUnit, * ) NeutrinoMomentumX2_Initial
    READ( FileUnit, * ) NeutrinoMomentumX2_OffGrid
    READ( FileUnit, * ) NeutrinoMomentumX3_Initial
    READ( FileUnit, * ) NeutrinoMomentumX3_OffGrid
    READ( FileUnit, * ) NeutrinoEnergy_PL
    READ( FileUnit, * ) NeutrinoMomentumX1_PL
    READ( FileUnit, * ) NeutrinoMomentumX2_PL
    READ( FileUnit, * ) NeutrinoMomentumX3_PL

    CLOSE( FileUnit )

    IF( amrex_parallel_ioprocessor() ) &
      WRITE(*,'(A)') '  Restored TwoMoment tally from ' // TRIM( FileName )

  END SUBROUTINE ReadTallyCheckpoint_TwoMoment

END MODULE MF_TwoMoment_TallyModule