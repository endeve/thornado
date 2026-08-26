MODULE MF_TwoMoment_TallyModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE UnitsModule, ONLY: &
    UnitsActive, &
    SpeedOfLight, &
    PlanckConstant, &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    swX, swE, nDOFX, nDOFE, nDOFZ, nNodesE, &
    iE_B0, iE_E0
  USE ReferenceElementModule, ONLY: &
    Weights_q
  USE MeshModule, ONLY: &
    MeshType, &
    CreateMesh, &
    DestroyMesh
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_SqrtGm
  USE GeometryFieldsModuleE, ONLY: &
    uGE, iGE_Ep2, iGE_Ep3
  USE RadiationFieldsModule, ONLY: &
    nSpecies, LeptonNumber, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, Zero, One, FourPi
  USE InputParsingModule, ONLY: &
    nLevels, &
    ProgramName, &
    eL, eR, nE, zoomE, &
    UseTiling
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_Z, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_TwoMoment_MF
  PUBLIC :: ComputeTally_TwoMoment_MF
  PUBLIC :: FinalizeTally_TwoMoment_MF
  PUBLIC :: IncrementOffGridTally_TwoMoment_MF
  PUBLIC :: IncrementPositivityLimiterTally_TwoMoment_MF


  LOGICAL, PARAMETER :: TallyLevel0Only = .TRUE.

  LOGICAL  :: SuppressTally
  REAL(DP) :: hc3

  CHARACTER(256) :: NeutrinoLeptonNumber_FileName
  REAL(DP), ALLOCATABLE :: NeutrinoLeptonNumber_Interior(:)
  REAL(DP), ALLOCATABLE :: NeutrinoLeptonNumber_Initial (:)
  REAL(DP), ALLOCATABLE :: NeutrinoLeptonNumber_OffGrid (:)
  REAL(DP), ALLOCATABLE :: NeutrinoLeptonNumber_Change  (:)

  CHARACTER(256) :: NeutrinoEnergy_FileName
  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_Interior(:)
  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_Initial (:)
  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_OffGrid (:)
  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_Change  (:)

  CHARACTER(256) :: Momentum_FileName
  REAL(DP), ALLOCATABLE :: Momentum_X1(:)
  REAL(DP), ALLOCATABLE :: Momentum_X2(:)
  REAL(DP), ALLOCATABLE :: Momentum_X3(:)

  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_PL   (:)
  REAL(DP), ALLOCATABLE :: NeutrinoMomentum1_PL(:)
  REAL(DP), ALLOCATABLE :: NeutrinoMomentum2_PL(:)
  REAL(DP), ALLOCATABLE :: NeutrinoMomentum3_PL(:)

CONTAINS


  SUBROUTINE InitializeTally_TwoMoment_MF &
    ( SuppressTally_Option, BaseFileName_Option )

    LOGICAL,          INTENT(in), OPTIONAL :: SuppressTally_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: BaseFileName_Option

    CHARACTER(256) :: BaseFileName
    INTEGER        :: FileUnit
    CHARACTER(256) :: TimeLabel
    CHARACTER(256) :: LeptonNumber_InteriorLabel, LeptonNumber_OffGridLabel
    CHARACTER(256) :: LeptonNumber_InitialLabel,  LeptonNumber_ChangeLabel
    CHARACTER(256) :: Energy_InteriorLabel, Energy_OffGridLabel
    CHARACTER(256) :: Energy_InitialLabel,  Energy_ChangeLabel
    CHARACTER(256) :: Momentum1Label, Momentum2Label, Momentum3Label

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

    IF( SuppressTally ) RETURN

    IF( UnitsActive )THEN
      hc3 = ( PlanckConstant * SpeedOfLight )**3
    ELSE
      hc3 = One
    END IF

    ALLOCATE( NeutrinoLeptonNumber_Interior(0:nLevels-1) )
    ALLOCATE( NeutrinoLeptonNumber_Initial (0:nLevels-1) )
    ALLOCATE( NeutrinoLeptonNumber_OffGrid (0:nLevels-1) )
    ALLOCATE( NeutrinoLeptonNumber_Change  (0:nLevels-1) )

    ALLOCATE( NeutrinoEnergy_Interior(0:nLevels-1) )
    ALLOCATE( NeutrinoEnergy_Initial (0:nLevels-1) )
    ALLOCATE( NeutrinoEnergy_OffGrid (0:nLevels-1) )
    ALLOCATE( NeutrinoEnergy_Change  (0:nLevels-1) )

    ALLOCATE( Momentum_X1(0:nLevels-1) )
    ALLOCATE( Momentum_X2(0:nLevels-1) )
    ALLOCATE( Momentum_X3(0:nLevels-1) )

    ALLOCATE( NeutrinoEnergy_PL   (0:nLevels-1) )
    ALLOCATE( NeutrinoMomentum1_PL(0:nLevels-1) )
    ALLOCATE( NeutrinoMomentum2_PL(0:nLevels-1) )
    ALLOCATE( NeutrinoMomentum3_PL(0:nLevels-1) )

    NeutrinoLeptonNumber_Interior = Zero
    NeutrinoLeptonNumber_Initial  = Zero
    NeutrinoLeptonNumber_OffGrid  = Zero
    NeutrinoLeptonNumber_Change   = Zero

    NeutrinoEnergy_Interior = Zero
    NeutrinoEnergy_Initial  = Zero
    NeutrinoEnergy_OffGrid  = Zero
    NeutrinoEnergy_Change   = Zero

    Momentum_X1 = Zero
    Momentum_X2 = Zero
    Momentum_X3 = Zero

    NeutrinoEnergy_PL    = Zero
    NeutrinoMomentum1_PL = Zero
    NeutrinoMomentum2_PL = Zero
    NeutrinoMomentum3_PL = Zero

    IF( amrex_parallel_ioprocessor() )THEN

      BaseFileName = ''
      IF( PRESENT( BaseFileName_Option ) ) &
        BaseFileName = TRIM( BaseFileName_Option )

      BaseFileName = TRIM( BaseFileName ) // TRIM( ProgramName )

      TimeLabel = 'Time [' // TRIM( UnitsDisplay % TimeLabel ) // ']'

      ! --- Neutrino Lepton Number ---

      NeutrinoLeptonNumber_FileName &
        = TRIM( BaseFileName ) // '.Tally_NeutrinoLeptonNumber.dat'

      LeptonNumber_InteriorLabel = 'LeptonNumber_Interior'
      LeptonNumber_OffGridLabel  = 'LeptonNumber_OffGrid'
      LeptonNumber_InitialLabel  = 'LeptonNumber_Initial'
      LeptonNumber_ChangeLabel   = 'LeptonNumber_Change'

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( NeutrinoLeptonNumber_FileName ) )
      WRITE( FileUnit, '(5(A25,x))' ) &
        TRIM( TimeLabel ), TRIM( LeptonNumber_InteriorLabel ), &
        TRIM( LeptonNumber_OffGridLabel ), TRIM( LeptonNumber_InitialLabel ), &
        TRIM( LeptonNumber_ChangeLabel )
      CLOSE( FileUnit )

      ! --- Neutrino Energy ---

      NeutrinoEnergy_FileName &
        = TRIM( BaseFileName ) // '.Tally_NeutrinoEnergy.dat'

      Energy_InteriorLabel &
        = 'Energy_Interior [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      Energy_OffGridLabel &
        = 'Energy_OffGrid ['  // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      Energy_InitialLabel &
        = 'Energy_Initial ['  // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      Energy_ChangeLabel &
        = 'Energy_Change ['   // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( NeutrinoEnergy_FileName ) )
      WRITE( FileUnit, '(5(A25,x))' ) &
        TRIM( TimeLabel ), TRIM( Energy_InteriorLabel ), &
        TRIM( Energy_OffGridLabel ), TRIM( Energy_InitialLabel ), &
        TRIM( Energy_ChangeLabel )
      CLOSE( FileUnit )

      ! --- Momentum ---

      Momentum_FileName = TRIM( BaseFileName ) // '.Tally_Momentum.dat'

      Momentum1Label = 'Momentum_1'
      Momentum2Label = 'Momentum_2'
      Momentum3Label = 'Momentum_3'

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( Momentum_FileName ) )
      WRITE( FileUnit, '(4(A25,x))' ) &
        TRIM( TimeLabel ), TRIM( Momentum1Label ), &
        TRIM( Momentum2Label ), TRIM( Momentum3Label )
      CLOSE( FileUnit )

    END IF

  END SUBROUTINE InitializeTally_TwoMoment_MF


  SUBROUTINE FinalizeTally_TwoMoment_MF

    IF( SuppressTally ) RETURN

    DEALLOCATE( NeutrinoLeptonNumber_Interior )
    DEALLOCATE( NeutrinoLeptonNumber_Initial  )
    DEALLOCATE( NeutrinoLeptonNumber_OffGrid  )
    DEALLOCATE( NeutrinoLeptonNumber_Change   )

    DEALLOCATE( NeutrinoEnergy_Interior )
    DEALLOCATE( NeutrinoEnergy_Initial  )
    DEALLOCATE( NeutrinoEnergy_OffGrid  )
    DEALLOCATE( NeutrinoEnergy_Change   )

    DEALLOCATE( Momentum_X1 )
    DEALLOCATE( Momentum_X2 )
    DEALLOCATE( Momentum_X3 )

    DEALLOCATE( NeutrinoEnergy_PL    )
    DEALLOCATE( NeutrinoMomentum1_PL )
    DEALLOCATE( NeutrinoMomentum2_PL )
    DEALLOCATE( NeutrinoMomentum3_PL )

  END SUBROUTINE FinalizeTally_TwoMoment_MF


  SUBROUTINE ComputeTally_TwoMoment_MF &
    ( MF_uGF, MF_uCR, Time, SetInitialValues_Option, Verbose_Option )

    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCR(0:nLevels-1)
    REAL(DP),             INTENT(in) :: Time
    LOGICAL,              INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL  :: SetInitialValues, Verbose
    INTEGER  :: iLevel
    INTEGER  :: iX_B0(3), iX_E0(3), iZ_B0(4), iZ_E0(4)
    INTEGER  :: iLo_GF(4), iLo_CR(4)
    REAL(DP) :: dN, dE, dG1, dG2, dG3

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(MeshType)     :: MeshX(3), MeshE

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U(:,:,:,:,:,:,:)

    IF( SuppressTally ) RETURN

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL CreateMesh( MeshE, nE, nNodesE, swE, eL, eR, zoomE )

    DO iLevel = 0, nLevels-1

      NeutrinoLeptonNumber_Interior(iLevel) = Zero
      NeutrinoEnergy_Interior      (iLevel) = Zero
      Momentum_X1                  (iLevel) = Zero
      Momentum_X2                  (iLevel) = Zero
      Momentum_X3                  (iLevel) = Zero

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )

        iLo_GF = LBOUND( uGF )
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

        CALL AllocateArray_Z &
               ( [ 1       , iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), &
                   1       , 1        ], &
                 [ nDOFZ   , iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4), &
                   nCR     , nSpecies ], &
                 U )

        CALL amrex2thornado_X &
               ( nGF, iX_B0, iX_E0, iLo_GF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B0, iZ_E0, iLo_CR, iZ_B0, iZ_E0, uCR, U )

        CALL ComputeTally_Interior &
               ( iZ_B0, iZ_E0, G, U, MeshX, MeshE, dN, dE, dG1, dG2, dG3 )

        NeutrinoLeptonNumber_Interior(iLevel) &
          = NeutrinoLeptonNumber_Interior(iLevel) + dN
        NeutrinoEnergy_Interior(iLevel) &
          = NeutrinoEnergy_Interior(iLevel) + dE
        Momentum_X1(iLevel) = Momentum_X1(iLevel) + dG1
        Momentum_X2(iLevel) = Momentum_X2(iLevel) + dG2
        Momentum_X3(iLevel) = Momentum_X3(iLevel) + dG3

        CALL DeallocateArray_Z &
               ( [ 1       , iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), &
                   1       , 1        ], &
                 [ nDOFZ   , iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4), &
                   nCR     , nSpecies ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO

    CALL DestroyMesh( MeshE )

    CALL amrex_parallel_reduce_sum( NeutrinoLeptonNumber_Interior, nLevels )
    CALL amrex_parallel_reduce_sum( NeutrinoEnergy_Interior      , nLevels )
    CALL amrex_parallel_reduce_sum( Momentum_X1                  , nLevels )
    CALL amrex_parallel_reduce_sum( Momentum_X2                  , nLevels )
    CALL amrex_parallel_reduce_sum( Momentum_X3                  , nLevels )

    IF( SetInitialValues )THEN

      DO iLevel = 0, nLevels-1
        NeutrinoLeptonNumber_Initial(iLevel) &
          = NeutrinoLeptonNumber_Interior(iLevel)
        NeutrinoEnergy_Initial(iLevel) &
          = NeutrinoEnergy_Interior(iLevel)
      END DO

    END IF

    DO iLevel = 0, nLevels-1

      NeutrinoLeptonNumber_Change(iLevel) &
        = NeutrinoLeptonNumber_Interior(iLevel) &
            - (   NeutrinoLeptonNumber_Initial(iLevel) &
                + NeutrinoLeptonNumber_OffGrid(iLevel) )

      NeutrinoEnergy_Change(iLevel) &
        = NeutrinoEnergy_Interior(iLevel) &
            - (   NeutrinoEnergy_Initial(iLevel) &
                + NeutrinoEnergy_OffGrid(iLevel) )

    END DO

    CALL WriteTally_TwoMoment( Time )

    IF( Verbose ) CALL DisplayTally( Time )

  END SUBROUTINE ComputeTally_TwoMoment_MF


  SUBROUTINE ComputeTally_Interior &
    ( iZ_B0, iZ_E0, G, U, MeshX, MeshE, N_Int, E_Int, G1_Int, G2_Int, G3_Int )

    INTEGER,        INTENT(in)  :: iZ_B0(4), iZ_E0(4)
    REAL(DP),       INTENT(in)  :: G(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(DP),       INTENT(in)  :: U(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):, &
                                     iZ_B0(4):,1:,1:)
    TYPE(MeshType), INTENT(in)  :: MeshX(3), MeshE
    REAL(DP),       INTENT(out) :: N_Int, E_Int, G1_Int, G2_Int, G3_Int

    INTEGER  :: iS, iZ1, iZ2, iZ3, iZ4, iNodeZ, iNodeE, iNodeX
    REAL(DP) :: d4Z, W, Wn, We

    N_Int  = Zero
    E_Int  = Zero
    G1_Int = Zero
    G2_Int = Zero
    G3_Int = Zero

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      d4Z =   MeshE    % Width(iZ1) * MeshX(1) % Width(iZ2) &
            * MeshX(2) % Width(iZ3) * MeshX(3) % Width(iZ4)

      DO iNodeZ = 1, nDOFZ

        iNodeE = MOD( ( iNodeZ - 1 )         , nDOFE ) + 1
        iNodeX = MOD( ( iNodeZ - 1 ) / nDOFE , nDOFX ) + 1

        W  = Weights_q(iNodeZ) * d4Z * G(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

        Wn = W * uGE(iNodeE,iZ1,iGE_Ep2)
        We = W * uGE(iNodeE,iZ1,iGE_Ep3)

        N_Int = N_Int + LeptonNumber(iS) * Wn &
                          * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)
        E_Int = E_Int + We * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)

        G1_Int = G1_Int + We * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)
        G2_Int = G2_Int + We * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)
        G3_Int = G3_Int + We * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    N_Int  = FourPi * N_Int  / hc3
    E_Int  = FourPi * E_Int  / hc3
    G1_Int = FourPi * G1_Int / hc3
    G2_Int = FourPi * G2_Int / hc3
    G3_Int = FourPi * G3_Int / hc3

  END SUBROUTINE ComputeTally_Interior


  SUBROUTINE IncrementOffGridTally_TwoMoment_MF( dM )

    REAL(DP), INTENT(in) :: dM(1:,0:)

    INTEGER :: iLevel

    IF( SuppressTally ) RETURN

    DO iLevel = 0, nLevels-1

      NeutrinoLeptonNumber_OffGrid(iLevel) &
        = NeutrinoLeptonNumber_OffGrid(iLevel) &
            + FourPi * dM(iCR_N,iLevel) / hc3

      NeutrinoEnergy_OffGrid(iLevel) &
        = NeutrinoEnergy_OffGrid(iLevel) &
            + FourPi * dM(nCR+iCR_N,iLevel) / hc3

    END DO

  END SUBROUTINE IncrementOffGridTally_TwoMoment_MF


  SUBROUTINE IncrementPositivityLimiterTally_TwoMoment_MF( dM )

    REAL(DP), INTENT(in) :: dM(1:,0:)

    INTEGER :: iLevel

    IF( SuppressTally ) RETURN

    DO iLevel = 0, nLevels-1

      NeutrinoEnergy_PL   (iLevel) = NeutrinoEnergy_PL   (iLevel) + dM(iCR_N ,iLevel)
      NeutrinoMomentum1_PL(iLevel) = NeutrinoMomentum1_PL(iLevel) + dM(iCR_G1,iLevel)
      NeutrinoMomentum2_PL(iLevel) = NeutrinoMomentum2_PL(iLevel) + dM(iCR_G2,iLevel)
      NeutrinoMomentum3_PL(iLevel) = NeutrinoMomentum3_PL(iLevel) + dM(iCR_G3,iLevel)

    END DO

  END SUBROUTINE IncrementPositivityLimiterTally_TwoMoment_MF


  SUBROUTINE WriteTally_TwoMoment( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER  :: FileUnit
    REAL(DP) :: N_Int, N_Off, N_Ini, N_Chg
    REAL(DP) :: E_Int, E_Off, E_Ini, E_Chg
    REAL(DP) :: P1, P2, P3

    IF( .NOT. amrex_parallel_ioprocessor() ) RETURN

    IF( TallyLevel0Only )THEN
      N_Int = NeutrinoLeptonNumber_Interior(0)
      N_Off = NeutrinoLeptonNumber_OffGrid (0)
      N_Ini = NeutrinoLeptonNumber_Initial (0)
      N_Chg = NeutrinoLeptonNumber_Change  (0)
      E_Int = NeutrinoEnergy_Interior(0)
      E_Off = NeutrinoEnergy_OffGrid (0)
      E_Ini = NeutrinoEnergy_Initial (0)
      E_Chg = NeutrinoEnergy_Change  (0)
      P1    = Momentum_X1(0)
      P2    = Momentum_X2(0)
      P3    = Momentum_X3(0)
    ELSE
      N_Int = SUM( NeutrinoLeptonNumber_Interior )
      N_Off = SUM( NeutrinoLeptonNumber_OffGrid  )
      N_Ini = SUM( NeutrinoLeptonNumber_Initial  )
      N_Chg = SUM( NeutrinoLeptonNumber_Change   )
      E_Int = SUM( NeutrinoEnergy_Interior )
      E_Off = SUM( NeutrinoEnergy_OffGrid  )
      E_Ini = SUM( NeutrinoEnergy_Initial  )
      E_Chg = SUM( NeutrinoEnergy_Change   )
      P1    = SUM( Momentum_X1 )
      P2    = SUM( Momentum_X2 )
      P3    = SUM( Momentum_X3 )
    END IF

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( NeutrinoLeptonNumber_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )
    WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, N_Int, N_Off, N_Ini, N_Chg
    CLOSE( FileUnit )

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( NeutrinoEnergy_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )
    WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, &
      E_Int / UnitsDisplay % EnergyGlobalUnit, &
      E_Off / UnitsDisplay % EnergyGlobalUnit, &
      E_Ini / UnitsDisplay % EnergyGlobalUnit, &
      E_Chg / UnitsDisplay % EnergyGlobalUnit
    CLOSE( FileUnit )

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( Momentum_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )
    WRITE( FileUnit, '(4(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, P1, P2, P3
    CLOSE( FileUnit )

  END SUBROUTINE WriteTally_TwoMoment


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    IF( .NOT. amrex_parallel_ioprocessor() ) RETURN

    IF( NeutrinoEnergy_Interior .NE. Zero ) &
      WRITE(*,'(6x,A40,2ES15.7E3)') 'Energy  Change | PL  (/Interior).: ', &
        NeutrinoEnergy_Change / NeutrinoEnergy_Interior, &
        NeutrinoEnergy_PL     / NeutrinoEnergy_Interior

    WRITE(*,*)
    WRITE(*,'(A8,A,ES8.2E2,x,A)') &
      '', 'TwoMoment Tally. t = ', &
      Time / UnitsDisplay % TimeUnit, TRIM( UnitsDisplay % TimeLabel )
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES14.7E2)') &
      '', 'Neutrino Lepton Number.: ', NeutrinoLeptonNumber_Interior(0)
    WRITE(*,'(A6,A40,ES14.7E2)') &
      '', 'Neutrino Lepton Number Change.: ', NeutrinoLeptonNumber_Change(0)
    WRITE(*,'(A6,A40,ES14.7E2)') &
      '', 'Neutrino Energy.: ', &
      NeutrinoEnergy_Interior(0) / UnitsDisplay % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES14.7E2)') &
      '', 'Neutrino Energy Change.: ', &
      NeutrinoEnergy_Change(0) / UnitsDisplay % EnergyGlobalUnit
    WRITE(*,*)

  END SUBROUTINE DisplayTally


END MODULE MF_TwoMoment_TallyModule
