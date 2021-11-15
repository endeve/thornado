MODULE Euler_TallyModule_NonRelativistic

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDimsX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    iGF_SqrtGm, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Phi_N
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE UnitsModule, ONLY: &
    UnitsDisplay

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_Euler_NonRelativistic
  PUBLIC :: ComputeTally_Euler_NonRelativistic
  PUBLIC :: IncrementOffGridTally_Euler_NonRelativistic
  PUBLIC :: FinalizeTally_Euler_NonRelativistic

  LOGICAL :: SuppressTally

  CHARACTER(256) :: BaryonicMass_FileName
  REAL(DP)       :: BaryonicMass_Interior
  REAL(DP)       :: BaryonicMass_Initial
  REAL(DP)       :: BaryonicMass_OffGrid
  REAL(DP)       :: BaryonicMass_Change

  CHARACTER(256) :: Energy_FileName
  REAL(DP)       :: Energy_Interior
  REAL(DP)       :: Energy_Initial
  REAL(DP)       :: Energy_OffGrid
  REAL(DP)       :: Energy_Change

  CHARACTER(256) :: InternalEnergy_FileName
  REAL(DP)       :: InternalEnergy_Interior
  REAL(DP)       :: InternalEnergy_Initial
  REAL(DP)       :: InternalEnergy_OffGrid
  REAL(DP)       :: InternalEnergy_Change

  CHARACTER(256) :: KineticEnergy_FileName
  REAL(DP)       :: KineticEnergy_Interior
  REAL(DP)       :: KineticEnergy_Initial
  REAL(DP)       :: KineticEnergy_OffGrid
  REAL(DP)       :: KineticEnergy_Change

  CHARACTER(256) :: GravitationalEnergy_FileName
  REAL(DP)       :: GravitationalEnergy_Interior
  REAL(DP)       :: GravitationalEnergy_Initial
  REAL(DP)       :: GravitationalEnergy_OffGrid
  REAL(DP)       :: GravitationalEnergy_Change


CONTAINS


  SUBROUTINE InitializeTally_Euler_NonRelativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, SuppressTally_Option, &
      BaseFileName_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)           :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressTally_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      BaseFileName_Option

    CHARACTER(256) :: BaseFileName
    INTEGER        :: FileUnit

    CHARACTER(256) :: TimeLabel
    CHARACTER(256) :: InteriorLabel, InitialLabel, OffGridLabel, ChangeLabel

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

    IF( SuppressTally ) RETURN

    BaseFileName = '../Output/'
    IF( PRESENT( BaseFileName_Option ) ) &
      BaseFileName = TRIM( BaseFileName_Option )

    BaseFileName = TRIM( BaseFileName ) // TRIM( ProgramName )

    ! --- Baryonic Mass ---

    BaryonicMass_FileName &
      = TRIM( BaseFileName ) // '_Tally_BaryonicMass.dat'

    BaryonicMass_Interior = Zero
    BaryonicMass_Initial  = Zero
    BaryonicMass_OffGrid  = Zero
    BaryonicMass_Change   = Zero

    TimeLabel     &
      = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
    InteriorLabel &
      = 'Interior [' // TRIM( UnitsDisplay % MassLabel ) // ']'
    OffGridLabel  &
      = 'Off Grid [' // TRIM( UnitsDisplay % MassLabel ) // ']'
    InitialLabel  &
      = 'Initial ['  // TRIM( UnitsDisplay % MassLabel ) // ']'
    ChangeLabel   &
      = 'Change ['   // TRIM( UnitsDisplay % MassLabel ) // ']'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( BaryonicMass_FileName ) )

    WRITE(FileUnit,'(5(A25,x))') &
      TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
      TRIM( InitialLabel ), TRIM( ChangeLabel )

    CLOSE( FileUnit )

    ! --- Energy ---

    Energy_FileName &
      = TRIM( BaseFileName ) // '_Tally_Energy.dat'

    Energy_Interior = Zero
    Energy_Initial  = Zero
    Energy_OffGrid  = Zero
    Energy_Change   = Zero

    TimeLabel     &
      = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
    InteriorLabel &
      = 'Interior [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    OffGridLabel  &
      = 'Off Grid [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    InitialLabel  &
      = 'Initial ['  // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    ChangeLabel   &
      = 'Change ['   // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( Energy_FileName ) )

    WRITE(FileUnit,'(5(A25,x))') &
      TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
      TRIM( InitialLabel ), TRIM( ChangeLabel )

    CLOSE( FileUnit )

    ! --- Internal Energy ---

    InternalEnergy_FileName &
      = TRIM( BaseFileName ) // '_Tally_InternalEnergy.dat'

    InternalEnergy_Interior = Zero
    InternalEnergy_Initial  = Zero
    InternalEnergy_OffGrid  = Zero
    InternalEnergy_Change   = Zero

    TimeLabel     &
      = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
    InteriorLabel &
      = 'Interior [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    OffGridLabel  &
      = 'Off Grid [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    InitialLabel  &
      = 'Initial ['  // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    ChangeLabel   &
      = 'Change ['   // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( InternalEnergy_FileName ) )

    WRITE(FileUnit,'(5(A25,x))') &
      TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
      TRIM( InitialLabel ), TRIM( ChangeLabel )

    CLOSE( FileUnit )

    ! --- Kinetic Energy ---

    KineticEnergy_FileName &
      = TRIM( BaseFileName ) // '_Tally_KineticEnergy.dat'

    KineticEnergy_Interior = Zero
    KineticEnergy_Initial  = Zero
    KineticEnergy_OffGrid  = Zero
    KineticEnergy_Change   = Zero

    TimeLabel     &
      = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
    InteriorLabel &
      = 'Interior [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    OffGridLabel  &
      = 'Off Grid [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    InitialLabel  &
      = 'Initial ['  // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    ChangeLabel   &
      = 'Change ['   // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( KineticEnergy_FileName ) )

    WRITE(FileUnit,'(5(A25,x))') &
      TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
      TRIM( InitialLabel ), TRIM( ChangeLabel )

    CLOSE( FileUnit )

    ! --- Gravitational Energy ---

    GravitationalEnergy_FileName &
      = TRIM( BaseFileName ) // '_Tally_GravitationalEnergy.dat'

    GravitationalEnergy_Interior = Zero
    GravitationalEnergy_Initial  = Zero
    GravitationalEnergy_OffGrid  = Zero
    GravitationalEnergy_Change   = Zero

    TimeLabel     &
      = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
    InteriorLabel &
      = 'Interior [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    OffGridLabel  &
      = 'Off Grid [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    InitialLabel  &
      = 'Initial ['  // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
    ChangeLabel   &
      = 'Change ['   // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( GravitationalEnergy_FileName ) )

    WRITE(FileUnit,'(5(A25,x))') &
      TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
      TRIM( InitialLabel ), TRIM( ChangeLabel )

    CLOSE( FileUnit )

  END SUBROUTINE InitializeTally_Euler_NonRelativistic


  SUBROUTINE ComputeTally_Euler_NonRelativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, Time, &
      SetInitialValues_Option, Verbose_Option )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      Time
    LOGICAL, INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: SetInitialValues
    LOGICAL :: Verbose

    IF( SuppressTally ) RETURN

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    CALL ComputeTally_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    IF( SetInitialValues )THEN

      BaryonicMass_Initial        = BaryonicMass_Interior
      Energy_Initial              = Energy_Interior
      InternalEnergy_Initial      = InternalEnergy_Interior
      KineticEnergy_Initial       = KineticEnergy_Interior
      GravitationalEnergy_Initial = GravitationalEnergy_Interior

    END IF

    BaryonicMass_Change &
      = BaryonicMass_Interior &
          - ( BaryonicMass_Initial + BaryonicMass_OffGrid )

    Energy_Change &
      = Energy_Interior &
          - ( Energy_Initial + Energy_OffGrid )

    InternalEnergy_Change &
      = InternalEnergy_Interior &
          - ( InternalEnergy_Initial + InternalEnergy_OffGrid )

    KineticEnergy_Change &
      = KineticEnergy_Interior &
          - ( KineticEnergy_Initial + KineticEnergy_OffGrid )

    GravitationalEnergy_Change &
      = GravitationalEnergy_Interior &
          - ( GravitationalEnergy_Initial + GravitationalEnergy_OffGrid )

    CALL WriteTally_Euler( Time )

    IF( Verbose )THEN

      CALL DisplayTally( Time )

    END IF

  END SUBROUTINE ComputeTally_Euler_NonRelativistic


  SUBROUTINE ComputeTally_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3

    REAL(DP) :: P(nDOFX,nPF)

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    BaryonicMass_Interior        = Zero
    Energy_Interior              = Zero
    InternalEnergy_Interior      = Zero
    KineticEnergy_Interior       = Zero
    GravitationalEnergy_Interior = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( U(iNX,iX1,iX2,iX3,iCF_D ), &
               U(iNX,iX1,iX2,iX3,iCF_S1), &
               U(iNX,iX1,iX2,iX3,iCF_S2), &
               U(iNX,iX1,iX2,iX3,iCF_S3), &
               U(iNX,iX1,iX2,iX3,iCF_E ), &
               U(iNX,iX1,iX2,iX3,iCF_Ne), &
               P(iNX,iPF_D ), &
               P(iNX,iPF_V1), &
               P(iNX,iPF_V2), &
               P(iNX,iPF_V3), &
               P(iNX,iPF_E ), &
               P(iNX,iPF_Ne), &
               G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      BaryonicMass_Interior &
        = BaryonicMass_Interior &
            + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U(iNX,iX1,iX2,iX3,iCF_D)

      Energy_Interior &
        = Energy_Interior &
            + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U(iNX,iX1,iX2,iX3,iCF_E)

      InternalEnergy_Interior &
        = InternalEnergy_Interior &
            + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * P(iNX,iPF_E)

      KineticEnergy_Interior &
        = KineticEnergy_Interior &
            + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * Half * ( U(iNX,iX1,iX2,iX3,iCF_S1) * P(iNX,iPF_V1) &
                           + U(iNX,iX1,iX2,iX3,iCF_S2) * P(iNX,iPF_V2) &
                           + U(iNX,iX1,iX2,iX3,iCF_S3) * P(iNX,iPF_V3) )

      GravitationalEnergy_Interior &
        = GravitationalEnergy_Interior &
            + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                * WeightsX_q(iNX) &
                * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                * Half * P(iNX,iPF_D) * G(iNX,iX1,iX2,iX3,iGF_Phi_N)

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTally_Euler


  SUBROUTINE IncrementOffGridTally_Euler_NonRelativistic( dM )

    REAL(DP), INTENT(in) :: dM(nCF)

    BaryonicMass_OffGrid &
      = dM(iCF_D)

    Energy_OffGrid &
      = dM(iCF_E)

    InternalEnergy_OffGrid &
      = Zero

    KineticEnergy_OffGrid &
      = Zero

    GravitationalEnergy_OffGrid &
      = Zero

  END SUBROUTINE IncrementOffGridTally_Euler_NonRelativistic


  SUBROUTINE WriteTally_Euler( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    ! --- Baryonic Mass ---

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( BaryonicMass_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, &
      BaryonicMass_Interior / UnitsDisplay % MassUnit, &
      BaryonicMass_OffGrid  / UnitsDisplay % MassUnit, &
      BaryonicMass_Initial  / UnitsDisplay % MassUnit, &
      BaryonicMass_Change   / UnitsDisplay % MassUnit

    CLOSE( FileUnit )

    ! --- Energy ---

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( Energy_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, &
      Energy_Interior / UnitsDisplay % EnergyGlobalUnit, &
      Energy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
      Energy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
      Energy_Change   / UnitsDisplay % EnergyGlobalUnit

    CLOSE( FileUnit )

    ! --- Internal Energy ---

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( InternalEnergy_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, &
      InternalEnergy_Interior / UnitsDisplay % EnergyGlobalUnit, &
      InternalEnergy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
      InternalEnergy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
      InternalEnergy_Change   / UnitsDisplay % EnergyGlobalUnit

    CLOSE( FileUnit )

    ! --- Kinetic Energy ---

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( KineticEnergy_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, &
      KineticEnergy_Interior / UnitsDisplay % EnergyGlobalUnit, &
      KineticEnergy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
      KineticEnergy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
      KineticEnergy_Change   / UnitsDisplay % EnergyGlobalUnit

    CLOSE( FileUnit )

    ! --- Gravitational Energy ---

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( GravitationalEnergy_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, &
      GravitationalEnergy_Interior / UnitsDisplay % EnergyGlobalUnit, &
      GravitationalEnergy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
      GravitationalEnergy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
      GravitationalEnergy_Change   / UnitsDisplay % EnergyGlobalUnit

    CLOSE( FileUnit )

  END SUBROUTINE WriteTally_Euler


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    WRITE(*,*)
    WRITE(*,'(A8,A,ES8.2E2,x,A)') &
      '', 'Euler Tally. t = ', &
      Time / UnitsDisplay % TimeUnit, &
      UnitsDisplay % TimeLabel
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Baryonic Mass Interior.: ', &
      BaryonicMass_Interior / UnitsDisplay % MassUnit, &
      UnitsDisplay % MassLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Baryonic Mass Initial..: ', &
      BaryonicMass_Initial  / UnitsDisplay % MassUnit, &
      UnitsDisplay % MassLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Baryonic Mass Off Grid.: ', &
      BaryonicMass_OffGrid  / UnitsDisplay % MassUnit, &
      UnitsDisplay % MassLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Baryonic Mass Change...: ', &
      BaryonicMass_Change   / UnitsDisplay % MassUnit, &
      UnitsDisplay % MassLabel
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Energy Interior.: ', &
      Energy_Interior / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Energy Initial..: ', &
      Energy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Energy Off Grid.: ', &
      Energy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Energy Change...: ', &
      Energy_Change  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Internal Energy Interior.: ', &
      InternalEnergy_Interior / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Internal Energy Initial..: ', &
      InternalEnergy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Internal Energy Off Grid.: ', &
      InternalEnergy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Internal Energy Change...: ', &
      InternalEnergy_Change   / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Kinetic Energy Interior.: ', &
      KineticEnergy_Interior / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Kinetic Energy Initial..: ', &
      KineticEnergy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Kinetic Energy Off Grid.: ', &
      KineticEnergy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Kinetic Energy Change...: ', &
      KineticEnergy_Change   / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Gravitational Energy Interior.: ', &
      GravitationalEnergy_Interior / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Gravitational Energy Initial..: ', &
      GravitationalEnergy_Initial  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Gravitational Energy Off Grid.: ', &
      GravitationalEnergy_OffGrid  / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'Gravitational Energy Change...: ', &
      GravitationalEnergy_Change   / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel

    WRITE(*,*)

  END SUBROUTINE DisplayTally


  SUBROUTINE FinalizeTally_Euler_NonRelativistic

  END SUBROUTINE FinalizeTally_Euler_NonRelativistic


END MODULE Euler_TallyModule_NonRelativistic
