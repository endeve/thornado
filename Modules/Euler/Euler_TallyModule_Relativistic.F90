MODULE Euler_TallyModule_Relativistic

  USE KindModule, ONLY: &
    DP, &
    Zero, &
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
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_E
  USE UnitsModule, ONLY: &
    UnitsDisplay

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE ADM_Mass_Module, ONLY: &
    Calc_ADM_Mass

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_Euler_Relativistic
  PUBLIC :: ComputeTally_Euler_Relativistic
  PUBLIC :: IncrementOffGridTally_Euler_Relativistic
  PUBLIC :: FinalizeTally_Euler_Relativistic

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

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  CHARACTER(256) :: ADMMass_FileName
  REAL(DP)       :: ADMMass_Interior
  REAL(DP)       :: ADMMass_Initial
  REAL(DP)       :: ADMMass_OffGrid
  REAL(DP)       :: ADMMass_Change

#endif


CONTAINS


  SUBROUTINE InitializeTally_Euler_Relativistic &
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

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ! --- ADM Mass ---

    ADMMass_FileName &
      = TRIM( BaseFileName ) // '_Tally_ADMMass.dat'

    ADMMass_Interior = Zero
    ADMMass_Initial  = Zero
    ADMMass_OffGrid  = Zero
    ADMMass_Change   = Zero

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

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( ADMMass_FileName ) )

    WRITE(FileUnit,'(5(A25,x))') &
      TRIM( TimeLabel ), TRIM( InteriorLabel ), TRIM( OffGridLabel ), &
      TRIM( InitialLabel ), TRIM( ChangeLabel )

    CLOSE( FileUnit )

#endif

  END SUBROUTINE InitializeTally_Euler_Relativistic


  SUBROUTINE ComputeTally_Euler_Relativistic &
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

      BaryonicMass_Initial = BaryonicMass_Interior
      Energy_Initial       = Energy_Interior

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      ADMMass_Initial      = ADMMass_Interior

#endif


    END IF

    BaryonicMass_Change &
      = BaryonicMass_Interior &
          - ( BaryonicMass_Initial + BaryonicMass_OffGrid )

    Energy_Change &
      = Energy_Interior &
          - ( Energy_Initial + Energy_OffGrid )

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ADMMass_Change &
      = ADMMass_Interior &
          - ( ADMMass_Initial + ADMMass_OffGrid )

#endif

    CALL WriteTally_Euler( Time )

    IF( Verbose )THEN

      CALL DisplayTally( Time )

    END IF

  END SUBROUTINE ComputeTally_Euler_Relativistic


  SUBROUTINE ComputeTally_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    BaryonicMass_Interior = Zero
    Energy_Interior       = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

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

    END DO
    END DO
    END DO
    END DO

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL Calc_ADM_Mass( ADMMass_Interior )

#endif

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTally_Euler


  SUBROUTINE IncrementOffGridTally_Euler_Relativistic( dM )

    REAL(DP), INTENT(in) :: dM(nCF)

    BaryonicMass_OffGrid &
      = BaryonicMass_OffGrid + dM(iCF_D)

    Energy_OffGrid &
      = Energy_OffGrid + dM(iCF_E)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ADMMass_OffGrid &
      = Zero

#endif

  END SUBROUTINE IncrementOffGridTally_Euler_Relativistic


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

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ! --- ADM Mass ---

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( ADMMass_FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(5(ES25.16E3,1x))' ) &
      Time / UnitsDisplay % TimeUnit, &
      ADMMass_Interior / UnitsDisplay % MassUnit, &
      ADMMass_OffGrid  / UnitsDisplay % MassUnit, &
      ADMMass_Initial  / UnitsDisplay % MassUnit, &
      ADMMass_Change   / UnitsDisplay % MassUnit

    CLOSE( FileUnit )

#endif

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
      Energy_Change   / UnitsDisplay % EnergyGlobalUnit, &
      UnitsDisplay % EnergyGlobalLabel

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    WRITE(*,*)

    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'ADM Mass Interior.: ', &
      ADMMass_Interior / UnitsDisplay % MassUnit, &
      UnitsDisplay % MassLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'ADM Mass Initial..: ', &
      ADMMass_Initial  / UnitsDisplay % MassUnit, &
      UnitsDisplay % MassLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'ADM Mass Off Grid.: ', &
      ADMMass_OffGrid  / UnitsDisplay % MassUnit, &
      UnitsDisplay % MassLabel
    WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
      '', 'ADM Mass Change...: ', &
      ADMMass_Change   / UnitsDisplay % MassUnit, &
      UnitsDisplay % MassLabel

#endif

    WRITE(*,*)

  END SUBROUTINE DisplayTally


  SUBROUTINE FinalizeTally_Euler_Relativistic

  END SUBROUTINE FinalizeTally_Euler_Relativistic


END MODULE Euler_TallyModule_Relativistic
