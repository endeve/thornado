#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_EOS
#endif

MODULE EquationOfStateModule_TABLE

#ifdef MICROPHYSICS_WEAKLIB

  ! --- weaklib modules --------------------------

  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType
  USE wlEOSInversionModule, ONLY: &
    InitializeEOSInversion, &
    ComputeTemperatureWith_DEY_Single_Guess, &
    ComputeTemperatureWith_DEY_Single_NoGuess, &
    ComputeTemperatureWith_DPY_Single_NoGuess, &
    DescribeEOSInversionError
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_3D_Custom_Point, &
    LogInterpolateDifferentiateSingleVariable_3D_Custom_Point

  ! ----------------------------------------------

#endif

  USE, INTRINSIC :: ISO_C_BINDING
  USE KindModule, ONLY: &
    DP, Zero, One
  USE UnitsModule, ONLY: &
    AtomicMassUnit, &
    BoltzmannConstant, &
    Gram, &
    Centimeter, &
    Kelvin, &
    Dyne, &
    Erg, &
    MeV

  IMPLICIT NONE
  PRIVATE

  CHARACTER(256) :: &
    EquationOfStateTableName
  INTEGER :: &
    iD_T, iT_T, iY_T, &
    iP_T, iS_T, iE_T, iMe_T, iMp_T, iMn_T, &
    iXp_T, iXn_T, iXa_T, iXh_T, iGm_T
  REAL(DP) :: &
    UnitD, UnitT, UnitY, &
    UnitP, UnitS, UnitE, UnitMe, UnitMp, UnitMn, &
    UnitXp, UnitXn, UnitXa, UnitXh, UnitGm
  REAL(DP), PUBLIC :: &
    OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, &
    OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm
  REAL(DP), PARAMETER :: &
    BaryonMass = AtomicMassUnit
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    D_T, T_T, Y_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    P_T, S_T, E_T, Me_T, Mp_T, Mn_T, &
    Xp_T, Xn_T, Xa_T, Xh_T, Gm_T
#ifdef MICROPHYSICS_WEAKLIB
  LOGICAL :: UsingExternalEOS
  TYPE(EquationOfStateTableType), POINTER :: EOS
#endif

  PUBLIC :: InitializeEquationOfState_TABLE
  PUBLIC :: FinalizeEquationOfState_TABLE
  PUBLIC :: ApplyEquationOfState_TABLE
  PUBLIC :: ComputeTemperatureFromPressure_TABLE
  PUBLIC :: ComputeTemperatureFromSpecificInternalEnergy_TABLE
  PUBLIC :: ComputeThermodynamicStates_Primitive_TABLE
  PUBLIC :: ComputeThermodynamicStates_Auxiliary_TABLE
  PUBLIC :: ComputePressureFromPrimitive_TABLE
  PUBLIC :: ComputePressureFromSpecificInternalEnergy_TABLE
  PUBLIC :: ComputeSoundSpeedFromPrimitive_TABLE
  PUBLIC :: ComputeAuxiliary_Fluid_TABLE
  PUBLIC :: ComputePressure_TABLE
  PUBLIC :: ComputeSpecificInternalEnergy_TABLE
  PUBLIC :: ComputeElectronChemicalPotential_TABLE
  PUBLIC :: ComputeProtonChemicalPotential_TABLE
  PUBLIC :: ComputeNeutronChemicalPotential_TABLE
  PUBLIC :: ComputeProtonMassFraction_TABLE
  PUBLIC :: ComputeNeutronMassFraction_TABLE

  REAL(DP), PUBLIC :: Min_D, Min_T, Min_Y
  REAL(DP), PUBLIC :: Max_D, Max_T, Max_Y

  INTERFACE ApplyEquationOfState_TABLE
    MODULE PROCEDURE ApplyEquationOfState_Table_Scalar
    MODULE PROCEDURE ApplyEquationOfState_Table_Vector
  END INTERFACE ApplyEquationOfState_TABLE

  INTERFACE ComputePressure_TABLE
    MODULE PROCEDURE ComputePressure_TABLE_Scalar
    MODULE PROCEDURE ComputePressure_TABLE_Vector
  END INTERFACE ComputePressure_TABLE

  INTERFACE ComputeSpecificInternalEnergy_TABLE
    MODULE PROCEDURE ComputeSpecificInternalEnergy_TABLE_Scalar
    MODULE PROCEDURE ComputeSpecificInternalEnergy_TABLE_Vector
  END INTERFACE ComputeSpecificInternalEnergy_TABLE

  INTERFACE ComputePressureFromPrimitive_TABLE
    MODULE PROCEDURE ComputePressureFromPrimitive_TABLE_Scalar
    MODULE PROCEDURE ComputePressureFromPrimitive_TABLE_Vector
  END INTERFACE ComputePressureFromPrimitive_TABLE

  INTERFACE ComputePressureFromSpecificInternalEnergy_TABLE
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_TABLE_Scalar
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_TABLE_Vector
  END INTERFACE ComputePressureFromSpecificInternalEnergy_TABLE

  INTERFACE ComputeSoundSpeedFromPrimitive_TABLE
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_TABLE_Scalar
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_TABLE_Vector
  END INTERFACE ComputeSoundSpeedFromPrimitive_TABLE

  INTERFACE ComputeThermodynamicStates_Primitive_TABLE
    MODULE PROCEDURE ComputeThermodynamicStates_Primitive_TABLE_Scalar
    MODULE PROCEDURE ComputeThermodynamicStates_Primitive_TABLE_Vector
  END INTERFACE ComputeThermodynamicStates_Primitive_TABLE

  INTERFACE ComputeThermodynamicStates_Auxiliary_TABLE
    MODULE PROCEDURE ComputeThermodynamicStates_Auxiliary_TABLE_Scalar
    MODULE PROCEDURE ComputeThermodynamicStates_Auxiliary_TABLE_Vector
  END INTERFACE ComputeThermodynamicStates_Auxiliary_TABLE

  INTERFACE ComputeAuxiliary_Fluid_TABLE
    MODULE PROCEDURE ComputeAuxiliary_Fluid_TABLE_Scalar
    MODULE PROCEDURE ComputeAuxiliary_Fluid_TABLE_Vector
  END INTERFACE ComputeAuxiliary_Fluid_TABLE

  INTERFACE ComputeTemperatureFromSpecificInternalEnergy_TABLE
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector
  END INTERFACE ComputeTemperatureFromSpecificInternalEnergy_TABLE

  INTERFACE ComputeTemperatureFromPressure_TABLE
    MODULE PROCEDURE ComputeTemperatureFromPressure_TABLE_Scalar
    MODULE PROCEDURE ComputeTemperatureFromPressure_TABLE_Vector
  END INTERFACE ComputeTemperatureFromPressure_TABLE

  INTERFACE ComputeElectronChemicalPotential_TABLE
    MODULE PROCEDURE ComputeElectronChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeElectronChemicalPotential_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeProtonChemicalPotential_TABLE
    MODULE PROCEDURE ComputeProtonChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeProtonChemicalPotential_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeNeutronChemicalPotential_TABLE
    MODULE PROCEDURE ComputeNeutronChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeNeutronChemicalPotential_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeProtonMassFraction_TABLE
    MODULE PROCEDURE ComputeProtonMassFraction_TABLE_Scalar
    MODULE PROCEDURE ComputeProtonMassFraction_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeNeutronMassFraction_TABLE
    MODULE PROCEDURE ComputeNeutronMassFraction_TABLE_Scalar
    MODULE PROCEDURE ComputeNeutronMassFraction_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeDependentVariable_TABLE
    MODULE PROCEDURE ComputeDependentVariable_TABLE_Scalar
    MODULE PROCEDURE ComputeDependentVariable_TABLE_Vector
  END INTERFACE ComputeDependentVariable_TABLE

  INTERFACE ComputeDependentVariableAndDerivatives_TABLE
    MODULE PROCEDURE ComputeDependentVariableAndDerivatives_TABLE_Scalar
    MODULE PROCEDURE ComputeDependentVariableAndDerivatives_TABLE_Vector
  END INTERFACE ComputeDependentVariableAndDerivatives_TABLE

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( D_T, T_T, Y_T, &
  !$OMP   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, UnitMe, UnitMp, UnitMn, &
  !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, &
  !$OMP   OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$OMP   OS_Xa, OS_Xh, OS_Gm, &
  !$OMP   P_T, S_T, E_T, Me_T, Mp_T, Mn_T, Xp_T, Xn_T, &
  !$OMP   Xa_T, Xh_T, Gm_T, &
  !$OMP   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( D_T, T_T, Y_T, &
  !$ACC   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, UnitMe, UnitMp, UnitMn, &
  !$ACC   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, &
  !$ACC   OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$ACC   OS_Xa, OS_Xh, OS_Gm, &
  !$ACC   P_T, S_T, E_T, Me_T, Mp_T, Mn_T, Xp_T, Xn_T, &
  !$ACC   Xa_T, Xh_T, Gm_T, &
  !$ACC   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
#endif

CONTAINS


  SUBROUTINE InitializeEquationOfState_TABLE &
    ( EquationOfStateTableName_Option, Verbose_Option, External_EOS )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option
#ifdef MICROPHYSICS_WEAKLIB
    TYPE(EquationOfStateTableType), POINTER, &
                      INTENT(in), OPTIONAL :: External_EOS
#else
    INTEGER,          INTENT(in), OPTIONAL :: External_EOS
#endif

    LOGICAL :: Verbose

    IF( PRESENT( EquationOfStateTableName_Option ) )THEN
       EquationOfStateTableName = TRIM( EquationOfStateTableName_Option )
    ELSE
       EquationOfStateTableName = 'EquationOfStateTable.h5'
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
       Verbose = Verbose_Option
    ELSE
       Verbose = .FALSE.
    END IF

    IF( Verbose )THEN
       WRITE(*,*)
       WRITE(*,'(A7,A12,A)') &
            '', 'Table Name: ', TRIM( EquationOfStateTableName )
    END IF

#ifdef MICROPHYSICS_WEAKLIB

    IF( .NOT. PRESENT( External_EOS ) ) THEN

       ALLOCATE( EOS )
       UsingExternalEOS = .FALSE.

       CALL InitializeHDF( )

       CALL ReadEquationOfStateTableHDF &
            ( EOS, TRIM( EquationOfStateTableName ) )

       CALL FinalizeHDF( )

    ELSE

       EOS => External_EOS
       UsingExternalEOS = .TRUE.

    END IF

    ! --- Thermodynamic State Indices ---

    iD_T = EOS % TS % Indices % iRho
    iT_T = EOS % TS % Indices % iT
    iY_T = EOS % TS % Indices % iYe

    ! --- Units ---

    UnitD = Gram / Centimeter**3
    UnitT = Kelvin
    UnitY = One

    UnitP  = Dyne / Centimeter**2
    UnitS  = BoltzmannConstant
    UnitE  = Erg / Gram
    UnitMe = MeV
    UnitMp = MeV
    UnitMn = MeV
    UnitXp = One
    UnitXn = One
    UnitXa = One
    UnitXh = One
    UnitGm = One

    ! --- Thermodynamic States ---

    ALLOCATE( D_T(EOS % TS % nPoints(iD_T)) )
    D_T = EOS % TS % States(iD_T) % Values

    Min_D = MINVAL( D_T ) * Gram / Centimeter**3
    Max_D = MAXVAL( D_T ) * Gram / Centimeter**3

    ALLOCATE( T_T(EOS % TS % nPoints(iT_T)) )
    T_T = EOS % TS % States(iT_T) % Values

    Min_T = MINVAL( T_T ) * Kelvin
    Max_T = MAXVAL( T_T ) * Kelvin

    ALLOCATE( Y_T(EOS % TS % nPoints(iY_T)) )
    Y_T = EOS % TS % States(iY_T) % Values

    Min_Y = MINVAL( Y_T )
    Max_Y = MAXVAL( Y_T )

    ! --- Dependent Variables Indices ---

    iP_T  = EOS % DV % Indices % iPressure
    iS_T  = EOS % DV % Indices % iEntropyPerBaryon
    iE_T  = EOS % DV % Indices % iInternalEnergyDensity
    iMe_T = EOS % DV % Indices % iElectronChemicalPotential
    iMp_T = EOS % DV % Indices % iProtonChemicalPotential
    iMn_T = EOS % DV % Indices % iNeutronChemicalPotential
    iXp_T = EOS % DV % Indices % iProtonMassFraction
    iXn_T = EOS % DV % Indices % iNeutronMassFraction
    iXa_T = EOS % DV % Indices % iAlphaMassFraction
    iXh_T = EOS % DV % Indices % iHeavyMassFraction
    iGm_T = EOS % DV % Indices % iGamma1

    ! --- Dependent Variables Offsets ---

    OS_P  = EOS % DV % Offsets(iP_T)
    OS_S  = EOS % DV % Offsets(iS_T)
    OS_E  = EOS % DV % Offsets(iE_T)
    OS_Me = EOS % DV % Offsets(iMe_T)
    OS_Mp = EOS % DV % Offsets(iMp_T)
    OS_Mn = EOS % DV % Offsets(iMn_T)
    OS_Xp = EOS % DV % Offsets(iXp_T)
    OS_Xn = EOS % DV % Offsets(iXn_T)
    OS_Xa = EOS % DV % Offsets(iXa_T)
    OS_Xh = EOS % DV % Offsets(iXh_T)
    OS_Gm = EOS % DV % Offsets(iGm_T)

    ALLOCATE &
      ( P_T (1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( S_T (1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( E_T (1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Me_T(1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Mp_T(1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Mn_T(1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Xp_T(1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Xn_T(1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Xa_T(1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Xh_T(1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Gm_T(1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )

    P_T  = EOS % DV % Variables(iP_T ) % Values
    S_T  = EOS % DV % Variables(iS_T ) % Values
    E_T  = EOS % DV % Variables(iE_T ) % Values
    Me_T = EOS % DV % Variables(iMe_T) % Values
    Mp_T = EOS % DV % Variables(iMp_T) % Values
    Mn_T = EOS % DV % Variables(iMn_T) % Values
    Xp_T = EOS % DV % Variables(iXp_T) % Values
    Xn_T = EOS % DV % Variables(iXn_T) % Values
    Xa_T = EOS % DV % Variables(iXa_T) % Values
    Xh_T = EOS % DV % Variables(iXh_T) % Values
    Gm_T = EOS % DV % Variables(iGm_T) % Values

    CALL InitializeEOSInversion &
           ( D_T, T_T, Y_T, &
             10.0d0**( E_T ) - OS_E, &
             10.0d0**( P_T ) - OS_P, &
             10.0d0**( S_T ) - OS_S, &
             Verbose_Option = Verbose )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO &
    !$OMP ( D_T, T_T, Y_T, &
    !$OMP   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMe, UnitMp, UnitMn, &
    !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, OS_P, OS_S, OS_E, OS_Me, &
    !$OMP   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm, P_T, S_T, &
    !$OMP   E_T, Me_T, Mp_T, Mn_T, Xp_T, Xn_T, Xa_T, Xh_T, Gm_T, &
    !$OMP   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( D_T, T_T, Y_T, &
    !$ACC   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMe, UnitMp, UnitMn, &
    !$ACC   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, OS_P, OS_S, OS_E, OS_Me, &
    !$ACC   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm, P_T, S_T, &
    !$ACC   E_T, Me_T, Mp_T, Mn_T, Xp_T, Xn_T, Xa_T, Xh_T, Gm_T, &
    !$ACC   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
#endif

#endif

  END SUBROUTINE InitializeEquationOfState_TABLE


  SUBROUTINE FinalizeEquationOfState_TABLE

#ifdef MICROPHYSICS_WEAKLIB

    DEALLOCATE( D_T, T_T, Y_T )

    DEALLOCATE( P_T  )
    DEALLOCATE( S_T  )
    DEALLOCATE( E_T  )
    DEALLOCATE( Me_T )
    DEALLOCATE( Mp_T )
    DEALLOCATE( Mn_T )
    DEALLOCATE( Xp_T )
    DEALLOCATE( Xn_T )
    DEALLOCATE( Xa_T )
    DEALLOCATE( Xh_T )
    DEALLOCATE( Gm_T )

    IF ( .NOT. UsingExternalEOS ) THEN
       DEALLOCATE( EOS )
    END IF

#endif

  END SUBROUTINE FinalizeEquationOfState_TABLE


  SUBROUTINE ApplyEquationOfState_TABLE_Scalar &
    ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm

    ! --- Interpolate Pressure ----------------------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, P, P_T, OS_P, Units_V = UnitP )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, S, S_T, OS_S, Units_V = UnitS )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, E, E_T, OS_E, Units_V = UnitE )

    ! --- Interpolate Electron Chemical Potential ---------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Me, Me_T, OS_Me, Units_V = UnitMe )

    ! --- Interpolate Proton Chemical Potential -----------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Mp, Mp_T, OS_Mp, Units_V = UnitMp )

    ! --- Interpolate Neutron Chemical Potential ----------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Mn, Mn_T, OS_Mn, Units_V = UnitMn )

    ! --- Interpolate Proton Mass Fraction ----------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Xp, Xp_T, OS_Xp, Units_V = UnitXp )

    ! --- Interpolate Neutron Mass Fraction ---------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Xn, Xn_T, OS_Xn, Units_V = UnitXn )

    ! --- Interpolate Alpha Mass Fraction -----------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Xa, Xa_T, OS_Xa, Units_V = UnitXa )

    ! --- Interpolate Heavy Mass Fraction -----------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Xh, Xh_T, OS_Xh, Units_V = UnitXh )

    ! --- Gamma1 ------------------------------------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Gm, Gm_T, OS_Gm, Units_V = UnitGm )

  END SUBROUTINE ApplyEquationOfState_TABLE_Scalar


  SUBROUTINE ApplyEquationOfState_TABLE_Vector &
    ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out) :: P(1:), S(1:), E(1:), Me(1:), Mp(1:), Mn(1:)
    REAL(DP), INTENT(out) :: Xp(1:), Xn(1:), Xa(1:), Xh(1:), Gm(1:)

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ApplyEquationOfState_TABLE_Scalar &
             ( D (iP), T (iP), Y (iP), P (iP), S (iP), E (iP), Me(iP), &
               Mp(iP), Mn(iP), Xp(iP), Xn(iP), Xa(iP), Xh(iP), Gm(iP) )

    END DO

  END SUBROUTINE ApplyEquationOfState_TABLE_Vector


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
    ( D, E, Y, T, Guess_Option, Error_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)            :: D, E, Y
    REAL(DP), INTENT(out)           :: T
    REAL(DP), INTENT(in),  OPTIONAL :: Guess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    REAL(DP) :: D_P, E_P, Y_P, T_Lookup, T_Guess
    INTEGER  :: Error

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / ( Gram / Centimeter**3 )
    E_P = E / ( Erg / Gram )
    Y_P = Y

    IF ( PRESENT( Guess_Option ) ) THEN

      T_Guess = Guess_Option / Kelvin

      CALL ComputeTemperatureWith_DEY_Single_Guess &
             ( D_P, E_P, Y_P, D_T, T_T, Y_T, E_T, OS_E, T_Lookup, T_Guess, &
               Error_Option = Error )

    ELSE

      CALL ComputeTemperatureWith_DEY_Single_NoGuess &
             ( D_P, E_P, Y_P, D_T, T_T, Y_T, E_T, OS_E, T_Lookup, &
               Error_Option = Error )

    END IF
    T = T_Lookup * Kelvin

#endif

    IF ( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector &
    ( D, E, Y, T, Guess_Option, Error_Option )

    REAL(DP), INTENT(in )           :: D(1:), E(1:), Y(1:)
    REAL(DP), INTENT(out)           :: T(1:)
    REAL(DP), INTENT(in ), OPTIONAL :: Guess_Option(1:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(1:)

    INTEGER  :: iP, nP
    INTEGER  :: Error(SIZE(D))
    REAL(DP) :: D_P, E_P, Y_P, T_Lookup, T_Guess

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

    IF( PRESENT( Guess_Option ) )THEN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, E_P, Y_P, T_Lookup, T_Guess ) &
    !$OMP MAP( to: D, E, Y, D_T, T_T, Y_T, E_T, Guess_Option ) &
    !$OMP MAP( from: T, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, E_P, Y_P, T_Lookup, T_Guess ) &
    !$ACC COPYIN( D, E, Y, D_T, T_T, Y_T, E_T, Guess_Option ) &
    !$ACC COPYOUT( T, Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, E_P, Y_P, T_Lookup, T_Guess )
#endif
      DO iP = 1, nP

        D_P = D(iP) / ( Gram / Centimeter**3 )
        E_P = E(iP) / ( Erg / Gram )
        Y_P = Y(iP)

        T_Guess = Guess_Option(iP) / Kelvin

        CALL ComputeTemperatureWith_DEY_Single_Guess &
               ( D_P, E_P, Y_P, D_T, T_T, Y_T, E_T, OS_E, T_Lookup, T_Guess, &
                 Error_Option = Error(iP) )

        T(iP) = T_Lookup * Kelvin

      END DO

    ELSE

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, E_P, Y_P, T_Lookup ) &
    !$OMP MAP( to: D, E, Y, D_T, T_T, Y_T, E_T ) &
    !$OMP MAP( from: T, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, E_P, Y_P, T_Lookup ) &
    !$ACC COPYIN( D, E, Y, D_T, T_T, Y_T, E_T ) &
    !$ACC COPYOUT( T, Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, E_P, Y_P, T_Lookup )
#endif
      DO iP = 1, nP

        D_P = D(iP) / ( Gram / Centimeter**3 )
        E_P = E(iP) / ( Erg / Gram )
        Y_P = Y(iP)

        CALL ComputeTemperatureWith_DEY_Single_NoGuess &
               ( D_P, E_P, Y_P, D_T, T_T, Y_T, E_T, OS_E, T_Lookup, &
                 Error_Option = Error(iP) )

        T(iP) = T_Lookup * Kelvin

      END DO

    END IF

    IF ( ANY( Error > 0 ) ) THEN
      DO iP = 1, nP
        IF ( Error(iP) > 0 ) THEN
          CALL DescribeEOSInversionError( Error(iP) )
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM &
          !$OMP ( D(iP), E(iP), Y(iP) )
#elif defined(THORNADO_OACC)
          !$ACC UPDATE HOST &
          !$ACC ( D(iP), E(iP), Y(iP) )
#endif
          D_P = D(iP) / ( Gram / Centimeter**3 )
          E_P = E(iP) / ( Erg / Gram )
          Y_P = Y(iP)
          WRITE(*,*)                 '[ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector] Error'
          WRITE(*,'(a,i5,3es23.15)') '  iP, D, E, Y : ', iP, D_P, E_P, Y_P
        END IF
      END DO
      STOP
    END IF

#endif

    IF( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector


  SUBROUTINE ComputePressureFromPrimitive_TABLE_Scalar( D, Ev, Ne, P )

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
  !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: P

    REAL(DP) :: Em, T, Y

    Em = Ev / D              ! --- Internal Energy per Mass
    Y  = Ne / D * BaryonMass ! --- Electron Fraction

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
           ( D, Em, Y, T )

    CALL ComputePressure_TABLE_Scalar &
           ( D, T, Y, P )

  END SUBROUTINE ComputePressureFromPrimitive_TABLE_Scalar


  SUBROUTINE ComputePressureFromPrimitive_TABLE_Vector( D, Ev, Ne, P )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(1:), INTENT(out) :: P

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputePressureFromPrimitive_TABLE_Scalar &
             ( D(iP), Ev(iP), Ne(iP), P(iP) )

    END DO

  END SUBROUTINE ComputePressureFromPrimitive_TABLE_Vector


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Scalar &
    ( D, Em, Y, P )

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
  !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Em, Y
    REAL(DP), INTENT(out) :: P

    REAL(DP) :: D_P, E_P, Y_P, T_P, T

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D  / UnitD
    E_P = Em / UnitE
    Y_P = Y  / UnitY

    CALL ComputeTemperatureWith_DEY_Single_NoGuess &
           ( D_P, E_P, Y_P, D_T, T_T, Y_T, E_T, OS_E, T_P )

    T = T_P * UnitT

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, P, P_T, OS_P, Units_V = UnitP )

#endif

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Vector &
    ( D, Em, Y, P )

    REAL(DP), INTENT(in)  :: D(:), Em(:), Y(:)
    REAL(DP), INTENT(out) :: P(:)

    INTEGER iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputePressureFromSpecificInternalEnergy_TABLE_Scalar &
             ( D(iP), Em(iP), Y(iP), P(iP) )

    END DO

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Vector


  SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Scalar( D, Ev, Ne, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: Cs

    REAL(DP) :: P, T, Y, Em, Gm

    Em = Ev / D
    Y  = Ne * ( BaryonMass / D )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
           ( D, Em, Y, T )

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, P, P_T, OS_P, Units_V = Dyne / Centimeter**2 )

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Gm, Gm_T, OS_Gm, Units_V = 1.0_DP )

    Cs = SQRT( Gm * P / D )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Scalar


  SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Vector( D, Ev, Ne, Cs )

    REAL(DP), INTENT(in)  :: D(1:), Ev(1:), Ne(1:)
    REAL(DP), INTENT(out) :: Cs(1:)

    REAL(DP), DIMENSION(SIZE(D)) :: P, T, Y, Em, Gm

    Em = Ev / D
    Y  = Ne * ( BaryonMass / D )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector &
           ( D, Em, Y, T )

    CALL ComputeDependentVariable_TABLE_Vector &
           ( D, T, Y, P, P_T, OS_P, Units_V = Dyne / Centimeter**2 )

    CALL ComputeDependentVariable_TABLE_Vector &
           ( D, T, Y, Gm, Gm_T, OS_Gm, Units_V = 1.0_DP )

    Cs = SQRT( Gm * P / D )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Vector


  SUBROUTINE ComputeTemperatureFromPressure_TABLE_Scalar( D, P, Y, T )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, P, Y
    REAL(DP), INTENT(out) :: T

    INTEGER  :: Error
    REAL(DP) :: D_P, P_P, Y_P, T_P

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    P_P = P / UnitP
    Y_P = Y / UnitY

    CALL ComputeTemperatureWith_DPY_Single_NoGuess &
           ( D_P, P_P, Y_P, D_T, T_T, Y_T, P_T, OS_P, T_P, &
             Error_Option = Error )

    T = T_P * UnitT

#endif

  END SUBROUTINE ComputeTemperatureFromPressure_TABLE_Scalar


  SUBROUTINE ComputeTemperatureFromPressure_TABLE_Vector( D, P, Y, T )

    REAL(DP), INTENT(in)  :: D(1:), P(1:), Y(1:)
    REAL(DP), INTENT(out) :: T(1:)

    INTEGER :: iP, nP

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputeTemperatureFromPressure_TABLE_Scalar &
             ( D(iP), P(iP), Y(iP), T(iP) )

    END DO

#endif

  END SUBROUTINE ComputeTemperatureFromPressure_TABLE_Vector


  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Scalar &
    ( D, T, Y, Ev, Em, Ne )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: Ev, Em, Ne

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Em, E_T, OS_E, Units_V = UnitE )

      Ev = Em * D              ! --- Internal Energy per Unit Volume
      Ne = Y  * D / BaryonMass ! --- Electrons per Unit Volume

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Scalar


  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Vector &
    ( D, T, Y, Ev, Em, Ne )

    REAL(DP), INTENT(in)  :: D (1:), T (1:), Y (1:)
    REAL(DP), INTENT(out) :: Ev(1:), Em(1:), Ne(1:)

    INTEGER :: iP, nP

    nP = SIZE( D )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: D, T, Y ) &
    !$OMP MAP( from: Em, Ev, Ne )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( D, T, Y ) &
    !$ACC COPYOUT( Em, Ev, Ne )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeThermodynamicStates_Primitive_TABLE_Scalar &
             ( D(iP), T(iP), Y(iP), Ev(iP), Em(iP), Ne(iP) )

    END DO

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Vector


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Scalar &
    ( D, Ev, Ne, T, Em, Y )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: T, Em, Y

    Em = Ev / D              ! --- Internal Energy per Mass
    Y  = Ne / D * BaryonMass ! --- Electron Fraction

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
           ( D, Em, Y, T )

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Scalar


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Vector &
    ( D, Ev, Ne, T, Em, Y )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(1:), INTENT(out) :: T, Em, Y

    INTEGER :: iP, nP, Error(SIZE(D))
    REAL(DP) :: D_P, E_P, Y_P

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: D, Ev, Ne ) &
    !$OMP MAP( from: T, Em, Y, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( D, Ev, Ne ) &
    !$ACC COPYOUT( T, Em, Y, Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      Em(iP) = Ev(iP) / D(iP)              ! --- Internal Energy per Mass
      Y (iP) = Ne(iP) / D(iP) * BaryonMass ! --- Electron Fraction

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
             ( D(iP), Em(iP), Y(iP), T(iP), Error_Option = Error(iP) )

    END DO

    IF ( ANY( Error > 0 ) ) THEN
      DO iP = 1, nP
        IF ( Error(iP) > 0 ) THEN
          CALL DescribeEOSInversionError( Error(iP) )
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM &
          !$OMP ( D(iP), Em(iP), Y(iP) )
#elif defined(THORNADO_OACC)
          !$ACC UPDATE HOST &
          !$ACC ( D(iP), Em(iP), Y(iP) )
#endif
          D_P = D(iP)  / ( Gram / Centimeter**3 )
          E_P = Em(iP) / ( Erg / Gram )
          Y_P = Y(iP)
          WRITE(*,*)                 '[ComputeThermodynamicStates_Auxiliary_TABLE_Vector] Error'
          WRITE(*,'(a,i5,3es23.15)') '  iP, D, E, Y : ', iP, D_P, E_P, Y_P
        END IF
      END DO
      STOP
    END IF

#endif

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Vector


  SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Scalar &
    ( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: P, T, Y, S, Em, Gm, Cs

    Em = Ev / D
    Y  = Ne * ( BaryonMass / D )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
           ( D, Em, Y, T )

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, P, P_T, OS_P, Units_V = UnitP )

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, S, S_T, OS_S, Units_V = UnitS )

    CALL ComputeDependentVariable_TABLE_Scalar &
           ( D, T, Y, Gm, Gm_T, OS_Gm, Units_V = UnitGm )

    Cs = SQRT( Gm * P / D )

  END SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Scalar


  SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Vector &
    ( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

    REAL(DP), INTENT(in)  :: D(1:), Ev(1:), Ne(1:)
    REAL(DP), INTENT(out) :: P(1:), T (1:), Y (1:), S(1:), Em(1:), Gm(1:), Cs(1:)

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputeAuxiliary_Fluid_TABLE_Scalar &
             ( D(iP), Ev(iP), Ne(iP), &
               P(iP), T (iP), Y (iP), S(iP), Em(iP), Gm(iP), Cs(iP) )

    END DO

  END SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Vector


  SUBROUTINE ComputePressure_TABLE_Scalar &
    ( D, T, Y, P, dPdD_Option, dPdT_Option, dPdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: P
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dPdD_Local, dPdT_Local, dPdY_Local
    REAL(DP), POINTER :: dPdD      , dPdT      , dPdY

    ComputeDerivatives &
      =      PRESENT( dPdD_Option ) &
        .OR. PRESENT( dPdT_Option ) &
        .OR. PRESENT( dPdY_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dPdD_Option ) )THEN
        dPdD => dPdD_Option
      ELSE
        dPdD => dPdD_Local
      END IF

      IF( PRESENT( dPdT_Option ) )THEN
        dPdT => dPdT_Option
      ELSE
        dPdT => dPdT_Local
      END IF

      IF( PRESENT( dPdY_Option ) )THEN
        dPdY => dPdY_Option
      ELSE
        dPdY => dPdY_Local
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Scalar &
             ( D, T, Y, P, dPdD, dPdT, dPdY, P_T, OS_P, Units_V = UnitP )

    ELSE

      CALL ComputeDependentVariable_TABLE_Scalar &
             ( D, T, Y, P, P_T, OS_P, Units_V = UnitP )

    END IF

  END SUBROUTINE ComputePressure_TABLE_Scalar


  SUBROUTINE ComputePressure_TABLE_Vector &
    ( D, T, Y, P, dPdD_Option, dPdT_Option, dPdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out)                   :: P(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dPdD_Local, dPdT_Local, dPdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dPdD      , dPdT      , dPdY

    ComputeDerivatives &
      =      PRESENT( dPdD_Option ) &
        .OR. PRESENT( dPdT_Option ) &
        .OR. PRESENT( dPdY_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )

      IF( PRESENT( dPdD_Option ) )THEN
        dPdD(1:nP) => dPdD_Option(:)
      ELSE
        dPdD(1:nP) => dPdD_Local(:)
      END IF

      IF( PRESENT( dPdT_Option ) )THEN
        dPdT(1:nP) => dPdT_Option(:)
      ELSE
        dPdT(1:nP) => dPdT_Local(:)
      END IF

      IF( PRESENT( dPdY_Option ) )THEN
        dPdY(1:nP) => dPdY_Option(:)
      ELSE
        dPdY(1:nP) => dPdY_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Vector &
             ( D, T, Y, P, dPdD, dPdT, dPdY, P_T, OS_P, Units_V = UnitP )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Y, P, P_T, OS_P, Units_V = UnitP )

    END IF

  END SUBROUTINE ComputePressure_TABLE_Vector


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Scalar &
    ( D, T, Y, E, dEdD_Option, dEdT_Option, dEdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: E
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dEdD_Local, dEdT_Local, dEdY_Local
    REAL(DP), POINTER :: dEdD,       dEdT,       dEdY

    ComputeDerivatives &
      =      PRESENT( dEdD_Option ) &
        .OR. PRESENT( dEdT_Option ) &
        .OR. PRESENT( dEdY_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dEdD_Option ) )THEN
        dEdD => dEdD_Option
      ELSE
        dEdD => dEdD_Local
      END IF

      IF( PRESENT( dEdT_Option ) )THEN
        dEdT => dEdT_Option
      ELSE
        dEdT => dEdT_Local
      END IF

      IF( PRESENT( dEdY_Option ) )THEN
        dEdY => dEdY_Option
      ELSE
        dEdY => dEdY_Local
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Scalar &
             ( D, T, Y, E, dEdD, dEdT, dEdY, E_T, OS_E, Units_V = UnitE )

    ELSE

      CALL ComputeDependentVariable_TABLE_Scalar &
             ( D, T, Y, E, E_T, OS_E, Units_V = UnitE )

    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Vector &
    ( D, T, Y, E, dEdD_Option, dEdT_Option, dEdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out)                   :: E(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dEdD_Local, dEdT_Local, dEdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dEdD      , dEdT      , dEdY

    ComputeDerivatives &
      =      PRESENT( dEdD_Option ) &
        .OR. PRESENT( dEdT_Option ) &
        .OR. PRESENT( dEdY_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )

      IF( PRESENT( dEdD_Option ) )THEN
        dEdD(1:nP) => dEdD_Option(:)
      ELSE
        dEdD(1:nP) => dEdD_Local(:)
      END IF

      IF( PRESENT( dEdT_Option ) )THEN
        dEdT(1:nP) => dEdT_Option(:)
      ELSE
        dEdT(1:nP) => dEdT_Local(:)
      END IF

      IF( PRESENT( dEdY_Option ) )THEN
        dEdY(1:nP) => dEdY_Option(:)
      ELSE
        dEdY(1:nP) => dEdY_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Vector &
             ( D, T, Y, E, dEdD, dEdT, dEdY, E_T, OS_E, Units_V = UnitE )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Y, E, E_T, OS_E, Units_V = UnitE )

    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Vector


  SUBROUTINE ComputeElectronChemicalPotential_TABLE_Scalar &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdY

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dMdD_Option ) )THEN
        dMdD => dMdD_Option
      ELSE
        dMdD => dMdD_Local
      END IF

      IF( PRESENT( dMdT_Option ) )THEN
        dMdT => dMdT_Option
      ELSE
        dMdT => dMdT_Local
      END IF

      IF( PRESENT( dMdY_Option ) )THEN
        dMdY => dMdY_Option
      ELSE
        dMdY => dMdY_Local
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Scalar &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Me_T, OS_Me, Units_V = UnitMe )

    ELSE

      CALL ComputeDependentVariable_TABLE_Scalar &
             ( D, T, Y, M, Me_T, OS_Me, Units_V = UnitMe )

    END IF

  END SUBROUTINE ComputeElectronChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeElectronChemicalPotential_TABLE_Vector &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out)                   :: M(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdY

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF ( ComputeDerivatives ) THEN

      nP = SIZE( D )

      IF( PRESENT( dMdD_Option ) )THEN
        dMdD(1:nP) => dMdD_Option(:)
      ELSE
        dMdD(1:nP) => dMdD_Local(:)
      END IF

      IF( PRESENT( dMdT_Option ) )THEN
        dMdT(1:nP) => dMdT_Option(:)
      ELSE
        dMdT(1:nP) => dMdT_Local(:)
      END IF

      IF( PRESENT( dMdY_Option ) )THEN
        dMdY(1:nP) => dMdY_Option(:)
      ELSE
        dMdY(1:nP) => dMdY_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Vector &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Me_T, OS_Me, Units_V = UnitMe )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Y, M, Me_T, OS_Me, Units_V = UnitMe )

    END IF

  END SUBROUTINE ComputeElectronChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeProtonChemicalPotential_TABLE_Scalar &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdY

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF ( ComputeDerivatives ) THEN

      IF( PRESENT( dMdD_Option ) )THEN
        dMdD => dMdD_Option
      ELSE
        dMdD => dMdD_Local
      END IF

      IF( PRESENT( dMdT_Option ) )THEN
        dMdT => dMdT_Option
      ELSE
        dMdT => dMdT_Local
      END IF

      IF( PRESENT( dMdY_Option ) )THEN
        dMdY => dMdY_Option
      ELSE
        dMdY => dMdY_Local
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Scalar &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mp_T, OS_Mp, Units_V = UnitMp )

    ELSE

      CALL ComputeDependentVariable_TABLE_Scalar &
             ( D, T, Y, M, Mp_T, OS_Mp, Units_V = UnitMp )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeProtonChemicalPotential_TABLE_Vector &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out)                   :: M(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdY

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )

      IF( PRESENT( dMdD_Option ) )THEN
        dMdD(1:nP) => dMdD_Option(:)
      ELSE
        dMdD(1:nP) => dMdD_Local(:)
      END IF

      IF( PRESENT( dMdT_Option ) )THEN
        dMdT(1:nP) => dMdT_Option(:)
      ELSE
        dMdT(1:nP) => dMdT_Local(:)
      END IF

      IF( PRESENT( dMdY_Option ) )THEN
        dMdY(1:nP) => dMdY_Option(:)
      ELSE
        dMdY(1:nP) => dMdY_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Vector &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mp_T, OS_Mp, Units_V = UnitMp )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Y, M, Mp_T, OS_Mp, Units_V = UnitMp )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Scalar &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdY

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dMdD_Option ) )THEN
        dMdD => dMdD_Option
      ELSE
        dMdD => dMdD_Local
      END IF

      IF( PRESENT( dMdT_Option ) )THEN
        dMdT => dMdT_Option
      ELSE
        dMdT => dMdT_Local
      END IF

      IF( PRESENT( dMdY_Option ) )THEN
        dMdY => dMdY_Option
      ELSE
        dMdY => dMdY_Local
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Scalar &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mn_T, OS_Mn, Units_V = UnitMn )

    ELSE

      CALL ComputeDependentVariable_TABLE_Scalar &
             ( D, T, Y, M, Mn_T, OS_Mn, Units_V = UnitMn )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Vector &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out)                   :: M(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdY

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )

      IF( PRESENT( dMdD_Option ) )THEN
        dMdD(1:nP) => dMdD_Option(:)
      ELSE
        dMdD(1:nP) => dMdD_Local(:)
      END IF

      IF( PRESENT( dMdT_Option ) )THEN
        dMdT(1:nP) => dMdT_Option(:)
      ELSE
        dMdT(1:nP) => dMdT_Local(:)
      END IF

      IF( PRESENT( dMdY_Option ) )THEN
        dMdY(1:nP) => dMdY_Option(:)
      ELSE
        dMdY(1:nP) => dMdY_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Vector &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mn_T, OS_Mn, Units_V = UnitMn )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Y, M, Mn_T, OS_Mn, Units_V = UnitMn )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeProtonMassFraction_TABLE_Scalar &
    ( D, T, Y, X, dXdD_Option, dXdT_Option, dXdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: X
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dXdD_Local, dXdT_Local, dXdY_Local
    REAL(DP), POINTER :: dXdD      , dXdT      , dXdY

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdY_Option )

    IF ( ComputeDerivatives ) THEN

      IF( PRESENT( dXdD_Option ) )THEN
        dXdD => dXdD_Option
      ELSE
        dXdD => dXdD_Local
      END IF

      IF( PRESENT( dXdT_Option ) )THEN
        dXdT => dXdT_Option
      ELSE
        dXdT => dXdT_Local
      END IF

      IF( PRESENT( dXdY_Option ) )THEN
        dXdY => dXdY_Option
      ELSE
        dXdY => dXdY_Local
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Scalar &
             ( D, T, Y, X, dXdD, dXdT, dXdY, Xp_T, OS_Xp, Units_V = UnitXp )

    ELSE

      CALL ComputeDependentVariable_TABLE_Scalar &
             ( D, T, Y, X, Xp_T, OS_Xp, Units_V = UnitXp )

    END IF

  END SUBROUTINE ComputeProtonMassFraction_TABLE_Scalar


  SUBROUTINE ComputeProtonMassFraction_TABLE_Vector &
    ( D, T, Y, X, dXdD_Option, dXdT_Option, dXdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out)                   :: X(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dXdD_Local, dXdT_Local, dXdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dXdD      , dXdT      , dXdY

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdY_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )

      IF( PRESENT( dXdD_Option ) )THEN
        dXdD(1:nP) => dXdD_Option(:)
      ELSE
        dXdD(1:nP) => dXdD_Local(:)
      END IF

      IF( PRESENT( dXdT_Option ) )THEN
        dXdT(1:nP) => dXdT_Option(:)
      ELSE
        dXdT(1:nP) => dXdT_Local(:)
      END IF

      IF( PRESENT( dXdY_Option ) )THEN
        dXdY(1:nP) => dXdY_Option(:)
      ELSE
        dXdY(1:nP) => dXdY_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Vector &
             ( D, T, Y, X, dXdD, dXdT, dXdY, Xp_T, OS_Xp, Units_V = UnitXp )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Y, X, Xp_T, OS_Xp, Units_V = UnitXp )

    END IF

  END SUBROUTINE ComputeProtonMassFraction_TABLE_Vector


  SUBROUTINE ComputeNeutronMassFraction_TABLE_Scalar &
    ( D, T, Y, X, dXdD_Option, dXdT_Option, dXdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: X
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dXdD_Local, dXdT_Local, dXdY_Local
    REAL(DP), POINTER :: dXdD      , dXdT      , dXdY

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdY_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dXdD_Option ) )THEN
        dXdD => dXdD_Option
      ELSE
        dXdD => dXdD_Local
      END IF

      IF( PRESENT( dXdT_Option ) )THEN
        dXdT => dXdT_Option
      ELSE
        dXdT => dXdT_Local
      END IF

      IF( PRESENT( dXdY_Option ) )THEN
        dXdY => dXdY_Option
      ELSE
        dXdY => dXdY_Local
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Scalar &
             ( D, T, Y, X, dXdD, dXdT, dXdY, Xn_T, OS_Xn, Units_V = UnitXn )

    ELSE

      CALL ComputeDependentVariable_TABLE_Scalar &
             ( D, T, Y, X, Xn_T, OS_Xn, Units_V = UnitXn )

    END IF

  END SUBROUTINE ComputeNeutronMassFraction_TABLE_Scalar


  SUBROUTINE ComputeNeutronMassFraction_TABLE_Vector &
    ( D, T, Y, X, dXdD_Option, dXdT_Option, dXdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out)                   :: X(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dXdD_Local, dXdT_Local, dXdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dXdD      , dXdT      , dXdY

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdY_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )

      IF( PRESENT( dXdD_Option ) )THEN
        dXdD(1:nP) => dXdD_Option(:)
      ELSE
        dXdD(1:nP) => dXdD_Local(:)
      END IF

      IF( PRESENT( dXdT_Option ) )THEN
        dXdT(1:nP) => dXdT_Option(:)
      ELSE
        dXdT(1:nP) => dXdT_Local(:)
      END IF

      IF( PRESENT( dXdY_Option ) )THEN
        dXdY(1:nP) => dXdY_Option(:)
      ELSE
        dXdY(1:nP) => dXdY_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivatives_TABLE_Vector &
             ( D, T, Y, X, dXdD, dXdT, dXdY, Xn_T, OS_Xn, Units_V = UnitXn )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Y, X, Xn_T, OS_Xn, Units_V = UnitXn )

    END IF

  END SUBROUTINE ComputeNeutronMassFraction_TABLE_Vector


  SUBROUTINE ComputeDependentVariable_TABLE_Scalar &
    ( D, T, Y, V, V_T, OS_V, Units_V )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: V
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Y_P, V_P

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Y_P = Y / UnitY

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Y_P, D_T, T_T, Y_T, OS_V, V_T, V_P )

    V = V_P * Units_V

#else

    V = Zero

#endif

  END SUBROUTINE ComputeDependentVariable_TABLE_Scalar


  SUBROUTINE ComputeDependentVariable_TABLE_Vector &
    ( D, T, Y, V, V_T, OS_V, Units_V )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out) :: V(1:)
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    INTEGER  :: iP, nP
    REAL(DP) :: D_P, T_P, Y_P, V_P

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( D_P, T_P, Y_P, V_P ) &
      !$OMP MAP( to: D, T, Y, D_T, T_T, Y_T, OS_V, V_T ) &
      !$OMP MAP( from: V )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( D_P, T_P, Y_P, V_P ) &
      !$ACC COPYIN( D, T, Y, D_T, T_T, Y_T, OS_V, V_T ) &
      !$ACC COPYOUT( V )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( D_P, T_P, Y_P, V_P )
#endif
    DO iP = 1, nP

      D_P = D(iP) / UnitD
      T_P = T(iP) / UnitT
      Y_P = Y(iP) / UnitY

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Y_P, D_T, T_T, Y_T, OS_V, V_T, V_P )

      V(iP) = V_P * Units_V

    END DO

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: V )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( V )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP
      V(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeDependentVariable_TABLE_Vector


  SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE_Scalar &
    ( D, T, Y, V, dVdD, dVdT, dVdY, V_T, OS_V, Units_V )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: V, dVdD, dVdT, dVdY
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Y_P, V_P, dV_P(3)

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Y_P = Y / UnitY

    CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Y_P, D_T, T_T, Y_T, OS_V, V_T, V_P, dV_P )

    V = V_P * Units_V

    dVdD = dV_P(1) * Units_V / UnitD
    dVdT = dV_P(2) * Units_V / UnitT
    dVdY = dV_P(3) * Units_V / UnitY

#else

    V    = Zero
    dVdD = Zero
    dVdT = Zero
    dVdY = Zero

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE_Vector &
    ( D, T, Y, V, dVdD, dVdT, dVdY, V_T, OS_V, Units_V )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Y(1:)
    REAL(DP), INTENT(out) :: V(1:), dVdD(1:), dVdT(1:), dVdY(1:)
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    INTEGER :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP MAP( to: D, T, Y, OS_V, V_T ) &
      !$OMP MAP( from: V, dVdD, dVdT, dVdY )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Y, OS_V, V_T ) &
      !$ACC COPYOUT( V, dVdD, dVdT, dVdY )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeDependentVariableAndDerivatives_TABLE_Scalar &
             ( D(iP), T(iP), Y(iP), V(iP), dVdD(iP), dVdT(iP), dVdY(iP), &
               V_T, OS_V, Units_V )

    END DO

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: V, dVdD, dVdT, dVdY )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( V, dVdD, dVdT, dVdY )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP
      V   (iP) = Zero
      dVdD(iP) = Zero
      dVdT(iP) = Zero
      dVdY(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE_Vector


END MODULE EquationOfStateModule_TABLE
