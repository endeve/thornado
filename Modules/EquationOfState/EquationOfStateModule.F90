MODULE EquationOfStateModule

  USE KindModule, ONLY: &
    DP, Zero
  USE EquationOfStateModule_IDEAL, ONLY: &
    InitializeEquationOfState_IDEAL, &
    ComputeInternalEnergyDensityFromPressure_IDEAL, &
    ComputePressureFromPrimitive_IDEAL, &
    ComputePressureFromSpecificInternalEnergy_IDEAL, &
    ComputeSoundSpeedFromPrimitive_IDEAL, &
    ComputeAuxiliary_Fluid_IDEAL
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    ApplyEquationOfState_TABLE, &
    ComputePressure_TABLE, &
    ComputeTemperatureFromPressure_TABLE, &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputePressureFromPrimitive_TABLE, &
    ComputePressureFromSpecificInternalEnergy_TABLE, &
    ComputeSoundSpeedFromPrimitive_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE
  USE UnitsModule, ONLY: &
    AtomicMassUnit

  IMPLICIT NONE
  PRIVATE

  CHARACTER(5), PUBLIC :: &
    EquationOfState
  REAL(DP), PUBLIC, PARAMETER :: &
    BaryonMass = AtomicMassUnit

  PUBLIC :: ApplyEquationOfState
  PUBLIC :: ComputePressure
  PUBLIC :: ComputeSpecificInternalEnergy
  PUBLIC :: ComputePressureFromPrimitive
  PUBLIC :: ComputePressureFromSpecificInternalEnergy
  PUBLIC :: ComputeInternalEnergyDensityFromPressure
  PUBLIC :: ComputeSoundSpeedFromPrimitive
  PUBLIC :: ComputeAuxiliary_Fluid
  PUBLIC :: ComputeThermodynamicStates_Primitive
  PUBLIC :: ComputeThermodynamicStates_Auxiliary
  PUBLIC :: ComputeTemperatureFromSpecificInternalEnergy
  PUBLIC :: ComputeTemperatureFromPressure
  PUBLIC :: ComputeElectronChemicalPotential
  PUBLIC :: ComputeProtonChemicalPotential
  PUBLIC :: ComputeNeutronChemicalPotential

  INTERFACE ApplyEquationOfState
    MODULE PROCEDURE ApplyEquationOfState_Scalar
    MODULE PROCEDURE ApplyEquationOfState_Vector
  END INTERFACE ApplyEquationOfState

  INTERFACE ComputePressure
    MODULE PROCEDURE ComputePressure_Scalar
    MODULE PROCEDURE ComputePressure_Vector
  END INTERFACE ComputePressure

  INTERFACE ComputeSpecificInternalEnergy
    MODULE PROCEDURE ComputeSpecificInternalEnergy_Scalar
    MODULE PROCEDURE ComputeSpecificInternalEnergy_Vector
  END INTERFACE ComputeSpecificInternalEnergy

  INTERFACE ComputePressureFromPrimitive
    MODULE PROCEDURE ComputePressureFromPrimitive_Scalar
    MODULE PROCEDURE ComputePressureFromPrimitive_Vector
  END INTERFACE ComputePressureFromPrimitive

  INTERFACE ComputePressureFromSpecificInternalEnergy
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_Scalar
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_Vector
  END INTERFACE ComputePressureFromSpecificInternalEnergy

  INTERFACE ComputeInternalEnergyDensityFromPressure
    MODULE PROCEDURE ComputeInternalEnergyDensityFromPressure_Scalar
    MODULE PROCEDURE ComputeInternalEnergyDensityFromPressure_Vector
  END INTERFACE ComputeInternalEnergyDensityFromPressure

  INTERFACE ComputeSoundSpeedFromPrimitive
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_Scalar
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_Vector
  END INTERFACE ComputeSoundSpeedFromPrimitive

  INTERFACE ComputeAuxiliary_Fluid
    MODULE PROCEDURE ComputeAuxiliary_Fluid_Scalar
    MODULE PROCEDURE ComputeAuxiliary_Fluid_Vector
  END INTERFACE ComputeAuxiliary_Fluid

  INTERFACE ComputeThermodynamicStates_Primitive
    MODULE PROCEDURE ComputeThermodynamicStates_Primitive_Scalar
    MODULE PROCEDURE ComputeThermodynamicStates_Primitive_Vector
  END INTERFACE ComputeThermodynamicStates_Primitive

  INTERFACE ComputeThermodynamicStates_Auxiliary
    MODULE PROCEDURE ComputeThermodynamicStates_Auxiliary_Scalar
    MODULE PROCEDURE ComputeThermodynamicStates_Auxiliary_Vector
  END INTERFACE ComputeThermodynamicStates_Auxiliary

  INTERFACE ComputeTemperatureFromSpecificInternalEnergy
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_Scalar
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_Vector
  END INTERFACE ComputeTemperatureFromSpecificInternalEnergy

  INTERFACE ComputeTemperatureFromPressure
    MODULE PROCEDURE ComputeTemperatureFromPressure_Scalar
    MODULE PROCEDURE ComputeTemperatureFromPressure_Vector
  END INTERFACE ComputeTemperatureFromPressure

  INTERFACE ComputeElectronChemicalPotential
    MODULE PROCEDURE ComputeElectronChemicalPotential_Scalar
    MODULE PROCEDURE ComputeElectronChemicalPotential_Vector
  END INTERFACE ComputeElectronChemicalPotential

  INTERFACE ComputeProtonChemicalPotential
    MODULE PROCEDURE ComputeProtonChemicalPotential_Scalar
    MODULE PROCEDURE ComputeProtonChemicalPotential_Vector
  END INTERFACE ComputeProtonChemicalPotential

  INTERFACE ComputeNeutronChemicalPotential
    MODULE PROCEDURE ComputeNeutronChemicalPotential_Scalar
    MODULE PROCEDURE ComputeNeutronChemicalPotential_Vector
  END INTERFACE ComputeNeutronChemicalPotential

  PUBLIC :: InitializeEquationOfState
  PUBLIC :: FinalizeEquationOfState

CONTAINS


  SUBROUTINE InitializeEquationOfState &
    ( EquationOfState_Option, EquationOfStateTableName_Option, &
      Gamma_IDEAL_Option, Verbose_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfState_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Gamma_IDEAL_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( EquationOfState_Option ) )THEN
      EquationOfState = EquationOfState_Option
    ELSE
      EquationOfState = 'IDEAL'
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A19,A)') &
        '', 'Equation Of State: ', TRIM( EquationOfState )
      WRITE(*,'(A5,A19)') &
        '', '------------------ '
    END IF

    SELECT CASE ( TRIM( EquationOfState ) )

      CASE ( 'IDEAL' )

        CALL InitializeEquationOfState_IDEAL &
               ( Gamma_IDEAL_Option &
                   = Gamma_IDEAL_Option, &
                 Verbose_Option = Verbose )

      CASE ( 'TABLE' )

        CALL InitializeEquationOfState_TABLE &
               ( EquationOfStateTableName_Option &
                   = EquationOfStateTableName_Option, &
                 Verbose_Option = Verbose )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A27,A5)') &
          ' ', 'Invalid Equation of State: ', EquationOfState
        STOP

    END SELECT

  END SUBROUTINE InitializeEquationOfState


  SUBROUTINE FinalizeEquationOfState

    SELECT CASE ( EquationOfState )

      CASE ( 'IDEAL' )

      CASE ( 'TABLE' )

    END SELECT

  END SUBROUTINE FinalizeEquationOfState


  ! --- ApplyEquationOfState ---


  SUBROUTINE ApplyEquationOfState_Scalar &
    ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm

#ifdef MICROPHYSICS_WEAKLIB

    CALL ApplyEquationOfState_TABLE &
           ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

#else

#endif

  END SUBROUTINE ApplyEquationOfState_Scalar


  SUBROUTINE ApplyEquationOfState_Vector &
    ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

    REAL(DP), INTENT(in)  :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out) :: P(:), S(:), E(:), Me(:), Mp(:), Mn(:)
    REAL(DP), INTENT(out) :: Xp(:), Xn(:), Xa(:), Xh(:), Gm(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ApplyEquationOfState_TABLE &
           ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

#else

#endif

  END SUBROUTINE ApplyEquationOfState_Vector


  ! --- ComputePressure ---


  SUBROUTINE ComputePressure_Scalar &
    ( D, T, Y, P, dPdD_Option, dPdT_Option, dPdY_Option )

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: P
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dPdD_Local, dPdT_Local, dPdY_Local
    REAL(DP), POINTER :: dPdD      , dPdT      , dPdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputePressure_TABLE &
             ( D, T, Y, P, dPdD, dPdT, dPdY )

    ELSE

      CALL ComputePressure_TABLE &
             ( D, T, Y, P )

    END IF

#else

#endif

  END SUBROUTINE ComputePressure_Scalar


  SUBROUTINE ComputePressure_Vector &
    ( D, T, Y, P, dPdD_Option, dPdT_Option, dPdY_Option )

    REAL(DP), INTENT(in)                    :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out)                   :: P(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdD_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdT_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdY_Option(:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dPdD_Local, dPdT_Local, dPdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dPdD      , dPdT      , dPdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputePressure_TABLE &
             ( D, T, Y, P, dPdD, dPdT, dPdY )

    ELSE

      CALL ComputePressure_TABLE &
             ( D, T, Y, P )

    END IF

#else

#endif

  END SUBROUTINE ComputePressure_Vector


  ! --- ComputeSpecificInternalEnergy ---


  SUBROUTINE ComputeSpecificInternalEnergy_Scalar &
    ( D, T, Y, E, dEdD_Option, dEdT_Option, dEdY_Option )

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: E
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dEdD_Local, dEdT_Local, dEdY_Local
    REAL(DP), POINTER :: dEdD      , dEdT      , dEdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( D, T, Y, E, dEdD, dEdT, dEdY )

    ELSE

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( D, T, Y, E )

    END IF

#else

#endif

  END SUBROUTINE ComputeSpecificInternalEnergy_Scalar


  SUBROUTINE ComputeSpecificInternalEnergy_Vector &
    ( D, T, Y, E, dEdD_Option, dEdT_Option, dEdY_Option )

    REAL(DP), INTENT(in)                    :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out)                   :: E(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdD_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdT_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdY_Option(:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dEdD_Local, dEdT_Local, dEdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dEdD      , dEdT      , dEdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( D, T, Y, E, dEdD, dEdT, dEdY )

    ELSE

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( D, T, Y, E )

    END IF

#else

#endif

  END SUBROUTINE ComputeSpecificInternalEnergy_Vector


  ! --- ComputePressureFromPrimitive ---


  SUBROUTINE ComputePressureFromPrimitive_Scalar &
    ( D, Ev, Ne, P )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: P

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputePressureFromPrimitive_TABLE( D, Ev, Ne, P )

#else

    CALL ComputePressureFromPrimitive_IDEAL( D, Ev, Ne, P )

#endif

  END SUBROUTINE ComputePressureFromPrimitive_Scalar


  SUBROUTINE ComputePressureFromPrimitive_Vector &
    ( D, Ev, Ne, P )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: P(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputePressureFromPrimitive_TABLE( D, Ev, Ne, P )

#else

    CALL ComputePressureFromPrimitive_IDEAL( D, Ev, Ne, P )

#endif

  END SUBROUTINE ComputePressureFromPrimitive_Vector


  ! --- ComputePressureFromSpecificInternalEnergy ---


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_Scalar &
    ( D, Em, Y, P )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Em, Y
    REAL(DP), INTENT(out) :: P

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputePressureFromSpecificInternalEnergy_TABLE( D, Em, Y, P )

#else

    CALL ComputePressureFromSpecificInternalEnergy_IDEAL( D, Em, Y, P )

#endif

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_Scalar


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_Vector &
    ( D, Em, Y, P )

    REAL(DP), INTENT(in)  :: D(:), Em(:), Y(:)
    REAL(DP), INTENT(out) :: P(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputePressureFromSpecificInternalEnergy_TABLE( D, Em, Y, P )

#else

    CALL ComputePressureFromSpecificInternalEnergy_IDEAL( D, Em, Y, P )

#endif

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_Vector


  ! --- ComputeInternalEnergyDensityFromPressure ---


  SUBROUTINE ComputeInternalEnergyDensityFromPressure_Scalar &
    ( D, P, Y, Ev )

    REAL(DP), INTENT(in)  :: D, P, Y
    REAL(DP), INTENT(out) :: Ev

#ifdef MICROPHYSICS_WEAKLIB

#else

    CALL ComputeInternalEnergyDensityFromPressure_IDEAL( D, P, Y, Ev )

#endif

  END SUBROUTINE ComputeInternalEnergyDensityFromPressure_Scalar


  SUBROUTINE ComputeInternalEnergyDensityFromPressure_Vector &
    ( D, P, Y, Ev )

    REAL(DP), INTENT(in)  :: D(:), P(:), Y(:)
    REAL(DP), INTENT(out) :: Ev(:)

#ifdef MICROPHYSICS_WEAKLIB

#else

    CALL ComputeInternalEnergyDensityFromPressure_IDEAL( D, P, Y, Ev )

#endif

  END SUBROUTINE ComputeInternalEnergyDensityFromPressure_Vector


  ! --- ComputeSoundSpeedFromPrimitive ---


  SUBROUTINE ComputeSoundSpeedFromPrimitive_Scalar &
    ( D, Ev, Ne, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: Cs

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeSoundSpeedFromPrimitive_TABLE( D, Ev, Ne, Cs )

#else

    CALL ComputeSoundSpeedFromPrimitive_IDEAL( D, Ev, Ne, Cs )

#endif

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_Scalar


  SUBROUTINE ComputeSoundSpeedFromPrimitive_Vector &
    ( D, Ev, Ne, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: Cs(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeSoundSpeedFromPrimitive_TABLE( D, Ev, Ne, Cs )

#else

    CALL ComputeSoundSpeedFromPrimitive_IDEAL( D, Ev, Ne, Cs )

#endif

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_Vector


  ! --- ComputeAuxiliary_Fluid ---


  SUBROUTINE ComputeAuxiliary_Fluid_Scalar &
    ( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: P, T, Y, S, Em, Gm, Cs

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeAuxiliary_Fluid_TABLE( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

#else

    CALL ComputeAuxiliary_Fluid_IDEAL( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

#endif

  END SUBROUTINE ComputeAuxiliary_Fluid_Scalar


  SUBROUTINE ComputeAuxiliary_Fluid_Vector &
    ( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: P(:), T (:), Y (:), S(:), Em(:), Gm(:), Cs(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeAuxiliary_Fluid_TABLE( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

#else

    CALL ComputeAuxiliary_Fluid_IDEAL( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

#endif

  END SUBROUTINE ComputeAuxiliary_Fluid_Vector


  ! --- ComputeThermodynamicStates_Primitive ---


  SUBROUTINE ComputeThermodynamicStates_Primitive_Scalar &
    ( D, T, Y, Ev, Em, Ne )

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: Ev, Em, Ne

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeThermodynamicStates_Primitive_TABLE( D, T, Y, Ev, Em, Ne )

#else

#endif

  END SUBROUTINE ComputeThermodynamicStates_Primitive_Scalar


  SUBROUTINE ComputeThermodynamicStates_Primitive_Vector &
    ( D, T, Y, Ev, Em, Ne )

    REAL(DP), INTENT(in)  :: D (:), T (:), Y (:)
    REAL(DP), INTENT(out) :: Ev(:), Em(:), Ne(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeThermodynamicStates_Primitive_TABLE( D, T, Y, Ev, Em, Ne )

#else

#endif

  END SUBROUTINE ComputeThermodynamicStates_Primitive_Vector


  ! --- ComputeThermodynamicStates_Auxiliary ---


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_Scalar &
    ( D, Ev, Ne, T, Em, Y )

    REAL(DP), INTENT(in)  :: D, Ev, Ne
    REAL(DP), INTENT(out) :: T, Em, Y

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeThermodynamicStates_Auxiliary_TABLE( D, Ev, Ne, T, Em, Y )

#else

#endif

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_Scalar


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_Vector &
    ( D, Ev, Ne, T, Em, Y )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: T(:), Em(:), Y (:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeThermodynamicStates_Auxiliary_TABLE( D, Ev, Ne, T, Em, Y )

#else

#endif

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_Vector


  ! --- ComputeTemperatureFromSpecificInternalEnergy ---


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_Scalar &
    ( D, E, Y, T, Guess_Option, Error_Option )

    REAL(DP), INTENT(in ) :: D, E, Y
    REAL(DP), INTENT(out) :: T
    REAL(DP), INTENT(in ), OPTIONAL :: Guess_Option
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    INTEGER  :: Error

#ifdef MICROPHYSICS_WEAKLIB

    IF( PRESENT( Guess_Option ) )THEN

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D, E, Y, T, Guess_Option, Error )

    ELSE

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D, E, Y, T, Error )

    END IF

    IF( PRESENT( Error_Option ) ) Error_Option = Error

#else

    T = Zero

#endif

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_Scalar


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_Vector &
    ( D, E, T, Y, Guess_Option, Error_Option )

    REAL(DP), INTENT(in ) :: D(:), E(:), Y(:)
    REAL(DP), INTENT(out) :: T(:)
    REAL(DP), INTENT(in ), OPTIONAL :: Guess_Option(:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(:)

    INTEGER :: Error(SIZE(D))

#ifdef MICROPHYSICS_WEAKLIB

    IF( PRESENT( Guess_Option ) )THEN

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D, E, Y, T, Guess_Option, Error )

    ELSE

       CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D, E, Y, T, Error_Option = Error )

    END IF

    IF( PRESENT( Error_Option ) ) Error_Option = Error

#else

    T = Zero

#endif

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_Vector


  ! --- ComputeTemperatureFromPressure ---


  SUBROUTINE ComputeTemperatureFromPressure_Scalar( D, P, Y, T )

    REAL(DP), INTENT(in)  :: D, P, Y
    REAL(DP), INTENT(out) :: T

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeTemperatureFromPressure_TABLE( D, P, Y, T )

#else

    T = Zero

#endif

  END SUBROUTINE ComputeTemperatureFromPressure_Scalar


  SUBROUTINE ComputeTemperatureFromPressure_Vector( D, P, Y, T )

    REAL(DP), INTENT(in)  :: D(:), P(:), Y(:)
    REAL(DP), INTENT(out) :: T(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeTemperatureFromPressure_TABLE( D, P, Y, T )

#else

    T = Zero

#endif

  END SUBROUTINE ComputeTemperatureFromPressure_Vector


  ! --- ComputeElectronChemicalPotential ---


  SUBROUTINE ComputeElectronChemicalPotential_Scalar &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputeElectronChemicalPotential_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY )

    ELSE

      CALL ComputeElectronChemicalPotential_TABLE &
             ( D, T, Y, M )

    END IF

#else

#endif

  END SUBROUTINE ComputeElectronChemicalPotential_Scalar


  SUBROUTINE ComputeElectronChemicalPotential_Vector &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out)                   :: M(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option(:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputeElectronChemicalPotential_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY )

    ELSE

      CALL ComputeElectronChemicalPotential_TABLE &
             ( D, T, Y, M )

    END IF

#else

#endif

  END SUBROUTINE ComputeElectronChemicalPotential_Vector


  ! --- ComputeProtonChemicalPotential ---


  SUBROUTINE ComputeProtonChemicalPotential_Scalar &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputeProtonChemicalPotential_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY )

    ELSE

      CALL ComputeProtonChemicalPotential_TABLE &
             ( D, T, Y, M )

    END IF

#else

#endif

  END SUBROUTINE ComputeProtonChemicalPotential_Scalar


  SUBROUTINE ComputeProtonChemicalPotential_Vector &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out)                   :: M(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option(:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputeProtonChemicalPotential_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY )

    ELSE

      CALL ComputeProtonChemicalPotential_TABLE &
             ( D, T, Y, M )

    END IF

#else

#endif

  END SUBROUTINE ComputeProtonChemicalPotential_Vector


  ! --- ComputeNeutronChemicalPotential ---


  SUBROUTINE ComputeNeutronChemicalPotential_Scalar &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputeNeutronChemicalPotential_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY )

    ELSE

      CALL ComputeNeutronChemicalPotential_TABLE &
             ( D, T, Y, M )

    END IF

#else

#endif

  END SUBROUTINE ComputeNeutronChemicalPotential_Scalar


  SUBROUTINE ComputeNeutronChemicalPotential_Vector &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out)                   :: M(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option(:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdY

#ifdef MICROPHYSICS_WEAKLIB

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

      CALL ComputeNeutronChemicalPotential_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY )

    ELSE

      CALL ComputeNeutronChemicalPotential_TABLE &
             ( D, T, Y, M )

    END IF

#else

#endif

  END SUBROUTINE ComputeNeutronChemicalPotential_Vector


END MODULE EquationOfStateModule
