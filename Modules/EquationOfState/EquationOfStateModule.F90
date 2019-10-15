MODULE EquationOfStateModule

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule_IDEAL, ONLY: &
    InitializeEquationOfState_IDEAL, &
    ComputeInternalEnergyDensityFromPressure_IDEAL, &
    ComputePressureFromPrimitive_IDEAL, &
    ComputePressureFromSpecificInternalEnergy_IDEAL, &
    ComputeSoundSpeedFromPrimitive_IDEAL, &
    ComputeAuxiliary_Fluid_IDEAL, &
    Auxiliary_Fluid_IDEAL
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    ApplyEquationOfState_TABLE, &
    ComputeTemperatureFromPressure_TABLE, &
    ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputePressureFromPrimitive_TABLE, &
    ComputePressureFromSpecificInternalEnergy_TABLE, &
    ComputeSoundSpeedFromPrimitive_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    Auxiliary_Fluid_TABLE, &
    ComputeSpecificInternalEnergy_Points_TABLE, &
    ComputeElectronChemicalPotentialPoints_TABLE, &
    ComputeProtonChemicalPotentialPoints_TABLE, &
    ComputeNeutronChemicalPotentialPoints_TABLE
  USE UnitsModule, ONLY: &
    AtomicMassUnit

  IMPLICIT NONE
  PRIVATE

  CHARACTER(5), PUBLIC :: &
    EquationOfState
  REAL(DP), PUBLIC, PARAMETER :: &
    BaryonMass = AtomicMassUnit

  PUBLIC :: ComputePressureFromPrimitive
  PUBLIC :: ComputeSoundSpeedFromPrimitive
  PUBLIC :: ComputeTemperatureFromSpecificInternalEnergy
  PUBLIC :: ComputeTemperatureFromPressure

  INTERFACE ComputePressureFromPrimitive
    MODULE PROCEDURE ComputePressureFromPrimitive_Scalar
    MODULE PROCEDURE ComputePressureFromPrimitive_Vector
  END INTERFACE ComputePressureFromPrimitive

  INTERFACE ComputeSoundSpeedFromPrimitive
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_Scalar
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_Vector
  END INTERFACE ComputeSoundSpeedFromPrimitive

  INTERFACE ComputeTemperatureFromSpecificInternalEnergy
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_Scalar
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_Vector
  END INTERFACE ComputeTemperatureFromSpecificInternalEnergy

  INTERFACE ComputeTemperatureFromPressure
    MODULE PROCEDURE ComputeTemperatureFromPressure_Scalar
    MODULE PROCEDURE ComputeTemperatureFromPressure_Vector
  END INTERFACE ComputeTemperatureFromPressure

  ! ---
  ! --- Interfaces for Various Equation of State Functions and Subroutines ---
  ! ---

  INTERFACE
    FUNCTION EosFunction( PF )
      USE KindModule, ONLY: DP
      USE FluidFieldsModule, ONLY: nPF, nAF
      REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
      REAL(DP), DIMENSION(nAF)             :: EosFunction
    END FUNCTION EosFunction
  END INTERFACE

  INTERFACE
    SUBROUTINE EosSubroutine_1( X, Y, Z, V1 )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)  :: X, Y, Z
      REAL(DP), DIMENSION(:), INTENT(out) :: V1
    END SUBROUTINE EosSubroutine_1
  END INTERFACE

  INTERFACE
    SUBROUTINE EosSubroutine_3( X, Y, Z, V1, V2, V3 )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)  :: X, Y, Z
      REAL(DP), DIMENSION(:), INTENT(out) :: V1, V2, V3
    END SUBROUTINE EosSubroutine_3
  END INTERFACE

  INTERFACE
    SUBROUTINE EosSubroutine_6( X, Y, Z, V1, V2, V3, V4, V5, V6 )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)  :: X, Y, Z
      REAL(DP), DIMENSION(:), INTENT(out) :: V1, V2, V3, V4, V5, V6
    END SUBROUTINE EosSubroutine_6
  END INTERFACE

  INTERFACE
    SUBROUTINE EosSubroutine_7( X, Y, Z, V1, V2, V3, V4, V5, V6, V7 )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)  :: X, Y, Z
      REAL(DP), DIMENSION(:), INTENT(out) :: V1, V2, V3, V4, V5, V6, V7
    END SUBROUTINE EosSubroutine_7
  END INTERFACE

  INTERFACE
    SUBROUTINE EosSubroutine_11 &
                 ( X, Y, Z, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11 )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)  :: X, Y, Z
      REAL(DP), DIMENSION(:), INTENT(out) :: V1, V2, V3, V4, V5, V6
      REAL(DP), DIMENSION(:), INTENT(out) :: V7, V8, V9, V10, V11
    END SUBROUTINE EosSubroutine_11
  END INTERFACE

  INTERFACE
    SUBROUTINE EosSubroutine_1_3( X, Y, Z, V, dVdX, dVdY, dVdZ )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)                    :: X, Y, Z
      REAL(DP), DIMENSION(:), INTENT(out)                   :: V
      REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dVdX
      REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dVdY
      REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dVdZ
    END SUBROUTINE EosSubroutine_1_3
  END INTERFACE

  INTERFACE
    SUBROUTINE EosSubroutine_1_2 &
                 ( X, Y, Z, V1, Guess_Option_V1, Error_Option )
      USE KindModule, ONLY: DP
        REAL(DP), INTENT(in)            :: X(:), Y(:), Z(:)
        REAL(DP), INTENT(out)           :: V1(:)
        REAL(DP), INTENT(in),  OPTIONAL :: Guess_Option_V1(:)
        INTEGER,  INTENT(out), OPTIONAL :: Error_Option(:)
    END SUBROUTINE EosSubroutine_1_2
  END INTERFACE

  ! ---
  ! --- Declaration of Equation of State Functions and Subroutines ---
  ! ---

  PROCEDURE (EosFunction),       POINTER, PUBLIC :: &
    Auxiliary_Fluid                              => NULL()
  PROCEDURE (EosSubroutine_1),   POINTER, PUBLIC :: &
    ComputeInternalEnergyDensityFromPressure     => NULL(), &
    ComputePressureFromSpecificInternalEnergy    => NULL()
  PROCEDURE (EosSubroutine_3),   POINTER, PUBLIC :: &
    ComputeThermodynamicStates_Primitive         => NULL(), &
    ComputeThermodynamicStates_Auxiliary         => NULL()
  PROCEDURE (EosSubroutine_7),   POINTER, PUBLIC :: &
    ComputeAuxiliary_Fluid                       => NULL()
  PROCEDURE (EosSubroutine_11),  POINTER, PUBLIC :: &
    ApplyEquationOfState                         => NULL()
  PROCEDURE (EosSubroutine_1_3), POINTER, PUBLIC :: &
    ComputeSpecificInternalEnergy                => NULL(), &
    ComputeElectronChemicalPotential             => NULL(), &
    ComputeProtonChemicalPotential               => NULL(), &
    ComputeNeutronChemicalPotential              => NULL()

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
                   = Gamma_IDEAL_Option )

        ComputeInternalEnergyDensityFromPressure &
          => ComputeInternalEnergyDensityFromPressure_IDEAL
        ComputePressureFromSpecificInternalEnergy &
          => ComputePressureFromSpecificInternalEnergy_IDEAL
        ComputeAuxiliary_Fluid &
          => ComputeAuxiliary_Fluid_IDEAL
        Auxiliary_Fluid &
          => Auxiliary_Fluid_IDEAL

      CASE ( 'TABLE' )

        CALL InitializeEquationOfState_TABLE &
               ( EquationOfStateTableName_Option &
                   = EquationOfStateTableName_Option, &
                 Verbose_Option = .TRUE. )

        ApplyEquationOfState &
          => ApplyEquationOfState_TABLE
        ComputeThermodynamicStates_Primitive &
          => ComputeThermodynamicStates_Primitive_TABLE
        ComputeThermodynamicStates_Auxiliary &
          => ComputeThermodynamicStates_Auxiliary_TABLE
        ComputePressureFromSpecificInternalEnergy &
          => ComputePressureFromSpecificInternalEnergy_TABLE
        ComputeAuxiliary_Fluid &
          => ComputeAuxiliary_Fluid_TABLE
        Auxiliary_Fluid &
          => Auxiliary_Fluid_TABLE
        ComputeSpecificInternalEnergy &
          => ComputeSpecificInternalEnergy_Points_TABLE
        ComputeElectronChemicalPotential &
          => ComputeElectronChemicalPotentialPoints_TABLE
        ComputeProtonChemicalPotential &
          => ComputeProtonChemicalPotentialPoints_TABLE
        ComputeNeutronChemicalPotential &
          => ComputeNeutronChemicalPotentialPoints_TABLE

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

        NULLIFY( ComputeInternalEnergyDensityFromPressure )
        NULLIFY( ComputePressureFromSpecificInternalEnergy )
        NULLIFY( ComputeAuxiliary_Fluid )
        NULLIFY( Auxiliary_Fluid )

      CASE ( 'TABLE' )

        NULLIFY( ApplyEquationOfState )
        NULLIFY( ComputeThermodynamicStates_Primitive )
        NULLIFY( ComputeThermodynamicStates_Auxiliary )
        NULLIFY( ComputePressureFromSpecificInternalEnergy )
        NULLIFY( ComputeAuxiliary_Fluid )
        NULLIFY( Auxiliary_Fluid )
        NULLIFY( ComputeSpecificInternalEnergy )
        NULLIFY( ComputeElectronChemicalPotential )
        NULLIFY( ComputeProtonChemicalPotential )
        NULLIFY( ComputeNeutronChemicalPotential )

    END SELECT

  END SUBROUTINE FinalizeEquationOfState


  ! --- ComputePressureFromPrimitive ---


  SUBROUTINE ComputePressureFromPrimitive_Scalar &
    ( D, Ev, Ne, P )

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


  ! --- ComputeSoundSpeedFromPrimitive ---


  SUBROUTINE ComputeSoundSpeedFromPrimitive_Scalar &
    ( D, Ev, Ne, Cs )

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
             ( D, E, Y, T, Error_Option = Error )

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


END MODULE EquationOfStateModule
