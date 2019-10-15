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

  PUBLIC :: ApplyEquationOfState
  PUBLIC :: ComputePressureFromPrimitive
  PUBLIC :: ComputePressureFromSpecificInternalEnergy
  PUBLIC :: ComputeSoundSpeedFromPrimitive
  PUBLIC :: ComputeThermodynamicStates_Primitive
  PUBLIC :: ComputeTemperatureFromSpecificInternalEnergy
  PUBLIC :: ComputeTemperatureFromPressure

  INTERFACE ApplyEquationOfState
    MODULE PROCEDURE ApplyEquationOfState_Scalar
    MODULE PROCEDURE ApplyEquationOfState_Vector
  END INTERFACE ApplyEquationOfState

  INTERFACE ComputePressureFromPrimitive
    MODULE PROCEDURE ComputePressureFromPrimitive_Scalar
    MODULE PROCEDURE ComputePressureFromPrimitive_Vector
  END INTERFACE ComputePressureFromPrimitive

  INTERFACE ComputePressureFromSpecificInternalEnergy
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_Scalar
    MODULE PROCEDURE ComputePressureFromSpecificInternalEnergy_Vector
  END INTERFACE ComputePressureFromSpecificInternalEnergy

  INTERFACE ComputeSoundSpeedFromPrimitive
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_Scalar
    MODULE PROCEDURE ComputeSoundSpeedFromPrimitive_Vector
  END INTERFACE ComputeSoundSpeedFromPrimitive

  INTERFACE ComputeThermodynamicStates_Primitive
    MODULE PROCEDURE ComputeThermodynamicStates_Primitive_Scalar
    MODULE PROCEDURE ComputeThermodynamicStates_Primitive_Vector
  END INTERFACE ComputeThermodynamicStates_Primitive

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
    ComputeInternalEnergyDensityFromPressure     => NULL()
  PROCEDURE (EosSubroutine_3),   POINTER, PUBLIC :: &
    ComputeThermodynamicStates_Auxiliary         => NULL()
  PROCEDURE (EosSubroutine_7),   POINTER, PUBLIC :: &
    ComputeAuxiliary_Fluid                       => NULL()
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
        ComputeAuxiliary_Fluid &
          => ComputeAuxiliary_Fluid_IDEAL
        Auxiliary_Fluid &
          => Auxiliary_Fluid_IDEAL

      CASE ( 'TABLE' )

        CALL InitializeEquationOfState_TABLE &
               ( EquationOfStateTableName_Option &
                   = EquationOfStateTableName_Option, &
                 Verbose_Option = .TRUE. )

        ComputeThermodynamicStates_Auxiliary &
          => ComputeThermodynamicStates_Auxiliary_TABLE
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
        NULLIFY( ComputeAuxiliary_Fluid )
        NULLIFY( Auxiliary_Fluid )

      CASE ( 'TABLE' )

        NULLIFY( ComputeThermodynamicStates_Auxiliary )
        NULLIFY( ComputeAuxiliary_Fluid )
        NULLIFY( Auxiliary_Fluid )
        NULLIFY( ComputeSpecificInternalEnergy )
        NULLIFY( ComputeElectronChemicalPotential )
        NULLIFY( ComputeProtonChemicalPotential )
        NULLIFY( ComputeNeutronChemicalPotential )

    END SELECT

  END SUBROUTINE FinalizeEquationOfState


  ! --- ApplyEquationOfState ---


  SUBROUTINE ApplyEquationOfState_Scalar &
    ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

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


  ! --- ComputePressureFromSpecificInternalEnergy ---


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_Scalar &
    ( D, Em, Ne, P )

    REAL(DP), INTENT(in)  :: D, Em, Ne
    REAL(DP), INTENT(out) :: P

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputePressureFromSpecificInternalEnergy_TABLE( D, Em, Ne, P )

#else

    CALL ComputePressureFromSpecificInternalEnergy_IDEAL( D, Em, Ne, P )

#endif

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_Scalar


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_Vector &
    ( D, Em, Ne, P )

    REAL(DP), INTENT(in)  :: D(:), Em(:), Ne(:)
    REAL(DP), INTENT(out) :: P(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputePressureFromSpecificInternalEnergy_TABLE( D, Em, Ne, P )

#else

    CALL ComputePressureFromSpecificInternalEnergy_IDEAL( D, Em, Ne, P )

#endif

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_Vector


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

    REAL(DP), INTENT(in)  :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out) :: Ev(:), Em(:), Ne(:)

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeThermodynamicStates_Primitive_TABLE( D, T, Y, Ev, Em, Ne )

#else

#endif

  END SUBROUTINE ComputeThermodynamicStates_Primitive_Vector


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
