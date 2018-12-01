MODULE EquationOfStateModule

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule_IDEAL, ONLY: &
    InitializeEquationOfState_IDEAL, &
    ComputeInternalEnergyDensityFromPressure_IDEAL, &
    ComputePressureFromPrimitive_IDEAL, &
    ComputePressureFromSpecificInternalEnergy_IDEAL, &
    ComputeSoundSpeedFromPrimitive_IDEAL, &
    ComputeSoundSpeedFromPrimitive_GR_IDEAL, &
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
      REAL(DP), DIMENSION(:), INTENT(in)            :: X, Y, Z
      REAL(DP), DIMENSION(:), INTENT(out)           :: V
      REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dVdX
      REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dVdY
      REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dVdZ
    END SUBROUTINE EosSubroutine_1_3
  END INTERFACE

  ! ---
  ! --- Declaration of Equation of State Functions and Subroutines ---
  ! ---

  PROCEDURE (EosFunction),       POINTER, PUBLIC :: &
    Auxiliary_Fluid                              => NULL()
  PROCEDURE (EosSubroutine_1),   POINTER, PUBLIC :: &
    ComputeInternalEnergyDensityFromPressure     => NULL(), &
    ComputePressureFromPrimitive                 => NULL(), &
    ComputePressureFromSpecificInternalEnergy    => NULL(), &
    ComputeTemperatureFromPressure               => NULL(), &
    ComputeTemperatureFromSpecificInternalEnergy => NULL(), &
    ComputeSoundSpeedFromPrimitive               => NULL(), &
    ComputeSoundSpeedFromPrimitive_GR            => NULL()
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
        ComputePressureFromPrimitive &
          => ComputePressureFromPrimitive_IDEAL
        ComputePressureFromSpecificInternalEnergy &
          => ComputePressureFromSpecificInternalEnergy_IDEAL
        ComputeSoundSpeedFromPrimitive &
          => ComputeSoundSpeedFromPrimitive_IDEAL
        ComputeSoundSpeedFromPrimitive_GR &
          => ComputeSoundSpeedFromPrimitive_GR_IDEAL
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
        ComputeTemperatureFromPressure &
          => ComputeTemperatureFromPressure_TABLE
        ComputeTemperatureFromSpecificInternalEnergy &
          => ComputeTemperatureFromSpecificInternalEnergy_TABLE
        ComputeThermodynamicStates_Primitive &
          => ComputeThermodynamicStates_Primitive_TABLE
        ComputeThermodynamicStates_Auxiliary &
          => ComputeThermodynamicStates_Auxiliary_TABLE
        ComputePressureFromPrimitive &
          => ComputePressureFromPrimitive_TABLE
        ComputePressureFromSpecificInternalEnergy &
          => ComputePressureFromSpecificInternalEnergy_TABLE
        ComputeSoundSpeedFromPrimitive &
          => ComputeSoundSpeedFromPrimitive_TABLE
        ComputeAuxiliary_Fluid &
          => ComputeAuxiliary_Fluid_TABLE
        Auxiliary_Fluid &
          => Auxiliary_Fluid_TABLE
        ComputeSpecificInternalEnergy &
          => ComputeSpecificInternalEnergy_TABLE
        ComputeElectronChemicalPotential &
          => ComputeElectronChemicalPotential_TABLE
        ComputeProtonChemicalPotential &
          => ComputeProtonChemicalPotential_TABLE
        ComputeNeutronChemicalPotential &
          => ComputeNeutronChemicalPotential_TABLE

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
        NULLIFY( ComputePressureFromPrimitive )
        NULLIFY( ComputePressureFromSpecificInternalEnergy )
        NULLIFY( ComputeSoundSpeedFromPrimitive )
        NULLIFY( ComputeAuxiliary_Fluid )
        NULLIFY( Auxiliary_Fluid )

      CASE ( 'TABLE' )

        NULLIFY( ApplyEquationOfState )
        NULLIFY( ComputeTemperatureFromPressure )
        NULLIFY( ComputeTemperatureFromSpecificInternalEnergy )
        NULLIFY( ComputeThermodynamicStates_Primitive )
        NULLIFY( ComputeThermodynamicStates_Auxiliary )
        NULLIFY( ComputePressureFromPrimitive )
        NULLIFY( ComputePressureFromSpecificInternalEnergy )
        NULLIFY( ComputeSoundSpeedFromPrimitive )
        NULLIFY( ComputeAuxiliary_Fluid )
        NULLIFY( Auxiliary_Fluid )
        NULLIFY( ComputeSpecificInternalEnergy )
        NULLIFY( ComputeElectronChemicalPotential )
        NULLIFY( ComputeProtonChemicalPotential )
        NULLIFY( ComputeNeutronChemicalPotential )

    END SELECT

  END SUBROUTINE FinalizeEquationOfState


END MODULE EquationOfStateModule
