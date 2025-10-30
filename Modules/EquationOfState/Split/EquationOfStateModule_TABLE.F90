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
    EquationOfStateCompOSETableType
#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE
  USE wlEOSComponentsSeparateInversionModule, ONLY: &
    InitializeEOSComponentsInversion, &
    ComputeTemperatureWith_DEYpYl_Single_Guess_Error, &
    ComputeTemperatureWith_DEYpYl_Single_Guess_NoError, &
    ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error, &
    ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError, &
    ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error, &
    DescribeEOSComponentsInversionError
#else
  USE wlEOSComponentsCombinedInversionModule, ONLY: &
    InitializeEOSComponentsInversion, &
    ComputeTemperatureWith_DEYpYl_Single_Guess_Error, &
    ComputeTemperatureWith_DEYpYl_Single_Guess_NoError, &
    ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error, &
    ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError, &
    ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error, &
    DescribeEOSComponentsInversionError
#endif
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_3D_Custom_Point, &
    LogInterpolateDifferentiateSingleVariable_3D_Custom_Point
  USE wlLeptonEOSModule, ONLY: &
    HelmTableType, MuonTableType
  USE wlElectronPhotonEOS, ONLY: &
    ElectronPhotonEOS, ElectronPhotonStateType
  USE wlMuonEOS, ONLY: &
    FullMuonEOS, MuonStateType
  USE wlHelmMuonIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF, ReadMuonTableHDF
  USE wlSoundSpeedModule, ONLY: &
    CalculateSoundSpeed
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, Index1D_Log, &
    GetIndexAndDelta_Log, GetIndexAndDelta_Lin, &
    LinearInterpDeriv_Array_Point
    
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
    Second, &
    Kelvin, &
    Dyne, &
    Erg, &
    MeV

  IMPLICIT NONE
  PRIVATE

  CHARACTER(256) :: &
    EquationOfStateTableName
  INTEGER :: &
    iD_T, iT_T, iYp_T, &
    iP_T, iS_T, iE_T, iMe_T, iMp_T, iMn_T, &
    iXp_T, iXn_T, iXa_T, iXh_T, iAh_T, iGm_T, &
    iEmp_T, iEmn_T, iSep_T, iSen_T
  REAL(DP) :: &
    UnitD, UnitT, UnitY, &
    UnitP, UnitS, UnitE, UnitMl, UnitMp, UnitMn, &
    UnitXp, UnitXn, UnitXa, UnitXh, UnitAh, UnitGm, &
    UnitEmp, UnitEmn, UnitSep, UnitSen
  REAL(DP), PUBLIC :: &
    OS_P, OS_S, OS_E, OS_Mp, OS_Mn, &
    OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Ah, OS_Gm, &
    OS_Emp, OS_Emn, OS_Sep, OS_Sen
  REAL(DP), PARAMETER :: &
    BaryonMass = AtomicMassUnit
  REAL(DP) :: minvar, OS_loc
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    D_T, T_T, Yp_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    P_T, S_T, E_T, Mp_T, Mn_T, &
    Xp_T, Xn_T, Xa_T, Xh_T, Ah_T, Gm_T, &
    Emp_T, Emn_T, Sep_T, Sen_T
#ifdef MICROPHYSICS_WEAKLIB
  LOGICAL :: UsingExternalEOS
  TYPE(EquationOfStateCompOSETableType), POINTER :: EOS
  TYPE(HelmTableType), POINTER :: HelmTable
  TYPE(MuonTableType), POINTER :: MuonTable
#endif
  REAL(DP), PUBLIC :: Eos_MinD ! This is not handled consistently everywhere, 
                               ! but should not matter

  PUBLIC :: InitializeEquationOfState_TABLE
  PUBLIC :: FinalizeEquationOfState_TABLE
  PUBLIC :: ApplyChemicalPotentialShift_TABLE
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
  PUBLIC :: ComputeMuonChemicalPotential_TABLE
  PUBLIC :: ComputeProtonChemicalPotential_TABLE
  PUBLIC :: ComputeNeutronChemicalPotential_TABLE
  PUBLIC :: ComputeProtonMassFraction_TABLE
  PUBLIC :: ComputeNeutronMassFraction_TABLE
  PUBLIC :: ComputeHeavyMassFraction_TABLE
  PUBLIC :: ComputeHeavyMassNumber_TABLE
  PUBLIC :: ComputeElectronNeutrinoChemicalPotential_TABLE
  PUBLIC :: ComputeMuonNeutrinoChemicalPotential_TABLE

  ! These are not present in the old weaklib table
  PUBLIC :: ComputeProtonEffectiveMass_TABLE
  PUBLIC :: ComputeNeutronEffectiveMass_TABLE
  PUBLIC :: ComputeProtonSelfEnergy_TABLE
  PUBLIC :: ComputeNeutronSelfEnergy_TABLE

  REAL(DP), PUBLIC :: Min_D, Min_T, Min_Y
  REAL(DP), PUBLIC :: Max_D, Max_T, Max_Y

  INTERFACE ApplyEquationOfState_TABLE
    MODULE PROCEDURE ApplyEquationOfState_TABLE_Scalar
    MODULE PROCEDURE ApplyEquationOfState_TABLE_Vector
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
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_NG_NE
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_NG_E
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_G_NE
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_G_E
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector
  END INTERFACE ComputeTemperatureFromSpecificInternalEnergy_TABLE

  INTERFACE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_NG_NE
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_NG_E
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_G_NE
    MODULE PROCEDURE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_G_E
  END INTERFACE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar

  INTERFACE ComputeTemperatureFromPressure_TABLE
    MODULE PROCEDURE ComputeTemperatureFromPressure_TABLE_Scalar
    MODULE PROCEDURE ComputeTemperatureFromPressure_TABLE_Vector
  END INTERFACE ComputeTemperatureFromPressure_TABLE

  INTERFACE ComputeElectronChemicalPotential_TABLE
    MODULE PROCEDURE ComputeElectronChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeElectronChemicalPotential_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeMuonChemicalPotential_TABLE
    MODULE PROCEDURE ComputeMuonChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeMuonChemicalPotential_TABLE_Vector
  END INTERFACE
  
  INTERFACE ComputeProtonChemicalPotential_TABLE
    MODULE PROCEDURE ComputeProtonChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeProtonChemicalPotential_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeNeutronChemicalPotential_TABLE
    MODULE PROCEDURE ComputeNeutronChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeNeutronChemicalPotential_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeElectronNeutrinoChemicalPotential_TABLE
    MODULE PROCEDURE ComputeElectronNeutrinoChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeElectronNeutrinoChemicalPotential_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeMuonNeutrinoChemicalPotential_TABLE
    MODULE PROCEDURE ComputeMuonNeutrinoChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeMuonNeutrinoChemicalPotential_TABLE_Vector
  END INTERFACE
  
  INTERFACE ComputeProtonMassFraction_TABLE
    MODULE PROCEDURE ComputeProtonMassFraction_TABLE_Scalar
    MODULE PROCEDURE ComputeProtonMassFraction_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeNeutronMassFraction_TABLE
    MODULE PROCEDURE ComputeNeutronMassFraction_TABLE_Scalar
    MODULE PROCEDURE ComputeNeutronMassFraction_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeHeavyMassFraction_TABLE
      MODULE PROCEDURE ComputeHeavyMassFraction_TABLE_Scalar
    MODULE PROCEDURE ComputeHeavyMassFraction_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeHeavyMassNumber_TABLE
    MODULE PROCEDURE ComputeHeavyMassNumber_TABLE_Scalar
    MODULE PROCEDURE ComputeHeavyMassNumber_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeDependentVariableBaryons
    MODULE PROCEDURE ComputeDependentVariableBaryons_TABLE_Scalar
    MODULE PROCEDURE ComputeDependentVariableBaryons_TABLE_Vector
  END INTERFACE ComputeDependentVariableBaryons
  
  INTERFACE ComputeDependentVariableAndDerivativesBaryons_TABLE
    MODULE PROCEDURE ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar
    MODULE PROCEDURE ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector
  END INTERFACE ComputeDependentVariableAndDerivativesBaryons_TABLE

  INTERFACE ComputeDependentVariableTotal
    MODULE PROCEDURE ComputeDependentVariableTotal_TABLE_Scalar
    MODULE PROCEDURE ComputeDependentVariableTotal_TABLE_Vector
  END INTERFACE ComputeDependentVariableTotal

  INTERFACE ComputeTotalPES_TABLE
    MODULE PROCEDURE ComputeTotalPES_TABLE_Scalar
    MODULE PROCEDURE ComputeTotalPES_TABLE_Vector
  END INTERFACE ComputeTotalPES_TABLE

  INTERFACE ComputeDependentVariableAndDerivativesTotal_TABLE
    MODULE PROCEDURE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar
    MODULE PROCEDURE ComputeDependentVariableAndDerivativesTotal_TABLE_Vector
  END INTERFACE ComputeDependentVariableAndDerivativesTotal_TABLE

  INTERFACE ComputeProtonEffectiveMass_TABLE
    MODULE PROCEDURE ComputeProtonEffectiveMass_TABLE_Scalar
    MODULE PROCEDURE ComputeProtonEffectiveMass_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeNeutronEffectiveMass_TABLE
    MODULE PROCEDURE ComputeNeutronEffectiveMass_TABLE_Scalar
    MODULE PROCEDURE ComputeNeutronEffectiveMass_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeProtonSelfEnergy_TABLE
    MODULE PROCEDURE ComputeProtonSelfEnergy_TABLE_Scalar
    MODULE PROCEDURE ComputeProtonSelfEnergy_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeNeutronSelfEnergy_TABLE
    MODULE PROCEDURE ComputeNeutronSelfEnergy_TABLE_Scalar
    MODULE PROCEDURE ComputeNeutronSelfEnergy_TABLE_Vector
  END INTERFACE
  ! Define local constants

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( D_T, T_T, Yp_T, &
  !$OMP   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, &
  !$OMP   UnitMl, UnitMp, UnitMn, UnitXp, UnitXn, UnitXa, UnitXh, &
  !$OMP   UnitAh, UnitGm, UnitEmp, UnitEmn, UnitSep, UnitSen, &
  !$OMP   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$OMP   OS_Xa, OS_Xh, OS_Ah, OS_Gm, OS_Emp, OS_Emn, OS_Sep, OS_Sen, &
  !$OMP   P_T, S_T, E_T, Mp_T, Mn_T, Xp_T, Xn_T, &
  !$OMP   Xa_T, Xh_T, Ah_T, Gm_T, Emp_T, Emn_T, Sep_T, Sen_T, &
  !$OMP   HelmTable, MuonTable, &
  !$OMP   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y, Eos_MinD )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( D_T, T_T, Yp_T, &
  !$ACC   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, &
  !$ACC   UnitMl, UnitMp, UnitMn, UnitXp, UnitXn, UnitXa, UnitXh, &
  !$ACC   UnitAh, UnitGm, UnitEmp, UnitEmn, UnitSep, UnitSen, &
  !$ACC   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$ACC   OS_Xa, OS_Xh, OS_Ah, OS_Gm, OS_Emp, OS_Emn, OS_Sep, OS_Sen, &
  !$ACC   P_T, S_T, E_T, Mp_T, Mn_T, Xp_T, Xn_T, &
  !$ACC   Xa_T, Xh_T, Ah_T, Gm_T, Emp_T, Emn_T, Sep_T, Sen_T, &
  !$ACC   HelmTable, MuonTable, &
  !$ACC   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y, Eos_MinD )
#endif

CONTAINS

  SUBROUTINE InitializeEquationOfState_TABLE &
    ( EquationOfStateTableName_Option, UseChemicalPotentialShift_Option, &
      Eos_MinD_Option , Verbose_Option, External_EOS, External_Helm, External_Muon )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    LOGICAL,          INTENT(in), OPTIONAL :: UseChemicalPotentialShift_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Eos_MinD_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option
#ifdef MICROPHYSICS_WEAKLIB
    ! Notice that you use EquationOfState/Split only if EOSMODE == COMPOSE
    TYPE(EquationOfStateCompOSETableType), POINTER, &
                      INTENT(in), OPTIONAL :: External_EOS
    type(HelmTableType), pointer, &
                      intent(in), optional :: External_Helm
    type(MuonTableType), pointer, &
                      intent(in), optional :: External_Muon
#else
    INTEGER,          INTENT(in), OPTIONAL :: External_EOS
    INTEGER,          INTENT(in), OPTIONAL :: External_Helm
    INTEGER,          INTENT(in), OPTIONAL :: External_Muon
#endif
    
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Ps, Ss, Es
    LOGICAL :: UseChemicalPotentialShift, Verbose
    
    REAL(DP) :: Eele, Pele, Sele

    ! Added variables to handle Combined EOSs
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState

    IF( PRESENT( EquationOfStateTableName_Option ) )THEN
       EquationOfStateTableName = TRIM( EquationOfStateTableName_Option )
    ELSE
       EquationOfStateTableName = 'EquationOfStateTable.h5'
    END IF

    IF( PRESENT( UseChemicalPotentialShift_Option ) )THEN
       UseChemicalPotentialShift = UseChemicalPotentialShift_Option
    ELSE
       UseChemicalPotentialShift = .FALSE.
    END IF

    IF( PRESENT( Eos_MinD_Option ) )THEN
       Eos_MinD = Eos_MinD_Option
    ELSE
       Eos_MinD = Zero
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
       ALLOCATE( HelmTable)
       ALLOCATE( MuonTable )
       UsingExternalEOS = .FALSE.

       CALL InitializeHDF( )
        
       ! read in helmholtz table
       CALL ReadHelmholtzTableHDF( HelmTable, TRIM( EquationOfStateTableName ) )

       ! read in helmholtz table
       CALL ReadMuonTableHDF( MuonTable, TRIM( EquationOfStateTableName ) )

       ! read in baryon table -------------------------------
       CALL ReadEquationOfStateTableHDF( EOS, TRIM( EquationOfStateTableName ) )
        
       CALL FinalizeHDF( )

    ELSE

       EOS => External_EOS
       HelmTable => External_Helm
       MuonTable => External_Muon
       UsingExternalEOS = .TRUE.

    END IF

    ! --- Thermodynamic State Indices ---

    iD_T  = EOS % TS % Indices % iRho
    iT_T  = EOS % TS % Indices % iT
    iYp_T = EOS % TS % Indices % iYe

    ! --- Units ---

    UnitD = Gram / Centimeter**3
    UnitT = Kelvin
    UnitY = One

    UnitP  = Dyne / Centimeter**2
    UnitS  = BoltzmannConstant
    UnitE  = Erg / Gram
    UnitMl = MeV
    UnitMp = MeV
    UnitMn = MeV
    UnitXp = One
    UnitXn = One
    UnitXa = One
    UnitXh = One
    UnitAh = One
    UnitGm = One
    UnitEmp = One
    UnitEmn = One
    UnitSep = One
    UnitSen = One

    ! --- Thermodynamic States ---

    ALLOCATE( D_T(EOS % TS % nPoints(iD_T)) )
    D_T = EOS % TS % States(iD_T) % Values

    Min_D = MINVAL( D_T ) * Gram / Centimeter**3
    Max_D = MAXVAL( D_T ) * Gram / Centimeter**3

    ALLOCATE( T_T(EOS % TS % nPoints(iT_T)) )
    T_T = EOS % TS % States(iT_T) % Values

    Min_T = MINVAL( T_T ) * Kelvin
    Max_T = MAXVAL( T_T ) * Kelvin

    ALLOCATE( Yp_T(EOS % TS % nPoints(iYp_T)) )
    Yp_T = EOS % TS % States(iYp_T) % Values

    Min_Y = MINVAL( Yp_T )
    Max_Y = MAXVAL( Yp_T )

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
    iAh_T = EOS % DV % Indices % iHeavyMassNumber
    iGm_T = EOS % DV % Indices % iGamma1
    iEmp_T = EOS % DV % Indices % iProtonEffMass
    iEmn_T = EOS % DV % Indices % iNeutronEffMass
    iSep_T = EOS % DV % Indices % iProtonSelfEnergy
    iSen_T = EOS % DV % Indices % iNeutronSelfEnergy

    ! --- Dependent Variables Offsets ---

    OS_P  = EOS % DV % Offsets(iP_T)
    OS_S  = EOS % DV % Offsets(iS_T)
    OS_E  = EOS % DV % Offsets(iE_T)
    OS_Mp = EOS % DV % Offsets(iMp_T)
    OS_Mn = EOS % DV % Offsets(iMn_T)
    OS_Xp = EOS % DV % Offsets(iXp_T)
    OS_Xn = EOS % DV % Offsets(iXn_T)
    OS_Xa = EOS % DV % Offsets(iXa_T)
    OS_Xh = EOS % DV % Offsets(iXh_T)
    OS_Ah = EOS % DV % Offsets(iAh_T)
    OS_Gm = EOS % DV % Offsets(iGm_T)
    OS_Emp = EOS % DV % Offsets(iEmp_T)
    OS_Emn = EOS % DV % Offsets(iEmn_T)
    OS_Sep = EOS % DV % Offsets(iSep_T)
    OS_Sen = EOS % DV % Offsets(iSen_T)

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
      ( Ah_T(1:EOS % DV % nPoints(1), &
            1:EOS % DV % nPoints(2), &
            1:EOS % DV % nPoints(3)) )

    ALLOCATE &
      ( Gm_T(1:EOS % DV % nPoints(1), &
            1:EOS % DV % nPoints(2), &
            1:EOS % DV % nPoints(3)) )

    ALLOCATE &
      ( Emp_T(1:EOS % DV % nPoints(1), &
            1:EOS % DV % nPoints(2), &
            1:EOS % DV % nPoints(3)) )

    ALLOCATE &
      ( Emn_T(1:EOS % DV % nPoints(1), &
            1:EOS % DV % nPoints(2), &
            1:EOS % DV % nPoints(3)) ) 

    ALLOCATE &
      ( Sep_T(1:EOS % DV % nPoints(1), &
            1:EOS % DV % nPoints(2), &
            1:EOS % DV % nPoints(3)) ) 

    ALLOCATE &
      ( Sen_T(1:EOS % DV % nPoints(1), &
            1:EOS % DV % nPoints(2), &
            1:EOS % DV % nPoints(3)) ) 

    P_T  = EOS % DV % Variables(iP_T ) % Values
    S_T  = EOS % DV % Variables(iS_T ) % Values
    E_T  = EOS % DV % Variables(iE_T ) % Values
    Mp_T = EOS % DV % Variables(iMp_T) % Values
    Mn_T = EOS % DV % Variables(iMn_T) % Values
    Xp_T = EOS % DV % Variables(iXp_T) % Values
    Xn_T = EOS % DV % Variables(iXn_T) % Values
    Xa_T = EOS % DV % Variables(iXa_T) % Values
    Xh_T = EOS % DV % Variables(iXh_T) % Values
    Ah_T = EOS % DV % Variables(iAh_T) % Values
    Gm_T = EOS % DV % Variables(iGm_T) % Values
    Emp_T = EOS % DV % Variables(iEmp_T) % Values
    Emn_T = EOS % DV % Variables(iEmn_T) % Values
    Sep_T = EOS % DV % Variables(iSep_T) % Values
    Sen_T = EOS % DV % Variables(iSen_T) % Values

    IF ( UseChemicalPotentialShift ) CALL ApplyChemicalPotentialShift_TABLE( Mp_T, Mn_T, OS_Mp, OS_Mn )

    ALLOCATE &
      ( Ps  (1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Ss  (1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Es  (1:EOS % DV % nPoints(1), &
             1:EOS % DV % nPoints(2), &
             1:EOS % DV % nPoints(3)) )
             
    ! Build full EOS without muons to determine the bounds of the EOS.
    ! Muons should not be so important to cause P, S, and E to go above the 
    ! bounds already calculated assuming Ym = 0.
    DO iD_T=1,EOS % DV % nPoints(1)
      DO iT_T=1,EOS % DV % nPoints(2)
        DO iYp_T=1,EOS % DV % nPoints(3)
        
          ! Now add electron component
          ! Initialize temperature, DensitY, Ye
          ElectronPhotonState % t   = T_T (iT_T)
          ElectronPhotonState % rho = D_T (iD_T)
          ElectronPhotonState % ye  = Yp_T(iYp_T)

          ! calculate electron quantities
          CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

          Eele = ElectronPhotonState % e
          Pele = ElectronPhotonState % p
          Sele = ElectronPhotonState % s

          Es(iD_T,iT_T,iYp_T) = 10.0d0**( E_T(iD_T,iT_T,iYp_T) ) + Eele - OS_E
          Ps(iD_T,iT_T,iYp_T) = ( P_T(iD_T,iT_T,iYp_T) ) + Pele - OS_P
          Ss(iD_T,iT_T,iYp_T) = 10.0d0**( S_T(iD_T,iT_T,iYp_T) ) + Sele - OS_S

        ENDDO
      ENDDO
    ENDDO
    
    CALL InitializeEOSComponentsInversion &
         ( D_T, T_T, Yp_T, Es, Ps, Ss, &
         HelmTable, MuonTable, Verbose_Option = Verbose )
             
    DEALLOCATE( Ps, Ss, Es )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( always, to: D_T, T_T, Yp_T, &
    !$OMP   UnitD, UnitT, UnitY, UnitP, UnitE, &
    !$OMP   UnitMl, UnitMp, UnitMn, UnitXp, UnitXn, UnitXa, UnitXh, &
    !$OMP   UnitAh, UnitGm, UnitEmp, UnitEmn, UnitSep, UnitSen, &
    !$OMP   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, &
    !$OMP   OS_Ah, OS_Gm, OS_Emp, OS_Emn, OS_Sep, OS_Sen, &
    !$OMP   P_T, S_T, E_T, Mp_T, Mn_T, Xp_T, Xn_T, Xa_T, &
    !$OMP   Xh_T, Ah_T, Gm_T, Emp_T, Emn_T, Sep_T, Sen_T, &
    !$OMP   HelmTable, MuonTable, &
    !$OMP   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y, Eos_MinD )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( D_T, T_T, Yp_T, &
    !$ACC   UnitD, UnitT, UnitY, UnitP, UnitE, &
    !$ACC   UnitMl, UnitMp, UnitMn, UnitXp, UnitXn, UnitXa, UnitXh, &
    !$ACC   UnitAh, UnitGm, UnitEmp, UnitEmn, UnitSep, UnitSen, &
    !$ACC   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, &
    !$ACC   OS_Ah, OS_Gm, OS_Emp, OS_Emn, OS_Sep, OS_Sen, &
    !$ACC   P_T, S_T, E_T, Mp_T, Mn_T, Xp_T, Xn_T, Xa_T, &
    !$ACC   Xh_T, Ah_T, Gm_T, Emp_T, Emn_T, Sep_T, Sen_T, &
    !$ACC   HelmTable, MuonTable, &
    !$ACC   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y, Eos_MinD )
#endif

#endif

  END SUBROUTINE InitializeEquationOfState_TABLE


  SUBROUTINE FinalizeEquationOfState_TABLE

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: D_T, T_T, Yp_T, &
    !$OMP   UnitD, UnitT, UnitY, UnitP, UnitE, &
    !$OMP   UnitMl, UnitMp, UnitMn, UnitXp, UnitXn, UnitXa, UnitXh, &
    !$OMP   UnitAh, UnitGm, UnitEmp, UnitEmn, UnitSep, UnitSen, &
    !$OMP   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, &
    !$OMP   OS_Ah, OS_Gm, OS_Emp, OS_Emn, OS_Sep, OS_Sen, &
    !$OMP   P_T, S_T, E_T, Mp_T, Mn_T, Xp_T, Xn_T, Xa_T, &
    !$OMP   Xh_T, Ah_T, Gm_T, Emp_T, Emn_T, Sep_T, Sen_T, &
    !$OMP   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y, Eos_MinD, &
    !$OMP   HelmTable, MuonTable )
#endif

    DEALLOCATE( D_T, T_T, Yp_T )

    DEALLOCATE( P_T  )
    DEALLOCATE( S_T  )
    DEALLOCATE( E_T  )
    DEALLOCATE( Mp_T )
    DEALLOCATE( Mn_T )
    DEALLOCATE( Xp_T )
    DEALLOCATE( Xn_T )
    DEALLOCATE( Xa_T )
    DEALLOCATE( Xh_T )
    DEALLOCATE( Ah_T )
    DEALLOCATE( Gm_T )
    DEALLOCATE( Emp_T )
    DEALLOCATE( Emn_T )
    DEALLOCATE( Sep_T )
    DEALLOCATE( Sen_T )

    IF ( .NOT. UsingExternalEOS ) THEN
      IF (ASSOCIATED(MuonTable))      DEALLOCATE( MuonTable )
      IF (ASSOCIATED(HelmTable)) DEALLOCATE( HelmTable )
      IF (ASSOCIATED(EOS))            DEALLOCATE( EOS )
    END IF
  
#endif

  END SUBROUTINE FinalizeEquationOfState_TABLE


  SUBROUTINE ApplyChemicalPotentialShift_TABLE( Mp, Mn, OS_Mp_loc, OS_Mn_loc )

    !For SFHo tables from
    !https://code.ornl.gov/astro/weaklib-tables/-/tree/master/SFHo/LowRes
    !up until git commit hash a36240ed
    !the neutron and proton chemical potentials
    !are tabulated subtracting the neutron-proton-mass difference
    !in order to use the same conventions as used in Chimera
    !to account for this, and get detailed balance factors
    !correct, we add the mass difference back in and then
    !add the SFHo reference masses for the chemical potential
    !(mn for Mun, mp for Mup)
    !For this renomalisation to the original SFHo tables,
    !we need to recalculate the offsets first

    REAL(dp), INTENT(inout), DIMENSION(:,:,:) :: Mp, Mn
    REAL(dp), INTENT(inout) :: OS_Mp_loc, OS_Mn_loc

    REAL(DP), PARAMETER :: neutron_mass = 939.56542052d0
    REAL(DP), PARAMETER :: proton_mass = 938.2720813d0
    REAL(DP), PARAMETER :: dmnp = 1.29333922d0
    REAL(DP) :: min_M, OS_M_new

    ! Apply the shift for proton chemical potential
    IF ( OS_Mp_loc > 0.0d0 ) THEN
      min_M = -0.5d0 * OS_Mp_loc
    ELSE
      min_M = MINVAL( 10.0d0**Mp )
    ENDIF
    OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + proton_mass + dmnp )
    Mp     = LOG10( 10.0d0**Mp - OS_Mp_loc + proton_mass + dmnp + OS_M_new)
    OS_Mp_loc    = OS_M_new

    ! Apply the shift for neutron chemical potential
    IF ( OS_Mn_loc > 0.0d0 ) THEN
      min_M = -0.5d0 * OS_Mn_loc
    ELSE
      min_M = MINVAL( 10.0d0**Mn )
    ENDIF
    OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + proton_mass + dmnp )
    Mn     = LOG10( 10.0d0**Mn - OS_Mn_loc + proton_mass + dmnp + OS_M_new)
    OS_Mn_loc    = OS_M_new

  END SUBROUTINE ApplyChemicalPotentialShift_TABLE


  SUBROUTINE ApplyEquationOfState_TABLE_Scalar &
    ( D, T, Ye, Ym, P, S, E, Mue, Mum, Mup, Mun, Xp, Xn, Xa, Xh, Gm )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: P, S, E, Mue, Mum, Mup, Mun, Xp, Xn, Xa, Xh, Gm

    REAL(DP) :: D_P, T_P, Ye_P, Ym_P, &
                Pbary, Sbary, Ebary, &
                Pele, Sele, Eele, &
                P_mu, S_mu, E_mu

    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = T / UnitY
    Ym_P = T / UnitY
    
#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE

    ! Calculate Electron Quantities
    ! Initialize Electron and Photon state
    ElectronPhotonState % t   = T_P
    ElectronPhotonState % rho = D_P
    ElectronPhotonState % ye  = Ye_P
    CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

    Eele = ElectronPhotonState % e
    Pele = ElectronPhotonState % p
    Sele = ElectronPhotonState % s
    Mue  = ElectronPhotonState % mue * UnitMl

    ! Calculate Muon Quantities
    MuonState % t     = T_P
    MuonState % rhoym = D_P * Ym_P
    CALL FullMuonEOS(MuonTable, MuonState)

    E_mu = MuonState % e
    P_mu = MuonState % p
    S_mu = MuonState % s
    Mum  = MuonState % mu * UnitMl
    
    ! --- Interpolate Pressure ----------------------------------------

    CALL ComputePressureBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, P )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, S, S_T, OS_S, One )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, E, E_T, OS_E, One )
           
    E = ( E + Eele + E_mu ) * UnitE
    P = ( P + Pele + P_mu ) * UnitP
    S = ( S + Sele + S_mu ) * UnitS

#else

  CALL ComputeTotalPES_TABLE_Scalar( D, T, Ye, Ym, P, E, S )

  ! Now take care of checmical potentials
  ElectronPhotonState % t   = T_P
  ElectronPhotonState % rho = D_P
  ElectronPhotonState % ye  = Ye_P
  CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)
  Mue  = ElectronPhotonState % mue * UnitMl

  ! Calculate Muon Quantities
  MuonState % t     = T_P
  MuonState % rhoym = D_P * Ym_P
  CALL FullMuonEOS(MuonTable, MuonState)
  Mum  = MuonState % mu * UnitMl

#endif

    ! --- Interpolate Proton Chemical Potential -----------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, Mup, Mp_T, OS_Mp, UnitMp )

    ! --- Interpolate Neutron Chemical Potential ----------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, Mun, Mn_T, OS_Mn, UnitMn )

    ! --- Interpolate Proton Mass Fraction ----------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, Xp, Xp_T, OS_Xp, UnitXp )

    ! --- Interpolate Neutron Mass Fraction ---------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, Xn, Xn_T, OS_Xn, UnitXn )

    ! --- Interpolate Alpha Mass Fraction -----------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, Xa, Xa_T, OS_Xa, UnitXa )

    ! --- Interpolate Heavy Mass Fraction -----------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, Xh, Xh_T, OS_Xh, UnitXh )

  END SUBROUTINE ApplyEquationOfState_TABLE_Scalar


  SUBROUTINE ApplyEquationOfState_TABLE_Vector &
    ( D, T, Ye, Ym, P, S, E, Mue, Mum, Mup, Mun, Xp, Xn, Xa, Xh, Gm )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: P(1:), S(1:), E(1:), Mue(1:), Mum(1:), Mup(1:), Mun(1:)
    REAL(DP), INTENT(out) :: Xp(1:), Xn(1:), Xa(1:), Xh(1:), Gm(1:)

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ApplyEquationOfState_TABLE_Scalar &
             ( D (iP), T (iP), Ye (iP), Ym (iP), P (iP), S (iP), E (iP), &
             Mue(iP), Mum(iP), Mup(iP), Mun(iP), Xp(iP), Xn(iP), Xa(iP), Xh(iP), Gm(iP) )

    END DO

  END SUBROUTINE ApplyEquationOfState_TABLE_Vector


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_G_E &
    ( D, E, Ye, Ym, T, Guess, Error )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)            :: D, E, Ye, Ym
    REAL(DP), INTENT(out)           :: T
    REAL(DP), INTENT(in)            :: Guess
    INTEGER,  INTENT(out)           :: Error

    REAL(DP) :: D_P, E_P, Ye_P, Ym_P, T_Lookup, T_Guess

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    E_P = E / UnitE
    Ye_P = Ye
    Ym_P = Ym

    T_Guess = Guess / Kelvin

    IF ( D_P >= Eos_MinD ) THEN
      CALL ComputeTemperatureWith_DEYpYl_Single_Guess_Error &
            ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, & 
                T_Lookup, T_Guess, Error, HelmTable, MuonTable )
    ELSE
      T_Lookup = T_Guess
      Error = 0
    END IF
    T = T_Lookup * Kelvin

#endif

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_G_E


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_G_NE &
    ( D, E, Ye, Ym, T, Guess )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)            :: D, E, Ye, Ym
    REAL(DP), INTENT(out)           :: T
    REAL(DP), INTENT(in)            :: Guess

    REAL(DP) :: D_P, E_P, Ye_P, Ym_P, T_Lookup, T_Guess

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    E_P = E / UnitE
    Ye_P = Ye
    Ym_P = Ym
    
    T_Guess = Guess / Kelvin

    IF ( D_P >= Eos_MinD ) THEN
      CALL ComputeTemperatureWith_DEYpYl_Single_Guess_NoError &
            ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
            T_Lookup, T_Guess, HelmTable, MuonTable )
    ELSE
      T_Lookup = T_Guess
    END IF

    T = T_Lookup * Kelvin

#endif

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_G_NE


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_NG_E &
    ( D, E, Ye, Ym, T, Error )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)            :: D, E, Ye, Ym
    REAL(DP), INTENT(out)           :: T
    INTEGER,  INTENT(out)           :: Error

    REAL(DP) :: D_P, E_P, Ye_P, Ym_P, T_Lookup

#ifdef MICROPHYSICS_WEAKLIB

    D_P  = D / UnitD
    E_P  = E / UnitE
    Ye_P = Ye
    Ym_P = Ym
    
    IF ( D_P >= Eos_MinD ) THEN
      CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error &
            ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
            T_Lookup, Error, HelmTable, MuonTable )
    ELSE
      T_Lookup = Zero
      Error = 0
    END IF

    T = T_Lookup * Kelvin

#endif

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_NG_E


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_NG_NE &
    ( D, E, Ye, Ym, T )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)            :: D, E, Ye, Ym
    REAL(DP), INTENT(out)           :: T

    REAL(DP) :: D_P, E_P, Ye_P, Ym_P, T_Lookup

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    E_P = E / UnitE
    Ye_P = Ye
    Ym_P = Ym
    
    IF ( D_P >= Eos_MinD ) THEN
      CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError &
            ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
            T_Lookup, HelmTable, MuonTable )
    ELSE
      T_Lookup = Zero
    END IF

    T = T_Lookup * Kelvin

#endif

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar_NG_NE


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector &
    ( D, E, Ye, Ym, T, Guess_Option, Error_Option )

    REAL(DP), INTENT(in )           :: D(1:), E(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)           :: T(1:)
    REAL(DP), INTENT(in ), OPTIONAL :: Guess_Option(1:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(1:)

    INTEGER  :: iP, nP
    INTEGER  :: Error(SIZE(D))
    REAL(DP) :: D_P, E_P, Ye_P, Ym_P, T_Lookup, T_Guess

    Error = 0

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

    IF( PRESENT( Guess_Option ) )THEN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup, T_Guess ) &
    !$OMP MAP( to: D, E, Ye, Ym, D_T, T_T, Yp_T, E_T, Guess_Option ) &
    !$OMP MAP( from: T, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup, T_Guess ) &
    !$ACC COPYIN( D, E, Ye, Ym, D_T, T_T, Yp_T, E_T, Guess_Option ) &
    !$ACC COPYOUT( T, Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup, T_Guess )
#endif
      DO iP = 1, nP

        D_P = D(iP) / UnitD
        E_P = E(iP) / UnitE
        Ye_P = Ye(iP)
        Ym_P = Ym(iP)
        
        T_Guess = Guess_Option(iP) / Kelvin

        IF ( D_P >= Eos_MinD ) THEN
          CALL ComputeTemperatureWith_DEYpYl_Single_Guess_Error &
                ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
                T_Lookup, T_Guess, Error(iP), HelmTable, MuonTable )
        ELSE
          T_Lookup = T_Guess
          Error(iP) = 0
        END IF

        T(iP) = T_Lookup * Kelvin

      END DO

    ELSE

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup ) &
    !$OMP MAP( to: D, E, Ye, Ym, D_T, T_T, Yp_T, E_T ) &
    !$OMP MAP( from: T, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup ) &
    !$ACC COPYIN( D, E, Ye, Ym, D_T, T_T, Yp_T, E_T ) &
    !$ACC COPYOUT( T, Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup )
#endif
      DO iP = 1, nP

        D_P = D(iP) / UnitD
        E_P = E(iP) / UnitE
        Ye_P = Ye(iP)
        Ym_P = Ym(iP)
        
        IF ( D_P >= Eos_MinD ) THEN
          CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error &
                ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, T_Lookup, &
                  Error(iP), HelmTable, MuonTable )
        ELSE
          T_Lookup = Zero
          Error(iP) = 0
        END IF

        T(iP) = T_Lookup * Kelvin

      END DO

    END IF

    IF ( ANY( Error > 0 ) ) THEN
      DO iP = 1, nP
        IF ( Error(iP) > 0 ) THEN
          CALL DescribeEOSComponentsInversionError( Error(iP) )
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM &
          !$OMP ( D(iP), E(iP), Ye(iP), Ym(iP) )
#elif defined(THORNADO_OACC)
          !$ACC UPDATE HOST &
          !$ACC ( D(iP), E(iP), Ye(iP), Ym(iP) )
#endif
          D_P = D(iP) / UnitD
          E_P = E(iP) / UnitE
          Ye_P = Ye(iP)
          Ym_P = Ym(iP)
          
          WRITE(*,*) '[ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector] Error'
          WRITE(*,'(a,i5,4es23.15)') '  iP, D, E, Ye, Ym : ', iP, D_P, E_P, Ye_P, Ym_P
        END IF
      END DO
      IF( .NOT. PRESENT( Error_Option ) ) STOP
    END IF

#endif

    IF( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector

  SUBROUTINE ComputePressureFromPrimitive_TABLE_Scalar( D, Ev, Ne, Nm, P )

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
  !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne, Nm
    REAL(DP), INTENT(out) :: P

    REAL(DP) :: Em, T, Ye, Ym

    Em = Ev / D              ! --- Internal Energy per Mass
    Ye = Ne / D * BaryonMass ! --- Electron Fraction
    Ym = Nm / D * BaryonMass ! --- Muon Fraction
  
    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
           ( D, Em, Ye, Ym, T )

    CALL ComputePressure_TABLE_Scalar &
           ( D, T, Ye, Ym, P )

  END SUBROUTINE ComputePressureFromPrimitive_TABLE_Scalar

  SUBROUTINE ComputePressureFromPrimitive_TABLE_Vector( D, Ev, Ne, Nm, P )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: D, Ev, Ne, Nm
    REAL(DP), DIMENSION(1:), INTENT(out) :: P

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputePressureFromPrimitive_TABLE_Scalar &
             ( D(iP), Ev(iP), Ne(iP), Nm(iP), P(iP) )

    END DO

  END SUBROUTINE ComputePressureFromPrimitive_TABLE_Vector


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Scalar &
    ( D, Em, Ye, Ym, P )

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
  !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Em, Ye, Ym
    REAL(DP), INTENT(out) :: P

    REAL(DP) :: D_P, E_P, Ye_P, Ym_P, T_P, T

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE
  TYPE(ElectronPhotonStateType) :: ElectronPhotonState
  TYPE(MuonStateType) :: MuonState
  REAL(DP) :: Pele, P_mu
#endif

#ifdef MICROPHYSICS_WEAKLIB

    D_P  = D  / UnitD
    E_P  = Em / UnitE
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    
    IF ( D_P >= Eos_MinD ) THEN
      CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError &
            ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
            T_P, HelmTable, MuonTable )
    ELSE
      P = Zero
    END IF

    T = T_P * UnitT

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE

    CALL ComputePressureBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, P )

    ! Calculate Electron Quantities
    ElectronPhotonState % t   = T_P
    ElectronPhotonState % rho = D_P
    ElectronPhotonState % ye  = Ye_P
    
    CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)
    Pele = ElectronPhotonState % p

    ! Calculate Muon Quantities
    MuonState % t     = T_P
    MuonState % rhoym = D_P * Ym_P

    CALL FullMuonEOS(MuonTable, MuonState)
    P_mu = MuonState % p
           
    P = ( P + Pele + P_mu ) * UnitP
#else
  
  CALL ComputeDependentVariableTotal_TABLE_Scalar &
    ( D, T, Ye, Ym, P, P_T, OS_P, &
    UnitP, 1, 0, 0 )

#endif
#endif

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Vector &
    ( D, Em, Ye, Ym, P )

    REAL(DP), INTENT(in)  :: D(:), Em(:), Ye(:), Ym(:)
    REAL(DP), INTENT(out) :: P(:)

    INTEGER iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputePressureFromSpecificInternalEnergy_TABLE_Scalar &
             ( D(iP), Em(iP), Ye(iP), Ym(iP), P(iP) )

    END DO

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Vector

  ! DECIDE HOW TO HANDLE SOUND SPEED BEFORE DOING THIS. SOUND SPEED FROM
  ! EACH Combined COMPONENT MIGHT NOT BE DOABLE, SO FIND A SOLUTION.
  ! CAREFUL BECAUSE YOU NEED DERIVATIVES AND WHEN YOU VARY BOTH YE AND Ym 
  ! IT MIGHT NOT BE STRAIGHTFORWARD. CURRENTLY BOLLIG'S FORMULA IS USED
  SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Scalar( D, Ev, Ne, Nm, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne, Nm
    REAL(DP), INTENT(out) :: Cs

    REAL(DP) :: P, T, Ye, Ym, Em, Gm

    Em = Ev / D
    Ye = Ne * ( BaryonMass / D )
    Ym = Nm * ( BaryonMass / D )
    
    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
           ( D, Em, Ye, Ym, T )

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE
    CALL CalculateSoundSpeed( D / UnitD, T / UnitT, Ye / UnitY, Ym / UnitY, D_T, T_T, Yp_T, &
        P_T, OS_P, E_T, OS_E, S_T, OS_S, HelmTable, MuonTable, Gm, Cs, .TRUE.)
#else
    CALL CalculateSoundSpeed( D / UnitD, T / UnitT, Ye / UnitY, Ym / UnitY, D_T, T_T, Yp_T, &
        P_T, OS_P, E_T, OS_E, S_T, OS_S, HelmTable, MuonTable, Gm, Cs, .FALSE.)
#endif

    Gm = Gm * UnitGm
    Cs = Cs * Centimeter / Second

    ! Temporary sound speed calculation to test accuracy of solver
    ! CALL ComputeDependentVariableBaryons_TABLE_Scalar &
    !     ( D, T, Ye, Ym, Gm, Gm_T, OS_Gm, Units_V = 1.0_DP )

    ! CALL ComputePressure_TABLE_Scalar( D, T, Ye, Ym, P)
    ! Cs = SQRT(Gm * P / D)

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Scalar

  SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Vector( D, Ev, Ne, Nm, Cs )

    REAL(DP), INTENT(in)  :: D(1:), Ev(1:), Ne(1:), Nm(1:)
    REAL(DP), INTENT(out) :: Cs(1:)

    REAL(DP), DIMENSION(SIZE(D)) :: P, T, Ye, Ym, Em, Gm
    
    INTEGER :: iP, nP
    
    nP = SIZE(D)
    
    Em = Ev / D
    Ye = Ne * ( BaryonMass / D )
    Ym = Nm * ( BaryonMass / D )
      
    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector &
           ( D, Em, Ye, Ym, T )

    DO iP=1,nP

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE
      CALL CalculateSoundSpeed( D(iP) / UnitD, T(iP) / UnitT, &
          Ye(iP) / UnitY, Ym(iP) / UnitY, &
          D_T, T_T, Yp_T, P_T, OS_P, E_T, OS_E, S_T, OS_S, &
          HelmTable, MuonTable, Gm(iP), Cs(iP), .TRUE.)
#else
      CALL CalculateSoundSpeed( D(iP) / UnitD, T(iP) / UnitT, &
          Ye(iP) / UnitY, Ym(iP) / UnitY, &
          D_T, T_T, Yp_T, P_T, OS_P, E_T, OS_E, S_T, OS_S, &
          HelmTable, MuonTable, Gm(iP), Cs(iP), .FALSE.)
#endif
      Gm(iP) = Gm(iP) * UnitGm
      Cs(iP) = Cs(iP) * Centimeter / Second

    ENDDO

    ! Temporary sound speed calculation to test accuracy of solver
    ! CALL ComputeDependentVariableBaryons_TABLE_Vector &
    !     ( D, T, Ye, Ym, Gm, Gm_T, OS_Gm, Units_V = 1.0_DP )

    ! CALL ComputePressure_TABLE_Vector( D, T, Ye, Ym, P)
    ! Cs = SQRT( Gm * P / D )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Vector


  SUBROUTINE ComputeTemperatureFromPressure_TABLE_Scalar( D, P, Ye, Ym, T )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, P, Ye, Ym
    REAL(DP), INTENT(out) :: T

    INTEGER  :: Error
    REAL(DP) :: D_P, P_P, Ye_P, Ym_P, T_P

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    P_P = P / UnitP
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY

    IF ( D_P >= Eos_MinD ) THEN
      CALL ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error &
            ( D_P, P_P, Ye_P, Ym_P, D_T, T_T, Yp_T, P_T, OS_P, T_P, &
              Error, HelmTable, MuonTable )
    ELSE
      T_P = Zero
      Error = 0
    END IF

    T = T_P * UnitT

#endif

  END SUBROUTINE ComputeTemperatureFromPressure_TABLE_Scalar


  SUBROUTINE ComputeTemperatureFromPressure_TABLE_Vector( D, P, Ye, Ym, T )

    REAL(DP), INTENT(in)  :: D(1:), P(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: T(1:)

    INTEGER :: iP, nP

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputeTemperatureFromPressure_TABLE_Scalar &
             ( D(iP), P(iP), Ye(iP), Ym(iP), T(iP) )

    END DO

#endif

  END SUBROUTINE ComputeTemperatureFromPressure_TABLE_Vector

  
  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Scalar &
    ( D, T, Ye, Ym, Ev, Em, Ne, Nm )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: Ev, Em, Ne, Nm
    
#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
    REAL(DP) :: Eele, E_mu
#endif
    ! --- Interpolate Specific Internal Energy ------------------------

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, Em, E_T, OS_E, One )
           
    ! Calculate Electron Quantities
    ! Initialize Electron and Photon state
    ElectronPhotonState % t   = T  / UnitT
    ElectronPhotonState % rho = D  / UnitD
    ElectronPhotonState % ye  = Ye / UnitY
    CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

    Eele = ElectronPhotonState % e

    ! Calculate Muon Quantities
    MuonState % t = T / UnitT
    MuonState % rhoym = D * Ym / UnitD / UnitY
    CALL FullMuonEOS(MuonTable, MuonState)

    E_mu = MuonState % e
           
    Em = ( Em + Eele + E_mu ) * UnitE
#else
  
  CALL ComputeDependentVariableTotal_TABLE_Scalar &
          ( D, T, Ye, Ym, Em, E_T, OS_E, &
          UnitE, 0, 1, 0 )

#endif

    Ev = Em * D              ! --- Internal Energy per Unit Volume
    Ne = Ye * D / BaryonMass ! --- Electrons per Unit Volume
    Nm = Ym * D / BaryonMass ! --- Muons per Unit Volume

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Scalar

  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Vector &
    ( D, T, Ye, Ym, Ev, Em, Ne, Nm )

    REAL(DP), INTENT(in)  :: D (1:), T (1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: Ev(1:), Em(1:), Ne(1:), Nm(1:)

    INTEGER :: iP, nP

    nP = SIZE( D )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: D, T, Ye, Ym ) &
    !$OMP MAP( from: Em, Ev, Ne )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( D, T, Ye, Ym ) &
    !$ACC COPYOUT( Em, Ev, Ne )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeThermodynamicStates_Primitive_TABLE_Scalar &
             ( D(iP), T(iP), Ye(iP), Ym(iP), Ev(iP), Em(iP), Ne(iP), Nm(iP) )

    END DO

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Vector

  SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Scalar &
    ( D, Ev, Ne, Nm, T, Em, Ye, Ym )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne, Nm
    REAL(DP), INTENT(out) :: T, Em, Ye, Ym

    Em  = Ev  / D              ! --- Internal Energy per Mass
    Ye  = Ne  / D * BaryonMass ! --- Electron Fraction
    Ym = Nm / D * BaryonMass ! --- Muon Fraction

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
           ( D, Em, Ye, Ym, T )

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Scalar

  SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Vector &
    ( D, Ev, Ne, Nm, T, Em, Ye, Ym )

    REAL(DP), DIMENSION(1:), INTENT(in)  :: D, Ev, Ne, Nm
    REAL(DP), DIMENSION(1:), INTENT(out) :: T, Em, Ye, Ym

    INTEGER :: iP, nP, Error(SIZE(D))
    REAL(DP) :: D_P, E_P, Ye_P, Ym_P

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( to: D, Ev, Ne, Nm ) &
    !$OMP MAP( from: T, Em, Ye, Ym, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( D, Ev, Ne, Nm ) &
    !$ACC COPYOUT( T, Em, Ye, Ym, Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      Em(iP) = Ev(iP) / D(iP)              ! --- Internal Energy per Mass
      Ye(iP) = Ne(iP) / D(iP) * BaryonMass ! --- Electron Fraction
      Ym(iP) = Nm(iP) / D(iP) * BaryonMass ! --- Muon Fraction

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
             ( D(iP), Em(iP), Ye(iP), Ym(iP), T(iP), Error(iP) )

    END DO

    IF ( ANY( Error > 0 ) ) THEN
      DO iP = 1, nP
        IF ( Error(iP) > 0 ) THEN
          CALL DescribeEOSComponentsInversionError( Error(iP) )
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM &
          !$OMP ( D(iP), Em(iP), Ye(iP), Ym(iP) )
#elif defined(THORNADO_OACC)
          !$ACC UPDATE HOST &
          !$ACC ( D(iP), Em(iP), Ye(iP), Ym(iP) )
#endif
          D_P  = D(iP)  / ( Gram / Centimeter**3 )
          E_P  = Em(iP) / ( Erg / Gram )
          Ye_P = Ye(iP)
          Ym_P = Ym(iP)
          WRITE(*,*)                 '[ComputeThermodynamicStates_Auxiliary_TABLE_Vector] Error'
          WRITE(*,'(a,i5,5es23.15)') '  iP, D, E, Ye, Ym : ', iP, D_P, E_P, Ye_P, Ym_P
        END IF
      END DO
      STOP
    END IF

#endif

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Vector

  SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Scalar &
    ( D, Ev, Ne, Nm, P, T, Ye, Ym, S, Em, Gm, Cs )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, Ev, Ne, Nm
    REAL(DP), INTENT(out) :: P, T, Ye, Ym, S, Em, Gm, Cs

    Em  = Ev / D
    Ye  = Ne * ( BaryonMass / D )
    Ym = Nm * ( BaryonMass / D )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
           ( D, Em, Ye, Ym, T )

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE
    CALL CalculateSoundSpeed( D / UnitD, T / UnitT, Ye / UnitY, Ym / UnitY, &
        D_T, T_T, Yp_T, P_T, OS_P, E_T, OS_E, S_T, OS_S, &
        HelmTable, MuonTable, Gm, Cs, .TRUE.)
#else
    CALL CalculateSoundSpeed( D / UnitD, T / UnitT, Ye / UnitY, Ym / UnitY, &
        D_T, T_T, Yp_T, P_T, OS_P, E_T, OS_E, S_T, OS_S, &
        HelmTable, MuonTable, Gm, Cs, .FAlSE.)
#endif

    Gm = Gm * UnitGm
    Cs = Cs * Centimeter / Second

    ! Temporary sound speed calculation to test accuracy of solver
    ! Cs = SQRT(Gm * P / D)

  END SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Scalar

  SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Vector &
    ( D, Ev, Ne, Nm, P, T, Ye, Ym, S, Em, Gm, Cs )

    REAL(DP), INTENT(in)  :: D(1:), Ev(1:), Ne(1:), Nm(1:)
    REAL(DP), INTENT(out) :: P(1:), T (1:), Ye(1:), Ym(1:), &
                             S(1:), Em(1:), Gm(1:), Cs(1:)

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputeAuxiliary_Fluid_TABLE_Scalar &
             ( D(iP), Ev(iP), Ne(iP), Nm(iP), &
               P(iP), T (iP), Ye(iP), Ym(iP), &
               S(iP), Em(iP), Gm(iP), Cs(iP) )

    END DO

  END SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Vector

  SUBROUTINE ComputePressure_TABLE_Scalar &
    ( D, T, Ye, Ym, P, dPdD_Option, dPdT_Option, dPdYe_Option, dPdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: P
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dPdD_Local, dPdT_Local, dPdYe_Local, dPdYm_Local
    REAL(DP), POINTER :: dPdD      , dPdT      , dPdYe      , dPdYm

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
    REAL(DP) :: Pele, P_mu
#endif

    ComputeDerivatives &
      =      PRESENT( dPdD_Option ) &
        .OR. PRESENT( dPdT_Option ) &
        .OR. PRESENT( dPdYe_Option ) &
        .OR. PRESENT( dPdYm_Option )

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

      IF( PRESENT( dPdYe_Option ) )THEN
        dPdYe => dPdYe_Option
      ELSE
        dPdYe => dPdYe_Local
      END IF

      IF( PRESENT( dPdYm_Option ) )THEN
        dPdYm => dPdYm_Option
      ELSE
        dPdYm => dPdYm_Local
      END IF
      
#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE

      WRITE(*,*) 'Still needs to be implemented, can be a bit tedious, maybe I can avoid doing this'
      STOP

#else
  
  CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
        ( D, T, Ye, Ym, P, dPdD, dPdT, dPdYe, dPdYm, P_T, OS_P, &
        UnitP, 1, 0, 0 )

#endif

    ELSE

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE

      CALL ComputePressureBaryons_TABLE_Scalar &
          ( D, T, Ye, Ym, P )

      ! Calculate Electron Quantities
      ! Initialize Electron and Photon state
      ElectronPhotonState % t   = T  / UnitT
      ElectronPhotonState % rho = D  / UnitD
      ElectronPhotonState % ye  = Ye / UnitY
      CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

      Pele = ElectronPhotonState % p

      ! Calculate Muon Quantities
      MuonState % t = T / UnitT
      MuonState % rhoym = D * Ym / UnitD / UnitY
      CALL FullMuonEOS(MuonTable, MuonState)

      P_mu = MuonState % p
      
      P = ( P + Pele + P_mu ) * UnitP
#else

  CALL ComputeDependentVariableTotal_TABLE_Scalar &
      ( D, T, Ye, Ym, P, P_T, OS_P, &
      UnitP, 1, 0, 0 )

#endif

    END IF

  END SUBROUTINE ComputePressure_TABLE_Scalar

  SUBROUTINE ComputePressure_TABLE_Vector &
    ( D, T, Ye, Ym, P, dPdD_Option, dPdT_Option, dPdYe_Option, dPdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: P(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dPdD_Local, dPdT_Local, dPdYe_Local, dPdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dPdD      , dPdT      , dPdYe, dPdYm

    ComputeDerivatives &
      =      PRESENT( dPdD_Option ) &
        .OR. PRESENT( dPdT_Option ) &
        .OR. PRESENT( dPdYe_Option ) &
        .OR. PRESENT( dPdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dPdD_Local(nP), dPdT_Local(nP), dPdYe_Local(nP), dPdYm_Local(nP) )

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

      IF( PRESENT( dPdYe_Option ) )THEN
        dPdYe(1:nP) => dPdYe_Option(:)
      ELSE
        dPdYe(1:nP) => dPdYe_Local(:)
      END IF

      IF( PRESENT( dPdYm_Option ) )THEN
        dPdYm(1:nP) => dPdYm_Option(:)
      ELSE
        dPdYm(1:nP) => dPdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Vector &
             ( D, T, Ye, Ym, P, dPdD, dPdT, dPdYe, dPdYm, P_T, OS_P, &
             UnitP, 1, 0, 0 )

    ELSE

      CALL ComputeDependentVariableTotal_TABLE_Vector &
             ( D, T, Ye, Ym, P, P_T, OS_P, &
             UnitP, 1, 0, 0 )

    END IF

  END SUBROUTINE ComputePressure_TABLE_Vector


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Scalar &
    ( D, T, Ye, Ym, E, dEdD_Option, dEdT_Option, dEdYe_Option, dEdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: E
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dEdD_Local, dEdT_Local, dEdYe_Local, dEdYm_Local
    REAL(DP), POINTER :: dEdD,       dEdT,       dEdYe, dEdYm

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
    REAL(DP) :: Eele, E_mu
#endif

    ComputeDerivatives &
      =      PRESENT( dEdD_Option ) &
        .OR. PRESENT( dEdT_Option ) &
        .OR. PRESENT( dEdYe_Option ) &
        .OR. PRESENT( dEdYm_Option )

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

      IF( PRESENT( dEdYe_Option ) )THEN
        dEdYe => dEdYe_Option
      ELSE
        dEdYe => dEdYe_Local
      END IF

      IF( PRESENT( dEdYm_Option ) )THEN
        dEdYm => dEdYm_Option
      ELSE
        dEdYm => dEdYm_Local
      END IF

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE

      WRITE(*,*) 'Still needs to be implemented, can be a bit tedious, maybe I can avoid doing this'
      STOP
#else
  
  CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
        ( D, T, Ye, Ym, E, dEdD, dEdT, dEdYe, dEdYm, E_T, OS_E, &
        UnitE, 0, 1, 0 )

#endif

    ELSE

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, E, E_T, OS_E, One )
           
    ! Calculate Electron Quantities
    ! Initialize Electron and Photon state
    ElectronPhotonState % t   = T  / UnitT
    ElectronPhotonState % rho = D  / UnitD
    ElectronPhotonState % ye  = Ye / UnitY

    CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)
    Eele = ElectronPhotonState % e

    ! Calculate Muon Quantities
    MuonState % t = T / UnitT
    MuonState % rhoym = D * Ym / UnitD / UnitY
    
    CALL FullMuonEOS(MuonTable, MuonState)
    E_mu = MuonState % e
           
    E = ( E + Eele + E_mu ) * UnitE

#else
  
  CALL ComputeDependentVariableTotal_TABLE_Scalar &
      ( D, T, Ye, Ym, E, E_T, OS_E, &
      UnitE, 0, 1, 0 )

#endif
    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Vector &
    ( D, T, Ye, Ym, E, dEdD_Option, dEdT_Option, dEdYe_Option, dEdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: E(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dEdD_Local, dEdT_Local, dEdYe_Local, dEdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dEdD      , dEdT      , dEdYe, dEdYm

    ComputeDerivatives &
      =      PRESENT( dEdD_Option ) &
        .OR. PRESENT( dEdT_Option ) &
        .OR. PRESENT( dEdYe_Option ) &
        .OR. PRESENT( dEdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dEdD_Local(nP), dEdT_Local(nP), dEdYe_Local(nP), dEdYm_Local(nP) )

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

      IF( PRESENT( dEdYe_Option ) )THEN
        dEdYe(1:nP) => dEdYe_Option(:)
      ELSE
        dEdYe(1:nP) => dEdYe_Local(:)
      END IF

      IF( PRESENT( dEdYm_Option ) )THEN
        dEdYm(1:nP) => dEdYm_Option(:)
      ELSE
        dEdYm(1:nP) => dEdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Vector &
             ( D, T, Ye, Ym, E, dEdD, dEdT, dEdYe, dEdYm, E_T, OS_E, &
             UnitE, 0, 1, 0 )

    ELSE

#ifdef INTERPOLATION_SPLIT_TABLE_SEPARATE

      WRITE(*,*) 'Not yet implemented in ComputeSpecificInternalEnergy_TABLE_Vector'
      STOP
      
#else

      CALL ComputeDependentVariableTotal_TABLE_Vector &
          ( D, T, Ye, Ym, E, E_T, OS_E, &
          UnitE, 0, 1, 0 )

#endif

    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Vector

  ! The Mue and Mum subroutines are special since you only need electron (muon) info
  SUBROUTINE ComputeElectronChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdYe_Option, dMdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym ! Ym is a dummy variable
    REAL(DP), INTENT(out) :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYm_Option
    
    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdYe_Local, dMdYm_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdYe, dMdYm
    
    REAL(DP) :: dD, dT
	  REAL(DP) :: aD, aT
    
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    
    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdYe_Option ) &
        .OR. PRESENT( dMdYm_Option )

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

      IF( PRESENT( dMdYe_Option ) )THEN
        dMdYe => dMdYe_Option
      ELSE
        dMdYe => dMdYe_Local
      END IF

      IF( PRESENT( dMdYm_Option ) )THEN
        dMdYm => dMdYm_Option
      ELSE
        dMdYm => dMdYm_Local
      END IF
    
      ! ADD STUFF
      dMdD = 0.0_dp
      dMdT = 0.0_dp
      dMdYe = 0.0_dp
      dMdYm = 0.0_dp

      ! Calculate Electron Quantities
      ! Initialize Electron and Photon state
      ElectronPhotonState % t    = T  / UnitT
      ElectronPhotonState % rho  = D  / UnitD
      ElectronPhotonState % ye   = Ye / UnitY

      CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

      M = ElectronPhotonState % mue * UnitMl
      
    ELSE
      
      ! Calculate Electron Quantities
      ! Initialize Electron and Photon state
      ElectronPhotonState % t    = T  / UnitT
      ElectronPhotonState % rho  = D  / UnitD
      ElectronPhotonState % ye   = Ye / UnitY
      
      CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

      M = ElectronPhotonState % mue * UnitMl
      
    ENDIF
      
  END SUBROUTINE ComputeElectronChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeElectronChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdYe_Option, dMdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:) ! Ym is a dummy variable
    REAL(DP), INTENT(out)                   :: M(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYm_Option(1:)
    
    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dMdD_Local, dMdT_Local, dMdYe_Local, dMdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdYe      , dMdYm

    INTEGER :: iP, nP

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdYe_Option ) &
        .OR. PRESENT( dMdYm_Option )

    nP = SIZE( D )

    IF( ComputeDerivatives )THEN

      ALLOCATE( dMdD_Local(nP), dMdT_Local(nP), dMdYe_Local(nP), dMdYm_Local(nP) )

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

      IF( PRESENT( dMdYe_Option ) )THEN
        dMdYe(1:nP) => dMdYe_Option(:)
      ELSE
        dMdYe(1:nP) => dMdYe_Local(:)
      END IF

      IF( PRESENT( dMdYm_Option ) )THEN
        dMdYm(1:nP) => dMdYm_Option(:)
      ELSE
        dMdYm(1:nP) => dMdYm_Local(:)
      END IF

      DO iP=1,nP
          CALL ComputeElectronChemicalPotential_TABLE_Scalar &
                (D(iP), T(iP), Ye(iP), Ym(iP), M(iP), dMdD(iP), dMdT(iP), dMdYe(iP), dMdYm(iP))
      END DO

    ELSE 
    
      DO iP=1,nP
        CALL ComputeElectronChemicalPotential_TABLE_Scalar &
                (D(iP), T(iP), Ye(iP), Ym(iP), M(iP))
      END DO
      
    END IF

  END SUBROUTINE ComputeElectronChemicalPotential_TABLE_Vector

  SUBROUTINE ComputeMuonChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdYe_Option, dMdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym ! Ye is a dummy variable
    REAL(DP), INTENT(out) :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYm_Option
    
    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdYe_Local, dMdYm_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdYe, dMdYm
    
    REAL(DP) :: aD, aT, dD, dT
    REAL(DP) :: D_P, T_P, Ym_P
    INTEGER  :: i, iD, iT

    TYPE(MuonStateType) :: MuonState
    
    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdYe_Option ) &
        .OR. PRESENT( dMdYm_Option )

    D_P  = D  / UnitD
    Ym_P = Ym / UnitY
    T_P  = T  / UnitT

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

      IF( PRESENT( dMdYe_Option ) )THEN
        dMdYe => dMdYe_Option
      ELSE
        dMdYe => dMdYe_Local
      END IF

      IF( PRESENT( dMdYm_Option ) )THEN
        dMdYm => dMdYm_Option
      ELSE
        dMdYm => dMdYm_Local
      END IF
    
      ! ADD STUFF
      dMdD = 0.0_dp
      dMdT = 0.0_dp
      dMdYe = 0.0_dp
      dMdYm = 0.0_dp
      
      CALL GetIndexAndDelta_Log( T_P,   MuonTable % t(:), iT, dT )

      DO i = 1, MuonTable % nPointsDen
        IF (MuonTable % rhoym(iT,i) >= D_P*Ym_P) THEN
          iD = i
          dD = LOG10( D_P*Ym_P / MuonTable % rhoym(iT,i) ) / LOG10( MuonTable % rhoym(iT,i+1) / MuonTable % rhoym(iT,i) )
          EXIT
        ENDIF 
      END DO

      aD = 1.0_dp / ( D_P * LOG10( MuonTable % rhoym(iT,iD+1) / MuonTable % rhoym(iT,iD) ) )
      aT = 1.0_dp / ( T_P * LOG10( MuonTable % t(iT+1) / MuonTable % t(iT) ) )

      CALL LinearInterpDeriv_Array_Point &
             ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % mu), M, &
               dMdD, dMdT )
      
      dMdYm = dMdD * D_P  * UnitMl / UnitD ! make sure the derivative is wr2 ym, not rhoym
      dMdD  = dMdD * Ym_P * UnitMl / UnitY ! make sure the derivative is wr2 rho, not rhoym
      
      ! THE WAY THE ABOVE IS DONE I THINK IS WRONG, the dD and dT are technically LOG
      ! BUT YOU ARE GIVING THE NON LOGGED VALUES. BOTTOM LINE RECHECK EVERYTHING
      dMdYm = 0.0_dp
      dMdD  = 0.0_dp
      dMdYe = 0.0_dp
      
    ELSE

      MuonState % t = T_P
      MuonState % rhoym = D_P * Ym_P
      
      CALL FullMuonEOS(MuonTable, MuonState)

      M = MuonState % mu

    ENDIF

    M = M * UnitMl
    
  END SUBROUTINE ComputeMuonChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeMuonChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdYe_Option, dMdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:) ! Ym is a dummy variable
    REAL(DP), INTENT(out)                   :: M(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYm_Option(1:)
    
    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dMdD_Local, dMdT_Local, dMdYe_Local, dMdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdYe      , dMdYm
    INTEGER :: iP, nP
    
    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdYe_Option ) &
        .OR. PRESENT( dMdYm_Option )

    nP = SIZE( D )

    IF( ComputeDerivatives )THEN

      ALLOCATE( dMdD_Local(nP), dMdT_Local(nP), dMdYe_Local(nP), dMdYm_Local(nP) )

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

      IF( PRESENT( dMdYe_Option ) )THEN
        dMdYe(1:nP) => dMdYe_Option(:)
      ELSE
        dMdYe(1:nP) => dMdYe_Local(:)
      END IF

      IF( PRESENT( dMdYm_Option ) )THEN
        dMdYm(1:nP) => dMdYm_Option(:)
      ELSE
        dMdYm(1:nP) => dMdYm_Local(:)
      END IF

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP MAP( to: D, T, Ye, Ym ) &
      !$OMP MAP( from: M, dMdD, dMdT, dMdYe, dMdYm )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Ye, Ym ) &
      !$ACC COPYOUT( M, dMdD, dMdT, dMdYe, dMdYm )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iP=1,nP
          CALL ComputeMuonChemicalPotential_TABLE_Scalar &
                (D(iP), T(iP), Ye(iP), Ym(iP), M(iP), dMdD(iP), dMdT(iP), dMdYe(iP), dMdYm(iP))
      END DO

    ELSE 
    
#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP MAP( to: D, T, Ye, Ym ) &
      !$OMP MAP( from: M )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Ye, Ym ) &
      !$ACC COPYOUT( M )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
      DO iP=1,nP
        CALL ComputeMuonChemicalPotential_TABLE_Scalar &
                (D(iP), T(iP), Ye(iP), Ym(iP), M(iP))
      END DO
      
    END IF

  END SUBROUTINE ComputeMuonChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeProtonChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdYe_Option, dMdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdYe_Local, dMdYm_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdYe, dMdYm
    
    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdYe_Option ) &
        .OR. PRESENT( dMdYm_Option )

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

      IF( PRESENT( dMdYe_Option ) )THEN
        dMdYe => dMdYe_Option
      ELSE
        dMdYe => dMdYe_Local
      END IF

      IF( PRESENT( dMdYm_Option ) )THEN
        dMdYm => dMdYm_Option
      ELSE
        dMdYm => dMdYm_Local
      END IF
      
      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, M, dMdD, dMdT, dMdYe, dMdYm, Mp_T, OS_Mp, UnitMp )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, M, Mp_T, OS_Mp, UnitMp )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE_Scalar


  ! THIS SUBROUTINE IS UNTOUCHED, NOT SURE HOW TO HANDLE IT
  SUBROUTINE ComputeProtonChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdYe_Option, dMdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: M(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dMdD_Local, dMdT_Local, dMdYe_Local, dMdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdYe      , dMdYm

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdYe_Option ) &
        .OR. PRESENT( dMdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dMdD_Local(nP), dMdT_Local(nP), dMdYe_Local(nP), dMdYm_Local(nP) )

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

      IF( PRESENT( dMdYe_Option ) )THEN
        dMdYe(1:nP) => dMdYe_Option(:)
      ELSE
        dMdYe(1:nP) => dMdYe_Local(:)
      END IF

      IF( PRESENT( dMdYm_Option ) )THEN
        dMdYm(1:nP) => dMdYm_Option(:)
      ELSE
        dMdYm(1:nP) => dMdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, M, dMdD, dMdT, dMdYe, dMdYm, Mp_T, OS_Mp, UnitMp )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, M, Mp_T, OS_Mp, UnitMp )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdYe_Option, dMdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdYe_Local, dMdYm_Local
    REAL(DP), POINTER :: dMdD      , dMdT      , dMdYe, dMdYm

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdYe_Option ) &
        .OR. PRESENT( dMdYm_Option )

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

      IF( PRESENT( dMdYe_Option ) )THEN
        dMdYe => dMdYe_Option
      ELSE
        dMdYe => dMdYe_Local
      END IF

      IF( PRESENT( dMdYm_Option ) )THEN
        dMdYm => dMdYm_Option
      ELSE
        dMdYm => dMdYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, M, dMdD, dMdT, dMdYe, dMdYm, Mn_T, OS_Mn, UnitMn )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, M, Mn_T, OS_Mn, UnitMn )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdYe_Option, dMdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: M(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dMdD_Local, dMdT_Local, dMdYe_Local, dMdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dMdD      , dMdT      , dMdYe      , dMdYm

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdYe_Option ) &
        .OR. PRESENT( dMdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dMdD_Local(nP), dMdT_Local(nP), dMdYe_Local(nP), dMdYm_Local(nP) )

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

      IF( PRESENT( dMdYe_Option ) )THEN
        dMdYe(1:nP) => dMdYe_Option(:)
      ELSE
        dMdYe(1:nP) => dMdYe_Local(:)
      END IF

      IF( PRESENT( dMdYm_Option ) )THEN
        dMdYm(1:nP) => dMdYm_Option(:)
      ELSE
        dMdYm(1:nP) => dMdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, M, dMdD, dMdT, dMdYe, dMdYm, Mn_T, OS_Mn, UnitMn )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, M, Mn_T, OS_Mn, UnitMn )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeProtonMassFraction_TABLE_Scalar &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdYe_Option, dXdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: X
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dXdD_Local, dXdT_Local, dXdYe_Local, dXdYm_Local
    REAL(DP), POINTER :: dXdD      , dXdT      , dXdYe      , dXdYm

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdYe_Option ) &
        .OR. PRESENT( dXdYm_Option )

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

      IF( PRESENT( dXdYe_Option ) )THEN
        dXdYe => dXdYe_Option
      ELSE
        dXdYe => dXdYe_Local
      END IF

      IF( PRESENT( dXdYm_Option ) )THEN
        dXdYm => dXdYm_Option
      ELSE
        dXdYm => dXdYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, X, dXdD, dXdT, dXdYe, dXdYm, Xp_T, OS_Xp, UnitXp )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, X, Xp_T, OS_Xp, UnitXp )

    END IF

  END SUBROUTINE ComputeProtonMassFraction_TABLE_Scalar


  ! THIS SUBROUTINE IS UNTOUCHED, NOT SURE HOW TO HANDLE IT
  SUBROUTINE ComputeProtonMassFraction_TABLE_Vector &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdYe_Option, dXdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: X(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dXdD_Local, dXdT_Local, dXdYe_Local, dXdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dXdD      , dXdT      , dXdYe      , dXdYm

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdYe_Option ) &
        .OR. PRESENT( dXdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dXdD_Local(nP), dXdT_Local(nP), dXdYe_Local(nP), dXdYm_Local(nP) )

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

      IF( PRESENT( dXdYe_Option ) )THEN
        dXdYe(1:nP) => dXdYe_Option(:)
      ELSE
        dXdYe(1:nP) => dXdYe_Local(:)
      END IF

      IF( PRESENT( dXdYm_Option ) )THEN
        dXdYm(1:nP) => dXdYm_Option(:)
      ELSE
        dXdYm(1:nP) => dXdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, X, dXdD, dXdT, dXdYe, dXdYm, Xp_T, OS_Xp, UnitXp )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, X, Xp_T, OS_Xp, UnitXp )

    END IF

  END SUBROUTINE ComputeProtonMassFraction_TABLE_Vector


  SUBROUTINE ComputeNeutronMassFraction_TABLE_Scalar &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdYe_Option, dXdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: X
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dXdD_Local, dXdT_Local, dXdYe_Local, dXdYm_Local
    REAL(DP), POINTER :: dXdD      , dXdT      , dXdYe      , dXdYm

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdYe_Option ) &
        .OR. PRESENT( dXdYm_Option )

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

      IF( PRESENT( dXdYe_Option ) )THEN
        dXdYe => dXdYe_Option
      ELSE
        dXdYe => dXdYe_Local
      END IF

      IF( PRESENT( dXdYm_Option ) )THEN
        dXdYm => dXdYm_Option
      ELSE
        dXdYm => dXdYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, X, dXdD, dXdT, dXdYe, dXdYm, Xn_T, OS_Xn, UnitXn )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, X, Xn_T, OS_Xn, UnitXn )

    END IF

  END SUBROUTINE ComputeNeutronMassFraction_TABLE_Scalar

  SUBROUTINE ComputeNeutronMassFraction_TABLE_Vector &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdYe_Option, dXdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: X(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dXdD_Local, dXdT_Local, dXdYe_Local, dXdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dXdD      , dXdT      , dXdYe      , dXdYm

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdYe_Option ) &
        .OR. PRESENT( dXdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dXdD_Local(nP), dXdT_Local(nP), dXdYe_Local(nP), dXdYm_Local(nP) )

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

      IF( PRESENT( dXdYe_Option ) )THEN
        dXdYe(1:nP) => dXdYe_Option(:)
      ELSE
        dXdYe(1:nP) => dXdYe_Local(:)
      END IF

      IF( PRESENT( dXdYm_Option ) )THEN
        dXdYm(1:nP) => dXdYm_Option(:)
      ELSE
        dXdYm(1:nP) => dXdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, X, dXdD, dXdT, dXdYe, dXdYm, Xn_T, OS_Xn, UnitXn )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, X, Xn_T, OS_Xn, UnitXn )

    END IF

  END SUBROUTINE ComputeNeutronMassFraction_TABLE_Vector


  SUBROUTINE ComputeHeavyMassFraction_TABLE_Scalar &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdYe_Option, dXdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: X
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dXdD_Local, dXdT_Local, dXdYe_Local, dXdYm_Local
    REAL(DP), POINTER :: dXdD      , dXdT      , dXdYe      , dXdYm

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdYe_Option ) &
        .OR. PRESENT( dXdYm_Option )

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

      IF( PRESENT( dXdYe_Option ) )THEN
        dXdYe => dXdYe_Option
      ELSE
        dXdYe => dXdYe_Local
      END IF

      IF( PRESENT( dXdYm_Option ) )THEN
        dXdYm => dXdYm_Option
      ELSE
        dXdYm => dXdYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, X, dXdD, dXdT, dXdYe, dXdYm, Xh_T, OS_Xh, UnitXh )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, X, Xh_T, OS_Xh, UnitXh )

    END IF

  END SUBROUTINE ComputeHeavyMassFraction_TABLE_Scalar

  SUBROUTINE ComputeHeavyMassFraction_TABLE_Vector &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdYe_Option, dXdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: X(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dXdD_Local, dXdT_Local, dXdYe_Local, dXdYm_Local
    REAL(DP), DIMENSION(:), POINTER :: &
      dXdD      , dXdT      , dXdYe      , dXdYm

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdYe_Option ) &
        .OR. PRESENT( dXdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dXdD_Local(nP), dXdT_Local(nP), dXdYe_Local(nP), dXdYm_Local(nP) )

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

      IF( PRESENT( dXdYe_Option ) )THEN
        dXdYe(1:nP) => dXdYe_Option(:)
      ELSE
        dXdYe(1:nP) => dXdYe_Local(:)
      END IF

      IF( PRESENT( dXdYm_Option ) )THEN
        dXdYm(1:nP) => dXdYm_Option(:)
      ELSE
        dXdYm(1:nP) => dXdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, X, dXdD, dXdT, dXdYe, dXdYm, Xh_T, OS_Xh, UnitXh )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, X, Xh_T, OS_Xh, UnitXh )

    END IF

  END SUBROUTINE ComputeHeavyMassFraction_TABLE_Vector

  SUBROUTINE ComputeHeavyMassNumber_TABLE_Scalar &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdYe_Option, dXdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: X
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dXdD_Local, dXdT_Local, dXdYe_Local, dXdYm_Local
    REAL(DP), POINTER :: dXdD      , dXdT      , dXdYe      , dXdYm

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdYe_Option ) &
        .OR. PRESENT( dXdYm_Option )

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

      IF( PRESENT( dXdYe_Option ) )THEN
        dXdYe => dXdYe_Option
      ELSE
        dXdYe => dXdYe_Local
      END IF

      IF( PRESENT( dXdYm_Option ) )THEN
        dXdYm => dXdYm_Option
      ELSE
        dXdYm => dXdYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, X, dXdD, dXdT, dXdYe, dXdYm, Ah_T, OS_Ah, UnitAh )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, X, Ah_T, OS_Ah, UnitAh )

    END IF

  END SUBROUTINE ComputeHeavyMassNumber_TABLE_Scalar

  SUBROUTINE ComputeHeavyMassNumber_TABLE_Vector &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdYe_Option, dXdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: X(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dXdD_Local, dXdT_Local, dXdYe_Local, dXdYm_Local
    REAL(DP), DIMENSION(:), POINTER :: &
      dXdD      , dXdT      , dXdYe      , dXdYm

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdYe_Option ) &
        .OR. PRESENT( dXdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dXdD_Local(nP), dXdT_Local(nP), dXdYe_Local(nP), dXdYm_Local(nP) )

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

      IF( PRESENT( dXdYe_Option ) )THEN
        dXdYe(1:nP) => dXdYe_Option(:)
      ELSE
        dXdYe(1:nP) => dXdYe_Local(:)
      END IF

      IF( PRESENT( dXdYm_Option ) )THEN
        dXdYm(1:nP) => dXdYm_Option(:)
      ELSE
        dXdYm(1:nP) => dXdYm_Local(:)
      END IF
      
      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, X, dXdD, dXdT, dXdYe, dXdYm, Ah_T, OS_Ah, UnitAh )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, X, Ah_T, OS_Ah, UnitAh )

    END IF

  END SUBROUTINE ComputeHeavyMassNumber_TABLE_Vector


  SUBROUTINE ComputeElectronNeutrinoChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, Ym, Munue )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: Munue

    REAL(DP) :: D_P
    REAL(DP) :: Mue, Mup, Mun

#ifdef MICROPHYSICS_WEAKLIB

    D_P  = D / UnitD
    
    IF ( D_P >= Eos_MinD ) THEN
      CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
        ( D, T, Ye, Ym, Mun )
        
      CALL ComputeProtonChemicalPotential_TABLE_Scalar &
        ( D, T, Ye, Ym, Mup )
        
      CALL ComputeElectronChemicalPotential_TABLE_Scalar &
        ( D, T, Ye, Ym, Mue )

      Mue  = Mue
      Mup  = Mup
      Mun  = Mun

      Munue  = ( Mue  + Mup ) - Mun
    ELSE
      Munue = Zero
    END IF

#else

    Munue  = Zero

#endif

  END SUBROUTINE ComputeElectronNeutrinoChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeElectronNeutrinoChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, Munue )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: Munue(1:)

    REAL(DP) :: D_P, T_P, Ye_P, Ym_P
    REAL(DP) :: Mue, Mup, Mun

    INTEGER  :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, T_P, Ye_P, Ym_P, Mue, Mup, Mun ) &
    !$OMP MAP( to: D, T, Ye, Ym, D_T, T_T, Yp_T, OS_Mp, OS_Mn, Mp_T, Mn_T ) &
    !$OMP MAP( from: Munue )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, T_P, Ye_P, Ym_P, Mue, Mup, Mun ) &
    !$ACC COPYIN( D, T, Ye, Ym, D_T, T_T, Yp_T, OS_Mp, OS_Mn, Mp_T, Mn_T ) &
    !$ACC COPYOUT( Munue )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, T_P, Ye_P, Ym_P, Mue, Mup, Mun )
#endif
    DO iP = 1, nP

      D_P  = D(iP)
      T_P  = T(iP)
      Ye_P = Ye(iP)
      Ym_P = Ym(iP)

      IF ( D_P / UnitD >= Eos_MinD ) THEN
        CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
          ( D_P, T_P, Ye_P, Ym_P, Mun )
          
        CALL ComputeProtonChemicalPotential_TABLE_Scalar &
          ( D_P, T_P, Ye_P, Ym_P, Mup )
          
        CALL ComputeElectronChemicalPotential_TABLE_Scalar &
          ( D_P, T_P, Ye_P, Ym_P, Mue )

        Munue(iP)  = ( Mue  + Mup ) - Mun
      ELSE
        Munue(iP) = Zero
      END IF

    END DO

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: Munue )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( Munue )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP
      Munue(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeElectronNeutrinoChemicalPotential_TABLE_Vector

  SUBROUTINE ComputeMuonNeutrinoChemicalPotential_TABLE_Scalar &
      ( D, T, Ye, Ym, Munum )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: Munum

    REAL(DP) :: D_P
    REAL(DP) :: Mum, Mup, Mun

#ifdef MICROPHYSICS_WEAKLIB

    D_P  = D / UnitD

    IF ( D_P >= Eos_MinD ) THEN
      CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
        ( D, T, Ye, Ym, Mun )
        
      CALL ComputeProtonChemicalPotential_TABLE_Scalar &
        ( D, T, Ye, Ym, Mup )
        
      CALL ComputeMuonChemicalPotential_TABLE_Scalar &
        ( D, T, Ye, Ym, Mum )

      Munum  = ( Mum  + Mup ) - Mun
    ELSE
      Munum = Zero
    END IF

#else

    Munum  = Zero

#endif

  END SUBROUTINE ComputeMuonNeutrinoChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeMuonNeutrinoChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, Munum )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: Munum(1:)

    REAL(DP) :: D_P, T_P, Ye_P, Ym_P
    REAL(DP) :: Mup, Mun, Mum

    INTEGER  :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, T_P, Ye_P, Ym_P, Mum, Mup, Mun ) &
    !$OMP MAP( to: D, T, Ye, Ym, D_T, T_T, Yp_T, OS_Mp, OS_Mn, Mp_T, Mn_T ) &
    !$OMP MAP( from: Munum )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, T_P, Ye_P, Ym_P, Mum, Mup, Mun ) &
    !$ACC COPYIN( D, T, Ye, Ym, D_T, T_T, Yp_T, OS_Mp, OS_Mn, Mp_T, Mn_T ) &
    !$ACC COPYOUT( Munum )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, T_P, Ye_P, Ym_P, Mum, Mup, Mun )
#endif
    DO iP = 1, nP

      D_P  = D(iP)
      T_P  = T(iP)
      Ye_P = Ye(iP)
      Ym_P = Ym(iP)

      IF ( D_P >= Eos_MinD ) THEN
        CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
          ( D_P, T_P, Ye_P, Ym_P, Mun )
          
        CALL ComputeProtonChemicalPotential_TABLE_Scalar &
          ( D_P, T_P, Ye_P, Ym_P, Mup )
          
        CALL ComputeMuonChemicalPotential_TABLE_Scalar &
          ( D_P, T_P, Ye_P, Ym_P, Mum )
          
        Munum(iP)  = ( Mum  + Mup ) - Mun
      ELSE
        Munum(iP) = Zero
      END IF

    END DO

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: Munum )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( Munum )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP
      Munum(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeMuonNeutrinoChemicalPotential_TABLE_Vector

  SUBROUTINE ComputeProtonEffectiveMass_TABLE_Scalar &
    ( D, T, Ye, Ym, Emp, dEmpdD_Option, dEmpdT_Option, dEmpdYe_Option, dEmpdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: Emp
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmpdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmpdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmpdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmpdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dEmpdD_Local, dEmpdT_Local, dEmpdYe_Local, dEmpdYm_Local
    REAL(DP), POINTER :: dEmpdD      , dEmpdT      , dEmpdYe, dEmpdYm

    ComputeDerivatives &
      =      PRESENT( dEmpdD_Option ) &
        .OR. PRESENT( dEmpdT_Option ) &
        .OR. PRESENT( dEmpdYe_Option ) &
        .OR. PRESENT( dEmpdYm_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dEmpdD_Option ) )THEN
        dEmpdD => dEmpdD_Option
      ELSE
        dEmpdD => dEmpdD_Local
      END IF

      IF( PRESENT( dEmpdT_Option ) )THEN
        dEmpdT => dEmpdT_Option
      ELSE
        dEmpdT => dEmpdT_Local
      END IF

      IF( PRESENT( dEmpdYe_Option ) )THEN
        dEmpdYe => dEmpdYe_Option
      ELSE
        dEmpdYe => dEmpdYe_Local
      END IF

      IF( PRESENT( dEmpdYm_Option ) )THEN
        dEmpdYm => dEmpdYm_Option
      ELSE
        dEmpdYm => dEmpdYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, Emp, dEmpdD, dEmpdT, dEmpdYe, dEmpdYm, Emp_T, OS_Emp, UnitEmp )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, Emp, Emp_T, OS_Emp, UnitEmp )

    END IF

  END SUBROUTINE ComputeProtonEffectiveMass_TABLE_Scalar

  SUBROUTINE ComputeProtonEffectiveMass_TABLE_Vector &
    ( D, T, Ye, Ym, Emp, dEmpdD_Option, dEmpdT_Option, dEmpdYe_Option, dEmpdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: Emp(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmpdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmpdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmpdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmpdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dEmpdD_Local, dEmpdT_Local, dEmpdYe_Local, dEmpdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dEmpdD      , dEmpdT      , dEmpdYe      , dEmpdYm

    ComputeDerivatives &
      =      PRESENT( dEmpdD_Option ) &
        .OR. PRESENT( dEmpdT_Option ) &
        .OR. PRESENT( dEmpdYe_Option ) &
        .OR. PRESENT( dEmpdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dEmpdD_Local(nP), dEmpdT_Local(nP), dEmpdYe_Local(nP), dEmpdYm_Local(nP) )

      IF( PRESENT( dEmpdD_Option ) )THEN
        dEmpdD(1:nP) => dEmpdD_Option(:)
      ELSE
        dEmpdD(1:nP) => dEmpdD_Local(:)
      END IF

      IF( PRESENT( dEmpdT_Option ) )THEN
        dEmpdT(1:nP) => dEmpdT_Option(:)
      ELSE
        dEmpdT(1:nP) => dEmpdT_Local(:)
      END IF

      IF( PRESENT( dEmpdYe_Option ) )THEN
        dEmpdYe(1:nP) => dEmpdYe_Option(:)
      ELSE
        dEmpdYe(1:nP) => dEmpdYe_Local(:)
      END IF

      IF( PRESENT( dEmpdYm_Option ) )THEN
        dEmpdYm(1:nP) => dEmpdYm_Option(:)
      ELSE
        dEmpdYm(1:nP) => dEmpdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, Emp, dEmpdD, dEmpdT, dEmpdYe, dEmpdYm, Emp_T, OS_Emp, UnitEmp )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, Emp, Emp_T, OS_Emp, UnitEmp )

    END IF

  END SUBROUTINE ComputeProtonEffectiveMass_TABLE_Vector

  
  SUBROUTINE ComputeNeutronEffectiveMass_TABLE_Scalar &
    ( D, T, Ye, Ym, Emn, dEmndD_Option, dEmndT_Option, dEmndYe_Option, dEmndYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: Emn
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmndD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmndT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmndYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmndYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dEmndD_Local, dEmndT_Local, dEmndYe_Local, dEmndYm_Local
    REAL(DP), POINTER :: dEmndD      , dEmndT      , dEmndYe, dEmndYm

    ComputeDerivatives &
      =      PRESENT( dEmndD_Option ) &
        .OR. PRESENT( dEmndT_Option ) &
        .OR. PRESENT( dEmndYe_Option ) &
        .OR. PRESENT( dEmndYm_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dEmndD_Option ) )THEN
        dEmndD => dEmndD_Option
      ELSE
        dEmndD => dEmndD_Local
      END IF

      IF( PRESENT( dEmndT_Option ) )THEN
        dEmndT => dEmndT_Option
      ELSE
        dEmndT => dEmndT_Local
      END IF

      IF( PRESENT( dEmndYe_Option ) )THEN
        dEmndYe => dEmndYe_Option
      ELSE
        dEmndYe => dEmndYe_Local
      END IF

      IF( PRESENT( dEmndYm_Option ) )THEN
        dEmndYm => dEmndYm_Option
      ELSE
        dEmndYm => dEmndYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, Emn, dEmndD, dEmndT, dEmndYe, dEmndYm, Emn_T, OS_Emn, UnitEmn )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, Emn, Emn_T, OS_Emn, UnitEmn )

    END IF

  END SUBROUTINE ComputeNeutronEffectiveMass_TABLE_Scalar

  SUBROUTINE ComputeNeutronEffectiveMass_TABLE_Vector &
    ( D, T, Ye, Ym, Emn, dEmndD_Option, dEmndT_Option, dEmndYe_Option, dEmndYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: Emn(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmndD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmndT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmndYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEmndYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dEmndD_Local, dEmndT_Local, dEmndYe_Local, dEmndYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dEmndD      , dEmndT      , dEmndYe      , dEmndYm

    ComputeDerivatives &
      =      PRESENT( dEmndD_Option ) &
        .OR. PRESENT( dEmndT_Option ) &
        .OR. PRESENT( dEmndYe_Option ) &
        .OR. PRESENT( dEmndYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dEmndD_Local(nP), dEmndT_Local(nP), dEmndYe_Local(nP), dEmndYm_Local(nP) )

      IF( PRESENT( dEmndD_Option ) )THEN
        dEmndD(1:nP) => dEmndD_Option(:)
      ELSE
        dEmndD(1:nP) => dEmndD_Local(:)
      END IF

      IF( PRESENT( dEmndT_Option ) )THEN
        dEmndT(1:nP) => dEmndT_Option(:)
      ELSE
        dEmndT(1:nP) => dEmndT_Local(:)
      END IF

      IF( PRESENT( dEmndYe_Option ) )THEN
        dEmndYe(1:nP) => dEmndYe_Option(:)
      ELSE
        dEmndYe(1:nP) => dEmndYe_Local(:)
      END IF

      IF( PRESENT( dEmndYm_Option ) )THEN
        dEmndYm(1:nP) => dEmndYm_Option(:)
      ELSE
        dEmndYm(1:nP) => dEmndYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, Emn, dEmndD, dEmndT, dEmndYe, dEmndYm, Emn_T, OS_Emn, UnitEmn )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, Emn, Emn_T, OS_Emn, UnitEmn )

    END IF

  END SUBROUTINE ComputeNeutronEffectiveMass_TABLE_Vector

  SUBROUTINE ComputeProtonSelfEnergy_TABLE_Scalar &
    ( D, T, Ye, Ym, Sep, dSepdD_Option, dSepdT_Option, dSepdYe_Option, dSepdYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: Sep
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSepdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSepdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSepdYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSepdYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dSepdD_Local, dSepdT_Local, dSepdYe_Local, dSepdYm_Local
    REAL(DP), POINTER :: dSepdD      , dSepdT      , dSepdYe, dSepdYm

    ComputeDerivatives &
      =      PRESENT( dSepdD_Option ) &
        .OR. PRESENT( dSepdT_Option ) &
        .OR. PRESENT( dSepdYe_Option ) &
        .OR. PRESENT( dSepdYm_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dSepdD_Option ) )THEN
        dSepdD => dSepdD_Option
      ELSE
        dSepdD => dSepdD_Local
      END IF

      IF( PRESENT( dSepdT_Option ) )THEN
        dSepdT => dSepdT_Option
      ELSE
        dSepdT => dSepdT_Local
      END IF

      IF( PRESENT( dSepdYe_Option ) )THEN
        dSepdYe => dSepdYe_Option
      ELSE
        dSepdYe => dSepdYe_Local
      END IF

      IF( PRESENT( dSepdYm_Option ) )THEN
        dSepdYm => dSepdYm_Option
      ELSE
        dSepdYm => dSepdYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, Sep, dSepdD, dSepdT, dSepdYe, dSepdYm, Sep_T, OS_Sep, UnitSep )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, Sep, Sep_T, OS_Sep, UnitSep )

    END IF

  END SUBROUTINE ComputeProtonSelfEnergy_TABLE_Scalar

  SUBROUTINE ComputeProtonSelfEnergy_TABLE_Vector &
    ( D, T, Ye, Ym, Sep, dSepdD_Option, dSepdT_Option, dSepdYe_Option, dSepdYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: Sep(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSepdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSepdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSepdYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSepdYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dSepdD_Local, dSepdT_Local, dSepdYe_Local, dSepdYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dSepdD      , dSepdT      , dSepdYe      , dSepdYm

    ComputeDerivatives &
      =      PRESENT( dSepdD_Option ) &
        .OR. PRESENT( dSepdT_Option ) &
        .OR. PRESENT( dSepdYe_Option ) &
        .OR. PRESENT( dSepdYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dSepdD_Local(nP), dSepdT_Local(nP), dSepdYe_Local(nP), dSepdYm_Local(nP) )

      IF( PRESENT( dSepdD_Option ) )THEN
        dSepdD(1:nP) => dSepdD_Option(:)
      ELSE
        dSepdD(1:nP) => dSepdD_Local(:)
      END IF

      IF( PRESENT( dSepdT_Option ) )THEN
        dSepdT(1:nP) => dSepdT_Option(:)
      ELSE
        dSepdT(1:nP) => dSepdT_Local(:)
      END IF

      IF( PRESENT( dSepdYe_Option ) )THEN
        dSepdYe(1:nP) => dSepdYe_Option(:)
      ELSE
        dSepdYe(1:nP) => dSepdYe_Local(:)
      END IF

      IF( PRESENT( dSepdYm_Option ) )THEN
        dSepdYm(1:nP) => dSepdYm_Option(:)
      ELSE
        dSepdYm(1:nP) => dSepdYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, Sep, dSepdD, dSepdT, dSepdYe, dSepdYm, Sep_T, OS_Sep, UnitSep )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, Sep, Sep_T, OS_Sep, UnitSep )

    END IF

  END SUBROUTINE ComputeProtonSelfEnergy_TABLE_Vector

  SUBROUTINE ComputeNeutronSelfEnergy_TABLE_Scalar &
    ( D, T, Ye, Ym, Sen, dSendD_Option, dSendT_Option, dSendYe_Option, dSendYm_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: Sen
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSendD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSendT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSendYe_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSendYm_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dSendD_Local, dSendT_Local, dSendYe_Local, dSendYm_Local
    REAL(DP), POINTER :: dSendD      , dSendT      , dSendYe, dSendYm

    ComputeDerivatives &
      =      PRESENT( dSendD_Option ) &
        .OR. PRESENT( dSendT_Option ) &
        .OR. PRESENT( dSendYe_Option ) &
        .OR. PRESENT( dSendYm_Option )

    IF( ComputeDerivatives )THEN

      IF( PRESENT( dSendD_Option ) )THEN
        dSendD => dSendD_Option
      ELSE
        dSendD => dSendD_Local
      END IF

      IF( PRESENT( dSendT_Option ) )THEN
        dSendT => dSendT_Option
      ELSE
        dSendT => dSendT_Local
      END IF

      IF( PRESENT( dSendYe_Option ) )THEN
        dSendYe => dSendYe_Option
      ELSE
        dSendYe => dSendYe_Local
      END IF

      IF( PRESENT( dSendYm_Option ) )THEN
        dSendYm => dSendYm_Option
      ELSE
        dSendYm => dSendYm_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, Sen, dSendD, dSendT, dSendYe, dSendYm, Sen_T, OS_Sen, UnitSen )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, Sen, Sen_T, OS_Sen, UnitSen )

    END IF

  END SUBROUTINE ComputeNeutronSelfEnergy_TABLE_Scalar

  SUBROUTINE ComputeNeutronSelfEnergy_TABLE_Vector &
    ( D, T, Ye, Ym, Sen, dSendD_Option, dSendT_Option, dSendYe_Option, dSendYm_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: Sen(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSendD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSendT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSendYe_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dSendYm_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dSendD_Local, dSendT_Local, dSendYe_Local, dSendYm_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dSendD      , dSendT      , dSendYe      , dSendYm

    ComputeDerivatives &
      =      PRESENT( dSendD_Option ) &
        .OR. PRESENT( dSendT_Option ) &
        .OR. PRESENT( dSendYe_Option ) &
        .OR. PRESENT( dSendYm_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dSendD_Local(nP), dSendT_Local(nP), dSendYe_Local(nP), dSendYm_Local(nP) )

      IF( PRESENT( dSendD_Option ) )THEN
        dSendD(1:nP) => dSendD_Option(:)
      ELSE
        dSendD(1:nP) => dSendD_Local(:)
      END IF

      IF( PRESENT( dSendT_Option ) )THEN
        dSendT(1:nP) => dSendT_Option(:)
      ELSE
        dSendT(1:nP) => dSendT_Local(:)
      END IF

      IF( PRESENT( dSendYe_Option ) )THEN
        dSendYe(1:nP) => dSendYe_Option(:)
      ELSE
        dSendYe(1:nP) => dSendYe_Local(:)
      END IF

      IF( PRESENT( dSendYm_Option ) )THEN
        dSendYm(1:nP) => dSendYm_Option(:)
      ELSE
        dSendYm(1:nP) => dSendYm_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, Sen, dSendD, dSendT, dSendYe, dSendYm, Sen_T, OS_Sen, UnitSen )

    ELSE

      CALL ComputeDependentVariableBaryons_TABLE_Vector &
             ( D, T, Ye, Ym, Sen, Sen_T, OS_Sen, UnitSen )

    END IF

  END SUBROUTINE ComputeNeutronSelfEnergy_TABLE_Vector

  SUBROUTINE ComputeDependentVariableBaryons_TABLE_Scalar &
    ( D, T, Ye, Ym, V, V_T, OS_V, Units_V )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: V
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P

    IF ( D_P >= Eos_MinD ) THEN
      CALL LogInterpolateSingleVariable_3D_Custom_Point &
            ( D_P, T_P, Yp_P, D_T, T_T, Yp_T, OS_V, V_T, V_P )
    ELSE
      V_P = Zero
    END IF
    V = V_P * Units_V

#else

    V = Zero

#endif

  END SUBROUTINE ComputeDependentVariableBaryons_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableBaryons_TABLE_Vector &
    ( D, T, Ye, Ym, V, V_T, OS_V, Units_V )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: V(1:)
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    INTEGER  :: iP, nP
    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P ) &
      !$OMP MAP( to: D, T, Ye, Ym, D_T, T_T, Yp_T, OS_V, V_T ) &
      !$OMP MAP( from: V )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P ) &
      !$ACC COPYIN( D, T, Ye, Ym, D_T, T_T, Yp_T, OS_V, V_T ) &
      !$ACC COPYOUT( V )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P )
#endif
    DO iP = 1, nP

      D_P = D(iP) / UnitD
      T_P = T(iP) / UnitT
      Ye_P = Ye(iP) / UnitY
      Ym_P = Ym(iP) / UnitY
      Yp_P = Ye_P + Ym_P

      IF ( D_P >= Eos_MinD ) THEN
        CALL LogInterpolateSingleVariable_3D_Custom_Point &
              ( D_P, T_P, Yp_P, D_T, T_T, Yp_T, OS_V, V_T, V_P )

      V(iP) = V_P * Units_V
      ELSE
        V(iP) = Zero
      END IF

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

  END SUBROUTINE ComputeDependentVariableBaryons_TABLE_Vector

  SUBROUTINE ComputePressureBaryons_TABLE_Scalar &
    ( D, T, Ye, Ym, P )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: P

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P
    REAL(DP) :: LocalOffset
    REAL(DP) :: D_cube(2), T_cube(2), Yp_cube(2), P_cube(2,2,2)
    INTEGER  :: iD, iT, iYp
    INTEGER  :: SizeDs, SizeTs, SizeYps

#ifdef MICROPHYSICS_WEAKLIB

    D_P  = D  / UnitD
    T_P  = T  / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P

    ! Bracekt points and find the cube in the table
    SizeDs = SIZE( D_T )
    SizeTs = SIZE( T_T )
    SizeYps = SIZE( Yp_T )

    iD  = Index1D_Log( D_P , D_T  )
    iT  = Index1D_Log( T_P , T_T  )
    iYp = Index1D_Lin( Yp_P, Yp_T )

    iD  = MIN( MAX( 1, iD ) , SizeDs - 1  )
    iT  = MIN( MAX( 1, iT ) , SizeTs - 1  )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    P_cube  = P_T (iD:iD+1,iT:iT+1,iYp:iYp+1)
    D_cube  = D_T (iD:iD+1)
    T_cube  = T_T (iT:iT+1)
    Yp_cube = Yp_T(iYp:iYp+1)
    
    LocalOffset = MINVAL(P_cube)
    IF (LocalOffset .lt. 0.0_dp) THEN
        LocalOffset = -1.1d0*LocalOffset
    ELSE
        LocalOffset = 0.0_dp
    ENDIF

    IF ( D_P >= Eos_MinD ) THEN
      CALL LogInterpolateSingleVariable_3D_Custom_Point &
            ( D_P, T_P, Yp_P, D_cube, T_cube, Yp_cube, LocalOffset, &
            LOG10(P_cube + LocalOffset), P )
    ELSE
      P = 0.0d0
    ENDIF

#else

    P = Zero

#endif

  END SUBROUTINE ComputePressureBaryons_TABLE_Scalar

  SUBROUTINE ComputeDependentVariableTotal_TABLE_Scalar &
    ( D, T, Ye, Ym, V, V_T, OS_V, Units_V, &
    ReturnP, ReturnE, ReturnS )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: V
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    INTEGER,  INTENT(in)  :: ReturnP, ReturnE, ReturnS

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: Vtot(2,2,2), E_LeptPhot(2,2,2), &
                P_LeptPhot(2,2,2), S_LeptPhot(2,2,2)
    
    INTEGER  :: i, iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps

    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
    
#ifdef MICROPHYSICS_WEAKLIB
    
    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P

    IF ( D_P < Eos_MinD ) THEN
      V = Zero
      RETURN
    ENDIF

    Ye_over_Yp = Ye_P / Yp_P
    Ym_over_Yp = Ym_P / Yp_P

    ! Now bracket the points 
    SizeDs  = SIZE( D_T  )
    SizeTs  = SIZE( T_T  )
    SizeYps = SIZE( Yp_T )

    iD  = Index1D_Log( D_P , D_T  )
    iT  = Index1D_Log( T_P , T_T  )
    iYp = Index1D_Lin( Yp_P, Yp_T )

    iD  = MIN( MAX( 1, iD ) , SizeDs - 1  )
    iT  = MIN( MAX( 1, iT ) , SizeTs - 1  )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    DO iL_T=1,2
      DO iL_D=1,2
        DO iL_Y=1,2
          ElectronPhotonState % t   = T_T (iT +iL_T-1)
          ElectronPhotonState % rho = D_T (iD +iL_D-1)
          ElectronPhotonState % ye  = Yp_T(iYp+iL_Y-1) * Ye_over_Yp
          
          CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

          MuonState % t     = T_T(iT+iL_T-1)
          MuonState % rhoym = D_T(iD+iL_D-1) * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
          
          CALL FullMuonEOS(MuonTable, MuonState)
          
          E_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % e + MuonState % e 
          P_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % p + MuonState % p
          S_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % s + MuonState % s

        END DO
      END DO
    END DO
    
    IF ( ReturnP == 1 ) THEN
      Vtot(:,:,:) = V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + P_LeptPhot(:,:,:)
    ELSE IF ( ReturnE == 1 ) THEN
      Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_LeptPhot(:,:,:)
    ELSE
      Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_LeptPhot(:,:,:)
    ENDIF
    
    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
           OS_V, LOG10(Vtot), V_P )

    V = V_P * Units_V

#else

    V = Zero

#endif

  END SUBROUTINE ComputeDependentVariableTotal_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableTotal_TABLE_Vector &
    ( D, T, Ye, Ym, V, V_T, OS_V, Units_V, &
    ReturnP, ReturnE, ReturnS )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: V(1:)
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    INTEGER,  INTENT(in)  :: ReturnP, ReturnE, ReturnS

    INTEGER  :: iP, nP
    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: Vtot(2,2,2), E_LeptPhot(2,2,2), &
                P_LeptPhot(2,2,2), S_LeptPhot(2,2,2)
    
    INTEGER  :: i, iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps

    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
        
    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P ) &
      !$OMP MAP( to: D, T, Ye, Ym, D_T, T_T, Yp_T, OS_V, V_T ) &
      !$OMP MAP( from: V )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P ) &
      !$ACC COPYIN( D, T, Ye, Ym, D_T, T_T, Yp_T, OS_V, V_T ) &
      !$ACC COPYOUT( V )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P )
#endif
    DO iP = 1, nP

      D_P  = D(iP)  / UnitD
      T_P  = T(iP)  / UnitT
      Ye_P = Ye(iP) / UnitY
      Ym_P = Ym(iP) / UnitY
      Yp_P = Ye_P + Ym_P
      IF ( D_P < Eos_MinD ) THEN
        V(iP) = Zero
        CYCLE
      ENDIF

      Ye_over_Yp = Ye_P / Yp_P
      Ym_over_Yp = Ym_P / Yp_P
      
      ! Now bracket the points 
      SizeDs  = SIZE( D_T )
      SizeTs  = SIZE( T_T )
      SizeYps = SIZE( Yp_T )

      iD  = Index1D_Log( D_P,  D_T  )
      iT  = Index1D_Log( T_P,  T_T  )
      iYp = Index1D_Lin( Yp_P, Yp_T )

      iD  = MIN( MAX( 1, iD  ), SizeDs  - 1 )
      iT  = MIN( MAX( 1, iT  ), SizeTs  - 1 )
      iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )
      
      ! now calculate muon and electron contribution on that specific point
      DO iL_T=1,2
        DO iL_D=1,2
          DO iL_Y=1,2
            ElectronPhotonState % t   = T_T (iT +iL_T-1)
            ElectronPhotonState % rho = D_T (iD +iL_D-1)
            ElectronPhotonState % ye  = Yp_T(iYp+iL_Y-1) * Ye_over_Yp
            
            ! CALL FullHelmEOS(1, HelmTable, ElectronPhotonState, .false., .false.)
            CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

            MuonState % t     = T_T(iT+iL_T-1)
            MuonState % rhoym = D_T(iD+iL_D-1) * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
            
            CALL FullMuonEOS(MuonTable, MuonState)
            
            E_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % e + MuonState % e
            P_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % p + MuonState % p
            S_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % s + MuonState % s
            
          END DO
        END DO
      END DO

      IF ( ReturnP == 1 ) THEN
        Vtot(:,:,:) = V_T(iD:iD+1,iT:iT+1,iYp:iYp+1)       + P_LeptPhot(:,:,:)
      ELSE IF ( ReturnE == 1 ) THEN
        Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_LeptPhot(:,:,:)
      ELSE
        Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_LeptPhot(:,:,:)
      ENDIF

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_V, LOG10(Vtot), V_P )

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

  END SUBROUTINE ComputeDependentVariableTotal_TABLE_Vector

  SUBROUTINE ComputeTotalPES_TABLE_Scalar &
    ( D, T, Ye, Ym, P, E, S )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: P, E, S

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, P_P, E_P, S_P
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: Ptot(2,2,2), Etot(2,2,2), Stot(2,2,2), &
                E_LeptPhot(2,2,2), P_LeptPhot(2,2,2), S_LeptPhot(2,2,2)
    
    INTEGER  :: i, iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps

    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P

    IF ( D_P < Eos_MinD ) THEN
      P = Zero
      E = Zero
      S = Zero
      RETURN
    ENDIF

    Ye_over_Yp = Ye_P / Yp_P
    Ym_over_Yp = Ym_P / Yp_P

    ! Now bracket the points 
    SizeDs  = SIZE( D_T  )
    SizeTs  = SIZE( T_T  )
    SizeYps = SIZE( Yp_T )

    iD  = Index1D_Log( D_P , D_T  )
    iT  = Index1D_Log( T_P , T_T  )
    iYp = Index1D_Lin( Yp_P, Yp_T )

    iD  = MIN( MAX( 1, iD ) , SizeDs - 1  )
    iT  = MIN( MAX( 1, iT ) , SizeTs - 1  )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    DO iL_T=1,2
      DO iL_D=1,2
        DO iL_Y=1,2
          ElectronPhotonState % t   = T_T (iT +iL_T-1)
          ElectronPhotonState % rho = D_T (iD +iL_D-1)
          ElectronPhotonState % ye  = Yp_T(iYp+iL_Y-1) * Ye_over_Yp
          
          ! CALL FullHelmEOS(1, HelmTable, ElectronPhotonState, .false., .false.)
          CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

          MuonState % t     = T_T(iT+iL_T-1)
          MuonState % rhoym = D_T(iD+iL_D-1) * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
          
          CALL FullMuonEOS(MuonTable, MuonState)
          
          E_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % e + MuonState % e
          P_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % p + MuonState % p
          S_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % s + MuonState % s

        END DO
      END DO
    END DO

    Ptot(:,:,:) = P_T(iD:iD+1,iT:iT+1,iYp:iYp+1)          + P_LeptPhot(:,:,:)
    Etot(:,:,:) = 10.0_DP**E_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_LeptPhot(:,:,:)
    Stot(:,:,:) = 10.0_DP**S_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_LeptPhot(:,:,:)

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
           OS_P, LOG10(Ptot), P_P )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
           OS_E, LOG10(Etot), E_P )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
           OS_S, LOG10(Stot), S_P )

    P = P_P * UnitP
    E = E_P * UnitE
    S = S_P * UnitS

#else

    P = Zero
    E = Zero
    S = Zero

#endif

  END SUBROUTINE ComputeTotalPES_TABLE_Scalar


  SUBROUTINE ComputeTotalPES_TABLE_Vector &
    ( D, T, Ye, Ym, P, E, S )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: E(1:), P(1:), S(1:)

    INTEGER  :: iP, nP
    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, P_P, E_P, S_P
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: Ptot(2,2,2), Etot(2,2,2), Stot(2,2,2), &
                P_LeptPhot(2,2,2), E_LeptPhot(2,2,2), S_LeptPhot(2,2,2)
    
    INTEGER  :: i, iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps

    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
        
    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, P_P, E_P, S_P ) &
      !$OMP MAP( to: D, T, Ye, Ym, D_T, T_T, Yp_T ) &
      !$OMP MAP( from: P, E, S )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, P_P, E_P, S_P ) &
      !$ACC COPYIN( D, T, Ye, Ym, D_T, T_T, Yp_T ) &
      !$ACC COPYOUT( P, E, S )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, P_P, E_P, S_P )
#endif
    DO iP = 1, nP

      D_P  = D(iP)  / UnitD
      T_P  = T(iP)  / UnitT
      Ye_P = Ye(iP) / UnitY
      Ym_P = Ym(iP) / UnitY
      Yp_P = Ye_P + Ym_P

      IF ( D_P < Eos_MinD ) THEN
        P(iP) = Zero
        E(iP) = Zero
        S(iP) = Zero
        CYCLE
      ENDIF
      
      Ye_over_Yp = Ye_P / Yp_P
      Ym_over_Yp = Ym_P / Yp_P
      
      ! Now bracket the points 
      SizeDs  = SIZE( D_T )
      SizeTs  = SIZE( T_T )
      SizeYps = SIZE( Yp_T )

      iD  = Index1D_Log( D_P,  D_T  )
      iT  = Index1D_Log( T_P,  T_T  )
      iYp = Index1D_Lin( Yp_P, Yp_T )

      iD  = MIN( MAX( 1, iD  ), SizeDs  - 1 )
      iT  = MIN( MAX( 1, iT  ), SizeTs  - 1 )
      iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )
      
      ! now calculate muon and electron contribution on that specific point
      DO iL_T=1,2
        DO iL_D=1,2
          DO iL_Y=1,2
            ElectronPhotonState % t   = T_T (iT +iL_T-1)
            ElectronPhotonState % rho = D_T (iD +iL_D-1)
            ElectronPhotonState % ye  = Yp_T(iYp+iL_Y-1) * Ye_over_Yp
            
            CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

            MuonState % t     = T_T(iT+iL_T-1)
            MuonState % rhoym = D_T(iD+iL_D-1) * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
            
            CALL FullMuonEOS(MuonTable, MuonState)
            
            E_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % e + MuonState % e
            P_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % p + MuonState % p
            S_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % s + MuonState % s
            
          END DO
        END DO
      END DO

      Ptot(:,:,:) = P_T(iD:iD+1,iT:iT+1,iYp:iYp+1)          + P_LeptPhot(:,:,:)
      Etot(:,:,:) = 10.0_DP**E_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_LeptPhot(:,:,:)
      Stot(:,:,:) = 10.0_DP**S_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_LeptPhot(:,:,:)

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_P, LOG10(Ptot), P_P )

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_E, LOG10(Etot), E_P )

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_S, LOG10(Stot), S_P )

      P(iP) = P_P * UnitP
      E(iP) = E_P * UnitE
      S(iP) = S_P * UnitS

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
      P(iP) = Zero
      E(iP) = Zero
      S(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeTotalPES_TABLE_Vector

  SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
    ( D, T, Ye, Ym, V, dVdD, dVdT, dVdYe, dVdYm, V_T, OS_V, Units_V, &
    ReturnP, ReturnE, ReturnS)

  ! CHECK THIS I THINK IT'S WRONG BUT IT'S NOT REALLY USED

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: V, dVdD, dVdT, dVdYe, dVdYm
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    INTEGER , INTENT(in)  :: ReturnP, ReturnE, ReturnS

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P, dV_P(3), dVbary_P(3)
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: Vtot(2,2,2), E_LeptPhot(2,2,2), &
                P_LeptPhot(2,2,2), S_LeptPhot(2,2,2)
    
    INTEGER  :: i, iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps
    REAL(DP) :: dD, dT
	  REAL(DP) :: aD, aT
    REAL(DP) :: dVdrhoym, dVdummy, Vdummy, LocalOffset
    
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P

    IF ( D_P < Eos_MinD ) THEN
      V     = Zero
      dVdD  = Zero
      dVdT  = Zero
      dVdYe = Zero
      dVdYm = Zero
      RETURN
    ENDIF

    Ye_over_Yp = Ye_P/(Ye_P + Ym_P)
    Ym_over_Yp = Ym_P/(Ye_P + Ym_P)

    ! Now bracket the points 
    SizeDs = SIZE( D_T )
    SizeTs = SIZE( T_T )
    SizeYps = SIZE( Yp_T )

    iD = Index1D_Log( D_P, D_T )
    iT = Index1D_Log( T_P, T_T )
    iYp = Index1D_Lin( Yp_P, Yp_T )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )
    
    ! now calculate muon and electron contribution on that specific point
    DO iL_T=1,2
      DO iL_D=1,2
        DO iL_Y=1,2
          ElectronPhotonState % t   = T_T (iT +iL_T-1)
          ElectronPhotonState % rho = D_T (iD +iL_D-1)
          ElectronPhotonState % ye  = Yp_T(iYp+iL_Y-1) * Ye_over_Yp
          
          CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

          MuonState % t = T_T(iT+iL_T-1)
          MuonState % rhoym = D_T(iD+iL_D-1) * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
          
          CALL FullMuonEOS(MuonTable, MuonState)
          
          E_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % e + MuonState % e
          P_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % p + MuonState % p
          S_LeptPhot(iL_D,iL_T,iL_Y) = ElectronPhotonState % s + MuonState % s

        END DO
      END DO
    END DO

    ! Now you have to calculate dVbarydYp, and if you have pressure, you
    ! have to handle the offset accordingly
    IF (ReturnP .EQ. 1) THEN
      Vtot(:,:,:) = V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + P_LeptPhot(:,:,:)
   
      LocalOffset = MINVAL(V_T(iD:iD+1,iT:iT+1,iYp:iYp+1))
      IF (LocalOffset .lt. 0.0_dp) THEN
          LocalOffset = -1.1d0*LocalOffset
      ELSE
          LocalOffset = 0.0_dp
      ENDIF
      
      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             LocalOffset, LOG10(V_T + LocalOffset), V_P, dVbary_P )
             
      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_P, LOG10(Vtot), V_P, dV_P )

      ! Now calculate electron and muon derivatives
      ! Initialize temperature, density, ye
      ElectronPhotonState % t   = T_P
      ElectronPhotonState % rho = D_P
      ElectronPhotonState % ye  = Ye_P

      ! calculate electron quantities
      CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

      ! dVdYe = d(Vbary + Vele + Vphot Vmu)dYe = dVbarydYe + d(Vele + Vphot)dYe = dVbarydYp + dVeledYe
      dVdYe = dVbary_P(3) + D_P/Ye_P * ElectronPhotonState % dpdr
      
      ! NOW MUONS ---------------------------------------------
      IF ( T_P .lt. MuonTable % t(1) ) THEN
        
        dVdYm = 0.0_dp
      
      ELSE

        CALL GetIndexAndDelta_Log( T_P, MuonTable % t(:), iT, dT )
        IF ( D_P*Ym_P .lt. MuonTable % rhoym(iT,i) ) THEN
          
          dVdYm = 0.0_dp

        ELSE

          DO i = 1, MuonTable % nPointsDen
            IF (MuonTable % rhoym(iT,i) >= D_P*Ym_P) THEN
              iD = i
              dD = LOG10( D_P*Ym_P / MuonTable % rhoym(iT,i) ) / LOG10( MuonTable % rhoym(iT,i+1) / MuonTable % rhoym(iT,i) )
              EXIT
            ENDIF 
          END DO

          aD = 1.0_dp / ( D_P * LOG10( MuonTable % rhoym(iT,iD+1) / MuonTable % rhoym(iT,iD) ) )
          aT = 1.0_dp / ( T_P * LOG10( MuonTable % t(iT+1) / MuonTable % t(iT) ) )
          
          CALL LinearInterpDeriv_Array_Point &
                ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % p), Vdummy, &
                  dVdrhoym, dVdummy )
          
          dVdYm = dVbary_P(3) + dVdrhoym * D_P ! make sure the derivative is wr2 ym, not rhoym
          
        ENDIF

      END IF

    ELSE IF (ReturnE .EQ. 1) THEN

      Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_LeptPhot(:,:,:)

      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_V, LOG10(V_T), V_P, dVbary_P )
             
      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_V, LOG10(Vtot), V_P, dV_P )
             
      ! Now calculate electron and muon derivatives
      ! Initialize temperature, density, ye
      ElectronPhotonState % t   = T_P
      ElectronPhotonState % rho = D_P
      ElectronPhotonState % ye  = Ye_P

      ! calculate electron quantities
      CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

      dVdYe = dVbary_P(3) + D_P/Ye_P * ElectronPhotonState % dedr ! This is valid because dPdrho for photons is zero

      ! NOW MUONS ---------------------------------------------
      IF ( T_P .lt. MuonTable % t(1) ) THEN
        
        dVdYm = 0.0_dp
      
      ELSE

        CALL GetIndexAndDelta_Log( T_P, MuonTable % t(:), iT, dT )
        IF ( D_P*Ym_P .lt. MuonTable % rhoym(iT,i) ) THEN
          
          dVdYm = 0.0_dp

        ELSE

          DO i = 1, MuonTable % nPointsDen
            IF (MuonTable % rhoym(iT,i) >= D_P*Ym_P) THEN
              iD = i
              dD = LOG10( D_P*Ym_P / MuonTable % rhoym(iT,i) ) / LOG10( MuonTable % rhoym(iT,i+1) / MuonTable % rhoym(iT,i) )
              EXIT
            ENDIF 
          END DO

          CALL LinearInterpDeriv_Array_Point &
                ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % e), Vdummy, &
                  dVdrhoym, dVdummy )
          
          dVdYm = dVbary_P(3) + dVdrhoym * D_P ! make sure the derivative is wr2 ym, not rhoym
        END IF
      END IF
             
    ELSE

      Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_LeptPhot(:,:,:)

      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_V, LOG10(V_T), V_P, dVbary_P )
             
      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_V, LOG10(Vtot), V_P, dV_P )
             
      ! Now calculate electron and muon derivatives
      ! Initialize temperature, density, ye
      ElectronPhotonState % t   = T_P
      ElectronPhotonState % rho = D_P
      ElectronPhotonState % ye  = Ye_P

      ! calculate electron quantities
      CALL ElectronPhotonEOS(HelmTable, ElectronPhotonState)

      dVdYe = dVbary_P(3) + D_P/Ye_P * ElectronPhotonState % dsdr
      WRITE(*,*) 'Need to get this properly from Electron EOS, careful'
      STOP
      
      ! NOW MUONS ---------------------------------------------
      IF ( T_P .lt. MuonTable % t(1) ) THEN
        
        dVdYm = 0.0_dp
      
      ELSE

        CALL GetIndexAndDelta_Log( T_P, MuonTable % t(:), iT, dT )
        IF ( D_P*Ym_P .lt. MuonTable % rhoym(iT,i) ) THEN
          
          dVdYm = 0.0_dp

        ELSE

          DO i = 1, MuonTable % nPointsDen
            IF (MuonTable % rhoym(iT,i) >= D_P*Ym_P) THEN
              iD = i
              dD = LOG10( D_P*Ym_P / MuonTable % rhoym(iT,i) ) / LOG10( MuonTable % rhoym(iT,i+1) / MuonTable % rhoym(iT,i) )
              EXIT
            ENDIF 
          END DO

          CALL LinearInterpDeriv_Array_Point &
                ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % s), Vdummy, &
                  dVdrhoym, dVdummy )
          
          dVdYm = dVbary_P(3) + dVdrhoym * D_P ! make sure the derivative is wr2 ym, not rhoym
        END IF
      END IF
             
    END IF
    
    CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
       ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
       OS_V, Vtot, V_P, dV_P )
             
    V = V_P * Units_V

    dVdD = dV_P(1) * Units_V / UnitD
    dVdT = dV_P(2) * Units_V / UnitT

#else

    V     = Zero
    dVdD  = Zero
    dVdT  = Zero
    dVdYe = Zero
    dVdYm = Zero

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Vector &
    ( D, T, Ye, Ym, V, dVdD, dVdT, dVdYe, dVdYm, V_T, OS_V, Units_V, &
    ReturnP, ReturnE, ReturnS )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: V(1:), dVdD(1:), dVdT(1:), dVdYe(1:), dVdYm(1:)
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    INTEGER , INTENT(in)  :: ReturnP, ReturnE, ReturnS

    INTEGER :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP MAP( to: D, T, Ye, Ym, OS_V, V_T ) &
      !$OMP MAP( from: V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Ye, Ym, OS_V, V_T ) &
      !$ACC COPYOUT( V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
          ( D(iP), T(iP), Ye(iP), Ym(iP), V(iP), dVdD(iP), dVdT(iP), &
            dVdYe(iP), dVdYm(iP), V_T, OS_V, Units_V, ReturnP, ReturnE, ReturnS )

    END DO

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP
      V   (iP) = Zero
      dVdD(iP) = Zero
      dVdT(iP) = Zero
      dVdYe(iP) = Zero
      dVdYm(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Vector

  SUBROUTINE ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
    ( D, T, Ye, Ym, V, dVdD, dVdT, dVdYe, dVdYm, V_T, OS_V, Units_V )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: V, dVdD, dVdT, dVdYe, dVdYm
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P, dV_P(3)

    INTEGER  :: i, iD, iT, iYp, iL_D, iL_Y, iL_T
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P

    IF ( D_P >= Eos_MinD ) THEN

      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
            ( D_P, T_P, Yp_P, D_T, T_T, Yp_T, OS_V, V_T, V_P, dV_P )
              
      V = V_P * Units_V

      dVdD = dV_P(1) * Units_V / UnitD
      dVdT = dV_P(2) * Units_V / UnitT
      dVdYe = dV_P(3) * Units_V / UnitY
      dVdYm = dV_P(3) * Units_V / UnitY
    
    ELSE

      V    = Zero
      dVdD = Zero
      dVdT = Zero
      dVdYe = Zero
      dVdYm = Zero

    ENDIF

#else

    V    = Zero
    dVdD = Zero
    dVdT = Zero
    dVdYe = Zero
    dVdYm = Zero

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
    ( D, T, Ye, Ym, V, dVdD, dVdT, dVdYe, dVdYm, V_T, OS_V, Units_V )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: V(1:), dVdD(1:), dVdT(1:), dVdYe(1:), dVdYm(1:)
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    INTEGER :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP MAP( to: D, T, Ye, Ym, OS_V, V_T ) &
      !$OMP MAP( from: V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Ye, Ym, OS_V, V_T ) &
      !$ACC COPYOUT( V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D(iP), T(iP), Ye(iP), Ym(iP), V(iP), dVdD(iP), dVdT(iP), dVdYe(iP), &
               dVdYm(iP), V_T, OS_V, Units_V )

    END DO

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP
      V   (iP) = Zero
      dVdD(iP) = Zero
      dVdT(iP) = Zero
      dVdYe(iP) = Zero
      dVdYm(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector
  
END MODULE EquationOfStateModule_TABLE
