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
#ifdef INVERSION_COMBINED
  USE wlEOSComponentsCombinedInversionModule, ONLY: &
    InitializeEOSComponentsInversion, &
    ComputeTemperatureWith_DEYpYl_Single_Guess_Error, &
    ComputeTemperatureWith_DEYpYl_Single_Guess_NoError, &
    ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error, &
    ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError, &
    ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error, &
    DescribeEOSComponentsInversionError
#else
  USE wlEOSComponentsSeparateInversionModule, ONLY: &
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
    HelmholtzEOSType, MuonEOSType
  USE wlElectronEOS, ONLY: &
    FullHelmEOS, MinimalHelmEOS_rt, ElectronStateType
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
    Kelvin, &
    Dyne, &
    Erg, &
    MeV

  IMPLICIT NONE
  PRIVATE

  CHARACTER(256) :: &
    EquationOfStateTableName
  INTEGER :: &
    iRho, iTemp, iYp, &
    iP_T, iS_T, iE_T, iMe_T, iMp_T, iMn_T, &
    iXp_T, iXn_T, iXa_T, iXh_T, iAh_T
  REAL(DP) :: &
    UnitD, UnitT, UnitY, &
    UnitP, UnitS, UnitE, UnitMl, UnitMp, UnitMn, &
    UnitXp, UnitXn, UnitXa, UnitXh, UnitAh
  REAL(DP), PUBLIC :: &
    OS_P, OS_S, OS_E, OS_Mp, OS_Mn, &
    OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Ah
  REAL(DP), PARAMETER :: &
    BaryonMass = AtomicMassUnit
  REAL(DP) :: minvar, OS_loc
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    D_T, T_T, Yp_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    P_T, S_T, E_T, Mp_T, Mn_T, &
    Xp_T, Xn_T, Xa_T, Xh_T, Ah_T
#ifdef MICROPHYSICS_WEAKLIB
  LOGICAL :: UsingExternalEOS
  TYPE(EquationOfStateTableType), POINTER :: EOS
  TYPE(HelmholtzEOSType), POINTER :: HelmholtzTable
  TYPE(MuonEOSType), POINTER :: MuonTable
#endif

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

  INTERFACE ComputeDependentVariableAndDerivativesTotal_TABLE
    MODULE PROCEDURE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar
    MODULE PROCEDURE ComputeDependentVariableAndDerivativesTotal_TABLE_Vector
  END INTERFACE ComputeDependentVariableAndDerivativesTotal_TABLE

  ! Define local constants
  
  REAL(dp), PARAMETER  :: ergmev    = 1.602177d-6   ! ergs per MeV
  REAL(dp), PARAMETER  :: avn       = 6.022141d+23  ! avogadro's number
  REAL(dp), PARAMETER  :: rmu       = 1.0d0/avn     ! atomic mass unit
  REAL(dp), PARAMETER  :: mass_ele  = 0.510998d+00  ! electron mass [MeV]
  REAL(dp), PARAMETER  :: mass_mu   = 105.65837d+00 ! muon mass [MeV]

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( D_T, T_T, Yp_T, &
  !$OMP   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, UnitMl, UnitMp, UnitMn, &
  !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitAh, &
  !$OMP   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$OMP   OS_Xa, OS_Xh, OS_Ah, &
  !$OMP   P_T, S_T, E_T, Mp_T, Mn_T, Xp_T, Xn_T, &
  !$OMP   Xa_T, Xh_T, Ah_T, &
  !$OMP   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( D_T, T_T, Yp_T, &
  !$ACC   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, UnitMl, UnitMp, UnitMn, &
  !$ACC   UnitXp, UnitXn, UnitXa, UnitXh, UnitAh, &
  !$ACC   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$ACC   OS_Xa, OS_Xh, OS_Ah, &
  !$ACC   P_T, S_T, E_T, Mp_T, Mn_T, Xp_T, Xn_T, &
  !$ACC   Xa_T, Xh_T, Ah_T, &
  !$ACC   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
#endif

CONTAINS

  SUBROUTINE InitializeEquationOfState_TABLE &
    ( EquationOfStateTableName_Option, UseChemicalPotentialShift_Option, &
      Verbose_Option, External_EOS )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    LOGICAL,          INTENT(in), OPTIONAL :: UseChemicalPotentialShift_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option
#ifdef MICROPHYSICS_WEAKLIB
    TYPE(EquationOfStateTableType), POINTER, &
                      INTENT(in), OPTIONAL :: External_EOS
#else
    INTEGER,          INTENT(in), OPTIONAL :: External_EOS
#endif
    
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Ps, Ss, Es
    LOGICAL :: UseChemicalPotentialShift, Verbose
    
    REAL(DP) :: Eele, Pele, Sele

    ! Added variables to handle Combined EOSs
    TYPE(ElectronStateType) :: ElectronState

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
       ALLOCATE( HelmholtzTable)
       ALLOCATE( MuonTable )
       UsingExternalEOS = .FALSE.

       CALL InitializeHDF( )
        
       ! read in helmholtz table
       CALL ReadHelmholtzTableHDF( HelmholtzTable, TRIM( EquationOfStateTableName ) )

       ! read in helmholtz table
       CALL ReadMuonTableHDF( MuonTable, TRIM( EquationOfStateTableName ) )

       ! read in baryon table -------------------------------
       CALL ReadEquationOfStateTableHDF( EOS, TRIM( EquationOfStateTableName ) )
        
       CALL FinalizeHDF( )

    ELSE

       EOS => External_EOS
       UsingExternalEOS = .TRUE.

    END IF

    ! --- Thermodynamic State Indices ---

    iRho  = EOS % TS % Indices % iRho
    iTemp = EOS % TS % Indices % iT
    iYp   = EOS % TS % Indices % iYe

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

    ! --- Thermodynamic States ---

    ALLOCATE( D_T(EOS % TS % nPoints(iRho)) )
    D_T = EOS % TS % States(iRho) % Values

    Min_D = MINVAL( D_T ) * Gram / Centimeter**3
    Max_D = MAXVAL( D_T ) * Gram / Centimeter**3

    ALLOCATE( T_T(EOS % TS % nPoints(iTemp)) )
    T_T = EOS % TS % States(iTemp) % Values

    Min_T = MINVAL( T_T ) * Kelvin
    Max_T = MAXVAL( T_T ) * Kelvin

    ALLOCATE( Yp_T(EOS % TS % nPoints(iYp)) )
    Yp_T = EOS % TS % States(iYp) % Values

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
    DO iRho=1,EOS % DV % nPoints(1)
      DO iTemp=1,EOS % DV % nPoints(2)
        DO iYp=1,EOS % DV % nPoints(3)
        
          ! Now add electron component
          ! Initialize temperature, DensitY, Yp, Zbar and Abar
          ElectronState % t    = T_T (iTemp)
          ElectronState % rho  = D_T (iRho)
          ElectronState % Y_e  = Yp_T(iYp)
          ElectronState % abar = 1.0_dp
          ElectronState % zbar = 1.0_dp
          ! ElectronState % abar = (Xn_T(iRho, iTemp, iYp) + Xp_T(iRho, iTemp, iYp) + &
          !                         Xa_T(iRho, iTemp, iYp) + Xh_T(iRho, iTemp, iYp)) / &
          !                        (Xn_T(iRho, iTemp, iYp) + Xp_T(iRho, iTemp, iYp) + &
          !                         Xa_T(iRho, iTemp, iYp)/4.0d0 + Xh_T(iRho, iTemp, iYp)/Ah_T(iRho, iTemp, iYp))
          ! ElectronState % zbar =  Yp_T (iYp) * ElectronState % abar
          
          ! calculate electron quantities
          CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

          Eele = ElectronState % e + mass_ele / rmu * ergmev * ElectronState % Y_e ! add back mass to internal Energy!
          Pele = ElectronState % p
          Sele = ElectronState % s

          Es(iRho,iTemp,iYp) = 10.0d0**( E_T(iRho,iTemp,iYp) ) + Eele - OS_E
          Ps(iRho,iTemp,iYp) = ( P_T(iRho,iTemp,iYp) ) + Pele - OS_P
          Ss(iRho,iTemp,iYp) = 10.0d0**( S_T(iRho,iTemp,iYp) ) + Sele - OS_S

        ENDDO
      ENDDO
    ENDDO
    
    CALL InitializeEOSComponentsInversion &
         ( D_T, T_T, Yp_T,  Es, Ps, Ss, &
         EquationOfStateTableName, EquationOfStateTableName, &
         Verbose_Option = Verbose )
             
    DEALLOCATE( Ps, Ss, Es )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( always, to: D_T, T_T, Yp_T, &
    !$OMP   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMl, UnitMp, UnitMn, &
    !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitAh, OS_P, OS_S, OS_E, &
    !$OMP   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, P_T, S_T, &
    !$OMP   E_T, Mp_T, Mn_T, Xp_T, Xn_T, Xa_T, Xh_T, Ah_T, &
    !$OMP   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( D_T, T_T, Yp_T, &
    !$ACC   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMl, UnitMp, UnitMn, &
    !$ACC   UnitXp, UnitXn, UnitXa, UnitXh, UnitAh, OS_P, OS_S, OS_E, &
    !$ACC   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, P_T, S_T, &
    !$ACC   E_T, Mp_T, Mn_T, Xp_T, Xn_T, Xa_T, Xh_T, Ah_T, &
    !$ACC   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
#endif

#endif

  END SUBROUTINE InitializeEquationOfState_TABLE


  SUBROUTINE FinalizeEquationOfState_TABLE

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: D_T, T_T, Yp_T, &
    !$OMP   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMl, UnitMp, UnitMn, &
    !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitAh, OS_P, OS_S, OS_E, &
    !$OMP   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, P_T, S_T, &
    !$OMP   E_T, Mp_T, Mn_T, Xp_T, Xn_T, Xa_T, Xh_T, &
    !$OMP   Min_D, Min_T, Min_Y, Max_D, Max_T, Max_Y )
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

    IF ( .NOT. UsingExternalEOS ) THEN
      IF (ASSOCIATED(MuonTable))      DEALLOCATE( MuonTable )
      IF (ASSOCIATED(HelmholtzTable)) DEALLOCATE( HelmholtzTable )
      IF (ASSOCIATED(EOS))            DEALLOCATE( EOS )
    END IF
  
#endif

  END SUBROUTINE FinalizeEquationOfState_TABLE


  SUBROUTINE ApplyChemicalPotentialShift_TABLE( Mp_T, Mn_T, OS_Mp, OS_Mn )

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

    REAL(dp), INTENT(inout), DIMENSION(:,:,:) :: Mp_T, Mn_T
    REAL(dp), INTENT(inout) :: OS_Mp, OS_Mn

    REAL(DP), PARAMETER :: neutron_mass = 939.56542052d0
    REAL(DP), PARAMETER :: proton_mass = 938.2720813d0
    REAL(DP), PARAMETER :: dmnp = 1.29333922d0
    REAL(DP) :: min_M, OS_M_new

    ! Apply the shift for proton chemical potential
    IF ( OS_Mp > 0.0d0 ) THEN
      min_M = -0.5d0 * OS_Mp
    ELSE
      min_M = MINVAL( 10.0d0**Mp_T )
    ENDIF
    OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + proton_mass + dmnp )
    Mp_T     = LOG10( 10.0d0**Mp_T - OS_Mp + proton_mass + dmnp + OS_M_new)
    OS_Mp    = OS_M_new

    ! Apply the shift for neutron chemical potential
    IF ( OS_Mn > 0.0d0 ) THEN
      min_M = -0.5d0 * OS_Mn
    ELSE
      min_M = MINVAL( 10.0d0**Mn_T )
    ENDIF
    OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + proton_mass + dmnp )
    Mn_T     = LOG10( 10.0d0**Mn_T - OS_Mn + proton_mass + dmnp + OS_M_new)
    OS_Mn    = OS_M_new

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

    REAL(DP) :: Pbary, Sbary, Ebary, Pele, Sele, Eele, &
                P_mu, S_mu, E_mu

    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    
#ifdef INVERSION_COMBINED
    ! --- Interpolate Pressure ----------------------------------------

    CALL ComputeDependentVariableTotal_TABLE_Scalar &
           ( D, T, Ye, Ym, P, P_T, OS_P, &
           UnitP, 0.0_dp, 1.0_dp, 0.0_dp )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariableTotal_TABLE_Scalar &
           ( D, T, Ye, Ym, S, S_T, OS_S, &
           UnitS, 0.0_dp, 0.0_dp, 1.0_dp )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariableTotal_TABLE_Scalar &
           ( D, T, Ye, Ym, E, E_T, OS_E, &
           UnitE, 1.0_dp, 0.0_dp, 0.0_dp )
#else

    ! Calculate Electron Quantities
    ! Initialize Electron state (Abar and Zbar not needed!!!)
    ElectronState % t = T
    ElectronState % rho = D
    ElectronState % y_e = Ye
    ElectronState % abar = One
    ElectronState % zbar = ElectronState % y_e
    
    CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

    Eele = ElectronState % e + mass_ele / rmu * ergmev * ElectronState % Y_e ! add back mass to internal Energy!
    Pele = ElectronState % p
    Sele = ElectronState % s
    Mue = ElectronState % mu_e

    ! Calculate Muon Quantities
    MuonState % t = T
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)

    E_mu = MuonState % e + mass_mu / rmu * ergmev * Ym ! add back mass to internal Energy!
    P_mu = MuonState % p
    S_mu = MuonState % s
    Mum = MuonState % mu

    ! --- Interpolate Pressure ----------------------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, P, P_T, OS_P, UnitP )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, S, S_T, OS_S, UnitS )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, E, E_T, OS_E, UnitE )
           
    E = E + Eele + E_mu
    P = P + Pele + P_mu
    S = S + Sele + S_mu
    
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

    D_P = D / ( Gram / Centimeter**3 )
    E_P = E / ( Erg / Gram )
    Ye_P = Ye
    Ym_P = Ym

    T_Guess = Guess / Kelvin

    CALL ComputeTemperatureWith_DEYpYl_Single_Guess_Error &
           ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, & 
              T_Lookup, T_Guess, Error )

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

    D_P = D / ( Gram / Centimeter**3 )
    E_P = E / ( Erg / Gram )
    Ye_P = Ye
    Ym_P = Ym
    
    T_Guess = Guess / Kelvin

    CALL ComputeTemperatureWith_DEYpYl_Single_Guess_NoError &
           ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
           T_Lookup, T_Guess )

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

    D_P  = D / ( Gram / Centimeter**3 )
    E_P  = E / ( Erg / Gram )
    Ye_P = Ye
    Ym_P = Ym
    
    CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error &
           ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
           T_Lookup, Error )

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

    D_P = D / ( Gram / Centimeter**3 )
    E_P = E / ( Erg / Gram )
    Ye_P = Ye
    Ym_P = Ym
    
    CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError &
           ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
           T_Lookup )

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

        D_P = D(iP) / ( Gram / Centimeter**3 )
        E_P = E(iP) / ( Erg / Gram )
        Ye_P = Ye(iP)
        Ym_P = Ym(iP)
        
        T_Guess = Guess_Option(iP) / Kelvin

        CALL ComputeTemperatureWith_DEYpYl_Single_Guess_Error &
               ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
               T_Lookup, T_Guess, Error(iP) )

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

        D_P = D(iP) / ( Gram / Centimeter**3 )
        E_P = E(iP) / ( Erg / Gram )
        Ye_P = Ye(iP)
        Ym_P = Ym(iP)
        
        CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error &
               ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, T_Lookup, &
                 Error(iP) )

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
          D_P = D(iP) / ( Gram / Centimeter**3 )
          E_P = E(iP) / ( Erg / Gram )
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

    Em  = Ev / D              ! --- Internal Energy per Mass
    Ye  = Ne / D * BaryonMass ! --- Electron Fraction
    Ym = Ne / D * BaryonMass ! --- Muon Fraction
  
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

    REAL(DP) :: D_P, E_P, Yp_P, Ye_P, Ym_P, T_P, T

#ifdef INVERSION_COMBINED
#else
    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    REAL(DP) :: Pele, P_mu
#endif

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D  / UnitD
    E_P = Em / UnitE
    Ye_P = Ye  / UnitY
    Ym_P = Ym  / UnitY
    Yp_P = Ye_P + Ym_P
    
    CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError &
           ( D_P, E_P, Ye_P, Ym_P, D_T, T_T, Yp_T, E_T, OS_E, &
           T_P )

    T = T_P * UnitT

#ifdef INVERSION_COMBINED

    CALL ComputeDependentVariableTotal_TABLE_Scalar &
           ( D, T, Ye_P, Ym_P, P, P_T, OS_P, &
           UnitP, 0.0_dp, 1.0_dp, 0.0_dp )
#else
    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye_P, Ym_P, P, P_T, OS_P, UnitP )
           
    ! Calculate Electron Quantities
    ! Initialize Electron state (Abar and Zbar not needed!!!)
    ElectronState % t = T
    ElectronState % rho = D
    ElectronState % Y_e = Ye_P
    
    CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)
    Pele = ElectronState % p

    ! Calculate Muon Quantities
    MuonState % t = T
    MuonState % rhoym = D * Ym_P

    CALL FullMuonEOS(MuonTable, MuonState)
    P_mu = MuonState % p
           
    P = P + Pele + P_mu
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

    CALL CalculateSoundSpeed( D, T, Ye, Ym, D_T, T_T, Yp_T, &
        P_T, OS_P, E_T, OS_E, HelmholtzTable, MuonTable, Gm, Cs, .true.)

    !Cs = SQRT( Gm * P / D )

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

      CALL CalculateSoundSpeed( D(iP), T(iP), Ye(iP), Ym(iP), &
          D_T, T_T, Yp_T, P_T, OS_P, E_T, OS_E, &
          HelmholtzTable, MuonTable, Gm(iP), Cs(iP),.true.)

    ENDDO

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

    CALL ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error &
           ( D_P, P_P, Ye_P, Ym_P, D_T, T_T, Yp_T, P_T, OS_P, T_P, &
             Error )

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
    
#ifdef INVERSION_COMBINED
#else
    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    REAL(DP) :: Eele, E_mu
#endif
    ! --- Interpolate Specific Internal Energy ------------------------

#ifdef INVERSION_COMBINED

    CALL ComputeDependentVariableTotal_TABLE_Scalar &
           ( D, T, Ye, Ym, Em, E_T, OS_E, &
           UnitE, 1.0_dp, 0.0_dp, 0.0_dp )
#else

    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, Em, E_T, OS_E, UnitE )
           
    ! Calculate Electron Quantities
    ! Initialize Electron state (Abar and Zbar not needed!!!)
    ElectronState % t = T
    ElectronState % rho = D
    ElectronState % Y_e = Ye
    ElectronState % abar = One 
    ElectronState % zbar = Ye
    
    CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)
    Eele = ElectronState % e + mass_ele / rmu * ergmev * ElectronState % Y_e ! add back mass to internal Energy!

    ! Calculate Muon Quantities
    MuonState % t = T
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
    E_mu = MuonState % e + mass_mu / rmu * ergmev * Ym ! add back mass to internal Energy!
           
    Em = Em + Eele + E_mu
#endif

    Ev  = Em * D                ! --- Internal Energy per Unit Volume
    Ne  = Ye  * D / BaryonMass  ! --- Electrons per Unit Volume
    Nm = Ym  * D / BaryonMass ! --- Muons per Unit Volume

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Scalar

  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Vector &
    ( D, T, Ye, Ym, Ev, Em, Ne, Nm )

    REAL(DP), INTENT(in)  :: D (1:), T (1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: Ev(1:), Em(1:), Ne(1:), Nm(1:)

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

  ! ! Calculation of sound speed now is more complicated, maybe this subroutine is not needed
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

    CALL CalculateSoundSpeed( D, T, Ye, Ym, D_T, T_T, Yp_T, &
        P_T, OS_P, E_T, OS_E, HelmholtzTable, MuonTable, Gm, Cs, .TRUE.)

    !Cs = SQRT( Gm * P / D )

  END SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Scalar

  ! Calculation of sound speed now is more complicated, maybe this subroutine is not needed
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

#ifdef INVERSION_COMBINED
#else
    TYPE(ElectronStateType) :: ElectronState
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
      
      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
             ( D, T, Ye, Ym, P, dPdD, dPdT, dPdYe, dPdYm, P_T, OS_P, &
             UnitP, 1.0_dp, 0.0_dp, 0.0_dp )

    ELSE

#ifdef INVERSION_COMBINED
      CALL ComputeDependentVariableTotal_TABLE_Scalar &
             ( D, T, Ye, Ym, P, P_T, OS_P, &
             UnitP, 1.0_dp, 0.0_dp, 0.0_dp )
#else
      CALL ComputeDependentVariableBaryons_TABLE_Scalar &
             ( D, T, Ye, Ym, P, P_T, OS_P, UnitP )

      ! Calculate Electron Quantities
      ! Initialize Electron state (Abar and Zbar not needed!!!)
      ElectronState % t = T
      ElectronState % rho = D
      ElectronState % Y_e = Ye
      
      CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)
      Pele = ElectronState % p

      ! Calculate Muon Quantities
      MuonState % t = T
      MuonState % rhoym = D * Ym
      
      CALL FullMuonEOS(MuonTable, MuonState)
      P_mu = MuonState % p
      
      P = P + Pele + P_mu
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
             UnitP, 1.0_dp, 0.0_dp, 0.0_dp )

    ELSE

      CALL ComputeDependentVariableTotal_TABLE_Vector &
             ( D, T, Ye, Ym, P, P_T, OS_P, &
             UnitP, 1.0_dp, 0.0_dp, 0.0_dp )

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

#ifdef INVERSION_COMBINED
#else
    TYPE(ElectronStateType) :: ElectronState
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

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
             ( D, T, Ye, Ym, E, dEdD, dEdT, dEdYe, dEdYm, E_T, OS_E, &
             UnitE, 1.0_dp, 0.0_dp, 0.0_dp )

    ELSE

#ifdef INVERSION_COMBINED

    CALL ComputeDependentVariableTotal_TABLE_Scalar &
           ( D, T, Ye, Ym, E, E_T, OS_E, &
           UnitE, 1.0_dp, 0.0_dp, 0.0_dp )
#else
    
    CALL ComputeDependentVariableBaryons_TABLE_Scalar &
           ( D, T, Ye, Ym, E, E_T, OS_E, UnitE )
           
    ! Calculate Electron Quantities
    ! Initialize Electron state (Abar and Zbar not needed!!!)
    ElectronState % t = T
    ElectronState % rho = D
    ElectronState % Y_e = Ye
    
    CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)
    Eele = ElectronState % e + mass_ele / rmu * ergmev * ElectronState % Y_e ! add back mass to internal Energy!

    ! Calculate Muon Quantities
    MuonState % t = T
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
    E_mu = MuonState % e + mass_mu / rmu * ergmev * Ym ! add back mass to internal Energy!
           
    E = E + Eele + E_mu
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
             UnitE, 1.0_dp, 0.0_dp, 0.0_dp )

    ELSE

#ifdef INVERSION_COMBINED

      CALL ComputeDependentVariableTotal_TABLE_Vector &
             ( D, T, Ye, Ym, E, E_T, OS_E, &
             UnitE, 1.0_dp, 0.0_dp, 0.0_dp )

#else
      ! Not really used
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
    
    TYPE(ElectronStateType) :: ElectronState
    
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
      
    ELSE
      
      ! Calculate Electron Quantities
      ! Initialize Electron state (Abar and Zbar not needed!!!)
      ElectronState % t = T
      ElectronState % rho = D
      ElectronState % Y_e = Ye
      
      CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

      M = ElectronState % mu_e
      
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
    INTEGER  :: iD, iT

    TYPE(MuonStateType) :: MuonState
    
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
      
      CALL GetIndexAndDelta_Log( D * Ym, MuonTable % rhoym(:), iD, dD )
      CALL GetIndexAndDelta_Log( T, MuonTable % t(:), iT, dT )
      
      aD = 1.0_dp / ( D * LOG10( MuonTable % rhoym(iD+1) / MuonTable % rhoym(iD) ) )
      aT = 1.0_dp / ( T * LOG10( MuonTable % t(iT+1) / MuonTable % t(iT) ) )

      CALL LinearInterpDeriv_Array_Point &
             ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % mu), M, &
               dMdD, dMdT )
      
      dMdYm = dMdD * D ! make sure the derivative is wr2 ym, not rhoym
      dMdD = dMdD * Ym ! make sure the derivative is wr2 rho, not rhoym
      dMdYe = 0.0_dp
      
    ELSE

      MuonState % t = T
      MuonState % rhoym = D * Ym
      
      CALL FullMuonEOS(MuonTable, MuonState)

      M = MuonState % mu

    ENDIF
    
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

      DO iP=1,nP
          CALL ComputeMuonChemicalPotential_TABLE_Scalar &
                (D(iP), T(iP), Ye(iP), Ym(iP), M(iP), dMdD(iP), dMdT(iP), dMdYe(iP), dMdYm(iP))
      END DO

    ELSE 
    
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

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P
    REAL(DP) :: Mue, Mup, Mun

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P
    
    CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P, Ym_P, Mun )
      
    CALL ComputeProtonChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P, Ym_P, Mup )
      
    CALL ComputeElectronChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P, Ym_P, Mue )

    Mue  = Mue  * UnitMl
    Mup  = Mup  * UnitMp
    Mun  = Mun  * UnitMn

    Munue  = ( Mue  + Mup ) - Mun

#else

    Munue  = Zero

#endif

  END SUBROUTINE ComputeElectronNeutrinoChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeElectronNeutrinoChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, Munue )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: Munue(1:)

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P
    REAL(DP) :: Mue, Mup, Mun

    INTEGER  :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, iD, iT, iY, dD, dT, dY, Mue, Mup, Mun ) &
    !$OMP MAP( to: D, T, Y, D_T, T_T, Yp_T, OS_Mp, OS_Mn, Mp_T, Mn_T ) &
    !$OMP MAP( from: Munue )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, iD, iT, iY, dD, dT, dY, Mue, Mup, Mun ) &
    !$ACC COPYIN( D, T, Y, D_T, T_T, Yp_T, OS_Mp, OS_Mn, Mp_T, Mn_T ) &
    !$ACC COPYOUT( Munue )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, iD, iT, iY, dD, dT, dY, Mue, Mup, Mun )
#endif
    DO iP = 1, nP

      D_P = D(iP) / UnitD
      T_P = T(iP) / UnitT
      Ye_P = Ye(iP) / UnitY
      Ym_P = Ym(iP) / UnitY
      Yp_P = Ye_P + Ym_P

      CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P, Ym_P, Mun )
        
      CALL ComputeProtonChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P, Ym_P, Mup )
        
      CALL ComputeElectronChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P, Ym_P, Mue )

      Mue  = Mue  * UnitMl
      Mup  = Mup  * UnitMp
      Mun  = Mun  * UnitMn

      Munue(iP)  = ( Mue  + Mup ) - Mun

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

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P
    REAL(DP) :: Mum, Mup, Mun

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P
    
    CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P, Ym_P, Mun )
      
    CALL ComputeProtonChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P, Ym_P, Mup )
      
    CALL ComputeMuonChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P, Ym_P, Munum )

    Mum  = Mum  * UnitMl
    Mup  = Mup  * UnitMp
    Mun  = Mun  * UnitMn

    Munum  = ( Mum  + Mup ) - Mun

#else

    Munum  = Zero

#endif

  END SUBROUTINE ComputeMuonNeutrinoChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeMuonNeutrinoChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, Munum )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: Munum(1:)

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P
    REAL(DP) :: Mup, Mun, Mum

    INTEGER  :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, iD, iT, iY, dD, dT, dY, Mum, Mup, Mun ) &
    !$OMP MAP( to: D, T, Y, D_T, T_T, Yp_T, OS_Mp, OS_Mn, Mp_T, Mn_T ) &
    !$OMP MAP( from: Munum )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, iD, iT, iY, dD, dT, dY, Mum, Mup, Mun ) &
    !$ACC COPYIN( D, T, Y, D_T, T_T, Yp_T, OS_Mp, OS_Mn, Mp_T, Mn_T ) &
    !$ACC COPYOUT( Munum )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, iD, iT, iY, dD, dT, dY, Mum, Mup, Mun )
#endif
    DO iP = 1, nP

      D_P = D(iP) / UnitD
      T_P = T(iP) / UnitT
      Ye_P = Ye(iP) / UnitY
      Ym_P = Ym(iP) / UnitY
      Yp_P = Ye_P + Ym_P

      CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P, Ym_P, Mun )
        
      CALL ComputeProtonChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P, Ym_P, Mup )
        
      CALL ComputeMuonChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P, Ym_P, Mum )

      Mum  = Mum  * UnitMl
      Mup  = Mup  * UnitMp
      Mun  = Mun  * UnitMn

      Munum(iP)  = ( Mum  + Mup ) - Mun

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

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, D_T, T_T, Yp_T, OS_V, V_T, V_P )

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
      !$OMP MAP( to: D, T, Y, D_T, T_T, Yp_T, OS_V, V_T ) &
      !$OMP MAP( from: V )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P ) &
      !$ACC COPYIN( D, T, Y, D_T, T_T, Yp_T, OS_V, V_T ) &
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

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T, T_T, Yp_T, OS_V, V_T, V_P )

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

  END SUBROUTINE ComputeDependentVariableBaryons_TABLE_Vector

  SUBROUTINE ComputeDependentVariableTotal_TABLE_Scalar &
    ( D, T, Ye, Ym, V, V_T, OS_V, Units_V, &
    InputE, InputP, InputS )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: V
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    REAL(DP), INTENT(in)  :: InputE, InputP, InputS

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: Vtot(2,2,2), E_leptons(2,2,2), &
                P_leptons(2,2,2), S_leptons(2,2,2)
    
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps

    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P

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

    DO iL_T=1,2
      DO iL_D=1,2
        DO iL_Y=1,2
          ElectronState % t = T_T(iT+iL_T-1)
          ElectronState % rho = D_T(iD+iL_D-1)
          ElectronState % y_e = Yp_T(iYp+iL_Y-1) * Ye_over_Yp
          ElectronState % abar = One
          ElectronState % zbar = ElectronState % y_e
          
          ! CALL FullHelmEOS(1, HelmholtzTable, ElectronState, .false., .false.)
          CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

          MuonState % t = T_T(iT+iL_T-1)
          MuonState % rhoym = D_T(iD+iL_D-1) * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
          
          CALL FullMuonEOS(MuonTable, MuonState)
          
          E_leptons(iL_D,iL_T,iL_Y) = ElectronState % e + mass_ele / rmu * ergmev * ElectronState % y_e + &
                        MuonState % e + mass_mu / rmu * ergmev * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
          P_leptons(iL_D,iL_T,iL_Y) = ElectronState % p + MuonState % p
          S_leptons(iL_D,iL_T,iL_Y) = ElectronState % s + MuonState % s

        END DO
      END DO
    END DO
    
    IF (InputP .EQ. 1.0d0) THEN
      Vtot(:,:,:) = V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + P_leptons(:,:,:)
    ELSE IF (InputE .EQ. 1.0d0) THEN
      Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_leptons(:,:,:)
    ELSE
      Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_leptons(:,:,:)
    ENDIF
    
    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
           OS_V, Vtot, V_P )

    V = V_P * Units_V

#else

    V = Zero

#endif

  END SUBROUTINE ComputeDependentVariableTotal_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableTotal_TABLE_Vector &
    ( D, T, Ye, Ym, V, V_T, OS_V, Units_V, &
    InputE, InputP, InputS )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: V(1:)
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    REAL(DP), INTENT(in)  :: InputE, InputP, InputS

    INTEGER  :: iP, nP
    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: Vtot(2,2,2), E_leptons(2,2,2), &
                P_leptons(2,2,2), S_leptons(2,2,2)
    
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps

    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
        
    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P ) &
      !$OMP MAP( to: D, T, Y, D_T, T_T, Yp_T, OS_V, V_T ) &
      !$OMP MAP( from: V )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( D_P, T_P, Yp_P, Ye_P, Ym_P, V_P ) &
      !$ACC COPYIN( D, T, Y, D_T, T_T, Yp_T, OS_V, V_T ) &
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
            ElectronState % t = T_T(iT+iL_T-1)
            ElectronState % rho = D_T(iD+iL_D-1)
            ElectronState % y_e = Yp_T(iYp+iL_Y-1) * Ye_over_Yp
            ElectronState % abar = One
            ElectronState % zbar = ElectronState % y_e
            
            ! CALL FullHelmEOS(1, HelmholtzTable, ElectronState, .false., .false.)
            CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

            MuonState % t = T_T(iT+iL_T-1)
            MuonState % rhoym = D_T(iD+iL_D-1) * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
            
            CALL FullMuonEOS(MuonTable, MuonState)
            
            E_leptons(iL_D,iL_T,iL_Y) = ElectronState % e + mass_ele / rmu * ergmev * ElectronState % y_e + &
                          MuonState % e + mass_mu / rmu * ergmev * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
            P_leptons(iL_D,iL_T,iL_Y) = ElectronState % p + MuonState % p
            S_leptons(iL_D,iL_T,iL_Y) = ElectronState % s + MuonState % s
            
          END DO
        END DO
      END DO

      IF (InputP .EQ. 1.0d0) THEN
        Vtot(:,:,:) = V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + P_leptons(:,:,:)
      ELSE IF (InputE .EQ. 1.0d0) THEN
        Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_leptons(:,:,:)
      ELSE
        Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_leptons(:,:,:)
      ENDIF

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_V, Vtot, V_P )

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


  SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
    ( D, T, Ye, Ym, V, dVdD, dVdT, dVdYe, dVdYm, V_T, OS_V, Units_V, &
    InputE, InputP, InputS)

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: V, dVdD, dVdT, dVdYe, dVdYm
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    REAL(DP), INTENT(in)  :: InputE, InputP, InputS

    REAL(DP) :: D_P, T_P, Yp_P, Ye_P, Ym_P, V_P, dV_P(3)
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: Vtot(2,2,2), E_leptons(2,2,2), &
                P_leptons(2,2,2), S_leptons(2,2,2)
    
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps
    REAL(DP) :: dD, dT
	  REAL(DP) :: aD, aT
    REAL(DP) :: dVdrhoym, dVdummy, Vdummy, LocalOffset
    
    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P
    
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
          ElectronState % t = T_T(iT+iL_T-1)
          ElectronState % rho = D_T(iD+iL_D-1)
          ElectronState % y_e = Yp_T(iYp+iL_Y-1) * Ye_over_Yp
          ElectronState % abar = One
          ElectronState % zbar = ElectronState % y_e
          
          ! CALL FullHelmEOS(1, HelmholtzTable, ElectronState, .false., .false.)
          CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

          MuonState % t = T_T(iT+iL_T-1)
          MuonState % rhoym = D_T(iD+iL_D-1) * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
          
          CALL FullMuonEOS(MuonTable, MuonState)
          
          E_leptons(iL_D,iL_T,iL_Y) = ElectronState % e + mass_ele / rmu * ergmev * ElectronState % y_e + &
                        MuonState % e + mass_mu / rmu * ergmev * Yp_T(iYp+iL_Y-1) * Ym_over_Yp
          P_leptons(iL_D,iL_T,iL_Y) = ElectronState % p + MuonState % p
          S_leptons(iL_D,iL_T,iL_Y) = ElectronState % s + MuonState % s

        END DO
      END DO
    END DO

    ! Now you have to calculate dVbarydYp, and if you have pressure, you
    ! have to handle the offset accordingly
    
    IF (InputP .EQ. 1.0d0) THEN
      Vtot(:,:,:) = V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + P_leptons(:,:,:)
   
      LocalOffset = MINVAL(V_T(iD:iD+1,iT:iT+1,iYp:iYp+1))
      IF (LocalOffset .lt. 0.0_dp) THEN
          LocalOffset = -1.1d0*LocalOffset
      ELSE
          LocalOffset = 0.0_dp
      ENDIF
      
      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             LocalOffset, LOG10(V_T + LocalOffset), V_P, dV_P )
             
      ! Now calculate electron and muon derivatives
      ! ELECTRON PART IS EASY -----------------------------------------------------
      ! Initialize temperature, density, yp, Zbar and Abar
      ElectronState % t = T
      ElectronState % rho = D
      ElectronState % abar = 1.0d0 ! these are only used for ion contribution
      ElectronState % zbar = Ye ! these are only used for ion contribution
      ElectronState % y_e = Ye

      ! calculate electron quantities
      CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

      dVdYe = D/Ye * ElectronState % dpdr ! This is valid because dPdrho for photons is zero
      
      ! NOW MUONS ---------------------------------------------
      IF (( D * Ym .lt. MuonTable % rhoym(1) ) .or. (T .lt. MuonTable % t(1))) THEN
        
        dVdYm = 0.0_dp
      
      ELSE
      
        CALL GetIndexAndDelta_Log( D * Ym, MuonTable % rhoym(:), iD, dD )
        CALL GetIndexAndDelta_Log( T, MuonTable % t(:), iT, dT )
        
        aD = 1.0_dp / ( D * LOG10( MuonTable % rhoym(iD+1) / MuonTable % rhoym(iD) ) )
        aT = 1.0_dp / ( T * LOG10( MuonTable % t(iT+1) / MuonTable % t(iT) ) )

        CALL LinearInterpDeriv_Array_Point &
               ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % p), Vdummy, &
                 dVdrhoym, dVdummy )
        
        dVdYm = dVdrhoym * D ! make sure the derivative is wr2 ym, not rhoym
      END IF
    
    ELSE IF (InputE .EQ. 1.0d0) THEN

      Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + E_leptons(:,:,:)

      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_V, V_T, V_P, dV_P )
             
      ! Now calculate electron and muon derivatives
      ! ELECTRON PART IS EASY -----------------------------------------------------
      ! Initialize temperature, density, yp, Zbar and Abar
      ElectronState % t = T
      ElectronState % rho = D
      ElectronState % abar = 1.0d0 ! these are only used for ion contribution
      ElectronState % zbar = 1.0d0 ! these are only used for ion contribution
      ElectronState % y_e = Ye

      ! calculate electron quantities
      CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

      dVdYe = D/Ye * ElectronState % dedr ! This is valid because dPdrho for photons is zero
      WRITE(*,*) 'Need to get this from Electron EOS, careful'
      STOP
      
      ! NOW MUONS ---------------------------------------------
      IF (( D * Ym .lt. MuonTable % rhoym(1) ) .or. (T .lt. MuonTable % t(1))) THEN
        
        dVdYm = 0.0_dp
      
      ELSE
      
        CALL GetIndexAndDelta_Log( D * Ym, MuonTable % rhoym(:), iD, dD )
        CALL GetIndexAndDelta_Log( T, MuonTable % t(:), iT, dT )
        
        aD = 1.0_dp / ( D * LOG10( MuonTable % rhoym(iD+1) / MuonTable % rhoym(iD) ) )
        aT = 1.0_dp / ( T * LOG10( MuonTable % t(iT+1) / MuonTable % t(iT) ) )

        CALL LinearInterpDeriv_Array_Point &
               ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % e), Vdummy, &
                 dVdrhoym, dVdummy )
        
        dVdYm = dVdrhoym * D ! make sure the derivative is wr2 ym, not rhoym
      END IF
             
    ELSE

      Vtot(:,:,:) = 10.0**V_T(iD:iD+1,iT:iT+1,iYp:iYp+1) + S_leptons(:,:,:)

      ! V_P is a dummy variable below
      CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
             OS_V, V_T, V_P, dV_P )
             
      ! Now calculate electron and muon derivatives
      ! ELECTRON PART IS EASY -----------------------------------------------------
      ! Initialize temperature, density, yp, Zbar and Abar
      ElectronState % t = T
      ElectronState % rho = D
      ElectronState % abar = 1.0d0 ! these are only used for ion contribution
      ElectronState % zbar = 1.0d0 ! these are only used for ion contribution
      ElectronState % y_e = Ye

      ! calculate electron quantities
      CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

      dVdYe = D/Ye * ElectronState % dsdr ! This is valid because dPdrho for photons is zero
      WRITE(*,*) 'Need to get this properly from Electron EOS, careful'
      STOP
      
      ! NOW MUONS ---------------------------------------------
      IF (( D * Ym .lt. MuonTable % rhoym(1) ) .or. (T .lt. MuonTable % t(1))) THEN
        
        dVdYm = 0.0_dp
      
      ELSE
      
        CALL GetIndexAndDelta_Log( D * Ym, MuonTable % rhoym(:), iD, dD )
        CALL GetIndexAndDelta_Log( T, MuonTable % t(:), iT, dT )
        
        aD = 1.0_dp / ( D * LOG10( MuonTable % rhoym(iD+1) / MuonTable % rhoym(iD) ) )
        aT = 1.0_dp / ( T * LOG10( MuonTable % t(iT+1) / MuonTable % t(iT) ) )

        CALL LinearInterpDeriv_Array_Point &
               ( iD, iT, dD, dT, aD, aT, 0.0_dp, LOG10(MuonTable % s), Vdummy, &
                 dVdrhoym, dVdummy )
        
        dVdYm = dVdrhoym * D ! make sure the derivative is wr2 ym, not rhoym
      END IF
             
    END IF
    
    CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
       ( D_P, T_P, Yp_P, D_T(iD:iD+1), T_T(iT:iT+1), Yp_T(iYp:iYp+1), &
       OS_V, Vtot, V_P, dV_P )
             
    V = V_P * Units_V

    dVdD = dV_P(1) * Units_V / UnitD
    dVdT = dV_P(2) * Units_V / UnitT

#else

    V    = Zero
    dVdD = Zero
    dVdT = Zero
    dVdYp = Zero

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Vector &
    ( D, T, Ye, Ym, V, dVdD, dVdT, dVdYe, dVdYm, V_T, OS_V, Units_V, &
    InputE, InputP, InputS )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: V(1:), dVdD(1:), dVdT(1:), dVdYe(1:), dVdYm(1:)
    REAL(DP), INTENT(in)  :: V_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    REAL(DP), INTENT(in)  :: InputE, InputP, InputS

    INTEGER :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP MAP( to: D, T, Y, OS_V, V_T ) &
      !$OMP MAP( from: V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Y, OS_V, V_T ) &
      !$ACC COPYOUT( V, dVdD, dVdT, dVdYe, dVdYm )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
             ( D(iP), T(iP), Ye(iP), Ym(iP), V(iP), dVdD(iP), dVdT(iP), &
               dVdYe(iP), dVdYm(iP), V_T, OS_V, Units_V, InputE, InputP, InputS )

    END DO

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: V, dVdD, dVdT, dVdYp )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( V, dVdD, dVdT, dVdYp )
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

    INTEGER  :: iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY
    Yp_P = Ye_P + Ym_P

    CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, D_T, T_T, Yp_T, OS_V, V_T, V_P, dV_P )
             
    V = V_P * Units_V

    dVdD = dV_P(1) * Units_V / UnitD
    dVdT = dV_P(2) * Units_V / UnitT
    dVdYe = dV_P(3) * Units_V / UnitY
    dVdYm = dV_P(3) * Units_V / UnitY

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
      !$OMP MAP( to: D, T, Y, OS_V, V_T ) &
      !$OMP MAP( from: V, dVdD, dVdT, dVdYp )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Y, OS_V, V_T ) &
      !$ACC COPYOUT( V, dVdD, dVdT, dVdYp )
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
    !$OMP MAP( from: V, dVdD, dVdT, dVdYp )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( V, dVdD, dVdT, dVdYp )
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
