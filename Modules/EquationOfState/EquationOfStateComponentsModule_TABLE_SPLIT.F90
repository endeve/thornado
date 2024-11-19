#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_EOS
#endif

MODULE EquationOfStateComponentsModule_TABLE

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
  USE wlwlSoundSpeedModule, ONLY: &
    CalculatewlSoundSpeedModule
  USE wlEosConstantsModule, ONLY: &
    kmev, rmu, kmev_inv, ergmev, me, mmu, cvel
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, Index1D_Log
    
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
    iXp_T, iXn_T, iXa_T, iXh_T, iGm_T
  REAL(DP) :: &
    UnitD, UnitT, UnitY, &
    UnitP, UnitS, UnitE, UnitMl, UnitMp, UnitMn, &
    UnitXp, UnitXn, UnitXa, UnitXh, UnitGm
  REAL(DP), PUBLIC :: &
    OS_P, OS_S, OS_E, OS_Mp, OS_Mn, &
    OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Mue_placeholder
  REAL(DP), PARAMETER :: &
    BaryonMass = AtomicMassUnit
  REAL(DP) :: minvar, OS_loc
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Dbary_T, Tbary_T, Ypbary_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    Pbary_T, Sbary_T, Ebary_T, Mpbary_T, Mnbary_T, &
    Xpbary_T, Xnbary_T, Xabary_T, Xhbary_T, Mue_placeholder
#ifdef MICROPHYSICS_WEAKLIB
  LOGICAL :: UsingExternalEOS
  TYPE(EquationOfStateTableType), POINTER :: EOSBary
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
  PUBLIC :: ComputeProtonChemicalPotential_TABLE
  PUBLIC :: ComputeNeutronChemicalPotential_TABLE
  PUBLIC :: ComputeProtonMassFraction_TABLE
  PUBLIC :: ComputeNeutronMassFraction_TABLE
  PUBLIC :: ComputeHeavyMassFraction_TABLE
  PUBLIC :: ComputeHeavyMassNumber_TABLE
  PUBLIC :: ComputeNeutrinoChemicalPotential_TABLE

  REAL(DP), PUBLIC :: Min_D, Min_T, Min_Yp
  REAL(DP), PUBLIC :: Max_D, Max_T, Max_Yp

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

  INTERFACE ComputeHeavyMassFraction_TABLE
      MODULE PROCEDURE ComputeHeavyMassFraction_TABLE_Scalar
    MODULE PROCEDURE ComputeHeavyMassFraction_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeHeavyMassNumber_TABLE
    MODULE PROCEDURE ComputeHeavyMassNumber_TABLE_Scalar
    MODULE PROCEDURE ComputeHeavyMassNumber_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeNeutrinoChemicalPotential_TABLE
    MODULE PROCEDURE ComputeNeutrinoChemicalPotential_TABLE_Scalar
    MODULE PROCEDURE ComputeNeutrinoChemicalPotential_TABLE_Vector
  END INTERFACE

  INTERFACE ComputeDependentVariableTotal
    MODULE PROCEDURE ComputeDependentVariableTotal_Scalar
    MODULE PROCEDURE ComputeDependentVariableTotal_Vector
  END INTERFACE ComputeDependentVariableTotal

  INTERFACE ComputeDependentVariableAndDerivatives_TABLE
    MODULE PROCEDURE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar
    MODULE PROCEDURE ComputeDependentVariableAndDerivativesTotal_TABLE_Vector
  END INTERFACE ComputeDependentVariableAndDerivatives_TABLE

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( Dbary_T, Tbary_T, Ypbary_T, &
  !$OMP   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, UnitMl, UnitMp, UnitMn, &
  !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, &
  !$OMP   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$OMP   OS_Xa, OS_Xh, &
  !$OMP   Pbary_T, Sbary_T, Ebary_T, Mpbary_T, Mnbary_T, Xpbary_T, Xnbary_T, &
  !$OMP   Xabary_T, Xhbary_T, &
  !$OMP   Min_D, Min_T, Min_Yp, Max_D, Max_T, Max_Yp )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( Dbary_T, Tbary_T, Ypbary_T, &
  !$ACC   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, UnitMl, UnitMp, UnitMn, &
  !$ACC   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, &
  !$ACC   OS_P, OS_S, OS_E, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$ACC   OS_Xa, OS_Xh, &
  !$ACC   Pbary_T, Sbary_T, Ebary_T, Mpbary_T, Mnbary_T, Xpbary_T, Xnbary_T, &
  !$ACC   Xabary_T, Xhbary_T, &
  !$ACC   Min_D, Min_T, Min_Yp, Max_D, Max_T, Max_Yp )
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

       ALLOCATE( EOSBary )
       UsingExternalEOS = .FALSE.

       CALL InitializeHDF( )
        
       ! read in helmholtz table
       CALL ReadHelmholtzTableHDF( HelmholtzTable, TRIM( EquationOfStateTableName ) )

       ! read in helmholtz table
       CALL ReadMuonTableHDF( MuonTable, TRIM( EquationOfStateTableName ) )

       ! read in baryon table -------------------------------
       CALL ReadEquationOfStateTableHDF( EOSBary, TRIM( EquationOfStateTableName ) )
        
       CALL FinalizeHDF( )

    ELSE

       EOSBary => External_EOS
       UsingExternalEOS = .TRUE.

    END IF

    ! --- Thermodynamic State Indices ---

    iRho = EOSBary % TS % Indices % iRho
    iTemp = EOSBary % TS % Indices % iT
    iYp = EOSBary % TS % Indices % iYe

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
    UnitGm = One

    ! --- Thermodynamic States ---

    ALLOCATE( Dbary_T(EOSBary % TS % nPoints(iRho)) )
    Dbary_T = EOSBary % TS % States(iRho) % Values

    Min_D = MINVAL( Dbary_T ) * Gram / Centimeter**3
    Max_D = MAXVAL( Dbary_T ) * Gram / Centimeter**3

    ALLOCATE( Tbary_T(EOSBary % TS % nPoints(iTemp)) )
    Tbary_T = EOSBary % TS % States(iTemp) % Values

    Min_T = MINVAL( Tbary_T ) * Kelvin
    Max_T = MAXVAL( Tbary_T ) * Kelvin

    ALLOCATE( Ypbary_T(EOSBary % TS % nPoints(iYp)) )
    Ypbary_T = EOSBary % TS % States(iYp) % Values

    Min_Yp = MINVAL( Ypbary_T )
    Max_Yp = MAXVAL( Ypbary_T )

    ! --- Dependent Variables Indices ---

    iP_T  = EOSBary % DV % Indices % iPressure
    iS_T  = EOSBary % DV % Indices % iEntropyPerBaryon
    iE_T  = EOSBary % DV % Indices % iInternalEnergyDensity
    iMe_T = EOSBary % DV % Indices % iElectronChemicalPotential
    iMp_T = EOSBary % DV % Indices % iProtonChemicalPotential
    iMn_T = EOSBary % DV % Indices % iNeutronChemicalPotential
    iXp_T = EOSBary % DV % Indices % iProtonMassFraction
    iXn_T = EOSBary % DV % Indices % iNeutronMassFraction
    iXa_T = EOSBary % DV % Indices % iAlphaMassFraction
    iXh_T = EOSBary % DV % Indices % iHeavyMassFraction

    ! --- Dependent Variables Offsets ---

    OS_P  = EOSBary % DV % Offsets(iP_T)
    OS_S  = EOSBary % DV % Offsets(iS_T)
    OS_E  = EOSBary % DV % Offsets(iE_T)
    OS_Mp = EOSBary % DV % Offsets(iMp_T)
    OS_Mn = EOSBary % DV % Offsets(iMn_T)
    OS_Xp = EOSBary % DV % Offsets(iXp_T)
    OS_Xn = EOSBary % DV % Offsets(iXn_T)
    OS_Xa = EOSBary % DV % Offsets(iXa_T)
    OS_Xh = EOSBary % DV % Offsets(iXh_T)

    ALLOCATE &
      ( Pbary_T (1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Sbary_T (1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Ebary_T (1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Mpbary_T(1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Mnbary_T(1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Xpbary_T(1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Xnbary_T(1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Xabary_T(1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Xhbary_T(1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )

    Pbary_T  = EOSBary % DV % Variables(iP_T ) % Values
    Sbary_T  = EOSBary % DV % Variables(iS_T ) % Values
    Ebary_T  = EOSBary % DV % Variables(iE_T ) % Values
    Mpbary_T = EOSBary % DV % Variables(iMp_T) % Values
    Mnbary_T = EOSBary % DV % Variables(iMn_T) % Values
    Xpbary_T = EOSBary % DV % Variables(iXp_T) % Values
    Xnbary_T = EOSBary % DV % Variables(iXn_T) % Values
    Xabary_T = EOSBary % DV % Variables(iXa_T) % Values
    Xhbary_T = EOSBary % DV % Variables(iXh_T) % Values

    IF ( UseChemicalPotentialShift ) CALL ApplyChemicalPotentialShift_TABLE( Mpbary_T, Mnbary_T, OS_Mp, OS_Mn )

    ALLOCATE &
      ( Ps  (1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Ss  (1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
    ALLOCATE &
      ( Es  (1:EOSBary % DV % nPoints(1), &
             1:EOSBary % DV % nPoints(2), &
             1:EOSBary % DV % nPoints(3)) )
             
    ! Build full EOS without muons to determine the bounds of the EOS.
    ! Muons should not be so important to cause P, S, and E to go above the 
    ! bounds already calculated assuming Ym = 0.
    DO iRho=1,EOSBary % DV % nPoints(1)
      DO iTemp=1,EOSBary % DV % nPoints(2)
        DO iYp=1,EOSBary % DV % nPoints(3)
        
          ! Now add electron component
          ! Initialize temperature, DensitY, Yp, Zbar and Abar
          ElectronState % t = Tbary_T(iTemp)
          ElectronState % rho = Dbary_T(iRho)
          ElectronState % Y_e = Ypbary_T(iYp)
          
          ! calculate electron quantities
          CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

          Eele = ElectronState % e + me / rmu * ergmev * ElectronState % Y_e ! add back mass to internal Energy!
          Pele = ElectronState % p
          Sele = ElectronState % s

          Es(iRho,iTemp,iYp) = 10.0d0**( Ebary_T(iRho,iTemp,iYp) ) + Eele - OS_E
          Ps(iRho,iTemp,iYp) = 10.0d0**( Pbary_T(iRho,iTemp,iYp) ) + Pele - OS_P
          Ss(iRho,iTemp,iYp) = 10.0d0**( Sbary_T(iRho,iTemp,iYp) ) + Sele - OS_S

        ENDDO
      ENDDO
    ENDDO
    
    CALL InitializeEOSComponentsInversion &
         ( Dbary_T, Tbary_T, Ypbary_T,  Es, Ps, Ss, &
         EquationOfStateTableName, EquationOfStateTableName, &
         Verbose_Option = Verbose )
             
    DEALLOCATE( Ps, Ss, Es )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( always, to: Dbary_T, Tbary_T, Ypbary_T, &
    !$OMP   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMl, UnitMp, UnitMn, &
    !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, OS_P, OS_S, OS_E, &
    !$OMP   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, Pbary_T, Sbary_T, &
    !$OMP   Ebary_T, Mpbary_T, Mnbary_T, Xpbary_T, Xnbary_T, Xabary_T, Xhbary_T, &
    !$OMP   Min_D, Min_T, Min_Yp, Max_D, Max_T, Max_Yp )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( Dbary_T, Tbary_T, Ypbary_T, &
    !$ACC   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMl, UnitMp, UnitMn, &
    !$ACC   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, OS_P, OS_S, OS_E, &
    !$ACC   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, Pbary_T, Sbary_T, &
    !$ACC   Ebary_T, Mpbary_T, Mnbary_T, Xpbary_T, Xnbary_T, Xabary_T, Xhbary_T, &
    !$ACC   Min_D, Min_T, Min_Yp, Max_D, Max_T, Max_Yp )
#endif

#endif

  END SUBROUTINE InitializeEquationOfState_TABLE


  SUBROUTINE FinalizeEquationOfState_TABLE

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Dbary_T, Tbary_T, Ypbary_T, &
    !$OMP   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMl, UnitMp, UnitMn, &
    !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, OS_P, OS_S, OS_E, &
    !$OMP   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, Pbary_T, Sbary_T, &
    !$OMP   Ebary_T, Mpbary_T, Mnbary_T, Xpbary_T, Xnbary_T, Xabary_T, Xhbary_T, &
    !$OMP   Min_D, Min_T, Min_Yp, Max_D, Max_T, Max_Yp )
#endif

    DEALLOCATE( Dbary_T, Tbary_T, Ypbary_T )

    DEALLOCATE( Pbary_T  )
    DEALLOCATE( Sbary_T  )
    DEALLOCATE( Ebary_T  )
    DEALLOCATE( Mpbary_T )
    DEALLOCATE( Mnbary_T )
    DEALLOCATE( Xpbary_T )
    DEALLOCATE( Xnbary_T )
    DEALLOCATE( Xabary_T )
    DEALLOCATE( Xhbary_T )

    IF ( .NOT. UsingExternalEOS ) THEN
       DEALLOCATE( MuonTable )
       DEALLOCATE( HelmholtzTable )
       DEALLOCATE( EOSBary )
    END IF
  
#endif

  END SUBROUTINE FinalizeEquationOfState_TABLE


  SUBROUTINE ApplyChemicalPotentialShift_TABLE( Mpbary_T, Mnbary_T, OS_Mp, OS_Mn )

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

    REAL(dp), INTENT(inout), DIMENSION(:,:,:) :: Mpbary_T, Mnbary_T
    REAL(dp), INTENT(inout) :: OS_Mp, OS_Mn

    REAL(DP), PARAMETER :: neutron_mass = 939.56542052d0
    REAL(DP), PARAMETER :: proton_mass = 938.2720813d0
    REAL(DP), PARAMETER :: dmnp = 1.29333922d0
    REAL(DP) :: min_M, OS_M_new

    ! Apply the shift for proton chemical potential
    IF ( OS_Mp > 0.0d0 ) THEN
      min_M = -0.5d0 * OS_Mp
    ELSE
      min_M = MINVAL( 10.0d0**Mpbary_T )
    ENDIF
    OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + proton_mass + dmnp )
    Mpbary_T     = LOG10( 10.0d0**Mpbary_T - OS_Mp + proton_mass + dmnp + OS_M_new)
    OS_Mp    = OS_M_new

    ! Apply the shift for neutron chemical potential
    IF ( OS_Mn > 0.0d0 ) THEN
      min_M = -0.5d0 * OS_Mn
    ELSE
      min_M = MINVAL( 10.0d0**Mnbary_T )
    ENDIF
    OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + proton_mass + dmnp )
    Mnbary_T     = LOG10( 10.0d0**Mnbary_T - OS_Mn + proton_mass + dmnp + OS_M_new)
    OS_Mn    = OS_M_new

  END SUBROUTINE ApplyChemicalPotentialShift_TABLE


  SUBROUTINE ApplyEquationOfState_TABLE_Scalar &
    ( D, T, Ye, Ym, P, S, E, Mue, Mumu, Mup, Mun, Xp, Xn, Xa, Xh )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: P, S, E, Mue, Mumu, Mup, Mun, Xp, Xn, Xa, Xh

    REAL(DP) :: Pbary, Sbary, Ebary, Pele, Sele, Eele, &
                P_mu, S_mu, E_mu

    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState

    ! Calculate Electron Quantities
    ! Initialize Electron state (Abar and Zbar not needed!!!)
    ElectronState % t = T
    ElectronState % rho = D
    ElectronState % Y_e = Ye
    
    CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

    Eele = ElectronState % e + me / rmu * ergmev * ElectronState % Y_e ! add back mass to internal Energy!
    Pele = ElectronState % p
    Sele = ElectronState % s
    Mue = ElectronState % mu_e

    ! Calculate Muon Quantities
    MuonState % t = T
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)

    E_mu = MuonState % e + mmu / rmu * ergmev * Ym ! add back mass to internal Energy!
    P_mu = MuonState % p
    S_mu = MuonState % s
    Mumu = MuonState % mu
    
    ! Calculate Baryon quantities

#ifdef INVERSION_COMBINED
    ! --- Interpolate Pressure ----------------------------------------

    CALL ComputeDependentVariableTotal_Scalar &
           ( D, T, Ye, Ym, P, Pbary_T, OS_P, &
           UnitP, 0.0_dp, 1.0_dp, 0.0_dp )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariableTotal_Scalar &
           ( D, T, Ye, Ym, S, Sbary_T, OS_S, &
           UnitS, 0.0_dp, 0.0_dp, 1.0_dp )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariableTotal_Scalar &
           ( D, T, Ye, Ym, E, Ebary_T, OS_E, &
           UnitE, 1.0_dp, 0.0_dp, 0.0_dp )
#else
    ! --- Interpolate Pressure ----------------------------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, P, Pbary_T, OS_P, UnitP )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, S, Sbary_T, OS_S, UnitS )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, E, Ebary_T, OS_E, UnitE )
           
    E = E + Eele + E_mu
    P = P + Pele + P_mu
    S = S + Sele + S_mu
    
#endif

    ! --- Interpolate Proton Chemical Potential -----------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, Mup, Mpbary_T, OS_Mp, UnitMp )

    ! --- Interpolate Neutron Chemical Potential ----------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, Mun, Mnbary_T, OS_Mn, UnitMn )

    ! --- Interpolate Proton Mass Fraction ----------------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, Xp, Xpbary_T, OS_Xp, UnitXp )

    ! --- Interpolate Neutron Mass Fraction ---------------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, Xn, Xnbary_T, OS_Xn, UnitXn )

    ! --- Interpolate Alpha Mass Fraction -----------------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, Xa, Xabary_T, OS_Xa, UnitXa )

    ! --- Interpolate Heavy Mass Fraction -----------------------------

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, Xh, Xhbary_T, OS_Xh, UnitXh )

  END SUBROUTINE ApplyEquationOfState_TABLE_Scalar


  SUBROUTINE ApplyEquationOfState_TABLE_Vector &
    ( D, T, Ye, Ym, P, S, E, Mue, Mumu, Mup, Mun, Xp, Xn, Xa, Xh )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: P(1:), S(1:), E(1:), Mue(1:), Mumu(1:), Mup(1:), Mun(1:)
    REAL(DP), INTENT(out) :: Xp(1:), Xn(1:), Xa(1:), Xh(1:)

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ApplyEquationOfState_TABLE_Scalar &
             ( D (iP), T (iP), Ye (iP), Ym (iP), P (iP), S (iP), E (iP), &
             Mue(iP), Mumu(iP), Mup(iP), Mun(iP), Xp(iP), Xn(iP), Xa(iP), Xh(iP) )

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
           ( D_P, E_P, Ye_P, Ym_P, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, OS_E, & 
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
           ( D_P, E_P, Ye_P, Ym_P, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, OS_E, &
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

    D_P = D / ( Gram / Centimeter**3 )
    E_P = E / ( Erg / Gram )
    Ye_P = Ye
    Ym_P = Ym
    
    CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error &
           ( D_P, E_P, Ye_P, Ym_P, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, OS_E, &
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
           ( D_P, E_P, Ye_P, Ym_P, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, OS_E, &
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
    !$OMP MAP( to: D, E, Ye, Ym, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, Guess_Option ) &
    !$OMP MAP( from: T, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup, T_Guess ) &
    !$ACC COPYIN( D, E, Ye, Ym, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, Guess_Option ) &
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
               ( D_P, E_P, Ye_P, Ym_P, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, OS_E, &
               T_Lookup, T_Guess, Error(iP) )

        T(iP) = T_Lookup * Kelvin

      END DO

    ELSE

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup ) &
    !$OMP MAP( to: D, E, Ye, Ym, Dbary_T, Tbary_T, Ypbary_T, Ebary_T ) &
    !$OMP MAP( from: T, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, E_P, Ye_P, Ym_P, T_Lookup ) &
    !$ACC COPYIN( D, E, Ye, Ym, Dbary_T, Tbary_T, Ypbary_T, Ebary_T ) &
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
               ( D_P, E_P, Ye_P, Ym_P, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, OS_E, T_Lookup, &
                 Error(iP) )

        T(iP) = T_Lookup * Kelvin

      END DO

    END IF

    IF ( ANY( Error > 0 ) ) THEN
      DO iP = 1, nP
        IF ( Error(iP) > 0 ) THEN
          CALL DescribeEOSInversionError( Error(iP) )
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
           ( D, Em, Ye+Ym, Ye, Ym, T )

    CALL ComputePressure_TABLE_Scalar &
           ( D, T, Ye+Ym, Ye, Ym, P )

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

    CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError &
           ( D_P, E_P, Ye_P, Ym_P, Dbary_T, Tbary_T, Ypbary_T, Ebary_T, OS_E, &
           T_P )

    T = T_P * UnitT

#ifdef INVERSION_COMBINED

    CALL ComputeDependentVariableTotal_Scalar &
           ( D, T, Ye, Ym, P, Pbary_T, OS_P, &
           UnitP, 0.0_dp, 1.0_dp, 0.0_dp )
#else
    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, P, Pbary_T, OS_P, UnitP )
           
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

    CALL CalculatewlSoundSpeedModule( D, T, Ye, Ym, Dbary_T, Tbary_T, Ypbary_T, &
        Pbary_T, OS_P, Ebary_T, OS_E, 'Energy', HelmholtzTable, MuonTable, Gm, Cs)

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

      CALL CalculatewlSoundSpeedModule( D(iP), T(iP), Ye(iP), Ym(iP), &
          Dbary_T, Tbary_T, Ypbary_T, Pbary_T, OS_P, Ebary_T, OS_E, &
          'Energy', HelmholtzTable, MuonTable, Gm(iP), Cs(iP))

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

    CALL ComputeTemperatureWith_DPY_Single_NoGuess_Error &
           ( D_P, P_P, Ye_P, Ym_P, Dbary_T, Tbary_T, Ypbary_T, Pbary_T, OS_P, T_P, &
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

    CALL ComputeDependentVariableTotal_Scalar &
           ( D, T, Ye, Ym, Em, Ebary_T, OS_E, &
           UnitE, 1.0_dp, 0.0_dp, 0.0_dp )
#else

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, Em, Ebary_T, OS_E, UnitE )
           
    ! Calculate Electron Quantities
    ! Initialize Electron state (Abar and Zbar not needed!!!)
    ElectronState % t = T
    ElectronState % rho = D
    ElectronState % Y_e = Ye
    
    CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)
    Eele = ElectronState % e + me / rmu * ergmev * ElectronState % Y_e ! add back mass to internal Energy!

    ! Calculate Muon Quantities
    MuonState % t = T
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
    E_mu = MuonState % e + mmu / rmu * ergmev * Ym ! add back mass to internal Energy!
           
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

  ! THIS SUBROUTINE IS UNTOUCHED, NOT SURE IF IT'S NEEDED WITH BOTH
  ! MUONS AND ELECTRONS, MY GUESS IS NO
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
    !$OMP MAP( to: D, Ev, Ne ) &
    !$OMP MAP( from: T, Em, Ye, Ym, Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYIN( D, Ev, Ne ) &
    !$ACC COPYOUT( T, Em, Ye, Ym, Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      Em(iP)  = Ev(iP)  / D(iP)              ! --- Internal Energy per Mass
      Ye(iP)  = Ne(iP)  / D(iP) * BaryonMass ! --- Electron Fraction
      Ym(iP) = Nm(iP) / D(iP) * BaryonMass ! --- Muon Fraction

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar &
             ( D(iP), Em(iP), Ye(iP), Ym(iP), T(iP), Error(iP) )

    END DO

    IF ( ANY( Error > 0 ) ) THEN
      DO iP = 1, nP
        IF ( Error(iP) > 0 ) THEN
          CALL DescribeEOSInversionError( Error(iP) )
#if defined(THORNADO_OMP_OL)
          !$OMP TARGET UPDATE FROM &
          !$OMP ( D(iP), Em(iP), Ye(iP), Ym(iP) )
#elif defined(THORNADO_OACC)
          !$ACC UPDATE HOST &
          !$ACC ( D(iP), Em(iP), Ye(iP), Ym(iP) )
#endif
          D_P = D(iP)  / ( Gram / Centimeter**3 )
          E_P = Em(iP) / ( Erg / Gram )
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

    CALL CalculatewlSoundSpeedModule( D, T, Ye, Ym, Dbary_T, Tbary_T, Ypbary_T, &
        Pbary_T, OS_P, Ebary_T, OS_E, 'Energy', HelmholtzTable, MuonTable, Gm, Cs)

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
    ( D, T, Ye, Ym, P, dPdD_Option, dPdT_Option, dPdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: P
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dPdD_Local, dPdT_Local, dPdY_Local
    REAL(DP), POINTER :: dPdD      , dPdT      , dPdYp

#ifdef INVERSION_COMBINED
#else
    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    REAL(DP) :: Pele, P_mu
#endif

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
        dPdYp => dPdY_Option
      ELSE
        dPdYp => dPdY_Local
      END IF

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
             ( D, T, Ye, Ym, P, dPdD, dPdT, dPdYp, Pbary_T, OS_P, &
             UnitP, 1.0_dp, 0.0_dp, 0.0_dp )

    ELSE

#ifdef INVERSION_COMBINED
      CALL ComputeDependentVariableTotal_Scalar &
             ( D, T, Ye, Ym, P, Pbary_T, OS_P, &
             UnitP, 1.0_dp, 0.0_dp, 0.0_dp )
#else
      CALL ComputeDependentVariableBaryons_Scalar &
             ( D, T, Ye+Ym, P, Pbary_T, OS_P, UnitP )

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
    ( D, T, Ye, Ym, P, dPdD_Option, dPdT_Option, dPdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: P(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dPdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: &
      dPdD_Local, dPdT_Local, dPdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: &
      dPdD      , dPdT      , dPdYp

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
        dPdYp(1:nP) => dPdY_Option(:)
      ELSE
        dPdYp(1:nP) => dPdY_Local(:)
      END IF

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Vector &
             ( D, T, Ye, Ym, P, dPdD, dPdT, dPdYp, Pbary_T, OS_P, &
             UnitP, 1.0_dp, 0.0_dp, 0.0_dp )

    ELSE

      CALL ComputeDependentVariableTotal_Vector &
             ( D, T, Ye, Ym, P, Pbary_T, OS_P, &
             UnitP, 1.0_dp, 0.0_dp, 0.0_dp )

    END IF

  END SUBROUTINE ComputePressure_TABLE_Vector


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Scalar &
    ( D, T, Ye, Ym, E, dEdD_Option, dEdT_Option, dEdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
    REAL(DP), INTENT(out)                   :: E
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dEdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), TARGET  :: dEdD_Local, dEdT_Local, dEdY_Local
    REAL(DP), POINTER :: dEdD,       dEdT,       dEdY

#ifdef INVERSION_COMBINED
#else
    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    REAL(DP) :: Eele, E_mu
#endif

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

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
             ( D, T, Ye, Ym, E, dEdD, dEdT, dEdY, Ebary_T, OS_E, &
             UnitE, 1.0_dp, 0.0_dp, 0.0_dp )

    ELSE

#ifdef INVERSION_COMBINED

    CALL ComputeDependentVariableTotal_Scalar &
           ( D, T, Ye, Ym, Em, Ebary_T, OS_E, &
           UnitE, 1.0_dp, 0.0_dp, 0.0_dp )
#else

    CALL ComputeDependentVariableBaryons_Scalar &
           ( D, T, Ye+Ym, E, Ebary_T, OS_E, UnitE )
           
    ! Calculate Electron Quantities
    ! Initialize Electron state (Abar and Zbar not needed!!!)
    ElectronState % t = T
    ElectronState % rho = D
    ElectronState % Y_e = Ye
    
    CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)
    Eele = ElectronState % e + me / rmu * ergmev * ElectronState % Y_e ! add back mass to internal Energy!

    ! Calculate Muon Quantities
    MuonState % t = T
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)
    E_mu = MuonState % e + mmu / rmu * ergmev * Ym ! add back mass to internal Energy!
           
    E = E + Eele + E_mu
#endif
    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Vector &
    ( D, T, Ye, Ym, E, dEdD_Option, dEdT_Option, dEdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
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

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Vector &
             ( D, T, Ye, Ym, E, dEdD, dEdT, dEdY, Ebary_T, OS_E, &
             UnitE, 1.0_dp, 0.0_dp, 0.0_dp )

    ELSE

      CALL ComputeDependentVariableTotal_Vector &
             ( D, T, Ye, Ym, E, Ebary_T, OS_E, &
             UnitE, 1.0_dp, 0.0_dp, 0.0_dp )

    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Vector

  ! Notice that only E, S, and P have the option of calculating derivatives.
  ! I do not see why you would need derivatives of anything else. In case you
  ! do, then there is no electron or muon contribution for those, so you might have to make 
  ! new subroutines, and should not use ComputeDependentVariableAndDerivativesTotal_TABLE_Vector,
  ! but instead create a new one, somthing like ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector

  ! The Mue and Mumu subroutines are special since you only need electron (muon) info
  SUBROUTINE ComputeElectronChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, M )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye
    REAL(DP), INTENT(out) :: M

    TYPE(ElectronStateType) :: ElectronState
    
    ! Calculate Electron Quantities
    ! Initialize Electron state (Abar and Zbar not needed!!!)
    ElectronState % t = T
    ElectronState % rho = D
    ElectronState % Y_e = Ye
    
    CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

    M = ElectronState % mu_e
    
  END SUBROUTINE ComputeElectronChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeElectronChemicalPotential_TABLE_Vector &
    ( D, T, Ye, M )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:)
    REAL(DP), INTENT(out)                   :: M(1:)

    INTEGER :: iP, nP

    DO iP=1,nP
      CALL ComputeElectronChemicalPotential_TABLE_Scalar &
              (D(iP), T(iP), Ye(iP), M(iP))
    END DO

  END SUBROUTINE ComputeElectronChemicalPotential_TABLE_Vector

  SUBROUTINE ComputeMuonChemicalPotential_TABLE_Scalar &
    ( D, T, Ym, M )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ym
    REAL(DP), INTENT(out) :: M

    TYPE(MuonStateType) :: MuonState
    
    ! Calculate Muon Quantities
    ! Initialize Muon state (Abar and Zbar not needed!!!)
    MuonState % t = T
    MuonState % rhoym = D * Ym
    
    CALL FullMuonEOS(MuonTable, MuonState)

    M = MuonState % mu
    
  END SUBROUTINE ComputeMuonChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeMuonChemicalPotential_TABLE_Vector &
    ( D, T, Ym, M )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: M(1:)

    INTEGER :: iP, nP

    DO iP=1,nP
      CALL ComputeMuonChemicalPotential_TABLE_Scalar &
              (D(iP), T(iP), Ym(iP), M(iP))
    END DO

  END SUBROUTINE ComputeMuonChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeProtonChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye+Ym, M, dMdD, dMdT, dMdY, Mpbary_T, OS_Mp, UnitMp )

    ELSE

      CALL ComputeDependentVariableBaryons_Scalar &
             ( D, T, Ye+Ym, M, Mpbary_T, OS_Mp, UnitMp )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE_Scalar


  ! THIS SUBROUTINE IS UNTOUCHED, NOT SURE HOW TO HANDLE IT
  SUBROUTINE ComputeProtonChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye+Ym, M, dMdD, dMdT, dMdY, Mpbary_T, OS_Mp, UnitMp )

    ELSE

      CALL ComputeDependentVariableBaryons_Vector &
             ( D, T, Ye+Ym, M, Mpbary_T, OS_Mp, UnitMp )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye+Ym, M, dMdD, dMdT, dMdY, Mnbary_T, OS_Mn, UnitMn )

    ELSE

      CALL ComputeDependentVariableBaryons_Scalar &
             ( D, T, Ye+Ym, M, Mnbary_T, OS_Mn, UnitMn )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye+Ym, M, dMdD, dMdT, dMdY, Mnbary_T, OS_Mn, UnitMn )

    ELSE

      CALL ComputeDependentVariableBaryons_Vector &
             ( D, T, Ye+Ym, M, Mnbary_T, OS_Mn, UnitMn )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeProtonMassFraction_TABLE_Scalar &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye+Ym, X, dXdD, dXdT, dXdY, Xpbary_T, OS_Xp, UnitXp )

    ELSE

      CALL ComputeDependentVariableBaryons_Scalar &
             ( D, T, Ye+Ym, X, Xpbary_T, OS_Xp, UnitXp )

    END IF

  END SUBROUTINE ComputeProtonMassFraction_TABLE_Scalar


  ! THIS SUBROUTINE IS UNTOUCHED, NOT SURE HOW TO HANDLE IT
  SUBROUTINE ComputeProtonMassFraction_TABLE_Vector &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye+Ym, X, dXdD, dXdT, dXdY, Xpbary_T, OS_Xp, UnitXp )

    ELSE

      CALL ComputeDependentVariableBaryons_Vector &
             ( D, T, Ye+Ym, X, Xpbary_T, OS_Xp, UnitXp )

    END IF

  END SUBROUTINE ComputeProtonMassFraction_TABLE_Vector


  SUBROUTINE ComputeNeutronMassFraction_TABLE_Scalar &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye+Ym, X, dXdD, dXdT, dXdY, Xnbary_T, OS_Xn, UnitXn )

    ELSE

      CALL ComputeDependentVariableBaryons_Scalar &
             ( D, T, Ye+Ym, X, Xnbary_T, OS_Xn, UnitXn )

    END IF

  END SUBROUTINE ComputeNeutronMassFraction_TABLE_Scalar


  ! THIS SUBROUTINE IS UNTOUCHED, NOT SURE HOW TO HANDLE IT
  SUBROUTINE ComputeNeutronMassFraction_TABLE_Vector &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
             ( D, T, Ye+Ym, X, dXdD, dXdT, dXdY, Xnbary_T, OS_Xn, UnitXn )

    ELSE

      CALL ComputeDependentVariableBaryons_Vector &
             ( D, T, Ye+Ym, X, Xnbary_T, OS_Xn, UnitXn )

    END IF

  END SUBROUTINE ComputeNeutronMassFraction_TABLE_Vector


  SUBROUTINE ComputeHeavyMassFraction_TABLE_Scalar &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
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

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D, T, Ye+Ym, X, dXdD, dXdT, dXdY, Xhbary_T, OS_Xh, UnitXh )

    ELSE

      CALL ComputeDependentVariableBaryons_Scalar &
             ( D, T, Ye+Ym, X, Xhbary_T, OS_Xh, UnitXh )

    END IF

  END SUBROUTINE ComputeHeavyMassFraction_TABLE_Scalar

  SUBROUTINE ComputeHeavyMassFraction_TABLE_Vector &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: X(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dXdD_Local, dXdT_Local, dXdY_Local
    REAL(DP), DIMENSION(:), POINTER :: &
      dXdD      , dXdT      , dXdY

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdY_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dXdD_Local(nP), dXdT_Local(nP), dXdY_Local(nP) )

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
             ( D, T, Ye+Ym, X, dXdD, dXdT, dXdY, Xh_T, OS_Xh, Units_V = UnitXh )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Ye+Ym, X, Xh_T, OS_Xh, Units_V = UnitXh )

    END IF

  END SUBROUTINE ComputeHeavyMassFraction_TABLE_Vector

  SUBROUTINE ComputeHeavyMassNumber_TABLE_Scalar &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdY_Option )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)                    :: D, T, Ye, Ym
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
             ( D, T, Ye+Ym, X, dXdD, dXdT, dXdY, Ah_T, OS_Ah, Units_V = UnitAh )

    ELSE

      CALL ComputeDependentVariable_TABLE_Scalar &
             ( D, T, Ye+Ym, X, Ah_T, OS_Ah, Units_V = UnitAh )

    END IF

  END SUBROUTINE ComputeHeavyMassNumber_TABLE_Scalar

  SUBROUTINE ComputeHeavyMassNumber_TABLE_Vector &
    ( D, T, Ye, Ym, X, dXdD_Option, dXdT_Option, dXdY_Option )

    REAL(DP), INTENT(in)                    :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out)                   :: X(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdD_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdT_Option(1:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dXdY_Option(1:)

    LOGICAL :: ComputeDerivatives
    INTEGER :: nP
    REAL(DP), DIMENSION(:), TARGET, ALLOCATABLE  :: &
      dXdD_Local, dXdT_Local, dXdY_Local
    REAL(DP), DIMENSION(:), POINTER :: &
      dXdD      , dXdT      , dXdY

    ComputeDerivatives &
      =      PRESENT( dXdD_Option ) &
        .OR. PRESENT( dXdT_Option ) &
        .OR. PRESENT( dXdY_Option )

    IF( ComputeDerivatives )THEN

      nP = SIZE( D )
      ALLOCATE( dXdD_Local(nP), dXdT_Local(nP), dXdY_Local(nP) )

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
             ( D, T, Ye+Ym, X, dXdD, dXdT, dXdY, Ah_T, OS_Ah, Units_V = UnitAh )

    ELSE

      CALL ComputeDependentVariable_TABLE_Vector &
             ( D, T, Ye+Ym, X, Ah_T, OS_Ah, Units_V = UnitAh )

    END IF

  END SUBROUTINE ComputeHeavyMassNumber_TABLE_Vector


  SUBROUTINE ComputeNeutrinoChemicalPotential_TABLE_Scalar &
    ( D, T, Ye, Ym, Mnue, Mnumu )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: Mnue, Mnumu

    REAL(DP) :: D_P, T_P, Ye_P, Ym_P
    INTEGER  :: iD, iT, iY
    REAL(DP) :: dD, dT, dY
    REAL(DP) :: Mue, Mumu, Mup, Mun

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY

    CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P+Ym_P, Mun )
      
    CALL ComputeProtonChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P+Ym_P, Mup )
      
    CALL ComputeElectronChemicalPotential_TABLE_Scalar &
      ( D, T, Ye_P, Mue )
      
    CALL ComputeMuonChemicalPotential_TABLE_Scalar &
      ( D, T, Ym_P, Mumu )

    Mue  = Mue  * UnitMl
    Mumu = Mumu * UnitMl
    Mup  = Mup  * UnitMp
    Mun  = Mun  * UnitMn

    Mnue  = ( Mue  + Mup ) - Mun
    Mnumu = ( Mumu + Mup ) - Mun

#else

    Mnue  = Zero
    Mnumu = Zero

#endif

  END SUBROUTINE ComputeNeutrinoChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeNeutrinoChemicalPotential_TABLE_Vector &
    ( D, T, Ye, Ym, Mnue, Mnumu )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: Mnue(1:), Mnumu(1:)

    REAL(DP) :: D_P, T_P, Ye_P, Ym_P
    INTEGER  :: iD, iT, iY
    REAL(DP) :: dD, dT, dY
    REAL(DP) :: Mue, Mumu, Mup, Mun

    INTEGER  :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP PRIVATE( D_P, T_P, Yp_P, iD, iT, iY, dD, dT, dY, Mue, Mumu, Mup, Mun ) &
    !$OMP MAP( to: D, T, Y, Dbary_T, Tbary_T, Ypbary_T, OS_Mp, OS_Mn, Mpbary_T, Mnbary_T ) &
    !$OMP MAP( from: Mnue, Mnumu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC PRIVATE( D_P, T_P, Yp_P, iD, iT, iY, dD, dT, dY, Mue, Mumu, Mup, Mun ) &
    !$ACC COPYIN( D, T, Y, Dbary_T, Tbary_T, Ypbary_T, OS_Mp, OS_Mn, Mpbary_T, Mnbary_T ) &
    !$ACC COPYOUT( Mnue, Mnumu )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, T_P, Yp_P, iD, iT, iY, dD, dT, dY, Mue, Mumu, Mup, Mun )
#endif
    DO iP = 1, nP

      D_P = D(iP) / UnitD
      T_P = T(iP) / UnitT
      Ye_P = Ye(iP) / UnitY
      Ym_P = Ym(iP) / UnitY

      CALL ComputeNeutronChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P+Ym_P, Mun )
        
      CALL ComputeProtonChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P+Ym_P, Mup )
        
      CALL ComputeElectronChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ye_P, Mue )
        
      CALL ComputeMuonChemicalPotential_TABLE_Scalar &
        ( D_P, T_P, Ym_P, Mumu )

      Mue  = Mue  * UnitMl
      Mumu = Mumu * UnitMl
      Mup  = Mup  * UnitMp
      Mun  = Mun  * UnitMn

      Mnue(iP)  = ( Mue  + Mup ) - Mun
      Mnumu(iP) = ( Mumu + Mup ) - Mun

    END DO

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP MAP( from: Mnue, Mnumu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC COPYOUT( Mnue, Mnumu )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP
      Mnue(iP) = Zero
      Mnumu(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeDependentVariableBaryons_Scalar &
    ( D, T, Yp, V, Vbary_T, OS_V, Units_V )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Yp
    REAL(DP), INTENT(out) :: V
    REAL(DP), INTENT(in)  :: Vbary_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Yp_P, V_P

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Yp_P = Yp / UnitY

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, Dbary_T, Tbary_T, Ypbary_T, OS_V, Vbary_T, V_P )

    V = V_P * Units_V

#else

    V = Zero

#endif

  END SUBROUTINE ComputeDependentVariableBaryons_Scalar


  SUBROUTINE ComputeDependentVariableBaryons_Vector &
    ( D, T, Yp, V, Vbary_T, OS_V, Units_V )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Yp(1:)
    REAL(DP), INTENT(out) :: V(1:)
    REAL(DP), INTENT(in)  :: Vbary_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    INTEGER  :: iP, nP
    REAL(DP) :: D_P, T_P, Yp_P, V_P

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( D_P, T_P, Yp_P, V_P ) &
      !$OMP MAP( to: D, T, Y, Dbary_T, Tbary_T, Ypbary_T, OS_V, Vbary_T ) &
      !$OMP MAP( from: V )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( D_P, T_P, Yp_P, V_P ) &
      !$ACC COPYIN( D, T, Y, Dbary_T, Tbary_T, Ypbary_T, OS_V, Vbary_T ) &
      !$ACC COPYOUT( V )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( D_P, T_P, Yp_P, V_P )
#endif
    DO iP = 1, nP

      D_P = D(iP) / UnitD
      T_P = T(iP) / UnitT
      Yp_P = Yp(iP) / UnitY

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, Dbary_T, Tbary_T, Ypbary_T, OS_V, Vbary_T, V_P )

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

  END SUBROUTINE ComputeDependentVariableBaryons_Vector

  SUBROUTINE ComputeDependentVariableTotal_Scalar &
    ( D, T, Ye, Ym, V, Vbary_T, OS_V, Units_V, &
    InputE, InputP, InputS )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: V
    REAL(DP), INTENT(in)  :: Vbary_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    REAL(DP), INTENT(in)  :: InputE, InputP, InputS

    REAL(DP) :: D_P, T_P, Ye_P, Ym_P, V_P
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: E_leptons, P_leptons, S_leptons
    REAL(DP) :: Vtot(2,2,2)
    
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps

    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY

    Ye_over_Yp = Ye_P/(Ye_P + Ym_P)
    Ym_over_Yp = Ym_P/(Ye_P + Ym_P)

    ! Now bracket the points 
    SizeDs = SIZE( Dbary_T )
    SizeTs = SIZE( Tbary_T )
    SizeYps = SIZE( Ypbary_T )

    iD = Index1D_Log( D_P, Dbary_T )
    iT = Index1D_Log( T_P, Tbary_T )
    iYp = Index1D_Lin( Yp_P, Ypbary_T )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    DO iL_T=1,2
      DO iL_D=1,2
        DO iL_Y=1,2
          ElectronState % t = Tbary_T(iT+iL_T-1)
          ElectronState % rho = Dbary_T(iD+iL_D-1)
          ElectronState % y_e = Ypbary_T(iYp+iL_Y-1) * Ye_over_Yp
          
          ! CALL FullHelmEOS(1, HelmholtzTable, ElectronState, .false., .false.)
          CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

          MuonState % t = Tbary_T(iT+iL_T-1)
          MuonState % rhoym = Dbary_T(iD+iL_D-1) * Ypbary_T(iYp+iL_Y-1) * Ym_over_Yp
          
          CALL FullMuonEOS(MuonTable, MuonState)
          
          E_leptons = ElectronState % e + me / rmu * ergmev * ElectronState % y_e + &
                        MuonState % e + mmu / rmu * ergmev * Ypbary_T(iYp+iL_Y-1) * Ym_over_Yp
          P_leptons = ElectronState % p + MuonState % p
          S_leptons = ElectronState % s + MuonState % s
          
          Vtot(iL_D,iL_T,iL_Y) = Vbary_T(iD+iL_D-1,iT+iL_T-1,iYp+iL_Y-1) + &
              InputE*E_leptons + InputP*P_leptons + InputS*S_leptons
              
        END DO
      END DO
    END DO

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, Dbary_T(iD:iD+1), Tbary_T(iT:iT+1), Ypbary_T(iYp:iYp+1), &
           OS_V, Vtot, V_P )

    V = V_P * Units_V

#else

    V = Zero

#endif

  END SUBROUTINE ComputeDependentVariableTotal_Scalar


  SUBROUTINE ComputeDependentVariableTotal_Vector &
    ( D, T, Ye, Ym, V, Vbary_T, OS_V, Units_V, &
    InputE, InputP, InputS )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: V(1:)
    REAL(DP), INTENT(in)  :: Vbary_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    REAL(DP), INTENT(in)  :: InputE, InputP, InputS

    INTEGER  :: iP, nP
    REAL(DP) :: D_P, T_P, Ye_P, Ym_P, V_P
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: E_leptons, P_leptons, S_leptons
    REAL(DP) :: Vtot(2,2,2)
    
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps

    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
        
    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP PRIVATE( D_P, T_P, Yp_P, V_P ) &
      !$OMP MAP( to: D, T, Y, Dbary_T, Tbary_T, Ypbary_T, OS_V, Vbary_T ) &
      !$OMP MAP( from: V )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC PRIVATE( D_P, T_P, Yp_P, V_P ) &
      !$ACC COPYIN( D, T, Y, Dbary_T, Tbary_T, Ypbary_T, OS_V, Vbary_T ) &
      !$ACC COPYOUT( V )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO &
      !$OMP PRIVATE( D_P, T_P, Yp_P, V_P )
#endif
    DO iP = 1, nP

      D_P = D(iP) / UnitD
      T_P = T(iP) / UnitT
      Yp_P = Yp(iP) / UnitY
      Ye_P = Ye(iP) / UnitY
      Ym_P = Ym(iP) / UnitY

      Ye_over_Yp = Ye_P/(Ye_P + Ym_P)
      Ym_over_Yp = Ym_P/(Ye_P + Ym_P)
      
      ! Now bracket the points 
      SizeDs = SIZE( Dbary_T )
      SizeTs = SIZE( Tbary_T )
      SizeYps = SIZE( Ypbary_T )

      iD = Index1D_Log( D_P, Dbary_T )
      iT = Index1D_Log( T_P, Tbary_T )
      iYp = Index1D_Lin( Yp_P, Ypbary_T )

      iD = MIN( MAX( 1, iD ), SizeDs - 1 )
      iT = MIN( MAX( 1, iT ), SizeTs - 1 )
      iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )
      
      ! now calculate muon and electron contribution on that specific point
      DO iL_T=1,2
        DO iL_D=1,2
          DO iL_Y=1,2
            ElectronState % t = Tbary_T(iT+iL_T-1)
            ElectronState % rho = Dbary_T(iD+iL_D-1)
            ElectronState % y_e = Ypbary_T(iYp+iL_Y-1) * Ye_over_Yp
            
            ! CALL FullHelmEOS(1, HelmholtzTable, ElectronState, .false., .false.)
            CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

            MuonState % t = Tbary_T(iT+iL_T-1)
            MuonState % rhoym = Dbary_T(iD+iL_D-1) * Ypbary_T(iYp+iL_Y-1) * Ym_over_Yp
            
            CALL FullMuonEOS(MuonTable, MuonState)
            
            E_leptons = ElectronState % e + me / rmu * ergmev * ElectronState % y_e + &
                          MuonState % e + mmu / rmu * ergmev * Ypbary_T(iYp+iL_Y-1) * Ym_over_Yp
            P_leptons = ElectronState % p + MuonState % p
            S_leptons = ElectronState % s + MuonState % s
            
            Vtot(iL_D,iL_T,iL_Y) = Vbary_T(iD+iL_D-1,iT+iL_T-1,iYp+iL_Y-1) + &
                InputE*E_leptons + InputP*P_leptons + InputS*S_leptons
              
          END DO
        END DO
      END DO

      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( D_P, T_P, Yp_P, Dbary_T(iD:iD+1), Tbary_T(iT:iT+1), Ypbary_T(iYp:iYp+1), &
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

  END SUBROUTINE ComputeDependentVariableTotal_Vector


  SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
    ( D, T, Ye, Ym, V, dVdD, dVdT, dVdYp, Vbary_T, OS_V, Units_V, &
    InputE, InputP, InputS)

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Ye, Ym
    REAL(DP), INTENT(out) :: V, dVdD, dVdT, dVdYp
    REAL(DP), INTENT(in)  :: Vbary_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    REAL(DP), INTENT(in)  :: InputE, InputP, InputS

    REAL(DP) :: D_P, T_P, Ye_P, Ym_P, V_P, dV_P(3)
    REAL(DP) :: Ye_over_Yp, Ym_over_Yp
    REAL(DP) :: E_leptons, P_leptons, S_leptons
    REAL(DP) :: Vtot(2,2,2)
    
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps
    
    TYPE(ElectronStateType) :: ElectronState
    TYPE(MuonStateType) :: MuonState
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Ye_P = Ye / UnitY
    Ym_P = Ym / UnitY

    Ye_over_Yp = Ye_P/(Ye_P + Ym_P)
    Ym_over_Yp = Ym_P/(Ye_P + Ym_P)

    ! Now bracket the points 
    SizeDs = SIZE( Dbary_T )
    SizeTs = SIZE( Tbary_T )
    SizeYps = SIZE( Ypbary_T )

    iD = Index1D_Log( D_P, Dbary_T )
    iT = Index1D_Log( T_P, Tbary_T )
    iYp = Index1D_Lin( Yp_P, Ypbary_T )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )
    
    ! now calculate muon and electron contribution on that specific point
    DO iL_T=1,2
      DO iL_D=1,2
        DO iL_Y=1,2
          ElectronState % t = Tbary_T(iT+iL_T-1)
          ElectronState % rho = Dbary_T(iD+iL_D-1)
          ElectronState % y_e = Ypbary_T(iYp+iL_Y-1) * Ye_over_Yp
          
          ! CALL FullHelmEOS(1, HelmholtzTable, ElectronState, .false., .false.)
          CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

          MuonState % t = Tbary_T(iT+iL_T-1)
          MuonState % rhoym = Dbary_T(iD+iL_D-1) * Ypbary_T(iYp+iL_Y-1) * Ym_over_Yp
          
          CALL FullMuonEOS(MuonTable, MuonState)
          
          E_leptons = ElectronState % e + me / rmu * ergmev * ElectronState % y_e + &
                        MuonState % e + mmu / rmu * ergmev * Ypbary_T(iYp+iL_Y-1) * Ym_over_Yp
          P_leptons = ElectronState % p + MuonState % p
          S_leptons = ElectronState % s + MuonState % s
          
          Vtot(iL_D,iL_T,iL_Y) = Vbary_T(iD+iL_D-1,iT+iL_T-1,iYp+iL_Y-1) + &
              InputE*E_leptons + InputP*P_leptons + InputS*S_leptons

        END DO
      END DO
    END DO
    
    CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Ye_P+Ym_P, Dbary_T(iD:iD+1), Tbary_T(iT:iT+1), Ypbary_T(iYp:iYp+1), &
           OS_V, Vtot, V_P, dV_P )
             
    V = V_P * Units_V

    dVdD = dV_P(1) * Units_V / UnitD
    dVdT = dV_P(2) * Units_V / UnitT
    dVdYp = dV_P(3) * Units_V / UnitY

#else

    V    = Zero
    dVdD = Zero
    dVdT = Zero
    dVdYp = Zero

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Vector &
    ( D, T, Ye, Ym, V, dVdD, dVdT, dVdYp, Vbary_T, OS_V, Units_V, &
    InputE, InputP, InputS )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Ye(1:), Ym(1:)
    REAL(DP), INTENT(out) :: V(1:), dVdD(1:), dVdT(1:), dVdYp(1:)
    REAL(DP), INTENT(in)  :: Vbary_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V
    REAL(DP), INTENT(in)  :: InputE, InputP, InputS

    INTEGER :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP MAP( to: D, T, Y, OS_V, Vbary_T ) &
      !$OMP MAP( from: V, dVdD, dVdT, dVdYp )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Y, OS_V, Vbary_T ) &
      !$ACC COPYOUT( V, dVdD, dVdT, dVdYp )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeDependentVariableAndDerivativesTotal_TABLE_Scalar &
             ( D(iP), T(iP), Ye(iP), Ym(iP), V(iP), dVdD(iP), dVdT(iP), dVdYp(iP), &
               Vbary_T, OS_V, Units_V, InputE, InputP, InputS )

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
      dVdYp(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesTotal_TABLE_Vector

  SUBROUTINE ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
    ( D, T, Yp, V, dVdD, dVdT, dVdYp, Vbary_T, OS_V, Units_V )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Yp
    REAL(DP), INTENT(out) :: V, dVdD, dVdT, dVdYp
    REAL(DP), INTENT(in)  :: Vbary_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Yp_P, V_P, dV_P(3)

    INTEGER  :: iD, iT, iYp, iL_D, iL_Y, iL_T
    INTEGER  :: SizeDs, SizeTs, SizeYps
    
#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Yp_P = Yp / UnitY

    CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
           ( D_P, T_P, Yp_P, Dbary_T, Tbary_T, Ypbary_T, OS_V, Vbary_T, V_P, dV_P )
             
    V = V_P * Units_V

    dVdD = dV_P(1) * Units_V / UnitD
    dVdT = dV_P(2) * Units_V / UnitT
    dVdYp = dV_P(3) * Units_V / UnitY

#else

    V    = Zero
    dVdD = Zero
    dVdT = Zero
    dVdYp = Zero

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar


  SUBROUTINE ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector &
    ( D, T, Yp, V, dVdD, dVdT, dVdYp, Vbary_T, OS_V, Units_V )

    REAL(DP), INTENT(in)  :: D(1:), T(1:), Yp(1:)
    REAL(DP), INTENT(out) :: V(1:), dVdD(1:), dVdT(1:), dVdYp(1:)
    REAL(DP), INTENT(in)  :: Vbary_T(1:,1:,1:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    INTEGER :: iP, nP

    nP = SIZE( D )

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP MAP( to: D, T, Y, OS_V, Vbary_T ) &
      !$OMP MAP( from: V, dVdD, dVdT, dVdYp )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC COPYIN( D, T, Y, OS_V, Vbary_T ) &
      !$ACC COPYOUT( V, dVdD, dVdT, dVdYp )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeDependentVariableAndDerivativesBaryons_TABLE_Scalar &
             ( D(iP), T(iP), Yp(iP), V(iP), dVdD(iP), dVdT(iP), dVdYp(iP), &
               Vbary_T, OS_V, Units_V )

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
      dVdYp(iP) = Zero
    END DO

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivativesBaryons_TABLE_Vector
  
END MODULE EquationOfStateComponentsModule_TABLE
