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
    ComputeTemperatureWith_DEY, &
    ComputeTemperatureWith_DPY, &
    DescribeEOSInversionError
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateDifferentiateSingleVariable

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
  USE DeviceModule, ONLY: &
    QueryOnGpu

  IMPLICIT NONE
  PRIVATE

  CHARACTER(256) :: &
    EquationOfStateTableName
  INTEGER :: &
    iD_T, iT_T, iY_T, &
    iP_T, iS_T, iE_T, iMe_T, iMp_T, iMn_T, &
    iXp_T, iXn_T, iXa_T, iXh_T, iGm_T
  INTEGER, DIMENSION(3) :: &
    LogInterp
  REAL(DP) :: &
    UnitD, UnitT, UnitY, &
    UnitP, UnitS, UnitE, UnitMe, UnitMp, UnitMn, &
    UnitXp, UnitXn, UnitXa, UnitXh, UnitGm
  REAL(DP), PUBLIC :: &
    OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, &
    OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm
  REAL(DP), PARAMETER :: &
    BaryonMass = AtomicMassUnit
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    Ds_T, Ts_T, Ys_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: &
    Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, &
    Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T
#ifdef MICROPHYSICS_WEAKLIB
  LOGICAL :: UsingExternalEOS
  TYPE(EquationOfStateTableType), POINTER :: EOS
#endif

  REAL(DP), PUBLIC :: MinD, MinT, MinY
  REAL(DP), PUBLIC :: MaxD, MaxT, MaxY

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
  !$OMP ( Ds_T, Ts_T, Ys_T, &
  !$OMP   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, UnitMe, UnitMp, UnitMn, &
  !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, &
  !$OMP   OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$OMP   OS_Xa, OS_Xh, OS_Gm, &
  !$OMP   Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, &
  !$OMP   Xas_T, Xhs_T, Gms_T, &
  !$OMP   MinD, MinT, MinY, MaxD, MaxT, MaxY )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( Ds_T, Ts_T, Ys_T, &
  !$ACC   UnitD, UnitT, UnitY, UnitP, UnitS, UnitE, UnitMe, UnitMp, UnitMn, &
  !$ACC   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, &
  !$ACC   OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, &
  !$ACC   OS_Xa, OS_Xh, OS_Gm, &
  !$ACC   Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, &
  !$ACC   Xas_T, Xhs_T, Gms_T, &
  !$ACC   MinD, MinT, MinY, MaxD, MaxT, MaxY )
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

    ALLOCATE( Ds_T(EOS % TS % nPoints(iD_T)) )
    Ds_T = EOS % TS % States(iD_T) % Values

    MinD = MINVAL( Ds_T ) * Gram / Centimeter**3
    MaxD = MAXVAL( Ds_T ) * Gram / Centimeter**3

    ALLOCATE( Ts_T(EOS % TS % nPoints(iT_T)) )
    Ts_T = EOS % TS % States(iT_T) % Values

    MinT = MINVAL( Ts_T ) * Kelvin
    MaxT = MAXVAL( Ts_T ) * Kelvin

    ALLOCATE( Ys_T(EOS % TS % nPoints(iY_T)) )
    Ys_T = EOS % TS % States(iY_T) % Values

    MinY = MINVAL( Ys_T )
    MaxY = MAXVAL( Ys_T )

    LogInterp &
      = EOS % TS % LogInterp

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
      ( Ps_T (1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Ss_T (1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Es_T (1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Mes_T(1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Mps_T(1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Mns_T(1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Xps_T(1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Xns_T(1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Xas_T(1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Xhs_T(1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )
    ALLOCATE &
      ( Gms_T(1:EOS % DV % nPoints(1), &
              1:EOS % DV % nPoints(2), &
              1:EOS % DV % nPoints(3)) )

    Ps_T (:,:,:) = EOS % DV % Variables(iP_T ) % Values(:,:,:)
    Ss_T (:,:,:) = EOS % DV % Variables(iS_T ) % Values(:,:,:)
    Es_T (:,:,:) = EOS % DV % Variables(iE_T ) % Values(:,:,:)
    Mes_T(:,:,:) = EOS % DV % Variables(iMe_T) % Values(:,:,:)
    Mps_T(:,:,:) = EOS % DV % Variables(iMp_T) % Values(:,:,:)
    Mns_T(:,:,:) = EOS % DV % Variables(iMn_T) % Values(:,:,:)
    Xps_T(:,:,:) = EOS % DV % Variables(iXp_T) % Values(:,:,:)
    Xns_T(:,:,:) = EOS % DV % Variables(iXn_T) % Values(:,:,:)
    Xas_T(:,:,:) = EOS % DV % Variables(iXa_T) % Values(:,:,:)
    Xhs_T(:,:,:) = EOS % DV % Variables(iXh_T) % Values(:,:,:)
    Gms_T(:,:,:) = EOS % DV % Variables(iGm_T) % Values(:,:,:)

    CALL InitializeEOSInversion &
           ( Ds_T, Ts_T, Ys_T, &
             10.0d0**( Es_T ) - OS_E, &
             10.0d0**( Ps_T ) - OS_P, &
             10.0d0**( Ss_T ) - OS_S, &
             Verbose_Option = Verbose )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO &
    !$OMP ( Ds_T, Ts_T, Ys_T, &
    !$OMP   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMe, UnitMp, UnitMn, &
    !$OMP   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, OS_P, OS_S, OS_E, OS_Me, &
    !$OMP   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm, Ps_T, Ss_T, &
    !$OMP   Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T, &
    !$OMP   MinD, MinT, MinY, MaxD, MaxT, MaxY )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( Ds_T, Ts_T, Ys_T, &
    !$ACC   UnitD, UnitT, UnitY, UnitP, UnitE, UnitMe, UnitMp, UnitMn, &
    !$ACC   UnitXp, UnitXn, UnitXa, UnitXh, UnitGm, OS_P, OS_S, OS_E, OS_Me, &
    !$ACC   OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm, Ps_T, Ss_T, &
    !$ACC   Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T, &
    !$ACC   MinD, MinT, MinY, MaxD, MaxT, MaxY )
#endif

#endif

  END SUBROUTINE InitializeEquationOfState_TABLE


  SUBROUTINE FinalizeEquationOfState_TABLE

#ifdef MICROPHYSICS_WEAKLIB

    DEALLOCATE( Ds_T, Ts_T, Ys_T )

    DEALLOCATE( Ps_T  )
    DEALLOCATE( Ss_T  )
    DEALLOCATE( Es_T  )
    DEALLOCATE( Mes_T )
    DEALLOCATE( Mps_T )
    DEALLOCATE( Mns_T )
    DEALLOCATE( Xps_T )
    DEALLOCATE( Xns_T )
    DEALLOCATE( Xas_T )
    DEALLOCATE( Xhs_T )
    DEALLOCATE( Gms_T )

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

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, P, Ps_T, OS_P, Units_V = UnitP )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, S, Ss_T, OS_S, Units_V = UnitS )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, E, Es_T, OS_E, Units_V = UnitE )

    ! --- Interpolate Electron Chemical Potential ---------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Me, Mes_T, OS_Me, Units_V = UnitMe )

    ! --- Interpolate Proton Chemical Potential -----------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Mp, Mps_T, OS_Mp, Units_V = UnitMp )

    ! --- Interpolate Neutron Chemical Potential ----------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Mn, Mns_T, OS_Mn, Units_V = UnitMn )

    ! --- Interpolate Proton Mass Fraction ----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Xp, Xps_T, OS_Xp, Units_V = UnitXp )

    ! --- Interpolate Neutron Mass Fraction ---------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Xn, Xns_T, OS_Xn, Units_V = UnitXn )

    ! --- Interpolate Alpha Mass Fraction -----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Xa, Xas_T, OS_Xa, Units_V = UnitXa )

    ! --- Interpolate Heavy Mass Fraction -----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Xh, Xhs_T, OS_Xh, Units_V = UnitXh )

    ! --- Gamma1 ------------------------------------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Gm, Gms_T, OS_Gm, Units_V = UnitGm )

  END SUBROUTINE ApplyEquationOfState_TABLE_Scalar


  SUBROUTINE ApplyEquationOfState_TABLE_Vector &
    ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

    REAL(DP), INTENT(in)  :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out) :: P(:), S(:), E(:), Me(:), Mp(:), Mn(:)
    REAL(DP), INTENT(out) :: Xp(:), Xn(:), Xa(:), Xh(:), Gm(:)

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ApplyEquationOfState_TABLE &
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
    INTEGER :: Error

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / ( Gram / Centimeter**3 )
    E_P = E / ( Erg / Gram )
    Y_P = Y

    IF ( PRESENT( Guess_Option ) ) THEN

      T_Guess = Guess_Option / Kelvin

      CALL ComputeTemperatureWith_DEY &
             ( D_P, E_P, Y_P, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_Lookup, T_Guess, &
               Error_Option = Error )

    ELSE

      CALL ComputeTemperatureWith_DEY &
             ( D_P, E_P, Y_P, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_Lookup, &
               Error_Option = Error )

    END IF
    T = T_Lookup * Kelvin

#endif

    IF ( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE_Vector &
    ( D, E, Y, T, Guess_Option, Error_Option )

    REAL(DP), INTENT(in )           :: D(:), E(:), Y(:)
    REAL(DP), INTENT(out)           :: T(:)
    REAL(DP), INTENT(in ), OPTIONAL :: Guess_Option(:)
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option(:)

    LOGICAL  :: do_gpu
    INTEGER  :: iP, nP
    INTEGER  :: Error(SIZE(D))
    REAL(DP) :: D_P, E_P, Y_P, T_Lookup, T_Guess

#ifdef MICROPHYSICS_WEAKLIB

    do_gpu = QueryOnGPU( D, E, Y, T )
#if defined(THORNADO_DEBUG_EOS) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeTemperatureFromSpecificInternalEnergy_TABLE] Data not present on device'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeTemperatureFromSpecificInternalEnergy_TABLE]   D    missing'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeTemperatureFromSpecificInternalEnergy_TABLE]   E    missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeTemperatureFromSpecificInternalEnergy_TABLE]   Y    missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeTemperatureFromSpecificInternalEnergy_TABLE]   T    missing'
    END IF
#endif

    nP = SIZE( D )

    IF( PRESENT( Guess_Option ) )THEN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( D_P, E_P, Y_P, T_Lookup, T_Guess ) &
    !$OMP MAP( from: Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( D_P, E_P, Y_P, T_Lookup, T_Guess ) &
    !$ACC PRESENT( D, E, Y, T ) &
    !$ACC COPYOUT( Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, E_P, Y_P, T_Lookup, T_Guess )
#endif
      DO iP = 1, nP

        D_P = D(iP) / ( Gram / Centimeter**3 )
        E_P = E(iP) / ( Erg / Gram )
        Y_P = Y(iP)

        T_Guess = Guess_Option(iP) / Kelvin

        CALL ComputeTemperatureWith_DEY &
               ( D_P, E_P, Y_P, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_Lookup, T_Guess, &
                 Error_Option = Error(iP) )

        T(iP) = T_Lookup * Kelvin

      END DO

    ELSE

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( D_P, E_P, Y_P, T_Lookup ) &
    !$OMP MAP( from: Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( D_P, E_P, Y_P, T_Lookup ) &
    !$ACC PRESENT( D, E, Y, T ) &
    !$ACC COPYOUT( Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE( D_P, E_P, Y_P, T_Lookup )
#endif
      DO iP = 1, nP

        D_P = D(iP) / ( Gram / Centimeter**3 )
        E_P = E(iP) / ( Erg / Gram )
        Y_P = Y(iP)

        CALL ComputeTemperatureWith_DEY &
               ( D_P, E_P, Y_P, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_Lookup, &
                 Error_Option = Error(iP) )

        T(iP) = T_Lookup * Kelvin

      END DO

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

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D, Em, Y, T )

    CALL ComputePressure_TABLE &
           ( D, T, Y, P )

  END SUBROUTINE ComputePressureFromPrimitive_TABLE_Scalar


  SUBROUTINE ComputePressureFromPrimitive_TABLE_Vector( D, Ev, Ne, P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: P

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputePressureFromPrimitive_TABLE( D(iP), Ev(iP), Ne(iP), P(iP) )

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

    CALL ComputeTemperatureWith_DEY &
           ( D_P, E_P, Y_P, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_P )

    T = T_P * UnitT

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, P, Ps_T, OS_P, Units_V = UnitP )

#endif

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE_Vector &
    ( D, Em, Y, P )

    REAL(DP), INTENT(in)  :: D(:), Em(:), Y(:)
    REAL(DP), INTENT(out) :: P(:)

    INTEGER iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputePressureFromSpecificInternalEnergy_TABLE &
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

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D, Em, Y, T )

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, P, Ps_T, OS_P, Units_V = Dyne / Centimeter**2 )

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Gm, Gms_T, OS_Gm, Units_V = 1.0_DP )

    Cs = SQRT( Gm * P / D )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Scalar


  SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE_Vector( D, Ev, Ne, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: Cs(:)

    REAL(DP), DIMENSION(SIZE(D)) :: P, T, Y, Em, Gm

    Em = Ev / D
    Y  = Ne * ( BaryonMass / D )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D, Em, Y, T )

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, P, Ps_T, OS_P, Units_V = Dyne / Centimeter**2 )

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Gm, Gms_T, OS_Gm, Units_V = 1.0_DP )

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

    CALL ComputeTemperatureWith_DPY &
           ( D_P, P_P, Y_P, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, &
             Error_Option = Error )

    T = T_P * UnitT

#endif

  END SUBROUTINE ComputeTemperatureFromPressure_TABLE_Scalar


  SUBROUTINE ComputeTemperatureFromPressure_TABLE_Vector( D, P, Y, T )

    REAL(DP), INTENT(in)  :: D(:), P(:), Y(:)
    REAL(DP), INTENT(out) :: T(:)

    INTEGER :: iP, nP

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputeTemperatureFromPressure_TABLE &
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

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Em, Es_T, OS_E, Units_V = UnitE )

      Ev = Em * D              ! --- Internal Energy per Unit Volume
      Ne = Y  * D / BaryonMass ! --- Electrons per Unit Volume

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Scalar


  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE_Vector &
    ( D, T, Y, Ev, Em, Ne )

    REAL(DP), INTENT(in)  :: D (:), T (:), Y (:)
    REAL(DP), INTENT(out) :: Ev(:), Em(:), Ne(:)

    INTEGER :: iP, nP
    LOGICAL :: do_gpu

    nP = SIZE( D )

    do_gpu = QueryOnGPU( D, T, Y, Ev, Em, Ne )
#if defined(THORNADO_DEBUG_EOS) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeThermodynamicStates_Primitive_TABLE] Data not present on device'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Primitive_TABLE]   D    missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Primitive_TABLE]   T    missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Primitive_TABLE]   Y    missing'
      IF ( .not. QueryOnGPU( Ev ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Primitive_TABLE]   Ev   missing'
      IF ( .not. QueryOnGPU( Em ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Primitive_TABLE]   Em   missing'
      IF ( .not. QueryOnGPU( Ne ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Primitive_TABLE]   Ne   missing'
    END IF
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, T, Y, Em, Ev, Ne )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeThermodynamicStates_Primitive_TABLE &
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

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D, Em, Y, T )

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Scalar


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE_Vector &
    ( D, Ev, Ne, T, Em, Y )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: T, Em, Y

    INTEGER :: iP, nP, Error(SIZE(D))
    LOGICAL :: do_gpu

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

    do_gpu = QueryOnGPU( D, Ev, Ne, T, Em, Y )
#if defined(THORNADO_DEBUG_EOS) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeThermodynamicStates_Auxiliary_TABLE] Data not present on device'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Auxiliary_TABLE]   D    missing'
      IF ( .not. QueryOnGPU( Ev ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Auxiliary_TABLE]   Ev   missing'
      IF ( .not. QueryOnGPU( Ne ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Auxiliary_TABLE]   Ne   missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Auxiliary_TABLE]   T    missing'
      IF ( .not. QueryOnGPU( Em ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Auxiliary_TABLE]   Em   missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeThermodynamicStates_Auxiliary_TABLE]   Y    missing'
    END IF
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( from: Error )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, Ev, Ne, T, Em, Y ) &
    !$ACC COPYOUT( Error )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO &
    !$OMP MAP( from: Error )
#endif
    DO iP = 1, nP

      Em(iP) = Ev(iP) / D(iP)              ! --- Internal Energy per Mass
      Y (iP) = Ne(iP) / D(iP) * BaryonMass ! --- Electron Fraction

      CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
             ( D(iP), Em(iP), Y(iP), T(iP), Error_Option = Error(iP) )

    END DO

    IF ( ANY( Error(:) > 0 ) ) THEN
      DO iP = 1, nP
        IF ( Error(iP) > 0 ) CALL DescribeEOSInversionError( Error(iP) )
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

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D, Em, Y, T )

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, P, Ps_T, OS_P, Units_V = UnitP )

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, S, Ss_T, OS_S, Units_V = UnitS )

    CALL ComputeDependentVariable_TABLE &
           ( D, T, Y, Gm, Gms_T, OS_Gm, Units_V = UnitGm )

    Cs = SQRT( Gm * P / D )

  END SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Scalar


  SUBROUTINE ComputeAuxiliary_Fluid_TABLE_Vector &
    ( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

    REAL(DP), INTENT(in)  :: D(:), Ev(:), Ne(:)
    REAL(DP), INTENT(out) :: P(:), T (:), Y (:), S(:), Em(:), Gm(:), Cs(:)

    INTEGER :: iP, nP

    nP = SIZE( D )

    DO iP = 1, nP

      CALL ComputeAuxiliary_Fluid_TABLE &
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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, P, dPdD, dPdT, dPdY, Ps_T, OS_P, Units_V = UnitP )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, P, Ps_T, OS_P, Units_V = UnitP )

    END IF

  END SUBROUTINE ComputePressure_TABLE_Scalar


  SUBROUTINE ComputePressure_TABLE_Vector &
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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, P, dPdD, dPdT, dPdY, Ps_T, OS_P, Units_V = UnitP )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, P, Ps_T, OS_P, Units_V = UnitP )

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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, E, dEdD, dEdT, dEdY, Es_T, OS_E, Units_V = UnitE )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, E, Es_T, OS_E, Units_V = UnitE )

    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Scalar


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE_Vector &
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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, E, dEdD, dEdT, dEdY, Es_T, OS_E, Units_V = UnitE )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, E, Es_T, OS_E, Units_V = UnitE )

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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mes_T, OS_Me, Units_V = UnitMe )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mes_T, OS_Me, Units_V = UnitMe )

    END IF

  END SUBROUTINE ComputeElectronChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeElectronChemicalPotential_TABLE_Vector &
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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mes_T, OS_Me, Units_V = UnitMe )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mes_T, OS_Me, Units_V = UnitMe )

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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mps_T, Os_Mp, Units_V = UnitMp )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mps_T, Os_Mp, Units_V = UnitMp )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeProtonChemicalPotential_TABLE_Vector &
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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mps_T, OS_Mp, Units_V = UnitMp )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mps_T, OS_Mp, Units_V = UnitMp )

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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mns_T, OS_Mn, Units_V = UnitMn )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mns_T, OS_Mn, Units_V = UnitMn )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Scalar


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Vector &
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

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mns_T, OS_Mn, Units_V = UnitMn )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mns_T, OS_Mn, Units_V = UnitMn )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE_Vector


  SUBROUTINE ComputeDependentVariable_TABLE_Scalar &
    ( D, T, Y, V, V_T, OS_V, Units_V )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: V
    REAL(DP), INTENT(in)  :: V_T(:,:,:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Y_P, V_P

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Y_P = Y / UnitY

    CALL LogInterpolateSingleVariable &
           ( D_P, T_P, Y_P, Ds_T, Ts_T, Ys_T, OS_V, V_T, V_P )

    V = V_P * Units_V

#else

    V = Zero

#endif

  END SUBROUTINE ComputeDependentVariable_TABLE_Scalar


  SUBROUTINE ComputeDependentVariable_TABLE_Vector &
    ( D, T, Y, V, V_T, OS_V, Units_V )

    REAL(DP), INTENT(in)  :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out) :: V(:)
    REAL(DP), INTENT(in)  :: V_T(:,:,:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    INTEGER :: iP, nP
    LOGICAL :: do_gpu

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

    do_gpu = QueryOnGPU( D, T, Y, V )
#if defined(THORNADO_DEBUG_EOS) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeDependentVariable_TABLE] Data not present on device'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeDependentVariable_TABLE]   D    missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeDependentVariable_TABLE]   T    missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeDependentVariable_TABLE]   Y    missing'
      IF ( .not. QueryOnGPU( V ) ) &
        WRITE(*,*) '[ComputeDependentVariable_TABLE]   V    missing'
    END IF
#endif

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC IF( do_gpu ) &
      !$ACC PRESENT( D, T, Y, V, OS_V, V_T )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeDependentVariable_TABLE &
             ( D(iP), T(iP), Y(iP), V(iP), V_T, OS_V, Units_V )

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
    REAL(DP), INTENT(in)  :: V_T(:,:,:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Y_P, V_P, dV_P(3)

#ifdef MICROPHYSICS_WEAKLIB

    D_P = D / UnitD
    T_P = T / UnitT
    Y_P = Y / UnitY

    CALL LogInterpolateDifferentiateSingleVariable &
           ( D_P, T_P, Y_P, Ds_T, Ts_T, Ys_T, OS_V, V_T, V_P, dV_P )

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

    REAL(DP), INTENT(in)  :: D(:), T(:), Y(:)
    REAL(DP), INTENT(out) :: V(:), dVdD(:), dVdT(:), dVdY(:)
    REAL(DP), INTENT(in)  :: V_T(:,:,:)
    REAL(DP), INTENT(in)  :: OS_V, Units_V

    INTEGER :: iP, nP
    LOGICAL :: do_gpu

#ifdef MICROPHYSICS_WEAKLIB

    nP = SIZE( D )

    do_gpu = QueryOnGPU( D, T, Y, V, dVdD, dVdT, dVdY )
#if defined(THORNADO_DEBUG_EOS) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeDependentVariableAndDerivatives_TABLE] Data not present on device'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeDependentVariableAndDerivatives_TABLE]   D    missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeDependentVariableAndDerivatives_TABLE]   T    missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeDependentVariableAndDerivatives_TABLE]   Y    missing'
      IF ( .not. QueryOnGPU( V ) ) &
        WRITE(*,*) '[ComputeDependentVariableAndDerivatives_TABLE]   V    missing'
      IF ( .not. QueryOnGPU( dVdD ) ) &
        WRITE(*,*) '[ComputeDependentVariableAndDerivatives_TABLE]   dVdD missing'
      IF ( .not. QueryOnGPU( dVdT ) ) &
        WRITE(*,*) '[ComputeDependentVariableAndDerivatives_TABLE]   dVdT missing'
      IF ( .not. QueryOnGPU( dVdY ) ) &
        WRITE(*,*) '[ComputeDependentVariableAndDerivatives_TABLE]   dVdY missing'
    END IF
#endif

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR &
      !$ACC IF( do_gpu ) &
      !$ACC PRESENT( D, T, Y, V, dVdD, dVdT, dVdY, OS_V, V_T )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO
#endif
    DO iP = 1, nP

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D(iP), T(iP), Y(iP), V(iP), dVdD(iP), dVdT(iP), dVdY(iP), &
               V_T, OS_V, Units_V )

    END DO

#else

    V    = Zero
    dVdD = Zero
    dVdT = Zero
    dVdY = Zero

#endif

  END SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE_Vector


END MODULE EquationOfStateModule_TABLE
