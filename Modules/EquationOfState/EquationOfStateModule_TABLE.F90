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
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateDifferentiateSingleVariable, &
    ComputeTempFromIntEnergy_Lookup, &
    ComputeTempFromIntEnergy_Bisection, &
    ComputeTempFromIntEnergy_Secant, &
    ComputeTempFromPressure, &
    ComputeTempFromPressure_Bisection

  ! ----------------------------------------------

#endif

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    AtomicMassUnit, &
    BoltzmannConstant, &
    Gram, &
    Centimeter, &
    Kelvin, &
    Dyne, &
    Erg, &
    MeV
  USE FluidFieldsModule, ONLY: &
    ! --- Primitive Fluid Fields:
    iPF_D, iPF_E, iPF_Ne, nPF, &
    ! --- Auxiliary Fluid Fields:
    iAF_P, iAF_T, iAF_Ye, iAF_E, iAF_S, iAF_Gm, iAF_Cs, nAF

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
    OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, &
    OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm
  REAL(DP), PARAMETER :: &
    BaryonMass = AtomicMassUnit
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Ds_T, Ts_T, Ys_T
#ifdef MICROPHYSICS_WEAKLIB
  TYPE(EquationOfStateTableType), POINTER :: EOS
  LOGICAL :: UsingExternalEOS
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
  PUBLIC :: Auxiliary_Fluid_TABLE
  PUBLIC :: ComputePressure_TABLE
  PUBLIC :: ComputeSpecificInternalEnergy_TABLE
  PUBLIC :: ComputeElectronChemicalPotential_TABLE
  PUBLIC :: ComputeProtonChemicalPotential_TABLE
  PUBLIC :: ComputeNeutronChemicalPotential_TABLE

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

!    IF( .NOT. PRESENT( External_EOS ) ) THEN

       ALLOCATE( EOS )
       UsingExternalEOS = .false.       

       CALL InitializeHDF( )

       CALL ReadEquationOfStateTableHDF &
            ( EOS, TRIM( EquationOfStateTableName ) )

       CALL FinalizeHDF( )

    ! ELSE

    !    EOS => External_EOS
    !    UsingExternalEOS = .true.

    ! END IF

    ! --- Thermodynamic State Indices ---

    iD_T = EOS % TS % Indices % iRho
    iT_T = EOS % TS % Indices % iT
    iY_T = EOS % TS % Indices % iYe

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

#endif

  END SUBROUTINE InitializeEquationOfState_TABLE


  SUBROUTINE FinalizeEquationOfState_TABLE

    DEALLOCATE( Ds_T, Ts_T, Ys_T )

#ifdef MICROPHYSICS_WEAKLIB

    IF ( UsingExternalEOS ) THEN
       DEALLOCATE( EOS )
    END IF

#endif

  END SUBROUTINE FinalizeEquationOfState_TABLE


  SUBROUTINE ApplyEquationOfState_TABLE &
               ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: P, S, E, Me, Mp, Mn
    REAL(DP), DIMENSION(:), INTENT(out) :: Xp, Xn, Xa, Xh, Gm

    REAL(DP), DIMENSION(SIZE( D )) :: TMP

    ! --- Interpolate Pressure ----------------------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), P(:), iP_T, OS_P, &
             Units_V = Dyne / Centimeter**2 )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), S(:), iS_T, OS_S, &
             Units_V = BoltzmannConstant )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), E(:), iE_T, OS_E, &
             Units_V = Erg / Gram )

    ! --- Interpolate Electron Chemical Potential ---------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Me(:), iMe_T, OS_Me, &
             Units_V = MeV )

    ! --- Interpolate Proton Chemical Potential -----------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Mp(:), iMp_T, OS_Mp, &
             Units_V = MeV )

    ! --- Interpolate Neutron Chemical Potential ----------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Mn(:), iMn_T, OS_Mn, &
             Units_V = MeV )

    ! --- Interpolate Proton Mass Fraction ----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Xp(:), iXp_T, OS_Xp, &
             Units_V = 1.0_DP )

    ! --- Interpolate Neutron Mass Fraction ---------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Xn(:), iXn_T, OS_Xn, &
             Units_V = 1.0_DP )

    ! --- Interpolate Alpha Mass Fraction -----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Xa(:), iXa_T, OS_Xa, &
             Units_V = 1.0_DP )

    ! --- Interpolate Heavy Mass Fraction -----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Xh(:), iXh_T, OS_Xh, &
             Units_V = 1.0_DP )

    ! --- Gamma1 ------------------------------------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Gm(:), iGm_T, OS_Gm, &
             Units_V = 1.0_DP )

  END SUBROUTINE ApplyEquationOfState_TABLE


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE( D, E, Y, T )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, E, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: T

    REAL(DP), DIMENSION(SIZE(T)) :: T_Lookup

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeTempFromIntEnergy_Lookup     &
           ( D / ( Gram / Centimeter**3 ),   &
             E / ( Erg / Gram ),             &
             Y, Ds_T, Ts_T, Ys_T, LogInterp, &
             EOS % DV % Variables(iE_T) % Values, OS_E, &
             T_Lookup )

    T = T_Lookup * Kelvin

#endif

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE


  SUBROUTINE ComputePressureFromPrimitive_TABLE( D, Ev, Ne, P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: P

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'ComputePressureFromPrimitive_TABLE Not Implemented'
    WRITE(*,*)

    STOP

  END SUBROUTINE ComputePressureFromPrimitive_TABLE


  SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE( D, E, Y, P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, E, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: P

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'ComputePressureFromSpecificInternalEnergy_TABLE Not Implemented'
    WRITE(*,*)

    STOP

  END SUBROUTINE ComputePressureFromSpecificInternalEnergy_TABLE


  SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE( D, Ev, Ne, Cs )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: Cs

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'ComputeSoundSpeedFromPrimitive_TABLE Not Implemented'
    WRITE(*,*)

    STOP

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE


  SUBROUTINE ComputeTemperatureFromPressure_TABLE( D, P, Y, T )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, P, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: T

    REAL(DP), DIMENSION(SIZE(T)) :: T_Bisection

#ifdef MICROPHYSICS_WEAKLIB

    CALL ComputeTempFromPressure_Bisection  &
           ( D / ( Gram / Centimeter**3 ),   &
             P / ( Dyne / Centimeter**2 ),             &
             Y, Ds_T, Ts_T, Ys_T, LogInterp, &
             EOS % DV % Variables(iP_T) % Values, OS_P, &
             T_Bisection )

    T = T_Bisection * Kelvin

#endif

  END SUBROUTINE ComputeTemperatureFromPressure_TABLE


  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE( D, T, Y, Ev, Em, Ne )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Ev, Em, Ne

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Em(:), iE_T, OS_E, &
             Units_V = Erg / Gram )

    Ev(:) = Em(:) * D(:)              ! --- Internal Energy per Unit Volume
    Ne(:) = Y (:) * D(:) / BaryonMass ! --- Electrons per Unit Volume

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE( D, Ev, Ne, T, Em, Y )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: T, Em, Y

    Em(:) = Ev(:) / D(:)              ! --- Internal Energy per Mass
    Y (:) = Ne(:) / D(:) * BaryonMass ! --- Electron Fraction

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D(:), Em(:), Y(:), T(:) )

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE


  SUBROUTINE ComputeAuxiliary_Fluid_TABLE( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: P, T, Y, S, Em, Gm, Cs

    Em(:) = Ev(:) / D(:)
    Y (:) = Ne(:) * ( BaryonMass / D(:) )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D(:), Em(:), Y(:), T(:) )

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), P(:), iP_T, OS_P, &
             Units_V = Dyne / Centimeter**2 )

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), S(:), iS_T, OS_S, &
             Units_V = BoltzmannConstant )

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Gm(:), iGm_T, OS_Gm, &
             Units_V = 1.0_DP )

    Cs(:) = SQRT( Gm(:) * P(:) / D(:) )

  END SUBROUTINE ComputeAuxiliary_Fluid_TABLE


  FUNCTION Auxiliary_Fluid_TABLE( PF )

    REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
    REAL(DP), DIMENSION(nAF)             :: Auxiliary_Fluid_TABLE

    REAL(DP), DIMENSION(1) :: TMP

    Auxiliary_Fluid_TABLE(1:nAF) = 0.0_DP

    Auxiliary_Fluid_TABLE(iAF_E) &
      = PF(iPF_E) / PF(iPF_D)

    Auxiliary_Fluid_TABLE(iAF_Ye) &
      = PF(iPF_Ne) * ( BaryonMass / PF(iPF_D) )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE   &
           ( [ PF(iPF_D) ], [ Auxiliary_Fluid_TABLE(iAF_E) ], &
             [ Auxiliary_Fluid_TABLE(iAF_Ye) ], TMP )

    Auxiliary_Fluid_TABLE(iAF_T) = TMP(1)

    CALL ComputeDependentVariable_TABLE &
           ( [ PF(iPF_D) ], [ Auxiliary_Fluid_TABLE(iAF_T) ], &
             [ Auxiliary_Fluid_TABLE(iAF_Ye) ], TMP, iP_T, OS_P, &
             Units_V = Dyne / Centimeter**2 )

    Auxiliary_Fluid_TABLE(iAF_P) = TMP(1)

    CALL ComputeDependentVariable_TABLE &
           ( [ PF(iPF_D) ], [ Auxiliary_Fluid_TABLE(iAF_T) ], &
             [ Auxiliary_Fluid_TABLE(iAF_Ye) ], TMP, iS_T, OS_S, &
             Units_V = BoltzmannConstant )

    Auxiliary_Fluid_TABLE(iAF_S) = TMP(1)

    CALL ComputeDependentVariable_TABLE &
           ( [ PF(iPF_D) ], [ Auxiliary_Fluid_TABLE(iAF_T) ], &
             [ Auxiliary_Fluid_TABLE(iAF_Ye) ], TMP, iGm_T, OS_Gm, &
             Units_V = 1.0_DP )

    Auxiliary_Fluid_TABLE(iAF_Gm) = TMP(1)

    Auxiliary_Fluid_TABLE(iAF_Cs) &
      = SQRT( Auxiliary_Fluid_TABLE(iAF_Gm) &
              * Auxiliary_Fluid_TABLE(iAF_P) / PF(iPF_D) )

    RETURN
  END FUNCTION Auxiliary_Fluid_TABLE


  SUBROUTINE ComputePressure_TABLE & 
               ( D, T, Y, P, dPdD_Option, dPdT_Option, dPdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: P   
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dPdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dPdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dPdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(1:SIZE( D ),1:3) :: TMP

    ComputeDerivatives = .FALSE.
    IF( ANY( [ PRESENT( dPdD_Option ), PRESENT( dPdT_Option ), & 
               PRESENT( dPdY_Option ) ] ) ) ComputeDerivatives = .TRUE.

    IF( ComputeDerivatives )THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D(:), T(:), Y(:), P(:), TMP(:,1:3), iP_T, OS_P, &
               Units_V = Dyne / Centimeter**2 )
      
       IF( PRESENT( dPdD_Option ) ) dPdD_Option(:) = TMP(:,1)
       IF( PRESENT( dPdT_Option ) ) dPdT_Option(:) = TMP(:,2)
       IF( PRESENT( dPdY_Option ) ) dPdY_Option(:) = TMP(:,3) 

    ELSE

      CALL ComputeDependentVariable_TABLE & 
             ( D(:), T(:), Y(:), P(:), iP_T, OS_P, & 
               Units_V = Dyne / Centimeter**2 )

    END IF

  END SUBROUTINE ComputePressure_TABLE


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE &
               ( D, T, Y, E, dEdD_Option, dEdT_Option, dEdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: E
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dEdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dEdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dEdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(1:SIZE( D ),1:3) :: TMP

    ComputeDerivatives = .FALSE.
    IF( ANY( [ PRESENT( dEdD_Option ), PRESENT( dEdT_Option ), &
               PRESENT( dEdY_Option ) ] ) ) ComputeDerivatives = .TRUE.

    IF( ComputeDerivatives )THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D(:), T(:), Y(:), E(:), TMP(:,1:3), iE_T, OS_E, &
               Units_V = Erg / Gram )

      IF( PRESENT( dEdD_Option ) ) dEdD_Option(:) = TMP(:,1)
      IF( PRESENT( dEdT_Option ) ) dEdT_Option(:) = TMP(:,2)
      IF( PRESENT( dEdY_Option ) ) dEdY_Option(:) = TMP(:,3)

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D(:), T(:), Y(:), E(:), iE_T, OS_E, &
               Units_V = Erg / Gram )

    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE


  SUBROUTINE ComputeElectronChemicalPotential_TABLE &
               ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: M
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(1:SIZE( D ),1:3) :: TMP

    ComputeDerivatives = .FALSE.
    IF( ANY( [ PRESENT( dMdD_Option ), PRESENT( dMdT_Option ), &
               PRESENT( dMdY_Option ) ] ) ) ComputeDerivatives = .TRUE.

    IF( ComputeDerivatives )THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D(:), T(:), Y(:), M(:), TMP(:,1:3), iMe_T, OS_Me, &
               Units_V = MeV )

      IF( PRESENT( dMdD_Option ) ) dMdD_Option(:) = TMP(:,1)
      IF( PRESENT( dMdT_Option ) ) dMdT_Option(:) = TMP(:,2)
      IF( PRESENT( dMdY_Option ) ) dMdY_Option(:) = TMP(:,3)

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D(:), T(:), Y(:), M(:), iMe_T, OS_Me, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeElectronChemicalPotential_TABLE


  SUBROUTINE ComputeProtonChemicalPotential_TABLE &
               ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: M
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(1:SIZE( D ),1:3) :: TMP

    ComputeDerivatives = .FALSE.
    IF( ANY( [ PRESENT( dMdD_Option ), PRESENT( dMdT_Option ), &
               PRESENT( dMdY_Option ) ] ) ) ComputeDerivatives = .TRUE.

    IF( ComputeDerivatives )THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D(:), T(:), Y(:), M(:), TMP(:,1:3), iMp_T, OS_Mp, &
               Units_V = MeV )

      IF( PRESENT( dMdD_Option ) ) dMdD_Option(:) = TMP(:,1)
      IF( PRESENT( dMdT_Option ) ) dMdT_Option(:) = TMP(:,2)
      IF( PRESENT( dMdY_Option ) ) dMdY_Option(:) = TMP(:,3)

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D(:), T(:), Y(:), M(:), iMp_T, OS_Mp, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE &
               ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: M
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(1:SIZE( D ),1:3) :: TMP

    ComputeDerivatives = .FALSE.
    IF( ANY( [ PRESENT( dMdD_Option ), PRESENT( dMdT_Option ), &
               PRESENT( dMdY_Option ) ] ) ) ComputeDerivatives = .TRUE.

    IF( ComputeDerivatives )THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D(:), T(:), Y(:), M(:), TMP(:,1:3), iMn_T, OS_Mn, &
               Units_V = MeV )

      IF( PRESENT( dMdD_Option ) ) dMdD_Option(:) = TMP(:,1)
      IF( PRESENT( dMdT_Option ) ) dMdT_Option(:) = TMP(:,2)
      IF( PRESENT( dMdY_Option ) ) dMdY_Option(:) = TMP(:,3)

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D(:), T(:), Y(:), M(:), iMn_T, OS_Mn, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE


  SUBROUTINE ComputeDependentVariable_TABLE( D, T, Y, V, iV, OS_V, Units_V )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: V
    INTEGER,                INTENT(in)  :: iV
    REAL(DP),               INTENT(in)  :: OS_V, Units_V

    REAL(DP), DIMENSION(1:SIZE( D )) :: TMP

#ifdef MICROPHYSICS_WEAKLIB

    CALL LogInterpolateSingleVariable         &
           ( D(:) / ( Gram / Centimeter**3 ), &
             T(:) / Kelvin,                   &
             Y(:), Ds_T, Ts_T, Ys_T,          &
             LogInterp, OS_V,                 &
             EOS % DV % Variables(iV) % Values, TMP )

#endif

    V(:) = TMP(:) * Units_V

  END SUBROUTINE ComputeDependentVariable_TABLE


  SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE &
               ( D, T, Y, V, dV, iV, OS_V, Units_V )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:),   INTENT(out) :: V
    REAL(DP), DIMENSION(:,:), INTENT(out) :: dV
    INTEGER,                  INTENT(in)  :: iV
    REAL(DP),                 INTENT(in)  :: OS_V, Units_V

    REAL(DP), DIMENSION(1:SIZE( D ))     :: TMP
    REAL(DP), DIMENSION(1:SIZE( D ),1:3) :: dTMP

#ifdef MICROPHYSICS_WEAKLIB

    CALL LogInterpolateDifferentiateSingleVariable &
           ( D(:) / ( Gram / Centimeter**3 ),      &
             T(:) / Kelvin,                        &
             Y(:), Ds_T, Ts_T, Ys_T,               &
             LogInterp, OS_V,                      &
             EOS % DV % Variables(iV) % Values,    &
             TMP, dTMP )

#endif

    V(:) = TMP(:) * Units_V

    dV(:,1) = dTMP(:,1) * Units_V / ( Gram / Centimeter**3 )
    dV(:,2) = dTMP(:,2) * Units_V / Kelvin
    dV(:,3) = dTMP(:,3) * Units_V

  END SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE


END MODULE EquationOfStateModule_TABLE
