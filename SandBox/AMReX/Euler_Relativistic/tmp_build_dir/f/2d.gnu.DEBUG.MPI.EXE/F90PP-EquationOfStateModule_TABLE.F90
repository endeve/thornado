/Users/dunhamsj/Research/SN/thornado/Modules/EquationOfState/EquationOfStateModule_TABLE.F90:174:13: warning: empty character constant [-Winvalid-pp-token]
            '', 'Table Name: ', TRIM( EquationOfStateTableName )
            ^
/Users/dunhamsj/Research/SN/thornado/Modules/EquationOfState/EquationOfStateModule_TABLE.F90:531:7: warning: empty character constant [-Winvalid-pp-token]
      '', 'ComputePressureFromSpecificInternalEnergy_TABLE Not Implemented'
      ^
2 warnings generated.




MODULE EquationOfStateModule_TABLE



  USE, INTRINSIC :: ISO_C_BINDING
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
  USE DeviceModule, ONLY: &
    QueryOnGpu
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
  PUBLIC :: ComputeSpecificInternalEnergy_Point_TABLE
  PUBLIC :: ComputeSpecificInternalEnergy_Points_TABLE
  PUBLIC :: ComputeElectronChemicalPotential_TABLE
  PUBLIC :: ComputeElectronChemicalPotentialPoints_TABLE
  PUBLIC :: ComputeElectronChemicalPotentialPoint_TABLE
  PUBLIC :: ComputeProtonChemicalPotential_TABLE
  PUBLIC :: ComputeProtonChemicalPotentialPoints_TABLE
  PUBLIC :: ComputeProtonChemicalPotentialPoint_TABLE
  PUBLIC :: ComputeNeutronChemicalPotential_TABLE
  PUBLIC :: ComputeNeutronChemicalPotentialPoints_TABLE
  PUBLIC :: ComputeNeutronChemicalPotentialPoint_TABLE

  INTERFACE ComputeSpecificInternalEnergy_TABLE
    MODULE PROCEDURE ComputeSpecificInternalEnergy_Point_TABLE
    MODULE PROCEDURE ComputeSpecificInternalEnergy_Points_TABLE
  END INTERFACE ComputeSpecificInternalEnergy_TABLE

  INTERFACE ComputeElectronChemicalPotential_TABLE
    MODULE PROCEDURE ComputeElectronChemicalPotentialPoints_TABLE
    MODULE PROCEDURE ComputeElectronChemicalPotentialPoint_TABLE
  END INTERFACE

  INTERFACE ComputeProtonChemicalPotential_TABLE
    MODULE PROCEDURE ComputeProtonChemicalPotentialPoints_TABLE
    MODULE PROCEDURE ComputeProtonChemicalPotentialPoint_TABLE
  END INTERFACE

  INTERFACE ComputeNeutronChemicalPotential_TABLE
    MODULE PROCEDURE ComputeNeutronChemicalPotentialPoints_TABLE
    MODULE PROCEDURE ComputeNeutronChemicalPotentialPoint_TABLE
  END INTERFACE



CONTAINS


  SUBROUTINE InitializeEquationOfState_TABLE &
    ( EquationOfStateTableName_Option, Verbose_Option, External_EOS )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option




    INTEGER,          INTENT(in), OPTIONAL :: External_EOS


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



  END SUBROUTINE InitializeEquationOfState_TABLE


  SUBROUTINE FinalizeEquationOfState_TABLE



  END SUBROUTINE FinalizeEquationOfState_TABLE


  SUBROUTINE ApplyEquationOfState_TABLE &
               ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: P, S, E, Me, Mp, Mn
    REAL(DP), DIMENSION(:), INTENT(out) :: Xp, Xn, Xa, Xh, Gm

    REAL(DP), DIMENSION(SIZE( D )) :: TMP

    ! --- Interpolate Pressure ----------------------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), P(:), Ps_T, OS_P, &
             Units_V = Dyne / Centimeter**2 )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), S(:), Ss_T, OS_S, &
             Units_V = BoltzmannConstant )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), E(:), Es_T, OS_E, &
             Units_V = Erg / Gram )

    ! --- Interpolate Electron Chemical Potential ---------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Me(:), Mes_T, OS_Me, &
             Units_V = MeV )

    ! --- Interpolate Proton Chemical Potential -----------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Mp(:), Mps_T, OS_Mp, &
             Units_V = MeV )

    ! --- Interpolate Neutron Chemical Potential ----------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Mn(:), Mns_T, OS_Mn, &
             Units_V = MeV )

    ! --- Interpolate Proton Mass Fraction ----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Xp(:), Xps_T, OS_Xp, &
             Units_V = 1.0_DP )

    ! --- Interpolate Neutron Mass Fraction ---------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Xn(:), Xns_T, OS_Xn, &
             Units_V = 1.0_DP )

    ! --- Interpolate Alpha Mass Fraction -----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Xa(:), Xas_T, OS_Xa, &
             Units_V = 1.0_DP )

    ! --- Interpolate Heavy Mass Fraction -----------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Xh(:), Xhs_T, OS_Xh, &
             Units_V = 1.0_DP )

    ! --- Gamma1 ------------------------------------------------------

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Gm(:), Gms_T, OS_Gm, &
             Units_V = 1.0_DP )

  END SUBROUTINE ApplyEquationOfState_TABLE


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE( D, E, Y, T )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, E, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: T

    INTEGER :: iP, nP, Error(SIZE(D))
    LOGICAL :: do_gpu



  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergy_TABLE


  SUBROUTINE ComputeTemperatureFromSpecificInternalEnergyPoint_TABLE( D, E, Y, T, Error_Option )






    REAL(DP), INTENT(in)            :: D, E, Y
    REAL(DP), INTENT(out)           :: T
    INTEGER,  INTENT(out), OPTIONAL :: Error_Option

    REAL(DP) :: D_P, E_P, Y_P, T_Lookup
    INTEGER :: Error



    IF ( PRESENT( Error_Option ) ) Error_Option = Error

  END SUBROUTINE ComputeTemperatureFromSpecificInternalEnergyPoint_TABLE


  SUBROUTINE ComputePressureFromPrimitive_TABLE( D, Ev, Ne, P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: P

    REAL(DP), DIMENSION(SIZE(D)) :: Em, T, Y

    Em(:) = Ev(:) / D(:)              ! --- Internal Energy per Mass
    Y (:) = Ne(:) / D(:) * BaryonMass ! --- Electron Fraction

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D(:), Em(:), Y(:), T(:) )

    CALL ComputePressure_TABLE &
           ( D(:), T(:), Y(:), P(:) )

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

    REAL(DP), DIMENSION(SIZE(D)) :: P, T, Y, Em, Gm

    Em(:) = Ev(:) / D(:)
    Y (:) = Ne(:) * ( BaryonMass / D(:) )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D(:), Em(:), Y(:), T(:) )

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), P(:), Ps_T, OS_P, &
             Units_V = Dyne / Centimeter**2 )

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Gm(:), Gms_T, OS_Gm, &
             Units_V = 1.0_DP )

    Cs(:) = SQRT( Gm(:) * P(:) / D(:) )

  END SUBROUTINE ComputeSoundSpeedFromPrimitive_TABLE


  SUBROUTINE ComputeTemperatureFromPressure_TABLE( D, P, Y, T )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, P, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: T

    REAL(DP), DIMENSION(SIZE(T)) :: T_Bisection



  END SUBROUTINE ComputeTemperatureFromPressure_TABLE


  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE( D, T, Y, Ev, Em, Ne )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Ev, Em, Ne

    INTEGER :: iP, nP
    LOGICAL :: do_gpu

    nP = SIZE(D)

    do_gpu = QueryOnGPU( D, T, Y, Ev, Em, Ne )










    DO iP = 1, nP

      ! --- Interpolate Specific Internal Energy ------------------------

      CALL ComputeDependentVariablePoint_TABLE &
             ( D(iP), T(iP), Y(iP), Em(iP), Es_T, OS_E, &
               Units_V = Erg / Gram )

      Ev(iP) = Em(iP) * D(iP)              ! --- Internal Energy per Unit Volume
      Ne(iP) = Y (iP) * D(iP) / BaryonMass ! --- Electrons per Unit Volume

    END DO

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE( D, Ev, Ne, T, Em, Y )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: T, Em, Y

    INTEGER :: iP, nP, Error(SIZE(D))
    LOGICAL :: do_gpu



  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE


  SUBROUTINE ComputeAuxiliary_Fluid_TABLE( D, Ev, Ne, P, T, Y, S, Em, Gm, Cs )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: P, T, Y, S, Em, Gm, Cs

    Em(:) = Ev(:) / D(:)
    Y (:) = Ne(:) * ( BaryonMass / D(:) )

    CALL ComputeTemperatureFromSpecificInternalEnergy_TABLE &
           ( D(:), Em(:), Y(:), T(:) )

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), P(:), Ps_T, OS_P, &
             Units_V = Dyne / Centimeter**2 )

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), S(:), Ss_T, OS_S, &
             Units_V = BoltzmannConstant )

    CALL ComputeDependentVariable_TABLE &
           ( D(:), T(:), Y(:), Gm(:), Gms_T, OS_Gm, &
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
             [ Auxiliary_Fluid_TABLE(iAF_Ye) ], TMP, Ps_T, OS_P, &
             Units_V = Dyne / Centimeter**2 )

    Auxiliary_Fluid_TABLE(iAF_P) = TMP(1)

    CALL ComputeDependentVariable_TABLE &
           ( [ PF(iPF_D) ], [ Auxiliary_Fluid_TABLE(iAF_T) ], &
             [ Auxiliary_Fluid_TABLE(iAF_Ye) ], TMP, Ss_T, OS_S, &
             Units_V = BoltzmannConstant )

    Auxiliary_Fluid_TABLE(iAF_S) = TMP(1)

    CALL ComputeDependentVariable_TABLE &
           ( [ PF(iPF_D) ], [ Auxiliary_Fluid_TABLE(iAF_T) ], &
             [ Auxiliary_Fluid_TABLE(iAF_Ye) ], TMP, Gms_T, OS_Gm, &
             Units_V = 1.0_DP )

    Auxiliary_Fluid_TABLE(iAF_Gm) = TMP(1)

    Auxiliary_Fluid_TABLE(iAF_Cs) &
      = SQRT( Auxiliary_Fluid_TABLE(iAF_Gm) &
              * Auxiliary_Fluid_TABLE(iAF_P) / PF(iPF_D) )

    RETURN
  END FUNCTION Auxiliary_Fluid_TABLE


  SUBROUTINE ComputePressure_TABLE &
               ( D, T, Y, P, dPdD_Option, dPdT_Option, dPdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)                    :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)                   :: P
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dPdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dPdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dPdY_Option

    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: dPdD_Local, dPdT_Local, dPdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: dPdD, dPdT, dPdY
    INTEGER :: nP
    LOGICAL :: ComputeDerivatives

    nP = SIZE(D)
    ComputeDerivatives = .FALSE.

    IF ( PRESENT( dPdD_Option ) ) THEN
      dPdD(1:nP) => dPdD_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dPdD(1:nP) => dPdD_Local(:)
    END IF

    IF ( PRESENT( dPdT_Option ) ) THEN
      dPdT(1:nP) => dPdT_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dPdT(1:nP) => dPdT_Local(:)
    END IF

    IF ( PRESENT( dPdY_Option ) ) THEN
      dPdY(1:nP) => dPdY_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dPdY(1:nP) => dPdY_Local(:)
    END IF

    IF ( ComputeDerivatives ) THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, P, dPdD, dPdT, dPdY, Ps_T, OS_P, &
               Units_V = Dyne / Centimeter**2 )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, P, Ps_T, OS_P, &
               Units_V = Dyne / Centimeter**2 )

    END IF

  END SUBROUTINE ComputePressure_TABLE


  SUBROUTINE ComputeSpecificInternalEnergy_Point_TABLE &
    ( D, T, Y, E, dEdD_Option, dEdT_Option, dEdY_Option )







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

    IF ( PRESENT( dEdD_Option ) ) THEN
      dEdD => dEdD_Option
    ELSE
      dEdD => dEdD_Local
    END IF

    IF ( PRESENT( dEdT_Option ) ) THEN
      dEdT => dEdT_Option
    ELSE
      dEdT => dEdT_Local
    END IF

    IF ( PRESENT( dEdY_Option ) ) THEN
      dEdY => dEdY_Option
    ELSE
      dEdY => dEdY_Local
    END IF

    IF( ComputeDerivatives )THEN

      CALL ComputeDependentVariableAndDerivativesPoint_TABLE &
             ( D, T, Y, E, dEdD, dEdT, dEdY, Es_T, OS_E, Units_V = Erg / Gram )

    ELSE

      CALL ComputeDependentVariablePoint_TABLE &
             ( D, T, Y, E, Es_T, OS_E, Units_V = Erg / Gram )

    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_Point_TABLE


  SUBROUTINE ComputeSpecificInternalEnergy_Points_TABLE &
    ( D, T, Y, E, dEdD_Option, dEdT_Option, dEdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)                    :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)                   :: E
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dEdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dEdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dEdY_Option

    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: dEdD_Local, dEdT_Local, dEdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: dEdD, dEdT, dEdY
    INTEGER :: nP
    LOGICAL :: ComputeDerivatives

    nP = SIZE(D)
    ComputeDerivatives = .FALSE.

    IF ( PRESENT( dEdD_Option ) ) THEN
      dEdD(1:nP) => dEdD_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dEdD(1:nP) => dEdD_Local(:)
    END IF

    IF ( PRESENT( dEdT_Option ) ) THEN
      dEdT(1:nP) => dEdT_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dEdT(1:nP) => dEdT_Local(:)
    END IF

    IF ( PRESENT( dEdY_Option ) ) THEN
      dEdY(1:nP) => dEdY_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dEdY(1:nP) => dEdY_Local(:)
    END IF

    IF ( ComputeDerivatives ) THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, E, dEdD, dEdT, dEdY, Es_T, OS_E, &
               Units_V = Erg / Gram )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, E, Es_T, OS_E, &
               Units_V = Erg / Gram )

    END IF

  END SUBROUTINE ComputeSpecificInternalEnergy_Points_TABLE


  SUBROUTINE ComputeElectronChemicalPotentialPoints_TABLE &
    ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)                    :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)                   :: M
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: dMdD, dMdT, dMdY
    INTEGER :: nP
    LOGICAL :: ComputeDerivatives

    nP = SIZE(D)
    ComputeDerivatives = .FALSE.

    IF ( PRESENT( dMdD_Option ) ) THEN
      dMdD(1:nP) => dMdD_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdD(1:nP) => dMdD_Local(:)
    END IF

    IF ( PRESENT( dMdT_Option ) ) THEN
      dMdT(1:nP) => dMdT_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdT(1:nP) => dMdT_Local(:)
    END IF

    IF ( PRESENT( dMdY_Option ) ) THEN
      dMdY(1:nP) => dMdY_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdY(1:nP) => dMdY_Local(:)
    END IF

    IF ( ComputeDerivatives ) THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mes_T, OS_Me, &
               Units_V = MeV )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mes_T, OS_Me, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeElectronChemicalPotentialPoints_TABLE


  SUBROUTINE ComputeElectronChemicalPotentialPoint_TABLE &
               ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )






    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD, dMdT, dMdY
    LOGICAL :: ComputeDerivatives

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF ( PRESENT( dMdD_Option ) ) THEN
      dMdD => dMdD_Option
    ELSE
      dMdD => dMdD_Local
    END IF

    IF ( PRESENT( dMdT_Option ) ) THEN
      dMdT => dMdT_Option
    ELSE
      dMdT => dMdT_Local
    END IF

    IF ( PRESENT( dMdY_Option ) ) THEN
      dMdY => dMdY_Option
    ELSE
      dMdY => dMdY_Local
    END IF

    IF ( ComputeDerivatives ) THEN

      CALL ComputeDependentVariableAndDerivativesPoint_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mes_T, OS_Me, &
               Units_V = MeV )

    ELSE

      CALL ComputeDependentVariablePoint_TABLE &
             ( D, T, Y, M, Mes_T, OS_Me, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeElectronChemicalPotentialPoint_TABLE


  SUBROUTINE ComputeProtonChemicalPotentialPoints_TABLE &
               ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)                    :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)                   :: M
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: dMdD, dMdT, dMdY
    INTEGER :: nP
    LOGICAL :: ComputeDerivatives

    nP = SIZE(D)
    ComputeDerivatives = .FALSE.

    IF ( PRESENT( dMdD_Option ) ) THEN
      dMdD(1:nP) => dMdD_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdD(1:nP) => dMdD_Local(:)
    END IF

    IF ( PRESENT( dMdT_Option ) ) THEN
      dMdT(1:nP) => dMdT_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdT(1:nP) => dMdT_Local(:)
    END IF

    IF ( PRESENT( dMdY_Option ) ) THEN
      dMdY(1:nP) => dMdY_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdY(1:nP) => dMdY_Local(:)
    END IF

    IF ( ComputeDerivatives ) THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mps_T, OS_Mp, &
               Units_V = MeV )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mps_T, OS_Mp, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotentialPoints_TABLE


  SUBROUTINE ComputeProtonChemicalPotentialPoint_TABLE &
               ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )






    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD, dMdT, dMdY
    LOGICAL :: ComputeDerivatives

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF ( PRESENT( dMdD_Option ) ) THEN
      dMdD => dMdD_Option
    ELSE
      dMdD => dMdD_Local
    END IF

    IF ( PRESENT( dMdT_Option ) ) THEN
      dMdT => dMdT_Option
    ELSE
      dMdT => dMdT_Local
    END IF

    IF ( PRESENT( dMdY_Option ) ) THEN
      dMdY => dMdY_Option
    ELSE
      dMdY => dMdY_Local
    END IF

    IF ( ComputeDerivatives ) THEN

      CALL ComputeDependentVariableAndDerivativesPoint_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mps_T, Os_Mp, &
               Units_V = MeV )

    ELSE

      CALL ComputeDependentVariablePoint_TABLE &
             ( D, T, Y, M, Mps_T, Os_Mp, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeProtonChemicalPotentialPoint_TABLE


  SUBROUTINE ComputeNeutronChemicalPotentialPoints_TABLE &
               ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)                    :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)                   :: M
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    REAL(DP), DIMENSION(SIZE(D)), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), DIMENSION(:)      , POINTER :: dMdD, dMdT, dMdY
    INTEGER :: nP
    LOGICAL :: ComputeDerivatives

    nP = SIZE(D)
    ComputeDerivatives = .FALSE.

    IF ( PRESENT( dMdD_Option ) ) THEN
      dMdD(1:nP) => dMdD_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdD(1:nP) => dMdD_Local(:)
    END IF

    IF ( PRESENT( dMdT_Option ) ) THEN
      dMdT(1:nP) => dMdT_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdT(1:nP) => dMdT_Local(:)
    END IF

    IF ( PRESENT( dMdY_Option ) ) THEN
      dMdY(1:nP) => dMdY_Option(:)
      ComputeDerivatives = .TRUE.
    ELSE
      dMdY(1:nP) => dMdY_Local(:)
    END IF

    IF ( ComputeDerivatives ) THEN

      CALL ComputeDependentVariableAndDerivatives_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mns_T, OS_Mn, &
               Units_V = MeV )

    ELSE

      CALL ComputeDependentVariable_TABLE &
             ( D, T, Y, M, Mns_T, OS_Mn, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotentialPoints_TABLE


  SUBROUTINE ComputeNeutronChemicalPotentialPoint_TABLE &
               ( D, T, Y, M, dMdD_Option, dMdT_Option, dMdY_Option )






    REAL(DP), INTENT(in)                    :: D, T, Y
    REAL(DP), INTENT(out)                   :: M
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdD_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdT_Option
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: dMdY_Option

    REAL(DP), TARGET  :: dMdD_Local, dMdT_Local, dMdY_Local
    REAL(DP), POINTER :: dMdD, dMdT, dMdY
    LOGICAL :: ComputeDerivatives

    ComputeDerivatives &
      =      PRESENT( dMdD_Option ) &
        .OR. PRESENT( dMdT_Option ) &
        .OR. PRESENT( dMdY_Option )

    IF ( PRESENT( dMdD_Option ) ) THEN
      dMdD => dMdD_Option
    ELSE
      dMdD => dMdD_Local
    END IF

    IF ( PRESENT( dMdT_Option ) ) THEN
      dMdT => dMdT_Option
    ELSE
      dMdT => dMdT_Local
    END IF

    IF ( PRESENT( dMdY_Option ) ) THEN
      dMdY => dMdY_Option
    ELSE
      dMdY => dMdY_Local
    END IF

    IF ( ComputeDerivatives ) THEN

      CALL ComputeDependentVariableAndDerivativesPoint_TABLE &
             ( D, T, Y, M, dMdD, dMdT, dMdY, Mns_T, OS_Mn, &
               Units_V = MeV )

    ELSE

      CALL ComputeDependentVariablePoint_TABLE &
             ( D, T, Y, M, Mns_T, OS_Mn, &
               Units_V = MeV )

    END IF

  END SUBROUTINE ComputeNeutronChemicalPotentialPoint_TABLE


  SUBROUTINE ComputeDependentVariable_TABLE( D, T, Y, V, V_T, OS_V, Units_V )

    REAL(DP), DIMENSION(:),     INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:),     INTENT(out) :: V
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: V_T
    REAL(DP),                   INTENT(in)  :: OS_V, Units_V

    INTEGER :: iP, nP
    LOGICAL :: do_gpu



  END SUBROUTINE ComputeDependentVariable_TABLE


  SUBROUTINE ComputeDependentVariablePoint_TABLE( D, T, Y, V, V_T, OS_V, Units_V )






    REAL(DP),                   INTENT(in)  :: D, T, Y
    REAL(DP),                   INTENT(out) :: V
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: V_T
    REAL(DP),                   INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Y_P, V_P



    V = 0.0_DP



  END SUBROUTINE ComputeDependentVariablePoint_TABLE


  SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE &
               ( D, T, Y, V, dVdD, dVdT, dVdY, V_T, OS_V, Units_V )

    REAL(DP), DIMENSION(:),     INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:),     INTENT(out) :: V, dVdD, dVdT, dVdY
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: V_T
    REAL(DP),                   INTENT(in)  :: OS_V, Units_V

    INTEGER :: iP, nP
    LOGICAL :: do_gpu



  END SUBROUTINE ComputeDependentVariableAndDerivatives_TABLE


  SUBROUTINE ComputeDependentVariableAndDerivativesPoint_TABLE &
               ( D, T, Y, V, dVdD, dVdT, dVdY, V_T, OS_V, Units_V )






    REAL(DP),                   INTENT(in)  :: D, T, Y
    REAL(DP),                   INTENT(out) :: V, dVdD, dVdT, dVdY
    REAL(DP), DIMENSION(:,:,:), INTENT(in)  :: V_T
    REAL(DP),                   INTENT(in)  :: OS_V, Units_V

    REAL(DP) :: D_P, T_P, Y_P, V_P, dV_P(3)



    V = 0.0_DP
    dVdD = 0.0_DP
    dVdT = 0.0_DP
    dVdY = 0.0_DP



  END SUBROUTINE ComputeDependentVariableAndDerivativesPoint_TABLE


END MODULE EquationOfStateModule_TABLE

