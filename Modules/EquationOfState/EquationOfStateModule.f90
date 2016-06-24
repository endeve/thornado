MODULE EquationOfStateModule

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
    ComputeTempFromIntEnergy

  ! ----------------------------------------------

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    AtomicMassUnit, &
    Gram, &
    Centimeter, &
    Kelvin, &
    Dyne, &
    BoltzmannConstant, &
    Erg, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE UtilitiesModule, ONLY: &
    MapTo1D, &
    MapFrom1D
  USE FluidFieldsModule, ONLY: &
    uPF, nPF, & ! - Primitive Fluid Fields
    iPF_D, iPF_E, &
    uAF, nAF, & ! - Auxiliary Fluid Fields
    iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Gm, iAF_Cs

  IMPLICIT NONE
  PRIVATE

  CHARACTER(5) :: &
    EquationOfState &
      = 'IDEAL'
  CHARACTER(256) :: &
    EquationOfStateTableName &
      = 'EquationOfStateTable.h5'
  REAL(DP) :: &
    Gamma_IDEAL &
      = 5.0_DP / 3.0_DP
  TYPE(EquationOfStateTableType) :: &
    EOS

  REAL(DP), PUBLIC, PARAMETER :: &
    BaryonMass = AtomicMassUnit

  PROCEDURE ( ), POINTER, PUBLIC :: &
    ApplyEquationOfState => NULL()

  INTERFACE
    SUBROUTINE ComputeThermodynamicStates( IA1, IA2, IA3, OA1, OA2, OA3 )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)  :: IA1, IA2, IA3
      REAL(DP), DIMENSION(:), INTENT(out) :: OA1, OA2, OA3
    END SUBROUTINE ComputeThermodynamicStates
  END INTERFACE

  PROCEDURE (ComputeThermodynamicStates), POINTER, PUBLIC :: &
    ComputeThermodynamicStates_Primitive => NULL(), &
    ComputeThermodynamicStates_Auxiliary => NULL()

  INTERFACE
    SUBROUTINE ComputeChemicalPotentialsSubroutine &
                 ( D, T, Y, Me, Mp, Mn )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
      REAL(DP), DIMENSION(:), INTENT(out) :: Me, Mp, Mn
    END SUBROUTINE ComputeChemicalPotentialsSubroutine
  END INTERFACE

  PROCEDURE (ComputeChemicalPotentialsSubroutine), POINTER, PUBLIC :: &
    ComputeChemicalPotentials => NULL()

  INTERFACE
    SUBROUTINE ComputeEquationOfStateVariable( D, T, Y, V, dVdT, dVdY )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
      REAL(DP), DIMENSION(:), INTENT(out)           :: V
      REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dVdT
      REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dVdY
    END SUBROUTINE ComputeEquationOfStateVariable
  END INTERFACE

  PROCEDURE (ComputeEquationOfStateVariable), POINTER, PUBLIC :: &
    ComputeSpecificInternalEnergy    => NULL(), &
    ComputeElectronChemicalPotential => NULL(), &
    ComputeProtonChemicalPotential   => NULL(), &
    ComputeNeutronChemicalPotential  => NULL()

  INTERFACE
    SUBROUTINE AuxiliaryEosSubroutine( iB, iE )
      INTEGER, DIMENSION(3), INTENT(in) :: iB, iE
    END SUBROUTINE AuxiliaryEosSubroutine
  END INTERFACE

  PROCEDURE (AuxiliaryEosSubroutine), POINTER, PUBLIC :: &
    ComputeAuxiliary_Fluid => NULL()

  INTERFACE
    PURE FUNCTION AuxiliaryEosFunction_A( PF )
      USE KindModule, ONLY: DP
      USE FluidFieldsModule, ONLY: nPF, nAF
      REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
      REAL(DP), DIMENSION(nAF)             :: AuxiliaryEosFunction_A
    END FUNCTION AuxiliaryEosFunction_A
  END INTERFACE

  PROCEDURE (AuxiliaryEosFunction_A), POINTER, PUBLIC :: &
    Auxiliary_Fluid => NULL()

  INTERFACE
    PURE REAL(DP) FUNCTION AuxiliaryEosFunction_B( PF, AF )
      USE KindModule, ONLY: DP
      USE FluidFieldsModule, ONLY: nPF, nAF
      REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
      REAL(DP), DIMENSION(nAF), INTENT(in) :: AF
    END FUNCTION AuxiliaryEosFunction_B
  END INTERFACE

  PROCEDURE (AuxiliaryEosFunction_B), POINTER, PUBLIC :: &
    Pressure_Primitive => NULL()
  PROCEDURE (AuxiliaryEosFunction_B), POINTER, PUBLIC :: &
    InternalEnergy_Auxiliary => NULL()

  PUBLIC :: InitializeEquationOfState
  PUBLIC :: FinalizeEquationOfState

CONTAINS


  SUBROUTINE InitializeEquationOfState &
               ( EquationOfState_Option, EquationOfStateTableName_Option, &
                 Gamma_IDEAL_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfState_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Gamma_IDEAL_Option

    IF( PRESENT( EquationOfState_Option ) )THEN
      EquationOfState = EquationOfState_Option
    END IF

    IF( PRESENT( EquationOfStateTableName_Option ) )THEN
      EquationOfStateTableName = TRIM( EquationOfStateTableName_Option )
    END IF

    IF( PRESENT( Gamma_IDEAL_Option ) )THEN
      Gamma_IDEAL = Gamma_IDEAL_Option
    END IF

    WRITE(*,*)
    WRITE(*,'(A5,A19,A)') &
      '', 'Equation Of State: ', TRIM( EquationOfState )
    WRITE(*,'(A5,A19)') &
      '', '------------------ '

    SELECT CASE ( EquationOfState )
      CASE( 'IDEAL' )

        ApplyEquationOfState &
          => ApplyEquationOfState_IDEAL
        ComputeThermodynamicStates_Primitive &
          => ComputeThermodynamicStates_Primitive_IDEAL
        ComputeThermodynamicStates_Auxiliary &
          => ComputeThermodynamicStates_Auxiliary_IDEAL
        ComputeChemicalPotentials &
          => ComputeChemicalPotentials_IDEAL
        ComputeAuxiliary_Fluid &
          => ComputeAuxiliary_Fluid_IDEAL
        Auxiliary_Fluid &
          => Auxiliary_Fluid_IDEAL
        Pressure_Primitive &
          => Pressure_Primitive_IDEAL
        InternalEnergy_Auxiliary &
          => InternalEnergy_Auxiliary_IDEAL

      CASE( 'TABLE' )

        WRITE(*,*)
        WRITE(*,'(A7,A12,A)') &
          '', 'Table Name: ', TRIM( EquationOfStateTableName )

        CALL InitializeHDF( )

        CALL ReadEquationOfStateTableHDF &
               ( EOS, TRIM( EquationOfStateTableName ) )

        CALL FinalizeHDF( )

        ApplyEquationOfState &
          => ApplyEquationOfState_TABLE
        ComputeThermodynamicStates_Primitive &
          => ComputeThermodynamicStates_Primitive_TABLE
        ComputeThermodynamicStates_Auxiliary &
          => ComputeThermodynamicStates_Auxiliary_TABLE
        ComputeSpecificInternalEnergy &
          => ComputeSpecificInternalEnergy_TABLE
        ComputeChemicalPotentials &
          => ComputeChemicalPotentials_TABLE
        ComputeElectronChemicalPotential &
          => ComputeElectronChemicalPotential_TABLE
        ComputeProtonChemicalPotential &
          => ComputeProtonChemicalPotential_TABLE
        ComputeNeutronChemicalPotential &
          => ComputeNeutronChemicalPotential_TABLE
        ComputeAuxiliary_Fluid &
          => ComputeAuxiliary_Fluid_TABLE
        Auxiliary_Fluid &
          => Auxiliary_Fluid_TABLE
        Pressure_Primitive &
          => Pressure_Primitive_TABLE
        InternalEnergy_Auxiliary &
          => InternalEnergy_Auxiliary_TABLE

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A27,A5)') &
          ' ', 'Invalid Equation of State: ', EquationOfState
        STOP

    END SELECT

  END SUBROUTINE InitializeEquationOfState


  SUBROUTINE FinalizeEquationOfState

    NULLIFY( ApplyEquationOfState )
    NULLIFY( ComputeThermodynamicStates_Primitive )
    NULLIFY( ComputeThermodynamicStates_Auxiliary )
    NULLIFY( ComputeChemicalPotentials )
    NULLIFY( ComputeAuxiliary_Fluid )
    NULLIFY( Auxiliary_Fluid )
    NULLIFY( Pressure_Primitive )
    NULLIFY( InternalEnergy_Auxiliary )

  END SUBROUTINE FinalizeEquationOfState


  ! -- Ideal Equation of State --


  SUBROUTINE ApplyEquationOfState_IDEAL

    WRITE(*,'(A4,A)') '', 'ApplyEquationOfState_IDEAL'

  END SUBROUTINE ApplyEquationOfState_IDEAL


  SUBROUTINE ComputeThermodynamicStates_Primitive_IDEAL( D, T, Y, Ev, Em, Ne )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Ev, Em, Ne

    WRITE(*,'(A4,A)') '', 'ComputeThermodynamicStates_Primitive_IDEAL'

  END SUBROUTINE ComputeThermodynamicStates_Primitive_IDEAL


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_IDEAL( D, Ev, Ne, T, Em, Y )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: T, Em, Y

    WRITE(*,'(A4,A)') '', 'ComputeThermodynamicStates_Auxiliary_IDEAL'

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_IDEAL


  SUBROUTINE ComputeChemicalPotentials_IDEAL( D, T, Y, Me, Mp, Mn )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Me, Mp, Mn

    WRITE(*,'(A4,A)') '', 'ComputeChemicalPotentials_IDEAL'

  END SUBROUTINE ComputeChemicalPotentials_IDEAL


  SUBROUTINE ComputeElectronChemicalPotential_IDEAL( D, T, Y, Me )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Me

    WRITE(*,'(A4,A)') '', 'ComputeElectronChemicalPotential_IDEAL'

  END SUBROUTINE ComputeElectronChemicalPotential_IDEAL


  SUBROUTINE ComputeProtonChemicalPotential_IDEAL( D, T, Y, Mp )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Mp

    WRITE(*,'(A4,A)') '', 'ComputeProtonChemicalPotential_IDEAL'

  END SUBROUTINE ComputeProtonChemicalPotential_IDEAL


  SUBROUTINE ComputeNeutronChemicalPotential_IDEAL( D, T, Y, Mn )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Mn

    WRITE(*,'(A4,A)') '', 'ComputeNeutronChemicalPotential_IDEAL'

  END SUBROUTINE ComputeNeutronChemicalPotential_IDEAL


  SUBROUTINE ComputeAuxiliary_Fluid_IDEAL( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iX1, iX2, iX3

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          uAF(:,iX1,iX2,iX3,iAF_P) &
            = ( Gamma_IDEAL - 1.0_DP ) * uPF(:,iX1,iX2,iX3,iPF_E)

          uAF(:,iX1,iX2,iX3,iAF_Gm) &
            = Gamma_IDEAL

          uAF(:,iX1,iX2,iX3,iAF_Cs) &
            = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - 1.0_DP ) &
                      * uPF(:,iX1,iX2,iX3,iAF_E) &
                          / uPF(:,iX1,iX2,iX3,iPF_D) )

        END DO
      END DO
    END DO

  END SUBROUTINE ComputeAuxiliary_Fluid_IDEAL


  PURE FUNCTION Auxiliary_Fluid_IDEAL( PF )

    REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
    REAL(DP), DIMENSION(nAF)             :: Auxiliary_Fluid_IDEAL

    Auxiliary_Fluid_IDEAL(iAF_P)  &
      = ( Gamma_IDEAL - 1.0_DP ) * PF(iPF_E)
    Auxiliary_Fluid_IDEAL(iAF_Gm) &
      = Gamma_IDEAL
    Auxiliary_Fluid_IDEAL(iAF_E)  &
      = PF(iPF_E) / PF(iPF_D)
    Auxiliary_Fluid_IDEAL(iAF_Cs) &
      = SQRT( Gamma_IDEAL * ( Gamma_IDEAL - 1.0_DP ) &
                * PF(iPF_E) / PF(iPF_D) )

    RETURN
  END FUNCTION Auxiliary_Fluid_IDEAL


  PURE REAL(DP) FUNCTION Pressure_Primitive_IDEAL( PF, AF )

    REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
    REAL(DP), DIMENSION(nAF), INTENT(in) :: AF

    Pressure_Primitive_IDEAL &
      = ( Gamma_IDEAL - 1.0_DP ) * PF(iPF_E)

    RETURN
  END FUNCTION Pressure_Primitive_IDEAL


  PURE REAL(DP) FUNCTION InternalEnergy_Auxiliary_IDEAL( PF, AF )

    REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
    REAL(DP), DIMENSION(nAF), INTENT(in) :: AF

    InternalEnergy_Auxiliary_IDEAL &
      = AF(iAF_P) / ( Gamma_IDEAL - 1.0_DP )

    RETURN
  END FUNCTION InternalEnergy_Auxiliary_IDEAL


  ! -- Tabulated Equation of State (through weaklib) --


  SUBROUTINE ApplyEquationOfState_TABLE

    REAL(DP), DIMENSION(:), ALLOCATABLE :: D_1D, T_1D, Y_1D, TMP

    ASSOCIATE &
      ( iD_T  => EOS % TS % Indices % iRho,                       &
        iT_T  => EOS % TS % Indices % iT,                         &
        iY_T  => EOS % TS % Indices % iYe,                        &
        iP_T  => EOS % DV % Indices % iPressure,                  &
        iS_T  => EOS % DV % Indices % iEntropyPerBaryon,          &
        iE_T  => EOS % DV % Indices % iInternalEnergyDensity,     &
        iMe_T => EOS % DV % Indices % iElectronChemicalPotential, &
        iMp_T => EOS % DV % Indices % iProtonChemicalPotential,   &
        iMn_T => EOS % DV % Indices % iNeutronChemicalPotential )

    ASSOCIATE &
      ( D_T  => EOS % TS % States(iD_T) % Values,     &
        T_T  => EOS % TS % States(iT_T) % Values,     &
        Y_T  => EOS % TS % States(iY_T) % Values,     &
        Log  => EOS % TS % LogInterp,                 &
        P_T  => EOS % DV % Variables( iP_T) % Values, &
        S_T  => EOS % DV % Variables( iS_T) % Values, &
        E_T  => EOS % DV % Variables( iE_T) % Values, &
        Me_T => EOS % DV % Variables(iMe_T) % Values, &
        Mp_T => EOS % DV % Variables(iMp_T) % Values, &
        Mn_T => EOS % DV % Variables(iMn_T) % Values, &
        OS   => EOS % DV % Offsets )

    ALLOCATE( D_1D(nDOFX*PRODUCT( nX )) )
    ALLOCATE( T_1D(nDOFX*PRODUCT( nX )) )
    ALLOCATE( Y_1D(nDOFX*PRODUCT( nX )) )
    ALLOCATE( TMP (nDOFX*PRODUCT( nX )) )

    CALL MapTo1D( uPF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iPF_D ), D_1D )
    CALL MapTo1D( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_T ), T_1D )
    CALL MapTo1D( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Ye), Y_1D )

    ! --- Convert to weaklib units ---

    D_1D = D_1D / ( Gram / Centimeter**3 )
    T_1D = T_1D / ( Kelvin )
    Y_1D = Y_1D

    ! --- Interpolate Pressure ----------------------------------------

    CALL LogInterpolateSingleVariable &
           ( D_1D, T_1D, Y_1D, D_T, T_T, Y_T, Log, OS(iP_T), P_T, TMP )

    CALL MapFrom1D( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_P ), &
                    TMP * Dyne / Centimeter**2 )

    ! --- Interpolate Entropy Per Baryon ------------------------------

    CALL LogInterpolateSingleVariable &
           ( D_1D, T_1D, Y_1D, D_T, T_T, Y_T, Log, OS(iS_T), S_T, TMP )

    CALL MapFrom1D( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_S ), &
                    TMP * BoltzmannConstant )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL LogInterpolateSingleVariable &
           ( D_1D, T_1D, Y_1D, D_T, T_T, Y_T, Log, OS(iE_T), E_T, TMP )

    CALL MapFrom1D( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_E ), &
                    TMP * Erg / Gram )

    ! --- Interpolate Electron Chemical Potential ---------------------

    CALL LogInterpolateSingleVariable &
           ( D_1D, T_1D, Y_1D, D_T, T_T, Y_T, Log, OS(iMe_T), Me_T, TMP )

    CALL MapFrom1D( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Me), &
                    TMP * MeV )

    ! --- Interpolate Proton Chemical Potential -----------------------

    CALL LogInterpolateSingleVariable &
           ( D_1D, T_1D, Y_1D, D_T, T_T, Y_T, Log, OS(iMp_T), Mp_T, TMP )

    CALL MapFrom1D( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Mp), &
                    TMP * MeV )

    ! --- Interpolate Neutron Chemical Potential ----------------------

    CALL LogInterpolateSingleVariable &
           ( D_1D, T_1D, Y_1D, D_T, T_T, Y_T, Log, OS(iMn_T), Mn_T, TMP )

    CALL MapFrom1D( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Mn), &
                    TMP * MeV )

    DEALLOCATE( D_1D, T_1D, Y_1D, TMP )

    END ASSOCIATE ! D_T, etc.

    END ASSOCIATE ! iD_T, etc.

  END SUBROUTINE ApplyEquationOfState_TABLE


  SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE( D, T, Y, Ev, Em, Ne )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Ev, Em, Ne

    REAL(DP), DIMENSION(:), ALLOCATABLE :: TMP

    ALLOCATE( TMP(SIZE( D )) )

    ASSOCIATE &
      ( iD_T => EOS % TS % Indices % iRho, &
        iT_T => EOS % TS % Indices % iT,   &
        iY_T => EOS % TS % Indices % iYe,  &
        iE_T => EOS % DV % Indices % iInternalEnergyDensity )

    ASSOCIATE &
      ( D_T => EOS % TS % States(iD_T) % Values,    &
        T_T => EOS % TS % States(iT_T) % Values,    &
        Y_T => EOS % TS % States(iY_T) % Values,    &
        Log => EOS % TS % LogInterp,                &
        E_T => EOS % DV % Variables(iE_T) % Values, &
        OS  => EOS % DV % Offsets )

    ! --- Interpolate Specific Internal Energy ------------------------

    CALL LogInterpolateSingleVariable &
           ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, D_T, T_T, Y_T, &
             Log, OS(iE_T), E_T, TMP )

    Em(:) = TMP(:) * Erg / Gram      ! --- Internal Energy per Mass
    Ev(:) = Em(:) * D(:)             ! --- Internal Energy per Volume
    Ne(:) = Y(:) * D(:) / BaryonMass ! --- Electrons per Volume

    END ASSOCIATE ! D_T, etc.

    END ASSOCIATE ! iD_T, etc.

    DEALLOCATE( TMP )

  END SUBROUTINE ComputeThermodynamicStates_Primitive_TABLE


  SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE( D, Ev, Ne, T, Em, Y )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, Ev, Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: T, Em, Y

    INTEGER                :: iS
    REAL(DP), DIMENSION(1) :: TMP

    ASSOCIATE &
      ( iD_T => EOS % TS % Indices % iRho, &
        iT_T => EOS % TS % Indices % iT,   &
        iY_T => EOS % TS % Indices % iYe,  &
        iE_T => EOS % DV % Indices % iInternalEnergyDensity )

    ASSOCIATE &
      ( D_T => EOS % TS % States(iD_T) % Values,    &
        T_T => EOS % TS % States(iT_T) % Values,    &
        Y_T => EOS % TS % States(iY_T) % Values,    &
        Log => EOS % TS % LogInterp,                &
        E_T => EOS % DV % Variables(iE_T) % Values, &
        OS  => EOS % DV % Offsets )

    Em(:) = Ev(:) / D(:)              ! --- Internal Energy per Mass
    Y(:)  = Ne(:) / D(:) * BaryonMass ! --- Electron Fraction

    DO iS = 1, SIZE( D )

      CALL ComputeTempFromIntEnergy &
             ( D(iS) / ( Gram / Centimeter**3 ), Em(iS) / ( Erg / Gram ), &
               Y(iS), D_T, T_T, Y_T, Log, E_T, OS(iE_T), TMP )

      T(iS) = TMP(1) * Kelvin

    END DO

    END ASSOCIATE ! D_T, etc.

    END ASSOCIATE ! iD_T, etc.

  END SUBROUTINE ComputeThermodynamicStates_Auxiliary_TABLE


  SUBROUTINE ComputeSpecificInternalEnergy_TABLE &
               ( D, T, Y, E, dEdT_Option, dEdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: E
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dEdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dEdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: TMP
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dTMP

    ComputeDerivatives = .FALSE.
    IF( ALL( [ PRESENT( dEdT_Option ), PRESENT( dEdY_Option ) ] ) ) &
      ComputeDerivatives = .TRUE.

    ALLOCATE( TMP(SIZE( D )), dTMP(SIZE( D ), 3) )

    ASSOCIATE &
      ( iD_T => EOS % TS % Indices % iRho, &
        iT_T => EOS % TS % Indices % iT,   &
        iY_T => EOS % TS % Indices % iYe,  &
        iE_T => EOS % DV % Indices % iInternalEnergyDensity )

    ASSOCIATE &
      ( D_T => EOS % TS % States(iD_T) % Values,    &
        T_T => EOS % TS % States(iT_T) % Values,    &
        Y_T => EOS % TS % States(iY_T) % Values,    &
        Log => EOS % TS % LogInterp,                &
        E_T => EOS % DV % Variables(iE_T) % Values, &
        OS  => EOS % DV % Offsets )

    ! --- Interpolate Specific Internal Energy ------------------------

    IF( ComputeDerivatives )THEN

      CALL LogInterpolateDifferentiateSingleVariable &
             ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               D_T, T_T, Y_T, Log, OS(iE_T), E_T, TMP, dTMP )

      E(:) = TMP(:) * Erg / Gram

      dEdT_Option(:) &
        = dTMP(:,2) * Erg / Gram / Kelvin

      dEdY_Option(:) &
        = dTMP(:,3) * Erg / Gram

    ELSE

      CALL LogInterpolateSingleVariable &
             ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               D_T, T_T, Y_T, Log, OS(iE_T), E_T, TMP )

      E(:) = TMP(:) * Erg / Gram

    END IF

    END ASSOCIATE ! D_T, etc.

    END ASSOCIATE ! iD_T, etc.

    DEALLOCATE( TMP, dTMP )

  END SUBROUTINE ComputeSpecificInternalEnergy_TABLE


  SUBROUTINE ComputeChemicalPotentials_TABLE( D, T, Y, Me, Mp, Mn )

    REAL(DP), DIMENSION(:), INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Me, Mp, Mn

    REAL(DP), DIMENSION(:), ALLOCATABLE :: TMP

    ALLOCATE( TMP(SIZE( D )) )

    ASSOCIATE &
      ( iD_T  => EOS % TS % Indices % iRho,                       &
        iT_T  => EOS % TS % Indices % iT,                         &
        iY_T  => EOS % TS % Indices % iYe,                        &
        iMe_T => EOS % DV % Indices % iElectronChemicalPotential, &
        iMp_T => EOS % DV % Indices % iProtonChemicalPotential,   &
        iMn_T => EOS % DV % Indices % iNeutronChemicalPotential )

    ASSOCIATE &
      ( D_T  => EOS % TS % States(iD_T) % Values,     &
        T_T  => EOS % TS % States(iT_T) % Values,     &
        Y_T  => EOS % TS % States(iY_T) % Values,     &
        Log  => EOS % TS % LogInterp,                 &
        Me_T => EOS % DV % Variables(iMe_T) % Values, &
        Mp_T => EOS % DV % Variables(iMp_T) % Values, &
        Mn_T => EOS % DV % Variables(iMn_T) % Values, &
        OS   => EOS % DV % Offsets )

    ! --- Interpolate Electron Chemical Potential ---------------------

    CALL LogInterpolateSingleVariable &
           ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, D_T, T_T, Y_T, &
             Log, OS(iMe_T), Me_T, TMP )

    Me(:) = TMP(:) * MeV

    ! --- Interpolate Proton Chemical Potential -----------------------

    CALL LogInterpolateSingleVariable &
           ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, D_T, T_T, Y_T, &
             Log, OS(iMp_T), Mp_T, TMP )

    Mp(:) = TMP(:) * MeV

    ! --- Interpolate Neutron Chemical Potential ----------------------

    CALL LogInterpolateSingleVariable &
           ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, D_T, T_T, Y_T, &
             Log, OS(iMn_T), Mn_T, TMP )

    Mn(:) = TMP(:) * MeV

    END ASSOCIATE ! D_T, etc.

    END ASSOCIATE ! iD_T, etc.

    DEALLOCATE( TMP )

  END SUBROUTINE ComputeChemicalPotentials_TABLE


  SUBROUTINE ComputeElectronChemicalPotential_TABLE &
               ( D, T, Y, Me, dMedT_Option, dMedY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: Me
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMedT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMedY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: TMP
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dTMP

    ComputeDerivatives = .FALSE.
    IF( ALL( [ PRESENT( dMedT_Option ), PRESENT( dMedY_Option ) ] ) ) &
      ComputeDerivatives = .TRUE.

    ALLOCATE( TMP(SIZE( D )), dTMP(SIZE( D ), 3) )

    ASSOCIATE &
      ( iD_T  => EOS % TS % Indices % iRho, &
        iT_T  => EOS % TS % Indices % iT,   &
        iY_T  => EOS % TS % Indices % iYe,  &
        iMe_T => EOS % DV % Indices % iElectronChemicalPotential )

    ASSOCIATE &
      ( D_T  => EOS % TS % States(iD_T) % Values,     &
        T_T  => EOS % TS % States(iT_T) % Values,     &
        Y_T  => EOS % TS % States(iY_T) % Values,     &
        Log  => EOS % TS % LogInterp,                 &
        Me_T => EOS % DV % Variables(iMe_T) % Values, &
        OS   => EOS % DV % Offsets )

    ! --- Interpolate Electron Chemical Potential ---------------------

    IF( ComputeDerivatives )THEN

      CALL LogInterpolateDifferentiateSingleVariable &
             ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               D_T, T_T, Y_T, Log, OS(iMe_T), Me_T, TMP, dTMP )

      Me(:) = TMP(:) * MeV

      dMedT_Option(:) &
        = dTMP(:,2) * MeV / Kelvin

      dMedY_Option(:) &
        = dTMP(:,3) * MeV

    ELSE

      CALL LogInterpolateSingleVariable &
             ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               D_T, T_T, Y_T, Log, OS(iMe_T), Me_T, TMP )

      Me(:) = TMP(:) * MeV

    END IF

    END ASSOCIATE ! D_T, etc.

    END ASSOCIATE ! iD_T, etc.

    DEALLOCATE( TMP, dTMP )

  END SUBROUTINE ComputeElectronChemicalPotential_TABLE


  SUBROUTINE ComputeProtonChemicalPotential_TABLE &
               ( D, T, Y, Mp, dMpdT_Option, dMpdY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: Mp
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMpdT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMpdY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: TMP
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dTMP

    ComputeDerivatives = .FALSE.
    IF( ALL( [ PRESENT( dMpdT_Option ), PRESENT( dMpdY_Option ) ] ) ) &
      ComputeDerivatives = .TRUE.

    ALLOCATE( TMP(SIZE( D )), dTMP(SIZE( D ), 3) )

    ASSOCIATE &
      ( iD_T  => EOS % TS % Indices % iRho, &
        iT_T  => EOS % TS % Indices % iT,   &
        iY_T  => EOS % TS % Indices % iYe,  &
        iMp_T => EOS % DV % Indices % iProtonChemicalPotential )

    ASSOCIATE &
      ( D_T  => EOS % TS % States(iD_T) % Values,     &
        T_T  => EOS % TS % States(iT_T) % Values,     &
        Y_T  => EOS % TS % States(iY_T) % Values,     &
        Log  => EOS % TS % LogInterp,                 &
        Mp_T => EOS % DV % Variables(iMp_T) % Values, &
        OS   => EOS % DV % Offsets )

    ! --- Interpolate Proton Chemical Potential -----------------------

    IF( ComputeDerivatives )THEN

      CALL LogInterpolateDifferentiateSingleVariable &
             ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               D_T, T_T, Y_T, Log, OS(iMp_T), Mp_T, TMP, dTMP )

      Mp(:) = TMP(:) * MeV

      dMpdT_Option(:) &
        = dTMP(:,2) * MeV / Kelvin

      dMpdY_Option(:) &
        = dTMP(:,3) * MeV

    ELSE

      CALL LogInterpolateSingleVariable &
             ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               D_T, T_T, Y_T, Log, OS(iMp_T), Mp_T, TMP )

      Mp(:) = TMP(:) * MeV

    END IF

    END ASSOCIATE ! D_T, etc.

    END ASSOCIATE ! iD_T, etc.

    DEALLOCATE( TMP, dTMP )

  END SUBROUTINE ComputeProtonChemicalPotential_TABLE


  SUBROUTINE ComputeNeutronChemicalPotential_TABLE &
               ( D, T, Y, Mn, dMndT_Option, dMndY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: Mn
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMndT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dMndY_Option

    LOGICAL :: ComputeDerivatives
    REAL(DP), DIMENSION(:), ALLOCATABLE   :: TMP
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dTMP

    ComputeDerivatives = .FALSE.
    IF( ALL( [ PRESENT( dMndT_Option ), PRESENT( dMndY_Option ) ] ) ) &
      ComputeDerivatives = .TRUE.

    ALLOCATE( TMP(SIZE( D )), dTMP(SIZE( D ), 3) )

    ASSOCIATE &
      ( iD_T  => EOS % TS % Indices % iRho, &
        iT_T  => EOS % TS % Indices % iT,   &
        iY_T  => EOS % TS % Indices % iYe,  &
        iMn_T => EOS % DV % Indices % iNeutronChemicalPotential )

    ASSOCIATE &
      ( D_T  => EOS % TS % States(iD_T) % Values,     &
        T_T  => EOS % TS % States(iT_T) % Values,     &
        Y_T  => EOS % TS % States(iY_T) % Values,     &
        Log  => EOS % TS % LogInterp,                 &
        Mn_T => EOS % DV % Variables(iMn_T) % Values, &
        OS   => EOS % DV % Offsets )

    ! --- Interpolate Neutron Chemical Potential ----------------------

    IF( ComputeDerivatives )THEN

      CALL LogInterpolateDifferentiateSingleVariable &
             ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               D_T, T_T, Y_T, Log, OS(iMn_T), Mn_T, TMP, dTMP )

      Mn(:) = TMP(:) * MeV

      dMndT_Option(:) &
        = dTMP(:,2) * MeV / Kelvin

      dMndY_Option(:) &
        = dTMP(:,3) * MeV

    ELSE

      CALL LogInterpolateSingleVariable &
             ( D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               D_T, T_T, Y_T, Log, OS(iMn_T), Mn_T, TMP )

      Mn(:) = TMP(:) * MeV

    END IF

    END ASSOCIATE ! D_T, etc.

    END ASSOCIATE ! iD_T, etc.

    DEALLOCATE( TMP, dTMP )

  END SUBROUTINE ComputeNeutronChemicalPotential_TABLE


  SUBROUTINE ComputeAuxiliary_Fluid_TABLE( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

  END SUBROUTINE ComputeAuxiliary_Fluid_TABLE


  PURE FUNCTION Auxiliary_Fluid_TABLE( PF )

    REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
    REAL(DP), DIMENSION(nAF)             :: Auxiliary_Fluid_TABLE

    Auxiliary_Fluid_TABLE(1:nAF) &
      = 0.0_DP

    RETURN
  END FUNCTION Auxiliary_Fluid_TABLE


  PURE REAL(DP) FUNCTION Pressure_Primitive_TABLE( PF, AF )

    REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
    REAL(DP), DIMENSION(nAF), INTENT(in) :: AF

    Pressure_Primitive_TABLE &
      = 0.0_DP

    RETURN
  END FUNCTION Pressure_Primitive_TABLE


  PURE REAL(DP) FUNCTION InternalEnergy_Auxiliary_TABLE( PF, AF )

    REAL(DP), DIMENSION(nPF), INTENT(in) :: PF
    REAL(DP), DIMENSION(nAF), INTENT(in) :: AF

    InternalEnergy_Auxiliary_TABLE &
      = 0.0_DP

    RETURN
  END FUNCTION InternalEnergy_Auxiliary_TABLE


END MODULE EquationOfStateModule
