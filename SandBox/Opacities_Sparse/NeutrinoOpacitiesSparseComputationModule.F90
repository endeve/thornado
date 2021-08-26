#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_OPACITY
#endif

MODULE NeutrinoOpacitiesSparseComputationModule

  USE TasmanianSG
  USE KindModule, ONLY: &
    DP, Zero, One
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    Centimeter, &
    Gram, &
    Kelvin, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFE, nDOFX
  USE DeviceModule, ONLY: &
    QueryOnGpu
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &
    NodeNumbersX, &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE, &
    ComputeSpecificInternalEnergy_TABLE
  USE OpacityModule_TABLE, ONLY: &
#ifdef MICROPHYSICS_WEAKLIB
    OS_EmAb, OS_Iso, OS_NES, OS_Pair, &
    EmAb_T, Iso_T, NES_T, Pair_T, &
    NES_AT, Pair_AT, &
#endif
    LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar
  USE NeutrinoOpacitiesModule, ONLY: &
    f_EQ, opEC, opES, opIS, opPP

#ifdef MICROPHYSICS_WEAKLIB

  ! --- weaklib modules ---

  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateSingleVariable_1D3D_Custom, &
    LogInterpolateSingleVariable_1D3D_Custom_Point, &
    LogInterpolateSingleVariable_2D2D_Custom, &
    LogInterpolateSingleVariable_2D2D_Custom_Point, &
    LogInterpolateSingleVariable_2D2D_Custom_Aligned, &
    LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point, &
    LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Point, &
    LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P

  ! ----------------------------------------------

#endif

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeSparseGridOpacities
  PUBLIC :: FinalizeSparseGridOpacities
  PUBLIC :: ComputeEquilibriumDistributions_Points_SG
  PUBLIC :: ComputeEquilibriumDistributionAndDerivatives_Points_SG
  PUBLIC :: ComputeNeutrinoOpacities_EC_Points_SG
  PUBLIC :: ComputeNeutrinoOpacities_ES_Points_SG
  PUBLIC :: ComputeNeutrinoOpacities_NES_Point_SG
  PUBLIC :: ComputeNeutrinoOpacities_NES_Points_SG
  PUBLIC :: ComputeNeutrinoOpacitiesRates_NES_Points
  PUBLIC :: ComputeNeutrinoOpacities_Pair_Point_SG
  PUBLIC :: ComputeNeutrinoOpacities_Pair_Points_SG
  PUBLIC :: ComputeNeutrinoOpacitiesRates_Pair_Points
  PUBLIC :: FermiDirac
  PUBLIC :: dFermiDiracdT
  PUBLIC :: dFermiDiracdY

  REAL(DP), PARAMETER :: Log1d100 = LOG( 1.0d100 )
  REAL(DP), PARAMETER :: UnitD    = Gram / Centimeter**3
  REAL(DP), PARAMETER :: UnitT    = Kelvin
  REAL(DP), PARAMETER :: UnitY    = One
  REAL(DP), PARAMETER :: UnitE    = MeV
  REAL(DP), PARAMETER :: UnitEta  = One
  REAL(DP), PARAMETER :: UnitEC   = One / Centimeter
  REAL(DP), PARAMETER :: UnitES   = One / Centimeter
  REAL(DP), PARAMETER :: UnitNES  = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: UnitPair = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: cv       = 0.96d+00 ! weak interaction constant
  REAL(DP), PARAMETER :: ca       = 0.50d+00 ! weak interaction constant

  REAL(DP), PARAMETER :: C1(iNuE:iNuE_Bar) = [ ( cv + ca )**2, ( cv - ca )**2 ]
  REAL(DP), PARAMETER :: C2(iNuE:iNuE_Bar) = [ ( cv - ca )**2, ( cv + ca )**2 ]

  CHARACTER(256) :: &
    OpacityGridName_EmAb, &
    OpacityGridName_Iso,  &
    OpacityGridName_NES,  &
    OpacityGridName_Pair

  type(TasmanianSparseGrid), PUBLIC :: gridEmAb, gridIso, gridNES, gridPair

  INTERFACE FermiDirac
    MODULE PROCEDURE FermiDirac_Scalar
    MODULE PROCEDURE FermiDirac_Vector
  END INTERFACE

  INTERFACE dFermiDiracdT
    MODULE PROCEDURE dFermiDiracdT_Scalar
    MODULE PROCEDURE dFermiDiracdT_Vector
  END INTERFACE

  INTERFACE dFermiDiracdY
    MODULE PROCEDURE dFermiDiracdY_Scalar
    MODULE PROCEDURE dFermiDiracdY_Vector
  END INTERFACE

CONTAINS

    SUBROUTINE InitializeSparseGridOpacities &
      ( OpacityGridName_EmAb_Option, OpacityGridName_Iso_Option, &
        OpacityGridName_NES_Option, OpacityGridName_Pair_Option )

      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityGridName_EmAb_Option
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityGridName_Iso_Option
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityGridName_NES_Option
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityGridName_Pair_Option

      LOGICAL :: Include_EmAb
      LOGICAL :: Include_Iso
      LOGICAL :: Include_NES
      LOGICAL :: Include_Pair
      LOGICAL :: SGReadError

      CALL tsgAllocateGrid(gridEmAb)
      CALL tsgAllocateGrid(gridIso)
      CALL tsgAllocateGrid(gridNES)
      CALL tsgAllocateGrid(gridPair)

      IF( PRESENT( OpacityGridName_EmAb_Option ) &
          .AND. ( LEN( OpacityGridName_EmAb_Option ) > 1 ) )THEN
        OpacityGridName_EmAb = TRIM( OpacityGridName_EmAb_Option )
        Include_EmAb = .TRUE.
        SGReadError = tsgRead(gridEmAb, OpacityGridName_EmAb)
      ELSE
        OpacityGridName_EmAb = ''
        Include_EmAb = .FALSE.
      END IF

      IF( PRESENT( OpacityGridName_Iso_Option ) &
          .AND. ( LEN( OpacityGridName_Iso_Option ) > 1 ) )THEN
        OpacityGridName_Iso = TRIM( OpacityGridName_Iso_Option )
        Include_Iso = .TRUE.
        SGReadError = tsgRead(gridIso, OpacityGridName_Iso)
      ELSE
        OpacityGridName_Iso = ''
        Include_Iso = .FALSE.
      END IF

      IF( PRESENT( OpacityGridName_NES_Option ) &
          .AND. ( LEN( OpacityGridName_NES_Option ) > 1 ) )THEN
        OpacityGridName_NES = TRIM( OpacityGridName_NES_Option )
        Include_NES = .TRUE.
        SGReadError = tsgRead(gridNES, OpacityGridName_NES)
      ELSE
        OpacityGridName_NES = ''
        Include_NES = .FALSE.
      END IF

      IF( PRESENT( OpacityGridName_Pair_Option ) &
          .AND. ( LEN( OpacityGridName_Pair_Option ) > 1 ) )THEN
        OpacityGridName_Pair = TRIM( OpacityGridName_Pair_Option )
        Include_Pair = .TRUE.
        SGReadError = tsgRead(gridPair, OpacityGridName_Pair)
      ELSE
        OpacityGridName_Pair = ''
        Include_Pair = .FALSE.
      END IF
  END SUBROUTINE InitializeSparseGridOpacities

  SUBROUTINE FinalizeSparseGridOpacities

    CALL tsgDeallocateGrid(gridEmAb)
    CALL tsgDeallocateGrid(gridIso)
    CALL tsgDeallocateGrid(gridNES)
    CALL tsgDeallocateGrid(gridPair)

  END SUBROUTINE FinalizeSparseGridOpacities


  SUBROUTINE ComputeEquilibriumDistributions_Points_SG &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, f0_1, f0_2, iS_1, iS_2 )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    REAL(DP), INTENT(out) :: f0_1(:,:), f0_2(:,:)
    INTEGER,  INTENT(in)  :: iS_1, iS_2

    INTEGER  :: iX, iE, iSpecies
    REAL(DP) :: Me(iX_B:iX_E), Mp(iX_B:iX_E), Mn(iX_B:iX_E)
    REAL(DP) :: Mnu, kT, FD_Exp
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( f0_1, f0_2 )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeEquilibriumDistributions_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_Points]   Y missing'
      IF ( .not. QueryOnGPU( f0_1 ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_Points]   f0_1 missing'
      IF ( .not. QueryOnGPU( f0_2 ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_Points]   f0_2 missing'
    END IF
#endif

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E

        kT = BoltzmannConstant * T(iX)

        Mnu = ( Me(iX) + Mp(iX) ) - Mn(iX)

        IF ( iS_1 == iNuE ) THEN
          f0_1(iE,iX) = FermiDirac( E(iE), +Mnu, kT )
        ELSE IF ( iS_1 == iNuE_Bar ) THEN
          f0_1(iE,iX) = FermiDirac( E(iE), -Mnu, kT )
        ELSE
          f0_1(iE,iX) = FermiDirac( E(iE), Zero, kT )
        END IF

        IF ( iS_2 == iNuE ) THEN
          f0_2(iE,iX) = FermiDirac( E(iE), +Mnu, kT )
        ELSE IF ( iS_2 == iNuE_Bar ) THEN
          f0_2(iE,iX) = FermiDirac( E(iE), -Mnu, kT )
        ELSE
          f0_2(iE,iX) = FermiDirac( E(iE), Zero, kT )
        END IF

      END DO
    END DO

  END SUBROUTINE ComputeEquilibriumDistributions_Points_SG



  SUBROUTINE ComputeEquilibriumDistributionAndDerivatives_Points_SG &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, f0, df0dY, df0dU, iSpecies )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    REAL(DP), INTENT(out) :: f0   (:,:)
    REAL(DP), INTENT(out) :: df0dY(:,:)
    REAL(DP), INTENT(out) :: df0dU(:,:)
    INTEGER,  INTENT(in)  :: iSpecies

    REAL(DP), DIMENSION(iX_B:iX_E) :: Me,  dMedT , dMedY
    REAL(DP), DIMENSION(iX_B:iX_E) :: Mp,  dMpdT , dMpdY
    REAL(DP), DIMENSION(iX_B:iX_E) :: Mn,  dMndT , dMndY
    REAL(DP), DIMENSION(iX_B:iX_E) :: Mnu, dMnudT, dMnudY
    REAL(DP), DIMENSION(iX_B:iX_E) :: U,   dUdT,   dUdY, dUdD

    REAL(DP) :: kT, df0dT_Y, df0dY_T
    INTEGER  :: iX, iE
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( f0, df0dY, df0dU )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points_SG] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points_SG]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points_SG]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points_SG]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points_SG]   Y missing'
      IF ( .not. QueryOnGPU( f0 ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points_SG]   f0 missing'
      IF ( .not. QueryOnGPU( df0dY ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points_SG]   df0dY missing'
      IF ( .not. QueryOnGPU( df0dU ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points_SG]   df0dU missing'
    END IF
#endif

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me, dUdD, dMedT, dMedY )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp, dUdD, dMpdT, dMpdY )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn, dUdD, dMndT, dMndY )

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( D, T, Y, U,  dUdD, dUdT,  dUdY  )


    IF ( iSpecies == iNuE ) THEN

      DO iX = iX_B, iX_E
        Mnu   (iX) = ( Me   (iX) + Mp   (iX) ) - Mn   (iX)
        dMnudT(iX) = ( dMedT(iX) + dMpdT(iX) ) - dMndT(iX)
        dMnudY(iX) = ( dMedY(iX) + dMpdY(iX) ) - dMndY(iX)
      END DO

    ELSE IF ( iSpecies == iNuE_Bar ) THEN

      DO iX = iX_B, iX_E
        Mnu   (iX) = Mn   (iX) - ( Me   (iX) + Mp   (iX) )
        dMnudT(iX) = dMndT(iX) - ( dMedT(iX) + dMpdT(iX) )
        dMnudY(iX) = dMndY(iX) - ( dMedY(iX) + dMpdY(iX) )
      END DO

    ELSE

      DO iX = iX_B, iX_E
        Mnu   (iX) = Zero
        dMnudT(iX) = Zero
        dMnudY(iX) = Zero
      END DO

    END IF

    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E

        kT = BoltzmannConstant * T(iX)

        f0(iE,iX) = FermiDirac   ( E(iE), Mnu(iX), kT )
        df0dT_Y   = dFermiDiracdT( E(iE), Mnu(iX), kT, dMnudT(iX), T(iX) ) ! Constant T
        df0dY_T   = dFermiDiracdY( E(iE), Mnu(iX), kT, dMnudY(iX), T(iX) ) ! Constant Y

        df0dU(iE,iX) = df0dT_Y / dUdT(iX)
        df0dY(iE,iX) = df0dY_T - df0dU(iE,iX) * dUdY(iX)

      END DO
    END DO

  END SUBROUTINE ComputeEquilibriumDistributionAndDerivatives_Points_SG

  SUBROUTINE ComputeNeutrinoOpacities_EC_Points_SG &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iSpecies, opEC_Points )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---
    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iSpecies
    REAL(DP), INTENT(out) :: opEC_Points(:,:)

    INTEGER  :: iX, iE, nE, nX
    INTEGER  :: k
    REAL(DP) :: LogE_P(iE_B:iE_E), LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    REAL(DP), allocatable :: Points(:,:), Values(:,:)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opEC_Points )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points_SG] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points_SG]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points_SG]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points_SG]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points_SG]   Y missing'
      IF ( .not. QueryOnGPU( opEC_Points ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points_SG]   opEC_Points missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

    nE = iE_E - iE_B + 1
    nX = iX_E - iX_B + 1

    DO iX = iX_B, iX_E
      LogD_P(iX) = LOG10( D(iX) / UnitD )
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      Y_P(iX) = Y(iX) / UnitY
    END DO

    DO iE = iE_B, iE_E
      LogE_P(iE) = LOG10( E(iE) / UnitE )
    END DO

    allocate(Points(4, nX * nE))
    allocate(Values(2, nX * nE))

    k = 0
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E
        k = k + 1
        Points(1,k) = Y_P(iX)
        Points(2,k) = LogT_P(iX)
        Points(3,k) = LogD_P(iX)
        Points(4,k) = LogE_P(iE)
      END DO
    END DO


    CALL tsgEvaluateBatch(gridEmAb, Points,  nX * nE, Values)
    ! CALL tsgEvaluateBatch(gridEmAb, Points,  2, Values)

    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "Points(:,1) = ", Points(:,1)
    ! WRITE(*,*) "Points(:,2) = ", Points(:,2)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    !
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "Values(:,1) = ", Values(:,1)
    ! WRITE(*,*) "Values(:,2) = ", Values(:,2)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"



    k = 0
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E
        k = k + 1
        ! opEC_Points(iE,iX) = Values(iSpecies,k)
        ! opEC_Points(iE,iX) = 10.0d0**(Values(iSpecies,k)) - OS_EmAb(iSpecies)
        ! opEC_Points(iE,iX) = 10.0d0**(Values(iSpecies,k) - 20.0d0) - OS_EmAb(iSpecies)
        opEC_Points(iE,iX) = 10.0d0**(Values(iSpecies,k) - 20.0d0) - 0.0d0
      END DO
    END DO
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "Values(iSpecies,1) = ", Values(iSpecies,1)
    ! WRITE(*,*) "Values(iSpecies,2) = ", Values(iSpecies,2)
    ! WRITE(*,*) "OS_EmAb(iSpecies) = ", OS_EmAb(iSpecies)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "opEC_Points(1,1) = ", opEC_Points(1,1)
    ! WRITE(*,*) "opEC_Points(2,1) = ", opEC_Points(2,1)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"

    ! CALL LogInterpolateSingleVariable_1D3D_Custom &
    !        ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
    !          OS_EmAb(iSpecies), EmAb_T(:,:,:,:,iSpecies), opEC_Points, &
    !          GPU_Option = do_gpu )
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "opEC_Points(1,1) = ", opEC_Points(1,1)
    ! WRITE(*,*) "opEC_Points(2,1) = ", opEC_Points(2,1)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"

    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E
        opEC_Points(iE,iX) = opEC_Points(iE,iX) * UnitEC
      END DO
    END DO
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "UnitEC = ", UnitEC
    ! WRITE(*,*) "opEC_Points(1,1) = ", opEC_Points(1,1)
    ! WRITE(*,*) "opEC_Points(2,1) = ", opEC_Points(2,1)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"

#else

    opEC_Points = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC_Points_SG


  SUBROUTINE ComputeNeutrinoOpacities_ES_Points_SG &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iSpecies, iMoment, opES_Points )

    ! --- Elastic Scattering Opacities (Multiple D,T,Y) ---
    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES_Points(:,:)

    INTEGER  :: iX, iE, nE, nX
    INTEGER  :: k, iV
    REAL(DP) :: LogE_P(iE_B:iE_E), LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    REAL(DP), allocatable :: Points(:,:), Values(:,:)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opES_Points )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points_SG] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points_SG]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points_SG]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points_SG]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points_SG]   Y missing'
      IF ( .not. QueryOnGPU( opES_Points ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points_SG]   opES_Points missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

    nE = iE_E - iE_B + 1
    nX = iX_E - iX_B + 1

    DO iX = iX_B, iX_E
      LogD_P(iX) = LOG10( D(iX) / UnitD )
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      Y_P(iX) = Y(iX) / UnitY
    END DO

    DO iE = iE_B, iE_E
      LogE_P(iE) = LOG10( E(iE) / UnitE )
    END DO

    allocate(Points(4, nX * nE))
    allocate(Values(2, nX * nE))

    k = 0
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E
        k = k + 1
        Points(1,k) = Y_P(iX)
        Points(2,k) = LogT_P(iX)
        Points(3,k) = LogD_P(iX)
        Points(4,k) = LogE_P(iE)
      END DO
    END DO

    CALL tsgEvaluateBatch(gridIso, Points,  nX * nE, Values)
    ! iV = (iSpecies-1)*2 + iMoment
    iV = iSpecies

    k = 0
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E
        k = k + 1
        !opES_Points(iE,iX) = Values(iV,k)
        ! opES_Points(iE,iX) = 10.0d0**(Values(iV,k)) - OS_Iso(iSpecies, iMoment)
        ! opES_Points(iE,iX) = 10.0d0**(Values(iV,k) - 19.403037876305483) - OS_Iso(iSpecies, iMoment)
        opES_Points(iE,iX) = 10.0d0**(Values(iV,k) - 19.403037876305483) - 0.0d0
        ! IF ((iE == 12).and.(iX==2)) THEN
        !     WRITE(*,*) "Points at (12,2) = ", Points(:,k)
        ! END IF
        END DO
    END DO

    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28,157,60,13,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28,157,60,13+1,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28,157,60+1,13,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28,157,60+1,13+1,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28,157+1,60,13,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28,157+1,60,13+1,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28,157+1,60+1,13,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28,157+1,60+1,13+1,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28+1,157,60,13,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28+1,157,60,13+1,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28+1,157,60+1,13,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28+1,157,60+1,13+1,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28+1,157+1,60,13,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28+1,157+1,60,13+1,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28+1,157+1,60+1,13,iMoment,iSpecies)
    ! WRITE(*,*) "ISO_T around (12,2) = ", Iso_T(28+1,157+1,60+1,13+1,iMoment,iSpecies)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "opES_Points(1,1) = ", opES_Points(1,1)
    ! WRITE(*,*) "opES_Points(12,2) = ", opES_Points(12,2)
    ! WRITE(*,*) "opES_Points(11,1) = ", opES_Points(11,1)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    !
    ! CALL LogInterpolateSingleVariable_1D3D_Custom &
    !        ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
    !          OS_Iso(iSpecies,iMoment), Iso_T(:,:,:,:,iMoment,iSpecies), opES_Points, &
    !          GPU_Option = do_gpu )
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"
    ! WRITE(*,*) "opES_Points(1,1) = ", opES_Points(1,1)
    ! WRITE(*,*) "opES_Points(12,2) = ", opES_Points(12,2)
    ! WRITE(*,*) "opES_Points(11,1) = ", opES_Points(11,1)
    ! WRITE(*,*) "-------------------------------------------------------------------------------------------------"

    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E
        opES_Points(iE,iX) = opES_Points(iE,iX) * UnitES
      END DO
    END DO

#else

    opES_Points = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES_Points_SG

  SUBROUTINE ComputeNeutrinoOpacities_NES_Point_SG &
    ( iE_B, iE_E, E, D, T, Y, iS_1, iS_2, iMoment, &
      Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2, WORK1, WORK2 )

    ! --- Neutrino-Electron Scattering Opacities (Multiple D,T,Y) ---
    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D
    REAL(DP), INTENT(in)  :: T
    REAL(DP), INTENT(in)  :: Y
    INTEGER,  INTENT(in)  :: iS_1, iS_2
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In_1(:,:), Phi_Out_1(:,:)
    REAL(DP), INTENT(out) :: Phi_In_2(:,:), Phi_Out_2(:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK1(:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK2(:,:)

    REAL(DP), POINTER :: H1(:,:), H2(:,:)
    INTEGER  :: iE1, iE2, iH1, iH2, nE
    INTEGER  :: i, j, k, N
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P, LogEta_P, LogE_P(iE_B:iE_E)
    REAL(DP), allocatable :: Points(:,:), Values(:,:)

#ifdef MICROPHYSICS_WEAKLIB

    iH1 = ( iMoment - 1 ) * 2 + 1
    iH2 = ( iMoment - 1 ) * 2 + 2

    nE = iE_E - iE_B + 1

    IF ( PRESENT( WORK1 ) ) THEN
      H1 => WORK1(:,:)
    ELSE
      ALLOCATE( H1(iE_B:iE_E,iE_B:iE_E) )
    END IF
    IF ( PRESENT( WORK2 ) ) THEN
      H2 => WORK2(:,:)
    ELSE
      ALLOCATE( H2(iE_B:iE_E,iE_B:iE_E) )
    END IF

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    LogT_P = LOG10( T / UnitT )

    kT = BoltzmannConstant * T
    LogEta_P = LOG10( LogEta_P / kT / UnitEta )

    DO iE1 = iE_B, iE_E
      LogE_P(iE1) = LOG10( E(iE1) / UnitE )
    END DO

    allocate(Points(4, nE * nE))
    allocate(Values(2, nE * nE))

    k = 0
    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E
        k = k + 1
        Points(1,k) = LogEta_P
        Points(2,k) = LogT_P
        Points(3,k) = LogE_P(iE1)
        Points(4,k) = LogE_P(iE2)
      END DO
    END DO



    CALL tsgEvaluateBatch(gridNES, Points,  nE * nE, Values)



    k = 0
    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E
        k = k + 1
        ! H1(iE1,iE2) = Values(1,k)
        ! H2(iE1,iE2) = Values(2,k)
        ! H1(iE1,iE2) = 10.0d0**(Values(1,k)) - OS_NES(1,iH1)
        ! H2(iE1,iE2) = 10.0d0**(Values(2,k)) - OS_NES(1,iH2)
        H1(iE1,iE2) = 10.0d0**(Values(1,k) - 20.0d0) - 0.0d0
        H2(iE1,iE2) = 10.0d0**(Values(2,k) - 20.0d0) - 0.0d0
      END DO
    END DO


    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        kT = BoltzmannConstant * T
        DetBal = EXP( - ABS( E(iE2) - E(iE1) ) / kT )

        IF ( iE1 <= iE2 ) THEN
          Phi_Out_1(iE1,iE2) = ( C1(iS_1) * H1(iE1,iE2) + C2(iS_1) * H2(iE1,iE2) ) * UnitNES
          Phi_In_1 (iE1,iE2) = Phi_Out_1(iE1,iE2) * DetBal

          Phi_Out_2(iE1,iE2) = ( C1(iS_2) * H1(iE1,iE2) + C2(iS_2) * H2(iE1,iE2) ) * UnitNES
          Phi_In_2 (iE1,iE2) = Phi_Out_2(iE1,iE2) * DetBal
        ELSE
          Phi_In_1 (iE1,iE2) = ( C1(iS_1) * H1(iE2,iE1) + C2(iS_1) * H2(iE2,iE1) ) * UnitNES
          Phi_Out_1(iE1,iE2) = Phi_In_1(iE1,iE2) * DetBal

          Phi_In_2 (iE1,iE2) = ( C1(iS_2) * H1(iE2,iE1) + C2(iS_2) * H2(iE2,iE1) ) * UnitNES
          Phi_Out_2(iE1,iE2) = Phi_In_2(iE1,iE2) * DetBal
        END IF

      END DO
    END DO

    IF ( .NOT. PRESENT( WORK1 ) ) DEALLOCATE( H1 )
    IF ( .NOT. PRESENT( WORK2 ) ) DEALLOCATE( H2 )

#else

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E
        Phi_Out_1(iE1,iE2) = Zero
        Phi_In_1 (iE1,iE2) = Zero
        Phi_Out_2(iE1,iE2) = Zero
        Phi_In_2 (iE1,iE2) = Zero
      END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NES_Point_SG

  SUBROUTINE ComputeNeutrinoOpacities_NES_Points_SG &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iS_1, iS_2, iMoment, &
      Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2, WORK1, WORK2 )

    ! --- Neutrino-Electron Scattering Opacities (Multiple D,T,Y) ---
    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iS_1, iS_2
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In_1(:,:,:), Phi_Out_1(:,:,:)
    REAL(DP), INTENT(out) :: Phi_In_2(:,:,:), Phi_Out_2(:,:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK1(:,:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK2(:,:,:)

    REAL(DP), POINTER :: H1(:,:,:), H2(:,:,:)
    INTEGER  :: iX, iE1, iE2, iH1, iH2, nE, nX
    INTEGER  :: i, j, k, N
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E), LogE_P(iE_B:iE_E)
    REAL(DP), allocatable :: Points(:,:), Values(:,:)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( Phi_In_1, Phi_Out_1 ) &
       .AND. QueryOnGPU( Phi_In_2, Phi_Out_2 )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) &
        '[ComputeNeutrinoOpacities_NES_Points_SG] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points_SG]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points_SG]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points_SG]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points_SG]   Y missing'
      IF ( .not. QueryOnGPU( Phi_In_1 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points_SG]   Phi_In_1 missing'
      IF ( .not. QueryOnGPU( Phi_Out_1 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points_SG]   Phi_Out_1 missing'
      IF ( .not. QueryOnGPU( Phi_In_2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points_SG]   Phi_In_2 missing'
      IF ( .not. QueryOnGPU( Phi_Out_2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points_SG]   Phi_Out_2 missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

    iH1 = ( iMoment - 1 ) * 2 + 1
    iH2 = ( iMoment - 1 ) * 2 + 2

    nE = iE_E - iE_B + 1
    nX = iX_E - iX_B + 1

    IF ( PRESENT( WORK1 ) ) THEN
      H1 => WORK1(:,:,:)
    ELSE
      ALLOCATE( H1(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E) )
    END IF
    IF ( PRESENT( WORK2 ) ) THEN
      H2 => WORK2(:,:,:)
    ELSE
      ALLOCATE( H2(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E) )
    END IF

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )

      kT = BoltzmannConstant * T(iX)
      LogEta_P(iX) = LOG10( LogEta_P(iX) / kT / UnitEta )
    END DO

    DO iE1 = iE_B, iE_E
      LogE_P(iE1) = LOG10( E(iE1) / UnitE )
    END DO

    allocate(Points(4, nX * nE * nE))
    allocate(Values(2, nX * nE * nE))
    ! allocate(Values(2, nX * nE * nE))
    !
    k = 0
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          k = k + 1
          Points(1,k) = LogEta_P(iX)
          Points(2,k) = LogT_P(iX)
          Points(3,k) = LogE_P(iE1)
          Points(4,k) = LogE_P(iE2)
        END DO
      END DO
    END DO
    ! allocate(Points(4, 2))
    ! allocate(Values(4, 2))


    ! Points(1,:) = [LogEta_P(1), LogT_P(1), LogE_P(1), LogE_P(1)]
    ! Points(1,:) = LogEta_P(1:2)
    ! Points(2,:) = LogT_P(1:2)
    ! Points(3,:) = LogE_P(1:2)
    ! Points(4,:) = LogE_P(1:2)
    !
    ! WRITE(*,*) "LogEta_P = ", LogEta_P
    ! WRITE(*,*) "LogT_P = ", LogT_P
    ! WRITE(*,*) "LogE_P = ", LogE_P
    ! WRITE(*,*) "Points(:,1) = "
    ! WRITE(*,*) Points(:,1)
    ! WRITE(*,*) "Points(:,end) = "
    ! WRITE(*,*) Points(:,32*32)
    ! k = 0
    ! DO iX = iX_B, iX_E
    !   DO iE2 = iE_B, iE_E
    !     DO iE1 = iE_B, iE_E
    !       k = k + 1
    !       WRITE(*,*) Points(k,:)
    !     END DO
    !   END DO
    ! END DO

    ! N = tsgGetNumPoints(gridNES)
    ! WRITE(*,*) "N = ", N

    CALL tsgEvaluateBatch(gridNES, Points,  nX * nE * nE, Values)

    ! WRITE(*,*) "Values = ", Values
    ! k = 0
    ! DO iX = iX_B, iX_E
    !   DO iE2 = iE_B, iE_E
    !     DO iE1 = iE_B, iE_E
    !       k = k + 1
    !       WRITE(*,*) Values(k,:)
    !     END DO
    !   END DO
    ! END DO

    ! WRITE(*,*) "OS_NES(1,iH1) = "
    ! WRITE(*,*) OS_NES(1,iH1)
    ! WRITE(*,*) "OS_NES(1,iH2) = "
    ! WRITE(*,*) OS_NES(1,iH2)
    ! WRITE(*,*) "OS_NES(1,3) = "
    ! WRITE(*,*) OS_NES(1,3)
    ! WRITE(*,*) "OS_NES(1,4) = "
    ! WRITE(*,*) OS_NES(1,4)

    k = 0
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          k = k + 1
          ! H1(iE1,iE2,iX) = Values(1,k)
          ! H2(iE1,iE2,iX) = Values(2,k)
          ! H1(iE1,iE2,iX) = 10.0d0**(Values(1,k)) - OS_NES(1,iH1)
          ! H2(iE1,iE2,iX) = 10.0d0**(Values(2,k)) - OS_NES(1,iH2)
          H1(iE1,iE2,iX) = 10.0d0**(Values(1,k) - 20.0d0) - 0.0d0
          H2(iE1,iE2,iX) = 10.0d0**(Values(2,k) - 20.0d0) - 0.0d0
        END DO
      END DO
    END DO

    ! WRITE(*,*) "H1_SG = ["

    ! DO iE1 = iE_B, iE_E
    !   WRITE(*,*) H1(iE1,:,1), ";"
    ! END DO
    ! WRITE(*,*) H1(1,3,2)
    ! WRITE(*,*) H1(4,2,1)
    ! WRITE(*,*) " ];"

    ! ! --- Interpolate HI  (put result temporarily in Phi_In) ---
    !
    ! CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
    !        ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
    !          OS_NES(1,iH1), NES_AT(:,:,:,:,iH1,1), H1, &
    !          GPU_Option = do_gpu, ASYNC_Option = 1 )
    !
    ! WRITE(*,*) "H1 = ["
    ! ! DO iE1 = iE_B, iE_E
    ! !     WRITE(*,*) H1(iE1,:,1), ";"
    ! ! END DO
    ! WRITE(*,*) H1(1,3,2)
    ! WRITE(*,*) H1(4,2,1)
    ! WRITE(*,*) " ];"
    !
    ! ! --- Interpolate HII (put result temporarily in Phi_Out) ---
    !
    ! CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
    !        ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
    !          OS_NES(1,iH2), NES_AT(:,:,:,:,iH2,1), H2, &
    !          GPU_Option = do_gpu, ASYNC_Option = 1  )
    !
    !
    ! k = 0
    ! WRITE(*,*) "Points = "
    ! DO iX = iX_B, iX_B
    !  DO iE2 = iE_B, iE_E
    !    DO iE1 = iE_B, iE_E
    !      k = k + 1
    !      IF ( iE1 <= iE2 ) THEN
    !        WRITE(*,*) Points(1,k), Points(2,k), Points(3,k), Points(4,k)
    !      ENDIF
    !    END DO
    !  END DO
    ! END DO
    !
    ! k = 0
    ! WRITE(*,*) "Values = "
    ! DO iX = iX_B, iX_B
    !  DO iE2 = iE_B, iE_E
    !    DO iE1 = iE_B, iE_E
    !      k = k + 1
    !      IF ( iE1 <= iE2 ) THEN
    !        WRITE(*,*) Values(1,k), Values(2,k), LOG10(H1(iE1,iE2,iX)+OS_NES(1,iH1)), LOG10(H2(iE1,iE2,iX)+OS_NES(1,iH2))
    !      ! IF ( iE1 >= iE2 ) THEN
    !      !   WRITE(*,*) Values(1,k), Values(2,k), LOG10(H1(iE2,iE1,iX)+OS_NES(1,iH1)), LOG10(H2(iE2,iE1,iX)+OS_NES(1,iH2))
    !      END IF
    !    END DO
    !  END DO
    ! END DO


    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E

          kT = BoltzmannConstant * T(iX)
          DetBal = EXP( - ABS( E(iE2) - E(iE1) ) / kT )

          IF ( iE1 <= iE2 ) THEN
            Phi_Out_1(iE1,iE2,iX) = ( C1(iS_1) * H1(iE1,iE2,iX) + C2(iS_1) * H2(iE1,iE2,iX) ) * UnitNES
            Phi_In_1 (iE1,iE2,iX) = Phi_Out_1(iE1,iE2,iX) * DetBal

            Phi_Out_2(iE1,iE2,iX) = ( C1(iS_2) * H1(iE1,iE2,iX) + C2(iS_2) * H2(iE1,iE2,iX) ) * UnitNES
            Phi_In_2 (iE1,iE2,iX) = Phi_Out_2(iE1,iE2,iX) * DetBal
          ELSE
            Phi_In_1 (iE1,iE2,iX) = ( C1(iS_1) * H1(iE2,iE1,iX) + C2(iS_1) * H2(iE2,iE1,iX) ) * UnitNES
            Phi_Out_1(iE1,iE2,iX) = Phi_In_1(iE1,iE2,iX) * DetBal

            Phi_In_2 (iE1,iE2,iX) = ( C1(iS_2) * H1(iE2,iE1,iX) + C2(iS_2) * H2(iE2,iE1,iX) ) * UnitNES
            Phi_Out_2(iE1,iE2,iX) = Phi_In_2(iE1,iE2,iX) * DetBal
          END IF

        END DO
      END DO
    END DO

    IF ( .NOT. PRESENT( WORK1 ) ) DEALLOCATE( H1 )
    IF ( .NOT. PRESENT( WORK2 ) ) DEALLOCATE( H2 )

#else

    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          Phi_Out_1(iE1,iE2,iX) = Zero
          Phi_In_1 (iE1,iE2,iX) = Zero
          Phi_Out_2(iE1,iE2,iX) = Zero
          Phi_In_2 (iE1,iE2,iX) = Zero
        END DO
      END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NES_Points_SG


  SUBROUTINE ComputeNeutrinoOpacitiesRates_NES_Points &
    ( iE_B, iE_E, iX_B, iX_E, W2, J, Phi_In, Phi_Out, Eta, Chi )

    ! --- Neutrino-Electron Scattering Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2     (:)
    REAL(DP), INTENT(in)  :: J      (:,:)
    REAL(DP), INTENT(in)  :: Phi_In (:,:,:)
    REAL(DP), INTENT(in)  :: Phi_Out(:,:,:)
    REAL(DP), INTENT(out) :: Eta    (:,:)
    REAL(DP), INTENT(out) :: Chi    (:,:)

    REAL(DP) :: fEta(iE_B:iE_E)
    REAL(DP) :: fChi(iE_B:iE_E)
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE, iE1, iE2, nX, nE
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( W2 ) &
       .AND. QueryOnGPU( J, Eta, Chi ) &
       .AND. QueryOnGPU( Phi_In, Phi_Out )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_NES_Points] Data not present on device'
      IF ( .not. QueryOnGPU( W2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_NES_Points]   W2 missing'
      IF ( .not. QueryOnGPU( J ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_NES_Points]   J missing'
      IF ( .not. QueryOnGPU( Eta ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_NES_Points]   Eta missing'
      IF ( .not. QueryOnGPU( Chi ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_NES_Points]   Chi missing'
      IF ( .not. QueryOnGPU( Phi_In ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_NES_Points]   Phi_In missing'
      IF ( .not. QueryOnGPU( Phi_Out ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_NES_Points]   Phi_Out missing'
    END IF
#endif

    nX = iX_E - iX_B + 1
    nE = iE_E - iE_B + 1

    IF ( do_gpu ) THEN

      DO iX = iX_B, iX_E
        DO iE2 = iE_B, iE_E

          SUM1 = Zero
          DO iE1 = iE_B, iE_E
            SUM1 = SUM1 + Phi_In (iE1,iE2,iX) * W2(iE1) * J(iE1,iX)
          END DO
          Eta(iE2,iX) = SUM1

          SUM2 = Zero
          DO iE1 = iE_B, iE_E
            SUM2 = SUM2 + Phi_Out(iE1,iE2,iX) * W2(iE1) * ( One - J(iE1,iX) )
          END DO
          Chi(iE2,iX) = SUM1 + SUM2

        END DO
      END DO

    ELSE

      DO iX = iX_B, iX_E

        DO iE = iE_B, iE_E
          fEta(iE) = W2(iE) * J(iE,iX)
          fChi(iE) = W2(iE) * ( One - J(iE,iX) )
        END DO

        ! --- Emissivity ---

        CALL MatrixVectorMultiply &
          ( 'T', nE, nE, One, Phi_In(:,:,iX), nE, &
            fEta(iE_B), 1, Zero, Eta(:,iX), 1 )

        DO iE = iE_B, iE_E
          Chi(iE,iX) = Eta(iE,iX)
        END DO

        ! --- Absorptivity ---

        CALL MatrixVectorMultiply &
          ( 'T', nE, nE, One, Phi_Out(:,:,iX), nE, &
            fChi(iE_B), 1, One, Chi(:,iX), 1 )

      END DO

    END IF

  END SUBROUTINE ComputeNeutrinoOpacitiesRates_NES_Points

  SUBROUTINE ComputeNeutrinoOpacities_Pair_Point_SG &
    ( iE_B, iE_E, E, D, T, Y, iS_1, iS_2, iMoment, &
      Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2, WORK1, WORK2 )

    ! --- Pair Opacities (Multiple D,T,Y) ---
    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D
    REAL(DP), INTENT(in)  :: T
    REAL(DP), INTENT(in)  :: Y
    INTEGER,  INTENT(in)  :: iS_1, iS_2
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In_1(:,:), Phi_Out_1(:,:)
    REAL(DP), INTENT(out) :: Phi_In_2(:,:), Phi_Out_2(:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK1(:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK2(:,:)

    REAL(DP), POINTER :: J1(:,:), J2(:,:)
    INTEGER  :: iE1, iE2, iJ1, iJ2, nE
    INTEGER  :: i, j, k
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P, LogEta_P, LogE_P(iE_B:iE_E)
    REAL(DP), allocatable :: Points(:,:), Values(:,:)


#ifdef MICROPHYSICS_WEAKLIB

    iJ1 = ( iMoment - 1 ) * 2 + 1
    iJ2 = ( iMoment - 1 ) * 2 + 2

    nE = iE_E - iE_B + 1

    IF ( PRESENT( WORK1 ) ) THEN
      J1 => WORK1(:,:)
    ELSE
      ALLOCATE( J1(iE_B:iE_E,iE_B:iE_E) )
    END IF
    IF ( PRESENT( WORK2 ) ) THEN
      J2 => WORK2(:,:)
    ELSE
      ALLOCATE( J2(iE_B:iE_E,iE_B:iE_E) )
    END IF

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    LogT_P = LOG10( T / UnitT )
    LogEta_P = LOG10( LogEta_P / ( BoltzmannConstant * T ) / UnitEta )

    DO iE1 = iE_B, iE_E
      LogE_P(iE1) = LOG10( E(iE1) / UnitE )
    END DO

    allocate(Points(4, nE * nE))
    allocate(Values(2, nE * nE))
    ! allocate(Values(2, nX * nE * nE))
    !
    k = 0
    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E
        k = k + 1
        Points(1,k) = LogEta_P
        Points(2,k) = LogT_P
        Points(3,k) = LogE_P(iE1)
        Points(4,k) = LogE_P(iE2)
      END DO
    END DO

    CALL tsgEvaluateBatch(gridPair, Points,  nE * nE, Values)

    k = 0
    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E
        k = k + 1
        ! H1(iE1,iE2) = Values(1,k)
        ! H2(iE1,iE2) = Values(2,k)
        J1(iE1,iE2) = 10.0d0**(Values(1,k)) - OS_Pair(1,iJ1)
        J2(iE1,iE2) = 10.0d0**(Values(2,k)) - OS_Pair(1,iJ2)
        ! J1(iE1,iE2) = 10.0d0**(Values(1,k) - 20.0d0) - 9.391433233644253d-23
        ! J2(iE1,iE2) = 10.0d0**(Values(2,k) - 20.0d0) - 9.391433233644254d-23
      END DO
    END DO


    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        kT = BoltzmannConstant * T
        DetBal = EXP( - ABS( E(iE1) + E(iE2) ) / kT )

        IF ( iE1 <= iE2 ) THEN
          Phi_Out_1(iE1,iE2) = ( C1(iS_1) * J1(iE1,iE2) + C2(iS_1) * J2(iE1,iE2) ) * UnitPair
          Phi_Out_2(iE1,iE2) = ( C1(iS_2) * J1(iE1,iE2) + C2(iS_2) * J2(iE1,iE2) ) * UnitPair
        ELSE
          Phi_Out_1(iE1,iE2) = ( C1(iS_1) * J2(iE2,iE1) + C2(iS_1) * J1(iE2,iE1) ) * UnitPair
          Phi_Out_2(iE1,iE2) = ( C1(iS_2) * J2(iE2,iE1) + C2(iS_2) * J1(iE2,iE1) ) * UnitPair
        END IF
        Phi_In_1(iE1,iE2) = Phi_Out_1(iE1,iE2) * DetBal
        Phi_In_2(iE1,iE2) = Phi_Out_2(iE1,iE2) * DetBal

      END DO
    END DO

    IF ( .NOT. PRESENT( WORK1 ) ) DEALLOCATE( J1 )
    IF ( .NOT. PRESENT( WORK2 ) ) DEALLOCATE( J2 )

#else

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E
        Phi_Out_1(iE1,iE2) = Zero
        Phi_In_1 (iE1,iE2) = Zero
        Phi_Out_2(iE1,iE2) = Zero
        Phi_In_2 (iE1,iE2) = Zero
      END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Pair_Point_SG


  SUBROUTINE ComputeNeutrinoOpacities_Pair_Points_SG &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iS_1, iS_2, iMoment, &
      Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2, WORK1, WORK2 )

    ! --- Pair Opacities (Multiple D,T,Y) ---
    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iS_1, iS_2
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In_1(:,:,:), Phi_Out_1(:,:,:)
    REAL(DP), INTENT(out) :: Phi_In_2(:,:,:), Phi_Out_2(:,:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK1(:,:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK2(:,:,:)

    REAL(DP), POINTER :: J1(:,:,:), J2(:,:,:)
    INTEGER  :: iX, iE1, iE2, iJ1, iJ2, nE, nX
    INTEGER  :: i, j, k
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E), LogE_P(iE_B:iE_E)
    REAL(DP), allocatable :: Points(:,:), Values(:,:)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( Phi_In_1, Phi_Out_1 ) &
       .AND. QueryOnGPU( Phi_In_2, Phi_Out_2 )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG]   Y missing'
      IF ( .not. QueryOnGPU( Phi_In_1 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG]   Phi_In_1 missing'
      IF ( .not. QueryOnGPU( Phi_Out_1 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG]   Phi_Out_1 missing'
      IF ( .not. QueryOnGPU( Phi_In_2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG]   Phi_In_2 missing'
      IF ( .not. QueryOnGPU( Phi_Out_2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points_SG]   Phi_Out_2 missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

    iJ1 = ( iMoment - 1 ) * 2 + 1
    iJ2 = ( iMoment - 1 ) * 2 + 2

    nE = iE_E - iE_B + 1
    nX = iX_E - iX_B + 1

    IF ( PRESENT( WORK1 ) ) THEN
      J1 => WORK1(:,:,:)
    ELSE
      ALLOCATE( J1(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E) )
    END IF
    IF ( PRESENT( WORK2 ) ) THEN
      J2 => WORK2(:,:,:)
    ELSE
      ALLOCATE( J2(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E) )
    END IF

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

    DO iE1 = iE_B, iE_E
      LogE_P(iE1) = LOG10( E(iE1) / UnitE )
    END DO

    allocate(Points(4, nX * nE * nE))
    allocate(Values(2, nX * nE * nE))
    ! allocate(Values(2, nX * nE * nE))
    !
    k = 0
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          k = k + 1
          Points(1,k) = LogEta_P(iX)
          Points(2,k) = LogT_P(iX)
          Points(3,k) = LogE_P(iE1)
          Points(4,k) = LogE_P(iE2)
        END DO
      END DO
    END DO

    CALL tsgEvaluateBatch(gridPair, Points,  nX * nE * nE, Values)

    k = 0
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          k = k + 1
          ! H1(iE1,iE2,iX) = Values(1,k)
          ! H2(iE1,iE2,iX) = Values(2,k)
          ! J1(iE1,iE2,iX) = 10.0d0**(Values(1,k)) - OS_Pair(1,iJ1)
          ! J2(iE1,iE2,iX) = 10.0d0**(Values(2,k)) - OS_Pair(1,iJ2)
          J1(iE1,iE2,iX) = 10.0d0**(Values(1,k) - 20.0d0) - 9.391433233644253d-23
          J2(iE1,iE2,iX) = 10.0d0**(Values(2,k) - 20.0d0) - 9.391433233644254d-23
        END DO
      END DO
    END DO

    ! WRITE(*,*) "J1_SG = ["
    ! WRITE(*,*) J1(10,3,2)
    ! WRITE(*,*) J1(4,20,1)
    ! WRITE(*,*) " ];"

    ! ! --- Interpolate JI  (put result temporarily in Phi_In) ---
    !
    ! CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
    !        ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
    !          OS_Pair(1,iJ1), Pair_AT(:,:,:,:,iJ1,1), J1, &
    !          GPU_Option = do_gpu, ASYNC_Option = 1 )
    ! WRITE(*,*) "J1 = ["
    ! WRITE(*,*) J1(10,3,2)
    ! WRITE(*,*) J1(4,20,1)
    ! WRITE(*,*) " ];"
    !
    ! ! --- Interpolate JII (put result temporarily in Phi_Out) ---
    !
    ! CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
    !        ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
    !          OS_Pair(1,iJ2), Pair_AT(:,:,:,:,iJ2,1), J2, &
    !          GPU_Option = do_gpu, ASYNC_Option = 1 )

    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E

          kT = BoltzmannConstant * T(iX)
          DetBal = EXP( - ABS( E(iE1) + E(iE2) ) / kT )

          IF ( iE1 <= iE2 ) THEN
            Phi_Out_1(iE1,iE2,iX) = ( C1(iS_1) * J1(iE1,iE2,iX) + C2(iS_1) * J2(iE1,iE2,iX) ) * UnitPair
            Phi_Out_2(iE1,iE2,iX) = ( C1(iS_2) * J1(iE1,iE2,iX) + C2(iS_2) * J2(iE1,iE2,iX) ) * UnitPair
          ELSE
            Phi_Out_1(iE1,iE2,iX) = ( C1(iS_1) * J2(iE2,iE1,iX) + C2(iS_1) * J1(iE2,iE1,iX) ) * UnitPair
            Phi_Out_2(iE1,iE2,iX) = ( C1(iS_2) * J2(iE2,iE1,iX) + C2(iS_2) * J1(iE2,iE1,iX) ) * UnitPair
          END IF
          Phi_In_1(iE1,iE2,iX) = Phi_Out_1(iE1,iE2,iX) * DetBal
          Phi_In_2(iE1,iE2,iX) = Phi_Out_2(iE1,iE2,iX) * DetBal

        END DO
      END DO
    END DO

    IF ( .NOT. PRESENT( WORK1 ) ) DEALLOCATE( J1 )
    IF ( .NOT. PRESENT( WORK2 ) ) DEALLOCATE( J2 )

#else

    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          Phi_Out_1(iE1,iE2,iX) = Zero
          Phi_In_1 (iE1,iE2,iX) = Zero
          Phi_Out_2(iE1,iE2,iX) = Zero
          Phi_In_2 (iE1,iE2,iX) = Zero
        END DO
      END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Pair_Points_SG


  SUBROUTINE ComputeNeutrinoOpacitiesRates_Pair_Points &
    ( iE_B, iE_E, iX_B, iX_E, W2, J, Phi_In, Phi_Out, Eta, Chi )

    ! --- Pair Rates (Multiple J) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: W2     (:)
    REAL(DP), INTENT(in)  :: J      (:,:)
    REAL(DP), INTENT(in)  :: Phi_In (:,:,:)
    REAL(DP), INTENT(in)  :: Phi_Out(:,:,:)
    REAL(DP), INTENT(out) :: Eta    (:,:)
    REAL(DP), INTENT(out) :: Chi    (:,:)

    REAL(DP) :: fEta(iE_B:iE_E)
    REAL(DP) :: fChi(iE_B:iE_E)
    REAL(DP) :: SUM1, SUM2
    INTEGER  :: iX, iE, iE1, iE2, nX, nE
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( W2 ) &
       .AND. QueryOnGPU( J, Eta, Chi ) &
       .AND. QueryOnGPU( Phi_In, Phi_Out )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_Pair_Points] Data not present on device'
      IF ( .not. QueryOnGPU( W2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_Pair_Points]   W2 missing'
      IF ( .not. QueryOnGPU( J ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_Pair_Points]   J missing'
      IF ( .not. QueryOnGPU( Eta ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_Pair_Points]   Eta missing'
      IF ( .not. QueryOnGPU( Chi ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_Pair_Points]   Chi missing'
      IF ( .not. QueryOnGPU( Phi_In ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_Pair_Points]   Phi_In missing'
      IF ( .not. QueryOnGPU( Phi_Out ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacitiesRates_Pair_Points]   Phi_Out missing'
    END IF
#endif

    nX = iX_E - iX_B + 1
    nE = iE_E - iE_B + 1

    IF ( do_gpu ) THEN

      DO iX = iX_B, iX_E
        DO iE2 = iE_B, iE_E

          SUM1 = Zero
          DO iE1 = iE_B, iE_E
            SUM1 = SUM1 + Phi_In (iE1,iE2,iX) * W2(iE1) * ( One - J(iE1,iX) )
          END DO
          Eta(iE2,iX) = SUM1

          SUM2 = Zero
          DO iE1 = iE_B, iE_E
            SUM2 = SUM2 + Phi_Out(iE1,iE2,iX) * W2(iE1) * J(iE1,iX)
          END DO
          Chi(iE2,iX) = SUM1 + SUM2

        END DO
      END DO

    ELSE

      DO iX = iX_B, iX_E

        DO iE = iE_B, iE_E
          fEta(iE) = W2(iE) * ( One - J(iE,iX) )
          fChi(iE) = W2(iE) * J(iE,iX)
        END DO

        ! --- Emissivity ---

        CALL MatrixVectorMultiply &
          ( 'T', nE, nE, One, Phi_In(:,:,iX), nE, &
            fEta(iE_B), 1, Zero, Eta(:,iX), 1 )

        DO iE = iE_B, iE_E
          Chi(iE,iX) = Eta(iE,iX)
        END DO

        ! --- Absorptivity ---

        CALL MatrixVectorMultiply &
          ( 'T', nE, nE, One, Phi_Out(:,:,iX), nE, &
            fChi(iE_B), 1, One, Chi(:,iX), 1 )

      END DO

    END IF

  END SUBROUTINE ComputeNeutrinoOpacitiesRates_Pair_Points



  FUNCTION FermiDirac_Scalar( E, Mu, kT ) RESULT( FermiDirac )

    REAL(DP), INTENT(in) :: E, Mu, kT
    REAL(DP) :: FermiDirac

    REAL(DP) :: Exponent

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FermiDirac &
      = One / ( EXP( Exponent ) + One )

    RETURN
  END FUNCTION FermiDirac_Scalar


  FUNCTION FermiDirac_Vector( E, Mu, kT ) RESULT( FermiDirac )

    REAL(DP), INTENT(in) :: E(:), Mu, kT
    REAL(DP) :: FermiDirac(SIZE(E))

    REAL(DP) :: Exponent(SIZE(E))

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FermiDirac &
      = One / ( EXP( Exponent ) + One )

    RETURN
  END FUNCTION FermiDirac_Vector


  FUNCTION dFermiDiracdT_Scalar( E, Mu, kT, dMudT, T ) RESULT( dFermiDiracdT )

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudT, T
    REAL(DP) :: dFermiDiracdT

    REAL(DP) :: Exponent, FD

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdT &
      = ( FD * EXP( Exponent ) ) * FD * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT_Scalar


  FUNCTION dFermiDiracdT_Vector( E, Mu, kT, dMudT, T ) RESULT( dFermiDiracdT )

    REAL(DP), INTENT(in) :: E(:), Mu, kT, dMudT, T
    REAL(DP) :: dFermiDiracdT(SIZE(E))

    REAL(DP) :: Exponent(SIZE(E)), FD(SIZE(E))

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdT &
      = ( FD * EXP( Exponent ) ) * FD * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT_Vector


  FUNCTION dFermiDiracdY_Scalar( E, Mu, kT, dMudY, T ) RESULT( dFermiDiracdY )

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudY, T
    REAL(DP) :: dFermiDiracdY

    REAL(DP) :: Exponent, FD

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdY &
      = ( FD * EXP( Exponent ) ) * FD * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY_Scalar


  FUNCTION dFermiDiracdY_Vector( E, Mu, kT, dMudY, T ) RESULT( dFermiDiracdY )

    REAL(DP), INTENT(in) :: E(:), Mu, kT, dMudY, T
    REAL(DP) :: dFermiDiracdY(SIZE(E))

    REAL(DP) :: Exponent(SIZE(E)), FD(SIZE(E))

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdY &
      = ( FD * EXP( Exponent ) ) * FD * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY_Vector



  ! SUBROUTINE SGEvaluationBatch &
  !   (grid, Points, Values)
  !
  !   ! --- Wrapper for tsgEvaluationBatch ---
  !   ! --- Change the dimension of inputs ---
  !   type(TasmanianSparseGrid), INTENT(in) :: grid
  !   REAL(DP), INTENT(in)  :: Points(:,:)
  !   REAL(DP), INTENT(out) :: Values(:,:)
  !   INTEGER  :: iNumX, iNumDim, iNumOutputs
  !   REAL(DP), allocatable :: x(:), y(:)
  !   INTEGER :: inputshape(2)
  !
  !   inputshape = SHAPE(Points)
  !   iNumX = inputshape(1)
  !   iNumDim = inputshape(2)
  !   iNumOutputs = tsgGetNumOutputs(grid)
  !
  !   x = RESHAPE(Points, [iNumX*iNumDim])
  !   ALLOCATE(y(iNumX*iNumOutputs))
  !
  !   call tsgEvaluateBatch(grid, x, iNumX, y)
  !
  !   Values = RESHAPE(y, [iNumX, iNumOutputs])
  !
  ! END SUBROUTINE SGEvaluationBatch
END MODULE NeutrinoOpacitiesSparseComputationModule
