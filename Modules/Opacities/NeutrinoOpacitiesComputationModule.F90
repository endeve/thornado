#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_OPACITY
#endif

MODULE NeutrinoOpacitiesComputationModule

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

  PUBLIC :: ComputeNeutrinoOpacities
  PUBLIC :: ComputeNeutrinoOpacities_EC_Point
  PUBLIC :: ComputeNeutrinoOpacities_EC_Points
  PUBLIC :: ComputeNeutrinoOpacities_EC_Vector
  PUBLIC :: ComputeNeutrinoOpacities_ES_Point
  PUBLIC :: ComputeNeutrinoOpacities_ES_Points
  PUBLIC :: ComputeNeutrinoOpacities_ES_Vector
  PUBLIC :: ComputeNeutrinoOpacities_NES_Point
  PUBLIC :: ComputeNeutrinoOpacities_NES_Points
  PUBLIC :: ComputeNeutrinoOpacitiesAndDerivatives_NES_Point
  PUBLIC :: ComputeNeutrinoOpacitiesRates_NES_Points
  PUBLIC :: ComputeNeutrinoOpacities_Pair_Point
  PUBLIC :: ComputeNeutrinoOpacities_Pair_Points
  PUBLIC :: ComputeNeutrinoOpacitiesRates_Pair_Points
  PUBLIC :: ComputeNeutrinoOpacitiesAndDerivatives_Pair_Point
  PUBLIC :: ComputeEquilibriumDistributions_Point
  PUBLIC :: ComputeEquilibriumDistributions_Points
  PUBLIC :: ComputeEquilibriumDistributionAndDerivatives_Point
  PUBLIC :: ComputeEquilibriumDistributionAndDerivatives_Points
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

  INTERFACE ComputeEquilibriumDistributions_Point
    MODULE PROCEDURE ComputeEquilibriumDistributions_Point_1
    MODULE PROCEDURE ComputeEquilibriumDistributions_Point_2
    MODULE PROCEDURE ComputeEquilibriumDistributions_Point_3
  END INTERFACE

  INTERFACE ComputeEquilibriumDistributions_Points
    MODULE PROCEDURE ComputeEquilibriumDistributions_Points_1
    MODULE PROCEDURE ComputeEquilibriumDistributions_Points_2
  END INTERFACE

  INTERFACE ComputeNeutrinoOpacities_NES_Point
    MODULE PROCEDURE ComputeNeutrinoOpacities_NES_Point_1
    MODULE PROCEDURE ComputeNeutrinoOpacities_NES_Point_2
  END INTERFACE

  INTERFACE ComputeNeutrinoOpacities_NES_Points
    MODULE PROCEDURE ComputeNeutrinoOpacities_NES_Points_1
    MODULE PROCEDURE ComputeNeutrinoOpacities_NES_Points_2
  END INTERFACE

  INTERFACE ComputeNeutrinoOpacities_Pair_Point
    MODULE PROCEDURE ComputeNeutrinoOpacities_Pair_Point_1
    MODULE PROCEDURE ComputeNeutrinoOpacities_Pair_Point_2
  END INTERFACE

  INTERFACE ComputeNeutrinoOpacities_Pair_Points
    MODULE PROCEDURE ComputeNeutrinoOpacities_Pair_Points_1
    MODULE PROCEDURE ComputeNeutrinoOpacities_Pair_Points_2
  END INTERFACE

  INTERFACE ComputeNeutrinoOpacity_Point
    MODULE PROCEDURE ComputeNeutrinoOpacity_E_D_T_Y_Point
    MODULE PROCEDURE ComputeNeutrinoOpacity_E_E_T_Eta_Point
  END INTERFACE

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


  SUBROUTINE ComputeNeutrinoOpacities( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      T(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      Y(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):)

    INTEGER  :: iS
    REAL(DP) :: wTime

    DO iS = 1, nSpecies

      CALL ComputeEquilibriumDistributions &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y, iS )

      wTime = MPI_WTIME( )

      CALL ComputeNeutrinoOpacities_EC &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y, iS )

      wTime = MPI_WTIME( ) - wTime

      PRINT*
      PRINT*, "  ComputeNeutrinoOpacities: ", iS
      PRINT*
      PRINT*, "    EC: ", wTime

      wTime = MPI_WTIME( )

      CALL ComputeNeutrinoOpacities_ES &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y, iS )

      wTime = MPI_WTIME( ) - wTime

      PRINT*, "    ES: ", wTime

    END DO

  END SUBROUTINE ComputeNeutrinoOpacities


  SUBROUTINE ComputeEquilibriumDistributions &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y, iSpecies )

    ! --- Equilibrium Neutrino Distributions ---

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      T(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      Y(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):)
    INTEGER,  INTENT(in) :: &
      iSpecies

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, nZ(4)
    INTEGER  :: iOS_X, iOS_E, iNodeX, iNodeE
    REAL(DP) :: E(nDOFE), Me(nDOFX), Mp(nDOFX), Mn(nDOFX), Mnu(nDOFX)

    nZ = iZ_E0 - iZ_B0 + 1

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      ! --- Compute Chemical Potentials ---

      CALL ComputeElectronChemicalPotential_TABLE &
             ( D(:,iZ2,iZ3,iZ4), T(:,iZ2,iZ3,iZ4), Y(:,iZ2,iZ3,iZ4), Me )

      CALL ComputeProtonChemicalPotential_TABLE &
             ( D(:,iZ2,iZ3,iZ4), T(:,iZ2,iZ3,iZ4), Y(:,iZ2,iZ3,iZ4), Mp )

      CALL ComputeNeutronChemicalPotential_TABLE &
             ( D(:,iZ2,iZ3,iZ4), T(:,iZ2,iZ3,iZ4), Y(:,iZ2,iZ3,iZ4), Mn )

      IF( iSpecies .EQ. iNuE )THEN

        ! --- Electron Neutrinos ---

        Mnu = Me + Mp - Mn

      ELSEIF( iSpecies .EQ. iNuE_Bar )THEN

        ! --- Electron Antineutrinos ---

        Mnu = Mn - Me - Mp

      ELSE

        WRITE(*,*)
        WRITE(*,'(A4,A)') '', 'ERROR (ComputeEquilibriumDistributions)'
        WRITE(*,'(A4,A,I2.2)') '', 'Invalid Species: ', iSpecies
        WRITE(*,*)
        STOP

      END IF

      ! --- Offset (Position) ---

      iOS_X = ( (iZ4-1)*nZ(3)*nZ(2) + (iZ3-1)*nZ(2) + (iZ2-1) ) * nDOFX

      DO iZ1 = iZ_B0(1), iZ_E0(1)

        ! --- Offset (Energy) ---

        iOS_E = (iZ1-1) * nDOFE

        ! --- Energy Coordinates ---

        E = MeshE % Center(iZ1) + MeshE % Width(iZ1) * MeshE % Nodes(:)

        DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE

          f_EQ(iOS_E+iNodeE,iSpecies,iOS_X+iNodeX) &
            = FermiDirac( E(iNodeE), Mnu(iNodeX), T(iNodeX,iZ2,iZ3,iZ4) )

        END DO
        END DO

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeEquilibriumDistributions


  SUBROUTINE ComputeEquilibriumDistributions_Point_1 &
    ( iE_B, iE_E, E, D, T, Y, f0, iSpecies )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Equilibrium Neutrino Distributions (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: f0(:)
    INTEGER,  INTENT(in)  :: iSpecies

    INTEGER  :: iE
    REAL(DP) :: Me, Mp, Mn, Mnu, kT

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    kT = BoltzmannConstant * T

    Mnu = ( Me + Mp ) - Mn

    IF ( iSpecies == iNuE ) THEN
      f0(:) = FermiDirac( E(:), +Mnu, kT )
    ELSE IF ( iSpecies == iNuE_Bar ) THEN
      f0(:) = FermiDirac( E(:), -Mnu, kT )
    ELSE
      f0(:) = FermiDirac( E(:), Zero, kT )
    END IF

  END SUBROUTINE ComputeEquilibriumDistributions_Point_1


  SUBROUTINE ComputeEquilibriumDistributions_Point_2 &
    ( iE_B, iE_E, E, D, T, Y, f0_1, f0_2, iS_1, iS_2 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Equilibrium Neutrino Distributions (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: f0_1(:), f0_2(:)
    INTEGER,  INTENT(in)  :: iS_1, iS_2

    INTEGER  :: iE
    REAL(DP) :: Me, Mp, Mn, Mnu, kT

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    kT = BoltzmannConstant * T

    Mnu = ( Me + Mp ) - Mn

    IF ( iS_1 == iNuE ) THEN
      f0_1(:) = FermiDirac( E(:), +Mnu, kT )
    ELSE IF ( iS_1 == iNuE_Bar ) THEN
      f0_1(:) = FermiDirac( E(:), -Mnu, kT )
    ELSE
      f0_1(:) = FermiDirac( E(:), Zero, kT )
    END IF

    IF ( iS_2 == iNuE ) THEN
      f0_2(:) = FermiDirac( E(:), +Mnu, kT )
    ELSE IF ( iS_2 == iNuE_Bar ) THEN
      f0_2(:) = FermiDirac( E(:), -Mnu, kT )
    ELSE
      f0_2(:) = FermiDirac( E(:), Zero, kT )
    END IF

  END SUBROUTINE ComputeEquilibriumDistributions_Point_2


  SUBROUTINE ComputeEquilibriumDistributions_Point_3 &
    ( E, D, T, Y, f0, iSpecies )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Equilibrium Neutrino Distributions (Single E,D,T,Y) ---

    REAL(DP), INTENT(in)  :: E
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: f0
    INTEGER,  INTENT(in)  :: iSpecies

    INTEGER  :: iE
    REAL(DP) :: Me, Mp, Mn, Mnu, kT

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    kT = BoltzmannConstant * T

    Mnu = ( Me + Mp ) - Mn

    IF ( iSpecies == iNuE ) THEN
      f0 = FermiDirac( E, +Mnu, kT )
    ELSE IF ( iSpecies == iNuE_Bar ) THEN
      f0 = FermiDirac( E, -Mnu, kT )
    ELSE
      f0 = FermiDirac( E, Zero, kT )
    END IF

  END SUBROUTINE ComputeEquilibriumDistributions_Point_3


  SUBROUTINE ComputeEquilibriumDistributions_Points_1 &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, f0, iSpecies )

    ! --- Equilibrium Neutrino Distributions (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    REAL(DP), INTENT(out) :: f0(:,:)
    INTEGER,  INTENT(in)  :: iSpecies

    INTEGER  :: iX, iE
    REAL(DP) :: Me(iX_B:iX_E), Mp(iX_B:iX_E), Mn(iX_B:iX_E)
    REAL(DP) :: Mnu, kT
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( f0 )
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
      IF ( .not. QueryOnGPU( f0 ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributions_Points]   f0 missing'
    END IF
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: Me, Mp, Mn )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( Me, Mp, Mn )
#endif

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( Mnu, kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( Mnu, kT ) &
    !$ACC PRESENT( Me, Mp, Mn, E, T, f0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Mnu, kT ) &
#endif
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E

        kT = BoltzmannConstant * T(iX)

        Mnu = ( Me(iX) + Mp(iX) ) - Mn(iX)

        IF ( iSpecies == iNuE ) THEN
          f0(iE,iX) = FermiDirac( E(iE), +Mnu, kT )
        ELSE IF ( iSpecies == iNuE_Bar ) THEN
          f0(iE,iX) = FermiDirac( E(iE), -Mnu, kT )
        ELSE
          f0(iE,iX) = FermiDirac( E(iE), Zero, kT )
        END IF

      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: Me, Mp, Mn )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( Me, Mp, Mn )
#endif

  END SUBROUTINE ComputeEquilibriumDistributions_Points_1


  SUBROUTINE ComputeEquilibriumDistributions_Points_2 &
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: Me, Mp, Mn )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( Me, Mp, Mn )
#endif

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( Mnu, kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( Mnu, kT ) &
    !$ACC PRESENT( Me, Mp, Mn, E, T, f0_1, f0_2 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( Mnu, kT ) &
#endif
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: Me, Mp, Mn )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( Me, Mp, Mn )
#endif

  END SUBROUTINE ComputeEquilibriumDistributions_Points_2


  SUBROUTINE ComputeEquilibriumDistributionAndDerivatives_Point &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, f0, df0dY, df0dU )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Equilibrium Neutrino Distribution and Derivatives (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    REAL(DP), INTENT(out) :: f0   (iE_B:iE_E)
    REAL(DP), INTENT(out) :: df0dY(iE_B:iE_E)
    REAL(DP), INTENT(out) :: df0dU(iE_B:iE_E)

    INTEGER  :: iE
    REAL(DP) :: kT, Exponent
    REAL(DP) :: Me, dMedD, dMedT, dMedY
    REAL(DP) :: Mp, dMpdD, dMpdT, dMpdY
    REAL(DP) :: Mn, dMndD, dMndT, dMndY
    REAL(DP) :: M,  dMdD,  dMdT,  dMdY
    REAL(DP) :: U,  dUdD,  dUdT,  dUdY
    REAL(DP) :: df0dT_Y(iE_B:iE_E)
    REAL(DP) :: df0dY_T(iE_B:iE_E)

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me, dMedD, dMedT, dMedY )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp, dMpdD, dMpdT, dMpdY )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn, dMndD, dMndT, dMndY )

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( D, T, Y, U,  dUdD,  dUdT,  dUdY )

    IF( iSpecies .EQ. iNuE )THEN

      ! --- Electron Neutrinos ---

      M = ( Me + Mp ) - Mn

      dMdT = dMedT + dMpdT - dMndT
      dMdY = dMedY + dMpdY - dMndY

    ELSEIF( iSpecies .EQ. iNuE_Bar )THEN

      ! --- Electron Antineutrino ---

      M = Mn - ( Me + Mp )

      dMdT = dMndT - ( dMedT + dMpdT )
      dMdY = dMndY - ( dMedY + dMpdY )

    ELSE

      M = Zero

    END IF

    kT = BoltzmannConstant * T

    DO iE = iE_B, iE_E

      f0(iE) = FermiDirac( E(iE), M, kT )

      df0dT_Y(iE) = dFermiDiracdT( E(iE), M, kT, dMdT, T ) ! Constant T
      df0dY_T(iE) = dFermiDiracdY( E(iE), M, kT, dMdY, T ) ! Constant Y

    END DO

    df0dU = df0dT_Y / dUdT
    df0dY = df0dY_T - df0dT_Y * dUdY / dUdT

  END SUBROUTINE ComputeEquilibriumDistributionAndDerivatives_Point


  SUBROUTINE ComputeEquilibriumDistributionAndDerivatives_Points &
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
      WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points]   Y missing'
      IF ( .not. QueryOnGPU( f0 ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points]   f0 missing'
      IF ( .not. QueryOnGPU( df0dY ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points]   df0dY missing'
      IF ( .not. QueryOnGPU( df0dU ) ) &
        WRITE(*,*) '[ComputeEquilibriumDistributionAndDerivatives_Points]   df0dU missing'
    END IF
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: Me,  dMedT,  dMedY, &
    !$OMP             Mp,  dMpdT,  dMpdY, &
    !$OMP             Mn,  dMndT,  dMndY, &
    !$OMP             Mnu, dMnudT, dMnudY, &
    !$OMP             U,   dUdT,   dUdY, dUdD )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( Me,  dMedT,  dMedY, &
    !$ACC         Mp,  dMpdT,  dMpdY, &
    !$ACC         Mn,  dMndT,  dMndY, &
    !$ACC         Mnu, dMnudT, dMnudY, &
    !$ACC         U,   dUdT,   dUdY, dUdD )
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

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC IF( do_gpu ) &
      !$ACC PRESENT( Me,   dMedT,   dMedY, &
      !$ACC          Mp,   dMpdT,   dMpdY, &
      !$ACC          Mn,   dMndT,   dMndY, &
      !$ACC          Mnu,  dMnudT,  dMnudY )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO iX = iX_B, iX_E
        Mnu   (iX) = ( Me   (iX) + Mp   (iX) ) - Mn   (iX)
        dMnudT(iX) = ( dMedT(iX) + dMpdT(iX) ) - dMndT(iX)
        dMnudY(iX) = ( dMedY(iX) + dMpdY(iX) ) - dMndY(iX)
      END DO

    ELSE IF ( iSpecies == iNuE_Bar ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC IF( do_gpu ) &
      !$ACC PRESENT( Me,   dMedT,   dMedY, &
      !$ACC          Mp,   dMpdT,   dMpdY, &
      !$ACC          Mn,   dMndT,   dMndY, &
      !$ACC          Mnu,  dMnudT,  dMnudY )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO iX = iX_B, iX_E
        Mnu   (iX) = Mn   (iX) - ( Me   (iX) + Mp   (iX) )
        dMnudT(iX) = dMndT(iX) - ( dMedT(iX) + dMpdT(iX) )
        dMnudY(iX) = dMndY(iX) - ( dMedY(iX) + dMpdY(iX) )
      END DO

    ELSE

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
      !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG &
      !$ACC IF( do_gpu ) &
      !$ACC PRESENT( Mnu,  dMnudT,  dMnudY )
#elif defined(THORNADO_OMP)
      !$OMP PARALLEL DO SIMD
#endif
      DO iX = iX_B, iX_E
        Mnu   (iX) = Zero
        dMnudT(iX) = Zero
        dMnudY(iX) = Zero
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( kT ) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT ) &
    !$ACC PRESENT( Mnu,  dMnudT,  dMnudY, f0, df0dY, df0dU, E, T, dUdT, dUdY )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP PRIVATE( kT ) &
#endif
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: Me,  dMedT,  dMedY, &
    !$OMP               Mp,  dMpdT,  dMpdY, &
    !$OMP               Mn,  dMndT,  dMndY, &
    !$OMP               Mnu, dMnudT, dMnudY, &
    !$OMP               U,   dUdT,   dUdY, dUdD )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( Me,  dMedT,  dMedY, &
    !$ACC         Mp,  dMpdT,  dMpdY, &
    !$ACC         Mn,  dMndT,  dMndY, &
    !$ACC         Mnu, dMnudT, dMnudY, &
    !$ACC         U,   dUdT,   dUdY, dUdD )
#endif

  END SUBROUTINE ComputeEquilibriumDistributionAndDerivatives_Points


  SUBROUTINE ComputeNeutrinoOpacities_EC &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y, iSpecies )

    ! --- Electron Capture Opacities ---

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      T(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      Y(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):)
    INTEGER,  INTENT(in) :: &
      iSpecies

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, nZ(4)
    INTEGER  :: iNodeX, iNodeE, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iOS_X, iOS_E
    REAL(DP) :: D_K(1), T_K(1), Y_K(1), E_K(1), opEC_K(1,1)

    nZ = iZ_E0 - iZ_B0 + 1

#ifdef MICROPHYSICS_WEAKLIB

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
  
      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        D_K = D(iNodeX,iZ2,iZ3,iZ4)      
        T_K = T(iNodeX,iZ2,iZ3,iZ4)      
        Y_K = Y(iNodeX,iZ2,iZ3,iZ4)      
         
        ! --- Offset (Position) ---

        iOS_X = ( (iZ4-1)*nZ(3)*nZ(2) + (iZ3-1)*nZ(2) + (iZ2-1) ) * nDOFX

        DO iZ1 = iZ_B0(1), iZ_E0(1)
        DO iNodeE = 1, nDOFE

          E_K = NodeCoordinate( MeshE, iZ1, iNodeE )

          ! --- Offset (Energy) ---

          iOS_E = (iZ1-1) * nDOFE

          CALL LogInterpolateSingleVariable &
                 ( LOG10( E_K / UnitE ), LOG10( D_K / UnitD ), &
                   LOG10( T_K / UnitT ),      ( Y_K / UnitY ), &
                   LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                   OS_EmAb(iSpecies), EmAb_T(:,:,:,:,iSpecies), opEC_K )

          opEC(iOS_E+iNodeE,iSpecies,iOS_X+iNodeX) &
            = opEC_K(1,1) * UnitEC

        END DO ! iNodeE
        END DO ! iZ1
 
      END DO ! iNodeX

    END DO
    END DO
    END DO

#else

    opEC = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC


  SUBROUTINE ComputeNeutrinoOpacities_EC_Point &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, opEC_Point )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Electron Capture Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: opEC_Point(iE_B:iE_E)
    INTEGER,  INTENT(in)  :: iSpecies

    INTEGER :: iE

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL) || defined(THORNADO_OACC)

    DO iE = iE_B, iE_E

      CALL ComputeNeutrinoOpacity_Point &
             ( E(iE), D, T, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               opEC_Point(iE), EmAb_T(:,:,:,:,iSpecies), OS_EmAb(iSpecies), UnitEC )

    END DO

#else

    CALL LogInterpolateSingleVariable_1D3D_Custom_Point &
           ( LOG10( E / UnitE ), LOG10( D / UnitD ), &
             LOG10( T / UnitT ),      ( Y / UnitY ), &
             LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_EmAb(iSpecies), EmAb_T(:,:,:,:,iSpecies), opEC_Point )

    opEC_Point = opEC_Point * UnitEC

#endif

#else

    opEC_Point = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC_Point


  SUBROUTINE ComputeNeutrinoOpacities_EC_Points &
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

    INTEGER  :: iX, iE
    REAL(DP) :: LogE_P(iE_B:iE_E), LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opEC_Points )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   Y missing'
      IF ( .not. QueryOnGPU( opEC_Points ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   opEC_Points missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, LogD_P, T, LogT_P, Y, Y_P )
#endif
    DO iX = iX_B, iX_E
      LogD_P(iX) = LOG10( D(iX) / UnitD )
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      Y_P(iX) = Y(iX) / UnitY
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, LogE_P )
#endif
    DO iE = iE_B, iE_E
      LogE_P(iE) = LOG10( E(iE) / UnitE )
    END DO

    CALL LogInterpolateSingleVariable_1D3D_Custom &
           ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_EmAb(iSpecies), EmAb_T(:,:,:,:,iSpecies), opEC_Points, &
             GPU_Option = do_gpu )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( opEC_Points )
#endif
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E
        opEC_Points(iE,iX) = opEC_Points(iE,iX) * UnitEC
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#else

    opEC_Points = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC_Points

  SUBROUTINE ComputeNeutrinoOpacities_EC_Vector &
    ( iP_B, iP_E, E, D, T, Y, iSpecies, opEC_Points )

    ! --- Electron Capture Opacities (Multiple D,T,Y) ---
    ! --- Modified by Sherwood Richers to take in particle data ---

    INTEGER,  INTENT(in)  :: iP_B, iP_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iSpecies
    REAL(DP), INTENT(out) :: opEC_Points(:)

    INTEGER  :: iP
    REAL(DP) :: LogE_P(iP_B:iP_E), LogD_P(iP_B:iP_E), LogT_P(iP_B:iP_E), Y_P(iP_B:iP_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opEC_Points )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   Y missing'
      IF ( .not. QueryOnGPU( opEC_Points ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_EC_Points]   opEC_Points missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, LogD_P, T, LogT_P, Y, Y_P, E, LogE_P )
#endif
    DO iP = iP_B, iP_E
      LogD_P(iP) = LOG10( D(iP) / UnitD )
      LogT_P(iP) = LOG10( T(iP) / UnitT )
      Y_P(iP) = Y(iP) / UnitY
      LogE_P(iP) = LOG10( E(iP) / UnitE )
    END DO

    CALL LogInterpolateSingleVariable_4D_Custom &
           ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_EmAb(iSpecies), EmAb_T(:,:,:,:,iSpecies), opEC_Points)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( opEC_Points )
#endif
    DO iP = iP_B, iP_E
       opEC_Points(iP) = opEC_Points(iP) * UnitEC
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#else

    opEC_Points = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC_Vector


  SUBROUTINE ComputeNeutrinoOpacities_ES &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y, iSpecies )

    ! --- Elastic Scattering Opacities ---

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      T(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      Y(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):)
    INTEGER,  INTENT(in) :: &
      iSpecies

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, nZ(4)
    INTEGER  :: iNodeX, iNodeE, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iOS_X, iOS_E
    REAL(DP) :: D_K(1), T_K(1), Y_K(1), E_K(1), opES_K(1,1)

    nZ = iZ_E0 - iZ_B0 + 1

#ifdef MICROPHYSICS_WEAKLIB

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        D_K = D(iNodeX,iZ2,iZ3,iZ4)
        T_K = T(iNodeX,iZ2,iZ3,iZ4)
        Y_K = Y(iNodeX,iZ2,iZ3,iZ4)

        ! --- Offset (Position) ---

        iOS_X = ( (iZ4-1)*nZ(3)*nZ(2) + (iZ3-1)*nZ(2) + (iZ2-1) ) * nDOFX
 
        DO iZ1 = iZ_B0(1), iZ_E0(1)
        DO iNodeE = 1, nDOFE

          E_K = NodeCoordinate( MeshE, iZ1, iNodeE )

          ! --- Offset (Energy) ---

          iOS_E = (iZ1-1) * nDOFE

          CALL LogInterpolateSingleVariable &
                 ( LOG10( E_K / UnitE ), LOG10( D_K / UnitD ), &
                   LOG10( T_K / UnitT ),      ( Y_K / UnitY ), &
                   LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                   OS_Iso(iSpecies,1), &
                   Iso_T(:,:,:,:,1,iSpecies), &
                   opES_K )

          opES(iOS_E+iNodeE,iSpecies,iOS_X+iNodeX) &
            = opES_K(1,1) * UnitES

        END DO ! iNodeE
        END DO ! iZ1

      END DO ! iNodeX

    END DO
    END DO
    END DO

#else

    opES = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES


  SUBROUTINE ComputeNeutrinoOpacities_ES_Point &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, iMoment, opES_Point )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Elastic Scattering Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES_Point(iE_B:iE_E)

    INTEGER :: iE

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL) || defined(THORNADO_OACC)

    DO iE = iE_B, iE_E

      CALL ComputeNeutrinoOpacity_Point &
             ( E(iE), D, T, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               opES_Point(iE), Iso_T(:,:,:,:,iMoment,iSpecies), OS_Iso(iSpecies,iMoment), UnitES )

    END DO

#else

    CALL LogInterpolateSingleVariable_1D3D_Custom_Point &
           ( LOG10( E / UnitE ), LOG10( D / UnitD ), &
             LOG10( T / UnitT ),      ( Y / UnitY ), &
             LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_Iso(iSpecies,iMoment), Iso_T(:,:,:,:,iMoment,iSpecies), opES_Point )

    opES_Point = opES_Point * UnitES

#endif

#else

    opES_Point = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES_Point


  SUBROUTINE ComputeNeutrinoOpacities_ES_Points &
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

    INTEGER  :: iX, iE
    REAL(DP) :: LogE_P(iE_B:iE_E), LogD_P(iX_B:iX_E), LogT_P(iX_B:iX_E), Y_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opES_Points )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   Y missing'
      IF ( .not. QueryOnGPU( opES_Points ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   opES_Points missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, LogD_P, T, LogT_P, Y, Y_P )
#endif
    DO iX = iX_B, iX_E
      LogD_P(iX) = LOG10( D(iX) / UnitD )
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      Y_P(iX) = Y(iX) / UnitY
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, LogE_P )
#endif
    DO iE = iE_B, iE_E
      LogE_P(iE) = LOG10( E(iE) / UnitE )
    END DO

    CALL LogInterpolateSingleVariable_1D3D_Custom &
           ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_Iso(iSpecies,iMoment), Iso_T(:,:,:,:,iMoment,iSpecies), opES_Points, &
             GPU_Option = do_gpu )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( opES_Points )
#endif
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E
        opES_Points(iE,iX) = opES_Points(iE,iX) * UnitES
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#else

    opES_Points = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES_Points

  SUBROUTINE ComputeNeutrinoOpacities_ES_Vector &
    ( iP_B, ip_E, E, D, T, Y, iSpecies, iMoment, opES_Points )

    ! --- Elastic Scattering Opacities (Multiple D,T,Y) ---
    ! --- Modified by Sherwood Richers to take in particle data ---

    INTEGER,  INTENT(in)  :: ip_B, ip_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES_Points(:)

    INTEGER  :: iX, iE
    REAL(DP) :: LogE_P(ip_B:ip_E), LogD_P(ip_B:ip_E), LogT_P(ip_B:ip_E), Y_P(ip_B:ip_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( opES_Points )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   Y missing'
      IF ( .not. QueryOnGPU( opES_Points ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_ES_Points]   opES_Points missing'
    END IF
#endif

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( D, LogD_P, T, LogT_P, Y, Y_P, E, LogE_P )
#endif
    DO iP = iP_B, iP_E
      LogD_P(iP) = LOG10( D(iP) / UnitD )
      LogT_P(iP) = LOG10( T(iP) / UnitT )
      Y_P(iP) = Y(iP) / UnitY
      LogE_P(iE) = LOG10( E(iP) / UnitE )
    END DO

    CALL LogInterpolateSingleVariable_4D_Custom &
           ( LogE_P, LogD_P, LogT_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             OS_Iso(iSpecies,iMoment), Iso_T(:,:,:,:,iMoment,iSpecies), opES_Points)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( opES_Points )
#endif
    DO iP = iP_B, iP_E
       opES_Points(iP) = opES_Points(iP) * UnitES
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogD_P, LogT_P, Y_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogD_P, LogT_P, Y_P )
#endif

#else

    opES_Points = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES_Vector


  SUBROUTINE ComputeNeutrinoOpacities_NES_Point_1 &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Neutrino-Electron Scattering Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In (:,:)
    REAL(DP), INTENT(out) :: Phi_Out(:,:)

    INTEGER  :: iE1, iE2, iH1, iH2
    REAL(DP) :: H1(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: H2(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P, LogEta_P

#ifdef MICROPHYSICS_WEAKLIB

    iH1 = ( iMoment - 1 ) * 2 + 1
    iH2 = ( iMoment - 1 ) * 2 + 2

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    kT = BoltzmannConstant * T

    LogT_P = LOG10( T / UnitT )
    LogEta_P = LOG10( LogEta_P / kT / UnitEta )

    ! --- Interpolate HI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH1), NES_AT(:,:,:,:,iH1,1), H1 )

    ! --- Interpolate HII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH2), NES_AT(:,:,:,:,iH2,1), H2 )

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        DetBal = EXP( - ABS( E(iE2) - E(iE1) ) / kT )

        IF ( iE1 <= iE2 ) THEN
          Phi_Out(iE1,iE2) = ( C1(iSpecies) * H1(iE1,iE2) + C2(iSpecies) * H2(iE1,iE2) ) * UnitNES
          Phi_In (iE1,iE2) = Phi_Out(iE1,iE2) * DetBal
        ELSE
          Phi_In (iE1,iE2) = ( C1(iSpecies) * H1(iE2,iE1) + C2(iSpecies) * H2(iE2,iE1) ) * UnitNES
          Phi_Out(iE1,iE2) = Phi_In(iE1,iE2) * DetBal
        END IF
      END DO
    END DO

#else

    Phi_In  = Zero
    Phi_Out = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NES_Point_1


  SUBROUTINE ComputeNeutrinoOpacities_NES_Point_2 &
    ( iE_B, iE_E, E, D, T, Y, iS_1, iS_2, iMoment, &
      Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Neutrino-Electron Scattering Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iS_1, iS_2
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In_1(:,:), Phi_Out_1(:,:)
    REAL(DP), INTENT(out) :: Phi_In_2(:,:), Phi_Out_2(:,:)

    INTEGER  :: iE1, iE2, iH1, iH2
    REAL(DP) :: H1(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: H2(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P, LogEta_P

#ifdef MICROPHYSICS_WEAKLIB

    iH1 = ( iMoment - 1 ) * 2 + 1
    iH2 = ( iMoment - 1 ) * 2 + 2

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    kT = BoltzmannConstant * T

    LogT_P = LOG10( T / UnitT )
    LogEta_P = LOG10( LogEta_P / kT / UnitEta )

    ! --- Interpolate HI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH1), NES_AT(:,:,:,:,iH1,1), H1 )

    ! --- Interpolate HII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH2), NES_AT(:,:,:,:,iH2,1), H2 )

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

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

#else

    Phi_In_1  = Zero
    Phi_Out_1 = Zero
    Phi_In_2  = Zero
    Phi_Out_2 = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NES_Point_2


  SUBROUTINE ComputeNeutrinoOpacities_NES_Points_1 &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out, WORK1, WORK2 )

    ! --- Neutrino-Electron Scattering Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In (:,:,:)
    REAL(DP), INTENT(out) :: Phi_Out(:,:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK1(:,:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK2(:,:,:)

    REAL(DP), POINTER :: H1(:,:,:), H2(:,:,:)
    INTEGER  :: iX, iE1, iE2, iH1, iH2, nE, nX
    INTEGER  :: i, j, k
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( Phi_In, Phi_Out )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) &
        '[ComputeNeutrinoOpacities_NES_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   Y missing'
      IF ( .not. QueryOnGPU( Phi_In ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   Phi_In missing'
      IF ( .not. QueryOnGPU( Phi_Out ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   Phi_Out missing'
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogT_P, LogEta_P, H1, H2 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogT_P, LogEta_P, H1, H2 )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT ) &
    !$ACC PRESENT( T, LogT_P, LogEta_P )
#endif
    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )

      kT = BoltzmannConstant * T(iX)
      LogEta_P(iX) = LOG10( LogEta_P(iX) / kT / UnitEta )
    END DO

    ! --- Interpolate HI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH1), NES_AT(:,:,:,:,iH1,1), H1, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

    ! --- Interpolate HII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH2), NES_AT(:,:,:,:,iH2,1), H2, &
             GPU_Option = do_gpu, ASYNC_Option = 1  )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( kT, DetBal ) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) ASYNC(1) &
    !$ACC PRIVATE( kT, DetBal ) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, T, Phi_In, Phi_Out, H1, H2 )
#endif
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E

          kT = BoltzmannConstant * T(iX)
          DetBal = EXP( - ABS( E(iE2) - E(iE1) ) / kT )

          IF ( iE1 <= iE2 ) THEN
            Phi_Out(iE1,iE2,iX) = ( C1(iSpecies) * H1(iE1,iE2,iX) + C2(iSpecies) * H2(iE1,iE2,iX) ) * UnitNES
            Phi_In (iE1,iE2,iX) = Phi_Out(iE1,iE2,iX) * DetBal
          ELSE
            Phi_In (iE1,iE2,iX) = ( C1(iSpecies) * H1(iE2,iE1,iX) + C2(iSpecies) * H2(iE2,iE1,iX) ) * UnitNES
            Phi_Out(iE1,iE2,iX) = Phi_In(iE1,iE2,iX) * DetBal
          END IF
        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogT_P, LogEta_P, H1, H2 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogT_P, LogEta_P, H1, H2 )

    !$ACC WAIT(1)
#endif

    IF ( .NOT. PRESENT( WORK1 ) ) DEALLOCATE( H1 )
    IF ( .NOT. PRESENT( WORK2 ) ) DEALLOCATE( H2 )

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( Phi_In, Phi_Out )
#endif
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          Phi_Out(iE1,iE2,iX) = Zero
          Phi_In (iE1,iE2,iX) = Zero
        END DO
      END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NES_Points_1


  SUBROUTINE ComputeNeutrinoOpacities_NES_Points_2 &
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
    INTEGER  :: i, j, k
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( Phi_In_1, Phi_Out_1 ) &
       .AND. QueryOnGPU( Phi_In_2, Phi_Out_2 )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) &
        '[ComputeNeutrinoOpacities_NES_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   Y missing'
      IF ( .not. QueryOnGPU( Phi_In_1 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   Phi_In_1 missing'
      IF ( .not. QueryOnGPU( Phi_Out_1 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   Phi_Out_1 missing'
      IF ( .not. QueryOnGPU( Phi_In_2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   Phi_In_2 missing'
      IF ( .not. QueryOnGPU( Phi_Out_2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points]   Phi_Out_2 missing'
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogT_P, LogEta_P, H1, H2 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogT_P, LogEta_P, H1, H2 )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT ) &
    !$ACC PRESENT( T, LogT_P, LogEta_P )
#endif
    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )

      kT = BoltzmannConstant * T(iX)
      LogEta_P(iX) = LOG10( LogEta_P(iX) / kT / UnitEta )
    END DO

    ! --- Interpolate HI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH1), NES_AT(:,:,:,:,iH1,1), H1, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

    ! --- Interpolate HII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_NES(1,iH2), NES_AT(:,:,:,:,iH2,1), H2, &
             GPU_Option = do_gpu, ASYNC_Option = 1  )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( kT, DetBal ) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) ASYNC(1) &
    !$ACC PRIVATE( kT, DetBal ) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, T, Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2, H1, H2 )
#endif
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogT_P, LogEta_P, H1, H2 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogT_P, LogEta_P, H1, H2 )

    !$ACC WAIT(1)
#endif

    IF ( .NOT. PRESENT( WORK1 ) ) DEALLOCATE( H1 )
    IF ( .NOT. PRESENT( WORK2 ) ) DEALLOCATE( H2 )

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2 )
#endif
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

  END SUBROUTINE ComputeNeutrinoOpacities_NES_Points_2


  SUBROUTINE ComputeNeutrinoOpacitiesAndDerivatives_NES_Point &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out, &
      dPhi_In_dY, dPhi_In_dE, dPhi_Out_dY, dPhi_Out_dE )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Neutrino-Electron Scattering Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In     (iE_B:iE_E,iE_B:iE_E)
    REAL(DP), INTENT(out) :: Phi_Out    (iE_B:iE_E,iE_B:iE_E)
    REAL(DP), INTENT(out) :: dPhi_In_dY (iE_B:iE_E,iE_B:iE_E)
    REAL(DP), INTENT(out) :: dPhi_In_dE (iE_B:iE_E,iE_B:iE_E)
    REAL(DP), INTENT(out) :: dPhi_Out_dY(iE_B:iE_E,iE_B:iE_E)
    REAL(DP), INTENT(out) :: dPhi_Out_dE(iE_B:iE_E,iE_B:iE_E)

    INTEGER  :: iE1, iE2, iH1, iH2
    REAL(DP) :: kT, LogT, LogEta, C_Eta, C_T
    REAL(DP) :: M, dMdD, dMdT, dMdY
    REAL(DP) :: U, dUdD, dUdT, dUdY
    REAL(DP) :: H1    (iE_B:iE_E,iE_B:iE_E), H2      (iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: dH1dT (iE_B:iE_E,iE_B:iE_E), dH1dEta (iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: dH2dT (iE_B:iE_E,iE_B:iE_E), dH2dEta (iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: dPhidT(iE_B:iE_E,iE_B:iE_E), dPhidEta(iE_B:iE_E,iE_B:iE_E)

#ifdef MICROPHYSICS_WEAKLIB

    iH1 = ( iMoment - 1 ) * 2 + 1
    iH2 = ( iMoment - 1 ) * 2 + 2

    ! --- Compute Electron Chemical Potential and Derivatives ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, M, dMdD, dMdT, dMdY )

    kT = BoltzmannConstant * T

    LogT   = LOG10( T / UnitT )
    LogEta = LOG10( M / kT )

    ! --- Interpolate HI  ---

    CALL LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P &
           ( LogT, LogEta, LogTs_T, LogEtas_T, &
             OS_NES(1,iH1), NES_AT(:,:,:,:,iH1,1), H1, dH1dT, dH1dEta )

    ! --- Interpolate HII ---

    CALL LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P &
           ( LogT, LogEta, LogTs_T, LogEtas_T, &
             OS_NES(1,iH2), NES_AT(:,:,:,:,iH2,1), H2, dH2dT, dH2dEta )

    ! --- Phi_Out ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE2
      Phi_Out(iE1,iE2) = ( C1(iSpecies) * H1(iE1,iE2) + C2(iSpecies) * H2(iE1,iE2) ) * UnitNES
    END DO
    END DO

    DO iE2 = iE_B,  iE_E
    DO iE1 = iE2+1, iE_E
      Phi_Out(iE1,iE2) = Phi_Out(iE2,iE1) * EXP( ( E(iE2) - E(iE1) ) / kT )
    END DO
    END DO

    ! --- Phi_In ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      Phi_In(iE1,iE2) = Phi_Out(iE2,iE1)
    END DO
    END DO

    ! --- Derivative of Phi_Out wrt. Temperature ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE2
      dPhidT(iE1,iE2) &
        = ( C1(iSpecies) * dH1dT(iE1,iE2) + C2(iSpecies) * dH2dT(iE1,iE2) ) * UnitNES / UnitT
    END DO
    END DO

    DO iE2 = iE_B,  iE_E
    DO iE1 = iE2+1, iE_E
      dPhidT(iE1,iE2) &
        = ( dPhidT(iE2,iE1) &
            - ( Phi_Out(iE2,iE1) / T ) * ( E(iE2) - E(iE1) ) / kT ) &
          * EXP( ( E(iE2) - E(iE1) ) / kT )
    END DO
    END DO

    ! --- Derivative of Phi_Out wrt. Degeneracy Parameter ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE2
      dPhidEta(iE1,iE2) &
        = ( C1(iSpecies) * dH1dEta(iE1,iE2) + C2(iSpecies) * dH2dEta(iE1,iE2) ) * UnitNES
    END DO
    END DO

    DO iE2 = iE_B,  iE_E
    DO iE1 = iE2+1, iE_E
      dPhidEta(iE1,iE2) &
        = dPhidEta(iE2,iE1) * EXP( ( E(iE2) - E(iE1) ) / kT )
    END DO
    END DO

    ! --- Compute Specific Internal Energy and Derivatives ---

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( D, T, Y, U, dUdD, dUdT, dUdY )

    ! --- Derivative of Phi_Out wrt. Electron Fraction ---

    C_Eta = ( dMdY - dUdY * ( dMdT - M / T ) / dUdT ) / kT
    C_T   = - dUdY / dUdT

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      dPhi_Out_dY(iE1,iE2) &
        = C_Eta * dPhidEta(iE1,iE2) + C_T * dPhidT(iE1,iE2)
    END DO
    END DO

    ! --- Derivative of Phi_In wrt. Electron Fraction ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      dPhi_In_dY(iE1,iE2) = dPhi_Out_dY(iE2,iE1)
    END DO
    END DO

    ! --- Derivative of Phi_Out wrt. Specific Internal Energy ---

    C_Eta = ( dMdT - M / T ) / ( dUdT * kT )
    C_T   = One / dUdT

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      dPhi_Out_dE(iE1,iE2) &
        = C_Eta * dPhidEta(iE1,iE2) + C_T * dPhidT(iE1,iE2)
    END DO
    END DO

    ! --- Derivative of Phi_In wrt. Specific Internal Energy ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      dPhi_In_dE(iE1,iE2) = dPhi_Out_dE(iE2,iE1)
    END DO
    END DO

#else

    Phi_In      = Zero
    Phi_Out     = Zero
    dPhi_In_dY  = Zero
    dPhi_In_dE  = Zero
    dPhi_Out_dY = Zero
    dPhi_Out_dE = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacitiesAndDerivatives_NES_Point


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

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( SUM1, SUM2 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( SUM1, SUM2 ) &
      !$ACC PRESENT( Phi_In, Phi_Out, Eta, Chi, W2, J )
#endif
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


  SUBROUTINE ComputeNeutrinoOpacities_Pair_Point_1 &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Pair Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In (:,:)
    REAL(DP), INTENT(out) :: Phi_Out(:,:)

    INTEGER  :: iE1, iE2, iJ1, iJ2
    REAL(DP) :: J1(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: J2(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P, LogEta_P

#ifdef MICROPHYSICS_WEAKLIB

    iJ1 = ( iMoment - 1 ) * 2 + 1
    iJ2 = ( iMoment - 1 ) * 2 + 2

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    kT = BoltzmannConstant * T

    LogT_P = LOG10( T / UnitT )
    LogEta_P = LOG10( LogEta_P / kT / UnitEta )

    ! --- Interpolate JI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ1), Pair_AT(:,:,:,:,iJ1,1), J1 )

    ! --- Interpolate JII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ2), Pair_AT(:,:,:,:,iJ2,1), J2 )

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        DetBal = EXP( - ABS( E(iE1) + E(iE2) ) / kT )

        IF ( iE1 <= iE2 ) THEN
          Phi_Out(iE1,iE2) = ( C1(iSpecies) * J1(iE1,iE2) + C2(iSpecies) * J2(iE1,iE2) ) * UnitPair
        ELSE
          Phi_Out(iE1,iE2) = ( C1(iSpecies) * J2(iE2,iE1) + C2(iSpecies) * J1(iE2,iE1) ) * UnitPair
        END IF
        Phi_In(iE1,iE2) = Phi_Out(iE1,iE2) * DetBal

      END DO
    END DO

#else

    Phi_In  = Zero
    Phi_Out = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Pair_Point_1


  SUBROUTINE ComputeNeutrinoOpacities_Pair_Point_2 &
    ( iE_B, iE_E, E, D, T, Y, iS_1, iS_2, iMoment, &
      Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2 )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Pair Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iS_1, iS_2
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In_1(:,:), Phi_Out_1(:,:)
    REAL(DP), INTENT(out) :: Phi_In_2(:,:), Phi_Out_2(:,:)

    INTEGER  :: iE1, iE2, iJ1, iJ2
    REAL(DP) :: J1(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: J2(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P, LogEta_P

#ifdef MICROPHYSICS_WEAKLIB

    iJ1 = ( iMoment - 1 ) * 2 + 1
    iJ2 = ( iMoment - 1 ) * 2 + 2

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    kT = BoltzmannConstant * T

    LogT_P = LOG10( T / UnitT )
    LogEta_P = LOG10( LogEta_P / kT / UnitEta )

    ! --- Interpolate JI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ1), Pair_AT(:,:,:,:,iJ1,1), J1 )

    ! --- Interpolate JII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ2), Pair_AT(:,:,:,:,iJ2,1), J2 )

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

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

#else

    Phi_In_1  = Zero
    Phi_Out_1 = Zero
    Phi_In_2  = Zero
    Phi_Out_2 = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Pair_Point_2


  SUBROUTINE ComputeNeutrinoOpacities_Pair_Points_1 &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out, WORK1, WORK2 )

    ! --- Pair Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D(:)
    REAL(DP), INTENT(in)  :: T(:)
    REAL(DP), INTENT(in)  :: Y(:)
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In (:,:,:)
    REAL(DP), INTENT(out) :: Phi_Out(:,:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK1(:,:,:)
    REAL(DP), INTENT(out), TARGET, OPTIONAL :: WORK2(:,:,:)

    REAL(DP), POINTER :: J1(:,:,:), J2(:,:,:)
    INTEGER  :: iX, iE1, iE2, iJ1, iJ2, nE, nX
    INTEGER  :: i, j, k
    REAL(DP) :: kT, DetBal
    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( Phi_In, Phi_Out )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   Y missing'
      IF ( .not. QueryOnGPU( Phi_In ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   Phi_In missing'
      IF ( .not. QueryOnGPU( Phi_Out ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   Phi_Out missing'
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogT_P, LogEta_P, J1, J2 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogT_P, LogEta_P, J1, J2 )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( T, LogT_P, LogEta_P )
#endif
    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

    ! --- Interpolate JI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ1), Pair_AT(:,:,:,:,iJ1,1), J1, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

    ! --- Interpolate JII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ2), Pair_AT(:,:,:,:,iJ2,1), J2, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( kT, DetBal )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT, DetBal ) &
    !$ACC PRESENT( E, T, Phi_In, Phi_Out, J1, J2 )
#endif
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E

          kT = BoltzmannConstant * T(iX)
          DetBal = EXP( - ABS( E(iE1) + E(iE2) ) / kT )

          IF ( iE1 <= iE2 ) THEN
            Phi_Out(iE1,iE2,iX) = ( C1(iSpecies) * J1(iE1,iE2,iX) + C2(iSpecies) * J2(iE1,iE2,iX) ) * UnitPair
          ELSE
            Phi_Out(iE1,iE2,iX) = ( C1(iSpecies) * J2(iE2,iE1,iX) + C2(iSpecies) * J1(iE2,iE1,iX) ) * UnitPair
          END IF
          Phi_In(iE1,iE2,iX) = Phi_Out(iE1,iE2,iX) * DetBal

        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogT_P, LogEta_P, J1, J2 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogT_P, LogEta_P, J1, J2 )

    !$ACC WAIT(1)
#endif

    IF ( .NOT. PRESENT( WORK1 ) ) DEALLOCATE( J1 )
    IF ( .NOT. PRESENT( WORK2 ) ) DEALLOCATE( J2 )

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( Phi_In, Phi_Out )
#endif
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          Phi_Out(iE1,iE2,iX) = Zero
          Phi_In (iE1,iE2,iX) = Zero
        END DO
      END DO
    END DO

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Pair_Points_1


  SUBROUTINE ComputeNeutrinoOpacities_Pair_Points_2 &
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
    REAL(DP) :: LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( Phi_In_1, Phi_Out_1 ) &
       .AND. QueryOnGPU( Phi_In_2, Phi_Out_2 )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points] Data not present on device'
      IF ( .not. QueryOnGPU( E ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   E missing'
      IF ( .not. QueryOnGPU( D ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   D missing'
      IF ( .not. QueryOnGPU( T ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   T missing'
      IF ( .not. QueryOnGPU( Y ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   Y missing'
      IF ( .not. QueryOnGPU( Phi_In_1 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   Phi_In_1 missing'
      IF ( .not. QueryOnGPU( Phi_Out_1 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   Phi_Out_1 missing'
      IF ( .not. QueryOnGPU( Phi_In_2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   Phi_In_2 missing'
      IF ( .not. QueryOnGPU( Phi_Out_2 ) ) &
        WRITE(*,*) '[ComputeNeutrinoOpacities_Pair_Points]   Phi_Out_2 missing'
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogT_P, LogEta_P, J1, J2 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogT_P, LogEta_P, J1, J2 )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( T, LogT_P, LogEta_P )
#endif
    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

    ! --- Interpolate JI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ1), Pair_AT(:,:,:,:,iJ1,1), J1, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

    ! --- Interpolate JII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
           ( LogT_P, LogEta_P, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ2), Pair_AT(:,:,:,:,iJ2,1), J2, &
             GPU_Option = do_gpu, ASYNC_Option = 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( kT, DetBal )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) ASYNC(1) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT, DetBal ) &
    !$ACC PRESENT( E, T, Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2, J1, J2 )
#endif
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogT_P, LogEta_P, J1, J2 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogT_P, LogEta_P, J1, J2 )

    !$ACC WAIT(1)
#endif

    IF ( .NOT. PRESENT( WORK1 ) ) DEALLOCATE( J1 )
    IF ( .NOT. PRESENT( WORK2 ) ) DEALLOCATE( J2 )

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( Phi_In_1, Phi_Out_1, Phi_In_2, Phi_Out_2 )
#endif
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

  END SUBROUTINE ComputeNeutrinoOpacities_Pair_Points_2


  SUBROUTINE ComputeNeutrinoOpacitiesAndDerivatives_Pair_Point &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out, &
      dPhi_In_dY, dPhi_In_dE, dPhi_Out_dY, dPhi_Out_dE )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Pair Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(:)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In     (:,:)
    REAL(DP), INTENT(out) :: Phi_Out    (:,:)
    REAL(DP), INTENT(out) :: dPhi_In_dY (:,:)
    REAL(DP), INTENT(out) :: dPhi_In_dE (:,:)
    REAL(DP), INTENT(out) :: dPhi_Out_dY(:,:)
    REAL(DP), INTENT(out) :: dPhi_Out_dE(:,:)

    INTEGER  :: iE1, iE2, iJ1, iJ2
    REAL(DP) :: kT, LogT, LogEta, C_Eta, C_T
    REAL(DP) :: M, dMdD, dMdT, dMdY
    REAL(DP) :: U, dUdD, dUdT, dUdY
    REAL(DP) :: J1    (iE_B:iE_E,iE_B:iE_E), J2      (iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: dJ1dT (iE_B:iE_E,iE_B:iE_E), dJ1dEta (iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: dJ2dT (iE_B:iE_E,iE_B:iE_E), dJ2dEta (iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: dPhidT(iE_B:iE_E,iE_B:iE_E), dPhidEta(iE_B:iE_E,iE_B:iE_E)

#ifdef MICROPHYSICS_WEAKLIB

    iJ1 = ( iMoment - 1 ) * 2 + 1
    iJ2 = ( iMoment - 1 ) * 2 + 2

    ! --- Compute Electron Chemical Potential and Derivatives ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, M, dMdD, dMdT, dMdY )

    kT = BoltzmannConstant * T

    LogT   = LOG10( T / UnitT )
    LogEta = LOG10( M / kT )

    ! --- Interpolate JI  ---

    CALL LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P &
           ( LogT, LogEta, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ1), Pair_AT(:,:,:,:,iJ1,1), J1, dJ1dT, dJ1dEta )

    ! --- Interpolate JII ---

    CALL LogInterpolateDifferentiateSingleVariable_2D2D_Custom_Aligned_P &
           ( LogT, LogEta, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ2), Pair_AT(:,:,:,:,iJ2,1), J2, dJ2dT, dJ2dEta )

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E

      IF( iE1 <= iE2 )THEN
        Phi_Out(iE1,iE2) &
          = ( C1(iSpecies) * J1(iE1,iE2) + C2(iSpecies) * J2(iE1,iE2) ) * UnitPair
      ELSE
        Phi_Out(iE1,iE2) &
          = ( C1(iSpecies) * J2(iE2,iE1) + C2(iSpecies) * J1(iE2,iE1) ) * UnitPair
      END IF

      Phi_In(iE1,iE2) &
        = Phi_Out(iE1,iE2) * EXP( - ( E(iE1) + E(iE2) ) / kT )

    END DO
    END DO

    ! --- Derivative of Phi_Out wrt. Temperature ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E

      IF( iE1 <= iE2 )THEN
        dPhidT(iE1,iE2) &
          = ( C1(iSpecies) * dJ1dT(iE1,iE2) + C2(iSpecies) * dJ2dT(iE1,iE2) ) * UnitPair / UnitT
      ELSE
        dPhidT(iE1,iE2) &
          = ( C1(iSpecies) * dJ2dT(iE2,iE1) + C2(iSpecies) * dJ1dT(iE2,iE1) ) * UnitPair / UnitT
      END IF

    END DO
    END DO

    ! --- Derivative of Phi_Out wrt. Degeneracy Parameter ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E

      IF( iE1 <= iE2 )THEN
        dPhidEta(iE1,iE2) &
          = ( C1(iSpecies) * dJ1dEta(iE1,iE2) + C2(iSpecies) * dJ2dEta(iE1,iE2) ) * UnitPair
      ELSE
        dPhidEta(iE1,iE2) &
          = ( C1(iSpecies) * dJ2dEta(iE2,iE1) + C2(iSpecies) * dJ1dEta(iE2,iE1) ) * UnitPair
      END IF

    END DO
    END DO

    ! --- Compute Specific Internal Energy and Derivatives ---

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( D, T, Y, U, dUdD, dUdT, dUdY )

    ! --- Derivative of Phi_Out wrt. Electron Fraction ---

    C_Eta = ( dMdY - dUdY * ( dMdT - M / T ) / dUdT ) / kT
    C_T   = - dUdY / dUdT

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      dPhi_Out_dY(iE1,iE2) &
        = C_Eta * dPhidEta(iE1,iE2) + C_T * dPhidT(iE1,iE2)
    END DO
    END DO

    ! --- Derivative of Phi_In wrt. Electron Fraction ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      dPhi_In_dY(iE1,iE2) &
        = ( C_T * ( dPhidT(iE1,iE2) &
                      + ( Phi_Out(iE1,iE2) / T ) * ( E(iE1) + E(iE2) ) / kT ) &
            + C_Eta * dPhidEta(iE1,iE2) ) * EXP( - ( E(iE1) + E(iE2) ) / kT )
    END DO
    END DO

    ! --- Derivative of Phi_Out wrt. Specific Internal Energy ---

    C_Eta = ( dMdT - M / T ) / ( dUdT * kT )
    C_T   = One / dUdT

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      dPhi_Out_dE(iE1,iE2) &
        = C_Eta * dPhidEta(iE1,iE2) + C_T * dPhidT(iE1,iE2)
    END DO
    END DO

    ! --- Derivative of Phi_In wrt. Specific Internal Energy ---

    DO iE2 = iE_B, iE_E
    DO iE1 = iE_B, iE_E
      dPhi_In_dE(iE1,iE2) &
        = ( C_T * ( dPhidT(iE1,iE2) &
                      + ( Phi_Out(iE1,iE2) / T ) * ( E(iE1) + E(iE2) ) / kT ) &
            + C_Eta * dPhidEta(iE1,iE2) ) * EXP( - ( E(iE1) + E(iE2) ) / kT )
    END DO
    END DO

#else

    Phi_In      = Zero
    Phi_Out     = Zero
    dPhi_In_dY  = Zero
    dPhi_In_dE  = Zero
    dPhi_Out_dY = Zero
    dPhi_Out_dE = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacitiesAndDerivatives_Pair_Point


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

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( SUM1, SUM2 )
#elif defined(THORNADO_OACC)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( SUM1, SUM2 ) &
      !$ACC PRESENT( Phi_In, Phi_Out, Eta, Chi, W2, J )
#endif
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


  SUBROUTINE ComputeNeutrinoOpacityE_Point &
    ( E, D, T, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, Op, Op_T, OS_Op, Units_Op )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), DIMENSION(:),       INTENT(in)  :: E
    REAL(DP),                     INTENT(in)  :: D, T, Y
    REAL(DP), DIMENSION(:),       INTENT(in)  :: LogEs_T, LogDs_T, LogTs_T, Ys_T
    REAL(DP), DIMENSION(:),       INTENT(out) :: Op
    REAL(DP), DIMENSION(:,:,:,:), INTENT(in)  :: Op_T
    REAL(DP),                     INTENT(in)  :: OS_Op, Units_Op

    REAL(DP), DIMENSION(SIZE(E)) :: E_P, Op_P
    REAL(DP) :: D_P, T_P, Y_P

#ifdef MICROPHYSICS_WEAKLIB

    E_P(:) = LOG10( E(:) / UnitE )
    D_P    = LOG10( D    / UnitD )
    T_P    = LOG10( T    / UnitT )
    Y_P    = Y / UnitY

    CALL LogInterpolateSingleVariable &
           ( E_P, D_P, T_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_Op, Op_T, Op_P )

    Op(:) = Op_P(:) * Units_Op

#else

    Op(:) = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacityE_Point


  SUBROUTINE ComputeNeutrinoOpacity_E_D_T_Y_Point &
    ( E, D, T, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, Op, Op_T, OS_Op, Units_Op )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP),                     INTENT(in)  :: E, D, T, Y
    REAL(DP), DIMENSION(:),       INTENT(in)  :: LogEs_T, LogDs_T, LogTs_T, Ys_T
    REAL(DP),                     INTENT(out) :: Op
    REAL(DP), DIMENSION(:,:,:,:), INTENT(in)  :: Op_T
    REAL(DP),                     INTENT(in)  :: OS_Op, Units_Op

    REAL(DP) :: E_P, D_P, T_P, Y_P, Op_P

#ifdef MICROPHYSICS_WEAKLIB

    E_P = LOG10( E / UnitE )
    D_P = LOG10( D / UnitD )
    T_P = LOG10( T / UnitT )
    Y_P = Y / UnitY

    CALL LogInterpolateSingleVariable &
           ( E_P, D_P, T_P, Y_P, LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_Op, Op_T, Op_P )

    Op = Op_P * Units_Op

#else

    Op = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacity_E_D_T_Y_Point


  SUBROUTINE ComputeNeutrinoOpacity_E_E_T_Eta_Point &
    ( E1, E2, T, Eta, LogEs_T, LogTs_T, LogEtas_T, Op, Op_T, OS_Op, Units_Op )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP),                     INTENT(in)  :: E1, E2, T, Eta
    REAL(DP), DIMENSION(:),       INTENT(in)  :: LogEs_T, LogTs_T, LogEtas_T
    REAL(DP),                     INTENT(out) :: Op
    REAL(DP), DIMENSION(:,:,:,:), INTENT(in)  :: Op_T
    REAL(DP),                     INTENT(in)  :: OS_Op, Units_Op

    REAL(DP) :: E1_P, E2_P, T_P, Eta_P, Op_P

#ifdef MICROPHYSICS_WEAKLIB

    E1_P  = LOG10( E1 / UnitE )
    E2_P  = LOG10( E2 / UnitE )
    T_P   = LOG10( T / UnitT )
    Eta_P = LOG10( Eta / UnitEta )

    CALL LogInterpolateSingleVariable &
           ( E1_P, E2_P, T_P, Eta_P, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, OS_Op, Op_T, Op_P )

    Op = Op_P * Units_Op

#else

    Op = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacity_E_E_T_Eta_Point


  FUNCTION FermiDirac_Scalar( E, Mu, kT ) RESULT( FermiDirac )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E, Mu, kT
    REAL(DP) :: FermiDirac

    REAL(DP) :: Exponent

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FermiDirac &
      = One / ( EXP( Exponent ) + One )

    RETURN
  END FUNCTION FermiDirac_Scalar


  FUNCTION FermiDirac_Vector( E, Mu, kT ) RESULT( FermiDirac )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E(:), Mu, kT
    REAL(DP) :: FermiDirac(SIZE(E))

    REAL(DP) :: Exponent(SIZE(E))

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FermiDirac &
      = One / ( EXP( Exponent ) + One )

    RETURN
  END FUNCTION FermiDirac_Vector


  FUNCTION dFermiDiracdT_Scalar( E, Mu, kT, dMudT, T ) RESULT( dFermiDiracdT )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

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
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: E(:), Mu, kT, dMudY, T
    REAL(DP) :: dFermiDiracdY(SIZE(E))

    REAL(DP) :: Exponent(SIZE(E)), FD(SIZE(E))

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdY &
      = ( FD * EXP( Exponent ) ) * FD * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY_Vector


END MODULE NeutrinoOpacitiesComputationModule
