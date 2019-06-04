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
    nDOFE, nDOFX
  USE DeviceModule, ONLY: &
    QueryOnGpu
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE MeshModule, ONLY: &
    MeshE
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeElectronChemicalPotential_TABLE, &
    ComputeProtonChemicalPotential_TABLE, &
    ComputeNeutronChemicalPotential_TABLE
  USE OpacityModule_TABLE, ONLY: &
#ifdef MICROPHYSICS_WEAKLIB
    OS_EmAb, OS_Iso, OS_NES, OS_Pair, &
    EmAb_T, Iso_T, NES_T, Pair_T, &
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
    LogInterpolateSingleVariable_2D2D_Custom_Point

  ! ----------------------------------------------

#endif

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeNeutrinoOpacities
  PUBLIC :: ComputeNeutrinoOpacities_EC_Point
  PUBLIC :: ComputeNeutrinoOpacities_EC_Points
  PUBLIC :: ComputeNeutrinoOpacities_ES_Point
  PUBLIC :: ComputeNeutrinoOpacities_ES_Points
  PUBLIC :: ComputeNeutrinoOpacities_NES_Point
  PUBLIC :: ComputeNeutrinoOpacities_NES_Points
  PUBLIC :: ComputeNeutrinoOpacities_Pair_Point
  PUBLIC :: ComputeNeutrinoOpacities_Pair_Points
  PUBLIC :: ComputeEquilibriumDistributions_Point
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

  INTERFACE ComputeNeutrinoOpacity_Point
    MODULE PROCEDURE ComputeNeutrinoOpacity_E_D_T_Y_Point
    MODULE PROCEDURE ComputeNeutrinoOpacity_E_E_T_Eta_Point
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


  SUBROUTINE ComputeEquilibriumDistributions_Point &
    ( iE_B, iE_E, E, D, T, Y, f_EQ_Point, iSpecies )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Equilibrium Neutrino Distributions (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: f_EQ_Point(iE_B:iE_E)
    INTEGER,  INTENT(in)  :: iSpecies

    INTEGER  :: iE
    REAL(DP) :: Me, Mp, Mn, Mnu, kT, FD_Exp

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( D, T, Y, Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( D, T, Y, Mn )

    IF( iSpecies .EQ. iNuE )THEN

      ! --- Electron Neutrinos ---

      Mnu = ( Me + Mp ) - Mn

    ELSEIF( iSpecies .EQ. iNuE_Bar )THEN

      ! --- Electron Antineutrino ---

      Mnu = Mn - ( Me + Mp )

    ELSE

      Mnu = Zero

    END IF

    kT = BoltzmannConstant * T
    DO iE = iE_B, iE_E
      FD_Exp = ( E(iE) - Mnu ) / kT
      FD_Exp = MIN( MAX( FD_Exp, - Log1d100 ), + Log1d100 )
      f_EQ_Point(iE) = One / ( EXP( FD_Exp ) + One )
    END DO

  END SUBROUTINE ComputeEquilibriumDistributions_Point


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
    INTEGER  :: iOS_X, iOS_E
    REAL(DP) :: D_K(1), T_K(1), Y_K(1), E_K(1), opEC_K(1,1)

    nZ = iZ_E0 - iZ_B0 + 1

#ifdef MICROPHYSICS_WEAKLIB

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      ! --- Use Cell Averaged D, T, Y ---

      D_K = DOT_PRODUCT( WeightsX_q(:), D(:,iZ2,iZ3,iZ4) )
      T_K = DOT_PRODUCT( WeightsX_q(:), T(:,iZ2,iZ3,iZ4) )
      Y_K = DOT_PRODUCT( WeightsX_q(:), Y(:,iZ2,iZ3,iZ4) )

      ! --- Offset (Position) ---

      iOS_X = ( (iZ4-1)*nZ(3)*nZ(2) + (iZ3-1)*nZ(2) + (iZ2-1) ) * nDOFX

      DO iZ1 = iZ_B0(1), iZ_E0(1)

        ! --- Use Cell Center Energy ---

        E_K = MeshE % Center(iZ1)

        ! --- Offset (Energy) ---

        iOS_E = (iZ1-1) * nDOFE

        CALL LogInterpolateSingleVariable &
               ( LOG10( E_K / UnitE ), LOG10( D_K / UnitD ), &
                 LOG10( T_K / UnitT ),      ( Y_K / UnitY ), &
                 LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                 OS_EmAb(iSpecies), EmAb_T(:,:,:,:,iSpecies), opEC_K )

        opEC(iOS_E+1:iOS_E+nDOFE,iSpecies,iOS_X+1:iOS_X+nDOFX) &
          = opEC_K(1,1) * UnitEC

      END DO

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

#if defined(THORNADO_OMP_XL) || defined(THORNADO_OACC)

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
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iSpecies
    REAL(DP), INTENT(out) :: opEC_Points(iE_B:iE_E,iX_B:iX_E)

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
    INTEGER  :: iOS_X, iOS_E
    REAL(DP) :: D_K(1), T_K(1), Y_K(1), E_K(1), opES_K(1,1)

    nZ = iZ_E0 - iZ_B0 + 1

#ifdef MICROPHYSICS_WEAKLIB

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      ! --- Use Cell Averaged D, T, Y ---

      D_K = DOT_PRODUCT( WeightsX_q(:), D(:,iZ2,iZ3,iZ4) )
      T_K = DOT_PRODUCT( WeightsX_q(:), T(:,iZ2,iZ3,iZ4) )
      Y_K = DOT_PRODUCT( WeightsX_q(:), Y(:,iZ2,iZ3,iZ4) )

      ! --- Offset (Position) ---

      iOS_X = ( (iZ4-1)*nZ(3)*nZ(2) + (iZ3-1)*nZ(2) + (iZ2-1) ) * nDOFX

      DO iZ1 = iZ_B0(1), iZ_E0(1)

        ! --- Use Cell Center Energy ---

        E_K = MeshE % Center(iZ1)

        ! --- Offset (Energy) ---

        iOS_E = (iZ1-1) * nDOFE

        CALL LogInterpolateSingleVariable &
               ( LOG10( E_K / UnitE ), LOG10( D_K / UnitD ), &
                 LOG10( T_K / UnitT ),      ( Y_K / UnitY ), &
                 LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                 OS_Iso(iSpecies,1), &
                 Iso_T(:,:,:,:,1,iSpecies), &
                 opES_K )

        opES(iOS_E+1:iOS_E+nDOFE,iSpecies,iOS_X+1:iOS_X+nDOFX) &
          = opES_K(1,1) * UnitES

      END DO

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

#if defined(THORNADO_OMP_XL) || defined(THORNADO_OACC)

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
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES_Points(iE_B:iE_E,iX_B:iX_E)

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


  SUBROUTINE ComputeNeutrinoOpacities_NES_Point &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out )
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
    REAL(DP), INTENT(out) :: Phi_In (iE_B:iE_E,iE_B:iE_E)
    REAL(DP), INTENT(out) :: Phi_Out(iE_B:iE_E,iE_B:iE_E)

    INTEGER  :: iE1, iE2, iH1, iH2
    REAL(DP) :: C1, C2, kT, Phi_Local
    REAL(DP) :: LogE_P(iE_B:iE_E), LogT_P, LogEta_P

#ifdef MICROPHYSICS_WEAKLIB

    iH1 = ( iMoment - 1 ) * 2 + 1
    iH2 = ( iMoment - 1 ) * 2 + 2

    IF(     iSpecies == iNuE     )THEN

      C1 = ( cv + ca )**2
      C2 = ( cv - ca )**2

    ELSEIF( iSpecies == iNuE_Bar )THEN

      C1 = ( cv - ca )**2
      C2 = ( cv + ca )**2

    END IF

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    kT = BoltzmannConstant * T

    LogT_P = LOG10( T / UnitT )
    LogEta_P = LOG10( LogEta_P / kT / UnitEta )

    DO iE1 = iE_B, iE_E
      LogE_P(iE1) = LOG10( E(iE1) / UnitE )
    END DO

    ! --- Interpolate HI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Point &
           ( LogE_P, LogT_P, LogEta_P, LogEs_T, LogTs_T, LogEtas_T, &
             OS_NES(1,iH1), NES_T(:,:,:,:,iH1,1), Phi_In )

    ! --- Interpolate HII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Point &
           ( LogE_P, LogT_P, LogEta_P, LogEs_T, LogTs_T, LogEtas_T, &
             OS_NES(1,iH2), NES_T(:,:,:,:,iH2,1), Phi_Out )

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE2
        Phi_Out(iE1,iE2) = C1 * Phi_In(iE1,iE2) + C2 * Phi_Out(iE1,iE2)
      END DO
    END DO

    DO iE2 = iE_B, iE_E
      DO iE1 = iE2+1, iE_E
        Phi_Out(iE1,iE2) = Phi_Out(iE2,iE1) * EXP( ( E(iE2) - E(iE1) ) / kT )
      END DO
    END DO

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E
        Phi_Out(iE1,iE2) = Phi_Out(iE1,iE2) * UnitNES
      END DO
    END DO

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E
        Phi_In(iE1,iE2) = Phi_Out(iE2,iE1)
      END DO
    END DO

#else

    Phi_In  = Zero
    Phi_Out = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_NES_Point


  SUBROUTINE ComputeNeutrinoOpacities_NES_Points &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out )

    ! --- Neutrino-Electron Scattering Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Phi_Out(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    INTEGER  :: iX, iE1, iE2, iH1, iH2
    REAL(DP) :: C1, C2, kT
    REAL(DP) :: LogE_P(iE_B:iE_E), LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
    LOGICAL  :: do_gpu

    do_gpu = QueryOnGPU( E, D, T, Y ) &
       .AND. QueryOnGPU( Phi_In, Phi_Out )
#if defined(THORNADO_DEBUG_OPACITY) && defined(THORNADO_GPU)
    IF ( .not. do_gpu ) THEN
      WRITE(*,*) '[ComputeNeutrinoOpacities_NES_Points] Data not present on device'
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

    IF(     iSpecies == iNuE     )THEN

      C1 = ( cv + ca )**2
      C2 = ( cv - ca )**2

    ELSEIF( iSpecies == iNuE_Bar )THEN

      C1 = ( cv - ca )**2
      C2 = ( cv + ca )**2

    END IF

    iH1 = ( iMoment - 1 ) * 2 + 1
    iH2 = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogT_P, LogEta_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogT_P, LogEta_P )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( kT )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT ) &
    !$ACC PRESENT( T, LogT_P, LogEta_P )
#endif
    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )

      kT = BoltzmannConstant * T(iX)
      LogEta_P(iX) = LOG10( LogEta_P(iX) / kT / UnitEta )
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, LogE_P )
#endif
    DO iE1 = iE_B, iE_E
      LogE_P(iE1) = LOG10( E(iE1) / UnitE )
    END DO

    ! --- Interpolate HI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom &
           ( LogE_P, LogT_P, LogEta_P, LogEs_T, LogTs_T, LogEtas_T, &
             OS_NES(1,iH1), NES_T(:,:,:,:,iH1,1), Phi_In, &
             GPU_Option = do_gpu )

    ! --- Interpolate HII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom &
           ( LogE_P, LogT_P, LogEta_P, LogEs_T, LogTs_T, LogEtas_T, &
             OS_NES(1,iH2), NES_T(:,:,:,:,iH2,1), Phi_Out, &
             GPU_Option = do_gpu )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP PRIVATE( kT ) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT ) &
    !$ACC PRESENT( E, T, Phi_In, Phi_Out )
#endif
    DO iX = iX_B, iX_E

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR
#endif
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE2
          Phi_Out(iE1,iE2,iX) = C1 * Phi_In(iE1,iE2,iX) + C2 * Phi_Out(iE1,iE2,iX)
        END DO
      END DO

      kT = BoltzmannConstant * T(iX)

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO SIMD
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR
#endif
      DO iE2 = iE_B, iE_E
        DO iE1 = iE2+1, iE_E
          Phi_Out(iE1,iE2,iX) = Phi_Out(iE2,iE1,iX) * EXP( ( E(iE2) - E(iE1) ) / kT )
        END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR COLLAPSE(2)
#endif
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          Phi_Out(iE1,iE2,iX) = Phi_Out(iE1,iE2,iX) * UnitNES
        END DO
      END DO

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR COLLAPSE(2)
#endif
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          Phi_In(iE1,iE2,iX) = Phi_Out(iE2,iE1,iX)
        END DO
      END DO

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogT_P, LogEta_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogT_P, LogEta_P )
#endif

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

  END SUBROUTINE ComputeNeutrinoOpacities_NES_Points


  SUBROUTINE ComputeNeutrinoOpacities_Pair_Point &
    ( iE_B, iE_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    ! --- Pair Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In (iE_B:iE_E,iE_B:iE_E)
    REAL(DP), INTENT(out) :: Phi_Out(iE_B:iE_E,iE_B:iE_E)

    INTEGER  :: iE1, iE2, iJ1, iJ2
    REAL(DP) :: C1, C2, kT, Phi_Local
    REAL(DP) :: LogE_P(iE_B:iE_E), LogT_P, LogEta_P

#ifdef MICROPHYSICS_WEAKLIB

    iJ1 = ( iMoment - 1 ) * 2 + 1
    iJ2 = ( iMoment - 1 ) * 2 + 2

    IF(     iSpecies == iNuE     )THEN

      C1 = ( cv + ca )**2
      C2 = ( cv - ca )**2

    ELSEIF( iSpecies == iNuE_Bar )THEN

      C1 = ( cv - ca )**2
      C2 = ( cv + ca )**2

    END IF

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

    kT = BoltzmannConstant * T

    LogT_P = LOG10( T / UnitT )
    LogEta_P = LOG10( LogEta_P / kT / UnitEta )

    DO iE1 = iE_B, iE_E
      LogE_P(iE1) = LOG10( E(iE1) / UnitE )
    END DO

    ! --- Interpolate JI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Point &
           ( LogE_P, LogT_P, LogEta_P, LogEs_T, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ1), Pair_T(:,:,:,:,iJ1,1), Phi_In )

    ! --- Interpolate JII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom_Point &
           ( LogE_P, LogT_P, LogEta_P, LogEs_T, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ2), Pair_T(:,:,:,:,iJ2,1), Phi_Out )

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        IF ( iE1 <= iE2 ) THEN
          Phi_Local = ( C1 * Phi_In (iE1,iE2) + C2 * Phi_Out(iE1,iE2) ) * UnitPair
        ELSE
          Phi_Local = ( C1 * Phi_Out(iE2,iE1) + C2 * Phi_In (iE2,iE1) ) * UnitPair
        END IF

        Phi_Out(iE1,iE2) = Phi_Local
        Phi_In (iE1,iE2) = Phi_Local * EXP( - ( E(iE1) + E(iE2) ) / kT )

      END DO
    END DO

#else

    Phi_In  = Zero
    Phi_Out = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_Pair_Point


  SUBROUTINE ComputeNeutrinoOpacities_Pair_Points &
    ( iE_B, iE_E, iX_B, iX_E, E, D, T, Y, iSpecies, iMoment, Phi_In, Phi_Out )

    ! --- Pair Opacities (Multiple D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    INTEGER,  INTENT(in)  :: iX_B, iX_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: T(iX_B:iX_E)
    REAL(DP), INTENT(in)  :: Y(iX_B:iX_E)
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In (iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)
    REAL(DP), INTENT(out) :: Phi_Out(iE_B:iE_E,iE_B:iE_E,iX_B:iX_E)

    INTEGER  :: iX, iE1, iE2, iJ1, iJ2
    REAL(DP) :: C1, C2, kT, Phi_Local(iE_B:iE_E,iE_B:iE_E)
    REAL(DP) :: LogE_P(iE_B:iE_E), LogT_P(iX_B:iX_E), LogEta_P(iX_B:iX_E)
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

    IF(     iSpecies == iNuE     )THEN

      C1 = ( cv + ca )**2
      C2 = ( cv - ca )**2

    ELSEIF( iSpecies == iNuE_Bar )THEN

      C1 = ( cv - ca )**2
      C2 = ( cv + ca )**2

    END IF

    iJ1 = ( iMoment - 1 ) * 2 + 1
    iJ2 = ( iMoment - 1 ) * 2 + 2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( alloc: LogE_P, LogT_P, LogEta_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( LogE_P, LogT_P, LogEta_P )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, LogEta_P )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( T, LogT_P, LogEta_P )
#endif
    DO iX = iX_B, iX_E
      LogT_P(iX) = LOG10( T(iX) / UnitT )
      LogEta_P(iX) = LOG10( LogEta_P(iX) / ( BoltzmannConstant * T(iX) ) / UnitEta )
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, LogE_P )
#endif
    DO iE1 = iE_B, iE_E
      LogE_P(iE1) = LOG10( E(iE1) / UnitE )
    END DO

    ! --- Interpolate JI  (put result temporarily in Phi_In) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom &
           ( LogE_P, LogT_P, LogEta_P, LogEs_T, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ1), Pair_T(:,:,:,:,iJ1,1), Phi_In, &
             GPU_Option = do_gpu )

    ! --- Interpolate JII (put result temporarily in Phi_Out) ---

    CALL LogInterpolateSingleVariable_2D2D_Custom &
           ( LogE_P, LogT_P, LogEta_P, LogEs_T, LogTs_T, LogEtas_T, &
             OS_Pair(1,iJ2), Pair_T(:,:,:,:,iJ2,1), Phi_Out, &
             GPU_Option = do_gpu )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( kT, Phi_Local )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( kT, Phi_Local ) &
    !$ACC PRESENT( E, T, Phi_In, Phi_Out )
#endif
    DO iX = iX_B, iX_E

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR COLLAPSE(2)
#endif
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          IF ( iE1 <= iE2 ) THEN
            Phi_Local(iE1,iE2) = ( C1 * Phi_In (iE1,iE2,iX) + C2 * Phi_Out(iE1,iE2,iX) ) * UnitPair
          ELSE
            Phi_Local(iE1,iE2) = ( C1 * Phi_Out(iE2,iE1,iX) + C2 * Phi_In (iE2,iE1,iX) ) * UnitPair
          END IF
        END DO
      END DO

      kT = BoltzmannConstant * T(iX)

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO SIMD COLLAPSE(2)
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR COLLAPSE(2)
#endif
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E
          Phi_Out(iE1,iE2,iX) = Phi_Local(iE1,iE2)
          Phi_In (iE1,iE2,iX) = Phi_Local(iE1,iE2) * EXP( - ( E(iE1) + E(iE2) ) / kT )
        END DO
      END DO

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: LogE_P, LogT_P, LogEta_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( LogE_P, LogT_P, LogEta_P )
#endif

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

  END SUBROUTINE ComputeNeutrinoOpacities_Pair_Points


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


  PURE ELEMENTAL REAL(DP) FUNCTION FermiDirac( E, Mu, kT )

    REAL(DP), INTENT(in) :: E, Mu, kT

    REAL(DP) :: Exponent

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FermiDirac &
      = One / ( EXP( Exponent ) + One )

    RETURN
  END FUNCTION FermiDirac


  PURE ELEMENTAL REAL(DP) FUNCTION dFermiDiracdT( E, Mu, kT, dMudT, T )

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudT, T

    REAL(DP) :: Exponent, FD

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdT &
      = ( FD * EXP( Exponent ) ) * FD * ( dMudT + ( E - Mu ) / T ) / kT

    RETURN
  END FUNCTION dFermiDiracdT


  PURE ELEMENTAL REAL(DP) FUNCTION dFermiDiracdY( E, Mu, kT, dMudY, T )

    REAL(DP), INTENT(in) :: E, Mu, kT, dMudY, T

    REAL(DP) :: Exponent, FD

    Exponent = MIN( MAX( ( E - Mu ) / kT, - Log1d100 ), + Log1d100 )

    FD = FermiDirac( E, Mu, kT )

    dFermiDiracdY &
      = ( FD * EXP( Exponent ) ) * FD * dMudY / kT

    RETURN
  END FUNCTION dFermiDiracdY


END MODULE NeutrinoOpacitiesComputationModule
