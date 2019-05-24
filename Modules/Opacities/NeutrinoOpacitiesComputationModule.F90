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
  REAL(DP), PARAMETER :: UnitEC   = One / Centimeter
  REAL(DP), PARAMETER :: UnitES   = One / Centimeter
  REAL(DP), PARAMETER :: UnitNES  = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: UnitPair = One / ( Centimeter * MeV**3 )
  REAL(DP), PARAMETER :: cv       = 0.96d+00 ! weak interaction constant
  REAL(DP), PARAMETER :: ca       = 0.50d+00 ! weak interaction constant

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

    ! --- Equilibrium Neutrino Distributions (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: f_EQ_Point(iE_B:iE_E)
    INTEGER,  INTENT(in)  :: iSpecies

    REAL(DP) :: Me, Mp, Mn, Mnu

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

      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'ERROR (ComputeEquilibriumDistributions_Point):'
      WRITE(*,'(A4,A,I2.2)') '', 'Invalid Species: ', iSpecies
      WRITE(*,*)
      STOP

    END IF

    f_EQ_Point = FermiDirac( E, Mnu, BoltzmannConstant * T )

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

    ! --- Electron Capture Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: opEC_Point(iE_B:iE_E)
    INTEGER,  INTENT(in)  :: iSpecies

    CALL ComputeNeutrinoOpacityE_Point &
           ( E, D, T, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             opEC_Point, EmAb_T(:,:,:,:,iSpecies), OS_EmAb(iSpecies), UnitEC )

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

    INTEGER :: iX, iE
    LOGICAL :: do_gpu

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, D, T, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, opEC_Points, OS_EmAb, EmAb_T )
#endif
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE), D(iX), T(iX), Y(iX), LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                 opEC_Points(iE,iX), EmAb_T(:,:,:,:,iSpecies), OS_EmAb(iSpecies), UnitEC )

      END DO
    END DO

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
                 Iso_T(:,1,:,:,:,iSpecies), &
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

    ! --- Elastic Scattering Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: opES_Point(iE_B:iE_E)

    CALL ComputeNeutrinoOpacityE_Point &
           ( E, D, T, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
             opES_Point, Iso_T(:,iMoment,:,:,:,iSpecies), OS_Iso(iSpecies,iMoment), UnitES )

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

    INTEGER :: iX, iE
    LOGICAL :: do_gpu

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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP IF( do_gpu )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC IF( do_gpu ) &
    !$ACC PRESENT( E, D, T, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, opES_Points, OS_Iso, Iso_T )
#endif
    DO iX = iX_B, iX_E
      DO iE = iE_B, iE_E

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE), D(iX), T(iX), Y(iX), LogEs_T, LogDs_T, LogTs_T, Ys_T, &
                 opES_Points(iE,iX), Iso_T(:,iMoment,:,:,:,iSpecies), OS_Iso(iSpecies,iMoment), UnitES )

      END DO
    END DO

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
    REAL(DP) :: C1, C2, H1, H2, H1_P, H2_P, Phi, Phi_P, kT, dE
    REAL(DP) :: Me, LogEta

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
           ( D, T, Y, Me )

    kT = BoltzmannConstant * T
    LogEta = LOG10( Me / kT )

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        ! --- Interpolate HI ---

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE1), E(iE2), T, LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                 H1, NES_T(:,:,iH1,:,:,1), OS_NES(1,iH1), UnitNES )

        ! --- Interpolate HI' ---

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE2), E(iE1), T, LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                 H1_P, NES_T(:,:,iH1,:,:,1), OS_NES(1,iH1), UnitNES )

        ! --- Interpolate HII ---

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE1), E(iE2), T, LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                 H2, NES_T(:,:,iH2,:,:,1), OS_NES(1,iH2), UnitNES )

        ! --- Interpolate HII' ---

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE2), E(iE1), T, LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                 H2_P, NES_T(:,:,iH2,:,:,1), OS_NES(1,iH2), UnitNES )

        ! --- Compute Phi ---

        Phi = C1 * H1 + C2 * H2

        ! --- Enforce Detailed Balance ---

        dE = E(iE2) - E(iE1)
        Phi_P = ( C1 * H1_P + C2 * H2_P ) * EXP( dE / kT )

        IF ( iE1 <= iE2 ) THEN

          Phi_Out(iE1,iE2) = Phi
          Phi_In (iE1,iE2) = Phi_P

        ELSE

          Phi_Out(iE1,iE2) = Phi_P
          Phi_In (iE1,iE2) = Phi

        END IF

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
    REAL(DP) :: C1, C2, H1, H2, H1_P, H2_P, Phi, Phi_P, kT, dE
    REAL(DP) :: Me(iX_B:iX_E), LogEta
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
    !$OMP MAP( alloc: Me )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( Me )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( H1, H2, H1_P, H2_P, Phi, Phi_P, LogEta, kT, dE )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( H1, H2, H1_P, H2_P, Phi, Phi_P, LogEta, kT, dE )
    !$ACC PRESENT( E, T, Me, LogEs_T, LogTs_T, LogEtas_T, OS_NES, NES_T, Phi_In, Phi_Out )
#endif
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E

          kT = BoltzmannConstant * T(iX)
          LogEta = LOG10( Me(iX) / kT )

          ! --- Interpolate HI ---

          CALL ComputeNeutrinoOpacity_Point &
                 ( E(iE1), E(iE2), T(iX), LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                   H1, NES_T(:,:,iH1,:,:,1), OS_NES(1,iH1), UnitNES )

          ! --- Interpolate HI' ---

          CALL ComputeNeutrinoOpacity_Point &
                 ( E(iE2), E(iE1), T(iX), LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                   H1_P, NES_T(:,:,iH1,:,:,1), OS_NES(1,iH1), UnitNES )

          ! --- Interpolate HII ---

          CALL ComputeNeutrinoOpacity_Point &
                 ( E(iE1), E(iE2), T(iX), LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                   H2, NES_T(:,:,iH2,:,:,1), OS_NES(1,iH2), UnitNES )

          ! --- Interpolate HII' ---

          CALL ComputeNeutrinoOpacity_Point &
                 ( E(iE2), E(iE1), T(iX), LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                   H2_P, NES_T(:,:,iH2,:,:,1), OS_NES(1,iH2), UnitNES )

          ! --- Compute Phi ---

          Phi = C1 * H1 + C2 * H2

          ! --- Enforce Detailed Balance ---

          dE = E(iE2) - E(iE1)
          Phi_P = ( C1 * H1_P + C2 * H2_P ) * EXP( dE / kT )

          IF ( iE1 <= iE2 ) THEN

            Phi_Out(iE1,iE2,iX) = Phi
            Phi_In (iE1,iE2,iX) = Phi_P

          ELSE

            Phi_Out(iE1,iE2,iX) = Phi_P
            Phi_In (iE1,iE2,iX) = Phi

          END IF

        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: Me )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( Me )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
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

    ! --- Pair Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    INTEGER,  INTENT(in)  :: iSpecies
    INTEGER,  INTENT(in)  :: iMoment
    REAL(DP), INTENT(out) :: Phi_In (iE_B:iE_E,iE_B:iE_E)
    REAL(DP), INTENT(out) :: Phi_Out(iE_B:iE_E,iE_B:iE_E)

    INTEGER  :: iE1, iE2, iJ1, iJ2
    REAL(DP) :: C1, C2, J1, J2, J1_P, J2_P, Phi, Phi_P, kT, E12
    REAL(DP) :: Me, LogEta

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
           ( D, T, Y, Me )

    kT = BoltzmannConstant * T
    LogEta = LOG10( Me / kT )

    DO iE2 = iE_B, iE_E
      DO iE1 = iE_B, iE_E

        ! --- Interpolate JI ---

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE1), E(iE2), T, LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                 J1, Pair_T(:,:,iJ1,:,:,1), OS_Pair(1,iJ1), UnitPair )

        ! --- Interpolate JI' ---

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE2), E(iE1), T, LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                 J1_P, Pair_T(:,:,iJ1,:,:,1), OS_Pair(1,iJ1), UnitPair )

        ! --- Interpolate JII ---

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE1), E(iE2), T, LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                 J2, Pair_T(:,:,iJ2,:,:,1), OS_Pair(1,iJ2), UnitPair )

        ! --- Interpolate JII' ---

        CALL ComputeNeutrinoOpacity_Point &
               ( E(iE2), E(iE1), T, LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                 J2_P, Pair_T(:,:,iJ2,:,:,1), OS_Pair(1,iJ2), UnitPair )

        ! --- Compute Phi ---

        Phi   = C1 * J1   + C2 * J2
        Phi_P = C1 * J2_P + C2 * J1_P

        E12 = E(iE2) + E(iE1)

        IF ( iE1 <= iE2 ) THEN

          Phi_Out(iE1,iE2) = Phi
          Phi_In (iE1,iE2) = Phi_P * EXP( - E12 / kT )

        ELSE

          Phi_Out(iE1,iE2) = Phi_P
          Phi_In (iE1,iE2) = Phi   * EXP( - E12 / kT )

        END IF

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
    REAL(DP) :: C1, C2, J1, J2, J1_P, J2_P, Phi, Phi_P, kT, E12
    REAL(DP) :: Me(iX_B:iX_E), LogEta
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
    !$OMP MAP( alloc: Me )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC IF( do_gpu ) &
    !$ACC CREATE( Me )
#endif

    ! --- Compute Electron Chemical Potential ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( D, T, Y, Me )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
    !$OMP IF( do_gpu ) &
    !$OMP PRIVATE( J1, J2, J1_P, J2_P, Phi, Phi_P, LogEta, kT, dE )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC IF( do_gpu ) &
    !$ACC PRIVATE( J1, J2, J1_P, J2_P, Phi, Phi_P, LogEta, kT, dE )
    !$ACC PRESENT( E, T, Me, LogEs_T, LogTs_T, LogEtas_T, OS_Pair, Pair_T, Phi_In, Phi_Out )
#endif
    DO iX = iX_B, iX_E
      DO iE2 = iE_B, iE_E
        DO iE1 = iE_B, iE_E

          kT = BoltzmannConstant * T(iX)
          LogEta = LOG10( Me(iX) / kT )

          ! --- Interpolate JI ---

          CALL ComputeNeutrinoOpacity_Point &
                 ( E(iE1), E(iE2), T(iX), LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                   J1, Pair_T(:,:,iJ1,:,:,1), OS_Pair(1,iJ1), UnitPair )

          ! --- Interpolate JI' ---

          CALL ComputeNeutrinoOpacity_Point &
                 ( E(iE2), E(iE1), T(iX), LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                   J1_P, Pair_T(:,:,iJ1,:,:,1), OS_Pair(1,iJ1), UnitPair )

          ! --- Interpolate JII ---

          CALL ComputeNeutrinoOpacity_Point &
                 ( E(iE1), E(iE2), T(iX), LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                   J2, Pair_T(:,:,iJ2,:,:,1), OS_Pair(1,iJ2), UnitPair )

          ! --- Interpolate JII' ---

          CALL ComputeNeutrinoOpacity_Point &
                 ( E(iE2), E(iE1), T(iX), LogEta, LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
                   J2_P, Pair_T(:,:,iJ2,:,:,1), OS_Pair(1,iJ2), UnitPair )

          ! --- Compute Phi ---

          Phi   = C1 * J1   + C2 * J2
          Phi_P = C1 * J2_P + C2 * J1_P

          E12 = E(iE2) + E(iE1)

          IF ( iE1 <= iE2 ) THEN

            Phi_Out(iE1,iE2,iX) = Phi
            Phi_In (iE1,iE2,iX) = Phi_P * EXP( - E12 / kT )

          ELSE

            Phi_Out(iE1,iE2,iX) = Phi_P
            Phi_In (iE1,iE2,iX) = Phi   * EXP( - E12 / kT )

          END IF

        END DO
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP IF( do_gpu ) &
    !$OMP MAP( release: Me )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC IF( do_gpu ) &
    !$ACC DELETE( Me )
#endif

#else

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
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


  SUBROUTINE ComputeNeutrinoOpacity_Point &
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

  END SUBROUTINE ComputeNeutrinoOpacity_Point


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
