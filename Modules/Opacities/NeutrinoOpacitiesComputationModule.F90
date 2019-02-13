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
    OPACITIES, &
#endif
    LogEs_T, LogDs_T, LogTs_T, Ys_T
  USE NeutrinoOpacitiesModule, ONLY: &
    f_EQ, opEC, opES, opIS, opPP

#ifdef MICROPHYSICS_WEAKLIB

  ! --- weaklib modules ---

  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_1D3D_Custom

  ! ----------------------------------------------

#endif

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeNeutrinoOpacities
  PUBLIC :: ComputeNeutrinoOpacities_EC_Point
  PUBLIC :: ComputeNeutrinoOpacities_ES_Point
  PUBLIC :: ComputeEquilibriumDistributions_Point
  PUBLIC :: ComputeNeutrinoOpacities_EC_Points
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

CONTAINS


  SUBROUTINE ComputeNeutrinoOpacities( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      T(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      Y(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):)

    REAL(DP) :: wTime

    CALL ComputeEquilibriumDistributions &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y )

    wTime = MPI_WTIME( )

    CALL ComputeNeutrinoOpacities_EC &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y )

    wTime = MPI_WTIME( ) - wTime

    PRINT*
    PRINT*, "  ComputeNeutrinoOpacities:"
    PRINT*
    PRINT*, "    EC: ", wTime

    wTime = MPI_WTIME( )

    CALL ComputeNeutrinoOpacities_ES &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y )

    wTime = MPI_WTIME( ) - wTime

    PRINT*, "    ES: ", wTime

  END SUBROUTINE ComputeNeutrinoOpacities


  SUBROUTINE ComputeEquilibriumDistributions &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y )

    ! --- Equilibrium Neutrino Distributions ---

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      T(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      Y(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):)

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

          ! --- Electron Neutrinos ---

          Mnu = Me + Mp - Mn

          ! --- Offset (Position) ---

          iOS_X = ( (iZ4-1)*nZ(3)*nZ(2) + (iZ3-1)*nZ(2) + (iZ2-1) ) * nDOFX

          DO iZ1 = iZ_B0(1), iZ_E0(1)

            ! --- Offset (Energy) ---

            iOS_E = (iZ1-1) * nDOFE

            ! --- Energy Coordinates ---

            E = MeshE % Center(iZ1) + MeshE % Width(iZ1) * MeshE % Nodes(:)

            DO iNodeX = 1, nDOFX
              DO iNodeE = 1, nDOFE

                f_EQ(iOS_E+iNodeE,1,iOS_X+iNodeX) &
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

    REAL(DP) :: Me(1), Mp(1), Mn(1), Mnu

    ! --- Compute Chemical Potentials ---

    CALL ComputeElectronChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], Me )

    CALL ComputeProtonChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], Mp )

    CALL ComputeNeutronChemicalPotential_TABLE &
           ( [ D ], [ T ], [ Y ], Mn )

    IF    ( iSpecies .EQ. 1 )THEN

      ! --- Electron Neutrinos ---

      Mnu = ( Me(1) + Mp(1) ) - Mn(1)

    ELSEIF( iSpecies .EQ. 2 )THEN

      ! --- Electron Antineutrino ---

      Mnu = Mn(1) - ( Me(1) + Mp(1) )

    ELSE

      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'ERROR (ComputeEquilibriumDistributions_Point):'
      WRITE(*,'(A4,A,I2.2)') '', 'Unknown Species: ', iSpecies
      WRITE(*,*)

    END IF

    f_EQ_Point = FermiDirac( E, Mnu, BoltzmannConstant * T )

  END SUBROUTINE ComputeEquilibriumDistributions_Point


  SUBROUTINE ComputeNeutrinoOpacities_EC( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y )

    ! --- Electron Capture Opacities ---

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      T(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      Y(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, nZ(4)
    INTEGER  :: iOS_X, iOS_E
    REAL(DP) :: D_K(1), T_K(1), Y_K(1), E_K(1), opEC_K(1,1)

    nZ = iZ_E0 - iZ_B0 + 1

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( opEC_T => OPACITIES % EmAb % Opacity(1) % Values, &
        OS     => OPACITIES % EmAb % Offsets(1) )

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

        ! --- Electron Neutrinos ---

        CALL LogInterpolateSingleVariable_1D3D_Custom &
               ( LOG10( E_K / UnitE ), LOG10( D_K / UnitD ), &
                 LOG10( T_K / UnitT ),      ( Y_K / UnitY ), &
                 LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, opEC_T, opEC_K )

        opEC(iOS_E+1:iOS_E+nDOFE,1,iOS_X+1:iOS_X+nDOFX) &
          = opEC_K(1,1) * UnitEC

      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE ! opEC_T, etc.

#else

!    opEC = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC


  SUBROUTINE ComputeNeutrinoOpacities_EC_Point &
    ( iE_B, iE_E, E, D, T, Y, opEC_Point, iSpecies )

    ! --- Electron Capture Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: opEC_Point(iE_B:iE_E)
    INTEGER,  INTENT(in)  :: iSpecies

    REAL(DP) :: tmp(iE_B:iE_E,1)

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( opEC_T => OPACITIES % EmAb &
                    % Opacity(iSpecies) % Values, &
        OS     => OPACITIES % EmAb &
                    % Offsets(iSpecies) )

    CALL LogInterpolateSingleVariable_1D3D_Custom &
           ( LOG10( [ E ] / UnitE ), LOG10( [ D ] / UnitD ), &
             LOG10( [ T ] / UnitT ),      ( [ Y ] / UnitY ), &
             LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, opEC_T, tmp )

    opEC_Point(:) = tmp(:,1) * UnitEC

    END ASSOCIATE ! opEC_T, etc.

#else

    opEC_Point(:) = Zero

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

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( opEC_T => OPACITIES % EmAb &
                    % Opacity(iSpecies) % Values, &
        OS     => OPACITIES % EmAb &
                    % Offsets(iSpecies) )

    CALL LogInterpolateSingleVariable_1D3D_Custom &
           ( LOG10( E / UnitE ), LOG10( D / UnitD ), &
             LOG10( T / UnitT ),      ( Y / UnitY ), &
             LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, opEC_T, &
             opEC_Points )

    opEC_Points = opEC_Points * UnitEC

    END ASSOCIATE ! opEC_T, etc.

#else

    opEC_Points = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_EC_Points


  SUBROUTINE ComputeNeutrinoOpacities_ES( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D, T, Y )

    ! --- Elastic Scattering Opacities ---

    INTEGER, INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      T(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):), &
      Y(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, nZ(4)
    INTEGER  :: iOS_X, iOS_E
    REAL(DP) :: D_K(1), T_K(1), Y_K(1), E_K(1), opES_K(1,1)

    nZ = iZ_E0 - iZ_B0 + 1

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( opES_T => OPACITIES % Scat_Iso % Kernel(1) % Values(:,1,:,:,:), &
        OS     => OPACITIES % Scat_Iso % Offsets(1,1) )

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

        ! --- Electron Neutrinos ---

        CALL LogInterpolateSingleVariable_1D3D_Custom &
               ( LOG10( E_K / UnitE ), LOG10( D_K / UnitD ), &
                 LOG10( T_K / UnitT ),      ( Y_K / UnitY ), &
                 LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, opES_T, opES_K )

        opES(iOS_E+1:iOS_E+nDOFE,1,iOS_X+1:iOS_X+nDOFX) &
          = opES_K(1,1) * UnitES

      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE ! opES_T, etc.

#else

!    opES = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES


  SUBROUTINE ComputeNeutrinoOpacities_ES_Point &
    ( iE_B, iE_E, E, D, T, Y, opES_Point, iSpecies )

    ! --- Elastic Scattering Opacities (Single D,T,Y) ---

    INTEGER,  INTENT(in)  :: iE_B, iE_E
    REAL(DP), INTENT(in)  :: E(iE_B:iE_E)
    REAL(DP), INTENT(in)  :: D, T, Y
    REAL(DP), INTENT(out) :: opES_Point(iE_B:iE_E)
    INTEGER,  INTENT(in)  :: iSpecies

    REAL(DP) :: tmp(iE_B:iE_E,1)

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( opES_T => OPACITIES % Scat_Iso &
                    % Kernel(iSpecies) % Values(:,1,:,:,:), &
        OS     => OPACITIES % Scat_Iso &
                    % Offsets(iSpecies,1) )

    CALL LogInterpolateSingleVariable_1D3D_Custom &
           ( LOG10( [ E ] / UnitE ), LOG10( [ D ] / UnitD ), &
             LOG10( [ T ] / UnitT ),      ( [ Y ] / UnitY ), &
             LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, opES_T, tmp )

    opES_Point(:) = tmp(:,1) * UnitES

    END ASSOCIATE ! opEC_T, etc.

#else

    opES_Point(:) = Zero

#endif

  END SUBROUTINE ComputeNeutrinoOpacities_ES_Point


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
