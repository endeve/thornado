MODULE NeutrinoOpacitiesComputationModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE UnitsModule, ONLY: &
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
    OPACITIES, &
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
      ( opEC_T => OPACITIES % ThermEmAb % Absorptivity(1) % Values, &
        OS     => OPACITIES % ThermEmAb % Offsets(1) )

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
                     LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, &
                     opEC_T, opEC_K )

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
      ( opES_T => OPACITIES % Scatt_Iso % Kernel(1) % Values(:,:,:,:,1), &
        OS     => OPACITIES % Scatt_Iso % Offsets(1,1) )

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
                     LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, &
                     opES_T, opES_K )

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


  PURE ELEMENTAL REAL(DP) FUNCTION FermiDirac( E, Mu, kT )

    REAL(DP), INTENT(in) :: E, Mu, kT

    REAL(DP) :: Exponent

    Exponent = ( E - Mu ) / kT
    Exponent = MAX( Exponent, - Log1d100 )
    Exponent = MIN( Exponent, + Log1d100 )

    FermiDirac &
      = One / ( EXP( Exponent ) + One )

    RETURN
  END FUNCTION FermiDirac


END MODULE NeutrinoOpacitiesComputationModule
