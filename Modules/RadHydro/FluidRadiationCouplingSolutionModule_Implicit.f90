MODULE FluidRadiationCouplingSolutionModule_Implicit

  USE KindModule, ONLY: &
    DP, Pi
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    Centimeter, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE FluidFieldsModule, ONLY: &
    iPF_D, nPF, &
    iAF_T, iAF_Ye, iAF_Me, iAF_Mp, iAF_Mn, nAF
  USE RadiationFieldsModule, ONLY: &
    iPR_D, nPR
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeWeights, &
    InitializeFluidFields, &
    FinalizeFluidFields, &
    InitializeRadiationFields, &
    FinalizeRadiationFields
  USE EquationOfStateModule, ONLY: &
    ComputeChemicalPotentials_TABLE
  USE OpacityModule, ONLY: &
    ComputeAbsorptionCoefficients

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nNodesX_G, nNodesE_G
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Chi, f_FD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N

  PUBLIC :: CoupleFluidRadiation_Implicit_EmissionAbsorption

CONTAINS


  SUBROUTINE CoupleFluidRadiation_Implicit_EmissionAbsorption &
               ( dt, iX_Begin, iX_End )

    REAL(DP),              INTENT(in) :: dt
    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    CALL InitializeFluidRadiationCoupling

    CALL CoupleFluidRadiation_EmissionAbsorption( dt )

    CALL FinalizeFluidRadiationCoupling

  END SUBROUTINE CoupleFluidRadiation_Implicit_EmissionAbsorption


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    CALL WriteVector( SIZE( E_N ), E_N / MeV, 'E.dat' )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    CALL WriteVector( SIZE( W2_N ), W2_N / MeV**3, 'W2.dat' )
    CALL WriteVector( SIZE( W3_N ), W3_N / MeV**4, 'W3.dat' )

    ALLOCATE( uPF_N(nPF, nNodesX_G), uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    CALL WriteVector( SIZE( uPR_N(:,1,1) ), uPR_N(:,1,1), 'N.dat' )

    ALLOCATE( Chi(nNodesE_G, nNodesX_G), f_FD(nNodesE_G, nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N, uPF_N, uAF_N, uPR_N, Chi, f_FD )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation_EmissionAbsorption( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER  :: iX, iE
    REAL(DP) :: Mnu, Gamma

    CALL ComputeChemicalPotentials_TABLE &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
             uAF_N(iAF_Me,:), uAF_N(iAF_Mp,:), uAF_N(iAF_Mn,:) )

    DO iX = 1, nNodesX_G

      CALL ComputeAbsorptionCoefficients &
             ( E_N, [uPF_N(iPF_D,iX)], [uAF_N(iAF_T,iX)], [uAF_N(iAF_Ye,iX)], &
               Chi(:,iX) )

      Mnu = uAF_N(iAF_Me,iX) + uAF_N(iAF_Mp,iX) - uAF_N(iAF_Mn,iX)

      f_FD(:,iX) = FermiDirac( E_N, Mnu, BoltzmannConstant * uAF_N(iAF_T,iX) )

      DO iE = 1, nNodesE_G

        Gamma = dt * Chi(iE,iX)

        uPR_N(iE,iPR_D,iX) &
          = ( uPR_N(iE,iPR_D,iX) +  Gamma * 4.0_DP * Pi * f_FD(iE,iX) ) &
              / ( 1.0_DP + Gamma )

      END DO

    END DO

    CALL WriteVector( SIZE( f_FD(:,1) ), f_FD(:,1), 'f_FD.dat' )

    CALL WriteVector &
           ( SIZE( Chi(:,1) ), Chi(:,1) / ( 1.0_DP / Centimeter ), 'Chi.dat' )

  END SUBROUTINE CoupleFluidRadiation_EmissionAbsorption


  PURE ELEMENTAL REAL(DP) FUNCTION FermiDirac( E, Mu, kT )

    REAL(DP), INTENT(in) :: E, Mu, kT

    FermiDirac = 1.0_DP / ( EXP( ( E - Mu ) / kT ) + 1.0_DP )

    RETURN
  END FUNCTION FermiDirac


END MODULE FluidRadiationCouplingSolutionModule_Implicit
