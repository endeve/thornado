MODULE FluidRadiationCouplingSolutionModule_ThermalReservoir

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE FluidFieldsModule, ONLY: &
    nPF, iPF_D, &
    nAF, iAF_T, iAF_Ye, iAF_Me, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputePrimitiveMoments, &
    ComputeConservedMoments
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeNodesX, &
    InitializeFluidFields, &
    InitializeRadiationFields, &
    FinalizeFluidFields, &
    FinalizeRadiationFields, &
    FermiDirac
  USE OpacityModule, ONLY: &
    ComputeAbsorptionOpacity

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nNodesX_G, nNodesE_G
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: X_N, uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Chi, FD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N

  PUBLIC :: CoupleFluidRadiation_ThermalReservoir

CONTAINS


  SUBROUTINE CoupleFluidRadiation_ThermalReservoir &
               ( dt, iX_B0, iX_E0, iX_B1, iX_E1, U_F, dU_F, U_R, dU_R, &
                 EvolveFluid_Option )

    REAL(DP), INTENT(in)  :: &
      dt
    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U_F (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      dU_F(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in)  :: &
      U_R (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    REAL(DP), INTENT(out) :: &
      dU_R(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      EvolveFluid_Option

    CALL ComputePrimitiveMoments &
           ( iX_Begin = iX_B0, iX_End = iX_E0 )

    CALL InitializeFluidRadiationCoupling

    CALL CoupleFluidRadiation( dt )

    CALL FinalizeFluidRadiationCoupling

    CALL ComputeConservedMoments &
           ( iX_Begin = iX_B0, iX_End = iX_E0 )

  END SUBROUTINE CoupleFluidRadiation_ThermalReservoir


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( X_N(nNodesX_G,3) )
    CALL InitializeNodesX( X_N )

    ALLOCATE( uPF_N(nPF, nNodesX_G) )
    ALLOCATE( uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    ALLOCATE &
      ( Chi(nNodesE_G, nNodesX_G), &
        FD (nNodesE_G, nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    DEALLOCATE( E_N, X_N )

    CALL FinalizeFluidFields( uPF_N, uAF_N )
    DEALLOCATE( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )
    DEALLOCATE( uPR_N )

    DEALLOCATE( Chi, FD )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER  :: iX, iE
    REAL(DP) :: Gamma

    CALL SetRates

    CALL SetEquilibrium

    DO iX = 1, nNodesX_G

      DO iE = 1, nNodesE_G

        Gamma = dt * Chi(iE,iX)

        ! --- Number Density ---

        uPR_N(iE,iPR_D,iX) &
          = ( uPR_N(iE,iPR_D,iX) &
              + Gamma * FourPi * FD(iE,iX) ) / ( 1.0_DP + Gamma )

        ! --- Number Flux Density (1) ---

        uPR_N(iE,iPR_I1,iX) &
          = uPR_N(iE,iPR_I1,iX) / ( 1.0_DP + Gamma )

        ! --- Number Flux Density (2) ---

        uPR_N(iE,iPR_I2,iX) &
          = uPR_N(iE,iPR_I2,iX) / ( 1.0_DP + Gamma )

        ! --- Number Flux Density (3) ---

        uPR_N(iE,iPR_I3,iX) &
          = uPR_N(iE,iPR_I3,iX) / ( 1.0_DP + Gamma )

      END DO

    END DO

  END SUBROUTINE CoupleFluidRadiation


  SUBROUTINE SetRates

    ASSOCIATE &
      ( D_N => uPF_N(iPF_D, 1:nNodesX_G), &
        T_N => uAF_N(iAF_T, 1:nNodesX_G), &
        Y_N => uAF_N(iAF_Ye,1:nNodesX_G) )

    CALL ComputeAbsorptionOpacity &
           ( E_N, D_N, T_N, Y_N, X_N(:,1), X_N(:,2), X_N(:,3), Chi )

    END ASSOCIATE ! D_N, etc.

  END SUBROUTINE SetRates


  SUBROUTINE SetEquilibrium

    INTEGER  :: iX
    REAL(DP) :: Mnu

    ASSOCIATE &
      ( T  => uAF_N(iAF_T, :), &
        Me => uAF_N(iAF_Me,:), &
        Mp => uAF_N(iAF_Mp,:), &
        Mn => uAF_N(iAF_Mn,:) )

    DO iX = 1, nNodesX_G

      Mnu = Me(iX) + Mp(iX) - Mn(iX)

      FD(:,iX) = FermiDirac( E_N, Mnu, BoltzmannConstant * T(iX) )

    END DO

    END ASSOCIATE ! T, etc.

  END SUBROUTINE SetEquilibrium


END MODULE FluidRadiationCouplingSolutionModule_ThermalReservoir
