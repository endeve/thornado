MODULE FluidRadiationCouplingSolutionModule_ElasticScattering

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE FluidFieldsModule, ONLY: &
    nPF, iPF_D, &
    nAF, iAF_T, iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputeConservedMoments, &
    ComputePrimitiveMoments
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeFluidFields, &
    InitializeRadiationFields, &
    FinalizeFluidFields, &
    FinalizeRadiationFields
  USE OpacityModule, ONLY: &
    ComputeScatteringOpacity_ES

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nNodesX_G, nNodesE_G
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Kappa
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N

  PUBLIC :: CoupleFluidRadiation_ElasticScattering

CONTAINS


  SUBROUTINE CoupleFluidRadiation_ElasticScattering &
               ( dt, iX_Begin, iX_End, EvolveFluid_Option )

    REAL(DP),              INTENT(in) :: dt
    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    LOGICAL,               INTENT(in), OPTIONAL :: EvolveFluid_Option

    CALL ComputePrimitiveMoments &
           ( iX_Begin = iX_Begin, iX_End = iX_End )

    CALL InitializeFluidRadiationCoupling

    CALL CoupleFluidRadiation( dt )

    CALL FinalizeFluidRadiationCoupling

    CALL ComputeConservedMoments &
           ( iX_Begin = iX_Begin, iX_End = iX_End )

  END SUBROUTINE CoupleFluidRadiation_ElasticScattering


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( uPF_N(nPF, nNodesX_G) )
    ALLOCATE( uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    ALLOCATE( Kappa(nNodesE_G, nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    DEALLOCATE( E_N )

    CALL FinalizeFluidFields( uPF_N, uAF_N )
    DEALLOCATE( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )
    DEALLOCATE( uPR_N )

    DEALLOCATE( Kappa )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER :: iX, iE

    CALL SetRates

    DO iX = 1, nNodesX_G

      DO iE = 1, nNodesE_G

        ! --- Number Flux Density (1) ---

        uPR_N(iE,iPR_I1,iX) &
          = uPR_N(iE,iPR_I1,iX) / ( 1.0_DP + dt * Kappa(iE,iX) )

        ! --- Number Flux Density (2) ---

        uPR_N(iE,iPR_I2,iX) &
          = uPR_N(iE,iPR_I2,iX) / ( 1.0_DP + dt * Kappa(iE,iX) )

        ! --- Number Flux Density (3) ---

        uPR_N(iE,iPR_I3,iX) &
          = uPR_N(iE,iPR_I3,iX) / ( 1.0_DP + dt * Kappa(iE,iX) )

      END DO

    END DO

  END SUBROUTINE CoupleFluidRadiation


  SUBROUTINE SetRates

    ASSOCIATE &
      ( D_N => uPF_N(iPF_D, 1:nNodesX_G), &
        T_N => uAF_N(iAF_T, 1:nNodesX_G), &
        Y_N => uAF_N(iAF_Ye,1:nNodesX_G) )

!!$      Kappa = 1.0d-1 * ( 1.0_DP / Centimeter )

    CALL ComputeScatteringOpacity_ES( E_N, D_N, T_N, Y_N, Kappa )

    END ASSOCIATE ! D_N, etc.

  END SUBROUTINE SetRates


END MODULE FluidRadiationCouplingSolutionModule_ElasticScattering
