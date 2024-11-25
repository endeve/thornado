MODULE FluidRadiationCouplingSolutionModule_ConstantOpacities

  USE KindModule, ONLY: &
    DP, Zero
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE, &
    nDOF
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, &
    uAF, iAF_T, iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    rhsCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputePrimitiveMoments
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeNodesX, &
    MapForward_FluidField, &
    MapForward_RadiationField, &
    MapBackward_RadiationField
  USE OpacityModule, ONLY: &
    ComputeEmissivity, &
    ComputeAbsorptionOpacity, &
    ComputeScatteringOpacity_ES

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nNodesX_G, nNodesE_G
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: E_N, D_N, T_N, Y_N
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: X_N
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: J_N, H1_N, H2_N, H3_N
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Eta, Chi, Sigma
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dJ_N, dH1_N, dH2_N, dH3_N

  PUBLIC :: CoupleFluidRadiation_ConstantOpacities

CONTAINS


  SUBROUTINE CoupleFluidRadiation_ConstantOpacities &
               ( dt, iX_B0, iX_E0, iX_B1, iX_E1, U_F, dU_F, U_R, dU_R, &
                 EvolveFluid_Option )

    REAL(DP), INTENT(in)  :: &
      dt
    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_B1(3), iX_E0(3), iX_E1(3)
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

  END SUBROUTINE CoupleFluidRadiation_ConstantOpacities


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( X_N(nNodesX_G,3) )
    CALL InitializeNodesX( X_N )

    ALLOCATE( D_N(nNodesX_G) )
    CALL MapForward_FluidField &
           ( uPF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iPF_D), D_N )

    ALLOCATE( T_N(nNodesX_G) )
    CALL MapForward_FluidField &
           ( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_T), T_N )

    ALLOCATE( Y_N(nNodesX_G) )
    CALL MapForward_FluidField &
           ( uAF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),iAF_Ye), Y_N )

    ALLOCATE( J_N(nNodesE_G, nNodesX_G) )
    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_D, 1), &
             J_N(1:nNodesE_G,1:nNodesX_G) )

    ALLOCATE( H1_N(nNodesE_G, nNodesX_G) )
    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I1, 1), &
             H1_N(1:nNodesE_G,1:nNodesX_G) )

    ALLOCATE( H2_N(nNodesE_G, nNodesX_G) )
    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I2, 1), &
             H2_N(1:nNodesE_G,1:nNodesX_G) )

    ALLOCATE( H3_N(nNodesE_G, nNodesX_G) )
    CALL MapForward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I3, 1), &
             H3_N(1:nNodesE_G,1:nNodesX_G) )

    ALLOCATE( Eta  (nNodesE_G, nNodesX_G) )
    ALLOCATE( Chi  (nNodesE_G, nNodesX_G) )
    ALLOCATE( Sigma(nNodesE_G, nNodesX_G) )
    ALLOCATE( dJ_N (nNodesE_G, nNodesX_G) )
    ALLOCATE( dH1_N(nNodesE_G, nNodesX_G) )
    ALLOCATE( dH2_N(nNodesE_G, nNodesX_G) )
    ALLOCATE( dH3_N(nNodesE_G, nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    DEALLOCATE( E_N, X_N )

    DEALLOCATE( D_N, T_N, Y_N )

    CALL MapBackward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_D, 1), &
             J_N(1:nNodesE_G,1:nNodesX_G) )

    CALL MapBackward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I1, 1), &
             H1_N(1:nNodesE_G,1:nNodesX_G) )

    CALL MapBackward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I2, 1), &
             H2_N(1:nNodesE_G,1:nNodesX_G) )

    CALL MapBackward_RadiationField &
           ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I3, 1), &
             H3_N(1:nNodesE_G,1:nNodesX_G) )

    DEALLOCATE( J_N, H1_N, H2_N, H3_N )

    DEALLOCATE( Eta, Chi, Sigma )

    CALL MapBackward_RadiationField &
           ( rhsCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_N, 1), &
             dJ_N(1:nNodesE_G,1:nNodesX_G) )

    CALL MapBackward_RadiationField &
           ( rhsCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_G1, 1), &
             dH1_N(1:nNodesE_G,1:nNodesX_G) )

    CALL MapBackward_RadiationField &
           ( rhsCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_G2, 1), &
             dH2_N(1:nNodesE_G,1:nNodesX_G) )

    CALL MapBackward_RadiationField &
           ( rhsCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_G3, 1), &
             dH3_N(1:nNodesE_G,1:nNodesX_G) )

    DEALLOCATE( dJ_N, dH1_N, dH2_N, dH3_N )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER :: iX, iE

    CALL SetRates

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G

        ! --- Number Density ---

        J_N(iE,iX) &
          = ( J_N(iE,iX) + dt * Eta(iE,iX) ) / ( 1.0_DP + dt * Chi(iE,iX) )

        ! --- Number Flux Density (1) ---

        H1_N(iE,iX) &
          = H1_N(iE,iX) / ( 1.0_DP + dt * ( Chi(iE,iX) + Sigma(iE,iX) ) )

        ! --- Number Flux Density (2) ---

        H2_N(iE,iX) &
          = H2_N(iE,iX) / ( 1.0_DP + dt * ( Chi(iE,iX) + Sigma(iE,iX) ) )

        ! --- Number Flux Density (3) ---

        H3_N(iE,iX) &
          = H3_N(iE,iX) / ( 1.0_DP + dt * ( Chi(iE,iX) + Sigma(iE,iX) ) )

        ! --- Increments ---

        dJ_N (iE,iX) = Eta(iE,iX) - Chi(iE,iX) * J_N(iE,iX)

        dH1_N(iE,iX) = - ( Chi(iE,iX) + Sigma(iE,iX) ) * H1_N(iE,iX)

        dH2_N(iE,iX) = - ( Chi(iE,iX) + Sigma(iE,iX) ) * H2_N(iE,iX)

        dH3_N(iE,iX) = - ( Chi(iE,iX) + Sigma(iE,iX) ) * H3_N(iE,iX)

      END DO
    END DO

  END SUBROUTINE CoupleFluidRadiation


  SUBROUTINE SetRates

    CALL ComputeEmissivity &
           ( E_N, D_N, T_N, Y_N, X_N(:,1), X_N(:,2), X_N(:,3), Eta )

    CALL ComputeAbsorptionOpacity &
           ( E_N, D_N, T_N, Y_N, X_N(:,1), X_N(:,2), X_N(:,3), Chi )

    CALL ComputeScatteringOpacity_ES &
           ( E_N, D_N, T_N, Y_N, X_N(:,1), X_N(:,2), X_N(:,3), Sigma )

  END SUBROUTINE SetRates


END MODULE FluidRadiationCouplingSolutionModule_ConstantOpacities
