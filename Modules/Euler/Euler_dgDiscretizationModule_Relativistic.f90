MODULE Euler_dgDiscretizationModule_Relativistic

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Pi, TwoPi, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, WeightsX_X1, &
    nDOFX_X2, WeightsX_X2, &
    WeightsX_q,  &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, dLXdX2_q, &
    LX_X1_Dn, LX_X1_Up, &
    LX_X2_Dn, LX_X2_Up
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF,                                 &
    iGF_h_1, iGF_h_2, iGF_h_3,                &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_SqrtGm,                               &
    iGF_Alpha,                                &
    iGF_Beta_1, iGF_Beta_2, iGF_Beta_3
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_P, iAF_Gm
  USE Euler_BoundaryConditionsModule_Relativistic, ONLY: &
    Euler_ApplyBoundaryConditions_Relativistic
  USE Euler_UtilitiesModule_Relativistic, ONLY:  &
    Euler_ComputePrimitive_Relativistic,      &
    Euler_Eigenvalues_Relativistic,           &
    Euler_AlphaMiddle_Relativistic,           &
    Euler_Flux_X1_Relativistic,               &
    Euler_Flux_X2_Relativistic,               &
    Euler_Flux_X3_Relativistic,               &
    Euler_StressTensor_Diagonal_Relativistic, &
    Euler_NumericalFlux_HLL_Relativistic,     &
    Euler_NumericalFlux_X1_HLLC_Relativistic, &
    Euler_NumericalFlux_X2_HLLC_Relativistic, &
    Euler_NumericalFlux_X3_HLLC_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromSpecificInternalEnergy, &
    ComputeSoundSpeedFromPrimitive_GR

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: Euler_ComputeIncrement_DG_Explicit_Relativistic

  LOGICAL, PARAMETER :: DisplayTimers = .FALSE.
  REAL(DP) :: Timer_RHS_GR
  REAL(DP) :: Timer_RHS_1_GR, dT_RHS_1_GR
  REAL(DP) :: Timer_RHS_2_GR, dT_RHS_2_GR
  REAL(DP) :: Timer_RHS_3_GR, dT_RHS_3_GR
  REAL(DP) :: Timer_INT_F_GR, dT_INT_F_GR
  REAL(DP) :: Timer_INT_G_GR, dT_INT_G_GR
  REAL(DP) :: Timer_FLX_N_GR, dT_FLX_N_GR
  REAL(DP) :: Timer_Geo

CONTAINS


  SUBROUTINE Euler_ComputeIncrement_DG_Explicit_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU, SuppressBC_Option )

    INTEGER, INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout)        :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(out)          :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    REAL(DP) :: ErrorL1, ErrorIn, Error, X1

    INTEGER  :: iX1, iX2, iX3, iCF, iNodeX, iNodeX1
    REAL(DP) :: dX1, dX2, dX3
    LOGICAL  :: SuppressBC

    dU = Zero

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL Euler_ApplyBoundaryConditions_Relativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL ComputeIncrement_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    ! --- Multiply Inverse Mass Matrix ---

    DO iCF = 1, nCF

      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        dX3 = MeshX(3) % Width(iX3)
        dX2 = MeshX(2) % Width(iX2)
        dX1 = MeshX(1) % Width(iX1)

        dU(:,iX1,iX2,iX3,iCF) &
          = dU(:,iX1,iX2,iX3,iCF) &
            / ( WeightsX_q * G(:,iX1,iX2,iX3,iGF_SqrtGm) * dX1 * dX2 * dX3 )

      END DO
      END DO
      END DO

    END DO

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

  END SUBROUTINE Euler_ComputeIncrement_DG_Explicit_Relativistic


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iAF, iNodeX, iNodeX_X1, iNodeX1
    REAL(DP) :: dX2, dX3, &
                Alpha, Beta_1, Beta_2, Beta_3, &
                AlphaMns, AlphaMdl, AlphaPls, &
                Cs_L(nDOFX_X1), Cs_R(nDOFX_X1), &
                Lambda_L(nCF,nDOFX_X1), Lambda_R(nCF,nDOFX_X1), &
                uCF_L(nDOFX_X1,nCF), uCF_R(nDOFX_X1,nCF), &
                uPF_L(nDOFX_X1,nPF), uPF_R(nDOFX_X1,nPF), &
                uAF_L(nDOFX_X1,nAF), uAF_R(nDOFX_X1,nAF), &
                G_F(nDOFX_X1,nGF), &
                P_P(nDOFX), P_K(nDOFX), &
                Flux_X1_L(nDOFX_X1,nCF), Flux_X1_R(nDOFX_X1,nCF), &
                EigVals_L(nDOFX_X1,nCF), EigVals_R(nDOFX_X1,nCF), &
                EigVals(nDOFX_X1,nCF), &
                NumericalFlux(nDOFX_X1,nCF), &
                uCF_P(nDOFX,nCF), uCF_K(nDOFX,nCF), &    
                uPF_P(nDOFX,nPF), uPF_K(nDOFX,nPF), &
                uAF_P(nDOFX,nAF), uAF_K(nDOFX,nAF), &
                G_P(nDOFX,nGF), G_K(nDOFX,nGF), &
                Flux_X1_q(nDOFX,nCF)

    Timer_RHS_1_GR = 0.0_DP
    Timer_RHS_2_GR = 0.0_DP
    Timer_RHS_3_GR = 0.0_DP
    Timer_INT_F_GR = 0.0_DP
    Timer_INT_G_GR = 0.0_DP
    Timer_FLX_N_GR = 0.0_DP

    CALL Timer_Start( Timer_RHS_GR )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1) + 1

      dX3 = MeshX(3) % Width(iX3)
      dX2 = MeshX(2) % Width(iX2)

      DO iCF = 1, nCF

        uCF_P(:,iCF) = U(:,iX1-1,iX2,iX3,iCF)
        uCF_K(:,iCF) = U(:,iX1,  iX2,iX3,iCF)

      END DO

      DO iGF = 1, nGF

        G_P(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
        G_K(:,iGF) = G(:,iX1,  iX2,iX3,iGF)

      END DO

      !--------------------
      ! --- Volume Term ---
      !--------------------

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        CALL Euler_ComputePrimitive_Relativistic &
               ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 G_K(:,iGF_Gm_dd_11),                               &
                 G_K(:,iGF_Gm_dd_22),                               &
                 G_K(:,iGF_Gm_dd_33),                               &
                 P_K )

        DO iNodeX = 1, nDOFX

          Flux_X1_q(iNodeX,:) &
            = Euler_Flux_X1_Relativistic    &
                ( uPF_K(iNodeX,iPF_D ),     &
                  uPF_K(iNodeX,iPF_V1),     &
                  uPF_K(iNodeX,iPF_V2),     &
                  uPF_K(iNodeX,iPF_V3),     &
                  uPF_K(iNodeX,iPF_E ),     &
                  uPF_K(iNodeX,iPF_Ne),     &
                  P_K(iNodeX),              &
                  G_K(iNodeX,iGF_Gm_dd_11), &
                  G_K(iNodeX,iGF_Gm_dd_22), &
                  G_K(iNodeX,iGF_Gm_dd_33), &
                  G_K(iNodeX,iGF_Alpha),    &
                  G_K(iNodeX,iGF_Beta_1) )

        END DO

        CALL Timer_Start( dT_RHS_1_GR )

        DO iCF = 1, nCF

          Flux_X1_q(:,iCF) &
            = dX2 * dX3 * WeightsX_q * G_K(:,iGF_Alpha) &
                * G_K(:,iGF_SqrtGm) * Flux_X1_q(:,iCF)

          CALL DGEMV &
                 ( 'T', nDOFX, nDOFX, One, dLXdX1_q, nDOFX, &
                   Flux_X1_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

        END DO

        CALL Timer_Stop( dT_RHS_1_GR )

        CALL Timer_Add( Timer_RHS_1_GR, dT_RHS_1_GR )

      END IF

      !------------------------
      ! --- Divergence Term ---
      !------------------------

      ! --- Interpolate Fluid Fields ---

      CALL Timer_Start( dT_INT_F_GR )

      DO iCF = 1, nCF

        ! --- Left States ---

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
                 uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )

        ! --- Right States ---

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                 uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )

      END DO

      ! --- Left State Pressure ---

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               P_P, 1, Zero, uAF_L(:,iAF_P), 1 )

      ! --- Right State Pressure ---

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               P_K, 1, Zero, uAF_R(:,iAF_P), 1 )

      CALL Timer_Stop( dT_INT_F_GR )

      CALL Timer_Add( Timer_INT_F_GR, dT_INT_F_GR )

      ! --- Interpolate Geometry Fields ---

      CALL Timer_Start( dT_INT_G_GR )

      ! --- Face States (Average of Left and Right States) ---

      G_F = Zero

      ! --- Scale Factors ---

      DO iGF = iGF_h_1, iGF_h_3

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )

        G_F(:,iGF) = MAX( G_F(:,iGF), SqrtTiny )

      END DO

      CALL ComputeGeometryX_FromScaleFactors( G_F )

      ! --- Lapse Function ---

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_P(:,iGF_Alpha), 1, Zero, G_F(:,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_K(:,iGF_Alpha), 1, Half, G_F(:,iGF_Alpha), 1 )

      G_F(:,iGF_Alpha) = MAX( G_F(:,iGF_Alpha), SqrtTiny )

      CALL Timer_Stop( dT_INT_G_GR )

      CALL Timer_Add( Timer_INT_G_GR, dT_INT_G_GR )

      ! --- Left State Primitive, etc. ---

      CALL Euler_ComputePrimitive_Relativistic &
             ( uCF_L(:,iCF_D ), uCF_L(:,iCF_S1), uCF_L(:,iCF_S2), &
               uCF_L(:,iCF_S3), uCF_L(:,iCF_E ), uCF_L(:,iCF_Ne), &
               uPF_L(:,iPF_D ), uPF_L(:,iPF_V1), uPF_L(:,iPF_V2), &
               uPF_L(:,iPF_V3), uPF_L(:,iPF_E ), uPF_L(:,iPF_Ne), &
               G_F(:,iGF_Gm_dd_11),                               &
               G_F(:,iGF_Gm_dd_22),                               &
               G_F(:,iGF_Gm_dd_33),                               &
               uAF_L(:,iAF_P) )

      CALL ComputeSoundSpeedFromPrimitive_GR &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L )

      DO iNodeX_X1 = 1, nDOFX_X1

        Lambda_L(:,iNodeX_X1) &
          = Euler_Eigenvalues_Relativistic     &
              ( uPF_L(iNodeX_X1,iPF_V1),       &
                Cs_L (iNodeX_X1),              &
                uPF_L(iNodeX_X1,iPF_V1),       &
                uPF_L(iNodeX_X1,iPF_V2),       &
                uPF_L(iNodeX_X1,iPF_V3),       &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                G_F  (iNodeX_X1,iGF_Alpha),    &
                G_F  (iNodeX_X1,iGF_Beta_1) )

        Flux_X1_L(iNodeX_X1,:) &
          = Euler_Flux_X1_Relativistic         &
              ( uPF_L(iNodeX_X1,iPF_D ),       &
                uPF_L(iNodeX_X1,iPF_V1),       &
                uPF_L(iNodeX_X1,iPF_V2),       &
                uPF_L(iNodeX_X1,iPF_V3),       &
                uPF_L(iNodeX_X1,iPF_E ),       &
                uPF_L(iNodeX_X1,iPF_Ne),       &
                uAF_L(iNodeX_X1,iAF_P ),       &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                G_F  (iNodeX_X1,iGF_Alpha),    &
                G_F  (iNodeX_X1,iGF_Beta_1) )

      END DO

      ! --- Right State Primitive, etc. ---

      CALL Euler_ComputePrimitive_Relativistic &
             ( uCF_R(:,iCF_D ), uCF_R(:,iCF_S1), uCF_R(:,iCF_S2), &
               uCF_R(:,iCF_S3), uCF_R(:,iCF_E ), uCF_R(:,iCF_Ne), &
               uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), uPF_R(:,iPF_V2), &
               uPF_R(:,iPF_V3), uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
               G_F(:,iGF_Gm_dd_11),                               &
               G_F(:,iGF_Gm_dd_22),                               &
               G_F(:,iGF_Gm_dd_33),                               &
               uAF_R(:,iAF_P) )

      CALL ComputeSoundSpeedFromPrimitive_GR &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R )

      DO iNodeX_X1 = 1, nDOFX_X1

        Lambda_R(:,iNodeX_X1) &
          = Euler_Eigenvalues_Relativistic     &
              ( uPF_R(iNodeX_X1,iPF_V1),       &
                Cs_R (iNodeX_X1),              &
                uPF_R(iNodeX_X1,iPF_V1),       &
                uPF_R(iNodeX_X1,iPF_V2),       &
                uPF_R(iNodeX_X1,iPF_V3),       &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                G_F  (iNodeX_X1,iGF_Alpha),    &
                G_F  (iNodeX_X1,iGF_Beta_1) )

        Flux_X1_R(iNodeX_X1,:) &
          = Euler_Flux_X1_Relativistic         &
              ( uPF_R(iNodeX_X1,iPF_D ),       &
                uPF_R(iNodeX_X1,iPF_V1),       &
                uPF_R(iNodeX_X1,iPF_V2),       &
                uPF_R(iNodeX_X1,iPF_V3),       &
                uPF_R(iNodeX_X1,iPF_E ),       &
                uPF_R(iNodeX_X1,iPF_Ne),       &
                uAF_R(iNodeX_X1,iAF_P ),       &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                G_F  (iNodeX_X1,iGF_Alpha),    &
                G_F  (iNodeX_X1,iGF_Beta_1) )

      END DO

      ! --- Numerical Flux ---

      CALL Timer_Start( dT_FLX_N_GR )

      DO iNodeX_X1 = 1, nDOFX_X1

        AlphaMns &
          = MAX( Zero, &
                 MAXVAL( - Lambda_L(:,iNodeX_X1) ), &
                 MAXVAL( - Lambda_R(:,iNodeX_X1) ) )

        AlphaPls &
          = MAX( Zero, &
                 MAXVAL( + Lambda_L(:,iNodeX_X1) ), &
                 MAXVAL( + Lambda_R(:,iNodeX_X1) ) )

        AlphaMdl &
          = Euler_AlphaMiddle_Relativistic         &
              ( uCF_L    (iNodeX_X1,iCF_D ),       &
                uCF_L    (iNodeX_X1,iCF_S1),       &
                uCF_L    (iNodeX_X1,iCF_E ),       &
                Flux_X1_L(iNodeX_X1,iCF_D ),       &
                Flux_X1_L(iNodeX_X1,iCF_S1),       &
                Flux_X1_L(iNodeX_X1,iCF_E ),       &
                uCF_R    (iNodeX_X1,iCF_D ),       &
                uCF_R    (iNodeX_X1,iCF_S1),       &
                uCF_R    (iNodeX_X1,iCF_E ),       &
                Flux_X1_R(iNodeX_X1,iCF_D ),       &
                Flux_X1_R(iNodeX_X1,iCF_S1),       &
                Flux_X1_R(iNodeX_X1,iCF_E ),       &
                AlphaPls, AlphaMns,                &
                G_F      (iNodeX_X1,iGF_Gm_dd_11), &
                G_F      (iNodeX_X1,iGF_Alpha),    &
                G_F      (iNodeX_X1,iGF_Beta_1) )

        NumericalFlux(iNodeX_X1,:) &
!          = Euler_NumericalFlux_HLLC_X2_Relativistic &
          = Euler_NumericalFlux_HLL_Relativistic     &
              ( uCF_L    (iNodeX_X1,:),              &
                uCF_R    (iNodeX_X1,:),              &
                Flux_X1_L(iNodeX_X1,:),              &
                Flux_X1_R(iNodeX_X1,:),              &
                AlphaPls, AlphaMns, AlphaMdl,        &
                G_F      (iNodeX_X1,iGF_Gm_dd_11),   &
                uPF_L    (iNodeX_X1,iPF_V1),         &
                uPF_R    (iNodeX_X1,iPF_V1),         &
                uAF_L    (iNodeX_X1,iAF_P),          &
                uAF_R    (iNodeX_X1,iAF_P),          &
                G_F      (iNodeX_X1,iGF_Alpha),      &
                G_F      (iNodeX_X1,iGF_Beta_1) )

      END DO

      DO iCF = 1, nCF

        NumericalFlux(:,iCF) &
          = dX2 * dX3 * WeightsX_X1 * G_F(:,iGF_Alpha) &
              * G_F(:,iGF_SqrtGm) * NumericalFlux(:,iCF)

      END DO

      CALL Timer_Stop( dT_FLX_N_GR )

      CALL Timer_Add( Timer_FLX_N_GR, dT_FLX_N_GR )

      ! --- Contribution to This Element ---

      CALL Timer_Start( dT_RHS_2_GR )

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        DO iCF = 1, nCF

          CALL DGEMV &
                 ( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Dn, nDOFX_X1, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

        END DO

      END IF

      CALL Timer_Stop( dT_RHS_2_GR )

      CALL Timer_Add( Timer_RHS_2_GR, dT_RHS_2_GR )

      ! --- Contribution to Previous Element ---

      CALL Timer_Start( dT_RHS_3_GR )

      IF( iX1 .GT. iX_B0(1) )THEN

        DO iCF = 1, nCF

          CALL DGEMV &
               ( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Up, nDOFX_X1, &
                 NumericalFlux(:,iCF), 1, One, dU(:,iX1-1,iX2,iX3,iCF), 1 )

        END DO

      END IF

      CALL Timer_Stop( dT_RHS_3_GR )

      CALL Timer_Add( Timer_RHS_3_GR, dT_RHS_3_GR )

    END DO
    END DO
    END DO

    CALL Timer_Stop( Timer_RHS_GR )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A24,ES10.4E2)') &
        '', 'ComputeRHS_Euler: ', Timer_RHS_GR
      WRITE(*,*)
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 1: ', Timer_RHS_1_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 2: ', Timer_RHS_2_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 3: ', Timer_RHS_3_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'INT F: ', Timer_INT_F_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'INT G: ', Timer_INT_G_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'FLX N: ', Timer_FLX_N_GR
      WRITE(*,*)
      WRITE(*,'(A6,A18,ES10.4E2)') &
           '', 'Sum: ', Timer_RHS_1_GR+Timer_RHS_2_GR+Timer_RHS_3_GR+&
           Timer_INT_F_GR &
        + Timer_INT_G_GR + Timer_FLX_N_GR

    END IF

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_Divergence_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iAF, iNodeX, iNodeX_X2, iNodeX2
    REAL(DP) :: dX1, dX3, &
                Alpha, Beta_1, Beta_2, Beta_3, &
                AlphaMns, AlphaMdl, AlphaPls, &
                Cs_L(nDOFX_X2), Cs_R(nDOFX_X2), &
                Lambda_L(nCF,nDOFX_X2), Lambda_R(nCF,nDOFX_X2), &
                uCF_L(nDOFX_X2,nCF), uCF_R(nDOFX_X2,nCF), &
                uPF_L(nDOFX_X2,nPF), uPF_R(nDOFX_X2,nPF), &
                uAF_L(nDOFX_X2,nAF), uAF_R(nDOFX_X2,nAF), &
                G_F(nDOFX_X2,nGF), &
                P_P(nDOFX) , P_K(nDOFX), & 
                Flux_X2_L(nDOFX_X2,nCF), Flux_X2_R(nDOFX_X2,nCF), &
                EigVals_L(nDOFX_X2,nCF), EigVals_R(nDOFX_X2,nCF), &
                EigVals(nDOFX_X2,nCF), &
                NumericalFlux(nDOFX_X2,nCF), &
                uCF_P(nDOFX,nCF), uCF_K(nDOFX,nCF), &
                uPF_P(nDOFX,nPF), uPF_K(nDOFX,nPF), &
                uAF_P(nDOFX,nAF), uAF_K(nDOFX,nAF), &
                G_P(nDOFX,nGF), G_K(nDOFX,nGF), &
                Flux_X2_q(nDOFX,nCF)

    IF( iX_E0(2) .EQ. iX_B0(2) ) RETURN

    Timer_RHS_1_GR = 0.0_DP
    Timer_RHS_2_GR = 0.0_DP
    Timer_RHS_3_GR = 0.0_DP
    Timer_INT_F_GR = 0.0_DP
    Timer_INT_G_GR = 0.0_DP
    Timer_FLX_N_GR = 0.0_DP

    CALL Timer_Start( Timer_RHS_GR )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2) + 1
    DO iX1 = iX_B0(1), iX_E0(1)

      dX3 = MeshX(3) % Width(iX3)
      dX1 = MeshX(1) % Width(iX1)

      DO iCF = 1, nCF

        uCF_P(:,iCF) = U(:,iX1,iX2-1,iX3,iCF)
        uCF_K(:,iCF) = U(:,iX1,iX2,  iX3,iCF)

      END DO

      DO iGF = 1, nGF

        G_P(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
        G_K(:,iGF) = G(:,iX1,iX2,  iX3,iGF)

      END DO

      !--------------------
      ! --- Volume Term ---
      !--------------------

      IF( iX2 .LT. iX_E0(2) + 1 )THEN

        CALL Euler_ComputePrimitive_Relativistic &
               ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 G_K(:,iGF_Gm_dd_11),                               &
                 G_K(:,iGF_Gm_dd_22),                               &
                 G_K(:,iGF_Gm_dd_33),                               &
                 P_K )

        DO iNodeX = 1, nDOFX

          Flux_X2_q(iNodeX,:) &
            = Euler_Flux_X2_Relativistic      &
                ( uPF_K(iNodeX,iPF_D ),       &
                  uPF_K(iNodeX,iPF_V1),       &
                  uPF_K(iNodeX,iPF_V2),       &
                  uPF_K(iNodeX,iPF_V3),       &
                  uPF_K(iNodeX,iPF_E ),       &
                  uPF_K(iNodeX,iPF_Ne),       &
                  P_K  (iNodeX),              &
                  G_K  (iNodeX,iGF_Gm_dd_11), &
                  G_K  (iNodeX,iGF_Gm_dd_22), &
                  G_K  (iNodeX,iGF_Gm_dd_33), &
                  G_K  (iNodeX,iGF_Alpha),    &
                  G_K  (iNodeX,iGF_Beta_2) )

        END DO

        CALL Timer_Start( dT_RHS_2_GR )

        DO iCF = 1, nCF

          Flux_X2_q(:,iCF) &
            = dX1 * dX3 * WeightsX_q * G_K(:,iGF_Alpha) &
                * G_K(:,iGF_SqrtGm) * Flux_X2_q(:,iCF)

          CALL DGEMV &
                 ( 'T', nDOFX, nDOFX, One, dLXdX2_q, nDOFX, &
                   Flux_X2_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

        END DO

        CALL Timer_Stop( dT_RHS_2_GR )

        CALL Timer_Add( Timer_RHS_2_GR, dT_RHS_2_GR )

      END IF ! --- End of volume term

      !------------------------
      ! --- Divergence Term ---
      !------------------------

      ! --- Interpolate Fluid Fields ---

      CALL Timer_Start( dT_INT_F_GR )

      DO iCF = 1, nCF

        ! --- Left States ---

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
                 uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )

        ! --- Right States ---

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
                 uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )

      END DO

      ! --- Left State Pressure ---

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               P_P, 1, Zero, uAF_L(:,iAF_P), 1 )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               P_K, 1, Zero, uAF_R(:,iAF_P), 1 )

      CALL Timer_Stop( dT_INT_F_GR )

      CALL Timer_Add( Timer_INT_F_GR, dT_INT_F_GR )

      ! --- Interpolate Geometry Fields ---

      CALL Timer_Start( dT_INT_G_GR )

      ! --- Face States (Average of Left and Right States ) ---

      G_F = Zero

      ! --- Scale Factors ---

      DO iGF = iGF_h_1, iGF_h_3

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
                 G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                 G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )

        G_F(1:nDOFX_X2,iGF) = MAX( G_F(:,iGF), SqrtTiny )

      END DO

      CALL ComputeGeometryX_FromScaleFactors( G_F )

      ! --- Lapse Function ---

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               G_P(:,iGF_Alpha), 1, Zero, G_F(:,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               G_K(:,iGF_Alpha), 1, Half, G_F(:,iGF_Alpha), 1 )

      G_F(:,iGF_Alpha) = MAX( G_F(:,iGF_Alpha), SqrtTiny )

      CALL Timer_Stop( dT_INT_G_GR )

      CALL Timer_Add( Timer_INT_G_GR, dT_INT_G_GR )

      ! --- Left State Primitive, etc. ---

      CALL Euler_ComputePrimitive_Relativistic &
             ( uCF_L(:,iCF_D ), uCF_L(:,iCF_S1), uCF_L(:,iCF_S2), &
               uCF_L(:,iCF_S3), uCF_L(:,iCF_E ), uCF_L(:,iCF_Ne), &
               uPF_L(:,iPF_D ), uPF_L(:,iPF_V1), uPF_L(:,iPF_V2), &
               uPF_L(:,iPF_V3), uPF_L(:,iPF_E ), uPF_L(:,iPF_Ne), &
               G_F(:,iGF_Gm_dd_11),                               &
               G_F(:,iGF_Gm_dd_22),                               &
               G_F(:,iGF_Gm_dd_33),                               &
               uAF_L(:,iAF_P) )

      CALL ComputeSoundSpeedFromPrimitive_GR &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L )

      DO iNodeX_X2 = 1, nDOFX_X2

        Lambda_L(:,iNodeX_X2) &
          = Euler_Eigenvalues_Relativistic                  &
              ( uPF_L(iNodeX_X2,iPF_V2),       &
                Cs_L (iNodeX_X2),              &
                uPF_L(iNodeX_X2,iPF_V1),       &
                uPF_L(iNodeX_X2,iPF_V2),       &
                uPF_L(iNodeX_X2,iPF_V3),       &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_33), &
                G_F  (iNodeX_X2,iGF_Alpha),    &
                G_F  (iNodeX_X2,iGF_Beta_2) )

        Flux_X2_L(iNodeX_X2,:) &
          = Euler_Flux_X2_Relativistic                      &
              ( uPF_L(iNodeX_X2,iPF_D ),       &
                uPF_L(iNodeX_X2,iPF_V1),       &
                uPF_L(iNodeX_X2,iPF_V2),       &
                uPF_L(iNodeX_X2,iPF_V3),       &
                uPF_L(iNodeX_X2,iPF_E ),       &
                uPF_L(iNodeX_X2,iPF_Ne),       &
                uAF_L(iNodeX_X2,iAF_P ),       &
                G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_33), &
                G_F  (iNodeX_X2,iGF_Alpha),    &
                G_F  (iNodeX_X2,iGF_Beta_2) )

      END DO

      ! --- Right State Primitive, etc. ---

      CALL Euler_ComputePrimitive_Relativistic &
             ( uCF_R(:,iCF_D ), uCF_R(:,iCF_S1), uCF_R(:,iCF_S2), &
               uCF_R(:,iCF_S3), uCF_R(:,iCF_E ), uCF_R(:,iCF_Ne), &
               uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), uPF_R(:,iPF_V2), &
               uPF_R(:,iPF_V3), uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
               G_F(:,iGF_Gm_dd_11),                               &
               G_F(:,iGF_Gm_dd_22),                               &
               G_F(:,iGF_Gm_dd_33),                               &
               uAF_R(:,iAF_P) )

      CALL ComputeSoundSpeedFromPrimitive_GR &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R )

      DO iNodeX_X2 = 1, nDOFX_X2

        Lambda_R(:,iNodeX_X2) &
          = Euler_Eigenvalues_Relativistic     &
              ( uPF_R(iNodeX_X2,iPF_V2),       &
                Cs_R (iNodeX_X2),              &
                uPF_R(iNodeX_X2,iPF_V1),       &
                uPF_R(iNodeX_X2,iPF_V2),       &
                uPF_R(iNodeX_X2,iPF_V3),       &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_33), &
                G_F  (iNodeX_X2,iGF_Alpha),    &
                G_F  (iNodeX_X2,iGF_Beta_2) )

        Flux_X2_R(iNodeX_X2,:) &
          = Euler_Flux_X2_Relativistic         &
              ( uPF_R(iNodeX_X2,iPF_D ),       &
                uPF_R(iNodeX_X2,iPF_V1),       &
                uPF_R(iNodeX_X2,iPF_V2),       &
                uPF_R(iNodeX_X2,iPF_V3),       &
                uPF_R(iNodeX_X2,iPF_E ),       &
                uPF_R(iNodeX_X2,iPF_Ne),       &
                uAF_R(iNodeX_X2,iAF_P ),       &
                G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_33), &
                G_F  (iNodeX_X2,iGF_Alpha),    &
                G_F  (iNodeX_X2,iGF_Beta_2) )

      END DO

      ! --- Numerical Flux ---

      CALL Timer_Start( dT_FLX_N_GR )

      DO iNodeX_X2 = 1, nDOFX_X2

        AlphaMns &
          = MAX( Zero, &
                 MAXVAL( - Lambda_L(:,iNodeX_X2) ), &
                 MAXVAL( - Lambda_R(:,iNodeX_X2) ) )

        AlphaPls &
          = MAX( Zero, &
                 MAXVAL( + Lambda_L(:,iNodeX_X2) ), &
                 MAXVAL( + Lambda_R(:,iNodeX_X2) ) )

        AlphaMdl &
          = Euler_AlphaMiddle_Relativistic         &
              ( uCF_L    (iNodeX_X2,iCF_D ),       &
                uCF_L    (iNodeX_X2,iCF_S2),       &
                uCF_L    (iNodeX_X2,iCF_E ),       &
                Flux_X2_L(iNodeX_X2,iCF_D ),       &
                Flux_X2_L(iNodeX_X2,iCF_S2),       &
                Flux_X2_L(iNodeX_X2,iCF_E ),       &
                uCF_R    (iNodeX_X2,iCF_D ),       &
                uCF_R    (iNodeX_X2,iCF_S2),       &
                uCF_R    (iNodeX_X2,iCF_E ),       &
                Flux_X2_R(iNodeX_X2,iCF_D ),       &
                Flux_X2_R(iNodeX_X2,iCF_S2),       &
                Flux_X2_R(iNodeX_X2,iCF_E ),       &
                AlphaPls, AlphaMns,                &
                G_F      (iNodeX_X2,iGF_Gm_dd_22), &
                G_F      (iNodeX_X2,iGF_Alpha),    &
                G_F      (iNodeX_X2,iGF_Beta_2) )

        NumericalFlux(iNodeX_X2,:) &
!          = Euler_NumericalFlux_X2_HLLC_Relativistic &
          = Euler_NumericalFlux_HLL_Relativistic     &
              ( uCF_L    (iNodeX_X2,1:nCF),          &
                uCF_R    (iNodeX_X2,1:nCF),          &
                Flux_X2_L(iNodeX_X2,1:nCF),          &
                Flux_X2_R(iNodeX_X2,1:nCF),          &
                AlphaPls, AlphaMns, AlphaMdl,        &
                G_F      (iNodeX_X2,iGF_Gm_dd_22),   &
                uPF_L    (iNodeX_X2,iPF_V2),         &
                uPF_R    (iNodeX_X2,iPF_V2),         &
                uAF_L    (iNodeX_X2,iAF_P),          &
                uAF_R    (iNodeX_X2,iAF_P),          &
                G_F      (iNodeX_X2,iGF_Alpha),      &
                G_F      (iNodeX_X2,iGF_Beta_2) )


      END DO

      DO iCF = 1, nCF

        NumericalFlux(:,iCF) &
          = dX1 * dX3 * WeightsX_X2 * G_F(:,iGF_Alpha) &
              * G_F(:,iGF_SqrtGm) * NumericalFlux(:,iCF)

      END DO

      CALL Timer_Stop( dT_FLX_N_GR )

      CALL Timer_Add( Timer_FLX_N_GR, dT_FLX_N_GR )

      ! --- Contribution to This Element ---

      CALL Timer_Start( dT_RHS_2_GR )

      IF( iX2 .LT. iX_E0(2) + 1 )THEN

        DO iCF = 1, nCF

          CALL DGEMV &
                 ( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Dn, nDOFX_X2, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

        END DO

      END IF

      CALL Timer_Stop( dT_RHS_2_GR )

      CALL Timer_Add( Timer_RHS_2_GR, dT_RHS_2_GR )

      ! --- Contribution to Previous Element ---

      CALL Timer_Start( dT_RHS_3_GR )

      IF( iX2 .GT. iX_B0(2) )THEN

        DO iCF = 1, nCF

          CALL DGEMV &
                 ( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Up, nDOFX_X2, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2-1,iX3,iCF), 1 )

        END DO

      END IF

      CALL Timer_Stop( dT_RHS_3_GR )

      CALL Timer_Add( Timer_RHS_3_GR, dT_RHS_3_GR )

    END DO
    END DO
    END DO

    CALL Timer_Stop( Timer_RHS_GR )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A24,ES10.4E2)') &
        '', 'ComputeRHS_Euler: ', Timer_RHS_GR
      WRITE(*,*)
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 1: ', Timer_RHS_1_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 2: ', Timer_RHS_2_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'RHS 3: ', Timer_RHS_3_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'INT F: ', Timer_INT_F_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'INT G: ', Timer_INT_G_GR
      WRITE(*,'(A6,A18,ES10.4E2)') &
        '', 'FLX N: ', Timer_FLX_N_GR
      WRITE(*,*)
      WRITE(*,'(A6,A18,ES10.4E2)') &
           '', 'Sum: ', Timer_RHS_1_GR+Timer_RHS_2_GR+Timer_RHS_3_GR+&
           Timer_INT_F_GR &
        + Timer_INT_G_GR + Timer_FLX_N_GR

    END IF

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    IF( iX_E0(3) .EQ. iX_B0(3) )THEN
      RETURN
    ELSE
      STOP 'ComputeIncrement_Divergence_X3 not yet implemented'
    END IF

  END SUBROUTINE ComputeIncrement_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX
    REAL(DP) :: dX1, dX2, dX3, &
                P_K(nDOFX), &
                dh1dX1(nDOFX), dh2dX1(nDOFX), dh3dX1(nDOFX), &
                dadx1(nDOFX), &
                Stress(nDOFX,3), &
                uCF_K(nDOFX,nCF), &
                uPF_K(nDOFX,nPF), &
                G_K(nDOFX,nGF), &
                G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF), &
                G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF)

    CALL Timer_Start( Timer_Geo )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX3 = MeshX(3) % Width(iX3)
      dX2 = MeshX(2) % Width(iX2)
      dX1 = MeshX(1) % Width(iX1)

      DO iCF = 1, nCF

        uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

      END DO

      DO iGF = 1, nGF

        G_P_X1(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
        G_K   (:,iGF) = G(:,iX1,  iX2,iX3,iGF)
        G_N_X1(:,iGF) = G(:,iX1+1,iX2,iX3,iGF)

      END DO

      CALL Euler_ComputePrimitive_Relativistic &
           ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
             uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
             uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
             uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
             G_K(:,iGF_Gm_dd_11),                               &
             G_K(:,iGF_Gm_dd_22),                               &
             G_K(:,iGF_Gm_dd_33),                               &
             P_K )

      DO iNodeX = 1, nDOFX

        Stress(iNodeX,:) &
          = Euler_StressTensor_Diagonal_Relativistic        &
              ( uCF_K(iNodeX,iCF_S1), uCF_K(iNodeX,iCF_S2), &
                uCF_K(iNodeX,iCF_S3), uPF_K(iNodeX,iPF_V1), &
                uPF_K(iNodeX,iPF_V2), uPF_K(iNodeX,iPF_V3), &
                P_K  (iNodeX) )

      END DO

      ! --- Scale Factor Derivatives wrt X1 ---

      ! --- Face States (Average of Left and Right States) ---

      DO iGF = iGF_h_1, iGF_h_3

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P_X1(:,iGF), 1, Zero, G_X1_Dn(:,iGF), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K   (:,iGF), 1, Half, G_X1_Dn(:,iGF), 1 )

        G_X1_Dn(:,iGF) = MAX( G_X1_Dn(:,iGF), SqrtTiny )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_K   (:,iGF), 1, Zero, G_X1_Up(:,iGF), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_N_X1(:,iGF), 1, Half, G_X1_Up(:,iGF), 1 )

        G_X1_Up(:,iGF) = MAX( G_X1_Up(:,iGF), SqrtTiny )

      END DO

      ! --- dh1dx1 ---

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_h_1), 1, Zero, dh1dX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_1), 1,  One, dh1dX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_1), 1,  One, dh1dX1, 1 )

      dh1dx1 = dh1dx1 / ( WeightsX_q * dX1 )

      ! --- dh2dx1 ---

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_2), 1,  One, dh2dX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_2), 1,  One, dh2dX1, 1 )

      dh2dx1 = dh2dx1 / ( WeightsX_q(:) * dX1 )

      ! --- dh3dx1 ---

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_3), 1,  One, dh3dX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_3), 1,  One, dh3dX1, 1 )

      dh3dx1 = dh3dx1 / ( WeightsX_q * dX1 )

      ! --- Lapse Function Derivative wrt X1 ---

      ! --- Face States (Average of Left and Right States) ---

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_P_X1(:,iGF_Alpha), 1, Zero, G_X1_Dn(:,iGF_Alpha), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_K   (:,iGF_Alpha), 1, Half, G_X1_Dn(:,iGF_Alpha), 1 )

      G_X1_Dn(:,iGF_Alpha) = MAX( G_X1_Dn(:,iGF_Alpha), SqrtTiny )

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_K   (:,iGF_Alpha), 1, Zero, G_X1_Up(:,iGF_Alpha), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_N_X1(:,iGF_Alpha), 1, Half, G_X1_Up(:,iGF_Alpha), 1 )

      G_X1_Up(:,iGF_Alpha) = MAX( G_X1_Up(:,iGF_Alpha), SqrtTiny )

      ! --- dadx1 ---

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_Alpha), 1, Zero,  &
                  dadX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_Alpha), 1,  One,  &
                  dadX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_Alpha), 1,  One,  &
                  dadX1, 1 )

      dadx1 = dadx1 / ( WeightsX_q * dX1 )

      ! --- Compute Increments ---

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1)                                &
            + G_K(:,iGF_Alpha)                                    &
                * ( ( Stress(:,1) * dh1dX1 ) / G_K(:,iGF_h_1)     &
                    + ( Stress(:,2) * dh2dX1 ) / G_K(:,iGF_h_2)   &
                    + ( Stress(:,3) * dh3dX1 ) / G_K(:,iGF_h_3) ) &
            - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx1

      dU(:,iX1,iX2,iX3,iCF_E) &
        = dU(:,iX1,iX2,iX3,iCF_E) &
            - ( uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11) ) * dadx1

    END DO
    END DO
    END DO

    CALL Timer_Stop( Timer_Geo )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES10.4E2)') &
        '', 'ComputeIncrement_Geometry: ', Timer_Geo

    END IF

  END SUBROUTINE ComputeIncrement_Geometry


  SUBROUTINE Timer_Start( Timer )

    REAL(DP) :: Timer

    Timer = MPI_WTIME( )

  END SUBROUTINE Timer_Start


  SUBROUTINE Timer_Stop( Timer )

    REAL(DP) :: Timer

    Timer = MPI_WTIME( ) - Timer

  END SUBROUTINE Timer_Stop


  SUBROUTINE Timer_Add( Timer, dT )

    REAL(DP) :: Timer, dT

    Timer = Timer + dT

  END SUBROUTINE Timer_Add


END MODULE Euler_dgDiscretizationModule_Relativistic
