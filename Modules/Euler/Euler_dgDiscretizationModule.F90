MODULE Euler_dgDiscretizationModule

  USE KindModule, ONLY: &
    DP,       &
    Zero,     &
    SqrtTiny, &
    Third,    &
    Half,     &
    One,      &
    Three,    &
    Pi,       &
    TwoPi
  USE TimersModule_Euler,  ONLY: &
    TimersStart_Euler,            &
    TimersStop_Euler,             &
    Timer_Euler_dgDiscretization, &
    Timer_Euler_Divergence_X1,    &
    Timer_Euler_Divergence_X2,    &
    Timer_Euler_Divergence_X3,    &
    Timer_Euler_Geometry,         &
    Timer_Euler_Gravity,          &
    Timer_Euler_SurfaceTerm,      &
    Timer_Euler_NumericalFlux,    &
    Timer_Euler_VolumeTerm,       &
    Timer_Euler_Increment,        &
    Timer_Euler_MatrixVectorMultiply
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX, &
    nNodesX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1,    &
    nDOFX_X2,    &
    nDOFX_X3,    &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3, &
    WeightsX_q,  &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    dLXdX2_q, &
    dLXdX3_q, &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE GeometryFieldsModule, ONLY: &
    nGF,          &
    iGF_h_1,      &
    iGF_h_2,      &
    iGF_h_3,      &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm,   &
    iGF_Alpha,    &
    iGF_Psi,      &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_Phi_N,    &
    CoordinateSystem
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE FluidFieldsModule, ONLY: &
    nCF,       &
    iCF_D,     &
    iCF_S1,    &
    iCF_S2,    &
    iCF_S3,    &
    iCF_E,     &
    iCF_Ne,    &
    nPF,       &
    iPF_D,     &
    iPF_V1,    &
    iPF_V2,    &
    iPF_V3,    &
    iPF_E,     &
    iPF_Ne,    &
    nAF,       &
    iAF_P,     &
    iDF_Sh_X1, &
    iDF_Sh_X2, &
    iDF_Sh_X3
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler,      &
    Eigenvalues_Euler,           &
    AlphaMiddle_Euler,           &
    Flux_X1_Euler,               &
    Flux_X2_Euler,               &
    Flux_X3_Euler,               &
    StressTensor_Diagonal_Euler, &
    NumericalFlux_Euler_X1,      &
    NumericalFlux_Euler_X2,      &
    NumericalFlux_Euler_X3
  USE Euler_DiscontinuityDetectionModule, ONLY: &
    DetectShocks_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive
  USE InputOutputModuleHDF, ONLY: &
    WriteSourceTermDiagnosticsHDF

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_Euler_DG_Explicit

  LOGICAL,  PUBLIC :: WriteSourceTerms
  REAL(DP), PUBLIC :: Time


CONTAINS


  SUBROUTINE ComputeIncrement_Euler_DG_Explicit &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, SuppressBC_Option )

    INTEGER,  INTENT(in)            :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)            :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout)         :: &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout)         :: &
      D (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(out)           :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
    LOGICAL,  INTENT(in),  OPTIONAL :: &
      SuppressBC_Option

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: dX1, dX2, dX3
    LOGICAL  :: SuppressBC

    CALL TimersStart_Euler( Timer_Euler_dgDiscretization )

    dU = Zero

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC )THEN

      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

      CALL DetectShocks_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    END IF

    CALL ComputeIncrement_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    CALL ComputeIncrement_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    CALL ComputeIncrement_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    ! --- Multiply Inverse Mass Matrix ---

    CALL TimersStart_Euler( Timer_Euler_Increment )

    DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)

        dX1 = MeshX(1) % Width(iX1)
        dX2 = MeshX(2) % Width(iX2)
        dX3 = MeshX(3) % Width(iX3)

        dU(:,iX1,iX2,iX3,iCF) &
          = dU(:,iX1,iX2,iX3,iCF) &
              / ( WeightsX_q * G(:,iX1,iX2,iX3,iGF_SqrtGm) * dX1 * dX2 * dX3 )

      END DO
      END DO
      END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Increment )

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_dgDiscretization )

  END SUBROUTINE ComputeIncrement_Euler_DG_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      D (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX, iNodeX_X1
    REAL(DP) :: dX2, dX3
    REAL(DP) :: AlphaPls, AlphaMns, AlphaMdl
    REAL(DP) :: G_F(nDOFX_X1,nGF)
    REAL(DP) :: uCF_L(nDOFX_X1,nCF), uCF_R(nDOFX_X1,nCF)
    REAL(DP) :: uPF_L(nDOFX_X1,nPF), uPF_R(nDOFX_X1,nPF)
    REAL(DP) :: P_L(nDOFX_X1), P_R(nDOFX_X1)
    REAL(DP) :: Cs_L(nDOFX_X1), Cs_R(nDOFX_X1)
    REAL(DP) :: G_P(nDOFX,nGF), G_K(nDOFX,nGF)
    REAL(DP) :: uCF_P(nDOFX,nCF), uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: EigVals_L(nDOFX_X1,nCF), EigVals_R(nDOFX_X1,nCF)
    REAL(DP) :: Flux_X1_L(nDOFX_X1,nCF), Flux_X1_R(nDOFX_X1,nCF)
    REAL(DP) :: Flux_X1_q(nDOFX,nCF)
    REAL(DP) :: NumericalFlux(nDOFX_X1,nCF)

    IF( iX_E0(1) .EQ. iX_B0(1) ) RETURN

    CALL TimersStart_Euler( Timer_Euler_Divergence_X1 )

    !$OMP PARALLEL DO PRIVATE &
    !$OMP& ( iX1, iX2, iX3, iCF, iGF, iNodeX, iNodeX_X1, dX2, dX3, &
    !$OMP&   uCF_P, uCF_K, uCF_L, uCF_R, uPF_K, uPF_L, uPF_R, P_K, &
    !$OMP&   P_L, P_R, Cs_L, Cs_R, G_P, G_K, G_F, Flux_X1_q, &
    !$OMP&   Flux_X1_L, Flux_X1_R, NumericalFlux )
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

      CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        CALL ComputePrimitive_Euler &
               ( uCF_K(:,iCF_D ),     &
                 uCF_K(:,iCF_S1),     &
                 uCF_K(:,iCF_S2),     &
                 uCF_K(:,iCF_S3),     &
                 uCF_K(:,iCF_E ),     &
                 uCF_K(:,iCF_Ne),     &
                 uPF_K(:,iPF_D ),     &
                 uPF_K(:,iPF_V1),     &
                 uPF_K(:,iPF_V2),     &
                 uPF_K(:,iPF_V3),     &
                 uPF_K(:,iPF_E ),     &
                 uPF_K(:,iPF_Ne),     &
                 G_K(:,iGF_Gm_dd_11), &
                 G_K(:,iGF_Gm_dd_22), &
                 G_K(:,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

        DO iNodeX = 1, nDOFX

          Flux_X1_q(iNodeX,:) &
            = Flux_X1_Euler &
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
                  G_K  (iNodeX,iGF_Beta_1) )

        END DO

        DO iCF = 1, nCF

          Flux_X1_q(:,iCF) &
            = dX2 * dX3 * WeightsX_q * G_K(:,iGF_Alpha) * G_K(:,iGF_SqrtGm) &
                * Flux_X1_q(:,iCF)

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX, nDOFX, One, dLXdX1_q, nDOFX, &
                   Flux_X1_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

      !---------------------
      ! --- Surface Term ---
      !---------------------

      CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

      ! --- Interpolate Fluid Fields ---

      DO iCF = 1, nCF

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
                 uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                 uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      END DO

      ! --- Interpolate Geometry Fields ---

      ! --- Face States (Average of Left and Right States) ---

      G_F = Zero

      DO iGF = iGF_h_1, iGF_h_3

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_F(:,iGF) = MAX( G_F(:,iGF), SqrtTiny )

      END DO

      CALL ComputeGeometryX_FromScaleFactors( G_F )

      ! --- Lapse Function ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_P(:,iGF_Alpha), 1, Zero, G_F(:,iGF_Alpha), 1 )


      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_K(:,iGF_Alpha), 1, Half, G_F(:,iGF_Alpha), 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      G_F(:,iGF_Alpha) = MAX( G_F(:,iGF_Alpha), SqrtTiny )

      ! --- Left State Primitive, etc. ---

      CALL ComputePrimitive_Euler &
             ( uCF_L(:,iCF_D ),       &
               uCF_L(:,iCF_S1),       &
               uCF_L(:,iCF_S2),       &
               uCF_L(:,iCF_S3),       &
               uCF_L(:,iCF_E ),       &
               uCF_L(:,iCF_Ne),       &
               uPF_L(:,iPF_D ),       &
               uPF_L(:,iPF_V1),       &
               uPF_L(:,iPF_V2),       &
               uPF_L(:,iPF_V3),       &
               uPF_L(:,iPF_E ),       &
               uPF_L(:,iPF_Ne),       &
               G_F  (:,iGF_Gm_dd_11), &
               G_F  (:,iGF_Gm_dd_22), &
               G_F  (:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L )

      DO iNodeX_X1 = 1, nDOFX_X1

        EigVals_L(iNodeX_X1,:) &
          = Eigenvalues_Euler &
              ( uPF_L(iNodeX_X1,iPF_V1),       &
                Cs_L (iNodeX_X1),              &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                uPF_L(iNodeX_X1,iPF_V1),       &
                uPF_L(iNodeX_X1,iPF_V2),       &
                uPF_L(iNodeX_X1,iPF_V3),       &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                G_F  (iNodeX_X1,iGF_Alpha),    &
                G_F  (iNodeX_X1,iGF_Beta_1) )

        Flux_X1_L(iNodeX_X1,:) &
          = Flux_X1_Euler &
              ( uPF_L(iNodeX_X1,iPF_D ),       &
                uPF_L(iNodeX_X1,iPF_V1),       &
                uPF_L(iNodeX_X1,iPF_V2),       &
                uPF_L(iNodeX_X1,iPF_V3),       &
                uPF_L(iNodeX_X1,iPF_E ),       &
                uPF_L(iNodeX_X1,iPF_Ne),       &
                P_L  (iNodeX_X1),              &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                G_F  (iNodeX_X1,iGF_Alpha),    &
                G_F  (iNodeX_X1,iGF_Beta_1) )

      END DO

      ! --- Right State Primitive, etc. ---

      CALL ComputePrimitive_Euler &
             ( uCF_R(:,iCF_D ),       &
               uCF_R(:,iCF_S1),       &
               uCF_R(:,iCF_S2),       &
               uCF_R(:,iCF_S3),       &
               uCF_R(:,iCF_E ),       &
               uCF_R(:,iCF_Ne),       &
               uPF_R(:,iPF_D ),       &
               uPF_R(:,iPF_V1),       &
               uPF_R(:,iPF_V2),       &
               uPF_R(:,iPF_V3),       &
               uPF_R(:,iPF_E ),       &
               uPF_R(:,iPF_Ne),       &
               G_F  (:,iGF_Gm_dd_11), &
               G_F  (:,iGF_Gm_dd_22), &
               G_F  (:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R )

      DO iNodeX_X1 = 1, nDOFX_X1

        EigVals_R(iNodeX_X1,:) &
          = Eigenvalues_Euler &
              ( uPF_R(iNodeX_X1,iPF_V1),       &
                Cs_R (iNodeX_X1),              &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                uPF_R(iNodeX_X1,iPF_V1),       &
                uPF_R(iNodeX_X1,iPF_V2),       &
                uPF_R(iNodeX_X1,iPF_V3),       &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                G_F  (iNodeX_X1,iGF_Alpha),    &
                G_F  (iNodeX_X1,iGF_Beta_1) )

        Flux_X1_R(iNodeX_X1,:) &
          = Flux_X1_Euler &
              ( uPF_R(iNodeX_X1,iPF_D ),       &
                uPF_R(iNodeX_X1,iPF_V1),       &
                uPF_R(iNodeX_X1,iPF_V2),       &
                uPF_R(iNodeX_X1,iPF_V3),       &
                uPF_R(iNodeX_X1,iPF_E ),       &
                uPF_R(iNodeX_X1,iPF_Ne),       &
                P_R  (iNodeX_X1),              &
                G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                G_F  (iNodeX_X1,iGF_Alpha),    &
                G_F  (iNodeX_X1,iGF_Beta_1) )

      END DO

      ! --- Numerical Flux ---

      CALL TimersStart_Euler( Timer_Euler_NumericalFlux )

      DO iNodeX_X1 = 1, nDOFX_X1

        AlphaMns &
          = MAX( Zero, &
                 MAXVAL( - EigVals_L(iNodeX_X1,:) ), &
                 MAXVAL( - EigVals_R(iNodeX_X1,:) ) )

        AlphaPls &
          = MAX( Zero, &
                 MAXVAL( + EigVals_L(iNodeX_X1,:) ), &
                 MAXVAL( + EigVals_R(iNodeX_X1,:) ) )

        AlphaMdl &
          = AlphaMiddle_Euler &
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
                G_F      (iNodeX_X1,iGF_Gm_dd_11), &
                AlphaPls, AlphaMns,                &
                G_F      (iNodeX_X1,iGF_Alpha),    &
                G_F      (iNodeX_X1,iGF_Beta_1) )

        NumericalFlux(iNodeX_X1,:) &
          = NumericalFlux_Euler_X1 &
              ( uCF_L    (iNodeX_X1,:),                 &
                uCF_R    (iNodeX_X1,:),                 &
                Flux_X1_L(iNodeX_X1,:),                 &
                Flux_X1_R(iNodeX_X1,:),                 &
                AlphaPls, AlphaMns, AlphaMdl,           &
                G_F      (iNodeX_X1,iGF_Gm_dd_11),      &
                uPF_L    (iNodeX_X1,iPF_V1),            &
                uPF_R    (iNodeX_X1,iPF_V1),            &
                P_L      (iNodeX_X1),                   &
                P_R      (iNodeX_X1),                   &
                G_F      (iNodeX_X1,iGF_Alpha),         &
                G_F      (iNodeX_X1,iGF_Beta_1),        &
                MAXVAL( D(:,iX1-1,iX2,iX3,iDF_Sh_X2) ), &
                MAXVAL( D(:,iX1  ,iX2,iX3,iDF_Sh_X2) ), &
                MAXVAL( D(:,iX1-1,iX2,iX3,iDF_Sh_X3) ), &
                MAXVAL( D(:,iX1  ,iX2,iX3,iDF_Sh_X3) ) )

      END DO

      CALL TimersStop_Euler( Timer_Euler_NumericalFlux )

      DO iCF = 1, nCF

        NumericalFlux(:,iCF) &
          = dX2 * dX3 * WeightsX_X1 * G_F(:,iGF_Alpha) * G_F(:,iGF_SqrtGm) &
              * NumericalFlux(:,iCF)

      END DO

      ! --- Contribution to This Element ---

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        DO iCF = 1, nCF

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Dn, nDOFX_X1, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      ! --- Contribution to Previous Element ---

      IF( iX1 .GT. iX_B0(1) )THEN

        DO iCF = 1, nCF

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Up, nDOFX_X1, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1-1,iX2,iX3,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    END DO
    END DO
    END DO
    !$OMP END PARALLEL DO

    CALL TimersStop_Euler( Timer_Euler_Divergence_X1 )

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_Divergence_X2 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      D (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX, iNodeX_X2
    REAL(DP) :: dX1, dX3
    REAL(DP) :: AlphaPls, AlphaMns, AlphaMdl
    REAL(DP) :: G_F(nDOFX_X2,nGF)
    REAL(DP) :: uCF_L(nDOFX_X2,nCF), uCF_R(nDOFX_X2,nCF)
    REAL(DP) :: uPF_L(nDOFX_X2,nPF), uPF_R(nDOFX_X2,nPF)
    REAL(DP) :: P_L(nDOFX_X2), P_R(nDOFX_X2)
    REAL(DP) :: Cs_L(nDOFX_X2), Cs_R(nDOFX_X2)
    REAL(DP) :: G_P(nDOFX,nGF), G_K(nDOFX,nGF)
    REAL(DP) :: uCF_P(nDOFX,nCF), uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: EigVals_L(nDOFX_X2,nCF), EigVals_R(nDOFX_X2,nCF)
    REAL(DP) :: Flux_X2_L(nDOFX_X2,nCF), Flux_X2_R(nDOFX_X2,nCF)
    REAL(DP) :: Flux_X2_q(nDOFX,nCF)
    REAL(DP) :: NumericalFlux(nDOFX_X2,nCF)

    IF( iX_E0(2) .EQ. iX_B0(2) ) RETURN

    CALL TimersStart_Euler( Timer_Euler_Divergence_X2 )

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

      CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

      IF( iX2 .LT. iX_E0(2) + 1 )THEN

        CALL ComputePrimitive_Euler &
               ( uCF_K(:,iCF_D ),     &
                 uCF_K(:,iCF_S1),     &
                 uCF_K(:,iCF_S2),     &
                 uCF_K(:,iCF_S3),     &
                 uCF_K(:,iCF_E ),     &
                 uCF_K(:,iCF_Ne),     &
                 uPF_K(:,iPF_D ),     &
                 uPF_K(:,iPF_V1),     &
                 uPF_K(:,iPF_V2),     &
                 uPF_K(:,iPF_V3),     &
                 uPF_K(:,iPF_E ),     &
                 uPF_K(:,iPF_Ne),     &
                 G_K(:,iGF_Gm_dd_11), &
                 G_K(:,iGF_Gm_dd_22), &
                 G_K(:,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

        DO iNodeX = 1, nDOFX

          Flux_X2_q(iNodeX,:) &
            = Flux_X2_Euler &
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

        DO iCF = 1, nCF

          Flux_X2_q(:,iCF) &
            = dX1 * dX3 * WeightsX_q * G_K(:,iGF_Alpha) * G_K(:,iGF_SqrtGm) &
                * Flux_X2_q(:,iCF)

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX, nDOFX, One, dLXdX2_q, nDOFX, &
                   Flux_X2_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

      !---------------------
      ! --- Surface Term ---
      !---------------------

      CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

      ! --- Interpolate Fluid Fields ---

      DO iCF = 1, nCF

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
                 uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
                 uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      END DO

      ! --- Interpolate Geometry Fields ---

      ! --- Face States (Average of Left and Right States) ---

      G_F = Zero

      DO iGF = iGF_h_1, iGF_h_3

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                 G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                 G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_F(:,iGF) = MAX( G_F(:,iGF), SqrtTiny )

      END DO

      CALL ComputeGeometryX_FromScaleFactors( G_F )

      ! --- Lapse Function ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               G_P(:,iGF_Alpha), 1, Zero, G_F(:,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               G_K(:,iGF_Alpha), 1, Half, G_F(:,iGF_Alpha), 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      G_F(:,iGF_Alpha) = MAX( G_F(:,iGF_Alpha), SqrtTiny )

      ! --- Left State Primitive, etc. ---

      CALL ComputePrimitive_Euler &
             ( uCF_L(:,iCF_D ),       &
               uCF_L(:,iCF_S1),       &
               uCF_L(:,iCF_S2),       &
               uCF_L(:,iCF_S3),       &
               uCF_L(:,iCF_E ),       &
               uCF_L(:,iCF_Ne),       &
               uPF_L(:,iPF_D ),       &
               uPF_L(:,iPF_V1),       &
               uPF_L(:,iPF_V2),       &
               uPF_L(:,iPF_V3),       &
               uPF_L(:,iPF_E ),       &
               uPF_L(:,iPF_Ne),       &
               G_F  (:,iGF_Gm_dd_11), &
               G_F  (:,iGF_Gm_dd_22), &
               G_F  (:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L )

      DO iNodeX_X2 = 1, nDOFX_X2

        EigVals_L(iNodeX_X2,:) &
          = Eigenvalues_Euler &
              ( uPF_L(iNodeX_X2,iPF_V2),       &
                Cs_L (iNodeX_X2),              &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                uPF_L(iNodeX_X2,iPF_V1),       &
                uPF_L(iNodeX_X2,iPF_V2),       &
                uPF_L(iNodeX_X2,iPF_V3),       &
                G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_33), &
                G_F  (iNodeX_X2,iGF_Alpha),    &
                G_F  (iNodeX_X2,iGF_Beta_2) )

        Flux_X2_L(iNodeX_X2,:) &
          = Flux_X2_Euler &
              ( uPF_L(iNodeX_X2,iPF_D ),       &
                uPF_L(iNodeX_X2,iPF_V1),       &
                uPF_L(iNodeX_X2,iPF_V2),       &
                uPF_L(iNodeX_X2,iPF_V3),       &
                uPF_L(iNodeX_X2,iPF_E ),       &
                uPF_L(iNodeX_X2,iPF_Ne),       &
                P_L  (iNodeX_X2),              &
                G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_33), &
                G_F  (iNodeX_X2,iGF_Alpha),    &
                G_F  (iNodeX_X2,iGF_Beta_2) )

      END DO

      ! --- Right State Primitive, etc. ---

      CALL ComputePrimitive_Euler &
             ( uCF_R(:,iCF_D ),       &
               uCF_R(:,iCF_S1),       &
               uCF_R(:,iCF_S2),       &
               uCF_R(:,iCF_S3),       &
               uCF_R(:,iCF_E ),       &
               uCF_R(:,iCF_Ne),       &
               uPF_R(:,iPF_D ),       &
               uPF_R(:,iPF_V1),       &
               uPF_R(:,iPF_V2),       &
               uPF_R(:,iPF_V3),       &
               uPF_R(:,iPF_E ),       &
               uPF_R(:,iPF_Ne),       &
               G_F  (:,iGF_Gm_dd_11), &
               G_F  (:,iGF_Gm_dd_22), &
               G_F  (:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), P_R )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R )

      DO iNodeX_X2 = 1, nDOFX_X2

        EigVals_R(iNodeX_X2,:) &
          = Eigenvalues_Euler &
              ( uPF_R(iNodeX_X2,iPF_V2),       &
                Cs_R (iNodeX_X2),              &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                uPF_R(iNodeX_X2,iPF_V1),       &
                uPF_R(iNodeX_X2,iPF_V2),       &
                uPF_R(iNodeX_X2,iPF_V3),       &
                G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_33), &
                G_F  (iNodeX_X2,iGF_Alpha),    &
                G_F  (iNodeX_X2,iGF_Beta_2) )

        Flux_X2_R(iNodeX_X2,:) &
          = Flux_X2_Euler &
              ( uPF_R(iNodeX_X2,iPF_D ),       &
                uPF_R(iNodeX_X2,iPF_V1),       &
                uPF_R(iNodeX_X2,iPF_V2),       &
                uPF_R(iNodeX_X2,iPF_V3),       &
                uPF_R(iNodeX_X2,iPF_E ),       &
                uPF_R(iNodeX_X2,iPF_Ne),       &
                P_R  (iNodeX_X2),              &
                G_F  (iNodeX_X2,iGF_Gm_dd_11), &
                G_F  (iNodeX_X2,iGF_Gm_dd_22), &
                G_F  (iNodeX_X2,iGF_Gm_dd_33), &
                G_F  (iNodeX_X2,iGF_Alpha),    &
                G_F  (iNodeX_X2,iGF_Beta_2) )

      END DO

      ! --- Numerical Flux ---

      CALL TimersStart_Euler( Timer_Euler_NumericalFlux )

      DO iNodeX_X2 = 1, nDOFX_X2

        AlphaMns &
          = MAX( Zero, &
                 MAXVAL( - EigVals_L(iNodeX_X2,:) ), &
                 MAXVAL( - EigVals_R(iNodeX_X2,:) ) )

        AlphaPls &
          = MAX( Zero, &
                 MAXVAL( + EigVals_L(iNodeX_X2,:) ), &
                 MAXVAL( + EigVals_R(iNodeX_X2,:) ) )

        AlphaMdl &
          = AlphaMiddle_Euler &
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
                G_F      (iNodeX_X2,iGF_Gm_dd_22), &
                AlphaPls, AlphaMns,                &
                G_F      (iNodeX_X2,iGF_Alpha),    &
                G_F      (iNodeX_X2,iGF_Beta_2) )

        NumericalFlux(iNodeX_X2,:) &
          = NumericalFlux_Euler_X2 &
              ( uCF_L    (iNodeX_X2,:),                 &
                uCF_R    (iNodeX_X2,:),                 &
                Flux_X2_L(iNodeX_X2,:),                 &
                Flux_X2_R(iNodeX_X2,:),                 &
                AlphaPls, AlphaMns, AlphaMdl,           &
                G_F      (iNodeX_X2,iGF_Gm_dd_22),      &
                uPF_L    (iNodeX_X2,iPF_V2),            &
                uPF_R    (iNodeX_X2,iPF_V2),            &
                P_L      (iNodeX_X2),                   &
                P_R      (iNodeX_X2),                   &
                G_F      (iNodeX_X2,iGF_Alpha),         &
                G_F      (iNodeX_X2,iGF_Beta_2),        &
                MAXVAL( D(:,iX1,iX2-1,iX3,iDF_Sh_X1) ), &
                MAXVAL( D(:,iX1,iX2  ,iX3,iDF_Sh_X1) ), &
                MAXVAL( D(:,iX1,iX2-1,iX3,iDF_Sh_X3) ), &
                MAXVAL( D(:,iX1,iX2  ,iX3,iDF_Sh_X3) ) )

      END DO

      CALL TimersStop_Euler( Timer_Euler_NumericalFlux )

      DO iCF = 1, nCF

        NumericalFlux(:,iCF) &
          = dX1 * dX3 * WeightsX_X2 * G_F(:,iGF_Alpha) * G_F(:,iGF_SqrtGm) &
              * NumericalFlux(:,iCF)

      END DO

      ! --- Contribution to This Element ---

      IF( iX2 .LT. iX_E0(2) + 1 )THEN

        DO iCF = 1, nCF

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Dn, nDOFX_X2, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      ! --- Contribution to Previous Element ---

      IF( iX2 .GT. iX_B0(2) )THEN

        DO iCF = 1, nCF

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Up, nDOFX_X2, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2-1,iX3,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Divergence_X2 )

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      D (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX, iNodeX_X3
    REAL(DP) :: dX1, dX2
    REAL(DP) :: AlphaPls, AlphaMns, AlphaMdl
    REAL(DP) :: G_F(nDOFX_X3,nGF)
    REAL(DP) :: uCF_L(nDOFX_X3,nCF), uCF_R(nDOFX_X3,nCF)
    REAL(DP) :: uPF_L(nDOFX_X3,nPF), uPF_R(nDOFX_X3,nPF)
    REAL(DP) :: P_L(nDOFX_X3), P_R(nDOFX_X3)
    REAL(DP) :: Cs_L(nDOFX_X3), Cs_R(nDOFX_X3)
    REAL(DP) :: G_P(nDOFX,nGF), G_K(nDOFX,nGF)
    REAL(DP) :: uCF_P(nDOFX,nCF), uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: EigVals_L(nDOFX_X3,nCF), EigVals_R(nDOFX_X3,nCF)
    REAL(DP) :: Flux_X3_L(nDOFX_X3,nCF), Flux_X3_R(nDOFX_X3,nCF)
    REAL(DP) :: Flux_X3_q(nDOFX,nCF)
    REAL(DP) :: NumericalFlux(nDOFX_X3,nCF)

    IF( iX_E0(3) .EQ. iX_B0(3) ) RETURN

    CALL TimersStart_Euler( Timer_Euler_Divergence_X3 )

    DO iX3 = iX_B0(3), iX_E0(3) + 1
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX2 = MeshX(2) % Width(iX2)
      dX1 = MeshX(1) % Width(iX1)

      DO iCF = 1, nCF

        uCF_P(:,iCF) = U(:,iX1,iX2,iX3-1,iCF)
        uCF_K(:,iCF) = U(:,iX1,iX2,iX3  ,iCF)

      END DO

      DO iGF = 1, nGF

        G_P(:,iGF) = G(:,iX1,iX2,iX3-1,iGF)
        G_K(:,iGF) = G(:,iX1,iX2,iX3  ,iGF)

      END DO

      !--------------------
      ! --- Volume Term ---
      !--------------------

      CALL TimersStart_Euler( Timer_Euler_VolumeTerm )

      IF( iX3 .LT. iX_E0(3) + 1 )THEN

        CALL ComputePrimitive_Euler &
               ( uCF_K(:,iCF_D ),     &
                 uCF_K(:,iCF_S1),     &
                 uCF_K(:,iCF_S2),     &
                 uCF_K(:,iCF_S3),     &
                 uCF_K(:,iCF_E ),     &
                 uCF_K(:,iCF_Ne),     &
                 uPF_K(:,iPF_D ),     &
                 uPF_K(:,iPF_V1),     &
                 uPF_K(:,iPF_V2),     &
                 uPF_K(:,iPF_V3),     &
                 uPF_K(:,iPF_E ),     &
                 uPF_K(:,iPF_Ne),     &
                 G_K(:,iGF_Gm_dd_11), &
                 G_K(:,iGF_Gm_dd_22), &
                 G_K(:,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

        DO iNodeX = 1, nDOFX

          Flux_X3_q(iNodeX,:) &
            = Flux_X3_Euler &
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
                  G_K  (iNodeX,iGF_Beta_3) )

        END DO

        DO iCF = 1, nCF

          Flux_X3_q(:,iCF) &
            = dX1 * dX2 * WeightsX_q * G_K(:,iGF_Alpha) * G_K(:,iGF_SqrtGm) &
                * Flux_X3_q(:,iCF)

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX, nDOFX, One, dLXdX3_q, nDOFX, &
                   Flux_X3_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      CALL TimersStop_Euler( Timer_Euler_VolumeTerm )

      !---------------------
      ! --- Surface Term ---
      !---------------------

      CALL TimersStart_Euler( Timer_Euler_SurfaceTerm )

      ! --- Interpolate Fluid Fields ---

      DO iCF = 1, nCF

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Up, nDOFX_X3, &
                 uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
                 uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      END DO

      ! --- Interpolate Geometry Fields ---

      ! --- Face States (Average of Left and Right States) ---

      G_F = Zero

      DO iGF = iGF_h_1, iGF_h_3

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                 G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                 G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_F(:,iGF) = MAX( G_F(:,iGF), SqrtTiny )

      END DO

      CALL ComputeGeometryX_FromScaleFactors( G_F )

      ! --- Lapse Function ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
               G_P(:,iGF_Alpha), 1, Zero, G_F(:,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
               G_K(:,iGF_Alpha), 1, Half, G_F(:,iGF_Alpha), 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      G_F(:,iGF_Alpha) = MAX( G_F(:,iGF_Alpha), SqrtTiny )

      ! --- Left State Primitive, etc. ---

      CALL ComputePrimitive_Euler &
             ( uCF_L(:,iCF_D ),       &
               uCF_L(:,iCF_S1),       &
               uCF_L(:,iCF_S2),       &
               uCF_L(:,iCF_S3),       &
               uCF_L(:,iCF_E ),       &
               uCF_L(:,iCF_Ne),       &
               uPF_L(:,iPF_D ),       &
               uPF_L(:,iPF_V1),       &
               uPF_L(:,iPF_V2),       &
               uPF_L(:,iPF_V3),       &
               uPF_L(:,iPF_E ),       &
               uPF_L(:,iPF_Ne),       &
               G_F  (:,iGF_Gm_dd_11), &
               G_F  (:,iGF_Gm_dd_22), &
               G_F  (:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L )

      DO iNodeX_X3 = 1, nDOFX_X3

        EigVals_L(iNodeX_X3,:) &
          = Eigenvalues_Euler &
              ( uPF_L(iNodeX_X3,iPF_V3),       &
                Cs_L (iNodeX_X3),              &
                G_F  (iNodeX_X3,iGF_Gm_dd_33), &
                uPF_L(iNodeX_X3,iPF_V1),       &
                uPF_L(iNodeX_X3,iPF_V2),       &
                uPF_L(iNodeX_X3,iPF_V3),       &
                G_F  (iNodeX_X3,iGF_Gm_dd_11), &
                G_F  (iNodeX_X3,iGF_Gm_dd_22), &
                G_F  (iNodeX_X3,iGF_Gm_dd_33), &
                G_F  (iNodeX_X3,iGF_Alpha),    &
                G_F  (iNodeX_X3,iGF_Beta_3) )

        Flux_X3_L(iNodeX_X3,:) &
          = Flux_X3_Euler &
              ( uPF_L(iNodeX_X3,iPF_D ),       &
                uPF_L(iNodeX_X3,iPF_V1),       &
                uPF_L(iNodeX_X3,iPF_V2),       &
                uPF_L(iNodeX_X3,iPF_V3),       &
                uPF_L(iNodeX_X3,iPF_E ),       &
                uPF_L(iNodeX_X3,iPF_Ne),       &
                P_L  (iNodeX_X3),              &
                G_F  (iNodeX_X3,iGF_Gm_dd_11), &
                G_F  (iNodeX_X3,iGF_Gm_dd_22), &
                G_F  (iNodeX_X3,iGF_Gm_dd_33), &
                G_F  (iNodeX_X3,iGF_Alpha),    &
                G_F  (iNodeX_X3,iGF_Beta_3) )

      END DO

      ! --- Right State Primitive, etc. ---

      CALL ComputePrimitive_Euler &
             ( uCF_R(:,iCF_D ),       &
               uCF_R(:,iCF_S1),       &
               uCF_R(:,iCF_S2),       &
               uCF_R(:,iCF_S3),       &
               uCF_R(:,iCF_E ),       &
               uCF_R(:,iCF_Ne),       &
               uPF_R(:,iPF_D ),       &
               uPF_R(:,iPF_V1),       &
               uPF_R(:,iPF_V2),       &
               uPF_R(:,iPF_V3),       &
               uPF_R(:,iPF_E ),       &
               uPF_R(:,iPF_Ne),       &
               G_F  (:,iGF_Gm_dd_11), &
               G_F  (:,iGF_Gm_dd_22), &
               G_F  (:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), P_R )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R )

      DO iNodeX_X3 = 1, nDOFX_X3

        EigVals_R(iNodeX_X3,:) &
          = Eigenvalues_Euler &
              ( uPF_R(iNodeX_X3,iPF_V3),       &
                Cs_R (iNodeX_X3),              &
                G_F  (iNodeX_X3,iGF_Gm_dd_33), &
                uPF_R(iNodeX_X3,iPF_V1),       &
                uPF_R(iNodeX_X3,iPF_V2),       &
                uPF_R(iNodeX_X3,iPF_V3),       &
                G_F  (iNodeX_X3,iGF_Gm_dd_11), &
                G_F  (iNodeX_X3,iGF_Gm_dd_22), &
                G_F  (iNodeX_X3,iGF_Gm_dd_33), &
                G_F  (iNodeX_X3,iGF_Alpha),    &
                G_F  (iNodeX_X3,iGF_Beta_3) )

        Flux_X3_R(iNodeX_X3,:) &
          = Flux_X3_Euler &
              ( uPF_R(iNodeX_X3,iPF_D ),       &
                uPF_R(iNodeX_X3,iPF_V1),       &
                uPF_R(iNodeX_X3,iPF_V2),       &
                uPF_R(iNodeX_X3,iPF_V3),       &
                uPF_R(iNodeX_X3,iPF_E ),       &
                uPF_R(iNodeX_X3,iPF_Ne),       &
                P_R  (iNodeX_X3),              &
                G_F  (iNodeX_X3,iGF_Gm_dd_11), &
                G_F  (iNodeX_X3,iGF_Gm_dd_22), &
                G_F  (iNodeX_X3,iGF_Gm_dd_33), &
                G_F  (iNodeX_X3,iGF_Alpha),    &
                G_F  (iNodeX_X3,iGF_Beta_3) )

      END DO

      ! --- Numerical Flux ---

      CALL TimersStart_Euler( Timer_Euler_NumericalFlux )

      DO iNodeX_X3 = 1, nDOFX_X3

        AlphaMns &
          = MAX( Zero, &
                 MAXVAL( - EigVals_L(iNodeX_X3,:) ), &
                 MAXVAL( - EigVals_R(iNodeX_X3,:) ) )

        AlphaPls &
          = MAX( Zero, &
                 MAXVAL( + EigVals_L(iNodeX_X3,:) ), &
                 MAXVAL( + EigVals_R(iNodeX_X3,:) ) )

        AlphaMdl &
          = AlphaMiddle_Euler &
              ( uCF_L    (iNodeX_X3,iCF_D ),       &
                uCF_L    (iNodeX_X3,iCF_S3),       &
                uCF_L    (iNodeX_X3,iCF_E ),       &
                Flux_X3_L(iNodeX_X3,iCF_D ),       &
                Flux_X3_L(iNodeX_X3,iCF_S3),       &
                Flux_X3_L(iNodeX_X3,iCF_E ),       &
                uCF_R    (iNodeX_X3,iCF_D ),       &
                uCF_R    (iNodeX_X3,iCF_S3),       &
                uCF_R    (iNodeX_X3,iCF_E ),       &
                Flux_X3_R(iNodeX_X3,iCF_D ),       &
                Flux_X3_R(iNodeX_X3,iCF_S3),       &
                Flux_X3_R(iNodeX_X3,iCF_E ),       &
                G_F      (iNodeX_X3,iGF_Gm_dd_33), &
                AlphaPls, AlphaMns,                &
                G_F      (iNodeX_X3,iGF_Alpha),    &
                G_F      (iNodeX_X3,iGF_Beta_3) )

        NumericalFlux(iNodeX_X3,:) &
          = NumericalFlux_Euler_X3 &
              ( uCF_L    (iNodeX_X3,:),                 &
                uCF_R    (iNodeX_X3,:),                 &
                Flux_X3_L(iNodeX_X3,:),                 &
                Flux_X3_R(iNodeX_X3,:),                 &
                AlphaPls, AlphaMns, AlphaMdl,           &
                G_F      (iNodeX_X3,iGF_Gm_dd_33),      &
                uPF_L    (iNodeX_X3,iPF_V3),            &
                uPF_R    (iNodeX_X3,iPF_V3),            &
                P_L      (iNodeX_X3),                   &
                P_R      (iNodeX_X3),                   &
                G_F      (iNodeX_X3,iGF_Alpha),         &
                G_F      (iNodeX_X3,iGF_Beta_3),        &
                MAXVAL( D(:,iX1,iX2,iX3-1,iDF_Sh_X1) ), &
                MAXVAL( D(:,iX1,iX2,iX3  ,iDF_Sh_X1) ), &
                MAXVAL( D(:,iX1,iX2,iX3-1,iDF_Sh_X2) ), &
                MAXVAL( D(:,iX1,iX2,iX3  ,iDF_Sh_X2) ) )

      END DO

      CALL TimersStop_Euler( Timer_Euler_NumericalFlux )

      DO iCF = 1, nCF

        NumericalFlux(:,iCF) &
          = dX1 * dX2 * WeightsX_X3 * G_F(:,iGF_Alpha) * G_F(:,iGF_SqrtGm) &
              * NumericalFlux(:,iCF)

      END DO

      ! --- Contribution to This Element ---

      IF( iX3 .LT. iX_E0(3) + 1 )THEN

        DO iCF = 1, nCF

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Dn, nDOFX_X3, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      ! --- Contribution to Previous Element ---

      IF( iX3 .GT. iX_B0(3) )THEN

        DO iCF = 1, nCF

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Up, nDOFX_X3, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3-1,iCF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

      END IF

      CALL TimersStop_Euler( Timer_Euler_SurfaceTerm )

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_Divergence_X3 )

  END SUBROUTINE ComputeIncrement_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    CALL TimersStart_Euler( Timer_Euler_Geometry )

#if defined HYDRO_RELATIVISTIC

    CALL ComputeIncrement_Geometry_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#else

    CALL ComputeIncrement_Geometry_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#endif

    CALL TimersStop_Euler( Timer_Euler_Geometry )

  END SUBROUTINE ComputeIncrement_Geometry


  SUBROUTINE ComputeIncrement_Geometry_NonRelativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX
    REAL(DP) :: dX1, dX2
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: dh2dX1(nDOFX), dh3dX1(nDOFX), dh3dX2(nDOFX)
    REAL(DP) :: Stress(nDOFX,3)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: G_K(nDOFX,nGF)
    REAL(DP) :: G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF)
    REAL(DP) :: G_P_X2(nDOFX,nGF), G_N_X2(nDOFX,nGF)
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF)
    REAL(DP) :: G_X2_Dn(nDOFX_X2,nGF), G_X2_Up(nDOFX_X2,nGF)

    IF( TRIM( CoordinateSystem ) == 'CARTESIAN' ) RETURN

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX2 = MeshX(2) % Width(iX2)
      dX1 = MeshX(1) % Width(iX1)

!      print*,"iX1, iX2, iX3 = ", iX1, iX2, iX3

      DO iCF = 1, nCF

        uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

      END DO

      DO iGF = 1, nGF

        G_K   (:,iGF) = G(:,iX1,  iX2,iX3,iGF)
        G_P_X1(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
        G_N_X1(:,iGF) = G(:,iX1+1,iX2,iX3,iGF)

      END DO

      IF( nDimsX .GT. 1 )THEN

        DO iGF = 1, nGF

          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)

        END DO

      END IF

      CALL ComputePrimitive_Euler &
             ( uCF_K(:,iCF_D ), &
               uCF_K(:,iCF_S1), &
               uCF_K(:,iCF_S2), &
               uCF_K(:,iCF_S3), &
               uCF_K(:,iCF_E ), &
               uCF_K(:,iCF_Ne), &
               uPF_K(:,iPF_D ), &
               uPF_K(:,iPF_V1), &
               uPF_K(:,iPF_V2), &
               uPF_K(:,iPF_V3), &
               uPF_K(:,iPF_E ), &
               uPF_K(:,iPF_Ne), &
               G_K(:,iGF_Gm_dd_11), &
               G_K(:,iGF_Gm_dd_22), &
               G_K(:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

      DO iNodeX = 1, nDOFX

        Stress(iNodeX,1:3) &
          = StressTensor_Diagonal_Euler &
              ( uCF_K(iNodeX,iCF_S1), &
                uCF_K(iNodeX,iCF_S2), &
                uCF_K(iNodeX,iCF_S3), &
                uPF_K(iNodeX,iPF_V1), &
                uPF_K(iNodeX,iPF_V2), &
                uPF_K(iNodeX,iPF_V3), &
                P_K  (iNodeX) )

      END DO

      ! --- Scale Factor Derivatives wrt X1 ---

      ! --- Face States (Average of Left and Right States) ---

      DO iGF = iGF_h_2, iGF_h_3

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P_X1(:,iGF), 1, Zero, G_X1_Dn(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K   (:,iGF), 1, Half, G_X1_Dn(:,iGF), 1 )

        G_X1_Dn(1:nDOFX_X1,iGF) &
          = MAX( G_X1_Dn(1:nDOFX_X1,iGF), SqrtTiny )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_K   (:,iGF), 1, Zero, G_X1_Up(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_N_X1(:,iGF), 1, Half, G_X1_Up(:,iGF), 1 )

        G_X1_Up(1:nDOFX_X1,iGF) &
          = MAX( G_X1_Up(1:nDOFX_X1,iGF), SqrtTiny )

      END DO

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Dn(:,iGF_h_2), 1,  One, dh2dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_2), 1,  One, dh2dX1, 1 )

      dh2dx1 = dh2dx1 / ( WeightsX_q(:) * dX1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1(:) * G_X1_Dn(:,iGF_h_3), 1,  One, dh3dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX1, 1 )

      dh3dx1 = dh3dx1 / ( WeightsX_q(:) * dX1 )

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1) &
            + ( Stress(:,2) * dh2dX1(:) ) / G_K(:,iGF_h_2)  &
            + ( Stress(:,3) * dh3dX1(:) ) / G_K(:,iGF_h_3)

      IF( nDimsX .GT. 1 )THEN

        ! --- Scale Factor Derivatives wrt X2 ---

        ! --- Face States (Average of Left and Right States) ---

        DO iGF = iGF_h_3, iGF_h_3

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_P_X2(:,iGF), 1, Zero, G_X2_Dn(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_K   (:,iGF), 1, Half, G_X2_Dn(:,iGF), 1 )

          G_X2_Dn(1:nDOFX_X2,iGF) &
            = MAX( G_X2_Dn(1:nDOFX_X2,iGF), SqrtTiny )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_K   (:,iGF), 1, Zero, G_X2_Up(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_N_X2(:,iGF), 1, Half, G_X2_Up(:,iGF), 1 )

          G_X2_Up(1:nDOFX_X2,iGF) &
            = MAX( G_X2_Up(1:nDOFX_X2,iGF), SqrtTiny )

        END DO

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                  WeightsX_X2(:) * G_X2_Up(:,iGF_h_3), 1, Zero, dh3dX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                  WeightsX_X2(:) * G_X2_Dn(:,iGF_h_3), 1,  One, dh3dX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                  WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX2, 1 )

        dh3dx2 = dh3dx2 / ( WeightsX_q(:) * dX2 )

        dU(:,iX1,iX2,iX3,iCF_S2) &
          = dU(:,iX1,iX2,iX3,iCF_S2) &
              + ( Stress(:,3) * dh3dX2(:) ) / G_K(:,iGF_h_3)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_Geometry_NonRelativistic


  SUBROUTINE ComputeIncrement_Geometry_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX, iDim, jDim
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: dh1dX1(nDOFX), dh2dX1(nDOFX), dh3dX1(nDOFX), &
                dh1dX2(nDOFX), dh2dX2(nDOFX), dh3dX2(nDOFX), &
                dh1dX3(nDOFX), dh2dX3(nDOFX), dh3dX3(nDOFX)
    REAL(DP) :: db1dX1(nDOFX), db2dX1(nDOFX), db3dX1(nDOFX), &
                db1dX2(nDOFX), db2dX2(nDOFX), db3dX2(nDOFX), &
                db1dX3(nDOFX), db2dX3(nDOFX), db3dX3(nDOFX)
    REAL(DP) :: dadx1(nDOFX), dadx2(nDOFX), dadx3(nDOFX)
    REAL(DP) :: Stress(nDOFX,3)
    REAL(DP) :: uCF_K(nDOFX,nCF), uPF_K(nDOFX,nPF), G_K(nDOFX,nGF)
    REAL(DP) :: G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF), &
                G_P_X2(nDOFX,nGF), G_N_X2(nDOFX,nGF), &
                G_P_X3(nDOFX,nGF), G_N_X3(nDOFX,nGF)
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF), &
                G_X2_Dn(nDOFX_X2,nGF), G_X2_Up(nDOFX_X2,nGF), &
                G_X3_Dn(nDOFX_X3,nGF), G_X3_Up(nDOFX_X3,nGF)

    REAL(DP) :: EnergyDensitySourceTerms(nDOFX,iX_B0(1):iX_E0(1), &
                                               iX_B0(2):iX_E0(2), &
                                               iX_B0(3):iX_E0(3),7)

    REAL(DP) :: DivGridVolume      (nDOFX)
    REAL(DP) :: PressureTensorTrace(nDOFX)
    REAL(DP) :: PressureTensor     (nDOFX,3,3)
    REAL(DP) :: Xij                (nDOFX,3,3)
    REAL(DP) :: Christoffel3D_X1   (nDOFX,3,3)
    REAL(DP) :: Christoffel3D_X2   (nDOFX,3,3)
    REAL(DP) :: Christoffel3D_X3   (nDOFX,3,3)
    REAL(DP) :: Christoffel_X1     (nDOFX,3,3)
    REAL(DP) :: Christoffel_X2     (nDOFX,3,3)
    REAL(DP) :: Christoffel_X3     (nDOFX,3,3)

    REAL(DP) :: GradPsi (nDOFX)
    REAL(DP) :: GradPsiF(nDOFX)
    REAL(DP) :: X1     (nDOFX)

    dadx1  = Zero
    dadx2  = Zero
    dadx3  = Zero
    dh1dX1 = Zero
    dh2dX1 = Zero
    dh3dX1 = Zero
    dh1dX2 = Zero
    dh2dX2 = Zero
    dh3dX2 = Zero
    dh1dX3 = Zero
    dh2dX3 = Zero
    dh3dX3 = Zero
    db1dX1 = Zero
    db2dX1 = Zero
    db3dX1 = Zero
    db1dX2 = Zero
    db2dX2 = Zero
    db3dX2 = Zero
    db1dX3 = Zero
    db2dX3 = Zero
    db3dX3 = Zero
    PressureTensor = Zero

    GradPsi = Zero

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

        G_K   (:,iGF) = G(:,iX1,  iX2,iX3,iGF)
        G_P_X1(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
        G_N_X1(:,iGF) = G(:,iX1+1,iX2,iX3,iGF)

      END DO

      IF     ( nDimsX .EQ. 2 )THEN

        DO iGF = 1, nGF

          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)

        END DO

      ELSE IF( nDimsX .EQ. 3 )THEN

        DO iGF = 1, nGF

          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)
          G_P_X3(:,iGF) = G(:,iX1,iX2,iX3-1,iGF)
          G_N_X3(:,iGF) = G(:,iX1,iX2,iX3+1,iGF)

        END DO

      END IF

      CALL ComputePrimitive_Euler &
           ( uCF_K(:,iCF_D ),     &
             uCF_K(:,iCF_S1),     &
             uCF_K(:,iCF_S2),     &
             uCF_K(:,iCF_S3),     &
             uCF_K(:,iCF_E ),     &
             uCF_K(:,iCF_Ne),     &
             uPF_K(:,iPF_D ),     &
             uPF_K(:,iPF_V1),     &
             uPF_K(:,iPF_V2),     &
             uPF_K(:,iPF_V3),     &
             uPF_K(:,iPF_E ),     &
             uPF_K(:,iPF_Ne),     &
             G_K(:,iGF_Gm_dd_11), &
             G_K(:,iGF_Gm_dd_22), &
             G_K(:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_K(:,iPF_D), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

      ! --- Compute P^{ij} ---

      DO iNodeX = 1, nDOFX

        PressureTensor(:,1,1) &
          = ( uPF_K(:,iPF_V1) * uCF_K(:,iCF_S1) + P_K ) / G_K(:,iGF_Gm_dd_11)

        PressureTensor(:,1,2) &
          = uPF_K(:,iPF_V1) * uCF_K(:,iCF_S2) / G_K(:,iGF_Gm_dd_22)

        PressureTensor(:,1,3) &
          = uPF_K(:,iPF_V1) * uCF_K(:,iCF_S3) / G_K(:,iGF_Gm_dd_33)

        PressureTensor(:,2,1) &
          = uPF_K(:,iPF_V2) * uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11)

        PressureTensor(:,2,2) &
          = ( uPF_K(:,iPF_V2) * uCF_K(:,iCF_S2) + P_K ) / G_K(:,iGF_Gm_dd_22)

        PressureTensor(:,2,3) &
          = uPF_K(:,iPF_V2) * uCF_K(:,iCF_S3) / G_K(:,iGF_Gm_dd_33)

        PressureTensor(:,3,1) &
          = uPF_K(:,iPF_V3) * uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11)

        PressureTensor(:,3,2) &
          = uPF_K(:,iPF_V3) * uCF_K(:,iCF_S2) / G_K(:,iGF_Gm_dd_22)

        PressureTensor(:,3,3) &
          = ( uPF_K(:,iPF_V3) * uCF_K(:,iCF_S3) + P_K ) / G_K(:,iGF_Gm_dd_33)

      END DO

      ! --- Redundant calculation. Will modify code so this can be
      !     removed later ---

      DO iNodeX = 1, nDOFX

        Stress(iNodeX,:) &
          = StressTensor_Diagonal_Euler &
              ( uCF_K(iNodeX,iCF_S1), &
                uCF_K(iNodeX,iCF_S2), &
                uCF_K(iNodeX,iCF_S3), &
                uPF_K(iNodeX,iPF_V1), &
                uPF_K(iNodeX,iPF_V2), &
                uPF_K(iNodeX,iPF_V3), &
                P_K  (iNodeX) )

      END DO

      ! --- Scale factor derivatives wrt X1 ---

      ! --- Interpolation ---

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_h_1), G_K(:,iGF_h_1), G_N_X1(:,iGF_h_1), &
               G_X1_Dn(:,iGF_h_1), G_X1_Up(:,iGF_h_1) )

      G_X1_Dn(:,iGF_h_1) = MAX( G_X1_Dn(:,iGF_h_1), SqrtTiny )
      G_X1_Up(:,iGF_h_1) = MAX( G_X1_Up(:,iGF_h_1), SqrtTiny )

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_h_2), G_K(:,iGF_h_2), G_N_X1(:,iGF_h_2), &
               G_X1_Dn(:,iGF_h_2), G_X1_Up(:,iGF_h_2) )

      G_X1_Dn(:,iGF_h_2) = MAX( G_X1_Dn(:,iGF_h_2), SqrtTiny )
      G_X1_Up(:,iGF_h_2) = MAX( G_X1_Up(:,iGF_h_2), SqrtTiny )

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_h_3), G_K(:,iGF_h_3), G_N_X1(:,iGF_h_3), &
               G_X1_Dn(:,iGF_h_3), G_X1_Up(:,iGF_h_3) )

      G_X1_Dn(:,iGF_h_3) = MAX( G_X1_Dn(:,iGF_h_3), SqrtTiny )
      G_X1_Up(:,iGF_h_3) = MAX( G_X1_Up(:,iGF_h_3), SqrtTiny )

      ! --- Differentiation ---

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_h_1), G_X1_Dn(:,iGF_h_1), G_K(:,iGF_h_1), &
               dh1dX1 )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_h_2), G_X1_Dn(:,iGF_h_2), G_K(:,iGF_h_2), &
               dh2dX1 )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_h_3), G_X1_Dn(:,iGF_h_3), G_K(:,iGF_h_3), &
               dh3dX1 )

      ! --- Shift vector derivative wrt X1 ---

      ! --- Interpolation ---

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_Beta_1), G_K(:,iGF_Beta_1), G_N_X1(:,iGF_Beta_1), &
               G_X1_Dn(:,iGF_Beta_1), G_X1_Up(:,iGF_Beta_1) )

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_Beta_2), G_K(:,iGF_Beta_2), G_N_X1(:,iGF_Beta_2), &
               G_X1_Dn(:,iGF_Beta_2), G_X1_Up(:,iGF_Beta_2) )

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_Beta_3), G_K(:,iGF_Beta_3), G_N_X1(:,iGF_Beta_3), &
               G_X1_Dn(:,iGF_Beta_3), G_X1_Up(:,iGF_Beta_3) )

      ! --- Differentiation ---

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Beta_1), G_X1_Dn(:,iGF_Beta_1), &
               G_K(:,iGF_Beta_1), &
               db1dX1 )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Beta_2), G_X1_Dn(:,iGF_Beta_2), &
               G_K(:,iGF_Beta_2), &
               db2dX1 )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Beta_3), G_X1_Dn(:,iGF_Beta_3), &
               G_K(:,iGF_Beta_3), &
               db3dX1 )

      ! --- Lapse function derivative wrt X1 ---

      ! --- Interpolation ---

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_Alpha), G_K(:,iGF_Alpha), G_N_X1(:,iGF_Alpha), &
               G_X1_Dn(:,iGF_Alpha), G_X1_Up(:,iGF_Alpha) )

      G_X1_Dn(:,iGF_Alpha) = MAX( G_X1_Dn(:,iGF_Alpha), SqrtTiny )
      G_X1_Up(:,iGF_Alpha) = MAX( G_X1_Up(:,iGF_Alpha), SqrtTiny )

      ! --- Diffentiation ---

      CALL ComputeDerivative &
             ( nDOFX_X1,dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Alpha), G_X1_Dn(:,iGF_Alpha), &
               G_K(:,iGF_Alpha), &
               dadX1 )

      ! --- Momentum density source term (S1) ---

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1)                                &
            + G_K(:,iGF_Alpha)                                    &
                * ( ( Stress(:,1) * dh1dX1 ) / G_K(:,iGF_h_1)     &
                    + ( Stress(:,2) * dh2dX1 ) / G_K(:,iGF_h_2)   &
                    + ( Stress(:,3) * dh3dX1 ) / G_K(:,iGF_h_3) ) &
            + uCF_K(:,iCF_S1) * db1dX1 &
                + uCF_K(:,iCF_S2) * db2dX1 &
                + uCF_K(:,iCF_S3) * db3dX1 &
            - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx1

      IF( nDimsX .GT. 1 )THEN

        ! --- Scale factor derivatives wrt X2 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2(:,iGF_h_1), G_K(:,iGF_h_1), G_N_X2(:,iGF_h_1), &
                 G_X2_Dn(:,iGF_h_1), G_X2_Up(:,iGF_h_1) )

        G_X2_Dn(:,iGF_h_1) = MAX( G_X2_Dn(:,iGF_h_1), SqrtTiny )
        G_X2_Up(:,iGF_h_1) = MAX( G_X2_Up(:,iGF_h_1), SqrtTiny )

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2(:,iGF_h_2), G_K(:,iGF_h_2), G_N_X2(:,iGF_h_2), &
                 G_X2_Dn(:,iGF_h_2), G_X2_Up(:,iGF_h_2) )

        G_X2_Dn(:,iGF_h_2) = MAX( G_X2_Dn(:,iGF_h_2), SqrtTiny )
        G_X2_Up(:,iGF_h_2) = MAX( G_X2_Up(:,iGF_h_2), SqrtTiny )

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2(:,iGF_h_3), G_K(:,iGF_h_3), G_N_X2(:,iGF_h_3), &
                 G_X2_Dn(:,iGF_h_3), G_X2_Up(:,iGF_h_3) )

        G_X2_Dn(:,iGF_h_3) = MAX( G_X2_Dn(:,iGF_h_3), SqrtTiny )
        G_X2_Up(:,iGF_h_3) = MAX( G_X2_Up(:,iGF_h_3), SqrtTiny )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_h_1), G_X1_Dn(:,iGF_h_1), G_K(:,iGF_h_1), &
                 dh1dX2 )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_h_2), G_X1_Dn(:,iGF_h_2), G_K(:,iGF_h_2), &
                 dh2dX2 )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_h_3), G_X1_Dn(:,iGF_h_3), G_K(:,iGF_h_3), &
                 dh3dX2 )

        ! --- Shift vector derivatives wrt X2 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2(:,iGF_Beta_1), G_K(:,iGF_Beta_1), &
                 G_N_X2(:,iGF_Beta_1), &
                 G_X2_Dn(:,iGF_Beta_1), G_X2_Up(:,iGF_Beta_1) )

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2(:,iGF_Beta_2), G_K(:,iGF_Beta_2), &
                 G_N_X2(:,iGF_Beta_2), &
                 G_X2_Dn(:,iGF_Beta_2), G_X2_Up(:,iGF_Beta_2) )

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2(:,iGF_Beta_3), G_K(:,iGF_Beta_3), &
                 G_N_X2(:,iGF_Beta_3), &
                 G_X2_Dn(:,iGF_Beta_3), G_X2_Up(:,iGF_Beta_3) )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Beta_1), G_X1_Dn(:,iGF_Beta_1), &
                 G_K(:,iGF_Beta_1), &
                 db1dX2 )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Beta_2), G_X1_Dn(:,iGF_Beta_2), &
                 G_K(:,iGF_Beta_2), &
                 db2dX2 )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Beta_3), G_X1_Dn(:,iGF_Beta_3), &
                 G_K(:,iGF_Beta_3), &
                 db3dX2 )

        ! --- Lapse function derivative wrt X2 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2(:,iGF_Alpha), G_K(:,iGF_Alpha), &
                 G_N_X2(:,iGF_Alpha), &
                 G_X2_Dn(:,iGF_Alpha), G_X2_Up(:,iGF_Alpha) )

        G_X2_Dn(:,iGF_Alpha) = MAX( G_X2_Dn(:,iGF_Alpha), SqrtTiny )
        G_X2_Up(:,iGF_Alpha) = MAX( G_X2_Up(:,iGF_Alpha), SqrtTiny )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Alpha), G_X1_Dn(:,iGF_Alpha), &
                 G_K(:,iGF_Alpha), &
                 dadX2 )

        ! --- Momentum density source term (S2) ---

        dU(:,iX1,iX2,iX3,iCF_S2) &
          = dU(:,iX1,iX2,iX3,iCF_S2)                                &
              + G_K(:,iGF_Alpha)                                    &
                  * ( ( Stress(:,1) * dh1dX2 ) / G_K(:,iGF_h_1)     &
                      + ( Stress(:,2) * dh2dX2 ) / G_K(:,iGF_h_2)   &
                      + ( Stress(:,3) * dh3dX2 ) / G_K(:,iGF_h_3) ) &
              + uCF_K(:,iCF_S1) * db1dX2 &
                  + uCF_K(:,iCF_S2) * db2dX2 &
                  + uCF_K(:,iCF_S3) * db3dX2 &
              - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx2

      END IF

      IF( nDimsX .GT. 2 )THEN

        ! --- Scale factor derivatives wrt X3 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3(:,iGF_h_1), G_K(:,iGF_h_1), G_N_X3(:,iGF_h_1), &
                 G_X3_Dn(:,iGF_h_1), G_X3_Up(:,iGF_h_1) )

        G_X3_Dn(:,iGF_h_1) = MAX( G_X3_Dn(:,iGF_h_1), SqrtTiny )
        G_X3_Up(:,iGF_h_1) = MAX( G_X3_Up(:,iGF_h_1), SqrtTiny )

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3(:,iGF_h_2), G_K(:,iGF_h_2), G_N_X3(:,iGF_h_2), &
                 G_X3_Dn(:,iGF_h_2), G_X3_Up(:,iGF_h_2) )

        G_X3_Dn(:,iGF_h_2) = MAX( G_X3_Dn(:,iGF_h_2), SqrtTiny )
        G_X3_Up(:,iGF_h_2) = MAX( G_X3_Up(:,iGF_h_2), SqrtTiny )

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3(:,iGF_h_3), G_K(:,iGF_h_3), G_N_X3(:,iGF_h_3), &
                 G_X3_Dn(:,iGF_h_3), G_X3_Up(:,iGF_h_3) )

        G_X3_Dn(:,iGF_h_3) = MAX( G_X3_Dn(:,iGF_h_3), SqrtTiny )
        G_X3_Up(:,iGF_h_3) = MAX( G_X3_Up(:,iGF_h_3), SqrtTiny )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_h_1), G_X3_Dn(:,iGF_h_1), G_K(:,iGF_h_1), &
                 dh1dX3 )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_h_2), G_X3_Dn(:,iGF_h_2), G_K(:,iGF_h_2), &
                 dh2dX3 )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_h_3), G_X3_Dn(:,iGF_h_3), G_K(:,iGF_h_3), &
                 dh3dX3 )

        ! --- Shift vector derivatives wrt X3 ---

        ! --- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3(:,iGF_Beta_1), G_K(:,iGF_Beta_1), &
                 G_N_X3(:,iGF_Beta_1), &
                 G_X3_Dn(:,iGF_Beta_1), G_X3_Up(:,iGF_Beta_1) )

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3(:,iGF_Beta_2), G_K(:,iGF_Beta_2), &
                 G_N_X3(:,iGF_Beta_2), &
                 G_X3_Dn(:,iGF_Beta_2), G_X3_Up(:,iGF_Beta_2) )

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3(:,iGF_Beta_3), G_K(:,iGF_Beta_3), &
                 G_N_X3(:,iGF_Beta_3), &
                 G_X3_Dn(:,iGF_Beta_3), G_X3_Up(:,iGF_Beta_3) )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Beta_1), G_X3_Dn(:,iGF_Beta_1), &
                 G_K(:,iGF_Beta_1), &
                 db1dX3 )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Beta_2), G_X3_Dn(:,iGF_Beta_2), &
                 G_K(:,iGF_Beta_2), &
                 db2dX3 )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Beta_3), G_X3_Dn(:,iGF_Beta_3), &
                 G_K(:,iGF_Beta_3), &
                 db3dX3 )

        ! --- Lapse function derivative wrt X3 ---

        ! -- Interpolation ---

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3(:,iGF_Alpha), G_K(:,iGF_Alpha), G_N_X3(:,iGF_Alpha), &
                 G_X3_Dn(:,iGF_Alpha), G_X3_Up(:,iGF_Alpha) )

        G_X3_Dn(:,iGF_Alpha) = MAX( G_X3_Dn(:,iGF_Alpha), SqrtTiny )
        G_X3_Up(:,iGF_Alpha) = MAX( G_X3_Up(:,iGF_Alpha), SqrtTiny )

        ! --- Differentiation ---

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Alpha), G_X3_Dn(:,iGF_Alpha), &
                 G_K(:,iGF_Beta_3), &
                 dadX3 )

        ! --- Momentum density source term (S3) ---

        dU(:,iX1,iX2,iX3,iCF_S3) &
          = dU(:,iX1,iX2,iX3,iCF_S3)                                &
              + G_K(:,iGF_Alpha)                                    &
                  * ( ( Stress(:,1) * dh1dX3 ) / G_K(:,iGF_h_1)     &
                      + ( Stress(:,2) * dh2dX3 ) / G_K(:,iGF_h_2)   &
                      + ( Stress(:,3) * dh3dX3 ) / G_K(:,iGF_h_3) ) &
              + uCF_K(:,iCF_S1) * db1dX3 &
                  + uCF_K(:,iCF_S2) * db2dX3 &
                  + uCF_K(:,iCF_S3) * db3dX3 &
              - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx3

      END IF

      ! --- Compute divergence of grid volume ---

      ! --- Interpolate SqrtGm to faces (do this with scale factors instead?)---

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_SqrtGm), G_K(:,iGF_SqrtGm), &
               G_N_X1(:,iGF_SqrtGm), &
               G_X1_Dn(:,iGF_SqrtGm), G_X1_Up(:,iGF_SqrtGm) )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Alpha) * G_X1_Up(:,iGF_SqrtGm) &
                 * G_X1_Up(:,iGF_Beta_1), &
               G_X1_Dn(:,iGF_Alpha) * G_X1_Dn(:,iGF_SqrtGm) &
                 * G_X1_Dn(:,iGF_Beta_1), &
               G_K(:,iGF_Alpha) * G_K(:,iGF_SqrtGm) * G_K(:,iGF_Beta_1), &
               DivGridVolume, &
               Alpha_Option = One, Beta_Option = Zero )

      IF( nDimsX .GT. 1 )THEN

        CALL InterpolateToFace &
               ( nDOFX_X2, LX_X2_Up, LX_X2_Dn, &
                 G_P_X2(:,iGF_SqrtGm), G_K(:,iGF_SqrtGm), &
                 G_N_X2(:,iGF_SqrtGm), &
                 G_X2_Dn(:,iGF_SqrtGm), G_X2_Up(:,iGF_SqrtGm) )

        CALL ComputeDerivative &
               ( nDOFX_X2, dX2, LX_X2_Up, LX_X2_Dn, dLXdX2_q, WeightsX_X2, &
                 G_X2_Up(:,iGF_Alpha) * G_X2_Up(:,iGF_SqrtGm) &
                   * G_X2_Up(:,iGF_Beta_2), &
                 G_X2_Dn(:,iGF_Alpha) * G_X2_Dn(:,iGF_SqrtGm) &
                   * G_X2_Dn(:,iGF_Beta_2), &
                 G_K(:,iGF_Alpha) * G_K(:,iGF_SqrtGm) * G_K(:,iGF_Beta_2), &
                 DivGridVolume, &
                 Alpha_Option = One, Beta_Option = One )

      END IF

      IF( nDimsX .GT. 2 )THEN

        CALL InterpolateToFace &
               ( nDOFX_X3, LX_X3_Up, LX_X3_Dn, &
                 G_P_X3(:,iGF_SqrtGm), G_K(:,iGF_SqrtGm), &
                 G_N_X3(:,iGF_SqrtGm), &
                 G_X3_Dn(:,iGF_SqrtGm), G_X3_Up(:,iGF_SqrtGm) )

        CALL ComputeDerivative &
               ( nDOFX_X3, dX3, LX_X3_Up, LX_X3_Dn, dLXdX3_q, WeightsX_X3, &
                 G_X3_Up(:,iGF_Alpha) * G_X3_Up(:,iGF_SqrtGm) &
                   * G_X3_Up(:,iGF_Beta_3), &
                 G_X3_Dn(:,iGF_Alpha) * G_X3_Dn(:,iGF_SqrtGm) &
                   * G_X3_Dn(:,iGF_Beta_3), &
                 G_K(:,iGF_Alpha) * G_K(:,iGF_SqrtGm) * G_K(:,iGF_Beta_3), &
                 DivGridVolume, &
                 Alpha_Option = One, Beta_Option = One )

      END IF

      DivGridVolume &
        = One / ( G_K(:,iGF_Alpha) * G_K(:,iGF_SqrtGm) ) * DivGridVolume

      ! --- Compute energy increment ---

      ! --- Extrinsic curvature term ---

      EnergyDensitySourceTerms(:,iX1,iX2,iX3,1) &
        =   PressureTensor(:,1,1) * db1dX1 &
          + PressureTensor(:,1,2) * db2dX1 &
          + PressureTensor(:,1,3) * db3dX1 &
          + PressureTensor(:,2,1) * db1dX2 &
          + PressureTensor(:,2,2) * db2dX2 &
          + PressureTensor(:,2,3) * db3dX2 &
          + PressureTensor(:,3,1) * db1dX1 &
          + PressureTensor(:,3,2) * db2dX3 &
          + PressureTensor(:,3,3) * db3dX3

      ! --- Need to add Christoffel symbol term ---

      Christoffel3D_X1(:,1,1) = One / G_K(:,iGF_h_1) * dh1dX1
      Christoffel3D_X1(:,1,2) = One / G_K(:,iGF_h_1) * dh1dX2
      Christoffel3D_X1(:,1,3) = One / G_K(:,iGF_h_1) * dh1dX3
      Christoffel3D_X1(:,2,1) = Christoffel3D_X1(:,1,2)
      Christoffel3D_X1(:,2,2) = -G_K(:,iGF_h_2) / G_K(:,iGF_h_1)**2 * dh2dX1
      Christoffel3D_X1(:,2,3) = Zero
      Christoffel3D_X1(:,3,1) = Christoffel3D_X1(:,1,3)
      Christoffel3D_X1(:,3,2) = Christoffel3D_X1(:,2,3)
      Christoffel3D_X1(:,3,3) = -G_K(:,iGF_h_3) / G_K(:,iGF_h_1)**2 * dh3dX1

      Christoffel3D_X2(:,1,1) = -G_K(:,iGF_h_1) / G_K(:,iGF_h_2)**2 * dh1dX2
      Christoffel3D_X2(:,1,2) = One / G_K(:,iGF_h_2) * dh2dX1
      Christoffel3D_X2(:,1,3) = Zero
      Christoffel3D_X2(:,2,1) = Christoffel3D_X2(:,1,2)
      Christoffel3D_X2(:,2,2) = One / G_K(:,iGF_h_2) * dh2dX2
      Christoffel3D_X2(:,2,3) = One / G_K(:,iGF_h_2) * dh2dX3
      Christoffel3D_X2(:,3,1) = Christoffel3D_X2(:,1,3)
      Christoffel3D_X2(:,3,2) = Christoffel3D_X2(:,2,3)
      Christoffel3D_X2(:,3,3) = -G_K(:,iGF_h_3) / G_K(:,iGF_h_2)**2 * dh3dX2

      Christoffel3D_X3(:,1,1) = -G_K(:,iGF_h_1) / G_K(:,iGF_h_3)**2 * dh1dX3
      Christoffel3D_X3(:,1,2) = Zero
      Christoffel3D_X3(:,1,3) = One / G_K(:,iGF_h_3) * dh3dX1
      Christoffel3D_X3(:,2,1) = Christoffel3D_X3(:,1,2)
      Christoffel3D_X3(:,2,2) = -G_K(:,iGF_h_2) / G_K(:,iGF_h_3)**2 * dh2dX3
      Christoffel3D_X3(:,2,3) = One / G_K(:,iGF_h_3) * dh3dX2
      Christoffel3D_X3(:,3,1) = Christoffel3D_X3(:,1,3)
      Christoffel3D_X3(:,3,2) = Christoffel3D_X3(:,2,3)
      Christoffel3D_X3(:,3,3) = One / G_K(:,iGF_h_3) * dh3dX3

      Xij(:,1,1) &
        = G_K(:,iGF_Alpha)**( -2 ) &
            * ( G_K(:,iGF_Gm_dd_11) * db1dX1 &
                  + G_K(:,iGF_h_1) * (   G_K(:,iGF_Beta_1) * dh1dX1   &
                                       + G_K(:,iGF_Beta_2) * dh1dX2   &
                                       + G_K(:,iGF_Beta_3) * dh1dX3 ) &
                  - Third * G_K(:,iGF_Gm_dd_11) * DivGridVolume )

      Xij(:,1,2) &
        = Half * G_K(:,iGF_Alpha)**( -2 ) &
            * ( G_K(:,iGF_Gm_dd_11) * db1dX2 + G_K(:,iGF_Gm_dd_22) * db2dX1  )

      Xij(:,1,3) &
        = Half * G_K(:,iGF_Alpha)**( -2 ) &
            * ( G_K(:,iGF_Gm_dd_11) * db1dX3 + G_K(:,iGF_Gm_dd_33) * db3dX1  )

      Xij(:,2,1) = Xij(:,1,2)

      Xij(:,2,2) &
        = G_K(:,iGF_Alpha)**( -2 ) &
            * ( G_K(:,iGF_Gm_dd_22) * db2dX2 &
                  + G_K(:,iGF_h_2) * (   G_K(:,iGF_Beta_1) * dh2dX1   &
                                       + G_K(:,iGF_Beta_2) * dh2dX2   &
                                       + G_K(:,iGF_Beta_3) * dh2dX3 ) &
                  - Third * G_K(:,iGF_Gm_dd_22) * DivGridVolume )

      Xij(:,2,3) &
        = Half * G_K(:,iGF_Alpha)**( -2 ) &
            * ( G_K(:,iGF_Gm_dd_22) * db2dX3 + G_K(:,iGF_Gm_dd_33) * db3dX2  )

      Xij(:,3,1) = Xij(:,1,3)

      Xij(:,3,2) = Xij(:,2,3)

      Xij(:,3,3) &
        = G_K(:,iGF_Alpha)**( -2 ) &
            * ( G_K(:,iGF_Gm_dd_33) * db3dX3 &
                  + G_K(:,iGF_h_3) * (   G_K(:,iGF_Beta_1) * dh3dX1   &
                                       + G_K(:,iGF_Beta_2) * dh3dX2   &
                                       + G_K(:,iGF_Beta_3) * dh3dX3 ) &
                  - Third * G_K(:,iGF_Gm_dd_33) * DivGridVolume )

      DO iDim = 1, 3

        DO jDim = 1, 3

          Christoffel_X1(:,iDim,jDim) &
            = G_K(:,iGF_Beta_1) * Xij(:,iDim,jDim) &
                + Christoffel3D_X1(:,iDim,jDim)

          Christoffel_X2(:,iDim,jDim) &
            = G_K(:,iGF_Beta_2) * Xij(:,iDim,jDim) &
                + Christoffel3D_X2(:,iDim,jDim)

          Christoffel_X3(:,iDim,jDim) &
            = G_K(:,iGF_Beta_3) * Xij(:,iDim,jDim) &
                + Christoffel3D_X3(:,iDim,jDim)

        END DO

      END DO

      EnergyDensitySourceTerms(:,iX1,iX2,iX3,2) = Zero

      DO iNodeX = 1, nDOFX

        PressureTensor(iNodeX,:,1) &
          = G_K(iNodeX,iGF_Gm_dd_11) * PressureTensor(iNodeX,:,1)
        PressureTensor(iNodeX,:,2) &
          = G_K(iNodeX,iGF_Gm_dd_22) * PressureTensor(iNodeX,:,2)
        PressureTensor(iNodeX,:,3) &
          = G_K(iNodeX,iGF_Gm_dd_33) * PressureTensor(iNodeX,:,3)

      END DO

      ! Get gradient of conformal factor on upper face

      CALL InterpolateToFace &
             ( nDOFX_X1, LX_X1_Up, LX_X1_Dn, &
               G_P_X1(:,iGF_Psi), G_K(:,iGF_Psi), G_N_X1(:,iGF_Psi), &
               G_X1_Dn(:,iGF_Psi), G_X1_Up(:,iGF_Psi) )

      CALL ComputeDerivative &
             ( nDOFX_X1, dX1, LX_X1_Up, LX_X1_Dn, dLXdX1_q, WeightsX_X1, &
               G_X1_Up(:,iGF_Psi), &
               G_X1_Dn(:,iGF_Psi), &
               G_K(:,iGF_Psi),     &
               GradPsi,            &
               Alpha_Option = One, Beta_Option = Zero )

      GradPsiF = Zero

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GradPsi, 1, Zero, GradPsiF(1:nDOFX_X1), 1 )

      DO iDim = 1, 3

        EnergyDensitySourceTerms(:,iX1,iX2,iX3,2) &
        = EnergyDensitySourceTerms(:,iX1,iX2,iX3,2) &
            + PressureTensor(:,iDim,1) &
                * (   Christoffel_X1(:,iDim,1) * G_K(:,iGF_Beta_1)   &
                    + Christoffel_X1(:,iDim,2) * G_K(:,iGF_Beta_2)   &
                    + Christoffel_X1(:,iDim,3) * G_K(:,iGF_Beta_3) ) &
            + PressureTensor(:,iDim,2) &
                * (   Christoffel_X2(:,iDim,1) * G_K(:,iGF_Beta_1)   &
                    + Christoffel_X2(:,iDim,2) * G_K(:,iGF_Beta_2)   &
                    + Christoffel_X2(:,iDim,3) * G_K(:,iGF_Beta_3) ) &
            + PressureTensor(:,iDim,3) &
                * (   Christoffel_X3(:,iDim,1) * G_K(:,iGF_Beta_1)   &
                    + Christoffel_X3(:,iDim,2) * G_K(:,iGF_Beta_2)   &
                    + Christoffel_X3(:,iDim,3) * G_K(:,iGF_Beta_3) )

      END DO

      PressureTensorTrace &
        =   uCF_K(:,iCF_S1) * uPF_K(:,iPF_V1) &
          + uCF_K(:,iCF_S2) * uPF_K(:,iPF_V2) &
          + uCF_K(:,iCF_S3) * uPF_K(:,iPF_V3) + Three * P_K

      EnergyDensitySourceTerms(:,iX1,iX2,iX3,3) &
        = -Third * PressureTensorTrace * DivGridVolume

      EnergyDensitySourceTerms(:,iX1,iX2,iX3,4) &
        = -(   uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11) * dadx1 &
             + uCF_K(:,iCF_S2) / G_K(:,iGF_Gm_dd_22) * dadx2 &
             + uCF_K(:,iCF_S3) / G_K(:,iGF_Gm_dd_33) * dadx3 )

      EnergyDensitySourceTerms(:,iX1,iX2,iX3,5) &
        = -uCF_K(:,iCF_E) * DivGridVolume

      DO iNodeX = 1, nNodesX(1)

        X1(iNodeX) = NodeCoordinate( MeshX(1), iX1, iNodeX )

      END DO

      EnergyDensitySourceTerms(:,iX1,iX2,iX3,6) &
        = G_K(:,iGF_Psi)**12 / G_K(:,iGF_Alpha)**2 * 2.0_DP / 3.0_DP &
            * ( db1dX1**2 - 2.0_DP / X1 * db1dX1 * G_K(:,iGF_Beta_1) &
                  + 1.0_DP / X1**2 * G_K(:,iGF_Beta_1)**2 )

      EnergyDensitySourceTerms(:,iX1,iX2,iX3,7) &
        = GradPsiF

      ! --- Add to increments ---

      dU(:,iX1,iX2,iX3,iCF_E) &
        = dU(:,iX1,iX2,iX3,iCF_E) &
            + EnergyDensitySourceTerms(:,iX1,iX2,iX3,1) &
            + EnergyDensitySourceTerms(:,iX1,iX2,iX3,2) &
            + EnergyDensitySourceTerms(:,iX1,iX2,iX3,3) &
            + EnergyDensitySourceTerms(:,iX1,iX2,iX3,4)

      DO iCF = 1, nCF

        dU(:,iX1,iX2,iX3,iCF) &
          = dU(:,iX1,iX2,iX3,iCF) &
              - U(:,iX1,iX2,iX3,iCF) * DivGridVolume

      END DO

    END DO
    END DO
    END DO

    IF( WriteSourceTerms )THEN

      CALL WriteSourceTermDiagnosticsHDF( Time, EnergyDensitySourceTerms )

    END IF

  END SUBROUTINE ComputeIncrement_Geometry_Relativistic


  SUBROUTINE ComputeIncrement_Gravity &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

   INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    CALL TimersStart_Euler( Timer_Euler_Gravity )

#if defined HYDRO_RELATIVISTIC

    CALL ComputeIncrement_Gravity_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#else

    CALL ComputeIncrement_Gravity_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#endif

    CALL TimersStop_Euler( Timer_Euler_Gravity )

  END SUBROUTINE ComputeIncrement_Gravity


  SUBROUTINE ComputeIncrement_Gravity_NonRelativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: Phi_P_X1(nDOFX)
    REAL(DP) :: Phi_K   (nDOFX)
    REAL(DP) :: Phi_N_X1(nDOFX)
    REAL(DP) :: dPhidX1(nDOFX)
    REAL(DP) :: Phi_X1_Dn(nDOFX_X1)
    REAL(DP) :: Phi_X1_Up(nDOFX_X1)
    REAL(DP) :: uCF_K(nDOFX,nCF)

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      DO iCF = 1, nCF

        uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

      END DO

      Phi_P_X1(:) = G(:,iX1-1,iX2,iX3,iGF_Phi_N)
      Phi_K   (:) = G(:,iX1,  iX2,iX3,iGF_Phi_N)
      Phi_N_X1(:) = G(:,iX1+1,iX2,iX3,iGF_Phi_N)

      ! --- Derivative of Potential wrt X1 ---

      ! --- Face States (Average of Left and Right States) ---

      ! --- Face at X1_L ---

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               Phi_P_X1(:), 1, Zero, Phi_X1_Dn(:), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               Phi_K   (:), 1, Half, Phi_X1_Dn(:), 1 )

      ! --- Face at X1_H ---

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               Phi_K   (:), 1, Zero, Phi_X1_Up(:), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               Phi_N_X1(:), 1, Half, Phi_X1_Up(:), 1 )

      ! --- dPhidX1 ---

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1(:) * Phi_X1_Up(:), 1, Zero, dPhidX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1(:) * Phi_X1_Dn(:), 1,  One, dPhidX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q (:) * Phi_K    (:), 1,  One, dPhidX1, 1 )

      dPhidX1 = dPhidX1 / ( WeightsX_q(:) * dX1 )

      ! --- Increments ---

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1) &
            - uCF_K(:,iCF_D) * dPhidX1(:)

      dU(:,iX1,iX2,iX3,iCF_E) &
        = dU(:,iX1,iX2,iX3,iCF_E) &
            - uCF_K(:,iCF_S1) * dPhidX1(:)

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_Gravity_NonRelativistic


  SUBROUTINE ComputeIncrement_Gravity_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)
  END SUBROUTINE ComputeIncrement_Gravity_Relativistic


  SUBROUTINE InterpolateToFace &
    ( nDOFX_X, LX_Up, LX_Dn, &
      InterpolantP, InterpolantK, InterpolantN, &
      Answer_Dn, Answer_Up )

    INTEGER,  INTENT(in)  :: nDOFX_X
    REAL(DP), INTENT(in)  :: LX_Up(:,:)
    REAL(DP), INTENT(in)  :: LX_Dn(:,:)
    REAL(DP), INTENT(in)  :: InterpolantP(:), InterpolantK(:), InterpolantN(:)
    REAL(DP), INTENT(out) :: Answer_Dn(:), Answer_Up(:)

    CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

    CALL DGEMV &
           ( 'N', nDOFX_X, nDOFX, One,  LX_Up, nDOFX_X, &
             InterpolantP, 1, Zero, Answer_Dn, 1 )

    CALL DGEMV &
           ( 'N', nDOFX_X, nDOFX, Half, LX_Dn, nDOFX_X, &
             InterpolantK, 1, Half, Answer_Dn, 1 )

    CALL DGEMV &
           ( 'N', nDOFX_X, nDOFX, One,  LX_Up, nDOFX_X, &
             InterpolantK, 1, Zero, Answer_Up, 1 )

    CALL DGEMV &
           ( 'N', nDOFX_X, nDOFX, Half, LX_Dn, nDOFX_X, &
             InterpolantN, 1, Half, Answer_Up, 1 )

    CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

  END SUBROUTINE InterpolateToFace


  SUBROUTINE ComputeDerivative &
    ( nDOFX_X, dX, LX_Up, LX_Dn, dLX, WeightsX_X, &
      UpperFaceValues, LowerFaceValues, VolumeValues, &
      Answer, Alpha_Option, Beta_Option )

    INTEGER,  INTENT(in)  :: nDOFX_X
    REAL(DP), INTENT(in)  :: dX
    REAL(DP), INTENT(in)  :: LX_Up(:,:)
    REAL(DP), INTENT(in)  :: LX_Dn(:,:)
    REAL(DP), INTENT(in)  :: dLX  (:,:)
    REAL(DP), INTENT(in)  :: WeightsX_X(:)
    REAL(DP), INTENT(in)  :: UpperFaceValues(:)
    REAL(DP), INTENT(in)  :: LowerFaceValues(:)
    REAL(DP), INTENT(in)  :: VolumeValues   (:)
    REAL(DP), INTENT(out) :: Answer(:)
    REAL(DP), INTENT(in), OPTIONAL :: Alpha_Option, Beta_Option

    REAL(DP) :: Alpha, Beta

    Alpha = One
    IF( PRESENT( Alpha_Option ) ) Alpha = Alpha_Option

    Beta = Zero
    IF( PRESENT( Beta_Option ) ) Beta = Beta_Option

    Answer = Zero

    CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

    CALL DGEMV( 'T', nDOFX_X, nDOFX, + Alpha, LX_Up, nDOFX_X, &
                WeightsX_X * UpperFaceValues, 1, Beta, Answer, 1 )

    CALL DGEMV( 'T', nDOFX_X, nDOFX, - Alpha, LX_Dn, nDOFX_X, &
                WeightsX_X * LowerFaceValues, 1, One , Answer, 1 )

    CALL DGEMV( 'T', nDOFX,    nDOFX, - Alpha, dLX, nDOFX,    &
                WeightsX_q  * VolumeValues  , 1, One , Answer, 1 )

    CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

    Answer = Answer / ( WeightsX_q * dX )

  END SUBROUTINE ComputeDerivative


END MODULE Euler_dgDiscretizationModule
