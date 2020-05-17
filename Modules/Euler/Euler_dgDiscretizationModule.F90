MODULE Euler_dgDiscretizationModule

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One, Pi, TwoPi
  USE TimersModule_Euler,  ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_dgDiscretization, &
    Timer_Euler_Divergence, &
    Timer_Euler_Geometry, &
    Timer_Euler_Gravity, &
    Timer_Euler_SurfaceTerm, &
    Timer_Euler_NumericalFlux, &
    Timer_Euler_VolumeTerm, &
    Timer_Euler_Increment, &
    Timer_Euler_ComputePrimitive, &
    Timer_Euler_MatrixVectorMultiply
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, WeightsX_X1, &
    nDOFX_X2, WeightsX_X2, &
    nDOFX_X3, WeightsX_X3, &
    WeightsX_q,            &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, LX_X1_Dn, LX_X1_Up, &
    dLXdX2_q, LX_X2_Dn, LX_X2_Up, &
    dLXdX3_q, LX_X3_Dn, LX_X3_Up
  USE GeometryFieldsModule, ONLY: &
    nGF,                                      &
    iGF_h_1,      iGF_h_2,      iGF_h_3,      &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_SqrtGm,                               &
    iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3, &
    iGF_Phi_N,                                &
    CoordinateSystem
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_P, &
    iDF_Sh_X1, iDF_Sh_X2, iDF_Sh_X3
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

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_Euler_DG_Explicit


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

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectShocks_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    CALL TimersStart_Euler( Timer_Euler_Divergence )

    CALL ComputeIncrement_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    CALL ComputeIncrement_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    CALL ComputeIncrement_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU )

    CALL TimersStop_Euler( Timer_Euler_Divergence )

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

    CALL TimersStart_Euler( Timer_Euler_Geometry )

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_Geometry )

    CALL TimersStart_Euler( Timer_Euler_Gravity )

    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    CALL TimersStop_Euler( Timer_Euler_Gravity )

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

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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
                MAXVAL( D(:,iX1-1,iX2,iX3,iDF_Sh_X1) ), &
                MAXVAL( D(:,iX1  ,iX2,iX3,iDF_Sh_X1) ) )

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

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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
                MAXVAL( D(:,iX1,iX2-1,iX3,iDF_Sh_X2) ), &
                MAXVAL( D(:,iX1,iX2  ,iX3,iDF_Sh_X2) ) )

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

        CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

        CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStart_Euler( Timer_Euler_ComputePrimitive )

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

      CALL TimersStop_Euler( Timer_Euler_ComputePrimitive )

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
                MAXVAL( D(:,iX1,iX2,iX3-1,iDF_Sh_X3) ), &
                MAXVAL( D(:,iX1,iX2,iX3  ,iDF_Sh_X3) ) )

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

#if defined HYDRO_NONRELATIVISTIC

    CALL ComputeIncrement_Geometry_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputeIncrement_Geometry_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#endif

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

    ! --- This subroutine currently assumes a stationary spacetime ---

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX
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

      ! --- Scale factor face states (average of left and right states) ---

      DO iGF = iGF_h_1, iGF_h_3

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P_X1(:,iGF), 1, Zero, G_X1_Dn(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K   (:,iGF), 1, Half, G_X1_Dn(:,iGF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_X1_Dn(:,iGF) = MAX( G_X1_Dn(:,iGF), SqrtTiny )

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_K   (:,iGF), 1, Zero, G_X1_Up(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_N_X1(:,iGF), 1, Half, G_X1_Up(:,iGF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_X1_Up(:,iGF) = MAX( G_X1_Up(:,iGF), SqrtTiny )

      END DO

      ! --- dh1dx1 ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_h_1), 1, Zero, dh1dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_1), 1, One,  dh1dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_1), 1, One,  dh1dX1, 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      dh1dx1 = dh1dx1 / ( WeightsX_q * dX1 )

      ! --- dh2dx1 ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_2), 1, One,  dh2dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_2), 1, One,  dh2dX1, 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      dh2dx1 = dh2dx1 / ( WeightsX_q * dX1 )

      ! --- dh3dx1 ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_3), 1, One,  dh3dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_3), 1, One,  dh3dX1, 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      dh3dx1 = dh3dx1 / ( WeightsX_q * dX1 )

      ! --- Shift vector derivative wrt X1 ---

      ! --- Shift vector face states (average of left and right states) ---

      DO iGF = iGF_Beta_1, iGF_Beta_3

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P_X1(:,iGF), 1, Zero, G_X1_Dn(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K   (:,iGF), 1, Half, G_X1_Dn(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_K   (:,iGF), 1, Zero, G_X1_Up(:,iGF), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_N_X1(:,iGF), 1, Half, G_X1_Up(:,iGF), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      END DO

      ! --- db1dx1 ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_Beta_1), 1, Zero, db1dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_Beta_1), 1, One,  db1dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_Beta_1), 1, One,  db1dX1, 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      db1dx1 = db1dx1 / ( WeightsX_q * dX1 )

      ! --- db2dx1 ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_Beta_2), 1, Zero, db2dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_Beta_2), 1, One,  db2dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_Beta_2), 1, One,  db2dX1, 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      db2dx1 = db2dx1 / ( WeightsX_q * dX1 )

      ! --- db3dx1 ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_Beta_3), 1, Zero, db3dX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_Beta_3), 1, One,  db3dX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_Beta_3), 1, One,  db3dX1, 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      db3dx1 = db3dx1 / ( WeightsX_q * dX1 )

      ! --- Lapse function derivative wrt X1 ---

      ! --- Lapse function face states (average of left and right states) ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_P_X1(:,iGF_Alpha), 1, Zero, G_X1_Dn(:,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_K   (:,iGF_Alpha), 1, Half, G_X1_Dn(:,iGF_Alpha), 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      G_X1_Dn(:,iGF_Alpha) = MAX( G_X1_Dn(:,iGF_Alpha), SqrtTiny )

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_K   (:,iGF_Alpha), 1, Zero, G_X1_Up(:,iGF_Alpha), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_N_X1(:,iGF_Alpha), 1, Half, G_X1_Up(:,iGF_Alpha), 1 )
      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      G_X1_Up(:,iGF_Alpha) = MAX( G_X1_Up(:,iGF_Alpha), SqrtTiny )

      ! --- dadx1 ---

      CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_Alpha), 1, Zero, dadX1, 1 )

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_Alpha), 1, One,  dadX1, 1 )

      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_Alpha), 1, One,  dadX1, 1 )

      CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

      dadx1 = dadx1 / ( WeightsX_q * dX1 )

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

        ! --- Scale factor face states (average of left and right states) ---

        DO iGF = iGF_h_1, iGF_h_3

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_P_X2(:,iGF), 1, Zero, G_X2_Dn(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_K   (:,iGF), 1, Half, G_X2_Dn(:,iGF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

          G_X2_Dn(:,iGF) = MAX( G_X2_Dn(:,iGF), SqrtTiny )

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_K   (:,iGF), 1, Zero, G_X2_Up(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_N_X2(:,iGF), 1, Half, G_X2_Up(:,iGF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

          G_X2_Up(:,iGF) = MAX( G_X2_Up(:,iGF), SqrtTiny )

        END DO

        ! --- dh1dx2 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_h_1), 1, Zero, dh1dX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_h_1), 1, One,  dh1dX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_1), 1, One,  dh1dX2, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        dh1dx2 = dh1dx2 / ( WeightsX_q * dX2 )

        ! --- dh2dx2 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_h_2), 1, Zero, dh2dX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_h_2), 1, One,  dh2dX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_2), 1, One,  dh2dX2, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        dh2dx2 = dh2dx2 / ( WeightsX_q * dX2 )

        ! --- dh3dx2 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_h_3), 1, Zero, dh3dX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_h_3), 1, One,  dh3dX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_3), 1, One,  dh3dX2, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        dh3dx2 = dh3dx2 / ( WeightsX_q * dX2 )

        ! --- Shift vector derivatives wrt X2 ---

        ! --- Shift vector face states (average of left and right states) ---

        DO iGF = iGF_Beta_1, iGF_Beta_3

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_P_X2(:,iGF), 1, Zero, G_X2_Dn(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_K   (:,iGF), 1, Half, G_X2_Dn(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_K   (:,iGF), 1, Zero, G_X2_Up(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_N_X2(:,iGF), 1, Half, G_X2_Up(:,iGF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

        ! --- db1dx2 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_Beta_1), 1, Zero, db1dX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_Beta_1), 1, One,  db1dX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Beta_1), 1, One,  db1dX2, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        db1dx2 = db1dx2 / ( WeightsX_q * dX2 )

        ! --- db2dx2 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_Beta_2), 1, Zero, db2dX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_Beta_2), 1, One,  db2dX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Beta_2), 1, One,  db2dX2, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        db2dx2 = db2dx2 / ( WeightsX_q * dX2 )

        ! --- db3dx2 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_Beta_3), 1, Zero, db3dX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_Beta_3), 1, One,  db3dX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Beta_3), 1, One,  db3dX2, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        db3dx2 = db3dx2 / ( WeightsX_q * dX2 )

        ! --- Lapse function derivative wrt X2 ---

        ! --- Face states (average of left and right states) ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                 G_P_X2(:,iGF_Alpha), 1, Zero, G_X2_Dn(:,iGF_Alpha), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                 G_K   (:,iGF_Alpha), 1, Half, G_X2_Dn(:,iGF_Alpha), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_X2_Dn(:,iGF_Alpha) = MAX( G_X2_Dn(:,iGF_Alpha), SqrtTiny )

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                 G_K   (:,iGF_Alpha), 1, Zero, G_X2_Up(:,iGF_Alpha), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                 G_N_X2(:,iGF_Alpha), 1, Half, G_X2_Up(:,iGF_Alpha), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_X2_Up(:,iGF_Alpha) = MAX( G_X2_Up(:,iGF_Alpha), SqrtTiny )

        ! --- dadx2 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_Alpha), 1, Zero, dadX2, 1 )

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_Alpha), 1, One,  dadX2, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Alpha), 1, One,  dadX2, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        dadx2 = dadx2 / ( WeightsX_q * dX2 )

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

        ! --- Scale factor face states (average of left and right states) ---

        DO iGF = iGF_h_1, iGF_h_3

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                   G_P_X3(:,iGF), 1, Zero, G_X3_Dn(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                   G_K   (:,iGF), 1, Half, G_X3_Dn(:,iGF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

          G_X3_Dn(:,iGF) = MAX( G_X3_Dn(:,iGF), SqrtTiny )

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                   G_K   (:,iGF), 1, Zero, G_X3_Up(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                   G_N_X3(:,iGF), 1, Half, G_X3_Up(:,iGF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

          G_X3_Up(:,iGF) = MAX( G_X3_Up(:,iGF), SqrtTiny )

        END DO

        ! --- dh1dx3 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_h_1), 1, Zero, dh1dX3, 1 )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_h_1), 1, One,  dh1dX3, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_1), 1, One,  dh1dX3, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        dh1dx3 = dh1dx3 / ( WeightsX_q * dX3 )

        ! --- dh2dx3 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_h_2), 1, Zero, dh2dX3, 1 )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_h_2), 1, One,  dh2dX3, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_2), 1, One,  dh2dX3, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        dh2dx3 = dh2dx3 / ( WeightsX_q * dX3 )

        ! --- dh3dx3 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_h_3), 1, Zero, dh3dX3, 1 )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                  WeightsX_X3 * G_X3_Dn(:,iGF_h_3), 1, One,  dh3dX3, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_3), 1, One,  dh3dX3, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        dh3dx3 = dh3dx3 / ( WeightsX_q * dX3 )

        ! --- Shift vector derivatives wrt X3 ---

        ! --- Shift vector face states (average of left and right states) ---

        DO iGF = iGF_Beta_1, iGF_Beta_3

          CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                   G_P_X3(:,iGF), 1, Zero, G_X3_Dn(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                   G_K   (:,iGF), 1, Half, G_X3_Dn(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                   G_K   (:,iGF), 1, Zero, G_X3_Up(:,iGF), 1 )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                   G_N_X3(:,iGF), 1, Half, G_X3_Up(:,iGF), 1 )

          CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        END DO

        ! --- db1dx3 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_Beta_1), 1, Zero, db1dX3, 1 )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_Beta_1), 1, One,  db1dX3, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Beta_1), 1, One,  db1dX3, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        db1dx3 = db1dx3 / ( WeightsX_q * dX3 )

        ! --- db2dx3 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_Beta_2), 1, Zero, db2dX3, 1 )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_Beta_2), 1, One,  db2dX3, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Beta_2), 1, One,  db2dX3, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        db2dx3 = db2dx3 / ( WeightsX_q * dX3 )

        ! --- db3dx3 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_Beta_3), 1, Zero, db3dX3, 1 )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_Beta_3), 1, One,  db3dX3, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Beta_3), 1, One,  db3dX3, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        db3dx3 = db3dx3 / ( WeightsX_q * dX3 )

        ! --- Lapse function derivative wrt X3 ---

        ! --- Lapse function face states (average of left and right states) ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                 G_P_X3(:,iGF_Alpha), 1, Zero, G_X3_Dn(:,iGF_Alpha), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                 G_K   (:,iGF_Alpha), 1, Half, G_X3_Dn(:,iGF_Alpha), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_X3_Dn(:,iGF_Alpha) = MAX( G_X3_Dn(:,iGF_Alpha), SqrtTiny )

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                 G_K   (:,iGF_Alpha), 1, Zero, G_X3_Up(:,iGF_Alpha), 1 )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                 G_N_X3(:,iGF_Alpha), 1, Half, G_X3_Up(:,iGF_Alpha), 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        G_X3_Up(:,iGF_Alpha) = MAX( G_X3_Up(:,iGF_Alpha), SqrtTiny )

        ! --- dadx3 ---

        CALL TimersStart_Euler( Timer_Euler_MatrixVectorMultiply )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_Alpha), 1, Zero, dadX3, 1 )

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_Alpha), 1, One,  dadX3, 1 )

        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Alpha), 1, One,  dadX3, 1 )

        CALL TimersStop_Euler( Timer_Euler_MatrixVectorMultiply )

        dadx3 = dadx3 / ( WeightsX_q * dX3 )

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

      ! --- Compute energy increment (missing time-dependent metric term) ---

      dU(:,iX1,iX2,iX3,iCF_E) &
        = dU(:,iX1,iX2,iX3,iCF_E) &
            + Stress(:,1) / G_K(:,iGF_h_1) &
                * ( G_K(:,iGF_h_1) * db1dX1 &
                    + G_K(:,iGF_Beta_1) * dh1dX1 &
                    + G_K(:,iGF_Beta_2) * dh1dX2 &
                    + G_K(:,iGF_Beta_3) * dh1dX3 ) &
                    + uCF_K(:,iCF_S2) * uPF_K(:,iPF_V1) * db2dX1 &
                    + uCF_K(:,iCF_S3) * uPF_K(:,iPF_V1) * db3dX1 &
            + Stress(:,2) / G_K(:,iGF_h_2) &
                * ( G_K(:,iGF_h_2) * db2dX2 &
                    + G_K(:,iGF_Beta_1) * dh2dX1 &
                    + G_K(:,iGF_Beta_2) * dh2dX2 &
                    + G_K(:,iGF_Beta_3) * dh2dX3 &
                    + uCF_K(:,iCF_S1) * uPF_K(:,iPF_V2) * db1dX2 &
                    + uCF_K(:,iCF_S3) * uPF_K(:,iPF_V2) * db3dX2 ) &
            + Stress(:,3) / G_K(:,iGF_h_3) &
                * ( G_K(:,iGF_h_3) * db3dX3 &
                    + G_K(:,iGF_Beta_1) * dh3dX1 &
                    + G_K(:,iGF_Beta_2) * dh3dX2 &
                    + G_K(:,iGF_Beta_3) * dh3dX3 &
                    + uCF_K(:,iCF_S1) * uPF_K(:,iPF_V3) * db1dX3 &
                    + uCF_K(:,iCF_S2) * uPF_K(:,iPF_V3) * db2dX3 ) &
            - (   uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11) * dadx1 &
                + uCF_K(:,iCF_S2) / G_K(:,iGF_Gm_dd_22) * dadx2 &
                + uCF_K(:,iCF_S3) / G_K(:,iGF_Gm_dd_33) * dadx3 )

    END DO
    END DO
    END DO

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

#if defined HYDRO_NONRELATIVISTIC

    CALL ComputeIncrement_Gravity_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputeIncrement_Gravity_Relativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#endif

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


END MODULE Euler_dgDiscretizationModule
