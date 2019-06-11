MODULE Euler_dgDiscretizationModule

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One, Pi, TwoPi
  USE TimersModule_Euler
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
    nAF, iAF_P
  USE Euler_BoundaryConditionsModule, ONLY: &
    Euler_ApplyBoundaryConditions
  USE Euler_UtilitiesModule, ONLY: &
    Euler_ComputePrimitive,      &
    Euler_Eigenvalues,           &
    Euler_AlphaMiddle,           &
    Euler_Flux_X1,               &
    Euler_Flux_X2,               &
    Euler_StressTensor_Diagonal, &
    Euler_NumericalFlux_HLL,     &
    Euler_NumericalFlux_X1_HLLC, &
    Euler_NumericalFlux_X2_HLLC, &
    Euler_NumericalFlux_X3_HLLC
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: Euler_ComputeIncrement_DG_Explicit

  LOGICAL :: TimeIt = .FALSE.


CONTAINS


  SUBROUTINE Euler_ComputeIncrement_DG_Explicit &
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

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: dX1, dX2, dX3
    LOGICAL  :: SuppressBC

    IF( TimeIt ) CALL InitializeTimers_Euler
    IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_Inc )

    dU = Zero

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL Euler_ApplyBoundaryConditions &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_Div_X1 )
    CALL ComputeIncrement_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )
    IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_Div_X1 )

    IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_Div_X2 )
    CALL ComputeIncrement_Divergence_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )
    IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_Div_X2 )

    IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_Div_X3 )
    CALL ComputeIncrement_Divergence_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )
    IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_Div_X3 )

    ! --- Multiply Inverse Mass Matrix ---

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

    IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_Geom )
    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )
    IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_Geom )

    IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_Grav )
    CALL ComputeIncrement_Gravity &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )
    IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_Grav )

    IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_Inc )
    IF( TimeIt ) CALL FinalizeTimers_Euler
    IF( TimeIt ) STOP 'Euler_ComputeIncrement_DG_Explicit'

  END SUBROUTINE Euler_ComputeIncrement_DG_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
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

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        IF( TimeIt) CALL TimersStart_Euler( Timer_Euler_CompPrim )
        CALL Euler_ComputePrimitive &
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
        IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_CompPrim )

        CALL ComputePressureFromPrimitive &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

        DO iNodeX = 1, nDOFX

          Flux_X1_q(iNodeX,:) &
            = Euler_Flux_X1 &
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
          IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
          CALL DGEMV &
                 ( 'T', nDOFX, nDOFX, One, dLXdX1_q, nDOFX, &
                   Flux_X1_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )
          IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

        END DO

      END IF

      !---------------------
      ! --- Surface Term ---
      !---------------------

      ! --- Interpolate Fluid Fields ---

      DO iCF = 1, nCF

        IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
                 uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                 uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )
        IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

      END DO

      ! --- Interpolate Geometry Fields ---

      ! --- Face States (Average of Left and Right States) ---

      G_F = Zero

      DO iGF = iGF_h_1, iGF_h_3

        IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                 G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                 G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )
        IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

        G_F(:,iGF) = MAX( G_F(:,iGF), SqrtTiny )

      END DO

      CALL ComputeGeometryX_FromScaleFactors( G_F )

      ! --- Lapse Function ---

      IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               G_P(:,iGF_Alpha), 1, Zero, G_F(:,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               G_K(:,iGF_Alpha), 1, Half, G_F(:,iGF_Alpha), 1 )
      IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

      G_F(:,iGF_Alpha) = MAX( G_F(:,iGF_Alpha), SqrtTiny )

      ! --- Left State Primitive, etc. ---

      IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_CompPrim )
      CALL Euler_ComputePrimitive &
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
      IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_CompPrim )

      CALL ComputePressureFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L )

      DO iNodeX_X1 = 1, nDOFX_X1

        EigVals_L(iNodeX_X1,:) &
          = Euler_Eigenvalues &
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
          = Euler_Flux_X1 &
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

      IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_CompPrim )
      CALL Euler_ComputePrimitive &
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
      IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_CompPrim )

      CALL ComputePressureFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), P_R  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R )

      DO iNodeX_X1 = 1, nDOFX_X1

        EigVals_R(iNodeX_X1,:) &
          = Euler_Eigenvalues &
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
          = Euler_Flux_X1 &
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
          = Euler_AlphaMiddle &
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
                G_F      (iNodeX_X1,iGF_Alpha),    &
                G_F      (iNodeX_X1,iGF_Beta_1),   &
                AlphaPls, AlphaMns )

        NumericalFlux(iNodeX_X1,:) &
!          = Euler_NumericalFlux_X1_HLLC &
          = Euler_NumericalFlux_HLL &
              ( uCF_L    (iNodeX_X1,:),            &
                uCF_R    (iNodeX_X1,:),            &
                Flux_X1_L(iNodeX_X1,:),            &
                Flux_X1_R(iNodeX_X1,:),            &
                G_F      (iNodeX_X1,iGF_Gm_dd_11), &
                uPF_L    (iNodeX_X1,iPF_V1),       &
                uPF_R    (iNodeX_X1,iPF_V1),       &
                P_L      (iNodeX_X1),              &
                P_R      (iNodeX_X1),              &
                G_F      (iNodeX_X1,iGF_Alpha),    &
                G_F      (iNodeX_X1,iGF_Beta_1),   &
                AlphaPls, AlphaMns, AlphaMdl )

      END DO

      DO iCF = 1, nCF

        NumericalFlux(:,iCF) &
          = dX2 * dX3 * WeightsX_X1 * G_F(:,iGF_Alpha) * G_F(:,iGF_SqrtGm) &
              * NumericalFlux(:,iCF)

      END DO

      ! --- Contribution to This Element ---

      IF( iX1 .LT. iX_E0(1) + 1 )THEN

        DO iCF = 1, nCF

          IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
          CALL DGEMV &
                 ( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Dn, nDOFX_X1, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )
          IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

        END DO

      END IF

      ! --- Contribution to Previous Element ---

      IF( iX1 .GT. iX_B0(1) )THEN

        DO iCF = 1, nCF

          IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
          CALL DGEMV &
                 ( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Up, nDOFX_X1, &
                   NumericalFlux(:,iCF), 1, One, dU(:,iX1-1,iX2,iX3,iCF), 1 )
          IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

        END DO

      END IF

    END DO
    END DO
    END DO
    !$OMP END PARALLEL DO

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

      IF( iX2 .LT. iX_E0(2) + 1 )THEN

        IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_CompPrim )
        CALL Euler_ComputePrimitive &
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
        IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_CompPrim )

        CALL ComputePressureFromPrimitive &
               ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

        DO iNodeX = 1, nDOFX

          Flux_X2_q(iNodeX,:) &
            = Euler_Flux_X2 &
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

          IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
          CALL DGEMV &
                 ( 'T', nDOFX, nDOFX, One, dLXdX2_q, nDOFX, &
                   Flux_X2_q(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )
          IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

        END DO

      END IF

      !---------------------
      ! --- Surface Term ---
      !---------------------

      ! --- Interpolate Fluid Fields ---

      DO iCF = 1, nCF

        IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
                 uCF_P(:,iCF), 1, Zero, uCF_L(:,iCF), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
                 uCF_K(:,iCF), 1, Zero, uCF_R(:,iCF), 1 )
        IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

      END DO

      ! --- Interpolate Geometry Fields ---

      ! --- Face States (Average of Left and Right States) ---

      G_F = Zero

      DO iGF = iGF_h_1, iGF_h_3

        IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                 G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                 G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )
        IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

        G_F(:,iGF) = MAX( G_F(:,iGF), SqrtTiny )

      END DO

      CALL ComputeGeometryX_FromScaleFactors( G_F )

      ! --- Lapse Function ---

      IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               G_P(:,iGF_Alpha), 1, Zero, G_F(:,iGF_Alpha), 1 )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               G_K(:,iGF_Alpha), 1, Half, G_F(:,iGF_Alpha), 1 )
      IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

      G_F(:,iGF_Alpha) = MAX( G_F(:,iGF_Alpha), SqrtTiny )

      ! --- Left State Primitive, etc. ---

      IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_CompPrim )
      CALL Euler_ComputePrimitive &
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
      IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_CompPrim )

      CALL ComputePressureFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), P_L  )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L )

      DO iNodeX_X2 = 1, nDOFX_X2

        EigVals_L(iNodeX_X2,:) &
          = Euler_Eigenvalues &
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
          = Euler_Flux_X2 &
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

      IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_CompPrim )
      CALL Euler_ComputePrimitive &
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
      IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_CompPrim )

      CALL ComputePressureFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), P_R )

      CALL ComputeSoundSpeedFromPrimitive &
             ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R )

      DO iNodeX_X2 = 1, nDOFX_X2

        EigVals_R(iNodeX_X2,:) &
          = Euler_Eigenvalues &
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
          = Euler_Flux_X2 &
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
          = Euler_AlphaMiddle &
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
                G_F      (iNodeX_X2,iGF_Alpha),    &
                G_F      (iNodeX_X2,iGF_Beta_2),   &
                AlphaPls, AlphaMns )

        NumericalFlux(iNodeX_X2,:) &
!          = Euler_NumericalFlux_X2_HLLC &
          = Euler_NumericalFlux_HLL &
              ( uCF_L    (iNodeX_X2,:),            &
                uCF_R    (iNodeX_X2,:),            &
                Flux_X2_L(iNodeX_X2,:),            &
                Flux_X2_R(iNodeX_X2,:),            &
                G_F      (iNodeX_X2,iGF_Gm_dd_22), &
                uPF_L    (iNodeX_X2,iPF_V2),       &
                uPF_R    (iNodeX_X2,iPF_V2),       &
                P_L      (iNodeX_X2),              &
                P_R      (iNodeX_X2),              &
                G_F      (iNodeX_X2,iGF_Alpha),    &
                G_F      (iNodeX_X2,iGF_Beta_2),   &
                AlphaPls, AlphaMns, AlphaMdl )

      END DO

      DO iCF = 1, nCF

        NumericalFlux(:,iCF) &
          = dX1 * dX3 * WeightsX_X2 * G_F(:,iGF_Alpha) * G_F(:,iGF_SqrtGm) &
              * NumericalFlux(:,iCF)

      END DO

      ! --- Contribution to This Element ---

      IF( iX2 .LT. iX_E0(2) + 1 )THEN

        DO iCF = 1, nCF

          IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
          CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Dn, nDOFX_X2, &
                      NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )
          IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

        END DO

      END IF

      ! --- Contribution to Previous Element ---

      IF( iX2 .GT. iX_B0(2) )THEN

        DO iCF = 1, nCF

          IF( TimeIt ) CALL TimersStart_Euler( Timer_Euler_MV )
          CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Up, nDOFX_X2, &
                      NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2-1,iX3,iCF), 1 )
          IF( TimeIt ) CALL TimersStop_Euler( Timer_Euler_MV )

        END DO

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    IF( iX_E0(3) == iX_B0(3) ) RETURN

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

#ifdef HYDRO_NONRELATIVISTIC

    CALL ComputeIncrement_Geometry_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#elif HYDRO_RELATIVISTIC

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

      IF( nDimsX .EQ. 2 )THEN
        DO iGF = 1, nGF
          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)
        END DO
      END IF

      CALL Euler_ComputePrimitive &
             ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
               uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
               uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
               uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
               G_K(:,iGF_Gm_dd_11), &
               G_K(:,iGF_Gm_dd_22), &
               G_K(:,iGF_Gm_dd_33) )

      CALL ComputePressureFromPrimitive &
             ( uPF_K(:,iPF_D ), uPF_K(:,iPF_E), uPF_K(:,iPF_Ne), P_K )

      DO iNodeX = 1, nDOFX

        Stress(iNodeX,1:3) &
          = Euler_StressTensor_Diagonal &
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

      ! --- Compute Increments ---

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1) &
            + ( Stress(:,2) * dh2dX1(:) ) / G_K(:,iGF_h_2)  &
            + ( Stress(:,3) * dh3dX1(:) ) / G_K(:,iGF_h_3)

      dU(:,iX1,iX2,iX3,iCF_S2) &
        = dU(:,iX1,iX2,iX3,iCF_S2) &
            + ( Stress(:,3) * dh3dX2(:) ) / G_K(:,iGF_h_3)

    END DO
    END DO
    END DO


  END SUBROUTINE ComputeIncrement_Geometry_NonRelativistic


  SUBROUTINE ComputeIncrement_Geometry_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:), &
      U (:,iX_B1(1):,iX_B1(2):,iX_B1(3):,:)
    REAL(DP), INTENT(inout) :: &
      dU(:,iX_B0(1):,iX_B0(2):,iX_B0(3):,:)

    ! --- This subroutine currently assumes that the shift-vector
    !     is identically zero, which implies that the extrinsic
    !     curvature is identically zero ---

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: dh1dX1(nDOFX), dh2dX1(nDOFX), dh3dX1(nDOFX), &
                dh1dX2(nDOFX), dh2dX2(nDOFX), dh3dX2(nDOFX), &
                dh1dX3(nDOFX), dh2dX3(nDOFX), dh3dX3(nDOFX)
    REAL(DP) :: dadx1(nDOFX), dadx2(nDOFX), dadx3(nDOFX)
    REAL(DP) :: Stress(nDOFX,3)
    REAL(DP) :: uCF_K(nDOFX,nCF), uPF_K(nDOFX,nPF), G_K(nDOFX,nGF)
    REAL(DP) :: G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF), &
                G_P_X2(nDOFX,nGF), G_N_X2(nDOFX,nGF), &
                G_P_X3(nDOFX,nGF), G_N_X3(nDOFX,nGF)
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF), &
                G_X2_Dn(nDOFX_X2,nGF), G_X2_Up(nDOFX_X2,nGF), &
                G_X3_Dn(nDOFX_X3,nGF), G_X3_Up(nDOFX_X3,nGF)

    dadx1 = Zero
    dadx2 = Zero
    dadx3 = Zero

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

      IF( nDimsX .EQ. 2 )THEN
        DO iGF = 1, nGF
          G_P_X2(:,iGF) = G(:,iX1,iX2-1,iX3,iGF)
          G_N_X2(:,iGF) = G(:,iX1,iX2+1,iX3,iGF)
        END DO
      END IF

      IF( nDimsX .EQ. 3 )THEN
        DO iGF = 1, nGF
          G_P_X3(:,iGF) = G(:,iX1,iX2,iX3-1,iGF)
          G_N_X3(:,iGF) = G(:,iX1,iX2,iX3+1,iGF)
        END DO
      END IF

      CALL Euler_ComputePrimitive &
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
          = Euler_StressTensor_Diagonal &
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
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_1), 1, One,  dh1dX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_1), 1, One,  dh1dX1, 1 )

      dh1dx1 = dh1dx1 / ( WeightsX_q * dX1 )

      ! --- dh2dx1 ---

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_2), 1, One,  dh2dX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_2), 1, One,  dh2dX1, 1 )

      dh2dx1 = dh2dx1 / ( WeightsX_q * dX1 )

      ! --- dh3dx1 ---

      CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_h_3), 1, One,  dh3dX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_h_3), 1, One,  dh3dX1, 1 )

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
                  WeightsX_X1 * G_X1_Up(:,iGF_Alpha), 1, Zero, dadX1, 1 )
      CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                  WeightsX_X1 * G_X1_Dn(:,iGF_Alpha), 1, One,  dadX1, 1 )
      CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                  WeightsX_q  * G_K    (:,iGF_Alpha), 1, One,  dadX1, 1 )

      dadx1 = dadx1 / ( WeightsX_q * dX1 )

      dU(:,iX1,iX2,iX3,iCF_S1) &
        = dU(:,iX1,iX2,iX3,iCF_S1)                                &
            + G_K(:,iGF_Alpha)                                    &
                * ( ( Stress(:,1) * dh1dX1 ) / G_K(:,iGF_h_1)     &
                    + ( Stress(:,2) * dh2dX1 ) / G_K(:,iGF_h_2)   &
                    + ( Stress(:,3) * dh3dX1 ) / G_K(:,iGF_h_3) ) &
            - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx1

      IF( nDimsX .GT. 1 )THEN

        ! --- Scale Factor Derivatives wrt X2 ---

        ! --- Face States (Average of Left and Right States) ---

        DO iGF = iGF_h_1, iGF_h_3

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_P_X2(:,iGF), 1, Zero, G_X2_Dn(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_K   (:,iGF), 1, Half, G_X2_Dn(:,iGF), 1 )

          G_X2_Dn(:,iGF) = MAX( G_X2_Dn(:,iGF), SqrtTiny )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G_K   (:,iGF), 1, Zero, G_X2_Up(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G_N_X2(:,iGF), 1, Half, G_X2_Up(:,iGF), 1 )

          G_X2_Up(:,iGF) = MAX( G_X2_Up(:,iGF), SqrtTiny )

        END DO

        ! --- dh1dx2 ---

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_h_1), 1, Zero, dh1dX2, 1 )
        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_h_1), 1, One,  dh1dX2, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_1), 1, One,  dh1dX2, 1 )

        dh1dx2 = dh1dx2 / ( WeightsX_q * dX2 )

        ! --- dh2dx2 ---

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_h_2), 1, Zero, dh2dX2, 1 )
        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_h_2), 1, One,  dh2dX2, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_2), 1, One,  dh2dX2, 1 )

        dh2dx2 = dh2dx2 / ( WeightsX_q * dX2 )

        ! --- dh3dx2 ---

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_h_3), 1, Zero, dh3dX2, 1 )
        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_h_3), 1, One,  dh3dX2, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_3), 1, One,  dh3dX2, 1 )

        dh3dx2 = dh3dx2 / ( WeightsX_q * dX2 )

        ! --- Lapse Function Derivative wrt X2 ---

        ! --- Face States (Average of Left and Right States) ---

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                 G_P_X2(:,iGF_Alpha), 1, Zero, G_X2_Dn(:,iGF_Alpha), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                 G_K   (:,iGF_Alpha), 1, Half, G_X2_Dn(:,iGF_Alpha), 1 )

        G_X2_Dn(:,iGF_Alpha) = MAX( G_X2_Dn(:,iGF_Alpha), SqrtTiny )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                 G_K   (:,iGF_Alpha), 1, Zero, G_X2_Up(:,iGF_Alpha), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                 G_N_X2(:,iGF_Alpha), 1, Half, G_X2_Up(:,iGF_Alpha), 1 )

        G_X2_Up(:,iGF_Alpha) = MAX( G_X2_Up(:,iGF_Alpha), SqrtTiny )

        ! --- dadx2 ---

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Up(:,iGF_Alpha), 1, Zero, dadX2, 1 )
        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2 * G_X2_Dn(:,iGF_Alpha), 1, One,  dadX2, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Alpha), 1, One,  dadX2, 1 )

        dadx2 = dadx2 / ( WeightsX_q * dX2 )

        dU(:,iX1,iX2,iX3,iCF_S2) &
          = dU(:,iX1,iX2,iX3,iCF_S2)                                &
              + G_K(:,iGF_Alpha)                                    &
                  * ( ( Stress(:,1) * dh1dX2 ) / G_K(:,iGF_h_1)     &
                      + ( Stress(:,2) * dh2dX2 ) / G_K(:,iGF_h_2)   &
                      + ( Stress(:,3) * dh3dX2 ) / G_K(:,iGF_h_3) ) &
              - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx2

      END IF

      IF( nDimsX .GT. 2 )THEN

        ! --- Scale Factor Derivatives wrt X3 ---

        ! --- Face States (Average of Left and Right States) ---

        DO iGF = iGF_h_1, iGF_h_3

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                   G_P_X3(:,iGF), 1, Zero, G_X3_Dn(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                   G_K   (:,iGF), 1, Half, G_X3_Dn(:,iGF), 1 )

          G_X3_Dn(:,iGF) = MAX( G_X3_Dn(:,iGF), SqrtTiny )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                   G_K   (:,iGF), 1, Zero, G_X3_Up(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                   G_N_X3(:,iGF), 1, Half, G_X3_Up(:,iGF), 1 )

          G_X3_Up(:,iGF) = MAX( G_X3_Up(:,iGF), SqrtTiny )

        END DO

        ! --- dh1dx3 ---

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_h_1), 1, Zero, dh1dX3, 1 )
        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_h_1), 1, One,  dh1dX3, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_1), 1, One,  dh1dX3, 1 )

        dh1dx3 = dh1dx3 / ( WeightsX_q * dX3 )

        ! --- dh2dx3 ---

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_h_2), 1, Zero, dh2dX3, 1 )
        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_h_2), 1, One,  dh2dX3, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_2), 1, One,  dh2dX3, 1 )

        dh2dx3 = dh2dx3 / ( WeightsX_q * dX3 )

        ! --- dh3dx3 ---

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_h_3), 1, Zero, dh3dX3, 1 )
        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                  WeightsX_X3 * G_X3_Dn(:,iGF_h_3), 1, One,  dh3dX3, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_h_3), 1, One,  dh3dX3, 1 )

        dh3dx3 = dh3dx3 / ( WeightsX_q * dX3 )

        ! --- Lapse Function Derivative wrt X3 ---

        ! --- Face States (Average of Left and Right States) ---

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                 G_P_X3(:,iGF_Alpha), 1, Zero, G_X3_Dn(:,iGF_Alpha), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                 G_K   (:,iGF_Alpha), 1, Half, G_X3_Dn(:,iGF_Alpha), 1 )

        G_X3_Dn(:,iGF_Alpha) = MAX( G_X3_Dn(:,iGF_Alpha), SqrtTiny )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                 G_K   (:,iGF_Alpha), 1, Zero, G_X3_Up(:,iGF_Alpha), 1 )
        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                 G_N_X3(:,iGF_Alpha), 1, Half, G_X3_Up(:,iGF_Alpha), 1 )

        G_X3_Up(:,iGF_Alpha) = MAX( G_X3_Up(:,iGF_Alpha), SqrtTiny )

        ! --- dadx3 ---

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Up(:,iGF_Alpha), 1, Zero, dadX3, 1 )
        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3 * G_X3_Dn(:,iGF_Alpha), 1, One,  dadX3, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q  * G_K    (:,iGF_Alpha), 1, One,  dadX3, 1 )

        dadx3 = dadx3 / ( WeightsX_q * dX3 )

        dU(:,iX1,iX2,iX3,iCF_S3) &
          = dU(:,iX1,iX2,iX3,iCF_S3)                                &
              + G_K(:,iGF_Alpha)                                    &
                  * ( ( Stress(:,1) * dh1dX3 ) / G_K(:,iGF_h_1)     &
                      + ( Stress(:,2) * dh2dX3 ) / G_K(:,iGF_h_2)   &
                      + ( Stress(:,3) * dh3dX3 ) / G_K(:,iGF_h_3) ) &
              - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx3

      END IF

      ! --- Compute Energy Increment (missing extrinsic curvature term) ---

      dU(:,iX1,iX2,iX3,iCF_E) &
        = dU(:,iX1,iX2,iX3,iCF_E) &
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

#ifdef HYDRO_NONRELATIVISTIC

    CALL ComputeIncrement_Gravity_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

#elif HYDRO_RELATIVISTIC

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
