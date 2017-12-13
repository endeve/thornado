MODULE dgDiscretizationModule_Euler_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Pi, TwoPi, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    WeightsX_q, &
    WeightsX_X1, &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    LX_X1_Dn, &
    LX_X1_Up
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE GeometryComputationModule_Beta, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    rhsCF, uAF, nAF, iAF_P, iAF_Gm
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY:  &
    ComputePrimitive_GR, &
    Eigenvalues_GR, &
    AlphaC_GR, &
    Flux_X1_GR, &
    NumericalFlux_X1_HLLC_GR
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromSpecificInternalEnergy, &
    ComputeSoundSpeedFromPrimitive_GR

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_Euler_GR_DG_Explicit

  LOGICAL, PARAMETER :: DisplayTimers = .TRUE.
  REAL(DP) :: Timer_RHS_GR
  REAL(DP) :: Timer_RHS_1_GR, dT_RHS_1_GR
  REAL(DP) :: Timer_RHS_2_GR, dT_RHS_2_GR
  REAL(DP) :: Timer_RHS_3_GR, dT_RHS_3_GR
  REAL(DP) :: Timer_INT_F_GR, dT_INT_F_GR
  REAL(DP) :: Timer_INT_G_GR, dT_INT_G_GR
  REAL(DP) :: Timer_FLX_N_GR, dT_FLX_N_GR

CONTAINS


  SUBROUTINE ComputeIncrement_Euler_GR_DG_Explicit &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)   :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    REAL(DP) :: ErrorL1, ErrorIn, Error, X1

    INTEGER  :: iX1, iX2, iX3, iCF, iNodeX, iNodeX1
    REAL(DP) :: dX1, dX2, dX3

    dU = Zero

    CALL ComputeIncrement_Divergence_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    ! --- Multiply Inverse Mass Matrix ---

    DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)

        dX3 = MeshX(3) % Width(iX3)

        DO iX2 = iX_B0(2), iX_E0(2)

          dX2 = MeshX(2) % Width(iX2)

          DO iX1 = iX_B0(1), iX_E0(1)

            dX1 = MeshX(1) % Width(iX1)

            dU(:,iX1,iX2,iX3,iCF) &
              = dU(:,iX1,iX2,iX3,iCF) &
                  / ( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                        * dX1 * dX2 * dX3 )

          END DO
        END DO
      END DO
    END DO

    ! --- Compute Error ---

    ErrorL1 = 0.0_DP
    ErrorIn = 0.0_DP
    
    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            Error &
              = ABS( - ( 1.0_DP - 0.01_DP )**( -0.5_DP ) * 0.1_DP &
                     * Pi * COS( TwoPi * X1 ) &
                     - dU(iNodeX,iX1,iX2,iX3,iCF_D) )

            ErrorL1 = ErrorL1 + Error
            ErrorIn = MAX( ErrorIn, Error )

          END DO

        END DO
      END DO
    END DO

    ErrorL1 = ErrorL1 / REAL( nDOFX*nX(1)*nX(2)*nX(3) )

    WRITE(*,*)
    WRITE(*,'(A6,A,ES10.4E2)') &
      '', 'ErrorL1: ', ErrorL1
    WRITE(*,'(A6,A,ES10.4E2)') &
      '', 'ErrorIn: ', ErrorIn
    WRITE(*,*)

    CALL ComputeIncrement_Geometry &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

  END SUBROUTINE ComputeIncrement_Euler_GR_DG_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iAF, iNodeX, iNodeX_X1, iNodeX1
    REAL(DP) :: dX2, dX3
    REAL(DP) :: Alpha, Beta_1, Beta_2, Beta_3
    REAL(DP) :: AlphaMns, AlphaMdl, AlphaPls
    REAL(DP), DIMENSION(nDOFX_X1)     :: Cs_L, Cs_R
    REAL(DP), DIMENSION(nCF,nDOFX_X1) :: Lambda_L, Lambda_R
    REAL(DP), DIMENSION(nDOFX_X1,nPF) :: uPF_L, uPF_R
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: uCF_L, uCF_R
    REAL(DP), DIMENSION(nDOFX_X1,nGF) :: G_F
    REAL(DP), DIMENSION(nDOFX)        :: P_K, P_P
    REAL(DP), DIMENSION(nDOFX_X1,nAF) :: uAF_L, uAF_R
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: Flux_X1_L, Flux_X1_R
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: EigVals_L, EigVals_R, EigVals
    REAL(DP), DIMENSION(nDOFX_X1,nCF) :: NumericalFlux
    REAL(DP), DIMENSION(nDOFX,   nPF) :: uPF_P, uPF_K
    REAL(DP), DIMENSION(nDOFX,   nCF) :: uCF_P, uCF_K    
    REAL(DP), DIMENSION(nDOFX,   nGF) :: G_P, G_K
    REAL(DP), DIMENSION(nDOFX,   nAF) :: uAF_P, uAF_K
    REAL(DP), DIMENSION(nDOFX,   nCF) :: Flux_X1_q

    Timer_RHS_1_GR = 0.0_DP
    Timer_RHS_2_GR = 0.0_DP
    Timer_RHS_3_GR = 0.0_DP
    Timer_INT_F_GR = 0.0_DP
    Timer_INT_G_GR = 0.0_DP
    Timer_FLX_N_GR = 0.0_DP

    CALL Timer_Start( Timer_RHS_GR )

    DO iX3 = iX_B0(3), iX_E0(3)

      dX3 = MeshX(3) % Width(iX3)

      DO iX2 = iX_B0(2), iX_E0(2)

        dX2 = MeshX(2) % Width(iX2)

        DO iX1 = iX_B0(1), iX_E0(1) + 1

          DO iCF = 1, nCF

            uCF_P(:,iCF) = U(:,iX1-1,iX2,iX3,iCF)
            uCF_K(:,iCF) = U(:,iX1,  iX2,iX3,iCF)

          END DO

          DO iGF = 1, nGF

            G_P(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
            G_K(:,iGF) = G(:,iX1,  iX2,iX3,iGF)

          END DO

          P_P(:) = uAF(:,iX1-1,iX2,iX3,iAF_P)
          P_K(:) = uAF(:,iX1,  iX2,iX3,iAF_P)

          !--------------------
          ! --- Volume Term ---
          !--------------------

          IF( iX1 < iX_E0(1) + 1 )THEN

            CALL ComputePrimitive_GR &
                   ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                     uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                     uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                     uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                     P_K(:), &
                     G_K(:,iGF_Gm_dd_11), &
                     G_K(:,iGF_Gm_dd_22), &
                     G_K(:,iGF_Gm_dd_33) )

            DO iNodeX = 1, nDOFX

              Flux_X1_q(iNodeX,1:nCF)           &
                = Flux_X1_GR                    &
                    ( uPF_K(iNodeX,iPF_D ),     &
                      uPF_K(iNodeX,iPF_V1),     &
                      uPF_K(iNodeX,iPF_V2),     &
                      uPF_K(iNodeX,iPF_V3),     &
                      uPF_K(iNodeX,iPF_E ),     &                    
                      uPF_K(iNodeX,iPF_Ne),     &
                      uAF_K(iNodeX,iAF_P ),     &                    
                      G_K(iNodeX,iGF_Gm_dd_11), &
                      G_K(iNodeX,iGF_Gm_dd_22), &
                      G_K(iNodeX,iGF_Gm_dd_33), &
                      G_K(iNodeX,iGF_Alpha),    &
                      G_K(iNodeX,iGF_Beta_1) )

            END DO

            CALL Timer_Start( dT_RHS_1_GR )

            DO iCF = 1, nCF

              Flux_X1_q(:,iCF) &
                = dX2 * dX3 * WeightsX_q(:) * G_K(:,iGF_Alpha) &
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
                   P_P(:), 1, Zero, uAF_L(:,iAF_P), 1 )

          ! --- Right State Pressure ---

          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                   P_K(:), 1, Zero, uAF_R(:,iAF_P), 1 )

          CALL Timer_Stop( dT_INT_F_GR )

          CALL Timer_Add( Timer_INT_F_GR, dT_INT_F_GR )

          ! --- Interpolate Geometry Fields ---

          CALL Timer_Start( dT_INT_G_GR )

          G_F = Zero

          ! --- Face States (Average of Left and Right States) ---

          ! --- Scale Factors ---

          DO iGF = iGF_h_1, iGF_h_3

            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                     G_P(:,iGF), 1, Zero, G_F(:,iGF), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                     G_K(:,iGF), 1, Half, G_F(:,iGF), 1 )

            G_F(1:nDOFX_X1,iGF) &
              = MAX( G_F(1:nDOFX_X1,iGF), SqrtTiny )

          END DO

          CALL ComputeGeometryX_FromScaleFactors( G_F(:,:) )

          ! --- Lapse Function ---

          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                   G_P(:,iGF_Alpha), 1, Zero, G_F(:,iGF_Alpha), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                   G_K(:,iGF_Alpha), 1, Half, G_F(:,iGF_Alpha), 1 )

          G_F(1:nDOFX_X1,iGF_Alpha) &
            = MAX( G_F(1:nDOFX_X1,iGF_Alpha), SqrtTiny )

          CALL Timer_Stop( dT_INT_G_GR )

          CALL Timer_Add( Timer_INT_G_GR, dT_INT_G_GR )

          ! --- Left State Primitive, etc. ---

          CALL ComputePrimitive_GR &
                 ( uCF_L(:,iCF_D ), uCF_L(:,iCF_S1), uCF_L(:,iCF_S2), &
                   uCF_L(:,iCF_S3), uCF_L(:,iCF_E ), uCF_L(:,iCF_Ne), &
                   uPF_L(:,iPF_D ), uPF_L(:,iPF_V1), uPF_L(:,iPF_V2), &
                   uPF_L(:,iPF_V3), uPF_L(:,iPF_E ), uPF_L(:,iPF_Ne), &
                   uAF_L(:,iAF_P), &
                   G_F(:,iGF_Gm_dd_11), &
                   G_F(:,iGF_Gm_dd_22), &
                   G_F(:,iGF_Gm_dd_33) )

          CALL ComputeSoundSpeedFromPrimitive_GR &
                 ( uPF_L(:,iPF_D), uPF_L(:,iPF_E), uPF_L(:,iPF_Ne), Cs_L(:) )

          DO iNodeX_X1 = 1, nDOFX_X1

            Lambda_L(:,iNodeX_X1) &
              = Eigenvalues_GR &
                  ( uPF_L(iNodeX_X1,iPF_V1),       &
                    uPF_L(iNodeX_X1,iPF_V2),       &
                    uPF_L(iNodeX_X1,iPF_V3),       &
                    uPF_L(iNodeX_X1,iPF_V1),       &
                    Cs_L (iNodeX_X1),              &
                    G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                    G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                    G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                    G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                    G_F  (iNodeX_X1,iGF_Alpha),    &
                    G_F  (iNodeX_X1,iGF_Beta_1) )


            Flux_X1_L(iNodeX_X1,1:nCF)          &
              = Flux_X1_GR                         &
                  ( uPF_L(iNodeX_X1,iPF_D ),       &
                    uPF_L(iNodeX_X1,iPF_V1),       &
                    uPF_L(iNodeX_X1,iPF_V2),       &
                    uPF_L(iNodeX_X1,iPF_V3),       &
                    uPF_L(iNodeX_X1,iPF_E ),       &
                    uPF_L(iNodeX_X1,iPF_Ne),       &
                    uAF_L(iNodeX_X1,iAF_P ),       &
                    G_F(iNodeX_X1,iGF_Gm_dd_11), &
                    G_F(iNodeX_X1,iGF_Gm_dd_22), &
                    G_F(iNodeX_X1,iGF_Gm_dd_33), &
                    G_F(iNodeX_X1,iGF_Alpha),    &
                    G_F(iNodeX_X1,iGF_Beta_1) )

          END DO

          ! --- Right State Primitive, etc. ---

          CALL ComputePrimitive_GR &
                 ( uCF_R(:,iCF_D ), uCF_R(:,iCF_S1), uCF_R(:,iCF_S2), &
                   uCF_R(:,iCF_S3), uCF_R(:,iCF_E ), uCF_R(:,iCF_Ne), &
                   uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), uPF_R(:,iPF_V2), &
                   uPF_R(:,iPF_V3), uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
                   uAF_R(:,iAF_P), &
                   G_F(:,iGF_Gm_dd_11), &
                   G_F(:,iGF_Gm_dd_22), &
                   G_F(:,iGF_Gm_dd_33) )

          CALL ComputeSoundSpeedFromPrimitive_GR &
                 ( uPF_R(:,iPF_D), uPF_R(:,iPF_E), uPF_R(:,iPF_Ne), Cs_R(:) )

          DO iNodeX_X1 = 1, nDOFX_X1


            Lambda_R(:,iNodeX_X1) &
              = Eigenvalues_GR &
                  ( uPF_R(iNodeX_X1,iPF_V1),       &
                    uPF_R(iNodeX_X1,iPF_V2),       &
                    uPF_R(iNodeX_X1,iPF_V3),       &
                    uPF_R(iNodeX_X1,iPF_V1),       &
                    Cs_R (iNodeX_X1),              &
                    G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                    G_F  (iNodeX_X1,iGF_Gm_dd_22), &
                    G_F  (iNodeX_X1,iGF_Gm_dd_33), &
                    G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                    G_F  (iNodeX_X1,iGF_Alpha),    &
                    G_F  (iNodeX_X1,iGF_Beta_1) )


            Flux_X1_R(iNodeX_X1,1:nCF) &
              = Flux_X1_GR                         &
                  ( uPF_R(iNodeX_X1,iPF_D ),       &
                    uPF_R(iNodeX_X1,iPF_V1),       &
                    uPF_R(iNodeX_X1,iPF_V2),       &
                    uPF_R(iNodeX_X1,iPF_V3),       &
                    uPF_R(iNodeX_X1,iPF_E ),       &
                    uPF_R(iNodeX_X1,iPF_Ne),       &
                    uAF_L(iNodeX_X1,iAF_P ),       &
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
              = AlphaC_GR &
                  ( uCF_L(iNodeX_X1,1:nCF), Flux_X1_L(iNodeX_X1,1:nCF), &
                    uCF_R(iNodeX_X1,1:nCF), Flux_X1_R(iNodeX_X1,1:nCF), &
                    G_F  (iNodeX_X1,iGF_Gm_dd_11), &
                    G_F  (iNodeX_X1,iGF_Beta_1),  &
                    AlphaPls, AlphaMns )

            NumericalFlux( iNodeX_X1,:)                         &
              = NumericalFlux_X1_HLLC_GR                           &
                  ( uCF_L( iNodeX_X1, 1:nCF ),                     &
                    uCF_R( iNodeX_X1,1:nCF ),                      &
                    Flux_X1_L( iNodeX_X1,1:nCF ),               &
                    Flux_X1_R( iNodeX_X1,1:nCF ),               &
                    AlphaPls,      &
                    AlphaMns, AlphaMdl, nCF, &
                    uPF_L( iNodeX_X1, iPF_V1 ),                    &
                    uPF_R( iNodeX_X1, iPF_V1 ),                    &
                    uAF_L( iNodeX_X1, iAF_P ),                     &
                    uAF_R( iNodeX_X1, iAF_P ),                     &
                    G_F( iNodeX_X1, iGF_Beta_1 ),                &
                    G_F( iNodeX_X1, iGF_Gm_dd_11 ) )


          END DO

          DO iCF = 1, nCF

            NumericalFlux(:,iCF) &
              = dX2 * dX3 * WeightsX_X1(:) * G_F(:,iGF_Alpha) &
                  * G_F(:,iGF_SqrtGm) * NumericalFlux(:,iCF)

          END DO

          CALL Timer_Stop( dT_FLX_N_GR )

          CALL Timer_Add( Timer_FLX_N_GR, dT_FLX_N_GR )

          ! --- Contribution to This Element ---

          CALL Timer_Start( dT_RHS_2_GR )

          IF( iX1 < iX_E0(1) + 1 )THEN

            DO iCF = 1, nCF

              CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Dn, nDOFX_X1, &
                          NumericalFlux(:,iCF), 1, One, &
                          dU(:,iX1,iX2,iX3,iCF), 1 )

            END DO

          END IF

          CALL Timer_Stop( dT_RHS_2_GR )

          CALL Timer_Add( Timer_RHS_2_GR, dT_RHS_2_GR )

          ! --- Contribution to Previous Element ---

          CALL Timer_Start( dT_RHS_3_GR )

          IF( iX1 > iX_B0(1) )THEN

            DO iCF = 1, nCF

              CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Up, nDOFX_X1, &
                          NumericalFlux(:,iCF), 1, &
                          One, dU(:,iX1-1,iX2,iX3,iCF), 1 )

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


  SUBROUTINE ComputeIncrement_Geometry &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

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


END MODULE dgDiscretizationModule_Euler_GR
