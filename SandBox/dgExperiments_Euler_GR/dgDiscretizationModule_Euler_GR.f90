MODULE dgDiscretizationModule_Euler_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Pi, TwoPi, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1,    &
    WeightsX_q,  &
    WeightsX_X1, &
    NodeNumberTableX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    LX_X1_Dn, &
    LX_X1_Up
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
    uAF, nAF, iAF_P, iAF_Gm
  USE BoundaryConditionsModule_Beta, ONLY: &
    ApplyBoundaryConditions_Fluid
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY:  &
    ComputePrimitive_GR,      &
    Eigenvalues_GR,           &
    AlphaC_GR,                &
    Flux_X1_GR,               &
    StressTensor_Diagonal,    &
    NumericalFlux_X1_HLLC_GR, &
    NumericalFlux_X1_LLF_GR
  USE SlopeLimiterModule_Euler_GR, ONLY: &
    ApplySlopeLimiter_Euler_GR
  USE PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromSpecificInternalEnergy, &
    ComputeSoundSpeedFromPrimitive_GR

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeIncrement_Euler_GR_DG_Explicit

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


  SUBROUTINE ComputeIncrement_Euler_GR_DG_Explicit &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out)   :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    REAL(DP) :: ErrorL1, ErrorIn, Error, X1

    INTEGER  :: iX1, iX2, iX3, iCF, iNodeX, iNodeX1
    REAL(DP) :: dX1, dX2, dX3

    dU = Zero

!    WRITE(*,*) 'CALL ApplyBoundaryConditions_Fluid 1'
    CALL ApplyBoundaryConditions_Fluid &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U )

!    WRITE(*,*) 'CALL ApplySlopeLimiter_Euler_GR'
    CALL ApplySlopeLimiter_Euler_GR &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

!    WRITE(*,*) 'CALL ApplyPositivityLimiter'
    CALL ApplyPositivityLimiter &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

!    WRITE(*,*) 'CALL ApplyBoundaryConditions_Fluid 2'
    CALL ApplyBoundaryConditions_Fluid &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U )

!    WRITE(*,*) 'CALL ComputeIncrement_Divergence_X1'
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

            dU(:,iX1,iX2,iX3,iCF)                                 &
              = dU(:,iX1,iX2,iX3,iCF)                             &
                  / ( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) &
                        * dX1 * dX2 * dX3 )
            
          END DO
        END DO
      END DO
    END DO

!    WRITE(*,*) 'CALL ComputeIncrement_Geometry'
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

!          WRITE(*,*) 'iX1,iX2,iX3' , iX1,iX2,iX3
           
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
                 P_K(:),                                            &
                 G_K(:,iGF_Gm_dd_11),                               &
                 G_K(:,iGF_Gm_dd_22),                               &
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
                      P_K(iNodeX),     &                    
                      G_K(iNodeX,iGF_Gm_dd_11), &
                      G_K(iNodeX,iGF_Gm_dd_22), &
                      G_K(iNodeX,iGF_Gm_dd_33), &
                      G_K(iNodeX,iGF_Alpha),    &
                      G_K(iNodeX,iGF_Beta_1) )

            END DO

            CALL Timer_Start( dT_RHS_1_GR )

            DO iCF = 1, nCF

              Flux_X1_q(:,iCF)                                 &
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
                   uAF_L(:,iAF_P),                                    &
                   G_F(:,iGF_Gm_dd_11),                               &
                   G_F(:,iGF_Gm_dd_22),                               &
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

            Flux_X1_L(iNodeX_X1,1:nCF)             &
              = Flux_X1_GR                         &
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

!          WRITE(*,*) 'Right State Primitive'
!          WRITE(*,*) '  CALL ComputePrimitive'
          CALL ComputePrimitive_GR &
               ( uCF_R(:,iCF_D ), uCF_R(:,iCF_S1), uCF_R(:,iCF_S2), &
                 uCF_R(:,iCF_S3), uCF_R(:,iCF_E ), uCF_R(:,iCF_Ne), &
                 uPF_R(:,iPF_D ), uPF_R(:,iPF_V1), uPF_R(:,iPF_V2), &
                 uPF_R(:,iPF_V3), uPF_R(:,iPF_E ), uPF_R(:,iPF_Ne), &
                 uAF_R(:,iAF_P),                                    &
                 G_F(:,iGF_Gm_dd_11),                               &
                 G_F(:,iGF_Gm_dd_22),                               &
                 G_F(:,iGF_Gm_dd_33) )

!          WRITE(*,*) 'CALL ComputeSoundSpeedFromPrimitive_GR'
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
              = AlphaC_GR &
                  ( uCF_L(iNodeX_X1,1:nCF), Flux_X1_L(iNodeX_X1,1:nCF), &
                    uCF_R(iNodeX_X1,1:nCF), Flux_X1_R(iNodeX_X1,1:nCF), &
                    AlphaPls, AlphaMns,                                 &
                    G_F  (iNodeX_X1,iGF_Gm_dd_11),                      &
                    G_F  (iNodeX_X1,iGF_Beta_1) )

            NumericalFlux( iNodeX_X1,:)                &
              = NumericalFlux_X1_HLLC_GR               &
                  ( uCF_L(iNodeX_X1,1:nCF),            &
                    uCF_R(iNodeX_X1,1:nCF),            &
                    Flux_X1_L(iNodeX_X1,1:nCF),        &
                    Flux_X1_R(iNodeX_X1,1:nCF),        &
                    AlphaPls, AlphaMns, AlphaMdl, nCF, &
                    uPF_L(iNodeX_X1,iPF_V1),           &
                    uPF_R(iNodeX_X1,iPF_V1),           &
                    uAF_L(iNodeX_X1,iAF_P),            &
                    uAF_R(iNodeX_X1,iAF_P),            &
                    G_F(iNodeX_X1,iGF_Beta_1),         &
                    G_F(iNodeX_X1,iGF_Gm_dd_11) )


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

              CALL DGEMV &
                     ( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Dn, nDOFX_X1, &
                       NumericalFlux(:,iCF), 1, One, dU(:,iX1,iX2,iX3,iCF), 1 )

            END DO

          END IF

          CALL Timer_Stop( dT_RHS_2_GR )

          CALL Timer_Add( Timer_RHS_2_GR, dT_RHS_2_GR )

          ! --- Contribution to Previous Element ---

          CALL Timer_Start( dT_RHS_3_GR )

          IF( iX1 > iX_B0(1) )THEN

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


  SUBROUTINE ComputeIncrement_Geometry &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
      U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iNodeX
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: P_K(nDOFX)
    REAL(DP) :: dh1dX1(nDOFX), dh2dX1(nDOFX), dh3dX1(nDOFX)
    REAL(DP) :: dadx1(nDOFX)
    REAL(DP) :: Stress(nDOFX,3)
    REAL(DP) :: uCF_K(nDOFX,nCF)
    REAL(DP) :: uPF_K(nDOFX,nPF)
    REAL(DP) :: G_K(nDOFX,nGF)
    REAL(DP) :: G_P_X1(nDOFX,nGF), G_N_X1(nDOFX,nGF)
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF)

    CALL Timer_Start( Timer_Geo )

    DO iX3 = iX_B0(3), iX_E0(3)

      dX3 = MeshX(3) % Width(iX3)

      DO iX2 = iX_B0(2), iX_E0(2)

        dX2 = MeshX(2) % Width(iX2)

        DO iX1 = iX_B0(1), iX_E0(1)

          dX1 = MeshX(1) % Width(iX1)

          DO iCF = 1, nCF

            uCF_K(:,iCF) = U(:,iX1,iX2,iX3,iCF)

          END DO

          DO iGF = 1, nGF

            G_P_X1(:,iGF) = G(:,iX1-1,iX2,iX3,iGF)
            G_K   (:,iGF) = G(:,iX1,  iX2,iX3,iGF)
            G_N_X1(:,iGF) = G(:,iX1+1,iX2,iX3,iGF)

          END DO

          P_K(:) = uAF(:,iX1,iX2,iX3,iAF_P)

          CALL ComputePrimitive_GR &
               ( uCF_K(:,iCF_D ), uCF_K(:,iCF_S1), uCF_K(:,iCF_S2), &
                 uCF_K(:,iCF_S3), uCF_K(:,iCF_E ), uCF_K(:,iCF_Ne), &
                 uPF_K(:,iPF_D ), uPF_K(:,iPF_V1), uPF_K(:,iPF_V2), &
                 uPF_K(:,iPF_V3), uPF_K(:,iPF_E ), uPF_K(:,iPF_Ne), &
                 P_K(:),                                            &
                 G_K(:,iGF_Gm_dd_11),                               &
                 G_K(:,iGF_Gm_dd_22),                               &
                 G_K(:,iGF_Gm_dd_33) )

          DO iNodeX = 1, nDOFX

            Stress(iNodeX,1:3)            &
              = StressTensor_Diagonal     &
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

          ! --- dh1dx1 ---

          CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Up(:,iGF_h_1), 1, Zero, dh1dX1, 1 )
          CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Dn(:,iGF_h_1), 1,  One, dh1dX1, 1 )
          CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                      WeightsX_q (:) * G_K    (:,iGF_h_1), 1,  One, dh1dX1, 1 )

          dh1dx1 = dh1dx1 / ( WeightsX_q(:) * dX1 )

          ! --- dh2dx1 ---

          CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )
          CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Dn(:,iGF_h_2), 1,  One, dh2dX1, 1 )
          CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                      WeightsX_q (:) * G_K    (:,iGF_h_2), 1,  One, dh2dX1, 1 )

          dh2dx1 = dh2dx1 / ( WeightsX_q(:) * dX1 )

          ! --- dh3dx1 ---

          CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )
          CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Dn(:,iGF_h_3), 1,  One, dh3dX1, 1 )
          CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                      WeightsX_q (:) * G_K    (:,iGF_h_3), 1,  One, dh3dX1, 1 )

          dh3dx1 = dh3dx1 / ( WeightsX_q(:) * dX1 )

          ! --- Lapse Function Derivative wrt X1 ---

          ! --- Face States (Average of Left and Right States) ---

          iGF = iGF_Alpha

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

          ! --- dadx1 ---

          CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Up(:,iGF_Alpha), 1, Zero,  &
                      dadX1, 1 )
          CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                      WeightsX_X1(:) * G_X1_Dn(:,iGF_Alpha), 1,  One,  &
                      dadX1, 1 )
          CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                      WeightsX_q (:) * G_K    (:,iGF_Alpha), 1,  One,  &
                      dadX1, 1 )

          dadx1 = dadx1 / ( WeightsX_q(:) * dX1 )

          ! --- Compute Increments ---

          dU(:,iX1,iX2,iX3,iCF_S1)                                       &
            = dU(:,iX1,iX2,iX3,iCF_S1)                                   &
                + G_K(:,iGF_Alpha)                                       &
                    * ( ( Stress(:,1) * dh1dX1(:) ) / G_K(:,iGF_h_1)     &
                        + ( Stress(:,2) * dh2dX1(:) ) / G_K(:,iGF_h_2)   &
                        + ( Stress(:,3) * dh3dX1(:) ) / G_K(:,iGF_h_3) ) &
                - ( uCF_K(:,iCF_D) + uCF_K(:,iCF_E) ) * dadx1(:)

          dU(:,iX1,iX2,iX3,iCF_E)     &
            = dU(:,iX1,iX2,iX3,iCF_E) &
                - ( uCF_K(:,iCF_S1) / G_K(:,iGF_Gm_dd_11) ) * dadx1(:)

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


END MODULE dgDiscretizationModule_Euler_GR
