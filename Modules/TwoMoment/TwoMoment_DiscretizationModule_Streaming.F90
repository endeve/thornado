#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_EXPLICIT
#endif

MODULE TwoMoment_DiscretizationModule_Streaming

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    TwoPi, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, &
    nNodesE, nDOFE, &
    nDOF
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_Explicit, &
    Timer_Ex_In, &
    Timer_Ex_Div, &
    Timer_Ex_Geometry, &
    Timer_Ex_Permute, &
    Timer_Ex_Interpolate, &
    Timer_Ex_Flux, &
    Timer_Ex_Increment, &
    Timer_Ex_Out
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_q, &
    WeightsX_X1, &
    WeightsX_X2
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, &
    dLXdX2_q, &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, &
    nDOF_X2, &
    nDOF_X3, &
    Weights_q, &
    Weights_X1, &
    Weights_X2, &
    Weights_X3, &
    NodeNumberTable, &
    OuterProduct1D3D, &
    NodeNumbersX
  USE ReferenceElementModule_Lagrange, ONLY: &
    dLdX1_q, &
    dLdX2_q, &
    dLdX3_q, &
    L_X1_Dn, &
    L_X1_Up, &
    L_X2_Dn, &
    L_X2_Up, &
    L_X3_Dn, &
    L_X3_Up
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    CoordinateSystem
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModuleE, ONLY: &
    nGE, &
    iGE_Ep0, &
    iGE_Ep2
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor
  USE TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment, &
    Flux_X1, &
    Flux_X2, &
    Flux_X3, &
    NumericalFlux_LLF, &
    StressTensor_Diagonal

  IMPLICIT NONE
  PRIVATE

  REAL(DP), PARAMETER :: Ones(16) = One

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit

CONTAINS

  SUBROUTINE ComputeIncrement_TwoMoment_Explicit &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out)   :: &
      dU(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iNodeX, iNode, iZ1, iZ2, iZ3, iZ4, iCR, iS
    REAL(DP) :: Tau

    CALL TimersStart( Timer_Explicit )

    ASSOCIATE ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
                dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    CALL ApplyBoundaryConditions_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U )

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B1(4), iZ_E1(4)
          DO iZ3 = iZ_B1(3), iZ_E1(3)
            DO iZ2 = iZ_B1(2), iZ_E1(2)
              DO iZ1 = iZ_B1(1), iZ_E1(1)
                DO iNode = 1, nDOF
                  dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = Zero
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStart( Timer_Ex_Div )

    CALL ComputeIncrement_Divergence_X1 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL ComputeIncrement_Divergence_X2 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL ComputeIncrement_Divergence_X3 &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL TimersStop( Timer_Ex_Div )

    ! --- Multiply Inverse Mass Matrix ---

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF

                  iNodeX = MOD( (iNode-1) / nDOFE, nDOFX ) + 1
                  Tau = One * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

                  dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                    = dU(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) &
                        / ( Weights_q(iNode) * Tau * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) )

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStart( Timer_Ex_Geometry )

    CALL ComputeIncrement_Geometry &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    CALL TimersStop( Timer_Ex_Geometry )

    END ASSOCIATE

#ifdef THORNADO_DEBUG_EXPLICIT
    WRITE(*,'(a20,7i4)')     'MAXLOC(dU)', MAXLOC(dU)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(dU)', MAXVAL(dU)
#endif

    CALL TimersStop( Timer_Explicit )

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit


  SUBROUTINE ComputeIncrement_Divergence_X1 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)    :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNode
    INTEGER  :: iGF, iCR
    REAL(DP) :: dZ(4)
    REAL(DP) :: FF, EF
    REAL(DP) :: GX_P(nDOFX,   nGF)
    REAL(DP) :: GX_K(nDOFX,   nGF)
    REAL(DP) :: GX_F(nDOFX_X1,nGF)
    REAL(DP) :: G_K(nDOF,nGF)
    REAL(DP) :: G_F(nDOF_X1,nGF)
    REAL(DP), DIMENSION(nDOF_X1)     :: absLambda_L
    REAL(DP), DIMENSION(nDOF_X1)     :: absLambda_R
    REAL(DP), DIMENSION(nDOF_X1)     :: alpha
    REAL(DP), DIMENSION(nDOF)        :: Tau
    REAL(DP), DIMENSION(nDOF_X1)     :: Tau_X1
    REAL(DP), DIMENSION(nDOF_X1,nCR) :: uCR_L, uCR_R
    REAL(DP), DIMENSION(nDOF_X1,nPR) :: uPR_L, uPR_R
    REAL(DP), DIMENSION(nDOF_X1,nCR) :: Flux_X1_L
    REAL(DP), DIMENSION(nDOF_X1,nCR) :: Flux_X1_R
    REAL(DP), DIMENSION(nDOF_X1,nCR) :: NumericalFlux
    REAL(DP), DIMENSION(nDOF   ,nCR) :: uCR_P, uCR_K
    REAL(DP), DIMENSION(nDOF   ,nPR) :: uPR_K
    REAL(DP), DIMENSION(nDOF   ,nCR) :: Flux_X1_q

    IF( iZ_E0(2) .EQ. iZ_B0(2) ) RETURN

    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)

        dZ(4) = MeshX(3) % Width(iZ4)

        DO iZ3 = iZ_B0(3), iZ_E0(3)

          dZ(3) = MeshX(2) % Width(iZ3)

          DO iZ2 = iZ_B0(2), iZ_E0(2) + 1

            ! --- Geometry Fields in Element Nodes ---

            DO iGF = 1, nGF

              GX_P(:,iGF) = GX(:,iZ2-1,iZ3,iZ4,iGF) ! --- Previous Element
              GX_K(:,iGF) = GX(:,iZ2,  iZ3,iZ4,iGF) ! --- This     Element

              G_K(1:nDOF,iGF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_K(1:nDOFX,iGF), nDOFX )

            END DO

            ! --- Interpolate Geometry Fields on Shared Face ---

            ! --- Face States (Average of Left and Right States) ---

            ! --- Scale Factors ---

            DO iGF = iGF_h_1, iGF_h_3

              CALL DGEMV &
                     ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                       GX_P(:,iGF), 1, Zero, GX_F(:,iGF), 1 )
              CALL DGEMV &
                     ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                       GX_K(:,iGF), 1, Half, GX_F(:,iGF), 1 )

              GX_F(1:nDOFX_X1,iGF) &
                = MAX( GX_F(1:nDOFX_X1,iGF), SqrtTiny )

              G_F(1:nDOF_X1,iGF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_F(1:nDOFX_X1,iGF), nDOFX_X1 )

            END DO

            CALL ComputeGeometryX_FromScaleFactors( GX_F(:,:) )

            CALL ComputeGeometryX_FromScaleFactors( G_F(:,:) )

            ! --- Lapse Function ---

            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                     GX_P(:,iGF_Alpha), 1, Zero, GX_F(:,iGF_Alpha), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                     GX_K(:,iGF_Alpha), 1, Half, GX_F(:,iGF_Alpha), 1 )

            GX_F(1:nDOFX_X1,iGF_Alpha) &
              = MAX( GX_F(1:nDOFX_X1,iGF_Alpha), SqrtTiny )

            G_F(1:nDOF_X1,iGF_Alpha) &
              = OuterProduct1D3D &
                  ( Ones(1:nDOFE), nDOFE, &
                    GX_F(1:nDOFX_X1,iGF_Alpha), nDOFX_X1 )

            DO iZ1 = iZ_B0(1), iZ_E0(1)

              dZ(1) = MeshE % Width(iZ1)

              ! --- Volume Jacobian in Energy-Position Element ---

              Tau(1:nDOF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_K(:,iGF_SqrtGm), nDOFX )

              Tau_X1(1:nDOF_X1) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_F(:,iGF_SqrtGm), nDOFX_X1 )

              DO iCR = 1, nCR

                uCR_P(:,iCR) = U(:,iZ1,iZ2-1,iZ3,iZ4,iCR,iS)
                uCR_K(:,iCR) = U(:,iZ1,iZ2,  iZ3,iZ4,iCR,iS)

              END DO

              !--------------------
              ! --- Volume Term ---
              !--------------------

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                CALL ComputePrimitive_TwoMoment &
                       ( uCR_K(:,iCR_N ), uCR_K(:,iCR_G1), &
                         uCR_K(:,iCR_G2), uCR_K(:,iCR_G3), &
                         uPR_K(:,iPR_D ), uPR_K(:,iPR_I1), &
                         uPR_K(:,iPR_I2), uPR_K(:,iPR_I3), &
                         G_K(:,iGF_Gm_dd_11), &
                         G_K(:,iGF_Gm_dd_22), &
                         G_K(:,iGF_Gm_dd_33) )

                DO iNode = 1, nDOF

                  FF = FluxFactor &
                         ( uPR_K(iNode,iPR_D ), uPR_K(iNode,iPR_I1), &
                           uPR_K(iNode,iPR_I2), uPR_K(iNode,iPR_I3), &
                           G_K(iNode,iGF_Gm_dd_11), &
                           G_K(iNode,iGF_Gm_dd_22), &
                           G_K(iNode,iGF_Gm_dd_33) )

                  EF = EddingtonFactor( uPR_K(iNode,iPR_D), FF )

                  Flux_X1_q(iNode,1:nCR) &
                    = Flux_X1 &
                        ( uPR_K(iNode,iPR_D ), uPR_K(iNode,iPR_I1), &
                          uPR_K(iNode,iPR_I2), uPR_K(iNode,iPR_I3), &
                          FF, EF, &
                          G_K(iNode,iGF_Gm_dd_11), &
                          G_K(iNode,iGF_Gm_dd_22), &
                          G_K(iNode,iGF_Gm_dd_33) )

                END DO

                DO iCR = 1, nCR

                  Flux_X1_q(:,iCR) &
                    = dZ(3) * dZ(4) * Weights_q(:) &
                        * G_K(:,iGF_Alpha) * Tau(:) * Flux_X1_q(:,iCR)

                  CALL DGEMV &
                         ( 'T', nDOF, nDOF, One, dLdX1_q, nDOF, &
                           Flux_X1_q(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2,iZ3,iZ4,iCR,iS), 1 )

                END DO

              END IF

              !---------------------
              ! --- Surface Term ---
              !---------------------

              ! --- Interpolate Radiation Fields ---

              DO iCR = 1, nCR

                ! --- Interpolate Left State ---

                CALL DGEMV &
                       ( 'N', nDOF_X1, nDOF, One, L_X1_Up, nDOF_X1, &
                         uCR_P(:,iCR), 1, Zero, uCR_L(:,iCR), 1 )

                ! --- Interpolate Right State ---

                CALL DGEMV &
                       ( 'N', nDOF_X1, nDOF, One, L_X1_Dn, nDOF_X1, &
                         uCR_K(:,iCR), 1, Zero, uCR_R(:,iCR), 1 )

              END DO

              ! --- Left State Primitive, etc. ---

              CALL ComputePrimitive_TwoMoment &
                     ( uCR_L(:,iCR_N ), uCR_L(:,iCR_G1), &
                       uCR_L(:,iCR_G2), uCR_L(:,iCR_G3), &
                       uPR_L(:,iPR_D ), uPR_L(:,iPR_I1), &
                       uPR_L(:,iPR_I2), uPR_L(:,iPR_I3), &
                       G_F(:,iGF_Gm_dd_11), &
                       G_F(:,iGF_Gm_dd_22), &
                       G_F(:,iGF_Gm_dd_33) )

              DO iNode = 1, nDOF_X1

                FF = FluxFactor &
                       ( uPR_L(iNode,iPR_D ), uPR_L(iNode,iPR_I1), &
                         uPR_L(iNode,iPR_I2), uPR_L(iNode,iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11), &
                         G_F(iNode,iGF_Gm_dd_22), &
                         G_F(iNode,iGF_Gm_dd_33) )

                EF = EddingtonFactor( uPR_L(iNode,iPR_D), FF )

                Flux_X1_L(iNode,1:nCR) &
                  = Flux_X1 &
                      ( uPR_L(iNode,iPR_D ), uPR_L(iNode,iPR_I1), &
                        uPR_L(iNode,iPR_I2), uPR_L(iNode,iPR_I3), &
                        FF, EF, &
                        G_F(iNode,iGF_Gm_dd_11), &
                        G_F(iNode,iGF_Gm_dd_22), &
                        G_F(iNode,iGF_Gm_dd_33) )

                absLambda_L(iNode) = 1.0_DP

              END DO

              ! --- Right State Primitive, etc. ---

              CALL ComputePrimitive_TwoMoment &
                     ( uCR_R(:,iCR_N ), uCR_R(:,iCR_G1), &
                       uCR_R(:,iCR_G2), uCR_R(:,iCR_G3), &
                       uPR_R(:,iPR_D ), uPR_R(:,iPR_I1), &
                       uPR_R(:,iPR_I2), uPR_R(:,iPR_I3), &
                       G_F(:,iGF_Gm_dd_11), &
                       G_F(:,iGF_Gm_dd_22), &
                       G_F(:,iGF_Gm_dd_33) )

              DO iNode = 1, nDOF_X1

                FF = FluxFactor &
                       ( uPR_R(iNode,iPR_D ), uPR_R(iNode,iPR_I1), &
                         uPR_R(iNode,iPR_I2), uPR_R(iNode,iPR_I3), &
                         G_F(iNode,iGF_Gm_dd_11), &
                         G_F(iNode,iGF_Gm_dd_22), &
                         G_F(iNode,iGF_Gm_dd_33) )

                EF = EddingtonFactor( uPR_R(iNode,iPR_D), FF )

                Flux_X1_R(iNode,1:nCR) &
                  = Flux_X1 &
                      ( uPR_R(iNode,iPR_D ), uPR_R(iNode,iPR_I1), &
                        uPR_R(iNode,iPR_I2), uPR_R(iNode,iPR_I3), &
                        FF, EF, &
                        G_F(iNode,iGF_Gm_dd_11), &
                        G_F(iNode,iGF_Gm_dd_22), &
                        G_F(iNode,iGF_Gm_dd_33) )

                absLambda_R(iNode) = 1.0_DP

              END DO

              ! --- Numerical Flux ---

              alpha = MAX( absLambda_L, absLambda_R )

              DO iCR = 1, nCR

                NumericalFlux(:,iCR) &
                  = NumericalFlux_LLF &
                      ( uCR_L    (:,iCR), &
                        uCR_R    (:,iCR), &
                        Flux_X1_L(:,iCR), &
                        Flux_X1_R(:,iCR), alpha(:) )

                NumericalFlux(:,iCR) &
                  = dZ(3) * dZ(4) * Weights_X1(:) &
                      * G_F(:,iGF_Alpha) * Tau_X1(:) * NumericalFlux(:,iCR)

              END DO

              ! --- Contribution to this Element ---

              IF( iZ2 < iZ_E0(2) + 1 )THEN

                DO iCR = 1, nCR

                  CALL DGEMV &
                         ( 'T', nDOF_X1, nDOF, + One, L_X1_Dn, &
                           nDOF_X1, NumericalFlux(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2  ,iZ3,iZ4,iCR,iS), 1 )

                END DO

              END IF

              ! --- Contribution to Previous Element ---

              IF( iZ2 > iZ_B0(2) )THEN

                DO iCR = 1, nCR

                  CALL DGEMV &
                         ( 'T', nDOF_X1, nDOF, - One, L_X1_Up, &
                           nDOF_X1, NumericalFlux(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2-1,iZ3,iZ4,iCR,iS), 1 )

                END DO

              END IF

            END DO ! iZ1
          END DO ! iZ2
        END DO ! iZ3
      END DO ! iZ4
    END DO ! iS

  END SUBROUTINE ComputeIncrement_Divergence_X1


  SUBROUTINE ComputeIncrement_Divergence_X2 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)    :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iNode
    INTEGER  :: iGF, iCR
    REAL(DP) :: dE, dX1, dX3
    REAL(DP) :: Tau   (nDOF)
    REAL(DP) :: Tau_X2(nDOF_X2)
    REAL(DP) :: FF(nDOF)
    REAL(DP) :: EF(nDOF)
    REAL(DP) :: GX_P(nDOFX,   nGF)
    REAL(DP) :: GX_K(nDOFX,   nGF)
    REAL(DP) :: GX_F(nDOFX_X2,nGF)
    REAL(DP) :: G_K (nDOF,    nGF)
    REAL(DP) :: G_F (nDOF_X2, nGF)
    REAL(DP) :: uCR_P(nDOF,nCR)
    REAL(DP) :: uCR_K(nDOF,nCR)
    REAL(DP) :: uPR_K(nDOF,nPR)
    REAL(DP) :: Flux_X2_q(nDOF,nCR)
    REAL(DP) :: FF_L(nDOF_X2), EF_L(nDOF_X2)
    REAL(DP) :: FF_R(nDOF_X2), EF_R(nDOF_X2)
    REAL(DP) :: absLambda_L(nDOF_X2)
    REAL(DP) :: absLambda_R(nDOF_X2)
    REAL(DP) :: alpha(nDOF_X2)
    REAL(DP) :: uCR_L(nDOF_X2,nCR), uCR_R(nDOF_X2,nCR)
    REAL(DP) :: uPR_L(nDOF_X2,nPR), uPR_R(nDOF_X2,nPR)
    REAL(DP) :: Flux_X2_L(nDOF_X2,nCR)
    REAL(DP) :: Flux_X2_R(nDOF_X2,nCR)
    REAL(DP) :: NumericalFlux(nDOF_X2,nCR)

    IF( iZ_E0(3) .EQ. iZ_B0(3) ) RETURN

    DO iS = 1, nSpecies

      DO iZ4 = iZ_B0(4), iZ_E0(4)

        dX3 = MeshX(3) % Width(iZ4)

        DO iZ3 = iZ_B0(3), iZ_E0(3) + 1

          DO iZ2 = iZ_B0(2), iZ_E0(2)

            dX1 = MeshX(1) % Width(iZ2)

            ! --- Geometry Fields in Element Nodes ---

            DO iGF = 1, nGF

              GX_P(:,iGF) = GX(:,iZ2,iZ3-1,iZ4,iGF) ! --- Previous Element
              GX_K(:,iGF) = GX(:,iZ2,iZ3,  iZ4,iGF) ! --- This     Element

              G_K(1:nDOF,iGF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_K(1:nDOFX,iGF), nDOFX )

            END DO

            ! --- Interpolate Geometry Fields on Shared Face ---

            ! --- Face States (Average of Left and Right States) ---

            DO iGF = iGF_h_1, iGF_h_3

              CALL DGEMV &
                     ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                       GX_P(:,iGF), 1, Zero, GX_F(:,iGF), 1 )
              CALL DGEMV &
                     ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                       GX_K(:,iGF), 1, Half, GX_F(:,iGF), 1 )

              GX_F(1:nDOFX_X2,iGF) &
                = MAX( GX_F(1:nDOFX_X2,iGF), SqrtTiny )

              G_F(1:nDOF_X2,iGF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_F(1:nDOFX_X2,iGF), nDOFX_X2 )

            END DO

            CALL ComputeGeometryX_FromScaleFactors( G_F(:,:) )

            ! --- Lapse Function ---

            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                     GX_P(:,iGF_Alpha), 1, Zero, GX_F(:,iGF_Alpha), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                     GX_K(:,iGF_Alpha), 1, Half, GX_F(:,iGF_Alpha), 1 )

            GX_F(1:nDOFX_X2,iGF_Alpha) &
              = MAX( GX_F(1:nDOFX_X2,iGF_Alpha), SqrtTiny )

            G_F(1:nDOF_X2,iGF_Alpha) &
              = OuterProduct1D3D &
                  ( Ones(1:nDOFE), nDOFE, &
                    GX_F(1:nDOFX_X2,iGF_Alpha), nDOFX_X2 )

            DO iZ1 = iZ_B0(1), iZ_E0(1)

              dE = MeshE % Width(iZ1)

              ! --- Volume Jacobian in Energy-Position Element ---

              Tau(1:nDOF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, G_K(:,iGF_SqrtGm), nDOFX )

              Tau_X2(1:nDOF_X2) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, G_F(:,iGF_SqrtGm), nDOFX_X2 )

              DO iCR = 1, nCR

                uCR_P(:,iCR) = U(:,iZ1,iZ2,iZ3-1,iZ4,iCR,iS)
                uCR_K(:,iCR) = U(:,iZ1,iZ2,iZ3,  iZ4,iCR,iS)

              END DO

              !--------------------
              ! --- Volume Term ---
              !--------------------

              IF( iZ3 < iZ_E0(3) + 1 )THEN

                CALL ComputePrimitive_TwoMoment &
                       ( uCR_K(:,iCR_N ), uCR_K(:,iCR_G1), &
                         uCR_K(:,iCR_G2), uCR_K(:,iCR_G3), &
                         uPR_K(:,iPR_D ), uPR_K(:,iPR_I1), &
                         uPR_K(:,iPR_I2), uPR_K(:,iPR_I3), &
                         G_K(:,iGF_Gm_dd_11), &
                         G_K(:,iGF_Gm_dd_22), &
                         G_K(:,iGF_Gm_dd_33) )

                FF = FluxFactor &
                       ( uPR_K(:,iPR_D ), uPR_K(:,iPR_I1), &
                         uPR_K(:,iPR_I2), uPR_K(:,iPR_I3), &
                         G_K(:,iGF_Gm_dd_11), &
                         G_K(:,iGF_Gm_dd_22), &
                         G_K(:,iGF_Gm_dd_33) )

                EF = EddingtonFactor &
                       ( uPR_K(:,iPR_D), FF(:) )

                DO iNode = 1, nDOF

                  Flux_X2_q(iNode,1:nCR) &
                    = Flux_X2 &
                        ( uPR_K(iNode,iPR_D ), uPR_K(iNode,iPR_I1), &
                          uPR_K(iNode,iPR_I2), uPR_K(iNode,iPR_I3), &
                          FF(iNode), EF(iNode), &
                          G_K(iNode,iGF_Gm_dd_11), &
                          G_K(iNode,iGF_Gm_dd_22), &
                          G_K(iNode,iGF_Gm_dd_33) )

                END DO

                DO iCR = 1, nCR

                  Flux_X2_q(:,iCR) &
                    = dX1 * dX3 * Weights_q(:) &
                        * G_K(:,iGF_Alpha) * Tau(:) * Flux_X2_q(:,iCR)

                  CALL DGEMV &
                         ( 'T', nDOF, nDOF, One, dLdX2_q, nDOF, &
                           Flux_X2_q(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2,iZ3,iZ4,iCR,iS), 1 )

                END DO

              END IF

              !---------------------
              ! --- Surface Term ---
              !---------------------

              ! --- Interpolate Radiation Fields ---

              DO iCR = 1, nCR

                ! --- Interpolate Left State ---

                CALL DGEMV &
                       ( 'N', nDOF_X2, nDOF, One, L_X2_Up, nDOF_X2, &
                         uCR_P(:,iCR), 1, Zero, uCR_L(:,iCR), 1 )

                ! --- Interpolate Right State ---

                CALL DGEMV &
                       ( 'N', nDOF_X2, nDOF, One, L_X2_Dn, nDOF_X2, &
                         uCR_K(:,iCR), 1, Zero, uCR_R(:,iCR), 1 )

              END DO

              ! --- Left State Primitive, etc. ---

              CALL ComputePrimitive_TwoMoment &
                     ( uCR_L(:,iCR_N ), uCR_L(:,iCR_G1), &
                       uCR_L(:,iCR_G2), uCR_L(:,iCR_G3), &
                       uPR_L(:,iPR_D ), uPR_L(:,iPR_I1), &
                       uPR_L(:,iPR_I2), uPR_L(:,iPR_I3), &
                       G_F(:,iGF_Gm_dd_11), &
                       G_F(:,iGF_Gm_dd_22), &
                       G_F(:,iGF_Gm_dd_33) )

              FF_L(:) &
                = FluxFactor &
                    ( uPR_L(:,iPR_D ), uPR_L(:,iPR_I1), &
                      uPR_L(:,iPR_I2), uPR_L(:,iPR_I3), &
                      G_F(:,iGF_Gm_dd_11), &
                      G_F(:,iGF_Gm_dd_22), &
                      G_F(:,iGF_Gm_dd_33) )

              EF_L(:) &
                = EddingtonFactor &
                    ( uPR_L(:,iPR_D), FF_L(:) )

              absLambda_L(:) = One

              DO iNode = 1, nDOF_X2

                Flux_X2_L(iNode,1:nCR) &
                  = Flux_X2 &
                      ( uPR_L(iNode,iPR_D ), uPR_L(iNode,iPR_I1), &
                        uPR_L(iNode,iPR_I2), uPR_L(iNode,iPR_I3), &
                        FF_L(iNode), EF_L(iNode), &
                        G_F(iNode,iGF_Gm_dd_11), &
                        G_F(iNode,iGF_Gm_dd_22), &
                        G_F(iNode,iGF_Gm_dd_33) )

              END DO

              ! --- Right State Primitive, etc. ---

              CALL ComputePrimitive_TwoMoment &
                     ( uCR_R(:,iCR_N ), uCR_R(:,iCR_G1), &
                       uCR_R(:,iCR_G2), uCR_R(:,iCR_G3), &
                       uPR_R(:,iPR_D ), uPR_R(:,iPR_I1), &
                       uPR_R(:,iPR_I2), uPR_R(:,iPR_I3), &
                       G_F(:,iGF_Gm_dd_11), &
                       G_F(:,iGF_Gm_dd_22), &
                       G_F(:,iGF_Gm_dd_33) )

              FF_R(:) &
                = FluxFactor &
                    ( uPR_R(:,iPR_D ), uPR_R(:,iPR_I1), &
                      uPR_R(:,iPR_I2), uPR_R(:,iPR_I3), &
                      G_F(:,iGF_Gm_dd_11), &
                      G_F(:,iGF_Gm_dd_22), &
                      G_F(:,iGF_Gm_dd_33) )

              EF_R(:) &
                = EddingtonFactor &
                    ( uPR_R(:,iPR_D), FF_R(:) )

              absLambda_R(:) = One

              DO iNode = 1, nDOF_X2

                Flux_X2_R(iNode,1:nCR) &
                  = Flux_X2 &
                      ( uPR_R(iNode,iPR_D ), uPR_R(iNode,iPR_I1), &
                        uPR_R(iNode,iPR_I2), uPR_R(iNode,iPR_I3), &
                        FF_R(iNode), EF_R(iNode), &
                        G_F(iNode,iGF_Gm_dd_11), &
                        G_F(iNode,iGF_Gm_dd_22), &
                        G_F(iNode,iGF_Gm_dd_33) )

              END DO

              ! --- Numerical Flux ---

              alpha = MAX( absLambda_L, absLambda_R )

              DO iCR = 1, nCR

                NumericalFlux(:,iCR) &
                  = NumericalFlux_LLF &
                      ( uCR_L    (:,iCR), &
                        uCR_R    (:,iCR), &
                        Flux_X2_L(:,iCR), &
                        Flux_X2_R(:,iCR), alpha(:) )

                NumericalFlux(:,iCR) &
                  = dX1 * dX3 * Weights_X2(:) &
                      * G_F(:,iGF_Alpha) * Tau_X2(:) * NumericalFlux(:,iCR)

              END DO

              ! --- Contribution to this Element ---

              IF( iZ3 < iZ_E0(3) + 1 )THEN

                DO iCR = 1, nCR

                  CALL DGEMV &
                         ( 'T', nDOF_X2, nDOF, + One, L_X2_Dn, &
                           nDOF_X2, NumericalFlux(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2,iZ3  ,iZ4,iCR,iS), 1 )

                END DO

              END IF

              ! --- Contribution to Previous Element ---

              IF( iZ3 > iZ_B0(3) )THEN

                DO iCR = 1, nCR

                  CALL DGEMV &
                         ( 'T', nDOF_X2, nDOF, - One, L_X2_Up, &
                           nDOF_X2, NumericalFlux(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2,iZ3-1,iZ4,iCR,iS), 1 )

                END DO

              END IF

            END DO ! iZ1
          END DO ! iZ2
        END DO ! iZ3
      END DO ! iZ4
    END DO ! iS

  END SUBROUTINE ComputeIncrement_Divergence_X2


  SUBROUTINE ComputeIncrement_Divergence_X3 &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(in)    :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      dU(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iNode
    INTEGER  :: iGF, iCR
    REAL(DP) :: dE, dX1, dX2
    REAL(DP) :: Tau   (nDOF)
    REAL(DP) :: Tau_X3(nDOF_X3)
    REAL(DP) :: FF(nDOF)
    REAL(DP) :: EF(nDOF)
    REAL(DP) :: GX_P(nDOFX,   nGF)
    REAL(DP) :: GX_K(nDOFX,   nGF)
    REAL(DP) :: GX_F(nDOFX_X3,nGF)
    REAL(DP) :: G_K (nDOF,    nGF)
    REAL(DP) :: G_F (nDOF_X3, nGF)
    REAL(DP) :: uCR_P(nDOF,nCR)
    REAL(DP) :: uCR_K(nDOF,nCR)
    REAL(DP) :: uPR_K(nDOF,nPR)
    REAL(DP) :: Flux_X3_q(nDOF,nCR)
    REAL(DP) :: FF_L(nDOF_X3), EF_L(nDOF_X3)
    REAL(DP) :: FF_R(nDOF_X3), EF_R(nDOF_X3)
    REAL(DP) :: alpha(nDOF_X3)
    REAL(DP) :: absLambda_L(nDOF_X3)
    REAL(DP) :: absLambda_R(nDOF_X3)
    REAL(DP) :: uCR_L(nDOF_X3,nCR), uCR_R(nDOF_X3,nCR)
    REAL(DP) :: uPR_L(nDOF_X3,nPR), uPR_R(nDOF_X3,nPR)
    REAL(DP) :: Flux_X3_L(nDOF_X3,nCR)
    REAL(DP) :: Flux_X3_R(nDOF_X3,nCR)
    REAL(DP) :: NumericalFlux(nDOF_X3,nCR)

    IF( iZ_E0(4) .EQ. iZ_B0(4) ) RETURN

    DO iS = 1, nSpecies

      DO iZ4 = iZ_B0(4), iZ_E0(4) + 1

        DO iZ3 = iZ_B0(3), iZ_E0(3)

          dX2 = MeshX(2) % Width(iZ3)

          DO iZ2 = iZ_B0(2), iZ_E0(2)

            dX1 = MeshX(1) % Width(iZ2)

            ! --- Geometry Fields in Element Nodes ---

            DO iGF = 1, nGF

              GX_P(:,iGF) = GX(:,iZ2,iZ3,iZ4-1,iGF) ! --- Previous Element
              GX_K(:,iGF) = GX(:,iZ2,iZ3,iZ4,  iGF) ! --- This     Element

              G_K(1:nDOF,iGF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_K(1:nDOFX,iGF), nDOFX )

            END DO

            ! --- Interpolate Geometry Fields on Shared Face ---

            ! --- Face States (Average of Left and Right States) ---

            DO iGF = iGF_h_1, iGF_h_3

              CALL DGEMV &
                     ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                       GX_P(:,iGF), 1, Zero, GX_F(:,iGF), 1 )
              CALL DGEMV &
                     ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                       GX_K(:,iGF), 1, Half, GX_F(:,iGF), 1 )

              GX_F(1:nDOFX_X3,iGF) &
                = MAX( GX_F(1:nDOFX_X3,iGF), SqrtTiny )

              G_F(1:nDOF_X3,iGF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, GX_F(1:nDOFX_X3,iGF), nDOFX_X3 )

            END DO

            CALL ComputeGeometryX_FromScaleFactors( G_F(:,:) )

            ! --- Lapse Function ---

            CALL DGEMV &
                   ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                     GX_P(:,iGF_Alpha), 1, Zero, GX_F(:,iGF_Alpha), 1 )
            CALL DGEMV &
                   ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                     GX_K(:,iGF_Alpha), 1, Half, GX_F(:,iGF_Alpha), 1 )

            GX_F(1:nDOFX_X3,iGF_Alpha) &
              = MAX( GX_F(1:nDOFX_X3,iGF_Alpha), SqrtTiny )

            G_F(1:nDOF_X3,iGF_Alpha) &
              = OuterProduct1D3D &
                  ( Ones(1:nDOFE), nDOFE, &
                    GX_F(1:nDOFX_X3,iGF_Alpha), nDOFX_X3 )

            DO iZ1 = iZ_B0(1), iZ_E0(1)

              dE = MeshE % Width(iZ1)

              ! --- Volume Jacobian in Energy-Position Element ---

              Tau(1:nDOF) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, G_K(:,iGF_SqrtGm), nDOFX )

              Tau_X3(1:nDOF_X3) &
                = OuterProduct1D3D &
                    ( Ones(1:nDOFE), nDOFE, G_F(:,iGF_SqrtGm), nDOFX_X3 )

              DO iCR = 1, nCR

                uCR_P(:,iCR) = U(:,iZ1,iZ2,iZ3,iZ4-1,iCR,iS)
                uCR_K(:,iCR) = U(:,iZ1,iZ2,iZ3,iZ4,  iCR,iS)

              END DO

              !--------------------
              ! --- Volume Term ---
              !--------------------

              IF( iZ4 < iZ_E0(4) + 1 )THEN

                CALL ComputePrimitive_TwoMoment &
                       ( uCR_K(:,iCR_N ), uCR_K(:,iCR_G1), &
                         uCR_K(:,iCR_G2), uCR_K(:,iCR_G3), &
                         uPR_K(:,iPR_D ), uPR_K(:,iPR_I1), &
                         uPR_K(:,iPR_I2), uPR_K(:,iPR_I3), &
                         G_K(:,iGF_Gm_dd_11), &
                         G_K(:,iGF_Gm_dd_22), &
                         G_K(:,iGF_Gm_dd_33) )

                FF = FluxFactor &
                       ( uPR_K(:,iPR_D ), uPR_K(:,iPR_I1), &
                         uPR_K(:,iPR_I2), uPR_K(:,iPR_I3), &
                         G_K(:,iGF_Gm_dd_11), &
                         G_K(:,iGF_Gm_dd_22), &
                         G_K(:,iGF_Gm_dd_33) )

                EF = EddingtonFactor &
                       ( uPR_K(:,iPR_D), FF(:) )

                DO iNode = 1, nDOF

                  Flux_X3_q(iNode,1:nCR) &
                    = Flux_X3 &
                        ( uPR_K(iNode,iPR_D ), uPR_K(iNode,iPR_I1), &
                          uPR_K(iNode,iPR_I2), uPR_K(iNode,iPR_I3), &
                          FF(iNode), EF(iNode), &
                          G_K(iNode,iGF_Gm_dd_11), &
                          G_K(iNode,iGF_Gm_dd_22), &
                          G_K(iNode,iGF_Gm_dd_33) )

                END DO

                DO iCR = 1, nCR

                  Flux_X3_q(:,iCR) &
                    = dX1 * dX2 * Weights_q(:) &
                        * G_K(:,iGF_Alpha) * Tau(:) * Flux_X3_q(:,iCR)

                  CALL DGEMV &
                         ( 'T', nDOF, nDOF, One, dLdX3_q, nDOF, &
                           Flux_X3_q(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2,iZ3,iZ4,iCR,iS), 1 )

                END DO

              END IF

              !---------------------
              ! --- Surface Term ---
              !---------------------

              ! --- Interpolate Radiation Fields ---

              DO iCR = 1, nCR

                ! --- Interpolate Left State ---

                CALL DGEMV &
                       ( 'N', nDOF_X3, nDOF, One, L_X3_Up, nDOF_X3, &
                         uCR_P(:,iCR), 1, Zero, uCR_L(:,iCR), 1 )

                ! --- Interpolate Right State ---

                CALL DGEMV &
                       ( 'N', nDOF_X3, nDOF, One, L_X3_Dn, nDOF_X3, &
                         uCR_K(:,iCR), 1, Zero, uCR_R(:,iCR), 1 )

              END DO

              ! --- Left State Primitive, etc. ---

              CALL ComputePrimitive_TwoMoment &
                     ( uCR_L(:,iCR_N ), uCR_L(:,iCR_G1), &
                       uCR_L(:,iCR_G2), uCR_L(:,iCR_G3), &
                       uPR_L(:,iPR_D ), uPR_L(:,iPR_I1), &
                       uPR_L(:,iPR_I2), uPR_L(:,iPR_I3), &
                       G_F(:,iGF_Gm_dd_11), &
                       G_F(:,iGF_Gm_dd_22), &
                       G_F(:,iGF_Gm_dd_33) )

              FF_L(:) &
                = FluxFactor &
                    ( uPR_L(:,iPR_D ), uPR_L(:,iPR_I1), &
                      uPR_L(:,iPR_I2), uPR_L(:,iPR_I3), &
                      G_F(:,iGF_Gm_dd_11), &
                      G_F(:,iGF_Gm_dd_22), &
                      G_F(:,iGF_Gm_dd_33) )

              EF_L(:) &
                = EddingtonFactor &
                    ( uPR_L(:,iPR_D), FF_L(:) )

              absLambda_L(:) = One

              DO iNode = 1, nDOF_X3

                Flux_X3_L(iNode,1:nCR) &
                  = Flux_X3 &
                      ( uPR_L(iNode,iPR_D ), uPR_L(iNode,iPR_I1), &
                        uPR_L(iNode,iPR_I2), uPR_L(iNode,iPR_I3), &
                        FF_L(iNode), EF_L(iNode), &
                        G_F(iNode,iGF_Gm_dd_11), &
                        G_F(iNode,iGF_Gm_dd_22), &
                        G_F(iNode,iGF_Gm_dd_33) )

              END DO

              ! --- Right State Primitive, etc. ---

              CALL ComputePrimitive_TwoMoment &
                     ( uCR_R(:,iCR_N ), uCR_R(:,iCR_G1), &
                       uCR_R(:,iCR_G2), uCR_R(:,iCR_G3), &
                       uPR_R(:,iPR_D ), uPR_R(:,iPR_I1), &
                       uPR_R(:,iPR_I2), uPR_R(:,iPR_I3), &
                       G_F(:,iGF_Gm_dd_11), &
                       G_F(:,iGF_Gm_dd_22), &
                       G_F(:,iGF_Gm_dd_33) )

              FF_R(:) &
                = FluxFactor &
                    ( uPR_R(:,iPR_D ), uPR_R(:,iPR_I1), &
                      uPR_R(:,iPR_I2), uPR_R(:,iPR_I3), &
                      G_F(:,iGF_Gm_dd_11), &
                      G_F(:,iGF_Gm_dd_22), &
                      G_F(:,iGF_Gm_dd_33) )

              EF_R(:) &
                = EddingtonFactor &
                    ( uPR_R(:,iPR_D), FF_R(:) )

              absLambda_R(:) = One

              DO iNode = 1, nDOF_X3

                Flux_X3_R(iNode,1:nCR) &
                  = Flux_X3 &
                      ( uPR_R(iNode,iPR_D ), uPR_R(iNode,iPR_I1), &
                        uPR_R(iNode,iPR_I2), uPR_R(iNode,iPR_I3), &
                        FF_R(iNode), EF_R(iNode), &
                        G_F(iNode,iGF_Gm_dd_11), &
                        G_F(iNode,iGF_Gm_dd_22), &
                        G_F(iNode,iGF_Gm_dd_33) )

              END DO

              ! --- Numerical Flux ---

              alpha = MAX( absLambda_L, absLambda_R )

              DO iCR = 1, nCR

                NumericalFlux(:,iCR) &
                  = NumericalFlux_LLF &
                      ( uCR_L    (:,iCR), &
                        uCR_R    (:,iCR), &
                        Flux_X3_L(:,iCR), &
                        Flux_X3_R(:,iCR), alpha(:) )

                NumericalFlux(:,iCR) &
                  = dX1 * dX2 * Weights_X3(:) &
                      * G_F(:,iGF_Alpha) * Tau_X3(:) * NumericalFlux(:,iCR)

              END DO

              ! --- Contribution to this Element ---

              IF( iZ4 < iZ_E0(4) + 1 )THEN

                DO iCR = 1, nCR

                  CALL DGEMV &
                         ( 'T', nDOF_X3, nDOF, + One, L_X3_Dn, &
                           nDOF_X3, NumericalFlux(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2,iZ3  ,iZ4,iCR,iS), 1 )

                END DO

              END IF

              ! --- Contribution to Previous Element ---

              IF( iZ4 > iZ_B0(4) )THEN

                DO iCR = 1, nCR

                  CALL DGEMV &
                         ( 'T', nDOF_X3, nDOF, - One, L_X3_Up, &
                           nDOF_X3, NumericalFlux(:,iCR), 1, One, &
                           dU(:,iZ1,iZ2,iZ3,iZ4-1,iCR,iS), 1 )

                END DO

              END IF

            END DO ! iZ1
          END DO ! iZ2
        END DO ! iZ3
      END DO ! iZ4
    END DO ! iS

  END SUBROUTINE ComputeIncrement_Divergence_X3


  SUBROUTINE ComputeIncrement_Geometry &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(in)    :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(inout) :: &
      dU(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iGF
    INTEGER  :: iNodeZ, iNodeX
    REAL(DP) :: PR_K(nDOF,nPR), FF(nDOF), EF(nDOF), Stress(3)
    REAL(DP) :: &
      h2_X1(nDOFX_X1,iZ_B0(2):iZ_E0(2)+1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
      h3_X1(nDOFX_X1,iZ_B0(2):iZ_E0(2)+1,iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
      h3_X2(nDOFX_X2,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3)+1,iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      dh2dX1(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
      dh3dX1(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
      dh3dX2(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      G(nDOF,nGF,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    IF( TRIM( CoordinateSystem ) == 'CARTESIAN' ) RETURN

    ! --- Derivative of Scale Factor h_2 wrt X1 ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX(:,iZ2-1,iZ3,iZ4,iGF_h_2), 1, Zero, h2_X1(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX(:,iZ2  ,iZ3,iZ4,iGF_h_2), 1, Half, h2_X1(:,iZ2,iZ3,iZ4), 1 )

      ! --- h_2 on X1 Faces ---

      h2_X1(:,iZ2,iZ3,iZ4) &
        = WeightsX_X1(:) * MAX( h2_X1(:,iZ2,iZ3,iZ4), SqrtTiny )

    END DO
    END DO
    END DO

    ASSOCIATE( dZ2 => MeshX(1) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      CALL DGEMV &
             ( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
               h2_X1(:,iZ2+1,iZ3,iZ4), 1, Zero, dh2dX1(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
               h2_X1(:,iZ2  ,iZ3,iZ4), 1,  One, dh2dX1(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX, WeightsX_q &
               * GX(:,iZ2,iZ3,iZ4,iGF_h_2), 1, One, dh2dX1(:,iZ2,iZ3,iZ4), 1 )

      dh2dX1(:,iZ2,iZ3,iZ4) &
        = dh2dX1(:,iZ2,iZ3,iZ4) / ( WeightsX_q(:) * dZ2(iZ2) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dZ2, etc.

    ! --- Derivative of Scale Factor h_3 wrt X1 ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2) + 1

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
               GX(:,iZ2-1,iZ3,iZ4,iGF_h_3), 1, Zero, h3_X1(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
               GX(:,iZ2  ,iZ3,iZ4,iGF_h_3), 1, Half, h3_X1(:,iZ2,iZ3,iZ4), 1 )

      ! --- h_3 on X1 Faces ---

      h3_X1(:,iZ2,iZ3,iZ4) &
        = WeightsX_X1(:) * MAX( h3_X1(:,iZ2,iZ3,iZ4), SqrtTiny )

    END DO
    END DO
    END DO

    ASSOCIATE( dZ2 => MeshX(1) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      CALL DGEMV &
             ( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
               h3_X1(:,iZ2+1,iZ3,iZ4), 1, Zero, dh3dX1(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
               h3_X1(:,iZ2  ,iZ3,iZ4), 1,  One, dh3dX1(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX, WeightsX_q &
               * GX(:,iZ2,iZ3,iZ4,iGF_h_3), 1, One, dh3dX1(:,iZ2,iZ3,iZ4), 1 )

      dh3dX1(:,iZ2,iZ3,iZ4) &
        = dh3dX1(:,iZ2,iZ3,iZ4) / ( WeightsX_q(:) * dZ2(iZ2) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dZ2, etc.

    ! --- Derivative of Scale Factor h_3 wrt X2 ---

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3) + 1
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
               GX(:,iZ2,iZ3-1,iZ4,iGF_h_3), 1, Zero, h3_X2(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
               GX(:,iZ2,iZ3  ,iZ4,iGF_h_3), 1, Half, h3_X2(:,iZ2,iZ3,iZ4), 1 )

      ! --- h_3 on X2 Faces ---

      h3_X2(:,iZ2,iZ3,iZ4) &
        = WeightsX_X2(:) * MAX( h3_X2(:,iZ2,iZ3,iZ4), SqrtTiny )

    END DO
    END DO
    END DO

    ASSOCIATE( dZ3 => MeshX(2) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      CALL DGEMV &
             ( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
               h3_X2(:,iZ2,iZ3+1,iZ4), 1, Zero, dh3dX2(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
               h3_X2(:,iZ2,iZ3  ,iZ4), 1,  One, dh3dX2(:,iZ2,iZ3,iZ4), 1 )
      CALL DGEMV &
             ( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX, WeightsX_q &
               * GX(:,iZ2,iZ3,iZ4,iGF_h_3), 1, One, dh3dX2(:,iZ2,iZ3,iZ4), 1 )

      dh3dX2(:,iZ2,iZ3,iZ4) &
        = dh3dX2(:,iZ2,iZ3,iZ4) / ( WeightsX_q(:) * dZ3(iZ3) )

    END DO
    END DO
    END DO

    END ASSOCIATE ! dZ2, etc.

    DO iGF = 1, nGF
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeZ = 1, nDOF

        G(iNodeZ,iGF,iZ2,iZ3,iZ4) &
          = GX(NodeNumbersX(iNodeZ),iZ2,iZ3,iZ4,iGF)

      END DO

    END DO
    END DO
    END DO
    END DO

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      CALL ComputePrimitive_TwoMoment &
             ( U(:,iZ1,iZ2,iZ3,iZ4,iCR_N, iS), &
               U(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
               U(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
               U(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
               PR_K(:,iPR_D ), PR_K(:,iPR_I1), &
               PR_K(:,iPR_I2), PR_K(:,iPR_I3), &
               G(:,iGF_Gm_dd_11,iZ2,iZ3,iZ4),  &
               G(:,iGF_Gm_dd_22,iZ2,iZ3,iZ4),  &
               G(:,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

      FF = FluxFactor &
             ( PR_K(:,iPR_D ), PR_K(:,iPR_I1), &
               PR_K(:,iPR_I2), PR_K(:,iPR_I3), &
               G(:,iGF_Gm_dd_11,iZ2,iZ3,iZ4),  &
               G(:,iGF_Gm_dd_22,iZ2,iZ3,iZ4),  &
               G(:,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

      EF = EddingtonFactor( PR_K(:,iPR_D), FF )

      DO iNodeZ = 1, nDOF

        iNodeX = NodeNumbersX(iNodeZ)

        Stress = StressTensor_Diagonal &
                   ( PR_K(iNodeZ,iPR_D ), PR_K(iNodeZ,iPR_I1), &
                     PR_K(iNodeZ,iPR_I2), PR_K(iNodeZ,iPR_I3), &
                     FF(iNodeZ), EF(iNodeZ), &
                     G(iNodeZ,iGF_Gm_dd_11,iZ2,iZ3,iZ4), &
                     G(iNodeZ,iGF_Gm_dd_22,iZ2,iZ3,iZ4),  &
                     G(iNodeZ,iGF_Gm_dd_33,iZ2,iZ3,iZ4) )

        ! --- Add to Increments ---

        dU(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) &
          = dU(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) &
              + Stress(2) * dh2dX1(iNodeX,iZ2,iZ3,iZ4) &
                  / G(iNodeZ,iGF_h_2,iZ2,iZ3,iZ4) &
              + Stress(3) * dh3dX1(iNodeX,iZ2,iZ3,iZ4) &
                  / G(iNodeZ,iGF_h_3,iZ2,iZ3,iZ4)

        dU(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) &
          = dU(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) &
              + Stress(3) * dh3dX2(iNodeX,iZ2,iZ3,iZ4) &
                  / G(iNodeZ,iGF_h_3,iZ2,iZ3,iZ4)

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

#ifdef THORNADO_DEBUG_EXPLICIT
    WRITE(*,'(a,4i4,es23.15)') 'MINLOC(dh2dX1), MINVAL(dh2dX1)', MINLOC(dh2dX1), MINVAL(dh2dX1)
    WRITE(*,'(a,4i4,es23.15)') 'MAXLOC(dh2dX1), MAXVAL(dh2dX1)', MAXLOC(dh2dX1), MAXVAL(dh2dX1)
    WRITE(*,'(a,4i4,es23.15)') 'MINLOC(dh3dX1), MINVAL(dh3dX1)', MINLOC(dh3dX1), MINVAL(dh3dX1)
    WRITE(*,'(a,4i4,es23.15)') 'MAXLOC(dh3dX1), MAXVAL(dh3dX1)', MAXLOC(dh3dX1), MAXVAL(dh3dX1)
    WRITE(*,'(a,4i4,es23.15)') 'MINLOC(dh3dX2), MINVAL(dh3dX2)', MINLOC(dh3dX2), MINVAL(dh3dX2)
    WRITE(*,'(a,4i4,es23.15)') 'MAXLOC(dh3dX2), MAXVAL(dh3dX2)', MAXLOC(dh3dX2), MAXVAL(dh3dX2)
#endif

  END SUBROUTINE ComputeIncrement_Geometry


END MODULE TwoMoment_DiscretizationModule_Streaming
