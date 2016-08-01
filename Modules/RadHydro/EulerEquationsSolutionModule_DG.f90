MODULE EulerEquationsSolutionModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, dL_X1
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryModule, ONLY: &
    CoordinateSystem, &
    a, b
  USE FluidFieldsModule, ONLY: &
    rhsCF, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, nCF, iAF_P, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, nPF, &
    uAF, iAF_P, iAF_Gm, iAF_Cs, nAF
  USE EquationOfStateModule, ONLY: &
    ComputeAuxiliary_Fluid, &
    Auxiliary_Fluid
  USE RiemannSolverModule, ONLY: &
    NumericalFlux_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputePrimitive, &
    Primitive, &
    AlphaMax, &
    AlphaP, &
    AlphaM, &
    AlphaC, &
    Flux_X1, &
    GeometrySources

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ComputeRHS_Euler_DG

CONTAINS


  SUBROUTINE ComputeRHS_Euler_DG( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    CALL ComputePrimitive( iX_Begin, iX_End )

    CALL ComputeAuxiliary_Fluid( iX_Begin, iX_End )

    CALL ComputeRHS_Euler_DG_X1( iX_Begin, iX_End )

    CALL ComputeRHS_Euler_DG_GeometrySources(iX_Begin, iX_End )

  END SUBROUTINE ComputeRHS_Euler_DG


  SUBROUTINE ComputeRHS_Euler_DG_X1( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iX1, iX2, iX3, iCF
    INTEGER :: iNodeX1, jNodeX1, iNodeX2, iNodeX3, iNodeX, jNodeX
    REAL(DP) :: Alpha, AlphaPls, AlphaMns, AlphaMdl
    REAL(DP), DIMENSION(1:nCF) :: VolumeTerm, Flux_L, Flux_R, Flux
    REAL(DP), DIMENSION(1:nCF) :: uCF_L, uCF_R
    REAL(DP), DIMENSION(1:nPF) :: uPF_L, uPF_R
    REAL(DP), DIMENSION(1:nAF) :: uAF_L, uAF_R
    REAL(DP), DIMENSION(nX(1)) :: a_X1_L, a_X1_R
    REAL(DP), DIMENSION(nX(1)) :: b_X1_L, b_X1_R
    REAL(DP), DIMENSION(nNodesX(1)) :: L_X1_L, L_X1_R
    REAL(DP), DIMENSION(nNodesX(1), nX(1)) :: a_X1_q
    REAL(DP), DIMENSION(nNodesX(1), nX(1)) :: b_X1_q
    REAL(DP), DIMENSION(nNodesX(1),nNodesX(1)) :: dL_X1_q

    ASSOCIATE &
      ( x_q => MeshX(1) % Nodes, &
        w_q => MeshX(1) % Weights, &
        X1C => MeshX(1) % Center(1:nX(1)), &
        dX1 => MeshX(1) % Width(1:nX(1)) )

    DO jNodeX1 = 1, nNodesX(1)
      L_X1_L(jNodeX1) &
        = L_X1(jNodeX1) % P( - 0.5_DP )
      L_X1_R(jNodeX1) &
        = L_X1(jNodeX1) % P( + 0.5_DP )
      DO iNodeX1 = 1, nNodesX(1)
        dL_X1_q(iNodeX1,jNodeX1) &
          = dL_X1(jNodeX1) % P( x_q(iNodeX1) )
      END DO
    END DO

    DO iX1 = iX_Begin(1), iX_End(1)
      a_X1_L(iX1) &
        = a( [ X1C(iX1) - 0.5_DP * dX1(iX1), 0.0_DP, 0.0_DP ] )
      a_X1_R(iX1) &
        = a( [ X1C(iX1) + 0.5_DP * dX1(iX1), 0.0_DP, 0.0_DP ] )
      b_X1_L(iX1) &
        = b( [ X1C(iX1) - 0.5_DP * dX1(iX1), 0.0_DP, 0.0_DP ] )
      b_X1_R(iX1) &
        = b( [ X1C(iX1) + 0.5_DP * dX1(iX1), 0.0_DP, 0.0_DP ] )
      DO iNodeX1 = 1, nNodesX(1)
        a_X1_q(iNodeX1,iX1) &
          = a( [ X1C(iX1) + dX1(iX1) * x_q(iNodeX1), 0.0_DP, 0.0_DP ] )
        b_X1_q(iNodeX1,iX1) &
          = b( [ X1C(iX1) + dX1(iX1) * x_q(iNodeX1), 0.0_DP, 0.0_DP ] )
      END DO
    END DO

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          ASSOCIATE &
            ( uCF_P => uCF(:,iX1-1,iX2,iX3,:), & ! Previous Element
              uCF_K => uCF(:,iX1  ,iX2,iX3,:), & ! This     Element
              uCF_N => uCF(:,iX1+1,iX2,iX3,:), & ! Next     Element
              uPF_K => uPF(:,iX1  ,iX2,iX3,:), & ! This     Element
              uAF_K => uAF(:,iX1  ,iX2,iX3,:) )  ! This     Element

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) = 0.0_DP

                ! -- Volume Term --

                VolumeTerm = 0.0_DP
                DO jNodeX1 = 1, nNodesX(1)

                  jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

                  VolumeTerm(1:nCF) &
                    = VolumeTerm(1:nCF) &
                        + w_q(jNodeX1) &
                            * dL_X1_q(jNodeX1,iNodeX1) &
                            * a_X1_q(jNodeX1, iX1) &
                            * b_X1_q(jNodeX1, iX1) &
                            * Flux_X1( uPF_K(jNodeX,iPF_D ), &
                                       uPF_K(jNodeX,iPF_V1), &
                                       uPF_K(jNodeX,iPF_V2), &
                                       uPF_K(jNodeX,iPF_V3), &
                                       uPF_K(jNodeX,iPF_E ), &
                                       uAF_K(jNodeX,iAF_P ) )  
                END DO

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      + VolumeTerm(1:nCF) / ( w_q(iNodeX1) &
                                                * a_X1_q(iNodeX1, iX1) &
                                                * b_X1_q(iNodeX1, iX1) & 
                                                * dX1(iX1) )
                ! -- Left Face -- 

                ! -- Left State -- 

                DO iCF = 1, nCF

                  uCF_L(iCF) = 0.0_DP
                  DO jNodeX1 = 1, nNodesX(1)

                    jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

                    uCF_L(iCF) &
                      = uCF_L(iCF) &
                          + L_X1_R(jNodeX1) * uCF_P(jNodeX,iCF)

                  END DO

                END DO

                uPF_L = Primitive( uCF_L )

                uAF_L = Auxiliary_Fluid( uPF_L )

                Flux_L &
                  = Flux_X1( uPF_L(iPF_D), uPF_L(iPF_V1), uPF_L(iPF_V2), &
                             uPF_L(iPF_V3), uPF_L(iPF_E), uAF_L(iAF_P) )

                ! -- Right State --

                DO iCF = 1, nCF

                  uCF_R(iCF) = 0.0_DP
                  DO jNodeX1 = 1, nNodesX(1)

                    jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

                    uCF_R(iCF) &
                      = uCF_R(iCF) &
                          + L_X1_L(jNodeX1) * uCF_K(jNodeX,iCF)

                  END DO

                END DO

                uPF_R = Primitive( uCF_R )

                uAF_R = Auxiliary_Fluid( uPF_R )

                Flux_R &
                  = Flux_X1( uPF_R(iPF_D), uPF_R(iPF_V1), uPF_R(iPF_V2), &
                             uPF_R(iPF_V3), uPF_R(iPF_E), uAF_R(iAF_P) )

                ! -- Numerical Flux --

                Alpha &
                  = MAX( AlphaMax( uPF_L(iPF_V1), uAF_L(iAF_Cs) ), &
                         AlphaMax( uPF_R(iPF_V1), uAF_R(iAF_Cs) ) )

                AlphaPls &
                  = AlphaP( uPF_L(iPF_V1), uAF_L(iAF_Cs), &
                            uPF_R(iPF_V1), uAF_R(iAF_Cs) )

                AlphaMns &
                  = AlphaM( uPF_L(iPF_V1), uAF_L(iAF_Cs), &
                            uPF_R(iPF_V1), uAF_R(iAF_Cs) )

                AlphaMdl &
                  = AlphaC &
                      ( [ uCF_L (iCF_D), uCF_L (iCF_S1) ], &
                        [ uCF_R (iCF_D), uCF_R (iCF_S1) ], &
                        [ Flux_L(iCF_D), Flux_L(iCF_S1) ], &
                        [ Flux_R(iCF_D), Flux_R(iCF_S1) ], &
                        AlphaPls, AlphaMns )

                Flux = NumericalFlux_Fluid &
                         ( uCF_L, uCF_R, Flux_L, Flux_R, &
                           Alpha, AlphaPls, AlphaMns, AlphaMdl, nCF )

                ! -- Contribution to Right-Hand Side --

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      + L_X1_L(iNodeX1) &
                          * Flux(1:nCF) &
                          * a_X1_L(iX1) &
                          * b_X1_L(iX1) / ( w_q(iNodeX1) * dX1(iX1) & 
                                            * a_X1_q(iNodeX1, iX1) &
                                            * b_X1_q(iNodeX1, iX1) )

                ! -- Right Face --

                ! -- Left State --

                DO iCF = 1, nCF

                  uCF_L(iCF) = 0.0_DP
                  DO jNodeX1 = 1, nNodesX(1)

                    jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

                    uCF_L(iCF) &
                      = uCF_L(iCF) &
                          + L_X1_R(jNodeX1) * uCF_K(jNodeX,iCF)

                  END DO

                END DO

                uPF_L = Primitive( uCF_L )

                uAF_L = Auxiliary_Fluid( uPF_L )

                Flux_L &
                  = Flux_X1( uPF_L(iPF_D), uPF_L(iPF_V1), uPF_L(iPF_V2), &
                             uPF_L(iPF_V3), uPF_L(iPF_E), uAF_L(iAF_P) )

                ! -- Right State --

                DO iCF = 1, nCF

                  uCF_R(iCF) = 0.0_DP
                  DO jNodeX1 = 1, nNodesX(1)

                    jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

                    uCF_R(iCF) &
                      = uCF_R(iCF) &
                          + L_X1_L(jNodeX1) * uCF_N(jNodeX,iCF)

                  END DO

                END DO

                uPF_R = Primitive( uCF_R )

                uAF_R = Auxiliary_Fluid( uPF_R )

                Flux_R &
                  = Flux_X1( uPF_R(iPF_D), uPF_R(iPF_V1), uPF_R(iPF_V2), &
                             uPF_R(iPF_V3), uPF_R(iPF_E), uAF_R(iAF_P) )

                ! -- Numerical Flux --

                Alpha &
                  = MAX( AlphaMax( uPF_L(iPF_V1), uAF_L(iAF_Cs) ), &
                         AlphaMax( uPF_R(iPF_V1), uAF_R(iAF_Cs) ) )

                AlphaPls &
                  = AlphaP( uPF_L(iPF_V1), uAF_L(iAF_Cs), &
                            uPF_R(iPF_V1), uAF_R(iAF_Cs) )

                AlphaMns &
                  = AlphaM( uPF_L(iPF_V1), uAF_L(iAF_Cs), &
                            uPF_R(iPF_V1), uAF_R(iAF_Cs) )

                AlphaMdl &
                  = AlphaC &
                      ( [ uCF_L (iCF_D), uCF_L (iCF_S1) ], &
                        [ uCF_R (iCF_D), uCF_R (iCF_S1) ], &
                        [ Flux_L(iCF_D), Flux_L(iCF_S1) ], &
                        [ Flux_R(iCF_D), Flux_R(iCF_S1) ], &
                        AlphaPls, AlphaMns )

                Flux = NumericalFlux_Fluid &
                         ( uCF_L, uCF_R, Flux_L, Flux_R, &
                           Alpha, AlphaPls, AlphaMns, AlphaMdl, nCF )

                ! -- Contribution to Right-Hand Side --

                rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                  = rhsCF(iNodeX,iX1,iX2,iX3,1:nCF) &
                      - L_X1_R(iNodeX1) * Flux(1:nCF) &
                          * a_X1_R(iX1) * b_X1_R(iX1) &
                          / ( w_q(iNodeX1) * a_X1_q(iNodeX1, iX1) &
                                * b_X1_q(iNodeX1, iX1) * dX1(iX1))
              END DO
            END DO
          END DO

          END ASSOCIATE ! uCF_P, etc.

        END DO
      END DO
    END DO

    END ASSOCIATE ! x_q, etc.

  END SUBROUTINE ComputeRHS_Euler_DG_X1

  SUBROUTINE ComputeRHS_Euler_DG_GeometrySources(iX_Begin, iX_End)

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iX1, iX2, iX3
    INTEGER :: iNodeX1, iNodeX2, iNodeX3
    INTEGER :: iNode
    REAL(DP) :: X1, X2, X3

    IF( TRIM( CoordinateSystem ) == 'CARTESIAN' ) RETURN

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)
          ASSOCIATE &
            (uCF_K => uCF(:,iX1, iX2, iX3, :))
          DO iNodeX3 = 1, nNodesX(3)
            X3 = NodeCoordinate(MeshX(3), iX3, iNodeX3)
            DO iNodeX2 = 1, nNodesX(2)
              X2 = NodeCoordinate(MeshX(2), iX2, iNodeX2)
              DO iNodeX1 = 1, nNodesX(1)
                X1 = NodeCoordinate(MeshX(1), iX1, iNodeX1)
                iNode = NodeNumberX(iNodeX1, iNodeX2, iNodeX3)
                rhsCF(iNode, iX1, iX2, iX3, 1:nCF) &
                = rhsCF(iNode, iX1, iX2, iX3, 1:nCF) &
                + GeometrySources(uCF_K(iNode, iCF_D), &
                                  uCF_K(iNode, iCF_S1), &
                                  uCF_K(iNode, iCF_S2), &
                                  uCF_K(iNode, iCF_S3), &
                                  uCF_K(iNode, iAF_P), &
                                  [X1, X2, X3 ] )
              END DO
            END DO
          END DO
          END ASSOCIATE 
        END DO
      END DO 
    END DO
  END SUBROUTINE ComputeRHS_Euler_DG_GeometrySources

END MODULE EulerEquationsSolutionModule_DG
