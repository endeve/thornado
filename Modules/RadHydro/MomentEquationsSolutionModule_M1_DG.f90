MODULE MomentEquationsSolutionModule_M1_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, dL_X1
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem, &
    a, b
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    rhsCR, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, nCR
  USE RiemannSolverModule, ONLY: &
    NumericalFlux_Radiation
  USE MomentEquationsUtilitiesModule, ONLY: &
    AlphaMax, &
    AlphaP, &
    AlphaM, &
    Flux_X1, &
    GeometrySources
  

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeRHS_M1_DG

CONTAINS


  SUBROUTINE ComputeRHS_M1_DG( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    CALL ComputeRHS_M1_DG_X1( iX_Begin, iX_End )

    CALL ComputeRHS_M1_DG_GeometrySources( iX_Begin, iX_End )

  END SUBROUTINE ComputeRHS_M1_DG


  SUBROUTINE ComputeRHS_M1_DG_X1( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER :: iS, iX1, iX2, iX3, iE
    INTEGER :: iNodeX1, jNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER :: iNode, jNode, iCR
    REAL(DP) :: Alpha, AlphaPls, AlphaMns
    REAL(DP), DIMENSION(1:nCR) :: VolumeTerm, Flux_L, Flux_R, Flux
    REAL(DP), DIMENSION(1:nCR) :: uCR_L, uCR_R
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
        dX1 => MeshX(1) % Width (1:nX(1)) )

    ! -- Precompute Lagrange Polynomials --

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

    ! -- Precompute Metric Functions --

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

    ! -- Compute Right-Hand Side for Moment Equations --

    DO iS = 1, nSpecies

      DO iX3 = iX_Begin(3), iX_End(3)
        DO iX2 = iX_Begin(2), iX_End(2)
          DO iX1 = iX_Begin(1), iX_End(1)
            DO iE = 1, nE

              ASSOCIATE &
                ( uCR_P => uCR(:,iE,iX1-1,iX2,iX3,:,iS), & ! Previous Element
                  uCR_K => uCR(:,iE,iX1,  iX2,iX3,:,iS), & ! This     Element
                  uCR_N => uCR(:,iE,iX1+1,iX2,iX3,:,iS) )  ! Next     Element

              DO iNodeX3 = 1, nNodesX(3)
                DO iNodeX2 = 1, nNodesX(2)
                  DO iNodeX1 = 1, nNodesX(1)
                    DO iNodeE = 1, nNodesE

                      iNode &
                        = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                      rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) = 0.0_DP

                      ! -- Volume Term --

                      VolumeTerm = 0.0_DP
                      DO jNodeX1 = 1, nNodesX(1)

                        jNode &
                          = NodeNumber( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                        VolumeTerm(1:nCR) &
                          = VolumeTerm(1:nCR) &
                              + w_q(jNodeX1) &
                                  * a_X1_q(jNodeX1,iX1) &
                                  * b_X1_q(jNodeX1,iX1) &
                                  * Flux_X1( uCR_K(jNode,iCR_N), &
                                             uCR_K(jNode,iCR_G1), &
                                             uCR_K(jNode,iCR_G2), &
                                             uCR_K(jNode,iCR_G3) ) &
                                  * dL_X1_q(jNodeX1,iNodeX1)

                      END DO

                      rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                        = rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                            + VolumeTerm(1:nCR) &
                                / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                      * b_X1_q(iNodeX1,iX1) * dX1(iX1) )

                      ! -- Left Face --

                      ! -- Left State --

                      DO iCR = 1, nCR
                        uCR_L(iCR) = 0.0_DP
                        DO jNodeX1 = 1, nNodesX(1)

                          jNode &
                            = NodeNumber( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                          uCR_L(iCR) &
                            = uCR_L(iCR) &
                                + L_X1_R(jNodeX1) * uCR_P(jNode,iCR)

                        END DO
                      END DO

                      Flux_L &
                        = Flux_X1( uCR_L(iCR_N),  uCR_L(iCR_G1), &
                                   uCR_L(iCR_G2), uCR_L(iCR_G3) )

                      ! -- Right State --

                      DO iCR = 1, nCR
                        uCR_R(iCR) = 0.0_DP
                        DO jNodeX1 = 1, nNodesX(1)

                          jNode &
                            = NodeNumber( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                          uCR_R(iCR) &
                            = uCR_R(iCR) &
                                + L_X1_L(jNodeX1) * uCR_K(jNode,iCR)

                        END DO
                      END DO

                      Flux_R &
                        = Flux_X1( uCR_R(iCR_N),  uCR_R(iCR_G1), &
                                   uCR_R(iCR_G2), uCR_R(iCR_G3) )

                      ! -- Numerical Flux --

                      Alpha &
                        = MAX( AlphaMax( uCR_L(iCR_N) , uCR_L(iCR_G1)  , &
                                         uCR_L(iCR_G2), uCR_L(iCR_G3) ), &
                               AlphaMax( uCR_R(iCR_N) , uCR_R(iCR_G1)  , &
                                         uCR_R(iCR_G2), uCR_R(iCR_G3) ) )

                      AlphaPls &
                        = AlphaP( uCR_L(iCR_N) , uCR_L(iCR_G1), &
                                  uCR_L(iCR_G2), uCR_L(iCR_G3), &
                                  uCR_R(iCR_N) , uCR_R(iCR_G1), &
                                  uCR_R(iCR_G2), uCR_R(iCR_G3) )

                      AlphaMns &
                        = AlphaM( uCR_L(iCR_N) , uCR_L(iCR_G1), &
                                  uCR_L(iCR_G2), uCR_L(iCR_G3), &
                                  uCR_R(iCR_N) , uCR_R(iCR_G1), &
                                  uCR_R(iCR_G2), uCR_R(iCR_G3) )

                      Flux = NumericalFlux_Radiation &
                               ( uCR_L, uCR_R, Flux_L, Flux_R, &
                                 Alpha, AlphaPls, AlphaMns, 1.0_DP, nCR )

                      ! -- Contribution to Right-Hand Side --

                      rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                        = rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                            + a_X1_L(iX1) * b_X1_L(iX1) &
                                * Flux(1:nCR) * L_X1_L(iNodeX1) &
                                    / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                          * b_X1_q(iNodeX1,iX1) * dX1(iX1) )

                      ! -- Right Face --

                      ! -- Left State --

                      DO iCR = 1, nCR
                        uCR_L(iCR) = 0.0_DP
                        DO jNodeX1 = 1, nNodesX(1)

                          jNode &
                            = NodeNumber( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                          uCR_L(iCR) &
                            = uCR_L(iCR) &
                                + L_X1_R(jNodeX1) * uCR_K(jNode,iCR)

                        END DO
                      END DO

                      Flux_L &
                        = Flux_X1( uCR_L(iCR_N),  uCR_L(iCR_G1), &
                                   uCR_L(iCR_G2), uCR_L(iCR_G3) )

                      ! -- Right State --

                      DO iCR = 1, nCR
                        uCR_R(iCR) = 0.0_DP
                        DO jNodeX1 = 1, nNodesX(1)

                          jNode &
                            = NodeNumber( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                          uCR_R(iCR) &
                            = uCR_R(iCR) &
                                + L_X1_L(jNodeX1) * uCR_N(jNode,iCR)

                        END DO
                      END DO

                      Flux_R &
                        = Flux_X1( uCR_R(iCR_N),  uCR_R(iCR_G1), &
                                   uCR_R(iCR_G2), uCR_R(iCR_G3) )

                      ! -- Numerical Flux --

                      Alpha &
                        = MAX( AlphaMax( uCR_L(iCR_N),  uCR_L(iCR_G1),   &
                                         uCR_L(iCR_G2), uCR_L(iCR_G3) ), &
                               AlphaMax( uCR_R(iCR_N),  uCR_R(iCR_G1),   &
                                         uCR_R(iCR_G2), uCR_R(iCR_G3) ) )

                      AlphaPls &
                        = AlphaP( uCR_L(iCR_N) , uCR_L(iCR_G1), &
                                  uCR_L(iCR_G2), uCR_L(iCR_G3), &
                                  uCR_R(iCR_N) , uCR_R(iCR_G1), &
                                  uCR_R(iCR_G2), uCR_R(iCR_G3) )

                      AlphaMns &
                        = AlphaM( uCR_L(iCR_N) , uCR_L(iCR_G1), &
                                  uCR_L(iCR_G2), uCR_L(iCR_G3), &
                                  uCR_R(iCR_N) , uCR_R(iCR_G1), &
                                  uCR_R(iCR_G2), uCR_R(iCR_G3) )

                      Flux = NumericalFlux_Radiation &
                               ( uCR_L, uCR_R, Flux_L, Flux_R, &
                                 Alpha, AlphaPls, AlphaMns, 1.0_DP, nCR )

                      ! -- Contribution to Right-Hand Side --

                      rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                        = rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                            -  a_X1_R(iX1) * b_X1_R(iX1) &
                                 * Flux(1:nCR) * L_X1_R(iNodeX1) &
                                     / ( w_q(iNodeX1) * a_X1_q(iNodeX1,iX1) &
                                           * b_X1_q(iNodeX1,iX1) * dX1(iX1) )

                    END DO ! iNodeE
                  END DO ! iNodeX1
                END DO ! iNodeX2
              END DO ! iNodeX3

              END ASSOCIATE ! uCR_P, etc.

            END DO ! iE
          END DO ! iX1
        END DO ! iX2
      END DO ! iX3

    END DO ! iS

    END ASSOCIATE ! x_q, etc.

  END SUBROUTINE ComputeRHS_M1_DG_X1


  SUBROUTINE ComputeRHS_M1_DG_GeometrySources( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER  :: iNode
    REAL(DP) :: X1, X2, X3

    IF( TRIM( CoordinateSystem ) == 'CARTSEIAN' ) RETURN

    DO iS = 1, nSpecies

      DO iX3 = iX_Begin(3), iX_End(3)
        DO iX2 = iX_Begin(2), iX_End(2)
          DO iX1 = iX_Begin(1), iX_End(1)
            DO iE = 1, nE

              ASSOCIATE &
                ( uCR_K => uCR(:,iE,iX1,iX2,iX3,:,iS) ) ! This Element

              DO iNodeX3 = 1, nNodesX(3)

                X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

                DO iNodeX2 = 1, nNodesX(2)

                  X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

                  DO iNodeX1 = 1, nNodesX(1)

                    X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                    DO iNodeE = 1, nNodesE

                      iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                      rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                        = rhsCR(iNode,iE,iX1,iX2,iX3,1:nCR,iS) &
                            + GeometrySources &
                                ( uCR_K(iNode,iCR_N),  uCR_K(iNode,iCR_G1), &
                                  uCR_K(iNode,iCR_G2), uCR_K(iNode,iCR_G3), &
                                  [ X1, X2, X3 ] )

                    END DO
                  END DO
                END DO
              END DO

              END ASSOCIATE ! uCR_K

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ComputeRHS_M1_DG_GeometrySources


END MODULE MomentEquationsSolutionModule_M1_DG
