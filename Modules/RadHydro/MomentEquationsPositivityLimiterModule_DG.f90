MODULE MomentEquationsPositivityLimiterModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOF, &
    nX, nNodesX, &
    nE, nNodesE
  USE PolynomialBasisModule_Lagrange, ONLY: &
    evalL
  USE PolynomialBasisModule_Legendre, ONLY: &
    evalP
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Radiation, &
    MapModalToNodal_Radiation
  USE GeometryFieldsModule, ONLY: &
    Vol, VolJac
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    WeightsR, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, nCR
  USE MomentEquationsLimiterUtilitiesModule_DG, ONLY: &
    SolveTheta_Bisection
  USE MomentEquationsLimiterModule_DG, ONLY: &
    ApplyPositivityLimiter, &
    nPositivePoints, &
    Points_E, Points_X1, Points_X2, Points_X3, &
    Lagrange

  IMPLICIT NONE
  PRIVATE

  LOGICAL,  PARAMETER :: Debug = .FALSE.
  REAL(DP), PARAMETER :: Tol_N = SQRT(TINY(1.0_DP))
  REAL(DP), PARAMETER :: Tol_G = SQRT(TINY(1.0_DP))

  PUBLIC :: ApplyPositivityLimiter_M1_DG

CONTAINS


  SUBROUTINE ApplyPositivityLimiter_M1_DG

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iCR, iPoint
    REAL(DP) :: Theta_1, Theta_2, Theta_P
    REAL(DP) :: dN, dG1, dG2, dG3
    REAL(DP), DIMENSION(nPositivePoints)     :: absG_P
    REAL(DP), DIMENSION(1:nCR)               :: uCR_K
    REAL(DP), DIMENSION(nDOF,nCR)            :: uCR_M
    REAL(DP), DIMENSION(nPositivePoints,nCR) :: uCR_P

    IF( nDOF == 1 ) RETURN

    IF( .NOT. ApplyPositivityLimiter ) RETURN

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iCR = 1, nCR
                DO iPoint = 1, nPositivePoints

                  uCR_P(iPoint,iCR) &
                    = DOT_PRODUCT &
                        ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), Lagrange(:,iPoint) )

                END DO
              END DO

              absG_P(:) &
                = SQRT( uCR_P(:,iCR_G1)**2 &
                        + uCR_P(:,iCR_G2)**2 &
                        + uCR_P(:,iCR_G3)**2 )

              IF( ANY( uCR_P(:,iCR_N) < Tol_N ) .OR. &
                  ANY( uCR_P(:,iCR_N) - absG_P(:) < Tol_G ) )THEN

                ! --- Limiting Using Modal Representation ---

                DO iCR = 1, nCR

                  CALL MapNodalToModal_Radiation &
                         ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), uCR_M(:,iCR) )

                  uCR_K(iCR) &
                    = SUM( WeightsR(:) * uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                           * VolJac(:,iE,iX1,iX2,iX3) ) / Vol(iE,iX1,iX2,iX3)

                END DO

                ! --- Ensure Positive Number Density ---

                Theta_1 = 1.0_DP
                DO iPoint = 1, nPositivePoints
                  IF( uCR_P(iPoint,iCR_N) < Tol_N ) &
                    Theta_1 &
                      = MIN( Theta_1, &
                             ( Tol_N - uCR_K(iCR_N) ) &
                               / ( uCR_P(iPoint,iCR_N) - uCR_K(iCR_N) ) )
                END DO

                IF( Theta_1 < 1.0_DP )THEN

                  IF( Debug )THEN
                    PRINT*
                    PRINT*, "iE,iX1,iX2,iX3 = ", iE, iX1, iX2, iX3
                    PRINT*, "Theta_1 = ", Theta_1
                    PRINT*, " N_K = ", uCR_K(iCR_N)
                    PRINT*, " N_p = ", uCR_P(:,iCR_N)
                    PRINT*
                  END IF

                  uCR_M(1,iCR_N) &
                    = Theta_1 * uCR_M(1,iCR_N) &
                        + (1.0_DP - Theta_1) * uCR_K(iCR_N)

                  uCR_M(2:nDOF,iCR_N) &
                    = Theta_1 * uCR_M(2:nDOF,iCR_N)

                  DO iPoint = 1, nPositivePoints

                    uCR_P(iPoint,iCR_N) &
                      = evalP( uCR_M(:,iCR_N), &
                               Points_E (iPoint), Points_X1(iPoint), &
                               Points_X2(iPoint), Points_X3(iPoint) )

                  END DO

                END IF

                ! --- Ensure Limited Number Flux ---

                absG_P(:) &
                  = SQRT( uCR_P(:,iCR_G1)**2 &
                          + uCR_P(:,iCR_G2)**2 &
                          + uCR_P(:,iCR_G3)**2 )

                Theta_2 = 1.0_DP
                DO iPoint = 1, nPositivePoints

                  IF( uCR_P(iPoint,iCR_N) - absG_P(iPoint) < Tol_G )THEN

                    dN  = uCR_P(iPoint,iCR_N)  - uCR_K(iCR_N)
                    dG1 = uCR_P(iPoint,iCR_G1) - uCR_K(iCR_G1)
                    dG2 = uCR_P(iPoint,iCR_G2) - uCR_K(iCR_G2)
                    dG3 = uCR_P(iPoint,iCR_G3) - uCR_K(iCR_G3)

                    CALL SolveTheta_Bisection &
                           ( uCR_K(iCR_N), uCR_K(iCR_G1), uCR_K(iCR_G2), &
                             uCR_K(iCR_G3), dN, dG1, dG2, dG3, &
                             Tol_G, Theta_P )

                    Theta_2 = MIN( Theta_2, Theta_P )

                  END IF

                END DO

                IF( Theta_2 < 1.0_DP )THEN

                  IF( Debug )THEN
                    PRINT*
                    PRINT*, "iE,iX1,iX2,iX3 = ", iE, iX1, iX2, iX3
                    PRINT*, "Theta_2 = ", Theta_2
                    PRINT*, "Tol_G   = ", Tol_G
                    PRINT*, " G_K = ", uCR_K(iCR_N)-SQRT( uCR_K(iCR_G1)**2 &
                                                          + uCR_K(iCR_G2)**2 &
                                                          + uCR_K(iCR_G3)**2 )
                    PRINT*, " G_p = ", uCR_P(:,iCR_N)-absG_P(:)
                    PRINT*
                  END IF

                  DO iCR = 1, nCR

                    uCR_M(1,iCR) &
                      = Theta_2 * uCR_M(1,iCR) &
                          + (1.0_DP - Theta_2) * uCR_K(iCR)

                    uCR_M(2:nDOF,iCR) &
                      = Theta_2 * uCR_M(2:nDOF,iCR)

                  END DO

                END IF

                ! --- Back to Nodal Representation ---

                DO iCR = 1, nCR

                  CALL MapModalToNodal_Radiation &
                         ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), uCR_M(:,iCR) )

                END DO

                IF( Debug )THEN

                  DO iCR = 1, nCR
                    DO iPoint = 1, nPositivePoints

                      uCR_P(iPoint,iCR) &
                        = evalL( uCR(:,iE,iX1,iX2,iX3,iCR,iS), &
                                 Points_E (iPoint), Points_X1(iPoint), &
                                 Points_X2(iPoint), Points_X3(iPoint) )

                    END DO

                    uCR_K(iCR) &
                      = SUM( WeightsR(:) * uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                           * VolJac(:,iE,iX1,iX2,iX3) ) / Vol(iE,iX1,iX2,iX3)

                  END DO

                  absG_P &
                    = SQRT( uCR_P(:,iCR_G1)**2 + uCR_P(:,iCR_G2)**2 &
                            + uCR_P(:,iCR_G3)**2 )

                  IF( ANY( uCR_P(:,iCR_N) < Tol_N ) .OR. &
                      ANY( uCR_P(:,iCR_N) - absG_P(:) < Tol_G ) )THEN

                    PRINT*
                    PRINT*, "ApplyPositivityLimiter_M1_DG"
                    PRINT*
                    PRINT*, "  Problem with Positivity Limiter!"
                    PRINT*
                    PRINT*, "    iE, iX1, iX2, iX3 = ", iE, iX1, iX2, iX3
                    PRINT*, "    Theta_1 = ", Theta_1
                    PRINT*, "    Theta_2 = ", Theta_2
                    PRINT*, "    Tol_N   = ", Tol_N
                    PRINT*, "    Tol_G   = ", Tol_G
                    PRINT*
                    PRINT*, "  Conserved Radiation Fields (Nodal):"
                    PRINT*, "  N     = ", uCR_P(:,iCR_N)
                    PRINT*, "  G1    = ", uCR_P(:,iCR_G1)
                    PRINT*, "  G2    = ", uCR_P(:,iCR_G2)
                    PRINT*, "  G3    = ", uCR_P(:,iCR_G3)
                    PRINT*, "  N-|G| = ", &
                      uCR_P(:,iCR_N) &
                      - SQRT( uCR_P(:,iCR_G1)**2 &
                              + uCR_P(:,iCR_G2)**2 &
                              + uCR_P(:,iCR_G3)**2 )
                    PRINT*
                    PRINT*, "  Cell-Averages:"
                    PRINT*, "  N_K       = ", uCR_K(iCR_N)
                    PRINT*, "  G1_K      = ", uCR_K(iCR_G1)
                    PRINT*, "  G2_K      = ", uCR_K(iCR_G2)
                    PRINT*, "  G3_K      = ", uCR_K(iCR_G3)
                    PRINT*, "  N_K-|G_K| = ", &
                      uCR_K(iCR_N) &
                      - SQRT( uCR_K(iCR_G1)**2 &
                              + uCR_K(iCR_G2)**2 &
                              + uCR_K(iCR_G3)**2 )
                    PRINT*

                  END IF

                END IF ! Debug

              END IF

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ApplyPositivityLimiter_M1_DG


END MODULE MomentEquationsPositivityLimiterModule_DG
