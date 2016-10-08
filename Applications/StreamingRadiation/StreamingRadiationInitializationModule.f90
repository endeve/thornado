MODULE StreamingRadiationInitializationModule

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeStreamingSineWave1D
  PUBLIC :: InitializeGaussianSphericalWave1D
  PUBLIC :: InitializeLineSource1D

CONTAINS


  SUBROUTINE InitializeStreamingSineWave1D( Epsilon_Option, Delta_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Epsilon_Option, Delta_Option

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iNode
    REAL(DP) :: X1, Epsilon, Delta

    Epsilon = 0.0_DP
    IF( PRESENT( Epsilon_Option ) ) &
      Epsilon = Epsilon_Option

    Delta = 0.0_DP
    IF( PRESENT( Delta_Option ) ) &
      Delta = Delta_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A5,A10,ES10.4E2)') '', 'Epsilon = ', Epsilon
    WRITE(*,'(A5,A10,ES10.4E2)') '', 'Delta = ', Delta
    WRITE(*,*)

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3) 
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iNodeX3 = 1, nNodesX(3)
                DO iNodeX2 = 1, nNodesX(2)
                  DO iNodeX1 = 1, nNodesX(1)
                    DO iNodeE = 1, nNodesE

                      X1    = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                      iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                        = 1.0_DP + ( 1.0_DP - Epsilon ) &
                                     * SIN( 2.0_DP * Pi * X1 )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                        = ( 1.0_DP - Delta ) &
                            * uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS)

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                        = 0.0_DP

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                        = 0.0_DP

                    END DO
                  END DO
                END DO
              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE InitializeStreamingSineWave1D


  SUBROUTINE InitializeGaussianSphericalWave1D

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iNode
    REAL(DP) :: X1

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3) 
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iNodeX3 = 1, nNodesX(3)
                DO iNodeX2 = 1, nNodesX(2)
                  DO iNodeX1 = 1, nNodesX(1)
                    DO iNodeE = 1, nNodesE

                      X1    = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                      iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                        = EXP( - X1**2 ) / X1**2

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                        = uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS)

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                        = 0.0_DP

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                        = 0.0_DP

                    END DO
                  END DO
                END DO
              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE InitializeGaussianSphericalWave1D


  SUBROUTINE InitializeLineSource1D( Sigma_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Sigma_Option

    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iNode
    REAL(DP) :: X1, Sigma

    Sigma = 0.03_DP
    IF( PRESENT( Sigma_Option ) ) &
      Sigma = Sigma_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A5,A8,ES10.4E2)') '', 'Sigma = ', Sigma
    WRITE(*,*)

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iNodeX3 = 1, nNodesX(3)
                DO iNodeX2 = 1, nNodesX(2)
                  DO iNodeX1 = 1, nNodesX(1)
                    DO iNodeE = 1, nNodesE

                      X1    = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                      iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_N,iS) &
                        = MAX( EXP( - 0.5_DP * ( X1 / Sigma )**2 ) &
                               / ( 2.0 * Sigma**2 ), 1.0d-4 )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,iS) &
                        = EPSILON( 1.0_DP )

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G2,iS) &
                        = 0.0_DP

                      uCR(iNode,iE,iX1,iX2,iX3,iCR_G3,iS) &
                        = 0.0_DP

                    END DO
                  END DO
                END DO
              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE InitializeLineSource1D


END MODULE StreamingRadiationInitializationModule
