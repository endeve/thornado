MODULE PolynomialBasisMappingModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX
  USE QuadratureModule, ONLY: &
    xG5, wG5
  USE UtilitiesModule, ONLY: &
    WriteMatrix
  USE PolynomialBasisModule_Lagrange, ONLY: &
    L_X1, L_X2, L_X3
  USE PolynomialBasisModule_Legendre, ONLY: &
    P_X1, P_X2, P_X3, MassPX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePolynomialBasisMapping
  PUBLIC :: MapNodalToModal
  PUBLIC :: MapModalToNodal

  INTERFACE MapNodalToModal
    MODULE PROCEDURE MapNodalToModal_Fluid
    MODULE PROCEDURE MapNodalToModal_Radiation
  END INTERFACE MapNodalToModal

  INTERFACE MapModalToNodal
    MODULE PROCEDURE MapModalToNodal_Fluid
    MODULE PROCEDURE MapModalToNodal_Radiation
  END INTERFACE MapModalToNodal

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Kij_X
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Pij_X

CONTAINS


  SUBROUTINE InitializePolynomialBasisMapping( Nodes_X1, Nodes_X2, Nodes_X3 )

    REAL(DP), DIMENSION(:), INTENT(in) :: Nodes_X1, Nodes_X2, Nodes_X3

    INTEGER :: iX1, iX2, iX3
    INTEGER :: jX1, jX2, jX3
    INTEGER :: qX1, qX2, qX3
    INTEGER :: i, j

    ALLOCATE( Kij_X(1:nDOFX,1:nDOFX), Pij_X(1:nDOFX,1:nDOFX) )

    i = 0
    DO iX1 = 1, nNodesX(1)
      DO iX2 = 1, nNodesX(2)
        DO iX3 = 1, nNodesX(3)

          i = i + 1

          j = 0
          DO jX1 = 1, nNodesX(1)
            DO jX2 = 1, nNodesX(2)
              DO jX3 = 1, nNodesX(3)

                j = j + 1

                Kij_X(i,j) = 0.0_DP
                DO qX3 = 1, SIZE( xG5 )
                  DO qX2 = 1, SIZE( xG5 )
                    DO qX1 = 1, SIZE( xG5 )

                      Kij_X(i,j) &
                        = Kij_X(i,j) &
                            + wG5(qX1) * wG5(qX2) * wG5(qX3) &
                                * P_X1(iX1) % P( xG5(qX1) ) &
                                    * P_X2(iX2) % P( xG5(qX2) ) &
                                        * P_X3(iX3) % P( xG5(qX3) ) &
                                * L_X1(jX1) % P( xG5(qX1) ) &
                                    * L_X2(jX2) % P( xG5(qX2) ) &
                                        * L_X3(jX3) % P( xG5(qX3) )

                    END DO
                  END DO
                END DO

                Pij_X(i,j) &
                  = P_X1(jX1) % P( Nodes_X1(iX1) ) &
                      * P_X2(jX2) % P( Nodes_X2(iX2) ) &
                          * P_X3(jX3) % P( Nodes_X3(iX3) )

              END DO
            END DO
          END DO

        END DO
      END DO
    END DO

  END SUBROUTINE InitializePolynomialBasisMapping


  SUBROUTINE MapNodalToModal_Fluid( uN, uM )

    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uN
    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uM

    INTEGER :: i

    uM = 0.0_DP
    DO i = 1, nDOFX
      uM(:) = uM(:) + Kij_X(:,i) * uN(i)
    END DO
    uM = MassPX * uM

  END SUBROUTINE MapNodalToModal_Fluid


  SUBROUTINE MapNodalToModal_Radiation

  END SUBROUTINE MapNodalToModal_Radiation


  SUBROUTINE MapModalToNodal_Fluid( uN, uM )

    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uN
    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uM

    INTEGER :: i

    uN = 0.0_DP
    DO i = 1, nDOFX
      uN(:) = uN(:) + Pij_X(:,i) * uM(i)
    END DO

  END SUBROUTINE MapModalToNodal_Fluid


  SUBROUTINE MapModalToNodal_Radiation

  END SUBROUTINE MapModalToNodal_Radiation


END MODULE PolynomialBasisMappingModule
