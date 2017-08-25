MODULE PolynomialBasisMappingModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX, nDOF
  USE QuadratureModule, ONLY: &
    xG5, wG5
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE PolynomialBasisModule_Lagrange, ONLY: &
    IndL_Q, IndLX_Q, L_E, L_X1, L_X2, L_X3
  USE PolynomialBasisModule_Legendre, ONLY: &
    IndP_Q, INDPX_Q, P_E, P_X1, P_X2, P_X3, &
    MassP, MassPX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePolynomialBasisMapping
  PUBLIC :: MapNodalToModal_Fluid
  PUBLIC :: MapNodalToModal_Radiation_X
  PUBLIC :: MapNodalToModal_Radiation
  PUBLIC :: MapModalToNodal_Fluid
  PUBLIC :: MapModalToNodal_Radiation_X
  PUBLIC :: MapModalToNodal_Radiation

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Kij_X, K_ij
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Pij_X, P_ij

CONTAINS


  SUBROUTINE InitializePolynomialBasisMapping &
               ( Nodes_E, Nodes_X1, Nodes_X2, Nodes_X3 )

    REAL(DP), DIMENSION(:), INTENT(in) :: Nodes_E, Nodes_X1, Nodes_X2, Nodes_X3

    INTEGER :: i, j
    INTEGER :: qE, qX1, qX2, qX3

    ALLOCATE( Kij_X(1:nDOFX,1:nDOFX), Pij_X(1:nDOFX,1:nDOFX) )

    Kij_X = 0.0_DP
    DO j = 1, nDOFX
      DO i = 1, nDOFX

        DO qX3 = 1, SIZE( xG5 )
          DO qX2 = 1, SIZE( xG5 )
            DO qX1 = 1, SIZE( xG5 )

              Kij_X(i,j) &
                = Kij_X(i,j) &
                    + wG5(qX1) * wG5(qX2) * wG5(qX3) &
                        * P_X1(IndPX_Q(1,i)) % P( xG5(qX1) ) &
                            * P_X2(IndPX_Q(2,i)) % P( xG5(qX2) ) &
                                * P_X3(IndPX_Q(3,i)) % P( xG5(qX3) ) &
                        * L_X1(IndLX_Q(1,j)) % P( xG5(qX1) ) &
                            * L_X2(IndLX_Q(2,j)) % P( xG5(qX2) ) &
                                * L_X3(IndLX_Q(3,j)) % P( xG5(qX3) )

            END DO
          END DO
        END DO

        Pij_X(i,j) &
          = P_X1(IndPX_Q(1,j)) % P( Nodes_X1(IndLX_Q(1,i)) ) &
              * P_X2(IndPX_Q(2,j)) % P( Nodes_X2(IndLX_Q(2,i)) ) &
                  * P_X3(IndPX_Q(3,j)) % P( Nodes_X3(IndLX_Q(3,i)) )

      END DO
    END DO

    ALLOCATE( K_ij(1:nDOF,1:nDOF), P_ij(1:nDOF,1:nDOF) )

    K_ij = 0.0_DP
    DO j = 1, nDOF
      DO i = 1, nDOF

        DO qX3 = 1, SIZE( xG5 )
          DO qX2 = 1, SIZE( xG5 )
            DO qX1 = 1, SIZE( xG5 )
              DO qE  = 1, SIZE( xG5 )

                K_ij(i,j) &
                  = K_ij(i,j) &
                      + wG5(qE) * wG5(qX1) * wG5(qX2) * wG5(qX3) &
                          * P_E(IndP_Q(0,i)) % P( xG5(qE) ) &
                              * P_X1(IndP_Q(1,i)) % P( xG5(qX1) ) &
                                  * P_X2(IndP_Q(2,i)) % P( xG5(qX2) ) &
                                      * P_X3(IndP_Q(3,i)) % P( xG5(qX3) ) &
                          * L_E(IndL_Q(0,j)) % P( xG5(qE) ) &
                              * L_X1(IndL_Q(1,j)) % P( xG5(qX1) ) &
                                  * L_X2(IndL_Q(2,j)) % P( xG5(qX2) ) &
                                      * L_X3(IndL_Q(3,j)) % P( xG5(qX3) )

              END DO
            END DO
          END DO
        END DO

        P_ij(i,j) &
          = P_E(IndP_Q(0,j)) % P( Nodes_E(IndL_Q(0,i)) ) &
              * P_X1(IndP_Q(1,j)) % P( Nodes_X1(IndL_Q(1,i)) ) &
                  * P_X2(IndP_Q(2,j)) % P( Nodes_X2(IndL_Q(2,i)) ) &
                      * P_X3(IndP_Q(3,j)) % P( Nodes_X3(IndL_Q(3,i)) )

      END DO
    END DO

  END SUBROUTINE InitializePolynomialBasisMapping


  SUBROUTINE MapNodalToModal_Fluid( uN, uM )

    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uN
    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uM

    uM = MassPX * MATMUL( Kij_X, uN )

  END SUBROUTINE MapNodalToModal_Fluid


  SUBROUTINE MapNodalToModal_Radiation_X( uN, uM )

    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uN
    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uM

    uM = MassPX * MATMUL( Kij_X, uN )

  END SUBROUTINE MapNodalToModal_Radiation_X


  SUBROUTINE MapNodalToModal_Radiation( uN, uM )

    REAL(DP), DIMENSION(nDOF), INTENT(in)  :: uN
    REAL(DP), DIMENSION(nDOF), INTENT(out) :: uM

    uM = MassP * MATMUL( K_ij, uN )

  END SUBROUTINE MapNodalToModal_Radiation


  SUBROUTINE MapModalToNodal_Fluid( uN, uM )

    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uN
    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uM

    uN = MATMUL( Pij_X, uM )

  END SUBROUTINE MapModalToNodal_Fluid


  SUBROUTINE MapModalToNodal_Radiation_X( uN, uM )

    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uN
    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uM

    uN = MATMUL( Pij_X, uM )

  END SUBROUTINE MapModalToNodal_Radiation_X


  SUBROUTINE MapModalToNodal_Radiation( uN, uM )

    REAL(DP), DIMENSION(nDOF), INTENT(out) :: uN
    REAL(DP), DIMENSION(nDOF), INTENT(in)  :: uM

    uN = MATMUL( P_ij, uM )

  END SUBROUTINE MapModalToNodal_Radiation


END MODULE PolynomialBasisMappingModule
