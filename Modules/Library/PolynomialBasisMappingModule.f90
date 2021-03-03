MODULE PolynomialBasisMappingModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOFZ
  USE QuadratureModule, ONLY: &
    xG5, wG5
  USE PolynomialBasisModule_Lagrange, ONLY: &
    IndL_Q, IndLX_Q, L_E, L_X1, L_X2, L_X3
  USE PolynomialBasisModule_Legendre, ONLY: &
    IndP_Q, IndPX_Q, P_E, P_X1, P_X2, P_X3, &
    MassP, MassPX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePolynomialBasisMapping
  PUBLIC :: MapNodalToModalX
  PUBLIC :: MapNodalToModalZ
  PUBLIC :: MapNodalToModal_Fluid
  PUBLIC :: MapNodalToModal_Radiation_X
  PUBLIC :: MapNodalToModal_Radiation
  PUBLIC :: MapModalToNodalX
  PUBLIC :: MapModalToNodalZ
  PUBLIC :: MapModalToNodal_Fluid
  PUBLIC :: MapModalToNodal_Radiation_X
  PUBLIC :: MapModalToNodal_Radiation

  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: Kij_X, Pij_X
  REAL(DP), DIMENSION(:,:), ALLOCATABLE         :: Kij_Z, Pij_Z


CONTAINS


  SUBROUTINE InitializePolynomialBasisMapping &
    ( Nodes_E, Nodes_X1, Nodes_X2, Nodes_X3 )

    REAL(DP), INTENT(in) :: Nodes_E (:)
    REAL(DP), INTENT(in) :: Nodes_X1(:)
    REAL(DP), INTENT(in) :: Nodes_X2(:)
    REAL(DP), INTENT(in) :: Nodes_X3(:)

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
        =   P_X1(IndPX_Q(1,j)) % P( Nodes_X1(IndLX_Q(1,i)) ) &
          * P_X2(IndPX_Q(2,j)) % P( Nodes_X2(IndLX_Q(2,i)) ) &
          * P_X3(IndPX_Q(3,j)) % P( Nodes_X3(IndLX_Q(3,i)) )

    END DO
    END DO

    ALLOCATE( Kij_Z(1:nDOFZ,1:nDOFZ), Pij_Z(1:nDOFZ,1:nDOFZ) )

    Kij_Z = 0.0_DP
    DO j = 1, nDOFZ
    DO i = 1, nDOFZ

      DO qX3 = 1, SIZE( xG5 )
      DO qX2 = 1, SIZE( xG5 )
      DO qX1 = 1, SIZE( xG5 )
      DO qE  = 1, SIZE( xG5 )

        Kij_Z(i,j) &
          = Kij_Z(i,j) &
              + wG5(qE) * wG5(qX1) * wG5(qX2) * wG5(qX3) &
                  * P_E (IndP_Q(0,i)) % P( xG5(qE ) ) &
                  * P_X1(IndP_Q(1,i)) % P( xG5(qX1) ) &
                  * P_X2(IndP_Q(2,i)) % P( xG5(qX2) ) &
                  * P_X3(IndP_Q(3,i)) % P( xG5(qX3) ) &
                  * L_E (IndL_Q(0,j)) % P( xG5(qE ) ) &
                  * L_X1(IndL_Q(1,j)) % P( xG5(qX1) ) &
                  * L_X2(IndL_Q(2,j)) % P( xG5(qX2) ) &
                  * L_X3(IndL_Q(3,j)) % P( xG5(qX3) )

      END DO
      END DO
      END DO
      END DO

      Pij_Z(i,j) &
        =   P_E (IndP_Q(0,j)) % P( Nodes_E (IndL_Q(0,i)) ) &
          * P_X1(IndP_Q(1,j)) % P( Nodes_X1(IndL_Q(1,i)) ) &
          * P_X2(IndP_Q(2,j)) % P( Nodes_X2(IndL_Q(2,i)) ) &
          * P_X3(IndP_Q(3,j)) % P( Nodes_X3(IndL_Q(3,i)) )

    END DO
    END DO

  END SUBROUTINE InitializePolynomialBasisMapping


  ! --- Nodal to Modal ---


  SUBROUTINE MapNodalToModalX( uN, uM )

    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uN
    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uM

    uM = MassPX * MATMUL( Kij_X, uN )

  END SUBROUTINE MapNodalToModalX


  SUBROUTINE MapNodalToModalZ( uN, uM )

    REAL(DP), DIMENSION(nDOFZ), INTENT(in)  :: uN
    REAL(DP), DIMENSION(nDOFZ), INTENT(out) :: uM

    uM = MassP * MATMUL( Kij_Z, uN )

  END SUBROUTINE MapNodalToModalZ


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

    REAL(DP), DIMENSION(nDOFZ), INTENT(in)  :: uN
    REAL(DP), DIMENSION(nDOFZ), INTENT(out) :: uM

    uM = MassP * MATMUL( Kij_Z, uN )

  END SUBROUTINE MapNodalToModal_Radiation


  ! --- Modal to Nodal ---


  SUBROUTINE MapModalToNodalX( uN, uM )

    REAL(DP), DIMENSION(nDOFX), INTENT(out) :: uN
    REAL(DP), DIMENSION(nDOFX), INTENT(in)  :: uM

    uN = MATMUL( Pij_X, uM )

  END SUBROUTINE MapModalToNodalX


  SUBROUTINE MapModalToNodalZ( uN, uM )

    REAL(DP), DIMENSION(nDOFZ), INTENT(out) :: uN
    REAL(DP), DIMENSION(nDOFZ), INTENT(in)  :: uM

    uN = MATMUL( Pij_Z, uM )

  END SUBROUTINE MapModalToNodalZ


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

    REAL(DP), DIMENSION(nDOFZ), INTENT(out) :: uN
    REAL(DP), DIMENSION(nDOFZ), INTENT(in)  :: uM

    uN = MATMUL( Pij_Z, uM )

  END SUBROUTINE MapModalToNodal_Radiation


END MODULE PolynomialBasisMappingModule
