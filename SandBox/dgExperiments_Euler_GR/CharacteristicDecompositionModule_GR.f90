MODULE CharacteristicDecompositionModule_GR

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One, Two, Four
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputePrimitive_GR, &
    Eigenvalues_GR
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive_GR

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: ComputeCharacteristicDecomposition_GR

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_GR ( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF), U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF), invR(nCF,nCF)

    REAL(DP) :: D(1), V1(1), V2(1), V3(1), E(1), Ne(1), P(1), Cs(1)
    REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33
    REAL(DP) :: W, h, LapseFunction
    REAL(DP) :: ShiftVector(3)
    REAL(DP) :: epsilon, kappa, chi, VSq, VdSq, K
    REAL(DP) :: EigVals(nCF)
    REAL(DP) :: LAMBDA_m, LAMBDA_p
    REAL(DP) :: Vm, Vp, Am, Ap, Nm, Np

    LOGICAL  :: DEBUG = .FALSE.
    REAL(DP) :: dFdU(nCF,nCF), LAMBDA(nCF,nCF)
    INTEGER  :: i

    Gmdd11         = G(iGF_Gm_dd_11)
    Gmdd22         = G(iGF_Gm_dd_22)
    Gmdd33         = G(iGF_Gm_dd_33)
    LapseFunction  = G(iGF_Alpha)
    ShiftVector(1) = G(iGF_Beta_1)
    ShiftVector(2) = G(iGF_Beta_2)
    ShiftVector(3) = G(iGF_Beta_3)

    CALL ComputePrimitive_GR &
           ( [U(iCF_D )], [U(iCF_S1)], [U(iCF_S2)], &
             [U(iCF_S3)], [U(iCF_E )], [U(iCF_Ne)], &
             D, V1, V2, V3, E, Ne, P, &
             [G(iGF_Gm_dd_11)], &
             [G(iGF_Gm_dd_22)], &
             [G(iGF_Gm_dd_33)] )

    Vu1 = V1(1)
    Vu2 = V2(1)
    Vu3 = V3(1)

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    VdSq = Vd1**2 + Vd2**2 + Vd3**2

    VSq = Vu1*Vd1 + Vu2*Vd2 + Vu3*Vd3

    ! --- Auxiliary quantities ---
    W = One / SQRT( One - VSq )
    h = One + ( E(1) + P(1) ) / D(1)

    epsilon = P(1) / ( ( Gamma_IDEAL - One ) * D(1) )
    chi     = ( Gamma_IDEAL - One ) * epsilon
    kappa   = ( Gamma_IDEAL - One ) * D(1)

    CALL ComputeSoundSpeedFromPrimitive_GR( D, E, Ne, Cs )

    ! --- Rezzolla, Eq. (7.244) ---
    K = ( Gamma_IDEAL - One ) / ( Gamma_IDEAL - One - Cs(1)**2 )

    SELECT CASE( iDim )

      CASE( 1 )

        EigVals = Eigenvalues_GR &
                    ( Vu1, Vu2, Vu3, Vu1, Cs(1), &
                      Gmdd11, Gmdd22, Gmdd33, Gmdd11, &
                      LapseFunction, ShiftVector(1) )

        ! --- lambda_m = EigVals(1) ---
        ! --- lambda_p = EigVals(3) ---

        ! --- Rezzolla, Eq. (7.245) ---
        LAMBDA_m = ( EigVals(1) + ShiftVector(1) ) / LapseFunction
        LAMBDA_p = ( EigVals(3) + ShiftVector(1) ) / LapseFunction

        Vm = ( Vu1 - LAMBDA_m ) / ( One / Gmdd11 - Vu1 * LAMBDA_m )
        Vp = ( Vu1 - LAMBDA_p ) / ( One / Gmdd11 - Vu1 * LAMBDA_p )

        Am = ( One / Gmdd11 - Vu1 * Vu1 ) &
             / ( One / Gmdd11 - Vu1 * LAMBDA_m )
        Ap = ( One / Gmdd11 - Vu1 * Vu1 ) &
             / ( One / Gmdd11 - Vu1 * LAMBDA_p )

        R(:,1) = [ One, &
                   h * W * ( Vd1 - Vm ), &
                   h * W * Vd2, &
                   h * W * Vd3, &
                   h * W * Am - One, &
                   Zero ]
        R(:,2) = [ K / ( h * W ), &
                   Vd1, &
                   Vd2, &
                   Vd3, &
                   One - K / ( h * W ), &
                   Zero ]
        R(:,3) = [ W * Vd2, &
                   h *            Two * W**2 * Vd1 * Vd2, &
                   h * ( Gmdd22 + Two * W**2 * Vd2 * Vd2 ), &
                   h *            Two * W**2 * Vd3 * Vd2, &
                   W * Vd2 * ( Two * h * W - One ), &
                   Zero ]
        R(:,4) = [ W * Vd3, &
                   h *            Two * W**2 * Vd1 * Vd3, &
                   h *            Two * W**2 * Vd2 * Vd3, &
                   h * ( Gmdd33 + Two * W**2 * Vd3 * Vd3 ), &
                   W * Vd3 * ( Two * h * W - One ), &
                   Zero ]
        R(:,5) = [ One, &
                   h * W * ( Vd1 - Vp ), &
                   h * W * Vd2, &
                   h * W * Vd3, &
                   h * W * Ap - One, &
                   Zero ]
        R(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

        invR = inv(R)

        IF( DEBUG )THEN

          WRITE(*,'(A)') 'Debugging CharacteristicDecompositionModule (X1)...'
          WRITE(*,*)

          CALL ComputeFluxJacConsMatrix( iDim, U, G, dFdU )

          WRITE(*,*) 'EigVals:'
          DO i = 1, nCF
            WRITE(*,*) EigVals(i)
          END DO
          WRITE(*,*)

          LAMBDA = MATMUL( invR, MATMUL( dFdU, R ) )
          WRITE(*,*) 'LAMBDA = inv(R) x dFdU x R:'
          DO i = 1, nCF
            WRITE(*,*) LAMBDA(i,:)
          END DO
          WRITE(*,*)

!!$          WRITE(*,*) 'GR: R(i,:)'
!!$          DO i = 1, nCF
!!$            WRITE(*,*) R(i,:)
!!$          END DO
!!$          WRITE(*,*)
!          STOP 'End of debugging CharacteristicDecompositionModule'

        END IF

      CASE( 2 )

        EigVals = Eigenvalues_GR &
                    ( Vu1, Vu2, Vu3, Vu2, Cs(1), &
                      Gmdd11, Gmdd22, Gmdd33, Gmdd22, &
                      LapseFunction, ShiftVector(2) )

        ! --- lambda_m = EigVals(1) ---
        ! --- lambda_p = EigVals(3) ---

        ! --- Rezzolla, Eq. (7.245) ---
        LAMBDA_m = ( EigVals(1) + ShiftVector(2) ) / LapseFunction
        LAMBDA_p = ( EigVals(3) + ShiftVector(2) ) / LapseFunction

        Vm = ( Vu2 - LAMBDA_m ) / ( One / Gmdd22 - Vu2 * LAMBDA_m )
        Vp = ( Vu2 - LAMBDA_p ) / ( One / Gmdd22 - Vu2 * LAMBDA_p )

        Am = ( One / Gmdd22 - Vu2 * Vu2 ) &
             / ( One / Gmdd22 - Vu2 * LAMBDA_m )
        Ap = ( One / Gmdd22 - Vu2 * Vu2 ) &
             / ( One / Gmdd22 - Vu2 * LAMBDA_p )

        R(:,1) = [ One, &
                   h * W * Vd1, &
                   h * W * ( Vd2 - Vm ), &
                   h * W * Vd3, &
                   h * W * Am - One, &
                   Zero ]
        R(:,2) = [ K / ( h * W ), &
                   Vd1, &
                   Vd2, &
                   Vd3, &
                   One - K / ( h * W ), &
                   Zero ]
        R(:,3) = [ W * Vd1, &
                   h * ( Gmdd11 + Two * W**2 * Vd1 * Vd1 ), &
                   h *            Two * W**2 * Vd2 * Vd1, &
                   h *            Two * W**2 * Vd3 * Vd1, &
                   W * Vd1 * ( Two * h * W - One ), &
                   Zero ]
        R(:,4) = [ W * Vd3, &
                   h *            Two * W**2 * Vd1 * Vd3, &
                   h *            Two * W**2 * Vd2 * Vd3, &
                   h * ( Gmdd33 + Two * W**2 * Vd3 * Vd3 ), &
                   W * Vd3 * ( Two * h * W - One ), &
                   Zero ]
        R(:,5) = [ One, &
                   h * W * Vd1, &
                   h * W * ( Vd2 - Vp ), &
                   h * W * Vd3, &
                   h * W * Ap - One, &
                   Zero ]
        R(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

        invR = inv(R)

        IF( DEBUG )THEN

          WRITE(*,'(A)') 'Debugging CharacteristicDecompositionModule (X2)...'
          WRITE(*,*)

          CALL ComputeFluxJacConsMatrix( iDim, U, G, dFdU )

          WRITE(*,*) 'EigVals:'
          DO i = 1, nCF
            WRITE(*,*) EigVals(i)
          END DO
          WRITE(*,*)

          LAMBDA = MATMUL( invR, MATMUL( dFdU, R ) )
          WRITE(*,*) 'LAMBDA = inv(R) x dFdU x R:'
          DO i = 1, nCF
            WRITE(*,*) LAMBDA(i,:)
          END DO
          WRITE(*,*)

!!$          WRITE(*,*) 'GR: R(i,:)'
!!$          DO i = 1, nCF
!!$            WRITE(*,*) R(i,:)
!!$          END DO
!!$          WRITE(*,*)
          STOP 'End of debugging CharacteristicDecompositionModule'

        END IF

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A,I2.2)') '', 'Dimension = ', iDim
        WRITE(*,'(A5,A)') &
          '', 'Characteristic limiting not implemented for dimension > 2'
        WRITE(*,'(A5,A)') '', 'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition_GR
  

  ! --- Find the inverse of a matrix, function definition from
  !     http://fortranwiki.org/fortran/show/Matrix+inversion ---

  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  FUNCTION inv(A) RESULT(Ainv)
    REAL(DP), DIMENSION(:,:), INTENT(in)     :: A
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Ainv

    REAL(DP), DIMENSION(SIZE(A,1)) :: work  ! work array for LAPACK
    INTEGER, DIMENSION(SIZE(A,1)) :: ipiv   ! pivot indices
    INTEGER :: n, info

    ! External procedures defined in LAPACK
    EXTERNAL DGETRF
    EXTERNAL DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = SIZE(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    CALL DGETRF(n, n, Ainv, n, ipiv, info)

    IF (info /= 0) THEN
       STOP 'Matrix is numerically singular!'
    END IF

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

    IF (info /= 0) THEN
       STOP 'Matrix inversion failed!'
    END IF
  END FUNCTION inv

  SUBROUTINE ComputeFluxJacConsMatrix( iDim, U, G, dFdU )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: U(nCF), G(nGF)
    REAL(DP), INTENT(out) :: dFdU(nCF,nCF)

    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33, LapseFunction, ShiftVector(3), eta
    REAL(DP) :: D(1), V1(1), V2(1), V3(1), E(1), Ne(1), P(1), Cs(1)
    REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq, W
    REAL(DP) :: rho, epsilon, kappa, chi, h
    REAL(DP) :: dUdV(nCF,nCF), dFdV(nCF,nCF)

    Gmdd11         = G(iGF_Gm_dd_11)
    Gmdd22         = G(iGF_Gm_dd_22)
    Gmdd33         = G(iGF_Gm_dd_33)
    LapseFunction  = G(iGF_Alpha)
    ShiftVector(1) = G(iGF_Beta_1)
    ShiftVector(2) = G(iGF_Beta_2)
    ShiftVector(3) = G(iGF_Beta_3)

    CALL ComputePrimitive_GR &
           ( [U(iCF_D )], [U(iCF_S1)], [U(iCF_S2)], &
             [U(iCF_S3)], [U(iCF_E )], [U(iCF_Ne)], &
             D, V1, V2, V3, E, Ne, P, &
             [G(iGF_Gm_dd_11)], &
             [G(iGF_Gm_dd_22)], &
             [G(iGF_Gm_dd_33)] )

    CALL ComputeSoundSpeedFromPrimitive_GR( D, E, Ne, Cs )

    rho     = D(1)
    epsilon = E(1) / D(1)
    h = One + epsilon + P(1) / rho

    chi   = ( Gamma_IDEAL - One ) * epsilon
    kappa = ( Gamma_IDEAL - One ) * rho
    
    Vu1 = V1(1)
    Vu2 = V2(1)
    Vu3 = V3(1)

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3
 
    VSq = Vu1*Vd1 + Vu2*Vd2 + Vu3*Vd3

    W = One / SQRT( One - VSq )

    ! --- dU/dV ---
    dUdV(:,1) = [ W, &
                  ( One + epsilon + chi ) * W**2 * Vd1, &
                  ( One + epsilon + chi ) * W**2 * Vd2, &
                  ( One + epsilon + chi ) * W**2 * Vd3, &
                  ( One + epsilon + chi ) * W**2 - chi - W, &
                  Zero ]
    dUdV(:,2) = [ rho * W**3 * Vd1, &
                  rho * h * W**2 * ( Two * W**2 * Vd1 * Vd1 + Gmdd11 ), &
                  rho * h * W**2 *   Two * W**2 * Vd2 * Vd1, &
                  rho * h * W**2 *   Two * W**2 * Vd3 * Vd1, &
                  rho * W**3 * Vd1 * ( Two * h * W - One ), &
                  Zero ]
    dUdV(:,3) = [ rho * W**3 * Vd2, &
                  rho * h * W**2 *   Two * W**2 * Vd1 * Vd2, &
                  rho * h * W**2 * ( Two * W**2 * Vd2 * Vd2 + Gmdd22 ), &
                  rho * h * W**2 *   Two * W**2 * Vd3 * Vd2, &
                  rho * W**3 * Vd2 * ( Two * h * W - One ), &
                  Zero ]
    dUdV(:,4) = [ rho * W**3 * Vd3, &
                  rho * h * W**2 *   Two * W**2 * Vd1 * Vd3, &
                  rho * h * W**2 *   Two * W**2 * Vd2 * Vd3, &
                  rho * h * W**2 * ( Two * W**2 * Vd3 * Vd3 + Gmdd33 ), &
                  rho * W**3 * Vd3 * ( Two * h * W - One ), &
                 Zero ]
    dUdV(:,5) = [ Zero, &
                  ( rho + kappa ) * W**2 * Vd1, &
                  ( rho + kappa ) * W**2 * Vd2, &
                  ( rho + kappa ) * W**2 * Vd3, &
                  ( rho + kappa ) * W**2 - kappa, &
                  Zero ]
    dUdV(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

    ! --- Flux-Jacobian of primitive variables ---
    SELECT CASE( iDim )

      CASE( 1 )

        eta = ShiftVector(1) / LapseFunction

        dFdV(:,1) = [ Vu1 * W - eta * W, &
                      ( One + epsilon + chi ) * W**2 * Vu1 * Vd1 + chi &
                        - eta * ( One + epsilon + chi ) * W**2 * Vd1, &
                      ( One + epsilon + chi ) * W**2 * Vu1 * Vd2       &
                        - eta * ( One + epsilon + chi ) * W**2 * Vd2, &
                      ( One + epsilon + chi ) * W**2 * Vu1 * Vd3       &
                        - eta * ( One + epsilon + chi ) * W**2 * Vd3, &
                      ( One + epsilon + chi ) * W**2 * Vu1 - W * Vu1 &
                        - eta * ( W * ( W * ( One + epsilon + VSq * chi ) &
                          - One ) ), &
                      Zero ]
        dFdV(:,2) = [ rho * W * ( W**2 * Vu1 * Vd1 + One ) &
                        - eta * rho * W**3 * Vd1, &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd1 * Vd1 &
                        + Vd1 + Gmdd11 * Vu1 ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd1 * Vd1 + Gmdd11 ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd2 * Vd1 &
                        + Vd2                ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd2 * Vd1          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd3 * Vd1 &
                        + Vd3                ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd3 * Vd1          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd1 + One ) &
                        - rho * W * ( W**2 * Vu1 * Vd1 + One ) &
                        - eta * rho * W**3 * Vd1 * ( Two * h * W - One ), &
                      Zero ]

        dFdV(:,3) = [ rho * W * ( W**2 * Vu1 * Vd2       ) &
                        - eta * rho * W**3 * Vd2, &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd1 * Vd2 &
                                             ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd1 * Vd2          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd2 * Vd2 &
                        + Gmdd22 * Vu1 ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd2 * Vd2 + Gmdd22 ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd3 * Vd2 &
                                             ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd3 * Vd2          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd2       ) &
                        - rho * W * ( W**2 * Vu1 * Vd2       ) &
                        - eta * rho * W**3 * Vd2 * ( Two * h * W - One ), &
                      Zero ]

        dFdV(:,4) = [ rho * W * ( W**2 * Vu1 * Vd3       ) &
                        - eta * rho * W**3 * Vd3, &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd1 * Vd3 &
                                             ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd1 * Vd3          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd2 * Vd3 &
                                             ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd2 * Vd3          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd3 * Vd3 &
                              + Gmdd33 * Vu1 ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd3 * Vd3 + Gmdd33 ), &
                      rho * h * W**2 * ( Two * W**2 * Vu1 * Vd3       ) &
                        - rho * W * ( W**2 * Vu1 * Vd3       ) &
                        - eta * rho * W**3 * Vd3 * ( Two * h * W - One ), &
                      Zero ]
        dFdV(:,5) = [ Zero, &
                      ( rho + kappa ) * W**2 * Vu1 * Vd1 + kappa &
                        - eta * ( rho + kappa ) * W**2 * Vd1, &
                      ( rho + kappa ) * W**2 * Vu1 * Vd2         &
                        - eta * ( rho + kappa ) * W**2 * Vd2, &
                      ( rho + kappa ) * W**2 * Vu1 * Vd3         &
                        - eta * ( rho + kappa ) * W**2 * Vd3, &
                      W**2 * ( rho + kappa ) * Vu1 &
                        - eta * ( W**2 * ( rho + VSq * kappa ) ), &
                      Zero ]
        dFdV(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

      CASE( 2 )

        eta = ShiftVector(2) / LapseFunction

        dFdV(:,1) = [ Vu2 * W - eta * W, &
                      ( One + epsilon + chi ) * W**2 * Vu2 * Vd1       &
                        - eta * ( One + epsilon + chi ) * W**2 * Vd1, &
                      ( One + epsilon + chi ) * W**2 * Vu2 * Vd2 + chi &
                        - eta * ( One + epsilon + chi ) * W**2 * Vd2, &
                      ( One + epsilon + chi ) * W**2 * Vu2 * Vd3       &
                        - eta * ( One + epsilon + chi ) * W**2 * Vd3, &
                      ( One + epsilon + chi ) * W**2 * Vu2 - W * Vu2 &
                        - eta * ( W * ( W * ( One + epsilon + VSq * chi ) &
                          - One ) ), &
                      Zero ]
        dFdV(:,2) = [ rho * W * ( W**2 * Vu2 * Vd1       ) &
                        - eta * rho * W**3 * Vd1, &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd1 * Vd1 &
                        +       Gmdd11 * Vu2 ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd1 * Vd1 + Gmdd11 ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd2 * Vd1 &
                                             ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd2 * Vd1          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd3 * Vd1 &
                                             ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd3 * Vd1          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd1       ) &
                        - rho * W * ( W**2 * Vu2 * Vd1       ) &
                        - eta * rho * W**3 * Vd1 * ( Two * h * W - One ), &
                      Zero ]
        dFdV(:,3) = [ rho * W * ( W**2 * Vu2 * Vd2 + One ) &
                        - eta * rho * W**3 * Vd2, &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd1 * Vd2 &
                        + Vd1                ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd1 * Vd2          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd2 * Vd2 &
                        + Vd2 + Gmdd22 * Vu2 ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd2 * Vd2 + Gmdd22 ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd3 * Vd2 &
                        + Vd3                    ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd3 * Vd2          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd2 + One ) &
                        - rho * W * ( W**2 * Vu2 * Vd2 + One ) &
                        - eta * rho * W**3 * Vd2 * ( Two * h * W - One ), &
                      Zero ]
        dFdV(:,4) = [ rho * W * ( W**2 * Vu2 * Vd3       ) &
                        - eta * rho * W**3 * Vd3, &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd1 * Vd3 &
                                             ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd1 * Vd3          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd2 * Vd3 &
                                             ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd2 * Vd3          ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd3 * Vd3 &
                        +       Gmdd33 * Vu2 ) &
                        - eta * rho * h * W**2 &
                        * ( Two * W**2 * Vd3 * Vd3 + Gmdd33 ), &
                      rho * h * W**2 * ( Two * W**2 * Vu2 * Vd3       ) &
                        - rho * W * ( W**2 * Vu2 * Vd3       ) &
                        - eta * rho * W**3 * Vd3 * ( Two * h * W - One ), &
                      Zero ]
        dFdV(:,5) = [ Zero, &
                      ( rho + kappa ) * W**2 * Vu2 * Vd1         &
                        - eta * ( rho + kappa ) * W**2 * Vd1, &
                      ( rho + kappa ) * W**2 * Vu2 * Vd2 + kappa &
                        - eta * ( rho + kappa ) * W**2 * Vd2, &
                      ( rho + kappa ) * W**2 * Vu2 * Vd3         &
                        - eta * ( rho + kappa ) * W**2 * Vd3, &
                      W**2 * ( rho + kappa ) * Vu2 &
                        - eta * ( W**2 * ( rho + VSq * kappa ) ), &
                      Zero ]
        dFdV(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A,I2.2)') '', 'Dimension = ', iDim
        WRITE(*,'(A5,A)') &
          '', 'Characteristic limiting not implemented for dimension > 2'
        WRITE(*,'(A5,A)') '', 'Stopping...'
        STOP

    END SELECT

    dFdU = MATMUL( dFdV, inv(dUdV) )

  END SUBROUTINE ComputeFluxJacConsMatrix


END MODULE CharacteristicDecompositionModule_GR
