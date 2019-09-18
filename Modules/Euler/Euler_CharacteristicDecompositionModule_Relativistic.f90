MODULE Euler_CharacteristicDecompositionModule_Relativistic

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
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    Euler_ComputePrimitive_Relativistic, &
    Euler_Eigenvalues_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive, &
    ComputePressureFromPrimitive

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: Euler_ComputeCharacteristicDecomposition_Relativistic


CONTAINS


  SUBROUTINE Euler_ComputeCharacteristicDecomposition_Relativistic &
    ( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF), U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF), invR(nCF,nCF)

    ! --- Expressions for right and left eigenvector matrices from
    !     Rezzolla & Zanotti, Relativistic Hydrodynamics, Eq. 7.240-7.248 ---

    REAL(DP) :: D(1), V1(1), V2(1), V3(1), E(1), Ne(1)
    REAL(DP) :: P(1), Cs(1)
    REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq
    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33, DetGm, LapseFunction, ShiftVector(3)
    REAL(DP) :: W, h, K
    REAL(DP) :: EigVals(nCF)
    REAL(DP) :: LAMBDA_m, LAMBDA_p
    REAL(DP) :: Vm, Vp, Am, Ap, Nm, Np, Cm, Cp
    REAL(DP) :: xi, DELTA
    REAL(DP) :: GAMMA_11, GAMMA_22

    LOGICAL  :: DEBUG = .FALSE.
    REAL(DP) :: dFdU(nCF,nCF), LAMBDA(nCF,nCF)
    INTEGER  :: i

    Gmdd11         = G(iGF_Gm_dd_11)
    Gmdd22         = G(iGF_Gm_dd_22)
    Gmdd33         = G(iGF_Gm_dd_33)
    DetGm          = G(iGF_SqrtGm)**2
    LapseFunction  = G(iGF_Alpha)
    ShiftVector(1) = G(iGF_Beta_1)
    ShiftVector(2) = G(iGF_Beta_2)
    ShiftVector(3) = G(iGF_Beta_3)

    CALL Euler_ComputePrimitive_Relativistic &
           ( [U(iCF_D )], [U(iCF_S1)], [U(iCF_S2)], &
             [U(iCF_S3)], [U(iCF_E )], [U(iCF_Ne)], &
             D, V1, V2, V3, E, Ne, &
             [G(iGF_Gm_dd_11)], &
             [G(iGF_Gm_dd_22)], &
             [G(iGF_Gm_dd_33)] )

    CALL ComputePressureFromPrimitive( D, E, Ne, P )

    Vu1 = V1(1)
    Vu2 = V2(1)
    Vu3 = V3(1)

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    VSq = Vu1*Vd1 + Vu2*Vd2 + Vu3*Vd3

    ! --- Auxiliary quantities ---
    W = One / SQRT( One - VSq )
    h = One + ( E(1) + P(1) ) / D(1)

    CALL ComputeSoundSpeedFromPrimitive( D, E, Ne, Cs )

    ! --- Rezzolla, Eq. (7.244) ---
    K = ( Gamma_IDEAL - One ) / ( Gamma_IDEAL - One - Cs(1)**2 )

    IF( DEBUG )THEN
      WRITE(*,*)
      WRITE(*,'(A)') 'Debugging CharacteristicDecompositionModule...'
      WRITE(*,*)

      ! --- Input for python program to compute eigenvectors ---
      WRITE(*,'(A)') '# Geometry'
      WRITE(*,'(A,ES24.16E3)') 'Gmdd11        = ', Gmdd11
      WRITE(*,'(A,ES24.16E3)') 'Gmdd22        = ', Gmdd22
      WRITE(*,'(A,ES24.16E3)') 'Gmdd33        = ', Gmdd33
      WRITE(*,'(A,ES24.16E3)') 'DetGm         = ', DetGm
      WRITE(*,'(A,ES24.16E3)') 'ShiftVector1  = ', ShiftVector(1)
      WRITE(*,'(A,ES24.16E3)') 'ShiftVector2  = ', ShiftVector(2)
      WRITE(*,'(A,ES24.16E3)') 'ShiftVector3  = ', ShiftVector(3)
      WRITE(*,'(A,ES24.16E3)') 'LapseFunction = ', LapseFunction
      WRITE(*,*)
      WRITE(*,'(A)') '# Fluid variables'
      WRITE(*,'(A,ES24.16E3)') 'Gamma  = ', Gamma_IDEAL
      WRITE(*,'(A,ES24.16E3)') 'rho    = ', D(1)
      WRITE(*,'(A,ES24.16E3)') 'Vu1    = ', Vu1
      WRITE(*,'(A,ES24.16E3)') 'Vu2    = ', Vu2
      WRITE(*,'(A,ES24.16E3)') 'Vu3    = ', Vu3
      WRITE(*,'(A,ES24.16E3)') 'e      = ', E(1)
      WRITE(*,*)
    END IF

    SELECT CASE( iDim )

      CASE( 1 )

        EigVals = Euler_Eigenvalues_Relativistic &
                    ( Vu1, Cs(1), Gmdd11, Vu1, Vu2, Vu3, &
                      Gmdd11, Gmdd22, Gmdd33, &
                      LapseFunction, ShiftVector(1) )

        ! --- lambda_m = EigVals(1) ---
        ! --- lambda_p = EigVals(3) ---

        ! --- Rezzolla, Eq. (7.245) ---
        LAMBDA_m = ( EigVals(1) + ShiftVector(1) ) / LapseFunction
        LAMBDA_p = ( EigVals(3) + ShiftVector(1) ) / LapseFunction

        Vm = ( Vu1 - LAMBDA_m ) / ( One / Gmdd11 - Vu1 * LAMBDA_m )
        Vp = ( Vu1 - LAMBDA_p ) / ( One / Gmdd11 - Vu1 * LAMBDA_p )

        Am = ( One / Gmdd11 - Vu1 * Vu1 ) / ( One / Gmdd11 - Vu1 * LAMBDA_m )
        Ap = ( One / Gmdd11 - Vu1 * Vu1 ) / ( One / Gmdd11 - Vu1 * LAMBDA_p )

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

        ! --- Rezzolla, Eq. (7.251) ---
        GAMMA_11 = Gmdd22 * Gmdd33
        xi = GAMMA_11 - DetGm * Vu1**2

        ! --- Rezzolla, Eq. (7.249) ---
        Nm = ( One - K ) * ( -DetGm * Vu1 + Vp * ( W**2 * xi - GAMMA_11 ) ) &
               - K * W**2 * Vp * xi
        Np = ( One - K ) * ( -DetGm * Vu1 + Vm * ( W**2 * xi - GAMMA_11 ) ) &
               - K * W**2 * Vm * xi

        ! --- Rezzolla, Eq. (7.250) ---
        Cm = Vd1 - Vm
        Cp = Vd1 - Vp
        DELTA = h**3 * W * ( K - One ) * ( Cp - Cm ) * xi

        ! --- Transpose of L ---
        invR(1,:) = + h**2 / DELTA &
                      * [ h * W * Vp * xi + Nm, &
                          Vp * Vu1 * ( Two * K - One ) &
                            * ( W**2 * xi - GAMMA_11 ) &
                            + GAMMA_11 * ( One - K * Ap ), &
                          Vp * Vu2 * ( Two * K - One ) * W**2 * xi, &
                          Vp * Vu3 * ( Two * K - One ) * W**2 * xi, &
                          Nm, &
                          Zero ]
        invR(2,:) = W / ( K - One ) &
                      * [ h - W, W * Vu1, W * Vu2, W * Vu3, -W, Zero ]
        invR(3,:) = One / ( h * xi ) &
                      * [ -Gmdd33 * Vd2, &
                          Vu1 * Gmdd33 * Vd2, &
                          Gmdd33 * ( One - Vd1 * Vu1 ), &
                          Zero, &
                          -Gmdd33 * Vd2, &
                          Zero ]
        invR(4,:) = One / ( h * xi ) &
                      * [ -Gmdd22 * Vd3, &
                          Vu1 * Gmdd22 * Vd3, &
                          Zero, &
                          Gmdd22 * ( One - Vd1 * Vu1 ), &
                          -Gmdd22 * Vd3, &
                          Zero ]
        invR(5,:) = - h**2 / DELTA &
                      * [ h * W * Vm * xi + Np, &
                          Vm * Vu1 * ( Two * K - One ) &
                            * ( W**2 * xi - GAMMA_11 ) &
                            + GAMMA_11 * ( One - K * Am ), &
                          Vm * Vu2 * ( Two * K - One ) * W**2 * xi, &
                          Vm * Vu3 * ( Two * K - One ) * W**2 * xi, &
                          Np, &
                          Zero ]
        invR(6,:) = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( DEBUG )THEN

          WRITE(*,*)
          WRITE(*,'(A)') 'X1'
          WRITE(*,*)

          CALL ComputeFluxJacConsMatrix( iDim, U, G, dFdU )

          WRITE(*,*) 'EigVals:'
          DO i = 1, nCF
            WRITE(*,'(ES11.2E3)') EigVals(i)
          END DO
          WRITE(*,*)

          LAMBDA = MATMUL( invR, MATMUL( dFdU, R ) )
!!$          WRITE(*,*) 'LAMBDA = inv(R) x dFdU x R:'
          DO i = 1, nCF
            WRITE(*,'(6ES11.2E3)') LAMBDA(i,1), LAMBDA(i,2), LAMBDA(i,3), &
                                   LAMBDA(i,4), LAMBDA(i,5), LAMBDA(i,6)
          END DO
          WRITE(*,*)

        END IF

      CASE( 2 )

        EigVals = Euler_Eigenvalues_Relativistic &
                    ( Vu2, Cs(1), Gmdd22, Vu1, Vu2, Vu3, &
                      Gmdd11, Gmdd22, Gmdd33, &
                      LapseFunction, ShiftVector(2) )

        ! --- lambda_m = EigVals(1) ---
        ! --- lambda_p = EigVals(3) ---

        ! --- Rezzolla, Eq. (7.245) ---
        LAMBDA_m = ( EigVals(1) + ShiftVector(2) ) / LapseFunction
        LAMBDA_p = ( EigVals(3) + ShiftVector(2) ) / LapseFunction

        Vm = ( Vu2 - LAMBDA_m ) / ( One / Gmdd22 - Vu2 * LAMBDA_m )
        Vp = ( Vu2 - LAMBDA_p ) / ( One / Gmdd22 - Vu2 * LAMBDA_p )

        Am = ( One / Gmdd22 - Vu2 * Vu2 ) / ( One / Gmdd22 - Vu2 * LAMBDA_m )
        Ap = ( One / Gmdd22 - Vu2 * Vu2 ) / ( One / Gmdd22 - Vu2 * LAMBDA_p )

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

        ! --- Rezzolla, Eq. (7.251) ---
        GAMMA_22 = Gmdd11 * Gmdd33
        xi = GAMMA_22 - DetGm * Vu2**2

        ! --- Rezzolla, Eq. (7.249) ---
        Nm = ( One - K ) * ( -DetGm * Vu2 + Vp * ( W**2 * xi - GAMMA_22 ) ) &
               - K * W**2 * Vp * xi
        Np = ( One - K ) * ( -DetGm * Vu2 + Vm * ( W**2 * xi - GAMMA_22 ) ) &
               - K * W**2 * Vm * xi

        ! --- Rezzolla, Eq. (7.250) ---
        Cm = Vd2 - Vm
        Cp = Vd2 - Vp
        DELTA = h**3 * W * ( K - One ) * ( Cp - Cm ) * xi

        ! --- Transpose of L, Eqs. (7.246-7.248) ---
        invR(1,:) = + h**2 / DELTA &
                      * [ h * W * Vp * xi + Nm, &
                          Vp * Vu1 * ( Two * K - One ) * W**2 * xi, &
                          Vp * Vu2 * ( Two * K - One ) &
                            * ( W**2 * xi - GAMMA_22 ) &
                            + GAMMA_22 * ( One - K * Ap ), &
                          Vp * Vu3 * ( Two * K - One ) * W**2 * xi, &
                          Nm, &
                          Zero ]
        invR(2,:) = W / ( K - One ) &
                      * [ h - W, W * Vu1, W * Vu2, W * Vu3, -W, Zero ]
        invR(3,:) = One / ( h * xi ) &
                      * [ -Gmdd33 * Vd1, &
                          Gmdd33 * ( One - Vd2 * Vu2 ), &
                          Vu2 * Gmdd33 * Vd1, &
                          Zero, &
                          -Gmdd33 * Vd1, &
                          Zero ]
        invR(4,:) = One / ( h * xi ) &
                      * [ -Gmdd11 * Vd3, &
                          Zero, &
                          Vu2 * Gmdd11 * Vd3, &
                          Gmdd11 * ( One - Vd2 * Vu2 ), &
                          -Gmdd11 * Vd3, &
                          Zero ]
        invR(5,:) = - h**2 / DELTA &
                      * [ h * W * Vm * xi + Np, &
                          Vm * Vu1 * ( Two * K - One ) * W**2 * xi, &
                          Vm * Vu2 * ( Two * K - One ) &
                            * ( W**2 * xi - GAMMA_22 ) &
                            + GAMMA_22 * ( One - K * Am ), &
                          Vm * Vu3 * ( Two * K - One ) * W**2 * xi, &
                          Np, &
                          Zero ]
        invR(6,:) = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( DEBUG )THEN

          WRITE(*,*)
          WRITE(*,'(A)') 'X2'
          WRITE(*,*)

          CALL ComputeFluxJacConsMatrix( iDim, U, G, dFdU )

          WRITE(*,*) 'EigVals:'
          DO i = 1, nCF
            WRITE(*,'(ES11.2E3)') EigVals(i)
          END DO
          WRITE(*,*)

          LAMBDA = MATMUL( invR, MATMUL( dFdU, R ) )
!!$          WRITE(*,*) 'LAMBDA = inv(R) x dFdU x R:'
          DO i = 1, nCF
            WRITE(*,'(6ES11.2E3)') LAMBDA(i,1), LAMBDA(i,2), LAMBDA(i,3), &
                                   LAMBDA(i,4), LAMBDA(i,5), LAMBDA(i,6)
          END DO
          WRITE(*,*)

        END IF

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A,I2.2)') '', 'Dimension = ', iDim
        WRITE(*,'(A5,A)') &
          '', 'Characteristic limiting not implemented for dimension > 2'
        WRITE(*,'(A5,A)') '', 'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE Euler_ComputeCharacteristicDecomposition_Relativistic


!!$  ! --- Find the inverse of a matrix, function definition from
!!$  !     http://fortranwiki.org/fortran/show/Matrix+inversion ---
!!$
!!$  ! Returns the inverse of a matrix calculated by finding the LU
!!$  ! decomposition.  Depends on LAPACK.
!!$  FUNCTION inv(A) RESULT(Ainv)
!!$    REAL(DP), DIMENSION(:,:), INTENT(in)     :: A
!!$    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Ainv
!!$
!!$    REAL(DP), DIMENSION(SIZE(A,1)) :: work  ! work array for LAPACK
!!$    INTEGER, DIMENSION(SIZE(A,1)) :: ipiv   ! pivot indices
!!$    INTEGER :: n, info
!!$
!!$    ! External procedures defined in LAPACK
!!$    EXTERNAL DGETRF
!!$    EXTERNAL DGETRI
!!$
!!$    ! Store A in Ainv to prevent it from being overwritten by LAPACK
!!$    Ainv = A
!!$    n = SIZE(A,1)
!!$
!!$    ! DGETRF computes an LU factorization of a general M-by-N matrix A
!!$    ! using partial pivoting with row interchanges.
!!$    CALL DGETRF(n, n, Ainv, n, ipiv, info)
!!$
!!$    IF (info /= 0) THEN
!!$       STOP 'Matrix is numerically singular!'
!!$    END IF
!!$
!!$    ! DGETRI computes the inverse of a matrix using the LU factorization
!!$    ! computed by DGETRF.
!!$    CALL DGETRI(n, Ainv, n, ipiv, work, n, info)
!!$
!!$    IF (info /= 0) THEN
!!$       STOP 'Matrix inversion failed!'
!!$    END IF
!!$  END FUNCTION inv

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

    CALL Euler_ComputePrimitive_Relativistic &
           ( [U(iCF_D )], [U(iCF_S1)], [U(iCF_S2)], &
             [U(iCF_S3)], [U(iCF_E )], [U(iCF_Ne)], &
             D, V1, V2, V3, E, Ne, &
             [G(iGF_Gm_dd_11)], &
             [G(iGF_Gm_dd_22)], &
             [G(iGF_Gm_dd_33)] )

    CALL ComputePressureFromPrimitive( D, E, Ne, P )

    CALL ComputeSoundSpeedFromPrimitive( D, E, Ne, Cs )

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

!!$    dFdU = MATMUL( dFdV, inv(dUdV) )

  END SUBROUTINE ComputeFluxJacConsMatrix


END MODULE Euler_CharacteristicDecompositionModule_Relativistic
