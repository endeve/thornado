MODULE CharacteristicDecompositionModule_GR

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Four
  USE GeometryFieldsModule, ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Beta_1,   &
    iGF_Beta_2,   &
    iGF_Beta_3,   &
    iGF_SqrtGm
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


  SUBROUTINE ComputeCharacteristicDecomposition_GR ( iDim, G, U, Rs, invRs )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF), U(nCF)
    REAL(DP), INTENT(out) :: Rs(nCF,nCF), invRs(nCF,nCF)

    REAL(DP) :: D(1), V1(1), V2(1), V3(1), E(1), Ne(1), P(1), Cs(1)
    REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33
    REAL(DP) :: W, h, LapseFunction
    REAL(DP) :: ShiftVector(3)
    REAL(DP) :: epsilon, kappa, chi, VSq, VdSq, K
    REAL(DP) :: EigVals(nCF)
    REAL(DP) :: LAMBDA_m, LAMBDA_p
    REAL(DP) :: Vm, Vp, Cm, Cp, Am, Ap, Nm, Np

    LOGICAL  :: DEBUG = .FALSE., TimeIt = .FALSE.
    REAL(DP) :: GammaXX, detGamma, DELTA, xi
    REAL(DP) :: IDENTITY(nCF,nCF)
    REAL(DP) :: dUdV(nCF,nCF)
    REAL(DP) :: R(nCF,nCF), LT(nCF,nCF), invR(nCF,nCF)
    REAL(DP) :: Timer
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

        Cm = Vd1 - Vm
        Cp = Vd1 - Vp

        Am = ( One / Gmdd11 - Vu1 * Vu1 ) &
                      / ( One / Gmdd11 - Vu1 * LAMBDA_m )
        Ap = ( One / Gmdd11 - Vu1 * Vu1 ) &
                      / ( One / Gmdd11 - Vu1 * LAMBDA_p )

        ! --- R* = Rs = MATMUL( dUdV, R ) ---

        IF( TimeIt )THEN
          WRITE(*,*)
          WRITE(*,'(A)') "Timing 'inv' function"
          WRITE(*,'(A)') "---------------------"
          Timer = 0.0_DP
          CALL Timer_Start( Timer )
        END IF

        Rs(:,1) = [ W * ( One &
                      + D * h * W**3 * ( Vd1 * Cm + Vd2**2 + Vd3**2 ) ), &
                    W**2 * ( ( One + epsilon + chi &
                      + ( D + kappa ) * ( h * W * Am - One ) ) * Vd1 &
                      + D * h**2 * W * ( Gmdd11 * Cm  + Two * W**2 * Vd1 &
                      * ( Vd1 * Cm + Vd2**2 + Vd3**2 ) ) ), &
                    W**2 * ( ( One + epsilon + chi &
                      + ( D + kappa ) * ( h * W * Am - One ) ) * Vd2 &
                      + D * h**2 * W * ( Gmdd22 * Vd2 + Two * W**2 * Vd2 &
                      * ( Vd1 * Cm + Vd2**2 + Vd3**2 ) ) ), &
                    W**2 * ( ( One + epsilon + chi &
                      + ( D + kappa ) * ( h * W * Am - One ) ) * Vd3 &
                      + D * h**2 * W * ( Gmdd33 * Vd3 + Two * W**2 * Vd3 &
                      * ( Vd1 * Cm + Vd2**2 + Vd3**2 ) ) ), &
                    W * ( W * ( One + epsilon + VSq * chi ) - One &
                      + D * h * W**3 * ( Two * H * W - One ) &
                      * ( Vd1 * Cm + Vd2**2 + Vd3**2 ) &
                      + W * ( D + VSq * kappa ) * ( h * W * Am - One ) ), &
                    Zero ]
        Rs(:,2) = [ K / h + D * W**3 * ( Vd1**2 + Vd2**2 + Vd3**2 ), &
                    W**2 * Vd1 * ( ( One + epsilon + chi ) * K / ( h * W ) &
                      + ( D + kappa ) * ( One - K / ( h * W ) ) &
                      + D * h * ( Gmdd11 + Two * W**2 * VdSq ) ), &
                    W**2 * Vd2 * ( ( One + epsilon + chi ) * K / ( h * W ) &
                      + ( D + kappa ) * ( One - K / ( h * W ) ) &
                      + D * h * ( Gmdd22 + Two * W**2 * VdSq ) ), &
                    W**2 * Vd3 * ( ( One + epsilon + chi ) * K / ( h * W ) &
                      + ( D + kappa ) * ( One - K / ( h * W ) ) &
                      + D * h * ( Gmdd33 + Two * W**2 * VdSq ) ), &
                    K / h * ( W * ( One + epsilon + VSq * chi ) - One ) &
                      + W**2 * ( ( D + VSq * kappa ) * ( One - K / ( h * W ) ) &
                      + D * W * ( Two * h * W - One ) * VdSq ), &
                    Zero ]
        Rs(:,3) = [ W**2 * Vd2 * ( One + D * h * W * ( Gmdd22 &
                      + Two * W**2 * ( Vd1**2 + Vd2**2 + Vd3**2 ) ) ), &
                    W**3 * Vd1 * Vd2 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( Two * h * W - One ) &
                      + Two * D * h**2 * W &
                      * ( Gmdd11 + Gmdd22 + Two * W**2 * VdSq ) ), &
                    W**2 * ( W * Vd2**2 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( Two * h * W - One ) ) &
                      + D * h**2 * ( ( Gmdd22 + Two * W**2 * Vd2**2 )**2 &
                      + Four * W**4 * Vd2**2 * ( Vd1**2 + Vd3**2 ) ) ), &
                    W**3 * Vd3 * Vd2 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( Two * h * W - One ) &
                      + Two * D * h**2 * W &
                      * ( Gmdd22 + Gmdd33 + Two * W**2 * VdSq ) ), &
                    W**2 * Vd2 * ( W * ( One + epsilon + VSq * chi ) - One &
                      + D * h * W * ( Two * h * W - One ) &
                      * ( Gmdd22 + Two * W**2 * ( Vd1**2 + Vd2**2 + Vd3**2 ) ) &
                      + W * ( D + VSq * kappa ) * ( Two * h * W - One ) ), &
                    Zero ]
        Rs(:,4) = [ W**2 * Vd3 * ( One &
                      + D * h * W * ( Gmdd33 + Two * W**2 * VdSq ) ), &
                    W**3 * Vd1 * Vd3 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( Two * h * W - One ) &
                      + Two * D * h**2 * W &
                      * ( Gmdd11 + Gmdd33 + Two * W**2 * VdSq ) ), &
                    W**3 * Vd2 * Vd3 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( Two * h * W - One ) &
                      + Two * D * h**2 * W &
                      * ( Gmdd22 + Gmdd33 + Two * W**2 * VdSq ) ), &
                    W**2 * ( W * Vd3**2 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( Two * h * W - One ) ) &
                      + D * h**2 * ( ( Gmdd33 + Two * W**2 * Vd3**2 )**2 &
                      + Four * W**4 * Vd3**2 * ( Vd1**2 + Vd2**2 ) ) ), &
                    W**2 * Vd3 * ( W * ( One + epsilon + VSq * chi ) - One &
                      + D * h * W * ( Two * h * W - One ) &
                      * ( Gmdd33 + Two * W**2 * VdSq ) &
                      + W * ( D + VSq * kappa ) * ( Two * h * W - One ) ), &
                    Zero ]
        Rs(:,5) = [ W * ( One &
                      + D * h * W**3 * ( Cp * Vd1 + Vd2**2 + Vd3**2 ) ), &
                    W**2 * ( Vd1 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( h * W * Ap - One ) ) &
                      + D * h**2 * W * ( Two * W**2 * Vd1 &
                      * ( Cp * Vd1 + Vd2**2 + Vd3**2 ) + Cp * Gmdd11 ) ), &
                    W**2 * ( Vd2 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( h * W * Ap - One ) ) &
                      + D * h**2 * W * ( Two * W**2 * Vd2 &
                      * ( Cp * Vd1 + Vd2**2 + Vd3**2 ) + Vd2 * Gmdd22 ) ), &
                    W**2 * ( Vd3 * ( One + epsilon + chi &
                      + ( D + kappa ) * ( h * W * Ap - One ) ) &
                      + D * h**2 * W * ( Two * W**2 * Vd3 &
                      * ( Cp * Vd1 + Vd2**2 + Vd3**2 ) + Vd3 * Gmdd33 ) ), &
                    W * ( W * ( One + epsilon + VSq * chi ) - One &
                      + D * h * W**3 * ( Two * h * W - One ) &
                      * ( Vd1 * Cp + Vd2**2 + Vd3**2 ) &
                      + W * ( D + VSq * kappa ) * ( h * W * Ap - One ) ), &
                    Zero ]
        Rs(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]


        IF( TimeIt )THEN

          CALL Timer_Stop( Timer )
          WRITE(*,'(A)') 'Time to read in elements of 6x6 matrix:'
          WRITE(*,'(ES24.16E3,A2)') Timer, ' s'

          Timer = 0.0_DP
          CALL Timer_Start( Timer )
          invRs = inv( Rs )
          CALL Timer_Stop( Timer )
          WRITE(*,'(A)') 'Time to compute inverse of 6x6 matrix:'
          WRITE(*,'(ES24.16E3,A2)') Timer, ' s'
          WRITE(*,*)

          STOP
        END IF

        IF( DEBUG )THEN

          DO i = 1, nCF
            WRITE(*,*) Rs(i,:)
          END DO
          WRITE(*,*)

         ! --- Right eigenvactor matrix ---
          R(:,1) = [ One, &
                     h * W * Cm, &
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
                     h * W * Cp, &
                     h * W * Vd2, &
                     h * W * Vd3, &
                     h * W * Ap - One, &
                     Zero ]
          R(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

          ! --- Transpose of left eigenvector matrix ---
          GammaXX  = Gmdd22 * Gmdd33
          detGamma = G(iGF_SqrtGm)**2
          xi       = GammaXX - detGamma * Vu1 * Vu1
          DELTA    = h**3 * W * ( K - One ) * ( Cp - Cm ) * xi

          Nm = ( One - K ) &
                        * ( -detGamma * Vu1 + Vp &
                        * ( W**2 * xi - GammaXX ) ) &
                        - K * W**2 * Vp * xi
          Np = ( One - K ) &
                        * ( -detGamma * Vu1 + Vm &
                        * ( W**2 * xi - GammaXX ) ) &
                        - K * W**2 * Vm * xi

          LT(1,:) = +h**2 / DELTA &
                      * [ h * W * Vp * xi + Nm, &
                          Vp * ( Two * K - One ) &
                            * ( W**2 * Vu1 * xi - GammaXX * Vu1 ) &
                            + GammaXX * ( One - K * Ap ), &
                          Vp * ( Two * K - One ) &
                            *   W**2 * Vu2 * xi, &
                          Vp * ( Two * K - One ) &
                            *   W**2 * Vu3 * xi, &
                          Nm, &
                          Zero ]
          LT(2,:) = W / ( K - One ) &
                      * [ h - W, W * Vu1, W * Vu2, W * Vu3, -W, Zero ]
          LT(3,:) = One / ( h * xi ) &
                      * [ -Gmdd33 * Vd2, &
                          Vu1 * Gmdd33 * Vd2, &
                          Gmdd33 * ( One - Vd1 * Vu1 ), &
                          Zero, &
                          -Gmdd33 * Vd2, &
                          Zero ]
          LT(4,:) = One / ( h * xi ) &
                      * [ -Gmdd22 * Vd3, &
                          Vu1 * Gmdd22 * Vd3, &
                          Zero, &
                          Gmdd22 * ( One - Vd1 * Vu1 ), &
                          -Gmdd22 * Vd3, &
                          Zero ]
          LT(5,:) = -h**2 / DELTA &
                      * [ h * W * Vm * xi + Np, &
                          Vm * ( Two * K - One ) &
                            * ( W**2 * Vu1 * xi - GammaXX * Vu1 ) &
                            + GammaXX * ( One - K * Am ), &
                          Vm * ( Two * K - One ) &
                            *   W**2 * Vu2 * xi, &
                          Vm * ( Two * K - One ) &
                            *   W**2 * Vu3 * xi, &
                          Np, &
                          Zero ]
          LT(6,:) = [ Zero, Zero, Zero, Zero, Zero, One ]

          ! --- dU/dV ---
          dUdV(:,1) = [ W, &
                        ( One + epsilon + chi ) * W**2 * Vd1, &
                        ( One + epsilon + chi ) * W**2 * Vd2, &
                        ( One + epsilon + chi ) * W**2 * Vd3, &
                        ( One + epsilon + chi ) * W**2 - chi - W, &
                        Zero ]
          dUdV(:,2) = [ D * W**3 * Vd1, &
                        D * h * W**2 * ( Two * W**2 * Vd1 * Vd1 + Gmdd11 ), &
                        D * h * W**2 *   Two * W**2 * Vd2 * Vd1, &
                        D * h * W**2 *   Two * W**2 * Vd3 * Vd1, &
                        D * W**3 * Vd1 * ( Two * h * W - One ), &
                        Zero ]
          dUdV(:,3) = [ D * W**3 * Vd2, &
                        D * h * W**2 *   Two * W**2 * Vd1 * Vd2, &
                        D * h * W**2 * ( Two * W**2 * Vd2 * Vd2 + Gmdd22 ), &
                        D * h * W**2 *   Two * W**2 * Vd3 * Vd2, &
                        D * W**3 * Vd2 * ( Two * h * W - One ), &
                        Zero ]
          dUdV(:,4) = [ D * W**3 * Vd3, &
                        D * h * W**2 *   Two * W**2 * Vd1 * Vd3, &
                        D * h * W**2 *   Two * W**2 * Vd2 * Vd3, &
                        D * h * W**2 * ( Two * W**2 * Vd3 * Vd3 + Gmdd33 ), &
                        D * W**3 * Vd3 * ( Two * h * W - One ), &
                        Zero ]
          dUdV(:,5) = [ Zero, &
                        ( D + kappa ) * W**2 * Vd1, &
                        ( D + kappa ) * W**2 * Vd2, &
                        ( D + kappa ) * W**2 * Vd3, &
                        ( D + kappa ) * W**2 - kappa, &
                        Zero ]
          dUdV(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

          Rs = MATMUL( dUdV, R )
          invRs = inv( Rs )

          DO i = 1, nCF
            WRITE(*,*) Rs(i,:)
          END DO
          WRITE(*,*)

          !STOP

        END IF

        invRs = inv( Rs )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A19,I2.2)') &
          '', 'Invalid dimension: ', iDim
        STOP

    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition_GR


  ! --- Find the inverse of a matrix; function definition from
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


  SUBROUTINE Timer_Start( Timer )

    REAL(DP) :: Timer

    Timer = MPI_WTIME( )

  END SUBROUTINE Timer_Start


  SUBROUTINE Timer_Stop( Timer )

    REAL(DP) :: Timer

    Timer = MPI_WTIME( ) - Timer

  END SUBROUTINE Timer_Stop


END MODULE CharacteristicDecompositionModule_GR
