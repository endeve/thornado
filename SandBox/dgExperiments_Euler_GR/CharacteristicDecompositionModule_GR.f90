MODULE CharacteristicDecompositionModule_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, SqrtTiny
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

  PUBLIC :: ComputeCharacteristicDecomposition_GR

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_GR ( iDim, G, U, Rs, invRs )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: Rs(nCF,nCF)
    REAL(DP), INTENT(out) :: invRs(nCF,nCF)

    INTEGER :: i, j
    REAL(DP), DIMENSION(1)   :: D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, E, Ne, P
    REAL(DP), DIMENSION(1)   :: Gmdd11, Gmdd22, Gmdd33
    REAL(DP), DIMENSION(1)   :: W, h
    REAL(DP), DIMENSION(nCF) :: EigVals
    REAL(DP), DIMENSION(1)   :: Cs, epsilon, chi, kappa, VSq
    REAL(DP), DIMENSION(1)   :: DELTA_m, DELTA_p
    REAL(DP)                 :: FluxJacPrim(nCF,nCF), FluxJacCons(nCF,nCF)
    REAL(DP)                 :: DIAG(nCF,nCF)
    REAL(DP)                 :: dUdV(nCF,nCF), dVdU(nCF,nCF), R(nCF,nCF)
    LOGICAL                  :: DEBUG = .FALSE.

    CALL ComputePrimitive_GR &
           ( [U(iCF_D )], [U(iCF_S1)], [U(iCF_S2)], &
             [U(iCF_S3)], [U(iCF_E )], [U(iCF_Ne)], &
             D, Vu1, Vu2, Vu3, E, Ne, P, &
             [G(iGF_Gm_dd_11)], &
             [G(iGF_Gm_dd_22)], &
             [G(iGF_Gm_dd_33)] )
    !WRITE(*,*) D   ! D   = 1
    !WRITE(*,*) Vu1 ! Vu1 = 0
    !WRITE(*,*) Vu2 ! Vu2 = 0
    !WRITE(*,*) Vu3 ! Vu3 = 0
    !WRITE(*,*) Ne  ! Ne  = 0
    !WRITE(*,*) P   ! P   = 1

    CALL ComputeSoundSpeedFromPrimitive_GR( D, E, Ne, Cs )
    !WRITE(*,*) Cs ! Cs = 0.5164

    Gmdd11 = G(iGF_Gm_dd_11)
    Gmdd22 = G(iGF_Gm_dd_22)
    Gmdd33 = G(iGF_Gm_dd_33)
    !WRITE(*,*) Gmdd11, Gmdd22, Gmdd33 ! Gmdd11 = Gmdd22 = Gmdd33 = 1

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    VSq = Vu1*Vd1 + Vu2*Vd2 + Vu3*Vd3

    W = One / SQRT( One - VSq )
    !WRITE(*,*) W ! W = 1

    h = One + ( E + P ) / D
    !WRITE(*,*) h ! h = 5

    SELECT CASE( iDim )

      CASE( 1 )

        EigVals = Eigenvalues_GR &
                    ( Vu1(1), Vu2(1), Vu3(1), Vu1(1), Cs(1), &
                      Gmdd11(1), Gmdd22(1), Gmdd33(1), Gmdd11(1), &
                      G(iGF_Alpha), G(iGF_Beta_1) )

        ! --- lambda_m = EigVals(1) ---
        ! --- lambda_p = EigVals(3) ---
        !WRITE(*,*) EigVals(1) ! EigVals(1) = -Cs ~ -0.5164
        !WRITE(*,*) EigVals(3) ! EigVals(3) = +Cs ~ +0.5164

        ! --- Right eigenvectors from Font et al., (1994) ---
        epsilon = h - One - P / D
        !WRITE(*,*) epsilon ! epsilon = 3

        chi     = ( Gamma_IDEAL - One ) * epsilon
        !WRITE(*,*) chi ! chi = 1
        kappa   = ( Gamma_IDEAL - One ) * D
        !WRITE(*,*) kappa ! kappa = 1/3

        DELTA_m = -VSq * W**2 * EigVals(1)**2 &
                    + Two * Vu1 * W**2 * EigVals(1) - ( One + Vu1**2 * W**2 )
        DELTA_p = -VSq * W**2 * EigVals(3)**2 &
                    + Two * Vu1 * W**2 * EigVals(3) - ( One + Vu1**2 * W**2 )
        !WRITE(*,*) DELTA_m, DELTA_p ! DELTA_m = DELTA_p = -1

        R(:,1) = [ One, &
                   ( Vu1 - EigVals(1) ) * ( One - EigVals(1) * Vu1 ) &
                     / ( D * DELTA_m ), &
                   -EigVals(1) * Vu2 * ( Vu1 - EigVals(1) ) &
                     / ( D * DELTA_m ), &
                   -EigVals(1) * Vu3 * ( Vu1 - EigVals(1) ) &
                     / ( D * DELTA_m ), &
                   -chi / kappa - h * ( Vu1 - EigVals(1) )**2 * W**2 &
                     / ( kappa * DELTA_m ), &
                   Zero ]
        R(:,2) = [ -kappa, Zero, Zero, Zero, chi, Zero ]
        R(:,3) = [ -kappa, Zero, One,  Zero, chi, Zero ]
        R(:,4) = [ -kappa, Zero, Zero, One,  chi, Zero ]
        R(:,5) = [ One, &
                   ( Vu1 - EigVals(3) ) * ( One - EigVals(3) * Vu1 ) &
                     / ( D * DELTA_p ), &
                   -EigVals(3) * Vu2 * ( Vu1 - EigVals(3) ) &
                     / ( D * DELTA_p ), &
                   -EigVals(3) * Vu3 * ( Vu1 - EigVals(3) ) &
                     / ( D * DELTA_p ), &
                   -chi / kappa - h * ( Vu1 - EigVals(3) )**2 * W**2 &
                     / ( kappa * DELTA_p ), &
                   Zero ]
        R(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( DEBUG )THEN
          WRITE(*,*) 'R'
          DO i = 1, nCF
            WRITE(*,*) R(i,:)
          END DO
          WRITE(*,*)
          !STOP
        END IF

        ! --- A^{0} from Font paper ---
        dUdV(:,1) = [ W, &
                      ( One + epsilon + chi ) * Vu1 * W**2, &
                      ( One + epsilon + chi ) * Vu2 * W**2, &
                      ( One + epsilon + chi ) * Vu3 * W**2, &
                      ( One + epsilon + chi ) * W**2 - chi - W, &
                      Zero ]
        dUdV(:,2) = [ D * Vu1 * W**3, &
                      D * h * W**2 * ( One + Two * Vu1**2 * W**2 ), &
                      Two * D * h * Vu1 * Vu2 * W**4, &
                      Two * D * h * Vu1 * Vu3 * W**4, &
                      Two * D * h * Vu1 * W**4 - D * Vu1 * W**3, &
                      Zero ]
        dUdV(:,3) = [ D * Vu2 * W**3, &
                      Two * D * h * Vu1 * Vu2 * W**4, &
                      D * h * W**2 * ( One + Two * Vu2**2 * W**2 ), &
                      Two * D * h * Vu2 * Vu3 * W**4, &
                      Two * D * h * Vu2 * W**4 - D * Vu2 * W**3, &
                      Zero ]
        dUdV(:,4) = [ D * Vu3 * W**3, &
                      Two * D * h * Vu1 * Vu3 * W**4, &
                      Two * D * h * Vu2 * Vu3 * W**4, &
                      D * h * W**2 * ( One + Two * Vu3**2 * W**2 ), &
                      Two * D * h * Vu3 * W**4 - D * Vu3 * W**3, &
                      Zero ]
        dUdV(:,5) = [ Zero, &
                      ( D + kappa ) * Vu1 * W**2, &
                      ( D + kappa ) * Vu2 * W**2, &
                      ( D + kappa ) * Vu3 * W**2, &
                      ( D + kappa ) * W**2 - kappa, &
                      Zero ]
        dUdV(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( DEBUG )THEN
          WRITE(*,*) 'dUdV (A^{0})'
          DO i = 1, nCF
            WRITE(*,*) dUdV(i,:)
          END DO
          WRITE(*,*)
          !STOP
        END IF

        ! --- Convert R to R* (Eq. (30) in Font (1994)) ---
        Rs = MATMUL( dUdV, R )

        IF( DEBUG )THEN
          WRITE(*,*) 'R*'
          DO i = 1, nCF
            WRITE(*,*) Rs(i,:)
          END DO
          WRITE(*,*)
          !STOP
        END IF

        invRs = inv( Rs )

        IF( DEBUG )THEN

          ! --- Flux-Jacobian A^{1} from Font et al., (1994) ---
          FluxJacPrim(:,1) = [ Vu1 * W, &
                              ( One + epsilon + chi ) &
                                * Vu1**2 * W**2 + chi, &
                              ( One + epsilon + chi ) &
                                * Vu1 * Vu2 * W**2, &
                              ( One + epsilon + chi ) &
                                * Vu1 * Vu3 * W**2, &
                              ( One + epsilon + chi ) &
                                * Vu1 * W**2 - Vu1 * W, &
                               Zero ]
          FluxJacPrim(:,2) = [ D * W * ( One + Vu1**2 * W**2 ), &
                               Two * D * h * Vu1 * W**2 &
                                 * ( One + Vu1**2 * W**2 ), &
                               D * h * Vu2 * W**2 &
                                 * ( One + Two * Vu1**2 * W**2 ), &
                               D * h * Vu3 * W**2 &
                                 * ( One + Two * Vu1**2 * W**2 ), &
                               D * h * W**2 * ( One + Two * Vu1**2 * W**2 ) &
                                 - D * W * ( One + Vu1**2 * W**2 ), &
                               Zero ]
          FluxJacPrim(:,3) = [ D * Vu1 * Vu2 * W**3, &
                               Two * D * h * Vu1**2 * Vu2 * W**4, &
                               D * h * Vu1 * W**2 &
                                 * ( One + Two * Vu2**2 * W**2 ), &
                               Two * D * h * Vu1 * Vu2 * Vu3 * W**4, &
                               Two * D * h * Vu1 * Vu2 * W**4 &
                                 - D * Vu1 * Vu2 * W**3, &
                               Zero ]
          FluxJacPrim(:,4) = [ D * Vu1 * Vu3 * W**3, &
                               Two * D * h * Vu1**2 * Vu3 * W**4, &
                               Two * D * h * Vu1 * Vu2 * Vu3 * W**4, &
                               D * h * Vu1 * W**2 &
                                 * ( One + Two * Vu3**2 * W**2 ), &
                               Two * D * h * Vu1 * Vu3 * W**4 &
                                 - D * Vu1 * Vu3 * W**3, &
                               Zero ]
          FluxJacPrim(:,5) = [ Zero, &
                               ( D + kappa ) * Vu1**2 * W**2 + kappa, &
                               ( D + kappa ) * Vu1 * Vu2 * W**2, &
                               ( D + kappa ) * Vu1 * Vu3 * W**2, &
                               ( D + kappa ) * Vu1 * W**2, &
                               Zero ]
          FluxJacPrim(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

          WRITE(*,*) 'FluxJacPrim'
          DO i = 1, nCF
            WRITE(*,*) FluxJacPrim(i,:)
          END DO
          WRITE(*,*)
          !STOP

          dVdU = inv( dUdV )
          WRITE(*,*) 'dVdU ( (A^{0})^{-1} )'
          DO i = 1, nCF
            WRITE(*,*) dVdU(i,:)
          END DO
          WRITE(*,*)
          !STOP
          
          ! --- Convert 'A' matrices to 'B' matrices (Eq. (19)) ---
          FluxJacCons = MATMUL( FluxJacPrim, dVdU )
          WRITE(*,*) 'FluxJacCons'
          DO i = 1, nCF
            WRITE(*,*) FluxJacCons(i,:)
          END DO
          !STOP

          ! --- Diagonalize the flux-Jacobian ---
          DIAG = MATMUL( inv( Rs ), MATMUL( FluxJacCons, Rs ) )
          WRITE(*,*) 'DIAG'
          DO i = 1, nCF
            WRITE(*,*) DIAG(i,:)
          END DO
          WRITE(*,*)
          !STOP

        END IF

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A19,I2.2)') &
          '', 'Invalid dimension: ', iDim
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

END MODULE CharacteristicDecompositionModule_GR
