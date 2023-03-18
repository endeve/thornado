MODULE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic, &
    Eigenvalues_Euler_Relativistic
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE EquationOfStateModule, ONLY: &
    ComputeSoundSpeedFromPrimitive, &
    ComputePressureFromPrimitive

  USE MPI

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL


CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
    ( iDim, G, U, R, invR )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)    :: iDim
    REAL(DP), INTENT(in)    :: G(8)
    REAL(DP), INTENT(inout) :: U(nCF)
    REAL(DP), INTENT(out)   :: R(nCF,nCF), invR(nCF,nCF)

    ! --- Expressions for right and left eigenvector matrices from
    !     Rezzolla & Zanotti, Relativistic Hydrodynamics, Eq. 7.240-7.248 ---

    REAL(DP) :: D, V1, V2, V3, E, Ne
    REAL(DP) :: P, Cs
    REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq
    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33, DetGm, LapseFunction, ShiftVector(3)
    REAL(DP) :: W, h, K
    REAL(DP) :: EigVals(nCF)
    REAL(DP) :: LAMBDA_m, LAMBDA_p
    REAL(DP) :: Vm, Vp, Am, Ap, Nm, Np, Cm, Cp
    REAL(DP) :: xi, DELTA
    REAL(DP) :: GAMMA_11, GAMMA_22, GAMMA_33

    LOGICAL  :: DEBUG    = .FALSE.
    LOGICAL  :: DEBUG_X1 = .FALSE.
    LOGICAL  :: DEBUG_X2 = .FALSE.
    LOGICAL  :: DEBUG_X3 = .FALSE.
    REAL(DP) :: dFdU(nCF,nCF), LAMBDA(nCF,nCF)
    INTEGER  :: iCF, iErr

    Gmdd11         = G(1)
    Gmdd22         = G(2)
    Gmdd33         = G(3)
    DetGm          = G(4)
    LapseFunction  = G(5)
    ShiftVector(1) = G(6)
    ShiftVector(2) = G(7)
    ShiftVector(3) = G(8)

    CALL ComputePrimitive_Euler_Relativistic &
           ( U(iCF_D), U(iCF_S1), U(iCF_S2), U(iCF_S3), U(iCF_E), U(iCF_Ne), &
             D, V1, V2, V3, E, Ne, &
             Gmdd11, Gmdd22, Gmdd33, iErr )

    CALL ComputePressureFromPrimitive( D, E, Ne, P )

    Vu1 = V1
    Vu2 = V2
    Vu3 = V3

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    VSq = Vu1*Vd1 + Vu2*Vd2 + Vu3*Vd3

    ! --- Auxiliary quantities ---

    W = One / SQRT( One - VSq )
    h = One + ( E + P ) / D

    CALL ComputeSoundSpeedFromPrimitive( D, E, Ne, Cs )

    ! --- Rezzolla, Eq. (7.244) ---

    K = ( Gamma_IDEAL - One ) / ( Gamma_IDEAL - One - Cs**2 )

#if !defined(THORNADO_OMP_OL) && !defined(THORNADO_OACC)
    IF( DEBUG )THEN

      PRINT*
      PRINT '(A)', 'Debugging CharacteristicDecompositionModule...'
      PRINT*

      PRINT '(A)', '# Geometry Fields'

      PRINT*

      PRINT '(A,ES24.16E3)', 'Gmdd11        = ', Gmdd11
      PRINT '(A,ES24.16E3)', 'Gmdd22        = ', Gmdd22
      PRINT '(A,ES24.16E3)', 'Gmdd33        = ', Gmdd33
      PRINT '(A,ES24.16E3)', 'DetGm         = ', DetGm
      PRINT '(A,ES24.16E3)', 'ShiftVector1  = ', ShiftVector(1)
      PRINT '(A,ES24.16E3)', 'ShiftVector2  = ', ShiftVector(2)
      PRINT '(A,ES24.16E3)', 'ShiftVector3  = ', ShiftVector(3)
      PRINT '(A,ES24.16E3)', 'LapseFunction = ', LapseFunction
      PRINT*
      PRINT '(A)', '# Fluid Fields'
      PRINT '(A,ES24.16E3)', 'Gamma  = ', Gamma_IDEAL
      PRINT '(A,ES24.16E3)', 'rho    = ', D
      PRINT '(A,ES24.16E3)', 'Vu1    = ', Vu1
      PRINT '(A,ES24.16E3)', 'Vu2    = ', Vu2
      PRINT '(A,ES24.16E3)', 'Vu3    = ', Vu3
      PRINT '(A,ES24.16E3)', 'e      = ', E
      PRINT*

    END IF
#endif

    SELECT CASE( iDim )

      CASE( 1 )

        EigVals = Eigenvalues_Euler_Relativistic &
                    ( Vu1, Cs, Gmdd11, Vu1, Vu2, Vu3, &
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

        ! --- Rezzolla, Eqs. (7.240-7.241) ---

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

        ! --- Transpose of L, i.e., inverse of R
        !     Rezzolla, Eqs. (7.246-7.248) ---

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

!!$        invR = inv( R )

#if !defined(THORNADO_OMP_OL) && !defined(THORNADO_OACC)
        IF( DEBUG_X1 )THEN

           PRINT*
           PRINT '(2x,A)', 'Debugging characteristic decomposition (X1)'
           PRINT '(2x,A)', '-------------------------------------------'
           PRINT*

           CALL ComputeFluxJacConsMatrix( iDim, U, G, dFdU )

           PRINT '(4x,A)', 'Eigenvalues'
           PRINT '(4x,A)', '-----------'
           DO iCF = 1, nCF
             PRINT '(6x,ES10.2E3)', EigVals(iCF)
           END DO
           PRINT*

           PRINT '(4x,A)', 'Diagonalized flux-jacobian'
           PRINT '(4x,A)', '--------------------------'
           LAMBDA = MATMUL( invR, MATMUL( dFdU, R ) )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', LAMBDA(iCF,:)
           END DO
           PRINT*

           PRINT '(4x,A)', 'MATMUL( invR, R )'
           PRINT '(4x,A)', '-----------------'
           LAMBDA = MATMUL( invR, R )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', LAMBDA(iCF,:)
           END DO
           PRINT*

           PRINT '(4x,A)', 'invR - inv( R )'
           PRINT '(4x,A)', '---------------'
           LAMBDA = inv( R )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', invR(iCF,:) - LAMBDA(iCF,:)
           END DO
           PRINT*

        END IF
#endif

      CASE( 2 )

        EigVals = Eigenvalues_Euler_Relativistic &
                    ( Vu2, Cs, Gmdd22, Vu1, Vu2, Vu3, &
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

        ! --- Rezzolla, Eqs. (7.240-7.241) ---

        R(:,1) = [ One, &
                   h * W * Vd1, &
                   h * W * ( Vd2 - Vm ), &
                   h * W * Vd3, &
                   h * W * Am - One, &
                   Zero ]
        R(:,2) = [ W * Vd1, &
                   h * ( Gmdd11 + Two * W**2 * Vd1 * Vd1 ), &
                   h *            Two * W**2 * Vd2 * Vd1, &
                   h *            Two * W**2 * Vd3 * Vd1, &
                   W * Vd1 * ( Two * h * W - One ), &
                   Zero ]
        R(:,3) = [ K / ( h * W ), &
                   Vd1, &
                   Vd2, &
                   Vd3, &
                   One - K / ( h * W ), &
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

        ! --- Transpose of L, i.e., inverse of R
        !     Rezzolla, Eqs. (7.246-7.248) ---

        invR(1,:) = + h**2 / DELTA &
                      * [ h * W * Vp * xi + Nm, &
                          Vp * Vu1 * ( Two * K - One ) * W**2 * xi, &
                          Vp * Vu2 * ( Two * K - One ) &
                            * ( W**2 * xi - GAMMA_22 ) &
                            + GAMMA_22 * ( One - K * Ap ), &
                          Vp * Vu3 * ( Two * K - One ) * W**2 * xi, &
                          Nm, &
                          Zero ]
        invR(2,:) = One / ( h * xi ) &
                      * [ -Gmdd33 * Vd1, &
                          Gmdd33 * ( One - Vd2 * Vu2 ), &
                          Vu2 * Gmdd33 * Vd1, &
                          Zero, &
                          -Gmdd33 * Vd1, &
                          Zero ]
        invR(3,:) = W / ( K - One ) &
                      * [ h - W, W * Vu1, W * Vu2, W * Vu3, -W, Zero ]
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

!!$        invR = inv( R )

#if !defined(THORNADO_OMP_OL) && !defined(THORNADO_OACC)
        IF( DEBUG_X2 )THEN

           PRINT*
           PRINT '(2x,A)', 'Debugging characteristic decomposition (X2)'
           PRINT '(2x,A)', '-------------------------------------------'
           PRINT*

           CALL ComputeFluxJacConsMatrix( iDim, U, G, dFdU )

           PRINT '(4x,A)', 'Eigenvalues'
           PRINT '(4x,A)', '-----------'
           DO iCF = 1, nCF
             PRINT '(6x,ES10.2E3)', EigVals(iCF)
           END DO
           PRINT*

           PRINT '(4x,A)', 'Diagonalized flux-jacobian'
           PRINT '(4x,A)', '--------------------------'
           LAMBDA = MATMUL( invR, MATMUL( dFdU, R ) )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', LAMBDA(iCF,:)
           END DO
           PRINT*

           PRINT '(4x,A)', 'MATMUL( invR, R )'
           PRINT '(4x,A)', '-----------------'
           LAMBDA = MATMUL( invR, R )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', LAMBDA(iCF,:)
           END DO
           PRINT*

           PRINT '(4x,A)', 'invR - inv( R )'
           PRINT '(4x,A)', '---------------'
           LAMBDA = inv( R )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', invR(iCF,:) - LAMBDA(iCF,:)
           END DO
           PRINT*


        END IF
#endif

      CASE( 3 )

        EigVals = Eigenvalues_Euler_Relativistic &
                    ( Vu3, Cs, Gmdd33, Vu1, Vu2, Vu3, &
                      Gmdd11, Gmdd22, Gmdd33, &
                      LapseFunction, ShiftVector(3) )

        ! --- lambda_m = EigVals(1) ---
        ! --- lambda_p = EigVals(3) ---

        ! --- Rezzolla, Eq. (7.245) ---

        LAMBDA_m = ( EigVals(1) + ShiftVector(2) ) / LapseFunction
        LAMBDA_p = ( EigVals(3) + ShiftVector(2) ) / LapseFunction

        Vm = ( Vu3 - LAMBDA_m ) / ( One / Gmdd33 - Vu3 * LAMBDA_m )
        Vp = ( Vu3 - LAMBDA_p ) / ( One / Gmdd33 - Vu3 * LAMBDA_p )

        Am = ( One / Gmdd33 - Vu3 * Vu3 ) / ( One / Gmdd33 - Vu3 * LAMBDA_m )
        Ap = ( One / Gmdd33 - Vu3 * Vu3 ) / ( One / Gmdd33 - Vu3 * LAMBDA_p )

        ! --- Rezzolla, Eqs. (7.240-7.241) ---

        R(:,1) = [ One, &
                   h * W * Vd1, &
                   h * W * Vd2, &
                   h * W * ( Vd3 - Vm ), &
                   h * W * Am - One, &
                   Zero ]
        R(:,2) = [ W * Vd1, &
                   h * ( Gmdd11 + Two * W**2 * Vd1 * Vd1 ), &
                   h *            Two * W**2 * Vd2 * Vd1, &
                   h *            Two * W**2 * Vd3 * Vd1, &
                   W * Vd1 * ( Two * h * W - One ), &
                   Zero ]
        R(:,3) = [ W * Vd2, &
                   h *            Two * W**2 * Vd1 * Vd2, &
                   h * ( Gmdd22 + Two * W**2 * Vd2 * Vd2 ), &
                   h *            Two * W**2 * Vd3 * Vd2, &
                   W * Vd2 * ( Two * h * W - One ), &
                   Zero ]
        R(:,4) = [ K / ( h * W ), &
                   Vd1, &
                   Vd2, &
                   Vd3, &
                   One - K / ( h * W ), &
                   Zero ]
        R(:,5) = [ One, &
                   h * W * Vd1, &
                   h * W * Vd2, &
                   h * W * ( Vd3 - Vp ), &
                   h * W * Ap - One, &
                   Zero ]
        R(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

        ! --- Rezzolla, Eq. (7.251) ---

        GAMMA_33 = Gmdd11 * Gmdd22
        xi = GAMMA_33 - DetGm * Vu3**2

        ! --- Rezzolla, Eq. (7.249) ---

        Nm = ( One - K ) * ( -DetGm * Vu3 + Vp * ( W**2 * xi - GAMMA_33 ) ) &
               - K * W**2 * Vp * xi
        Np = ( One - K ) * ( -DetGm * Vu3 + Vm * ( W**2 * xi - GAMMA_33 ) ) &
               - K * W**2 * Vm * xi

        ! --- Rezzolla, Eq. (7.250) ---

        Cm = Vd3 - Vm
        Cp = Vd3 - Vp
        DELTA = h**3 * W * ( K - One ) * ( Cp - Cm ) * xi

        ! --- Transpose of L, i.e., inverse of R
        !     Rezzolla, Eqs. (7.246-7.248) ---

        invR(1,:) = + h**2 / DELTA &
                      * [ h * W * Vp * xi + Nm, &
                          Vp * Vu1 * ( Two * K - One ) * W**2 * xi, &
                          Vp * Vu2 * ( Two * K - One ) * W**2 * xi, &
                          Vp * Vu3 * ( Two * K - One ) &
                            * ( W**2 * xi - GAMMA_33 ) &
                            + GAMMA_33 * ( One - K * Ap ), &
                          Nm, &
                          Zero ]
        invR(2,:) = One / ( h * xi ) &
                      * [ -Gmdd22 * Vd1, &
                          Gmdd22 * ( One - Vd3 * Vu3 ), &
                          Zero, &
                          Vu3 * Gmdd22 * Vd1, &
                          -Gmdd22 * Vd1, &
                          Zero ]
        invR(3,:) = One / ( h * xi ) &
                      * [ -Gmdd11 * Vd2, &
                          Zero, &
                          Gmdd11 * ( One - Vd3 * Vu3 ), &
                          Vu3 * Gmdd11 * Vd2, &
                          -Gmdd11 * Vd2, &
                          Zero ]
        invR(4,:) = W / ( K - One ) &
                      * [ h - W, W * Vu1, W * Vu2, W * Vu3, -W, Zero ]
        invR(5,:) = - h**2 / DELTA &
                      * [ h * W * Vm * xi + Np, &
                          Vm * Vu1 * ( Two * K - One ) * W**2 * xi, &
                          Vm * Vu2 * ( Two * K - One ) * W**2 * xi, &
                          Vm * Vu3 * ( Two * K - One ) &
                            * ( W**2 * xi - GAMMA_33 ) &
                            + GAMMA_33 * ( One - K * Am ), &
                          Np, &
                          Zero ]
        invR(6,:) = [ Zero, Zero, Zero, Zero, Zero, One ]

!!$        invR = inv( R )

#if !defined(THORNADO_OMP_OL) && !defined(THORNADO_OACC)
        IF( DEBUG_X3 )THEN

           PRINT*
           PRINT '(2x,A)', 'Debugging characteristic decomposition (X3)'
           PRINT '(2x,A)', '-------------------------------------------'
           PRINT*

           CALL ComputeFluxJacConsMatrix( iDim, U, G, dFdU )

           PRINT '(4x,A)', 'Eigenvalues'
           PRINT '(4x,A)', '-----------'
           DO iCF = 1, nCF
             PRINT '(6x,ES10.2E3)', EigVals(iCF)
           END DO
           PRINT*

           PRINT '(4x,A)', 'Diagonalized flux-jacobian'
           PRINT '(4x,A)', '--------------------------'
           LAMBDA = MATMUL( invR, MATMUL( dFdU, R ) )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', LAMBDA(iCF,:)
           END DO
           PRINT*

           PRINT '(4x,A)', 'MATMUL( invR, R )'
           PRINT '(4x,A)', '-----------------'
           LAMBDA = MATMUL( invR, R )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', LAMBDA(iCF,:)
           END DO
           PRINT*

           PRINT '(4x,A)', 'invR - inv( R )'
           PRINT '(4x,A)', '---------------'
           LAMBDA = inv( R )
           DO iCF = 1, nCF
             PRINT '(6x,6ES11.2E3)', invR(iCF,:) - LAMBDA(iCF,:)
           END DO
           PRINT*

        END IF
#endif

    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL

  REAL(DP) FUNCTION absX( X )

    REAL(DP), INTENT(in) :: X

    IF( ABS( X ) .LT. 1.0d-15 )THEN
      absX = Zero
    ELSE
      absX = X
    END IF

    RETURN
  END FUNCTION absX


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

    INTEGER,  INTENT(in)    :: iDim
    REAL(DP), INTENT(in)    :: G(8)
    REAL(DP), INTENT(inout) :: U(nCF)
    REAL(DP), INTENT(out)   :: dFdU(nCF,nCF)

    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33, LapseFunction, ShiftVector(3), eta
    REAL(DP) :: D, V1, V2, V3, E, Ne, P, Cs
    REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, VSq, W
    REAL(DP) :: rho, epsilon, kappa, chi, h
    REAL(DP) :: dUdV(nCF,nCF), dFdV(nCF,nCF)

    REAL(DP) :: Vui(3), Vdj(3), Vdk(3)

    INTEGER :: iErr

    Gmdd11         = G(1)
    Gmdd22         = G(2)
    Gmdd33         = G(3)
    LapseFunction  = G(5)
    ShiftVector(1) = G(6)
    ShiftVector(2) = G(7)
    ShiftVector(3) = G(8)

    CALL ComputePrimitive_Euler_Relativistic &
           ( U(iCF_D), U(iCF_S1), U(iCF_S2), U(iCF_S3), U(iCF_E), U(iCF_Ne), &
             D, V1, V2, V3, E, Ne, &
             Gmdd11, Gmdd22, Gmdd33, iErr )

    CALL ComputePressureFromPrimitive( D, E, Ne, P )

    CALL ComputeSoundSpeedFromPrimitive( D, E, Ne, Cs )

    rho     = D
    epsilon = E / D
    h = One + epsilon + P / rho

    chi   = ( Gamma_IDEAL - One ) * epsilon
    kappa = ( Gamma_IDEAL - One ) * rho

    Vu1 = V1
    Vu2 = V2
    Vu3 = V3

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    Vui = [ Vu1, Vu2, Vu3 ]
    Vdj = [ Vd1, Vd2, Vd3 ]
    Vdk = [ Vd1, Vd2, Vd3 ]

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


        dFdV(:,1) = [ Vui(1) * W - eta * W, &
                      ( One + epsilon + chi ) * W**2 * Vui(1) * Vdj(1) + chi   &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(1),       &
                      ( One + epsilon + chi ) * W**2 * Vui(1) * Vdj(2) + Zero  &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(2),       &
                      ( One + epsilon + chi ) * W**2 * Vui(1) * Vdj(3) + Zero  &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(3),       &
                      ( One + epsilon + chi ) * W**2 * Vui(1) - W * Vui(1)     &
                        - eta * ( W * ( W * ( One + epsilon + VSq * chi )      &
                          - One ) ), &
                      Zero ]
        dFdV(:,2) = [ rho * W * ( W**2 * Vui(1) * Vdk(1) + One  )              &
                        - eta * rho * W**3 * Vdk(1),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(1) * Vdk(1) &
                        + Vdj(1) + Gmdd11 * Vui(1) )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(1) + Gmdd11 ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(2) * Vdk(1) &
                        + Vdj(2) +     Zero        )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(1) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(3) * Vdk(1) &
                        + Vdj(3) +     Zero        )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(1) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdk(1) + One  ) &
                        - rho * W * ( W**2 * Vui(1) * Vdk(1) + One  )          &
                        - eta * rho * W**3 * Vdk(1) * ( Two * h * W - One ),   &
                      Zero ]
        dFdV(:,3) = [ rho * W * ( W**2 * Vui(1) * Vdk(2) + Zero )              &
                        - eta * rho * W**3 * Vdk(2),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(1) * Vdk(2) &
                        + Zero + Zero              )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(2) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(2) * Vdk(2) &
                        + Zero + Gmdd22 * Vui(1)   )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(2) + Gmdd22 ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(3) * Vdk(2) &
                        + Zero + Zero              )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(2) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdk(2) + Zero ) &
                        - rho * W * ( W**2 * Vui(1) * Vdk(2) + Zero )          &
                        - eta * rho * W**3 * Vdk(2) * ( Two * h * W - One ),   &
                      Zero ]
        dFdV(:,4) = [ rho * W * ( W**2 * Vui(1) * Vdk(3) + Zero )              &
                        - eta * rho * W**3 * Vdk(3),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(1) * Vdk(3) &
                        + Zero + Zero              )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(3) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(2) * Vdk(3) &
                        + Zero + Zero              )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(3) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdj(3) * Vdk(3) &
                        + Zero + Gmdd33 * Vui(1)   )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(3) + Gmdd33   ),         &
                      rho * h * W**2 * ( Two * W**2 * Vui(1) * Vdk(3) + Zero ) &
                        - rho * W * ( W**2 * Vui(1) * Vdk(3) + Zero )          &
                        - eta * rho * W**3 * Vdk(3) * ( Two * h * W - One ),   &
                      Zero ]
        dFdV(:,5) = [ Zero,                                                    &
                      ( rho + kappa ) * W**2 * Vui(1) * Vdj(1) + kappa         &
                        - eta * ( rho + kappa ) * W**2 * Vdj(1),               &
                      ( rho + kappa ) * W**2 * Vui(1) * Vdj(2) + Zero          &
                        - eta * ( rho + kappa ) * W**2 * Vdj(2),               &
                      ( rho + kappa ) * W**2 * Vui(1) * Vdj(3) + Zero          &
                        - eta * ( rho + kappa ) * W**2 * Vdj(3),               &
                      W**2 * ( rho + kappa ) * Vui(1)                          &
                        - eta * ( W**2 * ( rho + VSq * kappa ) ),              &
                      Zero ]
        dFdV(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

      CASE( 2 )

        eta = ShiftVector(2) / LapseFunction

        dFdV(:,1) = [ Vui(2) * W - eta * W, &
                      ( One + epsilon + chi ) * W**2 * Vui(2) * Vdj(1) + Zero  &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(1),       &
                      ( One + epsilon + chi ) * W**2 * Vui(2) * Vdj(2) + chi   &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(2),       &
                      ( One + epsilon + chi ) * W**2 * Vui(2) * Vdj(3) + Zero  &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(3),       &
                      ( One + epsilon + chi ) * W**2 * Vui(2) - W * Vui(2)     &
                        - eta * ( W * ( W * ( One + epsilon + VSq * chi )      &
                          - One ) ), &
                      Zero ]
        dFdV(:,2) = [ rho * W * ( W**2 * Vui(2) * Vdk(1) + Zero )              &
                        - eta * rho * W**3 * Vdk(1),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(1) * Vdk(1) &
                        + Zero  + Gmdd11 * Vui(2)  )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(1) + Gmdd11 ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(2) * Vdk(1) &
                        + Zero  +     Zero         )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(1) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(3) * Vdk(1) &
                        + Zero  +     Zero         )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(1) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdk(1) + Zero ) &
                        - rho * W * ( W**2 * Vui(2) * Vdk(1) + Zero )          &
                        - eta * rho * W**3 * Vdk(1) * ( Two * h * W - One ),   &
                      Zero ]
        dFdV(:,3) = [ rho * W * ( W**2 * Vui(2) * Vdk(2) + One  )              &
                        - eta * rho * W**3 * Vdk(2),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(1) * Vdk(2) &
                        + Vdj(1) + Zero            )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(2) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(2) * Vdk(2) &
                        + Vdj(2) + Gmdd22 * Vui(2) )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(2) + Gmdd22 ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(3) * Vdk(2) &
                        + Vdj(3) + Zero            )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(2) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdk(2) + One  ) &
                        - rho * W * ( W**2 * Vui(2) * Vdk(2) + One  )          &
                        - eta * rho * W**3 * Vdk(2) * ( Two * h * W - One ),   &
                      Zero ]
        dFdV(:,4) = [ rho * W * ( W**2 * Vui(2) * Vdk(3) + Zero )              &
                        - eta * rho * W**3 * Vdk(3),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(1) * Vdk(3) &
                        + Zero + Zero              )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(3) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(2) * Vdk(3) &
                        + Zero + Zero              )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(3) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdj(3) * Vdk(3) &
                        + Zero + Gmdd33 * Vui(2)   )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(3) + Gmdd33   ),         &
                      rho * h * W**2 * ( Two * W**2 * Vui(2) * Vdk(3) + Zero ) &
                        - rho * W * ( W**2 * Vui(2) * Vdk(3) + Zero )          &
                        - eta * rho * W**3 * Vdk(3) * ( Two * h * W - One ),   &
                      Zero ]
        dFdV(:,5) = [ Zero,                                                    &
                      ( rho + kappa ) * W**2 * Vui(2) * Vdj(1) + Zero          &
                        - eta * ( rho + kappa ) * W**2 * Vdj(1),               &
                      ( rho + kappa ) * W**2 * Vui(2) * Vdj(2) + kappa         &
                        - eta * ( rho + kappa ) * W**2 * Vdj(2),               &
                      ( rho + kappa ) * W**2 * Vui(2) * Vdj(3) + Zero          &
                        - eta * ( rho + kappa ) * W**2 * Vdj(3),               &
                      W**2 * ( rho + kappa ) * Vui(2)                          &
                        - eta * ( W**2 * ( rho + VSq * kappa ) ),              &
                      Zero ]
        dFdV(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

      CASE( 3 )

        eta = ShiftVector(3) / LapseFunction

        dFdV(:,1) = [ Vui(3) * W - eta * W, &
                      ( One + epsilon + chi ) * W**2 * Vui(3) * Vdj(1) + Zero  &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(1),       &
                      ( One + epsilon + chi ) * W**2 * Vui(3) * Vdj(2) + Zero  &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(2),       &
                      ( One + epsilon + chi ) * W**2 * Vui(3) * Vdj(3) + chi   &
                        - eta * ( One + epsilon + chi ) * W**2 * Vdj(3),       &
                      ( One + epsilon + chi ) * W**2 * Vui(3) - W * Vui(3)     &
                        - eta * ( W * ( W * ( One + epsilon + VSq * chi )      &
                          - One ) ), &
                      Zero ]
        dFdV(:,2) = [ rho * W * ( W**2 * Vui(3) * Vdk(1) + Zero )              &
                        - eta * rho * W**3 * Vdk(1),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(1) * Vdk(1) &
                        + Zero  + Gmdd11 * Vui(3)  )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(1) + Gmdd11 ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(2) * Vdk(1) &
                        + Zero  +     Zero         )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(1) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(3) * Vdk(1) &
                        + Zero  +     Zero         )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(1) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdk(1) + Zero ) &
                        - rho * W * ( W**2 * Vui(3) * Vdk(1) + Zero )          &
                        - eta * rho * W**3 * Vdk(1) * ( Two * h * W - One ),   &
                      Zero ]
        dFdV(:,3) = [ rho * W * ( W**2 * Vui(3) * Vdk(2) + Zero )              &
                        - eta * rho * W**3 * Vdk(2),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(1) * Vdk(2) &
                        + Zero  + Zero             )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(2) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(2) * Vdk(2) &
                        + Zero  + Gmdd22 * Vui(3)  )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(2) + Gmdd22 ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(3) * Vdk(2) &
                        + Zero  +     Zero         )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(2) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdk(2) + Zero ) &
                        - rho * W * ( W**2 * Vui(3) * Vdk(2) + Zero )          &
                        - eta * rho * W**3 * Vdk(2) * ( Two * h * W - One ),   &
                      Zero ]
        dFdV(:,4) = [ rho * W * ( W**2 * Vui(3) * Vdk(3) + One  )              &
                        - eta * rho * W**3 * Vdk(3),                           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(1) * Vdk(3) &
                        + Vdj(1) + Zero            )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(1) * Vdk(3) + Zero   ),           &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(2) * Vdk(3) &
                        + Vdj(2) + Zero            )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(2) * Vdk(3) + Zero  ),            &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdj(3) * Vdk(3) &
                        + Vdj(3) + Gmdd33 * Vui(3) )                           &
                        - eta * rho * h * W**2                                 &
                        * ( Two * W**2 * Vdj(3) * Vdk(3) + Gmdd33   ),         &
                      rho * h * W**2 * ( Two * W**2 * Vui(3) * Vdk(3) + One  ) &
                        - rho * W * ( W**2 * Vui(3) * Vdk(3) + One  )          &
                        - eta * rho * W**3 * Vdk(3) * ( Two * h * W - One ),   &
                      Zero ]

        dFdV(:,5) = [ Zero,                                                    &
                      ( rho + kappa ) * W**2 * Vui(3) * Vdj(1) + Zero          &
                        - eta * ( rho + kappa ) * W**2 * Vdj(1),               &
                      ( rho + kappa ) * W**2 * Vui(3) * Vdj(2) + Zero          &
                        - eta * ( rho + kappa ) * W**2 * Vdj(2),               &
                      ( rho + kappa ) * W**2 * Vui(3) * Vdj(3) + kappa         &
                        - eta * ( rho + kappa ) * W**2 * Vdj(3),               &
                      W**2 * ( rho + kappa ) * Vui(3)                          &
                        - eta * ( W**2 * ( rho + VSq * kappa ) ),              &
                      Zero ]

        dFdV(:,6) = [ Zero, Zero, Zero, Zero, Zero, One ]

    END SELECT

    dFdU = MATMUL( dFdV, inv(dUdV) )

  END SUBROUTINE ComputeFluxJacConsMatrix


END MODULE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL
