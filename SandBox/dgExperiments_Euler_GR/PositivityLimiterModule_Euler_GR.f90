MODULE PositivityLimiterModule_Euler_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodesX_q, &
    nDOFX_X1, NodesX1, &
    nDOFX_X2, NodesX2, &
    nDOFX_X3, NodesX3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, LX_X1_Up, &
    LX_X2_Dn, LX_X2_Up, &
    LX_X3_Dn, LX_X3_Up
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, iGF_h_2, iGF_h_3, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_SqrtGm

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializePositivityLimiter_Euler_GR
  PUBLIC :: FinalizePositivityLimiter_Euler_GR
  PUBLIC :: ApplyPositivityLimiter_Euler_GR

  LOGICAL               :: UsePositivityLimiter, DEBUG = .FALSE.
  INTEGER, PARAMETER    :: nPS = 7  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total number of Positive Points
  REAL(DP)              :: Min_1, Min_2, Den
  REAL(DP), ALLOCATABLE :: U_PP(:,:), G_PP(:,:)

CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler_GR &
    ( Min_1_Option, Min_2_Option, UsePositivityLimiter_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option

    INTEGER :: i

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    UsePositivityLimiter = .TRUE.
    IF( PRESENT( UsePositivityLimiter_Option ) ) &
      UsePositivityLimiter = UsePositivityLimiter_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'InitializePositivityLimiter'
    WRITE(*,*)
    WRITE(*,'(A6,A,L1)') &
      '', 'Use Positivity Limiter: ', UsePositivityLimiter 
    WRITE(*,*)
    WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_1 = ', Min_1
    WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_2 = ', Min_2

    nPP = 0
    nPP(1) = PRODUCT( nNodesX )

    DO i = 1, 3

      IF( nNodesX(i) > 1 )THEN

        nPP(2*i:2*i+1) &
          = PRODUCT( nNodesX, MASK = [1,2,3] .NE. i )

      END IF

    END DO

    nPT = SUM( nPP )

    ALLOCATE( U_PP(nPT,nCF) )
    ALLOCATE( G_PP(nPT,nGF) )

  END SUBROUTINE InitializePositivityLimiter_Euler_GR

  
  SUBROUTINE FinalizePositivityLimiter_Euler_GR

    DEALLOCATE( U_PP )
    DEALLOCATE( G_PP )

  END SUBROUTINE FinalizePositivityLimiter_Euler_GR


  SUBROUTINE ApplyPositivityLimiter_Euler_GR &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    LOGICAL  :: NegativeStates(2), UseGeometry = .TRUE.
    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iP
    REAL(DP) :: Min_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: U_q(nDOFX,nCF), G_q(nDOFX,nGF), U_K(nCF), G_K(nGF), q(nPT)
    REAL(DP) :: SSq, alpha
    REAL(DP) :: q_K(1), alphaMax

    IF( nDOFX == 1 ) RETURN
    
    IF( .NOT. UsePositivityLimiter ) RETURN

    IF( DEBUG )THEN
      WRITE(*,*)
      WRITE(*,'(A)') 'Entering PositivityLimiterModule...'
      WRITE(*,'(A)') '-----------------------------------'
    END IF

    alphaMax = -HUGE( One )

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          IF( DEBUG ) THEN
            WRITE(*,*)
            WRITE(*,'(A,3I4)') 'iX1,iX2,iX3: ',iX1,iX2,iX3
          END IF

          U_q(1:nDOFX,1:nCF) &
            = U(1:nDOFX,iX1,iX2,iX3,1:nCF)

          G_q(1:nDOFX,1:nGF) &
            = G(1:nDOFX,iX1,iX2,iX3,1:nGF)

          NegativeStates = .FALSE.

          IF( DEBUG ) WRITE(*,'(A)') 'CALL ComputePointValues_Fluid'
          DO iCF = 1, nCF
            CALL ComputePointValues( U_q(:,iCF), U_PP(:,iCF) )
          END DO
          IF( DEBUG ) WRITE(*,'(A)') 'CALL ComputePointValues_Geometry'
          DO iGF = iGF_h_1, iGF_h_3
            CALL ComputePointValues( G_q(:,iGF), G_PP(:,iGF) )
          END DO
          CALL ComputeGeometryX_FromScaleFactors( G_PP )

          ! --- Ensure Positive Density ---

          Min_K = MINVAL( U_PP(:,iCF_D) )

          IF( Min_K < Min_1 )THEN
            IF( DEBUG ) WRITE(*,*) 'Ensuring positive density'

            ! --- Cell Average of Density---

            IF( DEBUG ) WRITE(*,'(A)') 'Computing cell-average of density'
            IF( UseGeometry )THEN
              U_K(iCF_D) &
                = SUM( WeightsX_q(:) * U_q(:,iCF_D) &
                    * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                    / SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) )
            ELSE
              U_K(iCF_D) &
                = SUM( WeightsX_q(:) * U_q(:,iCF_D) ) / SUM( WeightsX_q(:) )
            END IF

            IF( U_K(iCF_D) .LE. Min_1 )THEN
              IF( DEBUG ) WRITE(*,'(A)') 'U_K(iCF_D) < Min_1'
              U_q(:,iCF_D) = 1.01_DP * Min_1
            ELSE
              Theta_1 &
                = MIN( One, ABS( (Min_1-U_K(iCF_D) ) / (Min_K-U_K(iCF_D)) ) )

              ! --- Limit Density Towards Cell Average ---

              U_q(:,iCF_D) &
                = Theta_1 * U_q(:,iCF_D) + ( One - Theta_1 ) * U_K(iCF_D)
            END IF

            ! --- Recompute Point Value for CF_D ---
            IF( DEBUG ) WRITE(*,'(A)') 'Re-compute point value for CF_D'
            CALL ComputePointValues( U_q(:,iCF_D), U_PP(:,iCF_D) )

            NegativeStates(1) = .TRUE.

          END IF ! End of ensuring positive density IF statement

          ! --- Ensure positive q(u) a la Qin et al. (2016), JCP, 315, 323 ---

          IF( DEBUG ) WRITE(*,'(A)') 'CALL Computeq'
          CALL Computeq( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), q(1:nPT) )

          IF( ANY( q(:) < Zero ) )THEN
            IF( DEBUG ) WRITE(*,'(A)') 'Ensuring positive q'

            ! --- Cell Average ---

            IF( DEBUG ) WRITE(*,'(A)') 'Computing cell-average, fluid'
            DO iCF = 1, nCF
              IF( UseGeometry )THEN
                U_K(iCF) &
                  = SUM( WeightsX_q(:) * U_q(:,iCF) &
                      * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                      / SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) )
              ELSE
                U_K(iCF) &
                  = SUM( WeightsX_q(:) * U_q(:,iCF) ) / SUM( WeightsX_q(:) )
              END IF          
            END DO


            IF( DEBUG ) WRITE(*,'(A)') 'Computing cell-average, geometry'
            DO iGF = 1, nGF
              IF( UseGeometry )THEN
                G_K(iGF) &
                  = SUM( WeightsX_q(:) * G_q(:,iGF) &
                      * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                      / SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) )
              ELSE
                G_K(iGF) &
                  = SUM( WeightsX_q(:) * G_q(:,iGF) ) / SUM( WeightsX_q(:) )
              END IF

            END DO

            ! --- Compute q using cell-averages ---
            CALL Computeq( 1, U_K, G_K, q_K )
            IF( q_K(1) .LT. Zero )THEN

              ! --- Ensure positive tau ---
              IF( U_K(iCF_E) .LT. Zero ) U_K(iCF_E) = MAX( SqrtTiny, U_K(iCF_E) )

              IF( DEBUG ) &
                WRITE(*,'(A,ES24.16E3)') 'q_K = ', q_K(1)

              SSq =  U_K(iCF_S1)**2 / G_K(iGF_Gm_dd_11) &
                   + U_K(iCF_S2)**2 / G_K(iGF_Gm_dd_22) &
                   + U_K(iCF_S3)**2 / G_K(iGF_Gm_dd_33)
              IF( DEBUG ) &
                WRITE(*,'(A,ES11.3E3)') 'SSq = ', SSq

              IF( DEBUG ) &
                WRITE(*,'(A,ES11.3E3)') 'tau = ', U_K(iCF_E)
              IF( DEBUG ) &
                WRITE(*,'(A,ES11.3E3)') 'D   = ', U_K(iCF_D)

              ! --- Demand that q_K > 0 by modifying tau (alpha > 1) ---
              alpha = ( Min_2 - U_K(iCF_D) &
                        + SQRT( U_K(iCF_D)**2 + SSq + Min_2 ) ) / U_K(iCF_E)
              U_K(iCF_E) = alpha * U_K(iCF_E)

              alphaMax = MAX( alphaMax, alpha )

              IF( DEBUG ) &
                   WRITE(*,'(A,ES10.3E3)') 'alpha - 1 = ', alpha - 1.0_DP

              IF( DEBUG )THEN
                WRITE(*,'(A)') 'CALL Computeq after alpha limiting'
                CALL Computeq( 1, U_K, G_K, q_K )
                WRITE(*,'(A,ES24.16E3)') 'q_K = ', q_K(1)
              END IF

            END IF

            Theta_2 = One
            DO iP = 1, nPT

              IF( q(iP) < Zero ) THEN

                IF( DEBUG ) WRITE(*,'(A)') 'CALL SolveTheta_Bisection'
                CALL SolveTheta_Bisection &
                       ( U_PP(iP,1:nCF), U_K(1:nCF), &
                           G_PP(iP,1:nGF), G_K(1:nGF), Theta_P )

                Theta_2 = MIN( Theta_2, Theta_P )

              END IF

            END DO

            ! --- Limit Towards Cell Average ----

            DO iCF = 1, nCF

              U_q(:,iCF) = Theta_2 * U_q(:,iCF) + ( One - Theta_2 ) * U_K(iCF)

            END DO

            NegativeStates(2) = .TRUE.
            NegativeStates(1) = .FALSE.

          END IF ! End of ensuring positive q IF statement

          IF( NegativeStates(1) )THEN

            U(1:nDOFX,iX1,iX2,iX3,iCF_D) &
              = U_q(1:nDOFX,iCF_D)

          ELSEIF( NegativeStates(2) )THEN

            U(1:nDOFX,iX1,iX2,iX3,1:nCF) &
              = U_q(1:nDOFX,1:nCF)

          END IF

          ! --- Check q value AFTER limiting ---
          IF( DEBUG )THEN
            U_q(1:nDOFX,1:nCF) &
              = U(1:nDOFX,iX1,iX2,iX3,1:nCF)
            G_q(1:nDOFX,1:nGF) &
              = G(1:nDOFX,iX1,iX2,iX3,1:nGF)
            DO iCF = 1, nCF
              CALL ComputePointValues( U_q(:,iCF), U_PP(:,iCF) )
            END DO
            DO iGF = iGF_h_1, iGF_h_3
              CALL ComputePointValues( G_q(:,iGF), G_PP(:,iGF) )
            END DO
            CALL ComputeGeometryX_FromScaleFactors( G_PP )
            CALL Computeq( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), q(1:nPT) )
            IF( ANY( q .LT. Zero ) )THEN
              WRITE(*,'(A)') 'Positivity Limiter:'
              WRITE(*,'(A,3I5.4)') 'iX1,iX2,iX3: ', iX1, iX2, iX3
              WRITE(*,'(A,6ES12.3E3)') 'U_K(1:nCF): ', U_K(1:nCF)
              WRITE(*,'(A,5ES12.3E3)') 'q(:): ', q(:)
              STOP 'q < 0 after limiting.'
            END IF
          END IF

        END DO
      END DO
    END DO

    IF( DEBUG )THEN
      IF( alphaMax .GT. 0.0d0  ) &
        WRITE(*,'(A,ES10.3E3)') 'alphaMax - 1 = ', alphaMax - 1.0_DP
    END IF

  END SUBROUTINE ApplyPositivityLimiter_Euler_GR


  SUBROUTINE ComputePointValues( X_Q, X_P )

    REAL(DP), INTENT(in)  :: X_Q(nDOFX)
    REAL(DP), INTENT(out) :: X_P(nPT)

    INTEGER :: iOS
    
    X_P(1:nDOFX) = X_Q(1:nDOFX)

    IF( SUM( nPP(2:3) ) > 0 )THEN

      ! --- X1 ---

      iOS = nPP(1)

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X1), 1 )

      iOS = iOS + nPP(2)

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X1), 1 )

    END IF

    IF( SUM( nPP(4:5) ) > 0 )THEN

      ! --- X2 ---

      iOS = SUM( nPP(1:3) )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X2), 1 )

      iOS = iOS + nPP(4)

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X2), 1 )

    END IF

    IF( SUM( nPP(6:7) ) > 0 )THEN

      ! --- X3 ---

      iOS = SUM( nPP(1:5) )

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X3), 1 )

      iOS = iOS + nPP(6)

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Up, nDOFX_X3, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X3), 1 )

    END IF

  END SUBROUTINE ComputePointValues


  SUBROUTINE Computeq( N, U, G, q )

    INTEGER, INTENT(in)   :: N
    REAL(DP), INTENT(in)  :: U(N,nCF), G(N,nGF)
    REAL(DP), INTENT(out) :: q(N)

    q = qFun( U(:,iCF_D), U(:,iCF_S1), U(:,iCF_S2), U(:,iCF_S3), U(:,iCF_E), &
         G(:,iGF_Gm_dd_11), G(:,iGF_Gm_dd_22), G(:,iGF_Gm_dd_33) )

    RETURN
  END SUBROUTINE Computeq


  PURE REAL(DP) ELEMENTAL FUNCTION qFun( D, S1, S2, S3, tau, Gm11, Gm22, Gm33 )
    
    REAL(DP), INTENT(in) :: D, S1, S2, S3, tau, Gm11, Gm22, Gm33

    qFun = tau + D - SQRT( D**2 + S1**2 / Gm11 + S2**2 / Gm22 + S3**2 / Gm33 &
             + Min_2 )

    RETURN
  END FUNCTION qFun

  
  SUBROUTINE SolveTheta_Bisection( U_Q, U_K, G_Q, G_K, Theta_P )

    REAL(DP), INTENT(in)  :: U_Q(nCF), U_K(nCF), G_Q(nGF), G_K(nGF)
    REAL(DP), INTENT(out) :: Theta_P

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = qFun &
            ( x_a * U_Q(iCF_D)        + ( One - x_a ) * U_K(iCF_D),         &
              x_a * U_Q(iCF_S1)       + ( One - x_a ) * U_K(iCF_S1),        &
              x_a * U_Q(iCF_S2)       + ( One - x_a ) * U_K(iCF_S2),        &
              x_a * U_Q(iCF_S3)       + ( One - x_a ) * U_K(iCF_S3),        &
              x_a * U_Q(iCF_E)        + ( One - x_a ) * U_K(iCF_E),         &
              x_a * G_Q(iGF_Gm_dd_11) + ( One - x_a ) * G_K(iGF_Gm_dd_11),  &
              x_a * G_Q(iGF_Gm_dd_22) + ( One - x_a ) * G_K(iGF_Gm_dd_22),  &
              x_a * G_Q(iGF_Gm_dd_33) + ( One - x_a ) * G_K(iGF_Gm_dd_33) )

    x_b = One
    f_b = qFun &
            ( x_b * U_Q(iCF_D)        + ( One - x_b ) * U_K(iCF_D),         &
              x_b * U_Q(iCF_S1)       + ( One - x_b ) * U_K(iCF_S1),        &
              x_b * U_Q(iCF_S2)       + ( One - x_b ) * U_K(iCF_S2),        &
              x_b * U_Q(iCF_S3)       + ( One - x_b ) * U_K(iCF_S3),        &
              x_b * U_Q(iCF_E)        + ( One - x_b ) * U_K(iCF_E),         &
              x_b * G_Q(iGF_Gm_dd_11) + ( One - x_b ) * G_K(iGF_Gm_dd_11),  &
              x_b * G_Q(iGF_Gm_dd_22) + ( One - x_b ) * G_K(iGF_Gm_dd_22),  &
              x_b * G_Q(iGF_Gm_dd_33) + ( One - x_b ) * G_K(iGF_Gm_dd_33) )

    IF( .NOT. f_a * f_b < 0 )THEN

      WRITE(*,'(A6,A)') &
        '', 'SolveTheta_Bisection (Euler_GR):'
      WRITE(*,'(A8,A,I3.3)') &
        '', 'Error: No Root in Interval'
      WRITE(*,'(A8,A,2ES15.6e3)') &
        '', 'x_a, x_b = ', x_a, x_b
      WRITE(*,'(A8,A,2ES15.6e3)') &
        '', 'f_a, f_b = ', f_a, f_b
      STOP

    END IF

    dx = x_b - x_a

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx

      f_c = qFun &
            ( x_c * U_Q(iCF_D)        + ( One - x_c ) * U_K(iCF_D),         &
              x_c * U_Q(iCF_S1)       + ( One - x_c ) * U_K(iCF_S1),        &
              x_c * U_Q(iCF_S2)       + ( One - x_c ) * U_K(iCF_S2),        &
              x_c * U_Q(iCF_S3)       + ( One - x_c ) * U_K(iCF_S3),        &
              x_c * U_Q(iCF_E)        + ( One - x_c ) * U_K(iCF_E),         &
              x_c * G_Q(iGF_Gm_dd_11) + ( One - x_c ) * G_K(iGF_Gm_dd_11),  &
              x_c * G_Q(iGF_Gm_dd_22) + ( One - x_c ) * G_K(iGF_Gm_dd_22),  &
              x_c * G_Q(iGF_Gm_dd_33) + ( One - x_c ) * G_K(iGF_Gm_dd_33) )

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx < dx_min ) CONVERGED = .TRUE.

      IF( ITERATION > MAX_IT .AND. .NOT. CONVERGED )THEN

        WRITE(*,'(A6,A)') &
          '', 'SolveTheta_Bisection (Euler_GR):'
        WRITE(*,'(A8,A,I3.3)') &
          '', 'ITERATION ', ITERATION
        WRITE(*,'(A8,A,4ES15.6e3)') &
          '', 'x_a, x_c, x_b, dx = ', x_a, x_c, x_b, dx
        WRITE(*,'(A8,A,3ES15.6e3)') &
          '', 'f_a, f_c, f_b     = ', f_a, f_c, f_b

        IF( ITERATION > MAX_IT + 3 ) STOP

      END IF

    END DO

    Theta_P = x_a

  END SUBROUTINE SolveTheta_Bisection


END MODULE PositivityLimiterModule_Euler_GR
