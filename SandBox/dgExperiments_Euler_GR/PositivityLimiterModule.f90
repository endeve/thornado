MODULE PositivityLimiterModule

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

  PUBLIC :: InitializePositivityLimiter
  PUBLIC :: FinalizePositivityLimiter
  PUBLIC :: ApplyPositivityLimiter

  LOGICAL               :: UsePositivityLimiter, DEBUG = .FALSE.
  INTEGER, PARAMETER    :: nPS = 7  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total number of Positive Points
  REAL(DP)              :: Min_1, Min_2, Den
  REAL(DP), ALLOCATABLE :: U_PP(:,:), G_PP(:,:)

CONTAINS


  SUBROUTINE InitializePositivityLimiter &
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

  END SUBROUTINE InitializePositivityLimiter

  
  SUBROUTINE FinalizePositivityLimiter

    DEALLOCATE( U_PP )
    DEALLOCATE( G_PP )

  END SUBROUTINE FinalizePositivityLimiter


  SUBROUTINE ApplyPositivityLimiter &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    LOGICAL  :: NegativeStates(2)
    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iP
    REAL(DP) :: Min_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: U_q(nDOFX,nCF), G_q(nDOFX,nGF), U_K(nCF), G_K(nGF), q(nPT), q1(nDOFX)

    IF( nDOFX == 1 ) RETURN
    
    IF( .NOT. UsePositivityLimiter ) RETURN

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          IF( DEBUG ) THEN
            WRITE(*,*)
            WRITE(*,'(A,3I4)') 'iX1,iX2,iX3',iX1,iX2,iX3
          END IF
          !WRITE(*,'(A,3I4)') 'iX1,iX2,iX3',iX1,iX2,iX3

          U_q(1:nDOFX,1:nCF) &
            = U(1:nDOFX,iX1,iX2,iX3,1:nCF)

          G_q(1:nDOFX,1:nGF) &
            = G(1:nDOFX,iX1,iX2,iX3,1:nGF)

          NegativeStates = .FALSE.

          IF( DEBUG ) WRITE(*,*) 'CALL ComputePointValues_Fluid (1)'
          CALL ComputePointValues_Fluid   ( U_q, U_PP )
          IF( DEBUG ) WRITE(*,*) 'CALL ComputePointValues_Geometry (1)'
          CALL ComputePointValues_Geometry( G_q, G_PP )

          ! --- Ensure Positive Density ---

          Min_K = MINVAL( U_PP(:,iCF_D) )

          IF( Min_K < Min_1 )THEN
             IF( DEBUG ) WRITE(*,*) 'Ensuring positive density'

            ! --- Cell Average ---

            IF( DEBUG ) WRITE(*,*) 'Compute Cell average'
            U_K(iCF_D) &
              = DOT_PRODUCT( WeightsX_q, U_q(:,iCF_D) )
!              = DOT_PRODUCT( WeightsX_q, U_q(:,iCF_D) * G_q(:,iGF_SqrtGm) ) &
!                  / DOT_PRODUCT( WeightsX_q, G_q(:,iGF_SqrtGm) )

!            Theta_1 &
!              = MIN( One, ABS( (Min_1-U_K(iCF_D) ) / (Min_K-U_K(iCF_D)) ) )
            Theta_1 = Zero

            ! --- Limit Density Towards Cell Average ---

            U_q(:,iCF_D) &
              = Theta_1 * U_q(:,iCF_D) + ( One - Theta_1 ) * U_K(iCF_D)

            ! --- Recompute Point Values ---

            IF( DEBUG ) WRITE(*,*) 'CALL ComputePointValues_Fluid (2)'
            CALL ComputePointValues_Fluid   ( U_q, U_PP )
            IF( DEBUG ) WRITE(*,*) 'CALL ComputePointValues_Geometry (2)'
            CALL ComputePointValues_Geometry( G_q, G_PP )

            NegativeStates(1) = .TRUE.

          END IF

          ! --- Ensure positive q(u) a la Wu & Tang (2015), JCP, 298, 539 ---

          IF( DEBUG ) WRITE(*,*) 'CALL Computeq'
          CALL Computeq( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), q(1:nPT) )
          !IF( DEBUG )THEN
            IF( ANY( q(:) < Min_2 ) ) WRITE(*,'(A3,4ES19.10E3)') 'q: ', q
          !END IF

          IF( ANY( q(:) < Min_2 ) )THEN
             IF( DEBUG ) WRITE(*,*) 'Ensuring positive q'

            ! --- Cell Average ---

            DO iCF = 1, nCF

              U_K(iCF) &
                = DOT_PRODUCT( WeightsX_q, U_q(:,iCF) )
!                = DOT_PRODUCT( WeightsX_q, U_q(:,iCF) * G_q(:,iGF_SqrtGm) ) &
!                    / DOT_PRODUCT( WeightsX_q, G_q(:,iGF_SqrtGm) )

            END DO

            DO iGF = 1, nGF

              G_K(iGF) &
                = DOT_PRODUCT( WeightsX_q, G_q(:,iGF) )
!                = DOT_PRODUCT( WeightsX_q, G_q(:,iGF) * G_q(:,iGF_SqrtGm) ) &
!                    / DOT_PRODUCT( WeightsX_q, G_q(:,iGF_SqrtGm) )

            END DO

!!$            Theta_2 = One
!!$            DO iP = 1, nPT
!!$
!!$              IF( q(iP) < Min_2 ) THEN
!!$
!!$                CALL SolveTheta_Bisection &
!!$                       ( U_PP(iP,1:nCF), U_K(1:nCF), &
!!$                           G_PP(iP,1:nGF), G_K(1:nGF), Min_2, Theta_P )
!!$
!!$                Theta_2 = MIN( Theta_2, Theta_P )
!!$
!!$              END IF
!!$
!!$            END DO
            Theta_2 = Zero

            ! --- Limit Towards Cell Average ----

            DO iCF = 1, nCF

              U_q(:,iCF) = Theta_2 * U_q(:,iCF) + ( One - Theta_2 ) * U_K(iCF)

            END DO

            NegativeStates(2) = .TRUE.
            NegativeStates(1) = .FALSE.

          END IF

          IF( NegativeStates(1) )THEN

            U(1:nDOFX,iX1,iX2,iX3,iCF_D) &
              = U_q(1:nDOFX,iCF_D)

          ELSEIF( NegativeStates(2) )THEN

            U(1:nDOFX,iX1,iX2,iX3,1:nCF) &
              = U_q(1:nDOFX,1:nCF)

          END IF

          IF( DEBUG )THEN
            IF( MINVAL( U(1:nDOFX,iX1,iX2,iX3,1) ) < Zero )THEN
              WRITE(*,*) 'Negative density'
            END IF
            IF( MINVAL( U(1:nDOFX,iX1,iX2,iX3,5) ) < Zero )THEN
              WRITE(*,*) 'Negative energy'
            END IF
          END IF

          ! --- Check q value AFTER limiting ---
          U_q(1:nDOFX,1:nCF) &
            = U(1:nDOFX,iX1,iX2,iX3,1:nCF)
          G_q(1:nDOFX,1:nGF) &
            = G(1:nDOFX,iX1,iX2,iX3,1:nGF)
          CALL ComputePointValues_Fluid   ( U_q, U_PP )
          CALL ComputePointValues_Geometry( G_q, G_PP )
          CALL Computeq( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), q(1:nPT) )
          IF( ANY( q .LT. Min_2 ) )THEN
            WRITE(*,'(A)') 'q<Min_2 after limiting.'
            WRITE(*,'(A3,4ES19.10E3)') 'q: ', q
            WRITE(*,'(A)') 'Stopping...'
            STOP
          END IF

        END DO
      END DO
    END DO
  
  END SUBROUTINE ApplyPositivityLimiter


  SUBROUTINE ComputePointValues_Fluid( U_Q, U_P )

    REAL(DP), INTENT(in)  :: U_Q(nDOFX,nCF)
    REAL(DP), INTENT(out) :: U_P(nPT, nCF)

    INTEGER :: iOS, iCF
    
    DO iCF = 1, nCF

      U_P(1:nDOFX,iCF) = U_Q(1:nDOFX,iCF)

      IF( SUM( nPP(2:3) ) > 0 )THEN

        ! --- X1 ---

        iOS = nPP(1)

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                 U_Q(1:nDOFX,iCF), 1, Zero, U_P(iOS+1:iOS+nDOFX_X1,iCF), 1 )

        iOS = iOS + nPP(2)

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
                 U_Q(1:nDOFX,iCF), 1, Zero, U_P(iOS+1:iOS+nDOFX_X1,iCF), 1 )

      END IF

      IF( SUM( nPP(4:5) ) > 0 )THEN

        ! --- X2 ---

        iOS = SUM( nPP(1:3) )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
                 U_Q(1:nDOFX,iCF), 1, Zero, U_P(iOS+1:iOS+nDOFX_X2,iCF), 1 )

        iOS = iOS + nPP(4)

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
                 U_Q(1:nDOFX,iCF), 1, Zero, U_P(iOS+1:iOS+nDOFX_X2,iCF), 1 )

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        ! --- X3 ---

        iOS = SUM( nPP(1:5) )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
                 U_Q(1:nDOFX,iCF), 1, Zero, U_P(iOS+1:iOS+nDOFX_X3,iCF), 1 )

        iOS = iOS + nPP(6)

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Up, nDOFX_X3, &
                 U_Q(1:nDOFX,iCF), 1, Zero, U_P(iOS+1:iOS+nDOFX_X3,iCF), 1 )

      END IF

    END DO

  END SUBROUTINE ComputePointValues_Fluid


  SUBROUTINE ComputePointValues_Geometry( G_Q, G_P )

    REAL(DP), INTENT(in)  :: G_Q(nDOFX,nGF)
    REAL(DP), INTENT(out) :: G_P(nPT,  nGF)

    INTEGER :: iOS, iGF

    ! --- Interpolate scale factors to cell faces ---

    DO iGF = iGF_h_1, iGF_h_3

      G_P(1:nDOFX,iGF) = G_Q(1:nDOFX,iGF)

      IF( SUM( nPP(2:3) ) > 0 )THEN

        ! --- X1 ---

        iOS = nPP(1)

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                 G_Q(1:nDOFX,iGF), 1, Zero, G_P(iOS+1:iOS+nDOFX_X1,iGF), 1 )

        iOS = iOS + nPP(2)

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
                 G_Q(1:nDOFX,iGF), 1, Zero, G_P(iOS+1:iOS+nDOFX_X1,iGF), 1 )

      END IF

      IF( SUM( nPP(4:5) ) > 0 )THEN

        ! --- X2 ---

        iOS = SUM( nPP(1:3) )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
                 G_Q(1:nDOFX,iGF), 1, Zero, G_P(iOS+1:iOS+nDOFX_X2,iGF), 1 )

        iOS = iOS + nPP(4)

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
                 G_Q(1:nDOFX,iGF), 1, Zero, G_P(iOS+1:iOS+nDOFX_X2,iGF), 1 )

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        ! --- X3 ---

        iOS = SUM( nPP(1:5) )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
                 G_Q(1:nDOFX,iGF), 1, Zero, G_P(iOS+1:iOS+nDOFX_X3,iGF), 1 )

        iOS = iOS + nPP(6)

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Up, nDOFX_X3, &
                 G_Q(1:nDOFX,iGF), 1, Zero, G_P(iOS+1:iOS+nDOFX_X3,iGF), 1 )

      END IF

    END DO

    CALL ComputeGeometryX_FromScaleFactors( G_P )

    RETURN
  END SUBROUTINE ComputePointValues_Geometry


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

    qFun = tau + D - SQRT( D**2 + S1**2 / Gm11 + S2**2 / Gm22 + S3**2 / Gm33 )

    RETURN
  END FUNCTION qFun

  
  SUBROUTINE SolveTheta_Bisection( U_Q, U_K, G_Q, G_K, Minq, Theta_P )

    REAL(DP), INTENT(in)  :: U_Q(nCF), U_K(nCF), G_Q(nGF), G_K(nGF), Minq
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
              x_a * G_Q(iGF_Gm_dd_33) + ( One - x_a ) * G_K(iGF_Gm_dd_33) ) &
          - Minq

    x_b = One
    f_b = qFun &
            ( x_b * U_Q(iCF_D)        + ( One - x_b ) * U_K(iCF_D),         &
              x_b * U_Q(iCF_S1)       + ( One - x_b ) * U_K(iCF_S1),        &
              x_b * U_Q(iCF_S2)       + ( One - x_b ) * U_K(iCF_S2),        &
              x_b * U_Q(iCF_S3)       + ( One - x_b ) * U_K(iCF_S3),        &
              x_b * U_Q(iCF_E)        + ( One - x_b ) * U_K(iCF_E),         &
              x_b * G_Q(iGF_Gm_dd_11) + ( One - x_b ) * G_K(iGF_Gm_dd_11),  &
              x_b * G_Q(iGF_Gm_dd_22) + ( One - x_b ) * G_K(iGF_Gm_dd_22),  &
              x_b * G_Q(iGF_Gm_dd_33) + ( One - x_b ) * G_K(iGF_Gm_dd_33) ) &
          - Minq

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
              x_c * G_Q(iGF_Gm_dd_33) + ( One - x_c ) * G_K(iGF_Gm_dd_33) ) &
            - Minq

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


END MODULE PositivityLimiterModule
