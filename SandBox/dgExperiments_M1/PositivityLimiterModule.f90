MODULE PositivityLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, nDOF
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_X1_Dn, L_X1_Up, &
    L_X2_Dn, L_X2_Up, &
    L_X3_Dn, L_X3_Up
  USE RadiationFieldsModule, ONLY: &
    nSpecies, nCR, &
    iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializePositivityLimiter
  PUBLIC :: FinalizePositivityLimiter
  PUBLIC :: ApplyPositivityLimiter

  LOGICAL               :: UsePositivityLimiter
  INTEGER, PARAMETER    :: nPS = 9  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Number of Positive Points
  REAL(DP)              :: Min_1, Max_1, Min_2
  REAL(DP)              :: Theta_FD
  REAL(DP), ALLOCATABLE :: U_PP(:,:)

CONTAINS


  SUBROUTINE InitializePositivityLimiter &
    ( Min_1_Option, Max_1_Option, Min_2_Option, UsePositivityLimiter_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Max_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option

    INTEGER :: i

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Max_1 = + HUGE( One )
    IF( PRESENT( Max_1_Option ) ) &
      Max_1 = Max_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    UsePositivityLimiter = .TRUE.
    IF( PRESENT( UsePositivityLimiter_Option ) ) &
      UsePositivityLimiter = UsePositivityLimiter_Option

    Theta_FD = One
    IF( Max_1 > One ) &
      Theta_FD = Zero

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'InitializePositivityLimiter'
    WRITE(*,*)
    WRITE(*,'(A6,A,L1)') &
      '', 'Use Positivity Limiter: ', UsePositivityLimiter 
    WRITE(*,*)
    WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_1 = ', Min_1
    WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_1 = ', Max_1
    WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_2 = ', Min_2
    WRITE(*,'(A6,A12,ES12.4E3)') '', 'Theta_FD = ', Theta_FD

    nPP = 0
    nPP(1) = PRODUCT( nNodesZ )

    DO i = 1, 4

      IF( nNodesZ(i) > 1 )THEN

        nPP(2*i:2*i+1) &
          = PRODUCT( nNodesZ, MASK = [1,2,3,4] .NE. i )

      END IF

    END DO

    nPT = SUM( nPP )

    ALLOCATE( U_PP(nPT,nCR) )

  END SUBROUTINE InitializePositivityLimiter


  SUBROUTINE FinalizePositivityLimiter

    DEALLOCATE( U_PP )

  END SUBROUTINE FinalizePositivityLimiter


  SUBROUTINE ApplyPositivityLimiter &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)    :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout) :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    LOGICAL  :: NegativeStates(2)
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iCR, iP
    REAL(DP) :: Min_K, Max_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: U_q(nDOF,nCR), U_K(nCR), Gamma(nPT)

    IF( .NOT. UsePositivityLimiter ) RETURN

    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              U_q(1:nDOF,1:nCR) &
                = U(1:nDOF,iZ1,iZ2,iZ3,iZ4,1:nCR,iS)

              NegativeStates = .FALSE.

              CALL ComputePointValues( U_q, U_PP )

              ! --- Ensure Bounded Density ---

              Min_K = MINVAL( U_PP(:,iCR_N) )
              Max_K = MAXVAL( U_PP(:,iCR_N) )

              IF( Min_K < Min_1 .OR. Max_K > Max_1 )THEN

                ! --- Cell Average ---

                U_K(iCR_N) = DOT_PRODUCT( Weights_q, U_q(:,iCR_N) )

                Theta_1 &
                  = MIN( One, &
                         ABS( (Min_1-U_K(iCR_N)) / (Min_K-U_K(iCR_N)) ), &
                         ABS( (Max_1-U_K(iCR_N)) / (Max_K-U_K(iCR_N)) ) )

                ! --- Limit Density Towards Cell Average ---

                U_q(:,iCR_N) &
                  = Theta_1 * U_q(:,iCR_N) + ( One - Theta_1 ) * U_K(iCR_N)

                ! --- Recompute Point Values ---

                CALL ComputePointValues( U_q, U_PP )

                NegativeStates(1) = .TRUE.

              END IF

              ! --- Ensure Positive Gamma ---

              CALL ComputeGamma( nPT, U_PP(1:nPT,1:nCR), Gamma(1:nPT) )

              IF( ANY( Gamma(:) < Min_2 ) )THEN

                ! --- Cell Average ---

                DO iCR = 1, nCR

                  U_K(iCR) = DOT_PRODUCT( Weights_q, U_q(:,iCR) )

                END DO

                Theta_2 = One
                DO iP = 1, nPT

                  IF( Gamma(iP) < Min_2 ) THEN

                    CALL SolveTheta_Bisection &
                           ( U_PP(iP,1:nCR), U_K(1:nCR), Min_2, Theta_P )

                    Theta_2 = MIN( Theta_2, Theta_P )

                  END IF

                END DO

                ! --- Limit Towards Cell Average ---

                DO iCR = 1, nCR

                  U_q(:,iCR) &
                    = Theta_2 * U_q(:,iCR) + ( One - Theta_2 ) * U_K(iCR)

                END DO

                NegativeStates(2) = .TRUE.
                NegativeStates(1) = .FALSE.

              END IF

              IF( NegativeStates(1) )THEN

                U(1:nDOF,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) &
                  = U_q(1:nDOF,iCR_N)

              ELSEIF( NegativeStates(2) )THEN

                U(1:nDOF,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) &
                  = U_q(1:nDOF,1:nCR)

              END IF

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ApplyPositivityLimiter


  SUBROUTINE ComputePointValues( U_Q, U_P )

    REAL(DP), INTENT(in)  :: U_Q(nDOF,nCR)
    REAL(DP), INTENT(out) :: U_P(nPT, nCR)

    INTEGER :: iOS, iCR

    DO iCR = 1, nCR

      U_P(1:nDOF,iCR) = U_Q(1:nDOF,iCR)

      IF( SUM( nPP(4:5) ) > 0 )THEN

        ! --- X1 ---

        iOS = SUM( nPP(1:3) )

        CALL DGEMV &
               ( 'N', nDOF_X1, nDOF, One, L_X1_Dn, nDOF_X1, &
                 U_Q(1:nDOF,iCR), 1, Zero, U_P(iOS+1:iOS+nDOF_X1,iCR), 1 )

        iOS = iOS + nPP(4)

        CALL DGEMV &
               ( 'N', nDOF_X1, nDOF, One, L_X1_Up, nDOF_X1, &
                 U_Q(1:nDOF,iCR), 1, Zero, U_P(iOS+1:iOS+nDOF_X1,iCR), 1 )

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        ! --- X2 ---

        iOS = SUM( nPP(1:5) )

        CALL DGEMV &
               ( 'N', nDOF_X2, nDOF, One, L_X2_Dn, nDOF_X2, &
                 U_Q(1:nDOF,iCR), 1, Zero, U_P(iOS+1:iOS+nDOF_X2,iCR), 1 )

        iOS = iOS + nPP(6)

        CALL DGEMV &
               ( 'N', nDOF_X2, nDOF, One, L_X2_Up, nDOF_X2, &
                 U_Q(1:nDOF,iCR), 1, Zero, U_P(iOS+1:iOS+nDOF_X2,iCR), 1 )

      END IF

      IF( SUM( nPP(8:9) ) > 0 )THEN

        ! --- X3 ---

        iOS = SUM( nPP(1:7) )

        CALL DGEMV &
               ( 'N', nDOF_X3, nDOF, One, L_X3_Dn, nDOF_X3, &
                 U_Q(1:nDOF,iCR), 1, Zero, U_P(iOS+1:iOS+nDOF_X3,iCR), 1 )

        iOS = iOS + nPP(8)

        CALL DGEMV &
               ( 'N', nDOF_X3, nDOF, One, L_X3_Up, nDOF_X3, &
                 U_Q(1:nDOF,iCR), 1, Zero, U_P(iOS+1:iOS+nDOF_X3,iCR), 1 )

      END IF

    END DO

  END SUBROUTINE ComputePointValues


  SUBROUTINE ComputeGamma( N, U, Gamma )

    INTEGER,  INTENT(in)  :: N
    REAL(DP), INTENT(in)  :: U(N,nCR)
    REAL(DP), INTENT(out) :: Gamma(N)

    Gamma = GammaFun( U(:,iCR_N), U(:,iCR_G1), U(:,iCR_G2), U(:,iCR_G3) )

  END SUBROUTINE ComputeGamma


  PURE REAL(DP) ELEMENTAL FUNCTION GammaFun( N, G1, G2, G3 )

    REAL(DP), INTENT(in) :: N, G1, G2, G3

    GammaFun = ( One - Theta_FD * N ) * N - SQRT( G1**2 + G2**2 + G3**2 )

    RETURN
  END FUNCTION GammaFun


  SUBROUTINE SolveTheta_Bisection( U_Q, U_K, MinGamma, Theta_P )

    REAL(DP), INTENT(in)  :: U_Q(nCR), U_K(nCR), MinGamma
    REAL(DP), INTENT(out) :: Theta_P

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = GammaFun &
            ( x_a * U_Q(iCR_N)  + ( One - x_a ) * U_K(iCR_N),   &
              x_a * U_Q(iCR_G1) + ( One - x_a ) * U_K(iCR_G1),  &
              x_a * U_Q(iCR_G2) + ( One - x_a ) * U_K(iCR_G2),  &
              x_a * U_Q(iCR_G3) + ( One - x_a ) * U_K(iCR_G3) ) &
          - MinGamma

    x_b = One
    f_b = GammaFun &
            ( x_b * U_Q(iCR_N)  + ( One - x_b ) * U_K(iCR_N),   &
              x_b * U_Q(iCR_G1) + ( One - x_b ) * U_K(iCR_G1),  &
              x_b * U_Q(iCR_G2) + ( One - x_b ) * U_K(iCR_G2),  &
              x_b * U_Q(iCR_G3) + ( One - x_b ) * U_K(iCR_G3) ) &
          - MinGamma

    IF( .NOT. f_a * f_b < 0 )THEN

      WRITE(*,'(A6,A)') &
        '', 'SolveTheta_Bisection (M1):'
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

      f_c = GammaFun &
              ( x_c * U_Q(iCR_N)  + ( One - x_c ) * U_K(iCR_N),   &
                x_c * U_Q(iCR_G1) + ( One - x_c ) * U_K(iCR_G1),  &
                x_c * U_Q(iCR_G2) + ( One - x_c ) * U_K(iCR_G2),  &
                x_c * U_Q(iCR_G3) + ( One - x_c ) * U_K(iCR_G3) ) &
            - MinGamma

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
          '', 'SolveTheta_Bisection (M1):'
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
