MODULE Euler_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX

  USE ReferenceElementModuleX, ONLY: &
    NodesX_q, WeightsX_q, &
    nDOFX_X1, NodesX1, &
    nDOFX_X2, NodesX2, &
    nDOFX_X3, NodesX3
    

  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, LX_X1_Up, &
    LX_X2_Dn, LX_X2_Up, &
    LX_X3_Dn, LX_X3_Up

  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm

  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_InitializePositivityLimiter
  PUBLIC :: Euler_FinalizePositivityLimiter
  PUBLIC :: Euler_ApplyPositivityLimiter

  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: Verbose
  INTEGER, PARAMETER    :: nPS = 7
  INTEGER               :: nPP(nPS)
  INTEGER               :: nPT
  REAL(DP)              :: Min_1, Min_2
  REAL(DP), ALLOCATABLE :: U_PP(:,:)

CONTAINS


  SUBROUTINE Euler_InitializePositivityLimiter &
    ( Min_1_Option, Min_2_Option, UsePositivityLimiter_Option, Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

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

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'Euler_InitializePositivityLimiter'
      WRITE(*,'(A2,A)') '',    '-------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A6,A,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter
      WRITE(*,*)
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_1 = ', Min_1
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_2 = ', Min_2
    END IF

    ! --- Compute Number of Positive Points per Set ---

    nPP = 0
    nPP(1) = PRODUCT( nNodesX )

    DO i = 1, 3

      IF( nNodesX(i) > 1 )THEN

        nPP(2*i:2*i+1) &
          = PRODUCT( nNodesX, MASK = [1,2,3] .NE. i )

      END IF

    END DO

    ! --- Total Number of Positive Points ---

    nPT = SUM( nPP )

    ! --- Conserved Variables in Positive Points ---

    ALLOCATE( U_PP(nPT,nCF) )

  END SUBROUTINE Euler_InitializePositivityLimiter


  SUBROUTINE Euler_FinalizePositivityLimiter

    DEALLOCATE( U_PP )

  END SUBROUTINE Euler_FinalizePositivityLimiter


  SUBROUTINE Euler_ApplyPositivityLimiter &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    LOGICAL  :: NegativeStates(2)
    INTEGER  :: iX1, iX2, iX3, iCF, iP
    REAL(DP) :: Min_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: U_q(nDOFX,nCF), U_K(nCF), IntE(nPT)

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          U_q(1:nDOFX,1:nCF) = U(1:nDOFX,iX1,iX2,iX3,1:nCF)

          NegativeStates = .FALSE.

          CALL ComputePointValues_Fluid( U_q, U_PP )

          ! --- Ensure Positive Mass Density ---

          Min_K = MINVAL( U_PP(:,iCF_D) )

          IF( Min_K < Min_1 )THEN

            ! --- Cell Average ---

            U_K(iCF_D) = SUM( WeightsX_q(:) * U_q(:,iCF_D) * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                           / SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) )

            Theta_1 = MIN( One, (U_K(iCF_D)-Min_1)/(U_K(iCF_D)-Min_K) )

            ! --- Limit Density Towards Cell Average ---

            U_q(:,iCF_D) = Theta_1 * U_q(:,iCF_D) &
                           + ( One - Theta_1 ) * U_K(iCF_D)

            ! --- Recompute Point Values ---

            CALL ComputePointValues_Fluid( U_q, U_PP )

            ! --- Flag for Negative Density ---

            NegativeStates(1) = .TRUE.

          END IF

          ! --- Ensure Positive Internal Energy Density ---

          CALL ComputeInternalEnergyDensity &
                 ( nPT, U_PP(1:nPT,1:nCF), IntE(1:nPT) )

          IF( ANY( IntE(:) < Min_2 ) )THEN

            ! --- Cell Average ---

            DO iCF = 1, nCF

              U_K(iCF) = SUM( WeightsX_q(:) * U_q(:,iCF) * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                           / SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) )

            END DO

            Theta_2 = One
            DO iP = 1, nPT

              IF( IntE(iP) < Min_2 )THEN

                CALL SolveTheta_Bisection &
                       ( U_PP(iP,1:nCF), U_K(1:nCF), Min_2, Theta_P )

                Theta_2 = MIN( Theta_2, Theta_P )

              END IF

            END DO

            ! --- Limit Towards Cell Average ---

            DO iCF = 1, nCF

              U_q(:,iCF) = Theta_2 * U_q(:,iCF) &
                           + ( One - Theta_2 ) * U_K(iCF)

            END DO

            ! --- Flag for Negative Internal Energy Density ---

            NegativeStates(1) = .FALSE.
            NegativeStates(2) = .TRUE.

          END IF

          IF( NegativeStates(1) )THEN

            U(1:nDOFX,iX1,iX2,iX3,iCF_D) &
              = U_q(1:nDOFX,iCF_D)

          ELSEIF( NegativeStates(2) )THEN

            U(1:nDOFX,iX1,iX2,iX3,1:nCF) &
              = U_q(1:nDOFX,1:nCF)

          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE Euler_ApplyPositivityLimiter


  SUBROUTINE ComputePointValues_Fluid( U_q, U_p )

    REAL(DP), INTENT(in)  :: U_q(nDOFX,nCF)
    REAL(DP), INTENT(out) :: U_p(nPT,  nCF)

    INTEGER :: iCF, iOS

    DO iCF = 1, nCF

      U_p(1:nDOFX,iCF) = U_q(1:nDOFX,iCF)

      IF( SUM( nPP(2:3) ) > 0 )THEN

        ! --- Points of Faces with Normal in X1 Direction ---

        iOS = nPP(1)

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
                 U_q(1:nDOFX,iCF), 1, Zero, U_p(iOS+1:iOS+nDOFX_X1,iCF), 1 )

        iOS = iOS + nPP(2)

        CALL DGEMV &
               ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
                 U_q(1:nDOFX,iCF), 1, Zero, U_p(iOS+1:iOS+nDOFX_X1,iCF), 1 )

      END IF

      IF( SUM( nPP(4:5) ) > 0 )THEN

        ! --- Points of Faces with Normal in X2 Direction ---

        iOS = SUM( nPP(1:3) )

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
                 U_q(1:nDOFX,iCF), 1, Zero, U_p(iOS+1:iOS+nDOFX_X2,iCF), 1 )

        iOS = iOS + nPP(4)

        CALL DGEMV &
               ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
                 U_q(1:nDOFX,iCF), 1, Zero, U_p(iOS+1:iOS+nDOFX_X2,iCF), 1 )

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        ! --- Points of Faces with Normal in X3 Direction ---

        iOS = SUM( nPP(1:5) )

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
                 U_q(1:nDOFX,iCF), 1, Zero, U_p(iOS+1:iOS+nDOFX_X3,iCF), 1 )

        iOS = iOS + nPP(6)

        CALL DGEMV &
               ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Up, nDOFX_X3, &
                 U_q(1:nDOFX,iCF), 1, Zero, U_p(iOS+1:iOS+nDOFX_X3,iCF), 1 )

      END IF

    END DO

  END SUBROUTINE ComputePointValues_Fluid


  SUBROUTINE ComputeInternalEnergyDensity( N, U, IntE )

    INTEGER,  INTENT(in)  :: N
    REAL(DP), INTENT(in)  :: U(N,nCF)
    REAL(DP), INTENT(out) :: IntE(N)

    IntE &
      = eFun( U(:,iCF_D), U(:,iCF_S1), U(:,iCF_S2), U(:,iCF_S3), U(:,iCF_E) )

  END SUBROUTINE ComputeInternalEnergyDensity


  PURE REAL(DP) ELEMENTAL FUNCTION eFun( D, S1, S2, S3, E )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, E

    eFun = E - Half * ( S1**2 + S2**2 + S3**2 ) / D

    RETURN
  END FUNCTION eFun


  SUBROUTINE SolveTheta_Bisection( U_Q, U_K, MinE, Theta_P )

    REAL(DP), INTENT(in)  :: U_Q(nCF), U_K(nCF), MinE
    REAL(DP), INTENT(out) :: Theta_P

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = eFun &
            ( x_a * U_Q(iCF_D)  + ( One - x_a ) * U_K(iCF_D),  &
              x_a * U_Q(iCF_S1) + ( One - x_a ) * U_K(iCF_S1), &
              x_a * U_Q(iCF_S2) + ( One - x_a ) * U_K(iCF_S2), &
              x_a * U_Q(iCF_S3) + ( One - x_a ) * U_K(iCF_S3), &
              x_a * U_Q(iCF_E)  + ( One - x_a ) * U_K(iCF_E) ) &
          - MinE

    x_b = One
    f_b = eFun &
            ( x_b * U_Q(iCF_D)  + ( One - x_b ) * U_K(iCF_D),  &
              x_b * U_Q(iCF_S1) + ( One - x_b ) * U_K(iCF_S1), &
              x_b * U_Q(iCF_S2) + ( One - x_b ) * U_K(iCF_S2), &
              x_b * U_Q(iCF_S3) + ( One - x_b ) * U_K(iCF_S3), &
              x_b * U_Q(iCF_E)  + ( One - x_b ) * U_K(iCF_E) ) &
          - MinE

    IF( .NOT. f_a * f_b < 0 )THEN

      WRITE(*,'(A6,A)') &
        '', 'SolveTheta_Bisection (Euler):'
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

      f_c = eFun &
              ( x_c * U_Q(iCF_D)  + ( One - x_c ) * U_K(iCF_D),  &
                x_c * U_Q(iCF_S1) + ( One - x_c ) * U_K(iCF_S1), &
                x_c * U_Q(iCF_S2) + ( One - x_c ) * U_K(iCF_S2), &
                x_c * U_Q(iCF_S3) + ( One - x_c ) * U_K(iCF_S3), &
                x_c * U_Q(iCF_E)  + ( One - x_c ) * U_K(iCF_E) ) &
            - MinE

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
          '', 'SolveTheta_Bisection (Euler):'
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


END MODULE Euler_PositivityLimiterModule
