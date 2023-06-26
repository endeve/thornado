MODULE Euler_PositivityLimiterModule_NonRelativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One
  USE UtilitiesModule, ONLY: &
    IsCornerCell
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_PositivityLimiter
  USE UnitsModule, ONLY: &
    UnitsDisplay

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_Euler_NonRelativistic_IDEAL
  PUBLIC :: FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL
  PUBLIC :: ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL


  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: Verbose
  INTEGER, PARAMETER    :: nPS = 7
  INTEGER               :: nPP(nPS)
  INTEGER               :: nPT
  REAL(DP)              :: Min_1, Min_2
  REAL(DP)              :: D_Min_Euler_PL, IntE_Min_Euler_PL
  REAL(DP), ALLOCATABLE :: U_PP(:,:), G_PP(:,:)

CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler_NonRelativistic_IDEAL &
    ( UsePositivityLimiter_Option, Verbose_Option, Min_1_Option, Min_2_Option, &
      D_Min_Euler_PL_Option, IntE_Min_Euler_PL_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    REAL(DP), INTENT(in), OPTIONAL :: D_Min_Euler_PL_Option
    REAL(DP), INTENT(in), OPTIONAL :: IntE_Min_Euler_PL_Option

    INTEGER :: i

    UsePositivityLimiter = .TRUE.
    IF( PRESENT( UsePositivityLimiter_Option ) ) &
      UsePositivityLimiter = UsePositivityLimiter_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    D_Min_Euler_PL = Zero
    IF( PRESENT( D_Min_Euler_PL_Option ) ) &
      D_Min_Euler_PL = D_Min_Euler_PL_Option

    IntE_Min_Euler_PL = Zero
    IF( PRESENT( IntE_Min_Euler_PL_Option ) ) &
      IntE_Min_Euler_PL = IntE_Min_Euler_PL_Option

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') &
        '  INFO: InitializePositivityLimiter_Euler_NonRelativistic_IDEAL:'
      WRITE(*,'(A)') &
        '  --------------------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A6,A,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter
      WRITE(*,*)
      WRITE(*,'(A6,A20,ES12.4E3)') '', 'Min_1 = ', Min_1
      WRITE(*,'(A6,A20,ES12.4E3)') '', 'Min_2 = ', Min_2
      WRITE(*,'(A6,A20,ES12.4E3,A,A)') '', 'D_Min_Euler_PL = ', &
        D_Min_Euler_PL / UnitsDisplay % MassDensityUnit, ' ', &
        UnitsDisplay % MassDensityLabel
      WRITE(*,'(A6,A20,ES12.4E3,A,A)') '', 'IntE_Min_Euler_PL = ', &
        IntE_Min_Euler_PL / UnitsDisplay % EnergyDensityUnit, ' ', &
        UnitsDisplay % EnergyDensityLabel
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

    ! --- Geometry in Positive Points ---

    ALLOCATE( G_PP(nPT,nGF) )

  END SUBROUTINE InitializePositivityLimiter_Euler_NonRelativistic_IDEAL


  SUBROUTINE FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL

    DEALLOCATE( U_PP )
    DEALLOCATE( G_PP )

  END SUBROUTINE FinalizePositivityLimiter_Euler_NonRelativistic_IDEAL


  SUBROUTINE ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL &
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
    REAL(DP) :: Min_D, Min_E
    REAL(DP) :: U_q(nDOFX,nCF), G_q(nDOFX,nGF), U_K(nCF), &
                IntE(nPT), IntE_K_P(nPT), IntE_K_P_Min

    REAL(DP), PARAMETER :: Alpha = 1.01_DP

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( IsCornerCell( iX_B1, iX_E1, iX1, iX2, iX3 ) ) CYCLE

      U_q(1:nDOFX,1:nCF) = U(1:nDOFX,iX1,iX2,iX3,1:nCF)
      G_q(1:nDOFX,1:nGF) = G(1:nDOFX,iX1,iX2,iX3,1:nGF)

      NegativeStates = .FALSE.

      DO iCF = 1, nCF
        CALL ComputePointValues( U_q(:,iCF), U_PP(:,iCF) )
      END DO

      DO iGF = iGF_h_1, iGF_h_3
        CALL ComputePointValues( G_q(:,iGF), G_PP(:,iGF) )
      END DO

      CALL ComputeGeometryX_FromScaleFactors( G_PP )

      ! --- Ensure Positive Mass Density ---

      Min_K = MINVAL( U_PP(:,iCF_D) )

      ! --- Cell Average ---

      U_K(iCF_D) &
        = SUM( WeightsX_q(:) * U_q(:,iCF_D) * G_q(:,iGF_SqrtGm) ) &
            / SUM( WeightsX_q(:) * G_q(:,iGF_SqrtGm) )

      Min_D = Min_1 * U_K(iCF_D)

      IF( Min_K < Min_D )THEN

        Theta_1 = MIN( One, ( U_K(iCF_D) - Min_D ) / ( U_K(iCF_D) - Min_K ) )

        ! --- Limit Density Towards Cell Average ---

        U_q(:,iCF_D) = Theta_1 * U_q(:,iCF_D) &
                       + ( One - Theta_1 ) * U_K(iCF_D)

        ! --- Recompute Point Values ---

        CALL ComputePointValues( U_q(1:nDOFX,iCF_D), U_PP(1:nPT,iCF_D) )

        ! --- Flag for Negative Density ---

        NegativeStates(1) = .TRUE.

      END IF

      ! --- Ensure Positive Internal Energy Density ---

      DO iCF = 1, nCF

        U_K(iCF) &
          = SUM( WeightsX_q(:) * U_q(:,iCF) * G_q(:,iGF_SqrtGm) ) &
              / SUM( WeightsX_q(:) * G_q(:,iGF_SqrtGm) )

      END DO

      DO iP = 1, nPT

        IntE_K_P(iP) &
          = U_K(iCF_E) &
              - Half * ( U_K(iCF_S1)**2 / G_PP(iP,iGF_Gm_dd_11) &
                       + U_K(iCF_S2)**2 / G_PP(iP,iGF_Gm_dd_22) &
                       + U_K(iCF_S3)**2 / G_PP(iP,iGF_Gm_dd_33) ) / U_K(iCF_D)

      END DO

      IntE_K_P_Min = MINVAL( IntE_K_P )

      IF( IntE_K_P_Min .LT. IntE_Min_Euler_PL )THEN

        PRINT *, 'iX1, iX2, iX3, IntE (Old), dIntE: ', &
                  iX1, iX2, iX3, IntE_K_P_Min, IntE_Min_Euler_PL - IntE_K_P_Min

        U_K(iCF_E) = U_K(iCF_E) + Alpha * ( IntE_Min_Euler_PL - IntE_K_P_Min )

        DO iCF = 1, nCF

          U(1:nDOFX,iX1,iX2,iX3,iCF) = U_K(iCF)

        END DO

        CYCLE

      END IF ! ANY( IntE_K_P .LT. IntE_Min_Euler_PL )

      CALL ComputeInternalEnergyDensity &
             ( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), IntE(1:nPT) )

      Min_E = Min_2 * U_K(iCF_E)

      IF( ANY( IntE(:) < Min_E ) )THEN

        Theta_2 = One
        DO iP = 1, nPT

          IF( IntE(iP) < Min_E )THEN

            CALL SolveTheta_Bisection &
                   ( U_PP(iP,1:nCF), G_PP(iP,1:nGF), U_K(1:nCF), &
                     Min_E, Theta_P )

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

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

  END SUBROUTINE ApplyPositivityLimiter_Euler_NonRelativistic_IDEAL


  SUBROUTINE ComputePointValues( X_q, X_p )

    REAL(DP), INTENT(in)  :: X_q(nDOFX)
    REAL(DP), INTENT(out) :: X_p(nPT)

    INTEGER :: iOS

    X_p(1:nDOFX) = X_q(1:nDOFX)

    IF( SUM( nPP(2:3) ) > 0 )THEN

      ! --- Points of Faces with Normal in X1 Direction ---

      iOS = nPP(1)

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               X_q(1:nDOFX), 1, Zero, X_p(iOS+1:iOS+nDOFX_X1), 1 )

      iOS = iOS + nPP(2)

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               X_q(1:nDOFX), 1, Zero, X_p(iOS+1:iOS+nDOFX_X1), 1 )

    END IF

    IF( SUM( nPP(4:5) ) > 0 )THEN

      ! --- Points of Faces with Normal in X2 Direction ---

      iOS = SUM( nPP(1:3) )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               X_q(1:nDOFX), 1, Zero, X_p(iOS+1:iOS+nDOFX_X2), 1 )

      iOS = iOS + nPP(4)

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               X_q(1:nDOFX), 1, Zero, X_p(iOS+1:iOS+nDOFX_X2), 1 )

    END IF

    IF( SUM( nPP(6:7) ) > 0 )THEN

      ! --- Points of Faces with Normal in X3 Direction ---

      iOS = SUM( nPP(1:5) )

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
               X_q(1:nDOFX), 1, Zero, X_p(iOS+1:iOS+nDOFX_X3), 1 )

      iOS = iOS + nPP(6)

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Up, nDOFX_X3, &
               X_q(1:nDOFX), 1, Zero, X_p(iOS+1:iOS+nDOFX_X3), 1 )

    END IF

  END SUBROUTINE ComputePointValues


  SUBROUTINE ComputeInternalEnergyDensity( N, U, G, IntE )

    INTEGER,  INTENT(in)  :: N
    REAL(DP), INTENT(in)  :: U(N,nCF), G(N,nGF)
    REAL(DP), INTENT(out) :: IntE(N)

    IntE &
      = eFun( U(:,iCF_D), U(:,iCF_S1), U(:,iCF_S2), U(:,iCF_S3), U(:,iCF_E), &
              G(:,iGF_Gm_dd_11), G(:,iGF_Gm_dd_22), G(:,iGF_Gm_dd_33) )

  END SUBROUTINE ComputeInternalEnergyDensity


  REAL(DP) ELEMENTAL FUNCTION eFun( D, S1, S2, S3, E, Gm11, Gm22, Gm33 )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, E
    REAL(DP), INTENT(in) :: Gm11, Gm22, Gm33

    eFun = E - Half * ( S1**2 / Gm11 + S2**2 / Gm22 + S3**2 / Gm33 ) / D

    RETURN
  END FUNCTION eFun


  SUBROUTINE SolveTheta_Bisection( U_Q, G_Q, U_K, MinE, Theta_P )

    REAL(DP), INTENT(in)  :: U_Q(nCF), G_Q(nGF), U_K(nCF), MinE
    REAL(DP), INTENT(out) :: Theta_P

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = eFun &
            ( x_a * U_Q(iCF_D)  + ( One - x_a ) * U_K(iCF_D),         &
              x_a * U_Q(iCF_S1) + ( One - x_a ) * U_K(iCF_S1),        &
              x_a * U_Q(iCF_S2) + ( One - x_a ) * U_K(iCF_S2),        &
              x_a * U_Q(iCF_S3) + ( One - x_a ) * U_K(iCF_S3),        &
              x_a * U_Q(iCF_E)  + ( One - x_a ) * U_K(iCF_E),         &
              G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33) ) &
          - MinE

    x_b = One
    f_b = eFun &
            ( x_b * U_Q(iCF_D)  + ( One - x_b ) * U_K(iCF_D),         &
              x_b * U_Q(iCF_S1) + ( One - x_b ) * U_K(iCF_S1),        &
              x_b * U_Q(iCF_S2) + ( One - x_b ) * U_K(iCF_S2),        &
              x_b * U_Q(iCF_S3) + ( One - x_b ) * U_K(iCF_S3),        &
              x_b * U_Q(iCF_E)  + ( One - x_b ) * U_K(iCF_E),         &
              G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33) ) &
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
              ( x_c * U_Q(iCF_D)  + ( One - x_c ) * U_K(iCF_D),         &
                x_c * U_Q(iCF_S1) + ( One - x_c ) * U_K(iCF_S1),        &
                x_c * U_Q(iCF_S2) + ( One - x_c ) * U_K(iCF_S2),        &
                x_c * U_Q(iCF_S3) + ( One - x_c ) * U_K(iCF_S3),        &
                x_c * U_Q(iCF_E)  + ( One - x_c ) * U_K(iCF_E),         &
                G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33) ) &
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


END MODULE Euler_PositivityLimiterModule_NonRelativistic_IDEAL
