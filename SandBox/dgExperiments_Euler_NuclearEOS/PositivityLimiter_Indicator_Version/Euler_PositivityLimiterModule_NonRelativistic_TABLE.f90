MODULE Euler_PositivityLimiterModule_NonRelativistic_TABLE

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE UnitsModule, ONLY: &
    Gram, Centimeter, Kelvin, AtomicMassUnit, Erg, Second
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
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    Theta1, Theta2, Theta3, E_Minimum
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_Euler_NonRelativistic_TABLE
  PUBLIC :: FinalizePositivityLimiter_Euler_NonRelativistic_TABLE
  PUBLIC :: ApplyPositivityLimiter_Euler_NonRelativistic_TABLE

  REAL(DP), PARAMETER   :: Unit_D     = Gram / Centimeter**3
  REAL(DP), PARAMETER   :: Unit_T     = Kelvin
  REAL(DP), PARAMETER   :: BaryonMass = AtomicMassUnit

  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: Verbose
  INTEGER, PARAMETER    :: nPS = 7
  INTEGER               :: nPP(nPS)
  INTEGER               :: nPT
  REAL(DP)              :: Min_D, Max_D, Min_T, Max_T, &
                           Min_Y, Max_Y
  REAL(DP), ALLOCATABLE :: U_PP(:,:)

CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
    ( Min_D_Option, Max_D_Option, Min_T_Option, Max_T_Option, &
      Min_Y_Option, Max_Y_Option, UsePositivityLimiter_Option, &
      Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_D_Option, Max_D_Option, &
                                      Min_T_Option, Max_T_Option, &
                                      Min_Y_Option, Max_Y_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: i

    Min_D = - HUGE( One )
    IF( PRESENT( Min_D_Option ) ) &
      Min_D = Min_D_Option

    Max_D = - HUGE( One )
    IF( PRESENT( Max_D_Option ) ) &
      Max_D = Max_D_Option

    Min_T = - HUGE( One )
    IF( PRESENT( Min_T_Option ) ) &
      Min_T = Min_T_Option

    Max_T = - HUGE( One )
    IF( PRESENT( Max_T_Option ) ) &
      Max_T = Max_T_Option

    Min_Y = - HUGE( One )
    IF( PRESENT( Min_Y_Option ) ) &
      Min_Y = Min_Y_Option

    Max_Y = - HUGE( One )
    IF( PRESENT( Max_Y_Option ) ) &
      Max_Y = Max_Y_Option

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
      WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'Euler_PositivityLimiterModule_NonRelativistic_TABLE'
      WRITE(*,'(A2,A)') '',    '-------------------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A6,A,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter
      WRITE(*,*)
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_D = ', Min_D / Unit_D
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_D = ', Max_D / Unit_D
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_T = ', Min_T / Unit_T
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_T = ', Max_T / Unit_T
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_Y = ', Min_Y
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_Y = ', Max_Y

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

  END SUBROUTINE InitializePositivityLimiter_Euler_NonRelativistic_TABLE


  SUBROUTINE FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

    DEALLOCATE( U_PP )

  END SUBROUTINE FinalizePositivityLimiter_Euler_NonRelativistic_TABLE


  SUBROUTINE ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, ResetIndicators_Option )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL, INTENT(in), OPTIONAL :: &
      ResetIndicators_Option

    LOGICAL  :: ResetIndicators
    LOGICAL  :: NegativeStates(3)
    INTEGER  :: iX1, iX2, iX3, iCF, iP
    REAL(DP) :: Min_D_K, Max_D_K, Min_N_K, Max_N_K, &
                Min_N, Max_N, Min_E(nPT),  &
                Theta_1, Theta_2, Theta_3, Theta_P
    REAL(DP) :: U_q(nDOFX,nCF), U_K(nCF), E_K, Min_E_K(1), Min_E_PP
    REAL(DP) :: Y_PP(nPT), E_PP(nPT), Y_K

    Min_N = Min_D * Min_Y / BaryonMass
    Max_N = Max_D * Max_Y / BaryonMass

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    IF( PRESENT( ResetIndicators_Option ) )THEN
      ResetIndicators = ResetIndicators_Option
    ELSE
      ResetIndicators = .FALSE.
    END IF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      U_q(1:nDOFX,1:nCF) = U(1:nDOFX,iX1,iX2,iX3,1:nCF)

      NegativeStates = .FALSE.

      IF( ResetIndicators )THEN

        Theta1(iX1,iX2,iX3) = One
        Theta2(iX1,iX2,iX3) = One
        Theta3(iX1,iX2,iX3) = One

      END IF

      CALL ComputePointValues_Fluid( U_q, U_PP )

      ! --- Ensure Positive Mass Density ---

      Min_D_K = MINVAL( U_PP(:,iCF_D) )
      Max_D_K = MAXVAL( U_PP(:,iCF_D) )

      IF( Min_D_K < Min_D .OR. Max_D_K > Max_D )THEN

        ! --- Cell Average ---

        U_K(iCF_D) &
          = SUM( WeightsX_q(:) * U_q(:,iCF_D) * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
              / SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) )

        Theta_1 = MIN( One, &
                       ABS( ( Min_D-U_K(iCF_D) ) / ( Min_D_K-U_K(iCF_D) ) ), &
                       ABS( ( Max_D-U_K(iCF_D) ) / ( Max_D_K-U_K(iCF_D) ) ) )

        ! --- Limit Density Towards Cell Average ---

        U_q(:,iCF_D) = Theta_1 * U_q(:,iCF_D) &
                       + ( One - Theta_1 ) * U_K(iCF_D)

        ! --- Recompute Point Values ---

        CALL ComputePointValues_Fluid( U_q, U_PP )

        ! --- Flag for Negative Density ---

        NegativeStates(1) = .TRUE.

        Theta1(iX1,iX2,iX3) = MIN( Theta1(iX1,iX2,iX3), Theta_1 )

      END IF

      ! --- Ensure Bounded Electron Density ---

      Min_N_K = MINVAL( U_PP(:,iCF_Ne) )
      Max_N_K = MAXVAL( U_PP(:,iCF_Ne) )

      IF( Min_N_K < Min_N  .OR. Max_N_K > Max_N )THEN

        ! --- Cell Average ---

        U_K(iCF_Ne) &
          = SUM( WeightsX_q(:) * U_q(:,iCF_Ne) * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
              / SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) )

        Theta_2 = MIN( One, &
                       ABS( ( Min_N-U_K(iCF_Ne) ) / ( Min_N_K-U_K(iCF_Ne) ) ), &
                       ABS( ( Max_N-U_K(iCF_Ne) ) / ( Max_N_K-U_K(iCF_Ne) ) ) )

        ! --- Limit Electron Density Towards Cell Average ---

        U_q(:,iCF_Ne) = Theta_2 * U_q(:,iCF_Ne) &
                        + ( One - Theta_2 ) * U_K(iCF_Ne)

        ! --- Recompute Point Values ---

        CALL ComputePointValues_Fluid( U_q, U_PP )

        ! --- Flag for Negative Electron Density --

        NegativeStates(2) = .TRUE.

        Theta2(iX1,iX2,iX3) = MIN( Theta2(iX1,iX2,iX3), Theta_2 )

      END IF

      ! --- Ensure Positive Specific Internal Energy ---


      DO iP = 1, nPT

         CALL ComputeSpecificInternalEnergyAndElectronFraction &
                ( U_PP(iP,1:nCF), E_PP(iP), Y_PP(iP) )

         CALL ComputeSpecificInternalEnergy_TABLE &
                ( [U_PP(iP,iCF_D)], [Min_T], [Y_PP(iP)], Min_E(iP:iP) )

      END DO

       DO iP = 1, nDOFX

        E_Minimum(iP,iX1,iX2,iX3) = Min_E(iP) / ( Erg / Gram )

      END DO

      IF( ANY( E_PP(:) < Min_E(:) ) )THEN

        ! --- Cell Average ---

        DO iCF = 1, nCF

          U_K(iCF) &
            = SUM( WeightsX_q(:) * U_q(:,iCF) * G(:,iX1,iX2,iX3,iGF_SqrtGm) ) &
                / SUM( WeightsX_q(:) * G(:,iX1,iX2,iX3,iGF_SqrtGm) )

        END DO

        CALL ComputeSpecificInternalEnergyAndElectronFraction &
               ( U_K(1:nCF), E_K, Y_K )

        CALL ComputeSpecificInternalEnergy_TABLE &
               ( [U_K(iCF_D)], [Min_T], [Y_K], Min_E_K )

        DO

          IF( ALL( E_PP(:) > Min_E(:) ) )EXIT

          Theta_3 = One
          DO iP = 1, nPT

            IF( E_PP(iP) < Min_E(iP) )THEN

              CALL SolveTheta_Bisection &
                     ( U_PP(iP,1:nCF), U_K(1:nCF), Min_E(iP), &
                       Min_E_K(1), Theta_P )

              Theta_3 = MIN( Theta_3, 0.99_DP * Theta_P )

            END IF

          END DO

          ! --- Limit Towards Cell Average ---

          DO iCF = 1, nCF

            U_q(:,iCF) = Theta_3 * U_q(:,iCF) &
                         + ( One - Theta_3 ) * U_K(iCF)

          END DO

          CALL ComputePointValues_Fluid( U_q, U_PP )

          DO iP = 1, nPT

            CALL ComputeSpecificInternalEnergyAndElectronFraction &
                   ( U_PP(iP,1:nCF), E_PP(iP), Y_PP(iP) )

            CALL ComputeSpecificInternalEnergy_TABLE &
                   ( [U_PP(iP,iCF_D)], [Min_T], [Y_PP(iP)], Min_E(iP:iP) )

          END DO

        END DO

        ! --- Flag for Negative Internal Energy Density ---

        NegativeStates(1) = .FALSE.
        NegativeStates(2) = .FALSE.
        NegativeStates(3) = .TRUE.

        Theta3(iX1,iX2,iX3) = MIN( Theta3(iX1,iX2,iX3), Theta_3 )

      END IF

      IF( NegativeStates(1) )THEN

        U(1:nDOFX,iX1,iX2,iX3,iCF_D) &
          = U_q(1:nDOFX,iCF_D)

      ELSEIF( NegativeStates(2) )THEN

        U(1:nDOFX,iX1,iX2,iX3,iCF_Ne) &
          = U_q(1:nDOFX,iCF_Ne)

      ELSEIF( NegativeStates(3) )THEN

        U(1:nDOFX,iX1,iX2,iX3,1:nCF) &
          = U_q(1:nDOFX,1:nCF)

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyPositivityLimiter_Euler_NonRelativistic_TABLE


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


  SUBROUTINE ComputeSpecificInternalEnergyAndElectronFraction( U, E, Y )

    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: E, Y

    E = eFun( U(iCF_D), U(iCF_S1), U(iCF_S2), U(iCF_S3), U(iCF_E) )
    Y = BaryonMass * U(iCF_Ne) / U(iCF_D)

  END SUBROUTINE ComputeSpecificInternalEnergyAndElectronFraction


  PURE REAL(DP) ELEMENTAL FUNCTION eFun( D, S1, S2, S3, E )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, E

    eFun = ( E - Half * ( S1**2 + S2**2 + S3**2 ) / D ) / D

    RETURN
  END FUNCTION eFun


  SUBROUTINE SolveTheta_Bisection( U_Q, U_K, MinE, MinEK, Theta_P )

    REAL(DP), INTENT(in)  :: U_Q(nCF), U_K(nCF), MinE, MinEK
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
          - ( x_a * MinE + ( One - x_a ) * MinEK )

    x_b = One
    f_b = eFun &
            ( x_b * U_Q(iCF_D)  + ( One - x_b ) * U_K(iCF_D),  &
              x_b * U_Q(iCF_S1) + ( One - x_b ) * U_K(iCF_S1), &
              x_b * U_Q(iCF_S2) + ( One - x_b ) * U_K(iCF_S2), &
              x_b * U_Q(iCF_S3) + ( One - x_b ) * U_K(iCF_S3), &
              x_b * U_Q(iCF_E)  + ( One - x_b ) * U_K(iCF_E) ) &
          - ( x_b * MinE + ( One - x_b ) * MinEK )

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
            - ( x_c * MinE + ( One - x_c ) * MinEK )

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


END MODULE Euler_PositivityLimiterModule_NonRelativistic_TABLE
