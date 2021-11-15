MODULE Euler_PositivityLimiterModule_NonRelativistic_TABLE

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    AtomicMassUnit
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
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    iDF_T1, &
    iDF_T2, &
    iDF_T3, &
    iDF_MinE, &
    iDF_MaxE
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_PositivityLimiter

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_Euler_NonRelativistic_TABLE
  PUBLIC :: FinalizePositivityLimiter_Euler_NonRelativistic_TABLE
  PUBLIC :: ApplyPositivityLimiter_Euler_NonRelativistic_TABLE

  INTEGER,  PARAMETER   :: nPS          = 7
  REAL(DP), PARAMETER   :: Unit_D       = Gram / Centimeter**3
  REAL(DP), PARAMETER   :: Unit_T       = Kelvin
  REAL(DP), PARAMETER   :: BaryonMass   = AtomicMassUnit
  REAL(DP), PARAMETER   :: SafetyFactor = 0.99_DP

  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: Verbose
  INTEGER               :: nPP(nPS)
  INTEGER               :: nPT
  REAL(DP)              :: Min_D, Max_D, &
                           Min_T, Max_T, &
                           Min_Y, Max_Y
  REAL(DP), ALLOCATABLE :: U_PP(:,:)

CONTAINS

  !> @param Min_1 Minimum table Density
  !> @param Min_2 Minimum table Temperature
  !> @param Min_3 Minimum table Electron Fraction
  SUBROUTINE InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
    ( UsePositivityLimiter_Option, Verbose_Option, &
      Min_1_Option, Min_2_Option, Min_3_Option, &
      Max_1_Option, Max_2_Option, Max_3_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Max_1_Option, &
                                      Min_2_Option, Max_2_Option, &
                                      Min_3_Option, Max_3_Option

    INTEGER :: i

    UsePositivityLimiter = .TRUE.
    IF( PRESENT( UsePositivityLimiter_Option ) ) &
      UsePositivityLimiter = UsePositivityLimiter_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    Min_D = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_D = Min_1_Option

    Max_D = - HUGE( One )
    IF( PRESENT( Max_1_Option ) ) &
      Max_D = Max_1_Option

    Min_T = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_T = Min_2_Option

    Max_T = - HUGE( One )
    IF( PRESENT( Max_2_Option ) ) &
      Max_T = Max_2_Option

    Min_Y = - HUGE( One )
    IF( PRESENT( Min_3_Option ) ) &
      Min_Y = Min_3_Option

    Max_Y = - HUGE( One )
    IF( PRESENT( Max_3_Option ) ) &
      Max_Y = Max_3_Option

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A2,A6,A)') '', &
        'INFO: ', 'Euler_PositivityLimiterModule_NonRelativistic_TABLE'
      WRITE(*,'(A2,A)')    '', &
        '---------------------------------------------------------'
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
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) ! Diagnostic Fluid

    LOGICAL  :: NegativeStates
    INTEGER  :: iX1, iX2, iX3, iCF, iQ, iP
    INTEGER  :: ITERATION
    REAL(DP) :: Min_D_K, Max_D_K, Min_N_K
    REAL(DP) :: Min_Y_K, Max_Y_K, Min_E_K
    REAL(DP) :: Min_N, Max_N, Min_E(nPT), Max_E(nPT)
    REAL(DP) :: Theta_1, Theta_2, Theta_3, Theta_P, Alpha
    REAL(DP) :: D_P, N_P
    REAL(DP) :: G_q(nDOFX,nGF)
    REAL(DP) :: U_q(nDOFX,nCF), U_K(nCF)
    REAL(DP) :: Y_PP(nPT), E_PP(nPT), Y_K, E_K

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

    Min_N = Min_D * Min_Y / BaryonMass
    Max_N = Max_D * Max_Y / BaryonMass

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      G_q(1:nDOFX,1:nGF) = G(1:nDOFX,iX1,iX2,iX3,1:nGF)
      U_q(1:nDOFX,1:nCF) = U(1:nDOFX,iX1,iX2,iX3,1:nCF)

      NegativeStates = .FALSE.

      ! --- Compute Fluid Cell Averages ---

      DO iCF = 1, nCF

        U_K(iCF) &
          = SUM( WeightsX_q(:) * U_q(:,iCF) * G_q(:,iGF_SqrtGm) ) &
              / SUM( WeightsX_q(:) * G_q(:,iGF_SqrtGm) )

      END DO

      Y_K = BaryonMass * U_K(iCF_Ne) / U_K(iCF_D)

      IF( U_K(iCF_D) < Min_D .OR. U_K(iCF_D) > Max_D )THEN

        ! --- Cell Average Density Outside Bounds ---

        PRINT*
        PRINT*, "Cell Average Density Outside Bounds"
        PRINT*
        PRINT*, "iX1,iX2,iX3 = ", iX1, iX2, iX3
        PRINT*
        PRINT*, "Min_D, Max_D, D_K = ", &
          Min_D / Unit_D, Max_D / Unit_D, U_K(iCF_D) / Unit_D
        DO iQ = 1, nDOFX
          PRINT*, "iQ, D(iQ) = ", iQ, U_q(iQ,iCF_D) / Unit_D
        END DO
        PRINT*
        PRINT*, "Y_K = ", Y_K
        PRINT*

        ! --- Force Cell Averages Inside Table Bounds ---

        IF( U_K(iCF_D) < Min_D )THEN

          U_q(1:nDOFX,iCF_D) = ( One + 1.0d-6 ) * Min_D

        ELSE

          U_q(1:nDOFX,iCF_D) = ( One - 1.0d-6 ) * Max_D

        END IF

        ! --- Reset Electron Density by Preserving Cell-Averaged Ye ---

        U_q(1:nDOFX,iCF_Ne) = Y_K * U_q(1:nDOFX,iCF_D) / BaryonMass

        ! --- Recompute Cell Averages ---

        U_K(iCF_D ) &
          = SUM( WeightsX_q(:) * U_q(:,iCF_D ) * G_q(:,iGF_SqrtGm) ) &
              / SUM( WeightsX_q(:) * G_q(:,iGF_SqrtGm) )

        U_K(iCF_Ne) &
          = SUM( WeightsX_q(:) * U_q(:,iCF_Ne) * G_q(:,iGF_SqrtGm) ) &
              / SUM( WeightsX_q(:) * G_q(:,iGF_SqrtGm) )

        NegativeStates = .TRUE.

        ! --- Flag That Cell Average Violated Bounds ---

        D(:,iX1,iX2,iX3,iDF_T1) = - One

      END IF

      CALL ComputePointValues_Fluid( U_q, U_PP )

      ! -------------------------------------------------------------------
      ! --- Step 1 --------------------------------------------------------
      ! --- Ensure Bounded Mass Density and Positive Electron Density -----
      ! -------------------------------------------------------------------

      Min_D_K = MINVAL( U_PP(:,iCF_D)  ) ! --- Minimum D  in Element
      Max_D_K = MAXVAL( U_PP(:,iCF_D)  ) ! --- Maximum D  in Element
      Min_N_K = MINVAL( U_PP(:,iCF_Ne) ) ! --- Minimum Ne in Element

      IF( ANY( [ Min_D_K - Min_D  , &
                 Max_D   - Max_D_K, &
                 Min_N_K - Min_N ] < Zero ) ) &
      THEN

        PRINT*, "Step 1"

        PRINT*, "Min_D, Min_D_K, D_K = ", Min_D, Min_D_K, U_K(iCF_D)
        PRINT*, "Max_D, Max_D_K, D_K = ", Max_D, Max_D_K, U_K(iCF_D)
        PRINT*, "Min_N, Min_N_K, N_K = ", Min_N, Min_N_K, U_K(iCF_Ne)

        Theta_1 = One

        DO iP = 1, nPT

          Theta_P = One

          D_P = U_PP(iP,iCF_D )
          N_P = U_PP(iP,iCF_Ne)
          DO WHILE( ANY( [D_P-Min_D,Max_D-D_P,N_P-Min_N] < Zero ) )

            IF( Theta_P > 1.0d-2 )THEN

              Theta_P = 0.95_DP * Theta_P ! --- Shrink Theta_P

            ELSE

              Theta_P = Zero

            END IF

            D_P = (One-Theta_P) * U_K(iCF_D ) + Theta_P * U_PP(iP,iCF_D )
            N_P = (One-Theta_P) * U_K(iCF_Ne) + Theta_P * U_PP(iP,iCF_Ne)

          END DO

          Theta_1 = MIN( Theta_1, SafetyFactor * Theta_P )

        END DO

        ! --- Limit Mass and Electron Density Towards Cell Average ---

        U_q(:,iCF_D ) &
          = Theta_1 * U_q(:,iCF_D ) + ( One - Theta_1 ) * U_K(iCF_D )

        U_q(:,iCF_Ne) &
          = Theta_1 * U_q(:,iCF_Ne) + ( One - Theta_1 ) * U_K(iCF_Ne)

        ! --- Recompute Point Values ---

        CALL ComputePointValues_Fluid( U_q, U_PP )

        NegativeStates = .TRUE.

        ! --- Set Diagnostic Fields Theta 1 ---

        D(:,iX1,iX2,iX3,iDF_T1) &
          = MIN( MINVAL( D(:,iX1,iX2,iX3,iDF_T1) ), Theta_1 )

        IF( ANY( U_PP(:,iCF_D) < Min_D ) )THEN

          PRINT*, "Failed Step 1"
          STOP

        END IF

        PRINT*, "Passed Step 1.  Theta_1 = ", Theta_1

      END IF

      ! -------------------------------------------------------------------
      ! --- Step 2 --------------------------------------------------------
      ! --- Ensure Bounded Electron Fraction ------------------------------
      ! -------------------------------------------------------------------

      Min_Y_K = BaryonMass * MINVAL( U_PP(:,iCF_Ne) / U_PP(:,iCF_D) )
      Max_Y_K = BaryonMass * MAXVAL( U_PP(:,iCF_Ne) / U_PP(:,iCF_D) )

      IF( Min_Y_K < Min_Y .OR. Max_Y_K > Max_Y )THEN

        PRINT*, "Step 2"

        PRINT*, "Min_Y, Min_Y_K, Y_K = ", Min_Y, Min_Y_K, Y_K
        PRINT*, "Max_Y, Max_Y_K, Y_K = ", Max_Y, Max_Y_K, Y_K

        Alpha = MIN( One, &
                     ABS( ( Min_Y - Y_K ) / ( Min_Y_K - Y_K ) ), &
                     ABS( ( Max_Y - Y_K ) / ( Max_Y_K - Y_K ) ) )

        PRINT*, "Alpha   = ", Alpha

        Theta_2 &
          = SafetyFactor * Alpha * U_K(iCF_D) &
            / ( Alpha * U_K(iCF_D) + (One-Alpha) * MAXVAL(U_PP(:,iCF_D)) )

        PRINT*, "Theta_2 = ", Theta_2

        ! --- Limit Mass and Electron Density Towards Cell Average ---

        U_q(:,iCF_D ) &
          = Theta_2 * U_q(:,iCF_D ) + ( One - Theta_2 ) * U_K(iCF_D )

        U_q(:,iCF_Ne) &
          = Theta_2 * U_q(:,iCF_Ne) + ( One - Theta_2 ) * U_K(iCF_Ne)

        ! --- Recompute Point Values ---

        CALL ComputePointValues_Fluid( U_q, U_PP )

        NegativeStates = .TRUE.

        D(:,iX1,iX2,iX3,iDF_T2) &
          = MIN( MINVAL( D(:,iX1,iX2,iX3,iDF_T2) ), Theta_2 )

      END IF

      ! -------------------------------------------------------------------
      ! --- Step 3 --------------------------------------------------------
      ! --- Ensure Bounded Specific Internal Energy -----------------------
      ! -------------------------------------------------------------------

      DO iP = 1, nPT

         CALL ComputeSpecificInternalEnergyAndElectronFraction &
                ( U_PP(iP,1:nCF), E_PP(iP), Y_PP(iP) )

         CALL ComputeSpecificInternalEnergy_TABLE &
                ( U_PP(iP,iCF_D), Min_T, Y_PP(iP), Min_E(iP) )

         CALL ComputeSpecificInternalEnergy_TABLE &
                ( U_PP(iP,iCF_D), Max_T, Y_PP(iP), Max_E(iP) )

      END DO

      ITERATION = 0
      DO WHILE( ANY( E_PP < Min_E ) )

        ITERATION = ITERATION + 1

        DO iP = 1, nDOFX ! --- Diagnostics (excludes points on interfaces)

          D(iP,iX1,iX2,iX3,iDF_MinE) = Min_E(iP)
          D(iP,iX1,iX2,iX3,iDF_MaxE) = Max_E(iP)

        END DO

        CALL ComputeSpecificInternalEnergyAndElectronFraction &
               ( U_K(1:nCF), E_K, Y_K )

        CALL ComputeSpecificInternalEnergy_TABLE &
               ( U_K(iCF_D), Min_T, Y_K, Min_E_K )

        IF( ITERATION .LT. 10 )THEN

          Theta_3 = One
          DO iP = 1, nPT

            IF( E_PP(iP) < Min_E(iP) )THEN

              CALL SolveTheta_Bisection &
                     ( U_PP(iP,1:nCF), U_K(1:nCF), Min_E(iP), Min_E_K, Theta_P )

              Theta_3 = MIN( Theta_3, SafetyFactor * Theta_P )

            END IF

          END DO

        ELSE

          Theta_3 = Zero

        END IF

        ! --- Limit Towards Cell Average ---

        DO iCF = 1, nCF

          U_q(:,iCF) = Theta_3 * U_q(:,iCF) + ( One - Theta_3 ) * U_K(iCF)

        END DO

        NegativeStates = .TRUE.

        D(:,iX1,iX2,iX3,iDF_T3) &
          = MIN( MINVAL( D(:,iX1,iX2,iX3,iDF_T3) ), Theta_3 )

        ! --- Recompute Point Values ---

        CALL ComputePointValues_Fluid( U_q, U_PP )

        DO iP = 1, nPT

          CALL ComputeSpecificInternalEnergyAndElectronFraction &
                 ( U_PP(iP,1:nCF), E_PP(iP), Y_PP(iP) )

          CALL ComputeSpecificInternalEnergy_TABLE &
                 ( U_PP(iP,iCF_D), Min_T, Y_PP(iP), Min_E(iP) )

          CALL ComputeSpecificInternalEnergy_TABLE &
                 ( U_PP(iP,iCF_D), Max_T, Y_PP(iP), Max_E(iP) )

        END DO

      END DO ! --- WHILE( ANY( E_PP < Min_E ) )

      IF( NegativeStates )THEN

        U(1:nDOFX,iX1,iX2,iX3,1:nCF) = U_q(1:nDOFX,1:nCF)

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

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


  SUBROUTINE SolveTheta_Bisection( U_Q, U_K, MinE, MinE_K, Theta_P )

    REAL(DP), INTENT(in)  :: U_Q(nCF), U_K(nCF), MinE, MinE_K
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
          - ( x_a * MinE + ( One - x_a ) * MinE_K )

    x_b = One
    f_b = eFun &
            ( x_b * U_Q(iCF_D)  + ( One - x_b ) * U_K(iCF_D),  &
              x_b * U_Q(iCF_S1) + ( One - x_b ) * U_K(iCF_S1), &
              x_b * U_Q(iCF_S2) + ( One - x_b ) * U_K(iCF_S2), &
              x_b * U_Q(iCF_S3) + ( One - x_b ) * U_K(iCF_S3), &
              x_b * U_Q(iCF_E)  + ( One - x_b ) * U_K(iCF_E) ) &
          - ( x_b * MinE + ( One - x_b ) * MinE_K )

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
            - ( x_c * MinE + ( One - x_c ) * MinE_K )

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
