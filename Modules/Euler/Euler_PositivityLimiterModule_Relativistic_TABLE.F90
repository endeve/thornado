MODULE Euler_PositivityLimiterModule_Relativistic_TABLE

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    SqrtTiny
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
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputePressureFromPrimitive_TABLE, &
    ComputeTemperatureFromPressure_TABLE, &
    ComputeSpecificInternalEnergy_TABLE
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    AtomicMassUnit
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_PositivityLimiter
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_Euler_Relativistic_TABLE
  PUBLIC :: FinalizePositivityLimiter_Euler_Relativistic_TABLE
  PUBLIC :: ApplyPositivityLimiter_Euler_Relativistic_TABLE

  LOGICAL               :: UsePositivityLimiter
  INTEGER, PARAMETER    :: nPS = 7  ! Number of Positive Point Sets
  REAL(DP), PARAMETER   :: Unit_D       = Gram / Centimeter**3
  REAL(DP), PARAMETER   :: Unit_T       = Kelvin
  REAL(DP), PARAMETER   :: BaryonMass   = AtomicMassUnit
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total number of Positive Points
  REAL(DP)              :: Min_D, Max_D, Min_T, Max_T, Min_Ye, Max_Ye
  REAL(DP), ALLOCATABLE :: G_PP(:,:), U_PP(:,:), P_PP(:,:)


CONTAINS


  !> @param Min_1 Minimum table Density
  !> @param Min_2 Minimum table Temperature
  !> @param Min_3 Minimum table Electron Fraction
  SUBROUTINE InitializePositivityLimiter_Euler_Relativistic_TABLE &
    ( UsePositivityLimiter_Option, Verbose_Option, &
      Min_1_Option, Min_2_Option, Min_3_Option, &
      Max_1_Option, Max_2_Option, Max_3_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option, &
                                      Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Max_1_Option, &
                                      Min_2_Option, Max_2_Option, &
                                      Min_3_Option, Max_3_Option

    INTEGER :: iDim
    LOGICAL :: Verbose

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

    Min_Ye = - HUGE( One )
    IF( PRESENT( Min_3_Option ) ) &
      Min_Ye = Min_3_Option

    Max_Ye = - HUGE( One )
    IF( PRESENT( Max_3_Option ) ) &
      Max_Ye = Max_3_Option

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Positivity Limiter (Euler, Relativistic, TABLE)'
      WRITE(*,'(A)') &
        '    -----------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A6,A,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter
      WRITE(*,*)
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_D = ', Min_D / Unit_D
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_D = ', Max_D / Unit_D
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_T = ', Min_T / Unit_T
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_T = ', Max_T / Unit_T
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_Y = ', Min_Ye
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_Y = ', Max_Ye

    END IF

    nPP(1:nPS) = 0
    nPP(1)     = PRODUCT( nNodesX(1:3) )

    DO iDim = 1, 3

      IF( nNodesX(iDim) > 1 )THEN

        nPP(2*iDim:2*iDim+1) &
          = PRODUCT( nNodesX(1:3), MASK = [1,2,3] .NE. iDim )

      END IF

    END DO

    nPT = SUM( nPP(1:nPS) )

    ALLOCATE( G_PP(nPT,nGF) )
    ALLOCATE( U_PP(nPT,nCF) )
    ALLOCATE( P_PP(nPT,nPF) )

  END SUBROUTINE InitializePositivityLimiter_Euler_Relativistic_TABLE


  SUBROUTINE FinalizePositivityLimiter_Euler_Relativistic_TABLE

    DEALLOCATE( P_PP )
    DEALLOCATE( U_PP )
    DEALLOCATE( G_PP )

  END SUBROUTINE FinalizePositivityLimiter_Euler_Relativistic_TABLE


  SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_TABLE &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF
    REAL(DP) :: U_q(nDOFX,nCF), U_K(nCF), SqrtGm(nDOFX)

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( U_q, SqrtGm, U_K, U_PP )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      U_q = U(:,iX1,iX2,iX3,:)

      SqrtGm = G(:,iX1,iX2,iX3,iGF_SqrtGm)

      DO iCF = 1, nCF

        U_K(iCF) &
          = SUM( WeightsX_q * U_q(:,iCF) * SqrtGm ) &
              / SUM( WeightsX_q * SqrtGm )

        CALL ComputePointValues( U_q(:,iCF), U_PP(:,iCF) )

      END DO

      IF( ANY( U_PP(:,iCF_E) .LT. Zero ) )THEN

        DO iCF = 1, nCF

          U(:,iX1,iX2,iX3,iCF) = U_K(iCF)

        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

  END SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_TABLE


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


  REAL(DP) ELEMENTAL FUNCTION qFun &
    ( D, S1, S2, S3, tau, Gm11, Gm22, Gm33, Min_ESq )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, tau, Gm11, Gm22, Gm33, Min_ESq

    qFun = tau + D &
             - SQRT( D**2 + ( S1**2 / Gm11 + S2**2 / Gm22 + S3**2 / Gm33 ) &
                       + Min_ESq )

    RETURN
  END FUNCTION qFun


END MODULE Euler_PositivityLimiterModule_Relativistic_TABLE

