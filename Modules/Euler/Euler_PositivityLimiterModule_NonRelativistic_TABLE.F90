MODULE Euler_PositivityLimiterModule_NonRelativistic_TABLE

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    SqrtTiny
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    AtomicMassUnit
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
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
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
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
    Timer_Euler_PositivityLimiter, &
    Timer_Euler_PL_LimitCells, &
    Timer_Euler_PL_CopyIn, &
    Timer_Euler_PL_CopyOut, &
    Timer_Euler_PL_Permute, &
    Timer_Euler_PL_Integrate

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
  INTEGER               :: nX_G
  INTEGER               :: nPP(nPS)
  INTEGER               :: nPT
  REAL(DP)              :: Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y, Min_N, Max_N
  REAL(DP), ALLOCATABLE :: InterpMat(:,:)
  REAL(DP), ALLOCATABLE :: U_PP(:,:)

#if   defined( THORNADO_OMP_OL )
  !$OMP DECLARE TARGET( nPT )
#elif defined( THORNADO_OACC   )
  !$ACC DECLARE CREATE( Min_D, Max_D, Min_T, Max_T, &
  !$ACC                 Min_Y, Max_Y, Min_N, Max_N, &
  !$ACC                 nPT, InterpMat )
#endif

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

    INTEGER :: i, iNodeX, iOS

    IF( PRESENT( UsePositivityLimiter_Option ) )THEN
      UsePositivityLimiter = UsePositivityLimiter_Option
    ELSE
      UsePositivityLimiter = .TRUE.
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( PRESENT( Min_1_Option ) )THEN
      Min_D = Min_1_Option
    ELSE
      Min_D = - HUGE( One )
    END IF

    IF( PRESENT( Max_1_Option ) )THEN
      Max_D = Max_1_Option
    ELSE
      Max_D = + HUGE( One )
    END IF

    IF( PRESENT( Min_2_Option ) )THEN
      Min_T = Min_2_Option
    ELSE
      Min_T = - HUGE( One )
    END IF

    IF( PRESENT( Max_2_Option ) )THEN
      Max_T = Max_2_Option
    ELSE
      Max_T = + HUGE( One )
    END IF

    IF( PRESENT( Min_3_Option ) )THEN
      Min_Y = Min_3_Option
    ELSE
      Min_Y = - HUGE( One )
    END IF

    IF( PRESENT( Max_3_Option ) )THEN
      Max_Y = Max_3_Option
    ELSE
      Max_Y = + HUGE( One )
    END IF

    Min_N = Min_D * Min_Y / BaryonMass
    Max_N = Max_D * Max_Y / BaryonMass

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

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE TO( Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE DEVICE   ( Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y, Min_N, Max_N )
#endif

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

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET UPDATE TO( nPT )
#elif defined( THORNADO_OACC   )
    !$ACC UPDATE DEVICE( nPT )
#endif

    ! --- Interpolation Matrix ---

    ALLOCATE( InterpMat(nPT,nDOFX) )

    InterpMat = Zero
    DO iNodeX = 1, nDOFX

      InterpMat(iNodeX,iNodeX) = One

      IF( SUM( nPP(2:3) ) > 0 )THEN

        iOS = nPP(1)
        InterpMat(iOS+1:iOS+nDOFX_X1,iNodeX) = LX_X1_Dn(1:nDOFX_X1,iNodeX)

        iOS = iOS + nPP(2)
        InterpMat(iOS+1:iOS+nDOFX_X1,iNodeX) = LX_X1_Up(1:nDOFX_X1,iNodeX)

      END IF

      IF( SUM( nPP(4:5) ) > 0 )THEN

        iOS = SUM( nPP(1:3) )
        InterpMat(iOS+1:iOS+nDOFX_X2,iNodeX) = LX_X2_Dn(1:nDOFX_X2,iNodeX)

        iOS = iOS + nPP(4)
        InterpMat(iOS+1:iOS+nDOFX_X2,iNodeX) = LX_X2_Up(1:nDOFX_X2,iNodeX)

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        iOS = SUM( nPP(1:5) )
        InterpMat(iOS+1:iOS+nDOFX_X3,iNodeX) = LX_X3_Dn(1:nDOFX_X3,iNodeX)

        iOS = iOS + nPP(6)
        InterpMat(iOS+1:iOS+nDOFX_X3,iNodeX) = LX_X3_Up(1:nDOFX_X3,iNodeX)

      END IF

    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: InterpMat )
#elif defined( THORNADO_OACC   )
    PRINT*,"InterpMat (CPU) = ", InterpMat
    !$ACC ENTER DATA &
    !$ACC COPYIN( InterpMat )
#endif

    ! --- Conserved Variables in Positive Points ---

    ALLOCATE( U_PP(nPT,nCF) )

  END SUBROUTINE InitializePositivityLimiter_Euler_NonRelativistic_TABLE


  SUBROUTINE FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

    IF( ALLOCATED( InterpMat ) )THEN

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: InterpMat )
#elif defined( THORNADO_OACC   )
      !$ACC EXIT DATA &
      !$ACC DELETE( InterpMat )
#endif

      DEALLOCATE( InterpMat )

    END IF

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
    INTEGER  :: iX1, iX2, iX3, iCF, iP, iNodeX
    INTEGER  :: ITERATION
    REAL(DP) :: Min_D_K, Max_D_K, Min_N_K
    REAL(DP) :: Min_Y_K, Max_Y_K, Min_E_K
    REAL(DP) :: Min_E(nPT), Max_E(nPT)
    REAL(DP) :: Min_Eps_P(nPT), Max_Eps_P(nPT), Min_Eps_K(nPT)
    REAL(DP) :: Theta_1, Theta_2, Theta_3, Theta_P, Alpha
    REAL(DP) :: D_P, N_P
    REAL(DP) :: G_q(nDOFX,nGF)
    REAL(DP) :: U_q(nDOFX,nCF), U_K(nCF)
    REAL(DP) :: Eps_P(nPT), Eps_K(nPT), Ye_P(nPT), Ye_K(nPT)
    REAL(DP) :: Y_PP(nPT), E_PP(nPT), Y_K, E_K

    REAL(DP) :: &
      G_Q_h_d_1 (nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      G_Q_h_d_2 (nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      G_Q_h_d_3 (nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      G_Q_SqrtGm(nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP) :: &
      G_P_Gm_dd_11(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      G_P_Gm_dd_22(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      G_P_Gm_dd_33(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP) :: &
      U_K_D (iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_K_S1(iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_K_S2(iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_K_S3(iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_K_E (iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_K_Ne(iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP) :: &
      U_Q_D (nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_Q_S1(nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_Q_S2(nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_Q_S3(nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_Q_E (nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_Q_Ne(nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP) :: &
      U_P_D (nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_P_S1(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_P_S2(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_P_S3(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_P_E (nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_P_Ne(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

    nX_G = PRODUCT( iX_E0 - iX_B0 + 1 )

    CALL TimersStart_Euler( Timer_Euler_PL_CopyIn )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( iX_B0, iX_E0, U ) &
    !$ACC CREATE( G_Q_h_d_1, G_Q_h_d_2, G_Q_h_d_3, G_Q_SqrtGm,   &
    !$ACC         G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33,      &
    !$ACC         U_K_D, U_K_S1, U_K_S2, U_K_S3, U_K_E, U_K_Ne,  &
    !$ACC         U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne,  &
    !$ACC         U_P_D, U_P_S1, U_P_S2, U_P_S3, U_P_E, U_P_Ne,  &
    !$ACC         Eps_K, Eps_P, Min_Eps_K, Min_Eps_P, Max_Eps_P, &
    !$ACC         Ye_K, Ye_P )
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_CopyIn )

    ! --- Geometry Scale Factors and Metric Determinant in Quad. Points ---

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, G, &
    !$ACC          G_Q_h_d_1, G_Q_h_d_2, G_Q_h_d_3, G_Q_SqrtGm )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        G_Q_h_d_1 (iNodeX,iX1,iX2,iX3) = G(iNodeX,iX1,iX2,iX3,iGF_h_1)
        G_Q_h_d_2 (iNodeX,iX1,iX2,iX3) = G(iNodeX,iX1,iX2,iX3,iGF_h_2)
        G_Q_h_d_3 (iNodeX,iX1,iX2,iX3) = G(iNodeX,iX1,iX2,iX3,iGF_h_3)
        G_Q_SqrtGm(iNodeX,iX1,iX2,iX3) = G(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_Permute )

    ! --- Interpolate Scale Factors to Positive Point Set (PPS) ---

#if defined( THORNADO_OACC )
    !$ACC UPDATE DEVICE( InterpMat )
#endif

    CALL TimersStart_Euler( Timer_Euler_PL_Integrate )

    CALL ComputePointValues( iX_B0, iX_E0, G_Q_h_d_1, G_P_Gm_dd_11 )
    CALL ComputePointValues( iX_B0, iX_E0, G_Q_h_d_2, G_P_Gm_dd_22 )
    CALL ComputePointValues( iX_B0, iX_E0, G_Q_h_d_3, G_P_Gm_dd_33 )

    CALL TimersStop_Euler( Timer_Euler_PL_Integrate )

    ! --- Compute Metric Components from Scale Factors in PPS ---

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, &
    !$ACC          G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33 )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iP = 1, nPT

        G_P_Gm_dd_11(iP,iX1,iX2,iX3) &
          = MAX( G_P_Gm_dd_11(iP,iX1,iX2,iX3)**2, SqrtTiny )
        G_P_Gm_dd_22(iP,iX1,iX2,iX3) &
          = MAX( G_P_Gm_dd_22(iP,iX1,iX2,iX3)**2, SqrtTiny )
        G_P_Gm_dd_33(iP,iX1,iX2,iX3) &
          = MAX( G_P_Gm_dd_33(iP,iX1,iX2,iX3)**2, SqrtTiny )

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_Permute )

    ! --- Copy Variables to Contiguous Data Blocks ---

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, U, &
    !$ACC          U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        U_Q_D (iNodeX,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_D )
        U_Q_S1(iNodeX,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_S1)
        U_Q_S2(iNodeX,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_S2)
        U_Q_S3(iNodeX,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_S3)
        U_Q_E (iNodeX,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_E )
        U_Q_Ne(iNodeX,iX1,iX2,iX3) = U(iNodeX,iX1,iX2,iX3,iCF_Ne)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_Permute )

    ! --- Compute Point Values ---

    CALL TimersStart_Euler( Timer_Euler_PL_Integrate )

    CALL ComputePointValues( iX_B0, iX_E0, U_Q_D , U_P_D  )
    CALL ComputePointValues( iX_B0, iX_E0, U_Q_S1, U_P_S1 )
    CALL ComputePointValues( iX_B0, iX_E0, U_Q_S2, U_P_S2 )
    CALL ComputePointValues( iX_B0, iX_E0, U_Q_S3, U_P_S3 )
    CALL ComputePointValues( iX_B0, iX_E0, U_Q_E , U_P_E  )
    CALL ComputePointValues( iX_B0, iX_E0, U_Q_Ne, U_P_Ne )

    CALL TimersStop_Euler( Timer_Euler_PL_Integrate )

    ! --- Compute Cell Averages ---

    CALL TimersStart_Euler( Timer_Euler_PL_Integrate )

    CALL ComputeCellAverage( iX_B0, iX_E0, G_Q_SqrtGm, U_Q_D , U_K_D  )
    CALL ComputeCellAverage( iX_B0, iX_E0, G_Q_SqrtGm, U_Q_S1, U_K_S1 )
    CALL ComputeCellAverage( iX_B0, iX_E0, G_Q_SqrtGm, U_Q_S2, U_K_S2 )
    CALL ComputeCellAverage( iX_B0, iX_E0, G_Q_SqrtGm, U_Q_S3, U_K_S3 )
    CALL ComputeCellAverage( iX_B0, iX_E0, G_Q_SqrtGm, U_Q_E , U_K_E  )
    CALL ComputeCellAverage( iX_B0, iX_E0, G_Q_SqrtGm, U_Q_Ne, U_K_Ne )

    CALL TimersStop_Euler( Timer_Euler_PL_Integrate )

    ! --- Check if Cell-Average Density is Outside Table Bounds ---

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Y_K ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K_D, U_Q_D, U_P_D, U_K_Ne, U_Q_Ne, U_P_Ne )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( Y_K )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( U_K_D(iX1,iX2,iX3) < Min_D .OR. U_K_D(iX1,iX2,iX3) > Max_D )THEN

        Y_K = BaryonMass * U_K_Ne(iX1,iX2,iX3) / U_K_D(iX1,iX2,iX3)

        ! --- Force Cell Average Inside Table Bounds ---

        IF( U_K_D(iX1,iX2,iX3) < Min_D )THEN

          U_K_D(iX1,iX2,iX3) = ( One + 1.0d-6 ) * Min_D

          DO iNodeX = 1, nDOFX
            U_Q_D(iNodeX,iX1,iX2,iX3) = U_K_D(iX1,iX2,iX3)
          END DO

        ELSE

          U_K_D(iX1,iX2,iX3) = ( One - 1.0d-6 ) * Max_D

          DO iNodeX = 1, nDOFX
            U_Q_D(iNodeX,iX1,iX2,iX3) = U_K_D(iX1,iX2,iX3)
          END DO

        END IF

        ! --- Reset Electron Density by Preserving Cell-Averaged Ye ---

        U_K_Ne(iX1,iX2,iX3) = Y_K * U_K_D(iX1,iX2,iX3) / BaryonMass

        DO iNodeX = 1, nDOFX
          U_Q_Ne(iNodeX,iX1,iX2,iX3) = U_K_Ne(iX1,iX2,iX3)
        END DO

        ! --- Recompute Point Values from Limited Solution ---

        CALL ComputePointValues_Single &
               ( U_Q_D (:,iX1,iX2,iX3), nDOFX, U_P_D (:,iX1,iX2,iX3), nPT )
        CALL ComputePointValues_Single &
               ( U_Q_Ne(:,iX1,iX2,iX3), nDOFX, U_P_Ne(:,iX1,iX2,iX3), nPT )

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

#if defined( THORNADO_OACC )
    !$ACC UPDATE HOST( U_Q_D, U_Q_Ne )
#endif

    PRINT*
    PRINT*,"1: MIN/MAX/SUM U_Q_D   = ", MINVAL( U_Q_D  ), MAXVAL( U_Q_D  ), SUM( U_Q_D  )
    PRINT*,"1: MIN/MAX/SUM U_Q_Ne  = ", MINVAL( U_Q_Ne ), MAXVAL( U_Q_Ne ), SUM( U_Q_Ne )

    ! -------------------------------------------------------------------
    ! --- Step 1 --------------------------------------------------------
    ! --- Ensure Bounded Mass Density and Positive Electron Density -----
    ! -------------------------------------------------------------------

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Min_D_K, Max_D_K, Min_N_K, Theta_1, Theta_P, D_P, N_P ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K_D, U_Q_D, U_P_D, U_K_Ne, U_Q_Ne, U_P_Ne )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( Min_D_K, Max_D_K, Min_N_K, Theta_1, Theta_P, D_P, N_P )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Min_D_K = U_P_D (1,iX1,iX2,iX3) ! --- Minimum D  in Element
      Max_D_K = U_P_D (1,iX1,iX2,iX3) ! --- Maximum D  in Element
      Min_N_K = U_P_Ne(1,iX1,iX2,iX3) ! --- Minimum Ne in Element
      DO iP = 2, nPT
        Min_D_K = MIN( Min_D_K, U_P_D (iP,iX1,iX2,iX3) )
        Max_D_K = MAX( Max_D_K, U_P_D (iP,iX1,iX2,iX3) )
        Min_N_K = MIN( Min_N_K, U_P_Ne(iP,iX1,iX2,iX3) )
      END DO

      IF( ANY( [ Min_D_K-Min_D, Max_D-Max_D_K, Min_N_K-Min_N ] < Zero ) )THEN

        Theta_1 = One

        DO iP = 1, nPT

          Theta_P = One

          D_P = U_P_D (iP,iX1,iX2,iX3)
          N_P = U_P_Ne(iP,iX1,iX2,iX3)
          DO WHILE( ANY( [ D_P-Min_D, Max_D-D_P, N_P-Min_N ] < Zero ) )

            IF( Theta_P > 1.0d-2 )THEN

              Theta_P = 0.95_DP * Theta_P

            ELSE

              Theta_P = Zero

            END IF

            D_P = (One-Theta_P) * U_K_D   (iX1,iX2,iX3) &
                      + Theta_P * U_P_D(iP,iX1,iX2,iX3)

            N_P = (One-Theta_P) * U_K_Ne   (iX1,iX2,iX3) &
                      + Theta_P * U_P_Ne(iP,iX1,iX2,iX3)

          END DO

          Theta_1 = MIN( Theta_1, SafetyFactor * Theta_P )

        END DO

        ! --- Limit Mass and Electron Density Towards Cell Average ---

        DO iNodeX = 1, nDOFX

          U_Q_D (iNodeX,iX1,iX2,iX3) &
            = Theta_1 * U_Q_D (iNodeX,iX1,iX2,iX3) &
                + (One-Theta_1) * U_K_D (iX1,iX2,iX3)

          U_Q_Ne(iNodeX,iX1,iX2,iX3) &
            = Theta_1 * U_Q_Ne(iNodeX,iX1,iX2,iX3) &
                + (One-Theta_1) * U_K_Ne(iX1,iX2,iX3)

        END DO

        ! --- Recompute Point Values from Limited Solution ---

        CALL ComputePointValues_Single &
               ( U_Q_D (1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_D (1:nPT,iX1,iX2,iX3), nPT )
        CALL ComputePointValues_Single &
               ( U_Q_Ne(1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_Ne(1:nPT,iX1,iX2,iX3), nPT )

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

#if defined( THORNADO_OACC )
    !$ACC UPDATE HOST( U_Q_D, U_Q_Ne )
#endif

    PRINT*
    PRINT*,"2: MIN/MAX/SUM U_Q_D   = ", MINVAL( U_Q_D  ), MAXVAL( U_Q_D  ), SUM( U_Q_D  )
    PRINT*,"2: MIN/MAX/SUM U_Q_Ne  = ", MINVAL( U_Q_Ne ), MAXVAL( U_Q_Ne ), SUM( U_Q_Ne )

    ! -------------------------------------------------------------------
    ! --- Step 2 --------------------------------------------------------
    ! --- Ensure Bounded Electron Fraction ------------------------------
    ! -------------------------------------------------------------------

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Min_Y_K, Max_Y_K, Y_K, Alpha, Max_D_K, Theta_2 ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K_D, U_Q_D, U_P_D, U_K_Ne, U_Q_Ne, U_P_Ne )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( Min_Y_K, MAX_Y_K, Y_K, Alpha, Max_D_K, Theta_2 )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Min_Y_K = BaryonMass * U_P_Ne(1,iX1,iX2,iX3) / U_P_D(1,iX1,iX2,iX3)
      Max_Y_K = BaryonMass * U_P_Ne(1,iX1,iX2,iX3) / U_P_D(1,iX1,iX2,iX3)
      DO iP = 2, nPT
        Min_Y_K = MIN( Min_Y_K, BaryonMass * U_P_Ne(iP,iX1,iX2,iX3) &
                                           / U_P_D (iP,iX1,iX2,iX3) )
        Max_Y_K = MAX( Max_Y_K, BaryonMass * U_P_Ne(iP,iX1,iX2,iX3) &
                                           / U_P_D (iP,iX1,iX2,iX3) )
      END DO

      IF( Min_Y_K < Min_Y .OR. Max_Y_K > Max_Y )THEN

        Y_K = BaryonMass * U_K_Ne(iX1,iX2,iX3) / U_K_D(iX1,iX2,iX3)

        Alpha = MIN( One, &
                     ABS( ( Min_Y - Y_K ) / ( Min_Y_K - Y_K ) ), &
                     ABS( ( Max_Y - Y_K ) / ( Max_Y_K - Y_K ) ) )

        Max_D_K = U_P_D(1,iX1,iX2,iX3)
        DO iP = 2, nPT
          Max_D_K = MAX( Max_D_K, U_P_D(iP,iX1,iX2,iX3) )
        END DO

        Theta_2 &
          = SafetyFactor * Alpha * U_K_D(iX1,iX2,iX3) &
            / ( Alpha * U_K_D(iX1,iX2,iX3) + (One-Alpha) * Max_D_K )

        ! --- Limit Mass and Electron Density Towards Cell Average ---

        DO iNodeX = 1, nDOFX

          U_Q_D (iNodeX,iX1,iX2,iX3) &
            = Theta_2 * U_Q_D (iNodeX,iX1,iX2,iX3) &
                + (One-Theta_2) * U_K_D (iX1,iX2,iX3)

          U_Q_Ne(iNodeX,iX1,iX2,iX3) &
            = Theta_2 * U_Q_Ne(iNodeX,iX1,iX2,iX3) &
                + (One-Theta_2) * U_K_Ne(iX1,iX2,iX3)

        END DO

        ! --- Recompute Point Values from Limited Solution ---

        CALL ComputePointValues_Single &
               ( U_Q_D (1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_D (1:nPT,iX1,iX2,iX3), nPT )
        CALL ComputePointValues_Single &
               ( U_Q_Ne(1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_Ne(1:nPT,iX1,iX2,iX3), nPT )

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

#if defined( THORNADO_OACC )
    !$ACC UPDATE HOST( U_Q_D, U_P_D, U_Q_Ne, InterpMat, &
    !$ACC              Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y, Min_N, Max_N )
#endif

    PRINT*
    PRINT*,"3: MIN/MAX/SUM U_Q_D   = ", MINVAL( U_Q_D  ), MAXVAL( U_Q_D  ), SUM( U_Q_D  )
    PRINT*,"3: MIN/MAX/SUM U_P_D   = ", MINVAL( U_P_D  ), MAXVAL( U_P_D  ), SUM( U_P_D  )
    PRINT*,"3: MIN/MAX/SUM U_Q_Ne  = ", MINVAL( U_Q_Ne ), MAXVAL( U_Q_Ne ), SUM( U_Q_Ne )
    PRINT*,"3: InterpMat = ", InterpMat
    PRINT*,"3: Min_D = ", Min_D
    PRINT*,"3: Max_D = ", Max_D
    PRINT*,"3: Min_T = ", Min_T
    PRINT*,"3: Max_T = ", Max_T
    PRINT*,"3: Min_Y = ", Min_Y
    PRINT*,"3: Max_Y = ", Max_Y
    PRINT*,"3: Min_N = ", Min_N
    PRINT*,"3: Max_N = ", Max_N

    ! -------------------------------------------------------------------
    ! --- Step 3 --------------------------------------------------------
    ! --- Ensure Bounded Specific Internal Energy -----------------------
    ! -------------------------------------------------------------------

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Eps_P, Eps_K, Ye_P, Ye_K, Min_Eps_P, Min_Eps_K, &
    !$ACC          Max_Eps_P, ITERATION, Theta_3, Theta_P ) &
    !$ACC PRESENT( iX_B0, iX_E0, &
    !$ACC          U_K_D, U_K_S1, U_K_S2, U_K_S3, U_K_E, U_K_Ne, &
    !$ACC          U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne, &
    !$ACC          U_P_D, U_P_S1, U_P_S2, U_P_S3, U_P_E, U_P_Ne, &
    !$ACC          G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33,     &
    !$ACC          Eps_K, Eps_P, Min_Eps_K, Min_Eps_P, Max_Eps_P, &
    !$ACC          Ye_K, Ye_P )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( Eps_P, Eps_K, Ye_P, Ye_K, Min_Eps_P, Min_Eps_K, &
    !$OMP          Max_Eps_P, ITERATION, Theta_3, Theta_P )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iP = 1, nPT

        CALL ComputeEpsAndYe &
               ( U_P_D       (iP,iX1,iX2,iX3), &
                 U_P_S1      (iP,iX1,iX2,iX3), &
                 U_P_S2      (iP,iX1,iX2,iX3), &
                 U_P_S3      (iP,iX1,iX2,iX3), &
                 U_P_E       (iP,iX1,iX2,iX3), &
                 U_P_Ne      (iP,iX1,iX2,iX3), &
                 G_P_Gm_dd_11(iP,iX1,iX2,iX3), &
                 G_P_Gm_dd_22(iP,iX1,iX2,iX3), &
                 G_P_Gm_dd_33(iP,iX1,iX2,iX3), &
                 Eps_P(iP), Ye_P(iP) )

        CALL ComputeSpecificInternalEnergy_TABLE &
               ( U_P_D(iP,iX1,iX2,iX3), Min_T, Ye_P(iP), Min_Eps_P(iP) )

        CALL ComputeSpecificInternalEnergy_TABLE &
               ( U_P_D(iP,iX1,iX2,iX3), Max_T, Ye_P(iP), Max_Eps_P(iP) )

      END DO

!      ITERATION = 0
      IF( ANY( Eps_P < Min_Eps_P ) )THEN
!      DO WHILE( ANY( Eps_P < Min_Eps_P ) )

!        ITERATION = ITERATION + 1

!        DO iP = 1, nPT

!          CALL ComputeEpsAndYe &
!                 ( U_K_D          (iX1,iX2,iX3), &
!                   U_K_S1         (iX1,iX2,iX3), &
!                   U_K_S2         (iX1,iX2,iX3), &
!                   U_K_S3         (iX1,iX2,iX3), &
!                   U_K_E          (iX1,iX2,iX3), &
!                   U_K_Ne         (iX1,iX2,iX3), &
!                   G_P_Gm_dd_11(iP,iX1,iX2,iX3), &
!                   G_P_Gm_dd_22(iP,iX1,iX2,iX3), &
!                   G_P_Gm_dd_33(iP,iX1,iX2,iX3), &
!                   Eps_K(iP), Ye_K(iP) )

!          CALL ComputeSpecificInternalEnergy_TABLE &
!                 ( U_K_D(iX1,iX2,iX3), Min_T, Ye_K(iP), Min_Eps_K(iP) )

!        END DO

!        IF( ITERATION .LT. 10 )THEN

!          Theta_3 = One
!          DO iP = 1, nPT
!
!            IF( Eps_P(iP) < Min_Eps_P(iP) )THEN
!
!              CALL SolveTheta_Bisection &
!                     ( U_P_D       (iP,iX1,iX2,iX3), &
!                       U_P_S1      (iP,iX1,iX2,iX3), &
!                       U_P_S2      (iP,iX1,iX2,iX3), &
!                       U_P_S3      (iP,iX1,iX2,iX3), &
!                       U_P_E       (iP,iX1,iX2,iX3), &
!                       U_K_D          (iX1,iX2,iX3), &
!                       U_K_S1         (iX1,iX2,iX3), &
!                       U_K_S2         (iX1,iX2,iX3), &
!                       U_K_S3         (iX1,iX2,iX3), &
!                       U_K_E          (iX1,iX2,iX3), &
!                       G_P_Gm_dd_11(iP,iX1,iX2,iX3), &
!                       G_P_Gm_dd_22(iP,iX1,iX2,iX3), &
!                       G_P_Gm_dd_33(iP,iX1,iX2,iX3), &
!                       Min_Eps_P(iP), Min_Eps_K(iP), &
!                       Theta_P )
!
!              Theta_3 = MIN( Theta_3, SafetyFactor * Theta_P )
!
!            END IF
!
!          END DO

          Theta_3 = Zero

!        ELSE

!          Theta_3 = Zero

!        END IF

        ! --- Limit All Conserved Fields Towards Cell Average ---

        DO iNodeX = 1, nDOFX

!          U_Q_D (iNodeX,iX1,iX2,iX3) &
!            = Theta_3 * U_Q_D (iNodeX,iX1,iX2,iX3) &
!                + (One-Theta_3) * U_K_D (iX1,iX2,iX3)

!          U_Q_S1(iNodeX,iX1,iX2,iX3) &
!            = Theta_3 * U_Q_S1(iNodeX,iX1,iX2,iX3) &
!                + (One-Theta_3) * U_K_S1(iX1,iX2,iX3)

!          U_Q_S2(iNodeX,iX1,iX2,iX3) &
!            = Theta_3 * U_Q_S2(iNodeX,iX1,iX2,iX3) &
!                + (One-Theta_3) * U_K_S2(iX1,iX2,iX3)

!          U_Q_S3(iNodeX,iX1,iX2,iX3) &
!            = Theta_3 * U_Q_S3(iNodeX,iX1,iX2,iX3) &
!                + (One-Theta_3) * U_K_S3(iX1,iX2,iX3)

!          U_Q_E (iNodeX,iX1,iX2,iX3) &
!            = Theta_3 * U_Q_E (iNodeX,iX1,iX2,iX3) &
!                + (One-Theta_3) * U_K_E (iX1,iX2,iX3)

!          U_Q_Ne(iNodeX,iX1,iX2,iX3) &
!            = Theta_3 * U_Q_Ne(iNodeX,iX1,iX2,iX3) &
!                + (One-Theta_3) * U_K_Ne(iX1,iX2,iX3)

          U_Q_D (iNodeX,iX1,iX2,iX3) = U_K_D (iX1,iX2,iX3)
          U_Q_S1(iNodeX,iX1,iX2,iX3) = U_K_S1(iX1,iX2,iX3)
          U_Q_S2(iNodeX,iX1,iX2,iX3) = U_K_S2(iX1,iX2,iX3)
          U_Q_S3(iNodeX,iX1,iX2,iX3) = U_K_S3(iX1,iX2,iX3)
          U_Q_E (iNodeX,iX1,iX2,iX3) = U_K_E (iX1,iX2,iX3)
          U_Q_Ne(iNodeX,iX1,iX2,iX3) = U_K_Ne(iX1,iX2,iX3)

        END DO

        ! --- Recompute Point Values from Limited Solution ---

!        CALL ComputePointValues_Single &
!               ( U_Q_D (1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_D (1:nPT,iX1,iX2,iX3), nPT )
!        CALL ComputePointValues_Single &
!               ( U_Q_S1(1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_S1(1:nPT,iX1,iX2,iX3), nPT )
!        CALL ComputePointValues_Single &
!               ( U_Q_S2(1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_S2(1:nPT,iX1,iX2,iX3), nPT )
!        CALL ComputePointValues_Single &
!               ( U_Q_S3(1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_S3(1:nPT,iX1,iX2,iX3), nPT )
!        CALL ComputePointValues_Single &
!               ( U_Q_E (1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_E (1:nPT,iX1,iX2,iX3), nPT )
!        CALL ComputePointValues_Single &
!               ( U_Q_Ne(1:nDOFX,iX1,iX2,iX3), nDOFX, U_P_Ne(1:nPT,iX1,iX2,iX3), nPT )

!        DO iP = 1, nPT

!          CALL ComputeEpsAndYe &
!                 ( U_P_D       (iP,iX1,iX2,iX3), &
!                   U_P_S1      (iP,iX1,iX2,iX3), &
!                   U_P_S2      (iP,iX1,iX2,iX3), &
!                   U_P_S3      (iP,iX1,iX2,iX3), &
!                   U_P_E       (iP,iX1,iX2,iX3), &
!                   U_P_Ne      (iP,iX1,iX2,iX3), &
!                   G_P_Gm_dd_11(iP,iX1,iX2,iX3), &
!                   G_P_Gm_dd_22(iP,iX1,iX2,iX3), &
!                   G_P_Gm_dd_33(iP,iX1,iX2,iX3), &
!                   Eps_P(iP), Ye_P(iP) )

!          CALL ComputeSpecificInternalEnergy_TABLE &
!                 ( U_P_D(iP,iX1,iX2,iX3), Min_T, Ye_P(iP), Min_Eps_P(iP) )

!          CALL ComputeSpecificInternalEnergy_TABLE &
!                 ( U_P_D(iP,iX1,iX2,iX3), Max_T, Ye_P(iP), Max_Eps_P(iP) )

!        END DO
       END IF
!      END DO ! --- WHILE( ANY( Eps_P < Min_Eps_P ) )

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

#if defined( THORNADO_OACC )
    !$ACC UPDATE HOST( U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne, &
    !$ACC              Eps_P, Ye_P )
#endif

    PRINT*
    PRINT*,"4: MIN/MAX/SUM U_Q_D  = ", MINVAL( U_Q_D  ), MAXVAL( U_Q_D  ), SUM( U_Q_D  )
    PRINT*,"4: MIN/MAX/SUM U_Q_S1 = ", MINVAL( U_Q_S1 ), MAXVAL( U_Q_S1 ), SUM( U_Q_S1 )
    PRINT*,"4: MIN/MAX/SUM U_Q_S2 = ", MINVAL( U_Q_S2 ), MAXVAL( U_Q_S2 ), SUM( U_Q_S2 )
    PRINT*,"4: MIN/MAX/SUM U_Q_S3 = ", MINVAL( U_Q_S3 ), MAXVAL( U_Q_S3 ), SUM( U_Q_S3 )
    PRINT*,"4: MIN/MAX/SUM U_Q_E  = ", MINVAL( U_Q_E  ), MAXVAL( U_Q_E  ), SUM( U_Q_E  )
    PRINT*,"4: MIN/MAX/SUM U_Q_Ne = ", MINVAL( U_Q_Ne ), MAXVAL( U_Q_Ne ), SUM( U_Q_Ne )

    PRINT*
    PRINT*,"4: Eps_P = ", Eps_P
    PRINT*,"4: Ye_P  = ", Ye_P

    ! --- Copy Back Fluid Variables ---

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, U, &
    !$ACC          U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        U(iNodeX,iX1,iX2,iX3,iCF_D ) = U_Q_D (iNodeX,iX1,iX2,iX3)
        U(iNodeX,iX1,iX2,iX3,iCF_S1) = U_Q_S1(iNodeX,iX1,iX2,iX3)
        U(iNodeX,iX1,iX2,iX3,iCF_S2) = U_Q_S2(iNodeX,iX1,iX2,iX3)
        U(iNodeX,iX1,iX2,iX3,iCF_S3) = U_Q_S3(iNodeX,iX1,iX2,iX3)
        U(iNodeX,iX1,iX2,iX3,iCF_E ) = U_Q_E (iNodeX,iX1,iX2,iX3)
        U(iNodeX,iX1,iX2,iX3,iCF_Ne) = U_Q_Ne(iNodeX,iX1,iX2,iX3)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_Permute )

    CALL TimersStart_Euler( Timer_Euler_PL_CopyOut )

#if   defined( THORNADO_OMP_OL )

#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U ) &
    !$ACC DELETE( iX_B0, iX_E0, &
    !$ACC         G_Q_h_d_1, G_Q_h_d_2, G_Q_h_d_3, G_Q_SqrtGm,   &
    !$ACC         G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33,      &
    !$ACC         U_K_D, U_K_S1, U_K_S2, U_K_S3, U_K_E, U_K_Ne,  &
    !$ACC         U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne,  &
    !$ACC         U_P_D, U_P_S1, U_P_S2, U_P_S3, U_P_E, U_P_Ne,  &
    !$ACC         Eps_K, Eps_P, Min_Eps_K, Min_Eps_P, Max_Eps_P, &
    !$ACC         Ye_K, Ye_P )
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_CopyOut )

!!$    DO iX3 = iX_B0(3), iX_E0(3)
!!$    DO iX2 = iX_B0(2), iX_E0(2)
!!$    DO iX1 = iX_B0(1), iX_E0(1)
!!$
!!$      G_q(1:nDOFX,1:nGF) = G(1:nDOFX,iX1,iX2,iX3,1:nGF)
!!$      U_q(1:nDOFX,1:nCF) = U(1:nDOFX,iX1,iX2,iX3,1:nCF)
!!$
!!$      NegativeStates = .FALSE.
!!$
!!$      ! --- Compute Fluid Cell Averages ---
!!$
!!$      DO iCF = 1, nCF
!!$
!!$        U_K(iCF) &
!!$          = SUM( WeightsX_q(:) * U_q(:,iCF) * G_q(:,iGF_SqrtGm) ) &
!!$              / SUM( WeightsX_q(:) * G_q(:,iGF_SqrtGm) )
!!$
!!$      END DO
!!$
!!$      Y_K = BaryonMass * U_K(iCF_Ne) / U_K(iCF_D)
!!$
!!$      IF( U_K(iCF_D) < Min_D .OR. U_K(iCF_D) > Max_D )THEN
!!$
!!$        ! --- Force Cell Averages Inside Table Bounds ---
!!$
!!$        IF( U_K(iCF_D) < Min_D )THEN
!!$
!!$          U_q(1:nDOFX,iCF_D) = ( One + 1.0d-6 ) * Min_D
!!$
!!$        ELSE
!!$
!!$          U_q(1:nDOFX,iCF_D) = ( One - 1.0d-6 ) * Max_D
!!$
!!$        END IF
!!$
!!$        ! --- Reset Electron Density by Preserving Cell-Averaged Ye ---
!!$
!!$        U_q(1:nDOFX,iCF_Ne) = Y_K * U_q(1:nDOFX,iCF_D) / BaryonMass
!!$
!!$        ! --- Recompute Cell Averages ---
!!$
!!$        U_K(iCF_D ) &
!!$          = SUM( WeightsX_q(:) * U_q(:,iCF_D ) * G_q(:,iGF_SqrtGm) ) &
!!$              / SUM( WeightsX_q(:) * G_q(:,iGF_SqrtGm) )
!!$
!!$        U_K(iCF_Ne) &
!!$          = SUM( WeightsX_q(:) * U_q(:,iCF_Ne) * G_q(:,iGF_SqrtGm) ) &
!!$              / SUM( WeightsX_q(:) * G_q(:,iGF_SqrtGm) )
!!$
!!$        NegativeStates = .TRUE.
!!$
!!$        ! --- Flag That Cell Average Violated Bounds ---
!!$
!!$        D(:,iX1,iX2,iX3,iDF_T1) = - One
!!$
!!$      END IF
!!$
!!$      CALL ComputePointValues_Fluid( U_q, U_PP )
!!$
!!$      ! -------------------------------------------------------------------
!!$      ! --- Step 1 --------------------------------------------------------
!!$      ! --- Ensure Bounded Mass Density and Positive Electron Density -----
!!$      ! -------------------------------------------------------------------
!!$
!!$      Min_D_K = MINVAL( U_PP(:,iCF_D)  ) ! --- Minimum D  in Element
!!$      Max_D_K = MAXVAL( U_PP(:,iCF_D)  ) ! --- Maximum D  in Element
!!$      Min_N_K = MINVAL( U_PP(:,iCF_Ne) ) ! --- Minimum Ne in Element
!!$
!!$      IF( ANY( [ Min_D_K - Min_D  , &
!!$                 Max_D   - Max_D_K, &
!!$                 Min_N_K - Min_N ] < Zero ) ) &
!!$      THEN
!!$
!!$        Theta_1 = One
!!$
!!$        DO iP = 1, nPT
!!$
!!$          Theta_P = One
!!$
!!$          D_P = U_PP(iP,iCF_D )
!!$          N_P = U_PP(iP,iCF_Ne)
!!$          DO WHILE( ANY( [D_P-Min_D,Max_D-D_P,N_P-Min_N] < Zero ) )
!!$
!!$            IF( Theta_P > 1.0d-2 )THEN
!!$
!!$              Theta_P = 0.95_DP * Theta_P ! --- Shrink Theta_P
!!$
!!$            ELSE
!!$
!!$              Theta_P = Zero
!!$
!!$            END IF
!!$
!!$            D_P = (One-Theta_P) * U_K(iCF_D ) + Theta_P * U_PP(iP,iCF_D )
!!$            N_P = (One-Theta_P) * U_K(iCF_Ne) + Theta_P * U_PP(iP,iCF_Ne)
!!$
!!$          END DO
!!$
!!$          Theta_1 = MIN( Theta_1, SafetyFactor * Theta_P )
!!$
!!$        END DO
!!$
!!$!        PRINT*, "iX1,iX2,iX3,Theta_1 (old) = ",iX1,iX2,iX3,Theta_1
!!$
!!$        ! --- Limit Mass and Electron Density Towards Cell Average ---
!!$
!!$        U_q(:,iCF_D ) &
!!$          = Theta_1 * U_q(:,iCF_D ) + ( One - Theta_1 ) * U_K(iCF_D )
!!$
!!$        U_q(:,iCF_Ne) &
!!$          = Theta_1 * U_q(:,iCF_Ne) + ( One - Theta_1 ) * U_K(iCF_Ne)
!!$
!!$        ! --- Recompute Point Values ---
!!$
!!$        CALL ComputePointValues_Fluid( U_q, U_PP )
!!$
!!$        NegativeStates = .TRUE.
!!$
!!$        ! --- Set Diagnostic Fields Theta 1 ---
!!$
!!$        D(:,iX1,iX2,iX3,iDF_T1) &
!!$          = MIN( MINVAL( D(:,iX1,iX2,iX3,iDF_T1) ), Theta_1 )
!!$
!!$      END IF
!!$
!!$      ! -------------------------------------------------------------------
!!$      ! --- Step 2 --------------------------------------------------------
!!$      ! --- Ensure Bounded Electron Fraction ------------------------------
!!$      ! -------------------------------------------------------------------
!!$
!!$      Min_Y_K = BaryonMass * MINVAL( U_PP(:,iCF_Ne) / U_PP(:,iCF_D) )
!!$      Max_Y_K = BaryonMass * MAXVAL( U_PP(:,iCF_Ne) / U_PP(:,iCF_D) )
!!$
!!$      IF( Min_Y_K < Min_Y .OR. Max_Y_K > Max_Y )THEN
!!$
!!$        Alpha = MIN( One, &
!!$                     ABS( ( Min_Y - Y_K ) / ( Min_Y_K - Y_K ) ), &
!!$                     ABS( ( Max_Y - Y_K ) / ( Max_Y_K - Y_K ) ) )
!!$
!!$!        PRINT*,"iX1,iX2,iX3,Alpha   (old) = ",iX1,iX2,iX3,Alpha
!!$
!!$        Theta_2 &
!!$          = SafetyFactor * Alpha * U_K(iCF_D) &
!!$            / ( Alpha * U_K(iCF_D) + (One-Alpha) * MAXVAL(U_PP(:,iCF_D)) )
!!$
!!$!        PRINT*,"iX1,iX2,iX3,Theta_2 (old) = ",iX1,iX2,iX3,Theta_2
!!$
!!$        ! --- Limit Mass and Electron Density Towards Cell Average ---
!!$
!!$        U_q(:,iCF_D ) &
!!$          = Theta_2 * U_q(:,iCF_D ) + ( One - Theta_2 ) * U_K(iCF_D )
!!$
!!$        U_q(:,iCF_Ne) &
!!$          = Theta_2 * U_q(:,iCF_Ne) + ( One - Theta_2 ) * U_K(iCF_Ne)
!!$
!!$        ! --- Recompute Point Values ---
!!$
!!$        CALL ComputePointValues_Fluid( U_q, U_PP )
!!$
!!$        NegativeStates = .TRUE.
!!$
!!$        D(:,iX1,iX2,iX3,iDF_T2) &
!!$          = MIN( MINVAL( D(:,iX1,iX2,iX3,iDF_T2) ), Theta_2 )
!!$
!!$      END IF
!!$
!!$      ! -------------------------------------------------------------------
!!$      ! --- Step 3 --------------------------------------------------------
!!$      ! --- Ensure Bounded Specific Internal Energy -----------------------
!!$      ! -------------------------------------------------------------------
!!$
!!$      DO iP = 1, nPT
!!$
!!$         CALL ComputeSpecificInternalEnergyAndElectronFraction &
!!$                ( U_PP(iP,1:nCF), E_PP(iP), Y_PP(iP) )
!!$
!!$         CALL ComputeSpecificInternalEnergy_TABLE &
!!$                ( U_PP(iP,iCF_D), Min_T, Y_PP(iP), Min_E(iP) )
!!$
!!$         CALL ComputeSpecificInternalEnergy_TABLE &
!!$                ( U_PP(iP,iCF_D), Max_T, Y_PP(iP), Max_E(iP) )
!!$
!!$      END DO
!!$
!!$      ITERATION = 0
!!$      DO WHILE( ANY( E_PP < Min_E ) )
!!$
!!$        ITERATION = ITERATION + 1
!!$
!!$        DO iP = 1, nDOFX ! --- Diagnostics (excludes points on interfaces)
!!$
!!$          D(iP,iX1,iX2,iX3,iDF_MinE) = Min_E(iP)
!!$          D(iP,iX1,iX2,iX3,iDF_MaxE) = Max_E(iP)
!!$
!!$        END DO
!!$
!!$        CALL ComputeSpecificInternalEnergyAndElectronFraction &
!!$               ( U_K(1:nCF), E_K, Y_K )
!!$
!!$        CALL ComputeSpecificInternalEnergy_TABLE &
!!$               ( U_K(iCF_D), Min_T, Y_K, Min_E_K )
!!$
!!$        IF( ITERATION .LT. 10 )THEN
!!$
!!$          Theta_3 = One
!!$          DO iP = 1, nPT
!!$
!!$            IF( E_PP(iP) < Min_E(iP) )THEN
!!$
!!$              CALL SolveTheta_Bisection_Old &
!!$                     ( U_PP(iP,1:nCF), U_K(1:nCF), Min_E(iP), Min_E_K, Theta_P )
!!$
!!$              Theta_3 = MIN( Theta_3, SafetyFactor * Theta_P )
!!$
!!$            END IF
!!$
!!$          END DO
!!$
!!$        ELSE
!!$
!!$          Theta_3 = Zero
!!$
!!$        END IF
!!$
!!$!        PRINT*,"iX1,iX2,iX3,Theta_3 (old) = ", iX1,iX2,iX3,Theta_3
!!$
!!$        ! --- Limit Towards Cell Average ---
!!$
!!$        DO iCF = 1, nCF
!!$
!!$          U_q(:,iCF) = Theta_3 * U_q(:,iCF) + ( One - Theta_3 ) * U_K(iCF)
!!$
!!$        END DO
!!$
!!$        NegativeStates = .TRUE.
!!$
!!$        D(:,iX1,iX2,iX3,iDF_T3) &
!!$          = MIN( MINVAL( D(:,iX1,iX2,iX3,iDF_T3) ), Theta_3 )
!!$
!!$        ! --- Recompute Point Values ---
!!$
!!$        CALL ComputePointValues_Fluid( U_q, U_PP )
!!$
!!$        DO iP = 1, nPT
!!$
!!$          CALL ComputeSpecificInternalEnergyAndElectronFraction &
!!$                 ( U_PP(iP,1:nCF), E_PP(iP), Y_PP(iP) )
!!$
!!$          CALL ComputeSpecificInternalEnergy_TABLE &
!!$                 ( U_PP(iP,iCF_D), Min_T, Y_PP(iP), Min_E(iP) )
!!$
!!$          CALL ComputeSpecificInternalEnergy_TABLE &
!!$                 ( U_PP(iP,iCF_D), Max_T, Y_PP(iP), Max_E(iP) )
!!$
!!$        END DO
!!$
!!$      END DO ! --- WHILE( ANY( E_PP < Min_E ) )
!!$
!!$      IF( NegativeStates )THEN
!!$
!!$        U(1:nDOFX,iX1,iX2,iX3,1:nCF) = U_q(1:nDOFX,1:nCF)
!!$
!!$      END IF
!!$
!!$    END DO
!!$    END DO
!!$    END DO

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


  SUBROUTINE ComputePointValues( iX_B0, iX_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFX, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: &
      U_P(nPT, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3))

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, nX_G, nDOFX, One, InterpMat, nPT, &
             U_Q, nDOFX, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues


  SUBROUTINE ComputePointValues_Single( U_Q, N_Q, U_P, N_P )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    INTEGER,  INTENT(in)  :: N_Q, N_P
    REAL(DP), INTENT(in)  :: U_Q(N_Q)
    REAL(DP), INTENT(out) :: U_P(N_P)

    INTEGER :: i, j

    U_P = Zero
    DO j = 1, N_Q
    DO i = 1, N_P

      U_P(i) = U_P(i) + InterpMat(i,j) * U_Q(j)

    END DO
    END DO

  END SUBROUTINE ComputePointValues_Single


  SUBROUTINE ComputeCellAverage( iX_B0, iX_E0, SqrtGm_Q, U_Q, U_K )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: &
      SqrtGm_Q &
        (nDOFX, &
         iX_B0(1):iX_E0(1), &
         iX_B0(2):iX_E0(2), &
         iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFX, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: &
      U_K(iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3))

    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: SUM1, SUM2

#if   defined( THORNADO_OMP_OL )
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
      !$OMP PRIVATE( SUM1, SUM2 )
#elif defined( THORNADO_OACC   )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRIVATE( SUM1, SUM2 ) &
      !$ACC PRESENT( U_K, WeightsX_Q, SqrtGm_Q, U_Q, iX_B0, iX_E0 )
#elif defined( THORNADO_OMP    )
      !$OMP PARALLEL DO COLLAPSE(3) &
      !$OMP PRIVATE( SUM1, SUM2 )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      SUM1 = Zero
      SUM2 = Zero
      DO iNodeX = 1, nDOFX

        SUM1 = SUM1 + WeightsX_Q(iNodeX) * SqrtGm_Q(iNodeX,iX1,iX2,iX3) &
                                         * U_Q     (iNodeX,iX1,iX2,iX3)
        SUM2 = SUM2 + WeightsX_Q(iNodeX) * SqrtGm_Q(iNodeX,iX1,iX2,iX3)

      END DO

      U_K(iX1,iX2,iX3) = SUM1 / SUM2

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeCellAverage


  SUBROUTINE ComputeEpsAndYe &
    ( D, S_d_1, S_d_2, S_d_3, E, Ne, Gm_dd_11, Gm_dd_22, Gm_dd_33, Eps, Ye )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D, S_d_1, S_d_2, S_d_3, E, Ne
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: Eps, Ye

    Eps = EpsFun( D, S_d_1, S_d_2, S_d_3, E, Gm_dd_11, Gm_dd_22, Gm_dd_33 )
    Ye  = BaryonMass * Ne / D

  END SUBROUTINE ComputeEpsAndYe


  SUBROUTINE ComputeSpecificInternalEnergyAndElectronFraction( U, E, Y )

    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: E, Y

    E = eFun( U(iCF_D), U(iCF_S1), U(iCF_S2), U(iCF_S3), U(iCF_E) )
    Y = BaryonMass * U(iCF_Ne) / U(iCF_D)

  END SUBROUTINE ComputeSpecificInternalEnergyAndElectronFraction


  FUNCTION EpsFun( D, S_d_1, S_d_2, S_d_3, E, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP)             :: EpsFun
    REAL(DP), INTENT(in) :: D, S_d_1, S_d_2, S_d_3, E
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    REAL(DP) :: S2

    S2 =   S_d_1**2 / Gm_dd_11 &
         + S_d_2**2 / Gm_dd_22 &
         + S_d_3**2 / Gm_dd_33

    EpsFun = ( E - Half * S2 / D ) / D

    RETURN
  END FUNCTION EpsFun


  PURE REAL(DP) ELEMENTAL FUNCTION eFun( D, S1, S2, S3, E )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, E

    eFun = ( E - Half * ( S1**2 + S2**2 + S3**2 ) / D ) / D

    RETURN
  END FUNCTION eFun


  SUBROUTINE SolveTheta_Bisection &
    ( D_P, S_d_1_P, S_d_2_P, S_d_3_P, E_P, &
      D_K, S_d_1_K, S_d_2_K, S_d_3_K, E_K, &
      Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P, &
      MinEps_P, MinEps_K, Theta_P )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: D_P, S_d_1_P, S_d_2_P, S_d_3_P, E_P
    REAL(DP), INTENT(in)  :: D_K, S_d_1_K, S_d_2_K, S_d_3_K, E_K
    REAL(DP), INTENT(in)  :: Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P
    REAL(DP), INTENT(in)  :: MinEps_P, MinEps_K
    REAL(DP), INTENT(out) :: Theta_P

    INTEGER,  PARAMETER :: MaxIterations = 19
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = EpsFun( D_K, S_d_1_K, S_d_2_K, S_d_3_K, E_K, &
                  Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P ) - MinEps_K

    x_b = One
    f_b = EpsFun( D_P, S_d_1_P, S_d_2_P, S_d_3_P, E_P, &
                  Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P ) - MinEps_P

    IF( .NOT. f_a * f_b < Zero )THEN
      Theta_P = Zero
      RETURN
    END IF

    dx = x_b - x_a

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED )

      dx = Half * dx
      x_c = x_a + dx

      f_c = EpsFun &
              ( x_c * D_P     + ( One - x_c ) * D_K    , &
                x_c * S_d_1_P + ( One - x_c ) * S_d_1_K, &
                x_c * S_d_2_P + ( One - x_c ) * S_d_2_K, &
                x_c * S_d_3_P + ( One - x_c ) * S_d_3_K, &
                x_c * E_P     + ( One - x_c ) * E_K    , &
                Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P ) &
            - ( x_c * MinEps_P + ( One - x_c ) * MinEps_K )

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx < dx_min ) CONVERGED = .TRUE.

      IF( ITERATION .GT. MaxIterations .AND. .NOT. CONVERGED )THEN

        x_a = Zero
        CONVERGED = .TRUE.

      END IF

    END DO

    Theta_P = x_a

  END SUBROUTINE SolveTheta_Bisection


  SUBROUTINE SolveTheta_Bisection_Old( U_Q, U_K, MinE, MinE_K, Theta_P )

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

  END SUBROUTINE SolveTheta_Bisection_Old


END MODULE Euler_PositivityLimiterModule_NonRelativistic_TABLE
