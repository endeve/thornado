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
    nDF, &
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
  REAL(DP), PARAMETER   :: Theta_1_Min  = 0.01_DP

  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: Verbose
  INTEGER               :: nX_G
  INTEGER               :: nPP(nPS)
  INTEGER               :: nPT
  REAL(DP)              :: Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y, Min_N, Max_N
  REAL(DP), ALLOCATABLE :: InterpMat(:,:)
  REAL(DP), ALLOCATABLE :: U_PP(:,:)

#if   defined( THORNADO_OMP_OL )
  !$OMP DECLARE TARGET( Min_D, Max_D, Min_T, Max_T, &
  !$OMP                 Min_Y, Max_Y, Min_N, Max_N, &
  !$OMP                 nPT )
#elif defined( THORNADO_OACC   )
  !$ACC DECLARE CREATE( Min_D, Max_D, Min_T, Max_T, &
  !$ACC                 Min_Y, Max_Y, Min_N, Max_N, &
  !$ACC                 nPT )
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
    !$OMP TARGET UPDATE TO( Min_D, Max_D, Min_T, Max_T, Min_Y, Max_Y, Min_N, Max_N )
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
    !$OMP MAP( always, to: InterpMat )
#elif defined( THORNADO_OACC   )
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
      G(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nCF)
    REAL(DP), INTENT(inout) :: &
      D(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nDF)

    LOGICAL  :: NegativeStates
    LOGICAL  :: DoStep_3
    INTEGER  :: iX1, iX2, iX3, iCF, iP, iNodeX
    INTEGER  :: ITERATION
    REAL(DP) :: Min_D_K, Max_D_K, Min_N_K
    REAL(DP) :: Min_Y_K, Max_Y_K, Min_E_K
    REAL(DP) :: Min_E(nPT), Max_E(nPT)
    REAL(DP) :: Theta_1, Theta_2, Theta_3, Theta_P, Alpha
    REAL(DP) :: D_P, N_P
    REAL(DP) :: G_q(nDOFX,nGF)
    REAL(DP) :: U_q(nDOFX,nCF), U_K(nCF)
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

    REAL(DP) :: &
      Eps_P    (nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Min_Eps_P(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Max_Eps_P(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Ye_P     (nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Eps_K    (nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Min_Eps_K(nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Ye_K     (nPT,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

    nX_G = PRODUCT( iX_E0 - iX_B0 + 1 )

    CALL TimersStart_Euler( Timer_Euler_PL_CopyIn )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, G, U ) &
    !$OMP MAP( alloc: G_Q_h_d_1, G_Q_h_d_2, G_Q_h_d_3, G_Q_SqrtGm,  &
    !$OMP             G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33,     &
    !$OMP             U_K_D, U_K_S1, U_K_S2, U_K_S3, U_K_E, U_K_Ne, &
    !$OMP             U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne, &
    !$OMP             U_P_D, U_P_S1, U_P_S2, U_P_S3, U_P_E, U_P_Ne, &
    !$OMP             Eps_P, Min_Eps_P, Max_Eps_P, Ye_P, &
    !$OMP             Eps_K, Min_Eps_K, Ye_K )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( iX_B0, iX_E0, G, U ) &
    !$ACC CREATE( G_Q_h_d_1, G_Q_h_d_2, G_Q_h_d_3, G_Q_SqrtGm,  &
    !$ACC         G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33,     &
    !$ACC         U_K_D, U_K_S1, U_K_S2, U_K_S3, U_K_E, U_K_Ne, &
    !$ACC         U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne, &
    !$ACC         U_P_D, U_P_S1, U_P_S2, U_P_S3, U_P_E, U_P_Ne, &
    !$ACC         Eps_P, Min_Eps_P, Max_Eps_P, Ye_P, &
    !$ACC         Eps_K, Min_Eps_K, Ye_K )
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_CopyIn )

    ! --- Geometry Scale Factors and Metric Determinant in Quad. Points ---

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
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

    CALL TimersStart_Euler( Timer_Euler_PL_Integrate )

    CALL ComputePointValues( iX_B0, iX_E0, G_Q_h_d_1, G_P_Gm_dd_11 )
    CALL ComputePointValues( iX_B0, iX_E0, G_Q_h_d_2, G_P_Gm_dd_22 )
    CALL ComputePointValues( iX_B0, iX_E0, G_Q_h_d_3, G_P_Gm_dd_33 )

    CALL TimersStop_Euler( Timer_Euler_PL_Integrate )

    ! --- Compute Metric Components from Scale Factors in PPS ---

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
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
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
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
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Y_K, iNodeX )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( Y_K, iNodeX ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K_D, U_Q_D, U_P_D, U_K_Ne, U_Q_Ne, U_P_Ne )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( Y_K, iNodeX )
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
               ( InterpMat, U_Q_D (:,iX1,iX2,iX3), U_P_D (:,iX1,iX2,iX3) )
        CALL ComputePointValues_Single &
               ( InterpMat, U_Q_Ne(:,iX1,iX2,iX3), U_P_Ne(:,iX1,iX2,iX3) )

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

    ! -------------------------------------------------------------------
    ! --- Step 1 --------------------------------------------------------
    ! --- Ensure Bounded Mass Density and Positive Electron Density -----
    ! -------------------------------------------------------------------

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iP, iNodeX, Min_D_K, Max_D_K, Min_N_K, Theta_1, Theta_P, D_P, N_P )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iP, iNodeX, Min_D_K, Max_D_K, Min_N_K, Theta_1, Theta_P, D_P, N_P ) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K_D, U_Q_D, U_P_D, U_K_Ne, U_Q_Ne, U_P_Ne )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iP, iNodeX, Min_D_K, Max_D_K, Min_N_K, Theta_1, Theta_P, D_P, N_P )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Min_D_K = Max_D
      Max_D_K = Min_D
      Min_N_K = Max_N
      DO iP = 1, nPT
        Min_D_K = MIN( Min_D_K, U_P_D (iP,iX1,iX2,iX3) ) ! --- Minimum D  in Element
        Max_D_K = MAX( Max_D_K, U_P_D (iP,iX1,iX2,iX3) ) ! --- Maximum D  in Element
        Min_N_K = MIN( Min_N_K, U_P_Ne(iP,iX1,iX2,iX3) ) ! --- Minimum Ne in Element
      END DO

      IF( ( Min_D_K < Min_D ) .OR. ( Max_D < Max_D_K ) .OR. ( Min_N_K < Min_N ) ) THEN

        ! --- This calculation for Theta_1 has been changed from the original backtracing
        ! --- algorithm described in 3.4.1 of Pochik et al. (2021)

        Theta_1 &
          = MIN( One, &
                 ABS( ( Min_D   - U_K_D (iX1,iX2,iX3) ) &
                    / ( Min_D_K - U_K_D (iX1,iX2,iX3) - SqrtTiny ) ), &
                 ABS( ( Max_D   - U_K_D (iX1,iX2,iX3) ) &
                    / ( Max_D_K - U_K_D (iX1,iX2,iX3) + SqrtTiny ) ), &
                 ABS( ( Min_N   - U_K_Ne(iX1,iX2,iX3) ) &
                    / ( Min_N_K - U_K_Ne(iX1,iX2,iX3) - SqrtTiny ) ) )

        Theta_1 = SafetyFactor * Theta_1
        IF ( Theta_1 < Theta_1_Min ) Theta_1 = Zero

        ! --- Limit Mass and Electron Density Towards Cell Average ---

        CALL ApplyLimiter( Theta_1, U_K_D (iX1,iX2,iX3), U_Q_D (:,iX1,iX2,iX3) )
        CALL ApplyLimiter( Theta_1, U_K_Ne(iX1,iX2,iX3), U_Q_Ne(:,iX1,iX2,iX3) )

        ! --- Recompute Point Values from Limited Solution ---

        CALL ComputePointValues_Single &
               ( InterpMat, U_Q_D (:,iX1,iX2,iX3), U_P_D (:,iX1,iX2,iX3) )
        CALL ComputePointValues_Single &
               ( InterpMat, U_Q_Ne(:,iX1,iX2,iX3), U_P_Ne(:,iX1,iX2,iX3) )

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

    ! -------------------------------------------------------------------
    ! --- Step 2 --------------------------------------------------------
    ! --- Ensure Bounded Electron Fraction ------------------------------
    ! -------------------------------------------------------------------

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_Y_K, Max_Y_K, Y_K, Alpha, Max_D_K, Theta_2 )
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

      Min_Y_K = Max_Y
      Max_Y_K = Min_Y
      DO iP = 1, nPT
        Min_Y_K = MIN( Min_Y_K, BaryonMass * U_P_Ne(iP,iX1,iX2,iX3) &
                                           / U_P_D (iP,iX1,iX2,iX3) )
        Max_Y_K = MAX( Max_Y_K, BaryonMass * U_P_Ne(iP,iX1,iX2,iX3) &
                                           / U_P_D (iP,iX1,iX2,iX3) )
      END DO

      IF( Min_Y_K < Min_Y .OR. Max_Y_K > Max_Y )THEN

        Y_K = BaryonMass * U_K_Ne(iX1,iX2,iX3) / U_K_D(iX1,iX2,iX3)

        Alpha = MIN( One, &
                     ABS( ( Min_Y - Y_K ) / ( Min_Y_K - Y_K - SqrtTiny ) ), &
                     ABS( ( Max_Y - Y_K ) / ( Max_Y_K - Y_K + SqrtTiny ) ) )

        Max_D_K = Min_D
        DO iP = 1, nPT
          Max_D_K = MAX( Max_D_K, U_P_D(iP,iX1,iX2,iX3) )
        END DO

        Theta_2 &
          = SafetyFactor * Alpha * U_K_D(iX1,iX2,iX3) &
            / ( Alpha * U_K_D(iX1,iX2,iX3) + (One-Alpha) * Max_D_K )

        ! --- Limit Mass and Electron Density Towards Cell Average ---

        CALL ApplyLimiter( Theta_2, U_K_D (iX1,iX2,iX3), U_Q_D (:,iX1,iX2,iX3) )
        CALL ApplyLimiter( Theta_2, U_K_Ne(iX1,iX2,iX3), U_Q_Ne(:,iX1,iX2,iX3) )

        ! --- Recompute Point Values from Limited Solution ---

        CALL ComputePointValues_Single &
               ( InterpMat, U_Q_D (:,iX1,iX2,iX3), U_P_D (:,iX1,iX2,iX3) )
        CALL ComputePointValues_Single &
               ( InterpMat, U_Q_Ne(:,iX1,iX2,iX3), U_P_Ne(:,iX1,iX2,iX3) )

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

    ! -------------------------------------------------------------------
    ! --- Step 3 --------------------------------------------------------
    ! --- Ensure Bounded Specific Internal Energy -----------------------
    ! -------------------------------------------------------------------

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( iP, iNodeX, ITERATION, Theta_3, Theta_P, DoStep_3 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE( iP, iNodeX, ITERATION, Theta_3, Theta_P, DoStep_3 ) &
    !$ACC PRESENT( iX_B0, iX_E0, &
    !$ACC          U_K_D, U_K_S1, U_K_S2, U_K_S3, U_K_E, U_K_Ne, &
    !$ACC          U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne, &
    !$ACC          U_P_D, U_P_S1, U_P_S2, U_P_S3, U_P_E, U_P_Ne, &
    !$ACC          G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33, &
    !$ACC          Eps_P, Min_Eps_P, Max_Eps_P, Ye_P, &
    !$ACC          Eps_K, Min_Eps_K, Ye_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(3) &
    !$OMP PRIVATE( iP, iNodeX, ITERATION, Theta_3, Theta_P, DoStep_3 )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DoStep_3 = .FALSE.
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
                 Eps_P       (iP,iX1,iX2,iX3), &
                 Ye_P        (iP,iX1,iX2,iX3) )

        CALL ComputeSpecificInternalEnergy_TABLE &
                 ( U_P_D(iP,iX1,iX2,iX3), Min_T, Ye_P(iP,iX1,iX2,iX3), &
                   Min_Eps_P(iP,iX1,iX2,iX3) )

        CALL ComputeSpecificInternalEnergy_TABLE &
                 ( U_P_D(iP,iX1,iX2,iX3), Max_T, Ye_P(iP,iX1,iX2,iX3),&
                   Max_Eps_P(iP,iX1,iX2,iX3) )

        DoStep_3 = DoStep_3 .OR. ( Eps_P(iP,iX1,iX2,iX3) < Min_Eps_P(iP,iX1,iX2,iX3) )

      END DO

      IF ( DoStep_3 ) THEN

        DO ITERATION = 1, 10

          DO iP = 1, nPT

            CALL ComputeEpsAndYe &
                   ( U_K_D          (iX1,iX2,iX3), &
                     U_K_S1         (iX1,iX2,iX3), &
                     U_K_S2         (iX1,iX2,iX3), &
                     U_K_S3         (iX1,iX2,iX3), &
                     U_K_E          (iX1,iX2,iX3), &
                     U_K_Ne         (iX1,iX2,iX3), &
                     G_P_Gm_dd_11(iP,iX1,iX2,iX3), &
                     G_P_Gm_dd_22(iP,iX1,iX2,iX3), &
                     G_P_Gm_dd_33(iP,iX1,iX2,iX3), &
                     Eps_K       (iP,iX1,iX2,iX3), &
                     Ye_K        (iP,iX1,iX2,iX3) )

            CALL ComputeSpecificInternalEnergy_TABLE &
                   ( U_K_D(iX1,iX2,iX3), Min_T, Ye_K(iP,iX1,iX2,iX3), &
                     Min_Eps_K(iP,iX1,iX2,iX3) )

          END DO

          IF( ITERATION .LT. 10 )THEN

            Theta_3 = One
            DO iP = 1, nPT

              IF( Eps_P(iP,iX1,iX2,iX3) < Min_Eps_P(iP,iX1,iX2,iX3) )THEN

                CALL SolveTheta_Bisection &
                       ( U_P_D       (iP,iX1,iX2,iX3), &
                         U_P_S1      (iP,iX1,iX2,iX3), &
                         U_P_S2      (iP,iX1,iX2,iX3), &
                         U_P_S3      (iP,iX1,iX2,iX3), &
                         U_P_E       (iP,iX1,iX2,iX3), &
                         U_K_D          (iX1,iX2,iX3), &
                         U_K_S1         (iX1,iX2,iX3), &
                         U_K_S2         (iX1,iX2,iX3), &
                         U_K_S3         (iX1,iX2,iX3), &
                         U_K_E          (iX1,iX2,iX3), &
                         G_P_Gm_dd_11(iP,iX1,iX2,iX3), &
                         G_P_Gm_dd_22(iP,iX1,iX2,iX3), &
                         G_P_Gm_dd_33(iP,iX1,iX2,iX3), &
                         Min_Eps_P   (iP,iX1,iX2,iX3), &
                         Min_Eps_K   (iP,iX1,iX2,iX3), &
                         Theta_P )

                Theta_3 = MIN( Theta_3, SafetyFactor * Theta_P )

              END IF

            END DO ! iP

          ELSE ! --- ITERATION .GE. 10

            Theta_3 = Zero

          END IF ! IF( ITERATION .LT. 10 )

          ! --- Limit All Conserved Fields Towards Cell Average ---

          CALL ApplyLimiter( Theta_3, U_K_D (iX1,iX2,iX3), U_Q_D (:,iX1,iX2,iX3) )
          CALL ApplyLimiter( Theta_3, U_K_S1(iX1,iX2,iX3), U_Q_S1(:,iX1,iX2,iX3) )
          CALL ApplyLimiter( Theta_3, U_K_S2(iX1,iX2,iX3), U_Q_S2(:,iX1,iX2,iX3) )
          CALL ApplyLimiter( Theta_3, U_K_S3(iX1,iX2,iX3), U_Q_S3(:,iX1,iX2,iX3) )
          CALL ApplyLimiter( Theta_3, U_K_E (iX1,iX2,iX3), U_Q_E (:,iX1,iX2,iX3) )
          CALL ApplyLimiter( Theta_3, U_K_Ne(iX1,iX2,iX3), U_Q_Ne(:,iX1,iX2,iX3) )

          ! --- Recompute Point Values from Limited Solution ---

          CALL ComputePointValues_Single &
                 ( InterpMat, U_Q_D (:,iX1,iX2,iX3), U_P_D (:,iX1,iX2,iX3) )
          CALL ComputePointValues_Single &
                 ( InterpMat, U_Q_S1(:,iX1,iX2,iX3), U_P_S1(:,iX1,iX2,iX3) )
          CALL ComputePointValues_Single &
                 ( InterpMat, U_Q_S2(:,iX1,iX2,iX3), U_P_S2(:,iX1,iX2,iX3) )
          CALL ComputePointValues_Single &
                 ( InterpMat, U_Q_S3(:,iX1,iX2,iX3), U_P_S3(:,iX1,iX2,iX3) )
          CALL ComputePointValues_Single &
                 ( InterpMat, U_Q_E (:,iX1,iX2,iX3), U_P_E (:,iX1,iX2,iX3) )
          CALL ComputePointValues_Single &
                 ( InterpMat, U_Q_Ne(:,iX1,iX2,iX3), U_P_Ne(:,iX1,iX2,iX3) )

          DoStep_3 = .FALSE.
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
                     Eps_P       (iP,iX1,iX2,iX3), &
                     Ye_P        (iP,iX1,iX2,iX3) )

            CALL ComputeSpecificInternalEnergy_TABLE &
                   ( U_P_D(iP,iX1,iX2,iX3), Min_T, Ye_P(iP,iX1,iX2,iX3), &
                     Min_Eps_P(iP,iX1,iX2,iX3) )

            CALL ComputeSpecificInternalEnergy_TABLE &
                   ( U_P_D(iP,iX1,iX2,iX3), Max_T, Ye_P(iP,iX1,iX2,iX3), &
                     Max_Eps_P(iP,iX1,iX2,iX3) )

            DoStep_3 = DoStep_3 .OR. ( Eps_P(iP,iX1,iX2,iX3) < Min_Eps_P(iP,iX1,iX2,iX3) )

          END DO

          IF ( .NOT. DoStep_3 ) EXIT
        END DO ! --- WHILE( DoStep_3 )

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

    ! --- Copy Back Fluid Variables ---

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
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
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: U ) &
    !$OMP MAP( release: iX_B0, iX_E0, G, &
    !$OMP               G_Q_h_d_1, G_Q_h_d_2, G_Q_h_d_3, G_Q_SqrtGm,   &
    !$OMP               G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33,      &
    !$OMP               U_K_D, U_K_S1, U_K_S2, U_K_S3, U_K_E, U_K_Ne,  &
    !$OMP               U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne,  &
    !$OMP               U_P_D, U_P_S1, U_P_S2, U_P_S3, U_P_E, U_P_Ne,  &
    !$OMP               Eps_P, Min_Eps_P, Max_Eps_P, Ye_P, &
    !$OMP               Eps_K, Min_Eps_K, Ye_K )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U ) &
    !$ACC DELETE( iX_B0, iX_E0, G, &
    !$ACC         G_Q_h_d_1, G_Q_h_d_2, G_Q_h_d_3, G_Q_SqrtGm,   &
    !$ACC         G_P_Gm_dd_11, G_P_Gm_dd_22, G_P_Gm_dd_33,      &
    !$ACC         U_K_D, U_K_S1, U_K_S2, U_K_S3, U_K_E, U_K_Ne,  &
    !$ACC         U_Q_D, U_Q_S1, U_Q_S2, U_Q_S3, U_Q_E, U_Q_Ne,  &
    !$ACC         U_P_D, U_P_S1, U_P_S2, U_P_S3, U_P_E, U_P_Ne,  &
    !$ACC         Eps_P, Min_Eps_P, Max_Eps_P, Ye_P, &
    !$ACC         Eps_K, Min_Eps_K, Ye_K )
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_CopyOut )

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

  END SUBROUTINE ApplyPositivityLimiter_Euler_NonRelativistic_TABLE


  SUBROUTINE ComputePointValues( iX_B0, iX_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFX, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: &
      U_P(nPT, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3))

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, nX_G, nDOFX, One, InterpMat, nPT, &
             U_Q, nDOFX, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues


  SUBROUTINE ComputePointValues_Single( InterpMat, U_Q, U_P )

#if   defined( THORNADO_OMP_OL )
    !$OMP DECLARE TARGET
#elif defined( THORNADO_OACC   )
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: InterpMat(nPT,nDOFX)
    REAL(DP), INTENT(in)  :: U_Q(nDOFX)
    REAL(DP), INTENT(out) :: U_P(nPT)

    REAL(DP) :: SUM1
    INTEGER  :: iNodeX, iP

    DO iP = 1, nPT
      SUM1 = Zero
      DO iNodeX = 1, nDOFX
        SUM1 = SUM1 + InterpMat(iP,iNodeX) * U_Q(iNodeX)
      END DO
      U_P(iP) = SUM1
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


  SUBROUTINE ApplyLimiter( Theta, U_K, U_Q )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: Theta, U_K
    REAL(DP), INTENT(inout) :: U_Q(nDOFX)

    INTEGER :: iNodeX

    DO iNodeX = 1, nDOFX
      U_Q(iNodeX) = Theta * U_Q(iNodeX) + ( One - Theta ) * U_K
    END DO

    RETURN
  END SUBROUTINE ApplyLimiter


END MODULE Euler_PositivityLimiterModule_NonRelativistic_TABLE
