!> Initialize, finalize, and apply positivity limiter from
!> Qin et al., (2016), JCP, 315, 323
MODULE Euler_PositivityLimiterModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE GeometryFieldsModule, ONLY: &
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
    iCF_E
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_PositivityLimiter, &
    Timer_Euler_PL_LimitCells, &
    Timer_Euler_PL_CopyIn, &
    Timer_Euler_PL_CopyOut, &
    Timer_Euler_PL_Permute, &
    Timer_Euler_PL_Integrate, &
    Timer_Euler_PL_ErrorCheck
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializePositivityLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: FinalizePositivityLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: ApplyPositivityLimiter_Euler_Relativistic_IDEAL

  LOGICAL               :: UsePositivityLimiter
  INTEGER, PARAMETER    :: nPS = 7  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total number of Positive Points
  REAL(DP)              :: Min_1, Min_2
  REAL(DP), ALLOCATABLE :: L_X(:,:)

  INTERFACE ComputePointValues
    MODULE PROCEDURE ComputePointValues_SingleField
    MODULE PROCEDURE ComputePointValues_ManyFields
  END INTERFACE ComputePointValues


CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler_Relativistic_IDEAL &
    ( UsePositivityLimiter_Option, Verbose_Option, Min_1_Option, Min_2_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option, &
                                      Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Min_2_Option

    INTEGER :: iDim, iNX, iOS
    LOGICAL :: Verbose

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

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Positivity Limiter (Euler, Relativistic, IDEAL)'
      WRITE(*,'(A)') &
        '    -----------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A6,A27,L1)') &
        '', 'UsePositivityLimiter: ', UsePositivityLimiter
      WRITE(*,*)
      WRITE(*,'(A6,A27,ES11.4E3)') &
        '', 'Min_1: ', Min_1
      WRITE(*,'(A6,A27,ES11.4E3)') &
        '', 'Min_2: ', Min_2
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

    ALLOCATE( L_X(nPT,nDOFX) )

    L_X = Zero
    DO iNX = 1, nDOFX

      L_X(iNX,iNX) = One

      IF( SUM( nPP(2:3) ) > 0 )THEN

        iOS = nPP(1)
        L_X(iOS+1:iOS+nDOFX_X1,iNX) = LX_X1_Dn(1:nDOFX_X1,iNX)

        iOS = iOS + nPP(2)
        L_X(iOS+1:iOS+nDOFX_X1,iNX) = LX_X1_Up(1:nDOFX_X1,iNX)

      END IF

      IF( SUM( nPP(4:5) ) > 0 )THEN

        iOS = SUM( nPP(1:3) )
        L_X(iOS+1:iOS+nDOFX_X2,iNX) = LX_X2_Dn(1:nDOFX_X2,iNX)

        iOS = iOS + nPP(4)
        L_X(iOS+1:iOS+nDOFX_X2,iNX) = LX_X2_Up(1:nDOFX_X2,iNX)

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        iOS = SUM( nPP(1:5) )
        L_X(iOS+1:iOS+nDOFX_X3,iNX) = LX_X3_Dn(1:nDOFX_X3,iNX)

        iOS = iOS + nPP(6)
        L_X(iOS+1:iOS+nDOFX_X3,iNX) = LX_X3_Up(1:nDOFX_X3,iNX)

      END IF

    END DO

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: L_X )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  L_X )
#endif

  END SUBROUTINE InitializePositivityLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE FinalizePositivityLimiter_Euler_Relativistic_IDEAL

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: L_X )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC EXIT DATA &
    !$ACC DELETE(       L_X )
#endif

    DEALLOCATE( L_X )

  END SUBROUTINE FinalizePositivityLimiter_Euler_Relativistic_IDEAL


  !> Iterate through the entire spatial domain and apply the positivity
  !> limiter from Qin et al., (2016), JCP, 315, 323 to each element.
  !> @param Theta_D minimum value to ensure physical density
  !> @param Theta_q minimum value to ensure physical internal energy-density
  !>        and velocity
  !> @todo Modify to accomodate GR
  SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iCF, iPT, nX_K, nCF_K
    REAL(DP) :: Min_D, Min_K, Min_ESq, Theta_D, Theta_P, q

    INTEGER :: iErr(              iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    LOGICAL :: NegativeStates(2  ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: Theta_q(          iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: SqrtGm(1:nDOFX   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    REAL(DP) :: U_Q(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: U_P(1:nPT  ,1:nCF,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(1:nCF        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    ! --- Scale factors ---

    REAL(DP) :: h1Q(1:nDOFX      ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: h2Q(1:nDOFX      ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: h3Q(1:nDOFX      ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    REAL(DP) :: h1P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: h2P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: h3P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    ! --- Metric coefficients ---

    REAL(DP) :: g1P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: g2P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: g3P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

    nX_K  = PRODUCT( iX_E0 - iX_B0 + 1 )
    nCF_K = nCF * nX_K

    CALL TimersStart_Euler( Timer_Euler_PL_CopyIn )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, G, U ) &
    !$OMP MAP( alloc: iErr, NegativeStates, Theta_q, SqrtGm, &
    !$OMP             U_Q, U_P, U_K, &
    !$OMP             h1Q, h2Q, h3Q, h1P, h2P, h3P, g1P, g2P, g3P )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, G, U ) &
    !$ACC CREATE(     iErr, NegativeStates, Theta_q, SqrtGm, &
    !$ACC             U_Q, U_P, U_K, &
    !$ACC             h1Q, h2Q, h3Q, h1P, h2P, h3P, g1P, g2P, g3P )
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( SqrtGm, h1Q, h2Q, h3Q, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      SqrtGm(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

      h1Q   (iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_h_1)
      h2Q   (iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_h_2)
      h3Q   (iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_h_3)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U_Q, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_Q(iNX,iCF,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_Permute )

    CALL TimersStart_Euler( Timer_Euler_PL_Integrate )

    CALL ComputePointValues( nCF_K, iX_B0, iX_E0, U_Q, U_P )

    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h1Q, h1P )
    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h2Q, h2P )
    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h3Q, h3P )

    CALL TimersStop_Euler( Timer_Euler_PL_Integrate )

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT(  iX_B0, iX_E0, g1P, g2P, g3P, h1P, h2P, h3P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iPT = 1, nPT

      g1P(iPT,iX1,iX2,iX3) = MAX( h1P(iPT,iX1,iX2,iX3)**2, SqrtTiny )
      g2P(iPT,iX1,iX2,iX3) = MAX( h2P(iPT,iX1,iX2,iX3)**2, SqrtTiny )
      g3P(iPT,iX1,iX2,iX3) = MAX( h3P(iPT,iX1,iX2,iX3)**2, SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( U_K, WeightsX_q, SqrtGm, U_Q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      U_K(iCF,iX1,iX2,iX3) &
        = SUM( WeightsX_q * SqrtGm(:,iX1,iX2,iX3) * U_Q(:,iCF,iX1,iX2,iX3) ) &
          / SUM( WeightsX_q * SqrtGm(:,iX1,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_Permute )

    ! --- Limit Mass-Density ---

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_D, Min_K, Theta_D )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( U_K, U_Q, U_P, NegativeStates ) &
    !$ACC PRIVATE( Min_D, Min_K, Theta_D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_D, Min_K, Theta_D )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Min_D = Min_1 * U_K(iCF_D,iX1,iX2,iX3)

      Min_K = HUGE( 1.0_DP )
      DO iPT = 1, nPT

        Min_K &
          = MIN( Min_K, U_P(iPT,iCF_D,iX1,iX2,iX3) )

      END DO

      NegativeStates(1,iX1,iX2,iX3) = .FALSE.

      IF( Min_K .LT. Min_D )THEN

        NegativeStates(1,iX1,iX2,iX3) = .TRUE.

        Theta_D &
          =   ( U_K(iCF_D,iX1,iX2,iX3) - Min_D ) &
            / ( U_K(iCF_D,iX1,iX2,iX3) - Min_K )

        DO iNX = 1, nDOFX

          U_Q(iNX,iCF_D,iX1,iX2,iX3) &
            = U_K(iCF_D,iX1,iX2,iX3) &
                + Theta_D * ( U_Q(iNX,iCF_D,iX1,iX2,iX3) &
                                - U_K(iCF_D,iX1,iX2,iX3) )

        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

    ! --- Recompute point values ---

    CALL TimersStart_Euler( Timer_Euler_PL_Integrate )

    CALL ComputePointValues( nCF_K, iX_B0, iX_E0, U_Q, U_P )

    CALL TimersStop_Euler( Timer_Euler_PL_Integrate )

    ! --- Limit q-function ---

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_ESq, Theta_P )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( U_K, U_Q, U_P, g1P, g2P, g3P, Theta_q, &
    !$ACC          iErr, NegativeStates ) &
    !$ACC PRIVATE( Min_ESq, Theta_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_ESq, Theta_P )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      iErr(iX1,iX2,iX3) = 0

      IF( U_K(iCF_E,iX1,iX2,iX3) .LT. Zero ) iErr(iX1,iX2,iX3) = 01

      Min_ESq = Min_2 * U_K(iCF_E,iX1,iX2,iX3)**2
      iErr(iX1,iX2,iX3) = 0

      Theta_q(iX1,iX2,iX3) = One

      NegativeStates(2,iX1,iX2,iX3) = .FALSE.

      DO iPT = 1, nPT

        q = Computeq &
              ( U_P(iPT,iCF_D ,iX1,iX2,iX3), &
                U_P(iPT,iCF_S1,iX1,iX2,iX3), &
                U_P(iPT,iCF_S2,iX1,iX2,iX3), &
                U_P(iPT,iCF_S3,iX1,iX2,iX3), &
                U_P(iPT,iCF_E ,iX1,iX2,iX3), &
                g1P(iPT      ,iX1,iX2,iX3), &
                g2P(iPT      ,iX1,iX2,iX3), &
                g3P(iPT      ,iX1,iX2,iX3), &
                Min_ESq )

        IF( q .LT. Zero )THEN

          NegativeStates(2,iX1,iX2,iX3) = .TRUE.

          CALL SolveTheta_Bisection &
                 ( U_P(iPT,iCF_D ,iX1,iX2,iX3), &
                   U_P(iPT,iCF_S1,iX1,iX2,iX3), &
                   U_P(iPT,iCF_S2,iX1,iX2,iX3), &
                   U_P(iPT,iCF_S3,iX1,iX2,iX3), &
                   U_P(iPT,iCF_E ,iX1,iX2,iX3), &
                   U_K(    iCF_D ,iX1,iX2,iX3), &
                   U_K(    iCF_S1,iX1,iX2,iX3), &
                   U_K(    iCF_S2,iX1,iX2,iX3), &
                   U_K(    iCF_S3,iX1,iX2,iX3), &
                   U_K(    iCF_E ,iX1,iX2,iX3), &
                   g1P(iPT      ,iX1,iX2,iX3), &
                   g2P(iPT      ,iX1,iX2,iX3), &
                   g3P(iPT      ,iX1,iX2,iX3), &
                   Min_ESq, Theta_P, &
                   iErr(iX1,iX2,iX3) )

          Theta_q(iX1,iX2,iX3) = MIN( Theta_q(iX1,iX2,iX3), Theta_P )

        END IF

      END DO

    END DO
    END DO
    END DO

    ! --- Limit all variables towards cell-average ---

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( U, U_K, U_Q, Theta_q, NegativeStates )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( NegativeStates(2,iX1,iX2,iX3) )THEN

        DO iCF = 1, nCF
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX3,iCF) &
            = U_K(iCF,iX1,iX2,iX3) &
                + Theta_q(iX1,iX2,iX3) * ( U_Q(iNX,iCF,iX1,iX2,iX3) &
                                             - U_K(iCF,iX1,iX2,iX3) )

        END DO
        END DO

      ELSE

        DO iCF = 1, nCF
        DO iNX = 1, nDOFX

          U(iNX,iX1,iX2,iX3,iCF) = U_Q(iNX,iCF,iX1,iX2,iX3)

        END DO
        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

    CALL TimersStart_Euler( Timer_Euler_PL_CopyOut )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U, iErr ) &
    !$OMP MAP( release: NegativeStates, Theta_q, iX_B0, iX_E0, SqrtGm, &
    !$OMP               G, U_Q, U_P, U_K, &
    !$OMP               h1Q, h2Q, h3Q, h1P, h2P, h3P, g1P, g2P, g3P )
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U, iErr ) &
    !$ACC DELETE(       NegativeStates, Theta_q, iX_B0, iX_E0, SqrtGm, &
    !$ACC               G, U_Q, U_P, U_K, &
    !$ACC               h1Q, h2Q, h3Q, h1P, h2P, h3P, g1P, g2P, g3P )
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_CopyOut )

    CALL TimersStart_Euler( Timer_Euler_PL_ErrorCheck )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( iErr(iX1,iX2,iX3) .NE. 0 ) &
        CALL DescribeError_Euler &
          ( iErr(iX1,iX2,iX3), &
            Int_Option  = [ iX1, iX2, iX3 ], &
            Real_Option = [ U(iNX,iX1,iX2,iX3,iCF_D ), &
                            U(iNX,iX1,iX2,iX3,iCF_S1), &
                            U(iNX,iX1,iX2,iX3,iCF_S2), &
                            U(iNX,iX1,iX2,iX3,iCF_S3), &
                            U(iNX,iX1,iX2,iX3,iCF_E ), &
                            G(iNX,iX1,iX2,iX3,iGF_h_1), &
                            G(iNX,iX1,iX2,iX3,iGF_h_2), &
                            G(iNX,iX1,iX2,iX3,iGF_h_3) ] )


    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_ErrorCheck )

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

  END SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ComputePointValues_SingleField( N, iX_B0, iX_E0, h_Q, h_P )

    INTEGER,  INTENT(in)  :: N, iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: h_Q(1:nDOFX,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3) )
    REAL(DP), INTENT(out) :: h_P(1:nPT  ,iX_B0(1):iX_E0(1), &
                                         iX_B0(2):iX_E0(2), &
                                         iX_B0(3):iX_E0(3) )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, N, nDOFX, One, L_X, nPT, &
             h_Q, nDOFX, Zero, h_P, nPT )

  END SUBROUTINE ComputePointValues_SingleField


  SUBROUTINE ComputePointValues_ManyFields( N, iX_B0, iX_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: N, iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: U_Q(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                               iX_B0(2):iX_E0(2), &
                                               iX_B0(3):iX_E0(3) )
    REAL(DP), INTENT(out) :: U_P(1:nPT  ,1:nCF,iX_B0(1):iX_E0(1), &
                                               iX_B0(2):iX_E0(2), &
                                               iX_B0(3):iX_E0(3) )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, N, nDOFX, One, L_X, nPT, &
             U_Q, nDOFX, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues_ManyFields


  FUNCTION Computeq &
    ( D, S1, S2, S3, tau, g1, g2, g3, Min_ESq ) RESULT( q )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, S1, S2, S3, tau, g1, g2, g3, Min_ESq

    REAL(DP) :: q

    q = tau + D &
          - SQRT( D**2 + ( S1**2 / g1 + S2**2 / g2 + S3**2 / g3 ) &
                  + Min_ESq )

    RETURN
  END FUNCTION Computeq


  SUBROUTINE SolveTheta_Bisection &
    ( D_P, S1_P, S2_P, S3_P, E_P, &
      D_K, S1_K, S2_K, S3_K, E_K, &
      g1_P, g2_P, g3_P, Min_ESq, Theta_P, iErr )

#if defined(THORNADO_OMP_OL) && !defined(THORNADO_EULER_NOGPU)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC) && !defined(THORNADO_EULER_NOGPU)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: D_P, S1_P, S2_P, S3_P, E_P, &
                               D_K, S1_K, S2_K, S3_K, E_K, &
                               g1_P, g2_P, g3_P, Min_ESq
    REAL(DP), INTENT(out)   :: Theta_P
    INTEGER , INTENT(inout) :: iErr

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0e-3_DP

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = Computeq &
            ( x_a *  D_P + ( One - x_a ) *  D_K, &
              x_a * S1_P + ( One - x_a ) * S1_K, &
              x_a * S2_P + ( One - x_a ) * S2_K, &
              x_a * S3_P + ( One - x_a ) * S3_K, &
              x_a *  E_P + ( One - x_a ) *  E_K, &
              g1_P, g2_P, g3_P, Min_ESq )

    x_b = One
    f_b = Computeq &
            ( x_b *  D_P + ( One - x_b ) *  D_K, &
              x_b * S1_P + ( One - x_b ) * S1_K, &
              x_b * S2_P + ( One - x_b ) * S2_K, &
              x_b * S3_P + ( One - x_b ) * S3_K, &
              x_b *  E_P + ( One - x_b ) *  E_K, &
              g1_P, g2_P, g3_P, Min_ESq )

    IF( .NOT. f_a * f_b < 0 ) iErr = 02

    dx = x_b - x_a

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx
      f_c = Computeq &
              ( x_c *  D_P + ( One - x_c ) *  D_K, &
                x_c * S1_P + ( One - x_c ) * S1_K, &
                x_c * S2_P + ( One - x_c ) * S2_K, &
                x_c * S3_P + ( One - x_c ) * S3_K, &
                x_c *  E_P + ( One - x_c ) *  E_K, &
                g1_P, g2_P, g3_P, Min_ESq )

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx .LT. dx_min ) CONVERGED = .TRUE.

      IF( ITERATION .GT. MAX_IT .AND. .NOT. CONVERGED )THEN

        IF( ITERATION .GT. MAX_IT + 3 ) iErr = 03

      END IF

    END DO

    Theta_P = x_a

  END SUBROUTINE SolveTheta_Bisection


END MODULE Euler_PositivityLimiterModule_Relativistic_IDEAL
