#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_POSITIVITYLIMITER
#endif
MODULE TwoMoment_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, nDOF, nDOFE, nDOFX
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_PositivityLimiter, &
    Timer_PL_In, &
    Timer_PL_Points, &
    Timer_PL_CellAverage, &
    Timer_PL_Theta_1, &
    Timer_PL_Theta_2, &
    Timer_PL_Out
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply, &
    MatrixVectorMultiply
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_X1_Dn, L_X1_Up, &
    L_X2_Dn, L_X2_Up, &
    L_X3_Dn, L_X3_Up
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE RadiationFieldsModule, ONLY: &
    nSpecies, nCR, &
    iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_TwoMoment
  PUBLIC :: FinalizePositivityLimiter_TwoMoment
  PUBLIC :: ApplyPositivityLimiter_TwoMoment
  PUBLIC :: TallyPositivityLimiter_TwoMoment

  CHARACTER(256)        :: TallyFileName
  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: UsePositivityLimiterTally
  INTEGER,    PARAMETER :: nPS = 9  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Number of Positive Points
  INTEGER               :: nR, nR_1, nZ(4)
  REAL(DP)              :: Min_1, Max_1, Min_2
  REAL(DP)              :: Theta_FD, MinTheta_1, MinTheta_2
  REAL(DP),   PARAMETER :: Theta_Eps = 1.0_DP - 10.0_DP*EPSILON(1.0_DP)
  REAL(DP), ALLOCATABLE :: L_X(:,:)

  INTEGER,    PARAMETER :: MAX_IT = 19

  INTERFACE GammaFun
    MODULE PROCEDURE GammaFun_Scalar
    MODULE PROCEDURE GammaFun_Vector
  END INTERFACE

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET( Min_1, Max_1, Min_2, Theta_FD )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE( Min_1, Max_1, Min_2, Theta_FD )
#endif

CONTAINS


  SUBROUTINE InitializePositivityLimiter_TwoMoment &
    ( Min_1_Option, Max_1_Option, Min_2_Option, &
      UsePositivityLimiter_Option, &
      UsePositivityLimiterTally_Option, Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Max_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiterTally_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose
    INTEGER :: i, FileUnit
    INTEGER :: iNode, iNodeX1, iNodeX2, iNodeX3, iOS

    IF( PRESENT( Min_1_Option ) )THEN
      Min_1 = Min_1_Option
    ELSE
      Min_1 = - HUGE( One )
    END IF

    IF( PRESENT( Max_1_Option ) )THEN
      Max_1 = Max_1_Option
    ELSE
      Max_1 = + HUGE( One )
    END IF

    IF( PRESENT( Min_2_Option ) )THEN
      Min_2 = Min_2_Option
    ELSE
      Min_2 = - HUGE( One )
    END IF

    IF( PRESENT( UsePositivityLimiter_Option ) )THEN
      UsePositivityLimiter = UsePositivityLimiter_Option
    ELSE
      UsePositivityLimiter = .TRUE.
    END IF

    IF( PRESENT( UsePositivityLimiterTally_Option ) )THEN
      UsePositivityLimiterTally = UsePositivityLimiterTally_Option
    ELSE
      UsePositivityLimiterTally = .FALSE.
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    Theta_FD = One
    IF( Max_1 > One ) &
      Theta_FD = Zero

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'InitializePositivityLimiter'
      WRITE(*,*)
      WRITE(*,'(A6,A,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter 
      WRITE(*,*)
      WRITE(*,'(A6,A12,ES23.15E3)') '', 'Min_1 = ', Min_1
      WRITE(*,'(A6,A12,ES23.15E3)') '', 'Max_1 = ', Max_1
      WRITE(*,'(A6,A12,ES23.15E3)') '', 'Min_2 = ', Min_2
      WRITE(*,'(A6,A12,ES23.15E3)') '', 'Theta_FD = ', Theta_FD

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO( Min_1, Max_1, Min_2, Theta_FD )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE( Min_1, Max_1, Min_2, Theta_FD )
#endif

    nPP = 0
    nPP(1) = PRODUCT( nNodesZ )

    DO i = 2, 4 ! --- Exclude energy dimension for now ---

      IF( nNodesZ(i) > 1 )THEN

        nPP(2*i:2*i+1) &
          = PRODUCT( nNodesZ, MASK = [1,2,3,4] .NE. i )

      END IF

    END DO

    nPT = SUM( nPP )

    ALLOCATE( L_X(nPT,nDOF) )
    L_X = Zero
    DO iNode = 1, nDOF

      L_X(iNode,iNode) = One

      IF ( SUM( nPP(4:5) ) > 0 ) THEN

        iOS = SUM( nPP(1:3) )
        L_X(iOS+1:iOS+nDOF_X1,iNode) = L_X1_Dn(1:nDOF_X1,iNode)

        iOS = iOS + nPP(4)
        L_X(iOS+1:iOS+nDOF_X1,iNode) = L_X1_Up(1:nDOF_X1,iNode)

      END IF

      IF ( SUM( nPP(6:7) ) > 0 ) THEN

        iOS = SUM( nPP(1:5) )
        L_X(iOS+1:iOS+nDOF_X2,iNode) = L_X2_Dn(1:nDOF_X2,iNode)

        iOS = iOS + nPP(6)
        L_X(iOS+1:iOS+nDOF_X2,iNode) = L_X2_Up(1:nDOF_X2,iNode)

      END IF

      IF ( SUM( nPP(8:9) ) > 0 ) THEN

        iOS = SUM( nPP(1:7) )
        L_X(iOS+1:iOS+nDOF_X3,iNode) = L_X3_Dn(1:nDOF_X3,iNode)

        iOS = iOS + nPP(8)
        L_X(iOS+1:iOS+nDOF_X3,iNode) = L_X3_Up(1:nDOF_X3,iNode)

      END IF

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: L_X )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( L_X )
#endif

    ! --- For Tally of Positivity Limiter ---

    IF( UsePositivityLimiterTally )THEN

      TallyFileName = '../Output/PositivityLimiterTally.dat'

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ) )

      WRITE( FileUnit, '(3(A14,x))' ) 'Time', 'MinTheta_1', 'MinTheta_2'

      CLOSE( FileUnit )

    END IF

  END SUBROUTINE InitializePositivityLimiter_TwoMoment


  SUBROUTINE FinalizePositivityLimiter_TwoMoment

    IF ( ALLOCATED( L_X ) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: L_X )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA &
      !$ACC DELETE( L_X )
#endif
      DEALLOCATE( L_X )

    END IF

  END SUBROUTINE FinalizePositivityLimiter_TwoMoment


  SUBROUTINE ApplyPositivityLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iNode, iZ1, iZ2, iZ3, iZ4, iS, iCR, iP
    INTEGER  :: nNeg_1, nNeg_2
    REAL(DP) :: Min_K, Max_K, Theta_1, Min_Gam, Theta_2, Gam, Theta_P

    REAL(DP) :: Min_K_S  (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: Max_K_S  (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: Theta_1_S(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: Min_Gam_S(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: Theta_2_S(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    INTEGER  :: iError(nPT,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    LOGICAL  :: NegativeStates(2,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    REAL(DP) :: U_Q_N  (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_Q_G1 (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_Q_G2 (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_Q_G3 (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    REAL(DP) :: U_P_N  (nPT ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_P_G1 (nPT ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_P_G2 (nPT ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_P_G3 (nPT ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    REAL(DP) :: U_K_N  (     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_K_G1 (     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_K_G2 (     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_K_G3 (     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart( Timer_PositivityLimiter )

    nZ = iZ_E0 - iZ_B0 + 1
    nR_1 = nSpecies * PRODUCT( nZ )
    nR = nR_1 * nCR

    NegativeStates = .FALSE.
    iError = 0

    MinTheta_1 = One
    MinTheta_2 = One

    nNeg_1 = 0
    nNeg_2 = 0

    CALL TimersStart( Timer_PL_In )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, NegativeStates, &
    !$OMP          U, iError, MinTheta_1, MinTheta_2 ) &
    !$OMP MAP( alloc: U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, &
    !$OMP             U_P_N, U_P_G1, U_P_G2, U_P_G3, &
    !$OMP             U_K_N, U_K_G1, U_K_G2, U_K_G3, &
    !$OMP             Min_K_S, Max_K_S, Theta_1_S, Min_Gam_S, Theta_2_S )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( iZ_B0, iZ_E0, NegativeStates, &
    !$ACC         U, iError, MinTheta_1, MinTheta_2 ) &
    !$ACC CREATE( U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, &
    !$ACC         U_P_N, U_P_G1, U_P_G2, U_P_G3, &
    !$ACC         U_K_N, U_K_G1, U_K_G2, U_K_G3, &
    !$ACC         Min_K_S, Max_K_S, Theta_1_S, Min_Gam_S, Theta_2_S )
#endif

    CALL TimersStop( Timer_PL_In )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, U, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              U_Q_N (:,iZ1,iZ2,iZ3,iZ4,iS) = U(:,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)
              U_Q_G1(:,iZ1,iZ2,iZ3,iZ4,iS) = U(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)
              U_Q_G2(:,iZ1,iZ2,iZ3,iZ4,iS) = U(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)
              U_Q_G3(:,iZ1,iZ2,iZ3,iZ4,iS) = U(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)

            END DO
          END DO
        END DO
      END DO
    END DO

    ! --- Point Values ---

    CALL TimersStart( Timer_PL_Points )

    CALL ComputePointValues( iZ_B0, iZ_E0, U_Q_N , U_P_N  )
    CALL ComputePointValues( iZ_B0, iZ_E0, U_Q_G1, U_P_G1 )
    CALL ComputePointValues( iZ_B0, iZ_E0, U_Q_G2, U_P_G2 )
    CALL ComputePointValues( iZ_B0, iZ_E0, U_Q_G3, U_P_G3 )

    CALL TimersStop( Timer_PL_Points )

    ! --- Cell Average ---

    CALL TimersStart( Timer_PL_CellAverage )

    CALL MatrixVectorMultiply &
      ( 'T', nDOF, nR_1, One, U_Q_N , nDOF, Weights_q, 1, Zero, U_K_N , 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nDOF, nR_1, One, U_Q_G1, nDOF, Weights_q, 1, Zero, U_K_G1, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nDOF, nR_1, One, U_Q_G2, nDOF, Weights_q, 1, Zero, U_K_G2, 1 )
    CALL MatrixVectorMultiply &
      ( 'T', nDOF, nR_1, One, U_Q_G3, nDOF, Weights_q, 1, Zero, U_K_G3, 1 )

    CALL TimersStop( Timer_PL_CellAverage )

    ! --- Ensure Bounded Density ---

    CALL TimersStart( Timer_PL_Theta_1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( Min_K, Max_K, Theta_1 ) &
    !$OMP REDUCTION( min: MinTheta_1 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( Min_K, Max_K, Theta_1 ) &
    !$ACC REDUCTION( min: MinTheta_1 ) &
    !$ACC PRESENT( U_P_N, U_K_N, U_Q_N, NegativeStates, Min_K_S, Max_K_S, &
    !$ACC          Theta_1_S, MinTheta_1, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( Min_K, Max_K, Theta_1 ) &
    !$OMP REDUCTION( min: MinTheta_1 )
#endif
    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              Min_K = One
              Max_K = Zero
              DO iP = 1, nPT
                Min_K = MIN( Min_K, U_P_N(iP,iZ1,iZ2,iZ3,iZ4,iS) )
                Max_K = MAX( Max_K, U_P_N(iP,iZ1,iZ2,iZ3,iZ4,iS) )
              END DO
              !Min_K = MINVAL( U_P_N(:,iZ1,iZ2,iZ3,iZ4,iS) )
              !Max_K = MAXVAL( U_P_N(:,iZ1,iZ2,iZ3,iZ4,iS) )

              IF ( Min_K < Min_1 .OR. Max_K > Max_1 ) THEN

                NegativeStates(1,iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

                ! --- Limit Density Towards Cell Average ---

                !U_K_N(iZ1,iZ2,iZ3,iZ4,iS) = DOT_PRODUCT( Weights_q, U_Q_N(:,iZ1,iZ2,iZ3,iZ4,iS) )

                Theta_1 &
                  = Theta_Eps * MIN( One, &
                    ABS( (Min_1-U_K_N(iZ1,iZ2,iZ3,iZ4,iS)) &
                          / (Min_K-U_K_N(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny) ), &
                    ABS( (Max_1-U_K_N(iZ1,iZ2,iZ3,iZ4,iS)) &
                          / (Max_K-U_K_N(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny) ) )

                U_Q_N(:,iZ1,iZ2,iZ3,iZ4,iS) &
                  =   U_Q_N(:,iZ1,iZ2,iZ3,iZ4,iS) * Theta_1 &
                    + U_K_N(iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_1 )

                Min_K_S(iZ1,iZ2,iZ3,iZ4,iS) = Min_K
                Max_K_S(iZ1,iZ2,iZ3,iZ4,iS) = Max_K
                Theta_1_S(iZ1,iZ2,iZ3,iZ4,iS) = Theta_1
                MinTheta_1 = MIN( Theta_1, MinTheta_1 )

              END IF

            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_PL_Theta_1 )

    ! --- Recompute Point Values ---

    CALL TimersStart( Timer_PL_Points )

    CALL ComputePointValues( iZ_B0, iZ_E0, U_Q_N , U_P_N  )

    CALL TimersStop( Timer_PL_Points )

#ifdef THORNADO_DEBUG_POSITIVITYLIMITER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( NegativeStates )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( NegativeStates )
#endif
    IF ( ANY( NegativeStates ) ) THEN
      DO iS = 1, nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)

                IF ( NegativeStates(1,iZ1,iZ2,iZ3,iZ4,iS) ) THEN
                  nNeg_1 = nNeg_1 + 1
                  IF ( nNeg_1 <= 5 ) THEN
#if defined(THORNADO_OMP_OL)
                    !$OMP TARGET UPDATE FROM( &
                    !$OMP   Min_K_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   Max_K_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   Theta_1_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_K_N(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_Q_N(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS) )
#elif defined(THORNADO_OACC)
                    !$ACC UPDATE HOST( &
                    !$ACC   Min_K_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   Max_K_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   Theta_1_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_K_N(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_Q_N(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS) )
#endif
                    Min_K = Min_K_S(iZ1,iZ2,iZ3,iZ4,iS)
                    Max_K = Max_K_S(iZ1,iZ2,iZ3,iZ4,iS)
                    Theta_1 = Theta_1_S(iZ1,iZ2,iZ3,iZ4,iS)
                    WRITE(*,'(a30,5i4)') &
                      'Neg. UQ(N) @ ', iZ1, iZ2, iZ3, iZ4, iS
                    WRITE(*,'(a30,3es23.15)') &
                      'Min_K, Max_K, Theta_1 : ', &
                      Min_K, Max_K, Theta_1
                    WRITE(*,'(a30,es23.15,i4)') &
                      '        MINVAL(UQ(N)) : ', &
                      MINVAL(U_Q_N(:,iZ1,iZ2,iZ3,iZ4,iS)), &
                      MINLOC(U_Q_N(:,iZ1,iZ2,iZ3,iZ4,iS))
                    WRITE(*,'(a30,es23.15)') &
                      '                  U_K : ', &
                      U_K_N(iZ1,iZ2,iZ3,iZ4,iS)
                  END IF
                END IF

              END DO
            END DO
          END DO
        END DO
      END DO
    END IF
#endif

    ! --- Ensure Positive Gamma ---

    CALL TimersStart( Timer_PL_Theta_2 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( Gam, Min_Gam, Theta_2, Theta_P ) &
    !$OMP REDUCTION( min: MinTheta_2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( Gam, Min_Gam, Theta_2, Theta_P ) &
    !$ACC REDUCTION( min: MinTheta_2 ) &
    !$ACC PRESENT( U_P_N, U_P_G1, U_P_G2, U_P_G3, &
    !$ACC          U_K_N, U_K_G1, U_K_G2, U_K_G3, &
    !$ACC          U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, &
    !$ACC          NegativeStates, Min_Gam_S, Theta_2_S, &
    !$ACC          MinTheta_2, iError, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( Gam, Min_Gam, Theta_2, Theta_P ) &
    !$OMP REDUCTION( min: MinTheta_2 )
#endif
    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              Min_Gam = One
              Theta_2 = One

              DO iP = 1, nPT

                Gam = GammaFun &
                      ( U_P_N (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                        U_P_G1(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                        U_P_G2(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                        U_P_G3(iP,iZ1,iZ2,iZ3,iZ4,iS) )

                Min_Gam = MIN( Min_Gam, Gam )

                IF ( Gam < Min_2 ) THEN

                  CALL SolveTheta_Bisection &
                        ( U_P_N  (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                          U_P_G1 (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                          U_P_G2 (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                          U_P_G3 (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                          U_K_N  (   iZ1,iZ2,iZ3,iZ4,iS), &
                          U_K_G1 (   iZ1,iZ2,iZ3,iZ4,iS), &
                          U_K_G2 (   iZ1,iZ2,iZ3,iZ4,iS), &
                          U_K_G3 (   iZ1,iZ2,iZ3,iZ4,iS), &
                          Theta_P, &
                          iError (iP,iZ1,iZ2,iZ3,iZ4,iS) )

                  Theta_2 = MIN( Theta_2, Theta_P )

                END IF

              END DO

              !Min_Gam = MINVAL( Gam(:) )

              IF ( Min_Gam < Min_2 ) THEN

                NegativeStates(2,iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

                !U_K_N (iZ1,iZ2,iZ3,iZ4,iS) = DOT_PRODUCT( Weights_q, U_Q_N (:,iZ1,iZ2,iZ3,iZ4,iS) )
                !U_K_G1(iZ1,iZ2,iZ3,iZ4,iS) = DOT_PRODUCT( Weights_q, U_Q_G1(:,iZ1,iZ2,iZ3,iZ4,iS) )
                !U_K_G2(iZ1,iZ2,iZ3,iZ4,iS) = DOT_PRODUCT( Weights_q, U_Q_G2(:,iZ1,iZ2,iZ3,iZ4,iS) )
                !U_K_G3(iZ1,iZ2,iZ3,iZ4,iS) = DOT_PRODUCT( Weights_q, U_Q_G3(:,iZ1,iZ2,iZ3,iZ4,iS) )

                ! --- Limit Towards Cell Average ---

                !Theta_2 = Theta_Eps * MINVAL( Theta_P(:) )
                Theta_2 = Theta_Eps * Theta_2

                U_Q_N (:,iZ1,iZ2,iZ3,iZ4,iS) &
                  =   U_Q_N (:,iZ1,iZ2,iZ3,iZ4,iS) * Theta_2 &
                    + U_K_N (iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_2 )

                U_Q_G1(:,iZ1,iZ2,iZ3,iZ4,iS) &
                  =   U_Q_G1(:,iZ1,iZ2,iZ3,iZ4,iS) * Theta_2 &
                    + U_K_G1(iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_2 )

                U_Q_G2(:,iZ1,iZ2,iZ3,iZ4,iS) &
                  =   U_Q_G2(:,iZ1,iZ2,iZ3,iZ4,iS) * Theta_2 &
                    + U_K_G2(iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_2 )

                U_Q_G3(:,iZ1,iZ2,iZ3,iZ4,iS) &
                  =   U_Q_G3(:,iZ1,iZ2,iZ3,iZ4,iS) * Theta_2 &
                    + U_K_G3(iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_2 )

                Min_Gam_S(iZ1,iZ2,iZ3,iZ4,iS) = Min_Gam
                Theta_2_S(iZ1,iZ2,iZ3,iZ4,iS) = Theta_2
                MinTheta_2 = MIN( Theta_2, MinTheta_2 )

              END IF

            END DO
          END DO
        END DO
      END DO
    END DO

    CALL TimersStop( Timer_PL_Theta_2 )

#ifdef THORNADO_DEBUG_POSITIVITYLIMITER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( NegativeStates )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( NegativeStates )
#endif
    IF ( ANY( NegativeStates ) ) THEN
      DO iS = 1, nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)

                IF ( NegativeStates(2,iZ1,iZ2,iZ3,iZ4,iS) ) THEN
                  nNeg_2 = nNeg_2 + 1
                  IF ( nNeg_2 <= 5 ) THEN
#if defined(THORNADO_OMP_OL)
                    !$OMP TARGET UPDATE FROM( &
                    !$OMP   Min_Gam_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   Theta_2_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_K_N (iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_K_G1(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_K_G2(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_K_G3(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_Q_N (1:nDOF,iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_Q_G1(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_Q_G2(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS), &
                    !$OMP   U_Q_G3(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS) )
#elif defined(THORNADO_OACC)
                    !$ACC UPDATE HOST( &
                    !$ACC   Min_Gam_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   Theta_2_S(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_K_N(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_K_G1(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_K_G2(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_K_G3(iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_Q_N(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_Q_G1(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_Q_G2(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS), &
                    !$ACC   U_Q_G3(1:nDOF,iZ1,iZ2,iZ3,iZ4,iS) )
#endif
                    Min_Gam = Min_Gam_S(iZ1,iZ2,iZ3,iZ4,iS)
                    Theta_2 = Theta_2_S(iZ1,iZ2,iZ3,iZ4,iS)
                    WRITE(*,'(a30,5i4)') &
                      'Neg. Gamma @ ', iZ1, iZ2, iZ3, iZ4, iS
                    WRITE(*,'(a30,2es23.15)') &
                      '   Min_Gamma, Theta_2 : ', Min_Gam, Theta_2
                    WRITE(*,'(a30,4es23.15)') &
                      '        MINVAL(UQ(:)) : ', &
                      MINVAL(U_Q_N (:,iZ1,iZ2,iZ3,iZ4,iS)), &
                      MINVAL(U_Q_G1(:,iZ1,iZ2,iZ3,iZ4,iS)), &
                      MINVAL(U_Q_G2(:,iZ1,iZ2,iZ3,iZ4,iS)), &
                      MINVAL(U_Q_G3(:,iZ1,iZ2,iZ3,iZ4,iS))
                    WRITE(*,'(a30,4es23.15)') &
                      '                  U_K : ', &
                      U_K_N (iZ1,iZ2,iZ3,iZ4,iS), &
                      U_K_G1(iZ1,iZ2,iZ3,iZ4,iS), &
                      U_K_G2(iZ1,iZ2,iZ3,iZ4,iS), &
                      U_K_G3(iZ1,iZ2,iZ3,iZ4,iS)
                  END IF
                END IF

              END DO
            END DO
          END DO
        END DO
      END DO
    END IF
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U, U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, NegativeStates, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              IF ( NegativeStates(2,iZ1,iZ2,iZ3,iZ4,iS) ) THEN

                U(:,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) = U_Q_N (:,iZ1,iZ2,iZ3,iZ4,iS)
                U(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) = U_Q_G1(:,iZ1,iZ2,iZ3,iZ4,iS)
                U(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) = U_Q_G2(:,iZ1,iZ2,iZ3,iZ4,iS)
                U(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) = U_Q_G3(:,iZ1,iZ2,iZ3,iZ4,iS)

              ELSE IF ( NegativeStates(1,iZ1,iZ2,iZ3,iZ4,iS) ) THEN

                U(:,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) = U_Q_N (:,iZ1,iZ2,iZ3,iZ4,iS)

              END IF

            END DO
          END DO
        END DO
      END DO
    END DO

#ifdef THORNADO_DEBUG_POSITIVITYLIMITER
#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( U )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( U )
#endif
    WRITE(*,'(a20,7i4)')     'MAXLOC(U)', MAXLOC(U)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(U)', MAXVAL(U)
#endif

    CALL TimersStart( Timer_PL_Out )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from: U, iError, MinTheta_1, MinTheta_2 ) &
    !$OMP MAP( release: iZ_B0, iZ_E0, NegativeStates, &
    !$OMP               U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, &
    !$OMP               U_P_N, U_P_G1, U_P_G2, U_P_G3, &
    !$OMP               U_K_N, U_K_G1, U_K_G2, U_K_G3, &
    !$OMP               Min_K_S, Max_K_S, Theta_1_S, Min_Gam_S, Theta_2_S )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U, iError, MinTheta_1, MinTheta_2 ) &
    !$ACC DELETE( iZ_B0, iZ_E0, NegativeStates, &
    !$ACC         U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, &
    !$ACC         U_P_N, U_P_G1, U_P_G2, U_P_G3, &
    !$ACC         U_K_N, U_K_G1, U_K_G2, U_K_G3, &
    !$ACC         Min_K_S, Max_K_S, Theta_1_S, Min_Gam_S, Theta_2_S )
#endif

    CALL TimersStop( Timer_PL_Out )

    IF ( ANY( iError > 0 ) ) THEN
      DO iS = 1, nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iP = 1, nPT

                  IF ( iError(iP,iZ1,iZ2,iZ3,iZ4,iS) > 0 ) THEN

                    WRITE(*,'(A6,A,6I4)') &
                      '', 'SolveTheta_Bisection (M1):', iP, iZ1, iZ2, iZ3, iZ4, iS
                    WRITE(*,'(A8,A,I3.3)') &
                      '', 'ITERATION ', iError(iP,iZ1,iZ2,iZ3,iZ4,iS)
                    IF( iError(iP,iZ1,iZ2,iZ3,iZ4,iS) > MAX_IT + 3 ) STOP

                  ELSE IF ( iError(iP,iZ1,iZ2,iZ3,iZ4,iS) < 0 ) THEN

                    WRITE(*,*)
                    WRITE(*,'(A6,A,6I4)') &
                      '', 'SolveTheta_Bisection (M1):', iP, iZ1, iZ2, iZ3, iZ4, iS
                    WRITE(*,'(A8,A,I3.3)') &
                      '', 'Error: No Root in Interval'
                    WRITE(*,*)
                    STOP

                  END IF

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END IF

    CALL TimersStop( Timer_PositivityLimiter )

  END SUBROUTINE ApplyPositivityLimiter_TwoMoment


  SUBROUTINE TallyPositivityLimiter_TwoMoment( Time )

    USE UnitsModule, ONLY: UnitsDisplay

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    IF( .NOT. UsePositivityLimiterTally ) RETURN

    OPEN( NEWUNIT=FileUnit, FILE=TRIM( TallyFileName ), POSITION='APPEND', ACTION='WRITE' )

    WRITE( FileUnit, '(3(ES14.5,x))' ) &
      Time / UnitsDisplay % TimeUnit, MinTheta_1, MinTheta_2

    CLOSE( FileUnit )

  END SUBROUTINE TallyPositivityLimiter_TwoMoment


  SUBROUTINE ComputePointValues_Single( U_Q, U_P )

    REAL(DP), INTENT(in)  :: U_Q(nDOF)
    REAL(DP), INTENT(out) :: U_P(nPT)

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, 1, nDOF, One, L_X, nPT, &
             U_Q, nDOF, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues_Single


  SUBROUTINE ComputePointValues( iZ_B0, iZ_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: U_Q(nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP), INTENT(out) :: U_P(nPT ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, nR_1, nDOF, One, L_X, nPT, &
             U_Q, nDOF, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues


  FUNCTION GammaFun_Scalar( N, G1, G2, G3 ) &
      RESULT( GammaFun )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: N, G1, G2, G3
    REAL(DP) :: GammaFun

    GammaFun = ( One - Theta_FD * N ) * N - SQRT( G1**2 + G2**2 + G3**2 )

    RETURN
  END FUNCTION GammaFun_Scalar

  FUNCTION GammaFun_Vector( N, G1, G2, G3 ) &
      RESULT( GammaFun )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: N(:), G1(:), G2(:), G3(:)
    REAL(DP) :: GammaFun(SIZE(N))

    GammaFun = ( One - Theta_FD * N ) * N - SQRT( G1**2 + G2**2 + G3**2 )

    RETURN
  END FUNCTION GammaFun_Vector


  SUBROUTINE SolveTheta_Bisection( U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, U_K_N, U_K_G1, U_K_G2, U_K_G3, Theta_P, iError )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3
    REAL(DP), INTENT(in)  :: U_K_N, U_K_G1, U_K_G2, U_K_G3
    REAL(DP), INTENT(out) :: Theta_P
    INTEGER,  INTENT(out) :: iError

    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    iError = 0

    x_a = Zero
    f_a = GammaFun &
            ( x_a * U_Q_N  + ( One - x_a ) * U_K_N,   &
              x_a * U_Q_G1 + ( One - x_a ) * U_K_G1,  &
              x_a * U_Q_G2 + ( One - x_a ) * U_K_G2,  &
              x_a * U_Q_G3 + ( One - x_a ) * U_K_G3 ) &
          - Min_2

    x_b = One
    f_b = GammaFun &
            ( x_b * U_Q_N  + ( One - x_b ) * U_K_N,   &
              x_b * U_Q_G1 + ( One - x_b ) * U_K_G1,  &
              x_b * U_Q_G2 + ( One - x_b ) * U_K_G2,  &
              x_b * U_Q_G3 + ( One - x_b ) * U_K_G3 ) &
          - Min_2

    IF( .NOT. f_a * f_b < 0 )THEN

      iError = -1

    END IF

    dx = x_b - x_a

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx

      f_c = GammaFun &
              ( x_c * U_Q_N  + ( One - x_c ) * U_K_N,   &
                x_c * U_Q_G1 + ( One - x_c ) * U_K_G1,  &
                x_c * U_Q_G2 + ( One - x_c ) * U_K_G2,  &
                x_c * U_Q_G3 + ( One - x_c ) * U_K_G3 ) &
            - Min_2

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx < dx_min ) CONVERGED = .TRUE.

      IF( ITERATION > MAX_IT .AND. .NOT. CONVERGED )THEN

        iError = ITERATION
        IF ( ITERATION > MAX_IT + 3 ) EXIT

      END IF

    END DO

    Theta_P = x_a

  END SUBROUTINE SolveTheta_Bisection


END MODULE TwoMoment_PositivityLimiterModule
