#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_POSITIVITYLIMITER
#endif
MODULE TwoMoment_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, nDOF, nDOFE, nDOFX
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_PositivityLimiter
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
  REAL(DP),   PARAMETER :: Theta_Eps = 1.0_DP - EPSILON(1.0_DP)
  REAL(DP), ALLOCATABLE :: L_X(:,:)

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
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_1 = ', Min_1
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_1 = ', Max_1
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_2 = ', Min_2
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Theta_FD = ', Theta_FD

    END IF

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

!#if defined(THORNADO_OMP_OL)
!      !$OMP TARGET ENTER DATA &
!      !$OMP MAP( to: L_X )
!#elif defined(THORNADO_OACC)
!      !$ACC ENTER DATA &
!      !$ACC COPYIN( L_X )
!#endif

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

!#if defined(THORNADO_OMP_OL)
!      !$OMP TARGET EXIT DATA &
!      !$OMP MAP( release: L_X )
!#elif defined(THORNADO_OACC)
!      !$ACC EXIT DATA &
!      !$ACC DELETE( L_X )
!#endif
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
    INTEGER  :: iPack_1, iPack_2, nPack_1, nPack_2, iter
    REAL(DP) :: Min_K, Max_K, Theta_1, Theta_2, Theta_P(nPT), Gam(nPT)

    LOGICAL  :: NegativeStates(2,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    REAL(DP) :: U_q    (nDOF,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_PP   (nPT ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_K    (     nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    INTEGER,  ALLOCATABLE :: NegativeTable_1(:,:)
    INTEGER,  ALLOCATABLE :: NegativeTable_2(:,:)

    REAL(DP), ALLOCATABLE :: U_K_Pack1(:)
    REAL(DP), ALLOCATABLE :: U_Q_Pack1(:,:)
    REAL(DP), ALLOCATABLE :: U_P_Pack1(:,:)
    REAL(DP), ALLOCATABLE :: Min_K_Pack(:)
    REAL(DP), ALLOCATABLE :: Max_K_Pack(:)

    REAL(DP), ALLOCATABLE :: U_K_Pack2(:,:)
    REAL(DP), ALLOCATABLE :: U_Q_Pack2(:,:,:)
    REAL(DP), ALLOCATABLE :: U_P_Pack2(:,:,:)
    REAL(DP), ALLOCATABLE :: Gam_Pack(:,:)
    REAL(DP), ALLOCATABLE :: Theta_2_Pack(:)

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart( Timer_PositivityLimiter )

    nZ = iZ_E0 - iZ_B0 + 1
    nR_1 = nSpecies * PRODUCT( nZ )
    nR = nR_1 * nCR

    ALLOCATE( Min_K_Pack(nR_1) )
    ALLOCATE( Max_K_Pack(nR_1) )
    ALLOCATE( U_K_Pack1(nR_1) )
    ALLOCATE( U_Q_Pack1(nDOF,nR_1) )
    ALLOCATE( U_P_Pack1(nPT,nR_1) )

    ALLOCATE( U_K_Pack2(nCR,nR_1) )
    ALLOCATE( U_Q_Pack2(nDOF,nCR,nR_1) )
    ALLOCATE( U_P_Pack2(nCR,nPT,nR_1) )
    ALLOCATE( Gam_Pack(nPT,nR_1) )
    ALLOCATE( Theta_2_Pack(nR_1) )

    ALLOCATE( NegativeTable_1(5,nR_1) )
    ALLOCATE( NegativeTable_2(5,nR_1) )

    MinTheta_1 = One
    MinTheta_2 = One

    NegativeStates = .FALSE.
    NegativeTable_1 = 0
    NegativeTable_2 = 0

!#if defined(THORNADO_OMP_OL)
!    !$OMP TARGET ENTER DATA &
!    !$OMP MAP( to: U, iZ_B0, iZ_E0, MinTheta_1, MinTheta_2, NegativeStates, &
!    !$OMP          NegativeTable_1, NegativeTable_2 ) &
!    !$OMP MAP( alloc: Min_K_Pack, Max_K_Pack, U_K_Pack1, U_Q_Pack1, U_P_Pack1, &
!    !$OMP             U_K_Pack2, U_Q_Pack2, U_P_Pack2, Gam_Pack, Theta_2_Pack, &
!    !$OMP             U_q, U_PP, U_K, Theta_P, Gam )
!#elif defined(THORNADO_OACC)
!    !$ACC ENTER DATA &
!    !$ACC COPYIN( U, iZ_B0, iZ_E0, MinTheta_1, MinTheta_2, NegativeStates, &
!    !$ACC         NegativeTable_1, NegativeTable_2 ) &
!    !$ACC CREATE( Min_K_Pack, Max_K_Pack, U_K_Pack1, U_Q_Pack1, U_P_Pack1, &
!    !$ACC         U_K_Pack2, U_Q_Pack2, U_P_Pack2, Gam_Pack, Theta_2_Pack, &
!    !$ACC         U_q, U_PP, U_K, Theta_P, Gam )
!#endif

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iNode = 1, nDOF

                  U_q(iNode,iCR,iZ1,iZ2,iZ3,iZ4,iS) = U(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS)

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    ! --- Ensure Bounded Density ---

    CALL ComputePointValues( iZ_B0, iZ_E0, U_q, U_PP )

    iPack_1 = 0
    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              Min_K = MINVAL( U_PP(:,iCR_N,iZ1,iZ2,iZ3,iZ4,iS) )
              Max_K = MAXVAL( U_PP(:,iCR_N,iZ1,iZ2,iZ3,iZ4,iS) )

              IF ( Min_K < Min_1 .OR. Max_K > Max_1 ) THEN

                iPack_1 = iPack_1 + 1

                NegativeTable_1(:,iPack_1) = [ iZ1, iZ2, iZ3, iZ4, iS ]

                NegativeStates(1,iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

                Min_K_Pack(iPack_1) = Min_K
                Max_K_Pack(iPack_1) = Max_K

                U_Q_Pack1(:,iPack_1) = U_q(:,iCR_N,iZ1,iZ2,iZ3,iZ4,iS)

              END IF

            END DO
          END DO
        END DO
      END DO
    END DO
    nPack_1 = iPack_1

    IF ( nPack_1 > 0 ) THEN

      ! --- Cell Average ---

      CALL DGEMV &
        ( 'T', nDOF, nPack_1, One, U_Q_Pack1, nDOF, &
          Weights_q, 1, Zero, U_K_Pack1, 1 )

      ! --- Limit Density Towards Cell Average ---

      DO iPack_1 = 1, nPack_1

        Min_K = Min_K_pack(iPack_1)
        Max_K = Max_K_pack(iPack_1)

        DO iNode = 1, nDOF

          Theta_1 &
            = Theta_Eps * MIN( One, &
                   ABS( (Min_1-U_K_Pack1(iPack_1)) / (Min_K-U_K_Pack1(iPack_1)) ), &
                   ABS( (Max_1-U_K_Pack1(iPack_1)) / (Max_K-U_K_Pack1(iPack_1)) ) )

          U_Q_Pack1(iNode,iPack_1) &
            = U_Q_Pack1(iNode,iPack_1) * Theta_1 + U_K_Pack1(iPack_1) * ( One - Theta_1 )

          MinTheta_1 = MIN( Theta_1, MinTheta_1 )

        END DO

      END DO

      ! --- Recompute Point Values ---

      CALL ComputePointValues_Pack( nPack_1, U_Q_Pack1, U_P_Pack1 )

      DO iPack_1 = 1, nPack_1

        iZ1 = NegativeTable_1(1,iPack_1)
        iZ2 = NegativeTable_1(2,iPack_1)
        iZ3 = NegativeTable_1(3,iPack_1)
        iZ4 = NegativeTable_1(4,iPack_1)
        iS  = NegativeTable_1(5,iPack_1)

        U_q (:,iCR_N,iZ1,iZ2,iZ3,iZ4,iS) = U_Q_Pack1(:,iPack_1)
        U_PP(:,iCR_N,iZ1,iZ2,iZ3,iZ4,iS) = U_P_Pack1(:,iPack_1)

      END DO

#ifdef THORNADO_DEBUG_POSITIVITYLIMITER
      DO iPack_1 = 1, MIN( nPack_1, 5 )
        Min_K = Min_K_pack(iPack_1)
        Max_K = Max_K_pack(iPack_1)
        Theta_1 &
          = MIN( One, &
                 ABS( (Min_1-U_K_Pack1(iPack_1)) / (Min_K-U_K_Pack1(iPack_1)) ), &
                 ABS( (Max_1-U_K_Pack1(iPack_1)) / (Max_K-U_K_Pack1(iPack_1)) ) )
        WRITE(*,'(a30,5i4)')      'Neg. UQ(N) @ ', NegativeTable_1(:,iPack_1)
        WRITE(*,'(a30,3es23.15)')    'Min_K, Max_K, Theta_1 : ', Min_K, Max_K, Theta_1
        WRITE(*,'(a30,es23.15,i4)')  '        MINVAL(UQ(N)) : ', MINVAL(U_Q_Pack1(:,iPack_1)), MINLOC(U_Q_Pack1(:,iPack_1))
        WRITE(*,'(a30,es23.15)')     '                  U_K : ', U_K_Pack1(iPack_1)
        WRITE(*,'(a30,es23.15,i4)')  '        MINVAL(UP(N)) : ', MINVAL(U_P_Pack1(:,iPack_1)), MINLOC(U_P_Pack1(:,iPack_1))
      END DO
#endif

    END IF

    ! --- Ensure Positive Gamma ---

    iPack_2 = 0
    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              CALL ComputeGamma &
                     ( nPT, U_PP(1:nPT,1:nCR,iZ1,iZ2,iZ3,iZ4,iS), Gam(1:nPT) )

              IF ( ANY( Gam(:) < Min_2 ) ) THEN

                iPack_2 = iPack_2 + 1

                NegativeTable_2(:,iPack_2) = [ iZ1, iZ2, iZ3, iZ4, iS ]

                NegativeStates(1,iZ1,iZ2,iZ3,iZ4,iS) = .FALSE.
                NegativeStates(2,iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

                DO iCR = 1, nCR
                  U_Q_Pack2(:,iCR,iPack_2) = U_q (:,iCR,iZ1,iZ2,iZ3,iZ4,iS)
                  U_P_Pack2(iCR,:,iPack_2) = U_PP(:,iCR,iZ1,iZ2,iZ3,iZ4,iS)
                END DO
                Gam_Pack(:,iPack_2) = Gam(:)

              END IF

            END DO
          END DO
        END DO
      END DO
    END DO
    nPack_2 = iPack_2

    IF ( nPack_2 > 0 ) THEN

      ! --- Cell Average ---

      CALL DGEMV &
        ( 'T', nDOF, nPack_2*nCR, One, U_Q_Pack2, nDOF, &
          Weights_q, 1, Zero, U_K_Pack2, 1 )

      ! --- Limit Towards Cell Average ---

      DO iPack_2 = 1, nPack_2

        Theta_P(:) = One
        DO iP = 1, nPT
          IF( Gam_Pack(iP,iPack_2) < Min_2 ) THEN

            CALL SolveTheta_Bisection &
                   ( U_P_Pack2(1:nCR,iP,iPack_2), U_K_Pack2(1:nCR,iPack_2), Min_2, Theta_P(iP) )

          END IF
        END DO

        Theta_2 = Theta_Eps * MINVAL( Theta_P(:) )
        Theta_2_Pack(iPack_2) = Theta_2

        MinTheta_2 = MIN( Theta_2, MinTheta_2 )

      END DO

      DO iPack_2 = 1, nPack_2
        DO iCR = 1, nCR
          DO iNode = 1, nDOF

            Theta_2 = Theta_2_Pack(iPack_2)

            U_Q_Pack2(iNode,iCR,iPack_2) &
              = U_Q_Pack2(iNode,iCR,iPack_2) * Theta_2 + U_K_Pack2(iCR,iPack_2) * ( One - Theta_2 )

          END DO
        END DO
      END DO

      DO iPack_2 = 1, nPack_2
        DO iCR = 1, nCR

          iZ1 = NegativeTable_2(1,iPack_2)
          iZ2 = NegativeTable_2(2,iPack_2)
          iZ3 = NegativeTable_2(3,iPack_2)
          iZ4 = NegativeTable_2(4,iPack_2)
          iS  = NegativeTable_2(5,iPack_2)

          U_q (:,iCR,iZ1,iZ2,iZ3,iZ4,iS) = U_Q_Pack2(:,iCR,iPack_2)

        END DO
      END DO

#ifdef THORNADO_DEBUG_POSITIVITYLIMITER
      DO iPack_2 = 1, MIN( nPack_2, 5 )
        Theta_2 = Theta_2_Pack(iPack_2)
        WRITE(*,'(a30,5i4)')      'Neg. Gamma @ ', NegativeTable_2(:,iPack_2)
        WRITE(*,'(a30,2es23.15)') '   Min_Gamma, Theta_2 : ', MINVAL(Gam_Pack(:,iPack_2)), Theta_2
        WRITE(*,'(a30,4es23.15)') '        MINVAL(UQ(:)) : ', MINVAL(U_Q_Pack2(:,:,iPack_2),DIM=1)
        WRITE(*,'(a30,4es23.15)') '                  U_K : ', U_K_Pack2(:,iPack_2)
      END DO
#endif

    END IF

    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              IF ( NegativeStates(1,iZ1,iZ2,iZ3,iZ4,iS) ) THEN

                U(:,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) = U_q(:,iCR_N,iZ1,iZ2,iZ3,iZ4,iS)

              ELSE IF ( NegativeStates(2,iZ1,iZ2,iZ3,iZ4,iS) ) THEN

                U(:,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) = U_q(:,1:nCR,iZ1,iZ2,iZ3,iZ4,iS)

              END IF

            END DO
          END DO
        END DO
      END DO
    END DO

!#if defined(THORNADO_OMP_OL)
!    !$OMP TARGET EXIT DATA &
!    !$OMP MAP( from: U, MinTheta_1, MinTheta_2 ) &
!    !$OMP MAP( release: iZ_B0, iZ_E0, NegativeTable_1, NegativeTable_2, &
!    !$OMP               Min_K_Pack, Max_K_Pack, U_K_Pack1, U_Q_Pack1, U_P_Pack1, &
!    !$OMP               U_K_Pack2, U_Q_Pack2, U_P_Pack2, Gam_Pack, Theta_2_Pack, &
!    !$OMP               U_q, U_PP, U_K, Theta_P, Gam, NegativeStates )
!#elif defined(THORNADO_OACC)
!    !$ACC EXIT DATA &
!    !$ACC COPYOUT( U, MinTheta_1, MinTheta_2 ) &
!    !$ACC DELETE( iZ_B0, iZ_E0, NegativeTable_1, NegativeTable_2, &
!    !$ACC         Min_K_Pack, Max_K_Pack, U_K_Pack1, U_Q_Pack1, U_P_Pack1, &
!    !$ACC         U_K_Pack2, U_Q_Pack2, U_P_Pack2, Gam_Pack, Theta_2_Pack, &
!    !$ACC         U_q, U_PP, U_K, Theta_P, Gam, NegativeStates )
!#endif

    CALL TimersStop( Timer_PositivityLimiter )

    DEALLOCATE( U_K_Pack1, U_Q_Pack1, U_P_Pack1, NegativeTable_1 )
    DEALLOCATE( Min_K_Pack, Max_K_Pack )
    DEALLOCATE( U_K_Pack2, U_Q_Pack2, U_P_Pack2, NegativeTable_2 )
    DEALLOCATE( Gam_Pack, Theta_2_Pack )

#ifdef THORNADO_DEBUG_POSITIVITYLIMITER
    WRITE(*,'(a20,7i4)')     'MAXLOC(U)', MAXLOC(U)
    WRITE(*,'(a20,es23.15)') 'MAXVAL(U)', MAXVAL(U)
#endif

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
    REAL(DP), INTENT(inout) :: U_P(nPT)

    CALL DGEMM &
           ( 'N', 'N', nPT, 1, nDOF, One, L_X, nPT, &
             U_Q, nDOF, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues_Single


  SUBROUTINE ComputePointValues_Pack( nPack, U_Q, U_P )

    INTEGER,  INTENT(in)  :: nPack
    REAL(DP), INTENT(in)  :: U_Q(nDOF,nPack)
    REAL(DP), INTENT(inout) :: U_P(nPT ,nPack)

    CALL DGEMM &
           ( 'N', 'N', nPT, nPack, nDOF, One, L_X, nPT, &
             U_Q, nDOF, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues_Pack


  SUBROUTINE ComputePointValues( iZ_B0, iZ_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: U_Q(nDOF,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP), INTENT(inout) :: U_P(nPT ,nCR,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    CALL DGEMM &
           ( 'N', 'N', nPT, nR, nDOF, One, L_X, nPT, &
             U_Q, nDOF, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues


  SUBROUTINE ComputeGamma( N, U, Gam )

    INTEGER,  INTENT(in)  :: N
    REAL(DP), INTENT(in)  :: U(N,nCR)
    REAL(DP), INTENT(out) :: Gam(N)

    Gam = GammaFun( U(:,iCR_N), U(:,iCR_G1), U(:,iCR_G2), U(:,iCR_G3) )

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

      WRITE(*,*)
      WRITE(*,'(A6,A)') &
        '', 'SolveTheta_Bisection (M1):'
      WRITE(*,'(A8,A,I3.3)') &
        '', 'Error: No Root in Interval'
      WRITE(*,'(A8,A,2ES15.6e3)') &
        '', 'x_a, x_b = ', x_a, x_b
      WRITE(*,'(A8,A,2ES15.6e3)') &
        '', 'f_a, f_b = ', f_a, f_b
      WRITE(*,*)
      WRITE(*,'(A8,A,ES20.12e3)') &
        '', 'N_K  = ', U_K(iCR_N)
      WRITE(*,'(A8,A,ES20.12e3)') &
        '', 'G1_K = ', U_K(iCR_G1)
      WRITE(*,'(A8,A,ES20.12e3)') &
        '', 'G2_K = ', U_K(iCR_G2)
      WRITE(*,'(A8,A,ES20.12e3)') &
        '', 'G3_K = ', U_K(iCR_G3)
      WRITE(*,*)

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


END MODULE TwoMoment_PositivityLimiterModule
