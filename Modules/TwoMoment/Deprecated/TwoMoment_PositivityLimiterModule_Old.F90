#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_POSITIVITYLIMITER
#endif
MODULE TwoMoment_PositivityLimiterModule_Old

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, nDOF, nDOFE, nDOFX
  USE TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_PositivityLimiter
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_X1_Dn, L_X1_Up, &
    L_X2_Dn, L_X2_Up, &
    L_X3_Dn, L_X3_Up
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_SqrtGm, nGF
  USE GeometryFieldsModuleE, ONLY: &
    iGE_Ep2, nGE
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
  LOGICAL               :: Debug = .FALSE.
  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: UsePositivityLimiterTally
  INTEGER,    PARAMETER :: nPS = 9  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Number of Positive Points
  REAL(DP)              :: Min_1, Max_1, Min_2
  REAL(DP)              :: Theta_FD, MinTheta_1, MinTheta_2
  REAL(DP),   PARAMETER :: Theta_Eps = 1.0_DP - EPSILON(1.0_DP)
  REAL(DP), ALLOCATABLE :: U_PP(:,:), GX_PP(:,:)
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
    INTEGER :: i, iNode, iOS, FileUnit

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
      WRITE(*,'(A2,A6,A)') '', 'INFO: ', 'InitializePositivityLimiter_TwoMoment'
      WRITE(*,*)
      WRITE(*,'(A6,A,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter 
      WRITE(*,*)
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_1 = ', Min_1
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Max_1 = ', Max_1
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Min_2 = ', Min_2
      WRITE(*,'(A6,A12,ES12.4E3)') '', 'Theta_FD = ', Theta_FD
    END IF

    ! --- Compute Number of Positive Points per Set ---

    nPP = 0
    nPP(1) = PRODUCT( nNodesZ )

    DO i = 2, 4 ! --- Exclude energy dimension for now ---

      IF( nNodesZ(i) > 1 )THEN

        nPP(2*i:2*i+1) &
          = PRODUCT( nNodesZ, MASK = [1,2,3,4] .NE. i )

      END IF

    END DO

    ! --- Total Number of Positive Points ---

    nPT = SUM( nPP )

    ! --- Conserved Variables in Positive Points ---

    ALLOCATE( U_PP(nPT,nCR) )

    ! --- Geometry in Positive Points ---

    ALLOCATE( GX_PP(nPT,nGF) )

    ! ---                             ---

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
 
    ! --- For Tally of Positivity Limiter ---

    IF( UsePositivityLimiterTally )THEN

      TallyFileName = '../Output/PositivityLimiterTally.dat'

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ) )

      WRITE( FileUnit, '(3(A14,x))' ) 'Time', 'MinTheta_1', 'MinTheta_2'

      CLOSE( FileUnit )

    END IF

  END SUBROUTINE InitializePositivityLimiter_TwoMoment


  SUBROUTINE FinalizePositivityLimiter_TwoMoment

    IF( ALLOCATED( U_PP ) )THEN
       DEALLOCATE( U_PP )
    END IF
    IF( ALLOCATED( GX_PP ) )THEN
       DEALLOCATE( GX_PP )
    END IF
    IF( ALLOCATED( L_X ) )THEN
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

    LOGICAL  :: NegativeStates(2)
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iCR, iP
    INTEGER  :: iNode, iNodeE, iNodeX
    INTEGER  :: nNeg_1, nNeg_2
    REAL(DP) :: Min_K, Max_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: U_q(nDOF,nCR), U_K(nCR), Gamma(nPT)
    REAL(DP) :: GX_q(nDOF,nGF)
    REAL(DP) :: Tau_q(nDOF), Tau_K
    REAL(DP), EXTERNAL :: DDOT

    IF( nDOF == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart( Timer_PositivityLimiter )

    MinTheta_1 = One
    MinTheta_2 = One

    nNeg_1 = 0
    nNeg_2 = 0

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      ! --- Copy Data ---

      U_q(1:nDOF,1:nCR) = U(1:nDOF,iZ1,iZ2,iZ3,iZ4,1:nCR,iS)

      DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE
          iNode = iNodeE+(iNodeX-1)*nDOFE
          ! --- Compute Tau and Geometry factor ---
          Tau_q(iNode) = &
            GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) * GE(iNodeE,iZ1,iGE_Ep2)
          GX_q(iNode,1:nGF) = GX(iNodeX,iZ2,iZ3,iZ4,1:nGF)
        END DO ! iNodeE
      END DO ! iNodeX

      Tau_K = DDOT( nDOF, Weights_q, 1, Tau_q, 1 )

      NegativeStates = .FALSE.

      ! --- Point Values ---

      CALL ComputePointValues( nCR, U_q, U_PP )
      CALL ComputePointValues( nGF, GX_q, GX_PP )
      CALL ComputeGeometryX_FromScaleFactors( GX_PP )


      ! --- Ensure Bounded Density ---

      Min_K = MINVAL( U_PP(:,iCR_N) )
      Max_K = MAXVAL( U_PP(:,iCR_N) )

      IF( Min_K < Min_1 .OR. Max_K > Max_1 )THEN

        IF( Debug ) WRITE(*,*) 'Positivity trigger 1128 MinMax'

        ! --- Cell Average ---

        !U_K(iCR_N) = DOT_PRODUCT( Weights_q, U_q(:,iCR_N) )
        !U_K(iCR_N) = DDOT( nDOF, Weights_q, 1, U_q(:,iCR_N), 1 )
        CALL DGEMV( 'T', nDOF, 1, One, U_q(:,iCR_N), nDOF, Weights_q*Tau_q, 1, Zero, U_K(iCR_N), 1 )
        U_K(iCR_N) = U_K(iCR_N) / Tau_K 

        Theta_1 &
          = Theta_Eps * MIN( One, &
                 ABS( (Min_1-U_K(iCR_N)) / (Min_K-U_K(iCR_N) + SqrtTiny) ), &
                 ABS( (Max_1-U_K(iCR_N)) / (Max_K-U_K(iCR_N) + SqrtTiny) ) )

        ! --- Limit Density Towards Cell Average ---

        U_q(:,iCR_N) = Theta_1 * U_q(:,iCR_N) + ( One - Theta_1 ) * U_K(iCR_N)

        ! --- Recompute Point Values ---

        CALL ComputePointValues( nCR, U_q, U_PP )

        NegativeStates(1) = .TRUE.

      END IF

#ifdef THORNADO_DEBUG_POSITIVITYLIMITER
      IF( NegativeStates(1) )THEN
        nNeg_1 = nNeg_1 + 1
        IF ( nNeg_1 <= 5 ) THEN
          WRITE(*,'(a30,5i4)')         'Neg. UQ(N) @ ', iZ1, iZ2, iZ3, iZ4, iS
          WRITE(*,'(a30,3es23.15)')    'Min_K, Max_K, Theta_1 : ', Min_K, Max_K, Theta_1
          WRITE(*,'(a30,es23.15,i4)')  '        MINVAL(UQ(N)) : ', MINVAL(U_q(:,iCR_N)), MINLOC(U_q(:,iCR_N))
          WRITE(*,'(a30,es23.15)')     '                  U_K : ', U_K(iCR_N)
        END IF
      END IF
#endif

      ! --- Ensure Positive Gamma ---

      Gamma(1:nPT) = GammaFun &
                     ( U_PP(1:nPT,iCR_N ), &
                       U_PP(1:nPT,iCR_G1), &
                       U_PP(1:nPT,iCR_G2), &
                       U_PP(1:nPT,iCR_G3), &
                       GX_PP(1:nPT,iGF_Gm_dd_11), &
                       GX_PP(1:nPT,iGF_Gm_dd_22), &
                       GX_PP(1:nPT,iGF_Gm_dd_33) )

      IF( ANY( Gamma(:) < Min_2 ) )THEN

        IF( Debug ) WRITE(*,*) 'Positivity trigger 1129 Gamma'

        ! --- Cell Average ---

        !DO iCR = 1, nCR
        !  U_K(iCR) = DOT_PRODUCT( Weights_q, U_q(:,iCR) )
        !  !U_K(iCR) = DDOT( nDOF, Weights_q, 1, U_q(:,iCR), 1 )
        !END DO
        CALL DGEMV( 'T', nDOF, nCR, One, U_q, nDOF, Weights_q*Tau_q, 1, Zero, U_K, 1 )
        U_K = U_K / Tau_K

        Theta_2 = One
        DO iP = 1, nPT

          IF( Gamma(iP) < Min_2 ) THEN

            CALL SolveTheta_Bisection &
                   ( U_PP(iP,1:nCR), U_K(1:nCR), &
                     GX_PP(iP,iGF_Gm_dd_11), &
                     GX_PP(iP,iGF_Gm_dd_22), &
                     GX_PP(iP,iGF_Gm_dd_33), &
                     Min_2, Theta_P )

            Theta_2 = MIN( Theta_2, Theta_P )

          END IF

        END DO

        Theta_2 = Theta_Eps * Theta_2

        ! --- Limit Towards Cell Average ---

        DO iCR = 1, nCR

          U_q(:,iCR) = Theta_2 * U_q(:,iCR) + ( One - Theta_2 ) * U_K(iCR)

        END DO

        NegativeStates(2) = .TRUE.

        MinTheta_2 = MIN( Theta_2, MinTheta_2 )

      END IF

#ifdef THORNADO_DEBUG_POSITIVITYLIMITER
      IF( NegativeStates(2) )THEN
        nNeg_2 = nNeg_2 + 1
        IF ( nNeg_2 <= 5 ) THEN
          WRITE(*,'(a30,5i4)')      'Neg. Gamma @ ', iZ1, iZ2, iZ3, iZ4, iS
          WRITE(*,'(a30,2es23.15)') '   Min_Gamma, Theta_2 : ', MINVAL(Gamma(:)), Theta_2
          WRITE(*,'(a30,4es23.15)') '        MINVAL(UQ(:)) : ', MINVAL(U_q(:,:),DIM=1)
          WRITE(*,'(a30,4es23.15)') '                  U_K : ', U_K(:)
        END IF
      END IF
#endif

      ! --- Update If Needed ---

      IF( NegativeStates(2) )THEN

        U(1:nDOF,iZ1,iZ2,iZ3,iZ4,1:nCR,iS) = U_q(1:nDOF,1:nCR)

      ELSEIF( NegativeStates(1) )THEN

        U(1:nDOF,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) = U_q(1:nDOF,iCR_N)

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_PositivityLimiter )

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


  SUBROUTINE ComputePointValues( nDim, U_Q, U_P )

    INTEGER,  INTENT(in)  :: nDim
    REAL(DP), INTENT(in)  :: U_Q(nDOF,nDim)
    REAL(DP), INTENT(out) :: U_P(nPT, nDim)

    ! --- Compute U_P = One * L_X * U_Q ---

    CALL DGEMM &
           ( 'N', 'N', nPT, nDim, nDOF, One, L_X, nPT, &
             U_Q, nDOF, Zero, U_P, nPT )


  END SUBROUTINE ComputePointValues


  PURE REAL(DP) ELEMENTAL FUNCTION GammaFun &
    ( N, G1, G2, G3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in) :: N, G1, G2, G3, Gm_dd_11, Gm_dd_22, Gm_dd_33

    GammaFun = ( One - Theta_FD * N ) * N &
               - SQRT( G1**2 / Gm_dd_11 + G2**2 / Gm_dd_22 + G3**2 / Gm_dd_33 )

    RETURN
  END FUNCTION GammaFun


  SUBROUTINE SolveTheta_Bisection &
    ( U_Q, U_K, Gm_dd_11, Gm_dd_22, Gm_dd_33, MinGamma, Theta_P )

    REAL(DP), INTENT(in)  :: U_Q(nCR), U_K(nCR), MinGamma
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
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
              x_a * U_Q(iCR_G3) + ( One - x_a ) * U_K(iCR_G3),  &
              Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
          - MinGamma

    x_b = One
    f_b = GammaFun &
            ( x_b * U_Q(iCR_N)  + ( One - x_b ) * U_K(iCR_N),   &
              x_b * U_Q(iCR_G1) + ( One - x_b ) * U_K(iCR_G1),  &
              x_b * U_Q(iCR_G2) + ( One - x_b ) * U_K(iCR_G2),  &
              x_b * U_Q(iCR_G3) + ( One - x_b ) * U_K(iCR_G3),  &
              Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
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
                x_c * U_Q(iCR_G3) + ( One - x_c ) * U_K(iCR_G3),  &
                Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
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


END MODULE TwoMoment_PositivityLimiterModule_Old
