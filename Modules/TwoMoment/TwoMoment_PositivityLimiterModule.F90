#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_POSITIVITYLIMITER
#endif
MODULE TwoMoment_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOF, nNodesZ, &
    nDOFE, nNodesE, &
    nDOFX, nNodesX
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
    MatrixMatrixMultiply
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, nDOF_X2, nDOF_X3, &
    Weights_q
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, LX_X1_Up, &
    LX_X2_Dn, LX_X2_Up, &
    LX_X3_Dn, LX_X3_Up
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_X1_Dn, L_X1_Up, &
    L_X2_Dn, L_X2_Up, &
    L_X3_Dn, L_X3_Up
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm, nGF, iGF_h_1, iGF_h_2, iGF_h_3
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
  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: UsePositivityLimiterTally
  INTEGER               :: N_R, N_G
  INTEGER,    PARAMETER :: nPS_Z = 9    ! Number of Positive Point Sets
  INTEGER,    PARAMETER :: nPS_X = 7    ! Number of Positive Point Sets
  INTEGER               :: nPP_Z(nPS_Z) ! Number of Positive Points Per Set
  INTEGER               :: nPP_X(nPS_X) ! Number of Positive Points Per Set
  INTEGER               :: nPT_Z        ! Number of Positive Points
  INTEGER               :: nPT_X        ! Number of Positive Points
  INTEGER               :: nR, nR_1, nZ(4)
  INTEGER,  ALLOCATABLE :: PointZ2X(:)
  REAL(DP)              :: Min_1, Max_1, Min_2
  REAL(DP)              :: Theta_FD, MinTheta_1, MinTheta_2
  REAL(DP),   PARAMETER :: Theta_Eps = 1.0_DP - 10.0_DP*EPSILON(1.0_DP)
  REAL(DP), ALLOCATABLE :: InterpMat_Z(:,:)
  REAL(DP), ALLOCATABLE :: InterpMat_X(:,:)

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
    INTEGER :: iNodeZ, iNodeX, iOS_Z, iOS_X, iP_Z, iP_X
    INTEGER :: iNodeE, iNodeX1, iNodeX2, iNodeX3

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

    ! --- Interpolation Matrix for Energy-Position Space Variables ---

    nPP_Z    = 0
    nPP_Z(1) = PRODUCT( nNodesZ )

    DO i = 2, 4  ! --- Exclude energy dimension for Order-1 ---

      IF( nNodesZ(i) > 1 )THEN

        nPP_Z(2*i:2*i+1) &
          = PRODUCT( nNodesZ, MASK = [1,2,3,4] .NE. i )

      END IF

    END DO

    ! --- Total Number of Positive Points ---

    nPT_Z = SUM( nPP_Z )

    ALLOCATE( InterpMat_Z(nPT_Z,nDOFZ) )

    InterpMat_Z = Zero
    DO iNodeZ = 1, nDOFZ

      InterpMat_Z(iNodeZ,iNodeZ) = One

      IF( SUM( nPP_Z(4:5) ) > 0 )THEN

        iOS_Z = SUM( nPP_Z(1:3) )
        InterpMat_Z(iOS_Z+1:iOS_Z+nDOF_X1,iNodeZ) = L_X1_Dn(1:nDOF_X1,iNodeZ)

        iOS_Z = iOS_Z + nPP_Z(4)
        InterpMat_Z(iOS_Z+1:iOS_Z+nDOF_X1,iNodeZ) = L_X1_Up(1:nDOF_X1,iNodeZ)

      END IF

      IF( SUM( nPP_Z(6:7) ) > 0 )THEN

        iOS_Z = SUM( nPP_Z(1:5) )
        InterpMat_Z(iOS_Z+1:iOS_Z+nDOF_X2,iNodeZ) = L_X2_Dn(1:nDOF_X2,iNodeZ)

        iOS_Z = iOS_Z + nPP_Z(6)
        InterpMat_Z(iOS_Z+1:iOS_Z+nDOF_X2,iNodeZ) = L_X2_Up(1:nDOF_X2,iNodeZ)

      END IF

      IF( SUM( nPP_Z(8:9) ) > 0 )THEN

        iOS_Z = SUM( nPP_Z(1:7) )
        InterpMat_Z(iOS_Z+1:iOS_Z+nDOF_X3,iNodeZ) = L_X3_Dn(1:nDOF_X3,iNodeZ)

        iOS_Z = iOS_Z + nPP_Z(8)
        InterpMat_Z(iOS_Z+1:iOS_Z+nDOF_X3,iNodeZ) = L_X3_Up(1:nDOF_X3,iNodeZ)

      END IF

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: InterpMat_Z )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( InterpMat_Z )
#endif

    ! --- Interpolation Matrix for Position Space Variables ---

    nPP_X    = 0
    nPP_X(1) = PRODUCT( nNodesX )

    DO i = 1, 3

      IF( nNodesX(i) > 1 )THEN

        nPP_X(2*i:2*i+1) = PRODUCT( nNodesX, MASK = [1,2,3] .NE. i )

      END IF

    END DO

    nPT_X = SUM( nPP_X )

    ALLOCATE( InterpMat_X(nPT_X,nDOFX) )

    InterpMat_X = Zero
    DO iNodeX = 1, nDOFX

      InterpMat_X(iNodeX,iNodeX) = One

      IF( SUM( nPP_X(2:3) ) > 0 )THEN

        iOS_X = nPP_X(1)
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X1,iNodeX) &
          = LX_X1_Dn(1:nDOFX_X1,iNodeX)

        iOS_X = iOS_X + nPP_X(2)
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X1,iNodeX) &
          = LX_X1_Up(1:nDOFX_X1,iNodeX)

      END IF

      IF( SUM( nPP_X(4:5) ) > 0 )THEN

        iOS_X = SUM( nPP_X(1:3) )
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X2,iNodeX) &
          = LX_X2_Dn(1:nDOFX_X2,iNodeX)

        iOS_X = iOS_X + nPP_X(4)
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X2,iNodeX) &
          = LX_X2_Up(1:nDOFX_X2,iNodeX)

      END IF

      IF( SUM( nPP_X(6:7) ) > 0 )THEN

        iOS_X = SUM( nPP_X(1:5) )
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X3,iNodeX) &
          = LX_X3_Dn(1:nDOFX_X3,iNodeX)

        iOS_X = iOS_X + nPP_X(6)
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X3,iNodeX) &
          = LX_X3_Up(1:nDOFX_X3,iNodeX)

      END IF

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: InterpMat_X )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( InterpMat_X )
#endif

    ! --- Index Map from Energy-Position to Position ---

    ALLOCATE( PointZ2X(nPT_Z) )

    iP_Z = 0
    iP_X = 0
    DO iNodeX3 = 1, nNodesX(3)
    DO iNodeX2 = 1, nNodesX(2)
    DO iNodeX1 = 1, nNodesX(1)

      iP_X = iP_X + 1

      DO iNodeE  = 1, nNodesE

        iP_Z = iP_Z + 1

        PointZ2X(iP_Z) = iP_X

      END DO

    END DO
    END DO
    END DO

    IF( nPP_Z(2) > 0 )THEN

      iP_Z = SUM( nPP_Z(1:1) )
      iP_X = 0
      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        iP_Z = iP_Z + 1
        iP_X = iP_X + 1

        PointZ2X(iP_Z) = iP_X

      END DO
      END DO
      END DO

    END IF

    IF( nPP_Z(3) > 0 )THEN

      iP_Z = SUM( nPP_Z(1:2) )
      iP_X = 0

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        iP_Z = iP_Z + 1
        iP_X = iP_X + 1

        PointZ2X(iP_Z) = iP_X

      END DO
      END DO
      END DO

    END IF

    IF( nPP_Z(4) > 0 )THEN

      iP_Z = SUM( nPP_Z(1:3) )
      iP_X = SUM( nPP_X(1:1) )

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)

        iP_X = iP_X + 1

        DO iNodeE = 1, nNodesE

          iP_Z = iP_Z + 1

          PointZ2X(iP_Z) = iP_X

        END DO

      END DO
      END DO

    END IF

    IF( nPP_Z(5) > 0 )THEN

      iP_Z = SUM( nPP_Z(1:4) )
      iP_X = SUM( nPP_X(1:2) )

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)

        iP_X = iP_X + 1

        DO iNodeE = 1, nNodesE

          iP_Z = iP_Z + 1

          PointZ2X(iP_Z) = iP_X

        END DO

      END DO
      END DO

    END IF

    IF( nPP_Z(6) > 0 )THEN

      iP_Z = SUM( nPP_Z(1:5) )
      iP_X = SUM( nPP_X(1:3) )

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX1 = 1, nNodesX(1)

        iP_X = iP_X + 1

        DO iNodeE = 1, nNodesE

          iP_Z = iP_Z + 1

          PointZ2X(iP_Z) = iP_X

        END DO

      END DO
      END DO

    END IF

    IF( nPP_Z(7) > 0 )THEN

      iP_Z = SUM( nPP_Z(1:6) )
      iP_X = SUM( nPP_X(1:4) )

      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX1 = 1, nNodesX(1)

        iP_X = iP_X + 1

        DO iNodeE = 1, nNodesE

          iP_Z = iP_Z + 1

          PointZ2X(iP_Z) = iP_X

        END DO

      END DO
      END DO

    END IF

    IF( nPP_Z(8) > 0 )THEN

      iP_Z = SUM( nPP_Z(1:7) )
      iP_X = SUM( nPP_X(1:5) )

      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        iP_X = iP_X + 1

        DO iNodeE = 1, nNodesE

          iP_Z = iP_Z + 1

          PointZ2X(iP_Z) = iP_X

        END DO

      END DO
      END DO

    END IF

    IF( nPP_Z(9) > 0 )THEN

      iP_Z = SUM( nPP_Z(1:8) )
      iP_X = SUM( nPP_X(1:6) )

      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        iP_X = iP_X + 1

        DO iNodeE = 1, nNodesE

          iP_Z = iP_Z + 1

          PointZ2X(iP_Z) = iP_X

        END DO

      END DO
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: PointZ2X )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( PointZ2X )
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

    IF ( ALLOCATED( InterpMat_Z ) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: InterpMat_Z )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA &
      !$ACC DELETE( InterpMat_Z )
#endif
      DEALLOCATE( InterpMat_Z )

    END IF

    IF ( ALLOCATED( InterpMat_X ) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: InterpMat_X )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA &
      !$ACC DELETE( InterpMat_X )
#endif
      DEALLOCATE( InterpMat_X )

    END IF

    IF ( ALLOCATED( PointZ2X ) ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: PointZ2X )
#elif defined(THORNADO_OACC)
      !$ACC EXIT DATA &
      !$ACC DELETE( PointZ2X )
#endif
      DEALLOCATE( PointZ2X )

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

    INTEGER  :: iNode, iZ1, iZ2, iZ3, iZ4, iS, iCR, iP_Z, iP_X, iNodeE, iNodeX
    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: nNeg_1, nNeg_2
    INTEGER  :: iX_B0(3), iX_E0(3)
    REAL(DP) :: Min_K, Max_K, Theta_1, Min_Gam, Theta_2, Gam, Theta_P

    REAL(DP) :: Min_K_S  (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: Max_K_S  (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: Theta_1_S(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: Min_Gam_S(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: Theta_2_S(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    INTEGER  :: iError(nPT_Z,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    LOGICAL  :: NegativeStates(2,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    REAL(DP) :: Tau_Q  (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: GX_Q_hd1(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: GX_Q_hd2(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: GX_Q_hd3(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: U_Q_N  (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_Q_G1 (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_Q_G2 (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_Q_G3 (nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    REAL(DP) :: GX_P_Gmdd11(nPT_X,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: GX_P_Gmdd22(nPT_X,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: GX_P_Gmdd33(nPT_X,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: U_P_N  (nPT_Z ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_P_G1 (nPT_Z ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_P_G2 (nPT_Z ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_P_G3 (nPT_Z ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    REAL(DP) :: U_K_N  (     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_K_G1 (     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_K_G2 (     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP) :: U_K_G3 (     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    IF( nDOF == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart( Timer_PositivityLimiter )

    nZ   = iZ_E0 - iZ_B0 + 1
    nR_1 = nSpecies * PRODUCT( nZ )
    N_R  = nR_1
    nR   = nR_1 * nCR

    iX_B0 = iZ_B0(2:4)
    iX_E0 = iZ_E0(2:4)
    N_G = PRODUCT( iX_E0 - iX_B0 + 1 )

    NegativeStates = .FALSE.
    iError = 0

    MinTheta_1 = One
    MinTheta_2 = One

    nNeg_1 = 0
    nNeg_2 = 0

    CALL TimersStart( Timer_PL_In )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, iX_B0, iX_E0, NegativeStates, &
    !$OMP          U, iError, MinTheta_1, MinTheta_2 )         &
    !$OMP MAP( alloc: Tau_Q, U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3,    &
    !$OMP             GX_Q_hd1, GX_Q_hd2, GX_Q_hd3,            &
    !$OMP             GX_P_Gmdd11, GX_P_Gmdd22, GX_P_Gmdd33,   &
    !$OMP             U_P_N, U_P_G1, U_P_G2, U_P_G3,           &
    !$OMP             U_K_N, U_K_G1, U_K_G2, U_K_G3,           &
    !$OMP             Min_K_S, Max_K_S, Theta_1_S, Min_Gam_S, Theta_2_S )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( iZ_B0, iZ_E0, iX_B0, iX_E0, NegativeStates, &
    !$ACC         U, iError, MinTheta_1, MinTheta_2 )         &
    !$ACC CREATE( Tau_Q, U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3,       &
    !$ACC         GX_Q_hd1, GX_Q_hd2, GX_Q_hd3,               &
    !$ACC         GX_P_Gmdd11, GX_P_Gmdd22, GX_P_Gmdd33,      &
    !$ACC         U_P_N, U_P_G1, U_P_G2, U_P_G3,              &
    !$ACC         U_K_N, U_K_G1, U_K_G2, U_K_G3,              &
    !$ACC         Min_K_S, Max_K_S, Theta_1_S, Min_Gam_S, Theta_2_S )
#endif

    CALL TimersStop( Timer_PL_In )

    ! --- Geometry Factor ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( GX_Q_hd1, GX_Q_hd2, GX_Q_hd3, GX, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

            GX_Q_hd1(iNodeX,iX1,iX2,iX3) = GX(iNodeX,iX1,iX2,iX3,iGF_h_1)
            GX_Q_hd2(iNodeX,iX1,iX2,iX3) = GX(iNodeX,iX1,iX2,iX3,iGF_h_2)
            GX_Q_hd3(iNodeX,iX1,iX2,iX3) = GX(iNodeX,iX1,iX2,iX3,iGF_h_3)

      END DO

    END DO
    END DO
    END DO

    ! --- Point Values ---

    CALL TimersStart( Timer_PL_Points )

    CALL ComputePointValuesX( iX_B0, iX_E0, GX_Q_hd1, GX_P_Gmdd11 )
    CALL ComputePointValuesX( iX_B0, iX_E0, GX_Q_hd2, GX_P_Gmdd22 )
    CALL ComputePointValuesX( iX_B0, iX_E0, GX_Q_hd3, GX_P_Gmdd33 )

    CALL TimersStop( Timer_PL_Points )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( GX_P_Gmdd11, GX_P_Gmdd22, GX_P_Gmdd33, iX_B0, iX_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(4)
#endif

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iP_X = 1, nPT_X

            GX_P_Gmdd11(iP_X,iX1,iX2,iX3) = MAX( GX_P_Gmdd11(iP_X,iX1,iX2,iX3)**2, SqrtTiny )
            GX_P_Gmdd22(iP_X,iX1,iX2,iX3) = MAX( GX_P_Gmdd22(iP_X,iX1,iX2,iX3)**2, SqrtTiny )
            GX_P_Gmdd33(iP_X,iX1,iX2,iX3) = MAX( GX_P_Gmdd33(iP_X,iX1,iX2,iX3)**2, SqrtTiny )

      END DO

   END DO
   END DO
   END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNode )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( iNode )                      &
    !$ACC PRESENT( Tau_Q, GX, GE, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNode )
#endif

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

          iNode = (iNodeX-1) * nDOFE + iNodeE

          Tau_Q(iNode,iZ1,iZ2,iZ3,iZ4) &
            = GE(iNodeE,iZ1,iGE_Ep2) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRESENT( iZ_B0, iZ_E0, U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6)
#endif

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iNode = 1, nDOF

      ! --- Copy Data ---

      U_Q_N (iNode,iZ1,iZ2,iZ3,iZ4,iS) = U(iNode,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)
      U_Q_G1(iNode,iZ1,iZ2,iZ3,iZ4,iS) = U(iNode,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)
      U_Q_G2(iNode,iZ1,iZ2,iZ3,iZ4,iS) = U(iNode,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)
      U_Q_G3(iNode,iZ1,iZ2,iZ3,iZ4,iS) = U(iNode,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Point Values ---

    CALL TimersStart( Timer_PL_Points )

    CALL ComputePointValuesZ( iZ_B0, iZ_E0, U_Q_N , U_P_N  )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, U_Q_G1, U_P_G1 )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, U_Q_G2, U_P_G2 )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, U_Q_G3, U_P_G3 )

    CALL TimersStop( Timer_PL_Points )

    ! --- Cell Average ---

    CALL TimersStart( Timer_PL_CellAverage )

    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, U_Q_N , U_K_N  )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, U_Q_G1, U_K_G1 )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, U_Q_G2, U_K_G2 )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, U_Q_G3, U_K_G3 )

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

              Min_K = Max_1
              Max_K = Min_1
              DO iP_Z = 1, nPT_Z
                Min_K = MIN( Min_K, U_P_N(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )
                Max_K = MAX( Max_K, U_P_N(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )
              END DO

              IF ( Min_K < Min_1 .OR. Max_K > Max_1 ) THEN

                NegativeStates(1,iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

                ! --- Limit Density Towards Cell Average ---

                Theta_1 &
                  = Theta_Eps * MIN( One, &
                    ABS( (Min_1-U_K_N(iZ1,iZ2,iZ3,iZ4,iS)) &
                          / (Min_K-U_K_N(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny) ), &
                    ABS( (Max_1-U_K_N(iZ1,iZ2,iZ3,iZ4,iS)) &
                          / (Max_K-U_K_N(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny) ) )

                DO iNode = 1, nDOF
                  U_Q_N(iNode,iZ1,iZ2,iZ3,iZ4,iS) &
                    =   U_Q_N(iNode,iZ1,iZ2,iZ3,iZ4,iS) * Theta_1 &
                      + U_K_N(iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_1 )
                END DO

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

    CALL ComputePointValuesZ( iZ_B0, iZ_E0, U_Q_N , U_P_N  )

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
    !$OMP PRIVATE( Gam, Min_Gam, Theta_2, Theta_P, iP_X ) &
    !$OMP REDUCTION( min: MinTheta_2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( Gam, Min_Gam, Theta_2, Theta_P, iP_X ) &
    !$ACC REDUCTION( min: MinTheta_2 ) &
    !$ACC PRESENT( U_P_N, U_P_G1, U_P_G2, U_P_G3, &
    !$ACC          GX_P_Gmdd11, GX_P_Gmdd22, GX_P_Gmdd33, &
    !$ACC          U_K_N, U_K_G1, U_K_G2, U_K_G3, &
    !$ACC          U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, &
    !$ACC          NegativeStates, Min_Gam_S, Theta_2_S, &
    !$ACC          MinTheta_2, iError, iZ_B0, iZ_E0, PointZ2X )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( Gam, Min_Gam, Theta_2, Theta_P, iP_X ) &
    !$OMP REDUCTION( min: MinTheta_2 )
#endif
    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              Min_Gam = Min_2
              Theta_2 = One

              DO iP_Z = 1, nPT_Z

                iP_X = PointZ2X(iP_Z)
!! Shaoping we got compile error without -O0 as  "Total size of kernel arguments exceeds limit!"               
                Gam = GammaFun &
                      ( U_P_N (iP_Z,iZ1,iZ2,iZ3,iZ4,iS),   &
                        U_P_G1(iP_Z,iZ1,iZ2,iZ3,iZ4,iS),   &
                        U_P_G2(iP_Z,iZ1,iZ2,iZ3,iZ4,iS),   &
                        U_P_G3(iP_Z,iZ1,iZ2,iZ3,iZ4,iS),   &
                        GX_P_Gmdd11(iP_X,iZ2,iZ3,iZ4), &
                        GX_P_Gmdd22(iP_X,iZ2,iZ3,iZ4), &
                        GX_P_Gmdd33(iP_X,iZ2,iZ3,iZ4) )

                Min_Gam = MIN( Min_Gam, Gam )

                IF ( Gam < Min_2 ) THEN

                  CALL SolveTheta_Bisection &
                        ( U_P_N  (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                          U_P_G1 (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                          U_P_G2 (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                          U_P_G3 (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                          U_K_N  (   iZ1,iZ2,iZ3,iZ4,iS),   &
                          U_K_G1 (   iZ1,iZ2,iZ3,iZ4,iS),   &
                          U_K_G2 (   iZ1,iZ2,iZ3,iZ4,iS),   &
                          U_K_G3 (   iZ1,iZ2,iZ3,iZ4,iS),   &
                          GX_P_Gmdd11(iP_X,iZ2,iZ3,iZ4),    &
                          GX_P_Gmdd22(iP_X,iZ2,iZ3,iZ4),    &
                          GX_P_Gmdd33(iP_X,iZ2,iZ3,iZ4),    &
                          Theta_P,                          &
                          iError (iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )

                  Theta_2 = MIN( Theta_2, Theta_P )

                END IF

              END DO

              IF ( Min_Gam < Min_2 ) THEN

                NegativeStates(2,iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

                ! --- Limit Towards Cell Average ---

                !Theta_2 = Theta_Eps * MINVAL( Theta_P(:) )
                Theta_2 = Theta_Eps * Theta_2

                DO iNode = 1, nDOF
                  U_Q_N (iNode,iZ1,iZ2,iZ3,iZ4,iS) &
                    =   U_Q_N (iNode,iZ1,iZ2,iZ3,iZ4,iS) * Theta_2 &
                      + U_K_N (iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_2 )

                  U_Q_G1(iNode,iZ1,iZ2,iZ3,iZ4,iS) &
                    =   U_Q_G1(iNode,iZ1,iZ2,iZ3,iZ4,iS) * Theta_2 &
                      + U_K_G1(iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_2 )

                  U_Q_G2(iNode,iZ1,iZ2,iZ3,iZ4,iS) &
                    =   U_Q_G2(iNode,iZ1,iZ2,iZ3,iZ4,iS) * Theta_2 &
                      + U_K_G2(iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_2 )

                  U_Q_G3(iNode,iZ1,iZ2,iZ3,iZ4,iS) &
                    =   U_Q_G3(iNode,iZ1,iZ2,iZ3,iZ4,iS) * Theta_2 &
                      + U_K_G3(iZ1,iZ2,iZ3,iZ4,iS) * ( One - Theta_2 )
                END DO

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
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)

              ! --- Update If Needed ---

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
    !$OMP MAP( release: iZ_B0, iZ_E0, iX_B0, iX_E0, NegativeStates, &
    !$OMP               Tau_Q, U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3,       &
    !$OMP               GX_Q_hd1, GX_Q_hd2, GX_Q_hd3,               &
    !$OMP               GX_P_Gmdd11, GX_P_Gmdd22, GX_P_Gmdd33,      &
    !$OMP               U_P_N, U_P_G1, U_P_G2, U_P_G3,              &
    !$OMP               U_K_N, U_K_G1, U_K_G2, U_K_G3,              &
    !$OMP               Min_K_S, Max_K_S, Theta_1_S, Min_Gam_S, Theta_2_S )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT( U, iError, MinTheta_1, MinTheta_2 )        &
    !$ACC DELETE( iZ_B0, iZ_E0, iX_B0, iX_E0, NegativeStates, &
    !$ACC         Tau_Q, U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3,       &
    !$ACC         GX_Q_hd1, GX_Q_hd2, GX_Q_hd3,               &
    !$ACC         GX_P_Gmdd11, GX_P_Gmdd22, GX_P_Gmdd33,      &
    !$ACC         U_P_N, U_P_G1, U_P_G2, U_P_G3,              &
    !$ACC         U_K_N, U_K_G1, U_K_G2, U_K_G3,              &
    !$ACC         Min_K_S, Max_K_S, Theta_1_S, Min_Gam_S, Theta_2_S )
#endif

    CALL TimersStop( Timer_PL_Out )

    IF ( ANY( iError > 0 ) ) THEN
      DO iS = 1, nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
            DO iZ2 = iZ_B0(2), iZ_E0(2)
              DO iZ1 = iZ_B0(1), iZ_E0(1)
                DO iP_Z = 1, nPT_Z

                  IF ( iError(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) > 0 ) THEN

                    WRITE(*,'(A6,A,6I4)') &
                      '', 'SolveTheta_Bisection (M1):', iP_Z, iZ1, iZ2, iZ3, iZ4, iS
                    WRITE(*,'(A8,A,I3.3)') &
                      '', 'ITERATION ', iError(iP_Z,iZ1,iZ2,iZ3,iZ4,iS)
                    IF( iError(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) > MAX_IT + 3 ) STOP

                  ELSE IF ( iError(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) < 0 ) THEN

                    WRITE(*,*)
                    WRITE(*,'(A6,A,6I4)') &
                      '', 'SolveTheta_Bisection (M1):', iP_Z, iZ1, iZ2, iZ3, iZ4, iS
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


  SUBROUTINE ComputePointValuesZ_Single( U_Q, U_P )

    REAL(DP), INTENT(in)  :: U_Q(nDOF)
    REAL(DP), INTENT(out) :: U_P(nPT_Z)

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT_Z, 1, nDOFZ, One, InterpMat_Z, nPT_Z, &
             U_Q, nDOFZ, Zero, U_P, nPT_Z )

  END SUBROUTINE ComputePointValuesZ_Single


  SUBROUTINE ComputePointValuesZ( iZ_B0, iZ_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: U_Q(nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP), INTENT(out) :: U_P(nPT_Z ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT_Z, N_R, nDOFZ, One, InterpMat_Z, nPT_Z, &
             U_Q, nDOFZ, Zero, U_P, nPT_Z )

  END SUBROUTINE ComputePointValuesZ


  SUBROUTINE ComputePointValuesX( iX_B0, iX_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: U_Q(nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: U_P(nPT_X,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT_X, N_G, nDOFX, One, InterpMat_X, nPT_X, &
             U_Q, nDOFX, Zero, U_P, nPT_X )

  END SUBROUTINE ComputePointValuesX


  SUBROUTINE ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, U_Q, U_K )

    INTEGER,  INTENT(in)  :: iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: Tau_Q(nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in)  :: U_Q(nDOF,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP), INTENT(out) :: U_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iNode
    REAL(DP) :: SUM1, SUM2

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( SUM1, SUM2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRIVATE( SUM1, SUM2 ) &
    !$ACC PRESENT( U_K, Weights_q, Tau_Q, U_Q, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( SUM1, SUM2 )
#endif

    DO iS = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
            DO iZ1 = iZ_B0(1), iZ_E0(1)
              SUM1 = Zero
              SUM2 = Zero
              DO iNode = 1, nDOF
                SUM1 = SUM1 + Weights_q(iNode) * Tau_Q(iNode,iZ1,iZ2,iZ3,iZ4) * U_Q(iNode,iZ1,iZ2,iZ3,iZ4,iS)
                SUM2 = SUM2 + Weights_q(iNode) * Tau_Q(iNode,iZ1,iZ2,iZ3,iZ4)
              END DO
              U_K(iZ1,iZ2,iZ3,iZ4,iS) = SUM1 / SUM2

            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ComputeCellAverage


  FUNCTION GammaFun_Scalar( N, G1, G2, G3, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
      RESULT( GammaFun )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: N, G1, G2, G3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: GammaFun

    GammaFun = ( One - Theta_FD * N ) * N &
               - SQRT( G1**2 / Gm_dd_11 + G2**2 / Gm_dd_22 + G3**2 / Gm_dd_33 )

    RETURN
  END FUNCTION GammaFun_Scalar

  FUNCTION GammaFun_Vector( N, G1, G2, G3, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
      RESULT( GammaFun )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: N(:), G1(:), G2(:), G3(:), &
                            Gm_dd_11(:), Gm_dd_22(:), Gm_dd_33(:)
    REAL(DP) :: GammaFun(SIZE(N))

    GammaFun = ( One - Theta_FD * N ) * N &
               - SQRT( G1**2 / Gm_dd_11 + G2**2 / Gm_dd_22 + G3**2 / Gm_dd_33 )

    RETURN
  END FUNCTION GammaFun_Vector


  SUBROUTINE SolveTheta_Bisection &
    ( U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3, U_K_N, U_K_G1, U_K_G2, U_K_G3, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Theta_P, iError )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: U_Q_N, U_Q_G1, U_Q_G2, U_Q_G3
    REAL(DP), INTENT(in)  :: U_K_N, U_K_G1, U_K_G2, U_K_G3
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
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
              x_a * U_Q_G3 + ( One - x_a ) * U_K_G3,  &
              Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
          - Min_2

    x_b = One
    f_b = GammaFun &
            ( x_b * U_Q_N  + ( One - x_b ) * U_K_N,   &
              x_b * U_Q_G1 + ( One - x_b ) * U_K_G1,  &
              x_b * U_Q_G2 + ( One - x_b ) * U_K_G2,  &
              x_b * U_Q_G3 + ( One - x_b ) * U_K_G3,  &
              Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
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
                x_c * U_Q_G3 + ( One - x_c ) * U_K_G3,  &
                Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
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
