MODULE TwoMoment_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny, FourPi
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nNodesZ, &
    nDOFE, nNodesE, &
    nDOFX, nNodesX
  USE TwoMoment_TimersModule, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_PL, &
    Timer_PL_Permute, &
    Timer_PL_PointValues, &
    Timer_PL_CellAverage, &
    Timer_PL_Theta_1, &
    Timer_PL_Theta_2, &
    Timer_PL_EnergyLimiter
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, LX_X1_Up, &
    LX_X2_Dn, LX_X2_Up, &
    LX_X3_Dn, LX_X3_Up
  USE ReferenceElementModule, ONLY: &
    nDOF_E, &
    nDOF_X1, &
    nDOF_X2, &
    nDOF_X3, &
    Weights_Q
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_E_Dn,  L_E_Up, &
    L_X1_Dn, L_X1_Up, &
    L_X2_Dn, L_X2_Up, &
    L_X3_Dn, L_X3_Up
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep2, iGE_Ep3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm, iGF_h_1, iGF_h_2, iGF_h_3, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    nDR, iDR_PL_Theta_1, iDR_PL_Theta_2, iDR_PL_dEnergy

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_TwoMoment
  PUBLIC :: FinalizePositivityLimiter_TwoMoment
  PUBLIC :: ApplyPositivityLimiter_TwoMoment

  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: UseEnergyLimiter
  LOGICAL               :: Verbose
  INTEGER               :: N_R, N_G
  INTEGER,    PARAMETER :: nPS_Z = 9    ! Number of Positive Point Sets
  INTEGER,    PARAMETER :: nPS_X = 7    ! Number of Positive Point Sets
  INTEGER               :: nPP_Z(nPS_Z) ! Number of Positive Points Per Set
  INTEGER               :: nPP_X(nPS_X) ! Number of Positive Points Per Set
  INTEGER               :: nPT_Z        ! Total Number of Positive Points
  INTEGER               :: nPT_X        ! Total Number of Positive Points
  INTEGER,  ALLOCATABLE :: PointZ2X(:)
  REAL(DP)              :: Min_1, Max_1, Min_2
  REAL(DP)              :: W_Factor
  REAL(DP),   PARAMETER :: One_EPS = One - 1.0d1 * EPSILON( One )
  REAL(DP), ALLOCATABLE :: InterpMat_Z(:,:)
  REAL(DP), ALLOCATABLE :: InterpMat_X(:,:)

  REAL(DP), PUBLIC :: dEnergyMomentum_PL_TwoMoment(nCR)

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE &
  !$OMP TARGET( Min_1, Max_1, Min_2, W_Factor, nPT_Z )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE &
  !$ACC CREATE( Min_1, Max_1, Min_2, W_Factor, nPT_Z )
#endif

CONTAINS


  SUBROUTINE InitializePositivityLimiter_TwoMoment &
    ( Min_1_Option, Max_1_Option, Min_2_Option, UsePositivityLimiter_Option, &
      UseEnergyLimiter_Option, Verbose_Option )

    USE UnitsModule, ONLY: &
      UnitsActive, &
      PlanckConstant, &
      SpeedOfLight

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Max_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UseEnergyLimiter_Option
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: i, iNodeZ, iNodeX, iOS_Z, iOS_X, iP_Z, iP_X
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

    IF( PRESENT( UseEnergyLimiter_Option ) )THEN
      UseEnergyLimiter = UseEnergyLimiter_Option
    ELSE
      UseEnergyLimiter = .FALSE.
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN

      WRITE(*,*)
      WRITE(*,'(A)') '  INFO: InitializePositivityLimiter_TwoMoment:'
      WRITE(*,'(A)') '  --------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A32,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter
      WRITE(*,'(A4,A32,L1)') &
        '', 'Use Energy Correction: ' , UseEnergyLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A32,ES11.3E3)') '', 'Min_1: ', Min_1
      WRITE(*,'(A4,A32,ES11.3E3)') '', 'Max_1: ', Max_1
      WRITE(*,'(A4,A32,ES11.3E3)') '', 'Min_2: ', Min_2

    END IF

    IF( UnitsActive )THEN

      W_Factor = FourPi / ( PlanckConstant * SpeedOfLight )**3

    ELSE

      W_Factor = FourPi

    END IF

    ! --- Interpolation Matrix for Energy-Position Space Variables ---

    nPP_Z    = 0
    nPP_Z(1) = PRODUCT( nNodesZ )

    DO i = 1, 4

      IF( nNodesZ(i) > 1 )THEN

        nPP_Z(2*i:2*i+1) = PRODUCT( nNodesZ, MASK = [1,2,3,4] .NE. i )

      END IF

    END DO

    nPT_Z = SUM( nPP_Z )

    ALLOCATE( InterpMat_Z(nPT_Z,nDOFZ) )

    InterpMat_Z = Zero
    DO iNodeZ = 1, nDOFZ

      InterpMat_Z(iNodeZ,iNodeZ) = One

      IF( SUM( nPP_Z(2:3) ) > 0 )THEN

        iOS_Z = nPP_Z(1)
        InterpMat_Z(iOS_Z+1:iOS_Z+nDOF_E,iNodeZ) = L_E_Dn(1:nDOF_E,iNodeZ)

        iOS_Z = iOS_Z + nPP_Z(2)
        InterpMat_Z(iOS_Z+1:iOS_Z+nDOF_E,iNodeZ) = L_E_Up(1:nDOF_E,iNodeZ)

      END IF

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
    !$OMP TARGET UPDATE &
    !$OMP TO( Min_1, Max_1, Min_2, W_Factor, nPT_Z )

    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: InterpMat_Z, InterpMat_X, PointZ2X )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE &
    !$ACC DEVICE( Min_1, Max_1, Min_2, W_Factor, nPT_Z )

    !$ACC ENTER DATA &
    !$ACC COPYIN( InterpMat_Z, InterpMat_X, PointZ2X )
#endif

    dEnergyMomentum_PL_TwoMoment = Zero

  END SUBROUTINE InitializePositivityLimiter_TwoMoment


  SUBROUTINE FinalizePositivityLimiter_TwoMoment

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: InterpMat_Z, InterpMat_X, PointZ2X )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE( InterpMat_Z, InterpMat_X, PointZ2X )
#endif

    DEALLOCATE( InterpMat_Z )
    DEALLOCATE( InterpMat_X )
    DEALLOCATE( PointZ2X )

  END SUBROUTINE FinalizePositivityLimiter_TwoMoment


  SUBROUTINE ApplyPositivityLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, uDR_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE (1:nDOFE, &
          iZ_B1(1):iZ_E1(1), &
          1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR, &
          1:nSpecies)
    REAL(DP), INTENT(out), OPTIONAL :: &
      uDR_Option &
         (iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nDR)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP_Z, iP_X
    INTEGER  :: iNodeZ, iNodeE, iNodeX
    REAL(DP) :: Min_K, Max_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: Gamma_P, Gamma_Min
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: EnergyMomentum_0(nCR)
    REAL(DP) :: EnergyMomentum_1(nCR)

    LOGICAL,  ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: RealizableCellAverage, LimiterApplied
    LOGICAL,  ALLOCATABLE, DIMENSION(:,:,:,:)     :: ApplyEnergyLimiter
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)       :: Theta_1_K, Theta_2_K
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:)     :: h_d_1_Q, h_d_1_P, h_d_2_Q, h_d_2_P, h_d_3_Q, h_d_3_P
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: Energy_K, dEnergy_K
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: Tau_Q, N_K, G1_K, G2_K, G3_K
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: N_Q, N_P, G1_Q, G1_P, G2_Q, G2_P, G3_Q, G3_P

    ALLOCATE( RealizableCellAverage(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( LimiterApplied(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( ApplyEnergyLimiter(iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( Theta_1_K(iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( Theta_2_K(iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( Energy_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( dEnergy_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( Tau_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( N_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( N_P(nPT_Z,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( N_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G1_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G1_P(nPT_Z,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G1_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G2_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G2_P(nPT_Z,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G2_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G3_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G3_P(nPT_Z,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( G3_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies) )
    ALLOCATE( h_d_1_Q(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( h_d_1_P(nPT_X,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( h_d_2_Q(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( h_d_2_P(nPT_X,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( h_d_3_Q(nDOFX,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )
    ALLOCATE( h_d_3_P(nPT_X,iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)) )

    IF( .NOT. UsePositivityLimiter .OR. nDOFZ == 1 ) RETURN

    CALL TimersStart( Timer_PL )

    N_R = nSpecies * PRODUCT( iZ_E0 - iZ_B0 + 1 )

    N_G = PRODUCT( iZ_E0(2:4) - iZ_B0(2:4) + 1 )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iZ_B0, iZ_E0, GE, GX, U_F, U_R ) &
    !$OMP MAP( alloc:  RealizableCellAverage, &
    !$OMP              LimiterApplied, ApplyEnergyLimiter, &
    !$OMP              Energy_K, dEnergy_K, &
    !$OMP              Tau_Q, N_Q, N_P, N_K, G1_Q, G1_P, G1_K, &
    !$OMP              G2_Q, G2_P, G2_K, G3_Q, G3_P, G3_K, h_d_1_Q, h_d_1_P, &
    !$OMP              h_d_2_Q, h_d_2_P, h_d_3_Q, h_d_3_P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA ASYNC &
    !$ACC COPYIN( iZ_B0, iZ_E0, GE, GX, U_F, U_R ) &
    !$ACC CREATE( RealizableCellAverage, &
    !$ACC         LimiterApplied, ApplyEnergyLimiter, &
    !$ACC         Energy_K, dEnergy_K, &
    !$ACC         Tau_Q, N_Q, N_P, N_K, G1_Q, G1_P, G1_K, &
    !$ACC         G2_Q, G2_P, G2_K, G3_Q, G3_P, G3_K, h_d_1_Q, h_d_1_P, &
    !$ACC         h_d_2_Q, h_d_2_P, h_d_3_Q, h_d_3_P )
#endif

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      Theta_1_K(iZ2,iZ3,iZ4) = One
      Theta_2_K(iZ2,iZ3,iZ4) = One

    END DO
    END DO
    END DO

    CALL ComputeGlobalEnergyMomentum & ! --- For Global Tally ---
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, EnergyMomentum_0 )

    CALL TimersStart( Timer_PL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, GX, h_d_1_Q, h_d_2_Q, h_d_3_Q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        h_d_1_Q(iNodeX,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_h_1)
        h_d_2_Q(iNodeX,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_h_2)
        h_d_3_Q(iNodeX,iZ2,iZ3,iZ4) = GX(iNodeX,iZ2,iZ3,iZ4,iGF_h_3)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop( Timer_PL_Permute )

    CALL ComputePointValuesX( iZ_B0, iZ_E0, h_d_1_Q, h_d_1_P )
    CALL ComputePointValuesX( iZ_B0, iZ_E0, h_d_2_Q, h_d_2_P )
    CALL ComputePointValuesX( iZ_B0, iZ_E0, h_d_3_Q, h_d_3_P )

    CALL TimersStart( Timer_PL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC PRIVATE( iNodeZ ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, Tau_Q, GX, GE )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( iNodeZ )
#endif
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        Tau_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4) &
          = GE(iNodeE,iZ1,iGE_Ep2) * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, N_Q, G1_Q, G2_Q, G3_Q, U_R )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(6)
#endif
    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        N_Q (iNodeZ,iZ1,iZ2,iZ3,iZ4,iS) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)
        G1_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4,iS) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)
        G2_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4,iS) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)
        G3_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4,iS) &
          = U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)

     END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_PL_Permute )

    CALL ComputePointValuesZ( iZ_B0, iZ_E0, N_Q , N_P  )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G1_Q, G1_P )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G2_Q, G2_P )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G3_Q, G3_P )

    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, N_Q , N_K  )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G1_Q, G1_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G2_Q, G2_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G3_Q, G3_K )

    CALL CheckCellAverageRealizability &
           ( iZ_B0, iZ_E0, N_K, G1_K, G2_K, G3_K, &
             h_d_1_P, h_d_2_P, h_d_3_P, RealizableCellAverage )

    ! --- Ensure Bounded Density ---

    CALL TimersStart( Timer_PL_Theta_1 )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(5) &
    !$OMP PRIVATE( Min_K, Max_K, Theta_1 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG COLLAPSE(5) ASYNC &
    !$ACC PRIVATE( Min_K, Max_K, Theta_1 ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, LimiterApplied, &
    !$ACC          RealizableCellAverage, N_P, N_K, N_Q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( Min_K, Max_K, Theta_1 )
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      LimiterApplied(iZ1,iZ2,iZ3,iZ4,iS) = .FALSE.

      IF( RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) )THEN

        Min_K = Min_1
        Max_K = Max_1

#if   defined( THORNADO_OMP_OL )
        !$OMP PARALLEL DO SIMD &
        !$OMP REDUCTION( min: Min_K ) &
        !$OMP REDUCTION( max: Max_K )
#elif defined( THORNADO_OACC   )
        !$ACC LOOP VECTOR &
        !$ACC REDUCTION( min: Min_K ) &
        !$ACC REDUCTION( max: Max_K )
#endif
        DO iP_Z = 1, nPT_Z

          Min_K = MIN( Min_K, N_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )
          Max_K = MAX( Max_K, N_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )

        END DO

        IF( Min_K < Min_1 .OR. Max_K > Max_1 )THEN

          Theta_1 &
            = MIN( One, &
                   ABS( ( Min_1 - N_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                      / ( Min_K - N_K(iZ1,iZ2,iZ3,iZ4,iS)-SqrtTiny ) ), &
                   ABS( ( Max_1 - N_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                      / ( Max_K - N_K(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny ) ) )

          Theta_1 = One_EPS * Theta_1

          CALL ApplyLimiter &
                 ( Theta_1, N_K(iZ1,iZ2,iZ3,iZ4,iS), &
                   N_Q(:,iZ1,iZ2,iZ3,iZ4,iS) )

          CALL ComputePointValuesZ_Single &
                 ( InterpMat_Z, &
                   N_Q(:,iZ1,iZ2,iZ3,iZ4,iS), &
                   N_P(:,iZ1,iZ2,iZ3,iZ4,iS) )

          LimiterApplied(iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

          Theta_1_K(iZ2,iZ3,iZ4) = MIN( Theta_1, Theta_1_K(iZ2,iZ3,iZ4) )

        END IF

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_PL_Theta_1 )

    ! --- Ensure Positive "Gamma" ---

    CALL TimersStart( Timer_PL_Theta_2 )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(5) &
    !$OMP PRIVATE( Theta_2, Gamma_Min )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG COLLAPSE(5) ASYNC &
    !$ACC PRIVATE( Theta_2, Gamma_Min ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, RealizableCellAverage, LimiterApplied, &
    !$ACC          h_d_1_P, h_d_2_P, h_d_3_P, N_P, G1_P, G2_P, G3_P, &
    !$ACC          N_K, G1_K, G2_K, G3_K, N_Q, G1_Q, G2_Q, G3_Q, PointZ2X )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( Theta_2, Gamma_Min, &
    !$OMP          iP_X, Gm_dd_11, Gm_dd_22, Gm_dd_33, Gamma_P, Theta_P )
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      IF( RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) )THEN

        Theta_2   = One
        Gamma_Min = Min_2

#if   defined( THORNADO_OMP_OL )
        !$OMP PARALLEL DO SIMD &
        !$OMP PRIVATE( iP_X, Gm_dd_11, Gm_dd_22, Gm_dd_33, Gamma_P, Theta_P ) &
        !$OMP REDUCTION( min: Gamma_Min, Theta_2 )
#elif defined( THORNADO_OACC   )
        !$ACC LOOP VECTOR &
        !$ACC PRIVATE( iP_X, Gm_dd_11, Gm_dd_22, Gm_dd_33, Gamma_P, Theta_P ) &
        !$ACC REDUCTION( min: Gamma_Min, Theta_2 )
#endif
        DO iP_Z = 1, nPT_Z

          iP_X = PointZ2X(iP_Z)

          Gm_dd_11 = MAX( h_d_1_P(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )
          Gm_dd_22 = MAX( h_d_2_P(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )
          Gm_dd_33 = MAX( h_d_3_P(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )

          Gamma_P &
            = GammaFun &
                ( N_P (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                  G1_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                  G2_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                  G3_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                  Gm_dd_11, Gm_dd_22, Gm_dd_33 )

          Gamma_Min = MIN( Gamma_P, Gamma_Min )

          IF( Gamma_P < Min_2 )THEN

            CALL SolveTheta_Bisection &
                   ( N_P (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                     G1_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                     G2_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                     G3_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                     N_K (     iZ1,iZ2,iZ3,iZ4,iS), &
                     G1_K(     iZ1,iZ2,iZ3,iZ4,iS), &
                     G2_K(     iZ1,iZ2,iZ3,iZ4,iS), &
                     G3_K(     iZ1,iZ2,iZ3,iZ4,iS), &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33 , &
                     Theta_P )

            Theta_2 = MIN( Theta_2, Theta_P )

          END IF

        END DO

        IF( Gamma_Min < Min_2 )THEN

          ! --- Limit Towards Cell Average ---

          Theta_2 = One_EPS * Theta_2

          CALL ApplyLimiter &
                 ( Theta_2, N_K (iZ1,iZ2,iZ3,iZ4,iS), &
                   N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) )
          CALL ApplyLimiter &
                 ( Theta_2, G1_K(iZ1,iZ2,iZ3,iZ4,iS), &
                   G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) )
          CALL ApplyLimiter &
                 ( Theta_2, G2_K(iZ1,iZ2,iZ3,iZ4,iS), &
                   G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) )
          CALL ApplyLimiter &
                 ( Theta_2, G3_K(iZ1,iZ2,iZ3,iZ4,iS), &
                   G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) )

          LimiterApplied(iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

          Theta_2_K(iZ2,iZ3,iZ4) = MIN( Theta_2, Theta_2_K(iZ2,iZ3,iZ4) )

        END IF

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_PL_Theta_2 )

    CALL RecoverRealizableCellAverage &
           ( iZ_B0, iZ_E0, N_K, G1_K, G2_K, G3_K, N_Q, G1_Q, G2_Q, G3_Q, &
             h_d_1_P, h_d_2_P, h_d_3_P, RealizableCellAverage )

    CALL TimersStart( Timer_PL_EnergyLimiter )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, LimiterApplied, ApplyEnergyLimiter )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4)
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      ApplyEnergyLimiter(iZ2,iZ3,iZ4,iS) &
        = ANY( LimiterApplied(:,iZ2,iZ3,iZ4,iS) )

    END DO
    END DO
    END DO
    END DO

    CALL ComputeLocalEnergy & ! --- For Energy Correction ---
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, &
             ApplyEnergyLimiter, Energy_K )

    CALL TimersStop( Timer_PL_EnergyLimiter )

    CALL TimersStart( Timer_PL_Permute )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, U_R, N_Q, G1_Q, G2_Q, G3_Q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) = N_Q (:,iZ1,iZ2,iZ3,iZ4,iS)
      U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) = G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS)
      U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) = G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS)
      U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) = G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_PL_Permute )

    CALL TimersStart( Timer_PL_EnergyLimiter )

    CALL ComputeLocalEnergy & ! --- For Energy Correction ---
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, &
             ApplyEnergyLimiter, dEnergy_K )

    ! --- Energy Change Due to Positivity Limiter ---

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRESENT( iZ_B0, iZ_E0, ApplyEnergyLimiter, dEnergy_K, Energy_K )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5)
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      IF( ApplyEnergyLimiter(iZ2,iZ3,iZ4,iS) )THEN

        dEnergy_K(iZ1,iZ2,iZ3,iZ4,iS) &
          = dEnergy_K(iZ1,iZ2,iZ3,iZ4,iS) - Energy_K(iZ1,iZ2,iZ3,iZ4,iS)

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_PL_EnergyLimiter )

    CALL ApplyEnergyLimiter_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, &
             ApplyEnergyLimiter, dEnergy_K )

    CALL ComputeGlobalEnergyMomentum & ! --- For Global Tally ---
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, EnergyMomentum_1 )

    IF( PRESENT( uDR_Option ) )THEN

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC  )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
      !$ACC PRESENT( uDR_Option, iZ_B0, iZ_E0 )
#elif defined(THORNADO_OMP   )
      !$OMP PARALLEL DO COLLAPSE(3)
#endif
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)

        uDR_Option(iZ2,iZ3,iZ4,iDR_PL_Theta_1) = Theta_1_K(iZ2,iZ3,iZ4)
        uDR_Option(iZ2,iZ3,iZ4,iDR_PL_Theta_2) = Theta_2_K(iZ2,iZ3,iZ4)
        uDR_Option(iZ2,iZ3,iZ4,iDR_PL_dEnergy) = Zero

      END DO
      END DO
      END DO

#if   defined(THORNADO_OMP_OL)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC  )
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
      !$ACC PRESENT( uDR_Option, iZ_B0, iZ_E0, ApplyEnergyLimiter, dEnergy_K )
#elif defined(THORNADO_OMP   )
      !$OMP PARALLEL DO COLLAPSE(4)
#endif
      DO iS  = 1, nSpecies
      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iZ1 = iZ_B0(1), iZ_E0(1)

        IF( ApplyEnergyLimiter(iZ2,iZ3,iZ4,iS) )THEN

          uDR_Option(iZ2,iZ3,iZ4,iDR_PL_dEnergy) &
            = uDR_Option(iZ2,iZ3,iZ4,iDR_PL_dEnergy) &
                + dEnergy_K(iZ1,iZ2,iZ3,iZ4,iS)

        END IF

      END DO
      END DO
      END DO
      END DO
      END DO

    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: iZ_B0, iZ_E0, GE, GX, U_F, U_R, &
    !$OMP               RealizableCellAverage, &
    !$OMP               LimiterApplied, ApplyEnergyLimiter, &
    !$OMP               Energy_K, dEnergy_K, Tau_Q, N_Q, N_P, N_K, G1_Q, G1_P, G1_K, &
    !$OMP               G2_Q, G2_P, G2_K, G3_Q, G3_P, G3_K, h_d_1_Q, h_d_1_P, &
    !$OMP               h_d_2_Q, h_d_2_P, h_d_3_Q, h_d_3_P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA ASYNC &
    !$ACC DELETE( iZ_B0, iZ_E0, GE, GX, U_F, U_R, &
    !$ACC         RealizableCellAverage, &
    !$ACC         LimiterApplied, ApplyEnergyLimiter, &
    !$ACC         Energy_K, dEnergy_K, Tau_Q, N_Q, N_P, N_K, G1_Q, G1_P, G1_K, &
    !$ACC         G2_Q, G2_P, G2_K, G3_Q, G3_P, G3_K, h_d_1_Q, h_d_1_P, &
    !$ACC         h_d_2_Q, h_d_2_P, h_d_3_Q, h_d_3_P )

    !$ACC WAIT
#endif

    dEnergyMomentum_PL_TwoMoment = EnergyMomentum_1 - EnergyMomentum_0

    CALL TimersStop( Timer_PL )

  END SUBROUTINE ApplyPositivityLimiter_TwoMoment


  SUBROUTINE ApplyEnergyLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, ApplyEnergyLimiter, DeltaE )

    USE MeshModule, ONLY: &
      MeshE, MeshX

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE (1:nDOFE, &
          iZ_B1(1):iZ_E1(1), &
          1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR, &
          1:nSpecies)
    LOGICAL, INTENT(in)     :: &
      ApplyEnergyLimiter &
        (iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies)
    REAL(DP), INTENT(in)    :: &
      DeltaE &
        (iZ_B0(1):iZ_E0(1), &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies)

    INTEGER  :: iNodeE, iNodeX, iNodeZ
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iK1, iK2
    REAL(DP) :: N_K1, N_K2, E_K1, E_K2, ResidualE
    REAL(DP) :: Theta_K1, Theta_K2
    REAL(DP) :: MinTheta_K1, MinTheta_K2

    REAL(DP), ALLOCATABLE, DIMENSION(:)   :: V_u_1, V_u_2, V_u_3
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: W2_K, W3_K

    ALLOCATE( V_u_1(1:nDOFX) )
    ALLOCATE( V_u_2(1:nDOFX) )
    ALLOCATE( V_u_3(1:nDOFX) )
    ALLOCATE( W2_K(1:nDOFZ,iZ_B0(1):iZ_E0(1)) )
    ALLOCATE( W3_K(1:nDOFZ,iZ_B0(1):iZ_E0(1)) )

    IF( .NOT. UseEnergyLimiter ) RETURN

    CALL TimersStart( Timer_PL_EnergyLimiter )

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(4) &
    !$OMP MAP( to: dZ1, dZ2, dZ3, dZ4 ) &
    !$OMP PRIVATE( ResidualE, iK1, iK2, N_K1, N_K2, E_K1, E_K2, &
    !$OMP          Theta_K1, Theta_K2, MinTheta_K1, MinTheta_K2, &
    !$OMP          V_u_1, V_u_2, V_u_3, W2_K, W3_K )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG COLLAPSE(4) ASYNC &
    !$ACC COPYIN( dZ1, dZ2, dZ3, dZ4 ) &
    !$ACC PRIVATE( ResidualE, iK1, iK2, N_K1, N_K2, E_K1, E_K2, &
    !$ACC          Theta_K1, Theta_K2, MinTheta_K1, MinTheta_K2, &
    !$ACC          V_u_1, V_u_2, V_u_3, W2_K, W3_K ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, GE, GX, U_F, U_R, DeltaE, &
    !$ACC          Weights_Q, ApplyEnergyLimiter )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(4) &
    !$OMP PRIVATE( ResidualE, iK1, iK2, N_K1, N_K2, E_K1, E_K2, &
    !$OMP          Theta_K1, Theta_K2, MinTheta_K1, MinTheta_K2, &
    !$OMP          V_u_1, V_u_2, V_u_3, W2_K, W3_K, &
    !$OMP          iNodeX, iNodeE )
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      IF( ApplyEnergyLimiter(iZ2,iZ3,iZ4,iS) )THEN

        ! --- Three-Velocity ---

#if   defined( THORNADO_OMP_OL )
        !$OMP PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
        !$ACC LOOP VECTOR
#endif
        DO iNodeX = 1, nDOFX

          V_u_1(iNodeX) &
            = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S1) &
                / ( U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D) &
                      * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) )

          V_u_2(iNodeX) &
            = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S2) &
                / ( U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D) &
                      * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) )

          V_u_3(iNodeX) &
            = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S3) &
                / ( U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D) &
                      * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

        END DO

        ! --- Integration Weights ---

#if   defined( THORNADO_OMP_OL )
        !$OMP PARALLEL DO SIMD COLLAPSE(2) &
        !$OMP PRIVATE( iNodeX, iNodeE )
#elif defined( THORNADO_OACC   )
        !$ACC LOOP VECTOR COLLAPSE(2) &
        !$ACC PRIVATE( iNodeX, iNodeE )
#endif
        DO iZ1 = iZ_B0(1), iZ_E0(1)
          DO iNodeZ = 1, nDOFZ

            iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
            iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

            W2_K(iNodeZ,iZ1) &
              = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                  * W_Factor &
                  * GE(iNodeE,iZ1,iGE_Ep2) &
                  * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) &
                  * Weights_Q(iNodeZ)

            W3_K(iNodeZ,iZ1) &
              = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                  * W_Factor &
                  * GE(iNodeE,iZ1,iGE_Ep3) &
                  * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) &
                  * Weights_Q(iNodeZ)

          END DO
        END DO

        ResidualE = DeltaE(iZ_B0(1),iZ2,iZ3,iZ4,iS)

        DO iZ1 = iZ_B0(1), iZ_E0(1) - 1 ! --- Forward Sweep

          iK1 = iZ1
          iK2 = iZ1 + 1

          N_K1 &
            = ElementNumber &
                ( W2_K(:,iK1), &
                  U_R (:,iK1,iZ2,iZ3,iZ4,iCR_N,iS) )

          N_K2 &
            = ElementNumber &
                ( W2_K(:,iK2), &
                  U_R (:,iK2,iZ2,iZ3,iZ4,iCR_N,iS) )

          E_K1 &
            = ElementEnergy &
                ( W3_K (:,iK1), &
                  U_R  (:,iK1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                  U_R  (:,iK1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                  U_R  (:,iK1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                  U_R  (:,iK1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                  V_u_1, &
                  V_u_2, &
                  V_u_3 )

          E_K2 &
            = ElementEnergy &
                ( W3_K (:,iK2), &
                  U_R  (:,iK2,iZ2,iZ3,iZ4,iCR_N ,iS), &
                  U_R  (:,iK2,iZ2,iZ3,iZ4,iCR_G1,iS), &
                  U_R  (:,iK2,iZ2,iZ3,iZ4,iCR_G2,iS), &
                  U_R  (:,iK2,iZ2,iZ3,iZ4,iCR_G3,iS), &
                  V_u_1, &
                  V_u_2, &
                  V_u_3 )

          MinTheta_K1 = MIN(0.0_DP, Min_1 / &
                            MINVAL(U_R(:,iK1,iZ2,iZ3,iZ4,iCR_N ,iS)) - One)
          MinTheta_K2 = MIN(0.0_DP, Min_1 / &
                            MINVAL(U_R(:,iK2,iZ2,iZ3,iZ4,iCR_N ,iS)) - One)

          Theta_K1 = Zero
          Theta_K2 = Zero
          CALL UpdateResidualEnergy &
              ( N_K1, N_K2, E_K1, E_K2, DeltaE(iK2,iZ2,iZ3,iZ4,iS), &
                MinTheta_K1, MinTheta_K2, ResidualE, Theta_K1, Theta_K2 )

          CALL LimitEnergy( Theta_K1, U_R(:,iK1,iZ2,iZ3,iZ4,iCR_N ,iS) )
          CALL LimitEnergy( Theta_K2, U_R(:,iK2,iZ2,iZ3,iZ4,iCR_N ,iS) )
          CALL LimitEnergy( Theta_K1, U_R(:,iK1,iZ2,iZ3,iZ4,iCR_G1,iS) )
          CALL LimitEnergy( Theta_K2, U_R(:,iK2,iZ2,iZ3,iZ4,iCR_G1,iS) )
          CALL LimitEnergy( Theta_K1, U_R(:,iK1,iZ2,iZ3,iZ4,iCR_G2,iS) )
          CALL LimitEnergy( Theta_K2, U_R(:,iK2,iZ2,iZ3,iZ4,iCR_G2,iS) )
          CALL LimitEnergy( Theta_K1, U_R(:,iK1,iZ2,iZ3,iZ4,iCR_G3,iS) )
          CALL LimitEnergy( Theta_K2, U_R(:,iK2,iZ2,iZ3,iZ4,iCR_G3,iS) )

        END DO

        IF( ABS( ResidualE ) > Zero )THEN

          DO iZ1 = iZ_E0(1) - 2, iZ_B0(1), - 1 ! Backward Sweep

            iK1 = iZ1
            iK2 = iZ1 + 1

            N_K1 &
              = ElementNumber &
                  ( W2_K(:,iK1), &
                    U_R (:,iK1,iZ2,iZ3,iZ4,iCR_N,iS) )

            N_K2 &
              = ElementNumber &
                  ( W2_K(:,iK2), &
                    U_R (:,iK2,iZ2,iZ3,iZ4,iCR_N,iS) )

            E_K1 &
              = ElementEnergy &
                  ( W3_K (:,iK1), &
                    U_R  (:,iK1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                    U_R  (:,iK1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                    U_R  (:,iK1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                    U_R  (:,iK1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                    V_u_1, &
                    V_u_2, &
                    V_u_3 )

            E_K2 &
              = ElementEnergy &
                  ( W3_K (:,iK2), &
                    U_R  (:,iK2,iZ2,iZ3,iZ4,iCR_N ,iS), &
                    U_R  (:,iK2,iZ2,iZ3,iZ4,iCR_G1,iS), &
                    U_R  (:,iK2,iZ2,iZ3,iZ4,iCR_G2,iS), &
                    U_R  (:,iK2,iZ2,iZ3,iZ4,iCR_G3,iS), &
                    V_u_1, &
                    V_u_2, &
                    V_u_3 )

            MinTheta_K1 = MIN(0.0_DP, Min_1 / &
                              MINVAL(U_R(:,iK1,iZ2,iZ3,iZ4,iCR_N ,iS)) - One)
            MinTheta_K2 = MIN(0.0_DP, Min_1 / &
                              MINVAL(U_R(:,iK2,iZ2,iZ3,iZ4,iCR_N ,iS)) - One)

            Theta_K1 = Zero
            Theta_K2 = Zero
            CALL UpdateResidualEnergy &
                ( N_K1, N_K2, E_K1, E_K2, Zero, &
                MinTheta_K1, MinTheta_K2, ResidualE, Theta_K1, Theta_K2 )

            CALL LimitEnergy( Theta_K1, U_R(:,iK1,iZ2,iZ3,iZ4,iCR_N ,iS) )
            CALL LimitEnergy( Theta_K2, U_R(:,iK2,iZ2,iZ3,iZ4,iCR_N ,iS) )
            CALL LimitEnergy( Theta_K1, U_R(:,iK1,iZ2,iZ3,iZ4,iCR_G1,iS) )
            CALL LimitEnergy( Theta_K2, U_R(:,iK2,iZ2,iZ3,iZ4,iCR_G1,iS) )
            CALL LimitEnergy( Theta_K1, U_R(:,iK1,iZ2,iZ3,iZ4,iCR_G2,iS) )
            CALL LimitEnergy( Theta_K2, U_R(:,iK2,iZ2,iZ3,iZ4,iCR_G2,iS) )
            CALL LimitEnergy( Theta_K1, U_R(:,iK1,iZ2,iZ3,iZ4,iCR_G3,iS) )
            CALL LimitEnergy( Theta_K2, U_R(:,iK2,iZ2,iZ3,iZ4,iCR_G3,iS) )

          END DO

        END IF ! --- ABS( ResidualE ) > Zero

      END IF ! --- ANY( LimiterApplied )

    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

    CALL TimersStop( Timer_PL_EnergyLimiter )

  END SUBROUTINE ApplyEnergyLimiter_TwoMoment


  SUBROUTINE ComputePointValuesZ( iZ_B0, iZ_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFZ, &
          iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies)
    REAL(DP), INTENT(inout) :: &
      U_P(nPT_Z, &
          iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies)

    CALL TimersStart( Timer_PL_PointValues )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT_Z, N_R, nDOFZ, One, InterpMat_Z, nPT_Z, &
             U_Q, nDOFZ, Zero, U_P, nPT_Z )

    CALL TimersStop( Timer_PL_PointValues )

  END SUBROUTINE ComputePointValuesZ


  SUBROUTINE ComputePointValuesZ_Single( InterpMat_Z, U_Q, U_P )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE VECTOR
#endif

    REAL(DP), INTENT(in)  :: &
     InterpMat_Z(nPT_Z,nDOFZ)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFZ)
    REAL(DP), INTENT(out) :: &
      U_P(nPT_Z)

    REAL(DP) :: SUM1
    INTEGER  :: iNodeZ, iP_Z

#if   defined( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD &
    !$OMP PRIVATE( SUM1 )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR &
    !$ACC PRIVATE( SUM1 )
#endif
    DO iP_Z = 1, nPT_Z
      SUM1 = Zero
      DO iNodeZ = 1, nDOFZ
        SUM1 = SUM1 + InterpMat_Z(iP_Z,iNodeZ) * U_Q(iNodeZ)
      END DO
      U_P(iP_Z) = SUM1
    END DO

  END SUBROUTINE ComputePointValuesZ_Single


  SUBROUTINE ComputePointValuesX( iZ_B0, iZ_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFX, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(inout) :: &
      U_P(nPT_X, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4))

    CALL TimersStart( Timer_PL_PointValues )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT_X, N_G, nDOFX, One, InterpMat_X, nPT_X, &
             U_Q, nDOFX, Zero, U_P, nPT_X )

    CALL TimersStop( Timer_PL_PointValues )

  END SUBROUTINE ComputePointValuesX


  SUBROUTINE ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, U_Q, U_K )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: &
      Tau_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                  iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in)  :: &
      U_Q  (nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                  iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP), INTENT(out) :: &
      U_K  (      iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                  iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iNodeZ
    REAL(DP) :: SUM1, SUM2

    CALL TimersStart( Timer_PL_CellAverage )

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( SUM1, SUM2 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRIVATE( SUM1, SUM2 ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, U_K, Weights_Q, Tau_Q, U_Q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( SUM1, SUM2 )
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      SUM1 = Zero
      SUM2 = Zero
      DO iNodeZ = 1, nDOFZ
        SUM1 = SUM1 + Weights_Q(iNodeZ) * Tau_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4) * U_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4,iS)
        SUM2 = SUM2 + Weights_Q(iNodeZ) * Tau_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4)
      END DO
      U_K(iZ1,iZ2,iZ3,iZ4,iS) = SUM1 / SUM2

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop( Timer_PL_CellAverage )

  END SUBROUTINE ComputeCellAverage


  SUBROUTINE CheckCellAverageRealizability &
    ( iZ_B0, iZ_E0, N_K, G1_K, G2_K, G3_K, h_d_1, h_d_2, h_d_3, &
      RealizableCellAverage )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in) :: &
      N_K (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(in) :: &
      G1_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(in) :: &
      G2_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(in) :: &
      G3_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(in) :: &
      h_d_1(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      h_d_2(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      h_d_3(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    LOGICAL, INTENT(out) :: &
      RealizableCellAverage &
          (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP_X
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33, Gamma_K

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( Gm_dd_11, Gm_dd_22, Gm_dd_33, Gamma_K )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRIVATE( Gm_dd_11, Gm_dd_22, Gm_dd_33, Gamma_K ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, RealizableCellAverage, &
    !$ACC          h_d_1, h_d_2, h_d_3, N_K, G1_K, G2_K, G3_K )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( Gm_dd_11, Gm_dd_22, Gm_dd_33, Gamma_K )
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

      ! --- Check for Negative Density ---

      IF( N_K(iZ1,iZ2,iZ3,iZ4,iS) < Min_1 )THEN

        RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) = .FALSE.

      ELSE

        ! --- Check for Negative "Gamma" ---

        DO iP_X = 1, nPT_X

          Gm_dd_11 = MAX( h_d_1(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )
          Gm_dd_22 = MAX( h_d_2(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )
          Gm_dd_33 = MAX( h_d_3(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )

          Gamma_K &
            = GammaFun &
                ( N_K (iZ1,iZ2,iZ3,iZ4,iS), &
                  G1_K(iZ1,iZ2,iZ3,iZ4,iS), &
                  G2_K(iZ1,iZ2,iZ3,iZ4,iS), &
                  G3_K(iZ1,iZ2,iZ3,iZ4,iS), &
                  Gm_dd_11, Gm_dd_22, Gm_dd_33 )

          IF(  Gamma_K < Min_2 )THEN

            RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) = .FALSE.

          END IF

        END DO

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE CheckCellAverageRealizability


  SUBROUTINE RecoverRealizableCellAverage &
    ( iZ_B0, iZ_E0, N_K, G1_K, G2_K, G3_K, N_Q, G1_Q, G2_Q, G3_Q, &
      h_d_1, h_d_2, h_d_3, RealizableCellAverage )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)    :: &
      N_K (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(in)    :: &
      G1_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(in)    :: &
      G2_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(in)    :: &
      G3_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(inout) :: &
      N_Q (nDOFZ, &
           iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(inout) :: &
      G1_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(inout) :: &
      G2_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(inout) :: &
      G3_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP), INTENT(in) :: &
      h_d_1(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      h_d_2(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      h_d_3(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    LOGICAL,  INTENT(in) :: &
      RealizableCellAverage &
          (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP_X
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33, absG_K

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
    !$OMP PRIVATE( absG_K, Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#elif defined( THORNADO_OACC   )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) ASYNC &
    !$ACC PRIVATE( absG_K, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, RealizableCellAverage, &
    !$ACC          h_d_1, h_d_2, h_d_3, &
    !$ACC          N_K, G1_K, G2_K, G3_K, N_Q, G1_Q, G2_Q, G3_Q )
#elif defined( THORNADO_OMP    )
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( absG_K, Gm_dd_11, Gm_dd_22, Gm_dd_33 )
#endif
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      IF( .NOT. RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) )THEN

        IF( N_K(iZ1,iZ2,iZ3,iZ4,iS) < Min_1 )THEN

          N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) = 1.1_DP * Min_1
          G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = Zero
          G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = Zero
          G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = Zero

        ELSE

          absG_K = Zero
          DO iP_X = 1, nPT_X

            Gm_dd_11 = MAX( h_d_1(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )
            Gm_dd_22 = MAX( h_d_2(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )
            Gm_dd_33 = MAX( h_d_3(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )

            absG_K &
              = MAX( absG_K,   G1_K(iZ1,iZ2,iZ3,iZ4,iS)**2 / Gm_dd_11 &
                             + G2_K(iZ1,iZ2,iZ3,iZ4,iS)**2 / Gm_dd_22 &
                             + G3_K(iZ1,iZ2,iZ3,iZ4,iS)**2 / Gm_dd_33 )

          END DO

          absG_K = MAX( SQRT( absG_K ), SqrtTiny )

          N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) &
            = N_K(iZ1,iZ2,iZ3,iZ4,iS)

          G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            = ( G1_K(iZ1,iZ2,iZ3,iZ4,iS) / absG_K ) &
                * 0.99_DP * N_K(iZ1,iZ2,iZ3,iZ4,iS)

          G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            = ( G2_K(iZ1,iZ2,iZ3,iZ4,iS) / absG_K ) &
                * 0.99_DP * N_K(iZ1,iZ2,iZ3,iZ4,iS)

          G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            = ( G3_K(iZ1,iZ2,iZ3,iZ4,iS) / absG_K ) &
                * 0.99_DP * N_K(iZ1,iZ2,iZ3,iZ4,iS)

        END IF

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE RecoverRealizableCellAverage


  REAL(DP) FUNCTION GammaFun( N, G1, G2, G3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: N, G1, G2, G3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    GammaFun &
      = N - SQRT( G1**2 / Gm_dd_11 + G2**2 / Gm_dd_22 + G3**2 / Gm_dd_33 )

    RETURN
  END FUNCTION GammaFun


  SUBROUTINE SolveTheta_Bisection &
    ( N_P, G1_P, G2_P, G3_P, N_K, G1_K, G2_K, G3_K, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Theta )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: N_P, G1_P, G2_P, G3_P
    REAL(DP), INTENT(in)  :: N_K, G1_K, G2_K, G3_K
    REAL(DP), INTENT(in)  :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP), INTENT(out) :: Theta

    INTEGER,  PARAMETER :: ITERATION_MAX = 12
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = GammaFun &
            ( N_K, G1_K, G2_K, G3_K, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) - Min_2

    x_b = One
    f_b = GammaFun &
            ( N_P, G1_P, G2_P, G3_P, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) - Min_2

    dx = One

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED .AND. ITERATION < ITERATION_MAX )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx

      f_c = GammaFun &
              ( x_c * N_P  + ( One - x_c ) * N_K,  &
                x_c * G1_P + ( One - x_c ) * G1_K, &
                x_c * G2_P + ( One - x_c ) * G2_K, &
                x_c * G3_P + ( One - x_c ) * G3_K, &
                Gm_dd_11, Gm_dd_22, Gm_dd_33 ) - Min_2

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx < dx_min ) CONVERGED = .TRUE.

    END DO

    IF( ITERATION >= ITERATION_MAX )THEN
      Theta = Zero
    ELSE
      Theta = x_a
    END IF

  END SUBROUTINE SolveTheta_Bisection


  REAL(DP) FUNCTION ElementNumber( W, N )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE VECTOR
#endif

    REAL(DP), INTENT(in) :: W(nDOFZ)
    REAL(DP), INTENT(in) :: N(nDOFZ)

    INTEGER :: iNodeZ

    ElementNumber = Zero

#if defined  ( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD &
    !$OMP REDUCTION( +: ElementNumber )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR &
    !$ACC REDUCTION( +: ElementNumber )
#endif
    DO iNodeZ = 1, nDOFZ

      ElementNumber = ElementNumber + W(iNodeZ) * N(iNodeZ)

    END DO

    RETURN
  END FUNCTION ElementNumber


  REAL(DP) FUNCTION ElementEnergy &
    ( W, N, G_d_1, G_d_2, G_d_3, V_u_1, V_u_2, V_u_3 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE VECTOR
#endif

    REAL(DP), INTENT(in) :: W    (nDOFZ)
    REAL(DP), INTENT(in) :: N    (nDOFZ)
    REAL(DP), INTENT(in) :: G_d_1(nDOFZ)
    REAL(DP), INTENT(in) :: G_d_2(nDOFZ)
    REAL(DP), INTENT(in) :: G_d_3(nDOFZ)
    REAL(DP), INTENT(in) :: V_u_1(nDOFX)
    REAL(DP), INTENT(in) :: V_u_2(nDOFX)
    REAL(DP), INTENT(in) :: V_u_3(nDOFX)

    INTEGER :: iNodeX, iNodeZ

    ElementEnergy = Zero

#if defined  ( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD &
    !$OMP PRIVATE( iNodeX ) &
    !$OMP REDUCTION( +: ElementEnergy )
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR &
    !$ACC PRIVATE( iNodeX ) &
    !$ACC REDUCTION( +: ElementEnergy )
#endif
    DO iNodeZ = 1, nDOFZ

      iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

      ElementEnergy &
        = ElementEnergy &
            + W(iNodeZ) &
                * ( N(iNodeZ) &
                    + V_u_1(iNodeX) * G_d_1(iNodeZ) &
                    + V_u_2(iNodeX) * G_d_2(iNodeZ) &
                    + V_u_3(iNodeX) * G_d_3(iNodeZ) )

    END DO

    RETURN
  END FUNCTION ElementEnergy


  SUBROUTINE UpdateResidualEnergy &
      ( N_K1, N_K2, E_K1, E_K2, DeltaE, MinTheta_K1, MinTheta_K2, &
      ResidualE, Theta_K1, Theta_K2 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: N_K1, N_K2, E_K1, E_K2, DeltaE
    REAL(DP), INTENT(in)    :: MinTheta_K1, MinTheta_K2
    REAL(DP), INTENT(inout) :: ResidualE
    REAL(DP), INTENT(out)   :: Theta_K1, Theta_K2

    REAL(DP), PARAMETER :: MinTheta_K = - 0.5_DP
    REAL(DP), PARAMETER :: MaxTheta_K =   1.0_DP

    REAL(DP) :: Det
    REAL(DP) :: MinTheta_K1_Both, MinTheta_K2_Both
    REAL(DP) :: RedTheta_K


    Det = ( N_K1 * E_K2 - N_K2 * E_K1 )

    IF( ABS( Det ) .GT. SqrtTiny )THEN

      Theta_K1 =   N_K2 * ( ResidualE + DeltaE ) / Det
      Theta_K2 = - N_K1 * ( ResidualE + DeltaE ) / Det

    ELSE

      Theta_K1 = Zero
      Theta_K2 = Zero

    END IF

    MinTheta_K1_Both = MAX(MinTheta_K1, MinTheta_K)
    MinTheta_K2_Both = MAX(MinTheta_K2, MinTheta_K)

    RedTheta_K = 1.0_DP

    IF( Theta_K1 < MinTheta_K1_Both )THEN

      RedTheta_K = MinTheta_K1_Both/Theta_K1

    END IF

    IF( Theta_K2 < MinTheta_K2_Both )THEN

      RedTheta_K &
        = MIN( MinTheta_K2_Both/Theta_K2, RedTheta_K )

    END IF

    IF( Theta_K1 > MaxTheta_K )THEN

      RedTheta_K &
        = MIN( MaxTheta_K/Theta_K1, RedTheta_K )

    END IF

    IF( Theta_K2 > MaxTheta_K )THEN

      RedTheta_K &
        = MIN( MaxTheta_K/Theta_K2, RedTheta_K )

    END IF

    Theta_K1 = RedTheta_K * Theta_K1
    Theta_K2 = RedTheta_K * Theta_K2

    ResidualE &
      = E_K1 * Theta_K1 + E_K2 * Theta_K2 &
          + ( ResidualE + DeltaE )

    RETURN
  END SUBROUTINE UpdateResidualEnergy


  SUBROUTINE LimitEnergy( Theta_K, U )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE VECTOR
#endif

    REAL(DP), INTENT(in)    :: Theta_K
    REAL(DP), INTENT(inout) :: U(nDOFZ)

    INTEGER :: iNodeZ

#if defined  ( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR
#endif
    DO iNodeZ = 1, nDOFZ
      U(iNodeZ) = ( One + Theta_K ) * U(iNodeZ)
    END DO

    RETURN
  END SUBROUTINE LimitEnergy


  SUBROUTINE ApplyLimiter( Theta, U_K, U_Q )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE VECTOR
#endif

    REAL(DP), INTENT(in)    :: Theta, U_K
    REAL(DP), INTENT(inout) :: U_Q(nDOFZ)

    INTEGER :: iNodeZ

#if defined  ( THORNADO_OMP_OL )
    !$OMP PARALLEL DO SIMD
#elif defined( THORNADO_OACC   )
    !$ACC LOOP VECTOR
#endif
    DO iNodeZ = 1, nDOFZ
      U_Q(iNodeZ) = Theta * U_Q(iNodeZ) + ( One - Theta ) * U_K
    END DO

    RETURN
  END SUBROUTINE ApplyLimiter


  SUBROUTINE ComputeLocalEnergy &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, ApplyEnergyLimiter, Energy )

    USE MeshModule, ONLY: &
      MeshE, MeshX

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GE (1:nDOFE, &
          iZ_B1(1):iZ_E1(1), &
          1:nGE)
    REAL(DP), INTENT(in)  :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in)  :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(in)  :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR, &
          1:nSpecies)
    LOGICAL, INTENT(in)   :: &
      ApplyEnergyLimiter &
        (iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies)
    REAL(DP), INTENT(out) :: &
      Energy &
        (iZ_B0(1):iZ_E0(1), &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         1:nSpecies)

    INTEGER  :: iNodeE, iNodeX, iNodeZ
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: V_u_1, V_u_2, V_u_3, SUM_E

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(5) &
    !$OMP MAP( to: dZ1, dZ2, dZ3, dZ4 ) &
    !$OMP PRIVATE( SUM_E )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG COLLAPSE(5) ASYNC &
    !$ACC COPYIN( dZ1, dZ2, dZ3, dZ4 ) &
    !$ACC PRIVATE( SUM_E ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, dZ1, dZ2, dZ3, dZ4, ApplyEnergyLimiter, &
    !$ACC          Weights_Q, GE, GX, U_F, U_R, Energy )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(5) &
    !$OMP PRIVATE( iNodeZ, V_u_1, V_u_2, V_u_3, SUM_E )
#endif
    DO iS  = 1       , nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      IF( ApplyEnergyLimiter(iZ2,iZ3,iZ4,iS) )THEN

      SUM_E = Zero

#if defined(THORNADO_OMP_OL)
      !$OMP PARALLEL DO SIMD COLLAPSE(2) &
      !$OMP PRIVATE( iNodeZ, V_u_1, V_u_2, V_u_3 ) &
      !$OMP REDUCTION( +: SUM_E )
#elif defined(THORNADO_OACC)
      !$ACC LOOP VECTOR COLLAPSE(2) &
      !$ACC PRIVATE( iNodeZ, V_u_1, V_u_2, V_u_3 ) &
      !$ACC REDUCTION( +: SUM_E )
#endif
      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        V_u_1 &
          = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S1) &
              / ( U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D) &
                    * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) )

        V_u_2 &
          = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S2) &
              / ( U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D) &
                    * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) )

        V_u_3 &
          = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S3) &
              / ( U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D) &
                    * GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

        SUM_E &
          = SUM_E &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                * W_Factor &
                * Weights_Q(iNodeZ) &
                * GE(iNodeE,iZ1,iGE_Ep3) &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm) &
                * (   U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) &
                    + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) * V_u_1 &
                    + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) * V_u_2 &
                    + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) * V_u_3 )

      END DO
      END DO

      Energy(iZ1,iZ2,iZ3,iZ4,iS) = SUM_E

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

  END SUBROUTINE ComputeLocalEnergy


  SUBROUTINE ComputeGlobalEnergyMomentum &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, EnergyMomentum )

    USE MeshModule, ONLY: &
      MeshE, MeshX

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE (1:nDOFE, &
          iZ_B1(1):iZ_E1(1), &
          1:nGE)
    REAL(DP), INTENT(in) :: &
      GX (1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nGF)
    REAL(DP), INTENT(in) :: &
      U_F(1:nDOFX, &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCF)
    REAL(DP), INTENT(in) :: &
      U_R(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCR, &
          1:nSpecies)
    REAL(DP), INTENT(out) :: &
      EnergyMomentum(nCR)

    INTEGER  :: iNodeE, iNodeX, iNodeZ
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: V_u_1, V_u_2, V_u_3
    REAL(DP) :: W3_K, SUM_N, SUM_G1, SUM_G2, SUM_G3

    SUM_N  = Zero
    SUM_G1 = Zero
    SUM_G2 = Zero
    SUM_G3 = Zero

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(7) &
    !$OMP MAP( to: dZ1, dZ2, dZ3, dZ4 ) &
    !$OMP MAP( tofrom: SUM_N, SUM_G1, SUM_G2, SUM_G3 ) &
    !$OMP PRIVATE( iNodeZ, V_d_1, V_d_2, V_d_3, V_u_1, V_u_2, V_u_3, W3_K ) &
    !$OMP REDUCTION( +: SUM_N, SUM_G1, SUM_G2, SUM_G3 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(7) ASYNC &
    !$ACC COPYIN( dZ1, dZ2, dZ3, dZ4 ) &
    !$ACC COPY( SUM_N, SUM_G1, SUM_G2, SUM_G3 ) &
    !$ACC PRIVATE( iNodeZ, V_d_1, V_d_2, V_d_3, V_u_1, V_u_2, V_u_3, W3_K ) &
    !$ACC REDUCTION( +: SUM_N, SUM_G1, SUM_G2, SUM_G3 ) &
    !$ACC PRESENT( iZ_B0, iZ_E0, Weights_Q, GE, GX, U_F, U_R )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(7) &
    !$OMP PRIVATE( iNodeZ, V_d_1, V_d_2, V_d_3, V_u_1, V_u_2, V_u_3, W3_K ) &
    !$OMP REDUCTION( +: SUM_N, SUM_G1, SUM_G2, SUM_G3 )
#endif
    DO iS  = 1       , nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        V_d_1 = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S1) / U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D)
        V_d_2 = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S2) / U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D)
        V_d_3 = U_F(iNodeX,iZ2,iZ3,iZ4,iCF_S3) / U_F(iNodeX,iZ2,iZ3,iZ4,iCF_D)

        V_u_1 = V_d_1 / GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)
        V_u_2 = V_d_2 / GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)
        V_u_3 = V_d_3 / GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)

        W3_K = dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
               * W_Factor &
               * Weights_Q(iNodeZ) &
               * GE(iNodeE,iZ1,iGE_Ep3) &
               * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

        SUM_N &
          = SUM_N  + W3_K * (   U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) &
                              + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) * V_u_1 &
                              + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) * V_u_2 &
                              + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) * V_u_3 )

        SUM_G1 &
          = SUM_G1 + W3_K * (   U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) &
                              + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) * V_d_1 )

        SUM_G2 &
          = SUM_G2 + W3_K * (   U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) &
                              + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) * V_d_2 )

        SUM_G3 &
          = SUM_G3 + W3_K * (   U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) &
                              + U_R(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS) * V_d_3 )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
#elif defined(THORNADO_OACC)
    !$ACC WAIT
#endif
    EnergyMomentum(iCR_N ) = SUM_N
    EnergyMomentum(iCR_G1) = SUM_G1
    EnergyMomentum(iCR_G2) = SUM_G2
    EnergyMomentum(iCR_G3) = SUM_G3

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeGlobalEnergyMomentum


END MODULE TwoMoment_PositivityLimiterModule
