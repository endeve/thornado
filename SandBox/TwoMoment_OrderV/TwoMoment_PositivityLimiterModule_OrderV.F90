MODULE TwoMoment_PositivityLimiterModule_OrderV

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nNodesZ, &
    nDOFE, nNodesE, &
    nDOFX, nNodesX
  USE TwoMoment_TimersModule_OrderV, ONLY: &
    TimersStart, &
    TimersStop, &
    Timer_PositivityLimiter
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
    nGE, iGE_Ep2
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm, iGF_h_1, iGF_h_2, iGF_h_3
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_TwoMoment
  PUBLIC :: FinalizePositivityLimiter_TwoMoment
  PUBLIC :: ApplyPositivityLimiter_TwoMoment

  LOGICAL               :: UsePositivityLimiter
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
  REAL(DP),   PARAMETER :: One_EPS = One - 1.0d1 * EPSILON( One )
  REAL(DP), ALLOCATABLE :: InterpMat_Z(:,:)
  REAL(DP), ALLOCATABLE :: InterpMat_X(:,:)

CONTAINS


  SUBROUTINE InitializePositivityLimiter_TwoMoment &
    ( Min_1_Option, Max_1_Option, Min_2_Option, UsePositivityLimiter_Option, &
      Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Max_1_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option
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
      WRITE(*,*)
      WRITE(*,'(A4,A32,ES11.3E3)') '', 'Min_1: ', Min_1
      WRITE(*,'(A4,A32,ES11.3E3)') '', 'Max_1: ', Max_1
      WRITE(*,'(A4,A32,ES11.3E3)') '', 'Min_2: ', Min_2

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

  END SUBROUTINE InitializePositivityLimiter_TwoMoment


  SUBROUTINE FinalizePositivityLimiter_TwoMoment

    DEALLOCATE( InterpMat_Z )
    DEALLOCATE( InterpMat_X )
    DEALLOCATE( PointZ2X )

  END SUBROUTINE FinalizePositivityLimiter_TwoMoment


  SUBROUTINE ApplyPositivityLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R )

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

    LOGICAL  :: RecomputePointValues
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP_Z, iP_X
    INTEGER  :: iNodeZ, iNodeE, iNodeX
    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3)
    REAL(DP) :: Min_K, Max_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: Gamma, Gamma_Min
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: &
      Tau_Q(nDOFZ, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      N_Q(nDOFZ, &
          iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies), &
      N_P(nPT_Z, &
          iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies), &
      N_K(iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies)
    REAL(DP) :: &
      G1_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      G1_P(nPT_Z, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      G1_K(iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP) :: &
      G2_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      G2_P(nPT_Z, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      G2_K(iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP) :: &
      G3_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      G3_P(nPT_Z, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      G3_K(iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP) :: &
      h_d_1_Q(nDOFX, &
              iZ_B0(2):iZ_E0(2), &
              iZ_B0(3):iZ_E0(3), &
              iZ_B0(4):iZ_E0(4)), &
      h_d_1_P(nPT_X, &
              iZ_B0(2):iZ_E0(2), &
              iZ_B0(3):iZ_E0(3), &
              iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      h_d_2_Q(nDOFX, &
              iZ_B0(2):iZ_E0(2), &
              iZ_B0(3):iZ_E0(3), &
              iZ_B0(4):iZ_E0(4)), &
      h_d_2_P(nPT_X, &
              iZ_B0(2):iZ_E0(2), &
              iZ_B0(3):iZ_E0(3), &
              iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      h_d_3_Q(nDOFX, &
              iZ_B0(2):iZ_E0(2), &
              iZ_B0(3):iZ_E0(3), &
              iZ_B0(4):iZ_E0(4)), &
      h_d_3_P(nPT_X, &
              iZ_B0(2):iZ_E0(2), &
              iZ_B0(3):iZ_E0(3), &
              iZ_B0(4):iZ_E0(4))

    IF( .NOT. UsePositivityLimiter .OR. nDOFZ == 1 ) RETURN

    CALL TimersStart( Timer_PositivityLimiter )

    N_R = nSpecies * PRODUCT( iZ_E0 - iZ_B0 + 1 )

    iX_B0 = iZ_B0(2:4)
    iX_E0 = iZ_E0(2:4)

    N_G = PRODUCT( iX_E0 - iX_B0 + 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        h_d_1_Q(iNodeX,iX1,iX2,iX3) = GX(iNodeX,iX1,iX2,iX3,iGF_h_1)
        h_d_2_Q(iNodeX,iX1,iX2,iX3) = GX(iNodeX,iX1,iX2,iX3,iGF_h_2)
        h_d_3_Q(iNodeX,iX1,iX2,iX3) = GX(iNodeX,iX1,iX2,iX3,iGF_h_3)

      END DO

    END DO
    END DO
    END DO

    CALL ComputePointValuesX( iX_B0, iX_E0, h_d_1_Q, h_d_1_P )
    CALL ComputePointValuesX( iX_B0, iX_E0, h_d_2_Q, h_d_2_P )
    CALL ComputePointValuesX( iX_B0, iX_E0, h_d_3_Q, h_d_3_P )

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

    CALL ComputePointValuesZ( iZ_B0, iZ_E0, N_Q , N_P  )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G1_Q, G1_P )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G2_Q, G2_P )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G3_Q, G3_P )

    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, N_Q , N_K  )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G1_Q, G1_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G2_Q, G2_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G3_Q, G3_K )

    ! --- Ensure Bounded Density ---

    RecomputePointValues = .FALSE.

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      Min_K = Min_1
      Max_K = Max_1

      DO iP_Z = 1, nPT_Z

        Min_K = MIN( Min_K, N_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )
        Max_K = MAX( Max_K, N_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )

      END DO

      IF( Min_K < Min_1 .OR. Max_K > Max_1 )THEN

        Theta_1 &
          = MIN( One, &
                 ABS( ( Min_1 - N_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                      / ( Min_K - N_K(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny ) ), &
                 ABS( ( Max_1 - N_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                      / ( Max_K - N_K(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny ) ) )

        Theta_1 = One_EPS * Theta_1

        N_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_1 * N_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_1 ) * N_K(iZ1,iZ2,iZ3,iZ4,iS)

        RecomputePointValues = .TRUE.

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    IF( RecomputePointValues )THEN

      CALL ComputePointValuesZ( iZ_B0, iZ_E0, N_Q , N_P  )

    END IF

    ! --- Ensure Positive "Gamma" ---

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      Theta_2   = One
      Gamma_Min = Min_2

      DO iP_Z = 1, nPT_Z

        iP_X = PointZ2X(iP_Z)

        Gm_dd_11 = MAX( h_d_1_P(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )
        Gm_dd_22 = MAX( h_d_2_P(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )
        Gm_dd_33 = MAX( h_d_3_P(iP_X,iZ2,iZ3,iZ4)**2, SqrtTiny )

        Gamma &
          = GammaFun &
              ( N_P (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                G1_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                G2_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                G3_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        Gamma_Min = MIN( Gamma, Gamma_Min )

        IF( Gamma < Min_2 )THEN

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

        N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_2 * N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_2 ) * N_K (iZ1,iZ2,iZ3,iZ4,iS)

        G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_2 * G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_2 ) * G1_K(iZ1,iZ2,iZ3,iZ4,iS)

        G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_2 * G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_2 ) * G2_K(iZ1,iZ2,iZ3,iZ4,iS)

        G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
          = Theta_2 * G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            + ( One - Theta_2 ) * G3_K(iZ1,iZ2,iZ3,iZ4,iS)

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iS = 1, nSpecies
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

    CALL TimersStop( Timer_PositivityLimiter )

  END SUBROUTINE ApplyPositivityLimiter_TwoMoment


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
    REAL(DP), INTENT(out) :: &
      U_P(nPT_Z, &
          iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies)

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT_Z, N_R, nDOFZ, One, InterpMat_Z, nPT_Z, &
             U_Q, nDOFZ, Zero, U_P, nPT_Z )

  END SUBROUTINE ComputePointValuesZ


  SUBROUTINE ComputePointValuesX( iX_B0, iX_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFX, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: &
      U_P(nPT_X, &
          iX_B0(1):iX_E0(1), &
          iX_B0(2):iX_E0(2), &
          iX_B0(3):iX_E0(3))

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT_X, N_G, nDOFX, One, InterpMat_X, nPT_X, &
             U_Q, nDOFX, Zero, U_P, nPT_X )

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

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      U_K(iZ1,iZ2,iZ3,iZ4,iS) &
        = SUM( Weights_Q(:) * Tau_Q(:,iZ1,iZ2,iZ3,iZ4) &
                 * U_Q(:,iZ1,iZ2,iZ3,iZ4,iS) ) &
          / SUM( Weights_Q(:) * Tau_Q(:,iZ1,iZ2,iZ3,iZ4) )

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeCellAverage


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


END MODULE TwoMoment_PositivityLimiterModule_OrderV
