MODULE TwoMoment_PositivityLimiterModule_FMC

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny, FourPi
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nNodesZ, &
    nDOFE, nNodesE, &
    nDOFX, nNodesX
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
    nPF!, iCF_D, iCF_S1, iCF_S2, iCF_S3
  USE TwoMoment_FieldsModule_FMC, ONLY: &
    nSpecies, &
    nCM, iCM_E, iCM_F1, iCM_F2, iCM_F3

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

  END SUBROUTINE InitializePositivityLimiter_TwoMoment

  SUBROUTINE FinalizePositivityLimiter_TwoMoment

    DEALLOCATE( InterpMat_Z )
    DEALLOCATE( InterpMat_X )
    DEALLOCATE( PointZ2X )

  END SUBROUTINE FinalizePositivityLimiter_TwoMoment

  SUBROUTINE ApplyPositivityLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_M )

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
          1:nPF)
    REAL(DP), INTENT(inout) :: &
      U_M(1:nDOFZ, &
          iZ_B1(1):iZ_E1(1), &
          iZ_B1(2):iZ_E1(2), &
          iZ_B1(3):iZ_E1(3), &
          iZ_B1(4):iZ_E1(4), &
          1:nCM, &
          1:nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP_Z, iP_X
    INTEGER  :: iNodeZ, iNodeE, iNodeX
    REAL(DP) :: Min_K, Max_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: Gamma_P, Gamma_Min
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    ! REAL(DP) :: EnergyMomentum_0(nCM)
    ! REAL(DP) :: EnergyMomentum_1(nCM)
    LOGICAL  :: &
      RealizableCellAverage &
        (iZ_B0(1):iZ_E0(1), &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies)
    LOGICAL  :: &
      LimiterApplied &
        (iZ_B0(1):iZ_E0(1), &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies)
    LOGICAL  :: &
      ApplyEnergyLimiter &
        (iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies)
    REAL(DP) :: &
      Theta_1_K &
        (iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4)), &
      Theta_2_K &
        (iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      Energy_K &
        (iZ_B0(1):iZ_E0(1), &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies)
    REAL(DP) :: &
      Tau_Q(nDOFZ, &
            iZ_B0(1):iZ_E0(1), &
            iZ_B0(2):iZ_E0(2), &
            iZ_B0(3):iZ_E0(3), &
            iZ_B0(4):iZ_E0(4))
    REAL(DP) :: &
      E_Q(nDOFZ, &
          iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies), &
      E_P(nPT_Z, &
          iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies), &
      E_K(iZ_B0(1):iZ_E0(1), &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4), &
          nSpecies)
    REAL(DP) :: &
      F1_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      F1_P(nPT_Z, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      F1_K(iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP) :: &
      F2_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      F2_P(nPT_Z, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      F2_K(iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP) :: &
      F3_Q(nDOFZ, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      F3_P(nPT_Z, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           nSpecies), &
      F3_K(iZ_B0(1):iZ_E0(1), &
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

    Write(*,*)
    print *,'ApplyPositivityLimiter_TwoMoment'

    N_R = nSpecies * PRODUCT( iZ_E0 - iZ_B0 + 1 )

    N_G = PRODUCT( iZ_E0(2:4) - iZ_B0(2:4) + 1 )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      Theta_1_K(iZ2,iZ3,iZ4) = One
      Theta_2_K(iZ2,iZ3,iZ4) = One

    END DO
    END DO
    END DO

    ! CALL ComputeGlobalEnergyMomentum & ! --- For Global Tally ---
    !        ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, EnergyMomentum_0 )

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

    CALL ComputePointValuesX( iZ_B0, iZ_E0, h_d_1_Q, h_d_1_P )
    CALL ComputePointValuesX( iZ_B0, iZ_E0, h_d_2_Q, h_d_2_P )
    CALL ComputePointValuesX( iZ_B0, iZ_E0, h_d_3_Q, h_d_3_P )

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

        E_Q (iNodeZ,iZ1,iZ2,iZ3,iZ4,iS) &
          = U_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS)
        F1_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4,iS) &
          = U_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS)
        F2_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4,iS) &
          = U_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS)
        F3_Q(iNodeZ,iZ1,iZ2,iZ3,iZ4,iS) &
          = U_M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS)

     END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL ComputePointValuesZ( iZ_B0, iZ_E0, E_Q , E_P  )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, F1_Q, F1_P )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, F2_Q, F2_P )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, F3_Q, F3_P )

    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, E_Q , E_K  )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, F1_Q, F1_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, F2_Q, F2_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, F3_Q, F3_K )

    ! --- Ensure Bounded Density ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      LimiterApplied(iZ1,iZ2,iZ3,iZ4,iS) = .FALSE.

      Min_K = Min_1
      Max_K = Max_1

      DO iP_Z = 1, nPT_Z

        Min_K = MIN( Min_K, E_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )
        Max_K = MAX( Max_K, E_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS) )

      END DO

      IF( Min_K < Min_1 .OR. Max_K > Max_1 )THEN

        Theta_1 &
          = MIN( One, &
                  ABS( ( Min_1 - E_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                    / ( Min_K - E_K(iZ1,iZ2,iZ3,iZ4,iS)-SqrtTiny ) ), &
                  ABS( ( Max_1 - E_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                    / ( Max_K - E_K(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny ) ) )

        Theta_1 = One_EPS * Theta_1

        CALL ApplyLimiter &
                ( Theta_1, E_K(iZ1,iZ2,iZ3,iZ4,iS), &
                  E_Q(:,iZ1,iZ2,iZ3,iZ4,iS) )

        CALL ComputePointValuesZ_Single &
                ( InterpMat_Z, &
                  E_Q(:,iZ1,iZ2,iZ3,iZ4,iS), &
                  E_P(:,iZ1,iZ2,iZ3,iZ4,iS) )

        LimiterApplied(iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

        Theta_1_K(iZ2,iZ3,iZ4) = MIN( Theta_1, Theta_1_K(iZ2,iZ3,iZ4) )

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Ensure Positive "Gamma" ---

    DO iS  = 1, nSpecies
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

          Gamma_P &
            = GammaFun &
                ( E_P (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                  F1_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                  F2_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                  F3_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                  Gm_dd_11, Gm_dd_22, Gm_dd_33 )

          Gamma_Min = MIN( Gamma_P, Gamma_Min )

          IF( Gamma_P < Min_2 )THEN

            CALL SolveTheta_Bisection &
                   ( E_P (iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                     F1_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                     F2_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                     F3_P(iP_Z,iZ1,iZ2,iZ3,iZ4,iS), &
                     E_K (     iZ1,iZ2,iZ3,iZ4,iS), &
                     F1_K(     iZ1,iZ2,iZ3,iZ4,iS), &
                     F2_K(     iZ1,iZ2,iZ3,iZ4,iS), &
                     F3_K(     iZ1,iZ2,iZ3,iZ4,iS), &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33 , &
                     Theta_P )

            Theta_2 = MIN( Theta_2, Theta_P )

          END IF

        END DO

        IF( Gamma_Min < Min_2 )THEN

          ! --- Limit Towards Cell Average ---

          Theta_2 = One_EPS * Theta_2

          CALL ApplyLimiter &
                 ( Theta_2, E_K (iZ1,iZ2,iZ3,iZ4,iS), &
                   E_Q (:,iZ1,iZ2,iZ3,iZ4,iS) )
          CALL ApplyLimiter &
                 ( Theta_2, F1_K(iZ1,iZ2,iZ3,iZ4,iS), &
                   F1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) )
          CALL ApplyLimiter &
                 ( Theta_2, F2_K(iZ1,iZ2,iZ3,iZ4,iS), &
                   F2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) )
          CALL ApplyLimiter &
                 ( Theta_2, F3_K(iZ1,iZ2,iZ3,iZ4,iS), &
                   F3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) )

          LimiterApplied(iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

          Theta_2_K(iZ2,iZ3,iZ4) = MIN( Theta_2, Theta_2_K(iZ2,iZ3,iZ4) )

        END IF

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      U_M(:,iZ1,iZ2,iZ3,iZ4,iCM_E ,iS) = E_Q (:,iZ1,iZ2,iZ3,iZ4,iS)
      U_M(:,iZ1,iZ2,iZ3,iZ4,iCM_F1,iS) = F1_Q(:,iZ1,iZ2,iZ3,iZ4,iS)
      U_M(:,iZ1,iZ2,iZ3,iZ4,iCM_F2,iS) = F2_Q(:,iZ1,iZ2,iZ3,iZ4,iS)
      U_M(:,iZ1,iZ2,iZ3,iZ4,iCM_F3,iS) = F3_Q(:,iZ1,iZ2,iZ3,iZ4,iS)

    END DO
    END DO
    END DO
    END DO
    END DO

    ! IF( PRESENT( uDR_Option ) )THEN

    !   DO iZ4 = iZ_B0(4), iZ_E0(4)
    !   DO iZ3 = iZ_B0(3), iZ_E0(3)
    !   DO iZ2 = iZ_B0(2), iZ_E0(2)

    !     uDR_Option(iZ2,iZ3,iZ4,iDR_PL_Theta_1) = Theta_1_K(iZ2,iZ3,iZ4)
    !     uDR_Option(iZ2,iZ3,iZ4,iDR_PL_Theta_2) = Theta_2_K(iZ2,iZ3,iZ4)
    !     uDR_Option(iZ2,iZ3,iZ4,iDR_PL_dEnergy) = Zero

    !   END DO
    !   END DO
    !   END DO

    ! END IF

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

  SUBROUTINE ComputePointValuesZ_Single( InterpMat_Z, U_Q, U_P )

        REAL(DP), INTENT(in)  :: &
         InterpMat_Z(nPT_Z,nDOFZ)
        REAL(DP), INTENT(in)  :: &
          U_Q(nDOFZ)
        REAL(DP), INTENT(out) :: &
          U_P(nPT_Z)
    
        REAL(DP) :: SUM1
        INTEGER  :: iNodeZ, iP_Z
    
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
    REAL(DP), INTENT(out) :: &
      U_P(nPT_X, &
          iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3), &
          iZ_B0(4):iZ_E0(4))

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

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iNodeZ
    REAL(DP) :: SUM1, SUM2

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

  END SUBROUTINE ComputeCellAverage
  
  REAL(DP) FUNCTION GammaFun( E, F1, F2, F3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

    REAL(DP), INTENT(in) :: E, F1, F2, F3
    REAL(DP), INTENT(in) :: Gm_dd_11, Gm_dd_22, Gm_dd_33

    GammaFun &
      = E - SQRT( F1**2 / Gm_dd_11 + F2**2 / Gm_dd_22 + F3**2 / Gm_dd_33 )

    RETURN
  END FUNCTION GammaFun

  SUBROUTINE SolveTheta_Bisection &
    ( E_P, F1_P, F2_P, F3_P, E_K, F1_K, F2_K, F3_K, &
      Gm_dd_11, Gm_dd_22, Gm_dd_33, Theta )

    REAL(DP), INTENT(in)  :: E_P, F1_P, F2_P, F3_P
    REAL(DP), INTENT(in)  :: E_K, F1_K, F2_K, F3_K
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
            ( E_K, F1_K, F2_K, F3_K, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) - Min_2

    x_b = One
    f_b = GammaFun &
            ( E_P, F1_P, F2_P, F3_P, Gm_dd_11, Gm_dd_22, Gm_dd_33 ) - Min_2

    dx = One

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED .AND. ITERATION < ITERATION_MAX )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx

      f_c = GammaFun &
              ( x_c * E_P  + ( One - x_c ) * E_K,  &
                x_c * F1_P + ( One - x_c ) * F1_K, &
                x_c * F2_P + ( One - x_c ) * F2_K, &
                x_c * F3_P + ( One - x_c ) * F3_K, &
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

  SUBROUTINE ApplyLimiter( Theta, U_K, U_Q )

    REAL(DP), INTENT(in)    :: Theta, U_K
    REAL(DP), INTENT(inout) :: U_Q(nDOFZ)

    INTEGER :: iNodeZ

    DO iNodeZ = 1, nDOFZ
      U_Q(iNodeZ) = Theta * U_Q(iNodeZ) + ( One - Theta ) * U_K
    END DO

    RETURN
  END SUBROUTINE ApplyLimiter

END MODULE TwoMoment_PositivityLimiterModule_FMC