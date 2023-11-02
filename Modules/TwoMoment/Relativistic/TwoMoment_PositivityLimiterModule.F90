MODULE TwoMoment_PositivityLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nNodesZ, &
    nDOFE, nDOFX, nNodesX, nNodesE
  USE ReferenceElementModule, ONLY: &
    nDOF_E, &
    nDOF_X1, &
    nDOF_X2, &
    nDOF_X3, &
    Weights_Q
  USE ReferenceElementModuleX,          ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_q
  USE ReferenceElementModule_Lagrange, ONLY: &
    L_E_Dn,  L_E_Up, &
    L_X1_Dn, L_X1_Up, &
    L_X2_Dn, L_X2_Up, &
    L_X3_Dn, L_X3_Up
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep2
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_SqrtGm, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3, &
    iGF_h_1, iGF_h_2, iGF_h_3
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, &
    iCF_E, iCF_NE
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE MeshModule, ONLY: &
    MeshX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_TwoMoment
  PUBLIC :: FinalizePositivityLimiter_TwoMoment
  PUBLIC :: ApplyPositivityLimiter_TwoMoment

  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: Verbose
  INTEGER               :: N_R
  INTEGER,    PARAMETER :: nPS_Z = 9  ! Number of Positive Point Sets
  INTEGER               :: nPP_Z(nPS_Z) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total Number of Positive Points
  REAL(DP)              :: Min_1, Max_1, Min_2
  REAL(DP),   PARAMETER :: One_EPS = One - 1.0d1 * EPSILON( One )
  REAL(DP), ALLOCATABLE :: InterpMat(:,:)
  INTEGER,  ALLOCATABLE :: PointZ2X(:)

  INTEGER, PARAMETER    :: nPS_X = 7  ! Number of Positive Point Sets Euler
  INTEGER               :: nPP_X(nPS_X) ! Number of Positive Points Per Set Euler
  INTEGER               :: nPT_X      ! Total number of Positive Points Euler
  INTEGER               :: N_X
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

    INTEGER :: i, iNodeZ, iNodeX, iOS, iOS_X, iDim
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iP_X, iP_Z
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

    nPP_Z    = 0
    nPP_Z(1) = PRODUCT( nNodesZ )

    DO i = 1, 4

      IF( nNodesZ(i) > 1 )THEN

        nPP_Z(2*i:2*i+1) = PRODUCT( nNodesZ, MASK = [1,2,3,4] .NE. i )

      END IF

    END DO
    nPT = SUM( nPP_Z )



    ALLOCATE( InterpMat(nPT,nDOFZ) )


    InterpMat = Zero

    DO iNodeZ = 1, nDOFZ


      InterpMat(iNodeZ,iNodeZ) = One


      IF( SUM( nPP_Z(2:3) ) > 0 )THEN

        iOS = nPP_Z(1)

        InterpMat(iOS+1:iOS+nDOF_E,iNodeZ) = L_E_Dn(1:nDOF_E,iNodeZ)


        iOS = iOS + nPP_Z(2)

        InterpMat(iOS+1:iOS+nDOF_E,iNodeZ) = L_E_Up(1:nDOF_E,iNodeZ)


      END IF

      IF( SUM( nPP_Z(4:5) ) > 0 )THEN

        iOS = SUM( nPP_Z(1:3) )
        InterpMat(iOS+1:iOS+nDOF_X1,iNodeZ) = L_X1_Dn(1:nDOF_X1,iNodeZ)

        iOS = iOS + nPP_Z(4)
        InterpMat(iOS+1:iOS+nDOF_X1,iNodeZ) = L_X1_Up(1:nDOF_X1,iNodeZ)

      END IF

      IF( SUM( nPP_Z(6:7) ) > 0 )THEN

        iOS = SUM( nPP_Z(1:5) )
        InterpMat(iOS+1:iOS+nDOF_X2,iNodeZ) = L_X2_Dn(1:nDOF_X2,iNodeZ)

        iOS = iOS + nPP_Z(6)
        InterpMat(iOS+1:iOS+nDOF_X2,iNodeZ) = L_X2_Up(1:nDOF_X2,iNodeZ)

      END IF

      IF( SUM( nPP_Z(8:9) ) > 0 )THEN

        iOS = SUM( nPP_Z(1:7) )
        InterpMat(iOS+1:iOS+nDOF_X3,iNodeZ) = L_X3_Dn(1:nDOF_X3,iNodeZ)

        iOS = iOS + nPP_Z(8)
        InterpMat(iOS+1:iOS+nDOF_X3,iNodeZ) = L_X3_Up(1:nDOF_X3,iNodeZ)

      END IF

    END DO

    nPP_X(1:nPS_X) = 0
    nPP_X(1)     = PRODUCT( nNodesX(1:3) )

    DO iDim = 1, 3

      IF( nNodesX(iDim) > 1 )THEN

        nPP_X(2*iDim:2*iDim+1) &
          = PRODUCT( nNodesX(1:3), MASK = [1,2,3] .NE. iDim )

      END IF

    END DO

    nPT_X = SUM( nPP_X(1:nPS_X) )

    ALLOCATE( InterpMat_X(nPT_X,nDOFX) )
    InterpMat_X = Zero

    DO iNodeX = 1, nDOFX


      InterpMat_X(iNodeX,iNodeX) = One
      IF( SUM( nPP_X(2:3) ) > 0 )THEN

        ! --- X1 ---

        iOS_X = nPP_X(1)
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X1,iNodeX) = LX_X1_Dn(1:nDOFX_X1,iNodeX)

        iOS_X = iOS_X + nPP_X(2)
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X1,iNodeX) = LX_X1_Up(1:nDOFX_X1,iNodeX)

      END IF

      IF( SUM( nPP_X(4:5) ) > 0 )THEN

        ! --- X2 ---

        iOS_X = SUM( nPP_X(1:3) )
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X2,iNodeX) = LX_X2_Dn(1:nDOFX_X2,iNodeX)

        iOS_X = iOS_X + nPP_X(4)
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X2,iNodeX) = LX_X2_Up(1:nDOFX_X2,iNodeX)

      END IF

      IF( SUM( nPP_X(6:7) ) > 0 )THEN

        ! --- X3 ---

        iOS_X = SUM( nPP_X(1:5) )
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X3,iNodeX) = LX_X3_Dn(1:nDOFX_X3,iNodeX)

        iOS_X = iOS_X + nPP_X(6)
        InterpMat_X(iOS_X+1:iOS_X+nDOFX_X3,iNodeX) = LX_X3_Up(1:nDOFX_X3,iNodeX)

      END IF


    END DO
    ! --- Index Map from Energy-Position to Position ---

    ALLOCATE( PointZ2X(nPT) )

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

    PRINT*, "FinalizePositivityLimiter_TwoMoment"

    DEALLOCATE( InterpMat )

  END SUBROUTINE FinalizePositivityLimiter_TwoMoment


  SUBROUTINE ApplyPositivityLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U_F, U_R, Verbose_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE (1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX (1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                  iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U_F(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
                  iZ_B1(4):iZ_E1(4),1:nCF)
    REAL(DP), INTENT(inout) :: &
      U_R(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                  iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      Verbose_Option

    LOGICAL  :: RecomputePointValues
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP, iP_X, n, m, Change
    INTEGER  :: iX_B0(3), iX_E0(3)
    INTEGER  :: iNodeZ, iNodeE, iNodeX
    REAL(DP) :: Min_K, Max_K, Theta_1, Theta_2, Theta_P
    REAL(DP) :: Gamma, Gamma_Min
    REAL(DP) :: Tau_Q(nDOFZ, &
                      iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                      iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: N_Q(nDOFZ, &
                    iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    nSpecies), &
                N_P(nPT  , &
                    iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    nSpecies), &
                N_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                    iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                    nSpecies)

    REAL(DP) :: G1_Q(nDOFZ, &
                       iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                       nSpecies), &
                G1_P(nPT  , &
                       iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                       nSpecies), &
                G1_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                       iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                       nSpecies)
    REAL(DP) :: G2_Q(nDOFZ, &
                     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies), &
                G2_P(nPT  , &
                     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies), &
                G2_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies)
    REAL(DP) :: G3_Q(nDOFZ, &
                     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies), &
                G3_P(nPT  , &
                     iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies), &
                G3_K(iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
                     nSpecies)

    REAL(DP) :: D_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                D_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: S1_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                S1_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: S2_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                S2_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: S3_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                S3_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: E_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                E_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: NE_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                NE_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: DP_Q(nDOFX  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 V1_Q(nDOFX  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 V2_Q(nDOFX  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 V3_Q(nDOFX  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 EP_Q(nDOFX  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 NEP_Q(nDOFX  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: DP_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 V1_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 V2_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 V3_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 EP_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                 NEP_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: G_11_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                G_11_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: G_22_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                G_22_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: G_33_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                G_33_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: h_d_1_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                h_d_1_P(nPT_X, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: h_d_2_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                h_d_2_P(nPT_X, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: h_d_3_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                h_d_3_P(nPT_X, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: A_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                A_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    REAL(DP) :: B1_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                B1_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: B2_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                B2_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP) :: B3_Q(nDOFX, &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4)), &
                B3_P(nPT_X  , &
                     iZ_B0(2):iZ_E0(2), &
                     iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    LOGICAL :: &
      RealizableCellAverage &
        (iZ_B0(1):iZ_E0(1), &
         iZ_B0(2):iZ_E0(2), &
         iZ_B0(3):iZ_E0(3), &
         iZ_B0(4):iZ_E0(4), &
         nSpecies)

    INTEGER :: ITERATION_Q(nDOFX,iZ_B0(2):iZ_E0(2), &
                                 iZ_B0(3):iZ_E0(3), &
                                 iZ_B0(4):iZ_E0(4))
    INTEGER :: iErr_Q     (nDOFX,iZ_B0(2):iZ_E0(2), &
                                 iZ_B0(3):iZ_E0(3), &
                                 iZ_B0(4):iZ_E0(4))

    INTEGER :: ITERATION_P(nPT_X,iZ_B0(2):iZ_E0(2), &
                                 iZ_B0(3):iZ_E0(3), &
                                 iZ_B0(4):iZ_E0(4))
    INTEGER :: iErr_P     (nPT_X,iZ_B0(2):iZ_E0(2), &
                                 iZ_B0(3):iZ_E0(3), &
                                 iZ_B0(4):iZ_E0(4))

    INTEGER :: ErrorExists

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    END IF

    IF( .NOT. UsePositivityLimiter .OR. nDOFZ == 1 ) RETURN

    IF( Verbose )THEN

      PRINT*, "      ApplyPositivityLimiter_TwoMoment"

    END IF

    N_R = nSpecies * PRODUCT( iZ_E0 - iZ_B0 + 1 )

    iX_E0(1) = iZ_E0(2)
    iX_E0(2) = iZ_E0(3)
    iX_E0(3) = iZ_E0(4)

    iX_B0(1) = iZ_B0(2)
    iX_B0(2) = iZ_B0(3)
    iX_B0(3) = iZ_B0(4)

    N_X = PRODUCT( iX_E0 - iX_B0 + 1 )

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

      N_Q (:,iZ1,iZ2,iZ3,iZ4,iS) = U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)
      G1_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)
      G2_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)
      G3_Q(:,iZ1,iZ2,iZ3,iZ4,iS) = U_R(:,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)

    END DO
    END DO
    END DO
    END DO
    END DO

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      D_Q (:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_D )
      S1_Q(:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_S1)
      S2_Q(:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_S2)
      S3_Q(:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_S3)
      E_Q (:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_E )
      NE_Q(:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_NE)

      G_11_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Gm_dd_11)
      G_22_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Gm_dd_22)
      G_33_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Gm_dd_33)

      h_d_1_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_h_1)
      h_d_2_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_h_2)
      h_d_3_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_h_3)

      A_Q (:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Alpha )
      B1_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Beta_1)
      B2_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Beta_2)
      B3_Q(:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Beta_3)

    END DO
    END DO
    END DO

    ErrorExists = 0
    DO iZ4    = iZ_B0(4), iZ_E0(4)
    DO iZ3    = iZ_B0(3), iZ_E0(3)
    DO iZ2    = iZ_B0(2), iZ_E0(2)
    DO iNodeX = 1       , nDOFX
      ITERATION_Q(iNodeX,iZ2,iZ3,iZ4) = 0
      iErr_Q     (iNodeX,iZ2,iZ3,iZ4) = 0

      CALL ComputePrimitive_Euler_Relativistic &
             ( D_Q   (iNodeX,iZ2,iZ3,iZ4), &
               S1_Q  (iNodeX,iZ2,iZ3,iZ4), &
               S2_Q  (iNodeX,iZ2,iZ3,iZ4), &
               S3_Q  (iNodeX,iZ2,iZ3,iZ4), &
               E_Q   (iNodeX,iZ2,iZ3,iZ4), &
               NE_Q  (iNodeX,iZ2,iZ3,iZ4), &
               DP_Q  (iNodeX,iZ2,iZ3,iZ4), &
               V1_Q  (iNodeX,iZ2,iZ3,iZ4), &
               V2_Q  (iNodeX,iZ2,iZ3,iZ4), &
               V3_Q  (iNodeX,iZ2,iZ3,iZ4), &
               EP_Q  (iNodeX,iZ2,iZ3,iZ4), &
               NEP_Q (iNodeX,iZ2,iZ3,iZ4), &
               G_11_Q(iNodeX,iZ2,iZ3,iZ4), &
               G_22_Q(iNodeX,iZ2,iZ3,iZ4), &
               G_33_Q(iNodeX,iZ2,iZ3,iZ4), &
               ITERATION_Option = ITERATION_Q(iNodeX,iZ2,iZ3,iZ4), &
               iErr_Option      = iErr_Q     (iNodeX,iZ2,iZ3,iZ4) )

      ErrorExists = ErrorExists + iErr_Q(iNodeX,iZ2,iZ3,iZ4)

    END DO
    END DO
    END DO
    END DO

    IF( ErrorExists .NE. 0 )THEN

      DO iZ4    = iZ_B0(4), iZ_E0(4)
      DO iZ3    = iZ_B0(3), iZ_E0(3)
      DO iZ2    = iZ_B0(2), iZ_E0(2)
      DO iNodeX = 1       , nDOFX

        IF( iErr_Q(iNodeX,iZ2,iZ3,IZ4) .NE. 0 )THEN

          CALL DescribeError_Euler &
            ( iErr_Q(iNodeX,iZ2,iZ3,iZ4), &
              Int_Option = [ ITERATION_Q(iNodeX,iZ2,iZ3,iZ4), 99999999, &
                             iZ_B0(2), iZ_B0(3), iZ_B0(4), &
                             iZ_E0(2), iZ_E0(3), iZ_E0(4), &
                             iNodeX, iZ2, iZ3, iZ4 ], &
              Real_Option = [ MeshX(1) % Center(iZ2), &
                              MeshX(2) % Center(iZ3), &
                              MeshX(3) % Center(iZ4), &
                              MeshX(1) % Width (iZ2), &
                              MeshX(2) % Width (iZ3), &
                              MeshX(3) % Width (iZ4), &
                              D_Q   (iNodeX,iZ2,iZ3,iZ4), &
                              S1_Q  (iNodeX,iZ2,iZ3,iZ4), &
                              S2_Q  (iNodeX,iZ2,iZ3,iZ4), &
                              S3_Q  (iNodeX,iZ2,iZ3,iZ4), &
                              E_Q   (iNodeX,iZ2,iZ3,iZ4), &
                              Ne_Q  (iNodeX,iZ2,iZ3,iZ4), &
                              G_11_Q(iNodeX,iZ2,iZ3,iZ4), &
                              G_22_Q(iNodeX,iZ2,iZ3,iZ4), &
                              G_33_Q(iNodeX,iZ2,iZ3,iZ4) ], &
              Char_Option = [ 'NA' ], &
              Message_Option &
                = 'Calling from ApplyPositivityLimiter_TwoMoment (nDOFX)' )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF ! ErrorExists .NE. 0

    CALL ComputePointValuesZ( iZ_B0, iZ_E0, N_Q , N_P  )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G1_Q, G1_P )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G2_Q, G2_P )
    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G3_Q, G3_P )

    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, N_Q , N_K  )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G1_Q, G1_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G2_Q, G2_K )
    CALL ComputeCellAverage( iZ_B0, iZ_E0, Tau_Q, G3_Q, G3_K )

    CALL ComputePointValuesX( iX_B0, iX_E0, D_Q , D_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, S1_Q , S1_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, S2_Q , S2_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, S3_Q , S3_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, E_Q , E_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, NE_Q , NE_P  )

    CALL ComputePointValuesX( iX_B0, iX_E0, h_d_1_Q , h_d_1_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, h_d_2_Q , h_d_2_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, h_d_3_Q , h_d_3_P  )

    CALL ComputePointValuesX( iX_B0, iX_E0, A_Q , A_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, B1_Q , B1_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, B2_Q , B2_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, B3_Q , B3_P  )

    ErrorExists = 0

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iP  = 1       , nPT_X

      ITERATION_P(iP,iZ2,iZ3,iZ4) = 0
      iErr_P     (iP,iZ2,iZ3,iZ4) = 0


      G_11_P(iP,iZ2,iZ3,iZ4) = MAX( h_d_1_P(iP,iZ2,iZ3,iZ4)**2, SqrtTiny )
      G_22_P(iP,iZ2,iZ3,iZ4) = MAX( h_d_2_P(iP,iZ2,iZ3,iZ4)**2, SqrtTiny )
      G_33_P(iP,iZ2,iZ3,iZ4) = MAX( h_d_3_P(iP,iZ2,iZ3,iZ4)**2, SqrtTiny )


      CALL ComputePrimitive_Euler_Relativistic &
             ( D_P   (iP,iZ2,iZ3,iZ4), &
               S1_P  (iP,iZ2,iZ3,iZ4), &
               S2_P  (iP,iZ2,iZ3,iZ4), &
               S3_P  (iP,iZ2,iZ3,iZ4), &
               E_P   (iP,iZ2,iZ3,iZ4), &
               NE_P  (iP,iZ2,iZ3,iZ4), &
               DP_P  (iP,iZ2,iZ3,iZ4), &
               V1_P  (iP,iZ2,iZ3,iZ4), &
               V2_P  (iP,iZ2,iZ3,iZ4), &
               V3_P  (iP,iZ2,iZ3,iZ4), &
               EP_P  (iP,iZ2,iZ3,iZ4), &
               NEP_P (iP,iZ2,iZ3,iZ4), &
               G_11_P(iP,iZ2,iZ3,iZ4), &
               G_22_P(iP,iZ2,iZ3,iZ4), &
               G_33_P(iP,iZ2,iZ3,iZ4), &
               ITERATION_Option = ITERATION_P(iP,iZ2,iZ3,iZ4), &
               iErr_Option      = iErr_P     (iP,iZ2,iZ3,iZ4) )

      ErrorExists = ErrorExists + iErr_P(iP,iZ2,iZ3,iZ4)

    END DO
    END DO
    END DO
    END DO

    IF( ErrorExists .NE. 0 )THEN

      DO iZ4 = iZ_B0(4), iZ_E0(4)
      DO iZ3 = iZ_B0(3), iZ_E0(3)
      DO iZ2 = iZ_B0(2), iZ_E0(2)
      DO iP  = 1       , nPT_X

        IF( iErr_P(iP,iZ2,iZ3,IZ4) .NE. 0 )THEN

          CALL DescribeError_Euler &
            ( iErr_P(iP,iZ2,iZ3,iZ4), &
              Int_Option = [ ITERATION_P(iP,iZ2,iZ3,iZ4), 99999999, &
                             iZ_B0(2), iZ_B0(3), iZ_B0(4), &
                             iZ_E0(2), iZ_E0(3), iZ_E0(4), &
                             iNodeX, iZ2, iZ3, iZ4 ], &
              Real_Option = [ MeshX(1) % Center(iZ2), &
                              MeshX(2) % Center(iZ3), &
                              MeshX(3) % Center(iZ4), &
                              MeshX(1) % Width (iZ2), &
                              MeshX(2) % Width (iZ3), &
                              MeshX(3) % Width (iZ4), &
                              D_P   (iP,iZ2,iZ3,iZ4), &
                              S1_P  (iP,iZ2,iZ3,iZ4), &
                              S2_P  (iP,iZ2,iZ3,iZ4), &
                              S3_P  (iP,iZ2,iZ3,iZ4), &
                              E_P   (iP,iZ2,iZ3,iZ4), &
                              Ne_P  (iP,iZ2,iZ3,iZ4), &
                              G_11_P(iP,iZ2,iZ3,iZ4), &
                              G_22_P(iP,iZ2,iZ3,iZ4), &
                              G_33_P(iP,iZ2,iZ3,iZ4) ], &
              Char_Option = [ 'NA' ], &
              Message_Option &
                = 'Calling from ApplyPositivityLimiter_TwoMoment (nPT_X)' )

        END IF

      END DO
      END DO
      END DO
      END DO

    END IF ! ErrorExists .NE. 0

    CALL CheckCellAverageRealizability &
          (iZ_B0, iZ_E0, N_K, G1_K, G2_K, G3_K, &
           V1_P, V2_P, V3_P, G_11_P, G_22_P, G_33_P, &
           A_P, B1_P, B2_P, B3_P, RealizableCellAverage)

    ! --- Ensure Bounded Density ---

    RecomputePointValues = .FALSE.

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      IF( RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) )THEN

        Min_K = Min_1

        DO iP = 1, nPT

          Min_K = MIN( Min_K, N_P(iP,iZ1,iZ2,iZ3,iZ4,iS) )

        END DO

        IF( Min_K < Min_1 )THEN

          Theta_1 &
            = MIN( One, &
                   ABS( ( Min_1 - N_K(iZ1,iZ2,iZ3,iZ4,iS) ) &
                        / ( Min_K - N_K(iZ1,iZ2,iZ3,iZ4,iS)+SqrtTiny ) ) )

          Theta_1 = One_EPS * Theta_1

          N_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
            = Theta_1 * N_Q(:,iZ1,iZ2,iZ3,iZ4,iS) &
              + ( One - Theta_1 ) * N_K(iZ1,iZ2,iZ3,iZ4,iS)

          RecomputePointValues = .TRUE.

        END IF

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

n=0
m=0
    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      IF( RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) )THEN
        Gamma_Min = Min_2
        Theta_2   = One

        DO iP = 1, nPT

          !CALL PointsZtoPointsX( nNodesX(1), iP, iP_X )

          iP_X = PointZ2X(iP)

          Gamma &
            = GammaFun &
                ( N_P (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                  G1_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                  G2_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                  G3_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                  V1_P(iP_X,iZ2,iZ3,iZ4), &
                  V2_P(iP_X,iZ2,iZ3,iZ4), &
                  V3_P(iP_X,iZ2,iZ3,iZ4), &
                  G_11_P(iP_X,iZ2,iZ3,iZ4), &
                  G_22_P(iP_X,iZ2,iZ3,iZ4), &
                  G_33_P(iP_X,iZ2,iZ3,iZ4), &
                  A_P(iP_X,iZ2,iZ3,iZ4) , &
                  B1_P(iP_X,iZ2,iZ3,iZ4) , &
                  B2_P(iP_X,iZ2,iZ3,iZ4) , &
                  B3_P(iP_X,iZ2,iZ3,iZ4)  )

          Gamma_Min = MIN( Gamma, Gamma_Min )

          IF( Gamma_Min < Min_2 )THEN

            CALL SolveTheta_Bisection &
                 ( N_P (iP,iZ1,iZ2,iZ3,iZ4,iS), &
                   G1_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                   G2_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                   G3_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
                   N_K (   iZ1,iZ2,iZ3,iZ4,iS), &
                   G1_K(   iZ1,iZ2,iZ3,iZ4,iS), &
                   G2_K(   iZ1,iZ2,iZ3,iZ4,iS), &
                   G3_K(   iZ1,iZ2,iZ3,iZ4,iS), &
                   V1_P(iP_X,iZ2,iZ3,iZ4), &
                   V2_P(iP_X,iZ2,iZ3,iZ4), &
                   V3_P(iP_X,iZ2,iZ3,iZ4), &
                   G_11_P(iP_X,iZ2,iZ3,iZ4), &
                   G_22_P(iP_X,iZ2,iZ3,iZ4), &
                   G_33_P(iP_X,iZ2,iZ3,iZ4), &
                   A_P(iP_X,iZ2,iZ3,iZ4), &
                   B1_P(iP_X,iZ2,iZ3,iZ4), &
                   B2_P(iP_X,iZ2,iZ3,iZ4), &
                   B3_P(iP_X,iZ2,iZ3,iZ4), &
                   Theta_P )
            Theta_2 = MIN( Theta_2, Theta_P )

          END IF

        END DO

        m=m+1

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

          n=n+1

        END IF

      END IF

    END DO
    END DO
    END DO
    END DO
    END DO
!print*, n, m
    IF( .NOT. ALL( RealizableCellAverage ) )THEN
      CALL RecoverRealizableCellAverage &
             ( iZ_B0, iZ_E0, N_K, G1_K, G2_K, G3_K, N_Q, G1_Q, G2_Q, G3_Q, &
               V1_P, V2_P, V3_P, G_11_P, G_22_P, G_33_P, &
               A_P, B1_P, B2_P, B3_P, RealizableCellAverage )

    END IF


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

  END SUBROUTINE ApplyPositivityLimiter_TwoMoment


  SUBROUTINE ComputePointValuesX( iX_B0, iX_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFX,iX_B0(1):iX_E0(1), &
                iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(out) :: &
      U_P(nPT_X,iX_B0(1):iX_E0(1), &
                iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT_X, N_X, nDOFX, One, InterpMat_X, nPT_X, &
             U_Q, nDOFX, Zero, U_P, nPT_X )

  END SUBROUTINE ComputePointValuesX

  SUBROUTINE ComputePointValuesZ( iZ_B0, iZ_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4)
    REAL(DP), INTENT(in)  :: &
      U_Q(nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)
    REAL(DP), INTENT(out) :: &
      U_P(nPT  ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),nSpecies)

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, N_R, nDOFZ, One, InterpMat, nPT, &
             U_Q, nDOFZ, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValuesZ


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


  REAL(DP) FUNCTION GammaFun( N, G1, G2, G3, V1, V2, V3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B1, B2, B3 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: N, G1, G2, G3, V1, V2 ,V3, Gm_dd_11, Gm_dd_22, Gm_dd_33, alp, B1, B2, B3
    REAL(DP) :: absG, G_uu(0:3,0:3), G(0:3), B(3)
    INTEGER :: i, j

    B(1) = B1
    B(2) = B2
    B(3) = B3

    G_uu(0,0) = - 1.0_DP / alp**2
    G_uu(1:3,0) = B(1:3) / alp**2
    G_uu(0,1:3) = B(1:3) / alp**2
    G_uu(1,1) = 1.0_DP / Gm_dd_11 - B(1) * B(1) / alp**2
    G_uu(2,2) = 1.0_DP / Gm_dd_22 - B(2) * B(2) / alp**2
    G_uu(3,3) = 1.0_DP / Gm_dd_33 - B(3) * B(3) / alp**2
    G_uu(1,2) = - B(1) * B(2) / alp**2
    G_uu(1,3) = - B(1) * B(3) / alp**2
    G_uu(2,1) = - B(1) * B(2) / alp**2
    G_uu(3,1) = - B(1) * B(3) / alp**2
    G_uu(2,3) = - B(2) * B(3) / alp**2
    G_uu(3,2) = - B(2) * B(3) / alp**2

    G(0) = ( V1 * G1 + V2 * G2 + V3 * G3 ) / alp
    G(1) = G1
    G(2) = G2
    G(3) = G3

    absG = 0.0_DP

    DO i = 0,3
    DO j = 0,3

      absG = absG + G_uu(i,j) * G(i) * G(j)

    END DO
    END DO

    absG = SQRT(absG)

    GammaFun = N - absG


    RETURN
  END FUNCTION GammaFun


  SUBROUTINE SolveTheta_Bisection &
    ( N_P, G1_P, G2_P, G3_P, N_K, G1_K, G2_K, G3_K, &
      V1_P, V2_P, V3_P, Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P, &
      A_P, B1_P, B2_P, B3_P, Theta )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: N_P, G1_P, G2_P, G3_P
    REAL(DP), INTENT(in)  :: N_K, G1_K, G2_K, G3_K
    REAL(DP), INTENT(in)  :: V1_P, V2_P, V3_P, Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P
    REAL(DP), INTENT(in)  :: A_P, B1_P, B2_P, B3_P
    REAL(DP), INTENT(out) :: Theta

    INTEGER,  PARAMETER :: ITERATION_MAX = 12
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = GammaFun( N_K, G1_K, G2_K, G3_K, V1_P, V2_P, V3_P, &
                    Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P, A_P, B1_P, B2_P, B3_P ) - Min_2

    x_b = One
    f_b = GammaFun( N_P, G1_P, G2_P, G3_P, V1_P, V2_P, V3_P, &
                    Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P, A_P, B1_P, B2_P, B3_P ) - Min_2

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
                V1_P, V2_P, V3_P,                  &
                Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P, &
                A_P, B1_P, B2_P, B3_P  ) - Min_2

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

  SUBROUTINE CheckCellAverageRealizability(iZ_B0, iZ_E0, N_K, G1_K, G2_K, G3_K, &
              V_1_P, V_2_P, V_3_P, G_11_P, G_22_P, G_33_P, &
              A_P, B1_P, B2_P, B3_P, RealizableCellAverage)

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
      G_11_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      G_22_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      G_33_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      A_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      B1_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      B2_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      B3_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      V_1_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      V_2_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      V_3_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    LOGICAL, INTENT(out) :: &
      RealizableCellAverage &
          (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    REAL(DP) :: GammaOut, W
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP_X
    REAL(DP) :: Gamma_K

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) = .TRUE.

      IF( N_K(iZ1,iZ2,iZ3,iZ4,iS) < Min_1 )THEN

       ! PRINT*
       ! PRINT*, "  N_K < Min_1"
       ! PRINT*
       ! PRINT*, "  iZ1,iZ2,iZ3,iZ4,iS = ", iZ1,iZ2,iZ3,iZ4,iS
       ! PRINT*, "  N_K                = ", N_K (iZ1,iZ2,iZ3,iZ4,iS)
       ! PRINT*, "  G1_K               = ", G1_K(iZ1,iZ2,iZ3,iZ4,iS)
       ! PRINT*, "  G2_K               = ", G2_K(iZ1,iZ2,iZ3,iZ4,iS)
       ! PRINT*, "  G3_K               = ", G3_K(iZ1,iZ2,iZ3,iZ4,iS)
       ! PRINT*

        RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) = .FALSE.
      END IF

    END DO
    END DO
    END DO
    END DO
    END DO
    ! --- Check for Negative "Gamma" ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iP_X = 1, nPT_X


        Gamma_K &
          = GammaFun &
              ( N_K (iZ1,iZ2,iZ3,iZ4,iS), &
                G1_K(iZ1,iZ2,iZ3,iZ4,iS), &
                G2_K(iZ1,iZ2,iZ3,iZ4,iS), &
                G3_K(iZ1,iZ2,iZ3,iZ4,iS), &
                V_1_P(iP_X,iZ2,iZ3,iZ4), &
                V_2_P(iP_X,iZ2,iZ3,iZ4), &
                V_3_P(iP_X,iZ2,iZ3,iZ4), &
                G_11_P(iP_X,iZ2,iZ3,iZ4), &
                G_22_P(iP_X,iZ2,iZ3,iZ4), &
                G_33_P(iP_X,iZ2,iZ3,iZ4), &
                A_P(iP_X,iZ2,iZ3,iZ4), &
                B1_P(iP_X,iZ2,iZ3,iZ4), &
                B2_P(iP_X,iZ2,iZ3,iZ4), &
                B3_P(iP_X,iZ2,iZ3,iZ4) )

        IF(  Gamma_K < Min_2 )THEN

       !   PRINT*
       !   PRINT*, "  Gamma_K < Min_2"
       !   PRINT*
       !   PRINT*, "  iZ1,iZ2,iZ3,iZ4,iS = ", iZ1,iZ2,iZ3,iZ4,iS
       !   PRINT*, "  Gamma_K            = ", Gamma_K
       !   PRINT*, "  N_K                = ", N_K (iZ1,iZ2,iZ3,iZ4,iS)
       !   PRINT*, "  G1_K               = ", G1_K(iZ1,iZ2,iZ3,iZ4,iS)
       !   PRINT*, "  G2_K               = ", G2_K(iZ1,iZ2,iZ3,iZ4,iS)
       !   PRINT*, "  G3_K               = ", G3_K(iZ1,iZ2,iZ3,iZ4,iS)
       !   PRINT*, "  V1_P               = ", V_1_P(iP_X,iZ2,iZ3,iZ4)

          RealizableCellAverage(iZ1,iZ2,iZ3,iZ4,iS) = .FALSE.
        END IF

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO



  END SUBROUTINE CheckCellAverageRealizability


  SUBROUTINE RecoverRealizableCellAverage(iZ_B0, iZ_E0, N_K, G1_K, G2_K, G3_K, &
              N_Q, G1_Q, G2_Q, G3_Q, V_1_P, V_2_P, V_3_P, G_11_P, G_22_P, G_33_P, &
              A_P, B1_P, B2_P, B3_P, RealizableCellAverage)

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
      G_11_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      G_22_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      G_33_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      V_1_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      V_2_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      V_3_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      A_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      B1_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      B2_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))
    REAL(DP), INTENT(in) :: &
      B3_P(nPT_X, &
            iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4))

    LOGICAL, INTENT(in) :: &
      RealizableCellAverage &
          (iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4), &
           nSpecies)
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iP_X, i, j
    REAL(DP) :: absG_K, G_uu(0:3,0:3), G(0:3), absG, B(3), alp, Gm_dd_11, Gm_dd_22, Gm_dd_33, V1, V2, V3

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



            Gm_dd_11 = G_11_P(iP_X,iZ2,iZ3,iZ4)
            Gm_dd_22 = G_22_P(iP_X,iZ2,iZ3,iZ4)
            Gm_dd_33 = G_33_P(iP_X,iZ2,iZ3,iZ4)

            V1 = V_1_P(iP_X,iZ2,iZ3,iZ4)
            V2 = V_2_P(iP_X,iZ2,iZ3,iZ4)
            V3 = V_3_P(iP_X,iZ2,iZ3,iZ4)

            alp = A_P(iP_X,iZ2,iZ3,iZ4)
            B(1) = B1_P(iP_X,iZ2,iZ3,iZ4)
            B(2) = B2_P(iP_X,iZ2,iZ3,iZ4)
            B(3) = B3_P(iP_X,iZ2,iZ3,iZ4)


            G_uu(0,0) = - 1.0_DP / alp**2
            G_uu(1:3,0) = B(1:3) / alp**2
            G_uu(0,1:3) = B(1:3) / alp**2
            G_uu(1,1) = 1.0_DP / Gm_dd_11 - B(1) * B(1) / alp**2
            G_uu(2,2) = 1.0_DP / Gm_dd_22 - B(2) * B(2) / alp**2
            G_uu(3,3) = 1.0_DP / Gm_dd_33 - B(3) * B(3) / alp**2
            G_uu(1,2) = - B(1) * B(2) / alp**2
            G_uu(1,3) = - B(1) * B(3) / alp**2
            G_uu(2,1) = - B(1) * B(2) / alp**2
            G_uu(3,1) = - B(1) * B(3) / alp**2
            G_uu(2,3) = - B(2) * B(3) / alp**2
            G_uu(3,2) = - B(2) * B(3) / alp**2

            G(1) = G1_K(iZ1,iZ2,iZ3,iZ4,iS)
            G(2) = G2_K(iZ1,iZ2,iZ3,iZ4,iS)
            G(3) = G3_K(iZ1,iZ2,iZ3,iZ4,iS)
            G(0) = ( V1 * G(1) + V2 * G(2) + V3 * G(3) ) / alp

            absG = 0.0_DP

            DO i = 0,3
            DO j = 0,3

              absG = absG + G_uu(i,j) * G(i) * G(j)

            END DO
            END DO


            absG_K &
              = MAX( absG_K, absG )



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






END MODULE TwoMoment_PositivityLimiterModule
