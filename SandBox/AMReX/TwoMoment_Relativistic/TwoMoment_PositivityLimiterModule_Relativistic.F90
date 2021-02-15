MODULE TwoMoment_PositivityLimiterModule_Relativistic

  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nNodesZ, &
    nDOFE, nDOFX, nNodesX
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
    nGF, iGF_SqrtGm, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, &
    iCF_E, iCF_NE
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_TwoMoment
  PUBLIC :: FinalizePositivityLimiter_TwoMoment
  PUBLIC :: ApplyPositivityLimiter_TwoMoment

  LOGICAL               :: UsePositivityLimiter
  LOGICAL               :: Verbose
  INTEGER               :: N_R
  INTEGER,    PARAMETER :: nPS = 9  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total Number of Positive Points
  REAL(DP)              :: Min_1, Max_1, Min_2
  REAL(DP),   PARAMETER :: One_EPS = One - 1.0d1 * EPSILON( One )
  REAL(DP), ALLOCATABLE :: InterpMat(:,:)

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

    nPP    = 0
    nPP(1) = PRODUCT( nNodesZ )

    DO i = 1, 4

      IF( nNodesZ(i) > 1 )THEN

        nPP(2*i:2*i+1) = PRODUCT( nNodesZ, MASK = [1,2,3,4] .NE. i )

      END IF

    END DO
    nPT = SUM( nPP )



    ALLOCATE( InterpMat(nPT,nDOFZ) )


    InterpMat = Zero

    DO iNodeZ = 1, nDOFZ


      InterpMat(iNodeZ,iNodeZ) = One


      IF( SUM( nPP(2:3) ) > 0 )THEN

        iOS = nPP(1)

        InterpMat(iOS+1:iOS+nDOF_E,iNodeZ) = L_E_Dn(1:nDOF_E,iNodeZ)


        iOS = iOS + nPP(2)

        InterpMat(iOS+1:iOS+nDOF_E,iNodeZ) = L_E_Up(1:nDOF_E,iNodeZ)


      END IF

      IF( SUM( nPP(4:5) ) > 0 )THEN

        iOS = SUM( nPP(1:3) )
        InterpMat(iOS+1:iOS+nDOF_X1,iNodeZ) = L_X1_Dn(1:nDOF_X1,iNodeZ)

        iOS = iOS + nPP(4)
        InterpMat(iOS+1:iOS+nDOF_X1,iNodeZ) = L_X1_Up(1:nDOF_X1,iNodeZ)

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        iOS = SUM( nPP(1:5) )
        InterpMat(iOS+1:iOS+nDOF_X2,iNodeZ) = L_X2_Dn(1:nDOF_X2,iNodeZ)

        iOS = iOS + nPP(6)
        InterpMat(iOS+1:iOS+nDOF_X2,iNodeZ) = L_X2_Up(1:nDOF_X2,iNodeZ)

      END IF

      IF( SUM( nPP(8:9) ) > 0 )THEN

        iOS = SUM( nPP(1:7) )
        InterpMat(iOS+1:iOS+nDOF_X3,iNodeZ) = L_X3_Dn(1:nDOF_X3,iNodeZ)

        iOS = iOS + nPP(8)
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
    REAL(DP), INTENT(in) :: &
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


    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    END IF

    IF( .NOT. UsePositivityLimiter .OR. nDOFZ == 1 ) RETURN


    IF (Verbose) THEN

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

      D_Q (:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_D)
      S1_Q (:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_S1)
      S2_Q (:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_S2)
      S3_Q (:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_S3)
      E_Q (:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_E)
      NE_Q (:,iZ2,iZ3,iZ4) = U_F(:,iZ2,iZ3,iZ4,iCF_NE)

      G_11_Q (:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Gm_dd_11)
      G_22_Q (:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Gm_dd_22)
      G_33_Q (:,iZ2,iZ3,iZ4) = GX(:,iZ2,iZ3,iZ4,iGF_Gm_dd_33)

    END DO
    END DO
    END DO

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
                 (D_Q(iNodeX,iZ2,iZ3,iZ4), &
                 S1_Q(iNodeX,iZ2,iZ3,iZ4), &
                 S2_Q(iNodeX,iZ2,iZ3,iZ4), &
                 S3_Q(iNodeX,iZ2,iZ3,iZ4), &
                 E_Q(iNodeX,iZ2,iZ3,iZ4), &
                 NE_Q(iNodeX,iZ2,iZ3,iZ4), &
                 DP_Q(iNodeX,iZ2,iZ3,iZ4), &
                 V1_Q(iNodeX,iZ2,iZ3,iZ4), &
                 V2_Q(iNodeX,iZ2,iZ3,iZ4), &
                 V3_Q(iNodeX,iZ2,iZ3,iZ4), &
                 EP_Q(iNodeX,iZ2,iZ3,iZ4), &
                 NEP_Q(iNodeX,iZ2,iZ3,iZ4), &
                 G_11_Q(iNodeX,iZ2,iZ3,iZ4), &
                 G_22_Q(iNodeX,iZ2,iZ3,iZ4), &
                 G_33_Q(iNodeX,iZ2,iZ3,iZ4) )

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

    END DO
    END DO
    END DO
    END DO
    END DO

    IF( RecomputePointValues )THEN

      CALL ComputePointValuesZ( iZ_B0, iZ_E0, N_Q , N_P  )

    END IF

    ! --- Ensure Positive "Gamma" ---

    CALL ComputePointValuesX( iX_B0, iX_E0, D_Q , D_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, S1_Q , S1_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, S2_Q , S2_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, S3_Q , S3_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, E_Q , E_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, NE_Q , NE_P  )


    CALL ComputePointValuesX( iX_B0, iX_E0, G_11_Q , G_11_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, G_22_Q , G_22_P  )
    CALL ComputePointValuesX( iX_B0, iX_E0, G_33_Q , G_33_P  )

!



    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iP = 1, nPT_X

        CALL ComputePrimitive_Euler_Relativistic &
                 (D_P(iP,iZ2,iZ3,iZ4), &
                 S1_P(iP,iZ2,iZ3,iZ4), &
                 S2_P(iP,iZ2,iZ3,iZ4), &
                 S3_P(iP,iZ2,iZ3,iZ4), &
                 E_P(iP,iZ2,iZ3,iZ4), &
                 NE_P(iP,iZ2,iZ3,iZ4), &
                 DP_P(iP,iZ2,iZ3,iZ4), &
                 V1_P(iP,iZ2,iZ3,iZ4), &
                 V2_P(iP,iZ2,iZ3,iZ4), &
                 V3_P(iP,iZ2,iZ3,iZ4), &
                 EP_P(iP,iZ2,iZ3,iZ4), &
                 NEP_P(iP,iZ2,iZ3,iZ4), &
                 G_11_P(iP,iZ2,iZ3,iZ4), &
                 G_22_P(iP,iZ2,iZ3,iZ4), &
                 G_33_P(iP,iZ2,iZ3,iZ4) )


      END DO

    END DO
    END DO
    END DO
    Change = 0 
print*, "Cell Average"
    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iP = 1, nPT_X

        CALL CheckRealizability( N_K(iZ1,iZ2,iZ3,iZ4,iS) , &
                                 G1_K(iZ1,iZ2,iZ3,iZ4,iS) , &                                                                   
                                 G2_K(iZ1,iZ2,iZ3,iZ4,iS) , &                                                                   
                                 G3_K(iZ1,iZ2,iZ3,iZ4,iS) , &                                                                   
                                 V1_P(iP,iZ2,iZ3,iZ4), &
                                 V2_P(iP,iZ2,iZ3,iZ4), &
                                 V3_P(iP,iZ2,iZ3,iZ4), &
                                 G_11_P(iP,iZ2,iZ3,iZ4), &
                                 G_22_P(iP,iZ2,iZ3,iZ4), &
                                 G_33_P(iP,iZ2,iZ3,iZ4), &
                                 iZ1, iZ2, Min_1, Min_2 )

      END DO



    END DO
    END DO
    END DO
    END DO
    END DO

n=0
m=0
    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      Gamma_Min = Min_2
      Theta_2   = One

      DO iP = 1, nPT

        CALL PointsZtoPointsX( nNodesX(1), iP, iP_X )

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
                 G_33_P(iP_X,iZ2,iZ3,iZ4) )

        Gamma_Min = MIN( Gamma, Gamma_Min )

        IF( Gamma_Min < Min_2 )THEN
print*, iZ1, iZ2
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

    END DO
    END DO
    END DO
    END DO
    END DO
print*, n, m


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


!    CALL ComputePointValuesZ( iZ_B0, iZ_E0, N_Q , N_P  )
!    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G1_Q, G1_P )
!    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G2_Q, G2_P )
!    CALL ComputePointValuesZ( iZ_B0, iZ_E0, G3_Q, G3_P )
!!print*, "Point Values"
!    DO iS = 1, nSpecies
!    DO iZ4 = iZ_B0(4), iZ_E0(4)
!    DO iZ3 = iZ_B0(3), iZ_E0(3)
!    DO iZ2 = iZ_B0(2), iZ_E0(2)
!    DO iZ1 = iZ_B0(1), iZ_E0(1)
!
!      DO iP = 1, nPT
!
!        CALL PointsZtoPointsX( nNodesX(1), iP, iP_X )
!        CALL CheckRealizability( N_P (iP,iZ1,iZ2,iZ3,iZ4,iS), &
!                                 G1_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
!                                 G2_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
!                                 G3_P(iP,iZ1,iZ2,iZ3,iZ4,iS), &
!                                 V1_P(iP_X,iZ2,iZ3,iZ4), &
!                                 V2_P(iP_X,iZ2,iZ3,iZ4), &
!                                 V3_P(iP_X,iZ2,iZ3,iZ4), &
!                                 G_11_P(iP_X,iZ2,iZ3,iZ4), &
!                                 G_22_P(iP_X,iZ2,iZ3,iZ4), &
!                                 G_33_P(iP_X,iZ2,iZ3,iZ4), &
!                                 iZ1, iZ2, Min_1, Min_2 )
!
!      END DO
!
!
!
!    END DO
!    END DO
!    END DO
!    END DO
!    END DO
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


  REAL(DP) FUNCTION GammaFun( N, G1, G2, G3, V1, V2, V3, Gm_dd_11, Gm_dd_22, Gm_dd_33 )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: N, G1, G2, G3, V1, V2 ,V3, Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: G

    G = Gm_dd_11 * G1**2 + Gm_dd_22 * G2**2 + Gm_dd_33 * G3**2 &
      - ( V1**2 * Gm_dd_11**2 * G1**2 + V2**2 * Gm_dd_22**2 * G2**2 + V3**2 * Gm_dd_33**2 * G3**2 ) &
      - 2.0_DP * ( V1 * V2 * Gm_dd_11 * Gm_dd_22 * G1 * G2 + V1 * V3 * Gm_dd_11 * Gm_dd_33 * G1 * G3 &
      + V2 * V3 * Gm_dd_22 * Gm_dd_33 * G2 * G3 )
    G = SQRT(G)  
    GammaFun = N - G 
    RETURN
  END FUNCTION GammaFun


  SUBROUTINE SolveTheta_Bisection &
    ( N_P, G1_P, G2_P, G3_P, N_K, G1_K, G2_K, G3_K, &
      V1_P, V2_P, V3_P, Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P, &
      Theta )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)  :: N_P, G1_P, G2_P, G3_P
    REAL(DP), INTENT(in)  :: N_K, G1_K, G2_K, G3_K
    REAL(DP), INTENT(in)  :: V1_P, V2_P, V3_P, Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P
    REAL(DP), INTENT(out) :: Theta

    INTEGER,  PARAMETER :: ITERATION_MAX = 12
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = GammaFun( N_K, G1_K, G2_K, G3_K, V1_P, V2_P, V3_P, &
                    Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P ) - Min_2

    x_b = One
    f_b = GammaFun( N_P, G1_P, G2_P, G3_P, V1_P, V2_P, V3_P, &
                    Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P ) - Min_2

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
                Gm_dd_11_P, Gm_dd_22_P, Gm_dd_33_P ) - Min_2
!
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

  SUBROUTINE PointsZtoPointsX( nNodes, iPT_Z, iPT_X )


    INTEGER, INTENT(in)  :: nNodes, iPT_Z
    INTEGER, INTENT(out) :: iPT_X

    IF (nNodes .EQ. 2) THEN
      
      IF(iPT_Z .EQ. 9 .OR. iPT_Z .EQ. 10) THEN
        iPT_X = 1
      ELSE IF(iPT_Z .EQ. 1 .OR. iPT_Z .EQ. 2 .OR. iPT_Z .EQ. 5 .OR. iPT_Z .EQ. 7) THEN
        iPT_X = 2      
      ELSE IF(iPT_Z .EQ. 3 .OR. iPT_Z .EQ. 4 .OR. iPT_Z .EQ. 6 .OR. iPT_Z .EQ. 8) THEN
        iPT_X = 3
      ELSE IF(iPT_Z .EQ. 11 .OR. iPT_Z .EQ. 12) THEN
        iPT_X = 4
      END IF

    ELSE IF(nNodes .EQ. 3) THEN

      IF(iPT_Z .EQ. 16 .OR. iPT_Z .EQ. 17 .OR. iPT_Z .EQ. 18) THEN
        iPT_X = 1
      ELSE IF(iPT_Z .EQ. 1 .OR. iPT_Z .EQ. 2 .OR. iPT_Z .EQ. 3 .OR. iPT_Z .EQ. 10 .OR. iPT_Z .EQ. 13) THEN
        iPT_X = 2      
      ELSE IF(iPT_Z .EQ. 4 .OR. iPT_Z .EQ. 5 .OR. iPT_Z .EQ. 6 .OR. iPT_Z .EQ. 11 .OR. iPT_Z .EQ. 14) THEN
        iPT_X = 3
      ELSE IF(iPT_Z .EQ. 7 .OR. iPT_Z .EQ. 8 .OR. iPT_Z .EQ. 9 .OR. iPT_Z .EQ. 12 .OR. iPT_Z .EQ. 15) THEN
        iPT_X = 4
      ELSE IF(iPT_Z .EQ. 19 .OR. iPT_Z .EQ. 20 .OR. iPT_Z .EQ. 21) THEN
        iPT_X = 5
      END IF

    END IF

  END SUBROUTINE

  SUBROUTINE CheckRealizability(N, G1, G2, G3, V1, V2, V3, G_11, G_22, G_33, iZ1, iZ2, Min_1,Min_2)


    REAL(DP), INTENT(inout)  :: N, G1, G2, G3
    REAL(DP), INTENT(in)  :: V1, V2, V3, G_11, G_22, G_33
    INTEGER, INTENT(in)   :: iZ1, iZ2
    REAL(DP), INTENT(in)  :: Min_1, Min_2

    REAL(DP) :: GammaOut, W

    GammaOut = GammaFun( N, G1, G2, G3, V1, V2, V3, &
                    G_11, G_22, G_33 )

    IF ( N .LT. Min_1 ) THEN
!     print*, "N Unrealizable"
!     print*, iZ1,iZ2, V1
!     print*, "N = ", N
!     print*, "G1 = ", G1
!     print*, "G2 = ", G2
!     print*, "G3 = ", G3
    END IF


    IF ( GammaOut .LT. Min_2 ) THEN
      W = 1.0_DP / ( 1.0_DP - ( G_11 * V1**2 + G_22 * V2**2 + G_33 * V3**2 ) )
    
   !  print*, "N - G Unrealizable"
   !  print*, iZ1,iZ2, V1
   !  print*, "N = ", N
   !  print*, "G1 = ", G1
   !  print*, "G2 = ", G2
   !  print*, "G3 = ", G3
   !  print*, "Gamma", GammaOut
     G1 = 0.999_DP * W * ( N - Min_2 ) 
    END IF

  END SUBROUTINE CheckRealizability





END MODULE TwoMoment_PositivityLimiterModule_Relativistic
