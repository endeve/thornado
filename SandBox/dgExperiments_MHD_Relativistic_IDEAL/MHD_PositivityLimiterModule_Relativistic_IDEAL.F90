!> Initialize, finalize, and apply positivity limiter from
!> Wu and Tang, (2017)
MODULE MHD_PositivityLimiterModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three, &
    SqrtTiny
  USE UtilitiesModule, ONLY: &
    IsCornerCell
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
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_SqrtGm
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    iCM_D, &
    iCM_S1, &
    iCM_S2, &
    iCM_S3, &
    iCM_E, &
    iCM_B1, &
    iCM_B2, &
    iCM_B3

  USE MPI

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializePositivityLimiter_MHD_Relativistic_IDEAL
  PUBLIC :: FinalizePositivityLimiter_MHD_Relativistic_IDEAL
  PUBLIC :: ApplyPositivityLimiter_MHD_Relativistic_IDEAL

  LOGICAL               :: UsePositivityLimiter
  INTEGER, PARAMETER    :: nPS = 7  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total number of Positive Points
  REAL(DP)              :: Min_1, Min_2, Min_3
  REAL(DP), ALLOCATABLE :: L_X(:,:)

  INTERFACE ComputePointValues
    MODULE PROCEDURE ComputePointValues_SingleField
    MODULE PROCEDURE ComputePointValues_ManyFields
  END INTERFACE ComputePointValues

CONTAINS


  SUBROUTINE InitializePositivityLimiter_MHD_Relativistic_IDEAL &
    ( UsePositivityLimiter_Option, Verbose_Option, Min_1_Option, Min_2_Option, Min_3_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option, &
                                      Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Min_2_Option, Min_3_Option

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

    Min_3 = - HUGE( One )
    IF( PRESENT( Min_3_Option ) ) &
      Min_3 = Min_3_Option

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Positivity Limiter (MHD, Relativistic, IDEAL)'
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
      WRITE(*,'(A6,A27,ES11.4E3)') &
        '', 'Min_3: ', Min_3
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

  END SUBROUTINE InitializePositivityLimiter_MHD_Relativistic_IDEAL


  SUBROUTINE FinalizePositivityLimiter_MHD_Relativistic_IDEAL

    DEALLOCATE( L_X )

  END SUBROUTINE FinalizePositivityLimiter_MHD_Relativistic_IDEAL


  !> Iterate through the entire spatial domain and apply the positivity
  !> limiter from Wu and Tang, (2017) to each element.
  !> @param Theta_D minimum value to ensure physical density
  !> @param Theta_q minimum value to ensure physical internal energy-density
  !>        and velocity
  !> @param Theta_psi minimum value to ensure ...
  SUBROUTINE ApplyPositivityLimiter_MHD_Relativistic_IDEAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iCM, iGF, iPT, nX_K, nCM_K
    REAL(DP) :: Min_D, Min_D_K, Min_q, Min_q_K, Min_Psi, Min_Psi_K, &
                Theta_D, Theta_q, Theta_Psi_P
    REAL(DP) :: q_P, q_K, Psi_P, Psi_K

    LOGICAL :: NegativeStates(3  ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: SqrtGm(1:nDOFX   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    REAL(DP) :: U_Q(1:nDOFX,1:nCM,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: U_P(1:nPT  ,1:nCM,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(1:nCM        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: G_K(1:nGF        ,iX_B0(1):iX_E0(1), &
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

    REAL(DP) :: Theta_Psi(iX_B0(1):iX_E0(1), &
                          iX_B0(2):iX_E0(2), &
                          iX_B0(3):iX_E0(3))

    REAL(DP) :: D_P, S1_P, S2_P, S3_P, E_P, B1_P, B2_P, B3_P, &
                g1_P, g2_P, g3_P

    REAL(DP) :: D_K, S1_K, S2_K, S3_K, E_K, B1_K, B2_K, B3_K

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    nX_K  = PRODUCT( iX_E0 - iX_B0 + 1 )
    nCM_K = nCM * nX_K

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

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM
    DO iNX = 1, nDOFX

      U_Q(iNX,iCM,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCM)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL ComputePointValues( nCM_K, iX_B0, iX_E0, U_Q, U_P )

    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h1Q, h1P )
    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h2Q, h2P )
    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h3Q, h3P )

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

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCM = 1, nCM

      IF( IsCornerCell( iX_B1, iX_E1, iX1, iX2, iX3 ) ) CYCLE

      U_K(iCM,iX1,iX2,iX3) &
        = SUM( WeightsX_q * SqrtGm(:,iX1,iX2,iX3) * U_Q(:,iCM,iX1,iX2,iX3) ) &
          / SUM( WeightsX_q * SqrtGm(:,iX1,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iGF = 1, nGF

      IF( IsCornerCell( iX_B1, iX_E1, iX1, iX2, iX3 ) ) CYCLE

      G_K(iGF,iX1,iX2,iX3) &
        = SUM( WeightsX_q * SqrtGm(:,iX1,iX2,iX3) * G(:,iX1,iX2,iX3,iGF) ) &
          / SUM( WeightsX_q * SqrtGm(:,iX1,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    ! --- Limit Mass-Density ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( IsCornerCell( iX_B1, iX_E1, iX1, iX2, iX3 ) ) CYCLE

      Min_D = Min_1

      Min_D_K = HUGE( 1.0_DP )
      DO iPT = 1, nPT

        Min_D_K &
          = MIN( Min_D_K, U_P(iPT,iCM_D,iX1,iX2,iX3) )

      END DO

      NegativeStates(1,iX1,iX2,iX3) = .FALSE.

      IF( Min_D_K .LT. Min_D )THEN

        NegativeStates(1,iX1,iX2,iX3) = .TRUE.

        Theta_D &
          =   ( U_K(iCM_D,iX1,iX2,iX3) - Min_D ) &
            / ( U_K(iCM_D,iX1,iX2,iX3) - Min_D_K )

        DO iNX = 1, nDOFX

          U_Q(iNX,iCM_D,iX1,iX2,iX3) &
            = U_K(iCM_D,iX1,iX2,iX3) &
                + Theta_D * ( U_Q(iNX,iCM_D,iX1,iX2,iX3) &
                                - U_K(iCM_D,iX1,iX2,iX3) )

        END DO

      END IF

    END DO
    END DO
    END DO

    ! --- Recompute point values ---

    CALL ComputePointValues( nCM_K, iX_B0, iX_E0, U_Q, U_P )

    ! --- Limit q-function ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( IsCornerCell( iX_B1, iX_E1, iX1, iX2, iX3 ) ) CYCLE

      Min_q = Min_2

      Min_q_K = HUGE( 1.0_DP )
      DO iPT = 1, nPT

        q_P = Computeq( U_P(iPT,iCM_D ,iX1,iX2,iX3),  &
                        U_P(iPT,iCM_S1,iX1,iX2,iX3),  &
                        U_P(iPT,iCM_S2,iX1,iX2,iX3),  &
                        U_P(iPT,iCM_S3,iX1,iX2,iX3),  &
                        U_P(iPT,iCM_E ,iX1,iX2,iX3)   &
                        + U_P(iPT,iCM_D ,iX1,iX2,iX3), &
                        g1P(iPT,       iX1,iX2,iX3),  &
                        g2P(iPT,       iX1,iX2,iX3),  &
                        g3P(iPT,       iX1,iX2,iX3) )

      Min_q_K &
          = MIN( Min_q_K, q_P )

      END DO

      NegativeStates(2,iX1,iX2,iX3) = .FALSE.

      IF( Min_q_K .LT. Min_q )THEN

        NegativeStates(2,iX1,iX2,iX3) = .TRUE.

        q_K = Computeq( U_K(iCM_D,  iX1,iX2,iX3),    &
                        U_K(iCM_S1, iX1,iX2,iX3),    &
                        U_K(iCM_S2, iX1,iX2,iX3),    &
                        U_K(iCM_S3, iX1,iX2,iX3),    &
                        U_K(iCM_E,  iX1,iX2,iX3)     &
                        + U_K(iCM_D, iX1,iX2,iX3),   &
                        G_K(iGF_h_1,iX1,iX2,iX3)**2, &
                        G_K(iGF_h_2,iX1,iX2,iX3)**2, &
                        G_K(iGF_h_3,iX1,iX2,iX3)**2 )

        Theta_q &
          =   ( q_K - Min_q ) &
            / ( q_K - Min_q_K )

        DO iNX = 1, nDOFX

          U_Q(iNX,iCM_D,iX1,iX2,iX3) &
            = U_K(iCM_D,iX1,iX2,iX3) &
                + Theta_q * ( U_Q(iNX,iCM_D,iX1,iX2,iX3) &
                                - U_K(iCM_D,iX1,iX2,iX3) )

          U_Q(iNX,iCM_S1,iX1,iX2,iX3) &
            = U_K(iCM_S1,iX1,iX2,iX3)  &
                + Theta_q * ( U_Q(iNX,iCM_S1,iX1,iX2,iX3) &
                                - U_K(iCM_S1,iX1,iX2,iX3) )

          U_Q(iNX,iCM_S2,iX1,iX2,iX3) &
            = U_K(iCM_S2,iX1,iX2,iX3) &
                + Theta_q * ( U_Q(iNX,iCM_S2,iX1,iX2,iX3) &
                                - U_K(iCM_S2,iX1,iX2,iX3) )

          U_Q(iNX,iCM_S3,iX1,iX2,iX3) &
            = U_K(iCM_S3,iX1,iX2,iX3) &
                + Theta_q * ( U_Q(iNX,iCM_S3,iX1,iX2,iX3) &
                                - U_K(iCM_S3,iX1,iX2,iX3) )

          U_Q(iNX,iCM_E,iX1,iX2,iX3) &
            = U_K(iCM_E,iX1,iX2,iX3) + U_K(iCM_D,iX1,iX2,iX3) &
                + Theta_q * ( U_Q(iNX,iCM_E,iX1,iX2,iX3) &
                                + U_Q(iNX,iCM_D,iX1,iX2,iX3) &
                                - U_K(iCM_E,iX1,iX2,iX3) &
                                - U_K(iCM_D,iX1,iX2,iX3) ) &
                - U_Q(iNX,iCM_D,iX1,iX2,iX3)

        END DO

      END IF

    END DO
    END DO
    END DO

    ! --- Recompute point values ---

    CALL ComputePointValues( nCM_K, iX_B0, iX_E0, U_Q, U_P )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( IsCornerCell( iX_B1, iX_E1, iX1, iX2, iX3 ) ) CYCLE

      Min_Psi = Zero

      Min_Psi_K = HUGE( 1.0_DP )
      Theta_Psi(iX1,iX2,iX3) = HUGE( 1.0_DP )

        D_K = U_K(iCM_D, iX1,iX2,iX3)
       S1_K = U_K(iCM_S1,iX1,iX2,iX3)
       S2_K = U_K(iCM_S2,iX1,iX2,iX3)
       S3_K = U_K(iCM_S3,iX1,iX2,iX3)
        E_K = U_K(iCM_E, iX1,iX2,iX3) + U_K(iCM_D, iX1,iX2,iX3)
       B1_K = U_K(iCM_B1,iX1,iX2,iX3)
       B2_K = U_K(iCM_B2,iX1,iX2,iX3)
       B3_K = U_K(iCM_B3,iX1,iX2,iX3)

      DO iPT = 1, nPT

          D_P = U_P(iPT,iCM_D, iX1,iX2,iX3)
         S1_P = U_P(iPT,iCM_S1,iX1,iX2,iX3)
         S2_P = U_P(iPT,iCM_S2,iX1,iX2,iX3)
         S3_P = U_P(iPT,iCM_S3,iX1,iX2,iX3)
          E_P = U_P(iPT,iCM_E, iX1,iX2,iX3) + U_P(iPT,iCM_D, iX1,iX2,iX3)
         B1_P = U_P(iPT,iCM_B1,iX1,iX2,iX3)
         B2_P = U_P(iPT,iCM_B2,iX1,iX2,iX3)
         B3_P = U_P(iPT,iCM_B3,iX1,iX2,iX3)

         g1_P = g1P(iPT,iX1,iX2,iX3)
         g2_P = g2P(iPT,iX1,iX2,iX3)
         g3_P = g3P(iPT,iX1,iX2,iX3)

        Psi_P = ComputePsi( D_P, S1_P, S2_P, S3_P, E_P - Min_3, &
                            B1_P, B2_P, B3_P, g1_P, g2_P, g3_P )

        Min_Psi_K &
            = MIN( Min_Psi_K, Psi_P )

        CALL SolveTheta_Bisection &
               ( D_P, S1_P, S2_P, S3_P, E_P, &
                 B1_P, B2_P, B3_P, &
                 D_K, S1_K, S2_K, S3_K, E_K, &
                 B1_K, B2_K, B3_K, &
                 g1_P, g2_P, g3_P, Min_3, Theta_Psi_P )

        Theta_Psi(iX1,iX2,iX3) = MIN( Theta_Psi(iX1,iX2,iX3), Theta_Psi_P )

      END DO

      NegativeStates(3,iX1,iX2,iX3) = .FALSE.

      IF( Min_Psi_K .LT. Min_Psi )THEN

        NegativeStates(3,iX1,iX2,iX3) = .TRUE.

        DO iNX = 1, nDOFX

          U_Q(iNX,iCM_D,iX1,iX2,iX3) &
            = U_K(iCM_D,iX1,iX2,iX3) &
                + Theta_Psi(iX1,iX2,iX3) * ( U_Q(iNX,iCM_D,iX1,iX2,iX3) &
                                               - U_K(iCM_D,iX1,iX2,iX3) )

          U_Q(iNX,iCM_S1,iX1,iX2,iX3) &
            = U_K(iCM_S1,iX1,iX2,iX3)  &
                + Theta_Psi(iX1,iX2,iX3) * ( U_Q(iNX,iCM_S1,iX1,iX2,iX3) &
                                               - U_K(iCM_S1,iX1,iX2,iX3) )

          U_Q(iNX,iCM_S2,iX1,iX2,iX3) &
            = U_K(iCM_S2,iX1,iX2,iX3) &
                + Theta_Psi(iX1,iX2,iX3) * ( U_Q(iNX,iCM_S2,iX1,iX2,iX3) &
                                               - U_K(iCM_S2,iX1,iX2,iX3) )

          U_Q(iNX,iCM_S3,iX1,iX2,iX3) &
            = U_K(iCM_S3,iX1,iX2,iX3) &
                + Theta_Psi(iX1,iX2,iX3) * ( U_Q(iNX,iCM_S3,iX1,iX2,iX3) &
                                               - U_K(iCM_S3,iX1,iX2,iX3) )

          U_Q(iNX,iCM_E,iX1,iX2,iX3) &
            = U_K(iCM_E,iX1,iX2,iX3) + U_K(iCM_D,iX1,iX2,iX3) &
                + Theta_Psi(iX1,iX2,iX3) * ( U_Q(iNX,iCM_E,iX1,iX2,iX3) &
                                               + U_Q(iNX,iCM_D,iX1,iX2,iX3) &
                                               - U_K(iCM_E,iX1,iX2,iX3) &
                                               - U_K(iCM_D,iX1,iX2,iX3) ) &
                - U_Q(iNX,iCM_D,iX1,iX2,iX3)

          U_Q(iNX,iCM_B1,iX1,iX2,iX3) &
            = U_K(iCM_B1,iX1,iX2,iX3)  &
                + Theta_Psi(iX1,iX2,iX3) * ( U_Q(iNX,iCM_B1,iX1,iX2,iX3) &
                                               - U_K(iCM_B1,iX1,iX2,iX3) )

          U_Q(iNX,iCM_B2,iX1,iX2,iX3) &
            = U_K(iCM_B2,iX1,iX2,iX3) &
                + Theta_Psi(iX1,iX2,iX3) * ( U_Q(iNX,iCM_B2,iX1,iX2,iX3) &
                                               - U_K(iCM_B2,iX1,iX2,iX3) )

          U_Q(iNX,iCM_B3,iX1,iX2,iX3) &
            = U_K(iCM_B3,iX1,iX2,iX3) &
                + Theta_Psi(iX1,iX2,iX3) * ( U_Q(iNX,iCM_B3,iX1,iX2,iX3) &
                                               - U_K(iCM_B3,iX1,iX2,iX3) )

        END DO

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyPositivityLimiter_MHD_Relativistic_IDEAL


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
    REAL(DP), INTENT(in)  :: U_Q(1:nDOFX,1:nCM,iX_B0(1):iX_E0(1), &
                                               iX_B0(2):iX_E0(2), &
                                               iX_B0(3):iX_E0(3) )
    REAL(DP), INTENT(out) :: U_P(1:nPT  ,1:nCM,iX_B0(1):iX_E0(1), &
                                               iX_B0(2):iX_E0(2), &
                                               iX_B0(3):iX_E0(3) )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, N, nDOFX, One, L_X, nPT, &
             U_Q, nDOFX, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues_ManyFields


  FUNCTION Computeq &
    ( D, S1, S2, S3, E, g1, g2, g3 ) RESULT( q )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, E, g1, g2, g3

    REAL(DP) :: q

    q = E - D * SQRT( One + ( S1**2 / g1 + S2**2 / g2 + S3**2 / g3 ) / D**2 )

    RETURN
  END FUNCTION Computeq


  FUNCTION ComputePsi &
    ( D, S1, S2, S3, E, B1, B2, B3, g1, g2, g3 ) RESULT( Psi )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, E, B1, B2, B3, g1, g2, g3

    REAL(DP) :: SSq, BSq, SdotB, Phi, Psi

    SSq = S1**2 / g1 + S2**2 / g2 + S3**2 / g3

    BSq = g1 * B1**2 + g2 * B2**2 + g3 * B3**2

    SdotB = S1 * B1 + S2 * B2 + S3 * B3

    Phi = ComputePhi( D, SSq, E, BSq )

    !PRINT*, 'Phi + BSq - E: ', Phi + BSq - E

    Psi = ( Phi - Two * ( BSq - E ) ) * SQRT( Phi + BSq - E ) &
          - SQRT( ( 27.0_DP / Two ) * ( D**2 * BSq + SdotB**2 ) )

    RETURN
  END FUNCTION ComputePsi


   FUNCTION ComputePhi &
    ( D, SSq, E, BSq ) RESULT( Phi )

    REAL(DP), INTENT(in) :: D, SSq, E, BSq

    REAL(DP) :: Phi

    Phi = SQRT( ( BSq - E )**2 + Three * ( E**2 - D**2 - SSq ) )

    RETURN
  END FUNCTION ComputePhi


  SUBROUTINE SolveTheta_Bisection &
    ( D_P, S1_P, S2_P, S3_P, E_P, B1_P, B2_P, B3_P, &
      D_K, S1_K, S2_K, S3_K, E_K, B1_K, B2_K, B3_K, &
      g1_P, g2_P, g3_P, Min_3, Theta_Psi )

    REAL(DP), INTENT(in)    :: D_P, S1_P, S2_P, S3_P, E_P, B1_P, B2_P, B3_P, &
                               D_K, S1_K, S2_K, S3_K, E_K, B1_K, B2_K, B3_K, &
                               g1_P, g2_P, g3_P, Min_3
    REAL(DP), INTENT(out)   :: Theta_Psi

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0e-3_DP

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = ComputePsi &
            ( x_a *   D_P  + ( One - x_a ) *   D_K,         &
              x_a *  S1_P  + ( One - x_a ) *  S1_K,         &
              x_a *  S2_P  + ( One - x_a ) *  S2_K,         &
              x_a *  S3_P  + ( One - x_a ) *  S3_K,         &
              x_a *   E_P  + ( One - x_a ) *   E_K - Min_3, &
              x_a *  B1_P  + ( One - x_a ) *  B1_K,         &
              x_a *  B2_P  + ( One - x_a ) *  B2_K,         &
              x_a *  B3_P  + ( One - x_a ) *  B3_K,         &
              g1_P, g2_P, g3_P )

    x_b = One * ( One + 10.0_DP * dx_min )
    f_b = ComputePsi &
            ( x_b *   D_P  + ( One - x_b ) *   D_K,         &
              x_b *  S1_P  + ( One - x_b ) *  S1_K,         &
              x_b *  S2_P  + ( One - x_b ) *  S2_K,         &
              x_b *  S3_P  + ( One - x_b ) *  S3_K,         &
              x_b *   E_P  + ( One - x_b ) *   E_K - Min_3, &
              x_b *  B1_P  + ( One - x_b ) *  B1_K,         &
              x_b *  B2_P  + ( One - x_b ) *  B2_K,         &
              x_b *  B3_P  + ( One - x_b ) *  B3_K,         &
              g1_P, g2_P, g3_P )

    IF( .NOT. f_a * f_b < 0 )THEN
      !PRINT*, 'Cannot perform bisection in positivity limiter.'
      !PRINT*, 'x_a: ', x_a
      !PRINT*, 'x_b: ', x_b
      !PRINT*, 'f_a: ', f_a
      !PRINT*, 'f_b: ', f_b
      Theta_Psi = One
      RETURN
    END IF

    dx = x_b - x_a

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx
      f_c = ComputePsi &
              ( x_c *   D_P  + ( One - x_c ) *   D_K,         &
                x_c *  S1_P  + ( One - x_c ) *  S1_K,         &
                x_c *  S2_P  + ( One - x_c ) *  S2_K,         &
                x_c *  S3_P  + ( One - x_c ) *  S3_K,         &
                x_c *  E_P   + ( One - x_c ) *   E_K - Min_3, &
                x_c *  B1_P  + ( One - x_c ) *  B1_K,         &
                x_c *  B2_P  + ( One - x_c ) *  B2_K,         &
                x_c *  B3_P  + ( One - x_c ) *  B3_K,         &
                g1_P, g2_P, g3_P )

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx .LT. dx_min ) CONVERGED = .TRUE.

      IF( ITERATION .GT. MAX_IT .AND. .NOT. CONVERGED )THEN

        STOP

      END IF

    END DO

    Theta_Psi = x_c

  END SUBROUTINE SolveTheta_Bisection


END MODULE MHD_PositivityLimiterModule_Relativistic_IDEAL
