!> Initialize, finalize, and apply positivity limiter from
!> Qin et al., (2016), JCP, 315, 323
MODULE Euler_PositivityLimiterModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, Zero, Half, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesX, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodesX_q, &
    nDOFX_X1, NodesX1, &
    nDOFX_X2, NodesX2, &
    nDOFX_X3, NodesX3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, LX_X1_Up, &
    LX_X2_Dn, LX_X2_Up, &
    LX_X3_Dn, LX_X3_Up
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, iGF_h_2, iGF_h_3, &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_SqrtGm

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: Euler_InitializePositivityLimiter_Relativistic
  PUBLIC :: Euler_FinalizePositivityLimiter_Relativistic
  PUBLIC :: Euler_ApplyPositivityLimiter_Relativistic

  LOGICAL               :: UsePositivityLimiter
  INTEGER, PARAMETER    :: nPS = 7  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total number of Positive Points
  REAL(DP)              :: Min_1, Min_2
  REAL(DP), ALLOCATABLE :: U_PP(:,:), G_PP(:,:)

CONTAINS


  SUBROUTINE Euler_InitializePositivityLimiter_Relativistic &
    ( Min_1_Option, Min_2_Option, UsePositivityLimiter_Option, Verbose_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Min_2_Option
    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option, &
                                      Verbose_Option

    INTEGER :: iDim
    LOGICAL :: Verbose

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    UsePositivityLimiter = .TRUE.
    IF( PRESENT( UsePositivityLimiter_Option ) ) &
      UsePositivityLimiter = UsePositivityLimiter_Option

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,'(A)') '    INFO: Euler_InitializePositivityLimiter_Relativistic'
      WRITE(*,'(A)') '    ----------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)') &
        '', 'Use Positivity Limiter: ', UsePositivityLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES10.3E3)') '', 'Min_1: ', Min_1
      WRITE(*,'(A4,A27,ES10.3E3)') '', 'Min_2: ', Min_2
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

    ALLOCATE( U_PP(nPT,nCF) )
    ALLOCATE( G_PP(nPT,nGF) )

  END SUBROUTINE Euler_InitializePositivityLimiter_Relativistic


  SUBROUTINE Euler_FinalizePositivityLimiter_Relativistic

    DEALLOCATE( U_PP )
    DEALLOCATE( G_PP )

  END SUBROUTINE Euler_FinalizePositivityLimiter_Relativistic


  !> Iterate through the entire spatial domain and apply the positivity
  !> limiter from Qin et al., (2016), JCP, 315, 323 to each element.
  !> @param Theta_D minimum value to ensure physical density
  !> @param Theta_q minimum value to ensure physical internal energy-density
  !>        and velocity
  !> @todo Modify to accomodate GR
  SUBROUTINE Euler_ApplyPositivityLimiter_Relativistic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iGF, iP
    REAL(DP) :: Theta_D, Theta_q, Theta_P
    REAL(DP) :: Min_K, Min_D, Min_ESq
    REAL(DP) :: U_q(nDOFX,nCF), G_q(nDOFX,nGF), U_K(nCF), q(nPT), SSq(nPT)
    REAL(DP) :: chi

    ! --- For de-bugging ---
    REAL(DP) :: q_K(nPT)

    IF( nDOFX == 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      U_q(1:nDOFX,1:nCF) = U(1:nDOFX,iX1,iX2,iX3,1:nCF)
      G_q(1:nDOFX,1:nGF) = G(1:nDOFX,iX1,iX2,iX3,1:nGF)

      ! --- Check if cell-average of density is physical ---

      U_K(iCF_D) = SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) &
                          * U_q(1:nDOFX,iCF_D) ) &
                     / SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) )

      IF( U_K(iCF_D) .GT. Zero )THEN

        Min_D = Min_1 * U_K(iCF_D)

        CALL ComputePointValues( U_q(1:nDOFX,iCF_D), U_PP(1:nPT,iCF_D) )

        Min_K = MINVAL( U_PP(1:nPT,iCF_D) )

        IF( Min_K .LT. Min_D )THEN

          ! --- Limit point-values of density towards cell-average ---
          Theta_D &
            = ( U_K(iCF_D) - Min_D ) / ( U_K(iCF_D) - Min_K )
          U_q(1:nDOFX,iCF_D) &
            = U_K(iCF_D) + Theta_D * ( U_q(1:nDOFX,iCF_D) - U_K(iCF_D) )

        END IF

        ! --- Modify point-values of q if necessary ---

        DO iCF = 1, nCF
          CALL ComputePointValues( U_q(1:nDOFX,iCF), U_PP(1:nPT,iCF) )
        END DO

        DO iGF = iGF_h_1, iGF_h_3
          CALL ComputePointValues( G_q(1:nDOFX,iGF), G_PP(1:nPT,iGF) )
        END DO
        CALL ComputeGeometryX_FromScaleFactors( G_PP(1:nPT,1:nGF) )

        U_K(iCF_E) &
          = ( SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) &
                     * U_q(1:nDOFX,iCF_E) ) &
                / SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) ) )

        IF( U_K(iCF_E) .LT. Zero ) STOP 'U_K(iCF_E) < 0'

        Min_ESq = Min_2 * U_K(iCF_E)**2

        CALL Computeq &
               ( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), Min_ESq, q(1:nPT) )

        IF( ANY( q(1:nPT) .LT. Zero ) )THEN

          DO iCF = 1, nCF
            U_K(iCF) = SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) &
                              * U_q(1:nDOFX,iCF) ) &
                         / SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) )
          END DO

          ! --- Compute cell-average of q ---
          DO iP = 1, nPT
            CALL Computeq( 1, U_K(1:nCF), G_PP(iP,1:nGF), Min_ESq, q_K(iP) )
          END DO

          ! --- Artificially inject some thermal energy if any q_K < 0 ---
          IF( ANY( q_K(1:nPT) .LT. Zero ) )THEN

            WRITE(*,*)
            WRITE(*,'(A)') 'Adding internal energy'
            WRITE(*,'(A,3I5.4)') 'iX1, iX2, iX3 = ', iX1, iX2, iX3

            ! --- Use maximum value of SSq(1:nPT) to determine how much
            !     internal energy to inject ---
            DO iP = 1, nPT
              SSq(iP) = U_K(iCF_S1)**2 / G_PP(iP,iGF_Gm_dd_11) &
                           + U_K(iCF_S2)**2 / G_PP(iP,iGF_Gm_dd_22) &
                           + U_K(iCF_S3)**2 / G_PP(iP,iGF_Gm_dd_33)
            END DO

            WRITE(*,'(A,ES24.16E3)') 'tau_K (old):', U_K(iCF_E)
            chi = ( Min_2 * U_K(iCF_E) - U_K(iCF_D) &
                      + SQRT( U_K(iCF_D)**2 + MAXVAL( SSq ) &
                                + Min_ESq ) ) / U_K(iCF_E) - One
            U_K(iCF_E) = U_K(iCF_E) * ( One + chi )
            WRITE(*,'(A,ES13.6E3)') &
              'Fractional amount of internal energy added: ', chi

            DO iCF = 1, nCF
              U_q(1:nDOFX,iCF) = U_K(iCF)
            END DO

          ELSE

            ! --- Solve for Theta_q such that all point-values
            !     of q are positive ---
            Theta_q = One
            DO iP = 1, nPT
              IF( q(iP) .LT. Zero )THEN
                CALL SolveTheta_Bisection &
                  ( U_PP(iP,1:nCF), U_K(1:nCF), G_PP(iP,1:nGF), Min_ESq, &
                    Theta_P, iX1, iX2, iX3, iP )
                Theta_q = MIN( Theta_q, Theta_P )
              END IF
            END DO

            ! --- Limit all variables towards cell-average ---
            DO iCF = 1, nCF
              U_q(1:nDOFX,iCF) &
                = U_K(iCF) + Theta_q * ( U_q(1:nDOFX,iCF) - U_K(iCF) )
            END DO

          END IF

        END IF ! q < 0

        U(1:nDOFX,iX1,iX2,iX3,1:nCF) = U_q(1:nDOFX,1:nCF)

      ELSE

        WRITE(*,'(A)') &
          'Warning: Euler_ApplyPositivityLimiterModule_Relativistic'
        WRITE(*,'(A)') 'Cell-average of density <= 0'
        WRITE(*,'(A)') 'Setting all values to cell-average'

        DO iCF = 1, nCF
          U_K(iCF) &
            = SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) &
                     * U_q(1:nDOFX,iCF) ) &
              / SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) )
          U(1:nDOFX,iX1,iX2,iX3,iCF) = U_K(iCF)
        END DO

      END IF

      ! --- Ensure that all point-values of q are positive after limiting ---
      DO iCF = 1, nCF
        CALL ComputePointValues( U(1:nDOFX,iX1,iX2,iX3,iCF), U_PP(1:nPT,iCF) )
      END DO
      DO iGF = iGF_h_1, iGF_h_3
        CALL ComputePointValues( G(1:nDOFX,iX1,iX2,iX3,iGF), G_PP(1:nPT,iGF) )
      END DO
      CALL ComputeGeometryX_FromScaleFactors( G_PP(1:nPT,1:nGF) )

      CALL Computeq &
             ( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), Min_ESq, q(1:nPT) )
      IF( ANY( q(1:nPT) .LT. Zero ) )THEN
        WRITE(*,*)
        WRITE(*,'(A)') &
            'FATAL ERROR: Euler_ApplyPositivityLimiterModule_Relativistic'
        WRITE(*,'(A)') &
            '------------------------------------------------------------'
        WRITE(*,'(A)') 'q < 0 after all limiting'
        WRITE(*,*)
        WRITE(*,'(A,3I5.4)') 'iX1, iX2, iX3 = ', iX1, iX2, iX3
        WRITE(*,*)
        DO iP = 1, SIZE( q(1:nPT) )
          WRITE(*,'(A,I2.2,ES25.16E3)') &
            'iP, q(iP) = ', iP, q(iP)
          WRITE(*,'(A,6ES12.3E3)') &
            'U_PP = ', U_PP(iP,:)
          WRITE(*,'(A,ES23.16E3)') &
            'SSq = ', &
            U_PP(iP,iCF_S1)**2 / G_PP(iP,iGF_Gm_dd_11) &
            + U_PP(iP,iCF_S2)**2 / G_PP(iP,iGF_Gm_dd_22) &
            + U_PP(iP,iCF_S3)**2 / G_PP(iP,iGF_Gm_dd_33)
          WRITE(*,'(A,3ES11.3E3)') &
            'G_PP(Gm11:Gm33) = ', G_PP(iP,iGF_Gm_dd_11:iGF_Gm_dd_33)
          WRITE(*,*)
        END DO
        ! --- Compute cell-average of q ---
        DO iP = 1, nPT
          CALL Computeq( 1, U_K(1:nCF), G_PP(iP,1:nGF), Min_ESq, q_K(iP) )
        END DO
        IF( ANY( q_K(1:nPT) .LT. Zero ) )THEN
          WRITE(*,'(A)') 'At least one q_K(1:nPT) < 0'
          DO iP = 1, nPT
            WRITE(*,'(A,I2.2,ES25.16E3)') 'iP, q_K(iP) = ', iP, q_K(iP)
          END DO
        WRITE(*,*)
        END IF
        WRITE(*,'(A,F18.16)')    'Theta_q = ', Theta_q
        WRITE(*,*)
        WRITE(*,'(A)') 'U(1:nDOFX,iX1,iX2,iX3,iCF)'
        WRITE(*,'(A)') '--------------------------'
        DO iCF = 1, nCF
          WRITE(*,'(A,I1)') 'iCF = ', iCF
          WRITE(*,*) U(1:nDOFX,iX1,iX2,iX3,iCF)
          WRITE(*,*)
        END DO
        STOP ''
      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE Euler_ApplyPositivityLimiter_Relativistic


  SUBROUTINE ComputePointValues( X_Q, X_P )

    REAL(DP), INTENT(in)  :: X_Q(nDOFX)
    REAL(DP), INTENT(out) :: X_P(nPT)

    INTEGER :: iOS

    X_P(1:nDOFX) = X_Q(1:nDOFX)

    IF( SUM( nPP(2:3) ) > 0 )THEN

      ! --- X1 ---

      iOS = nPP(1)

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X1), 1 )

      iOS = iOS + nPP(2)

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X1), 1 )

    END IF

    IF( SUM( nPP(4:5) ) > 0 )THEN

      ! --- X2 ---

      iOS = SUM( nPP(1:3) )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X2), 1 )

      iOS = iOS + nPP(4)

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X2), 1 )

    END IF

    IF( SUM( nPP(6:7) ) > 0 )THEN

      ! --- X3 ---

      iOS = SUM( nPP(1:5) )

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X3), 1 )

      iOS = iOS + nPP(6)

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Up, nDOFX_X3, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X3), 1 )

    END IF

  END SUBROUTINE ComputePointValues


  SUBROUTINE Computeq( N, U, G, Min_ESq, q )

    INTEGER,  INTENT(in)  :: N
    REAL(DP), INTENT(in)  :: U(N,nCF), G(N,nGF), Min_ESq
    REAL(DP), INTENT(out) :: q(N)

    q = qFun( U(1:N,iCF_D ), &
              U(1:N,iCF_S1), &
              U(1:N,iCF_S2), &
              U(1:N,iCF_S3), &
              U(1:N,iCF_E ), &
              G(1:N,iGF_Gm_dd_11), &
              G(1:N,iGF_Gm_dd_22), &
              G(1:N,iGF_Gm_dd_33), &
              Min_ESq )

  END SUBROUTINE Computeq


  REAL(DP) ELEMENTAL FUNCTION qFun &
    ( D, S1, S2, S3, tau, Gm11, Gm22, Gm33, Min_ESq )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, tau, Gm11, Gm22, Gm33, Min_ESq

    qFun = tau + D &
             - SQRT( D**2 + ( S1**2 / Gm11 + S2**2 / Gm22 + S3**2 / Gm33 ) &
                       + Min_ESq )

    RETURN
  END FUNCTION qFun


  SUBROUTINE SolveTheta_Bisection &
    ( U_Q, U_K, G_Q, Min_ESq, Theta_P, iX1, iX2, iX3, iP )

    REAL(DP), INTENT(in)  :: U_Q(nCF), U_K(nCF), G_Q(nGF), Min_ESq
    REAL(DP), INTENT(out) :: Theta_P

    ! --- For de-bugging ---
    INTEGER,  INTENT(in) :: iX1, iX2, iX3, iP
    INTEGER :: iCF, iGF

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

!!$    WRITE(*,*)
!!$    WRITE(*,'(A,I3)') 'iP = ', iP
!!$    WRITE(*,'(A)') 'U_Q = np.array( ['
!!$    DO iCF = 1, nCF-1
!!$      WRITE(*,'(ES24.16E3,A)') U_Q(iCF), ','
!!$    END DO
!!$    WRITE(*,'(ES24.16E3,A)') U_Q(iCF), ' ] )'
!!$    WRITE(*,*)
!!$    WRITE(*,'(A)') 'U_K = np.array( ['
!!$    DO iCF = 1, nCF-1
!!$      WRITE(*,'(ES24.16E3,A)') U_K(iCF), ','
!!$    END DO
!!$    WRITE(*,'(ES24.16E3,A)') U_K(iCF), ' ] )'
!!$    WRITE(*,*)
!!$    WRITE(*,'(A)') 'G_Q = np.array( ['
!!$    DO iGF = iGF_Gm_dd_11, iGF_Gm_dd_33-1
!!$      WRITE(*,'(ES24.16E3,A)') G_Q(iCF), ','
!!$    END DO
!!$    WRITE(*,'(ES24.16E3,A)') G_Q(iGF), ' ] )'
!!$    WRITE(*,*)

    x_a = Zero
    f_a = qFun &
            ( x_a * U_Q(iCF_D)  + ( One - x_a ) * U_K(iCF_D),         &
              x_a * U_Q(iCF_S1) + ( One - x_a ) * U_K(iCF_S1),        &
              x_a * U_Q(iCF_S2) + ( One - x_a ) * U_K(iCF_S2),        &
              x_a * U_Q(iCF_S3) + ( One - x_a ) * U_K(iCF_S3),        &
              x_a * U_Q(iCF_E)  + ( One - x_a ) * U_K(iCF_E),         &
              G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33), Min_ESq )

    x_b = One
    f_b = qFun &
            ( x_b * U_Q(iCF_D)  + ( One - x_b ) * U_K(iCF_D),         &
              x_b * U_Q(iCF_S1) + ( One - x_b ) * U_K(iCF_S1),        &
              x_b * U_Q(iCF_S2) + ( One - x_b ) * U_K(iCF_S2),        &
              x_b * U_Q(iCF_S3) + ( One - x_b ) * U_K(iCF_S3),        &
              x_b * U_Q(iCF_E)  + ( One - x_b ) * U_K(iCF_E),         &
              G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33), Min_ESq )

    IF( .NOT. f_a * f_b < 0 )THEN

      WRITE(*,'(A6,A)') &
        '', 'SolveTheta_Bisection (Euler_PositivityLimiterModule_Relativistic_IDEAL):'

      WRITE(*,*) 'iP = ', iP
      WRITE(*,*) 'iX1, iX2, iX3 = ', iX1, iX2, iX3
      WRITE(*,*) 'U_Q: ', U_Q
      WRITE(*,*) 'U_K: ', U_K
      WRITE(*,*) 'G_Q: ', G_Q(iGF_Gm_dd_11:iGF_Gm_dd_33)
      WRITE(*,'(A8,A,I3.3)') &
        '', 'Error: No Root in Interval'
      WRITE(*,'(A8,A,2ES15.6e3)') &
        '', 'x_a, x_b = ', x_a, x_b
      WRITE(*,'(A8,A,2ES15.6e3)') &
        '', 'f_a, f_b = ', f_a, f_b
      STOP

    END IF

    dx = x_b - x_a

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx
      f_c = qFun &
            ( x_c * U_Q(iCF_D)  + ( One - x_c ) * U_K(iCF_D),         &
              x_c * U_Q(iCF_S1) + ( One - x_c ) * U_K(iCF_S1),        &
              x_c * U_Q(iCF_S2) + ( One - x_c ) * U_K(iCF_S2),        &
              x_c * U_Q(iCF_S3) + ( One - x_c ) * U_K(iCF_S3),        &
              x_c * U_Q(iCF_E)  + ( One - x_c ) * U_K(iCF_E),         &
              G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33), Min_ESq )

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx < dx_min ) CONVERGED = .TRUE.

      IF( ITERATION > MAX_IT .AND. .NOT. CONVERGED )THEN

        WRITE(*,'(6x,A)') &
          'SolveTheta_Bisection (Euler_ApplyPositivityLimiter_Relativistic):'
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


END MODULE Euler_PositivityLimiterModule_Relativistic_IDEAL
