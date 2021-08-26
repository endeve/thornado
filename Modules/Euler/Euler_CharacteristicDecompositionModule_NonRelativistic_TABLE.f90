MODULE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_SL_CharDecomp

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: Debug = .FALSE.

  REAL(DP), PARAMETER :: dCs_Threshold = 0.1_DP
  REAL(DP), PARAMETER :: D_Threshold   = 1.0d16 * Gram / Centimeter**3
  REAL(DP), PARAMETER, DIMENSION(6,6) :: &
    I_6x6 = RESHAPE( (/1, 0, 0, 0, 0, 0, &
                       0, 1, 0, 0, 0, 0, &
                       0, 0, 1, 0, 0, 0, &
                       0, 0, 0, 1, 0, 0, &
                       0, 0, 0, 0, 1, 0, &
                       0, 0, 0, 0, 0, 1/), &
                       (/6,6/) )

  PUBLIC :: ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
    ( iDim, G, U, R, invR, FJ, cs_T, UseAnalytic_Option )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)
    LOGICAL,  INTENT(in) , OPTIONAL :: UseAnalytic_Option
    REAL(DP), INTENT(out), OPTIONAL :: cs_T
    REAL(DP), INTENT(out), OPTIONAL :: FJ(nCF,nCF)

    LOGICAL  :: UseAnalytic
    REAL(DP) :: Gmdd11, Gmdd22, Gmdd33
    REAL(DP) :: D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
    REAL(DP) :: E, Ne, P, Cs, Cs_table
    REAL(DP) :: K, H, Tau, T, Y, Vsq, CsSq, W, Em, Gm, S
    REAL(DP) :: dPdD, dPdT, dPdY
    REAL(DP) :: dEdD, dEdT, dEdY
    REAL(DP) :: dPdE, dPdDe, dPdTau
    REAL(DP) :: dFdU(nCF,nCF)

    IF ( PRESENT ( UseAnalytic_Option ) ) THEN
      UseAnalytic = UseAnalytic_Option
    ELSE
      UseAnalytic = .TRUE.
    END IF

    CALL TimersStart_Euler( Timer_Euler_SL_CharDecomp )

    CALL ComputePrimitive_Euler_NonRelativistic &
           ( U(iCF_D ), U(iCF_S1), U(iCF_S2), &
             U(iCF_S3), U(iCF_E ), U(iCF_Ne), &
             D, Vu1, Vu2, Vu3, E, Ne, &
             G(iGF_Gm_dd_11), &
             G(iGF_Gm_dd_22), &
             G(iGF_Gm_dd_33) )

    CALL ComputeAuxiliary_Fluid_TABLE &
          ( D, E, Ne, P, T, Y, S, Em, Gm, Cs_table )

    CALL ComputeSpecificInternalEnergy_TABLE &
          ( D, T, Y, Em, dEdD, dEdT, dEdY )

    CALL ComputePressure_TABLE &
          ( D, T, Y, P, dPdD, dPdT, dPdY )

    Gmdd11 = G(iGF_Gm_dd_11)
    Gmdd22 = G(iGF_Gm_dd_22)
    Gmdd33 = G(iGF_Gm_dd_33)

    ! --- Distinguish co-/contra-variant velocities.

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    Tau = 1.0_DP / D

    dPdE   = dPdT / dEdT
    dPdDe  = ( Tau ) * ( dPdY - dEdY * dPdE )
    dPdTau = (dPdDe * Y + dEdD * dPdE - dPdD) / (Tau**2)

    Vsq = Vu1 * Vd1 + Vu2 * Vd2 + Vu3 * Vd3

    CsSq = Tau**2 * ( P * dPdE - dPdTau ) + Y * dPdDe

    IF ( CsSq .LT. Zero .OR. D .GT. D_Threshold ) THEN
      R    = I_6x6
      invR = I_6x6

      IF ( PRESENT( cs_T ) )THEN
        cs_T = Cs_table
      END IF

      IF( PRESENT( FJ ) ) THEN
        FJ = I_6x6
      END IF

      RETURN
    ELSE
      Cs = SQRT( CsSq )
      IF ( ABS( Cs - Cs_table ) / Cs_table .GT. dCs_Threshold ) THEN
        R    = I_6x6
        invR = I_6x6

        IF ( PRESENT( cs_T ) )THEN
          cs_T = Cs_table
        END IF

        IF( PRESENT( FJ ) ) THEN
          FJ = I_6x6
        END IF

        RETURN
      END IF
    END IF

    K = ( ( - ( Y / Tau ) * dPdDe + dPdE * ( &
          Half * Vsq + Em ) + dPdTau * Tau ) / ( dPdE ) )
    H = ( Cs**2 / ( dPdE * Tau ) ) + K
    W = Tau * ( dPdE * ( Vsq - 2.0_DP * Em ) &
               - 2.0_DP * dPdTau * Tau )

    ! --- Compute the flux Jacobian for debugging/use in numeric routine ---

    SELECT CASE( iDim )
      CASE(1)
        CALL ComputeFluxJacobian_X1 &
               ( Gmdd11, Tau, T, Y, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, Vsq, Em, H, dPdTau, dPdE, dPdDe, dFdU )
      CASE(2)
        CALL ComputeFluxJacobian_X2 &
               ( Gmdd22, Tau, T, Y, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, Vsq, Em, H, dPdTau, dPdE, dPdDe, dFdU )
      CASE(3)
        CALL ComputeFluxJacobian_X3 &
               ( Gmdd33, Tau, T, Y, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, Vsq, Em, H, dPdTau, dPdE, dPdDe, dFdU )
    END SELECT

    IF( PRESENT( FJ ) ) THEN
      FJ = dFdU
    END IF

    IF ( UseAnalytic ) THEN
      CALL ComputeCharacteristicDecomposition_Analytic &
             ( iDim, Gmdd11, Gmdd22, Gmdd33, &
               D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, E, &
               Ne, P, Cs, K, H, Tau, T, Y, VSq, W, Em, &
               dPdD, dPdT, dPdY, dEdD, dEdT, dEdY, dPdE, &
               dPdDe, dPdTau, R, invR )

    ELSE
      CALL ComputeCharacteristicDecomposition_Numeric( R, invR, dFdU )
    END IF

    CALL TimersStop_Euler( Timer_Euler_SL_CharDecomp )

    ! -- Begin debugging statements. ---

    IF ( PRESENT( cs_T ) )THEN
      cs_T = Cs_table
    END IF

  END SUBROUTINE ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE

  SUBROUTINE ComputeCharacteristicDecomposition_Analytic &
               ( iDim, Gmdd11, Gmdd22, Gmdd33, &
                 D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3, E, &
                 Ne, P, Cs, K, H, Tau, T, Y, VSq, W, Em, &
                 dPdD, dPdT, dPdY, dEdD, dEdT, dEdY, dPdE, &
                 dPdDe, dPdTau, R, invR )

    INTEGER,  INTENT(in)  :: iDim

    REAL(DP), INTENT(in) :: Gmdd11, Gmdd22, Gmdd33

    REAL(DP), INTENT(in) :: D, Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
    REAL(DP), INTENT(in) :: E, Ne, P, Cs
    REAL(DP), INTENT(in) :: K, H, Tau, T, Y, Vsq, W, Em

    REAL(DP), INTENT(in) :: dPdD, dPdT, dPdY
    REAL(DP), INTENT(in) :: dEdD, dEdT, dEdY
    REAL(DP), INTENT(in) :: dPdE, dPdDe, dPdTau

    INTEGER :: i

    REAL(DP) :: X, Alpha, B, Delta, Zero2, invCsSq
    REAL(DP), DIMENSION(3) :: Phi_u, Phi_d

    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

    Phi_u(1) = dPdE * Tau * Vu1
    Phi_u(2) = dPdE * Tau * Vu2
    Phi_u(3) = dPdE * Tau * Vu3

    Phi_d(1) = dPdE * Tau * Vd1
    Phi_d(2) = dPdE * Tau * Vd2
    Phi_d(3) = dPdE * Tau * Vd3

    invCsSq = 1.0_DP / Cs**2

    SELECT CASE( iDim )

      CASE( 1 )

        Delta = Vu1 * Vd1 - Vu2 * Vd2 - Vu3 * Vd3

        B = 0.5_DP * (Delta + 2.0_DP * Em + &
                     (2.0_DP * dPdTau * Tau)/dPdE)
        X = (dPdE * ( Delta + 2.0_DP * Em) + 2.0_DP * dPdTau * Tau )

        Alpha = 2.0_DP * (Y) * dPdDe - X * Tau

        R(:,1) = [ One, Vd1 - Cs * SQRT( Gmdd11 ), Vd2, &
                   Vd3, H - Cs * SQRT( Gmdd11 ) * Vu1, Y ]
        R(:,2) = [ Zero, Zero, One, Zero, Vu2, Zero ]
        R(:,3) = [ One, Vd1, Zero, Zero, B, Zero ]
        R(:,4) = [ One, Vd1, Zero, Zero, Zero, &
                   (Tau * X) / (2.0_DP * dPdDe) ]
        R(:,5) = [ Zero, Zero, Zero, One, Vu3, Zero ]
        R(:,6) = [ One, Vd1 + Cs * SQRT( Gmdd11 ), Vd2, &
                   Vd3, H + Cs * SQRT( Gmdd11 ) * Vu1, Y ]

        invR(:,1) = invCsSq * &
            [ + 0.25_DP * (W + 2.0_DP * Cs * SQRT( Gmdd11 ) * Vu1), &
              - Half * Vd2 * W, &
              + (2.0_DP * Cs**2 * X + Alpha * W / Tau)/(2.0_DP * X), &
              - (Y) * dPdDe * W / (X * Tau), &
              - Half * Vd3 * W, &
              + 0.25_DP * (W - 2.0_DP * Cs * SQRT( Gmdd11 ) * Vu1) ]

        invR(:,2) = invCsSq * &
            [ - Half * ( ( Cs / SQRT( Gmdd11 ) ) + Phi_u(1) ), &
              + Phi_u(1) * Vd2, &
              - Phi_u(1) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi_u(1) / (X * Tau), &
              + Phi_u(1) * Vd3, &
              + Half * ( ( Cs / SQRT( Gmdd11 ) ) - Phi_u(1) ) ]

        invR(:,3) = invCsSq * &
            [ - Half * Phi_u(2), &
              + Cs**2 + Phi_u(2) * Vd2, &
              - Phi_u(2) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi_u(2) / (X * Tau), &
              + Phi_u(2) * Vd3, &
              - Half * Phi_u(2) ]

        invR(:,4) = invCsSq * &
            [ - Half * Phi_u(3), &
              + Phi_u(3) * Vd2, &
              - Phi_u(3) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi_u(3) / (X * Tau), &
              + Cs**2 + Phi_u(3) * Vd3, &
              - Half * Phi_u(3) ]

        invR(:,5) = invCsSq * &
            [ + Half * dPdE * Tau, &
              - Phi_d(2), &
              + dPdE * Alpha  / X, &
              - ((2.0_DP * Y * dPdDe * dPdE)/ X), &
              - Phi_d(3), &
              + Half * dPdE * Tau ]

        invR(:,6) = invCsSq * &
            [ + Half * dPdDe, &
              - Vd2 * dPdDe, &
              + (dPdDe * (-2.0_DP * Cs**2 + Alpha))/(Tau * X), &
              + 2.0_DP * dPdDe * (Cs**2 - Y * dPdDe)/(Tau * X), &
              - Vd3 * dPdDe, &
              + Half * dPdDe ]

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X1):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

      CASE(2)

        Delta = Vu2 * Vd2 - Vu1 * Vd1 - Vu3 * Vd3

        B = 0.5_DP * (Delta + 2.0_DP * Em + &
                     (2.0_DP * dPdTau * Tau)/dPdE)
        X = (dPdE * ( Delta + 2.0_DP * Em) + 2.0_DP * dPdTau * Tau )

        Alpha = 2.0_DP * (Y) * dPdDe - X * Tau

        R(:,1) = [ One, Vd1, Vd2 - Cs * SQRT( Gmdd22 ), &
                   Vd3, H - Cs * SQRT( Gmdd22 ) * Vu2, Y ]
        R(:,2) = [ Zero, One, Zero, Zero, Vu1, Zero ]
        R(:,3) = [ One, Zero, Vd2, Zero, B, Zero ]
        R(:,4) = [ One, Zero, Vd2, Zero, Zero, &
                   (Tau * X) / (2.0_DP * dPdDe) ]
        R(:,5) = [ Zero, Zero, Zero, One, Vu3, Zero ]
        R(:,6) = [ One, Vd1, Vd2 + Cs * SQRT( Gmdd22 ), &
                   Vd3, H + Cs * SQRT( Gmdd22 ) * Vu2, Y ]

        invR(:,1) = invCsSq * &
            [ + 0.25_DP * (W + 2.0_DP * Cs * SQRT( Gmdd22 ) * Vu2), &
              - Half * Vd1 * W, &
              + (2.0_DP * Cs**2 * X + Alpha * W / Tau)/(2.0_DP * X), &
              - (Y) * dPdDe * W / (X * Tau), &
              - Half * Vd3 * W, &
              + 0.25_DP * (W - 2.0_DP * Cs * SQRT( Gmdd22 ) * Vu2) ]

        invR(:,2) = invCsSq * &
            [ - Half * (Phi_u(1)), &
              + Cs**2 + Phi_u(1) * Vd1, &
              - Phi_u(1) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi_u(1) / (X * Tau), &
              + Phi_u(1) * Vd3, &
              - Half * (Phi_u(1)) ]

        invR(:,3) = invCsSq * &
            [ - Half * ( ( Cs / Gmdd22 ) + Phi_u(2)), &
              + Phi_u(2) * Vd1, &
              - Phi_u(2) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi_u(2) / (X * Tau), &
              + Phi_u(2) * Vd3, &
              + Half * ( ( Cs / Gmdd22 ) - Phi_u(2)) ]

        invR(:,4) = invCsSq * &
            [ - Half * Phi_u(3), &
              + Phi_u(3) * Vd1, &
              - Phi_u(3) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi_u(3) / (X * Tau), &
              + Cs**2 + Phi_u(3) * Vd3, &
              - Half * Phi_u(3) ]

        invR(:,5) = invCsSq * &
            [ + Half * dPdE * Tau, &
              - Phi_d(1), &
              + dPdE * Alpha  / X, &
              - ((2.0_DP * Y * dPdDe * dPdE)/ X), &
              - Phi_d(3), &
              + Half * dPdE * Tau ]

        invR(:,6) = invCsSq * &
            [ + Half * dPdDe, &
              - Vd1 * dPdDe, &
              + (dPdDe * (-2.0_DP * Cs**2 + Alpha))/(Tau * X), &
              + 2.0_DP * dPdDe * (Cs**2 - Y * dPdDe)/(Tau * X), &
              - Vd3 * dPdDe, &
              + Half * dPdDe ]

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X2):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition_Analytic


  SUBROUTINE ComputeCharacteristicDecomposition_Numeric( R, invR, dFdU )

    REAL(DP), INTENT(in)  :: dFdU(nCF,nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)
    REAL(DP)              :: Lambda(nCF,nCF)

    INTEGER  :: i, INFO, LWORK
    INTEGER  :: IPIV(nCF)

    REAL(DP)              :: WR(nCF)
    REAL(DP)              :: WI(nCF)
    REAL(DP)              :: TEMP(1)
    REAL(DP)              :: dFdU_Copy(nCF,nCF)
    REAL(DP)              :: invR_Copy(nCF,nCF)
    REAL(DP), ALLOCATABLE :: WORK1(:), WORK2(:)

    ! --- Copy to avoid overwriting dFdU ---

    dFdU_Copy = dFdU

    ! --- Necessary workplace query to get LWORK. ----

    CALL DGEEV( 'V', 'N', nCF, dFdU_Copy, nCF, WR, &
                WI, invR, nCF, 0, nCF, TEMP, &
                -1, INFO )

    LWORK = TEMP(1)
    ALLOCATE( WORK1(LWORK) )

    Lambda(:,:) = Zero

    ! --- Actual computation of eigedecomposition. ---

    CALL DGEEV( 'V', 'N', nCF, dFdU_Copy, nCF, WR, &
                WI, invR, nCF, 0, nCF, WORK1,  &
                LWORK, INFO )

    invR = TRANSPOSE( invR )

    invR_Copy = invR

    CALL DGETRF( nCF, nCF, invR_Copy, nCF, IPIV, INFO )

    LWORK = -1
    CALL DGETRI( nCF, invR_Copy, nCF, IPIV, TEMP, LWORK, INFO )

    LWORK = TEMP(1)
    ALLOCATE( WORK2(LWORK) )

    CALL DGETRI( nCF, invR_Copy, nCF, IPIV, WORK2, LWORK, INFO )

    R = invR_Copy

    IF ( ( INFO .NE. 0 ) .OR. ( ANY( ABS( WI ) > 1d-15 ) ) ) THEN

      PRINT*, 'INFO: ', INFO
      PRINT*, 'WR: ', WR
      PRINT*, 'WI: ', WI

      DO i = 1, 6

        PRINT*, 'Lambda(i,:) : ', Lambda(i,:)

      END DO

    END IF

  END SUBROUTINE ComputeCharacteristicDecomposition_Numeric

  SUBROUTINE ComputeFluxJacobian_X1( Gmdd11, Tau, T, Y, Vu1, Vu2, Vu3, &
                                     Vd1, Vd2, Vd3, Vsq, Em, H, dPdTau, &
                                     dPdE, dPdDe, dFdU_X1 )

    REAL(DP),  INTENT(in) :: Gmdd11
    REAL(DP),  INTENT(in) :: Tau, T, Y
    REAL(DP),  INTENT(in) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
    REAL(DP),  INTENT(in) :: Vsq, Em, H
    REAL(DP),  INTENT(in) :: dPdTau, dPdE, dPdDe
    REAL(DP), INTENT(out) :: dFdU_X1(nCF,nCF)

    dFdU_X1(1,:) = [ 0.0_DP, 1.0_DP / Gmdd11, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X1(2,1) = - Vu1 * Vd1 - Tau**2 * dPdTau &
                   - Tau * dPdE * ( Em - 0.5_DP * Vsq )
    dFdU_X1(2,2) = Vu1 * ( 2.0_DP - Tau * dPdE )
    dFdU_X1(2,3) = - dPdE * Tau * Vu2
    dFdU_X1(2,4) = - dPdE * Tau * Vu3
    dFdU_X1(2,5) = dPdE * Tau
    dFdU_X1(2,6) = dPdDe

    dFdU_X1(3,:) = [ - Vu1 * Vd2, Vd2 / Gmdd11, Vu1, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X1(4,:) = [ - Vu1 * Vd3, Vd3 / Gmdd11, 0.0_DP, Vu1, 0.0_DP, 0.0_DP ]

    dFdU_X1(5,1) = Vu1 * ( - H - dPdTau * Tau**2 &
                          - Tau * dPdE * ( Em - 0.5_DP * Vsq ) )
    dFdU_X1(5,2) = ( H / Gmdd11 ) - dPdE * Tau * Vu1**2
    dFdU_X1(5,3) =   - dPdE * Tau * Vu1 * Vu2
    dFdU_X1(5,4) =   - dPdE * Tau * Vu1 * Vu3
    dFdU_X1(5,5) = Vu1 * ( 1.0_DP + dPdE * Tau )
    dFdU_X1(5,6) = Vu1 * dPdDe

    dFdU_X1(6,:) = [ - Vu1 * Y, Y / Gmdd11, 0.0_DP, 0.0_DP, 0.0_DP, Vu1 ]

  END SUBROUTINE ComputeFluxJacobian_X1

  SUBROUTINE ComputeFluxJacobian_X2( Gmdd22, Tau, T, Y, Vu1, Vu2, Vu3, &
                                     Vd1, Vd2, Vd3, Vsq, Em, H, dPdTau, &
                                     dPdE, dPdDe, dFdU_X2 )

    REAL(DP),  INTENT(in) :: Gmdd22
    REAL(DP),  INTENT(in) :: Tau, T, Y
    REAL(DP),  INTENT(in) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
    REAL(DP),  INTENT(in) :: Vsq, Em, H
    REAL(DP),  INTENT(in) :: dPdTau, dPdE, dPdDe
    REAL(DP), INTENT(out) :: dFdU_X2(nCF,nCF)

    dFdU_X2(1,:) = [ 0.0_DP, 0.0_DP, 1.0_DP / Gmdd22, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X2(2,:) = [ - Vu2 * Vd1, Vu2, Vd1 / Gmdd22, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X2(3,1) = - Vu2 * Vd2 - Tau**2 * dPdTau &
                   - Tau * dPdE * ( Tau * Em - 0.5_DP * Vsq )
    dFdU_X2(3,2) = - dPdE * Tau * Vu1
    dFdU_X2(3,3) = Vu2 * ( 2.0_DP - Tau * dPdE )
    dFdU_X2(3,4) = - dPdE * Tau * Vu3
    dFdU_X2(3,5) = dPdE * Tau
    dFdU_X2(3,6) = dPdDe

    dFdU_X2(4,:) = [ - Vu2 * Vd3, 0.0_DP, Vd3 / Gmdd22, Vu2, 0.0_DP, 0.0_DP ]

    dFdU_X2(5,1) = Vu2 * ( - H - dPdTau * Tau**2 &
                          - Tau * dPdE * ( Em - 0.5_DP * Vsq ) )
    dFdU_X2(5,2) =   - dPdE * Tau * Vu1 * Vu2
    dFdU_X2(5,3) = ( H / Gmdd22 ) - dPdE * Tau * Vu2**2
    dFdU_X2(5,4) =   - dPdE * Tau * Vu2 * Vu3
    dFdU_X2(5,5) = Vu2 * ( 1.0_DP + dPdE * Tau )
    dFdU_X2(5,6) = Vu2 * dPdDe

    dFdU_X2(6,:) = [ - Vu2 * Y, 0.0_DP, Y / Gmdd22, 0.0_DP, 0.0_DP, Vu2 ]

  END SUBROUTINE ComputeFluxJacobian_X2

  SUBROUTINE ComputeFluxJacobian_X3( Gmdd33, Tau, T, Y, Vu1, Vu2, Vu3, &
                                     Vd1, Vd2, Vd3, Vsq, Em, H, dPdTau, &
                                     dPdE, dPdDe, dFdU_X3 )

    REAL(DP),  INTENT(in) :: Gmdd33
    REAL(DP),  INTENT(in) :: Tau, T, Y
    REAL(DP),  INTENT(in) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
    REAL(DP),  INTENT(in) :: Vsq, Em, H
    REAL(DP),  INTENT(in) :: dPdTau, dPdE, dPdDe
    REAL(DP), INTENT(out) :: dFdU_X3(nCF,nCF)

    dFdU_X3(1,:) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP / Gmdd33, 0.0_DP, 0.0_DP ]

    dFdU_X3(2,:) = [ - Vu3 * Vd1, Vu3, 0.0_DP, Vd1 / Gmdd33, 0.0_DP, 0.0_DP ]

    dFdU_X3(3,:) = [ - Vu3 * Vd2, 0.0_DP, Vu3, Vd2 / Gmdd33, 0.0_DP, 0.0_DP ]

    dFdU_X3(4,1) = - Vu3 * Vd3 - Tau**2 * dPdTau &
                   - Tau * dPdE * ( Tau * Em - 0.5_DP * Vsq )
    dFdU_X3(4,2) = - dPdE * Tau * Vu1
    dFdU_X3(4,3) = - dPdE * Tau * Vu2
    dFdU_X3(4,4) = Vu3 * ( 2.0_DP - Tau * dPdE )
    dFdU_X3(4,5) = dPdE * Tau
    dFdU_X3(4,6) = dPdDe

    dFdU_X3(5,1) = Vu3 * ( - H - dPdTau * Tau**2 &
                          - Tau * dPdE * ( Em - 0.5_DP * Vsq ) )
    dFdU_X3(5,2) =   - dPdE * Tau * Vu3 * Vu1
    dFdU_X3(5,3) = ( H / Gmdd33 ) - dPdE * Tau * Vu3**2
    dFdU_X3(5,4) =   - dPdE * Tau * Vu3 * Vu2
    dFdU_X3(5,5) = Vu3 * ( 1.0_DP + dPdE * Tau )
    dFdU_X3(5,6) = Vu3 * dPdDe

    dFdU_X3(6,:) = [ - Vu3 * Y, 0.0_DP, 0.0_DP, Y / Gmdd33, 0.0_DP, Vu3 ]

  END SUBROUTINE ComputeFluxJacobian_X3

END MODULE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE

