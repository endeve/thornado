MODULE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE UnitsModule, ONLY: &
    AtomicMassUnit, Gram, Centimeter, MeV
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, TimersStop_Euler, &
    Timer_Euler_CharacteristicDecomposition

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: Debug = .FALSE.

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

    LOGICAL :: UseAnalytic

    INTEGER  :: i
    REAL(DP) :: D, V1, V2, V3, E, Ne, P, Cs, Cs_table
    REAL(DP) :: K, H, Tau, T, Y, Vsq, W, Em, Gm, S

    REAL(DP) :: dPdD, dPdT, dPdY
    REAL(DP) :: dEdD, dEdT, dEdY
    REAL(DP) :: dPdE, dPdDe, dPdTau

    REAL(DP) :: X, Alpha, B, Delta, Zero2
    REAL(DP), DIMENSION(3) :: Phi

    REAL(DP) :: dFdU(nCF,nCF)

    IF ( PRESENT ( UseAnalytic_Option ) ) THEN
      UseAnalytic = UseAnalytic_Option
    ELSE
      UseAnalytic = .FALSE.
    END IF

    CALL TimersStart_Euler( Timer_Euler_CharacteristicDecomposition )

    CALL ComputePrimitive_Euler_NonRelativistic &
           ( U(iCF_D ), U(iCF_S1), U(iCF_S2), &
             U(iCF_S3), U(iCF_E ), U(iCF_Ne), &
             D, V1, V2, V3, E, Ne, &
             G(iGF_Gm_dd_11), &
             G(iGF_Gm_dd_22), &
             G(iGF_Gm_dd_33) )

    CALL ComputeAuxiliary_Fluid_TABLE &
          ( D, E, Ne, P, T, Y, S, Em, Gm, Cs_table )

    CALL ComputeSpecificInternalEnergy_TABLE &
          ( D, T, Y, Em, dEdD, dEdT, dEdY )

    CALL ComputePressure_TABLE &
          ( D, T, Y, P, dPdD, dPdT, dPdY )

    Tau = 1.0_DP / D

    dPdE   = dPdT / dEdT
    dPdDe  = ( Tau ) * ( dPdY - dEdY * dPdE )
    dPdTau = (dPdDe * Y + dEdD * dPdE - dPdD) / (Tau**2)

    Vsq = V1**2 + V2**2 + V3**2

    IF ( Tau**2 * ( P * dPdE - dPdTau ) + Y * dPdDe .GT. Zero) THEN
      ! Cs = SQRT( Tau**2 * ( P * dPdE - dPdTau) + Y * dPdDe   )
      Cs = Cs_table
    ELSE
      Cs = Cs_table
      ! write(*,*) "Switching Sound Speed to TABLE value"
    END IF

    K = ( (- (Y/Tau) * dPdDe + dPdE * ( &
          Half * Vsq + Em) + dPdTau * Tau )/(dPdE) )
    H = ( Cs**2 / (dPdE * Tau) ) + K
    W = Tau * (dPdE*(Vsq - 2.0_DP * Em) &
               - 2.0_DP * dPdTau * Tau )

    ! --- Compute the flux Jacobian for debugging/use in numeric routine ---

    SELECT CASE( iDim )
      CASE(1)
        CALL ComputeFluxJacobian_X1( Tau, T, Y, V1, V2, V3, Vsq, Em, H, dPdTau, dPdE, dPdDe, dFdU )
      CASE(2)
        CALL ComputeFluxJacobian_X2( Tau, T, Y, V1, V2, V3, Vsq, Em, H, dPdTau, dPdE, dPdDe, dFdU )
      CASE(3)
        CALL ComputeFluxJacobian_X3( Tau, T, Y, V1, V2, V3, Vsq, Em, H, dPdTau, dPdE, dPdDe, dFdU )
    END SELECT

    IF( PRESENT( FJ ) ) THEN
      FJ = dFdU
    END IF

    IF ( UseAnalytic ) THEN
      CALL ComputeCharacteristicDecomposition_Analytic( iDim, D, V1, V2, V3, E, Ne, P, Cs,  &
                                                        K, H, Tau, T, Y, VSq, W, Em, dPdD,  &
                                                        dPdT, dPdY, dEdD, dEdT, dEdY, dPdE, &
                                                        dPdDe, dPdTau, R, invR )
    ELSE
      CALL ComputeCharacteristicDecomposition_Numeric( R, invR, dFdU )
    END IF

    CALL TimersStop_Euler( Timer_Euler_CharacteristicDecomposition )

    ! -- Begin debugging statements. ---

    IF ( PRESENT( cs_T ) )THEN
      cs_T = Cs
    END IF

  END SUBROUTINE ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE

  SUBROUTINE ComputeCharacteristicDecomposition_Analytic &
               ( iDim, D, V1, V2, V3, E, Ne, P, Cs,  &
                 K, H, Tau, T, Y, VSq, W, Em, dPdD,  &
                 dPdT, dPdY, dEdD, dEdT, dEdY, dPdE, &
                 dPdDe, dPdTau, R, invR )

    INTEGER,  INTENT(in)  :: iDim

    REAL(DP), INTENT(in) :: D, V1, V2, V3, E, Ne, P, Cs
    REAL(DP), INTENT(in) :: K, H, Tau, T, Y, Vsq, W, Em

    REAL(DP), INTENT(in) :: dPdD, dPdT, dPdY
    REAL(DP), INTENT(in) :: dEdD, dEdT, dEdY
    REAL(DP), INTENT(in) :: dPdE, dPdDe, dPdTau

    INTEGER :: i

    REAL(DP) :: X, Alpha, B, Delta, Zero2, invCsSq
    REAL(DP), DIMENSION(3) :: Phi

    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

    Phi(1) = dPdE * Tau * V1
    Phi(2) = dPdE * Tau * V2
    Phi(3) = dPdE * Tau * V3

    invCsSq = 1.0_DP / Cs**2

    SELECT CASE( iDim )

      CASE( 1 )

        Delta = V1**2 - V2**2 - V3**2

        B = 0.5_DP * (Delta + 2.0_DP * Em + &
                     (2.0_DP * dPdTau * Tau)/dPdE)
        X = (dPdE * ( Delta + 2.0_DP * Em) + 2.0_DP * dPdTau * Tau )

        Alpha = 2.0_DP * (Y) * dPdDe - X * Tau

        R(:,1) = [ One, V1 - Cs, V2, V3, H - Cs * V1, Y]
        R(:,2) = [ Zero, Zero, One, Zero, V2, Zero ]
        R(:,3) = [ One, V1, Zero, Zero, B, Zero ]
        R(:,4) = [ One, V1, Zero, Zero, Zero, (Tau * X) / (2.0_DP * dPdDe) ]
        R(:,5) = [ Zero, Zero, Zero, One, V3, Zero ]
        R(:,6) = [ One, V1 + Cs, V2, V3, H + Cs * V1, Y]

        invR(:,1) = invCsSq * &
            [ + (2.0_DP * Cs * V1 + W) * 0.25_DP, &
              - Half * V2 * W, &
              + (2.0_DP * Cs**2 * X + Alpha * W / Tau)/(2.0_DP * X), &
              - (Y) * dPdDe * W / (X * Tau), &
              - Half * V3 * W, &
              + 0.25_DP * (W - 2.0_DP * Cs * V1) ]

        invR(:,2) = invCsSq * &
            [ - Half * (Cs + Phi(1)), &
              + Phi(1) * V2, &
              - Phi(1) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(1) / (X * Tau), &
              + Phi(1) * V3, &
              + Half * (Cs - Phi(1)) ]

        invR(:,3) = invCsSq * &
            [ - Half * Phi(2), &
              + Cs**2 + Phi(2) * V2, &
              - Phi(2) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(2) / (X * Tau), &
              + Phi(2) * V3, &
              - Half * Phi(2) ]

        invR(:,4) = invCsSq * &
            [ - Half * Phi(3), &
              + Phi(3) * V2, &
              - Phi(3) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(3) / (X * Tau), &
              + Cs**2 + Phi(3) * V3, &
              - Half * Phi(3) ]

        invR(:,5) = invCsSq * &
            [ + Half * dPdE * Tau, &
              - Phi(2), &
              + dPdE * Alpha  / X, &
              - ((2.0_DP * Y * dPdDe * dPdE)/ X), &
              - Phi(3), &
              + Half * dPdE * Tau ]

        invR(:,6) = invCsSq * &
            [ + Half * dPdDe, &
              - V2 * dPdDe, &
              + (dPdDe * (-2.0_DP * Cs**2 + Alpha))/(Tau * X), &
              + 2.0_DP * dPdDe * (Cs**2 - Y * dPdDe)/(Tau * X), &
              - dPdDe * V3, &
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

        Delta = V2**2 - V3**2 - V1**2

        B = 0.5_DP * (Delta + 2.0_DP * Em + &
                     (2.0_DP * dPdTau * Tau)/dPdE)
        X = (dPdE * ( Delta + 2.0_DP * Em) + 2.0_DP * dPdTau * Tau )

        Alpha = 2.0_DP * (Y) * dPdDe - X * Tau

        R(:,1) = [ One, V1, V2 - Cs, V3, H - Cs * V2, Y]
        R(:,2) = [ Zero, One, Zero, Zero, V1, Zero ]
        R(:,3) = [ One, Zero, V2, Zero, B, Zero ]
        R(:,4) = [ One, Zero, V2, Zero, Zero, (Tau * X) / (2.0_DP * dPdDe) ]
        R(:,5) = [ Zero, Zero, Zero, One, V3, Zero ]
        R(:,6) = [ One, V1, V2 + Cs, V3, H + Cs * V2, Y]

        invR(:,1) = invCsSq * &
            [ + (2.0_DP * Cs * V2 + W) * 0.25_DP, &
              - Half * V1 * W, & ! CHANGE V2 TO V1
              + (2.0_DP * Cs**2 * X + Alpha * W / Tau)/(2.0_DP * X), &
              - (Y) * dPdDe * W / (X * Tau), &
              - Half * V3 * W, &
              + 0.25_DP * (W - 2.0_DP * Cs * V2) ] ! CHANGE TO V2

        invR(:,2) = invCsSq * &
            [ - Half * (Phi(1)), &
              + Cs**2 + Phi(1) * V1, &
              - Phi(1) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(1) / (X * Tau), &
              + Phi(1) * V3, &
              - Half * (Phi(1)) ]

        invR(:,3) = invCsSq * &
            [ - Half * (Cs + Phi(2)), &
              + Phi(2) * V1, &
              - Phi(2) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(2) / (X * Tau), &
              + Phi(2) * V3, &
              + Half * (Cs - Phi(2)) ]

        invR(:,4) = invCsSq * &
            [ - Half * Phi(3), &
              + Phi(3) * V1, &
              - Phi(3) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(3) / (X * Tau), &
              + Cs**2 + Phi(3) * V3, &
              - Half * Phi(3) ]

        invR(:,5) = invCsSq * &
            [ + Half * dPdE * Tau, &
              - Phi(1), &
              + dPdE * Alpha  / X, &
              - ((2.0_DP * Y * dPdDe * dPdE)/ X), &
              - Phi(3), &
              + Half * dPdE * Tau ]

        invR(:,6) = invCsSq * &
            [ + Half * dPdDe, &
              - V1 * dPdDe, &
              + (dPdDe * (-2.0_DP * Cs**2 + Alpha))/(Tau * X), &
              + 2.0_DP * dPdDe * (Cs**2 - Y * dPdDe)/(Tau * X), &
              - dPdDe * V3, &
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

  SUBROUTINE ComputeFluxJacobian_X1( Tau, T, Y, V1, V2, V3, Vsq, Em, H, &
                                     dPdTau, dPdE, dPdDe, dFdU_X1 )

    REAL(DP),  INTENT(in) :: Tau, T, Y, V1, V2, V3, VSq, Em, H
    REAL(DP),  INTENT(in) :: dPdTau, dPdE, dPdDe
    REAL(DP), INTENT(out) :: dFdU_X1(nCF,nCF)

    dFdU_X1(1,:) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X1(2,1) = - V1**2 - Tau**2 * dPdTau &
                   - Tau * dPdE * ( Em - 0.5_DP * Vsq )
    dFdU_X1(2,2) = V1 * ( 2.0_DP - Tau * dPdE )
    dFdU_X1(2,3) = - dPdE * Tau * V2
    dFdU_X1(2,4) = - dPdE * Tau * V3
    dFdU_X1(2,5) = dPdE * Tau
    dFdU_X1(2,6) = dPdDe

    dFdU_X1(3,:) = [ - V1 * V2, V2, V1, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X1(4,:) = [ - V1 * V3, V3, 0.0_DP, V1, 0.0_DP, 0.0_DP ]

    dFdU_X1(5,1) = V1 * ( - H - dPdTau * Tau**2 &
                          - Tau * dPdE * ( Em - 0.5_DP * Vsq ) )
    dFdU_X1(5,2) = H - dPdE * Tau * V1**2
    dFdU_X1(5,3) =   - dPdE * Tau * V1 * V2
    dFdU_X1(5,4) =   - dPdE * Tau * V1 * V3
    dFdU_X1(5,5) = V1 * ( 1.0_DP + dPdE * Tau )
    dFdU_X1(5,6) = V1 * dPdDe

    dFdU_X1(6,:) = [ - V1 * Y, Y, 0.0_DP, 0.0_DP, 0.0_DP, V1 ]

  END SUBROUTINE ComputeFluxJacobian_X1

  SUBROUTINE ComputeFluxJacobian_X2( Tau, T, Y, V1, V2, V3, Vsq, E, H, &
                                     dPdTau, dPdE, dPdDe, dFdU_X2 )

    REAL(DP),  INTENT(in) :: Tau, T, Y, V1, V2, V3, VSq, E, H
    REAL(DP),  INTENT(in) :: dPdTau, dPdE, dPdDe
    REAL(DP), INTENT(out) :: dFdU_X2(nCF,nCF)

    dFdU_X2(1,:) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X2(2,:) = [ - V1 * V2, V2, V1, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X2(3,1) = - V2**2 - Tau**2 * dPdTau &
                   - Tau * dPdE * ( Tau * E - 0.5_DP * Vsq )
    dFdU_X2(3,2) = - dPdE * Tau * V1
    dFdU_X2(3,3) = V2 * ( 2.0_DP - Tau * dPdE )
    dFdU_X2(3,4) = - dPdE * Tau * V3
    dFdU_X2(3,5) = dPdE * Tau
    dFdU_X2(3,6) = dPdDe

    dFdU_X2(4,:) = [ - V2 * V3, 0.0_DP, V3, V2, 0.0_DP, 0.0_DP ]

    dFdU_X2(5,1) = V2 * ( - H - dPdTau * Tau**2 &
                          - Tau * dPdE * ( Tau * E - 0.5_DP * Vsq ) )
    dFdU_X2(5,2) =   - dPdE * Tau * V1 * V2
    dFdU_X2(5,3) = H - dPdE * Tau * V2**2
    dFdU_X2(5,4) =   - dPdE * Tau * V2 * V3
    dFdU_X2(5,5) = V1 * ( 1.0_DP + dPdE * Tau )
    dFdU_X2(5,6) = V1 * dPdDe

    dFdU_X2(6,:) = [ - V2 * Y, 0.0_DP, Y, 0.0_DP, 0.0_DP, V2 ]

  END SUBROUTINE ComputeFluxJacobian_X2

  SUBROUTINE ComputeFluxJacobian_X3( Tau, T, Y, V1, V2, V3, Vsq, E, H, &
                                     dPdTau, dPdE, dPdDe, dFdU_X3 )

    REAL(DP),  INTENT(in) :: Tau, T, Y, V1, V2, V3, VSq, E, H
    REAL(DP),  INTENT(in) :: dPdTau, dPdE, dPdDe
    REAL(DP), INTENT(out) :: dFdU_X3(nCF,nCF)

    dFdU_X3(1,:) = [ 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]

    dFdU_X3(2,:) = [ - V1 * V3, V3, 0.0_DP, V1, 0.0_DP, 0.0_DP ]

    dFdU_X3(3,:) = [ - V3 * V2, 0.0_DP, V3, V2, 0.0_DP, 0.0_DP ]

    dFdU_X3(4,1) = - V3**2 - Tau**2 * dPdTau &
                   - Tau * dPdE * ( Tau * E - 0.5_DP * Vsq )
    dFdU_X3(4,2) = - dPdE * Tau * V1
    dFdU_X3(4,3) = - dPdE * Tau * V2
    dFdU_X3(4,4) = V3 * ( 2.0_DP - Tau * dPdE )
    dFdU_X3(4,5) = dPdE * Tau
    dFdU_X3(4,6) = dPdDe

    dFdU_X3(5,1) = V2 * ( - H - dPdTau * Tau**2 &
                          - Tau * dPdE * ( Tau * E - 0.5_DP * Vsq ) )
    dFdU_X3(5,2) =   - dPdE * Tau * V1 * V2
    dFdU_X3(5,3) = H - dPdE * Tau * V2**2
    dFdU_X3(5,4) =   - dPdE * Tau * V2 * V3
    dFdU_X3(5,5) = V1 * ( 1.0_DP + dPdE * Tau )
    dFdU_X3(5,6) = V1 * dPdDe

    dFdU_X3(6,:) = [ - V2 * Y, 0.0_DP, Y, 0.0_DP, 0.0_DP, V2 ]

  END SUBROUTINE ComputeFluxJacobian_X3

END MODULE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE

