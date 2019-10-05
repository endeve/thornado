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
    Euler_ComputePrimitive_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: Debug = .FALSE.

  PUBLIC :: ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
    ( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

    INTEGER :: i
    REAL(DP), DIMENSION(1) :: D, V1, V2, V3, E, Ne, P, Cs, invCsSq, Cs_table
    REAL(DP), DIMENSION(1) :: K, H, Tau, T, Y, Vsq, W, Em, Gm, S

    REAL(DP), DIMENSION(1) :: dPdD, dPdT, dPdY
    REAL(DP), DIMENSION(1) :: dEdD, dEdT, dEdY
    REAL(DP), DIMENSION(1) :: dPdE, dPdDe, dPdTau

    REAL(DP), DIMENSION(1) :: X, Alpha, B, Delta, Zero2
    REAL(DP), DIMENSION(3) :: Phi

    CALL Euler_ComputePrimitive_NonRelativistic &
           ( [ U(iCF_D ) ], [ U(iCF_S1) ], [ U(iCF_S2) ], &
             [ U(iCF_S3) ], [ U(iCF_E ) ], [ U(iCF_Ne) ], &
             D, V1, V2, V3, E, Ne, &
             [ G(iGF_Gm_dd_11) ], &
             [ G(iGF_Gm_dd_22) ], &
             [ G(iGF_Gm_dd_33) ] )

    CALL ComputeAuxiliary_Fluid_TABLE &
          ( D, E, Ne, P, T, Y, S, Em, Gm, Cs_table )

    CALL ComputeSpecificInternalEnergy_TABLE &
          ( D, T, Y, Em, dEdD, dEdT, dEdY )

    CALL ComputePressure_TABLE &
          ( D, T, Y, P, dPdD, dPdT, dPdY )

    Tau = 1.0_DP / D

    dPdE   = dPdT / dEdT
    dPdDe  = ( Tau ) * ( dPdY - dEdY * dPdE )
    dPdTau = (dPdDe * (Y - 1) + dEdD * dPdE)/ (Tau**2)


    Vsq = V1**2 + V2**2 + V3**2

    IF ( Tau(1)**2 * ( P(1) * dPdE(1) - dPdTau(1) ) + Y(1) * dPdDe(1) .GT. Zero) THEN
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

    Phi(1) = dPdE(1) * Tau(1) * V1(1)
    Phi(2) = dPdE(1) * Tau(1) * V2(1)
    Phi(3) = dPdE(1) * Tau(1) * V3(1)

    invCsSq = 1.0_DP / Cs**2

    SELECT CASE( iDim )

      CASE( 1 )

        Delta = V1**2 - V2**2 - V3**2

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

        invR(:,1) = invCsSq(1) * &
            [ + (2.0_DP * Cs * V1 + W) * 0.25_DP, &
              - Half * V2 * W, &
              + (2.0_DP * Cs**2 * X + Alpha * W / Tau)/(2.0_DP * X), &
              - (Y) * dPdDe * W / (X * Tau), &
              - Half * V3 * W, &
              + 0.25_DP * (W - 2.0_DP * Cs * V1) ]

        invR(:,2) = invCsSq(1) * &
            [ - Half * (Cs + Phi(1)), &
              + Phi(1) * V2, &
              - Phi(1) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(1) / (X * Tau), &
              + Phi(1) * V3, &
              + Half * (Cs - Phi(1)) ]

        invR(:,3) = invCsSq(1) * &
            [ - Half * Phi(2), &
              + Cs**2 + Phi(2) * V2, &
              - Phi(2) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(2) / (X * Tau), &
              + Phi(2) * V3, &
              - Half * Phi(2) ]

        invR(:,4) = invCsSq(1) * &
            [ - Half * Phi(3), &
              + Phi(3) * V2, &
              - Phi(3) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(3) / (X * Tau), &
              + Cs**2 + Phi(3) * V3, &
              - Half * Phi(3) ]

        invR(:,5) = invCsSq(1) * &
            [ + Half * dPdE * Tau, &
              - Phi(2), &
              + dPdE * Alpha  / X, &
              - ((2.0_DP * Y * dPdDe * dPdE)/ X), &
              - Phi(3), &
              + Half * dPdE * Tau ]

        invR(:,6) = invCsSq(1) * &
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

        invR(:,1) = invCsSq(1) * &
            [ + (2.0_DP * Cs * V2 + W) * 0.25_DP, &
              - Half * V1 * W, & ! CHANGE V2 TO V1
              + (2.0_DP * Cs**2 * X + Alpha * W / Tau)/(2.0_DP * X), &
              - (Y) * dPdDe * W / (X * Tau), &
              - Half * V3 * W, &
              + 0.25_DP * (W - 2.0_DP * Cs * V2) ] ! CHANGE TO V2

        invR(:,2) = invCsSq(1) * &
            [ - Half * (Phi(1)), &
              + Cs**2 + Phi(1) * V1, &
              - Phi(1) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(1) / (X * Tau), &
              + Phi(1) * V3, &
              - Half * (Phi(1)) ]

        invR(:,3) = invCsSq(1) * &
            [ - Half * (Cs + Phi(2)), &
              + Phi(2) * V1, &
              - Phi(2) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(2) / (X * Tau), &
              + Phi(2) * V3, &
              + Half * (Cs - Phi(2)) ]

        invR(:,4) = invCsSq(1) * &
            [ - Half * Phi(3), &
              + Phi(3) * V1, &
              - Phi(3) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdDe * Phi(3) / (X * Tau), &
              + Cs**2 + Phi(3) * V3, &
              - Half * Phi(3) ]

        invR(:,5) = invCsSq(1) * &
            [ + Half * dPdE * Tau, &
              - Phi(1), &
              + dPdE * Alpha  / X, &
              - ((2.0_DP * Y * dPdDe * dPdE)/ X), &
              - Phi(3), &
              + Half * dPdE * Tau ]

        invR(:,6) = invCsSq(1) * &
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

  END SUBROUTINE ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE


END MODULE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE
