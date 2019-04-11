MODULE Euler_CharacteristicDecompositionModule_TABLE

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE UnitsModule, ONLY: &
    AtomicMassUnit
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: Debug = .FALSE.

  PUBLIC :: ComputeCharacteristicDecomposition

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

    INTEGER :: i
    REAL(DP), DIMENSION(1) :: D, V1, V2, V3, E, Ne, P, Cs, Alpha, invCsSq, HJ, Cs2
    REAL(DP), DIMENSION(1) :: K, H, Tau, T, Y, Vsq, B, Delta, X, W, Em, Gm, S

    REAL(DP), DIMENSION(1) :: dPdD, dPdT, dPdY
    REAL(DP), DIMENSION(1) :: dEdD, dEdT, dEdY
    REAL(DP), DIMENSION(1) :: dPdE, dPdN, dPdTau

    REAL(DP), DIMENSION(3) :: Phi

    REAL(DP), DIMENSION(6,6) :: dFdU, DD, Matrix ! DELETE

    CALL ComputePrimitive &
           ( [ U(iCF_D ) ], [ U(iCF_S1) ], [ U(iCF_S2) ], &
             [ U(iCF_S3) ], [ U(iCF_E ) ], [ U(iCF_Ne) ], &
             D, V1, V2, V3, E, Ne, &
             [ G(iGF_Gm_dd_11) ], &
             [ G(iGF_Gm_dd_22) ], &
             [ G(iGF_Gm_dd_33) ] )

    CALL ComputeAuxiliary_Fluid_TABLE &
          ( D, E, Ne, P, T, Y, S, Em, Gm, Cs2 )

    CALL ComputeSpecificInternalEnergy_TABLE &
          ( D, T, Y, Em, dEdD, dEdT, dEdY )

    CALL ComputePressure_TABLE &
          ( D, T, Y, P, dPdD, dPdT, dPdY )

    Tau = 1.0_DP / D

    dPdE   = dPdT / dEdT
    dPdN   =  ( Tau ) * ( dPdY - dEdY * dPdE ) ! Got Rid of * AtomicMassUnit
    !dPdTau = - ( dPdD - Y * dPdN / AtomicMassUnit - dEdD * dPdE ) / Tau**2
    dPdTau = (dPdN * (Y - 1) + dEdD * dPdE)/ (Tau**2)


    Vsq = V1**2 + V2**2 + V3**2
  !  Cs = SQRT(Tau * ( (dPdN) * Ne + P * dPdE * Tau - dPdTau * Tau ))
    Cs = SQRT( Tau**2 * ( P * dPdE - dPdTau) + Y * dPdN   )

    K = ( (- (Y/Tau) * dPdN + dPdE * ( &
          Half * Vsq + Em) + dPdTau * Tau )/(dPdE) )  !CHANGES. Ne -> D Y

    H = ( Cs**2 / (dPdE * Tau) ) + K
    Delta = V1**2 - V2**2 - V3**2
    B = 0.5_DP * (Delta + 2.0_DP * Em + &
                 (2.0_DP * dPdTau * Tau)/dPdE)

    X = (dPdE * ( Delta + 2.0_DP * Em) + 2.0_DP * dPdTau * Tau )

    Alpha = 2.0_DP * (Y) * dPdN - X * Tau !Ne -> Y
    W = Tau * (dPdE*(Vsq - 2.0_DP * Em) &
               - 2.0_DP * dPdTau * Tau )

    Phi(1) = dPdE(1) * Tau(1) * V1(1)
    Phi(2) = dPdE(1) * Tau(1) * V2(1)
    Phi(3) = dPdE(1) * Tau(1) * V3(1)

    invCsSq = 1.0_DP / Cs**2

    ! -------
    ! TEMPORARY
    ! ComputeFluxJacobian
    ! --------

    DD = 0.0d0
    DD(1,1) = V1(1) - Cs(1)
    DD(6,6) = Cs(1) + V1(1)
    DO i = 2, 5
      DD(i,i) = V1(1)
    END DO
    HJ   = Em + 0.5_DP * VSq + P * Tau
    dFdU(1,:) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU(2,1) = - V1(1)**2 - Tau(1)**2 * dPdTau(1) &
                - Tau(1) * dPdE(1) * ( Em(1) -  0.5_DP * VSq(1) )
    dFdU(2,2) = V1(1) * ( 2.0_DP - Tau(1) * dPdE(1) )
    dFdU(2,3) = - dPdE(1) * Tau(1) * V2(1)
    dFdU(2,4) = - dPdE(1) * Tau(1) * V3(1)
    dFdU(2,5) = + dPdE(1) * Tau(1)
    dFdU(2,6) = + dPdN(1)

    dFdU(3,:) = [ - V1(1) * V2(1), V2(1), V1(1), 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU(4,:) = [ - V1(1) * V3(1), V3(1), 0.0_DP, V1(1), 0.0_DP, 0.0_DP ]

    dFdU(5,1) = V1(1) * ( - HJ(1) - dPdTau(1) * Tau(1)**2 &
                          - Tau(1) * dPdE(1) * ( Em(1) - 0.5_DP * VSq(1) ) )
    dFdU(5,2) = HJ(1) - dPdE(1) * Tau(1) * V1(1)**2
    dFdU(5,3) =      - dPdE(1) * Tau(1) * V1(1) * V2(1)
    dFdU(5,4) =      - dPdE(1) * Tau(1) * V1(1) * V3(1)
    dFdU(5,5) = V1(1) * ( 1.0_DP + dPdE(1) * Tau(1) )
    dFdU(5,6) = V1(1) * dPdN(1)

    dFdU(6,:) = [ - V1(1) * Y, & !Y(1) / AtomicMassUnit
                            Y, &
                  0.0_DP, 0.0_DP, 0.0_DP, V1(1) ]

    SELECT CASE( iDim )

      CASE( 1 )

        R(:,1) = [ One, V1 - Cs, V2, V3, H - Cs * V1, Y]! Ne Tau -> Y
        R(:,2) = [ Zero, Zero, One, Zero, V2, Zero ]
        R(:,3) = [ One, V1, Zero, Zero, B, Zero ]
        R(:,4) = [ One, V1, Zero, Zero, Zero, (Tau * X) / (2.0_DP * dPdN) ]
        R(:,5) = [ Zero, Zero, Zero, One, V3, Zero ]
        R(:,6) = [ One, V1 + Cs, V2, V3, H + Cs * V1, Y] ! Ne Tau -> Y

        invR(:,1) = invCsSq(1) * &
            [ + (2.0_DP * Cs * V1 + W) * 0.25_DP, &
              - Half * V2 * W, &
              + (2.0_DP * Cs**2 * X + Alpha * W / Tau)/(2.0_DP * X), &
              - (Y) * dPdN * W / (X * Tau), &
              - Half * V3 * W, &
              + 0.25_DP * (W - 2.0_DP * Cs * V1) ]

        invR(:,2) = invCsSq(1) * &
            [ - Half * (Cs + Phi(1)), &
              + Phi(1) * V2, &
              - Phi(1) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdN * Phi(1) / (X * Tau), &
              + Phi(1) * V3, &
              + Half * (Cs - Phi(1)) ]

        invR(:,3) = invCsSq(1) * &
            [ - Half * Phi(2), &
              + Cs**2 + Phi(2) * V2, &
              - Phi(2) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdN * Phi(2) / (X * Tau), &
              + Phi(2) * V3, &
              - Half * Phi(2) ]

        invR(:,4) = invCsSq(1) * &
            [ - Half * Phi(3), &
              + Phi(3) * V2, &
              - Phi(3) * Alpha / (X * Tau), &
              + 2.0_DP * Y * dPdN * Phi(3) / (X * Tau), &
              + Cs**2 + Phi(3) * V3, &
              - Half * Phi(3) ]

        invR(:,5) = invCsSq(1) * &
            [ + Half * dPdE * Tau, &
              - Phi(2), &
              + dPdE * Alpha  / X, &   ! Get rid of * Tau
              - ((2.0_DP * Y * dPdN * dPdE)/ X), &
              - Phi(3), &
              + Half * dPdE * Tau ]

        invR(:,6) = invCsSq(1) * &
            [ + Half * dPdN, &
              - V2 * dPdN, &
              + (dPdN * (-2.0_DP * Cs**2 + Alpha))/(Tau * X), &
              + 2.0_DP * dPdN * (Cs**2 - Y * dPdN)/(Tau * X), &
              - dPdN * V3, &
              + Half * dPdN ]

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X1):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

        ! Matrix = MATMUL( TRANSPOSE(invR), dFdU ) - MATMUL( DD, TRANSPOSE(invR) )
        !
        ! WRITE(*,*)
        ! WRITE(*,'(A8,A)') '', 'L^T dFdU - D L^T:'
        ! WRITE(*,*)
        ! DO i = 1, 6
        !   WRITE(*,'(A8,6ES20.10E3)') '', Matrix(i,:)
        ! END DO
        !
        ! Matrix = MATMUL( dFdU, R ) - MATMUL( R, DD )
        !
        ! WRITE(*,*)
        ! WRITE(*,'(A8,A)') '', 'dFdU R - R D:'
        ! WRITE(*,*)
        ! DO i = 1, 6
        !   WRITE(*,'(A8,6ES20.10E3)') '', Matrix(i,:)
        ! END DO

        ! Matrix =  dFdU  - MATMUL( MATMUL( R, DD ), invR )
        !
        ! WRITE(*,*)
        ! WRITE(*,'(A8,A)') '', 'dFdU - R D invR:'
        ! WRITE(*,*)
        ! DO i = 1, 6
        !   WRITE(*,'(A8,6ES20.10E3)') '', Matrix(i,:)
        ! END DO


    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition


END MODULE Euler_CharacteristicDecompositionModule_TABLE
