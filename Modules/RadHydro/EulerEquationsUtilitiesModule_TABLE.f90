MODULE EulerEquationsUtilitiesModule_TABLE

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    AtomicMassUnit
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputePressure_TABLE, &
    ComputeSpecificInternalEnergy_TABLE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeEigenvectors_L
  PUBLIC :: ComputeEigenvectors_R

CONTAINS


  SUBROUTINE ComputeEigenvectors_L &
               ( D, T, Y, V1, V2, V3, lambda, VL, A0, Componentwise )

    REAL(DP), DIMENSION(1),   INTENT(in)  :: D, T, Y, V1, V2, V3
    REAL(DP), DIMENSION(6),   INTENT(out) :: lambda
    REAL(DP), DIMENSION(6,6), INTENT(out) :: VL, A0
    LOGICAL,                  INTENT(in)  :: Componentwise

    REAL(DP), DIMENSION(1) :: dPdD, dPdT, dPdY
    REAL(DP), DIMENSION(1) :: dEdD, dEdT, dEdY
    REAL(DP), DIMENSION(1) :: dPdE, dPdN, dPdTau
    REAL(DP), DIMENSION(1) :: Tau, P, E, N, H, VSq, w, X, Alpha, Delta
    REAL(DP), DIMENSION(3) :: Phi
    ! E is specific internal energy: epsilon
    ! Y is Ye

    INTEGER                             :: INFO, LWORK, i
    REAL(DP), DIMENSION(1)              :: TEMP, Cs
    REAL(DP), DIMENSION(6)              :: WR, WI
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WORK
    REAL(DP), DIMENSION(6,6)            :: A

    IF( Componentwise )THEN

      A(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      A(:,5) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      A(:,6) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

      A0 = A
      VL = A

    ELSE

      Tau = 1.0_DP / D

      CALL ComputePressure_TABLE &
             ( D, T, Y, P, dPdD, dPdT, dPdY )

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( D, T, Y, E, dEdD, dEdT, dEdY )

      dPdE   = dPdT / dEdT
      dPdN   = ( AtomicMassUnit * Tau ) * ( dPdY - dEdY * dPdE )
      dPdTau = - ( dPdD - Y * dPdN / AtomicMassUnit - dEdD * dPdE ) / Tau**2

      Cs = SQRT(Tau * (dPdN * N + P * dPdE * Tau - dPdTau * Tau ))

      Phi(1) = dPdE(1) * Tau(1) * V1(1)
      Phi(2) = dPdE(1) * Tau(1) * V2(1)
      Phi(3) = dPdE(1) * Tau(1) * V3(1)

      w = Tau(1) * (dPdE(1)*(Delta + 2.0_DP * E(1)) &
                 - 2.0_DP * dPdTau(1) * Tau(1))

      Delta = V1(1)**2 - V2(1)**2 - V3(1)**2

      X = (dPdE(1) * ( Delta + 2.0_DP * E) + 2.0_DP * dPdTau(1) * Tau(1))

      Alpha = 2.0_DP * N * dPdN(1) - X

      lambda(1) = V1(1) - Cs(1)
      lambda(2) = V1(1)
      lambda(3) = V1(1)
      lambda(4) = V1(1)
      lambda(5) = V1(1)
      lambda(6) = V1(1) + Cs(1)


      VL(1,:) = (1.0_DP / Cs(1)**2) * &
                [ (2.0_DP * Cs * V1(1) + w)/ 4.0_DP, -(Cs + Phi(1))/2.0_DP, &
                - 0.5_DP * Phi(2), - 0.5_DP * Phi(3), 0.5_DP * dPdE(1) * Tau(1),  &
                + 0.5_DP * dPdN(1) ]


      VL(2,:) = (1.0_DP / Cs(1)**2) * &
                [ - 0.5_DP * V2(1) * w, Phi(1) * V2(1), Cs**2 + Phi(2) * V2(1), &
                  + Phi(2) * V3(1), - Phi(2), - dPdN(1) * V2(1) ]

      VL(3,:) = (1.0_DP / Cs(1)**2) * &
                [ Cs**2 + (Alpha * W +  2.0_DP * Cs**2 * X)/(2.0_DP * X), &
                - (Phi(1) * Alpha)/X, - (Phi(2) * Alpha)/X, - (Phi(3) * Alpha)/X, &
                (dPdE(1) * Tau(1) * Alpha)/X, dPdN(1) * (-2.0_DP * Cs**2 &
                + Tau(1) * Alpha)/(Tau(1) * X)]

      VL(4,:) = (1.0_DP / Cs(1)**2) * &
                [ - (N * dPdN(1) * w)/X, (2.0_DP * N * dPdN(1) * Phi(1))/X, &
                  + (2.0_DP * N * dPdN(1) * Phi(2))/X, &
                  + (2.0_DP * N * dPdN(1) * Phi(3))/X, &
                  - (2.0_DP * N * dPdN(1) * dPdE(1) * Tau(1) )/X, &
                  + (2.0_DP * dPdN(1) * (Cs**2 * N * - N * dPdN(1) * Tau(1)) )/( &
                     Tau(1) * X) ]

      VL(5,:) = (1.0_DP / Cs(1)**2) * &
                [ - 0.5_DP * V3(1) * w, Phi(3) * V1(1), Phi(3) * V2(1), &
                  + Cs**2 + Phi(3) * V3(1), - Phi(3), - dPdN(1) * V3(1) ]

      VL(6,:) = (1.0_DP / Cs(1)**2) * &
                [ (- 2.0_DP * Cs * V1(1) + w)/(4.0_DP), 0.5_DP  * (Cs - Phi(1)), &
                  - 0.5_DP * Phi(2), - 0.5_DP * Phi(3), 0.5_DP * dPdE(1) * Tau(1), &
                  + 0.5_DP * dPdN(1) ]


      ! CALL ComputeFluxJacobian( D, T, Y, V1, V2, V3, A )
      !
      ! A0 = A
      !
      ! ! --- Workspace Query ---
      !
      ! LWORK = - 1
      ! CALL DGEEV( 'V', 'N', 6, A, 6, WR, WI, VL, 6, 0, 6, TEMP, LWORK, INFO )
      !
      ! ! --- Compute Eigenvectors ---
      !
      ! LWORK = TEMP(1)
      ! ALLOCATE( WORK(LWORK) )
      ! CALL DGEEV( 'V', 'N', 6, A, 6, WR, WI, VL, 6, 0, 6, WORK, LWORK, INFO )
      !
      ! DO i = 1, 6
      !   lambda(i) = WR(i)
      ! END DO
      !
      ! VL = TRANSPOSE( VL )


    END IF

  END SUBROUTINE ComputeEigenvectors_L


  SUBROUTINE ComputeEigenvectors_R &
               ( D, T, Y, V1, V2, V3, lambda, VR, A0, Componentwise )

    REAL(DP), DIMENSION(1),   INTENT(in)  :: D, T, Y, V1, V2, V3
    REAL(DP), DIMENSION(6),   INTENT(out) :: lambda
    REAL(DP), DIMENSION(6,6), INTENT(out) :: VR, A0
    LOGICAL,                  INTENT(in)  :: Componentwise

    INTEGER                             :: INFO, LWORK, i
    INTEGER,  DIMENSION(6)              :: IPIV
    REAL(DP), DIMENSION(1)              :: dPdD, dPdT, dPdY
    REAL(DP), DIMENSION(1)              :: dEdD, dEdT, dEdY
    REAL(DP), DIMENSION(1)              :: dPdE, dPdN, dPdTau
    REAL(DP), DIMENSION(1)              :: TEMP, h, k, Cs
    REAL(DP), DIMENSION(1)              :: Tau, P, E, N
    REAL(DP), DIMENSION(1)              :: VSq, X, Alpha, Delta
    REAL(DP), DIMENSION(6)              :: WR, WI
    REAL(DP), DIMENSION(:), ALLOCATABLE :: WORK
    REAL(DP), DIMENSION(6,6)            :: A, VL

    IF( Componentwise )THEN

      A(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      A(:,5) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      A(:,6) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

      A0 = A
      VR = A

    ELSE

      Tau = 1.0_DP / D

      CALL ComputePressure_TABLE &
             ( D, T, Y, P, dPdD, dPdT, dPdY )

      CALL ComputeSpecificInternalEnergy_TABLE &
             ( D, T, Y, E, dEdD, dEdT, dEdY )

      dPdE   = dPdT / dEdT
      dPdN   = ( AtomicMassUnit * Tau ) * ( dPdY - dEdY * dPdE )
      dPdTau = - ( dPdD - Y * dPdN / AtomicMassUnit - dEdD * dPdE ) / Tau**2

      Cs = SQRT(Tau * (dPdY * Y + P * dPdE * Tau - dPdTau * Tau ))

      VSq = V1(1)**2 + V2(1)**2 + V3(1)**2

      Delta = V1(1)**2 - V2(1)**2 - V3(1)**2

      X = (dPdE(1) * ( Delta + 2.0_DP * E) + 2.0_DP * dPdTau(1) * Tau(1))

      k = ( (Y(1) / AtomicMassUnit * Tau(1)) * dPdN(1) + dPdE(1) * ( &
            0.5_DP * Vsq + E) + dPdTau(1) * Tau(1))/(dPdE(1))

      h = Cs**2 / (dPdE(1) * Tau(1)) + k

      lambda(1) = V1(1) - Cs(1)
      lambda(2) = V1(1)
      lambda(3) = V1(1)
      lambda(4) = V1(1)
      lambda(5) = V1(1)
      lambda(6) = V1(1) + Cs(1)

      VR(1,:) = [ 1.0_DP, 0.0_DP, 1.0_DP, 1.0_DP, 0.0_DP, 1.0_DP ]
      VR(2,:) = [ V1(1) - Cs, 0.0_DP, V1(1), V1(1), 0.0_DP, V1(1) + Cs]
      VR(3,:) = [ V2(1), 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, V2(1)]
      VR(4,:) = [ V3(1), 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, V3(1)]
      VR(5,:) = [ h - Cs * V1(1), V2(1), 0.5_DP * (Delta + 2.0_DP * E + &
                   (2.0_DP * dPdTau(1) * Tau(1))/dPdE(1)), 0.0_DP, V3(1), &
                     h + Cs * V1(1) ]
      VR(6,:) = [(Y(1) / AtomicMassUnit), 0.0_DP, 0.0_DP, (Tau(1) * X)/(&
                  2.0_DP * dPdN(1)), 0.0_DP, (Y(1) / AtomicMassUnit) ]

      ! CALL ComputeEigenvectors_L &
      !        ( D, T, Y, V1, V2, V3, lambda, VL, A, Componentwise )
      !
      ! A0 = A
      !
      ! CALL DGETRF( 6, 6, VL, 6, IPIV, INFO )
      !
      ! ! --- Workspace Query ---
      !
      ! LWORK = -1
      ! CALL DGETRI( 6, VL, 6, IPIV, TEMP, LWORK, INFO )
      !
      ! ! --- Compute Inverse of VL ---
      !
      ! LWORK = TEMP(1)
      ! ALLOCATE( WORK(LWORK) )
      ! CALL DGETRI( 6, VL, 6, IPIV, WORK, LWORK, INFO )
      !
      ! VR = VL

    END IF

  END SUBROUTINE ComputeEigenvectors_R


  SUBROUTINE ComputeFluxJacobian( D, T, Y, V1, V2, V3, dFdU )

    REAL(DP), DIMENSION(1),   INTENT(in)  :: D, T, Y, V1, V2, V3
    REAL(DP), DIMENSION(6,6), INTENT(out) :: dFdU

    REAL(DP), DIMENSION(1) :: dPdD, dPdT, dPdY
    REAL(DP), DIMENSION(1) :: dEdD, dEdT, dEdY
    REAL(DP), DIMENSION(1) :: dPdE, dPdN, dPdTau
    REAL(DP), DIMENSION(1) :: Tau, P, E, N, H, VSq

    CALL ComputePressure_TABLE &
           ( D, T, Y, P, dPdD, dPdT, dPdY )

    CALL ComputeSpecificInternalEnergy_TABLE &
           ( D, T, Y, E, dEdD, dEdT, dEdY )

    Tau = 1.0_DP / D

    dPdE   = dPdT / dEdT
    dPdN   = ( AtomicMassUnit * Tau ) * ( dPdY - dEdY * dPdE )
    dPdTau = - ( dPdD - Y * dPdN / AtomicMassUnit - dEdD * dPdE ) / Tau**2

    N   = D * Y / AtomicMassUnit
    H   = E + 0.5_DP * V1**2 + V2**2 + V3**2 + P * Tau
    VSq = V1(1)**2 + V2(1)**2 + V3(1)**2

    dFdU(1,:) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU(2,1) = - V1(1)**2 - Tau(1)**2 * dPdTau(1) &
                - Tau(1) * dPdE(1) * ( E(1) -  0.5_DP * VSq(1) )
    dFdU(2,2) = V1(1) * ( 2.0_DP - Tau(1) * dPdE(1) )
    dFdU(2,3) = - dPdE(1) * Tau(1) * V2(1)
    dFdU(2,4) = - dPdE(1) * Tau(1) * V3(1)
    dFdU(2,5) = + dPdE(1) * Tau(1)
    dFdU(2,6) = + dPdN(1)

    dFdU(3,:) = [ - V1(1) * V2(1), V2(1), V1(1), 0.0_DP, 0.0_DP, 0.0_DP ]

    dFdU(4,:) = [ - V1(1) * V3(1), V3(1), 0.0_DP, V1(1), 0.0_DP, 0.0_DP ]

    dFdU(5,1) = V1(1) * ( - H(1) - dPdTau(1) * Tau(1)**2 &
                          - Tau(1) * dPdE(1) * ( E(1) - 0.5_DP * VSq(1) ) )
    dFdU(5,2) = H(1) - dPdE(1) * Tau(1) * V1(1)**2
    dFdU(5,3) =      - dPdE(1) * Tau(1) * V1(1) * V2(1)
    dFdU(5,4) =      - dPdE(1) * Tau(1) * V1(1) * V3(1)
    dFdU(5,5) = V1(1) * ( 1.0_DP + dPdE(1) * Tau(1) )
    dFdU(5,6) = V1(1) * dPdN(1)

    dFdU(6,:) = [ - V1(1) * Y(1) / AtomicMassUnit, &
                            Y(1) / AtomicMassUnit, &
                  0.0_DP, 0.0_DP, 0.0_DP, V1(1) ]

  END SUBROUTINE ComputeFluxJacobian


END MODULE EulerEquationsUtilitiesModule_TABLE
