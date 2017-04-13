MODULE UtilitiesModule_NuclearEOS

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    ComputePressure_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    InitializeEquationOfState_TABLE
  USE UnitsModule, ONLY: &
    Gram, Centimeter, Kelvin, &
    AtomicMassUnit, Dyne, Erg, &
    Second, MeV, Meter
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightMKS
  USE UtilitiesModule, ONLY: &
    WriteMatrix

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeEigenvectors_R
  PUBLIC :: ComputeEigenvectors_L

CONTAINS

  SUBROUTINE ComputeEigenvectors_R(D, T, Y, V1, V2, V3, lambda, VR, A0, Componentwise)
    
    REAL(DP), DIMENSION(1),       INTENT(in)  :: D, T, Y, V1, V2, V3
    LOGICAL,                      INTENT(in)  :: Componentwise
    REAL(DP), DIMENSION(6),       INTENT(out) :: lambda
    REAL(DP), DIMENSION(6,6),     INTENT(out) :: VR, A0

    REAL(DP), DIMENSION(1)                    :: dPdE, dPdN, dPdTau, dEdY, dEdD, dEdT, dPdY, dPdT, dPdD
    REAL(DP), DIMENSION(6,6)                  :: A
    REAL(DP), DIMENSION(6)                    :: WR, WI
    REAL(DP), DIMENSION(1)                    :: E, P, Cs, Tau, TEMP, N
    REAL(DP), ALLOCATABLE, DIMENSION(:)       :: WORK
    INTEGER                                   :: INFO, LWORK

    IF( Componentwise )THEN

      A(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      A(:,5) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      A(:,6) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

    ELSE
      CALL ComputePressure_TABLE( D, T, Y, P, dPdD, dPdT, dPdY )
      CALL ComputeSpecificInternalEnergy_TABLE( D, T, Y, E, dEdD, dEdT, dEdY )    

      Tau(1) = 1/D(1)

      dPdE(1) = (1/dEdT(1)) * (dPdT(1)) 
      dPdN(1) = ( AtomicMassUnit * Tau(1) ) * ( dPdY(1) - dEdY(1) * dPdE(1) )
      dPdTau(1) = (-Tau(1)**(-2)) * (dPdD(1) - (Y(1)/AtomicMassUnit)*(dPdN(1)) - dEdD(1) * dPdE(1) )      
      N(1) = ( ( D(1) / AtomicMassUnit ) * Y(1) )
      H(1) =  (E(1) + 0.5*V1(1)**2 + V2(1)**2 + V3(1)**2 + P(1) *Tau(1) )

      A(1,1) = 0.0_DP
      A(1,2) = 1.0_DP
      A(1,3) = 0.0_DP
      A(1,4) = 0.0_DP
      A(1,5) = 0.0_DP
      A(1,6) = 0.0_DP
      A(2,1) = -V1(1)**2 - (Tau(1)**2)*dPdTau(1) - Tau(1) * dPdE(1) &
               * (E(1) - 0.5*(V1(1)**2 + V2(1)**2 + V3(1)**2 ) ) 
      A(2,2) = V1(1)*( 2 - Tau(1) * dPdE(1) )
      A(2,3) = - dPdE(1) * V2(1) * Tau(1)
      A(2,4) = - dPdE(1) * V3(1) * Tau(1)
      A(2,5) = dPdE(1) * Tau(1)
      A(2,6) = dPdN(1)
      A(3,1) = - V1(1) * V2(1)
      A(3,2) = V2(1)
      A(3,3) = V1(1)
      A(3,4) = 0
      A(3,5) = 0
      A(3,6) = 0
      A(4,1) = - V1(1) * V3(1)
      A(4,2) = V3(1)
      A(4,3) = 0
      A(4,4) = V1(1)
      A(4,5) = 0
      A(4,6) = 0
      A(5,1) = V1(1) * ( -H(1) - dPdTau(1) * (Tau(1)**2) - Tau(1) * dPdE(1) &
               * ( E(1)  - 0.5 * (V1(1)**2 + V2(1)**2 + V3(1)**2)) )
      A(5,2) = H(1) - dPdE(1) * V1(1)**2 * Tau(1)
      A(5,3) = - dPdE(1) * V1(1) * V2(1) * Tau(1)
      A(5,4) = - dPdE(1) * V1(1) * V3(1) * Tau(1)
      A(5,5) = V1(1) * ( 1 + dPdE(1) * Tau(1) )
      A(5,6) = V1(1) * dPdN(1)
      A(6,1) = - ( V1(1) / AtomicMassUnit ) * Y(1)
      A(6,2) = Y(1) / (AtomicMassUnit)
      A(6,3) = 0.0_DP
      A(6,4) = 0.0_DP
      A(6,5) = 0.0_dp
      A(6,6) = V1(1)

      A0 = A

      LWORK = -1

      CALL DGEEV('N', 'V', 6, A, 6, WR, WI, 0, 6, VR, 6, TEMP, LWORK, INFO)

      LWORK = TEMP(1)
      ALLOCATE(WORK(LWORK))

      CALL DGEEV('N', 'V', 6, A, 6, WR, WI, 0, 6, VR, 6, WORK, LWORK, INFO) 

    END IF


  END SUBROUTINE ComputeEigenvectors_R


  SUBROUTINE ComputeEigenvectors_L(D, T, Y, V1, V2, V3, lambda, VL, A0, Componentwise)

    REAL(DP), DIMENSION(1),   INTENT(in)      :: D, T, Y, V1, V2, V3
    REAL(DP), DIMENSION(6),   INTENT(out)     :: lambda
    REAL(DP), DIMENSION(6,6), INTENT(out)     :: VL, A0

    REAL(DP), DIMENSION(1)                    :: dPdE, dPdN, dPdTau, dEdY, dEdD, dEdT, dPdY, dPdT, dPdD
    REAL(DP), DIMENSION(1)                    :: Tau, TEMP, N
    REAL(DP), DIMENSION(6,6)                  :: A
    REAL(DP), DIMENSION(6)                    :: WR, WI

    LOGICAL,                      INTENT(in)  :: Componentwise

    REAL(DP), DIMENSION(1)                    :: k, h, B, E, P, Cs
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: WORK
    INTEGER                                   :: INFO, LWORK

    IF( Componentwise )THEN

      A(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      A(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      A(:,5) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      A(:,6) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]
 
    ELSE

      CALL ComputePressure_TABLE( D, T, Y, P, dPdD, dPdT, dPdY )
      CALL ComputeSpecificInternalEnergy_TABLE( D, T, Y, E, dEdD, dEdT, dEdY )

      Tau(1) = 1/D(1)

      dPdE(1) = (1/dEdT(1)) * (dPdT(1))
      dPdN(1) = ( AtomicMassUnit * Tau(1) ) * ( dPdY(1) - dEdY(1) * dPdE(1) )
      dPdTau(1) = (-Tau(1)**(-2)) * (dPdD(1) - (Y(1)/AtomicMassUnit)*(dPdN(1)) - dEdD(1) * dPdE(1) )
      N(1) = ( ( D(1) / AtomicMassUnit ) * Y(1) )
      H(1) =  (E(1) + 0.5*V1(1)**2 + V2(1)**2 + V3(1)**2 + P(1) *Tau(1) ) 

      A(1,1) = 0.0_DP
      A(1,2) = 1.0_DP
      A(1,3) = 0.0_DP
      A(1,4) = 0.0_DP
      A(1,5) = 0.0_DP
      A(1,6) = 0.0_DP
      A(2,1) = -V1(1)**2 - (Tau(1)**2)*dPdTau(1) - Tau(1) * dPdE(1) &
               * (E(1) - 0.5*(V1(1)**2 + V2(1)**2 + V3(1)**2 ) )
      A(2,2) = V1(1)*( 2 - Tau(1) * dPdE(1) )
      A(2,3) = - dPdE(1) * V2(1) * Tau(1)
      A(2,4) = - dPdE(1) * V3(1) * Tau(1)
      A(2,5) = dPdE(1) * Tau(1)
      A(2,6) = dPdN(1)
      A(3,1) = - V1(1) * V2(1)
      A(3,2) = V2(1)
      A(3,3) = V1(1)
      A(3,4) = 0
      A(3,5) = 0
      A(3,6) = 0
      A(4,1) = - V1(1) * V3(1)
      A(4,2) = V3(1)
      A(4,3) = 0
      A(4,4) = V1(1)
      A(4,5) = 0
      A(4,6) = 0
      A(5,1) = V1(1) * ( -H(1) - dPdTau(1) * (Tau(1)**2) - Tau(1) * dPdE(1) &
               * ( E(1)  - 0.5 * (V1(1)**2 + V2(1)**2 + V3(1)**2)) )
      A(5,2) = H(1) - dPdE(1) * V1(1)**2 * Tau(1)
      A(5,3) = - dPdE(1) * V1(1) * V2(1) * Tau(1)
      A(5,4) = - dPdE(1) * V1(1) * V3(1) * Tau(1)
      A(5,5) = V1(1) * ( 1 + dPdE(1) * Tau(1) )
      A(5,6) = V1(1) * dPdN(1)
      A(6,1) = - ( V1(1) / AtomicMassUnit ) * Y(1)
      A(6,2) = Y(1) / (AtomicMassUnit)
      A(6,3) = 0.0_DP
      A(6,4) = 0.0_DP
      A(6,5) = 0.0_dp
      A(6,6) = V1(1)

      LWORK = -1

      CALL DGEEV('V', 'N', 6, A, 6, WR, WI, VL, 6, 0, 6, TEMP, LWORK, INFO)

      LWORK = TEMP(1)
      ALLOCATE(WORK(LWORK))

      CALL DGEEV('V', 'N', 6, A, 6, WR, WI, VL, 6, 0, 6, WORK, LWORK, INFO)

    END IF

  END SUBROUTINE ComputeEigenvectors_L 

END MODULE UtilitiesModule_NuclearEOS
