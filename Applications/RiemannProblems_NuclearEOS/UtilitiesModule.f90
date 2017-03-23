MODULE UtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE FluidFieldsModule, ONLY: &
    nCF!, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
!    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne        
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

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeEigenvectors_R
  PUBLIC :: ComputeEigenvectors_L

CONTAINS

  SUBROUTINE ComputeEigenvectors_R(D, T, Ye, V1, V2, V3, VR, Componentwise)
    
    REAL(DP), DIMENSION(1),       INTENT(in)  :: D, T, Ye, V1, V2, V3
    REAL(DP), DIMENSION(6,6),     INTENT(out) :: VR
    REAL(DP), DIMENSION(1)                    :: dPdE, dPdN, dPdTau, dEdY, dEdD, dEdT, dPdY, dPdT, dPdD
    REAL(DP), DIMENSION(6,6)                  :: A
    ! What should the dimensions be???

    LOGICAL,                      INTENT(in)  :: Componentwise

    REAL(DP), DIMENSION(1)                    :: k, h, B, E, P, Cs

    INTEGER                                   :: INFO, LWORK

    IF( Componentwise )THEN

      R1(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,5) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      R1(:,6) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]

    ELSE

      CALL ComputePressure_TABLE( [D], [T], [Y], [P], dPdD, dPdT, dPdY )
      CALL ComputeSpecificInternalEnergy_TABLE( [D], [T], [Y], [E], dEdD, dEdT, dEdY )    

! compute Cs

      Tau(1) = 1/D(1)

      dPdE(1) = (1/dEdT(1)) * (dPdT(1)) 
      dPdN(1) = ( BaryonMass * Tau(1) ) * ( dPdY(1) - dEdY(1) * dPdE(1) ) 
      dPdTau(1) = (-Tau(1)**(-2)) * (dPdD(1) - (Y(1)/BaryonMass)*(dPdN(1)) - dEdD(1) * dPdE(1) ) 
      N(1) = ( ( D(1) / BaryonMass ) * Y(1) ) 

      k(1) = ( - N(1) * dPdN(1) + ( dPdTau(1) + dPdE(1) * (E(1) + (V1(1)**2 + V2(1)**2 + V3(1)**2) / 2 ) ) )
      h(1) = ( Cs(1)**2 ) / (dPdE(1) * Tau(1)) + k(1)
      B(1) = dPdTau(1) * Tau(1) - dPdE(1) * ( 0.5 * (-V1(1)**2 + V2(1)**2 + V3(1)**2 ) )

      A(:,1) = [ 0, 1, 0, 0, 0, 0]
      A(2,1) = -V(1)**2 - (Tau(1)**2)*dPdTau(1) - Tau(1) * dPdE(1) &
               * (E(1) - 0.5*(V(1)**2 + V(2)**2 + V(3)**2 ) ) 
      A(2,2) = V(1)*( 2 - Tau(1) * dPdE(1) )
      A(2,3) = - dPdE(1) * V(2) * Tau(1)
      A(2,4) = - dPdE(1) * V(3) * Tau(1)
      A(2,5) = dPdE(1) * Tau(1)
      A(2,6) = dPdN(1)
      A(:,3) = [ -V(1)*V(2), V(2), V(1), 0, 0, 0]
      A(:,4) = [ -V(1)*V(3), V(3), 0, V(1), 0, 0]
      A(5,1) = v(1) * ( -H(1) - dPdTau(1) * (Tau(1)**2) - Tau(1) * dPdE(1) &
               * ( E(1)  - 0.5 * (v(1)**2 + v(2)**2 + v(3)**2)) )
      A(5,2) = H(1) - dPdE(1) * V(1)**2 * Tau(1)
      A(5,3) = - dPdE(1) * V(1) * V(2) * Tau(1)
      A(5,4) = - dPdE(1) * V(1) * V(3) * Tau(1)
      A(5,5) = V(1) * ( 1 + dPdE(1) * Tau(1) )
      A(5,6) = V(1) * dPdN(1)
      A(:,6) = [ - ( V(1) / BaryonMass ) * Y(1), Y(1) / BaryonMass, 0, 0, 0, V(1) ]

               
      LWORK = -1

      CALL DGEEV('N', 'V', 6, A, 6, WR, WI, 6, VR, 6, TEMP, LWORK, INFO)

      LWORK = TEMP(1)
      ALLOCATE(WORK(LWORK))

      CALL DGEEV('N', 'V', 6, A, 6, WR, WI, 6, VR, 6, WORK, LWORK, INFO)  


  END SUBROUTINE ComputeEigenvectors_R


  SUBROUTINE ComputeEigenvectors_L(D, T, Ye, V1, V2, V3, R1, Componentwise)

    REAL(DP),                     INTENT(in)  :: D, T, Ye, V1, V2, V3
    REAL(DP), DIMENSION(nCF,nCF), INTENT(out) :: VR, VL
    REAL(DP), DIMENSION(1)                    :: dPdE, dPdN, dPdTau, dEdY, dEdD, dEdT, dPdY, dPdT, dPdD
    REAL(DP), DIMENSION(6,6)                  :: A, VL, VR, DD
    ! What should the dimensions be???

    LOGICAL,                      INTENT(in)  :: Componentwise

    REAL(DP), DIMENSION(1)                    :: k, h, B, E, P, Cs

    INTEGER                                   :: INFO, LWORK

    IF( Componentwise )THEN

      R1(:,1) = [ 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,2) = [ 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,3) = [ 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,4) = [ 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 0.0_DP ]
      R1(:,5) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP ]
      R1(:,6) = [ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP ]
 
    ELSE

      CALL ComputePressure_TABLE( [D], [T], [Y], [P], dPdD, dPdT, dPdY )
      CALL ComputeSpecificInternalEnergy_TABLE( [D], [T], [Y], [E], dEdD, dEdT, dEdY )

! compute Cs

      Tau(1) = 1/D(1)

      dPdE(1) = (1/dEdT(1)) * (dPdT(1))
      dPdN(1) = ( BaryonMass * Tau(1) ) * ( dPdY(1) - dEdY(1) * dPdE(1) )
      dPdTau(1) = (-Tau(1)**(-2)) * (dPdD(1) - (Y(1)/BaryonMass)*(dPdN(1)) - dEdD(1) * dPdE(1) )
      N(1) = ( ( D(1) / BaryonMass ) * Y(1) )

      k(1) = ( - N(1) * dPdN(1) + ( dPdTau(1) + dPdE(1) * (E(1) + (V1(1)**2 + V2(1)**2 + V3(1)**2) / 2 ) ) )
      h(1) = ( Cs(1)**2 ) / (dPdE(1) * Tau(1)) + k(1)
      B(1) = dPdTau(1) * Tau(1) - dPdE(1) * ( 0.5 * (-V1(1)**2 + V2(1)**2 + V3(1)**2 ) )

      A(:,1) = [ 0, 1, 0, 0, 0, 0]
      A(2,1) = -V(1)**2 - (Tau(1)**2)*dPdTau(1) - Tau(1) * dPdE(1) &
               * (E(1) - 0.5*(V(1)**2 + V(2)**2 + V(3)**2 ) )
      A(2,2) = V(1)*( 2 - Tau(1) * dPdE(1) )
      A(2,3) = - dPdE(1) * V(2) * Tau(1)
      A(2,4) = - dPdE(1) * V(3) * Tau(1)
      A(2,5) = dPdE(1) * Tau(1)
      A(2,6) = dPdN(1)
      A(:,3) = [ -V(1)*V(2), V(2), V(1), 0, 0, 0]
      A(:,4) = [ -V(1)*V(3), V(3), 0, V(1), 0, 0]
      A(5,1) = v(1) * ( -H(1) - dPdTau(1) * (Tau(1)**2) - Tau(1) * dPdE(1) &
               * ( E(1)  - 0.5 * (v(1)**2 + v(2)**2 + v(3)**2)) )
      A(5,2) = H(1) - dPdE(1) * V(1)**2 * Tau(1)
      A(5,3) = - dPdE(1) * V(1) * V(2) * Tau(1)
      A(5,4) = - dPdE(1) * V(1) * V(3) * Tau(1)
      A(5,5) = V(1) * ( 1 + dPdE(1) * Tau(1) )
      A(5,6) = V(1) * dPdN(1)
      A(:,6) = [ - ( V(1) / BaryonMass ) * Y(1), Y(1) / BaryonMass, 0, 0, 0, V(1) ]


      LWORK = -1

      CALL DGEEV('V', 'N', 6, A, 6, WR, WI, VL, 6, 6, TEMP, LWORK, INFO)

      LWORK = TEMP(1)
      ALLOCATE(WORK(LWORK))

      CALL DGEEV('V', 'N', 6, A, 6, WR, WI, VL, 6, 6, WORK, LWORK, INFO)


  END SUBROUTINE ComputeEigenvectors_L 

END MODULE UtilitiesModule      
