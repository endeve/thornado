PROGRAM ComputeEigensystem_NuclearEOS

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
    Second, MeV

  IMPLICIT NONE

  REAL(DP) :: BaryonMass = AtomicMassUnit

  REAL(DP), DIMENSION(1)   :: D, T, Y, P, E, H, Tau, N
  REAL(DP), DIMENSION(1)   :: TEMP
  REAL(DP), DIMENSION(1)   :: dPdD
  REAL(DP), DIMENSION(1)   :: dPdT
  REAL(DP), DIMENSION(1)   :: dPdY
  REAL(DP), DIMENSION(1)   :: dEdT
  REAL(DP), DIMENSION(1)   :: dPdE, dPdN, dPdTau, dEdY, dEdD
  REAL(DP), DIMENSION(6,6) :: J, VL, VR
  REAL(DP), DIMENSION(3)   :: v
  REAL(DP), DIMENSION(6)   :: WR, WI

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: WORK

  INTEGER                  :: INFO, LWORK, i, k

  
  D = [ 1.20d11 * Gram / Centimeter**3 ]
  T = [ 7.6 * MeV ]
  Y = [ 0.15 ]
  v = [ 0, 0, 0 ]
  Tau(1) = 1/D(1)

  CALL InitializeEquationOfState_TABLE( 'EquationOfStateTable.h5' ) 
  CALL ComputePressure_TABLE( D, T, Y, P, dPdD, dPdT, dPdY )
  CALL ComputeSpecificInternalEnergy_TABLE( D, T, Y, E, dEdD, dEdT, dEdY )

  dPdE(1) = (1/dEdT(1)) * (dPdT(1)) !* ( (Dyne / Centimeter**2) / (Erg / Gram) )
  dPdN(1) = ( BaryonMass * Tau(1) ) * ( dPdY(1) - dEdY(1) * dPdE(1) ) !* ( (Dyne / Centimeter**2) / Centimeter**3 )
  dPdTau(1) = (-Tau(1)**(-2)) * (dPdD(1) - (Y(1)/BaryonMass)*(dPdN(1)) - dEdD(1) * dPdE(1) ) !& 
          !* ( (Dyne / Centimeter**2) / (Centimeter**3 / Gram) )
  N(1) = ( ( D(1) / BaryonMass ) * Y(1) ) !* (Centimeter**(-3)) 

  H(1) = ( E(1) + 0.5*(v(1)**2 + v(2)**2 + v(3)**2) + P(1) * Tau(1) )

  print*,"D = ", D(1) / ( Gram / Centimeter**3 )
  print*,"P = ", P(1) / (Dyne / Centimeter**2 )
  print*,"tau = ",(1 / D(1)) / ( Centimeter**3 / Gram )
  print*,"Y = ", Y(1)
  print*,"T = ", T(1) / MeV
  print*,"v = (", v(1),", ", v(2),", ", v(3),")"
  print*,"E = ", E(1) / ( Erg / Gram )
  print*,"n = ", N(1) * Centimeter**3 
  print*,"dPdD = ", dPdD
  print*,"dPdT = ", dPdT
  print*,"dPdY = ", dPdY
  print*,"dEdT = ", dEdT
  print*, " "
  print*,"dPdE = ", (1/dEdT(1)) * (dPdT(1))
  print*,"dPdN = ", dPdN
  print*,"dPdTau =", - dPdD(1) * (1/D(1))**(-2) 
  print*,"H = ", H(1)
  print*," "
 ! Note to Self: E is epsilon.

  J(1,1) = 0
  J(1,2) = 1
  J(1,3) = 0
  J(1,4) = 0
  J(1,5) = 0
  J(1,6) = 0
  J(2,1) = -v(1)**2 - (Tau(1)**2)*dPdTau(1) - Tau(1) * dPdE(1) & 
          * (E(1) - 0.5*(v(1)**2 + v(2)**2 + v(3)**2 ) )
  J(2,2) = v(1)*( 2 - Tau(1) * dPdE(1) )   
  J(2,3) = - dPdE(1) * v(2) * Tau(1)
  J(2,4) = - dPdE(1) * v(3) * Tau(1)
  J(2,5) = dPdE(1) * Tau(1)
  J(2,6) = dPdN(1)
  J(3,1) = - v(1) * v(2)
  J(3,2) = v(2)
  J(3,3) = v(1)
  J(3,4) = 0
  J(3,5) = 0
  J(3,6) = 0
  J(4,1) = - v(1) * v(3)
  J(4,2) = v(3)
  J(4,3) = 0
  J(4,4) = v(1)
  J(4,5) = 0
  J(4,6) = 0
  J(5,1) = v(1) * ( -H(1) - dPdTau(1) * (Tau(1)**2) - Tau(1) * dPdE(1) & 
          * ( E(1)  - 0.5 * (v(1)**2 + v(2)**2 + v(3)**2)) )
  J(5,2) = H(1) - dPdE(1) * v(1)**2 * Tau(1)
  J(5,3) = - dPdE(1) * v(1) * v(2) * Tau(1)
  J(5,4) = - dPdE(1) * v(1) * v(3) * Tau(1)
  J(5,5) = v(1) * ( 1 + dPdE(1) * Tau(1) )
  J(5,6) = v(1) * dPdN(1)
  J(6,1) = - ( v(1) / BaryonMass ) * Y(1)
  J(6,2) = Y(1) / BaryonMass
  J(6,3) = 0
  J(6,4) = 0
  J(6,5) = 0
  J(6,6) = v(1)

!  DO i = 2, 6
!    DO k = 1, 6
!     print*, "J(",i,",",k,") = ", J(i,k)
!    END DO
!  END DO
!
!  print*," "

  LWORK = -1

  CALL DGEEV('V', 'V', 6, J, 6, WR, WI, VL, 6, VR, 6, TEMP, LWORK, INFO)

  LWORK = TEMP(1)
  ALLOCATE(WORK(LWORK))

  CALL DGEEV('V', 'V', 6, J, 6, WR, WI, VL, 6, VR, 6, WORK, LWORK, INFO)

  DO i = 1, 6
    print*,"Eigenvalue",i,"is: ", WR(i) / (Centimeter / Second )
  END DO

  print*," "
!  print*,"Analytical Sound Speed = ", sqrt( Tau(1) * ( N(1) * dPdN(1) * ( Dyne / Centimeter**2) & 
!          + P(1)*dPdE(1)*Tau(1)*( Dyne / Centimeter**2) - dPdTau(1)*Tau(1)*(Dyne/Centimeter**2)) )  / (Centimeter / Second )
  print*,"Analytical Sound Speed = ", sqrt( Tau(1) * ( N(1) * dPdN(1) + P(1) * dPdE(1) * Tau(1) - dPdTau(1)*Tau(1)) ) & 
          / ( Centimeter / Second )

END PROGRAM ComputeEigensystem_NuclearEOS
