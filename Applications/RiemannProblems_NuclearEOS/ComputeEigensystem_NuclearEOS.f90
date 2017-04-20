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
    Second, MeV, Meter
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightMKS
  USE UtilitiesModule_NuclearEOS, ONLY: &
    ComputeEigenvectors_L, &
    ComputeEigenvectors_R
  USE UtilitiesModule, ONLY: &
    WriteMatrix

  IMPLICIT NONE

  REAL(DP) :: BaryonMass = AtomicMassUnit

  REAL(DP), DIMENSION(1)   :: D, T, Y, P, E, H, Tau, N
  REAL(DP), DIMENSION(1)   :: TEMP
  REAL(DP), DIMENSION(6,6) :: A, R, L, DD
  REAL(DP), DIMENSION(3)   :: v
  REAL(DP), DIMENSION(6)   :: lambda

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: WORK

  INTEGER                  :: INFO, LWORK, i, k

  
  D = [ 1.20d10 * Gram / Centimeter**3 ]
  T = [ 3.1 * MeV ]
  Y = [ 0.26 ]
  v = [ 0, 0, 0 ]
  Tau = 1/D

  CALL InitializeEquationOfState_TABLE( 'EquationOfStateTable.h5' ) 
  CALL ComputeEigenvectors_R(D, T, Y, v(1), v(2), v(3), lambda, R, A, .FALSE.)

  DD = 0.0d0

  DO i = 1, 6
    DD(i,i) = lambda(i)
  END DO

  print*,"AR - RD = ", matmul(A,R) - matmul(R,DD)

  CALL WriteMatrix(6, 6, A, 'A.dat')
  CALL WriteMatrix(6, 6, R, 'R.dat')

!  DO i = 1, 6
!   print*,"Eigenvalue",i,"is: ", lambda(i) / ( Meter / Second  ) / SpeedOfLightMKS
!  END DO

  CALL ComputeEigenvectors_L(D, T, Y, v(1), v(2), v(3), lambda, L, A, .FALSE.)

!  DO i = 1, 6
!    print*,"Eigenvalue",i,"is: ", lambda(i) / ( Meter / Second  ) / SpeedOfLightMKS
!  END DO


!  print*,"Analytical Sound Speed = ", sqrt( Tau(1) * ( N(1) * dPdN(1) + P(1) * dPdE(1) * Tau(1) - dPdTau(1)*Tau(1))) &
!                    / ( Meter / Second  ) / SpeedOfLightMKS


  DD = 0.0d0

  DO i = 1, 6
    DD(i,i) = lambda(i)
  END DO

  print*," "
  
  print*,"L'J - DL' = ", matmul(transpose(L),A) - matmul(DD,transpose(L))

!  DO i = 1, 6
!    DO k = 1, 6
!     print*, "A(",i,",",k,") = ", A(i,k)
!    END DO
!  END DO


END PROGRAM ComputeEigensystem_NuclearEOS
