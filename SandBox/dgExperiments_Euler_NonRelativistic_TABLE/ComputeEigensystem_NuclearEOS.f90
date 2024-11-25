PROGRAM ComputeEigensystem_NuclearEOS

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, Centimeter, Kelvin, &
    AtomicMassUnit, Dyne, Erg, &
    Second, MeV, Meter
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightMKS
  USE UtilitiesModule, ONLY: &
    WriteMatrix
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    ApplyEquationOfState_TABLE
  USE EulerEquationsUtilitiesModule_TABLE, ONLY: &
    ComputeEigenvectors_L, &
    ComputeEigenvectors_R

  IMPLICIT NONE

  REAL(DP) :: BaryonMass = AtomicMassUnit

  REAL(DP), DIMENSION(1)   :: D, T, Y, P, S, E, Me, Mp, Mn
  REAL(DP), DIMENSION(1)   :: Xp, Xn, Xa, Xh, Gm, Tau, N
  REAL(DP), DIMENSION(6,6) :: dFdU, L, R, Matrix, DD, RL, LR
  REAL(DP), DIMENSION(3)   :: V
  REAL(DP), DIMENSION(6)   :: lambda

  INTEGER :: i, k

  D   = [ 1.20d10 * Gram / Centimeter**3 ]
  T   = [ 3.1_DP * MeV ]
  Y   = [ 0.26_DP ]
  v   = [ 0.1_DP, 0.1_DP, 0.1_DP ]
  Tau = 1.0_DP / D

  CALL InitializeEquationOfState_TABLE &
         ( 'EquationOfStateTable.h5' )

  CALL ApplyEquationOfState_TABLE &
         ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm )

  CALL ComputeEigenvectors_L &
         ( D, T, Y, v(1), v(2), v(3), lambda, L, dFdU, .FALSE. )

  DD = 0.0d0
  DO i = 1, 6
    DD(i,i) = lambda(i)
  END DO

  WRITE(*,*)
  WRITE(*,'(A8,A)') '', 'dFdU:'
  WRITE(*,*)
  DO i = 1, 6
    WRITE(*,'(A8,6ES20.10E3)') '', DFDU(i,:)
  END DO

  Matrix = MATMUL( TRANSPOSE(L), dFdU ) - MATMUL( DD, TRANSPOSE(L) )

  WRITE(*,*)
  WRITE(*,'(A8,A)') '', 'L^T dFdU - D L^T:'
  WRITE(*,*)
  DO i = 1, 6
    WRITE(*,'(A8,6ES20.10E3)') '', Matrix(i,:)
  END DO

  WRITE(*,*)
  WRITE(*,'(A8,A)') '', 'Eigenvalues (L):'
  WRITE(*,*)
  DO i = 1, 6
    WRITE(*,'(A11,I2.2,A1,ES20.10E3)') '', i, '', lambda(i)
  END DO

  WRITE(*,*)
  WRITE(*,'(A8,A)') '', 'Sound Speed (Table):'
  WRITE(*,*)
  WRITE(*,'(A8,ES20.10E3)') '', SQRT( Gm * P / D )
  WRITE(*,*)

  CALL ComputeEigenvectors_R &
         ( D, T, Y, v(1), v(2), v(3), lambda, R, dFdU, .FALSE. )

  DD = 0.0d0
  DO i = 1, 6
    DD(i,i) = lambda(i)
  END DO

  WRITE(*,*)
  WRITE(*,'(A8,A)') '', 'Eigenvalues (R):'
  WRITE(*,*)
  DO i = 1, 6
    WRITE(*,'(A11,I2.2,A1,ES20.10E3)') '', i, '', DD(i,i)
  END DO

  Matrix = MATMUL( dFdU, R ) - MATMUL( R, DD )

  WRITE(*,*)
  WRITE(*,'(A8,A)') '', 'dFdU R - R D:'
  WRITE(*,*)
  DO i = 1, 6
    WRITE(*,'(A8,6ES20.10E3)') '', Matrix(i,:)
  END DO

  RL = MATMUL( R, L )
  LR = MATMUL( L, R )

  WRITE(*,*)
  WRITE(*,'(A8,A)') '', 'R L: '
  WRITE(*,*)
  DO i = 1, 6
    WRITE(*,'(A8,6ES20.10E3)') '', RL(i,:)
  END DO

  WRITE(*,*)
  WRITE(*,'(A8,A)') '', 'L R: '
  WRITE(*,*)
  DO i = 1, 6
    WRITE(*,'(A8,6ES20.10E3)') '', LR(i,:)
  END DO
  WRITE(*,*)

END PROGRAM ComputeEigensystem_NuclearEOS
