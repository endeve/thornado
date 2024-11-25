PROGRAM EigensystemTest

  USE KindModule, ONLY: &
    DP, One
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteRank3Tensor
  USE UnitsModule, ONLY: &
    Centimeter, Gram, &
    Kelvin, Erg, &
    AtomicMassUnit
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightCGS
  USE EquationOfStateModule_TABLE, ONLY: &
    Min_D, Min_T, Min_Y, &
    Max_D, Max_T, Max_Y, &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE, &
    ApplyEquationOfState_TABLE, &
    ComputePressure_TABLE, &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_Table
  USE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE, ONLY: &
    ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE
  USE GeometryFieldsModule, ONLY: &
      nGF
  USE FluidFieldsModule, ONLY: &
      nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, &
      iCF_E, iCF_Ne

  IMPLICIT NONE

  INTEGER :: &
    n_rndm, &
    iPoint, iV
  INTEGER, PARAMETER :: &
    iD = 1, iT = 2, iY = 3, &
    nPoints = 2**8, &
    nVariables = 15

  ! --- Set custom table bounds here. ---

  !REAL(DP), PARAMETER :: &
    !Min_D = 0.15d14 * Gram / Centimeter**3, &
    !Max_D = 0.35d14 * Gram / Centimeter**3, &
    !Min_T = 0.70d11 * Kelvin, &
    !Max_T = 0.85d11 * Kelvin, &
    !Min_Y = 0.15, &
    !Max_Y = 0.15

  REAL(DP), DIMENSION(nPoints) :: &
    D, T, Y, Ne, P, E,            &
    S, Me, Mp, Mn,            &
    Xp, Xn, Xa, Xh,           &
    Gm_A, Gm_T,               &
    dPdD, dPdT, dPdY, dPdDe,  &
    dEdD, dEdT, dEdY,         &
    alpha, beta, cs_T,        &
    rndm_D,      &
    rndm_T,      &
    rndm_Y

  REAL(DP), DIMENSION(nGF) :: G
  REAL(DP), DIMENSION(nCF) :: U
  REAL(DP), DIMENSION(nCF,nCF,nPoints) :: R_1_a, invR_1_a, dFdU_1_a
  REAL(DP), DIMENSION(nCF,nCF,nPoints) :: R_1_n, invR_1_n, dFdU_1_n

  CALL InitializeEquationOfState_TABLE &
         ("wl-EOS-SFHo-25-50-100.h5")

  WRITE(*,*)
  WRITE(*,'(A4,A10,I10.10)') '', 'nPoints = ', nPoints
  WRITE(*,*)

  n_rndm = nPoints

  ! --- Initialize Density Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )

  D(:) = 10**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  CALL WriteVector(nPoints, D / Gram * Centimeter**3, 'D.dat')

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  CALL WriteVector(nPoints, T / Kelvin, 'T.dat')

  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Y )

  Y(:) = minY + ( maxY - minY ) * rndm_Y

  CALL WriteVector(nPoints, Y, 'Y.dat')

  ! --- Get derivatives with respect to D, T, and Y. ---

  CALL ApplyEquationOfState_TABLE &
         ( D, T, Y, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh, Gm_T )
  CALL ComputePressure_TABLE &
         ( D, T, Y, P, dPdD, dPdT, dPdY )
  CALL ComputeSpecificInternalEnergy_TABLE &
         ( D, T, Y, E, dEdD, dEdT, dEdY )

  ! --- Get sound speed from characteristic decomposition module. ---

  G(:) = 1.0_DP

  DO iPoint = 1, nPoints

    U(iCF_D) = D(iPoint)
    U(iCF_S1) = 0.0_DP
    U(iCF_S2) = 0.0_DP
    U(iCF_S3) = 0.0_DP
    U(iCF_E)  = D(iPoint) * E(iPoint)
    U(iCF_Ne) = D(iPoint) * Y(iPoint) / AtomicMassUnit

    E(iPoint)  = D(iPoint) * E(iPoint)
    Ne(iPoint) = D(iPoint) * Y(iPoint) / AtomicMassUnit

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
           ( 1, G, U, R_1_n(:,:,iPoint), invR_1_n(:,:,iPoint), dFdU_1_n(:,:,iPoint), cs_T(iPoint), UseAnalytic_Option = .FALSE. )

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
           ( 1, G, U, R_1_a(:,:,iPoint), invR_1_a(:,:,iPoint), dFdU_1_a(:,:,iPoint), cs_T(iPoint), UseAnalytic_Option = .TRUE. )

  END DO

  alpha(1:nPoints) = dPdT / ( dEdT * D )
  dPdDe(1:nPoints) = dPdY / D - dEdY * alpha
  beta(1:nPoints)  = dPdDe * Y + dEdD * alpha * D - dPdD
  Gm_A(1:nPoints)  = ( dPdDe * Y - beta ) * D / P  + alpha

  CALL WriteVector(nPoints, P / Erg * Centimeter**3, 'P.dat')
  CALL WriteVector(nPoints, E / Erg * Gram, 'E.dat')
  CALL WriteVector(nPoints, cs_T * SpeedOfLightCGS, 'cs_T.dat')

  CALL WriteVector(nPoints, alpha, 'alpha.dat')
  CALL WriteVector(nPoints, beta *  ( Gram / Erg ), 'beta.dat')
  CALL WriteVector(nPoints, dPdDe * ( Gram / Erg ), 'dPdDe.dat')
  CALL WriteVector(nPoints, Gm_A, 'Gm_A.dat')
  CALL WriteVector(nPoints, Gm_T, 'Gm_T.dat')

  CALL WriteRank3Tensor(6,6,nPoints, R_1_n, 'R_1_n.dat')
  CALL WriteRank3Tensor(6,6,nPoints, invR_1_n, 'invR_1_n.dat')
  CALL WriteRank3Tensor(6,6,nPoints, dFdU_1_n, 'dFdU_1_n.dat')

  CALL WriteRank3Tensor(6,6,nPoints, R_1_a, 'R_1_a.dat')
  CALL WriteRank3Tensor(6,6,nPoints, invR_1_a, 'invR_1_a.dat')
  CALL WriteRank3Tensor(6,6,nPoints, dFdU_1_a, 'dFdU_1_a.dat')

END PROGRAM EigensystemTest
