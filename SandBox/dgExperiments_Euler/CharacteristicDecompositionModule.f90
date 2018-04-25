MODULE CharacteristicDecompositionModule

  USE KindModule, ONLY: &
    DP, One
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE EulerEquationsUtilitiesModule_Beta, ONLY: &
    ComputePrimitive
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeCharacteristicDecomposition

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition( iDim, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

    REAL(DP) :: P(nPF,1), Pressure(1)

    CALL ComputePrimitive &
           ( [U(iCF_D )], [U(iCF_S1)], [U(iCF_S2)], &
             [U(iCF_S3)], [U(iCF_E )], [U(iCF_Ne)], &
             P(iPF_D, :), P(iPF_V1,:), P(iPF_V2,:), &
             P(iPF_V3,:), P(iPF_E, :), P(iPF_Ne,:), &
             [One], [One], [One] )

    CALL ComputePressureFromPrimitive &
           ( P(iPF_D,:), P(iPF_E,:), P(iPF_Ne,:), Pressure )

    SELECT CASE( iDim )

      CASE( 1 )

      CASE( 2 )

      CASE( 3 )

    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition


END MODULE CharacteristicDecompositionModule
