MODULE Euler_CharacteristicDecompositionModule

  USE KindModule,           ONLY: &
    DP
  USE FluidFieldsModule,    ONLY: &
    nCF
  USE GeometryFieldsModule, ONLY: &
    nGF







  USE Euler_CharacteristicDecompositionModule_Relativistic



  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_ComputeCharacteristicDecomposition


CONTAINS


  SUBROUTINE Euler_ComputeCharacteristicDecomposition( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)








    CALL Euler_ComputeCharacteristicDecomposition_Relativistic &
           ( iDim, G, U, R, invR )



  END SUBROUTINE Euler_ComputeCharacteristicDecomposition


END MODULE Euler_CharacteristicDecompositionModule

