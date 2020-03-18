MODULE Euler_CharacteristicDecompositionModule

  USE KindModule,           ONLY: &
    DP
  USE FluidFieldsModule,    ONLY: &
    nCF
  USE GeometryFieldsModule, ONLY: &
    nGF

  USE Euler_CharacteristicDecompositionModule_NonRelativistic_IDEAL
  USE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE
  USE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeCharacteristicDecomposition_Euler


CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_Euler &
    ( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
           ( iDim, G, U, R, invR )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
           ( iDim, G, U, R, invR )

#else

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_IDEAL &
           ( iDim, G, U, R, invR )

#endif

  END SUBROUTINE ComputeCharacteristicDecomposition_Euler


END MODULE Euler_CharacteristicDecompositionModule
