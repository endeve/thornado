MODULE Euler_CharacteristicDecompositionModule

  USE KindModule, ONLY: &
    DP
  USE FluidFieldsModule, ONLY: &
    nCF
  USE GeometryFieldsModule, ONLY: &
    nGF

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

!!$  USE Euler_CharacteristicDecompositionModule_Relativistic_TABLE

#else

  USE Euler_CharacteristicDecompositionModule_NonRelativistic_TABLE

#endif

#else

#ifdef HYDRO_RELATIVISTIC

  USE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL

#else

  USE Euler_CharacteristicDecompositionModule_NonRelativistic_IDEAL

#endif

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeCharacteristicDecomposition_Euler


CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_Euler &
    ( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)    :: iDim
    REAL(DP), INTENT(in)    :: G(nGF)
    REAL(DP), INTENT(inout) :: U(nCF)
    REAL(DP), INTENT(out)   :: R(nCF,nCF)
    REAL(DP), INTENT(out)   :: invR(nCF,nCF)

#ifdef MICROPHYSICS_WEAKLIB

#ifdef HYDRO_RELATIVISTIC

#else

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
           ( iDim, G, U, R, invR )

#endif

#else

#ifdef HYDRO_RELATIVISTIC

    CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
           ( iDim, G, U, R, invR )

#else

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_IDEAL &
           ( iDim, G, U, R, invR )

#endif

#endif

  END SUBROUTINE ComputeCharacteristicDecomposition_Euler


END MODULE Euler_CharacteristicDecompositionModule
