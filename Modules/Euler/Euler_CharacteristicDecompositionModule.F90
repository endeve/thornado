MODULE Euler_CharacteristicDecompositionModule

  USE KindModule,           ONLY: &
    DP
  USE FluidFieldsModule,    ONLY: &
    nCF
  USE GeometryFieldsModule, ONLY: &
    nGF

#if defined HYDRO_NONRELATIVISTIC

  USE Euler_CharacteristicDecompositionModule_NonRelativistic_IDEAL

#elif defined HYDRO_RELATIVISTIC

  USE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL

#else

  USE Euler_CharacteristicDecompositionModule_NonRelativistic_IDEAL

#endif

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

#if defined HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeCharacteristicDecomposition_NonRelativistic &
           ( iDim, G, U, R, invR )

#elif defined HYDRO_RELATIVISTIC

    CALL Euler_ComputeCharacteristicDecomposition_Relativistic &
           ( iDim, G, U, R, invR )

#else

    CALL Euler_ComputeCharacteristicDecomposition_NonRelativistic &
           ( iDim, G, U, R, invR )

#endif

  END SUBROUTINE Euler_ComputeCharacteristicDecomposition


END MODULE Euler_CharacteristicDecompositionModule
