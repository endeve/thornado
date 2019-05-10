MODULE Euler_CharacteristicDecompositionModule_TABLE

  USE KindModule,           ONLY: &
    DP
  USE FluidFieldsModule,    ONLY: &
    nCF
  USE GeometryFieldsModule, ONLY: &
    nGF

#ifdef HYDRO_NONRELATIVISTIC

  USE Euler_UtilitiesModule_NonRelativistic
  USE Euler_CharacteristicDecompositionModule_TABLE_NonRelativistic

#elif HYDRO_RELATIVISTIC

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Euler_ComputeCharacteristicDecomposition_TABLE


CONTAINS


  SUBROUTINE Euler_ComputeCharacteristicDecomposition_TABLE &
    ( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

#ifdef HYDRO_NONRELATIVISTIC

    CALL Euler_ComputeCharacteristicDecomposition_TABLE_NonRelativistic &
           ( iDim, G, U, R, invR )

#elif HYDRO_RELATIVISTIC

#endif

  END SUBROUTINE Euler_ComputeCharacteristicDecomposition_TABLE

END MODULE Euler_CharacteristicDecompositionModule_TABLE
