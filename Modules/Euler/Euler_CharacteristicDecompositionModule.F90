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
    ( iDim, G, U, R, invR, iErr_Option )

    INTEGER,  INTENT(in)            :: iDim
    REAL(DP), INTENT(in)            :: G(nGF)
    REAL(DP), INTENT(in)            :: U(nCF)
    REAL(DP), INTENT(out)           :: R(nCF,nCF)
    REAL(DP), INTENT(out)           :: invR(nCF,nCF)
    INTEGER,  INTENT(out), OPTIONAL :: iErr_Option

    INTEGER :: iErr

    iErr = 0

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_TABLE &
           ( iDim, G, U, R, invR )

#elif defined HYDRO_NONRELATIVISTIC

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_IDEAL &
           ( iDim, G, U, R, invR )

#elif defined HYDRO_RELATIVISTIC

    CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
           ( iDim, G, U, R, invR, iErr )

#else

    CALL ComputeCharacteristicDecomposition_Euler_NonRelativistic_IDEAL &
           ( iDim, G, U, R, invR )

#endif

    IF( PRESENT( iErr_Option ) ) &
      iErr_Option = iErr

  END SUBROUTINE ComputeCharacteristicDecomposition_Euler


END MODULE Euler_CharacteristicDecompositionModule
