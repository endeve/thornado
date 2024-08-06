MODULE MHD_CharacteristicDecompositionModule

  USE KindModule, ONLY: &
    DP
  USE MagnetofluidFieldsModule, ONLY: &
    nCM
  USE GeometryFieldsModule, ONLY: &
    nGF

#ifdef HYDRO_RELATIVISTIC

  USE MHD_CharacteristicDecompositionModule_Relativistic_IDEAL

#endif


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeCharacteristicDecomposition_MHD


CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_MHD &
    ( iDim, G, U, EvolveOnlyMagnetic, R, invR )

    INTEGER,  INTENT(in)    :: iDim
    REAL(DP), INTENT(in)    :: G(nGF)
    REAL(DP), INTENT(inout) :: U(nCM)
    LOGICAL,  INTENT(in)    :: EvolveOnlyMagnetic
    REAL(DP), INTENT(out)   :: R(nCM,nCM)
    REAL(DP), INTENT(out)   :: invR(nCM,nCM)

#ifdef HYDRO_RELATIVISTIC

    CALL ComputeCharacteristicDecomposition_MHD_Relativistic_IDEAL &
           ( iDim, G, U, EvolveOnlyMagnetic, R, invR )

#endif

  END SUBROUTINE ComputeCharacteristicDecomposition_MHD


END MODULE MHD_CharacteristicDecompositionModule
