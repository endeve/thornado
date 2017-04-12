MODULE OpacityModule_IDEAL

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeOpacities_IDEAL
  PUBLIC :: FinalizeOpacities_IDEAL
  PUBLIC :: ComputeAbsorptionOpacity_IDEAL
  PUBLIC :: ComputeScatteringOpacity_ES_IDEAL

CONTAINS


  SUBROUTINE InitializeOpacities_IDEAL

  END SUBROUTINE InitializeOpacities_IDEAL


  SUBROUTINE FinalizeOpacities_IDEAL

  END SUBROUTINE FinalizeOpacities_IDEAL


  SUBROUTINE ComputeAbsorptionOpacity_IDEAL( E, D, T, Y, Chi )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Chi

    Chi = 0.0_DP

  END SUBROUTINE ComputeAbsorptionOpacity_IDEAL


  SUBROUTINE ComputeScatteringOpacity_ES_IDEAL( E, D, T, Y, Sigma )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Sigma

    Sigma = 0.0_DP

  END SUBROUTINE ComputeScatteringOpacity_ES_IDEAL


END MODULE OpacityModule_IDEAL
