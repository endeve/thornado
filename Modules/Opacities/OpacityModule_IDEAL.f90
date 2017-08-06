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


  SUBROUTINE ComputeAbsorptionOpacity_IDEAL &
               ( E, D, T, Y, X1, X2, X3, Chi )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Chi

    INTEGER :: iX, iE

    DO iX = 1, SIZE( D )
      DO iE = 1, SIZE( E )

        Chi(iE,iX) = 0.0_DP

      END DO
    END DO

  END SUBROUTINE ComputeAbsorptionOpacity_IDEAL


  SUBROUTINE ComputeScatteringOpacity_ES_IDEAL &
               ( E, D, T, Y, X1, X2, X3, Sigma )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Sigma

    INTEGER :: iX, iE

    DO iX = 1, SIZE( D )
      DO iE = 1, SIZE( E )

        Sigma(iE,iX) = 1.0_DP / X1(iX)**1.5_DP

      END DO
    END DO

  END SUBROUTINE ComputeScatteringOpacity_ES_IDEAL


END MODULE OpacityModule_IDEAL
