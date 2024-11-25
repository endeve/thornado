MODULE OpacityModule_IDEAL

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  REAL(DP) :: Eta_0, Chi_0, Sig_0

  PUBLIC :: InitializeOpacities_IDEAL
  PUBLIC :: FinalizeOpacities_IDEAL
  PUBLIC :: ComputeEmissivity_IDEAL
  PUBLIC :: ComputeAbsorptionOpacity_IDEAL
  PUBLIC :: ComputeScatteringOpacity_ES_IDEAL

CONTAINS


  SUBROUTINE InitializeOpacities_IDEAL &
               ( Eta_Option, Chi_Option, Sig_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Eta_Option, Chi_Option, Sig_Option

    Eta_0 = 1.0_DP
    IF( PRESENT( Eta_Option ) ) Eta_0 = Eta_Option

    Chi_0 = 1.0_DP
    IF( PRESENT( Chi_Option ) ) Chi_0 = Chi_Option

    Sig_0 = 1.0_DP
    IF( PRESENT( Sig_Option ) ) Sig_0 = Sig_Option

  END SUBROUTINE InitializeOpacities_IDEAL


  SUBROUTINE FinalizeOpacities_IDEAL

  END SUBROUTINE FinalizeOpacities_IDEAL


  SUBROUTINE ComputeEmissivity_IDEAL( E, D, T, Y, X1, X2, X3, Eta )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Eta

    INTEGER :: iX, iE

    DO iX = 1, SIZE( D )
      DO iE = 1, SIZE( E )

        Eta(iE,iX) = Eta_0

      END DO
    END DO

  END SUBROUTINE ComputeEmissivity_IDEAL


  SUBROUTINE ComputeAbsorptionOpacity_IDEAL( E, D, T, Y, X1, X2, X3, Chi )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Chi

    INTEGER :: iX, iE

    DO iX = 1, SIZE( D )
      DO iE = 1, SIZE( E )

        Chi(iE,iX) = Chi_0

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

        Sigma(iE,iX) = Sig_0

      END DO
    END DO

  END SUBROUTINE ComputeScatteringOpacity_ES_IDEAL


END MODULE OpacityModule_IDEAL
