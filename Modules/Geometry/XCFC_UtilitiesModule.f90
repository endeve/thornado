MODULE XCFC_UtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE GeometryFieldsModule, ONLY: &
    iGF_Psi

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MultiplyWithPsi6

CONTAINS


  SUBROUTINE MultiplyWithPsi6( iX_B1, iX_E1, G, U, Power )

    INTEGER,  INTENT(in)    :: iX_B1(3), iX_E1(3), Power
    REAL(DP), INTENT(in)    :: G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iNX, iX1, iX2, iX3, iFF, nFF

    nFF = SIZE( U, DIM = 5 )

    DO iFF = 1       , nFF
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iNX = 1       , nDOFX

      U(iNX,iX1,iX2,iX3,iFF) &
        = U(iNX,iX1,iX2,iX3,iFF) * G(iNX,iX1,iX2,iX3,iGF_Psi)**( 6 * Power )

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MultiplyWithPsi6


END MODULE XCFC_UtilitiesModule
