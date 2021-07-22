MODULE TaggingModule

  USE ISO_C_BINDING

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TagError_MF_uCF

CONTAINS


  SUBROUTINE TagError_MF_uCF &
    ( iX_B0, iX_E0, iLo, iHi, uCF, uCFerr, &
      SetTag, ClearTag, Tag, TagLo, TagHi )

    INTEGER,  INTENT(in) :: iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(3), TagHi(3)
    REAL(DP), INTENT(in) :: uCF(iLo(1):iHi(1),iLo(2):iHi(2), &
                                iLo(3):iHi(3),iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: uCFerr
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3))

    REAL(DP) :: U(1:nDOFX,iX_B0(1):iX_E0(1), &
                          iX_B0(2):iX_E0(2), &
                          iX_B0(3):iX_E0(3), &
                  1:nCF)

    INTEGER :: iX1, iX2, iX3

    REAL(DP) :: uCFerr_this

    uCFerr_this = uCFerr

    CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo, iX_B0, iX_E0, uCF, U )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( ANY( U(:,iX1,iX2,iX3,iCF_D) .GE. uCFerr_this ) )THEN

        Tag(iX1,iX2,iX3) = SetTag

      ELSE

        Tag(iX1,iX2,iX3) = ClearTag

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagError_MF_uCF

END MODULE TaggingModule
