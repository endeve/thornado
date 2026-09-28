MODULE TaggingModule

  USE ISO_C_BINDING

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE MeshModule, ONLY: &
    MeshX
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer
  USE FluidFieldsModule, ONLY: &
    iDF_Sh_X1

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    StepNo

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TagElements

CONTAINS


  SUBROUTINE TagElements &
    ( iLevel, iX_B0, iX_E0, iLoC, iHiC, uCF, iLoD, iHiD, uDF, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), &
                            iLoC(4), iHiC(4), iLoD(4), iHiD(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCF(iLoC(1):iHiC(1),iLoC(2):iHiC(2), &
                                iLoC(3):iHiC(3),iLoC(4):iHiC(4))
    REAL(DP), INTENT(in) :: uDF(iLoD(1):iHiD(1),iLoD(2):iHiD(2), &
                                iLoD(3):iHiD(3),iLoD(4):iHiD(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))

    INTEGER :: iX1, iX2, iX3, indLo

    REAL(DP) :: TagCriteria_this

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      indLo = 1 + nDOFX * ( iDF_Sh_X1 - 1 )

      IF( ANY( uDF(iX1,iX2,iX3,indLo:indLo+nDOFX-1) .GT. 0.1_DP ) &
          .OR. ( ANY( uCF(iX1,iX2,iX3,1:nDOFX) &
                    .GT. 3.0e8_DP * ( Gram / Centimeter**3 ) ) ) &
        )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements


END MODULE TaggingModule
