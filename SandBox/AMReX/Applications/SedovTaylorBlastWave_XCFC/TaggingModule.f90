MODULE TaggingModule

  USE ISO_C_BINDING

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE FluidFieldsModule, ONLY: &
    iDF_Sh_X1, &
    iDF_Sh_X2, &
    iDF_Sh_X3, &
    iCF_E

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
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCF, uDF, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCF(iLo(1):iHi(1),iLo(2):iHi(2), &
                                iLo(3):iHi(3),iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: uDF(iLo(1):iHi(1),iLo(2):iHi(2), &
                                iLo(3):iHi(3),iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))

    INTEGER :: iX1, iX2, iX3, indLoX1, indLoX2, indLoX3, indLoE

    REAL(DP) :: TagCriteria_this

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      indLoX1 = 1 + nDOFX * ( iDF_Sh_X1 - 1 )
      indLoX2 = 1 + nDOFX * ( iDF_Sh_X2 - 1 )
      indLoX3 = 1 + nDOFX * ( iDF_Sh_X3 - 1 )

      IF( uDF(iX1,iX2,iX3,indLoX1) .GT. TagCriteria_this &
          .OR. uDF(iX1,iX2,iX3,indLoX2) .GT. TagCriteria_this &
          .OR. uDF(iX1,iX2,iX3,indLoX3) .GT. TagCriteria_this )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

      ! --- Ensure detonation region of initial conditions is refined ---
      indLoE = 1 + nDOFX * ( iCF_E - 1 )
      IF( uCF(iX1,iX2,iX3,indLoE) .GT. 0.01_DP .AND. StepNo(0) .LT. 1 ) &
        Tag(iX1,iX2,iX3,1) = SetTag

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements


END MODULE TaggingModule
