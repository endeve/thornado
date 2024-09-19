MODULE TaggingModule

  USE ISO_C_BINDING

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    iCF_E

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Two
  USE InputParsingModule, ONLY: &
    t_new

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TagElements

CONTAINS


  SUBROUTINE TagElements &
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCF, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCF(iLo(1):iHi(1),iLo(2):iHi(2), &
                                iLo(3):iHi(3),iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))

    INTEGER :: iX1, iX2, iX3, indLo, indHi

    REAL(DP) :: TagCriteria_this, Gradient(nDOFX)

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      indLo = 1 + nDOFX * ( iCF_E - 1 )
      indHi = nDOFX * iCF_E

      Gradient = ( uCF(iX1+1,iX2,iX3,indLo:indHi) &
                 - uCF(iX1-1,iX2,iX3,indLo:indHi) )**2 &
                / ( Two * MeshX(1) % Width(iX1) )**2

      Gradient = SQRT( Gradient )

      IF( ANY( Gradient .GT. TagCriteria_this ) )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

      IF( t_new(0) .LT. 0.1_DP .AND. MeshX(1) % Center(iX1) .LT. 0.1_DP ) &
        Tag(iX1,iX2,iX3,1) = SetTag

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements


END MODULE TaggingModule
