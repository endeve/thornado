MODULE TaggingModule

  USE ISO_C_BINDING

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE MeshModule, ONLY: &
    MeshX
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

  PUBLIC :: TagElements_uCF


CONTAINS


  SUBROUTINE TagElements_uCF &
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

    REAL(DP) :: U(1:nDOFX,iX_B0(1):iX_E0(1), &
                          iX_B0(2):iX_E0(2), &
                          iX_B0(3):iX_E0(3), &
                  1:nCF)

    INTEGER :: iX1, iX2, iX3

    REAL(DP) :: TagCriteria_this

    REAL(DP) :: Radius

    TagCriteria_this = TagCriteria

    CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo, iX_B0, iX_E0, uCF, U )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

!      IF( ANY( U(:,iX1,iX2,iX3,iCF_D) .GE. TagCriteria_this ) )THEN

!      Radius = SQRT( ( MeshX(1) % Center(iX1) - 0.5_DP )**2 &
!                   + ( MeshX(2) % Center(iX2) - 0.5_DP )**2 )

!      IF( Radius .LT. TagCriteria_this )THEN

!      IF( ( MeshX(1) % Center(iX1) .GT. TagCriteria_this ) &
!           .AND. ( MeshX(2) % Center(iX2) .GT. TagCriteria_this ) )THEN

!      IF( ABS( MeshX(1) % Center(iX1) - 0.5_DP ) .LT. TagCriteria_this )THEN

!      IF( ABS( MeshX(2) % Center(iX2) ) .GT. TagCriteria_this )THEN

      IF( MeshX(1) % Center(iX1) .GT. TagCriteria_this )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_uCF

END MODULE TaggingModule
