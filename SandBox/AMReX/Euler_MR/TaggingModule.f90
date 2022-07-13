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
    DP, &
    Two

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TagElements_Advection1D
  PUBLIC :: TagElements_RiemannProblem1D
  PUBLIC :: TagElements_Advection2D
  PUBLIC :: TagElements_KelvinHelmholtz2D
  PUBLIC :: TagElements_Advection3D
  PUBLIC :: TagElements_uCF

CONTAINS


  SUBROUTINE TagElements_Advection1D &
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

    REAL(DP) :: TagCriteria_this

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      indLo = 1 + nDOFX * ( iCF_D - 1 )
      indHi = nDOFX * iCF_D

!      IF( MeshX(1) % Center(iX1) .GT. TagCriteria_this )THEN
      IF( ANY( uCF(iX1,iX2,iX3,indLo:indHi) .GT. TagCriteria_this ) )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_Advection1D


  SUBROUTINE TagElements_RiemannProblem1D &
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

    REAL(DP) :: TagCriteria_this
    REAL(DP) :: GradD(nDOFX)

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      indLo = 1 + nDOFX * ( iCF_D - 1 )
      indHi = nDOFX * iCF_D

      GradD = ( (   uCF(iX1+1,iX2,iX3,indLo:indHi) &
                  - uCF(iX1-1,iX2,iX3,indLo:indHi) ) )**2! &
!                / ( Two * MeshX(1) % Width(iX1) ) )**2

      GradD = SQRT( GradD )

      !IF( MeshX(1) % Center(iX1) .GT. TagCriteria_this )THEN
      IF( ANY( GradD .GT. TagCriteria_this ) )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_RiemannProblem1D


  SUBROUTINE TagElements_Advection2D &
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

    INTEGER :: iX1, iX2, iX3

    REAL(DP) :: Radius

    REAL(DP) :: TagCriteria_this

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Radius = SQRT( ( MeshX(1) % Center(iX1) - 0.5_DP )**2 &
                   + ( MeshX(2) % Center(iX2) - 0.5_DP )**2 )

!!$      IF( ( MeshX(1) % Center(iX1) .GT. TagCriteria_this ) &
!!$           .AND. ( MeshX(2) % Center(iX2) .GT. TagCriteria_this ) )THEN

      IF( Radius .LT. TagCriteria_this )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_Advection2D


  SUBROUTINE TagElements_KelvinHelmholtz2D &
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

    REAL(DP) :: TagCriteria_this
    REAL(DP) :: GradD(nDOFX)

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      indLo = 1 + nDOFX * ( iCF_D - 1 )
      indHi = nDOFX * iCF_D

      GradD = ( (   uCF(iX1+1,iX2,iX3,indLo:indHi) &
                  - uCF(iX1-1,iX2,iX3,indLo:indHi) ))**2! &
!                / ( Two * MeshX(1) % Width(iX1) ) )**2

      GradD = GradD &
                + ( (   uCF(iX1,iX2+1,iX3,indLo:indHi) &
                      - uCF(iX1,iX2-1,iX3,indLo:indHi) ) )**2! &
!                    / ( Two * MeshX(2) % Width(iX2) ) )**2

      GradD = SQRT( GradD )

      IF( ANY( GradD .GT. TagCriteria_this ) )THEN
      !IF( MeshX(1) % Center(iX1) .GT. TagCriteria_this )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_KelvinHelmholtz2D


  SUBROUTINE TagElements_Advection3D &
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

    INTEGER :: iX1, iX2, iX3

    REAL(DP) :: Radius

    REAL(DP) :: TagCriteria_this

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Radius = SQRT( ( MeshX(1) % Center(iX1) - 0.5_DP )**2 &
                   + ( MeshX(2) % Center(iX2) - 0.5_DP )**2 &
                   + ( MeshX(3) % Center(iX3) - 0.5_DP )**2 )

!      IF( ( MeshX(1) % Center(iX1) .GT. TagCriteria_this ) &
!           .AND. ( MeshX(2) % Center(iX2) .GT. TagCriteria_this ) &
!           .AND. ( MeshX(3) % Center(iX3) .GT. TagCriteria_this ) )THEN

      IF( Radius .LT. TagCriteria_this )THEN

        Tag(iX1,iX2,iX3,1) = SetTag

      ELSE

        Tag(iX1,iX2,iX3,1) = ClearTag

      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_Advection3D


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

    INTEGER :: iX1, iX2, iX3

    REAL(DP) :: TagCriteria_this

    REAL(DP) :: Radius

    TagCriteria_this = TagCriteria

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

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
