MODULE TaggingModule

  USE ISO_C_BINDING

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    iE_B0, &
    iE_E0, &
    nDOFX, &
    nDOFE
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate, &
    MeshE
  USE RadiationFieldsModule, ONLY: &
    iCR_N, &
    iCR_G1, &
    iCR_G2, &
    nSpecies, &
    nCR   
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Pi, &
    TwoPi, One, Zero
  USE MF_FieldsModule_TwoMoment, ONLY: &
    CreateFields_TwoMoment_MF, &
    MF_uCR
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModuleE, ONLY: &
    WeightsE

    USE amrex_multifab_module, ONLY: &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_max
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TagElements, TagElements_Density, TagElements_ShadowCasting, TagElements_TransparentVortex_Spherical, &
            TagElements_TVSD, TagElements_Exterior, TagElements_TransparentVortex_patch, &
            TagElements_TransparentVortex_MovingWedge, TagElements_TransparentVortex_Wedge

CONTAINS


SUBROUTINE TagElements &
  ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
    SetTag, ClearTag, TagLo, TagHi, Tag )
  INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                          TagLo(4), TagHi(4)
  REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1),iLo(2):iHi(2), &
                              iLo(3):iHi(3),iLo(4):iHi(4))
  REAL(DP), INTENT(in) :: TagCriteria
  CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
  CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                               TagLo(2):TagHi(2), &
                                               TagLo(3):TagHi(3), &
                                               TagLo(4):TagHi(4))
  INTEGER :: iX1, iX2, iX3
  REAL(DP) :: TagCriteria_this

  TagCriteria_this = TagCriteria

  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)
  
    IF( MeshX(1) % Center(iX1) .LT. TagCriteria_this )THEN
      Tag(iX1,iX2,iX3,1) = SetTag
    ELSE
      Tag(iX1,iX2,iX3,1) = ClearTag
    END IF
  
  END DO
  END DO
  END DO

END SUBROUTINE TagElements


SUBROUTINE TagElements_Density &
  ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
    SetTag, ClearTag, TagLo, TagHi, Tag )
  
  INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                          TagLo(4), TagHi(4)
  REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1),iLo(2):iHi(2), &
                              iLo(3):iHi(3),iLo(4):iHi(4))
  REAL(DP), INTENT(in) :: TagCriteria
  CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
  CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                               TagLo(2):TagHi(2), &
                                               TagLo(3):TagHi(3), &
                                               TagLo(4):TagHi(4))
  
  INTEGER  :: iX1, iX2, iX3
  INTEGER  :: iS, iZ1, indLo, indHi, iNodeE, iNodeX
  REAL(DP) :: TagCriteria_this
  
  TagCriteria_this = TagCriteria
  

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Tag(iX1,iX2,iX3,1) = ClearTag
    
      DO iS = 1, nSpecies
      DO iZ1 = iE_B0, iE_E0
      DO iNodeE = 1    , nDOFE
      DO iNodeX = 1    , nDOFX
        
        indLo = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFE * nDOFX &
              + ( iCR_N - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFE * nDOFX &
              + ( iZ1 - iE_B0 ) * nDOFE * nDOFX + ( iNodeE - 1 ) * nDOFX + iNodeX
        !PRINT *, indLo
        !PRINT *, iNodeX
        indHi = indLo + nDOFE * nDOFX !nDOFE * nDOFX * iCR_N !
        
        IF( ANY( ABS( uCR(iX1,iX2,iX3,indLo:indHi) ) .GT. TagCriteria_this ) ) THEN
          Tag(iX1,iX2,iX3,1) = SetTag
      !ELSE

        !Tag(iX1,iX2,iX3,1) = ClearTag
        END IF

    END DO
    END DO    
    END DO
    END DO 
    
  END DO
  END DO
  END DO
  
END SUBROUTINE TagElements_Density


SUBROUTINE TagElements_Density2 &
  ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
    SetTag, ClearTag, TagLo, TagHi, Tag )

  INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                          TagLo(4), TagHi(4)
  REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1),iLo(2):iHi(2), &
                              iLo(3):iHi(3),iLo(4):iHi(4))
  REAL(DP), INTENT(in) :: TagCriteria
  CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
  CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                               TagLo(2):TagHi(2), &
                                               TagLo(3):TagHi(3), &
                                               TagLo(4):TagHi(4))

  INTEGER  :: iX1, iX2, iX3
  INTEGER  :: iS, iZ1, indLo, indHi
  REAL(DP) :: TagCriteria_this

  TagCriteria_this = TagCriteria

  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)

    Tag(iX1,iX2,iX3,1) = ClearTag

    species_energy: DO iS  = 1, nSpecies
                    DO iZ1 = iE_B0, iE_E0

      indLo = ( iS    - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFE * nDOFX &
            + ( iCR_N - 1 ) *       ( iE_E0 - iE_B0 + 1 ) * nDOFE * nDOFX &
            + ( iZ1 - iE_B0 ) * nDOFE * nDOFX + 1
      indHi = indLo + nDOFE * nDOFX - 1

      IF( ANY( uCR(iX1,iX2,iX3,indLo:indHi) .GE. TagCriteria_this ) ) THEN
        Tag(iX1,iX2,iX3,1) = SetTag
        EXIT species_energy
      END IF

    END DO
    END DO species_energy

  END DO
  END DO
  END DO

END SUBROUTINE TagElements_Density2


SUBROUTINE TagElements_ShadowCasting2 &
  ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
    SetTag, ClearTag, TagLo, TagHi, Tag )

  INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                          TagLo(4), TagHi(4)
  REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                              iLo(3):iHi(3), iLo(4):iHi(4))
  REAL(DP), INTENT(in) :: TagCriteria
  CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
  CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                TagLo(2):TagHi(2), &
                                                TagLo(3):TagHi(3), &
                                                TagLo(4):TagHi(4))

  REAL(DP), PARAMETER :: X1_SOURCE = 3.0_DP
  REAL(DP), PARAMETER :: X2_SOURCE = 0.0_DP
  !REAL(DP), PARAMETER :: TwoPi     = 2.0_DP * 3.141592653589793_DP

  INTEGER  :: iX1, iX2, iX3
  INTEGER  :: iS, iZ1, iNodeE, iNodeX
  INTEGER  :: iNodeX1, iNodeX2
  INTEGER  :: stride, base, ind_N, ind_G1, ind_G2, ind_G3
  REAL(DP) :: TagCriteria_this
  REAL(DP) :: X1, X2, r_from_source
  REAL(DP) :: G1, G2, G3, F_mag, L_node

  TagCriteria_this = TagCriteria
  stride = ( iE_E0 - iE_B0 + 1 ) * nDOFE * nDOFX

  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)

    Tag(iX1,iX2,iX3,1) = ClearTag

    DO iS     = 1, nSpecies
    DO iZ1    = iE_B0, iE_E0
    DO iNodeE = 1, nDOFE
    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1, iNodeX)
      iNodeX2 = NodeNumberTableX(2, iNodeX)

      X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

      r_from_source = SQRT( ( X1 - X1_SOURCE )**2 + ( X2 - X2_SOURCE )**2 )

      base = ( iS    - 1 ) * nCR * stride &
           + ( iZ1   - iE_B0 ) * nDOFE * nDOFX &
           + ( iNodeE - 1 ) * nDOFX &
           + iNodeX

      ind_N  = ( iCR_N  - 1 ) * stride + base
      ind_G1 = ( iCR_G1 - 1 ) * stride + base
      ind_G2 = ( iCR_G2 - 1 ) * stride + base

      G1 = uCR(iX1, iX2, iX3, ind_G1)
      G2 = uCR(iX1, iX2, iX3, ind_G2)

      F_mag  = SQRT( G1**2 + G2**2 )
      L_node = TwoPi * r_from_source * F_mag

      IF( L_node > TagCriteria_this ) THEN
        Tag(iX1, iX2, iX3, 1) = SetTag
      END IF

    END DO
    END DO
    END DO
    END DO

  END DO
  END DO
  END DO

END SUBROUTINE TagElements_ShadowCasting2


SUBROUTINE TagElements_ShadowCasting &
  ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
    SetTag, ClearTag, TagLo, TagHi, Tag )

  INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                          TagLo(4), TagHi(4)
  REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                              iLo(3):iHi(3), iLo(4):iHi(4))
  REAL(DP), INTENT(in) :: TagCriteria
  CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
  CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                TagLo(2):TagHi(2), &
                                                TagLo(3):TagHi(3), &
                                                TagLo(4):TagHi(4))

  REAL(DP), PARAMETER :: X1_SOURCE = 3.0_DP
  REAL(DP), PARAMETER :: X2_SOURCE = 0.0_DP


  REAL(DP), PARAMETER :: L_LO_TAG = 2.0e-2_DP !5.0e-2_DP
  REAL(DP), PARAMETER :: L_HI_TAG = 8.0e-2_DP

  INTEGER  :: iX1, iX2, iX3
  INTEGER  :: iS, iZ1, iNodeE, iNodeX
  INTEGER  :: iNodeX1, iNodeX2
  INTEGER  :: stride, base, ind_G1, ind_G2, ind_G3
  REAL(DP) :: X1, X2, r_from_source
  REAL(DP) :: G1, G2, G3, F_mag, L_node
  LOGICAL  :: in_band

  stride = ( iE_E0 - iE_B0 + 1 ) * nDOFE * nDOFX

  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)

    Tag(iX1,iX2,iX3,1) = ClearTag
    in_band = .FALSE.

    DO iS     = 1, nSpecies
    DO iZ1    = iE_B0, iE_E0
    DO iNodeE = 1, nDOFE
    DO iNodeX = 1, nDOFX

      iNodeX1 = NodeNumberTableX(1, iNodeX)
      iNodeX2 = NodeNumberTableX(2, iNodeX)

      X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
      X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

      r_from_source = SQRT( ( X1 - X1_SOURCE )**2 + ( X2 - X2_SOURCE )**2 )

      base = ( iS     - 1 ) * nCR * stride &
           + ( iZ1    - iE_B0 ) * nDOFE * nDOFX &
           + ( iNodeE - 1 ) * nDOFX &
           + iNodeX

      ind_G1 = ( iCR_G1 - 1 ) * stride + base
      ind_G2 = ( iCR_G2 - 1 ) * stride + base

      G1 = uCR(iX1, iX2, iX3, ind_G1)
      G2 = uCR(iX1, iX2, iX3, ind_G2)

      F_mag  = SQRT( G1**2 + G2**2 )
      L_node = TwoPi * r_from_source * F_mag

      IF( L_node >= L_LO_TAG .AND. L_node <= L_HI_TAG ) THEN
        in_band = .TRUE.
      END IF

    END DO
    END DO
    END DO
    END DO

    IF( in_band ) Tag(iX1, iX2, iX3, 1) = SetTag

  END DO
  END DO
  END DO

END SUBROUTINE TagElements_ShadowCasting



SUBROUTINE TagElements_TransparentVortex_Spherical1 &
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                                iLo(3):iHi(3), iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))

    REAL(DP), PARAMETER :: R0_VORTEX  = 6.0_DP
    REAL(DP), PARAMETER :: TH0_VORTEX = 0.25_DP * Pi
    REAL(DP), PARAMETER :: R_TAG      = 1.7_DP

    REAL(DP) :: X1_C, X2_C, d_from_vortex
    INTEGER  :: iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      X1_C = MeshX(1) % Center(iX1)
      X2_C = MeshX(2) % Center(iX2)

      d_from_vortex &
        = SQRT( X1_C**2 + R0_VORTEX**2 &
                - 2.0_DP * X1_C * R0_VORTEX * COS( X2_C - TH0_VORTEX ) )

      IF( d_from_vortex <= R_TAG )THEN
        Tag(iX1,iX2,iX3,1) = SetTag
      ELSE
        Tag(iX1,iX2,iX3,1) = ClearTag
      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_TransparentVortex_Spherical1

  SUBROUTINE TagElements_TransparentVortex_Spherical &
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                                iLo(3):iHi(3), iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))


    REAL(DP), PARAMETER :: R_MIN = 4.5_DP
    REAL(DP), PARAMETER :: R_MAX = 7.9_DP

    REAL(DP) :: X1_C
    INTEGER  :: iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      X1_C = MeshX(1) % Center(iX1)

      IF( X1_C >= R_MIN .AND. X1_C <= R_MAX )THEN
        Tag(iX1,iX2,iX3,1) = SetTag
      ELSE
        Tag(iX1,iX2,iX3,1) = ClearTag
      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_TransparentVortex_Spherical


    SUBROUTINE TagElements_TransparentVortex_Wedge &
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                                iLo(3):iHi(3), iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))

    REAL(DP), PARAMETER :: X2_MIN = 0.125_DP * Pi   ! pi/8   (22.5 deg)
    REAL(DP), PARAMETER :: X2_MAX = 0.375_DP * Pi   ! 3pi/8  (67.5 deg)

    REAL(DP) :: X2_C
    INTEGER  :: iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      X2_C = MeshX(2) % Center(iX2)

      DO iX1 = iX_B0(1), iX_E0(1)

        IF( X2_C >= X2_MIN .AND. X2_C <= X2_MAX )THEN
          Tag(iX1,iX2,iX3,1) = SetTag
        ELSE
          Tag(iX1,iX2,iX3,1) = ClearTag
        END IF

      END DO

    END DO
    END DO

  END SUBROUTINE TagElements_TransparentVortex_Wedge


  SUBROUTINE TagElements_TransparentVortex_patch &
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                                iLo(3):iHi(3), iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))


    REAL(DP), PARAMETER :: X1_MIN = -3.0_DP
    REAL(DP), PARAMETER :: X1_MAX = +3.0_DP
    REAL(DP), PARAMETER :: X2_MIN = -3.0_DP
    REAL(DP), PARAMETER :: X2_MAX = +3.0_DP

    REAL(DP) :: X1_C, X2_C
    INTEGER  :: iX1, iX2, iX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      X1_C = MeshX(1) % Center(iX1)
      X2_C = MeshX(2) % Center(iX2)

      IF( X1_C >= X1_MIN .AND. X1_C <= X1_MAX .AND. &
          X2_C >= X2_MIN .AND. X2_C <= X2_MAX )THEN
        Tag(iX1,iX2,iX3,1) = SetTag
      ELSE
        Tag(iX1,iX2,iX3,1) = ClearTag
      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_TransparentVortex_patch

SUBROUTINE TagElements_TVSD1 &
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                                iLo(3):iHi(3), iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))

    REAL(DP), PARAMETER :: R_Domain_Lo = 2.0_DP
    REAL(DP), PARAMETER :: R_Domain_Hi = 10.0_DP


    REAL(DP), PARAMETER :: R_Exclude_Lo = 2.5_DP
    REAL(DP), PARAMETER :: R_Exclude_Hi = 0.0_DP

    REAL(DP), PARAMETER :: RMS_Target = 15.6_DP

    LOGICAL , PARAMETER :: RefineInsideBand = .FALSE.

    REAL(DP), PARAMETER :: N_Beam_Min = 1.0e-4_DP

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iS, iZ1, iNodeE, iNodeX, iNodeZ
    INTEGER  :: nDOFZ_L, stride, ind_N
    REAL(DP) :: E_0, E, W2, N_node, R_Lo, R_Hi
    REAL(DP) :: Int3, Int5, N_Peak, RMS, Dev
    LOGICAL  :: exceeds

    nDOFZ_L = nDOFE * nDOFX
    stride  = ( iE_E0 - iE_B0 + 1 ) * nDOFZ_L

    E_0 = NodeCoordinate( MeshE, iE_B0, 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Tag(iX1,iX2,iX3,1) = ClearTag

      R_Lo = MeshX(1) % Center(iX1) - 0.5_DP * MeshX(1) % Width(iX1)
      R_Hi = MeshX(1) % Center(iX1) + 0.5_DP * MeshX(1) % Width(iX1)

      IF( R_Lo .LT. R_Domain_Lo + R_Exclude_Lo ) CYCLE
      IF( R_Exclude_Hi .GT. 0.0_DP .AND. &
          R_Hi .GT. R_Domain_Hi - R_Exclude_Hi ) CYCLE

      exceeds = .FALSE.

      DO iS     = 1, nSpecies
      DO iNodeX = 1, nDOFX

        Int3   = 0.0_DP
        Int5   = 0.0_DP
        N_Peak = 0.0_DP

        DO iZ1    = iE_B0, iE_E0
        DO iNodeE = 1, nDOFE

          E  = NodeCoordinate( MeshE, iZ1, iNodeE )
          W2 = WeightsE(iNodeE) * MeshE % Width(iZ1) * E**2

          iNodeZ = ( iNodeE - 1 ) * nDOFX + iNodeX

          ind_N = ( iS    - 1     ) * nCR * stride &
                + ( iCR_N - 1     ) * stride       &
                + ( iZ1   - iE_B0 ) * nDOFZ_L      &
                +   iNodeZ

          N_node = uCR(iX1,iX2,iX3,ind_N)

          N_Peak = MAX( N_Peak, N_node )

          Int3 = Int3 + W2 * ( E / E_0 )    * N_node
          Int5 = Int5 + W2 * ( E / E_0 )**3 * N_node

        END DO
        END DO

        IF( N_Peak .LT. N_Beam_Min ) CYCLE
        IF( Int3   .LE. 0.0_DP     ) CYCLE

        RMS = E_0 * SQRT( MAX( Int5 / Int3, 0.0_DP ) )

        Dev = ABS( RMS - RMS_Target )

        IF( RefineInsideBand )THEN
          IF( Dev .LE. TagCriteria ) exceeds = .TRUE.
        ELSE
          IF( Dev .GT. TagCriteria ) exceeds = .TRUE.
        END IF

      END DO
      END DO

      IF( exceeds ) Tag(iX1,iX2,iX3,1) = SetTag

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_TVSD1

SUBROUTINE TagElements_TVSD &
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                                iLo(3):iHi(3), iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))

    REAL(DP), PARAMETER :: R_Refine_Inner = 3.0_DP
    REAL(DP), PARAMETER :: RMS_Target     = 15.0_DP
    REAL(DP), PARAMETER :: N_Beam_Min     = 1.0e-4_DP

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iS, iZ1, iNodeE, iNodeX, iNodeZ
    INTEGER  :: nDOFZ_L, stride, ind_N
    REAL(DP) :: E_0, E, W2, N_node, R_Lo
    REAL(DP) :: Int3, Int5, N_Peak, RMS
    LOGICAL  :: exceeds

    nDOFZ_L = nDOFE * nDOFX
    stride  = ( iE_E0 - iE_B0 + 1 ) * nDOFZ_L

    E_0 = NodeCoordinate( MeshE, iE_B0, 1 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      !Tag(iX1,iX2,iX3,1) = ClearTag

      R_Lo = MeshX(1) % Center(iX1)

      IF( R_Lo .LT. R_Refine_Inner )THEN
        Tag(iX1,iX2,iX3,1) = SetTag
      END IF

      exceeds = .FALSE.

      DO iS     = 1, nSpecies
      DO iNodeX = 1, nDOFX

        Int3   = 0.0_DP
        Int5   = 0.0_DP
        N_Peak = 0.0_DP

        DO iZ1    = iE_B0, iE_E0
        DO iNodeE = 1, nDOFE

          E  = NodeCoordinate( MeshE, iZ1, iNodeE )
          W2 = WeightsE(iNodeE) * MeshE % Width(iZ1) * E**2

          iNodeZ = ( iNodeE - 1 ) * nDOFX + iNodeX

          ind_N = ( iS    - 1     ) * nCR * stride &
                + ( iCR_N - 1     ) * stride       &
                + ( iZ1   - iE_B0 ) * nDOFZ_L      &
                +   iNodeZ

          N_node = uCR(iX1,iX2,iX3,ind_N)

          N_Peak = MAX( N_Peak, N_node )

          Int3 = Int3 + W2 * ( E / E_0 )    * N_node
          Int5 = Int5 + W2 * ( E / E_0 )**3 * N_node

        END DO
        END DO

        IF( N_Peak .LT. N_Beam_Min ) CYCLE
        IF( .NOT. Int3 .GT. 0.0_DP ) CYCLE

        RMS = E_0 * SQRT( MAX( Int5 / Int3, 0.0_DP ) )

        IF( ABS( RMS - RMS_Target ) .GT. TagCriteria ) exceeds = .TRUE.

      END DO
      END DO

      IF( exceeds ) Tag(iX1,iX2,iX3,1) = SetTag

    END DO
    END DO
    END DO

  END SUBROUTINE TagElements_TVSD

    SUBROUTINE TagElements_TransparentVortex_MovingWedge &
    ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
      SetTag, ClearTag, TagLo, TagHi, Tag )

    INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                            TagLo(4), TagHi(4)
    REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1), iLo(2):iHi(2), &
                                iLo(3):iHi(3), iLo(4):iHi(4))
    REAL(DP), INTENT(in) :: TagCriteria
    CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
    CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                                 TagLo(2):TagHi(2), &
                                                 TagLo(3):TagHi(3), &
                                                 TagLo(4):TagHi(4))

    REAL(DP), PARAMETER :: X2_MIN       = 0.125_DP * Pi  
    REAL(DP), PARAMETER :: X2_MAX       = 0.375_DP * Pi  
    REAL(DP), PARAMETER :: WEDGE_WIDTH  = 4.0_DP         
    REAL(DP), PARAMETER :: R_LEAD_PAD   = 0.5_DP         
    REAL(DP), PARAMETER :: R_FLOOR      = 2.0_DP         
    REAL(DP), PARAMETER :: R_CEIL       = 10.0_DP

    REAL(DP) :: X1_C, X2_C, R_Lo_Tag, R_Hi_Tag
    INTEGER  :: iX1, iX2, iX3

    R_Hi_Tag = MIN( R_LEAD_PAD, R_CEIL )
    R_Lo_Tag = MAX( R_Hi_Tag - WEDGE_WIDTH, R_FLOOR )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)

      X2_C = MeshX(2) % Center(iX2)

      DO iX1 = iX_B0(1), iX_E0(1)

        Tag(iX1,iX2,iX3,1) = ClearTag

        IF( X2_C .LT. X2_MIN .OR. X2_C .GT. X2_MAX ) CYCLE

        X1_C = MeshX(1) % Center(iX1)

        IF( X1_C .GE. R_Lo_Tag .AND. X1_C .LE. R_Hi_Tag ) &
          Tag(iX1,iX2,iX3,1) = SetTag

      END DO

    END DO
    END DO

  END SUBROUTINE TagElements_TransparentVortex_MovingWedge


SUBROUTINE TagElements_Exterior &
  ( iLevel, iX_B0, iX_E0, iLo, iHi, uCR, TagCriteria, &
    SetTag, ClearTag, TagLo, TagHi, Tag )

  INTEGER,  INTENT(in) :: iLevel, iX_B0(3), iX_E0(3), iLo(4), iHi(4), &
                          TagLo(4), TagHi(4)
  REAL(DP), INTENT(in) :: uCR(iLo(1):iHi(1),iLo(2):iHi(2), &
                              iLo(3):iHi(3),iLo(4):iHi(4))
  REAL(DP), INTENT(in) :: TagCriteria
  CHARACTER(KIND=c_char), INTENT(in)    :: SetTag, ClearTag
  CHARACTER(KIND=c_char), INTENT(inout) :: Tag(TagLo(1):TagHi(1), &
                                               TagLo(2):TagHi(2), &
                                               TagLo(3):TagHi(3), &
                                               TagLo(4):TagHi(4))

  INTEGER  :: iX1, iX2, iX3
  REAL(DP) :: TagCriteria_this

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

END SUBROUTINE TagElements_Exterior


END MODULE TaggingModule
