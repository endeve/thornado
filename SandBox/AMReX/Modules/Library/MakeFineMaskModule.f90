MODULE MakeFineMaskModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_imultifab, &
    amrex_imultifab_build, &
    amrex_imultifab_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- Local Modules ---

  USE InputParsingModule, ONLY: &
    nLevels, &
    swX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MakeFineMask
  PUBLIC :: DestroyFineMask

  INTEGER, PUBLIC :: iLeaf_MFM  = 0
  INTEGER         :: iTrunk_MFM = 1

  INTERFACE

    SUBROUTINE amrex_fi_makefinemask_thornado &
      ( iMF_Mask, CrseBA, CrseDM, FineBA, iCoarse, iFine, swXX, geom ) BIND(c)

        IMPORT
        IMPLICIT NONE

        TYPE(c_ptr)           :: iMF_Mask
        TYPE(c_ptr)   , VALUE :: CrseBA
        TYPE(c_ptr)   , VALUE :: CrseDM
        TYPE(c_ptr)   , VALUE :: FineBA
        INTEGER(c_int), VALUE :: iCoarse
        INTEGER(c_int), VALUE :: iFine
        INTEGER(c_int), VALUE :: swXX
        TYPE(c_ptr)   , VALUE :: geom

    END SUBROUTINE amrex_fi_makefinemask_thornado


    SUBROUTINE amrex_fi_destroyfinemask_thornado( iMF_Mask ) BIND(c)

        IMPORT
        IMPLICIT NONE

        TYPE(c_ptr) :: iMF_Mask

    END SUBROUTINE amrex_fi_destroyfinemask_thornado

  END INTERFACE

CONTAINS


  SUBROUTINE MakeFineMask( iLevel, iMF_Mask, BA, DM )

    INTEGER              , INTENT(in)    :: iLevel
    TYPE(amrex_imultifab), INTENT(inout) :: iMF_Mask
    TYPE(amrex_boxarray) , INTENT(in)    :: BA(0:)
    TYPE(amrex_distromap), INTENT(in)    :: DM(0:)

    IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 )THEN

      ! This uses swX(1) because amrex uses IntVect(int), which creates
      ! a vector of dimension `nDimsX` and fills each value with `int`
      iMF_Mask % owner = .TRUE.
      iMF_Mask % nc    = 1
      CALL amrex_fi_makefinemask_thornado &
             ( iMF_Mask % p, BA(iLevel) % p, DM(iLevel) % p, BA(iLevel+1) % p, &
               iLeaf_MFM, iTrunk_MFM, swX(1), amrex_geom(iLevel) % p )

    ELSE

      CALL amrex_imultifab_build( iMF_Mask, BA(iLevel), DM(iLevel), 1, swX )
      CALL iMF_Mask % SetVal( iLeaf_MFM )

    END IF

  END SUBROUTINE MakeFineMask


  SUBROUTINE DestroyFineMask( iLevel, iMF_Mask )

    INTEGER              , INTENT(in)    :: iLevel
    TYPE(amrex_imultifab), INTENT(inout) :: iMF_Mask

    IF( nLevels .GT. 1 .AND. iLevel .LT. nLevels-1 )THEN

      CALL amrex_fi_destroyfinemask_thornado( iMF_Mask % p )
      CALL amrex_imultifab_destroy( iMF_Mask )

    ELSE

      CALL amrex_fi_destroyfinemask_thornado( iMF_Mask % p )
      CALL amrex_imultifab_destroy( iMF_Mask )

    END IF

  END SUBROUTINE DestroyFineMask

END MODULE MakeFineMaskModule
