MODULE MakeFineMaskModule

  USE ISO_C_BINDING

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_imultifab

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MakeFineMask

  INTERFACE

    SUBROUTINE amrex_fi_makefinemask &
      ( iMF_Mask, CrseBA, CrseDM, FineBA, iCoarse, iFine ) BIND(c)

        IMPORT
        IMPLICIT NONE

        TYPE(c_ptr)           :: iMF_Mask
        TYPE(c_ptr)   , VALUE :: CrseBA
        TYPE(c_ptr)   , VALUE :: CrseDM
        TYPE(c_ptr)   , VALUE :: FineBA
        INTEGER(c_int), VALUE :: iCoarse
        INTEGER(c_int), VALUE :: iFine

    END SUBROUTINE amrex_fi_makefinemask

  END INTERFACE


CONTAINS


  SUBROUTINE MakeFineMask &
    ( iMF_Mask, CrseBA, CrseDM, FineBA, iCoarse, iFine )

    TYPE(amrex_imultifab), INTENT(inout) :: iMF_Mask
    TYPE(amrex_boxarray) , INTENT(in)    :: CrseBA
    TYPE(amrex_distromap), INTENT(in)    :: CrseDM
    TYPE(amrex_boxarray) , INTENT(in)    :: FineBA
    INTEGER(c_int)       , INTENT(in)    :: iCoarse
    INTEGER(c_int)       , INTENT(in)    :: iFine

    iMF_Mask % owner = .TRUE.
    iMF_Mask % nc    = 1
    iMF_Mask % ng    = 0
    CALL amrex_fi_makefinemask &
           ( iMF_Mask % p, CrseBA % p, CrseDM % p, FineBA % p, &
             iCoarse, iFine )

  END SUBROUTINE MakeFineMask


END MODULE MakeFineMaskModule
