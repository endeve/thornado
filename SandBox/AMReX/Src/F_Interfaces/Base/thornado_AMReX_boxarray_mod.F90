MODULE thornado_amrex_boxarray_module

  USE ISO_C_BINDING
  USE amrex_box_module
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex_boxarray_issame

  INTERFACE
    PURE LOGICAL FUNCTION amrex_fi_boxarray_issame( baa, bab ) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(C_PTR), VALUE, INTENT(in) :: baa
      TYPE(C_PTR), VALUE, INTENT(in) :: bab
    END FUNCTION amrex_fi_boxarray_issame
  END INTERFACE

CONTAINS


  PURE FUNCTION amrex_boxarray_issame( baa, bab ) RESULT(r)
    CLASS(amrex_boxarray), INTENT(in) :: baa
    CLASS(amrex_boxarray), INTENT(in) :: bab
    LOGICAL :: r
    r = amrex_fi_boxarray_issame( baa%p, bab%p )
  END FUNCTION amrex_boxarray_issame

END MODULE thornado_amrex_boxarray_module
