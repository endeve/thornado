MODULE thornado_amrex_distromap_module

  USE iso_c_binding
  USE amrex_boxarray_module
  USE amrex_distromap_module, ONLY: &
    amrex_distromap

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex_distromap_issame

  INTERFACE
     PURE LOGICAL FUNCTION amrex_fi_distromap_issame( dma, dmb ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR), VALUE, INTENT(in) :: dma
       TYPE(C_PTR), VALUE, INTENT(in) :: dmb
     END FUNCTION amrex_fi_distromap_issame
  END INTERFACE

CONTAINS

 PURE FUNCTION amrex_distromap_issame( dma, dmb ) RESULT(r)
    TYPE(amrex_distromap), INTENT(in) :: dma
    TYPE(amrex_distromap), INTENT(in) :: dmb
    LOGICAL :: r
    r =  amrex_fi_distromap_issame( dma%p, dmb%p )
  END FUNCTION amrex_distromap_issame

END MODULE thornado_amrex_distromap_module
