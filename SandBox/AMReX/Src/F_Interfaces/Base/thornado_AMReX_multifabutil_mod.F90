MODULE thornado_amrex_multifabutil_module

  USE ISO_C_BINDING
  USE amrex_fort_module
  USE amrex_multifab_module
  USE amrex_geometry_module

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex_average_down_dg_conservative, &
            amrex_average_down_dg_pointwise, &
            amrex_average_down_cg

  INTERFACE

     SUBROUTINE amrex_fi_average_down_dg_conservative &
       ( FineMF, CrseMF, FineMF_G, CrseMF_G, nComp, RefRatio, &
         nDOFX, nFine, vpFineToCoarseProjectionMatrix ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE :: &
         FineMF, CrseMF, FineMF_G, CrseMF_G, vpFineToCoarseProjectionMatrix
       INTEGER(C_INT)  , VALUE :: nComp, RefRatio, nDOFX, nFine
     END SUBROUTINE amrex_fi_average_down_dg_conservative

     SUBROUTINE amrex_fi_average_down_dg_pointwise &
       ( FineMF, CrseMF, nComp, RefRatio, &
         nDOFX, nFine, vpFineToCoarseProjectionMatrix ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE :: &
         FineMF, CrseMF, vpFineToCoarseProjectionMatrix
       INTEGER(C_INT)  , VALUE :: nComp, RefRatio, nDOFX, nFine
     END SUBROUTINE amrex_fi_average_down_dg_pointwise

     SUBROUTINE amrex_fi_average_down_cg &
       ( FineMF, CrseMF, nComp, RefRatio, &
         nDOFX, nFine, G2L, L2G, F2C ) BIND(c)
       IMPORT
       IMPLICIT NONE
       TYPE(C_PTR)     , VALUE :: &
         FineMF, CrseMF, G2L, L2G, F2C
       INTEGER(C_INT)  , VALUE :: nComp, RefRatio, nDOFX, nFine
     END SUBROUTINE amrex_fi_average_down_cg

  END INTERFACE

CONTAINS


  SUBROUTINE amrex_average_down_dg_conservative &
    ( FineMF, CrseMF, FineMF_G, CrseMF_G, nComp, RefRatio, &
      nDOFX, nFine, vpFineToCoarseProjectionMatrix )

    TYPE(amrex_multifab), INTENT(in)    :: FineMF, FineMF_G, CrseMF_G
    TYPE(amrex_multifab), INTENT(inout) :: CrseMF
    INTEGER             , INTENT(in)    :: nComp, RefRatio, nDOFX, nFine
    TYPE(C_PTR)         , INTENT(in)    :: vpFineToCoarseProjectionMatrix

    CALL amrex_fi_average_down_dg_conservative &
           ( FineMF % p, CrseMF % p, FineMF_G % p, CrseMF_G % p, &
             nComp, RefRatio, nDOFX, nFine, vpFineToCoarseProjectionMatrix )

  END SUBROUTINE amrex_average_down_dg_conservative


  SUBROUTINE amrex_average_down_dg_pointwise &
    ( FineMF, CrseMF, nComp, RefRatio, &
      nDOFX, nFine, vpFineToCoarseProjectionMatrix )

    TYPE(amrex_multifab), INTENT(in)    :: FineMF
    TYPE(amrex_multifab), INTENT(inout) :: CrseMF
    INTEGER             , INTENT(in)    :: nComp, RefRatio, nDOFX, nFine
    TYPE(C_PTR)         , INTENT(in)    :: vpFineToCoarseProjectionMatrix

    CALL amrex_fi_average_down_dg_pointwise &
           ( FineMF % p, CrseMF % p, nComp, RefRatio, &
             nDOFX, nFine, vpFineToCoarseProjectionMatrix )

  END SUBROUTINE amrex_average_down_dg_pointwise


  SUBROUTINE amrex_average_down_cg &
    ( FineMF, CrseMF, nComp, RefRatio, &
      nDOFX, nFine, G2L, L2G, F2C )

    TYPE(amrex_multifab), INTENT(in)    :: FineMF
    TYPE(amrex_multifab), INTENT(inout) :: CrseMF
    INTEGER             , INTENT(in)    :: nComp, RefRatio, nDOFX, nFine
    TYPE(C_PTR)         , INTENT(in)    :: G2L, L2G, F2C

    CALL amrex_fi_average_down_cg &
           ( FineMF % p, CrseMF % p, nComp, RefRatio, &
             nDOFX, nFine, G2L, L2G, F2C )

  END SUBROUTINE amrex_average_down_cg

END MODULE thornado_amrex_multifabutil_module
