MODULE MF_FieldsModule_Geometry

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    nMaxLevels

  IMPLICIT NONE
  PRIVATE

  ! --- Geometry Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uGF(:)

  PUBLIC :: CreateFields_Geometry_MF
  PUBLIC :: DestroyFields_Geometry_MF

CONTAINS


  SUBROUTINE CreateFields_Geometry_MF

    ALLOCATE( MF_uGF(0:nMaxLevels-1) )

  END SUBROUTINE CreateFields_Geometry_MF


  SUBROUTINE DestroyFields_Geometry_MF

    INTEGER :: iLevel

    DO iLevel = 0, nMaxLevels-1

      CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    END DO

    DEALLOCATE( MF_uGF )

  END SUBROUTINE DestroyFields_Geometry_MF


END MODULE MF_FieldsModule_Geometry
