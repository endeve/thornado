MODULE FinalizationModule

  ! --- AMReX Modules ---
  USE amrex_amr_module,     ONLY: &
    amrex_geometry, &
    amrex_geometry_destroy, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram



 CONTAINS



  SUBROUTINE FinalizeProgram


    !TYPE(amrex_geometry),  INTENT(inout) :: GEOM(0:nLevels-1)
    !TYPE(MeshType),        INTENT(inout) :: MeshX(1:3)
    
    CALL amrex_amrcore_finalize

    CALL amrex_finalize


  END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
