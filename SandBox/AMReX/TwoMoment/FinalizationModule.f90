MODULE FinalizationModule

  ! --- AMReX Modules ---
  USE amrex_amr_module,     ONLY: &
    amrex_geometry, &
    amrex_geometry_destroy, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_finalize

  ! --- thornado Modules ---
  USE RadiationFieldsModule,                ONLY: &
    DestroyRadiationFields
  USE TwoMoment_TimersModule_Relativistic, ONLY: &
    FinalizeTimers

  ! --- Local Modules ---
  USE MyAmrModule,                 ONLY: &
    nLevels, MyAmrFinalize
  USE MF_TwoMoment_TimeSteppingModule_Relativistic,    ONLY: &
    MF_FinalizeField_IMEX_RK


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FinalizeProgram



 CONTAINS



  SUBROUTINE FinalizeProgram(GEOM)


    TYPE(amrex_geometry),  INTENT(inout) :: GEOM(0:nLevels-1)
    
    INTEGER :: iLevel

    CALL DestroyRadiationFields

    CALL MF_FinalizeField_IMEX_RK

    DO iLevel = 0, nLevels-1
      CALL amrex_geometry_destroy( GEOM(iLevel) )
    END DO  
 
    CALL MyAmrFinalize
 
    CALL FinalizeTimers
 
    CALL amrex_amrcore_finalize

    CALL amrex_finalize


  END SUBROUTINE FinalizeProgram


END MODULE FinalizationModule
