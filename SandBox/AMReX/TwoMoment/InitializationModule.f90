MODULE InitializationModule
  
  USE amrex_amr_module, ONLY: &
    amrex_init, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init


  USE MyAmrModule,                      ONLY: &
    MyAmrInit


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

CONTAINS


  SUBROUTINE InitializeProgram

    CALL amrex_init()

    CALL amrex_amrcore_init()

    
    CALL MyAmrInit


  END SUBROUTINE InitializeProgram  


END MODULE InitializationModule
