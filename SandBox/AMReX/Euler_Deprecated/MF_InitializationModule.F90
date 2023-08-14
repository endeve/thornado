MODULE MF_InitializationModule

  ! --- AMReX Modules ---

  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module,                         ONLY: &
    amrex_multifab

  ! --- Local Modules ---

  USE InputParsingModule,                            ONLY: &
    nLevels
  USE MF_InitializationModule_NonRelativistic_TABLE, ONLY: &
    MF_InitializeFields_NonRelativistic_TABLE
  USE MF_InitializationModule_NonRelativistic_IDEAL, ONLY: &
    MF_InitializeFields_NonRelativistic_IDEAL
  USE MF_InitializationModule_Relativistic_IDEAL,    ONLY: &
    MF_InitializeFields_Relativistic_IDEAL
  USE MF_InitializationModule_NonRelativistic_IDEAL, ONLY: &
    MF_InitializeFields_NonRelativistic_IDEAL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeFields


CONTAINS


  SUBROUTINE MF_InitializeFields( ProgramName, MF_uGF, MF_uCF, GEOM )

    CHARACTER(LEN=*),     INTENT(in)    :: ProgramName
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)

#if defined HYDRO_NONRELATIVISTIC && defined MICROPHYSICS_WEAKLIB

    CALL MF_InitializeFields_NonRelativistic_TABLE &
           ( ProgramName, MF_uGF, MF_uCF, GEOM )

#elif defined HYDRO_RELATIVISTIC

    CALL MF_InitializeFields_Relativistic_IDEAL &
           ( ProgramName, MF_uGF, MF_uCF, GEOM )

#else

    CALL MF_InitializeFields_NonRelativistic_IDEAL &
           ( ProgramName, MF_uGF, MF_uCF, GEOM )

#endif

  END SUBROUTINE MF_InitializeFields


END MODULE MF_InitializationModule
