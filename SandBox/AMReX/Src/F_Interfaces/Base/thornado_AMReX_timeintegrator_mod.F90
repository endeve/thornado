MODULE thornado_amrex_timeintegrator_module

  USE ISO_C_BINDING
  USE amrex_fort_module
  USE amrex_multifab_module
  USE amrex_geometry_module

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex_integrate_twomoment

  INTERFACE

    ! FIXME: declare interface
    SUBROUTINE amrex_fi_integrate_twomoment &
        ( GEOM, GE, MF_uGF, MF_uCF, MF_uCR, t, dt, &
          rhs_collision_p, rhs_streaming_p, &
          swX, nGF, nCF, nCR, nDOFX, nDOFE, nSpecies, iE_E0, iE_B0, iE_E1, iE_B1) BIND(c)
      IMPORT
      IMPLICIT NONE
      TYPE(c_ptr), VALUE                :: GEOM, MF_uGF, MF_uCF, MF_uCR
      INTEGER(c_int),     INTENT(in)    :: swX(3)
      INTEGER(c_int),     VALUE         :: nGF, nCF, nCR, nDOFX, nDOFE, nSpecies, iE_E0, iE_B0, iE_E1, iE_B1
      REAL(amrex_real),   INTENT(in)    :: GE(*)
      REAL(amrex_real), VALUE           :: t, dt
      TYPE(c_funptr), VALUE             :: rhs_collision_p, rhs_streaming_p
    END SUBROUTINE amrex_fi_integrate_twomoment

  END INTERFACE

CONTAINS

  SUBROUTINE amrex_integrate_twomoment &
    ( GEOM, GE, MF_uGF, MF_uCF, MF_uCR, t_new, dt, rhs_collision_p, rhs_streaming_p, &
      swX, nGF, nCF, nCR, nDOFX, nDOFE, nSpecies, iE_E0, iE_B0, iE_E1, iE_B1 )

    TYPE(amrex_geometry), INTENT(in)    :: GEOM
    REAL(amrex_real),     INTENT(in)    :: GE(*)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR
    REAL(amrex_real),     INTENT(in)    :: t_new, dt
    INTEGER(c_int),       INTENT(in)    :: swX(3), nGF, nCF, nCR, nDOFX, nDOFE, nSpecies, iE_E0, iE_B0, iE_E1, iE_B1
    TYPE(c_funptr), value               :: rhs_collision_p, rhs_streaming_p

    CALL amrex_fi_integrate_twomoment &
             ( GEOM%p, GE, &
               MF_uGF%p, MF_uCF%p, MF_uCR%p, &
               t_new, dt, &
               rhs_collision_p, &
               rhs_streaming_p, &
               swX, nGF, nCF, nCR, nDOFX, nDOFE, nSpecies, iE_E0, iE_B0, iE_E1, iE_B1 &
               )

  END SUBROUTINE amrex_integrate_twomoment

END MODULE thornado_amrex_timeintegrator_module
