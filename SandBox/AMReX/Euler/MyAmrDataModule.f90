
MODULE MyAmrDataModule

  USE ISO_C_BINDING
  USE amrex_amr_module
  USE amrex_fort_module, ONLY: &
    amrex_real

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_new, flux_reg, stepno_vec, dt_vec, do_reflux, &
            MF_uGF, MF_uCF, MF_uPF, MF_uAF
  PUBLIC :: amr_data_init, amr_data_finalize

  REAL(amrex_real), ALLOCATABLE :: t_new(:)

  TYPE(amrex_multifab), ALLOCATABLE :: MF_uGF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uCF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uPF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uAF(:)

  TYPE(amrex_fluxregister), ALLOCATABLE :: flux_reg(:)

  INTEGER, ALLOCATABLE, SAVE :: stepno_vec(:)
  REAL(amrex_real), ALLOCATABLE, SAVE :: dt_vec(:)
  LOGICAL :: do_reflux  = .TRUE.
  
CONTAINS

  SUBROUTINE amr_data_init()
    ALLOCATE(t_new(0:amrex_max_level))
    t_new = 0.0_amrex_real

    ALLOCATE(MF_uGF(0:amrex_max_level))
    ALLOCATE(MF_uCF(0:amrex_max_level))
    ALLOCATE(MF_uPF(0:amrex_max_level))
    ALLOCATE(MF_uAF(0:amrex_max_level))

    ALLOCATE(flux_reg(0:amrex_max_level))

    ALLOCATE(stepno_vec(0:amrex_max_level))
    stepno_vec = 0

    ALLOCATE(dt_vec(0:amrex_max_level))
    dt_vec = 1.0e-4_amrex_real
  END SUBROUTINE amr_data_init

  SUBROUTINE amr_data_finalize
    INTEGER :: iLevel
    DO iLevel = 0, amrex_max_level
       CALL amrex_multifab_destroy(MF_uGF(iLevel))
       CALL amrex_multifab_destroy(MF_uCF(iLevel))
       CALL amrex_multifab_destroy(MF_uPF(iLevel))
       CALL amrex_multifab_destroy(MF_uAF(iLevel))
    END DO
    DO iLevel = 1, amrex_max_level
       CALL amrex_fluxregister_destroy(flux_reg(iLevel))
    END DO
  END SUBROUTINE amr_data_finalize
  
END MODULE MyAmrDataModule
