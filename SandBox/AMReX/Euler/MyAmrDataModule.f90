
MODULE MyAmrDataModule

  USE ISO_C_BINDING
  USE amrex_amr_module
  USE amrex_fort_module, ONLY: &
    amrex_real

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_new, t_old, MF_uGF_new, MF_uCF_new, MF_uPF_new, MF_uAF_new, &
                          MF_uGF_old, MF_uCF_old, MF_uPF_old, MF_uAF_old, &
            flux_reg, StepNo_vec, dt_vec, do_reflux
  PUBLIC :: InitializeDataAMReX, FinalizeDataAMReX

  REAL(amrex_real),         ALLOCATABLE :: t_new(:)
  REAL(amrex_real),         ALLOCATABLE :: t_old(:)

  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uGF_old(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uCF_old(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uPF_old(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uAF_old(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uGF_new(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uCF_new(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uPF_new(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uAF_new(:)

  TYPE(amrex_fluxregister), ALLOCATABLE :: flux_reg(:)

  INTEGER,          ALLOCATABLE, SAVE :: StepNo_vec(:)
  REAL(amrex_real), ALLOCATABLE, SAVE :: dt_vec(:)
  REAL(amrex_real), ALLOCATABLE, SAVE :: t_vec(:)
  LOGICAL                             :: do_reflux  = .TRUE.
  
CONTAINS

  SUBROUTINE InitializeDataAMReX

    ALLOCATE( t_new(0:amrex_max_level) )
    t_new = 0.0_amrex_real

    ALLOCATE( t_old(0:amrex_max_level) )
    t_old = -1.0e100_amrex_real

    ALLOCATE( MF_uGF_old(0:amrex_max_level) )
    ALLOCATE( MF_uCF_old(0:amrex_max_level) )
    ALLOCATE( MF_uPF_old(0:amrex_max_level) )
    ALLOCATE( MF_uAF_old(0:amrex_max_level) )

    ALLOCATE( MF_uGF_new(0:amrex_max_level) )
    ALLOCATE( MF_uCF_new(0:amrex_max_level) )
    ALLOCATE( MF_uPF_new(0:amrex_max_level) )
    ALLOCATE( MF_uAF_new(0:amrex_max_level) )

    ALLOCATE( flux_reg(0:amrex_max_level) )

    ALLOCATE( StepNo_vec(0:amrex_max_level) )
    StepNo_vec = 0

    ALLOCATE( dt_vec(0:amrex_max_level) )
    dt_vec = 1.0e-4_amrex_real

    ALLOCATE( t_vec(0:amrex_max_level) )
    t_vec = 0.0e0_amrex_real

  END SUBROUTINE InitializeDataAMReX


  SUBROUTINE FinalizeDataAMReX

    INTEGER :: iLevel

    DO iLevel = 0, amrex_max_level
      CALL amrex_multifab_destroy( MF_uGF_old(iLevel) )
      CALL amrex_multifab_destroy( MF_uCF_old(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF_old(iLevel) )
      CALL amrex_multifab_destroy( MF_uAF_old(iLevel) )
      CALL amrex_multifab_destroy( MF_uGF_new(iLevel) )
      CALL amrex_multifab_destroy( MF_uCF_new(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF_new(iLevel) )
      CALL amrex_multifab_destroy( MF_uAF_new(iLevel) )
    END DO
    DEALLOCATE( MF_uGF_old )
    DEALLOCATE( MF_uCF_old )
    DEALLOCATE( MF_uPF_old )
    DEALLOCATE( MF_uAF_old )
    DEALLOCATE( MF_uGF_new )
    DEALLOCATE( MF_uCF_new )
    DEALLOCATE( MF_uPF_new )
    DEALLOCATE( MF_uAF_new )

    DO iLevel = 1, amrex_max_level
       CALL amrex_fluxregister_destroy( flux_reg(iLevel) )
    END DO
    DEALLOCATE( flux_reg )

    DEALLOCATE( StepNo_vec )
    DEALLOCATE( t_vec )
    DEALLOCATE( dt_vec )
    DEALLOCATE( t_old )
    DEALLOCATE( t_new )

  END SUBROUTINE FinalizeDataAMReX
  
END MODULE MyAmrDataModule
