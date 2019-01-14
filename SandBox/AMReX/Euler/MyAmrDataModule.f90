
MODULE MyAmrDataModule

  USE ISO_C_BINDING
  USE amrex_amr_module
  USE amrex_fort_module, ONLY: &
    amrex_real

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_new, MF_uGF, MF_uCF, MF_uPF, MF_uAF, &
            flux_reg, StepNo_vec, dt_vec, do_reflux
  PUBLIC :: InitializeDataAMReX, FinalizeDataAMReX

  REAL(amrex_real),         ALLOCATABLE :: t_new(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uGF(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uCF(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uPF(:)
  TYPE(amrex_multifab),     ALLOCATABLE :: MF_uAF(:)
  TYPE(amrex_fluxregister), ALLOCATABLE :: flux_reg(:)

  INTEGER,          ALLOCATABLE, SAVE :: StepNo_vec(:)
  REAL(amrex_real), ALLOCATABLE, SAVE :: dt_vec(:)
  LOGICAL                             :: do_reflux  = .TRUE.
  
CONTAINS

  SUBROUTINE InitializeDataAMReX

    WRITE(*,*)
    WRITE(*,'(A)') 'Calling InitializeDataAMReX...'
    WRITE(*,*)

    ALLOCATE( t_new(0:amrex_max_level) )
    t_new = 0.0_amrex_real

    ALLOCATE( MF_uGF(0:amrex_max_level) )
    ALLOCATE( MF_uCF(0:amrex_max_level) )
    ALLOCATE( MF_uPF(0:amrex_max_level) )
    ALLOCATE( MF_uAF(0:amrex_max_level) )

    ALLOCATE( flux_reg(0:amrex_max_level) )

    ALLOCATE( StepNo_vec(0:amrex_max_level) )
    StepNo_vec = 0

    ALLOCATE( dt_vec(0:amrex_max_level) )
    dt_vec = 1.0e-4_amrex_real

  END SUBROUTINE InitializeDataAMReX


  SUBROUTINE FinalizeDataAMReX

    INTEGER :: iLevel

    WRITE(*,*)
    WRITE(*,'(A)') 'Calling FinalizeDataAMReX...'
    WRITE(*,*)

    DO iLevel = 0, amrex_max_level
      CALL amrex_multifab_destroy( MF_uGF(iLevel) )
      CALL amrex_multifab_destroy( MF_uCF(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF(iLevel) )
      CALL amrex_multifab_destroy( MF_uAF(iLevel) )
    END DO
    DEALLOCATE( MF_uGF )
    DEALLOCATE( MF_uCF )
    DEALLOCATE( MF_uPF )
    DEALLOCATE( MF_uAF )

    DO iLevel = 1, amrex_max_level
       CALL amrex_fluxregister_destroy( flux_reg(iLevel) )
    END DO
    DEALLOCATE( flux_reg )

    DEALLOCATE( StepNo_vec )
    DEALLOCATE( dt_vec )
    DEALLOCATE( t_new )

  END SUBROUTINE FinalizeDataAMReX
  
END MODULE MyAmrDataModule
