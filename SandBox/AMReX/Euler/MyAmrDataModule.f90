
MODULE MyAmrDataModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy

  ! --- Local Modules ---
  USE amrex_amr_module, ONLY: &
    amrex_max_level

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: MF_uGF, MF_uCF, MF_uPF, MF_uAF

  PUBLIC :: InitializeDataAMReX, FinalizeDataAMReX

  TYPE(amrex_multifab), ALLOCATABLE :: MF_uGF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uCF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uPF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uAF(:)

  
CONTAINS


  SUBROUTINE InitializeDataAMReX( amrex_max_level )

    INTEGER, INTENT(in) :: amrex_max_level

    ALLOCATE( MF_uGF(0:amrex_max_level) )
    ALLOCATE( MF_uCF(0:amrex_max_level) )
    ALLOCATE( MF_uPF(0:amrex_max_level) )
    ALLOCATE( MF_uAF(0:amrex_max_level) )

  END SUBROUTINE InitializeDataAMReX


  SUBROUTINE FinalizeDataAMReX( amrex_max_level )

    INTEGER, INTENT(in) :: amrex_max_level

    INTEGER :: iLevel

    DO iLevel = 0, amrex_max_level
      CALL amrex_multifab_destroy( MF_uAF(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF(iLevel) )
      CALL amrex_multifab_destroy( MF_uCF(iLevel) )
      CALL amrex_multifab_destroy( MF_uGF(iLevel) )
    END DO

    DEALLOCATE( MF_uAF )
    DEALLOCATE( MF_uPF )
    DEALLOCATE( MF_uCF )
    DEALLOCATE( MF_uGF )

  END SUBROUTINE FinalizeDataAMReX
  
END MODULE MyAmrDataModule
