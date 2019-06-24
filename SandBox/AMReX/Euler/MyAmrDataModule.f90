MODULE MyAmrDataModule

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy

  IMPLICIT NONE
  PRIVATE

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uGF(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCF(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPF(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uAF(:)

  PUBLIC :: InitializeDataAMReX
  PUBLIC :: FinalizeDataAMReX

  
CONTAINS


  SUBROUTINE InitializeDataAMReX( nLevels )

    INTEGER, INTENT(in) :: nLevels

    ALLOCATE( MF_uGF(0:nLevels) )
    ALLOCATE( MF_uCF(0:nLevels) )
    ALLOCATE( MF_uPF(0:nLevels) )
    ALLOCATE( MF_uAF(0:nLevels) )

  END SUBROUTINE InitializeDataAMReX


  SUBROUTINE FinalizeDataAMReX( nLevels )

    INTEGER, INTENT(in) :: nLevels

    INTEGER :: iLevel

    DO iLevel = 0, nLevels
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
