
MODULE MyAmrDataModule

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy

  IMPLICIT NONE
  PRIVATE

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPR(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCR(:)

  PUBLIC :: InitializeDataAMReX
  PUBLIC :: FinalizeDataAMReX

CONTAINS


  SUBROUTINE InitializeDataAMReX( nLevels )

    INTEGER, INTENT(in) :: nLevels

    ALLOCATE( MF_uPR(0:nLevels-1) )
    ALLOCATE( MF_uCR(0:nLevels-1) )

  END SUBROUTINE InitializeDataAMReX


  SUBROUTINE FinalizeDataAMReX( nLevels )

    INTEGER, INTENT(in) :: nLevels

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1
      CALL amrex_multifab_destroy( MF_uPR(iLevel) )
      CALL amrex_multifab_destroy( MF_uCR(iLevel) )
    END DO

    DEALLOCATE( MF_uPR )
    DEALLOCATE( MF_uCR )

  END SUBROUTINE FinalizeDataAMReX



END MODULE MyAmrDataModule
