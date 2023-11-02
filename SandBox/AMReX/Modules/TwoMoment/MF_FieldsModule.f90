
MODULE MF_FieldsModule

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy

  IMPLICIT NONE
  PRIVATE

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCF(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPF(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uGF(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uAF(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uDF(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPR(:)
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCR(:)
  REAL(DP), ALLOCATABLE, PUBLIC :: MF_OffGridFlux_Twomoment(:,:)

  PUBLIC :: InitializeDataAMReX
  PUBLIC :: FinalizeDataAMReX

CONTAINS


  SUBROUTINE InitializeDataAMReX( nLevels )

    INTEGER, INTENT(in) :: nLevels

    ALLOCATE( MF_uCF(0:nLevels-1) )
    ALLOCATE( MF_uPF(0:nLevels-1) )
    ALLOCATE( MF_uGF(0:nLevels-1) )
    ALLOCATE( MF_uAF(0:nLevels-1) )
    ALLOCATE( MF_uDF(0:nLevels-1) )
    ALLOCATE( MF_uPR(0:nLevels-1) )
    ALLOCATE( MF_uCR(0:nLevels-1) )
    ALLOCATE( MF_OffGridFlux_TwoMoment(0:nLevels-1,2*nCR) )

  END SUBROUTINE InitializeDataAMReX


  SUBROUTINE FinalizeDataAMReX( nLevels )

    INTEGER, INTENT(in) :: nLevels

    INTEGER :: iLevel
      

    DEALLOCATE( MF_OffGridFlux_TwoMoment )
    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uCF(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF(iLevel) )
      CALL amrex_multifab_destroy( MF_uGF(iLevel) )
      CALL amrex_multifab_destroy( MF_uAF(iLevel) )
      CALL amrex_multifab_destroy( MF_uDF(iLevel) )
      CALL amrex_multifab_destroy( MF_uPR(iLevel) )
      CALL amrex_multifab_destroy( MF_uCR(iLevel) )

    END DO

    DEALLOCATE( MF_uCF )
    DEALLOCATE( MF_uPF )
    DEALLOCATE( MF_uGF )
    DEALLOCATE( MF_uAF )
    DEALLOCATE( MF_uDF )
    DEALLOCATE( MF_uPR )
    DEALLOCATE( MF_uCR )

  END SUBROUTINE FinalizeDataAMReX



END MODULE MF_FieldsModule
