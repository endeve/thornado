MODULE MF_FieldsModule

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy

  IMPLICIT NONE
  PRIVATE

  ! --- Geometry Fields ---
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uGF(:)

  ! --- Conserved Fluid Fields ---
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCF(:)

  ! --- Primitive Fluid Fields ---
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPF(:)

  ! --- Auxiliary Fluid Fields ---
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uAF(:)

  ! --- Diagnostic Fields ---
  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uDF(:)

  PUBLIC :: CreateFields_MF
  PUBLIC :: DestroyFields_MF


CONTAINS


  SUBROUTINE CreateFields_MF( nLevels )

    INTEGER, INTENT(in) :: nLevels

    ALLOCATE( MF_uGF(0:nLevels-1) )
    ALLOCATE( MF_uCF(0:nLevels-1) )
    ALLOCATE( MF_uPF(0:nLevels-1) )
    ALLOCATE( MF_uAF(0:nLevels-1) )
    ALLOCATE( MF_uDF(0:nLevels-1) )

  END SUBROUTINE CreateFields_MF


  SUBROUTINE DestroyFields_MF( nLevels )

    INTEGER, INTENT(in) :: nLevels

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uDF(iLevel) )
      CALL amrex_multifab_destroy( MF_uAF(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF(iLevel) )
      CALL amrex_multifab_destroy( MF_uCF(iLevel) )
      CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    END DO

    DEALLOCATE( MF_uDF )
    DEALLOCATE( MF_uAF )
    DEALLOCATE( MF_uPF )
    DEALLOCATE( MF_uCF )
    DEALLOCATE( MF_uGF )

  END SUBROUTINE DestroyFields_MF


END MODULE MF_FieldsModule
