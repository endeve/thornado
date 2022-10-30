MODULE MF_FieldsModule_MHD

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy

  ! --- thornado Modules ---

  USE MagnetofluidFieldsModule, ONLY: &
    nCM

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  ! --- Geometry Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uGF(:)

  ! --- Conserved Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCM(:)

  ! --- Primitive Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPM(:)

  ! --- Auxiliary Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uAM(:)

  ! --- Diagnostic Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uDM(:)

  PUBLIC :: CreateFields_MF
  PUBLIC :: DestroyFields_MF


CONTAINS


  SUBROUTINE CreateFields_MF( nLevels )

    INTEGER, INTENT(in) :: nLevels

    ALLOCATE( MF_uGF(0:nLevels-1) )
    ALLOCATE( MF_uCM(0:nLevels-1) )
    ALLOCATE( MF_uPM(0:nLevels-1) )
    ALLOCATE( MF_uAM(0:nLevels-1) )
    ALLOCATE( MF_uDM(0:nLevels-1) )

  END SUBROUTINE CreateFields_MF


  SUBROUTINE DestroyFields_MF( nLevels )

    INTEGER, INTENT(in) :: nLevels

    INTEGER :: iLevel

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uDM(iLevel) )
      CALL amrex_multifab_destroy( MF_uAM(iLevel) )
      CALL amrex_multifab_destroy( MF_uPM(iLevel) )
      CALL amrex_multifab_destroy( MF_uCM(iLevel) )
      CALL amrex_multifab_destroy( MF_uGF(iLevel) )

    END DO

    DEALLOCATE( MF_uDM )
    DEALLOCATE( MF_uAM )
    DEALLOCATE( MF_uPM )
    DEALLOCATE( MF_uCM )
    DEALLOCATE( MF_uGF )

  END SUBROUTINE DestroyFields_MF


END MODULE MF_FieldsModule_MHD
