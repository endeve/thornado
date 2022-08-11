MODULE MF_FieldsModule_Euler

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy
  USE amrex_fluxregister_module, ONLY: &
    amrex_fluxregister, &
    amrex_fluxregister_destroy

  ! --- thornado Modules ---

  USE FluidFieldsModule, ONLY: &
    nCF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    nMaxLevels

  IMPLICIT NONE
  PRIVATE

  ! --- Conserved Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCF(:)

  ! --- Primitive Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPF(:)

  ! --- Auxiliary Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uAF(:)

  ! --- Diagnostic Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uDF(:)

  TYPE(amrex_fluxregister), ALLOCATABLE, PUBLIC :: FluxRegister_Euler(:)

  REAL(DP), ALLOCATABLE, PUBLIC :: OffGridFlux_Euler_MF(:,:)

  PUBLIC :: CreateFields_Euler_MF
  PUBLIC :: DestroyFields_Euler_MF

CONTAINS


  SUBROUTINE CreateFields_Euler_MF

    ALLOCATE( MF_uCF(0:nMaxLevels-1) )
    ALLOCATE( MF_uPF(0:nMaxLevels-1) )
    ALLOCATE( MF_uAF(0:nMaxLevels-1) )
    ALLOCATE( MF_uDF(0:nMaxLevels-1) )

    ALLOCATE( FluxRegister_Euler(0:nMaxLevels-1) )

    ALLOCATE( OffGridFlux_Euler_MF(1:nCF,0:nMaxLevels-1) )

  END SUBROUTINE CreateFields_Euler_MF


  SUBROUTINE DestroyFields_Euler_MF

    INTEGER :: iLevel

    DEALLOCATE( OffGridFlux_Euler_MF )

    DO iLevel = 0, nMaxLevels-1

      CALL amrex_fluxregister_destroy( FluxRegister_Euler(iLevel) )

      CALL amrex_multifab_destroy( MF_uDF(iLevel) )
      CALL amrex_multifab_destroy( MF_uAF(iLevel) )
      CALL amrex_multifab_destroy( MF_uPF(iLevel) )
      CALL amrex_multifab_destroy( MF_uCF(iLevel) )

    END DO

    DEALLOCATE( FluxRegister_Euler )

    DEALLOCATE( MF_uDF )
    DEALLOCATE( MF_uAF )
    DEALLOCATE( MF_uPF )
    DEALLOCATE( MF_uCF )

  END SUBROUTINE DestroyFields_Euler_MF


END MODULE MF_FieldsModule_Euler
