MODULE MF_FieldsModule_TwoMoment

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy
  USE amrex_fluxregister_module, ONLY: &
    amrex_fluxregister, &
    amrex_fluxregister_destroy

  ! --- thornado Modules ---

  USE RadiationFieldsModule, ONLY: &
    nCR

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    nMaxLevels

  IMPLICIT NONE
  PRIVATE

  ! --- Conserved Radiation Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCR(:)

  ! --- Primitive Radiation Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPR(:)

  TYPE(amrex_fluxregister), ALLOCATABLE, PUBLIC :: FluxRegister_TwoMoment(:)

  REAL(DP), ALLOCATABLE, PUBLIC :: OffGridFlux_TwoMoment_MF(:,:)

  PUBLIC :: CreateFields_TwoMoment_MF
  PUBLIC :: DestroyFields_TwoMoment_MF

CONTAINS


  SUBROUTINE CreateFields_TwoMoment_MF

    ALLOCATE( MF_uCR(0:nMaxLevels-1) )
    ALLOCATE( MF_uPR(0:nMaxLevels-1) )

    ALLOCATE( FluxRegister_TwoMoment(0:nMaxLevels-1) )

    ALLOCATE( OffGridFlux_TwoMoment_MF(1:nCR,0:nMaxLevels-1) )

  END SUBROUTINE CreateFields_TwoMoment_MF


  SUBROUTINE DestroyFields_TwoMoment_MF

    INTEGER :: iLevel

    DEALLOCATE( OffGridFlux_TwoMoment_MF )

    DO iLevel = 0, nMaxLevels-1

      CALL amrex_fluxregister_destroy( FluxRegister_TwoMoment(iLevel) )

      CALL amrex_multifab_destroy( MF_uPR(iLevel) )
      CALL amrex_multifab_destroy( MF_uCR(iLevel) )

    END DO

    DEALLOCATE( FluxRegister_TwoMoment )

    DEALLOCATE( MF_uPR )
    DEALLOCATE( MF_uCR )

  END SUBROUTINE DestroyFields_TwoMoment_MF


END MODULE MF_FieldsModule_TwoMoment
