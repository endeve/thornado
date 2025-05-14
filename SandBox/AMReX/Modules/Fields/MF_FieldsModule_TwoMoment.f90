MODULE MF_FieldsModule_TwoMoment

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_destroy
  USE thornado_amrex_fluxregister_module, ONLY: &
    amrex_fluxregister, &
    amrex_fluxregister_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE RadiationFieldsModule, ONLY: &
    nCR, &
    DescribeRadiationFields_Conserved, &
    DescribeRadiationFields_Primitive, &
    SetUnitsRadiationFields

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    nMaxLevels

  IMPLICIT NONE
  PRIVATE

  ! --- Conserved Radiation Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCR(:)

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_Permute(:)

  ! --- Primitive Radiation Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPR(:)

  ! --- Auxilliary Radiation Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uAR(:)

  ! --- Integrated Variables ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uGR(:)

  TYPE(amrex_fluxregister), ALLOCATABLE, PUBLIC :: FluxRegister_TwoMoment(:)

  REAL(DP), ALLOCATABLE, PUBLIC :: OffGridFlux_TwoMoment_MF(:,:)

  PUBLIC :: CreateFields_TwoMoment_MF
  PUBLIC :: DestroyFields_TwoMoment_MF

CONTAINS


  SUBROUTINE CreateFields_TwoMoment_MF

    ALLOCATE( MF_uCR(0:nMaxLevels-1) )
    ALLOCATE( MF_Permute(0:nMaxLevels-1) )
    ALLOCATE( MF_uPR(0:nMaxLevels-1) )
    ALLOCATE( MF_uAR(0:nMaxLevels-1) )
    ALLOCATE( MF_uGR(0:nMaxLevels-1) )

    ALLOCATE( FluxRegister_TwoMoment(0:nMaxLevels-1) )

    ALLOCATE( OffGridFlux_TwoMoment_MF(1:2*nCR,0:nMaxLevels-1) )

    CALL SetUnitsRadiationFields

    CALL DescribeRadiationFields_Conserved( amrex_parallel_ioprocessor() )
    CALL DescribeRadiationFields_Primitive( amrex_parallel_ioprocessor() )

  END SUBROUTINE CreateFields_TwoMoment_MF


  SUBROUTINE DestroyFields_TwoMoment_MF

    INTEGER :: iLevel

    DEALLOCATE( OffGridFlux_TwoMoment_MF )

    DO iLevel = 0, nMaxLevels-1

      CALL amrex_fluxregister_destroy( FluxRegister_TwoMoment(iLevel) )

      CALL amrex_multifab_destroy( MF_uPR(iLevel) )
      CALL amrex_multifab_destroy( MF_uAR(iLevel) )
      CALL amrex_multifab_destroy( MF_uCR(iLevel) )
      CALL amrex_multifab_destroy( MF_uGR(iLevel) )
      CALL amrex_multifab_destroy( MF_Permute(iLevel) )

    END DO

    DEALLOCATE( FluxRegister_TwoMoment )

    DEALLOCATE( MF_uPR )
    DEALLOCATE( MF_uAR )
    DEALLOCATE( MF_uCR )
    DEALLOCATE( MF_uGR )
    DEALLOCATE( MF_Permute )

  END SUBROUTINE DestroyFields_TwoMoment_MF


END MODULE MF_FieldsModule_TwoMoment
