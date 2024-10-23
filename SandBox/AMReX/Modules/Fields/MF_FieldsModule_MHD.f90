MODULE MF_FieldsModule_MHD

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

  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem
  USE MagnetofluidFieldsModule, ONLY: &
    nCM, &
    DescribeFields_Conserved, &
    DescribeFields_Primitive, &
    DescribeFields_Auxiliary, &
    DescribeFields_Diagnostic, &
    SetUnitsFields

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    nMaxLevels

  IMPLICIT NONE
  PRIVATE

  ! --- Conserved Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uCM(:)

  ! --- Primitive Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uPM(:)

  ! --- Auxiliary Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uAM(:)

  ! --- Diagnostic Fluid Fields ---

  TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uDM(:)

  TYPE(amrex_fluxregister), ALLOCATABLE, PUBLIC :: FluxRegister_MHD(:)

  REAL(DP), ALLOCATABLE, PUBLIC :: OffGridFlux_MHD_MF(:,:)

  PUBLIC :: CreateFields_MHD_MF
  PUBLIC :: DestroyFields_MHD_MF

CONTAINS


  SUBROUTINE CreateFields_MHD_MF

    ALLOCATE( MF_uCM(0:nMaxLevels-1) )
    ALLOCATE( MF_uPM(0:nMaxLevels-1) )
    ALLOCATE( MF_uAM(0:nMaxLevels-1) )
    ALLOCATE( MF_uDM(0:nMaxLevels-1) )

    ALLOCATE( FluxRegister_MHD(0:nMaxLevels-1) )

    ALLOCATE( OffGridFlux_MHD_MF(1:nCM,0:nMaxLevels-1) )

    CALL SetUnitsFields( TRIM( CoordinateSystem ), &
                              Verbose_Option = amrex_parallel_ioprocessor() )

    CALL DescribeFields_Conserved ( amrex_parallel_ioprocessor() )
    CALL DescribeFields_Primitive ( amrex_parallel_ioprocessor() )
    CALL DescribeFields_Auxiliary ( amrex_parallel_ioprocessor() )
    CALL DescribeFields_Diagnostic( amrex_parallel_ioprocessor() )

  END SUBROUTINE CreateFields_MHD_MF


  SUBROUTINE DestroyFields_MHD_MF

    INTEGER :: iLevel

    DEALLOCATE( OffGridFlux_MHD_MF )

    DO iLevel = 0, nMaxLevels-1

      CALL amrex_fluxregister_destroy( FluxRegister_MHD(iLevel) )

      CALL amrex_multifab_destroy( MF_uDM(iLevel) )
      CALL amrex_multifab_destroy( MF_uAM(iLevel) )
      CALL amrex_multifab_destroy( MF_uPM(iLevel) )
      CALL amrex_multifab_destroy( MF_uCM(iLevel) )

    END DO

    DEALLOCATE( FluxRegister_MHD )

    DEALLOCATE( MF_uDM )
    DEALLOCATE( MF_uAM )
    DEALLOCATE( MF_uPM )
    DEALLOCATE( MF_uCM )

  END SUBROUTINE DestroyFields_MHD_MF


END MODULE MF_FieldsModule_MHD
