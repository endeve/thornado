MODULE ReGridModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor
  USE amrex_amrcore_module, ONLY: &
   amrex_regrid, &
   amrex_get_numlevels
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE XCFC_UtilitiesModule, ONLY: &
    nGS, &
    nMF, &
    swX_GS

  ! --- Local Modules ---

  USE FillPatchModule_MHD, ONLY: &
    FillPatch
  USE AverageDownModule_MHD, ONLY: &
    AverageDown
  USE MF_MHD_UtilitiesModule, ONLY: &
    ComputeFromConserved_MHD_MF
  USE MF_MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_MHD, ONLY: &
    MF_uCM, &
    MF_uPM, &
    MF_uAM, &
    MF_uDM
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_uCR
  USE MF_MHD_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_MHD_MF
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE MF_GravitySolutionModule, ONLY: &
    EvolveGravity
  USE MF_GravitySolutionModule_XCFC, ONLY: &
    ComputeConformalFactor_XCFC_MF, &
    ComputeLapseShiftCurvature_XCFC_MF
  USE MF_XCFC_UtilitiesModule, ONLY: &
    MultiplyWithPsi6_MF, &
    UpdateConformalFactorAndMetric_XCFC_MF, &
    UpdateLapseShiftCurvature_XCFC_MF, &
    ComputeConformalFactorSourcesAndMg_XCFC_MF, &
    ComputePressureTensorTrace_XCFC_MF
  USE MF_GravitySolutionModule_Newtonian_Poseidon, ONLY: &
    ComputeGravitationalPotential_Newtonian_MF_Poseidon
  USE InputParsingModule, ONLY: &
    DEBUG, &
    UseAMR, &
    StepNo, &
    iReGrid, &
    nLevels, &
    nMaxLevels, &
    t_old, &
    t_new
  USE MF_KindModule, ONLY: &
    Zero

  PUBLIC :: ReGrid

CONTAINS


  SUBROUTINE ReGrid

    TYPE(amrex_multifab) :: MF_uGS(0:nMaxLevels-1)
    TYPE(amrex_multifab) :: MF_uMF(0:nMaxLevels-1)

    INTEGER :: iLevel, iErr

    IF( .NOT. UseAMR ) RETURN

    IF( .NOT. PerformRegrid ) RETURN

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(6x,A,I2.2)') 'nLevels (before regrid): ', nLevels
        WRITE(*,'(6x,A)') 'Regridding'

      END IF

    END IF

    CALL MultiplyWithPsi6_MF &
           ( MF_uGF, MF_uCM, +1, OnlyLeafElements_Option = .FALSE. )

    CALL amrex_regrid( 0, t_new(0) )

    nLevels = amrex_get_numlevels()

    CALL MultiplyWithPsi6_MF &
           ( MF_uGF, MF_uCM, -1, OnlyLeafElements_Option = .FALSE. )

    CALL ApplyPositivityLimiter_MHD_MF &
           ( t_old(0), MF_uGF, MF_uCM, MF_uDM )

    CALL MultiplyWithPsi6_MF &
           ( MF_uGF, MF_uCM, +1, OnlyLeafElements_Option = .FALSE. )

    IF( EvolveGravity )THEN

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_build &
               ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
                 nDOFX * nGS, 0 )
        CALL MF_uGS(iLevel) % SetVal( Zero ) ! remove this after debugging

        CALL amrex_multifab_build &
               ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
                 nDOFX * nMF, swX_GS )
        CALL MF_uMF(iLevel) % SetVal( Zero ) ! remove this after debugging

      END DO

#ifndef THORNADO_NOTRANSPORT

      CALL ComputeConformalFactor &
             ( MF_uGF, MF_uCM, MF_uCR, MF_uGS, MF_uMF )

      CALL ComputeLapseShiftCurvature &
             ( MF_uGF, MF_uCM, MF_uCR, MF_uGS, MF_uMF )

#else

      CALL ComputeConformalFactor &
             ( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

      CALL ComputeLapseShiftCurvature &
             ( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

#endif

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_destroy( MF_uMF(iLevel) )
        CALL amrex_multifab_destroy( MF_uGS(iLevel) )

      END DO

#else

      CALL ComputeGravitationalPotential_Newtonian_MF_Poseidon( MF_uCM, MF_uGF )

#endif

    END IF ! EvolveGravity

    CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCM, -1 )

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )
    CALL AverageDown( MF_uGF, MF_uCM )

    CALL ApplyPositivityLimiter_MHD_MF &
           ( t_old(0), MF_uGF, MF_uCM, MF_uDM )

    CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

    CALL ApplyBoundaryConditions_MHD_MF( t_old(0), MF_uCM )

    ! --- nLevels <= nMaxLevels; entire arrays t_old(0:nMaxLevels-1) and
    !     t_new(0:nMaxLevels-1) must have valid data ---
    t_old = t_old(0)
    t_new = t_new(0)

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,*)
        WRITE(*,'(A)') 'CALL ComputeFromConserved_MHD_MF'
        WRITE(*,*)

      END IF

      CALL ComputeFromConserved_MHD_MF &
             ( MF_uGF, MF_uCM, MF_uPM, MF_uAM )

    END IF

  END SUBROUTINE ReGrid


  ! --- PRIVATE SUBROUTINES ---


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE ComputeConformalFactor &
    ( MF_uGF, MF_uCM, MF_uCR, MF_uGS, MF_uMF )

#else

  SUBROUTINE ComputeConformalFactor &
    ( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifndef THORNADO_NOTRANSPORT

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCM, MF_uCR, MF_uGS )

#else

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCM, MF_uGS )

#endif

    CALL ComputeConformalFactor_XCFC_MF &
           ( MF_uGS, MF_uMF )

    CALL UpdateConformalFactorAndMetric_XCFC_MF &
           ( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )

  END SUBROUTINE ComputeConformalFactor


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE ComputeLapseShiftCurvature &
    ( MF_uGF, MF_uCM, MF_uCR, MF_uGS, MF_uMF )

#else

  SUBROUTINE ComputeLapseShiftCurvature &
    ( MF_uGF, MF_uCM, MF_uGS, MF_uMF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifndef THORNADO_NOTRANSPORT

    CALL ComputePressureTensorTrace_XCFC_MF &
           ( MF_uGF, MF_uCM, MF_uCR, MF_uGS )

#else

    CALL ComputePressureTensorTrace_XCFC_MF &
           ( MF_uGF, MF_uCM, MF_uGS )
#endif

    CALL ComputeLapseShiftCurvature_XCFC_MF &
           ( MF_uGS, MF_uMF )

    CALL UpdateLapseShiftCurvature_XCFC_MF &
           ( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF, UpdateSpatialMetric_Option = .TRUE. )

    CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

  END SUBROUTINE ComputeLapseShiftCurvature


  LOGICAL FUNCTION PerformReGrid()

    PerformRegrid = .FALSE.

    IF( MOD( StepNo(0), iReGrid ) .EQ. 0 ) PerformRegrid = .TRUE.

    RETURN
  END FUNCTION PerformReGrid


END MODULE ReGridModule
