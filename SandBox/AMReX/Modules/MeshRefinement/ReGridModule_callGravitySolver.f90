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
    nMF

  ! --- Local Modules ---

  USE AverageDownModule, ONLY: &
    AverageDown
  USE MF_GeometryModule, ONLY: &
    ComputeGeometryX_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF, &
    MF_uDF
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_Euler_MF
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE MF_GravitySolutionModule_XCFC, ONLY: &
    ComputeConformalFactor_XCFC_MF, &
    ComputeLapseShiftCurvature_XCFC_MF, &
    EvolveGravity
  USE MF_XCFC_UtilitiesModule, ONLY: &
    swXX, &
    MultiplyWithPsi6_MF, &
    UpdateConformalFactorAndMetric_XCFC_MF, &
    UpdateLapseShiftCurvature_XCFC_MF, &
    ApplyBoundaryConditions_Geometry_XCFC_MF, &
    ComputeConformalFactorSourcesAndMg_XCFC_MF, &
    ComputePressureTensorTrace_XCFC_MF
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

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

    END IF

    IF( UseAMR )THEN

      IF( MOD( StepNo(0), iReGrid ) .EQ. 0 )THEN

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,*)
            WRITE(*,'(6x,A,I2.2)') 'nLevels (before regrid): ', nLevels
            WRITE(*,'(6x,A)') 'Regridding'

          END IF

        END IF

        CALL amrex_regrid( 0, t_new(0) )

        nLevels = amrex_get_numlevels()

        CALL ComputeGeometryX_MF( MF_uGF )

!        CALL AverageDown( MF_uGF )

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,'(6x,A,I2.2)') 'nLevels (after regrid): ', nLevels
            WRITE(*,*)
            WRITE(*,'(A)') 'CALL ApplyBoundaryConditions_Geometry_MF'

          END IF

        END IF

!        CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,*)
            WRITE(*,'(A)') 'CALL ApplyPositivityLimiter_Euler_MF'
            WRITE(*,*)

          END IF

        END IF

        CALL ApplyPositivityLimiter_Euler_MF &
               ( MF_uGF, MF_uCF, MF_uDF )

        IF( EvolveGravity )THEN

          DO iLevel = 0, nLevels-1

            CALL amrex_multifab_build &
                   ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
                     nDOFX * nGS, swXX )
            CALL MF_uGS(iLevel) % SetVal( Zero ) ! remove this after debugging

            CALL amrex_multifab_build &
                   ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
                     nDOFX * nMF, swXX )
            CALL MF_uMF(iLevel) % SetVal( Zero ) ! remove this after debugging

          END DO

          CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCF, +1 )

          CALL ComputeConformalFactor &
                 ( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

          CALL ComputeLapseShiftCurvature &
                 ( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

          CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCF, -1 )

          DO iLevel = 0, nLevels-1

            CALL amrex_multifab_destroy( MF_uMF(iLevel) )
            CALL amrex_multifab_destroy( MF_uGS(iLevel) )

          END DO

        END IF ! EvolveGravity

        ! --- nLevels <= nMaxLevels; entire arrays t_old(0:nMaxLevels-1) and
        !     t_new(0:nMaxLevels-1) must have valid data ---
        t_old = t_old(0)
        t_new = t_new(0)

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,*)
            WRITE(*,'(A)') 'CALL ComputeFromConserved_Euler_MF'
            WRITE(*,*)

          END IF

          CALL ComputeFromConserved_Euler_MF &
                 ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

        END IF

      END IF ! MOD( StepNo(0), 10 ) .EQ. 0

    END IF ! UseAMR

  END SUBROUTINE ReGrid


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE ComputeConformalFactor( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

    CALL ComputeConformalFactor_XCFC_MF &
           ( MF_uGS, MF_uMF )

    CALL UpdateConformalFactorAndMetric_XCFC_MF &
           ( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

  END SUBROUTINE ComputeConformalFactor


  SUBROUTINE ComputeLapseShiftCurvature( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

    CALL ComputePressureTensorTrace_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

    CALL ComputeLapseShiftCurvature_XCFC_MF &
           ( MF_uGS, MF_uMF )

    CALL UpdateLapseShiftCurvature_XCFC_MF &
           ( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

    CALL ApplyBoundaryConditions_Geometry_XCFC_MF( MF_uGF )

  END SUBROUTINE ComputeLapseShiftCurvature


END MODULE ReGridModule
