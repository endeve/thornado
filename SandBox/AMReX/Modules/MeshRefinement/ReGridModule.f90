MODULE ReGridModule

  ! --- AMReX Modules ---

  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor
  USE amrex_amrcore_module, ONLY: &
   amrex_regrid, &
   amrex_get_numlevels

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
  USE InputParsingModule, ONLY: &
    DEBUG, &
    UseAMR, &
    StepNo, &
    iReGrid, &
    nLevels, &
    nMaxLevels, &
    t_old, &
    t_new

  PUBLIC :: ReGrid

CONTAINS


  SUBROUTINE ReGrid

    INTEGER :: iErr, iLevel

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

        DO iLevel = 0, nLevels

          ! amrex_regrid operates on level `iLevel+1`, so don't operate
          ! on nMaxLevels
          IF( iLevel .LT. nMaxLevels-1 ) &
            CALL amrex_regrid( iLevel, t_new(iLevel) )

        END DO

        nLevels = amrex_get_numlevels()

        ! --- nLevels <= nMaxLevels; entire arrays t_old(0:nMaxLevels-1) and
        !     t_new(0:nMaxLevels-1) must have valid data ---
        t_old = t_old(0)
        t_new = t_new(0)

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,'(6x,A,I2.2)') 'nLevels (after regrid): ', nLevels
            WRITE(*,*)
            WRITE(*,'(A)') 'CALL ApplyBoundaryConditions_Geometry_MF'

          END IF

        END IF

        CALL ApplyBoundaryConditions_Geometry_MF( MF_uGF )

        IF( DEBUG )THEN

          CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

          IF( amrex_parallel_ioprocessor() )THEN

            WRITE(*,*)
            WRITE(*,'(A)') 'CALL ApplyPositivityLimiter_Euler_MF'
            WRITE(*,*)

          END IF

        END IF

        ! --- Regridding may cause some cells to be un-physical ---
        CALL ApplyPositivityLimiter_Euler_MF( MF_uGF, MF_uCF, MF_uDF )

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


END MODULE ReGridModule
