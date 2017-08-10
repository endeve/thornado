MODULE TimeSteppingModule_Penalization

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay, &
    Millisecond
  USE InputOutputModule, ONLY: &
    WriteFields1D
  USE FluidRadiationCouplingSolutionModule_Penalization, ONLY: &
    ComputeRHS_Penalization


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EvolveFields

CONTAINS


  SUBROUTINE EvolveFields &
               ( t_begin, t_end, dt_write, dt_fixed_Option, &
                 iDisplayCycle_Option )

    REAL(DP), INTENT(in) :: t_begin, t_end, dt_write
    REAL(DP), INTENT(in), OPTIONAL :: dt_fixed_Option
    INTEGER,  INTENT(in), OPTIONAL :: iDisplayCycle_Option

    LOGICAL  :: WriteOutput = .FALSE.
    LOGICAL  :: FixedTimeStep
    INTEGER  :: iCycle, iDisplayCycle
    REAL(DP) :: t, t_write, dt, dt_fixed

    FixedTimeStep = .FALSE.
    IF( PRESENT( dt_fixed_Option ) )THEN
      FixedTimeStep = .TRUE.
      dt_fixed = dt_fixed_Option
    END IF

    iDisplayCycle = 10
    IF( PRESENT( iDisplayCycle_Option ) ) &
      iDisplayCycle = iDisplayCycle_Option

    ASSOCIATE( U => UnitsDisplay )

    iCycle  = 0
    t       = t_begin
    t_write = t_begin + dt_write

    CALL WriteFields1D &
           ( Time = t, WriteFluidFields_Option = .TRUE., &
             WriteRadiationFields_Option = .TRUE. )

    DO WHILE( t < t_end )

      iCycle = iCycle + 1


      IF( FixedTimeStep )THEN
        dt = dt_fixed
      ELSE
        CALL ComputeTimeStep( dt )
      END IF

      IF( t + dt > t_end )THEN

        dt = t_end - t

      END IF

      IF( t + dt > t_write )THEN

        dt          = t_write - t
        t_write     = t_write + dt_write
        WriteOutput = .TRUE.

      END IF

      IF( MOD( iCycle, iDisplayCycle ) == 0 )THEN

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A2,A2,A5,ES12.6E2,A1,A2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
          '', 'dt = ', dt / U % TimeUnit, '', TRIM( U % TimeLabel )

      END IF
 
      CALL UpdateFields( t, dt )

      t = t + dt

      IF( WriteOutput )THEN

        CALL WriteFields1D &
               ( Time = t, WriteFluidFields_Option = .TRUE., &
                 WriteRadiationFields_Option = .TRUE. )

        WriteOutput = .FALSE.

      END IF

    END DO ! WHILE

    WRITE(*,*)
    WRITE(*,'(A6,A15,ES10.4E2,A1,A2,A6,I8.8,A7,A4,ES10.4E2,A2)') &
      '', 'Evolved to t = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
      ' with ', iCycle, ' cycles'
    WRITE(*,*)

    END ASSOCIATE ! U

    CALL WriteFields1D &
           ( Time = t, WriteFluidFields_Option = .TRUE., &
             WriteRadiationFields_Option = .TRUE. )

  END SUBROUTINE EvolveFields


  SUBROUTINE ComputeTimestep( dt )

    REAL(DP), INTENT(out) :: dt

    dt = 1.0d-2 * Millisecond

  END SUBROUTINE ComputeTimestep


  SUBROUTINE UpdateFields( t, dt )

    REAL(DP), INTENT(in) :: t, dt

!    CALL Initialized_Penalization

!    CALL ApplyBoundaryConditions_Fluid

!    CALL ComputeRHS_Penalization

!   IF( EvolveRadiation )THEN

!      CALL ApplyBoundaryConditions_Radiation( t )

!      CALL ComputeRHS_Radiation &
!             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

!    END IF

!    IF( EvolveFluid )THEN
!
!      CALL ApplyRHS_Fluid &
!             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
!               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )
!
!      CALL ApplyBoundaryConditions_Fluid
!
!      CALL ApplySlopeLimiter_Fluid
!
!      CALL ApplyPositivityLimiter_Fluid
!
!    END IF
!
!    IF( EvolveRadiation )THEN
!
!      CALL ApplyRHS_Radiation &
!             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
!               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )
!
!      CALL ApplyBoundaryConditions_Radiation &
!             ( t + dt, LimiterBC_Option = .TRUE. )
!
!      CALL ApplySlopeLimiter_Radiation
!
!      CALL ApplyPositivityLimiter_Radiation
!
!    END IF
!
!    CALL Finalize_SSP_RK
!

  END SUBROUTINE UpdateFields


END MODULE TimeSteppingModule_Penalization
