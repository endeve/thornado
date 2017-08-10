MODULE TimeSteppingModule_Penalization

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EvolveFields

CONTAINS


  SUBROUTINE EvolveFields &
               ( t_begin, t_end, dt_write, dt_fixed_Option )

    REAL(DP), INTENT(in) :: t_begin, t_end, dt_write
    REAL(DP), INTENT(in), OPTIONAL :: dt_fixed_Option

    LOGICAL  :: WriteOutput = .FALSE.
    LOGICAL  :: FixedTimeStep
    INTEGER  :: iCycle
    REAL(DP) :: t, t_write, dt, dt_fixed
   
    iCycle  = 0
    t       = t_begin
    t_write = t_begin + dt_write


    DO WHILE( t < t_end )

      iCycle = iCycle + 1


      IF( FixedTimeStep )THEN
        dt = dt_fixed
      ELSE
!        CALL ComputeTimeStep( dt )
      END IF

      IF( t + dt > t_end )THEN

        dt = t_end - t

      END IF

      IF( t + dt > t_write )THEN

        dt          = t_write - t
        t_write     = t_write + dt_write
        WriteOutput = .TRUE.

      END IF

      IF( MOD( iCycle, 100 ) == 0 )THEN

!        WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A2,A2,A5,ES12.6E2,A1,A2)') &
!          '', 'Cycle = ', iCycle, &
!          '', 't = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
!          '', 'dt = ', dt / U % TimeUnit, '', TRIM( U % TimeLabel )
!
!
      END IF
 
!      CALL UpdateFields( t, dt )

      t = t + dt

!      IF( WriteOutput )THEN
!
!        CALL WriteFields1D &
!               ( Time = t, &
!                 WriteGeometryFields_Option = EvolveGravity, &
!                 WriteFluidFields_Option = EvolveFluid, &
!                 WriteRadiationFields_Option = EvolveRadiation )
!
!        CALL WriteFieldsRestart1D &
!
!               ( Time = t, &
!                 WriteGeometryFields_Option = EvolveGravity, &
!                 WriteFluidFields_Option = EvolveFluid, &
!                 WriteRadiationFields_Option = EvolveRadiation )
!
!        WriteOutput = .FALSE.
!
!      END IF

    END DO ! WHILE
   


  END SUBROUTINE EvolveFields


END MODULE TimeSteppingModule_Penalization
