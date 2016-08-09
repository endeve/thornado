MODULE TimeSteppingModule_Explicit

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE InputOutputModule, ONLY: &
    WriteFields1D
  USE FluidRadiationCouplingSolutionModule_Explicit, ONLY: &
    ComputeRHS_FluidRadiationCoupling_Explicit_EmissionAbsorption

  IMPLICIT NONE
  PRIVATE

  INTERFACE
    SUBROUTINE TimeStepper( t, dt )
      USE KindModule, ONLY: DP
      REAL(DP), INTENT(in) :: t, dt
    END SUBROUTINE TimeStepper
  END INTERFACE

  PUBLIC :: EvolveFields
  PUBLIC :: ForwardEuler

CONTAINS


  SUBROUTINE EvolveFields( t_begin, t_end, dt_write, UpdateFields )

    REAL(DP),    INTENT(in) :: t_begin, t_end, dt_write
    PROCEDURE (TimeStepper) :: UpdateFields

    INTEGER  :: iCycle
    REAL(DP) :: t, t_write, dt
    REAL(DP), DIMENSION(0:1) :: WallTime

    ASSOCIATE( U => UnitsDisplay )

    WRITE(*,*)
    WRITE(*,'(A4,A21)') '', 'INFO: Evolving Fields'
    WRITE(*,*)
    WRITE(*,'(A6,A10,ES10.4E2,A1,2A2,A8,ES10.4E2,&
              &A1,2A2,A11,ES10.4E2,A1,A2)') &
      '', 't_begin = ',  t_begin  / U % TimeUnit, '', TRIM( U % TimeLabel ), &
      '', 't_end = ',    t_end    / U % TimeUnit, '', TRIM( U % TimeLabel ), &
      '', 'dt_write = ', dt_write / U % TimeUnit, '', TRIM( U % TimeLabel )
    WRITE(*,*)

    iCycle  = 0
    t       = t_begin
    t_write = dt_write

    CALL WriteFields1D( Time = t )

    CALL CPU_TIME( WallTime(0) )

    dt = 0.0_DP

    CALL UpdateFields( t, dt )

    CALL CPU_TIME( WallTime(1) )

    WRITE(*,*)
    WRITE(*,'(A6,A15,ES10.4E2,A1,A2,A6,I8.8,A7,A4,ES10.4E2,A2)') &
      '', 'Evolved to t = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
      ' with ', iCycle, ' cycles', &
      ' in ', WallTime(1)-WallTime(0), ' s'
    WRITE(*,*)

    CALL WriteFields1D( Time = t )

    END ASSOCIATE ! U

  END SUBROUTINE EvolveFields


  SUBROUTINE ForwardEuler( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL ComputeRHS_FluidRadiationCoupling_Explicit_EmissionAbsorption

  END SUBROUTINE ForwardEuler


END MODULE TimeSteppingModule_Explicit
