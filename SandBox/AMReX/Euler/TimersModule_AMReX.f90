MODULE TimersModule_AMReX

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
    I8=>INT64

  USE amrex_fort_module, ONLY: &
    amrex_real
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  IMPLICIT NONE
  PRIVATE

  LOGICAL,  PUBLIC :: TimeIt_AMReX

  REAL(amrex_real), PUBLIC :: Timer_AMReX_Program
  REAL(amrex_real), PUBLIC :: Timer_AMReX_InputOutput
  REAL(amrex_real), PUBLIC :: Timer_AMReX_ComputeInc
  REAL(amrex_real), PUBLIC :: Timer_AMReX_Initialize   
  REAL(amrex_real), PUBLIC :: Timer_AMReX_DataTransfer
  REAL(amrex_real), PUBLIC :: Timer_AMReX_InternalBC
  REAL(amrex_real), PUBLIC :: Timer_AMReX_CopyMF

  CHARACTER(17) :: OutputFMT = '(7x,A,F6.3,A)'

  PUBLIC :: InitializeTimers_AMReX
  PUBLIC :: FinalizeTimers_AMReX
  PUBLIC :: TimersStart_AMReX
  PUBLIC :: TimersStop_AMReX
  PUBLIC :: TimersWtime_AMReX

CONTAINS


  SUBROUTINE InitializeTimers_AMReX

    IF( .NOT. TimeIt_AMReX ) RETURN

    Timer_AMReX_InputOutput  = 0.0_amrex_real
    Timer_AMReX_Initialize   = 0.0_amrex_real
    Timer_AMReX_DataTransfer = 0.0_amrex_real
    Timer_AMReX_ComputeInc   = 0.0_amrex_real
    Timer_AMReX_InternalBC   = 0.0_amrex_real
    Timer_AMReX_CopyMF       = 0.0_amrex_real
    Timer_AMReX_Program      = 0.0_amrex_real

    CALL TimersStart_AMReX( Timer_AMReX_Program )

    RETURN
  END SUBROUTINE InitializeTimers_AMReX


  SUBROUTINE FinalizeTimers_AMReX

    IF( .NOT. TimeIt_AMReX ) RETURN

    CALL TimersStop_AMReX( Timer_AMReX_Program )

    IF( amrex_parallel_ioprocessor() )THEN
      WRITE(*,*)
      WRITE(*,'(5x,A)') 'Timers_AMReX Summary'
      WRITE(*,'(5x,A)') '--------------------'
      WRITE(*,*)

      WRITE(*,OutputFMT) &
        'Initialize:       ', &
        100.0_amrex_real * Timer_AMReX_Initialize / Timer_AMReX_Program, ' %'
      WRITE(*,OutputFMT) &
        'ComputeIncrement: ', &
        100.0_amrex_real * Timer_AMReX_ComputeInc / Timer_AMReX_Program, ' %'
      WRITE(*,OutputFMT) &
        'DataTransfer:     ', &
        100.0_amrex_real * Timer_AMReX_DataTransfer / Timer_AMReX_Program, ' %'
      WRITE(*,OutputFMT) &
        'InternalBC:       ', &
        100.0_amrex_real * Timer_AMReX_InternalBC / Timer_AMReX_Program, ' %'
      WRITE(*,OutputFMT) &
        'CopyMF:           ', &
        100.0_amrex_real * Timer_AMReX_CopyMF / Timer_AMReX_Program, ' %'
      WRITE(*,OutputFMT) &
        'InputOutput:      ', &
        100.0_amrex_real * Timer_AMReX_InputOutput / Timer_AMReX_Program, ' %'
      WRITE(*,*)
      WRITE(*,'(7x,A,ES13.6E3,A)') &
        'Program:          ', &
          Timer_AMReX_Program, ' s'
      WRITE(*,*)
    END IF

    RETURN
  END SUBROUTINE FinalizeTimers_AMReX


  SUBROUTINE TimersStart_AMReX( Timer )

    REAL(amrex_real), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX ) RETURN

    Timer = Timer - TimersWtime_AMReX()

    RETURN
  END SUBROUTINE TimersStart_AMReX


  SUBROUTINE TimersStop_AMReX( Timer )

    REAL(amrex_real), INTENT(inout) :: Timer

    IF( .NOT. TimeIt_AMReX ) RETURN

    Timer = Timer + TimersWtime_AMReX()

    RETURN
  END SUBROUTINE TimersStop_AMReX


  REAL(amrex_real) FUNCTION TimersWtime_AMReX()

    INTEGER(I8) :: clock_read
    INTEGER(I8) :: clock_rate
    INTEGER(I8) :: clock_max

    IF( .NOT. TimeIt_AMReX ) RETURN

    CALL SYSTEM_CLOCK( clock_read, clock_rate, clock_max )
    TimersWtime_AMReX = REAL( clock_read, amrex_real ) &
                          / REAL( clock_rate, amrex_real )

    RETURN
  END FUNCTION TimersWtime_AMReX


END MODULE TimersModule_AMReX
