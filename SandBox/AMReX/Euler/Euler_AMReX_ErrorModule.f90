MODULE Euler_AMReX_ErrorModule

  USE amrex_error_module, ONLY: &
    amrex_abort

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DescribeError_Euler_AMReX

CONTAINS

  SUBROUTINE DescribeError_Euler_AMReX( iErr, Message_Option )

    INTEGER,          INTENT(in)           :: iErr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Message_Option

    CHARACTER(LEN=128) :: Message

    Message = ''
    IF( PRESENT( Message_Option ) ) &
      Message = Message_Option

    SELECT CASE( iErr )

      CASE( 00 )

        RETURN

      CASE( 01 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_PositivityLimiterModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') 'U_K(iCF_E) < 0'

          CALL amrex_abort( TRIM( Message ) )

      CASE( 02 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_PositivityLimiterModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') 'SolveTheta_Bisection: No root in interval'

          CALL amrex_abort( TRIM( Message ) )

      CASE( 03 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_PositivityLimiterModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') 'SolveTheta_Bisection: Failure to converge'

          CALL amrex_abort( TRIM( Message ) )

      CASE( 04 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_PositivityLimiterModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') 'q < 0 after all limiting'

          CALL amrex_abort( TRIM( Message ) )

      CASE( 05 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: Euler_ApplyBC_X1'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Fluid X1'

          CALL amrex_abort( TRIM( Message ) )

      CASE( 06 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: Euler_ApplyBC_X2'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Fluid X2'

          CALL amrex_abort( TRIM( Message ) )

      CASE( 07 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: Euler_ApplyBC_X3'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Fluid X3'

          CALL amrex_abort( TRIM( Message ) )

      CASE DEFAULT

          WRITE(*,'(2x,A,I2.2)') 'Unknown error: ', iErr
          WRITE(*,'(2x,A)') 'Stopping...'

          CALL amrex_abort( TRIM( Message ) )

    END SELECT

  END SUBROUTINE DescribeError_Euler_AMReX

END MODULE Euler_AMReX_ErrorModule
