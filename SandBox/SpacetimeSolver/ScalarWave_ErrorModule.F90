MODULE ScalarWave_ErrorModule


  USE UtilitiesModule, ONLY: &
    thornado_abort

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DescribeError_ScalarWave


CONTAINS


  SUBROUTINE DescribeError_ScalarWave( iErr, Message_Option )

    INTEGER,          INTENT(in)           :: iErr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Message_Option

    CHARACTER(LEN=128) :: Message

    Message = ''
    IF( PRESENT( Message_Option ) ) &
      Message = TRIM( Message_Option )

    SELECT CASE( iErr )

      CASE( 00 )

        RETURN

      CASE( 01 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: ScalarWave_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyBC_ScalarWave_X1'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Scalar Field X1'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 02 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: ScalarWave_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyBC_ScalarWave_X2'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Scalar Field X2'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 03 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: ScalarWave_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyBC_ScalarWave_X3'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Scalar Field X3'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE DEFAULT

          WRITE(*,'(2x,A,I2.2)') 'Unknown error: ', iErr
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

    END SELECT

  END SUBROUTINE DescribeError_ScalarWave


END MODULE ScalarWave_ErrorModule
