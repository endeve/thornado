MODULE MF_Euler_ErrorModule

  ! --- thornado Modules ---

  USE UtilitiesModule, ONLY: &
    thornado_abort

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DescribeError_Euler_MF

CONTAINS


  SUBROUTINE DescribeError_Euler_MF &
    ( iErr, Message_Option, Int_Option, Real_Option )

    INTEGER,          INTENT(in)           :: iErr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Message_Option
    INTEGER,          INTENT(in), OPTIONAL :: Int_Option(:)
    REAL(DP),         INTENT(in), OPTIONAL :: Real_Option(:)

    CHARACTER(LEN=128) :: Message

    Message = ''
    IF( PRESENT( Message_Option ) ) &
      Message = TRIM( Message_Option )

    SELECT CASE( iErr )

      CASE( 00 )

        RETURN

      CASE( 01 )

        WRITE(*,*)
        WRITE(*,'(2x,A,A)') 'FATAL ERROR'
        WRITE(*,'(2x,A)') &
          'MODULE: MF_InitializationModule_Relativistic_IDEAL'
        WRITE(*,'(2x,A)') &
          'SUBROUTINE: InitializeFields_Euler_MF'
        WRITE(*,*)
        WRITE(*,'(2x,A)') TRIM( Message )
        WRITE(*,*)
        WRITE(*,'(2x,A)')   'Valid choices'
        WRITE(*,'(2x,A)')   '-------------'
        WRITE(*,'(2x,A)')   '  Advection1D'
        WRITE(*,'(2x,A)')   '  RiemannProblem1D'
        WRITE(*,'(2x,A)')   '  Advection2D'
        WRITE(*,'(2x,A)')   '  KelvinHelmholtz2D'
        WRITE(*,'(2x,A)')   '  Advection3D'
        WRITE(*,*)
        WRITE(*,'(A)') 'Stopping...'

        CALL thornado_abort

      CASE( 02 )

        WRITE(*,'(2x,A,A)') 'FATAL ERROR'
        WRITE(*,'(2x,A)') &
          'MODULE: InitializationModule'
        WRITE(*,'(2x,A)') &
          'SUBROUTINE: InitializeProgram'
        WRITE(*,*)
        WRITE(*,'(2x,A,1x,I3.2)') TRIM( Message ), Int_Option(1)
        WRITE(*,'(2x,A)') 'Valid choices'
        WRITE(*,'(2x,A)') '-------------'
        WRITE(*,'(2x,A)') '  0 (CARTESIAN)'
        WRITE(*,'(2x,A)') '  1 (CYLINDRICAL)'
        WRITE(*,'(2x,A)') '  2 (SPHERICAL)'
        WRITE(*,'(2x,A)')

        CALL thornado_abort

      CASE DEFAULT

          WRITE(*,'(2x,A,I2.2)') 'Unknown error: ', iErr
          WRITE(*,'(2x,A)') TRIM( Message )
          WRITE(*,*)

          CALL thornado_abort

    END SELECT

  END SUBROUTINE DescribeError_Euler_MF


END MODULE MF_Euler_ErrorModule
