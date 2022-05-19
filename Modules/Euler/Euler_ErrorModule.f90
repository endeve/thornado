MODULE Euler_ErrorModule

  USE KindModule, ONLY: &
    DP
  USE UtilitiesModule, ONLY: &
    thornado_abort

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DescribeError_Euler


CONTAINS


  SUBROUTINE DescribeError_Euler &
    ( iErr, Message_Option, Int_Option, Real_Option )

    INTEGER,          INTENT(in)           :: iErr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Message_Option
    INTEGER         , INTENT(in), OPTIONAL :: Int_Option(:)
    REAL(DP)        , INTENT(in), OPTIONAL :: Real_Option(:)

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
            'MODULE: Euler_PositivityLimiterModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') 'U_K(iCF_E) < 0'
          WRITE(*,*)
          WRITE(*,'(2x,A,I8.8)')      'iX1:       ', Int_Option(1)
          WRITE(*,'(2x,A,I8.8)')      'iX2:       ', Int_Option(2)
          WRITE(*,'(2x,A,I8.8)')      'iX3:       ', Int_Option(3)
          WRITE(*,'(2x,A,ES24.16E3)') 'uD  (iNX): ', Real_Option(1)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS1 (iNX): ', Real_Option(2)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS2 (iNX): ', Real_Option(3)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS3 (iNX): ', Real_Option(4)
          WRITE(*,'(2x,A,ES24.16E3)') 'uE  (iNX): ', Real_Option(5)
          WRITE(*,'(2x,A,ES24.16E3)') 'h1  (iNX): ', Real_Option(6)
          WRITE(*,'(2x,A,ES24.16E3)') 'h2  (iNX): ', Real_Option(7)
          WRITE(*,'(2x,A,ES24.16E3)') 'h3  (iNX): ', Real_Option(8)
          WRITE(*,'(2x,A)') TRIM( Message )


          CALL thornado_abort

      CASE( 02 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_PositivityLimiterModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') 'SolveTheta_Bisection: No root in interval'
          WRITE(*,'(2x,A,I8.8)')      'iX1:       ', Int_Option(1)
          WRITE(*,'(2x,A,I8.8)')      'iX2:       ', Int_Option(2)
          WRITE(*,'(2x,A,I8.8)')      'iX3:       ', Int_Option(3)
          WRITE(*,'(2x,A,ES24.16E3)') 'uD  (iNX): ', Real_Option(1)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS1 (iNX): ', Real_Option(2)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS2 (iNX): ', Real_Option(3)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS3 (iNX): ', Real_Option(4)
          WRITE(*,'(2x,A,ES24.16E3)') 'uE  (iNX): ', Real_Option(5)
          WRITE(*,'(2x,A,ES24.16E3)') 'h1  (iNX): ', Real_Option(6)
          WRITE(*,'(2x,A,ES24.16E3)') 'h2  (iNX): ', Real_Option(7)
          WRITE(*,'(2x,A,ES24.16E3)') 'h3  (iNX): ', Real_Option(8)
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 03 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_PositivityLimiterModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') 'SolveTheta_Bisection: Failure to converge'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 04 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_PositivityLimiterModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') 'q < 0 after all limiting'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 05 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyBC_Euler_X1'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Fluid X1'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 06 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyBC_Euler_X2'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Fluid X2'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 07 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_BoundaryConditionsModule'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyBC_Euler_X3'
           WRITE(*,'(2x,A)') &
            'Invalid Boundary Condition for Fluid X3'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 08 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_UtilitiesModule_Relativistic'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: SolveZ_Bisection'
           WRITE(*,'(2x,A)') &
            'No Root in Interval'
          WRITE(*,*)
          WRITE(*,'(2x,A,I8.8)')      'iNX:       ', Int_Option(1)
          WRITE(*,'(2x,A,ES24.16E3)') 'uD  (iNX): ', Real_Option(1)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS1 (iNX): ', Real_Option(2)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS2 (iNX): ', Real_Option(3)
          WRITE(*,'(2x,A,ES24.16E3)') 'uS3 (iNX): ', Real_Option(4)
          WRITE(*,'(2x,A,ES24.16E3)') 'uE  (iNX): ', Real_Option(5)
          WRITE(*,'(2x,A,ES24.16E3)') 'uNe (iNX): ', Real_Option(6)
          WRITE(*,'(2x,A,ES24.16E3)') 'Gm11(iNX): ', Real_Option(7)
          WRITE(*,'(2x,A,ES24.16E3)') 'Gm22(iNX): ', Real_Option(8)
          WRITE(*,'(2x,A,ES24.16E3)') 'Gm33(iNX): ', Real_Option(9)
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 09 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_UtilitiesModule_Relativistic'
          WRITE(*,'(2x,A)') &
            'FUNCTION: AlphaMiddle_Euler_Relativistic'
           WRITE(*,'(2x,A)') &
            'AlphaMiddle Undefined'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 10 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: MF_InitializationModule_Relativistic_IDEAL'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ApplyJumpConditions_LeftState'
           WRITE(*,'(2x,A)') &
            'Root not bracketed'
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE DEFAULT

          WRITE(*,'(2x,A,I2.2)') 'Unknown error: ', iErr
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

    END SELECT

  END SUBROUTINE DescribeError_Euler


END MODULE Euler_ErrorModule
