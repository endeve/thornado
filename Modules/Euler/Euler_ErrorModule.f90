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
    ( iErr, Message_Option, Int_Option, Real_Option, Char_Option )

    INTEGER,          INTENT(in)           :: iErr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Message_Option
    INTEGER         , INTENT(in), OPTIONAL :: Int_Option(:)
    REAL(DP)        , INTENT(in), OPTIONAL :: Real_Option(:)
    CHARACTER(*)    , INTENT(in), OPTIONAL :: Char_Option(:)

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

          CALL WriteOutput( Int_Option, Real_Option, Char_Option )

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

      CASE( 11 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_UtilitiesModule_Relativistic'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ComputePrimitive_Scalar'
           WRITE(*,'(2x,A)') &
            'Reason for failure: CF_D .LT. Min_D'
          WRITE(*,*)

          CALL WriteOutput( Int_Option, Real_Option, Char_Option )

          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 12 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_UtilitiesModule_Relativistic'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: SolveZ_Bisection_Scalar'
           WRITE(*,'(2x,A)') &
            'Reason for failure: ITERATION .EQ. MAX_IT'
          WRITE(*,*)

          CALL WriteOutput( Int_Option, Real_Option, Char_Option )

          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 13 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_UtilitiesModule_Relativistic'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ComputePrimitive_Scalar'
           WRITE(*,'(2x,A)') &
            'Reason for failure: NANS IN'
          WRITE(*,*)

          CALL WriteOutput( Int_Option, Real_Option, Char_Option )

          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE( 14 )

          WRITE(*,*)
          WRITE(*,'(2x,A)') 'FATAL ERROR'
          WRITE(*,'(2x,A)') &
            'MODULE: Euler_UtilitiesModule_Relativistic'
          WRITE(*,'(2x,A)') &
            'SUBROUTINE: ComputePrimitive_Scalar'
           WRITE(*,'(2x,A)') &
            'Reason for failure: NANS OUT'
          WRITE(*,*)

          CALL WriteOutput( Int_Option, Real_Option, Char_Option )

          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

      CASE DEFAULT

          WRITE(*,'(2x,A,I2.2)') 'Unknown error: ', iErr
          WRITE(*,'(2x,A)') TRIM( Message )

          CALL thornado_abort

    END SELECT

  END SUBROUTINE DescribeError_Euler


  SUBROUTINE WriteOutput( IntArray, RealArray, CharArray )

    INTEGER     , INTENT(in) :: IntArray(:)
    REAL(DP)    , INTENT(in) :: RealArray(:)
    CHARACTER(*), INTENT(in) :: CharArray(:)

    WRITE(*,'(2x,A,I8.8)') &
      'ITERATION:           ', IntArray(1)
    IF( IntArray(2) .NE. 99999999 ) &
      WRITE(*,'(2x,A,I9.8)') &
        'iNX (Packed array): ', IntArray(2)
    IF( TRIM( CharArray(1) ) .NE. 'NA' ) &
      WRITE(*,'(2x,A,A)') &
        'iDimX: ', TRIM( CharArray(1) )
    WRITE(*,'(2x,A,4I7.5)') &
      'iX1, iX2, iX3, iNX: ', &
       IntArray(10), IntArray(11), IntArray(12), IntArray(9)
    WRITE(*,'(2x,A,3I7.5)') &
      'iX_B0:              ', IntArray(3), IntArray(4), IntArray(5)
    WRITE(*,'(2x,A,3I7.5)') &
      'iX_E0:              ', IntArray(6), IntArray(7), IntArray(8)
    WRITE(*,'(2x,A,SP3ES15.6E3)') &
      'X1_C, X2_C, X3_C: ', &
       RealArray(1), RealArray(2), RealArray(3)
    WRITE(*,'(2x,A,SP3ES15.6E3)') &
      'dX1, dX2, dX3:    ', &
       RealArray(4), RealArray(5), RealArray(6)
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'U(iCF_D       ) = ', RealArray(7), '_DP'
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'U(iCF_S1      ) = ', RealArray(8), '_DP'
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'U(iCF_S2      ) = ', RealArray(9), '_DP'
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'U(iCF_S3      ) = ', RealArray(10), '_DP'
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'U(iCF_E       ) = ', RealArray(11), '_DP'
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'U(iCF_Ne      ) = ', RealArray(12), '_DP'
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'G(iGF_Gm_dd_11) = ', RealArray(13), '_DP'
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'G(iGF_Gm_dd_22) = ', RealArray(14), '_DP'
    WRITE(*,'(2x,A,SPES24.16E3,A)') &
      'G(iGF_Gm_dd_33) = ', RealArray(15), '_DP'

  END SUBROUTINE WriteOutput


END MODULE Euler_ErrorModule
