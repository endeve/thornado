MODULE Euler_ErrorModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE UtilitiesModule, ONLY: &
    thornado_abort
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem

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

          CALL WriteOutput_PL( Int_Option, Real_Option )

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
          WRITE(*,*)

          CALL WriteOutput_PL( Int_Option, Real_Option )

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
          WRITE(*,*)

          CALL WriteOutput_PL( Int_Option, Real_Option )

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
          WRITE(*,*)

          CALL WriteOutput_PL( Int_Option, Real_Option )

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
            'SUBROUTINE: SolveF_Bisection'
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
            'SUBROUTINE: SolveF_Bisection_Scalar'
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
    WRITE(*,'(2x,A,SPES15.6E3,1x,A)') &
      'X1_C: ', RealArray(1) / UnitsDisplay % LengthX1Unit, &
      TRIM( UnitsDisplay % LengthX1Label )
    WRITE(*,'(2x,A,SPES15.6E3,1x,A)') &
      'X2_C: ', RealArray(2) / UnitsDisplay % LengthX2Unit, &
      TRIM ( UnitsDisplay % LengthX2Label )
    WRITE(*,'(2x,A,SPES15.6E3,1x,A)') &
      'X3_C: ', RealArray(3) / UnitsDisplay % LengthX3Unit, &
      TRIM( UnitsDisplay % LengthX3Label )
    WRITE(*,'(2x,A,SPES15.6E3,1x,A)') &
      'dX1:  ', RealArray(4) / UnitsDisplay % LengthX1Unit, &
      TRIM( UnitsDisplay % LengthX1Label )
    WRITE(*,'(2x,A,SPES15.6E3,1x,A)') &
      'dX2:  ', RealArray(5) / UnitsDisplay % LengthX2Unit, &
      TRIM( UnitsDisplay % LengthX2Label )
    WRITE(*,'(2x,A,SPES15.6E3,1x,A)') &
      'dX3:  ', RealArray(6) / UnitsDisplay % LengthX3Unit, &
      TRIM( UnitsDisplay % LengthX3Label )

    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U(iCF_D       ) = ', &
      RealArray(7)  /       UnitsDisplay % MassDensityUnit, &
      '_DP',          TRIM( UnitsDisplay % MassDensityLabel )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U(iCF_S1      ) = ', &
      RealArray(8)  /       UnitsDisplay % MomentumDensityX1Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX1Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U(iCF_S2      ) = ', &
      RealArray(9)  /       UnitsDisplay % MomentumDensityX2Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX2Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U(iCF_S3      ) = ', &
      RealArray(10) /       UnitsDisplay % MomentumDensityX3Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX3Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U(iCF_E       ) = ', &
      RealArray(11) /       UnitsDisplay % EnergyDensityUnit, &
      '_DP',          TRIM( UnitsDisplay % EnergyDensityLabel )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U(iCF_Ne      ) = ', &
      RealArray(12) /       UnitsDisplay % ParticleDensityUnit, &
      '_DP',          TRIM( UnitsDisplay % ParticleDensityLabel )

    IF( TRIM( CoordinateSystem ) .EQ. 'CARTESIAN' )THEN

      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_11) = ', RealArray(13), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_22) = ', RealArray(14), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_33) = ', RealArray(15), '_DP'

    ELSE IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDRICAL' )THEN

      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_11) = ', RealArray(13), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_22) = ', RealArray(14), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_33) = ', &
        RealArray(15) /       ( UnitsDisplay % LengthX1Unit )**2, &
        '_DP',            TRIM( UnitsDisplay % LengthX1Label ), '^2'

    ELSE IF( TRIM( CoordinateSystem ) .EQ. 'SPHERICAL' )THEN

      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_11) = ', RealArray(13), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_22) = ', &
        RealArray(14) /       ( UnitsDisplay % LengthX1Unit )**2, &
        '_DP',            TRIM( UnitsDisplay % LengthX1Label ), '^2'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G(iGF_Gm_dd_33) = ', &
        RealArray(15) /       ( UnitsDisplay % LengthX1Unit )**2, &
        '_DP',            TRIM( UnitsDisplay % LengthX1Label ), '^2'

    ELSE

      WRITE(*,'(2x,A,A)') &
        'Invalid coordinate system: ', TRIM( CoordinateSystem )

    END IF

  END SUBROUTINE WriteOutput


  SUBROUTINE WriteOutput_PL( IntArray, RealArray )

    INTEGER , INTENT(in) :: IntArray(:)
    REAL(DP), INTENT(in) :: RealArray(:)

    WRITE(*,'(2x,A,4I7.5)') &
      'iX1, iX2, iX3, iPT: ', &
       IntArray(7), IntArray(8), IntArray(9), IntArray(10)
    WRITE(*,'(2x,A,3I7.5)') &
      'iX_B0:              ', IntArray(1), IntArray(2), IntArray(3)
    WRITE(*,'(2x,A,3I7.5)') &
      'iX_E0:              ', IntArray(4), IntArray(5), IntArray(6)

    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_K(iCF_D       ) = ', &
      RealArray(1)  /       UnitsDisplay % MassDensityUnit, &
      '_DP',          TRIM( UnitsDisplay % MassDensityLabel )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_K(iCF_S1      ) = ', &
      RealArray(2)  /       UnitsDisplay % MomentumDensityX1Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX1Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_K(iCF_S2      ) = ', &
      RealArray(3)  /       UnitsDisplay % MomentumDensityX2Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX2Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_K(iCF_S3      ) = ', &
      RealArray(4)  /       UnitsDisplay % MomentumDensityX3Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX3Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_K(iCF_E       ) = ', &
      RealArray(5)  /       UnitsDisplay % EnergyDensityUnit, &
      '_DP',          TRIM( UnitsDisplay % EnergyDensityLabel )

    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'q_K               = ', &
      RealArray(6)  /       UnitsDisplay % EnergyDensityUnit, &
      '_DP',          TRIM( UnitsDisplay % EnergyDensityLabel )

    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_P(iCF_D       ) = ', &
      RealArray(7)  /       UnitsDisplay % MassDensityUnit, &
      '_DP',          TRIM( UnitsDisplay % MassDensityLabel )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_P(iCF_S1      ) = ', &
      RealArray(8)  /       UnitsDisplay % MomentumDensityX1Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX1Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_P(iCF_S2      ) = ', &
      RealArray(9)  /       UnitsDisplay % MomentumDensityX2Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX2Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_P(iCF_S3      ) = ', &
      RealArray(10) /       UnitsDisplay % MomentumDensityX3Unit, &
      '_DP',          TRIM( UnitsDisplay % MomentumDensityX3Label )
    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'U_P(iCF_E       ) = ', &
      RealArray(11) /       UnitsDisplay % EnergyDensityUnit, &
      '_DP',         TRIM( UnitsDisplay % EnergyDensityLabel )

    WRITE(*,'(2x,A,SPES24.16E3,A,1x,A)') &
      'q_P               = ', &
      RealArray(12) /       UnitsDisplay % EnergyDensityUnit, &
      '_DP',          TRIM( UnitsDisplay % EnergyDensityLabel )

    IF( TRIM( CoordinateSystem ) .EQ. 'CARTESIAN' )THEN

      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_11) = ', RealArray(13), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_22) = ', RealArray(14), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_33) = ', RealArray(15), '_DP'

    ELSE IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDRICAL' )THEN

      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_11) = ', RealArray(13), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_22) = ', RealArray(14), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_33) = ', &
        RealArray(15) /     ( UnitsDisplay % LengthX1Unit  )**2, &
        '_DP',          TRIM( UnitsDisplay % LengthX1Label ), '^2'

    ELSE IF( TRIM( CoordinateSystem ) .EQ. 'SPHERICAL' )THEN

      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_11) = ', RealArray(13), '_DP'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_22) = ', &
        RealArray(14) /     ( UnitsDisplay % LengthX1Unit  )**2, &
        '_DP',          TRIM( UnitsDisplay % LengthX1Label ), '^2'
      WRITE(*,'(2x,A,SPES24.16E3,A,1x,A,A)') &
        'G_P(iGF_Gm_dd_33) = ', &
        RealArray(15) /     ( UnitsDisplay % LengthX1Unit  )**2, &
        '_DP',          TRIM( UnitsDisplay % LengthX1Label ), '^2'

    ELSE

      WRITE(*,'(2x,A,A)') &
        'Invalid coordinate system: ', TRIM( CoordinateSystem )

    END IF

  END SUBROUTINE WriteOutput_PL


END MODULE Euler_ErrorModule
