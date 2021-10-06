MODULE InitializationModule

  USE KindModule, ONLY: &
    DP,  &
    TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX, &
    iX_B0, &
    iX_E0
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    unitsPF, &
    nPF, &
    uPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    uCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    unitsAF, &
    nAF, &
    uAF, &
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_E
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UnitsModule, ONLY: &
    Kilometer, &
    SpeedOfLight, &
    UnitsDisplay

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields


CONTAINS


  SUBROUTINE InitializeFields &
               ( AdvectionProfile_Option, &
                 RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblemName_Option

    CHARACTER(LEN=64) :: AdvectionProfile   = 'SineWave'
    CHARACTER(LEN=64) :: RiemannProblemName = 'Sod'

    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    IF( PRESENT( RiemannProblemName_Option ) ) &
      RiemannProblemName = TRIM( RiemannProblemName_Option )

    WRITE(*,*)
    WRITE(*,'(4x,A,A)') 'INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection' )

        CALL InitializeFields_Advection &
               ( TRIM( AdvectionProfile ) )

      CASE( 'RiemannProblem' )

        CALL InitializeFields_RiemannProblem &
               ( TRIM( RiemannProblemName ) )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_Advection( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    REAL(DP) :: D_0 = 1.0e12_DP
    REAL(DP) :: Amp = 1.0e11_DP
    REAL(DP) :: L   = 1.0e02_DP

    D_0 = D_0 * unitsPF(iPF_D)
    Amp = Amp * unitsPF(iPF_D)
    L   = L   * Kilometer

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          CASE( 'SineWave' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
              = D_0 + Amp * SIN( TwoPi * X1 / L )

            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = 3.0e4_DP * unitsPF(iPF_V1)

            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = 0.0_DP   * unitsPF(iPF_V2)

            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = 0.0_DP   * unitsPF(iPF_V3)

            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) &
              = 0.01_DP  * uPF(iNodeX,iX1,iX2,iX3,iPF_D) * SpeedOfLight**2

            uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
              = 0.3_DP   * unitsAF(iAF_Ye)

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  SineWave'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeTemperatureFromPressure &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), &
               uAF(:,iX1,iX2,iX3,iAF_T ) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), &
               uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), &
               uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), &
               uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), &
               uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_RiemannProblem &
    ( RiemannProblemName )

    CHARACTER(LEN=*), INTENT(in)           :: RiemannProblemName

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1, XD

    REAL(DP) :: LeftStatePF(nPF), RightStatePF(nPF)
    REAL(DP) :: LeftStateAF(nAF), RightStateAF(nAF)

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )
    WRITE(*,*)

    SELECT CASE( TRIM( RiemannProblemName ) )

      CASE( 'Sod' )

        XD = 0.0_DP * Kilometer

        LeftStatePF(iPF_D ) = 1.00e12_DP * unitsPF(iPF_D )
        LeftStatePF(iPF_V1) = 0.0_DP     * unitsPF(iPF_V1)
        LeftStatePF(iPF_V2) = 0.0_DP     * unitsPF(iPF_V2)
        LeftStatePF(iPF_V3) = 0.0_DP     * unitsPF(iPF_V3)
        LeftStateAF(iAF_P ) = 1.00d32    * unitsAF(iAF_P )
        LeftStateAF(iAF_Ye) = 0.4_DP     * unitsAF(iAF_Ye)

        RightStatePF(iPF_D ) = 1.25e11_DP * unitsPF(iPF_D )
        RightStatePF(iPF_V1) = 0.0_DP     * unitsPF(iPF_V1)
        RightStatePF(iPF_V2) = 0.0_DP     * unitsPF(iPF_V2)
        RightStatePF(iPF_V3) = 0.0_DP     * unitsPF(iPF_V3)
        RightStateAF(iAF_P ) = 1.00d31    * unitsAF(iAF_P )
        RightStateAF(iAF_Ye) = 0.3_DP     * unitsAF(iAF_Ye)

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A,A)') &
          'Invalid choice for RiemannProblemName: ', RiemannProblemName
        WRITE(*,'(A)') 'Valid choices:'
        WRITE(*,'(A)') &
          "  'Sod' - &
          Sod's shock tube"
        WRITE(*,'(A)') 'Stopping...'
        STOP

    END SELECT

    WRITE(*,'(6x,A,F8.6)') 'XD = ', XD
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Right State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'PF_D  = ', &
      RightStatePF(iPF_D ) / unitsPF(iPF_D ), UnitsDisplay % MassDensityLabel
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'PF_V1 = ', &
      RightStatePF(iPF_V1) / unitsPF(iPF_V1), UnitsDisplay % VelocityX1Label
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'PF_V2 = ', &
      RightStatePF(iPF_V2) / unitsPF(iPF_V2), UnitsDisplay % VelocityX2Label
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'PF_V3 = ', &
      RightStatePF(iPF_V3) / unitsPF(iPF_V3), UnitsDisplay % VelocityX3Label
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'AF_P  = ', &
      RightStateAF(iAF_P ) / unitsAF(iAF_P ), UnitsDisplay % PressureLabel
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'AF_Ye = ', &
      RightStateAF(iAF_Ye) / unitsAF(iAF_Ye)
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Left State:'
    WRITE(*,*)
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'PF_D  = ', &
      LeftStatePF(iPF_D ) / unitsPF(iPF_D ), UnitsDisplay % MassDensityLabel
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'PF_V1 = ', &
      LeftStatePF(iPF_V1) / unitsPF(iPF_V1), UnitsDisplay % VelocityX1Label
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'PF_V2 = ', &
      LeftStatePF(iPF_V2) / unitsPF(iPF_V2), UnitsDisplay % VelocityX2Label
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'PF_V3 = ', &
      LeftStatePF(iPF_V3) / unitsPF(iPF_V3), UnitsDisplay % VelocityX3Label
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'AF_P  = ', &
      LeftStateAF(iAF_P ) / unitsAF(iAF_P ), UnitsDisplay % PressureLabel
    WRITE(*,'(8x,A,ES14.6E3,1x,A)') 'AF_Ye = ', &
      LeftStateAF(iAF_Ye) / unitsAF(iAF_Ye)

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LE. XD )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = LeftStatePF(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = LeftStatePF(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = LeftStatePF(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = LeftStatePF(iPF_V3)
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = LeftStateAF(iAF_P )
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = LeftStateAF(iAF_Ye)

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = RightStatePF(iPF_D )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = RightStatePF(iPF_V1)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = RightStatePF(iPF_V2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = RightStatePF(iPF_V3)
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = RightStateAF(iAF_P )
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = RightStateAF(iAF_Ye)

        END IF

      END DO

      CALL ComputeTemperatureFromPressure &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), &
               uAF(:,iX1,iX2,iX3,iAF_T ) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne) )

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), &
               uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), &
               uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), &
               uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), &
               uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), &
               uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), &
               uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblem


END MODULE InitializationModule
