MODULE InitializationModule

  USE KindModule, ONLY: &
    DP,  &
    TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFX,       &
    iX_B0,       &
    iX_E0
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iPF_Ne, &
    uCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iCF_Ne, &
    uAF,    &
    iAF_P,  &
    iAF_T,  &
    iAF_Ye, &
    iAF_E
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UnitsModule, ONLY: &
    Gram,       &
    Centimeter, &
    Kilometer,  &
    Second,     &
    SpeedOfLight

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields


CONTAINS


  SUBROUTINE InitializeFields &
               ( AdvectionProfile_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: AdvectionProfile_Option

    CHARACTER(LEN=64) :: AdvectionProfile = 'SineWave'

    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    WRITE(*,*)
    WRITE(*,'(4x,A,A)') 'INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Advection' )

        CALL InitializeFields_Advection &
               ( TRIM( AdvectionProfile ) )

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

    REAL(DP), PARAMETER :: D_0 = 1.0e12_DP * Gram / Centimeter**3
    REAL(DP), PARAMETER :: Amp = 1.0e11_DP * Gram / Centimeter**3
    REAL(DP), PARAMETER :: L   = 1.0e02_DP * Kilometer

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
              = 3.0e4_DP * Kilometer / Second

            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
              = 0.0_DP * Kilometer / Second

            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
              = 0.0_DP * Kilometer / Second

            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) &
              = 0.01_DP * D_0 * SpeedOfLight**2

            uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
              = 0.3_DP

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


END MODULE InitializationModule
