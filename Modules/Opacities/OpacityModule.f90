MODULE OpacityModule

  ! --- weaklib modules --------------------------

  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateDifferentiateSingleVariable

  ! ----------------------------------------------

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV

  IMPLICIT NONE
  PRIVATE

  CHARACTER(5) :: &
    Opacity &
      = 'IDEAL'
  CHARACTER(256) :: &
    OpacityTableName &
      = 'OpacityTable.h5'
  TYPE(OpacityTableType) :: &
    OPACITIES

  PROCEDURE (ComputeOpacity_A), POINTER, PUBLIC :: &
    ComputeAbsorptionCoefficients => NULL()

  INTERFACE
    SUBROUTINE ComputeOpacity_A( E, D, T, Y, Opacity )
      USE KindModule, ONLY: DP
      REAL(DP), DIMENSION(:), INTENT(in)  :: E, D, T, Y
      REAL(DP), DIMENSION(:), INTENT(out) :: Opacity
    END SUBROUTINE ComputeOpacity_A
  END INTERFACE

  PUBLIC :: InitializeOpacities
  PUBLIC :: FinalizeOpacities

CONTAINS


  SUBROUTINE InitializeOpacities &
               ( Opacity_Option, OpacityTableName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: Opacity_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Option

    IF( PRESENT( Opacity_Option ) )THEN
      Opacity = TRIM( Opacity_Option )
    END IF

    IF( PRESENT( OpacityTableName_Option ) )THEN
      OpacityTableName = TRIM( OpacityTableName_Option )
    END IF

    WRITE(*,*)
    WRITE(*,'(A5,A9,A)') &
      '', 'Opacity: ', TRIM( Opacity )
    WRITE(*,'(A5,A9)') &
      '', '-------- '

    SELECT CASE ( TRIM( Opacity ) )
      CASE( 'IDEAL' )

        ComputeAbsorptionCoefficients &
          => ComputeAbsorptionCoefficients_IDEAL

      CASE( 'TABLE' )

        WRITE(*,*)
        WRITE(*,'(A7,A12,A)') &
          '', 'Table Name: ', TRIM( OpacityTableName )

        CALL InitializeHDF( )

        CALL ReadOpacityTableHDF &
               ( OPACITIES, TRIM( OpacityTableName ) )

        CALL FinalizeHDF( )

        ComputeAbsorptionCoefficients &
          => ComputeAbsorptionCoefficients_TABLE

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A17,A5)') &
          '', 'Invalid Opacity: ', TRIM( Opacity )
        STOP

    END SELECT

  END SUBROUTINE InitializeOpacities


  SUBROUTINE FinalizeOpacities

    NULLIFY( ComputeAbsorptionCoefficients )

  END SUBROUTINE FinalizeOpacities


  ! --- Ideal Opacities ---


  SUBROUTINE ComputeAbsorptionCoefficients_IDEAL( E, D, T, Y, Chi )

    REAL(DP), DIMENSION(:), INTENT(in)  :: E, D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Chi

    Chi = 0.0_DP

  END SUBROUTINE ComputeAbsorptionCoefficients_IDEAL


  ! --- Tabulated Opacities (through weaklib) ---


  SUBROUTINE ComputeAbsorptionCoefficients_TABLE( E, D, T, Y, Chi )

    REAL(DP), DIMENSION(:), INTENT(in)  :: E, D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out) :: Chi

    INTEGER                :: iE
    REAL(DP), DIMENSION(1) :: TMP

    ASSOCIATE &
      ( EOS => OPACITIES % EOSTable )

    ASSOCIATE &
      ( iD_T => EOS % TS % Indices % iRho, &
        iT_T => EOS % TS % Indices % iT,   &
        iY_T => EOS % TS % Indices % iYe )

    ASSOCIATE &
      ( E_T => OPACITIES % EnergyGrid  % Values, &
        D_T => EOS % TS % States(iD_T) % Values, &
        T_T => EOS % TS % States(iT_T) % Values, &
        Y_T => EOS % TS % States(iY_T) % Values )

    ASSOCIATE &
      ( Chi_T => OPACITIES % ThermEmAb % Absorptivity(1) % Values, &
        OS    => OPACITIES % ThermEmAb % Offset )

    DO iE = 1, SIZE( E )

!      PRINT*
!      PRINT*, "Interpolating Opacity"
!      PRINT*, "  E = ", [ E(iE) ] / MeV
!      PRINT*, "  D = ", [ D ] / ( Gram / Centimeter**3 )
!      PRINT*, "  T = ", [ T ] / Kelvin
!      PRINT*, "  Y = ", [ Y ]

!      PRINT*
!      PRINT*, "  E_T = ", E_T
!      PRINT*, "  D_T = ", D_T
!      PRINT*, "  T_T = ", T_T
!      PRINT*, "  Y_T = ", Y_T

      CALL LogInterpolateSingleVariable &
             ( [ E(iE) ] / MeV, [ D ] / ( Gram / Centimeter**3 ), &
               [ T ] / Kelvin, [ Y ], E_T, D_T, T_T, Y_T, [ 1, 1, 1, 0 ], &
               OS, Chi_T, TMP )

      Chi(iE) = TMP(1) * ( 1.0_DP / Centimeter )

    END DO

    END ASSOCIATE ! Chi_T

    END ASSOCIATE ! E_T, etc.

    END ASSOCIATE ! iD_T, etc.

    END ASSOCIATE ! EOS

  END SUBROUTINE ComputeAbsorptionCoefficients_TABLE


END MODULE OpacityModule
