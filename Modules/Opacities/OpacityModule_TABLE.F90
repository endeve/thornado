MODULE OpacityModule_TABLE

#ifdef MICROPHYSICS_WEAKLIB

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

#endif

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV

  IMPLICIT NONE
  PRIVATE

  CHARACTER(256) :: &
    OpacityTableName
#ifdef MICROPHYSICS_WEAKLIB
  TYPE(OpacityTableType) :: &
    OPACITIES
#endif

  PUBLIC :: InitializeOpacities_TABLE
  PUBLIC :: FinalizeOpacities_TABLE
  PUBLIC :: ComputeAbsorptionCoefficients_TABLE

CONTAINS


  SUBROUTINE InitializeOpacities_TABLE( OpacityTableName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Option

    OpacityTableName = 'OpacityTable.h5'
    IF( PRESENT( OpacityTableName_Option ) ) &
      OpacityTableName = TRIM( OpacityTableName_Option )

    WRITE(*,*)
    WRITE(*,'(A7,A12,A)') &
      '', 'Table Name: ', TRIM( OpacityTableName )

#ifdef MICROPHYSICS_WEAKLIB

    CALL InitializeHDF( )

    CALL ReadOpacityTableHDF &
           ( OPACITIES, TRIM( OpacityTableName ) )

    CALL FinalizeHDF( )

#endif

  END SUBROUTINE InitializeOpacities_TABLE


  SUBROUTINE FinalizeOpacities_TABLE

  END SUBROUTINE FinalizeOpacities_TABLE


  SUBROUTINE ComputeAbsorptionCoefficients_TABLE &
               ( E, D, T, Y, Chi, dChidT_Option, dChidY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: E, D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: Chi
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dChidT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dChidY_Option

    LOGICAL                  :: ComputeDerivatives
    INTEGER                  :: iE
    REAL(DP), DIMENSION(1)   :: TMP
    REAL(DP), DIMENSION(1,4) :: dTMP

    ComputeDerivatives = .FALSE.
    IF( ALL( [ PRESENT( dChidT_Option ), PRESENT( dChidY_Option ) ] ) ) &
      ComputeDerivatives = .TRUE.

#ifdef MICROPHYSICS_WEAKLIB

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
      ( Chi_T => OPACITIES % ThermEmAb % Absorptivity(1) % Values )

    IF( ComputeDerivatives )THEN

      DO iE = 1, SIZE( E )

        CALL LogInterpolateDifferentiateSingleVariable &
               ( [ E(iE) ] / MeV, [ D ] / ( Gram / Centimeter**3 ), &
                 [ T ] / Kelvin, [ Y ], E_T, D_T, T_T, Y_T, &
                 [ 1, 1, 1, 0 ], 1.0d-100, Chi_T, TMP, dTMP, debug = .FALSE. )

        Chi(iE) &
          = TMP(1) * ( 1.0_DP / Centimeter )

        dChidT_Option(iE) &
          = dTMP(1,3) * ( 1.0_DP / Centimeter ) / Kelvin

        dChidY_Option(iE) &
          = dTMP(1,4) * ( 1.0_DP / Centimeter )

      END DO

    ELSE

      DO iE = 1, SIZE( E )

        CALL LogInterpolateSingleVariable &
               ( [ E(iE) ] / MeV, [ D ] / ( Gram / Centimeter**3 ), &
                 [ T ] / Kelvin, [ Y ], E_T, D_T, T_T, Y_T, &
                 [ 1, 1, 1, 0 ], 1.0d-100, Chi_T, TMP )

        Chi(iE) = TMP(1) * ( 1.0_DP / Centimeter )

      END DO

    END IF

    END ASSOCIATE ! Chi_T

    END ASSOCIATE ! E_T, etc.

    END ASSOCIATE ! iD_T, etc.

    END ASSOCIATE ! EOS

#else

   Chi(:) = 0.0_DP

   IF( ComputeDerivatives )THEN

     dChidT_Option(:) = 0.0_DP
     dChidY_Option(:) = 0.0_DP

   END IF

#endif

  END SUBROUTINE ComputeAbsorptionCoefficients_TABLE


END MODULE OpacityModule_TABLE
