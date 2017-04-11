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
  INTEGER :: &
    iD_T, iT_T, iY_T
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Es_T, Ds_T, Ts_T, Ys_T
#ifdef MICROPHYSICS_WEAKLIB
  TYPE(OpacityTableType) :: &
    OPACITIES
#endif

  PUBLIC :: InitializeOpacities_TABLE
  PUBLIC :: FinalizeOpacities_TABLE
  PUBLIC :: ComputeAbsorptionOpacity_TABLE
  PUBLIC :: ComputeScatteringOpacity_ES_TABLE

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

    ! --- Thermodynamic State Indices ---

    iD_T = OPACITIES % TS % Indices % iRho
    iT_T = OPACITIES % TS % Indices % iT
    iY_T = OPACITIES % TS % Indices % iYe

    ! --- Thermodynamic States ---

    ALLOCATE( Ds_T(OPACITIES % TS % nPoints(iD_T)) )
    Ds_T = OPACITIES % TS % States(iD_T) % Values

    ALLOCATE( Ts_T(OPACITIES % TS % nPoints(iT_T)) )
    Ts_T = OPACITIES % TS % States(iT_T) % Values

    ALLOCATE( Ys_T(OPACITIES % TS % nPoints(iY_T)) )
    Ys_T = OPACITIES % TS % States(iY_T) % Values

    ! --- Energy Grid ---

    ALLOCATE( Es_T(OPACITIES % EnergyGrid % nPoints) )
    Es_T = OPACITIES % EnergyGrid  % Values

#endif

  END SUBROUTINE InitializeOpacities_TABLE


  SUBROUTINE FinalizeOpacities_TABLE

    DEALLOCATE( Es_T, Ds_T, Ts_T, Ys_T )

  END SUBROUTINE FinalizeOpacities_TABLE


  SUBROUTINE ComputeAbsorptionOpacity_TABLE &
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
      ( Chi_T => OPACITIES % ThermEmAb % Absorptivity(1) % Values )

    IF( ComputeDerivatives )THEN

      DO iE = 1, SIZE( E )

        CALL LogInterpolateDifferentiateSingleVariable &
               ( [ E(iE) ] / MeV, [ D ] / ( Gram / Centimeter**3 ), &
                 [ T ] / Kelvin, [ Y ], Es_T, Ds_T, Ts_T, Ys_T, &
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
                 [ T ] / Kelvin, [ Y ], Es_T, Ds_T, Ts_T, Ys_T, &
                 [ 1, 1, 1, 0 ], 1.0d-100, Chi_T, TMP )

        Chi(iE) = TMP(1) * ( 1.0_DP / Centimeter )

      END DO

    END IF

    END ASSOCIATE ! Chi_T

#else

   Chi(:) = 0.0_DP

   IF( ComputeDerivatives )THEN

     dChidT_Option(:) = 0.0_DP
     dChidY_Option(:) = 0.0_DP

   END IF

#endif

  END SUBROUTINE ComputeAbsorptionOpacity_TABLE


  SUBROUTINE ComputeScatteringOpacity_ES_TABLE &
               ( E, D, T, Y, Sigma, dSigmadT_Option, dSigmadY_Option )

    REAL(DP), DIMENSION(:), INTENT(in)            :: E, D, T, Y
    REAL(DP), DIMENSION(:), INTENT(out)           :: Sigma
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dSigmadT_Option
    REAL(DP), DIMENSION(:), INTENT(out), OPTIONAL :: dSigmadY_Option

    LOGICAL                :: ComputeDerivatives
    INTEGER                :: iE
    REAL(DP), DIMENSION(1) :: TMP

    ComputeDerivatives = .FALSE.

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( Sigma_T => OPACITIES % Scatt_Iso % Kernel(1) % Values(:,:,:,:,1), &
        OS      => OPACITIES % Scatt_Iso % Offsets(1,1) )

    DO iE = 1, SIZE( E )

      CALL LogInterpolateSingleVariable &
             ( [ E(iE) ] / MeV, [ D ] / ( Gram / Centimeter**3 ), &
               [ T ] / Kelvin, [ Y ], Es_T, Ds_T, Ts_T, Ys_T, &
               [ 1, 1, 1, 0 ], OS, Sigma_T, TMP )

      Sigma(iE) = TMP(1) * ( 1.0_DP / Centimeter )

    END DO

    END ASSOCIATE ! Sigma_T, etc.

#else

    Sigma(:) = 0.0_DP

    IF( ComputeDerivatives )THEN

      dSigmadT_Option(:) = 0.0_DP
      dSigmadY_Option(:) = 0.0_DP

    END IF

#endif

  END SUBROUTINE ComputeScatteringOpacity_ES_TABLE


END MODULE OpacityModule_TABLE
