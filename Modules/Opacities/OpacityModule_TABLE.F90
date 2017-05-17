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
    LogInterpolateSingleVariable_1D3D, &
    LogInterpolateSingleVariable_2D2D

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
    Es_T, Ds_T, Ts_T, Ys_T, Etas_T
#ifdef MICROPHYSICS_WEAKLIB
  TYPE(OpacityTableType) :: &
    OPACITIES
#endif

  PUBLIC :: InitializeOpacities_TABLE
  PUBLIC :: FinalizeOpacities_TABLE
  PUBLIC :: ComputeAbsorptionOpacity_TABLE
  PUBLIC :: ComputeScatteringOpacity_ES_TABLE
  PUBLIC :: ComputeScatteringOpacity_NES_TABLE

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

    ! --- Eta Grid ---

    ALLOCATE( Etas_T(OPACITIES % EtaGrid % nPoints) )
    Etas_T = OPACITIES % EtaGrid  % Values

#endif

  END SUBROUTINE InitializeOpacities_TABLE


  SUBROUTINE FinalizeOpacities_TABLE

#ifdef MICROPHYSICS_WEAKLIB

    DEALLOCATE( Es_T, Ds_T, Ts_T, Ys_T, Etas_T )

#endif

  END SUBROUTINE FinalizeOpacities_TABLE


  SUBROUTINE ComputeAbsorptionOpacity_TABLE( E, D, T, Y, Chi )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Chi

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( Chi_T => OPACITIES % ThermEmAb % Absorptivity(1) % Values, &
        OS    => OPACITIES % ThermEmAb % Offsets(1) )

    CALL LogInterpolateSingleVariable_1D3D &
           ( E / MeV, D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
             Es_T, Ds_T, Ts_T, Ys_T, [ 1, 1, 1, 0 ], OS, Chi_T, Chi )

    Chi(:,:) = Chi(:,:) * ( 1.0_DP / Centimeter )

    END ASSOCIATE ! Chi_T, etc.

#else

    Chi(:,:) = 0.0_DP

#endif

  END SUBROUTINE ComputeAbsorptionOpacity_TABLE


  SUBROUTINE ComputeScatteringOpacity_ES_TABLE( E, D, T, Y, Sigma )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Sigma

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( Sigma_T => OPACITIES % Scatt_Iso % Kernel(1) % Values(:,:,:,:,1), &
        OS      => OPACITIES % Scatt_Iso % Offsets(1,1) )

    CALL LogInterpolateSingleVariable_1D3D &
           ( E / MeV, D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
             Es_T, Ds_T, Ts_T, Ys_T, [ 1, 1, 1, 0 ], OS, Sigma_T, Sigma )

    Sigma(:,:) = Sigma(:,:) * ( 1.0_DP / Centimeter )

    END ASSOCIATE ! Sigma_T, etc.

#else

    Sigma(:,:) = 0.0_DP

#endif

  END SUBROUTINE ComputeScatteringOpacity_ES_TABLE


  SUBROUTINE ComputeScatteringOpacity_NES_TABLE( E, T, Eta, R0_In, R0_Out )

    REAL(DP), DIMENSION(:),     INTENT(in)  :: E, T, Eta
    REAL(DP), DIMENSION(:,:,:), INTENT(out) :: R0_In, R0_Out

    INTEGER :: iX

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( R0_Out_T => OPACITIES % Scatt_NES % Kernel(1) % Values(:,:,:,:,1), &
        OS       => OPACITIES % Scatt_NES % Offsets(1,1) )

    CALL LogInterpolateSingleVariable_2D2D      &
           ( E / MeV, E / MeV, T / Kelvin, Eta, &
             Es_T, Es_T, Ts_T, Etas_T,          &
             [ 1, 1, 1, 1 ], OS, R0_Out_T, R0_Out )

    END ASSOCIATE ! R0_Out_T, etc.

    R0_Out = R0_Out * ( 1.0_DP / ( Centimeter * MeV**3 ) )

    DO iX = 1, SIZE( T )

      R0_In(:,:,iX) = TRANSPOSE( R0_Out(:,:,iX) )

    END DO

#else

  R0_In (:,:,:) = 0.0_DP
  R0_Out(:,:,:) = 0.0_DP

#endif

  END SUBROUTINE ComputeScatteringOpacity_NES_TABLE


END MODULE OpacityModule_TABLE
