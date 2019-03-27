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
    LogInterpolateSingleVariable_1D3D, &
    LogInterpolateSingleVariable_1D3D_Custom, &
    LogInterpolateSingleVariable_2D2D, &
    LogInterpolateSingleVariable_2D2D_Custom

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

  LOGICAL, PARAMETER :: InterpTest = .TRUE.
  CHARACTER(256) :: &
    OpacityTableName_EmAb, &
    OpacityTableName_Iso
  INTEGER :: &
    iD_T, iT_T, iY_T
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    Es_T, Ds_T, Ts_T, Ys_T, Etas_T, &
    LogEs_T, LogDs_T, LogTs_T, LogEtas_T
#ifdef MICROPHYSICS_WEAKLIB
  TYPE(OpacityTableType), PUBLIC :: &
    OPACITIES
#endif

  PUBLIC :: InitializeOpacities_TABLE
  PUBLIC :: FinalizeOpacities_TABLE
  PUBLIC :: ComputeAbsorptionOpacity_TABLE
  PUBLIC :: ComputeScatteringOpacity_ES_TABLE
  PUBLIC :: ComputeScatteringOpacity_NES_TABLE

CONTAINS


  SUBROUTINE InitializeOpacities_TABLE &
    ( OpacityTableName_EmAb_Option, OpacityTableName_Iso_Option, &
      Verbose_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_EmAb_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Iso_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( OpacityTableName_EmAb_Option ) )THEN
      OpacityTableName_EmAb = TRIM( OpacityTableName_EmAb_Option )
    ELSE
      OpacityTableName_EmAb = 'OpacityTable_EmAb.h5'
    END IF

    IF( PRESENT( OpacityTableName_Iso_Option ) )THEN
      OpacityTableName_Iso = TRIM( OpacityTableName_Iso_Option )
    ELSE
      OpacityTableName_Iso = 'OpacityTable_Iso.h5'
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A7,A20,A)') &
        '', 'Table Name (EmAb): ', TRIM( OpacityTableName_EmAb )
      WRITE(*,'(A7,A20,A)') &
        '', 'Table Name (Iso):  ', TRIM( OpacityTableName_Iso )
    END IF

#ifdef MICROPHYSICS_WEAKLIB

    CALL InitializeHDF( )

    CALL ReadOpacityTableHDF &
           ( OPACITIES, &
             TRIM( OpacityTableName_EmAb ), &
             TRIM( OpacityTableName_Iso  ) )

    CALL FinalizeHDF( )

    ! --- Thermodynamic State Indices ---

    iD_T = OPACITIES % TS % Indices % iRho
    iT_T = OPACITIES % TS % Indices % iT
    iY_T = OPACITIES % TS % Indices % iYe

    ! --- Thermodynamic States ---

    ! --- Density ---

    ALLOCATE( Ds_T(OPACITIES % TS % nPoints(iD_T)) )
    Ds_T = OPACITIES % TS % States(iD_T) % Values

    ALLOCATE( LogDs_T(SIZE( Ds_T )) )
    LogDs_T = LOG10( Ds_T )

    ! --- Temperature ---

    ALLOCATE( Ts_T(OPACITIES % TS % nPoints(iT_T)) )
    Ts_T = OPACITIES % TS % States(iT_T) % Values

    ALLOCATE( LogTs_T(SIZE( Ts_T )) )
    LogTs_T = LOG10( Ts_T )

    ! --- Electron Fraction ---

    ALLOCATE( Ys_T(OPACITIES % TS % nPoints(iY_T)) )
    Ys_T = OPACITIES % TS % States(iY_T) % Values

    ! --- Energy Grid ---

    ALLOCATE( Es_T(OPACITIES % EnergyGrid % nPoints) )
    Es_T = OPACITIES % EnergyGrid  % Values

    ALLOCATE( LogEs_T(SIZE( Es_T )) )
    LogEs_T = LOG10( Es_T )

    ! --- Eta Grid ---

    ALLOCATE( Etas_T(OPACITIES % EtaGrid % nPoints) )
    Etas_T = OPACITIES % EtaGrid  % Values

    ALLOCATE( LogEtas_T(SIZE( Etas_T )) )
    LogEtas_T = LOG10( Etas_T )

#endif

  END SUBROUTINE InitializeOpacities_TABLE


  SUBROUTINE FinalizeOpacities_TABLE

#ifdef MICROPHYSICS_WEAKLIB

    DEALLOCATE( Es_T, Ds_T, Ts_T, Ys_T, Etas_T )
    DEALLOCATE( LogEs_T, LogDs_T, LogTs_T, LogEtas_T )

#endif

  END SUBROUTINE FinalizeOpacities_TABLE


  SUBROUTINE ComputeAbsorptionOpacity_TABLE &
               ( E, D, T, Y, X1, X2, X3, Chi )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Chi

    REAL(DP), DIMENSION(SIZE(E)) :: LogE

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( Chi_T => OPACITIES % EmAb % Opacity(1) % Values, &
        OS    => OPACITIES % EmAb % Offsets(1) )

    IF( .NOT. InterpTest )THEN

      CALL LogInterpolateSingleVariable_1D3D &
             ( E / MeV, D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               Es_T, Ds_T, Ts_T, Ys_T, [ 1, 1, 1, 0 ], OS, Chi_T, Chi )

    ELSE

      LogE = LOG10( E / MeV )

      CALL LogInterpolateSingleVariable_1D3D_Custom           &
             ( LogE, LOG10( D / ( Gram / Centimeter**3 ) ), &
               LOG10( T / Kelvin ), Y, &
               LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, &
               Chi_T, Chi )

    END IF

    Chi(:,:) = Chi(:,:) * ( 1.0_DP / Centimeter )

    END ASSOCIATE ! Chi_T, etc.

#else

    Chi(:,:) = 0.0_DP

#endif

  END SUBROUTINE ComputeAbsorptionOpacity_TABLE


  SUBROUTINE ComputeScatteringOpacity_ES_TABLE &
               ( E, D, T, Y, X1, X2, X3, Sigma )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Sigma

    REAL(DP), DIMENSION(SIZE(E)) :: LogE

#ifdef MICROPHYSICS_WEAKLIB

    ASSOCIATE &
      ( Sigma_T => OPACITIES % Scat_Iso % Kernel(1) % Values(:,1,:,:,:), &
        OS      => OPACITIES % Scat_Iso % Offsets(1,1) )

    IF( .NOT. InterpTest )THEN

      CALL LogInterpolateSingleVariable_1D3D &
             ( E / MeV, D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               Es_T, Ds_T, Ts_T, Ys_T, [ 1, 1, 1, 0 ], OS, Sigma_T, Sigma )

    ELSE

      LogE = LOG10( E / MeV )

      CALL LogInterpolateSingleVariable_1D3D_Custom         &
             ( LogE, LOG10( D / ( Gram / Centimeter**3 ) ), &
               LOG10( T / Kelvin ), Y, &
               LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, &
               Sigma_T, Sigma )

    END IF

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
    REAL(DP), DIMENSION(SIZE(E)) :: LogE

    PRINT*, "ComputeScatteringOpacity_NES_TABLE Disabled"
    STOP

#ifdef MICROPHYSICS_WEAKLIB

!!$    IF( .NOT. InterpTest )THEN
!!$
!!$      CALL LogInterpolateSingleVariable_2D2D &
!!$             ( E / MeV, E / MeV, T / Kelvin, Eta, &
!!$               Es_T, Es_T, Ts_T, Etas_T, [ 1, 1, 1, 1 ], &
!!$               OPACITIES % Scatt_NES % Offsets(1,1), &
!!$               OPACITIES % Scatt_NES % Kernel(1) % Values(:,:,:,:,1), &
!!$               R0_Out )
!!$
!!$    ELSE
!!$
!!$      LogE = LOG10( E / MeV )
!!$
!!$      CALL LogInterpolateSingleVariable_2D2D_Custom &
!!$             ( LogE, LogE, LOG10( T / Kelvin ), LOG10( Eta ), &
!!$               LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
!!$               OPACITIES % Scatt_NES % Offsets(1,1), &
!!$               OPACITIES % Scatt_NES % Kernel(1) % Values(:,:,:,:,1), &
!!$               R0_Out )
!!$
!!$    END IF
!!$
!!$    R0_Out = R0_Out * ( 1.0_DP / ( Centimeter * MeV**3 ) )
!!$
!!$    DO iX = 1, SIZE( T )
!!$
!!$      R0_In(:,:,iX) = TRANSPOSE( R0_Out(:,:,iX) )
!!$
!!$    END DO

#else

  R0_In (:,:,:) = 0.0_DP
  R0_Out(:,:,:) = 0.0_DP

#endif

  END SUBROUTINE ComputeScatteringOpacity_NES_TABLE


END MODULE OpacityModule_TABLE
