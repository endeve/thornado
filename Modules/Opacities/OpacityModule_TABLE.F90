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
    OpacityTableName_Iso,  &
    OpacityTableName_NES,  &
    OpacityTableName_Pair
  INTEGER :: &
    iD_T, iT_T, iY_T
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    Es_T, Ds_T, Ts_T, Ys_T, Etas_T, &
    LogEs_T, LogDs_T, LogTs_T, LogEtas_T
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    OS_EmAb
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: &
    OS_Iso, OS_NES, OS_Pair
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: &
    EmAb_T
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, PUBLIC :: &
    Iso_T, NES_T, Pair_T
#ifdef MICROPHYSICS_WEAKLIB
  TYPE(OpacityTableType), PUBLIC :: &
    OPACITIES
#endif

  PUBLIC :: InitializeOpacities_TABLE
  PUBLIC :: FinalizeOpacities_TABLE
  PUBLIC :: ComputeAbsorptionOpacity_TABLE
  PUBLIC :: ComputeScatteringOpacity_ES_TABLE
  PUBLIC :: ComputeScatteringOpacity_NES_TABLE

#if defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
  !$ACC   OS_EmAb, OS_Iso, OS_NES, OS_Pair, &
  !$ACC   EmAb_T, Iso_T, NES_T, Pair_T )
#endif

CONTAINS


  SUBROUTINE InitializeOpacities_TABLE &
    ( OpacityTableName_EmAb_Option, OpacityTableName_Iso_Option, &
      OpacityTableName_NES_Option, OpacityTableName_Pair_Option, &
      Verbose_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_EmAb_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Iso_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_NES_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Pair_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iS, iM
    LOGICAL :: Include_NES
    LOGICAL :: Include_Pair
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

    IF( PRESENT( OpacityTableName_NES_Option ) )THEN
      OpacityTableName_NES = TRIM( OpacityTableName_NES_Option )
      Include_NES = .TRUE.
    ELSE
      Include_NES = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_Pair_Option ) )THEN
      OpacityTableName_Pair = TRIM( OpacityTableName_Pair_Option )
      Include_Pair = .TRUE.
    ELSE
      Include_Pair = .FALSE.
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
      IF( Include_NES )THEN
        WRITE(*,'(A7,A20,A)') &
          '', 'Table Name (NES):  ', TRIM( OpacityTableName_NES )
      END IF
      IF( Include_Pair )THEN
        WRITE(*,'(A7,A20,A)') &
          '', 'Table Name (Pair): ', TRIM( OpacityTableName_Pair )
      END IF
    END IF

#ifdef MICROPHYSICS_WEAKLIB

    CALL InitializeHDF( )

    IF( (.NOT.Include_NES) .AND. (.NOT.Include_Pair) )THEN

      CALL ReadOpacityTableHDF &
             ( OPACITIES, &
               FileName_EmAb_Option = TRIM( OpacityTableName_EmAb ), &
               FileName_Iso_Option  = TRIM( OpacityTableName_Iso  ) )

    ELSE IF( Include_NES .AND. (.NOT.Include_Pair) )THEN

      CALL ReadOpacityTableHDF &
             ( OPACITIES, &
               FileName_EmAb_Option = TRIM( OpacityTableName_EmAb ), &
               FileName_Iso_Option  = TRIM( OpacityTableName_Iso  ), &
               FileName_NES_Option  = TRIM( OpacityTableName_NES  ) )

    ELSE IF( (.NOT.Include_NES) .AND. Include_Pair )THEN

      CALL ReadOpacityTableHDF &
             ( OPACITIES, &
               FileName_EmAb_Option = TRIM( OpacityTableName_EmAb ), &
               FileName_Iso_Option  = TRIM( OpacityTableName_Iso  ), &
               FileName_Pair_Option = TRIM( OpacityTableName_Pair ) )

    ELSE

      CALL ReadOpacityTableHDF &
             ( OPACITIES, &
               FileName_EmAb_Option = TRIM( OpacityTableName_EmAb ), &
               FileName_Iso_Option  = TRIM( OpacityTableName_Iso  ), &
               FileName_NES_Option  = TRIM( OpacityTableName_NES  ), &
               FileName_Pair_Option = TRIM( OpacityTableName_Pair ) )

    END IF

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

    ALLOCATE( OS_EmAb(1:OPACITIES % EmAb % nOpacities) )
    OS_EmAb(:) = OPACITIES % EmAb % Offsets(:)

    ALLOCATE( OS_Iso(1:OPACITIES % Scat_Iso % nOpacities, &
                     1:OPACITIES % Scat_Iso % nMoments) )
    OS_Iso(:,:) = OPACITIES % Scat_Iso % Offsets(:,:)

    ALLOCATE( OS_NES(1:OPACITIES % Scat_NES % nOpacities, &
                     1:OPACITIES % Scat_NES % nMoments) )
    OS_NES(:,:) = OPACITIES % Scat_NES % Offsets(:,:)

    ALLOCATE( OS_Pair(1:OPACITIES % Scat_Pair % nOpacities, &
                      1:OPACITIES % Scat_Pair % nMoments) )
    OS_Pair(:,:) = OPACITIES % Scat_Pair % Offsets(:,:)

    ALLOCATE( EmAb_T(1:OPACITIES % EmAb % nPoints(1), &
                     1:OPACITIES % EmAb % nPoints(2), &
                     1:OPACITIES % EmAb % nPoints(3), &
                     1:OPACITIES % EmAb % nPoints(4), &
                     1:OPACITIES % EmAb % nOpacities) )
    DO iS = 1, OPACITIES % EmAb % nOpacities
      EmAb_T(:,:,:,:,iS) = OPACITIES % EmAb % Opacity(iS) % Values(:,:,:,:)
    END DO

    ALLOCATE( NES_T(1:OPACITIES % Scat_NES % nPoints(1), &
                    1:OPACITIES % Scat_NES % nPoints(2), &
                    1:OPACITIES % Scat_NES % nPoints(4), &
                    1:OPACITIES % Scat_NES % nPoints(5), &
                    1:OPACITIES % Scat_NES % nMoments, &
                    1:OPACITIES % Scat_NES % nOpacities) )
    DO iS = 1, OPACITIES % Scat_NES % nOpacities
      DO iM = 1, OPACITIES % Scat_NES % nMoments
        NES_T(:,:,:,:,iM,iS) = OPACITIES % Scat_NES % Kernel(iS) % Values(:,:,iM,:,:)
      END DO
    END DO

    ALLOCATE( Pair_T(1:OPACITIES % Scat_Pair % nPoints(1), &
                     1:OPACITIES % Scat_Pair % nPoints(2), &
                     1:OPACITIES % Scat_Pair % nPoints(4), &
                     1:OPACITIES % Scat_Pair % nPoints(5), &
                     1:OPACITIES % Scat_Pair % nMoments, &
                     1:OPACITIES % Scat_Pair % nOpacities) )
    DO iS = 1, OPACITIES % Scat_Pair % nOpacities
      DO iM = 1, OPACITIES % Scat_Pair % nMoments
        Pair_T(:,:,:,:,iM,iS) = OPACITIES % Scat_Pair % Kernel(iS) % Values(:,:,iM,:,:)
      END DO
    END DO

    ALLOCATE( Iso_T(1:OPACITIES % Scat_Iso % nPoints(1), &
                    1:OPACITIES % Scat_Iso % nPoints(3), &
                    1:OPACITIES % Scat_Iso % nPoints(4), &
                    1:OPACITIES % Scat_Iso % nPoints(5), &
                    1:OPACITIES % Scat_Iso % nMoments, &
                    1:OPACITIES % Scat_Iso % nOpacities) )
    DO iS = 1, OPACITIES % Scat_Iso % nOpacities
      DO iM = 1, OPACITIES % Scat_Iso % nMoments
        Iso_T(:,:,:,:,iM,iS) = OPACITIES % Scat_Iso % Kernel(iS) % Values(:,iM,:,:,:)
      END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$OMP          OS_EmAb, OS_Iso, OS_NES, OS_Pair, &
    !$OMP          EmAb_T, Iso_T, NES_T, Pair_T )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$ACC   OS_EmAb, OS_Iso, OS_NES, OS_Pair, &
    !$ACC   EmAb_T, Iso_T, NES_T, Pair_T )
#endif

#endif

  END SUBROUTINE InitializeOpacities_TABLE


  SUBROUTINE FinalizeOpacities_TABLE

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$OMP               OS_EmAb, OS_Iso, OS_NES, OS_Pair, &
    !$OMP               EmAb_T, Iso_T, NES_T, Pair_T )
#endif

    DEALLOCATE( Es_T, Ds_T, Ts_T, Ys_T, Etas_T )
    DEALLOCATE( LogEs_T, LogDs_T, LogTs_T, LogEtas_T )

    DEALLOCATE( OS_EmAb, OS_Iso, OS_NES, OS_Pair )
    DEALLOCATE( EmAb_T, Iso_T, NES_T, Pair_T )

#endif

  END SUBROUTINE FinalizeOpacities_TABLE


  SUBROUTINE ComputeAbsorptionOpacity_TABLE &
               ( E, D, T, Y, X1, X2, X3, Chi )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Chi

    REAL(DP), DIMENSION(SIZE(E)) :: LogE

#ifdef MICROPHYSICS_WEAKLIB

    IF( .NOT. InterpTest )THEN

      CALL LogInterpolateSingleVariable_1D3D &
             ( E / MeV, D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               Es_T, Ds_T, Ts_T, Ys_T, [ 1, 1, 1, 0 ], &
               OPACITIES % EmAb % Offsets(1), &
               OPACITIES % EmAb % Opacity(1) % Values, &
               Chi )

    ELSE

      LogE = LOG10( E / MeV )

      CALL LogInterpolateSingleVariable_1D3D_Custom           &
             ( LogE, LOG10( D / ( Gram / Centimeter**3 ) ), &
               LOG10( T / Kelvin ), Y, &
               LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               OPACITIES % EmAb % Offsets(1), &
               OPACITIES % EmAb % Opacity(1) % Values, &
               Chi )

    END IF

    Chi(:,:) = Chi(:,:) * ( 1.0_DP / Centimeter )

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

    IF( .NOT. InterpTest )THEN

      CALL LogInterpolateSingleVariable_1D3D &
             ( E / MeV, D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
               Es_T, Ds_T, Ts_T, Ys_T, [ 1, 1, 1, 0 ], &
               OPACITIES % Scat_Iso % Offsets(1,1), &
               OPACITIES % Scat_Iso % Kernel(1) % Values(:,1,:,:,:), &
               Sigma )

    ELSE

      LogE = LOG10( E / MeV )

      CALL LogInterpolateSingleVariable_1D3D_Custom         &
             ( LogE, LOG10( D / ( Gram / Centimeter**3 ) ), &
               LOG10( T / Kelvin ), Y, &
               LogEs_T, LogDs_T, LogTs_T, Ys_T, &
               OPACITIES % Scat_Iso % Offsets(1,1), &
               OPACITIES % Scat_Iso % Kernel(1) % Values(:,1,:,:,:), &
               Sigma )

    END IF

    Sigma(:,:) = Sigma(:,:) * ( 1.0_DP / Centimeter )

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
