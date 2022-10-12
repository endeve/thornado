MODULE OpacityModule_TABLE

#ifdef MICROPHYSICS_WEAKLIB

  ! --- weaklib modules --------------------------

  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DescribeOpacityTable
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_1D3D_Custom, &
    LogInterpolateSingleVariable_2D_Custom_Point
  USE wlGridModule, ONLY: &
    MakeLogGrid
  USE wlInterpolationUtilitiesModule, ONLY: &
    GetIndexAndDelta_Lin, &
    GetIndexAndDelta_Log

  ! ----------------------------------------------

#endif

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nE, &
    nNodesE
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE ReferenceElementModuleE, ONLY: &
    WeightsE

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: InterpTest = .TRUE.
  CHARACTER(256) :: &
    OpacityTableName_EmAb, &
    OpacityTableName_Iso,  &
    OpacityTableName_NES,  &
    OpacityTableName_Pair, &
    OpacityTableName_Brem
  INTEGER :: &
    iD_T, iT_T, iY_T, idxE1, idxE2
  REAL(DP) :: &
    dE1, dE2
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    Es_T, Ds_T, Ts_T, Ys_T, Etas_T, &
    LogEs_T, LogDs_T, LogTs_T, LogEtas_T,  &
    Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T
  REAL(DP), PUBLIC :: EC_dE
  INTEGER, PUBLIC :: EC_nE, EC_iE_max, EC_iNodeE_max
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    EC_kfmin, EC_kfmax
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    OS_EmAb
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: &
    OS_Iso, OS_NES, OS_Pair, OS_Brem
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: &
    EmAb_T
!EC table spectrum, integrated onto thornados energy elements
  INTEGER, PUBLIC :: use_EC_table
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    OS_EmAb_EC_spec
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    OS_EmAb_EC_rate
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: &
    EmAb_EC_spec_T
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: &
    EmAb_EC_rate_T
! Process_T(able), Process_A(ligned)T(able)
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, PUBLIC :: &
    Iso_T, NES_T, Pair_T, NES_AT, Pair_AT, Brem_T, Brem_AT
#ifdef MICROPHYSICS_WEAKLIB
  TYPE(OpacityTableType), PUBLIC :: &
    OPACITIES
#endif

  REAL(DP), DIMENSION(6), PUBLIC :: &
    C1, C2

  REAL(DP), PARAMETER :: cv       = 0.96d+00 ! weak interaction constant
  REAL(DP), PARAMETER :: ca       = 0.50d+00 ! weak interaction constant

  REAL(DP), PARAMETER :: C1_NuE     = ( cv + ca )**2, C2_NuE     = ( cv - ca )**2
  REAL(DP), PARAMETER :: C1_NuE_Bar = ( cv - ca )**2, C2_NuE_Bar = ( cv + ca )**2

  REAL(DP), PARAMETER :: C1_NuM     = ( cv + ca - 2.0d0 )**2, C2_NuM     = ( cv - ca         )**2
  REAL(DP), PARAMETER :: C1_NuM_Bar = ( cv - ca         )**2, C2_NuM_Bar = ( cv + ca - 2.0d0 )**2

  REAL(DP), PARAMETER :: C1_NuT     = ( cv + ca - 2.0d0 )**2, C2_NuT     = ( cv - ca         )**2
  REAL(DP), PARAMETER :: C1_NuT_Bar = ( cv - ca         )**2, C2_NuT_Bar = ( cv + ca - 2.0d0 )**2

  PUBLIC :: InitializeOpacities_TABLE
  PUBLIC :: FinalizeOpacities_TABLE
  PUBLIC :: ComputeAbsorptionOpacity_TABLE
  PUBLIC :: ComputeScatteringOpacity_ES_TABLE
  PUBLIC :: ComputeScatteringOpacity_NES_TABLE

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T,     &
  !$OMP   Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T,             &
  !$OMP   OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem,      &
  !$OMP   EmAb_T, Iso_T, NES_T, Pair_T, NES_AT, Pair_AT,  &
  !$OMP   Brem_T, Brem_AT, C1, C2,                        &
  !$OMP   use_EC_table, OS_EmAb_EC_rate, OS_EmAb_EC_spec, &
  !$OMP   EmAb_EC_rate_T, EmAb_EC_spec_T, EC_nE, EC_dE,   &
  !$OMP   EC_iE_max, EC_iNodeE_max, EC_kfmin, EC_kfmax )   
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T,     &
  !$ACC   Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T,             &
  !$ACC   OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem,      &
  !$ACC   EmAb_T, Iso_T, NES_T, Pair_T, NES_AT, Pair_AT,  & 
  !$ACC   Brem_T, Brem_AT, C1, C2,                        &
  !$ACC   use_EC_table, OS_EmAb_EC_rate, OS_EmAb_EC_spec, &
  !$ACC   EmAb_EC_rate_T, EmAb_EC_spec_T, EC_nE, EC_dE,   &
  !$ACC   EC_iE_max, EC_iNodeE_max, EC_kfmin, EC_kfmax )
#endif

CONTAINS


  SUBROUTINE InitializeOpacities_TABLE &
    ( OpacityTableName_EmAb_Option, OpacityTableName_Iso_Option, &
      OpacityTableName_NES_Option, OpacityTableName_Pair_Option, &
      OpacityTableName_Brem_Option,                              &
      EquationOfStateTableName_Option, Verbose_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_EmAb_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Iso_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_NES_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Pair_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Brem_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    CHARACTER(128)     :: EquationOfStateTableName
    REAL(DP) :: LogE1, LogE2
    INTEGER :: iS, iM, iEta, iD, iT, iN_E1, iN_E2, iE1, iE2, iNodeE1, iNodeE2, nPointsE
    LOGICAL :: Include_EmAb
    LOGICAL :: Include_Iso
    LOGICAL :: Include_NES
    LOGICAL :: Include_Pair
    LOGICAL :: Include_Brem
    LOGICAL :: Verbose

    IF( PRESENT( OpacityTableName_EmAb_Option ) &
        .AND. ( LEN( OpacityTableName_EmAb_Option ) > 1 ) )THEN
      OpacityTableName_EmAb = TRIM( OpacityTableName_EmAb_Option )
      Include_EmAb = .TRUE.
    ELSE
      OpacityTableName_EmAb = ''
      Include_EmAb = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_Iso_Option ) &
        .AND. ( LEN( OpacityTableName_Iso_Option ) > 1 ) )THEN
      OpacityTableName_Iso = TRIM( OpacityTableName_Iso_Option )
      Include_Iso = .TRUE.
    ELSE
      OpacityTableName_Iso = ''
      Include_Iso = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_NES_Option ) &
        .AND. ( LEN( OpacityTableName_NES_Option ) > 1 ) )THEN
      OpacityTableName_NES = TRIM( OpacityTableName_NES_Option )
      Include_NES = .TRUE.
    ELSE
      OpacityTableName_NES = ''
      Include_NES = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_Pair_Option ) &
        .AND. ( LEN( OpacityTableName_Pair_Option ) > 1 ) )THEN
      OpacityTableName_Pair = TRIM( OpacityTableName_Pair_Option )
      Include_Pair = .TRUE.
    ELSE
      OpacityTableName_Pair = ''
      Include_Pair = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_Brem_Option ) &
        .AND. ( LEN( OpacityTableName_Brem_Option ) > 1 ) )THEN
      OpacityTableName_Brem = TRIM( OpacityTableName_Brem_Option )
      Include_Brem = .TRUE.
    ELSE
      OpacityTableName_Brem = ''
      Include_Brem = .FALSE.
    END IF

    IF( PRESENT( EquationOfStateTableName_Option ) &
        .AND. ( LEN( EquationOfStateTableName_Option ) > 1 ) )THEN
       EquationOfStateTableName = TRIM( EquationOfStateTableName_Option )
    ELSE
       EquationOfStateTableName = 'EquationOfStateTable.h5'
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
      IF( Include_Brem )THEN
        WRITE(*,'(A7,A20,A)') &
          '', 'Table Name (Brem): ', TRIM( OpacityTableName_Brem )
      END IF
    END IF

#ifdef MICROPHYSICS_WEAKLIB

    CALL InitializeHDF( )

    CALL ReadOpacityTableHDF &
           ( OPACITIES, &
             FileName_EmAb_Option = TRIM( OpacityTableName_EmAb ), &
             FileName_Iso_Option  = TRIM( OpacityTableName_Iso  ), &
             FileName_NES_Option  = TRIM( OpacityTableName_NES  ), &
             FileName_Pair_Option = TRIM( OpacityTableName_Pair ), &
             FileName_Brem_Option = TRIM( OpacityTableName_Brem ), &
             EquationOfStateTableName_Option = EquationOfStateTableName )

    CALL FinalizeHDF( )

#ifdef THORNADO_DEBUG
    IF( Verbose ) &
      CALL DescribeOpacityTable( OPACITIES )
#endif

    nPointsE = nE * nNodesE

    C1 = [ C1_NuE, C1_NuE_Bar, &
           C1_NuM, C1_NuM_Bar, &
           C1_NuT, C1_NuT_Bar ]

    C2 = [ C2_NuE, C2_NuE_Bar, &
           C2_NuM, C2_NuM_Bar, &
           C2_NuT, C2_NuT_Bar ]

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
    OS_EmAb = OPACITIES % EmAb % Offsets

    ALLOCATE( OS_Iso(1:OPACITIES % Scat_Iso % nOpacities, &
                     1:OPACITIES % Scat_Iso % nMoments) )
    OS_Iso = OPACITIES % Scat_Iso % Offsets

    ALLOCATE( OS_NES(1:OPACITIES % Scat_NES % nOpacities, &
                     1:OPACITIES % Scat_NES % nMoments) )
    OS_NES = OPACITIES % Scat_NES % Offsets

    ALLOCATE( OS_Pair(1:OPACITIES % Scat_Pair % nOpacities, &
                      1:OPACITIES % Scat_Pair % nMoments) )
    OS_Pair = OPACITIES % Scat_Pair % Offsets

    ALLOCATE( OS_Brem(1:OPACITIES % Scat_Brem % nOpacities, &
                      1:OPACITIES % Scat_Brem % nMoments) )
    OS_Brem = OPACITIES % Scat_Brem % Offsets

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

    ALLOCATE( Brem_T(1:OPACITIES % Scat_Brem % nPoints(1), &
                     1:OPACITIES % Scat_Brem % nPoints(2), &
                     1:OPACITIES % Scat_Brem % nPoints(4), &
                     1:OPACITIES % Scat_Brem % nPoints(5), &
                     1:OPACITIES % Scat_Brem % nMoments, &
                     1:OPACITIES % Scat_Brem % nOpacities) )
    DO iS = 1, OPACITIES % Scat_Brem % nOpacities
      DO iM = 1, OPACITIES % Scat_Brem % nMoments
        Brem_T(:,:,:,:,iM,iS) = OPACITIES % Scat_Brem % Kernel(iS) % Values(:,:,iM,:,:)
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

    IF(OPACITIES % EmAb % nuclei_EC_table .gt. 0) THEN

      ALLOCATE( Ds_EC_T(OPACITIES % EmAb % EC_Table_nRho) )
      Ds_EC_T = OPACITIES % EmAb % EC_table_rho 
    
      ALLOCATE( Ts_EC_T(OPACITIES % EmAb % EC_Table_nT) )
      Ts_EC_T = OPACITIES % EmAb % EC_table_T 

      ALLOCATE( Ys_EC_T(OPACITIES % EmAb % EC_Table_nYe) )
      Ys_EC_T = OPACITIES % EmAb % EC_table_Ye 

      ALLOCATE( Es_EC_T(OPACITIES % EmAb % EC_Table_nE) )
      Es_EC_T = OPACITIES % EmAb % EC_table_E 

      EC_nE   = SIZE(Es_EC_T)
      EC_dE   = Es_EC_T(2) - Es_EC_T(1)

      ALLOCATE(EC_kfmin(nE))
      ALLOCATE(EC_kfmax(nE))

      ALLOCATE( OS_EmAb_EC_spec (OPACITIES % EmAb% EC_table_nOpacities) )
      OS_EmAb_EC_spec = OPACITIES % EmAb % EC_table_spec_Offsets(1)
      ALLOCATE( OS_EmAb_EC_rate (OPACITIES % EmAb% EC_table_nOpacities) )
      OS_EmAb_EC_rate = OPACITIES % EmAb % EC_table_rate_Offsets(1)

      ALLOCATE( EmAb_EC_spec_T (1:OPACITIES % EmAb % EC_Table_nRho, &
                                1:OPACITIES % EmAb % EC_Table_nT,   &
                                1:OPACITIES % EmAb % EC_Table_nYe,  &
                                1:OPACITIES % EmAb % EC_Table_nE) )
      EmAb_EC_spec_T(:,:,:,:) = OPACITIES % EmAB &
                              % EC_table_spec(1) % Values(:,:,:,:)      

      ALLOCATE( EmAb_EC_rate_T (1:OPACITIES % EmAb % EC_Table_nRho, &
                                1:OPACITIES % EmAb % EC_Table_nT,   &
                                1:OPACITIES % EmAb % EC_Table_nYe) )
      EmAb_EC_rate_T(:,:,:)  = OPACITIES % EmAb &
                              % EC_table_rate(1) % Values(:,:,:)

    ENDIF

    ALLOCATE( NES_AT(1:nPointsE, &
                     1:nPointsE, &
                     1:OPACITIES % Scat_NES % nPoints(4), &
                     1:OPACITIES % Scat_NES % nPoints(5), &
                     1:OPACITIES % Scat_NES % nMoments,   &
                     1:OPACITIES % Scat_NES % nOpacities) )
    NES_AT = 0.0d0

    ALLOCATE( Pair_AT(1:nPointsE, &
                      1:nPointsE, &
                      1:OPACITIES % Scat_Pair % nPoints(4), &
                      1:OPACITIES % Scat_Pair % nPoints(5), &
                      1:OPACITIES % Scat_Pair % nMoments,   &
                      1:OPACITIES % Scat_Pair % nOpacities) )
    Pair_AT = 0.0d0

    ALLOCATE( Brem_AT(1:nPointsE, &
                      1:nPointsE, &
                      1:OPACITIES % Scat_Brem % nPoints(4), &
                      1:OPACITIES % Scat_Brem % nPoints(5), &
                      1:OPACITIES % Scat_Brem % nMoments,   &
                      1:OPACITIES % Scat_Brem % nOpacities) )
    Brem_AT = 0.0d0

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP                         &
    !$OMP ( to:                                         &
    !$OMP   LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$OMP   OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem,  &
    !$OMP   EmAb_T, Iso_T, NES_T, Pair_T, Brem_T,       &
    !$OMP   NES_AT, Pair_AT, Brem_AT, C1, C2 )
    !$OMP TARGET UPDATE TO                              &
    !$OMP ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$OMP   OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem,  &
    !$OMP   EmAb_T, Iso_T, NES_T, Pair_T, Brem_T,       &
    !$OMP   NES_AT, Pair_AT, Brem_AT, C1, C2            )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$ACC   OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem,  &
    !$ACC   EmAb_T, Iso_T, NES_T, Pair_T, Brem_T,       &
    !$ACC   NES_AT, Pair_AT, Brem_AT, C1, C2            )
#endif

    Energy_Grid: BLOCK
    REAL(dp) :: CenterE(nE), WidthE(nE), NodesE(nNodesE)

    CenterE(:) = MeshE % Center(:)
    WidthE(:)  = MeshE % Width(:)
    NodesE(:)  = MeshE % Nodes(:)

    EC_TABLE: BLOCK

    use_EC_table = OPACITIES % EmAb % nuclei_EC_table

    IF( use_EC_table .gt. 0 ) THEN
write(*,*) 'lol EC table'

    build_kfmin_max: BLOCK

      REAL(dp), DIMENSION(nE)   :: E_cells, dE_cells
      REAL(dp), DIMENSION(nE+1) :: E_faces
      INTEGER                   :: iE1, iNodeE1
      INTEGER                   :: k, kk

      REAL(dp)                  :: EC_E_max

      REAL(dp)                  :: x

      EC_E_max = Es_EC_T(EC_nE)

      DO k = 1, nE
        E_cells(k)  = CenterE(k) / MeV
        dE_cells(k) = WidthE(k)  / MeV
        E_faces(k)  = (CenterE(k) - 0.5d0*WidthE(k)) / MeV
      ENDDO 
      E_faces(nE+1)  = (CenterE(nE) + 0.5d0*WidthE(nE)) / MeV

      EC_iE_max     = 1
      EC_iNodeE_max = 1
      DO kk = 1, nPointsE
        iE1     = MOD( (kk-1) / nNodesE, nE      ) + 1
        iNodeE1 = MOD( (kk-1)          , nNodesE ) + 1
        x = NodeCoordinate( CenterE(iE1), WidthE(iE1), NodesE(iNodeE1) ) / MeV
        if(E_faces(iE1)<=EC_E_max) EC_iE_max     = iE1
        if(x <= EC_E_max)          EC_iNodeE_max = kk
      ENDDO

      EC_kfmin = 1
      EC_kfmax = 1
 
      DO k = 1, EC_iE_max
        DO kk = 1, EC_nE
          IF(E_faces(k) .ge. Es_EC_T(kk)) THEN
            EC_kfmin(k) = kk
          ELSE
            EXIT
          ENDIF
        ENDDO
        DO kk = EC_nE, 1, -1
          IF(E_faces(k+1) .le. Es_EC_T(kk)) THEN
            EC_kfmax(k) = kk-1
          ELSE
            EXIT
          ENDIF
        ENDDO
        IF (EC_kfmax(k).gt.EC_nE) EC_kfmax(k)=EC_nE
        IF (EC_kfmin(k).gt.EC_nE) EC_kfmin(k)=EC_nE
        IF (EC_kfmax(k).lt.1)     EC_kfmax(k)=1
        IF (EC_kfmin(k).lt.1)     EC_kfmin(k)=1
      ENDDO

      EC_kfmax(EC_iE_max) = EC_nE

    END BLOCK build_kfmin_max

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA MAP                       & 
    !$OMP ( to:                                       &
    !$OMP   Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T,       &
    !$OMP   OS_EmAb_EC_rate, OS_EmAb_EC_spec,         &
    !$OMP   EmAb_EC_rate_T, EmAb_EC_spec_T,           &
    !$OMP   EC_nE, EC_dE, EC_iE_max, EC_iNodeE_max,   & 
    !$OMP   EC_kfmin, EC_kfmax, use_EC_table )
    !$OMP TARGET UPDATE TO                            &
    !$OMP ( Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T,       & 
    !$OMP   OS_EmAb_EC_rate, OS_EmAb_EC_spec,         &
    !$OMP   EmAb_EC_rate_T, EmAb_EC_spec_T,           &
    !$OMP   EC_nE, EC_dE, EC_iE_max, EC_iNodeE_max,   &
    !$OMP   EC_kfmin, EC_kfmax, use_EC_table )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE                               &
    !$ACC ( Ds_EC_T, Ts_EC_T, Ys_EC_T, Es_EC_T,       &
    !$ACC   OS_EmAb_EC_rate, OS_EmAb_EC_spec,         &
    !$ACC   EmAb_EC_rate_T, EmAb_EC_spec_T,           &
    !$ACC   EC_nE, EC_dE, EC_iE_max, EC_iNodeE_max,   &
    !$ACC   EC_kfmin, EC_kfmax, use_EC_table )
#endif

     ENDIF

    END BLOCK EC_TABLE

    NES: BLOCK

    INTEGER :: nSpecies, nMoments, nPointsEta, nPointsT

    nSpecies   = OPACITIES % Scat_NES % nOpacities
    nMoments   = OPACITIES % Scat_NES % nMoments
    nPointsEta = OPACITIES % Scat_NES % nPoints(5)
    nPointsT   = OPACITIES % Scat_NES % nPoints(4)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 ) &
    !$ACC PRESENT( LogEs_T, OS_NES, NES_T, NES_AT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#endif
    DO iS = 1, nSpecies !OPACITIES % Scat_NES % nOpacities
      DO iM = 1, nMoments !OPACITIES % Scat_NES % nMoments
        DO iEta = 1, nPointsEta !OPACITIES % Scat_NES % nPoints(5)
          DO iT = 1, nPointsT !OPACITIES % Scat_NES % nPoints(4)
            DO iN_E2 = 1, nPointsE
              DO iN_E1 = 1, nPointsE

                iE1     = MOD( (iN_E1-1) / nNodesE, nE      ) + 1
                iNodeE1 = MOD( (iN_E1-1)          , nNodesE ) + 1

                iE2     = MOD( (iN_E2-1) / nNodesE, nE      ) + 1
                iNodeE2 = MOD( (iN_E2-1)          , nNodesE ) + 1

                LogE1 = LOG10( NodeCoordinate( CenterE(iE1), WidthE(iE1), NodesE(iNodeE1) ) / MeV )
                LogE2 = LOG10( NodeCoordinate( CenterE(iE2), WidthE(iE2), NodesE(iNodeE2) ) / MeV )

                CALL LogInterpolateSingleVariable_2D_Custom_Point &
                       ( LogE1, LogE2, LogEs_T, LogEs_T, OS_NES(iS,iM), NES_T(:,:,iT,iEta,iM,iS), &
                         NES_AT(iN_E1,iN_E2,iT,iEta,iM,iS) )

                NES_AT(iN_E1,iN_E2,iT,iEta,iM,iS) &
                  = LOG10( NES_AT(iN_E1,iN_E2,iT,iEta,iM,iS) + OS_NES(iS,iM) )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    END BLOCK NES

    Pair: BLOCK

    INTEGER :: nSpecies, nMoments, nPointsEta, nPointsT

    nSpecies   = OPACITIES % Scat_Pair % nOpacities
    nMoments   = OPACITIES % Scat_Pair % nMoments
    nPointsEta = OPACITIES % Scat_Pair % nPoints(5)
    nPointsT   = OPACITIES % Scat_Pair % nPoints(4)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 ) &
    !$ACC PRESENT( LogEs_T, OS_Pair, Pair_T, Pair_AT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#endif
    DO iS = 1, nSpecies !OPACITIES % Scat_Pair % nOpacities
      DO iM = 1, nMoments !OPACITIES % Scat_Pair % nMoments
        DO iEta = 1, nPointsEta !OPACITIES % Scat_Pair % nPoints(5)
          DO iT = 1, nPointsT !OPACITIES % Scat_Pair % nPoints(4)
            DO iN_E2 = 1, nPointsE
              DO iN_E1 = 1, nPointsE

                iE1     = MOD( (iN_E1-1) / nNodesE, nE      ) + 1
                iNodeE1 = MOD( (iN_E1-1)          , nNodesE ) + 1

                iE2     = MOD( (iN_E2-1) / nNodesE, nE      ) + 1
                iNodeE2 = MOD( (iN_E2-1)          , nNodesE ) + 1

                LogE1 = LOG10( NodeCoordinate( CenterE(iE1), WidthE(iE1), NodesE(iNodeE1) ) / MeV )
                LogE2 = LOG10( NodeCoordinate( CenterE(iE2), WidthE(iE2), NodesE(iNodeE2) ) / MeV )

                CALL LogInterpolateSingleVariable_2D_Custom_Point &
                       ( LogE1, LogE2, LogEs_T, LogEs_T, OS_Pair(iS,iM), Pair_T(:,:,iT,iEta,iM,iS), &
                         Pair_AT(iN_E1,iN_E2,iT,iEta,iM,iS) )

                Pair_AT(iN_E1,iN_E2,iT,iEta,iM,iS) &
                  = LOG10( Pair_AT(iN_E1,iN_E2,iT,iEta,iM,iS) + OS_Pair(iS,iM) )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    END BLOCK Pair

    Brem: BLOCK

    INTEGER :: nSpecies, nMoments, nPointsD, nPointsT

    nSpecies = OPACITIES % Scat_Brem % nOpacities
    nMoments = OPACITIES % Scat_Brem % nMoments
    nPointsD = OPACITIES % Scat_Brem % nPoints(4)
    nPointsT = OPACITIES % Scat_Brem % nPoints(5)

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 ) &
    !$ACC PRESENT( LogEs_T, OS_Brem, Brem_T, Brem_AT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#endif
    DO iS = 1, nSpecies !OPACITIES % Scat_Brem % nOpacities
      DO iM = 1, nMoments !OPACITIES % Scat_Brem % nMoments
        DO iD = 1, nPointsD !OPACITIES % Scat_Brem % nPoints(4)
          DO iT = 1, nPointsT !OPACITIES % Scat_Brem % nPoints(5)
            DO iN_E2 = 1, nPointsE
              DO iN_E1 = 1, nPointsE 

                iE1     = MOD( (iN_E1-1) / nNodesE, nE      ) + 1
                iNodeE1 = MOD( (iN_E1-1)          , nNodesE ) + 1

                iE2     = MOD( (iN_E2-1) / nNodesE, nE      ) + 1
                iNodeE2 = MOD( (iN_E2-1)          , nNodesE ) + 1

                LogE1 = LOG10( NodeCoordinate( CenterE(iE1), WidthE(iE1), NodesE(iNodeE1) ) / MeV )
                LogE2 = LOG10( NodeCoordinate( CenterE(iE2), WidthE(iE2), NodesE(iNodeE2) ) / MeV )

                CALL LogInterpolateSingleVariable_2D_Custom_Point &
                       ( LogE1, LogE2, LogEs_T, LogEs_T, OS_Brem(iS,iM), Brem_T(:,:,iD,iT,iM,iS), &
                         Brem_AT(iN_E1,iN_E2,iD,iT,iM,iS) )

                Brem_AT(iN_E1,iN_E2,iD,iT,iM,iS) &
                  = LOG10( Brem_AT(iN_E1,iN_E2,iD,iT,iM,iS) + OS_Brem(iS,iM) )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    END BLOCK Brem

    END BLOCK Energy_Grid

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM &
    !$OMP ( NES_AT, Pair_AT, Brem_AT )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST &
    !$ACC ( NES_AT, Pair_AT, Brem_AT )
#endif

#endif

  END SUBROUTINE InitializeOpacities_TABLE


  SUBROUTINE FinalizeOpacities_TABLE

#ifdef MICROPHYSICS_WEAKLIB

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$OMP               OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem, &
    !$OMP               EmAb_T, Iso_T, NES_T, Pair_T, Brem_T, &
    !$OMP               NES_AT, Pair_AT, Brem_AT, C1, C2 )
#endif

    DEALLOCATE( Es_T, Ds_T, Ts_T, Ys_T, Etas_T )
    DEALLOCATE( LogEs_T, LogDs_T, LogTs_T, LogEtas_T )

    DEALLOCATE( OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem )
    DEALLOCATE( EmAb_T, Iso_T, NES_T, Pair_T, Brem_T )
    DEALLOCATE( NES_AT, Pair_AT, Brem_AT )

    IF(ALLOCATED(OS_EmAb_EC_spec)) DEALLOCATE(OS_EmAb_EC_spec)
    IF(ALLOCATED(OS_EmAb_EC_rate)) DEALLOCATE(OS_EmAb_EC_rate)
    IF(ALLOCATED(EmAb_EC_rate_T )) DEALLOCATE(EmAb_EC_rate_T)
    IF(ALLOCATED(EmAb_EC_spec_T )) DEALLOCATE(EmAb_EC_spec_T )
    IF(ALLOCATED(Ds_EC_T))         DEALLOCATE(Ds_EC_T)
    IF(ALLOCATED(Ts_EC_T))         DEALLOCATE(Ts_EC_T)
    IF(ALLOCATED(Ys_EC_T))         DEALLOCATE(Ys_EC_T)
    IF(ALLOCATED(Es_EC_T))         DEALLOCATE(Es_EC_T)
    IF(ALLOCATED(EC_kfmin))        DEALLOCATE(EC_kfmin)
    IF(ALLOCATED(EC_kfmax))        DEALLOCATE(EC_kfmax)

#endif

  END SUBROUTINE FinalizeOpacities_TABLE


  SUBROUTINE ComputeAbsorptionOpacity_TABLE &
               ( E, D, T, Y, X1, X2, X3, Chi )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Chi

    REAL(DP), DIMENSION(SIZE(E)) :: LogE

#ifdef MICROPHYSICS_WEAKLIB

!!$    IF( .NOT. InterpTest )THEN
!!$
!!$      CALL LogInterpolateSingleVariable_1D3D &
!!$             ( E / MeV, D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
!!$               Es_T, Ds_T, Ts_T, Ys_T, [ 1, 1, 1, 0 ], &
!!$               OPACITIES % EmAb % Offsets(1), &
!!$               OPACITIES % EmAb % Opacity(1) % Values, &
!!$               Chi )
!!$
!!$    ELSE
!!$
!!$      LogE = LOG10( E / MeV )
!!$
!!$      CALL LogInterpolateSingleVariable_1D3D_Custom           &
!!$             ( LogE, LOG10( D / ( Gram / Centimeter**3 ) ), &
!!$               LOG10( T / Kelvin ), Y, &
!!$               LogEs_T, LogDs_T, LogTs_T, Ys_T, &
!!$               OPACITIES % EmAb % Offsets(1), &
!!$               OPACITIES % EmAb % Opacity(1) % Values, &
!!$               Chi )
!!$
!!$    END IF
!!$
!!$    Chi(:,:) = Chi(:,:) * ( 1.0_DP / Centimeter )

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

!!$    IF( .NOT. InterpTest )THEN
!!$
!!$      CALL LogInterpolateSingleVariable_1D3D &
!!$             ( E / MeV, D / ( Gram / Centimeter**3 ), T / Kelvin, Y, &
!!$               Es_T, Ds_T, Ts_T, Ys_T, [ 1, 1, 1, 0 ], &
!!$               OPACITIES % Scat_Iso % Offsets(1,1), &
!!$               OPACITIES % Scat_Iso % Kernel(1) % Values(:,1,:,:,:), &
!!$               Sigma )
!!$
!!$    ELSE
!!$
!!$      LogE = LOG10( E / MeV )
!!$
!!$      CALL LogInterpolateSingleVariable_1D3D_Custom         &
!!$             ( LogE, LOG10( D / ( Gram / Centimeter**3 ) ), &
!!$               LOG10( T / Kelvin ), Y, &
!!$               LogEs_T, LogDs_T, LogTs_T, Ys_T, &
!!$               OPACITIES % Scat_Iso % Offsets(1,1), &
!!$               OPACITIES % Scat_Iso % Kernel(1) % Values(:,1,:,:,:), &
!!$               Sigma )
!!$
!!$    END IF
!!$
!!$    Sigma(:,:) = Sigma(:,:) * ( 1.0_DP / Centimeter )

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
