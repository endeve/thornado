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
    DeAllocateOpacityTable
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_1D3D_Custom, &
    LogInterpolateSingleVariable_2D_Custom_Point

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
    nOpacities_NES, nMoments_NES, nPointsT_NES, nPointsEta_NES
  INTEGER :: &
    nOpacities_Pair, nMoments_Pair, nPointsT_Pair, nPointsEta_Pair
  INTEGER :: &
    nOpacities_Brem, nMoments_Brem, nPointsD_Brem, nPointsT_Brem
  INTEGER :: &
    iD_T, iT_T, iY_T, idxE1, idxE2
  REAL(DP) :: &
    dE1, dE2
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    Es_T, Ds_T, Ts_T, Ys_T, Etas_T, &
    LogEs_T, LogDs_T, LogTs_T, LogEtas_T
  REAL(DP), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
    OS_EmAb
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: &
    OS_Iso, OS_NES, OS_Pair, OS_Brem
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: &
    EmAb_T
! Process_T(able), Process_A(ligned)T(able)
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, PUBLIC :: &
    Iso_T, NES_T, Pair_T, NES_AT, Pair_AT, Brem_T, Brem_AT
#ifdef MICROPHYSICS_WEAKLIB
  TYPE(OpacityTableType), PUBLIC :: &
    OPACITIES
#endif
  LOGICAL :: Use_OpacityTables

  REAL(DP), DIMENSION(6), PUBLIC :: &
    C1, C2

  REAL(DP) :: EmAb_Nucleon_MinD, EmAb_Nucleon_MaxD ! density cutoffs for EmAb (free nucleon) opacities
  REAL(DP) :: EmAb_Nuclei_MinD, EmAb_Nuclei_MaxD   ! density cutoffs for EmAb (nuclei) opacities
  REAL(DP) :: EmAb_MinD, EmAb_MaxD                 ! density cutoffs for all EmAb opacities
  REAL(DP) :: Iso_MinD, Iso_MaxD                   ! density cutoffs for all Iso opacities
  REAL(DP) :: NES_MinD, NES_MaxD                   ! density cutoffs for all NES opacities
  REAL(DP) :: Pair_MinD, Pair_MaxD                 ! density cutoffs for all Pair opacities
  REAL(DP) :: Brem_MinD, Brem_MaxD                 ! density cutoffs for all Brem opacities
  REAL(DP) :: NNS_MinD, NNS_MaxD                   ! density cutoffs for all NNS opacities
  REAL(DP) :: NuPair_MinD, NuPair_MaxD             ! density cutoffs for all NuPair opacities
  REAL(DP) :: Op_MinD, Op_MaxD                     ! density cutoffs for all opacities

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

  PUBLIC :: QueryOpacity_EmAb_Nucleon
  PUBLIC :: QueryOpacity_EmAb_Nuclei
  PUBLIC :: QueryOpacity_EmAb
  PUBLIC :: QueryOpacity_Iso
  PUBLIC :: QueryOpacity_NES
  PUBLIC :: QueryOpacity_Pair
  PUBLIC :: QueryOpacity_Brem
  PUBLIC :: QueryOpacity_NNS
  PUBLIC :: QueryOpacity_NuPair
  PUBLIC :: QueryOpacity

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T,    &
  !$OMP   OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem,     &
  !$OMP   EmAb_T, Iso_T, NES_T, Pair_T, NES_AT, Pair_AT, &
  !$OMP   Brem_T, Brem_AT, C1, C2,                       &
  !$OMP   EmAb_Nucleon_MinD, EmAb_Nucleon_MaxD,          &
  !$OMP   EmAb_Nuclei_MinD, EmAb_Nuclei_MaxD,            &
  !$OMP   EmAb_MinD, EmAb_MaxD,                          &
  !$OMP   Iso_MinD, Iso_MaxD,                            &
  !$OMP   NES_MinD, NES_MaxD,                            &
  !$OMP   Pair_MinD, Pair_MaxD,                          &
  !$OMP   Brem_MinD, Brem_MaxD,                          &
  !$OMP   NNS_MinD, NNS_MaxD,                            &
  !$OMP   NuPair_MinD, NuPair_MaxD,                      &
  !$OMP   Op_MinD, Op_MaxD )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T,    &
  !$ACC   OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem,     &
  !$ACC   EmAb_T, Iso_T, NES_T, Pair_T, NES_AT, Pair_AT, & 
  !$ACC   Brem_T, Brem_AT, C1, C2,                       &
  !$ACC   EmAb_Nucleon_MinD, EmAb_Nucleon_MaxD,          &
  !$ACC   EmAb_Nuclei_MinD, EmAb_Nuclei_MaxD,            &
  !$ACC   EmAb_MinD, EmAb_MaxD,                          &
  !$ACC   Iso_MinD, Iso_MaxD,                            &
  !$ACC   NES_MinD, NES_MaxD,                            &
  !$ACC   Pair_MinD, Pair_MaxD,                          &
  !$ACC   Brem_MinD, Brem_MaxD,                          &
  !$ACC   NNS_MinD, NNS_MaxD,                            &
  !$ACC   NuPair_MinD, NuPair_MaxD,                      &
  !$ACC   Op_MinD, Op_MaxD )
#endif

CONTAINS


  SUBROUTINE InitializeOpacities_TABLE &
    ( OpacityTableName_EmAb_Option, OpacityTableName_Iso_Option, &
      OpacityTableName_NES_Option, OpacityTableName_Pair_Option, &
      OpacityTableName_Brem_Option,                              &
      EmAb_Nucleon_MinD_Option, EmAb_Nucleon_MaxD_Option,        &
      EmAb_Nuclei_MinD_Option, EmAb_Nuclei_MaxD_Option,          &
      EmAb_MinD_Option, EmAb_MaxD_Option,                        &
      Iso_MinD_Option, Iso_MaxD_Option,                          &
      NES_MinD_Option, NES_MaxD_Option,                          &
      Pair_MinD_Option, Pair_MaxD_Option,                        &
      Brem_MinD_Option, Brem_MaxD_Option,                        &
      NNS_MinD_Option, NNS_MaxD_Option,                          &
      NuPair_MinD_Option, NuPair_MaxD_Option,                    &
      Op_MinD_Option, Op_MaxD_Option,                            &
      EquationOfStateTableName_Option, Verbose_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_EmAb_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Iso_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_NES_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Pair_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Brem_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: EquationOfStateTableName_Option
    REAL(DP),         INTENT(in), OPTIONAL :: EmAb_Nucleon_MinD_Option, EmAb_Nucleon_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: EmAb_Nuclei_MinD_Option, EmAb_Nuclei_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: EmAb_MinD_Option, EmAb_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Iso_MinD_Option, Iso_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: NES_MinD_Option, NES_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Pair_MinD_Option, Pair_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Brem_MinD_Option, Brem_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: NNS_MinD_Option, NNS_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: NuPair_MinD_Option, NuPair_MaxD_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Op_MinD_Option, Op_MaxD_Option
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    CHARACTER(128)     :: EquationOfStateTableName
    REAL(DP) :: LogE1, LogE2
    INTEGER :: iS, iM, iEta, iD, iT, iN_E1, iN_E2, iE1, iE2, iNodeE1, iNodeE2
    INTEGER :: nOpacities, nMoments, nPointsEta, nPointsD, nPointsT, nPointsE
    LOGICAL :: Include_EmAb
    LOGICAL :: Include_Iso
    LOGICAL :: Include_NES
    LOGICAL :: Include_Pair
    LOGICAL :: Include_Brem
    LOGICAL :: Verbose

    IF( PRESENT( OpacityTableName_EmAb_Option ) &
        .AND. ( LEN_TRIM( OpacityTableName_EmAb_Option ) > 1 ) )THEN
      OpacityTableName_EmAb = TRIM( OpacityTableName_EmAb_Option )
      Include_EmAb = .TRUE.
    ELSE
      OpacityTableName_EmAb = ''
      Include_EmAb = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_Iso_Option ) &
        .AND. ( LEN_TRIM( OpacityTableName_Iso_Option ) > 1 ) )THEN
      OpacityTableName_Iso = TRIM( OpacityTableName_Iso_Option )
      Include_Iso = .TRUE.
    ELSE
      OpacityTableName_Iso = ''
      Include_Iso = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_NES_Option ) &
        .AND. ( LEN_TRIM( OpacityTableName_NES_Option ) > 1 ) )THEN
      OpacityTableName_NES = TRIM( OpacityTableName_NES_Option )
      Include_NES = .TRUE.
    ELSE
      OpacityTableName_NES = ''
      Include_NES = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_Pair_Option ) &
        .AND. ( LEN_TRIM( OpacityTableName_Pair_Option ) > 1 ) )THEN
      OpacityTableName_Pair = TRIM( OpacityTableName_Pair_Option )
      Include_Pair = .TRUE.
    ELSE
      OpacityTableName_Pair = ''
      Include_Pair = .FALSE.
    END IF

    IF( PRESENT( OpacityTableName_Brem_Option ) &
        .AND. ( LEN_TRIM( OpacityTableName_Brem_Option ) > 1 ) )THEN
      OpacityTableName_Brem = TRIM( OpacityTableName_Brem_Option )
      Include_Brem = .TRUE.
    ELSE
      OpacityTableName_Brem = ''
      Include_Brem = .FALSE.
    END IF

    IF( PRESENT( EquationOfStateTableName_Option ) &
        .AND. ( LEN_TRIM( EquationOfStateTableName_Option ) > 1 ) )THEN
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

    Use_OpacityTables = .FALSE.

#ifdef MICROPHYSICS_WEAKLIB

    IF (   Include_EmAb &
      .OR. Include_Iso  &
      .OR. Include_NES  &
      .OR. Include_Pair &
      .OR. Include_Brem ) THEN
      Use_OpacityTables = .TRUE.
    ELSE
      Use_OpacityTables = .FALSE.
      RETURN
    END IF

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

    ! --- Thermodynamic State Indices ---

    iD_T = OPACITIES % TS % Indices % iRho
    iT_T = OPACITIES % TS % Indices % iT
    iY_T = OPACITIES % TS % Indices % iYe

    ! --- Opacity Density Cutoffs ---

    ! --- EmAb ---

    IF( PRESENT( EmAb_Nucleon_MinD_Option ) )THEN
      EmAb_Nucleon_MinD = EmAb_Nucleon_MinD_Option
    ELSE
      EmAb_Nucleon_MinD = OPACITIES % EOSTable % TS % minValues( iD_T )
    END IF

    IF( PRESENT( EmAb_Nucleon_MaxD_Option ) )THEN
      EmAb_Nucleon_MaxD = EmAb_Nucleon_MaxD_Option
    ELSE
      EmAb_Nucleon_MaxD = OPACITIES % EOSTable % TS % maxValues( iD_T )
    END IF

    IF( PRESENT( EmAb_Nuclei_MinD_Option ) )THEN
      EmAb_Nuclei_MinD = EmAb_Nuclei_MinD_Option
    ELSE
      EmAb_Nuclei_MinD = OPACITIES % EOSTable % TS % minValues( iD_T )
    END IF

    IF( PRESENT( EmAb_Nuclei_MaxD_Option ) )THEN
      EmAb_Nuclei_MaxD = EmAb_Nuclei_MaxD_Option
    ELSE
      EmAb_Nuclei_MaxD = OPACITIES % EOSTable % TS % maxValues( iD_T )
    END IF

    IF( PRESENT( EmAb_MinD_Option ) )THEN
      EmAb_MinD = EmAb_MinD_Option
    ELSE
      EmAb_MinD = MIN( EmAb_Nucleon_MinD, EmAb_Nuclei_MinD )
    END IF

    IF( PRESENT( EmAb_MaxD_Option ) )THEN
      EmAb_MaxD = EmAb_MaxD_Option
    ELSE
      EmAb_MaxD = MAX( EmAb_Nucleon_MaxD, EmAb_Nuclei_MaxD )
    END IF

    ! --- Iso ---

    IF( PRESENT( Iso_MinD_Option ) )THEN
      Iso_MinD = Iso_MinD_Option
    ELSE
      Iso_MinD = OPACITIES % EOSTable % TS % minValues( iD_T )
    END IF

    IF( PRESENT( Iso_MaxD_Option ) )THEN
      Iso_MaxD = Iso_MaxD_Option
    ELSE
      Iso_MaxD = OPACITIES % EOSTable % TS % maxValues( iD_T )
    END IF

    ! --- NES ---

    IF( PRESENT( NES_MinD_Option ) )THEN
      NES_MinD = NES_MinD_Option
    ELSE
      NES_MinD = OPACITIES % EOSTable % TS % minValues( iD_T )
    END IF

    IF( PRESENT( NES_MaxD_Option ) )THEN
      NES_MaxD = NES_MaxD_Option
    ELSE
      NES_MaxD = OPACITIES % EOSTable % TS % maxValues( iD_T )
    END IF

    ! --- Pair ---

    IF( PRESENT( Pair_MinD_Option ) )THEN
      Pair_MinD = Pair_MinD_Option
    ELSE
      Pair_MinD = OPACITIES % EOSTable % TS % minValues( iD_T )
    END IF

    IF( PRESENT( Pair_MaxD_Option ) )THEN
      Pair_MaxD = Pair_MaxD_Option
    ELSE
      Pair_MaxD = OPACITIES % EOSTable % TS % maxValues( iD_T )
    END IF

    ! --- Brem ---

    IF( PRESENT( Brem_MinD_Option ) )THEN
      Brem_MinD = Brem_MinD_Option
    ELSE
      Brem_MinD = OPACITIES % EOSTable % TS % minValues( iD_T )
    END IF

    IF( PRESENT( Brem_MaxD_Option ) )THEN
      Brem_MaxD = Brem_MaxD_Option
    ELSE
      Brem_MaxD = OPACITIES % EOSTable % TS % maxValues( iD_T )
    END IF

    ! --- NNS ---

    IF( PRESENT( NNS_MinD_Option ) )THEN
      NNS_MinD = NNS_MinD_Option
    ELSE
      NNS_MinD = OPACITIES % EOSTable % TS % minValues( iD_T )
    END IF

    IF( PRESENT( NNS_MaxD_Option ) )THEN
      NNS_MaxD = NNS_MaxD_Option
    ELSE
      NNS_MaxD = OPACITIES % EOSTable % TS % maxValues( iD_T )
    END IF

    ! --- NuPair ---

    IF( PRESENT( NuPair_MinD_Option ) )THEN
      NuPair_MinD = NuPair_MinD_Option
    ELSE
      NuPair_MinD = OPACITIES % EOSTable % TS % minValues( iD_T )
    END IF

    IF( PRESENT( NuPair_MaxD_Option ) )THEN
      NuPair_MaxD = NuPair_MaxD_Option
    ELSE
      NuPair_MaxD = OPACITIES % EOSTable % TS % maxValues( iD_T )
    END IF

    ! --- Cutoff For All Opacities ---

    IF( PRESENT( Op_MinD_Option ) )THEN
      Op_MinD = Op_MinD_Option
    ELSE
      !Op_MinD = MIN( EmAb_MinD, Iso_MinD, NES_MinD, Pair_MinD, Brem_MinD, NNS_MinD, NuPair_MinD )
      Op_MinD = MIN( EmAb_MinD, Iso_MinD, NES_MinD, Pair_MinD, Brem_MinD )
    END IF

    IF( PRESENT( Op_MaxD_Option ) )THEN
      Op_MaxD = Op_MaxD_Option
    ELSE
      !Op_MaxD = MAX( EmAb_MaxD, Iso_MaxD, NES_MaxD, Pair_MaxD, Brem_MaxD, NNS_MaxD, NuPair_MaxD )
      Op_MaxD = MAX( EmAb_MaxD, Iso_MaxD, NES_MaxD, Pair_MaxD, Brem_MaxD )
    END IF

    ! --- Make Cutoffs Consistent ---

    EmAb_MinD   = MAX( EmAb_MinD  , Op_MinD )
    Iso_MinD    = MAX( Iso_MinD   , Op_MinD )
    NES_MinD    = MAX( NES_MinD   , Op_MinD )
    Pair_MinD   = MAX( Pair_MinD  , Op_MinD )
    Brem_MinD   = MAX( Brem_MinD  , Op_MinD )
    NNS_MinD    = MAX( NNS_MinD   , Op_MinD )
    NuPair_MinD = MAX( NuPair_MinD, Op_MinD )

    EmAb_Nucleon_MinD = MAX( EmAb_Nucleon_MinD, EmAb_MinD )
    EmAb_Nuclei_MinD  = MAX( EmAb_Nuclei_MinD , EmAb_MinD )

    EmAb_MaxD   = MIN( EmAb_MaxD  , Op_MaxD )
    Iso_MaxD    = MIN( Iso_MaxD   , Op_MaxD )
    NES_MaxD    = MIN( NES_MaxD   , Op_MaxD )
    Pair_MaxD   = MIN( Pair_MaxD  , Op_MaxD )
    Brem_MaxD   = MIN( Brem_MaxD  , Op_MaxD )
    NNS_MaxD    = MIN( NNS_MaxD   , Op_MaxD )
    NuPair_MaxD = MIN( NuPair_MaxD, Op_MaxD )

    EmAb_Nucleon_MaxD = MIN( EmAb_Nucleon_MaxD, EmAb_MaxD )
    EmAb_Nuclei_MaxD  = MIN( EmAb_Nuclei_MaxD , EmAb_MaxD )

    nPointsE = nE * nNodesE

    nOpacities_NES  = OPACITIES % Scat_NES % nOpacities
    nMoments_NES    = OPACITIES % Scat_NES % nMoments
    nPointsT_NES    = OPACITIES % Scat_NES % nPoints(4)
    nPointsEta_NES  = OPACITIES % Scat_NES % nPoints(5)

    nOpacities_Pair = OPACITIES % Scat_Pair % nOpacities
    nMoments_Pair   = OPACITIES % Scat_Pair % nMoments
    nPointsT_Pair   = OPACITIES % Scat_Pair % nPoints(4)
    nPointsEta_Pair = OPACITIES % Scat_Pair % nPoints(5)

    nOpacities_Brem = OPACITIES % Scat_Brem % nOpacities
    nMoments_Brem   = OPACITIES % Scat_Brem % nMoments
    nPointsD_Brem   = OPACITIES % Scat_Brem % nPoints(4)
    nPointsT_Brem   = OPACITIES % Scat_Brem % nPoints(5)

    C1 = [ C1_NuE, C1_NuE_Bar, &
           C1_NuM, C1_NuM_Bar, &
           C1_NuT, C1_NuT_Bar ]

    C2 = [ C2_NuE, C2_NuE_Bar, &
           C2_NuM, C2_NuM_Bar, &
           C2_NuT, C2_NuT_Bar ]

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

    CALL DeAllocateOpacityTable( OPACITIES )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( always, to: LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$OMP                  OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem, &
    !$OMP                  EmAb_T, Iso_T, NES_T, Pair_T, Brem_T, &
    !$OMP                  NES_AT, Pair_AT, Brem_AT, C1, C2, &
    !$OMP                  EmAb_Nucleon_MinD, EmAb_Nucleon_MaxD, &
    !$OMP                  EmAb_Nuclei_MinD, EmAb_Nuclei_MaxD, &
    !$OMP                  EmAb_MinD, EmAb_MaxD, &
    !$OMP                  Iso_MinD, Iso_MaxD, &
    !$OMP                  NES_MinD, NES_MaxD, &
    !$OMP                  Pair_MinD, Pair_MaxD, &
    !$OMP                  Brem_MinD, Brem_MaxD, &
    !$OMP                  NNS_MinD, NNS_MaxD, &
    !$OMP                  NuPair_MinD, NuPair_MaxD, &
    !$OMP                  Op_MinD, Op_MaxD )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
    !$ACC   OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem, &
    !$ACC   EmAb_T, Iso_T, NES_T, Pair_T, Brem_T, &
    !$ACC   NES_AT, Pair_AT, Brem_AT, C1, C2, &
    !$ACC   EmAb_Nucleon_MinD, EmAb_Nucleon_MaxD, &
    !$ACC   EmAb_Nuclei_MinD, EmAb_Nuclei_MaxD, &
    !$ACC   EmAb_MinD, EmAb_MaxD, &
    !$ACC   Iso_MinD, Iso_MaxD, &
    !$ACC   NES_MinD, NES_MaxD, &
    !$ACC   Pair_MinD, Pair_MaxD, &
    !$ACC   Brem_MinD, Brem_MaxD, &
    !$ACC   NNS_MinD, NNS_MaxD, &
    !$ACC   NuPair_MinD, NuPair_MaxD, &
    !$ACC   Op_MinD, Op_MaxD )
#endif

    ASSOCIATE ( CenterE => MeshE % Center, &
                WidthE  => MeshE % Width, &
                NodesE  => MeshE % Nodes )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP MAP( to: CenterE, WidthE, NodesE ) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC COPYIN( CenterE, WidthE, NodesE ) &
    !$ACC PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 ) &
    !$ACC PRESENT( LogEs_T, OS_NES, NES_T, NES_AT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#endif
    DO iS = 1, nOpacities_NES
      DO iM = 1, nMoments_NES
        DO iEta = 1, nPointsEta_NES
          DO iT = 1, nPointsT_NES
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP MAP( to: CenterE, WidthE, NodesE ) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC COPYIN( CenterE, WidthE, NodesE ) &
    !$ACC PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 ) &
    !$ACC PRESENT( LogEs_T, OS_Pair, Pair_T, Pair_AT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#endif
    DO iS = 1, nOpacities_Pair
      DO iM = 1, nMoments_Pair
        DO iEta = 1, nPointsEta_Pair
          DO iT = 1, nPointsT_Pair
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

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP MAP( to: CenterE, WidthE, NodesE ) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC COPYIN( CenterE, WidthE, NodesE ) &
    !$ACC PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 ) &
    !$ACC PRESENT( LogEs_T, OS_Brem, Brem_T, Brem_AT )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO COLLAPSE(6) &
    !$OMP PRIVATE( LogE1, LogE2, iE1, iE2, iNodeE1, iNodeE2 )
#endif
    DO iS = 1, nOpacities_Brem
      DO iM = 1, nMoments_Brem
        DO iD = 1, nPointsD_Brem
          DO iT = 1, nPointsT_Brem
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

    END ASSOCIATE

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

  IF ( Use_OpacityTables ) THEN

#if defined(THORNADO_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
      !$OMP               OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem, &
      !$OMP               EmAb_T, Iso_T, NES_T, Pair_T, Brem_T, &
      !$OMP               NES_AT, Pair_AT, Brem_AT, C1, C2, &
      !$OMP               EmAb_Nucleon_MinD, EmAb_Nucleon_MaxD, &
      !$OMP               EmAb_Nuclei_MinD, EmAb_Nuclei_MaxD, &
      !$OMP               EmAb_MinD, EmAb_MaxD, &
      !$OMP               Iso_MinD, Iso_MaxD, &
      !$OMP               NES_MinD, NES_MaxD, &
      !$OMP               Pair_MinD, Pair_MaxD, &
      !$OMP               Brem_MinD, Brem_MaxD, &
      !$OMP               NNS_MinD, NNS_MaxD, &
      !$OMP               NuPair_MinD, NuPair_MaxD, &
      !$OMP               Op_MinD, Op_MaxD )
#endif

      DEALLOCATE( Es_T, Ds_T, Ts_T, Ys_T, Etas_T )
      DEALLOCATE( LogEs_T, LogDs_T, LogTs_T, LogEtas_T )

      DEALLOCATE( OS_EmAb, OS_Iso, OS_NES, OS_Pair, OS_Brem )
      DEALLOCATE( EmAb_T, Iso_T, NES_T, Pair_T, Brem_T )
      DEALLOCATE( NES_AT, Pair_AT, Brem_AT )

  END IF

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


  LOGICAL FUNCTION QueryOpacity_EmAb_Nuclei( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_EmAb_Nuclei &
        = ( D >= EmAb_Nuclei_MinD .AND. D <= EmAb_Nuclei_MaxD )

  END FUNCTION QueryOpacity_EmAb_Nuclei


  LOGICAL FUNCTION QueryOpacity_EmAb_Nucleon( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_EmAb_Nucleon &
        = ( D >= EmAb_Nucleon_MinD .AND. D <= EmAb_Nucleon_MaxD )

  END FUNCTION QueryOpacity_EmAb_Nucleon


  LOGICAL FUNCTION QueryOpacity_EmAb( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_EmAb &
        = ( D >= EmAb_MinD .AND. D <= EmAb_MaxD ) &
            .AND. (      QueryOpacity_EmAb_Nucleon( D ) &
                    .OR. QueryOpacity_EmAb_Nuclei ( D ) )

  END FUNCTION QueryOpacity_EmAb


  LOGICAL FUNCTION QueryOpacity_Iso( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_Iso &
        = ( D >= Iso_MinD .AND. D <= Iso_MaxD )

  END FUNCTION QueryOpacity_Iso


  LOGICAL FUNCTION QueryOpacity_NES( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_NES &
        = ( D >= NES_MinD .AND. D <= NES_MaxD )

  END FUNCTION QueryOpacity_NES


  LOGICAL FUNCTION QueryOpacity_Pair( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_Pair &
        = ( D >= Pair_MinD .AND. D <= Pair_MaxD )

  END FUNCTION QueryOpacity_Pair


  LOGICAL FUNCTION QueryOpacity_Brem( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_Brem &
        = ( D >= Brem_MinD .AND. D <= Brem_MaxD )

  END FUNCTION QueryOpacity_Brem


  LOGICAL FUNCTION QueryOpacity_NNS( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_NNS &
        = ( D >= NNS_MinD .AND. D <= NNS_MaxD )

  END FUNCTION QueryOpacity_NNS


  LOGICAL FUNCTION QueryOpacity_NuPair( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity_NuPair &
        = ( D >= NuPair_MinD .AND. D <= NuPair_MaxD )

  END FUNCTION QueryOpacity_NuPair


  LOGICAL FUNCTION QueryOpacity( D )
#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D

    QueryOpacity &
        = ( D >= Op_MinD .AND. D <= Op_MaxD ) &
            .AND. (      QueryOpacity_EmAb( D ) &
                    .OR. QueryOpacity_Iso ( D ) &
                    .OR. QueryOpacity_NES ( D ) &
                    .OR. QueryOpacity_Pair( D ) &
                    .OR. QueryOpacity_Brem( D ) )

  END FUNCTION QueryOpacity


END MODULE OpacityModule_TABLE
