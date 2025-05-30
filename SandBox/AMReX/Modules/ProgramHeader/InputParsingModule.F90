MODULE InputParsingModule

  ! --- AMReX modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    InitializeProgramHeader, &
    bcZ, &
    nDimsX
  USE UnitsModule, ONLY: &
    ActivateUnitsDisplay, &
    UnitsDisplay
  USE GeometryFieldsModule, ONLY: &
    CoordinateSystem
  USE RadiationFieldsModule, ONLY: &
    SetNumberOfSpecies

  ! --- Local modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF

  IMPLICIT NONE

  ! -- thornado ---

  CHARACTER(:), ALLOCATABLE :: ProgramName
  INTEGER     , ALLOCATABLE :: swX(:)
  INTEGER     , ALLOCATABLE :: bcX(:)
  INTEGER                   :: nNodes
  REAL(DP)                  :: t_end, t_wrt, t_chk, dt_wrt, dt_chk, dt_rel
  INTEGER                   :: iCycleW, iCycleChk, iCycleD, iRestart, iReGrid
  LOGICAL                   :: RwChkFields_uGF, RwChkFields_uCF, &
                               RwChkFields_uDF, RwChkFields_uCR
  LOGICAL                   :: UsePhysicalUnits
  LOGICAL                   :: DEBUG
  LOGICAL                   :: SolveGravity_NR

  ! --- Transport ---

  INTEGER, ALLOCATABLE :: bcZ_TwoMoment(:)
  INTEGER              :: nE, nSpecies, swE, bcE
  REAL(DP)             :: eL, eR, zoomE

  ! --- Opacity Tables ---

  CHARACTER(:), ALLOCATABLE :: OpacityTableName_AbEm
  CHARACTER(:), ALLOCATABLE :: OpacityTableName_Iso
  CHARACTER(:), ALLOCATABLE :: OpacityTableName_NES
  CHARACTER(:), ALLOCATABLE :: OpacityTableName_Pair
  CHARACTER(:), ALLOCATABLE :: OpacityTableName_Brem

  ! --- Non-Linear Solver Parameters ---

  INTEGER  ::  M_outer
  INTEGER  ::  MaxIter_outer
  REAL(DP) ::  Rtol_outer
  INTEGER  ::  M_inner
  INTEGER  ::  MaxIter_inner
  REAL(DP) ::  Rtol_inner
  LOGICAL  ::  Include_NES
  LOGICAL  ::  Include_Pair
  LOGICAL  ::  Include_Brem
  LOGICAL  ::  Include_LinCorr
  REAL(DP), ALLOCATABLE ::  wMatterRHS(:)

  ! --- geometry ---

  INTEGER               :: coord_sys
  REAL(DP), ALLOCATABLE :: xL(:), xR(:)

  ! --- amr ---

  INTEGER :: MaxGridSizeX1
  INTEGER :: MaxGridSizeX2
  INTEGER :: MaxGridSizeX3
  INTEGER :: BlockingFactorX1
  INTEGER :: BlockingFactorX2
  INTEGER :: BlockingFactorX3
  INTEGER :: MaxGridSizeX(3)
  INTEGER :: BlockingFactor(3)
  INTEGER :: MaxLevel
  INTEGER :: nLevels
  INTEGER :: nMaxLevels
  LOGICAL :: UseTiling
  LOGICAL :: UseFluxCorrection_Euler
  LOGICAL :: UseFluxCorrection_TwoMoment
  LOGICAL :: UseAMR
  LOGICAL :: IsPeriodic(3)
  INTEGER     , ALLOCATABLE :: IsPeriodicInt(:)
  INTEGER     , ALLOCATABLE :: nX(:)
  INTEGER     , ALLOCATABLE :: RefinementRatio(:)
  INTEGER     , ALLOCATABLE :: StepNo(:)
  INTEGER     , ALLOCATABLE :: nRefinementBuffer(:)
  REAL(DP)    , ALLOCATABLE :: TagCriteria(:)
  CHARACTER(:), ALLOCATABLE :: RefinementScheme

  REAL(DP), ALLOCATABLE :: dt   (:), dt_TM(:)
  REAL(DP), ALLOCATABLE :: t_old(:)
  REAL(DP), ALLOCATABLE :: t_new(:)
  CHARACTER(:), ALLOCATABLE :: PlotFileNameRoot
  INTEGER :: iOS_CPP(3)

  LOGICAL                   :: WriteNodalData
  CHARACTER(:), ALLOCATABLE :: NodalDataFileName

CONTAINS


  SUBROUTINE InitializeParameters

    TYPE(amrex_parmparse) :: PP
    INTEGER               :: iDimX

#if defined( THORNADO_AMREX_GIT_HASH )

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)

      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git hash:   ', THORNADO_AMREX_GIT_HASH
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git date:   ', THORNADO_AMREX_GIT_DATE
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git branch: ', THORNADO_AMREX_GIT_BRANCH
      WRITE(*,'(2x,A,A)') &
        'INFO: thornado git url:    ', THORNADO_AMREX_GIT_URL

    END IF

#endif

    ! --- debug Parameters debug.* ---

    DEBUG             = .FALSE.
    WriteNodalData    = .FALSE.
    NodalDataFileName = ''
    CALL amrex_parmparse_build( PP, 'debug' )
      CALL PP % query( 'DEBUG', &
                        DEBUG )
      CALL PP % query( 'WriteNodalData', &
                        WriteNodalData )
      CALL PP % query( 'NodalDataFileName', &
                        NodalDataFileName )
    CALL amrex_parmparse_destroy( PP )

    ! --- thornado Parameters thornado.* ---

    UsePhysicalUnits = .FALSE.
    PlotFileNameRoot = 'plt'
    iCycleD          = 10
    iCycleW          = -1
    iCycleChk        = -1
    RwChkFields_uGF  = .TRUE.
    RwChkFields_uCF  = .TRUE.
    RwChkFields_uDF  = .FALSE.
    RwChkFields_uCR  = .FALSE.
    iRestart         = -1
    dt_wrt           = -1.0_DP
    dt_chk           = -1.0_DP
    dt_rel           = 0.0_DP
    iReGrid          = 1
    SolveGravity_NR  = .FALSE.
    nE               = 1
    nSpecies         = 1
    swE              = 0
    bcE              = 0
    bcZ_TwoMoment    = [ 0, 0, 0, 0 ]
    eL               = 0.0_DP
    eR               = 1.0_DP
    zoomE            = 1.0_DP
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get   ( 'ProgramName', &
                         ProgramName )
      CALL PP % get   ( 'nNodes', &
                         nNodes )
      CALL PP % getarr( 'swX', &
                         swX )
      CALL PP % getarr( 'bcX', &
                         bcX )
      CALL PP % queryarr( 'bcZ_TwoMoment', &
                           bcZ_TwoMoment )
      CALL PP % get   ( 't_end', &
                         t_end )
      CALL PP % query ( 'iCycleD', &
                         iCycleD )
      CALL PP % query ( 'PlotFileNameRoot', &
                         PlotFileNameRoot )
      CALL PP % query ( 'iCycleW', &
                         iCycleW )
      CALL PP % query ( 'iCycleChk', &
                         iCycleChk )
      CALL PP % query ( 'RwChkFields_uGF', &
                         RwChkFields_uGF )
      CALL PP % query ( 'RwChkFields_uCF', &
                         RwChkFields_uCF )
      CALL PP % query ( 'RwChkFields_uDF', &
                         RwChkFields_uDF )
      CALL PP % query ( 'RwChkFields_uCR', &
                         RwChkFields_uCR )
      CALL PP % query ( 'iRestart', &
                         iRestart )
      CALL PP % query ( 'dt_wrt', &
                         dt_wrt )
      CALL PP % query ( 'dt_chk', &
                         dt_chk )
      CALL PP % query ( 'dt_rel', &
                         dt_rel )
      CALL PP % query ( 'UsePhysicalUnits', &
                         UsePhysicalUnits )
      CALL PP % query ( 'iReGrid', &
                         iReGrid )
      CALL PP % query ( 'SolveGravity_NR', &
                         SolveGravity_NR )
      CALL PP % query ( 'nE', &
                         nE )
      CALL PP % query ( 'nSpecies', &
                         nSpecies )
      CALL PP % query ( 'swE', &
                         swE )
      CALL PP % query ( 'bcE', &
                         bcE )
      CALL PP % query ( 'eL', &
                         eL )
      CALL PP % query ( 'eR', &
                         eR )
      CALL PP % query ( 'zoomE', &
                         zoomE )
    CALL amrex_parmparse_destroy( PP )

    CALL SetNumberOfSpecies &
           ( nSpecies, Verbose_Option = amrex_parallel_ioprocessor() )

    IF( iCycleW * dt_wrt .GT. Zero ) &
      CALL DescribeError_MF &
             ( 101, Int_Option = [ iCycleW ], Real_Option = [ dt_wrt ] )

    IF( iCycleChk * dt_chk .GT. Zero ) &
      CALL DescribeError_MF &
             ( 102, Int_Option = [ iCycleChk ], Real_Option = [ dt_chk ] )

    ! --- Parameters geometry.* ---

    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'coord_sys', &
                         coord_sys )
      CALL PP % getarr( 'prob_lo', &
                         xL )
      CALL PP % getarr( 'prob_hi', &
                         xR )
      CALL PP % getarr( 'is_periodic', &
                        IsPeriodicInt )
    CALL amrex_parmparse_destroy( PP )

    IsPeriodic = .FALSE.

    DO iDimX = 1, amrex_spacedim

      IF( IsPeriodicInt(iDimX) .EQ. 1 ) IsPeriodic(iDimX) = .TRUE.

    END DO

    IF     ( coord_sys .EQ. 0 )THEN

      CoordinateSystem = 'CARTESIAN'

    ELSE IF( coord_sys .EQ. 1 )THEN

      CoordinateSystem = 'CYLINDRICAL'

    ELSE IF( coord_sys .EQ. 2 )THEN

      CoordinateSystem = 'SPHERICAL'

    ELSE

      CALL DescribeError_MF( 103, Int_Option = [ coord_sys ] )

    END IF

    IF( UsePhysicalUnits )THEN

      CALL ActivateUnitsDisplay &
             ( CoordinateSystem_Option = TRIM( CoordinateSystem ) )

      t_end  = t_end  * UnitsDisplay % TimeUnit
      dt_wrt = dt_wrt * UnitsDisplay % TimeUnit
      dt_chk = dt_chk * UnitsDisplay % TimeUnit
      dt_rel = dt_rel * UnitsDisplay % TimeUnit

      xL(1) = xL(1) * UnitsDisplay % LengthX1Unit
      xR(1) = xR(1) * UnitsDisplay % LengthX1Unit
      xL(2) = xL(2) * UnitsDisplay % LengthX2Unit
      xR(2) = xR(2) * UnitsDisplay % LengthX2Unit
      xL(3) = xL(3) * UnitsDisplay % LengthX3Unit
      xR(3) = xR(3) * UnitsDisplay % LengthX3Unit

      eL = eL * UnitsDisplay % EnergyUnit
      eR = eR * UnitsDisplay % EnergyUnit

    END IF

    ! --- Opacity table parameters OP.* ---

    OpacityTableName_AbEm = ''
    OpacityTableName_Iso  = ''
    OpacityTableName_NES  = ''
    OpacityTableName_Pair = ''
    OpacityTableName_Brem = ''
    CALL amrex_parmparse_build( PP, 'OP' )
      CALL PP % query( 'OpacityTableName_AbEm', &
                        OpacityTableName_AbEm )
      CALL PP % query( 'OpacityTableName_Iso', &
                        OpacityTableName_Iso )
      CALL PP % query( 'OpacityTableName_NES', &
                        OpacityTableName_NES )
      CALL PP % query( 'OpacityTableName_Pair', &
                        OpacityTableName_Pair )
      CALL PP % query( 'OpacityTableName_Brem', &
                        OpacityTableName_Brem )
    CALL amrex_parmparse_destroy( PP )

    ! --- Non-Linear Solver parameters NL.* ---

    M_outer         = 2
    MaxIter_outer   = 100
    Rtol_outer      = 1.0d-8
    M_inner         = 2
    MaxIter_inner   = 100
    Rtol_inner      = 1.0d-8
    Include_NES     = .FALSE.
    Include_Pair    = .FALSE.
    Include_Brem    = .FALSE.
    Include_LinCorr = .FALSE.
    wMatterRHS      = [ 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP ]
    CALL amrex_parmparse_build( PP, 'NL' )
      CALL PP % query( 'M_outer', &
                        M_outer )
      CALL PP % query( 'MaxIter_outer', &
                        MaxIter_outer )
      CALL PP % query( 'Rtol_outer', &
                        Rtol_outer )
      CALL PP % query( 'M_inner', &
                        M_inner )
      CALL PP % query( 'MaxIter_inner', &
                        MaxIter_inner )
      CALL PP % query( 'Rtol_inner', &
                        Rtol_inner )
      CALL PP % query( 'Include_NES', &
                        Include_NES )
      CALL PP % query( 'Include_Pair', &
                        Include_Pair )
      CALL PP % query( 'Include_Brem', &
                        Include_Brem )
      CALL PP % query( 'Include_LinCorr', &
                        Include_LinCorr )
      CALL PP % queryarr( 'wMatterRHS', &
                           wMatterRHS )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters amr.* ---

    IF( amrex_spacedim .EQ. 1 )THEN

      MaxGridSizeX1 = 128
      MaxGridSizeX2 = 1
      MaxGridSizeX3 = 1
      BlockingFactorX1 = 8
      BlockingFactorX2 = 1
      BlockingFactorX3 = 1

    ELSE IF( amrex_spacedim .EQ. 2 )THEN

      MaxGridSizeX1 = 128
      MaxGridSizeX2 = 128
      MaxGridSizeX3 = 1
      BlockingFactorX1 = 8
      BlockingFactorX2 = 8
      BlockingFactorX3 = 1

    ELSE

      MaxGridSizeX1 = 32
      MaxGridSizeX2 = 32
      MaxGridSizeX3 = 32
      BlockingFactorX1 = 8
      BlockingFactorX2 = 8
      BlockingFactorX3 = 8

    END IF
    UseAMR                      = .FALSE.
    UseFluxCorrection_Euler     = .FALSE.
    UseFluxCorrection_TwoMoment = .FALSE.
    UseTiling                   = .FALSE.
    RefinementScheme            = ''
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr  ( 'n_cell', &
                           nX )
      CALL PP % query   ( 'max_grid_size_x', &
                           MaxGridSizeX1 )
      CALL PP % query   ( 'max_grid_size_y', &
                           MaxGridSizeX2 )
      CALL PP % query   ( 'max_grid_size_z', &
                           MaxGridSizeX3 )
      CALL PP % query   ( 'blocking_factor_x', &
                           BlockingFactorX1 )
      CALL PP % query   ( 'blocking_factor_y', &
                           BlockingFactorX2 )
      CALL PP % query   ( 'blocking_factor_z', &
                           BlockingFactorX3 )
      CALL PP % get     ( 'max_level', &
                           MaxLevel )
      IF( MaxLevel .GT. 0 )THEN
        CALL PP % query ( 'UseAMR', &
                           UseAMR )
        CALL PP % query ( 'UseFluxCorrection_Euler', &
                           UseFluxCorrection_Euler )
        CALL PP % query ( 'UseFluxCorrection_TwoMoment', &
                           UseFluxCorrection_TwoMoment )
        CALL PP % query ( 'RefinementScheme', &
                           RefinementScheme )
        CALL PP % getarr( 'TagCriteria', &
                           TagCriteria )
        CALL PP % getarr( 'n_error_buf', &
                           nRefinementBuffer )
      END IF
      CALL PP % getarr  ( 'ref_ratio', &
                           RefinementRatio )
      CALL PP % query   ( 'UseTiling', &
                           UseTiling )
    CALL amrex_parmparse_destroy( PP )

    MaxGridSizeX   = [ MaxGridSizeX1   , MaxGridSizeX2   , MaxGridSizeX3    ]
    BlockingFactor = [ BlockingFactorX1, BlockingFactorX2, BlockingFactorX3 ]

    nMaxLevels = MaxLevel + 1
    nLevels    = nMaxLevels

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option      = nNodes, &
             nX_Option          = nX, &
             nE_Option          = nE, &
             swX_Option         = swX, &
             swE_Option         = swE, &
             xL_Option          = xL, &
             xR_Option          = xR, &
             eL_option          = eL, &
             eR_option          = eR, &
             bcX_Option         = bcX, &
             bcE_Option         = bcE, &
             Verbose_Option     = amrex_parallel_ioprocessor() )

    bcZ = bcZ_TwoMoment

    IF( nDimsX .NE. amrex_spacedim ) &
      CALL DescribeError_MF &
             ( 104, Int_Option = [ nDimsX, amrex_spacedim ] )

    iOS_CPP = 0

    iOS_CPP(1) = 1

    IF( amrex_spacedim .GT. 1 ) iOS_CPP(2) = 1
    IF( amrex_spacedim .GT. 2 ) iOS_CPP(3) = 1

  END SUBROUTINE InitializeParameters


  SUBROUTINE DescribeProgramHeader_AMReX

    CHARACTER(32) :: RFMT, IFMT, MFMT, TagUnits
    INTEGER :: iLevel, iDimX
    REAL(DP) :: MeshWidths(3,0:nMaxLevels-1)

    TagUnits = ''
    IF( TRIM( RefinementScheme ) .EQ. 'Mesh' )THEN
      TagUnits = TRIM( UnitsDisplay % LengthX1Label )
    ELSE IF( TRIM( RefinementScheme ) .EQ. 'Density' )THEN
      TagUnits = TRIM( UnitsDisplay % MassDensityLabel )
    END IF

    IF( .NOT. ALLOCATED( TagCriteria ) )THEN
      ALLOCATE( TagCriteria(nMaxLevels) )
      TagCriteria = 0.0_DP
    END IF

    IF( .NOT. ALLOCATED( nRefinementBuffer ) )THEN
      ALLOCATE( nRefinementBuffer(nMaxLevels) )
      nRefinementBuffer = 1
    END IF

    WRITE(RFMT,'(A,I2.2,A)') '(4x,A26,1x,', SIZE( TagCriteria ), 'ES11.3E3,x,A)'
    WRITE(IFMT,'(A,I2.2,A)') '(4x,A26,1x,', SIZE( nRefinementBuffer ), 'I3.2)'
    WRITE(MFMT,'(A,I2.2,A)') '(4x,A26,1x,', nMaxLevels, 'ES11.3E3,x,A)'

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,'(4x,A)')            'INFO: AMReX'
      WRITE(*,'(4x,A)')            '-----------'
      WRITE(*,'(4x,A26,1x,3I4.3)')   'MaxGridSizeX:', &
                                      MaxGridSizeX
      WRITE(*,'(4x,A26,1x,3I4.3)')   'BlockingFactor:', &
                                      BlockingFactor
      WRITE(*,'(4x,A26,1x,I2.2)')    'nMaxLevels:', &
                                      nMaxLevels
      WRITE(*,'(4x,A26,1x,I2.2)')    'iReGrid:', &
                                      iReGrid
      WRITE(*,'(4x,A26,1x,L)')       'UseFluxCorrection_Euler:', &
                                      UseFluxCorrection_Euler
      WRITE(*,'(4x,A26,1x,L)')       'UseTiling:', &
                                      UseTiling
      WRITE(*,'(4x,A26,1x,L)')       'UseAMR:', &
                                      UseAMR
      WRITE(*,TRIM(IFMT))            'nRefinementBuffer:', &
                                      nRefinementBuffer
      WRITE(*,TRIM(RFMT))            'TagCriteria:', &
                                      TagCriteria, TRIM( TagUnits )

      DO iLevel = 0, nMaxLevels-1

        DO iDimX = 1, nDimsX

          MeshWidths(iDimX,iLevel) &
            = ( ( xR(iDimX) - xL(iDimX) ) / nX(iDimX) ) / 2**(iLevel)

        END DO

      END DO

      WRITE(*,TRIM(MFMT)) 'MeshWidths (X1):', &
                           MeshWidths(1,:) / UnitsDisplay % LengthX1Unit, &
                                       TRIM( UnitsDisplay % LengthX1Label )
      IF( nDimsX .GT. 1 ) &
        WRITE(*,TRIM(MFMT)) 'MeshWidths (X2):', &
                             MeshWidths(2,:) / UnitsDisplay % LengthX2Unit, &
                                         TRIM( UnitsDisplay % LengthX2Label )
      IF( nDimsX .GT. 2 ) &
        WRITE(*,TRIM(MFMT)) 'MeshWidths (X3):', &
                             MeshWidths(3,:) / UnitsDisplay % LengthX3Unit, &
                                         TRIM( UnitsDisplay % LengthX3Label )
      WRITE(*,*)

    END IF

  END SUBROUTINE DescribeProgramHeader_AMReX


subroutine myamrfinalize
end subroutine myamrfinalize

END MODULE InputParsingModule
