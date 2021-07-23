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
  USE amrex_init_module, ONLY: &
    amrex_init, &
    amrex_initialized
  USE amrex_amr_module, ONLY: &
    amrex_amrcore_init, &
    amrex_amrcore_initialized

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    InitializeProgramHeader

  ! --- Local modules ---

  USE MF_KindModule, ONLY: &
    DP

  IMPLICIT NONE

  ! -- thornado ---

  CHARACTER(:), ALLOCATABLE :: ProgramName
  INTEGER, ALLOCATABLE  :: swX(:)
  INTEGER, ALLOCATABLE  :: bcX(:)
  INTEGER               :: nNodes
  REAL(DP)              :: Gamma_IDEAL

  ! --- geometry ---

  INTEGER               :: CoordSys
  REAL(DP), ALLOCATABLE :: xL(:), xR(:)

  ! --- amr ---

  INTEGER, ALLOCATABLE :: nX(:)
  INTEGER :: MaxGridSizeX1
  INTEGER :: MaxGridSizeX2
  INTEGER :: MaxGridSizeX3
  INTEGER :: BlockingFactorX1
  INTEGER :: BlockingFactorX2
  INTEGER :: BlockingFactorX3
  INTEGER :: MaxGridSizeX(3)
  INTEGER :: MaxLevel
  INTEGER :: nLevels
  LOGICAL :: UseTiling
  LOGICAL :: do_reflux
  INTEGER, ALLOCATABLE :: RefinementRatio(:)
  INTEGER, ALLOCATABLE :: StepNo(:)

  REAL(DP), ALLOCATABLE :: dt   (:)
  REAL(DP), ALLOCATABLE :: t_old(:)
  REAL(DP), ALLOCATABLE :: t_new(:)
  CHARACTER(:), ALLOCATABLE :: PlotFileBaseName
  INTEGER :: iOS_CPP(3)

CONTAINS


  SUBROUTINE InitializeParameters

    TYPE(amrex_parmparse) :: PP

    IF( .NOT. amrex_initialized() ) &
      CALL amrex_init()

    IF( .NOT. amrex_amrcore_initialized() ) &
      CALL amrex_amrcore_init()

    ! --- thornado paramaters thornado.* ---

    Gamma_IDEAL = 4.0_DP / 3.0_DP
    PlotFileBaseName = 'plt'
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get   ( 'ProgramName', ProgramName )
      CALL PP % get   ( 'nNodes', nNodes )
      CALL PP % getarr( 'swX', swX )
      CALL PP % getarr( 'bcX', bcX )
      CALL PP % query ( 'Gamma_IDEAL', Gamma_IDEAL )
      CALL PP % query ( 'PlotFileBaseName', PlotFileBaseName )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters geometry.* ---

    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'coord_sys', CoordSys )
      CALL PP % getarr( 'prob_lo'  , xL       )
      CALL PP % getarr( 'prob_hi'  , xR       )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters amr.* ---

    MaxGridSizeX1    = 1
    MaxGridSizeX2    = 1
    MaxGridSizeX3    = 1
    BlockingFactorX1 = 1
    BlockingFactorX2 = 1
    BlockingFactorX3 = 1
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr( 'n_cell'           , nX               )
      CALL PP % query ( 'max_grid_size_x'  , MaxGridSizeX1    )
      CALL PP % query ( 'max_grid_size_y'  , MaxGridSizeX2    )
      CALL PP % query ( 'max_grid_size_z'  , MaxGridSizeX3    )
      CALL PP % query ( 'blocking_factor_x', BlockingFactorX1 )
      CALL PP % query ( 'blocking_factor_y', BlockingFactorX2 )
      CALL PP % query ( 'blocking_factor_z', BlockingFactorX3 )
      CALL PP % get   ( 'max_level'        , MaxLevel         )
      CALL PP % get   ( 'UseTiling'        , UseTiling        )
      CALL PP % get   ( 'do_reflux'        , do_reflux        )
      CALL PP % getarr( 'ref_ratio'        , RefinementRatio  )
    CALL amrex_parmparse_destroy( PP )

    MaxGridSizeX = [ MaxGridSizeX1, MaxGridSizeX2, MaxGridSizeX3 ]
    nLevels = MaxLevel + 1

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option      = nNodes,              &
             nX_Option          = nX,                  &
             swX_Option         = swX,                 &
             xL_Option          = xL,                  &
             xR_Option          = xR,                  &
             bcX_Option         = bcX,                 &
             Verbose_Option     = amrex_parallel_ioprocessor() )

    iOS_CPP = 0

    iOS_CPP(1) = 1

    IF( amrex_spacedim .GT. 1 ) iOS_CPP(2) = 1
    IF( amrex_spacedim .GT. 2 ) iOS_CPP(3) = 1

  END SUBROUTINE InitializeParameters

END MODULE InputParsingModule
