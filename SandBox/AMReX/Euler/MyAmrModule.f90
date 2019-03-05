MODULE MyAmrModule

  USE iso_c_binding
  USE amrex_base_module, ONLY: &
    amrex_init, &
    amrex_initialized, &
    amrex_parallel_ioprocessor
  USE amrex_amr_module, ONLY: &
    amrex_amrcore_init, &
    amrex_amrcore_initialized, &
    amrex_is_all_periodic, &
    amrex_spacedim
  USE amrex_bc_types_module, ONLY: &
    amrex_bc_int_dir, &
    amrex_bc_foextrap
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_fort_module, ONLY: &
    amrex_real

  USE ProgramHeaderModule, ONLY: &
    InitializeProgramHeader, nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF, nPF, nAF
  USE GeometryFieldsModule, ONLY: &
    nGF

  USE MyAmrDataModule

  IMPLICIT NONE

  REAL(amrex_real)                    :: t_end, dt_wrt, Gamma_IDEAL, CFL
  INTEGER                             :: nNodes, nStages, nLevels, coord_sys
  INTEGER                             :: iCycleD, iCycleW, iCycleChk
  INTEGER,          ALLOCATABLE       :: MaxGridSize(:), nX(:), swX(:), bcX(:)
  REAL(amrex_real), ALLOCATABLE       :: xL(:), xR(:), dt(:), t(:)
  CHARACTER(LEN=:), ALLOCATABLE       :: ProgramName
  INTEGER,          ALLOCATABLE, SAVE :: StepNo(:)
  CHARACTER(LEN=32),             SAVE :: Coordsys

  ! --- Boundary Conditions ---
  INTEGER, ALLOCATABLE, PUBLIC, SAVE :: bcAMReX(:)

  ! --- Slope limiter ---
  LOGICAL          :: UseSlopeLimiter
  LOGICAL          :: UseCharacteristicLimiting
  LOGICAL          :: UseTroubledCellIndicator
  LOGICAL          :: UseAMReX
  REAL(amrex_real) :: SlopeTolerance
  REAL(amrex_real) :: BetaTVD, BetaTVB
  REAL(amrex_real) :: LimiterThresholdParameter


CONTAINS


  SUBROUTINE MyAmrInit

    TYPE(amrex_parmparse) :: PP

    IF( .NOT. amrex_initialized() ) &
      CALL amrex_init()

    IF( .NOT. amrex_amrcore_initialized() ) &
      CALL amrex_amrcore_init()

    ! --- thornado paramaters thornado.* ---
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get   ( 'dt_wrt',      dt_wrt )
      CALL PP % get   ( 't_end',       t_end )
      CALL PP % get   ( 'nNodes',      nNodes )
      CALL PP % get   ( 'nStages',     nStages )
      CALL PP % get   ( 'CFL',         CFL )
      CALL PP % get   ( 'ProgramName', ProgramName )
      CALL PP % get   ( 'Gamma',       Gamma_IDEAL )
      CALL PP % getarr( 'bcX',         bcX )
      CALL PP % getarr( 'swX',         swX )
      CALL PP % get   ( 'iCycleD',     iCycleD )
      CALL PP % get   ( 'iCycleW',     iCycleW )
      CALL PP % get   ( 'iCycleChk',   iCycleChk )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters geometry.* ---
    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'coord_sys',        coord_sys )
      CALL PP % getarr( 'prob_lo',          xL )
      CALL PP % getarr( 'prob_hi',          xR )
      CALL PP % getarr( 'bcAMReX',          bcAMReX )
    CALL amrex_parmparse_destroy( PP )
    IF     ( coord_sys .EQ. 0 )THEN
      CoordSys = 'CARTESIAN'
    ELSE IF( coord_sys .EQ. 1 )THEN
      CoordSys = 'CYLINDRICAL'
    ELSE IF( coord_sys .EQ. 2 )THEN
      CoordSys = 'SPHERICAL'
    ELSE
      STOP 'Invalid choice for coord_sys'
    END IF

    ! --- Parameters amr.*
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr( 'n_cell',        nX )
      CALL PP % getarr( 'max_grid_size', MaxGridSize )
      CALL PP % get   ( 'max_level',     nLevels )
    CALL amrex_parmparse_destroy( PP )

    ! --- Slope limiter parameters SL.*
    CALL amrex_parmparse_build( PP, 'SL' )
      CALL PP % get( 'UseSlopeLimiter',           UseSlopeLimiter )
      CALL PP % get( 'UseCharacteristicLimiting', UseCharacteristicLimiting )
      CALL PP % get( 'UseTroubledCellIndicator',  UseTroubledCellIndicator )
      CALL PP % get( 'SlopeTolerance',            SlopeTolerance )
      CALL PP % get( 'BetaTVD',                   BetaTVD )
      CALL PP % get( 'BetaTVB',                   BetaTVB )
      CALL PP % get( 'LimiterThresholdParameter', LimiterThresholdParameter )
    CALL amrex_parmparse_destroy( PP )

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX, &
             xL_Option = xL, xR_Option = xR, bcX_Option = bcX, &
             Verbose_Option = amrex_parallel_ioprocessor() )

    ALLOCATE( StepNo(0:nLevels) )
    StepNo = 0

    ALLOCATE( dt(0:nLevels) )
    dt = 1.0e-4_amrex_real

    ALLOCATE( t(0:nLevels) )
    t = 0.0e0_amrex_real

    CALL InitializeDataAMReX

  END SUBROUTINE MyAmrInit


  SUBROUTINE MyAmrFinalize

    CALL FinalizeDataAMReX

    DEALLOCATE( t )
    DEALLOCATE( dt )
    DEALLOCATE( StepNo )

  END SUBROUTINE MyAmrFinalize


END MODULE MyAmrModule
