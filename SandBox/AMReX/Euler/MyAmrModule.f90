MODULE MyAmrModule

  USE iso_c_binding
  USE amrex_amr_module
  USE amrex_fort_module, ONLY: &
    amrex_real

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF, nPF, nAF
  USE GeometryFieldsModule, ONLY: &
    nGF

  USE MyAmrDataModule

  IMPLICIT NONE

  REAL(amrex_real)                    :: t_end, dt_wrt, Gamma_IDEAL, CFL
  INTEGER                             :: nNodes, nStages, nLevels
  INTEGER                             :: iCycleD, iCycleW, iCycleChk
  INTEGER,          ALLOCATABLE       :: MaxGridSize(:), nX(:), swX(:), bcX(:)
  REAL(amrex_real), ALLOCATABLE       :: xL(:), xR(:), dt(:), t(:)
  CHARACTER(LEN=:), ALLOCATABLE       :: ProgramName, CoordSys
  INTEGER,          ALLOCATABLE, SAVE :: StepNo(:), nSubSteps(:)
  CHARACTER(LEN=:), ALLOCATABLE, SAVE :: Restart

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
    INTEGER               :: iLevel

    IF( .NOT. amrex_amrcore_initialized() ) &
      CALL amrex_amrcore_init()

!!$    CALL amrex_init_virtual_functions ( MyMakeNewLevelFromScratch,     &
!!$                                        my_make_new_level_from_coarse, &
!!$                                        my_remake_level,               &
!!$                                        MyClearLevel,                  &
!!$                                        my_error_estimate )

    ALLOCATE( CHARACTER(LEN=0) :: Restart )

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
      CALL PP % get   ( 'CoordinateSystem', CoordSys )
      CALL PP % getarr( 'prob_lo',          xL )
      CALL PP % getarr( 'prob_hi',          xR )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters amr.*
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr( 'n_cell',      nX )
      CALL PP % getarr( 'MaxGridSize', MaxGridSize )
      CALL PP % get   ( 'max_level',   nLevels )
      CALL PP % query ( 'Restart',     Restart )
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

!!$    if (.not. amrex_is_all_periodic()) then
!!$       lo_bc = amrex_bc_foextrap
!!$       hi_bc = amrex_bc_foextrap
!!$    end if

    ALLOCATE( StepNo(0:nLevels) )
    StepNo = 0

    ALLOCATE( nSubSteps(0:nLevels) )
    nSubSteps(0) = 1
    DO iLevel = 1, nLevels
      nSubSteps(iLevel) = amrex_ref_ratio(iLevel-1)
    END DO

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
    DEALLOCATE( nSubSteps )
    DEALLOCATE( StepNo )

  END SUBROUTINE MyAmrFinalize


END MODULE MyAmrModule
