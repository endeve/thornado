MODULE InputParsingModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,      ONLY: &
    AR => amrex_real, &
    amrex_spacedim
  USE amrex_base_module,      ONLY: &
    amrex_init, &
    amrex_initialized, &
    amrex_parallel_ioprocessor
  USE amrex_amr_module,       ONLY: &
    amrex_amrcore_init, &
    amrex_amrcore_initialized, &
    amrex_is_all_periodic, &
    amrex_spacedim
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_boxarray_module,  ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap
  USE amrex_geometry_module,  ONLY: &
    amrex_geometry

  ! --- thornado Modules ---
  USE ProgramHeaderModule,   ONLY: &
    InitializeProgramHeader, nDimsX
  USE UnitsModule,            ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay, &
    Centimeter, &
    Kilometer, &
    MeV, &
    SolarMass
  ! --- Local Modules ---
  USE MF_FieldsModule, ONLY: &
    InitializeDataAMReX, &
    FinalizeDataAMReX
  ! --- thornado ---
  REAL(AR)                       :: t_end, t_wrt, dt_wrt, t_chk, dt_chk, dt_rel
  REAL(AR),          ALLOCATABLE :: t(:), dt(:)
  REAL(AR),          ALLOCATABLE :: V_0(:)
  REAL(AR)                       :: CFL
  INTEGER                        :: nNodes, nStages
  INTEGER                        :: nE, swE, bcE
  INTEGER                        :: iCycleD, iCycleW, iCycleChk, iRestart, nSpecies
  INTEGER,           ALLOCATABLE :: nX(:), swX(:), bcX(:)
  REAL(AR),          ALLOCATABLE :: xL(:), xR(:)
  REAL(AR)                       :: eL, eR, zoomE
  REAL(AR)                       :: D_0, Chi, Sigma
  CHARACTER(LEN=:),  ALLOCATABLE :: ProgramName
  CHARACTER(LEN=:),  ALLOCATABLE :: Direction
  CHARACTER(LEN=:),  ALLOCATABLE :: Scheme
  CHARACTER(LEN=32), SAVE        :: CoordSys
  LOGICAL,           SAVE        :: UsePhysicalUnits
  LOGICAL,           SAVE        :: DEBUG



  ! --- AMReX  ---
  INTEGER                                    :: MaxLevel, nLevels, coord_sys
  INTEGER                                    :: MaxGridSizeX1
  INTEGER                                    :: MaxGridSizeX2
  INTEGER                                    :: MaxGridSizeX3
  INTEGER                                    :: BlockingFactorX1
  INTEGER                                    :: BlockingFactorX2
  INTEGER                                    :: BlockingFactorX3
  INTEGER                                    :: MaxGridSizeX(3)
  INTEGER,               ALLOCATABLE, SAVE   :: StepNo(:)
  TYPE(amrex_boxarray),  ALLOCATABLE, PUBLIC :: BA(:)
  TYPE(amrex_distromap), ALLOCATABLE, PUBLIC :: DM(:)
  TYPE(amrex_geometry),  ALLOCATABLE, PUBLIC :: GEOM(:)
  LOGICAL                                    :: UseTiling

  ! --- Equation Of State ---
  REAL(AR)                      :: Gamma_IDEAL
  CHARACTER(LEN=:), ALLOCATABLE :: EquationOfState
  CHARACTER(LEN=:), ALLOCATABLE :: EosTableName

  CHARACTER(LEN=:), ALLOCATABLE :: OpacityTableName_AbEm
  CHARACTER(LEN=:), ALLOCATABLE :: OpacityTableName_Iso
  CHARACTER(LEN=:), ALLOCATABLE :: OpacityTableName_NES
  CHARACTER(LEN=:), ALLOCATABLE :: OpacityTableName_Pair
  CHARACTER(LEN=:), ALLOCATABLE :: OpacityTableName_Brem

  ! --- Positivity limiter ---
  LOGICAL  :: UsePositivityLimiter_TwoMoment
  REAL(AR) :: Min_1_TwoMoment, Min_2_TwoMoment

  LOGICAL  :: UseSlopeLimiter_TwoMoment
  REAL(AR) :: BetaTVD_TwoMoment


  REAL(AR) :: Mass, R0, kT, mu0, E0

CONTAINS

  SUBROUTINE MyAmrInit

    REAL(AR), PARAMETER :: Zero = 0.0_AR
    REAL(AR), PARAMETER :: One  = 1.0_AR
    REAL(AR), PARAMETER :: Two  = 2.0_AR
    REAL(AR), PARAMETER :: Pi     = ACOS( -1.0_AR )
    REAL(AR), PARAMETER :: TwoPi  = 2.0_AR * Pi

    TYPE(amrex_parmparse) :: PP

    IF( .NOT. amrex_initialized() ) &
      CALL amrex_init()

    IF( .NOT. amrex_amrcore_initialized() ) &
      CALL amrex_amrcore_init()
    DEBUG = .FALSE.
    CALL amrex_parmparse_build( PP )
      CALL PP % query( 'DEBUG', DEBUG )
    CALL amrex_parmparse_destroy( PP )

    UsePhysicalUnits = .FALSE.
    Chi = 0.0_AR
    Sigma = 0.0_AR
    D_0 = 0.0_AR
    dt_rel = 0.0_AR

    Direction = ' '

    ! --- thornado paramaters thornado.* ---
    CALL amrex_parmparse_build( PP, 'thornado' )
      CALL PP % get   ( 'dt_wrt',           dt_wrt )
      CALL PP % get   ( 'dt_chk',           dt_chk )
      CALL PP % query ( 'dt_rel',           dt_rel )
      CALL PP % get   ( 't_end',            t_end )
      CALL PP % get   ( 'nNodes',           nNodes )
      CALL PP % get   ( 'CFL',              CFL )
      CALL PP % get   ( 'ProgramName',      ProgramName )
      CALL PP % get   ( 'Scheme',           Scheme )
      CALL PP % getarr( 'bcX',              bcX )
      CALL PP % getarr( 'swX',              swX )
      CALL PP % getarr( 'V_0',              V_0 )
      CALL PP % query ( 'nE', nE )
      CALL PP % get   ( 'swE',              swE )
      CALL PP % get   ( 'bcE',              bcE )
      CALL PP % get   ( 'eL',  eL )
      CALL PP % get   ( 'eR',  eR )
      CALL PP % query   ( 'D_0',  D_0 )
      CALL PP % query   ( 'Chi',  Chi )
      CALL PP % query   ( 'Sigma',  Sigma )
      CALL PP % get   ( 'zoomE',  zoomE )
      CALL PP % get   ( 'nSpecies',        nSpecies )
      CALL PP % get   ( 'iCycleD',          iCycleD )
      CALL PP % get   ( 'iCycleW',          iCycleW )
      CALL PP % get   ( 'iCycleChk',        iCycleChk )
      CALL PP % get   ( 'iRestart',         iRestart )
      CALL PP % query ( 'UsePhysicalUnits', UsePhysicalUnits )
      CALL PP % query ( 'Direction', Direction )

    CALL amrex_parmparse_destroy( PP )


    CFL = CFL / ( DBLE( amrex_spacedim ) * ( Two * DBLE( nNodes ) - One ) )

    ! --- Parameters geometry.* ---
    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'coord_sys',  coord_sys )
      CALL PP % getarr( 'prob_lo',    xL )
      CALL PP % getarr( 'prob_hi',    xR )
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

    Mass = 0.0_AR
    R0 = 0.0_AR
    E0 = 0.0_AR
    mu0 = 0.0_AR
    kT = 0.0_AR
    CALL amrex_parmparse_build( PP, 'ST' )
      CALL PP % query( 'Mass', Mass )
      CALL PP % query( 'R0'               ,R0 )
      CALL PP % query( 'mu0'               ,mu0 )
      CALL PP % query( 'E0'               ,E0 )
      CALL PP % query( 'kT'               ,kT )
    CALL amrex_parmparse_destroy( PP )

    IF( UsePhysicalUnits )THEN

      CALL ActivateUnitsDisplay( CoordinateSystem_Option = TRIM( CoordSys ) )

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

      Chi = Chi * ( 1.0_AR / Centimeter )

      Mass = Mass * SolarMass
      E0 = E0 * MeV
      mu0 = mu0 * MeV
      kT = kT * MeV
      R0 = R0 * kilometer

    END IF



    ! --- Parameters amr.* ---
    MaxGridSizeX1    = 1
    MaxGridSizeX2    = 1
    MaxGridSizeX3    = 1
    BlockingFactorX1 = 1
    BlockingFactorX2 = 1
    BlockingFactorX3 = 1
    UseTiling        = .TRUE.
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr( 'n_cell',            nX )
      CALL PP % query ( 'max_grid_size_x',   MaxGridSizeX1 )
      CALL PP % query ( 'max_grid_size_y',   MaxGridSizeX2 )
      CALL PP % query ( 'max_grid_size_z',   MaxGridSizeX3 )
      CALL PP % query ( 'blocking_factor_x', BlockingFactorX1 )
      CALL PP % query ( 'blocking_factor_y', BlockingFactorX2 )
      CALL PP % query ( 'blocking_factor_z', BlockingFactorX3 )
      CALL PP % get   ( 'max_level',         MaxLevel )
      CALL PP % query ( 'UseTiling'        , UseTiling        )
    CALL amrex_parmparse_destroy( PP )

    ! --- Equation of state parameters EoS.* ---
    OpacityTableName_AbEm = ''
    OpacityTableName_Iso  = ''
    OpacityTableName_NES  = ''
    OpacityTableName_Pair  = ''
    OpacityTableName_Brem  = ''
    CALL amrex_parmparse_build( PP, 'OP' )
      CALL PP % query( 'OpacityTableName_AbEm',OpacityTableName_AbEm )
      CALL PP % query( 'OpacityTableName_Iso', OpacityTableName_Iso )
      CALL PP % query( 'OpacityTableName_NES', OpacityTableName_NES )
      CALL PP % query( 'OpacityTableName_Pair', OpacityTableName_Pair )
      CALL PP % query( 'OpacityTableName_Brem', OpacityTableName_Brem )
    CALL amrex_parmparse_destroy( PP )

    Gamma_IDEAL     = 5.0_AR / 3.0_AR
    EquationOfState = 'IDEAL'
    EosTableName    = ''
    CALL amrex_parmparse_build( PP, 'EoS' )
      CALL PP % query( 'Gamma',           Gamma_IDEAL )
      CALL PP % query( 'EquationOfState', EquationOfState )
      CALL PP % query( 'EosTableName',    EosTableName    )
    CALL amrex_parmparse_destroy( PP )
    ! --- Positivitiy limiter parameters PL.* ---
    UsePositivityLimiter_TwoMoment = .FALSE.
    Min_1_TwoMoment                = 1.0e-12_AR
    Min_2_TwoMoment                = 1.0e-12_AR
    CALL amrex_parmparse_build( PP, 'PL' )
      CALL PP % query( 'UsePositivityLimiter_TwoMoment', UsePositivityLimiter_TwoMoment )
      CALL PP % query( 'Min_1_TwoMoment'               , Min_1_TwoMoment                )
      CALL PP % query( 'Min_2_TwoMoment'               , Min_2_TwoMoment                )
    CALL amrex_parmparse_destroy( PP )

    ! --- Positivitiy limiter parameters PL.* ---
    UseSlopeLimiter_TwoMoment = .FALSE.
    BetaTVD_TwoMoment = 1.0_AR
    CALL amrex_parmparse_build( PP, 'SL' )
      CALL PP % query( 'UseSlopeLimiter_TwoMoment', UseSlopeLimiter_TwoMoment )
      CALL PP % query( 'BetaTVD_TwoMoment'               , BetaTVD_TwoMoment                )
    CALL amrex_parmparse_destroy( PP )



    MaxGridSizeX = [ MaxGridSizeX1, MaxGridSizeX2, MaxGridSizeX3 ]

    nLevels = MaxLevel + 1

    CALL InitializeProgramHeader &
           ( ProgramName_Option = TRIM( ProgramName ), &
             nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX, bcX_Option = bcX, &
             xL_Option = xL, xR_Option = xR, &
             nE_Option = nE, swE_Option = swE, bcE_Option = bcE, eL_option = eL, eR_option = eR, &
             Verbose_Option = amrex_parallel_ioprocessor() )


    IF( amrex_parallel_ioprocessor() ) &
      CALL DescribeUnitsDisplay



    IF( nDimsX .NE. amrex_spacedim )THEN
      WRITE(*,'(A)') 'ERROR'
      WRITE(*,'(A)') '-----'
      WRITE(*,'(A)') 'thornado nDimsX different from AMReX amrex_spacedim.'
      WRITE(*,'(A)') 'Check DIM parameter in GNUmakefile. Stopping...'
      STOP
    END IF

    ALLOCATE( StepNo(0:nLevels-1) )
    StepNo = 0

    ALLOCATE( dt(0:nLevels-1) )
    dt = -100.0e0_AR

    ALLOCATE( t(0:nLevels-1) )
    t = 0.0e0_AR


    CALL InitializeDataAMReX( nLevels )

  END SUBROUTINE MyAmrInit


  SUBROUTINE MyAmrFinalize

    CALL FinalizeDataAMReX( nLevels )

    DEALLOCATE( GEOM )
    DEALLOCATE( DM   )
    DEALLOCATE( BA   )

    DEALLOCATE( t      )
    DEALLOCATE( dt     )
    DEALLOCATE( StepNo )

  END SUBROUTINE MyAmrFinalize



END MODULE InputParsingModule
