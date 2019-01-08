PROGRAM main

  ! --- AMReX Modules ---

  USE amrex_base_module
  USE amrex_fort_module

  ! --- thornado Modules ---

  USE KindModule,                       ONLY: &
    DP
  USE ProgramHeaderModule,              ONLY: &
    InitializeProgramHeader,                  &
    DescribeProgramHeaderX,                   &
    nDOFX, nNodesX
  USE PolynomialBasisModuleX_Lagrange,  ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE PolynomialBasisModuleX_Legendre,  ONLY: &
    InitializePolynomialBasisX_Legendre
  USE ReferenceElementModuleX,          ONLY: &
    InitializeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE MeshModule,                       ONLY: &
    MeshX, CreateMesh, DestroyMesh
  USE EquationOfStateModule,            ONLY: &
    InitializeEquationOfState
  USE GeometryFieldsModule,             ONLY: &
    nGF, CoordinateSystem
  USE FluidFieldsModule,                ONLY: &
    nCF, nPF, nAF
  USE InputOutputModuleAMReX,           ONLY: &
    WriteFieldsAMReX_PlotFile

  ! --- Local Modules ---

  USE MF_GeometryModule,        ONLY: &
    MF_ComputeGeometryX
  USE MF_InitializationModule,  ONLY: &
    MF_InitializeFields
  USE MF_Euler_UtilitiesModule, ONLY: &
    MF_ComputeFromConserved

  ! --- Checkpoint ---
  USE MyRestartModule

  IMPLICIT NONE

  CHARACTER(LEN=:), ALLOCATABLE :: ProgramName
  INTEGER :: iCycle
  INTEGER :: iDim
  INTEGER :: nNodes
  INTEGER :: nStages
  INTEGER :: nX(3) = 1, swX(3) = 0
  INTEGER :: coord_sys
  INTEGER, PARAMETER :: nGhost = 2
  INTEGER, ALLOCATABLE  :: n_cell(:)
  INTEGER, ALLOCATABLE  :: max_grid_size(:)
  REAL(amrex_real)      :: t_end
  REAL(amrex_real)      :: dt_wrt
  REAL(amrex_real)      :: Gamma_IDEAL
  TYPE(amrex_parmparse) :: PP
  TYPE(amrex_box)       :: BX
  TYPE(amrex_boxarray)  :: BA
  TYPE(amrex_distromap) :: DM
  TYPE(amrex_geometry), ALLOCATABLE :: GEOM(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uGF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uCF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uPF(:)
  TYPE(amrex_multifab), ALLOCATABLE :: MF_uAF(:)

  ! --- Checkpoint ---
  INTEGER,              ALLOCATABLE :: StepNo(:)
  REAL(amrex_real),     ALLOCATABLE :: dt(:), t_new(:)
  INTEGER                           :: FinestLevel, iLevel, nLevels
  TYPE(amrex_multifab), ALLOCATABLE :: MF_Chk(:)

  ! --- Initialize AMReX ---

  CALL amrex_init( )

  ! --- Parse Parameter File ---

  CALL amrex_parmparse_build( PP )

  CALL PP % get( "t_end",       t_end )
  CALL PP % get( "dt_wrt",      dt_wrt )
  CALL PP % get( "nNodes",      nNodes )
  CALL PP % get( "nStages",     nStages )
  CALL PP % get( "ProgramName", ProgramName )
  CALL PP % get( "Gamma",       Gamma_IDEAL )

  CALL amrex_parmparse_destroy( PP )

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(A4,A6,A)') '', 'Name: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A4,A24,ES10.3E2)') '',   't_end =', t_end
    WRITE(*,'(A4,A24,ES10.3E2)') '',  'dt_wrt =', dt_wrt
    WRITE(*,'(A4,A24,I7.6)')     '',  'nNodes =', nNodes
    WRITE(*,'(A4,A24,I7.6)')     '', 'nStages =', nStages
    WRITE(*,'(A4,A24,I3.2)')     '',  'nDimsX =', amrex_spacedim
    WRITE(*,'(A4,A24,ES10.3E2)') '',   'Gamma =', Gamma_IDEAL

  END IF

  CALL amrex_parmparse_build( PP, "geometry" )

  CALL PP % get( "coord_sys", coord_sys )

  CALL amrex_parmparse_destroy( PP )

  SELECT CASE ( coord_sys )
    CASE ( 0 )
      CoordinateSystem = 'CARTESIAN'
    CASE ( 1 )
      CoordinateSystem = 'CYLINDRICAL'
    CASE ( 2 )
      CoordinateSystem = 'SPHERICAL'
    CASE DEFAULT
      CoordinateSystem = 'CARTESIAN'
  END SELECT

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(A4,A24,A)') &
      '', 'CoordinateSystem =', CoordinateSystem

  END IF

  CALL amrex_parmparse_build( PP, "amr" )

  CALL PP % getarr( "n_cell", n_cell )
  CALL PP % getarr( "max_grid_size", max_grid_size )
  CALL PP % get   ( "nLevels", nLevels )

  CALL amrex_parmparse_destroy( PP )

  FinestLevel = nLevels
  ALLOCATE( StepNo(0:nLevels) )
  StepNo = 0
  ALLOCATE( dt    (0:nLevels) )
  dt = 1.0e-4_amrex_real
  ALLOCATE( t_new (0:nLevels) )
  t_new = 0.0_amrex_real

  ALLOCATE( GEOM  (0:nLevels) )
  ALLOCATE( MF_uGF(0:nLevels) )
  ALLOCATE( MF_uCF(0:nLevels) )
  ALLOCATE( MF_uPF(0:nLevels) )
  ALLOCATE( MF_uAF(0:nLevels) )
  ALLOCATE( MF_Chk(0:nLevels) )

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(A4,A24,3I7.6)') '',        'n_cell =', n_cell
    WRITE(*,'(A4,A24,3I7.6)') '', 'max_grid_size =', max_grid_size

  END IF

  DO iDim = 1, amrex_spacedim
    nX (iDim) = n_cell(iDim)
    swX(iDim) = nGhost
  END DO

  CALL InitializeProgramHeader &
         ( ProgramName_Option = TRIM( ProgramName ), &
           nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX )

  IF( amrex_parallel_ioprocessor() )THEN

    CALL DescribeProgramHeaderX

  END IF

  CALL InitializePolynomialBasisX_Lagrange

  CALL InitializePolynomialBasisX_Legendre

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma_IDEAL )

  BX = amrex_box( [ 1, 1, 1 ], [ n_cell(1), n_cell(2), n_cell(3) ] )

  CALL amrex_boxarray_build( BA, BX )

  CALL BA % maxSize( max_grid_size )

  DO iLevel = 0, nLevels
    CALL amrex_geometry_build( GEOM(iLevel), BX )
  END DO

  CALL amrex_distromap_build( DM, BA )

  DO iLevel = 0, nLevels
    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, nGhost )
    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, nGhost )
    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, nGhost )
    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, nGhost )
    CALL amrex_multifab_build( MF_Chk(iLevel), BA, DM, &
                                 nDOFX * ( nGF + nCF + nPF + nAF ), nGhost )
  END DO

  CALL amrex_distromap_destroy( DM )

  CALL amrex_boxarray_destroy( BA )

  DO iDim = 1, 3

    CALL CreateMesh &
           ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
             amrex_problo(iDim), amrex_probhi(iDim) )

  END DO

  DO iLevel = 0, nLevels
    CALL MF_ComputeGeometryX( MF_uGF(iLevel) )
    CALL MF_InitializeFields &
           ( TRIM( ProgramName ), MF_uGF(iLevel), MF_uCF(iLevel) )
    CALL MF_ComputeFromConserved &
           ( MF_uGF(iLevel), MF_uCF(iLevel), MF_uPF(iLevel), MF_uAF(iLevel) )
  END DO

  CALL WriteFieldsAMReX_PlotFile &
         ( 0.0_DP, nLevels, GEOM, StepNo, &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF )

  CALL WriteCheckpointFile &
         ( StepNo, FinestLevel, dt, t_new, &
           MF_uGF(0:FinestLevel) % BA % P, &
           MF_uGF(0:FinestLevel) % P, &
           MF_uCF(0:FinestLevel) % P, &
           MF_uPF(0:FinestLevel) % P, &
           MF_uAF(0:FinestLevel) % P )

  iCycle = 0

  ! --- Evolution goes here
!!$  CALL ReadCheckpointFile()
  ! --- END of evolution

  DO iDim = 1, 3

    CALL DestroyMesh( MeshX(iDim) )

  END DO

  DO iLevel = 0, nLevels
    CALL amrex_geometry_destroy( GEOM  (iLevel) )
    CALL amrex_multifab_destroy( MF_uGF(iLevel) )
    CALL amrex_multifab_destroy( MF_uCF(iLevel) )
    CALL amrex_multifab_destroy( MF_uPF(iLevel) )
    CALL amrex_multifab_destroy( MF_uAF(iLevel) )
  END DO

  DEALLOCATE( StepNo )
  DEALLOCATE( dt )
  DEALLOCATE( t_new )

  ! --- Finalize AMReX ---

  CALL amrex_finalize( )

END PROGRAM main

