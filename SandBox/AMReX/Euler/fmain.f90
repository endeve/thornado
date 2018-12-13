PROGRAM main

  ! --- AMReX Modules ---

  USE amrex_base_module
  USE amrex_fort_module, only: amrex_spacedim

  ! --- thornado Modules ---

  USE KindModule,                       ONLY: &
    DP
  USE ProgramHeaderModule,              ONLY: &
    InitializeProgramHeader,                  &
    DescribeProgramHeaderX,                   &
    nDOFX
  USE PolynomialBasisModuleX_Lagrange,  ONLY: &
    InitializePolynomialBasisX_Lagrange
  USE PolynomialBasisModuleX_Legendre,  ONLY: &
    InitializePolynomialBasisX_Legendre
  USE ReferenceElementModuleX,          ONLY: &
    InitializeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange
  USE EquationOfStateModule,            ONLY: &
    InitializeEquationOfState
  USE GeometryFieldsModule,             ONLY: &
    nGF
  USE FluidFieldsModule,                ONLY: &
    nCF, nPF, nAF
  USE InputOutputModuleAMReX,           ONLY: &
    WriteFieldsAMReX_PlotFile, &
    WriteFieldsAMReX_Checkpoint

  ! --- Local Modules ---

  USE MF_GeometryModule,        ONLY: &
    MF_ComputeGeometryX
  USE MF_InitializationModule,  ONLY: &
    MF_InitializeFields
  USE MF_Euler_UtilitiesModule, ONLY: &
    MF_ComputeFromConserved

  IMPLICIT NONE

  CHARACTER(LEN=:), ALLOCATABLE :: ProgramName
  INTEGER :: iCycle
  INTEGER :: iDim
  INTEGER :: nNodes
  INTEGER :: nStages
  INTEGER :: nX(3) = 1, swX(3) = 0
  INTEGER :: CoordinateSystem
  INTEGER, PARAMETER :: nGhost = 2
  INTEGER, ALLOCATABLE  :: n_cell(:)
  INTEGER, ALLOCATABLE  :: max_grid_size(:)
  REAL(amrex_real)      :: t_end
  REAL(amrex_real)      :: dt_wrt
  REAL(amrex_real)      :: Gamma_IDEAL
  TYPE(amrex_parmparse) :: PP
  TYPE(amrex_box)       :: BX
  TYPE(amrex_boxarray)  :: BA
  TYPE(amrex_geometry)  :: GEOM
  TYPE(amrex_distromap) :: DM
  TYPE(amrex_multifab)  :: MF_uGF
  TYPE(amrex_multifab)  :: MF_uCF
  TYPE(amrex_multifab)  :: MF_uPF
  TYPE(amrex_multifab)  :: MF_uAF

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

  CALL PP % get( "coord_sys", CoordinateSystem )

  CALL amrex_parmparse_destroy( PP )

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(A4,A24,I3.2)') &
      '', 'CoordinateSystem =', CoordinateSystem

  END IF

  CALL amrex_parmparse_build( PP, "amr" )

  CALL PP % getarr( "n_cell", n_cell )
  CALL PP % getarr( "max_grid_size", max_grid_size )

  CALL amrex_parmparse_destroy( PP )

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

  CALL amrex_geometry_build( GEOM, BX )

  CALL amrex_distromap_build( DM, BA )

  CALL amrex_multifab_build( MF_uGF, BA, DM, nDOFX * nGF, nGhost )

  CALL amrex_multifab_build( MF_uCF, BA, DM, nDOFX * nCF, nGhost )

  CALL amrex_multifab_build( MF_uPF, BA, DM, nDOFX * nPF, nGhost )

  CALL amrex_multifab_build( MF_uAF, BA, DM, nDOFX * nAF, nGhost )

  CALL amrex_distromap_destroy( DM )

  CALL amrex_boxarray_destroy( BA )

  CALL MF_ComputeGeometryX( MF_uGF )

  CALL MF_InitializeFields( TRIM( ProgramName ), MF_uGF, MF_uCF )

  CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

  CALL WriteFieldsAMReX_PlotFile &
         ( 0.0_DP, GEOM, &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF )

  iCycle = 0
  CALL WriteFieldsAMReX_Checkpoint( iCycle, 0.0_DP, GEOM, &
                                    MF_uGF, MF_uCF, MF_uPF, MF_uAF )

  CALL amrex_multifab_destroy( MF_uGF )

  CALL amrex_multifab_destroy( MF_uCF )

  CALL amrex_multifab_destroy( MF_uPF )

  CALL amrex_multifab_destroy( MF_uAF )

  CALL amrex_geometry_destroy( GEOM )

  ! --- Finalize AMReX ---

  CALL amrex_finalize( )

END PROGRAM main

