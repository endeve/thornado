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
  USE GeometryFieldsModule,             ONLY: &
    nGF
  USE FluidFieldsModule,                ONLY: &
    nCF, nPF, nAF

  ! --- Local Modules ---

  USE MF_GeometryModule,    ONLY: &
    MF_ComputeGeometryX
  USE InitializationModule, ONLY: &
    InitializeFields

  IMPLICIT NONE

  CHARACTER(LEN=:), ALLOCATABLE :: ProgramName
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

  CALL PP % get( "t_end",   t_end )
  CALL PP % get( "dt_wrt",  dt_wrt )
  CALL PP % get( "nNodes",  nNodes )
  CALL PP % get( "nStages", nStages )
  CALL PP % get( "ProgramName", ProgramName )

  CALL amrex_parmparse_destroy( PP )

  PRINT*, "t_end   = ", t_end
  PRINT*, "dt_wrt  = ", dt_wrt
  PRINT*, "nNodes  = ", nNodes
  PRINT*, "nStages = ", nStages
  PRINT*, "nDimsX  = ", amrex_spacedim
  PRINT*, "Name    = ", ProgramName, LEN( ProgramName )

  CALL amrex_parmparse_build( PP, "geometry" )

  CALL PP % get( "coord_sys", CoordinateSystem )

  CALL amrex_parmparse_destroy( PP )

  PRINT*, "CoordinateSystem = ", CoordinateSystem

  CALL amrex_parmparse_build( PP, "amr" )

  CALL PP % getarr( "n_cell", n_cell )
  CALL PP % getarr( "max_grid_size", max_grid_size )

  CALL amrex_parmparse_destroy( PP )

  PRINT*, "n_cell = ", n_cell
  PRINT*, "max_grid_size = ", max_grid_size

  DO iDim = 1, amrex_spacedim
    nX (iDim) = n_cell(iDim)
    swX(iDim) = nGhost
  END DO

  CALL InitializeProgramHeader &
         ( ProgramName_Option = TRIM( ProgramName ), &
           nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX )

  CALL DescribeProgramHeaderX

  CALL InitializePolynomialBasisX_Lagrange

  CALL InitializePolynomialBasisX_Legendre

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  BX = amrex_box( [ 1, 1, 1 ], [ n_cell(1), n_cell(2), n_cell(3) ] )

  PRINT*, "BX % Lo = ", BX % Lo
  PRINT*, "BX % Hi = ", BX % Hi
  PRINT*, "BX % nodal = ", BX % nodal

  CALL amrex_print( BX )

  CALL amrex_boxarray_build( BA, BX )

  CALL BA % maxSize( max_grid_size )

  CALL amrex_geometry_build( GEOM, BX )

  PRINT*, "GEOM % dx = ", GEOM % dx

  CALL amrex_distromap_build( DM, BA )

  CALL amrex_multifab_build( MF_uGF, BA, DM, nDOFX * nGF, nGhost )

  CALL amrex_multifab_build( MF_uCF, BA, DM, nDOFX * nCF, nGhost )

  CALL amrex_multifab_build( MF_uPF, BA, DM, nDOFX * nPF, nGhost )

  CALL amrex_multifab_build( MF_uAF, BA, DM, nDOFX * nAF, nGhost )

  CALL amrex_distromap_destroy( DM )

  CALL amrex_boxarray_destroy( BA )

  CALL MF_ComputeGeometryX( MF_uGF )

  CALL InitializeFields( TRIM( ProgramName ), GEOM, MF_uCF )

  CALL amrex_multifab_destroy( MF_uGF )

  CALL amrex_multifab_destroy( MF_uCF )

  CALL amrex_multifab_destroy( MF_uPF )

  CALL amrex_multifab_destroy( MF_uAF )

  CALL amrex_geometry_destroy( GEOM )

  ! --- Finalize AMReX ---

  CALL amrex_finalize( )

END PROGRAM main

