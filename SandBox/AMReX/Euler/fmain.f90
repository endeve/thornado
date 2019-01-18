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
    WriteFieldsAMReX_PlotFile, &
    ReadCheckpointFile

  ! --- Local Modules ---

  USE MF_GeometryModule,        ONLY: &
    MF_ComputeGeometryX
  USE MF_InitializationModule,  ONLY: &
    MF_InitializeFields
  USE MF_Euler_UtilitiesModule, ONLY: &
    MF_ComputeFromConserved

  ! --- Checkpoint ---
  USE InputOutputModuleAMReX ! Needed to import WriteFieldsAMReX_Checkpoint
  USE amrex_amr_module ! To call amrex_amrcore_init
  USE MyAmrDataModule
  USE MyAmrModule

  IMPLICIT NONE

  INTEGER :: iCycle
  INTEGER :: iDim
  TYPE(amrex_box)       :: BX
  TYPE(amrex_boxarray)  :: BA
  TYPE(amrex_distromap) :: DM
  TYPE(amrex_geometry), ALLOCATABLE :: GEOM(:)

  ! --- Checkpoint ---
  INTEGER :: iLevel

  ! --- Initialize AMReX ---

  CALL amrex_init( )

  CALL amrex_amrcore_init() ! Gets refinement ratio

  ! --- Parse Parameter File ---
  CALL MyAmrInit

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

  SELECT CASE ( CoordSys )
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
    WRITE(*,'(A5,A24,A)') &
      '', 'CoordinateSystem = ', CoordinateSystem

  END IF

  IF( amrex_parallel_ioprocessor() )THEN

    WRITE(*,*)
    WRITE(*,'(A4,A24,3I7.6)') '', 'nX          =', nX
    WRITE(*,'(A4,A24,3I7.6)') '', 'swX         =', swX
    WRITE(*,'(A4,A24,3I7.6)') '', 'MaxGridSize =', MaxGridSize

  END IF

  CALL InitializeProgramHeader &
         ( ProgramName_Option = TRIM( ProgramName ), &
           nNodes_Option = nNodes, nX_Option = nX, swX_Option = swX, &
           xL_Option = xL, xR_Option = xR )

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

  BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )

  CALL amrex_boxarray_build( BA, BX )

  CALL BA % maxSize( MaxGridSize )

  ALLOCATE( GEOM(0:nLevels) )

  DO iLevel = 0, nLevels
    CALL amrex_geometry_build( GEOM(iLevel), BX )
  END DO

  CALL amrex_distromap_build( DM, BA )

  DO iLevel = 0, nLevels
    CALL amrex_multifab_build( MF_uGF(iLevel), BA, DM, nDOFX * nGF, swX(1) )
    CALL amrex_multifab_build( MF_uCF(iLevel), BA, DM, nDOFX * nCF, swX(1) )
    CALL amrex_multifab_build( MF_uPF(iLevel), BA, DM, nDOFX * nPF, swX(1) )
    CALL amrex_multifab_build( MF_uAF(iLevel), BA, DM, nDOFX * nAF, swX(1) )
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

  CALL WriteFieldsAMReX_Checkpoint &
         ( StepNo, nLevels, dt, t_new, &
           MF_uGF % BA % P, &
           MF_uGF % P, &
           MF_uCF % P, &
           MF_uPF % P, &
           MF_uAF % P )

  iCycle = 0

  ! --- Evolution goes here

  WRITE(*,*) 'Reading from checkpoint file...'
  CALL MyAmrFinalize
  CALL ReadCheckpointFile()

  ! --- END of evolution

  DO iDim = 1, 3

    CALL DestroyMesh( MeshX(iDim) )

  END DO

  ! --- Finalize AMReX ---

  CALL MyAmrFinalize

  DO iLevel = 0, nLevels
    CALL amrex_geometry_destroy( GEOM(iLevel) )
  END DO
  DEALLOCATE( GEOM )

  CALL amrex_amrcore_finalize()
  CALL amrex_finalize( )

END PROGRAM main

