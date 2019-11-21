MODULE InitializationModule
  
  ! --- AMReX Modules ---
  USE amrex_fort_module, ONLY: &
    AR => amrex_real, &
    amrex_spacedim
  USE amrex_amr_module, ONLY: &
    amrex_init, &
    amrex_finalize
  USE amrex_amrcore_module, ONLY: &
    amrex_amrcore_init
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray,         &
    amrex_boxarray_build,   &
    amrex_boxarray_destroy, &
    amrex_print
  USE amrex_distromap_module, ONLY: &
    amrex_distromap,       &
    amrex_distromap_build, &
    amrex_distromap_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry, &
    amrex_geometry_build
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse,       &
    amrex_parmparse_build, &
    amrex_parmparse_destroy

  ! --- thornado Modules ---
  USE ProgramHeaderModule,              ONLY: &
    DescribeProgramHeaderX, &
    nDOFX,                  &
    nNodesX,                &
    nDOFE,                  &
    nNodesE,                &
    nDOF
  USE RadiationFieldsModule,            ONLY: &
    nCR, &
    nPR, &
    CreateRadiationFields
  ! --- Local modules ---
  USE MyAmrDataModule,                  ONLY: &
    MF_uPF, &
    MF_uCF
  USE MyAmrModule,                      ONLY: &
    t_end,                     &
    t,                         &
    dt,                        &
    t_wrt,                     &
    dt_wrt,                    &
    t_chk,                     &
    dt_chk,                    &
    CFL,                       &
    nNodes,                    &
    nStages,                   &
    nX,                        &
    nE,                        &
    swX,                       &
    swE,                       &
    bcX,                       &
    xL,                        &
    xR,                        &
    ProgramName,               &
    CoordSys,                  &
    StepNo,                    &
    nLevels,                   &
    iRestart,                  &
    MaxGridSizeX,              &
    BA,                        &
    DM,                        &
    GEOM,                      &
    MyAmrInit


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram

  REAL(AR), PARAMETER :: Zero = 0.0_AR
  REAL(AR), PARAMETER :: One  = 1.0_AR
  REAL(AR), PARAMETER :: Two  = 2.0_AR

CONTAINS


  SUBROUTINE InitializeProgram

    INTEGER               :: iLevel, iDim
    TYPE(amrex_parmparse) :: PP
    TYPE(amrex_box)       :: BX
    REAL(AR)              :: Mass

    ! --- Initialize AMReX ---
    CALL amrex_init()

    CALL amrex_amrcore_init()

    ! --- Parse parameter file ---
    CALL MyAmrInit

    IF( iRestart .LT. 0 )THEN
      BX = amrex_box( [ 1, 1, 1 ], [ nX(1), nX(2), nX(3) ] )
      ALLOCATE( BA(0:nLevels-1) )
      DO iLevel = 0, nLevels-1
        CALL amrex_boxarray_build( BA(iLevel), BX )
      END DO
      DO iLevel = 0, nLevels-1
        CALL BA(iLevel) % maxSize( MaxGridSizeX )
      END DO
      ALLOCATE( GEOM(0:nLevels-1) )
      ALLOCATE( DM  (0:nLevels-1) )
      DO iLevel = 0, nLevels-1
        CALL amrex_geometry_build ( GEOM(iLevel), BX )
        CALL amrex_distromap_build( DM  (iLevel), BA(iLevel) )
      END DO
    
      DO iLevel = 0, nLevels-1
        CALL amrex_multifab_build &
               ( MF_uPF(iLevel), BA(iLevel), DM(iLevel), nDOF * nPR, swX(1) )
        CALL MF_uPF(iLevel) % SetVal( Zero )
        CALL amrex_multifab_build &
               ( MF_uCF(iLevel), BA(iLevel), DM(iLevel), nDOF * nCR, swX(1) )
        CALL MF_uCF(iLevel) % SetVal( Zero )
      END DO 
  
      t     = Zero
      dt    = Zero
      t_wrt = dt_wrt
      t_chk = dt_chk

    ELSE
      print*, 'else'
    END IF

    CALL CreateRadiationFields( nX, swX, nE, swE )
    print*, nDOFX,nDOFE,nDOF
  END SUBROUTINE InitializeProgram  


END MODULE InitializationModule
