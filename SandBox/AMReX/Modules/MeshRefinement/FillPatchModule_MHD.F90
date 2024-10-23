MODULE FillPatchModule_MHD

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_amr_module, ONLY: &
    amrex_geom, &
    amrex_ref_ratio
#if defined( THORNADO_USE_MESHREFINEMENT )
  USE thornado_amrex_interpolater_module, ONLY: &
    amrex_interp_dg
#endif
  USE thornado_amrex_fillpatch_module, ONLY: &
    thornado_amrex_fillpatch, &
    thornado_amrex_fillcoarsepatch
  USE amrex_geometry_module, ONLY: &
    amrex_geometry, &
    amrex_is_all_periodic
  USE amrex_filcc_module, ONLY: &
    amrex_filcc
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor
  USE amrex_bc_types_module, ONLY: &
    amrex_bc_bogus

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    swX, &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshType, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm, &
    CoordinateSystem
  USE MHD_MeshRefinementModule, ONLY: &
    nFine, &
    vpCoarseToFineProjectionMatrix

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE InputParsingModule, ONLY: &
    UseTiling, &
    DEBUG
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF, &
    UpdateSpatialMetric_MF
  USE MF_MHD_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_MHD_MF
  USE MF_TimersModule_MHD, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_FillPatch

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FillPatch, FillCoarsePatch

  INTERFACE FillPatch
    MODULE PROCEDURE FillPatch_PointWise_Scalar ! Only called by RemakeLevel
    MODULE PROCEDURE FillPatch_PointWise_Vector
    MODULE PROCEDURE FillPatch_Conservative_Scalar ! Only called by RemakeLevel
    MODULE PROCEDURE FillPatch_Conservative_Vector
  END INTERFACE FillPatch

  ! Only called by FillCoarsePatch
  INTERFACE FillCoarsePatch
    MODULE PROCEDURE FillCoarsePatch_PointWise
    MODULE PROCEDURE FillCoarsePatch_Conservative
  END INTERFACE FillCoarsePatch

CONTAINS


  SUBROUTINE FillPatch_PointWise_Scalar &
    ( FineLevel, MF, MF_dst, &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst
    LOGICAL             , INTENT(in), OPTIONAL :: &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option

    LOGICAL :: ApplyBoundaryConditions_MHD, &
               ApplyBoundaryConditions_Geometry

    INTEGER :: iErr

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. They only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t_old_fine = Zero
    REAL(DP), PARAMETER :: t_new_fine = Zero
    REAL(DP), PARAMETER :: t          = Zero

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_MHD = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_MHD_Option ) ) &
      ApplyBoundaryConditions_MHD = ApplyBoundaryConditions_MHD_Option

    ApplyBoundaryConditions_Geometry = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Geometry_Option ) ) &
      ApplyBoundaryConditions_Geometry = ApplyBoundaryConditions_Geometry_Option

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL FillPatch_PointWise_Scalar, FineLevel: ', FineLevel

      END IF

    END IF

    ! Assume MF_old_crse = MF_new_crse = MF_old_fine = MF_new_fine = MF

    IF( FineLevel .EQ. 0 )THEN

      CALL thornado_amrex_fillpatch &
             ( MF_dst, &
               t_old_fine, MF(FineLevel), &
               t_new_fine, MF(FineLevel), &
               amrex_geom(FineLevel), FillPhysicalBC_Dummy, &
               t, 1, 1, MF_dst % nComp() )

    ELSE

#if defined( THORNADO_USE_MESHREFINEMENT )

      ALLOCATE( lo_bc(1:nDimsX,MF(FineLevel)%ncomp()) )
      ALLOCATE( hi_bc(1:nDimsX,MF(FineLevel)%ncomp()) )

      lo_bc = amrex_bc_bogus
      hi_bc = amrex_bc_bogus

      CALL thornado_amrex_fillpatch &
             ( MF_dst, &
               t_old_crse, MF(FineLevel-1), &
               t_new_crse, MF(FineLevel-1), &
               amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
               t_old_fine, MF(FineLevel  ), &
               t_new_fine, MF(FineLevel  ), &
               amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
               t, 1, 1, MF_dst % nComp(), &
               amrex_ref_ratio(FineLevel-1), &
               amrex_interp_dg, &
               lo_bc, hi_bc, &
               nFine, nDOFX, vpCoarseToFineProjectionMatrix )

      DEALLOCATE( hi_bc )
      DEALLOCATE( lo_bc )

#endif

    END IF

    IF( ApplyBoundaryConditions_Geometry ) &
      CALL ApplyBoundaryConditions_Geometry_MF( FineLevel, MF_dst )

    IF( ApplyBoundaryConditions_MHD ) &
      CALL ApplyBoundaryConditions_MHD_MF( t, FineLevel, MF_dst )

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_PointWise_Scalar


  SUBROUTINE FillPatch_PointWise_Vector &
    ( FineLevel, MF, &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    LOGICAL             , INTENT(in), OPTIONAL :: &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option

    LOGICAL :: ApplyBoundaryConditions_MHD, &
               ApplyBoundaryConditions_Geometry

    INTEGER :: iErr

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t_old_fine = Zero
    REAL(DP), PARAMETER :: t_new_fine = Zero
    REAL(DP), PARAMETER :: t          = Zero

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_MHD = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_MHD_Option ) ) &
      ApplyBoundaryConditions_MHD = ApplyBoundaryConditions_MHD_Option

    ApplyBoundaryConditions_Geometry = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Geometry_Option ) ) &
      ApplyBoundaryConditions_Geometry = ApplyBoundaryConditions_Geometry_Option

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      WRITE(*,'(4x,A,I3.3)') &
        'CALL FillPatch_PointWise_Vector, FineLevel: ', FineLevel

    END IF

    ! Assume MF_old_crse = MF_new_crse = MF_old_fine = MF_new_fine = MF
    ! Assume t_old_crse  = t_new_crse  = t_old_fine  = t_new_fine  = t

    IF( FineLevel .EQ. 0 )THEN

      CALL thornado_amrex_fillpatch &
             ( MF(FineLevel), &
               t_old_fine, MF(FineLevel), &
               t_new_fine, MF(FineLevel), &
               amrex_geom(FineLevel), FillPhysicalBC_Dummy, &
               t, 1, 1, MF(FineLevel) % nComp() )

    ELSE

#if defined( THORNADO_USE_MESHREFINEMENT )

      ALLOCATE( lo_bc(1:nDimsX,MF(FineLevel)%ncomp()) )
      ALLOCATE( hi_bc(1:nDimsX,MF(FineLevel)%ncomp()) )

      lo_bc = amrex_bc_bogus
      hi_bc = amrex_bc_bogus

      CALL thornado_amrex_fillpatch &
             ( MF(FineLevel), &
               t_old_crse, MF(FineLevel-1), &
               t_new_crse, MF(FineLevel-1), &
               amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
               t_old_fine, MF(FineLevel  ), &
               t_new_fine, MF(FineLevel  ), &
               amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
               t, 1, 1, MF(FineLevel) % nComp(), &
               amrex_ref_ratio(FineLevel-1), &
               amrex_interp_dg, &
               lo_bc, hi_bc, &
               nFine, nDOFX, vpCoarseToFineProjectionMatrix )

      DEALLOCATE( hi_bc )
      DEALLOCATE( lo_bc )

#endif

    END IF

    IF( ApplyBoundaryConditions_Geometry ) &
      CALL ApplyBoundaryConditions_Geometry_MF( FineLevel, MF(FineLevel) )

    IF( ApplyBoundaryConditions_MHD ) &
      CALL ApplyBoundaryConditions_MHD_MF( t, FineLevel, MF(FineLevel) )

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_PointWise_Vector


  SUBROUTINE FillCoarsePatch_PointWise &
    ( FineLevel, MF, &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option, &
      UpdateSpatialMetric_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    LOGICAL             , INTENT(in), OPTIONAL :: &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option, &
      UpdateSpatialMetric_Option

#if defined( THORNADO_USE_MESHREFINEMENT )

    LOGICAL :: ApplyBoundaryConditions_MHD, &
               ApplyBoundaryConditions_Geometry, &
               UpdateSpatialMetric

    INTEGER :: iErr

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t          = Zero

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_MHD = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_MHD_Option ) ) &
      ApplyBoundaryConditions_MHD = ApplyBoundaryConditions_MHD_Option

    ApplyBoundaryConditions_Geometry = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Geometry_Option ) ) &
      ApplyBoundaryConditions_Geometry = ApplyBoundaryConditions_Geometry_Option

    UpdateSpatialMetric = .FALSE.
    IF( PRESENT( UpdateSpatialMetric_Option ) ) &
      UpdateSpatialMetric = UpdateSpatialMetric_Option

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      WRITE(*,'(4x,A,I3.3)') &
        'CALL FillCoarsePatch_PointWise, FineLevel: ', FineLevel

    END IF

    ! Assume MF_old_crse = MF_new_crse = MF

    ALLOCATE( lo_bc(1:nDimsX,MF(FineLevel)%ncomp()) )
    ALLOCATE( hi_bc(1:nDimsX,MF(FineLevel)%ncomp()) )

    lo_bc = amrex_bc_bogus
    hi_bc = amrex_bc_bogus

    CALL thornado_amrex_fillcoarsepatch &
           ( MF(FineLevel), &
             t_old_crse, MF(FineLevel-1), &
             t_new_crse, MF(FineLevel-1), &
             amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
             amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
             t, MF(FineLevel) % nComp(), &
             amrex_ref_ratio(FineLevel-1), &
             amrex_interp_dg, lo_bc, hi_bc, &
             nFine, nDOFX, vpCoarseToFineProjectionMatrix )

    DEALLOCATE( hi_bc )
    DEALLOCATE( lo_bc )

    IF( UpdateSpatialMetric )THEN

      CALL UpdateSpatialMetric_MF( FineLevel, MF(FineLevel) )

    END IF

    IF( ApplyBoundaryConditions_Geometry ) &
      CALL ApplyBoundaryConditions_Geometry_MF( FineLevel, MF(FineLevel) )

    IF( ApplyBoundaryConditions_MHD ) &
      CALL ApplyBoundaryConditions_MHD_MF( t, FineLevel, MF(FineLevel) )

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

#endif

  END SUBROUTINE FillCoarsePatch_PointWise


  SUBROUTINE FillPatch_Conservative_Scalar &
    ( FineLevel, MF_uGF, MF_uGF_tmp, MF_src, MF_dst, &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:), MF_uGF_tmp
    TYPE(amrex_multifab), INTENT(inout) :: MF_src(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst
    LOGICAL             , INTENT(in), OPTIONAL :: &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel), SqrtGm_tmp

    INTEGER :: iErr

    LOGICAL :: ApplyBoundaryConditions_MHD, &
               ApplyBoundaryConditions_Geometry

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t_old_fine = Zero
    REAL(DP), PARAMETER :: t_new_fine = Zero
    REAL(DP), PARAMETER :: t          = Zero

    ! Assume MF_old_crse = MF_new_crse = MF_old_fine = MF_new_fine = MF

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_MHD = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_MHD_Option ) ) &
      ApplyBoundaryConditions_MHD = ApplyBoundaryConditions_MHD_Option

    ApplyBoundaryConditions_Geometry = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Geometry_Option ) ) &
      ApplyBoundaryConditions_Geometry = ApplyBoundaryConditions_Geometry_Option

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL FillPatch_Conservative_Scalar, FineLevel: ', FineLevel

      END IF

    END IF

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel-1), MF_uGF(FineLevel-1) % BA, &
                                    MF_uGF(FineLevel-1) % DM, nDOFX, swX )
      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ) % BA, &
                                    MF_uGF(FineLevel  ) % DM, nDOFX, swX )
      CALL amrex_multifab_build &
             ( SqrtGm_tmp         , MF_uGF_tmp          % BA, &
                                    MF_uGF_tmp          % DM, nDOFX, swX )

      CALL PopulateWithFlatSpaceMetric_MF &
             ( FineLevel-1, SqrtGm(FineLevel-1), swX )
      CALL PopulateWithFlatSpaceMetric_MF &
             ( FineLevel  , SqrtGm(FineLevel  ), swX )
      CALL PopulateWithFlatSpaceMetric_MF &
             ( FineLevel  , SqrtGm_tmp         , swX )

    END IF

    IF( FineLevel .EQ. 0 )THEN

      CALL thornado_amrex_fillpatch &
             ( MF_dst, &
               t_old_fine, MF_src(FineLevel), &
               t_new_fine, MF_src(FineLevel), &
               amrex_geom(FineLevel), FillPhysicalBC_Dummy, &
               t, 1, 1, MF_dst % nComp() )

    ELSE

#if defined( THORNADO_USE_MESHREFINEMENT )

      ALLOCATE( lo_bc(1:nDimsX,MF_src(FineLevel)%ncomp()) )
      ALLOCATE( hi_bc(1:nDimsX,MF_src(FineLevel)%ncomp()) )

      lo_bc = amrex_bc_bogus
      hi_bc = amrex_bc_bogus

      CALL thornado_amrex_fillpatch &
             ( MF_dst, SqrtGm_tmp, &
               t_old_crse, MF_src(FineLevel-1), SqrtGm(FineLevel-1), &
               t_new_crse, MF_src(FineLevel-1), SqrtGm(FineLevel-1), &
               amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
               t_old_fine, MF_src(FineLevel  ), SqrtGm(FineLevel  ), &
               t_new_fine, MF_src(FineLevel  ), SqrtGm(FineLevel  ), &
               amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
               t, 1, 1, MF_dst % nComp(), &
               amrex_ref_ratio(FineLevel-1), &
               amrex_interp_dg, &
               lo_bc, hi_bc, &
               nFine, nDOFX, vpCoarseToFineProjectionMatrix )

      DEALLOCATE( hi_bc )
      DEALLOCATE( lo_bc )

#endif

    END IF

    IF( ApplyBoundaryConditions_Geometry ) &
      CALL ApplyBoundaryConditions_Geometry_MF( FineLevel, MF_dst )

    IF( ApplyBoundaryConditions_MHD ) &
      CALL ApplyBoundaryConditions_MHD_MF( t, FineLevel, MF_dst )

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_destroy( SqrtGm_tmp )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_Conservative_Scalar


  SUBROUTINE FillPatch_Conservative_Vector &
    ( FineLevel, MF_uGF, MF, &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)
    LOGICAL             , INTENT(in)   , OPTIONAL :: &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: iErr

    LOGICAL :: ApplyBoundaryConditions_MHD, &
               ApplyBoundaryConditions_Geometry

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t_old_fine = Zero
    REAL(DP), PARAMETER :: t_new_fine = Zero
    REAL(DP), PARAMETER :: t          = Zero

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_MHD = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_MHD_Option ) ) &
      ApplyBoundaryConditions_MHD = ApplyBoundaryConditions_MHD_Option

    ApplyBoundaryConditions_Geometry = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Geometry_Option ) ) &
      ApplyBoundaryConditions_Geometry = ApplyBoundaryConditions_Geometry_Option

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      WRITE(*,'(4x,A,I3.3)') &
        'CALL FillPatch_Conservative_Vector, FineLevel: ', FineLevel

    END IF

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel-1), MF_uGF(FineLevel-1) % BA, &
                                    MF_uGF(FineLevel-1) % DM, nDOFX, swX )

      CALL SqrtGm(FineLevel-1) % COPY &
             ( MF_uGF(FineLevel-1), nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ) % BA, &
                                    MF_uGF(FineLevel  ) % DM, nDOFX, swX )

      CALL SqrtGm(FineLevel  ) % COPY &
             ( MF_uGF(FineLevel  ), nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

    END IF

    ! Assume MF_old_crse = MF_new_crse = MF_old_fine = MF_new_fine = MF
    ! Assume t_old_crse  = t_new_crse  = t_old_fine  = t_new_fine  = t

    IF( FineLevel .EQ. 0 )THEN

      CALL thornado_amrex_fillpatch &
             ( MF(FineLevel), &
               t_old_fine, MF(FineLevel), &
               t_new_fine, MF(FineLevel), &
               amrex_geom(FineLevel), FillPhysicalBC_Dummy, &
               t, 1, 1, MF(FineLevel) % nComp() )

    ELSE

#if defined( THORNADO_USE_MESHREFINEMENT )

      ALLOCATE( lo_bc(1:nDimsX,MF(FineLevel)%ncomp()) )
      ALLOCATE( hi_bc(1:nDimsX,MF(FineLevel)%ncomp()) )

      lo_bc = amrex_bc_bogus
      hi_bc = amrex_bc_bogus

      CALL thornado_amrex_fillpatch &
             ( MF(FineLevel), SqrtGm(FineLevel), &
               t_old_crse, MF(FineLevel-1), SqrtGm(FineLevel-1), &
               t_new_crse, MF(FineLevel-1), SqrtGm(FineLevel-1), &
               amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
               t_old_fine, MF(FineLevel  ), SqrtGm(FineLevel  ), &
               t_new_fine, MF(FineLevel  ), SqrtGm(FineLevel  ), &
               amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
               t, 1, 1, MF(FineLevel) % nComp(), &
               amrex_ref_ratio(FineLevel-1), &
               amrex_interp_dg, &
               lo_bc, hi_bc, &
               nFine, nDOFX, vpCoarseToFineProjectionMatrix )

      DEALLOCATE( hi_bc )
      DEALLOCATE( lo_bc )

#endif

    END IF

    IF( ApplyBoundaryConditions_Geometry ) &
      CALL ApplyBoundaryConditions_Geometry_MF( FineLevel, MF(FineLevel) )

    IF( ApplyBoundaryConditions_MHD ) &
      CALL ApplyBoundaryConditions_MHD_MF( t, FineLevel, MF(FineLevel) )

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_Conservative_Vector


  SUBROUTINE FillCoarsePatch_Conservative &
    ( FineLevel, MF_uGF, MF, &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)
    LOGICAL             , INTENT(in)   , OPTIONAL :: &
      ApplyBoundaryConditions_MHD_Option, &
      ApplyBoundaryConditions_Geometry_Option

#if defined( THORNADO_USE_MESHREFINEMENT )

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: iErr

    LOGICAL :: ApplyBoundaryConditions_MHD, &
               ApplyBoundaryConditions_Geometry

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t          = Zero

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_MHD = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_MHD_Option ) ) &
      ApplyBoundaryConditions_MHD = ApplyBoundaryConditions_MHD_Option

    ApplyBoundaryConditions_Geometry = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Geometry_Option ) ) &
      ApplyBoundaryConditions_Geometry = ApplyBoundaryConditions_Geometry_Option

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      WRITE(*,'(4x,A,I3.3)') &
        'CALL FillCoarsePatch_Conservative, FineLevel: ', FineLevel

    END IF

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel-1), MF_uGF(FineLevel-1) % BA, &
                                    MF_uGF(FineLevel-1) % DM, nDOFX, swX )

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ) % BA, &
                                    MF_uGF(FineLevel  ) % DM, nDOFX, swX )

      CALL PopulateWithFlatSpaceMetric_MF &
             ( FineLevel-1, SqrtGm(FineLevel-1), swX )
      CALL PopulateWithFlatSpaceMetric_MF &
             ( FineLevel  , SqrtGm(FineLevel  ), swX )

    END IF

    ! Assume MF_old_crse = MF_new_crse = MF
    ! Assume t_old_crse  = t_new_crse  = t

    ALLOCATE( lo_bc(1:nDimsX,MF(FineLevel)%ncomp()) )
    ALLOCATE( hi_bc(1:nDimsX,MF(FineLevel)%ncomp()) )

    lo_bc = amrex_bc_bogus
    hi_bc = amrex_bc_bogus

    CALL thornado_amrex_fillcoarsepatch &
           ( MF(FineLevel), SqrtGm(FineLevel), &
             t_old_crse, MF(FineLevel-1), SqrtGm(FineLevel-1), &
             t_new_crse, MF(FineLevel-1), SqrtGm(FineLevel-1), &
             amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
             amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
             t, MF(FineLevel) % nComp(), &
             amrex_ref_ratio(FineLevel-1), &
             amrex_interp_dg, &
             lo_bc, hi_bc, &
             nFine, nDOFX, vpCoarseToFineProjectionMatrix )

    DEALLOCATE( hi_bc )
    DEALLOCATE( lo_bc )

    IF( ApplyBoundaryConditions_Geometry ) &
      CALL ApplyBoundaryConditions_Geometry_MF( FineLevel, MF(FineLevel) )

    IF( ApplyBoundaryConditions_MHD ) &
      CALL ApplyBoundaryConditions_MHD_MF( t, FineLevel, MF(FineLevel) )

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

#endif

  END SUBROUTINE FillCoarsePatch_Conservative


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE FillPhysicalBC_Dummy( pMF, sComp, nComp, Time, pGEOM ) BIND(c)

    ! --- No INTENT here because amrex source code doesn't have it ---

    TYPE(c_ptr)   , VALUE :: pMF, pGEOM
    INTEGER(c_int), VALUE :: sComp, nComp
    REAL(DP)      , VALUE :: Time

    TYPE(amrex_geometry) :: GEOM
    TYPE(amrex_multifab) :: MF
    TYPE(amrex_mfiter)   :: MFI

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    REAL(DP), CONTIGUOUS, POINTER, DIMENSION(:,:,:,:) :: p

    INTEGER :: plo(4), phi(4)

    IF( .NOT. amrex_is_all_periodic() )THEN

      GEOM = pGEOM
      MF   = pMF

      CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )

      DO WHILE( MFI % next() )

        p => mf%dataptr(mfi)

        ! Part of this box is outside the domain
        IF( .NOT. GEOM % domain % CONTAINS(p) )THEN

          plo = LBOUND(p)
          phi = UBOUND(p)

          ALLOCATE( lo_bc(1:nDimsX,plo(4):phi(4)) )
          ALLOCATE( hi_bc(1:nDimsX,plo(4):phi(4)) )

          lo_bc = amrex_bc_bogus
          hi_bc = amrex_bc_bogus

          CALL amrex_filcc &
                 ( p, plo, phi, &
                   GEOM % domain % lo, GEOM % domain % hi, &
                   GEOM % dx, &
                   GEOM % get_physical_location( plo ), &
                   lo_bc, hi_bc )

          ! amrex_filcc doesn't fill EXT_DIR (see amrex_bc_types_module for a list of bc types
          ! In that case, the user needs to fill it.

          DEALLOCATE( hi_bc )
          DEALLOCATE( lo_bc )

        END IF

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END IF

  END SUBROUTINE FillPhysicalBC_Dummy


  SUBROUTINE PopulateWithFlatSpaceMetric_MF( iLevel, MF_SqrtGm, swXX )

    INTEGER             , INTENT(in)  :: iLevel, swXX(3)
    TYPE(amrex_multifab), INTENT(inout) :: MF_SqrtGm

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    TYPE(MeshType) :: MeshXX(3)
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B(3), iX_E(3)
    INTEGER  :: iX1, iX2, iX3, iNX, iNX1, iNX2
    REAL(DP) :: X1, X2, h1, h2, h3

    CALL CreateMesh_MF( iLevel, MeshXX )

#if defined( THORNADO_OMP )
    !$OMP PARALLEL &
    !$OMP PRIVATE( BX, MFI, uGF, &
    !$OMP          iX_B0, iX_E0, iX_B, iX_E, iNX1, iNX2, X1, X2, h1, h2, h3 )
#endif

    CALL amrex_mfiter_build( MFI, MF_SqrtGm, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_SqrtGm % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B  = iX_B0 - swXX
      iX_E  = iX_E0 + swXX

      DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
      DO iX1 = iX_B(1), iX_E(1)
      DO iNX = 1      , nDOFX

        iNX1 = NodeNumberTableX(1,iNX)
        iNX2 = NodeNumberTableX(2,iNX)

        X1 = NodeCoordinate( MeshXX(1), iX1, iNX1 )
        X2 = NodeCoordinate( MeshXX(2), iX2, iNX2 )

        IF( TRIM( CoordinateSystem ) .EQ. 'CYLINDRICAL' )THEN

          h1  = One
          h2  = One
          h3  = ABS( X1 )

        ELSE IF( TRIM( CoordinateSystem ) .EQ. 'SPHERICAL' )THEN

          h1  = One
          h2  = ABS( X1 )
          h3  = ABS( X1 * SIN( X2 ) )

        ELSE

          h1  = One
          h2  = One
          h3  = One

        END IF

        uGF(iX1,iX2,iX3,iNX) = h1 * h2 * h3

      END DO
      END DO
      END DO
      END DO

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
    !$OMP END PARALLEL
#endif

    CALL DestroyMesh_MF( MeshXX )

  END SUBROUTINE PopulateWithFlatSpaceMetric_MF


END MODULE FillPatchModule_MHD
