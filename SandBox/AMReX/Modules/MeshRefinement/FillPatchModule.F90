MODULE FillPatchModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_amr_module, ONLY: &
    amrex_geom, &
    amrex_ref_ratio
#if defined( THORNADO_USE_MESHREFINEMENT )
  USE amrex_interpolater_module, ONLY: &
    amrex_interp_dg
#endif
  USE amrex_fillpatch_module, ONLY: &
    amrex_fillpatch, &
    amrex_fillcoarsepatch
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
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm
  USE Euler_MeshRefinementModule, ONLY: &
    nFine, &
    vpCoarseToFineProjectionMatrix

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE InputParsingModule, ONLY: &
    UseTiling, &
    DEBUG
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_GeometryModule, ONLY: &
    ApplyBoundaryConditions_Geometry_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE MF_TimersModule, ONLY: &
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
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst
    LOGICAL             , INTENT(in), OPTIONAL :: &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option

    LOGICAL :: ApplyBoundaryConditions_Euler, &
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

    ApplyBoundaryConditions_Euler = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Euler_Option ) ) &
      ApplyBoundaryConditions_Euler = ApplyBoundaryConditions_Euler_Option

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

      CALL amrex_fillpatch &
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

      CALL amrex_fillpatch &
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

    IF( ApplyBoundaryConditions_Euler ) &
      CALL ApplyBoundaryConditions_Euler_MF( FineLevel, MF_dst )

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_PointWise_Scalar


  SUBROUTINE FillPatch_PointWise_Vector &
    ( FineLevel, MF, &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    LOGICAL             , INTENT(in), OPTIONAL :: &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option

    LOGICAL :: ApplyBoundaryConditions_Euler, &
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

    ApplyBoundaryConditions_Euler = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Euler_Option ) ) &
      ApplyBoundaryConditions_Euler = ApplyBoundaryConditions_Euler_Option

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

      CALL amrex_fillpatch &
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

      CALL amrex_fillpatch &
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

    IF( ApplyBoundaryConditions_Euler ) &
      CALL ApplyBoundaryConditions_Euler_MF( FineLevel, MF(FineLevel) )

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_PointWise_Vector


  SUBROUTINE FillCoarsePatch_PointWise &
    ( FineLevel, MF, &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)
    LOGICAL             , INTENT(in), OPTIONAL :: &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option

#if defined( THORNADO_USE_MESHREFINEMENT )

    LOGICAL :: ApplyBoundaryConditions_Euler, &
               ApplyBoundaryConditions_Geometry

    INTEGER :: iErr

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t          = Zero

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_Euler = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Euler_Option ) ) &
      ApplyBoundaryConditions_Euler = ApplyBoundaryConditions_Euler_Option

    ApplyBoundaryConditions_Geometry = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Geometry_Option ) ) &
      ApplyBoundaryConditions_Geometry = ApplyBoundaryConditions_Geometry_Option

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

    CALL amrex_fillcoarsepatch &
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

    IF( ApplyBoundaryConditions_Geometry ) &
      CALL ApplyBoundaryConditions_Geometry_MF( FineLevel, MF(FineLevel) )

    IF( ApplyBoundaryConditions_Euler ) &
      CALL ApplyBoundaryConditions_Euler_MF( FineLevel, MF(FineLevel) )

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

#endif

  END SUBROUTINE FillCoarsePatch_PointWise


  SUBROUTINE FillPatch_Conservative_Scalar &
    ( FineLevel, MF_uGF, MF_uGF_tmp, MF_src, MF_dst, &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:), MF_uGF_tmp
    TYPE(amrex_multifab), INTENT(inout) :: MF_src(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst
    LOGICAL             , INTENT(in)   , OPTIONAL :: &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel), SqrtGm_tmp

    INTEGER :: iErr

    LOGICAL :: ApplyBoundaryConditions_Euler, &
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

    ApplyBoundaryConditions_Euler = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Euler_Option ) ) &
      ApplyBoundaryConditions_Euler = ApplyBoundaryConditions_Euler_Option

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
             ( SqrtGm_tmp, MF_uGF_tmp % BA, &
                           MF_uGF_tmp % DM, nDOFX, swX )

      CALL SqrtGm(FineLevel-1) % COPY &
             ( MF_uGF(FineLevel-1), nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

      CALL SqrtGm(FineLevel  ) % COPY &
             ( MF_uGF(FineLevel  ), nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

      CALL SqrtGm_tmp % COPY &
             ( MF_uGF_tmp         , nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

    END IF

    IF( FineLevel .EQ. 0 )THEN

      CALL amrex_fillpatch &
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

      CALL amrex_fillpatch &
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

    IF( ApplyBoundaryConditions_Euler ) &
      CALL ApplyBoundaryConditions_Euler_MF( FineLevel, MF_dst )

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_destroy( SqrtGm_tmp )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_Conservative_Scalar


  SUBROUTINE FillPatch_Conservative_Vector &
    ( FineLevel, MF_uGF, MF, &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)
    LOGICAL             , INTENT(in)   , OPTIONAL :: &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: iErr

    LOGICAL :: ApplyBoundaryConditions_Euler, &
               ApplyBoundaryConditions_Geometry

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t_old_fine = Zero
    REAL(DP), PARAMETER :: t_new_fine = Zero
    REAL(DP), PARAMETER :: t          = Zero

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_Euler = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Euler_Option ) ) &
      ApplyBoundaryConditions_Euler = ApplyBoundaryConditions_Euler_Option

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

      CALL amrex_fillpatch &
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

      CALL amrex_fillpatch &
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

    IF( ApplyBoundaryConditions_Euler ) &
      CALL ApplyBoundaryConditions_Euler_MF( FineLevel, MF(FineLevel) )

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_Conservative_Vector


  SUBROUTINE FillCoarsePatch_Conservative &
    ( FineLevel, MF_uGF, MF, &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option )

    INTEGER             , INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)
    LOGICAL             , INTENT(in)   , OPTIONAL :: &
      ApplyBoundaryConditions_Euler_Option, &
      ApplyBoundaryConditions_Geometry_Option

#if defined( THORNADO_USE_MESHREFINEMENT )

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: iErr

    LOGICAL :: ApplyBoundaryConditions_Euler, &
               ApplyBoundaryConditions_Geometry

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = Zero
    REAL(DP), PARAMETER :: t_new_crse = Zero
    REAL(DP), PARAMETER :: t          = Zero

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    ApplyBoundaryConditions_Euler = .FALSE.
    IF( PRESENT( ApplyBoundaryConditions_Euler_Option ) ) &
      ApplyBoundaryConditions_Euler = ApplyBoundaryConditions_Euler_Option

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

      CALL SqrtGm(FineLevel-1) % COPY &
             ( MF_uGF(FineLevel-1), nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ) % BA, &
                                    MF_uGF(FineLevel  ) % DM, nDOFX, swX )

      CALL SqrtGm(FineLevel) % COPY &
             ( MF_uGF(FineLevel  ), nDOFX*(iGF_SqrtGm-1)+1, 1, nDOFX, swX )

    END IF

    ! Assume MF_old_crse = MF_new_crse = MF
    ! Assume t_old_crse  = t_new_crse  = t

    ALLOCATE( lo_bc(1:nDimsX,MF(FineLevel)%ncomp()) )
    ALLOCATE( hi_bc(1:nDimsX,MF(FineLevel)%ncomp()) )

    lo_bc = amrex_bc_bogus
    hi_bc = amrex_bc_bogus

    CALL amrex_fillcoarsepatch &
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

    IF( ApplyBoundaryConditions_Euler ) &
      CALL ApplyBoundaryConditions_Euler_MF( FineLevel, MF(FineLevel) )

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


END MODULE FillPatchModule
