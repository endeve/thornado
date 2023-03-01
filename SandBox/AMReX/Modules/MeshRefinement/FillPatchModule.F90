MODULE FillPatchModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
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
  USE amrex_amr_module, ONLY: &
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
    nNodes, &
    nDOFX, &
    nDimsX
  USE GeometryFieldsModule, ONLY: &
    iGF_SqrtGm

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_UtilitiesModule, ONLY: &
    MultiplyWithMetric
  USE InputParsingModule, ONLY: &
    UseTiling, &
    t_old, &
    t_new, &
    swX, &
    DEBUG
  USE MF_TimersModule, ONLY: &
    TimersStart_AMReX, &
    TimersStop_AMReX, &
    Timer_AMReX_FillPatch

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FillPatch, FillCoarsePatch

  INTERFACE FillPatch
    MODULE PROCEDURE FillPatch_Scalar_WithMetric
    MODULE PROCEDURE FillPatch_Vector_WithMetric
  END INTERFACE FillPatch

CONTAINS


  SUBROUTINE FillPatch_Scalar_WithMetric &
    ( FineLevel, MF_uGF, MF_src, MF_dst )

    INTEGER,              INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_src(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: nF, iErr

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    nF = MF_dst % nComp() / nDOFX

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL FillPatch_Scalar_WithMetric, FineLevel: ', FineLevel

      END IF

    END IF

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel-1), MF_uGF(FineLevel-1) % BA, &
                                    MF_uGF(FineLevel-1) % DM, nDOFX, swX )
      CALL SqrtGm(FineLevel-1) % COPY &
             ( MF_uGF(FineLevel-1), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ) % BA, &
                                    MF_uGF(FineLevel  ) % DM, nDOFX, swX )
      CALL SqrtGm(FineLevel  ) % COPY &
             ( MF_uGF(FineLevel  ), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

      CALL MultiplyWithMetric &
             ( FineLevel-1, SqrtGm(FineLevel-1), MF_src, nF, +1 )
      CALL MultiplyWithMetric &
             ( FineLevel  , SqrtGm(FineLevel  ), MF_src, nF, +1 )

    END IF

    CALL FillPatch_Scalar( FineLevel, MF_src, MF_dst )

    IF( FineLevel .GT. 0 )THEN

      CALL MultiplyWithMetric &
             ( FineLevel-1, SqrtGm(FineLevel-1), MF_src, nF, -1 )
      CALL MultiplyWithMetric &
             ( FineLevel  , SqrtGm(FineLevel  ), MF_src, nF, -1 )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_Scalar_WithMetric


  SUBROUTINE FillPatch_Vector_WithMetric( FineLevel, MF_uGF, MF )

    INTEGER,              INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: nF, iErr

    CALL TimersStart_AMReX( Timer_AMReX_FillPatch )

    nF = MF(0) % nComp() / nDOFX

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      WRITE(*,'(4x,A,I3.3)') &
        'CALL FillPatch_Vector_WithMetric, FineLevel: ', FineLevel

    END IF

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel-1), MF_uGF(FineLevel-1) % BA, &
                                    MF_uGF(FineLevel-1) % DM, nDOFX, swX )
      CALL SqrtGm(FineLevel-1) % COPY &
             ( MF_uGF(FineLevel-1), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ) % BA, &
                                    MF_uGF(FineLevel  ) % DM, nDOFX, swX )
      CALL SqrtGm(FineLevel  ) % COPY &
             ( MF_uGF(FineLevel  ), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

      CALL MultiplyWithMetric &
             ( FineLevel-1, SqrtGm(FineLevel-1), MF, nF, +1 )

      CALL MultiplyWithMetric &
             ( FineLevel  , SqrtGm(FineLevel  ), MF, nF, +1 )

    END IF

    CALL FillPatch_Vector( FineLevel, MF )

    IF( FineLevel .GT. 0 )THEN

      CALL MultiplyWithMetric &
             ( FineLevel-1, SqrtGm(FineLevel-1), MF, nF, -1 )
      CALL MultiplyWithMetric &
             ( FineLevel  , SqrtGm(FineLevel  ), MF, nF, -1 )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

    CALL TimersStop_AMReX( Timer_AMReX_FillPatch )

  END SUBROUTINE FillPatch_Vector_WithMetric


  SUBROUTINE FillCoarsePatch( FineLevel, MF_uGF, MF )

    INTEGER,              INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: nF

    nF = MF(0) % nComp() / nDOFX

    IF( FineLevel .GT. 0 )THEN

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel-1), MF_uGF(FineLevel-1) % BA, &
                                 MF_uGF(FineLevel-1) % DM, nDOFX, swX )
      CALL SqrtGm(FineLevel-1) % COPY &
             ( MF_uGF(FineLevel-1), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

      CALL amrex_multifab_build &
             ( SqrtGm(FineLevel  ), MF_uGF(FineLevel  ) % BA, &
                                 MF_uGF(FineLevel  ) % DM, nDOFX, swX )
      CALL SqrtGm(FineLevel  ) % COPY &
             ( MF_uGF(FineLevel  ), 1+nDOFX*(iGF_SqrtGm-1), 1, nDOFX, swX )

      CALL MultiplyWithMetric &
             ( FineLevel-1, SqrtGm(FineLevel-1), MF, nF, +1 )
      CALL MultiplyWithMetric &
             ( FineLevel  , SqrtGm(FineLevel  ), MF, nF, +1 )

    END IF

    CALL FillCoarsePatch_Vector( FineLevel, MF )

    IF( FineLevel .GT. 0 )THEN

      CALL MultiplyWithMetric &
             ( FineLevel-1, SqrtGm(FineLevel-1), MF, nF, -1 )
      CALL MultiplyWithMetric &
             ( FineLevel  , SqrtGm(FineLevel  ), MF, nF, -1 )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

  END SUBROUTINE FillCoarsePatch


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE FillPatch_Scalar( FineLevel, MF_src, MF_dst )

    INTEGER,              INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(in)    :: MF_src(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = 0.0_DP
    REAL(DP), PARAMETER :: t_new_crse = 0.0_DP
    REAL(DP), PARAMETER :: t_old_fine = 0.0_DP
    REAL(DP), PARAMETER :: t_new_fine = 0.0_DP
    REAL(DP), PARAMETER :: t          = 0.0_DP

    ! Assume MF_old_crse = MF_new_crse = MF_old_fine = MF_new_fine = MF
    ! Assume t_old_crse  = t_new_crse  = t_old_fine  = t_new_fine  = t

    IF( FineLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF_dst, &
                            t_old_fine, MF_src(FineLevel), &
                            t_new_fine, MF_src(FineLevel), &
                            amrex_geom(FineLevel), FillPhysicalBC_Dummy, &
                            t, sComp, dComp, MF_dst % nComp() )

    ELSE

#if defined( THORNADO_USE_MESHREFINEMENT )

      ALLOCATE( lo_bc(1:nDimsX,MF_src(FineLevel)%ncomp()) )
      ALLOCATE( hi_bc(1:nDimsX,MF_src(FineLevel)%ncomp()) )

      lo_bc = amrex_bc_bogus
      hi_bc = amrex_bc_bogus

      CALL amrex_fillpatch( MF_dst, &
                            t_old_crse, MF_src(FineLevel-1), &
                            t_new_crse, MF_src(FineLevel-1), &
                            amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
                            t_old_fine, MF_src(FineLevel  ), &
                            t_new_fine, MF_src(FineLevel  ), &
                            amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
                            t, sComp, dComp, MF_dst % nComp(), &
                            amrex_ref_ratio(FineLevel-1), &
                            amrex_interp_dg, &
                            lo_bc, hi_bc )

      DEALLOCATE( hi_bc )
      DEALLOCATE( lo_bc )

#endif

    END IF

  END SUBROUTINE FillPatch_Scalar


  SUBROUTINE FillPatch_Vector( FineLevel, MF )

    INTEGER,              INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = 0.0_DP
    REAL(DP), PARAMETER :: t_new_crse = 0.0_DP
    REAL(DP), PARAMETER :: t_old_fine = 0.0_DP
    REAL(DP), PARAMETER :: t_new_fine = 0.0_DP
    REAL(DP), PARAMETER :: t          = 0.0_DP

    ! Assume MF_old_crse = MF_new_crse = MF_old_fine = MF_new_fine = MF
    ! Assume t_old_crse  = t_new_crse  = t_old_fine  = t_new_fine  = t

    IF( FineLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF(FineLevel), &
                            t_old_fine, MF(FineLevel), &
                            t_new_fine, MF(FineLevel), &
                            amrex_geom(FineLevel), FillPhysicalBC_Dummy, &
                            t, sComp, dComp, MF(FineLevel) % nComp() )

    ELSE

#if defined( THORNADO_USE_MESHREFINEMENT )

      ALLOCATE( lo_bc(1:nDimsX,MF(FineLevel)%ncomp()) )
      ALLOCATE( hi_bc(1:nDimsX,MF(FineLevel)%ncomp()) )

      lo_bc = amrex_bc_bogus
      hi_bc = amrex_bc_bogus

      CALL amrex_fillpatch( MF(FineLevel), &
                            t_old_crse, MF(FineLevel-1), &
                            t_new_crse, MF(FineLevel-1), &
                            amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
                            t_old_fine, MF(FineLevel  ), &
                            t_new_fine, MF(FineLevel  ), &
                            amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
                            t, sComp, dComp, MF(FineLevel) % nComp(), &
                            amrex_ref_ratio(FineLevel-1), &
                            amrex_interp_dg, &
                            lo_bc, hi_bc )

      DEALLOCATE( hi_bc )
      DEALLOCATE( lo_bc )

#endif

    END IF

  END SUBROUTINE FillPatch_Vector


  SUBROUTINE FillCoarsePatch_Vector( FineLevel, MF )

    INTEGER,              INTENT(in)    :: FineLevel
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    INTEGER, ALLOCATABLE :: lo_bc(:,:), hi_bc(:,:)

    ! Dummy variables. Only matter when interpolating in time
    REAL(DP), PARAMETER :: t_old_crse = 0.0_DP
    REAL(DP), PARAMETER :: t_new_crse = 0.0_DP
    REAL(DP), PARAMETER :: t          = 0.0_DP

    ! Assume MF_old_crse = MF_new_crse = MF
    ! Assume t_old_crse  = t_new_crse  = t

#if defined( THORNADO_USE_MESHREFINEMENT )

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
             t, sComp, dComp, MF(FineLevel) % nComp(), &
             amrex_ref_ratio(FineLevel-1), &
             amrex_interp_dg, lo_bc, hi_bc )

    DEALLOCATE( hi_bc )
    DEALLOCATE( lo_bc )

#endif

  END SUBROUTINE FillCoarsePatch_Vector


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

    END IF

  END SUBROUTINE FillPhysicalBC_Dummy


END MODULE FillPatchModule
