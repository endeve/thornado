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
    amrex_ref_ratio, &
    amrex_interp_dg
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

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodes, &
    nDOFX
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
    lo_bc, &
    hi_bc, &
    swX, &
    DEBUG

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FillPatch, FillCoarsePatch

  INTERFACE FillPatch
    MODULE PROCEDURE FillPatch_WithMetric_Scalar
    MODULE PROCEDURE FillPatch_WithMetric_Vector
  END INTERFACE FillPatch

CONTAINS


  SUBROUTINE FillPatch_WithMetric_Scalar &
    ( FineLevel, Time, MF_uGF, MF_src, MF_dst )

    INTEGER,              INTENT(in)    :: FineLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_src(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: nF, iErr

    nF = MF_dst % nComp() / nDOFX

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      IF( amrex_parallel_ioprocessor() )THEN

        WRITE(*,'(4x,A,I3.3)') &
          'CALL FillPatch_WithMetric_Scalar, FineLevel: ', FineLevel

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
             ( SqrtGm(FineLevel-1), MF_src(FineLevel-1), nF, +1 )
      CALL MultiplyWithMetric &
             ( SqrtGm(FineLevel  ), MF_src(FineLevel  ), nF, +1 )

    END IF

    CALL FillPatch_Scalar( FineLevel, Time, MF_src, MF_dst )

    IF( FineLevel .GT. 0 )THEN

      CALL MultiplyWithMetric &
             ( SqrtGm(FineLevel-1), MF_src(FineLevel-1), nF, -1 )
      CALL MultiplyWithMetric &
             ( SqrtGm(FineLevel  ), MF_src(FineLevel  ), nF, -1 )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

  END SUBROUTINE FillPatch_WithMetric_Scalar


  SUBROUTINE FillPatch_WithMetric_Vector( FineLevel, Time, MF_uGF, MF )

    INTEGER,              INTENT(in)    :: FineLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF    (0:)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    TYPE(amrex_multifab) :: SqrtGm(FineLevel-1:FineLevel)

    INTEGER :: nF, iErr

    nF = MF(0) % nComp() / nDOFX

    IF( DEBUG )THEN

      CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

      WRITE(*,'(4x,A,I3.3)') &
        'CALL FillPatch_WithMetric_Vector, FineLevel: ', FineLevel

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

      CALL MultiplyWithMetric( SqrtGm(FineLevel-1), MF(FineLevel-1), nF, +1 )
      CALL MultiplyWithMetric( SqrtGm(FineLevel  ), MF(FineLevel  ), nF, +1 )

    END IF

    CALL FillPatch_Vector( FineLevel, Time, MF )

    IF( FineLevel .GT. 0 )THEN

      CALL MultiplyWithMetric( SqrtGm(FineLevel-1), MF(FineLevel-1), nF, -1 )
      CALL MultiplyWithMetric( SqrtGm(FineLevel  ), MF(FineLevel  ), nF, -1 )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

  END SUBROUTINE FillPatch_WithMetric_Vector


  SUBROUTINE FillCoarsePatch( FineLevel, Time, MF_uGF, MF )

    INTEGER,              INTENT(in)    :: FineLevel
    REAL(DP),             INTENT(in)    :: Time
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

      CALL MultiplyWithMetric( SqrtGm(FineLevel-1), MF(FineLevel-1), nF, +1 )
      CALL MultiplyWithMetric( SqrtGm(FineLevel  ), MF(FineLevel  ), nF, +1 )

    END IF

    CALL FillCoarsePatch_Vector( FineLevel, Time, MF )

    IF( FineLevel .GT. 0 )THEN

      CALL MultiplyWithMetric( SqrtGm(FineLevel-1), MF(FineLevel-1), nF, -1 )
      CALL MultiplyWithMetric( SqrtGm(FineLevel  ), MF(FineLevel  ), nF, -1 )

      CALL amrex_multifab_destroy( SqrtGm(FineLevel-1) )
      CALL amrex_multifab_destroy( SqrtGm(FineLevel  ) )

    END IF

  END SUBROUTINE FillCoarsePatch


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE FillPatch_Scalar( FineLevel, Time, MF_src, MF_dst )

    INTEGER,              INTENT(in)    :: FineLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(in)    :: MF_src(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    IF( FineLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF_dst, &
                            Time, MF_src(FineLevel), &
                            Time, MF_src(FineLevel), &
                            amrex_geom(FineLevel), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, MF_dst % nComp() )

    ELSE

      CALL amrex_fillpatch( MF_dst, &
                            Time, MF_src(FineLevel-1), &
                            Time, MF_src(FineLevel-1), &
                            amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
                            Time, MF_src(FineLevel  ), &
                            Time, MF_src(FineLevel  ), &
                            amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, MF_dst % nComp(), &
                            amrex_ref_ratio(FineLevel-1), &
                            amrex_interp_dg, &
                            lo_bc, hi_bc )

    END IF

  END SUBROUTINE FillPatch_Scalar


  SUBROUTINE FillPatch_Vector( FineLevel, Time, MF )

    INTEGER,              INTENT(in)    :: FineLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    ! Assume MF_old = MF_new = MF and t_old = t_new = t

    IF( FineLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF(FineLevel), &
                            Time, MF(FineLevel), &
                            Time, MF(FineLevel), &
                            amrex_geom(FineLevel), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, MF(FineLevel) % nComp() )

    ELSE

      CALL amrex_fillpatch( MF(FineLevel), &
                            Time, MF(FineLevel-1), &
                            Time, MF(FineLevel-1), &
                            amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
                            Time, MF(FineLevel  ), &
                            Time, MF(FineLevel  ), &
                            amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, MF(FineLevel) % nComp(), &
                            amrex_ref_ratio(FineLevel-1), &
                            amrex_interp_dg, &
                            lo_bc, hi_bc )

    END IF

  END SUBROUTINE FillPatch_Vector


  SUBROUTINE FillCoarsePatch_Vector( FineLevel, Time, MF )

    INTEGER,              INTENT(in)    :: FineLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    ! Assume t_old = t_new = t and MF_old = MF_new = MF

    CALL amrex_fillcoarsepatch &
           ( MF(FineLevel), &
             Time, MF(FineLevel-1), &
             Time, MF(FineLevel-1), &
             amrex_geom(FineLevel-1), FillPhysicalBC_Dummy, &
             amrex_geom(FineLevel  ), FillPhysicalBC_Dummy, &
             Time, sComp, dComp, MF(FineLevel) % nComp(), &
             amrex_ref_ratio(FineLevel-1), &
             amrex_interp_dg, lo_bc, hi_bc )

  END SUBROUTINE FillCoarsePatch_Vector


  SUBROUTINE FillPhysicalBC_Dummy( pMF, sComp, nComp, Time, pGEOM ) BIND(c)

    ! --- No INTENT here because amrex source code doesn't have it ---

    TYPE(c_ptr),    VALUE :: pMF, pGEOM
    INTEGER(c_int), VALUE :: sComp, nComp
    REAL(DP),       VALUE :: Time

    RETURN
!!$
!!$    TYPE(amrex_geometry) :: GEOM
!!$    TYPE(amrex_multifab) :: MF
!!$    TYPE(amrex_mfiter)   :: MFI
!!$    INTEGER              :: pLo(4), pHi(4)
!!$    REAL(DP), CONTIGUOUS, POINTER :: p(:,:,:,:)
!!$
!!$    IF( .NOT. amrex_is_all_periodic() )THEN
!!$
!!$      GEOM = pGEOM
!!$      MF   = pMF
!!$
!!$      !$OMP PARALLEL PRIVATE(MFI,p,pLo,pHi)
!!$      CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )
!!$
!!$      DO WHILE( MFI % next() )
!!$
!!$        p => MF % DataPtr( MFI )
!!$
!!$        ! Check if part of this box is outside the domain
!!$        IF( .NOT. GEOM % DOMAIN % CONTAINS( p ) )THEN
!!$
!!$          pLo = LBOUND( p )
!!$          pHi = UBOUND( p )
!!$
!!$          CALL amrex_filcc &
!!$                 ( p, pLo, pHi, &
!!$                   GEOM % DOMAIN % lo, GEOM % DOMAIN % hi, &
!!$                   GEOM % dX, &
!!$                   GEOM % get_physical_location( pLo ), &
!!$                   lo_bc, hi_bc )
!!$
!!$        END IF
!!$
!!$      END DO
!!$      !$OMP END PARALLEL
!!$
!!$      CALL amrex_mfiter_destroy( MFI )
!!$
!!$    END IF

  END SUBROUTINE FillPhysicalBC_Dummy


END MODULE FillPatchModule
