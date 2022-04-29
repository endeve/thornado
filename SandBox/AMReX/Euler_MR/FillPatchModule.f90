MODULE FillPatchModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
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


  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodes

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    t_old, &
    t_new, &
    lo_bc, &
    hi_bc, &
    lo_bc_uCF, &
    hi_bc_uCF

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FillPatch, FillCoarsePatch

  INTERFACE FillPatch
    MODULE PROCEDURE FillPatch_MultiTime
    MODULE PROCEDURE FillPatch_SingleTime
    MODULE PROCEDURE FillPatch_SingleTime_Copy
  END INTERFACE FillPatch

  INTERFACE FillCoarsePatch
    MODULE PROCEDURE FillCoarsePatch_MultiTime
    MODULE PROCEDURE FillCoarsePatch_SingleTime
  END INTERFACE FillCoarsePatch

CONTAINS


  SUBROUTINE FillPatch_MultiTime( iLevel, t_old, t_new, t, MF_old, MF_new, MF )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: t_old(0:nLevels-1), &
                                           t_new(0:nLevels-1), &
                                           t
    TYPE(amrex_multifab), INTENT(in)    :: MF_old(0:nLevels-1), &
                                           MF_new(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    IF( iLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF, &
                            t_old(iLevel), MF_old(iLevel), &
                            t_new(iLevel), MF_new(iLevel), &
                            amrex_geom(iLevel), FillPhysicalBC_Dummy, &
                            t, sComp, dComp, MF % nComp() )

    ELSE

      CALL amrex_fillpatch( MF, &
                            t_old(iLevel-1), MF_old(iLevel-1), &
                            t_new(iLevel-1), MF_new(iLevel-1), &
                            amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
                            t_old(iLevel  ), MF_old(iLevel  ), &
                            t_new(iLevel  ), MF_new(iLevel  ), &
                            amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
                            t, sComp, dComp, MF % nComp(), &
                            amrex_ref_ratio(iLevel-1), &
                            amrex_interp_dg, &
                            lo_bc, hi_bc )

    END IF

  END SUBROUTINE FillPatch_MultiTime


  SUBROUTINE FillPatch_SingleTime( iLevel, Time, MF )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:nLevels-1)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    ! Assume MF_old = MF_new = MF and t_old = t_new = t

    IF( iLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF(iLevel), &
                            Time, MF(iLevel), &
                            Time, MF(iLevel), &
                            amrex_geom(iLevel), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, MF(iLevel) % nComp() )

    ELSE

      CALL amrex_fillpatch( MF(iLevel), &
                            Time, MF(iLevel-1), &
                            Time, MF(iLevel-1), &
                            amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
                            Time, MF(iLevel  ), &
                            Time, MF(iLevel  ), &
                            amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, MF(iLevel) % nComp(), &
                            amrex_ref_ratio(iLevel-1), &
                            amrex_interp_dg, &
                            lo_bc, hi_bc )

    END IF

  END SUBROUTINE FillPatch_SingleTime


  SUBROUTINE FillPatch_SingleTime_Copy( iLevel, Time, MF_src, MF_dst )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(in)    :: MF_src(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_dst

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    IF( iLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF_dst, &
                            Time, MF_src(iLevel), &
                            Time, MF_src(iLevel), &
                            amrex_geom(iLevel), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, MF_dst % nComp() )

    ELSE

      CALL amrex_fillpatch( MF_dst, &
                            Time, MF_src(iLevel-1), &
                            Time, MF_src(iLevel-1), &
                            amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
                            Time, MF_src(iLevel  ), &
                            Time, MF_src(iLevel  ), &
                            amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, MF_dst % nComp(), &
                            amrex_ref_ratio(iLevel-1), &
                            amrex_interp_dg, &
                            lo_bc, hi_bc )

    END IF

  END SUBROUTINE FillPatch_SingleTime_Copy


  SUBROUTINE FillCoarsePatch_MultiTime &
    ( iLevel, t_old, t_new, t, MF_old, MF_new, MF )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: t_old(0:nLevels-1), &
                                           t_new(0:nLevels-1), &
                                           t
    TYPE(amrex_multifab), INTENT(in)    :: MF_old(0:nLevels-1), &
                                           MF_new(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF

    CALL amrex_fillcoarsepatch &
           ( MF, t_old(iLevel-1), MF_old(iLevel-1), &
                 t_new(iLevel-1), MF_new(iLevel-1), &
             amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
             amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
             t, 1, 1, MF % nComp(), amrex_ref_ratio(iLevel-1), &
             amrex_interp_dg, lo_bc, hi_bc )

  END SUBROUTINE FillCoarsePatch_MultiTime


  SUBROUTINE FillCoarsePatch_SingleTime( iLevel, Time, MF )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF(0:nLevels-1)

    INTEGER, PARAMETER :: sComp = 1, dComp = 1

    ! Assume t_old = t_new = t and MF_old = MF_new = MF

    CALL amrex_fillcoarsepatch &
           ( MF(iLevel), &
             Time, MF(iLevel-1), &
             Time, MF(iLevel-1), &
             amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
             amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
             Time, sComp, dComp, MF(iLevel) % nComp(), &
             amrex_ref_ratio(iLevel-1), &
             amrex_interp_dg, lo_bc, hi_bc )

  END SUBROUTINE FillCoarsePatch_SingleTime


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE FillPhysicalBC_uCF( pMF, sComp, nComp, Time, pGEOM ) BIND(c)

    ! --- No INTENT here because amrex source code doesn't have it ---

    TYPE(c_ptr),    VALUE :: pMF, pGEOM
    INTEGER(c_int), VALUE :: sComp, nComp
    REAL(DP),       VALUE :: Time

    TYPE(amrex_geometry) :: GEOM
    TYPE(amrex_multifab) :: MF
    TYPE(amrex_mfiter)   :: MFI
    INTEGER              :: pLo(4), pHi(4)
    REAL(DP), CONTIGUOUS, POINTER :: p(:,:,:,:)

    IF( .NOT. amrex_is_all_periodic() )THEN

      GEOM = pGEOM
      MF   = pMF

      !$OMP PARALLEL PRIVATE(MFI,p,pLo,pHi)
      CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )

      DO WHILE( MFI % next() )

        p => MF % DataPtr( MFI )

        ! Check if part of this box is outside the domain
        IF( .NOT. GEOM % DOMAIN % CONTAINS( p ) )THEN

          pLo = LBOUND( p )
          pHi = UBOUND( p )

          CALL amrex_filcc &
                 ( p, pLo, pHi, &
                   GEOM % DOMAIN % lo, GEOM % DOMAIN % hi, &
                   GEOM % dX, &
                   GEOM % get_physical_location( pLo ), &
                   lo_bc_uCF, hi_bc_uCF)

        END IF

      END DO
      !$OMP END PARALLEL

      CALL amrex_mfiter_destroy( MFI )

    END IF

  END SUBROUTINE FillPhysicalBC_uCF


  SUBROUTINE FillPhysicalBC_Dummy( pMF, sComp, nComp, Time, pGEOM ) BIND(c)

    ! --- No INTENT here because amrex source code doesn't have it ---

    TYPE(c_ptr),    VALUE :: pMF, pGEOM
    INTEGER(c_int), VALUE :: sComp, nComp
    REAL(DP),       VALUE :: Time

    TYPE(amrex_geometry) :: GEOM
    TYPE(amrex_multifab) :: MF
    TYPE(amrex_mfiter)   :: MFI
    INTEGER              :: pLo(4), pHi(4)
    REAL(DP), CONTIGUOUS, POINTER :: p(:,:,:,:)

    IF( .NOT. amrex_is_all_periodic() )THEN

      GEOM = pGEOM
      MF   = pMF

      !$OMP PARALLEL PRIVATE(MFI,p,pLo,pHi)
      CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )

      DO WHILE( MFI % next() )

        p => MF % DataPtr( MFI )

        ! Check if part of this box is outside the domain
        IF( .NOT. GEOM % DOMAIN % CONTAINS( p ) )THEN

          pLo = LBOUND( p )
          pHi = UBOUND( p )

          CALL amrex_filcc &
                 ( p, pLo, pHi, &
                   GEOM % DOMAIN % lo, GEOM % DOMAIN % hi, &
                   GEOM % dX, &
                   GEOM % get_physical_location( pLo ), &
                   lo_bc, hi_bc )

        END IF

      END DO
      !$OMP END PARALLEL

      CALL amrex_mfiter_destroy( MFI )

    END IF

  END SUBROUTINE FillPhysicalBC_Dummy


END MODULE FillPatchModule
