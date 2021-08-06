MODULE FillPatchModule

  USE ISO_C_BINDING

  ! --- AMReX Modules ---

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
    amrex_geometry

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE InputParsingModule, ONLY: &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FillPatch
  PUBLIC :: FillCoarsePatch


CONTAINS


  SUBROUTINE FillPatch( iLevel, Time, MF_uCF )

    USE MF_FieldsModule, ONLY: &
      MF_uCF_old, &
      MF_uCF_new
    USE InputParsingModule, ONLY: &
      t_old, &
      t_new
    USE AMReX_BoundaryConditionsModule, ONLY: &
      lo_bc, &
      hi_bc

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    INTEGER, PARAMETER :: sComp = 1, dComp = 1
    INTEGER :: nCompCF

    nCompCF = MF_uCF_old(iLevel) % nComp()

    IF( iLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF_uCF, t_old (iLevel), MF_uCF_old(iLevel), &
                                    t_new (iLevel), MF_uCF_new(iLevel), &
                            amrex_geom(iLevel), FillPhysicalBC, &
                            Time, sComp, dComp, nCompCF )

    ELSE

      CALL amrex_fillpatch( MF_uCF, t_old (iLevel-1), MF_uCF_old(iLevel-1), &
                                    t_new (iLevel-1), MF_uCF_new(iLevel-1), &
                            amrex_geom(iLevel-1), FillPhysicalBC, &
                                t_old (iLevel  ), MF_uCF_old(iLevel  ), &
                                t_new (iLevel  ), MF_uCF_new(iLevel  ), &
                            amrex_geom(iLevel  ), FillPhysicalBC, &
                            Time, sComp, dComp, nCompCF, &
                            amrex_ref_ratio(iLevel-1), amrex_interp_dg, &
                            lo_bc, hi_bc )

    END IF

  END SUBROUTINE FillPatch


  SUBROUTINE FillCoarsePatch( iLevel, Time, MF_uCF )

    USE MF_FieldsModule, ONLY: &
      MF_uCF_old, &
      MF_uCF_new
    USE InputParsingModule, ONLY: &
      t_old, &
      t_new
    USE AMReX_BoundaryConditionsModule, ONLY: &
      lo_bc, &
      hi_bc

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF

    INTEGER :: nCompCF

    nCompCF = MF_uCF % nComp()

    CALL amrex_fillcoarsepatch &
           ( MF_uCF, t_old(iLevel-1), MF_uCF_old(iLevel-1), &
                     t_new(iLevel-1), MF_uCF_new(iLevel-1), &
             amrex_geom(iLevel-1), FillPhysicalBC, &
             amrex_geom(iLevel  ), FillPhysicalBC, &
             Time, 1, 1, nCompCF, amrex_ref_ratio(iLevel-1), &
             amrex_interp_dg, lo_bc, hi_bc )

  END SUBROUTINE FillCoarsePatch


  SUBROUTINE FillPhysicalBC( pMF, sComp, nComp, Time, pGEOM ) BIND(c)

    USE amrex_geometry_module, ONLY: &
      amrex_is_all_periodic
    USE amrex_filcc_module, ONLY: &
      amrex_filcc
    USE AMReX_BoundaryConditionsModule, ONLY: &
      lo_bc, &
      hi_bc

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

        IF( .NOT. GEOM % DOMAIN % CONTAINS( p ) )THEN ! part of this box is outside the domain

          pLo = LBOUND( p )
          pHi = UBOUND( p )

          CALL amrex_filcc &
                 ( p, pLo, pHi, & ! fortran array and bounds
                   GEOM % DOMAIN % lo, GEOM % DOMAIN % hi, & ! index extent of whole problem domain
                   GEOM % dX, & ! cell size in real
                   GEOM % get_physical_location( pLo ), & ! physical location of lower left corner
                   lo_bc, hi_bc) ! bc types for each component

             ! amrex_filcc doesn't fill EXT_DIR (see amrex_bc_types_module for a list of bc types
             ! In that case, the user needs to fill it.
        END IF

      END DO
      !$OMP END PARALLEL

      CALL amrex_mfiter_destroy( MFI )

    END IF

  END SUBROUTINE FillPhysicalBC

END MODULE FillPatchModule
