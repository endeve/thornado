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
    amrex_interp_dg_order1, &
    amrex_interp_dg_order2, &
    amrex_interp_dg_order3
  USE amrex_fillpatch_module, ONLY: &
    amrex_fillpatch, &
    amrex_fillcoarsepatch
  USE amrex_geometry_module, ONLY: &
    amrex_geometry

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nNodes

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF
  USE InputParsingModule, ONLY: &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: FillPatch_uCF, FillCoarsePatch_uCF
  PUBLIC :: FillPatch_uGF, FillCoarsePatch_uGF


CONTAINS


  SUBROUTINE FillPatch_uGF( iLevel, Time, MF_uGF )

    USE MF_FieldsModule, ONLY: &
      MF_uGF_old, &
      MF_uGF_new
    USE InputParsingModule, ONLY: &
      t_old, &
      t_new
    USE AMReX_BoundaryConditionsModule, ONLY: &
      lo_bc, &
      hi_bc

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF

    INTEGER, PARAMETER :: sComp = 1, dComp = 1
    INTEGER :: nCompGF

    nCompGF = MF_uGF_old(iLevel) % nComp()

    IF( iLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF_uGF, t_old(iLevel), MF_uGF_old(iLevel), &
                                    t_new(iLevel), MF_uGF_new(iLevel), &
                            amrex_geom(iLevel), FillPhysicalBC, &
                            Time, sComp, dComp, nCompGF )

    ELSE

      IF( nNodes .EQ. 1 )THEN

        CALL amrex_fillpatch( MF_uGF, t_old(iLevel-1), MF_uGF_old(iLevel-1), &
                                      t_new(iLevel-1), MF_uGF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC, &
                                      t_old(iLevel  ), MF_uGF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_uGF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC, &
                              Time, sComp, dComp, nCompGF, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order1, &
                              lo_bc, hi_bc )

      ELSE IF( nNodes .EQ. 2 )THEN

        CALL amrex_fillpatch( MF_uGF, t_old(iLevel-1), MF_uGF_old(iLevel-1), &
                                      t_new(iLevel-1), MF_uGF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC, &
                                      t_old(iLevel  ), MF_uGF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_uGF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC, &
                              Time, sComp, dComp, nCompGF, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order2, &
                              lo_bc, hi_bc )

      ELSE IF( nNodes .EQ. 3 )THEN

        CALL amrex_fillpatch( MF_uGF, t_old(iLevel-1), MF_uGF_old(iLevel-1), &
                                      t_new(iLevel-1), MF_uGF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC, &
                                      t_old(iLevel  ), MF_uGF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_uGF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC, &
                              Time, sComp, dComp, nCompGF, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order3, &
                              lo_bc, hi_bc )

      ELSE

        CALL DescribeError_Euler_MF &
               ( 04, Message_Option = 'uGF', Int_Option = [ nNodes ] )

      END IF

    END IF

  END SUBROUTINE FillPatch_uGF


  SUBROUTINE FillCoarsePatch_uGF( iLevel, Time, MF_uGF )

    USE MF_FieldsModule, ONLY: &
      MF_uGF_old, &
      MF_uGF_new
    USE InputParsingModule, ONLY: &
      t_old, &
      t_new
    USE AMReX_BoundaryConditionsModule, ONLY: &
      lo_bc, &
      hi_bc

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF

    INTEGER :: nCompGF

    nCompGF = MF_uGF % nComp()

    IF( nNodes .EQ. 1 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF_uGF, t_old(iLevel-1), MF_uGF_old(iLevel-1), &
                       t_new(iLevel-1), MF_uGF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC, &
               amrex_geom(iLevel  ), FillPhysicalBC, &
               Time, 1, 1, nCompGF, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order1, lo_bc, hi_bc )

    ELSE IF( nNodes .EQ. 2 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF_uGF, t_old(iLevel-1), MF_uGF_old(iLevel-1), &
                       t_new(iLevel-1), MF_uGF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC, &
               amrex_geom(iLevel  ), FillPhysicalBC, &
               Time, 1, 1, nCompGF, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order2, lo_bc, hi_bc )

    ELSE IF( nNodes .EQ. 3 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF_uGF, t_old(iLevel-1), MF_uGF_old(iLevel-1), &
                       t_new(iLevel-1), MF_uGF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC, &
               amrex_geom(iLevel  ), FillPhysicalBC, &
               Time, 1, 1, nCompGF, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order3, lo_bc, hi_bc )

    ELSE

      CALL DescribeError_Euler_MF &
             ( 05, Message_Option = 'uGF', Int_Option = [ nNodes ] )

    END IF

  END SUBROUTINE FillCoarsePatch_uGF


  SUBROUTINE FillPatch_uCF( iLevel, Time, MF_uCF )

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

      CALL amrex_fillpatch( MF_uCF, t_old(iLevel), MF_uCF_old(iLevel), &
                                    t_new(iLevel), MF_uCF_new(iLevel), &
                            amrex_geom(iLevel), FillPhysicalBC, &
                            Time, sComp, dComp, nCompCF )

    ELSE

      IF( nNodes .EQ. 1 )THEN

        CALL amrex_fillpatch( MF_uCF, t_old(iLevel-1), MF_uCF_old(iLevel-1), &
                                      t_new(iLevel-1), MF_uCF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC, &
                                      t_old(iLevel  ), MF_uCF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_uCF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC, &
                              Time, sComp, dComp, nCompCF, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order1, &
                              lo_bc, hi_bc )

      ELSE IF( nNodes .EQ. 2 )THEN

        CALL amrex_fillpatch( MF_uCF, t_old(iLevel-1), MF_uCF_old(iLevel-1), &
                                      t_new(iLevel-1), MF_uCF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC, &
                                      t_old(iLevel  ), MF_uCF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_uCF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC, &
                              Time, sComp, dComp, nCompCF, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order2, &
                              lo_bc, hi_bc )

      ELSE IF( nNodes .EQ. 3 )THEN

        CALL amrex_fillpatch( MF_uCF, t_old(iLevel-1), MF_uCF_old(iLevel-1), &
                                      t_new(iLevel-1), MF_uCF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC, &
                                      t_old(iLevel  ), MF_uCF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_uCF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC, &
                              Time, sComp, dComp, nCompCF, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order3, &
                              lo_bc, hi_bc )

      ELSE

        CALL DescribeError_Euler_MF &
               ( 04, Message_Option = 'uCF', Int_Option = [ nNodes ] )

      END IF

    END IF

  END SUBROUTINE FillPatch_uCF


  SUBROUTINE FillCoarsePatch_uCF( iLevel, Time, MF_uCF )

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

    IF( nNodes .EQ. 1 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF_uCF, t_old(iLevel-1), MF_uCF_old(iLevel-1), &
                       t_new(iLevel-1), MF_uCF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC, &
               amrex_geom(iLevel  ), FillPhysicalBC, &
               Time, 1, 1, nCompCF, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order1, lo_bc, hi_bc )

    ELSE IF( nNodes .EQ. 2 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF_uCF, t_old(iLevel-1), MF_uCF_old(iLevel-1), &
                       t_new(iLevel-1), MF_uCF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC, &
               amrex_geom(iLevel  ), FillPhysicalBC, &
               Time, 1, 1, nCompCF, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order2, lo_bc, hi_bc )

    ELSE IF( nNodes .EQ. 3 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF_uCF, t_old(iLevel-1), MF_uCF_old(iLevel-1), &
                       t_new(iLevel-1), MF_uCF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC, &
               amrex_geom(iLevel  ), FillPhysicalBC, &
               Time, 1, 1, nCompCF, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order3, lo_bc, hi_bc )

    ELSE

      CALL DescribeError_Euler_MF &
             ( 05, Message_Option = 'uCF', Int_Option = [ nNodes ] )

    END IF

  END SUBROUTINE FillCoarsePatch_uCF


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
