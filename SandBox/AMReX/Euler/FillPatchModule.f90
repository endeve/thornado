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
    amrex_interp_dg_order1, &
    amrex_interp_dg_order2, &
    amrex_interp_dg_order3
  USE amrex_amrcore_module, ONLY: &
    amrex_max_level
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
  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF
  USE InputParsingModule, ONLY: &
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


CONTAINS


  SUBROUTINE FillPatch( iLevel, Time, MF_old, MF_new, MF )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(in)    :: MF_old(0:amrex_max_level), &
                                           MF_new(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF

    INTEGER, PARAMETER :: sComp = 1, dComp = 1
    INTEGER :: nComp

    nComp = MF % nComp()

    IF( iLevel .EQ. 0 )THEN

      CALL amrex_fillpatch( MF, t_old(iLevel), MF_old(iLevel), &
                                t_new(iLevel), MF_new(iLevel), &
                            amrex_geom(iLevel), FillPhysicalBC_Dummy, &
                            Time, sComp, dComp, nComp )

    ELSE

      IF( nNodes .EQ. 1 )THEN

        CALL amrex_fillpatch( MF, t_old(iLevel-1), MF_old(iLevel-1), &
                                  t_new(iLevel-1), MF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
                                      t_old(iLevel  ), MF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
                              Time, sComp, dComp, nComp, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order1, &
                              lo_bc, hi_bc )

      ELSE IF( nNodes .EQ. 2 )THEN

        CALL amrex_fillpatch( MF, t_old(iLevel-1), MF_old(iLevel-1), &
                                  t_new(iLevel-1), MF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
                                      t_old(iLevel  ), MF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
                              Time, sComp, dComp, nComp, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order2, &
                              lo_bc, hi_bc )

      ELSE IF( nNodes .EQ. 3 )THEN

        CALL amrex_fillpatch( MF, t_old(iLevel-1), MF_old(iLevel-1), &
                                  t_new(iLevel-1), MF_new(iLevel-1), &
                              amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
                                      t_old(iLevel  ), MF_old(iLevel  ), &
                                      t_new(iLevel  ), MF_new(iLevel  ), &
                              amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
                              Time, sComp, dComp, nComp, &
                              amrex_ref_ratio(iLevel-1), &
                              amrex_interp_dg_order3, &
                              lo_bc, hi_bc )

      ELSE

        CALL DescribeError_Euler_MF &
               ( 04, Message_Option = 'MF', Int_Option = [ nNodes ] )

      END IF

    END IF

  END SUBROUTINE FillPatch


  SUBROUTINE FillCoarsePatch( iLevel, Time, MF_old, MF_new, MF )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(in)    :: MF_old(0:amrex_max_level), &
                                           MF_new(0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF

    INTEGER :: nComp

    nComp = MF % nComp()

    IF( nNodes .EQ. 1 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF, t_old(iLevel-1), MF_old(iLevel-1), &
                   t_new(iLevel-1), MF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
               amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
               Time, 1, 1, nComp, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order1, lo_bc, hi_bc )

    ELSE IF( nNodes .EQ. 2 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF, t_old(iLevel-1), MF_old(iLevel-1), &
                   t_new(iLevel-1), MF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
               amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
               Time, 1, 1, nComp, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order2, lo_bc, hi_bc )

    ELSE IF( nNodes .EQ. 3 )THEN

      CALL amrex_fillcoarsepatch &
             ( MF, t_old(iLevel-1), MF_old(iLevel-1), &
                   t_new(iLevel-1), MF_new(iLevel-1), &
               amrex_geom(iLevel-1), FillPhysicalBC_Dummy, &
               amrex_geom(iLevel  ), FillPhysicalBC_Dummy, &
               Time, 1, 1, nComp, amrex_ref_ratio(iLevel-1), &
               amrex_interp_dg_order3, lo_bc, hi_bc )

    ELSE

      CALL DescribeError_Euler_MF &
             ( 05, Message_Option = 'MF', Int_Option = [ nNodes ] )

    END IF

  END SUBROUTINE FillCoarsePatch


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
