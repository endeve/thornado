!> Module to apply slope-limiter to AMReX MultiFabs.
!> @todo Fix issue of multiple grids giving different results.
MODULE MF_MHD_SlopeLimiterModule

  ! --- AMReX Modules ---

  USE amrex_box_module,                  ONLY: &
    amrex_box
  USE amrex_geometry_module,             ONLY: &
    amrex_geometry
  USE amrex_multifab_module,             ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule,               ONLY: &
    swX, &
    nDOFX
  USE MagnetofluidFieldsModule,                 ONLY: &
    nCM, &
    nDM
  USE GeometryFieldsModule,              ONLY: &
    nGF
  USE MHD_SlopeLimiterModule,          ONLY: &
    ApplySlopeLimiter_MHD

  ! --- Local Modules ---

  USE MF_KindModule,                     ONLY: &
    DP
  USE MF_UtilitiesModule,                ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
  USE InputParsingModule,                ONLY: &
    nLevels,         &
    UseSlopeLimiter, &
    UseTiling,       &
    DEBUG
  USE MF_MHD_BoundaryConditionsModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap, &
    MF_ApplyBoundaryConditions_MHD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ApplySlopeLimiter_MHD


CONTAINS


  SUBROUTINE MF_ApplySlopeLimiter_MHD( MF_uGF, MF_uCM, MF_uDM, GEOM )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCM(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDM(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCM(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDM(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER       :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                     iLo_MF(4), iApplyBC(3)
    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCM(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uDM(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCM => MF_uCM(iLevel) % DataPtr( MFI )
        uDM => MF_uDM(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCM) )

        ALLOCATE( D(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nDM) )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCM, U )

        CALL amrex2thornado_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

        ! --- Apply boundary conditions to physical boundaries ---

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL MF_ApplyBoundaryConditions_MHD'

        CALL MF_ApplyBoundaryConditions_MHD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL ApplySlopeLimiter_MHD'

        CALL Edge_Map % MHD_GetBC( iApplyBC )

        CALL ApplySlopeLimiter_MHD &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
                 SuppressBC_Option = .TRUE., iApplyBC_Option = iApplyBC )

        CALL thornado2amrex_X( nCM, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCM, U )

        CALL thornado2amrex_X( nDM, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDM, D )

        DEALLOCATE( D )

        DEALLOCATE( U )

        DEALLOCATE( G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ApplySlopeLimiter_MHD


END MODULE MF_MHD_SlopeLimiterModule
