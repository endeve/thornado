!> Module to apply slope-limiter to AMReX MultiFabs.
!> @todo Fix issue of multiple grids giving different results.
MODULE MF_Euler_SlopeLimiterModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,                 ONLY: &
    AR => amrex_real
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
  USE FluidFieldsModule,                 ONLY: &
    nCF, &
    nDF
  USE GeometryFieldsModule,              ONLY: &
    nGF
  USE Euler_SlopeLimiterModule,          ONLY: &
    ApplySlopeLimiter_Euler

  ! --- Local Modules ---

  USE MF_UtilitiesModule,                ONLY: &
    amrex2thornado_Euler, &
    thornado2amrex_Euler
  USE InputParsingModule,                ONLY: &
    nLevels,         &
    UseSlopeLimiter, &
    DEBUG
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap, &
    MF_ApplyBoundaryConditions_Euler
  USE TimersModule_AMReX_Euler,          ONLY: &
    TimersStart_AMReX_Euler,      &
    TimersStop_AMReX_Euler,       &
    Timer_AMReX_Euler_InteriorBC, &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ApplySlopeLimiter_Euler


CONTAINS


  SUBROUTINE MF_ApplySlopeLimiter_Euler( MF_uGF, MF_uCF, MF_uDF, GEOM )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uDF(:,:,:,:)

    REAL(AR), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iApplyBC(3)
    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uDF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uDF => MF_uDF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF) )

        ALLOCATE( D(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nDF) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_Euler( nGF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_Euler( nCF, iX_B1, iX_E1, uCF, U )

        CALL amrex2thornado_Euler( nDF, iX_B1, iX_E1, uDF, D )

        ! --- Apply boundary conditions to physical boundaries ---

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL MF_ApplyBoundaryConditions_Euler'

        CALL MF_ApplyBoundaryConditions_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL ApplySlopeLimiter_Euler'

        CALL Edge_Map % Euler_GetBC( iApplyBC )

        CALL ApplySlopeLimiter_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
                 SuppressBC_Option = .TRUE., iApplyBC_Option = iApplyBC )

        CALL thornado2amrex_Euler( nCF, iX_B1, iX_E1, uCF, U )

        CALL thornado2amrex_Euler( nDF, iX_B1, iX_E1, uDF, D )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( D )

        DEALLOCATE( U )

        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ApplySlopeLimiter_Euler


END MODULE MF_Euler_SlopeLimiterModule
