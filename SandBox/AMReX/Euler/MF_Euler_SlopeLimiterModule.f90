!> Module to apply slope-limiter to AMReX MultiFabs.
!> @todo Fix issue of multiple grids giving different results.
MODULE MF_Euler_SlopeLimiterModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX
  USE FluidFieldsModule,        ONLY: &
    nCF
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE Euler_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_Euler

  ! --- Local Modules ---
  USE MF_UtilitiesModule,                ONLY: &
    AMReX2thornado, &
    thornado2AMReX
  USE MyAmrModule,                       ONLY: &
    nLevels, UseSlopeLimiter, DEBUG
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap, ConstructEdgeMap, &
    MF_ApplyBoundaryConditions_Euler
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler, TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_InteriorBC, &
    Timer_AMReX_Euler_DataTransfer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ApplySlopeLimiter_Euler


CONTAINS


  SUBROUTINE MF_ApplySlopeLimiter_Euler( MF_uGF, MF_uCF, GEOM )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    DO iLevel = 0, nLevels

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )
      ! --- Apply boundary conditions to interior domains ---
      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )
      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )
        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )
        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF) )

        CALL AMReX2thornado &
               ( nGF, iX_B1, iX_E1, &
                 uGF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nGF), &
                 G(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nGF) )

        CALL AMReX2thornado &
               ( nCF, iX_B1, iX_E1, &
                 uCF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nCF), &
                 U(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nCF) )
        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

        ! --- Apply boundary conditions to physical boundaries ---
        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL MF_ApplyBoundaryConditions_Euler'
        CALL MF_ApplyBoundaryConditions_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1,  &
                 U(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nCF), Edge_Map )

        IF( DEBUG ) WRITE(*,'(A)') '    CALL ApplySlopeLimiter_Euler'
        CALL ApplySlopeLimiter_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G (1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF), &
                 U (1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF), &
                 SuppressBC_Option = .TRUE. )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )
        CALL thornado2AMReX &
               ( nCF, iX_B1, iX_E1, &
                 uCF(      iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nDOFX*nCF), &
                 U(1:nDOFX,iX_B1(1):iX_E1(1), &
                           iX_B1(2):iX_E1(2), &
                           iX_B1(3):iX_E1(3),1:nCF) )

        DEALLOCATE( U )
        DEALLOCATE( G )
        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ApplySlopeLimiter_Euler


END MODULE MF_Euler_SlopeLimiterModule
