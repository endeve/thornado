!> Module to apply slope-limiter to AMReX MultiFabs.
!> @todo Fix issue of multiple grids giving different results.
MODULE MF_Euler_SlopeLimiterModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_communicator, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nDF
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE Euler_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
  USE InputParsingModule, ONLY: &
    nLevels, &
    swX, &
    UseSlopeLimiter, &
    UseTiling, &
    DEBUG
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap, &
    ApplyBoundaryConditions_Euler_MF
  USE FillPatchModule, ONLY: &
    FillPatch
!!$  USE AverageDownModule, ONLY: &
!!$    AverageDown
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_InteriorBC, &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplySlopeLimiter_Euler_MF

  INTERFACE ApplySlopeLimiter_Euler_MF
    MODULE PROCEDURE ApplySlopeLimiter_Euler_MF_MultipleLevels
    MODULE PROCEDURE ApplySlopeLimiter_Euler_MF_SingleLevel
  END INTERFACE ApplySlopeLimiter_Euler_MF

CONTAINS


  SUBROUTINE ApplySlopeLimiter_Euler_MF_MultipleLevels &
    ( Time, MF_uGF, MF_uCF, MF_uDF )

    REAL(DP),             INTENT(in)    :: Time  (0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)

    INTEGER :: iLevel, iErr

    DO iLevel = 0, nLevels-1

      IF( DEBUG )THEN

        CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )

        IF( amrex_parallel_ioprocessor() )THEN

          WRITE(*,'(2x,A,I3.3)') &
            'CALL ApplySlopeLimiter_Euler_MF_SingleLevel, iLevel: ', &
            iLevel

        END IF

      END IF

      CALL ApplySlopeLimiter_Euler_MF_SingleLevel &
             ( iLevel, Time(iLevel), MF_uGF, MF_uCF, MF_uDF )

    END DO

    ! --- Ensure underlying coarse cells are consistent with
    !     cells on refined level ---

!!$    CALL AverageDown( MF_uGF, MF_uCF )

  END SUBROUTINE ApplySlopeLimiter_Euler_MF_MultipleLevels


  SUBROUTINE ApplySlopeLimiter_Euler_MF_SingleLevel &
    ( iLevel, Time, MF_uGF, MF_uCF, MF_uDF )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uDF(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D(:,:,:,:,:)

    INTEGER       :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                     iLo_MF(4), iApplyBC(3)
    TYPE(EdgeMap) :: Edge_Map

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    ! --- Apply boundary conditions to interior domains ---

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

    CALL FillPatch( iLevel, Time, MF_uGF, MF_uGF )
    CALL FillPatch( iLevel, Time, MF_uGF, MF_uCF )
    CALL FillPatch( iLevel, Time, MF_uGF, MF_uDF )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF(iLevel) % DataPtr( MFI )
      uCF => MF_uCF(iLevel) % DataPtr( MFI )
      uDF => MF_uDF(iLevel) % DataPtr( MFI )

      iLo_MF = LBOUND( uGF )

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

      CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

      CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

      CALL amrex2thornado_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

      ! --- Apply boundary conditions to physical boundaries ---

      CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

      CALL ApplyBoundaryConditions_Euler_MF &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

      CALL Edge_Map % Euler_GetBC( iApplyBC )

      CALL ApplySlopeLimiter_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
               SuppressBC_Option = .TRUE., iApplyBC_Option = iApplyBC )

      CALL thornado2amrex_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, U )

      CALL thornado2amrex_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      DEALLOCATE( D )

      DEALLOCATE( U )

      DEALLOCATE( G )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

    END DO

    CALL amrex_mfiter_destroy( MFI )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ApplySlopeLimiter_Euler_MF_SingleLevel


END MODULE MF_Euler_SlopeLimiterModule
