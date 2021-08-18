MODULE  MF_Euler_dgDiscretizationModule

  ! --- AMReX Modules ---

  USE amrex_amrcore_module, ONLY: &
    amrex_max_level, &
    amrex_geom
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum

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
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit, &
    OffGridFlux_Euler
  USE Euler_DiscontinuityDetectionModule, ONLY: &
    DetectShocks_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
  USE MF_FieldsModule, ONLY: &
    MF_OffGridFlux_Euler
  USE InputParsingModule, ONLY: &
    UseTiling, &
    swX
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap, &
    ApplyBoundaryConditions_Euler_MF
  USE FillPatchModule, ONLY: &
    FillPatch
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_InteriorBC, &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_Euler_MF


CONTAINS


  SUBROUTINE ComputeIncrement_Euler_MF( MF_uGF, MF_uCF, MF_uDF, MF_duCF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout)    :: MF_uCF (0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uDF (0:amrex_max_level)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCF(0:amrex_max_level)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: duCF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU(:,:,:,:,:)

    INTEGER :: iLevel, iCF
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    TYPE(EdgeMap) :: Edge_Map

    MF_OffGridFlux_Euler = Zero

    DO iLevel = 0, amrex_max_level

      ! --- Apply boundary conditions to interior domains ---

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

      CALL MF_uGF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )

      CALL MF_uDF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )

      CALL FillPatch( iLevel, 0.0e0_DP, MF_uCF(iLevel) ) ! replace 0.0 with Time

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF(iLevel) % DataPtr( MFI )
        uCF  => MF_uCF(iLevel) % DataPtr( MFI )
        uDF  => MF_uDF(iLevel) % DataPtr( MFI )

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

        CALL DetectShocks_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )

        CALL thornado2amrex_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( D  )

        DEALLOCATE( U  )

        DEALLOCATE( G  )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iLevel = 0, amrex_max_level

      ! --- Maybe don't need to apply boudnary conditions since
      !     they're applied in the shock detector ---

      ! --- Apply boundary conditions to interior domains ---

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

      CALL MF_uGF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )

      CALL MF_uDF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )

      CALL FillPatch( iLevel, 0.0e0_DP, MF_uCF(iLevel) ) ! replace 0.0 with Time

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

      CALL MF_duCF(iLevel) % SetVal( Zero )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCF  => MF_uCF (iLevel) % DataPtr( MFI )
        uDF  => MF_uDF (iLevel) % DataPtr( MFI )
        duCF => MF_duCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( G (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCF) )

        ALLOCATE( D (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nDF) )

        ALLOCATE( dU(1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCF) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )

        CALL amrex2thornado_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )

        ! --- Apply boundary conditions to physical boundaries ---

        CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

        CALL ApplyBoundaryConditions_Euler_MF &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

        CALL ComputeIncrement_Euler_DG_Explicit &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dU, &
                 SuppressBC_Option = .TRUE. )

        MF_OffGridFlux_Euler(iLevel,:) &
          = MF_OffGridFlux_Euler(iLevel,:) + OffGridFlux_Euler

        CALL thornado2amrex_X &
               ( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, duCF, dU )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( dU )

        DEALLOCATE( D  )

        DEALLOCATE( U  )

        DEALLOCATE( G  )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO

    DO iCF = 1, nCF

      CALL amrex_parallel_reduce_sum &
             ( MF_OffGridFlux_Euler(:,iCF), amrex_max_level+1 )

    END DO

  END SUBROUTINE ComputeIncrement_Euler_MF


END MODULE MF_Euler_dgDiscretizationModule
