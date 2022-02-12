MODULE  MF_Euler_dgDiscretizationModule

  ! --- AMReX Modules ---

  USE amrex_amrcore_module, ONLY: &
    amrex_max_level, &
    amrex_geom, &
    amrex_get_finest_level
  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_X1, &
    WeightsX_X2, &
    WeightsX_X3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nDF
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Psi, &
    iGF_SqrtGm
  USE Euler_dgDiscretizationModule, ONLY: &
    ComputeIncrement_Euler_DG_Explicit, &
    OffGridFlux_Euler
  USE Euler_DiscontinuityDetectionModule, ONLY: &
    DetectShocks_Euler
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE Euler_MeshRefinementModule, ONLY: &
    LX_X1_Refined_C, &
    LX_X2_Refined_C, &
    LX_X3_Refined_C

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    thornado2amrex_X_F, &
    amrex2thornado_X_F
  USE MF_FieldsModule, ONLY: &
    MF_OffGridFlux_Euler, &
    FluxRegister
  USE InputParsingModule, ONLY: &
    UseTiling, &
    swX, &
    do_reflux
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap, &
    ConstructEdgeMap, &
    ApplyBoundaryConditions_Euler_MF
  USE FillPatchModule, ONLY: &
    FillPatch_uGF, &
    FillPatch_uCF
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_InteriorBC, &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_Euler_MF


CONTAINS


  SUBROUTINE ComputeIncrement_Euler_MF &
    ( iLevel, Time, MF_uGF, MF_uCF, MF_uDF, MF_duCF, &
      UseXCFC_Option )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF
    TYPE(amrex_multifab), INTENT(in)    :: MF_uDF
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCF
    LOGICAL,              INTENT(in), OPTIONAL :: UseXCFC_Option

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF  (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF  (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uDF  (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: duCF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X1(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X2(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X3(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G   (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U   (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: D   (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: dU  (:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X1(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X2(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X3(:,:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    LOGICAL :: UseXCFC

    TYPE(amrex_multifab) :: SurfaceFluxes(1:nDimsX)
    INTEGER              :: iDimX, nGhost(nDimsX), nDOFX_X(3)
    LOGICAL              :: Nodal(nDimsX)

    TYPE(EdgeMap) :: Edge_Map

    UseXCFC = .FALSE.
    IF( PRESENT( UseXCFC_Option ) ) &
      UseXCFC = UseXCFC_Option

!    DO iLevel = 0, amrex_max_level
!
!      ! --- Apply boundary conditions to interior domains ---
!
!      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )
!
!      CALL MF_uGF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )
!
!      CALL MF_uCF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )
!
!      CALL MF_uDF(iLevel) % Fill_Boundary( amrex_geom(iLevel) )
!
!      CALL FillPatch_uGF( iLevel, Time(iLevel), MF_uGF(iLevel) )
!
!      CALL FillPatch_uCF( iLevel, Time(iLevel), MF_uCF(iLevel) )
!
!      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )
!
!      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )
!
!      DO WHILE( MFI % next() )
!
!        uGF  => MF_uGF(iLevel) % DataPtr( MFI )
!        uCF  => MF_uCF(iLevel) % DataPtr( MFI )
!        uDF  => MF_uDF(iLevel) % DataPtr( MFI )
!
!        iLo_MF = LBOUND( uGF )
!
!        BX = MFI % tilebox()
!
!        iX_B0 = BX % lo
!        iX_E0 = BX % hi
!        iX_B1 = BX % lo - swX
!        iX_E1 = BX % hi + swX
!
!        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )
!
!        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
!                            iX_B1(2):iX_E1(2), &
!                            iX_B1(3):iX_E1(3),1:nGF) )
!
!        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
!                            iX_B1(2):iX_E1(2), &
!                            iX_B1(3):iX_E1(3),1:nCF) )
!
!        ALLOCATE( D(1:nDOFX,iX_B1(1):iX_E1(1), &
!                            iX_B1(2):iX_E1(2), &
!                            iX_B1(3):iX_E1(3),1:nDF) )
!
!        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )
!
!        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )
!
!        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U )
!
!        CALL amrex2thornado_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )
!
!        ! --- Apply boundary conditions to physical boundaries ---
!
!        CALL ConstructEdgeMap( iLevel, BX, Edge_Map )
!
!        CALL ApplyBoundaryConditions_Euler_MF &
!               ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )
!
!        CALL DetectShocks_Euler( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D )
!
!        CALL thornado2amrex_X( nDF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uDF, D )
!
!        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )
!
!        DEALLOCATE( D  )
!
!        DEALLOCATE( U  )
!
!        DEALLOCATE( G  )
!
!        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )
!
!      END DO
!
!      CALL amrex_mfiter_destroy( MFI )
!
!    END DO

    ! --- Maybe don't need to apply boundary conditions since
    !     they're applied in the shock detector ---

    ! --- Apply boundary conditions to interior domains ---

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

!    CALL MF_uGF % Fill_Boundary( amrex_geom(iLevel) )
    CALL FillPatch_uGF( iLevel, Time, MF_uGF )

!    CALL MF_uCF % Fill_Boundary( amrex_geom(iLevel) )
    CALL FillPatch_uCF( iLevel, Time, MF_uCF )

    CALL MF_uDF % Fill_Boundary( amrex_geom(iLevel) )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InteriorBC )

    CALL MF_duCF % SetVal( Zero )

    CALL CreateMesh_MF( iLevel, MeshX )

    nDOFX_X(1) = nDOFX_X1
    nDOFX_X(2) = nDOFX_X2
    nDOFX_X(3) = nDOFX_X3

    nGhost = 0

    DO iDimX = 1, nDimsX

      Nodal        = .FALSE.
      Nodal(iDimX) = .TRUE.

      CALL amrex_multifab_build &
             ( SurfaceFluxes(iDimX), &
               MF_uGF % BA, MF_uGF % DM, &
               nDOFX_X(iDimX) * nCF, nGhost, Nodal )

      CALL SurfaceFluxes(iDimX) % SetVal( Zero )

    END DO

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF  => MF_uGF  % DataPtr( MFI )
      uCF  => MF_uCF  % DataPtr( MFI )
      uDF  => MF_uDF  % DataPtr( MFI )
      duCF => MF_duCF % DataPtr( MFI )

      uSurfaceFlux_X1 => SurfaceFluxes(1) % DataPtr( MFI )
      IF( nDimsX .GT. 1 ) uSurfaceFlux_X2 => SurfaceFluxes(2) % DataPtr( MFI )
      IF( nDimsX .GT. 2 ) uSurfaceFlux_X3 => SurfaceFluxes(3) % DataPtr( MFI )

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

      ALLOCATE( SurfaceFlux_X1(1:nDOFX_X1,iX_B0(1):iX_E0(1)+1, &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( SurfaceFlux_X2(1:nDOFX_X2,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2)+1, &
                                          iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( SurfaceFlux_X3(1:nDOFX_X3,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3)+1,1:nCF) )

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
               SuppressBC_Option = .TRUE., &
               SurfaceFlux_X1_Option = SurfaceFlux_X1, &
               SurfaceFlux_X2_Option = SurfaceFlux_X2, &
               SurfaceFlux_X3_Option = SurfaceFlux_X3 )

      CALL thornado2amrex_X &
             ( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, duCF, dU )

      iLo_MF = LBOUND( uSurfaceFlux_X1 )

      CALL thornado2amrex_X_F &
             ( nDOFX_X1, nCF, [ iX_B0(1)  , iX_B0(2), iX_B0(3) ], &
                              [ iX_E0(1)+1, iX_E0(2), iX_E0(3) ], iLo_MF, &
                              [ iX_B0(1)  , iX_B0(2), iX_B0(3) ], &
                              [ iX_E0(1)+1, iX_E0(2), iX_E0(3) ], &
                              uSurfaceFlux_X1, SurfaceFlux_X1 )

      iLo_MF = LBOUND( uSurfaceFlux_X2 )

      IF( nDimsX .GT. 1 ) &
        CALL thornado2amrex_X_F &
               ( nDOFX_X2, nCF, [ iX_B0(1), iX_B0(2)  , iX_B0(3) ], &
                                [ iX_E0(1), iX_E0(2)+1, iX_E0(3) ], iLo_MF, &
                                [ iX_B0(1), iX_B0(2)  , iX_B0(3) ], &
                                [ iX_E0(1), iX_E0(2)+1, iX_E0(3) ], &
                                uSurfaceFlux_X2, SurfaceFlux_X2 )

      iLo_MF = LBOUND( uSurfaceFlux_X3 )

      IF( nDimsX .GT. 2 ) &
        CALL thornado2amrex_X_F &
               ( nDOFX_X3, nCF, [ iX_B0(1), iX_B0(2), iX_B0(3)   ], &
                                [ iX_E0(1), iX_E0(2), iX_E0(3)+1 ], iLo_MF, &
                                [ iX_B0(1), iX_B0(2), iX_B0(3)   ], &
                                [ iX_E0(1), iX_E0(2), iX_E0(3)+1 ], &
                                uSurfaceFlux_X3, SurfaceFlux_X3 )

      MF_OffGridFlux_Euler(iLevel,:) &
        = MF_OffGridFlux_Euler(iLevel,:) + OffGridFlux_Euler

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      DEALLOCATE( SurfaceFlux_X3 )
      DEALLOCATE( SurfaceFlux_X2 )
      DEALLOCATE( SurfaceFlux_X1 )

      DEALLOCATE( dU )

      DEALLOCATE( D  )

      DEALLOCATE( U  )

      DEALLOCATE( G  )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

    END DO ! MFI

    CALL amrex_mfiter_destroy( MFI )

    IF( do_reflux )THEN

      IF( iLevel .GT. 0 )THEN

        CALL FluxRegister( iLevel ) &
               % FineAdd_DG &
                   ( SurfaceFluxes, +One, nCF,&
                     nDOFX_X1, nDOFX_X2, nDOFX_X3, &
                     WeightsX_X1, WeightsX_X2, WeightsX_X3, &
                     LX_X1_Refined_C, LX_X2_Refined_C, LX_X3_Refined_C )

      END IF

      IF( iLevel .LT. amrex_get_finest_level() )THEN

        CALL FluxRegister( iLevel+1 ) &
               % CrseInit_DG &
                   ( SurfaceFluxes, -One, nCF, &
                     nDOFX_X1, nDOFX_X2, nDOFX_X3, &
                     WeightsX_X1, WeightsX_X2, WeightsX_X3 )

      END IF

    END IF ! do_reflux

    DO iDimX = 1, nDimsX

      CALL amrex_multifab_destroy( SurfaceFluxes(iDimX) )

    END DO

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ComputeIncrement_Euler_MF


END MODULE MF_Euler_dgDiscretizationModule
