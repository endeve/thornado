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
    LX_X1_Dn_Refined, &
    LX_X1_Up_Refined, &
    LX_X2_Dn_Refined, &
    LX_X2_Up_Refined, &
    LX_X3_Dn_Refined, &
    LX_X3_Up_Refined

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
      FluxIncrement, UseXCFC_Option )

    INTEGER,              INTENT(in)    :: iLevel
    REAL(DP),             INTENT(in)    :: Time

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF
    TYPE(amrex_multifab), INTENT(in)    :: MF_uDF
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCF
    TYPE(amrex_multifab), INTENT(inout) :: FluxIncrement(2*nDimsX)
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

!        CALL ComputeFluxIncrement_Fine &
!               ( MF_uGF, UseXCFC, SurfaceFluxes, FluxIncrement )

        CALL FluxRegister( iLevel ) % FineAdd( FluxIncrement, +One )

      END IF

      IF( iLevel .LT. amrex_get_finest_level() )THEN

        CALL ComputeFluxIncrement_Coarse &
               ( MF_uGF, UseXCFC, SurfaceFluxes, FluxIncrement )

        CALL FluxRegister( iLevel+1 ) % CrseInit_DG( FluxIncrement, -One )

      END IF

    END IF ! do_reflux

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ComputeIncrement_Euler_MF


  !  --- PRIVATE SUBROUTINES ---


  SUBROUTINE ComputeFluxIncrement_Coarse &
    ( MF_uGF, UseXCFC, SurfaceFluxes, FluxIncrement )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF
    LOGICAL             , INTENT(in)    :: UseXCFC
    TYPE(amrex_multifab), INTENT(inout) :: SurfaceFluxes(nDimsX)
    TYPE(amrex_multifab), INTENT(inout) :: FluxIncrement(2*nDimsX)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X1(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X2(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uSurfaceFlux_X3(:,:,:,:)

    REAL(DP), CONTIGUOUS, POINTER :: uFluxIncrement_X1_L(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uFluxIncrement_X1_U(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uFluxIncrement_X2_L(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uFluxIncrement_X2_U(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uFluxIncrement_X3_L(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uFluxIncrement_X3_U(:,:,:,:)

    REAL(DP), ALLOCATABLE :: FluxIncrement_X1_L_P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X1_U_P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X2_L_P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X2_U_P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X3_L_P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X3_U_P(:,:,:,:,:)

    REAL(DP), ALLOCATABLE :: FluxIncrement_X1_L(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X1_U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X2_L(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X2_U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X3_L(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: FluxIncrement_X3_U(:,:,:,:,:)

    REAL(DP), ALLOCATABLE :: SurfaceFlux_X1(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X2(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X3(:,:,:,:,:)

    REAL(DP), ALLOCATABLE :: SurfaceFlux_X1_P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X2_P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: SurfaceFlux_X3_P(:,:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: tau(:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4), iDimX
    INTEGER :: iNX, iX1, iX2, iX3, iCF, nX_K

    CALL amrex_mfiter_build( MFI, MF_uGF, tiling = UseTiling )

    DO WHILE( MFI % next() )

      uGF => MF_uGF % DataPtr( MFI )

      uSurfaceFlux_X1     => SurfaceFluxes(1) % DataPtr( MFI )
      uFluxIncrement_X1_L => FluxIncrement(1) % DataPtr( MFI )
      uFluxIncrement_X1_U => FluxIncrement(2) % DataPtr( MFI )

      IF( nDimsX .GT. 1 )THEN

        uSurfaceFlux_X2     => SurfaceFluxes(2) % DataPtr( MFI )
        uFluxIncrement_X2_L => FluxIncrement(3) % DataPtr( MFI )
        uFluxIncrement_X2_U => FluxIncrement(4) % DataPtr( MFI )

      END IF

      IF( nDimsX .GT. 2 )THEN

        uSurfaceFlux_X3     => SurfaceFluxes(3) % DataPtr( MFI )
        uFluxIncrement_X3_L => FluxIncrement(5) % DataPtr( MFI )
        uFluxIncrement_X3_U => FluxIncrement(6) % DataPtr( MFI )

      END IF

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX

      CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      ALLOCATE( FluxIncrement_X1_L(1:nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( FluxIncrement_X1_U(1:nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( FluxIncrement_X2_L(1:nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( FluxIncrement_X2_U(1:nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( FluxIncrement_X3_L(1:nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( FluxIncrement_X3_U(1:nDOFX,iX_B0(1):iX_E0(1), &
                                           iX_B0(2):iX_E0(2), &
                                           iX_B0(3):iX_E0(3),1:nCF) )

      ALLOCATE( FluxIncrement_X1_L_P(1:nDOFX,1:nCF,iX_B0(2):iX_E0(2), &
                                                   iX_B0(3):iX_E0(3), &
                                                   iX_B0(1):iX_E0(1)) )
      ALLOCATE( FluxIncrement_X1_U_P(1:nDOFX,1:nCF,iX_B0(2):iX_E0(2), &
                                                   iX_B0(3):iX_E0(3), &
                                                   iX_B0(1):iX_E0(1)) )
      ALLOCATE( FluxIncrement_X2_L_P(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                                   iX_B0(3):iX_E0(3), &
                                                   iX_B0(2):iX_E0(2)) )
      ALLOCATE( FluxIncrement_X2_U_P(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                                   iX_B0(3):iX_E0(3), &
                                                   iX_B0(2):iX_E0(2)) )
      ALLOCATE( FluxIncrement_X3_L_P(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                                   iX_B0(2):iX_E0(2), &
                                                   iX_B0(3):iX_E0(3)) )
      ALLOCATE( FluxIncrement_X3_U_P(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                                   iX_B0(2):iX_E0(2), &
                                                   iX_B0(3):iX_E0(3)) )

      ALLOCATE( SurfaceFlux_X1(1:nDOFX_X1,iX_B0(1):iX_E0(1)+1, &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( SurfaceFlux_X2(1:nDOFX_X2,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2)+1, &
                                          iX_B0(3):iX_E0(3),1:nCF) )
      ALLOCATE( SurfaceFlux_X3(1:nDOFX_X3,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3)+1,1:nCF) )

      ALLOCATE( SurfaceFlux_X1_P(1:nDOFX_X1,1:nCF,iX_B0(2):iX_E0(2), &
                                                  iX_B0(3):iX_E0(3), &
                                                  iX_B0(1):iX_E0(1)+1) )
      ALLOCATE( SurfaceFlux_X2_P(1:nDOFX_X2,1:nCF,iX_B0(1):iX_E0(1), &
                                                  iX_B0(3):iX_E0(3), &
                                                  iX_B0(2):iX_E0(2)+1) )
      ALLOCATE( SurfaceFlux_X3_P(1:nDOFX_X3,1:nCF,iX_B0(1):iX_E0(1), &
                                                  iX_B0(2):iX_E0(2), &
                                                  iX_B0(3):iX_E0(3)+1) )

      ALLOCATE( tau(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3)) )
      ALLOCATE( G  (1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),nGF) )

      iLo_MF = LBOUND( uGF )

      CALL amrex2thornado_X &
             ( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

      iLo_MF = LBOUND( uSurfaceFlux_X1 )

      CALL amrex2thornado_X_F &
             ( nDOFX_X1, nCF, [ iX_B0(1)  , iX_B0(2), iX_B0(3) ], &
                              [ iX_E0(1)+1, iX_E0(2), iX_E0(3) ], iLo_MF, &
                              [ iX_B0(1)  , iX_B0(2), iX_B0(3) ], &
                              [ iX_E0(1)+1, iX_E0(2), iX_E0(3) ], &
                              uSurfaceFlux_X1, SurfaceFlux_X1 )

      iLo_MF = LBOUND( uSurfaceFlux_X2 )

      IF( nDimsX .GT. 1 ) &
        CALL amrex2thornado_X_F &
               ( nDOFX_X2, nCF, [ iX_B0(1), iX_B0(2)  , iX_B0(3) ], &
                                [ iX_E0(1), iX_E0(2)+1, iX_E0(3) ], iLo_MF, &
                                [ iX_B0(1), iX_B0(2)  , iX_B0(3) ], &
                                [ iX_E0(1), iX_E0(2)+1, iX_E0(3) ], &
                                uSurfaceFlux_X2, SurfaceFlux_X2 )

      iLo_MF = LBOUND( uSurfaceFlux_X3 )

      IF( nDimsX .GT. 2 ) &
        CALL amrex2thornado_X_F &
               ( nDOFX_X3, nCF, [ iX_B0(1), iX_B0(2), iX_B0(3)   ], &
                                [ iX_E0(1), iX_E0(2), iX_E0(3)+1 ], iLo_MF, &
                                [ iX_B0(1), iX_B0(2), iX_B0(3)   ], &
                                [ iX_E0(1), iX_E0(2), iX_E0(3)+1 ], &
                                uSurfaceFlux_X3, SurfaceFlux_X3 )

      ASSOCIATE( dX1 => MeshX(1) % Width, &
                 dX2 => MeshX(2) % Width, &
                 dX3 => MeshX(3) % Width )

      nX_K = PRODUCT( iX_E0 - iX_B0 + 1 )

      DO iX1 = iX_B0(1), iX_E0(1)+1
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iCF = 1, nCF
      DO iNX = 1, nDOFX_X1

        SurfaceFlux_X1_P(iNX,iCF,iX2,iX3,iX1) &
          = SurfaceFlux_X1(iNX,iX1,iX2,iX3,iCF) &
              * dX2(iX2) * dX3(iX3) * WeightsX_X1(iNX)

      END DO
      END DO
      END DO
      END DO
      END DO

      ! --- Contribution to Upper Element (X1) ---

      CALL MatrixMatrixMultiply &
             ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X1, + One, &
               LX_X1_Dn, nDOFX_X1, &
               SurfaceFlux_X1_P(1,1,iX_B0(2),iX_B0(3),iX_B0(1)  ), &
               nDOFX_X1, Zero, FluxIncrement_X1_U_P, nDOFX )

      ! --- Contribution to Lower Element (X1) ---

      CALL MatrixMatrixMultiply &
             ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X1, + One, &
               LX_X1_Up, nDOFX_X1, &
               SurfaceFlux_X1_P(1,1,iX_B0(2),iX_B0(3),iX_B0(1)+1), &
               nDOFX_X1, Zero, FluxIncrement_X1_L_P, nDOFX )

      IF( nDimsX .GT. 1 )THEN

        DO iX2 = iX_B0(2), iX_E0(2)+1
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iCF = 1, nCF
        DO iNX = 1, nDOFX_X2

          SurfaceFlux_X2_P(iNX,iCF,iX1,iX3,iX2) &
            = SurfaceFlux_X2(iNX,iX1,iX2,iX3,iCF) &
                * dX1(iX1) * dX3(iX3) * WeightsX_X2(iNX)

        END DO
        END DO
        END DO
        END DO
        END DO

        ! --- Contribution to Upper Element (X2) ---

        CALL MatrixMatrixMultiply &
               ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X2, + One, &
                 LX_X2_Dn, nDOFX_X2, &
                 SurfaceFlux_X2_P(1,1,iX_B0(1),iX_B0(3),iX_B0(2)  ), &
                 nDOFX_X2, Zero, FluxIncrement_X2_U_P, nDOFX )

        ! --- Contribution to Lower Element (X2) ---

        CALL MatrixMatrixMultiply &
               ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X2, + One, &
                 LX_X2_Up, nDOFX_X2, &
                 SurfaceFlux_X2_P(1,1,iX_B0(1),iX_B0(3),iX_B0(2)+1), &
                 nDOFX_X2, Zero, FluxIncrement_X2_L_P, nDOFX )

      END IF

      IF( nDimsX .GT. 2 )THEN

        DO iX3 = iX_B0(3), iX_E0(3)+1
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iCF = 1, nCF
        DO iNX = 1, nDOFX_X3

          SurfaceFlux_X3_P(iNX,iCF,iX1,iX2,iX3) &
            = SurfaceFlux_X3(iNX,iX1,iX2,iX3,iCF) &
                * dX1(iX1) * dX2(iX2) * WeightsX_X3(iNX)

        END DO
        END DO
        END DO
        END DO
        END DO

        ! --- Contribution to Upper Element (X3) ---

        CALL MatrixMatrixMultiply &
               ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X3, + One, &
                 LX_X3_Dn, nDOFX_X3, &
                 SurfaceFlux_X3_P(1,1,iX_B0(1),iX_B0(2),iX_B0(3)  ), &
                 nDOFX_X3, Zero, FluxIncrement_X3_U_P, nDOFX )

        ! --- Contribution to Lower Element (X3) ---

        CALL MatrixMatrixMultiply &
               ( 'T', 'N', nDOFX, nX_K*nCF, nDOFX_X3, + One, &
                 LX_X3_Up, nDOFX_X3, &
                 SurfaceFlux_X3_P(1,1,iX_B0(1),iX_B0(2),iX_B0(3)+1), &
                 nDOFX_X3, Zero, FluxIncrement_X3_L_P, nDOFX )

      END IF

      IF( UseXCFC )THEN

        DO iX3 = iX_B1(3), iX_E1(3)
        DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)
        DO iNX = 1, nDOFX

          tau(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Psi)**6

        END DO
        END DO
        END DO
        END DO

      ELSE

        DO iX3 = iX_B1(3), iX_E1(3)
        DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)
        DO iNX = 1, nDOFX

          tau(iNX,iX1,iX2,iX3) = One

        END DO
        END DO
        END DO
        END DO

      END IF

      ! --- Multiply with inverse mass matrix ---

      DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
      DO iX2 = iX_B0(2), iX_E0(2)
      DO iX1 = iX_B0(1), iX_E0(1)
      DO iNX = 1, nDOFX

        FluxIncrement_X1_U(iNX,iX1,iX2,iX3,iCF) &
          = FluxIncrement_X1_U_P(iNX,iCF,iX2,iX3,iX1) &
              * tau(iNX,iX1,iX2,iX3) &
              / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * dX1(iX1) * dX2(iX2) * dX3(iX3) )

        FluxIncrement_X1_L(iNX,iX1,iX2,iX3,iCF) &
          = FluxIncrement_X1_L_P(iNX,iCF,iX2,iX3,iX1) &
              * tau(iNX,iX1,iX2,iX3) &
              / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                    * dX1(iX1) * dX2(iX2) * dX3(iX3) )

      END DO
      END DO
      END DO
      END DO
      END DO

      IF( nDimsX .GT. 1 )THEN

        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          FluxIncrement_X2_U(iNX,iX1,iX2,iX3,iCF) &
            = FluxIncrement_X2_U_P(iNX,iCF,iX1,iX3,iX2) &
                * tau(iNX,iX1,iX2,iX3) &
                / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * dX1(iX1) * dX2(iX2) * dX3(iX3) )

          FluxIncrement_X2_L(iNX,iX1,iX2,iX3,iCF) &
            = FluxIncrement_X2_L_P(iNX,iCF,iX1,iX3,iX2) &
                * tau(iNX,iX1,iX2,iX3) &
                / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * dX1(iX1) * dX2(iX2) * dX3(iX3) )

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      IF( nDimsX .GT. 2 )THEN

        DO iCF = 1, nCF
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1, nDOFX

          FluxIncrement_X3_U(iNX,iX1,iX2,iX3,iCF) &
            = FluxIncrement_X3_U_P(iNX,iCF,iX1,iX2,iX3) &
                * tau(iNX,iX1,iX2,iX3) &
                / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * dX1(iX1) * dX2(iX2) * dX3(iX3) )

          FluxIncrement_X3_L(iNX,iX1,iX2,iX3,iCF) &
            = FluxIncrement_X3_L_P(iNX,iCF,iX1,iX2,iX3) &
                * tau(iNX,iX1,iX2,iX3) &
                / ( WeightsX_q(iNX) * G(iNX,iX1,iX2,iX3,iGF_SqrtGm) &
                      * dX1(iX1) * dX2(iX2) * dX3(iX3) )

        END DO
        END DO
        END DO
        END DO
        END DO

      END IF

      END ASSOCIATE ! dX1, dX2, dX3

      iLo_MF = LBOUND( uFluxIncrement_X1_L )

      CALL thornado2amrex_X &
             ( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, &
               uFluxIncrement_X1_U, FluxIncrement_X1_U )

      CALL thornado2amrex_X &
             ( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, &
               uFluxIncrement_X1_L, FluxIncrement_X1_L )

      IF( nDimsX .GT. 1 )THEN

        CALL thornado2amrex_X &
               ( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, &
                 uFluxIncrement_X2_L, FluxIncrement_X2_L )

        CALL thornado2amrex_X &
               ( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, &
                 uFluxIncrement_X2_L, FluxIncrement_X2_L )

      END IF

      IF( nDimsX .GT. 2 )THEN

        CALL thornado2amrex_X &
               ( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, &
                 uFluxIncrement_X3_L, FluxIncrement_X3_L )

        CALL thornado2amrex_X &
               ( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, &
                 uFluxIncrement_X3_L, FluxIncrement_X3_L )

      END IF

      DEALLOCATE( tau )
      DEALLOCATE( G )

      DEALLOCATE( SurfaceFlux_X3_P )
      DEALLOCATE( SurfaceFlux_X2_P )
      DEALLOCATE( SurfaceFlux_X1_P )

      DEALLOCATE( SurfaceFlux_X3 )
      DEALLOCATE( SurfaceFlux_X2 )
      DEALLOCATE( SurfaceFlux_X1 )

      DEALLOCATE( FluxIncrement_X3_U_P )
      DEALLOCATE( FluxIncrement_X3_L_P )
      DEALLOCATE( FluxIncrement_X2_U_P )
      DEALLOCATE( FluxIncrement_X2_L_P )
      DEALLOCATE( FluxIncrement_X1_U_P )
      DEALLOCATE( FluxIncrement_X1_L_P )

      DEALLOCATE( FluxIncrement_X3_U )
      DEALLOCATE( FluxIncrement_X3_L )
      DEALLOCATE( FluxIncrement_X2_U )
      DEALLOCATE( FluxIncrement_X2_L )
      DEALLOCATE( FluxIncrement_X1_U )
      DEALLOCATE( FluxIncrement_X1_L )

      CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

    END DO ! MFI

    CALL amrex_mfiter_destroy( MFI )

  END SUBROUTINE ComputeFluxIncrement_Coarse


END MODULE MF_Euler_dgDiscretizationModule
