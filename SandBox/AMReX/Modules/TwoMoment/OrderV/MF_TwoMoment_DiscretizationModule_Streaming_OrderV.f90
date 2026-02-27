MODULE  MF_TwoMoment_DiscretizationModule_Streaming_OrderV

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab, &
    amrex_multifab_build, amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, amrex_mfiter_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_get_finest_level


  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX, nDOFZ, swE, nDOFE, iE_B0, iE_E0, iE_B1, iE_E1, nDimsX
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE GeometryFieldsModuleE,     ONLY: &
    nGE, uGE
  USE RadiationFieldsModule,            ONLY: &
    nCR, nSpecies
  USE FluidFieldsModule,            ONLY: &
    nCF
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit, &
    OffGridFlux_TwoMoment  
  USE MeshModule, ONLY: &
    MeshX
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3


  ! --- Local Modules ---
  USE MF_KindModule, ONLY: &
    Zero
  USE MF_UtilitiesModule,                ONLY: &
    amrex2thornado_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z
  USE InputParsingModule,                       ONLY: &
    nLevels, &
    UseTiling, &
    nE, &
    UseFluxCorrection_TwoMoment
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE FillPatchModule, ONLY: &
    FillPatch
  USE AverageDownModule, ONLY: &
    AverageDown
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_Permute, &
    OffGridFlux_TwoMoment_MF
  USE MF_UtilitiesModule, ONLY: &
    amrex2amrex_permute_Z, &
    amrex_permute2amrex_Z, &
    MF_amrex2amrex_permute_Z_Level, &
    MF_amrex_permute2amrex_Z_Level
  USE MF_FieldsModule_TwoMoment, ONLY: &
    FluxRegister_TwoMoment
  USE Euler_MeshRefinementModule, ONLY: &
    FaceRatio, &
    WeightsX_X1c, &
    WeightsX_X2c, &
    WeightsX_X3c, &
    nFineX_X1, &
    nFineX_X2, &
    nFineX_X3, &
    vpLX_X1_Refined, &
    vpLX_X2_Refined, &
    vpLX_X3_Refined



  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit_MF


CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Explicit_MF( Time, GEOM, MF_uGF, MF_uCF, MF_uCR, MF_duCR, Verbose_Option )


    REAL(amrex_real)    , INTENT(in)    :: Time(0:)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout)    :: MF_uCF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCR(0:nLevels-1)
    LOGICAL,  INTENT(in), OPTIONAL      :: Verbose_Option

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCR (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: duCR(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: C (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U (:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: dU(:,:,:,:,:,:,:)


    REAL(amrex_real), CONTIGUOUS, POINTER :: uSurfaceFlux_X1(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uSurfaceFlux_X2(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uSurfaceFlux_X3(:,:,:,:)


    REAL(amrex_real), ALLOCATABLE :: SurfaceFlux_X1(:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: SurfaceFlux_X2(:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: SurfaceFlux_X3(:,:,:,:,:,:,:)


    TYPE(amrex_multifab) :: SurfaceFluxes(1:nDimsX)
    INTEGER              :: iDimX, nGhost(nDimsX), nDOFX_X(3)
    LOGICAL              :: Nodal(3)
    INTEGER              :: nComp_Flux, iFd, iNodeE, iNodeZ, iS, nFields, nS, nDOFZ_X

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), i, j
    INTEGER :: iZ1, iZ2, iZ3, iZ4, iCR, iNodeX, iComp, iField

    LOGICAL :: Verbose

    TYPE(EdgeMap) :: Edge_Map

    Verbose = .TRUE.

    OffGridFlux_TwoMoment_MF = Zero

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL FillPatch &
           ( iLevel, MF_uGF,MF_uCR)
      
      CALL FillPatch( iLevel, MF_uGF, MF_uCF )


      CALL FillPatch &
           ( iLevel, MF_uGF, &
             ApplyBoundaryConditions_Geometry_Option = .TRUE. )

      CALL MF_duCR(iLevel) % setval( 0.0_amrex_real )

      CALL CreateMesh_MF( iLevel, MeshX )

      
      nDOFX_X(1) = nDOFX_X1
      nDOFX_X(2) = nDOFX_X2
      nDOFX_X(3) = nDOFX_X3

      nGhost = 0

      nComp_Flux = ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies * nDOFE * nDOFX_X1

      DO iDimX = 1, nDimsX

        Nodal        = .FALSE.
        Nodal(iDimX) = .TRUE.

        CALL amrex_multifab_build &
              ( SurfaceFluxes(iDimX), &
                MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
                nComp_Flux, nGhost, Nodal )

        CALL SurfaceFluxes(iDimX) % SetVal( Zero )

      END DO

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCF  => MF_uCF (iLevel) % DataPtr( MFI )
        uCR  => MF_uCR (iLevel) % DataPtr( MFI )
        duCR => MF_duCR(iLevel) % DataPtr( MFI )
        uSurfaceFlux_X1 => SurfaceFluxes(1) % DataPtr( MFI )
        IF( nDimsX .GT. 1 ) uSurfaceFlux_X2 => SurfaceFluxes(2) % DataPtr( MFI )
        IF( nDimsX .GT. 2 ) uSurfaceFlux_X3 => SurfaceFluxes(3) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX



        iZ_B0(1)=iE_B0
        iZ_E0(1)=iE_E0
        iZ_B1(1)=iE_B1
        iZ_E1(1)=iE_E1


        iZ_B0(2:4)=iX_B0(1:3)
        iZ_E0(2:4)=iX_E0(1:3)
        iZ_B1(2:4)=iX_B1(1:3)
        iZ_E1(2:4)=iX_E1(1:3)


        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 C )

        CALL AllocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nCR     , &
                   nSpecies ], &
                 U )

        CALL AllocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nCR     , &
                   nSpecies ], &
                 dU )

! SurfaceFlux_X1: faces in X1 direction
CALL AllocateArray_Z &
       ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
         [ nDOFE * nDOFX_X1, iZ_E0(1), iZ_E0(2)+1, iZ_E0(3), iZ_E0(4), nCR, nSpecies ], &
         SurfaceFlux_X1 )

! SurfaceFlux_X2: faces in X2 direction
CALL AllocateArray_Z &
       ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
         [ nDOFX_X2, iZ_E0(1), iZ_E0(2), iZ_E0(3)+1, iZ_E0(4), nComp_Flux, nSpecies ], &
         SurfaceFlux_X2 )

! SurfaceFlux_X3: faces in X3 direction
CALL AllocateArray_Z &
       ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
         [ nDOFX_X3, iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4)+1, nComp_Flux, nSpecies ], &
         SurfaceFlux_X3 )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, C )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B1, iZ_E1, uCR, U )

        CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

        CALL ApplyBoundaryConditions_TwoMoment_MF &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )

        CALL ApplyBoundaryConditions_Euler_MF &
               ( iX_B0, iX_E0, iX_B1, iX_E1, C, Edge_Map )

       CALL ComputeIncrement_TwoMoment_Explicit &
              ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, G, C, U, dU, SuppressBC_Option = .TRUE., SurfaceFlux_X1_Option = SurfaceFlux_X1 )

 DO i=1,nCR
        OffGridFlux_TwoMoment_MF(i,iLevel) &
          = OffGridFlux_TwoMoment_MF(i,iLevel) &
          + OffGridFlux_TwoMoment(i)  
        OffGridFlux_TwoMoment_MF(i+nCR,iLevel) &
          = OffGridFlux_TwoMoment_MF(i+nCR,iLevel) &
          + OffGridFlux_TwoMoment(i+nCR)  
END DO

        CALL thornado2amrex_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, duCR, dU )


iLo_MF = LBOUND( uSurfaceFlux_X1 )

!PRINT *, '=== DEBUG Surface Flux X1 ==='
!PRINT *, 'iZ_B0 = ', iZ_B0
!PRINT *, 'iZ_E0 = ', iZ_E0
!PRINT *, 'iLo_MF = ', iLo_MF
!PRINT *, 'SHAPE(uSurfaceFlux_X1) = ', SHAPE(uSurfaceFlux_X1)
!PRINT *, 'LBOUND(uSurfaceFlux_X1) = ', LBOUND(uSurfaceFlux_X1)
!PRINT *, 'UBOUND(uSurfaceFlux_X1) = ', UBOUND(uSurfaceFlux_X1)
!PRINT *, 'Loop X1 faces: ', iZ_B0(2), ' to ', iZ_E0(2)+1
!PRINT *, 'nDOFX_X1 = ', nDOFX_X1
!PRINT *, 'nE = ', iZ_E0(1) - iZ_B0(1) + 1
!PRINT *, 'nCR = ', nCR
!PRINT *, 'nSpecies = ', nSpecies
!PRINT *, 'nComp_Flux = ', nComp_Flux
!PRINT *, '=== END DEBUG ==='


!DO iS = 1, nSpecies
!DO iCR = 1, nCR
!DO iZ4 = iZ_B0(4), iZ_E0(4)
!DO iZ3 = iZ_B0(3), iZ_E0(3)
!DO iZ2 = iZ_B0(2), iZ_E0(2)+1
!DO iZ1 = iZ_B0(1), iZ_E0(1)
!DO iNodeX = 1, nDOFZ

!  iComp = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFX_X1 &
!        + ( iCR - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFX_X1 &
!        + ( iZ1 - iE_B0 ) * nDOFX_X1 &
!        + iNodeX

!uSurfaceFlux_X1( iZ2, &       ! NO +1 (LBOUND=0, iZ_B0(2)=0)
!                 iZ3 - iZ_B0(3) + 1, &   ! +1 (LBOUND=1, iZ_B0(3)=1)
!                 iZ4 - iZ_B0(4) + 1, &   ! +1 (LBOUND=1, iZ_B0(4)=1)
!                 iComp ) &
!  = SurfaceFlux_X1( iNodeX, iZ1, iZ2, iZ3, iZ4, iCR, iS )


!END DO
!END DO
!END DO
!END DO
!END DO
!END DO
!END DO

!DO iS = 1, nSpecies
!DO iCR = 1, nCR 
!DO iZ4 = iZ_B0(4), iZ_E0(4)
!DO iZ3 = iZ_B0(3), iZ_E0(3)
!DO iZ2 = iZ_B0(2), iZ_E0(2)+1
!DO iZ1 = iZ_B0(1), iZ_E0(1)
!DO iNodeE = 1, nDOFE
!DO iNodeX = 1, nDOFX_X1

DO iS = 1, nSpecies
DO iZ1 = iZ_B0(1), iZ_E0(1)
DO iNodeE = 1, nDOFE
DO iCR = 1, nCR                 
DO iNodeX = 1, nDOFX_X1
DO iZ4 = iZ_B0(4), iZ_E0(4)
DO iZ3 = iZ_B0(3), iZ_E0(3)
DO iZ2 = iZ_B0(2), iZ_E0(2)+1



  iNodeZ = iNodeE + ( iNodeX - 1 ) * nDOFE

  iField = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFE &
         + ( iCR - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFE &
         + ( iZ1 - iE_B0 ) * nDOFE &
         + iNodeE
  
  iComp = ( iField - 1 ) * nDOFX_X1 + iNodeX

  uSurfaceFlux_X1( iZ2, iZ3, iZ4, iComp ) &
    = SurfaceFlux_X1( iNodeZ, iZ1, iZ2, iZ3, iZ4, iCR, iS )

END DO
END DO
END DO
END DO
END DO
END DO
END DO
END DO


!PRINT *, '=== Surface Flux Debug ==='
!PRINT *, 'iLevel = ', iLevel
!PRINT *, 'SurfaceFlux_X1 min/max = ', MINVAL(SurfaceFlux_X1), MAXVAL(SurfaceFlux_X1)
!PRINT *, 'uSurfaceFlux_X1 min/max = ', MINVAL(uSurfaceFlux_X1), MAXVAL(uSurfaceFlux_X1)
!PRINT *, 'SurfaceFluxes(1) ncomp = ', SurfaceFluxes(1) % ncomp()
!PRINT *, '=========================='


CALL DeallocateArray_Z &
       ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
         [ nDOFX_X3, iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4), nCR, nSpecies ], &
         SurfaceFlux_X3 )

CALL DeallocateArray_Z &
       ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
         [ nDOFX_X2, iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4), nCR, nSpecies ], &
         SurfaceFlux_X2 )

CALL DeallocateArray_Z &
       ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
         [ nDOFE * nDOFX_X1, iZ_E0(1), iZ_E0(2)+1, iZ_E0(3), iZ_E0(4), nCR, nSpecies ], &
         SurfaceFlux_X1 )


        CALL DeallocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nCR     , &
                   nSpecies ], &
                 dU )

        CALL DeallocateArray_Z &
               ( [ 1       , &
                   iZ_B1(1), &
                   iZ_B1(2), &
                   iZ_B1(3), &
                   iZ_B1(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E1(1), &
                   iZ_E1(2), &
                   iZ_E1(3), &
                   iZ_E1(4), &
                   nCR     , &
                   nSpecies ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 C )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL amrex_parallel_reduce_sum( OffGridFlux_TwoMoment_MF(:,iLevel), 2 * nCR )


IF( UseFluxCorrection_TwoMoment )THEN


  !PRINT *, '=== Flux Register Debug ==='
  !PRINT *, 'iLevel = ', iLevel
  !PRINT *, 'SurfaceFluxes(1) % ncomp() = ', SurfaceFluxes(1) % ncomp()
  !PRINT *, 'nDOFX_X1 * nComp_Flux = ', nDOFX_X1 * nComp_Flux
  !PRINT *, '==========================='

  IF( iLevel .GT. 0 ) THEN
    !PRINT *, 'Calling FineAdd_DG at level ', iLevel
    CALL FluxRegister_TwoMoment(iLevel) &
           % FineAdd_DG &
              ( SurfaceFluxes, nComp_Flux, FaceRatio, &
                nDOFX_X1, nDOFX_X2, nDOFX_X3, &
                nFineX_X1, nFineX_X2, nFineX_X3, &
                WeightsX_X1c, WeightsX_X2c, WeightsX_X3c, &
                vpLX_X1_Refined, vpLX_X2_Refined, vpLX_X3_Refined )
  END IF
  IF( iLevel .LT. amrex_get_finest_level() ) THEN
  !PRINT *, '=== Before CrseInit_DG ==='
  !PRINT *, 'Accessing FluxRegister at level = ', iLevel+1
  !PRINT *, 'finest_level = ', amrex_get_finest_level()
  !PRINT *, 'SurfaceFluxes ncomp = ', SurfaceFluxes(1) % ncomp()
  !PRINT *, 'nDOFX_X1, nDOFX_X2, nDOFX_X3', nDOFX_X1, nDOFX_X2, nDOFX_X3
  !PRINT *, '==========================='
    CALL FluxRegister_TwoMoment(iLevel+1) &
           % CrseInit_DG &
               ( SurfaceFluxes, nComp_Flux, &
                 nDOFX_X1, nDOFX_X2, nDOFX_X3, &
                 WeightsX_X1c, WeightsX_X2c, WeightsX_X3c )

END IF

END IF


DO iDimX = 1, nDimsX
  CALL amrex_multifab_destroy( SurfaceFluxes(iDimX) )
END DO

      CALL DestroyMesh_MF( MeshX )

    END DO

    DO i = 0, nLevels-1
        CALL amrex_multifab_build &
               ( MF_Permute(i), MF_uGF(i) % BA, &
                 MF_uGF(i) % DM, nDOFX * nDOFE * ( iE_E0 - iE_B0 + 1 ) * nCR * nSpecies, swX )


      CALL MF_amrex2amrex_permute_Z_Level(i,nCR,MF_uGF(i),MF_uCR(i),MF_Permute(i))

    END DO


    CALL AverageDown( MF_uGF, MF_Permute )

    DO i = 0, nLevels-1

      CALL MF_amrex_permute2amrex_Z_Level(i,nCR,MF_uGF(i),MF_uCR(i),MF_Permute(i))

      CALL amrex_multifab_destroy( MF_Permute(i) )


    END DO
    


  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit_MF


END MODULE  MF_TwoMoment_DiscretizationModule_Streaming_OrderV
