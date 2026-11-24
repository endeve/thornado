MODULE  MF_TwoMoment_DiscretizationModule_Streaming_OrderV

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum, amrex_parallel_ioprocessor
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, amrex_imultifab, &
    amrex_multifab_build, amrex_multifab_destroy, &
    amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy
  USE MaskModule, ONLY: &
    CreateFineMask, DestroyFineMask, IsNotLeafElement
  USE amrex_amrcore_module, ONLY: &
    amrex_get_finest_level, &
    amrex_geom


  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX, nDOFZ, swE, nDOFE, iE_B0, iE_E0, iE_B1, iE_E1, nDimsX
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE GeometryFieldsModuleE,     ONLY: &
    nGE, uGE, iGE_Ep2, iGE_Ep3
  USE RadiationFieldsModule,            ONLY: &
    nCR, nSpecies, LeptonNumber, iCR_N
  USE FluidFieldsModule,            ONLY: &
    nCF
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit, &
    OffGridFlux_TwoMoment  
  USE MeshModule, ONLY: &
    MeshX, MeshE
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3
  USE ReferenceElementModule, ONLY: &
    nDOF_X1, nDOF_X2, &
    Weights_X1, Weights_X2


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
    UseFluxCorrection_TwoMoment, &
    IsPeriodic
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMapFromBool
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
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment_MF
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

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit_MF_SinglePatch

  REAL(amrex_real) :: OffGridFlux_TwoMoment_X1_Inner(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X1_Outer(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X2_Inner(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X2_Outer(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X3_Inner(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X3_Outer(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_E       (2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X1_Inner_All(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X1_Outer_All(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X2_Inner_All(2*nCR)
  REAL(amrex_real) :: OffGridFlux_TwoMoment_X2_Outer_All(2*nCR)


CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Explicit_MF_SinglePatch( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
      lo, hi, &
      uGF, uCF, uCR, duCR, Verbose_Option ) BIND(C, name="ComputeIncrement_TwoMoment_Explicit_MF_SinglePatch")
    USE ISO_C_BINDING

    INTEGER(C_INT), INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(C_DOUBLE), INTENT(in)    :: &
      uGF (iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nDOFX*nGF)
    REAL(C_DOUBLE), INTENT(inout) :: &
      uCF (iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nDOFX*nCF)
    REAL(C_DOUBLE), INTENT(inout) :: &
      uCR (iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nDOFZ*(iZ_E0(1)-iZ_B0(1)+1)*nCR*nSpecies)
    REAL(C_DOUBLE), INTENT(out)   :: &
      duCR(iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           1:nDOFZ*(iZ_E0(1)-iZ_B0(1)+1)*nCR*nSpecies)
    LOGICAL(C_BOOL), INTENT(in) :: lo(3), hi(3)
    LOGICAL(C_BOOL), INTENT(in), OPTIONAL :: Verbose_Option

    REAL(amrex_real), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: C (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U (:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: dU(:,:,:,:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: SurfaceFlux_X1(:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: SurfaceFlux_X2(:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: SurfaceFlux_X3(:,:,:,:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    TYPE(EdgeMap) :: Edge_Map

    iE_B0=iZ_B0(1)
    iE_E0=iZ_E0(1)
    iE_B1=iZ_B1(1)
    iE_E1=iZ_E1(1)


    iX_B0(1:3)=iZ_B0(2:4)
    iX_E0(1:3)=iZ_E0(2:4)
    iX_B1(1:3)=iZ_B1(2:4)
    iX_E1(1:3)=iZ_E1(2:4)


    CALL AllocateArray_X &
           ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
             [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
             G )

    CALL AllocateArray_X &
           ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
             [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
             C )

    CALL AllocateArray_Z &
           ( [ 1       , iZ_B1(1), iZ_B1(2), iZ_B1(3), iZ_B1(4), 1       , 1        ], &
             [ nDOFZ   , iZ_E1(1), iZ_E1(2), iZ_E1(3), iZ_E1(4), nCR     , nSpecies ], &
             U )

    CALL AllocateArray_Z &
           ( [ 1       , iZ_B1(1), iZ_B1(2), iZ_B1(3), iZ_B1(4), 1       , 1        ], &
             [ nDOFZ   , iZ_E1(1), iZ_E1(2), iZ_E1(3), iZ_E1(4), nCR     , nSpecies ], &
             dU )

    CALL AllocateArray_Z &
           ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
             [ nDOFE * nDOFX_X1, iZ_E0(1), iZ_E0(2)+1, iZ_E0(3), iZ_E0(4), nCR, nSpecies ], &
             SurfaceFlux_X1 )

    CALL AllocateArray_Z &
           ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
             [ nDOFE *nDOFX_X2, iZ_E0(1), iZ_E0(2), iZ_E0(3)+1, iZ_E0(4), nCR, nSpecies ], &
             SurfaceFlux_X2 )

    CALL AllocateArray_Z &
           ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
             [ nDOFE *nDOFX_X3, iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4)+1, nCR, nSpecies ], &
             SurfaceFlux_X3 )

    CALL amrex2thornado_X( nGF, iX_B1, iX_E1, LBOUND(uGF), iX_B1, iX_E1, uGF, G )

    CALL amrex2thornado_X( nCF, iX_B1, iX_E1, LBOUND(uCF), iX_B1, iX_E1, uCF, C )

    CALL amrex2thornado_Z &
           ( nCR, nSpecies, nE, iE_B0, iE_E0, &
             iZ_B1, iZ_E1, LBOUND(uCR), iZ_B1, iZ_E1, uCR, U )

    CALL ConstructEdgeMapFromBool( lo, hi, Edge_Map )

    CALL ApplyBoundaryConditions_TwoMoment_MF &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )

    CALL ApplyBoundaryConditions_Euler_MF &
           ( iX_B0, iX_E0, iX_B1, iX_E1, C, Edge_Map )

    CALL ComputeIncrement_TwoMoment_Explicit &
          ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, G, C, U, dU, &
            SuppressBC_Option = .TRUE., &
            SurfaceFlux_X1_Option = SurfaceFlux_X1, &
            SurfaceFlux_X2_Option = SurfaceFlux_X2, &
            SurfaceFlux_X3_Option = SurfaceFlux_X3 )

    CALL thornado2amrex_Z &
           ( nCR, nSpecies, nE, iE_B0, iE_E0, &
             iZ_B1, iZ_E1, LBOUND(duCR), iZ_B0, iZ_E0, duCR, dU )

    CALL DeallocateArray_Z &
           ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
             [ nDOFE *nDOFX_X3, iZ_E0(1), iZ_E0(2), iZ_E0(3), iZ_E0(4)+1, nCR, nSpecies ], &
             SurfaceFlux_X3 )

    CALL DeallocateArray_Z &
           ( [ 1, iZ_B0(1), iZ_B0(2), iZ_B0(3), iZ_B0(4), 1, 1 ], &
             [ nDOFE *nDOFX_X2, iZ_E0(1), iZ_E0(2), iZ_E0(3)+1, iZ_E0(4), nCR, nSpecies ], &
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

  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit_MF_SinglePatch


  SUBROUTINE ComputeOffGridFlux_MF &
    ( iZ_B0, iZ_E0, FineMask, SurfaceFlux_X1, SurfaceFlux_X2 )
 
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4)
    INTEGER, INTENT(in) :: FineMask(iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(amrex_real), INTENT(in) :: &
      SurfaceFlux_X1(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)
    REAL(amrex_real), INTENT(in) :: &
      SurfaceFlux_X2(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)
 
    OffGridFlux_TwoMoment_X1_Inner = Zero
    OffGridFlux_TwoMoment_X1_Outer = Zero
    OffGridFlux_TwoMoment_X2_Inner = Zero
    OffGridFlux_TwoMoment_X2_Outer = Zero
    OffGridFlux_TwoMoment_X3_Inner = Zero
    OffGridFlux_TwoMoment_X3_Outer = Zero
 
    OffGridFlux_TwoMoment_X1_Inner_All = Zero
    OffGridFlux_TwoMoment_X1_Outer_All = Zero
    OffGridFlux_TwoMoment_X2_Inner_All = Zero
    OffGridFlux_TwoMoment_X2_Outer_All = Zero
 
    CALL FaceSum_X1 &
           ( iZ_B0, iZ_E0, iZ_B0(2)  , iZ_B0(2), FineMask, SurfaceFlux_X1, &
             OffGridFlux_TwoMoment_X1_Inner_All, &
             OffGridFlux_TwoMoment_X1_Inner )
 
    CALL FaceSum_X1 &
           ( iZ_B0, iZ_E0, iZ_E0(2)+1, iZ_E0(2), FineMask, SurfaceFlux_X1, &
             OffGridFlux_TwoMoment_X1_Outer_All, &
             OffGridFlux_TwoMoment_X1_Outer )
 
    IF( nDimsX .GT. 1 )THEN
 
      CALL FaceSum_X2 &
             ( iZ_B0, iZ_E0, iZ_B0(3)  , iZ_B0(3), FineMask, SurfaceFlux_X2, &
               OffGridFlux_TwoMoment_X2_Inner_All, &
               OffGridFlux_TwoMoment_X2_Inner )
 
      CALL FaceSum_X2 &
             ( iZ_B0, iZ_E0, iZ_E0(3)+1, iZ_E0(3), FineMask, SurfaceFlux_X2, &
               OffGridFlux_TwoMoment_X2_Outer_All, &
               OffGridFlux_TwoMoment_X2_Outer )
 
    END IF
 
 
    OffGridFlux_TwoMoment_E &
      = OffGridFlux_TwoMoment &
          - ( OffGridFlux_TwoMoment_X1_Inner_All &
                - OffGridFlux_TwoMoment_X1_Outer_All ) &
          - ( OffGridFlux_TwoMoment_X2_Inner_All &
                - OffGridFlux_TwoMoment_X2_Outer_All )

 
  END SUBROUTINE ComputeOffGridFlux_MF
 
 
  SUBROUTINE FaceSum_X1 &
    ( iZ_B0, iZ_E0, iZ2_Face, iZ2_Cell, FineMask, SurfaceFlux, F_All, F_Leaf )
 
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ2_Face, iZ2_Cell
    INTEGER, INTENT(in) :: FineMask(iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(amrex_real), INTENT(in) :: &
      SurfaceFlux(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)
    REAL(amrex_real), INTENT(out) :: F_All(1:2*nCR), F_Leaf(1:2*nCR)
 
    INTEGER          :: iCR, iS, iZ1, iZ3, iZ4, iNodeZ_X, iNodeE
    LOGICAL          :: Leaf
    REAL(amrex_real) :: w, Flux, dN, dE
 
    F_All  = Zero
    F_Leaf = Zero
 
    ASSOCIATE &
      ( dZ1 => MeshE % Width, dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )
 
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
 
      Leaf = .NOT. IsNotLeafElement( FineMask(iZ2_Cell,iZ3,iZ4,1) )
 
      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ1 = iZ_B0(1), iZ_E0(1)
      DO iNodeZ_X = 1, nDOF_X1
 
        iNodeE = MOD( iNodeZ_X - 1, nDOFE ) + 1
 
        w    = dZ1(iZ1) * dZ3(iZ3) * dZ4(iZ4) * Weights_X1(iNodeZ_X)
        Flux = SurfaceFlux(iNodeZ_X,iZ1,iZ2_Face,iZ3,iZ4,iCR,iS)
 
        dN = LeptonNumber(iS) * w * uGE(iNodeE,iZ1,iGE_Ep2) * Flux
        dE =                    w * uGE(iNodeE,iZ1,iGE_Ep3) * Flux
 
        F_All(iCR)     = F_All(iCR)     + dN
        F_All(iCR+nCR) = F_All(iCR+nCR) + dE
 
        IF( Leaf )THEN
          F_Leaf(iCR)     = F_Leaf(iCR)     + dN
          F_Leaf(iCR+nCR) = F_Leaf(iCR+nCR) + dE
        END IF
 
      END DO
      END DO
      END DO
      END DO
 
    END DO
    END DO
 
    END ASSOCIATE
 
  END SUBROUTINE FaceSum_X1
 
 
  SUBROUTINE FaceSum_X2 &
    ( iZ_B0, iZ_E0, iZ3_Face, iZ3_Cell, FineMask, SurfaceFlux, F_All, F_Leaf )
 
    INTEGER, INTENT(in) :: iZ_B0(4), iZ_E0(4), iZ3_Face, iZ3_Cell
    INTEGER, INTENT(in) :: FineMask(iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(amrex_real), INTENT(in) :: &
      SurfaceFlux(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)
    REAL(amrex_real), INTENT(out) :: F_All(1:2*nCR), F_Leaf(1:2*nCR)
 
    INTEGER          :: iCR, iS, iZ1, iZ2, iZ4, iNodeZ_X, iNodeE
    LOGICAL          :: Leaf
    REAL(amrex_real) :: w, Flux, dN, dE
 
    F_All  = Zero
    F_Leaf = Zero
 
    ASSOCIATE &
      ( dZ1 => MeshE % Width, dZ2 => MeshX(1) % Width, dZ4 => MeshX(3) % Width )
 
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
 
      Leaf = .NOT. IsNotLeafElement( FineMask(iZ2,iZ3_Cell,iZ4,1) )
 
      DO iS  = 1, nSpecies
      DO iCR = 1, nCR
      DO iZ1 = iZ_B0(1), iZ_E0(1)
      DO iNodeZ_X = 1, nDOF_X2
 
        iNodeE = MOD( iNodeZ_X - 1, nDOFE ) + 1
 
        w    = dZ1(iZ1) * dZ2(iZ2) * dZ4(iZ4) * Weights_X2(iNodeZ_X)
        Flux = SurfaceFlux(iNodeZ_X,iZ1,iZ2,iZ3_Face,iZ4,iCR,iS)
 
        dN = LeptonNumber(iS) * w * uGE(iNodeE,iZ1,iGE_Ep2) * Flux
        dE =                    w * uGE(iNodeE,iZ1,iGE_Ep3) * Flux
 
        F_All(iCR)     = F_All(iCR)     + dN
        F_All(iCR+nCR) = F_All(iCR+nCR) + dE
 
        IF( Leaf )THEN
          F_Leaf(iCR)     = F_Leaf(iCR)     + dN
          F_Leaf(iCR+nCR) = F_Leaf(iCR+nCR) + dE
        END IF
 
      END DO
      END DO
      END DO
      END DO
 
    END DO
    END DO
 
    END ASSOCIATE
 
  END SUBROUTINE FaceSum_X2


  SUBROUTINE IncrementOffGridTally_TwoMoment( iLevel, iX_B0, iX_E0 )


    INTEGER, INTENT(in) :: iLevel
    INTEGER, INTENT(in) :: iX_B0(3), iX_E0(3)

    OffGridFlux_TwoMoment_MF(:,iLevel) &
      = OffGridFlux_TwoMoment_MF(:,iLevel) + OffGridFlux_TwoMoment_E

    IF( .NOT. IsPeriodic(1) )THEN
      IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo(1) ) &
        OffGridFlux_TwoMoment_MF(:,iLevel) &
          = OffGridFlux_TwoMoment_MF(:,iLevel) + OffGridFlux_TwoMoment_X1_Inner
      IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi(1) ) &
        OffGridFlux_TwoMoment_MF(:,iLevel) &
          = OffGridFlux_TwoMoment_MF(:,iLevel) - OffGridFlux_TwoMoment_X1_Outer
    END IF

    IF( .NOT. IsPeriodic(2) )THEN
      IF( iX_B0(2) .EQ. amrex_geom(iLevel) % domain % lo(2) ) &
        OffGridFlux_TwoMoment_MF(:,iLevel) &
          = OffGridFlux_TwoMoment_MF(:,iLevel) + OffGridFlux_TwoMoment_X2_Inner
      IF( iX_E0(2) .EQ. amrex_geom(iLevel) % domain % hi(2) ) &
        OffGridFlux_TwoMoment_MF(:,iLevel) &
          = OffGridFlux_TwoMoment_MF(:,iLevel) - OffGridFlux_TwoMoment_X2_Outer
    END IF

    IF( .NOT. IsPeriodic(3) )THEN
      IF( iX_B0(3) .EQ. amrex_geom(iLevel) % domain % lo(3) ) &
        OffGridFlux_TwoMoment_MF(:,iLevel) &
          = OffGridFlux_TwoMoment_MF(:,iLevel) + OffGridFlux_TwoMoment_X3_Inner
      IF( iX_E0(3) .EQ. amrex_geom(iLevel) % domain % hi(3) ) &
        OffGridFlux_TwoMoment_MF(:,iLevel) &
          = OffGridFlux_TwoMoment_MF(:,iLevel) - OffGridFlux_TwoMoment_X3_Outer
    END IF

  END SUBROUTINE IncrementOffGridTally_TwoMoment


END MODULE  MF_TwoMoment_DiscretizationModule_Streaming_OrderV
