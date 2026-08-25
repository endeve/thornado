MODULE MF_TwoMoment_OpacityModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE FillPatchModule, ONLY: &
    FillPatch
  USE amrex_boxarray_module, ONLY: &
    amrex_boxarray
  USE amrex_distromap_module, ONLY: &
    amrex_distromap
  USE amrex_multifab_module, ONLY: &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_box_module, ONLY: &
    amrex_box
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    CreateMesh, &
    NodeCoordinate
  USE ReferenceElementModuleZ, ONLY: &
    NodeNumberTableZ
  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    swX, &
    swE, &
    zoomE, &
    nDOFX, &
    nDOFZ, &
    nDOFE, &
    iE_B0, &
    iE_E0, &
    iE_B1, &
    iE_E1, &
    iZ_B1, &
    iZ_E1, &
    iZ_B0, &
    iZ_E0, &
    nNodesE, &
    eL, &
    eR, &
    nE, &
    nX, &
    DescribeProgramHeaderX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    nSpecies
  USE FluidFieldsModule, ONLY: &
    nCF
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    InitializeSlopeLimiter_TwoMoment, &
    FinalizeSlopeLimiter_TwoMoment, &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
    InitializeTroubledCellIndicator_TwoMoment, &
    FinalizeTroubledCellIndicator_TwoMoment
  USE TwoMoment_OpacityModule, ONLY: SetOpacities, &
    nOP, iOP_D0, iOP_Chi, iOP_Sigma, &
    uOP

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Half, &
    Three
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    nMaxLevels, &
    nE
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_TwoMoment_MF
  USE MF_EdgeMapModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF

     IMPLICIT NONE
      PRIVATE
     
      TYPE(amrex_multifab), ALLOCATABLE, PUBLIC :: MF_uOP(:)
     
      PUBLIC :: InitializeOpacities_MF
      PUBLIC :: FinalizeOpacities_MF
      PUBLIC :: BuildOpacities_MF_Level
      PUBLIC :: ClearOpacities_MF_Level
      PUBLIC :: SetOpacities_MF_Level
      PUBLIC :: LoadOpacities_Box
      PUBLIC :: LoadOpacities_GlobalSlot
      PUBLIC :: FinalizeOpacities_GlobalSlot
      !PUBLIC :: nComp_uOP
     
      REAL(DP), SAVE :: D_0_Save    = Zero
      REAL(DP), SAVE :: Chi_Save    = Zero
      REAL(DP), SAVE :: Sigma_Save  = Zero
      LOGICAL,  SAVE :: Initialized = .FALSE.
     
    CONTAINS
     
     
      PURE INTEGER FUNCTION nComp_uOP()
        nComp_uOP = nDOFZ * ( iE_E1 - iE_B1 + 1 ) * nOP * nSpecies
      END FUNCTION nComp_uOP
     
     
      PURE INTEGER FUNCTION iComp_OP( iNodeZ, iE, iOP_, iS )
        INTEGER, INTENT(in) :: iNodeZ, iE, iOP_, iS
        INTEGER :: nE_loc
        nE_loc = iE_E1 - iE_B1 + 1
        iComp_OP = iNodeZ &
                 + nDOFZ * ( ( iE - iE_B1 ) &
                           + nE_loc * ( ( iOP_ - 1 ) &
                                      + nOP * ( iS - 1 ) ) )
      END FUNCTION iComp_OP
     
     
      SUBROUTINE InitializeOpacities_MF( D_0, Chi, Sigma )
     
        REAL(DP), INTENT(in) :: D_0, Chi, Sigma
     
        IF( Initialized ) RETURN
     
        ALLOCATE( MF_uOP(0:nMaxLevels-1) )
     
        D_0_Save    = D_0
        Chi_Save    = Chi
        Sigma_Save  = Sigma
        Initialized = .TRUE.
     
        IF( amrex_parallel_ioprocessor() )THEN
          WRITE(*,*)
          WRITE(*,'(A5,A)') '', 'MF_TwoMoment_OpacityModule initialised'
          WRITE(*,'(A7,A)') '', '(D0/Chi/Sigma below apply to DEFAULT ProgramName branch;'
          WRITE(*,'(A7,A)') '',  ' named cases override these in native SetOpacities)'
          WRITE(*,'(A7,A8,ES10.4E2)') '',    'D0 = ', D_0
          WRITE(*,'(A7,A8,ES10.4E2)') '',   'Chi = ', Chi
          WRITE(*,'(A7,A8,ES10.4E2)') '', 'Sigma = ', Sigma
          WRITE(*,'(A7,A,I0)') '', 'nMaxLevels = ', nMaxLevels
          WRITE(*,'(A7,A,I0)') '', 'nComp_uOP  = ', nComp_uOP()
        END IF
     
      END SUBROUTINE InitializeOpacities_MF
     
     
      SUBROUTINE FinalizeOpacities_MF
        INTEGER :: iLevel
        IF( .NOT. Initialized ) RETURN
        DO iLevel = 0, SIZE( MF_uOP ) - 1
          IF( MF_uOP(iLevel) % owner ) &
            CALL amrex_multifab_destroy( MF_uOP(iLevel) )
        END DO
        DEALLOCATE( MF_uOP )
        IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )
        Initialized = .FALSE.
      END SUBROUTINE FinalizeOpacities_MF
     
     
      SUBROUTINE ClearOpacities_MF_Level( iLevel )
        INTEGER, INTENT(in) :: iLevel
        IF( .NOT. Initialized ) RETURN
        IF( MF_uOP(iLevel) % owner ) &
          CALL amrex_multifab_destroy( MF_uOP(iLevel) )
      END SUBROUTINE ClearOpacities_MF_Level
     
     
      SUBROUTINE BuildOpacities_MF_Level( iLevel, BA, DM )
     
        INTEGER,               INTENT(in) :: iLevel
        TYPE(amrex_boxarray),  INTENT(in) :: BA
        TYPE(amrex_distromap), INTENT(in) :: DM
     
        CALL amrex_multifab_build( MF_uOP(iLevel), BA, DM, nComp_uOP(), swX )
        CALL MF_uOP(iLevel) % SetVal( Zero )
     
        CALL SetOpacities_MF_Level( iLevel )
     
      END SUBROUTINE BuildOpacities_MF_Level
     
     
      SUBROUTINE SetOpacities_MF_Level( iLevel )
     

        INTEGER, INTENT(in) :: iLevel
     
        TYPE(amrex_mfiter) :: MFI
        TYPE(amrex_box)    :: BX
        REAL(DP), CONTIGUOUS, POINTER :: pOP(:,:,:,:)
     
        INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
        INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
        INTEGER :: iNodeZ, iE, iOP_, iS
        INTEGER :: iX1, iX2, iX3, iC
     
        IF( .NOT. Initialized ) RETURN
     
        CALL CreateMesh_MF( iLevel, MeshX )
     
        CALL amrex_mfiter_build( MFI, MF_uOP(iLevel), Tiling = .FALSE. )
     
        DO WHILE( MFI % next() )
     
          BX  =  MFI % TileBox()
          pOP => MF_uOP(iLevel) % DataPtr( MFI )
     
          iX_B0 = BX % lo
          iX_E0 = BX % hi
          iX_B1 = BX % lo - swX
          iX_E1 = BX % hi + swX
     
          iZ_B0(1)   = iE_B0
          iZ_E0(1)   = iE_E0
          iZ_B1(1)   = iE_B1
          iZ_E1(1)   = iE_E1
          iZ_B0(2:4) = iX_B0
          iZ_E0(2:4) = iX_E0
          iZ_B1(2:4) = iX_B1
          iZ_E1(2:4) = iX_E1
     
          IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )
          ALLOCATE( uOP( 1:nDOFZ, &
                         iZ_B1(1):iZ_E1(1), &
                         iZ_B1(2):iZ_E1(2), &
                         iZ_B1(3):iZ_E1(3), &
                         iZ_B1(4):iZ_E1(4), &
                         1:nOP, 1:nSpecies ) )
          uOP = Zero
     
          CALL SetOpacities &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
                   D_0_Save, Chi_Save, Sigma_Save, &
                   Verbose_Option = .TRUE. )
     
          DO iS     = 1, nSpecies
          DO iOP_   = 1, nOP
          DO iE     = iE_B1, iE_E1
          DO iNodeZ = 1, nDOFZ
            iC = iComp_OP( iNodeZ, iE, iOP_, iS )
            DO iX3 = iZ_B1(4), iZ_E1(4)
            DO iX2 = iZ_B1(3), iZ_E1(3)
            DO iX1 = iZ_B1(2), iZ_E1(2)
              pOP(iX1, iX2, iX3, iC) &
                = uOP(iNodeZ, iE, iX1, iX2, iX3, iOP_, iS)
            END DO
            END DO
            END DO
          END DO
          END DO
          END DO
          END DO
     
        END DO
     
        CALL amrex_mfiter_destroy( MFI )
     
        IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )
     
        CALL DestroyMesh_MF( MeshX )
     
      END SUBROUTINE SetOpacities_MF_Level
     
     
      SUBROUTINE LoadOpacities_Box( iLevel, MFI, uOP_Box )
     
        INTEGER,            INTENT(in)     :: iLevel
        TYPE(amrex_mfiter), INTENT(in)     :: MFI
        REAL(DP), ALLOCATABLE, INTENT(out) :: uOP_Box(:,:,:,:,:,:,:)
     
        REAL(DP), CONTIGUOUS, POINTER :: pOP(:,:,:,:)
        INTEGER :: lo1, lo2, lo3, hi1, hi2, hi3
        INTEGER :: iNodeZ, iE, iOP_, iS
        INTEGER :: iX1, iX2, iX3, iC
     
        pOP => MF_uOP(iLevel) % DataPtr( MFI )
     
        lo1 = LBOUND( pOP, 1 );  hi1 = UBOUND( pOP, 1 )
        lo2 = LBOUND( pOP, 2 );  hi2 = UBOUND( pOP, 2 )
        lo3 = LBOUND( pOP, 3 );  hi3 = UBOUND( pOP, 3 )
     
        ALLOCATE( uOP_Box(1:nDOFZ, iE_B1:iE_E1, &
                          lo1:hi1, lo2:hi2, lo3:hi3, &
                          1:nOP,   1:nSpecies) )
     
        DO iS    = 1, nSpecies
        DO iOP_  = 1, nOP
        DO iE    = iE_B1, iE_E1
        DO iNodeZ = 1, nDOFZ
          iC = iComp_OP( iNodeZ, iE, iOP_, iS )
          DO iX3 = lo3, hi3
          DO iX2 = lo2, hi2
          DO iX1 = lo1, hi1
            uOP_Box(iNodeZ, iE, iX1, iX2, iX3, iOP_, iS) &
              = pOP(iX1, iX2, iX3, iC)
          END DO
          END DO
          END DO
        END DO
        END DO
        END DO
        END DO
     
      END SUBROUTINE LoadOpacities_Box
     
     
      SUBROUTINE LoadOpacities_GlobalSlot( iLevel, MFI, iZ_B1, iZ_E1 )
     
        INTEGER,            INTENT(in) :: iLevel
        TYPE(amrex_mfiter), INTENT(in) :: MFI
        INTEGER,            INTENT(in) :: iZ_B1(4), iZ_E1(4)
     
        REAL(DP), CONTIGUOUS, POINTER :: pOP(:,:,:,:)
        INTEGER :: iNodeZ, iE, iOP_, iS, iX1, iX2, iX3, iC
     
        pOP => MF_uOP(iLevel) % DataPtr( MFI )
     
        IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )
     
        ALLOCATE( uOP( 1:nDOFZ, &
                       iZ_B1(1):iZ_E1(1), &
                       iZ_B1(2):iZ_E1(2), &
                       iZ_B1(3):iZ_E1(3), &
                       iZ_B1(4):iZ_E1(4), &
                       1:nOP, 1:nSpecies ) )
     
        uOP = Zero
     
        DO iS     = 1, nSpecies
        DO iOP_   = 1, nOP
        DO iE     = iZ_B1(1), iZ_E1(1)
        DO iNodeZ = 1, nDOFZ
          iC = iComp_OP( iNodeZ, iE, iOP_, iS )
          DO iX3 = iZ_B1(4), iZ_E1(4)
          DO iX2 = iZ_B1(3), iZ_E1(3)
          DO iX1 = iZ_B1(2), iZ_E1(2)
            uOP(iNodeZ, iE, iX1, iX2, iX3, iOP_, iS) &
              = pOP(iX1, iX2, iX3, iC)
          END DO
          END DO
          END DO
        END DO
        END DO
        END DO
        END DO
     
      END SUBROUTINE LoadOpacities_GlobalSlot
     
     
      SUBROUTINE FinalizeOpacities_GlobalSlot
        IF( ALLOCATED( uOP ) ) DEALLOCATE( uOP )
      END SUBROUTINE FinalizeOpacities_GlobalSlot
     
     
    END MODULE MF_TwoMoment_OpacityModule