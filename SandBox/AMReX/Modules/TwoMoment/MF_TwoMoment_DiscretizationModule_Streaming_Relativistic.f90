MODULE  MF_TwoMoment_DiscretizationModule_Streaming_Relativistic

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter,   &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX, nDOFZ, swE, nDOFE, iE_B0, iE_E0, iE_B1, iE_E1
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE GeometryFieldsModuleE,     ONLY: &
    nGE, uGE
  USE RadiationFieldsModule,            ONLY: &
    nCR
  USE FluidFieldsModule,            ONLY: &
    nCF
  USE TwoMoment_DiscretizationModule_Streaming_Relativistic, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
  USE MeshModule, ONLY: &
    MeshX

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
    nSpecies, &
    UseTiling, &
    nE
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap, &
    ApplyBoundaryConditions_TwoMoment_MF
  USE FillPatchModule, ONLY: &
    FillPatch
  USE AverageDownModule, ONLY: &
    AverageDown
  USE MF_FieldsModule_TwoMoment, ONLY: &
    MF_Permute
  USE MF_UtilitiesModule, ONLY: &
    amrex2amrex_permute_Z, &
    amrex_permute2amrex_Z, &
    MF_amrex2amrex_permute_Z_Level, &
    MF_amrex_permute2amrex_Z_Level
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Explicit_MF


CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Explicit_MF( Time, GEOM, MF_uGF, MF_uCF, MF_uCR, MF_duCR, Verbose_Option )


    REAL(amrex_real),     INTENT(in)    :: Time(0:)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCR(0:nLevels-1)
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

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

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), i, j


    LOGICAL :: Verbose

    TYPE(EdgeMap) :: Edge_Map

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option


    DO i = 0, nLevels-1

      IF (i .NE. 0) THEN
       
        DO j = i - 1, i

          !CALL MF_amrex2amrex_permute_Z_Level(j,nCR,MF_uGF(j),MF_uCR(j),MF_Permute(j))

        END DO

        CALL FillPatch( i, Time(i), MF_uGF, MF_Permute )
      ELSE

        CALL FillPatch( i, Time(i), MF_uGF, MF_uCR )

      END IF

      IF (i .NE. 0) THEN

        DO j = i - 1, i

          !CALL MF_amrex_permute2amrex_Z_Level(j,nCR,MF_uGF(j),MF_uCR(j),MF_Permute(j))

        END DO

      END IF
    END DO


    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uCR(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_duCR(iLevel) % setval( 0.0_amrex_real )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCF  => MF_uCF (iLevel) % DataPtr( MFI )
        uCR  => MF_uCR (iLevel) % DataPtr( MFI )
        duCR => MF_duCR(iLevel) % DataPtr( MFI )

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


        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, C )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B1, iZ_E1, uCR, U )

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

        CALL ApplyBoundaryConditions_TwoMoment_MF &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U, Edge_Map )



       CALL ComputeIncrement_TwoMoment_Explicit &
              ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, G, C, U, dU, &
                Verbose_Option = Verbose, &
                SuppressBC_Option = .TRUE.  )
        CALL thornado2amrex_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, duCR, dU )

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

      CALL DestroyMesh_MF( MeshX )

    END DO

    DO i = 0, nLevels-1

     ! CALL MF_amrex2amrex_permute_Z_Level(i,nCR,MF_uGF(i),MF_uCR(i),MF_Permute(i))

    END DO


   ! CALL AverageDown( MF_uGF, MF_Permute )

    DO i = 0, nLevels-1

    !  CALL MF_amrex_permute2amrex_Z_Level(i,nCR,MF_uGF(i),MF_uCR(i),MF_Permute(i))

    END DO


  END SUBROUTINE ComputeIncrement_TwoMoment_Explicit_MF


END MODULE  MF_TwoMoment_DiscretizationModule_Streaming_Relativistic
