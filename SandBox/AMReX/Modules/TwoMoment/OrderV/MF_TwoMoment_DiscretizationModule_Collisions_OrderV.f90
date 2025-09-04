MODULE  MF_TwoMoment_DiscretizationModule_Collisions_OrderV

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
  USE TwoMoment_DiscretizationModule_Collisions, ONLY: &
    ComputeIncrement_TwoMoment_Implicit
  USE MeshModule, ONLY: &
    MeshX

  ! --- Local Modules ---
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
    ApplyBoundaryConditions_TwoMoment_MF
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Euler_MF
  USE MF_EdgeMapModule, ONLY: &
    ConstructEdgeMap, &
    EdgeMap

USE MF_KindModule, ONLY: &
    Zero
USE MF_UtilitiesModule
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_TwoMoment_Implicit_MF


CONTAINS


  SUBROUTINE ComputeIncrement_TwoMoment_Implicit_MF &
    ( GEOM, MF_uGF, MF_uCF, MF_uCR, MF_duCR, MF_duCF, dt, Verbose_Option )

    TYPE(amrex_geometry), INTENT(in)    :: GEOM   (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF (0:nLevels-1)
    REAL(amrex_real),     INTENT(in)    :: dt
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_duCF(0:nLevels-1)
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCR (:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: duCR(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: duCF(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: GX  (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U_F (:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: dU_F(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U_R (:,:,:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: dU_R(:,:,:,:,:,:,:)

    INTEGER :: iLevel, iZ1, iZ2, iZ3, iZ4, iNodeZ
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4), i

    LOGICAL :: Verbose

    TYPE(EdgeMap) :: Edge_Map

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---
      CALL MF_uCR(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_duCR(iLevel) % setval( 0.0_amrex_real )
      CALL MF_duCF(iLevel) % setval( 0.0_amrex_real )

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF (iLevel) % DataPtr( MFI )
        uCF  => MF_uCF (iLevel) % DataPtr( MFI )
        uCR  => MF_uCR (iLevel) % DataPtr( MFI )
        duCR => MF_duCR(iLevel) % DataPtr( MFI )
        duCF => MF_duCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        i=1

        DO WHILE (i<=4)

          IF (i==1) THEN

            iZ_B0(i)=iE_B0
            iZ_E0(i)=iE_E0
            iZ_B1(i)=iE_B1
            iZ_E1(i)=iE_E1

          ELSE

            iZ_B0(i)=iX_B0(i-1)
            iZ_E0(i)=iX_E0(i-1)
            iZ_B1(i)=iX_B1(i-1)
            iZ_E1(i)=iX_E1(i-1)

          END IF
          i = i + 1
        END DO

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 GX )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U_F )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 dU_F )

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
                 U_R )

!PRINT *, SHAPE(U_R)

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
                 dU_R )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, GX )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uCF, U_F )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B1, iZ_E1, uCR, U_R )

        !PRINT *, 'U_R', U_R(:,:,:,:,:,:,:)
        !PRINT *, '****************************************************'
        !PRINT *, 'F', U_F (:,:,:,:,1)

        CALL ConstructEdgeMap( iLevel, BX, Edge_Map )

        CALL ApplyBoundaryConditions_Euler_MF &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U_F, Edge_Map )

        CALL ApplyBoundaryConditions_TwoMoment_MF &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, U_R, Edge_Map )

        !PRINT *, 'U_R', U_R(:,:,:,:,:,1,:)
        !PRINT *, '****************************************************'
        !PRINT *, 'F', U_F (:,:,:,:,1)
        CALL ComputeIncrement_TwoMoment_Implicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, uGE, GX, U_F, dU_F, U_R, dU_R )
        PRINT *, 'MADE IT PAST IMPLICIT'
        CALL thornado2amrex_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, duCR, dU_R )

        CALL thornado2amrex_X &
               ( nCF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, duCF, dU_F )

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
                 dU_R )


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
                 U_R )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 dU_F )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U_F )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 GX )

      END DO

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO

  END SUBROUTINE ComputeIncrement_TwoMoment_Implicit_MF


END MODULE  MF_TwoMoment_DiscretizationModule_Collisions_OrderV
