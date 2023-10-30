MODULE MF_XCFC_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFZ, &
    nNodesX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodeNumberTableX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Up
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Psi, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33, &
    nGF
  USE GeometryBoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting, &
    ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    nSp => nSpecies
  USE XCFC_UtilitiesModule, ONLY: &
    iMF_Psi, &
    iMF_Alpha, &
    iMF_Beta_1, &
    iMF_Beta_2, &
    iMF_Beta_3, &
    iMF_K_dd_11, &
    iMF_K_dd_12, &
    iMF_K_dd_13, &
    iMF_K_dd_22, &
    iMF_K_dd_23, &
    iMF_K_dd_33, &
    nMF, &
    nGS, &
    MultiplyWithPsi6, &
    ComputeGravitationalMass, &
    UpdateConformalFactorAndMetric_XCFC

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    SqrtTiny, &
    One, &
    Two, &
    Pi
  USE MF_UtilitiesModule, ONLY: &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z, &
    amrex2thornado_X, &
    thornado2amrex_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE InputParsingModule, ONLY: &
    nLevels, &
    swX, &
    UseTiling

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MultiplyWithPsi6_MF
  PUBLIC :: ApplyBoundaryConditions_Geometry_XCFC_MF
  PUBLIC :: ComputeGravitationalMass_MF
  PUBLIC :: UpdateConformalFactorAndMetric_XCFC_MF
  PUBLIC :: UpdateGeometry_MF

  PUBLIC :: PopulateMF_uMF
  PUBLIC :: UpdateGeometryFields_MF

  INTEGER, PUBLIC :: swXX(3)

  INTERFACE MultiplyWithPsi6_MF
    MODULE PROCEDURE MultiplyWithPsi6_MF_X
    MODULE PROCEDURE MultiplyWithPsi6_MF_Z
  END INTERFACE MultiplyWithPsi6_MF

CONTAINS


  SUBROUTINE MultiplyWithPsi6_MF_X &
    ( MF_uGF, MF_uCF, Power )

    INTEGER             , INTENT(in)    :: Power
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER  :: iLevel
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G, U, uGF, uCF, uFM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

        CALL MultiplyWithPsi6 &
               ( iX_B1, iX_E1, G, U, Power, Mask_Option = uFM )

        CALL thornado2amrex_X &
               ( nCF, iX_B1, iX_E1, LBOUND( uCF ), iX_B1, iX_E1, uCF, U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE MultiplyWithPsi6_MF_X


  SUBROUTINE MultiplyWithPsi6_MF_Z &
    ( iE_B0, iE_E0, iE_B1, iE_E1, MF_uGF, MF_uCR, Power )

    INTEGER             , INTENT(in)    :: iE_B0, iE_E0, iE_B1, iE_E1, Power
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)

    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI
    TYPE(amrex_imultifab) :: iMF_FineMask

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: uFM(:,:,:,:)

    INTEGER  :: iLevel
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), nE
    INTEGER  :: iZ_B1(4), iZ_E1(4)

    nE = iE_E0 - iE_B0 + 1

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G, U, uGF, uCR, uFM, &
      !$OMP          iX_B0, iX_E0, iX_B1, iX_E1, iZ_B1, iZ_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uFM => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        iZ_B1(1)   = iE_B1
        iZ_B1(2:4) = iX_B1
        iZ_E1(1)   = iE_E1
        iZ_E1(2:4) = iX_E1

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_Z &
               ( [ 1    , iZ_B1(1), iZ_B1(2), iZ_B1(3), iZ_B1(4), 1  , 1   ], &
                 [ nDOFZ, iZ_E1(1), iZ_E1(2), iZ_E1(3), iZ_E1(4), nCR, nSp ], &
                 U )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_Z &
               ( nCR, nSp, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, LBOUND( uCR ), iZ_B1, iZ_E1, uCR, U )

        CALL MultiplyWithPsi6 &
               ( iE_B1, iE_E1, iX_B1, iX_E1, G, U, Power, Mask_Option = uFM )

        CALL thornado2amrex_Z &
               ( nCR, nSp, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, LBOUND( uCR ), iZ_B1, iZ_E1, uCR, U )

        CALL DeallocateArray_Z &
               ( [ 1    , iZ_B1(1), iZ_B1(2), iZ_B1(3), iZ_B1(4), 1  , 1   ], &
                 [ nDOFZ, iZ_E1(1), iZ_E1(2), iZ_E1(3), iZ_E1(4), nCR, nSp ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE MultiplyWithPsi6_MF_Z


  SUBROUTINE UpdateConformalFactorAndMetric_XCFC_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1) ! Metric Fields
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), ALLOCATABLE :: M(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uMF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER  :: iLevel
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, MF, GF, uMF, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swXX
        iX_E1 = iX_E0 + swXX

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nMF ], &
                 M )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL amrex2thornado_X &
               ( nMF, iX_B0, iX_E0, LBOUND( uMF ), iX_B0, iX_E0, uMF, M )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL UpdateConformalFactorAndMetric_XCFC &
               ( iX_B0, iX_E0, iX_B1, iX_E1, M, G )

        CALL thornado2amrex_X &
               ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nMF ], &
                 M )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      CALL DestroyMesh_MF( MeshX )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE UpdateConformalFactorAndMetric_XCFC_MF


  SUBROUTINE UpdateGeometry_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uMF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uMF, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swXX
        iX_E1 = iX_E0 + swXX

        DO iX3 = iX_B1(3), iX_E1(3)
        DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)
        DO iNX = 1       , nDOFX

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Alpha-1)+iNX)

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_1-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_2-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_3-1)+iNX)

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_11-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_11-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_12-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_12-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_13-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_13-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_22-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_22-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_23-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_23-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_33-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_33-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE UpdateGeometry_MF


  SUBROUTINE PopulateMF_uMF( MF_uGF, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

    TYPE(amrex_imultifab) :: iMF_FineMask
    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uMF     (:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uGF, uMF, FineMask, iX_B0, iX_E0 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uMF      => MF_uMF(iLevel) % DataPtr( MFI )
        FineMask => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Psi-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Alpha-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_1-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_2-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_3-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_11-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_11-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_12-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_12-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_13-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_13-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_22-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_22-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_23-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_23-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_33-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_33-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE PopulateMF_uMF


  SUBROUTINE ComputeGravitationalMass_MF( MF_uGS, GravitationalMass )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:nLevels-1)
    REAL(DP)            , INTENT(out) :: GravitationalMass

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), ALLOCATABLE :: GS(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)

    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLevel
    REAL(DP) :: GravitationalMass_OMP

    ! --- Assuming 1D spherical symmetry ---

    GravitationalMass = Zero

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

      GravitationalMass_OMP = Zero

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, GS, uGS, iX_B0, iX_E0, iX_B1, iX_E1 ) &
      !$OMP REDUCTION( +:GravitationalMass_OMP )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGS(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGS ], &
                 GS )

        CALL amrex2thornado_X &
               ( nGS, iX_B0, iX_E0, LBOUND( uGS ), iX_B0, iX_E0, uGS, GS )

        CALL ComputeGravitationalMass &
               ( iX_B0, iX_E0, iX_B1, iX_E1, GS, GravitationalMass_OMP )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGS ], &
                 GS )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

      GravitationalMass = GravitationalMass_OMP

      CALL DestroyMesh_MF( MeshX )

    END DO ! iLevel = 0, nLevels-1

    CALL amrex_parallel_reduce_sum( GravitationalMass )

  END SUBROUTINE ComputeGravitationalMass_MF


  SUBROUTINE ApplyBoundaryConditions_Geometry_XCFC_MF( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    CALL ApplyBoundaryConditions_X1_Inner( MF_uGF )
    CALL ApplyBoundaryConditions_X1_Outer( MF_uGF )

  END SUBROUTINE ApplyBoundaryConditions_Geometry_XCFC_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE ApplyBoundaryConditions_X1_Inner( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    ! --- Inner Boundary: Reflecting ---

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

          CALL AllocateArray_X &
                 ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                   [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                   G )

          CALL amrex2thornado_X &
                 ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

          CALL ApplyBoundaryConditions_Geometry_X1_Inner_Reflecting &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G )

          CALL thornado2amrex_X &
                 ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

          CALL DeallocateArray_X &
                 ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                   [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                   G )

        END IF

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyBoundaryConditions_X1_Inner


  SUBROUTINE ApplyBoundaryConditions_X1_Outer( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLevel
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    ! --- Outer Boundary: Extrapolate fields to face ---

    DO iLevel = 0, nLevels-1

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, G, uGF, iX_B0, iX_E0, iX_B1, iX_E1 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swX
        iX_E1 = iX_E0 + swX

        IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) )THEN

          CALL AllocateArray_X &
                 ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                   [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                   G )

          CALL amrex2thornado_X &
                 ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

          CALL ApplyBoundaryConditions_Geometry_X1_Outer_ExtrapolateToFace &
                 ( iX_B0, iX_E0, iX_B1, iX_E1, G )

          CALL thornado2amrex_X &
                 ( nGF, iX_B1, iX_E1, LBOUND( uGF ), iX_B1, iX_E1, uGF, G )

          CALL DeallocateArray_X &
                 ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                   [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                   G )

        END IF

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE ApplyBoundaryConditions_X1_Outer


  SUBROUTINE UpdateGeometryFields_MF( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, iNX1, iNX2
    INTEGER  :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP) :: X1, X2, Psi, h1, h2, h3

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      CALL CreateMesh_MF( iLevel, MeshX )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = iX_B0 - swXX
        iX_E1 = iX_E0 + swXX

        DO iX3 = iX_B1(3), iX_E1(3)
        DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)
        DO iNX = 1       , nDOFX

          iNX1 = NodeNumberTableX(1,iNX)
          iNX2 = NodeNumberTableX(2,iNX)

          X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

          Psi = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)
          h1  = Psi**2
          h2  = Psi**2 * X1
          h3  = Psi**2 * X1 * SIN( X2 )

!          uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX) = Psi
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) = h1
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) = h2
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX) = h3

          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) = MAX( h1**2, SqrtTiny )
          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) = MAX( h2**2, SqrtTiny )
          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) = MAX( h3**2, SqrtTiny )

          uGF(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX) = h1 * h2 * h3

        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL DestroyMesh_MF( MeshX )

      CALL amrex_mfiter_destroy( MFI )

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE UpdateGeometryFields_MF


END MODULE MF_XCFC_UtilitiesModule
