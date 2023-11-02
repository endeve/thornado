MODULE MF_TwoMoment_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    swE, &
    nDOFX, &
    nDOFZ, &
    nDOFE, &
    iE_B0, &
    iE_E0, &
    iE_B1, &
    iE_E1
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    iCR_N, &
    iCR_G1, &
    iCR_G2, &
    iCR_G3, &
    nPR, &
    iPR_D, &
    iPR_I1,  &
    iPR_I2,  &
    iPR_I3,  &
    iGR_N,   &
    iGR_D,   &
    iGR_I1,  &
    iGR_I2,  &
    iGR_I3,  &
    iGR_J,   &
    iGR_H1,  &
    iGR_H2,  &
    iGR_H3,  &
    iGR_RMS, &
    iGR_F,   &
    iGR_K,   &
    iGR_Q,   &
    nGR,     &
    LeptonNumber
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nAF, &
    nDF
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3
  USE GeometryFieldsModuleE,     ONLY: &
    uGE, &
    nGE, &
    iGE_Ep2, &
    iGE_Ep3
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment_Vector_Richardson
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic, &
    ComputeFromConserved_Euler_Relativistic
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE UnitsModule, ONLY: &
    UnitsActive, &
    AtomicMassUnit, &
    SpeedOfLight, &
    PlanckConstant


  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two, &
    FourPi
  USE InputParsingModule, ONLY: &
    nLevels, &
    nSpecies, &
    nE, &
    UseTiling
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X, &
    amrex2thornado_Z, &
    thornado2amrex_Z, &
    amrex2thornado_Integrated, &
    thornado2amrex_Integrated, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z, &
    AllocateArray_Integrated, &
    DeallocateArray_Integrated

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeTimeStep_TwoMoment_Fancy_MF
  PUBLIC :: ComputeTimeStep_TwoMoment_MF
  PUBLIC :: ComputeFromConserved_TwoMoment_MF
  PUBLIC :: ComputeGray_TwoMoment_MF


CONTAINS


  SUBROUTINE ComputeTimeStep_TwoMoment_Fancy_MF &
    ( MF_uGF, nX, nNodes, xR, xL, CFL, TimeStepMin )

    TYPE(amrex_multifab),  INTENT(in)  :: MF_uGF(0:nLevels-1)
    INTEGER             ,  INTENT(in)  :: nX(:), nNodes
    REAL(DP)            ,  INTENT(in)  :: xR(:), xL(:)
    REAL(DP)            ,  INTENT(in)  :: CFL
    REAL(DP)            ,  INTENT(out) :: TimeStepMin(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)

    REAL(DP) :: TimeStep(0:nLevels-1)
    INTEGER  :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER  :: iLo_MF(4)

    TimeStepMin = HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0(1:3) = BX % lo(1:3)
        iX_E0(1:3) = BX % hi(1:3)
        iX_B1(1:3) = BX % lo(1:3) - swX(1:3)
        iX_E1(1:3) = BX % hi(1:3) + swX(1:3)

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL CalculateTimeStep &
               ( iX_B1, iX_E1, iX_B0, iX_E0, CFL, G, TimeStep( iLevel) )

        TimeStepMin( iLevel ) = MIN( TimeStepMin( iLevel ), TimeStep( iLevel ) )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStepMin, SIZE(TimeStepMin) )

  END SUBROUTINE ComputeTimeStep_TwoMoment_Fancy_MF


  SUBROUTINE ComputeTimeStep_TwoMoment_MF &
    ( nX, xR, xL, nNodes, CFL, TimeStepMin )

    INTEGER , INTENT(in)  :: nX(:), nNodes
    REAL(DP), INTENT(in)  :: xR(:), xL(:), CFL
    REAL(DP), INTENT(out) :: TimeStepMin(0:nLevels-1)

    INTEGER  :: iLevel
    REAL(DP) :: TimeStep(0:nLevels-1)

    TimeStepMin = HUGE( One )

    DO iLevel = 0, nLevels-1

      TimeStep( iLevel ) = CFL * MINVAL( (xR-xL) / DBLE(nX) ) &
                           / ( Two * DBLE(nNodes-1) + One )

      TimeStepMin( iLevel ) = MIN( TimeStepMin( iLevel ), TimeStep( iLevel ) )

    END DO

  END SUBROUTINE ComputeTimeStep_TwoMoment_MF


  SUBROUTINE ComputeFromConserved_TwoMoment_MF &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uPR )

    TYPE(amrex_multifab), INTENT(in)    :: &
      MF_uGF(0:nLevels-1), MF_uCR(0:nLevels-1), MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uPR(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPR(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CF(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PF(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CR(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR(:,:,:,:,:,:,:)
    REAL(DP) :: PR_D(1), PR_I1(1), PR_I2(1), PR_I3(1)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iLo_MF(4), iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    INTEGER :: iNX, iX1, iX2, iX3
    INTEGER :: iNZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER :: ErrorExists
    INTEGER :: nIter(1)

    INTEGER, ALLOCATABLE :: ITERATION_Euler(:,:,:,:)
    INTEGER, ALLOCATABLE :: iErr_Euler     (:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uPR => MF_uPR(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iZ_B0(1) = iE_B0
        iZ_E0(1) = iE_E0
        iZ_B1(1) = iE_B1
        iZ_E1(1) = iE_E1

        iZ_B0(2:4) = iX_B0
        iZ_E0(2:4) = iX_E0
        iZ_B1(2:4) = iX_B1
        iZ_E1(2:4) = iX_E1

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 CF )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 PF )

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
                 CR )

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
                   nPR     , &
                   nSpecies ], &
                 PR )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X &
               ( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, CF )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uCR, CR )

        CALL amrex2thornado_Z &
               ( nPR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uPR, PR )

        ALLOCATE( ITERATION_Euler(1:nDOFX,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3)) )
        ALLOCATE( iErr_Euler     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                          iX_B0(2):iX_E0(2), &
                                          iX_B0(3):iX_E0(3)) )

        ErrorExists = 0

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          ITERATION_Euler(iNX,iX1,iX2,iX3) = 0
          iErr_Euler     (iNX,iX1,iX2,iX3) = 0

          CALL ComputePrimitive_Euler_Relativistic &
               ( CF(iNX,iX1,iX2,iX3,iCF_D ), &
                 CF(iNX,iX1,iX2,iX3,iCF_S1), &
                 CF(iNX,iX1,iX2,iX3,iCF_S2), &
                 CF(iNX,iX1,iX2,iX3,iCF_S3), &
                 CF(iNX,iX1,iX2,iX3,iCF_E ), &
                 CF(iNX,iX1,iX2,iX3,iCF_Ne), &
                 PF(iNX,iX1,iX2,iX3,iPF_D ), &
                 PF(iNX,iX1,iX2,iX3,iPF_V1), &
                 PF(iNX,iX1,iX2,iX3,iPF_V2), &
                 PF(iNX,iX1,iX2,iX3,iPF_V3), &
                 PF(iNX,iX1,iX2,iX3,iPF_E ), &
                 PF(iNX,iX1,iX2,iX3,iPF_Ne), &
                 G (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G (iNX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 ITERATION_Option = ITERATION_Euler(iNX,iX1,iX2,iX3), &
                 iErr_Option      = iErr_Euler     (iNX,iX1,iX2,iX3) )

          ErrorExists = ErrorExists + iErr_Euler(iNX,iX1,iX2,iX3)

        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .NE. 0 )THEN

          CALL CreateMesh_MF( iLevel, MeshX )

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            IF( iErr_Euler(iNX,iX1,iX2,iX3) .NE. 0 )THEN

            CALL DescribeError_Euler &
              ( iErr_Euler(iNX,iX1,iX2,iX3), &
                Int_Option = [ ITERATION_Euler(iNX,iX1,iX2,iX3), 99999999, &
                               iX_B0(1), iX_B0(2), iX_B0(3), &
                               iX_E0(1), iX_E0(2), iX_E0(3), &
                               iNX, iX1, iX2, iX3 ], &
                Real_Option = [ MeshX(1) % Center(iX1), &
                                MeshX(2) % Center(iX2), &
                                MeshX(3) % Center(iX3), &
                                MeshX(1) % Width (iX1), &
                                MeshX(2) % Width (iX2), &
                                MeshX(3) % Width (iX3), &
                                CF(iNX,iX1,iX2,iX3,iCF_D ), &
                                CF(iNX,iX1,iX2,iX3,iCF_S1), &
                                CF(iNX,iX1,iX2,iX3,iCF_S2), &
                                CF(iNX,iX1,iX2,iX3,iCF_S3), &
                                CF(iNX,iX1,iX2,iX3,iCF_E ), &
                                CF(iNX,iX1,iX2,iX3,iCF_Ne), &
                                G (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                                G (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                                G (iNX,iX1,iX2,iX3,iGF_Gm_dd_33) ], &
                Char_Option = [ 'NA' ], &
                Message_Option &
                  = 'Calling from ComputeFromConserved_TwoMoment_MF (Fluid)' )


            END IF

          END DO
          END DO
          END DO
          END DO

          CALL DestroyMesh_MF( MeshX )

        END IF

        DO iS  = 1       , nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
        DO iZ1 = iZ_B0(1), iZ_E0(1)
        DO iNZ = 1       , nDOFZ

          iNX = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

          CALL ComputePrimitive_TwoMoment_Vector_Richardson &
               ( [CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)], &
                 [CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)], &
                 [CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)], &
                 [CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)], &
                 PR_D, PR_I1, PR_I2, PR_I3, &
!                 [PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS)], &
!                 [PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS)], &
!                 [PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS)], &
!                 [PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS)], &
                 [PF(iNX    ,iZ2,iZ3,iZ4,iPF_V1)], &
                 [PF(iNX    ,iZ2,iZ3,iZ4,iPF_V2)], &
                 [PF(iNX    ,iZ2,iZ3,iZ4,iPF_V3)], &
                 [G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11)], &
                 [G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22)], &
                 [G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33)], &
                 !Zero, Zero, Zero, &
                 [G (iNX    ,iZ2,iZ3,iZ4,iGF_Alpha )], &
                 [G (iNX    ,iZ2,iZ3,iZ4,iGF_Beta_1)], &
                 [G (iNX    ,iZ2,iZ3,iZ4,iGF_Beta_2)], &
                 [G (iNX    ,iZ2,iZ3,iZ4,iGF_Beta_3)], &
                 [1], &
                 nIterations_Option = nIter )

          PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = PR_D (1)
          PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = PR_I1(1)
          PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = PR_I2(1)
          PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = PR_I3(1)

          IF( nIter(1) .GT. 990 )THEN

            WRITE(*,*) 'ERROR: ComputeFromConserved_TwoMoment_MF (Radiation)'
            WRITE(*,*) 'iZ_B0: ', iZ_B0
            WRITE(*,*) 'iZ_E0: ', iZ_E0

            WRITE(*,*) 'iNZ, iZ1, iZ2, iZ3, iZ4, iS: ', &
                        iNZ, iZ1, iZ2, iZ3, iZ4, iS

            iNX = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

            WRITE(*,'(2x,A,ES24.16E3)') &
              ' CR_N: ', CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'CR_G1: ', CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'CR_G2: ', CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'CR_G3: ', CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)
            WRITE(*,*)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'PF_V1: ', PF(iNX,iZ2,iZ3,iZ4,iPF_V1)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'PF_V2: ', PF(iNX,iZ2,iZ3,iZ4,iPF_V2)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'PF_V3: ', PF(iNX,iZ2,iZ3,iZ4,iPF_V3)
            WRITE(*,*)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'GF_g1: ', G (iNX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'GF_g2: ', G (iNX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'GF_g3: ', G (iNX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'GF_Al: ', G (iNX,iZ2,iZ3,iZ4,iGF_Alpha)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'GF_b1: ', G (iNX,iZ2,iZ3,iZ4,iGF_Beta_1)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'GF_b2: ', G (iNX,iZ2,iZ3,iZ4,iGF_Beta_2)
            WRITE(*,'(2x,A,ES24.16E3)') &
              'GF_b3: ', G (iNX,iZ2,iZ3,iZ4,iGF_Beta_3)

            CALL DescribeError_Euler( 99 )

          END IF

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO

        CALL thornado2amrex_Z &
               ( nPR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uPR, PR )

        DEALLOCATE( ITERATION_Euler )
        DEALLOCATE( iErr_Euler      )

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
                   nPR     , &
                   nSpecies ], &
                 PR )

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
                 CR )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 PF )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nCF ], &
                 CF )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE ComputeFromConserved_TwoMoment_MF


  SUBROUTINE CalculateTimeStep( iX_B1, iX_E1, iX_B0, iX_E0, CFL, G, dt)

    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iX_B0(3), iX_E0(3)
    REAL(DP),     INTENT(in)  :: CFL
    REAL(DP),     INTENT(in) ::  &
      G (1:nDOFX,iX_B1(1):iX_E1(1), iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3),1:nGF)
    REAL(DP),     INTENT(out) :: dt


    REAL(DP) ::   dX(3)
    INTEGER :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: dt_min(3), dt_s

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )


    dt_min = 100000.0_DP
    dt_s = 0.0_DP

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        dX(1) = dX1(iX1)

        dt_s = ( dX(1) * G(iNodeX,iX1,iX2,iX3,iGF_h_1) * CFL ) / G(iNodeX,iX1,iX2,iX3,iGF_Alpha)

        IF ( dt_s .LT. dt_min(1) ) THEN

          dt_min(1) = dt_s

        END IF

      END DO


    END DO
    END DO
    END DO

    IF( iX_B0(2) .LT. iX_E0(2) )THEN
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        dX(2) = dX2(iX2)
        dt_s = ( dX(2) * G(iNodeX,iX1,iX2,iX3,iGF_h_2) * CFL ) / G(iNodeX,iX1,iX2,iX3,iGF_Alpha)

        IF ( dt_s .LT. dt_min(2) ) THEN

          dt_min(2) = dt_s

        END IF

      END DO


    END DO
    END DO
    END DO
    END IF

    IF( iX_B0(3) .LT. iX_E0(3) )THEN
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        dX(3) = dX3(iX3)

        dt_s = ( dX(3) * G(iNodeX,iX1,iX2,iX3,iGF_h_3) * CFL ) / G(iNodeX,iX1,iX2,iX3,iGF_Alpha)
        IF ( dt_s .LT. dt_min(3) ) THEN

          dt_min(3) = dt_s

        END IF

      END DO


    END DO
    END DO
    END DO
    END IF
    dt = MINVAL( dt_min )
    END ASSOCIATE
  END SUBROUTINE CalculateTimeStep

  SUBROUTINE ComputeGray_TwoMoment_MF &
    ( MF_uGF, MF_uPF, MF_uCR, MF_uPR, MF_uGR )


    TYPE(amrex_multifab), INTENT(in)    :: &
      MF_uGF(0:nLevels-1), MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: &
      MF_uCR(0:nLevels-1), MF_uPR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uGR(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGR(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PF(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: CR(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: PR(:,:,:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: GR(:,:,:,:,:,:)

    INTEGER :: iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iLevel, iLo_MF(4)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )


        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uPR => MF_uPR(iLevel) % DataPtr( MFI )
        uGR => MF_uGR(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iZ_B0(1) = iE_B0
        iZ_E0(1) = iE_E0
        iZ_B1(1) = iE_B1
        iZ_E1(1) = iE_E1

        iZ_B0(2:4) = iX_B0
        iZ_E0(2:4) = iX_E0
        iZ_B1(2:4) = iX_B1
        iZ_E1(2:4) = iX_E1

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

        CALL AllocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 PF )

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
                 CR )

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
                   nPR     , &
                   nSpecies ], &
                 PR )

        CALL AllocateArray_Integrated &
               ( [ 1       , &
                   iX_B1(1), &
                   iX_B1(2), &
                   iX_B1(3), &
                   1       , &
                   1        ], &
                 [ nDOFX   , &
                   iX_E1(1), &
                   iX_E1(2), &
                   iX_E1(3), &
                   nGR     , &
                   nSpecies ], &
                 GR )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X &
               ( nPF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uPF, PF )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uCR, CR )

        CALL amrex2thornado_Z &
               ( nPR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uPR, PR )

        CALL amrex2thornado_Integrated &
               ( nGR, nSpecies, iX_B1, iX_E1, &
                 iLo_MF, iX_B0, iX_E0, uGR, GR )

        CALL ComputeGray_TwoMoment &
              ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
                G, PF, CR, PR, GR )

        CALL thornado2amrex_Integrated &
               ( nGR, nSpecies, iX_B1, iX_E1, &
                 iLo_MF, iX_B0, iX_E0, uGR, GR )

        CALL DeallocateArray_Integrated &
               ( [ 1       , &
                   iX_B1(1), &
                   iX_B1(2), &
                   iX_B1(3), &
                   1       , &
                   1        ], &
                 [ nDOFX   , &
                   iX_E1(1), &
                   iX_E1(2), &
                   iX_E1(3), &
                   nGR     , &
                   nSpecies ], &
                 GR )

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
                   nPR     , &
                   nSpecies ], &
                 PR )

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
                 CR )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nPF ], &
                 PF )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B1(1), iX_B1(2), iX_B1(3), 1   ], &
                 [ nDOFX, iX_E1(1), iX_E1(2), iX_E1(3), nGF ], &
                 G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO


  END SUBROUTINE ComputeGray_TwoMoment_MF


  SUBROUTINE ComputeGray_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GX, PF, CR, PR, GR )

    INTEGER,  INTENT(in)  :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)  :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in)  :: &
      PF(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nPF)
    REAL(DP), INTENT(in)  :: &
      CR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(in)  :: &
      PR(1:nDOFZ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      GR(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGR,1:nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iGR, iS, iNodeZ, iNodeE, iNodeX
    REAL(DP) :: hc3, E_0, E, RMS_Int3, RMS_Int5
    REAL(DP) :: W2(1:nDOFE,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: W3(1:nDOFE,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: W3_RMS(1:nDOFE,iZ_B0(1):iZ_E0(1))
    REAL(DP) :: W5_RMS(1:nDOFE,iZ_B0(1):iZ_E0(1))

    IF( UnitsActive )THEN

      hc3 = ( PlanckConstant * SpeedOfLight )**3

    ELSE

      hc3 = One

    END IF

    ! --- Integration Weights ---

    ASSOCIATE( dZ1 => MeshE % Width )

    E_0 = NodeCoordinate( MeshE, iZ_B0(1), 1 )

    DO iZ1    = iZ_B0(1), iZ_E0(1)
    DO iNodeE = 1, nDOFE

      E = NodeCoordinate( MeshE, iZ1, iNodeE )

      W2(iNodeE,iZ1) = FourPi * WeightsE(iNodeE) * ( dZ1(iZ1) * E**2 / hc3 )

      W3(iNodeE,iZ1) = W2(iNodeE,iZ1) * E

      W3_RMS(iNodeE,iZ1) = W2(iNodeE,iZ1) * ( E / E_0 )

      W5_RMS(iNodeE,iZ1) = W2(iNodeE,iZ1) * ( E / E_0 )**3



    END DO
    END DO
    END ASSOCIATE ! dZ1

    ! --- Initialize Gray Radiation Fields ---
    DO iS  = 1, nSpecies
    DO iGR = 1, nGR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        GR(iNodeX,iZ2,iZ3,iZ4,iGR,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Integrate Over Energy ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        RMS_Int3 = Zero
        RMS_Int5 = Zero

        DO iZ1    = iZ_B0(1), iZ_E0(1)
        DO iNodeE = 1, nDOFE

          iNodeZ = (iNodeX-1) * nDOFE + iNodeE

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_N,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_N,iS) &
                + W2(iNodeE,iZ1) * CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_D,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_D,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_I1,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_I1,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_I2,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_I2,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_I3,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_I3,iS) &
                + W2(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_J,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_J,iS) &
                + W3(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_H1,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_H1,iS) &
                + W3(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_H2,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_H2,iS) &
                + W3(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS)

          GR(iNodeX,iZ2,iZ3,iZ4,iGR_H3,iS) &
            = GR(iNodeX,iZ2,iZ3,iZ4,iGR_H3,iS) &
                + W3(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS)

          RMS_Int3 &
            = RMS_Int3 &
                + W3_RMS(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)

          RMS_Int5 &
            = RMS_Int5 &
                + W5_RMS(iNodeE,iZ1) * PR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)


        END DO
        END DO

        GR(iNodeX,iZ2,iZ3,iZ4,iGR_RMS,iS) &
          = E_0 * SQRT( RMS_Int5 / RMS_Int3 )


      END DO

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeGray_TwoMoment


END MODULE MF_TwoMoment_UtilitiesModule
