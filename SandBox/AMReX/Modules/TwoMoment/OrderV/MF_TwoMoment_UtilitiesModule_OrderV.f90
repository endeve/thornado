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
    iGF_Gm_dd_33
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputePrimitive_TwoMoment
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic, &
    ComputeFromConserved_Euler_NonRelativistic
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
  USE TwoMoment_ClosureModule, ONLY: &
    FluxFactor, &
    EddingtonFactor, &
    HeatFluxFactor

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


  PUBLIC :: ComputeFromConserved_TwoMoment_MF
  


CONTAINS




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
    REAL(DP) :: PR_D, PR_I1, PR_I2, PR_I3

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iLo_MF(4), iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    INTEGER :: iNX, iX1, iX2, iX3
    INTEGER :: iNZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER :: ErrorExists
    INTEGER :: nIterations

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

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          CALL ComputePrimitive_Euler_NonRelativistic &
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
                 G (iNX,iX1,iX2,iX3,iGF_Gm_dd_33) )

        END DO
        END DO
        END DO
        END DO

        DO iS  = 1       , nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
        DO iZ1 = iZ_B0(1), iZ_E0(1)
        DO iNZ = 1       , nDOFZ

          iNX = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

          CALL ComputePrimitive_TwoMoment &
               ( CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 PR_D, PR_I1, PR_I2, PR_I3, &
                 PF(iNX    ,iZ2,iZ3,iZ4,iPF_V1), &
                 PF(iNX    ,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNX    ,iZ2,iZ3,iZ4,iPF_V3), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33), &
                 nIterations_Option = nIterations )

          PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) = PR_D
          PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) = PR_I1
          PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) = PR_I2
          PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) = PR_I3

          IF( nIterations .GE. 100 )THEN

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

 
END MODULE MF_TwoMoment_UtilitiesModule
