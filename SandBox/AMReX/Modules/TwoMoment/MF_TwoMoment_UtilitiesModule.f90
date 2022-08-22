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
    iPR_I1, &
    iPR_I2, &
    iPR_I3
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
  USE TwoMoment_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_TwoMoment
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic, &
    ComputeFromConserved_Euler_Relativistic
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE MeshModule, ONLY: &
    MeshX

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
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
    thornado2amrex_Z

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeTimeStep_TwoMoment_Fancy_MF
  PUBLIC :: ComputeTimeStep_TwoMoment_MF
  PUBLIC :: ComputeFromConserved_TwoMoment_MF

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

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        CALL amrex2thornado_X &
               ( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL CalculateTimeStep &
               ( iX_B1, iX_E1, iX_B0, iX_E0, CFL, G, TimeStep( iLevel) )

        TimeStepMin( iLevel ) = MIN( TimeStepMin( iLevel ), TimeStep( iLevel ) )

        DEALLOCATE( G )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStepMin, nLevels )

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

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER :: iLo_MF(4), iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    INTEGER :: iNX, iX1, iX2, iX3
    INTEGER :: iNZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER :: ErrorExists

    INTEGER, ALLOCATABLE :: iErr_Euler    (:,:,:,:)
    INTEGER, ALLOCATABLE :: iErr_TwoMoment(:,:,:,:,:,:)

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

        ALLOCATE( G (1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( CF(1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nCF) )

        ALLOCATE( PF(1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nPF) )

        ALLOCATE( CR(1:nDOFZ,iZ_B1(1):iZ_E1(1), &
                             iZ_B1(2):iZ_E1(2), &
                             iZ_B1(3):iZ_E1(3), &
                             iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )

        ALLOCATE( PR(1:nDOFZ,iZ_B1(1):iZ_E1(1), &
                             iZ_B1(2):iZ_E1(2), &
                             iZ_B1(3):iZ_E1(3), &
                             iZ_B1(4):iZ_E1(4),1:nPR,1:nSpecies) )

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

        ALLOCATE( iErr_Euler(1:nDOFX,iX_B0(1):iX_E0(1), &
                                     iX_B0(2):iX_E0(2), &
                                     iX_B0(3):iX_E0(3)) )

        ErrorExists = 0

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          iErr_Euler(iNX,iX1,iX2,iX3) = 0

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
                 iErr = iErr_Euler(iNX,iX1,iX2,iX3) )

          ErrorExists = ErrorExists + iErr_Euler(iNX,iX1,iX2,iX3)

        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .NE. 0 )THEN

          WRITE(*,*) 'ERROR: ComputeFromConserved_TwoMoment_MF'
          WRITE(*,*) 'iX_B0: ', iX_B0
          WRITE(*,*) 'iX_E0: ', iX_E0

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            IF( iErr_Euler(iNX,iX1,iX2,iX3) .NE. 0 )THEN

              WRITE(*,*) 'iNX, iX1, iX2, iX3 = ', iNX, iX1, iX2, iX3

              CALL DescribeError_Euler &
                ( iErr_Euler(iNX,iX1,iX2,iX3), &
                  Int_Option = [ iNX ], &
                  Real_Option = [ CF(iNX,iX1,iX2,iX3,iCF_D ), &
                                  CF(iNX,iX1,iX2,iX3,iCF_S1), &
                                  CF(iNX,iX1,iX2,iX3,iCF_S2), &
                                  CF(iNX,iX1,iX2,iX3,iCF_S3), &
                                  CF(iNX,iX1,iX2,iX3,iCF_E ), &
                                  CF(iNX,iX1,iX2,iX3,iCF_Ne), &
                                  G (iNX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                                  G (iNX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                                  G (iNX,iX1,iX2,iX3,iGF_Gm_dd_33) ] )

            END IF

          END DO
          END DO
          END DO
          END DO

        END IF

        ALLOCATE( iErr_TwoMoment(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                                         iZ_B0(2):iZ_E0(2), &
                                         iZ_B0(3):iZ_E0(3), &
                                         iZ_B0(4):iZ_E0(4),1:nSpecies) )
        ErrorExists = 0

        DO iS  = 1       , nSpecies
        DO iZ4 = iZ_B0(4), iZ_E0(4)
        DO iZ3 = iZ_B0(3), iZ_E0(3)
        DO iZ2 = iZ_B0(2), iZ_E0(2)
        DO iZ1 = iZ_B0(1), iZ_E0(1)
        DO iNZ = 1       , nDOFZ

          iErr_TwoMoment(iNZ,iZ1,iZ2,iZ3,iZ4,iS) = 0

          iNX = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

          CALL ComputePrimitive_TwoMoment &
               ( CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 CR(iNZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 PR(iNZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 PF(iNX    ,iZ2,iZ3,iZ4,iPF_V1), &
                 PF(iNX    ,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNX    ,iZ2,iZ3,iZ4,iPF_V3), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33), &
                 Zero, Zero, Zero, &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Alpha), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Beta_1), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Beta_2), &
                 G (iNX    ,iZ2,iZ3,iZ4,iGF_Beta_3), &
                 iErr = iErr_TwoMoment(iNZ,iZ1,iZ2,iZ3,iZ4,iS) )

          ErrorExists = ErrorExists + iErr_TwoMoment(iNZ,iZ1,iZ2,iZ3,iZ4,iS)

        END DO
        END DO
        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .NE. 0 )THEN

          WRITE(*,*) 'ERROR: ComputeFromConserved_TwoMoment_MF'
          WRITE(*,*) 'iZ_B0: ', iZ_B0
          WRITE(*,*) 'iZ_E0: ', iZ_E0

          DO iS  = 1       , nSpecies
          DO iZ4 = iZ_B0(4), iZ_E0(4)
          DO iZ3 = iZ_B0(3), iZ_E0(3)
          DO iZ2 = iZ_B0(2), iZ_E0(2)
          DO iZ1 = iZ_B0(1), iZ_E0(1)
          DO iNZ = 1       , nDOFZ

            IF( iErr_TwoMoment(iNZ,iZ1,iZ2,iZ3,iZ4,iS) .NE. 0 )THEN

              WRITE(*,*) 'iNZ, iZ1, iZ2, iZ3, iZ4, iS = ', &
                          iNZ, iZ1, iZ2, iZ3, iZ4, iS

              iNX = MOD( (iNZ-1) / nDOFE, nDOFX ) + 1

              WRITE(*,*) ' PF_D: ', PF(iNX,iZ2,iZ3,iZ4,iPF_D )
              WRITE(*,*) 'PF_V1: ', PF(iNX,iZ2,iZ3,iZ4,iPF_V1)
              WRITE(*,*) 'PF_V2: ', PF(iNX,iZ2,iZ3,iZ4,iPF_V2)
              WRITE(*,*) 'PF_V3: ', PF(iNX,iZ2,iZ3,iZ4,iPF_V3)
              WRITE(*,*) ' PF_E: ', PF(iNX,iZ2,iZ3,iZ4,iPF_E )
              WRITE(*,*) 'PF_Ne: ', PF(iNX,iZ2,iZ3,iZ4,iPF_Ne)

              CALL DescribeError_Euler &
                     ( iErr_TwoMoment(iNZ,iZ1,iZ2,iZ3,iZ4,iS) )

            END IF

          END DO
          END DO
          END DO
          END DO
          END DO
          END DO

        END IF

        CALL thornado2amrex_Z &
               ( nPR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B0, iZ_E0, uPR, PR )

        DEALLOCATE( iErr_TwoMoment )
        DEALLOCATE( iErr_Euler     )
        DEALLOCATE( G )
        DEALLOCATE( CF )
        DEALLOCATE( PF )
        DEALLOCATE( CR )
        DEALLOCATE( PR )

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

    IF( iX_B0(2) .GT. iX_E0(2) )THEN
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

    IF( iX_B0(3) .GT. iX_E0(3) )THEN
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


END MODULE MF_TwoMoment_UtilitiesModule
