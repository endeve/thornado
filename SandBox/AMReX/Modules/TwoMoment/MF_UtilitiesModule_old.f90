MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    AR => amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_parallel_module,             ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_geometry_module,             ONLY: &
    amrex_geometry
  ! --- thornado Modules ---
  USE ProgramHeaderModule, ONLY: &
    nDOFZ,                 &
    nDOFX,                 &
    nNodesX,               &
    nX,                     &
    nE,                     &
    swX,                    &
    nNodesX,                &
    nDOFE,                  &
    swE,                    &
    nNodesE,                &
    iE_E0,                  &
    iE_B0,                  &
    iE_E1,                  &
    iE_B1
  USE RadiationFieldsModule,       ONLY: &
    nCR,    &
    nPR,    &
    iPR_D,  &
    iPR_I1, &
    iPR_I2, &
    iPR_I3, &
    iCR_N,  &
    iCR_G1, &
    iCR_G2, &
    iCR_G3
  USE MeshModule,                        ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule,              ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Alpha
  USE FluidFieldsModule,       ONLY: &
    nCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iPF_Ne, &
    nPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iCF_Ne, &
    nAF,    &
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_S, &
    iAF_E, &
    iAF_Gm, &
    iAF_Cs, &
    iAF_Me, &
    iAF_Mp, &
    iAF_Mn, &
    iAF_Xp, &
    iAF_Xn, &
    iAF_Xa, &
    iAF_Xh, &
    unitsPF, &
    unitsAF
  USE TwoMoment_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_TwoMoment
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE GeometryFieldsModuleE,     ONLY: &
    nGE, uGE
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ComputeThermodynamicStates_Auxiliary_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ApplyEquationOfState_TABLE
  USE UnitsModule, ONLY: &
    MeV
  ! --- Local Modules ---
  USE InputParsingModule, ONLY: &
    nLevels, &
    nSpecies, &
    UseTiling

  USE MF_TwoMoment_BoundaryConditionsModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap, &
    ApplyBoundaryConditions_TwoMoment_MF


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex2thornado_Z
  PUBLIC :: thornado2amrex_Z
  PUBLIC :: amrex2thornado_X
  PUBLIC :: thornado2amrex_X
  PUBLIC :: WriteNodalDataToFile
  PUBLIC :: WriteEulertoFile

  CONTAINS

  SUBROUTINE amrex2thornado_Z &
    ( nVars, nS, nE, iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF, &
      iZ_B, iZ_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nVars, nS, nE
    INTEGER,  INTENT(in)  :: iE_B0, iE_E0, iZ_B1(4), iZ_E1(4), iLo_MF(4), &
                             iZ_B(4), iZ_E(4)
    REAL(AR), INTENT(in)  :: &
      Data_amrex   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(AR), INTENT(out) :: &
      Data_thornado(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iVar, iD, iNodeZ

    DO iS = 1, nS
    DO iVar = 1, nVars
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)

      DO iZ1 = iE_B0, iE_E0 ! always want iZ1 to not include ghost cells
      DO iNodeZ = 1, nDOFZ

        iD = ( iS - 1 ) * nVars * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                + ( iVar - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                + ( iZ1 - 1 ) * nDOFZ + iNodeZ

        Data_thornado(iNodeZ,iZ1,iZ2,iZ3,iZ4,iVar,iS) &
          = Data_amrex(iZ2,iZ3,iZ4,iD)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE amrex2thornado_Z


  SUBROUTINE thornado2amrex_Z &
    ( nVars, nS, nE, iE_B0, iE_E0, iZ_B1, iZ_E1, iLo_MF, &
      iZ_B, iZ_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nVars, nS, nE
    INTEGER,  INTENT(in)  :: iE_B0, iE_E0, iZ_B1(4), iZ_E1(4), iLo_MF(4), &
                             iZ_B(4), iZ_E(4)
    REAL(AR), INTENT(out) :: &
      Data_amrex   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(AR), INTENT(in)  :: &
      Data_thornado(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    INTEGER :: iZ1, iZ2, iZ3, iZ4, iS, iVar, iNodeZ, iD

    DO iS = 1, nS
    DO iVar = 1, nVars
    DO iZ4 = iZ_B(4), iZ_E(4)
    DO iZ3 = iZ_B(3), iZ_E(3)
    DO iZ2 = iZ_B(2), iZ_E(2)

      DO iZ1 = iE_B0, iE_E0
      DO iNodeZ = 1, nDOFZ

        iD = ( iS - 1 ) * nVars * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
               + ( iVar - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
               + ( iZ1 - 1 ) * nDOFZ + iNodeZ

        Data_amrex(iZ2,iZ3,iZ4,iD) &
          = Data_thornado(iNodeZ,iZ1,iZ2,iZ3,iZ4,iVar,iS)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE thornado2amrex_Z


  SUBROUTINE amrex2thornado_X &
    ( nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(AR), INTENT(in)  :: &
      Data_amrex   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(AR), INTENT(out) :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_thornado(1:nDOFX,iX1,iX2,iX3,iFd) &
        = Data_amrex(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE amrex2thornado_X


  SUBROUTINE thornado2amrex_X &
    ( nFields, iX_B1, iX_E1, iLo_MF, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B1(3), iX_E1(3), iLo_MF(4), iX_B(3), iX_E(3)
    REAL(AR), INTENT(out) :: &
      Data_amrex   (   iLo_MF(1):,iLo_MF(2):,iLo_MF(3):,iLo_MF(4):)
    REAL(AR), INTENT(in)  :: &
      Data_thornado(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3, iFd

    DO iFd = 1, nFields
    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      Data_amrex(iX1,iX2,iX3,nDOFX*(iFd-1)+1:nDOFX*iFd) &
        = Data_thornado(1:nDOFX,iX1,iX2,iX3,iFd)

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE thornado2amrex_X


  SUBROUTINE WriteNodalDataToFile( GEOM, MF_uGF, MF_uCF, MF_uCR, FileNameBase )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCR(0:nLevels-1)
    CHARACTER(LEN=*)    , INTENT(in) :: FileNameBase

    INTEGER            :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3),     &
                          iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4),     &
                          iZ_B (4), iZ_E (4), iE_B, iE_E,             &
                          iX_B (3), iX_E (3), iEL, iER, iLo_MF(4),    &
                          iLevel, nCompGF, nCompCF, nCompCR,          &
                          iX1, iX2, iX3, iCF, iGF, iCR, i,            &
                          iZ1, iZ2, iZ3, iZ4, iS, iE,                 &
                          iNodeZ,                                     &
                          iNodeX, iNodeX1, iNodeX2, iNodeX3, iNodeE
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(EdgeMap)      :: Edge_Map

    REAL(AR), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF (:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCR (:,:,:,:)

    REAL(AR), ALLOCATABLE :: G (:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: CF (:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: PF (:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: CR (:,:,:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: PR (:,:,:,:,:,:,:)

    REAL(AR) :: D, I1, I2, I3, N, G1, G2, G3

    ALLOCATE( G(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nGF) )

    ALLOCATE( CF(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nCF) )

    ALLOCATE( PF(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nPF) )
    ALLOCATE( CR(1:nDOFZ,1-swE:nE+swE, &
                        1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nCR,1:nSpecies) )
    ALLOCATE( PR(1:nDOFZ,1-swE:nE+swE, &
                        1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nPR,1:nSpecies) )
    G = 0.0_AR
    CF = 0.0_AR
    PF = 0.0_AR
    CR = 0.0_AR
    PR = 0.0_AR

    ! --- Convert AMReX MultiFabs to thornado arrays ---

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCR(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF     => MF_uGF(iLevel) % DataPtr( MFI )
        nCompGF =  MF_uGF(iLevel) % nComp()

        uCF     => MF_uCF(iLevel) % DataPtr( MFI )
        nCompCF =  MF_uCF(iLevel) % nComp()

        uCR     => MF_uCR(iLevel) % DataPtr( MFI )
        nCompCR =  MF_uCr(iLevel) % nComp()

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iX_B = iX_B0
        iX_E = iX_E0

        i=1

        DO WHILE (i<=4)

          IF (i==1) THEN

            iZ_B0(i)=iE_B0
            iZ_E0(i)=iE_E0
          ELSE

            iZ_B0(i)=iX_B0(i-1)
            iZ_E0(i)=iX_E0(i-1)
          END IF
          i = i + 1
        END DO

        IF( iX_B0(1) .EQ. 1     ) iX_B(1) = 1     - swX(1)
        IF( iX_B0(2) .EQ. 1     ) iX_B(2) = 1     - swX(2)
        IF( iX_B0(3) .EQ. 1     ) iX_B(3) = 1     - swX(3)
        IF( iX_E0(1) .EQ. nX(1) ) iX_E(1) = nX(1) + swX(1)
        IF( iX_E0(2) .EQ. nX(2) ) iX_E(2) = nX(2) + swX(2)
        IF( iX_E0(3) .EQ. nX(3) ) iX_E(3) = nX(3) + swX(3)
        IF( iE_B0 .EQ. 1     ) iE_B = 1     - swE
        IF( iE_E0 .EQ. nE ) iE_E = nE + swE


        i=1

        DO WHILE (i<=4)

          IF (i==1) THEN

            iZ_B(i)=iE_B
            iZ_E(i)=iE_E
            iZ_B1(i)=iE_B1
            iZ_E1(i)=iE_E1
          ELSE

            iZ_B(i)=iX_B(i-1)
            iZ_E(i)=iX_E(i-1)
            iZ_B1(i)=iX_B1(i-1)
            iZ_E1(i)=iX_E1(i-1)
          END IF
          i = i + 1
        END DO

        CALL amrex2thornado_X( nGF, iX_B, iX_E, iLo_MF, iX_B, iX_E, uGF, &
                               G(1:nDOFX,iX_B(1):iX_E(1), &
                                         iX_B(2):iX_E(2), &
                                         iX_B(3):iX_E(3),1:nGF) )

        CALL amrex2thornado_X( nCF, iX_B, iX_E, iLo_MF, iX_B, iX_E, uCF, &
                               CF(1:nDOFX,iX_B(1):iX_E(1), &
                                          iX_B(2):iX_E(2), &
                                          iX_B(3):iX_E(3),1:nCF) )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B1, iZ_E1, iLo_MF, iZ_B1, iZ_E1, uCR, &
                 CR(1:nDOFZ,iZ_B1(1):iZ_E1(1), &
                            iZ_B1(2):iZ_E1(2), &
                            iZ_B1(3):iZ_E1(3), &
                            iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies) )

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )



        CALL ApplyBoundaryConditions_TwoMoment_MF &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, CR, Edge_Map )

      END DO

    END DO

    DO iX3    = 1-swX(3), nX(3)+swX(3)
    DO iX2    = 1-swX(2), nX(2)+swX(2)
    DO iX1    = 1-swX(1), nX(1)+swX(1)
    DO iNodeX = 1       , nDOFX

      DO iGF = 1, nGF

        CALL amrex_parallel_reduce_sum( G(iNodeX,iX1,iX2,iX3,iGF) )

      END DO

      DO iCF = 1, nCF

        CALL amrex_parallel_reduce_sum( CF(iNodeX,iX1,iX2,iX3,iCF) )

      END DO

    END DO
    END DO
    END DO
    END DO

    DO iS     = 1, nSpecies
    DO iZ4    = 1-swX(3), nX(3)+swX(3)
    DO iZ3    = 1-swX(2), nX(2)+swX(2)
    DO iZ2    = 1-swX(1), nX(1)+swX(1)
    DO iZ1    = 1-swE, nE+swE
    DO iNodeZ = 1       , nDOFZ


      DO iCR = 1, nCR

        CALL amrex_parallel_reduce_sum( CR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR,iS) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
      DO iX1 = 1, nX(1)

        CALL ComputePrimitive_Euler_Relativistic &
               ( CF(:,iX1,iX2,iX3,iCF_D ),       &
                 CF(:,iX1,iX2,iX3,iCF_S1),       &
                 CF(:,iX1,iX2,iX3,iCF_S2),       &
                 CF(:,iX1,iX2,iX3,iCF_S3),       &
                 CF(:,iX1,iX2,iX3,iCF_E ),       &
                 CF(:,iX1,iX2,iX3,iCF_Ne),       &
                 PF(:,iX1,iX2,iX3,iPF_D ),       &
                 PF(:,iX1,iX2,iX3,iPF_V1),       &
                 PF(:,iX1,iX2,iX3,iPF_V2),       &
                 PF(:,iX1,iX2,iX3,iPF_V3),       &
                 PF(:,iX1,iX2,iX3,iPF_E ),       &
                 PF(:,iX1,iX2,iX3,iPF_Ne),       &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )


      END DO
      END DO
      END DO

    END IF

    IF( amrex_parallel_ioprocessor() )THEN

      DO iS = 1, nSpecies
      DO iZ4 = 1, nX(3)
      DO iZ3 = 1, nX(2)
      DO iZ2 = 1, nX(1)
      DO iZ1 = 1, nE

        DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE

          iNodeZ = (iNodeX-1) * nDOFE + iNodeE
          CALL ComputePrimitive_TwoMoment &
               ( CR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_N, iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G1, iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G2, iS), &
                 CR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iCR_G3, iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iPR_D, iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iPR_I1, iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iPR_I2, iS), &
                 PR(iNodeZ,iZ1,iZ2,iZ3,iZ4, iPR_I3, iS), &
                 PF(iNodeX,iZ2,iZ3,iZ4, iPF_V1), &
                 PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2), &
                 PF(iNodeX,iZ2,iZ3,iZ4, iPF_V3), &
                 G (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11), &
                 G (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22), &
                 G (iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33), &
                 0.0_AR, 0.0_AR, 0.0_AR,                &
                 G(iNodeX ,iZ2,iZ3,iZ4,iGF_Alpha  ), &
                 G(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_1  ), &
                 G(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_2  ), &
                 G(iNodeX  ,iZ2,iZ3,iZ4,iGF_Beta_3  ) )

        END DO
        END DO

      END DO
      END DO
      END DO
      END DO
      END DO

      OPEN( UNIT = 100, FILE = TRIM( FileNameBase ) // 'r.dat'      )
      OPEN( UNIT = 101, FILE = TRIM( FileNameBase ) // 'D.dat'      )
      OPEN( UNIT = 102, FILE = TRIM( FileNameBase ) // 'I1.dat'  )
      OPEN( UNIT = 103, FILE = TRIM( FileNameBase ) // 'I2.dat'    )
      OPEN( UNIT = 104, FILE = TRIM( FileNameBase ) // 'I3.dat' )
      OPEN( UNIT = 105, FILE = TRIM( FileNameBase ) // 'N.dat'      )
      OPEN( UNIT = 106, FILE = TRIM( FileNameBase ) // 'G1.dat'      )
      OPEN( UNIT = 107, FILE = TRIM( FileNameBase ) // 'G2.dat'      )
      OPEN( UNIT = 108, FILE = TRIM( FileNameBase ) // 'G3.dat'      )
      OPEN( UNIT = 109, FILE = TRIM( FileNameBase ) // 'D_middle.dat'      )
      OPEN( UNIT = 110, FILE = TRIM( FileNameBase ) // 'D_spatial.dat'      )


      ! --- Hacked to work only for 1D and 2D problems ---
      DO iS = 1, nSpecies
      DO iX3     = 1, nX(3)
      DO iNodeX3 = 1, nNodesX(3)

        DO iX2     = 1, nX(2)
        DO iNodeX2 = 1, nNodesX(2)

        DO iX1     = 1, nX(1)
        DO iNodeX1 = 1, nNodesX(1)
          IF     ( iNodeX1 .EQ. 1 )THEN

            iEL = 1
            iER = nNodesE

          ELSE IF( iNodeX1 .EQ. 2 )THEN

            iEL = nNodesE + 1
            iER = 2 * nNodesE

          ELSE IF( iNodeX1 .EQ. 3 )THEN

            iEL = 2 * nNodesE + 1
            iER = 3 * nNodesE

          END IF

          DO iE = 1, nE
          DO iNodeE = iEL, iER
            D = PR(iNodeE,iE,iX1,iX2,iX3,iPR_D,iS)
            I1 = PR(iNodeE,iE,iX1,iX2,iX3,iPR_I1,iS)
            I2 = PR(iNodeE,iE,iX1,iX2,iX3,iPR_I2,iS)
            I3 = PR(iNodeE,iE,iX1,iX2,iX3,iPR_I3,iS)
            N = CR(iNodeE,iE,iX1,iX2,iX3,iCR_N,iS)
            G1 = CR(iNodeE,iE,iX1,iX2,iX3,iCR_G1,iS)
            G2 = CR(iNodeE,iE,iX1,iX2,iX3,iCR_G2,iS)
            G3 = CR(iNodeE,iE,iX1,iX2,iX3,iCR_G3,iS)

            WRITE(101,'(ES24.16E3,1x)',ADVANCE='NO') &
              D
            WRITE(102,'(ES24.16E3,1x)',ADVANCE='NO') &
              I1
            WRITE(103,'(ES24.16E3,1x)',ADVANCE='NO') &
              I2
            WRITE(104,'(ES24.16E3,1x)',ADVANCE='NO') &
              I3
            WRITE(105,'(ES24.16E3,1x)',ADVANCE='NO') &
              N
            WRITE(106,'(ES24.16E3,1x)',ADVANCE='NO') &
              G1
            WRITE(107,'(ES24.16E3,1x)',ADVANCE='NO') &
              G2
            WRITE(108,'(ES24.16E3,1x)',ADVANCE='NO') &
              G3
            IF (iX1 .EQ. nX(1)/2 .AND. iNodeX1 .EQ. 1) THEN
              WRITE(109,'(ES24.16E3,1x)',ADVANCE='NO') &
                D
            END IF
            IF (iE .EQ. 1 .AND. iNodeE .EQ. 1) THEN
              WRITE(110,'(ES24.16E3,1x)',ADVANCE='NO') &
                D
            END IF
          END DO
          END DO
        WRITE(101,*)
        WRITE(102,*)
        WRITE(103,*)
        WRITE(104,*)
        WRITE(105,*)
        WRITE(106,*)
        WRITE(107,*)
        WRITE(108,*)
        END DO
        END DO
        END DO
        END DO


      END DO
      END DO
      END DO

      CLOSE( 110 )
      CLOSE( 109 )
      CLOSE( 108 )
      CLOSE( 107 )
      CLOSE( 106 )
      CLOSE( 105 )
      CLOSE( 104 )
      CLOSE( 103 )
      CLOSE( 102 )
      CLOSE( 101 )
      CLOSE( 100 )

    END IF
    DEALLOCATE( CF )
    DEALLOCATE( PF )
    DEALLOCATE( CR )
    DEALLOCATE( PR )
    DEALLOCATE( G )

    print*, "Writing Nodal Values"
  END SUBROUTINE WriteNodalDataToFile


  SUBROUTINE WriteEulertoFile( MF_uCF, MF_uGF, n )

    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    INTEGER, INTENT(in) :: n

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)


    REAL(AR), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: AF(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: PF(:,:,:,:,:)

    REAL(AR) :: Mu(1:nDOFX), num

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)
    INTEGER :: iX1, iX2, iX3, l

    character(len=64) :: nm
    CHARACTER(LEN=16) :: FMT

    WRITE(FMT,'(A3,I3.3,A10)') '(SP', nDOFX, 'ES25.16E3)'
    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX


        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF) )

        ALLOCATE( PF(1:nDOFX,iX_B1(1):iX_E1(1), &
                             iX_B1(2):iX_E1(2), &
                             iX_B1(3):iX_E1(3),1:nPF) )

        ALLOCATE( AF(1:nDOFX,iX_B1(1):iX_E1(1), &
                              iX_B1(2):iX_E1(2), &
                              iX_B1(3):iX_E1(3),1:nAF) )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, U )


        DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          CALL ComputePrimitive_Euler_Relativistic &
               ( U(:,iX1,iX2,iX3,iCF_D ),       &
                 U(:,iX1,iX2,iX3,iCF_S1),       &
                 U(:,iX1,iX2,iX3,iCF_S2),       &
                 U(:,iX1,iX2,iX3,iCF_S3),       &
                 U(:,iX1,iX2,iX3,iCF_E ),       &
                 U(:,iX1,iX2,iX3,iCF_Ne),       &
                 PF(:,iX1,iX2,iX3,iPF_D ),       &
                 PF(:,iX1,iX2,iX3,iPF_V1),       &
                 PF(:,iX1,iX2,iX3,iPF_V2),       &
                 PF(:,iX1,iX2,iX3,iPF_V3),       &
                 PF(:,iX1,iX2,iX3,iPF_E ),       &
                 PF(:,iX1,iX2,iX3,iPF_Ne),       &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

           CALL ComputeAuxiliary_Fluid_TABLE &
               ( PF(:,iX1,iX2,iX3,iPF_D ), &
                 PF(:,iX1,iX2,iX3,iPF_E ), &
                 PF(:,iX1,iX2,iX3,iPF_Ne), &
                 AF(:,iX1,iX2,iX3,iAF_P ), &
                 AF(:,iX1,iX2,iX3,iAF_T ), &
                 AF(:,iX1,iX2,iX3,iAF_Ye), &
                 AF(:,iX1,iX2,iX3,iAF_S ), &
                 AF(:,iX1,iX2,iX3,iAF_E ), &
                 AF(:,iX1,iX2,iX3,iAF_Gm), &
                 AF(:,iX1,iX2,iX3,iAF_Cs) )



          CALL ApplyEquationOfState_TABLE &
               ( PF(:,iX1,iX2,iX3,iPF_D ), &
                 AF(:,iX1,iX2,iX3,iAF_T ), &
                 AF(:,iX1,iX2,iX3,iAF_Ye), &
                 AF(:,iX1,iX2,iX3,iAF_P ), &
                 AF(:,iX1,iX2,iX3,iAF_S ), &
                 AF(:,iX1,iX2,iX3,iAF_E ), &
                 AF(:,iX1,iX2,iX3,iAF_Me), &
                 AF(:,iX1,iX2,iX3,iAF_Mp), &
                 AF(:,iX1,iX2,iX3,iAF_Mn), &
                 AF(:,iX1,iX2,iX3,iAF_Xp), &
                 AF(:,iX1,iX2,iX3,iAF_Xn), &
                 AF(:,iX1,iX2,iX3,iAF_Xa), &
                 AF(:,iX1,iX2,iX3,iAF_Xh), &
                 AF(:,iX1,iX2,iX3,iAF_Gm) )


        END DO
        END DO
        END DO
        WRITE(nm,*) n

        nm = ADJUSTL(nm)

        OPEN( UNIT = 101, FILE = 'DF.dat'      )
        OPEN( UNIT = 102, FILE = 'V1.dat'      )
        OPEN( UNIT = 103, FILE = 'V2.dat'      )
        OPEN( UNIT = 104, FILE = 'V3.dat'      )
        OPEN( UNIT = 105, FILE = 'Ye.dat'      )
        OPEN( UNIT = 106, FILE = 'Mu'//trim(nm) //'.dat'      )
        OPEN( UNIT = 107, FILE = 'E.dat'      )
        OPEN( UNIT = 108, FILE = 'T'//trim(nm) //'.dat'      )


        DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          Mu(:) = AF(:,iX1,iX2,iX3,iAF_Me) + AF(:,iX1,iX2,iX3,iAF_Mp) - AF(:,iX1,iX2,iX3,iAF_Mn)


          WRITE(101,FMT) &
            PF(:,iX1,iX2,iX3,iPF_D) &
              / unitsPF(iPF_D)

          WRITE(102,FMT) &
            PF(:,iX1,iX2,iX3,iPF_V1) &
              / unitsPF(iPF_V1)

          WRITE(103,FMT) &
            PF(:,iX1,iX2,iX3,iPF_V2) &
              / unitsPF(iPF_V2)

          WRITE(104,FMT) &
            PF(:,iX1,iX2,iX3,iPF_V3) &
              / unitsPF(iPF_V3)

          WRITE(105,FMT) &
            AF(:,iX1,iX2,iX3,iAF_Ye) &
              / unitsAF(iAF_Ye)

          WRITE(106,FMT) &
            Mu(:) &
              / unitsAF(iAF_Me)

          WRITE(107,FMT) &
            PF(:,iX1,iX2,iX3,iPF_E) &
              / unitsPF(iPF_E)

          WRITE(108,FMT) &
            AF(:,iX1,iX2,iX3,iAF_T) &
              / MeV
        END DO
        END DO
        END DO
        DEALLOCATE( U )
        DEALLOCATE( G )
        DEALLOCATE( PF )
        DEALLOCATE( AF )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO







  END SUBROUTINE WriteEulertoFile


END MODULE MF_UtilitiesModule
