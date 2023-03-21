MODULE MF_TwoMoment_TallyModule

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
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE UnitsModule, ONLY: &
    UnitsActive, &
    SpeedOfLight, &
    PlanckConstant, &
    MeV,            &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX, &
    nDimsX
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE ReferenceElementModule, ONLY: &
    Weights_q
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE MeshModule, ONLY: &
    MeshType, &
    CreateMesh, &
    DestroyMesh
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_SqrtGm, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    CoordinateSystem
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_E
  USE GeometryFieldsModuleE,     ONLY: &
    nGE, uGE, iGE_Ep2, iGE_Ep3
  USE RadiationFieldsModule, ONLY: &
    nSpecies, LeptonNumber, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE ProgramHeaderModule,      ONLY: &
    swX, nDOFX, nDOFZ, swE, nDOFE, iE_B0, iE_E0, iE_B1, iE_E1, nNodesE
  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    FourPi, &
    One,   &
    Half
  USE InputParsingModule, ONLY: &
    nX, &
    nLevels, &
    ProgramName, &
    xL, &
    xR, &
    eL, &
    eR, &
    nE, &
    UseTiling, &
    zoomE
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_Z, &
    AllocateArray_X, &
    DeallocateArray_X, &
    AllocateArray_Z, &
    DeallocateArray_Z
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF




  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally_TwoMoment_MF
  PUBLIC :: ComputeTally_TwoMoment_MF
  PUBLIC :: FinalizeTally_TwoMoment_MF
  PUBLIC :: IncrementOffGridTally_TwoMoment_MF

  LOGICAL :: SuppressTally


  REAL(DP) :: hc3


  CHARACTER(256) :: NeutrinoLeptonNumber_FileName
  REAL(DP), ALLOCATABLE :: NeutrinoLeptonNumber_Interior(:)
  REAL(DP), ALLOCATABLE :: NeutrinoLeptonNumber_Initial(:)
  REAL(DP), ALLOCATABLE :: NeutrinoLeptonNumber_OffGrid(:)
  REAL(DP), ALLOCATABLE :: NeutrinoLeptonNumber_Change(:)

  CHARACTER(256) :: NeutrinoEnergy_FileName
  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_Interior(:)
  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_Initial(:)
  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_Offgrid(:)
  REAL(DP), ALLOCATABLE :: NeutrinoEnergy_Change(:)


  CHARACTER(256) :: Momentum_FileName
  REAL(DP), ALLOCATABLE :: Momentum_X1(:)
  REAL(DP), ALLOCATABLE :: Momentum_X2(:)
  REAL(DP), ALLOCATABLE :: Momentum_X3(:)

CONTAINS


  SUBROUTINE InitializeTally_TwoMoment_MF &
    ( SuppressTally_Option, BaseFileName_Option )

    LOGICAL,  INTENT(in),         OPTIONAL :: &
      SuppressTally_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      BaseFileName_Option


    CHARACTER(256) :: BaseFileName
    INTEGER        :: FileUnit

    CHARACTER(256) :: TimeLabel
    CHARACTER(256) :: LeptonNumber_InteriorLabel
    CHARACTER(256) :: LeptonNumber_InitialLabel
    CHARACTER(256) :: LeptonNumber_OffgridLabel
    CHARACTER(256) :: LeptonNumber_ChangeLabel

    CHARACTER(256) :: Energy_InteriorLabel
    CHARACTER(256) :: Energy_InitialLabel
    CHARACTER(256) :: Energy_OffGridLabel
    CHARACTER(256) :: Energy_ChangeLabel

    CHARACTER(256) :: Momentum1Label
    CHARACTER(256) :: Momentum2Label
    CHARACTER(256) :: Momentum3Label

    SuppressTally = .FALSE.
    IF( PRESENT( SuppressTally_Option ) ) &
      SuppressTally = SuppressTally_Option

    IF( SuppressTally ) RETURN


    IF( UnitsActive )THEN

      hc3 = ( PlanckConstant * SpeedOfLight )**3

    ELSE

      hc3 = One

    END IF



    ALLOCATE(NeutrinoLeptonNumber_Interior(0:nLevels-1) )
    ALLOCATE(NeutrinoLeptonNumber_Initial (0:nLevels-1) )
    ALLOCATE(NeutrinoLeptonNumber_OffGrid (0:nLevels-1) )
    ALLOCATE(NeutrinoLeptonNumber_Change  (0:nLevels-1) )

    ALLOCATE( NeutrinoEnergy_Interior(0:nLevels-1) )
    ALLOCATE( NeutrinoEnergy_Initial (0:nLevels-1) )
    ALLOCATE( NeutrinoEnergy_Offgrid (0:nLevels-1) )
    ALLOCATE( NeutrinoEnergy_Change  (0:nLevels-1) )

    ALLOCATE( Momentum_X1(0:nLevels-1) )
    ALLOCATE( Momentum_X2(0:nLevels-1) )
    ALLOCATE( Momentum_X3(0:nLevels-1) )

    IF( amrex_parallel_ioprocessor() )THEN

      BaseFileName = ''
      IF( PRESENT( BaseFileName_Option ) ) &
        BaseFileName = TRIM( BaseFileName_Option )

      BaseFileName = TRIM( BaseFileName ) // TRIM( ProgramName )

      ! --- Neutrino Lepton Number ---

      NeutrinoLeptonNumber_FileName &
        = TRIM( BaseFileName ) // '.Tally_NeutrinoLeptonNumber.dat'

      TimeLabel     &
        = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
      LeptonNumber_InteriorLabel &
        = 'LeptonNumber_Interior'
      LeptonNumber_OffgridLabel &
        = 'LeptonNumber_OffGrid'
      LeptonNumber_InitialLabel &
        = 'LeptonNumber_Initial'
      LeptonNumber_ChangeLabel &
        = 'LeptonNumber_Change'


      OPEN( NEWUNIT = FileUnit, FILE = TRIM(NeutrinoLeptonNumber_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( LeptonNumber_InteriorLabel ), TRIM( LeptonNumber_OffGridLabel ), &
        TRIM( LeptonNumber_InitialLabel ), TRIM( LeptonNumber_ChangeLabel )
      CLOSE( FileUnit )

      ! --- Neutrino Energy  ---

      NeutrinoEnergy_FileName &
        = TRIM( BaseFileName ) // '.Tally_NeutrinoEnergy.dat'

      TimeLabel     &
        = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
      Energy_InteriorLabel &
        = 'Energy_Interior [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      Energy_OffGridLabel &
        = 'Energy_OffGrid [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      Energy_InitialLabel &
        = 'Energy_Initial [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'
      Energy_ChangeLabel &
        = 'Energy_Change [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'


      OPEN( NEWUNIT = FileUnit, FILE = TRIM(NeutrinoEnergy_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( Energy_InteriorLabel ), TRIM( Energy_OffGridLabel ), &
        TRIM( Energy_InitialLabel ), TRIM( Energy_ChangeLabel )
      CLOSE( FileUnit )


      ! --- Momentum ---

      Momentum_FileName &
        = TRIM( BaseFileName ) // '.Tally_Momentum.dat'
      Momentum1Label &
        = 'Momentum_1'
      Momentum2Label &
        = 'Momentum_2'
      Momentum3Label &
        = 'Momentum_3'

      OPEN( NEWUNIT = FileUnit, FILE = TRIM(Momentum_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), &
        TRIM( Momentum1Label ), TRIM( Momentum2Label ), TRIM( Momentum3Label )

      CLOSE( FileUnit )

    END IF

    NeutrinoLeptonNumber_Interior = Zero
    NeutrinoLeptonNumber_Initial  = Zero
    NeutrinoLeptonNumber_OffGrid  = Zero
    NeutrinoLeptonNumber_Change   = Zero


    NeutrinoEnergy_Interior       = Zero
    NeutrinoEnergy_Initial        = Zero
    NeutrinoEnergy_OffGrid        = Zero
    NeutrinoEnergy_Change         = Zero



    Momentum_X1                   = Zero
    Momentum_X2                   = Zero
    Momentum_X3                   = Zero

  END SUBROUTINE InitializeTally_TwoMoment_MF



  SUBROUTINE FinalizeTally_TwoMoment_MF

    IF( SuppressTally ) RETURN

    DEALLOCATE( NeutrinoLeptonNumber_Interior )
    DEALLOCATE( NeutrinoLeptonNumber_Initial )
    DEALLOCATE( NeutrinoLeptonNumber_OffGrid )
    DEALLOCATE( NeutrinoLeptonNumber_Change  )

    DEALLOCATE( NeutrinoEnergy_Interior )
    DEALLOCATE( NeutrinoEnergy_Initial  )
    DEALLOCATE( NeutrinoEnergy_Offgrid  )
    DEALLOCATE( NeutrinoEnergy_Change   )

    DEALLOCATE( Momentum_X1 )
    DEALLOCATE( Momentum_X2 )
    DEALLOCATE( Momentum_X3 )

  END SUBROUTINE FinalizeTally_TwoMoment_MF


  SUBROUTINE ComputeTally_TwoMoment_MF &
    ( GEOM, MF_uGF, MF_uCF, MF_uCR, Time, SetInitialValues_Option, Verbose_Option )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCR(0:nLevels-1)


    REAL(DP),             INTENT(in) :: Time
    LOGICAL,              INTENT(in), OPTIONAL :: SetInitialValues_Option
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option


    LOGICAL :: SetInitialValues
    LOGICAL :: Verbose

    INTEGER                       :: iX_B0(3), iX_E0(3), iZ_B0(4), iZ_E0(4)
    INTEGER                       :: iLevel, iLo_MF(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: UF(:,:,:,:,:)
    REAL(DP), ALLOCATABLE         :: U (:,:,:,:,:,:,:)


    IF( SuppressTally ) RETURN


    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) ) &
      SetInitialValues = SetInitialValues_Option



    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option



    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      NeutrinoLeptonNumber_Interior(iLevel) = Zero

      NeutrinoEnergy_Interior(iLevel)       = Zero

      Momentum_X1(iLevel)           = Zero
      Momentum_X2(iLevel)           = Zero
      Momentum_X3(iLevel)           = Zero

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        iZ_B0(1) = iE_B0
        iZ_E0(1) = iE_E0


        iZ_B0(2:4) = iX_B0
        iZ_E0(2:4) = iX_E0

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

        CALL AllocateArray_Z &
               ( [ 1       , &
                   iZ_B0(1), &
                   iZ_B0(2), &
                   iZ_B0(3), &
                   iZ_B0(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E0(1), &
                   iZ_E0(2), &
                   iZ_E0(3), &
                   iZ_E0(4), &
                   nCR     , &
                   nSpecies ], &
                 U )

        CALL AllocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCF ], &
                 UF )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uCF, UF )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B0, iZ_E0, iLo_MF, iZ_B0, iZ_E0, uCR, U )

        CALL ComputeTally_TwoMoment( iZ_B0, iZ_E0, G, UF, U, iLevel )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nCF ], &
                 UF )

        CALL DeallocateArray_Z &
               ( [ 1       , &
                   iZ_B0(1), &
                   iZ_B0(2), &
                   iZ_B0(3), &
                   iZ_B0(4), &
                   1       , &
                   1        ], &
                 [ nDOFZ   , &
                   iZ_E0(1), &
                   iZ_E0(2), &
                   iZ_E0(3), &
                   iZ_E0(4), &
                   nCR     , &
                   nSpecies ], &
                 U )

        CALL DeallocateArray_X &
               ( [ 1    , iX_B0(1), iX_B0(2), iX_B0(3), 1   ], &
                 [ nDOFX, iX_E0(1), iX_E0(2), iX_E0(3), nGF ], &
                 G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO




    CALL amrex_parallel_reduce_sum( NeutrinoLeptonNumber_Interior, nLevels   )
    CALL amrex_parallel_reduce_sum( NeutrinoEnergy_Interior, nLevels         )

    CALL amrex_parallel_reduce_sum( Momentum_X1, nLevels         )
    CALL amrex_parallel_reduce_sum( Momentum_X2, nLevels         )
    CALL amrex_parallel_reduce_sum( Momentum_X3, nLevels         )



    IF( SetInitialValues )THEN

      DO iLevel = 0, nLevels-1

        NeutrinoLeptonNumber_Initial(iLevel) = NeutrinoLeptonNumber_Interior(iLevel)
        NeutrinoEnergy_Initial      (iLevel) = NeutrinoEnergy_Interior      (iLevel)

      END DO

    END IF


    DO iLevel = 0, nLevels-1

      NeutrinoLeptonNumber_Change(iLevel) &
        = NeutrinoLeptonNumber_Interior(iLevel) &
            - ( NeutrinoLeptonNumber_Initial(iLevel) + NeutrinoLeptonNumber_OffGrid(iLevel) )

      NeutrinoEnergy_Change(iLevel) &
        = NeutrinoEnergy_Interior(iLevel) &
            - ( NeutrinoEnergy_Initial(iLevel)       + NeutrinoEnergy_OffGrid(iLevel)       )

    END DO











    CALL WriteTally_TwoMoment( Time )

    IF( Verbose )THEN

      CALL DisplayTally( Time )

    END IF




  END SUBROUTINE ComputeTally_TwoMoment_MF

  SUBROUTINE IncrementOffGridTally_TwoMoment_MF( dM )

    REAL(DP), INTENT(in) :: dM(1:,0:)

    INTEGER :: iLevel

    IF( SuppressTally ) RETURN

    DO iLevel = 0, nLevels-1

      NeutrinoLeptonNumber_OffGrid(iLevel) &
        = NeutrinoLeptonNumber_OffGrid(iLevel) + FourPi * dM(iCR_N,iLevel) / hc3

      NeutrinoEnergy_OffGrid &
        = NeutrinoEnergy_OffGrid(iLevel)       + FourPi * dM(nCR+iCR_N,iLevel ) / hc3


    END DO
  END SUBROUTINE IncrementOffGridTally_TwoMoment_MF





  SUBROUTINE ComputeTally_TwoMoment( iZ_B0, iZ_E0, G, UF, U, iLevel )


    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iLevel
    REAL(DP), INTENT(in) :: &
      G(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)
    REAL(DP), INTENT(inout) :: &
      UF(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)

    TYPE(MeshType) :: MeshE
    TYPE(MeshType) :: MeshX(3)
    INTEGER        :: iNodeZ, iNodeX, iNodeE, iZ1, iZ2, iZ3, iZ4, iDim, iS
    REAL(DP)       :: W, vsq, dX
    REAL(DP) :: &
      PF(1:nDOFX, &
        iZ_B0(2):iZ_E0(2), &
        iZ_B0(3):iZ_E0(3), &
        iZ_B0(4):iZ_E0(4), &
        1:nPF)
    REAL(DP) :: d4Z(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                    iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                    iZ_B0(4):iZ_E0(4))

    CALL CreateMesh_MF( iLevel, MeshX )

    CALL CreateMesh &
           ( MeshE, nE, nNodesE, swE, eL, eR, zoomOption = zoomE )

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_Relativistic &
               ( UF (iNodeX,iZ2,iZ3,iZ4,iCF_D ),        &
                 UF (iNodeX,iZ2,iZ3,iZ4,iCF_S1),        &
                 UF (iNodeX,iZ2,iZ3,iZ4,iCF_S2),        &
                 UF (iNodeX,iZ2,iZ3,iZ4,iCF_S3),        &
                 UF (iNodeX,iZ2,iZ3,iZ4,iCF_E ),        &
                 UF (iNodeX,iZ2,iZ3,iZ4,iCF_Ne),        &
                 PF (iNodeX,iZ2,iZ3,iZ4,iPF_D ),        &
                 PF (iNodeX,iZ2,iZ3,iZ4,iPF_V1),        &
                 PF (iNodeX,iZ2,iZ3,iZ4,iPF_V2),        &
                 PF (iNodeX,iZ2,iZ3,iZ4,iPF_V3),        &
                 PF (iNodeX,iZ2,iZ3,iZ4,iPF_E ),        &
                 PF (iNodeX,iZ2,iZ3,iZ4,iPF_Ne),        &
                 G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)


      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE


        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        d4Z(iNodeZ,iZ1,iZ2,iZ3,iZ4)                             &
            =   FourPi * dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
            * Weights_q(iNodeZ)                                &
            * ( uGE(iNodeE,iZ1,iGE_Ep2) / hc3 )                &
            * G(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

      END DO
      END DO




    END DO
    END DO
    END DO
    END DO


    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)


      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        NeutrinoLeptonNumber_Interior(iLevel)             &
          = NeutrinoLeptonNumber_Interior(iLevel)        &
              + d4Z(iNodeZ,iZ1,iZ2,iZ3,iZ4)               &
                * LeptonNumber(iS)                        &
                * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS)

      END DO
      END DO




    END DO
    END DO
    END DO
    END DO
    END DO

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)


      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE


        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        d4Z(iNodeZ,iZ1,iZ2,iZ3,iZ4)                             &
            =   FourPi * dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
            * Weights_q(iNodeZ)                                &
            * ( uGE(iNodeE,iZ1,iGE_Ep3) / hc3 )                &
            * G(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)

      END DO
      END DO




    END DO
    END DO
    END DO
    END DO


    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)


      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE


        vsq = PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1)**2 * G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11) &
            + PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2)**2 * G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22) &
            + PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3)**2 * G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)
        W = 1.0_DP / SQRT( 1.0_DP - vsq )

        NeutrinoEnergy_Interior                                       &
          = NeutrinoEnergy_Interior                                   &
              + d4Z(iNodeZ,iZ1,iZ2,iZ3,iZ4)               &
                * ( W * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS)    &
                    + PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1)           &
                        * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) &
                    + PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2)           &
                        * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) &
                    + PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3)           &
                        * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) )


        Momentum_X1                                           &
          = Momentum_X1                                       &
              + d4Z(iNodeZ,iZ1,iZ2,iZ3,iZ4)               &
                * ( U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)       &
                    + W * G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)  &
                        * PF(iNodeX,iZ2,iZ3,iZ4,iPF_V1)       &
                        * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )


        Momentum_X2                                           &
          = Momentum_X2                                       &
              + d4Z(iNodeZ,iZ1,iZ2,iZ3,iZ4)               &
                * ( U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)       &
                    + W * G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)  &
                        * PF(iNodeX,iZ2,iZ3,iZ4,iPF_V2)       &
                        * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )


        Momentum_X3                                           &
          = Momentum_X3                                       &
              + d4Z(iNodeZ,iZ1,iZ2,iZ3,iZ4)               &
                * ( U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)       &
                    + W * G(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)  &
                        * PF(iNodeX,iZ2,iZ3,iZ4,iPF_V3)       &
                        * U(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )
      END DO
      END DO




    END DO
    END DO
    END DO
    END DO
    END DO

    END ASSOCIATE

    CALL DestroyMesh( MeshE )

    CALL DestroyMesh_MF( MeshX )

  END SUBROUTINE ComputeTally_TwoMoment


  SUBROUTINE WriteTally_Twomoment( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    IF( amrex_parallel_ioprocessor() )THEN

      ! --- Neutrino Lepton Number ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( NeutrinoLeptonNumber_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' )                  &
        Time / UnitsDisplay % TimeUnit, &
        NeutrinoLeptonNumber_Interior(0), &
        NeutrinoLeptonNumber_OffGrid (0), &
        NeutrinoLeptonNumber_Initial (0), &
        NeutrinoLeptonNumber_Change  (0)


      CLOSE( FileUnit )

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( NeutrinoEnergy_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' )                  &
        Time / UnitsDisplay % TimeUnit, &
        NeutrinoEnergy_Interior(0) / UnitsDisplay % EnergyGlobalUnit, &
        NeutrinoEnergy_OffGrid (0) / UnitsDisplay % EnergyGlobalUnit, &
        NeutrinoEnergy_Initial (0) / UnitsDisplay % EnergyGlobalUnit, &
        NeutrinoEnergy_Change  (0) / UnitsDisplay % EnergyGlobalUnit


      CLOSE( FileUnit )

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( Momentum_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' )                  &
        Time / UnitsDisplay % TimeUnit,                       &
        Momentum_X1(0),                                       &
        Momentum_X2(0),                                       &
        Momentum_X3(0)

      CLOSE( FileUnit )

    END IF

  END SUBROUTINE WriteTally_TwoMoment


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(A8,A,ES8.2E2,x,A)') &
        '', 'TwoMoment Tally. t = ', &
        Time / UnitsDisplay % TimeUnit, &
        UnitsDisplay % TimeLabel
      WRITE(*,*)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Neutrino Lepton Number.: ', &
        NeutrinoLeptonNumber_Interior(0)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Neutrino Energy.: ', &
        NeutrinoEnergy_Interior(0) / UnitsDisplay % EnergyGlobalUnit
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Neutrino Momentum1.: ', &
        Momentum_X1(0)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Neutrino Momentum2.: ', &
        Momentum_X2(0)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Neutrino Momentum3.: ', &
        Momentum_X3(0)


      WRITE(*,*)

    END IF

  END SUBROUTINE DisplayTally



END MODULE MF_TwoMoment_TallyModule
