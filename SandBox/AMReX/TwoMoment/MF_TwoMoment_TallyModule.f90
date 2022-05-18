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
  USE MyAmrModule, ONLY: &
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
    amrex2thornado_Z




  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_InitializeTally_TwoMoment
  PUBLIC :: MF_ComputeTally_TwoMoment
  PUBLIC :: MF_FinalizeTally_TwoMoment

  LOGICAL :: SuppressTally


  REAL(DP) :: hc3    


  CHARACTER(256) :: Neutrino_FileName
  REAL(DP), ALLOCATABLE :: Neutrino_LeptonNumber(:)
  REAL(DP), ALLOCATABLE :: Neutrino_Energy(:)

  
  CHARACTER(256) :: Momentum_FileName
  REAL(DP), ALLOCATABLE :: Momentum_X1(:)
  REAL(DP), ALLOCATABLE :: Momentum_X2(:)
  REAL(DP), ALLOCATABLE :: Momentum_X3(:)



CONTAINS


  SUBROUTINE MF_InitializeTally_TwoMoment &
    ( SuppressTally_Option, BaseFileName_Option )

    LOGICAL,  INTENT(in),         OPTIONAL :: &
      SuppressTally_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      BaseFileName_Option


    CHARACTER(256) :: BaseFileName
    INTEGER        :: FileUnit

    CHARACTER(256) :: TimeLabel
    CHARACTER(256) :: LeptonNumberLabel
    CHARACTER(256) :: EnergyLabel
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


    ALLOCATE( Neutrino_LeptonNumber(0:nLevels-1) )
    ALLOCATE( Neutrino_Energy(0:nLevels-1) )
    ALLOCATE( Momentum_X1(0:nLevels-1) )
    ALLOCATE( Momentum_X2(0:nLevels-1) )
    ALLOCATE( Momentum_X3(0:nLevels-1) )


    IF( amrex_parallel_ioprocessor() )THEN

      BaseFileName = ''
      IF( PRESENT( BaseFileName_Option ) ) &
        BaseFileName = TRIM( BaseFileName_Option )

      BaseFileName = TRIM( BaseFileName ) // TRIM( ProgramName )

      ! --- Neutrino ---

      Neutrino_FileName &
        = TRIM( BaseFileName ) // '_Tally_Neutrino.dat'

      TimeLabel     &
        = 'Time ['     // TRIM( UnitsDisplay % TimeLabel ) // ']'
      LeptonNumberLabel &
        = 'LeptonNumber'
      EnergyLabel &
        = 'Energy [' // TRIM( UnitsDisplay % EnergyGlobalLabel ) // ']'



      OPEN( NEWUNIT = FileUnit, FILE = TRIM(Neutrino_FileName ) )

      WRITE(FileUnit,'(5(A25,x))') &
        TRIM( TimeLabel ), TRIM( LeptonNumberLabel ), TRIM( EnergyLabel )
      CLOSE( FileUnit )

      ! --- Momentum ---

      Momentum_FileName &
        = TRIM( BaseFileName ) // '_Tally_Momentum.dat'
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

    Neutrino_LeptonNumber = Zero
    Neutrino_Energy       = Zero

    Momentum_X1           = Zero
    Momentum_X2           = Zero
    Momentum_X3           = Zero

  END SUBROUTINE MF_InitializeTally_TwoMoment



  SUBROUTINE MF_FinalizeTally_TwoMoment

    IF( SuppressTally ) RETURN

    DEALLOCATE( Neutrino_LeptonNumber )
    DEALLOCATE( Neutrino_Energy )
    DEALLOCATE( Momentum_X1 )
    DEALLOCATE( Momentum_X2 )
    DEALLOCATE( Momentum_X3 )

  END SUBROUTINE MF_FinalizeTally_TwoMoment


  SUBROUTINE MF_ComputeTally_TwoMoment &
    ( GEOM, MF_uGF, MF_uCF, MF_uCR, Time, Verbose_Option )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCR(0:nLevels-1)


    REAL(DP),             INTENT(in) :: Time
    LOGICAL,              INTENT(in), OPTIONAL :: Verbose_Option



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


    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option



    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      Neutrino_LeptonNumber(iLevel) = Zero
      Neutrino_Energy(iLevel)       = Zero

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
 
        
  

        ALLOCATE( G(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nGF) )


        ALLOCATE( U (1:nDOFZ,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                             iZ_B0(3):iZ_E0(3), &
                             iZ_B0(4):iZ_E0(4),1:nCR,1:nSpecies) )

        ALLOCATE( UF (1:nDOFX,iZ_B0(2):iZ_E0(2), &
                             iZ_B0(3):iZ_E0(3), &
                             iZ_B0(4):iZ_E0(4),1:nCF) )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uCF, UF )

        CALL amrex2thornado_Z &
               ( nCR, nSpecies, nE, iE_B0, iE_E0, &
                 iZ_B0, iZ_E0, iLo_MF, iZ_B0, iZ_E0, uCR, U )


        CALL ComputeTally_TwoMoment( iZ_B0, iZ_E0, G, UF, U, iLevel )



        DEALLOCATE( U )
        DEALLOCATE( UF)
        DEALLOCATE( G )


      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    CALL WriteTally_TwoMoment( Time )

    IF( Verbose )THEN

      CALL DisplayTally( Time )

    END IF




  END SUBROUTINE MF_ComputeTally_TwoMoment


  SUBROUTINE ComputeTally_TwoMoment( iZ_B0, iZ_E0, G, UF, U, iLevel )


    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iLevel
    REAL(DP), INTENT(in) :: &
      G(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)
    REAL(DP), INTENT(in) :: &
      U(1:,iZ_B0(1):,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:,1:)
    REAL(DP), INTENT(in) :: &
      UF(1:,iZ_B0(2):,iZ_B0(3):,iZ_B0(4):,1:)

    TYPE(MeshType) :: MeshE
    TYPE(MeshType) :: MeshX(3)
    INTEGER        :: iNodeZ, iNodeX, iNodeE, iZ1, iZ2, iZ3, iZ4, iDim, iS
    REAL(DP)       :: W, vsq
    REAL(DP) :: &
      PF(1:nDOFX, &
        iZ_B0(2):iZ_E0(2), &
        iZ_B0(3):iZ_E0(3), &
        iZ_B0(4):iZ_E0(4), &
        1:nPF)
    REAL(DP) :: d4Z(1:nDOFZ,iZ_B0(1):iZ_E0(1), &
                    iZ_B0(2):iZ_E0(2),iZ_B0(3):iZ_E0(3), &
                    iZ_B0(4):iZ_E0(4))


    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), 0, &
               xL(iDim), xR(iDim) )

    END DO




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
        
        Neutrino_LeptonNumber(iLevel)                     &
          = Neutrino_LeptonNumber(iLevel)                 &
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

        Neutrino_Energy                                       &
          = Neutrino_Energy                                   &
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




    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

    CALL DestroyMesh( MeshE )

    END ASSOCIATE

  END SUBROUTINE ComputeTally_TwoMoment


  SUBROUTINE WriteTally_Twomoment( Time )

    REAL(DP), INTENT(in) :: Time

    INTEGER :: FileUnit

    IF( amrex_parallel_ioprocessor() )THEN

      ! --- Neutrino ---

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( Neutrino_FileName ), &
            POSITION = 'APPEND', ACTION = 'WRITE' )

      WRITE( FileUnit, '(5(ES25.16E3,1x))' )                  &
        Time / UnitsDisplay % TimeUnit,                       &
        Neutrino_LeptonNumber(0),                             & 
        Neutrino_Energy(0) / UnitsDisplay % EnergyGlobalUnit

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
        Neutrino_LeptonNumber(0)
      WRITE(*,'(A6,A40,ES14.7E2,x,A)') &
        '', 'Neutrino Energy.: ', &
        Neutrino_Energy(0) / UnitsDisplay % EnergyGlobalUnit
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
