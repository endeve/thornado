MODULE MF_AccretionShockUtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum, &
    amrex_parallel_reduce_min

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX, &
    nDimsX, &
    nNodesX
  USE UtilitiesModule, ONLY: &
    IsCornerCell
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1
  USE MeshModule, ONLY: &
    MeshType, &
    CreateMesh, &
    DestroyMesh
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
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
    iAF_P
  USE AccretionShockUtilitiesModule, ONLY: &
    ComputeAccretionShockDiagnostics, &
    ComputePowerInLegendreModes, &
    ComputeAngleAveragedShockRadius
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE UnitsModule, ONLY: &
    UnitsDisplay

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_X_Global
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    nX, &
    xL, &
    xR, &
    DEBUG, &
    UsePhysicalUnits
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeAccretionShockDiagnostics
  PUBLIC :: WriteNodal1DICToFile_SAS

  LOGICAL,          PUBLIC              :: WriteNodal1DIC_SAS
  LOGICAL,          PUBLIC              :: WriteAccretionShockDiagnostics
  CHARACTER(LEN=:), PUBLIC, ALLOCATABLE :: FileName_Nodal1DIC_SAS
  CHARACTER(LEN=:), PUBLIC, ALLOCATABLE :: AccretionShockDiagnosticsFileName


CONTAINS


  SUBROUTINE MF_ComputeAccretionShockDiagnostics( Time, MF_uGF, MF_uCF )

    REAL(DP),             INTENT(in) :: Time  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: A(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iLo_MF(4), &
               iNX, iNX1, iX1, iX2, iX3, iDim

    TYPE(MeshType) :: MeshX(3)

    INTEGER, PARAMETER    :: nLegModes = 3
    INTEGER               :: iLegMode
    REAL(DP)              :: Power(0:nLegModes-1)
    REAL(DP)              :: AngleAveragedShockRadius
    REAL(DP), ALLOCATABLE :: PowerIntegrand(:,:,:,:)
    REAL(DP), ALLOCATABLE :: ShockRadius   (:,:,:,:)

    INTEGER :: FileUnit
    LOGICAL :: IsFile

    CHARACTER(256) :: TimeLabel, PowerLabel(0:nLegModes-1), ShockRadiusLabel

    IF( .NOT. WriteAccretionShockDiagnostics ) RETURN

    IF( nDimsX .EQ. 1 ) RETURN

    IF( amrex_parallel_ioprocessor() )THEN

      INQUIRE( FILE = TRIM( AccretionShockDiagnosticsFileName ), &
               EXIST = IsFile )

      IF( .NOT. IsFile )THEN

        TimeLabel = 'Time [' // TRIM( UnitsDisplay % TimeLabel ) // ']'

        IF( UsePhysicalUnits )THEN

          PowerLabel(0) = 'P0 (Entropy) [cgs]'
          PowerLabel(1) = 'P1 (Entropy) [cgs]'
          PowerLabel(2) = 'P2 (Entropy) [cgs]'

        ELSE

          PowerLabel(0) = 'P0 (Entropy) []'
          PowerLabel(1) = 'P1 (Entropy) []'
          PowerLabel(2) = 'P2 (Entropy) []'

        END IF

        ShockRadiusLabel &
          = 'Shock Radius [' // TRIM( UnitsDisplay % LengthX1Label ) // ']'

        OPEN( FileUnit, FILE = TRIM( AccretionShockDiagnosticsFileName ) )

        WRITE( FileUnit, '(5(A25,1x))' ) &
          TRIM( TimeLabel ), &
          TRIM( PowerLabel(0) ), TRIM( PowerLabel(1) ), TRIM( PowerLabel(2) ), &
          TRIM( ShockRadiusLabel )

        CLOSE( FileUnit )

      END IF

    END IF

    DO iDim = 1, 3

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim) )

    END DO

    ALLOCATE( PowerIntegrand(0:nLevels-1,0:nLegModes, &
                             1:nNodesX(1),1:nX(1)) )

    ALLOCATE( ShockRadius(0:nLevels-1,nDOFX_X1,nX(2),nX(3)) )

    PowerIntegrand = Zero
    ShockRadius    = HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( G(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nCF) )

        ALLOCATE( P(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nPF) )

        ALLOCATE( A(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nAF) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_X( nGF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uGF, G )
        CALL amrex2thornado_X( nCF, iX_B0, iX_E0, iLo_MF, iX_B0, iX_E0, uCF, U )

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)

          CALL ComputePrimitive_Euler &
                 ( U   (:,iX1,iX2,iX3,iCF_D ), &
                   U   (:,iX1,iX2,iX3,iCF_S1), &
                   U   (:,iX1,iX2,iX3,iCF_S2), &
                   U   (:,iX1,iX2,iX3,iCF_S3), &
                   U   (:,iX1,iX2,iX3,iCF_E ), &
                   U   (:,iX1,iX2,iX3,iCF_Ne), &
                   P   (:,iX1,iX2,iX3,iPF_D ), &
                   P   (:,iX1,iX2,iX3,iPF_V1), &
                   P   (:,iX1,iX2,iX3,iPF_V2), &
                   P   (:,iX1,iX2,iX3,iPF_V3), &
                   P   (:,iX1,iX2,iX3,iPF_E ), &
                   P   (:,iX1,iX2,iX3,iPF_Ne), &
                   G   (:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   G   (:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   G   (:,iX1,iX2,iX3,iGF_Gm_dd_33) )

          CALL ComputePressureFromPrimitive &
                 ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_E), &
                   P(:,iX1,iX2,iX3,iPF_Ne), A(:,iX1,iX2,iX3,iAF_P) )

        END DO
        END DO
        END DO

        CALL ComputeAccretionShockDiagnostics &
               ( iX_B0, iX_E0, P, A, MeshX, &
                 PowerIntegrand(iLevel,:,:,iX_B0(1):iX_E0(1)), &
                 ShockRadius(iLevel,:,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)) )

        DEALLOCATE( A )
        DEALLOCATE( P )
        DEALLOCATE( U )
        DEALLOCATE( G )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iX1 = 1, nX(1)
    DO iNX1 = 1, nNodesX(1)
    DO iLegMode = 0, nLegModes-1

      CALL amrex_parallel_reduce_sum &
             ( PowerIntegrand(:,iLegMode,iNX1,iX1), nLevels )

    END DO
    END DO
    END DO

    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iNX1 = 1, nDOFX_X1

      CALL amrex_parallel_reduce_min( ShockRadius(:,iNX1,iX2,iX3), nLevels )

    END DO
    END DO
    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      CALL ComputePowerInLegendreModes &
             ( 1, nX(1), MeshX(1), PowerIntegrand(0,:,:,:), Power )

      CALL ComputeAngleAveragedShockRadius &
             ( [ 1, 1, 1 ], nX, ShockRadius(0,:,:,:), AngleAveragedShockRadius )

      OPEN( FileUnit, FILE = TRIM( AccretionShockDiagnosticsFileName ), &
            POSITION = 'APPEND' )

      WRITE( FileUnit, '(SPES25.16E3,1x)', ADVANCE = 'NO' ) &
        Time(0) / UnitsDisplay % TimeUnit

      WRITE( FileUnit, '(3(SPES25.16E3,1x))', ADVANCE = 'NO' ) &
        Power

      WRITE( FileUnit, '(SPES25.16E3,1x)', ADVANCE = 'NO' ) &
        AngleAveragedShockRadius / UnitsDisplay % LengthX1Unit

      WRITE( FileUnit, * )

      CLOSE( FileUnit )

    END IF

    DEALLOCATE( ShockRadius    )
    DEALLOCATE( PowerIntegrand )

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

  END SUBROUTINE MF_ComputeAccretionShockDiagnostics


  SUBROUTINE WriteNodal1DICToFile_SAS( GEOM, MF_uGF, MF_uCF )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)

    INTEGER           :: iLo(3), iHi(3), iNX, iX1, iX2, iX3
    CHARACTER(LEN=16) :: FMT

    REAL(DP) :: P(1:nDOFX,1:nPF)
    REAL(DP) :: A(1:nDOFX,1:nAF)
    REAL(DP) :: G(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                          1-swX(2):nX(2)+swX(2), &
                          1-swX(3):nX(3)+swX(3), &
                  1:nGF)
    REAL(DP) :: U(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                          1-swX(2):nX(2)+swX(2), &
                          1-swX(3):nX(3)+swX(3), &
                  1:nCF)

    IF( .NOT. WriteNodal1DIC_SAS ) RETURN

    IF( .NOT. nDimsX .EQ. 1 ) RETURN

    CALL amrex2thornado_X_Global &
           ( GEOM, MF_uGF, nGF, G, ApplyBC_Option = .FALSE. )

    CALL amrex2thornado_X_Global &
           ( GEOM, MF_uCF, nCF, U, ApplyBC_Option = .TRUE. )

    IF( amrex_parallel_ioprocessor() )THEN

      iLo = 1  - swX
      iHi = nX + swX

      OPEN( UNIT = 101, FILE = TRIM( FileName_Nodal1DIC_SAS ) )

      WRITE(FMT,'(A3,I3.3,A10)') '(SP', nDOFX, 'ES25.16E3)'

      WRITE(101,'(A)') FMT

      DO iX3 = iLo(3), iHi(3)
      DO iX2 = iLo(2), iHi(2)
      DO iX1 = iLo(1), iHi(1)

        IF( IsCornerCell( iLo, iHi, iX1, iX2, iX3 ) ) CYCLE

        CALL ComputePrimitive_Euler &
               ( U   (:,iX1,iX2,iX3,iCF_D ), &
                 U   (:,iX1,iX2,iX3,iCF_S1), &
                 U   (:,iX1,iX2,iX3,iCF_S2), &
                 U   (:,iX1,iX2,iX3,iCF_S3), &
                 U   (:,iX1,iX2,iX3,iCF_E ), &
                 U   (:,iX1,iX2,iX3,iCF_Ne), &
                 P   (:,iPF_D ), &
                 P   (:,iPF_V1), &
                 P   (:,iPF_V2), &
                 P   (:,iPF_V3), &
                 P   (:,iPF_E ), &
                 P   (:,iPF_Ne), &
                 G   (:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G   (:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G   (:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( P(:,iPF_D ), P(:,iPF_E ), P(:,iPF_Ne), A(:,iAF_P) )

        WRITE(101,FMT) P(:,iPF_D )
        WRITE(101,FMT) P(:,iPF_V1)
        WRITE(101,FMT) A(:,iAF_P )

      END DO
      END DO
      END DO

      CLOSE( 101 )

    END IF

  END SUBROUTINE WriteNodal1DICToFile_SAS


END MODULE MF_AccretionShockUtilitiesModule
