!> Module for operations on MultiFabs
MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,                 ONLY: &
    AR => amrex_real
  USE amrex_box_module,                  ONLY: &
    amrex_box
  USE amrex_parallel_module,             ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum
  USE amrex_geometry_module,             ONLY: &
    amrex_geometry
  USE amrex_multifab_module,             ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy

  ! --- thornado Modules ---

  USE ProgramHeaderModule,               ONLY: &
    nDOFX, &
    nNodesX
  USE MeshModule,                        ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule,              ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha,    &
    iGF_Psi,      &
    iGF_SqrtGm,   &
    unitsGF
  USE FluidFieldsModule,                 ONLY: &
    nCF,     &
    iCF_D,   &
    iCF_S1,  &
    iCF_S2,  &
    iCF_S3,  &
    iCF_E,   &
    iCF_Ne,  &
    nPF,     &
    iPF_D,   &
    iPF_V1,  &
    iPF_V2,  &
    iPF_V3,  &
    iPF_E,   &
    iPF_Ne,  &
    unitsPF, &
    nDF,     &
    iDF_TCI, &
    unitsAF, &
    nAF,     &
    iAF_P
  USE Euler_UtilitiesModule,             ONLY: &
    ComputePrimitive_Euler
  USE EquationOfStateModule,             ONLY: &
    ComputePressureFromPrimitive

  ! --- Local Modules ---

  USE InputParsingModule,                ONLY: &
    nLevels, &
    nX,      &
    swX
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap, &
    MF_ApplyBoundaryConditions_Euler
  USE TimersModule_AMReX_Euler,           ONLY: &
    TimersStart_AMReX_Euler,      &
    TimersStop_AMReX_Euler,       &
    Timer_AMReX_Euler_DataTransfer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: amrex2thornado_Euler
  PUBLIC :: thornado2amrex_Euler
  PUBLIC :: ShowVariableFromMultiFab
  PUBLIC :: WriteNodalDataToFile_SAS
  PUBLIC :: WriteNodalDataToFile
  PUBLIC :: CombineGridData

  REAL(AR), PARAMETER, PUBLIC :: Zero = 0.0_AR


CONTAINS


  SUBROUTINE amrex2thornado_Euler &
    ( nFields, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B(3), iX_E(3)
    REAL(AR), INTENT(in)  :: &
      Data_amrex   (   iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(AR), INTENT(out) :: &
      Data_thornado(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)

    INTEGER :: iX1, iX2, iX3, iField

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      DO iField = 1, nFields

        Data_thornado(1:nDOFX,iX1,iX2,iX3,iField) &
          = Data_amrex(iX1,iX2,iX3,nDOFX*(iField-1)+1:nDOFX*iField)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE amrex2thornado_Euler


  SUBROUTINE thornado2amrex_Euler &
    ( nFields, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nFields
    INTEGER,  INTENT(in)  :: iX_B(3), iX_E(3)
    REAL(AR), INTENT(out) :: &
      Data_amrex   (   iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(AR), INTENT(in)  :: &
      Data_thornado(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)

    INTEGER :: iX1, iX2, iX3, iField

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      DO iField = 1, nFields

        Data_amrex(iX1,iX2,iX3,nDOFX*(iField-1)+1:nDOFX*iField) &
          = Data_thornado(1:nDOFX,iX1,iX2,iX3,iField)

      END DO

    END DO
    END DO
    END DO

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

  END SUBROUTINE thornado2amrex_Euler


  SUBROUTINE ShowVariableFromMultiFab( MF, swXX, iComp )

    TYPE(amrex_multifab), INTENT(in) :: MF(0:nLevels-1)
    INTEGER,              INTENT(in) :: swXX(3)
    INTEGER,              INTENT(in) :: iComp

    INTEGER                       :: iX1, iX2, iX3, iLevel
    INTEGER                       :: lo(4), hi(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: U(:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        U => MF(iLevel) % DataPtr( MFI )
        BX = MFI % tilebox()

        lo = LBOUND( U ); hi = UBOUND( U )

        DO iX3 = BX % lo(3) - swXX(3), BX % hi(3) + swXX(3)
        DO iX2 = BX % lo(2) - swXX(2), BX % hi(2) + swXX(2)
        DO iX1 = BX % lo(1) - swXX(1), BX % hi(1) + swXX(1)

          WRITE(*,'(A,3I4.3,ES10.1E3)') &
            'iX1, iX2, iX3, Data: ',iX1, iX2, iX3, U(iX1,iX2,iX3,iComp)

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    WRITE(*,*)

  END SUBROUTINE ShowVariableFromMultiFab


  SUBROUTINE WriteNodalDataToFile_SAS( GEOM, MF_uGF, MF_uCF, FileNameBase )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    CHARACTER(LEN=*)    , INTENT(in) :: FileNameBase

    INTEGER            :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3),   &
                          iX_B (3), iX_E (3),                       &
                          iX1, iX2, iX3, iNodeX, iCF, iGF,          &
                          iLevel, nCompGF, nCompCF
    CHARACTER(LEN=16)  :: FMT
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(EdgeMap)      :: Edge_Map

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR), ALLOCATABLE         :: G  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: U  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: P  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: A  (:,:,:,:,:)

    ALLOCATE( G(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nGF) )

    ALLOCATE( U(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nCF) )

    ALLOCATE( P(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nPF) )

    ALLOCATE( A(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nAF) )

    G = 0.0_AR
    U = 0.0_AR
    P = 0.0_AR
    A = 0.0_AR

    ! --- Convert AMReX MultiFabs to thornado arrays ---

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF     => MF_uGF(iLevel) % DataPtr( MFI )
        nCompGF =  MF_uGF(iLevel) % nComp()

        uCF     => MF_uCF(iLevel) % DataPtr( MFI )
        nCompCF =  MF_uCF(iLevel) % nComp()

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iX_B = iX_B0
        iX_E = iX_E0

        ! --- Ensure exchange cells are excluded ---

        IF( iX_B0(1) .EQ. 1     ) iX_B(1) = 1     - swX(1)
        IF( iX_B0(2) .EQ. 1     ) iX_B(2) = 1     - swX(2)
        IF( iX_B0(3) .EQ. 1     ) iX_B(3) = 1     - swX(3)
        IF( iX_E0(1) .EQ. nX(1) ) iX_E(1) = nX(1) + swX(1)
        IF( iX_E0(2) .EQ. nX(2) ) iX_E(2) = nX(2) + swX(2)
        IF( iX_E0(3) .EQ. nX(3) ) iX_E(3) = nX(3) + swX(3)

        CALL amrex2thornado_Euler( nGF, iX_B, iX_E, &
                                   uGF(iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nCompGF), &
                                   G  (1:nDOFX, &
                                       iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nGF) )

        CALL amrex2thornado_Euler( nCF, iX_B, iX_E, &
                                   uCF(iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nCompCF), &
                                   U  (1:nDOFX, &
                                       iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nCF) )

        ! --- Apply boundary conditions ---

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

        CALL MF_ApplyBoundaryConditions_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                  U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF), Edge_Map )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    ! --- Combine data from all processes ---

    DO iX3    = 1-swX(3), nX(3)+swX(3)
    DO iX2    = 1-swX(2), nX(2)+swX(2)
    DO iX1    = 1-swX(1), nX(1)+swX(1)
    DO iNodeX = 1       , nDOFX

      DO iGF = 1, nGF

        CALL amrex_parallel_reduce_sum( G(iNodeX,iX1,iX2,iX3,iGF) )

      END DO

      DO iCF = 1, nCF

        CALL amrex_parallel_reduce_sum( U(iNodeX,iX1,iX2,iX3,iCF) )

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Compute primitive and write to file ---

    IF( amrex_parallel_ioprocessor() )THEN

      DO iX3 = 1-swX(3), nX(3)+swX(3)
      DO iX2 = 1-swX(2), nX(2)+swX(2)
      DO iX1 = 1-swX(1), nX(1)+swX(1)

        CALL ComputePrimitive_Euler &
               ( U(:,iX1,iX2,iX3,iCF_D ),       &
                 U(:,iX1,iX2,iX3,iCF_S1),       &
                 U(:,iX1,iX2,iX3,iCF_S2),       &
                 U(:,iX1,iX2,iX3,iCF_S3),       &
                 U(:,iX1,iX2,iX3,iCF_E ),       &
                 U(:,iX1,iX2,iX3,iCF_Ne),       &
                 P(:,iX1,iX2,iX3,iPF_D ),       &
                 P(:,iX1,iX2,iX3,iPF_V1),       &
                 P(:,iX1,iX2,iX3,iPF_V2),       &
                 P(:,iX1,iX2,iX3,iPF_V3),       &
                 P(:,iX1,iX2,iX3,iPF_E ),       &
                 P(:,iX1,iX2,iX3,iPF_Ne),       &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_E ), &
                 P(:,iX1,iX2,iX3,iPF_Ne), A(:,iX1,iX2,iX3,iAF_P) )

      END DO
      END DO
      END DO

      OPEN( UNIT = 101, FILE = TRIM( FileNameBase ) // '_D.dat' )
      OPEN( UNIT = 102, FILE = TRIM( FileNameBase ) // '_V.dat' )
      OPEN( UNIT = 103, FILE = TRIM( FileNameBase ) // '_P.dat' )

      WRITE(FMT,'(A3,I3.3,A10)') '(SP', nDOFX, 'ES25.16E3)'

      WRITE(101,'(A)') FMT
      WRITE(102,'(A)') FMT
      WRITE(103,'(A)') FMT

      DO iX3 = 1-swX(3), nX(3)+swX(3)
      DO iX2 = 1-swX(2), nX(2)+swX(2)
      DO iX1 = 1-swX(1), nX(1)+swX(1)

        WRITE(101,FMT) P(:,iX1,iX2,iX3,iPF_D )
        WRITE(102,FMT) P(:,iX1,iX2,iX3,iPF_V1)
        WRITE(103,FMT) A(:,iX1,iX2,iX3,iAF_P )

      END DO
      END DO
      END DO

      CLOSE( 103 )
      CLOSE( 102 )
      CLOSE( 101 )

    END IF

    DEALLOCATE( A )
    DEALLOCATE( P )
    DEALLOCATE( U )
    DEALLOCATE( G )

  END SUBROUTINE WriteNodalDataToFile_SAS


  SUBROUTINE WriteNodalDataToFile( GEOM, MF_uGF, MF_uCF, MF_uDF, FileNameBase )

    TYPE(amrex_geometry), INTENT(in) :: GEOM  (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uDF(0:nLevels-1)
    CHARACTER(LEN=*)    , INTENT(in) :: FileNameBase

    INTEGER            :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                          iX_B (3), iX_E (3), iXL(3), iXR(3),     &
                          iLevel, nCompGF, nCompCF, nCompDF,      &
                          iX1, iX2, iX3, iCF, iGF,                &
                          iNodeX, iNodeX1, iNodeX2, iNodeX3
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    TYPE(EdgeMap)      :: Edge_Map

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uDF(:,:,:,:)
    REAL(AR), ALLOCATABLE         :: G  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: U  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: P  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: A  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: D  (:,:,:,:,:)

    REAL(AR) :: Alpha, Psi, SqrtGm, rho, v, e, pr

    ALLOCATE( G(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nGF) )

    ALLOCATE( U(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nCF) )

    ALLOCATE( P(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nPF) )

    ALLOCATE( A(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nAF) )

    ALLOCATE( D(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nDF) )

    G = 0.0_AR
    U = 0.0_AR
    P = 0.0_AR
    D = 0.0_AR

    ! --- Convert AMReX MultiFabs to thornado arrays ---

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uDF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF     => MF_uGF(iLevel) % DataPtr( MFI )
        nCompGF =  MF_uGF(iLevel) % nComp()

        uCF     => MF_uCF(iLevel) % DataPtr( MFI )
        nCompCF =  MF_uCF(iLevel) % nComp()

        uDF     => MF_uDF(iLevel) % DataPtr( MFI )
        nCompDF =  MF_uDF(iLevel) % nComp()

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        iX_B = iX_B0
        iX_E = iX_E0

        ! --- Ensure exchange cells are excluded ---

        IF( iX_B0(1) .EQ. 1     ) iX_B(1) = 1     - swX(1)
        IF( iX_B0(2) .EQ. 1     ) iX_B(2) = 1     - swX(2)
        IF( iX_B0(3) .EQ. 1     ) iX_B(3) = 1     - swX(3)
        IF( iX_E0(1) .EQ. nX(1) ) iX_E(1) = nX(1) + swX(1)
        IF( iX_E0(2) .EQ. nX(2) ) iX_E(2) = nX(2) + swX(2)
        IF( iX_E0(3) .EQ. nX(3) ) iX_E(3) = nX(3) + swX(3)

        CALL amrex2thornado_Euler( nGF, iX_B, iX_E, &
                                   uGF(iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nCompGF), &
                                   G  (1:nDOFX, &
                                       iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nGF) )

        CALL amrex2thornado_Euler( nCF, iX_B, iX_E, &
                                   uCF(iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nCompCF), &
                                   U  (1:nDOFX, &
                                       iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nCF) )

        CALL amrex2thornado_Euler( nDF, iX_B, iX_E, &
                                   uDF(iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nCompDF), &
                                   D  (1:nDOFX, &
                                       iX_B(1):iX_E(1), &
                                       iX_B(2):iX_E(2), &
                                       iX_B(3):iX_E(3),1:nDF) )

        ! --- Apply boundary conditions ---

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

        CALL MF_ApplyBoundaryConditions_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                  U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF), Edge_Map )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    ! --- Combine data from all processes ---

    DO iX3    = 1-swX(3), nX(3)+swX(3)
    DO iX2    = 1-swX(2), nX(2)+swX(2)
    DO iX1    = 1-swX(1), nX(1)+swX(1)
    DO iNodeX = 1       , nDOFX

      DO iGF = 1, nGF

        CALL amrex_parallel_reduce_sum( G(iNodeX,iX1,iX2,iX3,iGF) )

      END DO

      DO iCF = 1, nCF

        CALL amrex_parallel_reduce_sum( U(iNodeX,iX1,iX2,iX3,iCF) )

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Compute primitive and write to file ---

    IF( amrex_parallel_ioprocessor() )THEN

      DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
      DO iX1 = 1, nX(1)

        CALL ComputePrimitive_Euler &
               ( U(:,iX1,iX2,iX3,iCF_D ),       &
                 U(:,iX1,iX2,iX3,iCF_S1),       &
                 U(:,iX1,iX2,iX3,iCF_S2),       &
                 U(:,iX1,iX2,iX3,iCF_S3),       &
                 U(:,iX1,iX2,iX3,iCF_E ),       &
                 U(:,iX1,iX2,iX3,iCF_Ne),       &
                 P(:,iX1,iX2,iX3,iPF_D ),       &
                 P(:,iX1,iX2,iX3,iPF_V1),       &
                 P(:,iX1,iX2,iX3,iPF_V2),       &
                 P(:,iX1,iX2,iX3,iPF_V3),       &
                 P(:,iX1,iX2,iX3,iPF_E ),       &
                 P(:,iX1,iX2,iX3,iPF_Ne),       &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_E ), &
                 P(:,iX1,iX2,iX3,iPF_Ne), A(:,iX1,iX2,iX3,iAF_P) )

      END DO
      END DO
      END DO

      OPEN( UNIT = 100, FILE = TRIM( FileNameBase ) // 'r.dat'      )
      OPEN( UNIT = 101, FILE = TRIM( FileNameBase ) // 'Alpha.dat'  )
      OPEN( UNIT = 102, FILE = TRIM( FileNameBase ) // 'Psi.dat'    )
      OPEN( UNIT = 103, FILE = TRIM( FileNameBase ) // 'SqrtGm.dat' )
      OPEN( UNIT = 104, FILE = TRIM( FileNameBase ) // 'D.dat'      )
      OPEN( UNIT = 105, FILE = TRIM( FileNameBase ) // 'V.dat'      )
      OPEN( UNIT = 106, FILE = TRIM( FileNameBase ) // 'E.dat'      )
      OPEN( UNIT = 107, FILE = TRIM( FileNameBase ) // 'P.dat'      )

      ! --- Hacked to work only for 1D and 2D problems ---

      DO iX3     = 1, nX(3)
      DO iNodeX3 = 1, nNodesX(3)

        DO iX2     = 1, nX(2)
        DO iNodeX2 = 1, nNodesX(2)

          IF     ( iNodeX2 .EQ. 1 )THEN

            iXL(1) = 1
            iXR(1) = nNodesX(1)

          ELSE IF( iNodeX2 .EQ. 2 )THEN

            iXL(1) = nNodesX(1) + 1
            iXR(1) = 2 * nNodesX(1)

          ELSE IF( iNodeX2 .EQ. 3 )THEN

            iXL(1) = 2 * nNodesX(1) + 1
            iXR(1) = 3 * nNodesX(1)

          END IF

          DO iX1     = 1, nX(1)
          DO iNodeX1 = iXL(1), iXR(1)

            Alpha  = G(iNodeX1,iX1,iX2,iX3,iGF_Alpha  )
            Psi    = G(iNodeX1,iX1,iX2,iX3,iGF_Psi    )
            SqrtGm = G(iNodeX1,iX1,iX2,iX3,iGF_SqrtGm )
            rho    = P(iNodeX1,iX1,iX2,iX3,iPF_D      )
            v      = P(iNodeX1,iX1,iX2,iX3,iPF_V1     )
            e      = P(iNodeX1,iX1,iX2,iX3,iPF_E      )
            pr     = A(iNodeX1,iX1,iX2,iX3,iAF_P      )

            WRITE(100,'(ES24.16E3,1x)',ADVANCE='NO') &
              NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            WRITE(101,'(ES24.16E3,1x)',ADVANCE='NO') &
              Alpha &
                / unitsGF(iGF_Alpha)

            WRITE(102,'(ES24.16E3,1x)',ADVANCE='NO') &
              Psi &
                / unitsGF(iGF_Psi)

            WRITE(103,'(ES24.16E3,1x)',ADVANCE='NO') &
              SqrtGm &
                / unitsGF(iGF_SqrtGm)

            WRITE(104,'(ES24.16E3,1x)',ADVANCE='NO') &
              rho &
                / unitsPF(iPF_D)

            WRITE(105,'(ES24.16E3,1x)',ADVANCE='NO') &
              v &
                / unitsPF(iPF_V1)

            WRITE(106,'(ES24.16E3,1x)',ADVANCE='NO') &
              e &
                / unitsPF(iPF_E)

            WRITE(107,'(ES24.16E3,1x)',ADVANCE='NO') &
              pr &
                / unitsAF(iAF_P)

          END DO
          END DO

        WRITE(100,*)
        WRITE(101,*)
        WRITE(102,*)
        WRITE(103,*)
        WRITE(104,*)
        WRITE(105,*)
        WRITE(106,*)
        WRITE(107,*)

        END DO
        END DO

      END DO
      END DO

      CLOSE( 107 )
      CLOSE( 106 )
      CLOSE( 105 )
      CLOSE( 104 )
      CLOSE( 103 )
      CLOSE( 102 )
      CLOSE( 101 )
      CLOSE( 100 )

    END IF

    DEALLOCATE( D )
    DEALLOCATE( P )
    DEALLOCATE( U )
    DEALLOCATE( G )

  END SUBROUTINE WriteNodalDataToFile


  SUBROUTINE CombineGridData( MF, nF, iField, U )

    TYPE(amrex_multifab), INTENT(in)  :: MF(0:nLevels-1)
    INTEGER,              INTENT(in)  :: nF, iField
    REAL(AR),             INTENT(out) :: U(1:,1-swX(1):,1-swX(2):,1-swX(3):)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(AR), CONTIGUOUS, POINTER :: U_P(:,:,:,:)
    REAL(AR)                      :: U_K(nDOFX,nF)

    REAL(AR) :: U0(nDOFX,1-swX(1):nX(1)+swX(1), &
                         1-swX(2):nX(2)+swX(2), &
                         1-swX(3):nX(3)+swX(3),0:nLevels-1)

    INTEGER  :: iNodeX, iX1, iX2, iX3, iLevel, iX_B(3), iX_E(3), lo(4), hi(4)

    U0 = Zero

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        U_P => MF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        lo = LBOUND( U_P )
        hi = UBOUND( U_P )

        iX_B = BX % lo
        iX_E = BX % hi

        ! -- Get physical ghost cells right ---

        IF( BX % lo(1) .EQ. 1     ) iX_B(1) = iX_B(1) - swX(1)
        IF( BX % lo(2) .EQ. 1     ) iX_B(2) = iX_B(2) - swX(2)
        IF( BX % lo(3) .EQ. 1     ) iX_B(3) = iX_B(3) - swX(3)
        IF( BX % hi(1) .EQ. nX(1) ) iX_E(1) = iX_E(1) + swX(1)
        IF( BX % hi(2) .EQ. nX(2) ) iX_E(2) = iX_E(2) + swX(2)
        IF( BX % hi(3) .EQ. nX(3) ) iX_E(3) = iX_E(3) + swX(3)

        DO iX3 = iX_B(3), iX_E(3)
        DO iX2 = iX_B(2), iX_E(2)
        DO iX1 = iX_B(1), iX_E(1)

          U_K = RESHAPE( U_P(iX1,iX2,iX3,lo(4):hi(4)), [ nDOFX, nF ] )

          U0(:,iX1,iX2,iX3,iLevel) = U_K(:,iField)

        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    ! --- Combine data from different grids ---

    DO iX3 = 1-swX(3), nX(3)+swX(3)
    DO iX2 = 1-swX(2), nX(2)+swX(2)
    DO iX1 = 1-swX(1), nX(1)+swX(1)

      DO iNodeX = 1, nDOFX

        CALL amrex_parallel_reduce_sum( U0(iNodeX,iX1,iX2,iX3,:), nLevels )

        U(iNodeX,iX1,iX2,iX3) = U0(iNodeX,iX1,iX2,iX3,0)

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE CombineGridData


END MODULE MF_UtilitiesModule
