MODULE MF_AccretionShockUtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    AR => amrex_real
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
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDOFX
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
    ComputeAccretionShockDiagnostics
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive

  ! --- Local Modules ---

  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    amrex2thornado_X_Global
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    nX, &
    DEBUG
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeAccretionShockDiagnostics
  PUBLIC :: WriteNodalDataToFile_SAS


CONTAINS


  SUBROUTINE MF_ComputeAccretionShockDiagnostics( MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uAF(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(AR), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    REAL(AR), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: A(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    INTEGER, PARAMETER :: nLegModes = 3
    INTEGER            :: iLegMode
    REAL(AR)           :: Power_Legendre(0:nLevels-1,0:nLegModes-1)

    Power_Legendre = 0.0_AR

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uPF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uAF => MF_uAF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uPF )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( P(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nPF ) )

        ALLOCATE( A(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nAF ) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_X( nPF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uPF, P )

        CALL amrex2thornado_X( nAF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uAF, A )

        CALL ComputeAccretionShockDiagnostics &
               ( iX_B0, iX_E0, iX_B1, iX_E1, P, A, &
                 Power_Legendre(iLevel,0:nLegModes-1) )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( A )

        DEALLOCATE( P )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

    DO iLegMode = 0, nLegModes-1

      CALL amrex_parallel_reduce_sum( Power_Legendre(:,iLegMode), nLevels )

    END DO

  END SUBROUTINE MF_ComputeAccretionShockDiagnostics


  SUBROUTINE WriteNodalDataToFile_SAS( GEOM, MF_uGF, MF_uCF, FileNameBase )

    TYPE(amrex_geometry), INTENT(in) :: GEOM(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)
    CHARACTER(LEN=*)    , INTENT(in) :: FileNameBase

    INTEGER           :: iLo(3), iHi(3), iX1, iX2, iX3
    CHARACTER(LEN=16) :: FMT

    REAL(AR) :: P(1:nDOFX,1:nPF)
    REAL(AR) :: A(1:nDOFX,1:nAF)
    REAL(AR) :: G(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                          1-swX(2):nX(2)+swX(2), &
                          1-swX(3):nX(3)+swX(3), &
                  1:nGF)
    REAL(AR) :: U(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                          1-swX(2):nX(2)+swX(2), &
                          1-swX(3):nX(3)+swX(3), &
                  1:nCF)

    CALL amrex2thornado_X_Global &
           ( GEOM, MF_uGF, nGF, G, ApplyBC_Option = .FALSE. )

    CALL amrex2thornado_X_Global &
           ( GEOM, MF_uCF, nCF, U, ApplyBC_Option = .TRUE. )

    IF( amrex_parallel_ioprocessor() )THEN

      iLo = 1  - swX
      iHi = nX + swX

      OPEN( UNIT = 101, FILE = TRIM( FileNameBase ) // '_D.dat' )
      OPEN( UNIT = 102, FILE = TRIM( FileNameBase ) // '_V.dat' )
      OPEN( UNIT = 103, FILE = TRIM( FileNameBase ) // '_P.dat' )

      WRITE(FMT,'(A3,I3.3,A10)') '(SP', nDOFX, 'ES25.16E3)'

      WRITE(101,'(A)') FMT
      WRITE(102,'(A)') FMT
      WRITE(103,'(A)') FMT

      DO iX3 = iLo(3), iHi(3)
      DO iX2 = iLo(2), iHi(2)
      DO iX1 = iLo(1), iHi(1)

        CALL ComputePrimitive_Euler &
               ( U(:,iX1,iX2,iX3,iCF_D ), &
                 U(:,iX1,iX2,iX3,iCF_S1), &
                 U(:,iX1,iX2,iX3,iCF_S2), &
                 U(:,iX1,iX2,iX3,iCF_S3), &
                 U(:,iX1,iX2,iX3,iCF_E ), &
                 U(:,iX1,iX2,iX3,iCF_Ne), &
                 P(:,iPF_D ), &
                 P(:,iPF_V1), &
                 P(:,iPF_V2), &
                 P(:,iPF_V3), &
                 P(:,iPF_E ), &
                 P(:,iPF_Ne), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( P(:,iPF_D ), P(:,iPF_E ), P(:,iPF_Ne), A(:,iAF_P) )

        WRITE(101,FMT) P(:,iPF_D )
        WRITE(102,FMT) P(:,iPF_V1)
        WRITE(103,FMT) A(:,iAF_P )

      END DO
      END DO
      END DO

      CLOSE( 103 )
      CLOSE( 102 )
      CLOSE( 101 )

    END IF

  END SUBROUTINE WriteNodalDataToFile_SAS


END MODULE MF_AccretionShockUtilitiesModule
