MODULE MF_AccretionShockUtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_multifab_module, ONLY: &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    swX, &
    nDimsX, &
    nNodesX
  USE FluidFieldsModule, ONLY: &
    nPF, &
    iPF_D, &
    iPF_V1, &
    nAF, &
    iAF_P

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero
  USE MF_FieldsModule_Geometry, ONLY: &
    MF_uGF
  USE MF_FieldsModule_Euler, ONLY: &
    MF_uCF, &
    MF_uPF, &
    MF_uAF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeFromConserved_Euler_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    nX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WriteNodal1DICToFile_SAS

  LOGICAL,          PUBLIC              :: WriteNodal1DIC_SAS
  CHARACTER(LEN=:), PUBLIC, ALLOCATABLE :: FileName_Nodal1DIC_SAS

CONTAINS


  SUBROUTINE WriteNodal1DICToFile_SAS

    INTEGER           :: iNX, iX1, iX2, iX3
    CHARACTER(LEN=16) :: FMT

    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: iX_B1(3), iX_E1(3)
    INTEGER :: iX_B (3), iX_E (3)
    INTEGER :: iLo  (3), iHi  (3)
    INTEGER :: iLevel

    REAL(DP) :: &
      P(1:nNodesX(1),0-swX(1):nX(1)-1+swX(1),1:1,1:1,1:nAF), &
      A(1:nNodesX(1),0-swX(1):nX(1)-1+swX(1),1:1,1:1,1:nPF)

    TYPE(amrex_parmparse) :: PP

    IF( .NOT. WriteNodal1DIC_SAS ) RETURN

    IF( .NOT. nDimsX .EQ. 1 ) RETURN

    IF( nLevels .GT. 1 )THEN

      IF( amrex_parallel_ioprocessor() ) &
        WRITE(*,*) &
          'WARNING: WriteNodal1DICToFile_SAS untested with multi-level mesh'

    END IF

    CALL ComputeFromConserved_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uPF, MF_uAF, &
             swXX_Option = swX )

    iLevel = 0

    P = Zero
    A = Zero

    CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

    DO WHILE( MFI % next() )

      uPF => MF_uPF(iLevel) % DataPtr( MFI )
      uAF => MF_uAF(iLevel) % DataPtr( MFI )

      BX = MFI % tilebox()

      iX_B0 = BX % lo
      iX_E0 = BX % hi
      iX_B1 = BX % lo - swX
      iX_E1 = BX % hi + swX

      iX_B = iX_B0
      iX_E = iX_E0

      IF( BX % lo(1) .EQ. amrex_geom(iLevel) % domain % lo(1) ) &
        iX_B(1) = iX_B1(1)
      IF( BX % hi(1) .EQ. amrex_geom(iLevel) % domain % hi(1) ) &
        iX_E(1) = iX_E1(1)

      DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
      DO iX1 = iX_B(1), iX_E(1)
      DO iNX = 1     , nNodesX(1)

        P(iNX,iX1,iX2,iX3,iPF_D ) &
          = uPF(iX1,iX2,iX3,nNodesX(1)*(iPF_D -1)+iNX)

        P(iNX,iX1,iX2,iX3,iPF_V1) &
          = uPF(iX1,iX2,iX3,nNodesX(1)*(iPF_V1-1)+iNX)

        A(iNX,iX1,iX2,iX3,iAF_P ) &
          = uAF(iX1,iX2,iX3,nNodesX(1)*(iAF_P -1)+iNX)

      END DO
      END DO
      END DO
      END DO

    END DO ! WHILE( MFI % next() )

    CALL amrex_mfiter_destroy( MFI )

    iLo = [ -1   , 1, 1 ]
    iHi = [ nX(1), 1, 1 ]

    DO iX3 = iLo(3), iHi(3)
    DO iX2 = iLo(2), iHi(2)
    DO iX1 = iLo(1), iHi(1)
    DO iNX = 1     , nNodesX(1)

      CALL amrex_parallel_reduce_sum( P(iNX,iX1,iX2,iX3,iPF_D ) )
      CALL amrex_parallel_reduce_sum( P(iNX,iX1,iX2,iX3,iPF_V1) )
      CALL amrex_parallel_reduce_sum( A(iNX,iX1,iX2,iX3,iAF_P ) )

    END DO
    END DO
    END DO
    END DO

    CALL amrex_parmparse_build( PP, 'SAS' )
      CALL PP % get  ( 'FileName_Nodal1DIC_SAS', &
                        FileName_Nodal1DIC_SAS )
    CALL amrex_parmparse_destroy( PP )

    IF( amrex_parallel_ioprocessor() )THEN

      OPEN( UNIT = 101, FILE = TRIM( FileName_Nodal1DIC_SAS ) )

      WRITE(FMT,'(A3,I3.3,A10)') '(SP', nNodesX(1), 'ES25.16E3)'

      WRITE(101,'(ES24.16E3)') Zero
      WRITE(101,'(ES24.16E3)') Zero
      WRITE(101,'(A)') FMT

      DO iX3 = iLo(3), iHi(3)
      DO iX2 = iLo(2), iHi(2)
      DO iX1 = iLo(1), iHi(1)

        WRITE(101,FMT) P(:,iX1,iX2,iX3,iPF_D )
        WRITE(101,FMT) P(:,iX1,iX2,iX3,iPF_V1)
        WRITE(101,FMT) A(:,iX1,iX2,iX3,iAF_P )

      END DO
      END DO
      END DO

      CLOSE( 101 )

    END IF

  END SUBROUTINE WriteNodal1DICToFile_SAS


END MODULE MF_AccretionShockUtilitiesModule
