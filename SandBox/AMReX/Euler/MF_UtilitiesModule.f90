!> Module for operations on MultiFabs
!> @todo Modify data transfer subroutine to allow specification of looping over
!>       ghost cells
MODULE MF_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,     ONLY: &
    AR => amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_sum
  USE amrex_geometry_module, ONLY: &
    amrex_geometry
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build

  ! --- thornado Modules ---

  USE ProgramHeaderModule,    ONLY: &
    nDOFX
  USE GeometryFieldsModule,   ONLY: &
    nGF,          &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule,      ONLY: &
    nCF,    &
    iCF_D,  &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E,  &
    iCF_Ne, &
    nPF,    &
    iPF_D,  &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E,  &
    iPF_Ne
  USE Euler_UtilitiesModule, ONLY: &
    ComputePrimitive_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive

  ! --- Local Modules ---

  USE MyAmrModule,                       ONLY: &
    nLevels, &
    nX,      &
    swX
  USE MF_Euler_BoundaryConditionsModule, ONLY: &
    EdgeMap,          &
    ConstructEdgeMap, &
    MF_ApplyBoundaryConditions_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: AMReX2thornado
  PUBLIC :: thornado2AMReX
  PUBLIC :: ShowVariableFromMultiFab
  PUBLIC :: WriteRawDataToFile


CONTAINS


  SUBROUTINE AMReX2thornado( nVars, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nVars
    INTEGER,  INTENT(in)  :: iX_B(3), iX_E(3)
    REAL(AR), INTENT(in)  :: &
      Data_amrex   (   iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(AR), INTENT(out) :: &
      Data_thornado(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)

    INTEGER :: iX1, iX2, iX3, iVar

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      DO iVar = 1, nVars
        Data_thornado(1:nDOFX,iX1,iX2,iX3,iVar) &
          = Data_amrex(iX1,iX2,iX3,nDOFX*(iVar-1)+1:nDOFX*iVar)
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE AMReX2thornado


  SUBROUTINE thornado2AMReX( nVars, iX_B, iX_E, Data_amrex, Data_thornado )

    INTEGER,  INTENT(in)  :: nVars
    INTEGER,  INTENT(in)  :: iX_B(3), iX_E(3)
    REAL(AR), INTENT(out) :: &
      Data_amrex   (   iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(AR), INTENT(in)  :: &
      Data_thornado(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)

    INTEGER :: iX1, iX2, iX3, iVar

    DO iX3 = iX_B(3), iX_E(3)
    DO iX2 = iX_B(2), iX_E(2)
    DO iX1 = iX_B(1), iX_E(1)

      DO iVar = 1, nVars
        Data_amrex(iX1,iX2,iX3,nDOFX*(iVar-1)+1:nDOFX*iVar) &
          = Data_thornado(1:nDOFX,iX1,iX2,iX3,iVar)
      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE thornado2AMReX


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

    END DO

    WRITE(*,*)

  END SUBROUTINE ShowVariableFromMultiFab


  SUBROUTINE WriteRawDataToFile( GEOM, MF_uGF, MF_uCF )

    TYPE(amrex_geometry), INTENT(in) :: GEOM   (0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)

    INTEGER            :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    INTEGER            :: iLevel, nCompGF, nCompCF
    INTEGER            :: iX1, iX2, iX3, iNode, iCF, iGF, nLines
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI
    CHARACTER(LEN=16)  :: FMT

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR), ALLOCATABLE         :: G  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: U  (:,:,:,:,:)
    TYPE(EdgeMap)                 :: Edge_Map

    REAL(AR)                      :: P       (nDOFX,nPF)
    REAL(AR)                      :: Pressure(nDOFX)

    ALLOCATE( G(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nGF) )

    ALLOCATE( U(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nCF) )

    G = 0.0_AR
    U = 0.0_AR

    ! --- Convert AMReX MultiFabs to thornado arrays ---

    DO iLevel = 0, nLevels-1

      ! --- Apply boundary conditions to interior domains ---

      CALL MF_uGF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL MF_uCF(iLevel) % Fill_Boundary( GEOM(iLevel) )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF    => MF_uGF(iLevel) % DataPtr( MFI )
        nCompGF = MF_uGF(iLevel) % nComp()

        uCF    => MF_uCF(iLevel) % DataPtr( MFI )
        nCompCF = MF_uCF(iLevel) % nComp()

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL AMReX2thornado( nGF, iX_B1, iX_E1, &
                             uGF(iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1:nCompGF), &
                             G  (1:nDOFX, &
                                 iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1:nGF) )

        CALL AMReX2thornado( nCF, iX_B1, iX_E1, &
                             uCF(iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1:nCompCF), &
                             U  (1:nDOFX, &
                                 iX_B1(1):iX_E1(1), &
                                 iX_B1(2):iX_E1(2), &
                                 iX_B1(3):iX_E1(3),1:nCF) )

        ! --- Apply boundary conditions ---

        CALL ConstructEdgeMap( GEOM(iLevel), BX, Edge_Map )

        CALL MF_ApplyBoundaryConditions_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                  U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF), Edge_Map )

      END DO

    END DO

    ! --- Combine data from all processes ---

    DO iX3   = 1-swX(3), nX(3)+swX(3)
    DO iX2   = 1-swX(2), nX(2)+swX(2)
    DO iX1   = 1-swX(1), nX(1)+swX(1)
    DO iNode = 1       , nDOFX

      DO iGF = 1, nGF

        CALL amrex_parallel_reduce_sum( G(iNode,iX1,iX2,iX3,iGF) )

      END DO

      DO iCF = 1, nCF

        CALL amrex_parallel_reduce_sum( U(iNode,iX1,iX2,iX3,iCF) )

      END DO

    END DO
    END DO
    END DO
    END DO

    ! --- Compute primitive and write to files ---

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(FMT,'(A3,I3.3,A10)') '(SP', nDOFX, 'ES25.16E3)'

      OPEN( UNIT = 100, FILE = 'D.dat' )
      OPEN( UNIT = 101, FILE = 'V.dat' )
      OPEN( UNIT = 102, FILE = 'P.dat' )

      WRITE(100,*) FMT
      WRITE(101,*) FMT
      WRITE(102,*) FMT

      nLines = 0

      DO iX3 = 1-swX(3), nX(3)+swX(3)
      DO iX2 = 1-swX(2), nX(2)+swX(2)
      DO iX1 = 1-swX(1), nX(1)+swX(1)

        nLines = nLines + 1

        CALL ComputePrimitive_Euler &
               ( U(:,iX1,iX2,iX3,iCF_D ),       &
                 U(:,iX1,iX2,iX3,iCF_S1),       &
                 U(:,iX1,iX2,iX3,iCF_S2),       &
                 U(:,iX1,iX2,iX3,iCF_S3),       &
                 U(:,iX1,iX2,iX3,iCF_E ),       &
                 U(:,iX1,iX2,iX3,iCF_Ne),       &
                 P(:,iPF_D ),                   &
                 P(:,iPF_V1),                   &
                 P(:,iPF_V2),                   &
                 P(:,iPF_V3),                   &
                 P(:,iPF_E ),                   &
                 P(:,iPF_Ne),                   &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

        CALL ComputePressureFromPrimitive &
               ( P(:,iPF_D ), P(:,iPF_E ), P(:,iPF_Ne), Pressure(:) )

        WRITE(100,TRIM(FMT)) P(:,iPF_D )
        WRITE(101,TRIM(FMT)) P(:,iPF_V1)
        WRITE(102,TRIM(FMT)) Pressure(:)

      END DO
      END DO
      END DO

      WRITE(100,*) nLines
      WRITE(101,*) nLines
      WRITE(102,*) nLines

      CLOSE( 102 )
      CLOSE( 101 )
      CLOSE( 100 )

    END IF

    DEALLOCATE( U )
    DEALLOCATE( G )

  END SUBROUTINE WriteRawDataToFile


END MODULE MF_UtilitiesModule
