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
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build

  ! --- thornado Modules ---

  USE ProgramHeaderModule,   ONLY: &
    nDOFX
  USE GeometryFieldsModule,  ONLY: &
    nGF
  USE FluidFieldsModule,     ONLY: &
    nCF

  ! --- Local Modules ---

  USE MyAmrModule,           ONLY: &
    nLevels, &
    nX,      &
    swX

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


  SUBROUTINE WriteRawDataToFile( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(in) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in) :: MF_uCF(0:nLevels-1)

    INTEGER            :: iX_B(3), iX_E(3)
    INTEGER            :: iLevel, nCompGF, nCompCF
    INTEGER            :: iX1, iX2, iX3, iNode, iCF, iGF
    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR), ALLOCATABLE         :: G  (:,:,:,:,:)
    REAL(AR), ALLOCATABLE         :: U  (:,:,:,:,:)

    ALLOCATE( G(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nGF) )

    ALLOCATE( U(1:nDOFX,1-swX(1):nX(1)+swX(1), &
                        1-swX(2):nX(2)+swX(2), &
                        1-swX(3):nX(3)+swX(3),1:nCF) )
    G = 0.0_AR
    U = 0.0_AR

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF  => MF_uGF(iLevel) % DataPtr( MFI )
        uCF  => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B = BX % lo
        iX_E = BX % hi

        IF( iX_B(1) .EQ. 1     ) iX_B(1) = iX_B(1) - swX(1)
        IF( iX_B(2) .EQ. 1     ) iX_B(2) = iX_B(2) - swX(2)
        IF( iX_B(3) .EQ. 1     ) iX_B(3) = iX_B(3) - swX(3)

        IF( iX_E(1) .EQ. nX(1) ) iX_E(1) = iX_E(1) + swX(1)
        IF( iX_E(2) .EQ. nX(2) ) iX_E(2) = iX_E(2) + swX(2)
        IF( iX_E(3) .EQ. nX(3) ) iX_E(3) = iX_E(3) + swX(3)

        nCompGF = MF_uGF(iLevel) % nComp()
        nCompCF = MF_uCF(iLevel) % nComp()

        CALL AMReX2thornado( nGF, iX_B, iX_E, &
                             uGF(iX_B(1):iX_E(1), &
                                 iX_B(2):iX_E(2), &
                                 iX_B(3):iX_E(3),1:nCompGF), &
                             G  (1:nDOFX, &
                                 iX_B(1):iX_E(1), &
                                 iX_B(2):iX_E(2), &
                                 iX_B(3):iX_E(3),1:nGF) )

        CALL AMReX2thornado( nCF, iX_B, iX_E, &
                             uCF(iX_B(1):iX_E(1), &
                                 iX_B(2):iX_E(2), &
                                 iX_B(3):iX_E(3),1:nCompCF), &
                             U  (1:nDOFX, &
                                 iX_B(1):iX_E(1), &
                                 iX_B(2):iX_E(2), &
                                 iX_B(3):iX_E(3),1:nCF) )

      END DO

    END DO

    DO iX3 = 1-swX(3), nX(3)+swX(3)
    DO iX2 = 1-swX(2), nX(2)+swX(2)
    DO iX1 = 1-swX(1), nX(1)+swX(1)

      DO iGF   = 1, nGF
      DO iNode = 1, nDOFX

        CALL amrex_parallel_reduce_sum( G(iNode,iX1,iX2,iX3,iGF) )

      END DO
      END DO

      DO iCF   = 1, nCF
      DO iNode = 1, nDOFX

        CALL amrex_parallel_reduce_sum( U(iNode,iX1,iX2,iX3,iCF) )

      END DO
      END DO

    END DO
    END DO
    END DO

    IF( amrex_parallel_ioprocessor() )THEN

      OPEN( UNIT = 100, FILE = 'G.dat' )
      OPEN( UNIT = 101, FILE = 'U.dat' )

      DO iX3 = 1-swX(3), nX(3)+swX(3)
      DO iX2 = 1-swX(2), nX(2)+swX(2)
      DO iX1 = 1-swX(1), nX(1)+swX(1)

        DO iGF   = 1, nGF
        DO iNode = 1, nDOFX

          WRITE(100,'(ES25.16E3)',ADVANCE='NO') G(iNode,iX1,iX2,iX3,iGF)

        END DO
        END DO

        DO iCF   = 1, nCF
        DO iNode = 1, nDOFX

          WRITE(101,'(ES25.16E3)',ADVANCE='NO') U(iNode,iX1,iX2,iX3,iCF)

        END DO
        END DO

      END DO
      END DO
      END DO

      CLOSE( 101 )
      CLOSE( 100 )

    END IF

    DEALLOCATE( U )
    DEALLOCATE( G )

  END SUBROUTINE WriteRawDataToFile


END MODULE MF_UtilitiesModule
