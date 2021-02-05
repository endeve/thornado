MODULE MF_AccretionShockDiagnosticsModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,               ONLY: &
    AR => amrex_real
  USE amrex_box_module,                ONLY: &
    amrex_box
  USE amrex_multifab_module,           ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,           ONLY: &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule,             ONLY: &
    swX, &
    nDOFX
  USE FluidFieldsModule,               ONLY: &
    nPF, &
    nAF
  USE AccretionShockDiagnosticsModule, ONLY: &
    ComputeAccretionShockDiagnostics

  ! --- Local Modules ---

  USE MF_UtilitiesModule,              ONLY: &
    amrex2thornado_Euler
  USE InputParsingModule,              ONLY: &
    nLevels, &
    DEBUG
  USE TimersModule_AMReX_Euler,        ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler,  &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeAccretionShockDiagnostics


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

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    INTEGER, PARAMETER :: nLegModes = 3
    INTEGER            :: iLegMode
    REAL(AR)           :: Power_Legendre(0:nLevels-1,0:nLegModes-1)

    Power_Legendre = 0.0_AR

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uPF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uAF => MF_uAF(iLevel) % DataPtr( MFI )

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

        CALL amrex2thornado_Euler( nPF, iX_B1, iX_E1, uPF, P )

        CALL amrex2thornado_Euler( nAF, iX_B1, iX_E1, uAF, A )

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


END MODULE MF_AccretionShockDiagnosticsModule
