MODULE MF_Euler_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_fort_module,        ONLY: &
    AR => amrex_real
  USE amrex_box_module,         ONLY: &
    amrex_box
  USE amrex_multifab_module,    ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,    ONLY: &
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule,      ONLY: &
    nDOFX, &
    swX
  USE GeometryFieldsModule,     ONLY: &
    nGF
  USE FluidFieldsModule,        ONLY: &
    nCF, &
    nPF, &
    nAF
  USE Euler_UtilitiesModule,    ONLY: &
    ComputeTimeStep_Euler, &
    ComputeFromConserved_Euler
  USE EquationOfStateModule,    ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  ! --- Local Modules ---

  USE InputParsingModule,       ONLY: &
    nLevels
  USE MF_UtilitiesModule,       ONLY: &
    amrex2thornado_Euler, &
    thornado2amrex_Euler
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler,        &
    TimersStop_AMReX_Euler,         &
    Timer_AMReX_Euler_Allocate, &
    Timer_AMReX_ComputeTimeStep_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeFromConserved
  PUBLIC :: MF_ComputeTimeStep


CONTAINS


  SUBROUTINE MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(in)    :: &
      MF_uGF(0:nLevels-1), MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uPF(0:nLevels-1), MF_uAF(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    REAL(AR), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: A(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uAF => MF_uAF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF) )

        ALLOCATE( P(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nPF) )

        ALLOCATE( A(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nAF) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_Euler( nGF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_Euler( nCF, iX_B1, iX_E1, uCF, U )

        CALL amrex2thornado_Euler( nPF, iX_B1, iX_E1, uPF, P )

        CALL amrex2thornado_Euler( nAF, iX_B1, iX_E1, uAF, A )

        CALL ComputeFromConserved_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

        CALL thornado2amrex_Euler( nPF, iX_B1, iX_E1, uPF, P )

        CALL thornado2amrex_Euler( nAF, iX_B1, iX_E1, uAF, A )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( A )

        DEALLOCATE( P )

        DEALLOCATE( U )

        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ComputeFromConserved


  SUBROUTINE MF_ComputeTimeStep( MF_uGF, MF_uCF, CFL, TimeStepMin )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:nLevels-1), &
                                         MF_uCF(0:nLevels-1)
    REAL(AR),             INTENT(in)  :: CFL
    REAL(AR),             INTENT(out) :: TimeStepMin(0:nLevels-1)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(AR), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(AR), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    REAL(AR) :: TimeStep(0:nLevels-1)

    CALL TimersStart_AMReX_Euler( Timer_AMReX_ComputeTimeStep_Euler )

    TimeStepMin = HUGE( 1.0e0_AR )

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi
        iX_B1 = BX % lo - swX
        iX_E1 = BX % hi + swX

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        ALLOCATE( G(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nGF) )

        ALLOCATE( U(1:nDOFX,iX_B1(1):iX_E1(1), &
                            iX_B1(2):iX_E1(2), &
                            iX_B1(3):iX_E1(3),1:nCF) )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_Euler( nGF, iX_B1, iX_E1, uGF, G )

        CALL amrex2thornado_Euler( nCF, iX_B1, iX_E1, uCF, U )

        CALL ComputeTimeStep_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep(iLevel) )

        TimeStepMin(iLevel) = MIN( TimeStepMin(iLevel), TimeStep(iLevel) )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( U )

        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStepMin, nLevels )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_ComputeTimeStep_Euler )

  END SUBROUTINE MF_ComputeTimeStep


END MODULE MF_Euler_UtilitiesModule
