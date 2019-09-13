MODULE MF_Euler_UtilitiesModule

  ! --- AMReX Modules ---
  USE amrex_fort_module,     ONLY: &
    amrex_real
  USE amrex_box_module,      ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_reduce_min, &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---
  USE ProgramHeaderModule,   ONLY: &
    nDOFX, swX
  USE GeometryFieldsModule,  ONLY: &
    nGF
  USE FluidFieldsModule,     ONLY: &
    nCF, nPF, nAF
  USE Euler_UtilitiesModule, ONLY: &
    Euler_ComputeTimeStep, &
    Euler_ComputeFromConserved
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  ! --- Local Modules ---
  USE MyAmrModule, ONLY: &
    nLevels
  USE MF_UtilitiesModule, ONLY: &
    AMReX2thornado, &
    thornado2AMReX
  USE TimersModule_AMReX_Euler, ONLY: &
    TimersStart_AMReX_Euler, TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_DataTransfer, &
    Timer_AMReX_Euler_ComputeTimeStep

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeFromConserved
  PUBLIC :: MF_ComputeTimeStep


CONTAINS


  SUBROUTINE MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(in   ) :: &
      MF_uGF(0:nLevels), MF_uCF(0:nLevels)
    TYPE(amrex_multifab), INTENT(inout) :: &
      MF_uPF(0:nLevels), MF_uAF(0:nLevels)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: A(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3)

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uAF => MF_uAF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )
        ALLOCATE( G(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3), 1:nGF ) )
        ALLOCATE( U(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3), 1:nCF ) )
        ALLOCATE( P(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3), 1:nPF ) )
        ALLOCATE( A(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3), 1:nAF ) )

        CALL AMReX2thornado &
               ( nGF, iX_B0, iX_E0, &
                 uGF(      iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nDOFX*nGF), &
                 G(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nGF) )
        CALL AMReX2thornado &
               ( nCF, iX_B0, iX_E0, &
                 uCF(      iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nDOFX*nCF), &
                 U(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nCF) )
        CALL AMReX2thornado &
               ( nPF, iX_B0, iX_E0, &
                 uPF(      iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nDOFX*nPF), &
                 P(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nPF) )
        CALL AMReX2thornado &
               ( nAF, iX_B0, iX_E0, &
                 uAF(      iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nDOFX*nAF), &
                 A(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nAF) )
        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

        CALL Euler_ComputeFromConserved &
               ( iX_B0, iX_E0, &
                 G(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nGF), &
                 U(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nCF), &
                 P(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nPF), &
                 A(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nAF) )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )
        CALL thornado2AMReX &
               ( nPF, iX_B0, iX_E0, &
                 uPF(      iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nDOFX*nPF), &
                 P(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nPF) )
        CALL thornado2AMReX &
               ( nAF, iX_B0, iX_E0, &
                 uAF(      iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nDOFX*nAF), &
                 A(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nAF) )

        DEALLOCATE( A )
        DEALLOCATE( P )
        DEALLOCATE( U )
        DEALLOCATE( G )
        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ComputeFromConserved


  SUBROUTINE MF_ComputeTimeStep( MF_uGF, MF_uCF, CFL, TimeStepMin )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:nlevels), MF_uCF(0:nLevels)
    REAL(amrex_real),     INTENT(in)  :: CFL
    REAL(amrex_real),     INTENT(out) :: TimeStepMin(0:nLevels)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(amrex_real), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(amrex_real), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(amrex_real), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(amrex_real), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3)

    REAL(amrex_real) :: TimeStep(0:nLevels)

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_ComputeTimeStep )

    TimeStepMin = HUGE( 1.0e0_amrex_real )

    DO iLevel = 0, nLevels

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = .TRUE. )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )
        ALLOCATE( G(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nGF) )
        ALLOCATE( U(1:nDOFX,iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3),1:nCF) )

        CALL AMReX2thornado &
               ( nGF, iX_B0, iX_E0, &
                 uGF(      iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nDOFX*nGF), &
                 G(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nGF) )
        CALL AMReX2thornado &
               ( nCF, iX_B0, iX_E0, &
                 uCF(      iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nDOFX*nCF), &
                 U(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nCF) )
        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

        CALL Euler_ComputeTimeStep &
               ( iX_B0, iX_E0, &
                 G(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nGF), &
                 U(1:nDOFX,iX_B0(1):iX_E0(1), &
                           iX_B0(2):iX_E0(2), &
                           iX_B0(3):iX_E0(3),1:nCF), &
                 CFL, TimeStep(iLevel) )

        TimeStepMin(iLevel) = MIN( TimeStepMin(iLevel), TimeStep(iLevel) )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )
        DEALLOCATE( U )
        DEALLOCATE( G )
        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_DataTransfer )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStepMin, nLevels+1 )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_ComputeTimeStep )

  END SUBROUTINE MF_ComputeTimeStep


END MODULE MF_Euler_UtilitiesModule
