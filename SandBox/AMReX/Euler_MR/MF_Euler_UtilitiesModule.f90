MODULE MF_Euler_UtilitiesModule

  ! --- AMReX Modules ---

  USE amrex_fort_module, ONLY: &
    amrex_spacedim
  USE amrex_box_module, ONLY: &
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

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF, &
    nPF, &
    nAF
  USE Euler_UtilitiesModule, ONLY: &
    ComputeTimeStep_Euler, &
    ComputeFromConserved_Euler
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    One
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    swX, &
    nX, &
    xL, &
    xR
  USE MF_UtilitiesModule, ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_TimersModule, ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler, &
    Timer_AMReX_Euler_Allocate, &
    Timer_AMReX_ComputeTimeStep_Euler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeFromConserved_Euler_MF
  PUBLIC :: ComputeTimeStep_Euler_MF

CONTAINS


  SUBROUTINE ComputeFromConserved_Euler_MF( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:), MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uPF(0:), MF_uAF(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uPF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uAF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: P(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: A(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uPF => MF_uPF(iLevel) % DataPtr( MFI )
        uAF => MF_uAF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

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

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, U )

        CALL amrex2thornado_X( nPF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uPF, P )

        CALL amrex2thornado_X( nAF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uAF, A )

        CALL ComputeFromConserved_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

        CALL thornado2amrex_X( nPF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uPF, P )

        CALL thornado2amrex_X( nAF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uAF, A )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( A )

        DEALLOCATE( P )

        DEALLOCATE( U )

        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE ComputeFromConserved_Euler_MF


  SUBROUTINE ComputeTimeStep_Euler_MF &
    ( MF_uGF, MF_uCF, CFL, TimeStepMin )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:), MF_uCF(0:)
    REAL(DP),             INTENT(in)  :: CFL
    REAL(DP),             INTENT(out) :: TimeStepMin(0:)

    TYPE(amrex_mfiter) :: MFI
    TYPE(amrex_box)    :: BX

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    REAL(DP), ALLOCATABLE :: G(:,:,:,:,:)
    REAL(DP), ALLOCATABLE :: U(:,:,:,:,:)

    INTEGER :: iLevel, iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), iLo_MF(4)

    REAL(DP) :: TimeStep(0:nLevels-1)

    CALL TimersStart_AMReX_Euler( Timer_AMReX_ComputeTimeStep_Euler )

    TimeStepMin = HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL CreateMesh_MF( iLevel, MeshX )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        iLo_MF = LBOUND( uGF )

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

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uGF, G )

        CALL amrex2thornado_X( nCF, iX_B1, iX_E1, iLo_MF, iX_B0, iX_E0, uCF, U )

        CALL ComputeTimeStep_Euler &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, CFL, TimeStep(iLevel) )

        TimeStepMin(iLevel) = MIN( TimeStepMin(iLevel), TimeStep(iLevel) )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( U )

        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO ! --- Loop over grids (boxes) ---

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyMesh_MF( MeshX )

    END DO ! --- Loop over levels ---

    CALL amrex_parallel_reduce_min( TimeStepMin, SIZE( TimeStepMin ) )

    CALL TimersStop_AMReX_Euler( Timer_AMReX_ComputeTimeStep_Euler )

  END SUBROUTINE ComputeTimeStep_Euler_MF


END MODULE MF_Euler_UtilitiesModule
