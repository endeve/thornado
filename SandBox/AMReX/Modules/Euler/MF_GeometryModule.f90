MODULE MF_GeometryModule

  ! --- AMReX Modules ---

  USE amrex_box_module,                          ONLY: &
    amrex_box
  USE amrex_multifab_module,                     ONLY: &
    amrex_multifab,     &
    amrex_mfiter,       &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module,                     ONLY: &
    amrex_parallel_ioprocessor

  ! --- thornado Modules ---

  USE ProgramHeaderModule,                       ONLY: &
    nDOFX, &
    swX
  USE GeometryFieldsModule,                      ONLY: &
    nGF, &
    iGF_Phi_N
  USE GeometryComputationModule,                 ONLY: &
    ComputeGeometryX
  USE GravitySolutionModule_Newtonian_PointMass, ONLY: &
    ComputeGravitationalPotential

  ! --- Local Modules ---

  USE MF_KindModule,                             ONLY: &
    DP
  USE InputParsingModule,                        ONLY: &
    nLevels, &
    UseTiling
  USE MF_UtilitiesModule,                        ONLY: &
    amrex2thornado_X, &
    thornado2amrex_X
  USE TimersModule_AMReX_Euler,                  ONLY: &
    TimersStart_AMReX_Euler, &
    TimersStop_AMReX_Euler,  &
    Timer_AMReX_Euler_Allocate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MF_ComputeGeometryX
  PUBLIC :: MF_ComputeGravitationalPotential


CONTAINS


  SUBROUTINE MF_ComputeGeometryX( MF_uGF, Mass )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    REAL(DP),             INTENT(in)    :: Mass

    INTEGER                       :: iLevel
    INTEGER                       :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                                     iLo_MF(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G(:,:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

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

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL ComputeGeometryX &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass_Option = Mass )

        CALL thornado2amrex_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ComputeGeometryX


  SUBROUTINE MF_ComputeGravitationalPotential( MF_uGF, Mass )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    REAL(DP),             INTENT(in)    :: Mass

    INTEGER                       :: iLevel
    INTEGER                       :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3), &
                                     iLo_MF(4)
    TYPE(amrex_box)               :: BX
    TYPE(amrex_mfiter)            :: MFI
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), ALLOCATABLE         :: G(:,:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

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

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        CALL amrex2thornado_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL ComputeGravitationalPotential &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, Mass )

        CALL thornado2amrex_X( nGF, iX_B1, iX_E1, iLo_MF, iX_B1, iX_E1, uGF, G )

        CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_Allocate )

        DEALLOCATE( G )

        CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_Allocate )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MF_ComputeGravitationalPotential


END MODULE MF_GeometryModule
