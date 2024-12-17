MODULE GravitySolutionModule_XCFC

  USE KindModule, ONLY: &
    DP, &
    Zero
  USE ProgramHeaderModule, ONLY: &
    swX
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE XCFC_UtilitiesModule, ONLY: &
    iMF_Psi, &
    iMF_Alpha, &
    iMF_Beta_1, &
    iMF_Beta_2, &
    iMF_Beta_3, &
    iMF_K_dd_11, &
    iMF_K_dd_12, &
    iMF_K_dd_13, &
    iMF_K_dd_22, &
    iMF_K_dd_23, &
    iMF_K_dd_33, &
    swX_GS
#ifdef GRAVITY_SOLVER_POSEIDON_XCFC
  USE GravitySolutionModule_XCFC_Poseidon, ONLY: &
    InitializeGravitySolver_XCFC_Poseidon, &
    FinalizeGravitySolver_XCFC_Poseidon, &
    ComputeConformalFactor_XCFC_Poseidon, &
    ComputeLapseShiftCurvature_XCFC_Poseidon
#endif
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler,  &
    Timer_GravitySolver

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_XCFC
  PUBLIC :: FinalizeGravitySolver_XCFC
  PUBLIC :: ComputeConformalFactor_XCFC
  PUBLIC :: ComputeLapseShiftCurvature_XCFC

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, FillGhostCells_Option )

    INTEGER , INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL , INTENT(in), OPTIONAL :: FillGhostCells_Option

    LOGICAL :: FillGhostCells

    ! --- Initialize to flat space ---
    CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    FillGhostCells = .FALSE.
    IF( PRESENT( FillGhostCells_Option ) ) &
      FillGhostCells = FillGhostCells_Option

    swX_GS = 0
    IF( FillGhostCells ) &
      swX_GS = swX

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL InitializeGravitySolver_XCFC_Poseidon

#endif

  END SUBROUTINE InitializeGravitySolver_XCFC


  SUBROUTINE FinalizeGravitySolver_XCFC

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL FinalizeGravitySolver_XCFC_Poseidon

#endif

  END SUBROUTINE FinalizeGravitySolver_XCFC


  SUBROUTINE ComputeConformalFactor_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout) :: &
      M (1:,iX_B0(1)-swX_GS(1):,iX_B0(2)-swX_GS(2):,iX_B0(3)-swX_GS(3):,1:)

    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  CALL ComputeConformalFactor_XCFC_Poseidon &
         ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

#else

    M(:,:,:,:,iMF_Psi) = Zero

#endif

    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeConformalFactor_XCFC


  SUBROUTINE ComputeLapseShiftCurvature_XCFC &
    ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      GS(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout) :: &
      M (1:,iX_B0(1)-swX_GS(1):,iX_B0(2)-swX_GS(2):,iX_B0(3)-swX_GS(3):,1:)

    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputeLapseShiftCurvature_XCFC_Poseidon &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GS, M )

#else

    M(:,:,:,:,iMF_Alpha)               = Zero
    M(:,:,:,:,iMF_Beta_1:iMF_Beta_3)   = Zero
    M(:,:,:,:,iMF_K_dd_11:iMF_K_dd_33) = Zero

#endif

    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeLapseShiftCurvature_XCFC


END MODULE GravitySolutionModule_XCFC
