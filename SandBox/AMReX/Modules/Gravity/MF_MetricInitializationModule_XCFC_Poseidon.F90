MODULE MF_MetricInitializationModule_XCFC_Poseidon

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_imultifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_max

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDOFZ, &
    iE_E0, &
    iE_B0, &
    iE_E1, &
    iE_B1, &
    nE, &
    xR, &
    swX
  USE GeometryFieldsModule, ONLY: &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33, &
    iGF_Psi
  USE GeometryComputationModule, ONLY: &
    LapseFunction, &
    ConformalFactor
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR, &
    nSpecies
  USE XCFC_UtilitiesModule, ONLY: &
    nGS, &
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
    nMF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler_MF
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling
  USE AverageDownModule, ONLY: &
    AverageDown
  USE MF_XCFC_UtilitiesModule, ONLY: &
    swXX, &
    MultiplyWithPsi6_MF, &
    UpdateConformalFactorAndMetric_XCFC_MF, &
    UpdateLapseShiftCurvature_XCFC_MF, &
    ApplyBoundaryConditions_Geometry_XCFC_MF, &
    ComputeGravitationalMass_MF
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE MF_GravitySolutionModule_XCFC, ONLY: &
    ComputeConformalFactor_XCFC_MF, &
    ComputeLapseShiftCurvature_XCFC_MF
  USE MF_XCFC_UtilitiesModule, ONLY: &
    ComputeConformalFactorSourcesAndMg_XCFC_MF, &
    ComputePressureTensorTrace_XCFC_MF


#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  ! --- Poseidon Modules ---

  USE Poseidon_Interface_Boundary_Conditions, ONLY : &
    Poseidon_Set_Uniform_Boundary_Conditions
  USE Poseidon_Interface_Source_Input, ONLY: &
    Poseidon_Input_Sources_Part1
  USE Poseidon_Interface_Run, ONLY: &
    Poseidon_XCFC_Run_Part1
  USE Poseidon_Interface_Initial_Guess, ONLY: &
    Poseidon_Input_Initial_Guess

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeMetric_XCFC_MF_Poseidon
  PUBLIC :: InitializeMetricFromCheckpoint_XCFC_MF_Poseidon

CONTAINS


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE InitializeMetric_XCFC_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uPF, MF_uAF )

#else

  SUBROUTINE InitializeMetric_XCFC_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
#endif
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAF(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    TYPE(amrex_multifab) :: MF_uGS(0:nLevels-1)
    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1)
    TYPE(amrex_multifab) :: LF1   (0:nLevels-1)
    TYPE(amrex_multifab) :: LF2   (0:nLevels-1)
    TYPE(amrex_multifab) :: dLF   (0:nLevels-1)
    TYPE(amrex_multifab) :: CF1   (0:nLevels-1)
    TYPE(amrex_multifab) :: CF2   (0:nLevels-1)
    TYPE(amrex_multifab) :: dCF   (0:nLevels-1)

    LOGICAL  :: CONVERGED
    INTEGER  :: iLevel, ITER, iNX
    REAL(DP) :: MaxLF, MaxCF

    REAL(DP), PARAMETER :: TOLERANCE = 1.0e-13_DP
    INTEGER , PARAMETER :: MAX_ITER = 10

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nGS, swXX )
      CALL MF_uGS(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nMF, swXX )
      CALL MF_uMF(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( LF1(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, swXX )
      CALL LF1(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( LF2(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, swXX )
      CALL LF2(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( dLF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, swXX )
      CALL dLF(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( CF1(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, swXX )
      CALL CF1(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( CF2(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, swXX )
      CALL CF2(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( dCF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, swXX )
      CALL dCF(iLevel) % SetVal( Zero )

    END DO ! iLevel = 0, nLevels-1

    ! --- Iterate to incorporate gravity in initial conditions ---

    CONVERGED = .FALSE.
    ITER = 0

    DO WHILE( .NOT. CONVERGED )

      MaxLF = -HUGE( One )
      MaxCF = -HUGE( One )

      ITER = ITER + 1

      DO iLevel = 0, nLevels - 1

        CALL LF1(iLevel) % COPY &
              ( MF_uGF(iLevel), nDOFX*(iGF_Alpha-1)+1, 1, nDOFX, 0 )

        CALL CF1(iLevel) % COPY &
              ( MF_uGF(iLevel), nDOFX*(iGF_Psi  -1)+1, 1, nDOFX, 0 )

      END DO

      CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCF, +1 )

#ifndef THORNADO_NOTRANSPORT

      CALL MultiplyWithPsi6_MF( iE_B0, iE_E0, iE_B1, iE_E1, MF_uGF, MF_uCR, +1 )

      CALL ComputeConformalFactor &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

      CALL ComputeLapseShiftCurvature &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

      CALL MultiplyWithPsi6_MF( iE_B0, iE_E0, iE_B1, iE_E1, MF_uGF, MF_uCR, -1 )

#else

      CALL ComputeConformalFactor( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

      CALL ComputeLapseShiftCurvature( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

#endif

      CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCF, -1 )

      DO iLevel = 0, nLevels-1

        CALL LF2(iLevel) % COPY &
              ( MF_uGF(iLevel), nDOFX*(iGF_Alpha-1)+1, 1, nDOFX, 0 )

        CALL CF2(iLevel) % COPY &
              ( MF_uGF(iLevel), nDOFX*(iGF_Psi  -1)+1, 1, nDOFX, 0 )

        CALL dLF(iLevel) &
               % LinComb( +One, LF2(iLevel), 1, &
                          -One, LF1(iLevel), 1, 1, &
                          nDOFX, 0 )

        CALL dCF(iLevel) &
               % LinComb( +One, CF2(iLevel), 1, &
                          -One, CF1(iLevel), 1, 1, &
                          nDOFX, 0 )

        DO iNX = 1, nDOFX

          MaxLF = MAX( MaxLF, dLF(iLevel) % Norm0( iNX ) )
          MaxCF = MAX( MaxCF, dCF(iLevel) % Norm0( iNX ) )

        END DO

      END DO ! iLevel = 0, nLevels-1

      CALL amrex_parallel_reduce_max( MaxLF )
      CALL amrex_parallel_reduce_max( MaxCF )

      CALL ComputeConserved_Euler_MF( MF_uGF, MF_uPF, MF_uAF, MF_uCF )

      IF( amrex_parallel_ioprocessor() )THEN

        IF( ITER .GT. MAX_ITER - 3 ) &
          WRITE(*,'(4x,A,I2.2,1x,ES24.16E3)') &
            'ITER, MAX( MaxLF, MaxCF ): ', ITER, MAX( MaxLF, MaxCF )

      END IF

      IF( MAX( MaxLF, MaxCF ) .LT. TOLERANCE ) CONVERGED = .TRUE.

      IF( ITER .EQ. MAX_ITER ) &
        CALL DescribeError_MF &
               ( 902, Real_Option = [ MAX( MaxLF, MaxCF ) ] )

    END DO ! WHILE( .NOT. CONVERGED )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( dCF   (iLevel) )
      CALL amrex_multifab_destroy( CF2   (iLevel) )
      CALL amrex_multifab_destroy( CF1   (iLevel) )
      CALL amrex_multifab_destroy( dLF   (iLevel) )
      CALL amrex_multifab_destroy( LF2   (iLevel) )
      CALL amrex_multifab_destroy( LF1   (iLevel) )
      CALL amrex_multifab_destroy( MF_uMF(iLevel) )
      CALL amrex_multifab_destroy( MF_uGS(iLevel) )

    END DO

#endif

  END SUBROUTINE InitializeMetric_XCFC_MF_Poseidon


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE InitializeMetricFromCheckpoint_XCFC_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uCR )

#else

  SUBROUTINE InitializeMetricFromCheckpoint_XCFC_MF_Poseidon &
    ( MF_uGF, MF_uCF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:)
#endif

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    TYPE(amrex_multifab)  :: MF_uGS    (0:nLevels-1), &
                             MF_uMF    (0:nLevels-1), &
                             MF_uCF_tmp(0:nLevels-1)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab)  :: MF_uCR_tmp(0:nLevels-1)
#endif

    INTEGER          :: iLevel
    REAL(DP)         :: GravitationalMass, Psi_xR, AlphaPsi_xR, Beta_u_xR(3)
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uCF_tmp(iLevel), MF_uCF(iLevel) % BA, MF_uCF(iLevel) % DM, &
               nDOFX * nCF, swX )
      CALL MF_uCF_tmp(iLevel) % SetVal( Zero )

      CALL MF_uCF_tmp(iLevel) % COPY &
             ( MF_uCF(iLevel), 1, 1, nDOFX * nCF, swX )

#ifndef THORNADO_NOTRANSPORT

      CALL amrex_multifab_build &
             ( MF_uCR_tmp(iLevel), MF_uCR(iLevel) % BA, MF_uCR(iLevel) % DM, &
               nDOFZ * nCR * nE * nSpecies, swX )
      CALL MF_uCR_tmp(iLevel) % SetVal( Zero )

      CALL MF_uCR_tmp(iLevel) % COPY &
             ( MF_uCR(iLevel), 1, 1, nDOFX * nCR, nDOFZ * nCR * nE * nSpecies )

#endif

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nGS, swXX )
      CALL MF_uGS(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nMF, swXX )
      CALL MF_uMF(iLevel) % SetVal( Zero )

    END DO

    CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCF_tmp, +1 )

#ifndef THORNADO_NOTRANSPORT

    CALL MultiplyWithPsi6_MF &
           ( iE_B0, iE_E0, iE_B1, iE_E1, MF_uGF, MF_uCR_tmp, +1 )

    CALL ComputeConformalFactor &
           ( MF_uGF, MF_uCF_tmp, MF_uCR_tmp, MF_uGS, MF_uMF )

    CALL ComputeLapseShiftCurvature &
           ( MF_uGF, MF_uCF_tmp, MF_uCR_tmp, MF_uGS, MF_uMF )

#else

    CALL ComputeConformalFactor( MF_uGF, MF_uCF_tmp, MF_uGS, MF_uMF )

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCF_tmp, MF_uGS )

#endif

    CALL Poseidon_Input_Sources_Part1( MF_uGS )

    ! --- Set Boundary Values ---

    CALL ComputeGravitationalMass_MF( MF_uGS, GravitationalMass )

    ! --- Approximate outer boundary with isotropic expressions ---

    Psi_xR       = ConformalFactor( xR(1), GravitationalMass )
    AlphaPsi_xR  = LapseFunction  ( xR(1), GravitationalMass ) * Psi_xR
    Beta_u_xR(1) = Zero
    Beta_u_xR(2) = Zero
    Beta_u_xR(3) = Zero

    INNER_BC_TYPES = [ 'N', 'N', 'N', 'N', 'N' ] ! Neumann
    OUTER_BC_TYPES = [ 'D', 'D', 'D', 'D', 'D' ] ! Dirichlet

    INNER_BC_VALUES &
      = [ Zero  , Zero       , Zero        , Zero        , Zero ]
    OUTER_BC_VALUES &
      = [ Psi_xR, AlphaPsi_xR, Beta_u_xR(1), Beta_u_xR(2), Beta_u_xR(3) ]

    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( 'I', INNER_BC_TYPES, INNER_BC_VALUES )
    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( 'O', OUTER_BC_TYPES, OUTER_BC_VALUES)

    CALL PopulateMF_uMF( MF_uGF, MF_uMF )

    CALL Poseidon_Input_Initial_Guess( MF_uMF )

    CALL Poseidon_XCFC_Run_Part1()

    DO iLevel = 0, nLevels-1

#ifndef THORNADO_NOTRANSPORT
      CALL amrex_multifab_destroy( MF_uCR_tmp(iLevel) )
#endif
      CALL amrex_multifab_destroy( MF_uCF_tmp(iLevel) )
      CALL amrex_multifab_destroy( MF_uMF    (iLevel) )
      CALL amrex_multifab_destroy( MF_uGS    (iLevel) )

    END DO

#endif

  END SUBROUTINE InitializeMetricFromCheckpoint_XCFC_MF_Poseidon


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE ComputeConformalFactor &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

#else

  SUBROUTINE ComputeConformalFactor &
    ( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifndef THORNADO_NOTRANSPORT

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

#else

    CALL ComputeConformalFactorSourcesAndMg_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

#endif

    CALL ComputeConformalFactor_XCFC_MF &
           ( MF_uGS, MF_uMF )

    CALL UpdateConformalFactorAndMetric_XCFC_MF &
           ( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

  END SUBROUTINE ComputeConformalFactor


#ifndef THORNADO_NOTRANSPORT

  SUBROUTINE ComputeLapseShiftCurvature &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

#else

  SUBROUTINE ComputeLapseShiftCurvature &
    ( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

#endif

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
#ifndef THORNADO_NOTRANSPORT
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
#endif
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifndef THORNADO_NOTRANSPORT

    CALL ComputePressureTensorTrace_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

#else

    CALL ComputePressureTensorTrace_XCFC_MF &
           ( MF_uGF, MF_uCF, MF_uGS )
#endif

    CALL ComputeLapseShiftCurvature_XCFC_MF &
           ( MF_uGS, MF_uMF )

    CALL UpdateLapseShiftCurvature_XCFC_MF &
           ( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

    CALL ApplyBoundaryConditions_Geometry_XCFC_MF( MF_uGF )

  END SUBROUTINE ComputeLapseShiftCurvature


  SUBROUTINE PopulateMF_uMF( MF_uGF, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

    TYPE(amrex_imultifab) :: iMF_FineMask
    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uMF     (:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3)

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

#if defined( THORNADO_OMP )
      !$OMP PARALLEL &
      !$OMP PRIVATE( BX, MFI, uGF, uMF, FineMask, iX_B0, iX_E0 )
#endif

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uMF      => MF_uMF(iLevel) % DataPtr( MFI )
        FineMask => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Psi-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Alpha-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_1-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_2-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_Beta_3-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX)

          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_11-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_11-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_12-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_12-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_13-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_13-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_22-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_22-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_23-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_23-1)+iNX)
          uMF    (iX1,iX2,iX3,nDOFX*(iMF_K_dd_33-1)+iNX) &
            = uGF(iX1,iX2,iX3,nDOFX*(iGF_K_dd_33-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyFineMask( iMF_FineMask )

#if defined( THORNADO_OMP )
      !$OMP END PARALLEL
#endif

    END DO ! iLevel = 0, nLevels-1

  END SUBROUTINE PopulateMF_uMF


END MODULE MF_MetricInitializationModule_XCFC_Poseidon
