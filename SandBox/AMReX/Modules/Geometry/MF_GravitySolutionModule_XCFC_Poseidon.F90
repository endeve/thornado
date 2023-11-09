MODULE MF_GravitySolutionModule_XCFC_Poseidon

  ! --- AMReX Modules ---

  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_max

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX, &
    iE_E0, &
    iE_B0, &
    iE_E1, &
    iE_B1
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    iGF_Alpha, &
    iGF_Psi
  USE GeometryComputationModule, ONLY: &
    LapseFunction, &
    ConformalFactor
  USE FluidFieldsModule, ONLY: &
    nCF
  USE XCFC_UtilitiesModule, ONLY: &
    nGS, &
    nMF

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler_MF
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    xR, &
    swX, &
    iRestart
  USE AverageDownModule, ONLY: &
    AverageDown
  USE MF_XCFC_UtilitiesModule, ONLY: &
    swXX, &
    MultiplyWithPsi6_MF, &
    UpdateConformalFactorAndMetric_XCFC_MF, &
    UpdateLapseShiftCurvature_XCFC_MF, &
    ApplyBoundaryConditions_Geometry_XCFC_MF, &
    ComputeGravitationalMass_MF, &
    PopulateMF_uMF

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

  USE MF_Euler_XCFC_UtilitiesModule, ONLY: &
    ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF, &
    ComputePressureTensorTrace_XCFC_Euler_MF

#ifndef THORNADO_NOTRANSPORT
  USE MF_TwoMoment_XCFC_UtilitiesModule, ONLY: &
    ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF, &
    ComputePressureTensorTrace_XCFC_TwoMoment_MF
#endif

  ! --- Poseidon Modules ---

  USE Poseidon_Interface_Initialization, ONLY: &
    Initialize_Poseidon
  USE Poseidon_Interface_Boundary_Conditions, ONLY : &
    Poseidon_Set_Uniform_Boundary_Conditions
  USE Poseidon_Interface_Source_Input, ONLY: &
    Poseidon_Input_Sources_Part1, &
    Poseidon_Input_Sources_Part2
  USE Poseidon_Interface_Run, ONLY: &
    Poseidon_XCFC_Run_Part1, &
    Poseidon_XCFC_Run_Part2
  USE Poseidon_Interface_Return_Routines, ONLY: &
    Poseidon_Return_Conformal_Factor, &
    Poseidon_Return_ALL
  USE Poseidon_Interface_Close, ONLY: &
    Poseidon_Close
  USE Poseidon_Interface_Initial_Guess, ONLY: &
    Poseidon_Input_Initial_Guess

#endif

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: FillGhostCells

  PUBLIC :: InitializeGravitySolver_XCFC_MF_Poseidon
  PUBLIC :: FinalizeGravitySolver_XCFC_MF_Poseidon
  PUBLIC :: ComputeConformalFactor_XCFC_MF_Poseidon
  PUBLIC :: ComputeLapseShiftCurvature_XCFC_MF_Poseidon

  PUBLIC :: InitializeMetric_Euler_MF_Poseidon
  PUBLIC :: InitializeMetric_TwoMoment_MF_Poseidon

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC_MF_Poseidon( MF_uGF, MF_uCF )

    TYPE(amrex_multifab), INTENT(inout), OPTIONAL :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)   , OPTIONAL :: MF_uCF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    TYPE(amrex_multifab)  :: MF_uGS    (0:nLevels-1), &
                             MF_uMF    (0:nLevels-1), &
                             MF_uCF_tmp(0:nLevels-1)
    TYPE(amrex_parmparse) :: PP

    INTEGER          :: iLevel
    REAL(DP)         :: GravitationalMass, Psi_xR, AlphaPsi_xR, Beta_u_xR(3)
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    FillGhostCells = .FALSE.
    CALL amrex_parmparse_build( PP, 'poseidon' )
      CALL PP % query( 'FillGhostCells', FillGhostCells )
    CALL amrex_parmparse_destroy( PP )

    swXX = 0
    IF( FillGhostCells ) swXX = swX

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)') &
        'INFO: Gravity Solver (Poseidon, XCFC)'
      WRITE(*,'(4x,A)') &
        '-------------------------------------'
      WRITE(*,*)
      WRITE(*,'(6x,A,L)') 'FillGhostCells: ', FillGhostCells
      WRITE(*,'(6x,A)') 'Only implemented for 1D spherical symmetry.'
      WRITE(*,*)

    END IF

    CALL Initialize_Poseidon &
           ( Source_NQ                    = nNodesX,          &
             Source_xL                    = [ -Half, +Half ], &
             Source_RQ_xlocs              = MeshX(1) % Nodes, &
             Source_TQ_xlocs              = MeshX(2) % Nodes, &
             Source_PQ_xlocs              = MeshX(3) % Nodes, &
             Source_Units                 = 'G',              &
             Source_Radial_Boundary_Units = 'km',             &
             Flat_Guess_Option            = .TRUE.,           &
             Verbose_Option               = .FALSE.,          &
             Print_Setup_Option           = .TRUE.,           &
             Print_Results_Option         = .FALSE. )

    IF( iRestart .GE. 0 )THEN

      DO iLevel = 0, nLevels-1

        CALL amrex_multifab_build &
               ( MF_uCF_tmp(iLevel), MF_uCF(iLevel) % BA, MF_uCF(iLevel) % DM, &
                 nDOFX * nCF, swX )
        CALL MF_uCF_tmp(iLevel) % SetVal( Zero )

        CALL MF_uCF_tmp(iLevel) % COPY &
               ( MF_uCF(iLevel), 1, 1, nDOFX * nCF, swX )

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

      CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF &
             ( MF_uGF, MF_uCF_tmp, MF_uGS )

      CALL MultiplyWithPsi6_MF( MF_uGF, MF_uCF_tmp, -1 )

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

        CALL amrex_multifab_destroy( MF_uCF_tmp(iLevel) )
        CALL amrex_multifab_destroy( MF_uMF    (iLevel) )
        CALL amrex_multifab_destroy( MF_uGS    (iLevel) )

      END DO

    END IF ! iRestart .GE. 0

#endif

  END SUBROUTINE InitializeGravitySolver_XCFC_MF_Poseidon


  SUBROUTINE FinalizeGravitySolver_XCFC_MF_Poseidon

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_XCFC_MF_Poseidon


  SUBROUTINE ComputeConformalFactor_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

    REAL(DP)         :: GravitationalMass, Psi_xR, AlphaPsi_xR, Beta_u_xR(3)
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    ! --- Set Boundary Values ---

    CALL ComputeGravitationalMass_MF( MF_uGS, GravitationalMass )

    ! --- Approximate outer boundary with isotropic expressions ---

    Psi_xR       = ConformalFactor( xR(1), GravitationalMass )
    AlphaPsi_xR  = LapseFunction  ( xR(1), GravitationalMass ) * Psi_xR
    Beta_u_xR(1) = Zero
    Beta_u_xR(2) = Zero
    Beta_u_xR(3) = Zero

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

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

    ! --- Set XCFC sources with current conformal factor ---
    CALL Poseidon_Input_Sources_Part1( MF_uGS )

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

    CALL Poseidon_Return_Conformal_Factor &
           ( MF_uMF, FillGhostCells_Option = FillGhostCells )

#endif

  END SUBROUTINE ComputeConformalFactor_XCFC_MF_Poseidon


  SUBROUTINE ComputeLapseShiftCurvature_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:nLevels-1) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:nLevels-1)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ! --- Set gravity sources with updated conformal factor ---

    CALL Poseidon_Input_Sources_Part2( MF_uGS )

    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()

    CALL Poseidon_Return_ALL &
           ( MF_uMF, FillGhostCells_Option = FillGhostCells )

#endif

  END SUBROUTINE ComputeLapseShiftCurvature_XCFC_MF_Poseidon


  SUBROUTINE InitializeMetric_Euler_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAF(0:nLevels-1)

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

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

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

      CALL ComputeConformalFactor_Euler( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

      CALL ComputeLapseShiftCurvature_Euler( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

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

  END SUBROUTINE InitializeMetric_Euler_MF_Poseidon


  SUBROUTINE InitializeMetric_TwoMoment_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAF(0:nLevels-1)

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

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nGS, 0 )
      CALL MF_uGS(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nMF, 0 )
      CALL MF_uMF(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( LF1(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )
      CALL LF1(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( LF2(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )
      CALL LF2(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( dLF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )
      CALL dLF(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( CF1(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )
      CALL CF1(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( CF2(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )
      CALL CF2(iLevel) % SetVal( Zero )

      CALL amrex_multifab_build &
             ( dCF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )
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

      CALL MultiplyWithPsi6_MF( iE_B0, iE_E0, iE_B1, iE_E1, MF_uGF, MF_uCR, +1 )

      CALL ComputeConformalFactor_TwoMoment &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

      CALL ComputeLapseShiftCurvature_TwoMoment &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

      CALL ComputeConserved_Euler_MF( MF_uGF, MF_uPF, MF_uAF, MF_uCF )

      DO iLevel = 0, nLevels - 1

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

      IF( amrex_parallel_ioprocessor() )THEN

        IF( ITER .GT. MAX_ITER - 3 ) &
          WRITE(*,'(4x,A,I2.2,1x,ES24.16E3)') &
            'ITER, MAX( MaxLF, MaxCF ): ', ITER, MAX( MaxLF, MaxCF )

      END IF

      IF( MAX( MaxLF, MaxCF ) .LT. TOLERANCE ) CONVERGED = .TRUE.

      IF( ITER .EQ. MAX_ITER ) &
        CALL DescribeError_MF &
               ( 903, Real_Option = [ MAX( MaxLF, MaxCF ) ] )

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

  END SUBROUTINE InitializeMetric_TwoMoment_MF_Poseidon


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE ComputeConformalFactor_Euler( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

    CALL ComputeConformalFactor_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

    CALL UpdateConformalFactorAndMetric_XCFC_MF( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

#endif

  END SUBROUTINE ComputeConformalFactor_Euler


  SUBROUTINE ComputeLapseShiftCurvature_Euler( MF_uGF, MF_uCF, MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    CALL ComputePressureTensorTrace_XCFC_Euler_MF &
           ( MF_uGF, MF_uCF, MF_uGS )

    CALL ComputeLapseShiftCurvature_XCFC_MF_Poseidon &
           ( MF_uGS, MF_uMF )

    CALL UpdateLapseShiftCurvature_XCFC_MF( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

    CALL ApplyBoundaryConditions_Geometry_XCFC_MF( MF_uGF )

#endif

  END SUBROUTINE ComputeLapseShiftCurvature_Euler


  SUBROUTINE ComputeConformalFactor_TwoMoment &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC
#ifndef THORNADO_NOTRANSPORT

    CALL ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    CALL ComputeConformalFactor_XCFC_MF_Poseidon( MF_uGS, MF_uMF )

    CALL UpdateConformalFactorAndMetric_XCFC_MF( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

#endif
#endif

  END SUBROUTINE ComputeConformalFactor_TwoMoment


  SUBROUTINE ComputeLapseShiftCurvature_TwoMoment &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS, MF_uMF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uMF(0:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC
#ifndef THORNADO_NOTRANSPORT

    CALL ComputePressureTensorTrace_XCFC_TwoMoment_MF &
           ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    CALL ComputeLapseShiftCurvature_XCFC_MF_Poseidon &
           ( MF_uGS, MF_uMF )

    CALL AverageDown( MF_uGF )

    CALL ApplyBoundaryConditions_Geometry_XCFC_MF( MF_uGF )

#endif
#endif

  END SUBROUTINE ComputeLapseShiftCurvature_TwoMoment


END MODULE MF_GravitySolutionModule_XCFC_Poseidon
