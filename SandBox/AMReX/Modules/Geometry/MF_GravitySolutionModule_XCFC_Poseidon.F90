MODULE MF_GravitySolutionModule_XCFC_Poseidon

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_parmparse_module, ONLY: &
    amrex_parmparse, &
    amrex_parmparse_build, &
    amrex_parmparse_destroy
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_max

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX, &
    nDOFZ, &
    iE_E0, &
    iE_B0, &
    nDOFE, &
    nE
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE
  USE GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Psi
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nCF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    nPF
  USE GeometryFieldsModuleE, ONLY: &
    iGE_Ep3
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    nCR,    &
    iCR_N,  &
    iCR_G1, &
    iCR_G2, &
    iCR_G3
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler
  USE ReferenceElementModuleE, ONLY: &
    WeightsE

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three, &
    FourPi
  USE MaskModule, ONLY: &
    CreateFineMask, &
    DestroyFineMask, &
    IsNotLeafElement
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler_MF
  USE MF_ErrorModule, ONLY: &
    DescribeError_MF
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    xR, &
    swX, &
    iRestart
  USE AverageDownModule, ONLY: &
    AverageDown
  USE XCFC_UtilitiesModule, ONLY: &
    iGS_E, &
    iGS_S1, &
    iGS_S2, &
    iGS_S3, &
    iGS_S, &
    iGS_Mg, &
    nGS, &
    nMF, &
    swXX, &
    MultiplyWithPsi6_MF, &
    UpdateConformalFactorAndMetric_MF, &
    UpdateGeometry_MF, &
    ApplyBoundaryConditions_Geometry, &
    ComputeGravitationalMass_MF, &
    PopulateMF_uMF

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

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
  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon
  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF_Poseidon
  PUBLIC :: ComputeConformalFactor_MF_Poseidon
  PUBLIC :: ComputeGeometry_MF_Poseidon
  PUBLIC :: ComputePressureTensorTrace_XCFC_Euler_MF_Poseidon
  PUBLIC :: ComputePressureTensorTrace_XCFC_TwoMoment_MF_Poseidon
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
    REAL(DP)         :: Psi_xR, AlphaPsi_xR, Beta_u_xR(3)
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

      CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_uCF_tmp )

      CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon &
             ( MF_uGF, MF_uCF_tmp, MF_uGS )

      CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_uCF_tmp )

      CALL Poseidon_Input_Sources_Part1( MF_uGS, nGS )

      ! --- Set Boundary Values ---

      CALL SetPoseidonBoundaryConditions_Outer &
             ( MF_uGS, Psi_xR, AlphaPsi_xR, Beta_u_xR )

      INNER_BC_TYPES = [ 'N', 'N', 'N', 'N', 'N' ] ! Neumann
      OUTER_BC_TYPES = [ 'D', 'D', 'D', 'D', 'D' ] ! Dirichlet

      INNER_BC_VALUES = [ Zero  , Zero       , Zero, Zero, Zero ]
      OUTER_BC_VALUES = [ Psi_xR, AlphaPsi_xR, &
                          Beta_u_xR(1), Beta_u_xR(2), Beta_u_xR(3) ]

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


  SUBROUTINE ComputeConformalFactor_MF_Poseidon( MF_uGS, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:nLevels-1) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    REAL(DP)         :: Psi_xR, AlphaPsi_xR, Beta_u_xR(3)
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1) ! Metric Fields

    INTEGER :: iLevel

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nMF, swXX )
      CALL MF_uMF(iLevel) % SetVal( Zero )

    END DO

    ! --- Set Boundary Values ---

    CALL SetPoseidonBoundaryConditions_Outer &
           ( MF_uGS, Psi_xR, AlphaPsi_xR, Beta_u_xR )

    INNER_BC_TYPES = [ 'N', 'N', 'N', 'N', 'N' ] ! Neumann
    OUTER_BC_TYPES = [ 'D', 'D', 'D', 'D', 'D' ] ! Dirichlet

    INNER_BC_VALUES = [ Zero  , Zero       , Zero, Zero, Zero ]
    OUTER_BC_VALUES = [ Psi_xR, AlphaPsi_xR, &
                        Beta_u_xR(1), Beta_u_xR(2), Beta_u_xR(3) ]

    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( 'I', INNER_BC_TYPES, INNER_BC_VALUES )
    CALL Poseidon_Set_Uniform_Boundary_Conditions &
           ( 'O', OUTER_BC_TYPES, OUTER_BC_VALUES)

    ! --- Set XCFC sources with current conformal factor ---
    CALL Poseidon_Input_Sources_Part1( MF_uGS, nGS )

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

    CALL Poseidon_Return_Conformal_Factor &
           ( MF_uMF, FillGhostCells_Option = FillGhostCells )

    CALL UpdateConformalFactorAndMetric_MF( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uMF(iLevel) )

    END DO

#endif

  END SUBROUTINE ComputeConformalFactor_MF_Poseidon


  SUBROUTINE ComputeGeometry_MF_Poseidon( MF_uGS, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:nLevels-1) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1) ! Metric Fields

    INTEGER :: iLevel

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nMF, swXX )
      CALL MF_uMF(iLevel) % SetVal( Zero )

    END DO

    ! --- Set gravity sources with updated conformal factor ---

    CALL Poseidon_Input_Sources_Part2( MF_uGS, nGS )

    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()

    CALL Poseidon_Return_ALL &
           ( MF_uMF, FillGhostCells_Option = FillGhostCells )

    ! --- Copy data from Poseidon to thornado ---

    CALL UpdateGeometry_MF( MF_uMF, MF_uGF )

    CALL AverageDown( MF_uGF )

    CALL ApplyBoundaryConditions_Geometry( MF_uGF )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uMF(iLevel) )

    END DO

#endif

  END SUBROUTINE ComputeGeometry_MF_Poseidon


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:)

    TYPE(amrex_imultifab) :: iMF_FineMask
    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS     (:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, ErrorExists
    INTEGER  :: iX_B0(3), iX_E0(3)
    REAL(DP) :: Psi6
    REAL(DP) :: uCF_K(nCF), uPF_K(nPF), &
                LorentzFactor, BetaDotV, Enthalpy, Pressure

    INTEGER, ALLOCATABLE :: ITERATION(:,:,:,:)
    INTEGER, ALLOCATABLE :: iErr     (:,:,:,:)

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uCF      => MF_uCF(iLevel) % DataPtr( MFI )
        uGS      => MF_uGS(iLevel) % DataPtr( MFI )
        FineMask => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ErrorExists = 0

        ALLOCATE( ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )
        ALLOCATE( iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          ITERATION(iNX,iX1,iX2,iX3) = 0
          iErr     (iNX,iX1,iX2,iX3) = 0

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          uGS       (iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) &
            =  ( uCF(iX1,iX2,iX3,nDOFX*(iCF_E-1)+iNX) &
               + uCF(iX1,iX2,iX3,nDOFX*(iCF_D-1)+iNX) )

          uGS       (iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) &
            = uCF   (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX)

          uGS       (iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) &
            = uCF   (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX)

          uGS       (iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) &
            = uCF   (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX)

          ! Assume Psi^(iStage) ~ Psi^(iStage+1) for Poseidon BCs

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          uCF_K(iCF_D ) = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6
          uCF_K(iCF_S1) = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6
          uCF_K(iCF_S2) = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6
          uCF_K(iCF_S3) = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6
          uCF_K(iCF_E ) = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6
          uCF_K(iCF_Ne) = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6

          CALL ComputePrimitive_Euler_Relativistic &
                 ( uCF_K(iCF_D ), &
                   uCF_K(iCF_S1), &
                   uCF_K(iCF_S2), &
                   uCF_K(iCF_S3), &
                   uCF_K(iCF_E ), &
                   uCF_K(iCF_Ne), &
                   uPF_K(iPF_D ), &
                   uPF_K(iPF_V1), &
                   uPF_K(iPF_V2), &
                   uPF_K(iPF_V3), &
                   uPF_K(iPF_E ), &
                   uPF_K(iPF_Ne), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
                   iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

          ErrorExists = ErrorExists + iErr(iNX,iX1,iX2,iX3)

          CALL ComputePressureFromPrimitive &
                 ( uPF_K(iPF_D), uPF_K(iPF_E), uPF_K(iPF_Ne), Pressure )

          LorentzFactor &
            = One / SQRT( One                              &
                - ( uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                      * uPF_K(iPF_V1)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                      * uPF_K(iPF_V2)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                      * uPF_K(iPF_V3)**2 ) )

          BetaDotV =   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX) &
                         * uPF_K(iPF_V1) &
                     + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX) &
                         * uPF_K(iPF_V2) &
                     + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX) &
                         * uPF_K(iPF_V3)

          Enthalpy = uPF_K(iPF_D) + uPF_K(iPF_E) + Pressure

          uGS(iX1,iX2,iX3,nDOFX*(iGS_Mg-1)+iNX) &
            = ( Enthalpy * ( Two * LorentzFactor**2 &
                  * ( One - BetaDotV &
                              / uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) ) &
                      - One ) &
                + Two * Pressure ) &
               * uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha -1)+iNX) &
               * uGF(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX)

        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .GT. 0 )THEN

          WRITE(*,*) 'ERROR'
          WRITE(*,*) '-----'
          WRITE(*,*) '    MODULE: MF_GravitySolutionModule_XCFC_Poseidon'
          WRITE(*,*) 'SUBROUTINE: ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon'

          CALL CreateMesh_MF( iLevel, MeshX )

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

            CALL DescribeError_Euler &
                   ( iErr(iNX,iX1,iX2,iX3), &
                     Int_Option &
                       = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                           iX_B0(1), iX_B0(2), iX_B0(3), &
                           iX_E0(1), iX_E0(2), iX_E0(3), &
                           iNX, iX1, iX2, iX3 ], &
                     Real_Option &
                       = [ MeshX(1) % Center(iX1), &
                           MeshX(2) % Center(iX2), &
                           MeshX(3) % Center(iX3), &
                           MeshX(1) % Width (iX1), &
                           MeshX(2) % Width (iX2), &
                           MeshX(3) % Width (iX3), &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ], &
                     Char_Option = [ 'NA' ], &
                     Message_Option &
                       = 'Calling from ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon' )

          END DO
          END DO
          END DO
          END DO

          CALL DestroyMesh_MF( MeshX )

        END IF

        DEALLOCATE( ITERATION )
        DEALLOCATE( iErr      )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

#endif

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon


  SUBROUTINE ComputePressureTensorTrace_XCFC_Euler_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_imultifab) :: iMF_FineMask
    TYPE(amrex_box)       :: BX
    TYPE(amrex_mfiter)    :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF     (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS     (:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: FineMask(:,:,:,:)

    INTEGER :: iLevel, iNX, iX1, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: ErrorExists

    INTEGER, ALLOCATABLE :: ITERATION(:,:,:,:)
    INTEGER, ALLOCATABLE :: iErr     (:,:,:,:)

    REAL(DP) :: uPF(nPF), Pressure, Psi6

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    DO iLevel = 0, nLevels-1

      CALL CreateFineMask( iLevel, iMF_FineMask, MF_uGF % BA, MF_uGF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF      => MF_uGF(iLevel) % DataPtr( MFI )
        uCF      => MF_uCF(iLevel) % DataPtr( MFI )
        uGS      => MF_uGS(iLevel) % DataPtr( MFI )
        FineMask => iMF_FineMask   % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ErrorExists = 0

        ALLOCATE( ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )
        ALLOCATE( iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          ITERATION(iNX,iX1,iX2,iX3) = 0
          iErr     (iNX,iX1,iX2,iX3) = 0

          IF( IsNotLeafElement( FineMask(iX1,iX2,iX3,1) ) ) CYCLE

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          ! --- Compute trace of stress tensor ---

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6

          CALL ComputePrimitive_Euler_Relativistic &
                 ( uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX), &
                   uPF(iPF_D ), &
                   uPF(iPF_V1), &
                   uPF(iPF_V2), &
                   uPF(iPF_V3), &
                   uPF(iPF_E ), &
                   uPF(iPF_Ne), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
                   iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) * Psi6

          ErrorExists = ErrorExists + iErr(iNX,iX1,iX2,iX3)

          CALL ComputePressureFromPrimitive &
                 ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) &
            =   (  uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6 * uPF(iPF_V1) &
                 + uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6 * uPF(iPF_V2) &
                 + uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6 * uPF(iPF_V3) &
                 + Three * Pressure ) * Psi6

        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .GT. 0 )THEN

          WRITE(*,*) 'ERROR'
          WRITE(*,*) '-----'
          WRITE(*,*) '    MODULE: MF_GravitySolutionModule_XCFC_Poseidon'
          WRITE(*,*) 'SUBROUTINE: ComputePressureTensorTrace_XCFC_Euler_MF_Poseidon'

          CALL CreateMesh_MF( iLevel, MeshX )

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

            CALL DescribeError_Euler &
                   ( iErr(iNX,iX1,iX2,iX3), &
                     Int_Option &
                       = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                           iX_B0(1), iX_B0(2), iX_B0(3), &
                           iX_E0(1), iX_E0(2), iX_E0(3), &
                           iNX, iX1, iX2, iX3 ], &
                     Real_Option &
                       = [ MeshX(1) % Center(iX1), &
                           MeshX(2) % Center(iX2), &
                           MeshX(3) % Center(iX3), &
                           MeshX(1) % Width (iX1), &
                           MeshX(2) % Width (iX2), &
                           MeshX(3) % Width (iX3), &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ], &
                     Char_Option = [ 'NA' ], &
                     Message_Option &
                       = 'Calling from ComputePressureTensorTrace_XCFC_Euler_MF_Poseidon' )

          END DO
          END DO
          END DO
          END DO

          CALL DestroyMesh_MF( MeshX )

        END IF

        DEALLOCATE( ITERATION )
        DEALLOCATE( iErr      )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

      CALL DestroyFineMask( iMF_FineMask )

    END DO ! iLevel = 0, nLevels-1

#endif

  END SUBROUTINE ComputePressureTensorTrace_XCFC_Euler_MF_Poseidon


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, ErrorExists
    INTEGER  :: iX_B0(3), iX_E0(3)
    REAL(DP) :: Psi6
    REAL(DP) :: uPF(nPF), LorentzFactor, BetaDotV, Enthalpy, Pressure

    INTEGER, ALLOCATABLE :: ITERATION(:,:,:,:)
    INTEGER, ALLOCATABLE :: iErr     (:,:,:,:)

    REAL(DP) :: E, S_i(3), E_int, S_i_int(3)
    REAL(DP) :: N, G_d_1, G_d_2, G_d_3, vG
    REAL(DP) :: V_u_1, V_u_2, V_u_3
    REAL(DP) :: V_d_1, V_d_2, V_d_3
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    INTEGER  :: iD_N, iD_G1, iD_G2, iD_G3, iE, iN_E, iS, iN_Z

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ASSOCIATE &
      ( dZ1 => MeshE  % Width )

    E       = Zero
    E_int   = Zero
    S_i     = Zero
    S_i_int = Zero

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ErrorExists = 0

        ALLOCATE( ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )
        ALLOCATE( iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          uGS       (iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) &
            =  ( uCF(iX1,iX2,iX3,nDOFX*(iCF_E-1)+iNX) &
               + uCF(iX1,iX2,iX3,nDOFX*(iCF_D-1)+iNX) )

          uGS    (iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX)

          uGS    (iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX)

          uGS    (iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX)

          ! Assume Psi^(iStage) ~ Psi^(iStage+1) for Poseidon BCs

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          ITERATION(iNX,iX1,iX2,iX3) = 0
          iErr     (iNX,iX1,iX2,iX3) = 0

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6

          CALL ComputePrimitive_Euler_Relativistic &
                 ( uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX), &
                   uPF(iPF_D ), &
                   uPF(iPF_V1), &
                   uPF(iPF_V2), &
                   uPF(iPF_V3), &
                   uPF(iPF_E ), &
                   uPF(iPF_Ne), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
                   iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) * Psi6

          ErrorExists = ErrorExists + iErr(iNX,iX1,iX2,iX3)

          CALL ComputePressureFromPrimitive &
                 ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

          LorentzFactor &
            = One / SQRT( One                              &
                - ( uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                      * uPF(iPF_V1)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                      * uPF(iPF_V2)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                      * uPF(iPF_V3)**2 ) )

          BetaDotV =   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX) &
                         * uPF(iPF_V1) &
                     + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX) &
                         * uPF(iPF_V2) &
                     + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                         * uGF(iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX) &
                         * uPF(iPF_V3)

          Enthalpy = uPF(iPF_D) + uPF(iPF_E) + Pressure

          uGS(iX1,iX2,iX3,nDOFX*(iGS_Mg-1)+iNX) &
            = ( Enthalpy * ( Two * LorentzFactor**2 &
                  * ( One - BetaDotV &
                              / uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) ) &
                      - One ) &
                + Two * Pressure ) &
               * uGF(iX1,iX2,iX3,nDOFX*(iGF_Alpha -1)+iNX) &
               * uGF(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX)

          DO iS   = 1, nSpecies
          DO iE   = 1, nE
          DO iN_E = 1, nDOFE

            iN_Z = (iNX-1) * nDOFE + iN_E

            iD_N = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iCR_N - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G1 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iCR_G1 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G2 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iCR_G2 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G3 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iCR_G3 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                 + ( iE - 1 ) * nDOFZ + iN_Z

            N     = uCR(iX1,iX2,iX3,iD_N)
            G_d_1 = uCR(iX1,iX2,iX3,iD_G1)
            G_d_2 = uCR(iX1,iX2,iX3,iD_G2)
            G_d_3 = uCR(iX1,iX2,iX3,iD_G3)

            V_u_1 = uPF (iPF_V1)
            V_u_2 = uPF (iPF_V2)
            V_u_3 = uPF (iPF_V3)

            Gm_dd_11 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX)
            Gm_dd_22 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX)
            Gm_dd_33 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX)

            V_d_1 = V_u_1 * Gm_dd_11
            V_d_2 = V_u_2 * Gm_dd_22
            V_d_3 = V_u_3 * Gm_dd_33

            vG = V_u_1 * G_d_1 + V_u_2 * G_d_2 + V_u_3 * G_d_3

            E_int      = LorentzFactor * N + vG
            S_i_int(1) = LorentzFactor * V_d_1 * N + G_d_1
            S_i_int(2) = LorentzFactor * V_d_2 * N + G_d_2
            S_i_int(3) = LorentzFactor * V_d_3 * N + G_d_3

            E = E &
              + FourPi * dZ1(iE) * WeightsE(iN_E) &
              * uGE(iN_E,iE,iGE_Ep3) * E_int
            S_i(1) = S_i(1) &
                   + FourPi * dZ1(iE) * WeightsE(iN_E) &
                   * uGE(iN_E,iE,iGE_Ep3) * S_i_int(1)
            S_i(2) = S_i(2) &
                   + FourPi * dZ1(iE) * WeightsE(iN_E) &
                   * uGE(iN_E,iE,iGE_Ep3) * S_i_int(2)
            S_i(3) = S_i(3) &
                   + FourPi * dZ1(iE) * WeightsE(iN_E) &
                   * uGE(iN_E,iE,iGE_Ep3) * S_i_int(3)

          END DO
          END DO
          END DO

          uGS(iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) + E

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) + S_i(1)

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) + S_i(2)

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) + S_i(3)

          E   = Zero
          S_i = Zero

        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .GT. 0 )THEN

          WRITE(*,*) 'ERROR'
          WRITE(*,*) '-----'
          WRITE(*,*) &
            '    MODULE: MF_GravitySolutionModule_XCFC_Poseidon'
          WRITE(*,*) &
            'SUBROUTINE: ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF_Poseidon'

          CALL CreateMesh_MF( iLevel, MeshX )

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

            CALL DescribeError_Euler &
                   ( iErr(iNX,iX1,iX2,iX3), &
                     Int_Option &
                       = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                           iX_B0(1), iX_B0(2), iX_B0(3), &
                           iX_E0(1), iX_E0(2), iX_E0(3), &
                           iNX, iX1, iX2, iX3 ], &
                     Real_Option &
                       = [ MeshX(1) % Center(iX1), &
                           MeshX(2) % Center(iX2), &
                           MeshX(3) % Center(iX3), &
                           MeshX(1) % Width (iX1), &
                           MeshX(2) % Width (iX2), &
                           MeshX(3) % Width (iX3), &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ], &
                     Char_Option = [ 'NA' ], &
                     Message_Option &
                       = 'Calling from ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF_Poseidon' )

          END DO
          END DO
          END DO
          END DO

          CALL DestroyMesh_MF( MeshX )

        END IF

        DEALLOCATE( ITERATION )
        DEALLOCATE( iErr      )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

    END DO ! iLevel = 0, nLevels-1

    END ASSOCIATE

#endif

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF_Poseidon


  SUBROUTINE ComputePressureTensorTrace_XCFC_TwoMoment_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCR(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCR(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)

    INTEGER :: iLevel, iNX, iX1, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: ErrorExists
    INTEGER, ALLOCATABLE :: ITERATION(:,:,:,:)
    INTEGER, ALLOCATABLE :: iErr     (:,:,:,:)

    REAL(DP) :: uPF(nPF), Pressure, Psi6

    REAL(DP) :: S, S_int, N, G_d_1, G_d_2, G_d_3, vG
    REAL(DP) :: LorentzFactor, V_u_1, V_u_2, V_u_3
    INTEGER  :: iD_N, iD_G1, iD_G2, iD_G3, iE, iN_E, iS, iN_Z

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    ASSOCIATE &
      ( dZ1 => MeshE  % Width )

    S = Zero

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uCR => MF_uCR(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ALLOCATE( ITERATION(1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )
        ALLOCATE( iErr     (1:nDOFX,iX_B0(1):iX_E0(1), &
                                    iX_B0(2):iX_E0(2), &
                                    iX_B0(3):iX_E0(3)) )

        ErrorExists = 0

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          ITERATION(iNX,iX1,iX2,iX3) = 0
          iErr     (iNX,iX1,iX2,iX3) = 0

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          ! --- Compute trace of stress tensor ---

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6

          CALL ComputePrimitive_Euler_Relativistic &
                 ( uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX), &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX), &
                   uPF(iPF_D ), &
                   uPF(iPF_V1), &
                   uPF(iPF_V2), &
                   uPF(iPF_V3), &
                   uPF(iPF_E ), &
                   uPF(iPF_Ne), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   ITERATION_Option = ITERATION(iNX,iX1,iX2,iX3), &
                   iErr_Option      = iErr     (iNX,iX1,iX2,iX3) )

          uCF    (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) * Psi6
          uCF    (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) &
            = uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) * Psi6

          ErrorExists = ErrorExists + iErr(iNX,iX1,iX2,iX3)

          CALL ComputePressureFromPrimitive &
                 ( uPF(iPF_D), uPF(iPF_E), uPF(iPF_Ne), Pressure )

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) &
            =   (  uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6 * uPF(iPF_V1) &
                 + uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6 * uPF(iPF_V2) &
                 + uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6 * uPF(iPF_V3) &
                 + Three * Pressure ) * Psi6

          LorentzFactor &
            = One / SQRT( One                              &
                - ( uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                      * uPF(iPF_V1)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                      * uPF(iPF_V2)**2 &
                  + uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                      * uPF(iPF_V3)**2 ) )

          DO iS   = 1, nSpecies
          DO iE   = 1, nE
          DO iN_E = 1, nDOFE

            iN_Z = ( iNX - 1 ) * nDOFE + iN_E

            iD_N  = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iCR_N - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G1 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iCR_G1 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G2 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iCR_G2 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iE - 1 ) * nDOFZ + iN_Z
            iD_G3 = ( iS - 1 ) * nCR * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iCR_G3 - 1 ) * ( iE_E0 - iE_B0 + 1 ) * nDOFZ &
                      + ( iE - 1 ) * nDOFZ + iN_Z

            N     = uCR(iX1,iX2,iX3,iD_N )
            G_d_1 = uCR(iX1,iX2,iX3,iD_G1)
            G_d_2 = uCR(iX1,iX2,iX3,iD_G2)
            G_d_3 = uCR(iX1,iX2,iX3,iD_G3)

            V_u_1 = uPF (iPF_V1)
            V_u_2 = uPF (iPF_V2)
            V_u_3 = uPF (iPF_V3)

            vG = V_u_1 * G_d_1 + V_u_2 * G_d_2 + V_u_3 * G_d_3

            S_int = LorentzFactor * N + vG

            S &
              = S + FourPi * dZ1(iE) * WeightsE(iN_E) &
                      * uGE(iN_E,iE,iGE_Ep3) * S_int

          END DO
          END DO
          END DO

          uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) &
            = uGS(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) + S

          S = Zero

        END DO
        END DO
        END DO
        END DO

        IF( ErrorExists .NE. 0 )THEN

          WRITE(*,*) &
            'ERROR'
          WRITE(*,*) &
            '-----'
          WRITE(*,*) &
            '    MODULE: Poseidon_UtilitiesModule'
          WRITE(*,*) &
            'SUBROUTINE: ComputePressureTensorTrace_XCFC_TwoMoment_MF_Poseidon'

          CALL CreateMesh_MF( iLevel, MeshX )

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

            CALL DescribeError_Euler &
                   ( iErr(iNX,iX1,iX2,iX3), &
                     Int_Option &
                       = [ ITERATION(iNX,iX1,iX2,iX3), 99999999, &
                           iX_B0(1), iX_B0(2), iX_B0(3), &
                           iX_E0(1), iX_E0(2), iX_E0(3), &
                           iNX, iX1, iX2, iX3 ], &
                     Real_Option &
                       = [ MeshX(1) % Center(iX1), &
                           MeshX(2) % Center(iX2), &
                           MeshX(3) % Center(iX3), &
                           MeshX(1) % Width (iX1), &
                           MeshX(2) % Width (iX2), &
                           MeshX(3) % Width (iX3), &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                           uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                           uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ], &
                     Char_Option = [ 'NA' ], &
                     Message_Option &
                       = 'Calling from ComputePressureTensorTrace_XCFC_TwoMoment_MF_Poseidon' )

          END DO
          END DO
          END DO
          END DO

          CALL DestroyMesh_MF( MeshX )

        END IF

        DEALLOCATE( ITERATION )
        DEALLOCATE( iErr      )

      END DO ! WHILE( MFI % next() )

      CALL amrex_mfiter_destroy( MFI )

    END DO ! iLevel = 0, nLevels-1

    END ASSOCIATE

#endif

  END SUBROUTINE ComputePressureTensorTrace_XCFC_TwoMoment_MF_Poseidon


  SUBROUTINE InitializeMetric_Euler_MF_Poseidon &
    ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uPF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uAF(0:nLevels-1)

    TYPE(amrex_multifab) :: MF_uGS(0:nLevels-1)
    TYPE(amrex_multifab) :: LF1   (0:nLevels-1)
    TYPE(amrex_multifab) :: LF2   (0:nLevels-1)
    TYPE(amrex_multifab) :: dLF   (0:nLevels-1)
    TYPE(amrex_multifab) :: CF1   (0:nLevels-1)
    TYPE(amrex_multifab) :: CF2   (0:nLevels-1)
    TYPE(amrex_multifab) :: dCF   (0:nLevels-1)

    LOGICAL  :: CONVERGED
    INTEGER  :: iLevel, ITER, iNX
    REAL(DP) :: MinLF, MinCF

    INTEGER, PARAMETER :: MAX_ITER = 10

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nGS, swXX )
      CALL MF_uGS(iLevel) % SetVal( Zero )

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

    MinLF = HUGE( One )
    MinCF = HUGE( One )

    DO WHILE( .NOT. CONVERGED )

      ITER = ITER + 1

      DO iLevel = 0, nLevels - 1

        CALL LF1(iLevel) % COPY &
              ( MF_uGF(iLevel), nDOFX*(iGF_Alpha-1)+1, 1, nDOFX, 0 )

        CALL CF1(iLevel) % COPY &
              ( MF_uGF(iLevel), nDOFX*(iGF_Psi  -1)+1, 1, nDOFX, 0 )

      END DO

      CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_uCF )

      CALL ComputeConformalFactorSourcesAndMg_XCFC_Euler_MF_Poseidon &
             ( MF_uGF, MF_uCF, MF_uGS )

      CALL ComputeConformalFactor_MF_Poseidon( MF_uGS, MF_uGF )

      CALL ComputePressureTensorTrace_XCFC_Euler_MF_Poseidon &
             ( MF_uGF, MF_uCF, MF_uGS )

      CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_uCF )

      CALL ComputeGeometry_MF_Poseidon( MF_uGS, MF_uGF )

      CALL ComputeConserved_Euler_MF( MF_uGF, MF_uPF, MF_uAF, MF_uCF )

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

          MinLF = MIN( MinLF, dLF(iLevel) % Norm0( iNX ) )
          MinCF = MIN( MinCF, dCF(iLevel) % Norm0( iNX ) )

        END DO

      END DO ! iLevel = 0, nLevels-1

      CALL amrex_parallel_reduce_max( MinLF )
      CALL amrex_parallel_reduce_max( MinCF )

      IF( MAX( MinLF, MinCF ) .LT. 1.0e-13_DP ) CONVERGED = .TRUE.

      IF( ITER .EQ. MAX_ITER ) &
        CALL DescribeError_MF &
               ( 901, Real_Option = [ MAX( MinLF, MinCF ) ] )

    END DO ! WHILE( .NOT. CONVERGED )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( dCF   (iLevel) )
      CALL amrex_multifab_destroy( CF2   (iLevel) )
      CALL amrex_multifab_destroy( CF1   (iLevel) )
      CALL amrex_multifab_destroy( dLF   (iLevel) )
      CALL amrex_multifab_destroy( LF2   (iLevel) )
      CALL amrex_multifab_destroy( LF1   (iLevel) )
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
    TYPE(amrex_multifab) :: LF1   (0:nLevels-1)
    TYPE(amrex_multifab) :: LF2   (0:nLevels-1)
    TYPE(amrex_multifab) :: dLF   (0:nLevels-1)
    TYPE(amrex_multifab) :: CF1   (0:nLevels-1)
    TYPE(amrex_multifab) :: CF2   (0:nLevels-1)
    TYPE(amrex_multifab) :: dCF   (0:nLevels-1)

    LOGICAL  :: CONVERGED
    INTEGER  :: iLevel, ITER, iNX
    REAL(DP) :: MinLF, MinCF

    INTEGER, PARAMETER :: MAX_ITER = 10

#ifdef GRAVITY_SOLVER_POSEIDON_XCFC

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nGS, 0 )
      CALL MF_uGS(iLevel) % SetVal( Zero )

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

    MinLF = HUGE( One )
    MinCF = HUGE( One )

    DO WHILE( .NOT. CONVERGED )

      ITER = ITER + 1

      DO iLevel = 0, nLevels - 1

        CALL LF1(iLevel) % COPY &
              ( MF_uGF(iLevel), nDOFX*(iGF_Alpha-1)+1, 1, nDOFX, 0 )

        CALL CF1(iLevel) % COPY &
              ( MF_uGF(iLevel), nDOFX*(iGF_Psi  -1)+1, 1, nDOFX, 0 )

      END DO

      CALL MultiplyWithPsi6_MF( MF_uGF, +1, 1, 1, 1, 1, MF_uCF )

      CALL ComputeConformalFactorSourcesAndMg_XCFC_TwoMoment_MF_Poseidon &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

      CALL ComputeConformalFactor_MF_Poseidon( MF_uGS, MF_uGF )

      CALL ComputePressureTensorTrace_XCFC_TwoMoment_MF_Poseidon &
             ( MF_uGF, MF_uCF, MF_uCR, MF_uGS )

      CALL MultiplyWithPsi6_MF( MF_uGF, -1, 1, 1, 1, 1, MF_uCF )

      CALL ComputeGeometry_MF_Poseidon( MF_uGS, MF_uGF )

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

          MinLF = MIN( MinLF, dLF(iLevel) % Norm0( iNX ) )
          MinCF = MIN( MinCF, dCF(iLevel) % Norm0( iNX ) )

        END DO

      END DO ! iLevel = 0, nLevels-1

      CALL amrex_parallel_reduce_max( MinLF )
      CALL amrex_parallel_reduce_max( MinCF )

      IF( MAX( MinLF, MinCF ) .LT. 1.0e-13_DP ) CONVERGED = .TRUE.

      IF( ITER .EQ. MAX_ITER ) &
        CALL DescribeError_MF &
               ( 901, Real_Option = [ MAX( MinLF, MinCF ) ] )

    END DO ! WHILE( .NOT. CONVERGED )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( dCF   (iLevel) )
      CALL amrex_multifab_destroy( CF2   (iLevel) )
      CALL amrex_multifab_destroy( CF1   (iLevel) )
      CALL amrex_multifab_destroy( dLF   (iLevel) )
      CALL amrex_multifab_destroy( LF2   (iLevel) )
      CALL amrex_multifab_destroy( LF1   (iLevel) )
      CALL amrex_multifab_destroy( MF_uGS(iLevel) )

    END DO

#endif

  END SUBROUTINE InitializeMetric_TwoMoment_MF_Poseidon


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE SetPoseidonBoundaryConditions_Outer &
    ( MF_uGS, Psi_xR, AlphaPsi_xR, Beta_u_xR )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:nLevels-1)
    REAL(DP)            , INTENT(out) :: Psi_xR, AlphaPsi_xR, Beta_u_xR(3)

    REAL(DP) :: GravitationalMass
    REAL(DP) :: Alpha_xR

    CALL ComputeGravitationalMass_MF( MF_uGS, GravitationalMass )

    ! --- Approximate outer boundary with isotropic expressions ---

    Psi_xR   = One + Half * Gravitationalmass / xR(1)
    Alpha_xR =   ( One - Half * GravitationalMass / xR(1) ) &
               / ( One + Half * GravitationalMass / xR(1) )

    AlphaPsi_xR = Alpha_xR * Psi_xR
    Beta_u_xR   = Zero

  END SUBROUTINE SetPoseidonBoundaryConditions_Outer


END MODULE MF_GravitySolutionModule_XCFC_Poseidon
