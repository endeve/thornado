MODULE MF_GravitySolutionModule_XCFC_Poseidon

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_multifab_build, &
    amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy, &
    amrex_imultifab
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_max, &
    amrex_parallel_reduce_sum

  ! --- thornado Modules ---

  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nNodesX
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE LinearAlgebraModule, ONLY: &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX, &
    nDOFX_X1, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Up
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_Psi, &
    iGF_K_dd_11, &
    iGF_K_dd_12, &
    iGF_K_dd_13, &
    iGF_K_dd_22, &
    iGF_K_dd_23, &
    iGF_K_dd_33, &
    nGF
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
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Two, &
    Three, &
    Pi, &
    SqrtTiny
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler_MF
  USE MF_Euler_ErrorModule, ONLY: &
    DescribeError_Euler_MF
  USE MakeFineMaskModule, ONLY: &
    MakeFineMask, &
    iLeaf_MFM
  USE InputParsingModule, ONLY: &
    nLevels, &
    UseTiling, &
    MaxGridSizeX, &
    nX, &
    nNodes, &
    xL, &
    xR

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

  ! --- Poseidon Modules ---

  USE Initialization_AMReX, ONLY: &
    Initialize_Poseidon_with_AMReX
  USE Poseidon_Main_Module, ONLY: &
    Poseidon_Close
  USE Poseidon_Interface_BC_Input, ONLY : &
    Poseidon_Set_Uniform_Boundary_Conditions
  USE Poseidon_Source_Input_Module, ONLY: &
    Poseidon_Input_Sources_Part1, &
    Poseidon_Input_Sources_Part2
  USE Poseidon_XCFC_Interface_Module, ONLY: &
    Poseidon_XCFC_Run_Part1, &
    Poseidon_XCFC_Run_Part2
  USE Poseidon_Return_Routines_Module, ONLY: &
    Poseidon_Return_Conformal_Factor, &
    Poseidon_Return_ALL
  USE Poseidon_Initial_Guess_Module, ONLY: &
    Poseidon_Initialize_Flat_Guess

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_XCFC_Poseidon_MF
  PUBLIC :: FinalizeGravitySolver_XCFC_Poseidon_MF
  PUBLIC :: ComputeConformalFactor_Poseidon_MF
  PUBLIC :: ComputeGeometry_Poseidon_MF
  PUBLIC :: ComputeConformalFactorSourcesAndMg_XCFC_MF
  PUBLIC :: ComputePressureTensorTrace_XCFC_MF
  PUBLIC :: MultiplyWithPsi6_MF
  PUBLIC :: InitializeMetric_MF

  ! --- MF: Metric Fields ---

  INTEGER, PARAMETER :: iMF_Psi      = 1
  INTEGER, PARAMETER :: iMF_Alpha    = 2
  INTEGER, PARAMETER :: iMF_Beta_1   = 3
  INTEGER, PARAMETER :: iMF_Beta_2   = 4
  INTEGER, PARAMETER :: iMF_Beta_3   = 5
  INTEGER, PARAMETER :: iMF_K_dd_11  = 6
  INTEGER, PARAMETER :: iMF_K_dd_12  = 7
  INTEGER, PARAMETER :: iMF_K_dd_13  = 8
  INTEGER, PARAMETER :: iMF_K_dd_22  = 9
  INTEGER, PARAMETER :: iMF_K_dd_23  = 10
  INTEGER, PARAMETER :: iMF_K_dd_33  = 11
  INTEGER, PARAMETER :: nMF          = 11

  ! --- GS: Gravity/Geometry Sources ---

  INTEGER, PARAMETER :: iGS_E  = 1
  INTEGER, PARAMETER :: iGS_S1 = 2
  INTEGER, PARAMETER :: iGS_S2 = 3
  INTEGER, PARAMETER :: iGS_S3 = 4
  INTEGER, PARAMETER :: iGS_S  = 5
  INTEGER, PARAMETER :: iGS_Mg = 6
  INTEGER, PARAMETER, PUBLIC :: nGS = 6

CONTAINS


  SUBROUTINE InitializeGravitySolver_XCFC_Poseidon_MF

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    IF( amrex_parallel_ioprocessor() )THEN

      WRITE(*,*)
      WRITE(*,'(4x,A)') &
        'INFO: Gravity Solver (Poseidon, XCFC)'
      WRITE(*,'(4x,A)') &
        '-------------------------------------'
      WRITE(*,*)
      WRITE(*,'(6x,A)') 'Only implemented for 1D spherical symmetry.'
      WRITE(*,*)

    END IF

    CALL Initialize_Poseidon_with_AMReX &
           ( Source_NQ                    = nNodesX, &
             Source_xL                    = [ -Half, +Half ], &
             Source_RQ_xlocs              = MeshX(1) % Nodes, &
             Source_TQ_xlocs              = MeshX(2) % Nodes, &
             Source_PQ_xlocs              = MeshX(3) % Nodes, &
             Source_Units                 = 'G', &
             Source_Radial_Boundary_Units = 'km', &
             Verbose_Option               = .FALSE., &
             Print_Setup_Option           = .TRUE.,  &
             Print_Results_Option         = .FALSE.   )

#endif

  END SUBROUTINE InitializeGravitySolver_XCFC_Poseidon_MF


  SUBROUTINE FinalizeGravitySolver_XCFC_Poseidon_MF

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    CALL Poseidon_Close()

#endif

  END SUBROUTINE FinalizeGravitySolver_XCFC_Poseidon_MF


  SUBROUTINE ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:nLevels-1) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    REAL(DP)         :: Psi_xR, AlphaPsi_xR, Beta_u_xR(3)
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1) ! Metric Fields

    INTEGER :: iLevel

!!$    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nMF, 0 )
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

    CALL Poseidon_Initialize_Flat_Guess() ! Possibly move this to init call

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

    CALL Poseidon_Return_Conformal_Factor( MF_uMF )

    CALL UpdateConformalFactorAndMetric_MF( MF_uMF, MF_uGF )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uMF(iLevel) )

    END DO

#endif

!!$    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeConformalFactor_Poseidon_MF


  SUBROUTINE ComputeGeometry_Poseidon_MF( MF_uGS, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:nLevels-1) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1) ! Metric Fields

    INTEGER :: iLevel
    REAL(DP) :: GravitationalMass

!!$    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uMF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nMF, 0 )
      CALL MF_uMF(iLevel) % SetVal( Zero )

    END DO

    ! --- Set gravity sources with updated conformal factor ---

    CALL Poseidon_Input_Sources_Part2( MF_uGS, nGS )

    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()

    CALL Poseidon_Return_ALL( MF_uMF )

    ! --- Copy data from Poseidon to thornado ---

    CALL ComputeGeometryFromPoseidon_MF( MF_uMF, MF_uGF )

    CALL ComputeGravitationalMass( MF_uGS, GravitationalMass )

    CALL ApplyBoundaryConditions_Geometry( GravitationalMass, MF_uGF )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_uMF(iLevel) )

    END DO

#endif

!!$    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeGeometry_Poseidon_MF


  SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_MF &
    ( MF_uGF, MF_uCF, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, jErr
    INTEGER  :: iX_B0(3), iX_E0(3)
    REAL(DP) :: Psi6
    REAL(DP) :: uPF(nPF), LorentzFactor, BetaDotV, Enthalpy, Pressure

    TYPE(amrex_imultifab) :: iMF_Mask

    INTEGER, ALLOCATABLE :: iErr(:,:,:,:)

    DO iLevel = 0, nLevels-1

      CALL MakeFineMask( iLevel, iMF_Mask, MF_uGF % BA, MF_uGF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        jErr = 0

        ALLOCATE( iErr(1:nDOFX,iX_B0(1):iX_E0(1), &
                               iX_B0(2):iX_E0(2), &
                               iX_B0(3):iX_E0(3)) )

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

          uGS       (iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) &
            =  ( uCF(iX1,iX2,iX3,nDOFX*(iCF_E-1)+iNX) &
               + uCF(iX1,iX2,iX3,nDOFX*(iCF_D-1)+iNX) )

          uGS        (iX1,iX2,iX3,nDOFX*(iGS_S1-1)+iNX) &
            = uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX)

          uGS        (iX1,iX2,iX3,nDOFX*(iGS_S2-1)+iNX) &
            = uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX)

          uGS        (iX1,iX2,iX3,nDOFX*(iGS_S3-1)+iNX) &
            = uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX)

          ! Assume Psi^(iStage) ~ Psi^(iStage+1) for Poseidon BCs

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          iErr(iNX,iX1,iX2,iX3) = 0

          CALL ComputePrimitive_Euler_Relativistic &
                 ( uCF (iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                   uCF (iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                   uCF (iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                   uCF (iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                   uCF (iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                   uCF (iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                   uPF (iPF_D ), &
                   uPF (iPF_V1), &
                   uPF (iPF_V2), &
                   uPF (iPF_V3), &
                   uPF (iPF_E ), &
                   uPF (iPF_Ne), &
                   uGF (iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF (iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF (iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   iErr(iNX,iX1,iX2,iX3) )

          jErr = jErr + iErr(iNX,iX1,iX2,iX3)

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

        END DO
        END DO
        END DO
        END DO

        IF( jErr .GT. 0 )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1       , nDOFX

            IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

            Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

            CALL DescribeError_Euler &
              ( iErr(iNX,iX1,iX2,iX3), &
                Int_Option = [ iNX ], &
                Real_Option = [ uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                                uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                                uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                                uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                                uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                                uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                                uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                                uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                                uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ] )

          END DO
          END DO
          END DO
          END DO

        END IF

        DEALLOCATE( iErr )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE ComputeConformalFactorSourcesAndMg_XCFC_MF


  SUBROUTINE ComputePressureTensorTrace_XCFC_MF( MF_uGF, MF_uCF, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! Psi^6 * U
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)

    INTEGER :: iLevel, iNX, iX1, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: jErr
    INTEGER, ALLOCATABLE :: iErr(:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_Mask

    REAL(DP) :: uPF(nPF), Pressure, Psi6

!    CALL TimersStart_Euler( Timer_GS_ComputeSourceTerms )

    DO iLevel = 0, nLevels-1

      CALL MakeFineMask( iLevel, iMF_Mask, MF_uGF % BA, MF_uGF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ALLOCATE( iErr(1:nDOFX,iX_B0(1):iX_E0(1), &
                               iX_B0(2):iX_E0(2), &
                               iX_B0(3):iX_E0(3)) )

        jErr = 0

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

          iErr(iNX,iX1,iX2,iX3) = 0

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          ! --- Compute trace of stress tensor ---

          CALL ComputePrimitive_Euler_Relativistic &
                 ( uCF(iX1,iX2,iX3,nDOFX*(iCF_D -1)+iNX) / Psi6, &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) / Psi6, &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) / Psi6, &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) / Psi6, &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_E -1)+iNX) / Psi6, &
                   uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne-1)+iNX) / Psi6, &
                   uPF(iPF_D ), &
                   uPF(iPF_V1), &
                   uPF(iPF_V2), &
                   uPF(iPF_V3), &
                   uPF(iPF_E ), &
                   uPF(iPF_Ne), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                   uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX), &
                   iErr(iNX,iX1,iX2,iX3) )

          jErr = jErr + iErr(iNX,iX1,iX2,iX3)

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

        IF( jErr .NE. 0 )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1, nDOFX

            IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

            IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

              WRITE(*,*) 'ERROR'
              WRITE(*,*) '-----'
              WRITE(*,*) '    MODULE: Poseidon_UtilitiesModule'
              WRITE(*,*) 'SUBROUTINE: ComputePressureTensorTrace_Poseidon'
              WRITE(*,*) 'iX_B0: ', iX_B0
              WRITE(*,*) 'iX_E0: ', iX_E0
              WRITE(*,*) 'iX1, iX2, iX3: ', iX1, iX2, iX3

              Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

              CALL DescribeError_Euler &
                ( iErr(iNX,iX1,iX2,iX3), &
                  Int_Option = [ iNX ], &
                  Real_Option &
                    = [ uCF(iX1,iX2,iX3,nDOFX*(iCF_D       -1)+iNX) / Psi6, &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_S1      -1)+iNX) / Psi6, &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_S2      -1)+iNX) / Psi6, &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_S3      -1)+iNX) / Psi6, &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_E       -1)+iNX) / Psi6, &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne      -1)+iNX) / Psi6, &
                        uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX), &
                        uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX), &
                        uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) ] )

            END IF

          END DO
          END DO
          END DO
          END DO

        END IF

        DEALLOCATE( iErr )

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

!    CALL TimersStop_Euler( Timer_GS_ComputeSourceTerms )

  END SUBROUTINE ComputePressureTensorTrace_XCFC_MF


  SUBROUTINE MultiplyWithPsi6_MF( MF_uGF, Power, MF_uCF )

    INTEGER             , INTENT(in)    :: Power
    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, iCF
    INTEGER  :: iX_B0(3), iX_E0(3)

    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)
    TYPE(amrex_imultifab) :: iMF_Mask

    REAL(DP) :: Psi6

    DO iLevel = 0, nLevels-1

      CALL MakeFineMask( iLevel, iMF_Mask, MF_uGF % BA, MF_uGF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

          Psi6 = uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX)**6

          DO iCF = 1, nCF

            uCF(iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) &
              = uCF(iX1,iX2,iX3,nDOFX*(iCF-1)+iNX) * Psi6**( Power )

          END DO

        END DO
        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE MultiplyWithPsi6_MF


  SUBROUTINE InitializeMetric_MF( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

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

    END DO

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

      CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_uCF )

      CALL ComputeConformalFactorSourcesAndMg_XCFC_MF( MF_uGF, MF_uCF, MF_uGS )

      CALL ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

      CALL ComputePressureTensorTrace_XCFC_MF( MF_uGF, MF_uCF, MF_uGS )

      CALL MultiplyWithPsi6_MF( MF_uGF, -1, MF_uCF )

      CALL ComputeGeometry_Poseidon_MF( MF_uGS, MF_uGF )

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

      END DO

      CALL amrex_parallel_reduce_max( MinLF )
      CALL amrex_parallel_reduce_max( MinCF )

      IF( MAX( MinLF, MinCF ) .LT. 1.0e-13_DP ) CONVERGED = .TRUE.

      IF( ITER .EQ. MAX_ITER ) &
        CALL DescribeError_Euler_MF &
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

  END SUBROUTINE InitializeMetric_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE UpdateConformalFactorAndMetric_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1) ! Metric Fields
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uMF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, iNX1, iNX2
    INTEGER  :: iX_B0(3), iX_E0(3)
    REAL(DP) :: X1, X2, Psi, h1, h2, h3

    TYPE(amrex_imultifab) :: iMF_Mask

    DO iLevel = 0, nLevels-1

      CALL MakeFineMask( iLevel, iMF_Mask, MF_uGF % BA, MF_uGF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      CALL CreateMesh_MF( iLevel, MeshX )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

          iNX1 = NodeNumberTableX(1,iNX)
          iNX2 = NodeNumberTableX(2,iNX)

          X1 = NodeCoordinate( MeshX(1), iX1, iNX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNX2 )

          Psi = uMF(iX1,iX2,iX3,nDOFX*(iMF_Psi-1)+iNX)
          h1  = Psi**2
          h2  = Psi**2 * X1
          h3  = Psi**2 * X1 * SIN( X2 )

          uGF(iX1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX) = Psi
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) = h1
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) = h2
          uGF(iX1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX) = h3

          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) = MAX( h1**2, SqrtTiny )
          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) = MAX( h2**2, SqrtTiny )
          uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) = MAX( h3**2, SqrtTiny )

          uGF(iX1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX) = h1 * h2 * h3

        END DO
        END DO
        END DO
        END DO

      END DO

      CALL DestroyMesh_MF( MeshX )

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE UpdateConformalFactorAndMetric_MF


  SUBROUTINE ComputeGeometryFromPoseidon_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uMF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_Mask

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3)

    DO iLevel = 0, nLevels-1

      CALL MakeFineMask( iLevel, iMF_Mask, MF_uGF % BA, MF_uGF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Alpha-1)+iNX)

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_1-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_2-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_Beta_3-1)+iNX)

          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_11-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_11-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_12-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_12-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_13-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_13-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_22-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_22-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_23-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_23-1)+iNX)
          uGF    (iX1,iX2,iX3,nDOFX*(iGF_K_dd_33-1)+iNX) &
            = uMF(iX1,iX2,iX3,nDOFX*(iMF_K_dd_33-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE ComputeGeometryFromPoseidon_MF


  SUBROUTINE SetPoseidonBoundaryConditions_Outer &
    ( MF_uGS, Psi_xR, AlphaPsi_xR, Beta_u_xR )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:nLevels-1)
    REAL(DP)            , INTENT(out) :: Psi_xR, AlphaPsi_xR, Beta_u_xR(3)

    REAL(DP) :: GravitationalMass
    REAL(DP) :: Alpha_xR

    CALL ComputeGravitationalMass( MF_uGS, GravitationalMass )

    ! --- Approximate outer boundary with isotropic expressions ---

    Psi_xR   = One + Half * GravitationalMass / xR(1)
    Alpha_xR =   ( One - Half * GravitationalMass / xR(1) ) &
               / ( One + Half * GravitationalMass / xR(1) )

    AlphaPsi_xR = Alpha_xR * Psi_xR
    Beta_u_xR   = Zero

  END SUBROUTINE SetPoseidonBoundaryConditions_Outer


  SUBROUTINE ComputeGravitationalMass( MF_uGS, GravitationalMass )

    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:nLevels-1)
    REAL(DP)            , INTENT(out) :: GravitationalMass

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_Mask

    INTEGER  :: iNX, iX1, iX2, iX3, iX_B0(3), iX_E0(3), iLevel
    REAL(DP) :: d3X

    ! --- Assuming 1D spherical symmetry ---

    GravitationalMass = Zero

    DO iLevel = 0, nLevels-1

      CALL MakeFineMask( iLevel, iMF_Mask, MF_uGS % BA, MF_uGS % DM )

      CALL amrex_mfiter_build( MFI, MF_uGS(iLevel), tiling = UseTiling )

      CALL CreateMesh_MF( iLevel, MeshX )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )

        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          IF( Mask(iX1,iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

          d3X = Two / Pi * MeshX(1) % Width(iX1) &
                         * MeshX(2) % Width(iX2) &
                         * MeshX(3) % Width(iX3)

          GravitationalMass &
            = GravitationalMass + d3X &
                * WeightsX_q(iNX) * uGS(iX1,iX2,iX3,nDOFX*(iGS_Mg-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO

      CALL DestroyMesh_MF( MeshX )

      CALL amrex_mfiter_destroy( MFI )

    END DO

    CALL amrex_parallel_reduce_sum( GravitationalMass )

  END SUBROUTINE ComputeGravitationalMass


  SUBROUTINE ApplyBoundaryConditions_Geometry( GravitationalMass, MF_uGF )

    REAL(DP)            , INTENT(in)    :: GravitationalMass
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    CALL ApplyBoundaryConditions_X1_Inner( MF_uGF )
    CALL ApplyBoundaryConditions_X1_Outer( GravitationalMass, MF_uGF )

  END SUBROUTINE ApplyBoundaryConditions_Geometry


  SUBROUTINE ApplyBoundaryConditions_X1_Inner( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    INTEGER , CONTIGUOUS, POINTER :: Mask(:,:,:,:)

    TYPE(amrex_imultifab) :: iMF_Mask

    INTEGER :: iLevel, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: iNX1, iNX2, iNX3, iNX
    INTEGER :: jNX1, jNX

    DO iLevel = 0, nLevels-1

      CALL MakeFineMask( iLevel, iMF_Mask, MF_uGF % BA, MF_uGF % DM )

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        Mask => iMF_Mask % DataPtr( MFI )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ! --- Inner Boundary: Reflecting ---

        IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)

            IF( Mask(iX_B0(1),iX2,iX3,1) .NE. iLeaf_MFM ) CYCLE

            DO iNX3 = 1, nNodesX(3)
            DO iNX2 = 1, nNodesX(2)
            DO iNX1 = 1, nNodesX(1)

              jNX1 = ( nNodesX(1) - iNX1 ) + 1

              iNX = NodeNumberX( iNX1, iNX2, iNX3 )
              jNX = NodeNumberX( jNX1, iNX2, iNX3 )

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Alpha-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Alpha-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Psi-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Psi-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Beta_1-1)+iNX) &
                = -uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Beta_1-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Beta_2-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Beta_2-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Beta_3-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_Beta_3-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_h_1-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_h_2-1)+jNX)

              uGF     (iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX) &
                = +uGF(iX_B0(1)  ,iX2,iX3,nDOFX*(iGF_h_3-1)+jNX)

              uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX) &
                = MAX( uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX)**2, &
                       SqrtTiny )

              uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX) &
                = MAX( uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX)**2, &
                       SqrtTiny )

              uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX) &
                = MAX( uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX)**2, &
                       SqrtTiny )

              uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_SqrtGm-1)+iNX) &
                =   uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_1-1)+iNX) &
                  * uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_2-1)+iNX) &
                  * uGF(iX_B0(1)-1,iX2,iX3,nDOFX*(iGF_h_3-1)+iNX)

            END DO
            END DO
            END DO

          END DO
          END DO

        END IF

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE ApplyBoundaryConditions_X1_Inner


  SUBROUTINE ApplyBoundaryConditions_X1_Outer( GravitationalMass, MF_uGF )

    REAL(DP)            , INTENT(in)    :: GravitationalMass
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), ALLOCATABLE :: G_K(:,:,:,:)
    REAL(DP), ALLOCATABLE :: G_F(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLevel, iX2, iX3, iGF, nX1_X
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: iNX

    REAL(DP) :: X1

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      CALL CreateMesh_MF( iLevel, MeshX )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ! --- Outer Boundary ---

        IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) )THEN

          DO iNX = 1, nNodesX(1)

            X1 = NodeCoordinate( MeshX(1), iX_E0(1)+1, iNX )

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_Alpha-1)+iNX) &
              =   ( One - Half * GravitationalMass / X1 ) &
                / ( One + Half * GravitationalMass / X1 )

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_Psi-1)+iNX) &
              = One + Half * GravitationalMass / X1

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_Beta_1-1)+iNX) = Zero

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_Beta_2-1)+iNX) = Zero

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_Beta_3-1)+iNX) = Zero

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_K_dd_11-1)+iNX) = Zero

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_K_dd_12-1)+iNX) = Zero

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_K_dd_13-1)+iNX) = Zero

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_K_dd_22-1)+iNX) = Zero

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_K_dd_23-1)+iNX) = Zero

            uGF(iX_E0(1)+1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3), &
                nDOFX*(iGF_K_dd_33-1)+iNX) = Zero

          END DO

!!$          nX1_X = ( iX_E0(3) - iX_B0(3) + 1 ) * ( iX_E0(2) - iX_B0(2) + 1 )
!!$
!!$          ALLOCATE( G_K(1:nDOFX   ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nGF) )
!!$          ALLOCATE( G_F(1:nDOFX_X1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nGF) )
!!$
!!$          DO iGF = 1       , nGF
!!$          DO iX3 = iX_B0(3), iX_E0(3)
!!$          DO iX2 = iX_B0(2), iX_E0(2)
!!$          DO iNX = 1       , nDOFX
!!$
!!$            G_K(iNX,iX2,iX3,iGF) = uGF(iX_E0(1),iX2,iX3,nDOFX*(iGF-1)+iNX)
!!$
!!$          END DO
!!$          END DO
!!$          END DO
!!$          END DO
!!$
!!$          DO iGF = 1, nGF
!!$
!!$            CALL MatrixMatrixMultiply &
!!$                   ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, &
!!$                     nDOFX_X1,   G_K(1,iX_B0(2),iX_B0(3),iGF), &
!!$                     nDOFX, Zero,G_F(1,iX_B0(2),iX_B0(3),iGF), &
!!$                     nDOFX_X1 )
!!$
!!$          END DO
!!$
!!$          DO iGF = 1       , nGF
!!$          DO iX3 = iX_B0(3), iX_E0(3)
!!$          DO iX2 = iX_B0(2), iX_E0(2)
!!$          DO iNX = 1       , nDOFX
!!$
!!$            uGF(iX_E0(1)+1,iX2,iX3,nDOFX*(iGF-1)+iNX) = G_F(1,iX2,iX3,iGF)
!!$
!!$          END DO
!!$          END DO
!!$          END DO
!!$          END DO
!!$
!!$          DEALLOCATE( G_F )
!!$          DEALLOCATE( G_K )

        END IF

      END DO

      CALL DestroyMesh_MF( MeshX )

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE ApplyBoundaryConditions_X1_Outer


END MODULE MF_GravitySolutionModule_XCFC_Poseidon
