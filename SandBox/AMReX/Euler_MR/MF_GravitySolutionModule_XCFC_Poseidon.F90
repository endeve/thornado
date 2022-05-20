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
    amrex_mfiter_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom
  USE amrex_parallel_module, ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_reduce_min

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
    nDOFX_X1
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
    Three, &
    SqrtTiny
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
  USE MF_Euler_UtilitiesModule, ONLY: &
    ComputeConserved_Euler_MF
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
    Poseidon_Close, &
    Poseidon_CFA_Set_Uniform_Boundary_Conditions
  USE Poseidon_Source_Input_Module, ONLY: &
    Poseidon_Input_Sources_Part1, &
    Poseidon_Input_Sources_Part2
  USE Poseidon_XCFC_Interface_Module, ONLY: &
    Poseidon_XCFC_Run_Part1, &
    Poseidon_XCFC_Run_Part2
  USE Poseidon_Return_Routines_Module, ONLY: &
    Poseidon_Return_Conformal_Factor, &
    Poseidon_Return_ALL
  USE Initial_Guess_Module, ONLY: &
    Poseidon_Init_FlatGuess

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_XCFC_Poseidon_MF
  PUBLIC :: FinalizeGravitySolver_XCFC_Poseidon_MF
  PUBLIC :: ComputeConformalFactor_Poseidon_MF
  PUBLIC :: ComputeGeometry_Poseidon_MF
  PUBLIC :: ComputeConformalFactorSources_XCFC_MF
  PUBLIC :: ComputePressureTensorTrace_XCFC_MF
  PUBLIC :: MultiplyWithPsi6_MF
  PUBLIC :: InitializeMetric_MF

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

  INTEGER, PARAMETER :: iGS_E  = 1
  INTEGER, PARAMETER :: iGS_S1 = 2
  INTEGER, PARAMETER :: iGS_S2 = 3
  INTEGER, PARAMETER :: iGS_S3 = 4
  INTEGER, PARAMETER :: iGS_S  = 5
  INTEGER, PARAMETER :: iGS_Mg = 6
  INTEGER, PARAMETER :: nGS    = 6

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
           ( Source_NQ          = nNodesX, &
             Source_xL          = [ -Half, +Half ], &
             Source_RQ_xlocs    = MeshX(1) % Nodes, &
             Source_TQ_xlocs    = MeshX(2) % Nodes, &
             Source_PQ_xlocs    = MeshX(3) % Nodes, &
             Units_Option       = 'G', &
             Verbose_Option     = .TRUE., &
             Print_Setup_Option = .TRUE. )

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

    REAL(DP)         :: Psi_xR, AlphaPsi_xR
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    REAL(DP) :: GravitationalMass

    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1) ! Metric Fields

!!$    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! --- Set Boundary Values ---

    CALL SetBoundaryConditions_Outer &
           ( MF_uGF, Psi_xR_Option = Psi_xR, AlphaPsi_xR_Option = AlphaPsi_xR )

    INNER_BC_TYPES = [ 'N', 'N', 'N', 'N', 'N' ] ! Neumann
    OUTER_BC_TYPES = [ 'D', 'D', 'D', 'D', 'D' ] ! Dirichlet

    INNER_BC_VALUES = [ Zero  , Zero       , Zero, Zero, Zero ]
    OUTER_BC_VALUES = [ Psi_xR, AlphaPsi_xR, Zero, Zero, Zero ]

    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions &
           ( "I", INNER_BC_TYPES, INNER_BC_VALUES )
    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions &
           ( "O", OUTER_BC_TYPES, OUTER_BC_VALUES)

    ! --- Set matter sources with current conformal factor ---
    CALL Poseidon_Input_Sources_Part1 &
           ( MF_Src_Input  = MF_uGS, &
             MF_Src_nComps = nGS )

    CALL Poseidon_Init_FlatGuess() ! Possibly move this to init call

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

    CALL Poseidon_Return_Conformal_Factor( MF_Results = MF_uMF )

    CALL UpdateConformalFactorAndMetric_MF( MF_uMF, MF_uGF )

#endif

!!$    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeConformalFactor_Poseidon_MF


  SUBROUTINE ComputeGeometry_Poseidon_MF( MF_uGS, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGS(0:nLevels-1) ! Gravity Sources
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1) ! Metric Fields

!!$    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! --- Set gravity sources with updated conformal factor ---

    CALL Poseidon_Input_Sources_Part2 &
           ( MF_Src_Input  = MF_uGS, &
             MF_Src_nComps = nGS )

    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()

    CALL Poseidon_Return_ALL( MF_Results = MF_uMF )

    ! --- Copy data from Poseidon to thornado ---

    CALL ComputeGeometryFromPoseidon_MF( MF_uMF, MF_uGF )

    CALL SetBoundaryConditions_Inner( MF_uGF )
    CALL SetBoundaryConditions_Outer( MF_uGF )

#endif

!!$    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeGeometry_Poseidon_MF


  SUBROUTINE ComputeConformalFactorSources_XCFC_MF &
    ( MF_uGF, MF_uCF, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! This is U^{*}
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )
        uGS => MF_uGS(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

          uGS      (iX1,iX2,iX3,nDOFX*(iGS_E-1)+iNX) &
            =   uCF(iX1,iX2,iX3,nDOFX*(iCF_E-1)+iNX) &
              + uCF(iX1,iX2,iX3,nDOFX*(iCF_D-1)+iNX)

          uGS        (iX1,iX2,iX3,nDOFX*(iGS_S1      -1)+iNX) &
            = uCF    (iX1,iX2,iX3,nDOFX*(iCF_S1      -1)+iNX) &
                / uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_11-1)+iNX)

          uGS        (iX1,iX2,iX3,nDOFX*(iGS_S2      -1)+iNX) &
            = uCF    (iX1,iX2,iX3,nDOFX*(iCF_S2      -1)+iNX) &
                / uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_22-1)+iNX)

          uGS        (iX1,iX2,iX3,nDOFX*(iGS_S3      -1)+iNX) &
            = uCF    (iX1,iX2,iX3,nDOFX*(iCF_S3      -1)+iNX) &
                / uGF(iX1,iX2,iX3,nDOFX*(iGF_Gm_dd_33-1)+iNX)

        END DO
        END DO
        END DO
        END DO

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE ComputeConformalFactorSources_XCFC_MF


  SUBROUTINE ComputePressureTensorTrace_XCFC_MF( MF_uGF, MF_uCF, MF_uGS )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)    :: MF_uCF(0:nLevels-1) ! This is U^{*}
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGS(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uCF(:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGS(:,:,:,:)

    INTEGER :: iLevel, iNX, iX1, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: jErr
    INTEGER, ALLOCATABLE :: iErr(:,:,:,:)

    REAL(DP) :: uPF(nPF), Pressure, Psi6

!    CALL TimersStart_Euler( Timer_GS_ComputeSourceTerms )

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

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

          uGF(iX1,iX2,iX3,nDOFX*(iGS_S-1)+iNX) &
            =   uCF(iX1,iX2,iX3,nDOFX*(iCF_S1-1)+iNX) * uPF(iPF_V1) &
              + uCF(iX1,iX2,iX3,nDOFX*(iCF_S2-1)+iNX) * uPF(iPF_V2) &
              + uCF(iX1,iX2,iX3,nDOFX*(iCF_S3-1)+iNX) * uPF(iPF_V3) &
              + Psi6 * Three * Pressure

        END DO
        END DO
        END DO
        END DO

        IF( jErr .NE. 0 )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)
          DO iNX = 1, nDOFX

            IF( iErr(iNX,iX1,iX2,iX3) .NE. 0 )THEN

              WRITE(*,*) 'ERROR'
              WRITE(*,*) '-----'
              WRITE(*,*) '    MODULE: Poseidon_UtilitiesModule'
              WRITE(*,*) 'SUBROUTINE: ComputePressureTensorTrace_Poseidon'
              WRITE(*,*) 'iX_B0: ', iX_B0
              WRITE(*,*) 'iX_E0: ', iX_E0
              WRITE(*,*) 'iX1, iX2, iX3: ', iX1, iX2, iX3

              CALL DescribeError_Euler &
                ( iErr(iNX,iX1,iX2,iX3), &
                  Int_Option = [ iNX ], &
                  Real_Option &
                    = [ uCF(iX1,iX2,iX3,nDOFX*(iCF_D       -1)+iNX), &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_S1      -1)+iNX), &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_S2      -1)+iNX), &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_S3      -1)+iNX), &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_E       -1)+iNX), &
                        uCF(iX1,iX2,iX3,nDOFX*(iCF_Ne      -1)+iNX), &
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

    REAL(DP) :: Psi6

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )
        uCF => MF_uCF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

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
    TYPE(amrex_multifab) :: Al1   (0:nLevels-1)
    TYPE(amrex_multifab) :: Al2   (0:nLevels-1)
    TYPE(amrex_multifab) :: dAl   (0:nLevels-1)
    TYPE(amrex_multifab) :: CF1   (0:nLevels-1)
    TYPE(amrex_multifab) :: CF2   (0:nLevels-1)
    TYPE(amrex_multifab) :: dCF   (0:nLevels-1)

    LOGICAL  :: CONVERGED
    INTEGER  :: iLevel, ITER
    REAL(DP) :: MinAl, MinCF

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
             ( MF_uGS(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX * nGS, 0 )

      CALL amrex_multifab_build &
             ( Al1(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )

      CALL amrex_multifab_build &
             ( Al2(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )

      CALL amrex_multifab_build &
             ( dAl(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )

      CALL amrex_multifab_build &
             ( CF1(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )

      CALL amrex_multifab_build &
             ( CF2(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )

      CALL amrex_multifab_build &
             ( dCF(iLevel), MF_uGF(iLevel) % BA, MF_uGF(iLevel) % DM, &
               nDOFX, 0 )

    END DO

    ! --- Iterate to incorporate gravity in initial conditions ---

    CONVERGED = .FALSE.
    ITER = 0

    DO WHILE( .NOT. CONVERGED )

      ITER = ITER + 1

      DO iLevel = 0, nLevels - 1

        CALL Al1(iLevel) % COPY &
              ( MF_uGF(iLevel), 1+nDOFX*(iGF_Alpha-1), 1, nDOFX, 0 )

        CALL CF1(iLevel) % COPY &
              ( MF_uGF(iLevel), 1+nDOFX*(iGF_Psi  -1), 1, nDOFX, 0 )

      END DO

      CALL MultiplyWithPsi6_MF( MF_uGF, +1, MF_uCF )

      CALL ComputeConformalFactorSources_XCFC_MF( MF_uGF, MF_uCF, MF_uGS )

      CALL ComputeConformalFactor_Poseidon_MF( MF_uGS, MF_uGF )

      CALL ComputePressureTensorTrace_XCFC_MF( MF_uGF, MF_uCF, MF_uGS )

      CALL ComputeGeometry_Poseidon_MF( MF_uGS, MF_uGF )

      DO iLevel = 0, nLevels - 1

        CALL Al2(iLevel) % COPY &
              ( MF_uGF(iLevel), 1+nDOFX*(iGF_Alpha-1), 1, nDOFX, 0 )

        CALL CF2(iLevel) % COPY &
              ( MF_uGF(iLevel), 1+nDOFX*(iGF_Psi  -1), 1, nDOFX, 0 )

        CALL dAl(iLevel) &
               % LinComb( +One, Al2(iLevel), 1, &
                          -One, Al1(iLevel), 1, 1, &
                          nDOFX, 0 )

        CALL dCF(iLevel) &
               % LinComb( +One, CF2(iLevel), 1, &
                          -One, CF1(iLevel), 1, 1, &
                          nDOFX, 0 )

       MinAl = dAl(iLevel) % Norm0( nDOFX )
       MinCF = dCF(iLevel) % Norm0( nDOFX )

      END DO

      CALL amrex_parallel_reduce_min( MinAl )
      CALL amrex_parallel_reduce_min( MinCF )

      CALL ComputeConserved_Euler_MF( MF_uGF, MF_uPF, MF_uAF, MF_uCF )

      IF( MAX( MinAl, MinCF ) .LT. 1.0e-13_DP ) CONVERGED = .TRUE.

      IF( ITER .EQ. 10 )THEN

        WRITE(*,*) 'Could not initialize fields. Exiting...'
        STOP

      END IF

    END DO ! WHILE( .NOT. CONVERGED )

print*,'yay! got to MF_GravitySolutionModule_XCFC_Poseidon.F90, line 687'
stop 'MF_GravitySolutionModule_XCFC_Poseidon (line 688)'

!!$    CALL ApplySlopeLimiter_Euler_Relativistic_TABLE &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )
!!$
!!$    CALL ApplyPositivityLimiter_Euler_Relativistic_TABLE &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF )
!!$
!!$    CALL MultiplyByPsi6( iX_B1, iX_E1, uGF, uCF )
!!$
!!$    CALL ComputeMatterSources_Poseidon &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, E, Si, Mg )
!!$
!!$    CALL ComputeConformalFactor_Poseidon &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, E, Si, Mg, uGF )
!!$
!!$    CALL ComputePressureTensorTrace_Poseidon &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, S )
!!$
!!$    CALL ComputeGeometry_Poseidon &
!!$           ( iX_B0, iX_E0, iX_B1, iX_E1, E, S, Si, uGF )
!!$
!!$    CALL DivideByPsi6( iX_B1, iX_E1, uGF, uCF )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( dCF   (iLevel) )
      CALL amrex_multifab_destroy( CF2   (iLevel) )
      CALL amrex_multifab_destroy( CF1   (iLevel) )
      CALL amrex_multifab_destroy( dAl   (iLevel) )
      CALL amrex_multifab_destroy( Al2   (iLevel) )
      CALL amrex_multifab_destroy( Al1   (iLevel) )
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

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3, iNX1, iNX2
    INTEGER  :: iX_B0(3), iX_E0(3)
    REAL(DP) :: X1, X2, Psi, h1, h2, h3

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      CALL CreateMesh_MF( iLevel, MeshX )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

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

    INTEGER  :: iLevel, iNX, iX1, iX2, iX3
    INTEGER  :: iX_B0(3), iX_E0(3)

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uMF => MF_uMF(iLevel) % DataPtr( MFI )
        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
        DO iNX = 1       , nDOFX

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


  SUBROUTINE SetBoundaryConditions_Inner( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    CALL SetBoundaryConditions_X1_Inner( MF_uGF )

  END SUBROUTINE SetBoundaryConditions_Inner


  SUBROUTINE SetBoundaryConditions_Outer &
    ( MF_uGF, Psi_xR_Option, AlphaPsi_xR_Option )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    REAL(DP)            , INTENT(out), OPTIONAL :: &
      Psi_xR_Option, AlphaPsi_xR_Option

    REAL(DP) :: Psi_xR
    REAL(DP) :: AlphaPsi_xR

    CALL SetBoundaryConditions_X1_Outer( MF_uGF, Psi_xR, AlphaPsi_xR )

    IF( PRESENT( Psi_xR_Option      ) ) Psi_xR_Option      = Psi_xR
    IF( PRESENT( AlphaPsi_xR_Option ) ) AlphaPsi_xR_Option = AlphaPsi_xR

  END SUBROUTINE SetBoundaryConditions_Outer


  SUBROUTINE SetBoundaryConditions_X1_Inner( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)

    INTEGER :: iLevel, iX2, iX3
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: iNX1, iNX2, iNX3, iNX
    INTEGER :: jNX1, jNX

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ! --- Inner Boundary: Reflecting ---

        IF( iX_B0(1) .EQ. amrex_geom(iLevel) % domain % lo( 1 ) )THEN

          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)

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

  END SUBROUTINE SetBoundaryConditions_X1_Inner


  SUBROUTINE SetBoundaryConditions_X1_Outer( MF_uGF, Psi_xR, AlphaPsi_xR )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    REAL(DP)            , INTENT(out)   :: Psi_xR, AlphaPsi_xR

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF(:,:,:,:)
    REAL(DP), ALLOCATABLE :: G_K(:,:,:,:)
    REAL(DP), ALLOCATABLE :: G_F(:,:,:,:)

    INTEGER :: iLevel, iX2, iX3, iGF, nX1_X
    INTEGER :: iX_B0(3), iX_E0(3)
    INTEGER :: iNX1, iNX2, iNX3, iNX
    INTEGER :: jNX1, jNX

    Psi_xR      = -HUGE( One )
    AlphaPsi_xR = -HUGE( One )

    DO iLevel = 0, nLevels-1

      CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

      DO WHILE( MFI % next() )

        uGF => MF_uGF(iLevel) % DataPtr( MFI )

        BX = MFI % tilebox()

        iX_B0 = BX % lo
        iX_E0 = BX % hi

        ! --- Outer Boundary ---

        IF( iX_E0(1) .EQ. amrex_geom(iLevel) % domain % hi( 1 ) )THEN

          nX1_X = ( iX_E0(3) - iX_B0(3) + 1 ) * ( iX_E0(2) - iX_B0(2) + 1 )

          ALLOCATE( G_K(1:nDOFX   ,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nGF) )
          ALLOCATE( G_F(1:nDOFX_X1,iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),1:nGF) )

          DO iGF = 1       , nGF
          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iNX = 1       , nDOFX

            G_K(iNX,iX2,iX3,iGF) = uGF(iX_E0(1),iX2,iX3,nDOFX*(iGF-1)+iNX)

          END DO
          END DO
          END DO
          END DO

          DO iGF = 1, nGF

            CALL MatrixMatrixMultiply &
                   ( 'N', 'N', nDOFX_X1, nX1_X, nDOFX, One, LX_X1_Up, &
                     nDOFX_X1,   G_K(1,iX_B0(2),iX_B0(3),iGF), &
                     nDOFX, Zero,G_F(1,iX_B0(2),iX_B0(3),iGF), &
                     nDOFX_X1 )

          END DO

          DO iGF = 1       , nGF
          DO iX3 = iX_B0(3), iX_E0(3)
          DO iX2 = iX_B0(2), iX_E0(2)
          DO iNX = 1       , nDOFX

            uGF(iX_E0(1)+1,iX2,iX3,nDOFX*(iGF-1)+iNX) = G_F(1,iX2,iX3,iGF)

          END DO
          END DO
          END DO
          END DO

          Psi_xR      = G_F(1,1,1,iGF_Psi)
          AlphaPsi_xR = G_F(1,1,1,iGF_Alpha) * G_F(1,1,1,iGF_Psi)

          DEALLOCATE( G_F )
          DEALLOCATE( G_K )

        END IF

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE SetBoundaryConditions_X1_Outer


  SUBROUTINE ComputeGravitationalMass( MF_uGF, MF_uGS, GravitationalMass )
    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:nLevels-1)
    REAL(DP)            , INTENT(out) :: GravitationalMass
  END SUBROUTINE ComputeGravitationalMass


END MODULE MF_GravitySolutionModule_XCFC_Poseidon
