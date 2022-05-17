MODULE MF_GravitySolutionModule_XCFC_Poseidon

  ! --- AMReX Modules ---

  USE amrex_box_module, ONLY: &
    amrex_box
  USE amrex_multifab_module, ONLY: &
    amrex_multifab, &
    amrex_mfiter, &
    amrex_mfiter_build, &
    amrex_mfiter_destroy
  USE amrex_amrcore_module, ONLY: &
    amrex_geom

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

  ! --- Local Modules ---

  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    SqrtTiny
  USE MF_MeshModule, ONLY: &
    CreateMesh_MF, &
    DestroyMesh_MF
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
  USE Source_Input_AMReX, ONLY: &
    Poseidon_Input_Sources1_AMReX, &
    Poseidon_Input_Sources2_AMReX
  USE Poseidon_XCFC_Interface_Module, ONLY: &
    Poseidon_XCFC_Run_Part1, &
    Poseidon_XCFC_Run_Part2
  USE Return_Functions_AMReX, ONLY: &
    Poseidon_Return_ConFactor_AMReX, &
    Poseidon_Return_ALL_AMReX
  USE Initial_Guess_Module, ONLY: &
    Poseidon_Init_FlatGuess

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeGravitySolver_XCFC_Poseidon_MF
  PUBLIC :: FinalizeGravitySolver_XCFC_Poseidon_MF
  PUBLIC :: ComputeConformalFactor_Poseidon_MF
  PUBLIC :: ComputeGeometry_Poseidon_MF

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


  SUBROUTINE InitializeGravitySolver_XCFC_Poseidon_MF &
    ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    INTEGER,  INTENT(in)    :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(inout) :: uGF(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

!!$    CALL ComputeGeometryX( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    WRITE(*,*)
    WRITE(*,'(4x,A)') &
      'INFO: Gravity Solver (Poseidon, XCFC)'
    WRITE(*,'(4x,A)') &
      '-------------------------------------'
    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Only implemented for 1D spherical symmetry.'
    WRITE(*,*)

    CALL Initialize_Poseidon_with_AMReX &
        (   FEM_Degree_Option           = MAX( 1, nNodes - 1 ), &
            L_Limit_Option              = 0,                    &
            Units_Option                = 'G',                  &
            Domain_Edge_Option          = [ xL(1), xR(1) ],     &
            Coarse_NE_Option            = nX,                   &
            NQ_Option                   = [ nNodes, 1, 1 ],     &
            Convergence_Criteria_Option = 1.0e-08_DP,           &
            AMReX_Max_Levels_Option     = nLevels-1,            &
            AMReX_Max_Grid_Size_Option  = MaxGridSizeX,         &
            Verbose_Option              = .FALSE.,              &
            Print_Setup_Option          = .TRUE. )

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

    REAL(DP)         :: Psi_BC, AlphaPsi_BC
    CHARACTER(LEN=1) :: INNER_BC_TYPES (5), OUTER_BC_TYPES (5)
    REAL(DP)         :: INNER_BC_VALUES(5), OUTER_BC_VALUES(5)

    REAL(DP) :: GravitationalMass

    TYPE(amrex_multifab) :: MF_uMF(0:nLevels-1) ! Metric Fields

!!$    CALL TimersStart_Euler( Timer_GravitySolver )

#ifdef GRAVITY_SOLVER_POSEIDON_CFA

    ! --- Set Boundary Values ---

    ! May not need this with improved boundary conditions
    CALL ComputeGravitationalMass( MF_uGF, MF_uGS, GravitationalMass )

!    Psi_BC      = ConformalFactor( xR(1), GravitationalMass )
!    AlphaPsi_BC = LapseFunction  ( xR(1), GravitationalMass ) * Psi_BC

    INNER_BC_TYPES = [ 'N', 'N', 'N', 'N', 'N' ] ! Neumann
    OUTER_BC_TYPES = [ 'D', 'D', 'D', 'D', 'D' ] ! Dirichlet

    INNER_BC_VALUES = [ Zero  , Zero       , Zero, Zero, Zero ]
    OUTER_BC_VALUES = [ Psi_BC, AlphaPsi_BC, Zero, Zero, Zero ]

    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions &
           ( "I", INNER_BC_TYPES, INNER_BC_VALUES )
    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions &
           ( "O", OUTER_BC_TYPES, OUTER_BC_VALUES)

    CALL Poseidon_Init_FlatGuess() ! Possibly move this to init call

! --- Waiting for update to Poseidon that removes reference element info
!     from argument list ---

!!$    ! --- Set matter sources with current conformal factor ---
!!$    CALL Poseidon_Input_Sources1_AMReX &
!!$           ( MF_Src_Input  = MF_uGS,           &
!!$             MF_Src_nComps = nGS,              &
!!$             num_levels    = nLevels,          &
!!$             Input_NQ      = nNodesX,          &
!!$             Input_R_Quad  = MeshX(1) % Nodes, &
!!$             Input_T_Quad  = MeshX(2) % Nodes, &
!!$             Input_P_Quad  = MeshX(3) % Nodes, &
!!$             Input_xL      = [ -Half, +Half ] )

    ! --- Compute conformal factor ---

    CALL Poseidon_XCFC_Run_Part1()

! --- Waiting for update to Poseidon that removes reference element info
!     from argument list ---

!!$    CALL Poseidon_Return_ConFactor_AMReX &
!!$         ( NQ          = nNodesX,          &
!!$           RQ_Input    = MeshX(1) % Nodes, &
!!$           TQ_Input    = MeshX(2) % Nodes, &
!!$           PQ_Input    = MeshX(3) % Nodes, &
!!$           Left_Limit  = -Half,            &
!!$           Right_Limit = +Half,            &
!!$           nLevels     = nLevels,          &
!!$           MF_Results  = MF_uMF )


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

! --- Waiting for update to Poseidon that removes reference element info
!     from argument list ---

!!$    CALL Poseidon_Input_Sources2_AMReX &
!!$           ( MF_Src_Input  = MF_uGS,           &
!!$             MF_Src_nComps = nGS,              &
!!$             num_levels    = nLevels,          &
!!$             Input_NQ      = nNodesX,          &
!!$             Input_R_Quad  = MeshX(1) % Nodes, &
!!$             Input_T_Quad  = MeshX(2) % Nodes, &
!!$             Input_P_Quad  = MeshX(3) % Nodes, &
!!$             Input_xL      = [ -Half, +Half ] )

    ! --- Compute lapse and shift ---

    CALL Poseidon_XCFC_Run_Part2()

! --- Waiting for update to Poseidon that removes reference element info
!     from argument list ---

!!$    CALL Poseidon_Return_ALL_AMReX &
!!$         ( NQ               = nNodesX,          &
!!$           RQ_Input         = MeshX(1) % Nodes, &
!!$           TQ_Input         = MeshX(2) % Nodes, &
!!$           PQ_Input         = MeshX(3) % Nodes, &
!!$           Left_Limit       = -Half,            &
!!$           Right_Limit      = +Half,            &
!!$           nLevels          = nLevels,          &
!!$           MF_Results       = MF_uMF )

    ! --- Copy data from Poseidon to thornado ---

    CALL ComputeGeometryFromPoseidon_MF( MF_uMF, MF_uGF )

    CALL SetBoundaryConditions( MF_uGF )

#endif

!!$    CALL TimersStop_Euler( Timer_GravitySolver )

  END SUBROUTINE ComputeGeometry_Poseidon_MF


  ! --- PRIVATE SUBROUTINES ---


  SUBROUTINE UpdateConformalFactorAndMetric_MF( MF_uMF, MF_uGF )

    TYPE(amrex_multifab), INTENT(in)    :: MF_uMF(0:nLevels-1) ! Metric Fields
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uMF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)

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

    REAL(DP), CONTIGUOUS, POINTER :: uMF (:,:,:,:)
    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)

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


  SUBROUTINE SetBoundaryConditions( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    CALL SetBoundaryConditions_X1( MF_uGF )

  END SUBROUTINE SetBoundaryConditions


  SUBROUTINE SetBoundaryConditions_X1( MF_uGF )

    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)

    TYPE(amrex_box)    :: BX
    TYPE(amrex_mfiter) :: MFI

    REAL(DP), CONTIGUOUS, POINTER :: uGF (:,:,:,:)
    REAL(DP), ALLOCATABLE :: G_K(:,:,:,:)
    REAL(DP), ALLOCATABLE :: G_F(:,:,:,:)

    INTEGER :: iLevel, iX1, iX2, iX3, iGF, nX1_X
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

        ! --- Upper Boundary ---

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

          DEALLOCATE( G_F )
          DEALLOCATE( G_K )

        END IF

      END DO

      CALL amrex_mfiter_destroy( MFI )

    END DO

  END SUBROUTINE SetBoundaryConditions_X1


  SUBROUTINE ComputeGravitationalMass( MF_uGF, MF_uGS, GravitationalMass )
    TYPE(amrex_multifab), INTENT(in)  :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(in)  :: MF_uGS(0:nLevels-1)
    REAL(DP)            , INTENT(out) :: GravitationalMass
  END SUBROUTINE ComputeGravitationalMass


END MODULE MF_GravitySolutionModule_XCFC_Poseidon
