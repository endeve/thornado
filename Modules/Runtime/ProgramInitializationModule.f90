MODULE ProgramInitializationModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, swX, xL, xR, ZoomX, &
    nE, swE, eL, eR, ZoomE, &
    nDimsX, nDimsE, nDims, &
    nDOFX,  nDOFE,  nDOF, &
    nNodesX, nNodesE, nNodes, &
    InitializeProgramHeader
  USE QuadratureModule, ONLY: &
    InitializeQuadratures
  USE PolynomialBasisModule_Lagrange, ONLY: &
    InitializePolynomialBasis_Lagrange
  USE PolynomialBasisModule_Legendre, ONLY: &
    InitializePolynomialBasis_Legendre
  USE PolynomialBasisMappingModule, ONLY: &
    InitializePolynomialBasisMapping
  USE MeshModule, ONLY: &
    MeshX, MeshE, &
    CreateMesh, &
    DestroyMesh
  USE GeometryModule, ONLY: &
    InitializeGeometry, &
    FinalizeGeometry
  USE FluidFieldsModule, ONLY: &
    CreateFluidFields, &
    DestroyFluidFields
  USE RadiationFieldsModule, ONLY: &
    CreateRadiationFields, &
    DestroyRadiationFields
  USE RadiationFieldsUtilitiesModule, ONLY: &
    InitializeRadiationFieldsUtilities, &
    FinalizeRadiationFieldsUtilities
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE OpacityModule, ONLY: &
    InitializeOpacities, &
    FinalizeOpacities
  USE FluidEvolutionModule, ONLY: &
    InitializeFluidEvolution, &
    FinalizeFluidEvolution
  USE RadiationEvolutionModule, ONLY: &
    InitializeRadiationEvolution, &
    FinalizeRadiationEvolution
  USE RiemannSolverModule, ONLY: &
    InitializeRiemannSolvers, &
    FinalizeRiemannSolvers
  USE FluidRadiationCouplingModule, ONLY: &
    InitializeFluidRadiationCoupling, &
    FinalizeFluidRadiationCoupling
  USE TimeSteppingModule, ONLY: &
    InitializeTimeStepping, &
    FinalizeTimeStepping

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeProgram
  PUBLIC :: FinalizeProgram

CONTAINS


  SUBROUTINE InitializeProgram &
    ( ProgramName_Option, nX_Option, swX_Option, bcX_Option, xL_Option, &
      xR_Option, zoomX_Option, nE_Option, swE_Option, bcE_Option, &
      eL_Option, eR_Option, zoomE_Option, CoordinateSystem_Option, &
      nNodes_Option, ActivateUnits_Option, EquationOfState_Option, &
      EquationOfStateTableName_Option, Gamma_IDEAL_Option, &
      Opacity_Option, OpacityTableName_Option, FluidSolver_Option, &
      RadiationSolver_Option, FluidRiemannSolver_Option, &
      RadiationRiemannSolver_Option, FluidRadiationCoupling_Option, &
      EvolveFluid_Option, EvolveRadiation_Option, BetaTVB_Option, &
      nStages_SSP_RK_Option, nStages_SI_RK_Option )

    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: ProgramName_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: nX_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: swX_Option
    INTEGER,  DIMENSION(3), INTENT(in), OPTIONAL :: bcX_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: xL_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: xR_Option
    REAL(DP), DIMENSION(3), INTENT(in), OPTIONAL :: zoomX_Option
    INTEGER,                INTENT(in), OPTIONAL :: nE_Option
    INTEGER,                INTENT(in), OPTIONAL :: swE_Option
    INTEGER,                INTENT(in), OPTIONAL :: bcE_Option
    REAL(DP),               INTENT(in), OPTIONAL :: eL_Option
    REAL(DP),               INTENT(in), OPTIONAL :: eR_Option
    REAL(DP),               INTENT(in), OPTIONAL :: zoomE_Option
    INTEGER,                INTENT(in), OPTIONAL :: nNodes_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: CoordinateSystem_Option
    LOGICAL,                INTENT(in), OPTIONAL :: ActivateUnits_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: EquationOfState_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: &
      EquationOfStateTableName_Option
    REAL(DP),               INTENT(in), OPTIONAL :: Gamma_IDEAL_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: Opacity_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: OpacityTableName_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: FluidSolver_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: RadiationSolver_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: FluidRiemannSolver_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: &
      RadiationRiemannSolver_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: &
      FluidRadiationCoupling_Option
    LOGICAL,                INTENT(in), OPTIONAL :: EvolveFluid_Option
    LOGICAL,                INTENT(in), OPTIONAL :: EvolveRadiation_Option
    REAL(DP),               INTENT(in), OPTIONAL :: BetaTVB_Option
    INTEGER,                INTENT(in), OPTIONAL :: nStages_SSP_RK_Option
    INTEGER,                INTENT(in), OPTIONAL :: nStages_SI_RK_Option

    LOGICAL :: ActivateUnits
    INTEGER :: iDim

    CALL InitializeProgramHeader &
           ( ProgramName_Option = ProgramName_Option, &
             nX_Option = nX_Option, swX_Option = swX_Option, &
             bcX_Option = bcX_Option, xL_Option = xL_Option, &
             xR_Option = xR_Option, zoomX_Option = zoomX_Option, &
             nE_Option = nE_Option, swE_Option = swE_Option, &
             bcE_Option = bcE_Option, eL_Option = eL_Option, &
             eR_Option = eR_Option, zoomE_Option = zoomE_Option, &
             nNodes_Option = nNodes_Option )

    WRITE(*,*)
    WRITE(*,'(A2,A28,A)') &
      '', 'INFO: Initializing Program: ', TRIM( ProgramName )
    WRITE(*,*)

    ! --- Units ---

    ActivateUnits = .FALSE.
    IF( PRESENT( ActivateUnits_Option ) ) &
      ActivateUnits = ActivateUnits_Option

    IF( ActivateUnits )THEN

      CALL ActivateUnitsDisplay
      CALL DescribeUnitsDisplay

    END IF

    ! --- Problem Dimensionality ---

    nDimsX = 0
    DO iDim = 1, 3
      nNodesX(iDim) = 1
      IF( nX(iDim) > 1 )THEN
        nDimsX = nDimsX + 1
        nNodesX(iDim) = nNodes
      END IF
    END DO

    nDimsE  = 0
    nNodesE = 1
    IF( nE > 1 )THEN
      nDimsE  = 1
      nNodesE = nNodes
    END IF

    nDims = nDimsX + nDimsE

    WRITE(*,'(A5,A17,I1)') &
      '', 'Dimensionality = ', nDims
    WRITE(*,'(A5,A18)') &
      '', '------------------'
    WRITE(*,*)
    WRITE(*,'(A7,A9,I1,A2,A9,I1)') &
      '', 'nDimsX = ', nDimsX, '', 'nDimsE = ', nDimsE
    WRITE(*,*)

    WRITE(*,'(A7,A9,I4.4,A2,A9,I4.4,A2,A9,I4.4,A2,A6,I4.4)') &
      '', 'nX(1) = ', nX(1), '', 'nX(2) = ', nX(2), &
      '', 'nX(3) = ', nX(3), '', 'nE = ', nE
    WRITE(*,'(A7,A9,I4,A2,A9,I4,A2,A9,I4,A2,A6,I4)') &
      '', 'swX(1) = ', swX(1), '', 'swX(2) = ', swX(2), &
      '', 'swX(3) = ', swX(3), '', 'swE = ', swE
    WRITE(*,*)

    ! --- Degrees of Freedom Per Element Per Physical Field ---

    nDOFX = nNodes**nDimsX
    nDOFE = nNodes**nDimsE
    nDOF  = nNodes**nDims

    WRITE(*,'(A5,A36)') '', 'Degrees of Freedom / Element / Field'
    WRITE(*,*)
    WRITE(*,'(A7,A9,I2.2)') &
      '', 'nNodes = ', nNodes
    WRITE(*,*)
    DO iDim = 1, 3
      WRITE(*,'(A9,A4,I1,A10,I1,A4,I2.2)') &
        '', 'i = ', iDim, ', nNodesX(', iDim, ') = ', nNodesX(iDim)
    END DO
    WRITE(*,'(A16,A13,I2.2)') &
      '', 'nNodesE    = ', nNodesE
    WRITE(*,*)
    WRITE(*,'(A7,A8,I4.4,A2,A8,I4.4,A2,A7,I4.4)') &
      '', 'nDOFX = ', nDOFX, '', 'nDOFE = ', nDOFE, '', 'nDOF = ', nDOF
    WRITE(*,*)

    ! --- Quadratures ---

    CALL InitializeQuadratures

    ! --- Polynomial Basis ---

    CALL InitializePolynomialBasis_Lagrange

    CALL InitializePolynomialBasis_Legendre

    ! --- Spatial Grid ---

    WRITE(*,'(A5,A20)') '', 'Computational Domain'
    WRITE(*,'(A5,A20)') '', '--------------------'
    WRITE(*,*)
    DO iDim = 1, 3

      WRITE(*,'(A7,A3,I1,A4,ES8.2E2,A2,A3,I1,A4,&
                &ES8.2E2,A2,A6,I1,A4,ES10.4E2)') &
        '',   'xL(', iDim, ') = ', xL(iDim), &
        ', ', 'xR(', iDim, ') = ', xR(iDim), &
        ', ', 'ZoomX(', iDim, ') = ', ZoomX(iDim)
      WRITE(*,*)

      CALL CreateMesh( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
                       xL(iDim), xR(iDim), ZoomOption = ZoomX(iDim) )

      WRITE(*,'(A9,A11,I1,A4,ES8.2E2,A3,ES8.2E2)') &
        '', 'MIN/MAX dx(', iDim, ') = ', &
        MINVAL( MeshX(iDim) % Width ), ' / ', MAXVAL( MeshX(iDim) % Width )
      WRITE(*,*)

    END DO

    ! --- Spectral Grid ---

    WRITE(*,'(A7,A8,ES8.2E2,A2,A8,ES8.2E2,A2,A11,ES10.4E2)') &
      '', 'eL    = ', eL, ', ', 'eR    = ', eR, &
      ', ', 'ZoomE    = ', ZoomE
    WRITE(*,*)

    CALL CreateMesh( MeshE, nE, nNodesE, swE, eL, eR, ZoomOption = ZoomE )

    WRITE(*,'(A9,A16,ES8.2E2,A3,ES8.2E2)') &
      '', 'MIN/MAX de    = ', &
      MINVAL( MeshE % Width ), ' / ', MAXVAL( MeshE % Width )

    ! --- Geometry ---

    CALL InitializeGeometry &
           ( nX, nNodesX, swX, MeshX, nE, nNodesE, swE, MeshE, &
             CoordinateSystem_Option = CoordinateSystem_Option )

    ! --- Physical Fields ---

    CALL CreateFluidFields( nX, swX )

    CALL CreateRadiationFields( nX, swX, nE, swE )

    CALL InitializeRadiationFieldsUtilities &
           ( MeshE % Weights, MeshX(1) % Weights, MeshX(2) % Weights, &
             MeshX(3) % Weights )

    ! --- For Mapping Between Nodal and Modal Representations ---

    CALL InitializePolynomialBasisMapping &
           ( MeshE % Nodes, MeshX(1) % Nodes, MeshX(2) % Nodes, MeshX(3) % Nodes )

    ! --- Equation of State ---

    CALL InitializeEquationOfState &
           ( EquationOfState_Option &
               = EquationOfState_Option, &
             EquationOfStateTableName_Option &
               = EquationOfStateTableName_Option, &
             Gamma_IDEAL_Option = Gamma_IDEAL_Option )

    ! --- Opacities ---

    CALL InitializeOpacities &
           ( Opacity_Option = Opacity_Option, &
             OpacityTableName_Option = OpacityTableName_Option )

    ! --- Fluid Solver ---

    CALL InitializeFluidEvolution &
           ( FluidSolver_Option = FluidSolver_Option, &
             BetaTVB_Option = BetaTVB_Option )

    ! --- Radiation Solver ---

    CALL InitializeRadiationEvolution &
           ( RadiationSolver_Option = RadiationSolver_Option, &
             BetaTVB_Option = BetaTVB_Option )

    ! --- Riemann Solvers ---

    CALL InitializeRiemannSolvers &
           ( FluidRiemannSolver_Option &
               = FluidRiemannSolver_Option, &
             RadiationRiemannSolver_Option &
               = RadiationRiemannSolver_Option )

    ! --- Fluid-Radiation Solver ---

    CALL InitializeFluidRadiationCoupling &
           ( FluidRadiationCoupling_Option &
               = FluidRadiationCoupling_Option )

    ! --- Time Stepping ---

    CALL InitializeTimeStepping &
           ( EvolveFluid_Option = EvolveFluid_Option, &
             EvolveRadiation_Option = EvolveRadiation_Option, &
             nStages_SSP_RK_Option = nStages_SSP_RK_Option, &
             nStages_SI_RK_Option = nStages_SI_RK_Option )

  END SUBROUTINE InitializeProgram


  SUBROUTINE FinalizeProgram

    INTEGER :: iDim

    ! --- Spatial Grid ---

    DO iDim = 1, 3

      CALL DestroyMesh( MeshX(iDim) )

    END DO

    ! --- Spectral Grid ---

    CALL DestroyMesh( MeshE )

    ! --- Geometry ---

    CALL FinalizeGeometry

    ! --- Physical Fields ---

    CALL DestroyFluidFields

    CALL DestroyRadiationFields

    CALL FinalizeRadiationFieldsUtilities

    ! --- Equation of State ---

    CALL FinalizeEquationOfState

    ! --- Opacities ---

    CALL FinalizeOpacities

    ! --- Fluid Solver ---

    CALL FinalizeFluidEvolution

    ! --- Radiation Solver ---

    CALL FinalizeRadiationEvolution

    ! --- Riemann Solvers ---

    CALL FinalizeRiemannSolvers

    ! --- Fluid-Radiation Solver ---

    CALL FinalizeFluidRadiationCoupling

    ! --- Time Stepping ---

    CALL FinalizeTimeStepping

  END SUBROUTINE FinalizeProgram


END MODULE ProgramInitializationModule
