MODULE ProgramInitializationModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    ActivateUnitsDisplay, &
    DescribeUnitsDisplay, &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, swX, xL, xR, ZoomX, &
    nE, swE, eL, eR, ZoomE, &
    nDimsX, nDimsE, nDims, &
    nDOFX,  nDOFE,  nDOF, &
    nNodesX, nNodesE, nNodes, &
    InitializeProgramHeader
  USE UtilitiesModule, ONLY: &
    InitializeWeights
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
  USE GeometryFieldsModule, ONLY: &
    WeightsGX, &
    WeightsGX_X1, &
    WeightsGX_X2, &
    WeightsGX_X3, &
    WeightsG, &
    CreateGeometryFields, &
    DestroyGeometryFields
  USE GeometryInitializationModule, ONLY: &
    InitializeGeometry, &
    FinalizeGeometry
  USE FluidFieldsModule, ONLY: &
    WeightsF, &
    WeightsF_X1, &
    WeightsF_X2, &
    WeightsF_X3, &
    CreateFluidFields, &
    DestroyFluidFields
  USE RadiationFieldsModule, ONLY: &
    WeightsR, &
    CreateRadiationFields, &
    DestroyRadiationFields
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE OpacityModule, ONLY: &
    InitializeOpacities, &
    FinalizeOpacities
  USE GravitySolutionModule, ONLY: &
    InitializeGravitySolver, &
    FinalizeGravitySolver
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

  INCLUDE 'mpif.h'

  INTEGER :: mpierr

  PUBLIC :: InitializeProgram
  PUBLIC :: FinalizeProgram

CONTAINS


  SUBROUTINE InitializeProgram &
    ( ProgramName_Option, nX_Option, swX_Option, bcX_Option, xL_Option, &
      xR_Option, zoomX_Option, nE_Option, swE_Option, bcE_Option, &
      eL_Option, eR_Option, zoomE_Option, CoordinateSystem_Option, &
      nNodes_Option, ActivateUnits_Option, EquationOfState_Option, &
      EquationOfStateTableName_Option, Gamma_IDEAL_Option, &
      Opacity_Option, OpacityTableName_Option, GravitySolver_Option, &
      PointMass_Option, FluidSolver_Option, RadiationSolver_Option, &
      FluidRiemannSolver_Option, RadiationRiemannSolver_Option, &
      FluidRadiationCoupling_Option, EvolveGravity_Option, &
      EvolveFluid_Option, EvolveRadiation_Option, &
      ApplySlopeLimiter_Option, BetaTVB_Option, BetaTVD_Option, &
      ApplyPositivityLimiter_Option, nStages_SSP_RK_Option, &
      nStages_SI_RK_Option )

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
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: GravitySolver_Option
    REAL(DP),               INTENT(in), OPTIONAL :: PointMass_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: FluidSolver_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: RadiationSolver_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: FluidRiemannSolver_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: &
      RadiationRiemannSolver_Option
    CHARACTER(LEN=*),       INTENT(in), OPTIONAL :: &
      FluidRadiationCoupling_Option
    LOGICAL,                INTENT(in), OPTIONAL :: EvolveGravity_Option
    LOGICAL,                INTENT(in), OPTIONAL :: EvolveFluid_Option
    LOGICAL,                INTENT(in), OPTIONAL :: EvolveRadiation_Option
    LOGICAL,                INTENT(in), OPTIONAL :: ApplySlopeLimiter_Option
    REAL(DP),               INTENT(in), OPTIONAL :: BetaTVB_Option
    REAL(DP),               INTENT(in), OPTIONAL :: BetaTVD_Option
    LOGICAL,                INTENT(in), OPTIONAL :: &
      ApplyPositivityLimiter_Option
    INTEGER,                INTENT(in), OPTIONAL :: nStages_SSP_RK_Option
    INTEGER,                INTENT(in), OPTIONAL :: nStages_SI_RK_Option

    LOGICAL :: ActivateUnits
    INTEGER :: iDim

    CALL MPI_INIT( mpierr )

    CALL InitializeProgramHeader     &
           ( ProgramName_Option      &
               = ProgramName_Option, &
             nX_Option               &
               = nX_Option,          &
             swX_Option              &
               = swX_Option,         &
             bcX_Option              &
               = bcX_Option,         &
             xL_Option               &
               = xL_Option,          &
             xR_Option               &
               = xR_Option,          &
             zoomX_Option            &
               = zoomX_Option,       &
             nE_Option               &
               = nE_Option,          &
             swE_Option              &
               = swE_Option,         &
             bcE_Option              &
               = bcE_Option,         &
             eL_Option               &
               = eL_Option,          &
             eR_Option               &
               = eR_Option,          &
             zoomE_Option            &
               = zoomE_Option,       &
             nNodes_Option           &
               = nNodes_Option )

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

    ASSOCIATE( U => UnitsDisplay )

    ! --- Spatial Grid ---

    WRITE(*,'(A5,A)') '', 'Computational Domain'
    WRITE(*,'(A5,A)') '', '--------------------'
    WRITE(*,*)
    WRITE(*,'(A7,A)') '', 'Spatial Domain:'
    WRITE(*,'(A7,A)') '', '---------------'
    WRITE(*,*)
    DO iDim = 1, 3

      WRITE(*,'(A7,A3,I1,A4,ES8.2E2,A1,A,A2,A3,I1,A4,&
                &ES8.2E2,A1,A,A2,A6,I1,A4,ES10.4E2)') &
        '',   'xL(', iDim, ') = ', xL(iDim) / U % LengthUnit, &
        '', TRIM( U % LengthLabel ), &
        ', ', 'xR(', iDim, ') = ', xR(iDim) / U % LengthUnit, &
        '', TRIM( U % LengthLabel ), &
        ', ', 'ZoomX(', iDim, ') = ', ZoomX(iDim)
      WRITE(*,*)

      CALL CreateMesh &
             ( MeshX(iDim), nX(iDim), nNodesX(iDim), swX(iDim), &
               xL(iDim), xR(iDim), ZoomOption = ZoomX(iDim) )

      WRITE(*,'(A9,A11,I1,A4,ES8.2E2,A1,A,A3,ES8.2E2,A1,A)') &
        '', 'MIN/MAX dx(', iDim, ') = ', &
        MINVAL( MeshX(iDim) % Width(1:nX(iDim)) ) &
          / U % LengthUnit, '', TRIM( U % LengthLabel ), &
        ' / ', &
        MAXVAL( MeshX(iDim) % Width(1:nX(iDim)) ) &
          / U % LengthUnit, '', TRIM( U % LengthLabel )
      WRITE(*,*)

    END DO

    ! --- Spectral Grid ---

    WRITE(*,*)
    WRITE(*,'(A7,A)') '', 'Spectral Domain:'
    WRITE(*,'(A7,A)') '', '----------------'
    WRITE(*,*)

    WRITE(*,'(A7,A5,ES8.2E2,A1,A,A2,A5,ES8.2E2,A1,A,A2,A8,ES10.4E2)') &
      '', 'eL = ', eL / U % EnergyUnit, '', TRIM( U % EnergyLabel ), &
      ', ', &
          'eR = ', eR / U % EnergyUnit, '', TRIM( U % EnergyLabel ), &
      ', ', 'ZoomE = ', ZoomE
    WRITE(*,*)

    CALL CreateMesh &
           ( MeshE, nE, nNodesE, swE, eL, eR, ZoomOption = ZoomE )

    WRITE(*,'(A9,A13,ES8.2E2,A1,A,A3,ES8.2E2,A1,A)') &
      '', 'MIN/MAX de = ', &
      MINVAL( MeshE % Width(1:nE) ) / U % EnergyUnit, &
      '', TRIM( U % EnergyLabel ), &
      ' / ', &
      MAXVAL( MeshE % Width(1:nE) ) / U % EnergyUnit, &
      '', TRIM( U % EnergyLabel )

    END ASSOCIATE ! U

    ! --- Geometry ---

    CALL CreateGeometryFields( nX, swX, nE, swE )

    CALL InitializeWeights & ! --- For Integration in Elements
           ( MeshX(1) % Weights, &
             MeshX(2) % Weights, &
             MeshX(3) % Weights, &
             WeightsGX, WeightsGX_X1, WeightsGX_X2, WeightsGX_X3 )

    CALL InitializeWeights & ! --- For Integration in Elements
           ( MeshE    % Weights, MeshX(1) % Weights, &
             MeshX(2) % Weights, MeshX(3) % Weights, &
             WeightsG )

    CALL InitializeGeometry &
           ( nX, nNodesX, swX, nE, nNodesE, swE, &
             CoordinateSystem_Option &
               = CoordinateSystem_Option )

    ! --- Physical Fields ---

    ! --- Fluid Fields ---

    CALL CreateFluidFields( nX, swX )

    CALL InitializeWeights & ! --- For Integration in Elements
           ( MeshX(1) % Weights, &
             MeshX(2) % Weights, &
             MeshX(3) % Weights, &
             WeightsF, WeightsF_X1, WeightsF_X2, WeightsF_X3 )

    ! --- Radiation Fields ---

    CALL CreateRadiationFields( nX, swX, nE, swE )

    CALL InitializeWeights & ! --- For Integration in Elements
           ( MeshE    % Weights, MeshX(1) % Weights, &
             MeshX(2) % Weights, MeshX(3) % Weights, &
             WeightsR )

    ! --- For Mapping Between Nodal and Modal Representations ---

    CALL InitializePolynomialBasisMapping &
           ( MeshE    % Nodes, MeshX(1) % Nodes, &
             MeshX(2) % Nodes, MeshX(3) % Nodes )

    ! --- Equation of State ---

    CALL InitializeEquationOfState &
           ( EquationOfState_Option &
               = EquationOfState_Option, &
             EquationOfStateTableName_Option &
               = EquationOfStateTableName_Option, &
             Gamma_IDEAL_Option &
               = Gamma_IDEAL_Option )

    ! --- Opacities ---

    CALL InitializeOpacities &
           ( Opacity_Option &
               = Opacity_Option, &
             OpacityTableName_Option &
               = OpacityTableName_Option )

    ! --- Gravity Solver ---

    CALL InitializeGravitySolver &
           ( GravitySolver_Option &
               = GravitySolver_Option, &
             PointMass_Option &
               = PointMass_Option )

    ! --- Fluid Solver ---

    CALL InitializeFluidEvolution &
           ( FluidSolver_Option &
               = FluidSolver_Option, &
             ApplySlopeLimiter_Option &
               = ApplySlopeLimiter_Option, &
             BetaTVB_Option &
               = BetaTVB_Option, &
             BetaTVD_Option &
               = BetaTVD_Option, &
             ApplyPositivityLimiter_Option &
               = ApplyPositivityLimiter_Option )

    ! --- Radiation Solver ---

    CALL InitializeRadiationEvolution &
           ( RadiationSolver_Option &
               = RadiationSolver_Option, &
             ApplySlopeLimiter_Option &
               = ApplySlopeLimiter_Option, &
             BetaTVB_Option &
               = BetaTVB_Option, &
             BetaTVD_Option &
               = BetaTVD_Option, &
             ApplyPositivityLimiter_Option &
               = ApplyPositivityLimiter_Option )

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
           ( EvolveGravity_Option &
               = EvolveGravity_Option, &
             EvolveFluid_Option &
               = EvolveFluid_Option, &
             EvolveRadiation_Option &
               = EvolveRadiation_Option, &
             nStages_SSP_RK_Option &
               = nStages_SSP_RK_Option, &
             nStages_SI_RK_Option &
               = nStages_SI_RK_Option )

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

    CALL DestroyGeometryFields

    CALL FinalizeGeometry

    ! --- Physical Fields ---

    CALL DestroyFluidFields

    CALL DestroyRadiationFields

    ! --- Equation of State ---

    CALL FinalizeEquationOfState

    ! --- Opacities ---

    CALL FinalizeOpacities

    ! --- Gravity Solver ---

    CALL FinalizeGravitySolver

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

    CALL MPI_FINALIZE( mpierr )

  END SUBROUTINE FinalizeProgram


END MODULE ProgramInitializationModule
