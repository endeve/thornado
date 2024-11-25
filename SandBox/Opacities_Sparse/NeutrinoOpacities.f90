PROGRAM NeutrinoOpacities
  USE TasmanianSG
  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE MeshModule, ONLY: &
    MeshE, &
    NodeCoordinate
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC_Point, &
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Point, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Point, &
    ComputeNeutrinoOpacities_NES_Points, &
    ComputeNeutrinoOpacitiesAndDerivatives_NES_Point, &
    ComputeNeutrinoOpacities_Pair_Point, &
    ComputeNeutrinoOpacities_Pair_Points, &
    ComputeNeutrinoOpacitiesAndDerivatives_Pair_Point
  USE DeviceModule, ONLY: &
    InitializeDevice, &
    FinalizeDevice
  USE NeutrinoOpacitiesSparseComputationModule, ONLY: &
    ComputeNeutrinoOpacities_NES_Points_SG

  IMPLICIT NONE

  INCLUDE 'mpif.h'


  INTEGER, PARAMETER :: &
    nNodes   = 2, &
    nE       = 2**4, &
    nX1      = 2**0, &
    nPointsX = nNodes * nX1, &
    nPointsE = nNodes * nE, &
    nSpecies = 2
  REAL(DP), PARAMETER :: &
    Unit_D     = Gram / Centimeter**3, &
    Unit_T     = Kelvin, &
    Unit_Y     = 1.0_DP, &
    Unit_E     = MeV, &
    Unit_Chi   = 1.0_DP / Centimeter, &
    Unit_Sigma = 1.0_DP / Centimeter, &
    eL         = 0.0e0_DP * Unit_E, &
    eR         = 3.0e2_DP * Unit_E, &
    ZoomE      = 1.183081754893913_DP

  INTEGER :: &
    mpierr, iE, iX, iS, iNodeE, iN_E
  REAL(DP) :: &
    Timer_ReadEos, &
    Timer_ReadOpacities, &
    Timer_Compute_EC, &
    Timer_Compute_EC_Point, &
    Timer_Compute_ES, &
    Timer_Compute_ES_Point, &
    Timer_Compute_NES, &
    Timer_Compute_NES_Point, &
    Timer_Compute_NES_D_Point, &
    Timer_Compute_Pair, &
    Timer_Compute_Pair_Point, &
    Timer_Compute_Pair_D_Point, &
    Timer_Total
  REAL(DP), DIMENSION(nPointsX) :: &
    D, T, Y
  REAL(DP), DIMENSION(nPointsE) :: &
    E, dE
  REAL(DP), DIMENSION(nPointsE,nPointsX,nSpecies) :: &
    Chi, &     ! --- Absorption Opacity
    Sigma, &   ! --- Scattering Opacity (Isoenergetic)
    Chi_NES, & ! --- Integrated NES Opacity
    Chi_Pair   ! --- Integrated Pair Opacity
  REAL(DP), DIMENSION(nPointsE,nPointsE,nPointsX) :: &
    Phi_0_NES_In_1, Phi_0_NES_In_2, Phi_0_NES_Out_1, Phi_0_NES_Out_2, &
    Phi_0_Pair_In_1, Phi_0_Pair_In_2, Phi_0_Pair_Out_1, Phi_0_Pair_Out_2, &
    Phi_0_NES_In_1_SG, Phi_0_NES_In_2_SG, Phi_0_NES_Out_1_SG, Phi_0_NES_Out_2_SG, &
    Phi_0_Pair_In_1_SG, Phi_0_Pair_In_2_SG, Phi_0_Pair_Out_1_SG, Phi_0_Pair_Out_2_SG

  type(TasmanianSparseGrid) :: gridNES, grid1, grid2, grid3
  INTEGER :: dims, outs, level
  INTEGER :: N1, N2, N3
  INTEGER :: N, i, j, verm, vern
  DOUBLE PRECISION :: err1, err2, err3, err4, exact
  REAL :: cpuStart, cpuEnd, stages(2,3)
  INTEGER :: conformal(3)
  CHARACTER, pointer :: string(:)
  DOUBLE PRECISION, pointer :: points(:,:), weights(:)
  DOUBLE PRECISION :: x, integ
  DOUBLE PRECISION, allocatable :: transformA(:), transformB(:), values(:,:), tvalues(:,:)
  DOUBLE PRECISION, allocatable :: res(:), res2(:,:)
  DOUBLE PRECISION :: desired_x(2)
  DOUBLE PRECISION :: randPoints(4,1000), randPoints2(2,1000), randPoints3(3,1000)
  DOUBLE PRECISION :: PI = 4.D0*DATAN(1.D0)
  CHARACTER (LEN=40) :: ReducedNESGridName
  logical :: SGReadError
  DOUBLE PRECISION :: SGerrorNES

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'NeutrinoOpacities', &
           nX_Option &
             = [ nX1, 1, 1 ], &
           swX_Option &
             = [ 01, 00, 00 ], &
           bcX_Option &
             = [ 32, 00, 00 ], &
           nE_Option &
             = nE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           ZoomE_Option &
             = ZoomE, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'CARTESIAN', &
           ActivateUnits_Option &
             = .TRUE., &
           nSpecies_Option &
             = nSpecies, &
           BasicInitialization_Option &
             = .TRUE. )

  WRITE(*,*)
  WRITE(*,'(A4,A)') '', 'NeutrinoOpacities'
  WRITE(*,*)
  WRITE(*,'(A6,A,I8.8)') '', 'nPointsX = ', nPointsX
  WRITE(*,'(A6,A,I8.8)') '', 'nPointsE = ', nPointsE
  WRITE(*,*)

  ! --- Thermodynamic State ---

!  D = 1.3d14 * Unit_D
  D = 5.6d13 * Unit_D
  T = 3.0d11 * Unit_T
  Y = 2.9d-1 * Unit_Y

  ! --- Energy Grid ---

  !dE(1:3) = 2.0_DP * Unit_E
  !DO iE = 4, nPointsE
  !  dE(iE) = 1.095_DP * dE(iE-1)
  !END DO

  DO iN_E = 1, nPointsE
    iE      = MOD( (iN_E-1) / nNodes, nE     ) + 1
    iNodeE  = MOD( (iN_E-1)         , nNodes ) + 1
    E(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )
    WRITE(*,'(A6,A2,I3.3,A10,ES8.2E2)') '','E(',iN_E,') [MeV] = ', E(iN_E) / Unit_E
  END DO


  ! --- Initialize Equation of State ---

  Timer_ReadEos = MPI_WTIME()

  CALL InitializeEquationOfState_TABLE &
         ( EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5', &
           Verbose_Option = .TRUE. )

  Timer_ReadEos = MPI_WTIME() - Timer_ReadEos

  ! --- Initialize Opacities ---

  Timer_ReadOpacities = MPI_WTIME()

  ! CALL InitializeOpacities_TABLE &
  !        ( OpacityTableName_EmAb_Option &
  !            = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &
  !          OpacityTableName_Iso_Option  &
  !            = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5',  &
  !          OpacityTableName_NES_Option &
  !            = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5',  &
  !          OpacityTableName_Pair_Option &
  !            = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &
  !          Verbose_Option = .TRUE. )

  CALL InitializeOpacities_TABLE &
         ( OpacityTableName_NES_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5',  &
           Verbose_Option = .TRUE. )
  Timer_ReadOpacities = MPI_WTIME() - Timer_ReadOpacities

  ! --- Compute NES Opacities ---

  Timer_Compute_NES = MPI_WTIME()

  CALL ComputeNeutrinoOpacities_NES_Points &
         ( 1, nPointsE, 1, nPointsX, E, D, T, Y, 1, 2, 1, &
           Phi_0_NES_In_1(:,:,:), Phi_0_NES_Out_1(:,:,:), &
           Phi_0_NES_In_2(:,:,:), Phi_0_NES_Out_2(:,:,:) )

  Timer_Compute_NES = MPI_WTIME() - Timer_Compute_NES




  CALL tsgAllocateGrid(gridNES)
  ! ReducedNESGridName =  'reduced_NES_bin.grid'
  ReducedNESGridName =  'NES_log_bin.grid'
  ! ReducedNESGridName =  'reduced_NES_log_bin.grid'
  SGReadError = tsgRead(gridNES, ReducedNESGridName)


  CALL ComputeNeutrinoOpacities_NES_Points_SG &
         (gridNES, 1, nPointsE, 1, nPointsX, E, D, T, Y, 1, 2, 1, &
           Phi_0_NES_In_1_SG(:,:,:), Phi_0_NES_Out_1_SG(:,:,:), &
           Phi_0_NES_In_2_SG(:,:,:), Phi_0_NES_Out_2_SG(:,:,:) )

  SGerrorNES = MAXVAL(ABS(Phi_0_NES_In_1_SG - Phi_0_NES_In_1)) &
                / MAXVAL(ABS(Phi_0_NES_In_1))

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "NES SG ERROR = ", SGerrorNES
  WRITE(*,*) "-------------------------------------------------------------------------------------------------"

  CALL WriteMatrix &
         ( nPointsE, nPointsE, Phi_0_NES_In_1 (:,:,1), 'Phi_0_NES_In_1.dat'  )
  CALL WriteMatrix &
         ( nPointsE, nPointsE, Phi_0_NES_In_1_SG (:,:,1), 'Phi_0_NES_In_1_SG.dat'  )





  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  CALL tsgDeallocateGrid(gridNES)

END PROGRAM NeutrinoOpacities
