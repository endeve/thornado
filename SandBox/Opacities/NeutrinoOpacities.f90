PROGRAM NeutrinoOpacities

  USE KindModule, ONLY: &
    DP, FourPi
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

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: &
    nNodes   = 2, &
    nE       = 2**4, &
    nX1      = 2**11, &
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
    Sigma1, &  !
    Sigma2, &  !
    Chi_NES, & ! --- Integrated NES Opacity
    Chi_Pair   ! --- Integrated Pair Opacity
  REAL(DP), DIMENSION(nPointsE,nPointsE,nPointsX,nSpecies) :: &
    Phi_0_NES_In,   dPhi_0_NES_In_dY,   dPhi_0_NES_In_dE, &
    Phi_0_NES_Out,  dPhi_0_NES_Out_dY,  dPhi_0_NES_Out_dE, &
    Phi_0_Pair_In,  dPhi_0_Pair_In_dY,  dPhi_0_Pair_In_dE, &
    Phi_0_Pair_Out, dPhi_0_Pair_Out_dY, dPhi_0_Pair_Out_dE

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

  CALL InitializeOpacities_TABLE &
         ( OpacityTableName_EmAb_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &
           OpacityTableName_Iso_Option  &
             = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5',  &
           OpacityTableName_NES_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5',  &
           OpacityTableName_Pair_Option &
             = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &
           Verbose_Option = .TRUE. )

  Timer_ReadOpacities = MPI_WTIME() - Timer_ReadOpacities

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: E, D, T, Y ) &
  !$OMP MAP( alloc: J0, Chi, Chi_NES, Chi_Pair, Chi_Brem, Sigma, Sigma1, Sigma2, &
  !$OMP             Phi_0_NES_In, Phi_0_NES_Out, Phi_0_Pair_In, Phi_0_Pair_Out )
#elif defined(THORNADO_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( E, D, T, Y ) &
  !$ACC CREATE( J0, Chi, Chi_NES, Chi_Pair, Chi_Brem, Sigma, Sigma1, Sigma2, &
  !$ACC         Phi_0_NES_In, Phi_0_NES_Out, Phi_0_Pair_In, Phi_0_Pair_Out )
#endif

  ! --- Compute Electron Capture Opacities ---

  Timer_Compute_EC = MPI_WTIME()

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_EC_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, Chi(:,:,iS) )

  END DO

  Timer_Compute_EC = MPI_WTIME() - Timer_Compute_EC

  ! --- Compute Electron Capture Opacities (Point) ---

  Timer_Compute_EC_Point = MPI_WTIME()

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
  !$ACC PRESENT( E, D, T, Y, Chi )
#endif
  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacities_EC_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, Chi(:,iX,iS) )

  END DO
  END DO

  Timer_Compute_EC_Point = MPI_WTIME() - Timer_Compute_EC_Point

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( Chi )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( Chi )
#endif

  ! --- Compute Elastic Scattering Opacities ---

  Timer_Compute_ES = MPI_WTIME()

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 1, Sigma1(:,:,iS) )

    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 2, Sigma2(:,:,iS) )

  END DO

  Timer_Compute_ES = MPI_WTIME() - Timer_Compute_ES

  ! --- Compute Elastic Scattering Opacities (Point) ---

  Timer_Compute_ES_Point = MPI_WTIME()

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
  !$ACC PRESENT( E, D, T, Y, Sigma, Sigma1, Sigma2 )
#endif
  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacities_ES_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 1, Sigma1(:,iX,iS) )

    CALL ComputeNeutrinoOpacities_ES_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 2, Sigma2(:,iX,iS) )

  END DO
  END DO

  DO iE = 1, nPointsE
    Sigma(iE,:,:) = FourPi * E(iE)**2 * ( Sigma1(iE,:,:) - Sigma2(iE,:,:) / 3.0d0 )
                    ! (A41) in Bruenn 85
  END DO

  Timer_Compute_ES_Point = MPI_WTIME() - Timer_Compute_ES_Point

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( Sigma, Sigma1, Sigma2 )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( Sigma, Sigma1, Sigma2 )
#endif

  ! --- Compute NES Opacities ---

  Timer_Compute_NES = MPI_WTIME()

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 1, &
             Phi_0_NES_In(:,:,:,iS), Phi_0_NES_Out(:,:,:,iS) )

  END DO

  Timer_Compute_NES = MPI_WTIME() - Timer_Compute_NES

  ! --- Compute NES Opacities (Point) ---

  Timer_Compute_NES_Point = MPI_WTIME()

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
  !$ACC PRESENT( E, D, T, Y, Phi_0_NES_In, Phi_0_NES_Out )
#endif
  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacities_NES_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 1, &
             Phi_0_NES_In(:,:,iX,iS), Phi_0_NES_Out(:,:,iX,iS) )

  END DO
  END DO

  Timer_Compute_NES_Point = MPI_WTIME() - Timer_Compute_NES_Point

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( Phi_0_NES_In, Phi_0_NES_Out )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( Phi_0_NES_In, Phi_0_NES_Out )
#endif

  ! --- Compute NES Opacities and Derivatives (Point) ---

  Timer_Compute_NES_D_Point = MPI_WTIME()

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacitiesAndDerivatives_NES_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 1, &
             Phi_0_NES_In     (:,:,iX,iS), &
             Phi_0_NES_Out    (:,:,iX,iS), &
             dPhi_0_NES_In_dY (:,:,iX,iS), &
             dPhi_0_NES_In_dE (:,:,iX,iS), &
             dPhi_0_NES_Out_dY(:,:,iX,iS), &
             dPhi_0_NES_Out_dE(:,:,iX,iS) )

  END DO
  END DO

  Timer_Compute_NES_D_Point = MPI_WTIME() - Timer_Compute_NES_D_Point

  ! --- Integrated NES Opacity ---

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX
  DO iE = 1, nPointsE

    Chi_NES(iE,iX,iS) &
      = TRAPEZ( nPointsE, E, Phi_0_NES_Out(:,iE,iX,iS) * E**2 )

  END DO
  END DO
  END DO

  Chi_NES = Chi_NES * FourPi ! (A38) in Bruenn85

  ! --- Compute Pair Opacities ---

  Timer_Compute_Pair = MPI_WTIME()

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 1, &
             Phi_0_Pair_In(:,:,:,iS), Phi_0_Pair_Out(:,:,:,iS) )

  END DO

  Timer_Compute_Pair = MPI_WTIME() - Timer_Compute_Pair

  ! --- Compute Pair Opacities (Point) ---

  Timer_Compute_Pair_Point = MPI_WTIME()

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
  !$ACC PRESENT( E, D, T, Y, Phi_0_Pair_In, Phi_0_Pair_Out )
#endif
  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacities_Pair_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 1, &
             Phi_0_Pair_In(:,:,iX,iS), Phi_0_Pair_Out(:,:,iX,iS) )

  END DO
  END DO

  Timer_Compute_Pair_Point = MPI_WTIME() - Timer_Compute_Pair_Point

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( Phi_0_Pair_In, Phi_0_Pair_Out )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( Phi_0_Pair_In, Phi_0_Pair_Out )
#endif

  ! --- Compute Pair Opacities and Derivatives (Point) ---

  Timer_Compute_Pair_D_Point = MPI_WTIME()

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacitiesAndDerivatives_Pair_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 1, &
             Phi_0_Pair_In     (:,:,iX,iS), &
             Phi_0_Pair_Out    (:,:,iX,iS), &
             dPhi_0_Pair_In_dY (:,:,iX,iS), &
             dPhi_0_Pair_In_dE (:,:,iX,iS), &
             dPhi_0_Pair_Out_dY(:,:,iX,iS), &
             dPhi_0_Pair_Out_dE(:,:,iX,iS) )

  END DO
  END DO

  Timer_Compute_Pair_D_Point = MPI_WTIME() - Timer_Compute_Pair_D_Point

  ! --- Integrated Pair Opacity ---

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX
  DO iE = 1, nPointsE

    Chi_Pair(iE,iX,iS) &
      = TRAPEZ( nPointsE, E, Phi_0_Pair_In(:,iE,iX,iS) * E**2 )

  END DO
  END DO
  END DO

  Chi_Pair  = Chi_Pair * FourPi ! (A47) in Bruenn85

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: E, D, T, Y, &
  !$OMP               Chi, Chi_NES, Chi_Pair, Sigma, Sigma1, Sigma2, &
  !$OMP               Phi_0_NES_In, Phi_0_NES_Out, Phi_0_Pair_In, Phi_0_Pair_Out )
#elif defined(THORNADO_OACC)
  !$ACC EXIT DATA &
  !$ACC DELETE( E, D, T, Y, &
  !$ACC         Chi, Chi_NES, Chi_Pair, Sigma, Sigma1, Sigma2, &
  !$ACC         Phi_0_NES_In, Phi_0_NES_Out, Phi_0_Pair_In, Phi_0_Pair_Out )
#endif

  CALL WriteVector &
         ( nPointsE, E / Unit_E, 'E.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi(:,1,iNuE) / Unit_Chi, 'Chi_NuE.dat' )

  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi(:,1,iNuE_Bar) / Unit_Chi, 'Chi_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Sigma(:,1,iNuE) / Unit_Sigma, 'Sigma_NuE.dat' )

  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Sigma(:,1,iNuE_Bar) / Unit_Sigma, 'Sigma_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi_NES(:,1,iNuE) / Unit_Chi, 'Chi_NES_NuE.dat' )

  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_NES(:,1,iNuE_Bar) / Unit_Chi, 'Chi_NES_NuE_Bar.dat' )

  CALL WriteVector & ! --- NuE
         ( nPointsE, Chi_Pair(:,1,iNuE) / Unit_Chi,  'Chi_Pair_NuE.dat' )

  CALL WriteVector & ! --- NuE_Bar
         ( nPointsE, Chi_Pair(:,1,iNuE_Bar) / Unit_Chi, 'Chi_Pair_NuE_Bar.dat' )

  CALL WriteMatrix &
         ( nPointsE, nPointsE, Phi_0_NES_Out (:,:,1,1), 'Phi_0_NES_Out.dat'  )

  CALL WriteMatrix &
         ( nPointsE, nPointsE, Phi_0_Pair_Out(:,:,1,1), 'Phi_0_Pair_Out.dat' )

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  Timer_Total &
    = Timer_Compute_EC + Timer_Compute_EC_Point &
      + Timer_Compute_ES + Timer_Compute_ES_Point &
      + Timer_Compute_NES + Timer_Compute_NES_Point &
      + Timer_Compute_NES_D_Point + Timer_Compute_Pair &
      + Timer_Compute_Pair_Point + Timer_Compute_Pair_D_Point

  WRITE(*,*)
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'ReadEos = ',       &
    Timer_ReadEos
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'ReadOpacities = ', &
    Timer_ReadOpacities
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_EC = ',    &
    Timer_Compute_EC, Timer_Compute_EC / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_EC (P) = ',    &
    Timer_Compute_EC_Point, Timer_Compute_EC_Point / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_ES = ',    &
    Timer_Compute_ES, Timer_Compute_ES / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_ES (P) = ',    &
    Timer_Compute_ES_Point, Timer_Compute_ES_Point / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_NES = ',   &
    Timer_Compute_NES, Timer_Compute_NES / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_NES (P) = ',   &
    Timer_Compute_NES_Point, Timer_Compute_NES_Point / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_NES_D (P) = ',   &
    Timer_Compute_NES_D_Point, Timer_Compute_NES_D_Point / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_Pair = ',  &
    Timer_Compute_Pair, Timer_Compute_Pair / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_Pair (P) = ',  &
    Timer_Compute_Pair_Point, Timer_Compute_Pair_Point / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_Pair_D (P) = ',   &
    Timer_Compute_Pair_D_Point, Timer_Compute_Pair_D_Point / Timer_Total
  WRITE(*,*)

  CALL FinalizeDevice

  CALL MPI_FINALIZE( mpierr )

CONTAINS


  PURE REAL(DP) FUNCTION TRAPEZ( n, x, y )

    INTEGER,  INTENT(in) :: n
    REAL(DP), INTENT(in) :: x(n), y(n)

    INTEGER :: i

    TRAPEZ = 0.0_DP
    DO i = 1, n - 1
      TRAPEZ = TRAPEZ + 0.5_dp * ( x(i+1) - x(i) ) * ( y(i) + y(i+1) )
    END DO

    RETURN
  END FUNCTION TRAPEZ


END PROGRAM NeutrinoOpacities
