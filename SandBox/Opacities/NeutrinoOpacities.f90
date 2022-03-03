PROGRAM NeutrinoOpacities

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV, &
    BoltzmannConstant
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
    FinalizeOpacities_TABLE, &
    C1, C2
  USE RadiationFieldsModule, ONLY: &
    iNuE, iNuE_Bar
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC, &
    ComputeNeutrinoOpacities_ES, &
    ComputeNeutrinoOpacities_NES, &
    ComputeNeutrinoOpacities_Pair
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
    UnitNES    = 1.0_DP / ( Centimeter * MeV**3 ), &
    UnitPair   = 1.0_DP / ( Centimeter * MeV**3 ), &
    eL         = 0.0e0_DP * Unit_E, &
    eR         = 3.0e2_DP * Unit_E, &
    ZoomE      = 1.183081754893913_DP

  INTEGER :: &
    mpierr, iE, iX, iS, iNodeE, iN_E, iE1, iE2
  REAL(DP) :: &
    kT, DetBal, &
    Timer_ReadEos, &
    Timer_ReadOpacities, &
    Timer_Compute_EC, &
    Timer_Compute_ES, &
    Timer_Compute_NES, &
    Timer_Compute_Pair, &
    Timer_Total
  REAL(DP), DIMENSION(nPointsX) :: &
    D, T, Y
  REAL(DP), DIMENSION(nPointsE) :: &
    E, dE, &
    Phi_0_In, Phi_0_Out
  REAL(DP), DIMENSION(nPointsE,nPointsX,nSpecies) :: &
    Chi, &     ! --- Absorption Opacity
    Sigma, &   ! --- Scattering Opacity (Isoenergetic)
    Chi_NES, & ! --- Integrated NES Opacity
    Chi_Pair   ! --- Integrated Pair Opacity
  REAL(DP), DIMENSION(nPointsE,nPointsE,nPointsX) :: &
    H1, H2, &  ! --- NES  Scattering Functions
    J1, J2     ! --- Pair Scattering Functions

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
  !$OMP MAP( alloc: Chi, Chi_NES, Chi_Pair, Sigma, &
  !$OMP             H1, H2, J1, J2 )
#elif defined(THORNADO_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( E, D, T, Y ) &
  !$ACC CREATE( Chi, Chi_NES, Chi_Pair, Sigma, &
  !$ACC         H1, H2, J1, J2 )
#endif

  ! --- Compute Electron Capture Opacities ---

  Timer_Compute_EC = MPI_WTIME()

  CALL ComputeNeutrinoOpacities_EC &
         ( 1, nPointsE, 1, nPointsX, 1, nSpecies, E, D, T, Y, Chi )

  Timer_Compute_EC = MPI_WTIME() - Timer_Compute_EC

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( Chi )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( Chi )
#endif

  ! --- Compute Elastic Scattering Opacities ---

  Timer_Compute_ES = MPI_WTIME()

  CALL ComputeNeutrinoOpacities_ES &
         ( 1, nPointsE, 1, nPointsX, 1, nSpecies, E, D, T, Y, 1, Sigma )

  Timer_Compute_ES = MPI_WTIME() - Timer_Compute_ES

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( Sigma )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( Sigma )
#endif

  ! --- Compute NES Opacities ---

  Timer_Compute_NES = MPI_WTIME()

  CALL ComputeNeutrinoOpacities_NES &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 1, H1, H2 )

  Timer_Compute_NES = MPI_WTIME() - Timer_Compute_NES

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( H1, H2 )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( H1, H2 )
#endif

  ! --- Integrated NES Opacity ---

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    kT = BoltzmannConstant * T(iX)

    DO iE2 = 1, nPointsE

      DO iE1 = 1, nPointsE
        DetBal = EXP( - ABS( E(iE2) - E(iE1) ) / kT )
        IF ( iE1 <= iE2 ) THEN
          Phi_0_Out(iE1) = ( C1(iS) * H1(iE1,iE2,iX) + C2(iS) * H2(iE1,iE2,iX) ) * UnitNES
          Phi_0_In (iE1) = Phi_0_Out(iE1) * DetBal
        ELSE
          Phi_0_In (iE1) = ( C1(iS) * H1(iE2,iE1,iX) + C2(iS) * H2(iE2,iE1,iX) ) * UnitNES
          Phi_0_Out(iE1) = Phi_0_In (iE1) * DetBal
        END IF
      END DO

      Chi_NES(iE2,iX,iS) &
        = TRAPEZ( nPointsE, E, Phi_0_Out * E**2 )

    END DO

  END DO
  END DO

  ! --- Compute Pair Opacities ---

  Timer_Compute_Pair = MPI_WTIME()

  CALL ComputeNeutrinoOpacities_Pair &
         ( 1, nPointsE, 1, nPointsX, D, T, Y, 1, J1, J2 )

  Timer_Compute_Pair = MPI_WTIME() - Timer_Compute_Pair

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET UPDATE FROM( J1, J2 )
#elif defined(THORNADO_OACC)
  !$ACC UPDATE HOST( J1, J2 )
#endif

  ! --- Integrated Pair Opacity ---

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    kT = BoltzmannConstant * T(iX)

    DO iE2 = 1, nPointsE

      DO iE1 = 1, nPointsE
        DetBal = EXP( - ABS( E(iE2) + E(iE1) ) / kT )
        IF ( iE1 <= iE2 ) THEN
          Phi_0_Out(iE1) = ( C1(iS) * J1(iE1,iE2,iX) + C2(iS) * J2(iE1,iE2,iX) ) * UnitPair
          Phi_0_In (iE1) = Phi_0_Out(iE1) * DetBal
        ELSE
          Phi_0_In (iE1) = ( C1(iS) * J2(iE2,iE1,iX) + C2(iS) * J1(iE2,iE1,iX) ) * UnitPair
          Phi_0_Out(iE1) = Phi_0_In (iE1) * DetBal
        END IF
      END DO

      Chi_Pair(iE2,iX,iS) &
        = TRAPEZ( nPointsE, E, Phi_0_In * E**2 )

    END DO

  END DO
  END DO

#if defined(THORNADO_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: E, D, T, Y, &
  !$OMP               Chi, Chi_NES, Chi_Pair, Sigma, &
  !$OMP               H1, H2, J1, J2 )
#elif defined(THORNADO_OACC)
  !$ACC EXIT DATA &
  !$ACC DELETE( E, D, T, Y, &
  !$ACC         Chi, Chi_NES, Chi_Pair, Sigma, &
  !$ACC         H1, H2, J1, J2 )
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
         ( nPointsE, nPointsE, H1(:,:,1), 'H1.dat'  )

  CALL WriteMatrix &
         ( nPointsE, nPointsE, J1(:,:,1), 'J1.dat' )

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  Timer_Total &
    = Timer_Compute_EC + Timer_Compute_ES &
      + Timer_Compute_NES + Timer_Compute_Pair

  WRITE(*,*)
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'ReadEos = ',       &
    Timer_ReadEos
  WRITE(*,'(A4,A22,1ES10.2E2)') '', 'ReadOpacities = ', &
    Timer_ReadOpacities
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_EC = ',    &
    Timer_Compute_EC, Timer_Compute_EC / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_ES = ',    &
    Timer_Compute_ES, Timer_Compute_ES / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_NES = ',   &
    Timer_Compute_NES, Timer_Compute_NES / Timer_Total
  WRITE(*,'(A4,A22,2ES10.2E2)') '', 'Compute_Pair = ',  &
    Timer_Compute_Pair, Timer_Compute_Pair / Timer_Total
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
