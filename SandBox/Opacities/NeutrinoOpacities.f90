PROGRAM NeutrinoOpacities

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV
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
    ComputeNeutrinoOpacities_Pair_Point, &
    ComputeNeutrinoOpacities_Pair_Points

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: &
    nPointsX = 2**4, &
    nPointsE = 2**5, &
    nSpecies = 2
  REAL(DP), PARAMETER :: &
    Unit_D     = Gram / Centimeter**3, &
    Unit_T     = Kelvin, &
    Unit_Y     = 1.0_DP, &
    Unit_E     = MeV, &
    Unit_Chi   = 1.0_DP / Centimeter, &
    Unit_Sigma = 1.0_DP / Centimeter

  INTEGER :: &
    mpierr, iE, iX, iS
  REAL(DP) :: &
    Timer_ReadEos, &
    Timer_ReadOpacities, &
    Timer_Compute_EC, &
    Timer_Compute_ES, &
    Timer_Compute_NES, &
    Timer_Compute_Pair, &
    Timer_Compute_EC_Point, &
    Timer_Compute_ES_Point, &
    Timer_Compute_NES_Point, &
    Timer_Compute_Pair_Point, &
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
  REAL(DP), DIMENSION(nPointsE,nPointsE,nPointsX,nSpecies) :: &
    Phi_0_NES_In, &
    Phi_0_NES_Out, &
    Phi_0_Pair_In, &
    Phi_0_Pair_Out

  CALL MPI_INIT( mpierr )

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

  dE(1:3) = 2.0_DP * Unit_E
  DO iE = 4, nPointsE
    dE(iE) = 1.095_DP * dE(iE-1)
  END DO

  DO iE = 1, nPointsE
    E(iE) = SUM( dE(1:iE-1) ) + 0.5_DP * dE(iE)
    WRITE(*,'(A6,A2,I3.3,A10,ES8.2E2)') '','E(',iE,') [MeV] = ', E(iE) / Unit_E
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

  ! --- Compute Electron Capture Opacities ---

  Timer_Compute_EC = MPI_WTIME()

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_EC_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, Chi(:,:,iS) )

  END DO

  Timer_Compute_EC = MPI_WTIME() - Timer_Compute_EC

  ! --- Compute Electron Capture Opacities (Point) ---

  Timer_Compute_EC_Point = MPI_WTIME()

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacities_EC_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, Chi(:,iX,iS) )

  END DO
  END DO

  Timer_Compute_EC_Point = MPI_WTIME() - Timer_Compute_EC_Point

  ! --- Compute Elastic Scattering Opacities ---

  Timer_Compute_ES = MPI_WTIME()

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 1, Sigma(:,:,iS) )

  END DO

  Timer_Compute_ES = MPI_WTIME() - Timer_Compute_ES

  ! --- Compute Elastic Scattering Opacities (Point) ---

  Timer_Compute_ES_Point = MPI_WTIME()

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacities_ES_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 1, Sigma(:,iX,iS) )

  END DO
  END DO

  Timer_Compute_ES_Point = MPI_WTIME() - Timer_Compute_ES_Point

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

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacities_NES_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 1, &
             Phi_0_NES_In(:,:,iX,iS), Phi_0_NES_Out(:,:,iX,iS) )

  END DO
  END DO

  Timer_Compute_NES_Point = MPI_WTIME() - Timer_Compute_NES_Point

  ! --- Integrated NES Opacity ---

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX
  DO iE = 1, nPointsE

    Chi_NES(iE,iX,iS) &
      = TRAPEZ( nPointsE, E, Phi_0_NES_Out(:,iE,iX,iS) * E**2 )

  END DO
  END DO
  END DO

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

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX

    CALL ComputeNeutrinoOpacities_Pair_Point &
           ( 1, nPointsE, E, D(iX), T(iX), Y(iX), iS, 1, &
             Phi_0_Pair_In(:,:,iX,iS), Phi_0_Pair_Out(:,:,iX,iS) )

  END DO
  END DO

  Timer_Compute_Pair_Point = MPI_WTIME() - Timer_Compute_Pair_Point

  ! --- Integrated Pair Opacity ---

  DO iS = 1, nSpecies
  DO iX = 1, nPointsX
  DO iE = 1, nPointsE

    Chi_Pair(iE,iX,iS) &
      = TRAPEZ( nPointsE, E, Phi_0_Pair_In(:,iE,iX,iS) * E**2 )

  END DO
  END DO
  END DO

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
      + Timer_Compute_Pair + Timer_Compute_Pair_Point

  WRITE(*,*)
  WRITE(*,'(A4,A19,1ES10.2E2)') '', 'ReadEos = ',       &
    Timer_ReadEos
  WRITE(*,'(A4,A19,1ES10.2E2)') '', 'ReadOpacities = ', &
    Timer_ReadOpacities
  WRITE(*,'(A4,A19,2ES10.2E2)') '', 'Compute_EC = ',    &
    Timer_Compute_EC, Timer_Compute_EC / Timer_Total
  WRITE(*,'(A4,A19,2ES10.2E2)') '', 'Compute_EC (P) = ',    &
    Timer_Compute_EC_Point, Timer_Compute_EC_Point / Timer_Total
  WRITE(*,'(A4,A19,2ES10.2E2)') '', 'Compute_ES = ',    &
    Timer_Compute_ES, Timer_Compute_ES / Timer_Total
  WRITE(*,'(A4,A19,2ES10.2E2)') '', 'Compute_ES (P) = ',    &
    Timer_Compute_ES_Point, Timer_Compute_ES_Point / Timer_Total
  WRITE(*,'(A4,A19,2ES10.2E2)') '', 'Compute_NES = ',   &
    Timer_Compute_NES, Timer_Compute_NES / Timer_Total
  WRITE(*,'(A4,A19,2ES10.2E2)') '', 'Compute_NES (P) = ',   &
    Timer_Compute_NES_Point, Timer_Compute_NES_Point / Timer_Total
  WRITE(*,'(A4,A19,2ES10.2E2)') '', 'Compute_Pair = ',  &
    Timer_Compute_Pair, Timer_Compute_Pair / Timer_Total
  WRITE(*,'(A4,A19,2ES10.2E2)') '', 'Compute_Pair (P) = ',  &
    Timer_Compute_Pair_Point, Timer_Compute_Pair_Point / Timer_Total
  WRITE(*,*)

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
