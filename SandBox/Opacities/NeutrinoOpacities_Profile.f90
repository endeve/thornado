PROGRAM NeutrinoOpacities_Profile

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
    ComputeNeutrinoOpacities_EC_Points, &
    ComputeNeutrinoOpacities_ES_Points, &
    ComputeNeutrinoOpacities_NES_Points, &
    ComputeNeutrinoOpacities_Pair_Points
  USE DeviceModule, ONLY: &
    InitializeDevice, &
    FinalizeDevice

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(64) :: ProfileName = 'Chimera_S25_100ms_Boltztran_h.profile'
  INTEGER       :: nE, nX(3), nNodes, nSpecies, nPointsE, nPointsX
  REAL(DP)      :: eL, eR, ZoomE
  REAL(DP)      :: xL(3), xR(3), ZoomX(3)
  REAL(DP), PARAMETER :: &
    Unit_D     = Gram / Centimeter**3, &
    Unit_T     = Kelvin, &
    Unit_Y     = 1.0_DP, &
    Unit_E     = MeV, &
    Unit_Chi   = 1.0_DP / Centimeter, &
    Unit_Sigma = 1.0_DP / Centimeter, &
    Unit_Phi   = Unit_Chi / Unit_E**3

  INTEGER :: &
    mpierr, iE, iX, iS, iNodeE, iN_E
  REAL(DP), DIMENSION(:), ALLOCATABLE :: E, R, D, T, Y
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    Chi,     & ! --- Absorption Opacity
    Sigma,   & ! --- Scattering Opacity (Isoenergetic + 1st)
    Sigma_0, & ! --- Scattering Opacity (Isoenergetic)
    Sigma_1, & ! --- Scattering Opacity (1st)
    Chi_NES, & ! --- Integrated NES Opacity
    Chi_Pair   ! --- Integrated Pair Opacity
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: &
    Phi_0_NES_In,  &
    Phi_0_NES_Out,  &
    Phi_0_Pair_In,  &
    Phi_0_Pair_Out

  nNodes   = 2
  nSpecies = 2

  CALL ReadTextProfile( TRIM(ProfileName), R, D, T, Y )

  nPointsX = SIZE(R)
  nX = [ nPointsX/nNodes, 1, 1 ]

  D = D * Unit_D
  T = T * Unit_T
  Y = Y * Unit_Y

  nE       = 10
  eL       = 0.0d0 * MeV
  eR       = 3.0d2 * MeV
  ZoomE    = 1.526_DP
  nPointsE = nE * nNodes

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'NeutrinoOpacities_Profile', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 01, 1, 1 ], &
           bcX_Option &
             = [ 32, 0, 0 ], &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           ZoomX_Option &
             = ZoomX, &
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
             = 'SPHERICAL', &
           ActivateUnits_Option &
             = .TRUE., &
           nSpecies_Option &
             = nSpecies, &
           BasicInitialization_Option &
             = .TRUE. )

  WRITE(*,*)
  WRITE(*,'(A4,A)') '', 'NeutrinoOpacities_Profile'
  WRITE(*,*)

  ! --- Initialize Energy Grid ---

  ALLOCATE( E(nPointsE) )

  DO iN_E = 1, nPointsE
    iE      = MOD( (iN_E-1) / nNodes, nE     ) + 1
    iNodeE  = MOD( (iN_E-1)         , nNodes ) + 1
    E(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )
    WRITE(*,'(A6,A2,I3.3,A10,ES8.2E2)') '','E(',iN_E,') [MeV] = ', E(iN_E) / Unit_E
  END DO

  ! --- Initialize Equation of State ---

  CALL InitializeEquationOfState_TABLE &
         ( EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5', &
           Verbose_Option = .TRUE. )

  ! --- Initialize Opacities ---

  ALLOCATE( Chi(nPointsE,nPointsX,nSpecies) )
  ALLOCATE( Sigma(nPointsE,nPointsX,nSpecies) )
  ALLOCATE( Sigma_0(nPointsE,nPointsX,nSpecies) )
  ALLOCATE( Sigma_1(nPointsE,nPointsX,nSpecies) )
  ALLOCATE( Chi_NES(nPointsE,nPointsX,nSpecies) )
  ALLOCATE( Chi_Pair(nPointsE,nPointsX,nSpecies) )

  ALLOCATE( Phi_0_NES_In(nPointsE,nPointsE,nPointsX,nSpecies) )
  ALLOCATE( Phi_0_NES_Out(nPointsE,nPointsE,nPointsX,nSpecies) )
  ALLOCATE( Phi_0_Pair_In(nPointsE,nPointsE,nPointsX,nSpecies) )
  ALLOCATE( Phi_0_Pair_Out(nPointsE,nPointsE,nPointsX,nSpecies) )

  CALL InitializeOpacities_TABLE &
         ( OpacityTableName_EmAb_Option &
             = 'wl-Op-LS220-20-40-100-Lower-T-E40-B85-EmAb.h5', &
           OpacityTableName_Iso_Option  &
             = 'wl-Op-LS220-20-40-100-Lower-T-E40-B85-Iso.h5',  &
           OpacityTableName_NES_Option &
             = 'wl-Op-LS220-20-40-100-Lower-T-E40-B85-NES.h5',  &
           OpacityTableName_Pair_Option &
             = 'wl-Op-LS220-20-40-100-Lower-T-E40-B85-Pair.h5', &
           Verbose_Option = .TRUE. )

  ! --- Compute Electron Capture Opacities ---

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_EC_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, Chi(:,:,iS) )

  END DO

  ! --- Compute Elastic Scattering Opacities ---

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 1, Sigma_0(:,:,iS) )

    CALL ComputeNeutrinoOpacities_ES_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 2, Sigma_1(:,:,iS) )

    Sigma(:,:,iS) = Sigma_0(:,:,iS) - Sigma_1(:,:,iS) / 3.0d0

  END DO

  ! --- Compute NES Opacities ---

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_NES_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 1, &
             Phi_0_NES_In(:,:,:,iS), Phi_0_NES_Out(:,:,:,iS) )

  END DO

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

  DO iS = 1, nSpecies

    CALL ComputeNeutrinoOpacities_Pair_Points &
           ( 1, nPointsE, 1, nPointsX, E, D, T, Y, iS, 1, &
             Phi_0_Pair_In(:,:,:,iS), Phi_0_Pair_Out(:,:,:,iS) )

  END DO

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

  CALL WriteVector &
         ( nPointsX, R, 'R.dat' )

  CALL WriteFourOpacities( Chi, Sigma, Chi_NES, Chi_Pair, 'mfp.dat' )

  CALL Write4DOpacities( Phi_0_NES_Out, 'Phi_0_NES_Out.dat' )

  CALL Write4DOpacities( Phi_0_Pair_Out, 'Phi_0_Pair_Out.dat' )

  CALL FinalizeEquationOfState_TABLE

  CALL FinalizeOpacities_TABLE

  DEALLOCATE( E, R, D, T, Y )
  DEALLOCATE( &
    Chi, Sigma, Sigma_0, Sigma_1, Chi_NES, Chi_Pair, &
    Phi_0_NES_In,  &
    Phi_0_NES_Out, &
    Phi_0_Pair_In, &
    Phi_0_Pair_Out )

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


  SUBROUTINE ReadTextProfile &
    ( FileName, R, D, T, Y )

    CHARACTER(*), INTENT(in) :: FileName
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: &
      R, D, T, Y

    CHARACTER(2) :: a
    INTEGER      :: i, datasize, ipos, nvar_stored
    REAL(dp), DIMENSION(:), ALLOCATABLE :: database

    CHARACTER(LEN=07)  :: Format1 = '(A2,I5)'
    CHARACTER(LEN=06)  :: Format2 = '(4A15)'
    CHARACTER(LEN=09)  :: Format3 = '(4ES15.6)'
    CHARACTER(LEN=80)  :: Current_Line

    PRINT*
    WRITE(*,*) 'Using Text Profile :', FileName

    OPEN( 1, FILE = FileName, FORM = "formatted", &
            ACTION = 'read')
    READ( 1, Format1 ) a, datasize
    READ( 1, Format2 )

    ALLOCATE( database ( datasize * 4) )
    ALLOCATE( R   ( datasize ) )
    ALLOCATE( D ( datasize ) )
    ALLOCATE( T   ( datasize ) )
    ALLOCATE( Y  ( datasize ) )

    READ( 1, Format3 ) database
    CLOSE( 1, STATUS = 'keep')

    DO i = 1, datasize
      R(i) = database(i*4-3)
      D(i) = database(i*4-2)
      T(i) = database(i*4-1)
      Y(i) = database(i*4)
    END DO

  END SUBROUTINE ReadTextProfile


  SUBROUTINE WriteFourOpacities( Chi, Sigma, Chi_NES, Chi_Pair, FileName )

    CHARACTER(len=*), INTENT(in) :: FileName
    REAL(DP),         INTENT(in) :: Chi(:,:,:), Sigma(:,:,:), &
                                    Chi_NES(:,:,:), Chi_Pair(:,:,:)

    INTEGER :: FUNIT, iE, iX, iS, nE, nX, nS, nDim(3)

    nDim = SHAPE( Chi )
    nE   = nDim(1)
    nX   = nDim(2)
    nS   = nDim(3)

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    DO iS = 1, nS
    DO iX = 1, nX
    DO iE = 1, nE

    WRITE( FUNIT, '(3I5,4ES15.4)' ) iS, iX, iE, &
                                    Chi(iE,iX,iS) / Unit_Chi, &
                                    Sigma(iE,iX,iS) / Unit_Sigma, &
                                    Chi_NES(iE,iX,iS) / Unit_Chi, &
                                    Chi_Pair(iE,iX,iS) / Unit_Chi

    END DO
    END DO
    END DO

    CLOSE( FUNIT )

  END SUBROUTINE WriteFourOpacities


  SUBROUTINE Write4DOpacities( Opacity, FileName )

    CHARACTER(len=*), INTENT(in) :: FileName
    REAL(DP),         INTENT(in) :: Opacity(:,:,:,:)

    INTEGER :: FUNIT, iE1, iE2, iX, iS, nE, nX, nS, nDim(4)

    nDim = SHAPE( Opacity )
    nE   = nDim(1)
    nX   = nDim(3)
    nS   = nDim(4)

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    DO iS = 1, nS
    DO iX = 1, nX
    DO iE2 = 1, nE
    DO iE1 = 1, nE

    WRITE( FUNIT, '(4I5,ES15.4)' ) iS, iX, iE2, iE1, &
                                    Opacity(iE1,iE2,iX,iS) / Unit_Phi

    END DO
    END DO
    END DO
    END DO

    CLOSE( FUNIT )

  END SUBROUTINE Write4DOpacities

END PROGRAM NeutrinoOpacities_Profile
