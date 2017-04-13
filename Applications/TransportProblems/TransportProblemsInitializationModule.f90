MODULE TransportProblemsInitializationModule

  USE KindModule, ONLY: &
    DP, Pi, FourPi
  USE UnitsModule, ONLY: &
    Centimeter, &
    Gram, &
    Kelvin, &
    MeV, &
    Kilometer, &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Gm, iAF_Cs
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved
  USE RadiationFieldsModule, ONLY: &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule, ONLY: &
    ApplyEquationOfState, &
    ComputeThermodynamicStates_Primitive
  USE OpacityModule, ONLY: &
    ComputeScatteringOpacity_ES
  USE MomentEquationsUtilitiesModule, ONLY: &
    Conserved

  IMPLICIT NONE
  PRIVATE

  TYPE ProfileType
    INTEGER                             :: nPoints
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Radius
    REAL(DP), DIMENSION(:), ALLOCATABLE :: MassDensity
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Temperature
    REAL(DP), DIMENSION(:), ALLOCATABLE :: ElectronFraction
  END TYPE ProfileType

  PUBLIC :: InitializeHomogeneousSphere1D
  PUBLIC :: InitializeGaussianSphericalDiffusion1D
  PUBLIC :: InitializeDeleptonizationProblem1D

CONTAINS


  SUBROUTINE InitializeHomogeneousSphere1D( CentralConditions_Option )

    CHARACTER(2), INTENT(in), OPTIONAL :: CentralConditions_Option

    CHARACTER(2) :: CentralConditions
    INTEGER      :: iX1, iX2, iX3, iE
    INTEGER      :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER      :: iNodeX, iNode
    REAL(DP)     :: X1, E
    REAL(DP)     :: Radius

    CentralConditions = '02'
    IF( PRESENT( CentralConditions_Option ) ) &
      CentralConditions = CentralConditions_Option

    Radius = 1.0d2 * Kilometer

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    ! --- Initialize Fluid Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                IF( X1 <= Radius )THEN

                  SELECT CASE ( CentralConditions )

                    CASE ( '01' )

                      uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                        = 1.0d14 * Gram / Centimeter**3
                      uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                        = 21.0_DP * MeV
                      uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                        = 0.25_DP

                    CASE ( '02' )

                      uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                        = 1.0d13 * Gram / Centimeter**3
                      uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                        = 16.0_DP * MeV
                      uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                        = 0.14_DP

                    CASE ( '03' )

                      uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                        = 1.0d12 * Gram / Centimeter**3
                      uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                        = 8.0_DP * MeV
                      uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                        = 0.12_DP

                    CASE ( '04' )

                      uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                        = 1.0d11 * Gram / Centimeter**3
                      uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                        = 8.0_DP * MeV
                      uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                        = 0.15_DP

                    CASE ( '05' )

                      uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                        = 1.0d10 * Gram / Centimeter**3
                      uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                        = 3.0_DP * MeV
                      uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                        = 0.26_DP

                    CASE DEFAULT

                      uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                        = 1.0d13 * Gram / Centimeter**3
                      uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                        = 16.0_DP * MeV
                      uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                        = 0.14_DP

                  END SELECT


                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                    = 1.0d8 * Gram / Centimeter**3

                  uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                    = 0.2_DP * MeV

                  uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                    = 0.4643_DP

                END IF

              END DO
            END DO
          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Gm) )

        END DO
      END DO
    END DO

    ! --- Initialize Radiation Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                  DO iNodeE = 1, nNodesE

                    E = NodeCoordinate( MeshE, iE, iNodeE )

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_D,1) &
                      = 1.0d-6

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,1) &
                      = 0.0_DP

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,1) &
                      = 0.0_DP

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,1) &
                      = 0.0_DP

                    uCR(iNode,iE,iX1,iX2,iX3,1:nCR,1) &
                      = Conserved( uPR(iNode,iE,iX1,iX2,iX3,1:nPR,1) )

                  END DO

                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeHomogeneousSphere1D


  SUBROUTINE InitializeGaussianSphericalDiffusion1D &
               ( t_0, BackgroundConditions_Option )

    REAL(DP),     INTENT(in)           :: t_0
    CHARACTER(2), INTENT(in), OPTIONAL :: BackgroundConditions_Option

    CHARACTER(2)             :: BackgroundConditions
    INTEGER                  :: iX1, iX2, iX3, iE
    INTEGER                  :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER                  :: iNodeX, iNode
    REAL(DP)                 :: X1, E_0, E
    REAL(DP), DIMENSION(1,1) :: Kappa

    BackgroundConditions = '02'
    IF( PRESENT( BackgroundConditions_Option ) ) &
      BackgroundConditions = BackgroundConditions_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    ! --- Initialize Fluid Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                SELECT CASE ( BackgroundConditions )

                  CASE ( '01' )

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                      = 1.0d14 * Gram / Centimeter**3
                    uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                      = 21.0_DP * MeV
                    uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                      = 0.25_DP

                  CASE ( '02' )

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                      = 1.0d13 * Gram / Centimeter**3
                    uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                      = 16.0_DP * MeV
                    uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                      = 0.14_DP

                  CASE ( '03' )

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                      = 1.0d12 * Gram / Centimeter**3
                    uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                      = 8.0_DP * MeV
                    uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                      = 0.12_DP

                  CASE ( '04' )

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                      = 1.0d11 * Gram / Centimeter**3
                    uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                      = 8.0_DP * MeV
                    uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                      = 0.15_DP

                  CASE ( '05' )

                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                      = 1.0d10 * Gram / Centimeter**3
                    uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                      = 3.0_DP * MeV
                    uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                      = 0.26_DP

                  END SELECT

              END DO
            END DO
          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Gm) )

        END DO
      END DO
    END DO

    ! --- Initialize Radiation Fields ---

    E_0 = NodeCoordinate( MeshE, 1, 1 )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                  DO iNodeE = 1, nNodesE

                    E = NodeCoordinate( MeshE, iE, iNodeE )

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    CALL ComputeScatteringOpacity_ES &
                           ( [ E ], &
                             [ uPF(iNodeX,iX1,iX2,iX3,iPF_D ) ], &
                             [ uAF(iNodeX,iX1,iX2,iX3,iAF_T ) ], &
                             [ uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) ], &
                             Kappa )

                    WRITE(*,'(A6,A10,ES12.4E2,A16,ES12.4E2)') &
                      '', 'E [MeV] = ', E / MeV, &
                      ' Kappa [1/cm] = ', Kappa / ( 1.0_DP / Centimeter )

                    ! --- Number Density ---

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_D,1) &
                      = EXP( - 3.0_DP * Kappa(1,1) * ( X1 * E_0 / E )**2 &
                               / ( 4.0_DP * t_0 ) )

                    ! --- Number Flux Density (1) ---

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,1) &
                      = X1 * uPR(iNode,iE,iX1,iX2,iX3,iPR_D,1) &
                          / ( 2.0_DP * t_0 )

                    ! --- Number Flux Density (2) ---

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,1) &
                      = 0.0_DP

                    ! --- Number Flux Density (3) ---

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,1) &
                      = 0.0_DP

                    ! --- Compute Conserved Radiation Fields ---

                    uCR(iNode,iE,iX1,iX2,iX3,1:nCR,1) &
                      = Conserved( uPR(iNode,iE,iX1,iX2,iX3,1:nPR,1) )

                  END DO

                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeGaussianSphericalDiffusion1D


  SUBROUTINE InitializeDeleptonizationProblem1D &
               ( Temperature_Option, ElectronFraction_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Temperature_Option
    REAL(DP), INTENT(in), OPTIONAL :: ElectronFraction_Option

    INTEGER             :: iX1, iX2, iX3, iR, iE
    INTEGER             :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER             :: iNodeX, iNode
    REAL(DP)            :: MassDensity
    REAL(DP)            :: Temperature
    REAL(DP)            :: ElectronFraction
    REAL(DP)            :: X1, E, Mnu, kT
    REAL(DP)            :: R_2, X1_0, Theta
    REAL(DP), PARAMETER :: R_0   = 1.0d00 * Kilometer
    REAL(DP), PARAMETER :: D_1   = 4.5d14 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: R_1   = 1.0d01 * Kilometer
    REAL(DP), PARAMETER :: D_2   = 1.0d13 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: Alpha = - 4.0_DP

    Temperature = 1.0d11 * Kelvin
    IF( PRESENT( Temperature_Option ) ) &
      Temperature = Temperature_Option

    ElectronFraction = 0.3_DP
    IF( PRESENT( ElectronFraction_Option ) ) &
      ElectronFraction = ElectronFraction_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    R_2  = R_1 * SQRT( LOG( D_1 / D_2 ) )
    X1_0 = R_2

    ! --- Initialize Fluid Fields with Profile ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                Theta = 0.5_DP * ( 1.0_DP + TANH( ( X1 - X1_0 ) / R_0 ) )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                MassDensity &
                  = (1.0_DP - Theta) * D_1 * EXP( - ( X1 / R_1 )**2 ) &
                      + Theta * D_2 * ( X1 / R_2 )**Alpha

                uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                  = MassDensity

                uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                  = (1.0_DP - Theta) * Temperature &
                      + Theta * Temperature * ( X1_0 / X1 )**1

                uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                  = ElectronFraction

              END DO
            END DO
          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

          CALL ApplyEquationOfState &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
                   uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Gm) )

          uAF(:,iX1,iX2,iX3,iAF_Cs) &
            = SQRT( uAF(:,iX1,iX2,iX3,iAF_Gm) &
                    * uAF(:,iX1,iX2,iX3,iAF_P) &
                    / uPF(:,iX1,iX2,iX3,iPF_D) )

        END DO
      END DO
    END DO

    ! --- Initialize Radiation Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                  kT  = BoltzmannConstant &
                        * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

                  Mnu = uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                        + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                        - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

                  DO iNodeE = 1, nNodesE

                    E = NodeCoordinate( MeshE, iE, iNodeE )

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_D,1) &
                      = MAX( 1.0d-8, FourPi / ( EXP( (E-Mnu)/kT ) + 1.0_DP ) )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,1) &
                      = 0.0_DP

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,1) &
                      = 0.0_DP

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,1) &
                      = 0.0_DP

                    uCR(iNode,iE,iX1,iX2,iX3,1:nCR,1) &
                      = Conserved( uPR(iNode,iE,iX1,iX2,iX3,1:nPR,1) )

                  END DO

                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeDeleptonizationProblem1D


  SUBROUTINE CreateProfile( P, N )

    TYPE(ProfileType)   :: P
    INTEGER, INTENT(in) :: N

    P % nPoints = N
    ALLOCATE( P % Radius          ( N ) )
    ALLOCATE( P % MassDensity     ( N ) )
    ALLOCATE( P % Temperature     ( N ) )
    ALLOCATE( P % ElectronFraction( N ) )

  END SUBROUTINE CreateProfile


  SUBROUTINE ReadProfile( Profile, ProfileName, Reverse_Option )

    TYPE(ProfileType)                      :: Profile
    CHARACTER(LEN=*), INTENT(in)           :: ProfileName
    LOGICAL,          INTENT(in), OPTIONAL :: Reverse_Option

    LOGICAL                     :: Reverse
    CHARACTER(80), DIMENSION(5) :: Buffer
    INTEGER                     :: FUNIT, STAT, i, Increment
    REAL(DP), DIMENSION(5)      :: BufferDP

    WRITE(*,*)
    WRITE(*,'(A4,A17,A)') &
      '', 'Reading Profile: ', TRIM( ProfileName )
    WRITE(*,*)

    Reverse = .FALSE.
    IF( PRESENT( Reverse_Option ) ) &
      Reverse = Reverse_Option

    Increment = + 1
    IF( Reverse ) &
      Increment = - 1
    
    OPEN( NEWUNIT = FUNIT, FILE = TRIM( ProfileName ) )

    READ( FUNIT, * ) Buffer

    i = 0
    IF( Reverse ) i = Profile % nPoints + 1

    DO

      READ( FUNIT, * , IOSTAT = STAT ) BufferDP

      IF( STAT /= 0 ) EXIT

      i = i + Increment

      Profile % Radius(i) &
        = BufferDP(1) * Centimeter
      Profile % MassDensity(i) &
        = BufferDP(2) * Gram / Centimeter**3
      Profile % Temperature(i) &
        = BufferDP(3) * Kelvin
      Profile % ElectronFraction(i) &
        = BufferDP(4)

    END DO

    CLOSE( FUNIT )

  END SUBROUTINE ReadProfile


  INTEGER FUNCTION nLines( ProfileName )

    CHARACTER(LEN=*), INTENT(in) :: ProfileName

    CHARACTER(80) :: Buffer
    INTEGER       :: FUNIT, STAT

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( ProfileName ) )

    READ( FUNIT, * ) Buffer

    nLines = 0
    DO

      READ( FUNIT, * , IOSTAT = STAT ) Buffer

      IF( STAT /= 0 ) EXIT

      nLines = nLines + 1

    END DO

    CLOSE( FUNIT )

  END FUNCTION nLines


END MODULE TransportProblemsInitializationModule
