MODULE TransportProblemsInitializationModule

  USE KindModule, ONLY: &
    DP, Pi
  USE UnitsModule, ONLY: &
    Centimeter, &
    Gram, &
    Kelvin, &
    Kilometer, &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    Locate, &
    NodeNumberX, &
    NodeNumber, &
    Interpolate1D_Linear
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_E, iPF_Ne, &
    uAF, iAF_T, iAF_E, iAF_Ye, iAF_Me, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule, ONLY: &
    ApplyEquationOfState, &
    ComputeThermodynamicStates_Primitive
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

  PUBLIC :: InitializeCoolingProblem1D

CONTAINS


  SUBROUTINE InitializeCoolingProblem1D( ProfileName )

    CHARACTER(LEN=*), INTENT(in) :: ProfileName

    INTEGER           :: iX1, iX2, iX3, iR, iE
    INTEGER           :: iNodeX1, iNodeX2, iNodeX3, iNodeE
    INTEGER           :: iNodeX, iNode
    REAL(DP)          :: X1, E
    TYPE(ProfileType) :: P

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)

    CALL CreateProfile( P, nLines( ProfileName ) )

    CALL ReadProfile( P, ProfileName, Reverse_Option = .TRUE. )

    ! --- Initialize Fluid Fields with Profile ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX3 = 1, nNodesX(3)
            DO iNodeX2 = 1, nNodesX(2)
              DO iNodeX1 = 1, nNodesX(1)

                X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

                iR = Locate( X1, P % Radius, P % nPoints )

                iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                  = Interpolate1D_Linear &
                      ( X1, P % Radius(iR), P % Radius(iR+1), &
                        P % MassDensity(iR), P % MassDensity(iR+1) )

                uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
                  = 1.0d11 * Kelvin
!!$                  = Interpolate1D_Linear &
!!$                      ( X1, P % Radius(iR), P % Radius(iR+1), &
!!$                        P % Temperature(iR), P % Temperature(iR+1) )

                uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
                  = 0.3_DP
!!$                  = Interpolate1D_Linear &
!!$                      ( X1, P % Radius(iR), P % Radius(iR+1), &
!!$                        P % ElectronFraction(iR), P % ElectronFraction(iR+1) )

              END DO
            END DO
          END DO

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
                   uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E), &
                   uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

        END DO
      END DO
    END DO

    CALL ApplyEquationOfState

    ! --- Initialize Radiation Fields ---

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)

                  iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

                  ASSOCIATE &
                    ( kT  => BoltzmannConstant &
                               * uAF(iNodeX,iX1,iX2,iX3,iAF_T), &
                      Mnu => uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                               + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                               - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )

                  DO iNodeE = 1, nNodesE

                    E = NodeCoordinate( MeshE, iE, iNodeE )

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_D,1) &
                      = 4.0_DP * Pi &
                          / ( EXP( ( E - Mnu ) / kT ) + 1.0_DP )

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,1) &
                      = 0.0_DP

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,1) &
                      = 0.0_DP

                    uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,1) &
                      = 0.0_DP

                    uCR(iNode,iE,iX1,iX2,iX3,1:nCR,1) &
                      = Conserved( uPR(iNode,iE,iX1,iX2,iX3,1:nPR,1) )

                  END DO

                  END ASSOCIATE ! kT, etc.

                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeCoolingProblem1D


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
