MODULE TwoMoment_OpacityModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Three
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFZ, nDOFE
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    nSpecies
  USE UnitsModule,            ONLY: &
    UnitsDisplay, &
    Centimeter

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC :: iOP_D0    = 1
  INTEGER, PUBLIC :: iOP_Chi   = 2
  INTEGER, PUBLIC :: iOP_Sigma = 3
  INTEGER, PUBLIC :: nOP       = 3

  REAL(DP), PUBLIC, ALLOCATABLE :: uOP(:,:,:,:,:,:,:)

  PUBLIC :: SetOpacities
  PUBLIC :: CreateOpacities
  PUBLIC :: DestroyOpacities

CONTAINS


  SUBROUTINE SetOpacities( iZ_B1, iZ_E1, iOS_CPP, D0, Chi, Sigma, kT, E0, mu0, R0, Verbose_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B1(4), iZ_E1(4), iOS_CPP(3)
    REAL(DP), INTENT(in) :: &
      D0, Chi, Sigma, kT, E0, mu0, R0
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option


    LOGICAL :: Verbose

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Setting Opacities:'
      WRITE(*,*)
      WRITE(*,'(A7,A8,ES10.4E2)') '',    'D0 = ', D0
      WRITE(*,'(A7,A8,ES10.4E2)') '',   'Chi = ', Chi
      WRITE(*,'(A7,A8,ES10.4E2)') '', 'Sigma = ', Sigma
    END IF

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'HomogeneousSphere1D' )

        CALL SetOpacities_HomogeneousSphere1D &
               ( iZ_B1, iZ_E1, iOS_CPP, D0, Chi )

      CASE( 'HomogeneousSphere2D' )

        CALL SetOpacities_HomogeneousSphere2D &
               ( iZ_B1, iZ_E1, iOS_CPP, D0, Chi )

      CASE( 'HomogeneousSphereGR' )

        CALL SetOpacities_HomogeneousSphereGR &
               ( iZ_B1, iZ_E1, iOS_CPP, D0, Chi, kT, E0, mu0, R0 )
      CASE( 'ShadowCasting' )

        CALL SetOpacities_ShadowCasting &
               ( iZ_B1, iZ_E1, iOS_CPP )

      CASE DEFAULT

        uOP(:,:,:,:,:,iOP_D0   ,:) = D0
        uOP(:,:,:,:,:,iOP_Chi  ,:) = Chi
        uOP(:,:,:,:,:,iOP_Sigma,:) = Sigma


    END SELECT

  END SUBROUTINE SetOpacities


  SUBROUTINE SetOpacities_HomogeneousSphere1D &
    ( iZ_B1, iZ_E1, iOS_CPP, D0, Chi )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B1(4), iZ_E1(4), iOS_CPP(3)
    REAL(DP), INTENT(in) :: &
      D0, Chi

    REAL(DP), PARAMETER :: R_0 = 1.0d-00 ! --- Radius of Sphere
    REAL(DP), PARAMETER :: L_R = 1.0d-08 ! --- Smoothing Lenght

    INTEGER  :: iNodeZ, iNodeZ2, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: Radius

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B1(4)-iOS_CPP(3), iZ_E1(4)-iOS_CPP(3)
    DO iZ3 = iZ_B1(3)-iOS_CPP(2), iZ_E1(3)-iOS_CPP(2)
    DO iZ2 = iZ_B1(2)-iOS_CPP(1), iZ_E1(2)-iOS_CPP(1)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        iNodeZ2 = NodeNumberTable(2,iNodeZ)

        Radius = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )

        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_D0   ,iS) &
          = D0
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Chi  ,iS) &
          = Chi * Half * ( One - TANH( ( Radius - R_0 ) / L_R ) )
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS) &
          = Zero
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE SetOpacities_HomogeneousSphere1D


  SUBROUTINE SetOpacities_HomogeneousSphere2D &
    ( iZ_B1, iZ_E1, iOS_CPP, D0, Chi )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B1(4), iZ_E1(4), iOS_CPP(3)
    REAL(DP), INTENT(in) :: &
      D0, Chi

    REAL(DP), PARAMETER :: R_0 = 1.0d-00 ! --- Radius of Sphere
    REAL(DP), PARAMETER :: L_R = 1.0d-01 ! --- Smoothing Lenght

    INTEGER  :: iNodeZ, iNodeE, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: E, X1_C, X2_C, X3_C, R_C, Spectrum

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B1(4)-iOS_CPP(3), iZ_E1(4)-iOS_CPP(3)
    DO iZ3 = iZ_B1(3)-iOS_CPP(2), iZ_E1(3)-iOS_CPP(2)
    DO iZ2 = iZ_B1(2)-iOS_CPP(1), iZ_E1(2)-iOS_CPP(1)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      X1_C = MeshX(1) % Center(iZ2)
      X2_C = MeshX(2) % Center(iZ3)
      X3_C = MeshX(3) % Center(iZ4)

      R_C = SQRT( X1_C**2 + X2_C**2 + X3_C**2 )

      DO iNodeZ = 1, nDOFZ

        IF( iZ_E1(1) .EQ. iZ_B1(1) )THEN

          Spectrum = One

        ELSE

          iNodeE = MOD( (iNodeZ-1), nDOFE ) + 1

          E = NodeCoordinate( MeshE, iZ1, iNodeE )

          Spectrum = One / ( EXP( E / Three - Three ) + One )

        END IF

        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_D0   ,iS) &
          = D0 * Spectrum
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Chi  ,iS) &
          = Chi * Half * ( One - TANH( ( R_C - R_0 ) / L_R ) )
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS) &
          = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE SetOpacities_HomogeneousSphere2D

  SUBROUTINE SetOpacities_HomogeneousSphereGR &
    ( iZ_B1, iZ_E1, iOS_CPP, D0, Chi, kT, E0, mu0, R0 )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B1(4), iZ_E1(4), iOS_CPP(3)
    REAL(DP), INTENT(in) :: &
      D0, Chi, kT, E0, mu0, R0

    REAL(DP), PARAMETER :: L_R = 3000.0_DP!1.0d-04  --- Smoothing Lenght

    INTEGER  :: iNodeZ, iNodeZ2, iZ1, iZ2, iZ3, iZ4, iS, iNodeE, iNodeX
    REAL(DP) :: Radius, E, Spectrum


    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B1(4)-iOS_CPP(3), iZ_E1(4)-iOS_CPP(3)
    DO iZ3 = iZ_B1(3)-iOS_CPP(2), iZ_E1(3)-iOS_CPP(2)
    DO iZ2 = iZ_B1(2)-iOS_CPP(1), iZ_E1(2)-iOS_CPP(1)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        iNodeE = MOD( (iNodeZ-1), nDOFE ) + 1
 
        E = NodeCoordinate( MeshE, iZ1, iNodeE )
!print*, E / UnitsDisplay % EnergyUnit

        Spectrum = One / ( EXP( ( E - mu0 ) / kT  ) + One )

        iNodeZ2 = NodeNumberTable(2,iNodeZ)

        Radius = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )

        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_D0   ,iS) &
          = D0 * Spectrum
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Chi  ,iS) &
          = Chi * ( E / E0 )**2 * Half * ( One - TANH( ( Radius - R0 ) / L_R ) )
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS) &
          = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO
  END SUBROUTINE SetOpacities_HomogeneousSphereGR



  SUBROUTINE SetOpacities_ShadowCasting &
    ( iZ_B1, iZ_E1, iOS_CPP )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B1(4), iZ_E1(4), iOS_CPP(3)

    REAL(DP), PARAMETER :: R_A = 2.0d+00 ! --- Radius of Absorbing Region
    REAL(DP), PARAMETER :: R_S = 1.5d+00 ! --- Radius of Radiating Region

    INTEGER  :: iNodeZ, iNodeZ2, iNodeZ3, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: X1, X2, Radius, Distance_A, WindowFunction, D0, Chi

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B1(4)-iOS_CPP(3), iZ_E1(4)-iOS_CPP(3)
    DO iZ3 = iZ_B1(3)-iOS_CPP(2), iZ_E1(3)-iOS_CPP(2)
    DO iZ2 = iZ_B1(2)-iOS_CPP(1), iZ_E1(2)-iOS_CPP(1)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        iNodeZ2 = NodeNumberTable(2,iNodeZ)
        iNodeZ3 = NodeNumberTable(3,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )
        X2 = NodeCoordinate( MeshX(2), iZ3, iNodeZ3 )

        Radius = SQRT( X1**2 + X2**2 )

        ! --- Distance to Center of Absorbing Region (R,Z) = ( 8, 0 )

        Distance_A = SQRT( ( X1 - 8.0d0 )**2 + X2**2 )

        D0  = Zero
        Chi = Zero

        IF( Distance_A <= R_A )THEN ! --- Inside Absorbing Region

           Chi = 10.d0

        ELSE

          WindowFunction & ! --- Smoothing
            = MAX( 1.0d-16, Half * ( One - TANH((Radius-R_S)/0.05_DP) ) )

          Chi = 10.d0 * EXP( - 2.0d0 * ( Radius / R_S )**2 ) * WindowFunction
          D0  = 1.0d-1 * WindowFunction

        END IF
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_D0   ,iS) = D0
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Chi  ,iS) = Chi
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS) = Zero
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO


  END SUBROUTINE SetOpacities_ShadowCasting

  SUBROUTINE CreateOpacities( iZ_B1, iZ_E1, iOS_CPP, Verbose_Option )

    INTEGER, INTENT(in) :: iZ_B1(4), iZ_E1(4), iOS_CPP(3)
    LOGICAL, INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A29,I2.2)') &
        '', 'Simple Opacities, nSpecies = ', nSpecies
    END IF


    ALLOCATE &
      ( uOP(1:nDOFZ, &
            iZ_B1(1):iZ_E1(1), &
            iZ_B1(2)-iOS_CPP(1):iZ_E1(2)-iOS_CPP(1), &
            iZ_B1(3)-iOS_CPP(2):iZ_E1(3)-iOS_CPP(2), &
            iZ_B1(4)-iOS_CPP(3):iZ_E1(4)-iOS_CPP(3), &
            1:nOP,1:nSpecies) )

  END SUBROUTINE CreateOpacities


  SUBROUTINE DestroyOpacities

    DEALLOCATE( uOP )

  END SUBROUTINE DestroyOpacities


END MODULE TwoMoment_OpacityModule
