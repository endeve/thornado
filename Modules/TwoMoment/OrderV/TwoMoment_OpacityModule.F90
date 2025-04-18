MODULE TwoMoment_OpacityModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Three, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFZ, nDOFE
  USE ReferenceElementModuleZ, ONLY: &
    NodeNumberTableZ
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE RadiationFieldsModule, ONLY: &
    nSpecies

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


  SUBROUTINE SetOpacities &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D0, Chi, Sigma, Verbose_Option )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D0, Chi, Sigma
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A5,A)') '', 'Setting Opacities:'
      WRITE(*,*)
      WRITE(*,'(A7,A8,ES10.4E2)') '',    'D0 = ', D0
      WRITE(*,'(A7,A8,ES10.4E2)') '',   'Chi = ', Chi
      WRITE(*,'(A7,A8,ES10.4E2)') '', 'Sigma = ', Sigma
    END IF

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'ExpandingAtmosphere' )

        CALL SetOpacities_ExpandingAtmosphere &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

      CASE( 'HomogeneousSphere1D' )

        CALL SetOpacities_HomogeneousSphere1D &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D0, Chi )

      CASE( 'HomogeneousSphere2D' )

        CALL SetOpacities_HomogeneousSphere2D &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D0, Chi )

      CASE( 'ShadowCasting2D_Cartesian' )

        CALL SetOpacities_ShadowCasting2D_Cartesian &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

      CASE( 'ShadowCasting2D_Cylindrical' )

        CALL SetOpacities_ShadowCasting2D_Cylindrical &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

      CASE DEFAULT

        uOP(:,:,:,:,:,iOP_D0   ,:) = D0
        uOP(:,:,:,:,:,iOP_Chi  ,:) = Chi
        uOP(:,:,:,:,:,iOP_Sigma,:) = Sigma

    END SELECT

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: uOP )
#elif defined( THORNADO_OACC   )
    !$ACC ENTER DATA &
    !$ACC COPYIN( uOP )
#endif

  END SUBROUTINE SetOpacities


  SUBROUTINE SetOpacities_ExpandingAtmosphere &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)

    REAL(DP), PARAMETER :: Rmax  = 11.0_DP
    REAL(DP), PARAMETER :: E_0   = 3.0_DP
    REAL(DP), PARAMETER :: Delta = 0.2_DP
    REAL(DP), PARAMETER :: Alpha = 10.9989_DP

    INTEGER  :: iNodeZ, iNodeZ1, iNodeZ2
    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: E, R

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeZ1 = NodeNumberTableZ(1,iNodeZ)
        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)

        E = NodeCoordinate( MeshE   , iZ1, iNodeZ1 )
        R = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )

        IF( R < Rmax )THEN

          IF( E > E_0 )THEN

            uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Chi,iS) = 1.0d1 * Alpha / R**2

          ELSE

            uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Chi,iS) &
              = Alpha * ( 1.0d1 * EXP( - (E-E_0)**2/Delta**2 ) &
                           + (One-EXP( - (E-E_0)**2/Delta**2 )) ) / R**2 

          END IF

        ELSE

          uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Chi,iS) = Zero

        END IF

        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_D0   ,iS) = One / ( EXP(E) - One )
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS) = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE SetOpacities_ExpandingAtmosphere


  SUBROUTINE SetOpacities_HomogeneousSphere1D &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D0, Chi )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D0, Chi

    REAL(DP), PARAMETER :: R_0 = 1.0d-00 ! --- Radius of Sphere
    REAL(DP), PARAMETER :: L_R = 1.0d-08 ! --- Smoothing Lenght

    INTEGER  :: iNodeZ, iNodeZ2, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: Radius

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)

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
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, D0, Chi )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      D0, Chi

    REAL(DP), PARAMETER :: R_0 = 1.0d-00 ! --- Radius of Sphere
    REAL(DP), PARAMETER :: L_R = 1.0d-01 ! --- Smoothing Lenght

    INTEGER  :: iNodeZ, iNodeE, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: E, X1_C, X2_C, X3_C, R_C, Spectrum

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      X1_C = MeshX(1) % Center(iZ2)
      X2_C = MeshX(2) % Center(iZ3)
      X3_C = MeshX(3) % Center(iZ4)

      R_C = SQRT( X1_C**2 + X2_C**2 + X3_C**2 )

      DO iNodeZ = 1, nDOFZ

        IF( iZ_E0(1) .EQ. iZ_B0(1) )THEN

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


  SUBROUTINE SetOpacities_ShadowCasting2D_Cartesian &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)

    REAL(DP), PARAMETER :: R_0_A = 2.0d+00 ! --- Radius of Absorbing Region
    REAL(DP), PARAMETER :: R_0_S = 1.5d+00 ! --- Radius of Radiating Region

    INTEGER  :: iNodeZ, iNodeZ2, iNodeZ3, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: X1, X2, Distance_A, Distance_S, D0_loc, Chi_loc

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)
        iNodeZ3 = NodeNumberTableZ(3,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iZ2, iNodeZ2 )
        X2 = NodeCoordinate( MeshX(2), iZ3, iNodeZ3 )

        ! Distance to Center of Absorbing Region (x,y) = ( 11, 0 )
        Distance_A = SQRT( (X1 - 1.1d1 )**2 + X2**2 )
        Chi_loc = Zero
        D0_loc  = Zero

        IF( Distance_A <= R_0_A )THEN ! Inside Absorbing Region
           Chi_loc = 10.d0
        ELSE
          ! Distance to Center of Radiation Region (x,y) = ( 3, 0 )
           Distance_S = SQRT( (X1 - 3.0d0 )**2 + X2**2 )
           IF( Distance_S <= R_0_S )THEN ! Inside Radiating Region
             Chi_loc = 10.d0 * EXP( - 2.0d0 * ( Distance_S / R_0_S )**2 )
             D0_loc  = 1.0d-1
           END IF
        END IF

        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_D0   ,iS) &
          = D0_loc
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Chi  ,iS) &
          = Chi_loc
        uOP(iNodeZ,iZ1,iZ2,iZ3,iZ4,iOP_Sigma,iS) &
          = Zero

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE SetOpacities_ShadowCasting2D_Cartesian


  SUBROUTINE SetOpacities_ShadowCasting2D_Cylindrical &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1 )

    ! --- {Z1,Z2,Z3,Z4} = {E,X1,X2,X3} ---

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)

    REAL(DP), PARAMETER :: R_A = 2.0d+00 ! --- Radius of Absorbing Region
    REAL(DP), PARAMETER :: R_S = 1.5d+00 ! --- Radius of Radiating Region

    INTEGER  :: iNodeZ, iNodeZ2, iNodeZ3, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: X1, X2, Radius, Distance_A, WindowFunction, D0, Chi

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B1(1), iZ_E1(1)

      DO iNodeZ = 1, nDOFZ

        iNodeZ2 = NodeNumberTableZ(2,iNodeZ)
        iNodeZ3 = NodeNumberTableZ(3,iNodeZ)

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

  END SUBROUTINE SetOpacities_ShadowCasting2D_Cylindrical


  SUBROUTINE CreateOpacities( nX, swX, nE, swE, Verbose_Option )

    INTEGER, INTENT(in) :: nX(3), swX(3), nE, swE
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
            1-swE:nE+swE, &
            1-swX(1):nX(1)+swX(1), &
            1-swX(2):nX(2)+swX(2), &
            1-swX(3):nX(3)+swX(3), &
            1:nOP,1:nSpecies) )

    uOP = Zero

  END SUBROUTINE CreateOpacities


  SUBROUTINE DestroyOpacities

#if   defined( THORNADO_OMP_OL )
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: uOP )
#elif defined( THORNADO_OACC   )
    !$ACC EXIT DATA &
    !$ACC DELETE( uOP )
#endif

    DEALLOCATE( uOP )

  END SUBROUTINE DestroyOpacities


END MODULE TwoMoment_OpacityModule
