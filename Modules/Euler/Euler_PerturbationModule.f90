MODULE Euler_PerturbationModule

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    Pi, &
    FourPi
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer
  USE ProgramHeaderModule, ONLY: &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  Use GeometryFieldsModule, ONLY: &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nPF, &
    iPF_D, &
    iPF_V1, &
    iPF_V2, &
    iPF_V3, &
    iPF_E, &
    iPF_Ne, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    nAF, &
    iAF_P, &
    iAF_T, &
    iAF_Ye, &
    iAF_S, &
    iAF_E, &
    iAF_Gm, &
    iAF_Cs
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic, &
    ComputeConserved_Euler_NonRelativistic, &
    ComputeFromConserved_Euler_NonRelativistic
  USE EquationOfStateModule, ONLY: &
    ComputeAuxiliary_Fluid

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeShellPerturbations, InitializeRandPerturbations, &
            ApplyShellPerturbations, ApplyRandPerturbations, &
            FinalizeRandPerturbations

  CHARACTER(LEN=32), PUBLIC :: PerturbType

  LOGICAL,  PUBLIC :: UseUnits
  INTEGER,  PUBLIC :: NumShells
  INTEGER,  PUBLIC :: BD_Order
  INTEGER,  PUBLIC :: CO_V_Nodes
  INTEGER,  PUBLIC :: MJ_D_Nodes
  INTEGER,  PUBLIC :: MJ_V_RadNodes
  INTEGER,  PUBLIC :: MJ_V_AngNodes
  REAL(DP), PUBLIC :: ShellWidth
  REAL(DP), PUBLIC :: BD_Amplitude
  REAL(DP), PUBLIC :: CO_V_Amplitude
  REAL(DP), PUBLIC :: MJ_D_Amplitude
  REAL(DP), PUBLIC :: MJ_V_Amplitude
  REAL(DP), PUBLIC :: MJ_V_Ratio

  INTEGER,  PUBLIC :: Rand_Points
  INTEGER,  PUBLIC, ALLOCATABLE :: Rand_X1(:), Rand_X2(:), Rand_X3(:)
  INTEGER,  PUBLIC, ALLOCATABLE :: Rand_NodeX(:)
  REAL(DP), PUBLIC :: Rand_Amplitude

CONTAINS


  SUBROUTINE InitializeShellPerturbations &
               ( PerturbType_Option, NumShells_Option, &
                 ShellWidth_Option, BD_Order_Option,  &
                 BD_Amplitude_Option, CO_V_Nodes_Option, &
                 CO_V_Amplitude_Option, MJ_D_Nodes_Option, &
                 MJ_D_Amplitude_Option, MJ_V_RadNodes_Option, &
                 MJ_V_AngNodes_Option, MJ_V_Amplitude_Option, &
                 MJ_V_Ratio_Option, UseUnits_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL  :: PerturbType_Option

    LOGICAL,  INTENT(in), OPTIONAL :: UseUnits_Option
    INTEGER,  INTENT(in), OPTIONAL :: NumShells_Option
    INTEGER,  INTENT(in), OPTIONAL :: BD_Order_Option
    INTEGER,  INTENT(in), OPTIONAL :: CO_V_Nodes_Option
    INTEGER,  INTENT(in), OPTIONAL :: MJ_D_Nodes_Option
    INTEGER,  INTENT(in), OPTIONAL :: MJ_V_RadNodes_Option
    INTEGER,  INTENT(in), OPTIONAL :: MJ_V_AngNodes_Option
    REAL(DP), INTENT(in), OPTIONAL :: BD_Amplitude_Option
    REAL(DP), INTENT(in), OPTIONAL :: CO_V_Amplitude_Option
    REAL(DP), INTENT(in), OPTIONAL :: MJ_D_Amplitude_Option
    REAL(DP), INTENT(in), OPTIONAL :: MJ_V_Amplitude_Option
    REAL(DP), INTENT(in), OPTIONAL :: MJ_V_Ratio_Option

    REAL(DP), INTENT(in), OPTIONAL :: ShellWidth_Option

    IF( PRESENT( PerturbType_Option ) ) THEN
      PerturbType = PerturbType_Option
    ELSE
      PerturbType = 'None'
    END IF

    IF( PRESENT( NumShells_Option ) ) THEN
      NumShells = NumShells_Option
    ELSE
      NumShells = 1
    END IF

    IF( PRESENT( ShellWidth_Option ) ) THEN
      ShellWidth = ShellWidth_Option
    ELSE
      ShellWidth = 0.4_DP
    END IF

    IF( PRESENT( UseUnits_Option ) ) THEN
      UseUnits = UseUnits_Option
    ELSE
      UseUnits = .FALSE.
    END IF

    IF( UseUnits ) THEN
      ShellWidth = ShellWidth * Kilometer
    END IF

    SELECT CASE ( TRIM( PerturbType ) )

      CASE ( 'BasicDensity' )

        IF( PRESENT ( BD_Order_Option ) ) THEN
          BD_Order = BD_Order_Option
        ELSE
          BD_Order = 0
        END IF

        IF( PRESENT( BD_Amplitude_Option ) ) THEN
          BD_Amplitude = BD_Amplitude_Option
        ELSE
          BD_Amplitude = One
        END IF

        IF( UseUnits ) THEN
          BD_Amplitude = BD_Amplitude * Gram / Centimeter**3
        END IF

      CASE ( 'CouchOtt_Velocity' )

       IF( PRESENT ( CO_V_Nodes_Option ) ) THEN
         CO_V_Nodes = CO_V_Nodes_Option
       ELSE
         CO_V_Nodes = 1
       END IF

       IF( PRESENT ( CO_V_Amplitude_Option ) ) THEN
         CO_V_Amplitude = CO_V_Amplitude_Option
       ELSE
         CO_V_Amplitude = One
       END IF

     CASE ( 'MullerJanka_Density' )

       IF( PRESENT( MJ_D_Nodes_Option ) ) THEN
         MJ_D_Nodes = MJ_D_Nodes_Option
       ELSE
         MJ_D_Nodes = 1
       END IF

       IF( PRESENT( MJ_D_Amplitude_Option ) ) THEN
         MJ_D_Amplitude = MJ_D_Amplitude_Option
       ELSE
         MJ_D_Amplitude = One
       END IF

       IF ( UseUnits ) THEN
         MJ_D_Amplitude = MJ_D_Amplitude * Gram / Centimeter**3
       END IF

     CASE ( 'MullerJanka_Velocity' )

       IF( PRESENT( MJ_V_AngNodes_Option ) ) THEN
         MJ_V_AngNodes = MJ_V_AngNodes_Option
       ELSE
         MJ_V_AngNodes = 1
       END IF

       IF( PRESENT( MJ_V_RadNodes_Option ) ) THEN
         MJ_V_RadNodes = MJ_V_RadNodes_Option
       ELSE
         MJ_V_RadNodes = 1
       END IF

       IF( PRESENT ( MJ_V_Amplitude_Option ) ) THEN
         MJ_V_Amplitude = MJ_V_Amplitude_Option
       ELSE
         MJ_V_Amplitude = One
       END IF

       IF( PRESENT ( MJ_V_Ratio_Option ) ) THEN
         MJ_V_Ratio = MJ_V_Ratio_Option
       ELSE
         MJ_V_Ratio = Half
       END IF

      CASE DEFAULT

    END SELECT

  END SUBROUTINE InitializeShellPerturbations


  SUBROUTINE InitializeRandPerturbations &
               ( Size_X1, Size_X2, Size_X3, Size_NodeX, &
                 Rand_Points_Option, Rand_Amplitude_Option )

    INTEGER,  INTENT(in) :: Size_X1, Size_X2, Size_X3
    INTEGER,  INTENT(in) :: Size_NodeX
    INTEGER,  INTENT(in), OPTIONAL :: Rand_Points_Option
    REAL(DP), INTENT(in), OPTIONAL :: Rand_Amplitude_Option

    INTEGER  :: Size1, Size2, Size3, SizeN
    INTEGER  :: iPoint
    REAL(DP) :: Rand1, Rand2, Rand3
    REAL(DP) :: RandN

    IF( PRESENT( Rand_Points_Option ) ) THEN
      Rand_Points = Rand_Points_Option
    ELSE
      Rand_Points = 1
    END IF

    DO iPoint = 1, Rand_Points

      CALL RANDOM_SEED  ( SIZE = Size1 )
      CALL RANDOM_NUMBER( Rand1 )

      CALL RANDOM_SEED  ( SIZE = Size2 )
      CALL RANDOM_NUMBER( Rand2 )

      CALL RANDOM_SEED  ( SIZE = Size3 )
      CALL RANDOM_NUMBER( Rand3 )

      CALL RANDOM_SEED  ( SIZE = SizeN )
      CALL RANDOM_NUMBER( RandN )

      Rand_X1(iPoint) = CEILING( Rand1 * ( Size_X1 ) )
      Rand_X2(iPoint) = CEILING( Rand2 * ( Size_X2 ) )
      Rand_X3(iPoint) = CEILING( Rand3 * ( Size_X3 ) )

    END DO

  END SUBROUTINE InitializeRandPerturbations


  SUBROUTINE ApplyShellPerturbations &
               ( Time, iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    REAL(DP), INTENT(in) :: Time

    INTEGER,  INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    INTEGER  :: GhostNode_Min, GhostNode_Max
    INTEGER  :: iShell
    REAL(DP) :: X1, X2
    REAL(DP) :: GhostRadius_Min, GhostRadius_Max
    REAL(DP) :: &
      P(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nPF)
    REAL(DP) :: &
      A(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nAF)

    IF( PerturbType .EQ. TRIM( 'None' ) )THEN
      RETURN
    END IF

    GhostNode_Min = NodeNumberTableX(1,1)
    GhostNode_Max = NodeNumberTableX(1,nDOFX)

    GhostRadius_Min = NodeCoordinate( MeshX(1), iX_E0(1)+1, GhostNode_Min )
    GhostRadius_Max = NodeCoordinate( MeshX(1), iX_E1(1),   GhostNode_Max )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_E0(1)+1, iX_E1(1)

      CALL ComputePrimitive_Euler_NonRelativistic &
             ( U(:,iX1,iX2,iX3,iCF_D ), U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E ), U(:,iX1,iX2,iX3,iCF_Ne), &
               P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_V1), &
               P(:,iX1,iX2,iX3,iPF_V2), P(:,iX1,iX2,iX3,iPF_V3), &
               P(:,iX1,iX2,iX3,iPF_E ), P(:,iX1,iX2,iX3,iPF_Ne), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeAuxiliary_Fluid &
             ( P(:,iX1,iX2,iX3,iPF_D  ), &
               P(:,iX1,iX2,iX3,iPF_E  ), &
               P(:,iX1,iX2,iX3,iPF_Ne ), &
               A(:,iX1,iX2,iX3,iAF_P  ), &
               A(:,iX1,iX2,iX3,iAF_T  ), &
               A(:,iX1,iX2,iX3,iAF_Ye ), &
               A(:,iX1,iX2,iX3,iAF_S  ), &
               A(:,iX1,iX2,iX3,iAF_E  ), &
               A(:,iX1,iX2,iX3,iAF_Gm ), &
               A(:,iX1,iX2,iX3,iAF_Cs ) )

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        DO iShell = 1, NumShells

          IF ( ( X1 - GhostRadius_Min ) &
               + ABS(P(iNodeX,iX1,iX2,iX3,iPF_V1)) * Time <= ShellWidth ) THEN

            SELECT CASE( TRIM( PerturbType ) )

              CASE ( 'BasicDensity' )

                CALL ApplyShellPerturbations_BasicDensity &
                       ( Time, BD_Order, BD_Amplitude, &
                         X2, P(iNodeX,iX1,iX2,iX3,iPF_D) )

              CASE( 'CouchOtt_Velocity' )

                CALL ApplyShellPerturbations_CouchOtt_Velocity &
                       ( Time, GhostRadius_Min, &
                         ShellWidth, &
                         CO_V_Nodes, CO_V_Amplitude, &
                         X1, X2, &
                         P(iNodeX,iX1,iX2,iX3,iPF_V1), &
                         P(iNodeX,iX1,iX2,iX3,iPF_V2) )

              CASE( 'MullerJanka_Velocity' )

                CALL ApplyShellPerturbations_MullerJanka_Velocity &
                       ( Time, GhostRadius_Min, &
                         ShellWidth, &
                         MJ_V_RadNodes, MJ_V_AngNodes, &
                         MJ_V_Amplitude, MJ_V_Ratio,   &
                         X1, X2, &
                         P(iNodeX,iX1,iX2,iX3,iPF_D),  &
                         P(iNodeX,iX1,iX2,iX3,iPF_V1), &
                         P(iNodeX,iX1,iX2,iX3,iPF_V2) )

            END SELECT

          END IF

        END DO

      END DO

      CALL ComputeConserved_Euler_NonRelativistic &
             ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_V1), &
               P(:,iX1,iX2,iX3,iPF_V2), P(:,iX1,iX2,iX3,iPF_V3), &
               P(:,iX1,iX2,iX3,iPF_E ), P(:,iX1,iX2,iX3,iPF_Ne), &
               U(:,iX1,iX2,iX3,iCF_D ), U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E ), U(:,iX1,iX2,iX3,iCF_Ne), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyShellPerturbations


  SUBROUTINE ApplyShellPerturbations_BasicDensity &
               ( Time, Order, Amplitude, theta, D )

    INTEGER,  INTENT(in) :: Order
    REAL(DP), INTENT(in) :: Time
    REAL(DP), INTENT(in) :: Amplitude
    REAL(DP), INTENT(in) :: theta

    REAL(DP), INTENT(inout) :: D

    SELECT CASE( Order )

      CASE ( 0 )

        D = ( One + Amplitude ) * D

      CASE ( 1 )

        D = ( One + Amplitude * COS( theta ) ) * D

    END SELECT

  END SUBROUTINE ApplyShellPerturbations_BasicDensity


  SUBROUTINE ApplyShellPerturbations_CouchOtt_Velocity &
               ( Time, GhostRadius_Min, ShellWidth, n, &
                 Amplitude, r, theta, V1, V2 )

    INTEGER,  INTENT(in)  :: n
    REAL(DP), INTENT(in)  :: Time
    REAL(DP), INTENT(in)  :: GhostRadius_Min
    REAL(DP), INTENT(in)  :: ShellWidth
    REAL(DP), INTENT(in)  :: Amplitude
    REAL(DP), INTENT(in)  :: r, theta
    REAL(DP), INTENT(in)  :: V1
    REAL(DP), INTENT(out) :: V2

    REAL(DP) :: Zeta

    ! --- Implementation of the theta velocity ---
    ! --- perturbations in Couch & Ott (2013). ---

    Zeta = Pi * ( ( r - GhostRadius_Min ) + ABS(V1) * Time ) / ( ShellWidth )

    V2 = ( Amplitude * ( One / r ) * V1 &
             * SIN( ( n - 1 ) * theta ) &
             * SIN( ( n - 1 ) * Zeta ) )

  END SUBROUTINE ApplyShellPerturbations_CouchOtt_Velocity


  SUBROUTINE ApplyShellPerturbations_MullerJanka_Velocity &
               ( Time, GhostRadius_Min, ShellWidth, &
                 n, l, Amplitude, Ratio, r, theta, D, V1, V2 )

    ! --- Still in progress. ---

    INTEGER,  INTENT(in) :: n, l
    REAL(DP), INTENT(in) :: Time
    REAL(DP), INTENT(in) :: GhostRadius_Min
    REAL(DP), INTENT(in) :: r, theta
    REAL(DP), INTENT(in) :: Amplitude
    REAL(DP), INTENT(in) :: Ratio
    REAL(DP), INTENT(in) :: D
    REAL(DP), INTENT(in) :: ShellWidth

    REAL(DP), INTENT(inout) :: V1, V2

    REAL(DP) :: V1_Perturb, V2_Perturb

    ! --- Implementation of the solenoidal velocity ---
    ! --- perturbations in Muller & Janka (2015).   ---

    !PRINT*, 'Amplitude: ', Amplitude

    !PRINT*, 'Ratio: ', Ratio

    !PRINT*, 'Density: ', D

    !PRINT*, 'r: ', r

    !PRINT*, 'theta: ', theta

    !PRINT*

    !PRINT*

    !PRINT*, 'V1 before perturbation: ', V1

    !PRINT*, 'V2 before perturbation: ', V2

    V1_Perturb = Amplitude * V1 * SIN( n * Pi * ( r - ShellWidth ) &
                                       / ( ShellWidth ) ) &
                * ( ( 3.0_DP / 2.0_DP ) * ( COS( theta ) / SQRT( SIN ( theta ) ) ) &
                * SphericalHarmonic_l1( l, theta ) + SQRT( SIN( theta ) ) &
                * SphericalHarmonic_l1_Derivative( l, theta ) ) !/ ( D * r**2 )

    V2_Perturb = - Amplitude * ( V1 / Ratio ) * ( ( SQRT( SIN( theta ) ) ) &
               * COS( n * Pi * ( r - ShellWidth ) &
                      / ( ShellWidth ) ) &
               * ( ( n * Pi ) ) & !/ ( ShellBounds(2) - ShellBounds(1) ) ) &
               * ( SphericalHarmonic_l1( l, theta ) ) ) * ( One / r )

    V1 = V1 + V1_Perturb

    V2 = V2 + V2_Perturb

    !PRINT*, 'V1 after perturbation: ', V1

    !PRINT*, 'V2 after perturbation: ', V2

    !PRINT*, 'Ratio of radial and theta velocities: ', V1_Perturb / ( r * V2_Perturb )

    !PRINT*, 'Radial kinetic energy density of perturbation: ', Half * D * V1_Perturb**2

    !PRINT*, 'Theta kinetic energy of perturbation: ', Half * D * r**2 * V2_Perturb**2

    !PRINT*, 'Ratio of kinetic energies: ', ( Half * D * V1_Perturb**2 ) / ( Half * D * r**2 * V2_Perturb**2 )

  END SUBROUTINE ApplyShellPerturbations_MullerJanka_Velocity


  SUBROUTINE ApplyRandPerturbations( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    ! --- Still in progress. ---

    INTEGER, INTENT(in) :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER :: iX1, iX2, iX3

    REAL(DP) :: &
      P(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nPF)
    REAL(DP) :: &
      A(1:nDOFX,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),1:nAF)

    ! --- Applying the random radial velocity      ---
    ! --- perturbations from Muller & Janka (2017) ---

    CALL ComputeFromConserved_Euler_NonRelativistic &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, P, A )

    P(Rand_NodeX(:),Rand_X1(:),Rand_X2(:),Rand_X3(:),iPF_V1) &
      = ( One + Rand_Amplitude ) &
        * P(Rand_NodeX(:),Rand_X1(:),Rand_X2(:),Rand_X3(:),iPF_V1)

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      CALL ComputeConserved_Euler_NonRelativistic &
               ( P(:,iX1,iX2,iX3,iPF_D ), P(:,iX1,iX2,iX3,iPF_V1), &
                 P(:,iX1,iX2,iX3,iPF_V2), P(:,iX1,iX2,iX3,iPF_V3), &
                 P(:,iX1,iX2,iX3,iPF_E ), P(:,iX1,iX2,iX3,iPF_Ne), &
                 U(:,iX1,iX2,iX3,iCF_D ), U(:,iX1,iX2,iX3,iCF_S1), &
                 U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
                 U(:,iX1,iX2,iX3,iCF_E ), U(:,iX1,iX2,iX3,iCF_Ne), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE ApplyRandPerturbations


  INTEGER RECURSIVE FUNCTION FACTORIAL( n ) RESULT( FACT )

    INTEGER, INTENT(in) :: n

    IF( n .LT. 0 ) THEN
     STOP "n < 0!"
    ELSE IF( ( n .EQ. 0 ) .OR. ( n .EQ. 1 ) ) THEN
      FACT = 1
    ELSE
      FACT = n * FACTORIAL( n - 1 )
    END IF
    RETURN
  END FUNCTION FACTORIAL


  REAL(DP) RECURSIVE FUNCTION &
    AssociatedLegendrePolynomial( l, mu ) RESULT( ALP_l1 )

    INTEGER,  INTENT(in) :: l
    REAL(DP), INTENT(in) :: mu

    ! --- Calculating Associated Legendre Polynomial l1 based on relation in ---
    ! --- http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html     ---

    IF( l .LT. 0 ) THEN
      STOP "l < 0!"
    ELSE IF( l .EQ. 0 ) THEN
      ALP_l1 = Zero
    ELSE IF( l .EQ. 1 ) THEN
      ALP_l1 = SQRT( One - mu**2 )
    ELSE
      ALP_l1 = ( mu * ( 2.0_DP * l - One ) / ( l - One ) ) &
                  * AssociatedLegendrePolynomial( l - 1, mu ) &
                  - ( l / ( l - One ) ) &
                  * AssociatedLegendrePolynomial( l - 2, mu )
    END IF
    RETURN
  END FUNCTION AssociatedLegendrePolynomial


  REAL(DP) FUNCTION &
    AssociatedLegendrePolynomial_Derivative( l, mu ) RESULT ( DALP_l1 )

    INTEGER,  INTENT(in) :: l
    REAL(DP), INTENT(in) :: mu

    ! --- Calculates theta derivative of Associated Legendre Polynomial  ---
    ! --- for l, m = 1 based on recurrance relation found at             ---
    ! --- http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html ---

    IF ( l .LT. 1 ) THEN
       STOP "l < 1"
    ELSE
      DALP_l1 = ( l * mu * AssociatedLegendrePolynomial( l, mu ) &
                  - ( l + 1 ) * AssociatedLegendrePolynomial( l - 1, mu ) ) &
                / SQRT( One - mu**2 )

    END IF
    RETURN

  END FUNCTION AssociatedLegendrePolynomial_Derivative


  REAL(DP) FUNCTION &
    SphericalHarmonic_l1( l, theta ) RESULT( SH_l1 )

    ! --- Calculates the l, m = 1 Spherical Harmonic as ---
    ! --- given in Muller & Janka (2015) when phi = 0. ---

    INTEGER,  INTENT(in) :: l
    REAL(DP), INTENT(in) :: theta

    IF ( l .LT. 1 ) THEN
      STOP "l < 1"
    ELSE
      SH_l1 = SQRT( ( ( 2 * l + 1 ) / FourPi ) &
                    * FACTORIAL( l - 1 ) &
                    / FACTORIAL( l + 1 ) ) &
                    * AssociatedLegendrePolynomial( l, COS( theta ) )
    END IF
    RETURN
  END FUNCTION SphericalHarmonic_l1


  REAL(DP) FUNCTION &
    SphericalHarmonic_l1_Derivative( l, theta ) RESULT( DSH_l1 )

    ! --- Calculates the theta derivative of the l, m = 1 Spherical  ---
    ! --- Harmonic as given in Muller & Janka (2015) when phi = 0. ---

    INTEGER,  INTENT(in) :: l
    REAL(DP), INTENT(in) :: theta

    IF ( l .LT. 1 ) THEN
      STOP "l < 1"
    ELSE
      DSH_l1 = SQRT( ( ( 2 * l + 1 ) / FourPi ) &
                     * FACTORIAL( l - 1 ) &
                     / FACTORIAL( l + 1 ) ) &
                     * AssociatedLegendrePolynomial_Derivative( l, COS( theta ) )
    END IF
    RETURN
  END FUNCTION SphericalHarmonic_l1_Derivative


  SUBROUTINE FinalizeRandPerturbations

    DEALLOCATE( Rand_X1 )
    DEALLOCATE( Rand_X2 )
    DEALLOCATE( Rand_X3 )

  END SUBROUTINE FinalizeRandPerturbations


END MODULE Euler_PerturbationModule
