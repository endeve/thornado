MODULE InitializationModule_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Three, Pi, TwoPi, FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nDOFX, &
    iX_B0, iX_B1, iX_E0, iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uAF, iAF_P
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeConserved_GR
  USE UnitsModule, ONLY: &
    Meter, Kilogram, Second, Joule
  USE UtilitiesModule, ONLY: &
    Locate

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_GR
  PUBLIC :: ReadParameters

CONTAINS

  SUBROUTINE InitializeFields_GR &
               ( RiemannProblemName_Option, &
                 RiemannProblem2dName_Option, &
                 SphericalRiemannProblemName_Option, &
                 nDetCells_Option, Eblast_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblemName_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblem2dName_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: SphericalRiemannProblemName_Option
    INTEGER,  INTENT(in), OPTIONAL         :: nDetCells_Option
    REAL(DP), INTENT(in), OPTIONAL         :: Eblast_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'RiemannProblem' )

        CALL InitializeFields_GR_RiemannProblem &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option, &
                 nDetCells_Option &
                   = nDetCells_Option, &
                 Eblast_Option &
                   = Eblast_Option )

      CASE( 'RiemannProblem2d' )

        CALL InitializeFields_GR_RiemannProblem2d &
               ( RiemannProblem2dName_Option &
                   = RiemannProblem2dName_Option )

      CASE( 'SphericalRiemannProblem' )

        CALL InitializeFields_GR_SphericalRiemannProblem &
               ( SphericalRiemannProblemName_Option &
                   = SphericalRiemannProblemName_Option )

      CASE( 'SphericalSedov' )

        CALL InitializeFields_GR_SphericalSedov &
               ( nDetCells_Option = nDetCells_Option, &
                 Eblast_Option = Eblast_Option)
        
      CASE( 'StandingAccretionShock' )

        CALL InitializeFields_GR_StandingAccretionShock
        
    END SELECT 

  END SUBROUTINE InitializeFields_GR


  SUBROUTINE InitializeFields_GR_RiemannProblem &
               ( RiemannProblemName_Option, &
                 nDetCells_Option, Eblast_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
         RiemannProblemName_Option
    INTEGER,  INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP), INTENT(in), OPTIONAL :: Eblast_Option

    CHARACTER(32) :: RiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    INTEGER       :: nDetCells
    REAL(DP)      :: Eblast, X_D

    RiemannProblemName = 'Sod'
    IF( PRESENT( RiemannProblemName_Option ) ) &
       RiemannProblemName = TRIM( RiemannProblemName_Option )

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    IF( TRIM( RiemannProblemName ) .EQ. 'CartesianSedov' )THEN
      nDetCells = 1
      IF( PRESENT( nDetCells_Option ) ) nDetCells = nDetCells_Option

      Eblast = 1.0d0
      IF( PRESENT( Eblast_Option ) ) Eblast = Eblast_Option

      X_D = DBLE( nDetCells ) * MeshX(1) % Width(1)
      WRITE(*,*)
      WRITE(*,'(A,I4.4)')      '     nDetCells:              ', nDetCells
      WRITE(*,'(A,ES23.16E3)') '     Initial blast radius:   ', X_D
      WRITE(*,'(A,ES23.16E3)') '     Initial blast pressure: ', &
                                       ( Gamma_IDEAL - One ) &
                                         * Eblast / X_D**3
    END IF

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            SELECT CASE ( TRIM( RiemannProblemName ) )

              CASE( 'Sod' )

                IF( X1 .LE. Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.125_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.1_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                END IF

              CASE( 'MBProblem1' )

                IF( X1 .LE. Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.9_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 10.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                END IF

              CASE( 'MBProblem4' )

                IF( X1 .LE. Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d3
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-2
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                END IF

              CASE( 'PerturbedShockTube' )

                IF( X1 .LE. Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 5.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 50.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
                    = 2.0_DP + 0.3_DP * SIN( 50.0_DP * X1 )
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 5.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                END IF

              CASE( 'CartesianSedov' )

                IF( X1 .LE. X_D )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = Eblast / X_D**3
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                    = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-6
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

               END IF

              CASE( 'ShockReflection' )

                IF( X1 .LE. One )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.99999_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.01_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                END IF


              CASE DEFAULT

                WRITE(*,*)
                WRITE(*,'(A,A)') &
                  'Invalid choice for RiemannProblemName: ', RiemannProblemName
                WRITE(*,'(A)') 'Valid choices:'
                WRITE(*,'(A)') &
                  "  'Sod' - &
                  Sod's shock tube"
                WRITE(*,'(A)') &
                  "  'MBProblem1' - &
                  Mignone & Bodo (2005) MNRAS, 364, 126, Problem 1"
                WRITE(*,'(A)') &
                  "  'MBProblem4' - &
                  Mignone & Bodo (2005) MNRAS, 364, 126, Problem 4"
                WRITE(*,'(A)') &
                  "  'PerturbedShockTube' - &
                  Del Zanna & Bucciantini (2002) AA, 390, 1177, &
                  Sinusoidal density perturbation"
                WRITE(*,'(A)') &
                  "  'CartesianSedov' - &
                  ..."
                WRITE(*,'(A)') &
                  "  'ShockReflection' - &
                  Del Zanna & Bucciantini (2002) AA, 390, 1177, &
                  Planar shock reflection"
                WRITE(*,'(A)') 'Stopping...'
                STOP

              END SELECT

            END DO

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
                   uAF(:,iX1,iX2,iX3,iAF_P) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFields_GR_RiemannProblem



  SUBROUTINE InitializeFields_GR_RiemannProblem2d &
               ( RiemannProblem2dName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblem2dName_Option

    CHARACTER(32) :: RiemannProblem2dName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1, iNodeX2
    REAL(DP)      :: X1, X2

    RiemannProblem2dName = 'DzB2002'
    IF( PRESENT( RiemannProblem2dName_Option ) )&
      RiemannProblem2dName = TRIM( RiemannProblem2dName_Option )

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem 2D Name: ', TRIM( RiemannProblem2dName )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)
            iNodeX2 = NodeNumberTableX(2,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

            SELECT CASE ( TRIM( RiemannProblem2dName ) )

              CASE( 'DzB2002' )

                ! --- SW ---
                IF( X1 .LE. Half .AND. X2 .LE. Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.5_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                ! --- NW ---
                ELSE IF( X1 .LE. Half .AND. X2 .GT. Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.1_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.99_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                ! --- NE ---
                ELSE IF( X1 .GT. Half .AND. X2 .GT. Half )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.1_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.01_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                ! --- SE ---
                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.1_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.99_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                END IF


              CASE DEFAULT

                WRITE(*,*)
                WRITE(*,'(A,A)') &
                  'Invalid choice for RiemannProblem2dName: ', &
                    RiemannProblem2dName
                WRITE(*,'(A)') 'Valid choices:'
                WRITE(*,'(A)') &
                  "  'DzB2002' - &
                  Del-Zanna & Bucciantini, 2D Riemann Problem"
                WRITE(*,'(A)') 'Stopping...'
                STOP

            END SELECT

          END DO

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
                   uAF(:,iX1,iX2,iX3,iAF_P) )

        END DO
      END DO
    END DO


  END SUBROUTINE InitializeFields_GR_RiemannProblem2d



  SUBROUTINE InitializeFields_GR_SphericalRiemannProblem &
               ( SphericalRiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
         SphericalRiemannProblemName_Option

    CHARACTER(32) :: SphericalRiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    SphericalRiemannProblemName = 'SphericalSod'
    IF( PRESENT( SphericalRiemannProblemName_Option ) ) &
       SphericalRiemannProblemName = TRIM( SphericalRiemannProblemName_Option )

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Spherical Riemann Problem Name: ', TRIM( SphericalRiemannProblemName )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            SELECT CASE ( TRIM( SphericalRiemannProblemName ) )

              CASE( 'SphericalSod' )

                IF( X1 <= One )THEN

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.125_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.1_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

                END IF

             CASE DEFAULT

                WRITE(*,*)
                WRITE(*,*) &
                  'Invalid choice for SphericalRiemannProblemName: ', &
                  SphericalRiemannProblemName
                WRITE(*,*) 'Valid choices:'
                WRITE(*,*) &
                  "'SphericalSod' - ", &
                  "Spherical Sod's shock tube"
                STOP

              END SELECT

            END DO

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
                   uAF(:,iX1,iX2,iX3,iAF_P) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFields_GR_SphericalRiemannProblem


  SUBROUTINE InitializeFields_GR_SphericalSedov &
               ( nDetCells_Option, Eblast_Option )
    
    INTEGER,  INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP), INTENT(in), OPTIONAL :: Eblast_Option

    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    INTEGER       :: nDetCells
    REAL(DP)      :: Eblast, X_D

    nDetCells = 1
    IF( PRESENT( nDetCells_Option ) ) nDetCells = nDetCells_Option

    Eblast = 1.0d0
    IF( PRESENT( Eblast_Option ) ) Eblast = Eblast_Option

    X_D = DBLE( nDetCells ) * MeshX(1) % Width(1)
    WRITE(*,*)
    WRITE(*,'(A,I4.4)')      '     nDetCells:              ', nDetCells
    WRITE(*,'(A,ES23.16E3)') '     Initial blast radius:   ', X_D
    WRITE(*,'(A,ES23.16E3)') '     Blast energy:           ', Eblast
    WRITE(*,'(A,ES23.16E3)') '     Initial blast pressure: ', &
                                     ( Gamma_IDEAL - One ) &
                                       * Eblast / ( FourPi / Three * X_D**3 )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            IF( X1 <= X_D)THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = Eblast / ( FourPi / Three * X_D**3 )
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = 1.0d-5
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

            END IF

          END DO

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
                   uAF(:,iX1,iX2,iX3,iAF_P) )

        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFields_GR_SphericalSedov


  SUBROUTINE InitializeFields_GR_StandingAccretionShock

    REAL(DP) :: D, V(3), P
    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    INTEGER, PARAMETER :: i_r = 1, i_D = 2, i_V1 = 3, i_E = 4
    INTEGER  :: iL, nLines
    REAL(DP) :: X1
    REAL(DP), ALLOCATABLE :: FluidFieldData(:,:), FluidFieldParameters(:)

    CALL ReadParameters &
           ( '../StandingAccretionShock_Parameters.dat', FluidFieldParameters )
    CALL ReadData &
           ( '../StandingAccretionShock_Data.dat', nLines, FluidFieldData )

    ! --- Interpolate initial conditions onto grid ---

    ! --- Loop over all elements ---
    DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          ! --- Loop over all nodes in an element ---
          DO iNodeX = 1, nDOFX

            ! --- Isolate node in X1 direction ---
            iNodeX1 = NodeNumberTableX(1,iNodeX)

            ! --- Physical coordinate corresponding to iNodeX1 ---
            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            ! --- Get lower index of input array
            !     (FluidFieldData) corresponding to physical coordinate (X1) ---
            iL = Locate( X1, FluidFieldData(:,i_r), nLines )

            ! --- Interpolate to the physical point X1 ---

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_D, i_r, iL, X1, FluidFieldData )

            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_V1, i_r, iL, X1, FluidFieldData )

            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_E, i_r, iL, X1, FluidFieldData )

            ! --- Compute pressure from internal energy density ---
            uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
              = ( Gamma_IDEAL - 1.0_DP ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

          END DO ! --- Loop over nodes ---

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33),                      &
                   uAF(:,iX1,iX2,iX3,iAF_P ) )

        END DO
      END DO
    END DO

    OPEN(  100, FILE = 'SAS_IC_interp.dat' )
    iX2 = 1
    iX3 = 1
    DO iX1 = iX_B0(1), iX_E0(1)
      DO iNodeX = 1, nDOFX

        ! --- Isolate node in X1 direction ---
        iNodeX1 = NodeNumberTableX(1,iNodeX)
          
        ! --- Physical coordinate corresponding to iNodeX1 ---
        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        WRITE(100,'(ES22.16E2,1x,ES22.16E2,1x,ES23.16E2,1x,ES22.16E2)') &
             X1, uPF(iNodeX,iX1,iX2,iX3,iPF_D), &
             uPF(iNodeX,iX1,iX2,iX3,iPF_V1), uPF(iNodeX,iX1,iX2,iX3,iPF_E)

      END DO
    END DO
    CLOSE( 100 )
    STOP 'Wrote SAS_IC_interp.dat'

  END SUBROUTINE InitializeFields_GR_StandingAccretionShock


  ! --- Auxiliary functions/subroutines for SAS problem ---

  REAL(DP) FUNCTION InterpolateInitialConditionsOntoGrid &
                      (iVar, i_r, iL, X, FluidFieldData) RESULT( yInterp )

    INTEGER,  INTENT(in) :: iVar, i_r, iL
    REAL(DP), INTENT(in) :: X
    REAL(DP), INTENT(in) :: FluidFieldData(:,:)
    REAL(DP)             :: X1, X2, Y1, Y2, m

    LOGICAL :: DEBUG = .FALSE.

    X1 = FluidFieldData(iL,i_r)
    X2 = FLuidFieldData(iL+1,i_r)
    Y1 = FluidFieldData(iL,iVar)
    Y2 = FluidFieldData(iL+1,iVar)

    m = ( Y2 - Y1 ) / ( X2 - X1 )

    ! --- Using only lower limit for slope ---
    yInterp = m * ( X - X1 ) + Y1

    IF( DEBUG )THEN
      WRITE(*,'(A)') 'Debugging InterpolateInitialConditionsOntoGrid'
      WRITE(*,'(A,I1)') 'Variable: ', iVar
      WRITE(*,'(A,ES24.16E3)') 'Y1             = ', Y1
      WRITE(*,'(A,ES24.16E3)') 'Y2             = ', Y2
      WRITE(*,'(A,ES24.16E3)') 'Y2 - Y1        = ', Y2 - Y1
      WRITE(*,'(A,ES24.16E3)') 'X2 - X1        = ', X2 - X1
      WRITE(*,'(A,ES24.16E3)') 'm              = ', m
      WRITE(*,'(A,ES24.16E3)') 'm * ( X - X1 ) = ', m * ( X - X1 )
      WRITE(*,*)
    END IF

    ! --- Using average slope ---
    ! --- Only changes accuracy in 12th decimal place ---
    !yInterp = ( 2.0_DP * m * ( X - X1 ) * ( X2 - X ) + ( Y1 * X2 + Y2 * X1 ) &
    !            - X * ( Y1 + Y2 ) ) / ( X1 + X2 - 2.0_DP * X )

    RETURN
  END FUNCTION InterpolateInitialConditionsOntoGrid


  SUBROUTINE ReadParameters( FILEIN, FluidFieldParameters)

    CHARACTER( LEN = * ), INTENT(in)   :: FILEIN
    REAL(DP), INTENT(out), ALLOCATABLE :: FluidFieldParameters(:)
    INTEGER                            :: i, nParams

    ! --- Get number of parameters ---
    nParams = 0
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO
      READ( 100, *, END = 10 )
      nParams = nParams + 1
    END DO
    10 CLOSE( 100 )

    ! --- Allocate and read in parameters ---
    ALLOCATE( FluidFieldParameters(nParams) )
    
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO i = 1, nParams
       READ( 100, '(ES23.16E2)' ) FluidFieldParameters(i)
    END DO
    CLOSE( 100 )

    ! --- Convert from physical-units to code-units ---
    FluidFieldParameters(1) = FluidFieldParameters(1) * Kilogram
    FluidFieldParameters(2) = FluidFieldParameters(2)
    FluidFieldParameters(3) = FluidFieldParameters(3) * Meter
    FluidFieldParameters(4) = FluidFieldParameters(4) * Meter
    FluidFieldParameters(5) = FluidFieldParameters(5) * Meter
    FluidFieldParameters(6) = FluidFieldParameters(6) * Meter
    FluidFieldParameters(7) = FluidFieldParameters(7) * Kilogram / Second

  END SUBROUTINE ReadParameters

  
  SUBROUTINE ReadData( FILEIN, nLines, FluidFieldData )

    CHARACTER( LEN = * ), INTENT(in)   :: FILEIN
    INTEGER,  INTENT(inout)            :: nLines
    REAL(DP), INTENT(out), ALLOCATABLE :: FluidFieldData(:,:)
    INTEGER                            :: i

    ! --- Get number of lines in data file ---
    nLines = 0
    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO
      READ( 100, *, END = 10 )
      nLines = nLines + 1
    END DO
    10 CLOSE( 100 )

    ! --- Allocate and read in data ---
    ALLOCATE( FluidFieldData( 1:nLines, 4 ) )

    OPEN( 100, FILE = TRIM( FILEIN ) )
    READ( 100, * ) ! --- Skip the header ---
    DO i = 1, nLines
       READ( 100, '(ES22.16E2,1x,ES22.16E2,1x,ES23.16E2,1x,ES22.16E2)' ) &
         FluidFieldData(i,:)
    END DO
    CLOSE( 100 )

    ! --- Convert from physical-units to code-units ---
    FluidFieldData(:,1) = FluidFieldData(:,1) * Meter
    FluidFieldData(:,2) = FluidFieldData(:,2) * Kilogram / Meter**3
    FluidFieldData(:,3) = FluidFieldData(:,3) * Meter / Second
    FluidFieldData(:,4) = FluidFieldData(:,4) * Joule / Meter**3

    OPEN(  100, FILE = 'SAS_IC_raw.dat' )
      DO i = 1, nLines
        WRITE(100,'(ES22.16E2,1x,ES22.16E2,1x,ES23.16E2,1x,ES22.16E2)') &
              FluidFieldData(i,:)
      END DO
    CLOSE( 100 )
!!$    STOP 'Wrote SAS initial conditions to file: SAS_IC.dat'

  END SUBROUTINE ReadData


END MODULE InitializationModule_GR
