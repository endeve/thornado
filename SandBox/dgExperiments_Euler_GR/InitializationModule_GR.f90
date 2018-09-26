MODULE InitializationModule_GR

  USE KindModule, ONLY: &
<<<<<<< .mine
    DP, Zero, Half, One
||||||| .r960
    DP, Zero, Half, One, TwoPi, Pi
=======
    DP, Half, One
>>>>>>> .r1035
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
<<<<<<< .mine
    nX, nNodesX, &
    nDOFX, &
    iX_B0, iX_B1, iX_E0, iX_E1
||||||| .r960
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nNodesX, &
    nDOFX
=======
    nX, nNodesX, &
    nDOFX
>>>>>>> .r1035
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
<<<<<<< .mine
  USE UnitsModule, ONLY: &
    Meter, Kilogram, Second, Joule
  USE UtilitiesModule, ONLY: &
    Locate

||||||| .r960
  USE DataFileReaderModule, ONLY: &
    ReadData, ReadParameters
  USE UtilitiesModule, ONLY: &
    Locate
  
=======

>>>>>>> .r1035
  IMPLICIT NONE
  PRIVATE

<<<<<<< .mine
  PUBLIC :: InitializeFields_GR
  PUBLIC :: ReadParameters
||||||| .r960
  PUBLIC :: InitializeFields_RiemannProblem
  PUBLIC :: InitializeFields_StandingAccretionShock
=======
  PUBLIC :: InitializeFields_GR
>>>>>>> .r1035

CONTAINS

<<<<<<< .mine
  SUBROUTINE InitializeFields_GR &
               ( RiemannProblemName_Option, &
                 SphericalRiemannProblemName_Option, &
                 nDetCells_Option, nRel_Option )
||||||| .r960
  SUBROUTINE InitializeFields_RiemannProblem &
    ( D_L, V_L, P_L, D_R, V_R, P_R, X_D_Option )
=======
  SUBROUTINE InitializeFields_GR &
               ( RiemannProblemName_Option, nDetCells_Option, nRel_Option )
>>>>>>> .r1035

<<<<<<< .mine
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblemName_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: SphericalRiemannProblemName_Option
    INTEGER,  INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP), INTENT(in), OPTIONAL :: nRel_Option
||||||| .r960
    REAL(DP), INTENT(in) :: D_L, P_L, D_R, P_R
    REAL(DP), INTENT(in) :: V_L(3), V_R(3)
    REAL(DP), INTENT(in), OPTIONAL :: X_D_Option
=======
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblemName_Option
    INTEGER,  INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP), INTENT(in), OPTIONAL :: nRel_Option
>>>>>>> .r1035

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

<<<<<<< .mine
      CASE( 'RiemannProblem' )

        CALL InitializeFields_GR_RiemannProblem &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option )

      CASE( 'SphericalRiemannProblem' )

        CALL InitializeFields_GR_SphericalRiemannProblem &
               ( SphericalRiemannProblemName_Option &
                   = SphericalRiemannProblemName_Option )

      CASE( 'SedovBlastWave' )

        CALL InitializeFields_GR_SedovBlastWave &
               ( nDetCells_Option = nDetCells_Option, &
                 nRel_Option = nRel_Option)

      CASE( 'StandingAccretionShock' )

         CALL InitializeFields_GR_StandingAccretionShock

    END SELECT 

  END SUBROUTINE InitializeFields_GR


  SUBROUTINE InitializeFields_GR_RiemannProblem &
               ( RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
         RiemannProblemName_Option

    CHARACTER(32) :: RiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    RiemannProblemName = 'Sod'
    IF( PRESENT( RiemannProblemName_Option ) ) &
       RiemannProblemName = TRIM( RiemannProblemName_Option )

||||||| .r960
=======
      CASE( 'RiemannProblem' )

        CALL InitializeFields_GR_RiemannProblem &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option )

      CASE( 'SedovBlastWave' )

        CALL InitializeFields_GR_SedovBlastWave &
               ( nDetCells_Option = nDetCells_Option, &
                 nRel_Option = nRel_Option)
        
    END SELECT 

  END SUBROUTINE InitializeFields_GR


  SUBROUTINE InitializeFields_GR_RiemannProblem &
               ( RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
         RiemannProblemName_Option

    CHARACTER(32) :: RiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    RiemannProblemName = 'Sod'
    IF( PRESENT( RiemannProblemName_Option ) ) &
       RiemannProblemName = TRIM( RiemannProblemName_Option )

>>>>>>> .r1035
    WRITE(*,*)
<<<<<<< .mine
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            SELECT CASE ( TRIM( RiemannProblemName ) )

              CASE( 'Sod' )

                IF( X1 <= Half )THEN

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

                IF( X1 <= Half )THEN

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

                IF( X1 <= Half )THEN

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

                IF( X1 <= Half )THEN

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

              CASE DEFAULT

                WRITE(*,*)
                WRITE(*,*) &
                  'Invalid choice for RiemannProblemName: ', RiemannProblemName
                WRITE(*,*) 'Valid choices:'
                WRITE(*,*) &
                  "Sod - ", &
                    "Sod's shock tube"
                WRITE(*,*) &
                  "MBProblem1 - ", &
                    "Mignone & Bodo (2005) MNRAS, 364, 126, Problem 1"
                WRITE(*,*) &
                  "MBProblem4 - ", &
                    "Mignone & Bodo (2005) MNRAS, 364, 126, Problem 4"
                WRITE(*,*) &
                  "PerturbedShockTube - ", &
                    "Del Zanna & Bucciantini (2002) AA, 390, 1177, ", &
                    "sinusoidal density perturbation"
                WRITE(*,*) 'Stopping...'
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
      '', 'Spherical Riemann Problem Name: ', &
      TRIM( SphericalRiemannProblemName )
||||||| .r960
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2)') &
      '', 'X_D = ', X_D
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
      '', 'D_L = ', D_L, '', 'D_R = ', D_R
    WRITE(*,'(A7,A6,3ES10.3E2,A4,A6,3ES10.3E2)') &
      '', 'V_L = ', V_L, '', 'V_R = ', V_R
    WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
      '', 'P_L = ', P_L, '', 'P_R = ', P_R
=======
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )
>>>>>>> .r1035

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

<<<<<<< .mine
            SELECT CASE ( TRIM( SphericalRiemannProblemName ) )
||||||| .r960
!!$            ! --- Riemann Problems ---
!!$            IF( X1 <= X_D )THEN
!!$
!!$              ! -- Left State --
!!$
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = D_L
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_L(1)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_L(2)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_L(3)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = P_L / (Gamma_IDEAL-One)
!!$              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = P_L
!!$
!!$            ELSE
!!$
!!$              ! -- Right State --
!!$
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = D_R
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_R(1)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_R(2)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_R(3)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = P_R / (Gamma_IDEAL-One)
!!$              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = P_R
!!$
!!$            END IF
=======
            SELECT CASE ( TRIM( RiemannProblemName ) )
>>>>>>> .r1035

<<<<<<< .mine
              CASE( 'SphericalSod' )
||||||| .r960
!!$            ! --- Smooth solution, advected sine wave ---
!!$
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
!!$              = 1.0_DP + 0.1_DP * SIN( TwoPi * X1 )
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0d-6 / (Gamma_IDEAL-One)
!!$            uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-6
=======
              CASE( 'Sod' )
>>>>>>> .r1035

<<<<<<< .mine
                IF( X1 <= One )THEN
||||||| .r960
!!$            ! --- Contact discontinuity, top-hat ---
!!$
!!$            IF( ( X1 <= 0.25d0 ) .OR. ( X1 >= 0.75d0 ) )THEN
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0d0
!!$            ELSE
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 2.0d0
!!$            END IF
!!$
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0d-6 / (Gamma_IDEAL-One)
!!$            uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-6
=======
                IF( X1 <= Half )THEN
>>>>>>> .r1035

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

<<<<<<< .mine
                END IF

              CASE DEFAULT

                WRITE(*,*)
                WRITE(*,*) &
                  'Invalid choice for SphericalRiemannProblemName: ', &
                  SphericalRiemannProblemName  
                WRITE(*,*) 'Valid choices:'
                WRITE(*,*) &
                  "SphericalSod - ", &
                    "Sod's shock tube, spherical coordinates"
                WRITE(*,*) 'Stopping...'
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


  SUBROUTINE InitializeFields_GR_SedovBlastWave( nDetCells_Option, nRel_Option )
    
    INTEGER,  INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP), INTENT(in), OPTIONAL :: nRel_Option

    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    INTEGER       :: nDetCells
    REAL(DP)      :: nRel, X_D

    nDetCells = 1
    IF( PRESENT( nDetCells_Option ) ) nDetCells = nDetCells_Option

    nRel = 1.0d0
    IF( PRESENT( nRel_Option ) ) nRel = nRel_Option

    X_D = DBLE( nDetCells ) * MeshX(1) % Width(1)

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
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                = nRel * uPF(iNodeX,iX1,iX2,iX3,iPF_D)
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ELSE
||||||| .r960
            ELSE
=======
                END IF
>>>>>>> .r1035

<<<<<<< .mine
              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                = 1.0d-6
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )
||||||| .r960
              ! -- Right State --
=======
              CASE( 'MBProblem1' )
>>>>>>> .r1035

<<<<<<< .mine
||||||| .r960
              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 2.0d0 + 0.3d0 * SIN( 50.0d0 * X1 )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 5.0d0 / (Gamma_IDEAL-One)
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 5.0d0

=======
                IF( X1 <= Half )THEN

>>>>>>> .r1035
                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.9_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

<<<<<<< .mine
          END DO
||||||| .r960

!!$            ! --- Acoustic problem, spherical symmetry ---
!!$
!!$            IF( ( X1 .GE. 0.4d0 ) .AND. ( X1 .LE. 0.6d0 ) )THEN
!!$              eps = 1.0d-4
!!$            ELSE
!!$              eps = 0.0d0
!!$            END IF
!!$
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
!!$              = 1.0_DP + eps * SIN( Pi * X1 )**4 / X1
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
!!$              = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
!!$              = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
!!$              = 0.0_DP
!!$            uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
!!$              = 1.0_DP / Gamma_IDEAL + eps * SIN( Pi * X1 )**4 / X1
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
!!$              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - 1.0_DP )
            
          END DO
=======
                ELSE

                  uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 10.0_DP
                  uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                    = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )
>>>>>>> .r1035

                END IF

              CASE( 'MBProblem4' )

                IF( X1 <= Half )THEN

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

                IF( X1 <= Half )THEN

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

             CASE DEFAULT

                WRITE(*,*)
                WRITE(*,*) &
                  'Invalid choice for RiemannProblemName: ', RiemannProblemName
                WRITE(*,*) 'Valid choices:'
                WRITE(*,*) &
                  "Sod's shock tube: ", &
                    "'Sod'"
                WRITE(*,*) &
                  "Mignone & Bodo (2005) MNRAS, 364, 126, Problem 1: ", &
                    "'MBProblem1'"
                WRITE(*,*) &
                  "Mignone & Bodo (2005) MNRAS, 364, 126, Problem 4: ", &
                    "'MBProblem4'"
                WRITE(*,*) &
                  "Del Zanna & Bucciantini (2002) AA, 390, 1177, ", &
                  "sinusoidal density perturbation: ", &
                    "'PerturbedShockTube'"
                WRITE(*,*) 'Stopping...'
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

<<<<<<< .mine
  END SUBROUTINE InitializeFields_GR_SedovBlastWave
||||||| .r960
  END SUBROUTINE InitializeFields_RiemannProblem
=======
  END SUBROUTINE InitializeFields_GR_RiemannProblem
>>>>>>> .r1035


<<<<<<< .mine
  SUBROUTINE InitializeFields_GR_StandingAccretionShock
||||||| .r960
  SUBROUTINE InitializeFields_StandingAccretionShock
=======
>>>>>>> .r1035

<<<<<<< .mine
    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    INTEGER, PARAMETER :: i_r = 1, i_D = 2, i_V1 = 3, i_E = 4
    INTEGER  :: iL, nLines
    REAL(DP) :: X1
    REAL(DP), ALLOCATABLE :: FluidFieldData(:,:)
    REAL(DP), ALLOCATABLE :: r(:), rho(:), v(:), e(:)
    REAL(DP) :: X1_arr(0:nX(1)+1)
||||||| .r960
    REAL(DP) :: D, V(3), P
    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    INTEGER, PARAMETER :: i_r = 1, i_D = 2, i_V1 = 3, i_E = 4
    INTEGER  :: iL, nLines, iX1_Debug = 800
    REAL(DP) :: X1
    REAL(DP), ALLOCATABLE :: FluidFieldData(:,:), FluidFieldParameters(:)
=======
  SUBROUTINE InitializeFields_GR_SedovBlastWave( nDetCells_Option, nRel_Option )
    
    INTEGER,  INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP), INTENT(in), OPTIONAL :: nRel_Option
>>>>>>> .r1035

    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

<<<<<<< .mine
    CALL ReadData &
           ( '../StandingAccretionShock_Data.dat', nLines, FluidFieldData )
||||||| .r960
    CALL ReadParameters &
           ( '../StandingAccretionShock_Parameters.dat', FluidFieldParameters )
    CALL ReadData &
           ( '../StandingAccretionShock_Data.dat', nLines, FluidFieldData )
=======
    INTEGER       :: nDetCells
    REAL(DP)      :: nRel, X_D
>>>>>>> .r1035

<<<<<<< .mine
    ALLOCATE( r  (1:nLines) )
    ALLOCATE( rho(1:nLines) )
    ALLOCATE( v  (1:nLines) )
    ALLOCATE( e  (1:nLines) )

    r   = FluidFieldData(:,1)
    rho = FluidFieldData(:,2)
    v   = FluidFieldData(:,3)
    e   = FluidFieldData(:,4)

    ! --- Interpolate initial conditions onto grid ---
    WRITE(*,'(A)') 'Interpolating initial conditions onto grid'
||||||| .r960
    ! --- Interpolate initial conditions onto grid ---
    WRITE(*,'(A)') 'Interpolating initial conditions onto grid'
=======
    nDetCells = 1
    IF( PRESENT( nDetCells_Option ) ) nDetCells = nDetCells_Option
>>>>>>> .r1035

    nRel = 1.0d0
    IF( PRESENT( nRel_Option ) ) nRel = nRel_Option

    X_D = DBLE( nDetCells ) * MeshX(1) % Width(1)

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX)

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
            X1_arr(iX1) = X1

            IF( X1 <= X_D)THEN

<<<<<<< .mine
||||||| .r960
            IF( iX1 == iX1_Debug ) iL_Debug = iL

=======
              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                = nRel * uPF(iNodeX,iX1,iX2,iX3,iPF_D)
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

>>>>>>> .r1035
            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
                = 1.0d-6
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF

          END DO

<<<<<<< .mine
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_E, i_r, iL, X1, FluidFieldData )

            ! --- Compute pressure from internal energy density ---
            uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
              = ( Gamma_IDEAL - 1.0_DP ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

||||||| .r960
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_E, i_r, iL, X1, FluidFieldData )

            ! --- Compute pressure from internal energy density ---
            uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
              = ( Gamma_IDEAL - 1.0_DP ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

          END DO ! --- Loop over nodes ---

=======
>>>>>>> .r1035
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

          END DO ! --- Loop over nodes ---

        END DO
      END DO
    END DO

<<<<<<< .mine
    DEALLOCATE( r, rho, v, e )
||||||| .r960
!!$    WRITE(*,*) iL_Debug, iX1_Debug
!!$    WRITE(*,'(A,ES22.16E2)') 'rhoL: ', uPF(1,iX1_Debug,1,1,iPF_D )
!!$    WRITE(*,'(A,ES22.16E2)') 'rhoU: ', uPF(1,iX1_Debug,1,1,iPF_D )
!!$    WRITE(*,'(A)') 'Stopping...'
!!$    STOP
=======
  END SUBROUTINE InitializeFields_GR_SedovBlastWave
>>>>>>> .r1035

<<<<<<< .mine
    ! --- Write initial conditions to file ---
    OPEN( 100, FILE = 'SAS_InitialConditions.dat' )
    DO iX1 = 1, nX(1)
      WRITE( 100, '(ES24.16E3,1x,ES24.16E3,1x,ES24.16E3,1x,ES24.16E3)' ) &
        X1_arr(iX1), uPF(:,iX1,1,1,iPF_D ), uPF(:,iX1,1,1,iPF_V1 ), &
          uAF(:,iX1,1,1,iAF_P )
    END DO
    CLOSE( 100 )
!    STOP 'Wrote initial conditions to file'
||||||| .r960
  END SUBROUTINE InitializeFields_StandingAccretionShock
=======
>>>>>>> .r1035

<<<<<<< .mine
  END SUBROUTINE InitializeFields_GR_StandingAccretionShock

  ! --- Auxiliary functions/subroutines for SAS problem ---

  REAL(DP) FUNCTION InterpolateInitialConditionsOntoGrid &
                      (iVar, i_r, iL, X, FluidFieldData) RESULT( yInterp )

    INTEGER,  INTENT(in) :: iVar, i_r, iL
    REAL(DP), INTENT(in) :: X
    REAL(DP), INTENT(in) :: FluidFieldData(:,:)
    REAL(DP)             :: m, X1, X2, Y1, Y2

    X1 = FluidFieldData(iL,i_r)
    X2 = FLuidFieldData(iL+1,i_r)
    Y1 = FluidFieldData(iL,iVar)
    Y2 = FluidFieldData(iL+1,iVar)

    m = ( Y2 - Y1 ) / ( X2 - X1 )

    ! --- Using only lower limit for slope ---
    yInterp = m * ( X - X1 ) + Y1

    ! --- Using average slope ---
    ! --- Only changes accuracy in 12th decimal place ---
!    yInterp = ( 2.0_DP * m * ( X - X1 ) * ( X2 - X ) + ( Y1 * X2 + Y2 * X1 ) &
!                - X * ( Y1 + Y2 ) ) / ( X1 + X2 - 2.0_DP * X )

    RETURN
  END FUNCTION InterpolateInitialConditionsOntoGrid

  SUBROUTINE ReadParameters( FILEIN, FluidFieldParameters )

||||||| .r960

  REAL(DP) FUNCTION InterpolateInitialConditionsOntoGrid &
                      (iVar, i_r, iL, X, FluidFieldData) RESULT( yInterp )

    INTEGER,  INTENT(in) :: iL, iVar, i_r
    REAL(DP), INTENT(in) :: X
    REAL(DP), INTENT(in) :: FluidFieldData(:,:)
    REAL(DP)             :: m, X1, X2, Y1, Y2

    X1 = FluidFieldData(iL,i_r)
    X2 = FLuidFieldData(iL+1,i_r)
    Y1 = FluidFieldData(iL,iVar)
    Y2 = FluidFieldData(iL+1,iVar)

    m = ( Y2 - Y1 ) / ( X2 - X1 )

    ! --- Using only lower limit for slope ---
    yInterp = m * ( X - X1 ) + Y1

    ! --- Using average slope ---
    ! --- Only changes accuracy in 12th decimal place ---
!    yInterp = ( 2.0_DP * m * ( X - X1 ) * ( X2 - X ) + ( Y1 * X2 + Y2 * X1 ) &
!                - X * ( Y1 + Y2 ) ) / ( X1 + X2 - 2.0_DP * X )

    IF( iL == iL_Debug  .AND. iVar == iPF_D ) THEN
      WRITE(*,'(A,ES24.16E3)') 'X1: ', X1
      WRITE(*,'(A,ES24.16E3)') 'X2: ', X2
      WRITE(*,'(A,ES24.16E3)') 'Y1: ', Y1
      WRITE(*,'(A,ES24.16E3)') 'Y2: ', Y2
      WRITE(*,'(A,ES24.16E3)') 'm:  ', m
      WRITE(*,'(A,ES24.16E3)') 'yInterp: ', yInterp
    END IF

    RETURN
  END FUNCTION InterpolateInitialConditionsOntoGrid


=======
>>>>>>> .r1035
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
       READ( 100, '(4ES23.16E2)' ) FluidFieldData(i,:)
    END DO
    CLOSE( 100 )

    ! --- Convert from physical-units to code-units ---
    FluidFieldData(:,1) = FluidFieldData(:,1) * Meter
    FluidFieldData(:,2) = FluidFieldData(:,2) * Kilogram / Meter**3
    FluidFieldData(:,3) = FluidFieldData(:,3) * Meter / Second
    FluidFieldData(:,4) = FluidFieldData(:,4) * Joule / Meter**3

  END SUBROUTINE ReadData


END MODULE InitializationModule_GR
