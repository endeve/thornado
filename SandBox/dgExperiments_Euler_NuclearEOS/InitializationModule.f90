MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Three, &
    Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    Gram, Centimeter, &
    Kilometer, Erg, Second, Kelvin
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    nNodesX, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Phi_N, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uAF, iAF_P, iAF_Ye, iAF_T, iAF_E
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    Euler_ComputeConserved_NonRelativistic
  USE EquationOfStateModule, ONLY: &
    ComputeTemperatureFromPressure, &
    ComputeThermodynamicStates_Primitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields &
    ( RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE ( 'RiemannProblem' )

        CALL InitializeFields_RiemannProblem &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option )

      CASE ( 'RiemannProblemSpherical' )

        CALL InitializeFields_RiemannProblemSpherical &
               ( RiemannProblemName_Option &
                   = RiemannProblemName_Option )

      CASE( 'Implosion' )

        CALL InitializeFields_Implosion

      CASE( 'Jet' )

        CALL InitializeFields_Jet

    END SELECT

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_RiemannProblem &
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

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'Sod' )

            IF( X1 <= Zero )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.00d12 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.4_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.25d11 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d31 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

            END IF

        END SELECT

      END DO

      CALL ComputeTemperatureFromPressure &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

      CALL Euler_ComputeConserved_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblem


  SUBROUTINE InitializeFields_RiemannProblemSpherical &
    ( RiemannProblemName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      RiemannProblemName_Option

    CHARACTER(32) :: RiemannProblemName
    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    RiemannProblemName = 'SphericalSod'
    IF( PRESENT( RiemannProblemName_Option ) ) &
       RiemannProblemName = TRIM( RiemannProblemName_Option )

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'SphericalSod' )

            IF( X1 <= One * Kilometer )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.00d12 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.4_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.25d11 * Gram / Centimeter**3
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
              uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d30 * Erg / Centimeter**3
              uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

            END IF

        END SELECT

      END DO

      CALL ComputeTemperatureFromPressure &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

      CALL Euler_ComputeConserved_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblemSpherical


  SUBROUTINE InitializeFields_Implosion

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP), PARAMETER :: D_0 = 0.125_DP
    REAL(DP), PARAMETER :: E_0 = 0.350_DP
    REAL(DP), PARAMETER :: D_1 = 1.000_DP
    REAL(DP), PARAMETER :: E_1 = 2.500_DP

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        IF( X1 + X2 .LT. 0.15_DP )THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
            = D_0
          uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
            = E_0

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
            = D_1
          uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
            = E_1

        ENDIF

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = Zero

      END DO

      CALL Euler_ComputeConserved_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Implosion


  SUBROUTINE InitializeFields_Jet

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) !* Kilometer
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 ) !* Kilometer

        IF( X1 .LE. Half * Kilometer .AND. X2 .LE. Half * Kilometer )THEN

          !SW
          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 0.80d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

        ELSE IF ( X1 .LE. Half * Kilometer .AND. X2 .GT. Half * Kilometer )THEN

          !NW
          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 7.275d4 * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

        ELSE IF( X1 .GT. Half * Kilometer .AND. X2 .GT. Half * Kilometer )THEN

          !NE
          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 0.5313d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 0.4d32 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

        ELSE

          !SE
          uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0d12 * Gram / Centimeter**3
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 7.275d4 * Kilometer / Second
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP * Kilometer / Second
          uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0d32 * Erg / Centimeter**3
          uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3_DP

        END IF

      END DO

      CALL ComputeTemperatureFromPressure &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_P), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_T) )

      CALL ComputeThermodynamicStates_Primitive &
             ( uPF(:,iX1,iX2,iX3,iPF_D),  uAF(:,iX1,iX2,iX3,iAF_T), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ),&
               uAF(:,iX1,iX2,iX3,iAF_E),  uPF(:,iX1,iX2,iX3,iPF_Ne) )

      CALL Euler_ComputeConserved_NonRelativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Jet


END MODULE InitializationModule
