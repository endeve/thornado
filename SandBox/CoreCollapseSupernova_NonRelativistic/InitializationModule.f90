MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    nDOFX, nDOFE
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
    uAF, iAF_P, iAF_Ye, iAF_T, iAF_E, iAF_S, iAF_Me, &
    iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, iAF_Xa, iAF_Xh, iAF_Gm
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ApplyEquationOfState_TABLE
  USE UtilitiesModule, ONLY: &
    Locate, &
    Interpolate1D_Linear
  USE ProgenitorModule, ONLY: &
    ProgenitorType1D, &
    ReadProgenitor1D

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields( ProgenitorFileName )

    CHARACTER(LEN=*), INTENT(in) :: ProgenitorFileName

    CALL InitializeFields_Fluid( ProgenitorFileName )

    CALL InitializeFields_Radiation

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_Fluid( ProgenitorFileName )

    CHARACTER(LEN=*), INTENT(in) :: ProgenitorFileName

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1
    TYPE(ProgenitorType1D) :: P1D

    CALL ReadProgenitor1D( TRIM( ProgenitorFileName ), P1D )

    ASSOCIATE &
      ( R1D => P1D % Radius, &
        D1D => P1D % MassDensity, &
        V1D => P1D % RadialVelocity, &
        T1D => P1D % Temperature, &
        Y1D => P1D % ElectronFraction )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
          = Interpolate1D( R1D, D1D, SIZE( R1D ), X1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = Interpolate1D( R1D, V1D, SIZE( R1D ), X1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = Zero

        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = Zero

        uAF(iNodeX,iX1,iX2,iX3,iAF_T) &
          = Interpolate1D( R1D, T1D, SIZE( R1D ), X1 )

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
          = Interpolate1D( R1D, Y1D, SIZE( R1D ), X1 )

        CALL ComputeThermodynamicStates_Primitive_TABLE &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) )

        CALL ApplyEquationOfState_TABLE &
                ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_S ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_E ), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Me), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Mp), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Mn), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xp), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xn), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xa), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Xh), &
                  uAF(iNodeX,iX1,iX2,iX3,iAF_Gm) )

        CALL ComputeConserved_Euler_NonRelativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE ! R1D, etc.

  END SUBROUTINE InitializeFields_Fluid


  SUBROUTINE InitializeFields_Radiation

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeE, iNodeX, iNodeZ
    REAL(DP) :: Gm_dd_11(nDOFX)
    REAL(DP) :: Gm_dd_22(nDOFX)
    REAL(DP) :: Gm_dd_33(nDOFX)

    DO iS = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Gm_dd_11 = uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11)
      Gm_dd_22 = uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22)
      Gm_dd_33 = uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33)

      DO iE = iE_B0, iE_E0

        DO iNodeX = 1, nDOFX
        DO iNodeE = 1, nDOFE

          iNodeZ = (iNodeX-1) * nDOFE + iNodeE

          uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) = 1.0d-40
          uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Zero
          uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
          uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

          CALL ComputeConserved_TwoMoment &
                 ( uPR(iNodeZ:iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                   uPR(iNodeZ:iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                   uPR(iNodeZ:iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                   uPR(iNodeZ:iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                   uCR(iNodeZ:iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                   uCR(iNodeZ:iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                   uCR(iNodeZ:iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                   uCR(iNodeZ:iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                   Gm_dd_11(iNodeX:iNodeX), &
                   Gm_dd_22(iNodeX:iNodeX), &
                   Gm_dd_33(iNodeX:iNodeX) )

        END DO
        END DO

      END DO

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Radiation


  REAL(DP) FUNCTION Interpolate1D( x, y, n, xq )

    INTEGER,  INTENT(in) :: n
    REAL(DP), INTENT(in) :: x(n), y(n)
    REAL(DP), INTENT(in) :: xq

    INTEGER :: i

    i = Locate( xq, x, n )

    IF( i == 0 )THEN

      ! --- Extrapolate Left ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(1), x(2), y(1), y(2) )

    ELSE IF( i == n )THEN

      ! --- Extrapolate Right ---

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(n-1), x(n), y(n-1), y(n) )

    ELSE

      Interpolate1D &
        = Interpolate1D_Linear( xq, x(i), x(i+1), y(i), y(i+1) )

    END IF

    RETURN

  END FUNCTION Interpolate1D


END MODULE InitializationModule
