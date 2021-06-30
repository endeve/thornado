MODULE InitializationModule_Neutrinos

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer, &
    MeV, &
    Kelvin, &
    Erg, &
    SpeedOfLight, &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOFX, nDOFE
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeConserved_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeConserved_TwoMoment
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeTemperatureFromPressure_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE, &
    ApplyEquationOfState_TABLE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields

CONTAINS


  SUBROUTINE InitializeFields()

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

       CASE( 'Relaxation' )

         CALL InitializeFields_Relaxation

       CASE( 'Deleptonization' )

         CALL InitializeFields_Deleptonization

       CASE( 'EquilibriumAdvection' )

         CALL InitializeFields_EquilibriumAdvection

    END SELECT

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_Relaxation

    REAL(DP), PARAMETER :: D_0   = 1.032d12 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: T_0   = 7.588d0 * MeV
    REAL(DP), PARAMETER :: Y_0   = 0.1347_DP
    REAL(DP), PARAMETER :: V_u_1 = 0.0_DP * SpeedOfLight
    REAL(DP), PARAMETER :: V_u_2 = 0.0_DP * SpeedOfLight
    REAL(DP), PARAMETER :: V_u_3 = 0.0_DP * SpeedOfLight
    REAL(DP), PARAMETER :: Mu_0  = - 1.0_DP ! \in [-1,1]

    INTEGER  :: iE, iX1, iX2, iX3, iS, iNodeE, iNodeX, iNodeZ
    REAL(DP) :: kT, E, f_E

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = D_0
        uAF(iNodeX,iX1,iX2,iX3,iAF_T ) = T_0
        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = Y_0

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

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_u_1
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_u_2
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_u_3

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

    ! --- Radiation Fields ---

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        kT = BoltzmannConstant * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

        E = NodeCoordinate( MeshE, iE, iNodeE )

        f_E = MAX( 0.99_DP * EXP( - ( E - Two*kT )**2 &
                                    / ( Two*(1.0d1*MeV)**2 ) ), 1.0d-99 )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D,iS) &
          = f_E * 0.50_DP * ( One - Mu_0 )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) &
          = f_E * 0.25_DP * ( One - Mu_0**2 )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V1),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V2),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V3),        &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_11),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_22),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Relaxation


  SUBROUTINE InitializeFields_Deleptonization

    REAL(DP), PARAMETER :: D_L = 1.0d13 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: V_L = 0.0d00 * SpeedOfLight
    REAL(DP), PARAMETER :: P_L = 1.1d32 * Erg / Centimeter**3
    REAL(DP), PARAMETER :: Y_L = 0.3d00

    REAL(DP), PARAMETER :: D_R = 1.0d12 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: V_R = 0.0d00 * SpeedOfLight
    REAL(DP), PARAMETER :: P_R = 1.0d32 * Erg / Centimeter**3
    REAL(DP), PARAMETER :: Y_R = 0.4d00

    REAL(DP), PARAMETER :: L_0 = 1.0d1 * Kilometer

    INTEGER  :: iE, iX1, iX2, iX3, iS, iNodeE, iNodeX, iNodeX1, iNodeZ
    REAL(DP) :: X1, Theta, kT, Mnu, E

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        Theta = Half * ( One - TANH( X1 / L_0 ) )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
          = Theta * D_L + (One-Theta) * D_R

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = Theta * V_L + (One-Theta) * V_R

        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = Zero

        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = Zero

        uAF(iNodeX,iX1,iX2,iX3,iAF_P ) &
          = Theta * P_L + (One-Theta) * P_R

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) &
          = Theta * Y_L + (One-Theta) * Y_R

        CALL ComputeTemperatureFromPressure_TABLE &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_Ye), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_T ) )

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

    ! --- Radiation Fields ---

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        kT = BoltzmannConstant * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

        Mnu = uAF  (iNodeX,iX1,iX2,iX3,iAF_Me) &
              + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
              - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

        IF( iS == iNuE_Bar )THEN

          Mnu = - Mnu

        END IF

        E = NodeCoordinate( MeshE, iE, iNodeE )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D,iS) &
          = One / ( EXP( ( E - Mnu ) / kT ) + One )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Zero
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V1),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V2),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V3),        &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_11),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_22),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Deleptonization


  SUBROUTINE InitializeFields_EquilibriumAdvection

    REAL(DP), PARAMETER :: D_0   = 3.0d14 * Gram / Centimeter**3
    REAL(DP), PARAMETER :: T_0   = 1.2d11 * Kelvin
    REAL(DP), PARAMETER :: Y_0   = 0.3_DP
    REAL(DP), PARAMETER :: V_u_1 = 0.1_DP * SpeedOfLight
    REAL(DP), PARAMETER :: V_u_2 = 0.0_DP * SpeedOfLight
    REAL(DP), PARAMETER :: V_u_3 = 0.0_DP * SpeedOfLight

    INTEGER  :: iE, iX1, iX2, iX3, iS, iNodeE, iNodeX, iNodeZ
    REAL(DP) :: kT, Mnu, E, f_E

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = D_0
        uAF(iNodeX,iX1,iX2,iX3,iAF_T ) = T_0
        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = Y_0

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

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_u_1
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_u_2
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_u_3

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

    ! --- Radiation Fields ---

    DO iS  = 1       , nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        kT = BoltzmannConstant * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

        Mnu = uAF  (iNodeX,iX1,iX2,iX3,iAF_Me) &
              + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
              - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn)

        IF( iS == iNuE_Bar )THEN

          Mnu = - Mnu

        END IF

        E = NodeCoordinate( MeshE, iE, iNodeE )

        f_E = MAX( One / ( EXP( ( E - Mnu ) / kT ) + One ), 1.0d-99 )

        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) = f_E
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) = Zero
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
        uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_N ,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(iNodeZ,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V1),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V2),        &
                 uPF(iNodeX   ,iX1,iX2,iX3,iPF_V3),        &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_11),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_22),  &
                 uGF(iNodeX   ,iX1,iX2,iX3,iGF_Gm_dd_33) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_EquilibriumAdvection


END MODULE InitializationModule_Neutrinos
