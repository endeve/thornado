MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Three, Pi, TwoPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFZ, nDOFX, nDOFE, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    MeshE, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Alpha, iGF_Beta_1, iGF_Beta_2, iGF_Beta_3
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  USE EquationOfStateModule_IDEAL, ONLY: &
    ComputePressureFromPrimitive_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields


CONTAINS


  SUBROUTINE InitializeFields( V_0, LengthScale )

    REAL(DP), INTENT(in) :: V_0(3)
    REAL(DP), INTENT(in) :: LengthScale

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') '', 'INFO: ', TRIM( ProgramName )

    SELECT CASE( TRIM( ProgramName ) )

      CASE( 'TransparentShock' )

        CALL InitializeFields_TransparentShock( V_0, LengthScale )
    
      CASE( 'StreamingDopplerShift' )

        CALL InitializeFields_StreamingDopplerShift &
               ( V_0 )


    END SELECT

  END SUBROUTINE InitializeFields


  SUBROUTINE InitializeFields_TransparentShock( V_0, ShockWidth )

    REAL(DP), INTENT(in) :: V_0(3)
    REAL(DP), INTENT(in) :: ShockWidth

    REAL(DP), PARAMETER :: X_Shock = 1.0d0
    REAL(DP), PARAMETER :: Mu_0    = 0.9d0

    INTEGER  :: iNodeX, iX1, iX2, iX3, iNodeE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    REAL(DP) :: X1, X2, X3, Pressure, E
    REAL(DP) :: V_u_1, V_u_2, V_u_3
    REAL(DP) :: Gm_dd_11, Gm_dd_22, Gm_dd_33
    REAL(DP) :: V2, W

    WRITE(*,*)
    WRITE(*,'(A6,A13,3ES9.2E2)') '', 'V_0 = ', V_0
    WRITE(*,'(A6,A13,1ES9.2E2)') '', 'ShockWidth = ', ShockWidth
    WRITE(*,*)

    ! --- Fluid Fields ---

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
        X3 = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
          = One
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
          = Half * V_0(1) * ( One + TANH( (X1-X_Shock)/ShockWidth ) )
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
          = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
          = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) &
          = 1.0d-1
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) &
          = 0.0d-0

        CALL ComputePressureFromPrimitive_IDEAL &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne), &
                 Pressure )

        CALL ComputeConserved_Euler_Relativistic &
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
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 Pressure )

      END DO

    END DO
    END DO
    END DO

    ! --- Radiation Fields ---

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B1(4), iZ_E1(4)
    DO iZ3 = iZ_B1(3), iZ_E1(3)
    DO iZ2 = iZ_B1(2), iZ_E1(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1
        iNodeE = MOD( (iNodeZ-1)        , nDOFE ) + 1

        V_u_1 = uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V1)
        V_u_2 = uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V2)
        V_u_3 = uPF(iNodeX,iZ2,iZ3,iZ4,iPF_V3)

        Gm_dd_11 = uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)
        Gm_dd_22 = uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)
        Gm_dd_33 = uGF(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)

        V2 =   Gm_dd_11 * V_u_1 * V_u_1 &
             + Gm_dd_22 * V_u_2 * V_u_2 &
             + Gm_dd_33 * V_u_3 * V_u_3

        W  = One / SQRT( One - V2 )

        E = NodeCoordinate( MeshE, iZ1, iNodeE )

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 0.5_DP * ( One - Mu_0 ) / ( EXP( E / Three - Three ) + One )
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = 0.5_DP * ( One + Mu_0 ) * W * uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D,iS)
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = 0.0_DP
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = 0.0_DP

        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 V_u_1, V_u_2, V_u_3, &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33, &
                 Zero, Zero, Zero, One, Zero, Zero, Zero )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_TransparentShock

  SUBROUTINE InitializeFields_StreamingDopplerShift &
    ( V_0 )

    REAL(DP),      INTENT(in) :: V_0(3)

    REAL(DP), PARAMETER :: X_0 = 2.0_DP
    REAL(DP), PARAMETER :: X_1 = 3.5_DP
    REAL(DP), PARAMETER :: X_2 = 6.5_DP
    REAL(DP), PARAMETER :: X_3 = 8.0_DP
    REAL(DP), PARAMETER :: L_X = 6.0_DP

    INTEGER  :: iNodeX, iX1, iX2, iX3
    INTEGER  :: iNodeZ, iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3
    REAL(DP) :: X1, X2, X3, uAF(nDOFX)

    ! --- Fluid Fields ---

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = 1.0_DP


          IF( X1 .LT. X_0 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = 0.0_DP
          ELSEIF( X1 .GE. X_0 .AND. X1 .LT. X_1 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSEIF( X1 .GE. X_1 .AND. X1 .LT. X_2 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = V_0(1)
          ELSEIF( X1 .GE. X_2 .AND. X1 .LT. X_3 )THEN
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = V_0(1) * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
          ELSE
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = 0.0_DP
          END IF
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_0(2)
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_0(3)


        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = 0.1_DP
        uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = 0.0_DP

        CALL ComputePressureFromPrimitive_IDEAL &
                 ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                   uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &          
                   uPF(iNodeX,iX1,iX2,iX3,iPF_Ne ), &          
                   uAF(iNodeX) )         

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1 ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2 ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3 ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ), &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ), &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNodeX) )  

      END DO

    END DO
    END DO
    END DO
    
    ! --- Radiation Fields ---
    
    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeZ = 1, nDOFZ

        iNodeX = MOD( (iNodeZ-1) / nDOFE, nDOFX ) + 1

        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS) &
          = 1.0d-40
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS) &
          = Zero
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS) &
          = Zero
        uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS) &
          = Zero
        
        CALL ComputeConserved_TwoMoment &
               ( uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_D ,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I1,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I2,iS), &
                 uPR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iPR_I3,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N ,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS), &
                 uCR(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS), &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V1),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V2),        &
                 uPF(iNodeX    ,iZ2,iZ3,iZ4,iPF_V3),        &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Gm_dd_33), &
                 0.0_DP, 0.0_DP, 0.0_DP,                   & 
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Alpha),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Beta_1),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Beta_2),  &
                 uGF(iNodeX    ,iZ2,iZ3,iZ4,iGF_Beta_3) )
      
         END DO

    END DO
    END DO
    END DO
    END DO
    END DO


  END SUBROUTINE InitializeFields_StreamingDopplerShift

END MODULE InitializationModule
