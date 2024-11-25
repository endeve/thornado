MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, TwoPi
  USE UnitsModule, ONLY: &
    Gram, Centimeter, &
    Kilometer, Kelvin, &
    BoltzmannConstant, MeV
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0, &
    nDOF, nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &
    NodeNumbersX, &
    NodeNumberTable
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
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule_TABLE, ONLY: &
    ApplyEquationOfState_TABLE, &
    ComputeThermodynamicStates_Primitive_TABLE
  USE TwoMoment_UtilitiesModule, ONLY: &
    ComputeConserved_TwoMoment

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializeFields_Relaxation

CONTAINS


  SUBROUTINE InitializeFields_Relaxation( D_0, T_0, Y_0 )

    REAL(DP), INTENT(in) :: D_0, T_0, Y_0

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'InitializeFields_Relaxation'

    CALL InitializeFluidFields_Relaxation( D_0, T_0, Y_0 )

    CALL InitializeRadiationFields_Relaxation

  END SUBROUTINE InitializeFields_Relaxation


  SUBROUTINE InitializeFluidFields_Relaxation( D_0, T_0, Y_0 )

    REAL(DP), INTENT(in) :: D_0, T_0, Y_0

    INTEGER  :: iX1, iX2, iX3, iNodeX
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)
        iNodeX3 = NodeNumberTableX(3,iNodeX)

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) = D_0

        uAF(iNodeX,iX1,iX2,iX3,iAF_T ) = T_0

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = Y_0

      END DO

      CALL ComputeThermodynamicStates_Primitive_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uPF(:,iX1,iX2,iX3,iPF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne) )

      uPF(:,iX1,iX2,iX3,iPF_V1) = Zero
      uPF(:,iX1,iX2,iX3,iPF_V2) = Zero
      uPF(:,iX1,iX2,iX3,iPF_V3) = Zero

      CALL ApplyEquationOfState_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
               uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
               uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
               uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

      ! --- Conserved ---

      uCF(:,iX1,iX2,iX3,iCF_D ) = uPF(:,iX1,iX2,iX3,iPF_D)
      uCF(:,iX1,iX2,iX3,iCF_S1) = Zero
      uCF(:,iX1,iX2,iX3,iCF_S2) = Zero
      uCF(:,iX1,iX2,iX3,iCF_S3) = Zero
      uCF(:,iX1,iX2,iX3,iCF_E ) = uPF(:,iX1,iX2,iX3,iPF_E)
      uCF(:,iX1,iX2,iX3,iCF_Ne) = uPF(:,iX1,iX2,iX3,iPF_Ne)

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFluidFields_Relaxation


  SUBROUTINE InitializeRadiationFields_Relaxation

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNode, iNodeE
    REAL(DP) :: kT(nDOF)
    REAL(DP) :: Mnu(nDOF), E
    REAL(DP) :: Gm_dd_11(nDOF)
    REAL(DP) :: Gm_dd_22(nDOF)
    REAL(DP) :: Gm_dd_33(nDOF)

    DO iS = 1, nSpecies
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Gm_dd_11 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_11)

      Gm_dd_22 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_22)

      Gm_dd_33 &
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_33)

      kT = BoltzmannConstant &
             * uAF(NodeNumbersX,iX1,iX2,iX3,iAF_T)

      IF( iS .EQ. iNuE )THEN

        Mnu = + ( uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn) )

      ELSEIF( iS .EQ. iNuE_Bar )THEN

        Mnu = - ( uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn) )

      END IF

      DO iE = iE_B0, iE_E0

        DO iNode = 1, nDOF

          iNodeE = NodeNumberTable(1,iNode)

          E = NodeCoordinate( MeshE, iE, iNodeE )

          uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
            = MAX( 0.99_DP * EXP( - ( E - Two*kT(iNode) )**2 &
                                    / ( Two*(1.0d1*MeV)**2 ) ), 1.0d-99 )

!          uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
!            = MAX( One / ( EXP( (E-Mnu(iNode))/kT(iNode) ) + One ), 1.0d-99 )

          uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) = Zero
          uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero
          uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

        END DO

        CALL ComputeConserved_TwoMoment &
               ( uPR(:,iE,iX1,iX2,iX3,iPR_D, iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_N, iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      END DO

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeRadiationFields_Relaxation


END MODULE InitializationModule
