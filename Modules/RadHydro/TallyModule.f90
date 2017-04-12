MODULE TallyModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    PlanckConstant, &
    SpeedOfLight
  USE ProgramHeaderModule, ONLY: &
    nE, &
    nX
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    WeightsGX, VolJacX, &
    WeightsG,  VolJac, VolJacE,  &
    uGF, iGF_Phi_N
  USE FluidFieldsModule, ONLY: &
    WeightsF, &
    uCF, iCF_D, iCF_E, iCF_Ne
  USE RadiationFieldsModule, ONLY: &
    WeightsR, &
    nSpecies, &
    uCR, iCR_N

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: TallyGravity
  LOGICAL :: TallyFluid
  LOGICAL :: TallyRadiation

  ! --- Tallied Gravity Quantities ---

  REAL(DP), DIMENSION(0:1) :: GlobalEnergy_Gravity

  ! --- Tallied Fluid Quantites ---

  REAL(DP), DIMENSION(0:1) :: GlobalBaryonMass_Fluid
  REAL(DP), DIMENSION(0:1) :: GlobalEnergy_Fluid
  REAL(DP), DIMENSION(0:1) :: GlobalElectronNumber_Fluid

  ! --- Tallied Radiation Quantities ---

  REAL(DP), DIMENSION(0:1) :: GlobalNumber_Radiation
  REAL(DP), DIMENSION(0:1) :: GlobalEnergy_Radiation

  PUBLIC :: InitializeGlobalTally
  PUBLIC :: ComputeGlobalTally

CONTAINS


  SUBROUTINE InitializeGlobalTally &
               ( TallyGravity_Option, TallyFluid_Option, &
                 TallyRadiation_Option )

    LOGICAL, INTENT(in), OPTIONAL :: TallyGravity_Option
    LOGICAL, INTENT(in), OPTIONAL :: TallyFluid_Option
    LOGICAL, INTENT(in), OPTIONAL :: TallyRadiation_Option

    TallyGravity = .FALSE.
    IF( PRESENT( TallyGravity_Option ) ) &
      TallyGravity = TallyGravity_Option

    TallyFluid = .FALSE.
    IF( PRESENT( TallyFluid_Option ) ) &
      TallyFluid = TallyFluid_Option

    TallyRadiation = .FALSE.
    IF( PRESENT( TallyRadiation_Option ) ) &
      TallyRadiation = TallyRadiation_Option

    IF( TallyGravity )THEN

      CALL ComputeGlobalTally_Gravity &
             ( iState_Option = 0 )

    END IF

    IF( TallyFluid )THEN

      CALL ComputeGlobalTally_Fluid &
             ( iState_Option = 0 )

    END IF

    IF( TallyRadiation )THEN

      CALL ComputeGlobalTally_Radiation &
             ( iState_Option = 0 )

    END IF

  END SUBROUTINE InitializeGlobalTally


  SUBROUTINE ComputeGlobalTally

    IF( TallyGravity )THEN

      CALL ComputeGlobalTally_Gravity
      CALL DisplayGlobalTally_Gravity

    END IF

    IF( TallyFluid )THEN

      CALL ComputeGlobalTally_Fluid
      CALL DisplayGlobalTally_Fluid

    END IF

    IF( TallyRadiation )THEN

      CALL ComputeGlobalTally_Radiation
      CALL DisplayGlobalTally_Radiation

    END IF

  END SUBROUTINE ComputeGlobalTally


  SUBROUTINE ComputeGlobalTally_Gravity( iState_Option )

    INTEGER, INTENT(in), OPTIONAL :: iState_Option

    INTEGER :: iS, iX1, iX2, iX3

    iS = 1
    IF( PRESENT( iState_Option ) ) &
      iS = iState_Option

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(1:nX(1)), &
        dX2 => MeshX(2) % Width(1:nX(2)), &
        dX3 => MeshX(3) % Width(1:nX(3)) )

    GlobalEnergy_Gravity(iS) = 0.0_DP
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          ! --- Newtonian Potential Energy ---

          GlobalEnergy_Gravity(iS) &
            = GlobalEnergy_Gravity(iS) &
                + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                    * SUM( WeightsGX(:) * uCF(:,iX1,iX2,iX3,iCF_D) &
                             * uGF(:,iX1,iX2,iX3,iGF_Phi_N) &
                             * VolJacX(:,iX1,iX2,iX3) )

        END DO
      END DO
    END DO

    GlobalEnergy_Gravity(iS) &
      = 0.5_DP * GlobalEnergy_Gravity(iS)

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeGlobalTally_Gravity


  SUBROUTINE DisplayGlobalTally_Gravity

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'INFO: Gravitational Tally'
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Global Potential Energy = ', &
      GlobalEnergy_Gravity(1)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Change = ', &
      GlobalEnergy_Gravity(1) &
        - GlobalEnergy_Gravity(0)
    WRITE(*,*)

  END SUBROUTINE DisplayGlobalTally_Gravity


  SUBROUTINE ComputeGlobalTally_Fluid( iState_Option )

    INTEGER, INTENT(in), OPTIONAL :: iState_Option
 
    INTEGER :: iS, iX1, iX2, iX3

    iS = 1
    IF( PRESENT( iState_Option ) ) &
      iS = iState_Option

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(1:nX(1)), &
        dX2 => MeshX(2) % Width(1:nX(2)), &
        dX3 => MeshX(3) % Width(1:nX(3)) )

    GlobalBaryonMass_Fluid    (iS) = 0.0_DP
    GlobalEnergy_Fluid        (iS) = 0.0_DP
    GlobalElectronNumber_Fluid(iS) = 0.0_DP
    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          GlobalBaryonMass_Fluid(iS) &
            = GlobalBaryonMass_Fluid(iS) &
                + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                    * SUM( WeightsF(:) * uCF(:,iX1,iX2,iX3,iCF_D) &
                             * VolJacX(:,iX1,iX2,iX3) )

          GlobalEnergy_Fluid(iS) &
            = GlobalEnergy_Fluid(iS) &
                + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                    * SUM( WeightsF(:) * uCF(:,iX1,iX2,iX3,iCF_E) &
                             * VolJacX(:,iX1,iX2,iX3) )

          GlobalElectronNumber_Fluid(iS) &
            = GlobalElectronNumber_Fluid(iS) &
                + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                    * SUM( WeightsF(:) * uCF(:,iX1,iX2,iX3,iCF_Ne) &
                             * VolJacX(:,iX1,iX2,iX3) )

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeGlobalTally_Fluid


  SUBROUTINE DisplayGlobalTally_Fluid

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'INFO: Fluid Tally'
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Global Baryon Mass = ', &
      GlobalBaryonMass_Fluid(1)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Change = ', &
      GlobalBaryonMass_Fluid(1) &
        - GlobalBaryonMass_Fluid(0)
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Global Energy = ', &
      GlobalEnergy_Fluid(1)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Change = ', &
      GlobalEnergy_Fluid(1) &
        - GlobalEnergy_Fluid(0)
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Global Electron Number = ', &
      GlobalElectronNumber_Fluid(1)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Change = ', &
      GlobalElectronNumber_Fluid(1) &
        - GlobalElectronNumber_Fluid(0)
    WRITE(*,*)

  END SUBROUTINE DisplayGlobalTally_Fluid


  SUBROUTINE ComputeGlobalTally_Radiation( iState_Option )

    INTEGER, INTENT(in), OPTIONAL :: iState_Option

    INTEGER :: iState
    INTEGER :: iS, iX1, iX2, iX3, iE

    iState = 1
    IF( PRESENT( iState_Option ) ) &
      iState = iState_Option

    ASSOCIATE &
      ( dE  => MeshE    % Width(1:nE   ), &
        dX1 => MeshX(1) % Width(1:nX(1)), &
        dX2 => MeshX(2) % Width(1:nX(2)), &
        dX3 => MeshX(3) % Width(1:nX(3)) )

    ASSOCIATE &
      ( hc3 => ( PlanckConstant * SpeedOfLight )**3 )

    GlobalNumber_Radiation(iState) = 0.0_DP
    GlobalEnergy_Radiation(iState) = 0.0_DP
    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3) 
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              GlobalNumber_Radiation(iState) &
                = GlobalNumber_Radiation(iState) &
                    + dE(iE) * dX1(iX1) * dX2(iX2) * dX3(iX3) &
                        * SUM( WeightsR(:) * uCR(:,iE,iX1,iX2,iX3,iCR_N,iS) &
                                 * VolJac(:,iE,iX1,iX2,iX3) ) / hc3

              GlobalEnergy_Radiation(iState) &
                = GlobalEnergy_Radiation(iState) &
                    + dE(iE) * dX1(iX1) * dX2(iX2) * dX3(iX3) &
                        * SUM( WeightsR(:) * uCR(:,iE,iX1,iX2,iX3,iCR_N,iS) &
                                 * VolJacE(:,iE,iX1,iX2,iX3) ) / hc3

            END DO
          END DO
        END DO
      END DO

    END DO

    END ASSOCIATE ! hc3

    END ASSOCIATE ! dE, etc.

  END SUBROUTINE ComputeGlobalTally_Radiation


  SUBROUTINE DisplayGlobalTally_Radiation

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'INFO: Radiation Tally'
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Global Number = ', &
      GlobalNumber_Radiation(1)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Change = ', &
      GlobalNumber_Radiation(1) &
        - GlobalNumber_Radiation(0)
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Global Energy = ', &
      GlobalEnergy_Radiation(1)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Change = ', &
      GlobalEnergy_Radiation(1) &
        - GlobalEnergy_Radiation(0)
    WRITE(*,*)

  END SUBROUTINE DisplayGlobalTally_Radiation


END MODULE TallyModule
