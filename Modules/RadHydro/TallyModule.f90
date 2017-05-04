MODULE TallyModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    PlanckConstant, &
    SpeedOfLight, &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nE, nX, nDOFX
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    WeightsGX, VolJacX, &
    WeightsG,  VolJac, VolJacE,  &
    uGF, iGF_Phi_N
  USE FluidFieldsModule, ONLY: &
    WeightsF, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  USE RadiationFieldsModule, ONLY: &
    WeightsR, &
    nSpecies, &
    uCR, iCR_N

  IMPLICIT NONE
  PRIVATE

  CHARACTER(80) :: TallyFileName
  LOGICAL       :: TallyGravity
  LOGICAL       :: TallyFluid
  LOGICAL       :: TallyRadiation

  ! --- Tallied Gravity Quantities ---

  REAL(DP), DIMENSION(0:1) :: GlobalEnergy_Gravity

  ! --- Tallied Fluid Quantites ---

  REAL(DP)                 :: MaximumMassDensity
  REAL(DP), DIMENSION(0:1) :: GlobalBaryonMass_Fluid
  REAL(DP), DIMENSION(0:1) :: GlobalEnergy_Fluid
  REAL(DP), DIMENSION(0:1) :: GlobalInternalEnergy_Fluid
  REAL(DP), DIMENSION(0:1) :: GlobalKineticEnergy_Fluid
  REAL(DP), DIMENSION(0:1) :: GlobalElectronNumber_Fluid

  ! --- Tallied Radiation Quantities ---

  REAL(DP), DIMENSION(0:1) :: GlobalNumber_Radiation
  REAL(DP), DIMENSION(0:1) :: GlobalEnergy_Radiation

  PUBLIC :: InitializeGlobalTally
  PUBLIC :: ComputeGlobalTally

CONTAINS


  SUBROUTINE InitializeGlobalTally &
               ( Time, TallyGravity_Option, TallyFluid_Option, &
                 TallyRadiation_Option )

    REAL(DP), INTENT(in), OPTIONAL :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: TallyGravity_Option
    LOGICAL,  INTENT(in), OPTIONAL :: TallyFluid_Option
    LOGICAL,  INTENT(in), OPTIONAL :: TallyRadiation_Option

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
      CALL ComputeGlobalTally_Gravity &
             ( iState_Option = 1 )

    END IF

    IF( TallyFluid )THEN

      CALL ComputeGlobalTally_Fluid &
             ( iState_Option = 0 )
      CALL ComputeGlobalTally_Fluid &
             ( iState_Option = 1 )

    END IF

    IF( TallyRadiation )THEN

      CALL ComputeGlobalTally_Radiation &
             ( iState_Option = 0 )
      CALL ComputeGlobalTally_Radiation &
             ( iState_Option = 1 )

    END IF

    TallyFileName &
      = '../Output/' // TRIM( ProgramName ) // '_GlobalTally.dat'

    CALL WriteGlobalTally( Time )

  END SUBROUTINE InitializeGlobalTally


  SUBROUTINE ComputeGlobalTally( Time )

    REAL(DP), INTENT(in) :: Time

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

    CALL WriteGlobalTally( Time, Append_Option = .TRUE. )

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

    ASSOCIATE( U => UnitsDisplay )

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'INFO: Gravitational Tally'
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Global Potential Energy = ', &
      GlobalEnergy_Gravity(1) / U % EnergyGlobalUnit, &
      '', U % EnergyGlobalLabel
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Change = ', &
      ( GlobalEnergy_Gravity(1) &
        - GlobalEnergy_Gravity(0) ) / U % EnergyGlobalUnit, &
      '', U % EnergyGlobalLabel
    WRITE(*,*)

    END ASSOCIATE ! U

  END SUBROUTINE DisplayGlobalTally_Gravity


  SUBROUTINE ComputeGlobalTally_Fluid( iState_Option )

    INTEGER, INTENT(in), OPTIONAL :: iState_Option
 
    INTEGER :: iS, iX1, iX2, iX3
    REAL(DP), DIMENSION(nDOFX) :: Eint, Ekin

    iS = 1
    IF( PRESENT( iState_Option ) ) &
      iS = iState_Option

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width(1:nX(1)), &
        dX2 => MeshX(2) % Width(1:nX(2)), &
        dX3 => MeshX(3) % Width(1:nX(3)) )

    MaximumMassDensity             = 0.0_DP
    GlobalBaryonMass_Fluid    (iS) = 0.0_DP
    GlobalEnergy_Fluid        (iS) = 0.0_DP
    GlobalInternalEnergy_Fluid(iS) = 0.0_DP
    GlobalKineticEnergy_Fluid (iS) = 0.0_DP
    GlobalElectronNumber_Fluid(iS) = 0.0_DP

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          MaximumMassDensity &
            = MAX( MaximumMassDensity, &
                   MAXVAL( uCF(:,iX1,iX2,iX3,iCF_D) ) )

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

          Ekin(:) = 0.5_DP * ( uCF(:,iX1,iX2,iX3,iCF_S1)**2 &
                               + uCF(:,iX1,iX2,iX3,iCF_S2)**2 &
                               + uCF(:,iX1,iX2,iX3,iCF_S3)**2 ) &
                    / uCF(:,iX1,iX2,iX3,iCF_D)

          Eint(:) = uCF(:,iX1,iX2,iX3,iCF_E) - Ekin(:)

          GlobalInternalEnergy_Fluid(iS) &
            = GlobalInternalEnergy_Fluid(iS) &
                + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                    * SUM( WeightsF(:) * Eint(:) * VolJacX(:,iX1,iX2,iX3) )

          GlobalKineticEnergy_Fluid(iS) &
            = GlobalKineticEnergy_Fluid(iS) &
                + dX1(iX1) * dX2(iX2) * dX3(iX3) &
                    * SUM( WeightsF(:) * Ekin(:) * VolJacX(:,iX1,iX2,iX3) )

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

    ASSOCIATE( U => UnitsDisplay )

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'INFO: Fluid Tally'
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Maximum Mass Density = ', &
      MaximumMassDensity / U % MassDensityUnit, &
      '', U % MassDensityLabel
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Global Baryon Mass = ', &
      GlobalBaryonMass_Fluid(1) &
      / U % MassUnit, &
      '', U % MassLabel
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Change = ', &
      ( GlobalBaryonMass_Fluid(1) &
        - GlobalBaryonMass_Fluid(0) ) &
      / U % MassUnit, &
      '', U % MassLabel
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Global Energy = ', &
      GlobalEnergy_Fluid(1) &
      / U % EnergyGlobalUnit, &
      '', U % EnergyGlobalLabel
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', '  Internal Energy = ', &
      GlobalInternalEnergy_Fluid(1) &
      / U % EnergyGlobalUnit, &
      '', U % EnergyGlobalLabel
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', '  Kinetic Energy = ', &
      GlobalKineticEnergy_Fluid(1) &
      / U % EnergyGlobalUnit, &
      '', U % EnergyGlobalLabel
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Change = ', &
      ( GlobalEnergy_Fluid(1) &
        - GlobalEnergy_Fluid(0) ) &
      / U % EnergyGlobalUnit, &
      '', U % EnergyGlobalLabel
    WRITE(*,*)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Global Electron Number = ', &
      GlobalElectronNumber_Fluid(1)
    WRITE(*,'(A8,A26,ES18.10E3)') &
      '', 'Change = ', &
      GlobalElectronNumber_Fluid(1) &
        - GlobalElectronNumber_Fluid(0)
    WRITE(*,*)

    END ASSOCIATE ! U

  END SUBROUTINE DisplayGlobalTally_Fluid


  SUBROUTINE ComputeGlobalTally_Radiation( iState_Option )

    INTEGER, INTENT(in), OPTIONAL :: iState_Option

    INTEGER  :: iState
    INTEGER  :: iS, iX1, iX2, iX3, iE
    REAL(DP) :: Scale

    iState = 1
    IF( PRESENT( iState_Option ) ) &
      iState = iState_Option

    ASSOCIATE &
      ( dE  => MeshE    % Width(1:nE   ), &
        dX1 => MeshX(1) % Width(1:nX(1)), &
        dX2 => MeshX(2) % Width(1:nX(2)), &
        dX3 => MeshX(3) % Width(1:nX(3)) )

    Scale = 1.0_DP
    IF( UnitsDisplay % Active ) &
      Scale = ( PlanckConstant * SpeedOfLight )**3

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
                                 * VolJac(:,iE,iX1,iX2,iX3) ) / Scale

              GlobalEnergy_Radiation(iState) &
                = GlobalEnergy_Radiation(iState) &
                    + dE(iE) * dX1(iX1) * dX2(iX2) * dX3(iX3) &
                        * SUM( WeightsR(:) * uCR(:,iE,iX1,iX2,iX3,iCR_N,iS) &
                                 * VolJacE(:,iE,iX1,iX2,iX3) ) / Scale

            END DO
          END DO
        END DO
      END DO

    END DO

    END ASSOCIATE ! dE, etc.

  END SUBROUTINE ComputeGlobalTally_Radiation


  SUBROUTINE DisplayGlobalTally_Radiation

    ASSOCIATE( U => UnitsDisplay )

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
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Global Energy = ', &
      GlobalEnergy_Radiation(1) &
      / U % EnergyGlobalUnit, &
      '', U % EnergyGlobalLabel
    WRITE(*,'(A8,A26,ES18.10E3,A1,A)') &
      '', 'Change = ', &
      ( GlobalEnergy_Radiation(1) &
        - GlobalEnergy_Radiation(0) ) &
      / U % EnergyGlobalUnit, &
      '', U % EnergyGlobalLabel
    WRITE(*,*)

    END ASSOCIATE ! U

  END SUBROUTINE DisplayGlobalTally_Radiation


  SUBROUTINE WriteGlobalTally( Time, Append_Option )

    REAL(DP), INTENT(in)           :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: Append_Option

    CHARACTER(6) :: AccessOption = 'APPEND'
    LOGICAL      :: Append
    INTEGER      :: FileUnit

    Append = .FALSE.
    IF( PRESENT( Append_Option ) ) &
      Append = Append_Option

    IF( .NOT. Append )THEN

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ) )

      WRITE( FileUnit, '(16(A14,x))' ) &
        'Time',          'Max D',      &
        'Total M',       'Change M',   &
        'Total E_F',                   &
        'Total E_F (I)',               &
        'Total E_F (K)', 'Change E_F', &
        'Total N_F',     'Change N_F', &
        'Total E_G',     'Change E_G', &
        'Total N_R',     'Change N_R', &
        'Total E_R',     'Change E_R'

    ELSE

      OPEN( NEWUNIT = FileUnit, FILE = TRIM( TallyFileName ), &
            ACCESS = AccessOption )

    END IF

    ASSOCIATE( U => UnitsDisplay )

    WRITE( FileUnit,                       &
           '(16(ES14.5,x))' )              &
      Time / U % TimeUnit,                 &
      MaximumMassDensity                   &
      / U % MassDensityUnit,               &
      GlobalBaryonMass_Fluid(1)            &
      / U % MassUnit,                      &
      ( GlobalBaryonMass_Fluid(1)          &
        - GlobalBaryonMass_Fluid(0) )      &
      / U % MassUnit,                      &
      GlobalEnergy_Fluid(1)                &
      / U % EnergyGlobalUnit,              &
      GlobalInternalEnergy_Fluid(1)        &
      / U % EnergyGlobalUnit,              &
      GlobalKineticEnergy_Fluid(1)         &
      / U % EnergyGlobalUnit,              &
      ( GlobalEnergy_Fluid(1)              &
        - GlobalEnergy_Fluid(0) )          &
      / U % EnergyGlobalUnit,              &
      GlobalElectronNumber_Fluid(1),       &
      ( GlobalElectronNumber_Fluid(1)      &
        - GlobalElectronNumber_Fluid(0) ), &
      GlobalEnergy_Gravity(1)              &
      / U % EnergyGlobalUnit,              &
      ( GlobalEnergy_Gravity(1)            &
        - GlobalEnergy_Gravity(0) )        &
      / U % EnergyGlobalUnit,              &
      GlobalNumber_Radiation(1),           &
      ( GlobalNumber_Radiation(1)          &
        - GlobalNumber_Radiation(0) ),     &
      GlobalEnergy_Radiation(1)            &
      / U % EnergyGlobalUnit,              &
      ( GlobalEnergy_Radiation(1)          &
        - GlobalEnergy_Radiation(0) )      &
      / U % EnergyGlobalUnit

    END ASSOCIATE ! U

    CLOSE( FileUnit )

  END SUBROUTINE WriteGlobalTally


END MODULE TallyModule
