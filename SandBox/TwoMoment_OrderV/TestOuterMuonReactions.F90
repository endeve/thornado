PROGRAM TestOuterMuonReactions

  ! Standalone, fixed-thermodynamics, fixed-radiation CC response test.
  ! This is not an evolution solver and does not modify a saved state.

  USE KindModule, ONLY: DP, Zero, One
  USE UnitsModule, ONLY: &
    Gram, Centimeter, Kelvin, MeV, Second, BoltzmannConstant, AtomicMassUnit
  USE OpacityModule_TABLE, ONLY: &
    EmAb_Muon_Method, nE_Muon_2D, &
    MuonCC_QuadratureNodes, MuonCC_QuadratureWeights
  USE EquationOfStateModule_TABLE, ONLY: &
    InitializeEquationOfState_TABLE, FinalizeEquationOfState_TABLE, &
    ComputeMuonChemicalPotential_TABLE, &
    ComputeMuonNeutrinoChemicalPotential_TABLE, &
    ComputeMuonFractionFromDegeneracy_TABLE
  USE RadiationFieldsModule, ONLY: LeptonNumber
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    nOp, iOp_EmAb_Muon, ComputeEquilibriumDistributions_Point, &
    ComputeNeutrinoOpacities_EC
  USE wlSemiLeptonicOpacityModule2D, ONLY: gauleg
  USE MuonSourceArithmeticCheck, ONLY: MuonMaterialNumberSource_Point

  IMPLICIT NONE

  CHARACTER(1024) :: TableName, InputName, OutputName, SpectrumName
  INTEGER :: iInput, iOutput, iSpectrum, nE, nX, nOrders
  INTEGER :: iE, iX, iS, iOrder, iCase, Error, Coordinate, iReference
  INTEGER, ALLOCATABLE :: Order(:), Point(:), Sample(:)
  REAL(DP), ALLOCATABLE :: Energy(:), Weight(:), Radius(:), Rho(:), Temp(:)
  REAL(DP), ALLOCATABLE :: D(:), T(:), Ye(:), Ym(:), Mu(:), MuNu(:)
  REAL(DP), ALLOCATABLE :: Saved(:,:,:), Equilibrium(:,:,:), Chi(:,:,:)
  REAL(DP), ALLOCATABLE :: Distribution(:,:), Transfer(:,:,:)
  LOGICAL, ALLOCATABLE :: Mask(:,:)
  REAL(DP) :: Value, Derivative, Source, AbsoluteSource, Theta, X, Vacancy
  REAL(DP) :: Raw(2), NumberPerBaryon

  IF( COMMAND_ARGUMENT_COUNT() /= 4 ) &
    ERROR STOP 'usage: TestOuterMuonReactions EOS INPUT SUMMARY SPECTRUM'
  CALL GET_COMMAND_ARGUMENT( 1, TableName )
  CALL GET_COMMAND_ARGUMENT( 2, InputName )
  CALL GET_COMMAND_ARGUMENT( 3, OutputName )
  CALL GET_COMMAND_ARGUMENT( 4, SpectrumName )

  OPEN( NEWUNIT = iInput, FILE = TRIM(InputName), STATUS = 'OLD' )
  READ(iInput,*) nE, nX, nOrders
  ALLOCATE( Order(nOrders), Energy(nE), Weight(nE) )
  READ(iInput,*) Order
  DO iE = 1, nE
    READ(iInput,*) Energy(iE), Weight(iE)
  END DO
  ALLOCATE( Point(nX), Sample(nX), Radius(nX), Rho(nX), Temp(nX), &
            D(nX), T(nX), Ye(nX), Ym(nX), Mu(nX), MuNu(nX) )
  ALLOCATE( Saved(nE,3:4,nX), Equilibrium(nE,3:4,nX), Chi(nE,3:4,nX) )
  ALLOCATE( Distribution(nE,2), Transfer(nE,2,3), Mask(nOp,nX) )

  CALL InitializeEquationOfState_TABLE &
         ( EquationOfStateTableName_Option = TRIM(TableName), &
           UseChemicalPotentialShift_Option = .FALSE., Verbose_Option = .FALSE. )
  LeptonNumber = [ One, -One, One, -One, One, -One ]
  EmAb_Muon_Method = 2
  Transfer = Zero
  Mask = .FALSE.
  Mask(iOp_EmAb_Muon,:) = .TRUE.

  DO iX = 1, nX
    READ(iInput,*) Point(iX), Sample(iX), Radius(iX), Rho(iX), Temp(iX), &
                   Ye(iX), Coordinate, Value
    D(iX) = Rho(iX) * Gram / Centimeter**3
    T(iX) = Temp(iX) * Kelvin
    Ym(iX) = Zero
    IF( Coordinate == 1 )THEN
      Theta = BoltzmannConstant * T(iX) / MeV
      CALL ComputeMuonFractionFromDegeneracy_TABLE &
             ( D(iX), T(iX), Ye(iX), Value / Theta, Ym(iX), Derivative, Error )
      IF( Error /= 0 ) ERROR STOP 'Muon EOS forward recovery failed in response test.'
    ELSE
      Ym(iX) = Value
    END IF
    CALL ComputeMuonChemicalPotential_TABLE &
           ( D(iX), T(iX), Ye(iX), Ym(iX), Mu(iX) )
    CALL ComputeMuonNeutrinoChemicalPotential_TABLE &
           ( D(iX), T(iX), Ye(iX), Ym(iX), MuNu(iX), &
             ChargedCurrentActive_Option = .TRUE. )
    DO iE = 1, nE
      READ(iInput,*) Saved(iE,3,iX), Saved(iE,4,iX)
      DO iS = 3, 4
        CALL ComputeEquilibriumDistributions_Point &
               ( Energy(iE) * MeV, D(iX), T(iX), Ye(iX), Ym(iX), &
                 Equilibrium(iE,iS,iX), iS, &
                 MuonChargedCurrentMask_Option = .TRUE. )
      END DO
    END DO
    Mu(iX) = Mu(iX) / MeV
    MuNu(iX) = MuNu(iX) / MeV
  END DO
  CLOSE(iInput)

  OPEN( NEWUNIT = iOutput, FILE = TRIM(OutputName), STATUS = 'REPLACE' )
  OPEN( NEWUNIT = iSpectrum, FILE = TRIM(SpectrumName), STATUS = 'REPLACE' )
  WRITE(iOutput,'(A)') &
    '# nQ point sample scenario r_km rho_gcm3 T_K Ye Ym mu_MeV munu_MeV' // &
    ' source_Ym_s_inv absolute_source_Ym_s_inv'
  WRITE(iOutput,'(A)') '# scenarios: 1 saved final spectrum; 2 equal 1e-40; 3 empty'
  WRITE(iOutput,'(A)') '# scenario 4: represented CC equilibrium at mu_mu=0'
  WRITE(iOutput,'(A,ES26.16E3)') '# atomic_mass_g ', AtomicMassUnit / Gram
  WRITE(iOutput,'(A,ES26.16E3)') &
    '# boltzmann_MeV_K ', BoltzmannConstant * Kelvin / MeV
  WRITE(iSpectrum,'(A)') &
    '# nQ point sample iE E_MeV weight_cm3_inv Ym mu_MeV munu_MeV' // &
    ' chi_nu_s_inv chi_bar_s_inv f0_nu f0_bar raw_nu_s_inv raw_bar_s_inv'

  DO iOrder = 1, nOrders
    nE_Muon_2D = Order(iOrder)
    ALLOCATE( MuonCC_QuadratureNodes(nE_Muon_2D), &
              MuonCC_QuadratureWeights(nE_Muon_2D) )
    CALL gauleg( Zero, One, MuonCC_QuadratureNodes, &
                 MuonCC_QuadratureWeights, nE_Muon_2D )
    CALL ComputeNeutrinoOpacities_EC &
           ( 1, nE, 3, 4, 1, nX, Energy * MeV, D, T, Ye, Ym, &
             Equilibrium, Chi, Mask )
    ! Production opacities are inverse code length. c=1 in code units.
    Chi = Chi * Second
    DO iX = 1, nX
      ! Use the physical baryon-mass conversion from the production units.
      NumberPerBaryon = ( AtomicMassUnit / Gram ) / Rho(iX)
      DO iCase = 1, 4
        SELECT CASE( iCase )
        CASE( 1 )
          Distribution = Saved(:,:,iX)
        CASE( 2 )
          Distribution = 1.0e-40_DP
        CASE( 3 )
          Distribution = Zero
        CASE( 4 )
          iReference = 0
          DO iS = 1, nX
            IF( Point(iS) == Point(iX) .AND. Mu(iS) == Zero ) iReference = iS
          END DO
          IF( iReference == 0 ) ERROR STOP 'No zero-muon reference state in test.'
          Distribution = Equilibrium(:,:,iReference)
        END SELECT
        CALL MuonMaterialNumberSource_Point &
               ( Weight, Chi(:,:,iX), Equilibrium(:,:,iX), Distribution, &
                 Transfer, .TRUE., [.FALSE.,.FALSE.,.FALSE.], &
                 Source, AbsoluteSource )
        WRITE(iOutput,'(4I7,9ES26.16E3)') &
          Order(iOrder), Point(iX), Sample(iX), iCase, &
          Radius(iX), Rho(iX), Temp(iX), Ye(iX), Ym(iX), Mu(iX), MuNu(iX), &
          Source * NumberPerBaryon, AbsoluteSource * NumberPerBaryon
      END DO
      Theta = BoltzmannConstant * T(iX) / MeV
      DO iE = 1, nE
        DO iS = 1, 2
          X = ( Energy(iE) - REAL(3-2*iS,DP)*MuNu(iX) ) / Theta
          ! Independent stable complement; never form 1-f0 near saturation.
          IF( X >= Zero )THEN
            Vacancy = One / ( One + EXP(-X) )
          ELSE
            Vacancy = EXP(X) / ( One + EXP(X) )
          END IF
          ! Reconstruct the raw absorption coefficient from the production
          ! relaxation coefficient. This is NOT an independent kernel check.
          Raw(iS) = Chi(iE,iS+2,iX) * Vacancy
        END DO
        WRITE(iSpectrum,'(4I7,11ES26.16E3)') &
          Order(iOrder), Point(iX), Sample(iX), iE, Energy(iE), Weight(iE), &
          Ym(iX), Mu(iX), MuNu(iX), Chi(iE,3,iX), Chi(iE,4,iX), &
          Equilibrium(iE,3,iX), Equilibrium(iE,4,iX), Raw
      END DO
    END DO
    DEALLOCATE( MuonCC_QuadratureNodes, MuonCC_QuadratureWeights )
  END DO

  CLOSE(iOutput)
  CLOSE(iSpectrum)
  CALL FinalizeEquationOfState_TABLE

END PROGRAM TestOuterMuonReactions
