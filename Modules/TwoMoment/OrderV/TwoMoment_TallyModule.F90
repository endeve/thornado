MODULE TwoMoment_TallyModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, FourPi
  USE UnitsModule, ONLY: &
    UnitsActive, &
    SpeedOfLight, &
    PlanckConstant, &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nDOFE, &
    nDOFX, &
    nDOFZ
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE ReferenceElementModule, ONLY: &
    Weights_q
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModuleE, ONLY: &
    nGE, iGE_Ep2, iGE_Ep3
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, iGF_SqrtGm, iGF_Phi_N
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE RadiationFieldsModule, ONLY: &
    nSpecies, LeptonNumber, &
    nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeTally
  PUBLIC :: ComputeTally
  PUBLIC :: IncrementOffGridTally_TwoMoment
  PUBLIC :: IncrementPositivityLimiterTally_TwoMoment
  PUBLIC :: IncrementOffGridTally_Euler
  PUBLIC :: FinalizeTally

  REAL(DP) :: hc3

  CHARACTER(256) :: NeutrinoLeptonNumber_FileName
  REAL(DP)       :: NeutrinoLeptonNumber_Interior
  REAL(DP)       :: NeutrinoLeptonNumber_Initial
  REAL(DP)       :: NeutrinoLeptonNumber_OffGrid
  REAL(DP)       :: NeutrinoLeptonNumber_Change

  CHARACTER(256) :: NeutrinoEnergy_FileName
  REAL(DP)       :: NeutrinoEnergy_Interior
  REAL(DP)       :: NeutrinoEnergy_Initial
  REAL(DP)       :: NeutrinoEnergy_OffGrid
  REAL(DP)       :: NeutrinoEnergy_Change
  CHARACTER(256) :: NeutrinoEnergy_PL_FileName
  REAL(DP)       :: NeutrinoEnergy_PL

  CHARACTER(256) :: NeutrinoMomentum1_FileName
  REAL(DP)       :: NeutrinoMomentum1_Interior
  REAL(DP)       :: NeutrinoMomentum1_Initial
  REAL(DP)       :: NeutrinoMomentum1_OffGrid
  REAL(DP)       :: NeutrinoMomentum1_Change
  CHARACTER(256) :: NeutrinoMomentum1_PL_FileName
  REAL(DP)       :: NeutrinoMomentum1_PL

  CHARACTER(256) :: NeutrinoMomentum2_FileName
  REAL(DP)       :: NeutrinoMomentum2_Interior
  REAL(DP)       :: NeutrinoMomentum2_Initial
  REAL(DP)       :: NeutrinoMomentum2_OffGrid
  REAL(DP)       :: NeutrinoMomentum2_Change
  CHARACTER(256) :: NeutrinoMomentum2_PL_FileName
  REAL(DP)       :: NeutrinoMomentum2_PL

  CHARACTER(256) :: NeutrinoMomentum3_FileName
  REAL(DP)       :: NeutrinoMomentum3_Interior
  REAL(DP)       :: NeutrinoMomentum3_Initial
  REAL(DP)       :: NeutrinoMomentum3_OffGrid
  REAL(DP)       :: NeutrinoMomentum3_Change
  CHARACTER(256) :: NeutrinoMomentum3_PL_FileName
  REAL(DP)       :: NeutrinoMomentum3_PL

  CHARACTER(256) :: FluidLeptonNumber_FileName
  REAL(DP)       :: FluidLeptonNumber_Interior
  REAL(DP)       :: FluidLeptonNumber_Initial
  REAL(DP)       :: FluidLeptonNumber_OffGrid
  REAL(DP)       :: FluidLeptonNumber_Change

  CHARACTER(256) :: FluidEnergy_FileName
  REAL(DP)       :: FluidEnergy_Interior
  REAL(DP)       :: FluidEnergy_Initial
  REAL(DP)       :: FluidEnergy_OffGrid
  REAL(DP)       :: FluidEnergy_Change

  CHARACTER(256) :: FluidMomentum1_FileName
  REAL(DP)       :: FluidMomentum1_Interior
  REAL(DP)       :: FluidMomentum1_Initial
  REAL(DP)       :: FluidMomentum1_OffGrid
  REAL(DP)       :: FluidMomentum1_Change

  CHARACTER(256) :: FluidMomentum2_FileName
  REAL(DP)       :: FluidMomentum2_Interior
  REAL(DP)       :: FluidMomentum2_Initial
  REAL(DP)       :: FluidMomentum2_OffGrid
  REAL(DP)       :: FluidMomentum2_Change

  CHARACTER(256) :: FluidMomentum3_FileName
  REAL(DP)       :: FluidMomentum3_Interior
  REAL(DP)       :: FluidMomentum3_Initial
  REAL(DP)       :: FluidMomentum3_OffGrid
  REAL(DP)       :: FluidMomentum3_Change

  CHARACTER(256) :: GravitationalEnergy_FileName
  REAL(DP)       :: GravitationalEnergy_Interior
  REAL(DP)       :: GravitationalEnergy_Initial
  REAL(DP)       :: GravitationalEnergy_OffGrid
  REAL(DP)       :: GravitationalEnergy_Change

CONTAINS


  SUBROUTINE InitializeTally

    CHARACTER(256) :: BaseFileName

    BaseFileName = '../Output/' // TRIM( ProgramName )

    ! --- Neutrino Lepton Number ---

    NeutrinoLeptonNumber_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoLeptonNumber.dat'

    CALL WriteTally_Header( NeutrinoLeptonNumber_FileName )

    NeutrinoLeptonNumber_Interior = Zero
    NeutrinoLeptonNumber_Initial  = Zero
    NeutrinoLeptonNumber_OffGrid  = Zero
    NeutrinoLeptonNumber_Change   = Zero

    ! --- Neutrino Energy ---

    NeutrinoEnergy_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoEnergy.dat'

    CALL WriteTally_Header( NeutrinoEnergy_FileName )

    NeutrinoEnergy_PL_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoEnergy_PL.dat'

    CALL WriteLimiterTally_Header( NeutrinoEnergy_PL_FileName )

    NeutrinoEnergy_Interior = Zero
    NeutrinoEnergy_Initial  = Zero
    NeutrinoEnergy_OffGrid  = Zero
    NeutrinoEnergy_Change   = Zero
    NeutrinoEnergy_PL       = Zero

    ! --- Neutrino Momentum 1 ---

    NeutrinoMomentum1_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoMomentum1.dat'

    CALL WriteTally_Header( NeutrinoMomentum1_FileName )

    NeutrinoMomentum1_PL_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoMomentum1_PL.dat'

    CALL WriteLimiterTally_Header( NeutrinoMomentum1_PL_FileName )

    NeutrinoMomentum1_Interior = Zero
    NeutrinoMomentum1_Initial  = Zero
    NeutrinoMomentum1_OffGrid  = Zero
    NeutrinoMomentum1_Change   = Zero
    NeutrinoMomentum1_PL       = Zero

    ! --- Neutrino Momentum 2 ---

    NeutrinoMomentum2_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoMomentum2.dat'

    CALL WriteTally_Header( NeutrinoMomentum2_FileName )

    NeutrinoMomentum2_PL_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoMomentum2_PL.dat'

    CALL WriteLimiterTally_Header( NeutrinoMomentum2_PL_FileName )

    NeutrinoMomentum2_Interior = Zero
    NeutrinoMomentum2_Initial  = Zero
    NeutrinoMomentum2_OffGrid  = Zero
    NeutrinoMomentum2_Change   = Zero
    NeutrinoMomentum2_PL       = Zero

    ! --- Neutrino Momentum 3 ---

    NeutrinoMomentum3_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoMomentum3.dat'

    CALL WriteTally_Header( NeutrinoMomentum3_FileName )

    NeutrinoMomentum3_PL_FileName &
      = TRIM( BaseFileName ) // '_Tally_NeutrinoMomentum3_PL.dat'

    CALL WriteLimiterTally_Header( NeutrinoMomentum3_PL_FileName )

    NeutrinoMomentum3_Interior = Zero
    NeutrinoMomentum3_Initial  = Zero
    NeutrinoMomentum3_OffGrid  = Zero
    NeutrinoMomentum3_Change   = Zero
    NeutrinoMomentum3_PL       = Zero

    ! --- Fluid Lepton Number ---

    FluidLeptonNumber_FileName &
      = TRIM( BaseFileName ) // '_Tally_FluidLeptonNumber.dat'

    CALL WriteTally_Header( FluidLeptonNumber_FileName )

    FluidLeptonNumber_Interior = Zero
    FluidLeptonNumber_Initial  = Zero
    FluidLeptonNumber_OffGrid  = Zero
    FluidLeptonNumber_Change   = Zero

    ! --- Fluid Energy ---

    FluidEnergy_FileName &
      = TRIM( BaseFileName ) // '_Tally_FluidEnergy.dat'

    CALL WriteTally_Header( FluidEnergy_FileName )

    FluidEnergy_Interior = Zero
    FluidEnergy_Initial  = Zero
    FluidEnergy_OffGrid  = Zero
    FluidEnergy_Change   = Zero

    ! --- Fluid Momentum 1 ---

    FluidMomentum1_FileName &
      = TRIM( BaseFileName ) // '_Tally_FluidMomentum1.dat'

    CALL WriteTally_Header( FluidMomentum1_FileName )

    FluidMomentum1_Interior = Zero
    FluidMomentum1_Initial  = Zero
    FluidMomentum1_OffGrid  = Zero
    FluidMomentum1_Change   = Zero

    ! --- Fluid Momentum 2 ---

    FluidMomentum2_FileName &
      = TRIM( BaseFileName ) // '_Tally_FluidMomentum2.dat'

    CALL WriteTally_Header( FluidMomentum2_FileName )

    FluidMomentum2_Interior = Zero
    FluidMomentum2_Initial  = Zero
    FluidMomentum2_OffGrid  = Zero
    FluidMomentum2_Change   = Zero

    ! --- Fluid Momentum 3 ---

    FluidMomentum3_FileName &
      = TRIM( BaseFileName ) // '_Tally_FluidMomentum3.dat'

    CALL WriteTally_Header( FluidMomentum3_FileName )

    FluidMomentum3_Interior = Zero
    FluidMomentum3_Initial  = Zero
    FluidMomentum3_OffGrid  = Zero
    FluidMomentum3_Change   = Zero

    ! --- Gravitational Energy ---

    GravitationalEnergy_FileName &
      = TRIM( BaseFileName ) // '_Tally_GravitationalEnergy.dat'

    CALL WriteTally_Header( GravitationalEnergy_FileName )

    GravitationalEnergy_Interior = Zero
    GravitationalEnergy_Initial  = Zero
    GravitationalEnergy_OffGrid  = Zero
    GravitationalEnergy_Change   = Zero

    IF( UnitsActive )THEN

      hc3 = ( PlanckConstant * SpeedOfLight )**3

    ELSE

      hc3 = One

    END IF

  END SUBROUTINE InitializeTally


  SUBROUTINE ComputeTally &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, Time, GE, GX, U, M, &
      SetInitialValues_Option, Verbose_Option )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      Time
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE, &
         iZ_B1(1):iZ_E1(1), &
         1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in) :: &
      U (1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCF)
    REAL(DP), INTENT(in) :: &
      M (1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)
    LOGICAL, INTENT(in), OPTIONAL :: &
      SetInitialValues_Option
    LOGICAL, INTENT(in), OPTIONAL :: &
      Verbose_Option

    LOGICAL :: SetInitialValues
    LOGICAL :: Verbose
    INTEGER :: iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)

    iX_B0 = iZ_B0(2:4); iX_B1 = iZ_B1(2:4)
    iX_E0 = iZ_E0(2:4); iX_E1 = iZ_E1(2:4)

    SetInitialValues = .FALSE.
    IF( PRESENT( SetInitialValues_Option ) )THEN
      SetInitialValues = SetInitialValues_Option
    END IF

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    END IF

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE FROM( GE, GX, U, M )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE HOST( GE, GX, U, M )
#endif

    CALL ComputeTally_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, M )

    CALL ComputeTally_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, GX, U )

    IF( SetInitialValues )THEN

      ! --- Neutrinos ---

      NeutrinoLeptonNumber_Initial = NeutrinoLeptonNumber_Interior
      NeutrinoEnergy_Initial       = NeutrinoEnergy_Interior
      NeutrinoMomentum1_Initial    = NeutrinoMomentum1_Interior
      NeutrinoMomentum2_Initial    = NeutrinoMomentum2_Interior
      NeutrinoMomentum3_Initial    = NeutrinoMomentum3_Interior

      ! --- Fluid ---

      FluidLeptonNumber_Initial    = FluidLeptonNumber_Interior
      FluidEnergy_Initial          = FluidEnergy_Interior
      FluidMomentum1_Initial       = FluidMomentum1_Interior
      FluidMomentum2_Initial       = FluidMomentum2_Interior
      FluidMomentum3_Initial       = FluidMomentum3_Interior

      ! --- Gravitational ---

      GravitationalEnergy_Initial  = GravitationalEnergy_Interior

    END IF

    ! --- Neutrinos ---

    NeutrinoLeptonNumber_Change &
      = NeutrinoLeptonNumber_Interior &
          - ( NeutrinoLeptonNumber_Initial + NeutrinoLeptonNumber_OffGrid )

    NeutrinoEnergy_Change &
      = NeutrinoEnergy_Interior &
          - ( NeutrinoEnergy_Initial       + NeutrinoEnergy_OffGrid       )

    NeutrinoMomentum1_Change &
      = NeutrinoMomentum1_Interior &
          - ( NeutrinoMomentum1_Initial    + NeutrinoMomentum1_OffGrid    )

    NeutrinoMomentum2_Change &
      = NeutrinoMomentum2_Interior &
          - ( NeutrinoMomentum2_Initial    + NeutrinoMomentum2_OffGrid    )

    NeutrinoMomentum3_Change &
      = NeutrinoMomentum3_Interior &
          - ( NeutrinoMomentum3_Initial    + NeutrinoMomentum3_OffGrid    )

    ! --- Fluid ---

    FluidLeptonNumber_Change &
      = FluidLeptonNumber_Interior &
          - ( FluidLeptonNumber_Initial    + FluidLeptonNumber_OffGrid    )

    FluidEnergy_Change &
      = FluidEnergy_Interior &
          - ( FluidEnergy_Initial          + FluidEnergy_OffGrid          )

    FluidMomentum1_Change &
      = FluidMomentum1_Interior &
          - ( FluidMomentum1_Initial       + FluidMomentum1_OffGrid       )

    FluidMomentum2_Change &
      = FluidMomentum2_Interior &
          - ( FluidMomentum2_Initial       + FluidMomentum2_OffGrid       )

    FluidMomentum3_Change &
      = FluidMomentum3_Interior &
          - ( FluidMomentum3_Initial       + FluidMomentum3_OffGrid       )

    ! --- Gravitational ---

    GravitationalEnergy_Change &
      = GravitationalEnergy_Interior &
          - ( GravitationalEnergy_Initial  + GravitationalEnergy_OffGrid  )

    CALL WriteTally( Time )

    IF( Verbose )THEN

      CALL DisplayTally( Time )

    END IF

  END SUBROUTINE ComputeTally


  SUBROUTINE ComputeTally_TwoMoment( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, M )

    INTEGER,  INTENT(in) :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in) :: &
      GE(1:nDOFE, &
         iZ_B1(1):iZ_E1(1), &
         1:nGE)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(in) :: &
      U (1:nDOFX, &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCF)
    REAL(DP), INTENT(in) :: &
      M (1:nDOFZ, &
         iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2), &
         iZ_B1(3):iZ_E1(3), &
         iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS
    INTEGER  :: iNodeE, iNodeX, iNodeZ
    REAL(DP) :: &
      P(1:nDOFX, &
        iZ_B0(2):iZ_E0(2), &
        iZ_B0(3):iZ_E0(3), &
        iZ_B0(4):iZ_E0(4), &
        1:nPF)

    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      DO iNodeX = 1, nDOFX

        CALL ComputePrimitive_Euler_NonRelativistic &
               ( U (iNodeX,iZ2,iZ3,iZ4,iCF_D ),        &
                 U (iNodeX,iZ2,iZ3,iZ4,iCF_S1),        &
                 U (iNodeX,iZ2,iZ3,iZ4,iCF_S2),        &
                 U (iNodeX,iZ2,iZ3,iZ4,iCF_S3),        &
                 U (iNodeX,iZ2,iZ3,iZ4,iCF_E ),        &
                 U (iNodeX,iZ2,iZ3,iZ4,iCF_Ne),        &
                 P (iNodeX,iZ2,iZ3,iZ4,iPF_D ),        &
                 P (iNodeX,iZ2,iZ3,iZ4,iPF_V1),        &
                 P (iNodeX,iZ2,iZ3,iZ4,iPF_V2),        &
                 P (iNodeX,iZ2,iZ3,iZ4,iPF_V3),        &
                 P (iNodeX,iZ2,iZ3,iZ4,iPF_E ),        &
                 P (iNodeX,iZ2,iZ3,iZ4,iPF_Ne),        &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11),  &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22),  &
                 GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33) )

      END DO

    END DO
    END DO
    END DO

    ASSOCIATE &
      ( dZ1 => MeshE    % Width, dZ2 => MeshX(1) % Width, &
        dZ3 => MeshX(2) % Width, dZ4 => MeshX(3) % Width )

    NeutrinoLeptonNumber_Interior = Zero

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        NeutrinoLeptonNumber_Interior                     &
          = NeutrinoLeptonNumber_Interior                 &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4) &
                * Weights_q(iNodeZ)                       &
                * GE(iNodeE,iZ1,iGE_Ep2)                  &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)       &
                * LeptonNumber(iS)                        &
                * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS)

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    NeutrinoLeptonNumber_Interior &
      = FourPi * NeutrinoLeptonNumber_Interior / hc3

    NeutrinoEnergy_Interior    = Zero
    NeutrinoMomentum1_Interior = Zero
    NeutrinoMomentum2_Interior = Zero
    NeutrinoMomentum3_Interior = Zero

    DO iS  = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)

      DO iNodeX = 1, nDOFX
      DO iNodeE = 1, nDOFE

        iNodeZ = (iNodeX-1) * nDOFE + iNodeE

        NeutrinoEnergy_Interior                               &
          = NeutrinoEnergy_Interior                           &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4)     &
                * Weights_q(iNodeZ)                           &
                * GE(iNodeE,iZ1,iGE_Ep3)                      &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)           &
                * ( M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS)        &
                    + P(iNodeX,iZ2,iZ3,iZ4,iPF_V1)            &
                        * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS) &
                    + P(iNodeX,iZ2,iZ3,iZ4,iPF_V2)            &
                        * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS) &
                    + P(iNodeX,iZ2,iZ3,iZ4,iPF_V3)            &
                        * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS) )

        NeutrinoMomentum1_Interior                            &
          = NeutrinoMomentum1_Interior                        &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4)     &
                * Weights_q(iNodeZ)                           &
                * GE(iNodeE,iZ1,iGE_Ep3)                      &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)           &
                * ( M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G1,iS)       &
                    + GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_11)     &
                        * P(iNodeX,iZ2,iZ3,iZ4,iPF_V1)        &
                        * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )

        NeutrinoMomentum2_Interior                            &
          = NeutrinoMomentum2_Interior                        &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4)     &
                * Weights_q(iNodeZ)                           &
                * GE(iNodeE,iZ1,iGE_Ep3)                      &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)           &
                * ( M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G2,iS)       &
                    + GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_22)     &
                        * P(iNodeX,iZ2,iZ3,iZ4,iPF_V2)        &
                        * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )

        NeutrinoMomentum3_Interior                            &
          = NeutrinoMomentum3_Interior                        &
              + dZ1(iZ1) * dZ2(iZ2) * dZ3(iZ3) * dZ4(iZ4)     &
                * Weights_q(iNodeZ)                           &
                * GE(iNodeE,iZ1,iGE_Ep3)                      &
                * GX(iNodeX,iZ2,iZ3,iZ4,iGF_SqrtGm)           &
                * ( M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_G3,iS)       &
                    + GX(iNodeX,iZ2,iZ3,iZ4,iGF_Gm_dd_33)     &
                        * P(iNodeX,iZ2,iZ3,iZ4,iPF_V3)        &
                        * M(iNodeZ,iZ1,iZ2,iZ3,iZ4,iCR_N,iS) )

      END DO
      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

    NeutrinoEnergy_Interior &
      = FourPi * NeutrinoEnergy_Interior    / hc3

    NeutrinoMomentum1_Interior &
      = FourPi * NeutrinoMomentum1_Interior / hc3

    NeutrinoMomentum2_Interior &
      = FourPi * NeutrinoMomentum2_Interior / hc3

    NeutrinoMomentum3_Interior &
      = FourPi * NeutrinoMomentum3_Interior / hc3

    END ASSOCIATE ! dZ1, etc.

  END SUBROUTINE ComputeTally_TwoMoment


  SUBROUTINE IncrementOffGridTally_TwoMoment( dM )

    REAL(DP), INTENT(in) :: dM(2*nCR)

    NeutrinoLeptonNumber_OffGrid &
      = NeutrinoLeptonNumber_OffGrid + FourPi * dM(iCR_N     ) / hc3

    NeutrinoEnergy_OffGrid &
      = NeutrinoEnergy_OffGrid       + FourPi * dM(nCR+iCR_N ) / hc3

    NeutrinoMomentum1_OffGrid &
      = NeutrinoMomentum1_OffGrid    + FourPi * dM(nCR+iCR_G1) / hc3

    NeutrinoMomentum2_OffGrid &
      = NeutrinoMomentum2_OffGrid    + FourPi * dM(nCR+iCR_G2) / hc3

    NeutrinoMomentum3_OffGrid &
      = NeutrinoMomentum3_OffGrid    + FourPi * dM(nCR+iCR_G3) / hc3

  END SUBROUTINE IncrementOffGridTally_TwoMoment


  SUBROUTINE IncrementPositivityLimiterTally_TwoMoment( dM )

    REAL(DP), INTENT(in) :: dM(nCR)

    NeutrinoEnergy_PL    = NeutrinoEnergy_PL    + dM(iCR_N )
    NeutrinoMomentum1_PL = NeutrinoMomentum1_PL + dM(iCR_G1)
    NeutrinoMomentum2_PL = NeutrinoMomentum2_PL + dM(iCR_G2)
    NeutrinoMomentum3_PL = NeutrinoMomentum3_PL + dM(iCR_G3)

  END SUBROUTINE IncrementPositivityLimiterTally_TwoMoment


  SUBROUTINE ComputeTally_Euler( iX_B0, iX_E0, iX_B1, iX_E1, GX, U )

    INTEGER, INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in) :: &
      GX(1:nDOFX, &
         iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3), &
         1:nGF)
    REAL(DP), INTENT(in) :: &
      U (1:nDOFX, &
         iX_B1(1):iX_E1(1), &
         iX_B1(2):iX_E1(2), &
         iX_B1(3):iX_E1(3), &
         1:nCF)

    INTEGER :: iX1, iX2, iX3, iNodeX

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    FluidLeptonNumber_Interior = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        FluidLeptonNumber_Interior                  &
          = FluidLeptonNumber_Interior              &
              + dX1(iX1) * dX2(iX2) * dX3(iX3)      &
                * WeightsX_q(iNodeX)                &
                * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U (iNodeX,iX1,iX2,iX3,iCF_Ne    )

      END DO

    END DO
    END DO
    END DO

    FluidEnergy_Interior = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        FluidEnergy_Interior                        &
          = FluidEnergy_Interior                    &
              + dX1(iX1) * dX2(iX2) * dX3(iX3)      &
                * WeightsX_q(iNodeX)                &
                * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U (iNodeX,iX1,iX2,iX3,iCF_E     )

      END DO

    END DO
    END DO
    END DO

    FluidMomentum1_Interior = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        FluidMomentum1_Interior                     &
          = FluidMomentum1_Interior                 &
              + dX1(iX1) * dX2(iX2) * dX3(iX3)      &
                * WeightsX_q(iNodeX)                &
                * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U (iNodeX,iX1,iX2,iX3,iCF_S1    )

      END DO

    END DO
    END DO
    END DO

    FluidMomentum2_Interior = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        FluidMomentum2_Interior                     &
          = FluidMomentum2_Interior                 &
              + dX1(iX1) * dX2(iX2) * dX3(iX3)      &
                * WeightsX_q(iNodeX)                &
                * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U (iNodeX,iX1,iX2,iX3,iCF_S2    )

      END DO

    END DO
    END DO
    END DO

    FluidMomentum3_Interior = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        FluidMomentum3_Interior                     &
          = FluidMomentum3_Interior                 &
              + dX1(iX1) * dX2(iX2) * dX3(iX3)      &
                * WeightsX_q(iNodeX)                &
                * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U (iNodeX,iX1,iX2,iX3,iCF_S3    )

      END DO

    END DO
    END DO
    END DO

    GravitationalEnergy_Interior = Zero

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        GravitationalEnergy_Interior                &
          = GravitationalEnergy_Interior            &
              + dX1(iX1) * dX2(iX2) * dX3(iX3)      &
                * WeightsX_q(iNodeX) * Half         &
                * GX(iNodeX,iX1,iX2,iX3,iGF_SqrtGm) &
                * U (iNodeX,iX1,iX2,iX3,iCF_D     ) &
                * GX(iNodeX,iX1,iX2,iX3,iGF_Phi_N )

      END DO

    END DO
    END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTally_Euler


  SUBROUTINE IncrementOffGridTally_Euler( dU )

    REAL(DP), INTENT(in) :: dU(nCF)

    FluidMomentum1_Offgrid &
      = FluidMomentum1_OffGrid    + dU(iCF_S1)

    FluidMomentum2_Offgrid &
      = FluidMomentum2_OffGrid    + dU(iCF_S2)

    FluidMomentum3_Offgrid &
      = FluidMomentum3_OffGrid    + dU(iCF_S3)

    FluidEnergy_OffGrid &
      = FluidEnergy_OffGrid       + dU(iCF_E )

    FluidLeptonNumber_OffGrid &
      = FluidLeptonNumber_OffGrid + dU(iCF_Ne)

  END SUBROUTINE IncrementOffGridTally_Euler


  SUBROUTINE WriteTally( Time )

    REAL(DP), INTENT(in) :: Time

    ASSOCIATE( U => UnitsDisplay )

    ! --- Neutrino Lepton Number ---

    CALL WriteTally_Variable &
           ( Time / U % TimeUnit, &
             NeutrinoLeptonNumber_Interior, &
             NeutrinoLeptonNumber_OffGrid , &
             NeutrinoLeptonNumber_Initial , &
             NeutrinoLeptonNumber_Change,   &
             NeutrinoLeptonNumber_FileName )

    ! --- Neutrino Energy ---

    CALL WriteTally_Variable &
           ( Time / U % TimeUnit, &
             NeutrinoEnergy_Interior / U % EnergyGlobalUnit, &
             NeutrinoEnergy_OffGrid  / U % EnergyGlobalUnit, &
             NeutrinoEnergy_Initial  / U % EnergyGlobalUnit, &
             NeutrinoEnergy_Change   / U % EnergyGlobalUnit, &
             NeutrinoEnergy_FileName )

    CALL WriteLimiterTally_Variable &
           ( Time / U % TimeUnit, &
             NeutrinoEnergy_PL / U % EnergyGlobalUnit, &
             NeutrinoEnergy_PL_FileName )

    ! --- Neutrino Momentum 1 ---

    CALL WriteTally_Variable             &
           ( Time / U % TimeUnit,        &
             NeutrinoMomentum1_Interior, &
             NeutrinoMomentum1_OffGrid,  &
             NeutrinoMomentum1_Initial,  &
             NeutrinoMomentum1_Change,   &
             NeutrinoMomentum1_FileName )

    CALL WriteLimiterTally_Variable &
           ( Time / U % TimeUnit, &
             NeutrinoMomentum1_PL, &
             NeutrinoMomentum1_PL_FileName )

    ! --- Neutrino Momentum 2 ---

    CALL WriteTally_Variable             &
           ( Time / U % TimeUnit,        &
             NeutrinoMomentum2_Interior, &
             NeutrinoMomentum2_OffGrid,  &
             NeutrinoMomentum2_Initial,  &
             NeutrinoMomentum2_Change,   &
             NeutrinoMomentum2_FileName )

    CALL WriteLimiterTally_Variable &
           ( Time / U % TimeUnit, &
             NeutrinoMomentum2_PL, &
             NeutrinoMomentum2_PL_FileName )

    ! --- Neutrino Momentum 3 ---

    CALL WriteTally_Variable             &
           ( Time / U % TimeUnit,        &
             NeutrinoMomentum3_Interior, &
             NeutrinoMomentum3_OffGrid,  &
             NeutrinoMomentum3_Initial,  &
             NeutrinoMomentum3_Change,   &
             NeutrinoMomentum3_FileName )

    CALL WriteLimiterTally_Variable &
           ( Time / U % TimeUnit, &
             NeutrinoMomentum3_PL, &
             NeutrinoMomentum3_PL_FileName )

    ! --- Fluid Lepton Number ---

    CALL WriteTally_Variable             &
           ( Time / U % TimeUnit,        &
             FluidLeptonNumber_Interior, &
             FluidLeptonNumber_OffGrid,  &
             FluidLeptonNumber_Initial,  &
             FluidLeptonNumber_Change,   &
             FluidLeptonNumber_FileName )

    ! --- Fluid Energy ---

    CALL WriteTally_Variable &
           ( Time / U % TimeUnit, &
             FluidEnergy_Interior / U % EnergyGlobalUnit, &
             FluidEnergy_OffGrid  / U % EnergyGlobalUnit, &
             FluidEnergy_Initial  / U % EnergyGlobalUnit, &
             FluidEnergy_Change   / U % EnergyGlobalUnit, &
             FluidEnergy_FileName )

    ! --- Fluid Momentum 1 ---

    CALL WriteTally_Variable          &
           ( Time / U % TimeUnit,     &
             FluidMomentum1_Interior, &
             FluidMomentum1_OffGrid,  &
             FluidMomentum1_Initial,  &
             FluidMomentum1_Change,   &
             FluidMomentum1_FileName )

    ! --- Fluid Momentum 2 ---

    CALL WriteTally_Variable          &
           ( Time / U % TimeUnit,     &
             FluidMomentum2_Interior, &
             FluidMomentum2_OffGrid,  &
             FluidMomentum2_Initial,  &
             FluidMomentum2_Change,   &
             FluidMomentum2_FileName )

    ! --- Fluid Momentum 3 ---

    CALL WriteTally_Variable          &
           ( Time / U % TimeUnit,     &
             FluidMomentum3_Interior, &
             FluidMomentum3_OffGrid,  &
             FluidMomentum3_Initial,  &
             FluidMomentum3_Change,   &
             FluidMomentum3_FileName )

    ! --- Gravitational Energy ---

    CALL WriteTally_Variable &
           ( Time / U % TimeUnit, &
             GravitationalEnergy_Interior / U % EnergyGlobalUnit, &
             GravitationalEnergy_OffGrid  / U % EnergyGlobalUnit, &
             GravitationalEnergy_Initial  / U % EnergyGlobalUnit, &
             GravitationalEnergy_Change   / U % EnergyGlobalUnit, &
             GravitationalEnergy_FileName )

    END ASSOCIATE

  END SUBROUTINE WriteTally


  SUBROUTINE WriteTally_Header( FileName )

    CHARACTER(*), INTENT(in) :: FileName

    INTEGER :: FileUnit

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ) )

    WRITE( FileUnit, '(5(A20,x))' ) &
      'Time', 'Interior', 'Off Grid', 'Initial', 'Change'

    CLOSE( FileUnit )

  END SUBROUTINE WriteTally_Header


  SUBROUTINE WriteLimiterTally_Header( FileName )

    CHARACTER(*), INTENT(in) :: FileName

    INTEGER :: FileUnit

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ) )

    WRITE( FileUnit, '(2(A20,x))' ) 'Time', 'Change'

    CLOSE( FileUnit )

  END SUBROUTINE WriteLimiterTally_Header


  SUBROUTINE WriteTally_Variable &
    ( Time, Interior, OffGrid, Initial, Change, FileName )

    REAL(DP)    , INTENT(in) :: Time, Interior, OffGrid, Initial, Change
    CHARACTER(*), INTENT(in) :: FileName

    INTEGER :: FileUnit

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(5(ES20.12,x))' ) &
      Time, Interior, OffGrid , Initial , Change

    CLOSE( FileUnit )

  END SUBROUTINE WriteTally_Variable


  SUBROUTINE WriteLimiterTally_Variable( Time, Change, FileName )

    REAL(DP)    , INTENT(in) :: Time, Change
    CHARACTER(*), INTENT(in) :: FileName

    INTEGER :: FileUnit

    OPEN( NEWUNIT = FileUnit, FILE = TRIM( FileName ), &
          POSITION = 'APPEND', ACTION = 'WRITE' )

    WRITE( FileUnit, '(2(ES20.12,x))' ) Time, Change

    CLOSE( FileUnit )

  END SUBROUTINE WriteLimiterTally_Variable


  SUBROUTINE DisplayTally( Time )

    REAL(DP), INTENT(in) :: Time

    ASSOCIATE( U => UnitsDisplay )

    WRITE(*,*)
    WRITE(*,'(A8,A,ES8.2E2)') '', 'Two-Moment Tally O(v/c). t = ', &
      Time / U % TimeUnit
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Lepton Number Interior.: ', &
      NeutrinoLeptonNumber_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Lepton Number Initial..: ', &
      NeutrinoLeptonNumber_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Lepton Number Off Grid.: ', &
      NeutrinoLeptonNumber_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Lepton Number Change...: ', &
      NeutrinoLeptonNumber_Change
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Energy Interior........: ', &
      NeutrinoEnergy_Interior / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Energy Initial.........: ', &
      NeutrinoEnergy_Initial  / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Energy Off Grid........: ', &
      NeutrinoEnergy_OffGrid  / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Energy Change..........: ', &
      NeutrinoEnergy_Change   / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Energy PL..............: ', &
      NeutrinoEnergy_PL       / U % EnergyGlobalUnit
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 1 Interior....: ', &
      NeutrinoMomentum1_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 1 Initial.....: ', &
      NeutrinoMomentum1_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 1 Off Grid....: ', &
      NeutrinoMomentum1_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 1 Change......: ', &
      NeutrinoMomentum1_Change
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 1 PL..........: ', &
      NeutrinoMomentum1_PL
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 2 Interior....: ', &
      NeutrinoMomentum2_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 2 Initial.....: ', &
      NeutrinoMomentum2_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 2 Off Grid....: ', &
      NeutrinoMomentum2_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 2 Change......: ', &
      NeutrinoMomentum2_Change
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 2 PL..........: ', &
      NeutrinoMomentum2_PL
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 3 Interior....: ', &
      NeutrinoMomentum3_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 3 Initial.....: ', &
      NeutrinoMomentum3_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 3 Off Grid....: ', &
      NeutrinoMomentum3_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 3 Change......: ', &
      NeutrinoMomentum3_Change
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Neutrino Momentum 3 PL..........: ', &
      NeutrinoMomentum3_PL
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Lepton Number Interior....: ', &
      FluidLeptonNumber_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Lepton Number Initial.....: ', &
      FluidLeptonNumber_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Lepton Number Off Grid....: ', &
      FluidLeptonNumber_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Lepton Number Change......: ', &
      FluidLeptonNumber_Change
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Energy Interior...........: ', &
      FluidEnergy_Interior / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Energy Initial............: ', &
      FluidEnergy_Initial  / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Energy Off Grid...........: ', &
      FluidEnergy_OffGrid  / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Energy Change.............: ', &
      FluidEnergy_Change   / U % EnergyGlobalUnit
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 1 Interior.......: ', &
      FluidMomentum1_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 1 Initial........: ', &
      FluidMomentum1_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 1 Off Grid.......: ', &
      FluidMomentum1_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 1 Change.........: ', &
      FluidMomentum1_Change
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 2 Interior.......: ', &
      FluidMomentum2_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 2 Initial........: ', &
      FluidMomentum2_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 2 Off Grid.......: ', &
      FluidMomentum2_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 2 Change.........: ', &
      FluidMomentum2_Change
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 3 Interior.......: ', &
      FluidMomentum3_Interior
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 3 Initial........: ', &
      FluidMomentum3_Initial
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 3 Off Grid.......: ', &
      FluidMomentum3_OffGrid
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Fluid Momentum 3 Change.........: ', &
      FluidMomentum3_Change
    WRITE(*,*)
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Gravitational Energy Interior...: ', &
      GravitationalEnergy_Interior / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Gravitational Energy Initial....: ', &
      GravitationalEnergy_Initial  / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Gravitational Energy Off Grid...: ', &
      GravitationalEnergy_OffGrid  / U % EnergyGlobalUnit
    WRITE(*,'(A6,A40,ES16.8E3)') &
      '', 'Gravitational Energy Change.....: ', &
      GravitationalEnergy_Change   / U % EnergyGlobalUnit
    WRITE(*,*)

    END ASSOCIATE

  END SUBROUTINE DisplayTally


  SUBROUTINE FinalizeTally

  END SUBROUTINE FinalizeTally


END MODULE TwoMoment_TallyModule
