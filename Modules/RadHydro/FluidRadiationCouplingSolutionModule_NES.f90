MODULE FluidRadiationCouplingSolutionModule_NES

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE UtilitiesModule, ONLY: &
    WriteVector, &
    WriteMatrix
  USE FluidFieldsModule, ONLY: &
    nPF, iPF_D, &
    nAF, iAF_T, iAF_Me
  USE RadiationFieldsModule, ONLY: &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputeConservedMoments, &
    ComputePrimitiveMoments
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeWeights, &
    InitializeFluidFields, &
    InitializeRadiationFields, &
    FinalizeFluidFields, &
    FinalizeRadiationFields
  USE OpacityModule, ONLY: &
    ComputeScatteringOpacity_NES

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: EvolveFluid
  INTEGER :: nNodesX_G, nNodesE_G
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: R0_In, R0_Out

  PUBLIC :: CoupleFluidRadiation_NES

CONTAINS


  SUBROUTINE CoupleFluidRadiation_NES &
               ( dt, iX_Begin, iX_End, EvolveFluid_Option )

    REAL(DP),              INTENT(in) :: dt
    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    LOGICAL,               INTENT(in), OPTIONAL :: EvolveFluid_Option

    EvolveFluid = .TRUE.
    IF( PRESENT( EvolveFluid_Option ) ) &
      EvolveFluid = EvolveFluid_Option

    CALL ComputePrimitiveMoments &
           ( iX_Begin = iX_Begin, iX_End = iX_End )

    CALL InitializeFluidRadiationCoupling

    CALL CoupleFluidRadiation( dt )

    CALL FinalizeFluidRadiationCoupling

    CALL ComputeConservedMoments &
           ( iX_Begin = iX_Begin, iX_End = iX_End )

  END SUBROUTINE CoupleFluidRadiation_NES


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    ALLOCATE( uPF_N(nPF, nNodesX_G) )
    ALLOCATE( uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    ALLOCATE &
      ( R0_In (nNodesE_G,nNodesE_G,nNodesX_G), &
        R0_Out(nNodesE_G,nNodesE_G,nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( uPF_N, uAF_N, uPR_N )
    DEALLOCATE( R0_In, R0_Out )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation( dt )

    REAL(DP), INTENT(in) :: dt

    ! --- Compute Scattering Opacities ---

    CALL SetRates

  END SUBROUTINE CoupleFluidRadiation


  SUBROUTINE SetRates

    ASSOCIATE &
      ( kT => BoltzmannConstant * uAF_N(iAF_T,1:nNodesX_G) )

    ASSOCIATE &
      ( T_N   => uAF_N(iAF_T, 1:nNodesX_G), &
        Eta_N => uAF_N(iAF_Me,1:nNodesX_G) / kT )

    CALL ComputeScatteringOpacity_NES &
           ( E_N, T_N, Eta_N, R0_In, R0_Out )

!!$    CALL WriteVector( SIZE(E_N), E_N, 'E_N.dat' )
!!$    CALL WriteMatrix( SIZE(E_N), SIZE(E_N), R0_In (:,:,1), 'R0_In.dat' )
!!$    CALL WriteMatrix( SIZE(E_N), SIZE(E_N), R0_Out(:,:,1), 'R0_Out.dat' )

    END ASSOCIATE ! T_N, etc.
    END ASSOCIATE ! kT

  END SUBROUTINE SetRates


END MODULE FluidRadiationCouplingSolutionModule_NES
