MODULE FluidRadiationCouplingSolutionModule_Implicit

  USE KindModule, ONLY: &
    DP, Pi
  USE UnitsModule, ONLY: &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE UtilitiesModule, ONLY: &
    WriteVector
  USE FluidFieldsModule, ONLY: &
    nPF, nAF
  USE RadiationFieldsModule, ONLY: &
    nPR
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeWeights, &
    InitializeFluidFields, &
    FinalizeFluidFields, &
    InitializeRadiationFields, &
    FinalizeRadiationFields

  IMPLICIT NONE
  PRIVATE

  INTEGER :: nNodesX_G, nNodesE_G
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Chi, f_FD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N

  PUBLIC :: CoupleFluidRadiation_Implicit_EmissionAbsorption

CONTAINS


  SUBROUTINE CoupleFluidRadiation_Implicit_EmissionAbsorption &
               ( dt, iX_Begin, iX_End )

    REAL(DP),              INTENT(in) :: dt
    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    CALL InitializeFluidRadiationCoupling

    WRITE(*,'(A5,A)') &
      '', 'CoupleFluidRadiation_Implicit_EmissionAbsorption'

    CALL FinalizeFluidRadiationCoupling

    STOP

  END SUBROUTINE CoupleFluidRadiation_Implicit_EmissionAbsorption


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    CALL WriteVector( SIZE( E_N ), E_N / MeV, 'E.dat' )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    CALL WriteVector( SIZE( W2_N ), W2_N / MeV**3, 'W2.dat' )
    CALL WriteVector( SIZE( W3_N ), W3_N / MeV**4, 'W3.dat' )

    ALLOCATE( uPF_N(nPF, nNodesX_G), uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    CALL WriteVector( SIZE( uPR_N(:,1,1) ), uPR_N(:,1,1), 'N.dat' )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N, uPF_N, uAF_N, uPR_N )

  END SUBROUTINE FinalizeFluidRadiationCoupling


END MODULE FluidRadiationCouplingSolutionModule_Implicit
