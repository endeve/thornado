MODULE FluidRadiationCouplingModule

  USE KindModule, ONLY: &
    DP
  USE FluidRadiationCouplingSolutionModule_ConstantOpacities, ONLY: &
    CoupleFluidRadiation_ConstantOpacities
  USE FluidRadiationCouplingSolutionModule_ThermalReservoir, ONLY: &
    CoupleFluidRadiation_ThermalReservoir
  USE FluidRadiationCouplingSolutionModule_EmissionAbsorption, ONLY: &
    CoupleFluidRadiation_EmissionAbsorption
  USE FluidRadiationCouplingSolutionModule_ElasticScattering, ONLY: &
    CoupleFluidRadiation_ElasticScattering
  USE FluidRadiationCouplingSolutionModule_NES, ONLY: &
    CoupleFluidRadiation_NES
  USE FluidRadiationCouplingSolutionModule_B85, ONLY: &
    CoupleFluidRadiation_B85

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: FluidRadiationCoupling = 'Dummy'

  PROCEDURE (CouplingProcedure), POINTER, PUBLIC :: &
    CoupleFluidRadiation => NULL(), &
    ComputeImplicitIncrement_FluidRadiation => NULL()

  INTERFACE
    SUBROUTINE CouplingProcedure &
                 ( dt, iX_B0, iX_E0, iX_B1, iX_E1, U_F, dU_F, U_R, dU_R, &
                   EvolveFluid_Option )
      USE KindModule, ONLY: DP
      REAL(DP), INTENT(in)  :: &
        dt
      INTEGER,  INTENT(in)  :: &
        iX_B0(3), iX_B1(3), iX_E0(3), iX_E1(3)
      REAL(DP), INTENT(in)  :: &
        U_F (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      REAL(DP), INTENT(out) :: &
        dU_F(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
      REAL(DP), INTENT(in)  :: &
        U_R (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
      REAL(DP), INTENT(out) :: &
        dU_R(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)
      LOGICAL,  INTENT(in), OPTIONAL :: &
        EvolveFluid_Option
    END SUBROUTINE CouplingProcedure
  END INTERFACE

  PUBLIC :: InitializeFluidRadiationCoupling
  PUBLIC :: FinalizeFluidRadiationCoupling

CONTAINS


  SUBROUTINE InitializeFluidRadiationCoupling &
               ( FluidRadiationCoupling_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      FluidRadiationCoupling_Option

    IF( PRESENT( FluidRadiationCoupling_Option ) )THEN
      FluidRadiationCoupling = FluidRadiationCoupling_Option
    END IF

    SELECT CASE ( TRIM( FluidRadiationCoupling ) )

      CASE( 'ConstantOpacities' )

        CoupleFluidRadiation &
          => CoupleFluidRadiation_ConstantOpacities
        ComputeImplicitIncrement_FluidRadiation &
          => CoupleFluidRadiation_ConstantOpacities

      CASE( 'ThermalReservoir' )

        CoupleFluidRadiation &
          => CoupleFluidRadiation_ThermalReservoir

      CASE( 'EmissionAbsorption' )

        CoupleFluidRadiation &
          => CoupleFluidRadiation_EmissionAbsorption

      CASE( 'ElasticScattering' )

        CoupleFluidRadiation &
          => CoupleFluidRadiation_ElasticScattering

      CASE( 'NES' )

        CoupleFluidRadiation &
          => CoupleFluidRadiation_NES

      CASE( 'B85' )

        CoupleFluidRadiation &
          => CoupleFluidRadiation_B85

      CASE( 'Complete' )

      CASE DEFAULT

        CoupleFluidRadiation &
          => CoupleFluidRadiation_Dummy

    END SELECT

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    NULLIFY( CoupleFluidRadiation )
    NULLIFY( ComputeImplicitIncrement_FluidRadiation )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation_Dummy &
               ( dt, iX_B0, iX_E0, iX_B1, iX_E1, U_F, dU_F, U_R, dU_R, &
                 EvolveFluid_Option )

    REAL(DP), INTENT(in)  :: &
      dt
    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_B1(3), iX_E0(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U_F (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(out) :: &
      dU_F(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in)  :: &
      U_R (1:,1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:,1:)
    REAL(DP), INTENT(out) :: &
      dU_R(1:,1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      EvolveFluid_Option

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'FluidRadiationCouplingModule: CoupleFluidRadiation_Dummy'
    WRITE(*,*)

  END SUBROUTINE CoupleFluidRadiation_Dummy


END MODULE FluidRadiationCouplingModule
