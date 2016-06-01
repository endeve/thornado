MODULE RiemannSolverModule

  USE KindModule, ONLY: &
    DP
  USE FluidFieldsModule, ONLY: &
    iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: FluidRiemannSolver     = 'LLF'
  CHARACTER(32) :: RadiationRiemannSolver = 'LLF'

  PROCEDURE (RiemannSolver), POINTER, PUBLIC :: &
    NumericalFlux_Fluid     => NULL(), &
    NumericalFlux_Radiation => NULL()

  INTERFACE
    PURE FUNCTION RiemannSolver &
      ( u_L, u_R, Flux_L, Flux_R, a, aP, aM, aC, nF )
      USE KindModule, ONLY: DP
      INTEGER,                   INTENT(in) :: nF
      REAL(DP), DIMENSION(1:nF), INTENT(in) :: u_L, u_R, flux_L, flux_R
      REAL(DP),                  INTENT(in) :: a, aP, aM, aC
      REAL(DP), DIMENSION(1:nF)             :: RiemannSolver
    END FUNCTION RiemannSolver
  END INTERFACE

  PUBLIC :: InitializeRiemannSolvers
  PUBLIC :: FinalizeRiemannSolvers

CONTAINS


  SUBROUTINE InitializeRiemannSolvers &
               ( FluidRiemannSolver_Option, RadiationRiemannSolver_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      FluidRiemannSolver_Option, &
      RadiationRiemannSolver_Option

    IF( PRESENT( FluidRiemannSolver_Option ) )THEN
      FluidRiemannSolver = FluidRiemannSolver_Option
    END IF

    WRITE(*,*)
    WRITE(*,'(A5,A22,A)') &
      '', 'Fluid Riemann Solver: ', TRIM( FluidRiemannSolver )
    WRITE(*,'(A5,A21)') &
      '', '---------------------'

    SELECT CASE ( TRIM( FluidRiemannSolver ) )
      CASE ( 'LLF' )
        NumericalFlux_Fluid => NumericalFlux_LLF
      CASE ( 'HLL' )
        NumericalFlux_Fluid => NumericalFlux_HLL
      CASE ( 'HLLC' )
        NumericalFlux_Fluid => NumericalFlux_HLLC
      CASE DEFAULT
        NumericalFlux_Fluid => NumericalFlux_LLF
    END SELECT

    IF( PRESENT( RadiationRiemannSolver_Option ) )THEN
      RadiationRiemannSolver = RadiationRiemannSolver_Option
    END IF

    WRITE(*,*)
    WRITE(*,'(A5,A26,A)') &
      '', 'Radiation Riemann Solver: ', TRIM( RadiationRiemannSolver )
    WRITE(*,'(A5,A25)') &
      '', '-------------------------'

    SELECT CASE( TRIM( RadiationRiemannSolver ) )
      CASE ( 'LLF' )
        NumericalFlux_Radiation => NumericalFlux_LLF
      CASE ( 'HLL' )
        NumericalFlux_Radiation => NumericalFlux_HLL
      CASE DEFAULT
        NumericalFlux_Radiation => NumericalFlux_LLF
    END SELECT

  END SUBROUTINE InitializeRiemannSolvers


  SUBROUTINE FinalizeRiemannSolvers

    NULLIFY( NumericalFlux_Fluid )
    NULLIFY( NumericalFlux_Radiation )

  END SUBROUTINE FinalizeRiemannSolvers


  PURE FUNCTION NumericalFlux_LLF &
    ( u_L, u_R, Flux_L, Flux_R, alpha, alpha_P, alpha_M, alpha_C, nF )

    ! --- Local Lax-Friedrichs Flux ---

    INTEGER,                   INTENT(in) :: nF
    REAL(DP), DIMENSION(1:nF), INTENT(in) :: u_L, u_R, flux_L, flux_R
    REAL(DP),                  INTENT(in) :: alpha, alpha_P, alpha_M, alpha_C
    REAL(DP), DIMENSION(1:nF)             :: NumericalFlux_LLF

    NumericalFlux_LLF &
      = 0.5_DP * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF


  PURE FUNCTION NumericalFlux_HLL &
    ( u_L, u_R, Flux_L, Flux_R, alpha, alpha_P, alpha_M, alpha_C, nF )

    ! --- Harten-Lax-van Leer Flux ---

    INTEGER,                   INTENT(in) :: nF
    REAL(DP), DIMENSION(1:nF), INTENT(in) :: u_L, u_R, flux_L, flux_R
    REAL(DP),                  INTENT(in) :: alpha, alpha_P, alpha_M, alpha_C
    REAL(DP), DIMENSION(1:nF)             :: NumericalFlux_HLL

    NumericalFlux_HLL &
      = ( alpha_P * flux_L + alpha_M * flux_R &
            - alpha_P * alpha_M * ( u_R - u_L ) ) / ( alpha_P + alpha_M )

    RETURN
  END FUNCTION NumericalFlux_HLL


  PURE FUNCTION NumericalFlux_HLLC &
    ( u_L, u_R, Flux_L, Flux_R, alpha, alpha_P, alpha_M, alpha_C, nF )

    ! --- Harten-Lax-van Leer Contact Flux ---
    ! --- Non-Relativistic Hydrodynamics -----

    INTEGER,                   INTENT(in) :: nF
    REAL(DP), DIMENSION(1:nF), INTENT(in) :: u_L, u_R, flux_L, flux_R
    REAL(DP),                  INTENT(in) :: alpha, alpha_P, alpha_M, alpha_C
    REAL(DP), DIMENSION(1:nF)             :: NumericalFlux_HLLC

    REAL(DP)                  :: D, V1, V2, V3, P, E
    REAL(DP), DIMENSION(1:nF) :: TMP

    IF( alpha_M .EQ. 0.0_DP )THEN

      NumericalFlux_HLLC = Flux_L

    ELSEIF( alpha_P .EQ. 0.0_DP )THEN

      NumericalFlux_HLLC = Flux_R

    ELSE

      IF( alpha_C .GE. 0.0_DP )THEN

        TMP = Flux_L + alpha_M * u_L

        D  = TMP(iCF_D) / ( alpha_C + alpha_M )
        V1 = alpha_C
        V2 = TMP(iCF_S2) / TMP(iCF_D)
        V3 = TMP(iCF_S3) / TMP(iCF_D)
        P  = TMP(iCF_S1) - alpha_C * TMP(iCF_D)
        E  = ( TMP(iCF_E) - alpha_C * P ) / ( alpha_C + alpha_M )

      ELSE

        TMP = Flux_R - alpha_P * u_R

        D  = TMP(iCF_D) / ( alpha_C - alpha_P )
        V1 = alpha_C
        V2 = TMP(iCF_S2) / TMP(iCF_D)
        V3 = TMP(iCF_S3) / TMP(iCF_D)
        P  = TMP(iCF_S1) - alpha_C * TMP(iCF_D)
        E  = ( TMP(iCF_E) - alpha_C * P ) / ( alpha_C - alpha_P )

      END IF

      NumericalFlux_HLLC(iCF_D) &
        = D * V1
      NumericalFlux_HLLC(iCF_S1) &
        = D * V1 * V1 + P
      NumericalFlux_HLLC(iCF_S2) &
        = D * V2 * V1
      NumericalFlux_HLLC(iCF_S3) &
        = D * V3 * V1
      NumericalFlux_HLLC(iCF_E) &
        = ( E + P ) * V1

    END IF

    RETURN
  END FUNCTION NumericalFlux_HLLC


END MODULE RiemannSolverModule
