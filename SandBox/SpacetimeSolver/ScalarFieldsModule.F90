MODULE ScalarFieldsModule
  
  USE KindModule, ONLY: &
    DP, One, SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX
  USE MeshModule, ONLY: &
    MeshX

  IMPLICIT NONE
  PRIVATE

  
  INTEGER, PUBLIC, PARAMETER :: iSF_u    = 1  ! 
  INTEGER, PUBLIC, PARAMETER :: iSF_v    = 2  ! 
  INTEGER, PUBLIC, PARAMETER :: nSF      = 2  ! n Scalar Fields

  CHARACTER(32), DIMENSION(nSF), PUBLIC, PARAMETER :: &
    namesSF = [ 'U                                           ', &
                'V                                           ' ]

  REAL(DP), DIMENSION(nSF), PUBLIC :: unitsSF

  REAL(DP), ALLOCATABLE, PUBLIC :: uSF( :, :, :, :, : )

  PUBLIC :: CreateScalarFields
  PUBLIC :: DestroyScalarFields
  PUBLIC :: Flux_X1
  PUBLIC :: NumericalFlux_X1_Upwind
  PUBLIC :: Source_X1
  PUBLIC :: ComputeTimeStep_ScalarWave

Contains

  SUBROUTINE CreateScalarFields &
    ( nX, swX )
    
    INTEGER,  INTENT(in)   :: nX(3), swX(3)

    ALLOCATE &
         ( uSF ( 1 : nDOFX, &
                 1 - sWX(1) : nX(1) + swX(1), &
                 1 - sWX(2) : nX(2) + swX(2), & 
                 1 - sWX(3) : nX(3) + swX(3), &
                 1 : nSF ) )

    unitsSF = One

    !-- Initialize Fields --
    
    uSF( :, :, :, :, iSF_u )   = 0.0_DP
    uSF( :, :, :, :, iSF_v )   = 0.0_DP

  END SUBROUTINE CreateScalarFields

  SUBROUTINE DestroyScalarFields

    DEALLOCATE( uSF )

  END SUBROUTINE DestroyScalarFields

  PURE FUNCTION Flux_X1 ( U, V )

    REAL(DP)             :: Flux_X1(1:nSF)
    REAL(DP), INTENT(in) :: U, V 

    Flux_X1 (iSF_u) = 1.0_DP * U
    Flux_X1 (iSF_v) = - 1.0_DP * V

  END FUNCTION Flux_X1

  PURE FUNCTION NumericalFlux_X1_Upwind &
    (uL, uR, fL, fR, lambda )

    REAL(DP), INTENT(in) :: uL(nSF), uR(nSF), fL(nSF), fR(nSF), lambda (nSF)
    
    REAL(DP) :: NumericalFlux_X1_Upwind(nSF)
    integer  :: iSF

    DO iSF=1, nSF 
      NumericalFlux_X1_Upwind(iSF) = 0.5_DP*(fL(iSF) + fR(iSF) - ABS(lambda(iSF))*(uR(iSF) - uL(iSF)))
    END DO


  END FUNCTION NumericalFlux_X1_Upwind

  PURE FUNCTION Source_X1 ( U, V )

    REAL(DP)             :: Source_X1(1:nSF)
    REAL(DP), INTENT(in) :: U, V 

    Source_X1 (iSF_u) = V
    Source_X1 (iSF_v) = 0.0_DP

  END FUNCTION Source_X1

  SUBROUTINE ComputeTimeStep_ScalarWave &
               ( iX_B0, iX_E0, iX_B1, iX_E1, U, CFL, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)  :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

    INTEGER  :: iX1, iX2, iX3, iNodeX, iDimX
    REAL(DP) :: dX(3), dt

    ASSOCIATE &
      ( dX1 => MeshX(1) % Width, &
        dX2 => MeshX(2) % Width, &
        dX3 => MeshX(3) % Width )

    TimeStep = HUGE( One )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        dX(1) = dX1(iX1)
        dX(2) = dX2(iX2)
        dX(3) = dX3(iX3)

        DO iDimX = 1, nDimsX

          dt = dX(iDimX)

          TimeStep = MIN( TimeStep, dt )

        END DO

      END DO

    END DO
    END DO
    END DO

    TimeStep = MAX( CFL * TimeStep, SqrtTiny )

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTimeStep_ScalarWave

END MODULE ScalarFieldsModule
