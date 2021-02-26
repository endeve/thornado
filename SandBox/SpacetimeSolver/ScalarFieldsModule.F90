MODULE ScalarFieldsModule
  
  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX

  IMPLICIT NONE
  PRIVATE

  
  INTEGER, PUBLIC, PARAMETER :: iSF_u    = 1  ! 
  INTEGER, PUBLIC, PARAMETER :: iSF_v    = 2  ! 
  INTEGER, PUBLIC, PARAMETER :: nSF      = 2  ! n Scalar Fields

  CHARACTER(32), DIMENSION(nSF), PUBLIC, PARAMETER :: &
    namesSF = [ 'U                                           ', &
                'V                                           ' ]


  REAL(DP), ALLOCATABLE, PUBLIC :: uSF( :, :, :, :, : )

  PUBLIC :: CreateScalarFields
  PUBLIC :: DestroyScalarFields
  PUBLIC :: Flux_X1
  PUBLIC :: NumericalFlux_X1_Upwind

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

    !-- Set Initial Conditions --
    
    uSF( :, :, :, :, iSF_u )   = 0.0_DP
    uSF( :, :, :, :, iSF_v )   = 0.0_DP

  END SUBROUTINE CreateScalarFields

  SUBROUTINE DestroyScalarFields

    DEALLOCATE( uSF )

  END SUBROUTINE DestroyScalarFields

  PURE FUNCTION Flux_X1 ( U, V )

    REAL(DP)             :: Flux_X1(1:nSF)
    REAL(DP), INTENT(in) :: U, V 

    Flux_X1 (iSF_u) = 1.0_dp * U
    Flux_X1 (iSF_v) = - 1.0_dp * V

  END FUNCTION Flux_X1

  PURE FUNCTION NumericalFlux_X1_Upwind &
    (uL, uR, fL, fR, lambda )

    REAL(DP), INTENT(in) :: uL(nSF), uR(nSF), fL(nSF), fR(nSF), lambda (nSF)
    
    REAL(DP) :: NumericalFlux_X1_Upwind(nSF)

    DO iSF=1, nSF 
      NumericalFlux_X1_Upwind(iSF) = 0.5_dp*(fL(iSF) + fR(iSF) - ABS(lambda(iSF))*(uR(iSF) - uL(iSF)))
    END DO


  END FUNCTION NumericalFlux_X1_Upwind

END MODULE ScalarFieldsModule
