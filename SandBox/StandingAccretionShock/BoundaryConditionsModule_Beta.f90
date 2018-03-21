MODULE ApplicationBoundaryConditionsModule_Beta  

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    xL, xR, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber,  &
    WriteVector
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_S1, nCF, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, nPF, &
    uAF, iAF_P, iAF_T,  iAF_Ye, iAF_S,  iAF_E, iAF_Gm, iAF_Cs, nAF
  USE EquationOfStateModule, ONLY: &
    ComputeInternalEnergyDensityFromPressure, &
    ComputeAuxiliary_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    Primitive

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyApplicationBoundaryConditions_Fluid_X1

CONTAINS

  SUBROUTINE ApplyApplicationBoundaryConditions_Fluid_X1

    SELECT CASE ( TRIM( ProgramName ) )
      
      CASE ( 'ApplicationDriver' )

        CALL ApplyBC_Fluid_X1_ApplicationDriver

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A38)') &
          '', 'Application Boundary Condition Missing'
        WRITE(*,*)
        WRITE(*,'(A7,A)') &
          '', TRIM( ProgramName )
        WRITE(*,*)
        STOP

    END SELECT

  END SUBROUTINE ApplyApplicationBoundaryConditions_Fluid_X1

  SUBROUTINE ApplyBC_Fluid_X1_ApplicationDriver

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeX_0, iNodeX_1, iNodeX_q, i_q
    INTEGER, PARAMETER :: N_q = 1
    REAL(DP) :: X1_0, X1_1, Beta, Gamma, Numerator_D, Numerator_E, Denominator_D, Denominator_E, K_D, K_E
    REAL(DP), DIMENSION(N_q * nNodesX(1))   :: X1_q, D_q, E_q
    REAL(DP), DIMENSION(nCF) :: uCF_q

    Gamma = 4.0_DP / 3.0_DP

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        ! --- Inner Boundary ---

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
        ! Pull coeffecient calculations into separate subroutine. 
            i_q = 1
            Numerator_D = 0
            Denominator_D = 0

            Numerator_E = 0
            Denominator_E = 0

            K_D = 0
            K_E = 0 

            DO iX1 = 1, N_q
              DO iNodeX1 = 1, nNodesX(1) 
              
                 X1_q(i_q) = NodeCoordinate( MeshX(1), iX1, iNodeX1)
                
                 iNodeX_q = NodeNumberX(iNodeX1, iNodeX2, iNodeX3)
                
                 uCF_q = uCF(iNodeX_q, iX1, iX2, iX3,:)

                 D_q(i_q) = uCF_q(iCF_D)

                 E_q(i_q) = uCF_q(iCF_E)
              
                 Numerator_D = Numerator_D + D_q(i_q)/X1_q(i_q)**3
                 Denominator_D = Denominator_D + 1.0_DP/X1_q(i_q)**6
                                 
                 Numerator_E = Numerator_E + E_q(i_q)/X1_q(i_q)**Beta
                 Denominator_E = Denominator_E + 1.0_DP/X1_q(i_q)**(2 * Beta) 

                 i_q = i_q + 1

              END DO
            END DO
               
            K_D = Numerator_D/Denominator_D
            K_E = Numerator_P/Denominator_E

            DO iNodeX1 = 1, nNodesX(1)

              iNodeX_0 = NodeNumberX( iNodeX1, iNodeX2, iNodeX3)

              X1_0 = NodeCoordinate( MeshX(1), 0, iNodeX1 )
         
              uCF(iNodeX_0,0,iX2,iX3,iCF_D) &
                = K_D / X1_0**3

              uCF(iNodeX_0,0,iX2,iX3,iCF_E) &
                = K_E / X1_0**Beta

            END DO           
          END DO
        END DO

        CALL ComputeInternalEnergyDensityFromPressure &
               ( uPF(:,0,iX2,iX3,iPF_D ), uAF(:,0,iX2,iX3,iAF_P), &
                 uAF(:,0,iX2,iX3,iAF_Ye), uPF(:,0,iX2,iX3,iPF_E) )

        CALL ComputeAuxiliary_Fluid &
               ( uPF(:,0,iX2,iX3,iPF_D ), uPF(:,0,iX2,iX3,iPF_E ), &
                 uPF(:,0,iX2,iX3,iPF_Ne), uAF(:,0,iX2,iX3,iAF_P ), &
                 uAF(:,0,iX2,iX3,iAF_T ), uAF(:,0,iX2,iX3,iAF_Ye), &
                 uAF(:,0,iX2,iX3,iAF_S ), uAF(:,0,iX2,iX3,iAF_E ), &
                 uAF(:,0,iX2,iX3,iAF_Gm), uAF(:,0,iX2,iX3,iAF_Cs) )

      END DO
    END DO

  END SUBROUTINE ApplyBC_Fluid_X1_StandingAccretionShock1D

END MODULE ApplicationBoundaryConditionsModule
