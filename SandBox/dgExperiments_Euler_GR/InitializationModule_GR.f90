MODULE InitializationModule_GR

  USE KindModule, ONLY: &
    DP, Zero, Half, One, TwoPi, Pi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nNodesX, &
    nDOFX
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uAF, iAF_P, iAF_Gm
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputeConserved_GR
  USE DataFileReaderModule, ONLY: &
    ReadData, ReadParameters
  USE UtilitiesModule, ONLY: &
    Locate
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_RiemannProblem
  PUBLIC :: InitializeFields_StandingAccretionShock

  INTEGER :: iL_Debug
CONTAINS

  SUBROUTINE InitializeFields_RiemannProblem &
    ( D_L, V_L, P_L, D_R, V_R, P_R, X_D_Option )

    REAL(DP), INTENT(in) :: D_L, P_L, D_R, P_R
    REAL(DP), INTENT(in) :: V_L(3), V_R(3)
    REAL(DP), INTENT(in), OPTIONAL :: X_D_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X_D, X1, eps

    X_D = 0.5_DP
    IF( PRESENT( X_D_Option ) ) &
      X_D = X_D_Option

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2)') &
      '', 'X_D = ', X_D
    WRITE(*,*)
    WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
      '', 'D_L = ', D_L, '', 'D_R = ', D_R
    WRITE(*,'(A7,A6,3ES10.3E2,A4,A6,3ES10.3E2)') &
      '', 'V_L = ', V_L, '', 'V_R = ', V_R
    WRITE(*,'(A7,A6,ES10.3E2,A24,A6,ES10.3E2)') &
      '', 'P_L = ', P_L, '', 'P_R = ', P_R

    ! Loop over elements
    DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          ! Loop over all nodes in an element
          DO iNodeX = 1, nDOFX

            iNodeX1 = NodeNumberTableX(1,iNodeX) ! Particular node

            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

!!$            ! --- Riemann Problems ---
!!$            IF( X1 <= X_D )THEN
!!$
!!$              ! -- Left State --
!!$
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = D_L
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_L(1)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_L(2)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_L(3)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = P_L / (Gamma_IDEAL-One)
!!$              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = P_L
!!$
!!$            ELSE
!!$
!!$              ! -- Right State --
!!$
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = D_R
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V_R(1)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = V_R(2)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = V_R(3)
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = P_R / (Gamma_IDEAL-One)
!!$              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = P_R
!!$
!!$            END IF

!!$            ! --- Smooth solution, advected sine wave ---
!!$
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
!!$              = 1.0_DP + 0.1_DP * SIN( TwoPi * X1 )
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0d-6 / (Gamma_IDEAL-One)
!!$            uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-6

!!$            ! --- Contact discontinuity, top-hat ---
!!$
!!$            IF( ( X1 <= 0.25d0 ) .OR. ( X1 >= 0.75d0 ) )THEN
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0d0
!!$            ELSE
!!$              uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 2.0d0
!!$            END IF
!!$
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 1.0d-6 / (Gamma_IDEAL-One)
!!$            uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-6

            ! --- Perturbed shock tube ---
            IF( X1 <= X_D )THEN

              ! -- Left State --

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 5.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 50.0d0 / (Gamma_IDEAL-One)
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 50.0d0

            ELSE

              ! -- Right State --

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 2.0d0 + 0.3d0 * SIN( 50.0d0 * X1 )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0d0
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  = 5.0d0 / (Gamma_IDEAL-One)
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 5.0d0

            END IF


!!$            ! --- Acoustic problem, spherical symmetry ---
!!$
!!$            IF( ( X1 .GE. 0.4d0 ) .AND. ( X1 .LE. 0.6d0 ) )THEN
!!$              eps = 1.0d-4
!!$            ELSE
!!$              eps = 0.0d0
!!$            END IF
!!$
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
!!$              = 1.0_DP + eps * SIN( Pi * X1 )**4 / X1
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
!!$              = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
!!$              = 0.0_DP
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) &
!!$              = 0.0_DP
!!$            uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
!!$              = 1.0_DP / Gamma_IDEAL + eps * SIN( Pi * X1 )**4 / X1
!!$            uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
!!$              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - 1.0_DP )
            
          END DO

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33),                      &
                   uAF(:,iX1,iX2,iX3,iAF_P ) )
          
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblem


  SUBROUTINE InitializeFields_StandingAccretionShock

    REAL(DP) :: D, V(3), P
    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    INTEGER, PARAMETER :: i_r = 1, i_D = 2, i_V1 = 3, i_E = 4
    INTEGER  :: iL, nLines, iX1_Debug = 800
    REAL(DP) :: X1
    REAL(DP), ALLOCATABLE :: FluidFieldData(:,:), FluidFieldParameters(:)

    WRITE(*,*)
    WRITE(*,'(A2,A6,A)') &
      '', 'INFO: ', TRIM( ProgramName )

    CALL ReadParameters &
           ( '../StandingAccretionShock_Parameters.dat', FluidFieldParameters )
    CALL ReadData &
           ( '../StandingAccretionShock_Data.dat', nLines, FluidFieldData )

    ! --- Interpolate initial conditions onto grid ---
    WRITE(*,'(A)') 'Interpolating initial conditions onto grid'

    ! --- Loop over all elements ---
    DO iX3 = iX_B1(3), iX_E1(3)
      DO iX2 = iX_B1(2), iX_E1(2)
        DO iX1 = iX_B1(1), iX_E1(1)

          ! --- Loop over all nodes in an element ---
          DO iNodeX = 1, nDOFX

            ! --- Isolate node in X1 direction ---
            iNodeX1 = NodeNumberTableX(1,iNodeX)

            ! --- Physical coordinate corresponding to iNodeX1 ---
            X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

            ! --- Get lower index of input array
            !     (FluidFieldData) corresponding to physical coordinate (X1) ---
            iL = Locate( X1, FluidFieldData(:,i_r), nLines )

            IF( iX1 == iX1_Debug ) iL_Debug = iL

            ! --- Interpolate to the physical point X1 ---

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_D, i_r, iL, X1, FluidFieldData )

            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_V1, i_r, iL, X1, FluidFieldData )

            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_Ne) = Zero

            uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
              = InterpolateInitialConditionsOntoGrid &
                  ( i_E, i_r, iL, X1, FluidFieldData )

            ! --- Compute pressure from internal energy density ---
            uAF(iNodeX,iX1,iX2,iX3,iAF_P) &
              = ( Gamma_IDEAL - 1.0_DP ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

          END DO ! --- Loop over nodes ---

          CALL ComputeConserved_GR &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
                   uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
                   uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
                   uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
                   uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
                   uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22),                      &
                   uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33),                      &
                   uAF(:,iX1,iX2,iX3,iAF_P ) )

        END DO
      END DO
    END DO

!!$    WRITE(*,*) iL_Debug, iX1_Debug
!!$    WRITE(*,'(A,ES22.16E2)') 'rhoL: ', uPF(1,iX1_Debug,1,1,iPF_D )
!!$    WRITE(*,'(A,ES22.16E2)') 'rhoU: ', uPF(1,iX1_Debug,1,1,iPF_D )
!!$    WRITE(*,'(A)') 'Stopping...'
!!$    STOP

  END SUBROUTINE InitializeFields_StandingAccretionShock


  REAL(DP) FUNCTION InterpolateInitialConditionsOntoGrid &
                      (iVar, i_r, iL, X, FluidFieldData) RESULT( yInterp )

    INTEGER,  INTENT(in) :: iL, iVar, i_r
    REAL(DP), INTENT(in) :: X
    REAL(DP), INTENT(in) :: FluidFieldData(:,:)
    REAL(DP)             :: m, X1, X2, Y1, Y2

    X1 = FluidFieldData(iL,i_r)
    X2 = FLuidFieldData(iL+1,i_r)
    Y1 = FluidFieldData(iL,iVar)
    Y2 = FluidFieldData(iL+1,iVar)

    m = ( Y2 - Y1 ) / ( X2 - X1 )

    ! --- Using only lower limit for slope ---
    yInterp = m * ( X - X1 ) + Y1

    ! --- Using average slope ---
    ! --- Only changes accuracy in 12th decimal place ---
!    yInterp = ( 2.0_DP * m * ( X - X1 ) * ( X2 - X ) + ( Y1 * X2 + Y2 * X1 ) &
!                - X * ( Y1 + Y2 ) ) / ( X1 + X2 - 2.0_DP * X )

    IF( iL == iL_Debug  .AND. iVar == iPF_D ) THEN
      WRITE(*,'(A,ES24.16E3)') 'X1: ', X1
      WRITE(*,'(A,ES24.16E3)') 'X2: ', X2
      WRITE(*,'(A,ES24.16E3)') 'Y1: ', Y1
      WRITE(*,'(A,ES24.16E3)') 'Y2: ', Y2
      WRITE(*,'(A,ES24.16E3)') 'm:  ', m
      WRITE(*,'(A,ES24.16E3)') 'yInterp: ', yInterp
    END IF

    RETURN
  END FUNCTION InterpolateInitialConditionsOntoGrid


END MODULE InitializationModule_GR
