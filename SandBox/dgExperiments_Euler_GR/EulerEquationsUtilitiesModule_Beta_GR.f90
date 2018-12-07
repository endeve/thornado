MODULE EulerEquationsUtilitiesModule_Beta_GR

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One, Two, Three, Four
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX, nX
  USE MeshModule, ONLY: &
    MeshX
  USE ReferenceElementModuleX, ONLY:          &
    nDOFX_X1, nDOFX_X2, nDOFX_X3, WeightsX_q, &
    WeightsX_X1, WeightsX_X2, WeightsX_X3,    &
    WeightsLX1,  WeightsLX2,  WeightsLX3
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    dLXdX1_q, dLXdX2_q, dLXdX3_q,             &
    LX_X1_Dn, LX_X1_Up,                       &
    LX_X2_Dn, LX_X2_Up,                       &
    LX_X3_Dn, LX_X3_Up
  USE FluidFieldsModule, ONLY:                         &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_P, iAF_Cs
  USE GeometryFieldsModule, ONLY:             &
    iGF_h_1, iGF_h_2, iGF_h_3,                &
    iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, iGF_Beta_2, iGF_Beta_3, nGF
  USE EquationOfStateModule, ONLY:             &
    ComputePressureFromSpecificInternalEnergy, &
    ComputeSoundSpeedFromPrimitive_GR
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL

  IMPLICIT NONE

  PRIVATE :: ComputeFunP

  ! --- ComputePressureWithB*Method is only public for debugging ---
  PUBLIC :: ComputePressureWithBisectionMethod
  PUBLIC :: ComputePressureWithBrentsMethod

  PUBLIC :: ComputeFromConserved_GR
  PUBLIC :: ComputePrimitive_GR
  PUBLIC :: ComputeFunJacP
  PUBLIC :: ComputeConserved_GR
  PUBLIC :: ComputeTimeStep_GR
  PUBLIC :: Eigenvalues_GR
  PUBLIC :: Flux_X1_GR
  PUBLIC :: Flux_X2_GR
  PUBLIC :: StressTensor_Diagonal
  PUBLIC :: AlphaC_GR
  PUBLIC :: NumericalFlux_LLF_GR
  PUBLIC :: NumericalFlux_HLL_GR
  PUBLIC :: NumericalFlux_X1_HLLC_GR
  PUBLIC :: NumericalFlux_X2_HLLC_GR

  REAL(DP), PARAMETER :: TolP = 1.0d-8, TolFunP = 1.0d-6, MachineEPS = 1.0d-16
  LOGICAL :: DEBUG          = .FALSE.
  LOGICAL :: DEBUG_bisec    = .FALSE.
  LOGICAL :: DEBUG_FunJacP  = .FALSE.
  LOGICAL :: DEBUG_TimeStep = .FALSE.

CONTAINS


  SUBROUTINE ComputeFromConserved_GR( iX_B0, iX_E0, G, U, P, A )

    INTEGER, INTENT(in)  :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:), &
      U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(inout)  :: &
      P(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:), &
      A(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    INTEGER :: iX1, iX2, iX3

    ! --- Update primitive variables, pressure, and sound speed
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      CALL ComputePrimitive_GR &
             ( U(1:,iX1,iX2,iX3,iCF_D),        &
               U(1:,iX1,iX2,iX3,iCF_S1),       &
               U(1:,iX1,iX2,iX3,iCF_S2),       &
               U(1:,iX1,iX2,iX3,iCF_S3),       &
               U(1:,iX1,iX2,iX3,iCF_E),        &
               U(1:,iX1,iX2,iX3,iCF_Ne),       &
               P(1:,iX1,iX2,iX3,iPF_D),        &
               P(1:,iX1,iX2,iX3,iPF_V1),       &
               P(1:,iX1,iX2,iX3,iPF_V2),       &
               P(1:,iX1,iX2,iX3,iPF_V3),       &
               P(1:,iX1,iX2,iX3,iPF_E),        &
               P(1:,iX1,iX2,iX3,iPF_Ne),       &
               A(1:,iX1,iX2,iX3,iAF_P),        &
               G(1:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(1:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(1:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeSoundSpeedFromPrimitive_GR &
             ( P(1:,iX1,iX2,iX3,iPF_D), P(1:,iX1,iX2,iX3,iPF_E), &
                 P(1:,iX1,iX2,iX3,iPF_Ne), A(1:,iX1,iX2,iX3,iAF_Cs) )

    END DO
    END DO
    END DO

  END SUBROUTINE ComputeFromConserved_GR


  SUBROUTINE ComputePrimitive_GR &
              ( CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
                PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
                AF_P, GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33 )

    REAL(DP), DIMENSION(:), INTENT(in)  :: CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, AF_P
    REAL(DP), DIMENSION(:), INTENT(in)  :: GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33

    LOGICAL :: CONVERGED

    INTEGER :: i, ITERATION, nNodes
    INTEGER, PARAMETER :: MAX_IT = 100

    REAL(DP) :: SSq, Pold, vSq, W, h, Pnew, q, Pbisec

    REAL(DP) :: FunP, JacP, FunP1

    IF( DEBUG )THEN
      WRITE(*,*)
      WRITE(*,'(A)') '    ================================'
      WRITE(*,'(A)') '    EE: Entering ComputePrimitive_GR'
      WRITE(*,'(A)') '    ================================'
      WRITE(*,*)
    END IF

    ! --- Loop through all the nodes ---
    nNodes = SIZE( CF_D )
    DO i = 1, nNodes

      IF( DEBUG ) WRITE(*,'(A,I1)') '    EE: iNode = ', i
      IF( DEBUG ) WRITE(*,'(A)')    '    -------------'
      IF( DEBUG ) WRITE(*,'(A,ES24.16E3)')    '    EE: S1: ', CF_S1(i)
      IF( DEBUG ) WRITE(*,'(A,ES24.16E3)')    '    EE: S2: ', CF_S2(i)
      IF( DEBUG ) WRITE(*,'(A,ES24.16E3)')    '    EE: S3: ', CF_S3(i)
      IF( DEBUG ) WRITE(*,*)

      q = CF_E(i) + CF_D(i) - SQRT( CF_D(i)**2 &
                                  + CF_S1(i)**2 / GF_Gm_dd_11(i)  &
                                  + CF_S2(i)**2 / GF_Gm_dd_22(i)  &
                                  + CF_S3(i)**2 / GF_Gm_dd_33(i) )

      IF( q .LT. Zero )THEN
        WRITE(*,*)
        WRITE(*,'(A6,I1)')        'Node:    ', i
        WRITE(*,'(A9,ES18.10E3)') 'q:       ', q
        WRITE(*,'(A9,ES18.10E3)') 'Gm11(i): ', GF_Gm_dd_11(i)
        WRITE(*,'(A9,ES18.10E3)') 'Gm22(i): ', GF_Gm_dd_22(i)
        WRITE(*,'(A9,ES18.10E3)') 'Gm33(i): ', GF_Gm_dd_33(i)
        WRITE(*,'(A9,ES18.10E3)') 'D(i):    ', CF_D(i)
        WRITE(*,'(A9,ES18.10E3)') 'tau(i):  ', CF_E(i)
        WRITE(*,'(A9,ES18.10E3)') 'S1(i):   ', CF_S1(i)
        WRITE(*,'(A9,ES18.10E3)') 'S2(i):   ', CF_S2(i)
        WRITE(*,'(A9,ES18.10E3)') 'S3(i):   ', CF_S3(i)
        STOP 'q < 0'
      END IF
    
      SSq =   CF_S1(i)**2 / GF_Gm_dd_11(i) &
            + CF_S2(i)**2 / GF_Gm_dd_22(i) &
            + CF_S3(i)**2 / GF_Gm_dd_33(i)

      ! --- Find Pressure with Newton's Method ---

      ! --- Get initial guess for pressure with bisection method ---
      IF( DEBUG ) &
        WRITE(*,'(A)') '    EE: CALL ComputePressureWithBisectionMethod'
      CALL ComputePressureWithBisectionMethod( CF_D(i), CF_E(i), SSq, Pbisec )
      Pold = Pbisec
      IF( Pbisec .LT. SqrtTiny )THEN
        WRITE(*,'(A)') '    EE: Pbisec < SqrTiny. Stopping...'
        STOP
      END IF
     
      ! --- Get initial value for FunP ---
      IF( DEBUG ) &
        WRITE(*,'(A)') '    EE: CALL ComputeFunJacP for Pbisec'
      CALL ComputeFunJacP( CF_D(i), CF_E(i), SSq, Pold, FunP, JacP )
      IF( DEBUG ) &
        WRITE(*,'(A,ES24.16E3)') '    EE: FunPbisec / Pbisec: ', FunP / Pbisec

      CONVERGED = .FALSE.

      IF( ABS( FunP ) .LT. TolFunP * ABS( TolFunP ) )THEN
        Pnew = Pold
        CONVERGED = .TRUE. 
      ELSE
        IF( DEBUG )THEN
          WRITE(*,*)
          WRITE(*,'(A)') "    EE: Start finding root with Newton's method"
          WRITE(*,'(A)') "    -------------------------------------------"
          WRITE(*,*)
        END IF
      END IF

      ITERATION = 0

      DO WHILE ( .NOT. CONVERGED )

        ITERATION = ITERATION + 1

        IF( DEBUG ) WRITE(*,'(A)') '    EE: CALL ComputeFunJacP'
        CALL ComputeFunJacP( CF_D(i), CF_E(i), SSq, Pold, FunP, JacP )

        IF( ABS( FunP / JacP ) .LT. MachineEPS ) THEN
          IF( DEBUG )THEN
            WRITE(*,'(A)')           '    EE: ABS(FunP/JacP) < MachineEps'
            WRITE(*,'(A,I3)')        '    EE: ITERATION = ', ITERATION
            WRITE(*,'(A,ES24.16E3)') '    EE: Pold      = ', Pold
            WRITE(*,'(A,ES24.16E3)') '    EE: Pnew      = ', Pnew
            WRITE(*,'(A,ES24.16E3)') '    EE: Pbisec    = ', Pbisec
            WRITE(*,*)
          END IF
          Pnew = Pold
          CONVERGED = .TRUE.
          EXIT
        END IF

        Pnew = MAX( Pold - FunP / JacP, SqrtTiny )

        ! --- Check if Newton's method has converged ---
        IF( ( ABS( ( Pnew - Pold ) / Pold ) ) .LT. TolP )THEN

          CALL ComputeFunJacP( CF_D(i), CF_E(i), SSq, Pnew, FunP, JacP )

          IF( .NOT. ( ABS( FunP ) .LT. TolFunP * ABS( Pnew ) ) ) THEN
            WRITE(*,*)
            WRITE(*,'(A)') 'WARNING:'
            WRITE(*,'(A)') '--------'
            WRITE(*,'(A)') 'ABS( FunP / Pnew ) > TolFunP'
            WRITE(*,'(A,ES24.16E3)') 'TolFunP              = ', TolFunP
            WRITE(*,'(A,ES24.16E3)') 'Pold                 = ', Pold
            WRITE(*,'(A,ES24.16E3)') 'Pnew                 = ', Pnew
            WRITE(*,'(A,ES24.16E3)') &
              '|(Pnew - Pold)/Pold| = ', ABS( ( Pnew - Pold ) / Pold )
            WRITE(*,'(A,ES24.16E3)') &
              '|FunP(Pnew) / Pnew|  = ', ABS( FunP / Pnew )
          END IF

          CONVERGED = .TRUE.
          EXIT

        END IF

        ! --- STOP after MAX_IT iterations ---
        IF( ITERATION .GE. MAX_IT - 4 .OR. DEBUG )THEN

          IF( ITERATION .EQ. MAX_IT - 4 )THEN
            WRITE(*,*)
            WRITE(*,'(A)')          '    EE: Max iterations IF statement'
            WRITE(*,'(A)')          '    -------------------------------'
            WRITE(*,'(A,I1)')       '    EE: Node: ', i
            WRITE(*,'(A,ES7.1E2)')  '    EE: TolP    = ', TolP
            WRITE(*,'(A,ES7.1E2)')  '    EE: TolFunP = ', TolFunP
            WRITE(*,*)
            WRITE(*,'(A,ES24.16E3)') '  U(i,iCF_D)  = ', CF_D(i)
            WRITE(*,'(A,ES24.16E3)') '  U(i,iCF_S1) = ', CF_S1(i)
            WRITE(*,'(A,ES24.16E3)') '  U(i,iCF_S2) = ', CF_S2(i)
            WRITE(*,'(A,ES24.16E3)') '  U(i,iCF_S3) = ', CF_S3(i)
            WRITE(*,'(A,ES24.16E3)') '  U(i,iCF_E)  = ', CF_E(i)
            WRITE(*,*)
            WRITE(*,'(A,ES24.16E3)') '  G(i,iGF_Gm_dd_11) = ', GF_Gm_dd_11(i)
            WRITE(*,'(A,ES24.16E3)') '  G(i,iGF_Gm_dd_22) = ', GF_Gm_dd_22(i)
            WRITE(*,'(A,ES24.16E3)') '  G(i,iGF_Gm_dd_33) = ', GF_Gm_dd_33(i)
            WRITE(*,*)
            WRITE(*,'(A,ES24.16E3)') '    EE: |v|   = ', &
              ABS( CF_S1(i) / ( CF_E(i) + CF_D(i) + Pnew ) )
            WRITE(*,'(A,ES24.16E3)') '    EE: W     = ', &
              1.0d0 / SQRT( 1.0d0 - SSq / ( CF_E(i) + CF_D(i) + Pnew )**2 )
            WRITE(*,'(A,ES24.16E3)') '    EE: P/rho = ', &
              Pnew / ( CF_D(i) / ( 1.0d0 / SQRT( 1.0d0 &
                - SSq / ( CF_E(i) + CF_D(i) + Pnew )**2 ) ) )
            WRITE(*,*)
            WRITE(*,'(A,ES24.16E3)') '    EE: Pold   = ', Pold
            WRITE(*,*)
          END IF

          WRITE(*,'(A,I3)') '    EE: ITERATION: ', ITERATION
          WRITE(*,'(A)')    '    EE: --------------'
          WRITE(*,'(A,ES24.16E3)') '    EE: Pold      = ', Pold
          WRITE(*,'(A,ES24.16E3)') '    EE: Pnew      = ', Pnew
          WRITE(*,'(A,ES24.16E3)') '    EE: FunP      = ', FunP
          WRITE(*,'(A,ES24.16E3)') '    EE: JacP      = ', JacP
          WRITE(*,'(A,ES24.16E3)') &
            '    EE: dP/|Pold| = ', ( Pnew - Pold ) / ABS( Pold )
          WRITE(*,'(A,ES24.16E3)') '    EE: FunP/Pold = ', FunP / Pold
          WRITE(*,*)
          IF( ITERATION .EQ. MAX_IT )THEN
            STOP 'Max allowed iterations reached, no convergence.'
          END IF
        END IF

        Pold = Pnew

      END DO

      AF_P(i) = Pnew

      IF( DEBUG )THEN
        WRITE(*,'(A,ES24.16E3)') '    EE: Pfinal: ', Pnew
        WRITE(*,*)
      END IF

      vSq = SSq / ( CF_E(i) + Pnew + CF_D(i) )**2

      W = 1.0_DP / SQRT( 1.0_DP - vSq )

      h = ( CF_E(i) + Pnew + CF_D(i) ) / ( W * CF_D(i) )

      IF( DEBUG ) THEN
        IF( h .LT. 1.0_DP ) THEN
          WRITE(*,'(A,ES18.10E3)') '    EE: h:    ', h
        END IF
      END IF

      ! --- Recover Primitive Variables ---

      PF_D(i)  = CF_D(i) / W
      
      PF_V1(i) = CF_S1(i) / ( CF_D(i) * h * W * GF_Gm_dd_11(i) )

      PF_V2(i) = CF_S2(i) / ( CF_D(i) * h * W * GF_Gm_dd_22(i) )

      PF_V3(i) = CF_S3(i) / ( CF_D(i) * h * W * GF_Gm_dd_33(i) )

      PF_E(i)  = CF_D(i) * ( h - 1.0_DP ) / W - Pnew
      
      PF_Ne(i) = CF_Ne(i) / W


   END DO ! --- End of loop over nodes ---


  END SUBROUTINE ComputePrimitive_GR


  SUBROUTINE ComputeFunJacP( CF_D, CF_E, SSq, P, FunP, JacP, Pnorm_Option )

    REAL(DP), INTENT(in)  :: CF_D, CF_E, SSq, P
    REAL(DP), INTENT(in), OPTIONAL :: Pnorm_Option
    REAL(DP), INTENT(out) :: FunP, JacP

    REAL(DP) :: HSq, RHO, EPS, dRHO, dEPS, Pnorm
    REAL(DP), DIMENSION(1) :: Pbar

    Pnorm = One
    IF( PRESENT( Pnorm_Option ) ) Pnorm = Pnorm_Option
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: Pnorm:   ', Pnorm

    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: CF_D:    ', CF_D
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: CF_E:    ', CF_E
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: SSq:     ', SSq
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: P:       ', P

    HSq = ( CF_E + P + CF_D )**2
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: HSq:     ', HSq
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') &
      '      CFJP: H^2-S^2: ', HSq - SSq

    RHO = CF_D * SQRT( HSq - SSq ) / SQRT( HSq )
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: RHO:     ', RHO

    EPS = ( SQRT( HSq - SSq ) &
            - P * SQRT( HSq ) / SQRT( HSq - SSq ) - CF_D ) / CF_D
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: EPS:     ', EPS

    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') &
      '      CFJP: CALL ComputePressureFromSpecificInternalEnergy'
    CALL ComputePressureFromSpecificInternalEnergy &
         ( [ RHO ], [ EPS ], [ 0.0_DP ], Pbar )

    FunP = ( P - Pbar(1) ) / Pnorm
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: FunP:     ', FunP

    dRHO = CF_D * SSq / ( SQRT( HSq - SSq ) * HSq )
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: dRHO:    ', dRHO

    dEPS = P * SSq / ( ( HSq - SSq ) * SQRT( HSq ) * RHO )
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: dEPS:     ', dEPS

    JacP = ( 1.0_DP - Pbar(1) * ( dRHO / RHO + dEPS / EPS ) ) / Pnorm
    IF( DEBUG_FunJacP ) WRITE(*,'(A,ES24.16E3)') '      CFJP: JacP:     ', JacP
    IF( DEBUG_FunJacP ) WRITE(*,*)

  END SUBROUTINE ComputeFunJacP


  SUBROUTINE ComputePressureWithBisectionMethod( CF_D, CF_E, SSq, P )

    REAL(DP), INTENT(in)  :: CF_D, CF_E, SSq
    REAL(DP), INTENT(out) :: P

    INTEGER,  PARAMETER :: MAX_IT = 10

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: PA, PB, PC, DeltaP
    REAL(DP) :: FunPA, FunPB, FunPC

    ! --- Get upper and lower bounds on pressure, PA, PB ---
    PA = SqrtTiny
    PB = Three * CF_E
!    PB = Two * ( One - One / Gamma_IDEAL ) * CF_E
    IF( DEBUG_Bisec ) WRITE(*,'(A,ES24.16E3)') '      CPWBM: PA: ', PA
    IF( DEBUG_Bisec ) WRITE(*,'(A,ES24.16E3)') '      CPWBM: PB: ', PB

    ! --- Compute FunP for upper and lower bounds ---
    IF( DEBUG_Bisec ) WRITE(*,'(A,ES24.16E3)') '      CPWBM: CALL ComputeFunPA'
    CALL ComputeFunP( CF_D, CF_E, SSq, PA, FunPA )
    IF( DEBUG_Bisec ) WRITE(*,'(A,ES24.16E3)') '      CPWBM: CALL ComputeFunPB'
    CALL ComputeFunP( CF_D, CF_E, SSq, PB, FunPB )

    ! --- Check that sign of FunP changes across bounds ---
    IF( .NOT. FunPA * FunPB .LT. 0 )THEN

      WRITE(*,'(A6,A)') &
        '', 'ComputePressureWithBisectionMethod:'
      WRITE(*,'(A8,A)') &
        '', 'Error: No Root in Interval'
      WRITE(*,'(A8,A,2ES15.6E3)') &
        '', 'PA, PB = ', PA, PB
      WRITE(*,'(A8,A,2ES15.6E3)') &
        '', 'FunPA, FunPB = ', FunPA, FunPB
      STOP

    END IF

    DeltaP = PB - PA

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( ITERATION .LT. MAX_IT )

      ITERATION = ITERATION + 1

      ! --- Compute midpoint, PC ---
      DeltaP = Half * DeltaP
      PC     = PA + DeltaP

      ! --- Compute FunP for midpoint pressure ---
      CALL ComputeFunP( CF_D, CF_E, SSq, PC, FunPC )

      ! --- Change PC to PA or PB, depending on sign of FunP(PC) ---
      IF( FunPA * FunPC .LT. Zero )THEN

        PB = PC
        FunPB = FunPC

      ELSE

        PA = PC
        FunPA = FunPC

      END IF

    END DO

    P = PA

    IF( DEBUG_bisec ) WRITE(*,'(A,ES24.16E3)') '      CPWBM: Pbisec: ', P
    IF( DEBUG_bisec ) WRITE(*,*)
    
  END SUBROUTINE ComputePressureWithBisectionMethod

  
  SUBROUTINE ComputePressureWithBrentsMethod( CF_D, CF_E, SSq, Pbrent )

    REAL(DP), INTENT(in)  :: CF_D, CF_E, SSq
    REAL(DP), INTENT(out) :: Pbrent

    REAL(DP) :: PA, PB, PC, PD, PS, PSwap, Pbisec
    REAL(DP) :: FunPA, FunPB, FunPC, FunPD, FunPS
    REAL(DP) :: JacPA, JacPB, JacPC, JacPD, JacPS

    LOGICAL :: mflag, COND1, COND2, COND3, COND4, COND5

    INTEGER, PARAMETER :: MAX_IT = 100
    INTEGER :: ITERATION

    ! --- Get Pbisec to normalize FunP ---
    CALL ComputePressureWithBisectionMethod( CF_D, CF_E, SSq, Pbisec )

    ! --- Define values that bracket the root ---
    PA = Zero
    PB = Three * CF_E
!    PB = Two * ( One - One / Gamma_IDEAL ) * CF_E
    CALL ComputeFunJacP( CF_D, CF_E, SSq, PA, FunPA, JacPA )
    CALL ComputeFunJacP( CF_D, CF_E, SSq, PB, FunPB, JacPB )

    ! --- Exit if solution is not bracketed ---
    IF( .NOT. FunPA * FunPB .LT. 0.0_DP )THEN
      WRITE(*,'(A)') 'No root in interval. Exiting...'
      STOP
    END IF

    ! --- Conditionally swap PA and PB ---
    IF( ABS( FunPA ) .LT. ABS( FunPB ) )THEN
       PSwap = PB
       PB    = PA
       PA    = PSwap
       CALL ComputeFunJacP( CF_D, CF_E, SSq, PA, FunPA, JacPA )
       CALL ComputeFunJacP( CF_D, CF_E, SSq, PB, FunPB, JacPB )
    END IF

    ! --- Set PC = PA and compute FunPC ---
    PC = PA
    CALL ComputeFunJacP( CF_D, CF_E, SSq, PC, FunPC, JacPC )

    ! --- Set mflag ---
    mflag = .TRUE.

    ITERATION = 0
    DO WHILE(       ( ABS( PB - PA ) / ABS( Pbisec ) .GE. TolP ) &
              .AND. ( ABS(FunPB) .GE. TolFunP ) &
              .AND. ( ITERATION .LE. MAX_IT ) )

      ITERATION = ITERATION + 1

      IF( ( FunPA .NE. FunPC ) .AND. ( FunPB .NE. FunPC ) )THEN
        ! --- Inverse quadratic interpolation ---  
        PS =   PA * FunPB * FunPC / ( ( FunPA - FunPB ) * ( FunPA - FunPC ) ) &
             + PB * FunPA * FunPC / ( ( FunPB - FunPA ) * ( FunPB - FunPC ) ) &
             + PC * FunPA * FunPB / ( ( FunPC - FunPA ) * ( FunPC - FunPB ) )
      ELSE
        ! --- Secant method (replace with Newton's method?) ---
        PS = PB - FunPB * ( PB - PA ) / ( FunPB - FunPA )
      END IF

      ! --- Condition 1 ---
      IF( PB .GT. ( 3.0_DP * PA + PB ) / 4.0_DP )THEN
        IF( .NOT. ( ( PS .GT. ( 3.0_DP * PA + PB ) / 4.0_DP ) &
                    .AND. ( PS .LT. PB ) ) )THEN
          COND1 = .TRUE.
        ELSE
          COND1 = .FALSE.
        END IF
      ELSE
        IF( .NOT. ( ( PS .LT. ( 3.0_DP * PA + PB ) / 4.0_DP ) &
                    .AND. ( PS .GT. PB ) ) )THEN
          COND1 = .TRUE.
        ELSE
          COND1 = .FALSE.
        END IF
      END IF

      ! --- Condition 2 ---
      IF( mflag .AND. ( ABS( PS - PB ) &
          .GE. ( ABS( PB - PC ) / 2.0_DP ) ) )THEN
        COND2 = .TRUE.
      ELSE
        COND2 = .FALSE.
      END IF

      ! --- Condition 3 ---
      IF( .NOT. mflag .AND. ( ABS( PS - PB ) &
          .GE. ( ABS( PC - PD ) / 2.0_DP ) ) )THEN
        COND3 = .TRUE.
      ELSE
        COND3 = .FALSE.
      END IF

      ! --- Condition 4 ---
      IF( mflag .AND. ( ABS( PB - PC ) .LT. TolP ) )THEN
        COND4 = .TRUE.
      ELSE
        COND4 = .FALSE.
      END IF

      ! --- Condition 5 ---
      IF( .NOT. mflag .AND. ( ABS( PC - PD ) .LT. TolP ) )THEN
        COND5 = .TRUE.
      ELSE
        COND5 = .FALSE.
      END IF

      IF( COND1 .OR. COND2 .OR. COND3 .OR. COND4 .OR. COND5 )THEN
        PS = ( PA + PB ) / 2.0_DP
        mflag = .TRUE.
      ELSE
        mflag = .FALSE.
      END IF

      CALL ComputeFunJacP( CF_D, CF_E, SSq, PS, FunPS, JacPS )
      PD = PC
      PC = PB

      IF( FunPA * FunPS .LT. 0.0_DP )THEN
        PB = PS
        CALL ComputeFunJacP( CF_D, CF_E, SSq, PB, FunPB, JacPB )
      ELSE
        PA = PS
        CALL ComputeFunJacP( CF_D, CF_E, SSq, PA, FunPA, JacPA )
      END IF

      ! --- Conditionally swap PA and PB ---
      IF( ABS( FunPA ) .LT. ABS( FunPB ) )THEN
        PSwap = PB
        PB    = PA
        PA    = PSwap
        CALL ComputeFunJacP( CF_D, CF_E, SSq, PA, FunPA, JacPA )
        CALL ComputeFunJacP( CF_D, CF_E, SSq, PB, FunPB, JacPB )
      END IF

      Pbrent = PB

      IF( ITERATION .EQ. MAX_IT )THEN
        WRITE(*,*)
        WRITE(*,'(A)') 'Maximum number of iterations reached.'
        WRITE(*,'(A,ES24.16E3)') 'Pbrent =        ', Pbrent
        WRITE(*,'(A,ES24.16E3)') '|dP|/|Pbisec| = ', &
          ABS( PB - PA ) / ABS( Pbisec )
        WRITE(*,'(A,ES24.16E3)') '|Pbrent-Pbisec|/|Pbisec| = ', &
          ABS( Pbrent - Pbisec ) / ABS( Pbisec )
        STOP
      END IF

    END DO

    RETURN
  END SUBROUTINE ComputePressureWithBrentsMethod


  SUBROUTINE ComputeFunP( CF_D, CF_E, SSq, P, FunP )

    REAL(DP), INTENT(in)  :: CF_D, CF_E, SSq, P
    REAL(DP), INTENT(out) :: FunP

    REAL(DP) :: HSq, RHO, EPS
    REAL(DP), DIMENSION(1) :: Pbar

    HSq = ( CF_E + P + CF_D )**2

    RHO = CF_D * SQRT( HSq - SSq ) / SQRT( HSq )

    EPS = ( SQRT( HSq - SSq ) &
            - P * SQRT( HSq ) / SQRT( HSq - SSq ) - CF_D ) / CF_D

    CALL ComputePressureFromSpecificInternalEnergy &
         ( [ RHO ], [ EPS ], [ 0.0_DP ], Pbar )

    FunP = P - Pbar(1)

  END SUBROUTINE ComputeFunP

  SUBROUTINE ComputeTimeStep_GR &
               ( iX_B0, iX_E0, G, U, CFL, TimeStep, UseSourceTerm_Option )

    INTEGER,  INTENT(in)  :: &
      iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:), &
      U(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep
    LOGICAL, INTENT(in), OPTIONAL :: UseSourceTerm_Option

    INTEGER  :: iX1, iX2, iX3, iNodeX, iDimX, iGF
    REAL(DP) :: dX(3), dt_X(3), dt_S(3)
    REAL(DP) :: P(nDOFX,nPF)
    REAL(DP) :: A(nDOFX,nAF)
    REAL(DP) :: PressureTensor(3)
    REAL(DP) :: EigVals_X1(nCF,nDOFX), Max_X1, &
                EigVals_X2(nCF,nDOFX), Max_X2, &
                EigVals_X3(nCF,nDOFX), Max_X3
    REAL(DP) :: SourceTerm(3), PosRoot(3), NegRoot(3), a2, a1, a0
    REAL(DP) :: alpha, W, h
    REAL(DP) :: tau
    REAL(DP) :: epsilon = Half
    REAL(DP), DIMENSION(nDOFX) :: dh1dX1, dh2dX1, dh3dX1, &
                                  dh1dX2, dh2dX2, dh3dX2, &
                                  dh1dX3, dh2dX3, dh3dX3
    REAL(DP) :: G_X1_Dn(nDOFX_X1,nGF), G_X1_Up(nDOFX_X1,nGF), &
                G_X2_Dn(nDOFX_X2,nGF), G_X2_Up(nDOFX_X2,nGF), &
                G_X3_Dn(nDOFX_X3,nGF), G_X3_Up(nDOFX_X3,nGF)
    LOGICAL :: UseSourceTerm

    UseSourceTerm = .TRUE.
    IF( PRESENT( UseSourceTerm_Option) ) &
      UseSourceTerm = UseSourceTerm_Option

    tau = ( Gamma_IDEAL - 1.0d0 ) / Gamma_IDEAL

    TimeStep = HUGE( One )
    dt_X(:)  = HUGE( One )
    dt_S(:)  = HUGE( One )

    PosRoot = HUGE( One )
    NegRoot = HUGE( One )

    Max_X1 = -Huge( One )
    Max_X2 = -Huge( One )
    Max_X3 = -Huge( One )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      dX(1) = MeshX(1) % Width(iX1)
      dX(2) = MeshX(2) % Width(iX2)
      dX(3) = MeshX(3) % Width(iX3)

      CALL ComputePrimitive_GR &
             ( U(:,iX1,iX2,iX3,iCF_D ), U(:,iX1,iX2,iX3,iCF_S1), &
               U(:,iX1,iX2,iX3,iCF_S2), U(:,iX1,iX2,iX3,iCF_S3), &
               U(:,iX1,iX2,iX3,iCF_E ), U(:,iX1,iX2,iX3,iCF_Ne), &
               P(:,iPF_D), P(:,iPF_V1), P(:,iPF_V2), P(:,iPF_V3), &
               P(:,iPF_E), P(:,iPF_Ne), A(:,iAF_P), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )

      CALL ComputeSoundSpeedFromPrimitive_GR &
             ( P(:,iPF_D), P(:,iPF_E), P(:,iPF_Ne), A(:,iAF_Cs) )

      IF( UseSourceTerm )THEN

        ! --- Interpolate scale factors to faces ---
        DO iGF = iGF_h_1, iGF_h_3

          ! --- X1 ---
          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                   G(:,iX1-1,iX2,iX3,iGF), 1, Zero, G_X1_Dn(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                   G(:,iX1,  iX2,iX3,iGF), 1, Half, G_X1_Dn(:,iGF), 1 )

          G_X1_Dn(1:nDOFX_X1,iGF) &
            = MAX( G_X1_Dn(1:nDOFX_X1,iGF), SqrtTiny )

          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, One,  LX_X1_Up, nDOFX_X1, &
                   G(:,iX1,  iX2,iX3,iGF), 1, Zero, G_X1_Up(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X1, nDOFX, Half, LX_X1_Dn, nDOFX_X1, &
                   G(:,iX1+1,iX2,iX3,iGF), 1, Half, G_X1_Up(:,iGF), 1 )

          G_X1_Up(1:nDOFX_X1,iGF) &
            = MAX( G_X1_Up(1:nDOFX_X1,iGF), SqrtTiny )

          ! --- X2 ---
          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G(:,iX1,iX2-1,iX3,iGF), 1, Zero, G_X2_Dn(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G(:,iX1,iX2,  iX3,iGF), 1, Half, G_X2_Dn(:,iGF), 1 )

          G_X2_Dn(1:nDOFX_X2,iGF) &
            = MAX( G_X2_Dn(1:nDOFX_X2,iGF), SqrtTiny )

          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, One,  LX_X2_Up, nDOFX_X2, &
                   G(:,iX1,iX2  ,iX3,iGF), 1, Zero, G_X2_Up(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X2, nDOFX, Half, LX_X2_Dn, nDOFX_X2, &
                   G(:,iX1,iX2+1,iX3,iGF), 1, Half, G_X2_Up(:,iGF), 1 )

          G_X2_Up(1:nDOFX_X2,iGF) &
            = MAX( G_X2_Up(1:nDOFX_X2,iGF), SqrtTiny )

          ! --- X3 ---
          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                   G(:,iX1,iX2,iX3-1,iGF), 1, Zero, G_X3_Dn(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                   G(:,iX1,iX2,iX3,  iGF), 1, Half, G_X3_Dn(:,iGF), 1 )

          G_X3_Dn(1:nDOFX_X3,iGF) &
            = MAX( G_X3_Dn(1:nDOFX_X3,iGF), SqrtTiny )

          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, One,  LX_X3_Up, nDOFX_X3, &
                   G(:,iX1,iX2,iX3,  iGF), 1, Zero, G_X3_Up(:,iGF), 1 )
          CALL DGEMV &
                 ( 'N', nDOFX_X3, nDOFX, Half, LX_X3_Dn, nDOFX_X3, &
                   G(:,iX1,iX2,iX3+1,iGF), 1, Half, G_X3_Up(:,iGF), 1 )

          G_X3_Up(1:nDOFX_X3,iGF) &
            = MAX( G_X3_Up(1:nDOFX_X3,iGF), SqrtTiny )

        END DO ! Loop over scale factors

        ! --- Scale factor derivatives with respect to X1 ---
        ! --- dh1dX1 ---

        CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                   WeightsX_X1(:) * G_X1_Up(:,iGF_h_1), 1, Zero, dh1dX1, 1 )
        CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                   WeightsX_X1(:) * G_X1_Dn(:,iGF_h_1), 1,  One, dh1dX1, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_1), 1, &
                    One, dh1dX1, 1 )

        dh1dX1 = dh1dX1 / ( WeightsX_q(:) * dX(1) )

        ! --- dh2dX1 ---

        CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                   WeightsX_X1(:) * G_X1_Up(:,iGF_h_2), 1, Zero, dh2dX1, 1 )
        CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                   WeightsX_X1(:) * G_X1_Dn(:,iGF_h_2), 1,  One, dh2dX1, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_2), &
                    1,  One, dh2dX1, 1 )

        dh2dX1 = dh2dX1 / ( WeightsX_q(:) * dX(1) )

        ! --- dh3dX1 ---

        CALL DGEMV( 'T', nDOFX_X1, nDOFX, + One, LX_X1_Up, nDOFX_X1, &
                   WeightsX_X1(:) * G_X1_Up(:,iGF_h_3), 1, Zero, dh3dX1, 1 )
        CALL DGEMV( 'T', nDOFX_X1, nDOFX, - One, LX_X1_Dn, nDOFX_X1, &
                   WeightsX_X1(:) * G_X1_Dn(:,iGF_h_3), 1,  One, dh3dX1, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX1_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_3), &
                    1,  One, dh3dX1, 1 )

        dh3dX1 = dh3dX1 / ( WeightsX_q(:) * dX(1) )

        ! --- Scale factor derivatives with respect to X2 ---
        ! --- dh1dX2 ---

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                   WeightsX_X2(:) * G_X2_Up(:,iGF_h_1), 1, Zero, dh1dX2, 1 )
        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                   WeightsX_X2(:) * G_X2_Dn(:,iGF_h_1), 1,  One, dh1dX2, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_2), 1,  &
                    One, dh1dX2, 1 )

        dh1dX2 = dh1dX2 / ( WeightsX_q(:) * dX(2) )

        ! --- dh2dX2 ---

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                   WeightsX_X2(:) * G_X2_Up(:,iGF_h_2), 1, Zero, dh2dX2, 1 )
        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                   WeightsX_X2(:) * G_X2_Dn(:,iGF_h_2), 1,  One, dh2dX2, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_2), 1,  &
                    One, dh2dX2, 1 )

        dh2dX2 = dh2dX2 / ( WeightsX_q(:) * dX(2) )

        ! --- dh3dX2 ---

        CALL DGEMV( 'T', nDOFX_X2, nDOFX, + One, LX_X2_Up, nDOFX_X2, &
                   WeightsX_X2(:) * G_X2_Up(:,iGF_h_3), 1, Zero, dh3dX2, 1 )
        CALL DGEMV( 'T', nDOFX_X2, nDOFX, - One, LX_X2_Dn, nDOFX_X2, &
                    WeightsX_X2(:) * G_X2_Dn(:,iGF_h_3), 1, One, dh3dX2, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX2_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_3), 1,  &
                    One, dh3dX2, 1 )

        dh3dX2 = dh3dX2 / ( WeightsX_q(:) * dX(2) )

        ! --- Scale factor derivatives with respect to X3 ---
        ! --- dh1dX3 ---

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                   WeightsX_X3(:) * G_X3_Up(:,iGF_h_1), 1, Zero, dh1dX3, 1 )
        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                   WeightsX_X3(:) * G_X3_Dn(:,iGF_h_1), 1,  One, dh1dX3, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_1), 1,  &
                    One, dh1dX3, 1 )

        dh1dX3 = dh1dX3 / ( WeightsX_q(:) * dX(3) )

        ! --- dh2dX3 ---

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                   WeightsX_X3(:) * G_X3_Up(:,iGF_h_2), 1, Zero, dh2dX3, 1 )
        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                    WeightsX_X3(:) * G_X3_Dn(:,iGF_h_2), 1, One, dh2dX3, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_2), 1,  &
                    One, dh2dX3, 1 )

        dh2dX3 = dh2dX3 / ( WeightsX_q(:) * dX(3) )

        ! --- dh3dx3 ---

        CALL DGEMV( 'T', nDOFX_X3, nDOFX, + One, LX_X3_Up, nDOFX_X3, &
                   WeightsX_X3(:) * G_X3_Up(:,iGF_h_3), 1, Zero, dh3dX3, 1 )
        CALL DGEMV( 'T', nDOFX_X3, nDOFX, - One, LX_X3_Dn, nDOFX_X3, &
                   WeightsX_X3(:) * G_X3_Dn(:,iGF_h_3), 1,  One, dh3dX3, 1 )
        CALL DGEMV( 'T', nDOFX,    nDOFX, - One, dLXdX3_q, nDOFX,    &
                    WeightsX_q (:) * G(:,iX1,iX2,iX3,iGF_h_3), 1,  &
                    One, dh3dX3, 1 )

        dh3dX3 = dh3dX3 / ( WeightsX_q(:) * dX(3) )

      END IF ! --- End of UseSourceTerm IF statement ---

      DO iNodeX = 1, nDOFX

        IF( UseSourceTerm )THEN
          ! --- Contribution from source-term (FIX THIS FOR MULTI-D) ---

          PressureTensor(1) = ( U(iNodeX,iX1,iX2,iX3,iCF_S1) &
                                  * P(iNodeX,iPF_V1) + A(iNodeX,iAF_P) ) &
                                / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11)
          PressureTensor(2) = ( U(iNodeX,iX1,iX2,iX3,iCF_S2) &
                                  * P(iNodeX,iPF_V2) + A(iNodeX,iAF_P) ) &
                                / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22)
          PressureTensor(3) = ( U(iNodeX,iX1,iX2,iX3,iCF_S3) &
                                  * P(iNodeX,iPF_V3) + A(iNodeX,iAF_P) ) &
                                / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33)

          SourceTerm(1) &
            =   PressureTensor(1) * G(iNodeX,iX1,iX2,iX3,iGF_h_1) &
                  * dh1dx1(iNodeX) &
              + PressureTensor(2) * G(iNodeX,iX1,iX2,iX3,iGF_h_2) &
                  * dh2dx1(iNodeX) &
              + PressureTensor(3) * G(iNodeX,iX1,iX2,iX3,iGF_h_3) &
                  * dh3dx1(iNodeX)
          SourceTerm(2) &
            =   PressureTensor(1) * G(iNodeX,iX1,iX2,iX3,iGF_h_1) &
                  * dh1dx2(iNodeX) &
              + PressureTensor(2) * G(iNodeX,iX1,iX2,iX3,iGF_h_2) &
                  * dh2dx2(iNodeX) &
              + PressureTensor(3) * G(iNodeX,iX1,iX2,iX3,iGF_h_3) &
                  * dh3dx2(iNodeX)
          SourceTerm(3) &
            =   PressureTensor(1) * G(iNodeX,iX1,iX2,iX3,iGF_h_1) &
                  * dh1dx3(iNodeX) &
              + PressureTensor(2) * G(iNodeX,iX1,iX2,iX3,iGF_h_2) &
                  * dh2dx3(iNodeX) &
              + PressureTensor(3) * G(iNodeX,iX1,iX2,iX3,iGF_h_3) &
                  * dh3dx3(iNodeX)

          ! --- Coefficients in quadratic equation ---
          a2 = (   SourceTerm(1)**2 / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11)   &
                 + SourceTerm(2)**2 / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22)   &
                 + SourceTerm(3)**2 / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) ) &
               / ( One - epsilon )**2

          a1 = Two / ( One - epsilon ) &
                 * (  SourceTerm(1) * U(iNodeX,iX1,iX2,iX3,iCF_S1) &
                      / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
                   +  SourceTerm(2) * U(iNodeX,iX1,iX2,iX3,iCF_S2) &
                     / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
                   +  SourceTerm(3) * U(iNodeX,iX1,iX2,iX3,iCF_S3) &
                     / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )

          a0 =   U(iNodeX,iX1,iX2,iX3,iCF_S1)**2 &
                   / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) &
               + U(iNodeX,iX1,iX2,iX3,iCF_S2)**2 &
                   / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) &
               + U(iNodeX,iX1,iX2,iX3,iCF_S3)**2 &
                   / G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) &
               - U(iNodeX,iX1,iX2,iX3,iCF_E) &
                 * ( U(iNodeX,iX1,iX2,iX3,iCF_E) &
                      + Two * U(iNodeX,iX1,iX2,iX3,iCF_D) )

          IF( a1**2 - Four * a2 * a0 .LT. Zero )THEN
            WRITE(*,*) 'Negative discriminant'
            WRITE(*,'(A,ES24.16E3)') 'a = ', a2
            WRITE(*,'(A,ES24.16E3)') 'b = ', a1
            WRITE(*,'(A,ES24.16E3)') 'c = ', a0
            WRITE(*,'(A)') 'Stopping...'
            STOP
          END IF

          IF( ABS(a2) .GT. SqrtTiny )THEN
            NegRoot(1) &
              = ( -a1 - SQRT( a1**2 - Four * a2 * a0 ) ) / ( Two * a2 )
            PosRoot(1) &
              = ( -a1 + SQRT( a1**2 - Four * a2 * a0 ) ) / ( Two * a2 )
          ELSE
            NegRoot(1) = HUGE( One )
            PosRoot(1) = -a0 / a1
          END IF

          IF( NegRoot(1) .LT. Zero ) NegRoot(1) = HUGE( One )
          IF( PosRoot(1) .LT. Zero ) PosRoot(1) = HUGE( One )

          IF( SourceTerm(1) .GT. SqrtTiny )THEN
            dt_S(1) = MIN( dt_S(1), MIN( PosRoot(1), NegRoot(1) ) )
          ELSE
            dt_S(1) = HUGE( One )
          END IF

        END IF ! --- End of UseSourceTerm IF statement ---

        ! --- Contribution from numerical flux term ---

        ! --- Eq. (2.7) from Qin et al., (2016), JCP, 315, 323 ---
        W = 1.0d0 / SQRT( 1.0d0 &
            - ( G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11) * P(iNodeX,iPF_V1)**2 &
              + G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22) * P(iNodeX,iPF_V2)**2 &
              + G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) * P(iNodeX,iPF_V3)**2 ) )
        h = 1.0d0 + ( P(iNodeX,iPF_E) + A(iNodeX,iAF_P) ) / P(iNodeX,iPF_D)
        alpha = ( ABS( P(iNodeX,iPF_V1) ) &
                  * ( h + 1.0d0 - 2.0d0 * h * tau ) &
                  * W**2 + SQRT( tau**4 * ( h - 1.0d0 )**2 &
                  + tau**2 * ( h - 1.0d0 ) &
                  * ( h + 1.0d0 - 2.0d0 * h * tau ) ) ) &
                  / ( W**2 * ( h + 1.0d0 - 2.0d0 * h * tau ) &
                  + tau**2 * ( h - 1.0d0 ) )

        EigVals_X1(:,iNodeX) = Eigenvalues_GR &
                                 ( P(iNodeX,iPF_V1), &
                                   P(iNodeX,iPF_V2), &
                                   P(iNodeX,iPF_V3), &
                                   P(iNodeX,iPF_V1), &
                                   A(iNodeX,iAF_Cs), &
                                   G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                                   G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                                   G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                                   G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                                   G(iNodeX,iX1,iX2,iX3,iGF_Alpha),    &
                                   G(iNodeX,iX1,iX2,iX3,iGF_Beta_1) )

        ! --- Maximum wavespeed ---
        Max_X1 &
          = MAX( Max_X1, MAX( ABS( alpha ), &
              MAXVAL( ABS( EigVals_X1(:,iNodeX) ) ) ) )

        IF( nDimsX .GT. 1 )THEN
          EigVals_X2(:,iNodeX) = Eigenvalues_GR &
                                   ( P(iNodeX,iPF_V1), &
                                     P(iNodeX,iPF_V2), &
                                     P(iNodeX,iPF_V3), &
                                     P(iNodeX,iPF_V2), &
                                     A(iNodeX,iAF_Cs), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Alpha),    &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Beta_2) )
          Max_X2 &
            = MAX( Max_X2, MAX( ABS( alpha ), &
                MAXVAL( ABS( EigVals_X2(:,iNodeX) ) ) ) )
        END IF

        IF( nDimsX .GT. 2 )THEN
          EigVals_X3(:,iNodeX) = Eigenvalues_GR &
                                   ( P(iNodeX,iPF_V1), &
                                     P(iNodeX,iPF_V2), &
                                     P(iNodeX,iPF_V3), &
                                     P(iNodeX,iPF_V3), &
                                     A(iNodeX,iAF_Cs), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Alpha),    &
                                     G(iNodeX,iX1,iX2,iX3,iGF_Beta_3) )
          Max_X3 &
            = MAX( Max_X3, MAX( ABS( alpha ), &
                MAXVAL( ABS( EigVals_X3(:,iNodeX) ) ) ) )
        END IF

      END DO ! Loop over nNodesX

      dt_X(1) &
        = dX(1) * WeightsLX1(1) / Max_X1

      IF( nDimsX .GT. 1 ) &
        dt_X(2) &
          = dX(2) * WeightsLX2(1) / Max_X2

      IF( nDimsX .GT. 2 ) &
        dt_X(3) &
          = dX(3) * WeightsLX3(1) / Max_X3

      TimeStep = MIN( TimeStep, MIN( MINVAL( dt_X ), MINVAL( dt_S ) ) )
      IF( TimeStep .LT. SqrtTiny )THEN

        WRITE(*,*)
        WRITE(*,*) 'iX1, iX2, iX3       = ', iX1, iX2, iX3
        WRITE(*,*) 'dt_S                = ', dt_S
        IF( UseSourceTerm ) &
          WRITE(*,*) 'MinRoot1            = ', MIN( PosRoot(1), NegRoot(1) )
        WRITE(*,*) 'D                   = ', U(:,iX1,iX2,iX3,iCF_D)
        WRITE(*,*) 'S1                  = ', U(:,iX1,iX2,iX3,iCF_S1)
        WRITE(*,*) 'tau                 = ', U(:,iX1,iX2,iX3,iCF_E)
        WRITE(*,*) 'Gm11                = ', G(:,iX1,iX2,iX3,iGF_Gm_dd_11)
        WRITE(*,*) 'Gm22                = ', G(:,iX1,iX2,iX3,iGF_Gm_dd_22)
        WRITE(*,*) 'Gm33                = ', G(:,iX1,iX2,iX3,iGF_Gm_dd_33)
        WRITE(*,*) 'SQRT(tau*(tau+2*D)) = ', &
             SQRT( U(:,iX1,iX2,iX3,iCF_E) * ( U(:,iX1,iX2,iX3,iCF_E) &
             + Two * U(:,iX1,iX2,iX3,iCF_D) ) )
        STOP 'Timestep < SqrtTiny'

      END IF

    END DO
    END DO
    END DO

    TimeStep = MAX( CFL * epsilon * TimeStep, SqrtTiny )

  END SUBROUTINE ComputeTimeStep_GR

  
  SUBROUTINE ComputeConserved_GR( PF_D, PF_V1, PF_V2, PF_V3, PF_E, PF_Ne, &
                                  CF_D, CF_S1, CF_S2, CF_S3, CF_E, CF_Ne, &
                                  GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33,  &
                                  AF_P )

    REAL(DP), DIMENSION(:), INTENT(in)  :: PF_D, PF_V1, PF_V2, PF_V3, &
                                           PF_E, PF_Ne
    REAL(DP), DIMENSION(:), INTENT(out) :: CF_D, CF_S1, CF_S2, CF_S3, &
                                           CF_E, CF_Ne
    REAL(DP), DIMENSION(:), INTENT(in)  :: GF_Gm_dd_11, GF_Gm_dd_22, GF_Gm_dd_33
    REAL(DP), DIMENSION(:), INTENT(in)  :: AF_P

    REAL(DP), DIMENSION( 1 : SIZE(PF_D) ) :: vSq, W, h

    vSq = GF_Gm_dd_11 * PF_V1**2 &
          + GF_Gm_dd_22 * PF_V2**2 &
            + GF_Gm_dd_33 * PF_V3**2

    W = 1.0_DP / SQRT( 1.0_DP - vSq )
    h = 1.0_DP + ( PF_E + AF_P ) / PF_D
    
    CF_D   = W * PF_D
    CF_S1  = h * W**2 * PF_D * GF_Gm_dd_11 * PF_V1
    CF_S2  = h * W**2 * PF_D * GF_Gm_dd_22 * PF_V2
    CF_S3  = h * W**2 * PF_D * GF_Gm_dd_33 * PF_V3
    CF_E   = h * W**2 * PF_D - AF_P - W * PF_D
    CF_Ne  = W * PF_Ne

  END SUBROUTINE ComputeConserved_GR
  

  PURE FUNCTION Eigenvalues_GR &
    ( V_1, V_2, V_3, V_i, Cs, Gm_11, Gm_22, Gm_33, Gm_ii, Alpha, Beta_i )

    ! Alpha is the lapse function
    ! V_i is the contravariant component V^i
    ! Beta_1 is the contravariant component Beta^1

    REAL(DP)             :: Eigenvalues_GR(1:nCF)
    REAL(DP), INTENT(in) :: V_1, V_2, V_3, V_i, Cs
    REAL(DP), INTENT(in) :: Gm_11, Gm_22, Gm_33, Gm_ii, Alpha, Beta_i

    REAL(DP) :: VSq

    VSq = Gm_11 * V_1**2 + Gm_22 * V_2**2 + Gm_33 * V_3**2

    Eigenvalues_GR(1) &
      = Alpha / ( One - VSq * Cs**2 ) * ( V_i * ( One - Cs**2 ) &
        - Cs * SQRT( ( One - VSq ) * ( ( One - VSq * Cs**2 ) / Gm_ii &
           - V_i**2 * ( One - Cs**2 ) ) ) ) - Beta_i

    Eigenvalues_GR(2) &
      = Alpha * V_i - Beta_i

    Eigenvalues_GR(3) &
      = Alpha / ( One - VSq * Cs**2 ) * ( V_i * ( One - Cs**2 ) &
        + Cs * SQRT( ( One - VSq ) * ( ( One - VSq * Cs**2 ) / Gm_ii &
           - V_i**2 * ( One - Cs**2 ) ) ) ) - Beta_i

    Eigenvalues_GR(4) &
      = Alpha * V_i - Beta_i

    Eigenvalues_GR(5) &
      = Alpha * V_i - Beta_i

    Eigenvalues_GR(6) &
      = Alpha * V_i - Beta_i

    RETURN
  END FUNCTION Eigenvalues_GR


  FUNCTION Flux_X1_GR &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in)       :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in)       :: Gm11, Gm22, Gm33, Lapse, Shift
    REAL(DP), DIMENSION(1:nCF) :: Flux_X1_GR

    REAL(DP) :: vSq, W, h

    vSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    W   = 1.0d0 / SQRT( 1.0d0 - vSq )

    h   = 1.0d0 + ( E + P ) / D

    Flux_X1_GR(iCF_D)  &
      = D * W * ( V1 - Shift / Lapse )

    Flux_X1_GR(iCF_S1) &
      = D * h * W**2 * Gm11 * V1  * ( V1 - Shift / Lapse ) + P

    Flux_X1_GR(iCF_S2) &
      = D * h * W**2 * Gm22 * V2  * ( V1 - Shift / Lapse )

    Flux_X1_GR(iCF_S3) &
      = D * h * W**2 * Gm33 * V3  * ( V1 - Shift / Lapse )

    Flux_X1_GR(iCF_E)  &
      = D * W * ( h * W - 1.0d0 ) * ( V1 - Shift / Lapse ) + Shift / Lapse * P

    Flux_X1_GR(iCF_Ne) &
      = Ne * W * ( V1 - Shift / Lapse )

    RETURN
  END FUNCTION Flux_X1_GR


  FUNCTION Flux_X2_GR &
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Lapse, Shift )

    REAL(DP), INTENT(in)       :: D, V1, V2, V3, E, Ne, P
    REAL(DP), INTENT(in)       :: Gm11, Gm22, Gm33, Lapse, Shift
    REAL(DP), DIMENSION(1:nCF) :: Flux_X2_GR

    REAL(DP) :: vSq, W, h

    vSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    W   = 1.0d0 / SQRT( 1.0d0 - vSq )

    h   = 1.0d0 + ( E + P ) / D

    Flux_X2_GR(iCF_D)  &
      = D * W * ( V2 - Shift / Lapse )

    Flux_X2_GR(iCF_S1) &
      = D * h * W**2 * Gm11 * V1  * ( V2 - Shift / Lapse )

    Flux_X2_GR(iCF_S2) &
      = D * h * W**2 * Gm22 * V2  * ( V2 - Shift / Lapse ) + P

    Flux_X2_GR(iCF_S3) &
      = D * h * W**2 * Gm33 * V3  * ( V2 - Shift / Lapse )

    Flux_X2_GR(iCF_E)  &
      = D * W * ( h * W - 1.0d0 ) * ( V2 - Shift / Lapse ) + Shift / Lapse * P

    Flux_X2_GR(iCF_Ne) &
      = Ne * W * ( V2 - Shift / Lapse )

    RETURN
  END FUNCTION Flux_X2_GR


  PURE FUNCTION StressTensor_Diagonal( S_1, S_2, S_3, V_1, V_2, V_3, P )

    REAL(DP)             :: StressTensor_Diagonal(1:3)
    REAL(DP), INTENT(in) :: S_1, S_2, S_3, V_1, V_2, V_3, P

    StressTensor_Diagonal(1) = S_1 * V_1 + P
    StressTensor_Diagonal(2) = S_2 * V_2 + P
    StressTensor_Diagonal(3) = S_3 * V_3 + P

    RETURN
  END FUNCTION StressTensor_Diagonal


  REAL(DP) FUNCTION AlphaC_GR &
             ( U_L, U_R, F_L, F_R, Gmdd, Lapse, Shift, aP, aM )

    ! --- Middle Wavespeed as suggested by Mignone and Bodo (2005) ---

    REAL(DP), INTENT(in) :: U_L(nCF), U_R(nCF), &
                            F_L(nCF), F_R(nCF), &
                            Gmdd, Lapse, Shift, &
                            aP, aM

    REAL(DP) :: E_L, E_R, F_E_L, F_E_R, A_L, A_R, B_L, B_R
    REAL(DP) :: A, B, C, eps = SqrtTiny

    E_L   = U_L(iCF_E) + U_L(iCF_D)
    E_R   = U_R(iCF_E) + U_R(iCF_D)
    F_E_L = F_L(iCF_E) + F_L(iCF_D)
    F_E_R = F_R(iCF_E) + F_R(iCF_D)

    A_L = -aM * E_L - Lapse * F_E_L
    A_R = +aP * E_R - Lapse * F_E_R
    B_L = -aM * U_L(iCF_S1) - Lapse * F_L(iCF_S1)
    B_R = +aP * U_R(iCF_S1) - Lapse * F_R(iCF_S1)

    ! --- A, B, and C from quadratic equation ---
    A = Gmdd**2 &
          * ( A_R * ( -aM + Shift ) - A_L * ( aP + Shift ) )
    B = Gmdd &
          * ( ( Lapse * A_L - B_R * ( -aM + Shift ) ) &
            - ( Lapse * A_R - B_L * (  aP + Shift ) ) )
    C = Lapse * ( B_R - B_L )

    ! --- Accounting for special cases of the solution to a
    !     quadratic equation when A = 0 ---

    IF     ( ( ABS( A ) .LT. eps ) .AND. ( ABS( B ) .LT. eps ) &
            .AND. ( ABS( C ) .LT. eps ) )THEN
      WRITE(*,*) 'AlphaC is undefined'
      AlphaC_GR = 0.0_DP
    ELSE IF( ( ABS( A ) .LT. eps ) .AND. ( ABS( B ) .LT. eps ) )THEN
      WRITE(*,*) 'AlphaC is undefined'
      WRITE(*,*) 'C:', C
      AlphaC_GR = 0.0_DP
    ELSE IF( ( ABS( A ) .LT. eps ) .AND. ( ABS( C ) .LT. eps ) )THEN
      AlphaC_GR = 0.0_DP
    ELSE IF( ABS( A ) .LT. eps )THEN
      AlphaC_GR = -C / B
    ELSE
      AlphaC_GR = ( -B - SQRT( MAX( B**2 - 4.0_DP * A * C, eps ) ) ) &
                    / ( 2.0_DP * A )
    END IF

    RETURN
  END FUNCTION AlphaC_GR


  PURE FUNCTION NumericalFlux_LLF_GR &
      ( UL, UR, FL, FR, aP, aM, aC, nF, &
        vL, vR, pL, pR, Lapse, Shift, Gmdd )

    ! --- Local Lax-Friedrichs Flux ---

    INTEGER,  INTENT(in) :: nF
    REAL(DP), INTENT(in) :: UL(1:nF), UR(1:nF), FL(1:nF), FR(1:nF), &
                            aP, aM, aC, vL, vR, pL, pR, &
                            Lapse, Shift, Gmdd
    REAL(DP)             :: p, D, S1, S2, S3, E, Ne
    REAL(DP)             :: NumericalFlux_LLF_GR(1:nF)
    REAL(DP)             :: alpha

    alpha = MAX( aM, aP )

    NumericalFlux_LLF_GR &
      = 0.5_DP * ( FL + FR - alpha * ( UR - UL ) )

    RETURN
  END FUNCTION NumericalFlux_LLF_GR


  PURE FUNCTION NumericalFlux_HLL_GR &
      ( UL, UR, FL, FR, aP, aM, aC, nF, &
        vL, vR, pL, pR, Lapse, Shift, Gmdd )

    INTEGER,  INTENT(in) :: nF
    REAL(DP), INTENT(in) :: UL(1:nF), UR(1:nF), FL(1:nF), FR(1:nF), &
                            aP, aM, aC, vL, vR, pL, pR, &
                            Lapse, Shift, Gmdd
    REAL(DP)             :: p, D, S1, S2, S3, E, Ne
    REAL(DP)             :: NumericalFlux_HLL_GR(1:nF)

    NumericalFlux_HLL_GR &
      = ( aP * FL + aM * FR &
            - aP * aM * ( UR - UL ) ) / ( aP + aM )

    RETURN
  END FUNCTION NumericalFlux_HLL_GR


  PURE FUNCTION NumericalFlux_X1_HLLC_GR &
      ( UL, UR, FL, FR, aP, aM, aC, nF, &
        vL, vR, pL, pR, Lapse, Shift, Gmdd )

    INTEGER,  INTENT(in) :: nF
    REAL(DP), INTENT(in) :: UL(1:nF), UR(1:nF), FL(1:nF), FR(1:nF), &
                            aP, aM, aC, vL, vR, pL, pR, &
                            Lapse, Shift, Gmdd
    REAL(DP)             :: p, D, S1, S2, S3, E, Ne
    REAL(DP)             :: UE, FE, FS
    REAL(DP)             :: VelocityRatio
    REAL(DP)             :: NumericalFlux_X1_HLLC_GR(1:nF)

    IF( aM .EQ. 0.0d0 )THEN

      NumericalFlux_X1_HLLC_GR = FL

    ELSEIF( aP .EQ. 0.0d0 )THEN

      NumericalFlux_X1_HLLC_GR = FR

    ELSE

      ! --- From Mignone & Bodo (2005)
      !     Note the sign change on aM which is due to it being
      !     read in as positive but the formulae assuming it is negative ---

      IF( aC .GE. 0.0d0 )THEN

        VelocityRatio = ( -aM - Lapse * vL + Shift ) &
                      / ( -aM - Lapse * aC + Shift )

        ! --- UL_star ---
        UE = UL(iCF_E) + UL(iCF_D)
        FE = UL(iCF_S1) / Gmdd - Shift / Lapse * UE
        FS = UL(iCF_S1) * ( vL - Shift / Lapse ) + pL

        p = ( Gmdd * aC * ( -aM * UE - Lapse * FE ) &
              - ( -aM * UL(iCF_S1) - Lapse * FS ) ) &
            / ( Lapse - Gmdd * aC * ( -aM + Shift ) )

        D = UL(iCF_D)   * VelocityRatio

        S1 = UL(iCF_S1) * VelocityRatio &
               + Lapse * ( p - pL ) &
               / ( -aM - Lapse * aC + Shift )

        S2 = UL(iCF_S2) * VelocityRatio

        S3 = UL(iCF_S3) * VelocityRatio

        E = UE          * VelocityRatio &
              + Lapse * ( p * aC - pL * vL ) &
              / ( -aM - Lapse * aC + Shift )

        Ne = UL(iCF_Ne) * VelocityRatio

      ELSE

        VelocityRatio = ( aP - Lapse * vR + Shift ) &
                      / ( aP - Lapse * aC + Shift )

        ! --- UR_star ---
        UE = UR(iCF_E) + UR(iCF_D)
        FE = UR(iCF_S1) / Gmdd - Shift / Lapse * UE
        FS = UR(iCF_S1) * ( vR - Shift / Lapse ) + pR

        p = ( Gmdd * aC * ( aP * UE - Lapse * FE ) &
              - ( aP * UR(iCF_S1) - Lapse * FS ) ) &
              / ( Lapse - Gmdd * aC * ( aP + Shift ) )

        D = UR(iCF_D)   * VelocityRatio

        S1 = UR(iCF_S1) * VelocityRatio &
               + Lapse * ( p - pR ) &
               / ( aP - Lapse * aC + Shift )

        S2 = UR(iCF_S2) * VelocityRatio

        S3 = UR(iCF_S3) * VelocityRatio

        E = UE          * VelocityRatio &
              + Lapse * ( p * aC - pR * vR ) &
              / ( aP - Lapse * aC + Shift )

        Ne = UR(iCF_Ne) * VelocityRatio

      END IF

      NumericalFlux_X1_HLLC_GR(iCF_D)  &
        = D  * ( aC - Shift / Lapse )
      NumericalFlux_X1_HLLC_GR(iCF_S1) &
        = S1 * ( aC - Shift / Lapse ) + p
      NumericalFlux_X1_HLLC_GR(iCF_S2) &
        = S2 * ( aC - Shift / Lapse )
      NumericalFlux_X1_HLLC_GR(iCF_S3) &
        = S3 * ( aC - Shift / Lapse )
      NumericalFlux_X1_HLLC_GR(iCF_E)  &
        = S1 / Gmdd - D * aC - Shift / Lapse * ( E - D )
      NumericalFlux_X1_HLLC_GR(iCF_Ne) &
        = Ne * ( aC - Shift / Lapse )

    END IF


  END FUNCTION NumericalFlux_X1_HLLC_GR


  PURE FUNCTION NumericalFlux_X2_HLLC_GR &
      ( UL, UR, FL, FR, aP, aM, aC, nF, &
        vL, vR, pL, pR, Lapse, Shift, Gmdd )

    INTEGER,  INTENT(in) :: nF
    REAL(DP), INTENT(in) :: UL(1:nF), UR(1:nF), FL(1:nF), FR(1:nF), &
                            aP, aM, aC, vL, vR, pL, pR, &
                            Lapse, Shift, Gmdd
    REAL(DP)             :: p, D, S1, S2, S3, E, Ne
    REAL(DP)             :: UE, FE, FS
    REAL(DP)             :: VelocityRatio
    REAL(DP)             :: NumericalFlux_X2_HLLC_GR(1:nF)

    IF( aM .EQ. 0.0d0 )THEN

      NumericalFlux_X2_HLLC_GR = FL

    ELSEIF( aP .EQ. 0.0d0 )THEN

      NumericalFlux_X2_HLLC_GR = FR

    ELSE

      ! --- From Mignone & Bodo (2005)
      !     Note the sign change on aM which is due to it being
      !     read in as positive but the formulae assuming it is negative ---

      IF( aC .GE. 0.0d0 )THEN

        VelocityRatio = ( -aM - Lapse * vL + Shift ) &
                      / ( -aM - Lapse * aC + Shift )

        ! --- UL_star ---
        UE = UL(iCF_E) + UL(iCF_D)
        FE = UL(iCF_S2) / Gmdd - Shift / Lapse * UE
        FS = UL(iCF_S2) * ( vL - Shift / Lapse ) + pL

        p = ( Gmdd * aC * ( -aM * UE - Lapse * FE ) &
              - ( -aM * UL(iCF_S2) - Lapse * FS ) ) &
            / ( Lapse - Gmdd * aC * ( -aM + Shift ) )

        D = UL(iCF_D)   * VelocityRatio

        S1 = UL(iCF_S1) * VelocityRatio

        S2 = UL(iCF_S2) * VelocityRatio &
               + Lapse * ( p - pL ) &
               / ( -aM - Lapse * aC + Shift )

        S3 = UL(iCF_S3) * VelocityRatio

        E = UE          * VelocityRatio &
              + Lapse * ( p * aC - pL * vL ) &
              / ( -aM - Lapse * aC + Shift )

        Ne = UL(iCF_Ne) * VelocityRatio

      ELSE

        VelocityRatio = ( aP - Lapse * vR + Shift ) &
                      / ( aP - Lapse * aC + Shift )

        ! --- UR_star ---
        UE = UR(iCF_E) + UR(iCF_D)
        FE = UR(iCF_S2) / Gmdd - Shift / Lapse * UE
        FS = UR(iCF_S2) * ( vR - Shift / Lapse ) + pR

        p = ( Gmdd * aC * ( aP * UE - Lapse * FE ) &
              - ( aP * UR(iCF_S2) - Lapse * FS ) ) &
              / ( Lapse - Gmdd * aC * ( aP + Shift ) )

        D = UR(iCF_D)   * VelocityRatio

        S1 = UR(iCF_S1) * VelocityRatio

        S2 = UR(iCF_S2) * VelocityRatio &
               + Lapse * ( p - pR ) &
               / ( aP - Lapse * aC + Shift )

        S3 = UR(iCF_S3) * VelocityRatio

        E = UE          * VelocityRatio &
              + Lapse * ( p * aC - pR * vR ) &
              / ( aP - Lapse * aC + Shift )

        Ne = UR(iCF_Ne) * VelocityRatio

      END IF

      NumericalFlux_X2_HLLC_GR(iCF_D)  &
        = D  * ( aC - Shift / Lapse )
      NumericalFlux_X2_HLLC_GR(iCF_S1) &
        = S1 * ( aC - Shift / Lapse )
      NumericalFlux_X2_HLLC_GR(iCF_S2) &
        = S2 * ( aC - Shift / Lapse ) + p
      NumericalFlux_X2_HLLC_GR(iCF_S3) &
        = S3 * ( aC - Shift / Lapse )
      NumericalFlux_X2_HLLC_GR(iCF_E)  &
        = S2 / Gmdd - D * aC - Shift / Lapse * ( E - D )
      NumericalFlux_X2_HLLC_GR(iCF_Ne) &
        = Ne * ( aC - Shift / Lapse )

    END IF


  END FUNCTION NumericalFlux_X2_HLLC_GR


END MODULE EulerEquationsUtilitiesModule_Beta_GR
