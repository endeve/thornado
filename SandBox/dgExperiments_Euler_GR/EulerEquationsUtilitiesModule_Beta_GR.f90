MODULE EulerEquationsUtilitiesModule_Beta_GR

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One, Two
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDimsX
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    nAF, iAF_P, iAF_Cs
  USE GeometryFieldsModule, ONLY: &
       iGF_h_1, iGF_h_2, iGF_h_3, &
       iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33, &
       iGF_Alpha, &
       iGF_Beta_1, iGF_Beta_2, iGF_Beta_3, nGF
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromSpecificInternalEnergy, &
    ComputeSoundSpeedFromPrimitive_GR
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL

  IMPLICIT NONE
  PRIVATE :: ComputePressureWithBisectionMethod
  PRIVATE :: ComputeFunP

  PUBLIC :: ComputeFromConserved_GR
  PUBLIC :: ComputePrimitive_GR
  PUBLIC :: ComputeFunJacP
  PUBLIC :: ComputeConserved_GR
  PUBLIC :: ComputeTimeStep_GR
  PUBLIC :: Eigenvalues_GR
  PUBLIC :: Flux_X1_GR
  PUBLIC :: StressTensor_Diagonal
  PUBLIC :: AlphaC_GR
  PUBLIC :: NumericalFlux_X1_LLF_GR
  PUBLIC :: NumericalFlux_X1_HLL_GR
  PUBLIC :: NumericalFlux_X1_HLLC_GR

  LOGICAL :: DEBUG = .FALSE.
  INTEGER, PARAMETER :: MAX_IT = 1000

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

    REAL(DP) :: SSq, Pold, vSq, W, h, Pnew, q, FunP0

    REAL(DP) :: FunP, JacP
    REAL(DP), PARAMETER :: TolP = 1.0d-10, TolFunP = 1.0d-6, TolFunP0 = 1.0d-10

    ! --- Loop through all the nodes ---
    nNodes = SIZE( CF_D )
    DO i = 1, nNodes

      q = CF_E(i) + CF_D(i) - SQRT( CF_D(i)**2 &
                                  + CF_S1(i)**2 / GF_Gm_dd_11(i)  &
                                  + CF_S2(i)**2 / GF_Gm_dd_22(i)  &
                                  + CF_S3(i)**2 / GF_Gm_dd_33(i) )

      IF( q .LT. Zero )THEN
        WRITE(*,*)
        WRITE(*,'(A6,I1)')        'Node: ', i
        WRITE(*,'(A9,ES18.10E3)') 'q:       ', q
        WRITE(*,'(A9,ES18.10E3)') 'Gm11(i): ', GF_Gm_dd_11(i)
        WRITE(*,'(A9,ES18.10E3)') 'Gm22(i): ', GF_Gm_dd_22(i)
        WRITE(*,'(A9,ES18.10E3)') 'Gm33(i): ', GF_Gm_dd_33(i)
        WRITE(*,'(A9,ES18.10E3)') 'D(i):    ', CF_D(i)
        WRITE(*,'(A9,ES18.10E3)') 'E(i):    ', CF_E(i)
        WRITE(*,'(A9,ES18.10E3)') 'S1(i):   ', CF_S1(i)
        WRITE(*,'(A9,ES18.10E3)') 'S2(i):   ', CF_S2(i)
        WRITE(*,'(A9,ES18.10E3)') 'S3(i):   ', CF_S3(i)
        STOP 'q < 0'
      END IF
    
      SSq = CF_S1(i)**2 / GF_Gm_dd_11(i) &
            + CF_S2(i)**2 / GF_Gm_dd_22(i) &
            + CF_S3(i)**2 / GF_Gm_dd_33(i)

      ! --- Find Pressure with Newton's Method ---

      ! --- Approximation for pressure assuming h^2~=1 ---
      Pold = MAX( -q, SqrtTiny )
      !WRITE(*,*)
      !WRITE(*,'(A13,ES24.16E3)') 'Pold_approx: ', Pold
      CALL ComputePressureWithBisectionMethod( CF_D(i), CF_E(i), SSq, Pold )
      !WRITE(*,'(A13,ES24.16E3)') 'Pold_bisec: ', Pold

      IF( DEBUG ) WRITE(*,*) 'Initial guess for pressure: ', Pold
      ! --- Get initial value for FunP ---
      IF( DEBUG ) WRITE(*,'(A)') 'CALL ComputeFunJacP0'
      CALL ComputeFunJacP( CF_D(i), CF_E(i), SSq, Pold, FunP0, JacP )
      IF( DEBUG )THEN
        WRITE(*,'(A)') 'Finished ComputeFunJacP0'
        WRITE(*,'(A7,ES24.16E3)') 'FunP0: ', FunP0
        WRITE(*,*)
        WRITE(*,'(A)') "Start finding root with Newton's method"
        WRITE(*,'(A)') "---------------------------------------"
      END IF

      CONVERGED = .FALSE.

      IF( ABS( FunP0 ) .LT. TolFunP0 )THEN
        Pnew = Pold
        CONVERGED = .TRUE.
      END IF

      ITERATION = 0

      DO WHILE ( .NOT. CONVERGED )

        ITERATION = ITERATION + 1

        IF( DEBUG )THEN
          WRITE(*,*)
          WRITE(*,*) 'Iteration ', ITERATION
          WRITE(*,'(A)') 'ComputeFunJacP'
        END IF
        CALL ComputeFunJacP( CF_D(i), CF_E(i), SSq, Pold, FunP, JacP )
        IF( DEBUG ) WRITE(*,'(A)') 'Finished ComputeFunJacP'

        Pnew = Pold - FunP / JacP

        IF( DEBUG )THEN
          WRITE(*,'(A6,ES24.16E3)') 'Pnew: ', Pnew
          WRITE(*,*)
        END IF
        Pnew = MAX( Pnew, SqrtTiny )

        ! --- Check if Newton's method has converged ---
        IF( DEBUG ) WRITE(*,'(A)') 'Checking relative tolerance...'

        IF( ABS( Pnew - Pold ) / ABS( Pnew ) .LT. TolP )THEN

          IF( DEBUG ) WRITE(*,'(A)') 'Relative tolerance met!'
          IF( DEBUG ) WRITE(*,'(A15,ES23.16E3)') &
               '  |dP|/|Pnew|: ', ABS( Pnew - Pold ) / ABS( Pnew )

          IF( DEBUG ) WRITE(*,'(A)') 'Checking FunP tolerance...'
          CALL ComputeFunJacP( CF_D(i), CF_E(i), SSq, Pnew, FunP, JacP )

          IF( ABS( FunP / FunP0 ) .LT. TolFunP ) THEN
            CONVERGED = .TRUE.
            IF( DEBUG )THEN
              WRITE(*,'(A)') 'FunP tolerance met!' 
              WRITE(*,'(A16,ES24.16E3)') '  |FunP/FunP0|: ', &
                ABS( FunP / FunP0 )
              WRITE(*,'(A13,I3.3,A11)') &
                'Converged in ', ITERATION, ' iterations'
              WRITE(*,'(A)') &
                '---------------------------'
              WRITE(*,*)
            END IF

          ELSE
            IF( DEBUG )THEN
              WRITE(*,'(A)') 'FunP tolerance NOT met'
              WRITE(*,'(I3.3,A11)') ITERATION, ' iterations' 
              WRITE(*,'(A7,ES24.16E3)') 'Pold:         ', &
                Pold
              WRITE(*,'(A7,ES24.16E3)') 'Pnew:         ', &
                Pnew
              WRITE(*,'(A7,ES24.16E3)') '|dP/P|:       ', &
                ABS( Pnew - Pold ) / ABS( Pnew )
            END IF
            WRITE(*,*)
            WRITE(*,'(A14,ES24.16E3)') '|FunP/FunP0|: ', &
              ABS( FunP / FunP0 )
            STOP 'Stopping because ABS( FunP / FunP0 ) .GE. TolFunP'

          END IF

        ELSE
          IF( DEBUG )THEN
            WRITE(*,'(A)') 'Relative tolerance NOT met'
             WRITE(*,'(A13,ES23.16E3)') &
               '|dP|/|Pnew|: ', ABS( Pnew - Pold ) / ABS( Pnew )
          END IF

        END IF

        ! --- STOP after MAX_IT iterations ---
        IF( ITERATION .GE. MAX_IT - 4 )THEN

          IF( ITERATION .EQ. MAX_IT - 4 )THEN
            WRITE(*,*)
            WRITE(*,*) 'Max iterations IF statement'
            WRITE(*,*) '---------------------------'
            WRITE(*,'(A,I1)') '  Node: ', i
            WRITE(*,'(A,ES7.1E2)')  '  TolP    = ', TolP
            WRITE(*,'(A,ES7.1E2)')  '  TolFunP = ', TolFunP
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
          END IF

          WRITE(*,'(A,I4.4)') '  ITERATION: ', ITERATION
          WRITE(*,'(A)')      '  ---------------'
          WRITE(*,'(A,ES24.16E3)') '  FunP = ', FunP
          WRITE(*,'(A,ES24.16E3)') &
                  '  dP   = ', ( Pnew - Pold ) / ABS( Pnew )
          WRITE(*,*)
          IF( ITERATION == MAX_IT )THEN
            STOP 'Max allowed iterations reached, no convergence.'
          END IF
        END IF

        Pold = Pnew

      END DO

      AF_P(i) = Pnew

      !WRITE(*,'(A13,ES24.16E3)') 'P_final: ', AF_P(i)
      vSq = SSq / ( CF_E(i) + Pnew + CF_D(i) )**2

      W = 1.0_DP / SQRT( 1.0_DP - vSq )

      h = ( CF_E(i) + Pnew + CF_D(i) ) / ( W * CF_D(i) )

      IF( DEBUG )THEN
        IF( h .LT. 1.0_DP ) WRITE(*,'(A6,ES18.10E3)') 'h:    ', h
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


  SUBROUTINE ComputeFunJacP( CF_D, CF_E, SSq, P, FunP, JacP )

    REAL(DP), INTENT(in)  :: CF_D, CF_E, SSq, P
    REAL(DP), INTENT(out) :: FunP, JacP

    REAL(DP) :: HSq, RHO, EPS, dRHO, dEPS
    REAL(DP), DIMENSION(1) :: Pbar
    LOGICAL :: DEBUG_FunJacP = .FALSE.

    IF( DEBUG_FunJacP ) WRITE(*,'(A11,ES24.16E3)') '  CF_D:    ', CF_D
    IF( DEBUG_FunJacP ) WRITE(*,'(A11,ES24.16E3)') '  CF_E:    ', CF_E
    IF( DEBUG_FunJacP ) WRITE(*,'(A11,ES24.16E3)') '  SSq:     ', SSq
    IF( DEBUG_FunJacP ) WRITE(*,'(A11,ES24.16E3)') '  P:       ', P

    HSq = ( CF_E + P + CF_D )**2
    IF( DEBUG_FunJacP ) WRITE(*,'(A11,ES24.16E3)') '  HSq:     ', HSq
    IF( DEBUG_FunJacP ) WRITE(*,'(A11,ES24.16E3)') '  H^2-S^2: ', HSq - SSq

    RHO = CF_D * SQRT( HSq - SSq ) / SQRT( HSq )
    IF( DEBUG_FunJacP ) WRITE(*,'(A11,ES24.16E3)') '  RHO:     ', RHO

    EPS = ( SQRT( HSq - SSq ) &
            - P * SQRT( HSq ) / SQRT( HSq - SSq ) - CF_D ) / CF_D
    IF( DEBUG_FunJacP ) WRITE(*,'(A11,ES24.16E3)') '  EPS:     ', EPS

    CALL ComputePressureFromSpecificInternalEnergy &
         ( [ RHO ], [ EPS ], [ 0.0_DP ], Pbar )

    FunP = P - Pbar(1)
    dRHO = CF_D * SSq / ( SQRT( HSq - SSq ) * HSq )
    dEPS = P * SSq / ( ( HSq - SSq ) * SQRT( HSq ) * RHO )

    JacP = 1.0_DP - Pbar(1) * ( dRHO / RHO + dEPS / EPS )

  END SUBROUTINE ComputeFunJacP


  SUBROUTINE ComputePressureWithBisectionMethod( CF_D, CF_E, SSq, P )

    REAL(DP), INTENT(in)  :: CF_D, CF_E, SSq
    REAL(DP), INTENT(out) :: P

    INTEGER,  PARAMETER :: MAX_IT = 96
    REAL(DP), PARAMETER :: DeltaP_min = 1.0d-3

    LOGICAL :: CONVERGED
    INTEGER :: ITERATION
    REAL(DP) :: PA, PB, PC, DeltaP
    REAL(DP) :: FunPA, FunPB, FunPC

    ! --- Get upper and lower bounds on pressure, PA, PB ---
    PA = Zero
    PB = Two * ( One - One / Gamma_IDEAL ) * CF_E

    ! --- Compute FunP for upper and lower bounds ---
    CALL ComputeFunP( CF_D, CF_E, SSq, PA, FunPA )
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
    DO WHILE ( .NOT. CONVERGED )

      ITERATION = ITERATION + 1

      ! --- Compute midpoint, PC ---
      DeltaP = Half * DeltaP
      PC     = PA + DeltaP

      ! --- Compute FunP for midpoint pressure ---
      CALL ComputeFunP( CF_D, CF_E, SSq, PC, FunPC )

      ! --- Change PC = PA or PB, depending on sign of FunP(PC) ---
      IF( FunPA * FunPC .LT. Zero )THEN

        PB = PC
        FunPB = FunPC

      ELSE

        PA = PC
        FunPA = FunPC

      END IF

      IF( ( DeltaP .LT. DeltaP_min ) .OR. ( ITERATION .EQ. MAX_IT ) ) &
        CONVERGED = .TRUE.

    END DO

    P = PA
    
  END SUBROUTINE ComputePressureWithBisectionMethod

  
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

  SUBROUTINE ComputeTimeStep_GR( iX_B, iX_E, G, U, CFL, TimeStep )

    INTEGER,  INTENT(in)  :: &
      iX_B(3), iX_E(3)
    REAL(DP), INTENT(in)  :: &
      G(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:), &
      U(1:,iX_B(1):,iX_B(2):,iX_B(3):,1:)
    REAL(DP), INTENT(in)  :: &
      CFL
    REAL(DP), INTENT(out) :: &
      TimeStep

    INTEGER  :: iX1, iX2, iX3, iNodeX, iDimX
    REAL(DP) :: dX(3), dt_X(3)
    REAL(DP) :: P(nDOFX,nPF)
    REAL(DP) :: A(nDOFX,nAF)
    REAL(DP) :: EigVals_X1(nCF,nDOFX), Max_X1, &
                EigVals_X2(nCF,nDOFX), Max_X2, &
                EigVals_X3(nCF,nDOFX), Max_X3
    REAL(DP) :: epsilon = Half

    TimeStep = HUGE( One )
    dt_X(:)  = HUGE( One )

    Max_X1 = -Huge( One )
    Max_X2 = -Huge( One )
    Max_X3 = -Huge( One )

    DO iX3 = iX_B(3), iX_E(3)
      DO iX2 = iX_B(2), iX_E(2)
        DO iX1 = iX_B(1), iX_E(1)

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

          DO iNodeX = 1, nDOFX
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
            Max_X1 &
              = MAX( Max_X1, &
                       MAXVAL( ABS( EigVals_X1(:,iNodeX) ) &
                         * G(iNodeX,iX1,iX2,iX3,iGF_h_1) ) )

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
                = MAX( Max_X2, &
                       MAXVAL( ABS( EigVals_X2(:,iNodeX) ) &
                         * G(iNodeX,iX1,iX2,iX3,iGF_h_2) ) )

            END IF

            IF( nDimsX .GT. 2 )THEN
              EigVals_X3(:,iNodeX) = Eigenvalues_GR &
                                       ( P(iNodeX,iPF_V1),       &
                                         P(iNodeX,iPF_V2),       &
                                         P(iNodeX,iPF_V3),       &
                                         P(iNodeX,iPF_V3),       &
                                         A(iNodeX,iAF_Cs),       &
                                         G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                                         G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                                         G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                                         G(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                                         G(iNodeX,iX1,iX2,iX3,iGF_Alpha),    &
                                         G(iNodeX,iX1,iX2,iX3,iGF_Beta_3) )
              Max_X3 &
                   = MAX( Max_X3, &
                       MAXVAL( ABS( EigVals_X3(:,iNodeX) ) &
                         * G(iNodeX,iX1,iX2,iX3,iGF_h_3) ) )

            END IF

          END DO

          ! --- NEED TO FIX THIS FOR SOURCE TERM ---
          dt_X(1) &
            = dX(1) / ( Two * Max_X1 )

          IF( nDimsX .GT. 1 ) &
            dt_X(2) &
              = dX(2) / ( Two * Max_X2 )

          IF( nDimsX .GT. 2 ) &
            dt_X(3) &
              = dX(3) / ( Two * Max_X3 )

          TimeStep = MIN( TimeStep, MINVAL( dt_X ) )

        END DO
      END DO
    END DO

    TimeStep = epsilon * CFL * TimeStep

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
    ( D, V1, V2, V3, E, Ne, P, Gm11, Gm22, Gm33, Alpha, Beta1 )

    REAL(DP), INTENT(in)       :: D, V1, V2, V3, E, P, Ne
    REAL(DP), INTENT(in)       :: Alpha, Beta1, Gm11, Gm22, Gm33
    REAL(DP), DIMENSION(1:nCF) :: Flux_X1_GR

    REAL(DP) :: vSq, W, h

    vSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

    W   = 1.0_DP / SQRT( 1.0_DP - vSq )

    h   = 1.0_DP + ( E + P ) / D

    Flux_X1_GR(iCF_D)  &
      = D * W * ( V1 - Beta1 / Alpha )

    Flux_X1_GR(iCF_S1) &
      = D * h * W**2 * Gm11 * V1 * ( V1 - Beta1 / Alpha ) + P

    Flux_X1_GR(iCF_S2) &
      = D * h * W**2 * Gm22 * V2 * ( V1 - Beta1 / Alpha )

    Flux_X1_GR(iCF_S3) &
      = D * h * W**2 * Gm33 * V3 * ( V1 - Beta1 / Alpha )

    Flux_X1_GR(iCF_E)  &
      = D * W * ( h * W - One ) * ( V1 - Beta1 / Alpha ) + Beta1 / Alpha * P

    Flux_X1_GR(iCF_Ne) &
      = Ne * W * ( V1 - Beta1 / Alpha )

    RETURN
  END FUNCTION Flux_X1_GR


  PURE FUNCTION StressTensor_Diagonal( S_1, S_2, S_3, V_1, V_2, V_3, P )

    REAL(DP)             :: StressTensor_Diagonal(1:3)
    REAL(DP), INTENT(in) :: S_1, S_2, S_3, V_1, V_2, V_3, P

    StressTensor_Diagonal(1) = S_1 * V_1 + P
    StressTensor_Diagonal(2) = S_2 * V_2 + P
    StressTensor_Diagonal(3) = S_3 * V_3 + P

    RETURN
  END FUNCTION StressTensor_Diagonal


  REAL(DP) FUNCTION AlphaC_GR &
             ( U_L, U_R, F_L, F_R, Gmdd, LapseFunction, ShiftVector, aP, aM )

    ! --- Middle Wavespeed as suggested by Mignone and Bodo (2005) ---

    REAL(DP), INTENT(in) :: U_L(nCF), U_R(nCF), &
                            F_L(nCF), F_R(nCF), &
                            Gmdd, LapseFunction, ShiftVector, &
                            aP, aM

    REAL(DP) :: E_L, E_R, F_E_L, F_E_R, A_L, A_R, B_L, B_R
    REAL(DP) :: A, B, C, eps = SqrtTiny

    E_L   = U_L(iCF_E) + U_L(iCF_D)
    E_R   = U_R(iCF_E) + U_R(iCF_D)
    F_E_L = F_L(iCF_E) + F_L(iCF_D)
    F_E_R = F_R(iCF_E) + F_R(iCF_D)

    A_L = -aM * E_L - LapseFunction * F_E_L
    A_R = +aP * E_R - LapseFunction * F_E_R
    B_L = -aM * U_L(iCF_S1) - LapseFunction * F_L(iCF_S1)
    B_R = +aP * U_R(iCF_S1) - LapseFunction * F_R(iCF_S1)

    ! --- A, B, and C from quadratic equation ---
    A = Gmdd**2 &
          * ( A_R * ( -aM + ShiftVector ) - A_L * ( aP + ShiftVector ) )
    B = Gmdd &
          * ( ( LapseFunction * A_L - B_R * ( -aM + ShiftVector ) ) &
            - ( LapseFunction * A_R - B_L * (  aP + ShiftVector ) ) )
    C = LapseFunction * ( B_R - B_L )

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


  PURE FUNCTION NumericalFlux_X1_LLF_GR &
      ( u_L, u_R, Flux_L, Flux_R, aP, aM, aC, nF, &
        V1_L, V1_R, p_L, p_R, ShiftVector1, Gm11 )    

    ! --- Local Lax-Friedrichs Flux ---

    INTEGER,  INTENT(in)                   :: nF
    REAL(DP), DIMENSION(1:nF),  INTENT(in) :: u_L, u_R, Flux_L, Flux_R
    REAL(DP), INTENT(in)                   :: aP, aM, aC, V1_L, V1_R, &
                                              p_L, p_R, ShiftVector1, Gm11
    REAL(DP), DIMENSION(1:nF)              :: NumericalFlux_X1_LLF_GR
    REAL(DP) :: alpha

    alpha    = MAX( aM, aP )

    NumericalFlux_X1_LLF_GR &
      = 0.5_DP * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_X1_LLF_GR


  PURE FUNCTION NumericalFlux_X1_HLL_GR &
      ( u_L, u_R, Flux_L, Flux_R, aP, aM, aC, nF, &
        V1_L, V1_R, p_L, p_R, ShiftVector1, Gm11 )

    INTEGER,  INTENT(in)                   :: nF
    REAL(DP), DIMENSION(1:nF),  INTENT(in) :: u_L, u_R, Flux_L, Flux_R
    REAL(DP), INTENT(in)                   :: aP, aM, aC, V1_L, V1_R, &
                                              p_L, p_R, ShiftVector1, Gm11
    REAL(DP), DIMENSION(1:nF)              :: NumericalFlux_X1_HLL_GR

    NumericalFlux_X1_HLL_GR &
      = ( aP * flux_L + aM * flux_R &
            - aP * aM * ( u_R - u_L ) ) / ( aP + aM )

    RETURN
  END FUNCTION NumericalFlux_X1_HLL_GR


  PURE FUNCTION NumericalFlux_X1_HLLC_GR &
      ( u_L, u_R, Flux_L, Flux_R, aP, aM, aC, nF, &
        V1_L, V1_R, p_L, p_R, ShiftVector1, Gm11 )

    INTEGER,  INTENT(in)                  :: nF
    REAL(DP), DIMENSION(1:nF), INTENT(in) :: u_L, u_R, Flux_L, Flux_R
    REAL(DP), INTENT(in)                  :: aP, aM, aC, V1_L, V1_R, &
                                             p_L, p_R, ShiftVector1, Gm11
    REAL(DP)                              :: p, CF_D, S1, S2, S3, CF_E, Ne
    REAL(DP), DIMENSION(1:nF)             :: NumericalFlux_X1_HLLC_GR

    IF( aM .EQ. 0.0_DP )THEN

      NumericalFlux_X1_HLLC_GR = Flux_L

    ELSEIF( aP .EQ. 0.0_DP )THEN

      NumericalFlux_X1_HLLC_GR = Flux_R

    ELSE

      ! --- From Mignone & Bodo (2005)
      ! --- Note the sign change on aM which is due to it being
      ! --- read in as positive but the formulae assuming it is negative

      IF( aC .GE. 0.0_DP )THEN    

        ! -- UL_star

        p     = ( Gm11 * aC * ( u_L( iCF_E ) + u_L( iCF_D ) ) &
                * ( -aM + ShiftVector1 ) - u_L( iCF_S1 ) * ( aC  &
                - aM - V1_L + ShiftVector1 ) + p_L ) / ( 1.0_DP       &
                - Gm11 * aC * ( -aM + ShiftVector1 ) )

        CF_D  = u_L( iCF_D  ) *  ( -aM - V1_L    + ShiftVector1 ) &
                              /  ( -aM - aC + ShiftVector1 )

        S1    = u_L( iCF_S1 ) *  ( -aM - V1_L    + ShiftVector1 ) &
                              /  ( -aM - aC + ShiftVector1 ) &
                + ( p - p_L ) /  ( -aM - aC + ShiftVector1 )

        S2    = u_L( iCF_S2 ) *  ( -aM - V1_L    + ShiftVector1 ) &
                              /  ( -aM - aC + ShiftVector1 )

        S3    = u_L( iCF_S3 ) *  ( -aM - V1_L    + ShiftVector1 ) &
                              /  ( -aM - aC + ShiftVector1 )

        CF_E  = ( ( u_L( iCF_E ) + u_L( iCF_D ) ) * ( -aM + ShiftVector1 ) &
                + 1.0_DP / Gm11 * S1 - 1.0_DP / Gm11 * u_L( iCF_S1 ) ) &
                / ( -aM + ShiftVector1 )

        Ne    = u_L( iCF_Ne ) *  ( -aM - V1_L    + ShiftVector1 ) &
                              /  ( -aM - aC + ShiftVector1 )

      ELSE

        ! -- UR_star

        p     = ( Gm11 * aC * ( u_R( iCF_E ) + u_R( iCF_D ) ) &
                * ( aP + ShiftVector1 ) - u_R( iCF_S1 ) * ( aC   &
                + aP - V1_R + ShiftVector1 ) + p_R ) / ( 1.0_DP       &
                - Gm11 * aC * ( aP + ShiftVector1 ) )

        CF_D  = u_R( iCF_D  ) *  ( aP - V1_R    + ShiftVector1 ) &
                              /  ( aP - aC + ShiftVector1 )

        S1    = u_R( iCF_S1 ) *  ( aP - V1_R    + ShiftVector1 ) &
                              /  ( aP - aC + ShiftVector1 ) &
                + ( p - p_R ) /  ( aP - aC + ShiftVector1 )

        S2    = u_R( iCF_S2 ) *  ( aP - V1_R    + ShiftVector1 ) &
                              /  ( aP - aC + ShiftVector1 )

        S3    = u_R( iCF_S3 ) *  ( aP - V1_R    + ShiftVector1 ) &
                              /  ( aP - aC + ShiftVector1 )

        CF_E  = ( ( u_R( iCF_E ) + u_R( iCF_D ) ) * ( aP + ShiftVector1 )     &
                + 1.0_DP / Gm11 * S1 - 1.0_DP / Gm11 * u_R( iCF_S1 ) ) &
                / ( aP + ShiftVector1 )

        Ne     = u_R( iCF_Ne ) *  ( aP - V1_R    + ShiftVector1 ) &
                               /  ( aP - aC + ShiftVector1 )

      END IF

      NumericalFlux_X1_HLLC_GR( iCF_D  ) &
        = CF_D  * ( aC - ShiftVector1 )
      NumericalFlux_X1_HLLC_GR( iCF_S1 ) &
        = S1 * ( aC - ShiftVector1 ) + p
      NumericalFlux_X1_HLLC_GR( iCF_S2 ) &
        = S2 * ( aC - ShiftVector1 )
      NumericalFlux_X1_HLLC_GR( iCF_S3 ) &
        = S3 * ( aC - ShiftVector1 )
      NumericalFlux_X1_HLLC_GR( iCF_E  ) &
        = 1.0_DP / Gm11 * S1 - ShiftVector1 * ( CF_E - CF_D ) - aC * CF_D
      NumericalFlux_X1_HLLC_GR( iCF_Ne ) &
        = Ne * ( aC - ShiftVector1 )

    END IF

    RETURN
  END FUNCTION NumericalFlux_X1_HLLC_GR


END MODULE EulerEquationsUtilitiesModule_Beta_GR
