PROGRAM ComputePrimitiveTest

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState, &
    ComputePressureFromSpecificInternalEnergy
  USE EulerEquationsUtilitiesModule_Beta_GR, ONLY: &
    ComputePrimitive_GR, &
    ComputeFunJacP, &
    ComputePressureWithBisectionMethod

  IMPLICIT NONE

  INTEGER,  PARAMETER :: nNodes = 1
  INTEGER             :: i
  REAL(DP)            :: Pressure(nNodes), q, SSq, Pbisec, Pbrent
  REAL(DP)            :: U(nNodes,nCF), P(nNodes,nPF), G(nNodes,nGF)
  REAL(DP)            :: Pmin, Pmax, DeltaP
  REAL(DP), PARAMETER :: Gamma_IDEAL = 4.0_DP / 3.0_DP

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma_IDEAL )

  i = 1

  ! --- MAKE SURE TO CHANGE 'E' TO 'd' ---
  U(i,iCF_D)  =  1.0000577410694345d+000
  U(i,iCF_S1) = -2.0044113494280685d-005
  U(i,iCF_S2) =  0.0000000000000000d+000
  U(i,iCF_S3) =  0.0000000000000000d+000
  U(i,iCF_E)  =  5.2867016077136483d-009

  G(i,iGF_Gm_dd_11) =  1.0000000000000000d+000
  G(i,iGF_Gm_dd_22) =  1.0000000000000000d+000
  G(i,iGF_Gm_dd_33) =  1.0000000000000000d+000

  SSq =  U(i,iCF_S1)**2 / G(i,iGF_Gm_dd_11) &
       + U(i,iCF_S2)**2 / G(i,iGF_Gm_dd_22) &
       + U(i,iCF_S3)**2 / G(i,iGF_Gm_dd_33)

  q = U(i,iCF_D) + U(i,iCF_E) - SQRT( U(i,iCF_D)**2 + SSq )

  Pmin = MAX( -( U(i,iCF_D) + U(i,iCF_E) ) + SQRT( SSq ), SqrtTiny )
  Pmax = 2.0_DP * ( One - One / Gamma_IDEAL ) * U(i,iCF_E)
  
  WRITE(*,'(A,ES24.16E3)') 'Pmin   = ', Pmin
  WRITE(*,'(A,ES24.16E3)') 'Pmax   = ', Pmax

  CALL ComputePressureWithBisectionMethod( U(i,iCF_D), U(i,iCF_E), SSq, Pbisec )
  WRITE(*,'(A,ES24.16E3)') 'Pbisec = ', Pbisec

  CALL FindRootBrent( U(i,iCF_D), U(i,iCF_E), SSq, Pmin, Pmax, Pbrent )
  WRITE(*,'(A,ES24.16E3)') 'Pbrent = ', Pbrent
  WRITE(*,*)

  WRITE(*,'(A)') 'Creating FunP array...'
  CALL CreateFunParray( U(i,iCF_D), U(i,iCF_E), SSq, Pmin, Pmax, Pbisec )
  STOP 'Created FunP array'
  WRITE(*,'(A)') 'Created FunP array!'
  WRITE(*,*)

  WRITE(*,'(A)') 'CALL ComputePrimitive_GR:'
  WRITE(*,'(A)') '-------------------------'
  CALL ComputePrimitive_GR &
         ( U(:,iCF_D), &
           U(:,iCF_S1), U(:,iCF_S2), U(:,iCF_S3), &
           U(:,iCF_E), U(:,iCF_Ne), &
           P(:,iPF_D), &
           P(:,iPF_V1), P(:,iPF_V2), P(:,iPF_V3), &
           P(:,iPF_E), P(:,iPF_Ne), &
           Pressure, &
           G(:,iGF_Gm_dd_11), G(:,iGF_Gm_dd_22), G(:,iGF_Gm_dd_33) )

  WRITE(*,'(A9,ES24.16E3)') "P(out) = ", Pressure(i)
  
  CALL FinalizeEquationOfState


CONTAINS
  

    SUBROUTINE CreateFunParray( CF_D, CF_E, SSq, Pmin, Pmax, Pbisec )

      REAL(DP), INTENT(in) :: CF_D, CF_E, SSq, Pmin, Pmax, Pbisec
      INTEGER, PARAMETER   :: N = 10000
      REAL(DP)             :: Parr(N), FunP(N), JacP(N), DeltaP
      INTEGER              :: j

      OPEN( 100, FILE = 'FunP.dat' )

      DeltaP = ( Pmax - Pmin ) / N
      WRITE( 100, '(A)' ) '# P, FunP, JacP'
      WRITE( 100, '(ES24.16E3,A,A)') Pbisec, ' nan ', ' nan '
      DO j = 1, N
        Parr(j) = Pmin + (j-1) * DeltaP
        CALL ComputeFunJacP &
               ( CF_D, CF_E, SSq, Parr(j), FunP(j), JacP(j), Pbisec )
        WRITE( 100, '(ES24.16E3,1x,ES24.16E3,1x,ES24.16E3)') &
          Parr(j), FunP(j), JacP(j)
      END DO

      CLOSE( 100 )

    END SUBROUTINE CreateFunParray


  SUBROUTINE FindRootBrent( CF_D, CF_E, SSq, Pmin, Pmax, P )

    REAL(DP), INTENT(in)  :: CF_D, CF_E, SSq, Pmin, Pmax
    REAL(DP), INTENT(out) :: P

    REAL(DP) :: PA, PB, PC, PD, PS, PSwap
    REAL(DP) :: FunPA, FunPB, FunPC, FunPD, FunPS
    REAL(DP) :: JacPA, JacPB, JacPC, JacPD, JacPS

    REAL(DP), PARAMETER :: TolP = 1.0d-12, TolFunP = 1.0d-10
    LOGICAL :: mflag, COND1, COND2, COND3, COND4, COND5

    ! --- Define values that bracket the root ---
    PA = Pmin
    PB = Pmax
    CALL ComputeFunJacP( CF_D, CF_E, SSq, PA, FunPA, JacPA, Pbisec )
    CALL ComputeFunJacP( CF_D, CF_E, SSq, PB, FunPB, JacPB, Pbisec )

    ! --- Exit if solution is not bracketed ---
    IF( .NOT. FunPA * FunPB .LT. 0.0_DP )THEN
      WRITE(*,'(A)') 'No root in interval. Exiting...'
      STOP
    END IF

    ! --- Swap a and b ---
    IF( ABS( FunPA ) .LT. ABS( FunPB ) )THEN
       PSwap = PB
       PB    = PA
       PA    = PSwap
       CALL ComputeFunJacP( CF_D, CF_E, SSq, PA, FunPA, JacPA, Pbisec )
       CALL ComputeFunJacP( CF_D, CF_E, SSq, PB, FunPB, JacPB, Pbisec )
    END IF

    ! --- Set PC = PA and compute FunPC ---
    PC = PA
    CALL ComputeFunJacP( CF_D, CF_E, SSq, PC, FunPC, JacPC, Pbisec )

    ! --- Set mflag ---
    mflag = .TRUE.

    DO WHILE( ( ABS( PB - PA ) .GE. TolP ) .AND. ( ABS(FunPB) .GE. TolFunP ) )

      IF( ( FunPA .NE. FunPC ) .AND. ( FunPB .NE. FunPC ) )THEN
        PS =   PA * FunPB * FunPC / ( ( FunPA - FunPB ) * ( FunPA - FunPC ) ) &
             + PB * FunPA * FunPC / ( ( FunPB - FunPA ) * ( FunPB - FunPC ) ) &
             + PC * FunPA * FunPB / ( ( FunPC - FunPA ) * ( FunPC - FunPB ) )
      ELSE
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

      CALL ComputeFunJacP( CF_D, CF_E, SSq, PS, FunPS, JacPS, Pbisec )
      PD = PC
      PC = PB

      IF( FunPA * FunPS .LT. 0.0_DP )THEN
        PB = PS
        CALL ComputeFunJacP( CF_D, CF_E, SSq, PB, FunPB, JacPB, Pbisec )
      ELSE
        PA = PS
        CALL ComputeFunJacP( CF_D, CF_E, SSq, PA, FunPA, JacPA, Pbisec )
      END IF

      ! --- Swap a and b ---
      IF( ABS( FunPA ) .LT. ABS( FunPB ) )THEN
        PSwap = PB
        PB    = PA
        PA    = PSwap
        CALL ComputeFunJacP( CF_D, CF_E, SSq, PA, FunPA, JacPA, Pbisec )
        CALL ComputeFunJacP( CF_D, CF_E, SSq, PB, FunPB, JacPB, Pbisec )
      END IF

      P = PB
    END DO

    RETURN
  END SUBROUTINE FindRootBrent

END PROGRAM ComputePrimitiveTest
