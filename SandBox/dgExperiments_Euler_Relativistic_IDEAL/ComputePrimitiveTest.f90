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
    ComputePressureFromPrimitive
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputePrimitive_Euler_Relativistic, &
    ComputeFunJacP, &
    ComputePressureWithBisectionMethod, &
    ComputePressureWithBrentsMethod

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
  U(i,iCF_D)  =  9.9994521915916246d-001
  U(i,iCF_S1) = -1.4111396032719828d-005
  U(i,iCF_S2) =  0.0000000000000000d+000
  U(i,iCF_S3) =  0.0000000000000000d+000
  U(i,iCF_E)  =  1.7191523746929956d-007

  G(i,iGF_Gm_dd_11) =  1.0000000000000000d+000
  G(i,iGF_Gm_dd_22) =  1.9775390625000015d-004
  G(i,iGF_Gm_dd_33) =  1.9775390625000015d-004

  SSq =  U(i,iCF_S1)**2 / G(i,iGF_Gm_dd_11) &
       + U(i,iCF_S2)**2 / G(i,iGF_Gm_dd_22) &
       + U(i,iCF_S3)**2 / G(i,iGF_Gm_dd_33)

  q = U(i,iCF_D) + U(i,iCF_E) - SQRT( U(i,iCF_D)**2 + SSq )

  Pmin = MAX( -( U(i,iCF_D) + U(i,iCF_E) ) + SQRT( SSq ), SqrtTiny )
  Pmax = 2.0_DP * ( One - One / Gamma_IDEAL ) * U(i,iCF_E)
  
  WRITE(*,'(A,ES24.16E3)') 'Pmin   = ', Pmin
  WRITE(*,'(A,ES24.16E3)') 'Pmax   = ', Pmax

  CALL ComputePressureWithBisectionMethod &
    ( U(i,iCF_D), U(i,iCF_E), SSq, Pbisec )
  WRITE(*,'(A,ES24.16E3)') 'Pbisec = ', Pbisec

  CALL ComputePressureWithBrentsMethod &
    ( U(i,iCF_D), U(i,iCF_E), SSq, Pbrent )
  WRITE(*,'(A,ES24.16E3)') 'Pbrent = ', Pbrent
  WRITE(*,*)

  WRITE(*,'(A)') 'Creating FunP array...'
  CALL CreateFunParray( U(i,iCF_D), U(i,iCF_E), SSq, Pmin, Pmax, Pbisec )
  !STOP 'Created FunP array'
  WRITE(*,'(A)') 'Created FunP array!'
  WRITE(*,*)

  WRITE(*,'(A)') 'CALL ComputePrimitive_Euler_Relativistic:'
  WRITE(*,'(A)') '-----------------------------------------'
  CALL ComputePrimitive_Euler_Relativistic &
         ( U(:,iCF_D), &
           U(:,iCF_S1), U(:,iCF_S2), U(:,iCF_S3), &
           U(:,iCF_E), U(:,iCF_Ne), &
           P(:,iPF_D), &
           P(:,iPF_V1), P(:,iPF_V2), P(:,iPF_V3), &
           P(:,iPF_E), P(:,iPF_Ne), &
           G(:,iGF_Gm_dd_11), G(:,iGF_Gm_dd_22), G(:,iGF_Gm_dd_33) )
  CALL ComputePressureFromPrimitive &
         ( P(:,iPF_D), P(:,iPF_E), P(:,iPF_Ne), Pressure )

  WRITE(*,'(A9,ES24.16E3)') "P(out) = ", Pressure(i)
  
  CALL FinalizeEquationOfState


CONTAINS
  

    SUBROUTINE CreateFunParray( CF_D, CF_E, SSq, Pmin, Pmax, Pbisec )

      REAL(DP), INTENT(in) :: CF_D, CF_E, SSq, Pmin, Pmax, Pbisec
      INTEGER, PARAMETER   :: N = 1000000
      REAL(DP)             :: Parr(N), FunP(N), JacP(N), DeltaP
      INTEGER              :: j

      OPEN( 100, FILE = 'FunP.dat' )

      DeltaP = ( Pmax - Pmin ) / N
      WRITE( 100, '(A)' ) '# P, FunP, JacP'
      WRITE( 100, '(ES24.16E3,A,A)') Pbisec, ' nan ', ' nan '
      DO j = 1, N
        Parr(j) = Pmin + (j-1) * DeltaP
        CALL ComputeFunJacP &
               ( CF_D, CF_E, SSq, Parr(j), FunP(j), JacP(j) )
        WRITE( 100, '(ES24.16E3,1x,ES24.16E3,1x,ES24.16E3)') &
          Parr(j), FunP(j), JacP(j)
      END DO

      CLOSE( 100 )

    END SUBROUTINE CreateFunParray


END PROGRAM ComputePrimitiveTest
