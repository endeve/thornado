PROGRAM Driver

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, One
  USE RadiationFieldsModule, ONLY: &
    nCR
  USE CharacteristicDecompositionModule, ONLY: &
    ComputeCharacteristicDecomposition_RC
  
  IMPLICIT NONE

  INTEGER  :: iDim
  REAL(DP) :: G(3)
  REAL(DP) :: U(nCR)
  REAL(DP) :: R(nCR,nCR)
  REAL(DP) :: invR(nCR,nCR)
  REAL(DP) :: TestMatrix(nCR,nCR)

  REAL(DP) :: D, I_1, I_2, I_3, Gm_dd_11, Gm_dd_22, Gm_dd_33

  WRITE(*,*)
  WRITE(*,'(A4,A)')'', 'Running CharacteristicDecomposition Test'


  ! ComputeFluxJacobian test in x:

  iDim = 1
  D        = 5.d0 * SqrtTiny!5.0d-1
  I_1      = SqrtTiny!1.0d-1 
  I_2      = Zero
  I_3      = Zero
  Gm_dd_11 = One
  Gm_dd_22 = Zero
  Gm_dd_33 = Zero

  G = [ Gm_dd_11, Gm_dd_22, Gm_dd_33 ]
  U = [ D, I_1, I_2, I_3 ]

  WRITE(*,*)
  WRITE(*,'(A4,4A18)') '', ' J', 'H1', 'H2', 'H3'
  WRITE(*,'(A4,4ES18.3)') '', D, I_1, I_2, I_3

  CALL ComputeCharacteristicDecomposition_RC &
         ( iDim, G, U, R, invR )

  ! ComputeFluxJacobian test in y:
  
  iDim = 2
  D        = 5.d0 * SqrtTiny!5.0d-1
  I_1      = Zero
  I_2      = SqrtTiny
  I_3      = Zero
  Gm_dd_11 = Zero
  Gm_dd_22 = One
  Gm_dd_33 = Zero

  G = [ Gm_dd_11, Gm_dd_22, Gm_dd_33]
  U = [D, I_1, I_2, I_3]

  WRITE(*,*)
  WRITE(*,'(A4,4A18)') '', ' J', 'H1', 'H2', 'H3'
  WRITE(*,'(A4,4ES18.3)') '', D, I_1, I_2, I_3

  CALL ComputeCharacteristicDecomposition_RC &
         ( iDim, G, U, R, invR )

  ! ComputeFluxJacobian test in z:
  
  iDim = 3
  D        = 5.d0 * SqrtTiny!5.0d-1
  I_1      = Zero
  I_2      = Zero
  I_3      = SqrtTiny
  Gm_dd_11 = Zero
  Gm_dd_22 = Zero
  Gm_dd_33 = One

  G = [ Gm_dd_11, Gm_dd_22, Gm_dd_33 ]
  U = [ D, I_1, I_2, I_3 ]

  WRITE(*,*)
  WRITE(*,'(A4,4A18)') '', ' J', 'H1', 'H2', 'H3'
  WRITE(*,'(A4,4ES18.3)') '', D, I_1, I_2, I_3

  CALL ComputeCharacteristicDecomposition_RC &
         ( iDim, G, U, R, invR )

END PROGRAM Driver
