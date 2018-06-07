MODULE CharacteristicDecompositionModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne
  USE EulerEquationsUtilitiesModule_Beta, ONLY: &
    ComputePrimitive
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: Debug = .FALSE.

  PUBLIC :: ComputeCharacteristicDecomposition

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition( iDim, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

    INTEGER :: i
    REAL(DP), DIMENSION(1) :: D, V1, V2, V3, E, Ne, P, Cs
    REAL(DP), DIMENSION(1) :: G, K, H, J, M1, M2, M3

    CALL ComputePrimitive &
           ( [U(iCF_D )], [U(iCF_S1)], [U(iCF_S2)], &
             [U(iCF_S3)], [U(iCF_E )], [U(iCF_Ne)], &
             D, V1, V2, V3, E, Ne, [One], [One], [One] )

    CALL ComputePressureFromPrimitive &
           ( D, E, Ne, P )

    CALL ComputeSoundSpeedFromPrimitive &
           ( D, E, Ne, Cs )

    G  = P / ( E * Cs )
    K  = Half * ( V1**2 + V2**2 + V3**2 )
    H  = G * K / Cs
    J  = ( E / P ) * Cs**2 + K
    M1 = V1 / Cs
    M2 = V2 / Cs
    M3 = v3 / Cs

    SELECT CASE( iDim )

      CASE( 1 )

        R(:,1) = [  One, V1 - Cs,     V2,   V3, J - Cs * V1, Zero ]
        R(:,2) = [  One,      V1,     V2,   V3,           K, Zero ]
        R(:,3) = [  One, V1 + Cs,     V2,   V3, J + Cs * V1, Zero ]
        R(:,4) = [ Zero,    Zero,  - One, Zero,        - V2, Zero ]
        R(:,5) = [ Zero,    Zero,   Zero,  One,          V3, Zero ]
        R(:,6) = [ Zero,    Zero,   Zero, Zero,        Zero,  One ]

        invR(:,1) &
          = [ Half * ( H + M1 ), &
              One - H, &
              Half * ( H - M1 ), &
              V2, - V3, Zero ]

        invR(:,2) &
          = [ - Half * ( G * M1 + One / Cs ), &
                         G * M1, &
              - Half * ( G * M1 - One / Cs ), &
              Zero, Zero, Zero ]

        invR(:,3) &
          = [ - Half * G * M2, &
                       G * M2, &
              - Half * G * M2, &
              - One, Zero, Zero ]

        invR(:,4) &
          = [ - Half * G * M3, &
                       G * M3, &
              - Half * G * M3, &
              Zero, One, Zero ]

        invR(:,5) &
          = [ Half * G / Cs, &
                   - G / Cs, &
              Half * G / Cs, &
              Zero, Zero, Zero ]

        invR(:,6) &
          = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X1):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

      CASE( 2 )

        R(:,1) = [ One,    V1, V2 - Cs,    V3, J - Cs * V2, Zero ]
        R(:,2) = [ One,    V1,      V2,    V3,           K, Zero ]
        R(:,3) = [ One,    V1, V2 + Cs,    V3, J + Cs * V2, Zero ]
        R(:,4) = [ Zero,  One,    Zero,  Zero,          V1, Zero ]
        R(:,5) = [ Zero, Zero,    Zero, - One,        - V3, Zero ]
        R(:,6) = [ Zero, Zero,    Zero,  Zero,        Zero,  One ]

        invR(:,1) &
          = [ Half * ( H + M2 ), &
              One - H, &
              Half * ( H - M2 ), &
              - V1, V3, Zero ]

        invR(:,2) &
          = [ - Half * G * M1, &
                       G * M1, &
              - Half * G * M1, &
              One, Zero, Zero ]

        invR(:,3) &
          = [ - Half * ( G * M2 + One / Cs ), &
                         G * M2, &
              - Half * ( G * M2 - One / Cs ), &
              Zero, Zero, Zero ]

        invR(:,4) &
          = [ - Half * G * M3, &
                       G * M3, &
              - Half * G * M3, &
              Zero, - One, Zero ]

        invR(:,5) &
          = [ Half * G / Cs, &
                   - G / Cs, &
              Half * G / Cs, &
              Zero, Zero, Zero ]

        invR(:,6) &
          = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X2):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

      CASE( 3 )

    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition


END MODULE CharacteristicDecompositionModule
