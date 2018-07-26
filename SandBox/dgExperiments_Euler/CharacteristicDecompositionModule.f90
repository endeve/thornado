MODULE CharacteristicDecompositionModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
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


  SUBROUTINE ComputeCharacteristicDecomposition( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

    INTEGER :: i
    REAL(DP), DIMENSION(1) :: D, V1, V2, V3, Gmdd11, Gmdd22, Gmdd33, Vu1, Vu2, Vu3, E, Ne, P, Cs
    REAL(DP), DIMENSION(1) :: Q, K, H, J, M1, M2, M3, Mu1, Mu2, Mu3, Vd1, Vd2, Vd3

    CALL ComputePrimitive &
           ( [ U(iCF_D ) ], [ U(iCF_S1) ], [ U(iCF_S2) ], &
             [ U(iCF_S3) ], [ U(iCF_E ) ], [ U(iCF_Ne) ], &
             D, V1, V2, V3, E, Ne, &
             [ G(iGF_Gm_dd_11) ], &
             [ G(iGF_Gm_dd_22) ], &
             [ G(iGF_Gm_dd_33) ] )

    CALL ComputePressureFromPrimitive &
           ( D, E, Ne, P )

    CALL ComputeSoundSpeedFromPrimitive &
           ( D, E, Ne, Cs )
    
    ! --- Variables for new eigendecomposition. ---
 
    Gmdd11 = G(iGF_Gm_dd_11)
    Gmdd22 = G(iGF_Gm_dd_22)
    Gmdd33 = G(iGF_Gm_dd_33)

    Vd1 = V1
    Vd2 = V2
    Vd3 = V3

    Vu1 = V1 / Gmdd11
    Vu2 = V2 / Gmdd22
    Vu3 = V3 / Gmdd33

    Q = P / ( E * Cs)
    K = Half * (Vu1 * Vd1 + Vu2 * Vd2 + Vu3 * Vd3)
    H = Q * K / Cs
    J = ( E / P ) * Cs**2 + K
    Mu1 = Vu1 / Cs
    Mu2 = Vu2 / Cs
    Mu3 = Vu3 / Cs

    ! --- Variables for old eigendecomposition ---

    !Q  = P / ( E * Cs )
    !K  = Half * ( V1**2 + V2**2 + V3**2 )
    !H  = Q * K / Cs
    !J  = ( E / P ) * Cs**2 + K
    !M1 = V1 / Cs
    !M2 = V2 / Cs
    !M3 = V3 / Cs

    SELECT CASE( iDim )

      CASE( 1 )

        ! --- New Eigendecomposition ---

        R(:,1) = [  One, Vd1 - Cs * SQRT(Gmdd11),     Vd2,   Vd3, J - Cs * SQRT(Gmdd11) * Vu1, Zero ]
        R(:,2) = [  One,                     Vd1,     Vd2,   Vd3,                           K, Zero ]
        R(:,3) = [  One, Vd1 + Cs * SQRT(Gmdd11),     Vd2,   Vd3, J + Cs * SQRT(Gmdd11) * Vu1, Zero ]
        R(:,4) = [ Zero,                    Zero,   - One,  Zero,                       - Vu2, Zero ]
        R(:,5) = [ Zero,                    Zero,    Zero,   One,                         Vu3, Zero ]
        R(:,6) = [ Zero,                    Zero,    Zero,  Zero,                        Zero,  One ]

        invR(:,1) &
          = [ Half * ( H + SQRT(Gmdd11) * Mu1 ), &
              One - H, &
              Half * ( H - SQRT(Gmdd11) * Mu1 ), &
              Vd2, - Vd3, Zero ]

        invR(:,2) &
          = [ - Half * ( Q * Mu1 + One / ( SQRT(Gmdd11) * Cs ) ), &
                         Q * Mu1, &
              - Half * ( Q * Mu1 - One / ( SQRT(Gmdd11) * Cs ) ), &
              Zero, Zero, Zero ]

        invR(:,3) &
          = [ - Half * Q * Mu2, &
                       Q * Mu2, &
              - Half * Q * Mu2, &
              - One, Zero, Zero ]

        invR(:,4) &
          = [ - Half * Q * Mu3, &
                       Q * Mu3, &
              - Half * Q * Mu3, &
              Zero, One, Zero ]

        invR(:,5) &
          = [ Half * Q / Cs, &
                   - Q / Cs, &
              Half * Q / Cs, &
              Zero, Zero, Zero ]

        invR(:,6) &
          = [ Zero, Zero, Zero, Zero, Zero, One ]

        ! -- Old Eigendecomp ---

        !R(:,1) = [  One, V1 - Cs,     V2,   V3, J - Cs * V1, Zero ]
        !R(:,2) = [  One,      V1,     V2,   V3,           K, Zero ]
        !R(:,3) = [  One, V1 + Cs,     V2,   V3, J + Cs * V1, Zero ]
        !R(:,4) = [ Zero,    Zero,  - One, Zero,        - V2, Zero ]
        !R(:,5) = [ Zero,    Zero,   Zero,  One,          V3, Zero ]
        !R(:,6) = [ Zero,    Zero,   Zero, Zero,        Zero,  One ]

        !invR(:,1) &
        !  = [ Half * ( H + M1 ), &
        !      One - H, &
        !      Half * ( H - M1 ), &
        !      V2, - V3, Zero ]

        !invR(:,2) &
        !  = [ - Half * ( Q * M1 + One / Cs ), &
        !                 Q * M1, &
        !      - Half * ( Q * M1 - One / Cs ), &
        !      Zero, Zero, Zero ]

        !invR(:,3) &
        !  = [ - Half * Q * M2, &
        !               Q * M2, &
        !      - Half * Q * M2, &
        !      - One, Zero, Zero ]

        !invR(:,4) &
        !  = [ - Half * Q * M3, &
        !               Q * M3, &
        !      - Half * Q * M3, &
        !     Zero, One, Zero ]

        !invR(:,5) &
        !  = [ Half * Q / Cs, &
        !           - Q / Cs, &
        !      Half * Q / Cs, &
        !      Zero, Zero, Zero ]

        !invR(:,6) &
        !  = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X1):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

      CASE( 2 )

        ! --- New Eigendecomp ---

        R(:,1) = [  One,    Vd1, Vd2 - Cs * SQRT(Gmdd22),    Vd3, J - Cs * SQRT(Gmdd22) * Vu2, Zero ]
        R(:,2) = [  One,    Vd1,                     Vd2,    Vd3,                           K, Zero ]
        R(:,3) = [  One,    Vd1, Vd2 + Cs * SQRT(Gmdd22),    Vd3, J + Cs * SQRT(Gmdd22) * Vu2, Zero ]
        R(:,4) = [ Zero,  - One,                    Zero,   Zero,                       - Vu1, Zero ]
        R(:,5) = [ Zero,   Zero,                    Zero,    One,                         Vu3, Zero ]
        R(:,6) = [ Zero,   Zero,                    Zero,   Zero,                        Zero,  One ]

        invR(:,1) &
          = [ Half * ( H + SQRT(Gmdd22) * Mu2 ), &
              One - H, &
              Half * ( H - SQRT(Gmdd22) * Mu2 ), &
                Vd1, - Vd3, Zero ]

        invR(:,2) &
          = [ - Half * Q * Mu1, &
                       Q * Mu1, &
              - Half * Q * Mu1, &
              - One, Zero, Zero ]

        invR(:,3) &
          = [ - Half * ( Q * Mu2 + One / ( SQRT(Gmdd22) * Cs ) ), &
                         Q * Mu2, &
              - Half * ( Q * Mu2 - One / ( SQRT(Gmdd22) * Cs ) ), &
              Zero, Zero, Zero ]

        invR(:,4) &
          = [ - Half * Q * Mu3, &
                       Q * Mu3, &
              - Half * Q * Mu3, &
              Zero, One, Zero ]

        invR(:,5) &
          = [ Half * Q / Cs, &
                   - Q / Cs, &
              Half * Q / Cs, &
              Zero, Zero, Zero ]

        invR(:,6) &
          = [ Zero, Zero, Zero, Zero, Zero, One ]

        ! --- Old Eigendecomp ---

        !R(:,1) = [ One,    V1, V2 - Cs,    V3, J - Cs * V2, Zero ]
        !R(:,2) = [ One,    V1,      V2,    V3,           K, Zero ]
        !R(:,3) = [ One,    V1, V2 + Cs,    V3, J + Cs * V2, Zero ]
        !R(:,4) = [ Zero,  One,    Zero,  Zero,          V1, Zero ]
        !R(:,5) = [ Zero, Zero,    Zero, - One,        - V3, Zero ]
        !R(:,6) = [ Zero, Zero,    Zero,  Zero,        Zero,  One ]

        !invR(:,1) &
        !  = [ Half * ( H + M2 ), &
        !      One - H, &
        !      Half * ( H - M2 ), &
        !      - V1, V3, Zero ]

        !invR(:,2) &
        !  = [ - Half * Q * M1, &
        !               Q * M1, &
        !      - Half * Q * M1, &
        !      One, Zero, Zero ]

        !invR(:,3) &
        !  = [ - Half * ( Q * M2 + One / Cs ), &
        !                 Q * M2, &
        !      - Half * ( Q * M2 - One / Cs ), &
        !      Zero, Zero, Zero ]

        !invR(:,4) &
        !  = [ - Half * Q * M3, &
        !               Q * M3, &
        !      - Half * Q * M3, &
        !      Zero, - One, Zero ]

        !invR(:,5) &
        !  = [ Half * Q / Cs, &
        !           - Q / Cs, &
        !      Half * Q / Cs, &
        !      Zero, Zero, Zero ]

        !invR(:,6) &
        !  = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X2):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

      CASE( 3 )

        R(:,1) = [  One,   Vd1,     Vd2,   Vd3 - Cs * SQRT(Gmdd33), J - Cs * SQRT(Gmdd33) * Vu3, Zero ]
        R(:,2) = [  One,   Vd1,     Vd2,                       Vd3,                           K, Zero ]
        R(:,3) = [  One,   Vd1,     Vd2,   Vd3 + Cs * SQRT(Gmdd33), J + Cs * SQRT(Gmdd33) * Vu3, Zero ]
        R(:,4) = [ Zero, - One,    Zero,                      Zero,                       - Vu1, Zero ]
        R(:,5) = [ Zero,  Zero,     One,                      Zero,                         Vu2, Zero ]
        R(:,6) = [ Zero,  Zero,    Zero,                      Zero,                        Zero,  One ]

        invR(:,1) &
          = [ Half * ( H + SQRT(Gmdd33) * Mu3 ), &
              One - H, &
              Half * ( H - SQRT(Gmdd33) * Mu3 ), &
              Vd1, - Vd2, Zero ]

        invR(:,2) &
          = [ - Half * Q * Mu1, &
                       Q * Mu1, &
              - Half * Q * Mu1, &
              - One, Zero, Zero ]

        invR(:,3) &
          = [ - Half * Q * Mu2, &
                       Q * Mu2, &
              - Half * Q * Mu2, &
                Zero, One, Zero ]

        invR(:,4) &
          = [ - Half * ( Q * Mu3 + One / ( SQRT(Gmdd33) * Cs  ) ), &
                         Q * Mu3, &
              - Half * ( Q * Mu3 + One / ( SQRT(Gmdd33) * Cs  ) ), &
               Zero, Zero, Zero ]

        invR(:,5) &
          = [ Half * Q / Cs, &
                   - Q / Cs, &
              Half * Q / Cs, &
              Zero, Zero, Zero ]

        invR(:,6) &
          = [ Zero, Zero, Zero, Zero, Zero, One ]

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X3):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

    END SELECT

  END SUBROUTINE ComputeCharacteristicDecomposition


END MODULE CharacteristicDecompositionModule
