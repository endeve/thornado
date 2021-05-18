MODULE Euler_CharacteristicDecompositionModule_NonRelativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule, ONLY: &
    ComputePressureFromPrimitive, &
    ComputeSoundSpeedFromPrimitive
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_SL_CharDecomp

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER :: Debug = .FALSE.

  PUBLIC :: ComputeCharacteristicDecomposition_Euler_NonRelativistic_IDEAL

CONTAINS


  SUBROUTINE ComputeCharacteristicDecomposition_Euler_NonRelativistic_IDEAL &
    ( iDim, G, U, R, invR )

    INTEGER,  INTENT(in)  :: iDim
    REAL(DP), INTENT(in)  :: G(nGF)
    REAL(DP), INTENT(in)  :: U(nCF)
    REAL(DP), INTENT(out) :: R(nCF,nCF)
    REAL(DP), INTENT(out) :: invR(nCF,nCF)

    INTEGER :: i
    REAL(DP), DIMENSION(1) :: D, V1, V2, V3, Gmdd11, Gmdd22, Gmdd33, Vu1, Vu2, Vu3, E, Ne, P, Cs
    REAL(DP), DIMENSION(1) :: Q, K, H, J, Mu1, Mu2, Mu3, Vd1, Vd2, Vd3

    CALL TimersStart_Euler( Timer_Euler_SL_CharDecomp )

    CALL ComputePrimitive_Euler_NonRelativistic &
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

    Gmdd11 = G(iGF_Gm_dd_11)
    Gmdd22 = G(iGF_Gm_dd_22)
    Gmdd33 = G(iGF_Gm_dd_33)

    Vu1 = V1
    Vu2 = V2
    Vu3 = V3

    Vd1 = Gmdd11 * Vu1
    Vd2 = Gmdd22 * Vu2
    Vd3 = Gmdd33 * Vu3

    Q = P / ( E * Cs)
    K = Half * (Vu1 * Vd1 + Vu2 * Vd2 + Vu3 * Vd3)
    H = Q * K / Cs
    J = ( E / P ) * Cs**2 + K
    Mu1 = Vu1 / Cs
    Mu2 = Vu2 / Cs
    Mu3 = Vu3 / Cs

    SELECT CASE( iDim )

      CASE( 1 )

        R(:,1) = [  One, Vd1 - Cs * SQRT( Gmdd11 ),     Vd2,   Vd3, J - Cs * SQRT( Gmdd11 ) * Vu1, Zero ]
        R(:,2) = [  One,                       Vd1,     Vd2,   Vd3,                             K, Zero ]
        R(:,3) = [  One, Vd1 + Cs * SQRT( Gmdd11 ),     Vd2,   Vd3, J + Cs * SQRT( Gmdd11 ) * Vu1, Zero ]
        R(:,4) = [ Zero,                      Zero,   - One,  Zero,                       -   Vu2, Zero ]
        R(:,5) = [ Zero,                      Zero,    Zero,   One,                           Vu3, Zero ]
        R(:,6) = [ Zero,                      Zero,    Zero,  Zero,                          Zero,  One ]

        invR(:,1) &
          = [ Half * ( H + SQRT( Gmdd11 ) * Mu1 ), &
              One - H, &
              Half * ( H - SQRT( Gmdd11 ) * Mu1 ), &
              Vd2, - Vd3, Zero ]

        invR(:,2) &
          = [ - Half * ( Q * Mu1 + One / ( SQRT( Gmdd11 ) * Cs ) ), &
                         Q * Mu1, &
              - Half * ( Q * Mu1 - One / ( SQRT( Gmdd11 ) * Cs ) ), &
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

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X1):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

      CASE( 2 )

        R(:,1) = [  One,    Vd1, Vd2 - Cs * SQRT( Gmdd22 ),    Vd3, J - Cs * SQRT( Gmdd22 ) * Vu2, Zero ]
        R(:,2) = [  One,    Vd1,                       Vd2,    Vd3,                             K, Zero ]
        R(:,3) = [  One,    Vd1, Vd2 + Cs * SQRT( Gmdd22 ),    Vd3, J + Cs * SQRT( Gmdd22 ) * Vu2, Zero ]
        R(:,4) = [ Zero,  - One,                      Zero,   Zero,                         - Vu1, Zero ]
        R(:,5) = [ Zero,   Zero,                      Zero,    One,                           Vu3, Zero ]
        R(:,6) = [ Zero,   Zero,                      Zero,   Zero,                          Zero,  One ]

        invR(:,1) &
          = [ Half * ( H + SQRT( Gmdd22 ) * Mu2 ), &
              One - H, &
              Half * ( H - SQRT( Gmdd22 ) * Mu2 ), &
                Vd1, - Vd3, Zero ]

        invR(:,2) &
          = [ - Half * Q * Mu1, &
                       Q * Mu1, &
              - Half * Q * Mu1, &
              - One, Zero, Zero ]

        invR(:,3) &
          = [ - Half * ( Q * Mu2 + One / ( SQRT( Gmdd22 ) * Cs ) ), &
                         Q * Mu2, &
              - Half * ( Q * Mu2 - One / ( SQRT( Gmdd22 ) * Cs ) ), &
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

        IF( Debug )THEN

          WRITE(*,*)
          WRITE(*,'(A4,A)') '', 'invR * R (X2):'
          WRITE(*,*)
          DO i = 1, 6
            WRITE(*,'(A4,6ES16.7E2)') '', MATMUL( invR(i,:), R(:,:) )
          END DO

        END IF

      CASE( 3 )

        R(:,1) = [  One,   Vd1,     Vd2,   Vd3 - Cs * SQRT( Gmdd33 ), J - Cs * SQRT( Gmdd33 ) * Vu3, Zero ]
        R(:,2) = [  One,   Vd1,     Vd2,                         Vd3,                             K, Zero ]
        R(:,3) = [  One,   Vd1,     Vd2,   Vd3 + Cs * SQRT( Gmdd33 ), J + Cs * SQRT( Gmdd33 ) * Vu3, Zero ]
        R(:,4) = [ Zero, - One,    Zero,                        Zero,                         - Vu1, Zero ]
        R(:,5) = [ Zero,  Zero,     One,                        Zero,                           Vu2, Zero ]
        R(:,6) = [ Zero,  Zero,    Zero,                        Zero,                          Zero,  One ]

        invR(:,1) &
          = [ Half * ( H + SQRT( Gmdd33 ) * Mu3 ), &
              One - H, &
              Half * ( H - SQRT( Gmdd33 ) * Mu3 ), &
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
          = [ - Half * ( Q * Mu3 + One / ( SQRT( Gmdd33 ) * Cs  ) ), &
                         Q * Mu3, &
              - Half * ( Q * Mu3 + One / ( SQRT( Gmdd33 ) * Cs  ) ), &
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

    CALL TimersStop_Euler( Timer_Euler_SL_CharDecomp )

  END SUBROUTINE ComputeCharacteristicDecomposition_Euler_NonRelativistic_IDEAL


END MODULE Euler_CharacteristicDecompositionModule_NonRelativistic_IDEAL
