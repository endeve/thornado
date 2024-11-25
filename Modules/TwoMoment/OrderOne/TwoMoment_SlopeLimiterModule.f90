MODULE TwoMoment_SlopeLimiterModule

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    nNodesZ, nDOF, nDOFE, nDOFX, nDimsX
  USE ReferenceElementModule, ONLY: &
    NodeNumbersE
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q
  USE UtilitiesModule, ONLY: &
    MinModB, NodeNumberX
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE RadiationFieldsModule, ONLY: &
    nSpecies, nCR, &
    iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE TwoMoment_CharacteristicDecompositionModule, ONLY: &
    TwoMoment_ComputeCharacteristicDecomposition

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_TwoMoment
  PUBLIC :: FinalizeSlopeLimiter_TwoMoment
  PUBLIC :: ApplySlopeLimiter_TwoMoment

  LOGICAL  :: UseSlopeLimiter
  LOGICAL  :: UseCharacteristicLimiting
  LOGICAL  :: Verbose
  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: I_4x4(1:4,1:4)

CONTAINS


  SUBROUTINE InitializeSlopeLimiter_TwoMoment &
    ( BetaTVD_Option, BetaTVB_Option, SlopeTolerance_Option, &
      UseSlopeLimiter_Option, UseCharacteristicLimiting_Option, &
      Verbose_Option ) 

    REAL(DP), INTENT(in), OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option
    REAL(DP), INTENT(in), OPTIONAL :: &
      SlopeTolerance_Option
    LOGICAL,  INTENT(in), OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      Verbose_Option

    INTEGER :: i

    IF( PRESENT( BetaTVD_Option ) )THEN
      BetaTVD = BetaTVD_Option
    ELSE
      BetaTVD = One
    END IF

    IF( PRESENT( BetaTVB_Option ) )THEN
      BetaTVB = BetaTVB_Option
    ELSE
      BetaTVB = Zero
    END IF

    IF( PRESENT( SlopeTolerance_Option ) )THEN
      SlopeTolerance = SlopeTolerance_Option
    ELSE
      SlopeTolerance = 1.0d-6
    END IF

    IF( PRESENT( UseSlopeLimiter_Option ) )THEN
      UseSlopeLimiter = UseSlopeLimiter_Option
    ELSE
      UseSlopeLimiter = .TRUE.
    END IF

    IF( PRESENT( UseCharacteristicLimiting_Option ) )THEN
      UseCharacteristicLimiting = UseCharacteristicLimiting_Option
    ELSE
      UseCharacteristicLimiting = .FALSE.
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .TRUE.
    END IF

    I_4x4 = Zero
    DO i = 1, 4
      I_4x4(i,i) = One
    END DO

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') '  INFO: InitializeSlopeLimiter_TwoMoment:'
      WRITE(*,'(A)') '  ---------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'      ) '', 'UseSlopeLimiter: ' , &
        UseSlopeLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A27,ES9.3E2)' ) '', 'BetaTVD: ' , &
        BetaTVD
      WRITE(*,'(A4,A27,ES9.3E2)' ) '', 'BetaTVB: ' , &
        BetaTVB
      WRITE(*,'(A4,A27,ES9.3E2)' ) '', 'SlopeTolerance: ' , &
        SlopeTolerance
      WRITE(*,'(A4,A27,L1)'      ) '', 'UseCharacteristicLimiting: ' , &
        UseCharacteristicLimiting
      WRITE(*,*)
    END IF
    
  END SUBROUTINE InitializeSlopeLimiter_TwoMoment


  SUBROUTINE FinalizeSlopeLimiter_TwoMoment

  END SUBROUTINE FinalizeSlopeLimiter_TwoMoment


  SUBROUTINE ApplySlopeLimiter_TwoMoment &
    ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, SuppressBC_Option )

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX, &
         iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nGF)
    REAL(DP), INTENT(inout) :: &
      U (1:nDOF,iZ_B1(1):iZ_E1(1), &
         iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
         1:nCR,1:nSpecies)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option

    LOGICAL  :: LimitedCell(nCR,iZ_B1(2):iZ_E1(2), &
                                iZ_B1(3):iZ_E1(3), &
                                iZ_B1(4):iZ_E1(4))

    INTEGER  :: iZ1, iZ2, iZ3, iZ4, iS, iCR, iGF
    REAL(DP) :: dZ2, dZ3, dZ4
    INTEGER  :: iNodeZ1
    REAL(DP) :: SlopeDifference(nCR)
    REAL(DP) :: GX_K(nGF)
    REAL(DP) :: dU (nCR,nDimsX)
    REAL(DP) :: R_Z2(nCR,nCR), invR_Z2(nCR,nCR)
    REAL(DP) :: R_Z3(nCR,nCR), invR_Z3(nCR,nCR)
    REAL(DP) :: R_Z4(nCR,nCR), invR_Z4(nCR,nCR)
    REAL(DP) :: U_M(nCR,0:2*nDimsX,nDOFX)
    REAL(DP) :: bufferArray(nDOFX)

    IF( .NOT. UseSlopeLimiter ) RETURN

    IF( nDOFX == 1 ) RETURN

    DO iS = 1, nSpecies
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)

      dZ2 = MeshX(1) % Width(iZ2)
      dZ3 = MeshX(2) % Width(iZ3)
      dZ4 = MeshX(3) % Width(iZ4)

      ! --- Cell Average of Geometry (Spatial Metric Only) ---

      DO iGF = iGF_Gm_dd_11, iGF_Gm_dd_33

        GX_K(iGF) = DOT_PRODUCT( WeightsX_q(:), GX(:,iZ2,iZ3,iZ4,iGF) )

      END DO

      DO iZ1 = iZ_B0(1), iZ_E0(1) ! energy loop
      DO iNodeZ1 = 1, nNodesZ(1)

      ! --- Map to Modal Representation ---
        
        DO iCR = 1, nCR

          CALL MapNodalToModal_Fluid &
                 ( U(NodeNumbersE(:,iNodeZ1), &
                     iZ1,iZ2,iZ3,iZ4,iCR,iS), &
                   U_M(iCR,0,:) )    

          ! --- Cell Average of Neighbors in X1 Direction ---
  
          U_M(iCR,1,1) &
            = DOT_PRODUCT &
                ( WeightsX_q(:), &
                  U(NodeNumbersE(:,iNodeZ1),iZ1,iZ2-1,iZ3,iZ4,iCR,iS) ) 
         
          U_M(iCR,2,1) &
            = DOT_PRODUCT &
                ( WeightsX_q(:), &
                  U(NodeNumbersE(:,iNodeZ1),iZ1,iZ2+1,iZ3,iZ4,iCR,iS) )

          IF( nDimsX > 1 )THEN

           ! --- Cell Average of Neighbors in X2 Direction ---

            U_M(iCR,3,1) &
              = DOT_PRODUCT &
                  ( WeightsX_q(:), &
                    U(NodeNumbersE(:,iNodeZ1),iZ1,iZ2,iZ3-1,iZ4,iCR,iS) )
         
            U_M(iCR,4,1) &
              = DOT_PRODUCT &
                  ( WeightsX_q(:), &
                    U(NodeNumbersE(:,iNodeZ1),iZ1,iZ2,iZ3+1,iZ4,iCR,iS) )
           
          END IF

          IF( nDimsX > 2 )THEN

            ! --- Cell Average of Neighbors in X3 Direction ---

            U_M(iCR,5,1) &
              = DOT_PRODUCT &
                  ( WeightsX_q(:), &
                    U(NodeNumbersE(:,iNodeZ1),iZ1,iZ2,iZ3,iZ4-1,iCR,iS) )

            U_M(iCR,6,1) &
              = DOT_PRODUCT &
                  ( WeightsX_q(:), &
                    U(NodeNumbersE(:,iNodeZ1),iZ1,iZ2,iZ3,iZ4+1,iCR,iS) )

          END IF

        END DO ! iCR

      IF( UseCharacteristicLimiting ) THEN
        
        ! --- Compute Eigenvectors
        
        CALL TwoMoment_ComputeCharacteristicDecomposition &
               ( 1, GX_K, U_M(:,0,1), R_Z2, invR_Z2 )

        IF( nDimsX > 1 )THEN

          CALL TwoMoment_ComputeCharacteristicDecomposition &
                 ( 2, GX_K, U_M(:,0,1), R_Z3, invR_Z3 )

        END IF

        IF( nDimsX > 2 )THEN
 
          CALL TwoMoment_ComputeCharacteristicDecomposition &
                 ( 3, GX_K, U_M(:,0,1), R_Z4, invR_Z4 )

        END IF 
     
      ELSE
  
        ! --- Componentwise Limiting ---
  
        R_Z2 = I_4x4; invR_Z2 = I_4x4
        R_Z3 = I_4x4; invR_Z3 = I_4x4
        R_Z4 = I_4x4; invR_Z4 = I_4x4
  
      END IF ! UseCharacteristicLimiting
  
      ! --- Compute Limited Slopes ---
  
      dU(:,1) &
        = MinModB &
            ( MATMUL( invR_Z2, U_M(:,0,2) ), &
              BetaTVD * MATMUL( invR_Z2, (U_M(:,0,1)-U_M(:,1,1)) ), &
              BetaTVD * MATMUL( invR_Z2, (U_M(:,2,1)-U_M(:,0,1)) ), &
              dZ2, BetaTVB )

      IF( nDimsX > 1 )THEN

        dU(:,2) &
          = MinModB &
              ( MATMUL( invR_Z3, U_M(:,0,3) ), &
                BetaTVD * MATMUL( invR_Z3, (U_M(:,0,1)-U_M(:,3,1)) ), &
                BetaTVD * MATMUL( invR_Z3, (U_M(:,4,1)-U_M(:,0,1)) ), &
                dZ3, BetaTVB )

      END IF

      IF( nDimsX > 2 )THEN

        dU(:,3) &
          = MinModB &
              ( MATMUL( invR_Z4, U_M(:,0,4) ), &
                BetaTVD * MATMUL( invR_Z4, (U_M(:,0,1)-U_M(:,5,1)) ), &
                BetaTVD * MATMUL( invR_Z4, (U_M(:,6,1)-U_M(:,0,1)) ), &
                dZ4, BetaTVB )

      END IF

      IF( UseCharacteristicLimiting )THEN

        ! --- Transform Back from Characteristic Variables ---

        dU(:,1) = MATMUL( R_Z2, dU(:,1) )

        IF( nDimsX > 1 )THEN

          dU(:,2) = MATMUL( R_Z3, dU(:,2) )

        END IF

        IF( nDimsX > 2 )THEN

          dU(:,3) = MATMUL( R_Z4, dU(:,3) )

        END IF

      END IF

      ! --- Compare Limited Slopes to Original Slopes ---

      DO iCR = 1, nCR

        SlopeDifference(iCR) &
          = ABS( U_M(iCR,0,2) - dU(iCR,1) )

        IF( nDimsX > 1 )THEN

          SlopeDifference(iCR) &
            = MAX( SlopeDifference(iCR), &
                   ABS( U_M(iCR,0,3) - dU(iCR,2) ) )

        END IF

        IF( nDimsX > 2 )THEN

          SlopeDifference(iCR) &
            = MAX( SlopeDifference(iCR), &
                   ABS( U_M(iCR,0,4) - dU(iCR,3) ) )

        END IF

      END DO

      ! --- Replace Slopes and Discard High-Order Components ---
      ! --- if Limited Slopes Deviate too Much from Original ---

      DO iCR = 1, nCR

        IF( SlopeDifference(iCR) &
              > SlopeTolerance * ABS( U_M(iCR,0,1) ) )THEN

          U_M(iCR,0,2:nDOFX) = Zero

          U_M(:,0,2) = dU(:,1)

          IF( nDimsX > 1 )THEN
            U_M(iCR,0,3) = dU(iCR,2)
          END IF

          IF( nDimsX > 2 )THEN
            U_M(iCR,0,4) = dU(iCR,3)
          END IF

        END IF

      ! --- Map to Nodal Representation ---

        CALL MapModalToNodal_Fluid &
               ( bufferArray, U_M(iCR,0,:) )

        U(NodeNumbersE(:,iNodeZ1), &
            iZ1,iZ2,iZ3,iZ4,iCR,iS) = bufferArray

!!$     LimitedCell(iCR,iX1,iX2,iX3) = .TRUE.

      END DO ! iCR

    END DO ! iNodeZ1
    END DO ! iZ1

    END DO ! iZ2
    END DO ! iZ3
    END DO ! iZ4
    END DO ! iS

  END SUBROUTINE ApplySlopeLimiter_TwoMoment


END MODULE TwoMoment_SlopeLimiterModule
