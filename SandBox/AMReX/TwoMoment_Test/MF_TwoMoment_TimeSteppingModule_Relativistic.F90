MODULE MF_TwoMoment_TimeSteppingModule_Relativistic



  ! --- AMReX Modules ---
  USE amrex_fort_module,      ONLY: &
    AR => amrex_real
  USE amrex_box_module,       ONLY: &
    amrex_box
  USE amrex_geometry_module,  ONLY: &
    amrex_geometry
  USE amrex_multifab_module,  ONLY: &
    amrex_multifab, &
    amrex_multifab_build, amrex_multifab_destroy, &
    amrex_mfiter, &
    amrex_mfiter_build, amrex_mfiter_destroy
  USE amrex_boxarray_module,  ONLY: &
    amrex_boxarray, &
    amrex_boxarray_build, amrex_boxarray_destroy
  USE amrex_distromap_module, ONLY: &
    amrex_distromap, &
    amrex_distromap_build, amrex_distromap_destroy

    ! --- thornado Modules ---
  USE ProgramHeaderModule, ONLY: &
    nDOFZ, nDOFX, nDOFE, &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1, swX, nX
  USE GeometryFieldsModuleE, ONLY: &
    nGE
  USE GeometryFieldsModule, ONLY: &
    nGF
  USE FluidFieldsModule, ONLY: &
    nCF
  USE RadiationFieldsModule, ONLY: &
    nCR
  USE MF_TwoMoment_DiscretizationModule_Streaming_Relativistic, ONLY: &
    ComputeIncrement_TwoMoment_Explicit_MF
  USE MF_TwoMoment_DiscretizationModule_Collisions_Relativistic, ONLY: &
    ComputeIncrement_TwoMoment_Implicit_MF
  USE MF_TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment_MF
  USE MF_TwoMoment_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_TwoMoment_MF
  USE MF_FieldsModule_TwoMoment, ONLY: &
    OffGridFlux_TwoMoment_MF
  ! --- Local Modules ---
  USE InputParsingModule,                      ONLY: &
    nLevels, DEBUG, UseTiling, nSpecies
  USE MF_KindModule, ONLY: &
    DP, &
    Zero, &
    One
  USE MF_TwoMoment_TallyModule, ONLY: &
    IncrementOffGridTally_TwoMoment_MF


  IMPLICIT NONE
  PRIVATE

  TYPE :: StageDataType
    REAL(AR), ALLOCATABLE :: dU_IM(:,:,:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: dU_EX(:,:,:,:,:,:,:)
    REAL(AR), ALLOCATABLE :: dF_EX(:,:,:,:,:)
  END TYPE StageDataType

  INTEGER                          :: nStages
  REAL(AR),            ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(AR),            ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  TYPE(amrex_multifab), DIMENSION(:),   ALLOCATABLE :: MF_U
  TYPE(amrex_multifab), DIMENSION(:),   ALLOCATABLE :: MF_F
  TYPE(amrex_multifab), DIMENSION(:,:), ALLOCATABLE :: MF_DU_Im
  TYPE(amrex_multifab), DIMENSION(:,:), ALLOCATABLE :: MF_DU_Ex
  TYPE(amrex_multifab), DIMENSION(:,:), ALLOCATABLE :: MF_DF_Im
  REAL(AR),            ALLOCATABLE :: U0(:,:,:,:,:,:,:)
  REAL(AR),            ALLOCATABLE :: Ui(:,:,:,:,:,:,:)
  TYPE(StageDataType), ALLOCATABLE :: StageData(:)



  PUBLIC :: Initialize_IMEX_RK
  PUBLIC :: Finalize_IMEX_RK
  PUBLIC :: Initialize_IMEX_RK_MF
  PUBLIC :: Finalize_IMEX_RK_MF
  PUBLIC :: Update_IMEX_RK_MF

CONTAINS


  SUBROUTINE Update_IMEX_RK_MF( t, dt, uGE, MF_uGF, MF_uCF, MF_uCR, GEOM, Verbose_Option )


    REAL(AR),     INTENT(in)    :: t(0:nLevels-1), dt(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uGF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCF(0:nLevels-1)
    TYPE(amrex_multifab), INTENT(inout) :: MF_uCR(0:nLevels-1)
    TYPE(amrex_geometry), INTENT(in)    :: GEOM  (0:nLevels-1)
    REAL(AR),     INTENT(in)    :: uGE(:,:,:)
    LOGICAL,          INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: iLevel, iS, jS

    ! --- For physical boundary conditions ---
    TYPE(amrex_mfiter)                    :: MFI
    REAL(AR), CONTIGUOUS, POINTER :: uCR(:,:,:,:), U(:,:,:,:)
    REAL(AR), CONTIGUOUS, POINTER :: uCF(:,:,:,:), F(:,:,:,:)

    REAL(DP) :: dM_OffGrid_TwoMoment(1:2*nCR,0:nLevels-1)

    INTEGER :: i

    LOGICAL :: Verbose
    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    dM_OffGrid_TwoMoment = Zero
    ! --- Set temporary MultiFabs U and dU to zero --
    DO iLevel = 0, nLevels-1


      CALL MF_U(iLevel) % setval( 0.0_AR )
      CALL MF_F(iLevel) % setval( 0.0_AR )

      DO iS = 1, nStages

        CALL MF_DU_Ex(iLevel,iS) % setval( 0.0_AR )
        CALL MF_DU_Im(iLevel,iS) % setval( 0.0_AR )
        CALL MF_DF_Im(iLevel,iS) % setval( 0.0_AR )

      END DO

    END DO


    DO iS = 1, nStages

      ! --- Copy data from input MultiFab to temporary MultiFab ---

      DO iLevel = 0, nLevels-1


        CALL MF_U(iLevel) &
               % COPY( MF_uCR(iLevel), 1, 1, &
                       MF_uCR(iLevel) % nComp(), swX )


        CALL MF_F(iLevel) &
               % COPY( MF_uCF(iLevel), 1, 1, &
                       MF_uCF(iLevel) % nComp(), swX )
        ! --- Apply boundary conditions to interior domains ---


        CALL MF_U(iLevel) % Fill_Boundary( GEOM(iLevel) )
        CALL MF_F(iLevel) % Fill_Boundary( GEOM(iLevel) )


        ! --- Copy ghost data from physical boundaries ---

        CALL amrex_mfiter_build( MFI, MF_uGF(iLevel), tiling = UseTiling )

        DO WHILE( MFI % next() )

          uCR => MF_uCR(iLevel) % DataPtr( MFI )
          U   => MF_U  (iLevel) % DataPtr( MFI )
          U   =  uCR

          uCF => MF_uCF(iLevel) % DataPtr( MFI )
          F   => MF_F  (iLevel) % DataPtr( MFI )
          F   =  uCF
        END DO

        CALL amrex_mfiter_destroy( MFI )

      END DO



      DO iLevel = 0, nLevels-1


        DO jS = 1, iS - 1


          IF( a_EX(iS,jS) .NE. 0.0_AR )THEN
            CALL MF_U(iLevel) &
                   % LinComb( 1.0_AR,              MF_U(iLevel),    1, &
                              dt(iLevel) * a_EX(iS,jS), MF_DU_Ex(iLevel,jS), 1, &
                              1, MF_U(iLevel) % nComp(), swX )
          END IF

          IF( a_IM(iS,jS) .NE. 0.0_AR )THEN


            CALL MF_U(iLevel) &
                 % LinComb( 1.0_AR,              MF_U(iLevel),    1, &
                            dt(iLevel) * a_IM(iS,jS), MF_DU_Im(iLevel,jS), 1, &
                            1, MF_U(iLevel) % nComp(), swX )
          END IF

          IF( jS == iS - 1 )THEN


            ! --- Apply Limiters ---

            CALL ApplySlopeLimiter_TwoMoment_MF &
                   ( GEOM, MF_uGF, MF_uCF, MF_U, Verbose_Option = Verbose  )

            CALL ApplyPositivityLimiter_TwoMoment_MF &
                   ( GEOM, MF_uGF, MF_uCF, MF_U, Verbose_Option = Verbose  )

          END IF


        END DO ! jS = 1, iS - 1

        IF( ANY( a_IM(:,iS) .NE. 0.0_AR ) .OR. ( w_IM(iS) .NE. 0.0_AR ) )THEN

          IF (Verbose) THEN
            PRINT*, "    IMPLICIT: ", iS
          END IF
          CALL ComputeIncrement_TwoMoment_Implicit_MF &
               ( GEOM, MF_uGF, MF_uCF, MF_U, MF_DU_Im(:,iS), dt(iLevel) * a_IM(iS,iS), Verbose_Option = Verbose )

          CALL MF_U(iLevel) &
                 % LinComb( 1.0_AR,              MF_U(iLevel),    1, &
                            dt(iLevel) * a_IM(iS,iS), MF_DU_Im(iLevel,iS), 1, &
                            1, MF_U(iLevel) % nComp(), swX )


          CALL ApplyPositivityLimiter_TwoMoment_MF &
                  ( GEOM, MF_uGF, MF_uCF, MF_U, Verbose_Option = Verbose  )


        END IF

        IF( ANY( a_EX(:,iS) .NE. 0.0_AR ) .OR. ( w_EX(iS) .NE. 0.0_AR ) )THEN

          ! --- Explicit Solve ---
          IF (Verbose) THEN
            PRINT*, "    EXPLICIT: ", iS
          END IF

          CALL ComputeIncrement_TwoMoment_Explicit_MF &
               ( t, GEOM, MF_uGF, MF_uCF, MF_U, MF_DU_Ex(:,iS), Verbose_Option = Verbose )


          dM_OffGrid_TwoMoment(:,iLevel) &
            = dM_OffGrid_TwoMoment(:,iLevel) &
            + dt(iLevel) * w_EX(iS) * OffGridFlux_TwoMoment_MF(:,iLevel)

        END IF

      END DO

    END DO ! iS = 1, nStages

    ! --- Assembly Step ---

  DO iLevel = 0, nLevels-1
    IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) .OR. &
        ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

      U = uCR

      IF (Verbose) THEN
        PRINT*, "    ASSEMBLY:"
      END IF
        !set Ui to U0 again ask about this
      DO iS = 1, nStages

        IF( w_IM(iS) .NE. 0.0_AR )THEN

          CALL MF_U(iLevel) &
                 % LinComb( 1.0_AR,              MF_U(iLevel),    1, &
                            dt(iLevel) * w_IM(iS), MF_DU_Im(iLevel,iS), 1, &
                            1, MF_U(iLevel) % nComp(), swX )

        END IF

        IF( w_EX(iS) .NE. 0.0_AR )THEN

          CALL MF_U(iLevel) &
                 % LinComb( 1.0_AR,              MF_U(iLevel),    1, &
                            dt(iLevel) * w_EX(iS), MF_DU_Ex(iLevel,iS), 1, &
                            1, MF_U(iLevel) % nComp(), swX )

        END IF

      END DO

        CALL ApplySlopeLimiter_TwoMoment_MF &
                   ( GEOM, MF_uGF, MF_uCF, MF_U, Verbose_Option = Verbose  )

        CALL ApplyPositivityLimiter_TwoMoment_MF &
                   ( GEOM, MF_uGF, MF_uCF, MF_U, Verbose_Option = Verbose  )

    END IF

  END DO

  CALL IncrementOffGridTally_TwoMoment_MF( dM_OffGrid_TwoMoment )


  uCR = U

  END SUBROUTINE Update_IMEX_RK_MF

  SUBROUTINE Initialize_IMEX_RK_MF &
    ( Scheme, BA, DM, Verbose_Option )

    CHARACTER(LEN=*), INTENT(in) :: Scheme
    TYPE(amrex_boxarray),  INTENT(in)           :: BA(0:nLevels-1)
    TYPE(amrex_distromap), INTENT(in)           :: DM(0:nLevels-1)
    LOGICAL,               INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER         :: iS, iLevel
    TYPE(amrex_box) :: BX

    CALL Initialize_IMEX_RK( Scheme, Verbose_Option )

    ALLOCATE( MF_U(0:nLevels-1) )
    ALLOCATE( MF_F(0:nLevels-1) )
    ALLOCATE( MF_DU_Ex(0:nLevels-1,1:nStages) )
    ALLOCATE( MF_DU_Im(0:nLevels-1,1:nStages) )
    ALLOCATE( MF_DF_Im(0:nLevels-1,1:nStages) )

    BX = amrex_box( [ 0, 0, 0 ], [ nX(1)-1, nX(2)-1, nX(3)-1 ] )

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_build &
        ( MF_U(iLevel), BA(iLevel), DM(iLevel), &
          nDOFZ * nCR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX )

      CALL amrex_multifab_build &
        ( MF_F(iLevel), BA(iLevel), DM(iLevel), nDOFX * nCF, swX )

      DO iS = 1, nStages

        CALL amrex_multifab_build &
               ( MF_DU_Ex(iLevel,iS), BA(iLevel), DM(iLevel), &
                 nDOFZ * nCR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX )

        CALL amrex_multifab_build &
               ( MF_DU_Im(iLevel,iS), BA(iLevel), DM(iLevel), &
                 nDOFZ * nCR * ( iZ_E0( 1 ) - iZ_B0( 1 ) + 1 ) * nSpecies, swX )

        CALL amrex_multifab_build &
               ( MF_DF_Im(iLevel,iS), BA(iLevel), DM(iLevel), nDOFX * nCF, swX )

      END DO

    END DO


  END SUBROUTINE Initialize_IMEX_RK_MF


  SUBROUTINE Finalize_IMEX_RK_MF

    INTEGER :: iLevel, iS

    CALL Finalize_IMEX_RK

    DO iLevel = 0, nLevels-1

      CALL amrex_multifab_destroy( MF_U(iLevel) )

      DO iS = 1, nStages

        CALL amrex_multifab_destroy( MF_DU_Ex(iLevel,iS) )
        CALL amrex_multifab_destroy( MF_DU_Im(iLevel,iS) )
        CALL amrex_multifab_destroy( MF_DF_Im(iLevel,iS) )

      END DO

    END DO

    DEALLOCATE( MF_U )
    DEALLOCATE( MF_DU_Ex)
    DEALLOCATE( MF_DU_Im)
    DEALLOCATE( MF_DF_Im)

  END SUBROUTINE Finalize_IMEX_RK_MF

  SUBROUTINE Initialize_IMEX_RK( Scheme, Verbose_Option )

    CHARACTER(LEN=*), INTENT(in) :: Scheme
    LOGICAL,               INTENT(in), OPTIONAL :: Verbose_Option

    INTEGER :: i
    LOGICAL :: Verbose

    Verbose = .FALSE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF (Verbose) THEN
      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'IMEX-RK Scheme: ', TRIM( Scheme )
    END IF

    SELECT CASE ( TRIM( Scheme ) )

      CASE ( 'BackwardEuler' )

        nStages = 1

        CALL AllocateButcherTables

        a_IM(1,1) = 1.0_AR
        w_IM(1)   = 1.0_AR

      CASE ( 'SSPRK1' )

        nStages = 1

        CALL AllocateButcherTables

        a_EX(1,1) = 0.0_AR
        w_EX(1)   = 1.0_AR

      CASE ( 'SSPRK2' )

        nStages = 2

        CALL AllocateButcherTables

        a_EX(1,1:2) = [ 0.0_AR, 0.0_AR ]
        a_EX(2,1:2) = [ 1.0_AR, 0.0_AR ]
        w_EX(1:2)   = [ 0.5_AR, 0.5_AR ]

      CASE ( 'SSPRK3' )

        nStages = 3

        CALL AllocateButcherTables

        a_EX(1,1:3) = [ 0.00_AR, 0.00_AR, 0.00_AR ]
        a_EX(2,1:3) = [ 1.00_AR, 0.00_AR, 0.00_AR ]
        a_EX(3,1:3) = [ 0.25_AR, 0.25_AR, 0.00_AR ]
        w_EX(1:3)   = [ 1.00_AR, 1.00_AR, 4.00_AR ] / 6.0_AR

      CASE ( 'IMEX_ARS_111' )

        nStages = 2

        CALL AllocateButcherTables

        ! --- Coefficients from Ascher et al. (1997) ---

        a_EX(2,1) = 1.0_AR
        w_EX(1)   = a_EX(2,1)

        a_IM(2,2) = 1.0_AR
        w_IM(2)   = a_IM(2,2)

      CASE ( 'IMEX_PDARS' )

        nStages = 3

        CALL AllocateButcherTables

        a_EX(2,1) = 1.0_AR
        a_EX(3,1) = 0.5_AR
        a_EX(3,2) = 0.5_AR

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        a_IM(2,2) = 1.0_AR
        a_IM(3,2) = 0.5_AR
        a_IM(3,3) = 0.5_AR

        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A6,A,A)') &
          '', 'Unknown Time Stepping Scheme: ', TRIM( Scheme )
        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'Available Options:'
        WRITE(*,*)
        WRITE(*,'(A6,A)') '', 'BackwardEuler'
        WRITE(*,'(A6,A)') '', 'SSPRK1'
        WRITE(*,'(A6,A)') '', 'SSPRK2'
        WRITE(*,'(A6,A)') '', 'SSPRK3'
        WRITE(*,'(A6,A)') '', 'IMEX_ARS_111'
        WRITE(*,'(A6,A)') '', 'IMEX_PDARS'
        WRITE(*,*)
        STOP

    END SELECT

    DO i = 1, nStages
      c_IM(i) = SUM( a_IM(i,1:i) )
      c_EX(i) = SUM( a_EX(i,1:i-1) )
    END DO

    IF (Verbose) THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Implicit Butcher Table:'
      WRITE(*,'(A6,A)') '', '-----------------------'
      DO i = 1, nStages
        WRITE(*,'(A6,5ES14.4E3)') '', c_IM(i), a_IM(i,1:nStages)
      END DO
      WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_IM(1:nStages)

      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'Explicit Butcher Table:'
      WRITE(*,'(A6,A)') '', '-----------------------'
      DO i = 1, nStages
        WRITE(*,'(A6,5ES14.4E3)') '', c_EX(i), a_EX(i,1:nStages)
      END DO
      WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_EX(1:nStages)
    END IF

    CALL AllocateArray7D( U0 )
    CALL AllocateArray7D( Ui )

    ALLOCATE( StageData(nStages) )

    DO i = 1, nStages


      CALL AllocateArray7D( StageData(i) % dU_IM )
      CALL AllocateArray7D( StageData(i) % dU_EX )
    END DO




  END SUBROUTINE Initialize_IMEX_RK



  SUBROUTINE Finalize_IMEX_RK

    INTEGER :: i

    DEALLOCATE( c_IM, w_IM, a_IM )
    DEALLOCATE( c_EX, w_EX, a_EX )

    CALL DeallocateArray7D( U0 )
    CALL DeallocateArray7D( Ui )

    DO i = 1, nStages

      CALL DeallocateArray7D( StageData(i) % dU_IM )
      CALL DeallocateArray7D( StageData(i) % dU_EX )

    END DO

    DEALLOCATE( StageData )

  END SUBROUTINE Finalize_IMEX_RK


  SUBROUTINE AllocateButcherTables

    ! --- Implicit Coefficients ---

    ALLOCATE( c_IM(nStages) )
    ALLOCATE( w_IM(nStages) )
    ALLOCATE( a_IM(nStages,nStages) )

    c_IM = 0.0_AR
    w_IM = 0.0_AR
    a_IM = 0.0_AR

    ! --- Explicit Coefficients ---

    ALLOCATE( c_EX(nStages) )
    ALLOCATE( w_EX(nStages) )
    ALLOCATE( a_EX(nStages,nStages) )

    c_EX = 0.0_AR
    w_EX = 0.0_AR
    a_EX = 0.0_AR

  END SUBROUTINE AllocateButcherTables


  SUBROUTINE AllocateArray7D( Array7D )

    REAL(AR), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    ALLOCATE &
      ( Array7D(nDOFZ, &
                iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
                nCR,nSpecies) )

    Array7D = 0.0_AR

  END SUBROUTINE AllocateArray7D

  SUBROUTINE AllocateArray5D( Array5D )

    REAL(AR), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    ALLOCATE &
      ( Array5D(nDOFZ, &
                iZ_B1(2):iZ_E1(2), &
                iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4), &
                nCF) )

    Array5D = 0.0_AR

  END SUBROUTINE AllocateArray5D

  SUBROUTINE DeallocateArray7D( Array7D )

    REAL(AR), ALLOCATABLE, INTENT(inout) :: Array7D(:,:,:,:,:,:,:)

    DEALLOCATE( Array7D )

  END SUBROUTINE DeallocateArray7D

  SUBROUTINE DeallocateArray5D( Array5D )

    REAL(AR), ALLOCATABLE, INTENT(inout) :: Array5D(:,:,:,:,:)

    DEALLOCATE( Array5D )

  END SUBROUTINE DeallocateArray5D



END MODULE MF_TwoMoment_TimeSteppingModule_Relativistic
