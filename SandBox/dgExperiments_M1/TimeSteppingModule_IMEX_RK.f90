MODULE TimeSteppingModule_IMEX_RK

  USE KindModule, ONLY: &
    DP, Zero
  USE ProgramHeaderModule, ONLY: &
    iZ_B0, iZ_B1, iZ_E0, iZ_E1, &
    nDOF
  USE RadiationFieldsModule, ONLY: &
    nCR, nSpecies
  USE TwoMoment_SlopeLimiterModule, ONLY: &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    ApplyPositivityLimiter_TwoMoment

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  INTEGER               :: nStages
  REAL(DP)              :: alpha_IM, alpha_EX
  REAL(DP), ALLOCATABLE :: c_IM(:), w_IM(:), a_IM(:,:)
  REAL(DP), ALLOCATABLE :: c_EX(:), w_EX(:), a_EX(:,:)
  REAL(DP), ALLOCATABLE :: U_IMEX(:), dU_IM(:,:), dU_EX(:,:)
  REAL(DP), ALLOCATABLE :: U0(:,:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: dU(:,:,:,:,:,:,:)

  PUBLIC :: Initialize_IMEX_RK
  PUBLIC :: Finalize_IMEX_RK
  PUBLIC :: Update_IMEX_RK

  INTERFACE
    SUBROUTINE IncrementExplicit &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )
      USE KindModule, ONLY: DP
      USE ProgramHeaderModule, ONLY: &
        nDOF, nDOFE, nDOFX
      USE GeometryFieldsModuleE, ONLY: &
        nGE
      USE GeometryFieldsModule, ONLY: &
        nGF
      USE RadiationFieldsModule, ONLY: &
        nCR, nSpecies

    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      GE(1:nDOFE,iZ_B1(1):iZ_E1(1),1:nGE)
    REAL(DP), INTENT(in)    :: &
      GX(1:nDOFX,iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nGF)
    REAL(DP), INTENT(inout) :: &
      U (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out)   :: &
      dU(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2), &
                 iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    END SUBROUTINE IncrementExplicit
  END INTERFACE

  INTERFACE
    SUBROUTINE IncrementImplicit &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, GE, GX, U, dU )
      USE KindModule, ONLY: DP
      INTEGER, INTENT(in)     :: &
        iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
      REAL(DP), INTENT(in)    :: &
        dt
      REAL(DP), INTENT(in)    :: &
        GE(1:,iZ_B1(1):,1:)
      REAL(DP), INTENT(in)    :: &
        GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
      REAL(DP), INTENT(inout) :: &
        U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
      REAL(DP), INTENT(inout) :: &
        dU(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    END SUBROUTINE IncrementImplicit
  END INTERFACE

  INTERFACE
    SUBROUTINE IncrementCorrection &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, alpha_EX, alpha_IM, GE, GX, U, dU )
      USE KindModule, ONLY: DP
      INTEGER, INTENT(in)     :: &
        iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
      REAL(DP), INTENT(in)    :: &
        dt, &
        alpha_EX, &
        alpha_IM
      REAL(DP), INTENT(in)    :: &
        GE(1:,iZ_B1(1):,1:)
      REAL(DP), INTENT(in)    :: &
        GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
      REAL(DP), INTENT(inout) :: &
        U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
      REAL(DP), INTENT(inout) :: &
        dU(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    END SUBROUTINE IncrementCorrection
  END INTERFACE

CONTAINS


  SUBROUTINE Initialize_IMEX_RK( Scheme )

    CHARACTER(LEN=*), INTENT(in) :: Scheme

    INTEGER :: i, nDOF_IMEX

    WRITE(*,*)
    WRITE(*,'(A6,A,A)') '', 'IMEX-RK Scheme: ', TRIM( Scheme )

    SELECT CASE ( TRIM( Scheme ) )

      CASE ( 'SSPRK1' )

        nStages = 1
        CALL AllocateButcherTables( nStages )

        a_EX(1,1) = 0.0_DP
        w_EX(1)   = 1.0_DP

      CASE ( 'SSPRK2' )

        nStages = 2
        CALL AllocateButcherTables( nStages )

        a_EX(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_EX(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_EX(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 'SSPRK3' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_EX(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_EX(1:3)   = [ 1.0_DP / 6.0_DP, &
                        1.0_DP / 6.0_DP, &
                        2.0_DP / 3.0_DP ]

      CASE ( 'IMEX_ARS_111' )

        nStages = 2
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Ascher et al. (1997) ---

        a_EX(2,1) = 1.0_DP
        w_EX(1)   = a_EX(2,1)

        a_IM(2,2) = 1.0_DP
        w_IM(2)   = a_IM(2,2)

      CASE ( 'IMEX_P_A2' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Hu et al. (2017) ---
        ! --- arXiv:1708.06279v1, Section 2.6.1 ----

        a_EX(2,1) = 0.7369502715_DP
        a_EX(3,1) = 0.3215281691_DP
        a_EX(3,2) = 0.6784718309_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        a_IM(1,1) = 0.6286351712_DP
        a_IM(2,1) = 0.2431004655_DP
        a_IM(2,2) = 0.1959392570_DP
        a_IM(3,1) = 0.4803651051_DP
        a_IM(3,2) = 0.0746432814_DP
        a_IM(3,3) = 0.4449916135_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

        alpha_IM = 0.2797373792_DP

      CASE ( 'IMEX_P_A2_RC' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Ran Chu (2018) ---

        a_EX(2,1) = 0.909090909090909_DP
        a_EX(3,1) = 0.450000000000000_DP
        a_EX(3,2) = 0.550000000000000_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        a_IM(1,1) = 0.521932391842510_DP
        a_IM(2,1) = 0.479820781424967_DP
        a_IM(2,2) = 0.002234534340252_DP
        a_IM(3,1) = 0.499900000000000_DP
        a_IM(3,2) = 0.001100000000000_DP
        a_IM(3,3) = 0.499000000000000_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

        alpha_IM = 0.260444263529413_DP

      CASE ( 'IMEX_P_A2_RC2' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Ran Chu (2018) ---

        a_EX(2,1) = 1.063612014417258_DP
        a_EX(3,1) = 0.529903768270289_DP
        a_EX(3,2) = 0.470096231729711_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        alpha_EX = 0.377549414083044_DP

        a_IM(1,1) = 0.231157005614910_DP
        a_IM(2,1) = 0.226956591610883_DP
        a_IM(2,2) = 0.576089691643813_DP
        a_IM(3,1) = 0.138355033665218_DP
        a_IM(3,2) = 0.250022702377217_DP
        a_IM(3,3) = 0.611622263957565_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

        alpha_IM = 0.344383801124686_DP

      CASE ( 'IMEX_P_A2_RC3' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Ran Chu (2018) ---

        a_EX(2,1) = 0.702049509161748_DP
        a_EX(3,1) = 0.287799516309037_DP
        a_EX(3,2) = 0.712200483690963_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        alpha_EX = 0.449518316898217_DP

        a_IM(1,1) = 0.000000122692639_DP
        a_IM(2,1) = 0.000000022134039_DP
        a_IM(2,2) = 0.702049437447734_DP
        a_IM(3,1) = 0.000000017307164_DP
        a_IM(3,2) = 0.169429711804114_DP
        a_IM(3,3) = 0.830570270888722_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

        alpha_IM = 0.449518308497898_DP

      CASE ( 'IMEX_P_A2_RC4' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Ran Chu (2018) ---

        a_EX(2,1) = 0.975992486704880_DP
        a_EX(3,1) = 0.487700974330154_DP
        a_EX(3,2) = 0.512299025669846_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        alpha_EX = 0.445978165340310_DP

        a_IM(1,1) = 0.390570636268524_DP
        a_IM(2,1) = 0.104174290767753_DP
        a_IM(2,2) = 0.500000819154403_DP
        a_IM(3,1) = 0.052692998443186_DP
        a_IM(3,2) = 0.055350847885380_DP
        a_IM(3,3) = 0.891956153671434_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

        alpha_IM = 0.445978096205720_DP

      CASE ( 'IMEX_P_A2_RC5' )

        nStages = 3
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Ran Chu (2018) ---

        a_EX(2,1) = 0.726393726720026_DP
        a_EX(3,1) = 0.311668064291096_DP
        a_EX(3,2) = 0.688331935708904_DP

        w_EX(1)   = a_EX(3,1)
        w_EX(2)   = a_EX(3,2)

        alpha_EX = 0.295233192284502_DP

        a_IM(1,1) = 0.292971931709409_DP
        a_IM(2,1) = 0.286422870538316_DP
        a_IM(2,2) = 0.307316841293593_DP
        a_IM(3,1) = 0.170547969724994_DP
        a_IM(3,2) = 0.125065984709675_DP
        a_IM(3,3) = 0.704386045565331_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

        alpha_IM = 0.328608455426277_DP

      CASE ( 'IMEX_P_ARS2' )

        nStages = 4
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Hu et al. (2017) ---
        ! --- arXiv:1708.06279v1, Section 2.6.2 ----

        a_EX(3,1) = 1.0_DP
        a_EX(4,1) = 0.5_DP
        a_EX(4,3) = 0.5_DP

        w_EX(1)   = a_EX(4,1)
        w_EX(3)   = a_EX(4,3)

        a_IM(2,2) = 1.6_DP
        a_IM(3,2) = 0.3_DP
        a_IM(3,3) = 0.7_DP
        a_IM(4,2) = 0.5_DP
        a_IM(4,3) = 0.3_DP
        a_IM(4,4) = 0.2_DP

        w_IM(2)   = a_IM(4,2)
        w_IM(3)   = a_IM(4,3)
        w_IM(4)   = a_IM(4,4)

        alpha_IM = 0.8_DP

      CASE ( 'IMEX_P_ARS2_RC' )

        nStages = 4
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Ran Chu (2018) ---

        a_EX(3,1) = 1.0_DP
        a_EX(4,1) = 0.5_DP
        a_EX(4,3) = 0.5_DP

        w_EX(1)   = a_EX(4,1)
        w_EX(3)   = a_EX(4,3)

        a_IM(2,2) = 4.062494753449722_DP
        a_IM(3,2) = 0.001699788261156_DP
        a_IM(3,3) = 0.998300211738844_DP
        a_IM(4,2) = 0.500000000000000_DP
        a_IM(4,3) = 0.499110102279259_DP
        a_IM(4,4) = 0.000889897720741_DP

        w_IM(2)   = a_IM(4,2)
        w_IM(3)   = a_IM(4,3)
        w_IM(4)   = a_IM(4,4)

        alpha_IM = 2.031247376724861_DP

      CASE ( 'IMEX_P_ARS2_RC2' )

        nStages = 4
        CALL AllocateButcherTables( nStages )

        ! --- Coefficients from Ran Chu (2018) ---

        a_EX(2,1) = 0.025819295974486_DP
        a_EX(3,1) = 1.556487769284859_DP
        a_EX(3,2) = 0.040322522993018_DP
        a_EX(4,1) = 0.579878574268125_DP
        a_EX(4,2) = 0.108755694346525_DP
        a_EX(4,3) = 0.311365731385350_DP

        w_EX(1)   = a_EX(4,1)
        w_EX(2)   = a_EX(4,2)
        w_EX(3)   = a_EX(4,3)

        a_IM(2,2) = 0.959327206951534_DP
        a_IM(3,2) = 0.322105700620029_DP
        a_IM(3,3) = 0.948643336989347_DP
        a_IM(4,2) = 0.579878574268125_DP
        a_IM(4,3) = 0.108755694346525_DP
        a_IM(4,4) = 0.311365731385350_DP

        w_IM(2)   = a_IM(4,2)
        w_IM(3)   = a_IM(4,3)
        w_IM(4)   = a_IM(4,4)

        alpha_IM = 0.505860218334414_DP

      CASE ( 'IMEX_SSP2332' )

        ! --- Coefficients from Pareschi & Russo (2005) ---
        ! --- J. Sci. Comput. 25, 129-154 -----------------

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(2,1) = 0.5_DP
        a_EX(3,1) = 0.5_DP
        a_EX(3,2) = 0.5_DP

        w_EX(1)   = 1.0_DP / 3.0_DP
        w_EX(2)   = 1.0_DP / 3.0_DP
        w_EX(3)   = 1.0_DP / 3.0_DP

        a_IM(1,1) = 0.25_DP
        a_IM(2,2) = 0.25_DP
        a_IM(3,1) = 1.0_DP / 3.0_DP
        a_IM(3,2) = 1.0_DP / 3.0_DP
        a_IM(3,3) = 1.0_DP / 3.0_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

      CASE ( 'IMEX_SSP2322' )

        ! --- Coefficients from Pareschi & Russo (2005) ---
        ! --- J. Sci. Comput. 25, 129-154 -----------------

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(3,2) = 1.0_DP

        w_EX(2)   =  0.5_DP
        w_EX(3)   =  0.5_DP

        a_IM(1,1) =  0.5_DP
        a_IM(2,1) = -0.5_DP
        a_IM(2,2) =  0.5_DP
        a_IM(3,2) =  0.5_DP
        a_IM(3,3) =  0.5_DP

        w_IM(1)   = a_IM(3,1)
        w_IM(2)   = a_IM(3,2)
        w_IM(3)   = a_IM(3,3)

      CASE ( 'IMEX_RKCB2' )

        ! --- Coefficients from Cavaglieri & Bewley (2015) ---
        ! --- JCP, 286, 172-193 ------------------------------

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(2,1) = 0.4_DP
        a_EX(3,2) = 1.0_DP

        w_EX(2)   = 5.0_DP / 6.0_DP
        w_EX(3)   = 1.0_DP / 6.0_DP

        a_IM(2,2) = 0.4_DP
        a_IM(3,2) = 5.0_DP / 6.0_DP
        a_IM(3,3) = 1.0_DP / 6.0_DP

        w_IM(2)   = 5.0_DP / 6.0_DP
        w_IM(3)   = 1.0_DP / 6.0_DP

      CASE ( 'IMEX_SIRK2' )

        ! --- Scheme from Chertock et al. (2015) ---
        ! --- SIAM J. Numer. Anal. 53, 2008-2029 ---

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(2,1) = 1.0_DP
        a_EX(3,1) = 1.0_DP
        a_EX(3,2) = 1.0_DP

        w_EX(1)   = 0.5_DP
        w_EX(2)   = 0.5_DP

        alpha_EX = 1.0_DP

        a_IM(2,2) = 1.0_DP
        a_IM(3,2) = 1.0_DP
        a_IM(3,3) = 1.0_DP

        w_IM(2)   = 0.5_DP
        w_IM(3)   = 0.5_DP

        alpha_IM = 1.0_DP

      CASE ( 'IMEX_PDARS_3' )

        ! --- Scheme from Ran Chu (2018) ---

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(2,1) = 1.0_DP
        a_EX(3,1) = 0.5_DP
        a_EX(3,2) = 0.5_DP

        w_EX(1)   = 0.5_DP
        w_EX(2)   = 0.5_DP

        alpha_EX = 0.0_DP

        a_IM(2,2) = 1.0_DP
        a_IM(3,2) = 0.5_DP
        a_IM(3,3) = 0.5_DP

        w_IM(2)   = 0.5_DP
        w_IM(3)   = 0.5_DP

        alpha_IM = 0.0_DP

      CASE ( 'IMEX_PDARS_4' )

        ! --- Scheme from Ran Chu (2018) ---

        nStages = 4
        CALL AllocateButcherTables( nStages )

        a_EX(2,1) = 1.0_DP
        a_EX(3,1) = 0.25_DP
        a_EX(3,2) = 0.25_DP
        a_EX(4,1) = 1.0_DP / 6.0_DP
        a_EX(4,2) = 1.0_DP / 6.0_DP
        a_EX(4,3) = 2.0_DP / 3.0_DP

        w_EX(1)   = 1.0_DP / 6.0_DP
        w_EX(2)   = 1.0_DP / 6.0_DP
        w_EX(3)   = 2.0_DP / 3.0_DP

        alpha_EX = 0.0_DP

        a_IM(2,2) = 1.0_DP
        a_IM(3,2) = 0.25_DP
        a_IM(3,3) = 0.25_DP
        a_IM(4,2) = 1.0_DP / 6.0_DP - (1.d-6)/4.0_DP
        a_IM(4,3) = 1.0_DP / 6.0_DP - (1.d-3)/4.0_DP
        a_IM(4,4) = 2.0_DP / 3.0_DP + (1.d-3)/4.0_DP + (1.d-6)/4.0_DP

        w_IM(2)   = a_IM(4,2)
        w_IM(3)   = a_IM(4,3)
        w_IM(4)   = a_IM(4,4)

        alpha_IM = 0.0_DP

      CASE ( 'IMEX_PC2' )

        ! --- Scheme from McClarren et al. (2008) ---
        ! --- JCP, 227, 7561-7586 -------------------

        nStages = 3
        CALL AllocateButcherTables( nStages )

        a_EX(2,1) = 0.5_DP
        a_EX(3,2) = 1.0_DP

        w_EX(2)   = a_EX(3,2)

        a_IM(2,2) = 0.5_DP
        a_IM(3,3) = 1.0_DP

        w_IM(3)   = a_IM(3,3)

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A6,A,A)') &
          '', 'Unknown Time Stepping Scheme: ', TRIM( Scheme )
        WRITE(*,*)
        WRITE(*,'(A6,A)') &
          '', 'Available Options:'
        WRITE(*,*)
        WRITE(*,'(A6,A)') '', 'SSPRK1'
        WRITE(*,'(A6,A)') '', 'SSPRK2'
        WRITE(*,'(A6,A)') '', 'SSPRK3'
        WRITE(*,'(A6,A)') '', 'IMEX_ARS_111'
        WRITE(*,'(A6,A)') '', 'IMEX_P_A2'
        WRITE(*,'(A6,A)') '', 'IMEX_P_A2_RC'
        WRITE(*,'(A6,A)') '', 'IMEX_P_A2_RC2'
        WRITE(*,'(A6,A)') '', 'IMEX_P_A2_RC3'
        WRITE(*,'(A6,A)') '', 'IMEX_P_A2_RC4'
        WRITE(*,'(A6,A)') '', 'IMEX_P_A2_RC5'
        WRITE(*,'(A6,A)') '', 'IMEX_P_ARS2'
        WRITE(*,'(A6,A)') '', 'IMEX_P_ARS2_RC'
        WRITE(*,'(A6,A)') '', 'IMEX_P_ARS2_RC2'
        WRITE(*,'(A6,A)') '', 'IMEX_SSP2332'
        WRITE(*,'(A6,A)') '', 'IMEX_SSP2322'
        WRITE(*,'(A6,A)') '', 'IMEX_RKCB2'
        WRITE(*,'(A6,A)') '', 'IMEX_SIRK2'
        WRITE(*,'(A6,A)') '', 'IMEX_PDARS_3'
        WRITE(*,'(A6,A)') '', 'IMEX_PDARS_4'
        WRITE(*,'(A6,A)') '', 'IMEX_PC2'

        WRITE(*,*)
        STOP

    END SELECT

    DO i = 1, nStages
      c_IM(i) = SUM( a_IM(i,1:i) )
      c_EX(i) = SUM( a_EX(i,1:i-1) )
    END DO

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Implicit Butcher Table:'
    WRITE(*,'(A6,A)') '', '-----------------------'
    DO i = 1, nStages
      WRITE(*,'(A6,5ES14.4E3)') '', c_IM(i), a_IM(i,1:nStages)
    END DO
    WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_IM(1:nStages)
    WRITE(*,*)
    WRITE(*,'(A6,A8,ES11.4E3)') '', 'alpha = ', alpha_IM

    WRITE(*,*)
    WRITE(*,'(A6,A)') '', 'Explicit Butcher Table:'
    WRITE(*,'(A6,A)') '', '-----------------------'
    DO i = 1, nStages
      WRITE(*,'(A6,5ES14.4E3)') '', c_EX(i), a_EX(i,1:nStages)
    END DO
    WRITE(*,'(A6,A14,4ES14.4E3)') '', '', w_EX(1:nStages)
    WRITE(*,*)
    WRITE(*,'(A6,A8,ES11.4E3)') '', 'alpha = ', alpha_EX

    ALLOCATE &
      ( U0(1:nDOF, &
           iZ_B0(1):iZ_E0(1), &
           iZ_B0(2):iZ_E0(2), &
           iZ_B0(3):iZ_E0(3), &
           iZ_B0(4):iZ_E0(4), &
           1:nCR, 1:nSpecies) )

    ALLOCATE &
      ( dU(1:nDOF, &
           iZ_B1(1):iZ_E1(1), &
           iZ_B1(2):iZ_E1(2), &
           iZ_B1(3):iZ_E1(3), &
           iZ_B1(4):iZ_E1(4), &
           1:nCR, 1:nSpecies) )

    nDOF_IMEX = nDOF * PRODUCT( iZ_E0 ) * nCR * nSpecies

    ALLOCATE( U_IMEX(1:nDOF_IMEX) )
    ALLOCATE( dU_IM (1:nDOF_IMEX,1:nStages) )
    ALLOCATE( dU_EX (1:nDOF_IMEX,1:nStages) )

  END SUBROUTINE Initialize_IMEX_RK


  SUBROUTINE Finalize_IMEX_RK

    DEALLOCATE( c_IM, w_IM, a_IM )
    DEALLOCATE( c_EX, w_EX, a_EX )
    DEALLOCATE( U0, dU )
    DEALLOCATE( U_IMEX, dU_IM, dU_EX )

  END SUBROUTINE Finalize_IMEX_RK


  SUBROUTINE AllocateButcherTables( nStages )

    INTEGER, INTENT(in) :: nStages

    ! --- Implicit Coefficients ---

    alpha_IM = Zero

    ALLOCATE( c_IM(nStages) )
    ALLOCATE( w_IM(nStages) )
    ALLOCATE( a_IM(nStages,nStages) )

    c_IM = Zero
    w_IM = Zero
    a_IM = Zero

    ! --- Explicit Coefficients ---

    alpha_EX = Zero

    ALLOCATE( c_EX(nStages) )
    ALLOCATE( w_EX(nStages) )
    ALLOCATE( a_EX(nStages,nStages) )

    c_EX = Zero
    w_EX = Zero
    a_EX = Zero

  END SUBROUTINE AllocateButcherTables


  SUBROUTINE Update_IMEX_RK &
    ( dt, GE, GX, U, ComputeIncrement_Explicit, &
      ComputeIncrement_Implicit, ComputeCorrection_Implicit )

    REAL(DP), INTENT(in)           :: &
      dt
    REAL(DP), INTENT(in)           :: &
      GE(1:,iZ_B1(1):,1:)
    REAL(DP), INTENT(in)           :: &
      GX(1:,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:)
    REAL(DP), INTENT(inout)        :: &
      U (1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)
    PROCEDURE(IncrementExplicit)   :: &
      ComputeIncrement_Explicit
    PROCEDURE(IncrementImplicit)   :: &
      ComputeIncrement_Implicit
    PROCEDURE(IncrementCorrection) :: &
      ComputeCorrection_Implicit

    INTEGER :: iS, jS

    CALL InitializeStep_IMEX_RK( U )

    DO iS = 1, nStages

      CALL MapToStage( iZ_B0, U0, U_IMEX )

      DO jS = 1, iS - 1

        IF( a_IM(iS,jS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * a_IM(iS,jS) * dU_IM(:,jS)

        END IF

        IF( a_EX(iS,jS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * a_EX(iS,jS) * dU_EX(:,jS)

        END IF

        IF( jS == iS - 1 )THEN

          ! --- Apply Positivity and Slope Limiter ---

          CALL MapFromStage( iZ_B1, U, U_IMEX )

          CALL ApplySlopeLimiter_TwoMoment &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

          CALL ApplyPositivityLimiter_TwoMoment &
                 ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

          CALL MapToStage( iZ_B1, U, U_IMEX )
        
        END IF

      END DO ! jS = 1, iS - 1

      IF( ANY( a_IM(:,iS) .NE. Zero ) .OR. ( w_IM(iS) .NE. Zero ) )THEN

        ! --- Implicit Solve ---

        CALL MapFromStage( iZ_B1, U, U_IMEX )

        CALL ComputeIncrement_Implicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt * a_IM(iS,iS), GE, GX, U, dU )

        CALL MapToStage( iZ_B1, dU, dU_IM(:,iS) )

        U_IMEX(:) = U_IMEX(:) + dt * a_IM(iS,iS) * dU_IM(:,iS)

        ! --- Apply Positivity and Slope Limiter ---

        CALL MapFromStage( iZ_B1, U, U_IMEX )

        CALL ApplySlopeLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

        CALL ApplyPositivityLimiter_TwoMoment &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

        CALL MapToStage( iZ_B1, U, U_IMEX )

      END IF

      IF( ANY( a_EX(:,iS) .NE. Zero ) .OR. ( w_EX(iS) .NE. Zero ) )THEN

        ! --- Explicit Solve ---

        CALL MapFromStage( iZ_B1, U, U_IMEX )
        
        CALL ComputeIncrement_Explicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

        CALL MapToStage( iZ_B1, dU, dU_EX(:,iS) )

      END IF

    END DO ! iS

    ! --- Weighted Sum ---

    IF( ANY( a_IM(nStages,:) .NE. w_IM(:) ) &
          .OR. ANY( a_EX(nStages,:) .NE. w_EX(:) ) )THEN

      CALL MapToStage( iZ_B0, U0, U_IMEX )

      DO iS = 1, nStages

        IF( w_IM(iS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * w_IM(iS) * dU_IM(:,iS)

        END IF

        IF( w_EX(iS) .NE. Zero )THEN

          U_IMEX(:) = U_IMEX(:) + dt * w_EX(iS) * dU_EX(:,iS)

        END IF

      END DO

    END IF

    CALL MapFromStage( iZ_B1, U, U_IMEX )

    CALL ApplySlopeLimiter_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

    CALL ApplyPositivityLimiter_TwoMoment &
           ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

    ! --- Correction Step ---

    IF( ANY( [ alpha_IM, alpha_EX ] > Zero ) )THEN

      IF( alpha_EX > Zero )THEN

        CALL ComputeIncrement_Explicit &
               ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U, dU )

      END IF

      CALL ComputeCorrection_Implicit &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
               dt, alpha_EX, alpha_IM, GE, GX, U, dU )

      U(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
          iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),:,:) &
        = U(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
              iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),:,:) &
          + dt**2 * dU(:,iZ_B0(1):iZ_E0(1),iZ_B0(2):iZ_E0(2), &
                         iZ_B0(3):iZ_E0(3),iZ_B0(4):iZ_E0(4),:,:)

      CALL ApplySlopeLimiter_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

      CALL ApplyPositivityLimiter_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, GE, GX, U )

    END IF

  END SUBROUTINE Update_IMEX_RK


  SUBROUTINE InitializeStep_IMEX_RK( U )

    REAL(DP), INTENT(in) :: &
      U(1:,iZ_B1(1):,iZ_B1(2):,iZ_B1(3):,iZ_B1(4):,1:,1:)

    U0(1:nDOF, &
       iZ_B0(1):iZ_E0(1), &
       iZ_B0(2):iZ_E0(2), &
       iZ_B0(3):iZ_E0(3), &
       iZ_B0(4):iZ_E0(4), &
       1:nCR,1:nSpecies) &
    = U(1:nDOF, &
        iZ_B0(1):iZ_E0(1), &
        iZ_B0(2):iZ_E0(2), &
        iZ_B0(3):iZ_E0(3), &
        iZ_B0(4):iZ_E0(4), &
        1:nCR,1:nSpecies)

    U_IMEX = Zero
    dU_IM  = Zero
    dU_EX  = Zero

  END SUBROUTINE InitializeStep_IMEX_RK


  SUBROUTINE MapToStage( iZ_B, U_7D, U_1D )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4)
    REAL(DP), INTENT(in)  :: &
      U_7D(1:,iZ_B(1):,iZ_B(2):,iZ_B(3):,iZ_B(4):,1:,1:)
    REAL(DP), INTENT(out) :: &
      U_1D(1:)

    INTEGER :: i, iNode, iZ1, iZ2, iZ3, iZ4, iCR, iS

    i = 1
    DO iS = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iNode = 1, nDOF

      U_1D(i) = U_7D(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS)

      i = i + 1

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapToStage


  SUBROUTINE MapFromStage( iZ_B, U_7D, U_1D )

    INTEGER,  INTENT(in)  :: &
      iZ_B(4)
    REAL(DP), INTENT(out) :: &
      U_7D(1:,iZ_B(1):,iZ_B(2):,iZ_B(3):,iZ_B(4):,1:,1:)
    REAL(DP), INTENT(in)  :: &
      U_1D(1:)

    INTEGER :: i, iNode, iZ1, iZ2, iZ3, iZ4, iCR, iS

    i = 1
    DO iS = 1, nSpecies
    DO iCR = 1, nCR
    DO iZ4 = iZ_B0(4), iZ_E0(4)
    DO iZ3 = iZ_B0(3), iZ_E0(3)
    DO iZ2 = iZ_B0(2), iZ_E0(2)
    DO iZ1 = iZ_B0(1), iZ_E0(1)
    DO iNode = 1, nDOF

      U_7D(iNode,iZ1,iZ2,iZ3,iZ4,iCR,iS) = U_1D(i)

      i = i + 1

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE MapFromStage

END MODULE TimeSteppingModule_IMEX_RK
