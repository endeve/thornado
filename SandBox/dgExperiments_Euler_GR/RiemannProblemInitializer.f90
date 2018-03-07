MODULE RiemannProblemInitializer

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: RiemannProblemChoice

CONTAINS

  SUBROUTINE RiemannProblemChoice &
               ( D_L, V_L, P_L, D_R, V_R, P_R, &
                 xR, x_D, K, t, t_end, CFL, Gamma, iRP )
    
    REAL(DP), INTENT(out) :: D_L, V_L(3), P_L, D_R, V_R(3), P_R
    REAL(DP), INTENT(out) :: xR, x_D, CFL, t, t_end, Gamma
    INTEGER,  INTENT(in)  :: iRP
    INTEGER, INTENT(out)  :: K

    SELECT CASE( iRP )

      CASE( 0 )

        WRITE(*,*) 'Sods Shock Tube'

        D_L = 1.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0_DP

        D_R = 0.125_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 0.1_DP

        t     = 0.0_DP
        t_end = 0.2_DP
        CFL   = 0.1_DP
        xR    = 1.0d0
        x_D   = 0.5_DP
        K     = 128
        Gamma = 4.0_DP / 3.0_DP

      CASE( 1 )

        WRITE(*,*) 'Blast Wave 1 (Del Zanna & Bucciantini (2002))'
         
        D_L = 10.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.33d+1

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d-6

        t     = 0.0_DP
        t_end = 4.0d-1
        CFL   = 0.5_DP
        xR    = 1.0d0
        x_D   = 0.5_DP
        K     = 400
        Gamma = 5.0_DP / 3.0_DP

      CASE( 2 )

        WRITE(*,*) 'Blast Wave 2 (Del Zanna & Bucciantini (2002))'
         
        D_L = 1.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d+3

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d-2

        t     = 0.0_DP
        t_end = 3.5d-1
        CFL   = 0.5_DP
        xR    = 1.0d0
        x_D   = 0.5_DP
        K     = 400
        Gamma = 5.0_DP / 3.0_DP

      CASE( 3 )

        WRITE(*,*) 'Strong Blast Wave (Qin et al., (2016)))'
         
        D_L = 1.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d+4

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d-6

        t     = 0.0_DP
        t_end = 3.0d-1
        CFL   = 0.5_DP
        xR    = 1.0d0
        x_D   = 0.5_DP
        K     = 400
        Gamma = 4.0_DP / 3.0_DP

      CASE( 4 )

        WRITE(*,*) 'Shock Reflection (Del Zanna & Bucciantini (2002))'
         
        D_L = 1.0_DP
        V_L = [  0.99999_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d-2

        D_R = 1.0_DP
        V_R = [ -0.99999_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d-2

        t     = 0.0_DP
        t_end = 7.5d-1
        CFL   = 0.5_DP
        xR    = 2.0d0
        x_D   = 1.0_DP
        K     = 500
        Gamma = 5.0_DP / 3.0_DP

      CASE( 5 )

        WRITE(*,*) 'Contact Discontinuity 1 (Mignone & Bodo (2005))'
         
        D_L = 1.0_DP
        V_L = [ 0.9_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0d+0

        D_R = 1.0_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 1.0d+1

        t     = 0.0_DP
        t_end = 4.0d-1
        CFL   = 0.8_DP
        xR    = 1.0d0
        x_D   = 0.5_DP
        K     = 100
        Gamma = 4.0_DP / 3.0_DP

      CASE( 6 )

        WRITE(*,*) 'Contact Discontinuity 2 (Mignone & Bodo (2005))'
         
        D_L = 1.0d+0
        V_L = [ -0.6d+0, 0.0d+0, 0.0d+0 ]
        P_L = 1.0d+1

        D_R = 1.0d+1
        V_R = [ 0.5d+0, 0.0d+0, 0.0d+0 ]
        P_R = 2.0d+1

        t     = 0.0d+0
        t_end = 4.0d+0
        CFL   = 0.8d+0
        xR    = 1.0d0
        x_D   = 0.5d+0
        K     = 100
        Gamma = 4.0_DP / 3.0_DP

      CASE( 7 )

        WRITE(*,*) 'Non-zero transverse velocity 1 (Qin et al., (2016))'
         
        D_L = 1.0d+0
        V_L = [ 0.0d+0, 0.9d+0, 0.0d+0 ]
        P_L = 1.0d+3

        D_R = 1.0d+0
        V_R = [ 0.0d+0, 0.9d+0, 0.0d+0 ]
        P_R = 1.0d-2

        t     = 0.0d+0
        t_end = 6.0d-1
        CFL   = 0.15d+0
        xR    = 1.0d0
        x_D   = 0.5d+0
        K     = 400
        Gamma = 5.0_DP / 3.0_DP

      CASE( 8 )

        WRITE(*,*) 'Smooth Problem'

        D_L = 1.0_DP
        V_L = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_L = 1.0_DP

        D_R = 0.125_DP
        V_R = [ 0.0_DP, 0.0_DP, 0.0_DP ]
        P_R = 0.1_DP

        t     = 0.0_DP
        CFL   = 0.1_DP
        xR    = 1.0d0
        x_D   = 0.5_DP
        K     = 128
        t_end = 1.0d0 * CFL * xR / ( 1.0d0 * K )
        Gamma = 4.0_DP / 3.0_DP

     CASE DEFAULT

        WRITE(*,*) 'Invalid choice for iRP, exiting...'
        STOP

    END SELECT

  END SUBROUTINE RiemannProblemChoice

END MODULE RiemannProblemInitializer

