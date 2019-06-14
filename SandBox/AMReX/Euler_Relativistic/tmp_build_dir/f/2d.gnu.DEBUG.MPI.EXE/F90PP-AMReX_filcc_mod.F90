/Users/dunhamsj/Research/AMReX/amrex/Src/Base/AMReX_filcc_mod.F90:2721:50: warning: missing terminating ' character [-Winvalid-pp-token]
  !     Set these to invalid values, they shouldn't be used if not reset
                                                 ^
1 warning generated.

module amrex_filcc_module

  use amrex_fort_module, only : amrex_real, amrex_spacedim, amrex_get_loop_bounds
  use amrex_bc_types_module
  use amrex_constants_module

  implicit none

  interface amrex_filcc
     module procedure amrex_filcc_1
     module procedure amrex_filcc_n
  end interface amrex_filcc

  private
  public :: amrex_filcc, amrex_fab_filcc, amrex_filccn, amrex_hoextraptocc




  public :: amrex_hoextraptocc_2d






  public :: filccn


contains

  subroutine amrex_filcc_n(q,qlo,qhi,domlo,domhi,dx,xlo,bclo,bchi)
    integer, intent(in) :: qlo(4), qhi(4)
    integer, dimension(amrex_spacedim), intent(in) :: domlo, domhi
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bclo(amrex_spacedim,*), bchi(amrex_spacedim,*)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qlo(4):qhi(4))
    integer :: i, bc(amrex_spacedim,2)
    do i = qlo(4), qhi(4)
       bc(:,1) = bclo(:,i)
       bc(:,2) = bchi(:,i)
       call amrex_filccn(qlo(1:3), qhi(1:3), q(:,:,:,i), qlo(1:3), qhi(1:3), 1, &
            domlo, domhi, dx, xlo, bc);
    end do
  end subroutine amrex_filcc_n



  subroutine amrex_filcc_1(q,qlo1,qlo2,qhi1,qhi2,domlo,domhi,dx,xlo,bc)
    integer, intent(in) :: qlo1,qlo2,qhi1,qhi2,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim, 2)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1,qlo2:qhi2)
    integer :: q_lo(3), q_hi(3)
    q_lo = [qlo1,qlo2,0]
    q_hi = [qhi1,qhi2,0]
    call amrex_filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc);
  end subroutine amrex_filcc_1



  subroutine amrex_fab_filcc (q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc) &
       bind(c, name='amrex_fab_filcc')

    implicit none

    integer, intent(in) :: qlo(3), qhi(3), nq
    integer, dimension(amrex_spacedim), intent(in) :: domlo, domhi
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim,2,nq)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nq)

    integer :: lo(3), hi(3)

    call amrex_get_loop_bounds(lo, hi, qlo, qhi)

    call amrex_filccn(lo, hi, q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc)

  end subroutine amrex_fab_filcc


  subroutine filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)
    implicit none
    integer,          intent(in   ) :: lo(3), hi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: ncomp
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)
    call amrex_filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)
  end subroutine filccn


  subroutine amrex_filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)

    implicit none

    integer,          intent(in   ) :: lo(3), hi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: ncomp
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)

    integer :: ilo, ihi, jlo, jhi, klo, khi
    integer :: is, ie, js, je, ks, ke
    integer :: i, j, k, n
    integer :: imin, imax, jmin, jmax, kmin, kmax

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))
    ilo = domlo(1)
    ihi = domhi(1)


    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))
    jlo = domlo(2)
    jhi = domhi(2)









    do n = 1, ncomp

       if (lo(1) < ilo) then
          imin = lo(1)
          imax = ilo-1

          if (bc(1,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(1,1,n) .eq. amrex_bc_foextrap) then
             
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,1,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

                      if (i < ilo - 1) then
                         q(i,j,k,n) = q(ilo,j,k,n)
                      else if (i == ilo - 1) then
                         if (ilo+2 <= ie) then
                            q(i,j,k,n) = eighth * (15*q(ilo,j,k,n) - 10*q(ilo+1,j,k,n) + 3*q(ilo+2,j,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(ilo,j,k,n) - q(ilo+1,j,k,n))
                         end if
                      end if

                   end do
                end do
             end do
             
          else if (bc(1,1,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,1,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(1) > ihi) then
          imin = ihi+1
          imax = hi(1)

          if (bc(1,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(1,2,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,2,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

                      if (i > ihi + 1) then
                         q(i,j,k,n) = q(ihi,j,k,n)
                      else if (i == ihi + 1) then
                         if (ihi-2 >= is) then
                            q(i,j,k,n) = eighth * (15*q(ihi,j,k,n) - 10*q(ihi-1,j,k,n) + 3*q(ihi-2,j,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(ihi,j,k,n) - q(ihi-1,j,k,n))
                         end if
                      end if

                   end do
                end do
             end do

          else if (bc(1,2,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,2,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do
             
          end if

       end if



       if (lo(2) < jlo) then
          jmin = lo(2)
          jmax = jlo-1

          if (bc(2,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(2,1,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jlo,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,1,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)

                      if (j < jlo - 1) then
                         q(i,j,k,n) = q(i,jlo,k,n)
                      else if (j == jlo - 1) then
                         if (jlo+2 <= je) then
                            q(i,j,k,n) = eighth * (15*q(i,jlo,k,n) - 10*q(i,jlo+1,k,n) + 3*q(i,jlo+2,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(i,jlo,k,n) - q(i,jlo+1,k,n))
                         end if
                      end if

                   end do
                end do
             end do

          else if (bc(2,1,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jlo+(jlo-j)-1,k,n)
                   end do
                end do
             end do

          else if (bc(2,1,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,jlo+(jlo-j)-1,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(2) > jhi) then
          jmin = jhi+1
          jmax = hi(2)

          if (bc(2,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(2,2,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jhi,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,2,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)

                      if (j > jhi + 1) then
                         q(i,j,k,n) = q(i,jhi,k,n)
                      else if (j == jhi + 1) then
                         if (jhi-2 >= js) then
                            q(i,j,k,n) = eighth * (15*q(i,jhi,k,n) - 10*q(i,jhi-1,k,n) + 3*q(i,jhi-2,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(i,jhi,k,n) - q(i,jhi-1,k,n))
                         end if
                      end if
                      
                   end do
                end do
             end do

          else if (bc(2,2,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jhi-(j-jhi)+1,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,2,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,jhi-(j-jhi)+1,k,n)
                   end do
                end do
             end do
             
          end if

       end if










       ! Now take care of the higher contributions

       !
       ! First correct the i-j edges and all corners
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(2,1,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. lo(2) < jlo) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)

             i = ilo-1
             j = jlo-1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)
                      
                   if (jlo+2 <= je) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,jlo,k,n) - 10*q(ilo-1,jlo+1,k,n) + 3*q(ilo-1,jlo+2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,jlo,k,n) - q(ilo-1,jlo+1,k,n))
                   end if
                   
                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,jlo-1,k,n) + &
                           half * eighth * (15*q(ilo,jlo-1,k,n) - 10*q(ilo+1,jlo-1,k,n) + 3*q(ilo+2,jlo-1,k,n))
                   else
                      q(i,j,k,n) = q(ilo-1,jlo-1,k,n) + half * half * (3*q(ilo,jlo-1,k,n) - q(ilo+1,jlo-1,k,n))
                   end if
                   

                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(2,2,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. hi(2) > jhi) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)

             i = ilo-1
             j = jhi+1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)

                   if (jhi-2 >= js) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,jhi,k,n) - 10*q(ilo-1,jhi-1,k,n) + 3*q(ilo-1,jhi-2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,jhi,k,n) - q(ilo-1,jhi-1,k,n))
                   end if
                   
                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,jhi+1,k,n) + &
                           half * eighth * (15*q(ilo,jhi+1,k,n) - 10*q(ilo+1,jhi+1,k,n) + 3*q(ilo+2,jhi+1,k,n))
                   else
                      q(i,j,k,n) = q(ilo-1,jhi+1,k,n) + half * half * (3*q(ilo,jhi+1,k,n) - q(ilo+1,jhi+1,k,n))
                   end if



                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(2,1,n) .eq. amrex_bc_hoextrap) then

          if (hi(1) > ihi .and. lo(2) < jlo) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)

             i = ihi+1
             j = jlo-1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)

                   if (jlo+2 <= je) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,jlo,k,n) - 10*q(ihi+1,jlo+1,k,n) + 3*q(ihi+1,jlo+2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,jlo,k,n) - q(ihi+1,jlo+1,k,n))
                   end if
                   
                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,jlo-1,k,n) + &
                           half * eighth * (15*q(ihi,jlo-1,k,n) - 10*q(ihi-1,jlo-1,k,n) + 3*q(ihi-2,jlo-1,k,n))
                   else
                      q(i,j,k,n) = q(ihi+1,jlo-1,k,n) + half * half * (3*q(ihi,jlo-1,k,n) - q(ihi-1,jlo-1,k,n))
                   end if
                   


                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(2,2,n) .eq. amrex_bc_hoextrap) then

          if (hi(1) > ihi .and. hi(2) > jhi) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)

             i = ihi+1
             j = jhi+1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)
                   
                   if (jhi-2 >= js) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,jhi,k,n) - 10*q(ihi+1,jhi-1,k,n) + 3*q(ihi+1,jhi-2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,jhi,k,n) - q(ihi+1,jhi-1,k,n))
                   end if
                   
                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,jhi+1,k,n) + &
                           half * eighth * (15*q(ihi,jhi+1,k,n) - 10*q(ihi-1,jhi+1,k,n) + 3*q(ihi-2,jhi+1,k,n))
                   else
                      q(i,j,k,n) = q(ihi+1,jhi+1,k,n) + half * half * (3*q(ihi,jhi+1,k,n) - q(ihi-1,jhi+1,k,n))
                   end if
                   

                   
                end do
             end if
          end if
       end if




    end do

  end subroutine amrex_filccn



  subroutine amrex_hoextraptocc (q, qlo, qhi, domlo, domhi, dx, xlo) &
       bind(c,name='amrex_hoextraptocc')
    integer, intent(in) :: qlo(3), qhi(3), domlo(*), domhi(*)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    real(amrex_real), intent(in) :: dx(*), xlo(*)




    call amrex_hoextraptocc_2d(q,qlo(1),qlo(2),qhi(1),qhi(2),domlo,domhi,dx,xlo)

  end subroutine amrex_hoextraptocc





subroutine amrex_hoextraptocc_2d(q,q_l1,q_l2,q_h1,q_h2,domlo,domhi,dx,xlo)

  use amrex_fort_module
  use amrex_constants_module

  implicit none

  integer    q_l1, q_l2, q_h1, q_h2
  integer    domlo(2), domhi(2)
  real(amrex_real)     xlo(2), dx(2)
  real(amrex_real)     q(q_l1:q_h1,q_l2:q_h2)

  integer    nlft, nrgt, nbot, ntop
  integer    ilo, ihi, jlo, jhi
  integer    i, j
  integer    is, ie, js, je

  nlft = max(0,domlo(1)-q_l1)
  nrgt = max(0,q_h1-domhi(1))
  nbot = max(0,domlo(2)-q_l2)
  ntop = max(0,q_h2-domhi(2))

  is = max(q_l1,domlo(1))
  ie = min(q_h1,domhi(1))
  js = max(q_l2,domlo(2))
  je = min(q_h2,domhi(2))

  !
  !     Set these to invalid values, they shouldn't be used if not reset
  !
  ilo = -10
  jlo = -10
  ihi = 100000000
  jhi = 100000000

  !
  !     First fill sides.
  !
  if (nlft .gt. 0) then
     ilo = domlo(1)
     do i = 2, nlft
        do j = q_l2, q_h2
           q(ilo-i,j) = q(ilo,j) 
        end do
     end do
     if (ilo+2 .le. ie) then 
        do j = q_l2, q_h2
           q(ilo-1,j) = 3*q(ilo,j) - 3*q(ilo+1,j) + q(ilo+2,j)
        end do
     else 
        do j = q_l2, q_h2
           q(ilo-1,j) = 2*q(ilo,j) - q(ilo+1,j)
        end do
     end if
  end if

  if (nrgt .gt. 0) then
     ihi = domhi(1)
     do i = 2, nrgt
        do j = q_l2, q_h2
           q(ihi+i,j) = q(ihi,j)
        end do
     end do
     if (ihi-2 .ge. is) then
        do j = q_l2, q_h2
           q(ihi+1,j) = 3*q(ihi,j) - 3*q(ihi-1,j) + q(ihi-2,j)
        end do
     else
        do j = q_l2, q_h2
           q(ihi+1,j) = 2*q(ihi,j) - q(ihi-1,j)
        end do
     end if
  end if

  if (nbot .gt. 0) then
     jlo = domlo(2)
     do j = 2, nbot
        do i = q_l1, q_h1
           q(i,jlo-j) = q(i,jlo)
        end do
     end do
     if (jlo+2 .le. je) then
        do i = q_l1, q_h1
           q(i,jlo-1) = 3*q(i,jlo) - 3*q(i,jlo+1) + q(i,jlo+2)
        end do
     else
        do i = q_l1, q_h1
           q(i,jlo-1) = 2*q(i,jlo) - q(i,jlo+1)
        end do
     end if
  end if

  if (ntop .gt. 0) then
     jhi = domhi(2)
     do j = 2, ntop
        do i = q_l1, q_h1
           q(i,jhi+j) = q(i,jhi)
        end do
     end do
     if (jhi-2 .ge. js) then
        do i = q_l1, q_h1
           q(i,jhi+1) = 3*q(i,jhi) - 3*q(i,jhi-1) + q(i,jhi-2)
        end do
     else
        do i = q_l1, q_h1
           q(i,jhi+1) = 2*q(i,jhi) - q(i,jhi-1)
        end do
     end if
  end if

  if ((nlft .gt. 0) .and. (nbot .gt. 0)) then
     if (jlo+2 .le. je) then
        q(ilo-1,jlo-1) = half * &
             (3*q(ilo-1,jlo) - 3*q(ilo-1,jlo+1) + q(ilo-1,jlo+2))
     else
        q(ilo-1,jlo-1) = half * (2*q(ilo-1,jlo) - q(ilo-1,jlo+1))
     end if

     if (ilo+2 .le. ie) then 
        q(ilo-1,jlo-1) =  q(ilo-1,jlo-1) + half * &
             (3*q(ilo,jlo-1) - 3*q(ilo+1,jlo-1) + q(ilo+2,jlo-1)) 
     else
        q(ilo-1,jlo-1) =  q(ilo-1,jlo-1) + half * &
             (2*q(ilo,jlo-1) - q(ilo+1,jlo-1))
     end if
  end if

  if ((nlft .gt. 0) .and. (ntop .gt. 0)) then 
     if (jhi-2 .ge. js) then 
        q(ilo-1,jhi+1) = half * &
             (3*q(ilo-1,jhi) - 3*q(ilo-1,jhi-1) + q(ilo-1,jhi-2))
     else
        q(ilo-1,jhi+1) = half * (2*q(ilo-1,jhi) - q(ilo-1,jhi-1))
     end if

     if (ilo+2 .le. ie) then 
        q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + half * &
             (3*q(ilo,jhi+1) - 3*q(ilo+1,jhi+1) + q(ilo+2,jhi+1))
     else
        q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + half * &
             (2*q(ilo,jhi+1) - q(ilo+1,jhi+1))
     end if
  end if

  if ((nrgt .gt. 0) .and. (nbot .gt. 0)) then 
     if (jlo+2 .le. je) then 
        q(ihi+1,jlo-1) = half * &
             (3*q(ihi+1,jlo) - 3*q(ihi+1,jlo+1) + q(ihi+1,jlo+2))
     else
        q(ihi+1,jlo-1) = half * (2*q(ihi+1,jlo) - q(ihi+1,jlo+1))
     end if

     if (ihi-2 .ge. is) then 
        q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + half * &
             (3*q(ihi,jlo-1) - 3*q(ihi-1,jlo-1) + q(ihi-2,jlo-1))
     else
        q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + half * &
             (2*q(ihi,jlo-1) - q(ihi-1,jlo-1))
     end if
  end if

  if ((nrgt .gt. 0) .and. (ntop .gt. 0)) then 
     if (jhi-2 .ge. js) then 
        q(ihi+1,jhi+1) = half * &
             (3*q(ihi+1,jhi) - 3*q(ihi+1,jhi-1) + q(ihi+1,jhi-2))
     else
        q(ihi+1,jhi+1) = half * (2*q(ihi+1,jhi) - q(ihi+1,jhi-1))
     end if

     if (ihi-2 .ge. is) then 
        q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + half * &
             (3*q(ihi,jhi+1) - 3*q(ihi-1,jhi+1) + q(ihi-2,jhi+1))
     else
        q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + half * &
             (2*q(ihi,jhi+1) - q(ihi-1,jhi+1))
     end if
  end if

end subroutine amrex_hoextraptocc_2d


end module amrex_filcc_module

