module mod_csp_vdiff
   use mod_csp_basic
   use mod_mpi_test
   implicit none
contains

!==============================================================================
   subroutine csp_vdiff_implicit(gTracer, KappaRX, recip_hFac, dT)
!==============================================================================
      implicit none
      integer :: i, k
      real(wp) :: a(nlpb, nk), b(nlpb, nk), c(nlpb, nk)
      real(wp) :: bet(nlpb, nk), gam(nlpb, nk), locTr(nlpb, nk)
      real(wp), intent(in) :: dT
      real(wp), intent(in) :: KappaRX(nlpb, nkp1), recip_hFac(nlpb, nk)
      real(wp), intent(inout) :: gTracer(nlpb, nk)

      !$ACC kernels create(a, b, c, bet, gam, locTr),   &
      !$ACC        present(KappaRX, recip_hFac, gTracer, recip_drC, recip_drF,  &
      !$ACC                kSurfC, nl, nlpb, nk, dT)

      ! -- Initialise
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            locTr(i, k) = 0.0_wp
         end do
      end do

      ! -- Old aLower
      !$ACC loop independent
      do i = 1, nlpb
         a(i, 1) = 0.0_wp
      end do

      !$ACC loop collapse(2) independent
      do k = 2, nk
         do i = 1, nlpb
            a(i, k) = -dT*recip_hFac(i, k)*recip_drF(k)*KappaRX(i, k)*recip_drC(k)
            if (recip_hFac(i, k - 1) .eq. 0.0_wp) a(i, k) = 0.0_wp
         end do
      end do

      ! -- Old aUpper
      !$ACC loop collapse(2) independent
      do k = 1, nk - 1
         do i = 1, nlpb
            c(i, k) = -dT*recip_hFac(i, k)*recip_drF(k)*KappaRX(i, k + 1)*recip_drC(k + 1)
            if (recip_hFac(i, k + 1) .eq. 0.0_wp) c(i, k) = 0.0_wp
         end do
      end do

      !$ACC loop independent
      do i = 1, nlpb
         c(i, nk) = 0.0_wp
      end do

      ! -- Old aCenter
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            b(i, k) = 1.0_wp - (a(i, k) + c(i, k))
         end do
      end do

      ! -- Old and new gam, bet are the same
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            bet(i, k) = 1.0_wp
            gam(i, k) = 0.0_wp
         end do
      end do

      ! -- Beginning of forward sweep (top level)
      !$ACC loop independent
      do i = 1, nlpb
         if (b(i, 1) .ne. 0.0_wp) bet(i, 1) = 1.0_wp/b(i, 1)
      end do

      ! -- Middle of forward sweep
      !$ACC loop independent
      do i = 1, nlpb
         !$ACC loop seq
         do k = 2, nk
            gam(i, k) = c(i, k - 1)*bet(i, k - 1)
            if ((b(i, k) - a(i, k)*gam(i, k)) .ne. 0.0_wp) then
               bet(i, k) = 1.0_wp/(b(i, k) - a(i, k)*gam(i, k))
            end if
         end do
      end do

      !$ACC loop independent
      do i = 1, nlpb
         locTr(i, 1) = gTracer(i, 1)*bet(i, 1)
      end do

      !$ACC loop independent
      do i = 1, nlpb
         !$ACC loop seq
         do k = 2, nk
            locTr(i, k) = bet(i, k)*(gTracer(i, k) - a(i, k)*locTr(i, k - 1))
         end do
      end do

      ! -- Backward sweep
      !$ACC loop independent
      do i = 1, nlpb
         !$ACC loop seq
         do k = nk - 1, 1, -1
            locTr(i, k) = locTr(i, k) - gam(i, k + 1)*locTr(i, k + 1)
         end do
      end do

      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            gTracer(i, k) = locTr(i, k)
         end do
      end do
      !$ACC end kernels

   end subroutine csp_vdiff_implicit

!==============================================================================
   subroutine csp_vdiff_implicit_gls(d, KappaRX, workc)
!==============================================================================
      implicit none
      integer :: i, k, ibot, ibotp1
      real(wp) :: a(nlpb, nk), b(nlpb, nk), c(nlpb, nk)
      real(wp) :: u(nlpb, nk), l(nlpb, nk), y(nlpb, nk)
      real(wp), intent(in) :: KappaRX(nlpb, nk), workc(nlpb, nkp1)
      real(wp), intent(inout) :: d(nlpb, nkp1)

      !$ACC kernels create(a,b,c,u,l,y),   &
      !$ACC         present(KappaRX, workc, d,   &
      !$ACC                 nl,nlpb,nk,nkp1,dTtracer,recip_hFacC,recip_drC,recip_drF,kSurfC)

      ! -- Initialise
      ! -- Old aLower
      !$ACC loop independent
      do i = 1, nlpb
         a(i, 1) = 0.0_wp
      end do

      !$ACC loop collapse(2)
      do k = 2, nk
         do i = 1, nlpb
            a(i, k) = -dTtracer*recip_hFacC(i, k - 1)*recip_drC(k)*KappaRX(i, k - 1)*recip_drF(k - 1)
         end do
      end do

      ! -- Old aUpper
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            c(i, k) = -dTtracer*recip_hFacC(i, k)*recip_drC(k)*KappaRX(i, k)*recip_drF(k)
         end do
      end do

      ! -- Old aCenter
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            b(i, k) = 1.0_wp - (a(i, k) + c(i, k)) + workc(i, k)*dTtracer
         end do
      end do

      ! -- set surface special boundary
      ! -- at surface level and level below, tke and psi not allowed transport downward,
      ! -- but level below below allow recieve tke and psi from level below.
      ! -- and tke and psi at the two level only equal boundary condition value
      !$ACC loop independent
      do i = 1, nlpb
         a(i, nk) = 0.0_wp
         c(i, nk) = 0.0_wp
         b(i, nk) = 1.0_wp
      end do

      ! -- set bottom boundary, same purpose as surface, only level up bottom allow transport
      ! -- tke and psi upward, and bottom and bottom level up level only equal boundary condition value.
      !$ACC loop independent private(ibot,ibotp1)
      do i = 1, nlpb
         ibot = min(nk, kSurfC(i))
         ibotp1 = min(nk, ibot + 1)
         a(i, ibot) = 0.0_wp
         c(i, ibot) = 0.0_wp
         b(i, ibot) = 1.0_wp

         a(i, ibotp1) = 0.0_wp
         c(i, ibotp1) = 0.0_wp
         b(i, ibotp1) = 1.0_wp
      end do

      ! -- set LU decompose
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            u(i, k) = b(i, k)
            l(i, k) = 0.0_wp
         end do
      end do

      !$ACC loop collapse(2)
      do k = 2, nk
         do i = 1, nlpb
            l(i, k) = a(i, k)/u(i, k)
            u(i, k) = b(i, k) - l(i, k)*c(i, k - 1)
         end do
      end do

      ! -- solve Ly = d
      !$ACC loop collapse(2) independent
      do k = 1, nk
         do i = 1, nlpb
            y(i, k) = d(i, k)
         end do
      end do

      !$ACC loop independent
      do i = 1, nlpb
         !$ACC loop seq
         do k = 2, nk
            y(i, k) = d(i, k) - l(i, k)*y(i, k - 1)
         end do
      end do

      ! -- solve Ux = y
      !$ACC loop independent
      do i = 1, nlpb
         d(i, nk) = y(i, nk)/u(i, nk)
      end do

      !$ACC loop independent
      do i = 1, nlpb
         !$ACC loop seq
         do k = nk - 1, 1, -1
            d(i, k) = (y(i, k) - c(i, k)*d(i, k + 1))/u(i, k)
         end do
      end do
      !$ACC end kernels

   end subroutine csp_vdiff_implicit_gls

end module mod_csp_vdiff
