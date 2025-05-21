module mod_csp_thdyn
   use mod_csp_init
   use mod_csp_vdiff
   use mod_csp_thdyn_hdiff
   use mod_csp_thdyn_adv
   use mod_csp_thdyn_force
   use mod_csp_bdy
   implicit none
contains
!==============================================================================
   subroutine csp_thdyn
!==============================================================================
      implicit none
      integer :: i, i_chk, k, itracer, w, s, e, n

      do itracer = 1, ntracer
         !$ACC kernels present(tracer, tFld, gTracer), copyin(itracer)
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, nlpb
               tracer(i, k) = tFld(i, k, itracer)
               gTracer(i, k) = 0.0_wp
            end do
         end do
         !$ACC end kernels

         if (itracer .eq. 1) then
            call csp_thdyn_force_temp
         else if (itracer .eq. 2) then
            call csp_thdyn_force_salt
         else
            continue
         end if

         call csp_thdyn_adv(itracer)

         call csp_thdyn_hdiff(itracer)

         call csp_thdyn_forward(itracer)

         if (itracer .eq. 2) then
            !$ACC kernels present(tracer)
            !$ACC loop collapse(2)
            do k = 1, nk
               do i = 1, nlpb
                  if (tracer(i, k) .le. 0.1_wp) then
                     tracer(i, k) = 0.1_wp
                  end if
               end do
            end do
            !$ACC end kernels
         end if

         ! update boundary condition
         if (ln_bdy .and. mpi_check_bdy) then
            if (itracer .eq. 1) then
               call csp_bdy_orlanski_3d(tracer, bdy_t(:,:,4), tFld(:,:,1), maskC, tw, ts, tn, te, dTtracer, real(fbt, wp))
            else if (itracer .eq. 2) then
               call csp_bdy_orlanski_3d(tracer, bdy_s(:,:,4), tFld(:,:,2), maskC, tw, ts, tn, te, dTtracer, real(fbs, wp))
            else
               continue
            end if
         end if

         CALL mpi_data_exchange(mpi_sendbuf_2d_nk_f, tracer, nk)

         !$ACC kernels present(tFld,tracer)
         !$ACC loop collapse(2)
         do k = 1, nk
            do i = 1, nlpb
               tFld(i, k, itracer) = tracer(i, k)
            end do
         end do
         !$ACC end kernels

      end do

   end subroutine csp_thdyn

!==============================================================================
   subroutine csp_thdyn_forward(its)
!==============================================================================
      implicit none
      integer, intent(in) :: its
      integer :: i, k
      real(wp) :: temp(nlpb, nk)
      real(dp) :: sumtest1, sumtest2, dsumtest, smass

      !$ACC data present(nk,gTracerPast,rStarExpC,gTracer,tracer,dTtracer,  &
      !$ACC      maskC,kappaRT,recip_hFacC, tFld_vdiff, sFld_vdiff), &
      !$ACC      create(temp)
      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nk
         ! gTracerPast need reset by new hFac, need future revised
         do i = 1, nlpb
            gTracerPast(i, k, its) = gTracerPast(i, k, its)*rStarExpC(i)
         end do
      end do
      !$ACC end kernels

      do k = 1, nk
         call csp_misc_adams_bashforth2(gTracer(:, k), gTracerPast(:, k, its), FirstTracer(its))
      end do

      FirstTracer(its) = .false.

      ! get guess field
      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            gTracer(i, k) = (tracer(i, k) + dTtracer*gTracer(i, k))*maskC(i, k)
         end do
      end do
      !$ACC end kernels

# if defined (EXDIAG)
      !$ACC kernels
      !$ACC loop collapse(2) independent
      do k = 1, nk
         do i = 1, loc_nlpb
            temp(i, k) = gTracer(i, k)
         end do
      end do
      !$ACC end kernels
# endif

      call csp_vdiff_implicit(gTracer, kappaRT, recip_hFacC, dTtracer)

# if defined (EXDIAG)
      if (its .eq. 1) then
         !$ACC kernels
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, loc_nlpb
               tFld_vdiff(i, k) = tFld_vdiff(i, k) + (gTracer(i, k) - temp(i, k))
            end do
         end do
         !$ACC end kernels
      end if

      if (its .eq. 2) then
         !$ACC kernels
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, loc_nlpb
               sFld_vdiff(i, k) = sFld_vdiff(i, k) + (gTracer(i, k) - temp(i, k))
            end do
         end do
         !$ACC end kernels
      end if
# endif

      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            tracer(i, k) = gTracer(i, k)
         end do
      end do
      !$ACC end kernels

      !$ACC end data

   end subroutine csp_thdyn_forward

end module mod_csp_thdyn
