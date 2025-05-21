module mod_csp_thdyn_hdiff
   use mod_csp_init
   implicit none
contains

!==============================================================================
   subroutine csp_thdyn_hdiff(its)
!==============================================================================
      implicit none
      integer, intent(in) :: its
      integer :: i, k, w, s, uel, uei, uep, vnl, vni, vnp
      real(wp) :: xyA(nlpb, nk, 2), dTdl(nlpb, nk, 2), &
                  del2T(nlpb, nk), dfl(nlpb, nk, 2), temp(nlpb, nk)

      !$ACC data present(dyZ,dxZ,drF,hFacW,hFacS,maskS,maskW,recip_dxC,recip_dyC,  &
      !$ACC      tracer,ue,vn,recip_rAc,recip_drF,recip_hFacC,maskC,  &
      !$ACC      diffKhu,diffKhv,diffK4u,diffK4v,tw,ts, tFld_hdiff, sFld_hdiff, dTtracer),  &
      !$ACC      create(xyA,dTdl,del2T,dfl,temp)

      call csp_thdyn_hdiff_coef

      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            xyA(i, k, iu) = dyZ(i)*drF(k)*hFacW(i, k)*maskW(i, k)
            xyA(i, k, iv) = dxZ(i)*drF(k)*hFacS(i, k)*maskS(i, k)
         end do
      end do
      !$ACC end kernels

      !$ACC kernels
      !$ACC loop collapse(2) private(w,s)
      do k = 1, nk
         do i = 1, nlpb
            w = tw(i)
            s = ts(i)
            dTdl(i, k, iu) = xyA(i, k, iu)*recip_dxC(i)*(tracer(i, k) - tracer(w, k))
            dTdl(i, k, iv) = xyA(i, k, iv)*recip_dyC(i)*(tracer(i, k) - tracer(s, k))
         end do
      end do
      !$ACC end kernels

      !$ACC kernels
      !$ACC loop collapse(2) private(uel,uei,uep,vnl,vni,vnp)
      do k = 1, nk
         do i = 1, nlpb
            uel = ue(i, 1)
            uei = ue(i, 2)
            uep = ue(i, 3)
            vnl = vn(i, 1)
            vni = vn(i, 2)
            vnp = vn(i, 3)
            del2T(i, k) = recip_rAc(i)*recip_drF(k)*recip_hFacC(i, k)*maskC(i, k) &
                          *((dTdl(uel, k, uei)*uep - dTdl(i, k, iu)) + (dTdl(vnl, k, vni)*vnp - dTdl(i, k, iv)))
         end do
      end do
      !$ACC end kernels

      if (harmonicT) then
         !$ACC kernels loop collapse(2)
         do k = 1, nk
            do i = 1, nlpb
               dfl(i, k, iu) = -diffKhu(i, k)*dTdl(i, k, iu)
               dfl(i, k, iv) = -diffKhv(i, k)*dTdl(i, k, iv)
            end do
         end do
         !$ACC end kernels
      else
         !$ACC kernels loop collapse(2)
         do k = 1, nk
            do i = 1, nlpb
               dfl(i, k, iu) = 0.0_wp
               dfl(i, k, iv) = 0.0_wp
            end do
         end do
         !$ACC end kernels
      end if

      if (biharmonicT) then
         !$ACC kernels loop collapse(2)
         do k = 1, nk
            do i = 1, nlpb
               w = tw(i)
               s = ts(i)
               dfl(i, k, iu) = dfl(i, k, iu) + diffK4u(i, k)*xyA(i, k, iu)*recip_dxC(i)*(del2T(i, k) - del2T(w, k))
               dfl(i, k, iv) = dfl(i, k, iv) + diffK4v(i, k)*xyA(i, k, iv)*recip_dyC(i)*(del2T(i, k) - del2T(s, k))
            end do
         end do
         !$ACC end kernels
      end if

      !$ACC kernels loop collapse(2) private(uel,uei,uep,vnl,vni,vnp)
      do k = 1, nk
         do i = 1, nlpb
            uel = ue(i, 1)
            uei = ue(i, 2)
            uep = ue(i, 3)
            vnl = vn(i, 1)
            vni = vn(i, 2)
            vnp = vn(i, 3)
            temp(i, k) = - recip_rAc(i)*recip_drF(k)*recip_hFacC(i, k)*maskC(i, k) &
                            *((dfl(uel, k, uei)*uep - dfl(i, k, iu)) + (dfl(vnl, k, vni)*vnp - dfl(i, k, iv)))
         end do
      end do
      !$ACC end kernels

      !$ACC kernels loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            gTracer(i, k) = gTracer(i, k) + temp(i, k)
         end do
      end do
      !$ACC end kernels

# if defined (EXDIAG)
      if (its .eq. 1) then
         !$ACC kernels
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, loc_nlpb
               tFld_hdiff(i, k) = tFld_hdiff(i, k) + temp(i, k)*dTtracer
            end do
         end do
         !$ACC end kernels
      end if

      if (its .eq. 2) then
         !$ACC kernels
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, loc_nlpb
               sFld_hdiff(i, k) = sFld_hdiff(i, k) + temp(i, k)*dTtracer
            end do
         end do
         !$ACC end kernels
      end if
# endif

      !$ACC end data

   end subroutine csp_thdyn_hdiff

!==============================================================================
   subroutine csp_thdyn_hdiff_coef
!==============================================================================
      implicit none
      integer :: i, k
      real(wp) :: r1_12

      r1_12 = 1.0_wp/12.0_wp

      if (harmonicT) then
         !$ACC kernels present(uFld,vFld,dxC,dyC,diffKhu,diffKhv),copyin(r1_12)
         !$ACC loop collapse(2)
         do k = 1, nk
            do i = 1, nlpb
               diffKhu(i, k) = 0.1_wp*abs(uFld(i, k))*dxC(i)*r1_12
               diffKhv(i, k) = 0.1_wp*abs(vFld(i, k))*dyC(i)*r1_12
            end do
         end do
         !$ACC end kernels
      end if

      if (biharmonicT) then
         !$ACC kernels present(uFld,vFld,dxC,dyC,diffK4u,diffK4v),copyin(r1_12)
         !$ACC loop collapse(2)
         do k = 1, nk
            do i = 1, nlpb
               diffK4u(i, k) = abs(uFld(i, k))*r1_12*dxC(i)**3
               diffK4v(i, k) = abs(vFld(i, k))*r1_12*dyC(i)**3
            end do
         end do
         !$ACC end kernels
      end if

   end subroutine csp_thdyn_hdiff_coef

end module mod_csp_thdyn_hdiff
