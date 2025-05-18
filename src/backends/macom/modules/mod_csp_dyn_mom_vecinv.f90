module mod_csp_dyn_mom_vecinv
   use mod_csp_misc
   use mod_csp_dyn_phi
   use mod_csp_dyn_hdiss
   use mod_csp_dyn_drag
   use mod_csp_dyn_forward
   use mod_csp_dyn_vort
   implicit none

contains

!==============================================================================
   subroutine csp_dyn_mom_vecinv
!==============================================================================
      implicit none
      integer :: i, k, bj, k1, k2, dk
      real(wp) :: phiSurfX(nlpb), phiSurfY(nlpb)
      real(dp) :: sshtemp, sshgloavg, rhsMax, rhsMin, mpi_tmp

      !$ACC kernels present(uBarCouple, vBarCouple)
      !$ACC loop independent
      do i = 1, nlpb
         uBarCouple(i) = 0.0_wp
         vBarCouple(i) = 0.0_wp
      end do
      !$ACC end kernels

      call csp_dyn_drag_bottom

      if (ln_bous) then
         k1 = nk
         k2 = 1
         dk = -1
      else
         k1 = 1
         k2 = nk
         dk = 1
      end if

      do k = k1, k2, dk

         call csp_dyn_mom_update_hFacZ(k)

         call csp_dyn_dphi_hyd(k)

         call csp_misc_hdiv(uFld(1:nlpb, k), vFld(1:nlpb, k), hDiv(1:nlpb))

         call csp_dyn_hrelvort(uFld(1:nlpb, k), vFld(1:nlpb, k), vort3(1:nlpbz))

         !$ACC kernels present(vort3,maskZ,omega3,fCoriG),  &
         !$ACC         copyin(k)
         !$ACC loop independent
         do i = 1, nlpbz
            vort3(i) = vort3(i)*maskZ(i, k)
            omega3(i) = fCoriG(i) + vort3(i)
         end do
         !$ACC end kernels

         call csp_misc_del2uv(k)

         call csp_dyn_vort_een(k)

         call csp_dyn_mom_vecinv_bernoulli(k)

         call csp_dyn_hdiss_coef(k)

         call csp_dyn_hdiss(k)

         call csp_dyn_drag_side(k)

         call csp_dyn_mom_vecinv_vshear(k)

         call csp_dyn_forward_guess(k)

      end do

      !$ACC kernels present(uBar,vBar,uFld,vFld,drF,maskW,maskS,recip_RcolW,recip_RcolS,h0FacW,h0FacS)
      !$ACC loop independent
      do i = 1, nlpb
         uBar(i) = 0.0_wp
         vBar(i) = 0.0_wp
         !$ACC loop seq
         do k = 1, nk
            uBar(i) = uBar(i) + uFld(i, k)*drF(k)*maskW(i,k)*h0FacW(i,k)*recip_RcolW(i)
            vBar(i) = vBar(i) + vFld(i, k)*drF(k)*maskS(i,k)*h0FacS(i,k)*recip_RcolS(i)
         end do
      end do
      !$ACC end kernels

      if (ln_bdy .and. mpi_check_bdy) then
         !$ACC kernels present(bdy_uBar,bdy_vBar,bdy_pbt,pbt_base,mpi_idx_bdy, &
         !$ACC                 bdy_u,bdy_v,drF,recip_RcolW,recip_RcolS,bdy_uBcl, &
         !$ACC                 bdy_vBcl,tide_u,tide_v,tide_ssh,gravityRau0,h0FacW,h0FacS)
         !$ACC loop independent private(bj)
         do i = 1, nlbdy
            bj = mpi_idx_bdy(i)
            bdy_pbt(i,4) = bdy_pbt(i,4) - pbt_base(bj)
            bdy_uBar(i) = 0.0_wp
            bdy_vBar(i) = 0.0_wp
         end do

         !$ACC loop independent private(bj)
         do i = 1, nlbdy
            bj = mpi_idx_bdy(i)
            !$ACC loop seq
            do k = 1, nk
               bdy_uBar(i) = bdy_uBar(i) + bdy_u(i, k, 4)*drF(k)*h0FacW(bj,k)*recip_RcolW(bj)
               bdy_vBar(i) = bdy_vBar(i) + bdy_v(i, k, 4)*drF(k)*h0FacS(bj,k)*recip_RcolS(bj)
            end do
         end do

         !$ACC loop independent private(bj)
         do i = 1, nlbdy
            bj = mpi_idx_bdy(i)
            !$ACC loop seq
            do k = 1, nk
               bdy_uBcl(i,k) = bdy_u(i, k, 4) - bdy_uBar(i)*maskW(bj,k)
               bdy_vBcl(i,k) = bdy_v(i, k, 4) - bdy_vBar(i)*maskS(bj,k)
            end do
         end do

         if (ln_tide) then
            !$ACC loop independent private(bj)
            do i = 1, nlbdy
               bj = mpi_idx_bdy(i)
               bdy_uBar(i) = bdy_uBar(i) + tide_u(bj)
               bdy_vBar(i) = bdy_vBar(i) + tide_v(bj)

               bdy_pbt(i,4) = bdy_pbt(i,4) + tide_ssh(bj)*gravityRau0
            end do
         end if

         !$ACC end kernels
      end if

      !$ACC kernels present(uBar,vBar,uFld,vFld,uFldBcl,vFldBcl,uFldBclDash, &
      !$ACC                 vFldBclDash,gUnow,gVnow,drF,maskW,maskS,uBarCouple, &
      !$ACC                 vBarCouple,recip_RcolW,recip_RcolS,dTMom,h0FacW,h0FacS)
      !$ACC loop collapse(2) independent
      do i = 1, nlpb
         do k = 1, nk
            uFldBcl(i,k) = uFld(i,k) - uBar(i)*maskW(i, k)
            vFldBcl(i,k) = vFld(i,k) - vBar(i)*maskS(i, k)
            uFldBclDash(i, k) = uFldBcl(i, k) + dTMom*gUnow(i, k)*maskW(i, k)
            vFldBclDash(i, k) = vFldBcl(i, k) + dTMom*gVnow(i, k)*maskS(i, k)
         end do
      end do

      !$ACC loop independent
      do i = 1, nlpb
         !$ACC loop seq
         do k = 1, nk
            uBarCouple(i) = uBarCouple(i) + uFldBclDash(i, k)*drF(k)*h0FacW(i,k)*recip_RcolW(i)*maskW(i,k)/dTMom
            vBarCouple(i) = vBarCouple(i) + vFldBclDash(i, k)*drF(k)*h0FacS(i,k)*recip_RcolS(i)*maskS(i,k)/dTMom
         end do
      end do

      !$ACC loop collapse(2) independent
      do i = 1, nlpb
         do k = 1, nk
            uFldBcl(i, k) = uFldBclDash(i, k) - dTMom*uBarCouple(i)*maskW(i, k)
            vFldBcl(i, k) = vFldBclDash(i, k) - dTMom*vBarCouple(i)*maskS(i, k)
         end do
      end do
      !$ACC end kernels

      FirstMomU = .false.
      FirstMomV = .false.

   end subroutine csp_dyn_mom_vecinv

!==============================================================================
   subroutine csp_dyn_mom_vecinv_bernoulli(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, w, s, uel, uei, vnl, vni
      real(wp) :: vecFldLoc(nlpb, 2), rAws(nlpb, 2), hFacWS(nlpb, 2)

      !$ACC data present(dKEdx,dKEdy,uFld,vFld,rAw,rAs,hFacW,hFacS,ue,vn,tw,ts,  &
      !$ACC              KE,recip_rAc,recip_hFacC,recip_dxC,recip_dyC,maskW,maskS),  &
      !$ACC      copyin(k),  &
      !$ACC      create(vecFldLoc,rAws,hFacWS)
      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         vecFldLoc(i, iu) = uFld(i, k)
         vecFldLoc(i, iv) = vFld(i, k)
         rAws(i, iu) = rAw(i)
         rAws(i, iv) = rAs(i)
         hFacWS(i, iu) = hFacW(i, k)
         hFacWS(i, iv) = hFacS(i, k)
      end do

      ! KE is define on C grid same with T S
      !$ACC loop independent private(uel,uei,vnl,vni)
      do i = 1, nlpb
         uel = ue(i, 1)
         uei = ue(i, 2)
         vnl = vn(i, 1)
         vni = vn(i, 2)
         KE(i) = 0.25_wp*( &
                 vecFldLoc(i, iu)*vecFldLoc(i, iu)*rAws(i, iu)*hFacWS(i, iu) &
                 + vecFldLoc(uel, uei)*vecFldLoc(uel, uei)*rAws(uel, uei)*hFacWS(uel, uei) &
                 + vecFldLoc(i, iv)*vecFldLoc(i, iv)*rAws(i, iv)*hFacWS(i, iv) &
                 + vecFldLoc(vnl, vni)*vecFldLoc(vnl, vni)*rAws(vnl, vni)*hFacWS(vnl, vni) &
                 )*recip_rAc(i)*recip_hFacC(i, k)
      end do
      !$ACC end kernels

      !$ACC kernels
      !$ACC loop private(w,s)
      do i = 1, nlpb
         w = tw(i)
         dKEdx(i) = -recip_dxC(i)*(KE(i) - KE(w))*maskW(i, k)

         s = ts(i)
         dKEdy(i) = -recip_dyC(i)*(KE(i) - KE(s))*maskS(i, k)
      end do
      !$ACC end kernels
      !$ACC end data

   end subroutine csp_dyn_mom_vecinv_bernoulli

!==============================================================================
   subroutine csp_dyn_mom_vecinv_vshear(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, w, s
      integer :: kp1, km1
      real(wp) :: mask_kp1, mask_km1
      real(wp) :: wBarm, wBarp, Zm, Zp

      kp1 = min(k + 1, nk)
      mask_kp1 = 1.0_wp
      if (k .eq. nk) mask_kp1 = 0.0_wp
      km1 = max(k - 1, 1)
      mask_km1 = 1.0_wp
      if (k .eq. 1) mask_km1 = 0.0_wp

      !$ACC kernels present(tw,ts,uFld,vFld,wFld,rAc,maskC,recip_rAw,recip_rAs,  &
      !$ACC                 recip_hFacW,recip_hFacS,recip_drF,uShearTerm,vShearTerm), &
      !$ACC         copyin(k,kp1,km1,mask_kp1,mask_km1)
      !$ACC loop independent private(w,wBarm,wBarp,Zm,Zp)
      do i = 1, nlpb
         w = tw(i)
         wBarm = 0.5*(wFld(i, k)*rAc(i)*maskC(i, km1) + wFld(w, k)*rAc(w)*maskC(w, km1)) &
                 *mask_km1*recip_rAw(i)

         wBarp = 0.5*(wFld(i, kp1)*rAc(i)*maskC(i, kp1) + wFld(w, kp1)*rAc(w)*maskC(w, kp1)) &
                 *mask_kp1*recip_rAw(i)

         Zm = uFld(i, km1)*mask_km1 - uFld(i, k)

         Zp = uFld(i, k) - uFld(i, kp1)*mask_kp1

         uShearTerm(i) = -0.5*(wBarp*Zp + wBarm*Zm)*recip_hFacW(i, k)*recip_drF(k)
      end do

      !$ACC loop independent private(s,wBarm,wBarp,Zm,Zp)
      do i = 1, nlpb
         s = ts(i)
         wBarm = 0.5*(wFld(i, k)*rAc(i)*maskC(i, km1) + wFld(s, k)*rAc(s)*maskC(s, km1)) &
                 *mask_km1*recip_rAs(i)

         wBarp = 0.5*(wFld(i, kp1)*rAc(i)*maskC(i, kp1) + wFld(s, kp1)*rAc(s)*maskC(s, kp1)) &
                 *mask_kp1*recip_rAs(i)

         Zm = mask_km1*vFld(i, km1) - vFld(i, k)

         Zp = vFld(i, k) - mask_kp1*vFld(i, kp1)

         vShearTerm(i) = -0.5*(wBarp*Zp + wBarm*Zm)*recip_hFacS(i, k)*recip_drF(k)
      end do

      !$ACC end kernels

   end subroutine csp_dyn_mom_vecinv_vshear

!==============================================================================
   subroutine csp_dyn_mom_update_hFacZ(k)
!==============================================================================
      implicit none
      integer :: i, k
      integer :: i1, i2, i3, i4
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      real(wp) :: hFacZOpen, hFacWS(nlpb, 2)

      !$ACC kernels present(nlpbz,hFacZ,recip_hFacZ,hFacW,hFacS,z1,z2,z3,z4),  &
      !$ACC         copyin(k),  &
      !$ACC         create(hFacWS)
      !$ACC loop
      do i = 1, nlpbz
         hFacZ(i) = 0.0_wp
         recip_hFacZ(i) = 0.0_wp
      end do

      !$ACC loop
      do i = 1, nlpb
         hFacWS(i, iu) = hFacW(i, k)
         hFacWS(i, iv) = hFacS(i, k)
      end do

      !$ACC loop private(i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4,hFacZOpen)
      do i = 1, nlpbz
         i1 = z1(i, 1)
         i2 = z2(i, 1)
         i3 = z3(i, 1)
         i4 = z4(i, 1)
         j1 = z1(i, 2)
         j2 = z2(i, 2)
         j3 = z3(i, 2)
         j4 = z4(i, 2)
         k1 = z1(i, 3)
         k2 = z2(i, 3)
         k3 = z3(i, 3)
         k4 = z4(i, 3)
         hFacZOpen = min(hFacWS(i1, j1)*abs(k1), hFacWS(i2, j2)*abs(k2))
         hFacZOpen = min(hFacZOpen, hFacWS(i3, j3)*abs(k3))
         if (k4 .ne. 0) then
            hFacZOpen = min(hFacZOpen, hFacWS(i4, j4)*abs(k4))
         end if
         hFacZ(i) = hFacZOpen
      end do

      !$ACC loop
      do i = 1, nlpbz
         if (hFacZ(i) .ne. 0) then
            recip_hFacZ(i) = 1.0_wp/hFacZ(i)
         else
            recip_hFacZ(i) = 0.0_wp
         end if
      end do
      !$ACC end kernels

   end subroutine csp_dyn_mom_update_hFacZ

end module mod_csp_dyn_mom_vecinv
