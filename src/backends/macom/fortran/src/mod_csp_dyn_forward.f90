module mod_csp_dyn_forward
   use mod_csp_misc
   implicit none
contains

!==============================================================================
   subroutine csp_dyn_forward_guess(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, ibot

      !$ACC data present(gUnow,gVnow,gUpast,gVpast,uCorTerm,dKEdx,uDissip,      &
      !$ACC              uDragTerms,uShearTerm,vCorTerm,dKEdy,vDissip,nlpb,nk,  &
      !$ACC              vDragTerms,vShearTerm,guExt,gvExt,dPhiHydX,dPhiHydY,   &
      !$ACC              kSurfC,maskW,maskS,uDragBottom,vDragBottom),copyin(k)
      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         gUnow(i, k) = gUnow(i, k) + uCorTerm(i) + dKEdx(i) + uDissip(i) + uDragTerms(i) + uShearTerm(i)
         gVnow(i, k) = gVnow(i, k) + vCorTerm(i) + dKEdy(i) + vDissip(i) + vDragTerms(i) + vShearTerm(i)
      end do
      !$ACC end kernels

      if (k .eq. nk) then
         !$ACC kernels
         !$ACC loop independent
         do i = 1, nlpb
            gUnow(i, nk) = gUnow(i, nk) + guExt(i)
            gVnow(i, nk) = gVnow(i, nk) + gvExt(i)
         end do
         !$ACC end kernels
      end if

      if (k .eq. 1) then
         !$ACC kernels
         !$ACC loop independent private(ibot)
         do i = 1, nlpb
            ibot = min(nk, kSurfC(i))
            gUnow(i, ibot) = gUnow(i, ibot) + uDragBottom(i)
            gVnow(i, ibot) = gVnow(i, ibot) + vDragBottom(i)
         end do
         !$ACC end kernels
      end if

      call csp_misc_adams_bashforth2(gUnow(1:nlpb, k), gUpast(1:nlpb, k), FirstMomU)

      call csp_misc_adams_bashforth2(gVnow(1:nlpb, k), gVpast(1:nlpb, k), FirstMomV)

      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         gUnow(i, k) = (gUnow(i, k) - dPhiHydX(i))*maskW(i, k)
         gVnow(i, k) = (gVnow(i, k) - dPhiHydY(i))*maskS(i, k)
      end do
      !$ACC end kernels
      !$ACC end data

   end subroutine csp_dyn_forward_guess

!==============================================================================
   subroutine csp_dyn_forward_correct
!==============================================================================
      implicit none

      call csp_dyn_forward_update_vertical_velocity

      call csp_dyn_forward_correct_p_star

   end subroutine csp_dyn_forward_correct

!==============================================================================
   subroutine csp_dyn_forward_correct_velocity
!==============================================================================
      implicit none
      integer :: i, k

      !$ACC kernels present(uFld, vFld, uFldBcl, vFldBcl, uBar, vBar, maskW, maskS)
      !$ACC loop collapse(2) independent
      do k = 1, nk
         do i = 1, nlpb
            uFld(i, k) = (uFldBcl(i, k) + uBar(i))*maskW(i, k)
            vFld(i, k) = (vFldBcl(i, k) + vBar(i))*maskS(i, k)
         end do
      end do
      !$ACC end kernels

      CALL mpi_data_exchange(mpi_sendbuf_2d_nk_f, uFld, nk)
      CALL mpi_data_exchange(mpi_sendbuf_2d_nk_f, vFld, nk)

   end subroutine csp_dyn_forward_correct_velocity

!==============================================================================
   subroutine csp_dyn_forward_update_vertical_velocity
!==============================================================================
      implicit none
      integer :: i, k, uel, uei, uep, vnl, vni, vnp
      real(wp) :: hDivFlow(nlpb), vecTrans(nlpb, 2, nk), conv2d(nlpb)
      real(wp) :: rStarDhDt(nlpb), dEtaHdt(nlpb)
      real(dp) :: sshtemp, sshgloavg

      !$ACC data create(hDivFlow, vecTrans, conv2d, rStarDhDt, dEtaHdt),   &
      !$ACC     present(nlpb, nkp1, nk, gravityDynForce, gravityRau0, uFld,  &
      !$ACC             vFld, wFld, dxZ, dyZ, drF, hFacW, hFacS, maskC,  &
      !$ACC             addMass, EmPmR, recip_RcolC, recip_rAc, ue, vn,  &
      !$ACC             ssh, etaH, h0FacC, mpi_sendbuf_2d_nkp1)
      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         hDivFlow(i) = 0.0_wp
         conv2d(i) = 0.0_wp
      end do

      !$ACC loop collapse(2) independent
      do k = 1, nkp1
         do i = 1, nlpb
            wFld(i, k) = 0.0_wp
         end do
      end do

      !$ACC loop collapse(2) independent
      do k = nk, 1, -1
         do i = 1, nlpb
            vecTrans(i, iu, k) = uFld(i, k)*dyZ(i)*drF(k)*hFacW(i, k)
            vecTrans(i, iv, k) = vFld(i, k)*dxZ(i)*drF(k)*hFacS(i, k)
         end do
      end do
      !$ACC end kernels

      !$ACC kernels
      !$ACC loop independent private(uel, uei, uep, vnl, vni, vnp)
      do i = 1, nlpb
         !$ACC loop seq
         do k = nk, 1, -1
            uel = ue(i, 1)
            uei = ue(i, 2)
            uep = ue(i, 3)
            vnl = vn(i, 1)
            vni = vn(i, 2)
            vnp = vn(i, 3)
            conv2d(i) = -maskC(i, k) &
                        *(vecTrans(uel, uei, k)*uep - vecTrans(i, iu, k)  &
                        + vecTrans(vnl, vni, k)*vnp - vecTrans(i, iv, k)) &
                        + gravityDynForce*addMass(i, k)

            hDivFlow(i) = hDivFlow(i) - conv2d(i)
         end do
      end do
      !$ACC end kernels

      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         dEtaHdt(i) = -hDivFlow(i)*recip_rAc(i) - gravityDynForce*EmPmR(i)

         rStarDhDt(i) = dEtaHdt(i)*recip_RcolC(i)

         conv2d(i) = 0.0_wp

         ssh(i) = ssh(i) + etaH(i)/gravityRau0
      end do
      !$ACC end kernels

      !$ACC kernels
      !$ACC loop independent private(uel, uei, uep, vnl, vni, vnp)
      do i = 1, nlpb
         !$ACC loop seq
         do k = nk, 2, -1
            uel = ue(i, 1)
            uei = ue(i, 2)
            uep = ue(i, 3)
            vnl = vn(i, 1)
            vni = vn(i, 2)
            vnp = vn(i, 3)
            conv2d(i) = -maskC(i, k) &
                        *(vecTrans(uel, uei, k)*uep - vecTrans(i, iu, k)  &
                        + vecTrans(vnl, vni, k)*vnp - vecTrans(i, iv, k)) &
                        + gravityDynForce*addMass(i, k)

            if (k .eq. nk) then
               wFld(i, k) = conv2d(i)*recip_rAc(i) - rStarDhDt(i)*drF(k)*h0FacC(i, k)
               wFld(i, k) = wFld(i, k) - gravityDynForce*EmPmR(i)*maskC(i, k)
            else
               wFld(i, k) = wFld(i, k + 1) + conv2d(i)*recip_rAc(i) &
                            - rStarDhDt(i)*drF(k)*h0FacC(i, k)
            end if

            wFld(i, k) = wFld(i, k) * maskC(i, k) * maskC(i, k-1)
         end do
      end do
      !$ACC end kernels

      CALL mpi_data_exchange(mpi_sendbuf_2d_nkp1, wFld, nkp1)
      !$ACC end data

   end subroutine csp_dyn_forward_update_vertical_velocity

!==============================================================================
   subroutine csp_dyn_forward_correct_p_star
!==============================================================================
      implicit none
      integer :: i, k, w, s
      real(wp) :: tmpfldW, tmpfldS

      !$ACC data present(nlpb, loc_nlpb, nk, rStarExpC, rStarFacC, rStarFacW,  &
      !$ACC             rStarFacS, etaH, RcolC, RcolW, RcolS, tw, ts, rAc, &
      !$ACC             recip_rAw, recip_rAs, recip_RcolC, maskC, maskW, maskS, &
      !$ACC             latC, lonC, rStarFacLow, rStarFacUp, hFacC, hFacW, &
      !$ACC             hFacS, h0FacC, h0FacW, h0FacS, recip_hFacC, recip_hFacS, &
      !$ACC             recip_hFacW, pbt_base)

      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         rStarExpC(i) = rStarFacC(i)

         rStarFacC(i) = (etaH(i) + RcolC(i) + pbt_base(i))*recip_RcolC(i)

         rStarExpC(i) = rStarFacC(i)/rStarExpC(i)
      end do

      !$ACC loop independent private(w, tmpfldW)
      do i = 1, nlpb
         w = tw(i)
         tmpfldW = RcolW(i)
         !tmpfldW = 0.5_wp*(RcolC(w)*rAc(w) + RcolC(i)*rAc(i))*recip_rAw(i)
         rStarFacW(i) = (0.5_wp*((etaH(w)+pbt_base(w))*rAc(w) + (etaH(i)+pbt_base(i))*rAc(i)) &
                         *recip_rAw(i) + tmpfldW)/tmpfldW
      end do

      !$ACC loop independent private(s, tmpfldS)
      do i = 1, nlpb
         s = ts(i)
         tmpfldS = RcolS(i)
         !tmpfldS = 0.5_wp*(RcolC(s)*rAc(s) + RcolC(i)*rAc(i))*recip_rAs(i)
         rStarFacS(i) = (0.5_wp*((etaH(s)+pbt_base(s))*rAc(s) + (etaH(i)+pbt_base(i))*rAc(i)) &
                         *recip_rAs(i) + tmpfldS)/tmpfldS
      end do

      !$ACC loop independent
      do i = 1, nlpb
         if (maskC(i, nk) .eq. 0) rStarFacC(i) = 1.0_wp
         if (maskW(i, nk) .eq. 0) rStarFacW(i) = 1.0_wp
         if (maskS(i, nk) .eq. 0) rStarFacS(i) = 1.0_wp
      end do

      !$ACC loop independent
      do i = 1, loc_nlpb
         if (rStarFacC(i) .lt. 0.0_wp) then
            write (*, *) "Error: rStarFacC = ", rStarFacC(i), " at i = ", i, latC(i), lonC(i)
            stop
         end if
         if (rStarFacW(i) .lt. 0.0_wp) then
            write (*, *) "Error: rStarFacW = ", rStarFacW(i), " at i = ", i, latC(i), lonC(i)
            stop
         end if
         if (rStarFacS(i) .lt. 0.0_wp) then
            write (*, *) "Error: rStarFacS = ", rStarFacS(i), " at i = ", i, latC(i), lonC(i)
            stop
         end if
      end do
      !$ACC end kernels

      !$ACC kernels
      !$ACC loop collapse(2) independent
      do k = 1, nk
         do i = 1, nlpb
            hFacC(i, k) = h0FacC(i, k)*rStarFacC(i)
            hFacW(i, k) = h0FacW(i, k)*rStarFacW(i)
            hFacS(i, k) = h0FacS(i, k)*rStarFacS(i)
         end do
      end do

      !$ACC loop collapse(2) independent
      do k = 1, nk
         do i = 1, nlpb
            if (hFacC(i, k) .ne. 0.0_wp) then
               recip_hFacC(i, k) = 1.0_wp/hFacC(i, k)
            else
               recip_hFacC(i, k) = 0.0_wp
            end if

            if (hFacW(i, k) .ne. 0.0_wp) then
               recip_hFacW(i, k) = 1.0_wp/hFacW(i, k)
            else
               recip_hFacW(i, k) = 0.0_wp
            end if

            if (hFacS(i, k) .ne. 0.0_wp) then
               recip_hFacS(i, k) = 1.0_wp/hFacS(i, k)
            else
               recip_hFacS(i, k) = 0.0_wp
            end if
         end do
      end do
      !$ACC end kernels

      !$ACC end data

   end subroutine csp_dyn_forward_correct_p_star

end module mod_csp_dyn_forward
