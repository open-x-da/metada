module mod_csp_dyn_drag
   use mod_csp_init
   use mod_csp_misc
   implicit none
contains

!==============================================================================
   subroutine csp_dyn_drag_side(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, n, e, unl, uni, vel, vei
      real(wp) :: hFacZClosedS, hFacZClosedN, hFacZClosedW, hFacZClosedE
      real(wp) :: dxy(nlpb, 2)
      real(wp) :: r_dxy(nlpb, 2)

      !$ACC data present(dxS,dyW,recip_dyW,recip_dxS,zn,un,ze,ve,hFacW,hFacZ,  &
      !$ACC              maskW,maskZ,uDragTerms,recip_hFacW,recip_rAw,         &
      !$ACC              SideDragFac,viscAh_Z,uFld,vFld,viscA4_Z,del2u,del2v,  &
      !$ACC              hFacS,maskS,vDragTerms,recip_hFacS,recip_rAs),  &
      !$ACC      copyin(k),  &
      !$ACC      create(dxy,r_dxy)
      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         dxy(i, iu) = dxS(i)
         dxy(i, iv) = dyW(i)
         r_dxy(i, iu) = recip_dyW(i)
         r_dxy(i, iv) = recip_dxS(i)
      end do

      ! --  u drag at walls
      !$ACC loop independent private(n,unl,uni,hFacZClosedS,hFacZClosedN)
      do i = 1, nlpb
         n = zn(i)
         unl = un(i, 1)
         uni = un(i, 2)
         hFacZClosedS = (hFacW(i, k) - hFacZ(i))*abs(maskW(i, k) - maskZ(i, k))
         hFacZClosedN = (hFacW(i, k) - hFacZ(n))*abs(maskW(i, k) - maskZ(n, k))
         uDragTerms(i) = -recip_hFacW(i, k)*recip_rAw(i)*SideDragFac(i, k) &
                         *(hFacZClosedS*dxy(i, iu)*r_dxy(i, iu) &
                           *(viscAh_Z(i)*uFld(i, k) - viscA4_Z(i)*del2u(i)) &
                           + hFacZClosedN*dxy(unl, uni)*r_dxy(unl, uni) &
                           *(viscAh_Z(n)*uFld(i, k) - viscA4_Z(n)*del2u(i)))
      end do

      ! --  v drag at walls
      !$ACC loop independent private(e,vel,vei,hFacZClosedW,hFacZClosedE)
      do i = 1, nlpb
         e = ze(i)
         vel = ve(i, 1)
         vei = ve(i, 2)
         hFacZClosedW = (hFacS(i, k) - hFacZ(i))*abs(maskS(i, k) - maskZ(i, k))
         hFacZClosedE = (hFacS(i, k) - hFacZ(e))*abs(maskS(i, k) - maskZ(e, k))
         vDragTerms(i) = -recip_hFacS(i, k)*recip_rAs(i)*SideDragFac(i, k) &
                         *(hFacZClosedW*dxy(i, iv)*r_dxy(i, iv) &
                           *(viscAh_Z(i)*vFld(i, k) - viscA4_Z(i)*del2v(i)) &
                           + hFacZClosedE*dxy(vel, vei)*r_dxy(vel, vei) &
                           *(viscAh_Z(e)*vFld(i, k) - viscA4_Z(e)*del2v(i)))
      end do
      !$ACC end kernels
      !$ACC end data

   end subroutine csp_dyn_drag_side

!==============================================================================
   subroutine csp_dyn_drag_bottom
!==============================================================================
      implicit none
      integer :: i, w, s, ibot
      real(wp) :: recDrC, zCdu, zCdv, zm1_2dt
      real(wp) :: DragBottomC(nlpb)

      ! here uDragBottom is reused by model, it set first in mod_csp_vdiff_coef_gls
      !$ACC kernels present(nlpb,nl,nk,kSurfC,uDragBottom,vDragBottom,gravityRau0,   &
      !$ACC                 recip_hFacW,recip_hFacS,recip_drF,tw,ts,uFld,vFld),   &
      !$ACC         create(DragBottomC)
      !$ACC loop independent
      do i = 1, nlpb
         DragBottomC(i) = uDragBottom(i)
      end do

      !$ACC loop independent private(w,s,ibot,recDrC,zCdu,zCdv,zm1_2dt)
      do i = 1, nlpb
         zm1_2dt = -1.0_wp/dTmom

         ibot = min(nk, kSurfC(i))

         ! -- Bottom drag for U
         recDrC = recip_hFacW(i, ibot)*recip_drF(ibot)*gravityRau0
         w = tw(i)
         zCdu = 0.5_wp*(DragBottomC(i) + DragBottomC(w))*recDrC

         uDragBottom(i) = max(zCdu, zm1_2dt)*uFld(i, ibot)

         ! -- Bottom drag for V
         recDrC = recip_hFacS(i, ibot)*recip_drF(ibot)*gravityRau0
         s = ts(i)
         zCdv = 0.5_wp*(DragBottomC(i) + DragBottomC(s))*recDrC

         vDragBottom(i) = max(zCdv, zm1_2dt)*vFld(i, ibot)
      end do
      !$ACC end kernels

   end subroutine csp_dyn_drag_bottom

end module mod_csp_dyn_drag
