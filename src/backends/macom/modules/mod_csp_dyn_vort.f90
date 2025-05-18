module mod_csp_dyn_vort
   use mod_csp_init
   implicit none

contains

!==============================================================================
   subroutine csp_dyn_hrelvort(uFldLoc, vFldLoc, Vort3Loc)
!==============================================================================
      implicit none
      integer :: i
      integer :: i1, i2, i3, i4
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      real(wp), intent(in) :: uFldLoc(nlpb), vFldLoc(nlpb)
      real(wp), intent(out) :: Vort3Loc(nlpbz)
      real(wp) :: VecFld(nlpb, 2), dlc(nlpb, 2)

      !$ACC kernels present(uFldLoc,vFldLoc,Vort3Loc,dxC,dyC,z1,z2,z3,z4,recip_rAz),   &
      !$ACC         create(VecFld,dlc)
      !$ACC loop independent
      do i = 1, nlpb
         VecFld(i, iu) = uFldLoc(i)
         VecFld(i, iv) = vFldLoc(i)
         dlc(i, iu) = dxC(i)
         dlc(i, iv) = dyC(i)
      end do

      !$ACC loop independent private(i1, i2, i3, i4,j1, j2, j3, j4,k1, k2, k3, k4)
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
         Vort3Loc(i) = (-VecFld(i1, j1)*dlc(i1, j1)*k1 + VecFld(i2, j2)*dlc(i2, j2)*k2 &
                        + VecFld(i3, j3)*dlc(i3, j3)*k3 - VecFld(i4, j4)*dlc(i4, j4)*k4)*recip_rAz(i)

      end do

      !$ACC end kernels

   end subroutine csp_dyn_hrelvort

!==============================================================================
   subroutine csp_dyn_vort_een(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i
      integer :: i1, i2, i3, i4, i5, i6
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4
      integer :: p1, p2, p3, p4
      real(wp) :: oneThird, oneFour
      real(wp) :: vort3uw, vort3ue, vort3dw, vort3de
      real(wp) :: VecFld(nlpb,2), dlz(nlpbz,2), hFacWS(nlpb,2)
      real(wp) :: maskW_s, maskS_s
      integer(4) :: maskWS(nlpb,2)

      !$ACC kernels present(uFld,vFld,hFacW,hFacS,dyZ,dxZ,  &
      !$ACC         auz1,auz2,auz3,auz4,auz5,auz6,au1,au2,au3,au4,  &
      !$ACC         recip_hFacZ,omega3,uCorTerm,recip_dxC,maskW,  &
      !$ACC         avz1,avz2,avz3,avz4,avz5,avz6,av1,av2,av3,av4,  &
      !$ACC         vCorTerm,recip_dyC,maskS,maskZ),  &
      !$ACC         copyin(k),  &
      !$ACC         create(VecFld,dlz,hFacWS,maskWS)

      !$ACC loop independent
      do i = 1, nlpb
         VecFld(i, iu) = uFld(i, k)
         VecFld(i, iv) = vFld(i, k)
         hFacWS(i, iu) = hFacW(i, k)
         hFacWS(i, iv) = hFacS(i, k)
         maskWS(i, iu) = maskW(i, k)
         maskWS(i, iv) = maskS(i, k)
      end do

      !$ACC loop independent
      do i = 1, nlpbz
         dlz(i, iu) = dyZ(i)
         dlz(i, iv) = dxZ(i)
      end do

      !$ACC loop independent private(i1, i2, i3, i4, i5, i6, j1, j2, j3, j4,  &
      !$ACC                  k1, k2, k3, k4, p1, p2, p3, p4, oneThird, oneFour, &
      !$ACC                  vort3uw, vort3dw, vort3ue, vort3de)
      do i = 1, nlpb
         i1 = auz1(i)
         i2 = auz2(i)
         i3 = auz3(i)
         i4 = auz4(i)
         i5 = auz5(i)
         i6 = auz6(i)
         j1 = au1(i, 1)
         j2 = au2(i, 1)
         j3 = au3(i, 1)
         j4 = au4(i, 1)
         k1 = au1(i, 2)
         k2 = au2(i, 2)
         k3 = au3(i, 2)
         k4 = au4(i, 2)
         p1 = au1(i, 3)
         p2 = au2(i, 3)
         p3 = au3(i, 3)
         p4 = au4(i, 3)

         oneThird = maskZ(i1,k)+maskZ(i2,k)+maskZ(i5,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3uw = (recip_hFacZ(i1)*omega3(i1) &
                    + recip_hFacZ(i2)*omega3(i2) &
                    + recip_hFacZ(i5)*omega3(i5)) &
                   *oneThird*VecFld(j1, k1)*p1*dlz(j1, k1)*hFacWS(j1, k1)

         oneThird = maskZ(i3,k)+maskZ(i2,k)+maskZ(i5,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3ue = (recip_hFacZ(i3)*omega3(i3) &
                    + recip_hFacZ(i2)*omega3(i2) &
                    + recip_hFacZ(i5)*omega3(i5)) &
                   *oneThird*VecFld(j2, k2)*p2*dlz(j2, k2)*hFacWS(j2, k2)

         oneThird = maskZ(i4,k)+maskZ(i5,k)+maskZ(i2,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3dw = (recip_hFacZ(i4)*omega3(i4) &
                    + recip_hFacZ(i5)*omega3(i5) &
                    + recip_hFacZ(i2)*omega3(i2)) &
                   *oneThird*VecFld(j3, k3)*p3*dlz(j3, k3)*hFacWS(j3, k3)

         oneThird = maskZ(i6,k)+maskZ(i5,k)+maskZ(i2,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3de = (recip_hFacZ(i6)*omega3(i6) &
                    + recip_hFacZ(i5)*omega3(i5) &
                    + recip_hFacZ(i2)*omega3(i2)) &
                   *oneThird*VecFld(j4, k4)*p4*dlz(j4, k4)*hFacWS(j4, k4)

         oneFour = maskWS(j1,k1) + maskWS(j2,k2) + maskWS(j3,k3) + maskWS(j4,k4)
         if (oneFour .eq. 0.0_wp) then
            oneFour = 0.0_wp
         else
            oneFour = 1.0_wp/oneFour
         end if
         uCorTerm(i) = +(vort3uw + vort3ue + vort3dw + vort3de)*oneFour*recip_dxC(i)*maskW(i, k)
      end do

      !$ACC loop independent private(i1, i2, i3, i4, i5, i6, j1, j2, j3, j4,  &
      !$ACC                  k1, k2, k3, k4, p1, p2, p3, p4, oneThird, oneFour, &
      !$ACC                  vort3uw, vort3dw, vort3ue, vort3de)
      do i = 1, nlpb
         i1 = avz1(i)
         i2 = avz2(i)
         i3 = avz3(i)
         i4 = avz4(i)
         i5 = avz5(i)
         i6 = avz6(i)
         j1 = av1(i, 1)
         j2 = av2(i, 1)
         j3 = av3(i, 1)
         j4 = av4(i, 1)
         k1 = av1(i, 2)
         k2 = av2(i, 2)
         k3 = av3(i, 2)
         k4 = av4(i, 2)
         p1 = av1(i, 3)
         p2 = av2(i, 3)
         p3 = av3(i, 3)
         p4 = av4(i, 3)

         oneThird = maskZ(i1,k)+maskZ(i2,k)+maskZ(i5,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3uw = (recip_hFacZ(i1)*omega3(i1) &
                    + recip_hFacZ(i2)*omega3(i2) &
                    + recip_hFacZ(i5)*omega3(i5)) &
                   *oneThird*VecFld(j1, k1)*p1*dlz(j1, k1)*hFacWS(j1, k1)

         oneThird = maskZ(i3,k)+maskZ(i2,k)+maskZ(i5,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3dw = (recip_hFacZ(i3)*omega3(i3) &
                    + recip_hFacZ(i2)*omega3(i2) &
                    + recip_hFacZ(i5)*omega3(i5)) &
                   *oneThird*VecFld(j2, k2)*p2*dlz(j2, k2)*hFacWS(j2, k2)

         oneThird = maskZ(i4,k)+maskZ(i5,k)+maskZ(i2,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3ue = (recip_hFacZ(i4)*omega3(i4) &
                    + recip_hFacZ(i5)*omega3(i5) &
                    + recip_hFacZ(i2)*omega3(i2)) &
                   *oneThird*VecFld(j3, k3)*p3*dlz(j3, k3)*hFacWS(j3, k3)

         oneThird = maskZ(i6,k)+maskZ(i5,k)+maskZ(i2,k)
         if (oneThird .eq. 0.0_wp) then
            oneThird = 0.0_wp
         else
            oneThird = 1.0_wp/oneThird
         end if
         vort3de = (recip_hFacZ(i6)*omega3(i6) &
                    + recip_hFacZ(i5)*omega3(i5) &
                    + recip_hFacZ(i2)*omega3(i2)) &
                   *oneThird*VecFld(j4, k4)*p4*dlz(j4, k4)*hFacWS(j4, k4)
         
         oneFour = maskWS(j1,k1) + maskWS(j2,k2) + maskWS(j3,k3) + maskWS(j4,k4)
         if (oneFour .eq. 0.0_wp) then
            oneFour = 0.0_wp
         else
            oneFour = 1.0_wp/oneFour
         end if
         vCorTerm(i) = -(vort3uw + vort3dw + vort3ue + vort3de)*oneFour*recip_dyC(i)*maskS(i, k)
      end do

      !$ACC end kernels

   end subroutine csp_dyn_vort_een

end module mod_csp_dyn_vort
