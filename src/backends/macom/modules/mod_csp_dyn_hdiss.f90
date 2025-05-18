module mod_csp_dyn_hdiss
   use mod_csp_misc
   use mod_csp_dyn_vort
   implicit none
contains

!==============================================================================
   subroutine csp_dyn_hdiss(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, w, s, e, n
      real(wp) :: Dim, Dij, Dmj, Zip, Zij, Zpj, uD2, vD2, uD4, vD4
      real(wp) :: rp_hFacW, rp_hFacS, maskW_s, maskS_s

      !$ACC data present(uDissip,vDissip,tw,ts,hDiv,viscAh_D,ze,zn,  &
      !$ACC         hFacZ,vort3,viscAh_Z,recip_dxC,recip_dyC,recip_dxZ,  &
      !$ACC         recip_dyZ,recip_hFacW,recip_hFacS,maskW,maskS,  &
      !$ACC         dStar,zStar,maskZ,del2u,del2v,viscA4_D,hFacZ,viscA4_Z),  &
      !$ACC      copyin(k)

      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         uDissip(i) = 0.0_wp
         vDissip(i) = 0.0_wp
      end do
      !$ACC end kernels

      if (harmonic) then
         !$ACC kernels
         !$ACC loop independent private(w,s,e,n,Dij,Dim,Dmj,Zij,Zip,Zpj,uD2,vD2)
         do i = 1, nlpb
            w = tw(i)
            s = ts(i)
            Dij = hDiv(i)*viscAh_D(i)
            Dim = hDiv(s)*viscAh_D(s)
            Dmj = hDiv(w)*viscAh_D(w)

            e = ze(i)
            n = zn(i)
            Zij = hFacZ(i)*vort3(i)*viscAh_Z(i)
            Zip = hFacZ(n)*vort3(n)*viscAh_Z(n)
            Zpj = hFacZ(e)*vort3(e)*viscAh_Z(e)

            uD2 = ((Dij - Dmj)*recip_dxC(i) - (Zip - Zij)*recip_dyZ(i))*recip_hFacW(i, k)
            vD2 = ((Dij - Dim)*recip_dyC(i) + (Zpj - Zij)*recip_dxZ(i))*recip_hFacS(i, k)

            uDissip(i) = uD2*maskW(i, k)
            vDissip(i) = vD2*maskS(i, k)
         end do
         !$ACC end kernels
      end if

      if (biharmonic) then
         call csp_misc_hdiv(del2u, del2v, dStar)
         call csp_dyn_hrelvort(del2u, del2v, zStar)

         !$ACC kernels
         !$ACC loop independent
         do i = 1,nlpbz
            zStar(i) = zStar(i)*maskZ(i, k)
         end do

         !$ACC loop independent private(w,s,e,n,Dij,Dim,Dmj,Zij,Zip,Zpj,uD4,vD4)
         do i = 1, nlpb
            w = tw(i)
            s = ts(i)
            Dim = dStar(s)*viscA4_D(s)
            Dij = dStar(i)*viscA4_D(i)
            Dmj = dStar(w)*viscA4_D(w)

            e = ze(i)
            n = zn(i)
            Zip = hFacZ(n)*zStar(n)*viscA4_Z(n)
            Zij = hFacZ(i)*zStar(i)*viscA4_Z(i)
            Zpj = hFacZ(e)*zStar(e)*viscA4_Z(e)

            uD4 = ((Dij - Dmj)*recip_dxC(i) - (Zip - Zij)*recip_dyZ(i))*recip_hFacW(i, k)
            vD4 = ((Dij - Dim)*recip_dyC(i) + (Zpj - Zij)*recip_dxZ(i))*recip_hFacS(i, k)

            uDissip(i) = (uDissip(i) - uD4)*maskW(i, k)
            vDissip(i) = (vDissip(i) - vD4)*maskS(i, k)
         end do
         !$ACC end kernels
      end if
      !$ACC end data

   end subroutine csp_dyn_hdiss

!==============================================================================
   subroutine csp_dyn_hdiss_coef(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, uel, uei, uep, vnl, vni, vnp
      integer :: i1, i2, i3, i4
      integer :: j1, j2, j3, j4
      integer :: k1, k2, k3, k4, k0
      real(wp) :: zdb, r1_4, csmag, csmag_min, csmag_max, vec_abs
      real(wp) :: temp_max, temp_min, rtemp, dTtemp
      real(wp) :: dtensq(nlpb), dshesq(nlpbz)
      real(wp) :: VecFld(nlpb, 2), dlcz(nlpbz, 2), dlc(nlpb, 2)

      !$ACC data present(uFld,vFld,dxC,dyC,dyZ,dxZ,ue,vn,recip_rAc,maskC,  &
      !$ACC              z1,z2,z3,z4,recip_rAz,maskZ,auz2,auz3,auz5,auz6,  &
      !$ACC              viscAh_D,viscA4_D,viscAh_Z,viscA4_Z,lsmagt2,lsmagf2,dTmom),  &
      !$ACC      copyin(k),  &
      !$ACC      create(dtensq,dshesq,VecFld,dlcz,dlc)

      ! for Smagorinsky viscosity, we should not mask viscAh_Z and viscA4_Z since
      ! side drag need them in lane-sea point
      !$ACC kernels
      !$ACC loop independent
      do i = 1, nlpb
         VecFld(i, iu) = uFld(i, k)
         VecFld(i, iv) = vFld(i, k)
         dlc(i, iu) = dxC(i)
         dlc(i, iv) = dyC(i)
      end do

      !$ACC loop independent
      do i = 1, nlpbz
         dlcz(i, iu) = dyZ(i)
         dlcz(i, iv) = dxZ(i)
      end do

      !$ACC loop independent private(uel,uei,uep,vnl,vni,vnp,zdb)
      do i = 1, nlpb
         uel = ue(i, 1)
         uei = ue(i, 2)
         uep = ue(i, 3)
         vnl = vn(i, 1)
         vni = vn(i, 2)
         vnp = vn(i, 3)
         zdb = (VecFld(uel, uei)*uep*dlcz(uel, uei) - VecFld(i, iu)*dlcz(i, iu) &
                - VecFld(vnl, vni)*vnp*dlcz(vnl, vni) + VecFld(i, iv)*dlcz(i, iv))*recip_rAc(i)
         dtensq(i) = zdb*zdb*maskC(i, k)
      end do

      !$ACC loop independent private(i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4,zdb)
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
         zdb = (+VecFld(i1, j1)*dlc(i1, j1)*k1 + VecFld(i2, j2)*dlc(i2, j2)*k2 &
                - VecFld(i3, j3)*dlc(i3, j3)*k3 - VecFld(i4, j4)*dlc(i4, j4)*k4)*recip_rAz(i)
         dshesq(i) = zdb*zdb*maskZ(i, k)
      end do
      !$ACC end kernels

      ! T-point value
      !$ACC kernels
      !$ACC loop independent private(i1,i2,i3,i4,rtemp,r1_4,vec_abs,temp_max,temp_min)
      do i = 1, nlpb
         csmag = 3.5_wp
         csmag_min = 1.0_wp*csmag
         csmag_max = 1.0_wp*csmag
         i1 = auz2(i)
         i2 = auz3(i)
         i3 = auz5(i)
         i4 = auz6(i)
         rtemp = real(maskZ(i1, k) + maskZ(i2, k) + maskZ(i3, k) + maskZ(i4, k))
         rtemp = max(1.0_wp, rtemp)
         r1_4 = 1.0_wp/rtemp
         viscAh_D(i) = sqrt(dtensq(i) + r1_4*(dshesq(i1) + dshesq(i2) + dshesq(i3) + dshesq(i4)))
         viscAh_D(i) = csmag*viscAh_D(i)*lsmagt2(i)

         vec_abs = sqrt(uFld(i, k)**2 + vFld(i, k)**2)

         temp_max = 0.125_wp*csmag_max*lsmagt2(i)/dTmom
         temp_min = 0.5_wp*csmag_min*sqrt(lsmagt2(i))*vec_abs
         viscAh_D(i) = max(viscAh_D(i), temp_min)
         viscAh_D(i) = min(viscAh_D(i), temp_max)

         viscA4_D(i) = 0.125_wp*viscAh_D(i)*lsmagt2(i)
         temp_max = 0.015625_wp*csmag_max*lsmagt2(i)*lsmagt2(i)/dTmom
         temp_min = 0.083333_wp*csmag_min*sqrt(lsmagt2(i))*lsmagt2(i)*vec_abs
         viscA4_D(i) = max(viscA4_D(i), temp_min)
         viscA4_D(i) = min(viscA4_D(i), temp_max)
      end do

      ! F-point value
      ! this version is not a perfect one, because Z point at corner is only a approximate
      !$ACC loop independent private(i1,i2,i3,i4,j1,j2,j3,j4,rtemp,r1_4,vec_abs,temp_max,temp_min)
      do i = 1, nlpbz
         csmag = 3.5_wp
         csmag_min = 1.0_wp*csmag
         csmag_max = 1.0_wp*csmag
         i1 = z1(i, 1)
         i2 = z2(i, 1)
         i3 = z3(i, 1)
         i4 = z4(i, 1)
         j1 = z1(i, 2)
         j2 = z2(i, 2)
         j3 = z3(i, 2)
         j4 = z4(i, 2)
         rtemp = real(maskC(i1, k) + maskC(i2, k) + maskC(i3, k) + maskC(i4, k))
         rtemp = max(1.0_wp, rtemp)
         r1_4 = 1.0_wp/rtemp
         viscAh_Z(i) = sqrt(dshesq(i) + r1_4*(dtensq(i1) + dtensq(i2) + dtensq(i3) + dtensq(i4)))
         viscAh_Z(i) = csmag*viscAh_Z(i)*lsmagf2(i)

         vec_abs = 2.0_wp*r1_4*sqrt(VecFld(i1, j1)**2 + VecFld(i2, j2)**2 + VecFld(i3, j3)**2 + VecFld(i4, j4)**2)

         temp_max = 0.125_wp*csmag_max*lsmagf2(i)/dTmom
         temp_min = 0.5_wp*csmag_min*sqrt(lsmagf2(i))*vec_abs
         viscAh_Z(i) = max(viscAh_Z(i), temp_min)
         viscAh_Z(i) = min(viscAh_Z(i), temp_max)

         viscA4_Z(i) = 0.125_wp*viscAh_Z(i)*lsmagf2(i)
         temp_max = 0.015625_wp*csmag_max*lsmagf2(i)*lsmagf2(i)/dTmom
         temp_min = 0.083333_wp*csmag_min*sqrt(lsmagf2(i))*lsmagf2(i)*vec_abs
         viscA4_Z(i) = max(viscA4_Z(i), temp_min)
         viscA4_Z(i) = min(viscA4_Z(i), temp_max)
      end do
      !$ACC end kernels
      !$ACC end data

   end subroutine csp_dyn_hdiss_coef

end module mod_csp_dyn_hdiss
