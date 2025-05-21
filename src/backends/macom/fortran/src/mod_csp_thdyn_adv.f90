module mod_csp_thdyn_adv
   use mod_csp_init
   implicit none

contains

!==============================================================================
   subroutine csp_thdyn_adv(its)
!==============================================================================
      implicit none
      integer, intent(in) :: its

      call csp_thdyn_adv_dst3fl(its)

      ! call csp_thdyn_adv_ndsl

   end subroutine csp_thdyn_adv

!==============================================================================
   subroutine csp_thdyn_adv_ndsl
!==============================================================================
      implicit none

      call csp_thdyn_adv_ppm

   end subroutine csp_thdyn_adv_ndsl

!==============================================================================
   subroutine csp_thdyn_adv_ppm
!==============================================================================
      implicit none

   end subroutine csp_thdyn_adv_ppm

!==============================================================================
   subroutine csp_thdyn_adv_dst3fl(its)
!==============================================================================
      implicit none
      integer, intent(in) :: its
      integer :: i, k, uel, uei, uep, vnl, vni, vnp, i_max
      integer :: e, w, n, s, w2, s2, km1, km2, kp1
      real(wp) :: Rjp, Rj, Rjm
      real(wp) :: uCFL, vCFL, wCFL, d0, d1
      real(wp) :: thetaP, thetaM, psiP, psiM
      real(wp) :: af(nlpb, 2)
      real(wp) :: xyA_iu, xyA_iv
      real(wp) :: Trans(nlpb, 2)
      real(wp) :: tl(nlpb, nk)
      real(wp) :: TransKp1(nlpb), TransK(nlpb), afK(nlpb), afKp1(nlpb)

      !$ACC data present(nk,nlpb,gTracer,tracer,uFld,vFld,dyZ,dxZ,drF,hFacW,hFacS,  &
      !$ACC      te,tw,tw2,ue,maskW,maskS,dTtracer,recip_dxC,recip_dyC,avoid0,  &
      !$ACC      tn,ts,ts2,vn,recip_hFacC,recip_drF,recip_rAc,maskC,rAc,wFld,  &
      !$ACC      mpi_sendbuf_1d, tFld_adv, sFld_adv),  &
      !$ACC      create(tl,Trans,af,TransKp1,TransK,afK,afKp1)

      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            tl(i, k) = tracer(i, k)
         end do
      end do

      !$ACC loop
      do i = 1, nlpb
         TransKp1(i) = 0.0_wp
         TransK(i) = 0.0_wp
         afK(i) = 0.0_wp
         afKp1(i) = 0.0_wp
      end do
      !$ACC end kernels

      do k = nk, 2, -1
         ! Calculate "mass transports" through tracer cell faces.
         !$ACC kernels
         !$ACC loop
         do i = 1, nlpb
            xyA_iu = dyZ(i)*drF(k)*hFacW(i, k)
            xyA_iv = dxZ(i)*drF(k)*hFacS(i, k)
            Trans(i, iu) = uFld(i, k)*xyA_iu
            Trans(i, iv) = vFld(i, k)*xyA_iv

            ! advection at u direction
            af(i, iu) = 0.0_wp
            af(i, iv) = 0.0_wp
         end do

         !$ACC loop private(e,w,w2,Rjp,Rj,Rjm,uCFL,d0,d1,thetaP,thetaM,psiP,psiM,n,s,s2,vCFL)
         do i = 1, nlpb
            e = te(i)
            w = tw(i)
            w2 = tw2(i)

            Rjp = (tl(e, k) - tl(i, k))*maskW(e, k)
            Rj = (tl(i, k) - tl(w, k))*maskW(i, k)
            Rjm = (tl(w, k) - tl(w2, k))*maskW(w, k)

            uCFL = abs(uFld(i, k)*dTtracer*recip_dxC(i))

            d0 = (2.0_wp - uCFL)*(1.0_wp - uCFL)/6.0_wp
            d1 = (1.0_wp - uCFL*uCFL)/6.0_wp

            if (abs(Rj)*thetaMax .le. abs(Rjm)) then
               thetaP = sign(thetaMax, Rjm*Rj)
            else
               thetaP = Rjm/Rj
            end if

            if (abs(Rj)*thetaMax .le. abs(Rjp)) then
               thetaM = sign(thetaMax, Rjp*Rj)
            else
               thetaM = Rjp/Rj
            end if

            psiP = d0 + d1*thetaP
            psiP = max(0.0_wp, min(min(1.0_wp, psiP), thetaP*(1.0_wp - uCFL)/(uCFL + avoid0)))
            psiM = d0 + d1*thetaM
            psiM = max(0.0_wp, min(min(1.0_wp, psiM), thetaM*(1.0_wp - uCFL)/(uCFL + avoid0)))

            af(i, iu) = 0.5_wp*(Trans(i, iu) + abs(Trans(i, iu)))*(tl(w, k) + psiP*Rj) &
                        + 0.5_wp*(Trans(i, iu) - abs(Trans(i, iu)))*(tl(i, k) - psiM*Rj)

            n = tn(i)
            s = ts(i)
            s2 = ts2(i)

            Rjp = (tl(n, k) - tl(i, k))*maskS(n, k)
            Rj = (tl(i, k) - tl(s, k))*maskS(i, k)
            Rjm = (tl(s, k) - tl(s2, k))*maskS(s, k)

            vCFL = abs(vFld(i, k)*dTtracer*recip_dyC(i))
            d0 = (2.0_wp - vCFL)*(1.0_wp - vCFL)/6.0_wp
            d1 = (1.0_wp - vCFL*vCFL)/6.0_wp

            if (abs(Rj)*thetaMax .le. abs(Rjm)) then
               thetaP = sign(thetaMax, Rjm*Rj)
            else
               thetaP = Rjm/Rj
            end if

            if (abs(Rj)*thetaMax .LE. abs(Rjp)) then
               thetaM = sign(thetaMax, Rjp*Rj)
            else
               thetaM = Rjp/Rj
            end if

            psiP = d0 + d1*thetaP
            psiP = max(0.0_wp, min(min(1.0_wp, psiP), thetaP*(1.0_wp - vCFL)/(vCFL + avoid0)))
            psiM = d0 + d1*thetaM
            psiM = max(0.0_wp, min(min(1.0_wp, psiM), thetaM*(1.0_wp - vCFL)/(vCFL + avoid0)))

            af(i, iv) = 0.5*(Trans(i, iv) + abs(Trans(i, iv)))*(tl(s, k) + psiP*Rj) &
                        + 0.5*(Trans(i, iv) - abs(Trans(i, iv)))*(tl(i, k) - psiM*Rj)
         end do

         !$ACC loop private(uel,uei,uep)
         do i = 1, nlpb
            uel = ue(i, 1)
            uei = ue(i, 2)
            uep = ue(i, 3)
            tl(i, k) = tl(i, k) &
                       - dTtracer*recip_hFacC(i, k)*recip_drF(k)*recip_rAc(i)*maskC(i, k) &
                       *(af(uel, uei)*uep - af(i, iu) - tracer(i, k)*(Trans(uel, uei)*uep - Trans(i, iu)))

         end do
         !$ACC end kernels

         CALL mpi_data_exchange(mpi_sendbuf_1d, tl(:, k))

         ! advection at y direction
         !$ACC kernels
         !$ACC loop private(e,w,w2,Rjp,Rj,Rjm,uCFL,d0,d1,thetaP,thetaM,psiP,psiM,n,s,s2,vCFL)
         do i = 1, nlpb
            af(i, iu) = 0.0_wp
            af(i, iv) = 0.0_wp

            e = te(i)
            w = tw(i)
            w2 = tw2(i)

            Rjp = (tl(e, k) - tl(i, k))*maskW(e, k)
            Rj = (tl(i, k) - tl(w, k))*maskW(i, k)
            Rjm = (tl(w, k) - tl(w2, k))*maskW(w, k)

            uCFL = abs(uFld(i, k)*dTtracer*recip_dxC(i))

            d0 = (2.0_wp - uCFL)*(1.0_wp - uCFL)/6.0_wp
            d1 = (1.0_wp - uCFL*uCFL)/6.0_wp

            if (abs(Rj)*thetaMax .le. abs(Rjm)) then
               thetaP = sign(thetaMax, Rjm*Rj)
            else
               thetaP = Rjm/Rj
            end if

            if (abs(Rj)*thetaMax .le. abs(Rjp)) then
               thetaM = sign(thetaMax, Rjp*Rj)
            else
               thetaM = Rjp/Rj
            end if

            psiP = d0 + d1*thetaP
            psiP = max(0.0_wp, min(min(1.0_wp, psiP), thetaP*(1.0_wp - uCFL)/(uCFL + avoid0)))
            psiM = d0 + d1*thetaM
            psiM = max(0.0_wp, min(min(1.0_wp, psiM), thetaM*(1.0_wp - uCFL)/(uCFL + avoid0)))

            af(i, iu) = 0.5_wp*(Trans(i, iu) + abs(Trans(i, iu)))*(tl(w, k) + psiP*Rj) &
                        + 0.5_wp*(Trans(i, iu) - abs(Trans(i, iu)))*(tl(i, k) - psiM*Rj)

            n = tn(i)
            s = ts(i)
            s2 = ts2(i)

            Rjp = (tl(n, k) - tl(i, k))*maskS(n, k)
            Rj = (tl(i, k) - tl(s, k))*maskS(i, k)
            Rjm = (tl(s, k) - tl(s2, k))*maskS(s, k)

            vCFL = abs(vFld(i, k)*dTtracer*recip_dyC(i))
            d0 = (2.0_wp - vCFL)*(1.0_wp - vCFL)/6.0_wp
            d1 = (1.0_wp - vCFL*vCFL)/6.0_wp

            if (abs(Rj)*thetaMax .le. abs(Rjm)) then
               thetaP = sign(thetaMax, Rjm*Rj)
            else
               thetaP = Rjm/Rj
            end if

            if (abs(Rj)*thetaMax .LE. abs(Rjp)) then
               thetaM = sign(thetaMax, Rjp*Rj)
            else
               thetaM = Rjp/Rj
            end if

            psiP = d0 + d1*thetaP
            psiP = max(0.0_wp, min(min(1.0_wp, psiP), thetaP*(1.0_wp - vCFL)/(vCFL + avoid0)))
            psiM = d0 + d1*thetaM
            psiM = max(0.0_wp, min(min(1.0_wp, psiM), thetaM*(1.0_wp - vCFL)/(vCFL + avoid0)))

            af(i, iv) = 0.5_wp*(Trans(i, iv) + abs(Trans(i, iv)))*(tl(s, k) + psiP*Rj) &
                      + 0.5_wp*(Trans(i, iv) - abs(Trans(i, iv)))*(tl(i, k) - psiM*Rj)
         end do

         !$ACC loop private(vnl,vni,vnp)
         do i = 1, nlpb
            vnl = vn(i, 1)
            vni = vn(i, 2)
            vnp = vn(i, 3)
            tl(i, k) = tl(i, k) &
                       - dTtracer*recip_hFacC(i, k)*recip_drF(k)*recip_rAc(i)*maskC(i, k) &
                       *(af(vnl, vni)*vnp - af(i, iv) - tracer(i, k)*(Trans(vnl, vni)*vnp - Trans(i, iv)))

         end do
         !$ACC end kernels

         !$ACC kernels
         !$ACC loop private(km1,km2,kp1,Rjp,Rj,Rjm,wCFL,d0,d1,thetaP,thetaM,psiP,psiM)
         do i = 1, nlpb
            km2 = max( 1,k-2)
            km1 = max( 1,k-1)
            kp1 = min(nk,k+1)

            TransKp1(i) = TransK(i)
            TransK(i) = wFld(i, k)*rAc(i)*maskC(i, k-1)*maskC(i, k)

            Rjp = (tl(i,k) - tl(i,kp1))*maskC(i,kp1)
            Rj  = (tl(i,km1) - tl(i,k))*maskC(i,k)*maskC(i,km1)
            Rjm = (tl(i,km2) - tl(i,km1))*maskC(i,km1)

            wCFL = abs(wFld(i, k)*dTtracer*recip_drC(k))
            d0 = (2.0_wp -wCFL)*(1.0_wp -wCFL)/6.0_wp
            d1 = (1.0_wp -wCFL*wCFL)/6.0_wp

            if ( abs(Rj)*thetaMax .le. abs(Rjm) ) then
               thetaP = sign(thetaMax,Rjm*Rj)
            else
               thetaP = Rjm/Rj
            end if
            if ( abs(Rj)*thetaMax .le. abs(Rjp) ) then
               thetaM = sign(thetaMax,Rjp*Rj)
            else
               thetaM = Rjp/Rj
            end if

            psiP = d0+d1*thetaP
            psiP = max(0.0_wp,min(min(1.0_wp,psiP),thetaP*(1.0_wp -wCFL)/(wCFL+1.0e-20_wp) ))
            psiM = d0+d1*thetaM
            psiM = max(0.0_wp,min(min(1.0_wp,psiM),thetaM*(1.0_wp -wCFL)/(wCFL+1.0e-20_wp) ))

            afKp1(i) = afK(i)
            afK(i) = 0.5_wp*(TransK(i)+abs(TransK(i)))*( tl(i,  k) + psiM*Rj ) &
                   + 0.5_wp*(TransK(i)-abs(TransK(i)))*( tl(i,km1) - psiP*Rj )

            tl(i, k) = tl(i, k) &
                     - dTtracer*recip_hFacC(i, k)*recip_drF(k)*recip_rAc(i)*maskC(i, k) &
                     *(afK(i) - afKp1(i) - tracer(i, k)*(TransK(i) - TransKp1(i)))

         end do
         !$ACC end kernels
      end do

      ! calculate tracer tendency by advection
      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            gTracer(i, k) = gTracer(i, k) + (tl(i, k) - tracer(i, k))/dTtracer
         end do
      end do
      !$ACC end kernels

# if defined (EXDIAG)
      if (its .eq. 1) then
         !$ACC kernels
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, loc_nlpb
               tFld_adv(i, k) = tFld_adv(i, k) + (tl(i, k) - tracer(i, k))
            end do
         end do
         !$ACC end kernels
      end if

      if (its .eq. 2) then
         !$ACC kernels
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, loc_nlpb
               sFld_adv(i, k) = sFld_adv(i, k) + (tl(i, k) - tracer(i, k))
            end do
         end do
         !$ACC end kernels
      end if
# endif

      !$ACC end data

   end subroutine csp_thdyn_adv_dst3fl

end module mod_csp_thdyn_adv
