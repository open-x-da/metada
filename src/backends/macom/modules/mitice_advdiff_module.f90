module mitice_advdiff
use mitice_parameters
use mitice_vars 
use mitice_itd
use mod_misc_basic
use mod_csp_basic
use mod_mpi_interfaces

   implicit none
!  Trans    :: volume transport, x/y direction
!  xyA     :: "areas" of X and Y face of tracer cells
   real(wp), dimension(:,:), allocatable   ::  Trans
   real(wp), dimension(:,:), allocatable   ::  xyA
   real(wp), dimension(:), allocatable   ::  recip_heff
   !$acc declare create(Trans,xyA,recip_heff)

contains

   subroutine mitice_adv_diff(uc, vc)
!  *===========================================================*
!  |    call dt3fl advection routines of mod_csp_thdyn_adv
!  *===========================================================*

      implicit none

!     interface vars
!     uc/vc     :: current ice velocity on C-grid;
      real(wp)   :: uc(:), vc(:)

!     === Local variables ===
!     afx       :: horizontal advective flux, x direction
!     afy       :: horizontal advective flux, y direction
!     gFld      :: tendency of seaice field
      INTEGER i
#ifdef SEAICE_ITD
      INTEGER it
#endif /* SEAICE_ITD */
      real(wp)   ::  afx       (nlpb)
      real(wp)   ::  afy       (nlpb)
      real(wp)   ::  gFld      (nlpb)
      CHARACTER(lc) msgBuf

      allocate(Trans(nlpb,2),xyA(nlpb,2), recip_heff(nlpb))

!$ACC data present(Trans,xyA,recip_heff,dyZ,dxZ,SIMaskU,SIMaskV,uc,vc,    &
#ifdef SEAICE_ITD
!$acc              AREAITD,HEFFITD,HSNOWITD,opnWtrFrac,                &
#endif
#ifdef SEAICE_VARIABLE_SALINITY
!$acc              HSALT,SEAICEdiffKhSalt,              &
#endif
!$acc              AREA,HEFF,HSNOW,HEFFM),         &
!$ACC      create(afx,afy,gFld)

!$acc kernels
!$acc loop
      do i = 1,nlpb
        recip_heff(i) = 1.0_wp
      enddo

      ! Calculate "mass transports" through tracer cell faces.
!$ACC loop
      do i = 1, nlpb
         xyA(i,iu) = dyZ(i)*SIMaskU(i)
         xyA(i,iv) = dxZ(i)*SIMaskV(i)
         Trans(i, iu) = uc(i)*xyA(i,iu)
         Trans(i, iv) = vc(i)*xyA(i,iv)
      enddo
!$acc end kernels

#ifdef SEAICE_ITD
!--   Effective Thickness (Volume)
      IF ( SEAICEadvHeff ) THEN
       DO it=1,SEAICE_multDim
         CALL mitice_adv_dst3fl(HEFFITD(:,it), uc, vc, gFld)
!-    Add tendency due to diffusion
         IF ( SEAICEdiffKhHeff .GT. 0.0_wp ) THEN
           CALL mitice_diff(SEAICEdiffKhHeff,HEFFITD(:,it),gFld)
         ENDIF
!     now do the "explicit" time step
!$acc kernels loop
         do i = 1,nlpb
           HEFFITD(i,it) = HEFFM(i) * (HEFFITD(i,it) + SEAICE_deltaTtherm * gFld(i))
         enddo
!$acc end kernels
       ENDDO
      ENDIF

!--   Fractional area
      IF ( SEAICEadvArea ) THEN
       DO it=1,SEAICE_multDim
          CALL mitice_adv_dst3fl(AREAITD(:,it), uc, vc, gFld)
!-    Add tendency due to diffusion
         IF ( SEAICEdiffKhArea .GT. 0.0_wp ) THEN
            CALL mitice_diff(SEAICEdiffKhArea,AREAITD(:,it),gFld)
         ENDIF
!     now do the "explicit" time step
!$acc kernels loop
         do i = 1,nlpb
           AREAITD(i,it) = HEFFM(i) * (AREAITD(i,it) + SEAICE_deltaTtherm * gFld(i))
         enddo
!$acc end kernels
       ENDDO
!     open water fraction needs to be advected for the ridging scheme
       CALL mitice_adv_dst3fl(opnWtrFrac, uc, vc, gFld)
!--   Add tendency due to diffusion
       IF ( SEAICEdiffKhArea .GT. 0.0_wp ) THEN
         CALL mitice_diff(SEAICEdiffKhArea,opnWtrFrac,gFld)
       ENDIF
!$acc kernels loop
       do i = 1,nlpb
         opnWtrFrac(i) = HEFFM(i) * (opnWtrFrac(i) + SEAICE_deltaTtherm * gFld(i))
       enddo
!$acc end kernels
      ENDIF

!--   Effective Snow Thickness (Volume)
      IF ( SEAICEadvSnow ) THEN
       DO it=1,SEAICE_multDim
         CALL mitice_adv_dst3fl(HSNOWITD(:,it), uc, vc, gFld)
!--   Add tendency due to diffusion
         IF ( SEAICEdiffKhSnow .GT. 0.0_wp ) THEN
           CALL mitice_diff(SEAICEdiffKhSnow,HSNOWITD(:,it),gFld)
         ENDIF
!     now do the "explicit" time step
!$acc kernels loop
         do i = 1,nlpb
           HSNOWITD(i,it) = HEFFM(i) * (HSNOWITD(i,it) + SEAICE_deltaTtherm * gFld(i))
         enddo
!$acc end kernels
       ENDDO
      ENDIF

!C     update mean ice thickness HEFF and total ice concentration AREA
!C     to match single category values
!C     (necessary here because updated HEFF is used below for SItracer)
      call seaice_itd_sum

#else /* not SEAICE_ITD */
!--   Effective Thickness (Volume)
      IF ( SEAICEadvHeff ) THEN
        CALL mitice_adv_dst3fl(HEFF, uc, vc, gFld)
        IF ( SEAICEdiffKhHeff .GT. 0.0_wp ) THEN
!-    Add tendency due to diffusion
          CALL mitice_diff(SEAICEdiffKhHeff,HEFF,gFld)
        ENDIF
!$acc kernels loop
        do i = 1,nlpb
          HEFF(i) = HEFFM(i) * (HEFF(i) + SEAICE_deltaTtherm * gFld(i))
        enddo
!$acc end kernels
      ENDIF

!--   Fractional area
      IF ( SEAICEadvArea ) THEN
        CALL mitice_adv_dst3fl(AREA, uc, vc, gFld)
        IF ( SEAICEdiffKhArea .GT. 0.0_wp ) THEN
!-    Add tendency due to diffusion
          CALL mitice_diff(SEAICEdiffKhArea,AREA,gFld)
        ENDIF
!$acc kernels loop
        do i = 1,nlpb
          AREA(i) = HEFFM(i) * (AREA(i) + SEAICE_deltaTtherm * gFld(i))
        enddo
!$acc end kernels
      ENDIF

!--   Effective Snow Thickness (Volume)
      IF ( SEAICEadvSnow ) THEN
        CALL mitice_adv_dst3fl(HSNOW, uc, vc, gFld)
        IF ( SEAICEdiffKhSnow .GT. 0.0_wp ) THEN
!--   Add tendency due to diffusion
          CALL mitice_diff(SEAICEdiffKhSnow,HSNOW,gFld)
        ENDIF
!$acc kernels loop
        do i = 1,nlpb
          HSNOW(i) = HEFFM(i) * (HSNOW(i) + SEAICE_deltaTtherm * gFld(i))
        enddo
!$acc end kernels
      ENDIF
#endif /* SEAICE_ITD */

#ifdef SEAICE_VARIABLE_SALINITY
!--   Effective Sea Ice Salinity (Mass of salt)
      IF ( SEAICEadvSalt ) THEN
        CALL mitice_adv_dst3fl(HSALT, uc, vc, gFld)
        IF ( SEAICEdiffKhSalt .GT. 0.0_wp ) THEN
!--   Add tendency due to diffusion
          CALL mitice_diff(SEAICEdiffKhSalt,HSALT,gFld)
        ENDIF
!$acc kernels loop
        do i = 1,nlpb
          HSALT(i) = HEFFM(i) * (HSALT(i) + SEAICE_deltaTtherm * gFld(i))
        enddo
!$acc end kernels
      ENDIF
#endif /* SEAICE_VARIABLE_SALINITY */

!$acc end data

      deallocate(Trans,xyA, recip_heff)

   end subroutine mitice_adv_diff


!==============================================================================
   subroutine mitice_adv_dst3fl(tracer,uc, vc, gFld)
       ! tracer should be icetracer or iceFld
!==============================================================================
      implicit none
      integer :: i, k, uel, uei, uep, vnl, vni, vnp, i_max
      integer :: e, w, n, s, w2, s2, km1, km2, kp1
      !interface vars
      real(wp) :: uc(:),vc(:),gFld(:), tracer(:)
      !local vars
      real(wp) :: Rjp, Rj, Rjm
      real(wp) :: uCFL, vCFL, wCFL, d0, d1
      real(wp) :: thetaP, thetaM, psiP, psiM
      real(wp) :: af(nlpb, 2)
      real(wp) :: tl(nlpb)
      real(wp) :: maskLocW(nlpb), maskLocS(nlpb)
      
      !$ACC data present(nk,nlpb,gFld,tracer,Trans,xyA,uc,vc,dyZ,dxZ,drF,hFacW,hFacS,  &
      !$ACC      te,tw,tw2,ue,recip_dxC,recip_dyC,avoid0,recip_heff,  &
      !$ACC      tn,ts,ts2,vn,recip_hFacC,recip_drF,recip_rAc,maskC,rAc,wFld,  &
      !$ACC      mpi_sendbuf_1d,SIMaskU,SIMaskV),  &
      !$ACC      create(tl,af,maskLocW,maskLocS)

      !$ACC kernels
      !$ACC loop
         do i = 1, nlpb
            tl(i) = tracer(i)
         end do
      !$ACC end kernels

      ! Calculate "mass transports" through tracer cell faces.
      !$ACC kernels
      !$ACC loop
      do i = 1, nlpb
         ! advection at horizontal directions
         af(i, iu) = 0.0_wp
         af(i, iv) = 0.0_wp
      end do

!#ifdef ALLOW_OBCS
!         do i = 1, nlpb
!            maskLocW(i) = SIMaskU(i)*maskInW(i)
!            maskLocS(i) = SIMaskV(i)*maskInS(i)
!         enddo 
!#else /* ALLOW_OBCS */
      !$acc loop
      do i = 1, nlpb
         maskLocW(i) = SIMaskU(i)
         maskLocS(i) = SIMaskV(i)
      enddo 
!#endif /* ALLOW_OBCS */

      !$ACC loop private(e,w,w2,Rjp,Rj,Rjm,uCFL,d0,d1,thetaP,thetaM,psiP,psiM,n,s,s2,vCFL)
      do i = 1, nlpb
         e = te(i)
         w = tw(i)
         w2 = tw2(i)

         Rjp = (tl(e) - tl(i))*maskLocW(e)
         Rj = (tl(i) - tl(w))*maskLocW(i)
         Rjm = (tl(w) - tl(w2))*maskLocW(w)

         uCFL = abs(uc(i)*SEAICE_deltaTtherm*recip_dxC(i))

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

         af(i, iu) = 0.5_wp*(Trans(i, iu) + abs(Trans(i, iu)))*(tl(w) + psiP*Rj) &
                     + 0.5_wp*(Trans(i, iu) - abs(Trans(i, iu)))*(tl(i) - psiM*Rj)

         n = tn(i)
         s = ts(i)
         s2 = ts2(i)

         Rjp = (tl(n) - tl(i))*maskLocS(n)
         Rj  = (tl(i) - tl(s))*maskLocS(i)
         Rjm = (tl(s) - tl(s2))*maskLocS(s)

         vCFL = abs(vc(i)*SEAICE_deltaTtherm*recip_dyC(i))
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

         af(i, iv) = 0.5*(Trans(i, iv) + abs(Trans(i, iv)))*(tl(s) + psiP*Rj) &
                     + 0.5*(Trans(i, iv) - abs(Trans(i, iv)))*(tl(i) - psiM*Rj)
      end do

      !$ACC loop private(uel,uei,uep)
      do i = 1, nlpb
         uel = ue(i, 1)
         uei = ue(i, 2)
         uep = ue(i, 3)
         tl(i) = tl(i) &
         !YY: maskC==heffm? r_hFld?  maskC(i,nk) or HEFFM(i)
                 !- SEAICE_deltaTtherm*recip_heff(i)*recip_rAc(i)*HEFFM(i) &
                 !*(af(uel, uei)*uep - af(i, iu) - tracer(i)*(Trans(uel, uei)*uep - Trans(i, iu)))
                 - SEAICE_deltaTtherm*recip_heff(i)*recip_rAc(i)*HEFFM(i) &
                 !- SEAICE_deltaTtherm*recip_heff(i)*recip_rAc(i)*maskC(i,nk) &
                 *(af(uel, uei)*uep - af(i, iu) )

      end do
      !$ACC end kernels

      CALL mpi_data_exchange(mpi_sendbuf_1d, tl)

      ! advection at y direction
      !$ACC kernels
      !$ACC loop private(e,w,w2,Rjp,Rj,Rjm,uCFL,d0,d1,thetaP,thetaM,psiP,psiM,n,s,s2,vCFL)
      do i = 1, nlpb
         af(i, iu) = 0.0_wp
         af(i, iv) = 0.0_wp

         e = te(i)
         w = tw(i)
         w2 = tw2(i)

         Rjp = (tl(e) - tl(i))*maskLocW(e)
         Rj = (tl(i) - tl(w))*maskLocW(i)
         Rjm = (tl(w) - tl(w2))*maskLocW(w)

         uCFL = abs(uc(i)*SEAICE_deltaTtherm*recip_dxC(i))

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

         af(i, iu) = 0.5_wp*(Trans(i, iu) + abs(Trans(i, iu)))*(tl(w) + psiP*Rj) &
                     + 0.5_wp*(Trans(i, iu) - abs(Trans(i, iu)))*(tl(i) - psiM*Rj)

         n = tn(i)
         s = ts(i)
         s2 = ts2(i)

         Rjp = (tl(n) - tl(i))*maskLocS(n)
         Rj  = (tl(i) - tl(s))*maskLocS(i)
         Rjm = (tl(s) - tl(s2))*maskLocS(s)

         vCFL = abs(vc(i)*SEAICE_deltaTtherm*recip_dyC(i))
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

         af(i, iv) = 0.5_wp*(Trans(i, iv) + abs(Trans(i, iv)))*(tl(s) + psiP*Rj) &
                   + 0.5_wp*(Trans(i, iv) - abs(Trans(i, iv)))*(tl(i) - psiM*Rj)
      end do

      !$ACC loop private(vnl,vni,vnp)
      do i = 1, nlpb
         vnl = vn(i, 1)
         vni = vn(i, 2)
         vnp = vn(i, 3)
         tl(i) = tl(i) &
                 !- SEAICE_deltaTtherm*recip_heff(i)*recip_rAc(i)*HEFFM(i) &
                 !*(af(vnl, vni)*vnp - af(i, iv) - tracer(i)*(Trans(vnl, vni)*vnp - Trans(i, iv)))
                 - SEAICE_deltaTtherm*recip_heff(i)*recip_rAc(i)*HEFFM(i) &
                 !- SEAICE_deltaTtherm*recip_heff(i)*recip_rAc(i)*maskC(i,nk) &
                 *(af(vnl, vni)*vnp - af(i, iv))

      end do
      !$ACC end kernels

      ! calculate tracer tendency by advection
      !$ACC kernels loop
      do i = 1, nlpb
         gFld(i) = (tl(i) - tracer(i))/SEAICE_deltaTtherm
      end do
      !$ACC end kernels

      !$ACC end data

   end subroutine mitice_adv_dst3fl


   subroutine mitice_diff(diffKh, iceFld, gFld)
!     *==========================================================*
!     | o Add tendency from horizontal diffusion
!     *==========================================================*

      implicit none

!INTERFACE VARS
      real(wp),intent(  in )  :: iceFld (nlpb)
      real(wp),intent(inout)  :: gFld(nlpb)
      real(wp)  ::  diffKh

!     === Local variables ===
      integer :: i, k, w, s, uel, uei, uep, vnl, vni, vnp
      real(wp) :: dTdl(nlpb, 2), dfl(nlpb, 2)

!$ACC data present(diffKh,iceFld,gFld,xyA,recip_dxC,recip_dyC,ue,vn,recip_rAc,HEFFM,tw,ts),  &
!$ACC      create(dTdl,dfl)

      IF ( diffKh .GT. 0.0_wp ) THEN
         !$ACC kernels
         !$ACC loop private(w,s)
         do i = 1, nlpb
            w = tw(i)
            s = ts(i)
            dTdl(i, iu) = xyA(i, iu)*recip_dxC(i)*(iceFld(i) - iceFld(w))
            dTdl(i, iv) = xyA(i, iv)*recip_dyC(i)*(iceFld(i) - iceFld(s))
         end do
         !$ACC end kernels

         !$ACC kernels
         !$ACC loop
         do i = 1, nlpb
            dfl(i, iu) = -diffKh*dTdl(i, iu)
            dfl(i, iv) = -diffKh*dTdl(i, iv)
         end do
         !$ACC end kernels

         !$ACC kernels loop private(uel,uei,uep,vnl,vni,vnp)
         do i = 1, nlpb
            uel = ue(i, 1)
            uei = ue(i, 2)
            uep = ue(i, 3)
            vnl = vn(i, 1)
            vni = vn(i, 2)
            vnp = vn(i, 3)
            gFld(i) = gFld(i) - recip_rAc(i)*HEFFM(i) &
                         *((dfl(uel,uei)*uep - dfl(i,iu)) + (dfl(vnl,vni)*vnp - dfl(i,iv)))
         end do
         !$ACC end kernels

      ENDIF
      !$ACC end data

   end subroutine mitice_diff

end module mitice_advdiff
