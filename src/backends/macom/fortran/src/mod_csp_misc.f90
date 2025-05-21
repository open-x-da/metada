module mod_csp_misc
   use mod_misc_basic
   use mod_csp_basic
   use mod_mpi_variables
   use mod_mpi_interfaces
   implicit none

   real(wp), save, public :: eosMDJWFnum(0:11), eosMDJWFden(0:12) ! coefficient for sea water density calculate
   !$ACC declare create(eosMDJWFnum, eosMDJWFden)

   ! TEOS10/EOS80 parameters
   real(wp),parameter :: rdeltaS = 20.0_wp
   real(wp),parameter :: r1_S0 = 1.0_wp/40.0_wp
   real(wp),parameter :: r1_T0 = 1.0_wp/40.0_wp
   real(wp),parameter :: r1_Z0 = 1.0e-4_wp

   ! EOS parameters
   real(wp),parameter :: EOS000 = 9.5356891948e+02_wp
   real(wp),parameter :: EOS100 = 1.7136499189e+02_wp
   real(wp),parameter :: EOS200 = -3.7501039454e+02_wp
   real(wp),parameter :: EOS300 = 5.1856810420e+02_wp
   real(wp),parameter :: EOS400 = -3.7264470465e+02_wp
   real(wp),parameter :: EOS500 = 1.4302533998e+02_wp
   real(wp),parameter :: EOS600 = -2.2856621162e+01_wp
   real(wp),parameter :: EOS010 = 1.0087518651e+01_wp
   real(wp),parameter :: EOS110 = -1.3647741861e+01_wp
   real(wp),parameter :: EOS210 = 8.8478359933_wp
   real(wp),parameter :: EOS310 = -7.2329388377_wp
   real(wp),parameter :: EOS410 = 1.4774410611_wp
   real(wp),parameter :: EOS510 = 2.0036720553e-01_wp
   real(wp),parameter :: EOS020 = -2.5579830599e+01_wp
   real(wp),parameter :: EOS120 = 2.4043512327e+01_wp
   real(wp),parameter :: EOS220 = -1.6807503990e+01_wp
   real(wp),parameter :: EOS320 = 8.3811577084_wp
   real(wp),parameter :: EOS420 = -1.9771060192_wp
   real(wp),parameter :: EOS030 = 1.6846451198e+01_wp
   real(wp),parameter :: EOS130 = -2.1482926901e+01_wp
   real(wp),parameter :: EOS230 = 1.0108954054e+01_wp
   real(wp),parameter :: EOS330 = -6.2675951440e-01_wp
   real(wp),parameter :: EOS040 = -8.0812310102_wp
   real(wp),parameter :: EOS140 = 1.0102374985e+01_wp
   real(wp),parameter :: EOS240 = -4.8340368631_wp
   real(wp),parameter :: EOS050 = 1.2079167803_wp
   real(wp),parameter :: EOS150 = 1.1515380987e-01_wp
   real(wp),parameter :: EOS060 = -2.4520288837e-01_wp
   real(wp),parameter :: EOS001 = 1.0748601068e+01_wp
   real(wp),parameter :: EOS101 = -1.7817043500e+01_wp
   real(wp),parameter :: EOS201 = 2.2181366768e+01_wp
   real(wp),parameter :: EOS301 = -1.6750916338e+01_wp
   real(wp),parameter :: EOS401 = 4.1202230403_wp
   real(wp),parameter :: EOS011 = -1.5852644587e+01_wp
   real(wp),parameter :: EOS111 = -7.6639383522e-01_wp
   real(wp),parameter :: EOS211 = 4.1144627302_wp
   real(wp),parameter :: EOS311 = -6.6955877448e-01_wp
   real(wp),parameter :: EOS021 = 9.9994861860_wp
   real(wp),parameter :: EOS121 = -1.9467067787e-01_wp
   real(wp),parameter :: EOS221 = -1.2177554330_wp
   real(wp),parameter :: EOS031 = -3.4866102017_wp
   real(wp),parameter :: EOS131 = 2.2229155620e-01_wp
   real(wp),parameter :: EOS041 = 5.9503008642e-01_wp
   real(wp),parameter :: EOS002 = 1.0375676547_wp
   real(wp),parameter :: EOS102 = -3.4249470629_wp
   real(wp),parameter :: EOS202 = 2.0542026429_wp
   real(wp),parameter :: EOS012 = 2.1836324814_wp
   real(wp),parameter :: EOS112 = -3.4453674320e-01_wp
   real(wp),parameter :: EOS022 = -1.2548163097_wp
   real(wp),parameter :: EOS003 = 1.8729078427e-02_wp
   real(wp),parameter :: EOS103 = -5.7238495240e-02_wp
   real(wp),parameter :: EOS013 = 3.8306136687e-01_wp

   ! ALPHA parameters
   real(wp),parameter :: ALP000 = -2.5218796628e-01_wp
   real(wp),parameter :: ALP100 = 3.4119354654e-01_wp
   real(wp),parameter :: ALP200 = -2.2119589983e-01_wp
   real(wp),parameter :: ALP300 = 1.8082347094e-01_wp
   real(wp),parameter :: ALP400 = -3.6936026529e-02_wp
   real(wp),parameter :: ALP500 = -5.0091801383e-03_wp
   real(wp),parameter :: ALP001 = 3.9631611467e-01_wp
   real(wp),parameter :: ALP101 = 1.9159845880e-02_wp
   real(wp),parameter :: ALP201 = -1.0286156825e-01_wp
   real(wp),parameter :: ALP301 = 1.6738969362e-02_wp
   real(wp),parameter :: ALP002 = -5.4590812035e-02_wp
   real(wp),parameter :: ALP102 = 8.6134185799e-03_wp
   real(wp),parameter :: ALP003 = -9.5765341718e-03_wp
   real(wp),parameter :: ALP010 = 1.2789915300_wp
   real(wp),parameter :: ALP110 = -1.2021756164_wp
   real(wp),parameter :: ALP210 = 8.4037519952e-01_wp
   real(wp),parameter :: ALP310 = -4.1905788542e-01_wp
   real(wp),parameter :: ALP410 = 9.8855300959e-02_wp
   real(wp),parameter :: ALP011 = -4.9997430930e-01_wp
   real(wp),parameter :: ALP111 = 9.7335338937e-03_wp
   real(wp),parameter :: ALP211 = 6.0887771651e-02_wp
   real(wp),parameter :: ALP012 = 6.2740815484e-02_wp
   real(wp),parameter :: ALP020 = -1.2634838399_wp
   real(wp),parameter :: ALP120 = 1.6112195176_wp
   real(wp),parameter :: ALP220 = -7.5817155402e-01_wp
   real(wp),parameter :: ALP320 = 4.7006963580e-02_wp
   real(wp),parameter :: ALP021 = 2.6149576513e-01_wp
   real(wp),parameter :: ALP121 = -1.6671866715e-02_wp
   real(wp),parameter :: ALP030 = 8.0812310102e-01_wp
   real(wp),parameter :: ALP130 = -1.0102374985_wp
   real(wp),parameter :: ALP230 = 4.8340368631e-01_wp
   real(wp),parameter :: ALP031 = -5.9503008642e-02_wp
   real(wp),parameter :: ALP040 = -1.5098959754e-01_wp
   real(wp),parameter :: ALP140 = -1.4394226233e-02_wp
   real(wp),parameter :: ALP050 = 3.6780433255e-02_wp

   ! BETA parameters
   real(wp),parameter :: BET000 = 2.1420623987_wp
   real(wp),parameter :: BET100 = -9.3752598635_wp
   real(wp),parameter :: BET200 = 1.9446303907e+01_wp
   real(wp),parameter :: BET300 = -1.8632235232e+01_wp
   real(wp),parameter :: BET400 = 8.9390837485_wp
   real(wp),parameter :: BET500 = -1.7142465871_wp
   real(wp),parameter :: BET010 = -1.7059677327e-01_wp
   real(wp),parameter :: BET110 = 2.2119589983e-01_wp
   real(wp),parameter :: BET210 = -2.7123520642e-01_wp
   real(wp),parameter :: BET310 = 7.3872053057e-02_wp
   real(wp),parameter :: BET410 = 1.2522950346e-02_wp
   real(wp),parameter :: BET020 = 3.0054390409e-01_wp
   real(wp),parameter :: BET120 = -4.2018759976e-01_wp
   real(wp),parameter :: BET220 = 3.1429341406e-01_wp
   real(wp),parameter :: BET320 = -9.8855300959e-02_wp
   real(wp),parameter :: BET030 = -2.6853658626e-01_wp
   real(wp),parameter :: BET130 = 2.5272385134e-01_wp
   real(wp),parameter :: BET230 = -2.3503481790e-02_wp
   real(wp),parameter :: BET040 = 1.2627968731e-01_wp
   real(wp),parameter :: BET140 = -1.2085092158e-01_wp
   real(wp),parameter :: BET050 = 1.4394226233e-03_wp
   real(wp),parameter :: BET001 = -2.2271304375e-01_wp
   real(wp),parameter :: BET101 = 5.5453416919e-01_wp
   real(wp),parameter :: BET201 = -6.2815936268e-01_wp
   real(wp),parameter :: BET301 = 2.0601115202e-01_wp
   real(wp),parameter :: BET011 = -9.5799229402e-03_wp
   real(wp),parameter :: BET111 = 1.0286156825e-01_wp
   real(wp),parameter :: BET211 = -2.5108454043e-02_wp
   real(wp),parameter :: BET021 = -2.4333834734e-03_wp
   real(wp),parameter :: BET121 = -3.0443885826e-02_wp
   real(wp),parameter :: BET031 = 2.7786444526e-03_wp
   real(wp),parameter :: BET002 = -4.2811838287e-02_wp
   real(wp),parameter :: BET102 = 5.1355066072e-02_wp
   real(wp),parameter :: BET012 = -4.3067092900e-03_wp
   real(wp),parameter :: BET003 = -7.1548119050e-04_wp

   real(wp),parameter :: PEN000 = -5.3743005340_wp
   real(wp),parameter :: PEN100 = 8.9085217499_wp
   real(wp),parameter :: PEN200 = -1.1090683384e+01_wp
   real(wp),parameter :: PEN300 = 8.3754581690_wp
   real(wp),parameter :: PEN400 = -2.0601115202_wp
   real(wp),parameter :: PEN010 = 7.9263222935_wp
   real(wp),parameter :: PEN110 = 3.8319691761e-01_wp
   real(wp),parameter :: PEN210 = -2.0572313651_wp
   real(wp),parameter :: PEN310 = 3.3477938724e-01_wp
   real(wp),parameter :: PEN020 = -4.9997430930_wp
   real(wp),parameter :: PEN120 = 9.7335338937e-02_wp
   real(wp),parameter :: PEN220 = 6.0887771651e-01_wp
   real(wp),parameter :: PEN030 = 1.7433051009_wp
   real(wp),parameter :: PEN130 = -1.1114577810e-01_wp
   real(wp),parameter :: PEN040 = -2.9751504321e-01_wp
   real(wp),parameter :: PEN001 = -6.9171176978e-01_wp
   real(wp),parameter :: PEN101 = 2.2832980419_wp
   real(wp),parameter :: PEN201 = -1.3694684286_wp
   real(wp),parameter :: PEN011 = -1.4557549876_wp
   real(wp),parameter :: PEN111 = 2.2969116213e-01_wp
   real(wp),parameter :: PEN021 = 8.3654420645e-01_wp
   real(wp),parameter :: PEN002 = -1.4046808820e-02_wp
   real(wp),parameter :: PEN102 = 4.2928871430e-02_wp
   real(wp),parameter :: PEN012 = -2.8729602515e-01_wp

   ! ALPHA_PEN parameters
   real(wp),parameter :: APE000 = -1.9815805734e-01_wp
   real(wp),parameter :: APE100 = -9.5799229402e-03_wp
   real(wp),parameter :: APE200 = 5.1430784127e-02_wp
   real(wp),parameter :: APE300 = -8.3694846809e-03_wp
   real(wp),parameter :: APE010 = 2.4998715465e-01_wp
   real(wp),parameter :: APE110 = -4.8667669469e-03_wp
   real(wp),parameter :: APE210 = -3.0443885826e-02_wp
   real(wp),parameter :: APE020 = -1.3074788257e-01_wp
   real(wp),parameter :: APE120 = 8.3359333577e-03_wp
   real(wp),parameter :: APE030 = 2.9751504321e-02_wp
   real(wp),parameter :: APE001 = 3.6393874690e-02_wp
   real(wp),parameter :: APE101 = -5.7422790533e-03_wp
   real(wp),parameter :: APE011 = -4.1827210323e-02_wp
   real(wp),parameter :: APE002 = 7.1824006288e-03_wp

   ! BETA_PEN parameters
   real(wp),parameter :: BPE000 = 1.1135652187e-01_wp
   real(wp),parameter :: BPE100 = -2.7726708459e-01_wp
   real(wp),parameter :: BPE200 = 3.1407968134e-01_wp
   real(wp),parameter :: BPE300 = -1.0300557601e-01_wp
   real(wp),parameter :: BPE010 = 4.7899614701e-03_wp
   real(wp),parameter :: BPE110 = -5.1430784127e-02_wp
   real(wp),parameter :: BPE210 = 1.2554227021e-02_wp
   real(wp),parameter :: BPE020 = 1.2166917367e-03_wp
   real(wp),parameter :: BPE120 = 1.5221942913e-02_wp
   real(wp),parameter :: BPE030 = -1.3893222263e-03_wp
   real(wp),parameter :: BPE001 = 2.8541225524e-02_wp
   real(wp),parameter :: BPE101 = -3.4236710714e-02_wp
   real(wp),parameter :: BPE011 = 2.8711395266e-03_wp
   real(wp),parameter :: BPE002 = 5.3661089288e-04_wp

contains

!==============================================================================
   subroutine csp_misc_del2uv(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, w, n, e, s

      !$ACC data present(zn,ze,tw,ts,del2u,del2v,hDiv,recip_dxC,recip_dyC,hFacZ,  &
      !$ACC         vort3,recip_dyZ,recip_dxZ,recip_hFacW,recip_hFacS,maskW,maskS)

      !$ACC kernels
      !$ACC loop independent private(w,n,e,s)
      do i = 1, nlpb
         w = tw(i)
         n = zn(i)
         del2u(i) = ((hDiv(i) - hDiv(w))*recip_dxC(i) &
                     - (hFacZ(n)*vort3(n) - hFacZ(i)*vort3(i)) &
                     *recip_dyZ(i)*recip_hFacW(i, k))*maskW(i, k)

         e = ze(i)
         s = ts(i)
         del2v(i) = ((hDiv(i) - hDiv(s))*recip_dyC(i) &
                     + (hFacZ(e)*vort3(e) - hFacZ(i)*vort3(i)) &
                     *recip_dxZ(i)*recip_hFacS(i, k))*maskS(i, k)
      end do
      !$ACC end kernels
      !$ACC end data

   end subroutine csp_misc_del2uv

!==============================================================================
   subroutine csp_misc_hdiv(uFldLoc, vFldLoc, hDivLoc)
!==============================================================================
      implicit none
      integer :: i, uel, uei, uep, vnl, vni, vnp
      real(wp), intent(in) :: uFldLoc(nlpb), vFldLoc(nlpb)
      real(wp) :: vecFldLoc(nlpb, 2), dxyZ(nlpbz, 2)
      real(wp), intent(out) :: hDivLoc(nlpb)

      !$ACC kernels present(dyZ,dxZ,ue,vn,recip_rAc,uFldLoc,vFldLoc,hDivLoc),  &
      !$ACC         create(vecFldLoc,dxyz,hDivLoc)
      !$ACC loop independent
      do i = 1, nlpbz
         dxyZ(i, iu) = dyZ(i)
         dxyZ(i, iv) = dxZ(i)
      end do

      !$ACC loop independent
      do i = 1, nlpb
         vecFldLoc(i, iu) = uFldLoc(i)
         vecFldLoc(i, iv) = vFldLoc(i)
      end do

      !$ACC loop independent private(uel,uei,uep,vnl,vni,vnp)
      do i = 1, nlpb
         uel = ue(i, 1)
         uei = ue(i, 2)
         uep = ue(i, 3)
         vnl = vn(i, 1)
         vni = vn(i, 2)
         vnp = vn(i, 3)
         hDivLoc(i) = (vecFldLoc(uel, uei)* uep * dxyZ(uel, uei) - vecFldLoc(i, iu)* dxyZ(i, iu) &
                     + vecFldLoc(vnl, vni)* vnp * dxyZ(vnl, vni) - vecFldLoc(i, iv)* dxyZ(i, iv))* recip_rAc(i)
      end do
      !$ACC end kernels

   end subroutine csp_misc_hdiv

!==============================================================================
   subroutine csp_misc_adams_bashforth2(gnowLoc, gpastLoc, FirstStep)
!==============================================================================
      implicit none
      integer :: i
      logical :: FirstStep
      real(wp) :: AB_gTr(nlpb)
      real(wp), intent(inout) :: gnowLoc(nlpb), gpastLoc(nlpb)

      if (FirstStep) then
         !$ACC kernels present(gnowLoc,gpastLoc),  &
         !$ACC         create(AB_gTr)
         !$ACC loop independent
         do i = 1, nlpb
            gpastLoc(i) = gnowLoc(i)
         end do
         !$ACC end kernels
      end if

      !$ACC kernels present(abFac) create(AB_gTr)
      !$ACC loop independent
      do i = 1, nlpb
         AB_gTr(i) = abFac*(gnowLoc(i) - gpastLoc(i))
         gpastLoc(i) = gnowLoc(i)
         gnowLoc(i) = gnowLoc(i) + AB_gTr(i)
      end do
      !$ACC end kernels

      return
   end subroutine csp_misc_adams_bashforth2

!==============================================================================
   subroutine csp_misc_rho_init
!==============================================================================
      implicit none

      eosMDJWFnum(0) = 9.99843699e+02_wp
      eosMDJWFnum(1) = 7.35212840e+00_wp
      eosMDJWFnum(2) = -5.45928211e-02_wp
      eosMDJWFnum(3) = 3.98476704e-04_wp
      eosMDJWFnum(4) = 2.96938239e+00_wp
      eosMDJWFnum(5) = -7.23268813e-03_wp
      eosMDJWFnum(6) = 2.12382341e-03_wp
      eosMDJWFnum(7) = 1.04004591e-02_wp
      eosMDJWFnum(8) = 1.03970529e-07_wp
      eosMDJWFnum(9) = 5.18761880e-06_wp
      eosMDJWFnum(10) = -3.24041825e-08_wp
      eosMDJWFnum(11) = -1.23869360e-11_wp

      eosMDJWFden(0) = 1.00000000e+00_wp
      eosMDJWFden(1) = 7.28606739e-03_wp
      eosMDJWFden(2) = -4.60835542e-05_wp
      eosMDJWFden(3) = 3.68390573e-07_wp
      eosMDJWFden(4) = 1.80809186e-10_wp
      eosMDJWFden(5) = 2.14691708e-03_wp
      eosMDJWFden(6) = -9.27062484e-06_wp
      eosMDJWFden(7) = -1.78343643e-10_wp
      eosMDJWFden(8) = 4.76534122e-06_wp
      eosMDJWFden(9) = 1.63410736e-09_wp
      eosMDJWFden(10) = 5.30848875e-06_wp
      eosMDJWFden(11) = -3.03175128e-16_wp
      eosMDJWFden(12) = -1.27934137e-17_wp

      !$ACC update device(eosMDJWFnum, eosMDJWFden)

   end subroutine csp_misc_rho_init

!==============================================================================
   subroutine csp_misc_rhoInSitu
!==============================================================================
      implicit none
      integer :: i, k
      real(wp) :: rhoDen, rhoNum, locPres
      real(wp) :: t1, t2, s1, p1, sp5, gp1

      if (ln_bous) then
         gp1 = rau0*gravity
      else
         gp1 = 1.0_wp
      end if

      !$ACC kernels present(nk,nlpb,rC,tFld,eosMDJWFnum,eosMDJWFden,rhoInSitu,maskC),copyin(gp1)
      !$ACC loop collapse(2) independent private(rhoDen, rhoNum, locPres, t1, t2, s1, p1, sp5)
      do k = 1, nk
         do i = 1, nlpb
            locPres = rC(k)*gp1

            t1 = tFld(i, k, 1)
            t2 = t1*t1
            s1 = tFld(i, k, 2)

            if (s1 .gt. 0.0_wp) then
               sp5 = sqrt(s1)
            else
               s1 = 0.0_wp
               sp5 = 0.0_wp
            end if

            p1 = locPres*SItodBar

            rhoNum = eosMDJWFnum(0) &
                     + t1*(eosMDJWFnum(1) + t1*(eosMDJWFnum(2) + eosMDJWFnum(3)*t1)) &
                     + s1*(eosMDJWFnum(4) + eosMDJWFnum(5)*t1 + eosMDJWFnum(6)*s1) &
                     + p1*(eosMDJWFnum(7) + eosMDJWFnum(8)*t2 &
                           + eosMDJWFnum(9)*s1 + p1*(eosMDJWFnum(10) + eosMDJWFnum(11)*t2))

            rhoDen = eosMDJWFden(0) &
                     + t1*(eosMDJWFden(1) + t1*(eosMDJWFden(2) + t1*(eosMDJWFden(3) + t1*eosMDJWFden(4)))) &
                     + s1*(eosMDJWFden(5) + t1*(eosMDJWFden(6) + eosMDJWFden(7)*t2) &
                           + sp5*(eosMDJWFden(8) + eosMDJWFden(9)*t2)) &
                     + p1*(eosMDJWFden(10) + p1*t1*(eosMDJWFden(11)*t2 + eosMDJWFden(12)*p1))

            rhoDen = 1.0_wp/rhoDen

            rhoInSitu(i, k) = (rhoDen*rhoNum - rau0)*maskC(i, k)

         end do
      end do
      !$ACC end kernels

      return
   end subroutine csp_misc_rhoInSitu

!==============================================================================
   subroutine csp_misc_bn2
!==============================================================================
      implicit none
      integer :: i, k
      real(wp) :: zh, zt, zs
      real(wp) :: zn, zn0, zn1, zn2, zn3
      real(wp) :: zaw, zbw, zrw
      real(wp) :: pab(nlpb, nk, 2)                  ! thermal/haline expansion ratio

      !$ACC kernels create(pab),  &
      !$ACC         present(rC,rF,tFld,rn2,drF,maskC,nk,nl,nlpb, nkp1, gravityRau0)
      !$ACC loop collapse(2) independent private(zh,zt,zs,zn3,zn2,zn1,zn0,zn)
      do k = 1, nk
         do i = 1, nlpb
            pab(i, k, 1) = 0.0_wp
            pab(i, k, 2) = 0.0_wp
            zh = rC(k)*r1_Z0/gravityRau0              ! depth
            zt = tFld(i, k, 1)*r1_T0                            ! temperature
            zs = sqrt(abs(tFld(i, k, 2) + rdeltaS)*r1_S0)   ! square root salinity
            !
            ! alpha
            zn3 = ALP003
            !
            zn2 = ALP012*zt + ALP102*zs + ALP002
            !
            zn1 = ((ALP031*zt &
                    + ALP121*zs + ALP021)*zt &
                   + (ALP211*zs + ALP111)*zs + ALP011)*zt &
                  + ((ALP301*zs + ALP201)*zs + ALP101)*zs + ALP001

            zn0 = ((((ALP050*zt &
                      + ALP140*zs + ALP040)*zt &
                     + (ALP230*zs + ALP130)*zs + ALP030)*zt &
                    + ((ALP320*zs + ALP220)*zs + ALP120)*zs + ALP020)*zt &
                   + (((ALP410*zs + ALP310)*zs + ALP210)*zs + ALP110)*zs + ALP010)*zt &
                  + ((((ALP500*zs + ALP400)*zs + ALP300)*zs + ALP200)*zs + ALP100)*zs + ALP000

            zn = ((zn3*zh + zn2)*zh + zn1)*zh + zn0
            !
            pab(i, k, 1) = zn*r1_rau0
            !
            ! beta
            zn3 = BET003
            !
            zn2 = BET012*zt + BET102*zs + BET002
            !
            zn1 = ((BET031*zt &
                    + BET121*zs + BET021)*zt &
                   + (BET211*zs + BET111)*zs + BET011)*zt &
                  + ((BET301*zs + BET201)*zs + BET101)*zs + BET001

            zn0 = ((((BET050*zt &
                      + BET140*zs + BET040)*zt &
                     + (BET230*zs + BET130)*zs + BET030)*zt &
                    + ((BET320*zs + BET220)*zs + BET120)*zs + BET020)*zt &
                   + (((BET410*zs + BET310)*zs + BET210)*zs + BET110)*zs + BET010)*zt &
                  + ((((BET500*zs + BET400)*zs + BET300)*zs + BET200)*zs + BET100)*zs + BET000

            zn = ((zn3*zh + zn2)*zh + zn1)*zh + zn0
            !
            pab(i, k, 2) = zn/zs*r1_rau0
         end do
      end do

      !$ACC loop collapse(2)
      do k = 1, nkp1
         do i = 1, nlpb
            rn2(i, k) = 0.0_wp
         end do
      end do

      !$ACC loop collapse(2) independent private(zrw,zaw,zbw)
      do k = nk, 2, -1
         do i = 1, nlpb
            zrw = (rC(k - 1) - rF(k))/(rC(k - 1) - rC(k))

            zaw = pab(i, k - 1, 1)*(1.0_wp - zrw) + pab(i, k, 1)*zrw
            zbw = pab(i, k - 1, 2)*(1.0_wp - zrw) + pab(i, k, 2)*zrw

            rn2(i, k) = gravityRau0*gravity* &
                        (zaw*(tFld(i, k, 1) - tFld(i, k - 1, 1)) - zbw*(tFld(i, k, 2) - tFld(i, k - 1, 2))) &
                        /drF(k)*maskC(i, k - 1)
         end do
      end do
      !$ACC end kernels
      return
   end subroutine csp_misc_bn2

end module mod_csp_misc
