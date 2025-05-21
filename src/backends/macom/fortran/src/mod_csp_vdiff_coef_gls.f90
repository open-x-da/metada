module mod_csp_vdiff_coef_gls
   use mod_csp_basic
   use mod_csp_misc
   use mod_csp_vdiff
   use mod_mpi_test
   implicit none
   private

   public csp_vdiff_coef_gls
   public csp_vdiff_coef_gls_init

   ! parameters will move to namelist read in future
   real(wp), parameter :: rn_crban = 100.0_wp   ! Craig and Banner 1994 constant for wb tke flux
   real(wp), parameter :: rn_z0 = 3.0e-3_wp  ! bottom roughness [m]
   real(wp), parameter :: rn_frac_hs = 1.3_wp   ! Fraction of wave height as roughness
   real(wp), parameter :: rn_hsro = 0.02_wp  !  Minimum surface roughness [m]
   real(wp), parameter :: sigma_en = 1.0_wp
   real(wp), parameter :: sigma_psi = 1.2_wp
   real(wp), parameter :: c_mu_0 = 0.5268_wp
   real(wp), parameter :: rn_p = 3.0_wp
   real(wp), parameter :: rn_m = 1.5_wp
   real(wp), parameter :: rn_n = -1.0_wp
   real(wp), parameter :: c_psi1 = 1.44_wp
   real(wp), parameter :: c_psi2 = 1.92_wp
   real(wp), parameter :: c_psi3_s = -0.629_wp
   real(wp), parameter :: c_psi3_us = 1.0_wp
   real(wp), parameter :: rghmin = -0.28_wp
   real(wp), parameter :: rgh0 = 0.0329_wp
   real(wp), parameter :: rl1 = 0.107_wp
   real(wp), parameter :: rl2 = 0.0032_wp
   real(wp), parameter :: rl3 = 0.0864_wp
   real(wp), parameter :: rl4 = 0.12_wp
   real(wp), parameter :: rl5 = 11.9_wp
   real(wp), parameter :: rl6 = 0.4_wp
   real(wp), parameter :: rl7 = 0.0_wp
   real(wp), parameter :: rl8 = 0.48_wp
   real(wp), parameter :: r2_3 = 2.0_wp/3.0_wp
   real(wp), parameter :: rcm_sf = 0.731_wp
   real(wp) :: ra_sf = -2.0_wp   ! Must be negative -2 < ra_sf < -1
   real(wp) :: rl_sf = 0.2_wp
   real(wp) :: rtrans = 0.1_wp
   real(wp) :: sigma_psi0
   real(wp) :: rs0
   real(wp) :: rs1
   real(wp) :: rs2
   real(wp) :: rs4
   real(wp) :: rs5
   real(wp) :: rs6
   real(wp) :: r_wp
   real(wp) :: rd1
   real(wp) :: rd2
   real(wp) :: rd3
   real(wp) :: rd4
   real(wp) :: rd5
   real(wp) :: rc02
   real(wp) :: rc03
   real(wp) :: rc04
   real(wp) :: rc02r
   real(wp) :: rsbc_tke
   real(wp) :: rsbc_zs2
   real(wp) :: rf6
   real(wp) :: rc_diff
   real(wp) :: psimin1
   real(wp) :: rough_b                                ! bottom roughness

   ! ---- OPENACC
   !$ACC declare copyin(rl_sf,ra_sf,rtrans)
   !$ACC declare create(rs0,rs1,rs2,rs4,rs5,rs6,r_wp,rd1,rd2,rd3,rd4,rd5,   &
   !$ACC                rc02,rc03,rc04,rc02r,rsbc_tke,rsbc_zs2,rf6,rc_diff, &
   !$ACC                psimin1,sigma_psi0,rough_b)
   ! ----
contains

!==============================================================================
   subroutine csp_vdiff_coef_gls_init
!==============================================================================
      implicit none
      integer :: i, k
      real(wp) :: zcr

      ! start constant parameter set, these will move out for fast in future
      rs0 = 1.5_wp*rl1*rl5*rl5
      rs1 = -rl4*(rl6 + rl7) + 2.0_wp*rl4*rl5*(rl1 - (1.0_wp/3.0_wp)*rl2 - rl3) + 1.5_wp*rl1*rl5*rl8
      rs2 = -(3.0_wp/8.0_wp)*rl1*(rl6*rl6 - rl7*rl7)
      rs4 = 2.0_wp*rl5
      rs5 = 2.0_wp*rl4
      rs6 = (2.0_wp/3.0_wp)*rl5*(3.0_wp*rl3*rl3 - rl2*rl2) - 0.5_wp*rl5*rl1*(3.0_wp*rl3 - rl2) &
            + 0.75_wp*rl1*(rl6 - rl7)
      r_wp = 3.0_wp*rl5*rl5
      rd1 = rl5*(7.0_wp*rl4 + 3.0_wp*rl8)
      rd2 = rl5*rl5*(3.0_wp*rl3*rl3 - rl2*rl2) - 0.75_wp*(rl6*rl6 - rl7*rl7)
      rd3 = rl4*(4.0_wp*rl4 + 3.0_wp*rl8)
      rd4 = rl4*(rl2*rl6 - 3.0_wp*rl3*rl7 - rl5*(rl2*rl2 - rl3*rl3)) + rl5*rl8*(3.0_wp*rl3*rl3 - rl2*rl2)
      rd5 = 0.25_wp*(rl2*rl2 - 3.0_wp*rl3*rl3)*(rl6*rl6 - rl7*rl7)
      rf6 = 8.0_wp/(c_mu_0**6)
      rc_diff = sqrt(2.0_wp)/(c_mu_0**3)

      rc02 = c_mu_0*c_mu_0; rc02r = 1.0_wp/rc02
      rc03 = rc02*c_mu_0
      rc04 = rc03*c_mu_0

      sigma_psi0 = sigma_psi

      ra_sf = -4.0_wp*rn_n*sqrt(sigma_en)/((1.0_wp + 4.0_wp*rn_m)*sqrt(sigma_en) &
                                           - sqrt(sigma_en + 24.0_wp*sigma_psi0*c_psi2))
      if (rn_crban == 0.0_wp) then
         rl_sf = vkarmn
      else
         rl_sf = c_mu_0*sqrt(c_mu_0/rcm_sf)*sqrt(((1.0_wp + 4.0_wp*rn_m + 8.0_wp*rn_m**2)*sigma_en &
                                                  + 12.0_wp*sigma_psi0*c_psi2 - (1.0_wp + 4.0_wp*rn_m) &
                                                  *sqrt(sigma_en*(sigma_en + 24.0_wp*sigma_psi0*c_psi2)))/(12.0_wp*rn_n**2))
      end if

      rsbc_tke = -3.0_wp/2.0_wp*rn_crban*ra_sf*rl_sf
      zcr = max(rsmall, rsbc_tke**(1.0_wp/(-ra_sf*3.0_wp/2.0_wp)) - 1.0_wp)
      rtrans = 0.2_wp/zcr                                              ! Ad. inverse transition length between log and wave layer
      rsbc_zs2 = rn_frac_hs/0.85_wp/gravity*665.0_wp                  ! Rascle formula for surface roughness

      psimin1 = rn_psimin

      rough_b = rn_z0*gravityRau0  ! convert bottom roughness from m into Pa

      if (.not. restart_in) then
         do k = 1, nkp1
            do i = 1, nlpb
               hmxl(i, k) = 0.05_wp*gravityRau0
            end do
         end do
      end if

      ! ---- OPENACC
      !$ACC update device(rl_sf, ra_sf, rtrans,  &
      !$ACC               rs0,rs1,rs2,rs4,rs5,rs6,r_wp,rd1,rd2,rd3,rd4,rd5,   &
      !$ACC               rc02,rc03,rc04,rc02r,rsbc_tke,rsbc_zs2,rf6,rc_diff, &
      !$ACC               psimin1,sigma_psi0,rough_b)
      ! ----

   end subroutine csp_vdiff_coef_gls_init

!==============================================================================
   subroutine csp_vdiff_coef_gls
!==============================================================================
      !!----------------------------------------------------------------------
      !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
      !!              coefficients using the GLS turbulent closure scheme.
      !!----------------------------------------------------------------------
      implicit none
      integer  :: i, k, uel, uei, uep, vnl, vni, vnp
      integer  :: ibot, ibotp1                           ! local integers
      real(wp) :: prod_en(nlpb, nkp1)
      real(wp) :: buo_en(nlpb, nkp1)
      real(wp) :: diss_en(nlpb, nkp1)
      real(wp) :: psi(nlpb, nkp1)
      real(wp) :: u_fric_s(nlpb)                         ! friction velocity at surface ocean side
      real(wp) :: u_fric_b                               ! friction velocity at bottom ocean side
      real(wp) :: u_fric_sa(nlpb)                        ! friction velocity at surface air side
      real(wp) :: zpelc(nlpb, nkp1), zhlc(nlpb)           ! Max depth of Langmuir cell [Pa]
      real(wp) :: plc(nlpb, nkp1)                         ! TKE due to Langmuir cell [Pa^2/s^3]
      real(wp) :: eb(nlpb, nkp1), hmxlb(nlpb, nkp1)
      real(wp) :: c_mu, c_mu_rho
      real(wp) :: work_c(nlpb, nkp1)
      real(wp) :: zzz, zcd                               !   -         -
      real(wp) :: zut, zvt                               !   -         -
      real(wp) :: zdir, zdiss, zprod
      real(wp) :: wlc
      real(wp) :: c_psi_3
      real(wp) :: zcof, gh, gm, sm, sh, rcff
      real(wp) :: uv(nlpb, nk, 2)
      real(wp) :: rCdU_bot                               ! top/bottom drag coeff. at t-point (<0) [m/s]
      real(wp) :: wave_age
      real(wp) :: rough_s(nlpb)                          ! surface roughness
      real(wp) :: zdep, zkar
      real(wp) :: kappa(nlpb, nk)
      real(wp) :: temp_sfc, minrn2
      real(wp) :: tprod, tbuo
      real(wp) :: wmaskC(nlpb, nkp1)
      real(wp) :: sh2(nlpb, nkp1)                         ! shear production term (M^2) [1/s^2]

      !$ACC data create(prod_en, buo_en, diss_en, psi, u_fric_s, u_fric_sa, &
      !$ACC             zpelc, zhlc, plc, eb, hmxlb, work_c, uv, rough_s, &
      !$ACC             kappa, wmaskC, sh2),   &
      !$ACC     present(en,nl,nlpb,nkp1,nk,uFld,vFld,drC,kSurfC,maskC,hmxl,  &
      !$ACC             KappaRM,KappaRT,rn2,taum,  &
      !$ACC             hFacC,rF,drF,rC,dTtracer,ue,vn,BottomDragFac,  &
      !$ACC             BottomDrag,BottomDragMax,uDragBottom,  &
      !$ACC             KappaRMCon,KappaRTCon)

      !$ACC kernels
      !$ACC loop
      do i = 1, nlpb
         wmaskC(i, nkp1) = 0.0_wp
         wmaskC(i, kSurfC(i)) = 0.0_wp
      end do

      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            wmaskC(i, k) = maskC(i, k)
            sh2(i, k) = 0.0_wp
            uv(i, k, 1) = uFld(i, k)
            uv(i, k, 2) = vFld(i, k)
         end do
      end do

      ! compute M^2(sh2) at w grid,
      ! Notice both of them is zero at sea-land surface
      ! en and diss at top and bottom surface will calculate by boundary condition
      !$ACC loop collapse(2)
      do k = 2, nk
         do i = 1, nlpb
            sh2(i, k) = ((uFld(i, k) - uFld(i, k - 1))*gravityRau0/drC(k))**2 &
                        + ((vFld(i, k) - vFld(i, k - 1))*gravityRau0/drC(k))**2
            sh2(i, k) = sh2(i, k)*wmaskC(i, k)
         end do
      end do
      !$ACC end kernels

      ! compute N^2(rn2) at w grid
      call csp_misc_bn2

      ! calculate turbulent energy production due to shear P(prod_en) and buoyancy B(buo_en)
      ! Notice rn2 can be negative here, which represent produce en by buoyancy instability
      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nkp1
         do i = 1, nlpb
            prod_en(i, k) = KappaRM(i, k)*sh2(i, k)
            buo_en(i, k) = -KappaRT(i, k)*rn2(i, k)

            ! calculate dissipation of turbulent energy
            diss_en(i, k) = rc03*en(i, k)*sqrt(en(i, k))/hmxl(i, k)

            eb(i, k) = en(i, k)
            hmxlb(i, k) = hmxl(i, k)

            ! initial plc for Langmuir circulation source term
            plc(i, k) = 0.0_wp
         end do
      end do

      ! set Langmuir circulation source term
      !$ACC loop
      do i = 1, nlpb
         u_fric_sa(i) = 0.377_wp*sqrt(taum(i))*gravityRau0        ! friction velocity at surface air side [Pa/s]
         zhlc(i) = rF(nk)

         zpelc(i, nkp1) = 0.0_wp
         zpelc(i, nk) = rn2(i, nk)*rF(nk)*drF(nk)

         !$ACC loop seq
         do k = nk - 1, 1, -1
            zpelc(i, k) = zpelc(i, k + 1) + rn2(i, k)*rF(k)*drF(k)
         end do

         !$ACC loop seq
         do k = 1, nk
            if (zpelc(i, k) > 0.5_wp*u_fric_sa(i)*u_fric_sa(i)) then
               zhlc(i) = rF(k)
            end if
         end do
      end do

      !$ACC loop collapse(2) private(wlc)
      do k = 1, nk
         do i = 1, nlpb
            if (rC(k) < zhlc(i)) then
               wlc = 0.15_wp*u_fric_sa(i)*sin(rpi*rC(k)/zhlc(i))
               plc(i, k) = wlc*wlc*wlc/zhlc(i)
            end if
         end do
      end do

      ! commit this to close or open LC
      !$ACC loop collapse(2)
      do k = 1, nkp1
         do i = 1, nlpb
            plc(i, k) = 0.0_wp
         end do
      end do
      ! end Langmuir circulation source term

      ! now solver turbulent energy transport equation through TMDA method
      ! notice we take diss_en as part of the LHS of equation
      !$ACC loop collapse(2) private(zdiss,zprod,zdir)
      do k = 1, nk - 1
         do i = 1, nlpb
            ! prod_en + buo_en < 0, buo_en will reduce TKE as diss, zdir = 0;
            ! prod_en + buo_en > 0, buo_en will increase TKE as shear production, zdir = 1
            zdir = 0.5_wp + sign(0.5_wp, buo_en(i, k) + prod_en(i, k))

            zdiss = diss_en(i, k) - (1.0_wp - zdir)*buo_en(i, k)
            zprod = zdir*buo_en(i, k) + prod_en(i, k) + plc(i, k)

            work_c(i, k) = zdiss/en(i, k)

            en(i, k) = en(i, k) + dTtracer*zprod
         end do
      end do
      !$ACC end kernels

      ! set surface boundary friction velocity and roughness
      !$ACC kernels
      !$ACC loop independent private(wave_age)
      do i = 1, nlpb
         u_fric_s(i) = sqrt(abs(taum(i))*r1_rau0)  ! unit is m/s
         wave_age = 30.0_wp*tanh(0.6_wp/28.0_wp/u_fric_s(i))  ! unit is dimensionless
         rough_s(i) = rsbc_zs2*u_fric_s(i)*u_fric_s(i)*wave_age*sqrt(wave_age)  ! unit is m
         rough_s(i) = max(rough_s(i), rn_hsro) ! unit is m
         rough_s(i) = rough_s(i)*gravityRau0  ! convert unit to Pa

         ! set surface boundary for en
         en(i, nkp1) = rc02r*u_fric_s(i)*u_fric_s(i)*(1.0_wp + rsbc_tke)**r2_3 ! unit is m^2/s^2
         en(i, nkp1) = en(i, nkp1)*gravityRau0*gravityRau0*maskC(i, nk) ! convert unit to Pa^2/s^2
      end do
      !$ACC end kernels

      ! set wave break penetrate, now we set only add at level nk, in future this may change to nk -- nk-x
      do k = nk, nk
         !$ACC kernels
         !$ACC loop independent private(temp_sfc)
         do i = 1, nlpb
            temp_sfc = rc02r*u_fric_s(i)*u_fric_s(i)*(1.0_wp + rsbc_tke* &
                       ((rough_s(i) + rF(k)*hFacC(i, k))/rough_s(i))**(1.5_wp*ra_sf))**r2_3 ! unit is m^2/s^2
            temp_sfc = temp_sfc*gravityRau0*gravityRau0*maskC(i, k - 1) ! convert unit to Pa^2/s^2

            en(i, k) = temp_sfc
         end do
         !$ACC end kernels
      end do
      ! end wave break penetrate

      ! set bottom boundary for en
      !$ACC kernels
      !$ACC loop independent private(uel,uei,uep,vnl,vni,vnp,ibot,ibotp1,zut,zvt,zzz,zcd,rCdU_bot,u_fric_b)
      do i = 1, nlpb
         uel = ue(i, 1)
         uei = ue(i, 2)
         uep = ue(i, 3)
         vnl = vn(i, 1)
         vni = vn(i, 2)
         vnp = vn(i, 3)
         ibot = min(nk, kSurfC(i))          ! ocean bottom level at w-points
         zut = 0.5_wp*(uv(i, ibot, iu) + uv(uel, ibot, uei)*uep)     ! ?2 x velocity at t-point?
         zvt = 0.5_wp*(uv(i, ibot, iv) + uv(vnl, ibot, vni)*vnp)

         zzz = max(1.0_wp, 0.5_wp*drF(ibot)*hFacC(i, ibot))   ! unit is Pa

         zcd = (vkarmn/log(zzz/rough_b))**2  ! unit is dimensionless

         zcd = BottomDragFac(i)*min(max(BottomDrag, zcd), BottomDragMax)   ! unit is dimensionless

         rCdU_bot = zut*zut + zvt*zvt

         uDragBottom(i) = -zcd*sqrt(rCdU_bot)

         rCdU_bot = -zcd*rCdU_bot  ! we cancel out rho_0 here for fast

         u_fric_b = sqrt(abs(rCdU_bot))  ! unit is m/s

         en(i, ibot) = rc02r*u_fric_b*u_fric_b*gravityRau0*gravityRau0 ! unit is Pa^2/s^2

         ibotp1 = min(nkp1, ibot + 1)

         en(i, ibotp1) = en(i, ibot)
      end do
      ! end set bottom boundary for en

      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            kappa(i, k) = 0.5_wp*(KappaRM(i, k) + KappaRM(i, k + 1))/sigma_en
         end do
      end do
      !$ACC end kernels

      ! now solver TKE transport equation through LU decompose method
      call csp_vdiff_implicit_gls(en, kappa, work_c)

      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nkp1
         do i = 1, nlpb
            en(i, k) = max(en(i, k), rn_emin)

            ! set psi equal diss_en
            psi(i, k) = diss_en(i, k)
         end do
      end do

      !$ACC loop collapse(2) private(tprod,tbuo,zdiss,zprod,zdir,c_psi_3)
      do k = 1, nk - 1
         do i = 1, nlpb
            ! rn2 > 0, stable, zdir = 1, c_psi_3 set to c_psi3_s;
            ! rn2 < 0, unstable, zdir = 0, c_psi_3 set to c_psi3_us
            zdir = 0.5_wp + sign(0.5_wp, rn2(i, k))
            c_psi_3 = zdir*c_psi3_s + (1.0_wp - zdir)*c_psi3_us

            tprod = c_psi1*prod_en(i, k)
            tbuo = c_psi_3*buo_en(i, k)

            ! tprod + tbuo < 0, buo_en will reduce psi as diss, zdir = 0;
            ! tprod + tbuo > 0, buo_en will increase psi as shear production, zdir = 1
            zdir = 0.5_wp + sign(0.5_wp, tprod + tbuo)
            zdiss = c_psi2*diss_en(i, k) - (1.0_wp - zdir)*tbuo
            zprod = tprod + zdir*tbuo

            work_c(i, k) = zdiss/eb(i, k)
            psi(i, k) = psi(i, k) + dTtracer*zprod*psi(i, k)/eb(i, k)
         end do
      end do

      ! set psi surface boundary condition
      !$ACC loop independent private(zdep,zkar)
      do i = 1, nlpb
         ! boundary condition at surface (nkp1)
         zdep = rough_s(i)*rl_sf
         psi(i, nkp1) = c_mu_0**rn_p*en(i, nkp1)**rn_m*zdep**rn_n*maskC(i, nk)
         zkar = (rl_sf + (vkarmn - rl_sf)*(1.0_wp - exp(-rtrans*rF(nk)/rough_s(i))))

         ! boundary condition at level below surface (nk)
         zdep = (rough_s(i) + rF(nk))*zkar
         psi(i, nk) = c_mu_0**rn_p*en(i, nk)**rn_m*zdep**rn_n*maskC(i, nk - 1)
      end do

      ! set psi bottom boundary condition
      !$ACC loop independent private(ibot,ibotp1,zdep)
      do i = 1, nlpb
         ! boundary condition at bottom
         ibot = min(nk, kSurfC(i))          ! ocean bottom level at w-points
         zdep = vkarmn*rough_b
         psi(i, ibot) = c_mu_0**rn_p*en(i, ibot)**rn_m*zdep**rn_n

         ! boundary condition at level up bottom
         ibotp1 = min(nkp1, ibot + 1)
         zdep = vkarmn*(rough_b + drF(ibot))
         psi(i, ibotp1) = c_mu_0**rn_p*en(i, ibot)**rn_m*zdep**rn_n
      end do

      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            kappa(i, k) = 0.5_wp*(KappaRM(i, k) + KappaRM(i, k + 1))/sigma_psi
         end do
      end do
      !$ACC end kernels

      ! now solver psi transport equation through LU decompose method
      call csp_vdiff_implicit_gls(psi, kappa, work_c)

      !$ACC kernels
      !$ACC loop collapse(2)
      do k = 1, nkp1
         do i = 1, nlpb
            psi(i, k) = max(psi(i, k), psimin1)

            ! update diss_en with new psi
            diss_en(i, k) = psi(i, k)

            hmxl(i, k) = rc03*en(i, k)*sqrt(en(i, k))/diss_en(i, k)
            ! hmxl = min(  0.56_wp * sqrt( 2.0_wp * en / rn2 ), hmxl )
            ! hmxl = max(0.05_wp * gravity * rau0, hmxl)

            work_c(i, k) = hmxl(i, k)*sqrt(2.0_wp*en(i, k))
         end do
      end do

      ! stable function
      !$ACC loop collapse(2) private(c_mu,c_mu_rho,zcof,gh,gm,rcff,sm,sh,minrn2)
      do k = 2, nk
         do i = 1, nl
            ! zcof =  l²/q²
            zcof = hmxlb(i, k)*hmxlb(i, k)/(2.0_wp*eb(i, k))
            ! Gh = -N²l²/q²
            gh = -rn2(i, k)*zcof
            gh = min(gh, rgh0)
            gh = max(gh, rghmin)
            gh = gh*rf6
            ! Gm =  M²l²/q² Shear number
            gm = max(sh2(i, k)*zcof, 1.0e-10_wp)
            gm = gm*rf6
            gm = min((r_wp - rd1*gh + rd3*gh*gh)/(rd2 - rd4*gh), gm)
            ! Stability functions from Canuto
            rcff = r_wp - rd1*gh + rd2*gm + rd3*gh*gh - rd4*gh*gm + rd5*gm*gm
            sm = (rs0 - rs1*gh + rs2*gm)/rcff
            sh = (rs4 - rs5*gh + rs6*gm)/rcff
            !
            ! Store stability function in zstt and zstm
            c_mu = rc_diff*sh*maskC(i, k)
            c_mu_rho = rc_diff*sm*maskC(i, k)

            KappaRM(i, k) = c_mu*work_c(i, k)
            KappaRT(i, k) = c_mu_rho*work_c(i, k)

            KappaRM(i, k) = KappaRM(i, k)*wmaskC(i, k)
            KappaRT(i, k) = KappaRT(i, k)*wmaskC(i, k)

            KappaRM(i, k) = max(KappaRM(i, k), KappaRMCon)
            KappaRT(i, k) = max(KappaRT(i, k), KappaRTCon)
            
            minrn2 = min(min(rn2(i, k-1), rn2(i, k)), rn2(i, k+1))
            if (minrn2 .lt. 0.0_wp) then
               KappaRT(i, k) = max(KappaRT(i, k), KappaRTevd)
            end if
            if (minrn2 .lt. -0.00001_wp) then
               KappaRT(i, k) = max(KappaRT(i, k), KappaRTevd*10.0_wp)
            end if
            if (minrn2 .lt. -0.0001_wp) then
               KappaRT(i, k) = max(KappaRT(i, k), KappaRTevd*100.0_wp)
            end if
            if (minrn2 .lt. -0.001_wp) then
               KappaRT(i, k) = max(KappaRT(i, k), KappaRTevd*1000.0_wp)
            end if
         end do
      end do
      !$ACC end kernels
      !$ACC end data

   end subroutine csp_vdiff_coef_gls

end module mod_csp_vdiff_coef_gls
