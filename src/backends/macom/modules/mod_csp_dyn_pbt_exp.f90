module mod_csp_dyn_pbt_exp
   use mod_csp_bdy
   use mod_misc
# if defined (PTIDE)
   use mod_csp_dyn_ptide
# endif
   use mod_mpi_test
   implicit none

   integer, save, public :: iterSurf
!
contains

!==============================================================================
   subroutine csp_dyn_pbt_exp
!==============================================================================
      implicit none
      integer ::  i, k, w, s, bj
      real(wp) :: iterSurfs_r, beta_exp, eps_exp
      real(wp) :: etaH_avg(nlpb), uBar_avg(nlpb), vBar_avg(nlpb), etaH_star(nlpb)
      real(wp) :: uBar_star(nlpb), vBar_star(nlpb), uBar_temp(nlpb), vBar_temp(nlpb)
      real(wp) :: eta_temp(nlpb)
      real(dp) :: sshtemp, sshgloavg, rhsMax, rhsMin, mpi_tmp

# if defined (PTIDE)
      call csp_dyn_ptide
# endif

      iterSurfs_r = real(iterSurfs, wp)

      beta_exp = 1.0_wp/3.0_wp
      eps_exp = 2.0_wp/3.0_wp

      !$ACC data create(etaH_avg, uBar_avg, vBar_avg, etaH_star, uBar_star, &
      !$ACC             vBar_star, uBar_temp, vBar_temp, eta_temp),   &
      !$ACC     present(nlpb, dTsurf, etaN, etaH, uBar, vBar, tw, ts, dphiSurfX,  &
      !$ACC             dphiSurfY, recip_dxC, recip_dyC, Bo_surf, uBarCouple,  &
      !$ACC             vBarCouple, bdy_pbt, bdy_uBar, bdy_vBar, maskC, maskW, &
# if defined (PTIDE)
      !$ACC             BAREL, &
# endif
      !$ACC             maskS, tn, te, mpi_sendbuf_1d), &
      !$ACC      copyin(iterSurfs_r, beta_exp, eps_exp)

      !$ACC kernels
      !$ACC loop
      do i = 1, nlpb
         etaH_avg(i) = 0.0_wp
         uBar_avg(i) = 0.0_wp
         vBar_avg(i) = 0.0_wp
      end do
      !$ACC end kernels

      do iterSurf = 1, Itersurf_en !nIterSurf

         !$ACC kernels
         !$ACC loop
         do i = 1, nlpb
            etaN(i) = etaH(i)
         end do
         !$ACC end kernels

         call csp_dyn_pbt_exp_etaH(etaN, uBar, vBar, etaH_star)

         !$ACC kernels
         !$ACC loop
         do i = 1, nlpb
            eta_temp(i) = beta_exp*etaH_star(i) + (1.0_wp-beta_exp)*etaN(i)
         end do
         !$ACC end kernels

         CALL mpi_data_exchange(mpi_sendbuf_1d, eta_temp)
         
         !$ACC kernels
# if defined (PTIDE)
         !$ACC loop independent private(w,s)
         do i = 1, nlpb
            w = tw(i)
            dphiSurfX(i) = recip_dxC(i)*(Bo_surf(i)*(eta_temp(i)-BAREL(i)) - Bo_surf(w)*(eta_temp(w)-BAREL(w)))
            s = ts(i)
            dphiSurfY(i) = recip_dyC(i)*(Bo_surf(i)*(eta_temp(i)-BAREL(i)) - Bo_surf(s)*(eta_temp(s)-BAREL(s)))
         end do
# else
         !$ACC loop independent private(w,s)
         do i = 1, nlpb
            w = tw(i)
            dphiSurfX(i) = recip_dxC(i)*(Bo_surf(i)*eta_temp(i) - Bo_surf(w)*eta_temp(w))
            s = ts(i)
            dphiSurfY(i) = recip_dyC(i)*(Bo_surf(i)*eta_temp(i) - Bo_surf(s)*eta_temp(s))
         end do
# endif

         !$ACC loop independent
         do i = 1, nlpb
            if (maskBdy(i) .eq. 0) then 
               uBar_star(i) = uBar(i) + dTsurf*(uBarCouple(i) - dphiSurfX(i))
               vBar_star(i) = vBar(i) + dTsurf*(vBarCouple(i) - dphiSurfY(i))
               uBar_temp(i) = 0.5_wp*(uBar_star(i) + uBar(i))
               vBar_temp(i) = 0.5_wp*(vBar_star(i) + vBar(i))
            else
               uBar_temp(i) = uBar(i)
               vBar_temp(i) = vBar(i)
            end if
         end do
         !$ACC end kernels

         call csp_dyn_pbt_exp_etaH(etaN, uBar_temp, vBar_temp, etaH)

         !$ACC kernels
         !$ACC loop
         do i = 1, nlpb
            eta_temp(i) = 0.5_wp*(eps_exp*etaH(i) + (1.0_wp-eps_exp)*etaH_star(i)+etaN(i))
         end do
         !$ACC end kernels

         CALL mpi_data_exchange(mpi_sendbuf_1d, eta_temp)

         !$ACC kernels
# if defined (PTIDE)
         !$ACC loop independent private(w,s)
         do i = 1, nlpb
            w = tw(i)
            dphiSurfX(i) = recip_dxC(i)*(Bo_surf(i)*(eta_temp(i)-BAREL(i)) - Bo_surf(w)*(eta_temp(w)-BAREL(w)))
            s = ts(i)
            dphiSurfY(i) = recip_dyC(i)*(Bo_surf(i)*(eta_temp(i)-BAREL(i)) - Bo_surf(s)*(eta_temp(s)-BAREL(s)))
         end do
# else
         !$ACC loop independent private(w,s)
         do i = 1, nlpb
            w = tw(i)
            dphiSurfX(i) = recip_dxC(i)*(Bo_surf(i)*eta_temp(i) - Bo_surf(w)*eta_temp(w))
            s = ts(i)
            dphiSurfY(i) = recip_dyC(i)*(Bo_surf(i)*eta_temp(i) - Bo_surf(s)*eta_temp(s))
         end do
# endif
         !$ACC end kernels

         !$ACC kernels
         !$ACC loop independent
         do i = 1, nlpb
            ! store last uBar and vBar in uBar_temp and vBar_temp
            uBar_temp(i) = uBar(i)
            vBar_temp(i) = vBar(i)
            if (maskBdy(i) .eq. 0) then 
               uBar(i) = uBar_temp(i) + dTsurf*(uBarCouple(i) - dphiSurfX(i))
               vBar(i) = vBar_temp(i) + dTsurf*(vBarCouple(i) - dphiSurfY(i))
            end if
         end do
         !$ACC end kernels

         if (ln_bdy .and. mpi_check_bdy) then
            call csp_bdy_orlanski_2d(etaH, bdy_pbt(:,4), etaN, maskC(:,nk), tw, ts, tn, te, dTsurf, 3600.0_wp)
            call csp_bdy_orlanski_2d(uBar, bdy_uBar, uBar_temp, maskW(:,nk), tw, ts, tn, te, dTsurf, 3600.0_wp)
            call csp_bdy_orlanski_2d(vBar, bdy_vBar, vBar_temp, maskS(:,nk), tw, ts, tn, te, dTsurf, 3600.0_wp)
         end if

         CALL mpi_data_exchange(mpi_sendbuf_1d, etaH)
         CALL mpi_data_exchange(mpi_sendbuf_1d, uBar)
         CALL mpi_data_exchange(mpi_sendbuf_1d, vBar)

         if (iterSurf .ge. Itersurf_st) then
            !$ACC kernels
            !$ACC loop independent
            do i = 1, nlpb
               etaH_avg(i) = etaH_avg(i) + etaH(i)
               uBar_avg(i) = uBar_avg(i) + uBar(i)
               vBar_avg(i) = vBar_avg(i) + vBar(i)
            end do
            !$ACC end kernels
         end if

      end do

      !$ACC kernels
      !$ACC loop
      do i = 1, nlpb
         etaH(i) = etaH_avg(i) / iterSurfs_r
         uBar(i) = uBar_avg(i) / iterSurfs_r
         vBar(i) = vBar_avg(i) / iterSurfs_r
      end do
      !$ACC end kernels

      !$ACC end data

      !$ACC update self(etaH)
      sshtemp = 0.0d0
      do i = 1, loc_nlpb
         sshtemp = sshtemp + dble(etaH(i))* dble(rAc(i))* maskC(i, nk)
      end do
      CALL mpi_comp_dp_allsum(sshtemp)
      sshgloavg = sshtemp / globalArea
      rhsMax = maxval(etaH)
      CALL MPI_ALLREDUCE(rhsMax, mpi_tmp, 1, mpi_real_wp, MPI_MAX, mpi_comp_comm, mpi_err)
      rhsMax = mpi_tmp
      rhsMin = minval(etaH)
      CALL MPI_ALLREDUCE(rhsMin, mpi_tmp, 1, mpi_real_wp, MPI_MIN, mpi_comp_comm, mpi_err)
      rhsMin = mpi_tmp
      if (mpi_rank == 0) write (sol_out_unit, "(a,i8,a,e23.16,a,e23.16,a,e23.16)") &
         "avg of etaH after step ", myIter, " is ", sshgloavg, " MAX is ", rhsMax, " MIN is ", rhsMin

   end subroutine csp_dyn_pbt_exp

!==============================================================================
  subroutine csp_dyn_pbt_exp_etaH(etaH0, uBar0, vBar0, etaH1)
!==============================================================================
     implicit none
     integer :: i, w, s, uel, uei, uep, vnl, vni, vnp
     real(wp), intent(in) :: uBar0(nlpb), vBar0(nlpb), etaH0(nlpb)
     real(wp) :: vecTrans(nlpb, 2), conv2d(nlpb), dEtaHdt(nlpb)
     real(wp), intent(inout) :: etaH1(nlpb)
     real(wp) :: etaHw, etaHs, RcolWb, RcolSb
     real(dp) :: sshtemp

     !$ACC data create(vecTrans, conv2d, dEtaHdt), &
     !$ACC         present(uBar0, vBar0, etaH0, etaH1, dxZ, dyZ, RcolC, maskW, &
     !$ACC                 maskS, maskC, recip_rAc, EmPmR, ue, vn, tw, ts, &
     !$ACC                 dTsurf, gravityDynForce, RcolW, RcolS, maskBdy)

     !$ACC kernels
     !$ACC loop independent private(w, s, etaHw, etaHs, RcolWb, RcolSb)
     do i = 1, nlpb
        w = tw(i)
        s = ts(i)
        etaHw = 0.5_wp*(etaH0(i)+etaH0(w))
        etaHs = 0.5_wp*(etaH0(i)+etaH0(s))
        RcolWb = RcolW(i)
        RcolSb = RcolS(i)
        vecTrans(i, iu) = uBar0(i)*dyZ(i)*(RcolWb + etaHw)*maskW(i, nk)
        vecTrans(i, iv) = vBar0(i)*dxZ(i)*(RcolSb + etaHs)*maskS(i, nk)
     end do

     !$ACC loop independent private(uel,uei,uep,vnl,vni,vnp)
     do i = 1, nlpb
        if (maskBdy(i) .eq. 0) then 
           uel = ue(i, 1)
           uei = ue(i, 2)
           uep = ue(i, 3)
           vnl = vn(i, 1)
           vni = vn(i, 2)
           vnp = vn(i, 3)
           conv2d(i) = -maskC(i, nk) &
                       *(vecTrans(uel, uei)*uep - vecTrans(i, iu) + &
                         vecTrans(vnl, vni)*vnp - vecTrans(i, iv))
   
           dEtaHdt(i) = conv2d(i)*recip_rAc(i)
   
           etaH1(i) = etaH0(i) + dEtaHdt(i)*dTsurf*maskC(i, nk) &
                               - EmPmR(i)*dTsurf*gravityDynForce
        else
           etaH1(i) = etaH0(i)
        end if
     end do
     !$ACC end kernels

     !$ACC end data

  end subroutine csp_dyn_pbt_exp_etaH

end module mod_csp_dyn_pbt_exp
