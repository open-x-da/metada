module mod_csp_dyn_ab
   use mod_csp_dyn_mom_vecinv
   use mod_csp_dyn_forward
   use mod_csp_vdiff
   use mod_csp_dyn_force
   use mod_csp_dyn_pbt_exp
   use mod_mpi_test
   implicit none
contains
!==============================================================================
   subroutine csp_dyn_ab
!==============================================================================
      implicit none
      integer :: i, k, bj

      !$ACC kernels present(gUnow,gVnow)
      !$ACC loop collapse(2)
      do k = 1, nk
         do i = 1, nlpb
            gUnow(i, k) = 0.0_wp
            gVnow(i, k) = 0.0_wp
         end do
      end do
      !$ACC end kernels

      call csp_dyn_force

      if (ln_tide) then
         call csp_tide_update
      end if

      call csp_dyn_mom_vecinv

      call csp_dyn_pbt_exp

      call csp_vdiff_implicit(uFldBcl, KappaRM, recip_hFacW, dTmom)

      call csp_vdiff_implicit(vFldBcl, KappaRM, recip_hFacS, dTmom)

      if (ln_bdy .and. mpi_check_bdy) then
         call csp_bdy_orlanski_3d(uFldBcl, bdy_uBcl, uFldBcl_past, maskW, tw, ts, tn, te, dTmom, real(fbvec, wp))
         call csp_bdy_orlanski_3d(vFldBcl, bdy_vBcl, vFldBcl_past, maskS, tw, ts, tn, te, dTmom, real(fbvec, wp))
      end if
      if (ln_bdy) then
         !$ACC kernels present(uFldBcl_past, uFldBcl, vFldBcl_past, vFldBcl)
         !$ACC loop collapse(2) independent
         do k = 1, nk
            do i = 1, nlpb
               uFldBcl_past(i,k) = uFldBcl(i,k)
               vFldBcl_past(i,k) = vFldBcl(i,k)
            end do
         end do
         !$ACC end kernels
      end if

      call csp_dyn_forward_correct_velocity

      call csp_dyn_forward_correct

   end subroutine csp_dyn_ab

end module mod_csp_dyn_ab
