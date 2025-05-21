module mod_csp_dyn_force
   use mod_csp_misc
   implicit none
contains

!==============================================================================
   subroutine csp_dyn_force
!==============================================================================
      implicit none
      integer :: i
      real(dp) :: slpglosum, slpgloavg

      !$ACC kernels present(nlpb,nk,guExt,gvExt,utau,vtau,maskW,maskS,drF)
      !$ACC loop independent
      do i = 1, nlpb
      !!YY, ZY: phi0surf should put before sea ice model, So moved to mpi_csp_io_rdforce end
         !phi0surf(i) = slp(i, 4)
!!YY: sIceLoad how to add to surface pressure, here is example from MITgcm
!!phi0surf(i,j,bi,bj) = ( pLoad(i,j,bi,bj) +sIceLoad(i,j,bi,bj)*gravity*sIceLoadFac )*recip_rhoConst

         guExt(i) = utau(i)*maskW(i, nk)*gravityDynForce/drF(nk)

         gvExt(i) = vtau(i)*maskS(i, nk)*gravityDynForce/drF(nk)

      end do
      !$ACC end kernels

# if !defined (NOSLP)
      ! remove global mean slp effect in phi0surf
      !if (.not. ln_bdy) then

      !   slpglosum = 0.0d0
      !   !$ACC kernels present(loc_nlpb,phi0surf,rAc,maskC),copy(slpglosum)
      !   !$ACC loop reduction(+:slpglosum)
      !   do i = 1, loc_nlpb
      !      slpglosum = slpglosum + dble(phi0surf(i))* dble(rAc(i))*maskC(i, nk)
      !   end do
      !   !$ACC end kernels

      !   CALL mpi_comp_dp_allsum(slpglosum)

      !   slpgloavg = slpglosum / globalArea

         !!$ACC kernels present(nlpb,phi0surf),copyin(slpgloavg)
         !!$ACC loop independent private(slpglosum)
         !$ACC kernels present(nlpb, phi0surf)
         !$ACC loop independent
         do i = 1, nlpb
            !phi0surf(i) = phi0surf(i) - slpgloavg
            phi0surf(i) = phi0surf(i) - 101325.0_wp
         end do
         !$ACC end kernels
      !end if
# endif

   end subroutine csp_dyn_force

end module mod_csp_dyn_force
