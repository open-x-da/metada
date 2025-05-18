module mod_csp_dyn_phi
   use mod_csp_misc
   use mod_mpi_interfaces
   implicit none
contains

!==============================================================================
   subroutine csp_dyn_phi_hyd(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i
      real(wp) :: ddRloc

      if (ln_bous) then
         ! here alphaRho in or out cycle will effect result
         !$ACC kernels present(alphaRho,maskC,rhoInSitu,rhoLev,phiHydC,phiHydF,drF),  &
         !$ACC         copyin(k)
         !$ACC loop
         do i = 1, nlpb
            alphaRho(i) = maskC(i, k)*rhoInSitu(i, k)
            phiHydC(i) = phiHydF(i) + 0.5_wp*drF(k)*alphaRho(i)*gravity*rhoLev(k)
            phiHydF(i) = phiHydC(i) + 0.5_wp*drF(k)*alphaRho(i)*gravity*rhoLev(k)
         end do
         !$ACC end kernels
      else
         !$ACC kernels present(alphaRho,maskC,rhoInSitu,rhoLev,kSurfC,rSurfC,rC,phiHydC,phiHydF,drF),  &
         !$ACC         copyin(k)
         !$ACC loop independent private(ddRloc)
         do i = 1, nlpb
            alphaRho(i) = maskC(i, k)*(1.0_wp/(rhoInSitu(i, k) + rau0) - rhoLev(k))

            ! This discretization is the "finite volume" form with Part-Cell Bathy
            ! which has not been used to date since it does not
            ! conserve KE+PE exactly even though it is more natural
            ! recommand since it is clear in physics
            if (k .eq. kSurfC(i)) then
               ddRloc = rSurfC(i) - rC(k)
               phiHydC(i) = ddRloc*alphaRho(i)
            else
               phiHydC(i) = phiHydF(i) + 0.5_wp*drF(k)*alphaRho(i)
            end if
            phiHydF(i) = phiHydC(i) + 0.5_wp*drF(k)*alphaRho(i)
         end do
         !$ACC end kernels
      end if

   end subroutine csp_dyn_phi_hyd

!==============================================================================
   subroutine csp_dyn_dphi_hyd(k)
!==============================================================================
      implicit none
      integer, intent(in) :: k
      integer :: i, w, s
      real(wp) :: varLoc(nlpb)
      real(dp) :: sshtemp, sshgloavg, rhsMax, rhsMin, mpi_tmp

      call csp_dyn_phi_hyd(k)

      !$ACC data present(phiHydC,rStarFacC,phi0surf,Bo_surf,ssh,tw,ts,etaH,  &
      !$ACC             recip_RcolC,dPhiHydX,dPhiHydY,recip_dxC,recip_dyC,   &
      !$ACC             alphaRho,rC,maskW,maskS,pbt_base,gravityRau0),  &
      !$ACC      copyin(k),  &
      !$ACC      create(varLoc)

      if (ln_bous) then
         !$ACC kernels
         !$ACC loop
         do i = 1, nlpb
            varLoc(i) = phiHydC(i)*rStarFacC(i) + phi0surf(i)*r1_rau0
         end do

         !$ACC loop independent private(w,s)
         do i = 1, nlpb
            w = tw(i)
            dPhiHydX(i) = (varLoc(i) - varLoc(w))*recip_dxC(i)*maskW(i, k)

            s = ts(i)
            dPhiHydY(i) = (varLoc(i) - varLoc(s))*recip_dyC(i)*maskS(i, k)
         end do
         !$ACC end kernels

         !$ACC kernels
         !$ACC loop
         do i = 1, nlpb
            varLoc(i) = etaH(i)*(1.0_wp + rC(k)*recip_RcolC(i))
         end do

         !$ACC loop independent private(w,s)
         do i = 1, nlpb
            w = tw(i)
            dPhiHydX(i) = dPhiHydX(i) &
                        + (0.5_wp*gravity*r1_rau0*(alphaRho(w) + alphaRho(i)) &  ! these two line is used calculate
                        *(varLoc(i) - varLoc(w)))*recip_dxC(i)*maskW(i, k)      ! slope for z star coordinate

            s = ts(i)
            dPhiHydY(i) = dPhiHydY(i) &
                        + (0.5_wp*gravity*r1_rau0*(alphaRho(s) + alphaRho(i)) &  ! these two line is used calculate
                        *(varLoc(i) - varLoc(s)))*recip_dyC(i)*maskS(i, k)      ! slope for z star coordinate
         end do
         !$ACC end kernels
      else
         !$ACC kernels
         !$ACC loop independent
         do i = 1, nlpb
            varLoc(i) = phiHydC(i)*rStarFacC(i) + phi0surf(i)*Bo_surf(i) + pbt_base(i)*Bo_surf(i)

            if (k .eq. nk) then
               !YY: ask ZY, what ssh mean here? sea surface height contributed by rho (steric)?
               ssh(i) = varLoc(i)/gravity - phi0surf(i)/gravityRau0
            end if
         end do

         !$ACC loop independent private(w,s)
         do i = 1, nlpb
            w = tw(i)
            dPhiHydX(i) = recip_dxC(i)*(varLoc(i) - varLoc(w))*maskW(i, k) &
                          + 0.5_wp*(alphaRho(w) + alphaRho(i)) &  ! these two line is used calculate
                          *(rStarFacC(i) - rStarFacC(w))*rC(k)*recip_dxC(i)  ! slope for p star coordinate

            s = ts(i)
            dPhiHydY(i) = recip_dyC(i)*(varLoc(i) - varLoc(s))*maskS(i, k) &
                          + 0.5_wp*(alphaRho(s) + alphaRho(i)) &  ! these two line is used calculate
                          *(rStarFacC(i) - rStarFacC(s))*rC(k)*recip_dyC(i)  ! slope for p star coordinate
         end do
         !$ACC end kernels
      end if
      !$ACC end data

   end subroutine csp_dyn_dphi_hyd

end module mod_csp_dyn_phi
