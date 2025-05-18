module mod_csp_bdy
   use mod_misc_basic
   use mod_csp_basic
   use mod_csp_misc
   implicit none

contains

!==============================================================================
   subroutine csp_bdy_relax_2d(model2d, bdy2d)
!==============================================================================
      !!---------------------------------------------------------------------
      !!
      !!
      !!
      !!
      !!---------------------------------------------------------------------
      implicit none
      integer :: i, i_loc, d
      real(wp) :: d_r, alpha_d
      real(wp) :: model2d(nlpb), bdy2d(nlbdy)

      do i = 1, nlbdy
         i_loc = mpi_idx_bdy(i)
         d = maskBdy(i_loc)
         if (d .gt. 2) then
            d_r = (real(d, wp)-1)/2.0_wp
            alpha_d = 1.0_wp - tanh(d_r)
         else
            alpha_d = 1.0_wp
         end if
         model2d(i_loc) = alpha_d*bdy2d(i) + (1.0_wp-alpha_d)*model2d(i_loc)
      end do

   end subroutine csp_bdy_relax_2d

!==============================================================================
   subroutine csp_bdy_relax_3d(model3d, bdy3d)
!==============================================================================
      !!---------------------------------------------------------------------
      !!
      !!
      !!
      !!
      !!---------------------------------------------------------------------
      implicit none
      integer :: k
      real(wp) :: model3d(nlpb,nk), bdy3d(nlbdy,nk)

      do k = 1, nk
         call csp_bdy_relax_2d(model3d(:,k), bdy3d(:,k))
      end do

   end subroutine csp_bdy_relax_3d

!==============================================================================
   subroutine csp_bdy_orlanski_2d(model2d, bdy2d, model2d_past, mask, xw, xs, xn, xe, dTin, dTout)
!==============================================================================
      !!---------------------------------------------------------------------
      !!
      !!
      !!
      !!
      !!---------------------------------------------------------------------
      implicit none
      integer :: i, j
      integer :: bj, bm1j, bm2j, bjp1, bjm1, bm1jp1, bm1jm1
      integer :: normalIdx
      integer, intent(in) :: xw(nlpb), xs(nlpb), xn(nlpb), xe(nlpb)
      integer(i2), intent(in) :: mask(nlpb)
      real(wp), intent(in) :: bdy2d(nlbdy), model2d_past(nlpb)
      real(wp), intent(inout) :: model2d(nlpb)
      real(wp) :: cnormal, ctange, dphit, dphix, dphiy, sdphiy, party
      real(wp) :: dTratio
      real(wp), intent(in) :: dTin, dTout

      dTratio = dTin/dTout

      do j = bdywidth, 1, -1
         ! loop for noraml point
         !$ACC kernels present(model2d, bdy2d, model2d_past, mask, xw, xs, xn, &
         !$ACC                 xe, mpi_idx_bdy, normal_idx, maskBdy), &
         !$ACC          copyin(dTratio, j)
         !$ACC loop independent private(bj, normalIdx, bm1j, bm2j, bjp1, bjm1, &
         !$ACC      bm1jp1, bm1jm1, dphit, dphix, sdphiy, dphiy, sdphiy, cnormal, &
         !$ACC      ctange, party)
         do i = 1, nlbdy
            bj = mpi_idx_bdy(i)
            normalIdx = normal_idx(bj)
            if ((maskBdy(bj) .eq. j) .and. (normalIdx .le. 4)) then
               !if (mpi_rank == 1) write(*,*) i, j, normalIdx, maskBdy(idxBdy)
               select case (normalIdx)
               case(1)
               ! normal direction is north, tangential is east
                  bm1j = xs(bj)
                  bm2j = xs(bm1j)
                  bjp1 = xe(bj)
                  bjm1 = xw(bj)
                  bm1jp1 = xe(bm1j)
                  bm1jm1 = xw(bm1j)

               case(2)
               ! normal direction is east, tangential is south
                  bm1j = xw(bj)
                  bm2j = xw(bm1j)
                  bjp1 = xs(bj)
                  bjm1 = xn(bj)
                  bm1jp1 = xs(bm1j)
                  bm1jm1 = xn(bm1j)

               case(3)
               ! normal direction is south, tangential is west
                  bm1j = xn(bj)
                  bm2j = xn(bm1j)
                  bjp1 = xw(bj)
                  bjm1 = xe(bj)
                  bm1jp1 = xw(bm1j)
                  bm1jm1 = xe(bm1j)

               case(4)
               ! normal direction is west, tangential is north
                  bm1j = xe(bj)
                  bm2j = xe(bm1j)
                  bjp1 = xn(bj)
                  bjm1 = xs(bj)
                  bm1jp1 = xn(bm1j)
                  bm1jm1 = xs(bm1j)

               case default
                  continue
               end select

               if (normalIdx .eq. 0) then
                  model2d(bj) = bdy2d(i)
               else
                  dphit = model2d(bm1j) - model2d_past(bm1j)
                  dphix = (model2d(bm1j) - model2d(bm2j))*mask(bm1j)*mask(bm2j)
                  sdphiy = dphit*(model2d_past(bm1jp1) - model2d_past(bm1jm1))*mask(bm1jp1)*mask(bm1jm1)

                  dphiy = 0.0_wp
                  if (sdphiy .gt. 0) then
                     dphiy = (model2d_past(bm1j) - model2d_past(bm1jm1))*mask(bm1j)*mask(bm1jm1)
                  else
                     dphiy = (model2d_past(bm1jp1) - model2d_past(bm1j))*mask(bm1jp1)*mask(bm1j)
                  end if

                  sdphiy = dphix**2+dphiy**2
                  if (sdphiy .gt. 1.0e-5_wp) then
                     cnormal = -dphit*dphix/sdphiy
                     ctange = -dphit*dphiy/sdphiy
                  else
                     cnormal = 0.0_wp
                     ctange = 0.0_wp
                  end if

                  if (cnormal .gt. 0.0_wp) then
                     party = 0.0_wp
                     if (ctange .gt. 0.0_wp) then
                        party = - ctange*(model2d_past(bj) - model2d_past(bjm1))*mask(bj)*mask(bjm1)
                     else
                        party = - ctange*(model2d_past(bjp1) - model2d_past(bj))*mask(bj)*mask(bjp1)
                     end if

                     model2d(bj) = (model2d_past(bj) + cnormal*model2d(bm1j)*mask(bm1j) + party)/(1.0_wp+cnormal)

                     model2d(bj) = model2d(bj) - (model2d(bj) - bdy2d(i))*dTratio
                  else
                     model2d(bj) = bdy2d(i)
                  end if
               end if

               !if ((mpi_rank .eq. 0 ) .and. (i .eq. 500)) then
               !   write(*,*) model2d(bj), party, cnormal, ctange, bdy2d(i), sdphiy, dphiy, dphix, dphit
               !   stop
               !end if
            end if
         end do

         ! loop for corner point
         !$ACC loop independent private(bj, normalIdx)
         do i = 1, nlbdy
            bj = mpi_idx_bdy(i)
            normalIdx = normal_idx(bj)
            if ((maskBdy(bj) .eq. j) .and. (normalIdx .gt. 4)) then
               model2d(bj) = bdy2d(i)
            end if
         end do
         !$ACC end kernels
      end do ! end do j = 1, bdywidth

   end subroutine csp_bdy_orlanski_2d

!==============================================================================
   subroutine csp_bdy_orlanski_3d(model3d, bdy3d, model3d_past, mask, xw, xs, xn, xe, dTin, dTout)
!==============================================================================
      !!---------------------------------------------------------------------
      !!
      !!
      !!
      !!
      !!---------------------------------------------------------------------
      implicit none
      integer :: i, j, k
      integer :: bj, bm1j, bm2j, bjp1, bjm1, bm1jp1, bm1jm1
      integer :: normalIdx
      integer, intent(in) :: xw(nlpb), xs(nlpb), xn(nlpb), xe(nlpb)
      integer(i2), intent(in) :: mask(nlpb, nk)
      real(wp), intent(in) :: bdy3d(nlbdy,nk), model3d_past(nlpb,nk)
      real(wp), intent(inout) :: model3d(nlpb,nk)
      real(wp) :: cnormal, ctange, dphit, dphix, dphiy, sdphiy, party
      real(wp) :: dTratio
      real(wp), intent(in) :: dTin, dTout

      dTratio = dTin/dTout

      do j = bdywidth, 1, -1
         ! loop for noraml point
         !$ACC kernels present(model3d, bdy3d, model3d_past, mask, xw, xs, xn, &
         !$ACC                 xe, mpi_idx_bdy, normal_idx, maskBdy), &
         !$ACC          copyin(dTratio, j)
         !$ACC loop independent private(bj, normalIdx, bm1j, bm2j, bjp1, bjm1, &
         !$ACC      bm1jp1, bm1jm1, dphit, dphix, sdphiy, dphiy, sdphiy, cnormal, &
         !$ACC      ctange, party)
         do i = 1, nlbdy
            bj = mpi_idx_bdy(i)
            normalIdx = normal_idx(bj)
            if ((maskBdy(bj) .eq. j) .and. (normalIdx .le. 4)) then
               !if (mpi_rank == 1) write(*,*) i, j, normalIdx, maskBdy(idxBdy)
               select case (normalIdx)
               case(1)
               ! normal direction is north, tangential is east
                  bm1j = xs(bj)
                  bm2j = xs(bm1j)
                  bjp1 = xe(bj)
                  bjm1 = xw(bj)
                  bm1jp1 = xe(bm1j)
                  bm1jm1 = xw(bm1j)

               case(2)
               ! normal direction is east, tangential is south
                  bm1j = xw(bj)
                  bm2j = xw(bm1j)
                  bjp1 = xs(bj)
                  bjm1 = xn(bj)
                  bm1jp1 = xs(bm1j)
                  bm1jm1 = xn(bm1j)

               case(3)
               ! normal direction is south, tangential is west
                  bm1j = xn(bj)
                  bm2j = xn(bm1j)
                  bjp1 = xw(bj)
                  bjm1 = xe(bj)
                  bm1jp1 = xw(bm1j)
                  bm1jm1 = xe(bm1j)

               case(4)
               ! normal direction is west, tangential is north
                  bm1j = xe(bj)
                  bm2j = xe(bm1j)
                  bjp1 = xn(bj)
                  bjm1 = xs(bj)
                  bm1jp1 = xn(bm1j)
                  bm1jm1 = xs(bm1j)

               case default
                  continue
               end select

               !$ACC loop seq
               do k = 1, nk
                  if (normalIdx .eq. 0) then
                     model3d(bj,k) = bdy3d(i,k)
                  else
                     dphit = model3d(bm1j,k) - model3d_past(bm1j,k)
                     dphix = (model3d(bm1j,k) - model3d(bm2j,k))*mask(bm1j,k)*mask(bm2j,k)
                     sdphiy = dphit*(model3d_past(bm1jp1,k) - model3d_past(bm1jm1,k))*mask(bm1jm1,k)*mask(bm1jp1,k)

                     dphiy = 0.0_wp
                     if (sdphiy .gt. 0) then
                        dphiy = (model3d_past(bm1j,k) - model3d_past(bm1jm1,k))*mask(bm1j,k)*mask(bm1jm1,k)
                     else
                        dphiy = (model3d_past(bm1jp1,k) - model3d_past(bm1j,k))*mask(bm1j,k)*mask(bm1jp1,k)
                     end if

                     sdphiy = dphix**2+dphiy**2
                     if (sdphiy .gt. 1.0e-5_wp) then
                        cnormal = -dphit*dphix/sdphiy
                        ctange = -dphit*dphiy/sdphiy
                      else
                        cnormal = 0.0_wp
                        ctange = 0.0_wp
                     end if

                     if (cnormal .gt. 0.0_wp) then
                        party = 0.0_wp
                        if (ctange .gt. 0.0_wp) then
                           party = - ctange*(model3d_past(bj,k) - model3d_past(bjm1,k))*mask(bj,k)*mask(bjm1,k)
                        else
                           party = - ctange*(model3d_past(bjp1,k) - model3d_past(bj,k))*mask(bj,k)*mask(bjp1,k)
                        end if

                        model3d(bj,k) = (model3d_past(bj,k) + cnormal*model3d(bm1j,k)*mask(bm1j,k) + party)/(1.0_wp+cnormal)

                        model3d(bj,k) = model3d(bj,k) - (model3d(bj,k) - bdy3d(i,k))*dTratio
                     else
                        model3d(bj,k) = bdy3d(i,k)
                     end if
                  end if
               end do
            end if
         end do

         ! loop for corner point
         !$ACC loop independent private(bj, normalIdx)
         do i = 1, nlbdy
            bj = mpi_idx_bdy(i)
            normalIdx = normal_idx(bj)
            if ((maskBdy(bj) .eq. j) .and. (normalIdx .gt. 4)) then
               !$ACC loop seq
               do k = 1, nk
                  model3d(bj,k) = bdy3d(i,k)
               end do
            end if
         end do
         !$ACC end kernels

      end do ! end do j = 1, bdywidth

   end subroutine csp_bdy_orlanski_3d

end module mod_csp_bdy
