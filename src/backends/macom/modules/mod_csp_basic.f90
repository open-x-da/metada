module mod_csp_basic
   use mod_misc_basic
   implicit none

!YY, ZY, add a switch for seaice module. use ifdef later
   logical, public          :: mitice_on
   !logical, public          :: mitice_on=.true.
   !logical, public          :: mitice_on=.false.

   integer, save, public          :: ndomain
   integer, save, public          :: idom
   integer, save, public          :: parent
   integer, save, public          :: nRatio
   integer, save, public          :: nl
   integer, save, public          :: nlpb
   integer, save, public          :: loc_nlpb ! the number of only local grids, added by pangrb
   integer, save, public          :: nlz
   integer, save, public          :: nlpbz
   integer, save, public          :: nk
   integer, save, public          :: nlbdy
   integer, save, public          :: bdywidth
   integer, save, public          :: loc_nlbdy
   integer, save, public          :: sbct = 4 ! the vertical number for reading forcing input data, added by pangrb
   integer, save, public          :: nkp1
   integer, save, public          :: ni
   real(wp), save, public         :: dTtracer
   real(wp), save, public         :: dTmom
   real(wp), save, public         :: dTsurf
   integer, save, public          :: nIterSurf
   integer, save, public          :: iterSurf_st
   integer, save, public          :: iterSurf_en
   integer, save, public          :: iterSurfs
   integer, save, public          :: nIter0
   integer, save, public          :: nIterMax
   integer, save, public          :: nOutFreq
   character(lc), save, public    :: fOutFreq
   integer, save, public          :: nResFreq
   integer, save, public          :: startyear, startmon, startday, starthour
   integer, save, public          :: startmin, startsec
   integer, save, public          :: nowyearf, nowmonf, nowdayf, nowhourf
   integer, save, public          :: nowminf, nowsecf
   integer, save, public          :: nowyear, nowmon, nowday, nowhour
   integer, save, public          :: nowmin, nowsec
   integer(i8), save, public      :: nowjulian
   integer(i8), save, public      :: endjulian
   !YY
   integer(i8), save, public      :: startTime     !time in sec when modeling begins
   integer(i8), save, public      :: endTime       !time in sec when modeling ends
   integer(i8), save, public      :: nowTime       !current time in sec
   integer(i8), save, public      :: dateNodalJulian
   character(lc), save, public    :: fsbc_dir
   character(lc), save, public    :: fbdy_dir
   character(lc), save, public    :: fout_dir
   logical, save, public          :: restart_in = .false.
   logical, save, public          :: assim_in = .false.
   logical, save, public          :: restart_out = .false.
   logical, save, public          :: date_res = .true.
   real(wp), save, public         :: output_nmuber = 0.0_wp
   logical, save, public          :: ln_bdy = .false.
   logical, save, public          :: ln_tide = .false.
   logical, save, public          :: ln_bous = .false.
   integer, save, public          :: ln_pbt_base = 0

   integer, save, public          :: myIter = 0
   real(dp), save, public         :: globalArea
   real(wp), save, public         :: abFac
   logical, save, public          :: FirstMomU = .true.
   logical, save, public          :: FirstMomV = .true.
   real(wp), save, public         :: gammas = 0.0_wp
   real(wp), save, public         :: gammat = 0.0_wp
   integer, save, public          :: ffsfrs
   integer, save, public          :: ffsfrt
   character(lc), save, public    :: nffsfrs
   character(lc), save, public    :: nffsfrt
   logical, save, public          :: fsfrsavg = .true.
   logical, save, public          :: fsfrtavg = .true.
   logical, save, public          :: sbc_cli = .true.
   logical, save, public          :: sbc_cyc = .false.
   integer, save, public          :: ffuv
   integer, save, public          :: fftq
   integer, save, public          :: ffrad
   integer, save, public          :: ffprc
   integer, save, public          :: ffslp
   integer, save, public          :: ffrunoff
   character(lc), save, public    :: nffuv
   character(lc), save, public    :: nfftq
   character(lc), save, public    :: nffrad
   character(lc), save, public    :: nffprc
   character(lc), save, public    :: nffslp
   character(lc), save, public    :: nffrunoff
   logical, save, public          :: fuvavg
   logical, save, public          :: ftqavg
   logical, save, public          :: fradavg
   logical, save, public          :: fprcavg
   logical, save, public          :: fslpavg
   logical, save, public          :: frunoffavg
   logical, save, public          :: ln_qdew = .false.
   logical, save, public          :: ln_rain = .false.
   logical, save, public          :: ln_EmPmRevise = .false.
   integer(i8), save, public      :: nxtjultq
   integer(i8), save, public      :: nxtjuluv
   integer(i8), save, public      :: nxtjulrad
   integer(i8), save, public      :: nxtjulprc
   integer(i8), save, public      :: nxtjulslp
   integer(i8), save, public      :: nxtjulrunoff
   integer(i8), save, public      :: nxtjulsfrs
   integer(i8), save, public      :: nxtjulsfrt
   logical, save, public          :: rhoref
   integer, save, public          :: fbvec
   integer, save, public          :: fbt
   integer, save, public          :: fbs
   integer, save, public          :: fbpbt
   integer(i8), save, public      :: nxtjulVEC
   integer(i8), save, public      :: nxtjulT
   integer(i8), save, public      :: nxtjulS
   integer(i8), save, public      :: nxtjulPbt
   character(lc), save, public    :: nfbvec
   character(lc), save, public    :: nfbt
   character(lc), save, public    :: nfbs
   character(lc), save, public    :: nfbpbt
   logical, save, public          :: bvecavg
   logical, save, public          :: btavg
   logical, save, public          :: bsavg
   logical, save, public          :: bpbtavg
   integer, save, public          :: run_out_unit
   integer, save, public          :: sol_out_unit
   real(wp), save, public         :: zqt
   real(wp), save, public         :: zuv

   real(wp), save, public         :: ViscAhDCon
   real(wp), save, public         :: ViscAhZCon
   real(wp), save, public         :: ViscA4DCon
   real(wp), save, public         :: ViscA4ZCon
   logical, save, public          :: harmonic
   logical, save, public          :: biharmonic
   real(wp), save, public         :: BottomDrag
   real(wp), save, public         :: BottomDragMax
   real(wp), save, public         :: KappaRMCon
   real(wp), save, public         :: abEps
   real(wp), save, public         :: rStarFacLow
   real(wp), save, public         :: rStarFacUp
   real(wp), save, public         :: refMaxDep = 8000.0_wp
   integer, save, public          :: nkbar
   real(wp), save, public         :: rn_emin    ! minimum value of e [m^2/s^2]
   real(wp), save, public         :: rn_emin0   ! surface minimum value of tke [m^2/s^2]
   real(wp), save, public         :: rn_psimin
   real(wp), save, public         :: gravityRau0
   real(wp), save, public         :: gravityDynForce

   integer, save, public          :: ntracer
   integer, save, public          :: ntide = 0
   logical, save, public          :: tide2d = .true.
   integer, save, public          :: tideDim
   real(wp), save, public         :: rn_abs
   real(wp), save, public         :: rn_si0
   real(wp), save, public         :: KappaRTCon
   real(wp), save, public         :: KappaRTevd
   integer, save, public          :: nksr
   logical, save, public          :: harmonicT
   logical, save, public          :: biharmonicT
   real(wp), save, public         :: diffKhCon
   real(wp), save, public         :: diffK4Con

   real(wp), save, public         :: surf_out_ratio, mom_out_ratio, tracer_out_ratio
   real(wp), save, public         :: avoid0
   logical, save, public          :: l_zt_equal_zu = .false. ! if q and t are given at different height than U
   integer, save, public          :: myIterm

   !OPENACC CODE ----------------------------------------------------------------
   !$ACC DECLARE CREATE(                                                      &
   !$acc ndomain, idom, parent, nRatio, nl, nlpb, loc_nlpb, nlz, nlpbz,       &
   !$acc nk, nlbdy, bdywidth, loc_nlbdy, sbct, nkp1, ni, dTtracer, dTmom,     &
   !$acc dTsurf, nIterSurf, iterSurf_st, iterSurf_en, iterSurfs, nIter0,      &
   !$acc nIterMax, nOutFreq, fOutFreq, nResFreq, startyear, startmon,         &
   !$acc startday, starthour, startmin, startsec, nowyearf, nowmonf,          &
   !$acc nowdayf, nowhourf, nowminf, nowsecf, nowyear, nowmon, nowday,        &
   !$acc nowhour, nowmin, nowsec, nowjulian, endjulian, dateNodalJulian,      &
   !$acc restart_out, date_res, output_nmuber, ln_bdy, ln_tide, ln_bous,      &
   !$acc ln_pbt_base, myIter, globalArea, abFac, FirstMomU, FirstMomV,        &
   !$acc gammas, gammat, ffsfrs, ffsfrt, nffsfrs, nffsfrt,                    &
   !$acc fsfrsavg, fsfrtavg, sbc_cli, sbc_cyc, ffuv, fftq, ffrad,             &
   !$acc ffprc, ffslp, ffrunoff, nffuv, nfftq, nffrad, nffprc, nffslp,        &
   !$acc nffrunoff, fuvavg, ftqavg, fradavg, fprcavg, fslpavg, frunoffavg, ln_qdew, &
   !$acc nxtjultq, nxtjuluv, nxtjulrad, nxtjulprc, nxtjulslp, nxtjulrunoff,   &
   !$acc nxtjulsfrs, nxtjulsfrt, rhoref, fbvec, fbt, fbs, fbpbt, nxtjulVEC,   &
   !$acc nxtjulT, nxtjulS, nxtjulPbt, nfbvec, nfbt, nfbs, nfbpbt, bvecavg,    &
   !$acc btavg, bsavg, bpbtavg, zqt, zuv, ViscAhDCon, ViscAhZCon, ViscA4DCon, &
   !$acc ViscA4ZCon, harmonic, biharmonic, BottomDrag, BottomDragMax,         &
   !$acc KappaRMCon, abEps, rStarFacLow, rStarFacUp, refMaxDep, nkbar,        &
   !$acc rn_emin, rn_emin0, rn_psimin, gravityRau0, gravityDynForce, ntracer, &
   !$acc ntide, tide2d, tideDim, rn_abs, rn_si0, KappaRTCon, kappaRTevd, nksr, &
   !$acc harmonicT, biharmonicT, diffKhCon, diffK4Con, surf_out_ratio,        &
   !$acc mom_out_ratio, tracer_out_ratio, avoid0, l_zt_equal_zu)
   ! ----------------------------------------------------------------------------

   ! -- array for geometric
   real(wp), public, save, allocatable, target, dimension(:)       :: drF, drC  ! Unit is Pa
   !$ACC DECLARE CREATE(drF, drC)
   real(wp), public, save, allocatable, target, dimension(:)       :: rF, rC, rC_z  ! Unit is Pa
   !$ACC DECLARE CREATE(rF, rC, rC_z)
   !!YY, ZY, rF_z is not available in MaCOM now
   !!introduce it, or use rF/(rho*g) see mitice_init, and locate var swfracba for reference
   real(wp), public, save, allocatable, target, dimension(:)       :: rF_z  ! Unit is m
   !!YY,ZY: fcori at C-point, compute at mod_csp_init
   real(wp), public, save, allocatable, target, dimension(:)       :: fCoriG, fCori
   !$ACC DECLARE CREATE(fCoriG,fCori)
   real(wp), public, save, allocatable, target, dimension(:)       :: latC, lonC
   !$ACC DECLARE CREATE(latC, lonC)
   real(wp), public, save, allocatable, target, dimension(:)       :: dxC, dxZ, dxW, dxS
   !$ACC DECLARE CREATE(dxC, dxZ, dxW, dxS)
   real(wp), public, save, allocatable, target, dimension(:)       :: dyC, dyZ, dyW, dyS
   !$ACC DECLARE CREATE(dyC, dyZ, dyW, dyS)
   real(wp), public, save, allocatable, target, dimension(:)       :: rAc, rAs, rAw, rAz
   !$ACC DECLARE CREATE(rAc, rAs, rAw, rAz)
   integer, public, save, allocatable, dimension(:)        :: tw, te, ts, tn
   !$ACC DECLARE CREATE(tw, te, ts, tn)
   integer, public, save, allocatable, dimension(:)        :: tw2, te2, ts2, tn2
   !$ACC DECLARE CREATE(tw2, te2, ts2, tn2)
   integer, public, save, allocatable, dimension(:)        :: zw, ze, zs, zn
   !$ACC DECLARE CREATE(zw, ze, zs, zn)
   integer, public, save, allocatable, dimension(:, :)     :: uw, ue, us, un
   !$ACC DECLARE CREATE(uw, ue, us, un)
   integer, public, save, allocatable, dimension(:, :)     :: vw, ve, vs, vn
   !$ACC DECLARE CREATE(vw, ve, vs, vn)
   integer, public, save, allocatable, dimension(:, :)     :: z1, z2, z3, z4
   !$ACC DECLARE CREATE(z1, z2, z3, z4)
   integer, public, save, allocatable, dimension(:)        :: auz1, auz2, auz3, auz4, auz5, auz6
   !$ACC DECLARE CREATE(auz1, auz2, auz3, auz4, auz5, auz6)
   integer, public, save, allocatable, dimension(:)        :: avz1, avz2, avz3, avz4, avz5, avz6
   !$ACC DECLARE CREATE(avz1, avz2, avz3, avz4, avz5, avz6)
   integer, public, save, allocatable, dimension(:, :)     :: au1, au2, au3, au4
   !$ACC DECLARE CREATE(au1, au2, au3, au4)
   integer, public, save, allocatable, dimension(:, :)     :: av1, av2, av3, av4
   !$ACC DECLARE CREATE(av1, av2, av3, av4)
   integer(i2), public, save, allocatable, target, dimension(:)    :: kSurfC, kSurfW, kSurfS
   !$ACC DECLARE CREATE(kSurfC, kSurfW, kSurfS)
   integer(i2), public, save, allocatable, dimension(:)    :: kTopC, kTopW, kTopS
   !$ACC DECLARE CREATE(kTopC, kTopW, kTopS)
   real(wp), public, save, allocatable, dimension(:, :)    :: hFacC, hFacS, hFacW
   !$ACC DECLARE CREATE(hFacC, hFacS, hFacW)
   real(wp), public, save, allocatable, dimension(:)       :: hFacZ, recip_hFacZ
   !$ACC DECLARE CREATE(hFacZ, recip_hFacZ)
   real(wp), public, save, allocatable, target, dimension(:, :)    :: h0FacC, h0FacS, h0FacW
   !$ACC DECLARE CREATE(h0FacC, h0FacS, h0FacW)
   integer(i2), public, save, allocatable, target, dimension(:, :) :: maskC, maskS
   !$ACC DECLARE CREATE(maskC, maskS)
   integer(i2), public, save, allocatable, target, dimension(:, :) :: maskW, maskZ
   !$ACC DECLARE CREATE(maskW, maskZ)
   integer(i2), public, save, allocatable, dimension(:)    :: maskBdy
   !$ACC DECLARE CREATE(maskBdy)
   integer(i2), public, save, allocatable, dimension(:)    :: normal_idx
   !$ACC DECLARE CREATE(normal_idx)
   real(wp), public, save, allocatable, dimension(:)       :: recip_drF, recip_drC
   !$ACC DECLARE CREATE(recip_drF, recip_drC)
   real(wp), public, save, allocatable, dimension(:)       :: recip_dxC, recip_dxZ
   !$ACC DECLARE CREATE(recip_dxC, recip_dxZ)
   real(wp), public, save, allocatable, dimension(:)       :: recip_dyC, recip_dyZ
   !$ACC DECLARE CREATE(recip_dyC, recip_dyZ)
   real(wp), public, save, allocatable, dimension(:)       :: recip_dxW, recip_dxS
   !$ACC DECLARE CREATE(recip_dxW, recip_dxS)
   real(wp), public, save, allocatable, dimension(:)       :: recip_dyW, recip_dyS
   !$ACC DECLARE CREATE(recip_dyW, recip_dyS)
   real(wp), public, save, allocatable, dimension(:)       :: recip_rAc, recip_rAs
   !$ACC DECLARE CREATE(recip_rAc, recip_rAs)
   real(wp), public, save, allocatable, dimension(:)       :: recip_rAw, recip_rAz
   !$ACC DECLARE CREATE(recip_rAw, recip_rAz)
   real(wp), public, save, allocatable, dimension(:)       :: recip_RcolC, recip_RcolZ
   !$ACC DECLARE CREATE(recip_RcolC, recip_RcolZ)
   real(wp), public, save, allocatable, dimension(:)       :: recip_RcolW, recip_RcolS
   !$ACC DECLARE CREATE(recip_RcolW, recip_RcolS)
   real(wp), public, save, allocatable, dimension(:)       :: RcolC, RcolW, RcolS, RcolZ
   !$ACC DECLARE CREATE(RcolC, RcolW, RcolS, RcolZ)
   real(wp), public, save, allocatable, dimension(:, :)    :: recip_hFacC, recip_hFacS
   !$ACC DECLARE CREATE(recip_hFacC, recip_hFacS)
   real(wp), public, save, allocatable, dimension(:, :)    :: recip_hFacW
   !$ACC DECLARE CREATE(recip_hFacW)

   ! -- array for sea surface force
   real(wp), public, save, allocatable, dimension(:)       :: taum, utau, vtau
   !$ACC DECLARE CREATE(taum, utau, vtau)
   real(wp), public, save, allocatable, dimension(:)       :: utau_out, vtau_out
   !$ACC DECLARE CREATE(utau_out, vtau_out)
   real(wp), public, save, allocatable, dimension(:)       :: guExt, gvExt
   !$ACC DECLARE CREATE(guExt, gvExt)
   real(wp), public, save, allocatable, dimension(:)       :: phi0surf
   !$ACC DECLARE CREATE(phi0surf)
   real(wp), public, save, allocatable, dimension(:)       :: EmPmR, EmPmR_long, EmPmR_out ! Unit is kg/(m^2*s)
   !$ACC DECLARE CREATE(EmPmR, EmPmR_long, EmPmR_out)
   real(wp), public, save, allocatable, dimension(:, :)    :: addMass ! Unit is kg/s
   !$ACC DECLARE CREATE(addMass)

   ! -- array for time forward
   ! -- array for dynamics
   ! array for baroclinic-pressure-velocity coupling
   real(wp), public, save, allocatable, dimension(:)       :: omega3
   !$ACC DECLARE CREATE(omega3)
   real(wp), public, save, allocatable, dimension(:)       :: vort3
   !$ACC DECLARE CREATE(vort3)
   real(wp), public, save, allocatable, dimension(:)       :: hDiv
   !$ACC DECLARE CREATE(hDiv)
   real(wp), public, save, allocatable, dimension(:)       :: KE
   !$ACC DECLARE CREATE(KE)
   real(wp), public, save, allocatable, dimension(:)       :: dKEdx, dKEdy
   !$ACC DECLARE CREATE(dKEdx, dKEdy)
   real(wp), public, save, allocatable, dimension(:)       :: uCorTerm, vCorTerm
   !$ACC DECLARE CREATE(uCorTerm, vCorTerm)
   real(wp), public, save, allocatable, dimension(:)       :: phiHydC, phiHydF
   !$ACC DECLARE CREATE(phiHydC, phiHydF)
   real(wp), public, save, allocatable, dimension(:)       :: dPhiHydX, dPhiHydY
   !$ACC DECLARE CREATE(dPhiHydX, dPhiHydY)
   real(wp), public, save, allocatable, target, dimension(:)       :: rhoLev
   !$ACC DECLARE CREATE(rhoLev)
   real(wp), public, save, allocatable, dimension(:)       :: uShearTerm, vShearTerm
   !$ACC DECLARE CREATE(uShearTerm, vShearTerm)
   real(wp), public, save, allocatable, target, dimension(:, :)    :: uFld, vFld, wFld
   !$ACC DECLARE CREATE(uFld, vFld, wFld)
   real(wp), public, save, allocatable, dimension(:, :)    :: uFldBcl, vFldBcl
   !$ACC DECLARE CREATE(uFldBcl, vFldBcl)
   real(wp), public, save, allocatable, dimension(:, :)    :: uFldBcl_past, vFldBcl_past
   !$ACC DECLARE CREATE(uFldBcl_past, vFldBcl_past)
   real(wp), public, save, allocatable, dimension(:, :)    :: uFldBclDash, vFldBclDash
   !$ACC DECLARE CREATE(uFldBclDash, vFldBclDash)
   real(wp), public, save, allocatable, dimension(:, :)    :: uFld_out, vFld_out, wFld_out
   !$ACC DECLARE CREATE(uFld_out, vFld_out, wFld_out)
   real(wp), public, save, allocatable, dimension(:, :)    :: gUnow, gVnow
   !$ACC DECLARE CREATE(gUnow, gVnow)
   real(wp), public, save, allocatable, dimension(:, :)    :: gUpast, gVpast
   !$ACC DECLARE CREATE(gUpast, gVpast)
   real(wp), public, save, allocatable, dimension(:)       :: uBar, vBar
   !$ACC DECLARE CREATE(uBar, vBar)
   real(wp), public, save, allocatable, dimension(:)       :: uBarPast, vBarPast
   !$ACC DECLARE CREATE(uBarPast, vBarPast)
   real(wp), public, save, allocatable, dimension(:)       :: uBarCouple, vBarCouple
   !$ACC DECLARE CREATE(uBarCouple, vBarCouple)
   real(wp), public, save, allocatable, dimension(:)       :: gUBarnow, gVBarnow
   !$ACC DECLARE CREATE(gUBarnow, gVBarnow)
   real(wp), public, save, allocatable, dimension(:)       :: gUBarpast, gVBarpast
   !$ACC DECLARE CREATE(gUBarpast, gVBarpast)

   ! array for barotropic-pressure-velocity coupling
   real(wp), public, save, allocatable, dimension(:)       :: etaN, etaH, etaH_out, pbt_base  ! Unit is Pa
   !$ACC DECLARE CREATE(etaN, etaH, etaH_out, pbt_base)
   real(wp), public, save, allocatable, dimension(:)       :: ssh, ssh_out ! Unit is m
   !$ACC DECLARE CREATE(ssh, ssh_out)
   real(wp), public, save, allocatable, dimension(:)       :: dphiSurfX, dphiSurfY
   !$ACC DECLARE CREATE(dphiSurfX, dphiSurfY)
   real(wp), public, save, allocatable, dimension(:)       :: Bo_surf, rBo_surf
   !$ACC DECLARE CREATE(Bo_surf, rBo_surf)
   real(wp), public, save, allocatable, dimension(:)       :: rSurfC, rSurfW, rSurfS
   !$ACC DECLARE CREATE(rSurfC, rSurfW, rSurfS)
   real(wp), public, save, allocatable, dimension(:)       :: rLowC, rLowW, rLowS
   !$ACC DECLARE CREATE(rLowC, rLowW, rLowS)
   real(wp), public, save, allocatable, dimension(:)       :: rSurfC_min
   !$ACC DECLARE CREATE(rSurfC_min)
   real(wp), public, save, allocatable, dimension(:)       :: rStarFacC, rStarFacW, rStarFacS
   !$ACC DECLARE CREATE(rStarFacC, rStarFacW, rStarFacS)
   real(wp), public, save, allocatable, dimension(:)       :: rStarExpC
   !$ACC DECLARE CREATE(rStarExpC)

   ! array for horizontal momentum dissipation
   real(wp), public, save, allocatable, dimension(:)       :: del2u, del2v
   !$ACC DECLARE CREATE(del2u, del2v)
   real(wp), public, save, allocatable, dimension(:)       :: dStar, zStar
   !$ACC DECLARE CREATE(dStar, zStar)
   real(wp), public, save, allocatable, dimension(:)       :: uDissip, vDissip
   !$ACC DECLARE CREATE(uDissip, vDissip)
   real(wp), public, save, allocatable, dimension(:)       :: viscAh_D, viscAh_Z
   !$ACC DECLARE CREATE(viscAh_D, viscAh_Z)
   real(wp), public, save, allocatable, dimension(:)       :: viscA4_D, viscA4_Z
   !$ACC DECLARE CREATE(viscA4_D, viscA4_Z)
   real(wp), public, save, allocatable, dimension(:)       :: lsmagt2, lsmagf2
   !$ACC DECLARE CREATE(lsmagt2, lsmagf2)
   real(wp), public, save, allocatable, dimension(:, :)    :: KappaRM, KappaRM_out
   !$ACC DECLARE CREATE(KappaRM, KappaRM_out)

   ! array for momentum drag
   ! SideDragFac: free slip = 0, partial slip 0~2 (most time is 1), no slip = 2
   ! SideDragFac is an array, in the internal ocean which should always be 1,
   ! the value only at lane-sea boundary has meaning.
   real(wp), public, save, allocatable, dimension(:, :)    :: SideDragFac
   !$ACC DECLARE CREATE(SideDragFac)
   real(wp), public, save, allocatable, dimension(:)       :: BottomDragFac
   !$ACC DECLARE CREATE(BottomDragFac)
   real(wp), public, save, allocatable, dimension(:)       :: uDragTerms, vDragTerms
   !$ACC DECLARE CREATE(uDragTerms, vDragTerms)
   real(wp), public, save, allocatable, dimension(:)       :: uDragBottom, vDragBottom
   !$ACC DECLARE CREATE(uDragBottom, vDragBottom)

   ! -- array for thermodynamic
   ! array for tracer
   real(wp), public, save, allocatable, dimension(:)       :: alphaRho
   !$ACC DECLARE CREATE(alphaRho)
   real(wp), public, save, allocatable, dimension(:, :)    :: rhoInSitu
   !$ACC DECLARE CREATE(rhoInSitu)
   real(wp), public, save, allocatable, target, dimension(:, :, :) :: tFld, tFld_out
   !$ACC DECLARE CREATE(tFld, tFld_out)
   real(wp), public, save, allocatable, dimension(:, :)    :: gTracer
   !$ACC DECLARE CREATE(gTracer)
   real(wp), public, save, allocatable, dimension(:, :, :) :: gTracerPast
   !$ACC DECLARE CREATE(gTracerPast)
   real(wp), public, save, allocatable, dimension(:, :)    :: tracer
   !$ACC DECLARE CREATE(tracer)
   logical, public, save, allocatable, dimension(:)        :: FirstTracer
   !$ACC DECLARE CREATE(FirstTracer)
   real(wp), public, save, allocatable, dimension(:, :)    :: diffKhu, diffKhv
   !$ACC DECLARE CREATE(diffKhu, diffKhv)
   real(wp), public, save, allocatable, dimension(:, :)    :: diffK4u, diffK4v
   !$ACC DECLARE CREATE(diffK4u, diffK4v)
   ! array for horizontal tracer diffusion

   ! array for vertical tracer diffusion
   real(wp), public, save, allocatable, dimension(:, :)    :: rn2
   !$ACC DECLARE CREATE(rn2)
   real(wp), public, save, allocatable, dimension(:, :)    :: KappaRT, KappaRS
   !$ACC DECLARE CREATE(KappaRT, KappaRS)
   real(wp), public, save, allocatable, dimension(:, :)    :: KappaRT_out
   !$ACC DECLARE CREATE(KappaRT_out)
   real(wp), public, save, allocatable, dimension(:, :)    :: KappaRX
   !$ACC DECLARE CREATE(KappaRX)
   real(wp), public, save, allocatable, dimension(:, :)    :: hmxl   ! mixing length [Pa]
   !$ACC DECLARE CREATE(hmxl)
   real(wp), public, save, allocatable, dimension(:, :)    :: en     ! turbulent kinetic energy [Pa^2/s^2]
   !$ACC DECLARE CREATE(en)

   ! -- array for sea surface force of tracer
   real(wp), public, save, allocatable, dimension(:, :)    :: gtExt
   !$ACC DECLARE CREATE(gtExt)
   real(wp), public, save, allocatable, dimension(:)       :: qsr, qsr_out
   !$ACC DECLARE CREATE(qsr, qsr_out)
   real(wp), public, save, allocatable, dimension(:)       :: qns, qns_out
   !$ACC DECLARE CREATE(qns, qns_out)

   ! -- ZY, YY:definition refer to mod_csp_forcing, positive corresponding to upward 
   real(wp), public, save, allocatable, dimension(:)       :: Qnet
   real(wp), public, save, allocatable, dimension(:)       :: Qsw
   real(wp), public, save, allocatable, dimension(:)       :: saltflux
   !$ACC DECLARE CREATE(Qnet,Qsw,saltflux)


   ! -- array for sea surface force data
   real(wp), public, save, allocatable, dimension(:, :)    :: u10, v10
   !$ACC DECLARE CREATE(u10, v10)
   real(wp), public, save, allocatable, dimension(:, :)    :: t10, q10
   !$ACC DECLARE CREATE(t10, q10)
   real(wp), public, save, allocatable, dimension(:, :)    :: lwdn, swdn
   !$ACC DECLARE CREATE(lwdn, swdn)
   real(wp), public, save, allocatable, dimension(:, :)    :: prec, snow
   !$ACC DECLARE CREATE(prec, snow)
   real(wp), public, save, allocatable, dimension(:, :)    :: slp, runoff
   !$ACC DECLARE CREATE(slp, runoff)
   real(wp), public, save, allocatable, dimension(:, :)    :: sfrt, sfrs
   !$ACC DECLARE CREATE(sfrt, sfrs)
   !YY,ZY: add zevap, which was a local var in mod_csp_force
   real(wp), public, save, allocatable, dimension(:)    :: zevap
   !$ACC DECLARE CREATE(zevap)

   ! -- array for sea boundary data
   real(wp), public, save, allocatable, dimension(:)       :: bdy_uBar, bdy_vBar
   !$ACC DECLARE CREATE(bdy_uBar, bdy_vBar)
   real(wp), public, save, allocatable, dimension(:, :)    :: bdy_uBcl, bdy_vBcl
   !$ACC DECLARE CREATE(bdy_uBcl, bdy_vBcl)
   real(wp), public, save, allocatable, dimension(:, :)    :: bdy_pbt
   !$ACC DECLARE CREATE(bdy_pbt)
   real(wp), public, save, allocatable, dimension(:, :, :) :: bdy_t, bdy_s, bdy_u, bdy_v, bdy_w
   !$ACC DECLARE CREATE(bdy_t, bdy_s, bdy_u, bdy_v, bdy_w)

   ! array for tide
   integer, public, save, allocatable, dimension(:)        :: tide_idx
   character(len=4), public, save, allocatable, dimension(:) :: tide_name_use
   real(wp), public, save, allocatable, dimension(:)       :: tide_freq_use
   !$ACC DECLARE CREATE(tide_freq_use)
   real(wp), public, save, allocatable, dimension(:,:)     :: tide_amp_ssh, tide_pha_ssh
   !$ACC DECLARE CREATE(tide_amp_ssh, tide_pha_ssh)
   real(wp), public, save, allocatable, dimension(:,:)     :: tide_amp_u, tide_pha_u
   !$ACC DECLARE CREATE(tide_amp_u, tide_pha_u)
   real(wp), public, save, allocatable, dimension(:,:)     :: tide_amp_v, tide_pha_v
   !$ACC DECLARE CREATE(tide_amp_v, tide_pha_v)
   real(wp), public, save, allocatable, dimension(:,:)     :: tide_u_ssh, tide_v_ssh, tide_f_ssh
   !$ACC DECLARE CREATE(tide_u_ssh, tide_v_ssh, tide_f_ssh)
   real(wp), public, save, allocatable, dimension(:,:)     :: tide_u_curu, tide_v_curu, tide_f_curu
   !$ACC DECLARE CREATE(tide_u_curu, tide_v_curu, tide_f_curu)
   real(wp), public, save, allocatable, dimension(:,:)     :: tide_u_curv, tide_v_curv, tide_f_curv
   !$ACC DECLARE CREATE(tide_u_curv, tide_v_curv, tide_f_curv)
   real(wp), public, save, allocatable, dimension(:)       :: tide_ssh, tide_u, tide_v
   !$ACC DECLARE CREATE(tide_ssh, tide_u, tide_v)

   ! array for long time tracer diagnosis
   real(wp), public, save, allocatable, dimension(:, :)    :: tFld_adv, sFld_adv, tFld_vdiff, sFld_vdiff, &
                                                              tFld_hdiff, sFld_hdiff, tFld_force, sFld_force
   !$ACC DECLARE CREATE(tFld_adv, sFld_adv, tFld_vdiff, sFld_vdiff, tFld_hdiff, sFld_hdiff, tFld_force, sFld_force) 


   type, public :: nc_attr
      character(lc)  :: var_units = 'none'
      character(lc)  :: var_longname = 'none'
   end type nc_attr

end module mod_csp_basic
