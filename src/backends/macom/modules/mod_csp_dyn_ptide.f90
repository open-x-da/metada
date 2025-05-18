!******************************************************************************
! Copyright          : 2022
! File Name          : mod_csp_dyn_ptide.f90
! Description        : This program will calculate tide potential, and add it in
!   barotropic process (csp_dyn_pbt_exp in mod_csp_dyn_pbt_exp.f90) as . Method
!   used in this program can reference "B.K. Arbic et al., 2004. The accuracy of
!   surface elevations in forward global barotropic and baroclinic tide models".
!   First version is provided by Prof. Li Wei (TJU).
!
! Revision History   :
! Date            Author         Comments
!******************************************************************************
! 2022.10.12      Wei Li         First version provided by Li Wei
! 2022.10.14      Yu Zhang       Modified for MaCOM
!******************************************************************************

module mod_csp_dyn_ptide
   use mod_misc
   use mod_mpi_test
   implicit none

   integer, parameter :: NTIDE0=2, NTIDE1=4, NTIDE2=4
   real(wp) :: THOUR
   real(wp) :: WU0(NTIDE0), WV0(NTIDE0), WF0(NTIDE0),  &
               WU1(NTIDE1), WV1(NTIDE1), WF1(NTIDE1),  WU2(NTIDE2), &
               WV2(NTIDE2), WF2(NTIDE2)
   integer, save :: IU45(2,13)
   DATA IU45/-2,-1, &
             -2, 0, &
             -2, 1, &
              0,-2, &
              0,-1, &
              0, 0, &
              0, 1, &
              0, 2, &
              0, 3, &
              2,-1, &
              2, 0, &
              2, 1, &
              2, 2/
   real(wp), save :: W0(NTIDE0), AMP0(NTIDE0), W1(NTIDE1), AMP1(NTIDE1), &
                     W2(NTIDE2), AMP2(NTIDE2)
   real(wp), save :: COEFF(13,10)
   DATA COEFF/ 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0008_wp,-0.0657_wp, 1.0000_wp, &
              -0.0649_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0534_wp,-0.0218_wp,-0.0059_wp, &
              -0.0023_wp, 0.0432_wp,-0.0028_wp, 0.0000_wp, 0.0000_wp, 1.0000_wp, &
               0.4143_wp, 0.0387_wp,-0.0008_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0058_wp, 0.1885_wp, 1.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0002_wp,-0.0064_wp,-0.0010_wp, 0.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0008_wp,-0.0112_wp, 1.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0015_wp,-0.0003_wp, 0.0000_wp, &
               0.0002_wp, 0.0000_wp, 0.0000_wp, 0.0001_wp,-0.0198_wp, 1.0000_wp, &
               0.1356_wp,-0.0029_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0294_wp, 1.0000_wp, &
               0.1980_wp,-0.0047_wp, 0.0000_wp, 0.0000_wp,-0.0152_wp,-0.0098_wp,-0.0057_wp, &
              -0.0037_wp, 0.1496_wp, 0.0296_wp, 0.0000_wp, 0.0000_wp, 1.0000_wp, &
               0.6398_wp, 0.1342_wp, 0.0086_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0005_wp,-0.0373_wp, 1.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0006_wp, 0.0002_wp, 0.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0366_wp, 1.0000_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0047_wp,-0.2505_wp,-0.1102_wp,-0.0156_wp, &
               0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0128_wp, 1.0000_wp, &
               0.2980_wp, 0.0324_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp/
   !$ACC declare create(WU0, WV0, WF0, WU1, WV1, WF1, WU2, WV2, WF2, IU45, W0, &
   !$ACC                AMP0, W1, AMP1, W2, AMP2, COEFF, THOUR)

   real(wp), save, allocatable, dimension(:,:) :: ELM0, ELM1, ELM2
   !$acc declare create(ELM0, ELM1, ELM2)

   real(wp), public, save, allocatable, dimension(:) :: BAREL
   !$acc declare create(BAREL)
!
contains

!==============================================================================
   subroutine csp_dyn_ptide_init
!==============================================================================
      implicit none
      integer :: i, k

      allocate(ELM0(nlpb, NTIDE0), ELM1(nlpb, NTIDE1), ELM2(nlpb, NTIDE2))
      ELM0 = 0.0_wp
      ELM1 = 0.0_wp
      ELM2 = 0.0_wp

      allocate(BAREL(nlpb))
      BAREL = 0.0_wp

      !	Mm,Mf
      !	K1,O1,P1,Q1
      !	M2,S2,N2,K2
      W0(1)=0.00950112_wp
      W0(2)=0.01916424_wp
      W1(1)=0.26251617_wp
      W1(2)=0.24335188_wp
      W1(3)=0.26108260_wp
      W1(4)=0.23385075_wp
      W2(1)=0.50586805_wp
      W2(2)=0.52359878_wp
      W2(3)=0.49636692_wp
      W2(4)=0.52503234_wp
      AMP0(1)=0.015_wp
      AMP0(2)=0.029_wp
      AMP1(1)=0.104_wp     ! 0.104_wp Arbic_etal_DSR_2004 or 0.098_wp Ye Anle ?
      AMP1(2)=0.070_wp
      AMP1(3)=0.033_wp     ! 0.033_wp Arbic_etal_DSR_2004 or 0.032_wp Ye Anle ?
      AMP1(4)=0.013_wp
      AMP2(1)=0.168_wp
      AMP2(2)=0.078_wp
      AMP2(3)=0.032_wp
      AMP2(4)=0.021_wp

      !$ACC UPDATE DEVICE(BAREL, IU45, W0, AMP0, W1, AMP1, W2, AMP2, COEFF, ELM0, ELM1, ELM2)
      !$ACC kernels present(latC, AMP0, AMP1, AMP2, ELM0, ELM1, ELM2)
      !$ACC loop independent
      do i = 1, nlpb
          do k = 1, NTIDE0
            ELM0(i, k) = AMP0(k)*(0.5_wp-1.5_wp*sin(latC(i)*radPi)*sin(latC(i)*radPi))
          end do
          do k = 1, NTIDE1
            ELM1(i, k) = AMP1(k)*sin(2.0_wp*latC(i)*radPi)
          end do
          do k = 1, NTIDE2
            ELM2(i, k) = AMP2(k)*cos(latC(i)*radPi)*cos(latC(i)*radPi)
          end do
      end do
      !$ACC end kernels

   end subroutine csp_dyn_ptide_init

!==============================================================================
   SUBROUTINE csp_dyn_ptide
!==============================================================================
      implicit none
      integer :: I, KK
      real(wp) :: YEAR, FTIME, DHOUR, RLON
      integer(8) :: juliantemp

      call greg2jul(0, 0, 0, 1, 1, nowyear, juliantemp)
      juliantemp = nowjulian - juliantemp
      THOUR = real(juliantemp, wp)/3600.0_wp
      FTIME = THOUR/24.0_wp
      DHOUR = dTtracer/3600.0_wp
      YEAR = FLOAT(nowyear)

      CALL FUV0(WV0, YEAR, 0)
      CALL FUV1(WV1, YEAR, 0)
      CALL FUV2(WV2, YEAR, 0)
      CALL FU10_0(WF0, WU0, YEAR, FTIME)
      CALL FU10_1(WF1, WU1, YEAR, FTIME)
      CALL FU10_2(WF2, WU2, YEAR, FTIME)

      !$ACC UPDATE DEVICE(THOUR, WF0, W0, WU0, WV0, WF1, W1, WU1, WU1, WV1, &
      !$ACC               WF2, W2, WU2, WV2)
      !$ACC data present(lonC, BAREL, maskC, ELM0, WF0, W0, WU0, WV0, ELM1, WF1, &
      !$ACC              W1, WU1, WU1, WV1, ELM2, WF2, W2, WU2, WV2)
      !$ACC kernels
      !$ACC loop independent
      DO I = 1, nlpb
         RLON = lonC(I)*radPi

         BAREL(I) = 0.0_wp

         DO KK = 1, NTIDE0
            BAREL(I) = BAREL(I)+ELM0(I,KK)*WF0(KK)*cos(W0(KK)*THOUR+WU0(KK)+WV0(KK))  !-8.0_wp
         END DO

         DO KK = 1, NTIDE1
            BAREL(I) = BAREL(I)+ELM1(I,KK)*WF1(KK)*cos(W1(KK)*THOUR+RLON+WU1(KK)+WV1(KK))  !-8.0_wp
         END DO

         DO KK = 1, NTIDE2
           BAREL(I) = BAREL(I)+ELM2(I,KK)*WF2(KK)*cos(W2(KK)*THOUR+2.0_wp*RLON+WU2(KK)+WV2(KK))  !-8.0_wp
         END DO

         BAREL(I) = BAREL(I)*gravityRau0*maskC(I,nk)

      END DO
      !$ACC end kernels

      !$ACC end data

      RETURN

   END SUBROUTINE csp_dyn_ptide

!==============================================================================
   SUBROUTINE FUV0(UU,Y,ID)
!==============================================================================
      IMPLICIT NONE
      real(wp) :: UU(2)
      real(wp) :: S0,H0,P0
      real(wp) :: Y
      INTEGER :: I,ID,IYC

      IYC=INT((Y-1901.0_wp)/4.0_wp)
      S0=277.025_wp+129.3848_wp*(Y-1900.0_wp)+13.17640_wp*(ID+IYC)
      H0=280.190_wp-0.23872_wp*(Y-1900._wp)+0.98565_wp*(ID+IYC)
      P0=334.385_wp+40.66249_wp*(Y-1900._wp)+0.11140_wp*(ID+IYC)
      UU(1)=(S0-P0)*atan(1.0_wp)/45.0_wp
      UU(2)=2.0_wp*S0*atan(1.0_wp)/45.0_wp

      RETURN
   END SUBROUTINE FUV0

!==============================================================================
   SUBROUTINE FUV1(UU,Y,ID)
!==============================================================================
      IMPLICIT NONE
      real(wp) :: UU(4)
      real(wp) :: S0,H0,P0
      real(wp) :: Y
      INTEGER :: I,ID,IYC

      IYC=INT((Y-1901.0_wp)/4.0_wp)
      S0=277.025_wp+129.3848_wp*(Y-1900.0_wp)+13.17640_wp*(ID+IYC)
      H0=280.190_wp-0.23872_wp*(Y-1900.0_wp)+0.98565_wp*(ID+IYC)
      P0=334.385_wp+40.66249_wp*(Y-1900.0_wp)+0.11140_wp*(ID+IYC)
      UU(1)=(H0+90.0_wp)*atan(1.0_wp)/45.0_wp
      UU(2)=(-2.0_wp*S0+H0+270.0_wp)*atan(1.0_wp)/45.0_wp
      UU(3)=(-H0+270.0_wp)*atan(1.0_wp)/45.0_wp
      UU(4)=(-3.0_wp*S0+H0+P0+270.0_wp)*atan(1.0_wp)/45.0_wp

      RETURN
   END SUBROUTINE FUV1

!==============================================================================
   SUBROUTINE FUV2(UU,Y,ID)
!==============================================================================
      IMPLICIT NONE
      real(wp) :: UU(4)
      real(wp) :: S0,H0,P0
      real(wp) :: Y
      INTEGER :: I,ID,IYC

      IYC=INT((Y-1901.0_wp)/4.0_wp)
      S0=277.025_wp+129.3848_wp*(Y-1900.0_wp)+13.17640_wp*(ID+IYC)
      H0=280.190_wp-0.23872_wp*(Y-1900.0_wp)+0.98565_wp*(ID+IYC)
      P0=334.385_wp+40.66249_wp*(Y-1900.0_wp)+0.11140_wp*(ID+IYC)
      UU(1)=(-2.0_wp*S0+2.0_wp*H0)*atan(1.0_wp)/45.0_wp
      UU(2)=0.0_wp*atan(1.0_wp)/45.0_wp
      UU(3)=(-3.0_wp*S0+2.0_wp*H0+P0)*atan(1.0_wp)/45.0_wp
      UU(4)=2.0_wp*H0*atan(1.0_wp)/45.0_wp

      RETURN
   END SUBROUTINE FUV2

!==============================================================================
   SUBROUTINE FU10_0(F2,U2,YY,DD)
!==============================================================================
      IMPLICIT NONE
      real(wp) :: YY,DD
      real(wp) :: FF(11),UU(11)
      real(wp) :: F2(2),U2(2)

      CALL FU(FF, UU, YY, DD)
      F2(1)=FF(1)
      U2(1)=UU(1)
      F2(2)=FF(2)
      U2(2)=UU(2)

      RETURN
   END SUBROUTINE FU10_0

!==============================================================================
   SUBROUTINE FU10_1(F4,U4,YY,DD)
!==============================================================================
      IMPLICIT NONE
      real(wp) :: YY,DD
      real(wp) :: FF(11),UU(11)
      real(wp) :: F4(4),U4(4)

      CALL FU(FF, UU, YY, DD)
      F4(1)=FF(5)
      U4(1)=UU(5)
      F4(2)=FF(3)
      U4(2)=UU(3)
      F4(3)=FF(4)
      U4(3)=UU(4)
      F4(4)=FF(3)
      U4(4)=UU(3)
      RETURN
   END SUBROUTINE FU10_1

!==============================================================================
   SUBROUTINE FU10_2(F4,U4,YY,DD)
!==============================================================================
      IMPLICIT NONE
      real(wp) :: YY,DD
      real(wp) :: FF(11),UU(11)
      real(wp) :: F4(4),U4(4)

      CALL FU(FF, UU, YY, DD)
      F4(1)=FF(8)
      U4(1)=UU(8)
      F4(2)=1.0_wp
      U4(2)=0.0_wp
      F4(3)=FF(8)
      U4(3)=UU(8)
      F4(4)=FF(10)
      U4(4)=UU(10)
      RETURN
   END SUBROUTINE FU10_2

!==============================================================================
   SUBROUTINE FU(FF,UU,YY,DD)
!==============================================================================
      IMPLICIT NONE
   ! --------------------
      real(wp) :: FF(11),UU(11),YY,DD
   ! -----
      INTEGER :: I,J,K,IYC
      real(wp) :: P0,PN0,FCOSU,FSINU

   ! --------------------
      IYC=(YY-1901.0_wp)/4.0_wp
      P0=334.385_wp+40.66249_wp*(YY-1900)+0.11140_wp*(DD+IYC)
      PN0=100.84_wp+19.3282_wp*(YY-1900)+0.0530_wp*(DD+IYC)
      P0=P0*radPi
      PN0=PN0*radPi
      DO K=1,11
         FCOSU=0.0_wp
         FSINU=0.0_wp
         IF(K.LE.10)THEN
            DO I=1,13
               FCOSU=FCOSU+COEFF(I,K)*cos(IU45(1,I)*P0+IU45(2,I)*PN0)
               FSINU=FSINU+COEFF(I,K)*sin(IU45(1,I)*P0+IU45(2,I)*PN0)
            END DO
         END IF
         IF(K.EQ.11)THEN
            FCOSU = - 0.008_wp*cos(-P0-2.0_wp*PN0) + 0.094_wp*cos(-P0-PN0) &
                    + 0.510_wp*cos(-P0) - 0.041_wp*cos(P0-PN0) &
                    + 1.418_wp*cos(P0) + 0.284_wp*cos(P0+PN0) - 0.008_wp*cos(P0+2.0_wp*PN0)
            FSINU = - 0.008_wp*sin(-P0-2.0_wp*PN0) + 0.094_wp*sin(-P0-PN0) &
                    + 0.510_wp*sin(-P0) - 0.041_wp*sin(P0-PN0) &
                    + 1.418_wp*sin(P0) + 0.284_wp*sin(P0+PN0) - 0.008_wp*sin(P0+2.0_wp*PN0)
         END IF
         FF(K)=sqrt(FCOSU*FCOSU+FSINU*FSINU)
         IF(FCOSU.EQ.0.0_wp) THEN
            IF(FSINU.GT.0.0_wp)THEN
               UU(K)=90.0_wp
            ELSE
               UU(K)=270.0_wp
            END IF
         ELSE
            UU(K)=atan2(FSINU,FCOSU)/radPi
            IF(K.EQ.11.AND.UU(K).LT.0.0_wp)THEN
               UU(K)=UU(K)+360.0_wp
            END IF
         END IF
      END DO
      DO I=1,11
         UU(I)=UU(I)*atan(1.0_wp)/45.0_wp
      END DO
      RETURN

   END SUBROUTINE FU

end module mod_csp_dyn_ptide
