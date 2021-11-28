module NumOfOrbitalAndElectrons
   integer*8 :: N3d, N2p, N, Ndim
   integer*8 :: nup, ndn, ntol
endmodule NumOfOrbitalAndElectrons

module ModelParas
   double precision :: U, E_site(1:30)
   double precision :: Uc, t, tt, ttt, miu
   double precision :: Ud, Up, UQ, E_d, E_pxy, E_pz, tpd, tpp, tpdz, tpzdz, tppz
   double precision :: U_eg, U_t2g, Upp, E_eg, E_t2g, cfs
   double precision :: rA, rB, rC, F_0, F_2, G_1, G_3, Dq, Ds, Dt
   double precision :: U_rest(5,5,5,5), U_pddp(3,5,5,3), U_dpdp(5,3,5,3)
   double precision :: U_ext(5,5), J_ext(5,5)
   double precision :: phase(30,30)
   double complex :: U_SO_2p(3,2,3,2)
   double precision :: xi_2p, kx, ky
   double precision :: disp(0:15), lsf(0:15), lambda
   double complex :: prefac(4) ! for twisted average boundary condition
   integer*8 :: nbx, nby, ibx, iby ! number of boundaries in x and y direction 
endmodule ModelParas

module MPIParas
   integer  ::   comm, myid, nprocs
   integer  ::   source, dest, tag, ierr, rc
   integer*8 :: nloc, localstart, localend
endmodule MPIParas

module ConstantParas
   integer*8 :: HsizeEst=72000, SprSizeEst=72000*200
   integer*8 :: Hsize_3dEst=72000, SprSize_3dEst=72000*200
   integer*8 :: NNN=1000, niter_CL=100, niter_CFE = 100, niter_CG=100
   double precision :: tol_CG=0.0001, tol_CL=0.0001
   double precision :: pi=3.1415926535, sqrthalf = 0.707106781
   double precision :: epsilone_CFE=0.3d0,  epsilone_CG = 1.0d0
endmodule ConstantParas

module BettsCluster
   character*1 :: Nmore
   integer*8 :: rpoint(0:15,2), runit
   integer*8 :: kpoint(0:15,2), kunit
   integer*8 :: site(0:15,4), gk, gk_opt
   integer*8 :: nQpoints, Kmap(1:16)
   double complex :: QPhase(0:15,0:8)
   double precision :: Qpoint(0:15,2)
!   integer*8 :: N
endmodule BettsCluster

module PolarParas
   double precision :: k_theta_in(0:15), k_theta_out(0:15)
   double precision :: k_phi_in(0:15), k_phi_out(0:15)
   double precision :: p_theta_in(0:15), p_phi_in(0:15)
   double precision :: p_theta_out(0:15), p_phi_out(0:15)
   double precision :: Gn_in(3,0:15), Gn_out(3,0:15)
   character*1 :: PorS_in, PorS_out
endmodule PolarParas

module ScanRegion
   integer*8 :: divX, divY
   integer*8:: eneX, eneY
   double precision :: startX, endX, startY, endY
endmodule ScanRegion




Program Main
use NumOfOrbitalAndElectrons; use ModelParas; use BettsCluster
use ConstantParas; use PolarParas; use ScanRegion
implicit none
integer*8 :: ksize_0, SprSize_0
integer*8 :: ksize_f, SprSize_f
integer*8 :: ii, jj, kk, mm, nn
integer*8 :: Ediv, Eloss_div
integer*8 :: qq, qqk, mqq
integer*8, allocatable :: Hsp_0(:), Hsp_f(:)
integer*8, external :: factorial, ksubtract, ksum
integer*8, allocatable :: IndexI_0(:), IndexJ_0(:)
integer*8, allocatable :: IndexI_f(:), IndexJ_f(:)
double precision :: a, maxE, minE, maxE_loss, minE_loss
double precision :: E_0, E_0_opt
double precision, allocatable :: specX(:), specY(:)
double precision, allocatable :: sparseH_0(:)
double precision, allocatable :: sparseH_f(:)
double precision, allocatable :: tempv(:)
double precision, allocatable :: H_0(:)
double complex :: z
character*1 :: onechar, onecharp
character*2 :: twochar


open(unit=10101,file='input_SpinCrlt',Status='old');
read(10101,'(3I8)')   N3d, N2p, Ndim
read(10101,'(2I8)')   nup, ndn
read(10101,'(2F8.2)') maxE_loss, minE_loss
read(10101,'(2F8.2)') maxE, minE
read(10101,'(2I8)') Ediv, Eloss_div
read(10101,'(3I8)')   niter_CL, niter_CFE, niter_CG
read(10101,'(2F8.6)') tol_CL, tol_CG
read(10101,'(2F8.2)') epsilone_CG, epsilone_CFE
close(10101)


write(*,*) ''
write(*,*) '     N3d     N2p     Ndim '
write(*,'(3I8)') N3d, N2p, Ndim
write(*,*) ''
write(*,*) '     nup     ndn'
write(*,'(2I8)') nup, ndn
write(*,*) ''
write(*,*) 'maxE_loss minE_loss'
write(*,'(2F8.2)') maxE_loss,minE_loss
write(*,*) ''
write(*,*) '    maxE    minE'
write(*,'(2F8.2)') maxE,minE
write(*,*) ''
write(*,*) '    Ediv    Eloss_div'
write(*,'(2I8)') Ediv, ELoss_div
write(*,*) ''
write(*,*) 'niter_CL     CFE      CG'
write(*,'(3I8)') niter_CL, niter_CFE, niter_CG
write(*,*) ''
write(*,*) '  tol_CL  tol_CG'
write(*,'(2F8.2)') tol_CL, tol_CG
write(*,*) ''
write(*,*) 'epsilone_CG, CFE'
write(*,'(2F8.2)') epsilone_CG, epsilone_CFE
write(*,*) ''
write(*,*) '******************************'
write(*,*) ''

N2p=0
N=N3d+N2p
ntol = nup+ndn
startX = minE_loss; endX = maxE_loss; divX = Eloss_div-1
startY = minE; endY = maxE; divY = Ediv-1

!============== The Step 0: preparation ==============


allocate(Hsp_0(Hsize_3dEst))
allocate(Hsp_f(Hsize_3dEst))

allocate(H_0(Hsize_3dEst))

allocate(IndexI_0(SprSize_3dEst))
allocate(IndexJ_0(SprSize_3dEst))
allocate(sparseH_0(SprSize_3dEst))
SprSize_0=SprSize_3dEst

allocate(IndexI_f(SprSize_3dEst))
allocate(IndexJ_f(SprSize_3dEst))
allocate(sparseH_f(SprSize_3dEst))
SprSize_f=SprSize_3dEst

allocate(specX(divX+1))
allocate(specY(divX+1))

allocate(tempv(Hsize_3dEst))

E_0=0.0d0;
H_0=0.0d0;
Hsp_0=0
Hsp_f=0
specX=0.0d0;
specY=0.0d0
tempv=0.0d0
IndexI_0=0
IndexJ_0=0
sparseH_0=0.0d0


nbx = 10
nby = 10
do ibx = 1,nbx
do iby = 1,nby

   write(*,*) 'ibx and iby', ibx, iby
   call SetModelParameters

   !prefac(1) = dcmplx(cos(ibx*2*pi/nbx), sin(ibx*2*pi/nbx))
   !prefac(2) = dcmplx(cos(ibx*2*pi/nbx), -1*sin(ibx*2*pi/nbx))
   !prefac(3) = dcmplx(cos(iby*2*pi/nby), sin(iby*2*pi/nby))
   !prefac(4) = dcmplx(cos(iby*2*pi/nby), -1*sin(iby*2*pi/nby))

   !============== The Step 1: Diagonalize 3d system ==============
   
   ! 1 stands for (0,0) momentum
   
   ! | phi 0 > at k = 0
   gk_opt = 1
   E_0_opt = 100
   do ii=1, nQpoints
      gk = Kmap(ii); 
      call GenHsp_kspace(gk, Hsp_0, ksize_0)
      call GenMatrix_kspace(gk, Hsp_0, ksize_0, IndexI_0, IndexJ_0, sparseH_0, SprSize_0)
      call ED_ARPACK(ksize_0, SprSize_0, IndexI_0, IndexJ_0, sparseH_0, E_0, H_0)
      write(*,*) ' '
      write(*,*) 'momentum point = ', gk
      !write(*,*) 'ksize_0 = ', ksize_0
      !write(*,*) 'SprSize_0 = ', SprSize_0
      write(*,*) 'E_0 = ', E_0
      if(E_0 < E_0_opt) then
         E_0_opt = E_0
         gk_opt = gk
      endif
         
   enddo
   
   gk = gk_opt; 
   call GenHsp_kspace(gk, Hsp_0, ksize_0)
   call GenMatrix_kspace(gk, Hsp_0, ksize_0, IndexI_0, IndexJ_0, sparseH_0, SprSize_0)
   call ED_ARPACK(ksize_0, SprSize_0, IndexI_0, IndexJ_0, sparseH_0, E_0, H_0)
   write(*,*) ' '
   write(*,*) 'optimal momentum point = ', gk
   write(*,*) 'ksize_0 = ', ksize_0
   write(*,*) 'SprSize_0 = ', SprSize_0
   write(*,*) 'E_0 = ', E_0
   
   !do ii=1,SprSize_0
   !   write(*,'(2B28, F10.3)') Hsp_0(IndexI_0(ii)), Hsp_0(IndexJ_0(ii)), sparseH_0(ii)
   !enddo
   !call ComplexLanczos(Hsize_3d_spin0, SprSize_3d_spin0, IndexI_3d_spin0, &
   !                 IndexJ_3d_spin0, sparseH_3d_spin0, E_0_3d, H_0_3d)
   
   
   do ii=1, nQpoints
      
      qq = Kmap(ii) 
   
      write(*,*) '**********************'
      write(*,*) '     qq  =  ', kpoint(qq,1), kpoint(qq,2)
      write(*,*) '**********************'
   
      qqk = ksubtract(gk, qq) !qqk = gk - qq
      mqq = ksubtract(int(0,8), qq) !mqq = -qq
      write(*,*) 'qqk = ', kpoint(qqk,1), kpoint(qqk,2)
      write(*,*) 'mqq = ', kpoint(mqq,1), kpoint(mqq,2)
   
      call GenHsp_kspace(qqk, Hsp_f, ksize_f)
      SprSize_f=SprSize_3dEst
      call GenMatrix_kspace(qqk, Hsp_f, ksize_f, IndexI_f, IndexJ_f, sparseH_f, SprSize_f)
      write(*,*) 'ksize_f = ', ksize_f
      write(*,*) 'SprSize_f = ', SprSize_f
   
   
      call GenSpinCorrelatedState(qq, H_0, tempv, Hsp_0, Hsp_f, ksize_0, ksize_f)
      write(*,*) 'Static spin correlation', qq, dot_product(tempv, tempv)   
      call ContFracExpan(ksize_f, SprSize_f, tempv, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
      write(onechar,'(i1)') ii-1
      if(((ibx-1)*10+(iby-1)) > 9) then
         write(twochar,'(i2)') (ibx-1)*10+(iby-1)
         open (111,file='Sqw-'//trim(twochar)//'-'//trim(onechar)//'.dat')
      else if(((ibx-1)*10+(iby-1)) > 0) then
         write(onecharp,'(i1)') (ibx-1)*10+(iby-1)
         open (111,file='Sqw-0'//trim(onecharp)//'-'//trim(onechar)//'.dat')
      else
         open (111,file='Sqw-00-'//trim(onechar)//'.dat')
      endif
      do mm=1, divX+1
         write(111,*) specX(mm), specY(mm)
      enddo
      close(111)
   !
   !   write(*,*) 'GenSpinThreeMagnon'
   !   call GenSpinThreeMagnon(qq, H_0, tempv, Hsp_0, Hsp_f, ksize_0, ksize_f)
   !   call ContFracExpan(ksize_f, SprSize_f, tempv, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
   !   write(onechar,'(i1)') qq
   !   open (112,file='ThreeMagnon-'//trim(onechar)//'.dat')
   !   do mm=1, divX+1
   !      write(112,*) specX(mm), specY(mm)
   !   enddo
   !   close(112)
   !
   !
   !   write(*,*) 'GenSpinCharge'
   !   call GenSpinCharge(qq, H_0, tempv, Hsp_0, Hsp_f, ksize_0, ksize_f)
   !   call ContFracExpan(ksize_f, SprSize_f, tempv, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
   !   write(onechar,'(i1)') qq
   !   open (113,file='SpinCharge-'//trim(onechar)//'.dat')
   !   do mm=1, divX+1
   !      write(113,*) specX(mm), specY(mm)
   !   enddo
   !   close(113)
   !
   !
   !   write(*,*) 'GenChargeCorrelatedState'
   !   call GenChargeCorrelatedState(qq, H_0, tempv, Hsp_0, Hsp_f, ksize_0, ksize_f)
   !   call ContFracExpan(ksize_f, SprSize_f, tempv, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
   !   write(onechar,'(i1)') qq
   !   open (114,file='Nqw-'//trim(onechar)//'.dat')
   !   do mm=1, divX+1
   !      write(114,*) specX(mm), specY(mm)
   !   enddo
   !   close(114)
   !
   !
   !   write(*,*) 'GenChargeTilda'
   !   call GenChargeTilda(qq, H_0, tempv, Hsp_0, Hsp_f, ksize_0, ksize_f)
   !   call ContFracExpan(ksize_f, SprSize_f, tempv, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
   !   write(onechar,'(i1)') qq
   !   open (1141,file='Nqwtilda-'//trim(onechar)//'.dat')
   !   do mm=1, divX+1
   !      write(1141,*) specX(mm), specY(mm)
   !   enddo
   !   close(1141)
   !
   !
   !
   !   write(*,*) 'GenTwoMagnon'
   !   call GenTwoMagnon(qq, H_0, tempv, Hsp_0, Hsp_f, ksize_0, ksize_f)
   !   call ContFracExpan(ksize_f, SprSize_f, tempv, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
   !   write(onechar,'(i1)') qq
   !   open (115,file='TwoMagnon-'//trim(onechar)//'.dat')
   !   do mm=1, divX+1
   !      write(115,*) specX(mm), specY(mm)
   !   enddo
   !   close(115)
   
   
   !   write(*,*) 'GenChargeChargeTilda'
   !   call GenChargeChargeTilda(qq, H_0, tempv, Hsp_0, Hsp_f, ksize_0, ksize_f)
   !   call ContFracExpan(ksize_f, SprSize_f, tempv, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
   !   write(onechar,'(i1)') qq
   !   open (116,file='ChargeChargeTilda-'//trim(onechar)//'.dat')
   !   do mm=1, divX+1
   !      write(116,*) specX(mm), specY(mm)
   !   enddo
   !   close(116)
   
   
   enddo  !nQpoints for loop
enddo
enddo

end program





subroutine GenChargeChargeTilda(qq, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas
use ModelParas
implicit none

INTEGER*8, INTENT(IN) :: qq, ksize0, ksize
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double precision, DIMENSION(ksize0), INTENT(IN) :: groundH
double precision, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, kk, kkp, kk1, kk2, kk3, kk4
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp, kstat, ktempx, ktempy

midH=0.0d0
do jj=1,ksize0          

   write(*,*) jj

   groundstatup=listqpt0(jj)/2**N    
   groundstatdn=mod(listqpt0(jj),2**N)

   do kk=0,N-1
      ! code to get kkp
      kkp = ksum(kk, qq)
      !write(*,*) 'kk    = ', kpoint(kk,1), kpoint(kk,2)
      !write(*,*) 'kk+qq = ', kpoint(kkp,1), kpoint(kkp,2)
   

! next***               !ii  sum over n
                        !            -----         -----
                        !             \             \     +                 2
                        !             /    < phi0 | /    C  C     | phi n >
                        !            -----         -----  k  k-q
!********               !              n             k

      if(BTEST(groundstatup,kkp).eqv..true.) then
         midstatup=IBCLR(groundstatup,kkp)
         midstatdn=groundstatdn
         midsign=(-1)**(sumbeforebit(midstatup,kkp)+sumeverybit(midstatdn))

         if(BTEST(midstatup,kk).eqv..false.) then
            afterstatup=IBSET(midstatup,kk)
            afterstatdn=midstatdn
            aftersign=midsign*(-1)**(sumbeforebit(midstatup,kk)+sumeverybit(midstatdn))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)*lsf(kkp)
         endif
      endif

      if(BTEST(groundstatdn,kkp).eqv..true.) then
         midstatdn=IBCLR(groundstatdn,kkp)
         midstatup=groundstatup
         midsign=(-1)**(sumbeforebit(midstatdn,kkp))

         if(BTEST(midstatdn,kk).eqv..false.) then
            afterstatdn=IBSET(midstatdn,kk)
            afterstatup=midstatup
            aftersign=midsign*(-1)**(sumbeforebit(midstatdn,kk))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)*lsf(kkp)
         endif
      endif


      do kk1 = 0,N-1
      do kk2 = 0,N-1

         kk3 = ksubtract(ksum(kk2, kk1), kk)
   !      kkp = ksubtract(kk, qq)
         kk4 = kkp

         if(BTEST(groundstatup,kk4).eqv..true.) then
            midstatup=IBCLR(groundstatup,kk4)
            midstatdn=groundstatdn
            midsign=(-1)**(sumbeforebit(midstatup,kk4)+sumeverybit(midstatdn))

         if(BTEST(midstatup,kk3).eqv..true.) then
            midstatup=IBCLR(midstatup,kk3)
            midstatdn=midstatdn
            midsign=midsign*(-1)**(sumbeforebit(midstatup,kk3)+sumeverybit(midstatdn))

         if(BTEST(midstatup,kk2).eqv..false.) then
            midstatup=IBSET(midstatup,kk2)
            midstatdn=midstatdn
            midsign=midsign*(-1)**(sumbeforebit(midstatup,kk2)+sumeverybit(midstatdn))

         if(BTEST(midstatup,kk1).eqv..false.) then
            afterstatup=IBSET(midstatup,kk1)
            afterstatdn=midstatdn
            aftersign=midsign*(-1)**(sumbeforebit(midstatup,kk1)+sumeverybit(midstatdn))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'up, No overlap between groundstate and N-1 state'
               write(*,*) 'kk=', kpoint(kk,1), kpoint(kk,2)
               write(*,*) 'qq=', kpoint(qq,1), kpoint(qq,2)
               write(*,*) 'kk1=', kpoint(kk1,1), kpoint(kk1,2)
               write(*,*) 'kk2=', kpoint(kk2,1), kpoint(kk2,2)
               write(*,*) 'kk3=', kpoint(kk3,1), kpoint(kk3,2)
               write(*,*) 'kk4=', kpoint(kk4,1), kpoint(kk4,2)
               write(*,'(2B16)') groundstatup, groundstatdn
               write(*,'(2B16)') afterstatup, afterstatdn
               stop
            endif
            midH(l)=midH(l) - 1.0/N*aftersign*groundH(jj)*lsf(kkp)
         endif
         endif
         endif
         endif

         if(BTEST(groundstatdn,kk4).eqv..true.) then
            midstatdn=IBCLR(groundstatdn,kk4)
            midstatup=groundstatup
            midsign=midsign*(-1)**(sumbeforebit(midstatdn,kk4))

         if(BTEST(midstatdn,kk3).eqv..true.) then
            midstatdn=IBCLR(midstatdn,kk3)
            midstatup=midstatup
            midsign=midsign*(-1)**(sumbeforebit(midstatdn,kk3))

         if(BTEST(midstatdn,kk2).eqv..false.) then
            midstatdn=IBSET(midstatdn,kk2)
            midstatup=midstatup
            midsign=midsign*(-1)**(sumbeforebit(midstatdn,kk2))

         if(BTEST(midstatdn,kk1).eqv..false.) then
            afterstatdn=IBSET(midstatdn,kk1)
            afterstatup=midstatup
            aftersign=midsign*(-1)**(sumbeforebit(midstatdn,kk1))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'dn, No overlap between groundstate and N-1 state'
               write(*,*) 'kk=', kpoint(kk,1), kpoint(kk,2)
               write(*,*) 'qq=', kpoint(qq,1), kpoint(qq,2)
               write(*,*) 'kk1=', kpoint(kk1,1), kpoint(kk1,2)
               write(*,*) 'kk2=', kpoint(kk2,1), kpoint(kk2,2)
               write(*,*) 'kk3=', kpoint(kk3,1), kpoint(kk3,2)
               write(*,*) 'kk4=', kpoint(kk4,1), kpoint(kk4,2)
               write(*,'(2B16)') groundstatup, groundstatdn
               write(*,'(2B16)') afterstatup, afterstatdn
               stop
            endif
            midH(l)=midH(l) - 1.0/N*aftersign*groundH(jj)*lsf(kkp)
         endif
         endif
         endif
         endif

      enddo
      enddo

   enddo

enddo

end subroutine






subroutine GenChargeCharge(qq, H_0, H_f, Hsp_0, Hsp_f, ksize_0, ksize_f)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas
use ModelParas
implicit none

INTEGER*8, INTENT(IN) :: qq, ksize_0, ksize_f
integer*8, DIMENSION(ksize_0), INTENT(IN) :: Hsp_0
integer*8, DIMENSION(ksize_f), INTENT(IN) :: Hsp_f
double precision, DIMENSION(ksize_0), INTENT(IN) :: H_0
double precision, DIMENSION(ksize_f), INTENT(OUT) :: H_f

integer*8, external :: ksubtract, ksum
integer*8 :: jj, kk, qqp, qqt
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py, spin
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: qqt1l, qqt1r, qqt2, qqs1, qqs2
integer*8 :: ksize_s1, ksize_s2
integer*8, allocatable :: Hsp_s1(:), Hsp_s2(:)
double precision, allocatable :: tempv1(:), tempv2(:)

allocate(tempv1(Hsize_3dEst))
allocate(tempv2(Hsize_3dEst))
allocate(Hsp_s1(Hsize_3dEst))
allocate(Hsp_s2(Hsize_3dEst))

H_f = 0.0
   !write(*,*) 'gk   = ', kpoint(gk,1), kpoint(gk,2)
   !write(*,*) 'kk   = ', kpoint(kk,1), kpoint(kk,2)

do kk=0,N-1
   !write(*,*) 'qq   = ', kpoint(qq,1), kpoint(qq,2)
   !write(*,*) 'qqp  = ', kpoint(qqp,1), kpoint(qqp,2)
   !write(*,*) 'qq and qqp', qq, qqp

   tempv1 = 0.0
   qqt1l = qq
   qqt1r = ksum(qq, kk)  !transfered momentum
   !qqs1 = ksubtract(gk, ksum(qqp,kk))            !momentum in the excited state
   !call GenHsp_kspace(qqs1, Hsp_s1, ksize_s1)
   !write(*,*) 'qq = ', qq
   call GenTwoCharge(qqt1l, qqt1r, H_0, tempv1, Hsp_0, Hsp_f, ksize_0, ksize_f)

   H_f = H_f + tempv1 * lsf(qqt1r)  !lattice structure factor

enddo

deallocate(tempv1)
deallocate(tempv2)
deallocate(Hsp_s1)
deallocate(Hsp_s2)

end subroutine



subroutine GenTwoMagnon(kk, H_0, H_f, Hsp_0, Hsp_f, ksize_0, ksize_f)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas
use ModelParas
implicit none

INTEGER*8, INTENT(IN) :: kk, ksize_0, ksize_f
integer*8, DIMENSION(ksize_0), INTENT(IN) :: Hsp_0
integer*8, DIMENSION(ksize_f), INTENT(IN) :: Hsp_f
double precision, DIMENSION(ksize_0), INTENT(IN) :: H_0
double precision, DIMENSION(ksize_f), INTENT(OUT) :: H_f

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, qq, qqp, qqt
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py, spin
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: qqt1, qqt2, qqs1, qqs2
integer*8 :: ksize_s1, ksize_s2
integer*8, allocatable :: Hsp_s1(:), Hsp_s2(:)
double precision, allocatable :: tempv1(:), tempv2(:), H_i(:)


allocate(tempv1(Hsize_3dEst))
allocate(tempv2(Hsize_3dEst))
allocate(Hsp_s1(Hsize_3dEst))
allocate(Hsp_s2(Hsize_3dEst))
allocate(H_i(Hsize_3dEst))

H_f = 0.0

   do qq=0,N-1
      write(*,*) 'qq', qq

      H_i = 0.0d0
            
      ! Eq (13) line 3
      !qqt1 = ksubtract(int(0,8),ksum(kk, qq))  !transfered momentum
      qqt1 = ksubtract(kk, qq)              !qqt1 = kk - qq
      qqs1 = ksubtract(gk, qqt1)            !qqs1 = gk - kk + qq
      call GenHsp_kspace(qqs1, Hsp_s1, ksize_s1)
      call GenSpinCorrelatedState(qqt1, H_0, tempv1, Hsp_0, Hsp_s1, ksize_0, ksize_s1)
      !write(*,*) 'Eq 13 line 3 done'    

      ! Eq (13) line 2
      qqt2 = qq
      qqs2 = ksubtract(qqs1, qqt2)         !qqs2 = qqs1 - qq = gk - kk
      call GenHsp_kspace(qqs2, Hsp_s2, ksize_s2)
      call GenSpinCorrelatedState(qqt2, tempv1, tempv2, Hsp_s1, Hsp_f, ksize_s1, ksize_f)
      H_i = H_i + tempv2/4
      !write(*,*) 'Eq 13 line 2 done'    

      ! Eq (13) line 1         
      ! d^{+}_q2,s d_{q2+qqt2},sbar d^{+}_q1,sbar d_{q1+qqt1},s
      call GenFourParticle(qqt2, qqt1, H_0, tempv1, Hsp_0, Hsp_f, ksize_0, ksize_f)
      H_i = H_i + tempv1/2
      !write(*,*) 'Eq 13 line 1 done'    

      tempv2 = H_i * lsf(ksubtract(kk, qq))  !lattice structure factor

      H_f = H_f + tempv2

      write(*,*) 'gk - kk = ', ksubtract(gk, kk), '  qqs2 = ', qqs2


   enddo

deallocate(tempv1)
deallocate(tempv2)
deallocate(Hsp_s1)
deallocate(Hsp_s2)
deallocate(H_i)

end subroutine



subroutine GenSpinCharge(kk, H_0, H_f, Hsp_0, Hsp_f, ksize_0, ksize_f)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas
use ModelParas
implicit none

INTEGER*8, INTENT(IN) :: kk, ksize_0, ksize_f
integer*8, DIMENSION(ksize_0), INTENT(IN) :: Hsp_0
integer*8, DIMENSION(ksize_f), INTENT(IN) :: Hsp_f
double precision, DIMENSION(ksize_0), INTENT(IN) :: H_0
double precision, DIMENSION(ksize_f), INTENT(OUT) :: H_f

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, qq, qqp, qqt
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py, spin
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: qqt1l, qqt1r, qqt2, qqs1, qqs2
integer*8 :: ksize_s1, ksize_s2
integer*8, allocatable :: Hsp_s1(:), Hsp_s2(:)
double precision, allocatable :: tempv1(:), tempv2(:)


allocate(tempv1(Hsize_3dEst))
allocate(tempv2(Hsize_3dEst))
allocate(Hsp_s1(Hsize_3dEst))
allocate(Hsp_s2(Hsize_3dEst))

H_f = 0.0
   !write(*,*) 'gk   = ', kpoint(gk,1), kpoint(gk,2)
   !write(*,*) 'kk   = ', kpoint(kk,1), kpoint(kk,2)

   do qq=0,N-1
   do qqp=0,N-1
      !write(*,*) 'qq   = ', kpoint(qq,1), kpoint(qq,2)
      !write(*,*) 'qqp  = ', kpoint(qqp,1), kpoint(qqp,2)
      !write(*,*) 'qq and qqp', qq, qqp

      qqt1l = qq
      qqt1r = ksubtract(ksum(qq, kk), qqp)  !transfered momentum
      qqs1 = ksubtract(gk, ksubtract(kk, qqp))            !momentum in the excited state
      call GenHsp_kspace(qqs1, Hsp_s1, ksize_s1)
      call GenTwoCharge(qqt1l, qqt1r, H_0, tempv1, Hsp_0, Hsp_s1, ksize_0, ksize_s1)

      qqt2 = qqp
      qqs2 = ksubtract(qqs1, qqt2)
      !call GenHsp_kspace(qqs2, Hsp_s2, ksize_s2)
      call GenSpinCorrelatedState(qqt2, tempv1, tempv2, Hsp_s1, Hsp_f, ksize_s1, ksize_f)

      tempv1 = tempv2 * lsf(ksubtract(ksum(qq, kk), qqp))  !lattice structure factor
            
      H_f = H_f + tempv1

   enddo
   enddo

deallocate(tempv1)
deallocate(tempv2)
deallocate(Hsp_s1)
deallocate(Hsp_s2)

end subroutine



subroutine GenSpinThreeMagnon(kk, H_0, H_f, Hsp_0, Hsp_f, ksize_0, ksize_f)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas
use ModelParas
implicit none

INTEGER*8, INTENT(IN) :: kk, ksize_0, ksize_f
integer*8, DIMENSION(ksize_0), INTENT(IN) :: Hsp_0
integer*8, DIMENSION(ksize_f), INTENT(IN) :: Hsp_f
double precision, DIMENSION(ksize_0), INTENT(IN) :: H_0
double precision, DIMENSION(ksize_f), INTENT(OUT) :: H_f

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, qq, qqp, qqt
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py, spin
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: qqt1, qqt2, qqt3, qqs1, qqs2, qqs3
integer*8 :: ksize_s1, ksize_s2, ksize_s3
integer*8, allocatable :: Hsp_s1(:), Hsp_s2(:), Hsp_s3(:)
double precision, allocatable :: tempv1(:), tempv2(:), H_i(:)


allocate(tempv1(Hsize_3dEst))
allocate(tempv2(Hsize_3dEst))
allocate(Hsp_s1(Hsize_3dEst))
allocate(Hsp_s2(Hsize_3dEst))
allocate(Hsp_s3(Hsize_3dEst))
allocate(H_i(Hsize_3dEst))

H_f = 0.0
   !write(*,*) 'gk   = ', kpoint(gk,1), kpoint(gk,2)
   !write(*,*) 'kk   = ', kpoint(kk,1), kpoint(kk,2)

   do qq=0,N-1
   do qqp=0,N-1
      !write(*,*) 'qq   = ', kpoint(qq,1), kpoint(qq,2)
      !write(*,*) 'qqp  = ', kpoint(qqp,1), kpoint(qqp,2)
      write(*,*) 'qq and qqp', qq, qqp

      H_i = 0.0d0
            
      ! Eq (13) line 3
      qqt1 = ksubtract(kk, ksum(qqp, qq))   !transfered momentum
      qqs1 = ksubtract(gk, qqt1)            !momentum in the excited state
      call GenHsp_kspace(qqs1, Hsp_s1, ksize_s1)
      call GenSpinCorrelatedState(qqt1, H_0, tempv1, Hsp_0, Hsp_s1, ksize_0, ksize_s1)
      !write(*,*) 'Eq 13 line 3 done'    

      ! Eq (13) line 2
      qqt2 = qq
      qqs2 = ksubtract(qqs1, qqt2)
      call GenHsp_kspace(qqs2, Hsp_s2, ksize_s2)
      call GenSpinCorrelatedState(qqt2, tempv1, tempv2, Hsp_s1, Hsp_s2, ksize_s1, ksize_s2)
      H_i = H_i + tempv2/4
      !write(*,*) 'Eq 13 line 2 done'    

      ! Eq (13) line 1         
      ! d^{+}_q2,s d_{q2+qqt2},sbar d^{+}_q1,sbar d_{q1+qqt1},s
      call GenFourParticle(qqt2, qqt1, H_0, tempv1, Hsp_0, Hsp_s2, ksize_0, ksize_s2)
      H_i = H_i + tempv1/2
      !write(*,*) 'Eq 13 line 1 done'    

      ! S_{q'}^z
      qqt3 = qqp
      qqs3 = ksubtract(qqs2, qqt3)
      !call GenHsp_kspace(qqs3, Hsp_s3, ksize_s3)
      call GenSpinCorrelatedState(qqt3, H_i, tempv1, Hsp_s2, Hsp_f, ksize_s2, ksize_f)
      !write(*,*) 'Eq S_{q}^z done'    
            
      tempv2 = tempv1 * lsf(ksum(qqp, ksubtract(qq, kk)))  !lattice structure factor

      H_f = H_f + tempv2

   enddo
   enddo

deallocate(tempv1)
deallocate(tempv2)
deallocate(Hsp_s1)
deallocate(Hsp_s2)
deallocate(Hsp_s3)
deallocate(H_i)

end subroutine


! d^{+}_q2,s d_{q2+qqt2},sbar d^{+}_q1,sbar d_{q1+qqt1},s
subroutine GenFourParticle(qq2t, qq1t, groundH, midH, listqpt0, listqpt, ksize0, ksize)
use NumOfOrbitalAndElectrons; use BettsCluster
implicit none

INTEGER*8, INTENT(IN) :: qq1t, qq2t, ksize0, ksize
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double precision, DIMENSION(ksize0), INTENT(IN) :: groundH
double precision, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj
integer*8 :: qq1s, qq2s, qq1, qq2
integer*8 :: initstatup, initstatdn
integer*8 :: afterstatup,afterstatdn,aftersign,l


     midH(:)=0.0d0
     do jj= 1, ksize0  !loop B

        initstatup=listqpt0(jj)/2**N;
        initstatdn=mod(listqpt0(jj),2**N);

        do qq1=0,N-1
        do qq2=0,N-1
   
           qq1s = ksum(qq1, qq1t)
           qq2s = ksum(qq2, qq2t)

              if(BTEST(initstatdn,qq1s)) then   !dn
                 afterstatdn=IBCLR(initstatdn,qq1s)
                 aftersign=((-1)**sumbeforebit(initstatdn,qq1s))
              if(.not.BTEST(initstatup,qq1)) then  !up
                 afterstatup=IBSET(initstatup,qq1)
                 aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,qq1)+sumeverybit(afterstatdn)))
              if(BTEST(afterstatup,qq2s)) then   !up
                 afterstatup=IBCLR(afterstatup,qq2s)
                 aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,qq2s)+sumeverybit(afterstatdn)))
              if(.not.BTEST(afterstatdn,qq2)) then !dn
                 afterstatdn=IBSET(afterstatdn,qq2)
                 aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,qq2))

                 call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                 if(l.eq.-1) then
                    write(*,*) 'No overlap between groundstate and N-1 state'
                    write(*,*) 'qqt1 & qqt2 =',qq1t,qq2t
                    stop
                 endif
                 midH(l)=midH(l)+aftersign*groundH(jj)

              endif                 
              endif                 
              endif                 
              endif                 

              if(BTEST(initstatup,qq1s)) then  !up
                 afterstatup=IBCLR(initstatup,qq1s)
                 aftersign=((-1)**(sumbeforebit(afterstatup,qq1s)+sumeverybit(initstatdn)))
              if(.not.BTEST(initstatdn,qq1)) then   !dn
                 afterstatdn=IBSET(initstatdn,qq1)
                 aftersign=aftersign*((-1)**sumbeforebit(initstatdn,qq1))
              if(BTEST(afterstatdn,qq2s)) then !dn
                 afterstatdn=IBCLR(afterstatdn,qq2s)
                 aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,qq2s))
              if(.not.BTEST(afterstatup,qq2)) then   !up
                 afterstatup=IBSET(afterstatup,qq2)
                 aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,qq2)+sumeverybit(afterstatdn)))

                 call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                 if(l.eq.-1) then
                    write(*,*) 'No overlap between groundstate and N-1 state'
                    write(*,*) 'qqt1 & qqt2 =',qq1t,qq2t
                    stop
                 endif
                 midH(l)=midH(l)+aftersign*groundH(jj)

              endif                 
              endif                 
              endif                 
              endif  
               
           enddo
           enddo

        enddo !loop B

end subroutine


subroutine GenTwoCharge(qql, qqr, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize)
use NumOfOrbitalAndElectrons; use BettsCluster
implicit none

INTEGER*8, INTENT(IN) :: qql, qqr, ksize0, ksize
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double precision, DIMENSION(ksize0), INTENT(IN) :: groundH
double precision, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, kk, kkp
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp, kstat, ktempx, ktempy



midH=0.0d0
!do kk=0,N-1
kk = qql
kkp = qqr
   
   do jj=1,ksize0          

    groundstatup=listqpt0(jj)/2**N
    groundstatdn=mod(listqpt0(jj),2**N)

      if(BTEST(groundstatup,kkp).eqv..true.) then
         midstatup=IBCLR(groundstatup,kkp)
         midstatdn=groundstatdn
         midsign=(-1)**(sumbeforebit(midstatup,kkp)+sumeverybit(midstatdn))

         if(BTEST(midstatup,kk).eqv..false.) then
            afterstatup=IBSET(midstatup,kk)
            afterstatdn=midstatdn
            aftersign=midsign*(-1)**(sumbeforebit(midstatup,kk)+sumeverybit(midstatdn))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif

      if(BTEST(groundstatdn,kkp).eqv..true.) then
         midstatdn=IBCLR(groundstatdn,kkp)
         midstatup=groundstatup
         midsign=(-1)**(sumbeforebit(midstatdn,kkp))

         if(BTEST(midstatdn,kk).eqv..false.) then
            afterstatdn=IBSET(midstatdn,kk)
            afterstatup=midstatup
            aftersign=midsign*(-1)**(sumbeforebit(midstatdn,kk))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif

   enddo

end subroutine


subroutine GenChargeTilda(qq, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize)
use NumOfOrbitalAndElectrons; use BettsCluster
implicit none

INTEGER*8, INTENT(IN) :: qq, ksize0, ksize
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double precision, DIMENSION(ksize0), INTENT(IN) :: groundH
double precision, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, kk, kkp, kk1, kk2, kk3, kk4
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp, kstat, ktempx, ktempy

midH=0.0d0
do jj=1,ksize0          

   groundstatup=listqpt0(jj)/2**N    
   groundstatdn=mod(listqpt0(jj),2**N)

   do kk=0,N-1
      ! code to get kkp
      kkp = ksum(kk, qq)
      !write(*,*) 'kk    = ', kpoint(kk,1), kpoint(kk,2)
      !write(*,*) 'kk+qq = ', kpoint(kkp,1), kpoint(kkp,2)
   

! next***               !ii  sum over n
                        !            -----         -----
                        !             \             \     +                 2
                        !             /    < phi0 | /    C  C     | phi n >
                        !            -----         -----  k  k-q
!********               !              n             k

      if(BTEST(groundstatup,kkp).eqv..true.) then
         midstatup=IBCLR(groundstatup,kkp)
         midstatdn=groundstatdn
         midsign=(-1)**(sumbeforebit(midstatup,kkp)+sumeverybit(midstatdn))

         if(BTEST(midstatup,kk).eqv..false.) then
            afterstatup=IBSET(midstatup,kk)
            afterstatdn=midstatdn
            aftersign=midsign*(-1)**(sumbeforebit(midstatup,kk)+sumeverybit(midstatdn))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif

      if(BTEST(groundstatdn,kkp).eqv..true.) then
         midstatdn=IBCLR(groundstatdn,kkp)
         midstatup=groundstatup
         midsign=(-1)**(sumbeforebit(midstatdn,kkp))

         if(BTEST(midstatdn,kk).eqv..false.) then
            afterstatdn=IBSET(midstatdn,kk)
            afterstatup=midstatup
            aftersign=midsign*(-1)**(sumbeforebit(midstatdn,kk))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif


      do kk1 = 1,N
      do kk2 = 1,N

         kk3 = ksubtract(kk1, kk)
         kk4 = ksubtract(ksum(kk2, kk), qq)

         if(BTEST(groundstatup,kk2).eqv..true.) then
            midstatup=IBCLR(groundstatup,kk2)
            midstatdn=groundstatdn
            midsign=(-1)**(sumbeforebit(midstatup,kk2)+sumeverybit(midstatdn))

         if(BTEST(midstatdn,kk1).eqv..true.) then
            midstatdn=IBCLR(midstatdn,kk1)
            midstatup=midstatup
            midsign=midsign*(-1)**(sumbeforebit(midstatdn,kk1))

         if(BTEST(midstatdn,kk3).eqv..false.) then
            midstatdn=IBSET(midstatdn,kk3)
            midstatup=midstatup
            midsign=midsign*(-1)**(sumbeforebit(midstatdn,kk3))

         if(BTEST(midstatup,kk4).eqv..false.) then
            afterstatup=IBSET(midstatup,kk4)
            afterstatdn=midstatdn
            aftersign=midsign*(-1)**(sumbeforebit(midstatup,kk4)+sumeverybit(midstatdn))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=', kpoint(kk,1), kpoint(kk,2)
               write(*,*) 'qq=', kpoint(qq,1), kpoint(qq,2)
               write(*,*) 'kk1=', kpoint(kk1,1), kpoint(kk1,2)
               write(*,*) 'kk2=', kpoint(kk2,1), kpoint(kk2,2)
               write(*,*) 'kk3=', kpoint(kk3,1), kpoint(kk3,2)
               write(*,*) 'kk4=', kpoint(kk4,1), kpoint(kk4,2)
               write(*,'(2B16)') groundstatup, groundstatdn
               write(*,'(2B16)') afterstatup, afterstatdn
               stop
            endif
            midH(l)=midH(l) - 2.0/N*aftersign*groundH(jj)
         endif
         endif
         endif
         endif

      enddo
      enddo

   enddo

enddo

end subroutine





subroutine GenChargeCorrelatedState(qq, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize)
use NumOfOrbitalAndElectrons; use BettsCluster
implicit none

INTEGER*8, INTENT(IN) :: qq, ksize0, ksize
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double precision, DIMENSION(ksize0), INTENT(IN) :: groundH
double precision, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, kk, kkp
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp, kstat, ktempx, ktempy



midH=0.0d0
do kk=0,N-1
   ! code to get kkp
   kkp = ksum(kk, qq)
   !write(*,*) 'kk    = ', kpoint(kk,1), kpoint(kk,2)
   !write(*,*) 'kk+qq = ', kpoint(kkp,1), kpoint(kkp,2)
   
   do jj=1,ksize0          

    groundstatup=listqpt0(jj)/2**N
    groundstatdn=mod(listqpt0(jj),2**N)

! next***               !ii  sum over n
                        !            -----         -----
                        !             \             \     +                 2
                        !             /    < phi0 | /    C  C     | phi n >
                        !            -----         -----  k  k-q
!********               !              n             k

      if(BTEST(groundstatup,kkp).eqv..true.) then
         midstatup=IBCLR(groundstatup,kkp)
         midstatdn=groundstatdn
         midsign=(-1)**(sumbeforebit(midstatup,kkp)+sumeverybit(midstatdn))



         if(BTEST(midstatup,kk).eqv..false.) then
            afterstatup=IBSET(midstatup,kk)
            afterstatdn=midstatdn
            aftersign=midsign*(-1)**(sumbeforebit(midstatup,kk)+sumeverybit(midstatdn))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif

      if(BTEST(groundstatdn,kkp).eqv..true.) then
         midstatdn=IBCLR(groundstatdn,kkp)
         midstatup=groundstatup
         midsign=(-1)**(sumbeforebit(midstatdn,kkp))

         if(BTEST(midstatdn,kk).eqv..false.) then
            afterstatdn=IBSET(midstatdn,kk)
            afterstatup=midstatup
            aftersign=midsign*(-1)**(sumbeforebit(midstatdn,kk))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif

   enddo

enddo

end subroutine





! next***               !ii  sum over n
!                                    -----         -----
!                                     \             \     +                 2
!                 !   S(q,w)=         /    < phi0 | /    C  C     | phi n >
!                                    -----         -----  k  k+q
!********                              n             k


!
!                    -----
!                     \       +              +   
!                     /    ( C     C      - C    C      ) | phi 0 >
!                    -----    k,up  k+q,up   k,dn k+q,dn
!                       k

subroutine GenSpinCorrelatedState(qq, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize)
use NumOfOrbitalAndElectrons; use BettsCluster
implicit none

INTEGER*8, INTENT(IN) :: qq, ksize0, ksize
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double precision, DIMENSION(ksize0), INTENT(IN) :: groundH
double precision, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, kk, kkp
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp, kstat, ktempx, ktempy

midH(:)=0.0d0
do kk=0,N-1
   ! code to get kkp
   kkp = ksum(kk, qq)
   !write(*,*) 'kk    = ', kpoint(kk,1), kpoint(kk,2)
   !write(*,*) 'kk+qq = ', kpoint(kkp,1), kpoint(kkp,2)
   
   do jj=1,ksize0          

    groundstatup=listqpt0(jj)/2**N
    groundstatdn=mod(listqpt0(jj),2**N)

      if(BTEST(groundstatup,kkp).eqv..true.) then
         midstatup=IBCLR(groundstatup,kkp)
         midstatdn=groundstatdn
         midsign=(-1)**(sumbeforebit(midstatup,kkp)+sumeverybit(midstatdn))

         if(BTEST(midstatup,kk).eqv..false.) then
            afterstatup=IBSET(midstatup,kk)
            afterstatdn=midstatdn
            aftersign=midsign*(-1)**(sumbeforebit(midstatup,kk)+sumeverybit(midstatdn))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif

      if(BTEST(groundstatdn,kkp).eqv..true.) then
         midstatdn=IBCLR(groundstatdn,kkp)
         midstatup=groundstatup
         midsign=(-1)**(sumbeforebit(midstatdn,kkp))

         if(BTEST(midstatdn,kk).eqv..false.) then
            afterstatdn=IBSET(midstatdn,kk)
            afterstatup=midstatup
            aftersign=midsign*(-1)**(sumbeforebit(midstatdn,kk))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               stop
            endif
            midH(l)=midH(l)-aftersign*groundH(jj)
         endif
      endif

   enddo

enddo


end subroutine

 



subroutine ED_ARPACK(Hsize_3d, SprSize_3d, IndexI_3d, &
                  IndexJ_3d, sparseH_3d, E_0_3d, H_0_3d)
implicit none
INTEGER*8, INTENT(IN) :: Hsize_3d, SprSize_3d
INTEGER*8, DIMENSION(SprSize_3d), INTENT(IN) :: IndexI_3d, IndexJ_3d
DOUBLE PRECISION, DIMENSION(SprSize_3d), INTENT(IN) :: sparseH_3d
DOUBLE PRECISION, INTENT(OUT) :: E_0_3d
DOUBLE PRECISION, DIMENSION(Hsize_3d), INTENT(OUT) :: H_0_3d


integer*8 :: ii
integer*8 :: spr_size

!-----------------------------------from arpack------------------------------------------------
character        bmat*1, which*2

integer     ::   maxncv, ldv,nn
integer     ::   ido, nev, ncv, lworkl, info, ierr, j, nconv, maxitr, mode, ishfts
integer     ::   iparam(11), ipntr(11)
integer     ::   flag, tempuse1

double precision  ::   tol, sigma, zero, time1
double precision, allocatable   ::   v(:,:), workl(:), workd(:), d(:,:)
double precision, allocatable   ::   resid(:)

logical, allocatable :: select(:)
logical     ::   rvec



spr_size=SprSize_3d !***number of non-zero elements of the matrix

ldv=Hsize_3d
nn=Hsize_3d
nev=5 !***number of eigen-states requied.
ncv=20
maxncv=ncv
lworkl=ncv*(ncv+8)

allocate(v(Hsize_3d,ncv))
allocate(workd(3*Hsize_3d))
allocate(d(ncv,2))
allocate(resid(Hsize_3d))
allocate(select(ncv))
allocate(workl(lworkl))


call CPU_time(time1)

zero=0.0d0;

bmat = 'I'
which = 'SA'!Require the smallest algebraic eigen-values.


!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%

!     tol = zero
      tol = 1.0d-8 !The tolerance factor for the convergence of Lanczos processes. If set to be zero, the error can be < 1.0d-16.
        info = 0
      ido = 0

!allocate(workl(lworkl))


!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%

      ishfts = 1
      maxitr = 300 !Max numver of iterations of the Lanczos processes.
        mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
      ipntr(1:11)=0

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%

 10   continue

!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%

         call dsaupd ( ido, bmat, nn, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl,lworkl, info )

         if (ido .eq. -1 .or. ido .eq. 1) then

!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%

!*********************************************************************************
        workd(ipntr(2):ipntr(2)+Hsize_3d-1)=0
        do ii=1, spr_size
        workd(ipntr(2)+IndexI_3d(ii)-1)=workd(ipntr(2)+IndexI_3d(ii)-1)+workd(ipntr(1)+IndexJ_3d(ii)-1)*sparseH_3d(ii)
        end do
!*********************************************************************************


!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
            go to 10

         end if


!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%


      if ( info .lt. 0 ) then

!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%

         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '

      else

!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%

         rvec = .true.

         call dseupd ( rvec, 'All', select, d, v, ldv, sigma,bmat, nn, which, nev, tol, resid, ncv, v, ldv,iparam, ipntr, workd, workl, lworkl, ierr )

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%

       open(unit=16, file='eigenvalues.dat')
       do ii=1,nev
         write(16,*) d(ii,1)
       end do
       close(16)

       open(unit=17, file='eigenvector.dat')
!
                do ii=1,Hsize_3d
                        H_0_3d(ii)=v(ii,1)
                        write(17,*) H_0_3d(ii)
                enddo
        close(17)

                E_0_3d=d(1,1);

         if ( ierr .ne. 0) then

!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%

             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '

         else

 20          continue

!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%

!             call dmout(6, nconv, 2, d, maxncv, -6, 'Ritz values and relative residuals')
         end if

!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%


         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit',&
                   ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if

         print *, ' '
         print *, ' _SDRV1 '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', nn, spr_size
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',&
                 ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
!         print *, ' The number of converged Ritz values is ',&
!                   nconv
         print *, ' The number of Implicit Arnoldi update',&
                 ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '

        end if

!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%

end subroutine










  

subroutine ExactDiagonalization2(Hsize_3d, SprSize_3d, IndexI_3d, &
                  IndexJ_3d, sparseH_3d, E_0_3d, H_0_3d)
implicit none
INTEGER*8, INTENT(IN) :: Hsize_3d, SprSize_3d
INTEGER*8, DIMENSION(SprSize_3d), INTENT(IN) :: IndexI_3d, IndexJ_3d
DOUBLE COMPLEX, DIMENSION(SprSize_3d), INTENT(IN) :: sparseH_3d
double precision, INTENT(OUT) :: E_0_3d
DOUBLE COMPLEX, DIMENSION(Hsize_3d), INTENT(OUT) :: H_0_3d


integer*8 :: info, ii
double precision, allocatable :: WORK(:),Ht(:,:),W_ED(:)

allocate(W_ED(Hsize_3d))
write(*,*) 'Start allocate Ht'
allocate(Ht(Hsize_3d, Hsize_3d))
write(*,*) 'finish allocate Ht'
allocate(WORK(5*Hsize_3d))

Ht=0.0d0
write(*,*) 'start assign'
do ii = 1, SprSize_3d
   Ht(IndexI_3d(ii), IndexJ_3d(ii)) = real(sparseH_3d(ii))
enddo

write(*,*) 'assign OK'

If (maxval(transpose(Ht)-Ht).ne.0.0d0.and.minval(transpose(Ht)-Ht).ne.0.0d0) then
write(*,*) "The Matrix is not Hermitian!"
else
write(*,*) "Matrix check ok."
end if

info=0
call DSYEV('V','L',Hsize_3d,Ht,Hsize_3d,W_ED,WORK,5*Hsize_3d,info)
write(*,*) info
if (info.eq.0) write(*,*) 'LAPCAK routine ok.'

E_0_3d = W_ED(1)

do ii = 1, Hsize_3d
   H_0_3d(ii) = Ht(ii,1)
enddo

deallocate(Ht)

end subroutine


subroutine ExactDiagonalization(Hsize_3d, SprSize_3d, IndexI_3d, &
                  IndexJ_3d, sparseH_3d, E_0_3d, H_0_3d)
implicit none
INTEGER*8, INTENT(IN) :: Hsize_3d, SprSize_3d
INTEGER*8, DIMENSION(SprSize_3d), INTENT(IN) :: IndexI_3d, IndexJ_3d
DOUBLE COMPLEX, DIMENSION(SprSize_3d), INTENT(IN) :: sparseH_3d
double precision, INTENT(OUT) :: E_0_3d
DOUBLE COMPLEX, DIMENSION(Hsize_3d), INTENT(OUT) :: H_0_3d


integer*8 :: LWORK
integer*8 :: info, ii
DOUBLE PRECISION, allocatable ::  RWORK(:)
double complex, allocatable :: WORK(:)
double complex, allocatable :: Ht(:,:)

allocate(Ht(Hsize_3d, Hsize_3d))
Ht=cmplx(0.0d0,0.0d0)

do ii = 1, SprSize_3d
   Ht(IndexI_3d(ii), IndexJ_3d(ii)) = sparseH_3d(ii)
enddo

LWORK=3*Hsize_3d
allocate(WORK(LWORK))
allocate(RWORK(LWORK))
info=0
!call zheev('V', 'U', Hsize_3d, Ht, Hsize_3d, &
!           E_0_3d, WORK, 3*Hsize_3d, RWORK, info)
call ZHEEV('V', 'U', Hsize_3d, Ht, Hsize_3d, &
           E_0_3d, WORK, 3*Hsize_3d, RWORK, info)

if (info.eq.0) then
   write(*,*) 'ZHEEV exit successfully.'
else
   write(*,*) '#info.neq.0# info = ', info
endif

H_0_3d(1:Hsize_3d) = Ht(1:Hsize_3d,1)

deallocate(WORK, RWORK, Ht)

end subroutine




subroutine GenHsp_kspace(kfixed, listqpt, ksize)
use ConstantParas; use NumOfOrbitalAndElectrons
use BettsCluster
implicit none

INTEGER*8, INTENT(IN):: kfixed
INTEGER*8, DIMENSION(Hsize_3dEst), INTENT(OUT):: listqpt
INTEGER*8, INTENT(OUT):: ksize

integer*8 :: iiup,iidn, Hsizet
integer*8 :: listup(65536,2),listdn(65536,2),listversup(65536),listversdn(65536)
integer*8 :: ktempx,ktempy, px, py
integer*8 :: jj, ii, i, kstat, iq
integer*8 :: iktemp1,iktemp2,iktemp
integer*8 :: kstat1,kstat2
integer*8 :: temp,tempi
integer*8 :: iup,idn, jup, jdn
integer*8 :: initstatup, initstatdn
integer*8, external :: sumeverybit, factorial
integer*8 :: ktempx1,ktempx2,ktempy1,ktempy2



iiup=0;iidn=0
do ii=0,2**N-1
        listup(ii+1,1)=0;listup(ii+1,2)=0;listdn(ii+1,1)=0;listdn(ii+1,2)=0;
        listversup(ii+1)=0;listversdn(ii+1)=0
end do

ktempx=0
ktempy=0
iktemp=0
kstat=1
do iup=0,2**N-1
        if(sumeverybit(iup).eq.nup) then                !find one up state
                do i=0,N-1
                   if(BTEST(iup,i)) then
                      ktempx = ktempx + kpoint(i,1)
                      ktempy = ktempy + kpoint(i,2)
                   endif
                enddo
                ktempx = mod(ktempx,kunit)
                ktempy = mod(ktempy,kunit)
                do while (kstat.ne.0)
                   if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
                      kstat=0
                   else
                      iktemp=iktemp+1
                   endif
                   if(iktemp.eq.N) then
                      write(*,*) 'iktemp out of bounds, spin up'
                      write(*,*) iup
                      stop
                   endif
                enddo
                listup(iiup+1,1)=iup;
                listup(iiup+1,2)=iktemp;
                listversup(iup)=iiup;
                iiup=iiup+1;
                ktempx=0;
                ktempy=0;
                iktemp=0;
                kstat=1;

        endif
enddo

!********************************************************************************
! This loop is superfluous at Sz=0
!********************************************************************************

ktempx=0
ktempy=0
iktemp=0
kstat=1

do idn=0,2**N-1
        if(sumeverybit(idn).eq.ndn) then                !find one down state
                do i=0,N-1
                   if(BTEST(idn,i)) then
                      ktempx = ktempx + kpoint(i,1)
                      ktempy = ktempy + kpoint(i,2)
                   endif
                enddo
                ktempx = mod(ktempx,kunit)
                ktempy = mod(ktempy,kunit)
                do while (kstat.ne.0)
                   if(kpoint(iktemp,1).eq.ktempx.and.kpoint(iktemp,2).eq.ktempy) then
                      kstat=0
                   else
                      iktemp=iktemp+1
                   endif
                   if(iktemp.eq.N) then
                      write(*,*) 'iktemp out of bounds, spin down'
                      write(*,*) idn
                      stop
                   endif
                enddo
                listdn(iidn+1,1)=idn;
                listdn(iidn+1,2)=iktemp;
                listversdn(idn)=iidn;
                iidn=iidn+1;
                ktempx=0;
                ktempy=0;
                iktemp=0;
                kstat=1;
        endif
enddo

Hsizet = (factorial(N)/factorial(nup)/factorial(N-nup))*&
(factorial(N)/factorial(ndn)/factorial(N-ndn))

ktempx=0
ktempy=0
iktemp=0
kstat=1
iq=0
listqpt=0
do jj=0,Hsizet-1         !loop B
        jdn=mod(jj,iidn);
        jup=jj/iidn;
        initstatup=listup(jup+1,1);
        initstatdn=listdn(jdn+1,1);
        ktempx=mod(kpoint(listup(jup+1,2),1)+kpoint(listdn(jdn+1,2),1),kunit)
        ktempy=mod(kpoint(listup(jup+1,2),2)+kpoint(listdn(jdn+1,2),2),kunit)
        do while (kstat.ne.0)
           if(kpoint(iktemp,1).eq.ktempx.and.kpoint(iktemp,2).eq.ktempy) then
              kstat=0
           else
              iktemp=iktemp+1
           endif
           if(iktemp.eq.N) then
              write(*,*) 'iktemp out of bounds, column'
              write(*,*) jj
              stop
           endif
        enddo
        if(iktemp.eq.kfixed) then
           iq=iq+1
           listqpt(iq) = initstatup*(2**N)+initstatdn
        endif
        iktemp=0
        kstat=1
enddo !loop B

!write(*,*) 'listup&dn set up OK', myid
ksize = iq

!write(*,*) 'Inside GenHsp kfixed = ', kfixed
!write(*,*) 'px and py = ', kpoint(kfixed, 1), kpoint(kfixed, 2)

end subroutine




subroutine GenMatrix_kspace(kfixed, listqpt, ksize, IndexIt, IndexJt, sparseHt, SprSizet)
use NumOfOrbitalAndElectrons; use ModelParas; use MPIParas; use ConstantParas
use BettsCluster
implicit none

integer*8, INTENT(IN) :: kfixed, ksize
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
integer*8, INTENT(INOUT) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(OUT) :: IndexIt, IndexJt
DOUBLE PRECISION, DIMENSION(SprSizet), INTENT(OUT) :: sparseHt

integer*8 :: H_index(1:3000)
integer*8 :: jj, kk, l, kkp, qq, ii, pp
integer*8 :: temp, tempi
integer*8 :: initstatup,initstatdn
integer*8 :: afterstatup,afterstatdn,initsign,aftersign
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp1,iktemp2,iktemp
integer*8 :: kstat1,kstat2
integer*8 :: ktempx,ktempy
integer*8 :: ktempx1,ktempx2,ktempy1,ktempy2
integer*8, external :: sumeverybit, sumbeforebit
integer*8, allocatable :: H(:)
double precision :: H_value(1:3000)


!!!******************Change the parameters here****************
SprSizet = ksize*200

!if(mod(ksize,nprocs).eq.0) then
!   nloc = (ksize / nprocs)
!else
!   nloc = (ksize / nprocs)+1
!endif

allocate(H(ksize))

H=0
H_value=0
H_index=0
tempi=0

do jj= 1, ksize  !loop B

        temp=0
        initstatup=listqpt(jj)/2**N;
        initstatdn=mod(listqpt(jj),2**N);

        do kk=0,N-1  !loop A

           !**************************************************************************
           ! KE
           !**************************************************************************

           if(BTEST(initstatup,kk)) then
              l=jj
              if(H(l).eq.0) then
                 temp=temp+1
                 H(l)=temp
                 H_value(temp)=disp(kk)
                 H_index(temp)=l
              else
                 H_value(H(l))=H_value(H(l))+disp(kk)
              endif
           endif
           if(BTEST(initstatdn,kk)) then
              l=jj
              if(H(l).eq.0) then
                 temp=temp+1
                 H(l)=temp
                 H_value(temp)=disp(kk)
                 H_index(temp)=l
              else
                 H_value(H(l))=H_value(H(l))+disp(kk)
              endif
           endif


           !******************************************************************************
           !               \---
           !                \   +   +
           ! Hubbard U:  U  /  c   c    c  c
           !               /--- k+q k'-q k' k
           !               q,k,k'
           !******************************************************************************

           ktempx1=0
           ktempy1=0
           iktemp1=0
           kstat1=1
           ktempx2=0
           ktempy2=0
           iktemp2=0
           kstat2=1

           do kkp=0,N-1
              do qq=0,N-1
                 if(BTEST(initstatup,kk).and.BTEST(initstatdn,kkp)) then
                    ktempx1 = mod(kpoint(kk,1) + kpoint(qq,1),kunit)
                    ktempx2 = mod(kpoint(kkp,1) - kpoint(qq,1)+kunit,kunit)
                    ktempy1 = mod(kpoint(kk,2) + kpoint(qq,2),kunit)
                    ktempy2 = mod(kpoint(kkp,2) - kpoint(qq,2)+kunit,kunit)
                    do while (kstat1.ne.0)
                       if(kpoint(iktemp1,1).eq.ktempx1.and.kpoint(iktemp1,2).eq.ktempy1) then
                         kstat1=0
                       else
                         iktemp1=iktemp1+1
                       endif
                       if(iktemp1.eq.N) then
                          write(*,*) 'iktemp out of bounds, row k+q'
                          write(*,*) kkp, qq
                          stop
                       endif
                    enddo
                    do while (kstat2.ne.0)
                       if(kpoint(iktemp2,1).eq.ktempx2.and.kpoint(iktemp2,2).eq.ktempy2) then
                         kstat2=0
                       else
                         iktemp2=iktemp2+1
                       endif
                       if(iktemp2.eq.N) then
                          write(*,*) 'iktemp out of bounds, row kp-q'
                          write(*,*) jj
                          stop
                       endif
                    enddo
                    afterstatdn=IBCLR(initstatdn,kkp)
                    aftersign=((-1)**sumbeforebit(initstatdn,kkp))
                    afterstatup=IBCLR(initstatup,kk)
                    aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,kk)+sumeverybit(afterstatdn)))
                    if(.not.BTEST(afterstatup,iktemp1).and..not.BTEST(afterstatdn,iktemp2)) then
                       afterstatdn=IBSET(afterstatdn,iktemp2)
                       aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,iktemp2))
                       afterstatup=IBSET(afterstatup,iktemp1)
                       aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,iktemp1)+sumeverybit(afterstatdn)))
                       call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                       if(l.eq.-1) stop

                       if(H(l).eq.0) then
                          temp=temp+1
                          H(l)=temp
                          H_value(temp)=-aftersign*U/dfloat(N)
                          H_index(temp)=l
                       else
                          H_value(H(l))=H_value(H(l))-aftersign*U/dfloat(N)
                       endif

                    endif
                    iktemp1=0
                    kstat1=1
                    iktemp2=0
                    kstat2=1
                 endif
              enddo
           enddo

        enddo !loop A

        do ii=1,temp
           if(H_value(ii) > 0.00001 .or. H_value(ii) < -0.00001) then
              tempi = tempi + 1
              IndexJt(tempi)  = H_index(ii) !still on the old numbering
              IndexIt(tempi)  = jj
              sparseHt(tempi) = H_value(ii)
           endif
           H(H_index(ii))=0
        enddo

enddo   !loop B

SprSizet = tempi

write(*,*) 'spr_size in GenMatrix subroutine =',SprSizet

end subroutine





subroutine SetModelParameters
use ModelParas; use NumOfOrbitalAndElectrons; use BettsCluster
use ConstantParas
implicit none

integer*8 :: ii, kk
t=1.00d0; tt=-0.30d0; ttt=0.0d0
miu=0.0d0
rA=0.0d0; rB=0.1365d0; rC=0.5093;;
Uc=-2.0d0;
Ud=rA+4.0d0*rB+3.0d0*rC
cfs=1.8d0;
Ds= 0.25d0;
Dt= 0.10d0;
E_d=0.0d0;
U=8.0d0;
lambda=0.0d0
E_site(1:N3d)=0.0d0;
E_site(N3d+1:N)= -930.0d0;
do ii=1,N
   E_site(ii+N)=E_site(ii)
end do

!THe hoppping:
phase(:,:)=0.0d0;

!===========================
!Defining the 2p S-O matrix:
!U_SO_2p(3 by 2 by 3 by 2).
!===========================
!U_SO_2p=0.0d0;!Strenght = xi_2p*U_SO_2p

U_SO_2p(1,1,2,1)= dcmplx(0.0d0, 0.5d0);  U_SO_2p(1,1,3,2)= dcmplx(-0.5d0, 0.0d0);
U_SO_2p(1,2,2,2)= dcmplx(0.0d0,-0.5d0);  U_SO_2p(1,2,3,1)= dcmplx(0.5d0, 0.0d0);
U_SO_2p(2,1,1,1)= dcmplx(0.0d0,-0.5d0);  U_SO_2p(2,1,3,2)= dcmplx(0.0d0,-0.5d0);
U_SO_2p(2,2,1,2)= dcmplx(0.0d0,0.5d0);   U_SO_2p(2,2,3,1)= dcmplx(0.0d0,-0.5d0);
U_SO_2p(3,1,1,2)= dcmplx(0.5d0, 0.0d0);  U_SO_2p(3,1,2,2)= dcmplx(0.0d0,0.5d0);
U_SO_2p(3,2,1,1)= dcmplx(-0.5d0, 0.0d0); U_SO_2p(3,2,2,1)= dcmplx(0.0d0,0.5d0);

!!Ti4+:
!xi_2p= 3.78d0;
xi_2p= 13.00d0;
!xi_2p= 0.00d0;

!!Fe2O3:
!!xi_2p=8.199d0;

!!xi_2p=0.0d0;
selectcase(N3d)
case(4)

 if(Ndim.eq.2) then

site(0,1)=2;site(0,2)=2;site(0,3)=1;site(0,4)=1;
site(1,1)=3;site(1,2)=3;site(1,3)=0;site(1,4)=0;
site(2,1)=0;site(2,2)=0;site(2,3)=3;site(2,4)=3;
site(3,1)=1;site(3,2)=1;site(3,3)=2;site(3,4)=2;

  else if(Ndim.eq.1) then

site(0,1)=1;site(0,2)=3;site(0,3)=0;site(0,4)=0;
site(1,1)=2;site(1,2)=0;site(1,3)=1;site(1,4)=1;
site(2,1)=3;site(2,2)=1;site(2,3)=2;site(2,4)=2;
site(3,1)=0;site(3,2)=2;site(3,3)=3;site(3,4)=3;

  endif



case(8)
!************************************************************
! K-points for the 8-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 |     x     |     x
!                 |      (1,3)|      (3,3)
!                 |           |
!                 x-----------x------
!                 |(0,2)      |(2,2)
!                 |           |
!                 |     x     |     x
!                 |      (1,1)|      (3,1)
!                 |           |
!                 x-----------x------
!                 |(0,0)      |(2,0)
!************************************************************
rpoint(0,1)=0;rpoint(0,2)=0;
rpoint(1,1)=2;rpoint(1,2)=0;
rpoint(2,1)=0;rpoint(2,2)=2;
rpoint(3,1)=2;rpoint(3,2)=2;
rpoint(4,1)=1;rpoint(4,2)=1;
rpoint(5,1)=3;rpoint(5,2)=1;
rpoint(6,1)=1;rpoint(6,2)=3;
rpoint(7,1)=3;rpoint(7,2)=3;

runit=4

! site( ,1) up right
! site( ,2) down left
! site( ,3) left up
! site( ,4) right down

if(Ndim.eq.2) then
site(0,1)=2;site(0,2)=7;site(0,3)=3;site(0,4)=6;
site(1,1)=3;site(1,2)=6;site(1,3)=2;site(1,4)=7;
site(2,1)=5;site(2,2)=0;site(2,3)=4;site(2,4)=1;
site(3,1)=4;site(3,2)=1;site(3,3)=5;site(3,4)=0;
site(4,1)=6;site(4,2)=3;site(4,3)=7;site(4,4)=2;
site(5,1)=7;site(5,2)=2;site(5,3)=6;site(5,4)=3;
site(6,1)=1;site(6,2)=4;site(6,3)=0;site(6,4)=5;
site(7,1)=0;site(7,2)=5;site(7,3)=1;site(7,4)=4;

else if(Ndim.eq.1) then
site(0,1)=1;site(0,2)=7;site(0,3)=0;site(0,4)=0;
site(1,1)=2;site(1,2)=0;site(1,3)=1;site(1,4)=1;
site(2,1)=3;site(2,2)=1;site(2,3)=2;site(2,4)=2;
site(3,1)=4;site(3,2)=2;site(3,3)=3;site(3,4)=3;
site(4,1)=5;site(4,2)=3;site(4,3)=4;site(4,4)=4;
site(5,1)=6;site(5,2)=4;site(5,3)=5;site(5,4)=5;
site(6,1)=7;site(6,2)=5;site(6,3)=6;site(6,4)=6;
site(7,1)=0;site(7,2)=6;site(7,3)=7;site(7,4)=7;
endif


case(10)

  if(Ndim.eq.1) then
site(0,1)=1;site(0,2)=9;site(0,3)=0;site(0,4)=0;
site(1,1)=2;site(1,2)=0;site(1,3)=1;site(1,4)=1;
site(2,1)=3;site(2,2)=1;site(2,3)=2;site(2,4)=2;
site(3,1)=4;site(3,2)=2;site(3,3)=3;site(3,4)=3;
site(4,1)=5;site(4,2)=3;site(4,3)=4;site(4,4)=4;
site(5,1)=6;site(5,2)=4;site(5,3)=5;site(5,4)=5;
site(6,1)=7;site(6,2)=5;site(6,3)=6;site(6,4)=6;
site(7,1)=8;site(7,2)=6;site(7,3)=7;site(7,4)=7;
site(8,1)=9;site(8,2)=7;site(8,3)=8;site(8,4)=8;
site(9,1)=0;site(9,2)=8;site(9,3)=9;site(9,4)=9;
  else
write(*,*) 'input Betts cluster dimension is wrong'
  endif


case(12)

  if(Ndim.eq.1) then
site(0,1)=1;site(0,2)=11;site(0,3)=0;site(0,4)=0;
site(1,1)=2;site(1,2)=0;site(1,3)=1;site(1,4)=1;
site(2,1)=3;site(2,2)=1;site(2,3)=2;site(2,4)=2;
site(3,1)=4;site(3,2)=2;site(3,3)=3;site(3,4)=3;
site(4,1)=5;site(4,2)=3;site(4,3)=4;site(4,4)=4;
site(5,1)=6;site(5,2)=4;site(5,3)=5;site(5,4)=5;
site(6,1)=7;site(6,2)=5;site(6,3)=6;site(6,4)=6;
site(7,1)=8;site(7,2)=6;site(7,3)=7;site(7,4)=7;
site(8,1)=9;site(8,2)=7;site(8,3)=8;site(8,4)=8;
site(9,1)=10;site(9,2)=8;site(9,3)=9;site(9,4)=9;
site(10,1)=11;site(10,2)=9;site(10,3)=10;site(10,4)=10;
site(11,1)=0;site(11,2)=10;site(11,3)=11;site(11,4)=11;
  else
write(*,*) 'input Betts cluster dimension is wrong'
  endif

endselect




selectcase(N3d)
case(8)
!************************************************************
! K-points for the 8-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 |     x     |     x
!                 |      (1,3)|      (3,3)
!                 |           |
!                 x-----------x------
!                 |(0,2)      |(2,2)
!                 |           |
!                 |     x     |     x
!                 |      (1,1)|      (3,1)
!                 |           |
!                 x-----------x------
!                 |(0,0)      |(2,0)
!************************************************************
kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=2;kpoint(1,2)=0;
kpoint(2,1)=0;kpoint(2,2)=2;
kpoint(3,1)=2;kpoint(3,2)=2;
kpoint(4,1)=1;kpoint(4,2)=1;
kpoint(5,1)=3;kpoint(5,2)=1;
kpoint(6,1)=1;kpoint(6,2)=3;
kpoint(7,1)=3;kpoint(7,2)=3;

kunit=4

nQpoints=4
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 3
Kmap(4) = 4


case(10)
!************************************************************
! K-points for the 16B-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 p-----x-----x-----x-----x-----p
!                 |                    b  |
!                 |                       |
!                 x     x     p     x     x     x
!                 |  b                    |
!                 |                       |
!                 x     x     x     x     p     x
!                 |              b        |
!                 |                       |
!                 x     p     x     x     x     x
!                 |                       |  b
!                 |                       | 
!                 x     x     x     p     x     x
!                 |        b              |
!                 |                       |
!                 p-----x-----x-----x-----x-----x

!************************************************************

kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=3;kpoint(1,2)=1;
kpoint(2,1)=2;kpoint(2,2)=4;
kpoint(3,1)=5;kpoint(3,2)=5;
kpoint(4,1)=6;kpoint(4,2)=2;
kpoint(5,1)=9;kpoint(5,2)=3;
kpoint(6,1)=8;kpoint(6,2)=6;
kpoint(7,1)=1;kpoint(7,2)=7;
kpoint(8,1)=4;kpoint(8,2)=8;
kpoint(9,1)=7;kpoint(9,2)=9;

kunit=10

nQpoints=4
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 3




case(12)
!************************************************************
! K-points for the 12D-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 x-----x-----x-----x-----x-----x
!                 |                 |      
!                 |                 |      
!                 x     x     x     x     x     x
!                 |                 |      
!                 |                 |      
!                 B     x     x     x     x     x
!                 |        _        |      
!                 |                 |      
!                 x-----x-----x-----B-----x-----x
!                 |                 |      
!                 B                 |      
!                 x     x     x     x     x     x
!                 |                 |        
!                 |        P        |       
!                 x     x     x     x     x     x
!                 |                 B      
!                 |                 |      
!                 B-----x-----x-----x-----x-----x

!************************************************************

kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=0;kpoint(1,2)=4;
kpoint(2,1)=3;kpoint(2,2)=3;
kpoint(3,1)=6;kpoint(3,2)=2;
kpoint(4,1)=9;kpoint(4,2)=1;
kpoint(5,1)=0;kpoint(5,2)=8;
kpoint(6,1)=3;kpoint(6,2)=7;
kpoint(7,1)=6;kpoint(7,2)=6;
kpoint(8,1)=9;kpoint(8,2)=5;
kpoint(9,1)=3;kpoint(9,2)=11;
kpoint(10,1)=6;kpoint(10,2)=10;
kpoint(11,1)=9;kpoint(11,2)=9;

kunit=12

nQpoints=5
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 3
Kmap(5) = 7



case(16)
   selectcase(Nmore)
   case('B')
!************************************************************
! K-points for the 16B-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 x     x     x     x     x
!                 |                       |
!                 |                       |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!                 |(0,3)(1,3) (2,3) (3,3) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,2)(1,2) (2,2) (3,2) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,1)(1,1) (2,1) (3,1) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,0)(1,0) (2,0) (3,0) |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!************************************************************
kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=1;kpoint(1,2)=0;
kpoint(2,1)=2;kpoint(2,2)=0;
kpoint(3,1)=3;kpoint(3,2)=0;
kpoint(4,1)=0;kpoint(4,2)=1;
kpoint(5,1)=1;kpoint(5,2)=1;
kpoint(6,1)=2;kpoint(6,2)=1;
kpoint(7,1)=3;kpoint(7,2)=1;
kpoint(8,1)=0;kpoint(8,2)=2;
kpoint(9,1)=1;kpoint(9,2)=2;
kpoint(10,1)=2;kpoint(10,2)=2;
kpoint(11,1)=3;kpoint(11,2)=2;
kpoint(12,1)=0;kpoint(12,2)=3;
kpoint(13,1)=1;kpoint(13,2)=3;
kpoint(14,1)=2;kpoint(14,2)=3;
kpoint(15,1)=3;kpoint(15,2)=3;

kunit=4



   case('A')
!************************************************************
! K-points for the 16A-site cluster with (+) integer values
! with a factor of (pi/4) removed
!
!                 |  x     x     x     x     x
!                 |                       |
!                 |                       |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!                 |  (1,6)(3,6)(5,6)(7,6) |
!                 |  x     x     x     x  |  x
!                 |                       |
!                 |(0,4)(2,4)(4,4)(6,4)   |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(1,2)(3,2) (5,2) (7,2) |
!                 |  x     x     x     x  |  x
!                 |                       |
!                 |(0,0)(2,0) (4,0) (6,0) |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!************************************************************
kpoint(0,1)=0;   kpoint(0,2)=0;
kpoint(1,1)=2;   kpoint(1,2)=0;
kpoint(2,1)=4;   kpoint(2,2)=0;
kpoint(3,1)=6;   kpoint(3,2)=0;
kpoint(4,1)=1;   kpoint(4,2)=2;
kpoint(5,1)=3;   kpoint(5,2)=2;
kpoint(6,1)=5;   kpoint(6,2)=2;
kpoint(7,1)=7;   kpoint(7,2)=2;
kpoint(8,1)=0;   kpoint(8,2)=4;
kpoint(9,1)=2;   kpoint(9,2)=4;
kpoint(10,1)=4;  kpoint(10,2)=4;
kpoint(11,1)=6;  kpoint(11,2)=4;
kpoint(12,1)=1;  kpoint(12,2)=6;
kpoint(13,1)=3;  kpoint(13,2)=6;
kpoint(14,1)=5;  kpoint(14,2)=6;
kpoint(15,1)=7;  kpoint(15,2)=6;
!**********************************************************************
kunit=8

   endselect
endselect


do kk=0,N3d-1
   kx = real(kpoint(kk,1))*2.0d0*pi/kunit + ibx*2.0*pi/nbx
   ky = real(kpoint(kk,2))*2.0d0*pi/kunit + iby*2.0*pi/nby
   disp(kk) = -2*t*(cos(kx)+cos(ky))&
                -4*tt*(cos(kx)*cos(ky))&
                  -2*ttt*(cos(2*kx)+cos(2*ky))&
                    -miu
   !disp(kk) = -2*t*(cos(real(kpoint(kk,1))*2.0d0*pi/kunit)+cos(real(kpoint(kk,2))*2.0d0*pi/kunit))&
   !             -4*tt*(cos(real(kpoint(kk,1))*2.0d0*pi/kunit)*cos(real(kpoint(kk,2))*2.0d0*pi/kunit))&
   !               -2*ttt*(cos(2*real(kpoint(kk,1))*2.0d0*pi/kunit)+cos(2*real(kpoint(kk,2))*2.0d0*pi/kunit))&
   !                 -miu
   lsf(kk) = (cos(real(kpoint(kk,1))*2.0d0*pi/kunit)+cos(real(kpoint(kk,2))*2.0d0*pi/kunit))*0.5&
             -0.3*cos(real(kpoint(kk,1))*2.0d0*pi/kunit)*cos(real(kpoint(kk,2))*2.0d0*pi/kunit); ! lattice structure factor
   write(*,*) 'kk, cos', kk, cos(real(kpoint(kk,1))*2.0d0*pi/kunit)+cos(real(kpoint(kk,2))*2.0d0*pi/kunit)
   write(*,*) 'kk,disp', kk, disp(kk)
enddo

end subroutine

!*****************************************************


subroutine ChenCmplxContFracExpan(Hsize_3d, SprSize_3d, H_e, groundE, IndexI, IndexJ, sparseH, specX, specY)
use NumOfOrbitalAndElectrons; use ConstantParas; use ScanRegion
implicit none

integer*8, INTENT(IN) :: Hsize_3d, SprSize_3d
integer*8, DIMENSION(SprSize_3d), INTENT(IN) :: IndexI, IndexJ
DOUBLE PRECISION, INTENT(IN) :: groundE
DOUBLE PRECISION, DIMENSION(Hsize_3d), INTENT(IN) :: H_e
DOUBLE PRECISION, DIMENSION(SprSize_3d), INTENT(IN) :: sparseH
DOUBLE PRECISION, DIMENSION(divX+1), INTENT(OUT) :: specX
DOUBLE PRECISION, DIMENSION(divX+1), INTENT(OUT) :: specY

integer*8 :: ii,jj,kk, i, j, k
integer*8 :: Itemp1, Itemp2, Itemp3
integer*8 :: spr_size, iter, iter_num
double precision :: rtemp1, rtemp2, rtemp3
double precision :: b_n
double precision :: time1, time2, time3
!!double precision, allocatable :: sparseH(:)
!!for complex Hamiltonian matrix:
double precision :: a_n
double complex :: ctemp1, ctemp2, ctemp3, ctemp4
double complex :: con_ct1, con_ct2, con_ct3
double complex, allocatable :: phi_n1(:), phi_nm(:), phi_n0(:)
double complex, allocatable :: carray_1(:), carray_2(:), carray_3(:)
double precision, allocatable :: work_an(:)
double precision, allocatable ::work_bn(:)
double precision :: Emin, Emax
integer*8 :: energy_pt, Hsizet

double complex :: Gamma

Gamma=dcmplx(0.0d0,epsilone_CFE)
write(*,*) 'starts to in ChenCmplx'
call CPU_time(time1)

spr_size=SprSize_3d;
Hsizet = Hsize_3d
write(*,*) 'Hsizet = ', Hsizet
!!========================================================================================
!!Part I: Construct the tri-diagonal matrix and obtain the continued faction coefficients:
!!========================================================================================

allocate(carray_1(Hsizet)); carray_1=0.0d0;
allocate(phi_n0(Hsizet)); phi_n0=0.0d0;
allocate(phi_n1(Hsizet)); phi_n1=0.0d0;
allocate(phi_nm(Hsizet)); phi_nm=0.0d0;
write(*,*) 'allocation is done'

!!For CFE: using the O|GS> as the initial starting vector.
phi_n0(1:Hsizet) = H_e(1:Hsizet)
iter_num=niter_CFE;
allocate(work_an(iter_num)); work_an=0.0d0;
allocate(work_bn(iter_num)); work_bn=0.0d0;

!!Specify the order of continued fraction expansions.
!!Note: Usually iter_num=50 already gives truncation error < 10-6;
!!      For iter_num > 100 orthogonality may leak and result in "NaN".
!!************************************************

!phi = initial
!do
!  c(:)=c(:)*Matrix
!  ct = c(:)*phi

iter=0; j=0;
do while(j.eq.0)

   call CPU_time(time2)

   iter=iter+1

   !!Constructing H|phi_n0>: later want to do in one step so only 3 vectors are needed.
   carray_1=0.0d0;
   do ii=1, spr_size
      carray_1(IndexI(ii))=  carray_1(IndexI(ii))+sparseH(ii)*phi_n0(IndexJ(ii));
   enddo

   !!Calculate <phi_n0|H|phi_n0>:
   ctemp1=0.0d0;
   rtemp1=0.0d0;
   do ii=1, Hsizet
      ctemp1=ctemp1+conjg(carray_1(ii))*phi_n0(ii)
      rtemp1=rtemp1+conjg(phi_n0(ii))*phi_n0(ii)
   enddo

   a_n=ctemp1/rtemp1;


   if(iter.eq.1) b_n=0.0d0;
   phi_n1=carray_1-a_n*phi_n0-b_n*phi_nm;

   !write(*,*) iter, a_n, b_n, rtemp1
   work_an(iter)=a_n;
   work_bn(iter)=b_n;

   rtemp2=0.0d0;
   do ii=1, Hsizet
      rtemp2=rtemp2+phi_n1(ii)*conjg(phi_n1(ii))
   enddo
   b_n=rtemp2/rtemp1;

   phi_nm=phi_n0
   phi_n0=phi_n1;
   phi_n1=0.0d0;

   if(iter.eq.iter_num) j=3

   call CPU_time(time3)
   do ii=1,Hsizet
      phi_n0(ii)=phi_n0(ii)/sqrt(rtemp2)
      phi_nm(ii)=phi_nm(ii)/sqrt(rtemp2)
   enddo


enddo

!!=======================================================
!!Part II: Calculating the continued fraction expansion:
!!=======================================================

!!Note and Suggestion:
!!Specify Emin to be the ground energy;
!!Emax-Emin would be the energy range of interest (for spectra).

!Emin=   -9.2832;
Emin = startX
Emax = endX
!!How many points in the spectra; 1001 is goodl
!!this in principle does not affect the time of the code.
energy_pt=divX+1;

!!The Lorenztian broadening of the spectra:
!Gamma=(0.0d0, 0.3d0);

carray_1=0.0d0;
rtemp1=0.0d0;

!open(unit=123, file='OGS.txt',status='old')
!do ii=1, Hsizet
!   read(123,*) carray_1(ii)
!   rtemp1=rtemp1+carray_1(ii)*conjg(carray_1(ii))
!enddo
carray_1(1:Hsizet) = H_e(1:Hsizet)
do ii=1, Hsizet
   rtemp1 = rtemp1 + carray_1(ii)*conjg(carray_1(ii))
enddo

!open(unit=1111,file='check_CFE.txt')


!do ii=1, energy_pt
do eneX=0,divX
   !rtemp2= Emin+ (Emax-Emin)*(ii-1)/(energy_pt-1)
   !ctemp1= rtemp2 + Gamma
   rtemp2 = dble(eneX)/dble(divX)*(endX-startX)+startX
   ctemp1 = rtemp2 + groundE + Gamma
   j=0
   do jj=iter_num, 2, -1

      j=j+1
      if(j.eq.1) ctemp4=0.0d0;

      ctemp3=0.0d0;
      ctemp3= work_bn(jj)/(ctemp1-work_an(jj) -ctemp4);
      ctemp4= ctemp3;

   enddo
   con_ct1=rtemp1/(ctemp1-work_an(1)-ctemp4)


   j=0
   do jj=iter_num-5, 2, -1

      j=j+1
      if(j.eq.1) ctemp4=0.0d0;

      ctemp3=0.0d0;
      ctemp3= work_bn(jj)/(ctemp1-work_an(jj) -ctemp4);
      ctemp4= ctemp3;

   enddo
   con_ct2=rtemp1/(ctemp1-work_an(1)-ctemp4)

   write(*,*) rtemp2, imag(con_ct1)/(-pi), imag(con_ct2)/(-pi)
   specX(eneX+1) = rtemp2
   specY(eneX+1) = imag(con_ct1)/(-pi)

enddo
!close(1111)


call CPU_time(time2)
write(*,*) 'Complex CFE (secs):', time2-time1
deallocate(carray_1, phi_n0, phi_n1, phi_nm, work_an, work_bn)

end subroutine



subroutine ContFracExpan(Hsize, SprSize, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
use ScanRegion; use ConstantParas
implicit none
INTEGER*8, INTENT(IN) :: Hsize, SprSize
INTEGER*8, DIMENSION(SprSize), INTENT(IN) :: IndexI, IndexJ
DOUBLE PRECISION, INTENT(IN) :: E_0
DOUBLE PRECISION, DIMENSION(Hsize),INTENT(IN) :: H_0
DOUBLE PRECISION, DIMENSION(SprSize), INTENT(IN) :: sparseH
DOUBLE PRECISION, DIMENSION(divX+1), INTENT(OUT) :: specX, specY

integer*8 :: niter
integer*8 :: ii,jj,kk
double precision :: mysum, factor
double precision, allocatable :: alpha(:), betha(:)
double precision, allocatable:: phi(:), phil(:),phip(:), phipp(:)
double complex :: z
double complex :: Intensity(divX+1)

   niter = 200

   allocate(phi(Hsize), phipp(Hsize))
   allocate(phil(Hsize), phip(Hsize))
   allocate(alpha(niter),betha(niter))

   phi(1:Hsize) = H_0(1:Hsize)

   mysum = 0.0;
   do jj=1, Hsize
        mysum = mysum + phi(jj)*phi(jj)
   enddo
   factor = mysum
   do jj=1, Hsize
        phi(jj) = phi(jj)/sqrt(mysum)
   enddo
   write(*,*) 'sum=',mysum
   !--------------------------------------------------

   do ii = 1, niter

!       According to different paper, the index of betha 
!       is different 
!       here, betha runs from 2 to niter
!             alpha runs from 1 to niter

      phip=0.0d0
      do kk=1,SprSize
         phip(IndexI(kk))=phip(IndexI(kk))+sparseH(kk)*phi(IndexJ(kk))
      enddo

      if(ii.ne.1) then
         do kk=1, Hsize
           phip(kk) = phip(kk) - betha(ii)*phil(kk)
         enddo
      endif

      alpha(ii)=0.0d0
      do kk=1, Hsize
        alpha(ii) = alpha(ii) + phip(kk)*phi(kk)
      enddo
      write(*,*) 'alpha', alpha(ii)

      do kk=1, Hsize
        phipp(kk) = phip(kk) - alpha(ii)*phi(kk)
      enddo
   if(ii.ne.niter) then
      betha(ii+1) = 0.0d0
      do kk=1, Hsize
        betha(ii+1) = betha(ii+1) + phipp(kk)**2
      enddo
      betha(ii+1) = sqrt(betha(ii+1))
      write(*,*) 'betha', betha(ii+1)
      do kk=1, Hsize
        phil(kk) = phi(kk)
      enddo
      do kk=1, Hsize
        phi(kk) = phipp(kk) / betha(ii+1)
      enddo
   endif

   enddo
   deallocate(phi, phil, phip, phipp)

do ii=0, divX
   z = CMPLX(dble(ii)/divX*(endX-startX)+startX+E_0,epsilone_CFE)
   Intensity(ii+1)=z-alpha(niter)
   do jj=1,niter-1
      Intensity(ii+1)=z-alpha(niter-jj)-betha(niter-jj+1)**2/Intensity(ii+1)
   enddo
   specX(ii+1)=dble(ii)/divX*(endX-startX)+startX
   specY(ii+1)=-1*1/pi*AIMAG(factor/Intensity(ii+1))
enddo

deallocate(alpha,betha)

end subroutine




!********************************************************

subroutine BinarySearch(listqpt,ksize,statupdn,l)
implicit none
INTEGER*8, INTENT(IN) :: ksize
INTEGER*8, DIMENSION(1:ksize), INTENT(IN) :: listqpt
INTEGER*8, INTENT(IN) :: statupdn
INTEGER*8, INTENT(OUT):: l

integer*8:: head, tail, middle

head=1; tail=ksize; middle=(head+tail)/2

do while((listqpt(middle).ne.statupdn).and.(head.le.tail))
   if(statupdn.gt.listqpt(middle)) then
      head = middle+1
   else
      tail = middle-1
   endif
   middle = (head+tail)/2
enddo

if(listqpt(middle).eq.statupdn) then
   l=middle;
else
   l=-1;
endif

end

!*****************************************************

integer*8 function factorial(tempn)
implicit none
integer*8:: tempn,tempn1

factorial=1
tempn1=tempn
do while (tempn1.ne.0)
        factorial=factorial*tempn1;
        tempn1=tempn1-1;
enddo
return
end function

!********************************************************

INTEGER*8 function sumeverybit(tempa)
implicit none
integer*8::tempa,tempa1,tempsum

tempsum=0
tempa1=tempa
do while (tempa1.ne.0)
        tempsum=tempsum+mod(tempa1,2)
        tempa1=ISHFT(tempa1,-1)
enddo
sumeverybit=tempsum
return
end function

!********************************************************

integer*8 function sumbeforebit(temps,ks)
implicit none
integer*8::temps,temps1,is,ks

sumbeforebit=0;

if(ks.eq.0) return

temps1=temps
do is=0,ks-1
        if(BTEST(temps1,is).eqv..true.) then
                sumbeforebit=sumbeforebit+1
        endif
enddo
return
end function

!**********************************************************
integer*8 function ksubtract(kk, kkp)
use BettsCluster; use NumOfOrbitalAndElectrons
implicit none
integer*8 :: kk, kkp
integer*8 :: iktemp, kstat, ktempx, ktempy
integer*8 :: px, py

px = kpoint(kkp,1)
py = kpoint(kkp,2)
      ! code to get kkp
      iktemp=0
      kstat=1
      ktempx = mod(kpoint(kk,1)-px+2*kunit,kunit)
      ktempy = mod(kpoint(kk,2)-py+2*kunit,kunit)
      do while (kstat.ne.0)
         if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
            kstat=0
         else
            iktemp=iktemp+1
         endif
         if(iktemp.eq.N) then
            write(*,*) 'iktemp out of bounds, spin up'
            stop
         endif
      enddo
      ksubtract=iktemp
end function


!**********************************************************
integer*8 function ksum(kk, kkp)
use BettsCluster; use NumOfOrbitalAndElectrons
implicit none
integer*8 :: kk, kkp
integer*8 :: iktemp, kstat, ktempx, ktempy
integer*8 :: px, py

px = kpoint(kkp,1)
py = kpoint(kkp,2)
      ! code to get kkp
      iktemp=0
      kstat=1
      ktempx = mod(kpoint(kk,1)+px+kunit,kunit)
      ktempy = mod(kpoint(kk,2)+py+kunit,kunit)
      do while (kstat.ne.0)
         if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
            kstat=0
         else
            iktemp=iktemp+1
         endif
         if(iktemp.eq.N) then
            write(*,*) 'iktemp out of bounds, spin up'
            stop
         endif
      enddo
      ksum=iktemp
end function

