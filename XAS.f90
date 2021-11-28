module NumOfOrbitalAndElectrons
   integer*8 :: N3d, N2p, N, Ndim
   integer*8 :: nup, ndn, ntol
endmodule NumOfOrbitalAndElectrons

module ModelParas
   double precision :: U(1:30), E_site(1:30)
   double precision :: Uc, t, tt
   double precision :: Ud, Up, UQ, E_d, E_pxy, E_pz, tpd, tpp, tpdz, tpzdz, tppz
   double precision :: U_eg, U_t2g, Upp, E_eg, E_t2g, cfs
   double precision :: rA, rB, rC0, F_0, F_2, G_1, G_3, Dq, Ds, Dt
   double precision :: U_rest(5,5,5,5), U_pddp(3,5,5,3), U_dpdp(5,3,5,3)
   double precision :: U_ext(5,5), J_ext(5,5)
   double precision :: phase(30,30)
   double precision :: xi_2p
   double complex :: U_SO_2p(3,2,3,2)
   double complex :: prefac(4) ! for twisted average boundary condition
   integer*8 :: nbx, nby, ibx, iby ! number of boundaries in x and y direction 
endmodule ModelParas

module MPIParas
   integer  ::   comm, myid, nprocs
   integer  ::   source, dest, tag, ierr, rc
   integer*8 ::  nloc, localstart, localend
   integer*8 ::  nloc_3d0, localstart_3d0, localend_3d0
   integer*8 ::  nloc_3d1, localstart_3d1, localend_3d1
endmodule MPIParas

module ConstantParas
   integer*8 :: HsizeEst=15077600, SprSizeEst=1007760*50
   integer*8 :: Hsize_3dEst=880400, SprSize_3dEst=63504*50
   integer*8 :: NNN=1000, niter_CL=100, niter_CFE = 100, niter_CG=100
   double precision :: tol_CG=0.0001, tol_CL=0.0001
   double precision :: pi=3.1415926535, sqrthalf = 0.707106781
   double precision :: epsilone_CFE=0.3d0,  epsilone_CG = 1.0d0
endmodule ConstantParas

module BettsCluster
   integer*8 :: rpoint(0:15,2), runit
   integer*8 :: site(0:15,4)
   integer*8 :: nQpoints
   double complex :: QPhase(0:15,0:8)
   double precision :: Qpoint(0:15,2), QA(0:15,2)
endmodule BettsCluster

module ScanRegion
   integer*8 :: divX, divY
   integer*8:: eneX, eneY
   double precision :: startX, endX, startY, endY
endmodule ScanRegion

module PolarParas
   double precision :: k_theta_in(0:15), k_theta_out(0:15)
   double precision :: k_phi_in(0:15), k_phi_out(0:15)
   double precision :: p_theta_in(0:15), p_phi_in(0:15)
   double precision :: p_theta_out(0:15), p_phi_out(0:15)
   double precision :: Gn_in(3,0:15), Gn_out(3,0:15)
   double precision :: k0 !proportional to Ein, k0=0.58 for Ein=930eV
   character*1 :: PorS_in, PorS_out
endmodule PolarParas


Program Main
use NumOfOrbitalAndElectrons; use ModelParas; use BettsCluster
use ConstantParas; use ScanRegion; use PolarParas; use MPIParas
implicit none
include 'mpif.h'

!CHARACTER(LEN=30) :: Format
!      call MPI_INIT( ierr )
!      comm = MPI_COMM_WORLD
!      call MPI_COMM_RANK( comm, myid, ierr )
!      call MPI_COMM_SIZE( comm, nprocs, ierr )
integer  ::   status(MPI_STATUS_SIZE)

   integer*8 :: Hsize, Nkept, SprSize
   integer*8 :: Hsize_3d_spin0, SprSize_3d_spin0
   integer*8 :: Hsize_3d_spin1, SprSize_3d_spin1
   integer*8 :: SprSize_O, SprSize_OD

integer*8 :: ii, jj, kk, mm, nn, popo, toto
integer*8 :: i, j, k, m, itempuse, tempsign
integer*8 :: Ediv, Eloss_div
integer*8 :: flag, nqq
integer*8 :: thread_num
integer*8 :: time
integer*8 :: CHsite, pp
integer*8, allocatable :: Hsp(:), Hsp_3d_spin0(:), Hsp_3d_spin1(:)
integer*8, external :: factorial
integer*8, external :: obcount
integer*8, external :: omp_get_num_procs
integer*8, allocatable :: IndexI(:), IndexJ(:)
integer*8, allocatable :: IndexI_3d_spin0(:), IndexJ_3d_spin0(:)
integer*8, allocatable :: IndexI_3d_spin1(:), IndexJ_3d_spin1(:)
integer*8, allocatable :: IndexI_O(:), IndexJ_O(:)
integer*8, allocatable :: IndexI_OD(:), IndexJ_OD(:)
integer*8, allocatable :: To3dMap(:), From3dMap(:)
integer*8 :: spin_tol
double precision :: a, maxE, minE, maxE_loss, minE_loss
double precision :: tempuse1, tempuse2, tempuse3
double precision :: time1, time2, rtemp_1
double precision, external :: Lorentzian
double precision :: E_0_3d
double precision, allocatable :: E_0(:), E_e(:)
double precision, allocatable :: specX(:)
double precision, allocatable :: specRIXS(:)
double precision, allocatable :: w(:), wloss(:), RIXS(:,:)
double complex :: gamma, tempusecomp, z
double complex, external :: VectorProductComplex
double complex, allocatable :: sparseH(:)
double complex, allocatable :: sparseH_3d_spin0(:)
double complex, allocatable :: sparseH_3d_spin1(:)
double complex, allocatable :: sparseH_O(:)
double complex, allocatable :: sparseH_OD(:)
double complex, allocatable :: H_0(:,:), H_e(:,:), H_f(:,:), H_i(:,:), H_id(:,:)
double complex, allocatable :: H_g(:), H_g_OD(:,:), InProd_XAS(:)
double complex, allocatable :: RIXS_Matrix(:,:,:), RIXS_Matrix2(:,:,:)
double complex, allocatable :: ctempv(:), ctempv0(:), tempv0_3d(:)
double complex, allocatable :: pvector(:), H_0_3d(:)
character*2 :: twochar
character*1 :: resf

!Lapack subroutine 
integer*8 :: LWORK
integer*8 :: info
DOUBLE PRECISION, allocatable ::  RWORK(:)
double complex, allocatable :: WORK(:)


      call MPI_INIT( ierr )
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
      if(myid.eq.0) write(*,*) '# of processors = ', nprocs
!      pp = myid

!============== Read input file ==============

open(unit=10101,file='input_fullED',Status='old');
read(10101,'(3I8)') N3d, N2p, Ndim
read(10101,'(2I8)') nup, ndn
read(10101,'(I8)') pp
read(10101,'(2F8.2)') maxE_loss,minE_loss
read(10101,'(2F8.2)') maxE,minE
read(10101,'(2I8)') Ediv,Eloss_div
read(10101,'(2A1)') PorS_in, PorS_out
read(10101,'(3I8)') niter_CL, niter_CFE, niter_CG
read(10101,'(2F8.4)') tol_CL, tol_CG
read(10101,'(2F8.2)') epsilone_CG, epsilone_CFE
close(10101)

if(myid.eq.0) then

write(*,*) ''
write(*,*) '     N3d     N2p     Ndim '
write(*,'(3I8)') N3d, N2p, Ndim
write(*,*) ''
write(*,*) '     nup     ndn'
write(*,'(2I8)') nup, ndn
write(*,*) ''
write(*,*) '   pp'
write(*,'(I8)') pp
write(*,*) ''
write(*,*) 'maxE_loss minE_loss'
write(*,'(2F8.2)') maxE_loss,minE_loss
write(*,*) ''
write(*,*) '    maxE    minE'
write(*,'(2F8.2)') maxE,minE
write(*,*) ''
write(*,*) '    Ediv Eloss_div'
write(*,'(2I8)') Ediv,Eloss_div
write(*,*) ''
write(*,*) 'PorS_in PorS_out'
write(*,'(2A1)') PorS_in, PorS_out
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

endif

N=N3d+N2p
ntol = nup+ndn
startX = minE_loss; endX = maxE_loss; divX = Eloss_div-1
startY = minE; endY = maxE; divY = Ediv-1

!============== The Step 0: preparation ==============

call SetQPhase
call SetModelParameters
call SetPolarParas

allocate(Hsp(HsizeEst))
allocate(Hsp_3d_spin0(Hsize_3dEst))
allocate(Hsp_3d_spin1(Hsize_3dEst))
Hsp=0; 
Hsp_3d_spin0=0
Hsp_3d_spin1=0
call GenHsp2pfull(int(0,8), Hsp_3d_spin0, Hsize_3d_spin0)
call GenHsp2pfull(int(2,8), Hsp_3d_spin1, Hsize_3d_spin1)
call GenHsp2p1hole(Hsp, Hsize)
write(*,*) 'Hsize_3d_spin0 = ', Hsize_3d_spin0
write(*,*) 'Hsize_3d_spin1 = ', Hsize_3d_spin1
write(*,*) 'Hsize = ', Hsize

SprSizeEst = Hsize*80
SprSize_3dEst = Hsize_3d_spin0*60
write(*,*) 'SprSizeEst', SprSizeEst

!============== The Step 1: Diagonalize 3d system ==============

CHsite=0
allocate(H_0_3d(Hsize_3dEst))
E_0_3d=0.0d0; 
H_0_3d=dcmplx(0.0d0,0.0d0); 

nbx = int(sqrt(real(nprocs)))
nby = int(sqrt(real(nprocs)))
ibx = myid / nby
iby = mod(myid,nby)
write(*,*) 'ibx and iby', ibx, iby

   prefac(1) = dcmplx(cos(ibx*2*pi/nbx), sin(ibx*2*pi/nbx))
   prefac(2) = dcmplx(cos(ibx*2*pi/nbx), -1*sin(ibx*2*pi/nbx))
   prefac(3) = dcmplx(cos(iby*2*pi/nby), sin(iby*2*pi/nby))
   prefac(4) = dcmplx(cos(iby*2*pi/nby), -1*sin(iby*2*pi/nby))

   allocate(IndexI_3d_spin0(SprSize_3dEst))
   allocate(IndexJ_3d_spin0(SprSize_3dEst))
   allocate(sparseH_3d_spin0(SprSize_3dEst))
   SprSize_3d_spin0=SprSize_3dEst
   call GenMatrix(CHsite, Hsp_3d_spin0, Hsize_3d_spin0, &
              IndexI_3d_spin0, IndexJ_3d_spin0, sparseH_3d_spin0, SprSize_3d_spin0) 
   H_0_3d=dcmplx(0.0d0,0.0d0); 
   call ARPACK(Hsize_3d_spin0, SprSize_3d_spin0, IndexI_3d_spin0, &
               IndexJ_3d_spin0, sparseH_3d_spin0, E_0_3d, H_0_3d)
   deallocate(IndexI_3d_spin0)
   deallocate(IndexJ_3d_spin0)
   deallocate(sparseH_3d_spin0)

!===========   The Step 4: forming RIXS using BiCGS   ==============

   allocate(specX(divX+1))
   allocate(specRIXS(divX+1))
   specX=0.0d0
   specRIXS=0.0d0

   allocate(IndexI(SprSizeEst))
   allocate(IndexJ(SprSizeEst))
   allocate(sparseH(SprSizeEst))
   SprSize=SprSizeEst

   allocate(ctempv(HsizeEst))
   ctempv=dcmplx(0.0d0,0.0d0)
   allocate(ctempv0(Hsize_3dEst))
   ctempv0=dcmplx(0.0d0,0.0d0)

   do CHsite=0,0
      !tempv = (z-H)**-1 * tempv

  !-------  Matrix with core-hole
  ! 
      call GenMatrix(CHsite, Hsp, Hsize, IndexI, IndexJ, sparseH, SprSize) 
      call FillInSpinOrbitCoupling(Hsp, Hsize, IndexI, IndexJ, sparseH, SprSize)
  ! 
  !-------  Matrix with core-hole done

   do pp=0, 0
      call GenOStatesnew(CHsite, H_0_3d, ctempv, Hsp_3d_spin0, Hsp, &
           Hsize_3d_spin0, Hsize, Gn_in(:,pp))
      deallocate(Hsp)
      call CPU_time(time1)
      call ChenCmplxContFracExpan(Hsize, SprSize, &
           ctempv, E_0_3d, IndexI, IndexJ, sparseH, specX, specRIXS)
      call CPU_time(time2)

      write(twochar,'(i2)') ibx*10+iby
      open (111,file='XAS-'//trim(twochar)//'.dat', access='append')
      do mm=1, divX+1
         write(111,*) specX(mm), specRIXS(mm)
      enddo
      close(111)

   enddo !pp
   enddo !CHsite


 !MPI  enddo !pp

!enddo !eneY

2808 call MPI_Finalize(ierr)

end program

!***************************************

subroutine PrintSprMtrx(IndexIt, IndexJt, sparseHt, SprSizet, filename)
implicit none

INTEGER*8, INTENT(IN) :: SprSizet
INTEGER*8, DIMENSION(SprSizet), INTENT(IN) :: IndexIt, IndexJt
integer*8, INTENT(IN) :: filename
double precision, DIMENSION(SprSizet), INTENT(IN) :: sparseHt
character :: FILENAME_1*100, FILENAME_2*100, FILENAME_3*100
integer*8 :: ii

write(unit=FILENAME_1,fmt="('sparseH_',i3.3,'.txt')") filename
write(unit=FILENAME_2,fmt= "('IndexI_',i3.3,'.txt')") filename
write(unit=FILENAME_3,fmt= "('IndexJ_',i3.3,'.txt')") filename
open(unit=12, file=FILENAME_1, status='unknown')
open(unit=13, file=FILENAME_2, status='unknown')
open(unit=14, file=FILENAME_3, status='unknown')
do ii=1, SprSizet
   write(12, *) sparseHt(ii)
   write(13, *) IndexIt(ii)
   write(14, *) IndexJt(ii)
enddo
close(12)
close(13)
close(14)

end subroutine

!****************************************

subroutine SetModelParameters
use ModelParas; use NumOfOrbitalAndElectrons; use BettsCluster
implicit none

integer*8 :: ii
!The Hamiltonain:
!include 'Coulombint.f90'
!define U_ext(:,:), J_ext(:,:)
t=0.40d0; tt=-0.12d0
rA=0.0d0; rB=0.1365d0; rC0=0.5093;;
Uc=-2.4d0;
U(1:N3d)=3.2d0;
!t=0.35d0; tt=-1.2d0; U=3.5d0; PRL 80, 4245, 1998 C. Kim
!Uc = 3.5/2 Physica B 284-288(2000) 457-458 K. Tsutsui 


!include 'multiplet.f90'
!define E_site(:), U(:), phase(:)
Ud=rA+4.0d0*rB+3.0d0*rC0
cfs=1.8d0;
Ds= 0.25d0;
Dt= 0.10d0;
E_d=0.0d0;
!U(1:N3d)=Ud;
E_site=0.0d0;
E_site(1:N3d)=0.0d0;
!!Oh Symmetry:
!E_site(1:2)=E_d+0.6d0*cfs;
!E_site(3:5)=E_d-0.4d0*cfs;

!!Td Symmetry:
!!E_site(1:2)=E_d-0.6d0*cfs;
!!E_site(3:5)=E_d+0.4d0*cfs;

!!D4h c-axis distortion:
!E_site(1)=E_site(1) -2.0d0*Ds-6.0d0*Dt
!E_site(2)=E_site(2) +2.0d0*Ds-1.0d0*Dt
!E_site(3)=E_site(3) +2.0d0*Ds-1.0d0*Dt
!E_site(4)=E_site(4) -1.0d0*Ds+4.0d0*Dt

!!test:
!E_site(1)=E_site(1)+0.8d0;

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

!Not really used in this code.....
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
 
write(*,*) 'It is a 10-site cluster'
  
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

 else if(Ndim.eq.2) then
!10 B (3,1) (-1,3)
! site( ,1) right
! site( ,2) left
! site( ,3) up
! site( ,4) down

site(0,1)=7;site(0,2)=3;site(0,3)=1;site(0,4)=9;
site(1,1)=2;site(1,2)=6;site(1,3)=4;site(1,4)=0;
site(2,1)=3;site(2,2)=1;site(2,3)=5;site(2,4)=7;
site(3,1)=0;site(3,2)=2;site(3,3)=6;site(3,4)=8;
site(4,1)=5;site(4,2)=9;site(4,3)=7;site(4,4)=1;
site(5,1)=6;site(5,2)=4;site(5,3)=8;site(5,4)=2;
site(6,1)=1;site(6,2)=5;site(6,3)=9;site(6,4)=3;
site(7,1)=8;site(7,2)=0;site(7,3)=2;site(7,4)=4;
site(8,1)=9;site(8,2)=7;site(8,3)=3;site(8,4)=5;
site(9,1)=4;site(9,2)=8;site(9,3)=0;site(9,4)=6;

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

  else if(Ndim.eq.2) then
!12 D (1,3) (4,0)
! site( ,1) right
! site( ,2) left
! site( ,3) up
! site( ,4) down

site(0,1)=1;site(0,2)=3;site(0,3)=7;site(0,4)=8;
site(1,1)=2;site(1,2)=0;site(1,3)=4;site(1,4)=9;
site(2,1)=3;site(2,2)=1;site(2,3)=5;site(2,4)=10;
site(3,1)=0;site(3,2)=2;site(3,3)=6;site(3,4)=11;
site(4,1)=5;site(4,2)=7;site(4,3)=8;site(4,4)=1;
site(5,1)=6;site(5,2)=4;site(5,3)=9;site(5,4)=2;
site(6,1)=7;site(6,2)=5;site(6,3)=10;site(6,4)=3;
site(7,1)=4;site(7,2)=6;site(7,3)=11;site(7,4)=0;
site(8,1)=9;site(8,2)=11;site(8,3)=0;site(8,4)=4;
site(9,1)=10;site(9,2)=8;site(9,3)=1;site(9,4)=5;
site(10,1)=11;site(10,2)=9;site(10,3)=2;site(10,4)=6;
site(11,1)=8;site(11,2)=10;site(11,3)=3;site(11,4)=7;

  else
write(*,*) 'input Betts cluster dimension is wrong'
  endif



case(14)

write(*,*) 'It is a 14-site cluster'

  if(Ndim.eq.1) then
site(0,1)=1;site(0,2)=13;site(0,3)=0;site(0,4)=0;
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
site(11,1)=12;site(11,2)=10;site(11,3)=11;site(11,4)=11;
site(12,1)=13;site(12,2)=11;site(12,3)=12;site(12,4)=12;
site(13,1)=0;site(13,2)=12;site(13,3)=13;site(13,4)=13;

  else if(Ndim.eq.2) then
!14 B (4,2) (-1,3)
! site( ,1) right
! site( ,2) left
! site( ,3) up
! site( ,4) down

site(0,1)=8;site(0,2)=7;site(0,3)=1;site(0,4)=13;
site(1,1)=2;site(1,2)=11;site(1,3)=4;site(1,4)=0;
site(2,1)=3;site(2,2)=1;site(2,3)=5;site(2,4)=8;
site(3,1)=12;site(3,2)=2;site(3,3)=6;site(3,4)=9;
site(4,1)=5;site(4,2)=13;site(4,3)=8;site(4,4)=1;
site(5,1)=6;site(5,2)=4;site(5,3)=9;site(5,4)=2;
site(6,1)=7;site(6,2)=5;site(6,3)=10;site(6,4)=3;
site(7,1)=0;site(7,2)=6;site(7,3)=11;site(7,4)=12;
site(8,1)=9;site(8,2)=0;site(8,3)=2;site(8,4)=4;
site(9,1)=10;site(9,2)=8;site(9,3)=3;site(9,4)=5;
site(10,1)=11;site(10,2)=9;site(10,3)=12;site(10,4)=6;
site(11,1)=1;site(11,2)=10;site(11,3)=13;site(11,4)=7;
site(12,1)=13;site(12,2)=3;site(12,3)=7;site(12,4)=10;
site(13,1)=4;site(13,2)=12;site(13,3)=0;site(13,4)=11;
  else
write(*,*) 'input Betts cluster dimension is wrong'
endif


endselect

end


!******************************************************************************

subroutine FillInSpinOrbitCoupling(Hsp, Hsizet, &
                                   IndexI, IndexJ, sparseH, SprSizet)
use NumOfOrbitalAndElectrons; use ModelParas; use ConstantParas
implicit none

!integer*8, INTENT(IN) :: tag
integer*8, INTENT(IN) :: Hsizet
integer*8, DIMENSION(Hsizet), INTENT(IN) :: Hsp
integer*8, DIMENSION(SprSizeEst), INTENT(INOUT) :: IndexI, IndexJ
double complex, DIMENSION(SprSizeEst), INTENT(INOUT) :: sparseH
integer*8, INTENT(INOUT) :: SprSizet
!double complex, DIMENSION(Hsize,Hsize), INTENT(INOUT) :: H

integer*8 :: ii,jj, tempuse, tempsign
integer*8 :: tempuse_1, tempcount, sumbeforebit, tempsign_1
integer*8 :: tempuse_0, tempsign_2, tempuse_2
integer*8 :: kk, popo, ncount, obcount
double precision:: time1
!**************************************************************************************  
!!!Read in the matrix without SO couplings:
!**************************************************************************************  
!include 'H_SO_2p.f90'

!***************************************************************************************************************
!01/26/2010. Include The 2p spin-orbital interactions for L2,3 absorption edges:
!the interaction matrix is given directly below:
!***************************************************************************************************************
ncount = SprSizet

do ii=1, Hsizet
!do ii=lstart, lend

!if(mod(ii,10000).eq.0) then
!   call CPU_time(time1)
!   write(*,*) 'time to run ii ', ii, time1
!endif

   !1. c^\dagger_{m,dn} c_{mm', dn}
   do jj=1+N3d, N2p+N3d
      if (BTEST(Hsp(ii), jj-1)) then
         tempuse_1= IBCLR(Hsp(ii), jj-1)
         tempcount=obcount(Hsp(ii), jj-1+1)
         tempsign_1= (-1)**(tempcount)
      else
         cycle
      endif

      do kk=1+N3d, N2p+N3d
         if(abs(U_SO_2p(kk-N3d, 1, jj-N3d, 1)).lt.0.0001) cycle
         if (BTEST(tempuse_1, kk-1)) cycle
         tempuse_2=IBSET(tempuse_1, kk-1)
         tempcount=obcount(tempuse_1,kk-1+1)
         tempsign_2=(-1)**(tempcount);
         tempsign=tempsign_1*tempsign_2
         tempuse=tempuse_2

         call BinarySearch(Hsp,Hsizet,tempuse,popo)
         if(popo.eq.-1) then
             write(*,'(B30)') tempuse  
             write(*,*) 'A' 
             stop
         endif
         !do popo=1, Hsize
         !   if (tempuse.eq.Hsp(popo)) then
              ! H(ii,popo)=H(ii,popo)+dcmplx(xi_2p*U_SO_2p(kk-N3d, 1, jj-N3d, 1)*tempsign)
               ncount = ncount+1
               IndexI(ncount) = ii
               IndexJ(ncount) = popo
               sparseH(ncount) = dcmplx(xi_2p*U_SO_2p(kk-N3d, 1, jj-N3d, 1)*tempsign)
         !   endif
         !enddo

      enddo!kk
   enddo!jj

   !2. c^\dagger_{m,up} c_{mm', up}
   do jj=1+N3d+N, N2p+N3d+N
      if (BTEST(Hsp(ii), jj-1)) then
         tempuse_1= IBCLR(Hsp(ii), jj-1)
         tempcount=obcount(Hsp(ii), jj-1+1)
         tempsign_1= (-1)**(tempcount)
      else
         cycle
      endif

      do kk=1+N3d+N, N2p+N3d+N
         if(abs(U_SO_2p(kk-N3d-N, 2, jj-N3d-N, 2)).lt.0.00001) cycle
         if (BTEST(tempuse_1, kk-1)) cycle
         tempuse_2=IBSET(tempuse_1, kk-1)
         tempcount=obcount(tempuse_1,kk-1+1)
         tempsign_2=(-1)**(tempcount);
         tempsign=tempsign_1*tempsign_2
         tempuse=tempuse_2

         call BinarySearch(Hsp,Hsizet,tempuse,popo)
         if(popo.eq.-1) then
             write(*,*) tempuse, 'B' 
             stop
         endif
         !do popo=1, Hsize
         !   if (tempuse.eq.Hsp(popo)) then
               !H(ii,popo)=H(ii,popo)+dcmplx(xi_2p*U_SO_2p(kk-N3d-N, 2, jj-N3d-N, 2)*tempsign)
               ncount = ncount+1
               IndexI(ncount) = ii
               IndexJ(ncount) = popo
               sparseH(ncount) = dcmplx(xi_2p*U_SO_2p(kk-N3d-N, 2, jj-N3d-N, 2)*tempsign)
         !   endif
         !enddo

      enddo!kk
   enddo!jj


!   if (ii.le.(Hsize_m+Hsize_0)) then
      !3. c^\dagger_{m,up} c_{mm', dn}
      do jj=1+N3d, N2p+N3d
         if (BTEST(Hsp(ii), jj-1)) then
            tempuse_1= IBCLR(Hsp(ii), jj-1)
            tempcount=obcount(Hsp(ii), jj-1+1)
            tempsign_1= (-1)**(tempcount)
         else
            cycle
         endif

         do kk=1+N3d+N, N2p+N3d+N
            if(abs(U_SO_2p(kk-N3d-N, 2, jj-N3d, 1)).lt.0.0001) cycle
            if (BTEST(tempuse_1, kk-1)) cycle
            tempuse_2=IBSET(tempuse_1, kk-1)
            tempcount=obcount(tempuse_1,kk-1+1)
            tempsign_2=(-1)**(tempcount);
            tempsign=tempsign_1*tempsign_2
            tempuse=tempuse_2

         call BinarySearch(Hsp,Hsizet,tempuse,popo)
         if(popo.eq.-1) then
             write(*,*) tempuse, 'C' 
             stop
         endif
         !   do popo=1, Hsize
         !      if (tempuse.eq.Hsp(popo)) then
                  !H(ii,popo)=H(ii,popo)+dcmplx(xi_2p*U_SO_2p(kk-N3d-N, 2, jj-N3d, 1)*tempsign)
                  ncount = ncount+1
                  IndexI(ncount) = ii
                  IndexJ(ncount) = popo
                  sparseH(ncount) = dcmplx(xi_2p*U_SO_2p(kk-N3d-N, 2, jj-N3d, 1)*tempsign)
!                  write(*,*) '3           ', ii, popo
!                  write(*,*) '3 kk and jj ', kk-N3d-N, jj-N3d
!                  write(*,*) '3 ', U_SO_2p(kk-N3d-N, 2, jj-N3d, 1), popo
!                  write(*,*) '3 ', tempsign, dcmplx(xi_2p*U_SO_2p(kk-N3d-N, 2, jj-N3d, 1)*tempsign)
         !      endif
         !   enddo

         enddo!kk
      enddo!jj
!   endif


!   if (ii.ge.(Hsize_m+1)) then
      !4. c^\dagger_{m,dn} c_{mm', up}
      do jj=1+N3d+N, N2p+N3d+N  !up
         if (BTEST(Hsp(ii), jj-1)) then
            tempuse_1= IBCLR(Hsp(ii), jj-1)
            tempcount=obcount(Hsp(ii), jj-1+1)
            tempsign_1= (-1)**(tempcount)
         else
            cycle
         endif

         do kk=1+N3d, N2p+N3d   !dn
            if(abs(U_SO_2p(kk-N3d, 1, jj-N3d-N, 2)).lt.0.0001) cycle
            if (BTEST(tempuse_1, kk-1)) cycle
            tempuse_2=IBSET(tempuse_1, kk-1)
            tempcount=obcount(tempuse_1,kk-1+1)
            tempsign_2=(-1)**(tempcount);
            tempsign=tempsign_1*tempsign_2
            tempuse=tempuse_2

         call BinarySearch(Hsp,Hsizet,tempuse,popo)
         if(popo.eq.-1) then
             write(*,'(B30)') tempuse  
             write(*,*) 'D' 
             stop
         endif
         !   do popo=1, Hsize
         !      if (tempuse.eq.Hsp(popo)) then
                 ! H(ii,popo)=H(ii,popo)+dcmplx(xi_2p*U_SO_2p(kk-N3d, 1, jj-N3d-N, 2)*tempsign)
                  ncount = ncount+1
                  IndexI(ncount) = ii
                  IndexJ(ncount) = popo
                  sparseH(ncount) = dcmplx(xi_2p*U_SO_2p(kk-N3d, 1, jj-N3d-N, 2)*tempsign)
 !                 write(*,*) '4           ', ii, popo
 !                 write(*,*) '4 kk and jj ', kk-N3d, jj-N3d-N
 !                 write(*,*) '4 U_SO_2p   ', U_SO_2p(kk-N3d, 1, jj-N3d-N, 2)
 !                 write(*,*) '4 ', tempsign, dcmplx(xi_2p*U_SO_2p(kk-N3d, 1, jj-N3d-N, 2)*tempsign)
         !      endif
         !   enddo

         enddo!kk
      enddo!jj
!   endif

enddo!ii

SprSizet = ncount

end


!***********************************************************************************
subroutine GenHsp2pfull(spin, Hsp, Hsizet)
use NumOfOrbitalAndElectrons; use ConstantParas
implicit none

integer*8, INTENT(IN) :: spin
integer*8, DIMENSION(Hsize_3dEst), INTENT(OUT) :: Hsp
integer*8, INTENT(OUT) :: Hsizet

integer*8 :: istat, ii, flag2pfull, jj
integer*8, external :: sumeverybit
integer*8 :: tspin


flag2pfull = 7*(2**N3d) + 7*(2**(N3d+N)) !test is 2p is full
ii=0
do istat=0,2**(2*N)-1
    if(sumeverybit(istat).eq.ntol .and. &
      IAND(istat, flag2pfull).eq.flag2pfull) then   !total spin is safistied

        tspin=0
        do jj=0,N-1 
           if(BTEST(istat,jj)) tspin=tspin+1
           if(BTEST(istat,jj+N)) tspin=tspin-1
        enddo
!        write(*,*) 'tspin', tspin, spin
        if(tspin.eq.spin) then

           ii=ii+1
           Hsp(ii)=istat;
        endif

    endif
enddo
Hsizet = ii

end subroutine

!************************************************

subroutine GenHsp2p1hole(Hsp, Hsizet)
use NumOfOrbitalAndElectrons; use ConstantParas
implicit none

integer*8, DIMENSION(HsizeEst), INTENT(OUT) :: Hsp
integer*8, INTENT(OUT) :: Hsizet

integer*8 :: istat, ii
integer*8, external :: sumeverybit
integer*8, external :: sumbeforebit

ii=0
do istat=0,2**(2*N)-1
    if((sumeverybit(istat).eq.ntol) .and. &
       ((sumbeforebit(istat, N3d)+sumbeforebit(istat/2**N, N3d)).eq.(ntol-5))) then   
       ! with one hole at 2p 
        ii=ii+1
        Hsp(ii)=istat;
    endif
enddo
Hsizet = ii

end subroutine



!*******************************************************************************
subroutine PGenMatrix(CHsite, Hspt, Hsizet, nloct, lstart, lend, IndexIt, IndexJt, sparseHt, SprSizet) 
use NumOfOrbitalAndElectrons; use ModelParas; use MPIParas; use ConstantParas
use BettsCluster
implicit none

integer*8, INTENT(IN) :: CHsite, Hsizet
integer*8, INTENT(IN) :: nloct, lstart, lend
integer*8, DIMENSION(Hsizet), INTENT(IN) :: Hspt
integer*8, INTENT(INOUT) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(OUT) :: IndexIt, IndexJt
double complex, DIMENSION(SprSizet), INTENT(OUT) :: sparseHt

integer*8 :: stpt
integer*8 :: hcount, count
integer*8 :: jj, ii, kk, i, j, toto, mm
integer*8 :: tempuse_1, tempuse_2
integer*8 :: reduc_num, cb_count
integer*8, external :: factorial
integer*8, allocatable :: hh(:), ulim(:), output(:), ini_output(:)
integer*8 :: sumeverybit, sumbeforebit
integer*8 :: direction, direction2, initstatup, initstatdn, afterstatup, afterstatdn, l
integer*8 :: kkNear, aftersign, initstat, afterstat, obcount
double precision :: Uctemp, time1
double precision, allocatable :: H(:)

integer*8 :: iiup, iidn, countt, iup, idn, jup, jdn 
integer*8 :: listup(65536),listdn(65536)
integer*8 :: listversup(65536),listversdn(65536)


IndexIt(1:SprSizet)=0
IndexJt(1:SprSizet)=0
sparseHt(1:SprSizet)=dcmplx(0.0d0,0.0d0)

reduc_num = 300

call CPU_time(time1)
!write(*,*) 'GenEachSpinBlock starts at time ',time1 
!write(*,*) 'nup and ndn = ', nup, ndn
!write(*,*) 'Hsize = ', Hsizet

allocate(H(reduc_num))
allocate(hh(reduc_num))

cb_count = 0
!*****Generating the tight-binding Hamiltonian*****
do ii=lstart, lend! Go through all the Hilbert space states
!do ii=stpt, stpt+nloc-1
     initstat = Hspt(ii)

        !if(mod(ii,100000).eq.0) then
        !    call CPU_time(time1)
        !    write(*,*) ii, time1
        !endif
        H(:)=0.0d0;
        hh(:)=0;
        hcount=1
        !The site energy
        do jj=0, 2*N-1
                if(BTEST(initstat,jj)) then
                        H(hcount)=H(hcount)+E_site(jj+1)
                        hh(hcount)=ii
                        hcount=hcount+1
                endif
        enddo

        !The on-site U
        if(nup*ndn.ne.0) then
           do jj=0, N-1
              if(BTEST(initstat,jj).AND.BTEST(initstat,jj+N)) then
                 H(hcount)=H(hcount)+U(jj+1)
                 hh(hcount)=ii
                 hcount=hcount+1
              end if
           end do
        endif

        !The core-hole potential Uc
        Uctemp = 0.0d0;
        do jj=N3d, N-1  !loop over all the 2p states to obtain the effective Uc
              if(.not.BTEST(initstat,jj)) then
                 Uctemp = Uctemp + Uc
              endif
              if(.not.BTEST(initstat,jj+N)) then
                 Uctemp = Uctemp + Uc
              endif
        enddo
        !do jj=0,N3d-1 !particular site CHsite
        !do jj=CHsite !particular site CHsite
        jj=CHsite
              if(BTEST(initstat,jj)) then
                 H(hcount)=H(hcount)+Uctemp
                 hh(hcount)=ii
                 hcount=hcount+1
              endif
              if(BTEST(initstat,jj+N)) then
                 H(hcount)=H(hcount)+Uctemp
                 hh(hcount)=ii
                 hcount=hcount+1
              endif
        !enddo
        !3d Betts cluster hopping
        !initstatup = Hspt(ii)/(2**N)
        !initstatdn = mod(Hspt(ii), 2**N)
        initstat = Hspt(ii)
       if(Ndim.eq.2) then
        do kk = 0, N3d-1
           do direction = 1, 4  !nearest neighbor hopping
              kkNear=site(kk,direction)
              if(BTEST(initstat,kk).eqv..true.) then
                 if(BTEST(initstat,kkNear).eqv..false.) then
                    afterstat=IBCLR(initstat,kk)
                    aftersign=(-1)**(obcount(initstat,kk+1))
                    afterstat=IBSET(afterstat,kkNear)
                    aftersign=aftersign*(-1)**(obcount(afterstat,kkNear+1))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*t
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
              if(BTEST(initstat,kk+N).eqv..true.) then
                 if(BTEST(initstat,kkNear+N).eqv..false.) then
                    afterstat=IBCLR(initstat,kk+N)
                    aftersign=(-1)**(obcount(initstat,kk+1+N))
                    afterstat=IBSET(afterstat,kkNear+N)
                    aftersign=aftersign*((-1)**obcount(afterstat,kkNear+1+N))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*t
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
           enddo

           do direction = 1, 2  !next nearest neighbor hopping
            do direction2 = 3, 4  !next nearest neighbor hopping
              kkNear=site(site(kk,direction2),direction)
              if(BTEST(initstat,kk).eqv..true.) then
                 if(BTEST(initstat,kkNear).eqv..false.) then
                    afterstat=IBCLR(initstat,kk)
                    aftersign=(-1)**(obcount(initstat,kk+1))
                    afterstat=IBSET(afterstat,kkNear)
                    aftersign=aftersign*(-1)**(obcount(afterstat,kkNear+1))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*tt
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
              !if(BTEST(initstatdn,kk).eqv..true.) then
              if(BTEST(initstat,kk+N).eqv..true.) then
                 if(BTEST(initstat,kkNear+N).eqv..false.) then
                    afterstat=IBCLR(initstat,kk+N)
                    aftersign=(-1)**(obcount(initstat,kk+1+N))
                    afterstat=IBSET(afterstat,kkNear+N)
                    aftersign=aftersign*((-1)**obcount(afterstat,kkNear+1+N))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*tt
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
            enddo
           enddo

        enddo


       else if(Ndim.eq.1) then
        do kk = 0, N3d-1
           do direction = 1, 2  !nearest neighbor hopping
              kkNear=site(kk,direction)
              if(BTEST(initstat,kk).eqv..true.) then
                 if(BTEST(initstat,kkNear).eqv..false.) then
                    afterstat=IBCLR(initstat,kk)
                    aftersign=(-1)**(obcount(initstat,kk+1))
                    afterstat=IBSET(afterstat,kkNear)
                    aftersign=aftersign*(-1)**(obcount(afterstat,kkNear+1))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*t
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
              if(BTEST(initstat,kk+N).eqv..true.) then
                 if(BTEST(initstat,kkNear+N).eqv..false.) then
                    afterstat=IBCLR(initstat,kk+N)
                    aftersign=(-1)**(obcount(initstat,kk+1+N))
                    afterstat=IBSET(afterstat,kkNear+N)
                    aftersign=aftersign*((-1)**obcount(afterstat,kkNear+1+N))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*t
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
           enddo
        enddo

       endif

        !******Generating the Sparse Matrix*******
        if((hcount-1).gt.reduc_num) then
           write(*,*) 'reduc_num is too small for row:', ii
           !call MPI_FINALIZE()
           stop
        endif

        do jj=1, hcount-1
           !To combine contritubitions from the same indices:
           do kk=jj+1, hcount-1
              if(hh(kk).eq.hh(jj)) then
                 H(jj)=H(jj)+H(kk)
                 H(kk)=0.0d0;
              endif
           enddo

           if(abs(H(jj)).gt.1.0d-6) then
              cb_count=cb_count+1
              sparseHt(cb_count)= H(jj)
              IndexIt(cb_count)=  ii
              IndexJt(cb_count)=  hh(jj)
           endif
        enddo

end do!end do for ii=1, Hsize


SprSizet = cb_count

if(cb_count.ge.size(IndexIt)) then
   write(*,*) 'the temp_IndexI array size is small; change its allocation size.'
   stop
endif

deallocate(H,hh)

end subroutine


!*******************************************************************************
subroutine GenMatrix(CHsite, Hspt, Hsizet, IndexIt, IndexJt, sparseHt, SprSizet) 
use NumOfOrbitalAndElectrons; use ModelParas; use MPIParas; use ConstantParas
use BettsCluster
implicit none

integer*8, INTENT(IN) :: CHsite, Hsizet
integer*8, DIMENSION(Hsizet), INTENT(IN) :: Hspt
integer*8, INTENT(INOUT) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(OUT) :: IndexIt, IndexJt
double complex, DIMENSION(SprSizet), INTENT(OUT) :: sparseHt

integer*8 :: stpt
integer*8 :: hcount, count
integer*8 :: jj, ii, kk, i, j, toto, mm
integer*8 :: tempuse_1, tempuse_2
integer*8 :: reduc_num, cb_count
integer*8, external :: factorial
integer*8, allocatable :: hh(:), ulim(:), output(:), ini_output(:)
integer*8 :: sumeverybit, sumbeforebit
integer*8 :: direction, direction2, initstatup, initstatdn, afterstatup, afterstatdn, l
integer*8 :: kkNear, aftersign, initstat, afterstat, obcount
double precision :: Uctemp, time1
double precision, allocatable :: H(:)

integer*8 :: iiup, iidn, countt, iup, idn, jup, jdn 
integer*8 :: listup(65536),listdn(65536)
integer*8 :: listversup(65536),listversdn(65536)


IndexIt(1:SprSizet)=0
IndexJt(1:SprSizet)=0
sparseHt(1:SprSizet)=dcmplx(0.0d0,0.0d0)

reduc_num = 300

call CPU_time(time1)
!write(*,*) 'GenEachSpinBlock starts at time ',time1 
!write(*,*) 'nup and ndn = ', nup, ndn
!write(*,*) 'Hsize = ', Hsizet

allocate(H(reduc_num))
allocate(hh(reduc_num))

cb_count = 0
!*****Generating the tight-binding Hamiltonian*****
do ii=1, Hsizet! Go through all the Hilbert space states
!do ii=stpt, stpt+nloc-1
     initstat = Hspt(ii)

        !if(mod(ii,100000).eq.0) then
        !    call CPU_time(time1)
        !    write(*,*) ii, time1
        !endif
        H(:)=0.0d0;
        hh(:)=0;
        hcount=1
        !The site energy
        do jj=0, 2*N-1
                if(BTEST(initstat,jj)) then
                        H(hcount)=H(hcount)+E_site(jj+1)
                        hh(hcount)=ii
                        hcount=hcount+1
                endif
        enddo

        !The on-site U
        if(nup*ndn.ne.0) then
           do jj=0, N-1
              if(BTEST(initstat,jj).AND.BTEST(initstat,jj+N)) then
                 H(hcount)=H(hcount)+U(jj+1)
                 hh(hcount)=ii
                 hcount=hcount+1
              end if
           end do
        endif

        !The core-hole potential Uc
        Uctemp = 0.0d0;
        do jj=N3d, N-1  !loop over all the 2p states to obtain the effective Uc
              if(.not.BTEST(initstat,jj)) then
                 Uctemp = Uctemp + Uc
              endif
              if(.not.BTEST(initstat,jj+N)) then
                 Uctemp = Uctemp + Uc
              endif
        enddo
        !do jj=0,N3d-1 !particular site CHsite
        !do jj=CHsite !particular site CHsite
        jj=CHsite
              if(BTEST(initstat,jj)) then
                 H(hcount)=H(hcount)+Uctemp
                 hh(hcount)=ii
                 hcount=hcount+1
              endif
              if(BTEST(initstat,jj+N)) then
                 H(hcount)=H(hcount)+Uctemp
                 hh(hcount)=ii
                 hcount=hcount+1
              endif
        !enddo
        !3d Betts cluster hopping
        !initstatup = Hspt(ii)/(2**N)
        !initstatdn = mod(Hspt(ii), 2**N)
        initstat = Hspt(ii)
       if(Ndim.eq.2) then
        do kk = 0, N3d-1
           do direction = 1, 4  !nearest neighbor hopping
              kkNear=site(kk,direction)
              if(BTEST(initstat,kk).eqv..true.) then
                 if(BTEST(initstat,kkNear).eqv..false.) then
                    afterstat=IBCLR(initstat,kk)
                    aftersign=(-1)**(obcount(initstat,kk+1))
                    afterstat=IBSET(afterstat,kkNear)
                    aftersign=aftersign*(-1)**(obcount(afterstat,kkNear+1))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*t*prefac(direction)
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
              if(BTEST(initstat,kk+N).eqv..true.) then
                 if(BTEST(initstat,kkNear+N).eqv..false.) then
                    afterstat=IBCLR(initstat,kk+N)
                    aftersign=(-1)**(obcount(initstat,kk+1+N))
                    afterstat=IBSET(afterstat,kkNear+N)
                    aftersign=aftersign*((-1)**obcount(afterstat,kkNear+1+N))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*t*prefac(direction)
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
           enddo

           do direction = 1, 2  !next nearest neighbor hopping
            do direction2 = 3, 4  !next nearest neighbor hopping
              kkNear=site(site(kk,direction2),direction)
              if(BTEST(initstat,kk).eqv..true.) then
                 if(BTEST(initstat,kkNear).eqv..false.) then
                    afterstat=IBCLR(initstat,kk)
                    aftersign=(-1)**(obcount(initstat,kk+1))
                    afterstat=IBSET(afterstat,kkNear)
                    aftersign=aftersign*(-1)**(obcount(afterstat,kkNear+1))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*tt*prefac(direction)*prefac(direction2)
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
              !if(BTEST(initstatdn,kk).eqv..true.) then
              if(BTEST(initstat,kk+N).eqv..true.) then
                 if(BTEST(initstat,kkNear+N).eqv..false.) then
                    afterstat=IBCLR(initstat,kk+N)
                    aftersign=(-1)**(obcount(initstat,kk+1+N))
                    afterstat=IBSET(afterstat,kkNear+N)
                    aftersign=aftersign*((-1)**obcount(afterstat,kkNear+1+N))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*tt*prefac(direction)*prefac(direction2)
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
            enddo
           enddo

        enddo


       else if(Ndim.eq.1) then
        do kk = 0, N3d-1
           do direction = 1, 2  !nearest neighbor hopping
              kkNear=site(kk,direction)
              if(BTEST(initstat,kk).eqv..true.) then
                 if(BTEST(initstat,kkNear).eqv..false.) then
                    afterstat=IBCLR(initstat,kk)
                    aftersign=(-1)**(obcount(initstat,kk+1))
                    afterstat=IBSET(afterstat,kkNear)
                    aftersign=aftersign*(-1)**(obcount(afterstat,kkNear+1))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*t
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
              if(BTEST(initstat,kk+N).eqv..true.) then
                 if(BTEST(initstat,kkNear+N).eqv..false.) then
                    afterstat=IBCLR(initstat,kk+N)
                    aftersign=(-1)**(obcount(initstat,kk+1+N))
                    afterstat=IBSET(afterstat,kkNear+N)
                    aftersign=aftersign*((-1)**obcount(afterstat,kkNear+1+N))
                    call BinarySearch(Hspt,Hsizet,afterstat,l)
                    if(l.eq.-1) stop
                    H(hcount)=H(hcount)-aftersign*t
                    hh(hcount)=l
                    hcount=hcount+1
                 endif
              endif
           enddo
        enddo

       endif

        !******Generating the Sparse Matrix*******
        if((hcount-1).gt.reduc_num) then
           write(*,*) 'reduc_num is too small for row:', ii
           !call MPI_FINALIZE()
           stop
        endif

        do jj=1, hcount-1
           !To combine contritubitions from the same indices:
           do kk=jj+1, hcount-1
              if(hh(kk).eq.hh(jj)) then
                 H(jj)=H(jj)+H(kk)
                 H(kk)=0.0d0;
              endif
           enddo

           if(abs(H(jj)).gt.1.0d-6) then
              cb_count=cb_count+1
              sparseHt(cb_count)= H(jj)
              IndexIt(cb_count)=  ii
              IndexJt(cb_count)=  hh(jj)
           endif
        enddo

end do!end do for ii=1, Hsize

SprSizet = cb_count

if(cb_count.ge.size(IndexIt)) then
   write(*,*) 'the temp_IndexI array size is small; change its allocation size.'
   stop
endif

deallocate(H,hh)

end subroutine

!********************************************************************************

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
subroutine BinarySearch(listqpt,endd,statupdn,l)
implicit none
INTEGER*8, INTENT(IN) :: endd
INTEGER*8, DIMENSION(max(endd,1)), INTENT(IN) :: listqpt
INTEGER*8, INTENT(IN) :: statupdn
INTEGER*8, INTENT(OUT):: l

integer*8:: head, tail, middle

!write(*,*) "in Binary Loop", endd
if(endd.lt.1) then
   l=-1;
   goto 200
endif
!write(*,*) "Pass The first test", endd

head=1; tail=endd; middle=(head+tail)/2

do while((listqpt(middle).ne.statupdn).and.(head.le.tail))
   if(statupdn.gt.listqpt(middle)) then
      head = middle+1
   else
      tail = middle-1
   endif
   middle = (head+tail)/2
   if((middle.gt.endd).or.(middle.lt.1)) then
     l=-1
     goto 200
   endif
enddo

if(listqpt(middle).eq.statupdn) then
   l=middle;
else
   l=-1;
endif

200 tail = endd

end

!*********************************************************
!*******************************************************

integer*8 function sumeverybit(tempa)
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

!*********************************************************
!********************************************************
subroutine GenODaggerStatesnew(CHsite, H_0, H_e, Hsp_src, Hsp_dst, Hsize_src, Hsize_dst, Gn)
use NumOfOrbitalAndElectrons
implicit none

integer*8, INTENT(IN) :: CHsite, Hsize_src, Hsize_dst
integer*8, DIMENSION(Hsize_src), INTENT(IN) :: Hsp_src
integer*8, DIMENSION(Hsize_dst), INTENT(IN) :: Hsp_dst
double complex, DIMENSION(Hsize_src), INTENT(IN) :: H_0
double complex, DIMENSION(Hsize_dst), INTENT(OUT) :: H_e
double precision, DIMENSION(3), INTENT(IN) :: Gn

integer*8 :: ii, jj, interstate, afterstate, l
integer*8 :: aftersign, sumbeforebit, iikept
integer*8 :: obcount, CHsitei


H_e=dcmplx(0.0d0,0.0d0)
CHsitei = CHsite
   do ii=1,Hsize_src
      do jj=N3d, N-1
         if(.not.BTEST(Hsp_src(ii),jj).and.(BTEST(Hsp_src(ii),CHsitei))) then
            interstate = IBCLR(Hsp_src(ii), CHsitei)
            aftersign = (-1)**obcount(Hsp_src(ii), CHsitei+1)
            afterstate = IBSET(interstate, jj)         
            aftersign = aftersign * (-1)**obcount(interstate, jj+1)
            call BinarySearch(Hsp_dst, Hsize_dst, afterstate, l)

            if(l.ne.-1) H_e(l)=H_e(l)+H_0(ii)*aftersign*Gn(jj-N3d+1)

         endif
         if(.not.BTEST(Hsp_src(ii),jj+N).and.(BTEST(Hsp_src(ii),CHsitei+N))) then
            interstate = IBCLR(Hsp_src(ii), CHsitei+N)
            aftersign = (-1)**obcount(Hsp_src(ii), CHsitei+N+1)
            afterstate = IBSET(interstate, jj+N)         
            aftersign = aftersign * (-1)**obcount(interstate, jj+N+1)
            call BinarySearch(Hsp_dst, Hsize_dst, afterstate, l)
          
            if(l.ne.-1) H_e(l)=H_e(l)+H_0(ii)*aftersign*Gn(jj-N3d+1)

         endif
     enddo
   enddo

end subroutine




subroutine GenODaggerMatrix(CHsite, Hsp_src, Hsp_dst, Hsize_src, Hsize_dst, Gn, &
                            IndexIt, IndexJt, sparseHt, SprSizet) 

use NumOfOrbitalAndElectrons; use ModelParas; use PolarParas
implicit none

integer*8, INTENT(IN) :: CHsite, Hsize_src, Hsize_dst
integer*8, DIMENSION(Hsize_src), INTENT(IN) :: Hsp_src
integer*8, DIMENSION(Hsize_dst), INTENT(in):: Hsp_dst
double precision, DIMENSION(3), INTENT(IN) :: Gn
integer*8, INTENT(INOUT) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(OUT) :: IndexIt, IndexJt
double complex, DIMENSION(SprSizet), INTENT(OUT) :: sparseHt

integer*8 :: ii, jj, interstate, afterstate, l
integer*8 :: aftersign, sumbeforebit, iikept
integer*8 :: obcount, CHsitei,iitemp


IndexIt(1:SprSizet)=0
IndexJt(1:SprSizet)=0
sparseHt(1:SprSizet)=dcmplx(0.0d0,0.0d0)

iitemp=0
CHsitei = CHsite
   do ii=1,Hsize_src
      do jj=N3d, N-1


         if(.not.BTEST(Hsp_src(ii),jj).and.(BTEST(Hsp_src(ii),CHsitei))) then
            interstate = Hsp_src(ii)
            interstate = IBCLR(interstate, CHsitei)
            aftersign = (-1)**obcount(interstate, CHsitei+1)
            afterstate = IBSET(interstate, jj)         
            aftersign = aftersign * (-1)**obcount(interstate, jj+1)
            call BinarySearch(Hsp_dst, Hsize_dst, afterstate, l)
            if(l.eq.-1) then
               write(*,*) 'Gen wrong A'
               write(*,'(B30)') afterstate
               stop
            endif
            !do iikept=1,Nkept
            !   H_e(l)=H_e(l)+H_0(ii)*aftersign*Gn(jj-N3d+1)
            !enddo
            !do l=1,Hsize
            !   if(afterstate.eq.Hsp(l)) then 
            !      do iikept=1,Nkept
            !         H_e(l,iikept)=H_e(l,iikept)+H_0(ii,iikept)*aftersign
            !      enddo
            !   endif
            !enddo
           if(abs(Gn(jj-N3d+1)).gt.0.0001d0) then
            iitemp=iitemp+1
            IndexJt(iitemp)=ii
            IndexIt(iitemp)=l             
            sparseHt(iitemp)= aftersign * Gn(jj-N3d+1) !+ (jj-N3d+1)*10*dcmplx(0.0d0,1.0d0)
           endif
         endif

         if(.not.BTEST(Hsp_src(ii),jj+N).and.(BTEST(Hsp_src(ii),CHsitei+N))) then
            interstate = Hsp_src(ii)
            interstate = IBCLR(interstate, CHsitei+N)
            aftersign = (-1)**obcount(interstate, CHsitei+N+1)
            afterstate = IBSET(interstate, jj+N)         
            aftersign = aftersign * (-1)**obcount(interstate, jj+N+1)
            call BinarySearch(Hsp_dst, Hsize_dst, afterstate, l)
            if(l.eq.-1) then
               write(*,*) 'Gen wrong B '
               write(*,'(B30)') afterstate
               stop
            endif
            !do iikept=1,Nkept
            !   H_e(l)=H_e(l)+H_0(ii)*aftersign*Gn(jj-N3d+1)
            !enddo
            !do l=1,Hsize
            !   if(afterstate.eq.Hsp(l)) then 
            !      do iikept=1,Nkept
            !         H_e(l,iikept)=H_e(l,iikept)+H_0(ii,iikept)*aftersign
            !      enddo
            !   endif
            !enddo
           if(abs(Gn(jj-N3d+1)).gt.0.0001d0) then
            iitemp=iitemp+1
            IndexJt(iitemp)=ii
            IndexIt(iitemp)=l             
            sparseHt(iitemp)= aftersign * Gn(jj-N3d+1) !- (jj-N3d+1)*10*dcmplx(0.0d0,1.0d0)
           endif 

         endif

      !  endif !Gn

     enddo
   enddo

SprSizet=iitemp

!enddo
end subroutine



subroutine GenOMatrix(CHsite, Hsp_src, Hsp_dst, Hsize_src, Hsize_dst, Gn, &
                      IndexIt, IndexJt, sparseHt, SprSizet) 

use NumOfOrbitalAndElectrons; use ModelParas; use PolarParas
implicit none

integer*8, INTENT(IN) :: CHsite, Hsize_src, Hsize_dst
integer*8, DIMENSION(Hsize_src), INTENT(IN) :: Hsp_src
integer*8, DIMENSION(Hsize_dst), INTENT(in):: Hsp_dst
double precision, DIMENSION(3), INTENT(IN) :: Gn
integer*8, INTENT(INOUT) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(OUT) :: IndexIt, IndexJt
double complex, DIMENSION(SprSizet), INTENT(OUT) :: sparseHt

integer*8 :: ii, jj, interstate, afterstate, l
integer*8 :: aftersign, sumbeforebit, iikept
integer*8 :: obcount, CHsitei,iitemp


IndexIt(1:SprSizet)=0
IndexJt(1:SprSizet)=0
sparseHt(1:SprSizet)=dcmplx(0.0d0,0.0d0)

iitemp=0
CHsitei = CHsite
!do iikept = 1, Nkept
   do ii=1,Hsize_src
      do jj=N3d, N-1

 !       if(abs(Gn(jj-N3d+1)).gt.0.0001d0) then

        !do CHsitei=0,N3d-1
         if(BTEST(Hsp_src(ii),jj).and.(.not.BTEST(Hsp_src(ii),CHsitei))) then
            interstate = Hsp_src(ii)
            interstate = IBCLR(interstate, jj)
            aftersign = (-1)**obcount(interstate, jj+1)
            !aftersign = (-1)**sumbeforebit(Hsp(ii), jj)
            afterstate = IBSET(interstate, CHsitei)
            aftersign = aftersign * (-1)**obcount(interstate, CHsitei+1)
            !aftersign = aftersign * (-1)**sumbeforebit(interstate, CHsite)
            call BinarySearch(Hsp_dst, Hsize_dst, afterstate, l)
            if(l.eq.-1) stop
            !do iikept=1,Nkept
            !   H_e(l)=H_e(l)+H_0(ii)*aftersign * Gn(jj-N3d+1)
            !enddo
           if(abs(Gn(jj-N3d+1)).gt.0.0001d0) then
            iitemp=iitemp+1
            IndexJt(iitemp)=ii
            IndexIt(iitemp)=l             
            sparseHt(iitemp)= aftersign * Gn(jj-N3d+1) !+ (jj-N3d+1)*dcmplx(0.0d0,1.0d0)
           endif 

         endif

         if(BTEST(Hsp_src(ii),jj+N).and.(.not.BTEST(Hsp_src(ii),CHsitei+N))) then
            interstate = Hsp_src(ii)
            interstate = IBCLR(interstate, jj+N)
            aftersign = (-1)**obcount(interstate, jj+N+1)
            !aftersign = (-1)**sumbeforebit(Hsp(ii), jj+N)
            afterstate = IBSET(interstate, CHsitei+N)
            aftersign = aftersign * (-1)**obcount(interstate, CHsitei+N+1)
            !aftersign = aftersign * (-1)**sumbeforebit(interstate, CHsite+N)
            call BinarySearch(Hsp_dst, Hsize_dst, afterstate, l)
            if(l.eq.-1) stop
            !do iikept=1,Nkept
            !   H_e(l)=H_e(l)+H_0(ii)*aftersign * Gn(jj-N3d+1)
            !enddo
           if(abs(Gn(jj-N3d+1)).gt.0.0001d0) then
            iitemp=iitemp+1
            IndexJt(iitemp)=ii
            IndexIt(iitemp)=l             
            sparseHt(iitemp)= aftersign * Gn(jj-N3d+1) !- (jj-N3d+1)*dcmplx(0.0d0,1.0d0)
           endif

         endif


        !endif
       !enddo !CHsitei
     enddo
   enddo

SprSizet=iitemp

!enddo
end subroutine






!**********************************************
subroutine GenOStatesnew(CHsite, H_0, H_e, Hsp_src, Hsp_dst, Hsize_src, Hsize_dst, Gn)
use NumOfOrbitalAndElectrons; use PolarParas
implicit none

integer*8, INTENT(IN) :: CHsite, Hsize_src, Hsize_dst
integer*8, DIMENSION(Hsize_src), INTENT(IN) :: Hsp_src
integer*8, DIMENSION(Hsize_dst), INTENT(in):: Hsp_dst
double complex, DIMENSION(Hsize_src), INTENT(IN) :: H_0
double complex, DIMENSION(Hsize_dst), INTENT(OUT) :: H_e
double precision, DIMENSION(3), INTENT(IN) :: Gn

integer*8 :: ii, jj, interstate, afterstate, l
integer*8 :: aftersign, sumbeforebit, iikept
integer*8 :: obcount, CHsitei

H_e=dcmplx(0.0d0,0.0d0)
CHsitei = CHsite
!do iikept = 1, Nkept
   do ii=1,Hsize_src
      do jj=N3d, N-1
        !do CHsitei=0,N3d-1
         if(BTEST(Hsp_src(ii),jj).and.(.not.BTEST(Hsp_src(ii),CHsitei))) then
            interstate = IBCLR(Hsp_src(ii), jj)
            aftersign = (-1)**obcount(Hsp_src(ii), jj+1)
            !aftersign = (-1)**sumbeforebit(Hsp(ii), jj)
            afterstate = IBSET(interstate, CHsitei)
            aftersign = aftersign * (-1)**obcount(interstate, CHsitei+1)
            !aftersign = aftersign * (-1)**sumbeforebit(interstate, CHsite)
            call BinarySearch(Hsp_dst, Hsize_dst, afterstate, l)
            if(l.eq.-1) stop
            !do iikept=1,Nkept
               H_e(l)=H_e(l)+H_0(ii)*aftersign * Gn(jj-N3d+1)
            !enddo
         endif
         if(BTEST(Hsp_src(ii),jj+N).and.(.not.BTEST(Hsp_src(ii),CHsitei+N))) then
            interstate = IBCLR(Hsp_src(ii), jj+N)
            aftersign = (-1)**obcount(Hsp_src(ii), jj+N+1)
            !aftersign = (-1)**sumbeforebit(Hsp(ii), jj+N)
            afterstate = IBSET(interstate, CHsitei+N)
            aftersign = aftersign * (-1)**obcount(interstate, CHsitei+N+1)
            !aftersign = aftersign * (-1)**sumbeforebit(interstate, CHsite+N)
            call BinarySearch(Hsp_dst, Hsize_dst, afterstate, l)
            if(l.eq.-1) stop
            !do iikept=1,Nkept
               H_e(l)=H_e(l)+H_0(ii)*aftersign * Gn(jj-N3d+1)
            !enddo
         endif
       !enddo !CHsitei
     enddo
   enddo
!enddo
end subroutine



!************************************
double complex function VectorProductComplex(arrayA, arrayB, length)
implicit none
integer*8, INTENT(IN) :: length
double complex, DIMENSION(length), INTENT(IN) :: arrayA, arrayB

integer*8 :: ii

VectorProductComplex = dcmplx(0.0d0,0.0d0);
do ii=1,length
   VectorProductComplex=VectorProductComplex+conjg(arrayA(ii))*arrayB(ii)
enddo
return
end function
!*****Lorentzian*****
double precision function Lorentzian(temp1, temp2)
implicit none
double precision, parameter::pi=3.1415926535d0
double precision :: temp1, temp2
Lorentzian=(1.0d0/pi)*(temp2)/(temp1**2+temp2**2)
return
end function
!*******************

!*****Count how many occupied states in front of the ith orbital*****
integer*8 function obcount(tempob, tempcount)!tempcount is the labeling of the real orbital
use NumOfOrbitalAndElectrons
implicit none
integer*8:: tempcount, tempob
integer*8:: ii, i
i=0
do ii=tempcount, 2*N-1
if (BTEST(tempob,ii)) i=i+1
end do
obcount=i
return
end function
!*******************
subroutine LinearSearch(mylistqpt,myksize,mystatupdn,myl)
implicit none
INTEGER*8, INTENT(IN) :: myksize
INTEGER*8, DIMENSION(1:myksize), INTENT(IN) :: mylistqpt
INTEGER*8, INTENT(IN) :: mystatupdn
INTEGER*8, INTENT(OUT):: myl

integer*8:: temp
do temp=1, myksize
   if(mylistqpt(temp).eq.mystatupdn) then
      myl = temp
      return
   endif
enddo
myl = -1
return
end

!************************************************************************
subroutine PComplexLanczos(Hsize, nloct, lstart, lend, &
                           spr_size, IndexI, IndexJ, sparseH, E_0, H_0)
use ConstantParas; use MPIParas
implicit none
include 'mpif.h'

integer*8, INTENT(IN) :: Hsize, spr_size
integer*8, INTENT(IN) :: nloct, lstart, lend
integer*8, DIMENSION(spr_size), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(spr_size), INTENT(IN) :: sparseH
double precision, INTENT(OUT) :: E_0
double complex, DIMENSION(Hsize), INTENT(OUT) :: H_0

integer*8 :: ii,jj,kk, i, j, k
integer*8 :: Itemp1, Itemp2, Itemp3
integer*8 :: iter, iter_num
double precision :: rtemp1, rtemp2, rtemp3
double precision :: myrtemp1, myrtemp2, myrtemp3
double precision :: tol, ortho_tol
double precision :: time1, time2, time3
double precision, allocatable :: rarray_1(:), rarray_2(:), rarray_3(:)
double precision, allocatable :: E0(:)
double complex :: ctemp1, ctemp2, ctemp3
double complex :: myctemp1, myctemp2, myctemp3
double complex :: a_n, b_n
double complex, allocatable :: DM(:,:), DM_temp(:,:)
double complex, allocatable :: phi_n1(:), phi_nm(:), phi_n0(:)
double complex, allocatable :: carray_1(:), carray_2(:), carray_3(:)
double complex, allocatable :: mycarray_1(:)
!!For LAPACK:
CHARACTER*10 :: JOBZ, UPLO
integer :: LDA, INFO, LWORK
DOUBLE PRECISION, allocatable ::  RWORK(:), W(:)
double complex, allocatable :: A(:,:), WORK(: )
double complex, allocatable :: ini_vec(:)


call CPU_time(time1)

!!=============================================
!!Part I: Constructing the tri-diagonal matrix:
!!=============================================

allocate(carray_2(nloct)); carray_2=0.0d0;
allocate(carray_1(Hsize)); carray_1=0.0d0;
allocate(mycarray_1(Hsize)); mycarray_1=0.0d0;
allocate(phi_n0(nloct)); phi_n0=0.0d0;
allocate(phi_n1(nloct)); phi_n1=0.0d0;
allocate(phi_nm(nloct)); phi_nm=0.0d0;
allocate(ini_vec(nloct))

myrtemp1=0.0d0;
do ii=1, nloct
   call random_number(rtemp2)
   call random_number(rtemp3)
   phi_n0(ii)=dcmplx(rtemp2,rtemp3)
   myrtemp1=myrtemp1+abs(phi_n0(ii))**2
enddo
call MPI_Allreduce(myrtemp1, rtemp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, rc)

phi_n0=phi_n0/sqrt(rtemp1)

ini_vec(1:nloct)=phi_n0(1:nloct)
!!************************************************
iter_num=niter_CL;
tol=tol_CL
ortho_tol=tol_CL;
!!************************************************

!write(*,*) '*************************************************************'
!write(*,*) 'Maximum iteration nums.:', iter_num
!write(*,*) 'convergence   criterion:', tol
!write(*,*) 'Orthogonality criterion:', ortho_tol
!write(*,*) '*************************************************************'

allocate(DM(iter_num,iter_num)); DM=0.0d0;
allocate(E0(iter_num)); E0=0.0d0;
if(myid.eq.0) open(unit=135,file='E_0.txt',status='unknown')

iter=0; j=0;
do while(j.eq.0)
      
   call MPI_Barrier(comm, rc)

   call CPU_time(time2)
   iter=iter+1
   !write(*,*) 'ITER:', iter

   !!Constructing H|phi_n0>: later want to do in one step so only 3 vectors are needed.
   !!rarray_1=0.0d0;
   mycarray_1=0.0d0;
   do ii=1, spr_size
      mycarray_1(IndexJ(ii))= mycarray_1(IndexJ(ii))+sparseH(ii)*phi_n0(IndexI(ii)-lstart+1);
   enddo
   call MPI_Allreduce(mycarray_1, carray_1, Hsize, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)

   !!Calculate <phi_n0|H|phi_n0>:
   !!rtemp1=0.0d0;
   myctemp1=0.0d0;
   do ii=1, nloct
      !rtemp1=rtemp1+rarray_1(ii)*phi_n0(ii)
      myctemp1=myctemp1+conjg(carray_1(ii+lstart-1))*phi_n0(ii)
   enddo
   call MPI_Allreduce(myctemp1, ctemp1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)

   !!a_n=rtemp1;
   a_n=ctemp1;
   if(iter.eq.1) b_n=0.0d0;
   !!phi_n1=rarray_1-a_n*phi_n0-b_n*phi_nm;

   phi_n1(1:nloct)=carray_1(lstart:lend)-a_n*phi_n0(1:nloct)-b_n*phi_nm(1:nloct);

   myrtemp2=0.0d0;
   do ii=1, nloct
      !!rtemp2=rtemp2+phi_n1(ii)**2
      myrtemp2=myrtemp2+abs(phi_n1(ii))**2
   enddo
   call MPI_Allreduce(myrtemp2, rtemp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, rc)
   phi_n1=phi_n1/sqrt(rtemp2)
   b_n=sqrt(rtemp2)

   !!Defining the tri-diagonal matrix:
   DM(iter,iter)= a_n
   if (iter.lt.iter_num) then
      DM(iter,iter+1)= b_n
      DM(iter+1,iter)= b_n
   endif

   !write(*,*) 'ME:', a_n, b_n

   !!Calculate overlap:
   rtemp1=0.0d0;
   rtemp2=0.0d0;
   rtemp3=0.0d0;
   myctemp1=0.0d0
   myctemp2=0.0d0
   myctemp3=0.0d0

   do ii=1, nloct
      myctemp1=myctemp1+conjg(phi_nm(ii))*phi_n0(ii)
      myctemp2=myctemp2+conjg(phi_nm(ii))*phi_n1(ii)
      myctemp3=myctemp3+conjg(phi_n0(ii))*phi_n1(ii)
   enddo
   call MPI_Allreduce(myctemp1, ctemp1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   call MPI_Allreduce(myctemp2, ctemp2, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   call MPI_Allreduce(myctemp3, ctemp3, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   !write(*,*) 'ortho_test:', ctemp1, ctemp2, ctemp3

   phi_nm=phi_n0
   phi_n0=phi_n1;
   phi_n1=0.0d0;


  !if(myid.eq.0) then

   if (allocated(DM_temp)) deallocate(DM_temp)
   allocate(DM_temp(iter,iter)); DM_temp=0.0d0;
   DM_temp(1,1)=DM(1,1)
   do ii=2, iter
      DM_temp(ii,ii)=DM(ii,ii)
      DM_temp(ii,ii-1)=DM(ii,ii-1)
      DM_temp(ii-1,ii)=DM(ii-1,ii)
   enddo

   !!For lapack:
   JOBZ='V'
   UPLO='U'

   Itemp1=iter
   LDA=Itemp1;
   LWORK=3*Itemp1;
   if (allocated(work)) deallocate(work)
   if (allocated(W)) deallocate(W)
   allocate(work(lwork)); work=0.0d0;
   allocate(W(Itemp1)); W=0.0d0;

   if(allocated(RWORK)) deallocate(RWORK)
   allocate(RWORK(lwork)); RWORK=0.0d0;

   !!SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
   !!call DSYEV(JOBZ, UPLO, Itemp1, DM_temp, LDA, W, Work, LWORK, INFO)
   !!SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
   call zheev(JOBZ, UPLO, Itemp1, DM_temp, LDA, W, WORK, LWORK, RWORK, info)
   if (info.eq.0) then
      continue
   else
      !write(*,*) 'DSYEV failed.'
      write(*,*) 'ZHEEV failed.'
      write(*,*) 'INFO:', INFO
      stop
   endif

   E0(iter)=W(1)
   !write(*,*) 'eigenvalue:', W(1)
   E_0=W(1)   !I added
   !if (iter.ge.5) then
   !   write(135,*) iter, W(1), W(2), W(3), W(4), W(5)
   !endif

   if(iter.ge.2) then
      if(abs(E0(iter)-E0(iter-1)).le.tol) j=1
   endif

   if(iter.eq.iter_num) j=3

 !  call CPU_time(time3)
 !  write(*,*) 'Lanczos iter ',iter,' takes ',time3-time2,' secs'


  !endif !nmyid.eq.0

  call MPI_Barrier(comm, rc)

enddo

if(myid.eq.0) then
   write(135,*) iter, W(1), W(2), W(3), W(4), W(5)
   close(135)

   if (j.eq.1) write(*,*) 'convergence reached at iter:', iter
   if (j.eq.2) write(*,*) 'orthogonality leaks at iter:', iter
   if (j.eq.3) write(*,*) 'maximum iteration reached.', iter
   write(*,*) ''

endif


!do ii=1, iter
!   write(*,*) 'DM', DM_temp(ii,1), DM_temp(ii,2)
!enddo
!write(*,*) ''

!!====================================
!!Part II: Extracting the eigenvector:
!!====================================
phi_nm=0.0d0;
phi_n1=0.0d0;
phi_n0=0.0d0;

!open(unit=111,file='ini_vec.txt',status='old')
!do ii=1, Hsize
!   read(111,*) phi_n0(ii)
!enddo
!close(111)
phi_n0(1:nloct)=ini_vec(1:nloct)

carray_2=carray_2+DM_temp(1,1)*phi_n0;
Itemp1=size(DM_temp(1,:))-1
iter=1
do while(iter.le.Itemp1)


   call MPI_Barrier(comm, rc)

   !write(*,*) 'Constructing the GS wavefunction at iter:', iter, size(DM_temp(1,:))
   !!Constructing H|phi_n0>: later want to do in one step so only 3 vectors are needed.
   mycarray_1=0.0d0;
   do ii=1, spr_size
      mycarray_1(IndexJ(ii))= mycarray_1(IndexJ(ii))+sparseH(ii)*phi_n0(IndexI(ii)-lstart+1);
   enddo
   call MPI_Allreduce(mycarray_1, carray_1, Hsize, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)

   !!Calculate <phi_n0|H|phi_n0>:
   myctemp1=0.0d0;
   do ii=1, nloct
      myctemp1=myctemp1+conjg(carray_1(ii+lstart-1))*phi_n0(ii)
   enddo
   call MPI_Allreduce(myctemp1, ctemp1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)

   a_n=ctemp1;
   if(iter.eq.1) b_n=0.0d0;
   phi_n1(1:nloct)=carray_1(lstart:lend)-a_n*phi_n0(1:nloct)-b_n*phi_nm(1:nloct);

   myrtemp2=0.0d0;
   do ii=1, nloct
      myrtemp2=myrtemp2+abs(phi_n1(ii))**2
   enddo
   call MPI_Allreduce(myrtemp2, rtemp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, rc)
   phi_n1=phi_n1/sqrt(rtemp2)
   b_n=sqrt(rtemp2)

   carray_2=carray_2+DM_temp(iter+1,1)*phi_n1;

   phi_nm=phi_n0
   phi_n0=phi_n1;
   phi_n1=0.0d0;

   iter=iter+1

enddo

myrtemp1=0.0d0;
do ii=1, nloct
   myrtemp1=myrtemp1+abs(carray_2(ii))**2
enddo
call MPI_Allreduce(myrtemp1, rtemp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, rc)

write(*,*) 'H_0 normalization:', sqrt(rtemp1)
mycarray_1=0.0d0;
do ii=1, spr_size
   mycarray_1(IndexJ(ii))= mycarray_1(IndexJ(ii))+sparseH(ii)*carray_2(IndexI(ii)-lstart+1);
enddo
call MPI_Allreduce(mycarray_1, carray_1, Hsize, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)

rtemp2=0.0d0;
do ii=1, Hsize
   rtemp2=rtemp2+abs(carray_1(ii))**2
enddo

!open(unit=137,file='H_0.txt',status='unknown')
!do ii=1, Hsize
!   write(137,*) carray_2(ii)
!enddo
!close(137)

H_0(1:nloct)=carray_2(1:nloct)
call CPU_time(time2)
write(*,*) 'Lanczos finished ', time2, 'secs at iter ', iter

deallocate(carray_2, carray_1, mycarray_1, phi_n0, phi_n1, phi_nm, ini_vec)
deallocate(work, W, DM, DM_temp, E0, RWORK)

end subroutine



subroutine ARPACK(Hsize, spr_size, IndexI, IndexJ, sparseH, E_0, H_0)
use ConstantParas
implicit none

integer*8, INTENT(IN) :: Hsize, spr_size
integer*8, DIMENSION(spr_size), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(spr_size), INTENT(IN) :: sparseH
double precision, INTENT(OUT) :: E_0
double complex, DIMENSION(Hsize), INTENT(OUT) :: H_0


double precision, allocatable :: sparseHin(:)
character        bmat*1, which*2
integer     ::   maxncv, ldv
integer     ::   ido, nev, ncv, lworkl, info, ierr, j, nconv, maxitr, mode, ishfts
integer     ::   iparam(11), ipntr(11)
integer     ::   flag, tempuse1
double precision  ::   tol, sigma, zero
double precision, allocatable   ::   v(:,:), workl(:), workd(:), d(:,:)
double precision, allocatable   ::   resid(:)
logical, allocatable :: select(:)
logical     ::   rvec

integer*8 :: ii, nn

allocate(sparseHin(spr_size))
do ii=1,spr_size
   sparseHin(ii) = real(sparseH(ii)) 
enddo

nev=5 !***number of eigen-states requied.
ncv=20
maxncv=ncv
lworkl=ncv*(ncv+8)
ldv=Hsize
nn=Hsize
zero=0.0d0;

bmat = 'I'
which = 'SA'!Require the smallest algebraic eigen-values.

allocate(v(nn,ncv))
allocate(workd(3*nn))
allocate(d(ncv,2))
allocate(resid(nn))
allocate(select(ncv))
allocate(workl(lworkl))


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

      lworkl = ncv*(ncv+8)

!     tol = zero
      tol = 1.0d-15 !The tolerance factor for the convergence of Lanczos processes. If set to be zero, the error can be < 1.0d-16.
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
        workd(ipntr(2):ipntr(2)+Hsize-1)=0
        do ii=1, spr_size
        workd(ipntr(2)+IndexI(ii)-1)=workd(ipntr(2)+IndexI(ii)-1)+workd(ipntr(1)+IndexJ(ii)-1)*sparseHin(ii)
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

!       open(unit=16, file='eigenvalues.dat')
       do ii=1,nev
!         write(16,*) d(ii,1)
          write(*,*) d(ii,1)
       end do
!       close(16)

       !open(unit=17, file='eigenvectors.dat')

                do ii=1,Hsize
                        H_0(ii)=v(ii,1)
       !                 write(17,*) base(jj)
                enddo

                E_0=d(1,1);!sedenergy=d(2,1)

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
         print *, ' The number of converged Ritz values is ',&
                   nconv
         print *, ' The number of Implicit Arnoldi update',&
                 ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '

        end if

!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%

9000 continue

end subroutine

!************************************************************************
subroutine ComplexLanczos(Hsize, spr_size, IndexI, IndexJ, sparseH, E_0, H_0)
use ConstantParas
implicit none

integer*8, INTENT(IN) :: Hsize, spr_size
integer*8, DIMENSION(spr_size), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(spr_size), INTENT(IN) :: sparseH
double precision, INTENT(OUT) :: E_0
double complex, DIMENSION(Hsize), INTENT(OUT) :: H_0

integer*8 :: ii,jj,kk, i, j, k
integer*8 :: Itemp1, Itemp2, Itemp3
integer*8 :: iter, iter_num
double precision :: rtemp1, rtemp2, rtemp3
double precision :: tol, ortho_tol
double precision :: time1, time2, time3
double precision, allocatable :: rarray_1(:), rarray_2(:), rarray_3(:)
double precision, allocatable :: E0(:)
double complex :: ctemp1, ctemp2, ctemp3
double complex :: a_n, b_n
double complex, allocatable :: DM(:,:), DM_temp(:,:)
double complex, allocatable :: phi_n1(:), phi_nm(:), phi_n0(:)
double complex, allocatable :: carray_1(:), carray_2(:), carray_3(:)
!!For LAPACK:
CHARACTER*10 :: JOBZ, UPLO
integer :: LDA, INFO, LWORK
DOUBLE PRECISION, allocatable ::  RWORK(:), W(:)
double complex, allocatable :: A(:,:), WORK(: )
double complex, allocatable :: ini_vec(:)


call CPU_time(time1)

!!=============================================
!!Part I: Constructing the tri-diagonal matrix:
!!=============================================

allocate(carray_2(Hsize)); carray_2=0.0d0;
allocate(carray_1(Hsize)); carray_1=0.0d0;
allocate(phi_n0(Hsize)); phi_n0=0.0d0;
allocate(phi_n1(Hsize)); phi_n1=0.0d0;
allocate(phi_nm(Hsize)); phi_nm=0.0d0;
allocate(ini_vec(Hsize))
rtemp1=0.0d0;
do ii=1, Hsize
   call random_number(rtemp2)
   call random_number(rtemp3)
   phi_n0(ii)=dcmplx(rtemp2,rtemp3)
   rtemp1=rtemp1+abs(phi_n0(ii))**2
enddo
phi_n0=phi_n0/sqrt(rtemp1)

!!open(unit=111,file='ini_vec.txt',status='unknown')
!do ii=1,Hsize
!   write(111,*) phi_n0(ii)
!enddo
!close(111)
ini_vec(1:Hsize)=phi_n0(1:Hsize)
!!************************************************
iter_num=niter_CL;
tol=tol_CL
ortho_tol=tol_CL;
!!************************************************

!write(*,*) '*************************************************************'
!write(*,*) 'Maximum iteration nums.:', iter_num
!write(*,*) 'convergence   criterion:', tol
!write(*,*) 'Orthogonality criterion:', ortho_tol
!write(*,*) '*************************************************************'

allocate(DM(iter_num,iter_num)); DM=0.0d0;
allocate(E0(iter_num)); E0=0.0d0;
open(unit=135,file='E_0.txt',status='unknown')

iter=0; j=0;
do while(j.eq.0)

   call CPU_time(time2)
   iter=iter+1
   !write(*,*) 'ITER:', iter

   !!Constructing H|phi_n0>: later want to do in one step so only 3 vectors are needed.
   !!rarray_1=0.0d0;
   carray_1=0.0d0;
   do ii=1, spr_size
      carray_1(IndexI(ii))= carray_1(IndexI(ii))+sparseH(ii)*phi_n0(IndexJ(ii));
   enddo

   !!Calculate <phi_n0|H|phi_n0>:
   !!rtemp1=0.0d0;
   ctemp1=0.0d0;
   do ii=1, Hsize
      !rtemp1=rtemp1+rarray_1(ii)*phi_n0(ii)
      ctemp1=ctemp1+conjg(carray_1(ii))*phi_n0(ii)
   enddo

   !!a_n=rtemp1;
   a_n=ctemp1;
   if(iter.eq.1) b_n=0.0d0;
   !!phi_n1=rarray_1-a_n*phi_n0-b_n*phi_nm;

   phi_n1=carray_1-a_n*phi_n0-b_n*phi_nm;

   rtemp2=0.0d0;
   do ii=1, Hsize
      !!rtemp2=rtemp2+phi_n1(ii)**2
      rtemp2=rtemp2+abs(phi_n1(ii))**2
   enddo
   phi_n1=phi_n1/sqrt(rtemp2)
   b_n=sqrt(rtemp2)

   !!Defining the tri-diagonal matrix:
   DM(iter,iter)= a_n
   if (iter.lt.iter_num) then
      DM(iter,iter+1)= b_n
      DM(iter+1,iter)= b_n
   endif

   !write(*,*) 'ME:', a_n, b_n

   !!Calculate overlap:
   rtemp1=0.0d0;
   rtemp2=0.0d0;
   rtemp3=0.0d0;
   ctemp1=0.0d0
   ctemp2=0.0d0
   ctemp3=0.0d0

   do ii=1, Hsize
      ctemp1=ctemp1+conjg(phi_nm(ii))*phi_n0(ii)
      ctemp2=ctemp2+conjg(phi_nm(ii))*phi_n1(ii)
      ctemp3=ctemp3+conjg(phi_n0(ii))*phi_n1(ii)
   enddo
   !write(*,*) 'ortho_test:', ctemp1, ctemp2, ctemp3

   phi_nm=phi_n0
   phi_n0=phi_n1;
   phi_n1=0.0d0;

   if (allocated(DM_temp)) deallocate(DM_temp)
   allocate(DM_temp(iter,iter)); DM_temp=0.0d0;
   DM_temp(1,1)=DM(1,1)
   do ii=2, iter
      DM_temp(ii,ii)=DM(ii,ii)
      DM_temp(ii,ii-1)=DM(ii,ii-1)
      DM_temp(ii-1,ii)=DM(ii-1,ii)
   enddo

   !!For lapack:
   JOBZ='V'
   UPLO='U'

   Itemp1=iter
   LDA=Itemp1;
   LWORK=3*Itemp1;
   if (allocated(work)) deallocate(work)
   if (allocated(W)) deallocate(W)
   allocate(work(lwork)); work=0.0d0;
   allocate(W(Itemp1)); W=0.0d0;

   if(allocated(RWORK)) deallocate(RWORK)
   allocate(RWORK(lwork)); RWORK=0.0d0;

   !!SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
   !!call DSYEV(JOBZ, UPLO, Itemp1, DM_temp, LDA, W, Work, LWORK, INFO)
   !!SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
   call zheev(JOBZ, UPLO, Itemp1, DM_temp, LDA, W, WORK, LWORK, RWORK, info)
   if (info.eq.0) then
      continue
   else
      !write(*,*) 'DSYEV failed.'
      write(*,*) 'ZHEEV failed.'
      write(*,*) 'INFO:', INFO
      stop
   endif

   E0(iter)=W(1)
   !write(*,*) 'eigenvalue:', W(1)
   E_0=W(1)   !I added
   !if (iter.ge.5) then
   !   write(135,*) iter, W(1), W(2), W(3), W(4), W(5)
   !endif

   if(iter.ge.2) then
      if(abs(E0(iter)-E0(iter-1)).le.tol) j=1
   endif

   if(iter.eq.iter_num) j=3

 !  call CPU_time(time3)
 !  write(*,*) 'Lanczos iter ',iter,' takes ',time3-time2,' secs'

enddo
write(135,*) iter, W(1), W(2), W(3), W(4), W(5)
close(135)

if (j.eq.1) write(*,*) 'convergence reached at iter:', iter
if (j.eq.2) write(*,*) 'orthogonality leaks at iter:', iter
if (j.eq.3) write(*,*) 'maximum iteration reached.', iter
write(*,*) ''

!do ii=1, iter
!   write(*,*) 'DM', DM_temp(ii,1), DM_temp(ii,2)
!enddo
!write(*,*) ''

!!====================================
!!Part II: Extracting the eigenvector:
!!====================================
phi_nm=0.0d0;
phi_n1=0.0d0;
phi_n0=0.0d0;

!open(unit=111,file='ini_vec.txt',status='old')
!do ii=1, Hsize
!   read(111,*) phi_n0(ii)
!enddo
!close(111)
phi_n0(1:Hsize)=ini_vec(1:Hsize)

carray_2=carray_2+DM_temp(1,1)*phi_n0;
Itemp1=size(DM_temp(1,:))-1
iter=1
do while(iter.le.Itemp1)

   !write(*,*) 'Constructing the GS wavefunction at iter:', iter, size(DM_temp(1,:))
   !!Constructing H|phi_n0>: later want to do in one step so only 3 vectors are needed.
   carray_1=0.0d0;
   do ii=1, spr_size
      carray_1(IndexI(ii))= carray_1(IndexI(ii))+sparseH(ii)*phi_n0(IndexJ(ii));
   enddo

   !!Calculate <phi_n0|H|phi_n0>:
   ctemp1=0.0d0;
   do ii=1, Hsize
      ctemp1=ctemp1+conjg(carray_1(ii))*phi_n0(ii)
   enddo

   a_n=ctemp1;
   if(iter.eq.1) b_n=0.0d0;
   phi_n1=carray_1-a_n*phi_n0-b_n*phi_nm;

   rtemp2=0.0d0;
   do ii=1, Hsize
      rtemp2=rtemp2+abs(phi_n1(ii))**2
   enddo
   phi_n1=phi_n1/sqrt(rtemp2)
   b_n=sqrt(rtemp2)

   carray_2=carray_2+DM_temp(iter+1,1)*phi_n1;

   phi_nm=phi_n0
   phi_n0=phi_n1;
   phi_n1=0.0d0;

   iter=iter+1

enddo

rtemp1=0.0d0;
do ii=1, Hsize
   rtemp1=rtemp1+abs(carray_2(ii))**2
enddo
write(*,*) 'H_0 normalization:', sqrt(rtemp1)
carray_1=0.0d0;
do ii=1, spr_size
   carray_1(IndexI(ii))= carray_1(IndexI(ii))+sparseH(ii)*carray_2(IndexJ(ii));
enddo

rtemp2=0.0d0;
do ii=1, Hsize
   rtemp2=rtemp2+abs(carray_1(ii))**2
enddo

!open(unit=137,file='H_0.txt',status='unknown')
!do ii=1, Hsize
!   write(137,*) carray_2(ii)
!enddo
!close(137)

H_0(1:Hsize)=carray_2(1:Hsize)
call CPU_time(time2)
write(*,*) 'Lanczos finished ', time2, 'secs at iter ', iter

deallocate(carray_2, carray_1, phi_n0, phi_n1, phi_nm, ini_vec)
deallocate(work, W, DM, DM_temp, E0, RWORK)

end subroutine

!*******************************************************************

!!Cheng-Chien Chen

!!Date:02/27/2011.
!!Lanczos code serial version. Note: careful about renormalizaton of Lanczos coefficients.

!!Date:05/30/2011.
!!Modifying the Lacnzos code to perform continued fraction expansion (serial).
!!The serial version is capable of dealing with Hsize=10^7 and obtain spectra in hours.

!!Note:the input Hamiltonian matirx (sparseH.txt) and OGS.txt are now complex.

subroutine ChenCmplxContFracExpan(Hsize_3d, SprSize_3d, H_e, groundE, IndexI, IndexJ, sparseH, specX, specY)
use NumOfOrbitalAndElectrons; use ConstantParas; use ScanRegion
implicit none

double precision, INTENT(IN) :: groundE
integer*8, INTENT(IN) :: Hsize_3d, SprSize_3d
double complex, DIMENSION(Hsize_3d), INTENT(IN) :: H_e
integer*8, DIMENSION(SprSize_3d), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(SprSize_3d), INTENT(IN) :: sparseH
double precision, DIMENSION(divX+1), INTENT(OUT) :: specX
double precision, DIMENSION(divX+1), INTENT(OUT) :: specY

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

   write(*,*) iter, a_n, b_n, rtemp1
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

open(unit=1111,file='check_CFE.txt')

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

   write(1111,*) rtemp2, imag(con_ct1)/(-pi), imag(con_ct2)/(-pi)
   specX(eneX+1) = rtemp2
   specY(eneX+1) = imag(con_ct1)/(-pi)

enddo
close(1111)


call CPU_time(time2)
write(*,*) 'Complex CFE (secs):', time2-time1
deallocate(carray_1, phi_n0, phi_n1, phi_nm, work_an, work_bn)

end subroutine

!*****************************************************************************


subroutine PBiCGS(z, Hsizet, nloct, lstart, lend, &
                  SprSizet, groundH, IndexI, IndexJ, sparseH)
use NumOfOrbitalAndElectrons; use ConstantParas; use ScanRegion
use MPIParas
implicit none
include 'mpif.h'

double complex, INTENT(IN) :: z
integer*8, INTENT(IN) :: SprSizet, Hsizet
integer*8, INTENT(IN) :: nloct, lstart, lend
double complex, DIMENSION(Hsizet), INTENT(INOUT) :: groundH
integer*8, DIMENSION(SprSizet), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(SprSizet), INTENT(IN) :: sparseH

double complex, allocatable:: x0(:), xn(:)
double complex, allocatable:: r0(:), rn(:), r00(:)
double complex, allocatable:: p0(:), pn(:)
double complex, allocatable:: v0(:), vn(:)
double complex, allocatable:: s(:), at(:)
double complex, allocatable:: dvtemp(:), phi(:)
double complex :: alpha, rho0, rhon
double complex :: myrhon, mytsum, tsum, mytsum1, tsum1
double complex :: betha, omega0, omegan
integer*8 :: spr_size, ii, jj
double precision :: time1, time2
     !z = (0.0d0, 0.0d0)
   !****************        STEP 3.A        *****************
   !********              Initialization           **********

spr_size = SprSizet

allocate(x0(nloct),xn(nloct))
allocate(p0(nloct),pn(nloct))
allocate(v0(nloct),vn(nloct))
allocate(r0(nloct),rn(nloct),r00(nloct))
allocate(s(nloct),at(nloct))
allocate(dvtemp(nloct))

allocate(phi(Hsizet))

  x0(:)=0.0d0;
 ! r0(:)=groundH(:)-A(:,:)*x0(:);
  r0(1:nloct)=groundH(lstart:lend);
  r00(:)=r0(:);
  rho0=1.0d0;
  alpha=1.0d0;
  omega0=1.0d0;
  v0(:)=0.0d0;
  p0(:)=0.0d0;

!if( myid.eq.0) then
!do ii=1,nloc
!   write(*,*) groundH(ii)
!enddo
!endif

   !****************        STEP 3.B        *****************
   !*********         Conjugate Gradient     *********************

call CPU_time(time1)

  do ii=1,niter_CG  !CG loop

   myrhon = DOT_PRODUCT(r00,r0);
   rhon = (0.0d0,0.0d0)
   call MPI_Barrier(comm, rc)
   call MPI_Allreduce(myrhon, rhon, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   betha = (rhon/rho0)*(alpha/omega0);
   pn(1:nloc) = r0(1:nloc) + betha*(p0(1:nloc)-omega0*v0(1:nloc))
   call MPI_Barrier(comm, rc)
   !call MPI_Gather(pn,nloc,MPI_DOUBLE_COMPLEX,vtemp,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Gather(pn,nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(vtemp,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(phi,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)

   call MPI_Allgather(pn,nloct,MPI_DOUBLE_COMPLEX,phi,nloct,MPI_DOUBLE_COMPLEX,comm,rc)

   vn(:)=(0.0d0,0.0d0)
   do jj=1, spr_size
      vn(IndexI(jj)-lstart+1)=vn(IndexI(jj)-lstart+1)-phi(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      vn(rIndex(jj,CHsite)-localstart+1)=vn(rIndex(jj,CHsite)-localstart+1)-pn(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloct
      vn(jj)=vn(jj)+pn(jj)*z
   enddo
   mytsum = DOT_PRODUCT(r00,vn)
   call MPI_Barrier(comm, rc)
   call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   alpha = rhon/tsum
   s(1:nloct) = r0(1:nloct) - alpha*vn(1:nloct)
   at(:)=(0.0d0,0.0d0)
   call MPI_Barrier(comm, rc)
   call MPI_Allgather(s,nloct,MPI_DOUBLE_COMPLEX,phi,nloct,MPI_DOUBLE_COMPLEX,comm,rc)

   do jj=1, spr_size
      at(IndexI(jj)-lstart+1)=at(IndexI(jj)-lstart+1)-phi(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      at(rIndex(jj,CHsite)-localstart+1)=at(rIndex(jj,CHsite)-localstart+1)-s(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloct
      at(jj)=at(jj)+s(jj)*z
   enddo
   mytsum = DOT_PRODUCT(at,s)
   mytsum1 = DOT_PRODUCT(at,at)
   call MPI_Barrier(comm, rc)
   call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   call MPI_Allreduce(mytsum1, tsum1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   omegan = tsum/tsum1
   xn(:) = x0(:)+alpha*pn(:)+omegan*s(:)
   rn(:) = s(:) - omegan*at(:)

 ! xn=x0+alpha*pn+omegan*s;
  ! dvtemp(1:nloc)= x0(1:nloc)+alpha*vn(1:nloc)+omegan*at(1:nloc)-groundH(1:nloc);
   dvtemp(:)=0.0d0
   call MPI_Barrier(comm, rc)
   !call MPI_Gather(xn(1),nloc,MPI_DOUBLE_COMPLEX,vtemp,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Gather(xn(1),nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(vtemp,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(phi,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   call MPI_Allgather(xn(1),nloct,MPI_DOUBLE_COMPLEX,phi,nloct,MPI_DOUBLE_COMPLEX,comm,rc)

   do jj=1, spr_size
      dvtemp(IndexI(jj)-lstart+1)=dvtemp(IndexI(jj)-lstart+1)-phi(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      dvtemp(rIndex(jj,CHsite)-localstart+1)=dvtemp(rIndex(jj,CHsite)-localstart+1)-xn(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloct
      dvtemp(jj)=dvtemp(jj)+xn(jj)*z
   enddo
   do jj=1,nloct
      dvtemp(jj)=dvtemp(jj)-groundH(jj+lstart-1)
   enddo

   mytsum = (0.0d0,0.0d0)
   do jj=1,nloct
      mytsum = mytsum + CONJG(dvtemp(jj))*dvtemp(jj)
   enddo
   call MPI_Barrier(comm, rc)
   call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   if(myid.eq.0) write(*,*) ii, tsum
   if(real(tsum).lt.tol_CG) then
      exit
   endif

   x0(:)=xn(:)
   rho0=rhon;
   p0(:)=pn(:)
   r0(:)=rn(:)
   omega0=omegan;
   v0(:)=vn(:)

  enddo  !CG loop

call CPU_time(time2)
 !write(*,*) 'CG finishes at:', time2, 'secs for ', ii,' iters'
 if(ii.gt.niter_CG) write(*,*) 'BLOWUP TOO BAD', real(tsum)

!groundH(:)=xn(:)
call MPI_Allgather(xn(1),nloct,MPI_DOUBLE_COMPLEX,groundH,nloct,MPI_DOUBLE_COMPLEX,comm,rc)

deallocate(x0,xn,p0,pn,v0,vn,r0,rn,r00,s,at,dvtemp,phi)

end subroutine 





subroutine BiCGS(z, Hsizet, SprSizet, groundH, IndexI, IndexJ, sparseH)
use NumOfOrbitalAndElectrons; use ConstantParas; use ScanRegion
implicit none

double complex, INTENT(IN) :: z
integer*8, INTENT(IN) :: SprSizet, Hsizet
double complex, DIMENSION(Hsizet), INTENT(INOUT) :: groundH
integer*8, DIMENSION(SprSizet), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(SprSizet), INTENT(IN) :: sparseH

double complex, allocatable:: x0(:), xn(:)
double complex, allocatable:: r0(:), rn(:), r00(:)
double complex, allocatable:: p0(:), pn(:)
double complex, allocatable:: v0(:), vn(:)
double complex, allocatable:: s(:), at(:)
double complex, allocatable:: dvtemp(:)
double complex :: alpha, rho0, rhon
double complex :: myrhon, mytsum, tsum, mytsum1, tsum1
double complex :: betha, omega0, omegan
integer*8 :: spr_size, nloc, ii, jj
double precision :: time1, time2
     !z = (0.0d0, 0.0d0)
   !****************        STEP 3.A        *****************
   !********              Initialization           **********

nloc = Hsizet
spr_size = SprSizet

allocate(x0(nloc),xn(nloc))
allocate(p0(nloc),pn(nloc))
allocate(v0(nloc),vn(nloc))
allocate(r0(nloc),rn(nloc),r00(nloc))
allocate(s(nloc),at(nloc))
allocate(dvtemp(nloc))

  x0(:)=0.0d0;
 ! r0(:)=groundH(:)-A(:,:)*x0(:);
  r0(:)=groundH(:);
  r00(:)=r0(:);
  rho0=1.0d0;
  alpha=1.0d0;
  omega0=1.0d0;
  v0(:)=0.0d0;
  p0(:)=0.0d0;

!if( myid.eq.0) then
!do ii=1,nloc
!   write(*,*) groundH(ii)
!enddo
!endif

   !****************        STEP 3.B        *****************
   !*********         Conjugate Gradient     *********************

call CPU_time(time1)

  do ii=1,niter_CG  !CG loop
   rhon = DOT_PRODUCT(r00,r0);
   !rhon = (0.0d0,0.0d0)
   !call MPI_Barrier(comm, rc)
   !call MPI_Allreduce(myrhon, rhon, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   betha = (rhon/rho0)*(alpha/omega0);
   pn(1:nloc) = r0(1:nloc) + betha*(p0(1:nloc)-omega0*v0(1:nloc))
   !call MPI_Barrier(comm, rc)
   !call MPI_Gather(pn,nloc,MPI_DOUBLE_COMPLEX,vtemp,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Gather(pn,nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(vtemp,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(phi,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)

   !call MPI_Allgather(pn,nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,comm,rc)

   vn(:)=(0.0d0,0.0d0)
   do jj=1, spr_size
      vn(IndexI(jj))=vn(IndexI(jj))-pn(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      vn(rIndex(jj,CHsite)-localstart+1)=vn(rIndex(jj,CHsite)-localstart+1)-pn(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloc
      vn(jj)=vn(jj)+pn(jj)*z
   enddo
   tsum = DOT_PRODUCT(r00,vn)
   !call MPI_Barrier(comm, rc)
   !call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   alpha = rhon/tsum
   s(1:nloc) = r0(1:nloc) - alpha*vn(1:nloc)
   at(:)=(0.0d0,0.0d0)
   !call MPI_Allgather(s,nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,comm,rc)

   do jj=1, spr_size
      at(IndexI(jj))=at(IndexI(jj))-s(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      at(rIndex(jj,CHsite)-localstart+1)=at(rIndex(jj,CHsite)-localstart+1)-s(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloc
      at(jj)=at(jj)+s(jj)*z
   enddo
   tsum = DOT_PRODUCT(at,s)
   tsum1 = DOT_PRODUCT(at,at)
   !call MPI_Barrier(comm, rc)
   !call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   !call MPI_Allreduce(mytsum1, tsum1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   omegan = tsum/tsum1
   xn(:) = x0(:)+alpha*pn(:)+omegan*s(:)
   rn(:) = s(:) - omegan*at(:)

 ! xn=x0+alpha*pn+omegan*s;
  ! dvtemp(1:nloc)= x0(1:nloc)+alpha*vn(1:nloc)+omegan*at(1:nloc)-groundH(1:nloc);
   dvtemp(:)=0.0d0
   !call MPI_Barrier(comm, rc)
   !call MPI_Gather(xn(1),nloc,MPI_DOUBLE_COMPLEX,vtemp,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Gather(xn(1),nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(vtemp,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(phi,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Allgather(xn(1),nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,comm,rc)

   do jj=1, spr_size
      dvtemp(IndexI(jj))=dvtemp(IndexI(jj))-xn(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      dvtemp(rIndex(jj,CHsite)-localstart+1)=dvtemp(rIndex(jj,CHsite)-localstart+1)-xn(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloc
      dvtemp(jj)=dvtemp(jj)+xn(jj)*z
   enddo
   do jj=1,nloc
      dvtemp(jj)=dvtemp(jj)-groundH(jj)
   enddo

   tsum = (0.0d0,0.0d0)
   do jj=1,nloc
      tsum = tsum + CONJG(dvtemp(jj))*dvtemp(jj)
   enddo
   !call MPI_Barrier(comm, rc)
   !call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   write(*,*) ii, tsum
   if(real(tsum).lt.tol_CG) then
      exit
   endif

   x0(:)=xn(:)
   rho0=rhon;
   p0(:)=pn(:)
   r0(:)=rn(:)
   omega0=omegan;
   v0(:)=vn(:)

  enddo  !CG loop

call CPU_time(time2)
 write(*,*) 'CG finishes at:', time2, 'secs for ', ii,' iters'
 if(ii.gt.niter_CG) write(*,*) 'BLOWUP TOO BAD', real(tsum)

groundH(:)=xn(:)

deallocate(x0,xn,p0,pn,v0,vn,r0,rn,r00,s,at,dvtemp)

end subroutine 


!********************************************
! QPhse for 4-1D; 4-2D; 8-1D; 8-2D

!*********************************************
subroutine SetQPhase
use BettsCluster; use NumOfOrbitalAndElectrons
use ConstantParas
implicit none

double complex:: nsqrt
integer*8 :: ii,jj

selectcase(N3d)
case(4)

  if(Ndim.eq.2) then
   nQpoints = 3
!Q=(0,0)
   Qpoint(0,1)=0.0d0/1.0d0   
   Qpoint(0,2)=0.0d0/1.0d0   

   QPhase(0,0)=(1.0d0,0.0d0); QPhase(1,0)=(1.0d0,0.0d0)
   QPhase(2,0)=(1.0d0,0.0d0); QPhase(3,0)=(1.0d0,0.0d0)
!Q=(pi,pi)
   Qpoint(1,1)=1.0d0/1.0d0   
   Qpoint(1,2)=1.0d0/1.0d0   

   QPhase(0,1)=(1.0d0,0.0d0); QPhase(1,1)=(-1.0d0,0.0d0)
   QPhase(2,1)=(-1.0d0,0.0d0); QPhase(3,1)=(1.0d0,0.0d0)
!Q=(pi,0)
   Qpoint(2,1)=1.0d0/1.0d0   
   Qpoint(2,2)=0.0d0/1.0d0   

   QPhase(0,2)=(1.0d0,0.0d0); QPhase(1,2)=(1.0d0,0.0d0)
   QPhase(2,2)=(-1.0d0,0.0d0); QPhase(3,2)=(-1.0d0,0.0d0)
  else if(Ndim.eq.1) then
   nQpoints = 3
!Q=(0)
   Qpoint(0,1)=0.0d0/2.0d0   
   Qpoint(0,2)=0.0d0/2.0d0   

   QPhase(0,0)=(1.0d0,0.0d0); QPhase(1,0)=(1.0d0,0.0d0)
   QPhase(2,0)=(1.0d0,0.0d0); QPhase(3,0)=(1.0d0,0.0d0)
!Q=(pi/2)
   Qpoint(1,1)=1.0d0/2.0d0   
   Qpoint(1,2)=0.0d0/2.0d0   

   QPhase(0,1)=(1.0d0,0.0d0); QPhase(1,1)=(0.0d0,1.0d0)
   QPhase(2,1)=(-1.0d0,0.0d0); QPhase(3,1)=(0.0d0,-1.0d0)
!Q=(pi)
   Qpoint(2,1)=2.0d0/2.0d0   
   Qpoint(2,2)=0.0d0/2.0d0   

   QPhase(0,2)=(1.0d0,0.0d0); QPhase(1,2)=(-1.0d0,0.0d0)
   QPhase(2,2)=(1.0d0,0.0d0); QPhase(3,2)=(-1.0d0,0.0d0)
  endif

case(8)

  if(Ndim.eq.2) then
   nQpoints = 7
!Q=(0,0)
   Qpoint(0,1)=0.0d0/2.0d0   
   Qpoint(0,2)=0.0d0/2.0d0   

   QPhase(0,0)=(1.0d0,0.0d0); QPhase(1,0)=(1.0d0,0.0d0)
   QPhase(2,0)=(1.0d0,0.0d0); QPhase(3,0)=(1.0d0,0.0d0)
   QPhase(4,0)=(1.0d0,0.0d0); QPhase(5,0)=(1.0d0,0.0d0)
   QPhase(6,0)=(1.0d0,0.0d0); QPhase(7,0)=(1.0d0,0.0d0)
!Q=(pi,pi)
   Qpoint(1,1)=2.0d0/2.0d0   
   Qpoint(1,2)=2.0d0/2.0d0   

   QPhase(0,1)=(1.0d0,0.0d0); QPhase(1,1)=(1.0d0,0.0d0)
   QPhase(2,1)=(-1.0d0,0.0d0); QPhase(3,1)=(-1.0d0,0.0d0)
   QPhase(4,1)=(1.0d0,0.0d0); QPhase(5,1)=(1.0d0,0.0d0)
   QPhase(6,1)=(-1.0d0,0.0d0); QPhase(7,1)=(-1.0d0,0.0d0)
!Q=(pi,0)
   Qpoint(2,1)=2.0d0/2.0d0   
   Qpoint(2,2)=0.0d0/2.0d0   

   QPhase(0,2)=(1.0d0,0.0d0); QPhase(1,2)=(-1.0d0,0.0d0)
   QPhase(2,2)=(-1.0d0,0.0d0); QPhase(3,2)=(1.0d0,0.0d0)
   QPhase(4,2)=(-1.0d0,0.0d0); QPhase(5,2)=(1.0d0,0.0d0)
   QPhase(6,2)=(1.0d0,0.0d0); QPhase(7,2)=(-1.0d0,0.0d0)
!Q=(pi/2,pi/2)
   Qpoint(3,1)=1.0d0/2.0d0   
   Qpoint(3,2)=1.0d0/2.0d0   

   QPhase(0,3)=(1.0d0,0.0d0); QPhase(1,3)=(1.0d0,0.0d0)
   QPhase(2,3)=(0.0d0,1.0d0); QPhase(3,3)=(0.0d0,1.0d0)
   QPhase(4,3)=(-1.0d0,0.0d0); QPhase(5,3)=(-1.0d0,0.0d0)
   QPhase(6,3)=(0.0d0,-1.0d0); QPhase(7,3)=(0.0d0,-1.0d0)
!Q=(-pi,-pi)
   Qpoint(4,1)=-2.0d0/2.0d0   
   Qpoint(4,2)=-2.0d0/2.0d0   

   QPhase(0,4)=(1.0d0,0.0d0); QPhase(1,4)=(1.0d0,0.0d0)
   QPhase(2,4)=(-1.0d0,0.0d0); QPhase(3,4)=(-1.0d0,0.0d0)
   QPhase(4,4)=(1.0d0,0.0d0); QPhase(5,4)=(1.0d0,0.0d0)
   QPhase(6,4)=(-1.0d0,0.0d0); QPhase(7,4)=(-1.0d0,0.0d0)
!Q=(-pi/2,-pi/2)
   Qpoint(6,1)=-1.0d0/2.0d0   
   Qpoint(6,2)=-1.0d0/2.0d0   

   QPhase(0,6)=(1.0d0,0.0d0); QPhase(1,6)=(1.0d0,0.0d0)
   QPhase(2,6)=(0.0d0,-1.0d0); QPhase(3,6)=(0.0d0,-1.0d0)
   QPhase(4,6)=(-1.0d0,0.0d0); QPhase(5,6)=(-1.0d0,0.0d0)
   QPhase(6,6)=(0.0d0,1.0d0); QPhase(7,6)=(0.0d0,1.0d0)
!Q=(-pi,0)
   Qpoint(5,1)=-2.0d0/2.0d0   
   Qpoint(5,2)=-0.0d0/2.0d0   

   QPhase(0,5)=(1.0d0,0.0d0); QPhase(1,5)=(-1.0d0,0.0d0)
   QPhase(2,5)=(-1.0d0,0.0d0); QPhase(3,5)=(1.0d0,0.0d0)
   QPhase(4,5)=(-1.0d0,0.0d0); QPhase(5,5)=(1.0d0,0.0d0)
   QPhase(6,5)=(1.0d0,0.0d0); QPhase(7,5)=(-1.0d0,0.0d0)

  else if(Ndim.eq.1) then
   nQpoints = 5
!Q=(0)
   Qpoint(0,1)=0.0d0/4.0d0   
   Qpoint(0,2)=0.0d0/4.0d0   

   QPhase(0,0)=(1.0d0,0.0d0); QPhase(1,0)=(1.0d0,0.0d0)
   QPhase(2,0)=(1.0d0,0.0d0); QPhase(3,0)=(1.0d0,0.0d0)
   QPhase(4,0)=(1.0d0,0.0d0); QPhase(5,0)=(1.0d0,0.0d0)
   QPhase(6,0)=(1.0d0,0.0d0); QPhase(7,0)=(1.0d0,0.0d0)
!Q=(pi/4)
   Qpoint(1,1)=1.0d0/4.0d0   
   Qpoint(1,2)=0.0d0/4.0d0   

   QPhase(0,1)=(1.0d0,0.0d0); QPhase(1,1)=dcmplx(sqrthalf,sqrthalf)
   QPhase(2,1)=(0.0d0,1.0d0); QPhase(3,1)=dcmplx(-sqrthalf,sqrthalf)
   QPhase(4,1)=(-1.0d0,0.0d0); QPhase(5,1)=dcmplx(-sqrthalf,-sqrthalf)
   QPhase(6,1)=(0.0d0,-1.0d0); QPhase(7,1)=dcmplx(sqrthalf,-sqrthalf)
!Q=(pi/2)
   Qpoint(2,1)=2.0d0/4.0d0   
   Qpoint(2,2)=0.0d0/4.0d0   

   QPhase(0,2)=(1.0d0,0.0d0); QPhase(1,2)=(0.0d0,1.0d0)
   QPhase(2,2)=(-1.0d0,0.0d0); QPhase(3,2)=(0.0d0,-1.0d0)
   QPhase(4,2)=(1.0d0,0.0d0); QPhase(5,2)=(0.0d0,1.0d0)
   QPhase(6,2)=(-1.0d0,0.0d0); QPhase(7,2)=(0.0d0,-1.0d0)
!Q=(3*pi/4)
   Qpoint(3,1)=3.0d0/4.0d0   
   Qpoint(3,2)=0.0d0/4.0d0   

   QPhase(0,3)=(1.0d0,0.0d0); QPhase(1,3)=dcmplx(-sqrthalf,sqrthalf)
   QPhase(2,3)=(0.0d0,-1.0d0); QPhase(3,3)=dcmplx(sqrthalf,sqrthalf)
   QPhase(4,3)=(-1.0d0,0.0d0); QPhase(5,3)=dcmplx(sqrthalf,-sqrthalf)
   QPhase(6,3)=(0.0d0,1.0d0); QPhase(7,3)=dcmplx(-sqrthalf,-sqrthalf)
!Q=(pi)
   Qpoint(4,1)=4.0d0/4.0d0   
   Qpoint(4,2)=0.0d0/4.0d0   

   QPhase(0,4)=(1.0d0,0.0d0); QPhase(1,4)=(-1.0d0,0.0d0)
   QPhase(2,4)=(1.0d0,0.0d0); QPhase(3,4)=(-1.0d0,0.0d0)
   QPhase(4,4)=(1.0d0,0.0d0); QPhase(5,4)=(-1.0d0,0.0d0)
   QPhase(6,4)=(1.0d0,0.0d0); QPhase(7,4)=(-1.0d0,0.0d0)
  endif

case(10)
  
 if(Ndim.eq.1) then

   nQpoints = 6
   
   nsqrt = dcmplx(cos(2*pi/10), sin(2*pi/10))

!Q=(0)
   Qpoint(0,1)=0.0d0/5.0d0   
   Qpoint(0,2)=0.0d0/5.0d0   

   QPhase(0,0)=(1.0d0,0.0d0); QPhase(1,0)=(1.0d0,0.0d0)
   QPhase(2,0)=(1.0d0,0.0d0); QPhase(3,0)=(1.0d0,0.0d0)
   QPhase(4,0)=(1.0d0,0.0d0); QPhase(5,0)=(1.0d0,0.0d0)
   QPhase(6,0)=(1.0d0,0.0d0); QPhase(7,0)=(1.0d0,0.0d0)
   QPhase(8,0)=(1.0d0,0.0d0); QPhase(9,0)=(1.0d0,0.0d0)
!Q=(pi/5)
   Qpoint(1,1)=1.0d0/5.0d0   
   Qpoint(1,2)=0.0d0/5.0d0   

   QPhase(0,1)=(1.0d0,0.0d0); QPhase(1,1)=nsqrt
   QPhase(2,1)=nsqrt**2;      QPhase(3,1)=nsqrt**3
   QPhase(4,1)=nsqrt**4;      QPhase(5,1)=nsqrt**5
   QPhase(6,1)=nsqrt**6;      QPhase(7,1)=nsqrt**7
   QPhase(8,1)=nsqrt**8;      QPhase(9,1)=nsqrt**9
!Q=(2*pi/5)
   Qpoint(2,1)=2.0d0/5.0d0   
   Qpoint(2,2)=0.0d0/5.0d0   

   QPhase(0,2)=(1.0d0,0.0d0); QPhase(1,2)=nsqrt**2
   QPhase(2,2)=nsqrt**4;      QPhase(3,2)=nsqrt**6
   QPhase(4,2)=nsqrt**8;      QPhase(5,2)=(1.0d0,0.0d0)
   QPhase(6,2)=nsqrt**2;      QPhase(7,2)=nsqrt**4
   QPhase(8,2)=nsqrt**6;      QPhase(9,2)=nsqrt**8
!Q=(3*pi/5)
   Qpoint(3,1)=3.0d0/5.0d0   
   Qpoint(3,2)=0.0d0/5.0d0   

   QPhase(0,3)=(1.0d0,0.0d0); QPhase(1,3)=nsqrt**3
   QPhase(2,3)=nsqrt**6;      QPhase(3,3)=nsqrt**9
   QPhase(4,3)=nsqrt**2;      QPhase(5,3)=nsqrt**5
   QPhase(6,3)=nsqrt**8;      QPhase(7,3)=nsqrt**1
   QPhase(8,3)=nsqrt**4;      QPhase(9,3)=nsqrt**7
!Q=(4*pi/5)
   Qpoint(4,1)=4.0d0/5.0d0   
   Qpoint(4,2)=0.0d0/5.0d0   

   QPhase(0,4)=(1.0d0,0.0d0); QPhase(1,4)=nsqrt**4
   QPhase(2,4)=nsqrt**8;      QPhase(3,4)=nsqrt**2
   QPhase(4,4)=nsqrt**6;      QPhase(5,4)=(1.0d0,0.0d0)
   QPhase(6,4)=nsqrt**4;      QPhase(7,4)=nsqrt**8
   QPhase(8,4)=nsqrt**2;      QPhase(9,4)=nsqrt**6
!Q=(pi)
   Qpoint(5,1)=5.0d0/5.0d0   
   Qpoint(5,2)=0.0d0/5.0d0   

   QPhase(0,5)=(1.0d0,0.0d0); QPhase(1,5)=(-1.0d0,0.0d0)
   QPhase(2,5)=(1.0d0,0.0d0); QPhase(3,5)=(-1.0d0,0.0d0)
   QPhase(4,5)=(1.0d0,0.0d0); QPhase(5,5)=(-1.0d0,0.0d0)
   QPhase(6,5)=(1.0d0,0.0d0); QPhase(7,5)=(-1.0d0,0.0d0)
   QPhase(8,5)=(1.0d0,0.0d0); QPhase(9,5)=(-1.0d0,0.0d0)

 else if(Ndim.eq.2) then

write(*,*) 'set Qpoint for 10-site 2D cluster'

   nQpoints = 7

   nsqrt = dcmplx(cos(2*pi/10), sin(2*pi/10))

   !?????????????????????????????????????
   !The coordinates of all the ten points
   QA(0,1)=0.0d0; QA(0,2)=0.0d0;
   QA(1,1)=0.0d0; QA(1,2)=1.0d0;
   QA(2,1)=1.0d0; QA(2,2)=1.0d0;
   QA(3,1)=2.0d0; QA(3,2)=1.0d0;
   QA(4,1)=0.0d0; QA(4,2)=2.0d0;
   QA(5,1)=1.0d0; QA(5,2)=2.0d0;
   QA(6,1)=2.0d0; QA(6,2)=2.0d0;
   QA(7,1)=0.0d0; QA(7,2)=3.0d0;
   QA(8,1)=1.0d0; QA(8,2)=3.0d0;
   QA(9,1)=2.0d0; QA(9,2)=3.0d0;

!Q=(0,0)
   Qpoint(0,1)=0.0d0/5.0d0
   Qpoint(0,2)=0.0d0/5.0d0
!Q=(3pi/5, pi/5)
   Qpoint(1,1)=3.0d0/5.0d0
   Qpoint(1,2)=1.0d0/5.0d0
!Q=(2*pi/5, 4*pi/5)
   Qpoint(2,1)=2.0d0/5.0d0
   Qpoint(2,2)=4.0d0/5.0d0
!Q=(pi,pi)
   Qpoint(3,1)=5.0d0/5.0d0
   Qpoint(3,2)=5.0d0/5.0d0
!Q=(-3pi/5, -pi/5)
   Qpoint(4,1)=-3.0d0/5.0d0
   Qpoint(4,2)=-1.0d0/5.0d0
!Q=(-2*pi/5, -4*pi/5)
   Qpoint(5,1)=-2.0d0/5.0d0
   Qpoint(5,2)=-4.0d0/5.0d0
!Q=(-pi,-pi)
   Qpoint(6,1)=-5.0d0/5.0d0
   Qpoint(6,2)=-5.0d0/5.0d0


  do jj=0, nQpoints-1
   do ii=0,9
      QPhase(ii,jj) = exp(dcmplx(0.0d0,pi*(Qpoint(jj,1)*QA(ii,1)+ Qpoint(jj,2)*QA(ii,2))))
   enddo
  enddo

  else
   write(*,*) 'The input Betts cluster dimension is wrong!'
  endif



case(12)

 if(Ndim.eq.1) then

   nQpoints = 7

   nsqrt = dcmplx(cos(2*pi/12), sin(2*pi/12))

!Q=(0)
   Qpoint(0,1)=0   
   Qpoint(0,2)=0   

   QPhase(0,0)=(1.0d0,0.0d0); QPhase(1,0)=(1.0d0,0.0d0)
   QPhase(2,0)=(1.0d0,0.0d0); QPhase(3,0)=(1.0d0,0.0d0)
   QPhase(4,0)=(1.0d0,0.0d0); QPhase(5,0)=(1.0d0,0.0d0)
   QPhase(6,0)=(1.0d0,0.0d0); QPhase(7,0)=(1.0d0,0.0d0)
   QPhase(8,0)=(1.0d0,0.0d0); QPhase(9,0)=(1.0d0,0.0d0)
   QPhase(10,0)=(1.0d0,0.0d0);QPhase(11,0)=(1.0d0,0.0d0)
!Q=(pi/6)
   Qpoint(1,1)=1.0d0/6.0d0   
   Qpoint(1,2)=0.0d0/6.0d0   

   QPhase(0,1)=(1.0d0,0.0d0); QPhase(1,1)=nsqrt**1
   QPhase(2,1)=nsqrt**2;      QPhase(3,1)=nsqrt**3
   QPhase(4,1)=nsqrt**4;      QPhase(5,1)=nsqrt**5
   QPhase(6,1)=nsqrt**6;      QPhase(7,1)=nsqrt**7
   QPhase(8,1)=nsqrt**8;      QPhase(9,1)=nsqrt**9
   QPhase(10,1)=nsqrt**10;    QPhase(11,1)=nsqrt**11
!Q=(2*pi/6)
   Qpoint(2,1)=2.0d0/6.0d0   
   Qpoint(2,2)=0.0d0/6.0d0   

   QPhase(0,2)=(1.0d0,0.0d0); QPhase(1,2)=nsqrt**2
   QPhase(2,2)=nsqrt**4;      QPhase(3,2)=nsqrt**6
   QPhase(4,2)=nsqrt**8;      QPhase(5,2)=nsqrt**10
   QPhase(6,2)=nsqrt**0;      QPhase(7,2)=nsqrt**2
   QPhase(8,2)=nsqrt**4;      QPhase(9,2)=nsqrt**6
   QPhase(10,2)=nsqrt**8;     QPhase(11,2)=nsqrt**10
!Q=(3*pi/6)
   Qpoint(3,1)=3.0d0/6.0d0   
   Qpoint(3,2)=0.0d0/6.0d0   

   QPhase(0,3)=(1.0d0,0.0d0); QPhase(1,3)=nsqrt**3
   QPhase(2,3)=nsqrt**6;      QPhase(3,3)=nsqrt**9
   QPhase(4,3)=nsqrt**0;      QPhase(5,3)=nsqrt**3
   QPhase(6,3)=nsqrt**6;      QPhase(7,3)=nsqrt**9
   QPhase(8,3)=nsqrt**0;      QPhase(9,3)=nsqrt**3
   QPhase(10,3)=nsqrt**6;     QPhase(11,3)=nsqrt**9
!Q=(4*pi/6)
   Qpoint(4,1)=4.0d0/6.0d0   
   Qpoint(4,2)=0.0d0/6.0d0   

   QPhase(0,4)=(1.0d0,0.0d0); QPhase(1,4)=nsqrt**4
   QPhase(2,4)=nsqrt**8;      QPhase(3,4)=nsqrt**0
   QPhase(4,4)=nsqrt**4;      QPhase(5,4)=nsqrt**8
   QPhase(6,4)=nsqrt**0;      QPhase(7,4)=nsqrt**4
   QPhase(8,4)=nsqrt**8;      QPhase(9,4)=nsqrt**0
   QPhase(10,4)=nsqrt**4;     QPhase(11,4)=nsqrt**8
!Q=(5*pi/6)
   Qpoint(5,1)=5.0d0/6.0d0   
   Qpoint(5,2)=0.0d0/6.0d0   

   QPhase(0,5)=(1.0d0,0.0d0); QPhase(1,5)=nsqrt**5
   QPhase(2,5)=nsqrt**10;     QPhase(3,5)=nsqrt**3
   QPhase(4,5)=nsqrt**8;      QPhase(5,5)=nsqrt**1
   QPhase(6,5)=nsqrt**6;      QPhase(7,5)=nsqrt**11
   QPhase(8,5)=nsqrt**4;      QPhase(9,5)=nsqrt**9
   QPhase(10,5)=nsqrt**2;     QPhase(11,5)=nsqrt**7
!Q=(pi)
   Qpoint(6,1)=6.0d0/6.0d0   
   Qpoint(6,2)=0.0d0/6.0d0   

   QPhase(0,6)=(1.0d0,0.0d0); QPhase(1,6)=(-1.0d0,0.0d0)
   QPhase(2,6)=(1.0d0,0.0d0); QPhase(3,6)=(-1.0d0,0.0d0)
   QPhase(4,6)=(1.0d0,0.0d0); QPhase(5,6)=(-1.0d0,0.0d0)
   QPhase(6,6)=(1.0d0,0.0d0); QPhase(7,6)=(-1.0d0,0.0d0)
   QPhase(8,6)=(1.0d0,0.0d0); QPhase(9,6)=(-1.0d0,0.0d0)
   QPhase(10,6)=(1.0d0,0.0d0);QPhase(11,6)=(-1.0d0,0.0d0)


  else if(Ndim.eq.2) then
! 12 D cluster
! a = (4,0)
! b = (1,3)
   nQpoints = 5


!Q=(0.0d0,0.0d0)
   Qpoint(0,1)=0.0d0/6.0d0
   Qpoint(0,2)=0.0d0/6.0d0
!Q=(0.0d0,2.0d0/3*pi)
   Qpoint(1,1)=0.0d0/6.0d0
   Qpoint(1,2)=2.0d0/3.0d0
!Q=(-0.0d0,-2.0d0/3*pi)
   Qpoint(2,1)=-0.0d0/6.0d0
   Qpoint(2,2)=-2.0d0/3.0d0
!Q=(1/2*pi, 1/2*pi)
   Qpoint(3,1)=1.0d0/2.0d0
   Qpoint(3,2)=1.0d0/2.0d0
!Q=(-1/2*pi, -1/2*pi)
   Qpoint(4,1)=-1.0d0/2.0d0
   Qpoint(4,2)=-1.0d0/2.0d0
!Q=(1*pi, 1/3*pi)
   Qpoint(5,1)=1.0d0/1.0d0
   Qpoint(5,2)=1.0d0/3.0d0
!Q=(-1*pi, -1/3*pi)
   Qpoint(6,1)=-1.0d0/1.0d0
   Qpoint(6,2)=-1.0d0/3.0d0
!Q=(1*pi, 1*pi)
   Qpoint(7,1)=1.0d0/1.0d0
   Qpoint(7,2)=1.0d0/1.0d0
!Q=(-1*pi, -1*pi)
   Qpoint(8,1)=-1.0d0/1.0d0
   Qpoint(8,2)=-1.0d0/1.0d0

   !The coordinates of all the points
   QA(0,1)=0.0d0; QA(0,2)=0.0d0;
   QA(1,1)=1.0d0; QA(1,2)=0.0d0;
   QA(2,1)=2.0d0; QA(2,2)=0.0d0;
   QA(3,1)=3.0d0; QA(3,2)=0.0d0;
   QA(4,1)=1.0d0; QA(4,2)=1.0d0;
   QA(5,1)=2.0d0; QA(5,2)=1.0d0;
   QA(6,1)=3.0d0; QA(6,2)=1.0d0;
   QA(7,1)=4.0d0; QA(7,2)=1.0d0;
   QA(8,1)=1.0d0; QA(8,2)=2.0d0;
   QA(9,1)=2.0d0; QA(9,2)=2.0d0;
   QA(10,1)=3.0d0; QA(10,2)=2.0d0;
   QA(11,1)=4.0d0; QA(11,2)=2.0d0;


  do jj=0, nQpoints-1
   do ii=0,11
      QPhase(ii,jj) = exp(dcmplx(0.0d0,pi*(Qpoint(jj,1)*QA(ii,1)+ Qpoint(jj,2)*QA(ii,2))))
   enddo
  enddo


  else
   write(*,*) 'The input Betts cluster dimension is wrong!'
  endif


case(14)

 if(Ndim.eq.1) then

   nQpoints = 8

   nsqrt = dcmplx(cos(2*pi/14), sin(2*pi/14))

!Q=(0)
   Qpoint(0,1)=0

   QPhase(0,0)=(1.0d0,0.0d0); QPhase(1,0)=(1.0d0,0.0d0)
   QPhase(2,0)=(1.0d0,0.0d0); QPhase(3,0)=(1.0d0,0.0d0)
   QPhase(4,0)=(1.0d0,0.0d0); QPhase(5,0)=(1.0d0,0.0d0)
   QPhase(6,0)=(1.0d0,0.0d0); QPhase(7,0)=(1.0d0,0.0d0)
   QPhase(8,0)=(1.0d0,0.0d0); QPhase(9,0)=(1.0d0,0.0d0)
   QPhase(10,0)=(1.0d0,0.0d0);QPhase(11,0)=(1.0d0,0.0d0)
   QPhase(12,0)=(1.0d0,0.0d0);QPhase(13,0)=(1.0d0,0.0d0)
!Q=(pi/7)
   Qpoint(1,1)=1.0d0/7.0d0

   QPhase(0,1)=(1.0d0,0.0d0); QPhase(1,1)=nsqrt**1
   QPhase(2,1)=nsqrt**2;      QPhase(3,1)=nsqrt**3
   QPhase(4,1)=nsqrt**4;      QPhase(5,1)=nsqrt**5
   QPhase(6,1)=nsqrt**6;      QPhase(7,1)=nsqrt**7
   QPhase(8,1)=nsqrt**8;      QPhase(9,1)=nsqrt**9
   QPhase(10,1)=nsqrt**10;    QPhase(11,1)=nsqrt**11
   QPhase(12,1)=nsqrt**12;    QPhase(13,1)=nsqrt**13
!Q=(2*pi/7)
   Qpoint(2,1)=2.0d0/7.0d0

   QPhase(0,2)=(1.0d0,0.0d0); QPhase(1,2)=nsqrt**2
   QPhase(2,2)=nsqrt**4;      QPhase(3,2)=nsqrt**6
   QPhase(4,2)=nsqrt**8;      QPhase(5,2)=nsqrt**10
   QPhase(6,2)=nsqrt**12;      QPhase(7,2)=nsqrt**0
   QPhase(8,2)=nsqrt**2;      QPhase(9,2)=nsqrt**4
   QPhase(10,2)=nsqrt**6;     QPhase(11,2)=nsqrt**8
   QPhase(12,2)=nsqrt**10;     QPhase(13,2)=nsqrt**12
!Q=(3*pi/7)
   Qpoint(3,1)=3.0d0/7.0d0

   QPhase(0,3)=(1.0d0,0.0d0); QPhase(1,3)=nsqrt**3
   QPhase(2,3)=nsqrt**6;      QPhase(3,3)=nsqrt**9
   QPhase(4,3)=nsqrt**12;      QPhase(5,3)=nsqrt**1
   QPhase(6,3)=nsqrt**4;      QPhase(7,3)=nsqrt**7
   QPhase(8,3)=nsqrt**10;      QPhase(9,3)=nsqrt**13
   QPhase(10,3)=nsqrt**2;     QPhase(11,3)=nsqrt**5
   QPhase(12,3)=nsqrt**8;     QPhase(13,3)=nsqrt**11
!Q=(4*pi/7)
   Qpoint(4,1)=4.0d0/7.0d0

   QPhase(0,4)=(1.0d0,0.0d0); QPhase(1,4)=nsqrt**4
   QPhase(2,4)=nsqrt**8;      QPhase(3,4)=nsqrt**12
   QPhase(4,4)=nsqrt**2;      QPhase(5,4)=nsqrt**6
   QPhase(6,4)=nsqrt**10;      QPhase(7,4)=nsqrt**0
   QPhase(8,4)=nsqrt**4;      QPhase(9,4)=nsqrt**8
   QPhase(10,4)=nsqrt**12;     QPhase(11,4)=nsqrt**2
   QPhase(12,4)=nsqrt**6;     QPhase(14,4)=nsqrt**10
!Q=(5*pi/7)
   Qpoint(5,1)=5.0d0/7.0d0

   QPhase(0,5)=(1.0d0,0.0d0); QPhase(1,5)=nsqrt**5
   QPhase(2,5)=nsqrt**10;     QPhase(3,5)=nsqrt**1
   QPhase(4,5)=nsqrt**6;      QPhase(5,5)=nsqrt**11
   QPhase(6,5)=nsqrt**2;      QPhase(7,5)=nsqrt**7
   QPhase(8,5)=nsqrt**12;      QPhase(9,5)=nsqrt**3
   QPhase(10,5)=nsqrt**8;     QPhase(11,5)=nsqrt**13
   QPhase(12,5)=nsqrt**4;     QPhase(13,5)=nsqrt**9
!Q=(6*pi/7)
   Qpoint(6,1)=6.0d0/7.0d0

   QPhase(0,6)=(1.0d0,0.0d0); QPhase(1,6)=nsqrt**6
   QPhase(2,6)=nsqrt**12;     QPhase(3,6)=nsqrt**4
   QPhase(4,6)=nsqrt**10;      QPhase(5,6)=nsqrt**2
   QPhase(6,6)=nsqrt**8;      QPhase(7,6)=nsqrt**0
   QPhase(8,6)=nsqrt**6;      QPhase(9,6)=nsqrt**12
   QPhase(10,6)=nsqrt**4;     QPhase(11,6)=nsqrt**10
   QPhase(12,6)=nsqrt**2;     QPhase(13,6)=nsqrt**8
!Q=(pi)
   Qpoint(7,1)=7.0d0/7.0d0

   QPhase(0,7)=(1.0d0,0.0d0); QPhase(1,7)=(-1.0d0,0.0d0)
   QPhase(2,7)=(1.0d0,0.0d0); QPhase(3,7)=(-1.0d0,0.0d0)
   QPhase(4,7)=(1.0d0,0.0d0); QPhase(5,7)=(-1.0d0,0.0d0)
   QPhase(6,7)=(1.0d0,0.0d0); QPhase(7,7)=(-1.0d0,0.0d0)
   QPhase(8,7)=(1.0d0,0.0d0); QPhase(9,7)=(-1.0d0,0.0d0)
   QPhase(10,7)=(1.0d0,0.0d0);QPhase(11,7)=(-1.0d0,0.0d0)
   QPhase(12,7)=(1.0d0,0.0d0);QPhase(12,7)=(-1.0d0,0.0d0)

  else if(Ndim.eq.2) then
! 14 B cluster
! a = (4,0)
! b = (1,3)
   nQpoints = 11

!Q=(0.0d0,0.0d0)
   Qpoint(0,1)=0.0d0/7.0d0
   Qpoint(0,2)=0.0d0/7.0d0
!Q=(3/7pi,1/7pi)
   Qpoint(1,1)=3.0d0/7.0d0
   Qpoint(1,2)=1.0d0/7.0d0
!Q=(6/7*pi, 2/7*pi)
   Qpoint(2,1)=6.0d0/7.0d0
   Qpoint(2,2)=2.0d0/7.0d0
!Q=(1/7*pi, 5/7*pi)
   Qpoint(3,1)=1.0d0/7.0d0
   Qpoint(3,2)=5.0d0/7.0d0
!Q=(4/7*pi, 6/7*pi)
   Qpoint(4,1)=4.0d0/7.0d0
   Qpoint(4,2)=6.0d0/7.0d0
!Q=(1*pi, 1*pi)
   Qpoint(5,1)=7.0d0/7.0d0
   Qpoint(5,2)=7.0d0/7.0d0
!Q=(-3/7pi,-1/7pi)
   Qpoint(6,1)=-3.0d0/7.0d0
   Qpoint(6,2)=-1.0d0/7.0d0
!Q=(-6/7*pi, -2/7*pi)
   Qpoint(7,1)=-6.0d0/7.0d0
   Qpoint(7,2)=-2.0d0/7.0d0
!Q=(-1/7*pi, -5/7*pi)
   Qpoint(8,1)=-1.0d0/7.0d0
   Qpoint(8,2)=-5.0d0/7.0d0
!Q=(-4/7*pi, -6/7*pi)
   Qpoint(9,1)=-4.0d0/7.0d0
   Qpoint(9,2)=-6.0d0/7.0d0
!Q=(-1*pi, -1*pi)
   Qpoint(10,1)=-7.0d0/7.0d0
   Qpoint(10,2)=-7.0d0/7.0d0

   !The coordinates of all the fourteen points
   QA(0,1)=0.0d0; QA(0,2)=0.0d0;
   QA(1,1)=0.0d0; QA(1,2)=1.0d0;
   QA(2,1)=1.0d0; QA(2,2)=1.0d0;
   QA(3,1)=2.0d0; QA(3,2)=1.0d0;
   QA(4,1)=0.0d0; QA(4,2)=2.0d0;
   QA(5,1)=1.0d0; QA(5,2)=2.0d0;
   QA(6,1)=2.0d0; QA(6,2)=2.0d0;
   QA(7,1)=3.0d0; QA(7,2)=2.0d0;
   QA(8,1)=0.0d0; QA(8,2)=3.0d0;
   QA(9,1)=1.0d0; QA(9,2)=3.0d0;
   QA(10,1)=2.0d0; QA(10,2)=3.0d0;
   QA(11,1)=3.0d0; QA(11,2)=3.0d0;
   QA(12,1)=2.0d0; QA(12,2)=4.0d0;
   QA(13,1)=3.0d0; QA(13,2)=4.0d0;

  do jj=0, nQpoints-1
   do ii=0,13
      QPhase(ii,jj) = exp(dcmplx(0.0d0,pi*(Qpoint(jj,1)*QA(ii,1)+ Qpoint(jj,2)*QA(ii,2))))
   enddo
  enddo

  else
   write(*,*) 'The input Betts cluster dimension is wrong!'
  endif



endselect

end subroutine



subroutine SetPolarParas
use PolarParas; use BettsCluster; use ConstantParas

integer*8 :: qq
double precision :: alpha_p, beta_p, gamma_p
double precision :: tin, delta, theta0
double precision :: ke, a, b
!alpha_p=cos(theta_in)*cos(phi_in)
!beta_p=cos(theta_in)*sin(phi_in)
!gamma_p=-sin(theta_in)
!Gn_in(1)= sqrt(3.0d0)*alpha_p; 
!Gn_in(2)= -sqrt(3.0d0)*beta_p;
!Gn_in(3)=0.0d0

k0= 0.58d0 ! For Ein=930eV
a = 0.233d0 !sin(50)/2/(1+cos(50))
b = 3.28858 !2*(1+cos(50)) 
!propogate direction for incoming and outgoing photon


!***************************
!   Old configuration   

goto 1920

do qq=0, nQpoints-1

   tin=dasin(0.5* Qpoint(qq,1)/0.58d0+0.233153829d0*dsqrt(3.286557522-( Qpoint(qq,1)/0.58d0)**2.0))
   k_theta_in(qq)=pi-tin
   k_phi_in(qq)=pi/2
   k_theta_out(qq)=abs(50.0d0/180.0d0*pi-(pi-k_theta_in(qq)))
   if(tin.lt.50.0d0/180.0d0*pi) then
      k_phi_out(qq)=pi/2
   else 
      k_phi_out(qq)=-pi/2
   endif
   if(myid.eq.0) write(*,'(I3,3F10.4)') qq, Qpoint(qq,1), &
   0.5* Qpoint(qq,1)/0.58d0+0.233153829d0*dsqrt(3.286557522-( Qpoint(qq,1)/0.58d0)**2.0), tin
!   write(*,*) Qpoint(qq,1), k_phi_in(qq)*180/pi, k_phi_out(qq)*180/pi
!   write(*,*) '****************************************************'
enddo

!***************************



!***************************
!   new configuration for k_theta_in, k_phi_in  

1919 do qq=0, nQpoints-1

   ke=sqrt(Qpoint(qq,1)**2+Qpoint(qq,2)**2)

   tin=dacos(-a* (ke/k0) + 0.5*sqrt(b-(ke/k0)**2) )
   write(*,'(I3,4F10.4)') qq, Qpoint(qq,1), Qpoint(qq,2),tin*180/pi

   k_theta_in(qq)=pi-tin
   if(Qpoint(qq,2).eq.0.0d0) then
      k_phi_in(qq)=pi/2
   else
      k_phi_in(qq)=datan(Qpoint(qq,1)/Qpoint(qq,2))
   endif

   if(k_theta_in(qq).gt.130.0d0/180.0d0*pi) then
      k_theta_out(qq)=k_theta_in(qq)-130.0d0/180.0d0*pi
      k_phi_out(qq)=k_phi_in(qq)
   else
      k_theta_out(qq)=-k_theta_in(qq)+130.0d0/180.0d0*pi
      k_phi_out(qq)=k_phi_in(qq)+pi
   endif

enddo


1920 theta0=65*pi/180.0
 
do qq=0, nQpoints-1

   ke=sqrt(Qpoint(qq,1)**2+Qpoint(qq,2)**2)
   if((Qpoint(qq,1)+Qpoint(qq,2)).lt.0) ke=-ke

   delta=dasin(ke/(2*k0*sin(theta0)))

   if(delta.lt.(pi/2.0-theta0)) then
      k_theta_in(qq) = pi/2.0+theta0+delta
   else 
      k_theta_in(qq) = 2*pi-(pi/2.0+theta0+delta)
   endif

   if(delta.lt.(pi/2.0-theta0)) then
      if(Qpoint(qq,2).eq.0) then
         k_phi_in(qq)=pi/2
      else 
         k_phi_in(qq)=datan(Qpoint(qq,1)/Qpoint(qq,2))
      endif
   else 
      if(Qpoint(qq,2).eq.0) then
         k_phi_in(qq)=pi/2+pi
      else 
         k_phi_in(qq)=pi+datan(Qpoint(qq,1)/Qpoint(qq,2))
      endif
   endif
   
   if(delta.gt.(-pi/2.0+theta0)) then
      k_theta_out(qq) = pi/2.0-theta0+delta
   else 
      k_theta_out(qq) = -(pi/2.0-theta0+delta)
   endif

   if(delta.gt.(-pi/2.0+theta0)) then
      if(Qpoint(qq,2).eq.0) then
         k_phi_out(qq)=pi/2
      else 
         k_phi_out(qq)=datan(Qpoint(qq,1)/Qpoint(qq,2))
      endif
   else 
      if(Qpoint(qq,2).eq.0) then
         k_phi_out(qq)=pi/2+pi
      else 
         k_phi_out(qq)=pi+datan(Qpoint(qq,1)/Qpoint(qq,2))
      endif
   endif

enddo

!sigma-polarized
!polarization vector perpendicular to the incidence plane
if(PorS_in.eq.'s') then
   do qq=0, nQpoints-1
      p_theta_in(qq)=pi/2.0
      p_phi_in(qq)=k_phi_in(qq)+pi/2
   enddo
!pi-polarized: 
!polarization vector in the incidence plane
else if(PorS_in.eq.'p') then
   do qq=0, nQpoints-1
      p_theta_in(qq)=k_theta_in(qq)-pi/2
      p_phi_in(qq)= k_phi_in(qq) 
      !Even if there is 180 difference, it does not matter 
   enddo
else 
   write(*,*) 'wrong polarization input'
endif

!polarization direction for outgoing photon
if(PorS_out.eq.'s') then
   do qq=0, nQpoints
      p_theta_out(qq)=pi/2.0
      p_phi_out(qq)=k_phi_out(qq)+pi/2
   enddo
else if(PorS_out.eq.'p') then
   do qq=0, nQpoints
      p_theta_out(qq)=k_theta_out(qq)+pi/2
      p_phi_out(qq)= k_phi_out(qq)
   enddo
else 
   write(*,*) 'wrong polarization input'
endif


!Dipole operator prefactor
do qq=0,nQpoints-1        
   alpha_p=sin(p_theta_in(qq))*cos(p_phi_in(qq))
   beta_p=sin(p_theta_in(qq))*sin(p_phi_in(qq))
   !gamma_p=-cos(theta_in(qq))

   !The original Chen's congifuration
   !Gn_in(1,qq)= sqrt(3.0d0)*alpha_p; 
   !Gn_in(2,qq)= -sqrt(3.0d0)*beta_p;
   !Gn_in(3,qq)=0.0d0

   ! x <--> z'
   ! y <--> x'
   ! z <--> y'


   ! z',x',y'     x,y,z
   Gn_in(3,qq)= sqrt(3.0d0)*alpha_p; 
   Gn_in(1,qq)= -sqrt(3.0d0)*beta_p;
   Gn_in(2,qq)=0.0d0
   !Gn_in(1,qq)=1.0
   !Gn_in(2,qq)=1.0
   !Gn_in(3,qq)=1.0


   alpha_p=sin(p_theta_out(qq))*cos(p_phi_out(qq))
   beta_p=sin(p_theta_out(qq))*sin(p_phi_out(qq))
   !gamma_p=-sin(theta_out(qq))

   !The original Chen's congifuration
   !Gn_out(1,qq)= sqrt(3.0d0)*alpha_p; 
   !Gn_out(2,qq)= -sqrt(3.0d0)*beta_p;
   !Gn_out(3,qq)=0.0d0
   Gn_out(3,qq)= sqrt(3.0d0)*alpha_p; 
   Gn_out(1,qq)= -sqrt(3.0d0)*beta_p;
   Gn_out(2,qq)=0.0d0
enddo

if(myid.eq.0) then

write(*,*) ''
write(*,*) '   ky(pi) kx(pi) k_theta_in k_phi_in'  
do qq=0, nQpoints-1
   write(*,'(4F8.2)') Qpoint(qq,1), Qpoint(qq,2), k_theta_in(qq)*180/pi, k_phi_in(qq)*180/pi
enddo

write(*,*) ''
write(*,*) '   ky(pi) kx(pi) k_theta_out k_phi_out'  
do qq=0, nQpoints-1
   write(*,'(4F8.2)') Qpoint(qq,1), Qpoint(qq,2), k_theta_out(qq)*180/pi, k_phi_out(qq)*180/pi
enddo

write(*,*) ''
write(*,*) 'Polarization incoming photon: ', PorS_in
write(*,*) 'Polarization outgoing photon: ', PorS_out
write(*,*) ''
write(*,*) '   ky(pi) kx(pi) p_theta_in  p_phi_in'  
do qq=0, nQpoints-1
   write(*,'(4F8.2)') Qpoint(qq,1), Qpoint(qq,2), p_theta_in(qq)*180/pi, p_phi_in(qq)*180/pi
enddo

write(*,*) ''
write(*,*) '   ky(pi) kx(pi) p_theta_out p_phi_out'  
do qq=0, nQpoints-1
   write(*,'(4F8.2)') Qpoint(qq,1), Qpoint(qq,2), p_theta_out(qq)*180/pi, p_phi_out(qq)*180/pi
enddo
write(*,*) ''

endif

end subroutine


subroutine Setnloc(Hsizet, nprocst, myidt, nloct, localstartt, localendt)
implicit none
INTEGER*8, INTENT(IN) :: Hsizet
INTEGER, INTENT(IN) :: nprocst, myidt
INTEGER*8, INTENT(OUT) :: nloct, localstartt, localendt

if(mod(Hsizet,nprocst).eq.0) then
   nloct = Hsizet/nprocst
else
   nloct = Hsizet/nprocst + 1
endif
localstartt = 1+nloct*myidt
localendt = nloct*(myidt+1)
if(myidt.eq.nprocst-1) then
   localendt=Hsizet
   nloct=localendt-localstartt+1
endif

end subroutine

