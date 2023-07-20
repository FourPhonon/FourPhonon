!  FourPhonon: An extension module to ShengBTE for computing four-phonon scattering rates and thermal conductivity
!  Copyright (C) 2021 Zherui Han <zrhan@purdue.edu>
!  Copyright (C) 2021 Xiaolong Yang <xiaolongyang1990@gmail.com>
!  Copyright (C) 2021 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2021 Tianli Feng <Tianli.Feng2011@gmail.com>
!  Copyright (C) 2021 Xiulin Ruan <ruan@purdue.edu>
!
!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2017 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2017 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2017 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2017 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.


! Compute the number of allowed three-phonon / four-phonon processes, their
! scattering amplitudes and their phase-space volume.

module processes
  use iso_fortran_env
  use misc
  use data
  use config
  !use omp_lib
  implicit none

  real(kind=8),parameter :: hbarp=hbar*1e22

contains

  ! Compute one of the matrix elements involved in the calculation of Ind_plus.
  function Vp_plus(i,j,k,q,qprime,qdprime,realqprime,realqdprime,eigenvect,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
    implicit none

    integer(kind=4),intent(in) :: i
    integer(kind=4),intent(in) :: j
    integer(kind=4),intent(in) :: k
    integer(kind=4),intent(in) :: q
    integer(kind=4),intent(in) :: qprime
    integer(kind=4),intent(in) :: qdprime
    real(kind=8),intent(in) :: realqprime(3)
    real(kind=8),intent(in) :: realqdprime(3)    
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(in) :: Ntri
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri)
    real(kind=8),intent(in) :: R_j(3,Ntri)
    real(kind=8),intent(in) :: R_k(3,Ntri)
    integer(kind=4),intent(in) :: Index_i(Ntri)
    integer(kind=4),intent(in) :: Index_j(Ntri)
    integer(kind=4),intent(in) :: Index_k(Ntri)

    complex(kind=8) :: Vp_plus

    integer(kind=4) :: ll
    integer(kind=4) :: rr
    integer(kind=4) :: ss
    integer(kind=4) :: tt
    complex(kind=8) :: prefactor
    complex(kind=8) :: Vp0

    Vp_plus=0.d0
    
    do ll=1,Ntri
       prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*&
            masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
            phexp(dot_product(realqprime,R_j(:,ll)))*&
            phexp(-dot_product(realqdprime,R_k(:,ll)))
       Vp0=0.d0
       do rr=1,3
          do ss=1,3
             do tt=1,3
                Vp0=Vp0+Phi(tt,ss,rr,ll)*&
                     eigenvect(q,i,tt+3*(Index_i(ll)-1))*&
                     eigenvect(qprime,j,ss+3*(Index_j(ll)-1))*&
                     conjg(eigenvect(qdprime,k,rr+3*(Index_k(ll)-1)))
             end do
          end do
       end do
       Vp_plus=Vp_plus+prefactor*Vp0
    end do
  end function Vp_plus

  ! Compute one of the matrix elements involved in the calculation of Ind_minus.
  function Vp_minus(i,j,k,q,qprime,qdprime,realqprime,realqdprime,eigenvect,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
    implicit none

    integer(kind=4),intent(in) :: i
    integer(kind=4),intent(in) :: j
    integer(kind=4),intent(in) :: k
    integer(kind=4),intent(in) :: q
    integer(kind=4),intent(in) :: qprime
    integer(kind=4),intent(in) :: qdprime
    real(kind=8),intent(in) :: realqprime(3)
    real(kind=8),intent(in) :: realqdprime(3)    
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(in) :: Ntri
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri)
    real(kind=8),intent(in) :: R_j(3,Ntri)
    real(kind=8),intent(in) :: R_k(3,Ntri)
    integer(kind=4),intent(in) :: Index_i(Ntri)
    integer(kind=4),intent(in) :: Index_j(Ntri)
    integer(kind=4),intent(in) :: Index_k(Ntri)

    complex(kind=8) :: Vp_minus

    integer(kind=4) :: ll
    integer(kind=4) :: rr
    integer(kind=4) :: ss
    integer(kind=4) :: tt
    complex(kind=8) :: prefactor
    complex(kind=8) :: Vp0

    Vp_minus=0.d0

    do ll=1,Ntri
       prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*&
            masses(types(Index_j(ll)))*masses(types(Index_k(ll))))*&
            phexp(-dot_product(realqprime,R_j(:,ll)))*&
            phexp(-dot_product(realqdprime,R_k(:,ll)))
       Vp0=0.
       do rr=1,3
          do ss=1,3
             do tt=1,3
                Vp0=Vp0+Phi(tt,ss,rr,ll)*&
                     eigenvect(q,i,tt+3*(Index_i(ll)-1))*&
                     conjg(eigenvect(qprime,j,ss+3*(Index_j(ll)-1)))*&
                     conjg(eigenvect(qdprime,k,rr+3*(Index_k(ll)-1)))
             end do
          end do
       end do
       Vp_minus=Vp_minus+prefactor*Vp0
    end do
  end function Vp_minus

  ! Compute the matrix elements for 4ph ++ process, Ind_plusplus
  function Vp_pp(i,j,k,l,q,qprime,qdprime,qtprime,realqprime,realqdprime,realqtprime,eigenvect,&
         Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
      implicit none

      integer(kind=4),intent(in) :: i
      integer(kind=4),intent(in) :: j
      integer(kind=4),intent(in) :: k
      integer(kind=4),intent(in) :: l
      integer(kind=4),intent(in) :: q
      integer(kind=4),intent(in) :: qprime
      integer(kind=4),intent(in) :: qdprime
      integer(kind=4),intent(in) :: qtprime
      real(kind=8),intent(in) :: realqprime(3)
      real(kind=8),intent(in) :: realqdprime(3) 
      real(kind=8),intent(in) :: realqtprime(3)   
      complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
      integer(kind=4),intent(in) :: Ntri
      real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri)
      real(kind=8),intent(in) :: R_j(3,Ntri)
      real(kind=8),intent(in) :: R_k(3,Ntri)
      real(kind=8),intent(in) :: R_l(3,Ntri)
      integer(kind=4),intent(in) :: Index_i(Ntri)
      integer(kind=4),intent(in) :: Index_j(Ntri)
      integer(kind=4),intent(in) :: Index_k(Ntri)
      integer(kind=4),intent(in) :: Index_l(Ntri)

      complex(kind=8) :: Vp_pp

      integer(kind=4) :: ll
      integer(kind=4) :: rr
      integer(kind=4) :: ss
      integer(kind=4) :: tt
      integer(kind=4) :: uu
      complex(kind=8) :: prefactor
      complex(kind=8) :: Vp0

      Vp_pp=0.d0

      do ll=1,Ntri
         prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*masses(types(Index_j(ll)))*&
               masses(types(Index_k(ll)))*masses(types(Index_l(ll))))*&
               phexp(dot_product(realqprime,R_j(:,ll)))*&
               phexp(dot_product(realqdprime,R_k(:,ll)))*&
               phexp(-dot_product(realqtprime,R_l(:,ll)))
         Vp0=0.d0
         do rr=1,3
            do ss=1,3
               do tt=1,3
                  do uu=1,3
                     Vp0=Vp0+Phi(uu,tt,ss,rr,ll)*&
                        eigenvect(q,i,uu+3*(Index_i(ll)-1))*&
                        eigenvect(qprime,j,tt+3*(Index_j(ll)-1))*&
                        eigenvect(qdprime,k,ss+3*(Index_k(ll)-1))*&
                        conjg(eigenvect(qtprime,l,rr+3*(Index_l(ll)-1)))
                  end do
               end do
            end do
         end do
         Vp_pp=Vp_pp+prefactor*Vp0
      end do

  end function Vp_pp

  ! Compute the matrix elements for 4ph +- process, Ind_plusminus
  function Vp_pm(i,j,k,l,q,qprime,qdprime,qtprime,realqprime,realqdprime,realqtprime,eigenvect,&
         Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
      implicit none

      integer(kind=4),intent(in) :: i
      integer(kind=4),intent(in) :: j
      integer(kind=4),intent(in) :: k
      integer(kind=4),intent(in) :: l
      integer(kind=4),intent(in) :: q
      integer(kind=4),intent(in) :: qprime
      integer(kind=4),intent(in) :: qdprime
      integer(kind=4),intent(in) :: qtprime
      real(kind=8),intent(in) :: realqprime(3)
      real(kind=8),intent(in) :: realqdprime(3) 
      real(kind=8),intent(in) :: realqtprime(3)   
      complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
      integer(kind=4),intent(in) :: Ntri
      real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri)
      real(kind=8),intent(in) :: R_j(3,Ntri)
      real(kind=8),intent(in) :: R_k(3,Ntri)
      real(kind=8),intent(in) :: R_l(3,Ntri)
      integer(kind=4),intent(in) :: Index_i(Ntri)
      integer(kind=4),intent(in) :: Index_j(Ntri)
      integer(kind=4),intent(in) :: Index_k(Ntri)
      integer(kind=4),intent(in) :: Index_l(Ntri)

      complex(kind=8) :: Vp_pm

      integer(kind=4) :: ll
      integer(kind=4) :: rr
      integer(kind=4) :: ss
      integer(kind=4) :: tt
      integer(kind=4) :: uu
      complex(kind=8) :: prefactor
      complex(kind=8) :: Vp0

      Vp_pm=0.d0

      do ll=1,Ntri
         prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*masses(types(Index_j(ll)))*&
               masses(types(Index_k(ll)))*masses(types(Index_l(ll))))*&
               phexp(dot_product(realqprime,R_j(:,ll)))*&
               phexp(-dot_product(realqdprime,R_k(:,ll)))*&
               phexp(-dot_product(realqtprime,R_l(:,ll)))
         Vp0=0.d0
         do rr=1,3
            do ss=1,3
               do tt=1,3
                  do uu=1,3
                     Vp0=Vp0+Phi(uu,tt,ss,rr,ll)*&
                        eigenvect(q,i,uu+3*(Index_i(ll)-1))*&
                        eigenvect(qprime,j,tt+3*(Index_j(ll)-1))*&
                        conjg(eigenvect(qdprime,k,ss+3*(Index_k(ll)-1)))*&
                        conjg(eigenvect(qtprime,l,rr+3*(Index_l(ll)-1)))
                  end do
               end do
            end do
         end do
         Vp_pm=Vp_pm+prefactor*Vp0
      end do

  end function Vp_pm

  ! Compute the matrix elements for 4ph -- process, Ind_minusminus
  function Vp_mm(i,j,k,l,q,qprime,qdprime,qtprime,realqprime,realqdprime,realqtprime,eigenvect,&
      Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
   implicit none

   integer(kind=4),intent(in) :: i
   integer(kind=4),intent(in) :: j
   integer(kind=4),intent(in) :: k
   integer(kind=4),intent(in) :: l
   integer(kind=4),intent(in) :: q
   integer(kind=4),intent(in) :: qprime
   integer(kind=4),intent(in) :: qdprime
   integer(kind=4),intent(in) :: qtprime
   real(kind=8),intent(in) :: realqprime(3)
   real(kind=8),intent(in) :: realqdprime(3) 
   real(kind=8),intent(in) :: realqtprime(3)   
   complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
   integer(kind=4),intent(in) :: Ntri
   real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri)
   real(kind=8),intent(in) :: R_j(3,Ntri)
   real(kind=8),intent(in) :: R_k(3,Ntri)
   real(kind=8),intent(in) :: R_l(3,Ntri)
   integer(kind=4),intent(in) :: Index_i(Ntri)
   integer(kind=4),intent(in) :: Index_j(Ntri)
   integer(kind=4),intent(in) :: Index_k(Ntri)
   integer(kind=4),intent(in) :: Index_l(Ntri)

   complex(kind=8) :: Vp_mm

   integer(kind=4) :: ll
   integer(kind=4) :: rr
   integer(kind=4) :: ss
   integer(kind=4) :: tt
   integer(kind=4) :: uu
   complex(kind=8) :: prefactor
   complex(kind=8) :: Vp0

   Vp_mm=0.d0

   do ll=1,Ntri
      prefactor=1.d0/sqrt(masses(types(Index_i(ll)))*masses(types(Index_j(ll)))*&
            masses(types(Index_k(ll)))*masses(types(Index_l(ll))))*&
            phexp(-dot_product(realqprime,R_j(:,ll)))*&
            phexp(-dot_product(realqdprime,R_k(:,ll)))*&
            phexp(-dot_product(realqtprime,R_l(:,ll)))
      Vp0=0.d0
      do rr=1,3
         do ss=1,3
            do tt=1,3
               do uu=1,3
                  Vp0=Vp0+Phi(uu,tt,ss,rr,ll)*&
                     eigenvect(q,i,uu+3*(Index_i(ll)-1))*&
                     conjg(eigenvect(qprime,j,tt+3*(Index_j(ll)-1)))*&
                     conjg(eigenvect(qdprime,k,ss+3*(Index_k(ll)-1)))*&
                     conjg(eigenvect(qtprime,l,rr+3*(Index_l(ll)-1)))
               end do
            end do
         end do
      end do
      Vp_mm=Vp_mm+prefactor*Vp0
   end do
  end function Vp_mm


  ! Use the 1st BZ image of each q point to distinguish normal and umklapp  processes
  subroutine ws_cell(q,shortest)
    implicit none

    integer(kind=4),intent(in) :: q(3)
    real(kind=8) :: r(3),dnrm2,tmp1,tmp2
    integer(kind=4) :: ix1,iy1,iz1
    real(kind=8),intent(out) :: shortest(3)

    shortest=matmul(rlattvec,q/dble(ngrid))
    tmp1=dnrm2(3,shortest,1)
    do ix1=-3,3
       do iy1=-3,3
          do iz1=-3,3
             r=matmul(rlattvec,q/dble(ngrid))+ix1*rlattvec(:,1)+iy1*rlattvec(:,2)+&
                     iz1*rlattvec(:,3)
             tmp2=dnrm2(3,r,1)
             if(tmp2.lt.tmp1) then
               tmp1=tmp2
               shortest=r
             end if
          end do
       end do
    end do
  end subroutine ws_cell


  ! Scattering amplitudes of absorption processes.
  subroutine Ind_plus(mm,N_plus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
       Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus,WP3_plus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_plus,Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(out) :: Indof2ndPhonon_plus(N_plus),Indof3rdPhonon_plus(N_plus)
    real(kind=8),intent(out) :: Gamma_plus(N_plus),WP3_plus

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_plus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime
    real(kind=8) :: omega,omegap,omegadp
    real(kind=8) :: realqprime(3),realqdprime(3)
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3)
    real(kind=8) :: WP3
    complex(kind=8) :: Vp

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus_count=0
    WP3_plus=0.d0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(list(ll),i)
    ! Loop over all processes, detecting those that are allowed and
    ! computing their amplitudes.
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             realqprime=matmul(rlattvec,qprime/dble(ngrid))
             omegap=energy(ii,j)
             fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
                ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
                omegadp=energy(ss,k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(ii,j,:)-&
                        velocity(ss,k,:))
                   if(abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                   !-------- call ws_cell function to distinguish normal and umklapp processes
                   call ws_cell(q,shortest_q)
                   call ws_cell(qprime,shortest_qprime)
                   call ws_cell(qdprime,shortest_qdprime)     
                   !  write(*,"(1E20.10)") shortest_q(1)+shortest_qprime(1)-shortest_qdprime(1)
                   !  write(*,"(1E20.10)") shortest_qdprime(1)-realqdprime(1)
                   !----------------------Normal process condition
                   if (normal) then
                   if (abs(shortest_q(1)+shortest_qprime(1)-shortest_qdprime(1))<1e-10.and.&
                       abs(shortest_q(2)+shortest_qprime(2)-shortest_qdprime(2))<1e-10.and.&
                       abs(shortest_q(3)+shortest_qprime(3)-shortest_qdprime(3))<1e-10) then
                   !-----------------------------------------------------------
                      N_plus_count=N_plus_count+1
                      Indof2ndPhonon_plus(N_plus_count)=(ii-1)*Nbands+j
                      Indof3rdPhonon_plus(N_plus_count)=(ss-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      Vp=Vp_plus(i,j,k,list(ll),ii,ss,&
                           realqprime,realqdprime,eigenvect,&
                           Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                      WP3=(fBEprime-fBEdprime)*&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_plus=WP3_plus+WP3
                      Gamma_plus(N_plus_count)=hbarp*pi/4.d0*WP3*abs(Vp)**2
                      ! At this point, Gamma's units are
                      ! (1.d-34J*s)*(1.d12/s)^(-4)*1amu^(-3)*(ev/angstrom**3)^2,
                      ! that is, 5.60626442*1.d8 THz
                      Gamma_plus(N_plus_count)=Gamma_plus(N_plus_count)*5.60626442*1.d8/nptk ! THz
                   end if
                   !---------------------Umklapp process condition
                   else if (umklapp) then
                   if (abs(shortest_q(1)+shortest_qprime(1)-shortest_qdprime(1))>1e-10.or.&
                       abs(shortest_q(2)+shortest_qprime(2)-shortest_qdprime(2))>1e-10.or.&
                       abs(shortest_q(3)+shortest_qprime(3)-shortest_qdprime(3))>1e-10) then
                   !-----------------------------------------------------------
                      N_plus_count=N_plus_count+1
                      Indof2ndPhonon_plus(N_plus_count)=(ii-1)*Nbands+j
                      Indof3rdPhonon_plus(N_plus_count)=(ss-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      Vp=Vp_plus(i,j,k,list(ll),ii,ss,&
                           realqprime,realqdprime,eigenvect,&
                           Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                      WP3=(fBEprime-fBEdprime)*&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_plus=WP3_plus+WP3
                      Gamma_plus(N_plus_count)=hbarp*pi/4.d0*WP3*abs(Vp)**2
                      ! At this point, Gamma's units are
                      ! (1.d-34J*s)*(1.d12/s)^(-4)*1amu^(-3)*(ev/angstrom**3)^2,
                      ! that is, 5.60626442*1.d8 THz
                      Gamma_plus(N_plus_count)=Gamma_plus(N_plus_count)*5.60626442*1.d8/nptk ! THz
                   end if
                   !----------------------------total scattering process
                   else
                      N_plus_count=N_plus_count+1
                      Indof2ndPhonon_plus(N_plus_count)=(ii-1)*Nbands+j
                      Indof3rdPhonon_plus(N_plus_count)=(ss-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      Vp=Vp_plus(i,j,k,list(ll),ii,ss,&
                           realqprime,realqdprime,eigenvect,&
                           Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                      WP3=(fBEprime-fBEdprime)*&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_plus=WP3_plus+WP3
                      Gamma_plus(N_plus_count)=hbarp*pi/4.d0*WP3*abs(Vp)**2
                      ! At this point, Gamma's units are
                      ! (1.d-34J*s)*(1.d12/s)^(-4)*1amu^(-3)*(ev/angstrom**3)^2,
                      ! that is, 5.60626442*1.d8 THz
                      Gamma_plus(N_plus_count)=Gamma_plus(N_plus_count)*5.60626442*1.d8/nptk ! THz
                   end if
                   end if
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
       WP3_plus=WP3_plus/nptk
    end if
  end subroutine Ind_plus

  ! Scattering amplitudes of emission processes. See Ind_plus() for details.
  subroutine Ind_minus(mm,N_minus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
       Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus,WP3_minus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),N_minus,Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(out) :: Indof2ndPhonon_minus(N_minus),Indof3rdPhonon_minus(N_minus)
    real(kind=8),intent(out) :: Gamma_minus(N_minus),WP3_minus

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_minus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime
    real(kind=8) ::  omega,omegap,omegadp
    real(kind=8) :: realqprime(3),realqdprime(3)
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3)
    real(kind=8) :: WP3
    complex(kind=8) :: Vp

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus_count=0
    WP3_minus=0.d0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(list(ll),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             realqprime=matmul(rlattvec,qprime/dble(ngrid))
             omegap=energy(ii,j)
             fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                qdprime=q-qprime
                qdprime=modulo(qdprime,Ngrid)
                realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
                ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
                omegadp=energy(ss,k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(ii,j,:)-&
                        velocity(ss,k,:))
                   if (abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                   !-------- call ws_cell function to distinguish normal and umklapp processes
                   call ws_cell(q,shortest_q)
                   call ws_cell(qprime,shortest_qprime)
                   call ws_cell(qdprime,shortest_qdprime)
                   !----------------------Normal process condition
                   if (normal) then
                   if (abs(shortest_q(1)-shortest_qprime(1)-shortest_qdprime(1))<1e-10.and.&
                       abs(shortest_q(2)-shortest_qprime(2)-shortest_qdprime(2))<1e-10.and.&
                       abs(shortest_q(3)-shortest_qprime(3)-shortest_qdprime(3))<1e-10) then
                    !-----------------------------------------------------------------------------
                      N_minus_count=N_minus_count+1
                      Indof2ndPhonon_minus(N_minus_count)=(ii-1)*Nbands+j
                      Indof3rdPhonon_minus(N_minus_count)=(ss-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      Vp=Vp_minus(i,j,k,list(ll),ii,ss,&
                           realqprime,realqdprime,eigenvect,&
                           Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                      WP3=(fBEprime+fBEdprime+1)*&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_minus=WP3_minus+WP3
                      Gamma_minus(N_minus_count)=hbarp*pi/4.d0*WP3*abs(Vp)**2
                      Gamma_minus(N_minus_count)=Gamma_minus(N_minus_count)*5.60626442*1.d8/nptk
                   end if
                   !--------------------Umklapp process condition
                   else if (umklapp) then
                   if (abs(shortest_q(1)-shortest_qprime(1)-shortest_qdprime(1))>1e-10.or.&
                       abs(shortest_q(2)-shortest_qprime(2)-shortest_qdprime(2))>1e-10.or.&
                       abs(shortest_q(3)-shortest_qprime(3)-shortest_qdprime(3))>1e-10) then
                    !-----------------------------------------------------------------------------
                      N_minus_count=N_minus_count+1
                      Indof2ndPhonon_minus(N_minus_count)=(ii-1)*Nbands+j
                      Indof3rdPhonon_minus(N_minus_count)=(ss-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      Vp=Vp_minus(i,j,k,list(ll),ii,ss,&
                           realqprime,realqdprime,eigenvect,&
                           Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                      WP3=(fBEprime+fBEdprime+1)*&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_minus=WP3_minus+WP3
                      Gamma_minus(N_minus_count)=hbarp*pi/4.d0*WP3*abs(Vp)**2
                      Gamma_minus(N_minus_count)=Gamma_minus(N_minus_count)*5.60626442*1.d8/nptk
                   end if
                   !------------------total scattering process
                   else
                      N_minus_count=N_minus_count+1
                      Indof2ndPhonon_minus(N_minus_count)=(ii-1)*Nbands+j
                      Indof3rdPhonon_minus(N_minus_count)=(ss-1)*Nbands+k
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      Vp=Vp_minus(i,j,k,list(ll),ii,ss,&
                           realqprime,realqdprime,eigenvect,&
                           Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                      WP3=(fBEprime+fBEdprime+1)*&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_minus=WP3_minus+WP3
                      Gamma_minus(N_minus_count)=hbarp*pi/4.d0*WP3*abs(Vp)**2
                      Gamma_minus(N_minus_count)=Gamma_minus(N_minus_count)*5.60626442*1.d8/nptk
                   end if
                   end if
                end if
             end do ! k
             !--------END emission process-------------
          end do ! ii
       end do  ! j
       WP3_minus=WP3_minus*5.d-1/nptk
    end if
  end subroutine Ind_minus

  ! Wrapper around Ind_plus and Ind_minus that splits the work among processors.
  subroutine Ind_driver(energy,velocity,eigenvect,Nlist,List,IJK,N_plus,N_minus,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,&
       Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus,&
       Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus,rate_scatt,rate_scatt_plus,rate_scatt_minus,WP3_plus,WP3_minus)
    implicit none

    include "mpif.h"

    real(kind=8),intent(in) :: energy(nptk,nbands)
    real(kind=8),intent(in) :: velocity(nptk,nbands,3)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(in) :: NList
    integer(kind=4),intent(in) :: List(Nlist)
    integer(kind=4),intent(in) :: IJK(3,nptk)
    integer(kind=4),intent(in) :: N_plus(Nlist*Nbands)
    integer(kind=4),intent(in) :: N_minus(Nlist*Nbands)
    integer(kind=4),intent(in) :: Ntri
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri)
    real(kind=8),intent(in) :: R_j(3,Ntri)
    real(kind=8),intent(in) :: R_k(3,Ntri)
    integer(kind=4),intent(in) :: Index_i(Ntri)
    integer(kind=4),intent(in) :: Index_j(Ntri)
    integer(kind=4),intent(in) :: Index_k(Ntri)
    integer(kind=4),intent(out) :: Indof2ndPhonon_plus(:)
    integer(kind=4),intent(out) :: Indof3rdPhonon_plus(:)
    real(kind=8),intent(out) :: Gamma_plus(:)
    integer(kind=4),intent(out) :: Indof2ndPhonon_minus(:)
    integer(kind=4),intent(out) :: Indof3rdPhonon_minus(:)
    real(kind=8),intent(out) :: Gamma_minus(:)
    real(kind=8),intent(out) :: rate_scatt(Nbands,Nlist),rate_scatt_plus(Nbands,Nlist),rate_scatt_minus(Nbands,Nlist)
    real(kind=8),intent(out) :: WP3_plus(Nbands,Nlist)
    real(kind=8),intent(out) :: WP3_minus(Nbands,Nlist)

    integer(kind=4) :: i
    integer(kind=4) :: ll
    integer(kind=4) :: mm
    integer(kind=4) :: maxsize
    integer(kind=4) :: Ntotal_plus
    integer(kind=4) :: Ntotal_minus
    integer(kind=4) :: Naccum_plus(Nbands*Nlist)
    integer(kind=4) :: Naccum_minus(Nbands*Nlist)
    integer(kind=4),allocatable :: Indof2ndPhonon(:)
    integer(kind=4),allocatable :: Indof3rdPhonon(:)
    integer(kind=4),allocatable :: Indof2ndPhonon_plus_reduce(:)
    integer(kind=4),allocatable :: Indof3rdPhonon_plus_reduce(:)
    integer(kind=4),allocatable :: Indof2ndPhonon_minus_reduce(:)
    integer(kind=4),allocatable :: Indof3rdPhonon_minus_reduce(:)
    real(kind=8) :: rate_scatt_plus_reduce(Nbands,Nlist),rate_scatt_minus_reduce(Nbands,Nlist)
    real(kind=8),allocatable :: Gamma0(:)
    real(kind=8),allocatable :: Gamma_plus_reduce(:)
    real(kind=8),allocatable :: Gamma_minus_reduce(:)
    real(kind=8) :: WP3_plus_reduce(Nbands*Nlist)
    real(kind=8) :: WP3_minus_reduce(Nbands*Nlist)

    maxsize=max(maxval(N_plus),maxval(N_minus))
    allocate(Indof2ndPhonon(maxsize))
    allocate(Indof3rdPhonon(maxsize))
    allocate(Gamma0(maxsize))

    Naccum_plus(1)=0
    Naccum_minus(1)=0
    do mm=2,Nbands*Nlist
       Naccum_plus(mm)=Naccum_plus(mm-1)+N_plus(mm-1)
       Naccum_minus(mm)=Naccum_minus(mm-1)+N_minus(mm-1)
    end do
    Ntotal_plus=sum(N_plus)
    Ntotal_minus=sum(N_minus)

    allocate(Indof2ndPhonon_plus_reduce(Ntotal_plus))
    allocate(Indof3rdPhonon_plus_reduce(Ntotal_plus))
    allocate(Indof2ndPhonon_minus_reduce(Ntotal_minus))
    allocate(Indof3rdPhonon_minus_reduce(Ntotal_minus))
    allocate(Gamma_plus_reduce(Ntotal_plus))
    allocate(Gamma_minus_reduce(Ntotal_minus))

    Indof2ndPhonon_plus_reduce=0
    Indof3rdPhonon_plus_reduce=0
    Indof2ndPhonon_minus_reduce=0
    Indof3rdPhonon_minus_reduce=0
    Gamma_plus_reduce=0.d0
    Gamma_minus_reduce=0.d0
    rate_scatt_plus_reduce=0.d0
    rate_scatt_minus_reduce=0.d0
    WP3_plus=0.d0
    WP3_minus=0.d0
    WP3_plus_reduce=0.d0
    WP3_minus_reduce=0.d0

    do mm=myid+1,Nbands*NList,numprocs
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       if(N_plus(mm).ne.0) then
          call Ind_plus(mm,N_plus(mm),energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
               Indof2ndPhonon(1:N_plus(mm)),Indof3rdPhonon(1:N_plus(mm)),&
               Gamma0(1:N_plus(mm)),WP3_plus_reduce(mm))
          Indof2ndPhonon_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Indof2ndPhonon(1:N_plus(mm))
          Indof3rdPhonon_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Indof3rdPhonon(1:N_plus(mm))
          Gamma_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Gamma0(1:N_plus(mm))
          rate_scatt_plus_reduce(i,ll)=sum(Gamma0(1:N_plus(mm)))
       end if
       if(N_minus(mm).ne.0) then
          call Ind_minus(mm,N_minus(mm),energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
               Indof2ndPhonon(1:N_minus(mm)),Indof3rdPhonon(1:N_minus(mm)),&
               Gamma0(1:N_minus(mm)),WP3_minus_reduce(mm))
          Indof2ndPhonon_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Indof2ndPhonon(1:N_minus(mm))
          Indof3rdPhonon_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Indof3rdPhonon(1:N_minus(mm))
          Gamma_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Gamma0(1:N_minus(mm))
          rate_scatt_minus_reduce(i,ll)=sum(Gamma0(1:N_minus(mm)))*5.D-1
       end if
    end do

    deallocate(Gamma0)
    deallocate(Indof3rdPhonon)
    deallocate(Indof2ndPhonon)

    call MPI_ALLREDUCE(Indof2ndPhonon_plus_reduce,Indof2ndPhonon_plus,&
         Ntotal_plus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Indof3rdPhonon_plus_reduce,Indof3rdPhonon_plus,&
         Ntotal_plus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Indof2ndPhonon_minus_reduce,Indof2ndPhonon_minus,&
         Ntotal_minus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Indof3rdPhonon_minus_reduce,Indof3rdPhonon_minus,&
         Ntotal_minus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Gamma_plus_reduce,Gamma_plus,Ntotal_plus,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(Gamma_minus_reduce,Gamma_minus,Ntotal_minus,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(rate_scatt_plus_reduce,rate_scatt_plus,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(rate_scatt_minus_reduce,rate_scatt_minus,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(WP3_plus_reduce,WP3_plus,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ll)
    call MPI_ALLREDUCE(WP3_minus_reduce,WP3_minus,Nbands*Nlist,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ll)
    rate_scatt=rate_scatt_plus+rate_scatt_minus

    deallocate(Indof2ndPhonon_plus_reduce)
    deallocate(Indof3rdPhonon_plus_reduce)
    deallocate(Indof2ndPhonon_minus_reduce)
    deallocate(Indof3rdPhonon_minus_reduce)
    deallocate(Gamma_plus_reduce)
    deallocate(Gamma_minus_reduce)
  end subroutine Ind_driver

  ! Compute the number of allowed absorption processes and their contribution
  ! to phase space.
  subroutine NP_plus(mm,energy,velocity,Nlist,List,IJK,N_plus,P_plus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    integer(kind=4),intent(out) :: N_plus
    real(kind=8),intent(out) :: P_plus

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss
    real(kind=8) :: sigma
    real(kind=8) :: omega,omegap,omegadp
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3)

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_plus=0
    P_plus=0.d00
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(list(ll),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             omegap=energy(ii,j)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
                omegadp=energy(ss,k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(ii,j,:)-&
                        velocity(ss,k,:))
                   if(abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                   !-------- call ws_cell function to distinguish normal and
                   !umklapp processes
                   call ws_cell(q,shortest_q)
                   call ws_cell(qprime,shortest_qprime)
                   call ws_cell(qdprime,shortest_qdprime)
                   !----------------------Normal process condition
                   if (normal) then
                   if (abs(shortest_q(1)+shortest_qprime(1)-shortest_qdprime(1))<1e-10.and.&
                       abs(shortest_q(2)+shortest_qprime(2)-shortest_qdprime(2))<1e-10.and.&
                       abs(shortest_q(3)+shortest_qprime(3)-shortest_qdprime(3))<1e-10) then
                   !-----------------------------------------------------------
                      N_plus=N_plus+1
                      P_plus=P_plus+&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                   !---------------------Umklapp process condition
                   else if (umklapp) then
                   if (abs(shortest_q(1)+shortest_qprime(1)-shortest_qdprime(1))>1e-10.or.&
                       abs(shortest_q(2)+shortest_qprime(2)-shortest_qdprime(2))>1e-10.or.&
                       abs(shortest_q(3)+shortest_qprime(3)-shortest_qdprime(3))>1e-10) then
                   !-----------------------------------------------------------
                      N_plus=N_plus+1
                      P_plus=P_plus+&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                   !-------------------total scattering process
                   else
                   !-----------------------------------------------------------
                      N_plus=N_plus+1
                      P_plus=P_plus+&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                   end if
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
    end if
  end subroutine NP_plus

  ! Same as NP_plus, but for emission processes.
  subroutine NP_minus(mm,energy,velocity,Nlist,List,IJK,N_minus,P_minus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    integer(kind=4),intent(out) :: N_minus
    real(kind=8),intent(out) :: P_minus

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss
    real(kind=8) :: sigma
    real(kind=8) :: omega,omegap,omegadp
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3)

    do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus=0
    P_minus=0.d00
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(list(ll),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             omegap=energy(ii,j)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                qdprime=q-qprime
                qdprime=modulo(qdprime,Ngrid)
                ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
                omegadp=energy(ss,k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(ii,j,:)-&
                        velocity(ss,k,:))
                   if(abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                   !-------- call ws_cell function to distinguish normal and
                   !umklapp processes
                   call ws_cell(q,shortest_q)
                   call ws_cell(qprime,shortest_qprime)
                   call ws_cell(qdprime,shortest_qdprime)
                   !----------------------Normal process condition
                   if (normal) then
                   if (abs(shortest_q(1)-shortest_qprime(1)-shortest_qdprime(1))<1e-10.and.&
                       abs(shortest_q(2)-shortest_qprime(2)-shortest_qdprime(2))<1e-10.and.&
                       abs(shortest_q(3)-shortest_qprime(3)-shortest_qdprime(3))<1e-10) then
                    !-----------------------------------------------------------------------------
                      N_minus=N_minus+1
                      P_minus=P_minus+&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                   !--------------------Umklapp process condition
                   else if (umklapp) then
                   if (abs(shortest_q(1)-shortest_qprime(1)-shortest_qdprime(1))>1e-10.or.&
                       abs(shortest_q(2)-shortest_qprime(2)-shortest_qdprime(2))>1e-10.or.&
                       abs(shortest_q(3)-shortest_qprime(3)-shortest_qdprime(3))>1e-10) then
                    !-----------------------------------------------------------------------------
                      N_minus=N_minus+1
                      P_minus=P_minus+&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                   !--------------------total scattering process
                   else
                    !-----------------------------------------------------------------------------
                      N_minus=N_minus+1
                      P_minus=P_minus+&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/&
                           (sigma*sqrt(Pi)*nptk**2*nbands**3)
                   end if
                   end if
                end if
             end do ! k
             !--------END emission process-------------!
          end do ! ii
       end do  ! j
    end if
  end subroutine NP_minus

  ! Wrapper around NP_plus and NP_minus that splits the work among processors.
  subroutine NP_driver(energy,velocity,Nlist,List,IJK,&
       N_plus,Pspace_plus_total,N_minus,Pspace_minus_total)
    implicit none

    include "mpif.h"

    real(kind=8),intent(in) :: energy(nptk,nbands)
    real(kind=8),intent(in) :: velocity(nptk,nbands,3)
    integer(kind=4),intent(in) :: NList
    integer(kind=4),intent(in) :: List(Nlist)
    integer(kind=4),intent(in) :: IJK(3,nptk)
    integer(kind=4),intent(out) :: N_plus(Nlist*Nbands)
    integer(kind=4),intent(out) :: N_minus(Nlist*Nbands)
    real(kind=8),intent(out) :: Pspace_plus_total(Nbands,Nlist)
    real(kind=8),intent(out) :: Pspace_minus_total(Nbands,Nlist)

    integer(kind=4) :: mm
    integer(kind=4) :: N_plus_reduce(Nlist*Nbands)
    integer(kind=4) :: N_minus_reduce(Nlist*Nbands)
    real(kind=8) :: Pspace_plus_reduce(Nlist*Nbands)
    real(kind=8) :: Pspace_minus_reduce(Nlist*Nbands)

    Pspace_plus_total=0.d0
    Pspace_plus_reduce=0.d0
    Pspace_minus_total=0.d0
    Pspace_minus_reduce=0.d0
    N_plus=0
    N_minus=0
    N_plus_reduce=0
    N_minus_reduce=0

    do mm=myid+1,Nbands*Nlist,numprocs
       if (energy(List(int((mm-1)/Nbands)+1),modulo(mm-1,Nbands)+1).le.omega_max) then
          call NP_plus(mm,energy,velocity,Nlist,List,IJK,&
               N_plus_reduce(mm),Pspace_plus_reduce(mm))
          call NP_minus(mm,energy,velocity,Nlist,List,IJK,&
               N_minus_reduce(mm),Pspace_minus_reduce(mm))
       endif
    end do

    call MPI_ALLREDUCE(N_plus_reduce,N_plus,Nbands*Nlist,MPI_INTEGER,&
         MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(N_minus_reduce,N_minus,Nbands*Nlist,MPI_INTEGER,&
         MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(Pspace_plus_reduce,Pspace_plus_total,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(Pspace_minus_reduce,Pspace_minus_total,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
  end subroutine NP_driver

  ! RTA-only version of Ind_plus.
  subroutine RTA_plus(mm,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
       Gamma_plus,WP3_plus,Gamma_plus_N,Gamma_plus_U)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    real(kind=8),intent(out) :: Gamma_plus,WP3_plus
    real(kind=8),intent(out) :: Gamma_plus_N,Gamma_plus_U
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime
    real(kind=8) :: omega,omegap,omegadp
    real(kind=8) :: realq(3),realqprime(3),realqdprime(3)
    real(kind=8) :: WP3
    complex(kind=8) :: Vp

    Gamma_plus=0.d00
    Gamma_plus_N=0.d00
    Gamma_plus_U=0.d00

    WP3_plus=0.d00
  
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    realq=matmul(rlattvec,q/dble(ngrid))
    omega=energy(list(ll),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             realqprime=matmul(rlattvec,qprime/dble(ngrid))
             omegap=energy(ii,j)
             fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
             !--------BEGIN absorption process-----------
             do k=1,Nbands
                qdprime=q+qprime
                qdprime=modulo(qdprime,Ngrid)
                realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
                ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
                omegadp=energy(ss,k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(ii,j,:)-&
                        velocity(ss,k,:))
                   if(abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      WP3=(fBEprime-fBEdprime)*&
                           exp(-(omega+omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_plus=WP3_plus+WP3
                      if (.not.onlyharmonic) then
                      Vp=Vp_plus(i,j,k,list(ll),ii,ss,&
                           realqprime,realqdprime,eigenvect,&
                           Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                      Gamma_plus=Gamma_plus+hbarp*pi/4.d0*WP3*abs(Vp)**2
                      !======================Normal and Umklapp scattering===========================
                      !-------- call ws_cell function to distinguish normal and umklapp processes
                      call ws_cell(q,shortest_q)
                      call ws_cell(qprime,shortest_qprime)
                      call ws_cell(qdprime,shortest_qdprime)
                      !----------------------Normal process condition
                      if (abs(shortest_q(1)+shortest_qprime(1)-shortest_qdprime(1))<1e-10.and.&
                          abs(shortest_q(2)+shortest_qprime(2)-shortest_qdprime(2))<1e-10.and.&
                          abs(shortest_q(3)+shortest_qprime(3)-shortest_qdprime(3))<1e-10) then
                         Gamma_plus_N=Gamma_plus_N+hbarp*pi/4.d0*WP3*abs(Vp)**2
                      else
                         Gamma_plus_U=Gamma_plus_U+hbarp*pi/4.d0*WP3*abs(Vp)**2
                      endif
                      !==============================================================================================
                      endif
                   end if
                end if
             end do ! k
             !--------END absorption process-------------!
          end do ! ii
       end do  ! j
       WP3_plus=WP3_plus/nptk
    end if
    Gamma_plus=Gamma_plus*5.60626442*1.d8/nptk ! THz
    Gamma_plus_N=Gamma_plus_N*5.60626442*1.d8/nptk ! THz
    Gamma_plus_U=Gamma_plus_U*5.60626442*1.d8/nptk ! THz
  end subroutine RTA_plus

  ! RTA-only version of Ind_minus.
  subroutine RTA_minus(mm,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
       Gamma_minus,WP3_minus,Gamma_minus_N,Gamma_minus_U)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    real(kind=8),intent(out) :: Gamma_minus,WP3_minus
    real(kind=8),intent(out) :: Gamma_minus_N,Gamma_minus_U
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_minus_count
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime
    real(kind=8) ::  omega,omegap,omegadp
    real(kind=8) :: realqprime(3),realqdprime(3)
    real(kind=8) :: WP3
    complex(kind=8) :: Vp

    Gamma_minus=0.d00
    Gamma_minus_N=0.d00
    Gamma_minus_U=0.d00
    WP3_minus=0.d00

    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    N_minus_count=0
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(list(ll),i)
    if(omega.ne.0) then
       do j=1,Nbands
          do ii=1,nptk
             qprime=IJK(:,ii)
             realqprime=matmul(rlattvec,qprime/dble(ngrid))
             omegap=energy(ii,j)
             fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
             !--------BEGIN emission process-----------
             do k=1,Nbands
                qdprime=q-qprime
                qdprime=modulo(qdprime,Ngrid)
                realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
                ss=Index_N(qdprime(1),qdprime(2),qdprime(3))
                omegadp=energy(ss,k)
                if ((omegap.ne.0).and.(omegadp.ne.0)) then
                   sigma=scalebroad*base_sigma(&
                        velocity(ii,j,:)-&
                        velocity(ss,k,:))
                   if (abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                      fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                      WP3=(fBEprime+fBEdprime+1)*&
                           exp(-(omega-omegap-omegadp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                           (omega*omegap*omegadp)
                      WP3_minus=WP3_minus+WP3
                      if (.not.onlyharmonic) then
                      Vp=Vp_minus(i,j,k,list(ll),ii,ss,&
                           realqprime,realqdprime,eigenvect,&
                           Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
                      Gamma_minus=Gamma_minus+hbarp*pi/4.d0*WP3*abs(Vp)**2
                     !---------------------------Normal and Umklapp scattering rates---------------------------
                     !-------- call ws_cell function to distinguish normal and umklapp processes
                     call ws_cell(q,shortest_q)
                     call ws_cell(qprime,shortest_qprime)
                     call ws_cell(qdprime,shortest_qdprime)
                     !----------------------Normal process condition----------------------
                      if (abs(shortest_q(1)-shortest_qprime(1)-shortest_qdprime(1))<1e-10.and.&
                          abs(shortest_q(2)-shortest_qprime(2)-shortest_qdprime(2))<1e-10.and.&
                          abs(shortest_q(3)-shortest_qprime(3)-shortest_qdprime(3))<1e-10) then
                         Gamma_minus_N=Gamma_minus_N+hbarp*pi/4.d0*WP3*abs(Vp)**2
                      else
                         Gamma_minus_U=Gamma_minus_U+hbarp*pi/4.d0*WP3*abs(Vp)**2
                      endif
                     !==================================================                     
                      endif
                   end if
                end if
             end do ! k
             !--------END emission process-------------
          end do ! ii
       end do  ! j
       WP3_minus=WP3_minus*5.d-1/nptk

    end if
    Gamma_minus=Gamma_minus*5.60626442*1.d8/nptk
    Gamma_minus_N=Gamma_minus_N*5.60626442*1.d8/nptk
    Gamma_minus_U=Gamma_minus_U*5.60626442*1.d8/nptk

  end subroutine RTA_minus

  ! Wrapper around 3ph RTA_plus and RTA_minus that splits the work among processors.
  subroutine RTA_driver(energy,velocity,eigenvect,Nlist,List,IJK,&
       Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,rate_scatt,rate_scatt_plus,rate_scatt_minus,WP3_plus,WP3_minus,&
       rate_scatt_plus_N,rate_scatt_minus_N,rate_scatt_plus_U,rate_scatt_minus_U)
    implicit none

    include "mpif.h"

    real(kind=8),intent(in) :: energy(nptk,nbands)
    real(kind=8),intent(in) :: velocity(nptk,nbands,3)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(in) :: NList
    integer(kind=4),intent(in) :: List(Nlist)
    integer(kind=4),intent(in) :: IJK(3,nptk)
    integer(kind=4),intent(in) :: Ntri
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri)
    real(kind=8),intent(in) :: R_j(3,Ntri)
    real(kind=8),intent(in) :: R_k(3,Ntri)
    integer(kind=4),intent(in) :: Index_i(Ntri)
    integer(kind=4),intent(in) :: Index_j(Ntri)
    integer(kind=4),intent(in) :: Index_k(Ntri)
    real(kind=8),intent(out) :: rate_scatt(Nbands,Nlist),rate_scatt_plus(Nbands,Nlist),rate_scatt_minus(Nbands,Nlist)
    real(kind=8),intent(out) :: rate_scatt_plus_N(Nbands,Nlist),rate_scatt_minus_N(Nbands,Nlist),rate_scatt_plus_U(Nbands,Nlist),rate_scatt_minus_U(Nbands,Nlist)
    real(kind=8),intent(out) :: WP3_plus(Nbands,Nlist)
    real(kind=8),intent(out) :: WP3_minus(Nbands,Nlist)

    integer(kind=4) :: i
    integer(kind=4) :: ll
    integer(kind=4) :: mm
    real(kind=8) :: Gamma_plus,Gamma_minus
    real(kind=8) :: Gamma_plus_N,Gamma_minus_N,Gamma_plus_U,Gamma_minus_U
    real(kind=8) :: rate_scatt_plus_reduce(Nbands,Nlist),rate_scatt_minus_reduce(Nbands,Nlist)
    real(kind=8) :: rate_scatt_plus_reduce_N(Nbands,Nlist),rate_scatt_minus_reduce_N(Nbands,Nlist)
    real(kind=8) :: rate_scatt_plus_reduce_U(Nbands,Nlist),rate_scatt_minus_reduce_U(Nbands,Nlist) 
    real(kind=8) :: WP3_plus_reduce(Nbands*Nlist)
    real(kind=8) :: WP3_minus_reduce(Nbands*Nlist)


    rate_scatt=0.d00
    rate_scatt_plus_reduce=0.d00
    rate_scatt_plus_reduce_N=0.d00
    rate_scatt_plus_reduce_U=0.d00
    rate_scatt_minus_reduce=0.d00
    rate_scatt_minus_reduce_N=0.d00
    rate_scatt_minus_reduce_U=0.d00
    WP3_plus=0.d00
    WP3_minus=0.d00
    WP3_plus_reduce=0.d00
    WP3_minus_reduce=0.d00


    do mm=myid+1,Nbands*NList,numprocs
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       if (energy(List(ll),i).le.omega_max) then
          call RTA_plus(mm,energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
               Gamma_plus,WP3_plus_reduce(mm),&
               Gamma_plus_N,Gamma_plus_U)
          rate_scatt_plus_reduce(i,ll)=Gamma_plus
          rate_scatt_plus_reduce_N(i,ll)=Gamma_plus_N
          rate_scatt_plus_reduce_U(i,ll)=Gamma_plus_U
          call RTA_minus(mm,energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
               Gamma_minus,WP3_minus_reduce(mm),Gamma_minus_N,Gamma_minus_U)
          rate_scatt_minus_reduce(i,ll)=Gamma_minus*5.D-1
          rate_scatt_minus_reduce_N(i,ll)=Gamma_minus_N*5.D-1
          rate_scatt_minus_reduce_U(i,ll)=Gamma_minus_U*5.D-1
       endif
    end do

    call MPI_ALLREDUCE(rate_scatt_plus_reduce,rate_scatt_plus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_plus_reduce_N,rate_scatt_plus_N,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_plus_reduce_U,rate_scatt_plus_U,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_minus_reduce,rate_scatt_minus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_minus_reduce_N,rate_scatt_minus_N,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_minus_reduce_U,rate_scatt_minus_U,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(WP3_plus_reduce,WP3_plus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(WP3_minus_reduce,WP3_minus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

    rate_scatt=rate_scatt_plus+rate_scatt_minus
  end subroutine RTA_driver

  ! Subroutines for 4ph calculations

  ! Compute the number of allowed four-phonon ++ +- -- processes and
  ! their contribution to phase space.
  subroutine NP_plusplus(mm,energy,velocity,Nlist,List,IJK,N_plusplus,P_plusplus)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    integer(kind=8),intent(out) :: N_plusplus
    real(kind=8),intent(out) :: P_plusplus

    integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
    real(kind=8) :: sigma
    real(kind=8) :: omega,omegap,omegadp,omegatp

    do ii=0,Ngrid(1)-1       ! G1 direction
      do jj=0,Ngrid(2)-1     ! G2 direction
         do kk=0,Ngrid(3)-1  ! G3 direction
            Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
         end do
      end do
    end do
    N_plusplus=0
    P_plusplus=0.d00
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    omega=energy(list(ll),i)
    if (omega.ne.0) then
      do ii=1,nptk
         do jj=ii,nptk
            do ss=1,nptk
               qprime=IJK(:,ii)
               qdprime=IJK(:,jj)
               qtprime=IJK(:,ss)
               if (all(qtprime.eq.modulo(q+qprime+qdprime,ngrid))) then
                  do j=1,Nbands
                     startbranch=1
                     if(jj==ii) startbranch=j
                     do k=startbranch,Nbands
                        do l=1,Nbands
                           omegap=energy(ii,j)
                           omegadp=energy(jj,k)
                           omegatp=energy(ss,l)
                           if ((omegap.ne.0).and.(omegadp.ne.0).and.(omegatp.ne.0)) then
                              sigma=scalebroad*base_sigma(&
                              -velocity(jj,k,:)+velocity(ss,l,:))
                              if (abs(omega+omegap+omegadp-omegatp).le.sigma) then
                                 if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                 N_plusplus=N_plusplus+1
                                 P_plusplus=P_plusplus+&
                                       exp(-(omega+omegap+omegadp-omegatp)**2/(sigma**2))/&
                                       (sigma*sqrt(Pi)*dble(nptk)**3*nbands**4)
                              end if
                           end if
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
    end if
  end subroutine

  subroutine NP_plusminus(mm,energy,velocity,Nlist,List,IJK,N_plusminus,P_plusminus)
      implicit none

      integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
      real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
      integer(kind=8),intent(out) :: N_plusminus
      real(kind=8),intent(out) :: P_plusminus

      integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
      integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
      integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
      real(kind=8) :: sigma
      real(kind=8) :: omega,omegap,omegadp,omegatp

      do ii=0,Ngrid(1)-1     ! G1 direction
      do jj=0,Ngrid(2)-1     ! G2 direction
         do kk=0,Ngrid(3)-1  ! G3 direction
            Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
         end do
      end do
      end do
      N_plusminus=0
      P_plusminus=0.d00
      i=modulo(mm-1,Nbands)+1
      ll=int((mm-1)/Nbands)+1
      q=IJK(:,list(ll))
      omega=energy(list(ll),i)
      if (omega.ne.0) then
         do ii=1,nptk
            do jj=1,nptk
               do ss=jj,nptk
                  qprime=IJK(:,ii)
                  qdprime=IJK(:,jj)
                  qtprime=IJK(:,ss)
                  if (all(qtprime.eq.modulo(q+qprime-qdprime,ngrid))) then
                     do j=1,Nbands
                        do k=1,Nbands
                           startbranch=1
                           if(ss==jj) startbranch=k
                           do l=startbranch,Nbands
                              omegap=energy(ii,j)
                              omegadp=energy(jj,k)
                              omegatp=energy(ss,l)
                              if ((omegap.ne.0).and.(omegadp.ne.0).and.(omegatp.ne.0)) then
                                 sigma=scalebroad*base_sigma(&
                                 velocity(jj,k,:)-velocity(ss,l,:))
                                 if ((list(ll).ne.jj .or. i.ne.k).and.(list(ll).ne.ss .or. i.ne.l).and.&
                                    (ii.ne.jj .or. j.ne.k).and.(ii.ne.ss .or. j.ne.l) )then ! none of the two pairs can be the same
                                    if (abs(omega+omegap-omegadp-omegatp).le.sigma) then
                                       if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                       N_plusminus=N_plusminus+1
                                       P_plusminus=P_plusminus+&
                                             exp(-(omega+omegap-omegadp-omegatp)**2/(sigma**2))/&
                                             (sigma*sqrt(Pi)*dble(nptk)**3*nbands**4) 
                                    end if
                                 end if
                              end if
                           end do
                        end do
                     end do
                  end if
               end do
            end do
         end do
       end if
  end subroutine

  subroutine NP_minusminus(mm,energy,velocity,Nlist,List,IJK,N_minusminus,P_minusminus)
      implicit none

      integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
      real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
      integer(kind=8),intent(out) :: N_minusminus
      real(kind=8),intent(out) :: P_minusminus

      integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
      integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
      integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
      real(kind=8) :: sigma
      real(kind=8) :: omega,omegap,omegadp,omegatp

      do ii=0,Ngrid(1)-1        ! G1 direction
      do jj=0,Ngrid(2)-1     ! G2 direction
         do kk=0,Ngrid(3)-1  ! G3 direction
            Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
         end do
      end do
      end do
      N_minusminus=0
      P_minusminus=0.d00
      i=modulo(mm-1,Nbands)+1
      ll=int((mm-1)/Nbands)+1
      q=IJK(:,list(ll))
      omega=energy(list(ll),i)
      if (omega.ne.0) then
         do ii=1,nptk
            do jj=ii,nptk
               do ss=jj,nptk
                  qprime=IJK(:,ii)
                  qdprime=IJK(:,jj)
                  qtprime=IJK(:,ss)
                  if (all(qtprime.eq.modulo(q-qprime-qdprime,ngrid))) then
                     do j=1,Nbands
                        startbranch=1
                        if(jj==ii) startbranch=j
                        do k=startbranch,Nbands
                           startbranch=1
                           if(ss==jj) startbranch=k
                           do l=startbranch,Nbands
                              omegap=energy(ii,j)
                              omegadp=energy(jj,k)
                              omegatp=energy(ss,l)
                              if ((omegap.ne.0).and.(omegadp.ne.0).and.(omegatp.ne.0)) then
                                 sigma=scalebroad*base_sigma(&
                                 velocity(jj,k,:)-velocity(ss,l,:))
                                 if (abs(omega-omegap-omegadp-omegatp).le.sigma) then
                                    if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                    N_minusminus=N_minusminus+1
                                    P_minusminus=P_minusminus+&
                                          exp(-(omega-omegap-omegadp-omegatp)**2/(sigma**2))/&
                                          (sigma*sqrt(Pi)*dble(nptk)**3*nbands**4)
                                 end if
                              end if
                           end do
                        end do
                     end do
                  end if
               end do
            end do
         end do
       end if
  end subroutine

  ! Wrapper around NP_plusplus,NP_plusminus and NP_minusminus that splits the work among processors.
  subroutine NP_driver_4ph(energy,velocity,Nlist,List,IJK,&
         N_plusplus,Pspace_plusplus_total,N_plusminus,Pspace_plusminus_total,N_minusminus,Pspace_minusminus_total)
    implicit none

    include "mpif.h"

    real(kind=8),intent(in) :: energy(nptk,nbands)
    real(kind=8),intent(in) :: velocity(nptk,nbands,3)
    integer(kind=4),intent(in) :: NList
    integer(kind=4),intent(in) :: List(Nlist)
    integer(kind=4),intent(in) :: IJK(3,nptk)
    integer(kind=8),intent(out) :: N_plusplus(Nlist*Nbands)
    integer(kind=8),intent(out) :: N_plusminus(Nlist*Nbands)
    integer(kind=8),intent(out) :: N_minusminus(Nlist*Nbands)
    real(kind=8),intent(out) :: Pspace_plusplus_total(Nbands,Nlist)
    real(kind=8),intent(out) :: Pspace_plusminus_total(Nbands,Nlist)
    real(kind=8),intent(out) :: Pspace_minusminus_total(Nbands,Nlist)

    integer(kind=4) :: mm
    integer(kind=8) :: N_plusplus_reduce(Nlist*Nbands)
    integer(kind=8) :: N_plusminus_reduce(Nlist*Nbands)
    integer(kind=8) :: N_minusminus_reduce(Nlist*Nbands)
    real(kind=8) :: Pspace_plusplus_reduce(Nlist*Nbands)
    real(kind=8) :: Pspace_plusminus_reduce(Nlist*Nbands)
    real(kind=8) :: Pspace_minusminus_reduce(Nlist*Nbands)

    Pspace_plusplus_total=0.d0
    Pspace_plusplus_reduce=0.d0
    Pspace_plusminus_total=0.d0
    Pspace_plusminus_reduce=0.d0
    Pspace_minusminus_total=0.d0
    Pspace_minusminus_reduce=0.d0
    N_plusplus=0
    N_plusminus=0
    N_minusminus=0
    N_plusplus_reduce=0
    N_plusminus_reduce=0
    N_minusminus_reduce=0

    do mm=myid+1,Nbands*Nlist,numprocs
      if (energy(List(int((mm-1)/Nbands)+1),modulo(mm-1,Nbands)+1).le.omega_max) then
         call NP_plusplus(mm,energy,velocity,Nlist,List,IJK,N_plusplus_reduce(mm),Pspace_plusplus_reduce(mm))
         call NP_plusminus(mm,energy,velocity,Nlist,List,IJK,N_plusminus_reduce(mm),Pspace_plusminus_reduce(mm))
         call NP_minusminus(mm,energy,velocity,Nlist,List,IJK,N_minusminus_reduce(mm),Pspace_minusminus_reduce(mm))
      endif
    end do

    call MPI_ALLREDUCE(N_plusplus_reduce,N_plusplus,Nbands*Nlist,MPI_INTEGER8,&
         MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(N_plusminus_reduce,N_plusminus,Nbands*Nlist,MPI_INTEGER8,&
         MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(N_minusminus_reduce,N_minusminus,Nbands*Nlist,MPI_INTEGER8,&
         MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(Pspace_plusplus_reduce,Pspace_plusplus_total,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(Pspace_plusminus_reduce,Pspace_plusminus_total,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(Pspace_minusminus_reduce,Pspace_minusminus_total,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

  end subroutine NP_driver_4ph

  ! RTA-only version of 4ph ++ process, Ind_plusplus; Avoid double counting
  subroutine RTA_plusplus(mm,energy,velocity,eigenvect,Nlist,List,&
      Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
      Gamma_plusplus,WP4_plusplus,Gamma_plusplus_N,Gamma_plusplus_U)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    real(kind=8),intent(out) :: Gamma_plusplus,WP4_plusplus
    real(kind=8),intent(out) :: Gamma_plusplus_N,Gamma_plusplus_U
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3),shortest_qtprime(3)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime,fBEtprime
    real(kind=8) :: omega,omegap,omegadp,omegatp
    real(kind=8) :: realq(3),realqprime(3),realqdprime(3),realqtprime(3)
    real(kind=8) :: WP4
    complex(kind=8) :: Vp

    Gamma_plusplus=0.d00
    Gamma_plusplus_N=0.d00
    Gamma_plusplus_U=0.d00
    WP4_plusplus=0.d00
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
            Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    realq=matmul(rlattvec,q/dble(ngrid))
    omega=energy(list(ll),i)

    if (omega.ne.0) then
      do ii=1,nptk
         do jj=ii,nptk
            do ss=1,nptk
               qprime=IJK(:,ii)
               realqprime=matmul(rlattvec,qprime/dble(ngrid))
               qdprime=IJK(:,jj)
               realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
               qtprime=IJK(:,ss)
               realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
               if (all(qtprime.eq.modulo(q+qprime+qdprime,ngrid))) then
                  do j=1,Nbands
                     startbranch=1
                     if(jj==ii) startbranch=j
                     do k=startbranch,Nbands
                        do l=1,Nbands
                           omegap=energy(ii,j)
                           omegadp=energy(jj,k)
                           omegatp=energy(ss,l)
                           if ((omegap.ne.0).and.(omegadp.ne.0).and.(omegatp.ne.0)) then
                              sigma=scalebroad*base_sigma(&
                              -velocity(jj,k,:)+velocity(ss,l,:))
                              if (abs(omega+omegap+omegadp-omegatp).le.sigma) then
                                 if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                 fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
                                 fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                                 fBEtprime=1.d0/(exp(hbar*omegatp/Kb/T)-1.D0)
                                 WP4=(fBEprime*fBEdprime*(1+fBEtprime)-(1+fBEprime)*(1+fBEdprime)*fBEtprime)*&
                                    exp(-(omega+omegap+omegadp-omegatp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                                    (omega*omegap*omegadp*omegatp)
                                 if(omegap.le.1.25 .or. omegadp.le.1.25 .or.omegatp.le.1.25) WP4=0
                                 WP4_plusplus=WP4_plusplus+WP4
                                 if (.not.onlyharmonic) then
                                    Vp=Vp_pp(i,j,k,l,list(ll),ii,jj,ss,&
                                             realqprime,realqdprime,realqtprime,eigenvect,&
                                             Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
                                    Gamma_plusplus=Gamma_plusplus+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                    !========Normal and Umklapp scattering========
                                    call ws_cell(q,shortest_q)
                                    call ws_cell(qprime,shortest_qprime)
                                    call ws_cell(qdprime,shortest_qdprime)
                                    call ws_cell(qtprime,shortest_qtprime)
                                    if (abs(shortest_q(1)+shortest_qprime(1)+shortest_qdprime(1)-shortest_qtprime(1))<1e-10.and.&
                                       abs(shortest_q(2)+shortest_qprime(2)+shortest_qdprime(2)-shortest_qtprime(2))<1e-10.and.&
                                       abs(shortest_q(3)+shortest_qprime(3)+shortest_qdprime(3)-shortest_qtprime(3))<1e-10) then
                                       Gamma_plusplus_N=Gamma_plusplus_N+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                       else
                                       Gamma_plusplus_U=Gamma_plusplus_U+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                    end if
                                 end if
                              end if
                           end if
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
      WP4_plusplus=WP4_plusplus/nptk/nptk
    end if
    ! converting to THz
    Gamma_plusplus=Gamma_plusplus*3.37617087*1.d9/nptk/nptk
    Gamma_plusplus_N=Gamma_plusplus_N*3.37617087*1.d9/nptk/nptk
    Gamma_plusplus_U=Gamma_plusplus_U*3.37617087*1.d9/nptk/nptk

  end subroutine

  ! RTA-only version of 4ph +- process, Ind_plusminus
  subroutine RTA_plusminus(mm,energy,velocity,eigenvect,Nlist,List,&
      Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
      Gamma_plusminus,WP4_plusminus,Gamma_plusminus_N,Gamma_plusminus_U)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    real(kind=8),intent(out) :: Gamma_plusminus,WP4_plusminus
    real(kind=8),intent(out) :: Gamma_plusminus_N,Gamma_plusminus_U
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3),shortest_qtprime(3)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime,fBEtprime
    real(kind=8) :: omega,omegap,omegadp,omegatp
    real(kind=8) :: realq(3),realqprime(3),realqdprime(3),realqtprime(3)
    real(kind=8) :: WP4
    complex(kind=8) :: Vp

    Gamma_plusminus=0.d00
    Gamma_plusminus_N=0.d00
    Gamma_plusminus_U=0.d00
    WP4_plusminus=0.d00
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
            Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    realq=matmul(rlattvec,q/dble(ngrid))
    omega=energy(list(ll),i)

    if (omega.ne.0) then
      do ii=1,nptk
         do jj=1,nptk
            do ss=jj,nptk
               qprime=IJK(:,ii)
               realqprime=matmul(rlattvec,qprime/dble(ngrid))
               qdprime=IJK(:,jj)
               realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
               qtprime=IJK(:,ss)
               realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
               if (all(qtprime.eq.modulo(q+qprime-qdprime,ngrid))) then
                  do j=1,Nbands
                     do k=1,Nbands
                        startbranch=1
                        if(ss==jj) startbranch=k
                        do l=startbranch,Nbands
                           omegap=energy(ii,j)
                           omegadp=energy(jj,k)
                           omegatp=energy(ss,l)
                           if ((omegap.ne.0).and.(omegadp.ne.0).and.(omegatp.ne.0)) then
                              sigma=scalebroad*base_sigma(&
                              velocity(jj,k,:)-velocity(ss,l,:))
                              if ((list(ll).ne.jj .or. i.ne.k).and.(list(ll).ne.ss .or. i.ne.l).and.&
                                 (ii.ne.jj .or. j.ne.k).and.(ii.ne.ss .or. j.ne.l) )then ! none of the two pairs can be the same
                                 if (abs(omega+omegap-omegadp-omegatp).le.sigma) then
                                    if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                    fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
                                    fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                                    fBEtprime=1.d0/(exp(hbar*omegatp/Kb/T)-1.D0)
                                    WP4=(fBEprime*(1+fBEdprime)*(1+fBEtprime)-(1+fBEprime)*fBEdprime*fBEtprime)*&
                                       exp(-(omega+omegap-omegadp-omegatp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                                       (omega*omegap*omegadp*omegatp)
                                    if(omegap.le.1.25 .or. omegadp.le.1.25 .or.omegatp.le.1.25) WP4=0
                                    WP4_plusminus=WP4_plusminus+WP4
                                    if (.not.onlyharmonic) then
                                       Vp=Vp_pm(i,j,k,l,list(ll),ii,jj,ss,&
                                                realqprime,realqdprime,realqtprime,eigenvect,&
                                                Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
                                       Gamma_plusminus=Gamma_plusminus+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                       !========Normal and Umklapp scattering========
                                       call ws_cell(q,shortest_q)
                                       call ws_cell(qprime,shortest_qprime)
                                       call ws_cell(qdprime,shortest_qdprime)
                                       call ws_cell(qtprime,shortest_qtprime)
                                       if (abs(shortest_q(1)+shortest_qprime(1)-shortest_qdprime(1)-shortest_qtprime(1))<1e-10.and.&
                                          abs(shortest_q(2)+shortest_qprime(2)-shortest_qdprime(2)-shortest_qtprime(2))<1e-10.and.&
                                          abs(shortest_q(3)+shortest_qprime(3)-shortest_qdprime(3)-shortest_qtprime(3))<1e-10) then
                                          Gamma_plusminus_N=Gamma_plusminus_N+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                          else
                                          Gamma_plusminus_U=Gamma_plusminus_U+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                       end if
                                    end if
                                 end if
                              end if
                           end if
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
      WP4_plusminus=WP4_plusminus/nptk/nptk
    end if
    
    ! converting to THz
    Gamma_plusminus=Gamma_plusminus*3.37617087*1.d9/nptk/nptk
    Gamma_plusminus_N=Gamma_plusminus_N*3.37617087*1.d9/nptk/nptk
    Gamma_plusminus_U=Gamma_plusminus_U*3.37617087*1.d9/nptk/nptk

  end subroutine

  ! RTA-only version of 4ph -- process, Ind_minusminus
  subroutine RTA_minusminus(mm,energy,velocity,eigenvect,Nlist,List,&
      Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
      Gamma_minusminus,WP4_minusminus,Gamma_minusminus_N,Gamma_minusminus_U)
    implicit none

    integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
    integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
    real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
    real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    real(kind=8),intent(out) :: Gamma_minusminus,WP4_minusminus
    real(kind=8),intent(out) :: Gamma_minusminus_N,Gamma_minusminus_U
    real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3),shortest_qtprime(3)

    integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
    integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
    integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
    real(kind=8) :: sigma
    real(kind=8) :: fBEprime,fBEdprime,fBEtprime
    real(kind=8) :: omega,omegap,omegadp,omegatp
    real(kind=8) :: realq(3),realqprime(3),realqdprime(3),realqtprime(3)
    real(kind=8) :: WP4
    complex(kind=8) :: Vp

    Gamma_minusminus=0.d00
    Gamma_minusminus_N=0.d00
    Gamma_minusminus_U=0.d00
    WP4_minusminus=0.d00
    do ii=0,Ngrid(1)-1
       do jj=0,Ngrid(2)-1
          do kk=0,Ngrid(3)-1
            Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
    end do
    i=modulo(mm-1,Nbands)+1
    ll=int((mm-1)/Nbands)+1
    q=IJK(:,list(ll))
    realq=matmul(rlattvec,q/dble(ngrid))
    omega=energy(list(ll),i)

    if (omega.ne.0) then
      do ii=1,nptk
         do jj=ii,nptk
            do ss=jj,nptk
               qprime=IJK(:,ii)
               realqprime=matmul(rlattvec,qprime/dble(ngrid))
               qdprime=IJK(:,jj)
               realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
               qtprime=IJK(:,ss)
               realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
               if (all(qtprime.eq.modulo(q-qprime-qdprime,ngrid))) then
                  do j=1,Nbands
                     startbranch=1
                     if(jj==ii) startbranch=j
                     do k=startbranch,Nbands
                        startbranch=1
                        if(ss==jj) startbranch=k
                        do l=startbranch,Nbands
                           omegap=energy(ii,j)
                           omegadp=energy(jj,k)
                           omegatp=energy(ss,l)
                           if ((omegap.ne.0).and.(omegadp.ne.0).and.(omegatp.ne.0)) then
                              sigma=scalebroad*base_sigma(&
                              velocity(jj,k,:)-velocity(ss,l,:))
                              if (abs(omega-omegap-omegadp-omegatp).le.sigma) then
                                 if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                 fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
                                 fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                                 fBEtprime=1.d0/(exp(hbar*omegatp/Kb/T)-1.D0)
                                 WP4=((1+fBEprime)*(1+fBEdprime)*(1+fBEtprime)-fBEprime*fBEdprime*fBEtprime)*&
                                    exp(-(omega-omegap-omegadp-omegatp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                                    (omega*omegap*omegadp*omegatp)
                                 if(omegap.le.1.25 .or. omegadp.le.1.25 .or.omegatp.le.1.25) WP4=0
                                 WP4_minusminus=WP4_minusminus+WP4
                                 if (.not.onlyharmonic) then
                                    Vp=Vp_mm(i,j,k,l,list(ll),ii,jj,ss,&
                                             realqprime,realqdprime,realqtprime,eigenvect,&
                                             Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
                                    Gamma_minusminus=Gamma_minusminus+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                    !========Normal and Umklapp scattering========
                                    call ws_cell(q,shortest_q)
                                    call ws_cell(qprime,shortest_qprime)
                                    call ws_cell(qdprime,shortest_qdprime)
                                    call ws_cell(qtprime,shortest_qtprime)
                                    if (abs(shortest_q(1)-shortest_qprime(1)-shortest_qdprime(1)-shortest_qtprime(1))<1e-10.and.&
                                       abs(shortest_q(2)-shortest_qprime(2)-shortest_qdprime(2)-shortest_qtprime(2))<1e-10.and.&
                                       abs(shortest_q(3)-shortest_qprime(3)-shortest_qdprime(3)-shortest_qtprime(3))<1e-10) then
                                       Gamma_minusminus_N=Gamma_minusminus_N+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                       else
                                       Gamma_minusminus_U=Gamma_minusminus_U+hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                    end if
                                 end if
                              end if
                           end if
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
      WP4_minusminus=WP4_minusminus/nptk/nptk
    end if
    
    ! converting to THz
    Gamma_minusminus=Gamma_minusminus*3.37617087*1.d9/nptk/nptk
    Gamma_minusminus_N=Gamma_minusminus_N*3.37617087*1.d9/nptk/nptk
    Gamma_minusminus_U=Gamma_minusminus_U*3.37617087*1.d9/nptk/nptk

  end subroutine

  ! Wrapper around 4ph RTA subroutines with 3ph subroutines that splits the work among processors.
  subroutine RTA_driver_4ph(energy,velocity,eigenvect,Nlist,List,IJK,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,rate_scatt_4ph,rate_scatt_plusplus,rate_scatt_plusminus,rate_scatt_minusminus,&
       WP4_plusplus,WP4_plusminus,WP4_minusminus,&
       rate_scatt_plusplus_N,rate_scatt_plusminus_N,rate_scatt_minusminus_N,&
       rate_scatt_plusplus_U,rate_scatt_plusminus_U,rate_scatt_minusminus_U)
    implicit none

    include "mpif.h"

    real(kind=8),intent(in) :: energy(nptk,nbands)
    real(kind=8),intent(in) :: velocity(nptk,nbands,3)
    complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
    integer(kind=4),intent(in) :: NList
    integer(kind=4),intent(in) :: List(Nlist)
    integer(kind=4),intent(in) :: IJK(3,nptk)
    integer(kind=4),intent(in) :: Ntri
    real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri)
    real(kind=8),intent(in) :: R_j(3,Ntri)
    real(kind=8),intent(in) :: R_k(3,Ntri)
    real(kind=8),intent(in) :: R_l(3,Ntri)
    integer(kind=4),intent(in) :: Index_i(Ntri)
    integer(kind=4),intent(in) :: Index_j(Ntri)
    integer(kind=4),intent(in) :: Index_k(Ntri)
    integer(kind=4),intent(in) :: Index_l(Ntri)
    real(kind=8),intent(out) :: rate_scatt_4ph(Nbands,Nlist)
    real(kind=8),intent(out) :: rate_scatt_plusplus(Nbands,Nlist),rate_scatt_plusminus(Nbands,Nlist),rate_scatt_minusminus(Nbands,Nlist)
    real(kind=8),intent(out) :: rate_scatt_plusplus_N(Nbands,Nlist),rate_scatt_plusminus_N(Nbands,Nlist),rate_scatt_minusminus_N(Nbands,Nlist)
    real(kind=8),intent(out) :: rate_scatt_plusplus_U(Nbands,Nlist),rate_scatt_plusminus_U(Nbands,Nlist),rate_scatt_minusminus_U(Nbands,Nlist)
    real(kind=8),intent(out) :: WP4_plusplus(Nbands,Nlist),WP4_plusminus(Nbands,Nlist),WP4_minusminus(Nbands,Nlist)

    integer(kind=4) :: i
    integer(kind=4) :: ll
    integer(kind=4) :: mm
    real(kind=8) :: Gamma_plusplus,Gamma_plusminus,Gamma_minusminus
    real(kind=8) :: Gamma_plusplus_N,Gamma_plusminus_N,Gamma_minusminus_N,Gamma_plusplus_U,Gamma_plusminus_U,Gamma_minusminus_U
    real(kind=8) :: rate_scatt_plusplus_reduce(Nbands,Nlist),rate_scatt_plusminus_reduce(Nbands,Nlist),rate_scatt_minusminus_reduce(Nbands,Nlist)
    real(kind=8) :: rate_scatt_plusplus_reduce_N(Nbands,Nlist),rate_scatt_plusminus_reduce_N(Nbands,Nlist),rate_scatt_minusminus_reduce_N(Nbands,Nlist)
    real(kind=8) :: rate_scatt_plusplus_reduce_U(Nbands,Nlist),rate_scatt_plusminus_reduce_U(Nbands,Nlist),rate_scatt_minusminus_reduce_U(Nbands,Nlist) 
    real(kind=8) :: WP4_plusplus_reduce(Nbands*Nlist),WP4_plusminus_reduce(Nbands*Nlist),WP4_minusminus_reduce(Nbands*Nlist)

    rate_scatt_4ph=0.d00
    rate_scatt_plusplus_reduce=0.d00
    rate_scatt_plusplus_reduce_N=0.d00
    rate_scatt_plusplus_reduce_U=0.d00

    rate_scatt_plusminus_reduce=0.d00
    rate_scatt_plusminus_reduce_N=0.d00
    rate_scatt_plusminus_reduce_U=0.d00

    rate_scatt_minusminus_reduce=0.d00
    rate_scatt_minusminus_reduce_N=0.d00
    rate_scatt_minusminus_reduce_U=0.d00

    WP4_plusplus=0.d00
    WP4_plusminus=0.d00
    WP4_minusminus=0.d00

    WP4_plusplus_reduce=0.d00
    WP4_plusminus_reduce=0.d00
    WP4_minusminus_reduce=0.d00

    do mm=myid+1,Nbands*NList,numprocs
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       if (energy(List(ll),i).le.omega_max) then
          call RTA_plusplus(mm,energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
               Gamma_plusplus,WP4_plusplus_reduce(mm),Gamma_plusplus_N,Gamma_plusplus_U)
          rate_scatt_plusplus_reduce(i,ll)=Gamma_plusplus
          rate_scatt_plusplus_reduce_N(i,ll)=Gamma_plusplus_N
          rate_scatt_plusplus_reduce_U(i,ll)=Gamma_plusplus_U
          call RTA_plusminus(mm,energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
               Gamma_plusminus,WP4_plusminus_reduce(mm),Gamma_plusminus_N,Gamma_plusminus_U)
          rate_scatt_plusminus_reduce(i,ll)=Gamma_plusminus
          rate_scatt_plusminus_reduce_N(i,ll)=Gamma_plusminus_N
          rate_scatt_plusminus_reduce_U(i,ll)=Gamma_plusminus_U
          call RTA_minusminus(mm,energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
               Gamma_minusminus,WP4_minusminus_reduce(mm),Gamma_minusminus_N,Gamma_minusminus_U)
          rate_scatt_minusminus_reduce(i,ll)=Gamma_minusminus
          rate_scatt_minusminus_reduce_N(i,ll)=Gamma_minusminus_N
          rate_scatt_minusminus_reduce_U(i,ll)=Gamma_minusminus_U
       endif
    end do

    call MPI_ALLREDUCE(rate_scatt_plusplus_reduce,rate_scatt_plusplus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_plusplus_reduce_N,rate_scatt_plusplus_N,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_plusplus_reduce_U,rate_scatt_plusplus_U,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

    call MPI_ALLREDUCE(rate_scatt_plusminus_reduce,rate_scatt_plusminus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_plusminus_reduce_N,rate_scatt_plusminus_N,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_plusminus_reduce_U,rate_scatt_plusminus_U,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

    call MPI_ALLREDUCE(rate_scatt_minusminus_reduce,rate_scatt_minusminus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_minusminus_reduce_N,rate_scatt_minusminus_N,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(rate_scatt_minusminus_reduce_U,rate_scatt_minusminus_U,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

    call MPI_ALLREDUCE(WP4_plusplus_reduce,WP4_plusplus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(WP4_plusminus_reduce,WP4_plusminus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
    call MPI_ALLREDUCE(WP4_minusminus_reduce,WP4_minusminus,Nbands*Nlist,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

    rate_scatt_4ph=rate_scatt_plusplus+rate_scatt_plusminus+rate_scatt_minusminus
  end subroutine RTA_driver_4ph

end module processes
