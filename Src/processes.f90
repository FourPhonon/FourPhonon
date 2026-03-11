! Compute the number of allowed three-phonon / four-phonon processes, their
! scattering amplitudes and their phase-space volume.

#ifndef NUM_GANGS
#define NUM_GANGS 1024
#endif

#ifndef VEC_LEN
#define VEC_LEN 512
#endif

module processes
   use iso_fortran_env
   use misc
   use data
   use config
   use omp_lib
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

         ! prefactor = 1.d0 / sqrt(masses(types(Index_i(ll))) * masses(types(Index_j(ll))) * masses(types(Index_k(ll)))) * &
         !    cmplx(cos(dot_product(realqprime, R_j(:, ll))), sin(dot_product(realqprime, R_j(:, ll))), kind=8) * &
         !    cmplx(cos(-dot_product(realqdprime, R_k(:, ll))), sin(-dot_product(realqdprime, R_k(:, ll))), kind=8)

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

      ! prefactor = 1.d0 / sqrt(masses(types(Index_i(ll))) * masses(types(Index_j(ll))) * masses(types(Index_k(ll)))) * &
      !      cmplx(cos(-dot_product(realqprime, R_j(:, ll))), sin(-dot_product(realqprime, R_j(:, ll))), kind=8) * &
      !      cmplx(cos(-dot_product(realqdprime, R_k(:, ll))), sin(-dot_product(realqdprime, R_k(:, ll))), kind=8)

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
      ! data for output 
      integer(kind=4),intent(out) :: Indof2ndPhonon_plus(:)
      integer(kind=4),intent(out) :: Indof3rdPhonon_plus(:)
      real(kind=8),intent(out) :: Gamma_plus(:)
      integer(kind=4),intent(out) :: Indof2ndPhonon_minus(:)
      integer(kind=4),intent(out) :: Indof3rdPhonon_minus(:)
      real(kind=8),intent(out) :: Gamma_minus(:)
      real(kind=8),intent(out) :: rate_scatt(Nbands,Nlist),rate_scatt_plus(Nbands,Nlist),rate_scatt_minus(Nbands,Nlist)
      real(kind=8),intent(out) :: WP3_plus(Nbands,Nlist)
      real(kind=8),intent(out) :: WP3_minus(Nbands,Nlist)
      !  end data for output
   
      integer(kind=4) :: i
      integer(kind=4) :: ll
      integer(kind=4) :: mm
      integer(kind=4) :: maxsize
      integer(kind=4) :: Ntotal_plus
      integer(kind=4) :: Ntotal_minus
      integer(kind=4) :: Naccum_plus(Nbands*Nlist)
      integer(kind=4) :: Naccum_minus(Nbands*Nlist)
      !   temporary variables
      integer(kind=4),allocatable :: Indof2ndPhonon(:)
      integer(kind=4),allocatable :: Indof3rdPhonon(:)
      real(kind=8),allocatable :: Gamma0(:)
      real(kind=8) :: WP3_temp ! variable for parallel computing
      !   end temporary variables
      real(kind=8),allocatable :: WP3_plus_reduce(:) ! change to allocatable for reduce matrix
      real(kind=8),allocatable :: WP3_minus_reduce(:) 
      real(kind=8),allocatable :: rate_scatt_plus_reduce(:),rate_scatt_minus_reduce(:) 
      ! end reduced variables


      integer(kind=4) :: chunkstart, chunkend, chunkid ! for parallel computing

 
      maxsize=max(maxval(N_plus),maxval(N_minus))

      !   allocate temporary arrays
      allocate(Indof2ndPhonon(maxsize))
      allocate(Indof3rdPhonon(maxsize))
      allocate(Gamma0(maxsize))
      !   end allocate temporary arrays
   
      Naccum_plus(1)=0
      Naccum_minus(1)=0
      do mm=2,Nbands*Nlist
         Naccum_plus(mm)=Naccum_plus(mm-1)+N_plus(mm-1)
         Naccum_minus(mm)=Naccum_minus(mm-1)+N_minus(mm-1)
      end do
      Ntotal_plus=sum(N_plus)
      Ntotal_minus=sum(N_minus)
   
      !  allocate reduced output arrays
      allocate(rate_scatt_plus_reduce(Nbands*Nlist))
      allocate(rate_scatt_minus_reduce(Nbands*Nlist))
      allocate(WP3_plus_reduce(Nbands*Nlist))
      allocate(WP3_minus_reduce(Nbands*Nlist))
      !  end allocate reduced output arrays

      ! initialize the output arrays
      WP3_plus=0.d0
      WP3_minus=0.d0
      Indof2ndPhonon_plus=0
      Indof3rdPhonon_plus=0
      Indof2ndPhonon_minus=0
      Indof3rdPhonon_minus=0
      Gamma_plus=0.d0
      Gamma_minus=0.d0
      !   end initialize the output arrays

      ! initialize the reduced output arrays
      rate_scatt_plus_reduce=0.d0
      rate_scatt_minus_reduce=0.d0
      WP3_plus_reduce=0.d0
      WP3_minus_reduce=0.d0
      !   end initialize the reduced output arrays

      ! initialize for the temporary variables
      Indof2ndPhonon=0
      Indof3rdPhonon=0
      Gamma0=0.d0
      WP3_temp = 0.d0
      !   end initialize for the temporary variables
 

      do chunkid=myid+1,numchunk,numprocs
         chunkstart = (chunkid-1)*chunksize+1
         chunkend = chunkid*chunksize
         if (chunkend.gt.Nbands*Nlist) chunkend = Nbands*Nlist
         !$OMP PARALLEL DO default(shared) schedule(dynamic,1) private(mm,i,ll,Gamma0,Indof2ndPhonon,Indof3rdPhonon,WP3_temp)
         do mm=chunkstart,chunkend
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1
            if(N_plus(mm).ne.0) then
               call Ind_plus(mm,N_plus(mm),energy,velocity,eigenvect,Nlist,List,&
                     Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                     Indof2ndPhonon(1:N_plus(mm)),Indof3rdPhonon(1:N_plus(mm)),&
                     Gamma0(1:N_plus(mm)),WP3_temp)
                  WP3_plus_reduce(mm) = WP3_temp
                  Indof2ndPhonon_plus((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
                        Indof2ndPhonon(1:N_plus(mm))
                  Indof3rdPhonon_plus((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
                        Indof3rdPhonon(1:N_plus(mm))
                  Gamma_plus((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
                        Gamma0(1:N_plus(mm))
                  rate_scatt_plus_reduce(mm)=sum(Gamma0(1:N_plus(mm)))! change to mm since we allocate the array in 1d
            end if
            if(N_minus(mm).ne.0) then
               call Ind_minus(mm,N_minus(mm),energy,velocity,eigenvect,Nlist,List,&
                     Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                     Indof2ndPhonon(1:N_minus(mm)),Indof3rdPhonon(1:N_minus(mm)),&
                     Gamma0(1:N_minus(mm)),WP3_temp)
                  WP3_minus_reduce(mm) = WP3_temp
                  Indof2ndPhonon_minus((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
                        Indof2ndPhonon(1:N_minus(mm))
                  Indof3rdPhonon_minus((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
                        Indof3rdPhonon(1:N_minus(mm))
                  Gamma_minus((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
                        Gamma0(1:N_minus(mm))
                  rate_scatt_minus_reduce(mm)=sum(Gamma0(1:N_minus(mm)))*5.D-1
            end if
         end do
         !$OMP END PARALLEL DO
      end do

   
      call MPI_ALLREDUCE(MPI_IN_PLACE,Indof2ndPhonon_plus,&
            Ntotal_plus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Indof3rdPhonon_plus,&
            Ntotal_plus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Indof2ndPhonon_minus,&
            Ntotal_minus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Indof3rdPhonon_minus,&
            Ntotal_minus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Gamma_plus,Ntotal_plus,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ll)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Gamma_minus,Ntotal_minus,&
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

      !  deallocate temporary arrays
      deallocate(Gamma0)
      deallocate(Indof3rdPhonon)
      deallocate(Indof2ndPhonon)
      !  end deallocate temporary arrays
      ! deallocate reduced output arrays
      deallocate(rate_scatt_plus_reduce)
      deallocate(rate_scatt_minus_reduce)
      deallocate(WP3_plus_reduce)
      deallocate(WP3_minus_reduce)
      ! end deallocate reduced output arrays
   
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
         do ii=1,nptk
            do j=1,Nbands
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
                       N_plus=N_plus+1
                       P_plus=P_plus+&
                            exp(-(omega+omegap-omegadp)**2/(sigma**2))/&
                            (sigma*sqrt(Pi)*nptk**2*nbands**3)
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
         do ii=1,nptk
            do j=1,Nbands
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
                       N_minus=N_minus+1
                       P_minus=P_minus+&
                            exp(-(omega-omegap-omegadp)**2/(sigma**2))/&
                            (sigma*sqrt(Pi)*nptk**2*nbands**3)
                    end if
                 end if
              end do ! k
              !--------END emission process-------------!
           end do ! ii
        end do  ! j
     end if
   end subroutine NP_minus
 
 
   ! Compute the number of allowed absorption processes and their contribution
   ! to phase space using sampling method.
   subroutine NP_plus_sample(mm,energy,velocity,Nlist,List,IJK,N_plus,P_plus)
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
 
    ! ----------- sampling method add -----------
     integer(kind=4) :: nn, rand_num, iter
     real :: rand_matrix(num_sample_process_3ph_phase_space),rand
    ! ----------- end sampling method add -----------
 
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
       ! ----------- sampling method add -----------
       ! Note: 
       ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
       ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
       ! phonon3: ss is momentum [1:nptk], k is energy [1:Nband]
       call random_seed()
       call random_number(rand_matrix)
       do iter=1,num_sample_process_3ph_phase_space
          rand = rand_matrix(iter-1)
          rand_num = FLOOR(INT(Nbands*nptk*Nbands)*rand)+1  ! generate random integer in [1,Nbands*nptk*Nbands]
          nn = int((rand_num-1)/Nbands)+1 ! get a random phonon number for 2nd phonon
          k = modulo(rand_num-1,Nbands)+1 ! phonon3 branch [1:Nband]
          j=modulo(nn-1,Nbands)+1 ! phonon2 branch [1:Nband]
          ii=int((nn-1)/Nbands)+1 ! phonon2 momentum [1:nptk]
       ! ----------- end sampling method add -----------
          qprime=IJK(:,ii)
          omegap=energy(ii,j)
          !--------BEGIN absorption process-----------
          qdprime=q+qprime
          qdprime=modulo(qdprime,Ngrid)
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3)) ! phonon3 momentum [1:nptk]
          omegadp=energy(ss,k)
          if ((omegap.ne.0).and.(omegadp.ne.0)) then
             sigma=scalebroad*base_sigma(&
                velocity(ii,j,:)-&
                velocity(ss,k,:))
             if(abs(omega+omegap-omegadp).le.(2.d0*sigma)) then
                N_plus=N_plus+1
                P_plus=P_plus+&
                   exp(-(omega+omegap-omegadp)**2/(sigma**2))/&
                   (sigma*sqrt(Pi)*nptk**2*nbands**3)
             end if
          end if
          !--------END absorption process-------------!
       end do ! iter
     end if
    ! ----------- sampling method add -----------
     N_plus = INT(REAL(N_plus, 8)/num_sample_process_3ph_phase_space*Nbands*nptk*Nbands, 4)   ! must first do the division and then do the multiplication. Otherwise it will overflow the max of int32
     P_plus = P_plus*Nbands*nptk*Nbands/num_sample_process_3ph_phase_space
    ! ----------- end sampling method add -----------
 
   end subroutine NP_plus_sample
 
   ! Same as NP_plus, but for emission processes using sampling method.
   subroutine NP_minus_sample(mm,energy,velocity,Nlist,List,IJK,N_minus,P_minus)
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
 
    ! ----------- sampling method add -----------
     integer(kind=4) :: nn, rand_num, iter
     real :: rand_matrix(num_sample_process_3ph_phase_space),rand
    ! ----------- end sampling method add -----------
 
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
       ! ----------- sampling method add -----------
       ! Note: 
       ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
       ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
       ! phonon3: ss is momentum [1:nptk], k is energy [1:Nband]
       call random_seed()
       call random_number(rand_matrix)
       do iter=1,num_sample_process_3ph_phase_space
          rand = rand_matrix(iter-1)
          rand_num = FLOOR(INT(Nbands*nptk*Nbands)*rand)+1  ! generate random integer in [1,Nbands*nptk*Nbands]
          nn = int((rand_num-1)/Nbands)+1 ! get a random phonon number for 2nd phonon
          k = modulo(rand_num-1,Nbands)+1 ! phonon3 branch [1:Nband]
          j=modulo(nn-1,Nbands)+1 ! phonon2 branch [1:Nband]
          ii=int((nn-1)/Nbands)+1 ! phonon2 momentum [1:nptk]
       ! ----------- end sampling method add -----------
          qprime=IJK(:,ii)
          omegap=energy(ii,j)
          !--------BEGIN emission process-----------
          qdprime=q-qprime
          qdprime=modulo(qdprime,Ngrid)
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3)) ! phonon3 momentum [1:nptk]
          omegadp=energy(ss,k)
          if ((omegap.ne.0).and.(omegadp.ne.0)) then
             sigma=scalebroad*base_sigma(&
                velocity(ii,j,:)-&
                velocity(ss,k,:))
             if(abs(omega-omegap-omegadp).le.(2.d0*sigma)) then
                N_minus=N_minus+1
                P_minus=P_minus+&
                   exp(-(omega-omegap-omegadp)**2/(sigma**2))/&
                   (sigma*sqrt(Pi)*nptk**2*nbands**3)
             end if
          end if
       !--------END emission process-------------!
       end do ! iter
     end if
    ! ----------- end sampling method add -----------
     N_minus = INT(REAL(N_minus, 8)/num_sample_process_3ph_phase_space*Nbands*nptk*Nbands, 4)   ! must first do the division and then do the multiplication. Otherwise it will overflow the max of int32
     P_minus = P_minus*Nbands*nptk*Nbands/num_sample_process_3ph_phase_space
    ! ----------- end sampling method add -----------
 
   end subroutine NP_minus_sample
 
  
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
   
      integer(kind=4) :: mm,i,ll

      integer(kind=4) :: N_plus_reduce(Nlist*Nbands)
      integer(kind=4) :: N_minus_reduce(Nlist*Nbands)
      real(kind=8) :: Pspace_plus_reduce(Nlist*Nbands)
      real(kind=8) :: Pspace_minus_reduce(Nlist*Nbands)

      integer(kind=4) :: N_temp
      real(kind=8) :: Pspace_temp


      integer(kind=4) :: chunkstart, chunkend, chunkid

      Pspace_plus_total=0.d0
      Pspace_plus_reduce=0.d0
      Pspace_minus_total=0.d0
      Pspace_minus_reduce=0.d0
      Pspace_temp = 0.d0

      N_plus=0
      N_minus=0
      N_plus_reduce=0
      N_minus_reduce=0
      N_temp=0

      do chunkid=myid+1,numchunk,numprocs
         chunkstart = (chunkid-1)*chunksize+1
         chunkend = chunkid*chunksize
         if (chunkend.gt.Nbands*Nlist) chunkend = Nbands*Nlist
         !$OMP PARALLEL DO default(shared) schedule(dynamic,1) private(mm,i,ll,N_temp,Pspace_temp)
         do mm=chunkstart,chunkend
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1
            if (energy(List(int((mm-1)/Nbands)+1),modulo(mm-1,Nbands)+1).le.omega_max) then
               if (num_sample_process_3ph_phase_space==-1) then ! do not sample
                  call NP_plus(mm,energy,velocity,Nlist,List,IJK,&
                     N_temp,Pspace_temp)
                  N_plus_reduce(mm) = N_temp
                  Pspace_plus_reduce(mm) = Pspace_temp
                  call NP_minus(mm,energy,velocity,Nlist,List,IJK,&
                     N_temp,Pspace_temp)
                  N_minus_reduce(mm) = N_temp
                  Pspace_minus_reduce(mm) = Pspace_temp
               else ! do sampling method
                  call NP_plus_sample(mm,energy,velocity,Nlist,List,IJK,&
                     N_temp,Pspace_temp)
                  N_plus_reduce(mm) = N_temp
                  Pspace_plus_reduce(mm) = Pspace_temp
                  call NP_minus_sample(mm,energy,velocity,Nlist,List,IJK,&
                     N_temp,Pspace_temp)
                  N_minus_reduce(mm) = N_temp
                  Pspace_minus_reduce(mm) = Pspace_temp
               endif
            endif
         end do
         !$OMP END PARALLEL DO
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




   ! Compute the number of allowed absorption processes and their contribution
   ! to phase space.
   subroutine Ind_only_plus(mm, energy, velocity, Nlist, List, IJK, N_plus, &
         Indof2ndPhonon_plus, Indof3rdPhonon_plus)
      implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     integer(kind=4),intent(in) :: N_plus
     integer(kind=4),intent(out) :: Indof2ndPhonon_plus(N_plus),Indof3rdPhonon_plus(N_plus)

     integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss
     real(kind=8) :: sigma
     real(kind=8) :: omega,omegap,omegadp

     integer(kind=4) :: N_plus_count
 
     do ii=0,Ngrid(1)-1        ! G1 direction
        do jj=0,Ngrid(2)-1     ! G2 direction
           do kk=0,Ngrid(3)-1  ! G3 direction
              Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
           end do
        end do
     end do
     N_plus_count=0

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
                       N_plus_count=N_plus_count+1
                       Indof2ndPhonon_plus(N_plus_count)=(ii-1)*Nbands+j
                       Indof3rdPhonon_plus(N_plus_count)=(ss-1)*Nbands+k
                    end if
                 end if
              end do ! k
              !--------END absorption process-------------!
           end do ! ii
        end do  ! j
     end if
   end subroutine Ind_only_plus



      ! Same as NP_plus, but for emission processes.
   subroutine Ind_only_minus(mm, energy, velocity, Nlist, List, IJK, N_minus, &
            Indof2ndPhonon_minus, Indof3rdPhonon_minus)
            
      implicit none
      integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
      real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
      integer(kind=4),intent(in) :: N_minus
      integer(kind=4),intent(out) :: Indof2ndPhonon_minus(N_minus),Indof3rdPhonon_minus(N_minus)
   
      integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k
      integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
      integer(kind=4) :: ii,jj,kk,ll,ss
      real(kind=8) :: sigma
      real(kind=8) :: omega,omegap,omegadp

      integer(kind=4) :: N_minus_count

   
      do ii=0,Ngrid(1)-1        ! G1 direction
         do jj=0,Ngrid(2)-1     ! G2 direction
            do kk=0,Ngrid(3)-1  ! G3 direction
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
                        N_minus_count=N_minus_count+1
                        Indof2ndPhonon_minus(N_minus_count)=(ii-1)*Nbands+j
                        Indof3rdPhonon_minus(N_minus_count)=(ss-1)*Nbands+k
                     end if
                  end if
               end do ! k
               !--------END emission process-------------!
            end do ! ii
         end do  ! j
      end if
   end subroutine Ind_only_minus
   



   ! Wrapper around NP_plus and NP_minus that splits the work among processors.
   subroutine Ind_only_driver(energy,velocity,Nlist,List,IJK,&
        N_plus,N_minus,&
         Indof2ndPhonon_plus, Indof3rdPhonon_plus,&
        Indof2ndPhonon_minus, Indof3rdPhonon_minus)
      implicit none
   
      include "mpif.h"
   
      real(kind=8),intent(in) :: energy(nptk,nbands)
      real(kind=8),intent(in) :: velocity(nptk,nbands,3)
      integer(kind=4),intent(in) :: NList
      integer(kind=4),intent(in) :: List(Nlist)
      integer(kind=4),intent(in) :: IJK(3,nptk)
      integer(kind=4),intent(in) :: N_plus(Nlist*Nbands)
      integer(kind=4),intent(in) :: N_minus(Nlist*Nbands)
      integer(kind=4),intent(out) :: Indof2ndPhonon_plus(:)
      integer(kind=4),intent(out) :: Indof3rdPhonon_plus(:)
      integer(kind=4),intent(out) :: Indof2ndPhonon_minus(:)
      integer(kind=4),intent(out) :: Indof3rdPhonon_minus(:)

      integer(kind=4) :: mm,i,ll

      ! reduced variables
      integer(kind=4),allocatable :: Indof2ndPhonon_plus_reduce(:)
      integer(kind=4),allocatable :: Indof3rdPhonon_plus_reduce(:)
      integer(kind=4),allocatable :: Indof2ndPhonon_minus_reduce(:)
      integer(kind=4),allocatable :: Indof3rdPhonon_minus_reduce(:)

      !   temporary variables
      integer(kind=4),allocatable :: Indof2ndPhonon(:)
      integer(kind=4),allocatable :: Indof3rdPhonon(:)

      integer(kind=4) :: maxsize
      integer(kind=4) :: Ntotal_plus
      integer(kind=4) :: Ntotal_minus
      integer(kind=4) :: Naccum_plus(Nbands*Nlist)
      integer(kind=4) :: Naccum_minus(Nbands*Nlist)


      integer(kind=4) :: chunkstart, chunkend, chunkid

      maxsize=max(maxval(N_plus),maxval(N_minus))

      !   allocate temporary arrays
      allocate(Indof2ndPhonon(maxsize))
      allocate(Indof3rdPhonon(maxsize))
      Naccum_plus(1)=0
      Naccum_minus(1)=0
      do mm=2,Nbands*Nlist
         Naccum_plus(mm)=Naccum_plus(mm-1)+N_plus(mm-1)
         Naccum_minus(mm)=Naccum_minus(mm-1)+N_minus(mm-1)
      end do
      Ntotal_plus=sum(N_plus)
      Ntotal_minus=sum(N_minus)
   
      !  allocate reduced output arrays
      allocate(Indof2ndPhonon_plus_reduce(Ntotal_plus))
      allocate(Indof3rdPhonon_plus_reduce(Ntotal_plus))
      allocate(Indof2ndPhonon_minus_reduce(Ntotal_minus))
      allocate(Indof3rdPhonon_minus_reduce(Ntotal_minus))

      ! initialize the reduced output arrays
      Indof2ndPhonon_plus_reduce=0
      Indof3rdPhonon_plus_reduce=0
      Indof2ndPhonon_minus_reduce=0
      Indof3rdPhonon_minus_reduce=0

      ! initialize for the temporary variables
      Indof2ndPhonon=0
      Indof3rdPhonon=0

      do chunkid=myid+1,numchunk,numprocs
         chunkstart = (chunkid-1)*chunksize+1
         chunkend = chunkid*chunksize
         if (chunkend.gt.Nbands*Nlist) chunkend = Nbands*Nlist
         !$OMP PARALLEL DO default(shared) schedule(dynamic,1) private(mm,i,ll,Indof2ndPhonon,Indof3rdPhonon)
         do mm=chunkstart,chunkend
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1
            if (energy(List(int((mm-1)/Nbands)+1),modulo(mm-1,Nbands)+1).le.omega_max) then
               call Ind_only_plus(mm,energy,velocity,Nlist,List,IJK,&
                  N_plus(mm),&
                  Indof2ndPhonon(1:N_plus(mm)),Indof3rdPhonon(1:N_plus(mm)))
               Indof2ndPhonon_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
                     Indof2ndPhonon(1:N_plus(mm))
               Indof3rdPhonon_plus_reduce((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
                     Indof3rdPhonon(1:N_plus(mm))


               call Ind_only_minus(mm,energy,velocity,Nlist,List,IJK,&
                  N_minus(mm),&
                  Indof2ndPhonon(1:N_minus(mm)),Indof3rdPhonon(1:N_minus(mm)))
               Indof2ndPhonon_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
                     Indof2ndPhonon(1:N_minus(mm))
               Indof3rdPhonon_minus_reduce((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
                     Indof3rdPhonon(1:N_minus(mm))

            endif
         end do
         !$OMP END PARALLEL DO
      end do

      call MPI_ALLREDUCE(Indof2ndPhonon_plus_reduce,Indof2ndPhonon_plus,&
            Ntotal_plus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
      call MPI_ALLREDUCE(Indof3rdPhonon_plus_reduce,Indof3rdPhonon_plus,&
            Ntotal_plus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
      call MPI_ALLREDUCE(Indof2ndPhonon_minus_reduce,Indof2ndPhonon_minus,&
            Ntotal_minus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)
      call MPI_ALLREDUCE(Indof3rdPhonon_minus_reduce,Indof3rdPhonon_minus,&
            Ntotal_minus,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ll)

   end subroutine Ind_only_driver
 
 
   ! RTA-only version of Ind_plus.
   subroutine RTA_plus(mm,energy,velocity,eigenvect,Nlist,List,&
        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
        Gamma_plus,WP3_plus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_plus,WP3_plus
 
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
         do ii=1,nptk
            do j=1,Nbands
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
   end subroutine RTA_plus
 
   ! RTA-only version of Ind_minus.
   subroutine RTA_minus(mm,energy,velocity,eigenvect,Nlist,List,&
        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
        Gamma_minus,WP3_minus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_minus,WP3_minus
 
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
         do ii=1,nptk
            do j=1,Nbands
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
 
   end subroutine RTA_minus
 
   ! RTA-only version of Ind_plus using sampling method.
   subroutine RTA_plus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
        Gamma_plus,WP3_plus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_plus,WP3_plus
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss
     real(kind=8) :: sigma
     real(kind=8) :: fBEprime,fBEdprime
     real(kind=8) :: omega,omegap,omegadp
     real(kind=8) :: realq(3),realqprime(3),realqdprime(3)
     real(kind=8) :: WP3
     complex(kind=8) :: Vp
 
    ! ----------- sampling method add -----------
     integer(kind=4) :: nn, rand_num, iter
     real :: rand_matrix(num_sample_process_3ph),rand
    ! ----------- end sampling method add -----------
 
     Gamma_plus=0.d00
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
 
       ! ----------- sampling method add -----------
       ! Note: 
       ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
       ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
       ! phonon3: ss is momentum [1:nptk], k is energy [1:Nband]
       call random_seed()
       call random_number(rand_matrix)
       do iter=1,num_sample_process_3ph
          rand = rand_matrix(iter-1)
          rand_num = FLOOR(INT(Nbands*nptk*Nbands)*rand)+1  ! generate random integer in [1,Nbands*nptk*Nbands]
          nn = int((rand_num-1)/Nbands)+1 ! get a random phonon number for 2nd phonon
          k = modulo(rand_num-1,Nbands)+1 ! phonon3 branch [1:Nband]
          j=modulo(nn-1,Nbands)+1 ! phonon2 branch [1:Nband]
          ii=int((nn-1)/Nbands)+1 ! phonon2 momentum [1:nptk]
       ! ----------- end sampling method add -----------
 
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,qprime/dble(ngrid))
          omegap=energy(ii,j)
          fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
          !--------BEGIN absorption process-----------
          qdprime=q+qprime
          qdprime=modulo(qdprime,Ngrid)
          realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3)) ! phonon3 momentum [1:nptk]
          omegadp=energy(ss,k)
          if ((omegap.ne.0).and.(omegadp.ne.0)) then ! obey energy conservation
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
                endif
             end if
          end if ! end energy conservation
          !--------END absorption process-------------!
       end do ! iter
        WP3_plus=WP3_plus/nptk
     end if
    ! ----------- sampling method add -----------
     Gamma_plus = Gamma_plus*Nbands*nptk*Nbands/num_sample_process_3ph
     WP3_plus=WP3_plus*Nbands*nptk*Nbands/num_sample_process_3ph
    ! ----------- end sampling method add -----------
 
     Gamma_plus=Gamma_plus*5.60626442*1.d8/nptk ! THz
   end subroutine RTA_plus_sample
 
   ! RTA-only version of Ind_minus using sampling method.
   subroutine RTA_minus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
        Gamma_minus,WP3_minus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_minus,WP3_minus
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),i,j,k,N_minus_count
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss
     real(kind=8) :: sigma
     real(kind=8) :: fBEprime,fBEdprime
     real(kind=8) ::  omega,omegap,omegadp
     real(kind=8) :: realqprime(3),realqdprime(3)
     real(kind=8) :: WP3
     complex(kind=8) :: Vp
 
    ! ----------- sampling method add -----------
     integer(kind=4) :: nn, rand_num, iter
     real :: rand_matrix(num_sample_process_3ph),rand
    ! ----------- end sampling method add -----------
 
     Gamma_minus=0.d00
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
       ! ----------- sampling method add -----------
       ! Note: 
       ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
       ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
       ! phonon3: ss is momentum [1:nptk], k is energy [1:Nband]
       call random_seed()
       call random_number(rand_matrix)
       do iter=1,num_sample_process_3ph
          rand = rand_matrix(iter-1)
          rand_num = FLOOR(INT(Nbands*nptk*Nbands)*rand)+1  ! generate random integer in [1,Nbands*nptk*Nbands]
          nn = int((rand_num-1)/Nbands)+1 ! get a random phonon number for 2nd phonon
          k = modulo(rand_num-1,Nbands)+1 ! phonon3 branch [1:Nband]
          j=modulo(nn-1,Nbands)+1 ! phonon2 branch [1:Nband]
          ii=int((nn-1)/Nbands)+1 ! phonon2 momentum [1:nptk]
       ! ----------- end sampling method add -----------
 
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,qprime/dble(ngrid))
          omegap=energy(ii,j)
          fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
          !--------BEGIN emission process-----------
          qdprime=q-qprime
          qdprime=modulo(qdprime,Ngrid)
          realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
          ss=Index_N(qdprime(1),qdprime(2),qdprime(3)) ! phonon3 momentum [1:nptk]
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
                endif
             end if
          end if ! end energy conservation
          !--------END emission process-------------
       end do ! iter
        WP3_minus=WP3_minus*5.d-1/nptk
     end if
 
    ! ----------- sampling method add -----------
     Gamma_minus = Gamma_minus*Nbands*nptk*Nbands/num_sample_process_3ph
     WP3_minus = WP3_minus*Nbands*nptk*Nbands/num_sample_process_3ph
    ! ----------- end sampling method add -----------
 
     Gamma_minus=Gamma_minus*5.60626442*1.d8/nptk
   end subroutine RTA_minus_sample
 
   ! Wrapper around 3ph RTA_plus and RTA_minus that splits the work among processors.
   subroutine RTA_driver(energy,velocity,eigenvect,Nlist,List,IJK,&
        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,rate_scatt,rate_scatt_plus,rate_scatt_minus,WP3_plus,WP3_minus)
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
      real(kind=8),intent(out) :: WP3_plus(Nbands,Nlist)
      real(kind=8),intent(out) :: WP3_minus(Nbands,Nlist)
   
      integer(kind=4) :: i
      integer(kind=4) :: ll
      integer(kind=4) :: mm
      real(kind=8) :: Gamma_plus,Gamma_minus
      real(kind=8) :: rate_scatt_plus_reduce(Nbands,Nlist),rate_scatt_minus_reduce(Nbands,Nlist)
      real(kind=8) :: WP3_plus_reduce(Nbands*Nlist)
      real(kind=8) :: WP3_minus_reduce(Nbands*Nlist)

      real(kind=8) :: WP3_temp ! variable for parallel computing

      integer(kind=4) :: chunkstart, chunkend, chunkid

      rate_scatt=0.d00
      rate_scatt_plus=0.d00
      rate_scatt_minus=0.d00
      rate_scatt_plus_reduce=0.d00
      rate_scatt_minus_reduce=0.d00

      WP3_plus=0.d00
      WP3_minus=0.d00
      WP3_plus_reduce=0.d00
      WP3_minus_reduce=0.d00

      do chunkid=myid+1,numchunk,numprocs  ! when numchunk=numprocs, each mpi process only calculate one chunk
         chunkstart = (chunkid-1)*chunksize+1
         chunkend = chunkid*chunksize
         if (chunkend.gt.Nbands*Nlist) chunkend = Nbands*Nlist
         ! !$OMP PARALLEL DO default(none) schedule(dynamic,1)  shared(myid,chunkstart,chunkend,Nbands,NList) &
         ! !$OMP & shared(num_sample_process_3ph,eigenvect,Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k) &
         ! !$OMP & shared(energy,velocity,List,IJK,omega_max) &
         ! !$OMP & shared(WP3_plus_reduce,WP3_minus_reduce,rate_scatt_plus_reduce,rate_scatt_minus_reduce) &
         ! !$OMP & private(mm,i,ll,WP3_temp,Gamma_plus,Gamma_minus)
      !$OMP PARALLEL DO default(shared) schedule(dynamic,1) private(mm,i,ll,WP3_temp,Gamma_plus,Gamma_minus)
         do mm=chunkstart,chunkend
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1
            if (energy(List(ll),i).le.omega_max) then
               if (num_sample_process_3ph==-1) then ! do not sample
                  call RTA_plus(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                        Gamma_plus,WP3_temp)
                  WP3_plus_reduce(mm) = WP3_temp
                  rate_scatt_plus_reduce(i,ll)=Gamma_plus
                  call RTA_minus(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                        Gamma_minus, WP3_temp)
                  WP3_minus_reduce(mm) = WP3_temp
                  rate_scatt_minus_reduce(i,ll)=Gamma_minus*5.D-1
               else ! do sampling method
                  call RTA_plus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                        Gamma_plus,WP3_temp)
                  WP3_plus_reduce(mm) = WP3_temp
                  rate_scatt_plus_reduce(i,ll)=Gamma_plus
                  call RTA_minus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                        Gamma_minus, WP3_temp)
                  WP3_minus_reduce(mm) = WP3_temp
                  rate_scatt_minus_reduce(i,ll)=Gamma_minus*5.D-1
               endif
      
            endif
         end do
         !$OMP END PARALLEL DO
      end do
   
      call MPI_ALLREDUCE(rate_scatt_plus_reduce,rate_scatt_plus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
      call MPI_ALLREDUCE(rate_scatt_minus_reduce,rate_scatt_minus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
      call MPI_ALLREDUCE(WP3_plus_reduce,WP3_plus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
      call MPI_ALLREDUCE(WP3_minus_reduce,WP3_minus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

      rate_scatt=rate_scatt_plus+rate_scatt_minus
   end subroutine RTA_driver

   
   subroutine RTA_driver_GPU_using_Ind(energy, velocity, eigenvect, Nlist, List, IJK, &
      Ntri, Phi, R_j, R_k, Index_i, Index_j, Index_k, &
      N_plus, N_minus,&
      Indof2ndPhonon_plus, Indof2ndPhonon_minus, Indof3rdPhonon_plus, Indof3rdPhonon_minus,&
      rate_scatt, rate_scatt_plus, rate_scatt_minus, WP3_plus, WP3_minus)
      implicit none
      include "mpif.h"

      ! Input variables
      real(kind=8), intent(in)    :: energy(nptk, nbands)
      real(kind=8), intent(in)    :: velocity(nptk, nbands, 3)
      complex(kind=8), intent(in) :: eigenvect(nptk, Nbands, Nbands)
      integer(kind=4), intent(in) :: Nlist, List(Nlist)
      integer(kind=4), intent(in) :: IJK(3, nptk)
      integer(kind=4), intent(in) :: Ntri
      real(kind=8), intent(in)    :: Phi(3, 3, 3, Ntri)
      real(kind=8), intent(in)    :: R_j(3, Ntri)
      real(kind=8), intent(in)    :: R_k(3, Ntri)
      integer(kind=4), intent(in) :: Index_i(Ntri)
      integer(kind=4), intent(in) :: Index_j(Ntri)
      integer(kind=4), intent(in) :: Index_k(Ntri)

      integer(kind=4), intent(in) :: N_plus(Nbands * Nlist)
      integer(kind=4), intent(in) :: N_minus(Nbands * Nlist)

      integer(kind=4), intent(in) :: Indof2ndPhonon_plus(:), Indof3rdPhonon_plus(:)
      integer(kind=4), intent(in) :: Indof2ndPhonon_minus(:), Indof3rdPhonon_minus(:)

      ! Output variables
      real(kind=8), intent(out) :: rate_scatt(Nbands, Nlist)
      real(kind=8), intent(out) :: rate_scatt_plus(Nbands, Nlist)
      real(kind=8), intent(out) :: rate_scatt_minus(Nbands, Nlist)
      real(kind=8), intent(out) :: WP3_plus(Nbands, Nlist)  ! WP3_plus_array from plus driver
      real(kind=8), intent(out) :: WP3_minus(Nbands, Nlist) ! WP3_minus_array from minus driver

      ! Local variables
      integer(kind=4) :: i, ll, mm
      real(kind=8) :: Gamma_plus, Gamma_minus

      if (num_sample_process_3ph .eq. -1) then
         if (myid .eq. 0) then
            call RTA_plus_driver_GPU_using_Ind(energy, velocity, eigenvect, Nlist, List, &
                  Ntri, Phi, R_j, R_k, Index_i, Index_j, Index_k, IJK, &
                  N_plus, &
                  Indof2ndPhonon_plus, Indof3rdPhonon_plus, &
                  rate_scatt_plus, WP3_plus)
            call RTA_minus_driver_GPU_using_Ind(energy, velocity, eigenvect, Nlist, List, &
               Ntri, Phi, R_j, R_k, Index_i, Index_j, Index_k, IJK, &
               N_minus, &
               Indof2ndPhonon_minus, Indof3rdPhonon_minus, &
               rate_scatt_minus, WP3_minus)
            rate_scatt_minus=rate_scatt_minus*5.D-1
         end if
      else ! do sampling method
         print *, "Error: Sampling method is not implemented for GPU version."
      endif
      rate_scatt = rate_scatt_plus + rate_scatt_minus

   end subroutine RTA_driver_GPU_using_Ind


   subroutine RTA_plus_driver_GPU_using_Ind(energy, velocity, eigenvect, Nlist, List, &
      Ntri, Phi, R_j, R_k, Index_i, Index_j, Index_k, IJK, &
      N_plus, &
      Indof2ndPhonon_plus, Indof3rdPhonon_plus, &
      rate_scatt_plus, WP3_plus_array)
      implicit none
      include "mpif.h"

      ! Input variables
      integer(kind=4), intent(in)    :: Nlist, List(Nlist), IJK(3, nptk)
      integer(kind=4), intent(in)    :: Ntri
      integer(kind=4), intent(in)    :: Index_i(Ntri), Index_j(Ntri), Index_k(Ntri)
      real(kind=8),    intent(in)    :: energy(nptk, Nbands), velocity(nptk, Nbands, 3)
      complex(kind=8), intent(in)    :: eigenvect(nptk, Nbands, Nbands)
      real(kind=8),    intent(in)    :: Phi(3, 3, 3, Ntri), R_j(3, Ntri), R_k(3, Ntri)

      ! Precomputed scattering process data:
      integer(kind=4), intent(in) :: N_plus(Nbands * Nlist)
      integer(kind=4), intent(in) :: Indof2ndPhonon_plus(:), Indof3rdPhonon_plus(:)

      ! Output variables
      real(kind=8), intent(out) :: rate_scatt_plus(Nbands, Nlist)
      real(kind=8), intent(out) :: WP3_plus_array(Nbands, Nlist)

      ! Local variables
      integer(kind=4) :: mm, i, ll
      real(kind=8) :: Gamma_plus, WP3_plus, WP3_val
      real(kind=8) :: omega, omegap, omegadp, sigma
      real(kind=8) :: fBEprime, fBEdprime
      integer(kind=4) :: q(3), qprime(3), qdprime(3)
      real(kind=8) :: realq(3), realqprime(3), realqdprime(3)
      integer(kind=4) :: ii, j, ss, k
      integer(kind=8) :: Naccum_plus(Nbands * Nlist)
      integer(kind=8) :: indnow

      ! For base_sigma calculation
      real(kind=8) :: v_base_sigma(3), base_sigma_value

      ! For inline Vp_plus calculation
      integer(kind=4) :: ll_Vp, rr_Vp, ss_Vp, tt_Vp
      integer(kind=4) :: temp_index1, temp_index2, temp_index3
      complex(kind=8) :: prefactor, Vp0, Vp_plus_inline, Vp

      ! Array for mapping q-vectors to a single index
      integer(kind=4), allocatable :: Index_N(:,:,:)
      integer(kind=4) :: d1, d2, d3

      ! Build the Index_N mapping array over the 3D grid.
      d1 = Ngrid(1)
      d2 = Ngrid(2)
      d3 = Ngrid(3)
      allocate(Index_N(0:d1-1, 0:d2-1, 0:d3-1))
      do ii = 0, d1-1
         do j = 0, d2-1
            do k = 0, d3-1
               Index_N(ii, j, k) = (k*d2 + j)*d1 + ii + 1
            end do
         end do
      end do

      ! Build an accumulation array so that for each (band, list) pair (indexed by mm)
      ! we know where its allowed process data begins in the index arrays.
      Naccum_plus(1) = 0
      do mm = 2, Nbands * Nlist
         Naccum_plus(mm) = Naccum_plus(mm-1) + N_plus(mm-1)
      end do

      ! Initialize output arrays.
      rate_scatt_plus = 0.d0
      WP3_plus_array  = 0.d0


#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc data copyin(energy, velocity, IJK, List, Nlist, Index_N, Ntri, Phi, R_j, R_k) &
      !$acc&  copyin(Index_i, Index_j, Index_k)  &
      !$acc& copyin(eigenvect)  &
      !$acc& copyin(T, scalebroad, Ngrid, nptk, Nbands, rlattvec, masses, types) &
      !$acc& copyin(omega_max, onlyharmonic)  &
      !$acc copyin(N_plus, Naccum_plus, Indof2ndPhonon_plus, Indof3rdPhonon_plus)  &
      !$acc copy(rate_scatt_plus, WP3_plus_array)
      !$acc parallel loop gang default(present) &
      !$acc& vector_length(VEC_LEN) &
      !$acc& num_gangs(NUM_GANGS) &
      !$acc private(Gamma_plus, WP3_plus, i, ll)
#endif
      do mm = 1, Nbands * Nlist
         Gamma_plus = 0.d0
         WP3_plus   = 0.d0
         i  = modulo(mm-1, Nbands) + 1
         ll = int((mm-1)/Nbands) + 1

#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
      !$acc data copyin(energy, velocity, IJK, List, Nlist, Index_N, Ntri, Phi, R_j, R_k) &
      !$acc&  copyin(i,ll)  &
      !$acc&  copyin(Index_i, Index_j, Index_k)  &
      !$acc& copyin(eigenvect)  &
      !$acc& copyin(T, scalebroad, Ngrid, nptk, Nbands, rlattvec, masses, types) &
      !$acc& copyin(omega_max, onlyharmonic)  &
      !$acc& copyin(N_plus, Naccum_plus)  &
      !$acc& copyin(Indof2ndPhonon_plus(Naccum_plus(mm) + 1:Naccum_plus(mm) + N_plus(mm))) &
      !$acc& copyin(Indof3rdPhonon_plus(Naccum_plus(mm) + 1:Naccum_plus(mm) + N_plus(mm)))  &
      !$acc copy(rate_scatt_plus, WP3_plus_array)
      !$acc parallel loop gang default(present) reduction(+:Gamma_plus, WP3_plus) &
      !$acc& vector_length(VEC_LEN) &
      !$acc& num_gangs(NUM_GANGS) &
      !$acc private(Gamma_plus, WP3_plus) &
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
         !$acc loop vector reduction(+:Gamma_plus, WP3_plus) &
#endif
         !$acc private(ii, j, ss, k, q, qprime, realq, realqprime, qdprime, realqdprime, &
         !$acc     omega, omegap, omegadp, sigma, fBEprime, fBEdprime, WP3_val, &
         !$acc     v_base_sigma, base_sigma_value)
         do indnow = Naccum_plus(mm) + 1, Naccum_plus(mm) + N_plus(mm)
         
            ! Recover the indices for this scattering process:
            ii = (Indof2ndPhonon_plus(indnow) - 1) / Nbands + 1
            j  = mod(Indof2ndPhonon_plus(indnow) - 1, Nbands) + 1
            ss = (Indof3rdPhonon_plus(indnow) - 1) / Nbands + 1
            k  = mod(Indof3rdPhonon_plus(indnow) - 1, Nbands) + 1

            q = IJK(:, List(ll))
            realq = matmul(rlattvec, q / dble(Ngrid))
            omega = energy(List(ll), i)

            qprime = IJK(:, ii)
            realqprime = matmul(rlattvec, qprime / dble(Ngrid))

            qdprime = modulo(q + qprime, Ngrid)
            realqdprime = matmul(rlattvec, qdprime / dble(Ngrid))

            omegap = energy(ii, j)
            omegadp = energy(ss, k)

            ! ---------- Sigma calculation inline ----------
            v_base_sigma = velocity(ii, j, :) - velocity(ss, k, :)
            base_sigma_value = max( (dot_product(rlattvec(:,1), v_base_sigma) / Ngrid(1))**2, &
                                    (dot_product(rlattvec(:,2), v_base_sigma) / Ngrid(2))**2, &
                                    (dot_product(rlattvec(:,3), v_base_sigma) / Ngrid(3))**2 )
            base_sigma_value = sqrt(base_sigma_value / 2.d0)
            sigma = scalebroad * base_sigma_value
            ! ---------- End Sigma calculation inline ----------

            fBEprime = 1.d0 / (exp(hbar * omegap / (Kb*T)) - 1.d0)
            fBEdprime = 1.d0 / (exp(hbar * omegadp / (Kb*T)) - 1.d0)

            WP3_val = (fBEprime - fBEdprime) * &
                        exp(- (omega + omegap - omegadp)**2 / (sigma**2)) / &
                        (sigma * sqrt(Pi) * (omega * omegap * omegadp))
            WP3_plus = WP3_plus + WP3_val

            if (.not. onlyharmonic) then
               ! ---------- Vp_plus calculation inline ----------
               Vp_plus_inline = (0.d0, 0.d0)
#ifdef GPU_ALL_MODE_PARALLELIZATION
               !$acc loop reduction(+:Vp_plus_inline) &
#endif
#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
               !$acc loop vector reduction(+:Vp_plus_inline) &
#endif
               !$acc private(rr_Vp, ss_Vp, tt_Vp, ll_Vp, Vp0, prefactor, temp_index1, temp_index2, temp_index3)
               do ll_Vp = 1, Ntri
                  prefactor = 1.d0 / sqrt( masses(types(Index_i(ll_Vp))) * &
                                             masses(types(Index_j(ll_Vp))) * &
                                             masses(types(Index_k(ll_Vp))) ) * &
                              cmplx(cos(dot_product(realqprime, R_j(:, ll_Vp))), &
                                    sin(dot_product(realqprime, R_j(:, ll_Vp))), kind=8) * &
                              cmplx(cos(-dot_product(realqdprime, R_k(:, ll_Vp))), &
                                    sin(-dot_product(realqdprime, R_k(:, ll_Vp))), kind=8)
                  Vp0 = (0.d0, 0.d0)
                  do rr_Vp = 1, 3
                     do ss_Vp = 1, 3
                        do tt_Vp = 1, 3
                           temp_index1 = tt_Vp + 3*(Index_i(ll_Vp) - 1)
                           temp_index2 = ss_Vp + 3*(Index_j(ll_Vp) - 1)
                           temp_index3 = rr_Vp + 3*(Index_k(ll_Vp) - 1)
                           Vp0 = Vp0 + Phi(tt_Vp, ss_Vp, rr_Vp, ll_Vp) * &
                                       eigenvect(List(ll), i, temp_index1) * &
                                       eigenvect(ii, j, temp_index2) * &
                                       conjg(eigenvect(ss, k, temp_index3))
                        end do
                     end do
                  end do
                  Vp_plus_inline = Vp_plus_inline + prefactor * Vp0
               end do
               Vp = Vp_plus_inline
               ! ---------- End Vp_plus calculation inline ----------
               Gamma_plus = Gamma_plus + (hbarp * pi / 4.d0) * WP3_val * abs(Vp)**2
            end if
         end do  ! End inner loop reduction
#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
            !$acc end parallel loop
            !$acc end data
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
            !$acc end loop
#endif
         WP3_plus = WP3_plus / nptk 
         rate_scatt_plus(i, ll) = Gamma_plus * 5.60626442d0 * 1.d8 / nptk
         WP3_plus_array(i, ll) = WP3_plus
      end do  ! End outer loop
#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc end parallel loop
      !$acc end data
#endif
      deallocate(Index_N)
   end subroutine RTA_plus_driver_GPU_using_Ind


   subroutine RTA_minus_driver_GPU_using_Ind(energy, velocity, eigenvect, Nlist, List, &
         Ntri, Phi, R_j, R_k, Index_i, Index_j, Index_k, IJK, &
         N_minus, &
         Indof2ndPhonon_minus, Indof3rdPhonon_minus, &
         rate_scatt_minus, WP3_minus_array)
      implicit none
      include "mpif.h"

      ! Input variables
      integer(kind=4), intent(in)    :: Nlist, List(Nlist), IJK(3, nptk)
      integer(kind=4), intent(in)    :: Ntri
      integer(kind=4), intent(in)    :: Index_i(Ntri), Index_j(Ntri), Index_k(Ntri)
      real(kind=8),    intent(in)    :: energy(nptk, Nbands), velocity(nptk, Nbands, 3)
      complex(kind=8), intent(in)    :: eigenvect(nptk, Nbands, Nbands)
      real(kind=8),    intent(in)    :: Phi(3, 3, 3, Ntri), R_j(3, Ntri), R_k(3, Ntri)

      ! Precomputed scattering process data:
      integer(kind=4), intent(in) :: N_minus(Nbands * Nlist)
      integer(kind=4), intent(in) :: Indof2ndPhonon_minus(:), Indof3rdPhonon_minus(:)

      ! Output variables
      real(kind=8), intent(out) :: rate_scatt_minus(Nbands, Nlist)
      real(kind=8), intent(out) :: WP3_minus_array(Nbands, Nlist)

      ! Local variables
      integer(kind=4) :: mm, i, ll
      real(kind=8) :: Gamma_minus, WP3_minus, WP3_val
      real(kind=8) :: omega, omegap, omegadp, sigma
      real(kind=8) :: fBEprime, fBEdprime
      integer(kind=4) :: q(3), qprime(3), qdprime(3)
      real(kind=8) :: realq(3), realqprime(3), realqdprime(3)
      integer(kind=4) :: ii, j, ss, k
      integer(kind=8) :: Naccum_minus(Nbands * Nlist)
      integer(kind=8) :: indnow

      ! For base_sigma calculation
      real(kind=8) :: v_base_sigma(3), base_sigma_value


      ! For inline Vp_minus calculation
      integer(kind=4) :: ll_Vp, rr_Vp, ss_Vp, tt_Vp
      integer(kind=4) :: temp_index1, temp_index2, temp_index3
      complex(kind=8) :: prefactor, Vp0, Vp_minus_inline, Vp

      ! Array for mapping q-vectors to a single index
      integer(kind=4), allocatable :: Index_N(:,:,:)
      integer(kind=4) :: d1, d2, d3

      ! Build the Index_N mapping array over the 3D grid.
      d1 = Ngrid(1)
      d2 = Ngrid(2)
      d3 = Ngrid(3)
      allocate(Index_N(0:d1-1, 0:d2-1, 0:d3-1))
      do ii = 0, d1-1
         do j = 0, d2-1
            do k = 0, d3-1
               Index_N(ii, j, k) = (k*d2 + j)*d1 + ii + 1
            end do
         end do
      end do

      ! Build an accumulation array so that for each (band, list) pair (indexed by mm)
      ! we know where its allowed process data begins in the precomputed index arrays.
      Naccum_minus(1) = 0
      do mm = 2, Nbands * Nlist
         Naccum_minus(mm) = Naccum_minus(mm-1) + N_minus(mm-1)
      end do

      ! Initialize output arrays.
      rate_scatt_minus = 0.d0
      WP3_minus_array  = 0.d0


#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc data copyin(energy, velocity, IJK, List, Nlist, Index_N, Ntri, Phi, R_j, R_k) &
      !$acc&  copyin(Index_i, Index_j, Index_k)  &
      !$acc& copyin(eigenvect)  &
      !$acc& copyin(T, scalebroad, Ngrid, nptk, Nbands, rlattvec, masses, types) &
      !$acc& copyin(omega_max, onlyharmonic)  &
      !$acc copyin(N_minus, Naccum_minus, Indof2ndPhonon_minus, Indof3rdPhonon_minus)  &
      !$acc copy(rate_scatt_minus, WP3_minus_array)
      !$acc parallel loop gang default(present) &
      !$acc& vector_length(VEC_LEN) &
      !$acc& num_gangs(NUM_GANGS) &
      !$acc private(Gamma_minus, WP3_minus, i, ll)
#endif
      do mm = 1, Nbands * Nlist
         Gamma_minus = 0.d0
         WP3_minus   = 0.d0
         i  = modulo(mm-1, Nbands) + 1
         ll = int((mm-1) / Nbands) + 1

#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
         !$acc data copyin(energy, velocity, IJK, List, Nlist, Index_N, Ntri, Phi, R_j, R_k) &
         !$acc&  copyin(i,ll)  &
         !$acc&  copyin(Index_i, Index_j, Index_k)  &
         !$acc& copyin(eigenvect)  &
         !$acc& copyin(T, scalebroad, Ngrid, nptk, Nbands, rlattvec, masses, types) &
         !$acc& copyin(omega_max, onlyharmonic)  &
         !$acc& copyin(N_minus, Naccum_minus)  &
         !$acc& copyin(Indof2ndPhonon_minus(Naccum_minus(mm) + 1:Naccum_minus(mm) + N_minus(mm))) &
         !$acc& copyin(Indof3rdPhonon_minus(Naccum_minus(mm) + 1:Naccum_minus(mm) + N_minus(mm)))  &
         !$acc copy(rate_scatt_minus, WP3_minus_array)
         !$acc parallel loop gang default(present) reduction(+:Gamma_minus, WP3_minus) &
         !$acc& vector_length(VEC_LEN) &
         !$acc& num_gangs(NUM_GANGS) &
         !$acc private(Gamma_minus, WP3_minus) &
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
         !$acc loop vector reduction(+:Gamma_minus, WP3_minus) &
#endif
         !$acc private(ii, j, ss, k, q, qprime, realq, realqprime, qdprime, realqdprime, &
         !$acc     omega, omegap, omegadp, sigma, fBEprime, fBEdprime, WP3_val, &
         !$acc     v_base_sigma, base_sigma_value)
         do indnow = Naccum_minus(mm) + 1, Naccum_minus(mm) + N_minus(mm)

            ! Recover the indices for this scattering process:
            ii = (Indof2ndPhonon_minus(indnow) - 1) / Nbands + 1
            j  = mod(Indof2ndPhonon_minus(indnow) - 1, Nbands) + 1
            ss = (Indof3rdPhonon_minus(indnow) - 1) / Nbands + 1
            k  = mod(Indof3rdPhonon_minus(indnow) - 1, Nbands) + 1

            q = IJK(:, List(ll))
            realq = matmul(rlattvec, q / dble(Ngrid))
            omega = energy(List(ll), i)

            qprime = IJK(:, ii)
            realqprime = matmul(rlattvec, qprime / dble(Ngrid))

            qdprime = modulo(q - qprime, Ngrid)
            realqdprime = matmul(rlattvec, qdprime / dble(Ngrid))

            omegap = energy(ii, j)
            omegadp = energy(ss, k)

            ! ---------- Sigma calculation inline ----------
            v_base_sigma = velocity(ii, j, :) - velocity(ss, k, :)
            base_sigma_value = max( (dot_product(rlattvec(:,1), v_base_sigma) / Ngrid(1))**2, &
                                    (dot_product(rlattvec(:,2), v_base_sigma) / Ngrid(2))**2, &
                                    (dot_product(rlattvec(:,3), v_base_sigma) / Ngrid(3))**2 )
            base_sigma_value = sqrt(base_sigma_value / 2.d0)
            sigma = scalebroad * base_sigma_value
            ! ---------- End Sigma calculation inline ----------

            fBEprime = 1.d0 / (exp(hbar * omegap / (Kb*T)) - 1.d0)
            fBEdprime = 1.d0 / (exp(hbar * omegadp / (Kb*T)) - 1.d0)

            WP3_val = (fBEprime + fBEdprime + 1) * &
                        exp(- (omega - omegap - omegadp)**2 / (sigma**2)) / &
                        (sigma * sqrt(Pi) * (omega * omegap * omegadp))
            WP3_minus = WP3_minus + WP3_val

            if (.not. onlyharmonic) then
               ! ---------- Vp_minus calculation inline ----------
               Vp_minus_inline = (0.d0, 0.d0)
#ifdef GPU_ALL_MODE_PARALLELIZATION
               !$acc loop reduction(+:Vp_minus_inline) &
#endif
#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
               !$acc loop vector reduction(+:Vp_minus_inline) &
#endif
               !$acc  private(rr_Vp, ss_Vp, tt_Vp, ll_Vp, Vp0, prefactor, temp_index1, temp_index2, temp_index3)
               do ll_Vp = 1, Ntri
                  prefactor = 1.d0 / sqrt( masses(types(Index_i(ll_Vp))) * &
                                             masses(types(Index_j(ll_Vp))) * &
                                             masses(types(Index_k(ll_Vp))) ) * &
                              cmplx(cos(-dot_product(realqprime, R_j(:, ll_Vp))), &
                                    sin(-dot_product(realqprime, R_j(:, ll_Vp))), kind=8) * &
                              cmplx(cos(-dot_product(realqdprime, R_k(:, ll_Vp))), &
                                    sin(-dot_product(realqdprime, R_k(:, ll_Vp))), kind=8)
                  Vp0 = (0.d0, 0.d0)
                  do rr_Vp = 1, 3
                     do ss_Vp = 1, 3
                        do tt_Vp = 1, 3
                           temp_index1 = tt_Vp + 3 * (Index_i(ll_Vp) - 1)
                           temp_index2 = ss_Vp + 3 * (Index_j(ll_Vp) - 1)
                           temp_index3 = rr_Vp + 3 * (Index_k(ll_Vp) - 1)
                           Vp0 = Vp0 + Phi(tt_Vp, ss_Vp, rr_Vp, ll_Vp) * &
                                       eigenvect(list(ll), i, temp_index1) * &
                                       conjg(eigenvect(ii, j, temp_index2)) * &
                                       conjg(eigenvect(ss, k, temp_index3))
                        end do
                     end do
                  end do
                  Vp_minus_inline = Vp_minus_inline + prefactor * Vp0
               end do
               Vp = Vp_minus_inline
               ! ---------- End Vp_minus calculation inline ----------
               Gamma_minus = Gamma_minus + (hbarp * pi / 4.d0) * WP3_val * abs(Vp)**2
            end if
         end do  ! End inner loop reduction
#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
         !$acc end parallel loop
         !$acc end data
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
         !$acc end loop
#endif
         WP3_minus = WP3_minus*5.d-1 / nptk 
         rate_scatt_minus(i, ll) = Gamma_minus * 5.60626442d0 * 1.d8 / nptk 
         WP3_minus_array(i, ll) = WP3_minus
      end do  ! End outer loop
#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc end parallel loop
      !$acc end data
#endif
      deallocate(Index_N)
   end subroutine RTA_minus_driver_GPU_using_Ind



 
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
             ! do ss=1,nptk    
                qprime=IJK(:,ii)
                qdprime=IJK(:,jj)
                ! ----------- fix -----------
                qtprime=modulo(q+qprime+qdprime,ngrid)
                ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                ! qtprime=IJK(:,ss)
                if (ss.ge.1) then
                   ! if (all(qtprime.eq.modulo(q+qprime+qdprime,ngrid))) then
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
                   ! end if
                end if
             ! end do
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
                ! do ss=jj,nptk
                   qprime=IJK(:,ii)
                   qdprime=IJK(:,jj)
                ! ----------- fix -----------
                   qtprime=modulo(q+qprime-qdprime,ngrid)
                   ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                   ! qtprime=IJK(:,ss)
                   if (ss.ge.jj) then
                      ! if (all(qtprime.eq.modulo(q+qprime-qdprime,ngrid))) then
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
                      ! end if
                   end if
                ! end do
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
                ! do ss=jj,nptk
                   qprime=IJK(:,ii)
                   qdprime=IJK(:,jj)
                ! ----------- fix -----------
                   qtprime=modulo(q-qprime-qdprime,ngrid)
                   ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                   ! qtprime=IJK(:,ss)
                   if (ss.ge.jj) then
                      ! if (all(qtprime.eq.modulo(q-qprime-qdprime,ngrid))) then
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
                      ! end if
                   end if
                ! end do
             end do
          end do
        end if
   end subroutine

   ! Compute the number of allowed four-phonon ++ +- -- processes and
   ! their contribution to phase space using sampling method.
    subroutine NP_plusplus_sample(mm,energy,velocity,Nlist,List,IJK,N_plusplus,P_plusplus)
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
 
       ! ----------- sampling method add -----------
       integer(kind=8) :: rand_num, iter, total_process, matrix_iter, temp_int
       real :: rand_matrix(num_sample_process_4ph_phase_space*5)
       total_process = Nbands*nptk*Nbands*nptk*Nbands
       ! ----------- end sampling method add -----------
 
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
          ! ----------- sampling method add -----------
          ! Note: 
          ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
          ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
          ! phonon3: jj is momentum [ii:nptk], k is energy [startbranch:Nband]
          ! phonon4: ss is momentum [1:nptk], l is energy [1:Nband]
          ! total combination: Nbands*nptk*Nbands*nptk*Nbands/2
          ! j = n + FLOOR((m+1-n)*u)  ! equation to generate random number j in [n,m] from a uniform distribution u in [0,1]
 
          call random_seed()
          call random_number(rand_matrix)
          iter=0
          do while (iter<=num_sample_process_4ph_phase_space)
             matrix_iter = (iter-1)*5
             ii=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+1)) ! phonon2 momentum [1:nptk]
             j=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+2)) ! phonon2 branch [1:Nband]
             jj=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+3)) ! phonon3 momentum [ii:nptk]
             startbranch=1
             if(jj==ii) startbranch=j
             k=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+4)) ! phonon3 branch [startbranch:Nband]
 
             qprime=IJK(:,ii)
             qdprime=IJK(:,jj)
             qtprime=modulo(q+qprime+qdprime,ngrid)
             ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! phonon4 momentum [1:nptk]
             l=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+5)) ! phonon4 branch [1:Nband]
             ! ----------- end sampling method add -----------
 
             if (jj>=ii .AND. k>=startbranch) then ! make sure that the phonon3 is in [ii:nptk], phonon3 branch is in [startbranch:Nband]
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
             end if ! check momentum and energy conservation
             iter = iter + 1
          end do ! iter
       end if
       ! ----------- sampling method add -----------
       N_plusplus = INT(REAL(N_plusplus, 8)/num_sample_process_4ph_phase_space*total_process, 8)    ! must first do the division and then do the multiplication. 
                                                                                              ! Otherwise it will overflow the max of int32
       P_plusplus = P_plusplus*total_process/num_sample_process_4ph_phase_space
       ! ----------- end sampling method add -----------
    end subroutine
 
   subroutine NP_plusminus_sample(mm,energy,velocity,Nlist,List,IJK,N_plusminus,P_plusminus)
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
 
       ! ----------- sampling method add -----------
       integer(kind=8) :: rand_num, iter, total_process, matrix_iter, temp_int
       real :: rand_matrix(num_sample_process_4ph_phase_space*5)
       total_process = Nbands*nptk*Nbands*nptk*Nbands
       ! ----------- end sampling method add -----------
 
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
       ! ----------- sampling method add -----------
       ! Note: 
       ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
       ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
       ! phonon3: jj is momentum [1:nptk], k is energy [1:Nband]
       ! phonon4: ss is momentum [jj:nptk], l is energy [startbranch:Nband]
       ! total combination: Nbands*nptk*Nbands*nptk*Nbands/2
       ! j = n + FLOOR((m+1-n)*u)  ! equation to generate random number j in [n,m] from a uniform distribution u in [0,1]
 
          call random_seed()
          call random_number(rand_matrix)
          iter=0
          do while (iter<=num_sample_process_4ph_phase_space)
             matrix_iter = (iter-1)*5
             ii=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+1)) ! phonon2 momentum [1:nptk]
             j=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+2)) ! phonon2 branch [1:Nband]
             jj=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+3)) ! phonon3 momentum [1:nptk]
             k=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+4)) ! phonon3 branch [1:Nband]
 
             qprime=IJK(:,ii)
             qdprime=IJK(:,jj)
             qtprime=modulo(q+qprime-qdprime,ngrid)
             ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! phonon4 momentum [1:nptk]
             ! ----------- end sampling method add -----------
             startbranch=1
             if(ss==jj) startbranch=k
             l=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+5)) ! phonon4 branch [startbranch:Nband]
             if (ss>=jj .AND. l>=startbranch) then ! make sure that the phonon4 is in [jj:nptk], phonon4 branch is in [startbranch:Nband]
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
             end if ! check momentum and energy conservation
             iter = iter + 1
          end do ! iter
        end if
       N_plusminus = INT(REAL(N_plusminus, 8)/num_sample_process_4ph_phase_space*total_process, 8)    ! must first do the division and then do the multiplication. 
                                                                                                ! Otherwise it will overflow the max of int32
       P_plusminus = P_plusminus*total_process/num_sample_process_4ph_phase_space
   end subroutine
 
   subroutine NP_minusminus_sample(mm,energy,velocity,Nlist,List,IJK,N_minusminus,P_minusminus)
       implicit none
 
       integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
       real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
       integer(kind=8),intent(out) :: N_minusminus
       real(kind=8),intent(out) :: P_minusminus
 
       integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
       integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
       integer(kind=4) :: ii,jj,kk,ll,ss,startbranch1,startbranch2 ! Ziqi add
       real(kind=8) :: sigma
       real(kind=8) :: omega,omegap,omegadp,omegatp
 
       ! ----------- sampling method add -----------
       integer(kind=8) :: rand_num, iter, total_process, matrix_iter, temp_int
       real :: rand_matrix(num_sample_process_4ph_phase_space*5)
       total_process = Nbands*nptk*Nbands*nptk*Nbands
       ! ----------- end sampling method add -----------
 
 
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
          ! ----------- sampling method add -----------
          ! Note: 
          ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
          ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
          ! phonon3: jj is momentum [ii:nptk], k is energy [startbranch1:Nband]
          ! phonon4: ss is momentum [jj:nptk], l is energy [startbranch2:Nband]
          ! total combination: Nbands*nptk*Nbands*nptk*Nbands/2
          ! j = n + FLOOR((m+1-n)*u)  ! equation to generate random number j in [n,m] from a uniform distribution u in [0,1]
 
          call random_seed()
          call random_number(rand_matrix)
          iter=0
          do while (iter<=num_sample_process_4ph_phase_space)
             matrix_iter = (iter-1)*5
             ii=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+1)) ! phonon2 momentum [1:nptk]
             j=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+2)) ! phonon2 branch [1:Nband]
             jj=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+3)) ! phonon3 momentum [ii:nptk]
             startbranch1=1
             if(jj==ii) startbranch1=j
             k=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+4)) ! phonon3 branch [startbranch1:Nband]
 
             qprime=IJK(:,ii)
             qdprime=IJK(:,jj)
             qtprime=modulo(q-qprime-qdprime,ngrid)
             ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! phonon4 momentum [jj:nptk]
             startbranch2=1
             if(ss==jj) startbranch2=k
             l=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+5)) ! phonon4 branch [startbranch2:Nband]
             ! ----------- end sampling method add -----------
 
             if (jj>=ii .AND. k>=startbranch1 .AND. ss>=jj .AND. l>=startbranch2) then ! make sure that the phonon3 is in [ii:nptk], phonon3 branch is in [startbranch:Nband]
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
             end if ! check momentum and energy conservation
             iter = iter + 1
          end do ! iter
        end if
 
       ! ----------- sampling method add -----------
       N_minusminus = INT(REAL(N_minusminus, 8)/num_sample_process_4ph_phase_space*total_process, 8)    ! must first do the division and then do the multiplication. 
                                                                                                  ! Otherwise it will overflow the max of int32
       P_minusminus = P_minusminus*total_process/num_sample_process_4ph_phase_space
       ! ----------- end sampling method add -----------
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
 
      integer(kind=4) :: mm,i,ll
      integer(kind=8) :: N_plusplus_reduce(Nlist*Nbands)
      integer(kind=8) :: N_plusminus_reduce(Nlist*Nbands)
      integer(kind=8) :: N_minusminus_reduce(Nlist*Nbands)
      real(kind=8) :: Pspace_plusplus_reduce(Nlist*Nbands)
      real(kind=8) :: Pspace_plusminus_reduce(Nlist*Nbands)
      real(kind=8) :: Pspace_minusminus_reduce(Nlist*Nbands)

      integer(kind=8) :: N_temp
      real(kind=8) :: Pspace_temp

      integer(kind=4) :: chunkstart, chunkend, chunkid


      Pspace_plusplus_total=0.d0
      Pspace_plusminus_total=0.d0
      Pspace_minusminus_total=0.d0
      Pspace_plusplus_reduce=0.d0
      Pspace_plusminus_reduce=0.d0
      Pspace_minusminus_reduce=0.d0

      N_plusplus=0
      N_plusminus=0
      N_minusminus=0
      N_plusplus_reduce=0
      N_plusminus_reduce=0
      N_minusminus_reduce=0

      N_temp=0
      Pspace_temp = 0.d0

      do chunkid=myid+1,numchunk,numprocs
         chunkstart = (chunkid-1)*chunksize+1
         chunkend = chunkid*chunksize
         if (chunkend.gt.Nbands*Nlist) chunkend = Nbands*Nlist
         !$OMP PARALLEL DO default(shared) schedule(dynamic,1) private(mm,i,ll,N_temp,Pspace_temp)
         do mm=chunkstart,chunkend
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1
            if (energy(List(int((mm-1)/Nbands)+1),modulo(mm-1,Nbands)+1).le.omega_max) then
               if (num_sample_process_4ph_phase_space==-1) then ! do not sample
                  call NP_plusplus(mm,energy,velocity,Nlist,List,IJK,N_temp,Pspace_temp)
                  N_plusplus_reduce(mm) = N_temp
                  Pspace_plusplus_reduce(mm) = Pspace_temp
                  call NP_plusminus(mm,energy,velocity,Nlist,List,IJK,N_temp,Pspace_temp)
                  N_plusminus_reduce(mm) = N_temp
                  Pspace_plusminus_reduce(mm) = Pspace_temp
                  call NP_minusminus(mm,energy,velocity,Nlist,List,IJK,N_temp,Pspace_temp)
                  N_minusminus_reduce(mm) = N_temp
                  Pspace_minusminus_reduce(mm) = Pspace_temp
               else ! do sampling method
                  call NP_plusplus_sample(mm,energy,velocity,Nlist,List,IJK,N_temp,Pspace_temp)
                  N_plusplus_reduce(mm) = N_temp
                  Pspace_plusplus_reduce(mm) = Pspace_temp
                  call NP_plusminus_sample(mm,energy,velocity,Nlist,List,IJK,N_temp,Pspace_temp)
                  N_plusminus_reduce(mm) = N_temp
                  Pspace_plusminus_reduce(mm) = Pspace_temp
                  call NP_minusminus_sample(mm,energy,velocity,Nlist,List,IJK,N_temp,Pspace_temp)
                  N_minusminus_reduce(mm) = N_temp
                  Pspace_minusminus_reduce(mm) = Pspace_temp
               endif
            endif
         end do
         !$OMP END PARALLEL DO
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



   ! Subroutines for 4ph calculations
   ! Compute the number of allowed four-phonon ++ +- -- processes and
   ! their contribution to phase space.
   subroutine Ind_only_plusplus(mm,energy,velocity,Nlist,List,IJK,N_plusplus,&
               Indof2ndPhonon_plusplus,Indof3rdPhonon_plusplus,Indof4thPhonon_plusplus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     integer(kind=8),intent(in) :: N_plusplus
     integer(kind=8),intent(out) :: Indof2ndPhonon_plusplus(N_plusplus),Indof3rdPhonon_plusplus(N_plusplus),Indof4thPhonon_plusplus(N_plusplus)
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
     real(kind=8) :: sigma
     real(kind=8) :: omega,omegap,omegadp,omegatp
     integer(kind=8) :: N_plusplus_count

 
     do ii=0,Ngrid(1)-1       ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
     end do
     N_plusplus_count=0
     i=modulo(mm-1,Nbands)+1
     ll=int((mm-1)/Nbands)+1
     q=IJK(:,list(ll))
     omega=energy(list(ll),i)
     if (omega.ne.0) then
       do ii=1,nptk
          do jj=ii,nptk
             ! do ss=1,nptk    
                qprime=IJK(:,ii)
                qdprime=IJK(:,jj)
                ! ----------- fix -----------
                qtprime=modulo(q+qprime+qdprime,ngrid)
                ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                ! qtprime=IJK(:,ss)
                if (ss.ge.1) then
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
                                    N_plusplus_count=N_plusplus_count+1
                                    if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                       Indof2ndPhonon_plusplus(N_plusplus_count)=(ii-1)*Nbands+j
                                       Indof3rdPhonon_plusplus(N_plusplus_count)=(jj-1)*Nbands+k
                                       Indof4thPhonon_plusplus(N_plusplus_count)=(ss-1)*Nbands+l
                                  end if
                               end if
                            end do
                         end do
                      end do
                end if
          end do
       end do
     end if
   end subroutine



   subroutine Ind_only_plusminus(mm,energy,velocity,Nlist,List,IJK,N_plusminus,&
            Indof2ndPhonon_plusminus,Indof3rdPhonon_plusminus,Indof4thPhonon_plusminus)
       implicit none
 
       integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
       real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
       integer(kind=8),intent(in) :: N_plusminus
     integer(kind=8) :: N_plusminus_count
       integer(kind=8),intent(out) :: Indof2ndPhonon_plusminus(N_plusminus),Indof3rdPhonon_plusminus(N_plusminus),Indof4thPhonon_plusminus(N_plusminus)

 
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
      N_plusminus_count=0
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       q=IJK(:,list(ll))
       omega=energy(list(ll),i)
       if (omega.ne.0) then
          do ii=1,nptk
             do jj=1,nptk
                ! do ss=jj,nptk
                   qprime=IJK(:,ii)
                   qdprime=IJK(:,jj)
                ! ----------- fix -----------
                   qtprime=modulo(q+qprime-qdprime,ngrid)
                   ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                   ! qtprime=IJK(:,ss)
                   if (ss.ge.jj) then
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
                                          N_plusminus_count=N_plusminus_count+1
                                           if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                             Indof2ndPhonon_plusminus(N_plusminus_count)=(ii-1)*Nbands+j
                                             Indof3rdPhonon_plusminus(N_plusminus_count)=(jj-1)*Nbands+k
                                             Indof4thPhonon_plusminus(N_plusminus_count)=(ss-1)*Nbands+l
                                        end if
                                     end if
                                  end if
                               end do
                            end do
                         end do
                      ! end if
                   end if
                ! end do
             end do
          end do
        end if
   end subroutine
 



   subroutine Ind_only_minusminus(mm,energy,velocity,Nlist,List,IJK,N_minusminus,&
            Indof2ndPhonon_minusminus,Indof3rdPhonon_minusminus,Indof4thPhonon_minusminus)
       implicit none
 
       integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk)
       real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
       integer(kind=8),intent(in) :: N_minusminus
       integer(kind=8),intent(out) :: Indof2ndPhonon_minusminus(N_minusminus),Indof3rdPhonon_minusminus(N_minusminus),Indof4thPhonon_minusminus(N_minusminus)
 
       integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
       integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
       integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
       real(kind=8) :: sigma
       real(kind=8) :: omega,omegap,omegadp,omegatp
 
     integer(kind=8) :: N_minusminus_count

       do ii=0,Ngrid(1)-1        ! G1 direction
       do jj=0,Ngrid(2)-1     ! G2 direction
          do kk=0,Ngrid(3)-1  ! G3 direction
             Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
          end do
       end do
       end do
       N_minusminus_count=0
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       q=IJK(:,list(ll))
       omega=energy(list(ll),i)
       if (omega.ne.0) then
          do ii=1,nptk
             do jj=ii,nptk
                ! do ss=jj,nptk
                   qprime=IJK(:,ii)
                   qdprime=IJK(:,jj)
                ! ----------- fix -----------
                   qtprime=modulo(q-qprime-qdprime,ngrid)
                   ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                   ! qtprime=IJK(:,ss)
                   if (ss.ge.jj) then
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
                                       N_minusminus_count=N_minusminus_count+1
                                        if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                          Indof2ndPhonon_minusminus(N_minusminus_count)=(ii-1)*Nbands+j
                                          Indof3rdPhonon_minusminus(N_minusminus_count)=(jj-1)*Nbands+k
                                          Indof4thPhonon_minusminus(N_minusminus_count)=(ss-1)*Nbands+l
                                     end if
                                  end if
                               end do
                            end do
                         end do
                      ! end if
                   end if
                ! end do
             end do
          end do
        end if
   end subroutine




   ! Wrapper around NP_plusplus,NP_plusminus and NP_minusminus that splits the work among processors.
   subroutine Ind_only_driver_4ph(energy,velocity,Nlist,List,IJK,&
      N_plusplus,N_plusminus,N_minusminus,&
      Indof2ndPhonon_plusplus,Indof3rdPhonon_plusplus,Indof4thPhonon_plusplus,&
      Indof2ndPhonon_plusminus,Indof3rdPhonon_plusminus,Indof4thPhonon_plusminus,&
      Indof2ndPhonon_minusminus,Indof3rdPhonon_minusminus,Indof4thPhonon_minusminus)
      implicit none

      include "mpif.h"

      real(kind=8),intent(in) :: energy(nptk,nbands)
      real(kind=8),intent(in) :: velocity(nptk,nbands,3)
      integer(kind=4),intent(in) :: NList
      integer(kind=4),intent(in) :: List(Nlist)
      integer(kind=4),intent(in) :: IJK(3,nptk)
      integer(kind=8),intent(in) :: N_plusplus(Nlist*Nbands)
      integer(kind=8),intent(in) :: N_plusminus(Nlist*Nbands)
      integer(kind=8),intent(in) :: N_minusminus(Nlist*Nbands)

      integer(kind=8),intent(inout) :: Indof2ndPhonon_plusplus(:),Indof3rdPhonon_plusplus(:),Indof4thPhonon_plusplus(:) 
      integer(kind=8),intent(inout) :: Indof2ndPhonon_plusminus(:),Indof3rdPhonon_plusminus(:),Indof4thPhonon_plusminus(:)
      integer(kind=8),intent(inout) :: Indof2ndPhonon_minusminus(:),Indof3rdPhonon_minusminus(:),Indof4thPhonon_minusminus(:)

      integer(kind=4) :: mm,i,ll

      !   temporary variables
      integer(kind=8),allocatable :: Indof2ndPhonon(:)
      integer(kind=8),allocatable :: Indof3rdPhonon(:)
      integer(kind=8),allocatable :: Indof4thPhonon(:)
      !   end temporary variables

      ! size of the matrix
      integer(kind=8) :: maxsize
      integer(kind=8) :: Ntotal_plusplus
      integer(kind=8) :: Ntotal_plusminus
      integer(kind=8) :: Ntotal_minusminus
      integer(kind=8) :: Naccum_plusplus(Nbands*Nlist)
      integer(kind=8) :: Naccum_plusminus(Nbands*Nlist)
      integer(kind=8) :: Naccum_minusminus(Nbands*Nlist)
      ! end size of the matrix

      integer(kind=4) :: chunkstart, chunkend, chunkid ! for parallel computing

      ! chunk size for parallelization allreduce
      integer(kind=8), parameter :: chunk_size_allreduce = 100000000
      integer(kind=8) :: offset_allreduce, remaining_allreduce, this_chunk_allreduce

      maxsize=max(maxval(N_plusplus),maxval(N_plusminus),maxval(N_minusminus))

      !   allocate temporary arrays
      allocate(Indof2ndPhonon(maxsize))
      allocate(Indof3rdPhonon(maxsize))
      allocate(Indof4thPhonon(maxsize))


      Naccum_plusplus(1)=0
      Naccum_plusminus(1)=0
      Naccum_minusminus(1)=0
      do mm=2,Nbands*Nlist
         Naccum_plusplus(mm)=Naccum_plusplus(mm-1)+N_plusplus(mm-1)
         Naccum_plusminus(mm)=Naccum_plusminus(mm-1)+N_plusminus(mm-1)
         Naccum_minusminus(mm)=Naccum_minusminus(mm-1)+N_minusminus(mm-1)
      end do
      Ntotal_plusplus=sum(N_plusplus)
      Ntotal_plusminus=sum(N_plusminus)
      Ntotal_minusminus=sum(N_minusminus)

      ! initialize the output arrays
      Indof2ndPhonon_plusplus=0
      Indof3rdPhonon_plusplus=0
      Indof4thPhonon_plusplus=0
      Indof2ndPhonon_plusminus=0
      Indof3rdPhonon_plusminus=0
      Indof4thPhonon_plusminus=0
      Indof2ndPhonon_minusminus=0
      Indof3rdPhonon_minusminus=0
      Indof4thPhonon_minusminus=0

      ! initialize for the temporary variables
      Indof2ndPhonon=0
      Indof3rdPhonon=0
      Indof4thPhonon=0

      do chunkid=myid+1,numchunk,numprocs
         chunkstart = (chunkid-1)*chunksize+1
         chunkend = chunkid*chunksize
         if (chunkend.gt.Nbands*Nlist) chunkend = Nbands*Nlist
         !$OMP PARALLEL DO default(shared) schedule(dynamic,1) private(mm,i,ll,Indof2ndPhonon,Indof3rdPhonon,Indof4thPhonon)
         do mm=chunkstart,chunkend
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1
            if (energy(List(int((mm-1)/Nbands)+1),modulo(mm-1,Nbands)+1).le.omega_max) then
               if(N_plusplus(mm).ne.0) then
                  call Ind_only_plusplus(mm,energy,velocity,Nlist,List,IJK,N_plusplus(mm),&
                              Indof2ndPhonon(1:N_plusplus(mm)),Indof3rdPhonon(1:N_plusplus(mm)),Indof4thPhonon(1:N_plusplus(mm)))
                     Indof2ndPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                        Indof2ndPhonon(1:N_plusplus(mm))
                     Indof3rdPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                        Indof3rdPhonon(1:N_plusplus(mm))
                     Indof4thPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                        Indof4thPhonon(1:N_plusplus(mm))
               end if
               if(N_plusminus(mm).ne.0) then
                  call Ind_only_plusminus(mm,energy,velocity,Nlist,List,IJK,N_plusminus(mm),&
                              Indof2ndPhonon(1:N_plusminus(mm)),Indof3rdPhonon(1:N_plusminus(mm)),Indof4thPhonon(1:N_plusminus(mm)))
                     Indof2ndPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                        Indof2ndPhonon(1:N_plusminus(mm))
                     Indof3rdPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                        Indof3rdPhonon(1:N_plusminus(mm))
                     Indof4thPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                        Indof4thPhonon(1:N_plusminus(mm))
               end if
               if(N_minusminus(mm).ne.0) then
                  call Ind_only_minusminus(mm,energy,velocity,Nlist,List,IJK,N_minusminus(mm),&
                              Indof2ndPhonon(1:N_minusminus(mm)),Indof3rdPhonon(1:N_minusminus(mm)),Indof4thPhonon(1:N_minusminus(mm)))
                     Indof2ndPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                        Indof2ndPhonon(1:N_minusminus(mm))
                     Indof3rdPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                        Indof3rdPhonon(1:N_minusminus(mm))
                     Indof4thPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                        Indof4thPhonon(1:N_minusminus(mm))
               end if
            endif
         end do
         !$OMP END PARALLEL DO
      end do

      ! ------- MPI_ALLREDUCE done in chunks -------
      offset_allreduce = 0
      remaining_allreduce = Ntotal_plusplus
      do while (remaining_allreduce > 0)
         this_chunk_allreduce = min(chunk_size_allreduce, remaining_allreduce)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof2ndPhonon_plusplus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof3rdPhonon_plusplus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof4thPhonon_plusplus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         offset_allreduce = offset_allreduce + this_chunk_allreduce
         remaining_allreduce = remaining_allreduce - this_chunk_allreduce
      end do

      offset_allreduce = 0
      remaining_allreduce = Ntotal_plusminus
      do while (remaining_allreduce > 0)
         this_chunk_allreduce = min(chunk_size_allreduce, remaining_allreduce)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof2ndPhonon_plusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof3rdPhonon_plusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof4thPhonon_plusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         offset_allreduce = offset_allreduce + this_chunk_allreduce
         remaining_allreduce = remaining_allreduce - this_chunk_allreduce
      end do

      offset_allreduce = 0
      remaining_allreduce = Ntotal_minusminus
      do while (remaining_allreduce > 0)
         this_chunk_allreduce = min(chunk_size_allreduce, remaining_allreduce)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof2ndPhonon_minusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof3rdPhonon_minusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof4thPhonon_minusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ll)
         offset_allreduce = offset_allreduce + this_chunk_allreduce
         remaining_allreduce = remaining_allreduce - this_chunk_allreduce
      end do
      ! ------- end MPI_ALLREDUCE -------

      !  deallocate temporary arrays
      deallocate(Indof4thPhonon)
      deallocate(Indof3rdPhonon)
      deallocate(Indof2ndPhonon)

   end subroutine Ind_only_driver_4ph




 
 
   ! Scattering amplitudes of 4ph ++ processes.
   subroutine Ind_plusplus(mm,N_plusplus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Indof2ndPhonon_plusplus,Indof3rdPhonon_plusplus,Indof4thPhonon_plusplus,Gamma_plusplus,WP4_plusplus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=8),intent(in) :: N_plusplus
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     integer(kind=8),intent(out) :: Indof2ndPhonon_plusplus(N_plusplus),Indof3rdPhonon_plusplus(N_plusplus),Indof4thPhonon_plusplus(N_plusplus)
     real(kind=8),intent(out) :: Gamma_plusplus(N_plusplus),WP4_plusplus
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
     integer(kind=8) :: N_plusplus_count
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
     real(kind=8) :: sigma
     real(kind=8) :: fBEprime,fBEdprime,fBEtprime
     real(kind=8) :: omega,omegap,omegadp,omegatp
     real(kind=8) :: realqprime(3),realqdprime(3),realqtprime(3)
     real(kind=8) :: WP4
     complex(kind=8) :: Vp
 
     do ii=0,Ngrid(1)-1        ! G1 direction
        do jj=0,Ngrid(2)-1     ! G2 direction
           do kk=0,Ngrid(3)-1  ! G3 direction
              Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
           end do
        end do
     end do
     N_plusplus_count=0
     WP4_plusplus=0.d0
     i=modulo(mm-1,Nbands)+1
     ll=int((mm-1)/Nbands)+1
     q=IJK(:,list(ll))
     omega=energy(list(ll),i)
     ! Loop over all processes, detecting those that are allowed and
     ! computing their amplitudes. Compared to three-phonon, the new algorithm 
     ! avoids double-counting when selecting allowed processes.
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
                if (all(qtprime.eq.modulo(q+qprime+qdprime,ngrid))) then ! Momentum conservation law
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
                               -velocity(jj,k,:)+velocity(ss,l,:)) ! Unit of sigma is rad/ps
                               if (abs(omega+omegap+omegadp-omegatp).le.sigma) then
                                  N_plusplus_count=N_plusplus_count+1
                                  ! Filtering those abnormal scattering rates (they actually strictly obey the energy conservation) 
                                  ! and set the delta function to be 1
                                  if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                  fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
                                  fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                                  fBEtprime=1.d0/(exp(hbar*omegatp/Kb/T)-1.D0)
                                  Indof2ndPhonon_plusplus(N_plusplus_count)=(ii-1)*Nbands+j
                                  Indof3rdPhonon_plusplus(N_plusplus_count)=(jj-1)*Nbands+k
                                  Indof4thPhonon_plusplus(N_plusplus_count)=(ss-1)*Nbands+l
                                  WP4=(fBEprime*fBEdprime*(1+fBEtprime)-(1+fBEprime)*(1+fBEdprime)*fBEtprime)*&
                                     exp(-(omega+omegap+omegadp-omegatp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                                     (omega*omegap*omegadp*omegatp)
                                  if(omegap.le.1.25 .or. omegadp.le.1.25 .or.omegatp.le.1.25) WP4=0 ! Very-low frequency phonons
                                  WP4_plusplus=WP4_plusplus+WP4
                                  Vp=Vp_pp(i,j,k,l,list(ll),ii,jj,ss,&
                                           realqprime,realqdprime,realqtprime,eigenvect,&
                                           Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
                                  Gamma_plusplus(N_plusplus_count)=hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                  ! At this point, units of Gamma are
                                  ! (1.d-34J*s)^2*(1.d12/s)^(-5)*1amu^(-4)*(ev/angstrom**4)^2,
                                  ! that is, 3.37617087*1.d9 THz
                                  Gamma_plusplus(N_plusplus_count)=Gamma_plusplus(N_plusplus_count)*3.37617087*1.d9/nptk/nptk
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
   end subroutine Ind_plusplus
 
   ! Scattering amplitudes of 4ph +-/-+ processes.
   subroutine Ind_plusminus(mm,N_plusminus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Indof2ndPhonon_plusminus,Indof3rdPhonon_plusminus,Indof4thPhonon_plusminus,Gamma_plusminus,WP4_plusminus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=8),intent(in) :: N_plusminus
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     integer(kind=8),intent(out) :: Indof2ndPhonon_plusminus(N_plusminus),Indof3rdPhonon_plusminus(N_plusminus),Indof4thPhonon_plusminus(N_plusminus)
     real(kind=8),intent(out) :: Gamma_plusminus(N_plusminus),WP4_plusminus
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
     integer(kind=8) :: N_plusminus_count
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
     real(kind=8) :: sigma
     real(kind=8) :: fBEprime,fBEdprime,fBEtprime
     real(kind=8) :: omega,omegap,omegadp,omegatp
     real(kind=8) :: realqprime(3),realqdprime(3),realqtprime(3)
     real(kind=8) :: WP4
     complex(kind=8) :: Vp
 
     do ii=0,Ngrid(1)-1        ! G1 direction
        do jj=0,Ngrid(2)-1     ! G2 direction
           do kk=0,Ngrid(3)-1  ! G3 direction
              Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
           end do
        end do
     end do
     N_plusminus_count=0
     WP4_plusminus=0.d0
     i=modulo(mm-1,Nbands)+1
     ll=int((mm-1)/Nbands)+1
     q=IJK(:,list(ll))
     omega=energy(list(ll),i)
     ! Loop over all processes, detecting those that are allowed and
     ! computing their amplitudes.
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
                                     N_plusminus_count=N_plusminus_count+1
                                     if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                     fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
                                     fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                                     fBEtprime=1.d0/(exp(hbar*omegatp/Kb/T)-1.D0)
                                     Indof2ndPhonon_plusminus(N_plusminus_count)=(ii-1)*Nbands+j
                                     Indof3rdPhonon_plusminus(N_plusminus_count)=(jj-1)*Nbands+k
                                     Indof4thPhonon_plusminus(N_plusminus_count)=(ss-1)*Nbands+l
                                     fBEtprime=1.d0/(exp(hbar*omegatp/Kb/T)-1.D0)
                                     WP4=(fBEprime*(1+fBEdprime)*(1+fBEtprime)-(1+fBEprime)*fBEdprime*fBEtprime)*&
                                        exp(-(omega+omegap-omegadp-omegatp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                                        (omega*omegap*omegadp*omegatp)
                                     if(omegap.le.1.25 .or. omegadp.le.1.25 .or.omegatp.le.1.25) WP4=0
                                     WP4_plusminus=WP4_plusminus+WP4
                                     Vp=Vp_pm(i,j,k,l,list(ll),ii,jj,ss,&
                                                 realqprime,realqdprime,realqtprime,eigenvect,&
                                                 Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
                                     Gamma_plusminus(N_plusminus_count)=hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                     Gamma_plusminus(N_plusminus_count)=Gamma_plusminus(N_plusminus_count)*3.37617087*1.d9/nptk/nptk
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
   end subroutine Ind_plusminus
 
   ! Scattering amplitudes of 4ph -- processes.
   subroutine Ind_minusminus(mm,N_minusminus,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Indof2ndPhonon_minusminus,Indof3rdPhonon_minusminus,Indof4thPhonon_minusminus,Gamma_minusminus,WP4_minusminus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=8), intent(in) :: N_minusminus
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     integer(kind=8),intent(out) :: Indof2ndPhonon_minusminus(N_minusminus),Indof3rdPhonon_minusminus(N_minusminus),Indof4thPhonon_minusminus(N_minusminus)
     real(kind=8),intent(out) :: Gamma_minusminus(N_minusminus),WP4_minusminus
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
     integer(kind=8) :: N_minusminus_count
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
     real(kind=8) :: sigma
     real(kind=8) :: fBEprime,fBEdprime,fBEtprime
     real(kind=8) :: omega,omegap,omegadp,omegatp
     real(kind=8) :: realqprime(3),realqdprime(3),realqtprime(3)
     real(kind=8) :: WP4
     complex(kind=8) :: Vp
 
     do ii=0,Ngrid(1)-1        ! G1 direction
        do jj=0,Ngrid(2)-1     ! G2 direction
           do kk=0,Ngrid(3)-1  ! G3 direction
              Index_N(ii,jj,kk)=(kk*Ngrid(2)+jj)*Ngrid(1)+ii+1
           end do
        end do
     end do
     N_minusminus_count=0
     WP4_minusminus=0.d0
     i=modulo(mm-1,Nbands)+1
     ll=int((mm-1)/Nbands)+1
     q=IJK(:,list(ll))
     omega=energy(list(ll),i)
     ! Loop over all processes, detecting those that are allowed and
     ! computing their amplitudes.
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
                                  N_minusminus_count=N_minusminus_count+1
                                  if (sigma.le.0.001) sigma=1/sqrt(Pi)
                                  fBEprime=1.d0/(exp(hbar*omegap/Kb/T)-1.D0)
                                  fBEdprime=1.d0/(exp(hbar*omegadp/Kb/T)-1.D0)
                                  fBEtprime=1.d0/(exp(hbar*omegatp/Kb/T)-1.D0)
                                  Indof2ndPhonon_minusminus(N_minusminus_count)=(ii-1)*Nbands+j
                                  Indof3rdPhonon_minusminus(N_minusminus_count)=(jj-1)*Nbands+k
                                  Indof4thPhonon_minusminus(N_minusminus_count)=(ss-1)*Nbands+l
                                  WP4=((1+fBEprime)*(1+fBEdprime)*(1+fBEtprime)-fBEprime*fBEdprime*fBEtprime)*&
                                     exp(-(omega-omegap-omegadp-omegatp)**2/(sigma**2))/sigma/sqrt(Pi)/&
                                     (omega*omegap*omegadp*omegatp)
                                  if(omegap.le.1.25 .or. omegadp.le.1.25 .or.omegatp.le.1.25) WP4=0
                                  WP4_minusminus=WP4_minusminus+WP4
                                  Vp=Vp_mm(i,j,k,l,list(ll),ii,jj,ss,&
                                           realqprime,realqdprime,realqtprime,eigenvect,&
                                           Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
                                  Gamma_minusminus(N_minusminus_count)=hbarp*hbarp*pi/8.d0*WP4*abs(Vp)**2
                                  Gamma_minusminus(N_minusminus_count)=Gamma_minusminus(N_minusminus_count)*3.37617087*1.d9/nptk/nptk
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
   end subroutine Ind_minusminus
 
   ! Wrapper around four-phonon Ind_plusplus, Ind_plusminus and Ind_minusminus that splits the work among processors.
   subroutine Ind_driver_4ph(energy,velocity,eigenvect,Nlist,List,IJK,N_plusplus,N_plusminus,N_minusminus,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,&
       Indof2ndPhonon_plusplus,Indof3rdPhonon_plusplus,Indof4thPhonon_plusplus,Gamma_plusplus,&
       Indof2ndPhonon_plusminus,Indof3rdPhonon_plusminus,Indof4thPhonon_plusminus,Gamma_plusminus,&
       Indof2ndPhonon_minusminus,Indof3rdPhonon_minusminus,Indof4thPhonon_minusminus,Gamma_minusminus,&
       rate_scatt_4ph,rate_scatt_plusplus,rate_scatt_plusminus,rate_scatt_minusminus, WP4_plusplus,WP4_plusminus,WP4_minusminus)
 
      implicit none
   
      include "mpif.h"
   
      real(kind=8),intent(in) :: energy(nptk,nbands)
      real(kind=8),intent(in) :: velocity(nptk,nbands,3)
      complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
      integer(kind=4),intent(in) :: NList
      integer(kind=4),intent(in) :: List(Nlist)
      integer(kind=4),intent(in) :: IJK(3,nptk)
      integer(kind=8),intent(in) :: N_plusplus(Nlist*Nbands)
      integer(kind=8),intent(in) :: N_plusminus(Nlist*Nbands)
      integer(kind=8),intent(in) :: N_minusminus(Nlist*Nbands)
      integer(kind=4),intent(in) :: Ntri
      real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri)
      real(kind=8),intent(in) :: R_j(3,Ntri)
      real(kind=8),intent(in) :: R_k(3,Ntri)
      real(kind=8),intent(in) :: R_l(3,Ntri)
      integer(kind=4),intent(in) :: Index_i(Ntri)
      integer(kind=4),intent(in) :: Index_j(Ntri)
      integer(kind=4),intent(in) :: Index_k(Ntri)
      integer(kind=4),intent(in) :: Index_l(Ntri)
      ! data for output 
      integer(kind=8),intent(out) :: Indof2ndPhonon_plusplus(:),Indof3rdPhonon_plusplus(:),Indof4thPhonon_plusplus(:) 
      integer(kind=8),intent(out) :: Indof2ndPhonon_plusminus(:),Indof3rdPhonon_plusminus(:),Indof4thPhonon_plusminus(:)
      integer(kind=8),intent(out) :: Indof2ndPhonon_minusminus(:),Indof3rdPhonon_minusminus(:),Indof4thPhonon_minusminus(:)
      real(kind=8),intent(out) :: Gamma_plusplus(:),Gamma_plusminus(:),Gamma_minusminus(:)
      !  end data for rate_scatt_4ph
      ! data for output, they are not allocatable
      real(kind=8),intent(out) :: rate_scatt_4ph(Nbands,Nlist)
      real(kind=8),intent(out) :: rate_scatt_plusplus(Nbands,Nlist),rate_scatt_plusminus(Nbands,Nlist),rate_scatt_minusminus(Nbands,Nlist)
      real(kind=8),intent(out) :: WP4_plusplus(Nbands,Nlist),WP4_plusminus(Nbands,Nlist),WP4_minusminus(Nbands,Nlist)
      !  end data for output
   
      integer(kind=4) :: i,nthread
      integer(kind=4) :: ll
      integer(kind=4) :: mm
      integer(kind=8) :: maxsize
      integer(kind=8) :: Ntotal_plusplus
      integer(kind=8) :: Ntotal_plusminus
      integer(kind=8) :: Ntotal_minusminus
      integer(kind=8) :: Naccum_plusplus(Nbands*Nlist)
      integer(kind=8) :: Naccum_plusminus(Nbands*Nlist)
      integer(kind=8) :: Naccum_minusminus(Nbands*Nlist)
      !   temporary variables
      integer(kind=8),allocatable :: Indof2ndPhonon(:)
      integer(kind=8),allocatable :: Indof3rdPhonon(:)
      integer(kind=8),allocatable :: Indof4thPhonon(:)
      real(kind=8),allocatable :: Gamma0(:)
      real(kind=8) :: WP4_temp ! variable for parallel computing
      !   end temporary variables

      ! reduced variables
      real(kind=8),allocatable  :: WP4_plusplus_reduce(:),WP4_plusminus_reduce(:),WP4_minusminus_reduce(:) ! change to allocatable for reduce matrix
      real(kind=8),allocatable :: rate_scatt_plusplus_reduce(:),rate_scatt_plusminus_reduce(:),rate_scatt_minusminus_reduce(:)
      ! end reduced variables
   
      integer(kind=4) :: chunkstart, chunkend, chunkid ! for parallel computing
      integer(kind=4) :: ierr

      ! chunk size for parallelization allreduce
      integer(kind=8), parameter :: chunk_size_allreduce = 100000000
      integer(kind=8) :: offset_allreduce, remaining_allreduce, this_chunk_allreduce

   
      maxsize=max(maxval(N_plusplus),maxval(N_plusminus),maxval(N_minusminus))

      !   allocate temporary arrays
      allocate(Indof2ndPhonon(maxsize))
      allocate(Indof3rdPhonon(maxsize))
      allocate(Indof4thPhonon(maxsize))
      allocate(Gamma0(maxsize))
      !   end allocate temporary arrays
   
      Naccum_plusplus(1)=0
      Naccum_plusminus(1)=0
      Naccum_minusminus(1)=0
      do mm=2,Nbands*Nlist
         Naccum_plusplus(mm)=Naccum_plusplus(mm-1)+N_plusplus(mm-1)
         Naccum_plusminus(mm)=Naccum_plusminus(mm-1)+N_plusminus(mm-1)
         Naccum_minusminus(mm)=Naccum_minusminus(mm-1)+N_minusminus(mm-1)
      end do
      Ntotal_plusplus=sum(N_plusplus)
      Ntotal_plusminus=sum(N_plusminus)
      Ntotal_minusminus=sum(N_minusminus)


      !  allocate reduced output arrays
      allocate(WP4_plusplus_reduce(Nbands*Nlist))
      allocate(WP4_plusminus_reduce(Nbands*Nlist))
      allocate(WP4_minusminus_reduce(Nbands*Nlist))
      allocate(rate_scatt_plusplus_reduce(Nbands*Nlist))
      allocate(rate_scatt_plusminus_reduce(Nbands*Nlist))
      allocate(rate_scatt_minusminus_reduce(Nbands*Nlist))
      !  end allocate reduced output arrays

      ! initialize the output arrays
      Indof2ndPhonon_plusplus=0
      Indof3rdPhonon_plusplus=0
      Indof4thPhonon_plusplus=0
      Indof2ndPhonon_plusminus=0
      Indof3rdPhonon_plusminus=0
      Indof4thPhonon_plusminus=0
      Indof2ndPhonon_minusminus=0
      Indof3rdPhonon_minusminus=0
      Indof4thPhonon_minusminus=0
      Gamma_plusplus=0.d0
      Gamma_plusminus=0.d0
      Gamma_minusminus=0.d0
      rate_scatt_plusplus=0.d0
      rate_scatt_plusminus=0.d0
      rate_scatt_minusminus=0.d0
      WP4_plusplus=0.d0
      WP4_plusminus=0.d0
      WP4_minusminus=0.d0
      !   end initialize the output arrays

      ! initialize the reduced output arrays
      WP4_plusplus_reduce=0.d0
      WP4_plusminus_reduce=0.d0
      WP4_minusminus_reduce=0.d0
      rate_scatt_plusplus_reduce=0.d0
      rate_scatt_plusminus_reduce=0.d0
      rate_scatt_minusminus_reduce=0.d0
      !   end initialize the reduced output arrays

      ! initialize for the temporary variables
      Indof2ndPhonon=0
      Indof3rdPhonon=0
      Indof4thPhonon=0
      Gamma0=0.d0
      WP4_temp=0.d0
      !   end initialize for the temporary variables
            
      do chunkid=myid+1,numchunk,numprocs
         chunkstart = (chunkid-1)*chunksize+1
         chunkend = chunkid*chunksize
         if (chunkend.gt.Nbands*Nlist) chunkend = Nbands*Nlist
         !$OMP PARALLEL DO default(shared) schedule(dynamic,1) private(mm,i,ll,Gamma0,Indof2ndPhonon,Indof3rdPhonon,Indof4thPhonon,WP4_temp)
         do mm=chunkstart,chunkend
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1
            if(N_plusplus(mm).ne.0) then
               call Ind_plusplus(mm,N_plusplus(mm),energy,velocity,eigenvect,Nlist,List,&
                     Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                     Indof2ndPhonon(1:N_plusplus(mm)),Indof3rdPhonon(1:N_plusplus(mm)),Indof4thPhonon(1:N_plusplus(mm)),&
                     Gamma0(1:N_plusplus(mm)),WP4_temp)
                  WP4_plusplus_reduce(mm)=WP4_temp
                  Indof2ndPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                     Indof2ndPhonon(1:N_plusplus(mm))
                  Indof3rdPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                     Indof3rdPhonon(1:N_plusplus(mm))
                  Indof4thPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                     Indof4thPhonon(1:N_plusplus(mm))
                  Gamma_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                     Gamma0(1:N_plusplus(mm))
                  rate_scatt_plusplus_reduce(mm)=sum(Gamma0(1:N_plusplus(mm))) ! change to mm since we allocate the array in 1d
            end if
            if(N_plusminus(mm).ne.0) then
               call Ind_plusminus(mm,N_plusminus(mm),energy,velocity,eigenvect,Nlist,List,&
                     Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                     Indof2ndPhonon(1:N_plusminus(mm)),Indof3rdPhonon(1:N_plusminus(mm)),Indof4thPhonon(1:N_plusminus(mm)),&
                     Gamma0(1:N_plusminus(mm)),WP4_temp)
                  WP4_plusminus_reduce(mm)=WP4_temp
                  Indof2ndPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                     Indof2ndPhonon(1:N_plusminus(mm))
                  Indof3rdPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                     Indof3rdPhonon(1:N_plusminus(mm))
                  Indof4thPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                     Indof4thPhonon(1:N_plusminus(mm))
                  Gamma_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                     Gamma0(1:N_plusminus(mm))
                  rate_scatt_plusminus_reduce(mm)=sum(Gamma0(1:N_plusminus(mm)))
            end if
            if(N_minusminus(mm).ne.0) then
               call Ind_minusminus(mm,N_minusminus(mm),energy,velocity,eigenvect,Nlist,List,&
                     Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                     Indof2ndPhonon(1:N_minusminus(mm)),Indof3rdPhonon(1:N_minusminus(mm)),Indof4thPhonon(1:N_minusminus(mm)),&
                     Gamma0(1:N_minusminus(mm)),WP4_temp)
                  WP4_minusminus_reduce(mm)=WP4_temp
                  Indof2ndPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                     Indof2ndPhonon(1:N_minusminus(mm))
                  Indof3rdPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                     Indof3rdPhonon(1:N_minusminus(mm))
                  Indof4thPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                     Indof4thPhonon(1:N_minusminus(mm))
                  Gamma_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                     Gamma0(1:N_minusminus(mm))
                  rate_scatt_minusminus_reduce(mm)=sum(Gamma0(1:N_minusminus(mm)))
            end if
         end do
         !$OMP END PARALLEL DO
      end do

      ! ------- MPI_ALLREDUCE done in chunks -------
      offset_allreduce = 0
      remaining_allreduce = Ntotal_plusplus
      do while (remaining_allreduce > 0)
         this_chunk_allreduce = min(chunk_size_allreduce, remaining_allreduce)

         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof2ndPhonon_plusplus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof3rdPhonon_plusplus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof4thPhonon_plusplus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Gamma_plusplus(offset_allreduce+1), this_chunk_allreduce, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

         offset_allreduce = offset_allreduce + this_chunk_allreduce
         remaining_allreduce = remaining_allreduce - this_chunk_allreduce
      end do

      offset_allreduce = 0
      remaining_allreduce = Ntotal_plusminus
      do while (remaining_allreduce > 0)
         this_chunk_allreduce = min(chunk_size_allreduce, remaining_allreduce)

         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof2ndPhonon_plusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof3rdPhonon_plusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof4thPhonon_plusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Gamma_plusminus(offset_allreduce+1), this_chunk_allreduce, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

         offset_allreduce = offset_allreduce + this_chunk_allreduce
         remaining_allreduce = remaining_allreduce - this_chunk_allreduce
      end do

      offset_allreduce = 0
      remaining_allreduce = Ntotal_minusminus
      do while (remaining_allreduce > 0)
         this_chunk_allreduce = min(chunk_size_allreduce, remaining_allreduce)

         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof2ndPhonon_minusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof3rdPhonon_minusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Indof4thPhonon_minusminus(offset_allreduce+1), this_chunk_allreduce, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, Gamma_minusminus(offset_allreduce+1), this_chunk_allreduce, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

         offset_allreduce = offset_allreduce + this_chunk_allreduce
         remaining_allreduce = remaining_allreduce - this_chunk_allreduce
      end do

      call MPI_ALLREDUCE(WP4_plusplus_reduce,WP4_plusplus,Nbands*Nlist,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rate_scatt_plusplus_reduce,rate_scatt_plusplus,Nbands*Nlist,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(WP4_plusminus_reduce,WP4_plusminus,Nbands*Nlist,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rate_scatt_plusminus_reduce,rate_scatt_plusminus,Nbands*Nlist,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(WP4_minusminus_reduce,WP4_minusminus,Nbands*Nlist,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rate_scatt_minusminus_reduce,rate_scatt_minusminus,Nbands*Nlist,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      ! ------- end MPI_ALLREDUCE -------

      !  deallocate temporary arrays
      deallocate(Gamma0)
      deallocate(Indof4thPhonon)
      deallocate(Indof3rdPhonon)
      deallocate(Indof2ndPhonon)
      !  end deallocate temporary arrays
      ! deallocate reduced output arrays
      deallocate(WP4_plusplus_reduce)
      deallocate(WP4_plusminus_reduce)
      deallocate(WP4_minusminus_reduce)
      deallocate(rate_scatt_plusplus_reduce)
      deallocate(rate_scatt_plusminus_reduce)
      deallocate(rate_scatt_minusminus_reduce)
      ! end deallocate reduced output arrays

   
      rate_scatt_4ph=rate_scatt_plusplus+rate_scatt_plusminus+rate_scatt_minusminus

   
   end subroutine 
 
 
 
 
   ! RTA-only version of 4ph ++ process, Ind_plusplus; Avoid double counting
   subroutine RTA_plusplus(mm,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Gamma_plusplus,WP4_plusplus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_plusplus,WP4_plusplus
 
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
             ! do ss=1,nptk
                qprime=IJK(:,ii)
                realqprime=matmul(rlattvec,qprime/dble(ngrid))
                qdprime=IJK(:,jj)
                realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
                ! ----------- fix -----------
                qtprime=modulo(q+qprime+qdprime,ngrid)
                realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
                ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                ! qtprime=IJK(:,ss)
                ! realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
                if (ss.ge.1) then
                   ! if (all(qtprime.eq.modulo(q+qprime+qdprime,ngrid))) then
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
                                     end if
                                  end if
                               end if
                            end do
                         end do
                      end do
                   ! end if
                end if ! greater equal than 1
             ! end do
          end do
       end do
       WP4_plusplus=WP4_plusplus/nptk/nptk
     end if
     ! converting to THz
     Gamma_plusplus=Gamma_plusplus*3.37617087*1.d9/nptk/nptk
 
   end subroutine
 
   ! RTA-only version of 4ph +- process, Ind_plusminus
   subroutine RTA_plusminus(mm,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Gamma_plusminus,WP4_plusminus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_plusminus,WP4_plusminus
 
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
             ! do ss=jj,nptk
                qprime=IJK(:,ii)
                realqprime=matmul(rlattvec,qprime/dble(ngrid))
                qdprime=IJK(:,jj)
                realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
                ! ----------- fix -----------
                qtprime=modulo(q+qprime-qdprime,ngrid)
                realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
                ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                ! qtprime=IJK(:,ss)
                ! realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
                if (ss.ge.jj) then
                   ! if (all(qtprime.eq.modulo(q+qprime-qdprime,ngrid))) then
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
                                        end if
                                     end if
                                  end if
                               end if
                            end do
                         end do
                      end do
                   ! end if
                end if ! greater equal than jj
             ! end do
          end do
       end do
       WP4_plusminus=WP4_plusminus/nptk/nptk
     end if
     
     ! converting to THz
     Gamma_plusminus=Gamma_plusminus*3.37617087*1.d9/nptk/nptk
 
   end subroutine
 
   ! RTA-only version of 4ph -- process, Ind_minusminus
   subroutine RTA_minusminus(mm,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Gamma_minusminus,WP4_minusminus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_minusminus,WP4_minusminus
 
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
             ! do ss=jj,nptk
                qprime=IJK(:,ii)
                realqprime=matmul(rlattvec,qprime/dble(ngrid))
                qdprime=IJK(:,jj)
                realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
                ! ----------- fix -----------
                qtprime=modulo(q-qprime-qdprime,ngrid)
                realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
                ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! Instead of iteration, calculate momentem ss from other phonons can have higher efficiency.
                ! ----------- end fix -----------
                ! qtprime=IJK(:,ss)
                ! realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
                if (ss.ge.jj) then
                   ! if (all(qtprime.eq.modulo(q-qprime-qdprime,ngrid))) then
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
                                     end if
                                  end if
                               end if
                            end do
                         end do
                      end do
                   ! end if
                end if ! greater than jj
             ! end do
          end do
       end do
       WP4_minusminus=WP4_minusminus/nptk/nptk
     end if
     
     ! converting to THz
     Gamma_minusminus=Gamma_minusminus*3.37617087*1.d9/nptk/nptk
   
      end subroutine


   subroutine RTA_plusplus_driver_GPU_using_Ind(energy,velocity,eigenvect,Nlist,List, &
      Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK, &
      N_plusplus, &
      Indof2ndPhonon_plusplus,Indof3rdPhonon_plusplus,Indof4thPhonon_plusplus,&
      rate_scatt_plusplus,WP4_plusplus_array)
      implicit none
      include "mpif.h"

      ! Input variables
      integer(kind=4), intent(in) :: Nlist, List(Nlist), IJK(3,nptk)
      integer(kind=4), intent(in) :: Ntri
      integer(kind=4), intent(in) :: Index_i(Ntri), Index_j(Ntri), Index_k(Ntri), Index_l(Ntri)
      real(kind=8), intent(in)    :: energy(nptk,Nbands), velocity(nptk,Nbands,3)
      complex(kind=8), intent(in) :: eigenvect(nptk,Nbands,Nbands)
      real(kind=8), intent(in)    :: Phi(3,3,3,3,Ntri), R_j(3,Ntri), R_k(3,Ntri), R_l(3,Ntri)

      ! Precomputed scattering process data:
      integer(kind=8),intent(in) :: N_plusplus(Nlist*Nbands) 
      integer(kind=8),intent(in) :: Indof2ndPhonon_plusplus(:),Indof3rdPhonon_plusplus(:),Indof4thPhonon_plusplus(:)

      ! Output variables
      real(kind=8), intent(out)   :: rate_scatt_plusplus(Nbands,Nlist)
      real(kind=8), intent(out)   :: WP4_plusplus_array(Nbands,Nlist)

      ! Local variables
      integer(kind=4) :: mm, i, ll
      real(kind=8)  :: Gamma_plusplus, WP4_temp, WP4_val
      real(kind=8)  :: omega, omegap, omegadp, omegatp, sigma
      real(kind=8)  :: fBEprime, fBEdprime, fBEtprime
      integer(kind=4) :: q(3), qprime(3), qdprime(3), qtprime(3)
      real(kind=8)  :: realq(3), realqprime(3), realqdprime(3), realqtprime(3)
      integer(kind=4) :: ii, jj, j, k, l, ss
      integer(kind=8) :: Naccum_plusplus(Nbands*Nlist)
      integer(kind=8) indnow

      ! For base_sigma calculation
      real(kind=8)  :: v_base_sigma(3), base_sigma_value

      ! For inline Vp_pp calculation        
      integer(kind=4) :: ll_Vp, rr_Vp,ss_Vp,tt_Vp,uu_Vp
      complex(kind=8) :: prefactor, Vp0, Vp_pp_inline

      ! Array for mapping q-vectors to a single index
      integer(kind=4), allocatable :: Index_N(:,:,:)
      integer(kind=4) :: d1, d2, d3

      ! Build the Index_N mapping array over the 3D grid.
      d1 = Ngrid(1)
      d2 = Ngrid(2)
      d3 = Ngrid(3)
      allocate(Index_N(0:d1-1, 0:d2-1, 0:d3-1))
      do ii = 0, d1-1
         do jj = 0, d2-1
            do l = 0, d3-1
               Index_N(ii,jj,l) = (l*d2 + jj)*d1 + ii + 1
            end do
         end do
      end do

      ! Build an accumulation array so that for each (band, list) pair (indexed by mm)
      ! we know where its allowed process data begins in the index arrays.
      Naccum_plusplus(1) = 0
      do mm = 2, Nbands*Nlist
         Naccum_plusplus(mm) = Naccum_plusplus(mm-1) + N_plusplus(mm-1)
      end do

      ! Initialize the output arrays
      rate_scatt_plusplus = 0.d0
      WP4_plusplus_array   = 0.d0

#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc data copyin(energy,velocity,IJK,List,Nlist,Index_N,Ntri,Phi,R_j,R_k,R_l) &
      !$acc&  copyin(Index_i,Index_j,Index_k,Index_l)  &
      !$acc&  copyin(eigenvect)  &
      !$acc&  copyin(T,scalebroad,Ngrid,nptk,Nbands,rlattvec,masses,types,omega_max,onlyharmonic)  &
      !$acc&  copyin(N_plusplus,Naccum_plusplus)  &  ! add these index
      !$acc&  copyin(Indof2ndPhonon_plusplus, Indof3rdPhonon_plusplus, Indof4thPhonon_plusplus)  &
      !$acc copy(rate_scatt_plusplus, WP4_plusplus_array) 
      !$acc parallel loop gang default(present) &
      !$acc& num_gangs(NUM_GANGS) &
      !$acc& vector_length(VEC_LEN) &
      !$acc private(Gamma_plusplus, WP4_temp, i, ll)
#endif
      do mm = 1, Nbands*Nlist
         Gamma_plusplus=0.d00
         WP4_temp=0.d00
         i=modulo(mm-1,Nbands)+1
         ll=int((mm-1)/Nbands)+1

#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
      !$acc data copyin(energy,velocity,IJK,List,Nlist,Index_N,Ntri,Phi,R_j,R_k,R_l) &
      !$acc&  copyin(i,ll)  &
      !$acc&  copyin(Index_i,Index_j,Index_k,Index_l)  &
      !$acc&  copyin(eigenvect)  &
      !$acc&  copyin(T,scalebroad,Ngrid,nptk,Nbands,rlattvec,masses,types,omega_max,onlyharmonic)  &
      !$acc&  copyin(N_plusplus,Naccum_plusplus)  &  ! add these index
      !$acc& copyin(Indof2ndPhonon_plusplus(Naccum_plusplus(mm)+1:Naccum_plusplus(mm) + N_plusplus(mm))) &
      !$acc& copyin(Indof3rdPhonon_plusplus(Naccum_plusplus(mm)+1:Naccum_plusplus(mm) + N_plusplus(mm))) &
      !$acc& copyin(Indof4thPhonon_plusplus(Naccum_plusplus(mm)+1:Naccum_plusplus(mm) + N_plusplus(mm))) &
      !$acc copy(rate_scatt_plusplus, WP4_plusplus_array)
      !$acc parallel loop gang default(present) reduction(+:Gamma_plusplus,WP4_temp) &
      !$acc& num_gangs(NUM_GANGS) &
      !$acc& vector_length(VEC_LEN) &
      !$acc private(Gamma_plusplus, WP4_temp) &
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc loop vector reduction(+:Gamma_plusplus,WP4_temp) &
#endif
      !$acc private(ii,j,jj,k,ss,l,q,realq,qprime,realqprime,qdprime,realqdprime,qtprime,realqtprime, &
      !$acc         omega,omegap,omegadp,omegatp,sigma,fBEprime,fBEdprime,fBEtprime,WP4_val, &
      !$acc         v_base_sigma,base_sigma_value)
         do indnow = Naccum_plusplus(mm) + 1, Naccum_plusplus(mm) + N_plusplus(mm)

            ! Recover the indices for this scattering process:
            ii = (Indof2ndPhonon_plusplus(indnow)-1)/Nbands + 1
            j  = mod(Indof2ndPhonon_plusplus(indnow)-1, Nbands) + 1
            jj = (Indof3rdPhonon_plusplus(indnow)-1)/Nbands + 1
            k  = mod(Indof3rdPhonon_plusplus(indnow)-1, Nbands) + 1
            ss = (Indof4thPhonon_plusplus(indnow)-1)/Nbands + 1
            l  = mod(Indof4thPhonon_plusplus(indnow)-1, Nbands) + 1

            q=IJK(:,list(ll))
            realq=matmul(rlattvec,q/dble(ngrid))
            omega=energy(list(ll),i)

            qprime     = IJK(:,ii)
            realqprime = matmul(rlattvec, qprime/dble(ngrid))

            qdprime    = IJK(:,jj)
            realqdprime= matmul(rlattvec, qdprime/dble(ngrid))

            qtprime    = modulo(q + qprime + qdprime, ngrid)
            realqtprime= matmul(rlattvec, qtprime/dble(ngrid))

            omegap  = energy(ii,j)
            omegadp = energy(jj,k)
            omegatp = energy(ss,l)


            ! ---------- Sigma calculation inline ----------
            v_base_sigma = -velocity(jj,k,:) + velocity(ss,l,:)
            base_sigma_value = max( (dot_product(rlattvec(:,1),v_base_sigma)/ngrid(1))**2, &
                                    (dot_product(rlattvec(:,2),v_base_sigma)/ngrid(2))**2, &
                                    (dot_product(rlattvec(:,3),v_base_sigma)/ngrid(3))**2 )
            base_sigma_value = sqrt(base_sigma_value/2.)
            sigma = scalebroad * base_sigma_value
            if (sigma <= 0.001d0) sigma = 1.d0/sqrt(Pi)
            ! ---------- End Sigma calculation inline ----------

            fBEprime  = 1.d0/(exp(hbar*omegap / (Kb*T)) - 1.d0)
            fBEdprime = 1.d0/(exp(hbar*omegadp/(Kb*T)) - 1.d0)
            fBEtprime = 1.d0/(exp(hbar*omegatp/(Kb*T)) - 1.d0)

            WP4_val = (fBEprime*fBEdprime*(1.d0+fBEtprime) - &
                     (1.d0+fBEprime)*(1.d0+fBEdprime)*fBEtprime) * &
                     exp(- (omega+omegap+omegadp-omegatp)**2/(sigma**2)) / &
                     ( sigma * sqrt(Pi) * (omega*omegap*omegadp*omegatp) )

            if ((omegap <= 1.25d0) .or. (omegadp <= 1.25d0) .or. (omegatp <= 1.25d0)) then
               WP4_val = 0.d0
            end if
            WP4_temp = WP4_temp + WP4_val

            if (.not. onlyharmonic) then
               ! ---------- Vp_pp calculation inline ----------
               Vp_pp_inline=0.d0
#ifdef GPU_ALL_MODE_PARALLELIZATION
               !$acc loop reduction(+:Vp_pp_inline) &
#endif
#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
               !$acc loop vector reduction(+:Vp_pp_inline) &
#endif
               !$acc& private(Vp0,prefactor,rr_Vp,ss_Vp,tt_Vp,uu_Vp)
               do ll_Vp=1,Ntri
                  prefactor = 1.d0 / sqrt(masses(types(Index_i(ll_Vp))) * &
                                 masses(types(Index_j(ll_Vp))) * &
                                 masses(types(Index_k(ll_Vp))) * &
                                 masses(types(Index_l(ll_Vp))) ) * &
                     cmplx(cos(dot_product(realqprime, R_j(:, ll_Vp))), sin(dot_product(realqprime, R_j(:, ll_Vp))), kind=8) * &
                     cmplx(cos(dot_product(realqdprime, R_k(:, ll_Vp))), sin(dot_product(realqdprime, R_k(:, ll_Vp))), kind=8) * &
                     cmplx(cos(-dot_product(realqtprime, R_l(:, ll_Vp))), sin(-dot_product(realqtprime, R_l(:, ll_Vp))), kind=8)

                  Vp0=0.d0
                  do rr_Vp=1,3
                     do ss_Vp=1,3
                        do tt_Vp=1,3
                           do uu_Vp=1,3
                              Vp0=Vp0+Phi(uu_Vp,tt_Vp,ss_Vp,rr_Vp,ll_Vp)*&
                                 eigenvect(list(ll),i,uu_Vp+3*(Index_i(ll_Vp)-1))*&
                                 eigenvect(ii,j,tt_Vp+3*(Index_j(ll_Vp)-1))*&
                                 eigenvect(jj,k,ss_Vp+3*(Index_k(ll_Vp)-1))*&
                                 conjg(eigenvect(ss,l,rr_Vp+3*(Index_l(ll_Vp)-1)))
                           end do
                        end do
                     end do
                  end do
                  Vp_pp_inline=Vp_pp_inline+prefactor*Vp0
               end do
               ! ---------- End Vp_pp calculation inline ----------
               Gamma_plusplus = Gamma_plusplus + (hbarp*hbarp*pi/8.d0) * WP4_val * abs(Vp_pp_inline)**2
            end if
         end do ! indnow

#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
         !$acc end parallel loop
         !$acc end data
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
         !$acc end loop
#endif
         WP4_temp = WP4_temp / (nptk * nptk)
         rate_scatt_plusplus(i,ll) = Gamma_plusplus* 3.37617087d9 / (nptk*nptk)
         WP4_plusplus_array(i,ll)  = WP4_temp
      end do
#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc end parallel loop
      !$acc end data
#endif
      deallocate(Index_N)
   end subroutine RTA_plusplus_driver_GPU_using_Ind


   subroutine RTA_plusminus_driver_GPU_using_Ind(energy,velocity,eigenvect,Nlist,List, &
      Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK, &
      N_plusminus, &
      Indof2ndPhonon_plusminus,Indof3rdPhonon_plusminus,Indof4thPhonon_plusminus, &
      rate_scatt_plusminus, WP4_plusminus_array)
      implicit none
      include "mpif.h"

      ! Input variables
      integer(kind=4), intent(in)    :: Nlist, List(Nlist), IJK(3,nptk)
      integer(kind=4), intent(in)    :: Ntri
      integer(kind=4), intent(in)    :: Index_i(Ntri), Index_j(Ntri), Index_k(Ntri), Index_l(Ntri)
      real(kind=8), intent(in)       :: energy(nptk,Nbands), velocity(nptk,Nbands,3)
      complex(kind=8), intent(in)    :: eigenvect(nptk,Nbands,Nbands)
      real(kind=8), intent(in)       :: Phi(3,3,3,3,Ntri), R_j(3,Ntri), R_k(3,Ntri), R_l(3,Ntri)

      ! Precomputed scattering process data:
      integer(kind=8), intent(in)    :: N_plusminus(Nbands*Nlist)
      integer(kind=8), intent(in)    :: Indof2ndPhonon_plusminus(:), Indof3rdPhonon_plusminus(:), Indof4thPhonon_plusminus(:)

      ! Output variables
      real(kind=8), intent(out)      :: rate_scatt_plusminus(Nbands,Nlist)
      real(kind=8), intent(out)      :: WP4_plusminus_array(Nbands,Nlist)

      ! Local variables
      integer(kind=4) :: mm, i, ll
      real(kind=8)  :: Gamma_plusminus, WP4_temp, WP4_val
      real(kind=8)  :: omega, omegap, omegadp, omegatp, sigma
      real(kind=8)  :: fBEprime, fBEdprime, fBEtprime
      integer(kind=4) :: q(3), qprime(3), qdprime(3), qtprime(3)
      real(kind=8)  :: realq(3), realqprime(3), realqdprime(3), realqtprime(3)
      integer(kind=4) :: j, k, l, ii, jj, ss
      integer(kind=8) :: Naccum_plusminus(Nbands*Nlist)
      integer(kind=8) :: indnow
      
      ! For base_sigma calculation
      real(kind=8)  :: v_base_sigma(3), base_sigma_value
      
      ! For inline Vp_pm calculation
      integer(kind=4) :: ll_Vp, rr_Vp, ss_Vp, tt_Vp, uu_Vp
      complex(kind=8) :: prefactor, Vp0, Vp_pm_inline

      ! Array for mapping q-vectors to a single index
      integer(kind=4), allocatable :: Index_N(:,:,:)
      integer(kind=4) :: d1, d2, d3

      ! Build the Index_N mapping array over the 3D grid.
      d1 = Ngrid(1)
      d2 = Ngrid(2)
      d3 = Ngrid(3)
      allocate(Index_N(0:d1-1, 0:d2-1, 0:d3-1))
      do ii = 0, d1-1
         do jj = 0, d2-1
            do l = 0, d3-1
               Index_N(ii,jj,l) = (l*d2 + jj)*d1 + ii + 1
            end do
         end do
      end do


      ! Build an accumulation array so that for each (band, list) pair (indexed by mm)
      ! we know where its allowed process data begins in the index arrays.
      Naccum_plusminus(1) = 0
      do mm = 2, Nbands*Nlist
         Naccum_plusminus(mm) = Naccum_plusminus(mm-1) + N_plusminus(mm-1)
      end do

      ! Initialize the output arrays
      rate_scatt_plusminus = 0.d0
      WP4_plusminus_array  = 0.d0

#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc data copyin(energy,velocity,IJK,List,Nlist,Index_N,Ntri,Phi,R_j,R_k,R_l) &
      !$acc&  copyin(Index_i,Index_j,Index_k,Index_l)  &
      !$acc&  copyin(eigenvect)  &
      !$acc&  copyin(T,scalebroad,Ngrid,nptk,Nbands,rlattvec,masses,types,omega_max,onlyharmonic)  &
      !$acc&  copyin(N_plusminus,Naccum_plusminus)  &  ! add these index
      !$acc&  copyin(Indof2ndPhonon_plusminus, Indof3rdPhonon_plusminus, Indof4thPhonon_plusminus)  &
      !$acc copy(rate_scatt_plusminus, WP4_plusminus_array) 
      !$acc parallel loop gang default(present) &
      !$acc& vector_length(VEC_LEN) &
      !$acc& num_gangs(NUM_GANGS) &
      !$acc private(Gamma_plusminus, WP4_temp, i, ll)
#endif
      do mm = 1, Nbands*Nlist
         Gamma_plusminus = 0.d0
         WP4_temp        = 0.d0
         i  = modulo(mm-1, Nbands) + 1
         ll = int((mm-1)/Nbands) + 1


#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
      !$acc data copyin(energy,velocity,IJK,List,Nlist,Index_N,Ntri,Phi,R_j,R_k,R_l) &
      !$acc&  copyin(i,ll)  &
      !$acc&  copyin(Index_i,Index_j,Index_k,Index_l)  &
      !$acc&  copyin(eigenvect)  &
      !$acc&  copyin(T,scalebroad,Ngrid,nptk,Nbands,rlattvec,masses,types,omega_max,onlyharmonic)  &
      !$acc&  copyin(N_plusminus,Naccum_plusminus)  &  ! add these index
      !$acc& copyin(Indof2ndPhonon_plusminus(Naccum_plusminus(mm)+1:Naccum_plusminus(mm) + N_plusminus(mm))) &
      !$acc& copyin(Indof3rdPhonon_plusminus(Naccum_plusminus(mm)+1:Naccum_plusminus(mm) + N_plusminus(mm))) &
      !$acc& copyin(Indof4thPhonon_plusminus(Naccum_plusminus(mm)+1:Naccum_plusminus(mm) + N_plusminus(mm))) &
      !$acc copy(rate_scatt_plusminus, WP4_plusminus_array) 
      !$acc parallel loop gang default(present) reduction(+:Gamma_plusminus,WP4_temp) &
      !$acc& vector_length(VEC_LEN) &
      !$acc& num_gangs(NUM_GANGS) &
      !$acc private(Gamma_plusminus, WP4_temp) &
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
         !$acc loop vector reduction(+:Gamma_plusminus,WP4_temp) &
#endif
         !$acc private(ii,j,jj,k,ss,l,q,realq,qprime,realqprime,qdprime,realqdprime,qtprime,realqtprime, &
         !$acc         omega,omegap,omegadp,omegatp,sigma,fBEprime,fBEdprime,fBEtprime,WP4_val, &
         !$acc         v_base_sigma,base_sigma_value)
         do indnow = Naccum_plusminus(mm)+1, Naccum_plusminus(mm) + N_plusminus(mm)

            ! Recover the indices for this scattering process:
            ii = (Indof2ndPhonon_plusminus(indnow) - 1) / Nbands + 1
            j  = mod(Indof2ndPhonon_plusminus(indnow) - 1, Nbands) + 1
            jj = (Indof3rdPhonon_plusminus(indnow) - 1) / Nbands + 1
            k  = mod(Indof3rdPhonon_plusminus(indnow) - 1, Nbands) + 1
            ss = (Indof4thPhonon_plusminus(indnow) - 1) / Nbands + 1
            l  = mod(Indof4thPhonon_plusminus(indnow) - 1, Nbands) + 1

            q = IJK(:, List(ll))
            realq = matmul(rlattvec, q / dble(Ngrid))
            omega = energy(List(ll), i)

            qprime = IJK(:, ii)
            realqprime = matmul(rlattvec, qprime / dble(Ngrid))

            qdprime = IJK(:, jj)
            realqdprime = matmul(rlattvec, qdprime / dble(Ngrid))

            qtprime = modulo(q + qprime - qdprime, Ngrid)
            realqtprime = matmul(rlattvec, qtprime / dble(Ngrid))

            omegap  = energy(ii, j)
            omegadp = energy(jj, k)
            omegatp = energy(ss, l)

            ! ---------- Sigma calculation inline ----------
            v_base_sigma = velocity(jj, k, :) - velocity(ss, l, :)
            base_sigma_value = max( (dot_product(rlattvec(:,1), v_base_sigma)/Ngrid(1))**2, &
                                    (dot_product(rlattvec(:,2), v_base_sigma)/Ngrid(2))**2, &
                                    (dot_product(rlattvec(:,3), v_base_sigma)/Ngrid(3))**2 )
            base_sigma_value = sqrt(base_sigma_value/2.d0)
            sigma = scalebroad * base_sigma_value
            if (sigma <= 0.001d0) sigma = 1.d0/sqrt(Pi)
            ! ---------- End Sigma calculation inline ----------

            fBEprime  = 1.d0 / (exp(hbar*omegap/(Kb*T)) - 1.d0)
            fBEdprime = 1.d0 / (exp(hbar*omegadp/(Kb*T)) - 1.d0)
            fBEtprime = 1.d0 / (exp(hbar*omegatp/(Kb*T)) - 1.d0)

            WP4_val = ( fBEprime*(1.d0+fBEdprime)*(1.d0+fBEtprime) - &
                        (1.d0+fBEprime)*fBEdprime*fBEtprime ) * &
                        exp( - (omega+omegap-omegadp-omegatp)**2 / (sigma**2) ) / &
                        ( sigma * sqrt(Pi) * (omega*omegap*omegadp*omegatp) )

            if ((omegap <= 1.25d0) .or. (omegadp <= 1.25d0) .or. (omegatp <= 1.25d0)) then
               WP4_val = 0.d0
            end if

            WP4_temp = WP4_temp + WP4_val

            if (.not. onlyharmonic) then
               ! ---------- Vp_pm calculation inline ----------
               Vp_pm_inline = (0.d0, 0.d0)
#ifdef GPU_ALL_MODE_PARALLELIZATION
               !$acc loop reduction(+:Vp_pm_inline) &
#endif
#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
               !$acc loop vector reduction(+:Vp_pm_inline) &
#endif
               !$acc& private(Vp0,prefactor,rr_Vp,ss_Vp,tt_Vp,uu_Vp)
               do ll_Vp = 1, Ntri
                  prefactor = 1.d0 / sqrt( masses(types(Index_i(ll_Vp))) * &
                                             masses(types(Index_j(ll_Vp))) * &
                                             masses(types(Index_k(ll_Vp))) * &
                                             masses(types(Index_l(ll_Vp))) ) * &
                              cmplx(cos(dot_product(realqprime, R_j(:,ll_Vp))), &
                                    sin(dot_product(realqprime, R_j(:,ll_Vp))), kind=8) * &
                              cmplx(cos(-dot_product(realqdprime, R_k(:,ll_Vp))), &
                                    sin(-dot_product(realqdprime, R_k(:,ll_Vp))), kind=8) * &
                              cmplx(cos(-dot_product(realqtprime, R_l(:,ll_Vp))), &
                                    sin(-dot_product(realqtprime, R_l(:,ll_Vp))), kind=8)
                  Vp0 = (0.d0, 0.d0)
                  do rr_Vp = 1, 3
                     do ss_Vp = 1, 3
                        do tt_Vp = 1, 3
                           do uu_Vp = 1, 3
                              Vp0 = Vp0 + Phi(uu_Vp,tt_Vp,ss_Vp,rr_Vp,ll_Vp) * &
                                          eigenvect(List(ll), i, uu_Vp+3*(Index_i(ll_Vp)-1)) * &
                                          eigenvect(ii, j, tt_Vp+3*(Index_j(ll_Vp)-1)) * &
                                          conjg(eigenvect(jj, k, ss_Vp+3*(Index_k(ll_Vp)-1))) * &
                                          conjg(eigenvect(ss, l, rr_Vp+3*(Index_l(ll_Vp)-1)))
                           end do
                        end do
                     end do
                  end do
                  Vp_pm_inline = Vp_pm_inline + prefactor * Vp0
               end do
               ! ---------- End Vp_pm calculation inline ----------
               Gamma_plusminus = Gamma_plusminus + (hbarp*hbarp*pi/8.d0) * WP4_val * abs(Vp_pm_inline)**2
            end if
         end do  ! indnow

#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
         !$acc end parallel loop
         !$acc end data
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
         !$acc end loop
#endif
         WP4_temp = WP4_temp / (nptk*nptk)
         rate_scatt_plusminus(i,ll) = Gamma_plusminus * 3.37617087d9 / (nptk*nptk)
         WP4_plusminus_array(i,ll)  = WP4_temp
      end do  ! End loop over mm
#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc end parallel loop
      !$acc end data
#endif
      deallocate(Index_N)
   end subroutine RTA_plusminus_driver_GPU_using_Ind



   subroutine RTA_minusminus_driver_GPU_using_Ind(energy,velocity,eigenvect,Nlist,List, &
      Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK, &
      N_minusminus, &
      Indof2ndPhonon_minusminus,Indof3rdPhonon_minusminus,Indof4thPhonon_minusminus, &
      rate_scatt_minusminus, WP4_minusminus_array)
      implicit none
      include "mpif.h"

      ! Input variables
      integer(kind=4), intent(in) :: Nlist, List(Nlist), IJK(3,nptk)
      integer(kind=4), intent(in) :: Ntri
      integer(kind=4), intent(in) :: Index_i(Ntri), Index_j(Ntri), Index_k(Ntri), Index_l(Ntri)
      real(kind=8), intent(in)    :: energy(nptk,Nbands), velocity(nptk,Nbands,3)
      complex(kind=8), intent(in) :: eigenvect(nptk,Nbands,Nbands)
      real(kind=8), intent(in)    :: Phi(3,3,3,3,Ntri), R_j(3,Ntri), R_k(3,Ntri), R_l(3,Ntri)

      ! Precomputed scattering process data:
      integer(kind=8), intent(in) :: N_minusminus(Nbands*Nlist)
      integer(kind=8), intent(in) :: Indof2ndPhonon_minusminus(:), Indof3rdPhonon_minusminus(:), Indof4thPhonon_minusminus(:)

      ! Output variables
      real(kind=8), intent(out)   :: rate_scatt_minusminus(Nbands,Nlist)
      real(kind=8), intent(out)   :: WP4_minusminus_array(Nbands,Nlist)

      ! Local variables
      integer(kind=4) :: mm, i, ll
      real(kind=8)  :: Gamma_minusminus, WP4_temp, WP4_val
      real(kind=8)  :: omega, omegap, omegadp, omegatp, sigma
      real(kind=8)  :: fBEprime, fBEdprime, fBEtprime
      integer(kind=4) :: q(3), qprime(3), qdprime(3), qtprime(3)
      real(kind=8)  :: realq(3), realqprime(3), realqdprime(3), realqtprime(3)
      integer(kind=4) :: j, k, l, ii, jj, ss
      integer(kind=8) :: Naccum_minusminus(Nbands*Nlist)
      integer(kind=8) :: indnow

      ! For base_sigma calculation
      real(kind=8)  :: v_base_sigma(3), base_sigma_value
      
      ! For inline Vp calculation
      integer(kind=4) :: ll_Vp, rr_Vp, ss_Vp, tt_Vp, uu_Vp
      complex(kind=8) :: prefactor, Vp0, Vp_mm_inline

      ! Array for mapping q-vectors to a single index
      integer(kind=4), allocatable :: Index_N(:,:,:)
      integer(kind=4) :: d1, d2, d3

      ! Build the Index_N mapping array over the 3D grid.
      d1 = Ngrid(1)
      d2 = Ngrid(2)
      d3 = Ngrid(3)
      allocate(Index_N(0:d1-1, 0:d2-1, 0:d3-1))
      do ii = 0, d1-1
         do jj = 0, d2-1
            do l = 0, d3-1
               Index_N(ii,jj,l) = (l*d2 + jj)*d1 + ii + 1
            end do
         end do
      end do

      ! Build an accumulation array so that for each (band, list) pair (indexed by mm)
      ! we know where its allowed process data begins in the index arrays.
      Naccum_minusminus(1) = 0
      do mm = 2, Nbands*Nlist
         Naccum_minusminus(mm) = Naccum_minusminus(mm-1) + N_minusminus(mm-1)
      end do

      ! Initialize the output arrays
      rate_scatt_minusminus = 0.d0
      WP4_minusminus_array  = 0.d0

#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc data copyin(energy,velocity,IJK,List,Nlist,Index_N,Ntri,Phi,R_j,R_k,R_l) &
      !$acc&  copyin(Index_i,Index_j,Index_k,Index_l)  &
      !$acc&  copyin(eigenvect)  &
      !$acc&  copyin(T,scalebroad,Ngrid,nptk,Nbands,rlattvec,masses,types,omega_max,onlyharmonic)  &
      !$acc&  copyin(N_minusminus,Naccum_minusminus)  &  ! add these index
      !$acc&  copyin(Indof2ndPhonon_minusminus, Indof3rdPhonon_minusminus, Indof4thPhonon_minusminus)  &
      !$acc copy(rate_scatt_minusminus, WP4_minusminus_array) 
      !$acc parallel loop gang default(present) &
      !$acc& vector_length(VEC_LEN) &
      !$acc& num_gangs(NUM_GANGS) &
      !$acc private(Gamma_minusminus, WP4_temp, i, ll)
#endif
      do mm = 1, Nbands*Nlist
         Gamma_minusminus = 0.d0
         WP4_temp         = 0.d0
         i  = modulo(mm-1, Nbands) + 1
         ll = int((mm-1)/Nbands) + 1


#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
         !$acc data copyin(energy,velocity,IJK,List,Nlist,Index_N,Ntri,Phi,R_j,R_k,R_l) &
         !$acc&  copyin(i,ll)  &
         !$acc&  copyin(Index_i,Index_j,Index_k,Index_l)  &
         !$acc&  copyin(eigenvect)  &
         !$acc&  copyin(T,scalebroad,Ngrid,nptk,Nbands,rlattvec,masses,types,omega_max,onlyharmonic)  &
         !$acc&  copyin(N_minusminus,Naccum_minusminus)  &  ! add these index
         !$acc& copyin(Indof2ndPhonon_minusminus(Naccum_minusminus(mm)+1:Naccum_minusminus(mm) + N_minusminus(mm))) &
         !$acc& copyin(Indof3rdPhonon_minusminus(Naccum_minusminus(mm)+1:Naccum_minusminus(mm) + N_minusminus(mm))) &
         !$acc& copyin(Indof4thPhonon_minusminus(Naccum_minusminus(mm)+1:Naccum_minusminus(mm) + N_minusminus(mm))) &
         !$acc copy(rate_scatt_minusminus, WP4_minusminus_array)
         !$acc parallel loop gang default(present) reduction(+:Gamma_minusminus,WP4_temp) &
         !$acc& vector_length(VEC_LEN) &
         !$acc& num_gangs(NUM_GANGS) &
         !$acc private(Gamma_minusminus, WP4_temp) &
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
         !$acc loop vector reduction(+:Gamma_minusminus,WP4_temp) &
#endif
         !$acc private(ii,j,jj,k,ss,l,q,realq,qprime,realqprime,qdprime,realqdprime,qtprime,realqtprime, &
         !$acc         omega,omegap,omegadp,omegatp,sigma,fBEprime,fBEdprime,fBEtprime,WP4_val, &
         !$acc         v_base_sigma,base_sigma_value)


         do indnow = Naccum_minusminus(mm)+1, Naccum_minusminus(mm) + N_minusminus(mm)
            ! Recover the indices for this scattering process:
            ii = (Indof2ndPhonon_minusminus(indnow)-1)/Nbands + 1
            j  = mod(Indof2ndPhonon_minusminus(indnow)-1, Nbands) + 1
            jj = (Indof3rdPhonon_minusminus(indnow)-1)/Nbands + 1
            k  = mod(Indof3rdPhonon_minusminus(indnow)-1, Nbands) + 1
            ss = (Indof4thPhonon_minusminus(indnow)-1)/Nbands + 1
            l  = mod(Indof4thPhonon_minusminus(indnow)-1, Nbands) + 1

            q = IJK(:, List(ll))
            realq = matmul(rlattvec, q/dble(Ngrid))
            omega = energy(List(ll), i)

            qprime = IJK(:, ii)
            realqprime = matmul(rlattvec, qprime/dble(Ngrid))

            qdprime = IJK(:, jj)
            realqdprime = matmul(rlattvec, qdprime/dble(Ngrid))

            qtprime = modulo(q - qprime - qdprime, Ngrid)
            realqtprime = matmul(rlattvec, qtprime/dble(Ngrid))

            omegap  = energy(ii, j)
            omegadp = energy(jj, k)
            omegatp = energy(ss, l)


            ! ---------- Sigma calculation inline ----------
            v_base_sigma = velocity(jj,k,:) - velocity(ss,l,:)
            base_sigma_value = max( (dot_product(rlattvec(:,1), v_base_sigma)/Ngrid(1))**2, &
                                    (dot_product(rlattvec(:,2), v_base_sigma)/Ngrid(2))**2, &
                                    (dot_product(rlattvec(:,3), v_base_sigma)/Ngrid(3))**2 )
            base_sigma_value = sqrt(base_sigma_value/2.d0)
            sigma = scalebroad * base_sigma_value
            if (sigma <= 1.d-3) sigma = 1.d0/sqrt(Pi)
            ! ---------- End Sigma calculation inline ----------

            fBEprime  = 1.d0/(exp(hbar*omegap/(Kb*T))-1.d0)
            fBEdprime = 1.d0/(exp(hbar*omegadp/(Kb*T))-1.d0)
            fBEtprime = 1.d0/(exp(hbar*omegatp/(Kb*T))-1.d0)

            WP4_val = ((1.d0+fBEprime)*(1.d0+fBEdprime)*(1.d0+fBEtprime) - &
                        fBEprime*fBEdprime*fBEtprime) * &
                     exp( - (omega - omegap - omegadp - omegatp)**2/(sigma**2) ) / &
                     ( sigma*sqrt(Pi)*(omega*omegap*omegadp*omegatp) )
            if ((omegap <= 1.25d0) .or. (omegadp <= 1.25d0) .or. (omegatp <= 1.25d0)) then
               WP4_val = 0.d0
            end if
            WP4_temp = WP4_temp + WP4_val

            if (.not. onlyharmonic) then
               ! ---------- Vp_mm calculation inline ----------
               Vp_mm_inline = (0.d0, 0.d0)
#ifdef GPU_ALL_MODE_PARALLELIZATION
               !$acc loop reduction(+:Vp_mm_inline) &
#endif
#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
               !$acc loop vector reduction(+:Vp_mm_inline) &
#endif
               !$acc& private(Vp0,prefactor,rr_Vp,ss_Vp,tt_Vp,uu_Vp)
               do ll_Vp = 1, Ntri
                  prefactor = 1.d0 / sqrt( masses(types(Index_i(ll_Vp))) * &
                                             masses(types(Index_j(ll_Vp))) * &
                                             masses(types(Index_k(ll_Vp))) * &
                                             masses(types(Index_l(ll_Vp))) ) * &
                              cmplx(cos(-dot_product(realqprime, R_j(:,ll_Vp))), &
                                    sin(-dot_product(realqprime, R_j(:,ll_Vp))), kind=8) * &
                              cmplx(cos(-dot_product(realqdprime, R_k(:,ll_Vp))), &
                                    sin(-dot_product(realqdprime, R_k(:,ll_Vp))), kind=8) * &
                              cmplx(cos(-dot_product(realqtprime, R_l(:,ll_Vp))), &
                                    sin(-dot_product(realqtprime, R_l(:,ll_Vp))), kind=8)
                  Vp0 = (0.d0, 0.d0)
                  do rr_Vp = 1, 3
                     do ss_Vp = 1, 3
                        do tt_Vp = 1, 3
                           do uu_Vp = 1, 3
                              Vp0 = Vp0 + Phi(uu_Vp,tt_Vp,ss_Vp,rr_Vp,ll_Vp) * &
                                       eigenvect(List(ll), i, uu_Vp+3*(Index_i(ll_Vp)-1)) * &
                                       conjg(eigenvect(ii, j, tt_Vp+3*(Index_j(ll_Vp)-1))) * &
                                       conjg(eigenvect(jj, k, ss_Vp+3*(Index_k(ll_Vp)-1))) * &
                                       conjg(eigenvect(ss, l, rr_Vp+3*(Index_l(ll_Vp)-1)))
                           end do
                        end do
                     end do
                  end do
                  Vp_mm_inline = Vp_mm_inline + prefactor * Vp0
               end do
               ! ---------- End Vp_mm calculation inline ----------
               Gamma_minusminus = Gamma_minusminus + (hbarp*hbarp*pi/8.d0)*WP4_val*abs(Vp_mm_inline)**2
            end if
         end do  ! End inner loop reduction
#ifdef GPU_MODE_BY_MODE_PARALLELIZATION
      !$acc end parallel loop
      !$acc end data
#endif
#ifdef GPU_ALL_MODE_PARALLELIZATION
      !$acc end loop
#endif
         WP4_temp = WP4_temp/(nptk*nptk)
         rate_scatt_minusminus(i,ll) = Gamma_minusminus * 3.37617087d9/(nptk*nptk)
         WP4_minusminus_array(i,ll)  = WP4_temp
      end do  ! End outer loop
#ifdef GPU_ALL_MODE_PARALLELIZATION
   !$acc end parallel loop
   !$acc end data
#endif
      deallocate(Index_N)
   end subroutine RTA_minusminus_driver_GPU_using_Ind







   ! RTA-only version of 4ph ++ process using sampling method, Ind_plusplus; Avoid double counting
   subroutine RTA_plusplus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Gamma_plusplus,WP4_plusplus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_plusplus,WP4_plusplus
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
     real(kind=8) :: sigma
     real(kind=8) :: fBEprime,fBEdprime,fBEtprime
     real(kind=8) :: omega,omegap,omegadp,omegatp
     real(kind=8) :: realq(3),realqprime(3),realqdprime(3),realqtprime(3)
     real(kind=8) :: WP4
     complex(kind=8) :: Vp
 
    ! ----------- sampling method add -----------
    integer(kind=8) :: rand_num, iter, total_process, matrix_iter, temp_int
    real :: rand_matrix(num_sample_process_4ph*5)
    total_process = Nbands*nptk*Nbands*nptk*Nbands
    ! ----------- end sampling method add -----------
 
 
     Gamma_plusplus=0.d00
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
       ! ----------- sampling method add -----------
       ! Note: 
       ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
       ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
       ! phonon3: jj is momentum [ii:nptk], k is energy [startbranch:Nband]
       ! phonon4: ss is momentum [1:nptk], l is energy [1:Nband]
       ! total combination: Nbands*nptk*Nbands*nptk*Nbands/2
       ! j = n + FLOOR((m+1-n)*u)  ! equation to generate random number j in [n,m] from a uniform distribution u in [0,1]
 
       call random_seed()
       call random_number(rand_matrix)
       iter=0
       do while (iter<=num_sample_process_4ph)
          matrix_iter = (iter-1)*5
          ii=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+1)) ! phonon2 momentum [1:nptk]
          j=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+2)) ! phonon2 branch [1:Nband]
          jj=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+3)) ! phonon3 momentum [ii:nptk]
          startbranch=1
          if(jj==ii) startbranch=j
          k=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+4)) ! phonon3 branch [startbranch:Nband]
 
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,qprime/dble(ngrid))
 
          qdprime=IJK(:,jj)
          realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
 
          qtprime=modulo(q+qprime+qdprime,ngrid)
          realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
          ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! phonon4 momentum [1:nptk]
          l=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+5)) ! phonon4 branch [1:Nband]
          ! ----------- end sampling method add -----------
          if (jj>=ii .AND. k>=startbranch) then ! make sure that the phonon3 is in [ii:nptk], phonon3 branch is in [startbranch:Nband]
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
                   end if
                end if
             end if
          end if ! check momentum and energy conservation
          iter = iter + 1
       end do ! iter                        
       WP4_plusplus=WP4_plusplus/nptk/nptk
     end if
     ! converting to THz
     Gamma_plusplus=Gamma_plusplus*3.37617087*1.d9/nptk/nptk
 
    ! ----------- sampling method add -----------
    Gamma_plusplus = Gamma_plusplus/num_sample_process_4ph*total_process  
    WP4_plusplus = WP4_plusplus*total_process/num_sample_process_4ph
    ! ----------- end sampling method add -----------
 
 
   end subroutine
 
   ! RTA-only version of 4ph +- process using sampling method, Ind_plusminus
   subroutine RTA_plusminus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Gamma_plusminus,WP4_plusminus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_plusminus,WP4_plusminus
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss,startbranch
     real(kind=8) :: sigma
     real(kind=8) :: fBEprime,fBEdprime,fBEtprime
     real(kind=8) :: omega,omegap,omegadp,omegatp
     real(kind=8) :: realq(3),realqprime(3),realqdprime(3),realqtprime(3)
     real(kind=8) :: WP4
     complex(kind=8) :: Vp
 
    ! ----------- sampling method add -----------
    integer(kind=8) :: rand_num, iter, total_process, matrix_iter, temp_int
    real :: rand_matrix(num_sample_process_4ph*5)
    total_process = Nbands*nptk*Nbands*nptk*Nbands
    ! ----------- end sampling method add -----------
 
 
     Gamma_plusminus=0.d00
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
       ! ----------- sampling method add -----------
       ! Note: 
       ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
       ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
       ! phonon3: jj is momentum [1:nptk], k is energy [1:Nband]
       ! phonon4: ss is momentum [jj:nptk], l is energy [startbranch:Nband]
       ! total combination: Nbands*nptk*Nbands*nptk*Nbands/2
       ! j = n + FLOOR((m+1-n)*u)  ! equation to generate random number j in [n,m] from a uniform distribution u in [0,1]
 
       call random_seed()
       call random_number(rand_matrix)
       iter=0
       do while (iter<=num_sample_process_4ph)
          matrix_iter = (iter-1)*5
          ii=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+1)) ! phonon2 momentum [1:nptk]
          j=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+2)) ! phonon2 branch [1:Nband]
          jj=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+3)) ! phonon3 momentum [1:nptk]
          k=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+4)) ! phonon3 branch [1:Nband]
 
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,qprime/dble(ngrid))
          qdprime=IJK(:,jj)
          realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
          qtprime=modulo(q+qprime-qdprime,ngrid)
          ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! phonon4 momentum [1:nptk]
          realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
          ! ----------- end sampling method add -----------
          startbranch=1
          if(ss==jj) startbranch=k
          l=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+5)) ! phonon4 branch [startbranch:Nband]
          if (ss>=jj .AND. l>=startbranch) then ! make sure that the phonon4 is in [jj:nptk], phonon4 branch is in [startbranch:Nband]
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
                      end if
                   end if
                end if
             end if
          end if ! check momentum and energy conservation
          iter = iter + 1
       end do ! iter
       WP4_plusminus=WP4_plusminus/nptk/nptk
     end if
     
     ! converting to THz
     Gamma_plusminus=Gamma_plusminus*3.37617087*1.d9/nptk/nptk
 
    ! ----------- sampling method add -----------
    Gamma_plusminus = Gamma_plusminus/num_sample_process_4ph*total_process    
    WP4_plusminus = WP4_plusminus*total_process/num_sample_process_4ph
    ! ----------- end sampling method add -----------
 
 
   end subroutine
 
   ! RTA-only version of 4ph -- process using sampling method, Ind_minusminus
   subroutine RTA_minusminus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
       Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
       Gamma_minusminus,WP4_minusminus)
     implicit none
 
     integer(kind=4),intent(in) :: mm,NList,List(Nlist),IJK(3,nptk),Ntri
     integer(kind=4),intent(in) :: Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri)
     real(kind=8),intent(in) :: energy(nptk,Nbands),velocity(nptk,Nbands,3)
     real(kind=8),intent(in) :: Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri)
     complex(kind=8),intent(in) :: eigenvect(nptk,Nbands,Nbands)
     real(kind=8),intent(out) :: Gamma_minusminus,WP4_minusminus
 
     integer(kind=4) :: q(3),qprime(3),qdprime(3),qtprime(3),i,j,k,l
     integer(kind=4) :: Index_N(0:(Ngrid(1)-1),0:(Ngrid(2)-1),0:(Ngrid(3)-1))
     integer(kind=4) :: ii,jj,kk,ll,ss,startbranch1,startbranch2 ! Ziqi add
     real(kind=8) :: sigma
     real(kind=8) :: fBEprime,fBEdprime,fBEtprime
     real(kind=8) :: omega,omegap,omegadp,omegatp
     real(kind=8) :: realq(3),realqprime(3),realqdprime(3),realqtprime(3)
     real(kind=8) :: WP4
     complex(kind=8) :: Vp
 
    ! ----------- sampling method add -----------
    integer(kind=8) :: rand_num, iter, total_process, matrix_iter, temp_int
    real :: rand_matrix(num_sample_process_4ph*5)
    total_process = Nbands*nptk*Nbands*nptk*Nbands
    ! ----------- end sampling method add -----------
 
 
     Gamma_minusminus=0.d00
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
       ! ----------- sampling method add -----------
       ! Note: 
       ! phonon1: ll is momentum [1:nptk], i is energy [1:Nband], ll and i are calculated from mm, mm is phonon number [1, Nbands*NList]
       ! phonon2: ii is momentum [1:nptk], j is energy [1:Nband]
       ! phonon3: jj is momentum [ii:nptk], k is energy [startbranch1:Nband]
       ! phonon4: ss is momentum [jj:nptk], l is energy [startbranch2:Nband]
       ! total combination: Nbands*nptk*Nbands*nptk*Nbands/2
       ! j = n + FLOOR((m+1-n)*u)  ! equation to generate random number j in [n,m] from a uniform distribution u in [0,1]
 
       call random_seed()
       call random_number(rand_matrix)
       iter=0
       do while (iter<=num_sample_process_4ph)
          matrix_iter = (iter-1)*5
          ii=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+1)) ! phonon2 momentum [1:nptk]
          j=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+2)) ! phonon2 branch [1:Nband]
          jj=1+FLOOR(INT(nptk+1-1)*rand_matrix(matrix_iter+3)) ! phonon3 momentum [ii:nptk]
          startbranch1=1
          if(jj==ii) startbranch1=j
          k=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+4)) ! phonon3 branch [startbranch1:Nband]
 
          qprime=IJK(:,ii)
          realqprime=matmul(rlattvec,qprime/dble(ngrid))
          qdprime=IJK(:,jj)
          realqdprime=matmul(rlattvec,qdprime/dble(ngrid))
          qtprime=modulo(q-qprime-qdprime,ngrid)
          realqtprime=matmul(rlattvec,qtprime/dble(ngrid))
          ss=Index_N(qtprime(1),qtprime(2),qtprime(3)) ! phonon4 momentum [jj:nptk]
          startbranch2=1
          if(ss==jj) startbranch2=k
          l=1+FLOOR(INT(Nbands+1-1)*rand_matrix(matrix_iter+5)) ! phonon4 branch [startbranch2:Nband]
          ! ----------- end sampling method add -----------
 
          if (jj>=ii .AND. k>=startbranch1 .AND. ss>=jj .AND. l>=startbranch2) then ! make sure that the phonon3 is in [ii:nptk], phonon3 branch is in [startbranch:Nband]
 
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
                   end if
                end if
             end if
          end if ! check momentum and energy conservation
          iter = iter + 1
       end do ! iter
       WP4_minusminus=WP4_minusminus/nptk/nptk
     end if
     
     ! converting to THz
     Gamma_minusminus=Gamma_minusminus*3.37617087*1.d9/nptk/nptk
 
    ! ----------- sampling method add -----------
    Gamma_minusminus = Gamma_minusminus/num_sample_process_4ph*total_process  
    WP4_minusminus = WP4_minusminus*total_process/num_sample_process_4ph
    ! ----------- end sampling method add -----------
 
 
 
   end subroutine
 
   ! Wrapper around 4ph RTA subroutines with 3ph subroutines that splits the work among processors.
   subroutine RTA_driver_4ph(energy,velocity,eigenvect,Nlist,List,IJK,&
        Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,&
        rate_scatt_4ph,rate_scatt_plusplus,rate_scatt_plusminus,rate_scatt_minusminus,&
        WP4_plusplus,WP4_plusminus,WP4_minusminus)
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
      integer(kind=4) :: i
      integer(kind=4) :: ll
      integer(kind=4) :: mm
      ! data for output 
      real(kind=8),intent(out) :: rate_scatt_4ph(Nbands,Nlist)
      real(kind=8),intent(out) :: rate_scatt_plusplus(Nbands,Nlist),rate_scatt_plusminus(Nbands,Nlist),rate_scatt_minusminus(Nbands,Nlist)
      real(kind=8),intent(out) :: WP4_plusplus(Nbands,Nlist),WP4_plusminus(Nbands,Nlist),WP4_minusminus(Nbands,Nlist)
      !  end data for rate_scatt_4ph
   
      ! reduced variables
      real(kind=8) :: rate_scatt_plusplus_reduce(Nbands,Nlist),rate_scatt_plusminus_reduce(Nbands,Nlist),rate_scatt_minusminus_reduce(Nbands,Nlist)
      real(kind=8) :: WP4_plusplus_reduce(Nbands*Nlist),WP4_plusminus_reduce(Nbands*Nlist),WP4_minusminus_reduce(Nbands*Nlist)
      ! end reduced variables

         
      !   temporary variables
      real(kind=8) :: Gamma_plusplus,Gamma_plusminus,Gamma_minusminus
      real(kind=8) :: WP4_temp
      !   end temporary variables

      integer(kind=4) :: chunkstart, chunkend, chunkid ! for parallel computing

      ! initialize the output arrays
      rate_scatt_4ph=0.d00
      rate_scatt_plusplus=0.d00
      rate_scatt_plusminus=0.d00
      rate_scatt_minusminus=0.d00

      WP4_plusplus=0.d00
      WP4_plusminus=0.d00
      WP4_minusminus=0.d00
      !   end initialize the output arrays

      ! initialize the reduced output arrays
      rate_scatt_plusplus_reduce=0.d00
      rate_scatt_plusminus_reduce=0.d00
      rate_scatt_minusminus_reduce=0.d00

      WP4_plusplus_reduce=0.d00
      WP4_plusminus_reduce=0.d00
      WP4_minusminus_reduce=0.d00
      !   end initialize the reduced output arrays

      ! initialize for the temporary variables
      Gamma_plusplus=0.d0
      Gamma_plusminus=0.d0
      Gamma_minusminus=0.d0
      WP4_temp=0.d0
      !   end initialize for the temporary variables

      do chunkid=myid+1,numchunk,numprocs
         chunkstart = (chunkid-1)*chunksize+1
         chunkend = chunkid*chunksize
         if (chunkend.gt.Nbands*Nlist) chunkend = Nbands*Nlist
         !$OMP PARALLEL DO default(shared) schedule(dynamic,1) private(mm,i,ll,Wp4_temp,Gamma_plusplus) &
         !$OMP & private(Gamma_plusminus,Gamma_minusminus)
         do mm=chunkstart,chunkend
            i=modulo(mm-1,Nbands)+1
            ll=int((mm-1)/Nbands)+1
            if (energy(List(ll),i).le.omega_max) then
               if (num_sample_process_4ph==-1) then ! do not sample
                  call RTA_plusplus(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                        Gamma_plusplus,WP4_temp)
                     WP4_plusplus_reduce(mm)=WP4_temp
                     rate_scatt_plusplus_reduce(i,ll)=Gamma_plusplus
                  call RTA_plusminus(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                        Gamma_plusminus,WP4_temp)
                     WP4_plusminus_reduce(mm)=WP4_temp
                     rate_scatt_plusminus_reduce(i,ll)=Gamma_plusminus
                  call RTA_minusminus(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                        Gamma_minusminus,WP4_temp)
                     WP4_minusminus_reduce(mm)=WP4_temp
                     rate_scatt_minusminus_reduce(i,ll)=Gamma_minusminus
               else ! do sampling method
                  call RTA_plusplus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                        Gamma_plusplus,WP4_temp)
                     WP4_plusplus_reduce(mm)=WP4_temp
                  rate_scatt_plusplus_reduce(i,ll)=Gamma_plusplus
                  call RTA_plusminus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                        Gamma_plusminus,WP4_temp)
                     WP4_plusminus_reduce(mm)=WP4_temp
                  rate_scatt_plusminus_reduce(i,ll)=Gamma_plusminus
                  call RTA_minusminus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                        Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                        Gamma_minusminus,WP4_temp)
                     WP4_minusminus_reduce(mm)=WP4_temp
                  rate_scatt_minusminus_reduce(i,ll)=Gamma_minusminus
               endif
      
            endif
         end do
         !$OMP END PARALLEL DO
      end do
   
      call MPI_ALLREDUCE(rate_scatt_plusplus_reduce,rate_scatt_plusplus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
      call MPI_ALLREDUCE(rate_scatt_plusminus_reduce,rate_scatt_plusminus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
      call MPI_ALLREDUCE(rate_scatt_minusminus_reduce,rate_scatt_minusminus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

      call MPI_ALLREDUCE(WP4_plusplus_reduce,WP4_plusplus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
      call MPI_ALLREDUCE(WP4_plusminus_reduce,WP4_plusminus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)
      call MPI_ALLREDUCE(WP4_minusminus_reduce,WP4_minusminus,Nbands*Nlist,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mm)

      rate_scatt_4ph=rate_scatt_plusplus+rate_scatt_plusminus+rate_scatt_minusminus
   end subroutine RTA_driver_4ph



   subroutine RTA_driver_4ph_GPU_using_Ind(energy,velocity,eigenvect,Nlist,List,IJK,&
      Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l, &
      N_plusplus,N_plusminus,N_minusminus,&
      Indof2ndPhonon_plusplus,Indof3rdPhonon_plusplus,Indof4thPhonon_plusplus,&
      Indof2ndPhonon_plusminus,Indof3rdPhonon_plusminus,Indof4thPhonon_plusminus,&
      Indof2ndPhonon_minusminus,Indof3rdPhonon_minusminus,Indof4thPhonon_minusminus,&
      rate_scatt_4ph,&
      rate_scatt_plusplus,rate_scatt_plusminus,rate_scatt_minusminus, &
      WP4_plusplus,WP4_plusminus,WP4_minusminus)
      implicit none
      include "mpif.h"

      ! Input variables 
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

      integer(kind=8),intent(in) :: N_plusplus(Nlist*Nbands)
      integer(kind=8),intent(in) :: N_plusminus(Nlist*Nbands)
      integer(kind=8),intent(in) :: N_minusminus(Nlist*Nbands)

      integer(kind=8),intent(in) :: Indof2ndPhonon_plusplus(:),Indof3rdPhonon_plusplus(:),Indof4thPhonon_plusplus(:) 
      integer(kind=8),intent(in) :: Indof2ndPhonon_plusminus(:),Indof3rdPhonon_plusminus(:),Indof4thPhonon_plusminus(:)
      integer(kind=8),intent(in) :: Indof2ndPhonon_minusminus(:),Indof3rdPhonon_minusminus(:),Indof4thPhonon_minusminus(:)

      ! Output variables
      real(kind=8),intent(out) :: rate_scatt_4ph(Nbands,Nlist)
      real(kind=8),intent(out) :: rate_scatt_plusplus(Nbands,Nlist)
      real(kind=8),intent(out) :: rate_scatt_plusminus(Nbands,Nlist)
      real(kind=8),intent(out) :: rate_scatt_minusminus(Nbands,Nlist)
      real(kind=8),intent(out) :: WP4_plusplus(Nbands,Nlist)
      real(kind=8),intent(out) :: WP4_plusminus(Nbands,Nlist)
      real(kind=8),intent(out) :: WP4_minusminus(Nbands,Nlist)

      if (num_sample_process_4ph .eq. -1) then  ! do not use sampling
         if (myid .eq. 0) then
            call RTA_plusplus_driver_GPU_using_Ind(energy,velocity,eigenvect,Nlist,List, &
               Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK, &
                  N_plusplus, &
                  Indof2ndPhonon_plusplus,Indof3rdPhonon_plusplus,Indof4thPhonon_plusplus,&
               rate_scatt_plusplus,WP4_plusplus)
            call RTA_plusminus_driver_GPU_using_Ind(energy,velocity,eigenvect,Nlist,List, &
               Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK, &
               N_plusminus, &
               Indof2ndPhonon_plusminus,Indof3rdPhonon_plusminus,Indof4thPhonon_plusminus, &
               rate_scatt_plusminus,WP4_plusminus)
            call RTA_minusminus_driver_GPU_using_Ind(energy,velocity,eigenvect,Nlist,List, &
               Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK, &
               N_minusminus, &
               Indof2ndPhonon_minusminus,Indof3rdPhonon_minusminus,Indof4thPhonon_minusminus, &
               rate_scatt_minusminus,WP4_minusminus)
         end if

      else ! do sampling method
         print *, "Error: Sampling method is not implemented for GPU version."
      endif

      ! Combine the contributions from the three processes into the total 4ph rate.
      rate_scatt_4ph = rate_scatt_plusplus + rate_scatt_plusminus + rate_scatt_minusminus

   end subroutine RTA_driver_4ph_GPU_using_Ind

 
 end module processes
