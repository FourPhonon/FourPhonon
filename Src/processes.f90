! Compute the number of allowed three-phonon / four-phonon processes, their
! scattering amplitudes and their phase-space volume.

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
     integer(kind=4),intent(out), allocatable :: Indof2ndPhonon_plus(:),Indof3rdPhonon_plus(:)
     real(kind=8),intent(out), allocatable :: Gamma_plus(:),Gamma_minus(:)
     integer(kind=4),intent(out), allocatable :: Indof2ndPhonon_minus(:),Indof3rdPhonon_minus(:)
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
 
     real(kind=8) :: rate_scatt_plus_reduce(Nbands,Nlist),rate_scatt_minus_reduce(Nbands,Nlist)
     real(kind=8),allocatable :: Gamma0(:)
 
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
 
     allocate(Indof2ndPhonon_plus(Ntotal_plus))
     allocate(Indof3rdPhonon_plus(Ntotal_plus))
     allocate(Indof2ndPhonon_minus(Ntotal_minus))
     allocate(Indof3rdPhonon_minus(Ntotal_minus))
     allocate(Gamma_plus(Ntotal_plus))
     allocate(Gamma_minus(Ntotal_minus))
 
     WP3_plus=0.d0
     WP3_minus=0.d0
     rate_scatt_plus=0.d0
     rate_scatt_minus=0.d0
     Gamma_plus=0.d0
     Gamma_minus=0.d0
     Indof2ndPhonon_plus=0
     Indof3rdPhonon_plus=0
     Indof2ndPhonon_minus=0
     Indof3rdPhonon_minus=0
 
     !$OMP PARALLEL DO default(shared) private(mm,i,ll,Gamma0,Indof2ndPhonon,Indof3rdPhonon) 
     do mm=1,Nbands*NList
        i=modulo(mm-1,Nbands)+1
        ll=int((mm-1)/Nbands)+1
        if(N_plus(mm).ne.0) then
          call Ind_plus(mm,N_plus(mm),energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
               Indof2ndPhonon(1:N_plus(mm)),Indof3rdPhonon(1:N_plus(mm)),&
               Gamma0(1:N_plus(mm)),WP3_plus(i,ll))
          Indof2ndPhonon_plus((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Indof2ndPhonon(1:N_plus(mm))
          Indof3rdPhonon_plus((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Indof3rdPhonon(1:N_plus(mm))
          Gamma_plus((Naccum_plus(mm)+1):(Naccum_plus(mm)+N_plus(mm)))=&
               Gamma0(1:N_plus(mm))
          rate_scatt_plus(i,ll)=sum(Gamma0(1:N_plus(mm)))
       end if
       if(N_minus(mm).ne.0) then
          call Ind_minus(mm,N_minus(mm),energy,velocity,eigenvect,Nlist,List,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
               Indof2ndPhonon(1:N_minus(mm)),Indof3rdPhonon(1:N_minus(mm)),&
               Gamma0(1:N_minus(mm)),WP3_minus(i,ll))
          Indof2ndPhonon_minus((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Indof2ndPhonon(1:N_minus(mm))
          Indof3rdPhonon_minus((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Indof3rdPhonon(1:N_minus(mm))
          Gamma_minus((Naccum_minus(mm)+1):(Naccum_minus(mm)+N_minus(mm)))=&
               Gamma0(1:N_minus(mm))
          rate_scatt_minus(i,ll)=sum(Gamma0(1:N_minus(mm)))*5.D-1
       end if
     end do
     !$OMP END PARALLEL DO
 
     deallocate(Gamma0)
     deallocate(Indof3rdPhonon)
     deallocate(Indof2ndPhonon)
 
     rate_scatt=rate_scatt_plus+rate_scatt_minus
 
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
     real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3)
 
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
     real(kind=8) :: shortest_q(3),shortest_qprime(3),shortest_qdprime(3)
 
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
 
     Pspace_plus_total=0.d0
     Pspace_minus_total=0.d0
     N_plus=0
     N_minus=0
     
     !$OMP PARALLEL DO default(shared) private(mm,i,ll)
     do mm=1,Nbands*Nlist
        i=modulo(mm-1,Nbands)+1
        ll=int((mm-1)/Nbands)+1
        if (energy(List(int((mm-1)/Nbands)+1),modulo(mm-1,Nbands)+1).le.omega_max) then
          if (num_sample_process_3ph_phase_space==-1) then ! do not sample
             call NP_plus(mm,energy,velocity,Nlist,List,IJK,&
                N_plus(mm),Pspace_plus_total(i,ll))
             call NP_minus(mm,energy,velocity,Nlist,List,IJK,&
                N_minus(mm),Pspace_minus_total(i,ll))
          else ! do sampling method
             call NP_plus_sample(mm,energy,velocity,Nlist,List,IJK,&
                N_plus(mm),Pspace_plus_total(i,ll))
             call NP_minus_sample(mm,energy,velocity,Nlist,List,IJK,&
                N_minus(mm),Pspace_minus_total(i,ll))
          endif
        endif
     end do
     !$OMP END PARALLEL DO
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
 
   ! RTA-only version of Ind_plus using sampling method.
   subroutine RTA_plus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
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
 
    ! ----------- sampling method add -----------
     integer(kind=4) :: nn, rand_num, iter
     real :: rand_matrix(num_sample_process_3ph),rand
    ! ----------- end sampling method add -----------
 
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
          end if ! end energy conservation
          !--------END absorption process-------------!
       end do ! iter
        WP3_plus=WP3_plus/nptk
     end if
    ! ----------- sampling method add -----------
     Gamma_plus = Gamma_plus*Nbands*nptk*Nbands/num_sample_process_3ph
     Gamma_plus_N = Gamma_plus_N*Nbands*nptk*Nbands/num_sample_process_3ph
     Gamma_plus_U = Gamma_plus_U*Nbands*nptk*Nbands/num_sample_process_3ph
     WP3_plus=WP3_plus*Nbands*nptk*Nbands/num_sample_process_3ph
    ! ----------- end sampling method add -----------
 
     Gamma_plus=Gamma_plus*5.60626442*1.d8/nptk ! THz
     Gamma_plus_N=Gamma_plus_N*5.60626442*1.d8/nptk ! THz
     Gamma_plus_U=Gamma_plus_U*5.60626442*1.d8/nptk ! THz
   end subroutine RTA_plus_sample
 
   ! RTA-only version of Ind_minus using sampling method.
   subroutine RTA_minus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
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
 
    ! ----------- sampling method add -----------
     integer(kind=4) :: nn, rand_num, iter
     real :: rand_matrix(num_sample_process_3ph),rand
    ! ----------- end sampling method add -----------
 
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
          end if ! end energy conservation
          !--------END emission process-------------
       end do ! iter
        WP3_minus=WP3_minus*5.d-1/nptk
     end if
 
    ! ----------- sampling method add -----------
     Gamma_minus = Gamma_minus*Nbands*nptk*Nbands/num_sample_process_3ph
     Gamma_minus_N = Gamma_minus_N*Nbands*nptk*Nbands/num_sample_process_3ph
     Gamma_minus_U = Gamma_minus_U*Nbands*nptk*Nbands/num_sample_process_3ph
     WP3_minus = WP3_minus*Nbands*nptk*Nbands/num_sample_process_3ph
    ! ----------- end sampling method add -----------
 
     Gamma_minus=Gamma_minus*5.60626442*1.d8/nptk
     Gamma_minus_N=Gamma_minus_N*5.60626442*1.d8/nptk
     Gamma_minus_U=Gamma_minus_U*5.60626442*1.d8/nptk
 
   end subroutine RTA_minus_sample
 
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
    !  real(kind=8) :: rate_scatt_plus_reduce(Nbands,Nlist),rate_scatt_minus_reduce(Nbands,Nlist)
    !  real(kind=8) :: rate_scatt_plus_reduce_N(Nbands,Nlist),rate_scatt_minus_reduce_N(Nbands,Nlist)
    !  real(kind=8) :: rate_scatt_plus_reduce_U(Nbands,Nlist),rate_scatt_minus_reduce_U(Nbands,Nlist) 
    !  real(kind=8) :: WP3_plus_reduce(Nbands*Nlist)
    !  real(kind=8) :: WP3_minus_reduce(Nbands*Nlist)
 
 
     rate_scatt=0.d00
     rate_scatt_plus=0.d00
     rate_scatt_minus=0.d00
     rate_scatt_plus_N=0.d00
     rate_scatt_plus_U=0.d00
     rate_scatt_minus_N=0.d00
     rate_scatt_minus_U=0.d00
 
     WP3_plus=0.d00
     WP3_minus=0.d00
 
     !$OMP PARALLEL DO default(shared) private(mm,i,ll,Gamma_plus,Gamma_plus_N,Gamma_plus_U,Gamma_minus,Gamma_minus_N,Gamma_minus_U)
     do mm=1,Nbands*NList
        i=modulo(mm-1,Nbands)+1
        ll=int((mm-1)/Nbands)+1
        if (energy(List(ll),i).le.omega_max) then
          if (num_sample_process_3ph==-1) then ! do not sample
             call RTA_plus(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                   Gamma_plus,WP3_plus(i,ll),&
                   Gamma_plus_N,Gamma_plus_U)
             rate_scatt_plus(i,ll)=Gamma_plus
             rate_scatt_plus_N(i,ll)=Gamma_plus_N
             rate_scatt_plus_U(i,ll)=Gamma_plus_U
             call RTA_minus(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                   Gamma_minus,WP3_minus(i,ll),Gamma_minus_N,Gamma_minus_U)
             rate_scatt_minus(i,ll)=Gamma_minus*5.D-1
             rate_scatt_minus_N(i,ll)=Gamma_minus_N*5.D-1
             rate_scatt_minus_U(i,ll)=Gamma_minus_U*5.D-1
          else ! do sampling method
             call RTA_plus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                   Gamma_plus,WP3_plus(i,ll),&
                   Gamma_plus_N,Gamma_plus_U)
             rate_scatt_plus(i,ll)=Gamma_plus
             rate_scatt_plus_N(i,ll)=Gamma_plus_N
             rate_scatt_plus_U(i,ll)=Gamma_plus_U
             call RTA_minus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,IJK,&
                   Gamma_minus,WP3_minus(i,ll),Gamma_minus_N,Gamma_minus_U)
             rate_scatt_minus(i,ll)=Gamma_minus*5.D-1
             rate_scatt_minus_N(i,ll)=Gamma_minus_N*5.D-1
             rate_scatt_minus_U(i,ll)=Gamma_minus_U*5.D-1
 
          endif
 
        endif
     end do
     !$OMP END PARALLEL DO
 
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
 
     Pspace_plusplus_total=0.d0
     Pspace_plusminus_total=0.d0
     Pspace_minusminus_total=0.d0
     N_plusplus=0
     N_plusminus=0
     N_minusminus=0
 
     !$OMP PARALLEL DO default(shared) private(mm,i,ll)
     do mm=1,Nbands*NList
       i=modulo(mm-1,Nbands)+1
       ll=int((mm-1)/Nbands)+1
       if (energy(List(int((mm-1)/Nbands)+1),modulo(mm-1,Nbands)+1).le.omega_max) then
          if (num_sample_process_4ph_phase_space==-1) then ! do not sample
             call NP_plusplus(mm,energy,velocity,Nlist,List,IJK,N_plusplus(mm),Pspace_plusplus_total(i,ll))
             call NP_plusminus(mm,energy,velocity,Nlist,List,IJK,N_plusminus(mm),Pspace_plusminus_total(i,ll))
             call NP_minusminus(mm,energy,velocity,Nlist,List,IJK,N_minusminus(mm),Pspace_minusminus_total(i,ll))
          else ! do sampling method
             call NP_plusplus_sample(mm,energy,velocity,Nlist,List,IJK,N_plusplus(mm),Pspace_plusplus_total(i,ll))
             call NP_plusminus_sample(mm,energy,velocity,Nlist,List,IJK,N_plusminus(mm),Pspace_plusminus_total(i,ll))
             call NP_minusminus_sample(mm,energy,velocity,Nlist,List,IJK,N_minusminus(mm),Pspace_minusminus_total(i,ll))
          endif
       endif
     end do
     !$OMP END PARALLEL DO
   end subroutine NP_driver_4ph
 
 
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
     integer(kind=8),intent(out), allocatable :: Indof2ndPhonon_plusplus(:),Indof3rdPhonon_plusplus(:),Indof4thPhonon_plusplus(:) 
     integer(kind=8),intent(out), allocatable :: Indof2ndPhonon_plusminus(:),Indof3rdPhonon_plusminus(:),Indof4thPhonon_plusminus(:)
     integer(kind=8),intent(out), allocatable :: Indof2ndPhonon_minusminus(:),Indof3rdPhonon_minusminus(:),Indof4thPhonon_minusminus(:)
     real(kind=8),intent(out), allocatable :: Gamma_plusplus(:),Gamma_plusminus(:),Gamma_minusminus(:)
     real(kind=8),intent(out) :: rate_scatt_4ph(Nbands,Nlist)
     real(kind=8),intent(out) :: rate_scatt_plusplus(Nbands,Nlist),rate_scatt_plusminus(Nbands,Nlist),rate_scatt_minusminus(Nbands,Nlist)
     real(kind=8),intent(out) :: WP4_plusplus(Nbands,Nlist),WP4_plusminus(Nbands,Nlist),WP4_minusminus(Nbands,Nlist)
 
     integer(kind=4) :: i,nthread
     integer(kind=4) :: ll
     integer(kind=4) :: mm
     integer(kind=4) :: maxsize
     integer(kind=8) :: Ntotal_plusplus
     integer(kind=8) :: Ntotal_plusminus
     integer(kind=8) :: Ntotal_minusminus
     integer(kind=8) :: Naccum_plusplus(Nbands*Nlist)
     integer(kind=8) :: Naccum_plusminus(Nbands*Nlist)
     integer(kind=8) :: Naccum_minusminus(Nbands*Nlist)
     integer(kind=8),allocatable :: Indof2ndPhonon(:)
     integer(kind=8),allocatable :: Indof3rdPhonon(:)
     integer(kind=8),allocatable :: Indof4thPhonon(:)
     real(kind=8),allocatable :: Gamma0(:)
 
 
     maxsize=max(maxval(N_plusplus),maxval(N_plusminus),maxval(N_minusminus))
     allocate(Indof2ndPhonon(maxsize))
     allocate(Indof3rdPhonon(maxsize))
     allocate(Indof4thPhonon(maxsize))
     allocate(Gamma0(maxsize))
 
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
 
     allocate(Indof2ndPhonon_plusplus(Ntotal_plusplus))
     allocate(Indof3rdPhonon_plusplus(Ntotal_plusplus))
     allocate(Indof4thPhonon_plusplus(Ntotal_plusplus))
     allocate(Indof2ndPhonon_plusminus(Ntotal_plusminus))
     allocate(Indof3rdPhonon_plusminus(Ntotal_plusminus))
     allocate(Indof4thPhonon_plusminus(Ntotal_plusminus))
     allocate(Indof2ndPhonon_minusminus(Ntotal_minusminus))
     allocate(Indof3rdPhonon_minusminus(Ntotal_minusminus))
     allocate(Indof4thPhonon_minusminus(Ntotal_minusminus))
     allocate(Gamma_plusplus(Ntotal_plusplus))
     allocate(Gamma_plusminus(Ntotal_plusminus))
     allocate(Gamma_minusminus(Ntotal_minusminus))
 
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
 
     
     !$OMP PARALLEL DO default(shared) private(mm,i,ll,Gamma0,Indof2ndPhonon,Indof3rdPhonon,Indof4thPhonon)      
     do mm=1,Nbands*NList ! This loop is for pure Openmp parallel
       !  if ( mm == Nbands*NList/4 ) then
       !     print *, "We are at a quater of all irreducible modes."
       !  end if
        i=modulo(mm-1,Nbands)+1
        ll=int((mm-1)/Nbands)+1
        if(N_plusplus(mm).ne.0) then
           call Ind_plusplus(mm,N_plusplus(mm),energy,velocity,eigenvect,Nlist,List,&
                Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                Indof2ndPhonon(1:N_plusplus(mm)),Indof3rdPhonon(1:N_plusplus(mm)),Indof4thPhonon(1:N_plusplus(mm)),&
                Gamma0(1:N_plusplus(mm)),WP4_plusplus(i,ll))
             Indof2ndPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                Indof2ndPhonon(1:N_plusplus(mm))
             Indof3rdPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                Indof3rdPhonon(1:N_plusplus(mm))
             Indof4thPhonon_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                Indof4thPhonon(1:N_plusplus(mm))
             Gamma_plusplus((Naccum_plusplus(mm)+1):(Naccum_plusplus(mm)+N_plusplus(mm)))=&
                Gamma0(1:N_plusplus(mm))
             rate_scatt_plusplus(i,ll)=sum(Gamma0(1:N_plusplus(mm)))
        end if
        if(N_plusminus(mm).ne.0) then
           call Ind_plusminus(mm,N_plusminus(mm),energy,velocity,eigenvect,Nlist,List,&
                Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                Indof2ndPhonon(1:N_plusminus(mm)),Indof3rdPhonon(1:N_plusminus(mm)),Indof4thPhonon(1:N_plusminus(mm)),&
                Gamma0(1:N_plusminus(mm)),WP4_plusminus(i,ll))
             Indof2ndPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                Indof2ndPhonon(1:N_plusminus(mm))
             Indof3rdPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                Indof3rdPhonon(1:N_plusminus(mm))
             Indof4thPhonon_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                Indof4thPhonon(1:N_plusminus(mm))
             Gamma_plusminus((Naccum_plusminus(mm)+1):(Naccum_plusminus(mm)+N_plusminus(mm)))=&
                Gamma0(1:N_plusminus(mm))
             rate_scatt_plusminus(i,ll)=sum(Gamma0(1:N_plusminus(mm)))
        end if
        if(N_minusminus(mm).ne.0) then
           call Ind_minusminus(mm,N_minusminus(mm),energy,velocity,eigenvect,Nlist,List,&
                Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                Indof2ndPhonon(1:N_minusminus(mm)),Indof3rdPhonon(1:N_minusminus(mm)),Indof4thPhonon(1:N_minusminus(mm)),&
                Gamma0(1:N_minusminus(mm)),WP4_minusminus(i,ll))
             Indof2ndPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                Indof2ndPhonon(1:N_minusminus(mm))
             Indof3rdPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                Indof3rdPhonon(1:N_minusminus(mm))
             Indof4thPhonon_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                Indof4thPhonon(1:N_minusminus(mm))
             Gamma_minusminus((Naccum_minusminus(mm)+1):(Naccum_minusminus(mm)+N_minusminus(mm)))=&
                Gamma0(1:N_minusminus(mm))
             rate_scatt_minusminus(i,ll)=sum(Gamma0(1:N_minusminus(mm)))
        end if
     end do
     !$OMP END PARALLEL DO
 
     deallocate(Gamma0)
     deallocate(Indof4thPhonon)
     deallocate(Indof3rdPhonon)
     deallocate(Indof2ndPhonon)
 
     rate_scatt_4ph=rate_scatt_plusplus+rate_scatt_plusminus+rate_scatt_minusminus
 
   end subroutine Ind_driver_4ph
 
 
 
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
                   ! end if
                end if ! greater equal than 1
             ! end do
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
                   ! end if
                end if ! greater equal than jj
             ! end do
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
                   ! end if
                end if ! greater than jj
             ! end do
          end do
       end do
       WP4_minusminus=WP4_minusminus/nptk/nptk
     end if
     
     ! converting to THz
     Gamma_minusminus=Gamma_minusminus*3.37617087*1.d9/nptk/nptk
     Gamma_minusminus_N=Gamma_minusminus_N*3.37617087*1.d9/nptk/nptk
     Gamma_minusminus_U=Gamma_minusminus_U*3.37617087*1.d9/nptk/nptk
 
   end subroutine
   ! RTA-only version of 4ph ++ process using sampling method, Ind_plusplus; Avoid double counting
   subroutine RTA_plusplus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
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
 
    ! ----------- sampling method add -----------
    integer(kind=8) :: rand_num, iter, total_process, matrix_iter, temp_int
    real :: rand_matrix(num_sample_process_4ph*5)
    total_process = Nbands*nptk*Nbands*nptk*Nbands
    ! ----------- end sampling method add -----------
 
 
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
          end if ! check momentum and energy conservation
          iter = iter + 1
       end do ! iter                        
       WP4_plusplus=WP4_plusplus/nptk/nptk
     end if
     ! converting to THz
     Gamma_plusplus=Gamma_plusplus*3.37617087*1.d9/nptk/nptk
     Gamma_plusplus_N=Gamma_plusplus_N*3.37617087*1.d9/nptk/nptk
     Gamma_plusplus_U=Gamma_plusplus_U*3.37617087*1.d9/nptk/nptk
 
    ! ----------- sampling method add -----------
    Gamma_plusplus = Gamma_plusplus/num_sample_process_4ph*total_process  
    Gamma_plusplus_N = Gamma_plusplus_N*total_process/num_sample_process_4ph
    Gamma_plusplus_U = Gamma_plusplus_U*total_process/num_sample_process_4ph
    WP4_plusplus = WP4_plusplus*total_process/num_sample_process_4ph
    ! ----------- end sampling method add -----------
 
 
   end subroutine
 
   ! RTA-only version of 4ph +- process using sampling method, Ind_plusminus
   subroutine RTA_plusminus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
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
 
    ! ----------- sampling method add -----------
    integer(kind=8) :: rand_num, iter, total_process, matrix_iter, temp_int
    real :: rand_matrix(num_sample_process_4ph*5)
    total_process = Nbands*nptk*Nbands*nptk*Nbands
    ! ----------- end sampling method add -----------
 
 
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
          end if ! check momentum and energy conservation
          iter = iter + 1
       end do ! iter
       WP4_plusminus=WP4_plusminus/nptk/nptk
     end if
     
     ! converting to THz
     Gamma_plusminus=Gamma_plusminus*3.37617087*1.d9/nptk/nptk
     Gamma_plusminus_N=Gamma_plusminus_N*3.37617087*1.d9/nptk/nptk
     Gamma_plusminus_U=Gamma_plusminus_U*3.37617087*1.d9/nptk/nptk
 
    ! ----------- sampling method add -----------
    Gamma_plusminus = Gamma_plusminus/num_sample_process_4ph*total_process    
    Gamma_plusminus_N = Gamma_plusminus_N*total_process/num_sample_process_4ph
    Gamma_plusminus_U = Gamma_plusminus_U*total_process/num_sample_process_4ph
    WP4_plusminus = WP4_plusminus*total_process/num_sample_process_4ph
    ! ----------- end sampling method add -----------
 
 
   end subroutine
 
   ! RTA-only version of 4ph -- process using sampling method, Ind_minusminus
   subroutine RTA_minusminus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
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
          end if ! check momentum and energy conservation
          iter = iter + 1
       end do ! iter
       WP4_minusminus=WP4_minusminus/nptk/nptk
     end if
     
     ! converting to THz
     Gamma_minusminus=Gamma_minusminus*3.37617087*1.d9/nptk/nptk
     Gamma_minusminus_N=Gamma_minusminus_N*3.37617087*1.d9/nptk/nptk
     Gamma_minusminus_U=Gamma_minusminus_U*3.37617087*1.d9/nptk/nptk
 
    ! ----------- sampling method add -----------
    Gamma_minusminus = Gamma_minusminus/num_sample_process_4ph*total_process  
    Gamma_minusminus_N = Gamma_minusminus_N*total_process/num_sample_process_4ph
    Gamma_minusminus_U = Gamma_minusminus_U*total_process/num_sample_process_4ph
    WP4_minusminus = WP4_minusminus*total_process/num_sample_process_4ph
    ! ----------- end sampling method add -----------
 
 
 
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
    !  real(kind=8) :: rate_scatt_plusplus_reduce(Nbands,Nlist),rate_scatt_plusminus_reduce(Nbands,Nlist),rate_scatt_minusminus_reduce(Nbands,Nlist)
    !  real(kind=8) :: rate_scatt_plusplus_reduce_N(Nbands,Nlist),rate_scatt_plusminus_reduce_N(Nbands,Nlist),rate_scatt_minusminus_reduce_N(Nbands,Nlist)
    !  real(kind=8) :: rate_scatt_plusplus_reduce_U(Nbands,Nlist),rate_scatt_plusminus_reduce_U(Nbands,Nlist),rate_scatt_minusminus_reduce_U(Nbands,Nlist) 
    !  real(kind=8) :: WP4_plusplus_reduce(Nbands*Nlist),WP4_plusminus_reduce(Nbands*Nlist),WP4_minusminus_reduce(Nbands*Nlist)
 
     rate_scatt_4ph=0.d00
     rate_scatt_plusplus=0.d00
     rate_scatt_plusplus_N=0.d00
     rate_scatt_plusplus_U=0.d00
 
     rate_scatt_plusminus=0.d00
     rate_scatt_plusminus_N=0.d00
     rate_scatt_plusminus_U=0.d00
 
     rate_scatt_minusminus=0.d00
     rate_scatt_minusminus_N=0.d00
     rate_scatt_minusminus_U=0.d00
 
     WP4_plusplus=0.d00
     WP4_plusminus=0.d00
     WP4_minusminus=0.d00
 
     !$OMP PARALLEL DO default(shared) private(mm,i,ll,Gamma_plusplus,Gamma_plusplus_N,Gamma_plusplus_U) &
     !$OMP & private(Gamma_plusminus,Gamma_plusminus_N,Gamma_plusminus_U,Gamma_minusminus,Gamma_minusminus_N,Gamma_minusminus_U)
     do mm=1,Nbands*NList ! This loop is for pure Openmp parallel
        i=modulo(mm-1,Nbands)+1
        ll=int((mm-1)/Nbands)+1
        if (energy(List(ll),i).le.omega_max) then
          if (num_sample_process_4ph==-1) then ! do not sample
             call RTA_plusplus(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                   Gamma_plusplus,WP4_plusplus(i,ll),Gamma_plusplus_N,Gamma_plusplus_U)
             rate_scatt_plusplus(i,ll)=Gamma_plusplus
             rate_scatt_plusplus_N(i,ll)=Gamma_plusplus_N
             rate_scatt_plusplus_U(i,ll)=Gamma_plusplus_U
             call RTA_plusminus(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                   Gamma_plusminus,WP4_plusminus(i,ll),Gamma_plusminus_N,Gamma_plusminus_U)
             rate_scatt_plusminus(i,ll)=Gamma_plusminus
             rate_scatt_plusminus_N(i,ll)=Gamma_plusminus_N
             rate_scatt_plusminus_U(i,ll)=Gamma_plusminus_U
             call RTA_minusminus(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                   Gamma_minusminus,WP4_minusminus(i,ll),Gamma_minusminus_N,Gamma_minusminus_U)
             rate_scatt_minusminus(i,ll)=Gamma_minusminus
             rate_scatt_minusminus_N(i,ll)=Gamma_minusminus_N
             rate_scatt_minusminus_U(i,ll)=Gamma_minusminus_U
          else ! do sampling method
             call RTA_plusplus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                   Gamma_plusplus,WP4_plusplus(i,ll),Gamma_plusplus_N,Gamma_plusplus_U)
             rate_scatt_plusplus(i,ll)=Gamma_plusplus
             rate_scatt_plusplus_N(i,ll)=Gamma_plusplus_N
             rate_scatt_plusplus_U(i,ll)=Gamma_plusplus_U
             call RTA_plusminus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                   Gamma_plusminus,WP4_plusminus(i,ll),Gamma_plusminus_N,Gamma_plusminus_U)
             rate_scatt_plusminus(i,ll)=Gamma_plusminus
             rate_scatt_plusminus_N(i,ll)=Gamma_plusminus_N
             rate_scatt_plusminus_U(i,ll)=Gamma_plusminus_U
             call RTA_minusminus_sample(mm,energy,velocity,eigenvect,Nlist,List,&
                   Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l,IJK,&
                   Gamma_minusminus,WP4_minusminus(i,ll),Gamma_minusminus_N,Gamma_minusminus_U)
             rate_scatt_minusminus(i,ll)=Gamma_minusminus
             rate_scatt_minusminus_N(i,ll)=Gamma_minusminus_N
             rate_scatt_minusminus_U(i,ll)=Gamma_minusminus_U
          endif
 
        endif
     end do
 
     rate_scatt_4ph=rate_scatt_plusplus+rate_scatt_plusminus+rate_scatt_minusminus
   end subroutine RTA_driver_4ph
 
 end module processes