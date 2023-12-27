! Subroutines and functions related to Grüneisen parameters.
module gruneisen
  use misc
  use data
  use config
  use input
  implicit none

contains

  ! Subroutine to compute the mode Grüneisen parameters.
  subroutine mode_grun(omega,eigenvect,Ntri,Phi,&
       R_j,R_k,Index_i,Index_j,Index_k,grun)
    implicit none

    integer(kind=4),intent(in) :: Ntri,Index_i(Ntri),Index_j(Ntri),Index_k(Ntri)
    real(kind=8),intent(in) :: omega(nptk,nbands)
    real(kind=8),intent(in) :: Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri)
    complex(kind=8),intent(in) :: eigenvect(nptk,nbands,nbands)
    real(kind=8),intent(out) :: grun(nptk,nbands)

    real(kind=8),parameter :: unitfactor=9.6472d4 ! From nm*eV/(amu*A^3*THz^2) to 1.

    integer(kind=4) :: ik,ii,jj,kk,iband,itri,ialpha,ibeta
    real(kind=8) :: kspace(nptk,3)
    complex(kind=8) :: factor1,factor2,factor3,g(nbands)

    do ii=1,Ngrid(1)
       do jj=1,Ngrid(2)
          do kk=1,Ngrid(3)
             ik=((kk-1)*Ngrid(2)+(jj-1))*Ngrid(1)+ii
             kspace(ik,:)=rlattvec(:,1)*(ii-1.0)/ngrid(1)+&
                  rlattvec(:,2)*(jj-1.0)/ngrid(2)+&
                  rlattvec(:,3)*(kk-1.0)/ngrid(3)
          end do
       end do
    end do    
    do ik=1,nptk
       g=0.0
       do iband=1,nbands
          do itri=1,Ntri
             factor1=phexp(dot_product(kspace(ik,:),R_j(:,itri)))/&
                  sqrt(masses(types(Index_i(itri)))*&
                  masses(types(Index_j(itri))))
             do ialpha=1,3
                factor2=factor1*conjg(eigenvect(ik,iband,&
                     3*(Index_i(itri)-1)+ialpha))
                do ibeta=1,3
                   factor3=factor2*eigenvect(ik,iband,&
                        3*(Index_j(itri)-1)+ibeta)
                   g(iband)=g(iband)+factor3*dot_product(&
                        Phi(ialpha,ibeta,:,itri),&
                        (cartesian(:,Index_k(itri))+R_k(:,itri)))
                end do
             end do
          end do
       if (omega(ik,iband).eq.0) then
          grun(ik,iband)=0.d0
       else
          g(iband)=-unitfactor*g(iband)/6.d00/omega(ik,iband)**2
          grun(ik,iband)=real(g(iband))
       endif
       end do
    end do
  end subroutine mode_grun

  ! Obtain the total Grüneisen parameter as a weighted sum over modes.
  function total_grun(omega,grun)
    implicit none
    real(kind=8),intent(in) :: omega(nptk,nbands),grun(nptk,nbands)

    integer(kind=4) :: ii,jj
    real(kind=8) :: total_grun,weight,dBE,x

    total_grun=0.
    weight=0.
    do jj=1,nbands
       do ii=1,nptk
          x=hbar*omega(ii,jj)/(2.*kB*T)
          if(x.gt.1e-6) then
             dBE=(x/sinh(x))**2.
             weight=weight+dBE
             total_grun=total_grun+dBE*grun(ii,jj)
          end if
       end do
    end do
    total_grun=total_grun/weight
  end function total_grun
end module gruneisen
