! Compute the thermal conductivity and its cumulative version.
module conductivity
  use config
  use data
  implicit none

contains

  ! Straightforward implementation of the thermal conductivity as an integral
  ! over the whole Brillouin zone in terms of frequencies, velocities and F_n.
  subroutine TConduct(omega,velocity,F_n,ThConductivity,ThConductivityMode)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands,3),F_n(Nbands,nptk,3)
    real(kind=8),intent(out) :: ThConductivity(Nbands,3,3)
    real(kind=8),intent(out) :: ThConductivityMode(nptk,Nbands,3,3)

    real(kind=8) :: fBE,tmp(3,3)
    integer(kind=4) :: ii,jj,dir1,dir2

    ThConductivity=0.d0
    ThConductivityMode=0.d0
    do jj=1,Nbands
       do ii=2,nptk
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(ii,jj,dir1)*F_n(jj,ii,dir2)
             end do
          end do
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          ThConductivityMode(ii,jj,:,:)=fBE*(fBE+1)*omega(ii,jj)*tmp
          ThConductivity(jj,:,:)=ThConductivity(jj,:,:)+ThConductivityMode(ii,jj,:,:)
       end do
    end do
    ThConductivity=1e21*hbar**2*ThConductivity/(kB*T*T*V*nptk)
    ThConductivityMode=1e21*hbar**2*ThConductivityMode/(kB*T*T*V*nptk)
  end subroutine TConduct

  ! Specialized version of the above subroutine for those cases where kappa
  ! can be treated as a scalar.
  subroutine TConductScalar(omega,velocity,F_n,ThConductivity)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands),F_n(Nbands,nptk)
    real(kind=8),intent(out) :: ThConductivity(Nbands)

    real(kind=8) :: fBE
    integer(kind=4) :: ii,jj

    ThConductivity=0.d0
    do jj=1,Nbands
       do ii=2,nptk
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          ThConductivity(jj)=ThConductivity(jj)+fBE*(fBE+1)*omega(ii,jj)*&
               velocity(ii,jj)*F_n(jj,ii)
       end do
    end do
    ThConductivity=1e21*hbar**2*ThConductivity/(kB*T*T*V*nptk)
  end subroutine TConductScalar

  ! "Cumulative thermal conductivity": value of kappa obtained when
  ! only phonon with mean free paths below a threshold are considered.
  ! ticks is a list of the thresholds to be employed.
  subroutine CumulativeTConduct(omega,velocity,F_n,ticks,results)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands,3),F_n(Nbands,nptk,3)
    real(kind=8),intent(out) :: ticks(nticks),results(Nbands,3,3,Nticks)

    real(kind=8) :: fBE,tmp(3,3),lambda
    integer(kind=4) :: ii,jj,kk,dir1,dir2

    real(kind=8) :: dnrm2

    do ii=1,nticks
       ticks(ii)=10.**(8.*(ii-1.)/(nticks-1.)-2.)
    end do

    results=0.
    do jj=1,Nbands
       do ii=2,nptk
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(ii,jj,dir1)*F_n(jj,ii,dir2)
             end do
          end do
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          tmp=fBE*(fBE+1)*omega(ii,jj)*tmp
          lambda=dot_product(F_n(jj,ii,:),velocity(ii,jj,:))/(&
               (omega(ii,jj)+1d-12)*dnrm2(3,velocity(ii,jj,:),1))
          do kk=1,nticks
             if(ticks(kk).gt.lambda) then
                results(jj,:,:,kk)=results(jj,:,:,kk)+tmp
             end if
          end do
       end do
    end do
    results=1e21*hbar**2*results/(kB*T*T*V*nptk)
  end subroutine CumulativeTConduct

  ! "Cumulative thermal conductivity vs angular frequency": value of kappa obtained when
  ! only frequencies below a threshold are considered.
  ! ticks is a list of the thresholds to be employed.
  subroutine CumulativeTConductOmega(omega,velocity,F_n,ticks,results)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands,3),F_n(Nbands,nptk,3)
    real(kind=8),intent(out) :: ticks(nticks),results(Nbands,3,3,Nticks)

    real(kind=8) :: fBE,tmp(3,3),lambda
    integer(kind=4) :: ii,jj,kk,dir1,dir2

    REAL(kind=8)  EMIN,EMAX

    EMIN=1.d10
    EMAX=-1.d10
    DO ii=1,NPTK
    DO jj=1,Nbands
       emin = min(emin, omega(ii,jj))
       emax = max(emax, omega(ii,jj))
    ENDDO
    ENDDO
    emax=emax*1.1d0

    do ii=1,nticks
       ticks(ii)=emin+(emax-emin)*ii/nticks
    end do

    results=0.
    do jj=1,Nbands
       do ii=2,nptk
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(ii,jj,dir1)*F_n(jj,ii,dir2)
             end do
          end do
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          tmp=fBE*(fBE+1)*omega(ii,jj)*tmp
          lambda=omega(ii,jj)
          do kk=1,nticks
             if(ticks(kk).gt.lambda) then
                results(jj,:,:,kk)=results(jj,:,:,kk)+tmp
             end if
          end do
       end do
    end do
    results=1e21*hbar**2*results/(kB*T*T*V*nptk)
  end subroutine CumulativeTConductOmega
end module conductivity
