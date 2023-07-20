!  FourPhonon: An extension module to ShengBTE for computing four phonon anharmonicity
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

! Compute the thermal conductivity and its cumulative version.
module conductivity
  use config
  use data
  implicit none

contains

  ! Straightforward implementation of the thermal conductivity as an integral
  ! over the whole Brillouin zone in terms of frequencies, velocities and F_n.
  subroutine TConduct(omega,velocity,velocity_offdiag,F_n,Nlist,Nequi,ALLEquiList,ThConductivity,ThConductivityMode,ThConductivityCoh,ThConductivityCohMode,rate)
    implicit none

    integer(kind=4),intent(in) :: Nlist,Nequi(nptk),ALLEquiList(Nsymm_rot,nptk)
    real(kind=8),intent(in) :: omega(nptk,Nbands)
    real(kind=8),intent(in) :: velocity(nptk,Nbands,3),velocity_offdiag(nptk,Nbands,Nbands,3),F_n(Nbands,nptk,3)
    real(kind=8),intent(out) :: ThConductivity(Nbands,3,3), ThConductivityCoh(Nbands,Nbands,3,3)
    real(kind=8),intent(out) :: ThConductivityMode(nptk,Nbands,3,3),ThConductivityCohMode(nptk,Nbands,Nbands,3,3)
    real(kind=8),intent(out) :: rate(Nbands,nptk)

    integer(kind=4) :: ii,jj,kk,ll,dir1,dir2
    real(kind=8) :: fBE,fBE1,fBE2,tmp(3,3),tmp_coh(3,3)


    ThConductivity=0.d0
    ThConductivityMode=0.d0
    ThConductivityCoh=0.d0
    ThConductivityCohMode=0.d0
    rate=0.d-5

    ! Scattering rate
    do ii=1,nptk
       do jj=1,Nbands
          if (norm2(F_n(jj,ii,:)).ne.0 .and. omega(ii,jj).gt.0) then
             rate(jj,ii)=norm2(velocity(ii,jj,:))*omega(ii,jj)/norm2(F_n(jj,ii,:)) ! F_n has omega factor
          end if
       end do
    end do
    print *, "RATE avg. over q-point:", sum(rate,dim=2)/nptk

    ! Calculate thermal conductivity
    do jj=1,Nbands
       do ii=2,nptk
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(ii,jj,dir1)*F_n(jj,ii,dir2)
             end do
          end do
          ! The Wigner coherence term: Simoncelli, Marzari, & Mauri. Nature Physics 15:809-813 (2019)  
          do kk=jj,Nbands
             if (jj.eq.kk) then ! test jj.ne.kk to check off-diagonal formulation = diagonal formulation of normal BTE
                cycle
             end if
             do dir1=1,3
                do dir2=1,3
                   tmp_coh(dir1,dir2)=velocity_offdiag(ii,jj,kk,dir1)*velocity_offdiag(ii,kk,jj,dir2)
                end do
             end do
             fBE1=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
             fBE2=1.d0/(exp(hbar*omega(ii,kk)/Kb/T)-1.D0)
             ThConductivityCohMode(ii,jj,kk,:,:)=(fBE1*(fBE1+1)*omega(ii,jj)+fBE2*(fBE2+1)*omega(ii,kk))*tmp_coh
             ThConductivityCohMode(ii,jj,kk,:,:)=ThConductivityCohMode(ii,jj,kk,:,:)*(omega(ii,jj)+omega(ii,kk))/2*(rate(jj,ii)+rate(kk,ii))
             ThConductivityCohMode(ii,jj,kk,:,:)=ThConductivityCohMode(ii,jj,kk,:,:)/(4*(omega(ii,jj)-omega(ii,kk))**2+(rate(jj,ii)+rate(kk,ii))**2)
             ThConductivityCoh(jj,kk,:,:)=ThConductivityCoh(jj,kk,:,:)+ThConductivityCohMode(ii,jj,kk,:,:)
             if (maxval(ThConductivityCohMode(ii,jj,kk,:,:)).gt.10**5) then
                print *,ii,jj,kk,ThConductivityCohMode(ii,jj,kk,1,1),ThConductivityCohMode(ii,jj,kk,2,2), abs(omega(ii,jj)-omega(ii,kk)), rate(jj,ii)+rate(kk,ii)
             end if
             ThConductivityCohMode(ii,kk,jj,:,:)=ThConductivityCohMode(ii,jj,kk,:,:)
             ThConductivityCoh(kk,jj,:,:)=ThConductivityCoh(jj,kk,:,:)
          end do
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          ThConductivityMode(ii,jj,:,:)=fBE*(fBE+1)*omega(ii,jj)*tmp
          ThConductivity(jj,:,:)=ThConductivity(jj,:,:)+ThConductivityMode(ii,jj,:,:)
       end do
    end do
    do jj=1,Nbands
       print *, "Average for band ",jj,":", sum(sum(ThConductivityCohMode(:,jj,:,1,1),dim=1),dim=1)/nptk
    end do
    ThConductivity=1e21*hbar**2*ThConductivity/(kB*T*T*V*nptk)
    ThConductivityMode=1e21*hbar**2*ThConductivityMode/(kB*T*T*V*nptk)
    ThConductivityCoh=1e21*hbar**2*ThConductivityCoh/(kB*T*T*V*nptk)
    ThConductivityCohMode=1e21*hbar**2*ThConductivityCohMode/(kB*T*T*V*nptk)
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
  subroutine CumulativeTConduct(omega,rate,velocity,velocity_offdiag,F_n,ticks,results,results_coh)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),rate(Nbands,nptk),velocity(nptk,Nbands,3),velocity_offdiag(nptk,Nbands,Nbands,3),F_n(Nbands,nptk,3)
    real(kind=8),intent(out) :: ticks(nticks),results(Nbands,3,3,Nticks),results_coh(Nbands,3,3,Nticks)

    real(kind=8) :: fBE,fBE1,fBE2,tmp(3,3),tmp_coh(3,3),lambda,lambda_coh
    integer(kind=4) :: ii,jj,kk,mm,dir1,dir2

    real(kind=8) :: dnrm2

    do ii=1,nticks
       ticks(ii)=10.**(8.*(ii-1.)/(nticks-1.)-2.)
    end do

    results=0.
    results_coh=0.
    do jj=1,Nbands
       do ii=2,nptk
          lambda=dot_product(F_n(jj,ii,:),velocity(ii,jj,:))/(&
               (omega(ii,jj)+1d-12)*dnrm2(3,velocity(ii,jj,:),1))
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(ii,jj,dir1)*F_n(jj,ii,dir2)
             end do
          end do
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          tmp=fBE*(fBE+1)*omega(ii,jj)*tmp
          do kk=1,nticks
             if(ticks(kk).gt.lambda) then
                results(jj,:,:,kk)=results(jj,:,:,kk)+tmp
             end if
          end do
          ! The Wigner coherence term: Simoncelli, Marzari, & Mauri. Nature Physics 15:809-813 (2019)
          do mm=1,Nbands
             if (jj.eq.mm) then
                cycle
             end if
             lambda_coh=dot_product(velocity_offdiag(ii,jj,mm,:),velocity_offdiag(ii,mm,jj,:))/&
                  (rate(mm,ii)*dnrm2(3,velocity_offdiag(ii,jj,mm,:),1))
             do dir1=1,3
                do dir2=1,3
                   tmp_coh(dir1,dir2)=velocity_offdiag(ii,jj,mm,dir1)*velocity_offdiag(ii,mm,jj,dir2)
                end do
             end do
             fBE1=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
             fBE2=1.d0/(exp(hbar*omega(ii,mm)/Kb/T)-1.D0)
             tmp_coh=tmp_coh*(fBE1*(fBE1+1)*omega(ii,jj)+fBE2*(fBE2+1)*omega(ii,mm))
             tmp_coh=tmp_coh*(omega(ii,jj)+omega(ii,mm))/2*(rate(jj,ii)+rate(mm,ii))
             tmp_coh=tmp_coh/(4*(omega(ii,jj)-omega(ii,mm))**2+(rate(jj,ii)+rate(mm,ii))**2)
             do kk=1,nticks             
                if(ticks(kk).gt.lambda_coh) then
                   results_coh(jj,:,:,kk)=results_coh(jj,:,:,kk)+tmp_coh
                end if
             end do
          end do
       end do
    end do
    results=1e21*hbar**2*results/(kB*T*T*V*nptk)
    results_coh=1e21*hbar**2*results_coh/(kB*T*T*V*nptk)
  end subroutine CumulativeTConduct

  ! "Cumulative thermal conductivity vs angular frequency": value of kappa obtained when
  ! only frequencies below a threshold are considered.
  ! ticks is a list of the thresholds to be employed.
  subroutine CumulativeTConductOmega(omega,rate,velocity,velocity_offdiag,F_n,ticks,results,results_coh)
    implicit none

    real(kind=8),intent(in) :: omega(nptk,Nbands),rate(Nbands,nptk),velocity(nptk,Nbands,3),velocity_offdiag(nptk,Nbands,Nbands,3),F_n(Nbands,nptk,3)
    real(kind=8),intent(out) :: ticks(nticks),results(Nbands,3,3,Nticks),results_coh(Nbands,3,3,Nticks)

    real(kind=8) :: fBE,fBE1,fBE2,tmp(3,3),tmp_coh(3,3),lambda,lambda_coh
    integer(kind=4) :: ii,jj,kk,mm,dir1,dir2

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
    results_coh=0.
    do jj=1,Nbands
       do ii=2,nptk
          lambda=omega(ii,jj)
          do dir1=1,3
             do dir2=1,3
                tmp(dir1,dir2)=velocity(ii,jj,dir1)*F_n(jj,ii,dir2)
             end do
          end do
          fBE=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
          tmp=fBE*(fBE+1)*omega(ii,jj)*tmp
          do kk=1,nticks
             if(ticks(kk).gt.lambda) then
                results(jj,:,:,kk)=results(jj,:,:,kk)+tmp
             end if
          end do
          ! The Wigner coherence term: Simoncelli, Marzari, & Mauri. Nature Physics 15:809-813 (2019)  
          do mm=1,Nbands
             if (jj.eq.mm) then
                cycle
             end if
             lambda_coh=(omega(ii,jj)+omega(ii,mm))/2
             do dir1=1,3
                do dir2=1,3
                   tmp_coh(dir1,dir2)=velocity_offdiag(ii,jj,mm,dir1)*velocity_offdiag(ii,mm,jj,dir2)
                end do
             end do
             fBE1=1.d0/(exp(hbar*omega(ii,jj)/Kb/T)-1.D0)
             fBE2=1.d0/(exp(hbar*omega(ii,mm)/Kb/T)-1.D0)
             tmp_coh=tmp_coh*(fBE1*(fBE1+1)*omega(ii,jj)+fBE2*(fBE2+1)*omega(ii,mm))
             tmp_coh=tmp_coh*(omega(ii,jj)+omega(ii,mm))/2*(rate(jj,ii)+rate(mm,ii))
             tmp_coh=tmp_coh/(4*(omega(ii,jj)-omega(ii,mm))**2+(rate(jj,ii)+rate(mm,ii))**2)
             do kk=1,nticks
                if(ticks(kk).gt.lambda_coh) then
                   results_coh(jj,:,:,kk)=results_coh(jj,:,:,kk)+tmp_coh
                end if
             end do
          end do
       end do
    end do
    results=1e21*hbar**2*results/(kB*T*T*V*nptk)
    results_coh=1e21*hbar**2*results_coh/(kB*T*T*V*nptk)
  end subroutine CumulativeTConductOmega
end module conductivity
