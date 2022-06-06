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

! Lattice specific heat and small-grain-limit reduced thermal
! conductivity. Both are calculated as integrals over the BZ and contain
! no anharmonic information.
module integrals
  use config
  use data
  implicit none

contains

  ! Lattice specific heat per unit volume.
  function cv(omega)
    implicit none
    real(kind=8),intent(in) :: omega(nptk,nbands)

    integer(kind=4) :: ii,jj
    real(kind=8) :: cv,dBE,x

    cv=0.
    do jj=1,nbands
       do ii=1,nptk
          x=hbar*omega(ii,jj)/(2.*kB*T)
          if(x.eq.0.) then
             dBE=1.
          else
             dBE=(x/sinh(x))**2.
          end if
          cv=cv+dBE
       end do
    end do
    cv=kB*cv/(1e-27*V*nptk) ! J/(K m^3)
  end function cv

  ! Small-grain-limit reduced thermal conductivity. In contrast with cv,
  ! group velocities play a role in this integral.
  Subroutine kappasg(omega,velocity,velocity_offdiag,nruter,nruter_coh)
    implicit none
    real(kind=8),intent(in) :: omega(nptk,nbands),velocity(nptk,nbands,3),velocity_offdiag(nptk,nbands,nbands,3)
    real(kind=8),intent(out) :: nruter(3,3),nruter_coh(3,3)

    integer(kind=4) :: ii,jj,kk,dir1,dir2
    real(kind=8) :: x,x_coh,dBE,dBE_coh,tmp,tmp_coh

    real(kind=8) :: dnrm2

    nruter=0.
    do jj=1,nbands
       do ii=2,nptk
          x=hbar*omega(ii,jj)/(2.*kB*T)
          dBE=(x/sinh(x))**2.
          tmp=dnrm2(3,velocity(ii,jj,:),1)
          if(tmp.lt.1d-12) then
             cycle
          else
             dBE=dBE/tmp
          end if
          do dir1=1,3
             do dir2=1,3
                nruter(dir1,dir2)=nruter(dir1,dir2)+dBE*&
                     velocity(ii,jj,dir1)*velocity(ii,jj,dir2)
             end do
          end do
          ! Coherence term
          do kk=jj+1,nbands
             x_coh=hbar*(omega(ii,jj)+omega(ii,kk))/2/(2.*kB*T)
             dBE_coh=(x/sinh(x))**2.
             tmp_coh=dnrm2(3,velocity_offdiag(ii,jj,kk,:),1)
             if(tmp_coh.lt.1d-12) then
                cycle
             else
                dBE_coh=dBE_coh/tmp_coh
             end if
             do dir1=1,3
                do dir2=1,3
                   nruter_coh(dir1,dir2)=nruter_coh(dir1,dir2)+dBE_coh*&
                        velocity_offdiag(ii,jj,kk,dir1)*velocity_offdiag(ii,jj,kk,dir2)
                end do
             end do
          end do
       end do
    end do
    nruter=1e21*kB*nruter/(nptk*V)
    nruter_coh=1e21*kB*nruter_coh/(nptk*V)
  end subroutine kappasg
end module integrals
