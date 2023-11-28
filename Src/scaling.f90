!  FourPhonon: An extension module to ShengBTE for computing four phonon anharmonicity
!  Copyright (C) 2021-2023 Zherui Han <zrhan@purdue.edu>
!  Copyright (C) 2021 Xiaolong Yang <xiaolongyang1990@gmail.com>
!  Copyright (C) 2021 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2021 Tianli Feng <Tianli.Feng2011@gmail.com>
!  Copyright (C) 2021-2023 Xiulin Ruan <ruan@purdue.edu>
!  Copyright (C) 2023 Ziqi Guo <gziqi@purdue.edu>
!  Copyright (C) 2023 Guang Lin <guanglin@purdue.edu>
!
!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2017 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2017 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2017 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2017 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Correction to the bulk tau required for nanowires.

module scaling
  use config
  use data
  implicit none

contains
  ! Correction to tau for nanowires, based on a radial average. Note
  ! that a direction needs to be selected before calling this
  ! subroutine: velocity_z is the projection of velocity over that
  ! axis.
  subroutine ScalingOfTau(Nlist,Nequi,ALLEquiList,velocity_z,&
       velocity,tauzero_wedge,radnw,ffunc)
    implicit none
    integer(kind=4),intent(in) :: Nlist,Nequi(nptk),ALLEquiList(Nsymm_rot,nptk)
    real(kind=8),intent(in) :: velocity(nptk,Nbands,3),velocity_z(nptk,Nbands)
    real(kind=8),intent(in) :: tauzero_wedge(Nbands,Nlist),radnw
    real(kind=8),intent(out) :: ffunc(nptk,Nbands)

    real(kind=8) :: vrho,xsum,tauint,vmod,vaxial,tauzero(Nbands,nptk)
    integer(kind=4),parameter :: nsum=500
    integer(kind=4) :: iisum,ii,jj

    real(kind=8) :: dnrm2

    do ii=1,Nlist
      do jj=1,Nequi(ii)
          tauzero(:,ALLEquiList(jj,ii))=tauzero_wedge(:,ii)
       end do
    end do

    ffunc=0.d0
    do ii=1,nptk
       do jj=1,Nbands
          vmod=dnrm2(3,velocity(ii,jj,:),1)
          tauint=abs(tauzero(jj,ii))
          vaxial=velocity_z(ii,jj)
          if ((vmod**2-vaxial**2).le.vaxial**2*1.0d-8) then
             vaxial=vaxial*2.d0*sqrt(dble(cgrid**2))/sqrt(4.d0*(cgrid**2)+1.d0)
          end if
          vrho=tauint*sqrt(vmod**2-vaxial**2)
          if ((vrho.ne.0.)) then
             ffunc(ii,jj)=0.0d0
             do iisum=1,nsum
                xsum=dfloat(iisum)*radnw/dfloat(nsum)
                ffunc(ii,jj)=ffunc(ii,jj)+2.0d0*sqrt(abs(radnw**2-xsum**2))/vrho-&
                     1.0d0+exp(-2.0d0*dsqrt(abs(radnw**2-xsum**2))/vrho)
             end do
             ffunc(ii,jj)=(vrho*2.0d0/(pi*radnw*nsum))*ffunc(ii,jj)
          end if
       end do
    end do
  end subroutine ScalingOfTau
end module scaling
