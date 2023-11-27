!  FourPhonon: An extension module to ShengBTE for computing four phonon anharmonicity
!  Copyright (C) 2021-2023 Zherui Han <zrhan@purdue.edu>
!  Copyright (C) 2021 Xiaolong Yang <xiaolongyang1990@gmail.com>
!  Copyright (C) 2021 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2021 Tianli Feng <Tianli.Feng2011@gmail.com>
!  Copyright (C) 2021-2023 Xiulin Ruan <ruan@purdue.edu>
!  Copyright (C) 2023 Ziqi Guo <gziqi@purdue.edu>
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

! Routines to read FORCE_CONSTANTS_2ND and FORCE_CONSTANTS_3RD.
! Routines to read FORCE_CONSTANTS_4TH.
module input
  use iso_fortran_env
  use config
  implicit none

  ! The file names are defined here and can be changed for custom builds.
  ! This is especially relevant is the format needs to be changed.
  character(len=*),parameter :: filename_2fc="FORCE_CONSTANTS_2ND"
  character(len=*),parameter :: filename_3fc="FORCE_CONSTANTS_3RD"
  character(len=*),parameter :: filename_4fc="FORCE_CONSTANTS_4TH"
  real(kind=8),parameter :: unitfactor=9648.53336213 ! from eV/(A^2*amu) to THz^2

contains

  ! Read FORCE_CONSTANTS_2ND.
  subroutine read2fc(fc)
    implicit none

    include "mpif.h"

    real(kind=8),allocatable,intent(out) :: fc(:,:,:,:,:,:,:)
    
    integer(kind=4) :: ntot,atom1,atom2,i,j,ip,ierr
    integer(kind=4) :: ix1,iy1,iz1,ix2,iy2,iz2,iatom1,iatom2
    real(kind=8) :: mm(natoms,natoms)

    do i=1,natoms
       mm(i,i)=masses(types(i))
       do j=i+1,natoms
          mm(i,j)=sqrt(masses(types(i))*masses(types(j)))
          mm(j,i)=mm(i,j)
       end do
    end do

    allocate(fc(natoms,3,scell(1),scell(2),scell(3),natoms,3))

    ! Phonopy's 2nd-order format is quite straightforward. Each file
    ! is essentially a sequence of 3x3 blocks, one for each pair of atoms
    ! in the supercell. A single header line contains the number of atoms.
    open(1,file=filename_2fc,status="old")
    read(1,*) ntot
    if(ntot.ne.scell(1)*scell(2)*scell(3)*natoms) then
       if(myid.eq.0)write(error_unit,*) "Error: wrong number of force constants for the specified scell"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    do i=1,ntot
       do j=1,ntot
          read(1,*) atom1,atom2
          call split_index(atom1,scell(1),scell(2),scell(3),&
               ix1,iy1,iz1,iatom1)
          call split_index(atom2,scell(1),scell(2),scell(3),&
               ix2,iy2,iz2,iatom2)
          if(ix1.eq.1.and.iy1.eq.1.and.iz1.eq.1) then
             do ip=1,3
                read(1,*) fc(iatom1,ip,ix2,iy2,iz2,iatom2,:)
             end do
          else
             do ip=1,3
                read(1,*)
             end do
          end if
       end do
    end do
    close(1)

    ! After reading the force constants, they are reduced using the
    ! tensor product of the square root of the atomic masses. It is these
    ! reduced constants that enter the expression of the dynamical matrix.
    do iatom1=1,natoms
       do iatom2=1,natoms
          fc(iatom1,:,:,:,:,iatom2,:)=&
               fc(iatom1,:,:,:,:,iatom2,:)/mm(iatom1,iatom2)
       end do
    end do
    fc=unitfactor*fc
  end subroutine read2fc

  ! Convert a supercell index of the kind used by Phonopy into
  ! a set of unit cell+atom indices.
  subroutine split_index(index,nx,ny,nz,ix,iy,iz,iatom)
    implicit none

    integer(kind=4),intent(in) :: index,nx,ny,nz
    integer(kind=4),intent(out) :: ix,iy,iz,iatom

    integer(kind=4) :: tmp1,tmp2

    call divmod(index-1,nx,tmp1,ix)
    call divmod(tmp1,ny,tmp2,iy)
    call divmod(tmp2,nz,iatom,iz)

    ix=ix+1
    iy=iy+1
    iz=iz+1
    iatom=iatom+1
  end subroutine split_index

  ! Read FORCE_CONSTANTS_3RD.
  subroutine read3fc(Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)
    implicit none
    integer(kind=4),intent(out) :: Ntri
    integer(kind=4),allocatable,intent(out) :: Index_i(:),Index_j(:),Index_k(:)
    real(kind=8),allocatable,intent(out) :: Phi(:,:,:,:),R_j(:,:),R_k(:,:)

    real(kind=8) :: tmp(3,3)
    integer(kind=4) :: ii,jj,ll,mm,nn,ltem,mtem,ntem,info,P(3)

    ! The file is in a simple sparse format, described in detail in
    ! the user documentation. See Doc/ShengBTE.pdf.
    open(1,file=filename_3fc,status="old")
    read(1,*) Ntri
    allocate(Index_i(Ntri),Index_j(Ntri),Index_k(Ntri))
    allocate(Phi(3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri))
    do ii=1,Ntri
       read(1,*) jj
       read(1,*) R_j(:,ii)
       read(1,*) R_k(:,ii)
       read(1,*) Index_i(ii),Index_j(ii),Index_k(ii)
       do ll=1,3
          do mm=1,3
             do nn=1,3
                read(1,*) ltem,mtem,ntem,Phi(ll,mm,nn,ii)
             end do
          end do
       end do
    end do
    close(1)
    ! Each vector is rounded to the nearest lattice vector.
    tmp=lattvec
    call dgesv(3,Ntri,tmp,3,P,R_j,3,info)
    R_j=matmul(lattvec,anint(R_j/10.))
    tmp=lattvec
    call dgesv(3,Ntri,tmp,3,P,R_k,3,info)
    R_k=matmul(lattvec,anint(R_k/10.))
  end subroutine read3fc

  ! Read FORCE_CONSTANTS_4TH.
  subroutine read4fc(Ntri,Phi,R_j,R_k,R_l,Index_i,Index_j,Index_k,Index_l)
      implicit none
      integer(kind=4),intent(out) :: Ntri
      integer(kind=4),allocatable,intent(out) :: Index_i(:),Index_j(:),Index_k(:),Index_l(:)
      real(kind=8),allocatable,intent(out) :: Phi(:,:,:,:,:),R_j(:,:),R_k(:,:),R_l(:,:)

      real(kind=8) :: tmp(3,3)
      integer(kind=4) :: ii,jj,ll,mm,nn,pp,ltem,mtem,ntem,ptem,info,P(3)

      ! The file format is similar to 3rd-ifcs and should be self-documented.
      open(1,file=filename_4fc,status="old")
      read(1,*) Ntri
      allocate(Index_i(Ntri),Index_j(Ntri),Index_k(Ntri),Index_l(Ntri))
      allocate(Phi(3,3,3,3,Ntri),R_j(3,Ntri),R_k(3,Ntri),R_l(3,Ntri))
      do ii=1,Ntri
         read(1,*) jj
         read(1,*) R_j(:,ii)
         read(1,*) R_k(:,ii)
         read(1,*) R_l(:,ii)
         read(1,*) Index_i(ii),Index_j(ii),Index_k(ii),Index_l(ii)
         do ll=1,3
            do mm=1,3
               do nn=1,3
                  do pp=1,3
                     read(1,*) ltem,mtem,ntem,ptem,Phi(ll,mm,nn,pp,ii)
                  end do
               end do
            end do
         end do
      end do
      close(1)
      ! Each vector is rounded to the nearest lattice vector.
      tmp=lattvec
      call dgesv(3,Ntri,tmp,3,P,R_j,3,info)
      R_j=matmul(lattvec,anint(R_j/10.))
      tmp=lattvec
      call dgesv(3,Ntri,tmp,3,P,R_k,3,info)
      R_k=matmul(lattvec,anint(R_k/10.))
      tmp=lattvec
      call dgesv(3,Ntri,tmp,3,P,R_l,3,info)
      R_l=matmul(lattvec,anint(R_l/10.))
  end subroutine read4fc

end module input
