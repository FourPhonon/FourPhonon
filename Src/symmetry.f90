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

! Thin, specialized wrapper around spglib, a library by Atsushi Togo.

module symmetry
  use iso_c_binding
  implicit none
  ! Tolerance parameter passed to spglib.
  real(kind=C_DOUBLE),parameter :: symprec=1d-5

  ! Explicit interfaces to spglib.
  interface
     function spg_get_symmetry(rotations,translations,nops,lattice,&
          positions,types,natoms,symprec) bind (C, name="spg_get_symmetry")
       use iso_c_binding
       integer(kind=C_INT) :: spg_get_symmetry
       integer(kind=C_INT),dimension(3,3,nops) :: rotations
       real(kind=C_DOUBLE),dimension(3,nops) :: translations
       integer(kind=C_INT),value :: nops
       real(kind=C_DOUBLE),dimension(3,3) :: lattice
       real(kind=C_DOUBLE),dimension(3,natoms) :: positions
       integer(kind=C_INT),dimension(natoms) :: types
       integer(kind=C_INT),value :: natoms
       real(kind=C_DOUBLE),value :: symprec
     end function spg_get_symmetry
     function spg_get_international(symbol,lattice,&
          positions,types,natoms,symprec) bind (C, name="spg_get_international")
       use iso_c_binding
       integer(kind=C_INT) :: spg_get_international
       character(kind=C_CHAR),dimension(11) :: symbol
       real(kind=C_DOUBLE),dimension(3,3) :: lattice
       real(kind=C_DOUBLE),dimension(3,natoms) :: positions
       integer(kind=C_INT),dimension(natoms) :: types
       integer(kind=C_INT),value :: natoms
       real(kind=C_DOUBLE),value :: symprec
     end function spg_get_international
     function spg_get_multiplicity(lattice,&
          positions,types,natoms,symprec) bind (C, name="spg_get_multiplicity")
       use iso_c_binding
       integer(kind=C_INT) :: spg_get_multiplicity
       real(kind=C_DOUBLE),dimension(3,3) :: lattice
       real(kind=C_DOUBLE),dimension(3,natoms) :: positions
       integer(kind=C_INT),dimension(natoms) :: types
       integer(kind=C_INT),value :: natoms
       real(kind=C_DOUBLE),value :: symprec
     end function spg_get_multiplicity
  end interface
contains
  
  ! Return the number of symmetry operations. Useful for allocating
  ! memory for get_operations().
  function get_num_operations(lattice,natoms,types,positions)
    implicit none
    real(kind=8),dimension(3,3),intent(in) :: lattice
    integer(kind=4),intent(in) :: natoms
    integer(kind=4),dimension(natoms),intent(in) :: types
    real(kind=8),dimension(3,natoms),intent(in) :: positions

    integer(kind=4) :: get_num_operations

    ! Notice the explicit C-compatible types used through this module.
    real(kind=C_DOUBLE),dimension(3,3) :: clattice
    integer(kind=C_INT) :: cnatoms
    integer(kind=C_INT),dimension(natoms) :: ctypes
    real(kind=C_DOUBLE),dimension(3,natoms) :: cpositions
    integer(kind=C_INT) :: num

    ! This kind of transposition is needed for interoperability with
    ! C.
    clattice=transpose(lattice)
    cnatoms=natoms
    ctypes=types
    cpositions=positions

    num=spg_get_multiplicity(clattice,cpositions,ctypes,&
         cnatoms,symprec)
    get_num_operations=num
  end function get_num_operations

  ! Return the matrix and vector representations of the symmetry
  ! operations of the system.
  subroutine get_operations(lattice,natoms,types,positions,nops,&
       rotations,translations,international)
    implicit none
    real(kind=8),dimension(3,3),intent(in) :: lattice
    integer(kind=4),intent(in) :: natoms
    integer(kind=4),dimension(natoms),intent(in) :: types
    real(kind=8),dimension(3,natoms),intent(in) :: positions
    integer(kind=4),intent(inout) :: nops
    integer(kind=4),dimension(3,3,nops),intent(out) :: rotations
    real(kind=8),dimension(3,nops),intent(out) :: translations
    character(len=10),intent(out) :: international

    integer(kind=C_INT) :: i
    integer(kind=C_INT) ::  newnops
    real(kind=C_DOUBLE),dimension(3,3) :: clattice
    integer(kind=C_INT) :: cnatoms
    integer(kind=C_INT),dimension(natoms) :: ctypes
    real(kind=C_DOUBLE),dimension(3,natoms) :: cpositions
    integer(kind=C_INT) :: cnops
    integer(kind=C_INT),dimension(3,3,nops) :: crotations
    real(kind=C_DOUBLE),dimension(3,nops) :: ctranslations
    character(len=11,kind=C_CHAR) :: intertmp

    clattice=transpose(lattice)
    cnatoms=natoms
    ctypes=types
    cpositions=positions
    cnops=nops

    ! If nops changes value, something went wrong. Checking
    ! this condition is up to the user.
    newnops = spg_get_symmetry(crotations,ctranslations,cnops,&
         clattice,cpositions,ctypes,cnatoms,symprec)
    i=spg_get_international(intertmp,clattice, cpositions,&
         ctypes,cnatoms,symprec)
    international=intertmp(1:10)
    nops=newnops
    ! Transform from C to Fortran order.
    do i=1,nops
       rotations(:,:,i)=transpose(crotations(:,:,i))
    end do
    translations=ctranslations
  end subroutine get_operations

  ! Return the Cartesian components of the rotations and translations
  ! returned by get_operations().
  subroutine get_cartesian_operations(lattice,nops,rotations,translations,&
       crotations,ctranslations)
    implicit none
    real(kind=8),dimension(3,3),intent(in) :: lattice
    integer(kind=4),intent(in) :: nops
    integer(kind=4),dimension(3,3,nops),intent(in) :: rotations
    real(kind=8),dimension(3,nops),intent(in) :: translations
    real(kind=8),dimension(3,3,nops),intent(out) :: crotations
    real(kind=8),dimension(3,nops),intent(out) :: ctranslations

    integer(kind=4) :: i,info
    integer(kind=4),dimension(3) :: P
    real(kind=8),dimension(3,3) :: tmp1,tmp2

    ctranslations=matmul(lattice,translations)
    do i=1,nops
       tmp1=transpose(lattice)
       tmp2=transpose(matmul(lattice,rotations(:,:,i)))
       ! Rotations transform as tensors: both the lattice-vector matrix
       ! and its inverse are needed. Explicit inversions are avoided.
       call dgesv(3,3,tmp1,3,P,tmp2,3,info)
       crotations(:,:,i)=transpose(tmp2)
    end do
  end subroutine get_cartesian_operations
end module symmetry
