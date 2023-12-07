! Zherui HAN <zrhan@purdue.edu>
! 05/20/2021; Last edited: 12/07/2023
! Convert TDEP IFC-2nd to Phonopy format
! How to run: this program takes five arguments
! 1. filename of TDEP IFC-2nd 
! 2-4. supercell size/dimension
! 5. number of atoms in the unit cell
! An example: ./tdep2phonopy outfile.force_constant 10 10 1 2
! converts TDEP IFCs to a 10 10 1 supercell of a two-atom unitcell

! GNU GENERAL PUBLIC LICENSE; similar to FourPhonon

!> Main program
PROGRAM tdep2phonopy
    implicit none

    !input
    integer :: qx,qy,qz,flag=0
    integer :: uc_phonopy,ss_phonopy
    
    real(kind=8) :: ifc_tdep(3,3)=0.d0
    
    integer :: i,j,k,m,n,num_arg,temp,omitline
    integer :: atomtdep
    integer :: ix1,iy1,iz1,ix2,iy2,iz2,iatom1,iatom2 ! index of supercell and atom
    integer :: celltdep(3) ! relative index in TDEP
    real(kind=8),allocatable :: ifc_phonopy(:,:,:,:)
    integer, allocatable :: neighbor(:)
    
    CHARACTER(LEN=36) :: filename,inputname,chartemp

    CALL GET_COMMAND_ARGUMENT(1,inputname)
    ! dimension of supercell qx qy qz
    CALL GET_COMMAND_ARGUMENT(2,chartemp)
    read (chartemp,*) qx
    CALL GET_COMMAND_ARGUMENT(3,chartemp)
    read (chartemp,*) qy
    CALL GET_COMMAND_ARGUMENT(4,chartemp)
    read (chartemp,*) qz
    CALL GET_COMMAND_ARGUMENT(5,chartemp) ! number of atoms in the unitcell
    read (chartemp,*) uc_phonopy
    filename = TRIM(ADJUSTL(inputname))
    ss_phonopy = uc_phonopy*qx*qy*qz ! number of atoms in the supercell
    allocate(ifc_phonopy(ss_phonopy,ss_phonopy,3,3))
    ifc_phonopy = 0.d0 ! all elements are initialized to zero
    
    OPEN (2,FILE=filename,STATUS='old')
    read (2,*) temp
    if ( temp .ne. uc_phonopy ) then
        write(*,*) "Error: inconsistent number of atoms in the unitcell!"
        stop 0
    end if
    read (2,*)
    allocate(neighbor(uc_phonopy))
    ! get number of neighbors in TDEP
    do i = 1, uc_phonopy
        read (2,*) neighbor(i)
        do j = 1, 5*neighbor(i)
            read (2,*)
        end do
    end do
    close(2)

    do i = 1, ss_phonopy ! atom1
        do j = 1, ss_phonopy ! atom2
            call split_index(i,qx,qy,qz,ix1,iy1,iz1,iatom1)
            call split_index(j,qx,qy,qz,ix2,iy2,iz2,iatom2)
            ! now we are tring to find if there is a match in TDEP IFCs
            OPEN (2,FILE=filename,STATUS='old')
            ! first discard useless information and atoms
            read(2,*)
            read(2,*)
            if ( iatom1 .ne. 1 ) then
                do m = 1, iatom1-1
                    do n = 1, (5*neighbor(m)+1)
                        read (2,*) ! omitted lines
                    end do
                end do
            end if
            ! look for useful pairs
            read (2,*) ! number of neighbors that atom1 has
            do m = 1, neighbor(iatom1)
                read (2,*) atomtdep ! index in TDEP
                if ( atomtdep == iatom2 ) then
                    read (2,*) celltdep(:)
                    ! if it satisfies the index
                    if ( celltdep(1)==(ix2-ix1).and.celltdep(2)==(iy2-iy1)&
                        .and.celltdep(3)==(iz2-iz1)) then
                        do k = 1, 3
                            read (2,*) ifc_phonopy(i,j,k,:)
                        end do
                        exit 
                    else ! discard this block
                        do n = 1, 3
                            read (2,*)
                        end do
                    end if
                else ! discard this block
                    do n = 1, 4
                        read (2,*)
                    end do
                end if
            end do
            close (2)
        end do  
    end do

    CALL export_ifc(ifc_phonopy,ss_phonopy)


END PROGRAM tdep2phonopy

!> Subroutines

! Convert a supercell index of the kind used by Phonopy into
! a set of unit cell+atom indices. This is extracted from ShengBTE
subroutine split_index(index,nx,ny,nz,ix,iy,iz,iatom)
    implicit none

    integer(kind=4),intent(in) :: index,nx,ny,nz
    integer(kind=4),intent(out) :: ix,iy,iz,iatom

    integer(kind=4) :: tmp1,tmp2

    call divmod(index-1,nx,tmp1,ix)
    call divmod(tmp1,ny,tmp2,iy)
    call divmod(tmp2,nz,iatom,iz)

    ! now there are index of supercell ix iy iz
    ix=ix+1
    iy=iy+1
    iz=iz+1
    ! atom index
    iatom=iatom+1

    ! there is something we need to modify: periodic boundary conditions
    ! this operation would make ix/y/z go back to +-nxyz/2
    if ( ix .GT. ceiling(nx/2.) ) then
        ix = ix - nx
    end if
    if ( iy .GT. ceiling(ny/2.) ) then
        iy = iy - ny
    end if
    if ( iz .GT. ceiling(nz/2.) ) then
        iz = iz - nz
    end if
end subroutine split_index

! Subroutine to write out phonopy format force constant
subroutine export_ifc(ifc_phonopy,ss_phonopy)
    implicit NONE

    integer, intent(in) :: ss_phonopy
    real(kind=8),intent(in) :: ifc_phonopy(ss_phonopy,ss_phonopy,3,3)
    integer :: i,j,m,n

    OPEN (1,FILE='FORCE_CONSTANTS_2ND',STATUS='replace')
    write(1,'(1x,I8,1x,I8)')ss_phonopy,ss_phonopy

    do i = 1, ss_phonopy
        do j = 1, ss_phonopy
            write(1,'(1x,I8,1x,I8)') i,j
            write(1,"(<3>(' ',f21.15,''))") ((ifc_phonopy(i,j,m,n),n=1,3),m=1,3)
        end do
    
    end do
end subroutine

! Quotient and remainder of integer division.
subroutine divmod(a,b,q,r)
    implicit none

    integer(kind=4),intent(in) :: a,b
    integer(kind=4),intent(out) :: q,r

    q=a/b
    r=mod(a,b)
end subroutine divmod
