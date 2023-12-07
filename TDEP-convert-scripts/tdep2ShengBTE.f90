! Zherui HAN <zrhan@purdue.edu>
! 06/17/2021; Last edited: 12/07/2023
! Convert TDEP IFC-3rd/4th to ShengBTE format
! Input argument: 1 order of IFCs, 2 number of atoms in
! the unit cell, 3 filename to be converted
! We also need infile.ucposcar from TDEP at hand

! An example: ./tdep2ShengBTE 3 2 outfile.forceconstant_thirdorder
! converts TDEP 3rd-IFCs of a two-atom unitcell to ShengBTE format

! GNU GENERAL PUBLIC LICENSE; similar to FourPhonon


program tdep2ShengBTE
    implicit none

    integer          :: order, ucatom
    character(len=36):: filename, chartemp, inputname
    real(kind=8)     :: lfactor, r(3,3)

    ! command line read in
    CALL GET_COMMAND_ARGUMENT(1,chartemp)
    read (chartemp,*) order ! 3 or 4
    CALL GET_COMMAND_ARGUMENT(2,chartemp)
    read (chartemp,*) ucatom ! number of atoms in unit cell
    CALL GET_COMMAND_ARGUMENT(3,inputname)
    filename = TRIM(ADJUSTL(inputname))

    ! get the lattice vector
    open(01,file='infile.ucposcar',status='old')
    write(*,*) "Found infile.ucposcar"
    read(01,*)
    read(01,*) lfactor
    read(01,*) r(1,:)
    read(01,*) r(2,:)
    read(01,*) r(3,:)
    close(01)
    r = r*lfactor

    if ( order .eq. 3 ) then
        call threeconvert(ucatom,filename,r)
    else if ( order .eq. 4) then
        call fourconvert(ucatom,filename,r)
    else 
        write(*,*) "Error: order should be 3 or 4!"
        stop 0
    end if

end program tdep2ShengBTE

! 3rd-IFC conversion
subroutine threeconvert(ucatom,filename,r)
    implicit none 

    integer,intent(in)      :: ucatom
    character(len=36),intent(in)    :: filename
    real(kind=8),intent(in) :: r(3,3)

    integer,allocatable     :: neighbor(:)
    real(kind=8)            :: ifc_tdep(9,3), rtdep(2,3), cart_r(2,3)
    integer                 :: nblock, maxnb, indexblock, indexneighbor, indexatom, triplet(3), i, j, k

    allocate(neighbor(ucatom))
    ! read and let's first scan the file to get number of blocks
    open(2,file=filename,status='old')
    write(*,*) "Found 3rd-IFCs file", filename
    read(2,*)
    read(2,*)
    do i = 1, ucatom
        read(2,*) neighbor(i)
        do j = 1, neighbor(i)*15
            read(2,*)
        end do
    end do
    close(2)
    nblock = sum(neighbor)
    maxnb = maxval(neighbor)
    indexblock = 0
    ifc_tdep = 0.0d0
    open(2,file=filename,status='old')
    read(2,*)
    read(2,*)
    open(3,file="FORCE_CONSTANTS_3RD_tdep",status="replace")
    write(3,*) nblock
    write(3,*) ! a blank line
    do indexatom = 1, ucatom
        read(2,*)
        do indexneighbor = 1, neighbor(indexatom)
            read(2,*) triplet(1)
            read(2,*) triplet(2)
            read(2,*) triplet(3)
            read(2,*) 
            read(2,*) rtdep(1,:)
            read(2,*) rtdep(2,:)
            cart_r(1,:) = matmul(rtdep(1,:),r)
            cart_r(2,:) = matmul(rtdep(2,:),r)
            do i = 1,9
                read(2,*) ifc_tdep(i,:)
            end do
            indexblock = indexblock+1
            write(3,*) indexblock
            write(3,*) cart_r(1,:)
            write(3,*) cart_r(2,:)
            write(3,*) triplet(:)
            do i = 1,3
                do j = 1,3
                    do k = 1,3
                        write(3,"(3I2,f21.15)") i,j,k,ifc_tdep(3*i+j-3,k)                        
                    end do                   
                end do               
            end do
            write(3,*) ! a blank line
        end do        
    end do ! numberate ending
    if ( indexblock .eq. nblock ) then
        write(*,*) "That's a match in number of blocks!"
    else
        write(*,*) "Check the number of blocks!"
    end if



end subroutine threeconvert

! 4th-IFC conversion
subroutine fourconvert(ucatom,filename,r)
    implicit none 

    integer,intent(in)      :: ucatom
    character(len=36),intent(in)    :: filename
    real(kind=8),intent(in) :: r(3,3)

    integer,allocatable     :: neighbor(:)
    real(kind=8)            :: ifc_tdep(27,3), rtdep(3,3), cart_r(3,3)
    integer                 :: nblock, maxnb, indexblock, indexneighbor, indexatom, triplet(4), i, j, k, l

    allocate(neighbor(ucatom))
    ! read and let's first scan the file to get number of blocks
    open(2,file=filename,status='old')
    write(*,*) "Found 4th-IFCs file", filename
    read(2,*)
    read(2,*)
    do i = 1, ucatom
        read(2,*) neighbor(i)
        do j = 1, neighbor(i)*35
            read(2,*)
        end do
    end do
    close(2)
    nblock = sum(neighbor)
    maxnb = maxval(neighbor)
    indexblock = 0
    ifc_tdep = 0.0d0
    open(2,file=filename,status='old')
    read(2,*)
    read(2,*)
    open(3,file="FORCE_CONSTANTS_4TH_tdep",status="replace")
    write(3,*) nblock
    write(3,*) ! a blank line
    do indexatom = 1, ucatom
        read(2,*)
        do indexneighbor = 1, neighbor(indexatom)
            read(2,*) triplet(1)
            read(2,*) triplet(2)
            read(2,*) triplet(3)
            read(2,*) triplet(4)
            read(2,*) 
            read(2,*) rtdep(1,:)
            read(2,*) rtdep(2,:)
            read(2,*) rtdep(3,:)
            cart_r(1,:) = matmul(rtdep(1,:),r)
            cart_r(2,:) = matmul(rtdep(2,:),r)
            cart_r(3,:) = matmul(rtdep(3,:),r)
            do i = 1,27
                read(2,*) ifc_tdep(i,:)
            end do
            indexblock = indexblock+1
            write(3,*) indexblock
            write(3,*) cart_r(1,:)
            write(3,*) cart_r(2,:)
            write(3,*) cart_r(3,:)
            write(3,*) triplet(:)
            do i = 1,3
                do j = 1,3
                    do k = 1,3
                        do l = 1,3
                            write(3,"(4I2,f21.15)") i,j,k,l,ifc_tdep(9*i+3*j+k-12,l)  
                        end do                      
                    end do                   
                end do               
            end do
            write(3,*) ! a blank line
        end do        
    end do ! numberate ending
    if ( indexblock .eq. nblock ) then
        write(*,*) "That's a match in number of blocks!"
    else
        write(*,*) "Check the number of blocks!"
    end if



end subroutine fourconvert