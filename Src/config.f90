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

! Variables and routines used to read and store the configuration.
module config
  use iso_fortran_env
  use misc
  use data
  use symmetry
  implicit none

  integer(kind=4) :: nelements,natoms,ngrid(3),norientations
  namelist /allocations/ nelements,natoms,ngrid,norientations
  real(kind=8) :: lfactor,lattvec(3,3),epsilon(3,3)
  character(len=3),allocatable :: elements(:)
  integer(kind=4),allocatable :: types(:),orientations(:,:)
  integer(kind=4) :: scell(3)
  real(kind=8),allocatable :: positions(:,:),masses(:),gfactors(:),born(:,:,:)
  namelist /crystal/ lfactor,lattvec,elements,types,positions,masses,gfactors,&
       epsilon,born,scell,orientations
  integer(kind=4) :: maxiter,nticks
  real(kind=8) :: T,scalebroad,rmin,rmax,dr,eps
  real(kind=8) :: T_min,T_max,T_step,omega_max,Length
   ! ----------- sampling method add -----------
  integer(kind=4) :: num_sample_process_3ph,num_sample_process_3ph_phase_space, num_sample_process_4ph, num_sample_process_4ph_phase_space
  namelist /parameters/ T,scalebroad,rmin,rmax,dr,maxiter,nticks,eps,&
           T_min,T_max,T_step,omega_max,Length, &
           num_sample_process_3ph,num_sample_process_3ph_phase_space, num_sample_process_4ph,num_sample_process_4ph_phase_space
   ! ----------- end sampling method add -----------

  logical :: nonanalytic,convergence,isotopes,autoisotopes,nanowires,onlyharmonic,espresso,normal,umklapp,&
             four_phonon,four_phonon_iteration,nanolength
  namelist /flags/ nonanalytic,convergence,isotopes,autoisotopes,&
       nanowires,onlyharmonic,espresso,normal,umklapp,four_phonon,four_phonon_iteration,nanolength

  integer(kind=4) :: nbands,nptk,nwires
  real(kind=8) :: cgrid,V,rV,rlattvec(3,3),slattvec(3,3)
  real(kind=8),allocatable :: cartesian(:,:),uorientations(:,:)

  integer(kind=4) :: nsymm,nsymm_rot
  integer(kind=4),allocatable :: rotations(:,:,:)
  integer(kind=4),allocatable :: rotations_orig(:,:,:)
  real(kind=8),allocatable :: crotations(:,:,:),qrotations(:,:,:)
  real(kind=8),allocatable :: crotations_orig(:,:,:)
  real(kind=8),allocatable :: qrotations_orig(:,:,:)
  real(kind=8),allocatable :: symmetrizers(:,:,:)
  character(len=10) :: international

  ! MPI variables, assigned in ShengBTE.f90.
  integer(kind=4) :: myid,numprocs

contains

  subroutine read_config()
    implicit none

    include "mpif.h"

    integer(kind=4) :: i,j,k,ii,jj,kk,ll,info,ierr
    integer(kind=4) :: P(3)
    integer(kind=4),allocatable :: rtmp(:,:,:),ID_equi(:,:)
    logical,allocatable :: valid(:)
    real(kind=8) :: dnrm2,tmp1(3,3),tmp2(3,3),tmp3(3,3)
    real(kind=8),allocatable :: crtmp(:,:,:),qrtmp(:,:,:)
    real(kind=8),allocatable :: translations(:,:),ctranslations(:,:)

    ! Set the defaults and read the namelist.
    nelements=0
    natoms=0
    ngrid=0
    norientations=0
    open(1,file="CONTROL",status="old")
    read(1,nml=allocations)
    if(norientations.lt.0) then
       if(myid.eq.0)write(error_unit,*) "Error: norientations must be >=0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    if(nelements.lt.1.or.natoms.lt.1.or.natoms.lt.nelements) then
       if(myid.eq.0)write(error_unit,*) "Error: nelements,natoms must be >0, natoms must be >=nelements"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    if(any(ngrid.lt.1)) then
       if(myid.eq.0)write(error_unit,*) "Error: all components of ngrid must be must be >0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    allocate(elements(nelements),types(natoms),positions(3,natoms),&
         masses(nelements),gfactors(nelements),born(3,3,natoms),&
         cartesian(3,natoms))
    if(norientations.ne.0) then
       allocate(orientations(3,norientations))
       allocate(uorientations(3,norientations))
       orientations=0
    end if
    lfactor=1.
    scell=-1
    epsilon=0.
    epsilon(1,1)=1.
    epsilon(2,2)=1.
    epsilon(3,3)=1.
    born=0.
    gfactors=0.
    types=0
    read(1,nml=crystal)
    if(.not.all(scell.gt.0)) then
       if(myid.eq.0)write(error_unit,*) "Error: all supercell sizes must be >0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    if(.not.all(types.ne.0)) then
       if(myid.eq.0)write(error_unit,*) "Error: atom types must be initialized correctly"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    do i=1,norientations
       if(all(orientations(:,i).eq.0)) then
          if(myid.eq.0)write(error_unit,*) "Error: orientation uninitialized or zero"
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_FINALIZE(ierr)
       end if
    end do
    T=0
    T_min=0
    scalebroad=1.0
    rmin=5.0
    rmax=505.0
    dr=100.0
    maxiter=1000
    nticks=100
    eps=1e-5
    omega_max=1.d100
    Length=1
   ! ----------- sampling method add -----------
    num_sample_process_3ph = -1
    num_sample_process_3ph_phase_space = -1
    num_sample_process_4ph = -1
    num_sample_process_4ph_phase_space = -1
   ! ----------- end sampling method add -----------

    read(1,nml=parameters)
    if ((T.le.0.).and.(T_min.le.0)) then
       if(myid.eq.0)write(error_unit,*) "Error: T must be >0 K"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    elseif (T.gt.0) then
    T_min=T
    T_max=T
    T_step=T
    end if
    if(rmin.le.0.or.rmax.le.rmin.or.dr.le.0) then
       if(myid.eq.0)write(error_unit,*) "Error: rmin and dr must be >0, and rmax must be > rmin"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    if(maxiter.le.0) then
       if(myid.eq.0)write(error_unit,*) "Error: maxiter must be >0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
    nonanalytic=.true.
    convergence=.true.
    isotopes=.true.
    autoisotopes=.true.
    nanowires=.false.
    onlyharmonic=.false.
    espresso=.false.
    
    ! Four-phonon namelist
    four_phonon=.false.
    four_phonon_iteration=.false.
    read(1,nml=flags)
    if(four_phonon_iteration.and.convergence.eq. .false.) then
      if(myid.eq.0)write(error_unit,*) "Error: four_phonon_iteration=.TRUE. but convergence=.FALSE."
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
    end if
    if(four_phonon_iteration.and.four_phonon.eq. .false.) then
      if(myid.eq.0)write(error_unit,*) "Error: four_phonon_iteration=.TRUE. but four_phonon=.FALSE."
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
    end if
    if(nanowires.and.norientations.eq.0) then
       if(myid.eq.0)write(error_unit,*) "Error: nanowires=.TRUE. but norientations=0"
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_FINALIZE(ierr)
    end if
   ! ----------- sampling method add -----------
    if(((four_phonon.eq. .false.) .and. (num_sample_process_4ph.gt.0)) .or.((four_phonon.eq. .false.) .and. (num_sample_process_4ph_phase_space.gt.0)) ) then
      if(myid.eq.0)write(error_unit,*) "Error: four_phonon=.false. but num_sample_process_4ph.gt.0 or num_sample_process_4ph_phase_space.gt.0. "
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
    end if

    if(((convergence.eq. .true.) .and. (num_sample_process_3ph.gt.0)) ) then
      if(myid.eq.0)write(error_unit,*) "Error: convergence=.true. but num_sample_process_3ph.gt.0. "
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
    end if

    if(((four_phonon_iteration.eq. .true.) .and. (num_sample_process_4ph.gt.0)) ) then
      if(myid.eq.0)write(error_unit,*) "Error: four_phonon_iteration=.true. but num_sample_process_4ph.gt.0. "
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
    end if
   ! ----------- end sampling method add -----------
    close(1)


    nptk=product(ngrid)
    cgrid=nptk**(1./3.)
    nbands=3*natoms
    nwires=ceiling((rmax-rmin)/dr)

    lattvec=lfactor*lattvec

    do i=1,3
       j=mod(i,3)+1
       k=mod(j,3)+1
       call cross_product(lattvec(:,j),lattvec(:,k),rlattvec(:,i))
    end do

    V=abs(dot_product(lattvec(:,1),rlattvec(:,1)))
    rV=2.*pi/V
    rlattvec=rV*rlattvec

    cartesian=matmul(lattvec,positions)

    if(nanowires) then
       uorientations=matmul(lattvec,orientations)
       do i=1,norientations
          uorientations(:,i)=uorientations(:,i)/&
               dnrm2(3,uorientations(:,i),1)
       end do
    end if

    ! Compute the average masses and g-factors automatically, if
    ! requested.
    if(autoisotopes) then
       call data_fill_isotopes()
       do i=1,nelements
          call data_calc_mandg(elements(i),masses(i),gfactors(i))
       end do
       call data_free_isotopes()
    end if

    ! Find out the symmetries of the system.
    nsymm=get_num_operations(lattvec,natoms,types,positions)
    ! We need do double that number to take time reversal symetry into
    ! account.
    nsymm_rot=2*nsymm
    allocate(translations(3,nsymm),rotations(3,3,nsymm_rot),&
         ctranslations(3,nsymm),crotations(3,3,nsymm_rot),&
         qrotations(3,3,nsymm_rot))
    allocate(rotations_orig(3,3,nsymm),&
         crotations_orig(3,3,nsymm),&
         qrotations_orig(3,3,nsymm),symmetrizers(3,3,nptk))
    call get_operations(lattvec,natoms,types,positions,nsymm,&
         rotations_orig,translations,international)
    rotations(:,:,1:nsymm)=rotations_orig
    if(myid.eq.0)write(*,*) "Info: symmetry group ",trim(international)," detected"
    if(myid.eq.0)write(*,*) "Info: ",nsymm," symmetry operations"
    call get_cartesian_operations(lattvec,nsymm,rotations_orig,translations,&
         crotations_orig,ctranslations)
    crotations(:,:,1:nsymm)=crotations_orig
    deallocate(translations,ctranslations)

    ! Transform the rotation matrices to the reciprocal-space basis.
    do i=1,nsymm
       tmp1=matmul(transpose(lattvec),lattvec)
       tmp2=transpose(rotations_orig(:,:,i))
       tmp3=tmp1
       call dgesv(3,3,tmp1,3,P,tmp2,3,info)
       qrotations_orig(:,:,i)=transpose(matmul(tmp2,tmp3))
    end do
    qrotations(:,:,1:nsymm)=qrotations_orig

    ! Fill the second half of the rotation matrix list
    ! using time reversal symmetry.
    rotations(:,:,nsymm+1:2*nsymm)=-rotations_orig(:,:,1:nsymm)
    qrotations(:,:,nsymm+1:2*nsymm)=-qrotations_orig(:,:,1:nsymm)
    crotations(:,:,nsymm+1:2*nsymm)=-crotations_orig(:,:,1:nsymm)

    allocate(ID_Equi(nsymm_rot,nptk))
    call symmetry_map_notransl(ID_equi)
    ! Create the "symmetrizers", linear operators that extract the
    ! component of a vector compatible with the symmetries at each q point.
    symmetrizers=0
    do ii=1,nptk
       kk=0
       do jj=1,nsymm_rot
          if(ID_equi(jj,ii).eq.ii) then
             symmetrizers(:,:,ii)=symmetrizers(:,:,ii)+&
                  crotations(:,:,jj)
             kk=kk+1
          end if
       end do
       if(kk.gt.1) then
          symmetrizers(:,:,ii)=symmetrizers(:,:,ii)/kk
       end if
    end do
    ! Find rotations that are either duplicated or incompatible with
    ! the q-point grid.
    call symmetry_map(ID_equi)
    allocate(valid(nsymm_rot))
    valid=.TRUE.
    jj=0
    do ii=1,nsymm_rot
       if(valid(ii).and.any(ID_equi(ii,:).eq.-1)) then
          valid(ii)=.FALSE.
          jj=jj+1
       end if
    end do
    if(myid.eq.0.and.jj.ne.0) then
       write(*,*) "Info:",jj,&
            "rotations are incompatible with the q-point grid and will be discarded"
    end if
    ll=0
    do ii=2,nsymm_rot
       do i=1,ii-1
          if(.not.valid(i))cycle
          if(all(rotations(:,:,ii).eq.rotations(:,:,i))) then
             valid(ii)=.FALSE.
             ll=ll+1
             exit
          end if
       end do
    end do
    if(myid.eq.0.and.ll.ne.0) then
       write(*,*) "Info:",ll,"duplicated rotations will be discarded"
    end if
    ! Filter out those rotations through a series of move_alloc calls.
    ! Arrays to take into account: rotations,crotations,qrotations.
    if(ll+jj.ne.0) then
       allocate(rtmp(3,3,nsymm_rot-ll-jj))
       allocate(crtmp(3,3,nsymm_rot-ll-jj))
       allocate(qrtmp(3,3,nsymm_rot-ll-jj))
       kk=0
       do ii=1,nsymm_rot
          if(valid(ii)) then
             kk=kk+1
             rtmp(:,:,kk)=rotations(:,:,ii)
             crtmp(:,:,kk)=crotations(:,:,ii)
             qrtmp(:,:,kk)=qrotations(:,:,ii)
          end if
       end do
       nsymm_rot=nsymm_rot-ll-jj
       call move_alloc(rtmp,rotations)
       call move_alloc(crtmp,crotations)
       call move_alloc(qrtmp,qrotations)
    end if
    deallocate(ID_Equi,valid)
  end subroutine read_config

  ! Free the memory used by all config structures.
  subroutine free_config()
    deallocate(elements,types,positions,masses,gfactors,born,cartesian,&
         rotations,crotations,qrotations,symmetrizers)
    if(nanowires) then
       deallocate(orientations,uorientations)
    end if
  end subroutine free_config

  ! Compute all images of a reciprocal-space point using the
  ! rotational part of the symmetry operations. Everything is
  ! performed in lattice coordinates.
  subroutine symm(r_in,r_out)
    implicit none
    integer(kind=4),intent(in) :: r_in(3)
    real(kind=8),intent(out) :: r_out(3,nsymm_rot)

    integer(kind=4) :: ii

    do ii=1,nsymm_rot
       r_out(:,ii)=ngrid*matmul(qrotations(:,:,ii),dble(r_in)/ngrid)
    end do
  end subroutine symm

  ! Find the equivalences among points.
  subroutine symmetry_map(ID_equi)
    implicit none
    integer(kind=4),intent(out) :: ID_equi(nsymm_rot,nptk)

    integer(kind=4) :: Ind_cell(3,nptk)
    integer(kind=4) :: i,isym,ivec(3)
    real(kind=8) :: vec(3),vec_symm(3,nsymm_rot),dnrm2

    call Id2Ind(Ind_cell)
    do i=1,nptk
       call symm(Ind_cell(:,i),vec_symm)
       do isym=1,nsymm_rot
          vec=vec_symm(:,isym)
          ivec=nint(vec)
          if(dnrm2(3,abs(vec-dble(ivec)),1).gt.1e-2) then
             ID_equi(isym,i)=-1
          else
             ID_equi(isym,i)=Ind2Id(modulo(ivec,ngrid))
          end if
       end do
    end do
  end subroutine symmetry_map

  ! Find symmetry operations excluding translations that can bring q to itself.
  subroutine symmetry_map_notransl(ID_equi)
    implicit none
    integer(kind=4),intent(out) :: ID_equi(nsymm_rot,nptk)

    integer(kind=4) :: Ind_cell(3,nptk)
    integer(kind=4) :: i,isym,ivec(3)
    real(kind=8) :: vec(3),vec_symm(3,nsymm_rot),dnrm2

    call Id2Ind(Ind_cell)
    do i=1,nptk
       call symm(Ind_cell(:,i),vec_symm)
       do isym=1,nsymm_rot
          vec=vec_symm(:,isym)
          ivec=nint(vec)
          if(dnrm2(3,abs(dble(Ind_cell(:,i))-dble(vec)),1).le.1e-5) then
             ID_equi(isym,i)=i
          else
             ID_equi(isym,i)=-1
          endif
       end do
    end do
  end subroutine symmetry_map_notransl

  ! Create a table that can be used to demultiplex cell indices.
  subroutine Id2ind(Ind_cell)
    implicit none
    integer(kind=4),intent(out) :: Ind_cell(3,nptk)

    integer(kind=4) :: ii,tmp1

    do ii=1,nptk
       call divmod(ii-1,Ngrid(1),tmp1,Ind_cell(1,ii))
       call divmod(tmp1,Ngrid(2),Ind_cell(3,ii),Ind_cell(2,ii))
    end do
  end subroutine Id2ind

  ! Multiplex three cell indices into one.
  function Ind2Id(Ind_cell)
    implicit none
    integer(kind=4),intent(in) :: Ind_cell(3)

    integer(kind=4) :: Ind2Id

    Ind2Id=1+Ind_cell(1)+(Ind_cell(2)+Ind_cell(3)*Ngrid(2))*Ngrid(1)
  end function Ind2Id

  ! Return the base broadening (without prefactor) for a mode.
  function base_sigma(v)
    implicit none
    real(kind=8),intent(in) :: v(3)

    real(kind=8) :: base_sigma

    integer(kind=4) :: nu

    base_sigma=0.
   !  do nu=1,3
   !    base_sigma=base_sigma+(dot_product(rlattvec(:,nu),v)/ngrid(nu))**2
   !  end do

   !  base_sigma=sqrt(base_sigma/6.)

    base_sigma=DMAX1((dot_product(rlattvec(:,1),v)/ngrid(1))**2,(dot_product(rlattvec(:,2),v)/ngrid(2))**2,(dot_product(rlattvec(:,3),v)/ngrid(3))**2)
    base_sigma=sqrt(base_sigma/2.)

  end function base_sigma


  ! Force a real 3x3 Cartesian tensor to fulfill all the symmetries.
  subroutine symmetrize_tensor(tensor)
    implicit none

    real(kind=8),intent(inout) :: tensor(3,3)

    integer(kind=4) :: isym
    real(kind=8) :: tmp(3,3)

    tmp = 0.
    do isym=1,nsymm_rot
       tmp = tmp + matmul(crotations(:,:,isym),&
            matmul(tensor,transpose(crotations(:,:,isym))))
    end do

    tensor=tmp/nsymm_rot
  end subroutine symmetrize_tensor
end module config
