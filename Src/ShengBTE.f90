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



! Main program.
program ShengBTE
  use misc
  use config
  use input
  use processes
  use wedgetc
  use iterations
  use conductivity
  use scaling
  use phonon_routines
  use dos_routines
  use gruneisen
  use integrals
  ! use omp_lib
  implicit none

  include "mpif.h"
  

  real(kind=8) :: kappa_sg(3,3),kappa_old(3,3),relchange
  integer(kind=4) :: i,j,ii,jj,kk,ll,mm
  integer(kind=4) :: Tcounter
  real(kind=8),allocatable :: energy(:,:),q0(:,:),q0_reduced(:,:),velocity(:,:,:),velocity_z(:,:)
  complex(kind=8),allocatable :: eigenvect(:,:,:)

  real(kind=8),allocatable :: grun(:,:)

  real(kind=8),allocatable :: rate_scatt(:,:),rate_scatt_plus(:,:),rate_scatt_minus(:,:)
  real(kind=8),allocatable :: rate_scatt_plus_N(:,:),rate_scatt_minus_N(:,:),rate_scatt_plus_U(:,:),rate_scatt_minus_U(:,:)
   
  ! Variables related to four-phonon scatterings
  real(kind=8),allocatable :: rate_scatt_4ph(:,:),rate_scatt_plusplus(:,:),rate_scatt_plusminus(:,:),rate_scatt_minusminus(:,:)
  real(kind=8),allocatable :: rate_scatt_plusplus_N(:,:),rate_scatt_plusminus_N(:,:),rate_scatt_minusminus_N(:,:)
  real(kind=8),allocatable :: rate_scatt_plusplus_U(:,:),rate_scatt_plusminus_U(:,:),rate_scatt_minusminus_U(:,:)
  integer(kind=4) :: Ntri_4fc
  real(kind=8),allocatable :: Psi(:,:,:,:,:),R_s(:,:),R_t(:,:),R_u(:,:)
  integer(kind=4),allocatable :: Index_r(:),Index_s(:),Index_t(:),Index_u(:)
  integer(kind=8) :: Ntotal_plusplus,Ntotal_plusminus,Ntotal_minusminus
  integer(kind=8),allocatable :: N_plusplus(:),N_plusminus(:),N_minusminus(:)
  integer(kind=8),allocatable :: Indof2ndPhonon_plusplus(:),Indof3rdPhonon_plusplus(:),Indof4thPhonon_plusplus(:)
  integer(kind=8),allocatable :: Indof2ndPhonon_plusminus(:),Indof3rdPhonon_plusminus(:),Indof4thPhonon_plusminus(:)
  integer(kind=8),allocatable :: Indof2ndPhonon_minusminus(:),Indof3rdPhonon_minusminus(:),Indof4thPhonon_minusminus(:)
  real(kind=8),allocatable :: Pspace_plusplus_total(:,:),Pspace_plusminus_total(:,:),Pspace_minusminus_total(:,:)
  real(kind=8),allocatable :: Gamma_plusplus(:),Gamma_plusminus(:),Gamma_minusminus(:)

  real(kind=8),allocatable :: tau_zero(:,:),tau(:,:),tau_b(:,:),tau2(:,:),tau_b2(:,:)
  real(kind=8),allocatable :: dos(:),pdos(:,:),rate_scatt_isotope(:,:)
  real(kind=8),allocatable :: F_n(:,:,:),F_n_0(:,:,:),F_n_aux(:,:)
  real(kind=8),allocatable :: ThConductivity(:,:,:)
  real(kind=8),allocatable :: ThConductivityMode(:,:,:,:)
  real(kind=8),allocatable :: ticks(:),cumulative_kappa(:,:,:,:)


  integer(kind=4) :: Ntri
  ! IJK are vectors of product(ngrid) points in the lattice coordinate
  ! system, with each component ranging from 0 to Ngrid(:)-1.  Index_N
  ! is a mapping of 3 indices for an individual q-point into one index.
  integer(kind=4),allocatable :: Index_i(:),Index_j(:),Index_k(:),IJK(:,:),Index_N(:,:,:)
  real(kind=8),allocatable :: Phi(:,:,:,:),R_j(:,:),R_k(:,:) 

  integer(kind=4) :: nlist,Ntotal_plus,Ntotal_minus
  integer(kind=4),allocatable :: nequi(:),list(:)
  integer(kind=4),allocatable :: AllEquiList(:,:),TypeofSymmetry(:,:),eqclasses(:)
  integer(kind=4),allocatable :: N_plus(:),N_minus(:)
  integer(kind=4),allocatable :: Indof2ndPhonon_plus(:),Indof3rdPhonon_plus(:)
  integer(kind=4),allocatable :: Indof2ndPhonon_minus(:),Indof3rdPhonon_minus(:)
  real(kind=8) :: radnw,kappa_or_old
  real(kind=8),allocatable :: Gamma_plus(:),Gamma_minus(:),Gamma_plus_N(:),Gamma_minus_N(:),Gamma_plus_U(:),Gamma_minus_U(:)
  real(kind=8),allocatable :: Pspace_plus_total(:,:)
  real(kind=8),allocatable :: Pspace_minus_total(:,:)
  real(kind=8),allocatable :: ffunc(:,:),radnw_range(:),v_or(:,:),F_or(:,:)
  real(kind=8),allocatable :: kappa_or(:),kappa_wires(:,:),kappa_wires_reduce(:,:)

  integer(kind=4) :: iorient,ierr,REQUIRED,PROVIDED
  character(len=4) :: aux,aux2
  character(len=1024) :: path
  character(len=128) :: sorientation
  character(len=30) :: tempstring1,tempstring2,tempstring3

  real(kind=8) :: dnrm2

  call MPI_INIT_thread(MPI_THREAD_FUNNELED,PROVIDED,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  if (provided < MPI_THREAD_FUNNELED) then
      write(*,*) "Insufficient threading support"
      call MPI_FINALIZE(ierr)
      stop
  endif

  ! Read CONTROL and initialize variables such as the
  ! direct/reciprocal lattice vectors, the symmetry operations, the
  ! Pearson deviation coefficients related to isotopic disorder and so
  ! on and so forth.
  call read_config()

  allocate(energy(nptk,nbands))
  allocate(eigenvect(nptk,nbands,nbands))
  allocate(grun(nptk,nbands))
  allocate(q0(nptk,3),q0_reduced(nptk,3))
  allocate(velocity(nptk,nbands,3))
  allocate(velocity_z(nptk,nbands))
  allocate(IJK(3,nptk))
  allocate(Index_N(0:(ngrid(1)-1),0:(ngrid(2)-1),0:(ngrid(3)-1)))
  allocate(F_n(nbands,nptk,3))
  allocate(F_n_0(nbands,nptk,3),F_n_aux(nbands,nptk))
  allocate(ThConductivity(nbands,3,3))
  allocate(ThConductivityMode(nptk,nbands,3,3))
  allocate(kappa_wires(nbands,nwires))
  allocate(kappa_wires_reduce(nbands,nwires))
  allocate(Nequi(nptk))
  allocate(List(nptk))
  allocate(AllEquiList(nsymm_rot,nptk))
  allocate(TypeOfSymmetry(nsymm_rot,nptk))
  allocate(eqclasses(nptk))
  allocate(ffunc(nptk,nbands))
  allocate(v_or(nptk,nbands))
  allocate(F_or(nbands,nptk))
  allocate(kappa_or(nbands))

  ! Obtain the q-point equivalence classes defined by symmetry
  ! operations.
  call Id2Ind(IJK)
  call wedge(Nlist,Nequi,List,ALLEquiList,TypeofSymmetry)
  do ll=1,Nlist
     do kk=1,Nequi(ll)
        eqclasses(ALLEquiList(kk,ll))=ll !List(ll)
     end do
  end do


  do ii=1,Ngrid(1)
     do jj=1,Ngrid(2)
        do kk=1,Ngrid(3)
           ll=((kk-1)*Ngrid(2)+(jj-1))*Ngrid(1)+ii
           q0(ll,:)=rlattvec(:,1)*(ii-1.d0)/ngrid(1)+&
                rlattvec(:,2)*(jj-1.d0)/ngrid(2)+&
                rlattvec(:,3)*(kk-1.d0)/ngrid(3)
           q0_reduced(ll,1)=(ii-1.d0)/ngrid(1)
           q0_reduced(ll,2)=(jj-1.d0)/ngrid(2)
           q0_reduced(ll,3)=(kk-1.d0)/ngrid(3)
        end do
     end do
  end do

  if(myid.eq.0) then
     open(1,file="BTE.ReciprocalLatticeVectors",status="replace")
     write(1,"(3F20.10,A12)") rlattvec(:,1), "# nm-1,b1"
     write(1,"(3F20.10,A12)") rlattvec(:,2), "# nm-1,b2"
     write(1,"(3F20.10,A12)") rlattvec(:,3), "# nm-1,b3"
     close(1)
     open(1,file="BTE.qpoints",status="replace")
     do ll=1,Nlist
        write(1,"(I9,x,I9,x,I9,x,3(E20.10,x))") ll, List(ll),Nequi(ll),q0_reduced(List(ll),:)
     end do
     close(1)
     open(1,file="BTE.qpoints_full",status="replace")
     do ll=1,nptk
        write(1,"(I9,x,I9,x,3(E20.10,x))") ll,eqclasses(ll),q0_reduced(ll,:)
     end do
     close(1)
  end if

  ! Obtain the phonon spectrum from the 2nd-order force constants (in
  ! either of the two formats suported) and the dielectric parameters
  ! in CONTROL.
  if(myid.eq.0) then
     write(*,*) "Info: about to obtain the spectrum"
     if(espresso) then
        write(*,*) "Info: expecting Quantum Espresso 2nd-order format"
     else
        write(*,*) "Info: expecting Phonopy 2nd-order format"
     end if
  end if
  call eigenDM(energy,eigenvect,velocity)
  if(myid.eq.0)write(*,*) "Info: spectrum calculation finished"

  open(101,file="BTE.cvVsT")
  open(102,file="BTE.KappaTensorVsT_sg")
  if (myid.eq.0) print*, "Info: start calculating specific heat and kappa in the small-grain limit "
  do Tcounter=1,ceiling((T_max-T_min)/T_step)+1
     T=T_min+(Tcounter-1)*T_step
     if ((T.gt.T_max).and.(T.lt.(T_max+1.d0)))  exit
     if (T.gt.(T_max+1.d0)) T=T_max
     if (myid.eq.0) then
        print *, "Info: Temperature=", T
        write(aux2,"(I0)") NINT(T)
        path="T"//trim(adjustl(aux2))//"K"
        call create_directory(trim(adjustl(path))//C_NULL_CHAR)
        call change_directory(trim(adjustl(path))//C_NULL_CHAR)
        ! Compute the harmonic integrals: lattice specific heat and
        ! small-grain-limit reduced thermal conductivity. Write out this
        ! information, as well as the spectrum itself.
        open(1,file="BTE.cv",status="replace")
        write(1,*) cv(energy)
        close(1)
        call kappasg(energy,velocity,kappa_sg)
        open(1,file="BTE.kappa_sg",status="replace")
        write(1,"(9E14.5)") kappa_sg
        close(1)
        write(101,"(F7.1,E14.5)") T,cv(energy)
        flush(101)
        write(102,"(F7.1,9E14.5)") T,kappa_sg
        flush(102)
        call change_directory(".."//C_NULL_CHAR)
     end if
  end do
  close(101)
  close(102)

  write(aux,"(I0)") 3*Nbands
  if(myid.eq.0) then
     open(1,file="BTE.v",status="replace")
     do ii=1,Nbands
        do ll=1,Nlist
           write(1,"(3E20.10)") velocity(list(ll),ii,:)
        end do
     end do
     close(1)
     open(1,file="BTE.v_full",status="replace")
     do ii=1,Nbands
        do ll=1,nptk
           write(1,"(3E20.10)") velocity(ll,ii,:)
        end do
     end do
  end if
  write(aux,"(I0)") Nbands
  if(myid.eq.0) then
     open(1,file="BTE.omega",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") energy(list(ll),:)
     end do
     close(1)
  end if

  ! Compute the normalized boundary scattering rates.
  allocate(tau_b(Nbands,Nlist))
  allocate(tau_b2(Nbands,nptk))
  do ll=1,Nlist
     do ii=1,Nbands
        tau_b(ii,ll)=1.d00/dnrm2(3,velocity(List(ll),ii,:),1)
     end do
  end do
  do ll=1,nptk
     do ii=1,Nbands
        tau_b2(ii,ll)=1.d00/dnrm2(3,velocity(ll,ii,:),1)
     end do
  end do
  if (myid.eq.0) then
     write(aux,"(I0)") Nbands
     open(2,file="BTE.w_boundary",status="replace")
     do i=1,Nbands
        do ll = 1,Nlist
           write(2,"(2E20.10)") energy(list(ll),i),  1./tau_b(i,ll)
        end do
     end do
     close(2)
  end if
  ! Locally adaptive estimates of the total and projected densities of states.
  allocate(ticks(nticks))
  allocate(dos(nticks))
  allocate(pdos(nticks,natoms))
  call calc_dos(energy,velocity,eigenvect,ticks,dos,pdos)
  if(myid.eq.0) then
     open(1,file="BTE.dos",status="replace")
     do mm=1,nticks
        write(1,"(2E14.5)") ticks(mm),dos(mm)
     end do
     close(1)
     write(aux,"(I0)") Natoms
     open(1,file="BTE.pdos",status="replace")
     do mm=1,nticks
        write(1,"(E14.5,"//trim(adjustl(aux))//"E14.5)")&
             ticks(mm),pdos(mm,:)
     end do
     close(1)
  end if
  deallocate(ticks)


  allocate(rate_scatt_isotope(Nbands,Nlist))
  if (isotopes) then
     call calc_isotopescatt(energy,velocity,eigenvect,nlist,list,&
          rate_scatt_isotope)
    else
     rate_scatt_isotope=0.d00
  end if

  if(myid.eq.0) then
     write(aux,"(I0)") Nbands
     open(1,file="BTE.w_isotopic",status="replace")
     do i=1,Nbands
        do ll=1,Nlist
           write(1,"(2E20.10)") energy(list(ll),i),rate_scatt_isotope(i,ll)
        end do
     end do
     close(1)
  end if
 

  allocate(rate_scatt(Nbands,Nlist))
  allocate(rate_scatt_plus(Nbands,Nlist))
  allocate(rate_scatt_plus_N(Nbands,Nlist))
  allocate(rate_scatt_plus_U(Nbands,Nlist))

  allocate(rate_scatt_minus(Nbands,Nlist))
  allocate(rate_scatt_minus_N(Nbands,Nlist))
  allocate(rate_scatt_minus_U(Nbands,Nlist))

  allocate(rate_scatt_4ph(Nbands,Nlist))
  allocate(rate_scatt_plusplus(Nbands,Nlist))
  allocate(rate_scatt_plusplus_N(Nbands,Nlist))
  allocate(rate_scatt_plusplus_U(Nbands,Nlist))
  allocate(rate_scatt_plusminus(Nbands,Nlist))
  allocate(rate_scatt_plusminus_N(Nbands,Nlist))
  allocate(rate_scatt_plusminus_U(Nbands,Nlist))
  allocate(rate_scatt_minusminus(Nbands,Nlist))
  allocate(rate_scatt_minusminus_N(Nbands,Nlist))
  allocate(rate_scatt_minusminus_U(Nbands,Nlist))



  allocate(tau_zero(Nbands,Nlist))
  allocate(tau(Nbands,Nlist))
  allocate(tau2(Nbands,nptk))

  rate_scatt=0.d0
  rate_scatt_plus=0.d0
  rate_scatt_plus_N=0.d0
  rate_scatt_plus_U=0.d0
  rate_scatt_minus=0.d0
  rate_scatt_minus_N=0.d0
  rate_scatt_minus_U=0.d0
  rate_scatt_4ph=0.d00
  rate_scatt_plusplus=0.d00
  rate_scatt_plusplus_N=0.d00
  rate_scatt_plusplus_U=0.d00
  rate_scatt_plusminus=0.d00
  rate_scatt_plusminus_N=0.d00
  rate_scatt_plusminus_U=0.d00
  rate_scatt_minusminus=0.d00
  rate_scatt_minusminus_N=0.d00
  rate_scatt_minusminus_U=0.d00


  allocate(radnw_range(nwires))
  do ii=1,nwires
     radnw_range(ii)=rmin+(ii-1.0)*dr
  end do

  ! N_plus: number of allowed absorption processes
  ! N_minus: number of allowed emission processes.
  allocate(N_plus(Nlist*Nbands))
  allocate(N_minus(Nlist*Nbands))
    ! N_plusplus: number of allowed p1+p2+p3 => p4 processes.
    ! N_plusminus: number of allowed p1+p2-p3 => p4 processes.
    ! N_minusminus: number of allowed p1-p2-p3 => p4 processes.
    if (four_phonon) then
      allocate(N_plusplus(Nlist*Nbands))
      allocate(N_plusminus(Nlist*Nbands))
      allocate(N_minusminus(Nlist*Nbands))
      allocate(Pspace_plusplus_total(Nbands,Nlist))
      allocate(Pspace_plusminus_total(Nbands,Nlist))
      allocate(Pspace_minusminus_total(Nbands,Nlist))
    end if
  ! Phase space volume per mode and their sum.
  allocate(Pspace_plus_total(Nbands,Nlist))
  allocate(Pspace_minus_total(Nbands,Nlist))


  call NP_driver(energy,velocity,Nlist,List,IJK,&
       N_plus,Pspace_plus_total,N_minus,Pspace_minus_total)

  Ntotal_plus=sum(N_plus)
  Ntotal_minus=sum(N_minus)
  
  if(myid.eq.0)write(*,*) "Info: Ntotal_plus =",Ntotal_plus
  if(myid.eq.0)write(*,*) "Info: Ntotal_minus =",Ntotal_minus
  if (four_phonon) then
    if(myid.eq.0.and.four_phonon)write(*,*) "Info: Start calculating 4ph allowed processes"
    call NP_driver_4ph(energy,velocity,Nlist,List,IJK,&
         N_plusplus,Pspace_plusplus_total,N_plusminus,&
         Pspace_plusminus_total,N_minusminus,Pspace_minusminus_total)
    Ntotal_plusplus=sum(N_plusplus)
    Ntotal_plusminus=sum(N_plusminus)
    Ntotal_minusminus=sum(N_minusminus)
  end if
  if(myid.eq.0.and.four_phonon) write(*,*) "Info: Ntotal_plusplus,Ntotal_plusminus,Ntotal_minusminus="
  if(myid.eq.0.and.four_phonon) write(*,"(1X,I0,1X,I0,1X,I0)") Ntotal_plusplus,Ntotal_plusminus,Ntotal_minusminus
  if(myid.eq.0.and.four_phonon) then
    open(1,file="BTE.Numprocess_4ph",status="replace")
    do i=1,Nbands
       do ll=1,Nlist
          write(1,"(1X,I0,1X,I0,1X,I0)") N_plusplus(ll+Nlist*(i-1)),N_plusminus(ll+Nlist*(i-1)),N_minusminus(ll+Nlist*(i-1))
       end do
    end do
    close(1)
  end if
  
  if(myid.eq.0) then
     open(1,file="BTE.P3_plus",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_plus_total(:,ll)
     end do
     close(1)
     open(1,file="BTE.P3_minus",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_minus_total(:,ll)
     end do
     close(1)
     open(1,file="BTE.P3",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") 2.*(Pspace_plus_total(:,ll)+&
             Pspace_minus_total(:,ll)/2.)/3.
     end do
     close(1)
     do ii=1,Nlist
        Pspace_plus_total(:,ii)=Pspace_plus_total(:,ii)*Nequi(ii)
     end do
     open(1,file="BTE.P3_plus_total",status="replace")
     write(1,*) sum(Pspace_plus_total)
     close(1)
     do ii=1,Nlist
        Pspace_minus_total(:,ii)=Pspace_minus_total(:,ii)*Nequi(ii)
     end do
     open(1,file="BTE.P3_minus_total",status="replace")
     write(1,*) sum(Pspace_minus_total)
     close(1)
     open(1,file="BTE.P3_total",status="replace")
     write(1,*) 2.*(sum(Pspace_plus_total)+sum(Pspace_minus_total)/2.)/3.
     close(1)
     if (four_phonon) then
        open(1,file="BTE.P4_plusplus",status="replace")
        do ll=1,Nlist
            write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_plusplus_total(:,ll)
        end do
        close(1)
        open(1,file="BTE.P4_plusmius",status="replace")
        do ll=1,Nlist
            write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_plusminus_total(:,ll)
        end do
        close(1)
        open(1,file="BTE.P4_minusminus",status="replace")
        do ll=1,Nlist
            write(1,"("//trim(adjustl(aux))//"E20.10)") Pspace_minusminus_total(:,ll)
        end do
        close(1)
        open(1,file="BTE.P4",status="replace")
        do ll=1,Nlist
            write(1,"("//trim(adjustl(aux))//"E20.10)") 6.*(Pspace_plusplus_total(:,ll)+&
                Pspace_plusminus_total(:,ll)+Pspace_minusminus_total(:,ll))/7.
        end do
        close(1)
        do ii=1,Nlist
            Pspace_plusplus_total(:,ii)=Pspace_plusplus_total(:,ii)*Nequi(ii)
        end do
        open(1,file="BTE.P4_plusplus_total",status="replace")
        write(1,*) sum(Pspace_plusplus_total)
        close(1)
        do ii=1,Nlist
          Pspace_plusminus_total(:,ii)=Pspace_plusminus_total(:,ii)*Nequi(ii)
        end do
        open(1,file="BTE.P4_plusminus_total",status="replace")
        write(1,*) sum(Pspace_plusminus_total)
        close(1)
        do ii=1,Nlist
            Pspace_minusminus_total(:,ii)=Pspace_minusminus_total(:,ii)*Nequi(ii)
        end do
        open(1,file="BTE.P4_minusminus_total",status="replace")
        write(1,*) sum(Pspace_minusminus_total)
        close(1)
        open(1,file="BTE.P4_total",status="replace")
        write(1,*) 6.*(sum(Pspace_plusplus_total)+sum(Pspace_plusminus_total)+&
                    sum(Pspace_minusminus_total))/7.
        close(1)
     end if
  end if

  if(onlyharmonic) then
      ! weighted phase space (WP3/WP4) is calculated here if onlyharmonic=.true., 
      ! otherwise WP3/WP4 will be calculated later together with BTE.w_anharmonic(BTE.w_3ph/4ph)
      if (myid.eq.0) print*, "Info: start calculating weighted phase space"
      do Tcounter=1,CEILING((T_max-T_min)/T_step)+1
         T=T_min+(Tcounter-1)*T_step
         if ((T.gt.T_max).and.(T.lt.(T_max+1.d0)))  exit
         if (T.gt.(T_max+1.d0)) T=T_max
         if (myid.eq.0) then
            print *, "Info: Temperature=", T
            write(aux2,"(I0)") NINT(T)
            path="T"//trim(adjustl(aux2))//"K"
            call change_directory(trim(adjustl(path))//C_NULL_CHAR)
         end if
         call RTA_driver(energy,velocity,eigenvect,Nlist,List,IJK,&
               Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,rate_scatt,&
               rate_scatt_plus,rate_scatt_minus,Pspace_plus_total,Pspace_minus_total,&
               rate_scatt_plus_N,rate_scatt_minus_N,rate_scatt_plus_U,rate_scatt_minus_U)
         if (four_phonon) then ! four-phonon function
            call RTA_driver_4ph(energy,velocity,eigenvect,Nlist,List,IJK,&
                  Ntri_4fc,Psi,R_s,R_t,R_u,Index_r,Index_s,Index_t,Index_u,rate_scatt_4ph,&
                  rate_scatt_plusplus,rate_scatt_plusminus,rate_scatt_minusminus,&
                  Pspace_plusplus_total,Pspace_plusminus_total,Pspace_minusminus_total,&
                  rate_scatt_plusplus_N,rate_scatt_plusminus_N,rate_scatt_minusminus_N,&
                  rate_scatt_plusplus_U,rate_scatt_plusminus_U,rate_scatt_minusminus_U)
         end if
         if(myid.eq.0) then
            open(1,file="BTE.WP3_plus",status="replace")
            do i=1,Nbands
               do ll=1,Nlist
                  write(1,"(6E14.5)") energy(list(ll),i),Pspace_plus_total(i,ll)                
               end do
            end do
            close(1)
            open(1,file="BTE.WP3_minus",status="replace")
            do i=1,Nbands
               do ll=1,Nlist
                  write(1,"(6E14.5)") energy(list(ll),i),Pspace_minus_total(i,ll)
               end do
            end do
            close(1)
            open(1,file="BTE.WP3",status="replace")
            do i=1,Nbands
               do ll=1,Nlist
                  write(1,"(2E14.5)") energy(list(ll),i),&
                        Pspace_plus_total(i,ll)+Pspace_minus_total(i,ll)
               end do
            end do
            close(1)
            if (four_phonon) then
              open(1,file="BTE.WP4_plusplus",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(6E14.5)") energy(list(ll),i),Pspace_plusplus_total(i,ll)
                end do
              end do
              close(1)
              open(1,file="BTE.WP4_plusminus",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(6E14.5)") energy(list(ll),i),Pspace_plusminus_total(i,ll)
                end do
              end do
              close(1)
              open(1,file="BTE.WP4_minusminus",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(6E14.5)") energy(list(ll),i),Pspace_minusminus_total(i,ll)
                end do
              end do
              close(1)
              open(1,file="BTE.WP4",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(2E14.5)") energy(list(ll),i),Pspace_plusplus_total(i,ll)+&
                                        Pspace_plusminus_total(i,ll)+Pspace_minusminus_total(i,ll)
                end do
              end do
              close(1)
            end if
         end if
      end do
  endif  ! onlyharmonic
  deallocate(Pspace_plus_total)
  deallocate(Pspace_minus_total)
  if (four_phonon) then
    deallocate(Pspace_plusplus_total)
    deallocate(Pspace_plusminus_total)
    deallocate(Pspace_minusminus_total)
  end if


  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! Up to this point, no anharmonic information is used.
  if(onlyharmonic) then
     write(*,*) "Info: onlyharmonic=.true., stopping here"
     call MPI_FINALIZE(ierr)
     stop
  end if

  ! Load the 3rd-IFCs from FORCE_CONSTANTS_3RD.
  call read3fc(Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k)

  ! Load the 4th-IFCs from FORCE_CONSTANTS_4TH.
  if (four_phonon) then
      call read4fc(Ntri_4fc,Psi,R_s,R_t,R_u,Index_r,Index_s,Index_t,Index_u)
  end if

  call mode_grun(energy,eigenvect,Ntri,Phi,R_j,R_k,&
       Index_i,Index_j,Index_k,grun)
  write(aux,"(I0)") Nbands
  if(myid.eq.0) then
     open(1,file="BTE.gruneisen",status="replace")
     do ll=1,Nlist
        write(1,"("//trim(adjustl(aux))//"E20.10)") grun(list(ll),:)
     end do
     close(1)
  end if
  open(101,file="BTE.gruneisenVsT_total")
  do Tcounter=1,CEILING((T_max-T_min)/T_step)+1
     T=T_min+(Tcounter-1)*T_step
     if ((T.gt.T_max).and.(T.lt.(T_max+1.d0)))  exit
     if (T.gt.(T_max+1.d0)) T=T_max
     if (myid.eq.0) then
        write(aux2,"(I0)") NINT(T)
        path="T"//trim(adjustl(aux2))//"K"
        call change_directory(trim(adjustl(path))//C_NULL_CHAR)
        open(1,file="BTE.gruneisen_total",status="replace")
        write(1,*) total_grun(energy,grun)
        close(1)
        write(101,"(F7.1,E14.5)") T,total_grun(energy,grun)
        flush(101)
        call change_directory(".."//C_NULL_CHAR)
     end if
  end do
  close(101)
  deallocate(grun)

  ! This is the most expensive part of the calculation: obtaining the
  ! three-phonon / four-phonon scattering 
  ! amplitudes for all allowed processes.
  ! When the iterative solution to the full linearized BTE is not
  ! requested (i.e., when the relaxation-time approximation is
  ! enough) we use optimized routines with a much smaller memory footprint.
  ! weighted phase space volume per mode .
  allocate(Pspace_plus_total(Nbands,Nlist))
  allocate(Pspace_minus_total(Nbands,Nlist))
  if (four_phonon) then
    allocate(Pspace_plusplus_total(Nbands,Nlist))
    allocate(Pspace_plusminus_total(Nbands,Nlist))
    allocate(Pspace_minusminus_total(Nbands,Nlist))
  end if


  open(303,file="BTE.KappaTensorVsT_RTA")
  if(convergence) then
     open(403,file="BTE.KappaTensorVsT_CONV")
     allocate(Indof2ndPhonon_plus(Ntotal_plus))
     allocate(Indof3rdPhonon_plus(Ntotal_plus))
     allocate(Indof2ndPhonon_minus(Ntotal_minus))
     allocate(Indof3rdPhonon_minus(Ntotal_minus))
     allocate(Gamma_plus(Ntotal_plus))
     allocate(Gamma_minus(Ntotal_minus))
     if (four_phonon_iteration) then
       allocate(Indof2ndPhonon_plusplus(Ntotal_plusplus))
       allocate(Indof3rdPhonon_plusplus(Ntotal_plusplus))
       allocate(Indof4thPhonon_plusplus(Ntotal_plusplus))
       allocate(Indof2ndPhonon_plusminus(Ntotal_plusminus))
       allocate(Indof3rdPhonon_plusminus(Ntotal_plusminus))
       allocate(Indof4thPhonon_plusminus(Ntotal_plusminus))
       allocate(Indof2ndPhonon_minusminus(Ntotal_minusminus))
       allocate(Indof3rdPhonon_minusminus(Ntotal_minusminus))
       allocate(Indof4thPhonon_minusminus(Ntotal_minusminus))
       allocate(Gamma_plusplus(Ntotal_plusplus))
       allocate(Gamma_plusminus(Ntotal_plusminus))
       allocate(Gamma_minusminus(Ntotal_minusminus))
     end if
  end if
  if (myid.eq.0) print*, "Info: start calculating kappa"
  do Tcounter=1,CEILING((T_max-T_min)/T_step)+1
     T=T_min+(Tcounter-1)*T_step
     if ((T.gt.T_max).and.(T.lt.(T_max+1.d0)))  exit
     if (T.gt.(T_max+1.d0)) T=T_max
     if (myid.eq.0) then
        print *, "Info: Temperature=", T
        write(aux2,"(I0)") NINT(T)
        path="T"//trim(adjustl(aux2))//"K"
        call change_directory(trim(adjustl(path))//C_NULL_CHAR)
     end if
     if(convergence) then
        call Ind_driver(energy,velocity,eigenvect,Nlist,List,IJK,N_plus,N_minus,&
             Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,&
             Indof2ndPhonon_plus,Indof3rdPhonon_plus,Gamma_plus,&
             Indof2ndPhonon_minus,Indof3rdPhonon_minus,Gamma_minus,rate_scatt,&
             rate_scatt_plus,rate_scatt_minus,Pspace_plus_total,Pspace_minus_total)
        ! Will add fou-phonon iteration capabilities when we publish
        ! tag four_phonon_iteration is reserved here
        if (four_phonon.and.four_phonon_iteration.eq. .false.) then
          call RTA_driver_4ph(energy,velocity,eigenvect,Nlist,List,IJK,&
                Ntri_4fc,Psi,R_s,R_t,R_u,Index_r,Index_s,Index_t,Index_u,rate_scatt_4ph,&
                rate_scatt_plusplus,rate_scatt_plusminus,rate_scatt_minusminus,&
                Pspace_plusplus_total,Pspace_plusminus_total,Pspace_minusminus_total,&
                rate_scatt_plusplus_N,rate_scatt_plusminus_N,rate_scatt_minusminus_N,&
                rate_scatt_plusplus_U,rate_scatt_plusminus_U,rate_scatt_minusminus_U)
        end if
     else
        call RTA_driver(energy,velocity,eigenvect,Nlist,List,IJK,&
             Ntri,Phi,R_j,R_k,Index_i,Index_j,Index_k,rate_scatt,rate_scatt_plus,&
             rate_scatt_minus,Pspace_plus_total,Pspace_minus_total,&
             rate_scatt_plus_N,rate_scatt_minus_N,rate_scatt_plus_U,rate_scatt_minus_U)
        if (four_phonon) then
          call RTA_driver_4ph(energy,velocity,eigenvect,Nlist,List,IJK,&
              Ntri_4fc,Psi,R_s,R_t,R_u,Index_r,Index_s,Index_t,Index_u,rate_scatt_4ph,&
              rate_scatt_plusplus,rate_scatt_plusminus,rate_scatt_minusminus,&
              Pspace_plusplus_total,Pspace_plusminus_total,Pspace_minusminus_total,&
              rate_scatt_plusplus_N,rate_scatt_plusminus_N,rate_scatt_minusminus_N,&
              rate_scatt_plusplus_U,rate_scatt_plusminus_U,rate_scatt_minusminus_U)
        end if
     end if

     
     if(convergence) then
        if(myid.eq.0) then
            open(1,file="BTE.WP3_plus",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(2E14.5)") energy(list(ll),i),Pspace_plus_total(i,ll)
              end do
            end do
            close(1)
            open(1,file="BTE.WP3_minus",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(2E14.5)") energy(list(ll),i),Pspace_minus_total(i,ll)
              end do
            end do
            close(1)
            open(1,file="BTE.WP3",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(2E14.5)") energy(list(ll),i),Pspace_plus_total(i,ll)+Pspace_minus_total(i,ll)
              end do
            end do
            close(1)
            if (four_phonon) then
              if (four_phonon_iteration.eq. .false.) then
                write(*,*) "Info: four_phonon_iteration=.false., using RTA on 4ph calculations"
              end if
              open(1,file="BTE.WP4_plusplus",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(2E14.5)") energy(list(ll),i),Pspace_plusplus_total(i,ll)
                end do
              end do
              close(1)
              open(1,file="BTE.WP4_plusminus",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(2E14.5)") energy(list(ll),i),Pspace_plusminus_total(i,ll)
                end do
              end do
              close(1)
              open(1,file="BTE.WP4_minusminus",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(2E14.5)") energy(list(ll),i),Pspace_minusminus_total(i,ll)
                end do
              end do
              close(1)
              open(1,file="BTE.WP4",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(2E14.5)") energy(list(ll),i),Pspace_plusplus_total(i,ll)+&
                                        Pspace_plusminus_total(i,ll)+Pspace_minusminus_total(i,ll)
                end do
              end do
              close(1)
            end if
        end if

        write(aux,"(I0)") Nbands
        if(myid.eq.0) then
            open(1,file="BTE.w_3ph",status="replace")
            open(2,file="BTE.w_3ph_plus",status="replace")
            open(3,file="BTE.w_3ph_minus",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(2E20.10)") energy(list(ll),i),rate_scatt(i,ll)
                  write(2,"(2E20.10)") energy(list(ll),i),rate_scatt_plus(i,ll)
                  write(3,"(2E20.10)") energy(list(ll),i),rate_scatt_minus(i,ll)
              end do
            end do
            close(1)
            close(2)
            close(3)
            if (four_phonon) then
              open(1,file="BTE.w_4ph",status="replace")
              open(2,file="BTE.w_4ph_plusplus",status="replace")
              open(3,file="BTE.w_4ph_plusminus",status="replace")
              open(4,file="BTE.w_4ph_minusminus",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(2E20.10)") energy(list(ll),i),rate_scatt_4ph(i,ll)
                    write(2,"(2E20.10)") energy(list(ll),i),rate_scatt_plusplus(i,ll)
                    write(3,"(2E20.10)") energy(list(ll),i),rate_scatt_plusminus(i,ll)
                    write(4,"(2E20.10)") energy(list(ll),i),rate_scatt_minusminus(i,ll)
                end do
              end do
              close(1)
              close(2)
              close(3)
              close(4)
              if (four_phonon_iteration.eq. .false.) then
                write(*,*) "Info: Output four-phonon N/U scattering rates"
                open(5,file="BTE.w_4ph_normal",status="replace")
                open(6,file="BTE.w_4ph_Umklapp",status="replace")
                do i=1,Nbands
                  do ll=1,Nlist
                      write(5,"(5E20.10)") energy(list(ll),i),rate_scatt_plusplus_N(i,ll),rate_scatt_plusminus_N(i,ll),&
                                            rate_scatt_minusminus_N(i,ll),rate_scatt_plusplus_N(i,ll)+&
                                            rate_scatt_plusminus_N(i,ll)+rate_scatt_minusminus_N(i,ll)
                      write(6,"(5E20.10)") energy(list(ll),i),rate_scatt_plusplus_U(i,ll),rate_scatt_plusminus_U(i,ll),&
                                            rate_scatt_minusminus_U(i,ll),rate_scatt_plusplus_U(i,ll)+&
                                            rate_scatt_plusminus_U(i,ll)+rate_scatt_minusminus_U(i,ll)
                  end do
                end do
                close(5)
                close(6)
              end if
            end if
        end if
     else 
        if(myid.eq.0) then
          write(*,*) "Info: convergence=.false."
          open(1,file="BTE.WP3_plus",status="replace")
          do i=1,Nbands
            do ll=1,Nlist
                write(1,"(6E14.5)") energy(list(ll),i),Pspace_plus_total(i,ll)
            end do
          end do
          close(1)
          open(1,file="BTE.WP3_minus",status="replace")
          do i=1,Nbands
            do ll=1,Nlist
                write(1,"(6E14.5)") energy(list(ll),i),Pspace_minus_total(i,ll)
            end do
          end do
          close(1)
          open(1,file="BTE.WP3",status="replace")
          do i=1,Nbands
            do ll=1,Nlist
                write(1,"(2E14.5)") energy(list(ll),i),Pspace_plus_total(i,ll)+Pspace_minus_total(i,ll)
            end do
          end do
          close(1)
          if (four_phonon) then
            open(1,file="BTE.WP4_plusplus",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(2E14.5)") energy(list(ll),i),Pspace_plusplus_total(i,ll)
              end do
            end do
            close(1)
            open(1,file="BTE.WP4_plusminus",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(2E14.5)") energy(list(ll),i),Pspace_plusminus_total(i,ll)
              end do
            end do
            close(1)
            open(1,file="BTE.WP4_minusminus",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(2E14.5)") energy(list(ll),i),Pspace_minusminus_total(i,ll)
              end do
            end do
            close(1)
            open(1,file="BTE.WP4",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(2E14.5)") energy(list(ll),i),Pspace_plusplus_total(i,ll)+&
                                      Pspace_plusminus_total(i,ll)+Pspace_minusminus_total(i,ll)
              end do
            end do
            close(1)
          end if
        end if

        write(aux,"(I0)") Nbands
        if(myid.eq.0) then
            open(1,file="BTE.w_3ph",status="replace")
            open(2,file="BTE.w_3ph_plus",status="replace")
            open(3,file="BTE.w_3ph_minus",status="replace")
            open(4,file="BTE.w_3ph_normal",status="replace")
            open(5,file="BTE.w_3ph_Umklapp",status="replace")
            do i=1,Nbands
              do ll=1,Nlist
                  write(1,"(6E20.10)") energy(list(ll),i),rate_scatt(i,ll)
                  write(2,"(6E20.10)") energy(list(ll),i),rate_scatt_plus(i,ll)
                  write(3,"(6E20.10)") energy(list(ll),i),rate_scatt_minus(i,ll)
                  write(4,"(4E20.10)") energy(list(ll),i),rate_scatt_plus_N(i,ll),rate_scatt_minus_N(i,ll),rate_scatt_plus_N(i,ll)+rate_scatt_minus_N(i,ll)
                  write(5,"(4E20.10)") energy(list(ll),i),rate_scatt_plus_U(i,ll),rate_scatt_minus_U(i,ll),rate_scatt_plus_U(i,ll)+rate_scatt_minus_U(i,ll)
              end do
            end do
            close(1)
            close(2)
            close(3)
            close(4)
            close(5)
            if (four_phonon) then
              write(*,*) "Info: Applying RTA scheme on both 3ph and 4ph"
              write(*,*) "Info: Output four-phonon N/U scattering rates"
              open(1,file="BTE.w_4ph",status="replace")
              open(2,file="BTE.w_4ph_plusplus",status="replace")
              open(3,file="BTE.w_4ph_plusminus",status="replace")
              open(4,file="BTE.w_4ph_minusminus",status="replace")
              open(5,file="BTE.w_4ph_normal",status="replace")
              open(6,file="BTE.w_4ph_Umklapp",status="replace")
              do i=1,Nbands
                do ll=1,Nlist
                    write(1,"(6E20.10)") energy(list(ll),i),rate_scatt_4ph(i,ll)
                    write(2,"(6E20.10)") energy(list(ll),i),rate_scatt_plusplus(i,ll)
                    write(3,"(6E20.10)") energy(list(ll),i),rate_scatt_plusminus(i,ll)
                    write(4,"(6E20.10)") energy(list(ll),i),rate_scatt_minusminus(i,ll)
                    write(5,"(5E20.10)") energy(list(ll),i),rate_scatt_plusplus_N(i,ll),rate_scatt_plusminus_N(i,ll),&
                                          rate_scatt_minusminus_N(i,ll),rate_scatt_plusplus_N(i,ll)+&
                                          rate_scatt_plusminus_N(i,ll)+rate_scatt_minusminus_N(i,ll)
                    write(6,"(5E20.10)") energy(list(ll),i),rate_scatt_plusplus_U(i,ll),rate_scatt_plusminus_U(i,ll),&
                                          rate_scatt_minusminus_U(i,ll),rate_scatt_plusplus_U(i,ll)+&
                                          rate_scatt_plusminus_U(i,ll)+rate_scatt_minusminus_U(i,ll)
                end do
              end do
              close(1)
              close(2)
              close(3)
              close(4)
              close(5)
              close(6)
            end if
        end if
     end if
     
     ! Obtain the total scattering rates in the relaxation time approximation.
     ! Four-phonon scattering rates are denoted rate_scatt_4ph
     rate_scatt=rate_scatt+rate_scatt_isotope+rate_scatt_4ph

     if(myid.eq.0) then
        open(1,file="BTE.w",status="replace")
        do i=1,Nbands
           do ll = 1,Nlist
              write(1,"(2E20.10)") energy(list(ll),i),rate_scatt(i,ll)
           end do
        end do
        close(1)
     end if

     tau_zero=0.d0
     do ll = 1,Nlist
        do i=1,Nbands
           if(rate_scatt(i,ll).ne.0) then
              tau_zero(i,ll)=1./rate_scatt(i,ll)
           end if
        end do
     end do


     ! Set up everything to start the iterative process.
     call iteration0(Nlist,Nequi,ALLEquiList,energy,velocity,tau_zero,F_n)
     F_n_0=F_n
     if(myid.eq.0) then
        ! Open all output files.
        open(2001,file="BTE.kappa",status="replace")
        open(2002,file="BTE.kappa_tensor",status="replace")
        open(2003,file="BTE.kappa_scalar",status="replace")
        call TConduct(energy,velocity,F_n,ThConductivity,ThConductivityMode)
        do ll=1,nbands
           call symmetrize_tensor(ThConductivity(ll,:,:))
        end do
        write(aux,"(I0)") 9*Nbands
        write(2001,"(I9,"//trim(adjustl(aux))//"E20.10)") 0,ThConductivity
        write(2002,"(I9,9E20.10,E20.10)") 0,sum(ThConductivity,dim=1)
        write(2003,"(I9,E20.10,E20.10)") 0,&
             sum(sum(ThConductivity,dim=1),reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.
        write(303,"(F7.1,9E14.5)") T,sum(ThConductivity,dim=1)
        flush(303)

        ! Iterate to convergence if desired.
        if(convergence) then
           do ii=1,maxiter
              kappa_old=sum(ThConductivity,dim=1)
            ! will add four_phonon_iteration when we publish
            call iteration(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
                    Ntotal_plus,Ntotal_minus,Indof2ndPhonon_plus,Indof3rdPhonon_plus,&
                    Indof2ndPhonon_minus,Indof3rdPhonon_minus,energy,velocity,&
                    Gamma_plus,Gamma_minus,tau_zero,F_n)
              ! Correct F_n to prevent it drifting away from the symmetry of the system.
              do ll=1,nptk
                 F_n(:,ll,:)=transpose(matmul(symmetrizers(:,:,ll),transpose(F_n(:,ll,:))))
              end do
              call TConduct(energy,velocity,F_n,ThConductivity,ThConductivityMode)
              do ll=1,nbands
                 call symmetrize_tensor(ThConductivity(ll,:,:))
              end do
              write(2001,"(I9,"//trim(adjustl(aux))//"E20.10)") ii,ThConductivity
              flush(2001)
              write(2002,"(I9,9E20.10)") ii,sum(ThConductivity,dim=1)
              flush(2002)
              write(2003,"(I9,E20.10)") ii,&
                   sum(sum(ThConductivity,dim=1),reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.
              flush(2003)
              relchange=twonorm3x3(sum(ThConductivity,dim=1)-kappa_old)/&
                   twonorm3x3(kappa_old)
              write(*,*) "Info: Iteration",ii
              write(*,*) "Info:","Relative change","=",relchange
              if(relchange.lt.eps) exit
           end do
           write(403,"(F7.1,9E14.5,I6)") T,sum(ThConductivity,dim=1),ii
           flush(403)
        end if
        close(2001)
        close(2002)
        close(2003)

        ! Write out the converged scattering rates.
        do ll=1,Nlist
           do ii=1,Nbands
              tau(ii,ll)=dot_product(F_n(ii,List(ll),:),velocity(List(ll),ii,:))/&
                   (dot_product(velocity(List(ll),ii,:),velocity(List(ll),ii,:))*energy(List(ll),ii))
           end do
        end do
        write(aux,"(I0)") Nbands
        open(1,file="BTE.w_final",status="replace")
        do i=1,Nbands
           do ll = 1,Nlist
              write(1,"(2E20.10)") energy(list(ll),i),  1./tau(i,ll)
           end do
        end do
        close(1)
        do ll=1,nptk
           do ii=1,Nbands
              tau2(ii,ll)=dot_product(F_n(ii,ll,:),velocity(ll,ii,:))/&
                   (dot_product(velocity(ll,ii,:),velocity(ll,ii,:))*energy(ll,ii))
           end do
        end do
        write(aux,"(I0)") Nbands

        ! If results for nanowires have been requested, obtain a lower bound
        ! for the thermal conductivity along each crystallographic orientation
        ! by using the bulk RTA results.
        if(nanowires) then
           do iorient=1,norientations
              write(sorientation,"(I128)") iorient
              open(3001,file="BTE.kappa_nw_"//trim(adjustl(sorientation))//"_lower",status="replace")
              do ii=1,nptk
                 do jj=1,Nbands
                    v_or(ii,jj)=dot_product(velocity(ii,jj,:),uorientations(:,iorient))
                    F_or(jj,ii)=dot_product(F_n_0(jj,ii,:),uorientations(:,iorient))
                 end do
              end do
              do mm=1,Nwires
                 radnw=radnw_range(mm)
                 call ScalingOfTau(Nlist,Nequi,ALLEquiList,v_or,velocity,tau_zero,radnw,ffunc)
                 do ii=1,nptk
                    do jj=1,Nbands
                       F_n_aux(jj,ii)=F_or(jj,ii)*ffunc(ii,jj)
                    end do
                 end do
                 call TConductScalar(energy,v_or,F_n_aux,kappa_or)
                 write(3001,"(E30.20,"//trim(adjustl(aux))//"E20.10,E20.10)") 2.d0*radnw,&
                      kappa_or,sum(kappa_or)
              end do
              close(3001)
           end do
        end if
        allocate(ticks(nticks),cumulative_kappa(nbands,3,3,nticks))
        ! Cumulative thermal conductivity.
        call CumulativeTConduct(energy,velocity,F_n,ticks,cumulative_kappa)
        do ii=1,nticks
           do ll=1,nbands
              call symmetrize_tensor(cumulative_kappa(ll,:,:,ii))
           end do
        end do
        write(aux,"(I0)") 9*nbands+1
        open(2002,file="BTE.cumulative_kappa_tensor",status="replace")
        open(2003,file="BTE.cumulative_kappa_scalar",status="replace")
        do ii=1,nticks
           write(2002,"(10E20.10)") ticks(ii),&
                sum(cumulative_kappa(:,:,:,ii),dim=1)
           write(2003,"(2E20.10)") ticks(ii),&
                sum(sum(cumulative_kappa(:,:,:,ii),dim=1),&
                reshape((/((i==j,i=1,3),j=1,3)/),(/3,3/)))/3.
        end do
        close(2002)
        close(2003)
        deallocate(ticks,cumulative_kappa)

        ! Cumulative thermal conductivity vs angular frequency.
        allocate(ticks(nticks),cumulative_kappa(nbands,3,3,nticks))
        call CumulativeTConductOmega(energy,velocity,F_n,ticks,cumulative_kappa)
        do ii=1,nticks
           do ll=1,nbands
              call symmetrize_tensor(cumulative_kappa(ll,:,:,ii))
           end do
        end do
        write(aux,"(I0)") 9*nbands+1
        open(2002,file="BTE.cumulative_kappaVsOmega_tensor",status="replace")
        do ii=1,nticks
           write(2002,"(10E20.10)") ticks(ii),&
                sum(cumulative_kappa(:,:,:,ii),dim=1)
        end do
        close(2002)
        deallocate(ticks,cumulative_kappa)
     end if

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     ! If results for nanowires have been requested, repeat the iterative process
     ! for each desired orientation, introducing the appropriate scaling of
     ! relaxation times.
     if(nanowires) then
        kappa_wires=0.d00
        kappa_wires_reduce=0.d00
        kk=ceiling(float(nwires)/numprocs)
        do iorient=1,norientations
           if(myid.eq.0) then
              write(*,"(A,I0,A,3(x,I0))") "Info: nanowires with orientation ",&
                   iorient,":",orientations(:,iorient)
           end if
           write(sorientation,"(I128)") iorient
           do ii=1,nptk
              do jj=1,Nbands
                 v_or(ii,jj)=dot_product(velocity(ii,jj,:),uorientations(:,iorient))
                 F_or(jj,ii)=dot_product(F_n_0(jj,ii,:),uorientations(:,iorient))
              end do
           end do
           do mm=myid+1,nwires,numprocs
              radnw=radnw_range(mm)
              call ScalingOfTau(Nlist,Nequi,ALLEquiList,v_or,velocity,tau_zero,radnw,ffunc)
              do ii=1,nptk
                 do jj=1,Nbands
                    F_n_aux(jj,ii)=F_or(jj,ii)*ffunc(ii,jj)
                 end do
              end do
              call TConductScalar(energy,v_or,F_n_aux,kappa_or)
              if(convergence) then
                 do ii=1,maxiter
                    kappa_or_old=sum(kappa_or)
                    call iteration_scalar(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,&
                         Ntotal_plus,Ntotal_minus,Indof2ndPhonon_plus,Indof3rdPhonon_plus,&
                         Indof2ndPhonon_minus,Indof3rdPhonon_minus,energy,v_or,&
                         Gamma_plus,Gamma_minus,tau_zero,F_n_aux)
                    do ll=1,nptk
                       do jj=1,nbands
                          F_n_aux(jj,ll)=F_n_aux(jj,ll)*ffunc(ll,jj)
                       end do
                    end do
                    call TConductScalar(energy,v_or,F_n_aux,kappa_or)
                    relchange=abs((sum(kappa_or)-kappa_or_old)/kappa_or_old)
                    if(relchange.lt.eps)exit
                 end do
              end if
              kappa_wires_reduce(:,mm)=kappa_or
           end do
           call MPI_ALLREDUCE(kappa_wires_reduce,kappa_wires,Nbands*Nwires,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
           if(myid.eq.0) then
              write(aux,"(I0)") 3*nbands
              open(3001,file="BTE.kappa_nw_"//trim(adjustl(sorientation)),status="replace")
              do ii=1,Nwires
                 radnw=radnw_range(ii)
                 write(3001,"(E30.20,"//trim(adjustl(aux))//"E20.10,E20.10)") 2.d0*radnw,&
                      kappa_wires(:,ii),sum(kappa_wires(:,ii))
              end do
              close(3001)
           end if
        end do
     end if
     if (myid.eq.0) call change_directory(".."//C_NULL_CHAR)
     
  end do ! Tcounter


  if(myid.eq.0)write(*,*) "Info: normal exit"

  call free_config()
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
end program ShengBTE
