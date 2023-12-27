! Routines used to calculate the phonon spectrum.
module phonon_routines
  use misc
  use data
  use config
  use input
  implicit none

contains

  ! Create the q-point grid and compute all relevant properties.
  subroutine eigenDM(omega,eigenvect,velocity)
    implicit none

    include "mpif.h"

    real(kind=8),intent(out) :: omega(nptk,nbands),velocity(nptk,nbands,3)
    complex(kind=8),intent(out) :: eigenvect(nptk,Nbands,Nbands)

    real(kind=8),allocatable :: omega_reduce(:,:),velocity_reduce(:,:,:)
    complex(kind=8),allocatable :: eigenvect_reduce(:,:,:)
    real(kind=8) :: kspace(nptk,3)
    integer(kind=4) :: indexK,ii,jj,kk
    character(len=1) :: aux

    do ii=1,Ngrid(1)        ! rlattvec(:,1) direction
       do jj=1,Ngrid(2)     ! rlattvec(:,2) direction
          do kk=1,Ngrid(3)  ! rlattvec(:,3) direction
             indexK=((kk-1)*Ngrid(2)+(jj-1))*Ngrid(1)+ii
             kspace(indexK,:)=rlattvec(:,1)*(ii-1.0)/ngrid(1)+&
                  rlattvec(:,2)*(jj-1.0)/ngrid(2)+&
                  rlattvec(:,3)*(kk-1.0)/ngrid(3)
          end do
       end do
    end do
    allocate(omega_reduce(nptk,nbands),velocity_reduce(nptk,nbands,3),&
         eigenvect_reduce(nptk,Nbands,Nbands))
    omega_reduce=0.
    velocity_reduce=0.
    eigenvect_reduce=0.
    kk=ceiling(float(nptk)/numprocs)
    ii=min(nptk,kk*myid)+1
    jj=min(nptk,kk*(myid+1))
    ! The routine to be called depends on the input format, selected through a
    ! flag in the CONTROL file.
    if(espresso) then
       call phonon_espresso(kspace(ii:jj,:),omega_reduce(ii:jj,:),&
            velocity_reduce(ii:jj,:,:),eigenvect_reduce(ii:jj,:,:))
    else if (tdep) then
       call phonon_tdep(kspace(ii:jj,:),omega_reduce(ii:jj,:),&
            velocity_reduce(ii:jj,:,:),eigenvect_reduce(ii:jj,:,:))
    else
       call phonon_phonopy(kspace(ii:jj,:),omega_reduce(ii:jj,:),&
            velocity_reduce(ii:jj,:,:),eigenvect_reduce(ii:jj,:,:))
    end if
    call MPI_ALLREDUCE(omega_reduce,omega,nptk*nbands,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,kk)
    call MPI_ALLREDUCE(velocity_reduce,velocity,nptk*nbands*3,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,kk)
    call MPI_ALLREDUCE(eigenvect_reduce,eigenvect,nptk*nbands*nbands,&
         MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,kk)
    deallocate(omega_reduce,velocity_reduce,eigenvect_reduce)
    ! Make sure that group velocities have the right symmetry at each q point.
    ! This solves the problem of undefined components for degenerate modes.
    do ii=1,nptk
       velocity(ii,:,:)=transpose(&
            matmul(symmetrizers(:,:,ii),transpose(velocity(ii,:,:))))
    end do
    ! Make sure that acoustic frequencies and group velocities at Gamma
    ! are exactly zero.
    if(myid.eq.0) then
       write(*,*) "Info: about to set the acoustic frequencies at Gamma to zero"
       write(*,*) "Info: original values:"
       do ii=1,3
          write(aux,"(I1)") ii
          write(*,*) "Info: omega(1,"//aux//") =",omega(1,ii),"rad/ps"
       end do
    end if
    omega(1,1:3)=0.d0
    velocity(1,1:3,:)=0.
  end subroutine eigenDM

  ! Compute phonon dispersions, Phonopy style.
  subroutine phonon_phonopy(kpoints,omegas,velocities,eigenvect)
    implicit none

    real(kind=8),intent(in) :: kpoints(:,:)
    real(kind=8),intent(out) :: omegas(:,:),velocities(:,:,:)
    complex(kind=8),intent(out),optional :: eigenvect(:,:,:)

    real(kind=8),parameter :: prefactor=1745.91429109 ! THz^2 * amu * nm^3

    integer(kind=4) :: nk
    real(kind=8),allocatable :: mm(:,:)
    complex(kind=8),allocatable :: dyn_total(:,:),dyn_nac(:,:)
    complex(kind=8),allocatable :: ddyn_total(:,:,:),ddyn_nac(:,:,:)
    real(kind=8),allocatable :: fc_short(:,:,:,:,:,:,:)
    real(kind=8),allocatable :: fc_diel(:,:,:,:,:,:,:)
    real(kind=8),allocatable :: fc_total(:,:,:,:,:,:,:)

    integer(kind=4) :: i,j,ip,ik,neq
    integer(kind=4) :: ix1,iy1,iz1,iatom1,ix2,iy2,iz2,iatom2
    real(kind=8) :: tmp1,tmp2,tmp3,dmin,Rnorm
    real(kind=8) :: rcell(3),r(3),rl(3),rr(3,27),qr(27)
    complex(kind=8) :: ztmp,star

    real(kind=8), allocatable :: shortest(:,:)
    real(kind=8), allocatable :: omega2(:),rwork(:)
    complex(kind=8), allocatable :: work(:)
    integer(kind=4) :: nwork=1

    real(kind=8) :: dnrm2

    nk=size(kpoints,1)

    allocate(mm(natoms,natoms))
    allocate(omega2(nbands))
    allocate(rwork(max(1,9*natoms-2)))

    allocate(fc_diel(natoms,3,scell(1),scell(2),scell(3),natoms,3))
    allocate(fc_total(natoms,3,scell(1),scell(2),scell(3),natoms,3))

    do i=1,natoms
       mm(i,i)=masses(types(i))
       do j=i+1,natoms
          mm(i,j)=sqrt(masses(types(i))*masses(types(j)))
          mm(j,i)=mm(i,j)
       end do
    end do

    ! Read FORCE_CONSTANTS_2ND and reduce the constants using mm.
    call read2fc(fc_short)

    allocate(dyn_total(nbands,nbands))
    allocate(dyn_nac(nbands,nbands))
    allocate(ddyn_total(nbands,nbands,3))
    allocate(ddyn_nac(nbands,nbands,3))
    allocate(work(nwork))
    allocate(shortest(3,nk))

    ! Use the 1st BZ image of each q point to improve the behavior of
    ! the non-analytic correction.
    do ik=1,nk
       shortest(:,ik)=kpoints(ik,:)
       tmp1=dnrm2(3,shortest(:,ik),1)
       do ix1=-2,2
          do iy1=-2,2
             do iz1=-2,2
                r=kpoints(ik,:)+ix1*rlattvec(:,1)+iy1*rlattvec(:,2)+&
                     iz1*rlattvec(:,3)
                tmp2=dnrm2(3,r,1)
                if(tmp2.lt.tmp1) then
                   tmp1=tmp2
                   shortest(:,ik)=r
                end if
             end do
          end do
       end do
    end do

    do ik=1,nk
       dyn_total=0.
       dyn_nac=0.
       ddyn_total=0.
       ddyn_nac=0.
       fc_diel=0.
       ! If the nonanalytic flag is set to TRUE, add the electrostatic
       ! correction. No correction is applied exactly at \Gamma in
       ! order not to rely on guesses about directions.
       if(nonanalytic.and..not.all(shortest(:,ik).eq.0.)) then
          tmp3=dot_product(shortest(:,ik),matmul(epsilon,shortest(:,ik)))
          do iatom1=1,natoms
             do iatom2=1,natoms
                do i=1,3
                   do j=1,3
                      tmp1=dot_product(shortest(:,ik),born(:,i,iatom1))
                      tmp2=dot_product(shortest(:,ik),born(:,j,iatom2))
                      dyn_nac(3*(iatom1-1)+i,3*(iatom2-1)+j)=tmp1*tmp2/&
                           mm(iatom1,iatom2)
                      ! The derivatives of the nonanalytic correction
                      ! will be needed later to make group velocities
                      ! and frequencies completely consistent.
                      do ip=1,3
                         ddyn_nac(3*(iatom1-1)+i,3*(iatom2-1)+j,ip)=&
                              tmp1*born(ip,j,iatom2)+tmp2*born(ip,i,iatom1)-&
                              2.*tmp1*tmp2*dot_product(epsilon(ip,:),shortest(:,ik))/tmp3
                      end do
                      ddyn_nac(3*(iatom1-1)+i,3*(iatom2-1)+j,:)=&
                           ddyn_nac(3*(iatom1-1)+i,3*(iatom2-1)+j,:)/&
                           mm(iatom1,iatom2)
                   end do
                end do
             end do
          end do
          dyn_nac=prefactor*dyn_nac/tmp3/V
          ddyn_nac=prefactor*ddyn_nac/tmp3/V
          ! Transform back to real space to obtain a correction to the
          ! short-range force constants.
          do iatom1=1,natoms
             do iatom2=1,natoms
                do i=1,3
                   do j=1,3
                      fc_diel(iatom1,i,:,:,:,iatom2,j)=real(dyn_nac(3*(iatom1-1)+i,&
                           3*(iatom2-1)+j))
                   end do
                end do
             end do
          end do
          fc_diel=fc_diel/(scell(1)*scell(2)*scell(3))
       end if

       ! Force constants with long-range correction.
       fc_total=fc_short+fc_diel
       ! Build the dynamical matrix and its derivatives.
       do iatom1=1,natoms
          do iatom2=1,natoms
             do ix1=1,scell(1)
                do iy1=1,scell(2)
                   do iz1=1,scell(3)
                      rcell=matmul(lattvec,(/ix1,iy1,iz1/)-(/1,1,1/))
                      r=cartesian(:,iatom1)-cartesian(:,iatom2)+rcell
                      dmin=huge(dmin)
                      do ix2=-2,2
                         do iy2=-2,2
                            do iz2=-2,2
                               rl=ix2*scell(1)*lattvec(:,1)+iy2*scell(2)*lattvec(:,2)+&
                                    iz2*scell(3)*lattvec(:,3)
                               Rnorm=dnrm2(3,rl+r,1)
                               if(abs(Rnorm-dmin).gt.1e-5) then
                                  if(Rnorm.lt.dmin) then
                                     neq=1
                                     dmin=Rnorm
                                     qr(neq)=dot_product(kpoints(ik,:),rl+rcell)
                                     rr(:,neq)=rl+rcell
                                  endif
                               else
                                  neq=neq+1
                                  qr(neq)=dot_product(kpoints(ik,:),rl+rcell)
                                  rr(:,neq)=rl+rcell
                               endif
                            end do
                         end do
                      end do
                      star=0.
                      do ip=1,neq
                         ztmp=phexp(-qr(ip))/neq
                         star=star+ztmp
                         do i=1,3
                            do j=1,3
                               dyn_total(3*(iatom1-1)+i,3*(iatom2-1)+j)=&
                                    dyn_total(3*(iatom1-1)+i,3*(iatom2-1)+j)+&
                                    ztmp*fc_total(iatom2,j,ix1,iy1,iz1,iatom1,i)
                               ddyn_total(3*(iatom1-1)+i,3*(iatom2-1)+j,:)=&
                                    ddyn_total(3*(iatom1-1)+i,3*(iatom2-1)+j,:)-&
                                    iunit*ztmp*rr(:,ip)*fc_total(iatom2,j,ix1,iy1,iz1,iatom1,i)
                            end do
                         end do
                      end do
                      if(nonanalytic.and..not.all(kpoints(ik,:).eq.0)) then
                         do i=1,3
                            do j=1,3
                               ddyn_total(3*(iatom1-1)+i,3*(iatom2-1)+j,:)=&
                                    ddyn_total(3*(iatom1-1)+i,3*(iatom2-1)+j,:)+&
                                    star*ddyn_nac(3*(iatom1-1)+i,3*(iatom2-1)+j,:)/&
                                    (scell(1)*scell(2)*scell(3))
                            end do
                         end do
                      end if
                   end do
                end do
             end do
          end do
       end do

       ! Frequencies squared result from a diagonalization of the
       ! dynamical matrix. The first call to zheev serves to ensure that
       ! enough space has been allocated for this.
       call zheev("V","U",nbands,dyn_total,nbands,omega2,work,-1,rwork,i)

       if(real(work(1)).gt.nwork) then
          nwork=nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V","U",nbands,dyn_total,nbands,omega2,work,nwork,rwork,i)

       ! Eigenvectors are also returned if required.
       if(present(eigenvect)) then
          eigenvect(ik,:,:)=transpose(dyn_total)
       end if

       ! As is conventional, imaginary frequencies are returned as negative.
       omegas(ik,:)=sign(sqrt(abs(omega2)),omega2)

       ! Group velocities are obtained perturbatively. This is very
       ! advatageous with respect to finite differences.
       do i=1,nbands
          do ip=1,3
             velocities(ik,i,ip)=real(dot_product(dyn_total(:,i),&
                  matmul(ddyn_total(:,:,ip),dyn_total(:,i))))
          end do
          velocities(ik,i,:)=velocities(ik,i,:)/(2.*omegas(ik,i))
       end do
    end do
    deallocate(mm,omega2,rwork,fc_short,fc_diel,fc_total,&
         dyn_total,dyn_nac,ddyn_total,ddyn_nac,work,shortest)
  end subroutine phonon_phonopy

  ! Adapted from the code of Quantum Espresso (
  ! http://www.quantum-espresso.org/ ), licensed under the GPL.
  subroutine phonon_espresso(kpoints,omegas,velocities,eigenvect)
    implicit none

    real(kind=8),intent(in) :: kpoints(:,:)
    real(kind=8),intent(out) :: omegas(:,:),velocities(:,:,:)
    complex(kind=8),optional,intent(out) :: eigenvect(:,:,:)

    ! QE's 2nd-order files are in Ryd units.
    real(kind=8),parameter :: bohr2nm=0.052917721092,toTHz=20670.687,&
         massfactor=1.8218779*6.022e-4

    integer(kind=4) :: ir,nreq,ntype,nat,ibrav,qscell(3)
    integer(kind=4) :: i,j,ipol,jpol,iat,jat,idim,jdim,t1,t2,t3,m1,m2,m3,ik
    integer(kind=4) :: ndim,nk,nwork,ncell_g(3)
    integer(kind=8),allocatable :: tipo(:)
    character(len=1) :: polar_key
    character(len=5),allocatable :: label(:)
    real(kind=8) :: weight,total_weight,exp_g,ck
    real(kind=8) :: celldm(6),r_ws(3),rws(124,0:3),wscell(3,0:3),at(3,3)
    real(kind=8) :: alpha,geg,gmax,kt,gr,volume_r,dnrm2
    real(kind=8) :: cell_r(1:3,0:3),cell_g(1:3,0:3)
    real(kind=8) :: zig(3),zjg(3),dgeg(3),t(0:3),g(0:3),g_old(0:3)
    real(kind=8), allocatable :: omega2(:),rwork(:)
    real(kind=8),allocatable :: k(:,:),mass(:),r(:,:),eps(:,:),mm(:,:),rr(:,:,:)
    real(kind=8),allocatable :: eival(:,:),vels(:,:,:),zeff(:,:,:),fc_s(:,:,:,:,:,:,:)
    complex(kind=8) :: auxi(3)
    complex(kind=8),allocatable :: cauxiliar(:),eigenvectors(:,:),work(:)
    complex(kind=8),allocatable :: dyn(:,:),dyn_s(:,:,:),dyn_g(:,:,:)
    complex(kind=8),allocatable :: ddyn(:,:,:),ddyn_s(:,:,:,:),ddyn_g(:,:,:,:)

    ! Quantum Espresso's 2nd-order format contains information about
    ! lattice vectors, atomic positions, Born effective charges and so
    ! forth in its header. The information is read but completely
    ! ignored. It is the user's responsibility to ensure that
    ! it is consistent with the CONTROL file.
    nwork=1
    nk=size(kpoints,1)
    open(1,file="espresso.ifc2",status="old")
    read(1,*) ntype,nat,ibrav,celldm(1:6)
    if (ibrav==0) then
      read(1,*) ((at(i,j),i=1,3),j=1,3)
    end if
    ntype=nelements
    nat=natoms
    ndim=3*nat

    allocate(omega2(nbands))
    allocate(work(nwork))
    allocate(rwork(max(1,9*natoms-2)))
    allocate(k(nk,3))
    allocate(label(ntype))
    allocate(mass(ntype))
    allocate(tipo(nat))
    allocate(r(nat,3))
    allocate(eps(3,3))
    allocate(zeff(nat,3,3))
    allocate(fc_s(3,3,nat,nat,scell(1),scell(2),scell(3)))
    allocate(mm(nat,nat))
    allocate(rr(nat,nat,3))
    allocate(dyn(ndim,ndim))
    allocate(dyn_s(nk,ndim,ndim))
    allocate(dyn_g(nk,ndim,ndim))
    allocate(ddyn(ndim,ndim,3))
    allocate(ddyn_s(nk,ndim,ndim,3))
    allocate(ddyn_g(nk,ndim,ndim,3))
    allocate(eival(ndim,nk))
    allocate(vels(ndim,nk,3))
    allocate(eigenvectors(ndim,ndim))
    allocate(cauxiliar(ndim))

    do i=1,ntype
       read(1,*) j,label(i),mass(i)
    end do
    mass=masses/massfactor
    do i=1,nat
       read(1,*) j,tipo(i),r(i,1:3)
    end do
    tipo=types
    r=transpose(matmul(lattvec,positions))/bohr2nm
    read(1,*) polar_key
    if(polar_key.eq."T") then
       do i=1,3
          read(1,*) eps(i,1:3)
       end do
       do i=1,nat
          read(1,*)
          do j=1,3
             read(1,*) zeff(i,j,1:3)
          end do
       end do
    end if
    eps=transpose(epsilon)
    do i=1,nat
       zeff(i,:,:)=transpose(born(:,:,i))
    end do
    read(1,*) qscell(1:3)
    ! Read the force constants.
    do i=1,3*3*nat*nat
       read(1,*) ipol,jpol,iat,jat
       do j=1,scell(1)*scell(2)*scell(3)
          read(1,*) t1,t2,t3,fc_s(ipol,jpol,iat,jat,t1,t2,t3)
       end do
    end do
    ! Enforce the conservation of momentum in the simplest way possible.
    ! Note that this is not necessary for the Phonopy format.
    do i=1,3
       do j=1,3
          do iat=1,nat
             fc_s(i,j,iat,iat,1,1,1)=fc_s(i,j,iat,iat,1,1,1)-&
                  sum(fc_s(i,j,iat,:,:,:,:))
          end do
       end do
    end do
    close(1)

    ! Make sure operations are performed in consistent units.
    k=kpoints*bohr2nm
    cell_r(:,1:3)=transpose(lattvec)/bohr2nm
    volume_r=V/bohr2nm**3
    do i=1,3
       cell_r(i,0)=dnrm2(3,cell_r(i,1:3),1)
    end do
    cell_g(:,1:3)=transpose(rlattvec)*bohr2nm
    do i=1,3
       cell_g(i,0)=dnrm2(3,cell_g(i,1:3),1)
    end do

    ! The dynamical matrix is built in a way similar to the previous
    ! subroutine.
    wscell(1,1:3)=cell_r(1,1:3)*scell(1)
    wscell(2,1:3)=cell_r(2,1:3)*scell(2)
    wscell(3,1:3)=cell_r(3,1:3)*scell(3)

    j=1
    do m1=-2,2
       do m2=-2,2
          do m3=-2,2
             if(all((/m1,m2,m3/).eq.0)) then
                cycle
             end if
             do i=1,3
                rws(j,i)=wscell(1,i)*m1+wscell(2,i)*m2+wscell(3,i)*m3
             end do
             rws(j,0)=0.5*dot_product(rws(j,1:3),rws(j,1:3))
             j=j+1
          end do
       end do
    end do

    do i=1,nat
       mm(i,i)=mass(tipo(i))
       rr(i,i,:)=0
       do j=i+1,nat
          mm(i,j)=sqrt(mass(tipo(i))*mass(tipo(j)))
          rr(i,j,1:3)=r(i,1:3)-r(j,1:3)
          mm(j,i)=mm(i,j)
          rr(j,i,1:3)=-rr(i,j,1:3)
       end do
    end do

    gmax=14.
    alpha=(2.*pi*bohr2nm/dnrm2(3,lattvec(:,1),1))**2
    geg=gmax*4.*alpha
    ncell_g=int(sqrt(geg)/cell_g(:,0))+1

    dyn_s=0.
    ddyn_s=0.

    do iat=1,nat
       do jat=1,nat
          total_weight=0.0d0
          do m1=-2*scell(1),2*scell(1)
             do m2=-2*scell(2),2*scell(2)
                do m3=-2*scell(3),2*scell(3)
                   do i=1,3
                      t(i)=m1*cell_r(1,i)+m2*cell_r(2,i)+m3*cell_r(3,i)
                      r_ws(i)=t(i)+rr(iat,jat,i)
                   end do
                   weight=0.d0
                   nreq=1
                   j=0
                   Do ir=1,124
                      ck=dot_product(r_ws,rws(ir,1:3))-rws(ir,0)
                      if(ck.gt.1e-6) then
                         j=1
                         cycle
                      end if
                      if(abs(ck).lt.1e-6) then
                         nreq=nreq+1
                      end if
                   end do
                   if(j.eq.0) then
                      weight=1.d0/dble(nreq)
                   end if
                   if(weight.gt.0.d0) then
                      t1=mod(m1+1,scell(1))
                      if(t1.le.0) then
                         t1=t1+scell(1)
                      end if
                      t2=mod(m2+1,scell(2))
                      if(t2.Le.0) then
                         t2=t2+scell(2)
                      end if
                      t3=mod(m3+1,scell(3))
                      if(t3.le.0) then
                         t3=t3+scell(3)
                      end if
                      do ik=1,nk
                         kt=dot_product(k(ik,1:3),t(1:3))
                         do ipol=1,3
                            idim = (iat-1)*3+ipol
                            do jpol=1,3
                               jdim = (jat-1)*3+jpol
                               dyn_s(ik,idim,jdim)=dyn_s(ik,idim,jdim)+&
                                    fc_s(ipol,jpol,iat,jat,t1,t2,t3)*&
                                    phexp(-kt)*weight
                               ddyn_s(ik,idim,jdim,1:3)=ddyn_s(ik,idim,jdim,1:3)-&
                                    iunit*t(1:3)*&
                                    fc_s(ipol,jpol,iat,jat,t1,t2,t3)*&
                                    phexp(-kt)*weight
                            end do
                         end do
                      end do
                   end if
                   total_weight=total_weight+weight
                end do
             end do
          end do
       end do
    end do
    ! The nonanalytic correction has two components in this
    ! approximation. Results may differ slightly between this method
    ! and the one implemented in the previous subroutine.
    dyn_g=0.
    ddyn_g=0.
    if(nonanalytic) then
       do m1=-ncell_g(1),ncell_g(1)
          do m2=-ncell_g(2),ncell_g(2)
             do m3=-ncell_g(3),ncell_g(3)
                g(1:3)=m1*cell_g(1,1:3)+&
                     m2*cell_g(2,1:3)+m3*cell_g(3,1:3)
                geg=dot_product(g(1:3),matmul(eps,g(1:3)))
                if(geg.gt.0.0d0.and.geg/alpha/4.0d0.lt.gmax) then
                   exp_g=exp(-geg/alpha/4.0d0)/geg
                   do iat=1,nat
                      zig(1:3)=matmul(g(1:3),zeff(iat,1:3,1:3))
                      auxi(1:3)=0.
                      do jat=1,nat
                         gr=dot_product(g(1:3),rr(iat,jat,1:3))
                         zjg(1:3)=matmul(g(1:3),zeff(jat,1:3,1:3))
                         auxi(1:3)=auxi(1:3)+zjg(1:3)*phexp(gr)
                      end do
                      do ipol=1,3
                         idim=(iat-1)*3+ipol
                         do jpol=1,3
                            jdim=(iat-1)*3+jpol
                            dyn_g(1:nk,idim,jdim)=dyn_g(1:nk,idim,jdim)-&
                                 exp_g*zig(ipol)*auxi(jpol)
                         end do
                      end do
                   end do
                end if
                g_old(0:3)=g(0:3)
                do ik=1,nk
                   g(1:3)=g_old(1:3)+k(ik,1:3)
                   geg=dot_product(g(1:3),matmul(eps,g(1:3)))
                   if (geg.gt.0.0d0.and.geg/alpha/4.0d0.lt.gmax) then
                      exp_g=exp(-geg/alpha/4.0d0)/geg
                      dgeg=matmul(eps+transpose(eps),g(1:3))
                      do iat=1,nat
                         zig(1:3)=matmul(g(1:3),zeff(iat,1:3,1:3))
                         do jat=1,nat
                            gr=dot_product(g(1:3),rr(iat,jat,1:3))
                            zjg(1:3)=matmul(g(1:3),zeff(jat,1:3,1:3))
                            do ipol=1,3
                               idim=(iat-1)*3+ipol
                               do jpol=1,3
                                  jdim=(jat-1)*3+jpol
                                  dyn_g(ik,idim,jdim)=dyn_g(ik,idim,jdim)+&
                                       exp_g*zig(ipol)*zjg(jpol)*phexp(gr)
                                  do i=1,3
                                     ddyn_g(ik,idim,jdim,i)=ddyn_g(ik,idim,jdim,i)+&
                                          exp_g*phexp(gr)*&
                                          (zjg(jpol)*zeff(iat,i,ipol)+zig(ipol)*zeff(jat,i,jpol)+&
                                          zig(ipol)*zjg(jpol)*iunit*rr(iat,jat,i)-&
                                          zig(ipol)*zjg(jpol)*(dgeg(i)/alpha/4.0+dgeg(i)/geg))
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end if
                end do
             end do
          end do
       end do
       dyn_g=dyn_g*8.*pi/volume_r
       ddyn_g=ddyn_g*8.*pi/volume_r
    end if
    ! Once the dynamical matrix has been built, the frequencies and
    ! group velocities are extracted exactly like in the previous
    ! subroutine.
    do ik=1,nk
       dyn(:,:)=dyn_s(ik,:,:)+dyn_g(ik,:,:)
       ddyn(:,:,:)=ddyn_s(ik,:,:,:)+ddyn_g(ik,:,:,:)
       do ipol=1,3
          do jpol=1,3
             do iat=1,nat
                do jat=1,nat
                   idim=(iat-1)*3+ipol
                   jdim=(jat-1)*3+jpol
                   dyn(idim,jdim)=dyn(idim,jdim)/mm(iat,jat)
                   ddyn(idim,jdim,1:3)=ddyn(idim,jdim,1:3)/mm(iat,jat)
                end do
             end do
          end do
       end do

       call zheev("V","U",nbands,dyn(:,:),nbands,omega2,work,-1,rwork,i)
       if(real(work(1)).gt.nwork) then
          nwork=nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V","U",nbands,dyn(:,:),nbands,omega2,work,nwork,rwork,i)

       if(present(eigenvect)) then
          eigenvect(ik,:,:)=transpose(dyn(:,:))
       end if

       omegas(ik,:)=sign(sqrt(abs(omega2)),omega2)

       do i=1,nbands
          do j=1,3
             velocities(ik,i,j)=real(dot_product(dyn(:,i),&
                  matmul(ddyn(:,:,j),dyn(:,i))))
          end do
          velocities(ik,i,:)=velocities(ik,i,:)/(2.*omegas(ik,i))
       end do
    end do
    ! Return the result to the units used in the rest of ShengBTE.
    omegas=omegas*toTHz
    velocities=velocities*toTHz*bohr2nm
    deallocate(k)
    deallocate(label)
    deallocate(mass)
    deallocate(tipo)
    deallocate(r)
    deallocate(eps)
    deallocate(zeff)
    deallocate(fc_s)
    deallocate(mm)
    deallocate(rr)
    deallocate(dyn)
    deallocate(dyn_s)
    deallocate(dyn_g)
    deallocate(ddyn)
    deallocate(ddyn_s)
    deallocate(ddyn_g)
    deallocate(eival)
    deallocate(vels)
    deallocate(eigenvectors)
    deallocate(cauxiliar)
    deallocate(work)
    deallocate(rwork)
    deallocate(omega2)
  end subroutine phonon_espresso

  ! Compute phonon properties, TDEP (Temperature Dependent Effective Potential) style
  ! Parts of the codes are adapted from TDEP source codes https://ollehellman.github.io/ (MIT license)
  ! First version: 2021/06/11; Last edited: 2023/12/04
  subroutine phonon_tdep(kpoints,omegas,velocities,eigenvect)
    implicit none
 
    real(kind=8),intent(in) :: kpoints(:,:)
    real(kind=8),intent(out) :: omegas(:,:),velocities(:,:,:)
    complex(kind=8),optional,intent(out) :: eigenvect(:,:,:)
 
    integer(kind=4) :: ntype,nat
    integer(kind=4) :: i,j,ipol,jpol,iat,jat,idim,jdim,m1,m2,m3,ik
    integer(kind=4) :: ndim,nk,nwork,ncell_g(3)
    integer(kind=8),allocatable :: tipo(:)
    real(kind=8) :: alpha,geg,gmax,volume_r,dnrm2
    real(kind=8) :: cell_r(1:3,0:3),cell_g(1:3,0:3)
    real(kind=8) :: rdotq
    real(kind=8), allocatable :: omega2(:),rwork(:)
    real(kind=8),allocatable :: k(:,:),mass(:),r(:,:),eps(:,:),mm(:,:),rr(:,:,:),zeff(:,:,:)
    complex(kind=8) :: cfac
    complex(kind=8),allocatable :: eigenvectors(:,:),work(:)
    complex(kind=8),allocatable :: dyn(:,:),dyn_s(:,:,:),dyn_g(:,:,:)
    complex(kind=8),allocatable :: ddyn(:,:,:),ddyn_s(:,:,:,:),ddyn_g(:,:,:,:)
    complex(kind=8),allocatable :: dyn_loc(:,:,:,:,:),ddyn_loc(:,:,:,:,:,:),&
                                   dyn_2_l_tmp(:,:,:,:),ddyn_2_l_tmp(:,:,:,:,:)
    real(kind=8),allocatable :: ifc_tdep(:,:,:,:),lv2(:,:,:)
    integer,allocatable   :: a2(:,:),neighbor(:)
    integer ::atom1,atom2,ii,jj,kk,maxneighbor
 
    real(kind=8),parameter :: bohr2nm=0.052917721092,toTHz=20670.687,&
    massfactor=1.8218779*6.022e-4,lo_bohr_to_A= 0.52917721067d0,&
    lo_Hartree_to_eV= 27.21138602d0
    real(kind=8) :: lo_forceconstant_2nd_eVA_to_HartreeBohr,&
                     lo_eV_to_Hartree,lo_A_to_bohr
    lo_A_to_bohr =1/lo_bohr_to_A
    lo_eV_to_Hartree=1/lo_Hartree_to_eV
    lo_forceconstant_2nd_eVA_to_HartreeBohr=lo_eV_to_Hartree/(lo_A_to_bohr**2)
 
 
    !> Basic things, similar to QE format
    nwork=1
    nk=size(kpoints,1)
    ntype=nelements
    nat=natoms
    ndim=3*nat
    allocate(omega2(nbands))
    allocate(work(nwork))
    allocate(rwork(max(1,9*natoms-2)))
    allocate(k(nk,3))
    allocate(mass(ntype))
    allocate(tipo(nat))
    allocate(r(nat,3))
    allocate(eps(3,3))
    allocate(zeff(nat,3,3))
    allocate(mm(nat,nat))
    allocate(rr(nat,nat,3))
    allocate(dyn(ndim,ndim))
    allocate(dyn_s(nk,ndim,ndim))
    allocate(dyn_g(nk,ndim,ndim))
    allocate(ddyn(ndim,ndim,3))
    allocate(ddyn_s(nk,ndim,ndim,3))
    allocate(ddyn_g(nk,ndim,ndim,3))
    allocate(eigenvectors(ndim,ndim))
    allocate(dyn_loc(3,3,nat,nat,nk),ddyn_loc(3,3,nat,nat,nk,3))
    allocate(dyn_2_l_tmp(3,3,nat,nat),ddyn_2_l_tmp(3,3,nat,nat,3))
    allocate(neighbor(nat))
    mass=masses/massfactor
    tipo=types
    r=transpose(matmul(lattvec,positions))/bohr2nm
    eps=transpose(epsilon)
    do i=1,nat
      zeff(i,:,:)=transpose(born(:,:,i))
    end do
    k=kpoints*bohr2nm
    cell_r(:,1:3)=transpose(lattvec)/bohr2nm
    volume_r=V/bohr2nm**3
    do i=1,3
       cell_r(i,0)=dnrm2(3,cell_r(i,1:3),1)
    end do
    cell_g(:,1:3)=transpose(rlattvec)*bohr2nm
    do i=1,3
       cell_g(i,0)=dnrm2(3,cell_g(i,1:3),1)
    end do
    do i=1,nat
      mm(i,i)=mass(tipo(i))
      rr(i,i,:)=0 
      do j=i+1,nat
         mm(i,j)=sqrt(mass(tipo(i))*mass(tipo(j)))
         rr(i,j,1:3)=r(i,1:3)-r(j,1:3)
         mm(j,i)=mm(i,j)
         rr(j,i,1:3)=-rr(i,j,1:3)
      end do
    end do
    ! These parameters are for Ewald summation: not exactly how TDEP does it but gives satisfactory results
    gmax=14.
    alpha=(2.*pi*bohr2nm/dnrm2(3,lattvec(:,1),1))**2
    geg=gmax*4.*alpha
    ncell_g=int(sqrt(geg)/cell_g(:,0))+1
 
    !> Read in TDEP force constants and construct dyn
    dyn_loc = 0.0d0
    maxneighbor = 200
    allocate(ifc_tdep(3,3,maxneighbor,nat),a2(maxneighbor,nat),lv2(maxneighbor,3,nat))
    open(123,File='infile.forceconstant',status='old')
    read(123,*)
    read(123,*)
    do atom1 = 1, nat
       read(123,*)neighbor(atom1)
       do j =1, neighbor(atom1)
          read(123,*)a2(j,atom1)
          atom2=a2(j,atom1)
          read(123,*)lv2(j,1:3,atom1)
          lv2(j,1:3,atom1) = lv2(j,1,atom1)*cell_r(1,1:3)+lv2(j,2,atom1)*cell_r(2,1:3)+lv2(j,3,atom1)*cell_r(3,1:3)
          do ii=1, 3
             read(123,*)ifc_tdep(ii,1:3,j,atom1)
             ifc_tdep(ii,1:3,j,atom1)=ifc_tdep(ii,1:3,j,atom1)*lo_forceconstant_2nd_eVA_to_HartreeBohr*2
          end do
         end do
    end do
    close(123)
    ! Fourier Transformation
    do ik=1, nk
       do atom1=1, nat
          do j=1, neighbor(atom1)
             atom2=a2(j,atom1)
             rdotq = (k(ik,1)*lv2(j,1,atom1) + k(ik,2)*lv2(j,2,atom1) +k(ik,3)*lv2(j,3,atom1) )
             cfac = cmplx(cos(rdotq),sin(rdotq),kind=8)
             dyn_loc(1:3,1:3,atom1,atom2,ik) = dyn_loc(1:3,1:3,atom1,atom2,ik) + cfac*ifc_tdep(1:3,1:3,j,atom1)
             do ii = 1,3
             do jj = 1,3
                ddyn_loc(ii,jj,atom1,atom2,ik,1:3) = ddyn_loc(ii,jj,atom1,atom2,ik,1:3)+&
                cfac*ifc_tdep(ii,jj,j,atom1)*iunit*lv2(j,1:3,atom1)
             end do
             end do
          end do
       end do
    end do
    deallocate(ifc_tdep,a2,lv2)
    ddyn_s=0.
    dyn_s =0.
    !Transform the dimension of the matrix in TDEP format to ShengBTE
    do ik  = 1,nk
       do iat =1, nat
          do jat =1, nat
             do ipol=1,3
                do jpol=1,3
                   idim = (iat-1)*3+ipol
                   jdim = (jat-1)*3+jpol
                   dyn_s(ik,idim,jdim) = dyn_loc(ipol,jpol,iat,jat,ik)
                   ddyn_s(ik,idim,jdim,1:3)= ddyn_loc(ipol,jpol,iat,jat,ik,1:3)
                end do
             end do
          end do
       end do
    end do
 
    !> NAC part
    ! Make sure you input the force constants with long-range subtracted by TDEP
    dyn_g = 0.0d0
    ddyn_g = 0.0d0
    if ( nonanalytic ) then
      do ik = 1, nk
         dyn_2_l_tmp = (0.0d0,0.0d0)
         ddyn_2_l_tmp = (0.0d0,0.0d0)
         call rgd_blk(dyn_2_l_tmp,ddyn_2_l_tmp,k(ik,:),nat,rr,eps,zeff,Ewald,ncell_g,cell_g,volume_r,.TRUE.)
         ! Transfer back to ShengBTE 2 indices
         do iat =1, nat
            do jat =1, nat
               do ipol=1,3
                  do jpol=1,3
                     idim = (iat-1)*3+ipol
                     jdim = (jat-1)*3+jpol
                     dyn_g(ik,idim,jdim) = dyn_2_l_tmp(ipol,jpol,jat,iat)
                     ddyn_g(ik,idim,jdim,1:3)= ddyn_2_l_tmp(ipol,jpol,jat,iat,1:3)
                  end do
               end do
            end do
         end do
      end do      
    end if
 
    !> Solve eigen problems
    ! Once the dynamical matrix has been built, the frequencies and
    ! group velocities are extracted exactly like in the previous
    ! subroutine.
    do ik=1,nk
      dyn(:,:)=dyn_s(ik,:,:)+dyn_g(ik,:,:)
      !call hermi(dyn(:,:), nat)
      ddyn(:,:,:)=ddyn_s(ik,:,:,:)+ddyn_g(ik,:,:,:)
      do ipol=1,3
         do jpol=1,3
            do iat=1,nat
               do jat=1,nat
                  idim=(iat-1)*3+ipol
                  jdim=(jat-1)*3+jpol
                  dyn(idim,jdim)=dyn(idim,jdim)/mm(iat,jat)
                  ddyn(idim,jdim,1:3)=ddyn(idim,jdim,1:3)/mm(iat,jat)
               end do
            end do
         end do
      end do
 
      call zheev("V","U",nbands,dyn(:,:),nbands,omega2,work,-1,rwork,i)
      if(real(work(1)).gt.nwork) then
         nwork=nint(2*real(work(1)))
         deallocate(work)
         allocate(work(nwork))
      end if
      call zheev("V","U",nbands,dyn(:,:),nbands,omega2,work,nwork,rwork,i)
 
      if(present(eigenvect)) then
         eigenvect(ik,:,:)=transpose(dyn(:,:))
      end if
 
      omegas(ik,:)=sign(sqrt(abs(omega2)),omega2)
 
      do i=1,nbands
         do j=1,3
            velocities(ik,i,j)=real(dot_product(dyn(:,i),&
                 matmul(ddyn(:,:,j),dyn(:,i))))
         end do
         velocities(ik,i,:)=velocities(ik,i,:)/(2.*omegas(ik,i))
      end do
   end do
   ! Return the result to the units used in the rest of ShengBTE.
   omegas=omegas*toTHz
   velocities=velocities*toTHz*bohr2nm
   deallocate(k)
   deallocate(mass)
   deallocate(tipo)
   deallocate(r)
   deallocate(eps)
   deallocate(zeff)
   deallocate(mm)
   deallocate(rr)
   deallocate(dyn)
   deallocate(dyn_s)
   deallocate(dyn_g)
   deallocate(ddyn)
   deallocate(ddyn_s)
   deallocate(ddyn_g)
   deallocate(eigenvectors)
   deallocate(work)
   deallocate(rwork)
   deallocate(omega2)

 end subroutine phonon_tdep

 !> subroutine for polar materials, called by phonon_tdep
 recursive subroutine rgd_blk(dyn,ddyn,q,nat,rr,epsil,zeu,lambda,ncell_g,bg,volume_r,chgmult)
   implicit none 
 
   integer     , intent(in)               :: nat
   real(kind=8), dimension(3), intent(in) :: q
   integer(kind=4), intent(in)            :: ncell_g(3)
   real(KIND=8), intent(in)               :: epsil(3,3),zeu(nat,3,3),rr(nat,nat,3),lambda,bg(1:3,0:3),volume_r
   complex(kind=8), intent(out)           :: dyn(3,3,nat,nat),ddyn(3,3,nat,nat,3)
   logical, intent(in)                    :: chgmult
 
   logical :: mult
   complex(kind=8) :: expikr,Chi
   real(kind=8), dimension(3,3) :: kk
   real(kind=8), dimension(3) :: Gvec,Kvec,tauvec,Keps
   real(kind=8) :: inv4lambda2,knorm,expLambdaKnorm,ikr,invknorm,partialChi,f0,Ksum
   integer      :: m1,m2,m3,i,j,na,nb,ii,jj
   complex(kind=8) ::D(3,3,nat,nat),DL1(3,3,nat,nat),D0(3,3,nat,nat),&
                           D2(3,3,nat,nat),D3(3,3,nat,nat),ddyn_temp(3,3,nat,nat,3),&
                           DQLx1(3,3,nat,nat),DQLy1(3,3,nat,nat),DQLz1(3,3,nat,nat),&
                           Dx0(3,3,nat,nat),Dy0(3,3,nat,nat),Dz0(3,3,nat,nat)
   Real(kind=8) :: dynamical_reciprocal(3,3,nat,nat),polar_onsite_correction(3,3,nat),&
                     m(3,3),re,im
   Real(kind=8) :: dynmat(3,3,nat,nat),qdir(3),f0_G,Kx,Ky,Kz
   real(kind=8),dimension(3,3) :: kkx,kky,kkz
   complex(kind=8),dimension(3)::v0
 
   if ( chgmult ) then
      mult=.true.
   else
      mult=.false.
   end if
 
   inv4lambda2=1.0d0/(4.0d0*(lambda**2))
   DQLx1=0.0d0
   DQLy1=0.0d0
   DQLz1=0.0d0
   DL1=0.0d0
   do m1=-ncell_g(1),ncell_g(1)
      do m2=-ncell_g(2),ncell_g(2)
         do m3=-ncell_g(3),ncell_g(3)
            Gvec=m1*bg(1,1:3)+m2*bg(2,1:3)+m3*bg(3,1:3)
            Kvec=(Gvec+q)!*2*pi!/celldm(1)
            if((Kvec(1)*Kvec(1)+Kvec(2)*Kvec(2)+Kvec(3)*Kvec(3))**0.5d0 .lt. 1e-10)cycle
            knorm=dot_product(Kvec,matmul(epsil,Kvec))
            invknorm=1.0d0/knorm
            expLambdaKnorm=exp(-knorm*inv4lambda2)
            ! Half of the Chi-function
            partialChi=expLambdaKnorm*invknorm
            ! Kronecker product K \otimes K
            kk(:,1)=Kvec(1:3)*Kvec(1)
            kk(:,2)=Kvec(1:3)*Kvec(2)
            kk(:,3)=Kvec(1:3)*Kvec(3)
            Kx=Kvec(1)
            Ky=Kvec(2)
            Kz=Kvec(3)
            kkx(:,1) = [ 2*Kx,      Ky,       Kz     ]
            kkx(:,2) = [   Ky,      0.0d0,   0.0d0   ]
            kkx(:,3) = [   Kz, 0.0d0,   0.0d0   ]
 
            kky(:,1) = [   0.0d0,   Kx,      0.0d0   ]
            kky(:,2) = [   Kx,      2*Ky,    Kz      ]
            kky(:,3) = [   0.0d0,   Kz,      0.0d0   ]
 
            kkz(:,1) = [   0.0d0,   0.0d0,   Kx      ]
            kkz(:,2) = [   0.0d0,   0.0d0,   Ky      ]
            kkz(:,3) = [   Kx,      Ky,      2*Kz    ]
 
            Keps=2*matmul(Kvec,epsil) 
            Keps=Keps*(invknorm+inv4lambda2)
 
            do na=1,nat
            do nb=na,nat
               tauvec=rr(nb,na,:) ! or rr(nb,na,:)
               ikr=dot_product(Kvec,tauvec) !*celldm(1)
               !expikr=cmplx( cos(ikr) , sin(ikr) ,kind=8)
               expikr=phexp(ikr)
               chi=partialChi*expikr
               DL1(:,:,na,nb)=DL1(:,:,na,nb)+kk*chi
               v0=iunit*tauvec-Keps
               DQLx1(:,:,na,nb)=DQLx1(:,:,na,nb)+kkx*chi+kk*chi*v0(1)
               DQLy1(:,:,na,nb)=DQLy1(:,:,na,nb)+kky*chi+kk*chi*v0(2)
               DQLz1(:,:,na,nb)=DQLz1(:,:,na,nb)+kkz*chi+kk*chi*v0(3)
            end do
            end do
         end do
      end do
   end do
   DL1=DL1*4.*pi/volume_r
   !ddyn
   DQLx1=DQLx1*4.*pi/volume_r
   DQLy1=DQLy1*4.*pi/volume_r
   DQLz1=DQLz1*4.*pi/volume_r
   ! Fix the other half
      D0 =0.0d0
      Dx0=0.0d0
      Dy0=0.0d0
      Dz0=0.0d0
      do na=1,nat
      do nb=na,nat
         if(na .eq. nb)then
            D0(:,:,na,nb)=DL1(:,:,na,nb)
            Dx0(:,:,na,nb)=DQLx1(:,:,na,nb)
            Dy0(:,:,na,nb)=DQLy1(:,:,na,nb)
            Dz0(:,:,na,nb)=DQLz1(:,:,na,nb)
         else
            D0(:,:,na,nb)=DL1(:,:,na,nb)
            Dx0(:,:,na,nb)=DQLx1(:,:,na,nb)
            Dy0(:,:,na,nb)=DQLy1(:,:,na,nb)
            Dz0(:,:,na,nb)=DQLz1(:,:,na,nb)
            do i=1,3
            do j=1,3
               D0(j,i,nb,na)=conjg(D0(i,j,na,nb))
               Dx0(j,i,nb,na)=conjg(Dx0(i,j,na,nb))
               Dy0(j,i,nb,na)=conjg(Dy0(i,j,na,nb))
               Dz0(j,i,nb,na)=conjg(Dz0(i,j,na,nb))
            enddo
            enddo
         endif
      enddo
      enddo
 
      if ( mult ) then
         ! Multiply in the charges
         !
         ddyn=0.0d0
         D=0.0d0
         do na=1,nat
            do nb=1,nat
               do j=1,3
                  do i=1,3
                     do jj=1,3
                     do ii=1,3
                        D(i,j,nb,na)=D(i,j,nb,na)+zeu(na,i,ii)*zeu(nb,j,jj)*D0(ii,jj,nb,na)
                        ddyn(i,j,nb,na,1)=ddyn(i,j,nb,na,1)+zeu(nb,i,ii)*zeu(na,j,jj)*Dx0(ii,jj,nb,na)
                        ddyn(i,j,nb,na,2)=ddyn(i,j,nb,na,2)+zeu(nb,i,ii)*zeu(na,j,jj)*Dy0(ii,jj,nb,na)
                        ddyn(i,j,nb,na,3)=ddyn(i,j,nb,na,3)+zeu(nb,i,ii)*zeu(na,j,jj)*Dz0(ii,jj,nb,na)
                     enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         ! Add the on-site correction
         D2=0.0d0
         call rgd_blk(D2,ddyn_temp,[0.0d0,0.0d0,0.0d0],nat,rr,epsil,zeu,lambda,ncell_g,bg,volume_r,.False.)
 
         dynamical_reciprocal=real(D2)
         D3=0.0d0
         do na=1,nat
         do nb=1,nat
            do j=1,3
            do i=1,3
               do jj=1,3
               do ii=1,3
                  D3(i,j,na,nb)=D3(i,j,na,nb)+zeu(na,i,ii)*zeu(nb,j,jj)*dynamical_reciprocal(ii,jj,na,nb)
               end do
               end do
            end do
            end do
         end do
         end do
         do na=1,nat
            m=0.0d0
            do nb=1,nat
               m=m+D3(:,:,na,nb)
            end do
            polar_onsite_correction(:,:,na)=-m
         enddo
 
         do na=1,nat
            D(:,:,na,na)=D(:,:,na,na)+polar_onsite_correction(:,:,na)
         end do
      else
         D=D0
         dyn=D0
      end if
 
 
      !!!chop the ugly number 
      if (chgmult) then
         do na=1,nat
         do nb=1,nat
            do i=1,3
            do j=1,3
               re=real(D(i,j,na,nb))
               im=aimag(D(i,j,na,nb))
               if(abs( re ).lt. 1e-10) re=0.0d0
               if(abs( im ).lt. 1e-10) im=0.0d0
               D(i,j,na,nb)=cmplx(re,im,kind=8)
            end do
            end do
         end do
         end do
      end if
 
      !!Are we gamma
      ! We abandon the correction at Gamma point. In Phonopy subroutines we have a comment:
      ! No correction is applied exactly at \Gamma in
      ! order not to rely on guesses about directions.
      dynmat=0.0d0
      ! IF(abs( q(1)) .lt. 1e-10 .and. abs (q(2)) .lt. 1e-10 .and. abs (q(3)) .lt. 1e-10) then
      !    qdir(1)=1.0d0
      !    qdir(2)=0.0d0
      !    qdir(3)=0.0d0
      !    f0_G= 0.0d0
      !    do j=1,3
      !    Do i=1,3
      !       f0_G=f0_G+qdir(i)*epsil(i,j)*qdir(j)
      !    Enddo
      !    Enddo
      !    f0_G=1.0d0/f0_G!*4.0d0*pi/omega
      !    Do na=1,nat
      !    Do nb=1,nat
      !       Do i=1,3
      !       Do j=1,3
      !          dynmat(i,j,na,nb)=dot_product(qdir,zeu(na,:,i)*dot_product(qdir,zeu(nb,:,j)))*f0_G
      !       Enddo
      !       Enddo
      !    Enddo
      !    Enddo
      ! ENDIF
      !!add together
      if(chgmult)then
         dyn=(dynmat+D)*2.0d0
         ddyn=ddyn*2.0d0
      end if
 
  end subroutine rgd_blk
 
 
  subroutine hermi(dyn_2, nat)
    implicit none
    INTEGER, INTENT(IN)            :: nat
    COMPLEX(KIND=8), INTENT(INOUT) :: dyn_2(3*nat,3*nat) 
    INTEGER                        :: dir1, dir2, atom1, atom2, mu, nu
    COMPLEX(KIND=8)                :: dyn_hermi(3,3,nat,nat)
    COMPLEX(KIND=8)                :: dyn_tmp(3,3,nat,nat)
    !
    !
    DO atom1 = 1, nat
       DO atom2 = 1, nat
          DO dir1 = 1, 3
             DO dir2 = 1, 3
                !
                mu = (atom1-1)*3+dir1
                nu = (atom2-1)*3+dir2
                dyn_tmp(dir1,dir2,atom1,atom2) = dyn_2(mu,nu)
                !
             END DO
          END DO
       END DO
    END DO
    !
    DO atom1 = 1, nat
       DO atom2 = 1, nat
          DO dir1 = 1, 3
             DO dir2 = 1, 3
                !
                dyn_hermi(dir1,dir2,atom1,atom2) = 0.50d0 * ( dyn_tmp(dir1,dir2,atom1,atom2) + CONJG(dyn_tmp(dir2,dir1,atom2,atom1)) )
                !
             END DO
          END DO
       END DO
    END DO
    !
    DO atom1 = 1, nat
       DO atom2 = 1, nat
          DO dir1 = 1, 3
             DO dir2 = 1, 3
                !
                mu = (atom1-1)*3+dir1
                nu = (atom2-1)*3+dir2
                dyn_2(mu,nu) = dyn_hermi(dir1,dir2,atom1,atom2)
                !
             END DO
          END DO
       END DO
    END DO
    !
    !
    RETURN

 end subroutine




end module phonon_routines
