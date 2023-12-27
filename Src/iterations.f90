! Code used to initialize F_n for the zeroth-order iteration of the
! BTE and to perform successive iterations.
module iterations
  use data
  use config
  use wedgetc
  implicit none

contains

  ! Fill F_n with its initial values for the iteration, computed from the
  ! RTA approximation to tau.
  subroutine iteration0(Nlist,Nequi,ALLEquiList,omega,velocity,tau_zero,F_n)
    implicit none

    integer(kind=4),intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    real(kind=8),intent(in) :: omega(nptk,Nbands),velocity(nptk,Nbands,3),tau_zero(Nbands,Nlist)
    real(kind=8),intent(out) :: F_n(Nbands,nptk,3)

    integer(kind=4) :: ii,kk,ll

    do ll=1,Nlist
       do ii=1,Nbands
          do kk=1,Nequi(ll)
             F_n(ii,ALLEquiList(kk,ll),:)=tau_zero(ii,ll)*velocity(ALLEquiList(kk,ll),ii,:)*&
                  omega(ALLEquiList(kk,ll),ii)
          end do
       end do
    end do
  end subroutine iteration0

  ! Advance the algorithm one iteration. F_n is used both as the input
  ! and as the output.
  subroutine iteration(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,Ntotal_plus,Ntotal_minus,&
       Indof2ndPhonon_plus,Indof3rdPhonon_plus,Indof2ndPhonon_minus,Indof3rdPhonon_minus,omega,&
       velocity,Gamma_plus,Gamma_minus,tau_zero,F_n)
    implicit none

    integer(kind=4),intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    integer(kind=4),intent(in) :: TypeofSymmetry(Nsymm_rot,nptk)
    integer(kind=4),intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands),Ntotal_plus,Ntotal_minus
    integer(kind=4),intent(in) :: Indof2ndPhonon_plus(Ntotal_plus),Indof3rdPhonon_plus(Ntotal_plus)
    integer(kind=4),intent(in) :: Indof2ndPhonon_minus(Ntotal_minus),Indof3rdPhonon_minus(Ntotal_minus)
    real(kind=8),intent(in) :: omega(nptk,nbands),velocity(nptk,nbands,3)
    real(kind=8),intent(in) :: Gamma_plus(Ntotal_plus),Gamma_minus(Ntotal_minus),tau_zero(nbands,nlist)
    real(kind=8),intent(inout) :: F_n(Nbands,nptk,3)

    integer(kind=4) :: ID_equi(Nsymm_rot,nptk),Naccum_plus,Naccum_minus
    integer(kind=4) :: i,j,k,jj,kk,ll,mm,nn
    real(kind=8) :: DeltaF(Nbands,nptk,3)

    call symmetry_map(ID_equi)
    DeltaF=0.d0
    do ll=1,Nlist
       do i=1,Nbands
          if (((ll-1)*Nbands+i).eq.1) then
             Naccum_plus=0
             Naccum_minus=0
          else
             Naccum_plus=Naccum_plus+N_plus((ll-1)*Nbands+i-1)
             Naccum_minus=Naccum_minus+N_minus((ll-1)*Nbands+i-1)
          end if
          do kk=1,Nequi(ll)
             if ((N_plus((ll-1)*Nbands+i).ne.0)) then
                do jj=1,N_plus((ll-1)*Nbands+i)
                   j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                   mm=int((Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                   k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                   nn=int((Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                   DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                        Gamma_plus(Naccum_plus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:)-&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:))
                end do !jj
             end if
             if ((N_minus((ll-1)*Nbands+i).ne.0)) then
                do jj=1,N_minus((ll-1)*Nbands+i)
                   j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                   mm=int((Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                   k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                   nn=int((Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                   DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                        Gamma_minus(Naccum_minus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:)+&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:))*5.D-1
                end do !jj
             end if
             F_n(i,ALLEquiList(kk,ll),:)=tau_zero(i,ll)*velocity(ALLEquiList(kk,ll),i,:)*&
                  omega(ALLEquiList(kk,ll),i)+tau_zero(i,ll)*DeltaF(i,ALLEquiList(kk,ll),:)
             F_n(:,ALLEquiList(kk,ll),:)=transpose(matmul(symmetrizers(:,:,ALLEquiList(kk,ll)),transpose(F_n(:,ALLEquiList(kk,ll),:))))
          end do !kk
       end do
    end do
  end subroutine iteration

    ! Four-phonon iteration with Gauss-Seidel approach
  subroutine iteration_4ph(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,Ntotal_plus,Ntotal_minus,&
       N_plusplus,N_plusminus,N_minusminus,Ntotal_plusplus,Ntotal_plusminus,Ntotal_minusminus,&
       Indof2ndPhonon_plus,Indof3rdPhonon_plus,Indof2ndPhonon_minus,Indof3rdPhonon_minus,&
       Indof2ndPhonon_plusplus,Indof3rdPhonon_plusplus,Indof4thPhonon_plusplus,&
       Indof2ndPhonon_plusminus,Indof3rdPhonon_plusminus,Indof4thPhonon_plusminus,&
       Indof2ndPhonon_minusminus,Indof3rdPhonon_minusminus,Indof4thPhonon_minusminus,&
       omega,velocity,Gamma_plus,Gamma_minus,Gamma_plusplus,Gamma_plusminus,Gamma_minusminus,tau_zero,F_n)
 
    implicit none 
   
    integer(kind=4),intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    integer(kind=4),intent(in) :: TypeofSymmetry(Nsymm_rot,nptk)
    integer(kind=4),intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands),Ntotal_plus,Ntotal_minus
    integer(kind=4),intent(in) :: Indof2ndPhonon_plus(Ntotal_plus),Indof3rdPhonon_plus(Ntotal_plus)
    integer(kind=4),intent(in) :: Indof2ndPhonon_minus(Ntotal_minus),Indof3rdPhonon_minus(Ntotal_minus)
   
    integer(kind=8),intent(in) :: N_plusplus(Nlist*Nbands),N_plusminus(Nlist*Nbands),N_minusminus(Nlist*Nbands)
    integer(kind=8),intent(in) :: Ntotal_plusplus,Ntotal_plusminus,Ntotal_minusminus
    integer(kind=8),intent(in) :: Indof2ndPhonon_plusplus(Ntotal_plusplus),Indof3rdPhonon_plusplus(Ntotal_plusplus),Indof4thPhonon_plusplus(Ntotal_plusplus)
    integer(kind=8),intent(in) :: Indof2ndPhonon_plusminus(Ntotal_plusminus),Indof3rdPhonon_plusminus(Ntotal_plusminus),Indof4thPhonon_plusminus(Ntotal_plusminus)
    integer(kind=8),intent(in) :: Indof2ndPhonon_minusminus(Ntotal_minusminus),Indof3rdPhonon_minusminus(Ntotal_minusminus),Indof4thPhonon_minusminus(Ntotal_minusminus)
   
    real(kind=8),intent(in) :: omega(nptk,nbands),velocity(nptk,nbands,3)
    real(kind=8),intent(in) :: Gamma_plus(Ntotal_plus),Gamma_minus(Ntotal_minus),tau_zero(nbands,nlist)
    real(kind=8),intent(in) :: Gamma_plusplus(Ntotal_plusplus),Gamma_plusminus(Ntotal_plusminus),Gamma_minusminus(Ntotal_minusminus)
    real(kind=8),intent(inout) :: F_n(Nbands,nptk,3)
   
    integer(kind=4) :: ID_equi(Nsymm_rot,nptk),Naccum_plus,Naccum_minus
    integer(kind=8) :: Naccum_plusplus,Naccum_plusminus,Naccum_minusminus
    integer(kind=4) :: i,j,k,l,jj,kk,ll,mm,nn,ss
    real(kind=8) :: DeltaF(Nbands,nptk,3)
   
    call symmetry_map(ID_equi)
    DeltaF=0.d0
   
    do ll=1,Nlist
      do i=1,Nbands
         if (((ll-1)*Nbands+i).eq.1) then
            Naccum_plus=0
            Naccum_minus=0
            Naccum_plusplus=0
            Naccum_plusminus=0
            Naccum_minusminus=0
         else
            Naccum_plus=Naccum_plus+N_plus((ll-1)*Nbands+i-1)
            Naccum_minus=Naccum_minus+N_minus((ll-1)*Nbands+i-1)
            Naccum_plusplus=Naccum_plusplus+N_plusplus((ll-1)*Nbands+i-1)
            Naccum_plusminus=Naccum_plusminus+N_plusminus((ll-1)*Nbands+i-1)
            Naccum_minusminus=Naccum_minusminus+N_minusminus((ll-1)*Nbands+i-1)
         end if
         do kk=1,Nequi(ll)
            if ((N_plus((ll-1)*Nbands+i).ne.0)) then
               do jj=1,N_plus((ll-1)*Nbands+i)
                  j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                  mm=int((Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                  k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                  nn=int((Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                  DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                       Gamma_plus(Naccum_plus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:)-&
                       F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:))
               end do !jj
            end if
            if ((N_minus((ll-1)*Nbands+i).ne.0)) then
               do jj=1,N_minus((ll-1)*Nbands+i)
                  j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                  mm=int((Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                  k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                  nn=int((Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                  DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                       Gamma_minus(Naccum_minus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:)+&
                       F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:))*5.D-1
               end do !jj
            end if
            if ((N_plusplus((ll-1)*Nbands+i).ne.0)) then
               do jj=1,N_plusplus((ll-1)*Nbands+i)
                  j=modulo(Indof2ndPhonon_plusplus(Naccum_plusplus+jj)-1,Nbands)+1
                  mm=int((Indof2ndPhonon_plusplus(Naccum_plusplus+jj)-1)/Nbands)+1
                  k=modulo(Indof3rdPhonon_plusplus(Naccum_plusplus+jj)-1,Nbands)+1
                  nn=int((Indof3rdPhonon_plusplus(Naccum_plusplus+jj)-1)/Nbands)+1
                  l=modulo(Indof4thPhonon_plusplus(Naccum_plusplus+jj)-1,Nbands)+1
                  ss=int((Indof4thPhonon_plusplus(Naccum_plusplus+jj)-1)/Nbands)+1
                  DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                        Gamma_plusplus(Naccum_plusplus+jj)*(F_n(l,ID_equi(TypeofSymmetry(kk,ll),ss),:)-&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:)-F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:))
               end do
            end if 
            if ((N_plusminus((ll-1)*Nbands+i).ne.0)) then
               do jj=1,N_plusminus((ll-1)*Nbands+i)
                  j=modulo(Indof2ndPhonon_plusminus(Naccum_plusminus+jj)-1,Nbands)+1
                  mm=int((Indof2ndPhonon_plusminus(Naccum_plusminus+jj)-1)/Nbands)+1
                  k=modulo(Indof3rdPhonon_plusminus(Naccum_plusminus+jj)-1,Nbands)+1
                  nn=int((Indof3rdPhonon_plusminus(Naccum_plusminus+jj)-1)/Nbands)+1
                  l=modulo(Indof4thPhonon_plusminus(Naccum_plusminus+jj)-1,Nbands)+1
                  ss=int((Indof4thPhonon_plusminus(Naccum_plusminus+jj)-1)/Nbands)+1
                  DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                        Gamma_plusminus(Naccum_plusminus+jj)*(F_n(l,ID_equi(TypeofSymmetry(kk,ll),ss),:)-&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:)+F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:))
               end do
            end if
            if ((N_minusminus((ll-1)*Nbands+i).ne.0)) then
               do jj=1,N_minusminus((ll-1)*Nbands+i)
                  j=modulo(Indof2ndPhonon_minusminus(Naccum_minusminus+jj)-1,Nbands)+1
                  mm=int((Indof2ndPhonon_minusminus(Naccum_minusminus+jj)-1)/Nbands)+1
                  k=modulo(Indof3rdPhonon_minusminus(Naccum_minusminus+jj)-1,Nbands)+1
                  nn=int((Indof3rdPhonon_minusminus(Naccum_minusminus+jj)-1)/Nbands)+1
                  l=modulo(Indof4thPhonon_minusminus(Naccum_minusminus+jj)-1,Nbands)+1
                  ss=int((Indof4thPhonon_minusminus(Naccum_minusminus+jj)-1)/Nbands)+1
                  DeltaF(i,ALLEquiList(kk,ll),:)=DeltaF(i,ALLEquiList(kk,ll),:)+&
                        Gamma_minusminus(Naccum_minusminus+jj)*(F_n(l,ID_equi(TypeofSymmetry(kk,ll),ss),:)+&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm),:)+F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn),:))
               end do
            end if
            F_n(i,ALLEquiList(kk,ll),:)=tau_zero(i,ll)*velocity(ALLEquiList(kk,ll),i,:)*&
                  omega(ALLEquiList(kk,ll),i)+tau_zero(i,ll)*DeltaF(i,ALLEquiList(kk,ll),:)
            F_n(:,ALLEquiList(kk,ll),:)=transpose(matmul(symmetrizers(:,:,ALLEquiList(kk,ll)),transpose(F_n(:,ALLEquiList(kk,ll),:))))
         end do
      end do
    end do 
  
  
  end subroutine iteration_4ph


  ! Restricted variation of the above subroutine, limited to cases where kappa
  ! is an scalar. Used in the nanowire calculation.
  subroutine iteration_scalar(Nlist,Nequi,ALLEquiList,TypeofSymmetry,N_plus,N_minus,Ntotal_plus,&
       Ntotal_minus,Indof2ndPhonon_plus,Indof3rdPhonon_plus,Indof2ndPhonon_minus,&
       Indof3rdPhonon_minus,omega,velocity,Gamma_plus,Gamma_minus,tau_zero,F_n)
    implicit none

    integer(kind=4),intent(in) :: Nlist,Nequi(Nlist),ALLEquiList(Nsymm_rot,nptk)
    integer(kind=4),intent(in) :: TypeofSymmetry(Nsymm_rot,nptk)
    integer(kind=4),intent(in) :: N_plus(Nlist*Nbands),N_minus(Nlist*Nbands),Ntotal_plus,Ntotal_minus
    integer(kind=4),intent(in) :: Indof2ndPhonon_plus(Ntotal_plus),Indof3rdPhonon_plus(Ntotal_plus)
    integer(kind=4),intent(in) :: Indof2ndPhonon_minus(Ntotal_minus),Indof3rdPhonon_minus(Ntotal_minus)
    real(kind=8),intent(in) :: omega(nptk,nbands),velocity(nptk,nbands)
    real(kind=8),intent(in) :: Gamma_plus(Ntotal_plus),Gamma_minus(Ntotal_minus),tau_zero(nbands,nlist)
    real(kind=8),intent(inout) :: F_n(Nbands,nptk)

    integer(kind=4) :: ID_equi(Nsymm_Rot,nptk),Naccum_plus,Naccum_minus
    integer(kind=4) :: i,j,k,jj,kk,ll,mm,nn
    real(kind=8) :: DeltaF(Nbands,nptk)

    DeltaF=0.d0
    call symmetry_map(ID_equi)
    DeltaF=0.d0
    do ll=1,Nlist
       do i=1,Nbands
          if (((ll-1)*Nbands+i).eq.1) then
             Naccum_plus=0
             Naccum_minus=0
          else
             Naccum_plus=Naccum_plus+N_plus((ll-1)*Nbands+i-1)
             Naccum_minus=Naccum_minus+N_minus((ll-1)*Nbands+i-1)
          end if
          do kk=1,Nequi(ll)
             if ((N_plus((ll-1)*Nbands+i).ne.0)) then
                do jj=1,N_plus((ll-1)*Nbands+i)
                   j=modulo(Indof2ndPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                   mm=int((Indof2ndPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                   k=modulo(Indof3rdPhonon_plus(Naccum_plus+jj)-1,Nbands)+1
                   nn=int((Indof3rdPhonon_plus(Naccum_plus+jj)-1)/Nbands)+1
                   DeltaF(i,ALLEquiList(kk,ll))=DeltaF(i,ALLEquiList(kk,ll))+&
                        Gamma_plus(Naccum_plus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn))-&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm)))
                end do !jj
             end if
             if ((N_minus((ll-1)*Nbands+i).ne.0)) then
                do jj=1,N_minus((ll-1)*Nbands+i)
                   j=modulo(Indof2ndPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                   mm=int((Indof2ndPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                   k=modulo(Indof3rdPhonon_minus(Naccum_minus+jj)-1,Nbands)+1
                   nn=int((Indof3rdPhonon_minus(Naccum_minus+jj)-1)/Nbands)+1
                   DeltaF(i,ALLEquiList(kk,ll))=DeltaF(i,ALLEquiList(kk,ll))+&
                        Gamma_minus(Naccum_minus+jj)*(F_n(k,ID_equi(TypeofSymmetry(kk,ll),nn))+&
                        F_n(j,ID_equi(TypeofSymmetry(kk,ll),mm)))*5.D-1
                end do !jj
             end if
             F_n(i,ALLEquiList(kk,ll))=tau_zero(i,ll)*velocity(ALLEquiList(kk,ll),i)*&
                  omega(ALLEquiList(kk,ll),i)+tau_zero(i,ll)*DeltaF(i,ALLEquiList(kk,ll))
          end do !kk
       end do
    end do
  end subroutine iteration_scalar


end module iterations
