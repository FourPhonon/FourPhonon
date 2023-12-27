! Some miscellaneous math functions.
module misc
  use iso_fortran_env
  use iso_c_binding
  implicit none

  ! Create a directory (interface to the system's mkdir).
  interface
     function mkdir(path,mode) bind(c,name="mkdir")
       use iso_c_binding
       implicit none
       character(kind=c_char,len=1) :: path(*)
       integer(c_int16_t), value :: mode
       integer(c_int) :: mkdir
     end function mkdir

     function syschdir(path) bind(c,name="chdir")
       use iso_c_binding
       implicit none
       character(kind=c_char,len=1) :: path(*)
       integer(c_int) :: syschdir
     end function syschdir
  end interface

contains

  ! 3D cross product.
  subroutine cross_product(a,b,res)
    real(kind=8),intent(in) :: a(3),b(3)
    real(kind=8),intent(out) :: res(3)

    integer(kind=4) :: i,j,k

    do i=1,3
       j=mod(i,3)+1
       k=mod(j,3)+1
       res(i)=a(j)*b(k)-a(k)*b(j)
    end do
  end subroutine cross_product

  ! Quotient and remainder of integer division.
  subroutine divmod(a,b,q,r)
    implicit none

    integer(kind=4),intent(in) :: a,b
    integer(kind=4),intent(out) :: q,r

    q=a/b
    r=mod(a,b)
  end subroutine divmod

  ! 2-norm of a 3x3 matrix.
  function twonorm3x3(a)
    implicit none

    real(kind=8),intent(in) :: a(3,3)

    real(kind=8) :: twonorm3x3

    integer(kind=4) :: info,iwork(24)
    real(kind=8) :: b(3,3),S(3),work(100)

    b=a
    call dgesdd("N",3,3,b,3,S,b,3,b,3,work,100,iwork,info)
    twonorm3x3=S(1)
  end function twonorm3x3

  ! exp(iunit*x). Sometimes it is much faster than the general
  ! complex exponential.
  function phexp(x)
    implicit none

    real(kind=8),intent(in) :: x
    complex(kind=8) :: phexp

    phexp=cmplx(cos(x),sin(x),kind=8)
  end function phexp

  ! Create a directory with permissions set to 0777.
  subroutine create_directory(path,result)
    implicit none
    character(kind=c_char,len=1),intent(in) :: path(*)
    integer(kind=4),intent(out),optional :: result

    integer(kind=4) :: tmp
    tmp=int(mkdir(path, int(o"777",c_int16_t)),4)
    if (present(result)) then
       result = tmp
    end if
  end subroutine create_directory

  ! Change the working directory
  subroutine change_directory(path,result)
    implicit none
    character(kind=c_char,len=1),intent(in) :: path(*)
    integer(kind=4),intent(out),optional :: result

    integer(kind=4) :: tmp
    tmp=int(syschdir(path),4)
    if (present(result)) then
       result = tmp
    end if
  end subroutine change_directory
end module misc
