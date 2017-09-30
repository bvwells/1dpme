module linear_solvers

contains

  subroutine tridiagonal(upper,diag,lower,x,B,idim)

    implicit none
    !---------------------------------------------------------------------------------
    integer, intent(IN) :: idim
    double precision, dimension(1:idim) :: upper, diag, lower, B, X
    !---------------------------------------------------------------------------------
    double precision, dimension(1:idim) :: gamma
    double precision :: beta
    integer :: i
    !---------------------------------------------------------------------------------      
    gamma=0; beta=0

    if (dabs(diag(1))<1d-10) stop 'Rewrite equations: trisolve'

    beta = diag(1)
    x(1) = B(1)/beta

    ! Forward substitution

    do i = 1, idim
       gamma(i) = upper(i-1) / beta
       beta = diag(i) - lower(i)*gamma(i)

       if (dabs(beta)<1d-10) stop 'trisolve failed'

       x(i) = ( B(i) - lower(i)*x(i-1) ) / beta
    end do

    ! Backward substitution

    do i = idim-1, 1, -1
       x(i) = x(i) - gamma(i+1)*x(i+1)
    end do

    return

  end subroutine tridiagonal

  subroutine gaussian_elimination(mat,xsol,B,idim)

    implicit none
    !---------------------------------------------------------------------------------
    integer, intent(IN) :: idim
    double precision, intent(IN), dimension(1:idim,1:idim) :: mat
    double precision, intent(INOUT), dimension(1:idim) :: xsol
    double precision, intent(IN), dimension(1:idim) :: B
    !---------------------------------------------------------------------------------      
    double precision, dimension(1:idim,1:idim+1) :: a
    double precision :: pmax, aux
    integer :: i, j, k, idim1, ll
    !---------------------------------------------------------------------------------

    idim1= idim +1

    a=0.0d0

    do i=1,idim
       do j=1,idim1
          a(i,j)=mat(i,j)
       end do
    end do

    a(:,idim1) = B(:)

    do i=1,idim-1
       pmax = 0.0

       do j=i,idim
          if(abs(a(j,i)).gt.pmax) then
             pmax = abs(a(j,i))
             ll = j
          endif
       enddo

       if (ll.ne.i) then
          do j=i,idim1
             aux = a(i,j)
             a(i,j) = a(ll,j)
             a(ll,j) = aux
          enddo
       endif

       aux = a(i,i)

       do j =i+1,idim1
          a(i,j) = a(i,j)/aux
       end do

       do k=i+1,idim
          aux = a(k,i)

          do j=i+1,idim1
             a(k,j)=a(k,j) - aux*a(i,j)
          end do
       end do

    end do

    a(idim,idim1) = a(idim,idim1)/a(idim,idim)

    do i=idim-1,1,-1
       do j=i+1,idim
          a(i,idim1) = a(i,idim1) - a(i,j)*a(j,idim1)
       enddo
    enddo

    xsol = a(:,idim1)    

    return  

  end subroutine gaussian_elimination

end module linear_solvers
