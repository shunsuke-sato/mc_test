module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  integer,parameter :: ndim = 10
  real(8),parameter :: dt = 0.0d0
  real(8) :: x_ho(ndim)
  integer,parameter :: Niter_MC = 10000

end module global_variables
program main
  use global_variables
  implicit none
  integer,parameter :: Nx = 256
  integer,parameter :: nmax = 0
  real(8),parameter :: xmin = -6d0,xmax= 6d0,dx = (xmax-xmin)/Nx
  integer :: ix,iter
  real(8) :: xx
  complex(8) :: zs1,zs2
  complex(8) :: propagator_HO
  complex(8) :: zf, gg
  complex(8) :: zs
  real(8)    :: ss


  open(20,file='out.out')
  do ix = 0,Nx
    xx = xmin + dx*ix
    zs1 = propagator_HO(0d0,xx,nmax,dt)
    zs2 = propagator_HO(1d0,xx,nmax,dt)

    write(20,"(999e26.16e3)")xx,zs1,zs2

  end do
  close(20)
  
! Monte Carlo
  x_ho = 0d0
  zs = 0d0
  ss = 0d0
  do iter = 1,Niter_MC
    write(*,*)iter
    call Metropolis_update(x_ho,ndim,dt,nmax)

    zf  = exp(-0.5d0*(x_ho(1)/1d0)**2)*exp(-0.5d0*(x_ho(10)/1d0)**2)/(sqrt(pi)*1d0) &
      *propagator_HO(x_ho(2),x_ho( 1),nmax,dt) &
      * propagator_HO(x_ho(3),x_ho( 2),nmax,dt) &
      * propagator_HO(x_ho(4),x_ho( 3),nmax,dt) &
      * propagator_HO(x_ho(5),x_ho( 4),nmax,dt) &
      * propagator_HO(x_ho(7),x_ho( 6),nmax,dt) &
      * propagator_HO(x_ho(8),x_ho( 7),nmax,dt) &
      * propagator_HO(x_ho(9),x_ho( 8),nmax,dt) &
      * propagator_HO(x_ho(10),x_ho(9),nmax,dt)


    zs = zs + zf/abs(zf) &
      * exp(-0.5d0*(x_ho(5)/1d0)**2)*exp(-0.5d0*(x_ho(6)/1d0)**2)/(sqrt(pi)*1d0)

    ss = ss + exp(-0.5d0*sum(x_ho(:)**2))/abs(zf)/(sqrt(2d0*pi))**10

  end do
  zs = zs/Niter_MC
  ss = ss/Niter_MC

  write(*,*)"result",zs,ss
  write(*,*)"result",zs/ss
  
end program main

subroutine Metropolis_update(x_ho,ndim,dt,nmax)
  implicit none
  integer,parameter :: nskip = 1000
  integer :: ndim,nmax
  real(8) :: x_ho(ndim),x_ho_old(ndim)
  real(8) :: dt
  integer :: n_accept,niter
  real(8)  :: weight,weight_old,r

  niter = 0
  n_accept = 0

  call calc_weight(x_ho,ndim,weight,dt)
  weight_old = weight
  x_ho_old   = x_ho

  do 

    niter = niter + 1
    call renew_x_ho(x_ho,x_ho_old,ndim)
    call calc_weight(x_ho,ndim,weight,dt)
    call random_number(r)
    if(r <weight/weight_old)then
      n_accept = n_accept + 1
      x_ho_old = x_ho
      weight_old = weight
    end if

    if(niter >= nskip)exit
!    if(n_accept == 100)exit
  end do

!  write(*,*)"niter=",niter
!  write(*,*)"n_accept=",n_accept
!  write(*,*)"acceptance ratio=",100d0*dble(n_accept)/niter,"[%]"
!
  contains
    subroutine renew_x_ho(x_ho,x_ho_old,ndim)
      implicit none
      integer :: ndim
      real(8) :: x_ho(ndim),x_ho_old(ndim)
      real(8) :: x1,x2
      integer :: n

      do n = 1,ndim
        call gaussian_random_number(x1,x2)
        x_ho(n) = x_ho_old(n) + x1
      end do


      
    end subroutine renew_x_ho

    subroutine calc_weight(x_ho,ndim,weight,dt)
      implicit none
      integer :: ndim
      real(8) :: dt
      real(8) :: x_ho(ndim),weight
      complex(8) :: propagator_HO

      weight = exp(-0.5d0*(x_ho(1)/1d0)**2)*exp(-0.5d0*(x_ho(10)/1d0)**2)
      weight = weight * abs( &
                        propagator_HO(x_ho(2),x_ho( 1),nmax,dt) &
                      * propagator_HO(x_ho(3),x_ho( 2),nmax,dt) &
                      * propagator_HO(x_ho(4),x_ho( 3),nmax,dt) &
                      * propagator_HO(x_ho(5),x_ho( 4),nmax,dt) &
                      * propagator_HO(x_ho(6),x_ho( 7),nmax,dt) &
                      * propagator_HO(x_ho(7),x_ho( 8),nmax,dt) &
                      * propagator_HO(x_ho(8),x_ho( 9),nmax,dt) &
                      * propagator_HO(x_ho(9),x_ho(10),nmax,dt))



    end subroutine calc_weight

end subroutine Metropolis_update

function propagator_HO(x1,x2,nmax,dt) result(z)
  implicit  none
  complex(8),parameter :: zI = (0d0,1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8) :: x1,x2,dt
  integer :: nmax
  complex(8) :: z
  integer :: n
  real(8) :: ss,eps
  real(8) :: Hermite_polynomials

  z = 0d0
  ss = 1d0
  do n = 0,nmax
    eps = dble(n) + 0.5d0

    ss = ss*2d0**min(n,1)*max(n,1)

    z = z + Hermite_polynomials(x1,n)*Hermite_polynomials(x2,n)/ss &
      *sqrt(1d0/pi)*exp(-0.5d0*x1**2)*exp(-0.5d0*x2**2) &
      *exp(-zI*eps*dt)

  end do


end function propagator_HO

function Hermite_polynomials(x,n) result(y)
  implicit none
  real(8) :: x,y
  integer :: n

  select case(n)
  case(0)
    y = 1d0
  case(1)
    y = 2d0*x
  case(2)
    y = 4d0*x**2-2d0
  case(3)
    y = 8d0*x**3-12d0*x
  case(4)
    y = 16d0*x**4-48d0*x**2+12d0
  case(5)
    y = 32d0*x**5 - 160d0*x**3 + 120d0*x
  case(6)
    y = 64d0*x**6 - 480d0*x**4 + 720d0*x**2 -120d0
  case(7)
    y = 128d0*x**7 - 1344d0*x**5 + 3360d0*x**3 - 1680d0*x
  case(8)
    y = 256d0*x**8 - 3584d0*x**6 + 13440d0*x**4 - 13440d0*x**2 + 1680d0
  case(9)
    y = 512d0*x**9 - 9216d0*x**7 + 48384d0*x**5 - 80640d0*x**3 + 30240d0*x
  case(10)
    y = 1024d0*x**10 - 23040d0*x**8 + 161280d0*x**6 + 302400*x**2 - 30240d0
  case default
    stop 'Hermite polynomials is not implemented for n>10.'
  end select
      
end function Hermite_polynomials

subroutine gaussian_random_number(x1,x2)
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  integer :: len = 2
  real(8) :: r1,r2,tmp


  call random_number(r1)
  call random_number(r2)
!
  if(r1 == 0d0)then
    x1 = 0d0
    x2 = 0d0
  else 
    tmp = sqrt(-2d0*log(r1))
    x1 = tmp*cos(2d0*pi*r2)
    x2 = tmp*sin(2d0*pi*r2)
  end if

end subroutine gaussian_random_number
