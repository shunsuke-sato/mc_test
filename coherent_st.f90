program main
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0,1d0)
  integer,parameter :: ndim = 12
  complex(8) :: z1(ndim),z2(ndim)
  complex(8) :: z1_org(ndim),z2_org(ndim), zk(ndim)
  integer,parameter :: Nmc = 1000000
  integer :: imc
  complex(8) :: znorm,zs
  real(8)    :: sigma_norm
  real(8)    :: alpha
  integer,parameter :: nmax = 60
  real(8),parameter :: rmax = 20d0, dr = rmax/nmax
  real(8) :: sigma_r(0:nmax+1)
  integer(8) :: icount_r(0:nmax+1),ir



!Normal Monte-Carlo
  znorm = 0d0
  sigma_norm = 0d0
  sigma_r = 0d0
  icount_r = 0d0


  do imc = 1,nmc
    call init_coherent_st(z1,z2,ndim)
    ir = aint(sqrt(sum(abs(z1)**2))/dr)
    ir = min(ir ,nmax+1)
    icount_r(ir) = icount_r(ir) + 1

    zs = (4d0/3d0)**ndim*exp(-zI*aimag(sum(z1*conjg(z2))))
    znorm = znorm + zs
    sigma_norm = sigma_norm + abs(zs-1d0)**2
    sigma_r(ir) = sigma_r(ir) + abs(zs-1d0)**2
    
  end do
  znorm = znorm/nmc
  sigma_norm = sqrt(sigma_norm/nmc)

  write(*,"(A,2x,I9)")"nmc=",nmc
  write(*,"(A,2x,2e26.16e3)")"norm=",znorm
  write(*,"(A,2x,2e26.16e3)")"sigma=",sigma_norm

  open(20,file="sigma_r_correlated_gaussian_normal.out")
  do ir = 0,nmax+1
     if(icount_r(ir) /= 0)then
        sigma_r(ir) = sqrt(sigma_r(ir)/icount_r(ir))
     else
        sigma_r(ir) = 0d0
     end if
     write(20,"(999e16.6e3)")dble(ir)*dr,dble(icount_r(ir)),sigma_r(ir)
  end do
  close(20)

!Normal Monte-Carlo
  znorm = 0d0
  sigma_norm = 0d0
  sigma_r = 0d0
  icount_r = 0d0

  do imc = 1,nmc
    call init_coherent_st_double_gauss(z1,z2,ndim)
    z2 = z2 + 0.5d0*z1

    ir = aint(sqrt(sum(abs(z1)**2))/dr)
    ir = min(ir ,nmax+1)
    icount_r(ir) = icount_r(ir) + 1

    zs = exp(0.25d0*sum(abs(z1)**2)) *exp(-zI*aimag(sum(z1*conjg(z2))))
    znorm = znorm + zs
    sigma_norm = sigma_norm + abs(zs-1d0)**2
    sigma_r(ir) = sigma_r(ir) + abs(zs-1d0)**2
    
  end do
  znorm = znorm/nmc
  sigma_norm = sqrt(sigma_norm/nmc)

  write(*,"(A,2x,I9)")"nmc=",nmc
  write(*,"(A,2x,2e26.16e3)")"norm=",znorm
  write(*,"(A,2x,2e26.16e3)")"sigma=",sigma_norm

  open(20,file="sigma_r_independent_gaussian_normal.out")
  do ir = 0,nmax+1
     if(icount_r(ir) /= 0)then
        sigma_r(ir) = sqrt(sigma_r(ir)/icount_r(ir))
     else
        sigma_r(ir) = 0d0
     end if
     write(20,"(999e16.6e3)")dble(ir)*dr,dble(icount_r(ir)),sigma_r(ir)
  end do
  close(20)


! Replica-shift Monte-Carlo
  znorm = 0d0
  sigma_norm = 0d0
  sigma_r = 0d0
  icount_r = 0d0
  do imc = 1,nmc
    call init_coherent_st_double_gauss(z1,z2,ndim)
    z2 = z2 + 0.5d0*z1

    ir = aint(sqrt(sum(abs(z1)**2))/dr)
    ir = min(ir ,nmax+1)
    icount_r(ir) = icount_r(ir) + 1

    zs = 0.5d0 * exp(0.25d0*sum(abs(z1)**2)) *exp(-zI*aimag(sum(z1*conjg(z2))))

    z2_org = z2

    zk = aimag(z1) - zI* real(z1)
    alpha = pi/sqrt(sum(abs(zk)**2)) ; alpha = alpha*exp(-0.5d0*alpha**2)

    z2 = z2_org + alpha*zk/sqrt(sum(abs(zk)**2))
    zs = zs + 0.5d0 *0.5d0* exp(0.25d0*sum(abs(z1)**2)) *exp(-zI*aimag(sum(z1*conjg(z2)))) &
         *exp(sum(-abs(z2-0.5d0*z1)**2 +abs(z2_org -0.5d0*z1)**2))

    z2 = z2_org - alpha*zk/sqrt(sum(abs(zk)**2))
    zs = zs + 0.5d0 *0.5d0* exp(0.25d0*sum(abs(z1)**2)) *exp(-zI*aimag(sum(z1*conjg(z2)))) &
         *exp(sum(-abs(z2-0.5d0*z1)**2 +abs(z2_org -0.5d0*z1)**2))


    znorm = znorm + zs
    sigma_norm = sigma_norm + abs(zs-1d0)**2
    sigma_r(ir) = sigma_r(ir) + abs(zs-1d0)**2
    
  end do
  znorm = znorm/nmc
  sigma_norm = sqrt(sigma_norm/nmc)

  write(*,"(A,2x,I9)")"nmc=",nmc
  write(*,"(A,2x,2e26.16e3)")"norm=",znorm
  write(*,"(A,2x,2e26.16e3)")"sigma=",sigma_norm

  open(20,file="sigma_r_independent_gaussian_shift.out")
  do ir = 0,nmax+1
     if(icount_r(ir) /= 0)then
        sigma_r(ir) = sqrt(sigma_r(ir)/icount_r(ir))
     else
        sigma_r(ir) = 0d0
     end if
     write(20,"(999e16.6e3)")dble(ir)*dr,dble(icount_r(ir)),sigma_r(ir)
  end do
  close(20)

  if(.true.)then

! Replica shift Monte-Carlo
  znorm = 0d0
  sigma_norm = 0d0
  sigma_r = 0d0
  icount_r = 0d0
  do imc = 1,nmc
    call init_coherent_st(z1,z2,ndim)
    ir = aint(sqrt(sum(abs(z1)**2))/dr)
    ir = min(ir ,nmax+1)
    icount_r(ir) = icount_r(ir) + 1

    zs = 0.5d0*(4d0/3d0)**ndim*exp(-zI*aimag(sum(z1*conjg(z2))))
    z1_org = z1
    z2_org = z2

! shift for z1
    zk = -aimag(z2) + zI* real(z2)
    alpha = pi/sqrt(sum(abs(zk)**2)) ; alpha = alpha*exp(-0.3d0*alpha)
    z1 = z1_org + alpha*zk/sqrt(sum(abs(zk)**2))
    zs = zs + 0.5d0*0.25d0*(4d0/3d0)**ndim*exp(-zI*aimag(sum(z1*conjg(z2)))) &
      *exp(sum(-0.5d0*abs(z1)**2 -0.5d0*abs(z2)**2 -0.5d0*abs(z1-z2)**2 &
      +0.5d0*abs(z1_org)**2 +0.5d0*abs(z2_org)**2 +0.5d0*abs(z1_org-z2_org)**2))

    z1 = z1_org - alpha*zk/sqrt(sum(abs(zk)**2))
    zs = zs + 0.5d0*0.25d0*(4d0/3d0)**ndim*exp(-zI*aimag(sum(z1*conjg(z2)))) &
      *exp(sum(-0.5d0*abs(z1)**2 -0.5d0*abs(z2)**2 -0.5d0*abs(z1-z2)**2 &
      +0.5d0*abs(z1_org)**2 +0.5d0*abs(z2_org)**2 +0.5d0*abs(z1_org-z2_org)**2 ))

    z1 = z1_org

! shift for z2
    zk = aimag(z1) - zI* real(z1)
    alpha = pi/sqrt(sum(abs(zk)**2)) ; alpha = alpha*exp(-0.3d0*alpha)

    z2 = z2_org + alpha*zk/sqrt(sum(abs(zk)**2))
    zs = zs + 0.5d0*0.25d0*(4d0/3d0)**ndim*exp(-zI*aimag(sum(z1*conjg(z2)))) &
      *exp(sum(-0.5d0*abs(z1)**2 -0.5d0*abs(z2)**2 -0.5d0*abs(z1-z2)**2 &
      +0.5d0*abs(z1_org)**2 +0.5d0*abs(z2_org)**2 +0.5d0*abs(z1_org-z2_org)**2 ))


    z2 = z2_org - alpha*zk/sqrt(sum(abs(zk)**2))
    zs = zs + 0.5d0*0.25d0*(4d0/3d0)**ndim*exp(-zI*aimag(sum(z1*conjg(z2)))) &
      *exp(sum(-0.5d0*abs(z1)**2 -0.5d0*abs(z2)**2 -0.5d0*abs(z1-z2)**2 &
      +0.5d0*abs(z1_org)**2 +0.5d0*abs(z2_org)**2 +0.5d0*abs(z1_org-z2_org)**2 ))




    znorm = znorm + zs
    sigma_norm = sigma_norm + abs(zs-1d0)**2
    sigma_r(ir) = sigma_r(ir) + abs(zs-1d0)**2
    
  end do
  znorm = znorm/nmc
  sigma_norm = sqrt(sigma_norm/nmc)

  write(*,"(A,2x,I9)")"nmc=",nmc
  write(*,"(A,2x,2e26.16e3)")"norm=",znorm
  write(*,"(A,2x,2e26.16e3)")"sigma=",sigma_norm

  open(20,file="sigma_r_correlated_gaussian_shift.out")
  do ir = 0,nmax+1
     if(icount_r(ir) /= 0)then
        sigma_r(ir) = sqrt(sigma_r(ir)/icount_r(ir))
     else
        sigma_r(ir) = 0d0
     end if
     write(20,"(999e16.6e3)")dble(ir)*dr,dble(icount_r(ir)),sigma_r(ir)
  end do
  close(20)

end if

end program main
subroutine init_coherent_st(z1,z2,ndim)
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  integer, intent(in) :: ndim
  complex(8),intent(out) :: z1(ndim),z2(ndim)
  integer :: idim
  real(8) :: x1,x2,p1,p2

  do idim = 1,ndim

    call correlated_gaussian_random_number(x1,x2)
    call correlated_gaussian_random_number(p1,p2)
    z1(idim) = x1 + zI*p1
    z2(idim) = x2 + zI*p2

  end do

end subroutine init_coherent_st

subroutine init_coherent_st_double_gauss(z1,z2,ndim)
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  integer, intent(in) :: ndim
  complex(8),intent(out) :: z1(ndim),z2(ndim)
  integer :: idim
  real(8) :: x1,x2,p1,p2

  do idim = 1,ndim
     call gaussian_random_number(x1,p1)
     call gaussian_random_number(x2,p2)
     z1(idim) = (x1 + zI*p1)/sqrt(2d0)
     z2(idim) = (x2 + zI*p2)/sqrt(2d0)
  end do

end subroutine init_coherent_st_double_gauss

subroutine correlated_gaussian_random_number(x1,x2)
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  real(8) :: r1,r2,tmp, sigma
  
  sigma = 1d0
  
  if(sigma >= 1d0)then
    
    do 
      call random_number(r1)
      call random_number(r2)
      
      if(r1 == 0d0)then
        x1 = 0d0
        x2 = 0d0
      else 
        tmp = sqrt(-2d0*log(r1))
        x1 = tmp*cos(2d0*pi*r2)
        x2 = tmp*sin(2d0*pi*r2)
      end if
      
      tmp = x1 -x2
      r1 = exp(-0.5d0*tmp**2/sigma)
      call random_number(r2)
      if(r2 < r1)exit
    end do
    
  else
    
    do 
      call random_number(r1)
      call random_number(r2)
      
      if(r1 == 0d0)then
        x1 = 0d0
        x2 = 0d0
      else 
        tmp = sqrt(-2d0*log(r1))
        x1 = tmp*cos(2d0*pi*r2)
        x2 = tmp*sin(2d0*pi*r2)
      end if
      
      x2 = x1 + x2*sqrt(sigma)
      r1 = exp(-0.5d0*x2**2)
      call random_number(r2)
      if(r2 < r1)exit
    end do
  end if
  
end subroutine correlated_gaussian_random_number

subroutine gaussian_random_number(x1,x2)
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  integer :: len = 2
  real(8) :: rvec(2),tmp
  real(8) :: r1,r2

  call random_number(r1)
  call random_number(r2)
  rvec(1) = r1
  rvec(2) = r2
!
  if(rvec(1) == 0d0)then
    x1 = 0d0
    x2 = 0d0
  else 
    tmp = sqrt(-2d0*log(rvec(1)))
    x1 = tmp*cos(2d0*pi*rvec(2))
    x2 = tmp*sin(2d0*pi*rvec(2))
  end if

end subroutine gaussian_random_number
