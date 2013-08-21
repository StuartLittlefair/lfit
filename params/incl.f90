program incl
  
  implicit none

  real :: q, dphi, qerr, dphierr
  real :: qsave, dphisave
  real :: r2_a
  real :: r2_a_sq, sin_thing, cos_thing
  real :: cossq_i, i
  integer :: k
  real, dimension(1000) :: rstore, istore 
  real, parameter :: pi=3.141
  real, external :: gaussran
  real :: mean, var

  print*, '> Enter mass ratio and err: '
  read(*,*) q, qerr
  print*, '> Enter phase width of white dwarf ingress and error: '
  read(*, *) dphi, dphierr
  qsave = q
  dphisave = dphi
  do k=1,1000
     q = qsave + gaussran()*qerr
     dphi = dphisave + gaussran()*dphierr

     r2_a = 0.49*q**(2./3.) / ( 0.6*q**(2./3.) + log(1.0+q**(1./3.)) )
     r2_a_sq = r2_a **2.0
     sin_thing = (sin(pi*dphi)**2.0)
     cos_thing = (cos(pi*dphi)**2.0)

     cossq_i = (r2_a_sq-sin_thing)/cos_thing
     i = sqrt(abs(cossq_i))
     i = acos(i)
     i = i*180/pi
     
     rstore(k) = r2_a
     istore(k) = i
  end do

!  call pgopen('?')
!  call pgsubp(1,2)
!  call pghist(1000,istore,minval(istore),maxval(istore),100,0)
!  call pghist(1000,rstore,minval(rstore),maxval(rstore),100,0)
!  call pgclos

  mean = sum(rstore)/real(size(rstore))
  var = 0.0
  do k=1,1000
     var = var + (rstore(k)-mean)**2.0
  end do
  var = var/size(rstore)

  print*, 'R_2/a: ',  mean, '+/- ', sqrt(var)

  mean = sum(istore)/real(size(istore))
  var = 0.0
  do k=1,1000
     var = var + (istore(k)-mean)**2.0
  end do
  var = var/size(istore)

  print*, 'Inclination: ', mean, '+/- ', sqrt(var)

end program incl
  
real function gaussran()
  
  ! uses box-muller transformation to obtain a random number with a
  ! gaussian probability distribution, mean 0, and unit variance
  implicit none
  
  real :: v1, v2, rsq, fac
  
1 call random_number(v1)
  call random_number(v2)
  v1 = 2.0*v1-1.0
  v2 = 2.0*v2-1.0
  rsq = v1**2. + v2**2.
  if (rsq .ge. 1.0 .or. rsq .eq. 0.0) goto 1
  fac = sqrt(-2.0*log(rsq)/rsq)
  gaussran = v1*fac
  return
  
end function gaussran
