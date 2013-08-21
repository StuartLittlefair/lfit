module params

  double precision, parameter :: G = 6.670d-11   ! Grav. Const Nm**2/kg**2
  double precision, parameter :: sig = 5.67d-8   ! stefans J/m**2/K**4
  double precision, parameter :: Msun = 1.989d30 ! kg
  double precision, parameter :: Rsun = 6.9599d8 ! m
  double precision, parameter :: Lsun = 3.826d26 ! W
  double precision, parameter :: pi = 3.141818d0
  double precision, parameter :: day_in_secs = 86400.0d0 
  double precision, parameter :: year_in_secs = 3.1556925e7

end module params

program stellar_calcs

  ! performs common stellar calculations

  use params
  implicit none
  integer :: choice


  write(*,*) '*************************************'
  write(*,*) '     Common Stellar Calculations     '
  write(*,*) '                                     '
  write(*,*) ' Pick an option:                     '
  write(*,*) '                                     '
  write(*,*) ' 1) Calculate disc luminosity and    '
  write(*,*) '    boundary layer temp, given       '
  write(*,*) '    Mdot, Rstar and central mass     '
  write(*,*) '                                     '
  write(*,*) ' 2) Calculate stellar radius given   '
  write(*,*) '    log of luminosity and temp       '
  write(*,*) '                                     '
  write(*,*) ' 3) Calculate Rot. Per. given vsini  '
  write(*,*) '    and stellar radius (assumes      '
  write(*,*) '    inclination of 90 degrees).      '
  write(*,*) '                                     '
  write(*,*) ' 4) Calculate Kelvin-Helmhotz time   '
  write(*,*) '                                     '     
  write(*,*) '                                     '     
  write(*,*) ' 5) Calculate Kelvin-Helmhotz time   '
  write(*,*) '    given log g, M and L             '
  write(*,*) '                                     '
  write(*,*) '                                     '
  write(*,*) ' 6) Calculate log g given M and R    '
  write(*,*) '                                     '     
  write(*,*) '                                     '
  write(*,*) ' 7) Calculate R given M and log g    '
  write(*,*) '                                     '     
  write(*,*) '                                     '

  write(*,'(A)',advance='no') 'Option: '
  read(*,*) choice

  select case(choice)
  case(1) 
     call disc_calc
  case(2)
     call radius_calc
  case(3)
     call period_calc
  case(4)
     call kh_calc
  case(5)
     call kh_calc2
  case(6)
     call logg_calc
  case(7)
     call logg_calc2
  end select

end program stellar_calcs

subroutine disc_calc

  use params
  implicit none
  double precision :: mdot,m1,rstar,delta_r_bl,v
  double precision :: L,T_bl



  write(*,'(A)',advance='no') 'Input m_dot (solar masses/yr): '
  read(*,*) mdot
  mdot = mdot*Msun/year_in_secs ! convert to SI units.
  write(*,'(A)',advance='no') &
       'Input central object mass (solar masses): '
  read(*,*) m1
  m1 = m1*Msun ! convert to SI units.
  write(*,'(A)',advance='no') &
       'Input central object radius (solar radii): '
  read(*,*) rstar
  rstar = rstar*Rsun ! convert to SI units.


  print*, ''
  print*, ''

  L = G*m1*mdot/2./rstar/Lsun
  print*, 'Disc Luminosity (solar units): ',L


  delta_r_bl = 0.02d0*rstar
  T_bl = (G*m1*mdot/rstar/sig/4./pi/rstar/delta_r_bl)**(1./4.)
  print*, 'Boundary Layer temperature (K): ',T_bl

  v = sqrt(G*m1/rstar)
  v = v/1000.d0
  print*, 'Inner disc velocity (km/s): ',v

end subroutine disc_calc

subroutine period_calc
  use params

  implicit none
  double precision :: rstar
  double precision :: vsini, P


  write(*,'(A)',advance='no') 'Input stellar radius (solar units): '
  read(*,*) rstar
  rstar = rstar*rsun ! convert to SI units.
  write(*,'(A)',advance='no') &
       'Input vsini (km/s): '
  read(*,*) vsini
  vsini = vsini*1000.0d0 ! convert to SI units

  print*,  '...assuming sin i = 1...'
  print*,  ''

  P = 2.0d0*pi*rstar/vsini/day_in_secs
  print*, 'rotation period, P < ', P

end subroutine period_calc
      
subroutine radius_calc
  use params

  implicit none

  double precision :: rstar
  double precision :: L,T

  write(*,'(A)',advance='no') 'Input log10 luminosity (solar units): '
  read(*,*) L
  L = 10.0d0**L
  L = L*Lsun ! convert to SI units.
  write(*,'(A)',advance='no') &
       'Input log10 effective temperature (K): '
  read(*,*) T
  T = 10.0d0**T

  print*,  ''
  print*,  ''

  rstar = sqrt(L/4.0d0/pi/sig/T**4.0d0)
  print*, 'radius of star =', rstar/rsun, ' solar radii'

end subroutine radius_calc

subroutine kh_calc
  use params

  implicit none

  double precision :: rstar,M
  double precision :: L,T,tkh

  write(*,'(A)',advance='no') 'Input log10 luminosity (solar units): '
  read(*,*) L
  L = 10.0d0**L
  L = L*Lsun ! convert to SI units.
  write(*,'(A)',advance='no') 'Input effective temperature (K): '
  read(*,*) T
  write(*,'(A)',advance='no') 'Input mass: '
  read(*,*) M
  M = M*Msun

  print*,  ''
  print*,  ''

  rstar = L/(T**4.0d0)
  rstar = rstar/4.0d0/pi/sig
  rstar = sqrt(rstar)
 
  tkh = G*M*M/L/Rstar

  print*, 'Kelvin Helmholtz timescale =', tkh/year_in_secs, ' yr'

end subroutine kh_calc

subroutine kh_calc2

  use params
  implicit none
  double precision rstar, mstar, L, logg, tkh

  write(*, '(A)', advance='no') 'Input mass of star (solar units): '
  read(*,*) mstar
  mstar = mstar*MSun*1000.0d0

  write(*, '(A)', advance='no') 'Input log g of star (cgs units): '
  read(*,*) logg
  
  rstar = sqrt ( G * 1.0d3 * mstar / (10.0d0**logg) )
  
  ! convert back to solar units
  mstar = mstar/MSun/1000.0d0
  rstar = rstar/Rsun/100.0d0
  
  ! convert to SI units
  mstar = mstar*Msun
  rstar = rstar*Rsun

  write(*,'(A)',advance='no') 'Input log10 luminosity (solar units): '
  read(*,*) L
  L = 10.0d0**L
  L = L*Lsun ! convert to SI units
  print*,  ''
  print*,  ''

  tkh = G*mstar*mstar/L/rstar
  print*, 'Kelvin Helmholtz timescale =', tkh/year_in_secs, ' yr'

end subroutine kh_calc2

subroutine logg_calc2
  use params
  implicit none

  double precision :: rstar, mstar, logg
  write(*,'(A)',advance='no') 'Input log g (cgs units)'
  read(*,*) logg
  write(*, '(A)', advance='no') 'Input mass of star (solar units): '
  read(*,*) mstar
  mstar = mstar*MSun*1000

  rstar = sqrt( (10.0**logg)/G/mstar/1.0e3 )
  rstar = 1.0/rstar
  rstar = rstar/Rsun/100
  print*, 'Radius (solar units) = ', rstar

end subroutine logg_calc2

subroutine logg_calc
  use params
  implicit none

  double precision :: rstar, mstar
  write(*,'(A)',advance='no') 'Input radius of star (solar units)'
  read(*,*) rstar
  rstar = rstar*Rsun*100
  write(*, '(A)', advance='no') 'Input mass of star (solar units): '
  read(*,*) mstar
  mstar = mstar*MSun*1000

  print*, 'log g (cgs) = ', log10(G*1.0e3*mstar/rstar/rstar)

end subroutine logg_calc
