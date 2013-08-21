      program nauenberg

c     given Rw/a = y, q and Period, calculates Mw and therefore Rw
c     uses the Nauenberg 1972 relation taken from Cook & Warner 1984
c     this corresponds to Hamada-Salpeter relation (approx.)

      implicit none

      double precision grav, pi, msun, rsun, y, q, period, part1, part2
      double precision mw43, mw, rw, m2, a, part3, Kw, Kr, r2, inc
      double precision yerr, qerr, perr, incerr, mwerr, err, value
      double precision rwerr, r2err, Kwerr, Krerr, aerr, m2err, errsini
      double precision constant,correction
      integer i

c     input constants in SI units

      grav = 6.67e-11
      pi = 4.*atan(1.)
      msun = 1.9891e30
      rsun = 6.96e8

      write(*,'(a,$)') 'Enter Rw/a = y and error : '
      read(*,*) y, yerr
      write(*,'(a,$)') 'Enter q and error : '
      read(*,*) q, qerr
      write(*,'(a,$)') 'Enter period (days) and error : '
      read(*,*) period, perr
      write(*,'(a,$)') 'Input inclination (degrees) and error : '
      read(*,*) inc, incerr
      write(*,'(a,$)') 'Input correction factor in percent : '
      read(*,*) correction

      correction = 1+(correction/100.)

      period = period*24*60*60

      constant = 7.795e6

      do i=1,2

         if (i .eq. 1) then
            write(*,*) 'Nauenberg cold'
            write(*,*) ' '
         else if (i .eq. 2) then
            write(*,*) 'Nauenberg hot'
            write(*,*) ' '
         end if
          

         part1 = ((grav*(1+q)*(period**2.))/(4*(pi**2.)))**(2./3.)
         part2 = (y/CONSTANT)**2.
         part3 = 1.44*msun
         
c     calculate primary mass...
         mw43 = (part3**(2./3.))/((part2*part1)+(part3**(-2./3.)))
         mw = mw43**(3./4.)
         mw = mw/msun
         
c     and error
         value = part2*part1
         err = value*sqrt((2.*yerr/y)**2.+(2./3.*(qerr)/(1.+q))**2.)
         
         mwerr=mw*3./4.*((err)/(value+(part3**(-2./3.))))
         
         write(*,*) ' '
         write(*,*) 'Primary mass        = ',mw,' +/-',mwerr,' Msun'
         
c     calculate primary radius...
         rw = CONSTANT*sqrt(((1.44/mw)**(2./3.))-((mw/1.44)**(2./3.)))
         rw = rw/rsun
         
c     and error
         rwerr=rw*sqrt((yerr/y)**2.+(mwerr/3./mw)**2.+
     $        (qerr/3./(1+q))**2.)
         
         write(*,*) 'Primary radius      = ',rw,' +/-',rwerr,' Rsun'
         
c     check this is consistent...
         rw = y*(grav*mw*msun*(1+q)*(period**2.)/4/(pi**2.))**(1./3.)
         rw = rw/rsun
         
         write(*,*) 'Primary radius (2)  = ',rw,' Rsun'
         
c     so get m2 from q
         m2 = q*mw
         
c     and error
         m2err = m2*sqrt((qerr/q)**2.+(mwerr/mw)**2.)
         
         write(*,*) 'Secondary mass      = ',m2,' +/-',m2err,' Msun'
         
c     orbital separation...
         a=((grav*(1.+q)*mw*msun*period**2./(4.*pi**2.))**(1./3.))/rsun
         
c     and error
         aerr = a/3.*sqrt((qerr/(1.+q))**2.+(2.*perr/period)**2.+
     $        (mwerr/mw)**2.)
         write(*,*) 'Orbital separation  = ',a,' +/-',aerr,' Rsun'
         
c     now get radial velocities
         Kw = sin(inc*pi/180.)*2*pi*q*a*rsun/(1+q)/period/1000.
         Kr = sin(inc*pi/180.)*2*pi*a*rsun/(1+q)/period/1000.
         
c     get errors
         errsini = 0.5*abs(sin((inc+incerr)*pi/180.)-
     $        sin((inc-incerr)*pi/180.))
         Kwerr = Kw*sqrt((errsini/sin(inc*pi/180.))**2.+(qerr/q)**2.
     $        +(aerr/a)**2.+(qerr/(1+q))**2.+(perr/period)**2.)
         Krerr = Kw*sqrt((errsini/sin(inc*pi/180.))**2.
     $        +(aerr/a)**2.+(qerr/(1+q))**2.+(perr/period)**2.)
         
         
         write(*,*) 'Kw                  = ',Kw,' +/-',Kwerr,' km/s'
         write(*,*) 'Kr                  = ',Kr,' +/-',Krerr,' km/s'
         
c     companion radius
c     use from Warner 1995 (Eggleton etc)
c     radius of sphere with same radius as Roche lobe
         r2 =a*rsun*0.49*q**(2./3.)/(0.6*q**(2./3.)+log(1+q**(1./3.)))
     $        /rsun
         
c     get error
         r2err=r2*sqrt((aerr/a)**2.+(2./3.*qerr/q)**2.+
     $        ((2./3.*qerr/(q**(1./3.)))**2.+
     $        (qerr/(3.*(q**(2./3.))*(1.+q)))**2.)/(0.6*(q**(2./3.))+
     $        log(1.+(q**(1./3.))))**2.)
         
         write(*,*) 'Radius of secondary = ',r2,' +/-',r2err,' Rsun'
         write(*,*) ' '
         
         constant = 7.795e6*correction
      
      end do

      end
