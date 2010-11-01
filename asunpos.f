c   this program calculates the solar position given the year, day,
c   time, latitude, and longitude
c
c   it outputs solar azimuth, elevation, hour angle, declination, and
c   air mass
c
c
	program sunpos
	real lat,long, mn
	write(6,'(3x,"This program calculates sun position if you
     # specify time and location!")')
	write(6,'(/)')
	write(6,'(3x,"Type the latitude like this sxx.xx where s is a
     # sign!")')
	read(5,'(f8.4)')lat
	write(6,'(f8.4)')lat
	write(6,'(3x,"Type the longitude like this sxxx.xx!")')
	read(5,'(f9.4)')long
	write(6,'(f9.4)')long
	write(6,'(3x,"Type the year like this xxxx.!")')
	read(5,'(f5.0)')year
	write(6,'(f5.0)')year
	write(6,'(3x,"Indicate time zone,e.g., est=+5,mst=+7,etc.(The
     # sign is important)!")')
	read(5,'(i3)')itz
	write(6,'(i3)')itz
	tz=float(itz)
60      write(6,'(3x,"Type the day of year(Feb 1=32) like this xxx.!")')
	read(5,'(f4.0)')day
	write(6,'(f4.0)')day
50      write(6,'(3x,"Type the local standard time(not daylight savings
     # time) like this xx.xx.xx., e.g., 12.15.47.!")')
	read(5,'(3f3.0)')hour, mn,sec
	write(6,'(3f3.0)')hour, mn,sec
	tm=((tz+hour)*3600+ mn*60+sec)/3600.
	call sunae(year,day,tm,lat,long,az,el,ha,dec,soldst)
	write(6,'(//)')
	write(6,'("The solar azimuth and elevation are")')
	write(6,100)az,el
100     format(3x,f9.4,3x,f8.4,/)
	write(6,'("The solar hour angle and declination are")')
	write(6,100)ha,dec
	am=armass(el)
	write(6,'("The air mass is")')
	write(6,200)am
200     format(3x,f6.2)
      write(6,'("The solar distance is")')
      write(6,250)soldst
250     format(3x,f7.4)
	write(6,'(//)')
	write(6,'("If you want another time same day, type 9!")')
	read(5,'(i1)')i
	if(i.eq.9)go to 50
	write(6,'("If you want another day, type 99!")')
	read(5,'(i2)')i
	if(i.eq.99)go to 60
	end
c---------------------------------------------------------------------+
c                                                                     |
c   subroutine sunae(year,day,hour,lat,long,az,el,ha,dec,soldst)      |
c                                                                     |
c   this subroutine calculates the local azimuth and elevation of the |
c     sun at a specific location and time using an approximation to   |
c     equations used to generate tables in The Astronomical Almanac.  |
c     refraction correction is added so sun position is apparent one. |
c                                                                     |
c   The Astronomical Almanac, U.S. Gov't Printing Office, Washington, |
c     D.C. (1985).                                                    |
c                                                                     |
c   input parameters                                                  |
c     year=year, e.g., 1986                                           |
c     day=day of year, e.g., feb 1=32                                 |
c     hour=hours plus fraction in UT, e.g., 8:30 am eastern daylight  |
c       time is equal to 8.5 + 5(5 hours west of Greenwich) -1(for    |
c       daylight savings time correction)                             |
c     lat=latitude in degrees (north is positive)                     |
c     long=longitude in degrees (east is positive)                    |
c                                                                     |
c   output parameters                                                 |
c     a=sun azimuth angle (measured east from north, 0 to 360 degs)   |
c     e=sun elevation angle (degs)                                    |
c     plus others, but note the units indicated before return         |
c                                                                     |
c---------------------------------------------------------------------+
c
      subroutine sunae(year,day,hour,lat,long,az,el,ha,dec,soldst)
c  work with real variables and define some constants, including
c  one to change between degs and radians
      implicit real (a-z)
      data twopi,pi,rad/6.2831853,3.1415927,.017453293/
c
c   get the current julian date (actually add 2,400,000 for jd)
      delta=year-1949.
      leap=aint(delta/4.)
      jd=32916.5+delta*365.+leap+day+hour/24.
c   1st no. is mid. 0 jan 1949 minus 2.4e6; leap=leap days since 1949
c  the last yr of century is not leap yr unless divisible by 400
      if(amod(year,100.).eq.0.0.and.amod(year,400.).ne.0.0)jd=jd-1.
c
c   calculate ecliptic coordinates
      time=jd-51545.0
c   51545.0 + 2.4e6 = noon 1 jan 2000
c
c   force mean longitude between 0 and 360 degs
      mnlong=280.460+.9856474*time
      mnlong=mod(mnlong,360.)
      if(mnlong.lt.0.)mnlong=mnlong+360.
c
c   mean anomaly in radians between 0 and 2*pi
      mnanom=357.528+.9856003*time
      mnanom=mod(mnanom,360.)
      if(mnanom.lt.0.)mnanom=mnanom+360.
      mnanom=mnanom*rad
c
c   compute the ecliptic longitude and obliquity of ecliptic in radians
      eclong=mnlong+1.915*sin(mnanom)+.020*sin(2.*mnanom)
      eclong=mod(eclong,360.)
      if(eclong.lt.0.)eclong=eclong+360.
      oblqec=23.439-.0000004*time
      eclong=eclong*rad
      oblqec=oblqec*rad
c
c   calculate right ascension and declination
      num=cos(oblqec)*sin(eclong)
      den=cos(eclong)
      ra=atan(num/den)
c   force ra between 0 and 2*pi
      if(den.lt.0)then
          ra=ra+pi
      elseif(num.lt.0)then
          ra=ra+twopi
      endif
c
c   dec in radians
      dec=asin(sin(oblqec)*sin(eclong))
c
c   calculate Greenwich mean sidereal time in hours
      gmst=6.697375+.0657098242*time+hour 
c   hour not changed to sidereal time since 'time' includes
c   the fractional day 
      gmst=mod(gmst,24.)
      if(gmst.lt.0.)gmst=gmst+24.
c
c   calculate local mean sidereal time in radians 
      lmst=gmst+long/15.
      lmst=mod(lmst,24.)
      if(lmst.lt.0.)lmst=lmst+24.
      lmst=lmst*15.*rad
c
c   calculate hour angle in radians between -pi and pi
      ha=lmst-ra
      if(ha.lt.-pi)ha=ha+twopi
      if(ha.gt.pi)ha=ha-twopi
c
c   change latitude to radians
      lat=lat*rad
c
c   calculate azimuth and elevation
      el=asin(sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(ha))
      az=asin(-cos(dec)*sin(ha)/cos(el))
c
c   this puts azimuth between 0 and 2*pi radians
      if(sin(dec)-sin(el)*sin(lat).ge.0.)then
      if(sin(az).lt.0.)az=az+twopi
      else
      az=pi-az
      endif
cc   if az=90 degs, elcritical=asin(sin(dec)/sin(lat))
cc    elc=asin(sin(dec)/sin(lat))
cc    if(el.ge.elc)az=pi-az
cc    if(el.le.elc.and.ha.gt.0.)az=twopi+az
c
c   calculate refraction correction for US stan. atmosphere
c   need to have el in degs before calculating correction
      el=el/rad
c
      if(el.ge.19.225) then 
         refrac=.00452*3.51823/tan(el*rad)
      else if (el.gt.-.766.and.el.lt.19.225) then
         refrac=3.51823*(.1594+.0196*el+.00002*el**2)/
     1   (1.+.505*el+.0845*el**2)
      else if (el.le.-.766) then
         refrac=0.0
      end if
c
c   note that 3.51823=1013.25 mb/288 C
      el=el+refrac
c   elevation in degs
c
c   calculate distance to sun in A.U. & diameter in degs
      soldst=1.00014-.01671*cos(mnanom)-.00014*cos(2.*mnanom)
      soldia=.5332/soldst
c
c   convert az and lat to degs before returning
      az=az/rad
      lat=lat/rad
	 ha=ha/rad
	 dec=dec/rad
c   mnlong in degs, gmst in hours, jd in days if 2.4e6 added;
c   mnanom,eclong,oblqec,ra,and lmst in radians
      return
      end
c
c
c   this function calculates air mass using kasten's
c   approximation to bemporad's tables
c
c
	function armass(el)
	z=(90.-el)*3.141592654/180.
	armass=1./(cos(z)+.50572*(6.07995+el)**-1.6364)
	return
	end
