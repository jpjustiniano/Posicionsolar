Program jdcalc

 y =2010
 m = 10
 d = 15
 hour = 20
 jd =(1461*(y+4800+(m-14)/12))/4+(367*(m-2-12*((m-14)/12)))/12-(3*((y+4900+(m-14)/12)/100))/4+d-32075-.5+hour/24.0
JD2= d-32075+1461*(y+4800+(m-14)/12)/4+367*(m-2-(m-14)/12*12)/12-3*((y+4900+(m-14)/12)/100)/4
write(*,*) jd, jd2
read (*,*) 

end Program