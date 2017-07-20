	PROGRAM maincontrol
	IMPLICIT NONE

! This program moves one of the hydrogen atoms of a formaldehyde molecule 
!from equilibrium around the molecule on the oxygen side to the other 
!hydrogen.

	REAL(8), DIMENSION(36,3) :: xx
	REAL(8), DIMENSION(36) :: xmass
	REAL(8) :: v,v1,v2,s,r,theta,pi,c,a,b,e,l,r2,phi,psi,rxy,del,lim
	INTEGER :: natom,k,j,pts

	pi = 4*ATAN(1d0)
!	PRINT *,pi
	natom = 4

	xmass(1) = 1837.15  ! H
	xmass(2) = 29156.95   ! O
	xmass(3) = 21874.66   ! C
	xmass(4) = 1837.15  ! H 
C  *** Formaldehyde eqiuilibrium geometry ***
c	xx(1,1) = -1.12460784  
c	xx(1,2) = 1.79389857
c	xx(1,3) = 0.00000000
c	xx(2,1) = 2.28820140
c	xx(2,2) = 0.00000000
c	xx(2,3) = 0.00000000
c	xx(3,1) = 0.00000000
c	xx(3,2) = 0.d0
c	xx(3,3) = 0.d0
c	xx(4,1) = -1.12460784
c	xx(4,2) = -1.79389857
c	xx(4,3) = 0.00000000
cccccccccccccccccccccccccccccccccccccccccc

C  *** HCO equilibrium geometry ***
	xx(1,1) = -1.27500000  
	xx(1,2) = 1.666000000
	xx(1,3) = 0.00000000
	xx(2,1) = 2.21320000
	xx(2,2) = 0.00000000
	xx(2,3) = 0.00000000
	xx(3,1) = 0.00000000
	xx(3,2) = 0.d0
	xx(3,3) = 0.d0

	CALL getpot(natom,xmass,xx,v,v1,v2,s)
	PRINT *, v*219475d0

	OPEN(UNIT=25,FILE='data.txt',STATUS='REPLACE')
	r = SQRT(xx(4,1)*xx(4,1)+xx(4,2)*xx(4,2))
	theta = ATAN(xx(4,2)/xx(4,1))
	PRINT *, theta*180d0/pi
C elliptical orbit components--------------
C	c = 2.28820140d0/2d0
C	a = 0.5*(r+SQRT(r**2+(c/2d0)**2+2d0*r*c/2d0*COS(theta)))
C	b = SQRT(a**2-c**2)
C	e = SQRT(1d0-b**2/a**2)
C-------------------------------------------
!	theta = theta+60d0*(pi/180d0)
	l = r
	PRINT *,l
	phi = ACOS(xx(4,2)/xx(4,1))
	pts = 20
	lim = 7d0
	r = lim
c	r = r+l
	del = 30d0*pi/180d0
c	xx(4,3) = r*SIN(del)
c	rxy = r*COS(del)
	xx(4,1) = -r*COS(theta+del)
	xx(4,2) = r*SIN(theta+del)
	DO k=1,pts
	  phi = theta-45d0*pi/180d0+REAL(k-1)*90d0/REAL(pts-1)*pi/180d0
c	  xx(1,1) = -l*COS(phi)
c	  xx(1,2) = l*SIN(phi)
c	  xx(4,1) = -(l+r)*COS(theta)
c	  xx(4,2) = (l+r)*SIN(theta)
	  DO j = 1,pts
	    r = REAL(j-1)*(10.5d0)/REAL(pts-1)+1.5d0
c	    psi = pi/180*(-90d0+180/REAL(pts-1)*REAL(j-1))
c	    xx(1,3) = r*SIN(psi)
c	    rxy = r*COS(psi)
	    xx(1,1) = -r*COS(phi)
	    xx(1,2) = r*SIN(phi)
	    CALL Getpot(natom,xmass,xx,v,v1,v2,s)
	    v = v*219474.63067d0
	    WRITE(UNIT=25,FMT=*) REAL(r),
     $        REAL(180-(phi)*180d0/pi),REAL(v)
	  END DO
	END DO

	CLOSE(UNIT=25)

	END PROGRAM
