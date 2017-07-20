C**********************************************************************
C this getpot() subroutine use cartesian coordinates as input (in bohr)
C it will return potential value in wavenumber, combining CCSD(T) fit with MRCI
C fit with switch range 1.5 - 2.0
C**********************************************************************

      subroutine Getpot (natom,xmass,xx,v,v1,v2,s)

      implicit double precision (a-h,o-z)
      dimension xx(36,3),rr(36,36),r(6),rmrci(6)
      dimension xmass(36)
      data rech,reoh,rehh / 2.15d0 , 1.85d0, 1.45d0 /

      hcm=219474.63

c Switching range

      rl=1.5
      ru=2.0

      xmab=(xmass(1)*xx(1,1)+xmass(4)*xx(4,1))/(xmass(1)+xmass(4))
      ymab=(xmass(1)*xx(1,2)+xmass(4)*xx(4,2))/(xmass(1)+xmass(4))
      zmab=(xmass(1)*xx(1,3)+xmass(4)*xx(4,3))/(xmass(1)+xmass(4))
c     write(*,*) xmab,ymab,zmab

      xmcd=(xmass(2)*xx(2,1)+xmass(3)*xx(3,1))/(xmass(2)+xmass(3))
      ymcd=(xmass(2)*xx(2,2)+xmass(3)*xx(3,2))/(xmass(2)+xmass(3))
      zmcd=(xmass(2)*xx(2,3)+xmass(3)*xx(3,3))/(xmass(2)+xmass(3))
c     write(*,*) xmcd,ymcd,zmcd

      rbig=sqrt((xmcd-xmab)*(xmcd-xmab)+
     $          (ymcd-ymab)*(ymcd-ymab)+
     $          (zmcd-zmab)*(zmcd-zmab))
c     write(*,*) rbig

c     do i=1,4
c     write(*,*) (xx(i,j),j=1,3)
c     end do

c Calculating the internuclear distances

      Call BONDS (natom,xx,rr) 

      r(1)=rr(1,2)
      r(2)=rr(3,4)
      r(3)=rr(2,3)
      r(4)=rr(1,3)
      r(5)=rr(2,4)
      r(6)=rr(1,4)

      rmrci(1)=rr(2,3)
      rmrci(2)=rr(1,3)
      rmrci(3)=rr(3,4)
      rmrci(4)=rr(1,2)
      rmrci(5)=rr(2,4)
      rmrci(6)=rr(1,4)

c     do i=1,6
c      write(*,*) r(i)
c     end do

c calculating the switch function 

      scha = rr(1,3)/rech
      schb = rr(3,4)/rech
      soha = rr(1,2)/reoh
      sohb = rr(2,4)/reoh
      shh  = rr(1,4)/rehh

      Ra=dmin1(scha,soha,shh)
      Rb=dmin1(schb,sohb,shh)
      S0=dmax1(Ra,Rb)

      x=(S0-rl)/(ru-rl)

      if (x.LT.0.) then
        s=1.0
      end if
      if (x.GT.1.) then
        s=0.0
      end if
      if ((x.GE.0.).and.(x.LE.1.)) then
        s=1.0-(10.0*x**3-15.0*x**4+6.0*x**5)
      end if

c Getting potential from MRCI fit

        Call h2co_mrci(rmrci,en) 
        v1=(en+114.341634)

c Getting potential from CCSC(T) fit

        Call h2co_ccsdt(r,rbig,v2)
        v2=v2/hcm

c Combining the two fit with switch function s

      v=v2*s+v1*(1-s)

c     write(*,*) v1,v2,v
      return
      end

C*************************************************************************
C*************************************************************************
