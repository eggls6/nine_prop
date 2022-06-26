module getephem_m
use global_m
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!		reads planetary ephemeris from JPL binary ephemeris files
!		not usable for DE406 or versions<DE405
! 
! 		written by Siegfried Eggl 20150110
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 public::getdeversion
 public::getrv
 public::getrvbody
 public::getrvhel
 public::getrvsun
 public::getp430
 public::getp405
 public::getparam
 
 interface 
  subroutine Pleph(et,i,nctr,rvread)
  use global_m
  integer,intent(in)::i,nctr
  double precision,intent(out)::rvread(1:6)
  double precision,intent(in)::et
  end subroutine
 end interface
 
 
 contains 
!############################################################################# 
subroutine getrv(et,rv) 
!reads rv from JPL binary ephemeris, needs Pleph subroutine from readephem
!input: 
!et... time in JD
!np... number of massive bodies including the sun 
!rv... (1:np,1:6) distance and velocities to the solar system barycenter
!nctr... center id (11 Sun, 12 SolSys barycenter (does not coincide with calculated barycenter due to missing asteroids)
implicit none
integer,parameter::nctr=12,np=11
integer::i
real(kind=dp),intent(in)::et
real(kind=dp),intent(out)::rv(1:np,1:6)
double precision:: rvread(1:np,1:6),etin

etin=dble(et)

 do i=1,np
        CALL  Pleph ( etin, i, NCTR, rvread(i,1:6)) 
 end do

do i=1,np
 rv(i,1:6)=real(rvread(ord(i),1:6),kind=dp)
end do


return
end subroutine
!********************************************************************************************************
subroutine getrvhel(et,rv) 
!reads rv from JPL binary ephemeris, needs Pleph subroutine from readephem
!input: 
!et... time in JD
!np... number of massive bodies including the sun 
!rv... (1:np,1:6) distance and velocities to the SUN (HELIOCENTRIC DISTANCES AND VELOCITIES
!nctr... center id (11 Sun, 12 SolSys barycenter (does not coincide with calculated barycenter due to missing asteroids)
implicit none
integer,parameter::nctr=11,np=11
integer::i
real(kind=dp),intent(in)::et
real(kind=dp),intent(out)::rv(1:np,1:6)
double precision:: rvread(1:np,1:6),etin

etin=dble(et)

 do i=1,np
        CALL  Pleph (etin, i, NCTR, rvread(i,1:6))
 end do

do i=1,np
          rv(i,1:6)=real(rvread(ord(i),1:6),kind=dp)
end do


return
end subroutine
!###########################################################################
subroutine getrvsun(et,rv1) 
!reads rv from JPL binary ephemeris, needs Pleph subroutine from readephem
!input: 
!et... time in JD 
!rv1... (1:6) distance and velocities to the solar system barycenter of the Sun
implicit none
integer,parameter::nctr=12,np=11
real(kind=dp),intent(in)::et
real(kind=dp),intent(out)::rv1(1:6)
double precision::rvread(1:6),etin

        etin=dble(et)
!the Sun is body nr. 11
        CALL  Pleph ( etin , 11, NCTR, rvread(1:6)) 
        
        rv1(1:6)=real(rvread(1:6),kind=dp)
        

return
end subroutine
!###########################################################################
subroutine getrvbody(et,body,rv1) 
!reads rv from JPL binary ephemeris, needs Pleph subroutine from readephem
!input: 
!et... time in JD
!body ... JPLID of body to get ephemeris for
!rv1... (1:6) distance and velocities to the solar system barycenter of the Earth
!
!                1 = MERCURY           8 = NEPTUNE
!                2 = VENUS             9 = PLUTO
!                3 = EARTH            10 = MOON
!                4 = MARS             11 = SUN
!                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
!                6 = SATURN           13 = EARTH-MOON BARYCENTER
!                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)

implicit none
integer,parameter::nctr=12,np=11
integer,intent(in)::body
real(kind=dp),intent(in)::et
real(kind=dp),intent(out)::rv1(1:6)
double precision::rvread(1:6),etin

        etin=dble(et)

        CALL  Pleph ( etin, body, NCTR,rvread(1:6)) 
      
        rv1(1:6)=real(rvread(1:6),kind=dp)

return
end subroutine
!###############################################################################
subroutine getp430
!get start and end times of the ephemeris files and masses
implicit none
integer,parameter::np=11
integer::i,nvs
double precision::vals(1:400),TJD(1:3),m(1:np),emrat
character*6::nams(1:400)

call const(nams,vals,TJD,nvs)

 tephstart=real(TJD(1),kind=dp)
 tephend=real(TJD(2),kind=dp)
 aukm=real(vals(10),kind=dp)
 clight=real(vals(7),kind=dp)/aukm*86400._dp 
 cm2=clight**(-2._dp)
 kgc2=real(vals(21),kind=dp)
 kgc=sqrt(kgc2)
 emrat=real(vals(11),kind=dp)


!order ascending in distance from the sun
m(1)=vals(12)                       !Mercury
m(2)=vals(13)                      !Venus

m(4)=vals(15)                      !Mars
m(5)=vals(16)                      !Jupiter
m(6)=vals(17)                      !Saturn
m(7)=vals(18)                      !Uranus
m(8)=vals(19)                      !Neptue
m(9)=vals(20)                      !Pluto
m(10)=vals(14)/(emrat+1.d0)        !moon
m(11)=1.d0                             !Sun

m(3)=vals(14)-m(10)                 !Earth


!order ascending in mass
do i=1,10
   gmeph(i)=real(m(ord(i)),kind=dp)
   meph(i)=gmeph(i)/kgc2
end do
   gmeph(11)=kgc2
   meph(11)=1._dp


return
end subroutine
!***********************************************************
subroutine getp405 
!get start and end times of the ephemeris files and masses
implicit none
integer,parameter::np=11
integer::i,nvs
double precision::vals(1:400),TJD(1:3),m(1:np),emrat
character*6::nams(1:400)

call const(nams,vals,TJD,nvs)

 tephstart=real(TJD(1),kind=dp)
 tephend=real(TJD(2),kind=dp)
 aukm=real(vals(7),kind=dp)
 clight=real(vals(6),kind=dp)/aukm*86400._dp 
 cm2=clight**(-2._dp)
 kgc2=real(vals(18),kind=dp)
 kgc=sqrt(kgc2)
 emrat=real(vals(8),kind=dp)

 
!order ascending in distance from the sun
m(1)=vals(9)                       !Mercury
m(2)=vals(10)                      !Venus

m(4)=vals(12)                      !Mars
m(5)=vals(13)                      !Jupiter
m(6)=vals(14)                      !Saturn
m(7)=vals(15)                      !Uranus
m(8)=vals(16)                      !Neptue
m(9)=vals(17)                      !Pluto
m(10)=vals(11)/(emrat+1.d0)        !moon
m(11)=1.d0                             !Sun

m(3)=vals(11)-m(10)                 !Earth

!order ascending in mass
do i=1,10
   gmeph(i)=real(m(ord(i)),kind=dp)
   meph(i)=gmeph(i)/kgc2
end do
   gmeph(11)=kgc2
   meph(11)=1._dp
   

return
end subroutine

!*******************************************************************************
subroutine getdeversion
!Get JPL DExxx version
!This will determine the readout routines since headers changed from DE423 to DE430
implicit none
integer,parameter::np=11
integer::nvs
real*8::vals(1:400),TJD(1:3)
character*6::nams(1:400)

call const(nams,vals,TJD,nvs)

write(*,*)'using ephermeris ',nams(1),vals(1)

DE4XX=vals(1)

return
end subroutine
!############################################################################
end module


