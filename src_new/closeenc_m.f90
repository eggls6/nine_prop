 module closeenc_m
  use global_m
  use transform_m
  use output_m
 

!  interface 
!   real(kind=dp) function q(a,b,c,d,t,soi)
!   use global_m
!   real(kind=dp),intent(in)::a,b,c,d,t,soi
!   end function
! 
!   real(kind=dp) function qdot(b,c,d,t)
!    use global_m
!   real(kind=dp),intent(in)::b,c,d,t
!   end function
!  end interface
 
  public::cccheck
  public::tpcoo
  public::q
  public::qdot

  contains
!*************************************************************************************  
real(kind=dp) function q(a,b,c,d,t,soi)
!calculates cubic spline interpolation value - distance to be reached (soi) for Newton Raphson, 
!a,b,c,d... spline coefficients (real)
!t......... spline parameter/ time (real)
!soi......  target distance (real)
use global_m
implicit none
real(kind=dp),intent(in)::a,b,c,d,t,soi

q=a+b*t+c*t*t+d*t**3.d0-soi

return
end function
!***********************
real(kind=dp) function qdot(b,c,d,t)
!calculates fist derivative of cubic spline interpolation value  for Newton Raphson
!b,c,d...  dspline coefficients (real)
!t......... dspline parameter/ time (real)
use global_m
implicit none
real(kind=dp),intent(in)::b,c,d,t

qdot=b+2.d0*c*t+3.d0*d*t*t

return
end function
!*********************************************************************************
  subroutine cccheck(t,dt,rv,id)
  use global_m
  !this subroutine checks for close encounters of the non ephemeris bodies with the ephemeris body "id"
  implicit none
  integer(kind=idp)::i,ipm0
  integer,intent(in)::id
  real(kind=dp),intent(in)::t,dt
  real(kind=dp),intent(in)::rv(:,:)
  real(kind=dp),dimension(1:6)::rveph,rvepho,rvephsoi,rvint,rvag,rvago
  real(kind=dp)::b,tp(1:2),dm,tpp,n
  real(kind=dp)::dag,ddag,dago,ddago,tmin,rmin,vmin,tsoi,soi
  
  

  !get ephemeris for now
  call getrvbody(t,ord(id),rveph)
  !get ephemeris for previous timestep
  call getrvbody(t-dt,ord(id),rvepho)
  !output to file when body comes within cclim
  do i=1,inbody
  ipm0=i+m0cur
  !calculate geocentric scalar distance (dag) and velocity (ddag) for now
  rvag=rv(i,1:6)
  call bodykoo1(rvag,ord(id),t)
  dag=sqrt(dot_product(rvag(1:3),rvag(1:3)))
  ddag=sqrt(dot_product(rvag(4:6),rvag(4:6)))
  
   !calculate geocentric scalar distance (dag) and velocity (ddag) for the previous step
  rvago=rvo(i,1:6)
  call bodykoo1(rvago,ord(id),t-dt)
  dago=sqrt(dot_product(rvago(1:3),rvago(1:3)))
  ddago=sqrt(dot_product(rvago(4:6),rvago(4:6)))
  
  !do spline interpolation to see whether minimum distance during the previous timestep dt is smaller than cclim
  call  srmin(dago,dag,ddago,ddag,dt,tmin,rmin,vmin)
    
!       write(*,*)'rv'
!     write(*,*)rv(i,1:6)
!           write(*,*)'rvo'
!     write(*,*)rvo(i,1:6)
!     write(*,*)'delta rv'
!     write(*,*)rv(i,1:6)-rvo(i,1:6)
!     write(*,*)'delta rveph'
!     write(*,*)rveph-rvepho
!   write(*,*)'dt',dt  
!   write(*,*)'rvag'
!   write(*,*)rvag
!   write(*,*)'rvago'
!   write(*,*)rvago
!   write(*,*)'dag,dago,ddag,ddago'
!   write(*,*)dag,dago,ddag,ddago  
!     
!     
  

  
  !check whether asteroid is within close encounter limits
  if(rmin.le.cclim) then
  !construct b-plane at nsoi x sphere of influence 
  
     soi=meph(id)**(2._dp/5._dp)*nsoi
     
          
     if(ccenter(i)) then
     else
     ncc(ipm0)=ncc(ipm0)+1
     end if
     
     ccenter(i)=.true.

    !db
   !    write(*,*)'dt',dt  
   !   write(*,*)'rvag'
   !   write(*,*)rvag
!   write(*,*)'rvago'
!   write(*,*)rvago
!   write(*,*)'dag,dago,ddag,ddago,rmin,tmin'
!   write(*,*)dag,dago,ddag,ddago,rmin,tmin  
!     write(*,*)'dag<cclim',dag,dagold(i),soi
!   edb


!catalogue minimum distances as global minima?
  !only start cataloguing if t greater than catalogue limit cct0
  if (t.gt.cct0) then
  !db
  !write(*,*)'I am catalogueing'
  
  
    !calculated minimum distance (rmin) via spline interpolation
    if(abs(rmin).lt.dagminn(ipm0)) then
      dagminn(ipm0)=rmin
      tdagminn(ipm0)=t+tmin
      vcc(ipm0)=ddag
    end if  

    !calculate minimum distance (dm) via hyperbolic fit
     call dmin(id,i,t,dag,rvag(1:6),dm,tpp)
  
     if(dm.lt.dagmina(ipm0)) then
      dagmina(ipm0)=dm
      tdagmina(ipm0)=tpp
      vinf(ipm0)=sqrt(dot_product(rvag(4:6),rvag(4:6)))
     end if            
    
!calculate target plane coordinates when entering SOI
    !db
    !write(*,*)'dago',dago,'dag',dag,'soi',soi
    !if(dago.gt.soi.and.dag.le.soi) then 
    !projection of the relative distance between the asteroid and the Earth onto the asteroid's relative velocity vector (which is orthogonal to the b-plane)
    !gives the distance of the asteroid to the b-plane along the asymptote. if this distance is smaller than soi, then...
   
   !db
   !write(*,*)'dago',dago,'dag',dag,'soi',soi
   !write(*,*)'rvag/dag',Dot_product(rvag(4:6)/ddag,rvag(4:6)/ddag),'rvago/dago',Dot_product(rvago(4:6)/ddago,rvago(4:6)/ddago)
   !write(*,*)'proj',abs(dot_product(rvag(1:3),rvag(4:6)/ddag)),'projo',abs(dot_product(rvago(1:3),rvago(4:6)/ddago)),'soi',soi
   if(abs(dot_product(rvag(1:3),rvag(4:6)/ddag)).le.soi.and.abs(dot_product(rvago(1:3),rvago(4:6)/ddago)).gt.soi) then
    
    !sphere of influence entered, calculate time of entrance (tsoi)
     call srsoi(dago,dag,ddago,ddag,dt,soi,tsoi)
        ndcc(ipm0)=ndcc(ipm0)+1
    
     !interpolate asteroid position and velocity at time of entrance into soi
     !to get b-plane coordinates
     call srvint(rvo(i,1:6),rv(i,1:6),dt,tsoi,rvint(1:6))
  
       !get ephemeris at interpolated time
     call getrvbody(t+tsoi,ord(id),rvephsoi)
  
   
     call tpcoo(t+tsoi,rvint(1:6),i,rvephsoi(1:6),id,b,tp(1:2))
     call ccout(ipm0,id,t+tsoi,soi,tpp,dm,b,tp)
     
    ! write(*,*)'fbplane(i)=',fbplane(i)
     
     if(fbplane(i)) then
         !first b-plane
         fbpcoo(ipm0,1:2)=tp(1:2)*aukm
         fbpcoo(ipm0,3)=t+tsoi
         
         !local close encounter b-plane
         loccc(ipm0,ndcc(ipm0),4:5)=tp(1:2)*aukm
         loccc(ipm0,ndcc(ipm0),6)=t+tsoi
         
         fbplane(ipm0)=.false.
     end if
     
     
          
     !local close encounter minima via spline interpolation
    if(abs(rmin).lt.loccc(ipm0,ncc(ipm0),2)) then
     loccc(ipm0,ndcc(ipm0),2)=rmin
     loccc(ipm0,ndcc(ipm0),1)=t+tmin
     loccc(ipm0,ndcc(ipm0),3)=ddag
    end if
     
     
!   write(*,*)'dt,tsoi,t',dt,tsoi,t  
!   write(*,*)'rmin,dm'
!   write(*,*)rmin,dm
!   write(*,*)'rvepho'
!   write(*,*)rvepho
!   write(*,*)'rvephsoi'
!   write(*,*)rvephsoi
!   write(*,*)'rveph'
!   write(*,*)rveph
!   write(*,*)'dag,dago,ddag,ddago,rmin,tmin'
!   write(*,*)dag,dago,ddag,ddago,rmin,tmin  

    end if
    
   end if !catalogue
    
  else
     ccenter(i)=.false.
  end if
   
 
 !debug
!  if (i.eq.4) then
!   write(*,*)'i,dag,aukm',i,dag,aukm
!  end if 
 !end debug
  
  
!    if(dag.gt.cclim) then
!    ccenter(i)=.true.
!    end if
  end do
  

     
     
   !   !start saving previous stuff
!   nsavep1=nsave+1
!   !latest value is rvold(:,:,1)
!   do i=1,nsave-1
!   rvold(:,:,nsave-i)=rvold(:,:,nsavep1-i)
!   dold(:,nsave-i)=dold(:,nsavep1-i)
!   end do
   return
   end subroutine
!***********************************************************************************
   subroutine dmin(id,i,t,dag,rvag,dm,tpp)
   implicit none
   !input 
   integer,intent(in)::id
   integer(kind=idp),intent(in)::i
   real(kind=dp),intent(in)::t,dag
   real(kind=dp),dimension(1:6),intent(in)::rvag
   !output minimum geocentric distance (pericenter for hyperbolic and parabolic orbits)
   real(kind=dp),intent(out)::dm,tpp
   
   !local
   real(kind=dp)::mu,a,e,p,v2,h(1:3),vesc2,egy
   real(kind=dp)::sinhf,coshf,f
   
   mu=kgc2*(meph(id)+mass(i))
   v2=dot_product(rvag(4:6),rvag(4:6))
   
   !calculate specific orbital energy sign to determine whether the orbit is 
   !elliptic, parabolic or hyperbolic
   
   eps=v2/2._dp-mu/dag
   
   if(eps.lt.0._dp) then !elliptic (capture)
   
   elseif(eps.eq.0._dp) then !parabolic
   
   else !hyperbolic
   
   a=1._dp/(v2/mu-2._dp/dag)
   
   call crossp3d(rvag(1:3),rvag(4:6),h(1:3))
   
   p=dot_product(h,h)/mu
   
   e=sqrt(1._dp+p/a)
  
   end if
  
   dm=a*abs(1._dp-e)
  
  !calculate time of pericenter passage tpp
  coshf=(dag/a+1._dp)/e
  if (dot_product(rvag(1:3),rvag(4:6)).ge.0._dp) then
  sinhf=sqrt(coshf*coshf-1._dp)
  else
  sinhf=-sqrt(coshf*coshf-1._dp)
  end if
  
  f=log(coshf+sinhf)
  
  tpp=t-(e*sinhf-f)/sqrt(mu/a**3._dp)
  !db
    !db
  !  write(*,*)'tpp,t,dag,a,e,f',tpp,t,dag,dm,dag-dm,f   
    !edb
 !  write(*,*)'t,dm,dag,a,e,eps',t,dm,dag,a,e,eps
   !edb
   return
   end subroutine
!***********************************************************************************

subroutine tpcoo(time,rva,ida,rve,ide,b,tp)
!calculates target plane coordinates
implicit none
integer(kind=idp),parameter::gdim=3

integer(kind=idp),intent(in)::ida
integer,intent(in)::ide
real(kind=dp),intent(in)::time
real(kind=dp),dimension(1:6),intent(in)::rva,rve
real(kind=dp),dimension(1:2),intent(out)::tp
real(kind=dp),intent(out)::b

real(kind=dp)::phi,theta,cost,sint,costm1,cosp,sinp,cospm1,vag
real(kind=dp),dimension(1:6)::rvah,rveh,rvag,rvaxyz
real(kind=dp),dimension(1:3)::zeta,eta,xi,ratp,x,y,z
real(kind=dp),dimension(1:3,1:3)::roty,rotxi,tpbasis

real(kind=dp)::scal

!db
 integer::i

!end db
!heliocentric state vector of asteroid
rvah=rva
call bodykoo1(rvah,ord(idsun),time)
!geocentric state vector of asteroid
rvag=rva
call bodykoo1(rvag,ord(ide),time)
!heliocentric state vector of Earth
rveh=rve
call bodykoo1(rveh,ord(idsun),time)



!xi axis is VxU -> minus of V is incorporated in the cross product
call crossp3d(rveh(4:6),rvag(4:6),xi)

xi=xi/sqrt(dot_product(xi,xi))


!geocentric velocity of the asteroid
vag=sqrt(dot_product(rvag(4:6),rvag(4:6)))

eta=rvag(4:6)/vag

call crossp3d(xi,eta,zeta)

!target plane coordinates
! ratp(1)=dot_product(rvag(1:3),xi(1:3))
! ratp(2)=dot_product(rvag(1:3),eta(1:3))
! ratp(3)=dot_product(rvag(1:3),zeta(1:3))

!xbplane=(xi, eta, zeta)^T x
!tpbasisi is already the transposed matrix
tpbasis(1,:)=xi(:)
tpbasis(2,:)=eta(:)
tpbasis(3,:)=zeta(:)

ratp(1:3)=matmul(tpbasis,rvag(1:3))


! 
! !planetocentric basis following Valsecchi 2003
! x(1:3)=rveh(1:3)/sqrt(dot_product(rveh(1:3),rveh(1:3)))
! y(1:3)=rveh(4:6)/sqrt(dot_product(rveh(4:6),rveh(4:6)))
! call crossp3d(x,y,z)
! 
! z=z/sqrt(dot_product(z,z))
! 
! do i=1,3
!  tpbasis(i,1)=x(i)
!  tpbasis(i,2)=y(i)
!  tpbasis(i,3)=z(i)
! end do
!  
!  call gauss(gdim,tpbasis,rvag(1:3),rvaxyz(1:3))       
!  call gauss(gdim,tpbasis,rvag(4:6),rvaxyz(4:6))       
!  
!  vag=sqrt(dot_product(rvaxyz(4:6),rvaxyz(4:6)))
!  !calculate phi and theta angle rotations following Valsecchi 2003
!  phi=atan2(rvaxyz(4),rvaxyz(6))
!   
!  call crossp3d(rvaxyz(4:6),-rveh(4:6),xi(1:3))
! 
!  xi(1:3)=xi(1:3)/Sqrt(dot_product(xi,xi))
! 
! !rotation should be in -phi direction 
!  cosp=cos(-phi)
!  sinp=sin(-phi)
!  cospm1=1._dp-cosp
! 
!  cost=rvaxyz(5)/vag
!  sint=-sqrt(1._dp-cost*cost)
!  costm1=1._dp-cost
!  
!  !rotation should be in -phi direction 
! !  roty(1,1:3)=(/cosp+y(1)**2*cospm1, y(1)*y(2)*cospm1-y(3)*sinp,&
! !                y(1)*y(3)*cospm1+y(2)*sinp/)
! !  roty(2,1:3)=(/y(2)*y(1)*cospm1+y(3)*sinp,cosp+y(2)**2*cospm1,&
! !                y(2)*y(3)*cospm1-y(1)*sinp/)
! !  roty(3,1:3)=(/y(3)*y(1)*cospm1-y(2)*sinp,y(2)*y(3)*cospm1+y(1)*sinp,&
! !                cosp+y(3)**2*cospm1/)   
!  
! !  roty(1,1:3)=(/1._dp,0._dp,0._dp/)
! !  roty(2,1:3)=(/0._dp,cosp,sinp/)
! !  roty(3,1:3)=(/0._dp,-sinp,cosp/)
!  
!   roty(1,1:3)=(/cosp,0._dp,sinp/)
!  roty(2,1:3)=(/0._dp,1._dp,0._dp/)
!  roty(3,1:3)=(/-sinp,0._dp,cosp/)
!  
!  
!   !rotation should be in -theta direction 
!  rotxi(1,1:3)=(/cost+xi(1)**2*costm1, xi(1)*xi(2)*costm1-xi(3)*sint,&
!                xi(1)*xi(3)*costm1+xi(2)*sint/)
!  rotxi(2,1:3)=(/xi(2)*xi(1)*costm1+xi(3)*sint,cost+xi(2)**2*costm1,&
!                xi(2)*xi(3)*costm1-xi(1)*sint/)
!  rotxi(3,1:3)=(/xi(3)*xi(1)*costm1-xi(2)*sint,xi(2)*xi(3)*costm1+xi(1)*sint,&
!                cost+xi(3)**2*costm1/)                            
! !db
!  write(*,*)'idsun,ide',idsun,ide
!  write(*,*)'helio rva',rva
!  write(*,*)'helio rve',rve
!  write(*,*)'rvah',rvah
!  write(*,*)'rvag',rvag
!  write(*,*)'rveh',rveh
!  write(*,*)'epoch',time
!  STOP
! !edb
! 
! 
! ratp(1:3)=matmul(matmul(rotxi,roty),rvaxyz(1:3))


!db
!   write(*,*)'time',time
!   write(*,*)'rvag',rvag(1:3)
!   write(*,*)'xi',xi(:)
!    write(*,*)'eta',eta(:)
!    write(*,*)'zeta',zeta(:)
!  
!  write(*,*)'basis transformation matrix'
!  do i=1,3
!  write(*,*)tpbasis(i,1:3)
!  end do
!edb

!determine the asteroids position in tp coordinate basis
!yes it works with gauss, but projection is more efficient
! do i=1,3
! tpbasis(i,1)=xi(i)
! tpbasis(i,2)=eta(i)
! tpbasis(i,3)=zeta(i)
! end do
! call gauss(gdim,tpbasis,rvag(1:3),ratp)


! ratp(1)=dot_product(rvag(1:3),xi(1:3))
! ratp(3)=dot_product(rvag(1:3),zeta(1:3))

! !db
!  write(*,*)'vector in new basis',ratp(1),dot_product(rvag(1:3),xi(1:3)),ratp(3),dot_product(rvag(1:3),zeta(1:3))
! ! 
!  STOP
! !edb

!scal=(sqrt(1._dp+ve2(ide)/dot_product(rvaxyz(4:6),rvaxyz(4:6))))
!no scaling of axis
scal=1._dp
!db
write(*,*)'ve2,rvag',ve2(ide),rvag
write(*,*)'scale',scal,'aukm',aukm
!edb
!xi coordinate
tp(1)=ratp(1)/scal
!zeta coordinate
tp(2)=ratp(3)/scal
!calculate scaled impact parameter b=b'/sqrt(xi^2+zeta^2)/sqrt(1+ve^2/vinf^2)
!ve^2=2GM/r_earth is a global variable, and vinf=geocentric velocity 
b=sqrt(tp(1)**2+tp(2)**2)
!db
write(*,*)'xi,zeta,b [au]',tp,b
!edb
end subroutine
!********************************************************************************   

subroutine srmin(r1,r2,v1,v2,dt,tmin,rmin,vmin)
!calculates mimimum distance (rmin) and time of occurence (tmin) from cubic spline interpolation of
!two distances r1,r2 and two corresponding velocities v1=dot r1 and v2= dot r2
implicit none

real(kind=dp),intent(in)::r1,r2,v1,v2,dt
real(kind=dp),intent(out)::tmin,rmin,vmin
real(kind=dp)::a,b,c,d,sq,rddot,discr


if(dt.eq.0.d0) then
rmin=r1
vmin=v1
tmin=0.d0

else

 a=r1
 b=v1
 c=(3.d0*(r2-r1)/dt-2.d0*v1-v2)/dt
 d=(2.d0*(r1-r2)/dt+v1+v2)/(dt*dt)

!calculate t from solving the quadratic equation of the first derivative of the cubic spline

discr=c*c-3.d0*b*d 
!no extreme values in the interval, use bounary values
if (discr.lt.0.d0) then
  if(r1.lt.r2) then
    rmin=r1
    vmin=v1
  else
    rmin=r2
    vmin=v2
  end if
else  
!calculate sqrt used by both solutions
sq=sqrt(discr)
  
if(d.eq.0.d0) then
!quadratic polynomial only
 if (c.eq.0.d0) then
  !linear ploynomial only
  if(b.gt.0.d0) then
     tmin=0.d0
     rmin=a
  else
     tmin=dt
     rmin=a+b*tmin
  end if
 else 
    tmin=-b/(2.d0*c)
    rmin=a+b*tmin+c*tmin**2.d0
    vmin=b+2.d0*c*tmin
 !   write(*,*)"tmin,a,b*tmin,c*timin**2",tmin,a,b*tmin,c*tmin**2.d0
 end if
else
 
!positive solution first
tmin=(-c+sq)/(3.d0*d)

!check whether this is a maximum or minimum by using the second derivative test
 rddot=2.d0*c+6.d0*d*tmin
 !write(*,*)'tmin,rddot',tmin,rddot
 if(rddot.ge.0.d0) then
 !this is the minimum, all fine, calculate rmin from original spline
  rmin=a+b*tmin+c*tmin**2.d0+d*tmin**3.d0
  vmin=b+2.d0*c*tmin+3.d0*d*tmin**2.d0
 else
 !negative solution second
  tmin=(-c-sq)/(3.d0*d)
 
 !check whether this is a maximum or minimum by using the second derivative test
   rddot=2.d0*c+6.d0*d*tmin
  ! write(*,*)'tmin,rddot',tmin,rddot
   if(rddot.ge.0.d0) then
   !this is the minimum, all fine, calculate rmin from original spline
   rmin=a+b*tmin+c*tmin**2.d0+d*tmin**3.d0
    vmin=b+2.d0*c*tmin+3.d0*d*tmin**2.d0
   else
 write(*,*)'Error in routine srmin: no minimum distance could be identified'
 write(*,*)'rddot',rddot
 write(*,*)'tmin',tmin
 write(*,*)'c,d,sq',c,d,sq
 write(*,*)'discr',c*c-3.d0*b*d
 write(*,*)'a,b,c,d',a,b,c,d
 write(*,*)'r1,r2,v1,v2,dt',r1,r2,v1,v2,dt
 STOP
   end if
  end if 
 end if
end if
end if

if (rmin.lt.0.d0) then
 rmin=min(r1,r2)
end if
if (vmin.lt.0.d0) then
 vmin=min(v1,v2)
end if
if (tmin.lt.0.d0) then
tmin=0.d0
end if
if (tmin.gt.dt) then
tmin=dt
end if

return
end subroutine

!*********************************************************************
subroutine srsoi(r1,r2,v1,v2,dt,soi,tsoi)
!calculates the time (tsoi) when a cubic spline interpolated distance reaches the value soi
!dt... timestep (real)
!r1... old distance (real)
!r2... new distance after dt (real)
!v1... old velocity (dot r1) (real)
!v2... new velocity (dot r2) (real)
!soi.. target distance (real)
!tsoi. time when target distance has been reached (real)
!
implicit none
integer::i,j,l(1)
real(kind=dp),intent(in)::r1,r2,v1,v2,dt,soi
real(kind=dp),intent(out)::tsoi
real(kind=dp)::a,b,c,d,sq,tsoitest(1:3),qtest(1:3)

 a=r1
 b=v1
 c=(3.d0*(r2-r1)/dt-2.d0*v1-v2)/dt
 d=(2.d0*(r1-r2)/dt+v1+v2)/(dt*dt)

if(d.eq.0.d0) then
 if(c.eq.0.d0) then
  tsoi=(soi-a)/b
 else
  sq=sqrt(b*b-4.d0*a*c+4.d0*c*soi)
  tsoi=min((-b+sq)/(2.d0*c),(-b-sq)/(2.d0*c))
  if (tsoi.lt.0.d0) then
   tsoi=max((-b+sq)/(2.d0*c),(-b-sq)/(2.d0*c))
  end if
 end if
else

!use three newton raphson trials with 6 steps to identify the minimum
 do j=1,3
 tsoi=real(j)*dt/3.d0
  do i=1,6
   tsoi=tsoi-q(a,b,c,d,tsoi,soi)/qdot(b,c,d,tsoi)
  end do
  
  if(tsoi.le.0.d0) then
   tsoi=0.d0
  end if

  if(tsoi.ge.dt) then
    tsoi=dt
  end if
  
  tsoitest(j)=tsoi
  qtest(j)=q(a,b,c,d,tsoi,soi)
 end do
 
  l=minloc(abs(qtest(:)))
  tsoi=tsoitest(l(1))
 

!  if(tsoi.le.0.d0) then
!   tsoi=0.d0
!  end if
! 
!  if(tsoi.ge.dt) then
!    tsoi=dt
!  end if

!search for solution with netwon raphson
! tsoi=-c/(3.d0*d) - (2.d0**third*(-c**2.d0 + 3.d0*b*d))/ &
!      (3.d0*d*(-2.d0*c**3.d0 + 9.d0*b*c*d - 27.d0*a*d**2.d0 + 27.d0*d**2.d0*soi + &
!   Sqrt(4.d0*(-c**2.d0 + 3.d0*b*d)**3.d0 + (-2.d0*c**3.d0 + 9.d0*b*c*d - 27.d0*a*d**2.d0 + 27.d0*d**2.d0*soi)**2.d0))**third) &
!   + (-2.d0*c**3.d0 + 9.d0*b*c*d - 27.d0*a*d**2.d0 + 27.d0*d**2.d0*soi + &
!   Sqrt(4.d0*(-c**2.d0 + 3.d0*b*d)**3.d0 +(-2.d0*c**3.d0 + 9.d0*b*c*d - 27.d0*a*d**2.d0 + 27.d0*d**2.d0*soi)**2.d0))**third/ &
!   (3.d0*2.d0**third*d)
    
end if             
return
end subroutine 

!************************************************
subroutine srvint(rv1,rv2,dt,tint,rvint)
!interpolates 3D distance and velocity state vectors given at endpoints rv1 and rv2 
!to produce interpolated value rvint at tint
!dt is the timestep between rv1 and rv2
implicit none

real(kind=dp),intent(in),dimension(1:6)::rv1,rv2
real(kind=dp),intent(in)::tint,dt
real(kind=dp),intent(out),dimension(1:6)::rvint
real(kind=dp)::a(3),b(3),c(3),d(3),sq(3),tint2,tint3

tint2=tint*tint
tint3=tint2*tint

 a(1:3)=rv1(1:3)
 b(1:3)=rv1(4:6)
 c(1:3)=(3.d0*(rv2(1:3)-rv1(1:3))/dt-2.d0*rv1(4:6)-rv2(4:6))/dt
 d(1:3)=(2.d0*(rv1(1:3)-rv2(1:3))/dt+rv1(4:6)+rv2(4:6))/(dt*dt)
 
 
 rvint(1:3)=a(1:3)+b(1:3)*tint+c(1:3)*tint2+d(1:3)*tint3
 rvint(4:6)=b(1:3)+2.d0*c(1:3)*tint+3.d0*d(1:3)*tint2
 
 return
 end subroutine
 
  



! 
! subroutine tpcoo(time,rva,ida,rve,ide,b,tp)
! !calculates target plane coordinates
! implicit none
! integer(kind=idp),parameter::gdim=3
! 
! integer(kind=idp),intent(in)::ida
! integer,intent(in)::ide
! real(kind=dp),intent(in)::time
! real(kind=dp),dimension(1:6),intent(in)::rva,rve
! real(kind=dp),dimension(1:2),intent(out)::tp
! real(kind=dp),intent(out)::b
! 
! real(kind=dp),dimension(1:6)::rvah,rveh,rvag
! real(kind=dp),dimension(1:3)::zeta,eta,xi,ratp
! real(kind=dp),dimension(1:3,1:3)::tpbasis
! 
! real(kind=dp)::scal
! 
! !db
!  integer::i
! 
! !end db
! 
! rvah=rva
! call bodykoo1(rvah,ord(idsun),time)
! 
! rvag=rva
! call bodykoo1(rvag,ord(ide),time)
! 
! rveh=rve
! call bodykoo1(rveh,ord(idsun),time)
! 
! !db
! !  write(*,*)'idsun,ide',idsun,ide
! !  write(*,*)'helio rva',rva
! !  write(*,*)'helio rve',rve
! !  write(*,*)'rvah',rvah
! !  write(*,*)'rvag',rvag
! !  write(*,*)'rveh',rveh
! !edb
! 
! !fist basis vector is orthogonal to the b-plane, parallel to the asteroid geocentric velocity 
! eta(1:3)=rvag(4:6)
! !normalize first basis vector
! eta(1:3)=eta(1:3)/sqrt(dot_product(eta,eta))
! 
! !second basis vector is the negative heliocentric velocity of earth projected onto the b-plane
! zeta(1:3)=eta(1:3)*dot_product(rveh(4:6),eta(1:3))-rveh(4:6)
! !normalize second basis vector
! zeta(1:3)=zeta(1:3)/sqrt(dot_product(zeta,zeta))
! 
! !construct third basis vector as xi=eta x zeta
! call crossp3d(eta,zeta,xi)
! !check whether xi is normalized
! 
! xi(:)=xi(:)/sqrt(dot_product(xi,xi))
! 
! if(abs(sqrt(dot_product(xi,xi))-1.d0).gt.1.d-12) then
!  write(*,*)'target plane xi coordinate not normalized'
! end if
! if(abs(dot_product(xi,zeta)).gt.1.d-12) then
!  write(*,*)'target plane xi and zeta axes non orthogonal'
! end if
! 
! 
! 
! !db
! !   write(*,*)'time',time
! !   write(*,*)'rvag',rvag(1:3)
! !   write(*,*)'xi',xi(:)
! !    write(*,*)'eta',eta(:)
! !    write(*,*)'zeta',zeta(:)
! !  
! !  write(*,*)'basis transformation matrix'
! !  do i=1,3
! !  write(*,*)tpbasis(i,1:3)
! !  end do
! !edb
! 
! !determine the asteroids position in tp coordinate basis
! !yes it works with gauss, but projection is more efficient
! ! do i=1,3
! ! tpbasis(i,1)=xi(i)
! ! tpbasis(i,2)=eta(i)
! ! tpbasis(i,3)=zeta(i)
! ! end do
! ! call gauss(gdim,tpbasis,rvag(1:3),ratp)
! 
! 
! ! ratp(1)=dot_product(rvag(1:3),xi(1:3))
! ! ratp(3)=dot_product(rvag(1:3),zeta(1:3))
! 
! ! !db
! !  write(*,*)'vector in new basis',ratp(1),dot_product(rvag(1:3),xi(1:3)),ratp(3),dot_product(rvag(1:3),zeta(1:3))
! ! ! 
! !  STOP
! ! !edb
! 
! scal=(sqrt(1._dp+ve2(ide)/dot_product(rvag(4:6),rvag(4:6))))
! !db
! write(*,*)'ve2,rvag',ve2(ide),rvag
! write(*,*)'scale',scal,'aukm',aukm
! write(*,*)'rvxi,rvzeta',dot_product(rvag(1:3),xi(1:3)),dot_product(rvag(1:3),zeta(1:3))
! !edb
! tp(1)=dot_product(rvag(1:3),xi(1:3))/scal !xi
! tp(2)=dot_product(rvag(1:3),zeta(1:3))/scal !zeta
! !calculate scaled impact parameter b=b'/sqrt(xi^2+zeta^2)/sqrt(1+ve^2/vinf^2)
! !ve^2=2GM/r_earth is a global variable, and vinf=geocentric velocity 
! b=sqrt(tp(1)**2+tp(2)**2)
! !db
! write(*,*)'tp(1:2),b [au]',tp,b
! !edb
!end subroutine
  

!***********************************************
   end module