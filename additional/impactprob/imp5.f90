program imp4
!uses orbfit .fou and nine .dmi output files and calculates Impact probability
implicit none 

integer::nclones,eof,n,n2,i,j,nVIh,nVIs,a,mon,day,np
real*8::dumr(1:8),dumr2(11),pi,rad2grad,rearthau,sigstart
real*8::dsig,au,cursig,rearth,dimp,dmin,dmax,dd
real*8,dimension(:),allocatable::d,prob
real*8,dimension(:,:),allocatable::sigd
character(len=50)::dum,ini,infile,infile1,infile2,outfile1,outfile2
character(len=150)::out
character(len=10)::color

write(*,*)'name of asteroid'
read(*,*)infile
write(*,*)'name of massless.inn file'
read(*,*)infile1
write(*,*)'name of nine .dmi file'
read(*,*)infile2
write(*,*)'number of clones'
read(*,*)nclones
write(*,*)'number of distance bins for probability'
read(*,*)np
write(*,*)'# sigma of clone sampling'
read(*,*)sigstart

!sigma sampling
dsig=(2.d0*sigstart)/real(nclones)

allocate(d(np),prob(np),sigd(nclones,4))
d(:)=0.d0
prob(:)=0.d0

sigd(:,:)=0.d0

!get sigma LOV values for clones
!ini=" no  RMS ,lastcor,  magn,  MOID ,nod+,nod-, sigQ"

outfile1="imps.out"
outfile2="imph.out"
!out="n    sigma calc   sigma fou   dmin [Rearth]   tmin [JD]"

pi=4.d0*atan(1.d0)
rad2grad=180.d0/pi
au=149597870.7d0
rearth=6371.d0
rearthau=6371.d0/au 

!define impact distacne
dimp=rearthau

open(21,file=trim(infile1),status="old")
open(22,file=trim(infile2),status="old")
open(23,file=trim(outfile1),status='replace')
open(24,file=trim(outfile2),status='replace')
write(23,*)out
write(24,*)out


open(25,file='imp.gnu',status='replace')
call gnuhead(25,-sigstart,sigstart,infile,rearth,au)



eof=0
!number of virtual impactors due to hyperbolic min dist
nVIh=0
!number of virtual impactors due to spline min dist
nVIs=0

!read header
! do while (dum.ne.trim(ini))
!  read(21,'(A)',iostat=eof)dum
!  if(eof.ne.0) then
!   write(*,*)'no initialization statement detected to initiate extraction'
!   STOP
!  end if
! end do
read(21,*)dum
read(22,*)dum

!start comparing files
n=1
do while(n<nclones)
read(21,*,iostat=eof)dumr(:)

!  if(n.le.1) then
!   !sigstart=dumr(8)
    
!   write(*,*)'LOV sigma start,dsigma',sigstart,dsig
!  end if
 
 if(eof.ne.0) then
  write(*,*)'error in reading clone elements, possibly number of clones incorrect'
  STOP
 end if
 
read(22,*,iostat=eof)n2,dumr2(:) 
 if(eof.ne.0) then
  write(*,*)'error in reading clone elements, possibly number of clones incorrect'
  STOP
 end if
 
 if(n2.ne.n) then 
 write(*,*)'possible mismatch between input files asteroid ids'
 end if
 
 cursig=dumr(8)
 !sigma
 sigd(n,1)=cursig
 !spline d
 sigd(n,2)=dumr2(4)
 !hyp d
 sigd(n,3)=dumr2(7)
 !jd
 sigd(n,4)=dumr2(5)
 
 !spline
 if(abs(dumr2(4)).lt.dimp) then

  call jd2greg(dumr2(5),a,mon,day)
  if(a>2013) then
  color='red'
  call setlabgnp(25,2040,a,cursig,dumr2(4)/rearthau,color)
  end if
  write(23,*)n,cursig,dumr(8),dumr2(4)/rearthau, dumr2(5),a,mon,day
  nVIs=nVIs+1
 end if

 
 if(abs(dumr2(7)).lt.dimp) then 
  call jd2greg(dumr2(5),a,mon,day)
  if (a>2013) then
  color='blue'
    call setlabgnp(25,2015,a,cursig,dumr2(7)/rearthau,color)
  end if
  write(24,*)n,cursig,dumr(8),dumr2(7)/rearthau, dumr2(8),a,mon,day
  nVIh=nVIh+1
 end if
 
 n=n+1
end do
if(nVIs.eq.0) then
write(*,*)'no VIs found'
else
write(*,*)'nVIs/nclones',real(nVIs)/real(nclones),'1/',int(real(nclones)/real(nVIs))
end if

if(nVIh.eq.0) then
write(*,*)'no VIh found'
else
write(*,*)'nVIh/nclones',real(nVIh)/real(nclones),'1/',int(real(nclones)/real(nVIh))
end if

call gnutail(25,nclones,-sigstart,sigstart,infile2)
!######################################
!calculate cumulative probability vs distance
!######################################
!construct distance grid
dmin=1.d5
dmax=0.d0
do i=1,nclones
if(sigd(i,2).gt.1.d-12) then
dmin=min(dmin,sigd(i,2))
end if

dmax=max(sigd(i,2),dmax)

end do
dd=(dmax-dmin)/real(np)

d(1)=dmin
d(np)=dmax

do i=2,np-1
 d(i)=dmin+real(i)*dd
end do

!calculate cumulative probabilities for distances (only valid for gaussian clones and equidistribution of mitigation parameters)


do j=1,nclones
 do i=1,np
  if(sigd(j,2).lt.d(i)+dd.and.sigd(j,2).ge.d(i)) then
   prob(i)=prob(i)+0.5d0*(-erf(sigd(j,1)/sqrt(2.d0))+erf((sigd(j,1)+dsig)/sqrt(2.d0)))
 end if
 end do
end do

open(26,file='prob.dat',status='replace')
do j=1,np
write(26,*)j,d(j),sum(prob(1:j)),prob(j)
end do
close(26)

open(26,file='prob.gnu',status='replace')
call gnuprob(26,infile,rearth,au)
close(26)


close(21)
close(22)
close(23)
close(24)
close(25)
end program 
!********************************************************
subroutine gnuprob(un,tit,rearth,au)
!generates header for gnuplot plot file 
implicit none
integer::un
real*8::ss,se,rearth,au
character(len=50)::tit
!write(un,*)'set term postscript enhanced  color font "Helvetica, 16" landscape'
write(un,*)'set term postscript color enhanced'
write(un,*)'set size 0.9,1'
write(un,*)'rearth=',rearth
write(un,*)'au=',au
write(un,*)'rpau=rearth/au'
write(un,*)'set title "',trim(tit),'"'
!write(un,*)'set logscale y'
write(un,*)'set key bottom right'
write(un,*)'set output "d_prob.eps"'
write(un,*)'set xlabel "Probability"'
write(un,*)'set xlabel "Minimum distance [R_{Earth}]"'
write(un,*)'p[:][0:1] "prob.dat"  u ($2/rpau):3 notit w l lw 3'
write(un,*)'set output "dum.eps"'
write(un,*)'set term x11'

return
end subroutine
!*************************************************************
subroutine gnuhead(un,ss,se,tit,rearth,au)
!generates header for gnuplot plot file 
implicit none
integer::un
real*8::ss,se,rearth,au
character(len=50)::tit
!write(un,*)'set term postscript enhanced  color font "Helvetica, 16" landscape'
write(un,*)'set term postscript color enhanced'
write(un,*)'set size 0.9,1'
write(un,*)'rearth=',rearth
write(un,*)'au=',au
write(un,*)'rpau=rearth/au'
write(un,*)'set title "',trim(tit),'"'
!write(un,*)'set logscale y'
write(un,*)'set key bottom right'
write(un,*)'set output "rearth_sigma.eps"'
write(un,*)'set arrow from', ss,',1 to ',se,',1 nohead' 
write(un,*)'set xlabel "{/Symbol s}_{LOV}"'
write(un,*)'set ylabel "Min. distance [R_{Earth}]"'

return
end subroutine
!*************************************************************
subroutine gnutail(un,nc,ss,se,infile2)
!generates tail for gnuplot plot file 
implicit none
integer::un,nc,ncmid
real*8::ss,se
character(len=50)::infile2
ncmid=int(real(nc)/2.d0+4)

write(un,*)'p[',ss,':',se,'][:] "',trim(infile2),'"  u (',-abs(ss-se),'*(',ncmid,'-$1)/',nc,'):($8/rpau)', &
           'tit "hyperbolic" w lp ps 0.9, "" u (',-abs(ss-se),'*(',ncmid,'-$1)/',nc,'):($5/rpau) ', &
           'tit "cubic spline" w lp ps 0.5'
write(un,*)'set output "dum.eps"'
write(un,*)'set term x11'
return
end subroutine

!*************************************************************
subroutine jd2greg(jd,a,mon,day)
!calculates Gregorian year, month, day from JD following Richards
implicit none
integer,intent(out)::a,mon,day
real*8,intent(in)::jd

integer,parameter::y=4716,n=12,m=2
real*8,parameter::j=1401,r=4,p=1461
real*8,parameter::v=3,u=5,s=153,w=2,B=274277,C=-38
real*8::f,e,g,h

f=jd+j
f=f+real(int((int((4.d0*jd+B)/146097.d0)*3.d0)/4.d0))+C
e=r*f+v
g=real(int(mod(e,p)/r))
h=u*g+w
day=int(mod(h,s)/u)+1
mon=mod(int(real(h)/real(s))+m,n)+1
a=int(e/p)-y+int(real(n+m-mon)/real(n))

return
end subroutine
!****************************************************************
subroutine setlabgnp(unit,aa,a,sig,dmin,color)
!generates set label list for gnuplot containing location for labels and year of impact solution as label
implicit none
integer,intent(in)::unit,aa,a
!aa avoid impact date
real*8,intent(in)::sig,dmin
character(len=10)::color

if(a.ne.aa) then
write(unit,'((A),I4.4,(A),X,F15.10,(A),X,F15.10,3(A))')'set label "',a,'" at ',sig,',',&
           ! max(abs(dmin+100.*(rand()-0.5d0)),10.),&
            dmin+dmin*(rand()-0.5d0)!,' textcolor rgbcolor "',trim(color),'"' 
end if
return
end subroutine





