module output_m
  use global_m                      
  use transform_m
  use getephem_m

  implicit none

  public::OpenOutput
  public::CloseOutput
  public::Out
  public::ccout
  public::locccdmp
  public::mdistdmp
  
  character(len=32),private,dimension(1:3,1:11)::fname


!/////////////////////////////////////////////////////////////////////////////////////////////////
!                                           Ausgabe und Auswertung
!
!                                         by Siegfried Eggl  20080110 
!------------------------------------------------------------------------------
!**********************************************************************
contains
  
!**********************************************************************
  subroutine OpenOutput(intname)
  use global_m
                       implicit none
                        integer::bu,du,dc
                        integer(kind=idp)::tdigit,i,j
                        real(kind=dp)::dtd,nn
                        
                        character(len=5),intent(in)::intname
                        character(len=4),dimension(1:3)::filemid
                        character(len=4),dimension(1:11)::fileend
                        character(len=20)::acc,form

                          write(unit=*,fmt=*)'output file(s) will be:'
!------------------------------------------------------------------------
!For Binary output: use gfortran or ifort(not supported???)

      acc = 'stream'
      form = 'unformatted'
!---------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------
!defines base unit for output files and the delta unit to augment the unit integer for opening files
! the latter needs to be the same as the filemid dimension
     bu=22
     du=3
     dc=0
     
!--------------------------------------------------------------------------------------------------     
! here the file mid sections and endings are defined for massive bodies (_big) massless bodies (_sml)
! and ephemeris data (_eph)
!-----------------------------------------------------------------------------------------------------------
                        
                        filemid(1) ='_sml'
                        filemid(2) ='_big'
                        filemid(3)='_eph'
                        
                        
                        fileend(1) = '.bco'
                        fileend(2) = '.hco'
                        fileend(3) = '.hel'
                      
                        fileend(4) = '.bbc'
                        fileend(5) = '.bhc'
                        fileend(6) = '.bhe'
                        
                        fileend(7)='.cen'
                        fileend(8)='.dmi'
                        fileend(9)='.kin'
                        
                        fileend(10)='.egy'
                        fileend(11)='.loc'
                        
                        do i=1,11
                         do j=1,3
                          fname(j,i)=intname//filemid(j)//fileend(i)
                         end do 
                        end do
  
  
!-------------------------- ASCII output-------------------------------------------
!--------------- BCO --------------------------------------
do j=1,du
  uc(j,1)=bu+dc
  if(outbco) then
        open(unit=uc(j,1), file=fname(j,1), status="replace", action="write")
                        write(unit=uc(j,1),fmt=1001)
                        write(unit=*,fmt=*)fname(j,1)
  end if             
  dc=dc+1
end do    


!----------------- HCO -----------------------------------   
do j=1,du
!current unit
 uc(j,2)=bu+dc  
      if(outhco) then

        open(unit=uc(j,2), file=fname(j,2), status="replace", action="write")   
                        write(unit=uc(j,2),fmt=1001)
                        write(unit=*,fmt=*)fname(j,2) 
      end if 
  dc=dc+1                
end do
 

!----------------- HEL -----------------------------------        
do j=1,du
  uc(j,3)=bu+dc
  if(outhele) then
        open(unit=uc(j,3), file=fname(j,3), status="replace", action="write")   
                if(input.eq.'ye'.or. input.eq.'Ye'.or.input.eq.'YE') then
                        write(unit=uc(j,3),fmt=1007)
                else  
                        write(unit=uc(j,3),fmt=1002)
                end if
                        write(unit=*,fmt=*)fname(j,3)
      
  end if
  dc=dc+1
end do
               
      
!--------------------------Binary output-------------------------------------------
!----------------- BBC ----------------------------------- 
do j=1,du
  uc(j,4)=bu+dc
  if(outbinbco) then
        open(unit=uc(j,4), file=fname(j,4), status='replace',form=form,Access=acc)   
                    write(unit=*,fmt=*)fname(j,4)       
  end if     
 dc=dc+1
end do 
!----------------- BHC -----------------------------------   
do j=1,du
  uc(j,5)=bu+dc                 
if(outbinhco) then
        open(unit=uc(j,5), file=fname(j,5), status='replace',form=form,Access=acc)   
                        write(unit=*,fmt=*)fname(j,5)
end if 
  dc=dc+1
end do
    
!----------------- BHE -----------------------------------   
do j=1,du
  uc(j,6)=bu+dc    
if(outbinhele) then
        open(unit=uc(j,6), file=fname(j,6),  status='replace',form=form,Access=acc)
                    write(unit=*,fmt=*)fname(j,6)
                   
end if
 dc=dc+1
end do

     !--------------------------Close Encounter output-------------------------------------------         
 !----------------------  CEN ------------------------------
 do j=1,du-1
  uc(j,7)=bu+dc     
 if(outcc) then  
      open(unit=uc(j,7), file=fname(j,7), status='replace',action="write")   
                         write(unit=*,fmt=*)fname(j,7)   
                         write(unit=uc(j,7),fmt=1006)   
                                         
 end if 
  dc=dc+1
end do

!--------------------- DMI ----------------------------
 do j=1,du-1
  uc(j,8)=bu+dc     
 if(outcc) then  
       open(unit=uc(j,8), file=fname(j,8), status='replace',action="write")   
                         write(unit=*,fmt=*)fname(j,8)   
                         write(unit=uc(j,8),fmt=1007)     
                         
   end if                        
  dc=dc+1
end do                         

!-------------------- KIN ------------------------------
  uc(1,9)=bu+dc 
 if(kick.eq.'t'.or.kick.eq.'p') then
	open(unit=uc(1,9), file=fname(1,9), status='replace', action="write")
                         write(unit=uc(1,9),fmt=1008)
                         write(unit=*,fmt=*)fname(1,9)
end if
  dc=dc+1   	         
 
!-------------------- INTEGRALS OF MOTION ------------------------------------------- 
 do j=1,du
  uc(j,10)=bu+dc  
if(outegy) then
        open(unit=uc(j,10), file=fname(j,10), status="replace", action="write")
                        write(unit=uc(j,10),fmt=1003)
                        write(unit=*,fmt=*)fname(j,10)
end if   	         
  dc=dc+1
end do   	         

!-------------------- Local Close encounter data ------------------------------------------- 

  uc(1,11)=bu+dc  
        open(unit=uc(1,11), file=fname(1,11), status="replace", action="write")
                        write(unit=uc(1,11),fmt=1009)
                        write(unit=*,fmt=*)fname(1,11)  	         
  dc=dc+1
!-------------------- SHOW PROGRESS IN TERMINAL ------------------------------------------- 

if(showprog) then
write(unit=*,fmt=*)'progress:'
else
 write(unit=*,fmt=*)'working...'      
end if


!---------------------Fix time output format----------------

dtd=huge(dtd)
nn=1._dp
tdigit=1
if(tout.lt.1._dp) then
 do while(dtd.gt.1.d-13)
  tdigit=tdigit+1
  nn=nn*10._dp
  dtd=nn*tout-real(int(nn*tout,kind=idp),kind=dp)
 end do 
else
 do while(dtd.gt.1.d-13)
   tdigit=tdigit+1
   nn=nn/10._dp
  dtd=nn*tout-real(int(nn*tout,kind=idp),kind=dp)
 end do 
end if

 if (tdigit.gt.0.and.tdigit.lt.100) then
   write(tdigc,"(I2)")tdigit
 else 
   tdigc="00"
 end if
!define global output format vars tdigc, tendc
write(tdigc,"(I2)")int(log10(tend))+tdigit
write(tendc,"(I2)")int(log10(tend))+tdigit+7


!------------output format statements---------------------------------------

1001 format(10X,'time',15X,'body',20X,'x [AU]',20X,'y [AU]',20X,'z [AU]',20X, &
             'vx [AU/D]',20X,'vy [AU/D]',20X,'vz [AU/D]',20X,'mass [M_sun]',20X,'name')

1002 format(10X,'time',15X,'body',25X,'a',25X,'e',25X,'i',25X,'omega',25X,'OMEGA', &
              25X,'mean anomaly',20X,'mass [M_sun]',20X,'name')

1003 format( 18X, 'time',20X,'total energy',10X,'delta energy',10X,'kinetic energy',10X,'potential energy',& 
                             10X, 'sum of delta energy', &
                             10X,'total angular momentum', 10X,'delta t.a.m',10X,'sum of delta t.a.m.', &
                             20X,'t.a.m. x',20X,'t.a.m. y',20X,'t.a.m. z', & 
                             10X,'barycenter x',10X,'barycenter y',10X,'barycenter z',20X)
1004 format(10X,'time',15X,'body',25X,'x',25X,'y',25X,'z',25X,'vx',25X,'vy',25X,'vz',25X,'specific energy')

1005 format(10X,'time',15X,'body1',5X,'body2',10X,'x',15X,'y',15X,'z',15X,'v_relative/v_escape')

1006 format(15X,'body1',20X,'body2',5X, 'scaled impact param b [km]',5X, 'target plane xi [km]', &
                 5X,'target plane zeta [km]',5X,'time [JD]',10X, 'distance [au]',10X, &
                 'time of pericenter passage [JD]',10X,'pericenter dist [au]')
                 
1007 format(15X,'body',5X, 'first bplane xi [km]',5X, 'first bplane zeta [km]',5X, &
                'first bplane epoch [JD]',5X,'min dist spline [au]',5X, 'epoch of min dist [JD]', 5X, &
                 'vrel at closest encounter [km/s]',5X,  'pericenter dist hyperbolic fit [au]', &
                 5X,'epoch of pericenter passage [JD]',5X,'vinf [km/s]',5X,'number of close encounters',5X, &
                 'number of deep close encounters',5X,'date of closest encounter')         
1008 format(15X,'body id',5X, 'kick time [JD]', 10X, 'asteroid r [km]', 10X, 'asteroid v before impact [km/s]', &
                5X,'asteroid v after impact [km/s]', &
                5X,'change in asteroid velocity dV (3) [km/s]',5X,  'relative velocity S/C - asteroid (3) [km/s]', &
                5X,'ejecta momentum (3) [kg km/s]',5X,'S/C impact velocity in ICRF (3) [km/s]' , &
                'phi (angle between impact velocity and ejecta vector [deg]')
                
1009 format(5X,'number of close encounter',5X, 'body id',5X, 'time of closest encounter [JD]', 10X, 'dmin [au]', &
                5X,'v_min (3) [au/D]',5X,  'bplane xi [km]', &
                5X,'bplane zeta [km]',5X,'time of bplane [JD]' , &
                'time of closest encounter yr, month, day')                
                            
! 1007 format(10X,'time',15X,'body',25X,'a',25X,'e',25X,'i',25X,'omega',25X,'OMEGA', &
!               20X,'mean anomaly',5X,'Yark da/dt diurnal [au/Myr]',3X, 'Yark < da/dt > diurnal [au/Myr]', &
!               5X,'mass [M_sun]',8X,'name')                            
 

end subroutine 

!*************************************************
  subroutine CloseOutput
  implicit none
  integer::i,j
   do i=1,11
    do j=1,3
      close(uc(j,i))
     end do
   end do
                 
     !                  write(*,*)'all units closed'
  end subroutine 
!***************************************************
subroutine Out(tt,rvin,massin)
!routine that is responsible for file output
use global_m

 !output routine for all integrators
real(kind=dp),intent(in)::tt
!type(vec),dimension(:),intent(in)::r,v
real(kind=dp),dimension(:),intent(in)::massin
real(kind=dp),dimension(1:inbody,1:6),intent(in)::rvin
real(kind=dp),dimension(1:inbody,1:6)::rv,ele,rvh,rvb
real(kind=dp),dimension(1:neph,1:6)::rveph,eleph
real(kind=dp),dimension(1:inbody)::m
!real(kind=dp)::engy,ekin,epot,lpsum,lpsum2(1:3),rb(1:6)
!real(kind=dp), dimension(:,:),allocatable::lp
!real(kind=dp), dimension(1:p%NrKoerper,1:3)::lrl
real(kind=dp)::lperror,engyerror,t
!real(kind=dp)::omega1,omega2
!real(kind=dp)::spezen,spezenel
integer(kind=idp)::n,mp1,nephm1,i,j

mp1=mcount+1
!assuming the Sun is the last in the list of ephemeris bodies
nephm1=neph-1
t=tt+tstart
!------------------- BARYCENTRIC COORDINATES ---------------------------
 if(outbco) then   
         !make local copy to avoid overwriting in different output files
         rvb(:,:)=rvin(:,:)
         m(:)=massin(1:mcount+m0run)
         !massless
      do n=mp1,inbody
         !Barycentric
        write(unit=uc(1,1),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,(n-mcount)+m0cur      
        write(unit=uc(1,1),fmt="(7(f25.16,3X),A)",advance="yes") rvb(n,1:6),0.d0                     
      end do   
         !massive 
      if(massiveout)   then
      !massive file
        do n=1,mcount
         !Barycentric
        write(unit=uc(2,1),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n !(n+m0count)      
        write(unit=uc(2,1),fmt="(7(f25.16,3X),A)",advance="yes") rvb(n,1:6),m(n)                     
        end do
      end if
        
       !ephemeris 
       if(ephemout) then
        call getrv(t,rveph)
          do n =1,neph
            write(unit=uc(3,1),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n !(totbody+n)      
            write(unit=uc(3,1),fmt="(7(f25.16,3X),A)",advance="yes") rveph(n,1:6),meph(n)   
           end do
       end if       
 end if
  
!------------------ HELIOCENTRIC COORDINATES --------------------------------   
 if(outhco) then
       !make local copy to avoid overwriting in different output files
         rvh(:,:)=rvin(:,:)
         m(:)=massin(1:mcount+m0run)
      !massless
      !transformation into heliocentric coordinates 
      call hekoo(rvh(mp1:inbody,1:6),m0count,t)
      do n=mp1,inbody
        write(unit=uc(1,2),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,(n-mcount)   +m0cur   
        write(unit=uc(1,2),fmt="(7(f25.16,3X),A)",advance="yes") rvh(n,1:6),0.d0                     
      end do   
      !massive 
      if(massiveout)   then
       !transformation into heliocentric coordinates 
      call hekoo(rvh(1:mcount,1:6),mcount,t)
        do n=1,mcount
        write(unit=uc(2,2),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n !(n+m0count)      
        write(unit=uc(2,2),fmt="(7(f25.16,3X),A)",advance="yes") rvh(n,1:6),m(n)                     
        end do
      end if
       !ephemeris 
       if(ephemout) then
        call getrvhel(t,rveph)
        call icrs2he(rveph,nephm1)
          do n =1,nephm1
            write(unit=uc(3,2),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n !(totbody+n)      
            write(unit=uc(3,2),fmt="(7(f25.16,3X),A)",advance="yes") rveph(n,1:6),meph(n)   
           end do
       end if       
   end if
       
!------------------ KEPLERIAN ORBITAL ELEMENTS (HELIOCENTRIC) ---------------------------       
   if(outhele) then
         !make local copy to avoid overwriting in different output files
         rvh(:,:)=rvin(:,:)
         m(:)=massin(1:mcount+m0run)
      !massless
      !transformation into heliocentric elements
      call hekoo(rvh(mp1:inbody,1:6),m0run,t)
      call htrnsel(rvh(mp1:inbody,1:6),ele(mp1:inbody,1:6),m(mp1:inbody),m0run)
      do n=mp1,inbody                                          
         write(unit=uc(1,3),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,(n-mcount)+m0cur  
          write(unit=uc(1,3),fmt="(7(f25.16,3X),A)",advance="yes")ele (n,:),m(n)      
      end do       
     !massive 
      if(massiveout)   then
       !transformation into heliocentric coordinates 
      call hekoo(rvh(1:mcount,1:6),mcount,t)
      call htrnsel(rvh(1:mcount,1:6),ele(1:mcount,1:6),m(1:mcount),mcount)
        do n=1,mcount
        write(unit=uc(2,3),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n !(n+m0count)      
        write(unit=uc(2,3),fmt="(7(f25.16,3X),A)",advance="yes") ele(n,:),m(n)                     
        end do
      end if
       !ephemeris 
       if(ephemout) then
        call getrvhel(t,rveph)
        call icrs2he(rveph,nephm1)
        call htrnsel(rveph(1:nephm1,1:6),eleph(1:nephm1,1:6),meph(1:nephm1),nephm1)
          do n =1,nephm1
            write(unit=uc(3,3),fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n !(totbody+n)      
            write(unit=uc(3,3),fmt="(7(f25.16,3X),A)",advance="yes") eleph(n,:),meph(n)   
           end do
       end if          
   end if    
       
       
!------------------------------------------------
!-------- BINARY FILE OUTPUT --------------------
!------------------------------------------------

!----------BARYCENTRIC COORDINATES --------------
if(outbinbco) then

    rvb(:,:)=rvin(:,:)
    m(:)=massin(:)      
         !massless
      do n=mp1,inbody
         !Barycentric
        write(unit=uc(1,4)) t,(n-mcount)+m0cur  , rvb(n,1:6),0.d0                     
      end do   
         !massive 
      if(massiveout)   then
      !massive file
        do n=1,mcount
         !Barycentric
        write(unit=uc(2,4)) t,n,rvb(n,1:6),m(n)                     
        end do
      end if     
       !ephemeris 
       if(ephemout) then
        call getrv(t,rveph)
          do n =1,neph
            write(unit=uc(3,4)) t,n,rveph(n,1:6),meph(n)   
           end do
       end if        
   end if
   
!---------- HELIOCENTRIC COORDINATES --------------
if(outbinhco) then
    !make local copy to avoid overwriting in different output files
    rvh(:,:)=rvin(:,:)
    m(:)=massin(:)      
         !massless
      call hekoo(rvh(mp1:inbody,1:6),m0run,t)
      do n=mp1,inbody
         !Barycentric
        write(unit=uc(1,5)) t,(n-mcount)+m0cur  , rvh(n,:),0.d0                     
      end do   
         !massive 
      if(massiveout)   then
      call hekoo(rvh(1:mcount,1:6),mcount,t)
      !massive file
        do n=1,mcount
         !Barycentric
        write(unit=uc(2,5)) t,n,rvh(n,:),m(n)                     
        end do
      end if     
       !ephemeris 
       if(ephemout) then
        call getrvhel(t,rveph)
         call icrs2he(rveph,neph)
          do n =1,neph-1
            write(unit=uc(3,5)) t,n,rveph(n,1:6),meph(n)   
           end do
       end if        
   end if
   
!----------KEPLERIAN ORBITAL ELEMENTS (HELIOCENTRIC, OSCULATING) -------------- 
if(outbinhele) then
         !make local copy to avoid overwriting in different output files
         rvh(:,:)=rvin(:,:)
         m(:)=massin(:)
      !massless
      !transformation into heliocentric elements
      call hekoo(rvh(mp1:inbody,1:6),m0run,t)
      !transformation into heliocentric elements
      call htrnsel(rvh(mp1:inbody,1:6),ele(mp1:inbody,1:6),m(mp1:inbody),m0run)
      do n=mp1,inbody                                          
         write(unit=uc(1,6)) t,(n-mcount)+m0cur,ele (n,:),m(n)      
      end do       
     !massive 
      if(massiveout)   then
       !transformation into heliocentric coordinates 
      call hekoo(rvh(1:mcount,1:6),mcount,t)
       !transformation into heliocentric elements
      call htrnsel(rvh(1:mcount,1:6),ele(1:mcount,1:6),m(mp1:inbody),mcount)
        do n=1,mcount
        write(unit=uc(2,6)) t,n, ele(n,:),m(n)                     
        end do
      end if
       !ephemeris 
       if(ephemout) then
        call getrvhel(t,rveph)
          call icrs2he(rveph,neph)
        !transformation into heliocentric elements
        call htrnsel(rveph(1:neph,1:6),eleph(1:neph,1:6),meph(1:neph),neph)
          do n =1,neph-1
            write(unit=uc(3,6)) t,n, eleph(n,:),meph(n)   
           end do
       end if   
end if

!---------------  TERMINAL OUTPUT OF PROGRESS -----------------------
if(showprog) then
       call  progress(tt,tcalc)
end if
!---------------------------------------------------------------------
 return
  end subroutine Out
!*********************************************************
subroutine ccout(i,id,t,d,tpp,dm,b,tp)
!output in close encounter file for target plane analysis
implicit none
integer,intent(in)::id
integer(kind=idp),intent(in)::i
integer::yr,mon,day
real(kind=dp),intent(in)::t,tpp,d,dm,b,tp(1:2)

call jd2greg(t,yr,mon,day)

write(uc(1,7),*)i-mcount,id,b*aukm,tp*aukm,t,d,tpp,dm, yr,mon,day,aukm

return
end subroutine
!********************************************************
subroutine mdistdmp
!output in close encounter file for target plane analysis
implicit none
integer::i,yr,mon,day
real(kind=dp),parameter::daysec=86400._dp

do i=mcount+1+m0cur-m0run,mcount+m0cur
call jd2greg(tdagminn(i),yr,mon,day)

write(uc(1,8),*)i-mcount,fbpcoo(i,1:3),dagminn(i),tdagminn(i),vcc(i)*aukm/daysec,dagmina(i),tdagmina(i),&
           vinf(i)*aukm/daysec,ncc(i),ndcc(i),yr,mon,day
end do


return
end subroutine
!*********************************************************
subroutine locccdmp
!output of stats for each close encounter
implicit none
integer::i,j,yr,mon,day
real(kind=dp),parameter::daysec=86400._dp
real(kind=dp)::maxndcc
real(kind=dp),dimension(6)::lmin,lmax,lmed,lmeddev

! do j=1,maxval(ncc(:))
!   do i=mcount+1+m0cur-m0run,mcount+m0cur
!     call jd2greg(loccc(i,j,1),yr,mon,day)
! write(uc(1,11),*)i-mcount,loccc(i,j,1:6),ncc(i),ndcc(i),yr,mon,day
!   end do
! end do
! 

maxndcc=maxval(ndcc(:)) 

do j=1,maxndcc
  do i=mcount+1,totbody
    call jd2greg(loccc(i,j,1),yr,mon,day)
    write(uc(1,11),*)j,i-mcount,loccc(i,j,1:6),ncc(i),ndcc(i),yr,mon,day
  end do
   write(uc(1,11),*)' '
end do

return
end subroutine
!*********************************************************

subroutine progress(now,total)
!--------------------------------------
! view the progress of a calculation in %
! now[real] ... current number of iterations 
! total [real]... total number of iterations
!
! dependencies: none
!---------------------------------------
real(kind=dp)::now,total
real(kind=dp)::dum
character::back

!backspace character in ASCII Code
back=achar(08)
!percentage of progress
dum=now/total*100.d0
!repositioning the cursor
write(unit=*,fmt='(A)',advance='no')back//back//back//back//back//back//back
!write the result and suppress the new line ($)
write(unit=*,fmt='(F6.2,A)',advance='no')dum,'%'

!final carriage return
if (now.ge.total) then
   write(unit=*,fmt='(A)')' '
end if
end subroutine progress


end module output_m


