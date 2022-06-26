program nine_prop
 
!----------------------------------------------------------------------------
! Asteroid propagation package written by Siegfried Eggl
! last update: 20121004
!---------------------------------------------------------------------------- 

use global_m
use getephem_m
use input_m
use transform_m
use output_m
use radau_m
use symp_m
!use yark_m

implicit none
!#########################################################################
!file opening and allocation
integer::astat,um,um0,hlines,hlines0
!dummy
integer(kind=idp)::i,j,nrun
!location and velocities of the bodies, initial rv, keplerian orbital elements
real(kind=dp),dimension(:,:),allocatable::rv,ele,rvinit,eleinit
!name of input files for massive and massless bodies
character(len=50)::inm,inm0
logical::ex
!local variables
real(kind=dp)::JDcomp,rnd
real(kind=dp)::t,tstartbu,tendbu

!debug
real(kind=dp),dimension(11,1:6)::rveph

!-------------------------------------------------------------------------
write(unit=*,fmt='(A)')'gathering of initial conditions...'

!get ephemeris version
call getdeversion

!get ephemeris parameters, e.g. gravitational constant, speed of light, masses, etc.
if(nint(DE4XX).ge.430) then
  call getp430
 else
  call getp405
end if



!read in config.inn file
call ReadConfig(inm,inm0)
!number of massless bodies per integration
 
 !m0run=10
!if(m0count.lt.10) then
! m0run=1
!end if

if(mod(m0count,m0run).ne.0) then
 write(*,*)'the number of massless bodies must be divisible by',m0run
 STOP
end if


 inbody=mcount+m0run
 totbody=mcount+m0count

 

!allocate things
allocate(mass(totbody),gmass(totbody), rvinit(totbody,6),eleinit(totbody,6), &
          kicked(inbody),kicktime1(inbody) ,rv(inbody,6), rvo(inbody,6),&    !,rvinit(inbody,6)
         initep(totbody),dadt_yark(totbody,2),trojan(totbody), ve2(neph), &
         ncc(totbody),ndcc(totbody),loccc(totbody,100,6), stat=astat)
 !if close encounter output then allocate logical array
  if(outcc) then
    allocate(ccenter(inbody),fbplane(totbody),fbpcoo(totbody,3),vinf(totbody),vcc(totbody))
    allocate(dagminn(totbody),dagmina(totbody),tdagminn(totbody),tdagmina(totbody))
  end if

 if(astat.eq.0) then
  else
   write(*,*)'Error: Allocation problem in main.f90'
   STOP
 end if
!-------------------------------------------------------------------------------
!zero local vars: rv,rvinit,um,um0,hlines,hlines0
!zero global vars: mass, dadt_yark, back
 call init(rv,rvo,um,um0,hlines,hlines0)
!------------------------------------------------------------------------------  
 !read in initial conditions
 
 !open massive input file
 call fopen(inm,um,hlines)
 !open massless input file
 call fopen(inm0,um0,hlines0)
 
 !calculate escape velocity for all ephemeris bodies (required for target plane analysis)
 !ve^2=2GM/r
 do i=1,neph
        ve2(i)=2.d0*kgc2*meph(i)/(rearth/aukm)
 end do
 

 
 !see which input is given
 select case (input)

  case('he','el')
       allocate(ele(max(totbody,neph),6),stat=astat)      
       call ReadIniEL(um,um0,eleinit)
       call htrnsko(eleinit,totbody,rvinit)

  case('ye')
          allocate(ele(max(totbody,neph),6),yark_p(totbody,10),stat=astat)
          call ReadIniYEL(um,um0,eleinit)
          call htrnsko(eleinit,totbody,rvinit)
                                              
  case('rv','ve','co')
       call ReadIniRV(um,um0,rvinit)        
       
!     case('te','tr')
!        allocate(ele(inbody,6),stat=astat)            
!        call ReadIniELT(um,um0,ele)                                            
!        call trojtrnsko(ele,rv)     
       
  case default
     write(unit=*,fmt=*)' '
     write(unit=*,fmt=*)'initial condition error: unknown type of input (elements or rv-coordinates)'
  end select
  

  if(astat.eq.0) then
  else
   write(*,*)'Error: Allocation problem in main.f90'
   STOP
  end if

 
  
 call fclose(um)
 call fclose(um0)  
  
  !propagate everything to tstart
  !for this purpose all massively interacting bodies have to be given at the same epoch!!!!
  if(abs(sum(initep(1:mcount))-real(mcount)*initep(1)).gt.1.d-1) then
          write(*,*)'Massive asteroid initial conditions need to be given at the same inital epoch!'
           write(*,*)'Correct and restart, please!'
          STOP
  end if
  
  !check whether required dates are within ephemeris range
      do i=1,totbody
         if(initep(i).lt.tephstart.or.initep(i).gt.tephend) then
         write(*,*)'Initial epoch of body'
              write(*,*)i,initep(i)
         write(*,*)'out of ephemeris range! Ephemeris t_min',tephstart,' Ephemeris t_max', tephend
         STOP
        end if
      end do
      
      
select case (kick)
 case ('t')
 write(*,*)'KICK at time', kicktime 
 case('p')
 write(*,*)'KICK at pericenter passage after', kicktime
 case default
end select 
 
!##################################################################################################
! HERE THE LOOP OVER TRAILS OF ASTEROIDS INSTEAD OF THE COMPLETE CLONE SET BEGINS. CLONES TO BE COMPUTED 
! ARE SELECTED AND PROPAGATED

   inbody=mcount+m0run
   mcp1=mcount+1
  
   nrun=nint(real(m0count)/real(m0run))
   m0cur=0
   
   loccc=1.d99
   
   
   write(*,*)'Radau-15 working...'
   call OpenOutput('Radau')
   
   !here the loop over the nrun integration packages a m0run particles begins
   do j=1,nrun
  
    rv(1:mcount,1:6)=rvinit(1:mcount,1:6)
    rv(mcp1:inbody,:)=rvinit(mcp1+m0cur:mcount+m0cur+m0run,:)
    fbplane(:)=.true.
    
    
    !initiate kick time uncertainty
    do i=mcp1,inbody
    call random_number(rnd)
    kicktime1(i)=(rnd-0.5_dp)*dkicktime+kicktime
    end do
    
   
 !   tstart=tstartbu
 !   tend=tendbu

!    !db
! write(*,*)'INITIAL CONDITIONS'
!   do i=1,mcount+m0run
!   write(*,*)rv(i,:)
!   end do
!   
!   
! write(*,*)'RVINIT'
!   do i=mcount+1+m0cur,mcount+m0cur+m0run
!   write(*,*)i,rvinit(i,:)
!   end do  
!  !edb
!  STOP
  !#####################################################################################################
  !HERE A PROPAGATION OF THE MASSIVE ASTEROIDS TO THE INITAL EPOCH OF THE MASSLESS SHOULD BE PERFORMED, 
  !IN ORDER TO PROPAGATE THEM TOGETHER TOWARDS TINIT.   
    
  !Preliminary propagation needed?
  if(abs(initep(mcount+1)-initep(1)).gt.1.d-12) then
  !propagate massive bodies only, adapt global number of bodies to propagate!    
  inbody=mcount
  
  write(*,*)'iniep NEO,massive',initep(mcount+1),initep(1)
  
  !specific for propagation of massive bodies with common start time
    if (initep(mcount+1).lt.initep(1)) then
      call bakoo(rv(1:mcount,1:6),mcount,initep(1))
      call init_Radau_bwd(initep(1),initep(mcount+1),rv(1:mcount,1:6)) 
     
      write(*,*)inbody,'massiv bodies propagated from,', max(initep(mcount+1),initep(1)), &
               'to JD ',min(initep(mcount+1),initep(1))           
    else
      call bakoo(rv(1:mcount,1:6),mcount,initep(mcount+1))
      call init_Radau_fwd(initep(1),initep(mcount+1),rv(1:mcount,1:6)) 
        
      write(*,*)inbody,'massiv bodies propagated from,', min(initep(mcount+1),initep(1)), &
               'to JD ',max(initep(mcount+1),initep(1))   
    end if
    
  else !NO PROPAGATION NECESSARY
   call bakoo(rv(1:mcount,1:6),mcount,initep(1))
  end if

  
  
 !########################################################### 
 !PROPAGATE EVERYTHING TO TSTART
  inbody=mcount+m0run
  
  if(abs(initep(mcount+1)-tstart).gt.1.d-12) then 
    !specifically for propagation of massive bodies with common start time
    if (tstart.lt.initep(mcount+1)) then
      call bakoo(rv(mcount+1:inbody,1:6),m0run,initep(mcount+1))
      
      !propagate everyting backwards to tstart 
      call init_Radau_bwd(initep(mcount+1),tstart,rv(1:inbody,1:6)) 
       write(*,*)inbody,'bodies propagated from JD',max(initep(mcount+1),tstart), &
               'to JD ',min(initep(mcount+1),tstart)     
    else
     call bakoo(rv(mcount+1:inbody,1:6),m0run,tstart) 
       !propagate everyting forwards to tstart 
     call init_Radau_fwd(initep(mcount+1),tstart,rv(1:inbody,1:6)) 
     !rectify direction of velocities
      write(*,*)inbody,'bodies propagated from JD',min(initep(mcount+1),tstart), &
            'to JD ',max(initep(mcount+1),tstart)    
    end if
   
  else !NO PROPAGATION NECESSARY
    call bakoo(rv(mcount+1:inbody,1:6),m0run,initep(mcount+1))           
  end if 
  !#########################################################
  !START OF ACTUAL PROPAGATION
 

   
  
  
! do a barycentric correction of all bodies
!   call bakoo(rv,inbody,initep(1))    
 
!start propagation      



   
!  
!   
  !debug
!    write(*,*)'****************** EPHEM **************'
! !    do i =1,neph
! !     write(*,*)i,meph(i)
! !    end do
!   call getrvhel(tstart,rveph)
!   do i =1,neph
!     write(*,*)i,rveph(i,:),meph(i)
!    end do
! !    
!   write(*,*)'****************** EPHEM **************'
!  call htrnsel(rveph,ele,meph,neph)
!   do i =1,neph
!     write(*,*)i,ele(i,:),meph(i)
!    end do
! 
! !   
!     write(*,*)'****************** HEL COO **************'
!    do i =1,inbody
!     write(*,*)i,rv(i,:),mass(i)
!    end do
! !   
!     write(*,*)'****************** HEL **************'
!   ! call hekoo(rv,tstart)
!    call htrnsel (rv,ele,mass,inbody)
!    do i =1,inbody
!     write(*,*)i,ele(i,:),mass(i)
!    end do
! !    
! ! 
!      write(*,*)'****************** BA COO **************'
!         call bakoo(rv,inbody,tstart)   
!     do i =1,inbody
!      write(*,*)i,rv(i,:),mass(i)
!     end do
! !    
! !     
!        write(*,*)'****************** HEL 2 **************'
!         call hekoo(rv,inbody,tstart)
!         call htrnsel (rv,ele,mass,inbody)
!     do i =1,inbody
!      write(*,*)i,ele(i,:),mass(i)
!     end do 
!       
!       write(*,*)'mcount',mcount
!       write(*,*)'m0count',m0count
!       write(*,*)'inbody',inbody
      
   !end debug
  
  ! do a barycentric correction of all bodies
 !  call bakoo(rv,inbody,tstart)    
  !start RADAU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

       call Out(0._dp,rv,mass)
       call Radau(rv) 
       
       m0cur=m0cur+m0run
       kicked(:)=.false.
       write(*,*)'iteration',j,'of',nrun,'completed'
       if (outcc) then
         !output all min distances for the current run
         call mdistdmp
       end if
       ephemout=.false.
       massiveout=.false.
   end do    
       
   call locccdmp    
       
       
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       
         !start Candy
!      write(*,*)'Candy 4th order symplectic working...'
!        call OpenOutput('Candy')
!        call Out(0._dp,rv,mass)
!        call Candy(rv) 
!        call CloseOutput
!        
       
       
!        start Yoshida
!      write(*,*)'Yoshida 8th order symplectic working...'
!        call OpenOutput('Yoshi')
!        call Out(0._dp,rv,mass)
!        call Yoshida(rv) 
!        if (outcc) then
!        call mdistout
!        end if
!        call CloseOutput
! !        
  !db
!   do i=1,inbody
!   write(*,*)'dagminn',i,dagminn(i),tdagminn(i)
!   end do     
!    do i=1,inbody
!   write(*,*)'dagmina',i,dagmina(i),tdagmina(i)
!   end do     
  
       
 end program    
 !****************************************************

  