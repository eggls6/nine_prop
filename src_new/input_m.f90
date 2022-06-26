module input_m
  use global_m
 ! use yark_m
  implicit none

  public::ReadConfig
  public::ReadIniRV
  public::ReadIniEL
  public::ReadIniYEL

 contains

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!                                                          INPUT 
!
!                                         by Siegfried Eggl  20080110 
!--------------------------------------------------------------------------------------------------------------------
!***************************************************************************************************************  
  subroutine ReadConfig(inm,inm0)
      implicit none
      character(len=50),intent(out)::inm,inm0
      integer::i
      logical::infile
      !Which output files
      character(len=3),dimension(1:8)::output
                       
        inquire(file="config.inn",exist=infile)

        if(infile) then

        open(unit=21, file="config.inn", status="old", action="read")

        !start time, end time, output time intervall
        read(unit=21, fmt=*)tstart, tend, tout
        !putting the start time value into a global variable
       
        tcalc=tend-tstart
        if(tend<tstart) then
          write(*,*) 'integration start time > integration end time'
          STOP 
        end if   
        if(tend>tephend.or.tstart<tephstart) then
          write(*,*) 'integration time outside of ephemeris data'
          write(*,*)'t ephemeris start:',tephstart
           write(*,*)'t ephemeris end:',tephstart
          STOP
        end if  
          
                       !Anzahl der Koerper
       read(unit=21, fmt=*)mcount,m0count,m0run
      
       !minimum stepsize
       read(unit=21,fmt=*)dtmin
       !maximum stepsize
       read(unit=21,fmt=*)dtmax
                       !one step accuracy
       read(unit=21,fmt=*)eps
       !massive body input file name
       read(unit=21,fmt=*)inm
        !massless body input file name
       read(unit=21,fmt=*)inm0
                       !input format (keplerian elements, rv vectors)
       read(unit=21,fmt=*)input
                      !Cutoffradius (if = 0, no cuttoff), stop calcualtion when massive (1)/massless (2) particle is beyond cutoff (logical)
!                        read(unit=21,fmt=*)cutoff,cutoffstop(1:2)
!-------------------------------------------------------------------------------------------------------------------------
!                     !close encounter output?
!                       read(unit=21,fmt=*)cen,cend,cemlim 
!--------------------------------------------------------------------------------------------------------
           
!--------------------------------------------------------------------------------------------------------        
      !show progress in terminal?
      read(unit=21,fmt=*)showprog
!--------------------------------------------------------------------------------------------------------           
      !want massive bodies in outputfile as well?
      read(unit=21,fmt=*)massiveout
      !want ephemeris bodies in outputfile as well?
      read(unit=21,fmt=*)ephemout
!--------------------------------------------------------------------------------------------------------
       !limit distance to consider close encounter with Earth / maximum stepsize during close encounter
        read(unit=21,fmt=*)cclim,ccdtmax   
!--------------------------------------------------------------------------------------------------------    
       !distance to construct bplane (n x sphere of influence)
      read(unit=21,fmt=*)nsoi
!--------------------------------------------------------------------------------------
      !mitigation relevant variables
      read(unit=21,fmt=*)kick 
      !kick has the following values: 
      !'p'... kick at next pericenter passage after kicktime  
      !'t'... kick happens for all clones at kicktime irrespecitve of their position
      !'n'.... no kick 
      read(unit=21,fmt=*)kicktime,dkicktime !time for mitigation attept [JD], range of mitigation time [JD]
      read(unit=21,fmt=*)kickrnd  !should kick rnd be random (i.e. between min or max) or nominal (i.e. maximal)
!      read(unit=21,fmt=*)kickmax(1:3) !maximum size of kick [au/D]: (1) radial, (2) tangential, (3) out of plane (angular momentum direction)
!      read(unit=21,fmt=*)kickmin(1:3) !minimum size of kick [au/D]: (1) radial, (2) tangential, (3) out of plane (angular momentum direction)
      read(unit=21,fmt=*)betaf(1:3)  !min, nominal, max value of beta factor (momentum enhancement factor)
      read(unit=21,fmt=*)ejphi       !maximum angle between impact and ejecta vectors [deg]
      read(unit=21,fmt=*)vimph(1:3)  !relative impact velocity vector [km/s] ICRF
                         vimph=vimph/aukm*86400.d0 !change units to [au/D]
                         !call he2icrs(vimph,1)
      read(unit=21,fmt=*)msc(1:3)    !min, nominal, max value of S/C mass [kg]
      read(unit=21,fmt=*)mast(1:3)   !min, nominal, max value of asterid mass [kg] 
                         !no transformation into Msun necessary since the [kg] units will cancel when delta v is calculated
      read(unit=21,fmt=*)failureprob  !probability of mission success (1-failue probability) [0,1]
      read(unit=21,fmt=*)cct0     !start time for close encounter check (so as to avoid to get past close encounters)
!--------------------------------------------------------------------------
      !which output files do you want?
      read(unit=21,fmt=*)output(:)
                       
                       outcc=.false.
                       outegy=.false.
                       outhele=.false.
                       outbele=.false.
                       outhco=.false.
                       outbco=.false.
                       outdef=.false.
                       outbinhele=.false.
                       outbinhco=.false.
                       outbinbco=.false.
                     

                       do i=1,8
                          select case (output(i))
                          case ('egy','en')
                             outegy=.true.
                         case('hel')
                            outhele=.true.
                         case('bel')
                            outbele=.true.   
                         case('hco')
                            outhco=.true.
                          case('bco')
                             outbco=.true.
                          case('bbc')
                             outbinbco=.true.
                          case('bhe')
                             outbinhele=.true. 
                           case('bhc')
                             outbinhco=.true.   
                           case('cc','clo','enc','cen')
                             outcc=.true.
                         end select
                      end do

                      if (input.eq.'te') then      
                        outtroj=.true.
                      else
                        outtroj=.false.
                      end if

!no outputfiles given:
                      if (.not.outegy.and..not. outhele.and..not. outbele &
                           .and..not.outhco.and..not.outbco.and..not. &
                           outbinbco.and..not.outbinhco.and..not.outbinhele &
                           .and..not.outcc) then
                              write(unit=*,fmt=*)' '
                              write(unit=*,fmt=*) 'no specific outputfiles selected'
                              write(*,*)'stopped'
                              STOP
                       end if
                    
!generate standard config.inn file
                    else
write(unit=*,fmt='(A)')'No inputfile found (config.inn)!'
end if
 close(unit=21)
  return
  end subroutine ReadConfig
!*************************************************************
 subroutine ReadIniRV(um,um0,rv)
 !READ initial positions (rv(1:3)) and velocities (rv(4:6)) as well as masses (m) form input file 
        implicit none
        integer,intent(in)::um,um0
        real(kind=dp),intent(inout)::rv(:,:)
        integer(kind=idp)::j
        !----------global vars---------:
        ! inbody, m, mcount 
        !------------------------------
        
        !initial conditions
        !massive 
       if (mcount.ne.0) then    
        do j=1,mcount
           read(unit=um, fmt=*) rv(j,1:6),mass(j),initep(j)
           gmass(j)=mass(j)*kgc2
        end do
       end if
        
        !massless
        if (m0count.ne.0) then
         do j=mcount+1,totbody
                read(unit=um0, fmt=*) rv(j,1:6),initep(j)
         end do
        end if
        return
end subroutine ReadIniRV
!*************************************************************
 subroutine ReadIniEL(um,um0,ele)
  !READ initial Keplerian orbital elements as well as masses (m) form input file 
           implicit none   
           integer,intent(in)::um,um0     
           real(kind=dp),intent(inout)::ele(:,:)
           integer(kind=idp)::j                            
        !----------global vars---------:
        ! inbody, m, mcount 
        !------------------------------  
         
      !massive 
      if (mcount.ne.0) then  
      do j=1,mcount
       read(unit=um, fmt=*) ele(j,1:6),mass(j),initep(j)
       gmass(j)=mass(j)*kgc2
      end do
      end if
      
      !massless 
      if (m0count.ne.0) then
      do j=mcount+1,totbody
         read(unit=um0, fmt=*) ele(j,1:6),initep(j)
      end do
      end if
          return
 end subroutine ReadIniEL

!*************************************************************
 subroutine ReadIniYEL(um,um0,ele)
 !READ initial Keplerian orbital elements as well as masses (m) form input file 
   implicit none         
   integer,intent(in)::um,um0
   real(kind=dp),dimension(:,:),intent(out)::ele
   integer(kind=idp)::j
   character::dc
        !----------global vars---------:
        ! mcount, inbody, m, yark_p 
        !------------------------------ 
      
      !massive   
      !read data
        do j=1,mcount
                read(unit=um, fmt=*) ele(j,1:6),mass(j),initep(j),yark_p(j,2:10)
        end do
        
      !massless  
      !read data
        do j=mcount+1,totbody
                read(unit=um0, fmt=*) ele(j,1:6),initep(j),yark_p(j,2:10)
        end do
        
        
   do j=1,inbody
   ! call init_yark(yark_p(j,1:3))
    !given in km in the input file
    yark_p(j,7)=yark_p(j,7)*1000.d0     
    yark_p(j,1)=0.d0 
   end do         

   return
 end subroutine ReadIniYEL
!*************************************************************
!  subroutine ReadIniELT(ele,m,NrK)
!  use global_m
!    implicit none         
!    real(kind=dp),dimension(:,:),intent(out)::ele
!    real(kind=dp),dimension(:),intent(out)::m
!    integer(kind=idp),intent(in)::NrK
!    integer(kind=idp)::j
!  
!    
!    do j=1,NrK
!     read(unit=21,fmt=*)ele(j,1:6),mass(j),trojan(j)
!    end do
!                        
!    close(unit=21)
!    return
!  end subroutine ReadIniELT
 end module input_m
