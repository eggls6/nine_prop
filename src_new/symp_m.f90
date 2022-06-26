module symp_m
  use global_m
  use transform_m
  use output_m
  use getephem_m
  use closeenc_m
  
  implicit none

public::Candy
public::Yoshida
 contains
!***************************************************************** 
subroutine Candy(rv)
integer(kind=idp)::i
real(kind=dp),dimension(1:inbody,1:6)::rv,rv1
real(kind=dp),dimension(1:inbody,1:3)::acc
real(kind=dp),dimension(1:4)::coeffa,coeffb
real(kind=dp)::dt,t,tmod,ts

dt=0.0001_dp
t=0._dp
ts=tstart
tcalc=abs(tend-tstart)

!Coefficients (multiplied by timestep to save cpu resources)
 coeffa(1)=(2._dp+2._dp**(1._dp/3._dp)+2**(-1._dp/3._dp))/6._dp*dt
 coeffa(4)=coeffa(1)
 coeffa(2)=(1._dp-2._dp**(1._dp/3._dp)-2._dp**(-1._dp/3._dp))/6._dp*dt
 coeffa(3)=coeffa(2)

 coeffb(1)=0._dp
 coeffb(2)=1._dp/(2._dp-2._dp**(1._dp/3._dp))*dt
 coeffb(4)=coeffb(2)
 coeffb(3)=1._dp/(1._dp-2._dp**(2._dp/3._dp))*dt




prgs: do while(t.lt.tcalc)

! step1
                do i=1,inbody
                        rv1(i,1:3)=rv(i,1:3) + coeffa(1)* rv(i,4:6) 
                end do
       
                call accrho(ts+t+coeffa(1),rv1,acc) 
                                        
 !step 2      
          
                do i=1,inbody
                  rv1(i,4:6) = rv(i,4:6) + coeffb(2)* acc(i,1:3)
                  rv1(i,1:3) = rv1(i,1:3) +coeffa(2)* rv1(i,4:6)
                end do
              
                call accrho(ts+t+coeffa(2),rv1,acc)       

              
    !step 3        

               
                do i=1,inbody
                    rv1(i,4:6)= rv1(i,4:6) + coeffb(3) * acc(i,1:3)
                    rv1(i,1:3) = rv1(i,1:3) + coeffa(3) * rv1(i,4:6)
                end do
 
           
                call accrho(ts+t+coeffa(3),rv1,acc)       

 !step 4 
                !massive
                do i=1,inbody   
                        rv1(i,4:6) = rv1(i,4:6) + coeffb(4) * acc(i,1:3)
                        rv(i,1:3)= rv1(i,1:3) + coeffa(4) * rv1(i,4:6)     
                        rv(i,4:6) = rv1(i,4:6)
                end do

             

  
   t=(real(nint(t/dt,Kind=idp),kind=dp)+1._dp)*dt
 ! output

   tmod=modulo(t/tout,1._dp)
   if (tmod.lt.1.d-14 .or. abs(tmod-1._dp).lt.1.d-14.and.t.ge.tout)then
     call Out(t, rv,mass)    
   end if

   !Close Encounter check
   if(outcc) then
    call cccheck(t+tstart,dt,rv,idearth)
   end if
   
   rvo=rv
end do prgs
return
end subroutine 


!*************************************************************************
subroutine Yoshida(rv)
integer(kind=idp)::i,j
real(kind=dp),dimension(1:inbody,1:6)::rv,rv1
real(kind=dp),dimension(1:inbody,1:3)::acc
real(kind=dp),dimension(0:7)::w
real(kind=dp),dimension(1:16)::coeffc
real(kind=dp),dimension(1:15)::coeffd
real(kind=dp)::dt,t,tmod,ts

dt=1._dp
t=0._dp
ts=tstart
tcalc=abs(tend-tstart)

w(1)=-0.161582374150097E1
w(2)=-0.244699182370524E1
w(3)=-0.716989419708120E-2
w(4)=0.244002732616735E1
w(5)=0.157739928123617E0
w(6)= 0.182020630970714E1
w(7)=0.104242620869991E1

!Coefficients (multiplied by timestep to save cpu resources)

w(0)=1._dp-2._dp*sum(w(1:7))

do i=1,7
 coeffd(i)=w(8-i)
 coeffd(16-i)=coeffd(i)

 coeffc(i+1)=0.5_dp*(w(8-i)+w(7-i))
 coeffc(16-i)=coeffc(i+1)
end do

 coeffd(8)=w(0)
 coeffc(1)=0.5_dp*w(7)
 coeffc(16)=coeffc(1)
    
  !Coefficients (multiplied by timestep to save cpu resources)
        
   coeffc(:)=coeffc(:)*dt
   coeffd(:)=coeffd(:)*dt

prgs: do while(t.lt.tcalc)


! step1
                do i=1,inbody
                        rv1(i,1:3)=rv(i,1:3) + coeffc(1)*rv(i,4:6) 
                end do
    
                call accrho(ts+t+coeffc(1),rv1,acc) 
       
         
                 do i=1,inbody
                   rv1(i,4:6)=rv(i,4:6)
                 end do
                                   
 !steps 2-15      
do j=2,15
 !step j
             
                do i=1,inbody
                  rv1(i,4:6) = rv1(i,4:6) + coeffd(j-1)* acc(i,1:3)
                  rv1(i,1:3) = rv1(i,1:3) +coeffc(j)* rv1(i,4:6)
                end do
                
                 call accrho(ts+t+sum(coeffc(1:j)),rv1,acc) 
 end do

!step 16 

                do i=1,inbody   
                        rv(i,4:6) = rv1(i,4:6) + coeffd(15) * acc(i,1:3)
                        rv(i,1:3)= rv1(i,1:3) + coeffc(16) * rv(i,4:6)                            
                end do
          
  
   t=(real(nint(t/dt,Kind=idp),kind=dp)+1._dp)*dt
 ! output

   tmod=modulo(t/tout,1._dp)
   if (tmod.lt.1.d-14 .or. abs(tmod-1._dp).lt.1.d-14.and.t.ge.tout)then
     call Out(t, rv,mass)    
   end if

   !Close Encounter check
   if(outcc) then
    call cccheck(t+tstart,dt,rv,idearth)
   end if
   
   rvo=rv
end do prgs
return
end subroutine 
!******************************************************************************************
subroutine accrho(t,rv,acc) 
implicit none
integer(kind=idp)::i,j,l,k
real(kind=dp)::t !current time for ephemeris
real(kind=dp)::rv(1:inbody,1:6),acc(1:inbody,1:3)
real(kind=dp)::rij(1:inbody,1:mcount,1:3),rijeph(1:inbody,1:neph,1:3),rveph(1:neph,1:6)
real(kind=dp)::d(1:inbody,1:mcount),deph(1:inbody,1:neph), dv(1:3)

!------------- NEWTONIAN ACCELERATIONS -----------------------
call getrv(t,rveph)

!$omp parallel default(shared)   &
!$omp private(i,j,l)
!$omp do 
do i=1,inbody
acc(i,:)=0._dp
end do
!$omp end do nowait

!massive pairdistances

!$omp do
do i=1,mcount-1
   do j=i+1,mcount
      do l=1,3
         rij(i,j,l)=rv(j,l)-rv(i,l)
      end do
        rij(j,i,:)=-rij(i,j,:)
!norm of pair distances
        d(i,j)=Sqrt(Dot_Product(rij(i,j,:),rij(i,j,:)))
        d(j,i)=d(i,j)
    end do
end do
!$omp end do nowait

!$omp do
do i=1,mcount
   rij(i,i,:)=0._dp
end do
!$omp end do nowait

!massless distances to massive
!$omp do
do i=mcount+1,inbody
   do j=1,mcount
      do l=1,3
         rij(i,j,l)=rv(j,l)-rv(i,l)
      end do
!norm of pair distances
        d(i,j)=Sqrt(Dot_Product(rij(i,j,:),rij(i,j,:)))
    end do
end do
!$omp end do


!massive and massless distances to ephemeris
do i=1,inbody 
   do j=1,neph
       do l=1,3
         rijeph(i,j,l)=rveph(j,l)-rv(i,l)
       end do
         
!norm of pair distances
        deph(i,j)=Sqrt(Dot_Product(rijeph(i,j,:),rijeph(i,j,:)))
    end do
end do

!massive: accelerations

!$omp do
   do i=1,mcount
      do j=1,mcount!massive on massive
         if(i.eq.j) then
            else
          acc(i,:)=acc(i,:)+gmass(j)*rij(i,j,:)/d(i,j)**3._dp
         end if
      end do
 
      do j=1,neph !ephemeris on massive 
          acc(i,:)=acc(i,:)+gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp
      end do
   end do
!$omp end do nowait

!massive on massless
!$omp do
   do i=mcount+1,inbody
      !massive on massless
      do j=1,mcount
          acc(i,:)=acc(i,:)+gmass(j)*rij(i,j,:)/d(i,j)**3._dp
      end do
      
      !ephemeris on massless
      do j=1,neph
          acc(i,:)=acc(i,:)+gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp  
      end do
   end do  
!$omp end do
!$omp end parallel 


!-------------- RELATIVISTIC ACCELERATIONS DUE TO THE SUN -------------------
!write(*,*)'t -reph(earh),acc',t,rveph(idearth,:)

call accgr(acc,deph,rijeph,rv,rveph,idsun)

!-------------- RELATIVISTIC ACCELERATIONS DUE TO EARTH -------------------
!call accgr(acc,deph,rijeph,rv,rveph,idearth)

return
end subroutine accrho


!**************************************************************************
subroutine accgr(acc,d,rijeph,rv,rveph,id)
!calculates relativistic monopole accelerations due to body "id"
!following Beutler 2006
use global_m
implicit none
integer,intent(in)::id
real(kind=dp),intent(in)::d(:,:),rijeph(:,:,:),rveph(:,:),rv(:,:)
real(kind=dp),intent(inout)::acc(:,:)
integer(kind=idp)::i
integer::j
real(kind=dp)::drv(1:3)

!-------------- RELATIVISTIC ACCELERATIONS DUE TO SUN -------------------------
 do i=1,inbody
! dummy variable delta v sun particle
         do j=1,3
           drv(j)=rv(i,j+3)-rveph(id,j+3)
         end do
 acc(i,1:3)=acc(i,1:3)-cm2/d(i,id)**3*gmeph(id)*(rijeph(i,id,1:3)*(4._dp*gmeph(id)/d(i,id)-Dot_Product(drv(:),drv(:))) + &
             4._dp*Dot_Product(rijeph(i,id,1:3),drv(1:3))*drv(1:3))    
 !db
!           write(*,*)'acc,rel',-cm2/d(i,id)**3*gmeph(id)*(rijeph(i,id,1:3)*(4._dp*gmeph(id)/d(i,id)-Dot_Product(drv(:),drv(:))) + &
!              4._dp*Dot_Product(rijeph(i,id,1:3),drv(1:3))*drv(1:3))    
!              write(*,*)1.d0/sqrt(cm2),d(i,id),gmeph(id),drv
!               write(*,*)'A',cm2/d(i,id)**3*gmeph(id)
!               write(*,*)'B', -rijeph(i,id,1:3)*(4._dp*gmeph(id)/d(i,id)-Dot_Product(drv(:),drv(:)))
!               write(*,*)'C',4._dp*Dot_Product(-rijeph(i,id,1:3),drv(1:3))*drv(1:3)
!edb
 end do !i
return
end subroutine
!#################################################################
end module  
   