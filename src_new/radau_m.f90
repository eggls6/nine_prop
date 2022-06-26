 module radau_m
  use global_m
  use transform_m
  use output_m
  use getephem_m
  use closeenc_m
 

  implicit none

public::Radau  !Gauss Radau 15 (Everhart 1974, 1985)
public::init_radau_fwd !propagates different initial conditions to tstart
public::init_radau_bwd !propagates different initial conditions to tstart
public::kikickt !kinetic impactor kick at a given time 
!public::kikickp !kinetic impactor kick at pericenter after a given time
 contains

!***************************************************************************
subroutine Radau(rv)
implicit none
integer(kind=idp)::allocstat,i
integer(kind=idp)::cnt
integer(kind=idp)::stiffcnt
real(kind=dp),dimension(1:inbody,1:6)::rv,rvout
real(kind=dp)::t,tmod,tsub,ts,tpts
real(kind=dp)::dt,dtdmax,dtdmin,dtdid,dto,dtdido,sumdto
real(kind=dp)::c_bo(7,inbody,3),c_eo(7,inbody,3)


!initial stepsize
dt=1.d-8
!initialize time related vars
dtdid=0._dp
dtdmin=1.d11
dtdmax=0._dp
tsub=0._dp
t=0._dp
tkc=0._dp
tky=0._dp
tkt=0._dp
ts=tstart
tcalc=abs(tend-tstart)
!initialize counter vars
 cnt=1
stiffcnt=0
!initialize radau integrator constants
call radauconst
!initialize old state vector
 rvo(:,:)=rv(:,:)  
 ccenter=.false.

 
!--------------------- START INTEGRATION LOOP -----------------------------
prgs: do while(t.lt.tcalc)


 
      if(dt<dtmin) then !stepsize should not be smaller than minimum stepsize
         dt=dtmin     
      end if
      if(dt>dtmax) then !stepsize should not be larger than maximum stepsize
         dt=dtmax
      end if
      
      if(any(ccenter)) then !close encounter happening, reduce stepsize
         dt=min(ccdtmax,dt)
      end if


  !check for output / stepsize 
  !this part of the code is responsible for keeping the optimum stepsize of the actual integration 
  !while making substeps to ensure that the required output intervals are met correctly
   if(tout.ne.0._dp) then

      tmod=abs(tout-modulo(t,tout)) !calculate how much time is left before output    

      if(dt>tmod) then !Next step would be larger than permissible for output 
         dto=tmod      !make the step exactly as large as output
         rvout=rv      !copy rv and predictor values into dummy variable that will be propagated only for output reasons
         c_bo=c_b        
         c_eo=c_e
         
         !make some substeps if necessary
         
         if(dt.ge.(tout+tmod)) then !if the time step is much larger than the output interval make several substeps with constant step-size
          !sum of substeps (substep time)
          sumdto=0.d0
          !predictors
          dto=tmod
          
          do while(dt>sumdto) !sub-step loop
    
          !call intermediate steps with fixed stepsize (output interval)
          !fixed stepsize is permissible in this case since the true step would be much larger anyway
          !i.e. the local error condition is fulfilled anyway for all step sizes dto < dt
          call radau_fstep(tsub,rvout,dto,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
             
          sumdto=sumdto+dto      !keep track of substep time 
          tsub=t+sumdto
          
          if(tsub.ge.tcalc) then !substeps should not go beyond the total integration time 
             t=tcalc
             exit
          end if
          
          if(tsub.gt.t+dt) then  !do not duplicate the next integrator output
            exit
          end if
         
           !Output 
          call Out(tsub,rvout,mass)  
!             !Close Encounter check
!           if(outcc) then
!             call cccheck(tsub+tstart,dto,rvout,idearth)
!           end if
          
           !fixed stepsize since dto < dt permits it anyway
          dto=tout  
          end do
          
         else ! the timestep is only a little larger than the time to the next output
          dto=tmod
          rvout=rv
          c_bo=c_b
          c_eo=c_e
          
            !make one substep only 
           call radau_fstep(t,rvout,dto,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r) 
           !Output 
           call Out(t+dto,rvout,mass)
           
!           !Close Encounter check
!           if(outcc) then
!             call cccheck(tsub+tstart,dto,rvout,idearth)
!           end if
      end if           
    end if   
   end if 
   
  
 !backup old constants     
   c_bo=c_b
   c_eo=c_e
!INTEGRATION STEP
call radau_step(t,ts,rv,dt,dtdid,c_b,c_e,c_h,c_xc,c_vc,c_c,c_d,c_r)
!

!compensated summation for time
   tky=dtdid-tkc
   tkt=t+tky
   tkc=(tkt-t)-tky
   t=tkt !update
 
   if(dtdid.gt.dtdmax) then
      dtdmax=dtdid
    end if
    if (dtdid.le.dtdmin) then
       dtdmin=dtdid
       stiffcnt=stiffcnt+1
      else 
       stiffcnt=0 
   end if

   if (stiffcnt.gt.10000) then
     write(*,*)'possibly stiff region encountered, stopping'
     STOP
   end if 

   !Close Encounter check
   if(outcc) then
   tpts=t+tstart
     call cccheck(tpts,dtdid,rv,idearth)
   end if
 
   !db
  !write(*,*)dtdid,dtmax,tpts,cct0
   !edb
   
     !------------- Kick due to Kinetic impactor ----------------------------
 select case (kick)
 case ('t')
  call kikickt(t,dtdid,rvo,rv,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
 case ('p')
  !call kikickp(t,dtdid,rvo,rv,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
  write(*,*)'kick at pericenter deprecated'
  STOP
 case default
 
 end select
   
 rvo(:,:)=rv(:,:)  
 cnt=cnt+1
 

 
end do prgs

!deallocate radauconst constants
call dealloc

write(*,*)'time accomplished:',t,'with an average stepsize of:',t/dble(cnt)
write(*,*)'minimum stepsize',dtdmin
write(*,*)'maximum stepsize',dtdmax

return
end subroutine Radau
!**************************************************************************************

subroutine init_radau_fwd(ts,te,rv)
implicit none
!initial propagation: 

!t.... time to calculate (0...tcalc)
!rv... initial conditions of massive (and perhaps massless) bodies


integer(kind=idp)::allocstat
integer(kind=idp)::cnt
integer(kind=idp)::stiffcnt
real(kind=dp),dimension(1:inbody,1:6)::rv,rvout
real(kind=dp)::t,tmod,dtmin,tsub,ts,te
real(kind=dp)::dt,dtdmax,dtdmin,dtdid,dto,dtdido,sumdto
real(kind=dp)::c_bo(7,inbody,3),c_eo(7,inbody,3)

!initialize time related vars
dtmin=1.d-12
dt=1.d-4
dtdid=0._dp
dtdmin=1.d11
dtdmax=0._dp
t=0._dp
!initialize counter vars
 cnt=1
stiffcnt=0
!initialize radau integrator constants
call radauconst

!calculate integration time
tcalc=te-ts

Write(*,*)'tcalc,t,ts,te', tcalc,t,ts,te

!--------------------- START INTEGRATION LOOP -----------------------------
prgs: do while(t.lt.tcalc)
 
      if(t+dt.gt.tcalc) then !Next step would be larger than permissible for initial epoch
         dt=tcalc-t      !make the step exactly as large as output
      end if  
   
      if(dt<dtmin) then !stepsize should not be smaller than minimum stepsize
         dt=dtmin     
      end if
      if(dt>dtmax) then !stepsize should not be larger than maximum stepsize
         dt=dtmax
      end if
! 
!INTEGRATION STEP
call radau_step(t,ts,rv,dt,dtdid,c_b,c_e,c_h,c_xc,c_vc,c_c,c_d,c_r)
!
  
    t=t+dtdid !update timestep
  
   if(dtdid.gt.dtdmax) then
      dtdmax=dtdid
    end if
    if (dtdid.le.dtdmin) then
       dtdmin=dtdid
       stiffcnt=stiffcnt+1
      else 
       stiffcnt=0 
   end if

   if (stiffcnt.gt.10000) then
     write(*,*)'possibly stiff region encountered, stopping'
     STOP
   end if 
   
 cnt=cnt+1
end do prgs


write(*,*)'time accomplished:',t,'with an average stepsize of:',t/dble(cnt)
write(*,*)'minimum stepsize',dtdmin
write(*,*)'maximum stepsize',dtdmax

!deallocate radauconst constants
call dealloc

return
end subroutine
!******************************************************************************
subroutine init_radau_bwd(ts,te,rv)
implicit none
!initial propagation: 

!t.... time to calculate (0...tcalc)
!rv... initial conditions of massive (and perhaps massless) bodies

integer(kind=idp)::allocstat
integer(kind=idp)::cnt
integer(kind=idp)::stiffcnt
real(kind=dp),dimension(1:inbody,1:6)::rv,rvout
real(kind=dp)::t,tmod,dtmin,tsub,ts,te
real(kind=dp)::dt,dtdmax,dtdmin,dtdid,dto,dtdido,sumdto
real(kind=dp)::c_bo(7,inbody,3),c_eo(7,inbody,3)

!initialize time related vars
dtmin=1.d-12
dt=1.d-4
dtdid=0._dp
dtdmin=1.d11
dtdmax=0._dp
t=te
!initialize counter vars
 cnt=1
stiffcnt=0
!initialize radau integrator constants
call radauconst


!calculate integration time
tcalc=abs(ts-te)

!if start JD of massive bodies > JD of massless bodies, integrate backwards
!i.e. change velocities of massive bodies (only!) 
rv(1:inbody,4:6)=-rv(1:inbody,4:6)


!--------------------- START INTEGRATION LOOP -----------------------------
prgs: do while(t+te.gt.ts)
 
      if(t-dt.gt.ts) then !Next step would be larger than permissible for initial epoch
         dt=abs(ts-t)      !make the step exactly as large as output
      end if  
   
      if(dt<dtmin) then !stepsize should not be smaller than minimum stepsize
         dt=dtmin     
      end if
      if(dt>dtmax) then !stepsize should not be larger than maximum stepsize
         dt=dtmax
      end if
! 
!INTEGRATION STEP
call radau_step_bwd(t,te,rv,dt,dtdid,c_b,c_e,c_h,c_xc,c_vc,c_c,c_d,c_r)
!
    
    t=t-dtdid !update timestep
    
  
  
   if(dtdid.gt.dtdmax) then
      dtdmax=dtdid
    end if
    if (dtdid.le.dtdmin) then
       dtdmin=dtdid
       stiffcnt=stiffcnt+1
      else 
       stiffcnt=0 
   end if

   if (stiffcnt.gt.10000) then
     write(*,*)'possibly stiff region encountered, stopping'
     STOP
   end if 
   
 cnt=cnt+1
end do prgs


 !rectify direction of velocities
 rv(1:inbody,4:6)=-rv(1:inbody,4:6)

write(*,*)'time accomplished:',t,'with an average stepsize of:',t/dble(cnt)
write(*,*)'minimum stepsize',dtdmin
write(*,*)'maximum stepsize',dtdmax

!deallocate radauconst constants
call dealloc

return
end subroutine
!******************************************************************************
!
! Author: John E. Chambers
! Modified: Siegfried Eggl  201104
!
! Integrates bodies for one timestep DTDID using
! Everhart's RA15 integrator algorithm. The accelerations are calculated
! using the subroutine FROH. The accuracy of the step is approximately 
! determined by the tolerance parameter TOL.
!
! Based on RADAU by E. Everhart, Physics Department, University of Denver.
! Comments giving equation numbers refer to Everhart (1985) ``An Efficient
! Integrator that Uses Gauss-Radau Spacings'', in The Dynamics of Comets:
! Their Origin and Evolution, p185-202, eds. A. Carusi & G. B. Valsecchi,
! pub Reidel. (A listing of the original subroutine is also given in this 
! paper.)
! 
!------------------------------------------------------------------------------
      subroutine radau_step(t,ts,rv,dt,dtdid,b,e,h,xc,vc,c,d,r)
      implicit none
!
! Input/Output
      real(kind=dp):: t,ts,dt,dtdid,tol
      real(kind=dp):: rv(:,:)  
    
!
! Local
      integer(kind=idp):: j,k,n,i
      real(kind=dp),dimension(7,inbody,3),intent(inout)::b,e
      real(kind=dp),dimension(7,inbody,3)::g
      real(kind=dp),intent(in)::h(8),xc(8),vc(7),c(21),d(21),r(28)
      real(kind=dp)::s(9),q,q2,q3,q4,q5,q6,q7,temp,gk,gkmax
      real(kind=dp),dimension(inbody,3)::acc,acc1
      real(kind=dp)::rv1(inbody,6),temp3(3)
      logical::iter

      tol=eps


!
! Calculate accs at the start of the sequence
      call accrho(t+ts,rv,acc) 
!
! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
!
!massive

!$omp parallel default(NONE) &
!$omp shared(g,b,d) &
!$omp firstprivate(inbody) &
!$omp private(k)
!$omp do
      do  k= 1, inbody
        g(1,k,:) = b(7,k,:)*d(16) + b(6,k,:)*d(11) + b(5,k,:)*d(7) &
                 + b(4,k,:)*d(4)  + b(3,k,:)*d(2)  + b(2,k,:)*d(1)  + b(1,k,:)
        g(2,k,:) = b(7,k,:)*d(17) + b(6,k,:)*d(12) + b(5,k,:)*d(8) &
                 + b(4,k,:)*d(5)  + b(3,k,:)*d(3)  + b(2,k,:)
        g(3,k,:) = b(7,k,:)*d(18) + b(6,k,:)*d(13) + b(5,k,:)*d(9)  &
                 + b(4,k,:)*d(6)  + b(3,k,:)
        g(4,k,:) = b(7,k,:)*d(19) + b(6,k,:)*d(14) + b(5,k,:)*d(10) + b(4,k,:)
        g(5,k,:) = b(7,k,:)*d(20) + b(6,k,:)*d(15) + b(5,k,:)
        g(6,k,:) = b(7,k,:)*d(21) + b(6,k,:)
        g(7,k,:) = b(7,k,:)
     end do
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
! For each iteration (six for firsdtcall to subroutine, two otherwise)...
 
!   do n = 1, niter
    n=0
    iter=.true.
    do while(iter)
!    do while(n<6)
!    
    n=n+1    
! For each substep within a sequence...
        do j = 2, 8
!
!$omp parallel default(NONE) &
!$omp shared(h,rv1,rv,b,acc,tol,iter,j) &
!$omp firstprivate(inbody,dt) &
!$omp private(s,k)
! Calculate position predictors using Eqn. 9 of Everhart
          s(1) = dt* h(j)
          s(2) = s(1) * s(1) * .5_dp
          s(3) = s(2) * h(j) *  1._dp/3._dp
          s(4) = s(3) * h(j) * .5_dp
          s(5) = s(4) * h(j) * .6_dp
          s(6) = s(5) * h(j) *  2._dp/3._dp
          s(7) = s(6) * h(j) * 5._dp/7._dp
          s(8) = s(7) * h(j) * .75_dp
          s(9) = s(8) * h(j) *  7._dp/9._dp
!
!massive
!$omp do
          do k = 1, inbody           
            rv1(k,1:3) = s(9)*b(7,k,:) + s(8)*b(6,k,:) + s(7)*b(5,k,:) &
                       + s(6)*b(4,k,:) + s(5)*b(3,k,:) + s(4)*b(2,k,:) &
                       + s(3)*b(1,k,:) + s(2)*acc(k,:)  + s(1)*rv(k,4:6) + rv(k,1:3)
          end do
!$omp end do

! calculate velocity predictors, from Eqn. 10 of Everhart
            s(1) = dt* h(j)
            s(2) = s(1) * h(j) * .5_dp
            s(3) = s(2) * h(j) * 2._dp/3._dp
            s(4) = s(3) * h(j) * .75_dp
            s(5) = s(4) * h(j) * .8_dp
            s(6) = s(5) * h(j) * 5._dp/6._dp
            s(7) = s(6) * h(j) * 6._dp/7._dp
            s(8) = s(7) * h(j) * .875_dp
!massive
!$omp do
            do k = 1, inbody
              rv1(k,4:6) = s(8)*b(7,k,:) + s(7)*b(6,k,:) + s(6)*b(5,k,:) &
                         + s(5)*b(4,k,:) + s(4)*b(3,k,:) + s(3)*b(2,k,:)    &
                         + s(2)*b(1,k,:) + s(1)*acc(k,:)  + rv(k,4:6)
            end do
!$omp end do
!$omp end parallel

!
! Calculate accs at the current substep (ATTENTION: UPDATE EPHEMERIS TIME!!!)
      

          call accrho(ts+t+s(1),rv1,acc1) 
 
! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
 
  select case(j)
     case(2)
!$omp parallel default(NONE) &
!$omp shared(b,r,g,acc,acc1) &
!$omp firstprivate(inbody) &
!$omp private(k,temp3)
!massive
!$omp do
            do k = 1, inbody
              temp3 = g(1,k,:)
              g(1,k,:) = (acc1(k,:) - acc(k,:)) * r(1)
              b(1,k,:) = b(1,k,:) + g(1,k,:) - temp3
            end do
!$omp end do
!$omp end parallel

      case(3)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,inbody
              do i=1,3             
               temp = g(2,k,i)
               gk = acc1(k,i) - acc(k,i)
               g(2,k,i) = (gk*r(2) - g(1,k,i))*r(3)
               temp = g(2,k,i) - temp
               b(1,k,i) = b(1,k,i)  +  temp * c(1)
               b(2,k,i) = b(2,k,i)  +  temp
              end do
            end do
!$omp end do
!$omp end parallel
      case(4)
!$omp parallel default(NONE) &
!$omp shared(b,r,c,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,inbody
             do i=1,3
              temp = g(3,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(3,k,i) = ((gk*r(4) - g(1,k,i))*r(5) - g(2,k,i))*r(6)
              temp = g(3,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(2)
              b(2,k,i) = b(2,k,i)  +  temp * c(3)
              b(3,k,i) = b(3,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
       case(5)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(4,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(4,k,i) = (((gk*r(7) - g(1,k,i))*r(8) - g(2,k,i))*r(9) &
                    - g(3,k,i))*r(10)
              temp = g(4,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(4)
              b(2,k,i) = b(2,k,i)  +  temp * c(5)
              b(3,k,i) = b(3,k,i)  +  temp * c(6)
              b(4,k,i) = b(4,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
        case(6)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,j) &
!$omp firstprivate(inbody)&
!$omp private(i,k,gk,temp3,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(5,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(5,k,i) =  ((((gk*r(11) - g(1,k,i))*r(12) - g(2,k,i))*r(13) &
                    - g(3,k,i))*r(14) - g(4,k,i))*r(15)
              temp = g(5,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(7)
              b(2,k,i) = b(2,k,i)  +  temp * c(8)
              b(3,k,i) = b(3,k,i)  +  temp * c(9)
              b(4,k,i) = b(4,k,i)  +  temp * c(10)
              b(5,k,i) = b(5,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
         case(7)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,inbody,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(6,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(6,k,i) = (((((gk*r(16) - g(1,k,i))*r(17) - g(2,k,i))*r(18) &
                     - g(3,k,i))*r(19) - g(4,k,i))*r(20) - g(5,k,i))*r(21)
              temp = g(6,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(11)
              b(2,k,i) = b(2,k,i)  +  temp * c(12)
              b(3,k,i) = b(3,k,i)  +  temp * c(13)
              b(4,k,i) = b(4,k,i)  +  temp * c(14)
              b(5,k,i) = b(5,k,i)  +  temp * c(15)
              b(6,k,i) = b(6,k,i)  +  temp
            end do
          end do

!$omp end do
!$omp end parallel

         case(8)
           gkmax=0._dp
!massive
            do k = 1, inbody
             do i=1,3
              temp = g(7,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(7,k,i) = ((((((gk*r(22) - g(1,k,i))*r(23) - g(2,k,i))*r(24) &
                     - g(3,k,i))*r(25) - g(4,k,i))*r(26) - g(5,k,i))*r(27) &
                     - g(6,k,i))*r(28)
              temp = g(7,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(16)
              b(2,k,i) = b(2,k,i)  +  temp * c(17)
              b(3,k,i) = b(3,k,i)  +  temp * c(18)
              b(4,k,i) = b(4,k,i)  +  temp * c(19)
              b(5,k,i) = b(5,k,i)  +  temp * c(20)
              b(6,k,i) = b(6,k,i)  +  temp * c(21)
              b(7,k,i) = b(7,k,i)  +  temp
              gkmax=max(gkmax,abs(temp))             
             end do !i 1..3
            end do  !k 1..nbody
           
!more iterations necessary?
                if(gkmax.gt.tol) then
                  iter=.true.
               else
                  iter=.false.
               end if
!              write(*,*)gkmax,n
        end select
      end do      !j 2...8
!     rerun subroutine with smaller stepsize if things do not converge      
      if(n>6) then
         dt=dt/3.13_dp
         dtdid=0._dp
!         write(*,*)'rvint',rv(2,:)
       !  write(*,*)'no convergence, trying smaller initial stepsize!',t*3.13d0,'->',dt       
         return
      end if
     end do      !n... while(iter)
!
!------------------------------------------------------------------------------
!
!  END  OF  MAIN  LOOP
!
! Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16)
      temp = 0._dp

      do k = 1, inbody
       do i=1,3
        temp = max( temp, abs( b(7,k,i) ) )
       end do
      end do
!
      temp = temp / (72._dp * abs(dt)**7)
      dtdid= dt
    
      if (temp.eq.0.0) then
        dt= dtdid* 1.4_dp
      else
        dt= sign( (tol/temp)**(1._dp/9._dp), dtdid)
      end if
!
! If new sequence size is much bigger than the current one, reduce it
      if (abs(dt/dtdid).gt.1.4_dp) dt= dtdid* 1.4_dp
!
! Find new position and velocity values at end of the sequence (Eqs. 11, 12)

      temp = dtdid* dtdid

!$omp parallel default(NONE) &
!$omp shared(rv,xc,vc,b,acc,inbody,dt,dtdid,temp) &
!$omp private(k)
!massive
!$omp do
      do k = 1 , inbody
        rv(k,1:3) = (xc(8)*b(7,k,:) + xc(7)*b(6,k,:) + xc(6)*b(5,k,:) &
                  +  xc(5)*b(4,k,:) + xc(4)*b(3,k,:) + xc(3)*b(2,k,:) &
                  +  xc(2)*b(1,k,:) + xc(1)*acc(k,:))*temp + rv(k,4:6)*dtdid+ rv(k,1:3)
!
        rv(k,4:6) = (vc(7)*b(7,k,:) + vc(6)*b(6,k,:) + vc(5)*b(5,k,:) &
                  +  vc(4)*b(4,k,:) + vc(3)*b(3,k,:) + vc(2)*b(2,k,:)     &
                  +  vc(1)*b(1,k,:) + acc(k,:))*dtdid+ rv(k,4:6)
      end do

!$omp end do
!$omp end parallel

! Predicdtnew B values to use at the start of the next sequence. The predicted
! values from the last call are saved as E. The correction, BD, between the
! actual and predicted values of B is applied in advance as a correction.
      q = dt/ dtdid
      q2 = q  * q
      q3 = q  * q2
      q4 = q2 * q2
      q5 = q2 * q3
      q6 = q3 * q3
      q7 = q3 * q4
!
!massive
      do k = 1, inbody
       do i =1,3
        s(1) = b(1,k,i) - e(1,k,i)
        s(2) = b(2,k,i) - e(2,k,i)
        s(3) = b(3,k,i) - e(3,k,i)
        s(4) = b(4,k,i) - e(4,k,i)
        s(5) = b(5,k,i) - e(5,k,i)
        s(6) = b(6,k,i) - e(6,k,i)
        s(7) = b(7,k,i) - e(7,k,i)
!
! Estimate B values for the nexdtsequence (Eqs. 13 of Everhart).
        e(1,k,i) = q* (b(7,k,i)* 7._dp + b(6,k,i)* 6._dp + b(5,k,i)* 5._dp &
               +       b(4,k,i)* 4._dp + b(3,k,i)* 3._dp + b(2,k,i)*2._dp + b(1,k,i))
        e(2,k,i) = q2*(b(7,k,i)*21._dp + b(6,k,i)*15._dp + b(5,k,i)*10._dp &
               +       b(4,k,i)* 6._dp + b(3,k,i)* 3._dp + b(2,k,i))      
        e(3,k,i) = q3*(b(7,k,i)*35._dp + b(6,k,i)*20._dp + b(5,k,i)*10._dp  &
               +       b(4,k,i)*4._dp  + b(3,k,i))
        e(4,k,i) = q4*(b(7,k,i)*35._dp + b(6,k,i)*15._dp + b(5,k,i)*5._dp + b(4,k,i))
        e(5,k,i) = q5*(b(7,k,i)*21._dp + b(6,k,i)*6._dp  + b(5,k,i))
        e(6,k,i) = q6*(b(7,k,i)*7._dp  + b(6,k,i))
        e(7,k,i) = q7* b(7,k,i)
!
        b(1,k,i) = e(1,k,i) + s(1)
        b(2,k,i) = e(2,k,i) + s(2)
        b(3,k,i) = e(3,k,i) + s(3)
        b(4,k,i) = e(4,k,i) + s(4)
        b(5,k,i) = e(5,k,i) + s(5)
        b(6,k,i) = e(6,k,i) + s(6)
        b(7,k,i) = e(7,k,i) + s(7)
       end do
      end do

      return
      end subroutine 
!*************************************************************************************      
     subroutine radau_fstep(t,rv,dt,b,e,h,xc,vc,c,d,r)
     !radau step with fixed stepsize dt
      implicit none
!
! Input/Output
      real(kind=dp):: t,dt
      real(kind=dp):: rv(:,:)  
!
! Local
      integer(kind=idp):: j,k,n,i
      real(kind=dp),dimension(7,inbody,3),intent(inout)::b,e
      real(kind=dp),dimension(7,inbody,3)::g
      real(kind=dp),intent(in)::h(8),xc(8),vc(7),c(21),d(21),r(28)
      real(kind=dp)::s(9),q,q2,q3,q4,q5,q6,q7,temp,gk,gkmax
      real(kind=dp),dimension(inbody,3)::acc,acc1
      real(kind=dp)::rv1(inbody,6),temp3(3)
      logical::iter

!
! Calculate accs at the start of the sequence
      call accrho(t+tstart,rv,acc) 
!
! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
!
!massive

!$omp parallel default(NONE) &
!$omp shared(g,b,d) &
!$omp firstprivate(inbody) &
!$omp private(k)
!$omp do
      do  k= 1, inbody
        g(1,k,:) = b(7,k,:)*d(16) + b(6,k,:)*d(11) + b(5,k,:)*d(7) &
                 + b(4,k,:)*d(4)  + b(3,k,:)*d(2)  + b(2,k,:)*d(1)  + b(1,k,:)
        g(2,k,:) = b(7,k,:)*d(17) + b(6,k,:)*d(12) + b(5,k,:)*d(8) &
                 + b(4,k,:)*d(5)  + b(3,k,:)*d(3)  + b(2,k,:)
        g(3,k,:) = b(7,k,:)*d(18) + b(6,k,:)*d(13) + b(5,k,:)*d(9)  &
                 + b(4,k,:)*d(6)  + b(3,k,:)
        g(4,k,:) = b(7,k,:)*d(19) + b(6,k,:)*d(14) + b(5,k,:)*d(10) + b(4,k,:)
        g(5,k,:) = b(7,k,:)*d(20) + b(6,k,:)*d(15) + b(5,k,:)
        g(6,k,:) = b(7,k,:)*d(21) + b(6,k,:)
        g(7,k,:) = b(7,k,:)
     end do
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
! For each iteration (six for firsdtcall to subroutine, two otherwise)...
 

   n=0
    do while(n<6)
!    
    n=n+1    
! For each substep within a sequence...
        do j = 2, 8
!
!$omp parallel default(NONE) &
!$omp shared(h,rv1,rv,b,acc,tol,iter,j) &
!$omp firstprivate(inbody,dt) &
!$omp private(s,k)
! Calculate position predictors using Eqn. 9 of Everhart
          s(1) = dt* h(j)
          s(2) = s(1) * s(1) * .5_dp
          s(3) = s(2) * h(j) *  1._dp/3._dp
          s(4) = s(3) * h(j) * .5_dp
          s(5) = s(4) * h(j) * .6_dp
          s(6) = s(5) * h(j) *  2._dp/3._dp
          s(7) = s(6) * h(j) * 5._dp/7._dp
          s(8) = s(7) * h(j) * .75_dp
          s(9) = s(8) * h(j) *  7._dp/9._dp
!
!massive
!$omp do
          do k = 1, inbody           
            rv1(k,1:3) = s(9)*b(7,k,:) + s(8)*b(6,k,:) + s(7)*b(5,k,:) &
                       + s(6)*b(4,k,:) + s(5)*b(3,k,:) + s(4)*b(2,k,:) &
                       + s(3)*b(1,k,:) + s(2)*acc(k,:)  + s(1)*rv(k,4:6) + rv(k,1:3)
          end do
!$omp end do

! calculate velocity predictors, from Eqn. 10 of Everhart
            s(1) = dt* h(j)
            s(2) = s(1) * h(j) * .5_dp
            s(3) = s(2) * h(j) * 2._dp/3._dp
            s(4) = s(3) * h(j) * .75_dp
            s(5) = s(4) * h(j) * .8_dp
            s(6) = s(5) * h(j) * 5._dp/6._dp
            s(7) = s(6) * h(j) * 6._dp/7._dp
            s(8) = s(7) * h(j) * .875_dp
!massive
!$omp do
            do k = 1, inbody
              rv1(k,4:6) = s(8)*b(7,k,:) + s(7)*b(6,k,:) + s(6)*b(5,k,:) &
                         + s(5)*b(4,k,:) + s(4)*b(3,k,:) + s(3)*b(2,k,:)    &
                         + s(2)*b(1,k,:) + s(1)*acc(k,:)  + rv(k,4:6)
            end do
!$omp end do
!$omp end parallel

!
! Calculate accs at the current substep (ATTENTION: UPDATE EPHEMERIS TIME!!!)
      
          call accrho(tstart+t+s(1),rv1,acc1) 
 
! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
 
  select case(j)
     case(2)
!$omp parallel default(NONE) &
!$omp shared(b,r,g,acc,acc1) &
!$omp firstprivate(inbody) &
!$omp private(k,temp3)
!massive
!$omp do
            do k = 1, inbody
              temp3 = g(1,k,:)
              g(1,k,:) = (acc1(k,:) - acc(k,:)) * r(1)
              b(1,k,:) = b(1,k,:) + g(1,k,:) - temp3
            end do
!$omp end do
!$omp end parallel

      case(3)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,inbody
              do i=1,3             
               temp = g(2,k,i)
               gk = acc1(k,i) - acc(k,i)
               g(2,k,i) = (gk*r(2) - g(1,k,i))*r(3)
               temp = g(2,k,i) - temp
               b(1,k,i) = b(1,k,i)  +  temp * c(1)
               b(2,k,i) = b(2,k,i)  +  temp
              end do
            end do
!$omp end do
!$omp end parallel
      case(4)
!$omp parallel default(NONE) &
!$omp shared(b,r,c,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,inbody
             do i=1,3
              temp = g(3,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(3,k,i) = ((gk*r(4) - g(1,k,i))*r(5) - g(2,k,i))*r(6)
              temp = g(3,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(2)
              b(2,k,i) = b(2,k,i)  +  temp * c(3)
              b(3,k,i) = b(3,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
       case(5)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(4,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(4,k,i) = (((gk*r(7) - g(1,k,i))*r(8) - g(2,k,i))*r(9) &
                    - g(3,k,i))*r(10)
              temp = g(4,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(4)
              b(2,k,i) = b(2,k,i)  +  temp * c(5)
              b(3,k,i) = b(3,k,i)  +  temp * c(6)
              b(4,k,i) = b(4,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
        case(6)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,j) &
!$omp firstprivate(inbody)&
!$omp private(i,k,gk,temp3,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(5,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(5,k,i) =  ((((gk*r(11) - g(1,k,i))*r(12) - g(2,k,i))*r(13) &
                    - g(3,k,i))*r(14) - g(4,k,i))*r(15)
              temp = g(5,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(7)
              b(2,k,i) = b(2,k,i)  +  temp * c(8)
              b(3,k,i) = b(3,k,i)  +  temp * c(9)
              b(4,k,i) = b(4,k,i)  +  temp * c(10)
              b(5,k,i) = b(5,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
         case(7)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,inbody,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(6,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(6,k,i) = (((((gk*r(16) - g(1,k,i))*r(17) - g(2,k,i))*r(18) &
                     - g(3,k,i))*r(19) - g(4,k,i))*r(20) - g(5,k,i))*r(21)
              temp = g(6,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(11)
              b(2,k,i) = b(2,k,i)  +  temp * c(12)
              b(3,k,i) = b(3,k,i)  +  temp * c(13)
              b(4,k,i) = b(4,k,i)  +  temp * c(14)
              b(5,k,i) = b(5,k,i)  +  temp * c(15)
              b(6,k,i) = b(6,k,i)  +  temp
            end do
          end do

!$omp end do
!$omp end parallel

         case(8)
           gkmax=0._dp
!massive
            do k = 1, inbody
             do i=1,3
              temp = g(7,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(7,k,i) = ((((((gk*r(22) - g(1,k,i))*r(23) - g(2,k,i))*r(24) &
                     - g(3,k,i))*r(25) - g(4,k,i))*r(26) - g(5,k,i))*r(27) &
                     - g(6,k,i))*r(28)
              temp = g(7,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(16)
              b(2,k,i) = b(2,k,i)  +  temp * c(17)
              b(3,k,i) = b(3,k,i)  +  temp * c(18)
              b(4,k,i) = b(4,k,i)  +  temp * c(19)
              b(5,k,i) = b(5,k,i)  +  temp * c(20)
              b(6,k,i) = b(6,k,i)  +  temp * c(21)
              b(7,k,i) = b(7,k,i)  +  temp
                      
             end do !i 1..3
            end do  !k 1..nbody
           

        end select
      end do      !j 2...8
     end do      !n... while(iter)
!
!------------------------------------------------------------------------------
!
!  END  OF  MAIN  LOOP
!
!
! Find new position and velocity values at end of the sequence (Eqs. 11, 12)

      temp = dt* dt

!$omp parallel default(NONE) &
!$omp shared(rv,xc,vc,b,acc,inbody,dt,dtdid,temp) &
!$omp private(k)
!massive
!$omp do
      do k = 1 , inbody
        rv(k,1:3) = (xc(8)*b(7,k,:) + xc(7)*b(6,k,:) + xc(6)*b(5,k,:) &
                  +  xc(5)*b(4,k,:) + xc(4)*b(3,k,:) + xc(3)*b(2,k,:) &
                  +  xc(2)*b(1,k,:) + xc(1)*acc(k,:))*temp + rv(k,4:6)*dt+ rv(k,1:3)
!
        rv(k,4:6) = (vc(7)*b(7,k,:) + vc(6)*b(6,k,:) + vc(5)*b(5,k,:) &
                  +  vc(4)*b(4,k,:) + vc(3)*b(3,k,:) + vc(2)*b(2,k,:)     &
                  +  vc(1)*b(1,k,:) + acc(k,:))*dt+ rv(k,4:6)
      end do

!$omp end do
!$omp end parallel


! Predicdtnew B values to use at the start of the next sequence. The predicted
! values from the last call are saved as E. The correction, BD, between the
! actual and predicted values of B is applied in advance as a correction.
!
!massive
      do k = 1, inbody
       do i =1,3
        s(1) = b(1,k,i) - e(1,k,i)
        s(2) = b(2,k,i) - e(2,k,i)
        s(3) = b(3,k,i) - e(3,k,i)
        s(4) = b(4,k,i) - e(4,k,i)
        s(5) = b(5,k,i) - e(5,k,i)
        s(6) = b(6,k,i) - e(6,k,i)
        s(7) = b(7,k,i) - e(7,k,i)
!
! Estimate B values for the nexdtsequence (Eqs. 13 of Everhart).
        e(1,k,i) = (b(7,k,i)* 7._dp + b(6,k,i)* 6._dp + b(5,k,i)* 5._dp &
               +       b(4,k,i)* 4._dp + b(3,k,i)* 3._dp + b(2,k,i)*2._dp + b(1,k,i))
        e(2,k,i) = (b(7,k,i)*21._dp + b(6,k,i)*15._dp + b(5,k,i)*10._dp &
               +       b(4,k,i)* 6._dp + b(3,k,i)* 3._dp + b(2,k,i))      
        e(3,k,i) = (b(7,k,i)*35._dp + b(6,k,i)*20._dp + b(5,k,i)*10._dp  &
               +       b(4,k,i)*4._dp  + b(3,k,i))
        e(4,k,i) = (b(7,k,i)*35._dp + b(6,k,i)*15._dp + b(5,k,i)*5._dp + b(4,k,i))
        e(5,k,i) = (b(7,k,i)*21._dp + b(6,k,i)*6._dp  + b(5,k,i))
        e(6,k,i) = (b(7,k,i)*7._dp  + b(6,k,i))
        e(7,k,i) =  b(7,k,i)
!
        b(1,k,i) = e(1,k,i) + s(1)
        b(2,k,i) = e(2,k,i) + s(2)
        b(3,k,i) = e(3,k,i) + s(3)
        b(4,k,i) = e(4,k,i) + s(4)
        b(5,k,i) = e(5,k,i) + s(5)
        b(6,k,i) = e(6,k,i) + s(6)
        b(7,k,i) = e(7,k,i) + s(7)
       end do
      end do

      return
      end subroutine 
      
      
!******************************************************************************************************      
  subroutine radau_fstep1(t,id,rv,dt,b,e,h,xc,vc,c,d,r)
     !radau step with fixed stepsize dt for individual particle (id)
      implicit none
!
! Input/Output
      integer(kind=idp)::id
      real(kind=dp):: t,dt
      real(kind=dp):: rv(:,:)  
!
! Local
      integer(kind=idp):: j,k,n,i,mcp1
      real(kind=dp),dimension(7,inbody,3),intent(inout)::b,e
      real(kind=dp),dimension(7,inbody,3)::g
      real(kind=dp),intent(in)::h(8),xc(8),vc(7),c(21),d(21),r(28)
      real(kind=dp)::s(9),q,q2,q3,q4,q5,q6,q7,temp,gk,gkmax
      real(kind=dp),dimension(inbody,3)::acc,acc1
      real(kind=dp)::rv1(inbody,6),temp3(3)
      real(kind=dp)::rvex(6),bex(7,1,3)
      logical::iter

!switch position of particle to be calculated with first particle after massive to speed up indexing
     mcp1=mcount+1
     
     rvex(1:6)=rv(id,1:6)
     rv(id,1:6)=rv(mcp1,1:6)
     rv(mcp1,1:6)=rvex(1:6)
     
     
     
     bex(1:7,1,1:3)=b(1:7,id,1:3)
     b(1:7,id,1:3)=b(1:7,mcp1,1:3)
     b(1:7,mcp1,1:3)=bex(1:7,1,1:3)

! Calculate accs at the start of the sequence
     call accrho1(t+tstart,mcp1,rv(1:mcp1,:),acc(1:mcp1,:)) 
   
!  
! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
!
!massive

!$omp parallel default(NONE) &
!$omp shared(g,b,d) &
!$omp firstprivate(mcp1) &
!$omp private(k)
!$omp do
      do  k= 1, mcp1
        g(1,k,:) = b(7,k,:)*d(16) + b(6,k,:)*d(11) + b(5,k,:)*d(7) &
                 + b(4,k,:)*d(4)  + b(3,k,:)*d(2)  + b(2,k,:)*d(1)  + b(1,k,:)
        g(2,k,:) = b(7,k,:)*d(17) + b(6,k,:)*d(12) + b(5,k,:)*d(8) &
                 + b(4,k,:)*d(5)  + b(3,k,:)*d(3)  + b(2,k,:)
        g(3,k,:) = b(7,k,:)*d(18) + b(6,k,:)*d(13) + b(5,k,:)*d(9)  &
                 + b(4,k,:)*d(6)  + b(3,k,:)
        g(4,k,:) = b(7,k,:)*d(19) + b(6,k,:)*d(14) + b(5,k,:)*d(10) + b(4,k,:)
        g(5,k,:) = b(7,k,:)*d(20) + b(6,k,:)*d(15) + b(5,k,:)
        g(6,k,:) = b(7,k,:)*d(21) + b(6,k,:)
        g(7,k,:) = b(7,k,:)
     end do
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
! For each iteration (six for firsdtcall to subroutine, two otherwise)...
 

   n=0
    do while(n<6)
!    
    n=n+1    
! For each substep within a sequence...
        do j = 2, 8
!
!$omp parallel default(NONE) &
!$omp shared(h,rv1,rv,b,acc,tol,iter,j) &
!$omp firstprivate(mcp1,dt) &
!$omp private(s,k)
! Calculate position predictors using Eqn. 9 of Everhart
          s(1) = dt* h(j)
          s(2) = s(1) * s(1) * .5_dp
          s(3) = s(2) * h(j) *  1._dp/3._dp
          s(4) = s(3) * h(j) * .5_dp
          s(5) = s(4) * h(j) * .6_dp
          s(6) = s(5) * h(j) *  2._dp/3._dp
          s(7) = s(6) * h(j) * 5._dp/7._dp
          s(8) = s(7) * h(j) * .75_dp
          s(9) = s(8) * h(j) *  7._dp/9._dp
!
!massive
!$omp do
          do k = 1, mcp1          
            rv1(k,1:3) = s(9)*b(7,k,:) + s(8)*b(6,k,:) + s(7)*b(5,k,:) &
                       + s(6)*b(4,k,:) + s(5)*b(3,k,:) + s(4)*b(2,k,:) &
                       + s(3)*b(1,k,:) + s(2)*acc(k,:)  + s(1)*rv(k,4:6) + rv(k,1:3)
          end do
!$omp end do

! calculate velocity predictors, from Eqn. 10 of Everhart
            s(1) = dt* h(j)
            s(2) = s(1) * h(j) * .5_dp
            s(3) = s(2) * h(j) * 2._dp/3._dp
            s(4) = s(3) * h(j) * .75_dp
            s(5) = s(4) * h(j) * .8_dp
            s(6) = s(5) * h(j) * 5._dp/6._dp
            s(7) = s(6) * h(j) * 6._dp/7._dp
            s(8) = s(7) * h(j) * .875_dp
!massive
!$omp do
            do k = 1, mcp1
              rv1(k,4:6) = s(8)*b(7,k,:) + s(7)*b(6,k,:) + s(6)*b(5,k,:) &
                         + s(5)*b(4,k,:) + s(4)*b(3,k,:) + s(3)*b(2,k,:)    &
                         + s(2)*b(1,k,:) + s(1)*acc(k,:)  + rv(k,4:6)
            end do
!$omp end do
!$omp end parallel

!
! Calculate accs at the current substep (ATTENTION: UPDATE EPHEMERIS TIME!!!)
      
          call accrho1(tstart+t+s(1),mcp1,rv1(1:mcp1,:),acc1(1:mcp1,:)) 
               
  
 
! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
 
  select case(j)
     case(2)
!$omp parallel default(NONE) &
!$omp shared(b,r,g,acc,acc1) &
!$omp firstprivate(mcp1) &
!$omp private(k,temp3)
!massive
!$omp do
            do k = 1, mcp1
              temp3 = g(1,k,:)
              g(1,k,:) = (acc1(k,:) - acc(k,:)) * r(1)
              b(1,k,:) = b(1,k,:) + g(1,k,:) - temp3
            end do
!$omp end do
!$omp end parallel

      case(3)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1) &
!$omp firstprivate(mcp1)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,mcp1
              do i=1,3             
               temp = g(2,k,i)
               gk = acc1(k,i) - acc(k,i)
               g(2,k,i) = (gk*r(2) - g(1,k,i))*r(3)
               temp = g(2,k,i) - temp
               b(1,k,i) = b(1,k,i)  +  temp * c(1)
               b(2,k,i) = b(2,k,i)  +  temp
              end do
            end do
!$omp end do
!$omp end parallel
      case(4)
!$omp parallel default(NONE) &
!$omp shared(b,r,c,g,acc,acc1) &
!$omp firstprivate(mcp1)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,mcp1
             do i=1,3
              temp = g(3,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(3,k,i) = ((gk*r(4) - g(1,k,i))*r(5) - g(2,k,i))*r(6)
              temp = g(3,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(2)
              b(2,k,i) = b(2,k,i)  +  temp * c(3)
              b(3,k,i) = b(3,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
       case(5)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1) &
!$omp firstprivate(mcp1)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, mcp1
             do i=1,3
              temp = g(4,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(4,k,i) = (((gk*r(7) - g(1,k,i))*r(8) - g(2,k,i))*r(9) &
                    - g(3,k,i))*r(10)
              temp = g(4,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(4)
              b(2,k,i) = b(2,k,i)  +  temp * c(5)
              b(3,k,i) = b(3,k,i)  +  temp * c(6)
              b(4,k,i) = b(4,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
        case(6)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,j) &
!$omp firstprivate(mcp1)&
!$omp private(i,k,gk,temp3,temp)
!massive
!$omp do
            do k = 1, mcp1
             do i=1,3
              temp = g(5,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(5,k,i) =  ((((gk*r(11) - g(1,k,i))*r(12) - g(2,k,i))*r(13) &
                    - g(3,k,i))*r(14) - g(4,k,i))*r(15)
              temp = g(5,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(7)
              b(2,k,i) = b(2,k,i)  +  temp * c(8)
              b(3,k,i) = b(3,k,i)  +  temp * c(9)
              b(4,k,i) = b(4,k,i)  +  temp * c(10)
              b(5,k,i) = b(5,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
         case(7)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,mcp1,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, mcp1
             do i=1,3
              temp = g(6,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(6,k,i) = (((((gk*r(16) - g(1,k,i))*r(17) - g(2,k,i))*r(18) &
                     - g(3,k,i))*r(19) - g(4,k,i))*r(20) - g(5,k,i))*r(21)
              temp = g(6,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(11)
              b(2,k,i) = b(2,k,i)  +  temp * c(12)
              b(3,k,i) = b(3,k,i)  +  temp * c(13)
              b(4,k,i) = b(4,k,i)  +  temp * c(14)
              b(5,k,i) = b(5,k,i)  +  temp * c(15)
              b(6,k,i) = b(6,k,i)  +  temp
            end do
          end do

!$omp end do
!$omp end parallel

         case(8)
           gkmax=0._dp
!massive
            do k = 1, mcp1
             do i=1,3
              temp = g(7,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(7,k,i) = ((((((gk*r(22) - g(1,k,i))*r(23) - g(2,k,i))*r(24) &
                     - g(3,k,i))*r(25) - g(4,k,i))*r(26) - g(5,k,i))*r(27) &
                     - g(6,k,i))*r(28)
              temp = g(7,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(16)
              b(2,k,i) = b(2,k,i)  +  temp * c(17)
              b(3,k,i) = b(3,k,i)  +  temp * c(18)
              b(4,k,i) = b(4,k,i)  +  temp * c(19)
              b(5,k,i) = b(5,k,i)  +  temp * c(20)
              b(6,k,i) = b(6,k,i)  +  temp * c(21)
              b(7,k,i) = b(7,k,i)  +  temp
                      
             end do !i 1..3
            end do  !k 1..nbody
           

        end select
      end do      !j 2...8
     end do      !n... while(iter)
!
!------------------------------------------------------------------------------
!
!  END  OF  MAIN  LOOP
!
!
! Find new position and velocity values at end of the sequence (Eqs. 11, 12)

      temp = dt* dt

!$omp parallel default(NONE) &
!$omp shared(rv,xc,vc,b,acc,mcp1,dt,dtdid,temp) &
!$omp private(k)
!massive
!$omp do
      do k = 1 , mcp1
        rv(k,1:3) = (xc(8)*b(7,k,:) + xc(7)*b(6,k,:) + xc(6)*b(5,k,:) &
                  +  xc(5)*b(4,k,:) + xc(4)*b(3,k,:) + xc(3)*b(2,k,:) &
                  +  xc(2)*b(1,k,:) + xc(1)*acc(k,:))*temp + rv(k,4:6)*dt+ rv(k,1:3)
!
        rv(k,4:6) = (vc(7)*b(7,k,:) + vc(6)*b(6,k,:) + vc(5)*b(5,k,:) &
                  +  vc(4)*b(4,k,:) + vc(3)*b(3,k,:) + vc(2)*b(2,k,:)     &
                  +  vc(1)*b(1,k,:) + acc(k,:))*dt+ rv(k,4:6)
      end do

!$omp end do
!$omp end parallel


! Predicdtnew B values to use at the start of the next sequence. The predicted
! values from the last call are saved as E. The correction, BD, between the
! actual and predicted values of B is applied in advance as a correction.
!
!massive
      do k = 1, mcp1
       do i =1,3
        s(1) = b(1,k,i) - e(1,k,i)
        s(2) = b(2,k,i) - e(2,k,i)
        s(3) = b(3,k,i) - e(3,k,i)
        s(4) = b(4,k,i) - e(4,k,i)
        s(5) = b(5,k,i) - e(5,k,i)
        s(6) = b(6,k,i) - e(6,k,i)
        s(7) = b(7,k,i) - e(7,k,i)
!
! Estimate B values for the nexdtsequence (Eqs. 13 of Everhart).
        e(1,k,i) = (b(7,k,i)* 7._dp + b(6,k,i)* 6._dp + b(5,k,i)* 5._dp &
               +       b(4,k,i)* 4._dp + b(3,k,i)* 3._dp + b(2,k,i)*2._dp + b(1,k,i))
        e(2,k,i) = (b(7,k,i)*21._dp + b(6,k,i)*15._dp + b(5,k,i)*10._dp &
               +       b(4,k,i)* 6._dp + b(3,k,i)* 3._dp + b(2,k,i))      
        e(3,k,i) = (b(7,k,i)*35._dp + b(6,k,i)*20._dp + b(5,k,i)*10._dp  &
               +       b(4,k,i)*4._dp  + b(3,k,i))
        e(4,k,i) = (b(7,k,i)*35._dp + b(6,k,i)*15._dp + b(5,k,i)*5._dp + b(4,k,i))
        e(5,k,i) = (b(7,k,i)*21._dp + b(6,k,i)*6._dp  + b(5,k,i))
        e(6,k,i) = (b(7,k,i)*7._dp  + b(6,k,i))
        e(7,k,i) =  b(7,k,i)
!
        b(1,k,i) = e(1,k,i) + s(1)
        b(2,k,i) = e(2,k,i) + s(2)
        b(3,k,i) = e(3,k,i) + s(3)
        b(4,k,i) = e(4,k,i) + s(4)
        b(5,k,i) = e(5,k,i) + s(5)
        b(6,k,i) = e(6,k,i) + s(6)
        b(7,k,i) = e(7,k,i) + s(7)
       end do
      end do
      
      !switch stuff back
      
     
     rvex(1:6)=rv(id,1:6)
     rv(id,1:6)=rv(mcp1,1:6)
     rv(mcp1,1:6)=rvex(1:6)
     
     
     
     bex(1:7,1,1:3)=b(1:7,id,1:3)
     b(1:7,id,1:3)=b(1:7,mcp1,1:3)
     b(1:7,mcp1,1:3)=bex(1:7,1,1:3)
      

      return
      end subroutine       
!***********************************************************************      
     subroutine radau_step_bwd(t,te,rv,dt,dtdid,b,e,h,xc,vc,c,d,r)
   !  radau step back in time, i.e. forward in time with exchanged directions of velocities 
      implicit none
!
! Input/Output
      real(kind=dp):: t,dt,dtdid,tol
      real(kind=dp):: rv(:,:)  
    
!
! Local
      integer(kind=idp):: j,k,n,i
      real(kind=dp),dimension(7,inbody,3),intent(inout)::b,e
      real(kind=dp),dimension(7,inbody,3)::g
      real(kind=dp),intent(in)::h(8),xc(8),vc(7),c(21),d(21),r(28)
      real(kind=dp)::s(9),q,q2,q3,q4,q5,q6,q7,temp,gk,gkmax
      real(kind=dp),dimension(inbody,3)::acc,acc1
      real(kind=dp)::rv1(inbody,6),temp3(3),te
      logical::iter

tol=eps
!
! Calculate accs at the start of the sequence
      call accbwd(t+te,rv,acc) 
!
! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
!
!massive

!$omp parallel default(NONE) &
!$omp shared(g,b,d) &
!$omp firstprivate(inbody) &
!$omp private(k)
!$omp do
      do  k= 1, inbody
        g(1,k,:) = b(7,k,:)*d(16) + b(6,k,:)*d(11) + b(5,k,:)*d(7) &
                 + b(4,k,:)*d(4)  + b(3,k,:)*d(2)  + b(2,k,:)*d(1)  + b(1,k,:)
        g(2,k,:) = b(7,k,:)*d(17) + b(6,k,:)*d(12) + b(5,k,:)*d(8) &
                 + b(4,k,:)*d(5)  + b(3,k,:)*d(3)  + b(2,k,:)
        g(3,k,:) = b(7,k,:)*d(18) + b(6,k,:)*d(13) + b(5,k,:)*d(9)  &
                 + b(4,k,:)*d(6)  + b(3,k,:)
        g(4,k,:) = b(7,k,:)*d(19) + b(6,k,:)*d(14) + b(5,k,:)*d(10) + b(4,k,:)
        g(5,k,:) = b(7,k,:)*d(20) + b(6,k,:)*d(15) + b(5,k,:)
        g(6,k,:) = b(7,k,:)*d(21) + b(6,k,:)
        g(7,k,:) = b(7,k,:)
     end do
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
! For each iteration (six for firsdtcall to subroutine, two otherwise)...
 
!   do n = 1, niter
    n=0
    iter=.true.
    do while(iter)
!    do while(n<6)
!    
    n=n+1    
! For each substep within a sequence...
        do j = 2, 8
!
!$omp parallel default(NONE) &
!$omp shared(h,rv1,rv,b,acc,tol,iter,j) &
!$omp firstprivate(inbody,dt) &
!$omp private(s,k)
! Calculate position predictors using Eqn. 9 of Everhart
          s(1) = dt* h(j)
          s(2) = s(1) * s(1) * .5_dp
          s(3) = s(2) * h(j) *  1._dp/3._dp
          s(4) = s(3) * h(j) * .5_dp
          s(5) = s(4) * h(j) * .6_dp
          s(6) = s(5) * h(j) *  2._dp/3._dp
          s(7) = s(6) * h(j) * 5._dp/7._dp
          s(8) = s(7) * h(j) * .75_dp
          s(9) = s(8) * h(j) *  7._dp/9._dp
!
!massive
!$omp do
          do k = 1, inbody           
            rv1(k,1:3) = s(9)*b(7,k,:) + s(8)*b(6,k,:) + s(7)*b(5,k,:) &
                       + s(6)*b(4,k,:) + s(5)*b(3,k,:) + s(4)*b(2,k,:) &
                       + s(3)*b(1,k,:) + s(2)*acc(k,:)  + s(1)*rv(k,4:6) + rv(k,1:3)
          end do
!$omp end do

! calculate velocity predictors, from Eqn. 10 of Everhart
            s(1) = dt* h(j)
            s(2) = s(1) * h(j) * .5_dp
            s(3) = s(2) * h(j) * 2._dp/3._dp
            s(4) = s(3) * h(j) * .75_dp
            s(5) = s(4) * h(j) * .8_dp
            s(6) = s(5) * h(j) * 5._dp/6._dp
            s(7) = s(6) * h(j) * 6._dp/7._dp
            s(8) = s(7) * h(j) * .875_dp
!massive
!$omp do
            do k = 1, inbody
              rv1(k,4:6) = s(8)*b(7,k,:) + s(7)*b(6,k,:) + s(6)*b(5,k,:) &
                         + s(5)*b(4,k,:) + s(4)*b(3,k,:) + s(3)*b(2,k,:)    &
                         + s(2)*b(1,k,:) + s(1)*acc(k,:)  + rv(k,4:6)
            end do
!$omp end do
!$omp end parallel

!
! Calculate accs at the current substep (ATTENTION: UPDATE EPHEMERIS TIME!!!)
      

          call accbwd(te+t-s(1),rv1,acc1) 
 
! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
 
  select case(j)
     case(2)
!$omp parallel default(NONE) &
!$omp shared(b,r,g,acc,acc1) &
!$omp firstprivate(inbody) &
!$omp private(k,temp3)
!massive
!$omp do
            do k = 1, inbody
              temp3 = g(1,k,:)
              g(1,k,:) = (acc1(k,:) - acc(k,:)) * r(1)
              b(1,k,:) = b(1,k,:) + g(1,k,:) - temp3
            end do
!$omp end do
!$omp end parallel

      case(3)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,inbody
              do i=1,3             
               temp = g(2,k,i)
               gk = acc1(k,i) - acc(k,i)
               g(2,k,i) = (gk*r(2) - g(1,k,i))*r(3)
               temp = g(2,k,i) - temp
               b(1,k,i) = b(1,k,i)  +  temp * c(1)
               b(2,k,i) = b(2,k,i)  +  temp
              end do
            end do
!$omp end do
!$omp end parallel
      case(4)
!$omp parallel default(NONE) &
!$omp shared(b,r,c,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,inbody
             do i=1,3
              temp = g(3,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(3,k,i) = ((gk*r(4) - g(1,k,i))*r(5) - g(2,k,i))*r(6)
              temp = g(3,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(2)
              b(2,k,i) = b(2,k,i)  +  temp * c(3)
              b(3,k,i) = b(3,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
       case(5)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1) &
!$omp firstprivate(inbody)&
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(4,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(4,k,i) = (((gk*r(7) - g(1,k,i))*r(8) - g(2,k,i))*r(9) &
                    - g(3,k,i))*r(10)
              temp = g(4,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(4)
              b(2,k,i) = b(2,k,i)  +  temp * c(5)
              b(3,k,i) = b(3,k,i)  +  temp * c(6)
              b(4,k,i) = b(4,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
        case(6)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,j) &
!$omp firstprivate(inbody)&
!$omp private(i,k,gk,temp3,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(5,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(5,k,i) =  ((((gk*r(11) - g(1,k,i))*r(12) - g(2,k,i))*r(13) &
                    - g(3,k,i))*r(14) - g(4,k,i))*r(15)
              temp = g(5,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(7)
              b(2,k,i) = b(2,k,i)  +  temp * c(8)
              b(3,k,i) = b(3,k,i)  +  temp * c(9)
              b(4,k,i) = b(4,k,i)  +  temp * c(10)
              b(5,k,i) = b(5,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
         case(7)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,inbody,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, inbody
             do i=1,3
              temp = g(6,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(6,k,i) = (((((gk*r(16) - g(1,k,i))*r(17) - g(2,k,i))*r(18) &
                     - g(3,k,i))*r(19) - g(4,k,i))*r(20) - g(5,k,i))*r(21)
              temp = g(6,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(11)
              b(2,k,i) = b(2,k,i)  +  temp * c(12)
              b(3,k,i) = b(3,k,i)  +  temp * c(13)
              b(4,k,i) = b(4,k,i)  +  temp * c(14)
              b(5,k,i) = b(5,k,i)  +  temp * c(15)
              b(6,k,i) = b(6,k,i)  +  temp
            end do
          end do

!$omp end do
!$omp end parallel

         case(8)
           gkmax=0._dp
!massive
            do k = 1, inbody
             do i=1,3
              temp = g(7,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(7,k,i) = ((((((gk*r(22) - g(1,k,i))*r(23) - g(2,k,i))*r(24) &
                     - g(3,k,i))*r(25) - g(4,k,i))*r(26) - g(5,k,i))*r(27) &
                     - g(6,k,i))*r(28)
              temp = g(7,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(16)
              b(2,k,i) = b(2,k,i)  +  temp * c(17)
              b(3,k,i) = b(3,k,i)  +  temp * c(18)
              b(4,k,i) = b(4,k,i)  +  temp * c(19)
              b(5,k,i) = b(5,k,i)  +  temp * c(20)
              b(6,k,i) = b(6,k,i)  +  temp * c(21)
              b(7,k,i) = b(7,k,i)  +  temp
              gkmax=max(gkmax,abs(temp))             
             end do !i 1..3
            end do  !k 1..nbody
           
!more iterations necessary?
                if(gkmax.gt.tol) then
                  iter=.true.
               else
                  iter=.false.
               end if
!              write(*,*)gkmax,n
        end select
      end do      !j 2...8
!     rerun subroutine with smaller stepsize if things do not converge      
      if(n>6) then
         dt=dt/3.13_dp
         dtdid=0._dp
!         write(*,*)'rvint',rv(2,:)
       !  write(*,*)'no convergence, trying smaller initial stepsize!',t*3.13d0,'->',dt       
         return
      end if
     end do      !n... while(iter)
!
!------------------------------------------------------------------------------
!
!  END  OF  MAIN  LOOP
!
! Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16)
      temp = 0._dp

      do k = 1, inbody
       do i=1,3
        temp = max( temp, abs( b(7,k,i) ) )
       end do
      end do
!
      temp = temp / (72._dp * abs(dt)**7)
      dtdid= dt
    
      if (temp.eq.0.0) then
        dt= dtdid* 1.4_dp
      else
        dt= sign( (tol/temp)**(1._dp/9._dp), dtdid)
      end if
!
! If new sequence size is much bigger than the current one, reduce it
      if (abs(dt/dtdid).gt.1.4_dp) dt= dtdid* 1.4_dp
!
! Find new position and velocity values at end of the sequence (Eqs. 11, 12)

      temp = dtdid* dtdid

!$omp parallel default(NONE) &
!$omp shared(rv,xc,vc,b,acc,inbody,dt,dtdid,temp) &
!$omp private(k)
!massive
!$omp do
      do k = 1 , inbody
        rv(k,1:3) = (xc(8)*b(7,k,:) + xc(7)*b(6,k,:) + xc(6)*b(5,k,:) &
                  +  xc(5)*b(4,k,:) + xc(4)*b(3,k,:) + xc(3)*b(2,k,:) &
                  +  xc(2)*b(1,k,:) + xc(1)*acc(k,:))*temp + rv(k,4:6)*dtdid+ rv(k,1:3)
!
        rv(k,4:6) = (vc(7)*b(7,k,:) + vc(6)*b(6,k,:) + vc(5)*b(5,k,:) &
                  +  vc(4)*b(4,k,:) + vc(3)*b(3,k,:) + vc(2)*b(2,k,:)     &
                  +  vc(1)*b(1,k,:) + acc(k,:))*dtdid+ rv(k,4:6)
      end do

!$omp end do
!$omp end parallel

! Predicdtnew B values to use at the start of the next sequence. The predicted
! values from the last call are saved as E. The correction, BD, between the
! actual and predicted values of B is applied in advance as a correction.
      q = dt/ dtdid
      q2 = q  * q
      q3 = q  * q2
      q4 = q2 * q2
      q5 = q2 * q3
      q6 = q3 * q3
      q7 = q3 * q4
!
!massive
      do k = 1, inbody
       do i =1,3
        s(1) = b(1,k,i) - e(1,k,i)
        s(2) = b(2,k,i) - e(2,k,i)
        s(3) = b(3,k,i) - e(3,k,i)
        s(4) = b(4,k,i) - e(4,k,i)
        s(5) = b(5,k,i) - e(5,k,i)
        s(6) = b(6,k,i) - e(6,k,i)
        s(7) = b(7,k,i) - e(7,k,i)
!
! Estimate B values for the nexdtsequence (Eqs. 13 of Everhart).
        e(1,k,i) = q* (b(7,k,i)* 7._dp + b(6,k,i)* 6._dp + b(5,k,i)* 5._dp &
               +       b(4,k,i)* 4._dp + b(3,k,i)* 3._dp + b(2,k,i)*2._dp + b(1,k,i))
        e(2,k,i) = q2*(b(7,k,i)*21._dp + b(6,k,i)*15._dp + b(5,k,i)*10._dp &
               +       b(4,k,i)* 6._dp + b(3,k,i)* 3._dp + b(2,k,i))      
        e(3,k,i) = q3*(b(7,k,i)*35._dp + b(6,k,i)*20._dp + b(5,k,i)*10._dp  &
               +       b(4,k,i)*4._dp  + b(3,k,i))
        e(4,k,i) = q4*(b(7,k,i)*35._dp + b(6,k,i)*15._dp + b(5,k,i)*5._dp + b(4,k,i))
        e(5,k,i) = q5*(b(7,k,i)*21._dp + b(6,k,i)*6._dp  + b(5,k,i))
        e(6,k,i) = q6*(b(7,k,i)*7._dp  + b(6,k,i))
        e(7,k,i) = q7* b(7,k,i)
!
        b(1,k,i) = e(1,k,i) + s(1)
        b(2,k,i) = e(2,k,i) + s(2)
        b(3,k,i) = e(3,k,i) + s(3)
        b(4,k,i) = e(4,k,i) + s(4)
        b(5,k,i) = e(5,k,i) + s(5)
        b(6,k,i) = e(6,k,i) + s(6)
        b(7,k,i) = e(7,k,i) + s(7)
       end do
      end do

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
real(kind=dp)::acckt(3),acckc(3),accky(3) !compensated summation

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
      !kahan vars
    acckc(:)=0._dp
    acckt(:)=0._dp
    accky(:)=0._dp
      do j=1,mcount!massive on massive
         if(i.eq.j) then
            else
          !acc(i,:)=acc(i,:)+gmass(j)*rij(i,j,:)/d(i,j)**3._dp
            accky(:)=gmass(j)*rij(i,j,:)/d(i,j)**3._dp-acckc(:)
            acckt(:)=acc(i,:)+accky(:)
            acckc(:)=(acckt(:)-acc(i,:))-accky(:)
            acc(i,:)=acckt(:)
         end if
      end do
 
      do j=1,neph !ephemeris on massive 
          !acc(i,:)=acc(i,:)+gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp
            accky(:)=gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp-acckc(:)
            acckt(:)=acc(i,:)+accky(:)
            acckc(:)=(acckt(:)-acc(i,:))-accky(:)
            acc(i,:)=acckt(:)
      end do
   end do
!$omp end do nowait

!massive on massless
!$omp do
   do i=mcount+1,inbody
      
        !kahan vars
    acckc(:)=0._dp
    acckt(:)=0._dp
    accky(:)=0._dp
    
      !massive on massless  
      do j=1,mcount
          !acc(i,:)=acc(i,:)+gmass(j)*rij(i,j,:)/d(i,j)**3._dp
          accky(:)=gmass(j)*rij(i,j,:)/d(i,j)**3._dp-acckc(:)
          acckt(:)=acc(i,:)+accky(:)
          acckc(:)=(acckt(:)-acc(i,:))-accky(:)
          acc(i,:)=acckt(:)
      end do
      
      !ephemeris on massless
      do j=1,neph
          !acc(i,:)=acc(i,:)+gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp  
           accky(:)=gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp-acckc(:)
           acckt(:)=acc(i,:)+accky(:)
           acckc(:)=(acckt(:)-acc(i,:))-accky(:)
           acc(i,:)=acckt(:)
      end do
   end do  
!$omp end do
!$omp end parallel 


!-------------- RELATIVISTIC ACCELERATIONS DUE TO THE SUN -------------------
!write(*,*)'t -reph(earh),acc',t,rveph(idearth,:)

call accgr(acc,deph,inbody,rijeph,rv,rveph,idsun)

!-------------- RELATIVISTIC ACCELERATIONS DUE TO EARTH -------------------
!call accgr(acc,deph,rijeph,rv,rveph,idearth)

!-------------- ACCELERATIONS DUE TO GRAVITY TRACTOR -------------------
!call accgt(t,acc,rv)



return
end subroutine accrho
!********************************************************************************
subroutine accrho1(t,n,rv,acc) 
!calculates accelerations due to planets and massive bodies onto massive bodies and 1!!! testparticle
implicit none
integer(kind=idp)::i,j,l,k,n
real(kind=dp)::t !current time for ephemeris
real(kind=dp)::rv(1:n,1:6),acc(1:n,1:3)
real(kind=dp)::rij(1:n,1:mcount,1:3),rijeph(1:n,1:neph,1:3),rveph(1:neph,1:6)
real(kind=dp)::d(1:n,1:mcount),deph(1:n,1:neph), dv(1:3)
real(kind=dp)::acckt(3),acckc(3),accky(3) !compensated summation

!------------- NEWTONIAN ACCELERATIONS -----------------------
call getrv(t,rveph)

!$omp parallel default(shared)   &
!$omp private(i,j,l)
!$omp do 
do i=1,n
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
do i=mcount+1,n
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
do i=1,n 
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
   !kahan vars
    acckc(:)=0._dp
    acckt(:)=0._dp
    accky(:)=0._dp
      do j=1,mcount!massive on massive
         if(i.eq.j) then
            else
            !acc(i,:)=acc(i,:)+gmass(j)*rij(i,j,:)/d(i,j)**3._dp
            accky(:)=gmass(j)*rij(i,j,:)/d(i,j)**3._dp-acckc(:)
            acckt(:)=acc(i,:)+accky(:)
            acckc(:)=(acckt(:)-acc(i,:))-accky(:)
            acc(i,:)=acckt(:)
         end if
      end do
 
      do j=1,neph !ephemeris on massive 
          !acc(i,:)=acc(i,:)+gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp
            accky(:)=gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp-acckc(:)
            acckt(:)=acc(i,:)+accky(:)
            acckc(:)=(acckt(:)-acc(i,:))-accky(:)
            acc(i,:)=acckt(:)
      end do
   end do
!$omp end do nowait

!massive on massless
!$omp do
   do i=mcount+1,n
      !massive on massless
          acckc(:)=0._dp
          acckt(:)=0._dp
          accky(:)=0._dp
      do j=1,mcount
         ! acc(i,:)=acc(i,:)+gmass(j)*rij(i,j,:)/d(i,j)**3._dp
          accky(:)=gmass(j)*rij(i,j,:)/d(i,j)**3._dp-acckc(:)
          acckt(:)=acc(i,:)+accky(:)
          acckc(:)=(acckt(:)-acc(i,:))-accky(:)
          acc(i,:)=acckt(:)
      end do
      
      !ephemeris on massless
      do j=1,neph
          !acc(i,:)=acc(i,:)+gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp
           accky(:)=gmeph(j)*rijeph(i,j,:)/deph(i,j)**3._dp-acckc(:)
           acckt(:)=acc(i,:)+accky(:)
           acckc(:)=(acckt(:)-acc(i,:))-accky(:)
           acc(i,:)=acckt(:)
      end do
   end do  
!$omp end do
!$omp end parallel 


!-------------- RELATIVISTIC ACCELERATIONS DUE TO THE SUN -------------------
!write(*,*)'t -reph(earh),acc',t,rveph(idearth,:)

call accgr(acc,deph,n,rijeph,rv,rveph,idsun)

!-------------- RELATIVISTIC ACCELERATIONS DUE TO EARTH -------------------
!call accgr(acc,deph,rijeph,rv,rveph,idearth)

!-------------- ACCELERATIONS DUE TO GRAVITY TRACTOR -------------------
!call accgt(t,acc,rv)



return
end subroutine accrho1
!*********************************************************************
subroutine accbwd(t,rv,acc)
!accelerations for backward integration  
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
          acc(i,:)=acc(i,:)+mass(j)*rij(i,j,:)/d(i,j)**3._dp
         end if
      end do
 
      do j=1,neph !ephemeris on massive 
          acc(i,:)=acc(i,:)+meph(j)*rijeph(i,j,:)/deph(i,j)**3._dp
      end do
      acc(i,:)=acc(i,:)*kgc2
   end do
!$omp end do nowait

!massive on massless
!$omp do
   do i=mcount+1,inbody
      !massive on massless
      do j=1,mcount
          acc(i,:)=acc(i,:)+mass(j)*rij(i,j,:)/d(i,j)**3._dp
      end do
      
      !ephemeris on massless
      do j=1,neph
          acc(i,:)=acc(i,:)+meph(j)*rijeph(i,j,:)/deph(i,j)**3._dp
          
!           !db
!           write(*,*)t,i,j,meph(j),deph(i,j),rveph(j,1:3)
!           !edb
          
      end do
      acc(i,:)=acc(i,:)*kgc2
   end do
!$omp end do
!$omp end parallel 


!-------------- RELATIVISTIC ACCELERATIONS DUE TO THE SUN -------------------
!For backward integration (back... global logical) change sign of ephemeris velocities
rveph(idsun,4:6)=-rveph(idsun,4:6)
call accgr(acc,deph,inbody,rijeph,rv,rveph,idsun)


!-------------- RELATIVISTIC ACCELERATIONS DUE TO EARTH -------------------
!call accgr(acc,deph,rijeph,rv,rveph,idearth)

return
end subroutine accbwd


!*********************************************************************************************************
subroutine radauconst
!RA15 constants
implicit none
integer::allocstat,n,i,j

allocate(c_b(7,inbody,3),c_e(7,inbody,3),c_h(8),c_xc(8),c_vc(7),c_c(21),c_d(21),c_r(28), & 
         stat=allocstat)

  if (allocstat.ne.0) then
            write(*,*)'error in radau_m: allocation not possible, try reducing the number of particles'
  end if 


! Gauss-Radau spacings for substeps within a sequence, for the 15th order 
! integrator. The sum of the H values should be 3.733333333333333
!
  c_h(:)=(/          0._dp,.0562625605269221464656522_dp,.1802406917368923649875799_dp,  &
       .3526247171131696373739078_dp,.5471536263305553830014486_dp,.7342101772154105315232106_dp,  &
       .8853209468390957680903598_dp,.9775206135612875018911745_dp/)
!
! Constant coefficients used in series expansions for X and V
!  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
!  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8
 
!   c_xc(:)=(/.5_dp,.1666666666666667_dp,.08333333333333333_dp,.05_dp,  &
!        .03333333333333333_dp,.02380952380952381_dp,.01785714285714286_dp, &
!        .01388888888888889_dp/)
! 
!   c_vc(:)=(/.5d0,.3333333333333333_dp,.25_dp,.2_dp, &
!         .1666666666666667_dp,.1428571428571429_dp,.125_dp/)
  c_xc(1)=0.5_dp
  c_xc(2)=1._dp/6._dp
  c_xc(3)=1._dp/12._dp
  c_xc(4)=0.05_dp
  c_xc(5)=1._dp/30._dp
  c_xc(6)=1._dp/42._dp
  c_xc(7)=1._dp/56._dp
  c_xc(8)=1._dp/72._dp

 c_vc(1)=0.5_dp
 c_vc(2)=1._dp/3._dp
 c_vc(3)=0.25_dp
 c_vc(4)=0.2_dp
 c_vc(5)=1._dp/6._dp
 c_vc(6)=1._dp/7._dp
 c_vc(7)=0.125_dp

  c_b(:,:,:)=0._dp
  c_e(:,:,:)=0._dp
!
! set values of the constant arrays
! (R = R21, R31, R32, R41, R42, R43 in Everhart's paper.)
    
        n = 0
        do j = 2, 8
          do i = 1, j - 1
            n = n + 1
            c_r(n) = 1._dp / (c_h(j) - c_h(i))
          end do
        end do
!
! Constants to convert between B and G arrays (C = C21, C31, C32, C41, C42...)
        c_c(1) = - c_h(2)
        c_d(1) =   c_h(2)
        n = 1
        do j = 3, 7
          n = n + 1
          c_c(n) = -c_h(j) * c_c(n-j+2)
          c_d(n) =  c_h(2) * c_d(n-j+2)
          do i = 3, j - 1
            n = n + 1
            c_c(n) = c_c(n-j+1)  -  c_h(j) * c_c(n-j+2)
            c_d(n) = c_d(n-j+1)  +  c_h(i) * c_d(n-j+2)
          end do
          n = n + 1
          c_c(n) = c_c(n-j+1) - c_h(j)
          c_d(n) = c_d(n-j+1) + c_h(j)
        end do
return
end subroutine

!*************************************************************************
subroutine dealloc
implicit none

deallocate(c_b,c_e,c_h,c_xc,c_vc,c_c,c_d,c_r)
         
return
end subroutine

!**************************************************************************
subroutine accgr(acc,d,n,rijeph,rv,rveph,id)
!calculates relativistic monopole accelerations for n bodies due to body "id"
!following Beutler 2006
use global_m
implicit none
integer(kind=idp),intent(in)::n
integer,intent(in)::id
real(kind=dp),intent(in)::d(:,:),rijeph(:,:,:),rveph(:,:),rv(:,:)
real(kind=dp),intent(inout)::acc(:,:)
integer(kind=idp)::i
integer::j
real(kind=dp)::drv(1:3)

!-------------- RELATIVISTIC ACCELERATIONS DUE TO SUN -------------------------
 do i=1,n
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
!********************************************************************
subroutine kikickt(t,dt,rvo,rv,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
!kick at kicktime t
implicit none

real(kind=dp),dimension(7,inbody,3),intent(inout)::c_bo,c_eo
real(kind=dp),intent(in)::c_h(8),c_xc(8),c_vc(7),c_c(21),c_d(21),c_r(28)
real(kind=dp),intent(in)::t,dt
integer(kind=idp)::i,j
real(kind=dp)::dt1,dt2,pp,rnorm,rnd,tdum
real(kind=dp),dimension(1:inbody,1:6),intent(inout)::rv,rvo
real(kind=dp),dimension(1:3)::LRL,LRLn,L
real(kind=dp),dimension(1:3)::astt,asth,astr,astv
real(kind=dp),dimension(1:4)::rnd2
real(kind=dp),dimension(1:3)::impvp,ejdir,ejdirp,imp,impvrel,impvb,impvpol,dvast
real(kind=dp),dimension(1:3,1:3)::vmat
real(kind=dp)::kmax,kmin
real(kind=dp),dimension(1:3)::lkmin,lkmax
real(kind=dp),dimension(1:6)::rvsun
real(kind=dp)::mscrnd,mastrnd 
real(kind=dp)::mej !effective ejecta mass ~(beta-1)*||v_ast - v_S/C||s
real(kind=dp)::newphi !new angle between impact and ejecta vectors 
real(kind=dp)::aud2kms !conversion factor from au/D to km/s
!absolute values of deflection magnitudes
!kmin=sqrt(dot_product(kickmin,kickmin))
!kmax=sqrt(dot_product(kickmax,kickmax))

!if (kmin.gt.kmax) then
!write(*,*)'the minimum deflection in config.inn is smaller than the maximum deflection. Stopping.'
! STOP
!end if

aud2kms=aukm/86400._dp

tdum=t+tstart

if (tdum.gt.minval(kicktime1(mcp1:inbody))) then

 if(all(kicked)) then
 else
 
 do i=mcp1,inbody
 
  if(kicked(i)) then
  else
  
  !propagate partilce id to kicktime
  dt1=dt-(tdum-kicktime1(i))
  
  call radau_fstep1(t-dt,i,rvo,dt1,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
      
  !get rv vector of the sun at kicktime    
  call getrvsun(((dt1-dt)+t)+tstart,rvsun) 
      
!       !db 
!    call crossp3d(rv(i,1:3),rv(i,4:6),L)
! !   
!    call crossp3d(rv(i,4:6),L,LRL)
!    rnorm=Sqrt(Dot_Product(rv(i,1:3),rv(i,1:3)))
!    LRL(1:3)=LRL(1:3)/kgc2-rv(i,1:3)/rnorm
!    LRLn=LRL/dot_product(LRL,LRL)
!    write(*,*)LRLn(1:3),rv(i,1:3)/rnorm
!   if(abs(dot_product(LRLn(1:3),rv(i,1:3)/rnorm)-1.d0).lt.0.005d0) then
!edb
    
    
  !kick it
 !db!
  !rnd=0.d0
  !edb
 
 !check for complete failure
 call random_number(rnd)
 !rnd=rand()
 if(rnd.gt.failureprob) then
 vmat=0._dp
 else
 !actual mitigation 
   astr=rvo(i,1:3)/sqrt(dot_product(rvo(i,1:3),rvo(i,1:3)))
   astv=rvo(i,4:6)/sqrt(dot_product(rvo(i,4:6),rvo(i,4:6)))
   !calculate tangential, radial and out of plane directions 
  !GAUSS R S W SYSTEM
   !call crossp3d(astr,astv,asth)
   !asth=asth/sqrt(dot_product(asth,asth))
   !call crossp3d(asth,astr,astt)   

  !GAUSS T N W system 
!     call crossp3d(astr,astv,asth)
!    asth=asth/sqrt(dot_product(asth,asth))
!    astt=astv
!    call crossp3d(astt,asth,astr)
!    astr=astr/sqrt(dot_product(astr,astr))

   
   
    ! barycentric S/C velocity w.r.t barycentric NEO velocity to get relative impact vector
    do j=1,3
     impvb(j)=vimph(j)+rvo(i,3+j)!barycentric S/C velociy
     impvrel(j)=vimph(j)!-rvo(i,3+j) !S/C velocity relative to the NEO
    end do
   
   !should kick be a random process or nominal?
   if(kickrnd) then
   
    
!     if(kmax.le.kmin*1.d1.or.kmin.eq.0.d0) then !same order of magnitude -> uniform distribution of values
!       !define matrix that contains radial tangential and out of plane velocities
!        vmat(1,1:3)=astr(1:3)*rnd2(1)*(kickmax(1)-kickmin(1))+kickmin(1)
!        vmat(2,1:3)=astt(1:3)*rnd2(2)*(kickmax(2)-kickmin(2))+kickmin(2)
!        vmat(3,1:3)=asth(1:3)*rnd2(3)*(kickmax(3)-kickmin(3))+kickmin(3)
!     
!     else !use uniform logarithmic distribution
!           !define matrix that contains radial tangential and out of plane velocities
!       do j=1,3
!         if(kickmax(j).eq.0.d0) then
!          lkmax(j)=-98.d0
!         else
!          lkmax(j)=log10(abs(kickmax(j)))
!         end if       
!         if(kickmin(j).eq.0.d0) then
!          lkmin(j)=-99.d0
!         else
!          lkmin(j)=log10(abs(kickmin(j)))
!         end if
!       end do
!        
!        vmat(1,1:3)=astr(1:3)*kickmax(1)*10.d0**((lkmin(1)-lkmax(1))*rnd2(1))
!        vmat(2,1:3)=astt(1:3)*kickmax(2)*10.d0**((lkmin(2)-lkmax(2))*rnd2(2))
!        vmat(3,1:3)=asth(1:3)*kickmax(3)*10.d0**((lkmin(3)-lkmax(3))*rnd2(3))
!     end if !log distr


   newphi=380.d0 
   do while(newphi.gt.ejphi) !make sure the new ejecta vector falls within the cone spanned by 2*ejphi 
    call random_number(rnd2(1:2))
   
    !change to polar cooridnates to calculate ejecta direction variance 
    call topolar2(impvrel,impvpol)
    
    ejdirp(1)=1.d0
    ejdirp(2:3)=impvpol(2:3)+2._dp*(rnd2(1:2)-0.5_dp)*ejphi*grad
    
    call tocartesian2(ejdirp,ejdir)
   
    newphi=acos(dot_product(ejdir,impvrel)/&
               (sqrt(dot_product(impvrel,impvrel))))*180./pi
    !db write(*,*)'newphi,norms',newphi,sqrt(dot_product(ejdir,ejdir)),sqrt(dot_product(impvrel,impvrel))                    
   end do 
   
   call random_number(rnd2(1:4))
   !calculate the impact momentum from the mass of the sc
   mscrnd=((msc(3)-msc(1))*rnd2(3)+msc(1))
   imp=impvrel*mscrnd
   
   !calculate mass of asteroid
   mastrnd=((mast(3)-mast(1))*rnd2(3)+mast(1))
   
   
   !scale ejecta velocities via (beta-1)*|imp| and add them to the impact momentum to get the total transferred momentum 
!    dvast=(imp+(ejdir*sqrt(dot_product(imp,imp))*((betaf(3)-betaf(1))*rnd2(4)+betaf(1)-1._dp)))/(mastrnd+mscrnd)

  !make it simpler by using center of mass of asteroid, S/C and ejecta
  !all we need is to calculate the effective ejecta mass: m_E~(beta-1)*|v_ast-v_S/C|
  mej=((betaf(3)-betaf(1))*rnd2(4)+betaf(1)-1._dp)*sqrt(dot_product(impvrel,impvrel))
   
  astv(1:3)=rvo(i,4:6) !dbg
  !deflection change in asteroid's velocity (ma*va+msc*vsc-me*ve)/(ma+msc-me) 
  rvo(i,4:6)=(mastrnd*rvo(i,4:6)+(mscrnd*impvb(1:3)-mej*ejdir(1:3)))/(mastrnd+(mscrnd-mej))
    
   else !deterministic kick
   
   !dvast=betaf(2)*impvrel*msc(2)/(mast(2)+msc(2))
  
    ejdir=impvrel !ejecta and impact velocity are (anti) aligned --- sign comes later with the mass!
    mej=(betaf(2)-1._dp)
   
    astv(1:3)=rvo(i,4:6) !dbg 
  !deflection change in asteroid's velocity (ma*va+msc*vsc-me*ve)/(ma+msc-me) 
   rvo(i,4:6)=(mast(2)*rvo(i,4:6)+(msc(2)*impvb(1:3)-mej*ejdir(1:3)))/(mast(2)+(msc(2)-mej))
   end if !kickrnd
   
   !deflection change in asteroid's velocity
    dvast(1:3)=rvo(i,4:6)-astv(1:3)
 
!    write(*,*)'id,ast(6),kick(3),dvasc(3),pej(3),vsc,phi',i,((dt1-dt)+t)+tstart,rvo(i,1:3)*aukm, &
!                rvo(i,4:6)*aukm/86400.d0,dvast(1:3)*aukm/86400.d0, &
!                impvrel*aukm/86400.d0,mej*ejdir*aukm/86400.d0, &
!                (rvo(i,4:6)+impvrel(1:3))*aukm/86400.d0, &
!                acos(dot_product(mej*ejdir,impvrel)/&
!                (sqrt(dot_product(mej*ejdir,mej*ejdir)*dot_product(impvrel,impvrel))))*180./pi
!                
!   write(unit=uc(1,9),fmt=*)'id, time kick, ast rv(6), dV_a (3) [km/s], v_rel (3) [km/s],p_ej(3)[kg km/s],v_sc[km/S],phi[deg]', &
!                i+m0cur, ((dt1-dt)+t)+tstart,rvo(i,1:3)*aukm, &
!                rvo(i,4:6)*aukm/86400.d0,dvast(1:3)*aukm/86400.d0, &
!                impvrel*aukm/86400.d0,mej*ejdir*aukm/86400.d0, &
!                (rvo(i,4:6)+impvrel(1:3))*aukm/86400.d0, &
!                acos(dot_product(mej*ejdir,impvrel)/&
!                (sqrt(dot_product(mej*ejdir,mej*ejdir)*dot_product(impvrel,impvrel))))*180./pi

!calculate deflection signal to noise ratio directly


write(unit=uc(1,9),fmt=*)i,((dt1-dt)+t)+tstart,rvo(i,1:3)*aukm, &
               astv(1:3)*aud2kms,rvo(i,4:6)*aud2kms,dvast(1:3)*aud2kms, &
               impvrel*aud2kms,mej*ejdir*aud2kms, &
               (rvo(i,4:6)+impvrel(1:3))*aud2kms, &
               acos(dot_product(mej*ejdir,impvrel)/&
               (sqrt(dot_product(mej*ejdir,mej*ejdir)*dot_product(impvrel,impvrel))))*180.d0/pi
 end if                     
   
  
   
   
   kicked(i)=.true.

 
   !propagate particle id to end of original timestep
   dt2=dt-dt1
   call radau_fstep1(t-dt2,i,rvo,dt2,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
   
   !db
   !write(*,*)'*********  AFTER ******************'
   !write(*,*)rv(i,1:6)-rvo(i,1:6)
   !write(*,*)'**********************************'
   !edb
   rv(i,1:6)=rvo(i,1:6)
  
   !write(unit=uc(1,9),fmt=*)i+m0cur,t+dt1-dt,t+tstart+dt1-dt,&
   !sqrt(dot_product(vmat(1,:),vmat(1,:))),sqrt(dot_product(vmat(2,:),vmat(2,:))),&
   !sqrt(dot_product(vmat(3,:),vmat(3,:)))
   
   
   
 
 !db
 ! write(*,*)i+m0cur,'KICKED',tdum,kicktime,t+tstart-dt+dt1,dt,dt1,dt2,t,t+tstart-dt+dt1+dt2   
  !  write(*,*)'i,mcount',i,mcount,kicked              
  !edb
  end if
     
 end do
 end if !all kicked
end if !tmin
  
return
end subroutine
!*************************************************************************
! subroutine kikickp(t,dt,rvo,rv,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
! !kick at pericenter
! implicit none
! 
! real(kind=dp),dimension(7,inbody,3),intent(inout)::c_bo,c_eo
! real(kind=dp),intent(in)::c_h(8),c_xc(8),c_vc(7),c_c(21),c_d(21),c_r(28)
! real(kind=dp),intent(in)::t,dt
! integer(kind=idp)::i,j
! real(kind=dp)::dt1,dt2,pp,rnorm,rnd
! real(kind=dp),dimension(1:inbody,1:6)::rv,rvo
! real(kind=dp),dimension(1:3)::LRL,LRLn,L
! real(kind=dp),dimension(1:3)::astt,asth,astr,astv,rnd2,vscal
! real(kind=dp),dimension(1:6)::ele1
! real(kind=dp),dimension(1:3,1:3)::vmat
! 
! if (t+tstart.gt.kicktime) then
! 
!  if(all(kicked)) then
!  else
!  do i=mcp1,inbody
!   
!   if(kicked(i)) then
!   else 
!   !determine time to perihelion passage
!   
!     call htrnselm0 (rvo(i,1:6),ele1(1:6),pp)
!   
!     if (dt.gt.pp) then !pp
!   !propagate partilce id to kicktime
!       dt1=pp
!  
!   call radau_fstep1(t-dt,i,rvo,dt1,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
! 
!   !kick it
!   
!    !check for complete failure
!     call random_number(rnd)
!     !rnd=rand()
!      
!      if(rnd.gt.failureprob) then !failureprob
! 	vmat=0._dp
!       else
!  !actual mitigation 
!    !calculate tangential, radial and out of plane directions 
!    astr=rvo(i,1:3)/sqrt(dot_product(rvo(i,1:3),rvo(i,1:3)))
!    astv=rvo(i,4:6)/sqrt(dot_product(rvo(i,4:6),rvo(i,4:6)))
!  !GAUSS R S W system 
! !    call crossp3d(astr,astv,asth)
! !    asth=asth/sqrt(dot_product(asth,asth))
! !    call crossp3d(asth,astr,astt)
! !    astt=astt/sqrt(dot_product(astt,astt))
! 
!   !GAUSS T N W system 
!     call crossp3d(astr,astv,asth)
!    asth=asth/sqrt(dot_product(asth,asth))
!    astt=astv
!    call crossp3d(astt,asth,astr)
!    astr=astr/sqrt(dot_product(astr,astr))
! 
!    !should kick be a random process or nominal?
!    if(kickrnd) then
!     call random_number(rnd2(1:3))
!    else
!     rnd2(1:3)=1._dp
!     kickmin(:)=0._dp
!    end if
!    !deflection change in asteroid's velocity
!      !define matrix that contains radial tangential and out of plane velocities
!        vscal(1:3)=rnd2(1:3)*(kickmax(1:3)-kickmin(1:3))+kickmin(1:3)
!   
!        vmat(1,1:3)=astr(1:3)*vscal(1)
!        vmat(2,1:3)=astt(1:3)*vscal(2)
!        vmat(3,1:3)=asth(1:3)*vscal(3)
!    
!    !deflection change in asteroid's velocity
!    do j=1,3
!     rvo(i,4:6)=rvo(i,4:6)+vmat(j,1:3)
!    end do
!  end if
!       kicked(i)=.true.
!   
!    !propagate particle id to end of original timestep
!       dt2=dt-dt1
!    call radau_fstep1(t-dt2,i,rvo,dt2,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
!    
!       rv(i,1:6)=rvo(i,1:6)
!     write(unit=uc(1,9),fmt=*)i+m0cur,t+dt1-dt,t+tstart+dt1-dt,&
!     sign(sqrt(dot_product(vmat(1,:),vmat(1,:))),vscal(1)),&
!     sign(sqrt(dot_product(vmat(2,:),vmat(2,:))),vscal(2)),&
!     sign(sqrt(dot_product(vmat(3,:),vmat(3,:))),vscal(3))
!    
!       end if !failureprob
!      
!   end if !pp
!  end do
!  end if
! end if
!   
! return
! end subroutine
!********************************************************************************************
! subroutine kikickp(t,dt,rvo,rv,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
! !kick at pericenter
! implicit none
! 
! real(kind=dp),dimension(7,inbody,3),intent(inout)::c_bo,c_eo
! real(kind=dp),intent(in)::c_h(8),c_xc(8),c_vc(7),c_c(21),c_d(21),c_r(28)
! real(kind=dp),intent(in)::t,dt
! integer(kind=idp)::i
! real(kind=dp)::dt1,dt2,pp,rnorm,rnd
! real(kind=dp),dimension(1:inbody,1:6)::rv,rvo
! real(kind=dp),dimension(1:3)::LRL,LRLn,L
! real(kind=dp),dimension(1:6)::ele1
! 
! 
! if (t+tstart.gt.kicktime) then
! 
!  if(all(kicked)) then
!  else
!  do i=4,inbody
!   if(kicked(i)) then
!   else
! 
!  
!   
!   !determine time of perihelion passage
!   
!   call htrnselm0 (rvo(i,1:6),ele1(1:6),pp)
!   
!   if (dt.gt.pp) then
!   
!   !propagate partilce id to kicktime
!   dt1=pp
!   
!    !db
!   ! write(*,*)'******* BEFORE *********************'
!   ! write(*,*)rv(i,1:6)-rvo(i,1:6)
!   ! write(*,*)'**********************************'
!    !edb
!   
!   call radau_fstep1(t-dt,i,rvo,dt1,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
!       
!       
! !       !db 
! !    call crossp3d(rv(i,1:3),rv(i,4:6),L)
! ! !   
! !    call crossp3d(rv(i,4:6),L,LRL)
! !    rnorm=Sqrt(Dot_Product(rv(i,1:3),rv(i,1:3)))
! !    LRL(1:3)=LRL(1:3)/kgc2-rv(i,1:3)/rnorm
! !    LRLn=LRL/dot_product(LRL,LRL)
! !    write(*,*)LRLn(1:3),rv(i,1:3)/rnorm
! !   if(abs(dot_product(LRLn(1:3),rv(i,1:3)/rnorm)-1.d0).lt.0.005d0) then
! !edb
!     
!     
!   !kick it
!  !db!
!   !rnd=0.d0
!   !edb
!       rnd=rand()
!    
!    
!    rvo(i,4:6)=rvo(i,4:6)+rvo(i,4:6)/dot_product(rvo(i,4:6),rvo(i,4:6))*rnd*kickmax
!     kicked(i)=.true.
! 
!    !propagate particle id to end of original timestep
!    dt2=dt-dt1
!    call radau_fstep1(t-dt2,i,rvo,dt2,c_bo,c_eo,c_h,c_xc,c_vc,c_c,c_d,c_r)
!    
!    !db
!    !write(*,*)'*********  AFTER ******************'
!    !write(*,*)rv(i,1:6)-rvo(i,1:6)
!    !write(*,*)'**********************************'
!    !edb
!    rv(i,1:6)=rvo(i,1:6)
!    
!    
!               
!    write(unit=uc(1,9),fmt=*)i+m0cur,t,t+tstart,rnd*kickmax,pp,dt,dt1,dt2,ele1(6)
!    write(*,*)i+m0cur,'KICKED',pp,dt,dt1,dt2,t,t+tstart,&
!               ele1(6),rnd*kickmax
!                      
!     
!   end if
!   end if
!  end do
!  end if
! end if
!   
! return
! end subroutine
!***************************************************************************************! 

end module
