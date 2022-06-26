program genclones
!generates clones for nine-prop propagator using pseudo Monte Carlo approach with Gaussian Errors following Bancelin Thesis 2011
!input : needs keplerian or JPL cometary element covariance matrix, 
!        needs lapack routines
!compile with gfortran -llapack 
 
implicit none

integer::i,j,lwork,liwork,info,nclones
integer,dimension(:),allocatable::iwork
real*8,dimension(6,6)::covceke,covke,covceket,covce,eigenv
real*8,dimension(6)::ce,ke,elem
real*8,dimension(:),allocatable::w,work
real*8::epoc,mu,e,q,tp,rnd,pix2
real*8,parameter::k=0.0172020989d0
character::dumc,elec

 covceke=0.d0
 covceket=0.d0
 covke=0.d0
 covce=0.d0
 ce=0.d0
 
 !asteroids are considered as test-particles
 mu=k*k
 pix2=8.d0*atan(1.d0)
 
 liwork=400
 lwork=liwork
 allocate(w(lwork),iwork(liwork),work(lwork))
 
 
 
 open(21,file="genclones.inn",status='old')
 read(21,*)dumc
 read(21,*)elec
 read(21,*)dumc
 read(21,*)epoc
 read(21,*) dumc
 read(21,*)ce(:)
 read(21,*)dumc
 
 do i=1,6
         read(21,*)covce(i,:)
 end do
 read(21,*)dumc
 read(21,*)nclones
 close(21)

 
 select case (elec)
  case('K','k')
  covke=covce
  
  case default
 
 q=ce(2)
 e=ce(1)
 tp=ce(3)
 
 !transform cometary to keplerian elemens
 ke(1)=q/(1.d0-e)
 ke(2)=e
 ke(3)=ce(6)
 ke(4)=ce(5)
 ke(5)=ce(4)
 ke(6)=modulo((epoc-tp)*sqrt((1.d0-e)**3*mu/q**3),pix2)*360.d0/pix2
 

 
 !construct transition matrix from cometary to keplerian elements
 covceke(1,1)= q/(1.d0-e)**2.d0
 covceke(1,2)= 1.d0/(1.d0-e)
 covceke(2,1)= 1.d0
 covceke(3,6)= 1.d0
 covceke(4,5)= 1.d0
 covceke(5,4)= 1.d0
 covceke(6,1)= -3.d0*(1.d0-e)**2.d0*(epoc-tp)*mu/(2.d0*q**3.d0*sqrt((1.d0-e)**3.d0/q**3.d0*mu))
 covceke(6,2)= -3.d0*(1.d0-e)**3.d0*(epoc-tp)*mu/(2.d0*q**4.d0*sqrt((1.d0-e)**3.d0/q**3.d0*mu))
 covceke(6,3)= -sqrt((1.d0-e)**3.d0/q**3.d0*mu)
 
 covceket=transpose(covceke)
 
 covke=matmul(covceke,matmul(covce,covceket))
 
 end select
 
 do i=1,6
 write(*,*) covke(i,:)
 end do
 
 
 call DSYEVD('V','U',6,eigenv,6,w,work,lwork,iwork,liwork,info)
 
 if (info.ne.0) then
  write(*,*)'Error in Lapack subroutine DSYEVD'
  STOP
 end if
 
 
 do i=1,6
         write(*,*)'eigenvector #',i 
         write(*,*)'V',eigenv(i,:)
         write(*,*)'corresponding eigenvalue', w(i)
         write(*,*)'w*V',w(i)*eigenv(i,:)
         write(*,*)'******************************'
 end do
 
 !check wheter eigenvectors are of unit lenght
 if(abs(dot_product(eigenv(6,1:6),eigenv(6,1:6))-1.d0).lt.1d-12) then
   write(*,*)'normalizing largest eigenvector (EV)...'
   eigenv(6,1:6)=eigenv(6,1:6)/sqrt(dot_product(eigenv(6,1:6),eigenv(6,1:6)))
   write(*,*)'normalized EV',eigenv(6,1:6)
 end if
 

 open(22,file="massless.inn",status='replace')
 write(22,*) 'a [au]  e    i[°]  om[°] Om[°] M[°]  epoch'
 
 !output nominal orbit
 write(22,*)ke,epoc
 !output clones
 do i=1,nclones
   call random_number(rnd)
    elem(:)=ke(:)+(rnd-0.5d0)*6.d0*sqrt(w(6))*eigenv(6,:)
    
    !correct for negative angles
    do j=1,3
     if(elem(j+3).lt.0.d0) then 
        ke(j+3)=ke(j+3)+360.d0
     end if
    end do
     if(elem(2).lt.0.d0) then !do not consider retrograde orbits
       elem(2)=abs(elem(2))
     end if
   write(22,*)elem(:),epoc
 end do
 
 close(22)
 
  write(*,*)ke(6),((epoc-tp)*sqrt((1.d0-e)**3.d0*mu/q**3.d0)/pix2)*360.d0,(epoc-tp)*sqrt((1.d0-e)**3.d0*mu/q**3.d0)
 
end program