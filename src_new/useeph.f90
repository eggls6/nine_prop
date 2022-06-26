program useeph
!uses JPL planetary ephemeris  
use getephem

implicit none
integer::i,n,j
integer,parameter::np=11
real*8::time,tini,tfin,tend,k,clight
real*8::dt
real*8::rv(1:np,1:6),mass(1:np)


open(unit=21,file='planetdat.out',status='unknown')

call getparam(mass,tini,tfin,k,clight)


time=tini
dt=31.41552d0

do i=1,3000
        time=time+dt

        call getrv(time,rv)

        do j=1,np
        write(21,*)time,j,rv(j,1:6),mass(j)
        end do
end do 



write(*,*)'all done'
write(*,*)'parameters: k',k
write(*,*)'parameters: clight',clight
write(*,*)'tend',tfin
write(*,*)'tini',tini
close(21)

end program