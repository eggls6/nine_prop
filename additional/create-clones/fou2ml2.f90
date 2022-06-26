program fou2ml
!uses orbfit .fou output files and extracts LOV clones
implicit none 

integer::nclones,eof,n
real*8::ke(1:6),pi,epoch,rad2grad
character(len=50)::dum,ini,infile,outfile
character(len=150)::out
write(*,*)'inputfile name'
read(*,*)infile
write(*,*)'clone epoch'
read(*,*)epoch
write(*,*)'number of clones'
read(*,*)nclones
ini=" no       a      e      I    Omeg    omeg    mean.an"
outfile="massless.inn"
out="a               e                 i                w            Om                M              Epoch"


pi=4.d0*atan(1.d0)
rad2grad=180.d0/pi

open(21,file=trim(infile),status="old")
open(22,file=trim(outfile),status='replace')
write(22,*)out

eof=0

do while (dum.ne.trim(ini))
 read(21,'(A)',iostat=eof)dum
 if(eof.ne.0) then
  write(*,*)'no initialization statement detected to initiate extraction'
  STOP
 end if
end do



do while(n<nclones)
read(21,*,iostat=eof)n,ke
 if(eof.ne.0) then
  write(*,*)'error in reading clone elements, possibly number of clones incorrect'
  STOP
 end if
 
write(22,*) ke(1:2),ke(3)*rad2grad,ke(5)*rad2grad,ke(4)*rad2grad,ke(6)*rad2grad,epoch
 
end do

close(22)
close(21)
end program 
 


