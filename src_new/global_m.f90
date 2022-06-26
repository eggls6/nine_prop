module global_m
  implicit none
  
  
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!                                        global variables definition
!
!                                         by Siegfried Eggl  20080110 
!------------------------------------------------------------------------------
!********************************************************************
!REAL PRECISION (32bit,64bit,128bit) 
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)   !15,307
  integer, parameter :: qp = selected_real_kind(33, 4931)
!INTEGER PRECISION 
  integer, parameter :: idp= selected_int_kind(16)

!CONSTANTS
  ! PI
  real(kind=dp),parameter::PI=3.14159265358979323846264338328_dp
  ! Degree to Rad conversion
  real(kind=dp),parameter::grad=0.0174532925199432957692369076849_dp
  !rad to degree conversion
  real(kind=dp),parameter::invgrad=57.2957795130823208767981548141_dp
  !epoch in Julian Date
  real(kind=dp),parameter::J2000=2451545._dp
  ! JPL DE ephemeris version will be taken from ephemeris file
  real(kind=dp)::DE4XX
  ! Gaussian Gravitational Constant will be taken from ephemeris file
  real(kind=dp)::kgc,kgc2
  !speed of light in au taken from ephemeris file
  real(kind=dp)::clight
  !1/clight^2
  real(kind=dp)::cm2
  !au in [km] taken from ephemeris file
  real(kind=dp)::aukm
  !Earth's equatorial radius 6378.1 km in [km]
  real(kind=dp),parameter::rearth=6378.1d0
  !Sun's mass in [kg]
  real(kind=dp),parameter::msunkg=1.98855E30
  !Earth's obliquity on J2000 in degrees
!  real(kind=dp),parameter::obl=23.439280444444444444444444444444444445_dp
   real(kind=dp),parameter::obl=23.439291111111_dp
   !Earth's phi angle to get J2000 
!  real(kind=dp),parameter::phi= -0.0000139667_dp
    real(kind=dp),parameter::phi=0._dp
!NUMBER OF BODIES
! Number of bodies in ephemeris files
  integer(kind=idp),parameter::neph=11
! Number of Massless Bodies in input file massless.inn  
  integer(kind=idp)::m0count
! number of massive bodies in input file massive.inn / same + 1
  integer(kind=idp)::mcount,mcp1  
! number of objects to be calculated in one run
  integer(kind=idp)::inbody
! total number of bodies to be calculated (massive + massless)
  integer(kind=idp)::totbody
! points towards the last particle id that was calculated in the previous run
  integer(kind=idp)::lastid
! number of massless bodies to be calculated in one run
  integer(kind=idp)::m0run
! current index of massless bodies (will be calculated in packages of m0run)
  integer(kind=idp)::m0cur
  

!MASS OF BODIES
  !ephemeris masses (meph) and masses x gravitational constant (gmeph)
  real(kind=dp),dimension(1:neph)::meph,gmeph
  !masses of massive bodies (mass) and masses x gravitational constant (gmass)
  real(kind=dp),dimension(:),allocatable::mass,gmass
  !order the ephemeris in terms of ascending mass for summations
  integer,parameter::ord(1:neph)=(/9,10,1,4,2,3,7,8,6,5,11/)
  !1...Mercury    -> 3
  !2...Venus      -> 5
  !3...Earth      -> 6
  !10..Moon       -> 2
  !4...Mars       -> 4
  !5...Jupiter    -> 10
  !6...Saturn     -> 9
  !7...Uranus     -> 8
  !8...Neptue     -> 7
  !9...Pluto      -> 1
  !11..Sun        -> 11
  

  
  !which number is the Sun (after ord transform) ?
  integer,parameter::idsun=11
  !which number is the Earth (after ord transform)?
  integer,parameter::idearth=6
    !which number is the Jupiter (after ord transform)?
  integer,parameter::idjup=10

!INTEGRATION TIME RELATED
  !start and end times for integration [JD!]
    real(kind=dp)::tstart,tend
  !start and end times for ephemeris
    real(kind=dp)::tephstart,tephend
    !time of simulation
    real(kind=dp)::tcalc
    !time of output
    real(kind=dp)::tout
  !epoch of initial conditions [JD]
    real(kind=dp),dimension(:),allocatable::initep  
    !minimum/maximum stepsize:
    real(kind=dp)::dtmin,dtmax
    !compensated (kahan) summation dummy
    real(kind=dp)::tky,tkt,tkc
    
   

!YARKOVSKY (see Bottke 2006)
     !solar constant 1367 W/m^2
     real(kind=dp),parameter::solconst=1367._dp 
     !C_diurnal, C_seasonal, Sin(eps_lag), Cos(eps_lag)
     real(kind=dp),allocatable,dimension(:,:)::yark_p
!    for output drift in semimajor axis due to yarkovsky effect [AU/Myrs] 
!     dadt_yark(1): diurnal effect,     dadt_yark(2): seasonal effect.
     real(kind=dp),allocatable,dimension(:,:)::dadt_yark

!     !spin direction (stays constant)
!     real(kind=dp),allocatable,dimension(:,:)::yb
!     real(kind=dp)::yeps


!TROJAN MOTION
      integer,allocatable,dimension(:)::trojan


!!PATH TO EPMEMERIS  
!  character(len=80)::path2eph
  

!AUXILLIARY VARIABLES
! current unit, see output_m
  integer::uc(3,11)
  !print progress to stdout?
  logical::showprog
!type of input variable
 character(len=2)::input
 
 !output style variables
  character(len=2)::tdigc,tendc
 
 !output of massive bodies as well?
  logical::massiveout 
  !output of ephemeris as well?
  logical::ephemout
 !output file choices
  logical::outegy,outbco,outhco,outbele,outhele,outdef
  logical::outbinbco,outbinhele,outbinhco
  logical::outtroj,outcc
  
  !matrices for transforming J2000 to ICRS to Ecliptic coordinates at J2000 and vice versa
  real(kind=dp),dimension(3,3)::jrote,erotj

!RADAU INTEGRATOR

  !integrator one step accuracy
  real(kind=dp)::eps
  !radau integrator constants
  real(kind=dp),allocatable,dimension(:,:,:)::c_b,c_e
  real(kind=dp),allocatable,dimension(:)::c_h,c_vc,c_xc,c_c,c_d,c_r
  
  
  
  
!CLOSE ENCOUNTER CHECK
  !which ephemeris body should be checked for close encounters
  integer::ccid  
  
  !  construct bplane at a distance of n x sphere of influence
  integer::nsoi
  
  ! number of close and deep close encounters per body
  integer(kind=idp),allocatable,dimension(:)::ncc,ndcc
  
  !limit distance for /maximum stepsize during close encounter 
  real(kind=dp)::cclim,ccdtmax
  !state vectors of the previous step
  real(kind=dp),dimension(:,:),allocatable::rvo  
  !partilce entered sphere of cclim (important to determine the point of entry)?
  logical,dimension(:),allocatable::ccenter
  
  !first bplane?
  logical,dimension(:),allocatable::fbplane
  !first bplane coordinates for all minor planets (inbody,3) with second dimension: 1...tpxi, 2...tpzeta, 3...time 
  real(kind=dp),dimension(:,:),allocatable::fbpcoo
  
  !geocentric distances/minimum distances and times
  real(kind=dp),dimension(:),allocatable::dagminn,dagmina,tdagminn,tdagmina
  real(kind=dp),dimension(:),allocatable::vcc,vinf

  !geocentric distances/minimum distances and times for each close approach with the Earth
   real(kind=dp),dimension(:,:,:),allocatable::loccc  !loccc(clone,number of cc, 1...6): 1:tmin, 2:rmin, 3:vmin, 4: tp xi, 5:tp zeta, 6:tp t_enter 

  !escape velocities for all massive bodies
  real(kind=dp),dimension(:),allocatable::ve2
 
 !start time for close encounter check
  real(kind=dp)::cct0
 
  !kinetic impactor
  character::kick   !used to determine whether kicked at a given time ('t') or at pericenter ('p')
  logical,dimension(:),allocatable::kicked
  logical::kickrnd !should kick contain random number? i.e. should kick be nominal or random process between min and max
  !real(kind=dp),dimension(1:3)::kickmin,kickmax !kickmax(1)=radial, kickmax(2)=tangential, kickmax(3)=out of plane (+ angular momentum)
  real(kind=dp)::kicktime,dkicktime
  !individual kicktime due to uncertainties
  real(kind=dp),dimension(:),allocatable::kicktime1
  real(kind=dp)::failureprob
  real(kind=dp)::betaf(1:3)   !range of beta factors (min, nominal, max)
  real(kind=dp)::ejphi        !max angle between (negative) impact and ejecta direction [deg]
  real(kind=dp)::vimph(1:3)   !heliocentric impact velocity [km/s]
  real(kind=dp),dimension(1:3)::msc,mast     !mass range of spacecraft/asteroid (min, nominal, max)
 
  
 !************ USE THE FOLLOWING IDEA PERHAPS LATER FOR Modified Target plane *****************
  !how many values are saved to interpolate (minimum 3 to identify close encounter)
!  integer,parameter,nsave=3
  !dummy variables for interpolating hyperbolic bodycentric encounter
  !rvold saves rv(i,j,nsave) for bodies i with coo/v j and times nsave
  !dold saves d(i,nsave) for body i with respect to target body for times nsave  
  !ddold saves differentials of dold to see whether a minimum distance has been reached 
  !real(kind=dp),dimension,allocatable::rvold(:,:,:),dold(:,:),ddold(:,:)


 

!################################################################################
   public::fopen
   public::fclose
   public::init
   public::invert
   public::crossp3d
   public::median

!****************************************************
 contains
!**************************************************** 
  subroutine fopen(name,un,hl)
  !open file with filename "name" in unit "un" and read "hl" header lines
 implicit none
  integer,intent(in)::un,hl
  character(len=50),intent(in)::name
  integer::i
  logical::ex
  character::dc
  
 inquire(file=trim(name),exist=ex)
 if(ex) then
     open(unit=un, file=trim(name), status="old", action="read")
     !read header lines
     do i=1,hl
      read(unit=un,fmt=*) dc
     end do
 else 
    write(*,*)'No input file found:',trim(name)
    STOP
 end if
 return
 end subroutine
 !************************************************
 subroutine fclose(un)
 !close file with unit un
 implicit none
 integer,intent(in)::un 
 close(un)
 return
 end subroutine
 !************************************************
  subroutine init(rv,rvold,um,um0,hl,hl0)
 !initializes global variables as well as local rv and rvinit
  implicit none
  integer::i
  integer,intent(inout)::um,um0,hl,hl0
  real(kind=dp),intent(inout),dimension(:,:)::rv,rvold
  real(kind=dp)::oblg,phig
  real(kind=dp),dimension(3,3)::rx,rz
  
   !initialize local variables
  oblg=obl*grad
  phig=phi*grad
  rv(:,:)=0._dp
  rvold(:,:)=0._dp
 
  !initialize global variables
  mass(:)=0._dp
  dadt_yark(:,:)=0
  !units for massive and massless input files
  um=40
  um0=41
  !number of header lines in massive and massless input files
  hl=1
  hl0=1
  
  !close encounter measures
  ccenter(:)=.true.
  ncc(:)=0
  ndcc(:)=0
  dagminn(:)=1._dp
  dagmina(:)=1._dp
  tdagmina(:)=0._dp
  tdagminn(:)=0._dp
  vinf(:)=0._dp
  vcc(:)=0._dp
  fbplane(:)=.true.
  fbpcoo(:,:)=0._dp
  !matrices for transforming ICRS to Ecliptic coordinates at J2000 and vice versa
  !J2000->ICRS->Ecliptic
!   i2e(1,:)=(/1._dp,0._dp,0._dp/)
!   i2e(2,:)=(/0._dp,cos(oblg),sin(oblg)/)
!   i2e(3,:)=(/0._dp,-sin(oblg),cos(oblg)/)
!   !Following Hilton & Hohenkerk 2004
!   j2i(1,:)=(/0.9999999999999938_dp,7.1d-8,-8.6d-8/)
!   j2i(2,:)=(/-7.1d-8,0.9999999999999938_dp,-2.6d-8/)
!   j2i(3,:)=(/8.6d-8,2.6d-8,0.99999999999999598_dp/)
! 
!    !USE DIFFERENCE IN ROTATION BETWEEN ICRS AND J2000
!    jrote=matmul(i2e,j2i)
   !IGNORE DIFFERENCE BETWEEN ICRS AND J2000 (JPL only uses obliquity rotation matrix comment previous line and uncomment the following)
 !  jrote=i2e

!TOTAL MATRIX FOLLOWING SOUAMI 2012
  
   
 rx(1,:)=(/1._dp,0._dp,0._dp/)
 rx(2,:)=(/0._dp,cos(oblg),sin(oblg)/)
 rx(3,:)=(/0._dp,-sin(oblg),cos(oblg)/)
  
 rz(1,:)=(/cos(phig),sin(phig),0._dp/)
 rz(2,:)=(/-sin(phig),cos(phig),0._dp/)
 rz(3,:)=(/0._dp,0._dp,1._dp/)
 
 
!    j2i(1,:)=(/0.9999999999999938_dp,7.1d-8,-8.6d-8/)
!    j2i(2,:)=(/-7.1d-8,0.9999999999999938_dp,-2.6d-8/)
!    j2i(3,:)=(/8.6d-8,2.6d-8,0.99999999999999598_dp/)
 
 !J2000->ICRS ~ Ecliptic
! jrote=matmul(rx,rz)
 
  jrote=rx
 
 !Ecliptic -> ICRS-> J2000
  
 call invert(jrote,erotj)


!***************************************
! KINETIC IMPACTOR (2cm/s in au/D)
 !kickmax=1.1551d-8
 !kickmin=0.d0
! kicktime=2458849.500000
 
 !kicktime=2459580.500000 !2022,1.1.
 !kinetic=.true.
 kicked(:)=.false.
  
  return
 end subroutine
!************************************************
 subroutine invert(m,minv)
 !inverts 3x3 matrix
 implicit none
 real(kind=dp),dimension(3,3),intent(in)::m
 real(kind=dp),dimension(3,3),intent(out)::minv
 real(kind=dp),dimension(3)::x0,x1,x2,dum
 real(kind=dp)::deta
 
 x0(1:3)=m(1:3,1)
 x1(1:3)=m(1:3,2)
 x2(1:3)=m(1:3,3)
 
 call crossp3d(x1,x2,dum)
 
 deta=Dot_Product(x0,dum)
 
 if (deta.eq.0._dp) then
   write(*,*)'Error in subroutine invert, module global_m: could not invert martix' 
 end if
 
 call crossp3d(x1,x2,minv(1,1:3))
 call crossp3d(x2,x0,minv(2,1:3))
 call crossp3d(x0,x1,minv(3,1:3))
 
 minv=minv/deta
  
 return
 end subroutine
 
!************************************************************ 
subroutine crossp3d(a,b,c)
!--------------------------------
!  3 dimensional crossproduct
!  a,b .... input vectors
!  c   .... output vector
!  c = a x b
!------------------------------- 
implicit none
real(kind=dp),dimension(1:3),intent(in)::a,b
real(kind=dp),dimension(1:3),intent(out)::c

 c(1)=a(2)*b(3)-b(2)*a(3)
 c(2)=a(3)*b(1)-b(3)*a(1)
 c(3)=a(1)*b(2)-b(1)*a(2)

return
end subroutine crossp3d
!***************************************************************
!*****************************************************************************
subroutine median(cnt,a,med)
implicit none
integer(kind=idp),intent(in)::cnt
real(kind=dp), dimension(1:cnt),intent(in)::a
integer(kind=idp),dimension(1:cnt)::ind
real(kind=dp)::half,med
logical::gerade


call bubble_sort_a(cnt,a,ind)

half=real(cnt)/2.d0

if (abs(half-aint(half)).le.0.1d0) then
   gerade=.true.
else
   gerade=.false.
end if

!write(*,*)'half gerade?',half,gerade

if (gerade) then
   med=0.5d0*(a(int(half))+a(int(half)+1))
else
   med=a(nint(half))
end if

return
end subroutine median

!****************************************************************************

subroutine bubble_sort_a ( n, a, ind )

!*****************************************************************************
!
!! BUBBLE_SORT_A ascending sorts a real vector using bubble sort.
!
!  Discussion:
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!    An associated index vector is adjusted to correspond to the
!    changes in the real vector.
!
!  Modified:
!
!    09 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
!    Input/output, real IND(N).
!    On input, an integer array.
!    On output, the entries of this array have been shifted
!    in a way that corresponds to the changes to the real vector.
!
  implicit none

  integer(kind=idp):: n

  real(kind=dp):: a(n)
  integer(kind=idp):: i
  integer(kind=idp):: ind(n)
  integer(kind=idp):: j

  do i = 1, n-1
    do j = i+1, n
      if ( a(i) > a(j) ) then
        call r_swap ( a(i), a(j) )
        call i_swap ( ind(i), ind(j) )
      end if
    end do
  end do

  return
end subroutine bubble_sort_a
!****************************************************************************
subroutine r_swap ( x, y )

!*****************************************************************************80
!
!! R_SWAP swaps two R's.
!
!  Modified:
!
!    29.1. 2007 Eggl
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real(kind=dp):: x,y,z

  z = x
  x = y
  y = z

  return
end subroutine r_swap
!******************************************************************************
subroutine i_swap ( x, y )

!*****************************************************************************
!
!! Swaps two Integers
!
!  Modified:
!
!    29.1. 2007 Eggl
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  integer(kind=idp):: x,y,z

  z = x
  x = y
  y = z

  return
end subroutine i_swap



!////////////////////////////////////////////////////////////////////
end module global_m