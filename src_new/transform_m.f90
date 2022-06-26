module transform_m
  use global_m
  use getephem_m
  
  implicit none

  public::htrnsko
  public::htrnsel
  public::trojtrnsko
  public::trojtrnsel
  public::hekoo
  public::bakoo
  public::bodykoo
  public::bodykoo1
  public::topolar2
  public::tocartesian2

  contains

!/////////////////////////////////////////////////////////////////////////////////
!                              Transformationen fuer Eingabe und Ausgabe 
!
!                                         by Siegfried Eggl  20080110 
!--------------------------------------------------------------------------------
!*****************************************************************

!....BERECHNUNG DER HELIOZENTRISCHEN KOORDINATEN UND GESCHWINDIGKEITEN
!....AUS DEN BAHNELEMENTEN


        SUBROUTINE htrnsko(ele,n,rv)
         use global_m
         implicit none   
         integer(kind=idp)::i,n
         real(kind=dp),dimension(:,:)::ele,rv
          real(kind=dp)::xh1,xh2,xh3,vh1,vh2,vh3,a,e,incl,om,gom,m       
          real(kind=dp)::cosgom,singom,cosi,sini,px,py,pz,qx,qy,qz,ah,exan,cappaq
          real(kind=dp)::cosex,sinex,rh,bh,cosom,sinom,ea,grd,p,kk
  
         p=pi
         grd=grad
         kk=kgc

!$omp parallel do default(private)   &
!$omp firstprivate(kk,grd,p,n) &
!$omp shared(rv,mass,meph,ele) 
    do  i=1,n
          a=ele(i,1)
          e=ele(i,2)
          incl=ele(i,3)*grd
          om=ele(i,4)*grd
          gom=ele(i,5)*grd
          m=ele(i,6)*grd

          cappaq=meph(idsun)+mass(i)
          
         
          call nr(m,e,ea)
 !        cosea=cos(ea)
 !         cosphi=(cosea-e)/(1._dp-e*cosea)
         
          IF (E .EQ.0.D0)  then
             OM  = 0.D0
          end if
          IF (INCL .EQ.0.D0) then
             GOM  = 0.D0
          end if
          COSOM = COS(OM )
          SINOM = SIN(OM )
          COSGOM = COS(GOM)
          SINGOM = SIN(GOM )
          COSI = COS(INCL )
          SINI = SIN(INCL ) 
          PX= COSOM*COSGOM-SINOM*SINGOM*COSI
          PY= COSOM*SINGOM+SINOM*COSGOM*COSI
          PZ= SINOM*SINI
          QX=-SINOM*COSGOM-COSOM*SINGOM*COSI
          QY=-SINOM*SINGOM+COSOM*COSGOM*COSI
          QZ= COSOM*SINI
!          AH = SQRT((1-E )/(1+E ))*SQRT((1._dp-cosphi(n))/(1._dp+cosphi(n)))
!frueher... *DTAN(V /2._dp)
!          EXAN = 2.D0*DATAN(AH)
          EXAN=ea
 !         ta=2._dp*DATAN(SQRT((1._dp+e)/(1._dp-e))*DTAN(ea/2._dp))
 
          COSEX = COS(EXAN)
          SINEX =SIN(EXAN)

          AH = A *SQRT(1._dp-E *E )
          XH1  = A *PX*(COSEX-E )+AH*QX*SINEX
          XH2  = A *PY*(COSEX-E )+AH*QY*SINEX
          XH3  = A *PZ*(COSEX-E )+AH*QZ*SINEX
          RH = SQRT(XH1 *XH1 +XH2 *XH2 +XH3 *XH3 )

          if(A.eq.0._dp.or.RH.eq.0._dp) then
             AH=1._dp
          else
             AH = SQRT(cappaq )/(SQRT(A )*RH)
          end if
  
          BH = A *SQRT(1._dp-E *E )*COSEX
          VH1  = AH*(-A *PX*SINEX+BH*QX)
          VH2  = AH*(-A *PY*SINEX+BH*QY)
          VH3  = AH*(-A *PZ*SINEX+BH*QZ)

   
          rv(i,1)=  XH1
          rv(i,2)=  XH2
          rv(i,3)=  XH3
          rv(i,4)=VH1*kgc
          rv(i,5)=VH2*kgc
          rv(i,6)=VH3*kgc                   
end do
  !$omp end parallel do  


        RETURN
      END  subroutine htrnsko
!*******************************************************
subroutine he2icrs(rv,n)
!rotates heliocentric ecliptic rv vectors into ICRS system of JPL Ephemeris
implicit none
integer(kind=idp)::n,i,j
real(kind=dp),dimension(:,:)::rv


!rotate positions
do i=1,n
   rv(i,1:3)=matmul(erotj(:,:),rv(i,1:3))
end do
!rotate velocities
do i=1,n
   rv(i,4:6)=matmul(erotj(:,:),rv(i,4:6))
end do

return
end subroutine
 
!*********************************************************************

subroutine icrs2he(rv,n)
!rotates heliocentric ecliptic rv vectors into ICRS system of JPL Ephemeris at reference J2000
implicit none
integer(kind=idp)::i,j,n
real(kind=dp),dimension(:,:)::rv

!rotate positions
do i=1,n
   rv(i,1:3)=matmul(jrote(:,:),rv(i,1:3))
end do 

!rotate velocities
do i=1,n
   rv(i,4:6)=matmul(jrote(:,:),rv(i,4:6))
end do

return
end subroutine




!*****************************************
!...Calculates Keplerian Orbital Elements from heliocentric position and velocity vectors

SUBROUTINE htrnsel (rv,ele,mas,n)
use global_m
implicit none
integer(kind=idp)::i,j,n,dim 
real(kind=dp),dimension(:,:),intent(in)::rv     
real(kind=dp):: ele(:,:),mas(:)
real(kind=dp)::e2,eanom,Lprojxy,U(1:3,1:3)
real(kind=dp):: tanom,Hnorm,sinE,cosE
real(kind=dp),dimension(1:3)::LRL,L,r,v,H,node,rU
real(kind=dp)::cappaq,a,e,incl,om,gom,mmm,L2,rnorm,nnorm,igrad,p
 
p=pi
igrad=invgrad
       
dim=3

!$omp parallel do default(private)   &
!$omp firstprivate (igrad,p,inbody,dim)  &
!$omp shared(rv,mass,ele)
 do i=1,n    
       cappaq=meph(idsun)+mas(i)

!*******************************************************
! LETS BE SERIOUS HERE

          r=rv(i,1:3)
          v=rv(i,4:6)/kgc
!calculation of unit mass angular momentum vector L and its projection onto xy-plane
        call crossp3d(r,v,L)
        L2=Dot_Product(L,L)
        Lprojxy=SQRT( L(1)*L(1)+L(2)*L(2) )

 ! inclination "incl"
       incl=ATAN2(Lprojxy,L(3))*igrad

!argument of the ascending node "gom"
       gom=Atan2(L(1),-L(2))*igrad
       if (gom.lt.0._dp) then
           gom = gom+360._dp
        end if

! Laplace-Runge-Lenz vector
        call crossp3d(v,L,LRL)
        rnorm=Sqrt(Dot_Product(r(:),r(:)))

        LRL(:)=LRL(:)/cappaq-r(:)/rnorm
         
! eccentricity "e" = |LRL|
        e=Sqrt(Dot_Product(LRL,LRL))

        if(e.lt.1.d-13)then 
           e=0._dp
           e2=0._dp
           om=0._dp
         else
           e2=e*e
         end if
         
! semi-major axis "a"
     ! a=L2/(cappaq*(1._dp+e2))
       a=L2/(cappaq*(1._dp-e2))    
   
   !calculation of vector "node" pointing to ascending node    
    ! attention: will change node for retrograde orbits!!!!
         node(1)=-L(2)*L(3)
         node(2)=L(1)*L(3)
         node(3)=0._dp

         nnorm=Sqrt(Dot_Product(node,node))
 
       if(incl.lt.1.d-13) then
          gom=0._dp
         ! argument of pericenter "om" = angle of LRL to x-axis
          om=Acos(LRL(1)/Sqrt(Dot_Product(LRL,LRL)))
              if(LRL(2).lt.0._dp) then
                 om=360._dp-om*igrad
              else
                 om=om*igrad
              end if
  
        else

          call crossp3d(L,node,H)
          Hnorm=Sqrt(Dot_Product(H,H))
          om=Atan2(Dot_Product(LRL,H)*nnorm,Dot_Product(LRL,node)*Hnorm)*igrad
              if(om.lt.0._dp) om=om+360._dp 
  
        end if  

          if(e.lt.1.2d-13.and.incl.le.1.d-13) then
              gom=0._dp
              om=0._dp

              !mean anomaly 'mmm'
              mmm=Atan2(r(2),r(1))*igrad
           if(mmm.lt.0._dp) mmm=mmm+360._dp

           elseif (e.lt.1.2d-13.and.incl.gt.1.d-13) then             
     
              !transformation of positionvector r into coordinate-system U spanned by Angular Momentum vector, vector Sun-ascending Node, 
              !and their crossproduct H, corresponding to motion of rU in the orbital plane. 
              

              call crossp3d(L,node,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=node(:)/nnorm
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))

              !mean anomaly 'mmm'
              mmm=tanom*igrad     
              if(mmm.lt.0._dp) mmm=mmm+360._dp
      
   !  ??       om=0._dp

            elseif (incl.lt.1.d-13.and.e.gt.1.d-13) then
              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
                      call crossp3d(L,LRL,H)
                      Hnorm=Sqrt(Dot_Product(H,H))
                      tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm) 
     

              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._dp+e*Cos(tanom))
               sinE=Sqrt(1._dp-e2)*Sin(tanom)/(1._dp+e*Cos(tanom))
               eanom=ATAN2(sinE,cosE)  

              !mean anomaly 'mmm' via Kepler's equation
              mmm=(eanom-e*sinE)*igrad              

             if(mmm.lt.0._dp) mmm=mmm+360._dp

              else
  
             call crossp3d(L,LRL,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=LRL(:)/e
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))


              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._dp+e*Cos(tanom))
               sinE=Sqrt(1._dp-e2)*Sin(tanom)/(1._dp+e*Cos(tanom))
               eanom=Atan2(sinE,cosE)  
                       
              !mean anomaly 'mmm' via Kepler's equation
             mmm=(eanom-e*sinE)*igrad
             if(mmm.lt.0._dp) mmm=mmm+360._dp
           end if
      
          if(e.lt.0._dp.or.e.gt.1._dp) e=1._dp
          if(a.lt.0._dp.or.a.gt.huge(a)) a=0._dp 
          if(incl.gt.180._dp.or.incl.lt.0._dp) om=-100._dp

          ele(i,1)=a
          ele(i,2)=e
          ele(i,3)=incl
          ele(i,4)=om
          ele(i,5)=gom
          ele(i,6)=mmm

          do j=4,6
            if(ele(i,j).gt.360._dp.or.ele(i,j).lt.0._dp.or.ele(i,j).ne.ele(i,j)) then
              ele(i,j)=-100._dp
            end if
          end do

       end do
!$omp end parallel do
    
       RETURN
     END subroutine htrnsel
!*************************************************************************************************

!*****************************************
!...Calculates Keplerian Orbital Elements from heliocentric position and velocity vectors for 1 massless particle

SUBROUTINE htrnselm0 (rv1,ele1,ppt)
use global_m
implicit none
integer(kind=idp)::i,j,n
integer(kind=idp),parameter::gd=3
real(kind=dp),dimension(1:6),intent(in)::rv1     
real(kind=dp),intent(out):: ele1(1:6),ppt
real(kind=dp)::e2,eanom,Lprojxy,U(1:3,1:3)
real(kind=dp):: tanom,Hnorm,sinE,cosE
real(kind=dp),dimension(1:3)::LRL,L,r,v,H,node,rU
real(kind=dp)::cappaq,a,e,incl,om,gom,mm,mmm,L2,rnorm,nnorm,igrad,p
 
p=pi
igrad=invgrad
       

       cappaq=meph(idsun)

!*******************************************************
! LETS BE SERIOUS HERE

          r=rv1(1:3)
          v=rv1(4:6)/kgc
         
!calculation of unit mass angular momentum vector L and its projection onto xy-plane
        call crossp3d(r,v,L)
        L2=Dot_Product(L,L)
        Lprojxy=SQRT( L(1)*L(1)+L(2)*L(2) )

 ! inclination "incl"
       incl=ATAN2(Lprojxy,L(3))*igrad

!argument of the ascending node "gom"
       gom=Atan2(L(1),-L(2))*igrad
       if (gom.lt.0._dp) then
           gom = gom+360._dp
        end if

! Laplace-Runge-Lenz vector
        call crossp3d(v,L,LRL)
        rnorm=Sqrt(Dot_Product(r(:),r(:)))

        LRL(:)=LRL(:)/cappaq-r(:)/rnorm
         
! eccentricity "e" = |LRL|
        e=Sqrt(Dot_Product(LRL,LRL))

        if(e.lt.1.d-13)then 
           e=0._dp
           e2=0._dp
           om=0._dp
         else
           e2=e*e
         end if
         
! semi-major axis "a"
     ! a=L2/(cappaq*(1._dp+e2))
       a=L2/(cappaq*(1._dp-e2))    
   
   !calculation of vector "node" pointing to ascending node    
    ! attention: will change node for retrograde orbits!!!!
         node(1)=-L(2)*L(3)
         node(2)=L(1)*L(3)
         node(3)=0._dp

         nnorm=Sqrt(Dot_Product(node,node))
 
       if(incl.lt.1.d-13) then
          gom=0._dp
         ! argument of pericenter "om" = angle of LRL to x-axis
          om=Acos(LRL(1)/Sqrt(Dot_Product(LRL,LRL)))
              if(LRL(2).lt.0._dp) then
                 om=360._dp-om*igrad
              else
                 om=om*igrad
              end if
  
        else

          call crossp3d(L,node,H)
          Hnorm=Sqrt(Dot_Product(H,H))
          om=Atan2(Dot_Product(LRL,H)*nnorm,Dot_Product(LRL,node)*Hnorm)*igrad
              if(om.lt.0._dp) om=om+360._dp 
  
        end if  

          if(e.lt.1.2d-13.and.incl.le.1.d-13) then
              gom=0._dp
              om=0._dp

              !mean anomaly 'mmm'
              mmm=Atan2(r(2),r(1))*igrad
           if(mmm.lt.0._dp) mmm=mmm+360._dp

           elseif (e.lt.1.2d-13.and.incl.gt.1.d-13) then             
     
              !transformation of positionvector r into coordinate-system U spanned by Angular Momentum vector, vector Sun-ascending Node, 
              !and their crossproduct H, corresponding to motion of rU in the orbital plane. 
              

              call crossp3d(L,node,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=node(:)/nnorm
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(gd,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))

              !mean anomaly 'mmm'
              mmm=tanom*igrad     
              if(mmm.lt.0._dp) mmm=mmm+360._dp
      
   !  ??       om=0._dp

            elseif (incl.lt.1.d-13.and.e.gt.1.d-13) then
              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
                      call crossp3d(L,LRL,H)
                      Hnorm=Sqrt(Dot_Product(H,H))
                      tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm) 
     

              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._dp+e*Cos(tanom))
               sinE=Sqrt(1._dp-e2)*Sin(tanom)/(1._dp+e*Cos(tanom))
               eanom=ATAN2(sinE,cosE)  

              !mean anomaly 'mmm' via Kepler's equation
              mmm=(eanom-e*sinE)*igrad              

             if(mmm.lt.0._dp) mmm=mmm+360._dp

              else
  
             call crossp3d(L,LRL,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=LRL(:)/e
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(gd,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))


              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._dp+e*Cos(tanom))
               sinE=Sqrt(1._dp-e2)*Sin(tanom)/(1._dp+e*Cos(tanom))
               eanom=Atan2(sinE,cosE)  
                       
              !mean anomaly 'mmm' via Kepler's equation
             mmm=(eanom-e*sinE)*igrad
             if(mmm.lt.0._dp) mmm=mmm+360._dp
           end if
      
          if(e.lt.0._dp.or.e.gt.1._dp) e=1._dp
          if(a.lt.0._dp.or.a.gt.huge(a)) a=0._dp 
          if(incl.gt.180._dp.or.incl.lt.0._dp) om=-100._dp

          
          ele1(1)=a
          ele1(2)=e
          ele1(3)=incl
          ele1(4)=om
          ele1(5)=gom
          ele1(6)=mmm

          do j=4,6
            if(ele1(j).gt.360._dp.or.ele1(j).lt.0._dp.or.ele1(j).ne.ele1(j)) then
              ele1(j)=-100._dp
            end if
          end do
      
          
!mean motion
          mm=kgc*sqrt(cappaq)*a**(-1.5_dp)
!time until next pericenter passage          
          ppt=(360.d0-mmm)*p/180._dp/mm
    
       RETURN
     END subroutine
!*************************************************************************************************


subroutine gauss(dim,a,b,c)
!----------------------------------------------------------
! Gauß'sches Eliminationsverfahren zur Lösung eines 
! Systems Linearer Gleichungen der Form 
!                                    A.c=b
! mit Teilpivotisierung
!
! dim....[Int] Dimension des Gleichungssystems
! a...[Real] dim x dim Koeffizientenmatrix des Gleichungssystems
! b...[Real] dim Vektor der rechten Seite der Gleichung
! c...[Real] dim  Vektor der Unbekannten
! dependencies: none
!----------------------------------------------------------
implicit none
  integer(kind=idp),intent(in)::dim
  integer(kind=idp)::i,j,k,pivot(1:dim),amax(1:1),idum
  real(kind=dp)::a(1:dim,1:dim+1),b(1:dim),dum(1:dim+1)
  real(kind=dp)::c(1:dim)

  do i=1,dim
     pivot(i)=i
  end do

  dum(:)=0._dp
  c(:)=0._dp
  a(:,dim+1)=b(:)
 
 do i=1,dim
!-----------------------------------------------------------------------
!  Partial Pivoting (nur Zeilen werden vertauscht)
!-----------------------------------------------------------------------
    amax(:)=maxloc(abs(a(:,i))) !Welche Zeile hat den größten Wert?
    idum=amax(1)
   if (idum.gt.i) then
      pivot(i)=idum
       dum(:)=a(i,:)                                !Sichere die Zeile die ersetzt wird
      a(i,:)=a(pivot(i),:)                       !Ersetze Zeile i mit jener mit größtem Wert 
       a(pivot(i),:)=dum(:)                 !Beende die Vertauschung, der gesicherte Wert ersetzt die ersetzende Zeile
     end if
!------------------------------------------------------
!   Gauss Elimination
!-----------------------------------------------------     
      a(i,:)=a(i,:)/a(i,i)
      if(i+1.le.dim) then
         do j=i+1,dim
            a(j,:)=a(j,:)-a(j,i)*a(i,:)         
        end do  
      end if
  end do
  
!-----------------------------------------------
!   Rückwärtssubstitution
!-----------------------------------------------

do i=0,dim-1
   j=dim-i  
   dum(1)=0._dp
   if (j+1.le.dim) then
      do k=j+1,dim
         dum(1)=dum(1)+a(j,k)*c(k)
      end do
   end if
   c(j)=(a(j,dim+1)-dum(1))
end do

!-----------------------------------------------------------------------
!  Back Pivoting enfällt, da nur Zeilen vertauscht wurden, und sich
!  somit zwar die Reihenfolge der Gleichungen, nicht aber die Koeffizienten
!  der Variablen geändert haben.
!-----------------------------------------------------------------------

return
end subroutine gauss

!******************************************************************************************
!....BERECHNUNG DER HELIOZENTRISCHEN KOORDINATEN UND GESCHWINDIGKEITEN
!....AUS DEN BAHNELEMENTEN FUER TROJANER
        SUBROUTINE trojtrnsko(ele,rv)
        use global_m
         implicit none   
         integer(kind=idp)::i,n
         real(kind=dp)::rv(:,:),ele(:,:)
          real(kind=dp)::xh1,xh2,xh3,vh1,vh2,vh3,a,e,incl,om,gom,m       
          real(kind=dp)::cosgom,singom,cosi,sini,px,py,pz,qx,qy,qz,ah,exan,cappaq
          real(kind=dp)::cosex,sinex,rh,bh,cosom,sinom,ea,p
          real(kind=dp)::avv(3)
   
    p=pi  
    do  i=1,inbody

          a=ele(i,1)
          e=ele(i,2)
          incl=ele(i,3)*grad
          om=ele(i,4)*grad
          gom=ele(i,5)*grad
          m=ele(i,6)*grad

!trojan defined in global_m

          if (trojan(i).le.1) then
            cappaq=mass(1)+mass(i)
          else
            cappaq=mass(1)+mass(trojan(i))
          end if

          call nr(m,e,ea)
 !        cosea=cos(ea)
 !         cosphi=(cosea-e)/(1._dp-e*cosea)
         
          IF (E .EQ.0.D0)  then
             OM  = 0.D0
          end if
          IF (INCL .EQ.0.D0) then
             GOM  = 0.D0
          end if
          COSOM = COS(OM )
          SINOM = SIN(OM )
          COSGOM = COS(GOM )
          SINGOM = SIN(GOM )
          COSI = COS(INCL )
          SINI = SIN(INCL ) 
          PX= COSOM*COSGOM-SINOM*SINGOM*COSI
          PY= COSOM*SINGOM+SINOM*COSGOM*COSI
          PZ= SINOM*SINI
          QX=-SINOM*COSGOM-COSOM*SINGOM*COSI
          QY=-SINOM*SINGOM+COSOM*COSGOM*COSI
          QZ= COSOM*SINI
!          AH = SQRT((1-E )/(1+E ))*SQRT((1._dp-cosphi(n))/(1._dp+cosphi(n)))
!frueher... *DTAN(V /2._dp)
!          EXAN = 2.D0*DATAN(AH)
          EXAN=ea
 !         ta=2._dp*DATAN(SQRT((1._dp+e)/(1._dp-e))*DTAN(ea/2._dp))
 
          COSEX = COS(EXAN)
          SINEX =SIN(EXAN)

          AH = A *SQRT(1._dp-E *E )
          XH1  = A *PX*(COSEX-E )+AH*QX*SINEX
          XH2  = A *PY*(COSEX-E )+AH*QY*SINEX
          XH3  = A *PZ*(COSEX-E )+AH*QZ*SINEX
          RH = SQRT(XH1 *XH1 +XH2 *XH2 +XH3 *XH3 )
    
          
          if(A.eq.0._dp.or.RH.eq.0._dp) then
             AH=1._dp
          else
             AH = SQRT(cappaq )/(SQRT(A )*RH)
          end if
           
  
          BH = A *SQRT(1._dp-E *E )*COSEX
          VH1  = AH*(-A *PX*SINEX+BH*QX)
          VH2  = AH*(-A *PY*SINEX+BH*QY)
          VH3  = AH*(-A *PZ*SINEX+BH*QZ)

   
          rv(i,1)=  XH1
          rv(i,2)=  XH2
          rv(i,3)=  XH3

          rv(i,4)=VH1*kgc
          rv(i,5)=VH2*kgc
          rv(i,6)=VH3*kgc                   


!         now rescale velocities to make them have the same angular velocity as the trojan host
          ! dr/dt = avv x r  so it should rest in a corotating frame
          
      if(trojan(i).le.1) then
      else    
          if(ele(trojan(i),1).eq.a) then
          else
            call crossp3d(rv(trojan(i),1:3),rv(trojan(i),4:6),avv)
            avv=avv/dot_product(rv(trojan(i),1:3),rv(trojan(i),1:3))
  
           call crossp3d(avv,rv(i,1:3),rv(i,4:6))
          end if
     end if
end do
        RETURN
      END  subroutine trojtrnsko

!*****************************************
!...BERECHNUNG DER BAHNELEMENTE AUS ORTS- UND GESCHWINDIGKEITS-VEKTOREN
!....INCLUSIVE QUADRANTENKORREKTUR 

SUBROUTINE trojtrnsel (rv,ele)
use global_m
implicit none
integer(kind=idp)::i,j,dim  
real(kind=dp),dimension(:,:):: ele,rv
real(kind=dp)::e2,eanom,Lprojxy,U(1:3,1:3)
real(kind=dp):: tanom,Hnorm,sinE,cosE
real(kind=dp),dimension(1:3)::LRL,L,r,v,H,node,rU
real(kind=dp)::cappaq,a,e,incl,om,gom,mmm,L2,igrad,rnorm,nnorm,p
!real(kind=dp)::avv
       
p=pi
igrad=invgrad 
dim=3

!$omp parallel default(private)   &
!$omp firstprivate(p,n,igrad) &
!$omp shared(rv,mass,ele)
!$omp do
 do i=1,inbody    

         do j=1,3
          r(j)=rv(i,j)
          v(j)=rv(i,3+j)
         end do
         
         if (trojan(i).le.1) then
            cappaq=meph(idsun)+mass(i)
            
          else
            cappaq=meph(idsun)+meph(trojan(i))
          end if


          
          
!calculation of unit mass angular momentum vector L and its projection onto xy-plane
        call crossp3d(r,v,L)
        L2=Dot_Product(L,L)
       Lprojxy=SQRT( L(1)*L(1)+L(2)*L(2) )

 ! inclination "incl"
       incl=ATAN2(Lprojxy,L(3))*igrad

!argument of the ascending node "gom"
       gom=Atan2(L(1),-L(2))*igrad
       if (gom.lt.0._dp) then
           gom = gom+360._dp
        end if

! Laplace-Runge-Lenz vector
        call crossp3d(v,L,LRL)
        rnorm=Sqrt(Dot_Product(r(:),r(:)))

        LRL(:)=LRL(:)/cappaq-r(:)/rnorm
         
! eccentricity "e" = |LRL|
        e=Sqrt(Dot_Product(LRL,LRL))

        if(e.lt.1.d-13)then 
           e=0._dp
           e2=0._dp
           om=0._dp
         else
           e2=e*e
         end if
         
! semi-major axis "a"
     ! a=L2/(cappaq*(1._dp+e2))
       a=L2/(cappaq*(1._dp-e2))    
   
   !calculation of vector "node" pointing to ascending node    
    ! attention: will change node for retrograde orbits!!!!
         node(1)=-L(2)*L(3)
         node(2)=L(1)*L(3)
         node(3)=0._dp

         nnorm=Sqrt(Dot_Product(node,node))
 
       if(incl.lt.1.d-13) then
          gom=0._dp
         ! argument of pericenter "om" = angle of LRL to x-axis
          om=Acos(LRL(1)/Sqrt(Dot_Product(LRL,LRL)))
              if(LRL(2).lt.0._dp) then
                 om=360._dp-om*igrad
              else
                 om=om*igrad
              end if
  
        else
   

          call crossp3d(L,node,H)
          Hnorm=Sqrt(Dot_Product(H,H))
          om=Atan2(Dot_Product(LRL,H)*nnorm,Dot_Product(LRL,node)*Hnorm)*igrad
              if(om.lt.0._dp) om=om+360._dp 
  
        end if  

          if(e.lt.1.2d-13.and.incl.le.1.d-13) then
              gom=0._dp
              om=0._dp

              !mean anomaly 'mmm'
              mmm=Atan2(r(2),r(1))*igrad
             if(mmm.lt.0._dp) mmm=mmm+360._dp

           elseif (e.lt.1.2d-13.and.incl.gt.1.d-13) then
        
              !transformation of positionvector r into coordinate-system U spanned by Angular Momentum vector, vector Sun-ascending Node, 
              !and their crossproduct H, corresponding to motion of rU in the orbital plane. 
              

              call crossp3d(L,node,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=node(:)/nnorm
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))

              !mean anomaly 'mmm'
              mmm=tanom*igrad     
              if(mmm.lt.0._dp) mmm=mmm+360._dp
      

            elseif (incl.lt.1.d-13.and.e.gt.1.d-13) then
              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
                      call crossp3d(L,LRL,H)
                      Hnorm=Sqrt(Dot_Product(H,H))
                      tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm) 
 
              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._dp+e*Cos(tanom))
               sinE=Sqrt(1._dp-e2)*Sin(tanom)/(1._dp+e*Cos(tanom))
               eanom=ATAN2(sinE,cosE)  

              !mean anomaly 'mmm' via Kepler's equation
              mmm=(eanom-e*sinE)*igrad              
 
!                write(*,*)'hey'

                 if(mmm.lt.0._dp) mmm=mmm+360._dp

              else
        
  
             call crossp3d(L,LRL,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=LRL(:)/e
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))


              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._dp+e*Cos(tanom))
               sinE=Sqrt(1._dp-e2)*Sin(tanom)/(1._dp+e*Cos(tanom))
               eanom=Atan2(sinE,cosE)  
                           
              !mean anomaly 'mmm' via Kepler's equation
             mmm=(eanom-e*sinE)*igrad
             if(mmm.lt.0._dp) mmm=mmm+360._dp
           end if
      
          if(e.lt.0._dp.or.e.gt.1._dp) e=1._dp
          if(a.lt.0._dp.or.a.gt.huge(a)) a=0._dp 
          if(incl.gt.180._dp.or.incl.lt.0._dp) om=-100._dp

          ele(i,1)=a
          ele(i,2)=e
          ele(i,3)=incl
          ele(i,4)=om
          ele(i,5)=gom
          ele(i,6)=mmm

          do j=4,6
            if(ele(i,j).gt.360._dp.or.ele(i,j).lt.0._dp.or.ele(i,j).ne.ele(i,j)) then
              ele(i,j)=-100._dp
            end if
          end do

       end do
!$omp end do       
!$omp end parallel
     
       RETURN
     END subroutine trojtrnsel

!   *****************************************************************
!...BERECHNUNG DER BAHNELEMENTE AUS DEN KARTHESISCHEN KOORDINATEN
!...UND GESCHWINDIGKEITEN 

  subroutine xtoel(x,ele,mu)

        IMPLICIT none
        real(kind=dp)::mu,rsquare,dotrv,vsquare,crossi,crossj,crossk
        real(kind=dp)::rvet,rsobrea,semia,ecosu,esinu,esquare
        real(kind=dp)::excen,anomu,anomalia,nodo,sinimod,incli
        real(kind=dp)::psubz,qsubz,argumento
        real(kind=dp):: X(1:6),ELE(1:6),V(1:3)
!	
        V(1)=X(4)/SQRT(mu)
        V(2)=X(5)/SQRT(mu)
        V(3)=X(6)/SQRT(mu)
        RSQUARE=X(1)**2+X(2)**2+X(3)**2
        DOTRV=X(1)*V(1)+X(2)*V(2)+X(3)*V(3)
        VSQUARE=V(1)**2+V(2)**2+V(3)**2
        CROSSI=X(2)*V(3)-X(3)*V(2)
        CROSSJ=X(3)*V(1)-X(1)*V(3)
        CROSSK=X(1)*V(2)-X(2)*V(1)
        RVET=SQRT(RSQUARE)
        RSOBREA=2._dp-RVET*VSQUARE
        SEMIA=RVET/RSOBREA
!  write(*,*)'xtoel',semia,rvet,vsquare,rsobrea
        ECOSU=1.D0-RSOBREA
        ESINU=DOTRV/SQRT(SEMIA)
        ESQUARE=ECOSU**2._dp+ESINU**2._dp
        EXCEN=SQRT(ESQUARE)
        ANOMU=ATAN2(ESINU,ECOSU)
!		A operacao ATAN2 da resultado entre -PI e +PI.
!     write(*,*)'esinu,ecosu,anomu',esinu,ecosu,anomu
        ANOMALIA=ANOMU-EXCEN*SIN(ANOMU)
!	MODCROSS=SQRT(CROSSI**2+CROSSJ**2+CROSSK**2)
        NODO=ATAN2(CROSSI,-CROSSJ)
        SINIMOD=SQRT(CROSSI**2+CROSSJ**2)      
        INCLI=ATAN2(SINIMOD,CROSSK)
        PSUBZ=X(3)/RVET*COS(ANOMU)-V(3)*SQRT(SEMIA)*SIN(ANOMU)
        QSUBZ=X(3)/RVET*SIN(ANOMU)+V(3)*SQRT(SEMIA)* (COS(ANOMU)-EXCEN)
        QSUBZ=QSUBZ/SQRT(1-ESQUARE)
        ARGUMENTO=ATAN2(PSUBZ,QSUBZ)

 !      write(*,*)'ANOMU',ANOMU
        ELE(1)=SEMIA
        ELE(2)=EXCEN
        ELE(3)=INCLI
        ELE(5)=NODO
        ELE(4)=ARGUMENTO
        ELE(6)=ANOMALIA
 
 
        RETURN
END subroutine xtoel
!**************************************************************************************
subroutine nr(m,ecc,ea)
!***********************************************
! solves kepler's equation by applying  Newton Raphson Method
!with an accuracy limit of E-14
!
! m[real]...............Mean Anomaly (radian!)
!ecc[real]..............Eccentricity (numerical Eccentricity <1!)
!ea[real]................Eccentric Anomaly (radian!)
!step[int]...............Iteration steps
!
!written by Siegfried Eggl  20061222
!modified                   20111026
!dependencies: none
!**********************************************
implicit none
        integer(kind=idp)::i
        real(kind=dp)::m,ecc,ea,ea0,dea,deat

ea=0._dp
ea0=1.3421_dp
deat=1.d-12
dea=ABS(ea-ea0)

do i=1,100

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._dp-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

end do

if (ea.le.1.d-14) then
   ea=0._dp
end if
!if precision is not achieved try with different initial condition
if(dea>deat) then

ea=0._dp
ea0=0.3421_dp
deat=1.d-12
dea=ABS(ea-ea0)
 do i=1,100

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._dp-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

 end do
end if


if(dea>deat) then

ea=0._dp
ea0=2.3421_dp
deat=1.E-12
dea=ABS(ea-ea0)
 do i=1,100

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._dp-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

 end do
end if

if(dea>deat) then
   write(unit=*,fmt=*)'convergence-error in subroutine nr'
   write(unit=*,fmt=*)'target precision:',deat
   write(unit=*,fmt=*)'achieved precision:',dea
end if
return

end subroutine nr
!*****************************************************************
SUBROUTINE hekoo(rv,n,time)
!------------------------------------------------------
! Transformation ins Heliozentirsche Koordinatensystem
! Die Koordinaten und Geschwindigkeiten  der Sonne werden von denen der massiven Körper
! abgezogen. Anschließend wird die Sonne in den Koordinaten Mittelpunkt gesetzt.
! Würde man die Schleife j bei 1 starten, hätte die Sonne schon nach dem ersten
! Schritt den Nullvektor als Koordinaten, und man bekäme kein brauchbares Ergebnis
!-------------------------------------------------------
implicit none
!IO
integer(kind=idp),intent(in)::n
real(kind=dp),intent(in)::time
real(kind=dp),intent(inout)::rv(1:n,1:6)
!local
integer(kind=idp)::j
integer(kind=idp),parameter::one=1
!real(kind=dp),dimension(1:6)::rvsun
real(kind=dp),dimension(1,1:6)::rvsun

!get ephemeris of the Sun 
call getrvsun(time,rvsun(1,1:6))

!$omp parallel do default(shared)   &
!$omp private(j)
do j=1,n
   rv(j,1:6)=rv(j,1:6)-rvsun(1,1:6)
end do
!$omp end parallel do

!rotate ICRS to heliocentric ecliptic 
!(transformation from barycentric to heliocentric will be done in subroutine hekoo)
 call icrs2he(rv,n)


return
end subroutine hekoo
!*******************************************************
SUBROUTINE bakoo(rv,n,time)
!-----------------------------------------------------
! Using ephemeris barycentric transformation does not quite work the usual way, since not all contributing partilces are present
! We simply assume they knew what they did and the planet data was in fact barycentric at some point.
! Thus,we only perform the inverse heliocentric transformation to make the asteroids koordinates barycentric 
!-----------------------------------------------------
implicit none
!IO
integer(kind=idp),intent(in)::n
real(kind=dp),intent(in)::time
real(kind=dp),intent(inout)::rv(:,:)
!local
integer(kind=idp)::j
integer(kind=idp),parameter::one=1
real(kind=dp),dimension(1,1:6)::rvsun

!get ephemeris 
call getrvsun(time,rvsun(1,1:6))

!Rotate from ICRS to ecliptic
call icrs2he(rvsun,one)

!do barycentric correction

!$omp parallel do default(shared)   &
!$omp private(j)
do j=1,n 
   rv(j,1:6)=rv(j,1:6)+rvsun(1,1:6)
end do
!$omp end parallel do

!rotate all from barycentric ecliptic into ICRS system of Ephemeris
call he2icrs(rv(1:n,1:6),n)

return
end subroutine bakoo
!****************************
subroutine bodykoo(rv,nmin,nmax,cen,time)
!------------------------------------------------------
! Transformation to bodycentric coordinate frame
! rv....... barycentric distance and velocity vectors
! nmin..... start id number of bodies to transform
! nmax....  end id number of bodies to transform
! cen....   id number of body center
! time...   current ephemeris time
!-------------------------------------------------------
implicit none
!IO
integer::cen
integer(kind=idp),intent(in)::nmin,nmax
real(kind=dp),intent(in)::time
real(kind=dp),intent(inout)::rv(:,:)
!local
integer(kind=idp)::j
real(kind=dp),dimension(1:6)::rvcen

!get ephemeris of the Sun 
call getrvbody(time,cen,rvcen)
!$omp parallel do default(shared)   &
!$omp private(j)
do j=nmin,nmax
   rv(j,1:6)=rv(j,1:6)-rvcen(1:6)
end do
!$omp end parallel do

return
end subroutine bodykoo

!****************************
subroutine bodykoo1(rv,cen,time)
!------------------------------------------------------
! Transformation to bodycentric coordinate frame
! rv....... barycentric distance and velocity vectors
! cen....   id number of body center
! time...   current ephemeris time
!-------------------------------------------------------
implicit none
!IO
integer,intent(in)::cen
real(kind=dp),intent(in)::time
real(kind=dp),intent(inout)::rv(1:6)
!local
integer(kind=idp)::j
real(kind=dp),dimension(1:6)::rvcen

!get ephemeris for the body with JPL ID cen



call getrvbody(time,cen,rvcen)
   
   rv(1:6)=rv(1:6)-rvcen(1:6)

return
end subroutine bodykoo1



! SUBROUTINE bakoo(rv,time)
! !-----------------------------------------------------
! ! Transformation ins Baryzentrische Koordinatensystem
! !-----------------------------------------------------
! implicit none
! integer(kind=idp)::j
! real(kind=dp),intent(in)::time
! real(kind=dp)::rv(:,:)
! real(kind=dp)::rb(1:6),rbx,rby,rbz,rbvx,rbvy,rbvz,summass
! real(kind=dp),dimension(1:neph,1:6)::rveph
! !db
! real(kind=dp)::rb2(1:6)
! !edb
! 
! !get ephemeris 
! call getrv(time,rveph)
! 
! !sorry, but open mp reduction just works for scalars...
! rbx=0._dp
! rby=0._dp
! rbz=0._dp
! rbvx=0._dp
! rbvy=0._dp
! rbvz=0._dp
! 
! !$omp parallel do default(shared)   &
! !$omp private(j) &
! !$omp reduction(+:rbx,rby,rbz,rbvx,rbvy,rbvz)
! do j=1,neph
!    rbx=rbx+meph(j)*rveph(j,1)
!    rby=rby+meph(j)*rveph(j,2)
!    rbz=rbz+meph(j)*rveph(j,3)
!    rbvx=rbvx+meph(j)*rveph(j,4)
!    rbvy=rbvy+meph(j)*rveph(j,5)
!    rbvz=rbvz+meph(j)*rveph(j,6)
! end do
! !$omp end parallel do
! 
! !db
! summass=Sum(meph)
! rb2(1)=rbx/summass
! rb2(2)=rby/summass
! rb2(3)=rbz/summass
! rb2(4)=rbvx/summass
! rb2(5)=rbvy/summass
! rb2(6)=rbvz/summass
! !edb
! 
! !$omp parallel do default(shared)   &
! !$omp private(j) &
! !$omp reduction(+:rbx,rby,rbz,rbvx,rbvy,rbvz) 
! do j=1,inbody
!    rbx=rbx+mass(j)*rv(j,1)
!    rby=rby+mass(j)*rv(j,2)
!    rbz=rbz+mass(j)*rv(j,3)
!    rbvx=rbvx+mass(j)*rv(j,4)
!    rbvy=rbvy+mass(j)*rv(j,5)
!    rbvz=rbvz+mass(j)*rv(j,6)
! end do
! !$omp end parallel do 
! 
! summass=Sum(mass)+Sum(meph)
! 
! rb(1)=rbx/summass
! rb(2)=rby/summass
! rb(3)=rbz/summass
! rb(4)=rbvx/summass
! rb(5)=rbvy/summass
! rb(6)=rbvz/summass
! 
! !db
! write(*,*)'delta rb',rb-rb2
! write(*,*)'rb',rb
! !edb
! 
! !$omp parallel do default(NONE)   &
! !$omp shared(rv,inbody,rb) private(j)
! do j=1,inbody
!    rv(j,1:6)=rv(j,1:6)-rb(1:6)
! end do
! !$omp end parallel do
! 
! return
! end subroutine bakoo
! !*****************************************************************


!**********************************************************
! subroutine jd2greg(jd,a,mon,day)
! !calculates Gregorian year, month, day from JD following Richards
! implicit none
! integer,intent(out)::a,mon,day
! real(kind=dp),intent(in)::jd
! 
! integer,parameter::y=4716,n=12,m=2
! real(kind=dp),parameter::j=1401,r=4,p=1461
! real(kind=dp),parameter::v=3,u=5,s=153,w=2,B=274277,C=-38
! real(kind=dp)::f,e,g,h
! 
! f=jd+j
! f=f+real(int((int((4.d0*jd+B)/146097.d0)*3.d0)/4.d0))+C
! e=r*f+v
! g=real(int(mod(e,p)/r))
! h=u*g+w
! day=int(mod(h,s)/u)+1
! mon=mod(int(real(h)/real(s))+m,n)+1
! a=int(e/p)-y+int(real(n+m-mon)/real(n))
! 
! return
! end subroutine 


SUBROUTINE jd2greg (JDr, YEAR,MONTH,DAY)

!---COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY)
!   GIVEN THE JULIAN DATE (JD).
!
implicit none
    real(kind=dp),intent(in)::JDr
    INTEGER,intent(out):: YEAR,MONTH,DAY
    INTEGER::L,N,I,J,K,JD
    
    JD=nint(JDr)
!
    L= JD+68569
    N= 4*L/146097
    L= L-(146097*N+3)/4
    I= 4000*(L+1)/1461001
    L= L-1461*I/4+31
    J= 80*L/2447
    K= L-2447*J/80
    L= J/11
    J= J+2-12*L
    I= 100*(N-49)+I+L
!
    YEAR= I
    MONTH= J
    DAY= K
!
  RETURN
END SUBROUTINE

!*************************************************************    

        subroutine topolar2(x,r)
!----------------------------------------------------
!       Umwandlung von kartesischen in 
!       Polarkoordinaten -
!       entgegen der üblichen Konvention normiert
!       wird hier theta auf die Ekliptik normiert
!--------------------------------------------------
          implicit none
          real(kind=dp),intent(in)::x(1:3) !cartesian coordinates (x,y,z)
          real(kind=dp),intent(out)::r(1:3) !polar coordinates (r,phi,theta)
        
          r(1)=sqrt(Dot_Product(x,x)) 	!r
          r(2)=atan2(x(2),x(1)) 	!phi
          r(3)=asin(x(3)/r(1)) 		!theta
     
          return
        end subroutine

!*************************************************************        
        subroutine tocartesian2(r,x)
!----------------------------------------------------
!       Umwandlung von Polarkoordinaten in 
!       kartesische Koordinaten -
!       entgegen der üblichen Konvention normiert
!       wird hier theta auf die Ekliptik normiert
!--------------------------------------------------
        implicit none
        real(kind=dp),intent(in)::r(1:3)      !polar coordinates (r,phi,theta)
        real(kind=dp),intent(out)::x(1:3)   !cartesian coordinates (x,y,z)

        x(1) = r(1)*cos(r(3))*cos(r(2))  !x
        x(2) = r(1)*cos(r(3))*sin(r(2))  !y
        x(3) = r(1)*sin(r(3))            !z

        return
        end subroutine
!**********************************************************************
!       subroutine gs(kS1,kS2,c1,c2,c3,PI)
! !----------------------------------------------------------------
! !         gs steht fuer Gram-Schmidtsches Orthonomalisierungs-
! !         Verfahren
! !         es werden 2 verschiedene gs-Verfahren angewendet,  
! !         da fuer die Winkel alpha und beta jeweils Rechts-
! !         und Linkssysteme verwendet werden
! !---------------------------------------------------------------
!   
!        implicit none
!        real*8:: kS1(1:3),kS2(1:3),b1(1:3),b2(1:3),b3(1:3), &
!      c1(1:3),c2(1:3),c3(1:3),betr,sk1,sk2,sk3,skcheck(1:3),&
!      polr,polphi,poltheta,PI
!        integer::i
!       
! !------------------------------ c1 -----------------------------
! !         entspricht dem normierten Vektor b1 (S1S2)
! !
! !         die restilichen Basisvektoren b2 & b3 eigentlich 
! !         frei waehlbar, nur linear unabhaengig
! !         aber hier orthogonal aufeinander gewaehlt
! !---------------------------------------------------------------        
! 
!           b1(:)=(kS2(:)-kS1(:))
!       
!        call betrag(b1,betr)
! 
!           c1(:)=b1(:)/betr
!        
! !------------------------------ c2 -----------------------------
! !         nehme Vektor b2 mit  
! !         gleichen x,z Koordinaten mit auf c1 orthogonalem y 
! !---------------------------------------------------------------
!    
! !-----------------------------------------------------------
! !         Umwandlung in Polarkoordinaten, 
! !         wähle b2 im Rechten Winkel zu b1
! !        -> + 90°
! !----------------------------------------------------------
!           call topolar(b1(:),polr,polphi,poltheta)
!           
!           polphi=polphi+PI/2.d0
! !---------------------------------------------------------
! !        Rücktransformation in kartesische Koordinaten
! !--------------------------------------------------------
!           call tocartesian(polr,poltheta,polphi,b2(:))
! 
! 
! !---------------------------------------------------------------
! !         Gram Schmidt fuer Vektor c2
! !---------------------------------------------------------------
!        call skalarproduct(b2,c1,sk1)
!       
!      
!           c2(:)=b2(:)-sk1*c1(:)
!                 
!       call betrag(c2,betr)
! 
!           c2(:)=c2(:)/betr
! 
! !------------------------------ c3 ----------------------------- 
! !         b3 Kreuzprodukt von b1 mit b2
! !         die Orientierung bleibt nach oben
! !         -> Rechtssystem 
! !---------------------------------------------------------------
! 
!        b3(1)=(b1(2)*b2(3)-b1(3)*b2(2))
!        b3(2)=(b1(3)*b2(1)-b1(1)*b2(3))
!        b3(3)=(b1(1)*b2(2)-b1(2)*b2(1))
! 
! !---------------------------------------------------------------  
! !         Gram Schmidt fuer Vektor c3
! !---------------------------------------------------------------
! 
!        call skalarproduct(b3,c2,sk2)
!        call skalarproduct(b3,c1,sk3)
!    
!     
!        c3(:)=b3(:)-(sk2)*c2(:)-(sk3)*c1(:)
!   
!        call  betrag(c3,betr)
!        
!           c3(:)=c3(:)/betr
!        
!        
!        call skalarproduct(c1,c2,skcheck(1))
!        call skalarproduct(c1,c3,skcheck(2))
!        call skalarproduct(c2,c3,skcheck(3))
!   
! 
!        do i=1,3
!           if(abs(skcheck(i))>1D-12) then
!              write(*,*)'Fehler in gsalpha: '
!              write(*,*)'neue Basisvektoren nicht orthogonal',skcheck(:)
!           end if
!        end do
!        
!        return
! 
!      end subroutine gs
!////////////////////////////////////////////////////////////////
end module transform_m
