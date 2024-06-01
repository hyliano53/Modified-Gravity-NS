program TOV

implicit none 
      
! Constant variables
double precision pi,mpi4tokmm2,kmpersunsmass,cspeed,rs,ms,Nu_inf

! Variables for the Runge-Kutta iteration:
double precision ka1,ka2,ka3,ka4,kb1,kb2,kb3,kb4,kc1,kc2,kc3,kc4,h,kd1,kd2,kd3,kd4,ke1,ke2,ke3,ke4
double precision kf1,kf2,kf3,kf4,kg1,kg2,kg3,kg4
double precision alpha, Omega, V, derV, phi_a, phi_b, err, prec, sig, Omega_vel,Jtot

! Arrays to contain the dynamical functions (pressure, mass, density, r)
integer i,j,N,Nit,len,s,pos,k
parameter(N=13900)
parameter(Nit=500)
double precision lambda(N),Edens(N),P(N),r(N),Nu(N),derNu(N),phi(N),g(N),w(N),a(N),derlambda(N)
      
! Function for the equation of state
double precision Edens_interpol,Edensaux,badcs

external Edens_interpol
common/unitconversion/mpi4tokmm2
common/pinumber/pi
pi=3.14159265d0

! A file_to_save the functions
open(unit=7,file='Mass_of_Radius.dat')
open(unit=8,file='J_moment.dat')


! Unit conversion from MeV/fm3 to km**-2
mpi4tokmm2=1.3152d-6

! Unit coversion from s to km
cspeed=2.99792d5

! Unit conversion from km to solar masses
kmpersunsmass=1.47d0

! The energy density for which the speed of sound exceeds the light speed
! in the effective theory. For the free neutron gas this does not occur.
badcs=1d20*mpi4tokmm2

! Prepare a table rho(P) from a file
call init_eqofstate()

! Initialize the radial variable
h=0.005
r(1)=0.005
do i=2,N
  if(i.eq.4001) then
    h=0.016
  endif
  r(i)=r(i-1)+h
enddo

!Initial conditions
lambda(1)=0d0
g(1)=0d0
Nu(1)=1d0


! Main iteration: for each value of the central pressure we compute
! the mass and radius of the star
do j=1,Nit

!Interval that cotains the initial condition desired (for shoothing method)
phi_a=0d0
phi_b=1d0

rs=0
s=0
err=1d0
prec=5d-4

! The typical cetral pressure of a NS is about hundreds of MeV/fm^3. We need it in km**-2
P(1)=0d0*mpi4tokmm2+5*dble(j)*mpi4tokmm2
Edens(1)=Edens_interpol(P(1))

do while(prec.le.err)
  h=0.005
      s=s+1
      phi(1)=(phi_a+phi_b)/2d0

! Auxiliary iteration where we integrate out in r
	do i=2,N

    if(i.eq.4001) then
      h=0.016
    endif

		call auxders(r(i-1),Edens(i-1),lambda(i-1),Nu(i-1),P(i-1),g(i-1),phi(i-1),ka1,kb1,kc1,kd1,ke1)
            Edensaux=Edens_interpol(P(i-1)+h*kc1/2d0)
! Construct the Runge-Kutta derivatives
        call auxders(r(i-1)+h/2d0, Edensaux, lambda(i-1)+h/2d0*ka1, Nu(i-1)+h/2d0*kb1, &
        P(i-1)+h/2d0*kc1, g(i-1)+h/2d0*ke1, phi(i-1)+h/2d0*kd1, ka2, kb2, kc2, kd2, ke2)
        Edensaux=Edens_interpol(P(i-1)+h*kc2/2d0)
        
        call auxders(r(i-1)+h/2d0, Edensaux, lambda(i-1)+h/2d0*ka2, Nu(i-1)+h/2d0*kb2, & 
        P(i-1)+h/2d0*kc2, g(i-1)+h/2d0*ke2, phi(i-1)+h/2d0*kd2, ka3, kb3, kc3, kd3, ke3)
        Edensaux=Edens_interpol(P(i-1)+h*kc3)
        
        call auxders(r(i-1)+h, Edensaux, lambda(i-1)+h*ka3,  Nu(i-1)+h*kb3, &
        P(i-1)+h*kc3, g(i-1)+h*ke3, phi(i-1)+h*kd3, ka4, kb4, kc4, kd4, ke4)
! advance_ the pressure with Runge-Kutta
        P(i)=P(i-1)+(h/6d0)*(kc1+2d0*kc2+2d0*kc3+kc4)
! advance_the mass with Runge-Kutta
        lambda(i)=lambda(i-1)+(h/6d0)*(ka1+2d0*ka2+2d0*ka3+ka4)
! advance_the Nu fuction with Runge-Kutta
        Nu(i)=Nu(i-1)+(h/6d0)*(kb1+2d0*kb2+2d0*kb3+kb4)
! advance_the g fuction with Runge-Kutta
        g(i)=g(i-1)+(h/6d0)*(ke1+2d0*ke2+2d0*ke3+ke4)
! advance_the phi fuction with Runge-Kutta
        phi(i)=phi(i-1)+(h/6d0)*(kd1+2d0*kd2+2d0*kd3+kd4)

        Edens(i)=Edens_interpol(P(i))

        if(isnan(phi(i)).eqv..false.) then
            sig=phi(i-1)-phi(i-2)
            !print*, sig
        endif
	enddo
	
      err=abs(phi(i-1))

      if((isnan(err).eqv..true.).and.(sig.le.0d0)) then
            phi_a=phi(1)
            err=1d0
      elseif((isnan(err).eqv..true.).and.(sig.gt.0d0)) then
            phi_b=phi(1)
            err=1d0
      endif


      if((isnan(err).eqv..false.).and.(phi(i-1).le.0d0)) then
          phi_a=phi(1)
      elseif((isnan(err).eqv..false.).and.(phi(i-1).gt.0d0)) then
          phi_b=phi(1)
      endif


 if(s==60) then
      goto 150
 endif

enddo


150 continue

h=0.005
k=1d6
!We solve the system with the desired initial conditions
do i=2,N

  if(i.eq.4001) then
    h=0.016
  endif

  if(i.ge.k) then
    phi(i-1)=0d0
    g(i-1)=0d0
  endif

      call auxders(r(i-1),Edens(i-1),lambda(i-1),Nu(i-1),P(i-1),g(i-1),phi(i-1),ka1,kb1,kc1,kd1,ke1)
      Edensaux=Edens_interpol(P(i-1)+h*kc1/2d0)
      derNu=kb1
      derlambda=ka1
! Construct the Runge-Kutta derivatives
  call auxders(r(i-1)+h/2d0, Edensaux, lambda(i-1)+h/2d0*ka1, Nu(i-1)+h/2d0*kb1, &
  P(i-1)+h/2d0*kc1, g(i-1)+h/2d0*ke1, phi(i-1)+h/2d0*kd1, ka2, kb2, kc2, kd2, ke2)
  Edensaux=Edens_interpol(P(i-1)+h*kc2/2d0)
  
  call auxders(r(i-1)+h/2d0, Edensaux, lambda(i-1)+h/2d0*ka2, Nu(i-1)+h/2d0*kb2, & 
  P(i-1)+h/2d0*kc2, g(i-1)+h/2d0*ke2, phi(i-1)+h/2d0*kd2, ka3, kb3, kc3, kd3, ke3)
  Edensaux=Edens_interpol(P(i-1)+h*kc3)
  
  call auxders(r(i-1)+h, Edensaux, lambda(i-1)+h*ka3,  Nu(i-1)+h*kb3, &
  P(i-1)+h*kc3, g(i-1)+h*ke3, phi(i-1)+h*kd3, ka4, kb4, kc4, kd4, ke4)
! advance_ the pressure with Runge-Kutta
  P(i)=P(i-1)+(h/6d0)*(kc1+2d0*kc2+2d0*kc3+kc4)
! advance_the mass with Runge-Kutta
  lambda(i)=lambda(i-1)+(h/6d0)*(ka1+2d0*ka2+2d0*ka3+ka4)
! advance_the Nu fuction with Runge-Kutta
  Nu(i)=Nu(i-1)+(h/6d0)*(kb1+2d0*kb2+2d0*kb3+kb4)
  !derNu(i)=kb1
! advance_the g fuction with Runge-Kutta
  g(i)=g(i-1)+(h/6d0)*(ke1+2d0*ke2+2d0*ke3+ke4)
! advance_the phi fuction with Runge-Kutta
  phi(i)=phi(i-1)+(h/6d0)*(kd1+2d0*kd2+2d0*kd3+kd4)


  Edens(i)=Edens_interpol(P(i))

! Conditions to break from the loop: reached the_end of the star
! or reached an energy density for which causality breaks down
  if((P(i).le.5d-10).and.(rs.eq.0d0)) then
    rs=r(i)*dexp(-phi(i)/dsqrt(3d0))
    pos=i
    k=2*i
  endif
enddo


!Now, we reescale the nu function
Nu(N)=Nu_inf
do i=1,N
  Nu(i)=Nu(i)-Nu_inf
enddo

!Now, we compute the angular velocity
w(1)=0.004
a(1)=0d0
h=0.005
do i=2,N

  if(i.eq.4001) then
    h=0.016
  endif

  call auxders2(r(i-1),Nu(i-1),lambda(i-1),Edens(i-1),P(i-1),phi(i-1),derNu(i-1), &
  derlambda(i-1),a(i-1),w(i-1),kf1,kg1)
  call auxders2(r(i-1)+h/2d0,Nu(i-1),lambda(i-1),Edens(i-1),P(i-1),phi(i-1), &
  derNu(i-1),derlambda(i-1),a(i-1)+h/2d0*kf2,w(i-1)+h/2d0*kg2,kf2,kg2)
  call auxders2(r(i-1)+h/2d0,Nu(i-1),lambda(i-1),Edens(i-1),P(i-1),phi(i-1), &
  derNu(i-1),derlambda(i-1),a(i-1)+h/2d0*kf3,w(i-1)+h/2d0*kg3,kf3,kg3)
  call auxders2(r(i-1)+h,Nu(i-1),lambda(i-1),Edens(i-1),P(i-1),phi(i-1), &
  derNu(i-1),derlambda(i-1),a(i-1)+h*kf4,w(i-1)+h*kg4,kf4,kg4)

  a(i)=a(i-1)+(h/6d0)*(kf1+2d0*kf2+2d0*kf3+kf4)
  w(i)=w(i-1)+(h/6d0)*(kg1+2d0*kg2+2d0*kg3+kg4)
enddo

!The angular velocity of the star as seen by an observer at infinity
Omega_vel=w(N)

!Reescale the angular velocity
do i=1,N
  w(i)=w(i)*2000d0/(Omega_vel*cspeed)
  a(i)=a(i)*2000d0/(Omega_vel*cspeed)
enddo

Omega_vel=w(N)

!Angular momentum
Jtot=(dexp(-lambda(pos)-Nu(pos))*r(pos)**(4)*a(pos))/6d0
write(8,*) r(i-1)*(1d0-dexp(-2d0*lambda(i-1)))/(2d0*kmpersunsmass), Jtot/(Omega_vel*kmpersunsmass)

!150 continue

print*, j
write(7,*) rs, r(i-1)*(1d0-dexp(-2d0*lambda(i-1)))/(2d0*kmpersunsmass)

enddo

close(7)
close(8)
end




! SUBROUTINES
subroutine auxders(r,Edens,lambda,Nu,P,g,phi,derlambda,derNu,derP,derphi,derg)
!Calculates the derivatives of the desired functions
      implicit none
      double precision pi,r,Edens,lambda,Nu,P,g,phi,derlambda,derP,derNu,derg,derphi
      double precision Omega, V, derV, alpha
      common/pinumber/pi
      alpha=1d1
      Omega=dexp(phi/dsqrt(3d0))
      V=(1d0-Omega**(-2))**2/(4d0*alpha)
      derV=(1d0-Omega**(-2))*Omega**(-2)/(dsqrt(3d0)*alpha)

      derlambda=4d0*pi*r*Edens*Omega**(-4)*dexp(2d0*lambda)+(1d0-dexp(2d0*lambda))/(2d0*r)+ &
      r*g**2/2d0+r*dexp(2d0*lambda)*V/4d0
      derNu=4d0*pi*r*P*Omega**(-4)*dexp(2d0*lambda)+(dexp(2d0*lambda)-1d0)/(2d0*r)+ &
      r*g**2/2d0-r*dexp(2d0*lambda)*V/4d0
      derP=-(Edens+P)*(derNu-g/dsqrt(3d0))
      derphi=g
      derg=-(derNu-derlambda+2d0/r)*g+4d0*pi*dexp(2d0*lambda)*Omega**(-4)*(3d0*P-Edens)/dsqrt(3d0)+ &
      dexp(2d0*lambda)*derV/4d0

      return
      end

subroutine auxders2(r,Nu,lambda,Edens,P,phi,derNu,derlambda,a,w,dera,derw)
!Calculates the derivatives of the desired functions
      implicit none
      double precision pi,r,Edens,lambda,Nu,P,phi,derlambda,derP,derNu,dera,derw,a,w
      double precision Omega, V, derV, alpha
      common/pinumber/pi
      Omega=dexp(phi/dsqrt(3d0))

      dera=16d0*pi*Omega**(-4)*(Edens+P)*dexp(2d0*lambda)*w+(derNu+derlambda-4d0/r)*a
      derw=a
      return
      end


subroutine init_eqofstate()
      implicit none
      integer i,N
      parameter(N=121)
      double precision rhoinit(N),Pinit(N),mpi4tokmm2
      common/eqstate/rhoinit,Pinit
      common/unitconversion/mpi4tokmm2
      open(unit=5,file='P_of_rho.dat')

      do i=1,N
        read(5,*) rhoinit(i),Pinit(i)  
        !Convert to the correct units, from MeV/fm3 to km**-2
        rhoinit(i)=rhoinit(i)*mpi4tokmm2
        Pinit(i)=Pinit(i)*mpi4tokmm2
      enddo
      close(5)
      return
      end



double precision function Edens_interpol(P)
! This is to obtain the energy density from the pressure for a loaded
! equation of state in a table
      implicit none
      double precision rho,P,mpi4tokmm2
      integer N,j,jstart,jup,jm
      parameter(N=121)
      
      double precision rhoinit(N),Pinit(N)
      common/eqstate/rhoinit,Pinit
      common/unitconversion/mpi4tokmm2
      if(P.LT.Pinit(1)) then
! Linear interpolation between the last points
       rho=rhoinit(1)/Pinit(1)*P
       goto 200
      endif
      if(P.GT.Pinit(N)) then
! Relativistic equation of state for extrapolating
       rho=3d0*P
       goto 200
      endif
      jstart=1
      jup=N
 270  if((jup-jstart).gt.1) then
        jm=(jup+jstart)/2
        if(P.GT.Pinit(jm)) then
          jstart=jm
        else
          jup=jm
        endif   
        goto 270
      endif   
      j=jstart
      rho=rhoinit(j)+(P-Pinit(j))*(rhoinit(j+1)-rhoinit(j))/(Pinit(j+1)-Pinit(j))
 200  Edens_interpol=rho
      return 
      end


double precision function Const_density(P)
      implicit none
      double precision P
! This equation of state is just a constant density for neutron matter
! We employ the Gaussian system of units, so that energy density is erg/cm**3
      Const_density=1E15
      return
      end
