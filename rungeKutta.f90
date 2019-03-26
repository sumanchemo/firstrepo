! This program is to solve the differential equations by using different
! methods like Euler, Heun and RungeKutta.

program runga
  implicit none
  integer,parameter::nt=25,n=2    ! nt number of grid points, n is number of variable in equation.
  integer i
  real*8 x(n),y(n),z(n),dt,t1,t2,t
  
  t1=0.0; t2=2.0; dt=(t2-t1)/nt

  x(1)=0; x(2)=0  ! runga kutta
  y(1)=0; y(2)=0  ! Euler 
  z(1)=0; z(2)=0  ! Heun
  
  do i=0,nt
     t=i*dt
     write(7,*)t,x(1),x(2),y(1),y(2)
     call rk4(x,n,dt,t)
     call Euler(y,n,dt,t)
     call Heun(z,n,dt,t)
  end do                                                          
    
  
  
end program runga


subroutine Euler(y,n,dt,t)    ! Euler method Subroutine 
  implicit none
  integer n
  real*8 y(n),t,dt,f(n)
  call ff(y,n,t,f)
  y=y+f*dt

end subroutine Euler


subroutine Heun(y,n,dt,t)    ! Heun method Subroutine
  implicit none
  integer n
  real*8 y(n),t,dt,f(n),x(n),fx(n)
  call ff(y,n,t,f)
  x=y
  x=x+f*dt
  call ff(x,n,t,fx)
  y=y+(f+fx)*dt/2

end subroutine Heun


 ! RungeKutta Subroutines  RK1, RK2 and RK4 *******
subroutine ff(y,n,t,f)
  implicit none
  integer i,n
  real*8 y(n),t,f(n)
  
  f(1)=y(1)*y(2)+cos(t)-0.5*sin(2*t)
  f(2)=y(1)**2+y(2)**2-(1+sin(t))   ![ only change this for any DE ]
  
end subroutine ff


subroutine rk1(y,n,dt,t)           ! RK1  Runge-Kutta method  
  implicit none
  integer i,n
  real*8 y(n), dt, f(n),t

  call ff(y,n,t,f)
  y=y+dt*f
end subroutine rk1


subroutine rk2(y,n,dt,t)           ! RK2  Runge-Kutta method  
  implicit none
  integer i,n
  real*8 y(n), dt, f(n),t,y_var(n),t_var

  call ff(y,n,t,f)
  y_var=y+0.5*dt*f
  t_var=t+0.5*dt
  call ff(y_var, n, t_var, f)
  y=y+dt*f
  
end subroutine rk2

subroutine rk4(y,n,dt,t)           ! RK4  Runge-Kutta method  
  implicit none
  integer i,n
  real*8 y(n), dt, f(n), t, k1(n), k2(n), k3(n), k4(n), y_var(n), t_var

  call ff(y,n,t,f)
  k1=dt*f

  y_var=y+0.5*k1
  t_var=t+0.5*dt
  call ff(y_var,n,t_var,f)

  k2=dt*f
  y_var=y+0.5*k2
  t_var=t+0.5*dt
  call ff(y_var,n,t_var,f)

  k3=dt*f
  y_var=y+k3
  t_var=t+dt
  call ff(y_var,n,t_var,f)

  k4=dt*f

  y=y+(k1+2*k2+2*k3+k4)/6
   
end subroutine rk4



