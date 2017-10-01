program pme

  !*********************************************************************************
  !**                                                                             **
  !**     Program to solve the 1-D porous media equation using moving finite      **
  !**     elements.  The porous media equation is given as:-                      **
  !**                                                                             **
  !**                             u_t + (u^m u_x)_x =0                            **
  !**     Variables :-                                                            **
  !**                                                                             **
  !**             u (solution) - Double precision array                           **
  !**             x (mesh position)- Double precision array                       **
  !**             x_dot (mesh velocity) - Double precision array                  **
  !**                                                                             **                      
  !**             nodes - number of nodes in mesh. Integer                        **
  !**             m - power in equation. Double precision                         **
  !**             output_t - Time of solution output. Double precision            **
  !**                                                                             **
  !**             delta_t - Size of time-step. Double precision                   **
  !**             total_t - total time of numerical simulation. Double precision  **
  !**                                                                             **
  !**                                                                             **
  !*********************************************************************************
  implicit none
  !---------------------------------------------------------------------------------
  integer :: nodes, i, reports
  double precision, dimension(:), allocatable :: u, x, theta, x_dot
  double precision :: output_t, delta_t, total_t, t_init
  double precision :: report_step, report_time
  double precision :: Q, m

  ! Computational timing variables
  integer :: System_Time_Start,System_Time_Stop,System_Time_Rate
  real :: CPU_Time_Start,CPU_Time_Stop
  
  logical :: writesol
  character(LEN=10)::G
  parameter (G='(110F16.9)')
  !---------------------------------------------------------------------------------  

  ! Print simulator banner.
  write(6,*)
  write(6,*) '---------------------------------------------------------------------'
  write(6,*) '                      Porous Media Equation'
  write(6,*) '---------------------------------------------------------------------'
  write(6,*)

  ! Timing routine for calculation of time of numerical simulation.

   ! Start timing procedure
  call system_clock(System_Time_Start,System_Time_Rate)
  call cpu_time(CPU_Time_Start)

  ! Read variables into program.   
  open(unit=10,file='variables.data',status='old',form='formatted')
  read(10,*) nodes
  read(10,*) output_t
  read(10,*) m
  read(10,*) Q
  read(10,*) t_init     
  read(10,*) reports
  close(10)

  ! Allocate arrays.
  allocate(u(0:nodes),x(0:nodes),x_dot(0:nodes),theta(0:nodes))
  u=0.0d0; x=0.0d0; x_dot=0.0d0; theta=0.0d0

  ! Setup the initial conditions for the simulation
  call initial_conditions(u,x,nodes,m,Q,t_init)

  ! Calculate the theta vector when mesh not equidistributed                                          
  call theta_setup(theta,nodes,x,u)

  ! Setup the equispaced report times for the numerical simulation.
  total_t = t_init
  report_step = (output_t-t_init)/reports
  report_time = t_init +  report_step
  writesol=.false.                     

  ! Open the solution output file.
  open(unit=20,file='solution.out')
  write(20,G) u(:), total_t
  write(20,G) x(:), total_t

  ! Write time-stepping banner
  write(6,*)
  write(6,*) 'Advancing solution to time ', output_t
  write(6,*)
  write(6,*) '---------------------------------------------------------------------'
  write(6,*) '         Time                Time-step                               '
  write(6,*) '---------------------------------------------------------------------'


  do while (total_t<output_t)                                 

     ! Calculate the time-step using an adaptive time-stepping routine
      call adaptive_timestep(u,x,m,delta_t,total_t,nodes)

      if (total_t+delta_t>report_time) then
        delta_t=report_time-total_t 
        report_time = report_time +  report_step
        writesol=.true.
      endif

     if (total_t+delta_t>output_t) then
        delta_t=output_t-total_t
     endif
     
     total_t = total_t + delta_t 

     ! Write out the current time and time-step size
     write(6,'(2f20.8)') total_t, delta_t

     ! Mass monitor
     call mass_mesh_velocity(u,x,x_dot,nodes,delta_t,m)

     ! Recover the solution through conservation principle.
     call mass_u_calc(u,x,nodes,theta) 

     if (writesol) then
        write(20,G) u(:), total_t
        write(20,G) x(:), total_t 
        writesol=.false.
     endif
  end do

  close(20)

  ! Stop the timing
  call system_clock(System_Time_Stop)
  call cpu_time(CPU_Time_Stop)

  write(*,'(/)')
  write(*,'(1x,a)') 'Timing report:'
  write(6,*) '---------------------------------------------------------------------'
  write(*,'(1x,a)') '    Elapsed Time   CPU Time',&
                    '        (s)           (s)'
  write(6,*) '---------------------------------------------------------------------'
  write(*,'(1x,2e15.4)') dble(System_Time_Stop-System_Time_Start)/dble(System_Time_Rate)&
                       &,CPU_Time_Stop-CPU_Time_Start
  write(6,*) '---------------------------------------------------------------------'
  write(*,'(/)')

end program pme

subroutine initial_conditions(u,x,nodes,m,Q,t_init)

  use special_functions
  !*********************************************************************************
  !**                                                                             **
  !**  This function generates a equally spaced mesh to solve the porous medium   **
  !**  equation on. The solution is initialised to the self-similar solution      **
  !**  given by:                                                                  **
  !**                                                                             **
  !**       u(x,t) = (1/lambda)*(1-(x/(rzero*lambda)))^(1/m)                      **
  !**                                                                             **
  !**  where                                                                      **
  !**                                                                             **
  !**     rzero = Q*gamma(1/m + 3/2)/(SQRT(pi)*gamma(1/m + 1))                    **
  !**     tzero = ((rzero^2)*m)/(2*(m+2))                                         **
  !**     lambda = (tinit/tzero)^(1/(m+2))                                        **
  !**                                                                             **
  !*********************************************************************************
  implicit none
  !---------------------------------------------------------------------------------
  integer, intent(IN) :: nodes
  double precision, intent(INOUT), dimension(0:nodes) :: u, x
  double precision, intent(IN) :: Q, t_init, m
  !---------------------------------------------------------------------------------  
  double precision :: rzero, tzero, gamman, gammad, pi, lambda, delta_x, r
  double precision :: one_over_m, one_over_lambda
  integer :: i
  !---------------------------------------------------------------------------------
  rzero=0; tzero=0; gamman=0; gammad=0; lambda=0
  u=0; x=0; pi=dacos(-1.0d0); delta_x=0; r=0.0

  one_over_m = 1.0d0/m

  ! Calculate quantities for self-similar solution.     
  gamman = gamma(one_over_m+(3.0d0/2.0d0))
  gammad = gamma(one_over_m+1.0d0)

  rzero = Q*gamman/(SQRT(pi)*gammad)
  tzero = ((rzero**2)*m)/(2.0d0*(m+2.0d0))
  lambda = (t_init/tzero)**(1.0d0/(m+2.0d0))

  r = rzero*lambda
  one_over_lambda = 1.0d0/lambda
  delta_x = 2.0d0 * r/nodes

  ! Use equispaced mesh
  do i=0, nodes
     x(i) = -r + delta_x*i
     u(i) = one_over_lambda*(dabs(1.0d0 - ((x(i)/r)**2)))**one_over_m
  end do

  return

end subroutine initial_conditions

subroutine theta_setup(theta,nodes,x,u)

  implicit none
  !------------------------------------------------------------------------------   
  integer, intent(IN) :: nodes
  double precision, intent(INOUT), dimension(0:nodes) :: theta
  double precision, intent(IN), dimension(0:nodes) :: u, x
  !------------------------------------------------------------------------------
  integer :: i
  !------------------------------------------------------------------------------   

  ! Calculate theta for mass monitor
  theta(0) = (1.0d0/3.0d0)*(x(1)-x(0)) * u(0) + &
             (1.0d0/6.0d0)*(x(1)-x(0)) * u(1)

  do i=1, nodes-1
     theta(i) = (1.0d0/6.0d0)*(x(i)-x(i-1))   * u(i-1) + &
                (1.0d0/3.0d0)*(x(i+1)-x(i-1)) * u(i)   + &
                (1.0d0/6.0d0)*(x(i+1)-x(i))   * u(i+1)
  end do

  theta(nodes) = (1.0d0/3.0d0)*(x(nodes)-x(nodes-1)) * u(nodes) + &
                 (1.0d0/6.0d0)*(x(nodes)-x(nodes-1)) * u(nodes-1)
  return

end subroutine theta_setup

subroutine mass_mesh_velocity(u,x,x_dot,nodes,delta_t,m)
  use linear_solvers

  !*********************************************************************************
  !**                                                                             **
  !**     This subroutine solves the matrix system:-                              **
  !**                                                                             **
  !**                     K x_dot = E                                             **
  !**                                                                             **
  !*********************************************************************************
  implicit none
  !---------------------------------------------------------------------------------
  integer, intent(IN) :: nodes
  double precision, intent(INOUT), dimension(0:nodes) :: u, x, x_dot
  double precision, intent(IN) :: delta_t, m         
  !---------------------------------------------------------------------------------
  double precision, dimension(0:nodes,0:nodes) :: K
  double precision, dimension(0:nodes) :: psi, E
  integer :: i
  !--------------------------------------------------------------------------------
  x_dot=0.0d0; K=0.0d0; E=0.0d0; psi=0.0d0; 

  ! setup stiffness matrix and load vector
  K(0,0) = - 0.5d0*((u(1)+u(0))/(x(1)-x(0))) 
  K(1,0) =  0.5d0*((u(1)+u(0))/(x(1)-x(0)))
  E(0) =  -(1.0d0/(m+1.0d0))*(1.0d0/(x(1)-x(0)))*(u(1)**(m+1) - u(0)**(m+1))

  do i=1, nodes-1 
     K(i-1,i)=0.5d0*((u(i)+u(i-1))/(x(i)-x(i-1)))
     K(i,i)=-0.5d0*((u(i)+u(i-1))/(x(i)-x(i-1))) - 0.5d0*((u(i+1)+u(i))/(x(i+1)-x(i)))
     K(i+1,i) = 0.5d0*((u(i+1)+u(i))/(x(i+1)-x(i)))
     E(i) = (1.0d0/(m+1.0d0))*(1.0d0/(x(i)-x(i-1)))*(u(i)**(m+1) - u(i-1)**(m+1)) -&
          & (1.0d0/(m+1.0d0))*(1.0d0/(x(i+1)-x(i)))*(u(i+1)**(m+1) - u(i)**(m+1))
  end do

  K(nodes-1,nodes) = 0.5d0*((u(nodes)+u(nodes-1))/(x(nodes)-x(nodes-1)))
  K(nodes,nodes) = -0.5d0*((u(nodes)+u(nodes-1))/(x(nodes)-x(nodes-1)))    
  E(nodes) = (1.0d0/(m+1.0d0))*(1.0d0/(x(nodes)-x(nodes-1)))*(u(nodes)**(m+1) - u(nodes-1)**(m+1))

  ! Solve stiffess system for \Psi imposing dirichlet boundary conditions strongly (\Psi=0).    
  call gaussian_elimination(K(1:nodes-1,1:nodes-1),psi(1:nodes-1),E(1:nodes-1),nodes-1)
  psi(0)=0.0d0; psi(nodes)=0.0d0

  ! Recover the mesh velocity from psi
  call velocity_recovery(x,x_dot,psi,nodes)

  ! Time-step the mesh using forward Euler.
  x(:) = x(:) + delta_t*x_dot(:)

  return

end subroutine mass_mesh_velocity

subroutine velocity_recovery(x,x_dot,psi,nodes)
  use linear_solvers

  implicit none
  !---------------------------------------------------------------------------------
  integer, intent(IN) :: nodes
  double precision, intent(INOUT), dimension(0:nodes) :: x_dot
  double precision, intent(IN), dimension(0:nodes) :: x, psi
  !---------------------------------------------------------------------------------
  double precision, dimension(0:nodes,0:nodes) :: Mass
  double precision, dimension(0:nodes) :: f
  integer :: i
  !--------------------------------------------------------------------------------
  Mass=0.0d0; f=0.0d0; x_dot=0.0d0

  ! setup mass matrix and load vector for recovery of mesh velocity
  Mass(0,0) = (1.0d0/3.0d0)*(x(1)-x(0))
  Mass(1,0) = (1.0d0/6.0d0)*(x(1)-x(0))
  f(0) = 0.5d0*(psi(1)-psi(0))

  do i=1, nodes-1
     Mass(i-1,i) = (1.0d0/6.0d0)*(x(i)-x(i-1))
     Mass(i,i) = (1.0d0/3.0d0)*(x(i+1)-x(i-1))
     Mass(i+1,i) = (1.0d0/6.0d0)*(x(i+1)-x(i))
     f(i) = 0.5d0*(psi(i+1)-psi(i-1))
  end do

  Mass(nodes-1,nodes) = (1.0d0/6.0d0)*(x(nodes)-x(nodes-1))
  Mass(nodes,nodes) = (1.0d0/3.0d0)*(x(nodes)-x(nodes-1))   
  f(nodes) = 0.5d0*(psi(nodes)-psi(nodes-1))

  ! Solve for mesh velocity
  call gaussian_elimination(Mass,x_dot,f,nodes+1)

  return

end subroutine velocity_recovery

subroutine mass_u_calc(u,x,nodes,theta)
  use linear_solvers
  !*********************************************************************************
  !**                                                                             **
  !**     This subroutine solves the matrix system:-                              **
  !**                                                                             **
  !**                     M U = B                                                 **
  !**                                                                             **
  !*********************************************************************************
  implicit none
  !---------------------------------------------------------------------------------
  integer, intent(IN) :: nodes
  double precision, intent(INOUT), dimension(0:nodes) :: u
  double precision, intent(IN), dimension(0:nodes) :: x, theta
  !---------------------------------------------------------------------------------
  double precision, dimension(0:nodes,0:nodes) :: M
  integer :: i
  !---------------------------------------------------------------------------------
  u=0.0d0; M=0.0d0

  ! Setup mass matrix system to recover solution using mass monitor     
  M(0,0) = (1.0d0/3.0d0)*(x(1)-x(0))
  M(1,0) = (1.0d0/6.0d0)*(x(1)-x(0))

  do i=1, nodes-1
     M(i-1,i) = (1.0d0/6.0d0)*(x(i)-x(i-1))
     M(i,i) = (1.0d0/3.0d0)*(x(i+1)-x(i-1))
     M(i+1,i) = (1.0d0/6.0d0)*(x(i+1)-x(i))
  end do

  M(nodes-1,nodes) = (1.0d0/6.0d0)*(x(nodes)-x(nodes-1))
  M(nodes,nodes) = (1.0d0/3.0d0)*(x(nodes)-x(nodes-1))                    

  ! Solve for solution u 
  call gaussian_elimination(M(1:nodes-1,1:nodes-1),u(1:nodes-1),theta(1:nodes-1),nodes-1)
  u(0)=0; u(nodes)=0

  return

end subroutine mass_u_calc

subroutine adaptive_timestep(u,x,m,delta_t,t,nodes) 

  implicit none
  !------------------------------------------------------------------------------        
  integer, intent(IN) :: nodes
  double precision, intent(IN), dimension(0:nodes) :: u, x
  double precision, intent(IN) :: t, m                 
  double precision, intent(INOUT) :: delta_t
  !------------------------------------------------------------------------------        
  double precision :: min_delta_x, uomax
  integer :: i                                    
  !------------------------------------------------------------------------------                                       
  min_delta_x=1000.0d0

  do i=1, nodes
     min_delta_x = dmin1((x(i)-x(i-1))**2,min_delta_x)
  end do

  uomax = MAXVAL(u)
  delta_t = (0.1d0/(uomax**m))*(t**(1.0d0/(2.0d0+m)))*min_delta_x

  return

end subroutine adaptive_timestep
