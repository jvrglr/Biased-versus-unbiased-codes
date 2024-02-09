
!=============================================================================
! Author: Javier Aguilar
! Date: February 8, 2024
! Description:   Module for public variables to be used in the rest of the modules
!
! This program is distributed under the terms of the MIT license.
!
! If you use this code in your work, please cite the following article:
! Aguilar J et al. "Biased versus unbiased numerical methods for stochastic simulations: application to contagion processes."
! DOI:
!=============================================================================
module declarations_module
  implicit none
  public

  integer*8 :: ios,V,realiz,k_max,lambda,N0,kmax,jjcontrol,seeds,degree
  double precision :: t,mu,R0,beta,tmax,h,M,Dt,pM,pmu,e0,gamma,pgamma
  integer*8, dimension (:),allocatable :: S,I,DI,DS,Ex,DEx,R,DR,N
  integer*4, dimension (:,:),allocatable :: E


contains

  subroutine assingments()
    implicit none
    jjcontrol=1 !label outcome file, also used in seed for random number generator
    realiz=1000 !M in the paper, number of realizations, redefined in subroutines or main
    V=100 !Number of populations
    N0=10**3 !Number of initial inhabitants in every population
    degree=4 !Number of links of every population (in a regular lattice)
    tmax=7.5d0 !2.4d0 !7.5d0  !20.0d0 !maximum simulation time
    t=0.0d0 !current simulation time
    allocate(S(V));allocate(I(V));allocate(DS(V));allocate(DI(V))
    allocate(Ex(V));allocate(DEx(V));allocate(R(V));allocate(DR(V))
    seeds=int8(0.01*N0) !initial number of infected individuals

    R0=4.0d0 !Basic reproductive number
    mu=1.0d0 !Recovery rate
    gamma=1.0 !Rate Exposed->Infected
    M=1.0d0 ! Mobility rate
    beta=R0*mu
    lambda=V*N0*M !total mobility rate

    e0=0.045 !Scaling of errors binomial method
    print *, "parameters read and defined properly"



  end subroutine assingments

end module declarations_module
