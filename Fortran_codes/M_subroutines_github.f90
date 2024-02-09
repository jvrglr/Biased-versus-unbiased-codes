!=============================================================================
! Author: Javier Aguilar
! Date: February 8, 2024
! Description:  Module with public subroutines. Depends on declarations_module and functions.
!
! This program is distributed under the terms of the MIT license.
!
! If you use this code in your work, please cite the following article:
! Aguilar J et al. "Biased versus unbiased numerical methods for stochastic simulations: application to contagion processes."
! DOI: 
!=============================================================================
module Subroutines


  use declarations_module
  use functions

contains

  subroutine Compare_CPU_t_G_and_B(a,xG,eG,xB,eB,Mb,tB,tG)
    !Compute density of infected individuals with fixed precision
    !using both Gillespie and binomial methods
    !Input: eG, scaling of errors of the biased method (provided in declarations module)
    implicit none
    double precision, intent(inout):: a,xG,eG,xB,eB,tG,tB
    integer*8, intent(inout) :: Mb
    integer*8 :: ii
    double precision :: start,finish,dumt,dumI
    double precision :: xm,xm2

    tG=0.0d0
    tB=0.0d0

    !--------------Measure Gillespie---------
    xm=0.0d0;xm2=0.0d0
    call cpu_time(start)
    do ii = 1, realiz, 1
      !Initial condition---------------------
      S=N0;S(1)=N0-seeds 
      I=0;I(1)=seeds
      Ex=0;R=0
      t=0.0d0
      !--------------------------------------
      call Gillespie_SEIR(dumt,dumI)  !Unbiased Metapopulation SEIR 
      xm=xm+dumI;xm2=xm2+dumI*dumI  !measure
    enddo
    call cpu_time(finish)
    tG=(finish-start) !T CPU with unbiased method
    !----------------Prepare binomial--------
    xm=xm/realiz;xm2=xm2/realiz
    xG=xm
    eG=sqrt(abs(xm2-xm**2)/realiz)  !unbiased precision
    Dt=eG/e0/3.0d0 !optimal discretization time
    Mb=int8(2.25d0*realiz) !optimal number of realizations
    !----------------Measure binomial--------
    xm=0.0d0;xm2=0.0d0
    call cpu_time(start)
    do ii = 1, Mb, 1
      !Initial condition---------------------
      S=N0;S(1)=N0-seeds
      I=0;I(1)=seeds
      Ex=0;R=0
      t=0.0d0
      !--------------------------------------
      call Binomial_SEIR(dumt,dumI) !Biased Metapopulation SEIR
      xm=xm+dumI;xm2=xm2+dumI*dumI
    enddo
    call cpu_time(finish)
    xm=xm/Mb;xm2=xm2/Mb
    eB=sqrt(abs(xm2-xm**2)/Mb)
    tB=(finish-start) !T CPU with biased method
    xB=xm

    a=tB/tG !Numerical value for alpha
    print*, "alpha=",a,"tG=",tG,"tB=",tB,"eG=",eG,"Mb=",Mb,"dt=",Dt
    print*, "xt (Gillespie)=",xG,"pm",eG
    print*, "xt (Binomial )=",xm,"pm",eB !plus e0+Dt

  end subroutine Compare_CPU_t_G_and_B

  subroutine Binomial_SEIR(tf,xf)
    !Compute SEIR on meta-population topology using Binomial method
    !Outcomes: final time-->tf
    !          final density of infected individuals
    !IMPORTANT: THE USER CAN USE THREE DIFFERENT BINOMIAL GENERATOR, SEE APPENDIX D.
    implicit none
    integer*8 :: sumI,I_travel,S_travel,Ex_travel,R_travel,ii,k,j,dum,l
    integer*8 :: nEx,nI,nR
    double precision :: p,bdt
    double precision, intent(out) :: tf,xf
    integer*8,dimension(V) :: DI,DS
    integer :: ZBQLBIN !Binomial (ii) randgen
    integer :: ignbin !Binomial (i) kachitvichyanuku
    integer :: iran_bin !Binomial (ii) dranxor

    sumI=sum(I)

    !-------Jump probabilities-------
    pmu=1.0d0-exp(-mu*Dt)
    pgamma=1.0d0-exp(-gamma*Dt)
    pM=1.0d0-exp(-M*Dt)
    bdt=beta*Dt
    !--------------------------------
    
    do while (( sumI.gt.0 ).and.(t.lt.tmax))
      t=t+Dt
      DS=0;DEx=0;DI=0;DR=0 !Initialize increments

      do ii = 1, V, 1 !Move people
        ! Extract people: Binomial sampling ----------------------
        S_travel =ZBQLBIN(S(ii), pM)
        Ex_travel=ZBQLBIN(Ex(ii),pM)
        I_travel =ZBQLBIN(I(ii), pM)
        R_travel =ZBQLBIN(R(ii), pM)
        DI(ii)=DI(ii)-I_travel
        DS(ii)=DS(ii)-S_travel
        DEx(ii)=DEx(ii)-Ex_travel
        DR(ii)=DR(ii)-R_travel
        k=E(ii,0) !Degree of ii
        !---------------------------------------------------------
        do j = 1, k, 1
          !Distribute people extracted: multinomial sampling------
          l=E(ii,j) !Label of neighbor
          p=1.0d0/(k-j+1)
          dum=ZBQLBIN(I_travel,p); DI(l)=DI(l)+dum;  I_travel=I_travel-dum
          dum=ZBQLBIN(S_travel,p); DS(l)=DS(l)+dum;  S_travel=S_travel-dum
          dum=ZBQLBIN(Ex_travel,p);DEx(l)=DEx(l)+dum;Ex_travel=Ex_travel-dum
          dum=ZBQLBIN(R_travel,p); DR(l)=DR(l)+dum;  R_travel=R_travel-dum
          !-------------------------------------------------------
        end do
      end do

      I=I+DI;S=S+DS;Ex=Ex+DEx;R=R+DR

      do ii = 1, V, 1 !Epidemics
        nEx=ZBQLBIN(S(ii),1.0d0-exp(-bdt*I(ii)/dble(I(ii)+S(ii)+Ex(ii)+R(ii)))) !S->Ex
        nI=ZBQLBIN(Ex(ii),pgamma) !Ex->I
        nR=ZBQLBIN(I(ii),pmu) !I->R
        S(ii)=S(ii)-Nex
        Ex(ii)=Ex(ii)+Nex-nI
        I(ii)=I(ii)+nI-nR
        R(ii)=R(ii)+nR
      end do
      sumI=sum(I)
    end do !End do over time

    tf=t
    xf=dble(sumI)/dble(V*N0)
101  end subroutine Binomial_SEIR

  subroutine Gillespie_SEIR(tf,I_tf)
    !Compute SEIR on meta-population topology using Gillespie method
    !Outcomes: final time-->tf
    !          final density of infected individuals-->xf
    !IMPORTANT: if random generator for uniform real random variables can have 0 as outcome, program will fail (infinite jumping times)
    implicit none
    double precision, intent(out) :: tf,I_tf
    double precision :: W,pb,pmm,u,pg,pmob,pI,pR,pS
    double precision, dimension(V) :: q
    integer*8 :: sumI,sumEx,pos,count,k,count_mov
    integer*8 :: or,de
    double precision :: sumSIN
    double precision :: dran_u
    integer*8 :: i_dran

    N=S+Ex+I+R
    sumI=sum(I);sumSIN=sum(I*(S/dble(N)));sumEx=sum(Ex)
    count=0;count_mov=0

    do while (( sumI.gt.0 ).and.(t.lt.tmax))

      pb=beta*sumSIN;pmm=mu*sumI;pg=gamma*sumEx !Transition rates
      W=pmm+pb+lambda+pg !Total exit rate
      u=dran_u() 
      t=t-log(u)/W !Update time
      pmob=lambda/W
      u=dran_u()
      if ( u.lt.pmob ) then
        !Move agent
        count_mov=count_mov+1
        q=N/dble(N0*V)
        u=dran_u()
        call search_list_binary_algoritm(q,or,u) !Choose origin
        k=E(or,0)
        de=E(or,i_dran(k)) !choose destination
        u=dran_u()
        pI=dble(I(or))/N(or);pS=pI+dble(S(or))/N(or);pR=pS+dble(R(or))/N(or)
        if ( u.lt.pI ) then !Chose compartment
          I(or)=I(or)-1
          I(de)=I(de)+1
        else if ( u.lt.pS ) then 
          S(or)=S(or)-1
          S(de)=S(de)+1
        else if ( u.lt.pR ) then 
          R(or)=R(or)-1
          R(de)=R(de)+1
        else
          Ex(or)=Ex(or)-1
          Ex(de)=Ex(de)+1
        end if
        N(or)=N(or)-1
        N(de)=N(de)+1
      else
        !If there is not movement, then update epidemics
        pb=pb/W;pg=pg/W
        if ( u.lt.(pb+pmob) ) then
          !Create E due to contact with infected
          q=S*I/dble(N)/dble(sumSIN)
          u=dran_u()
          call search_list_binary_algoritm(q,pos,u)
          Ex(pos)=Ex(pos)+1;S(pos)=S(pos)-1
        else if (u.lt.(pmob+pb+pg)) then
          q=Ex/dble(sumEx)
          u=dran_u()
          call search_list_binary_algoritm(q,pos,u)
          Ex(pos)=Ex(pos)-1;I(pos)=I(pos)+1 !R=R+1
        else
          !Destroy one infected
          q=I/dble(sumI)
          u=dran_u()
          call search_list_binary_algoritm(q,pos,u)
          I(pos)=I(pos)-1;R(pos)=R(pos)+1 !R=R+1
        end if

        sumI=sum(I);sumSIN=sum(I*(S/dble(N)));sumEx=sum(Ex)

        if ( t.gt.count ) then
          count=count+1
        end if

      end if !This is the end if that discerns between update due to epidemics or mobiliy

    end do !Main do over time and I=0


    tf=t
    I_tf=dble(sumI)/dble(V*N0)
  end subroutine Gillespie_SEIR

  subroutine search_list_binary_algoritm(list,position,p)
    !Rafle event:
    !Given a list of probabilities called "list" such that sum(list)=1 and a probability p.
    !Look for "position" such that C(position)>=p and C(j)<p for all j in [1,position[.
    !Where C is the cumulative of list: C(i)=list(1)+list(2)+...+list(i)-
    !REFERENCE:Brainerd, W. S. (2015). Guide to Fortran 2008 programming (p. 141). Berlin: Springer.

    implicit none
    double precision, dimension(:), intent (in) :: list
    double precision, intent (in) :: p
    integer*8, intent(out) :: position
    double precision, dimension(size(list)) :: C
    integer*8 ii,N,first,last,half


    N=size(list) !It would be cool to define this as a parameter (constant), I don't know how...
    C(1)=list(1)
    do ii = 2, N, 1 !Compute cumulative of list
      c(ii)=C(ii-1)+list(ii)
    end do

    first=1;last=N
    do while ( first.ne.last )
      half=(first+last)/2
      if ( p>C(half) ) then
        first=half+1
      else
        last=half
      end if
    end do

    position=first

  end subroutine search_list_binary_algoritm

  subroutine read_neigbors(filename)
    !Read edge list from file
    implicit none
    character (len=*),intent(in) :: filename
    Allocate(E(V,0:degree))
    open(unit=2002, file=trim(filename), iostat=ios, status="old", action="read")
    if ( ios /= 0 ) stop "Error opening file "

    read(2002,*) E

    close(2992)
    return

  end subroutine read_neigbors

  subroutine read_xF_tf()
    !Example of subroutine, read data file
    implicit none
    integer*8 :: ios,i,dum,datapoints
    double precision :: dummy

    open(unit=3003, file="data/F_Target_test3.dat", iostat=ios, status="old", action="read")
    if ( ios /= 0 ) stop "Error opening file "
    do i = 1, datapoints, 1
      read(3003,*) dummy,dummy,dummy,dum
    end do


  end subroutine read_xF_tf

end module subroutines
