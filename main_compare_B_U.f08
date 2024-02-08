!=============================================================================
! Author: Javier Aguilar
! Date: February 8, 2024
! Description: Main program, depends on the rest of modules (declarations, functions and subroutines), 
!also requires subrutines for pseudo-random number generation)
!
! This program is distributed under the terms of the MIT license.
!
! If you use this code in your work, please cite the following article:
! Aguilar J et al. "Biased versus unbiased numerical methods for stochastic simulations: application to contagion processes."
! DOI: 
!=============================================================================

program YourProgram
  ! Your Fortran code goes here
  ...

end program YourProgram

program Main_SIS
  !
  !
  use declarations_module
  use functions
  use subroutines
  implicit none
  integer*8 :: ii,Mb
  double precision :: start,finish,tB,tG
  double precision :: alpha,xG,eG,xB,eB,max_dec,min_dec
  integer*8 :: metarealiz
  character(8)  :: date
  character(10) :: clock
  character(5)  :: zone
  integer,dimension(8) :: values
  character (len=100) :: addname


  call assingments()

  call read_neigbors(filename="networks/lattice_2D_"//trim(str(int(V)))//".dat") !Read edge list
  
  !Check compativility of computer memory with array size----------------
  k_max=maxval(E(:,0)) !maximum degree of populations
  print *, "kmax=",k_max
  print *,"Memory to store arrays=",(dble(3)*8*V+V*(kmax+1)*8)/10**6,"Mb"
  ii=N0*N0
  if (( ii<0 ).or.(ii.ne.ii)) then
    print*, "My integers are too big!!",ii
    goto 101
  end if
  !----------------------------------------------------------------------

  call dran_ini(time()+jjcontrol) !Initialize random variables, jjcontrol is a global variable

  !With addname I modify the file name and path for data outcome

  ! addname="compare_ATA_N0_"
  ! addname="compare_metapop_N0_"
  ! addname="compare_Raul_SEIR_"
  addname="compare_ZBQLBIN_SEIR_"
  ! addname="compare_ignbin_SEIR_"
  ! addname="compare_SEIR_ATA_"

  open(unit=9012, file="data/"//trim(addname)// &
  "realiz_"//trim(str(int(realiz)))// &
  "_Lattice_V"//trim(str(int(V)))// &
  "_R0"//trim(str(nint(R0*10)))// &
  "_n"//trim(str(int(jjcontrol)))//".dat", &
   status="unknown", action="write") ! writes h,b and final iteration STATISTICS OF THE WHOLE

  if ( ios /= 0 ) stop "Error opening file "

  print*, "--------------------- PROGRAM STARTS------------------"
  
  metarealiz=12 !number of datapoints
  max_dec=4.0d0 !maximum precision
  min_dec=1.0d0 !minimum precision
  call cpu_time(start)
  do ii = 1, metarealiz
    
    realiz=int8(1.0d0+10.0d0**(max_dec*(dble(ii)/dble(metarealiz))))

    print *, "realiz=",realiz,"N0=",N0
    call Compare_CPU_t_G_and_B(alpha,xG,eG,xB,eB,Mb,tB,tG)

    write(9012,*) realiz,tmax,N0,R0,alpha,xG,eG,xB,eB,Mb,Dt,tB,tG !Save data
    write(*,*)    realiz,tmax,N0,R0,alpha,xG,eG,xB,eB,Mb,Dt,tB,tG
  end do

  call cpu_time(finish)
  call date_and_time(DATE=date,TIME=clock,ZONE=zone,VALUES=values)

  write(*,*) "This took:",finish-start,"seconds"
  print *, "Clock says:", values(5),":",values(6)

  !Deallocate and close files
  close(9012)

 
  deallocate(S);deallocate(I);deallocate(DI);deallocate(DS);deallocate(E) 


101 end program Main_SIS
