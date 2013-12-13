program hard_sp_nvt
  use const
  use input_data
  use var
  implicit none
  integer  :: i, j, ii
  !
  ! get input data
  call read_inp()
  ! start with random orientation for initial structure
  call init()
  ! check approximate acceptance rate
  call rate_accept()
  ! equilibration
  do i=1,iEQ
    x_avg = 0.d0; y_avg = 0.d0; z_avg = 0.d0
    call mcmove()
    ! check averaged x for equilibration
    do ii = 1,N
      x_avg = x_avg + r(1,ii)
      y_avg = y_avg + r(2,ii)
      z_avg = z_avg + r(3,ii)
    end do
    ! 
    open(unit=25,file='out/r_av',action='write',access='append')
    write (25,"(I6,3F14.8)") i, x_avg/N, y_avg/N, z_avg/N
    close(25)
    !
    if (mod(i,int(0.1*iEQ)).eq.0) write(*,'("Equil: Step = ",I8," / ",I8)') i,iEQ
  end do
  ! production
  do j=1,iPd
    call mcmove()
    ! harvest pair correlation function
    call g_2(j)
    !
    if (mod(j,int(0.1*iPd)).eq.0) write(*,'("Prod : Step = ",I8," / ",I8)') j,iPd
  end do
  ! 
  if (tRes) then
    if (tGetRes) then
      write (*,*) 'Successfully get targeting restart file'
    else
      write (*,*) 'Did not get targeting restart file'
    end if
  end if
  !
  ! finalize
  call finalize()
  !
  stop
end program hard_sp_nvt
