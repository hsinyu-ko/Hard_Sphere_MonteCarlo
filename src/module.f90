module const
  implicit none
  real(8), parameter  :: pi = dacos(-1.d0)
end module const

module input_data
  implicit none
  integer(8)   ::  N        ! number of particle in the unit cell
  integer(8)   ::  iTst     ! number of particle in the unit cell
  integer(8)   ::  iEQ      ! number of particle in the unit cell
  integer(8)   ::  iPd      ! number of particle in the unit cell
  integer(8)   ::  nbin     ! number of bin in g_2
  real(8)      ::  eta      ! the reduced density
  real(8)      ::  res_eta  ! restart file writing threathold
  real(8)      ::  Dr       ! the displacement of MC move
  character(3) ::  ibeg     ! the reduced density
  logical      ::  tRes     ! want to print restart file

end module input_data

module var
  implicit none
  real(8), allocatable    :: r(:,:)   ! r(xyz,Iparticle)
  real(8), allocatable    :: r_t(:,:) ! trial r
  integer(8), allocatable :: gr(:)    ! temp g_2 
  real(8)                 :: D        ! the hard sphere radius
  real(8)                 :: D_res    ! the restart D
  real(8)                 :: x_avg    ! the averaged x used to check equilibration
  real(8)                 :: y_avg    ! the averaged y used to check equilibration
  real(8)                 :: z_avg    ! the averaged z used to check equilibration
  logical                 :: tOverlap ! check overlap
  logical                 :: tGetRes  ! if get restart file
  integer(8)              :: iter     ! number of iterations
  integer(8)              :: Nacc     ! number of accepted moves
  ! un init_

end module var
