subroutine read_inp()
    use input_data
    implicit none
    !
    read (*,*) N, iTst, iEQ, iPd
    read (*,*) nbin
    read (*,*) eta, res_eta, dr
    read (*,*) ibeg
    read (*,*) tRes
    !
    return
end subroutine read_inp

subroutine init()
    use const
    use input_data
    use var
    implicit none
    !
    real(8) init_space
    integer  :: i,j,n_fcc,i_atom, ii,jj,kk
    !
    ! greeting
    write(*,*) 
    write(*,*) '================================================='
    write(*,*) '== Hard Sphere RDF Simulation with Monte Carlo =='
    write(*,*) '== ------------------------------------------- =='
    write(*,*) '== Final Exam Project for StatMech, 2013 Fall  =='
    write(*,*) '==                                             =='
    write(*,*) '==   Author:        Hsin-Yu Ko                 =='
    write(*,*) '==   Email:         hsinyu@princeton.edu       =='
    write(*,*) '==   Department of Chemistry, Princeton Univ.  =='
    write(*,*) '==   Ver. 1.0   Wed Nov 20 17:04:03 EST 2013   ==' 
    write(*,*) '================================================='
    write(*,*) 
    call RANDOM_SEED()
    ! allocate and initialize the variables
    allocate(r(3,N));   r = 0.d0
    allocate(r_t(3,N)); r_t = 0.d0
    allocate(gr(nbin)); gr=0
    ! init
    tOverlap=.true.
    tGetRes =.false.
    iter  = 0
    ! setting the L = 1 => rho = N
    ! calculate eta=rho*pi*D**3/6
    D     = (6.d0*eta/(pi*N))**(1.d0/3.d0)
    D_res = (6.d0*res_eta/(pi*N))**(1.d0/3.d0)
    ! initial position
    if (ibeg.eq.'cub') then
      i = 1+int(N**(1.d0/3.d0))
      init_space = 1.d0/dble(i)
      if (init_space.le.D) then
        write(*,*) 'not able to init with grid point initialization...STOP!!'
        stop
      end if
      do j=1,N
        r(1,j) = mod(j,i)
        r(2,j) = mod(int(j/i),i)
        r(3,j) = mod(int(j/i/i),i)
      end do
      r(:,:) = init_space*r(:,:) 
    else if (ibeg.eq.'res') then
      open(unit=20,file='out/restart.dat',action='read')
      do i=1,N
        read(20,*) r(1,i),r(2,i),r(3,i)
      end do
      close(20)
    else if (ibeg.eq.'fcc') then
      init_space = 1.1d0*dsqrt(2.d0)*D ! make some space
      n_fcc = INT((N/4.D0)**(1.D0/3.D0))+1
      i_atom = 0
      DO ii = 0, n_fcc - 1
        DO jj = 0, n_fcc - 1
          DO kk = 0, n_fcc - 1
            i_atom = i_atom + 1
            r(1,i_atom) = kk
            r(2,i_atom) = jj
            r(3,i_atom) = ii
            if (i_atom.eq.N) goto 123
            !
            i_atom = i_atom + 1
            r(1,i_atom) = (kk+0.5D0)
            r(2,i_atom) = (jj+0.5D0)
            r(3,i_atom) =  ii        
            if (i_atom.eq.N) goto 123
            !
            i_atom = i_atom + 1
            r(1,i_atom) = kk
            r(2,i_atom) = (jj+0.5D0)
            r(3,i_atom) = (ii+0.5D0)
            if (i_atom.eq.N) goto 123
            !
            i_atom = i_atom + 1
            r(1,i_atom) = (kk+0.5D0)
            r(2,i_atom) = jj    
            r(3,i_atom) = (ii+0.5D0)
            if (i_atom.eq.N) goto 123
          END DO
        END DO
      END DO
      123   continue
      r = r*init_space
      !     debug
      !       do i_atom = 1,N
      !         write(10,*) (r(j,i_atom),j=1,3)
      !       end do
    else
      write (*,*) "don't know how to do with, ibeg= ",ibeg
      stop
    end if
    !
    r_t = r
    call check_overlap()
    if (tOverlap) write (*,*) 'wrong initial packing'
    !
    open(unit=30,file='out/eta',action='write')
    write(30,*) eta
    close(30)
    ! 
    return
end subroutine init

subroutine finalize()
    use var
    implicit none
    !
    if (allocated(r))   deallocate(r)
    if (allocated(r_t)) deallocate(r_t)
    if (allocated(gr))  deallocate(gr)
    write(*,*) 
    write(*,*) '================================================='
    write(*,*) '==                   JOB DONE                  ==' 
    write(*,*) '================================================='
    write(*,*) 
    !
    return
end subroutine finalize

subroutine check_overlap()
    use input_data
    use var
    implicit none
    !
    integer :: i, j
    real(8) :: min_img_dist, dis,min_spa
    !
    min_spa = 1.d0
    tOverlap = .false.
    do i=1,N-1
      do j=i+1,N
        dis = min_img_dist(r_t(:,i),r_t(:,j))
        if (dis.lt.D) then
          tOverlap = .true.
          exit
        else
          tOverlap = .false.
          min_spa = min(min_spa,dis)
        end if
      end do
      if (tOverlap) exit
    end do
    !
    if (.not.tOverlap) then
      if (iter.gt.0.99*iEQ) then
        if ((min_spa.gt.D_res).and.(.not.tGetRes).and.(tRes)) then
          ! write restart file
          open(unit=22,file='out/restart.dat',action="write")
          do i=1,N
            write(22,"(3F14.8)") r_t(1,i),r_t(2,i),r_t(3,i)
          end do
          close(22)
          tGetRes=.true.
        end if
      end if
    end if
    ! number of iter grows anyway
    iter = iter + 1
    !
    return
end subroutine check_overlap

function min_img_dist(v1,v2)
    ! this is in reduced box L = 1
    implicit none
    real(8) :: min_img_dist, v1(3), v2(3), dv(3)
    integer :: i
    min_img_dist = 0.d0
    do i=1,3
      dv(i)=v1(i)-v2(i)
      dv(i)=dv(i)-nint(dv(i))
      min_img_dist = min_img_dist + dv(i)*dv(i)
    end do
    min_img_dist = dsqrt(min_img_dist)
    return
end function min_img_dist

subroutine mcmove()
    use input_data
    use var
    implicit none
    !
    integer   :: o, i
    real(8)   :: ranf(4) 
    ! store trial
    r_t = r
    ! select a particle at random
    call RANDOM_NUMBER(ranf(:))
    o = INT(N*ranf(1)) + 1
    !give particle a random displacement
    do i = 1,3
      r_t(i,o) = r(i,o) + (ranf(1+i)-0.5D0)*Dr
    end do
    !acceptance test
    call check_overlap()
    !
    if (.not.tOverlap) then
      ! accepted
      Nacc = Nacc + 1
      do i = 1,3
        IF (r_t(i,o).LT.0.d0) r_t(i,o) = r_t(i,o) + 1.d0
        IF (r_t(i,o).GT.1.d0) r_t(i,o) = r_t(i,o) - 1.d0
        r(i,o) = r_t(i,o)
      end do
    end if
    !
    return
end subroutine mcmove

subroutine rate_accept()
    use input_data
    use var
    implicit none
    !
    integer :: i 
    !
    do i = 1, iTst
      call mcmove()
    end do
    !
    write (*,*) 'the acceptance ratio is ',dble(Nacc)/dble(iter)
    ! auto adjust the move step (optional)
    ! call adjust_dr()
    ! manual adjust
    if (dble(Nacc)/dble(iter).lt.0.04) then
      write (*,*) 'too small acceptance ratio => smaller Dr ...STOP!!' 
    else if (dble(Nacc)/dble(iter).gt.0.9) then
      write (*,*) 'too large acceptance ratio => larger Dr ...STOP!!' 
    else
      write (*,*) 'Dr looks good'
    end if
    !
    return
end subroutine rate_accept

subroutine g_2(Iconf)
    use var
    use input_data
    use const
    implicit none
    ! declare variables
    character(2)             :: atom1, atom2
    real(8)                  :: cel(3,3), time, r_max, delta_r,norm,&
        min_img_dist, delta_V,f_t_pi_r_cub, omega, omg_nat_conf,det,&
        start_percent
    real(8),allocatable      :: g_pair(:,:)
    character(256)           :: fname, cConf
    !
    integer :: i,j,Iconf,ig,at1,at2
    !
    !!!initialization!!!
    allocate(g_pair(2,nbin)); g_pair=0.D0
    !!!initialization!!! END
    r_max = 0.5d0
    !! init g_pair
    delta_r = r_max / nbin
    !! init g_pair END
    omega = 1.d0
    do at1 = 1, N
      do at2 = 1, N
        if (at1.ne.at2) then
          ig = int(min_img_dist(r(:,at1),r(:,at2))/delta_r)+1
          if (ig.le.nbin) then
            gr(ig) = gr(ig) + 1
          end if
        end if
      end do
    end do
    ! normalization
    if (Iconf.eq.iPd) then
      ! normalization factor
      f_t_pi_r_cub = pi*4.D0/3.D0*delta_r*delta_r*delta_r
      do i = 1,nbin
        g_pair(1,i) = (i-0.5D0)*delta_r/D ! use D as unit
      end do
      omg_nat_conf = omega/(N*N*iPd)
      do i = 1,nbin
        j = i-1
        delta_V = f_t_pi_r_cub*(i*i*i-j*j*j)
        g_pair(2,i) = g_pair(2,i)+gr(i)/delta_V*omg_nat_conf
      end do
      ! normalization END
      !!! start g_pair calculation !!! END
      !!! output !!!
      open(unit=24,file='out/g_2.dat',action='write')
      do i = 1,nbin
        write(24,'(2F14.8)') g_pair(1,i),g_pair(2,i)
      end do
      close(24)
      !!! output !!! END
    end if
    !
    return
    !
end subroutine g_2

