!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
!                           Heat exchange sub-model of ELCIRC				!
!                       	Version 5 (Sept. 05, 2003)                              !
!                                                                                       !
!                 Center for Coastal and Land-Margin Research                           !
!             Department of Environmental Science and Engineering                       !
!                   OGI School of Science and Engineering,                              !
!                     Oregon Health & Science University                                !
!                       Beaverton, Oregon 97006, USA                                    !
!                                                                                       !
!                   Scientific direction: Antonio Baptista                              !
!                   Code development: Mike A. Zulauf                  			!
!                                                                                       !
!               Copyright 1999-2003 Oregon Health and Science University                !
!                              All Rights Reserved                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  

!-----------------------------------------------------------------------
! Note: the following global variables (from module global) are used in
!       this code.  This list does not include variables passed in as
!       arguments. . .
!
!       mnp
!       np
!       x
!       y
!       kfp
!       uu2
!       vv2
!       tnd
!       snd
!-----------------------------------------------------------------------
      subroutine get_wind (time, u_air_node, v_air_node, p_air_node, &
     &                     t_air_node, q_air_node)

!       implicit none
        use global
        implicit real*8(a-h,o-z),integer(i-n)

! define some new names for things in header file
        integer max_nodes
        parameter (max_nodes = mnp)

! input/output variables
        real*8 u_air_node(max_nodes), v_air_node(max_nodes)
        real*8 t_air_node(max_nodes), q_air_node(max_nodes)
        real*8 p_air_node(max_nodes), time

! local variables
        integer max_ni, max_nj, max_times, max_files
        integer max_elems_in, max_nodes_in
        parameter (max_ni = 1024)
        parameter (max_nj = 1024)
        parameter (max_elems_in = (max_ni-1) * (max_nj-1) * 2)
        parameter (max_nodes_in = max_ni * max_nj)
        parameter (max_times = 1000)
        parameter (max_files = 999)
        integer num_nodes
        integer in_elem_to_out_node_1(max_nodes)
        integer in_elem_to_out_node_2(max_nodes)
        integer ni_1, nj_1, num_times_1, num_nodes_in_1, num_elems_in_1
        integer node_i_1(max_nodes_in), node_j_1(max_nodes_in)
        integer node_num_in_1(max_nodes_in)
        integer elem_nodes_in_1(max_elems_in,3)
        integer ni_2, nj_2, num_times_2, num_nodes_in_2, num_elems_in_2
        integer node_i_2(max_nodes_in), node_j_2(max_nodes_in)
        integer node_num_in_2(max_nodes_in)
        integer elem_nodes_in_2(max_elems_in,3)
        integer num_files_1, num_files_2
        integer max_rank, rank
        parameter (max_rank = 3)
        integer dim_sizes(max_rank)
        real*8 weight_wind_node_1(max_nodes,3)
        real*8 weight_wind_node_2(max_nodes,3)
        real*8 x_in_1(max_ni,max_nj), y_in_1(max_ni,max_nj)
        real*8 x_in_2(max_ni,max_nj), y_in_2(max_ni,max_nj)
        real*8 start_day_1
        real*8 relative_weight_1, relative_weight_2
        real*8 max_window_1, max_window_2
        parameter (relative_weight_1 = 1.0)
        parameter (relative_weight_2 = 2.0)
        parameter (max_window_1 = 24.0)
        parameter (max_window_2 = 2.0)
        real*8 frac_day, secs_per_day, utc_start
        parameter (secs_per_day = 86400.0)
        parameter (utc_start = 8.0)
        real*8 temp_arr_1(max_ni,max_nj)
        real*4 temp_arr_2(max_ni,max_nj), temp_arr_3(max_ni,max_nj)
        real*4 temp_arr_4(max_ni,max_nj), temp_arr_5(max_ni,max_nj)
        real*8 temp_arr_6(max_elems_in)
        real*8 temp_arr_8(max_nodes), temp_arr_9(max_nodes)
        real*4 temp_sca
        real*4 wind_times_1(max_times), wind_times_2(max_times)
        character wind_set_1*50, wind_set_2*50, start_day_file*50
        character wind_time_files_1(max_times)*50
        character wind_time_files_2(max_times)*50
        parameter (wind_set_1 = 'hdf/wind_file_1')
        parameter (wind_set_2 = 'hdf/wind_file_2')
        parameter (start_day_file = 'hdf/start_day.txt')
        logical first_call, have_wind_2, have_start_day_file
        data first_call/.true./

! retain the values of some local variables between calls
        save first_call, start_day_1, &
     &    in_elem_to_out_node_1, weight_wind_node_1, &
     &    in_elem_to_out_node_2, weight_wind_node_2, &
     &    num_nodes, have_wind_2, &
     &    wind_times_1, num_times_1, num_files_1, ni_1, nj_1, &
     &    num_nodes_in_1, num_elems_in_1, node_i_1, node_j_1, &
     &    node_num_in_1, elem_nodes_in_1, &
     &    wind_times_2, num_times_2, num_files_2, ni_2, nj_2, &
     &    num_nodes_in_2, num_elems_in_2, node_i_2, node_j_2, &
     &    node_num_in_2, elem_nodes_in_2, &
     &    wind_time_files_1, wind_time_files_2

        open(39,file='fort.39')
        rewind(39)
        write(39,*)
        write(39,*) 'enter get_wind'
        write(39,*) 'first_call = ', first_call

! if this is the first call to this routine then get some things ready
        if (first_call) then

! define the local variables num_nodes
          num_nodes = np

! check to see if start_day_file exists
          call file_exst (start_day_file, have_start_day_file, .false.)

! if the start day file does exist, get start_day from it, otherwise
! use the first start_day in wind_set_1
          if (have_start_day_file) then
            open (unit=77, file=start_day_file, status='old')
            read(77,*) temp_sca
            close (unit=77)
            start_day_1 = temp_sca
          else
            wind_time_files_1(1) = 'hdf/wind_file_1.001.hdf'
            call read_scalar(wind_time_files_1(1), temp_sca, &
     &                       'start_day           ', 0.0)
            start_day_1 = temp_sca
          endif

! check to see if _any_ wind_file_2 exists (use first possible name)
          wind_time_files_2(1) = 'hdf/wind_file_2.001.hdf'
          call file_exst (wind_time_files_2(1), have_wind_2, .false.)
          if (.not. have_wind_2) then
            write(39,*)
            write(39,*) wind_time_files_2(1), ' not exist. . .'
          endif

! get the times of the data available in wind_set_1
          call get_times(wind_times_1, wind_set_1, &
     &                   'u                   ', &
     &                   wind_time_files_1, num_times_1, &
     &                   num_files_1, max_times, max_files)

! get the dimensions of the datasets in wind_set_1 (use first dataset)
          call get_dims(wind_time_files_1(1), 'u                   ', &
     &                  wind_times_1(1), rank, dim_sizes)
          ni_1 = dim_sizes(1)
          nj_1 = dim_sizes(2)

! check the dimensions of wind_set_1, to ensure they don't exceed the
! maximums
          if (ni_1 .gt. max_ni .or. nj_1 .gt. max_nj) then
            write(*,*)
            write(*,*) 'wind_file_1: max dimensions exceeded!'
            write(11,*)
            write(11,*) 'wind_file_1: max dimensions exceeded!'
            stop
          endif

! calculate the total number of nodes and elements for wind_set_1
          num_nodes_in_1 = ni_1 * nj_1
          num_elems_in_1 = (ni_1-1) * (nj_1-1) * 2

! check the elems/nodes of wind_set_1, to ensure they don't exceed the
! maximums
          if (num_elems_in_1 .gt. max_elems_in .or. &
     &        num_nodes_in_1 .gt. max_nodes_in) then
            write(*,*)
            write(*,*) 'wind_file_1: max elems/nodes exceeded!'
            write(11,*)
            write(11,*) 'wind_file_1: max elems/nodes exceeded!'
            stop
          endif

! create list of all nodes for wind_set_1
          call list_nodes (node_i_1, node_j_1, node_num_in_1, &
     &                     num_nodes_in_1, ni_1, nj_1)

! now create the list of all the elements (and the nodes defining them)
! for wind_set_1
          call list_elems (elem_nodes_in_1, node_num_in_1, &
     &                     ni_1, nj_1, num_elems_in_1)

! do the same as above for wind_set_2 (if it exists)
          if (have_wind_2) then

! get the times of the data available in wind_set_2
            call get_times(wind_times_2, wind_set_2, &
     &                     'u                   ', &
     &                     wind_time_files_2, num_times_2, &
     &                     num_files_2, max_times, max_files)

! get the dimensions of the datasets in wind_set_2 (use first dataset)
            call get_dims(wind_time_files_2(1), 'u                   ', &
     &                    wind_times_2(1), rank, dim_sizes)
            ni_2 = dim_sizes(1)
            nj_2 = dim_sizes(2)

! check the dimensions of wind_set_2, to ensure they don't exceed the
! maximums
            if (ni_2 .gt. max_ni .or. nj_2 .gt. max_nj) then
              write(*,*)
              write(*,*) 'wind_file_2: max dimensions exceeded!'
              write(11,*)
              write(11,*) 'wind_file_2: max dimensions exceeded!'
              stop
            endif

! calculate the total number of nodes and elements for wind_set_2
            num_nodes_in_2 = ni_2 * nj_2
            num_elems_in_2 = (ni_2-1) * (nj_2-1) * 2

! check the elems/nodes of wind_set_2, to ensure they don't exceed the
! maximums
            if (num_elems_in_2 .gt. max_elems_in .or. &
     &          num_nodes_in_2 .gt. max_nodes_in) then
              write(*,*)
              write(*,*) 'wind_file_2: max elems/nodes exceeded!'
              write(11,*)
              write(11,*) 'wind_file_2: max elems/nodes exceeded!'
              stop
            endif

! create list of all nodes for wind_set_2
            call list_nodes (node_i_2, node_j_2, node_num_in_2, &
     &                       num_nodes_in_2, ni_2, nj_2)

! now create the list of all the elements (and the nodes defining them)
! for wind_set_2
            call list_elems (elem_nodes_in_2, node_num_in_2, &
     &                       ni_2, nj_2, num_elems_in_2)

          endif ! end of have_wind_2 block

! read in the x and y values for wind_set_1, and copy to full size
! real*8 arrays
          call read_2d_arr(wind_time_files_1(1), temp_arr_2, &
     &                     'x                   ', 0.0, &
     &                     ni_1, nj_1)
          call read_2d_arr(wind_time_files_1(1), temp_arr_3, &
     &                     'y                   ', 0.0, &
     &                     ni_1, nj_1)
          call copy_arr(temp_arr_2, ni_1, nj_1, x_in_1, &
     &                    max_ni, max_nj)
          call copy_arr(temp_arr_3, ni_1, nj_1, y_in_1, &
     &                    max_ni, max_nj)

! calculate the weightings from wind_set_1 to elcirc nodes
! (this is slow)
          write(*,*)
          write(*,*) &
     &      'begin calculating grid weightings for wind_file_1'
          write(16,*)
          write(16,*) &
     &      'begin calculating grid weightings for wind_file_1'
          call get_weight (x_in_1, y_in_1, x, y, &
     &                     elem_nodes_in_1, node_i_1, node_j_1, &
     &                     max_ni, max_nj, &
     &                     num_elems_in_1, num_nodes_in_1, &
     &                     num_nodes, &
     &                     max_nodes, &
     &                     in_elem_to_out_node_1, &
     &                     temp_arr_6, weight_wind_node_1)
          write(*,*) &
     &      'done calculating grid weightings for wind_file_1'
          write(16,*) &
     &      'done calculating grid weightings for wind_file_1'

! do the same but for wind_set_2 (if it exists)
          if (have_wind_2) then

! read in the x and y values for wind_set_2, and copy to full size
! real*8 arrays
            call read_2d_arr(wind_time_files_2(1), temp_arr_2, &
     &                       'x                   ', 0.0, &
     &                       ni_2, nj_2)
            call read_2d_arr(wind_time_files_2(1), temp_arr_3, &
     &                       'y                   ', 0.0, &
     &                       ni_2, nj_2)
            call copy_arr(temp_arr_2, ni_2, nj_2, x_in_2, &
     &                    max_ni, max_nj)
            call copy_arr(temp_arr_3, ni_2, nj_2, y_in_2, &
     &                    max_ni, max_nj)

! calculate the weightings from wind_set_2 to elcirc nodes
! (this is slow)
            write(*,*)
            write(*,*) &
     &        'begin calculating grid weightings for wind_file_2'
            write(16,*)
            write(16,*) &
     &        'begin calculating grid weightings for wind_file_2'
            call get_weight (x_in_2, y_in_2, x, y, &
     &                       elem_nodes_in_2, node_i_2, node_j_2, &
     &                       max_ni, max_nj, &
     &                       num_elems_in_2, num_nodes_in_2, &
     &                       num_nodes, &
     &                       max_nodes, &
     &                       in_elem_to_out_node_2, &
     &                       temp_arr_6, weight_wind_node_2)
            write(*,*) &
     &        'done calculating grid weightings for wind_file_2'
            write(16,*) &
     &        'done calculating grid weightings for wind_file_2'

          endif

! output starting date and time
          write(*,*)
          write(*,*) 'wind file starting Julian date: ', start_day_1
          write(*,*) 'wind file assumed UTC starting time: ', utc_start
          write(16,*)
          write(16,*) 'wind file starting Julian date: ', start_day_1
          write(16,*) 'wind file assumed UTC starting time: ', utc_start

        endif ! (end of first_call block)

! define frac_day - the fractional Julian date
! include offset for starting time in UTC (in hours)
        frac_day = start_day_1 + time/secs_per_day + utc_start/24.0

! output info to debug file
        write(39,*) 'num_nodes = ', num_nodes
        write(39,*) 'num_files_1 = ', num_files_1
        write(39,*) 'num_times_1 = ', num_times_1
        if (have_wind_2) then
          write(39,*) 'num_files_2 = ', num_files_2
          write(39,*) 'num_times_2 = ', num_times_2
        endif
        write(39,*)
        write(39,*) 'start_day = ', start_day_1
        write(39,*) 'time      = ', time
        write(39,*) 'frac_day  = ', frac_day

! now retrieve data from the wind files

        call combine_set (real(frac_day), 'u                   ', &
     &                    wind_set_1, wind_set_2, &
     &                    temp_arr_2, temp_arr_3, &
     &                    temp_arr_4, temp_arr_5, temp_arr_1, &
     &                    wind_times_1, num_times_1, ni_1, nj_1, &
     &                    wind_times_2, num_times_2, ni_2, nj_2, &
     &                    num_files_1, num_files_2, max_files, &
     &                    temp_arr_8, temp_arr_9, u_air_node, &
     &                    in_elem_to_out_node_1, weight_wind_node_1, &
     &                    in_elem_to_out_node_2, weight_wind_node_2, &
     &                    elem_nodes_in_1, node_i_1, node_j_1, &
     &                    elem_nodes_in_2, node_i_2, node_j_2, &
     &                    max_ni, max_nj, max_nodes_in, &
     &                    num_nodes, max_nodes, &
     &                    num_nodes_in_1, num_elems_in_1, &
     &                    num_nodes_in_2, num_elems_in_2, &
     &                    max_elems_in, have_wind_2, &
     &                    relative_weight_1, relative_weight_2, &
     &                    max_window_1, max_window_2, &
     &                    wind_time_files_1, wind_time_files_2)

        call combine_set (real(frac_day), 'v                   ', &
     &                    wind_set_1, wind_set_2, &
     &                    temp_arr_2, temp_arr_3, &
     &                    temp_arr_4, temp_arr_5, temp_arr_1, &
     &                    wind_times_1, num_times_1, ni_1, nj_1, &
     &                    wind_times_2, num_times_2, ni_2, nj_2, &
     &                    num_files_1, num_files_2, max_files, &
     &                    temp_arr_8, temp_arr_9, v_air_node, &
     &                    in_elem_to_out_node_1, weight_wind_node_1, &
     &                    in_elem_to_out_node_2, weight_wind_node_2, &
     &                    elem_nodes_in_1, node_i_1, node_j_1, &
     &                    elem_nodes_in_2, node_i_2, node_j_2, &
     &                    max_ni, max_nj, max_nodes_in, &
     &                    num_nodes, max_nodes, &
     &                    num_nodes_in_1, num_elems_in_1, &
     &                    num_nodes_in_2, num_elems_in_2, &
     &                    max_elems_in, have_wind_2, &
     &                    relative_weight_1, relative_weight_2, &
     &                    max_window_1, max_window_2, &
     &                    wind_time_files_1, wind_time_files_2)

        call combine_set (real(frac_day), 'p_msl               ', &
     &                    wind_set_1, wind_set_2, &
     &                    temp_arr_2, temp_arr_3, &
     &                    temp_arr_4, temp_arr_5, temp_arr_1, &
     &                    wind_times_1, num_times_1, ni_1, nj_1, &
     &                    wind_times_2, num_times_2, ni_2, nj_2, &
     &                    num_files_1, num_files_2, max_files, &
     &                    temp_arr_8, temp_arr_9, p_air_node, &
     &                    in_elem_to_out_node_1, weight_wind_node_1, &
     &                    in_elem_to_out_node_2, weight_wind_node_2, &
     &                    elem_nodes_in_1, node_i_1, node_j_1, &
     &                    elem_nodes_in_2, node_i_2, node_j_2, &
     &                    max_ni, max_nj, max_nodes_in, &
     &                    num_nodes, max_nodes, &
     &                    num_nodes_in_1, num_elems_in_1, &
     &                    num_nodes_in_2, num_elems_in_2, &
     &                    max_elems_in, have_wind_2, &
     &                    relative_weight_1, relative_weight_2, &
     &                    max_window_1, max_window_2, &
     &                    wind_time_files_1, wind_time_files_2)

        call combine_set (real(frac_day), 'T                   ', &
     &                    wind_set_1, wind_set_2, &
     &                    temp_arr_2, temp_arr_3, &
     &                    temp_arr_4, temp_arr_5, temp_arr_1, &
     &                    wind_times_1, num_times_1, ni_1, nj_1, &
     &                    wind_times_2, num_times_2, ni_2, nj_2, &
     &                    num_files_1, num_files_2, max_files, &
     &                    temp_arr_8, temp_arr_9, t_air_node, &
     &                    in_elem_to_out_node_1, weight_wind_node_1, &
     &                    in_elem_to_out_node_2, weight_wind_node_2, &
     &                    elem_nodes_in_1, node_i_1, node_j_1, &
     &                    elem_nodes_in_2, node_i_2, node_j_2, &
     &                    max_ni, max_nj, max_nodes_in, &
     &                    num_nodes, max_nodes, &
     &                    num_nodes_in_1, num_elems_in_1, &
     &                    num_nodes_in_2, num_elems_in_2, &
     &                    max_elems_in, have_wind_2, &
     &                    relative_weight_1, relative_weight_2, &
     &                    max_window_1, max_window_2, &
     &                    wind_time_files_1, wind_time_files_2)

        call combine_set (real(frac_day), 'qv                  ', &
     &                    wind_set_1, wind_set_2, &
     &                    temp_arr_2, temp_arr_3, &
     &                    temp_arr_4, temp_arr_5, temp_arr_1, &
     &                    wind_times_1, num_times_1, ni_1, nj_1, &
     &                    wind_times_2, num_times_2, ni_2, nj_2, &
     &                    num_files_1, num_files_2, max_files, &
     &                    temp_arr_8, temp_arr_9, q_air_node, &
     &                    in_elem_to_out_node_1, weight_wind_node_1, &
     &                    in_elem_to_out_node_2, weight_wind_node_2, &
     &                    elem_nodes_in_1, node_i_1, node_j_1, &
     &                    elem_nodes_in_2, node_i_2, node_j_2, &
     &                    max_ni, max_nj, max_nodes_in, &
     &                    num_nodes, max_nodes, &
     &                    num_nodes_in_1, num_elems_in_1, &
     &                    num_nodes_in_2, num_elems_in_2, &
     &                    max_elems_in, have_wind_2, &
     &                    relative_weight_1, relative_weight_2, &
     &                    max_window_1, max_window_2, &
     &                    wind_time_files_1, wind_time_files_2)

! set first_call to false, so subsequent calls will know that they're
! not the first call
        first_call = .false.

        close(39)

      return
      end
!-----------------------------------------------------------------------
!
! Note that net downward heat flux (except for solar) is given by:
!
! net_sfc_flux_d = - (sen_flux + lat_flux + (longwave_u - longwave_d))
!
      subroutine surf_fluxes (time, u_air, v_air, p_air, &
     &                   t_air, q_air, shortwave_d, &
     &                   sen_flux, lat_flux, longwave_u, longwave_d, &
     &                   tau_xz, tau_yz)

!       implicit none
        use global
        implicit real*8(a-h,o-z),integer(i-n)

! define some new names for things in header file
        integer max_nodes
        parameter (max_nodes = mnp)

! input/output variables
        real*8 u_air(max_nodes), v_air(max_nodes)
        real*8 t_air(max_nodes), q_air(max_nodes)
        real*8 p_air(max_nodes), time
        real*8 sen_flux(max_nodes), lat_flux(max_nodes)
        real*8 longwave_d(max_nodes), longwave_u(max_nodes)
        real*8 shortwave_d(max_nodes)
        real*8 tau_xz(max_nodes), tau_yz(max_nodes)

! local variables
        integer num_nodes, i_node
        real*4 temp_sca, window
        real*8 weight_rad_node(max_nodes,3)
        real*8 start_day, one_third, frac_day, secs_per_day, utc_start
        real*8 max_window
        parameter (max_window = 6.0)
        parameter (secs_per_day = 86400.0)
        parameter (utc_start = 8.0)
        integer in_elem_to_out_node(max_nodes)
        character flux_set*50, start_day_file*50
        parameter (flux_set = 'hdf/flux_file_1')
        parameter (start_day_file = 'hdf/start_day.txt')
        logical first_call, time_exst, have_start_day_file
        data first_call/.true./
        integer max_ni, max_nj, ni, nj, max_times, num_times
        integer num_files, max_files
        integer num_nodes_in, num_elems_in
        integer max_elems_in, max_nodes_in
        parameter (max_ni = 1024)
        parameter (max_nj = 1024)
        parameter (max_elems_in = (max_ni-1) * (max_nj-1) * 2)
        parameter (max_nodes_in = max_ni * max_nj)
        parameter (max_times = 1000)
        parameter (max_files = 999)
        real*4 temp_arr_2(max_ni,max_nj), temp_arr_3(max_ni,max_nj)
        real*4 temp_arr_4(max_ni,max_nj), temp_arr_5(max_ni,max_nj)
        real*8 temp_arr_6(max_elems_in)
        real*8 temp_arr_8(max_nodes)
        real*4 flux_times(max_times)
        real*8 shortwave_in(max_ni,max_nj), longwave_in(max_ni,max_nj)
        real*8 x_in(max_ni,max_nj), y_in(max_ni,max_nj)
        real*8 stefan, t_freeze, emissivity
        character flux_time_files(max_times)*50
        integer node_i(max_nodes_in), node_j(max_nodes_in)
        integer max_rank, rank
        integer node_num_in(max_nodes_in), elem_nodes_in(max_elems_in,3)
        parameter (max_rank = 3)
        integer dim_sizes(max_rank)
        parameter (stefan = 5.67e-8)
        parameter (t_freeze = 273.15)
        parameter (emissivity = 1.0)
        integer printit

! retain the values of some local variables between calls
        save first_call, start_day, in_elem_to_out_node, &
     &       weight_rad_node, num_nodes, &
     &       flux_times, num_times, num_files, ni, nj, &
     &       num_nodes_in, num_elems_in, node_i, node_j, node_num_in, &
     &       elem_nodes_in, flux_time_files

        printit = 100

        open(38,file='fort.38')
        rewind(38)
        write(38,*)
        write(38,*) 'enter surf_fluxes'
        write(38,*) 'first_call = ', first_call

! define some constants, etc
        one_third = 1.0 / 3.0

! if this is the first call to this routine then get some things ready
        if (first_call) then

! define the local variables num_nodes
          num_nodes = np

! check to see if start_day_file exists
          call file_exst (start_day_file, have_start_day_file, .false.)

! if the start day file does exist, get start_day from it, otherwise
! use the first start_day in flux_set
          if (have_start_day_file) then
            open (unit=77, file=start_day_file, status='old')
            read(77,*) temp_sca
            close (unit=77)
            start_day = temp_sca
          else
            flux_time_files(1) = 'hdf/flux_file_1.001.hdf'
            call read_scalar(flux_time_files(1), temp_sca, &
     &                       'start_day           ', 0.0)
            start_day = temp_sca
          endif

! get the times of the fluxes available in the flux_set
          call get_times(flux_times, flux_set, &
     &                   'shortwave_flux      ', &
     &                   flux_time_files, num_times, &
     &                   num_files, max_times, max_files)

! get the dimensions of the flux datasets (use first dataset)
          call get_dims(flux_time_files(1), 'shortwave_flux      ', &
     &                  flux_times(1), rank, dim_sizes)
          ni = dim_sizes(1)
          nj = dim_sizes(2)

! check the dimensions of the flux datasets, to ensure they don't exceed
! the maximums
          if (ni .gt. max_ni .or. nj .gt. max_nj) then
            write(*,*)
            write(*,*) 'Radiative flux data: max dimensions exceeded!'
            write(11,*)
            write(11,*) 'Radiative flux data: max dimensions exceeded!'
            stop
          endif

! calculate the total number of nodes and elements for the radiative
! fluxes input grid
          num_nodes_in = ni * nj
          num_elems_in = (ni-1) * (nj-1) * 2

! check the elems/nodes of the input dataset, to ensure they don't
! exceed the maximums
          if (num_elems_in .gt. max_elems_in .or. &
     &        num_nodes_in .gt. max_nodes_in) then
            write(*,*)
            write(*,*) 'Radiative flux data: max elems/nodes exceeded!'
            write(11,*)
            write(11,*) 'Radiative flux data: max elems/nodes exceeded!'
            stop
          endif

! create list of all nodes for input grid
          call list_nodes (node_i, node_j, node_num_in, &
     &                     num_nodes_in, ni, nj)

! now create the list of all the elements (and the nodes defining them)
          call list_elems (elem_nodes_in, node_num_in, &
     &                     ni, nj, num_elems_in)

! read in the x and y values for the input grid, and copy to full size
! real*8 arrays
          call read_2d_arr(flux_time_files(1), temp_arr_4, &
     &                     'x                   ', 0.0, &
     &                     ni, nj)
          call read_2d_arr(flux_time_files(1), temp_arr_5, &
     &                     'y                   ', 0.0, &
     &                     ni, nj)
          call copy_arr(temp_arr_4, ni, nj, x_in, max_ni, max_nj)
          call copy_arr(temp_arr_5, ni, nj, y_in, max_ni, max_nj)

! calculate the weightings from rad grid to elcirc nodes (this is slow)
          write(*,*)
          write(*,*) 'calculating grid weightings for rad fluxes'
          write(16,*)
          write(16,*) 'calculating grid weightings for rad fluxes'
          call get_weight (x_in, y_in, x, y, &
     &                     elem_nodes_in, node_i, node_j, &
     &                     max_ni, max_nj, &
     &                     num_elems_in, num_nodes_in, &
     &                     num_nodes, &
     &                     max_nodes, &
     &                     in_elem_to_out_node, &
     &                     temp_arr_6, weight_rad_node)

! output starting date and time
          write(*,*)
          write(*,*) 'rad fluxes file starting Julian date: ', &
     &               start_day
          write(*,*) 'rad fluxes file assumed UTC starting time: ', &
     &               utc_start
          write(16,*)
          write(16,*) 'rad fluxes file starting Julian date: ', &
     &                start_day
          write(16,*) 'rad fluxes file assumed UTC starting time: ', &
     &                utc_start

        endif ! (end of first_call block)

! define frac_day - the fractional Julian date
! include offset for starting time in UTC (in hours)
        frac_day = start_day + time/secs_per_day + utc_start/24.0

! output info to debug file
        write(38,*) 'num_nodes = ', num_nodes
        write(38,*) 'num_files = ', num_files
        write(38,*) 'num_times = ', num_times
        write(38,*)
        write(38,*) 'start_day = ', start_day
        write(38,*) 'time      = ', time
        write(38,*) 'frac_day  = ', frac_day

! convert air temperatures to Celcius
        do i_node = 1, num_nodes
          t_air(i_node) = t_air(i_node) - t_freeze

          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*)
            write(38,*) 'i_node, sfc u, v, T = ', i_node, &
     &                  uu2(i_node, kfp(i_node)), &
     &                  vv2(i_node, kfp(i_node)), &
     &                  tnd(i_node, kfp(i_node))
            write(38,*) 'u, v, p, T, q (air) = ', u_air(i_node), &
     &                  v_air(i_node), p_air(i_node), t_air(i_node), &
     &                  q_air(i_node)
          endif

        enddo

! calculate the turbulent fluxes at the nodes
        write(38,*) 'above turb_fluxes'
        call turb_fluxes (num_nodes, &
     &                    u_air, v_air, p_air, t_air, q_air, &
     &                    sen_flux, lat_flux, tau_xz, tau_yz)
        write(38,*) 'below turb_fluxes'

! now we must retrieve the radiative fluxes, and copy from the filled
! real*4 arrays of reduced size to full-size real*8 arrays

! first shortwave
        call get_set (real(frac_day), 'shortwave_flux      ', &
     &                flux_set, time_exst, &
     &                temp_arr_4, temp_arr_2, temp_arr_3, &
     &                flux_times, num_times, ni, nj, &
     &                num_files, max_files, window, &
     &                flux_time_files)
        call copy_arr(temp_arr_4, ni, nj, shortwave_in, max_ni, max_nj)

! next longwave
        call get_set (real(frac_day), 'longwave_flux       ', &
     &                flux_set, time_exst, &
     &                temp_arr_4, temp_arr_2, temp_arr_3, &
     &                flux_times, num_times, ni, nj, &
     &                num_files, max_files, window, &
     &                flux_time_files)
        call copy_arr(temp_arr_4, ni, nj, longwave_in, max_ni, max_nj)
!       write(*,*) 'time_exst, window = ', time_exst, window

! if the desired time is outside of the times supplied by the dataset,
! or if the data interpolation window is greater than the maximum, then
! quit with an error message
        if ( (.not. time_exst) .or. (window .gt. max_window) ) then
          write(*,*)
          write(*,*) 'Radiative flux data: time not present in file!'
          write(11,*)
          write(11,*) 'Radiative flux data: time not present in file!'
          stop
        endif

        write(38,*) 'rad fluxes obtained'

! now the input radiative fluxes must be put on a triangular grid

! use the weightings to interpolate from the input grid to the node
! points
        write(38,*) 'interpolating longwave'
        call interp_data (longwave_d, longwave_in, weight_rad_node, &
     &                    in_elem_to_out_node, elem_nodes_in, &
     &                    node_i, node_j, &
     &                    num_nodes_in, num_nodes, num_elems_in, &
     &                    max_ni, max_nj, max_nodes)

        write(38,*) 'interpolating shortwave'
        call interp_data (shortwave_d, shortwave_in, weight_rad_node, &
     &                    in_elem_to_out_node, elem_nodes_in, &
     &                    node_i, node_j, &
     &                    num_nodes_in, num_nodes, num_elems_in, &
     &                    max_ni, max_nj, max_nodes)

! now calculate upwards longwave flux at the surface, using black-body
! equation
        write(38,*) 'calculating longwave_u'
        do i_node = 1, num_nodes
          longwave_u(i_node) = emissivity * stefan * &
     &                      ( t_freeze + tnd(i_node, kfp(i_node)) ) ** 4
        enddo

! get the albedo (store in temp_arr_8)
        write(38,*) 'calculating albedo'
        call get_albedo (temp_arr_8, num_nodes, max_nodes)

! reduce the downwards shortwave flux at the surface by the albedo
! (and check for possible negative values from interpolation)
        write(38,*) 'reducing shortwave'
        do i_node = 1, num_nodes
          shortwave_d(i_node) = max( (1.0 - temp_arr_8(i_node)) * &
     &                               shortwave_d(i_node), &
     &                               0.0 )
        enddo

        do i_node = 1, num_nodes
          if (mod(i_node-1,printit) .eq. 0) then
            if (kfp(i_node) .ne. 0) then
              write(38,*)
              write(38,*) 'i_node = ', i_node
              write(38,*) 'net_sfc_flux_d = ', &
     &                     - sen_flux(i_node) - lat_flux(i_node) &
     &                     - ( longwave_u(i_node) - longwave_d(i_node) )
              write(38,*) 'shortwave_d = ', shortwave_d(i_node)
              write(38,*) 'longwave_d, longwave_u = ', &
     &                     longwave_d(i_node), longwave_u(i_node)
              write(38,*) 'sen_flux, lat_flux = ', &
     &                     sen_flux(i_node), lat_flux(i_node)
            else
              write(38,*)
              write(38,*) 'i_node = ', i_node
              write(38,*) 'dry!'
            endif
          endif
        enddo

! set first_call to false, so subsequent calls will know that they're
! not the first call
        first_call = .false.

        close(38)

      return
      end
!-----------------------------------------------------------------------
!
! Calculate bulk aerodynamic surface fluxes over water using method of
! Zeng et al., J Clim, v11, p 2628, Oct 1998.
!
! Note: this subroutine uses actual temperatures instead of potential
!       temperatures.  Since the temperature height is only 2m, the
!       difference should be negligible. . .
!
      subroutine turb_fluxes (num_nodes, &
     &                        u_air, v_air, p_air, t_air, q_air, &
     &                        sen_flux, lat_flux, tau_xz, tau_yz)

!       implicit none
        use global
        implicit real*8(a-h,o-z),integer(i-n)

! define some new names for things in header file
        integer max_nodes
        parameter (max_nodes = mnp)

! input/output variables
        integer num_nodes
        real*8 u_air(max_nodes), v_air(max_nodes), p_air(max_nodes)
        real*8 t_air(max_nodes), q_air(max_nodes)
        real*8 sen_flux(max_nodes), lat_flux(max_nodes)
        real*8 tau_xz(max_nodes), tau_yz(max_nodes)

! local variables
        integer i_node, iter, max_iter
        parameter (max_iter = 10)
        real*8 u_star, theta_star, q_star, z_0, monin
        real*8 zeta_u, zeta_t
        real*8 one_third, two_thirds, z_t, z_u
        real*8 a1, a2, b1, b2, nu, g, beta, w_star, mix_ratio
        real*8 t_v, z_i, speed, karman, zeta_m, psi_m, zeta_h, psi_h
        real*8 re, z_0_t, e_sfc, t_freeze, q_sfc, epsilon_r
        real*8 rho_air, c_p_air, latent, r_air, rb
        real*8 theta_air, theta_v_air, delta_theta, delta_q
        real*8 delta_theta_v, theta_v_star, speed_res, tau
        real*4 esat_flat_r
        parameter (z_t = 2.0)
        parameter (z_u = 10.0)
        parameter (a1 = 0.013)
        parameter (a2 = 0.11)
        parameter (b1 = 2.67)
        parameter (b2 = -2.57)
        parameter (nu = 1.46e-5)
        parameter (beta = 1.0)
        parameter (g = 9.81)
        parameter (z_i = 1000.0)
        parameter (karman = 0.4)
        parameter (zeta_m = -1.574)
        parameter (zeta_h = -0.465)
        parameter (t_freeze = 273.15)
        parameter (epsilon_r = 0.6220)
        parameter (c_p_air = 1004.0)
        parameter (latent = 2.501e6)
        parameter (r_air = 287.0)
        integer printit
        logical converged
        printit = 100

        write(38,*) 'enter turb_fluxes'

! precalculate constants
        one_third = 1.0 / 3.0
        two_thirds = 2.0 / 3.0

! now loop over all points
        do i_node = 1, num_nodes
          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*)
            write(38,*) 'i_node = ', i_node
          endif

! if this point isn't dry, then calculate fluxes (if dry, then skip)
        if (kfp(i_node) .ne. 0) then

! calculate q_sfc from e_sfc
! (e_sfc reduced for salinity using eqn from Smithsonian Met Tables)
          e_sfc = (1.0 - 0.000537 * snd(i_node, kfp(i_node))) &
     &          * esat_flat_r(real(tnd(i_node, kfp(i_node)) + t_freeze))
          q_sfc = epsilon_r * e_sfc &
     &          / ( p_air(i_node) - e_sfc * (1.0 - epsilon_r) )

! calculate the water vapor mixing ratio of the air
          mix_ratio = q_air(i_node) / (1.0 - q_air(i_node))

! calculate theta_air, theta_v_air, delta_theta, delta_q,
! and delta_theta_v
          theta_air = (t_air(i_node) + t_freeze) + 0.0098*z_t
          theta_v_air = theta_air * (1.0 + 0.608 * mix_ratio)
          delta_theta = theta_air - &
     &                  (tnd(i_node, kfp(i_node)) + t_freeze)
          delta_q = q_air(i_node) - q_sfc
          delta_theta_v = delta_theta * (1.0 + 0.608 * mix_ratio) &
     &                  + 0.608 * theta_air * delta_q

! calculate the air virtual temperature and density
          t_v = (t_air(i_node) + t_freeze) * (1.0 + 0.608 * mix_ratio)
          rho_air = p_air(i_node) / (r_air * t_v)
          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*) 'e_sfc, q_sfc, mix_ratio = ', &
     &                   e_sfc, q_sfc, mix_ratio
            write(38,*) 'theta_air, theta_v_air, delta_theta = ', &
     &                   theta_air, theta_v_air, delta_theta
            write(38,*) 'delta_q, delta_theta_v, t_v = ', &
     &                   delta_q, delta_theta_v, t_v
            write(38,*) 'rho_air = ', rho_air
          endif

! begin with initial values of u_star, w_star, and speed
          u_star = 0.06
          w_star = 0.5
          if (delta_theta_v .ge. 0) then                    ! stable
            speed = &
     &        max( sqrt( &
     &               (u_air(i_node) - uu2(i_node, kfp(i_node)))**2 + &
     &               (v_air(i_node) - vv2(i_node, kfp(i_node)))**2 ), &
     &             0.1)
          else                                              ! unstable
            speed = &
     &        sqrt( (u_air(i_node) - uu2(i_node, kfp(i_node)))**2 + &
     &              (v_air(i_node) - vv2(i_node, kfp(i_node)))**2 + &
     &              (beta * w_star)**2 )
          endif

! now loop to obtain good initial values for u_star and z_0
          do iter = 1, 5
            z_0 = a1 * u_star * u_star / g + a2 * nu / u_star
            u_star = karman * speed / log(z_u/z_0)
          enddo

! calculate rb (some stability parameter from Xubin's code?)
          rb = g * z_u * delta_theta_v / (theta_v_air * speed * speed)

! calculate initial values for zeta_u, monin, zeta_t
          if (rb .ge. 0) then                      ! neutral or stable
            zeta_u = rb * log(z_u/z_0) &
     &             / (1.0 - 0.5*min(rb,0.19))
          else
            zeta_u = rb * log(z_u/z_0)
          endif
          monin = z_u / zeta_u
          zeta_t = z_t / monin
          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*) 'speed, z_0, u_star = ', &
     &                   speed, z_0, u_star
            write(38,*) 'zeta_u, zeta_t, monin = ', &
     &                   zeta_u, zeta_t, monin
          endif

! iterate a maximum of max_iter times
          iter = 0
          converged = .false.
100       continue
            iter = iter + 1

! Calculate the roughness lengths
            z_0 = a1 * u_star * u_star / g + a2 * nu / u_star
            re = u_star * z_0 / nu
            z_0_t = z_0 / exp(b1 * (re**0.25) + b2)
            if (mod(i_node-1,printit) .eq. 0) then
              write(38,*) 're, z_0, z_0_t = ', &
     &                     re, z_0, z_0_t
            endif

! calculate the zetas
            zeta_u = z_u / monin
            zeta_t = z_t / monin

! apply asymptotic limit to stable conditions
            if (zeta_t .gt. 2.5) then
              converged = .true.
              zeta_t = 2.5
              monin = z_t / zeta_t
              zeta_u = z_u / monin
              if (mod(i_node-1,printit) .eq. 0) then
                write(38,*) 'limiting zeta_u, zeta_t, monin!'
              endif
            endif

! caulculate u_star, depending on zeta
            if (zeta_u .lt. zeta_m) then                ! very unstable

              u_star = speed * karman &
     &               / ( log(zeta_m*monin/z_0) &
     &                   - psi_m(zeta_m) &
! extra term?
     &                   + psi_m(z_0/monin) &
     &                   + 1.14 * ((-zeta_u)**(one_third) - &
     &                           (-zeta_m)**(one_third)) )

            else if (zeta_u .lt. 0.0) then              ! unstable

              u_star = speed * karman &
     &               / ( log(z_u/z_0) &
     &                   - psi_m(zeta_u) &
! extra term?
     &                   + psi_m(z_0/monin) &
     &                 )

            else if (zeta_u .le. 1.0) then              ! neutral/stable

              u_star = speed * karman &
     &               / ( log(z_u/z_0) + 5.0*zeta_u &
! extra term?
     &                   - 5.0*z_0/monin &
     &                 )

            else                                        ! very stable

              u_star = speed * karman &
     &               / ( log(monin/z_0) + 5.0 &
     &                   + 5.0*log(zeta_u) &
! extra term?
     &                   - 5.0*z_0/monin &
     &                   + zeta_u - 1.0 )

            endif

! caulculate theta_star and q_star, depending on zeta
            if (zeta_t .lt. zeta_h) then                ! very unstable

              theta_star = karman * delta_theta &
     &                   / ( log(zeta_h*monin/z_0_t) &
     &                       - psi_h(zeta_h) &
! extra term?
     &                       + psi_h(z_0_t/monin) &
     &                       + 0.8 * ((-zeta_h)**(-one_third) - &
     &                                (-zeta_t)**(-one_third)) )

              q_star = karman * delta_q &
     &                   / ( log(zeta_h*monin/z_0_t) &
     &                       - psi_h(zeta_h) &
! extra term?
     &                       + psi_h(z_0_t/monin) &
     &                       + 0.8 * ((-zeta_h)**(-one_third) - &
     &                                (-zeta_t)**(-one_third)) )

            else if (zeta_t .lt. 0.0) then              ! unstable

              theta_star = karman * delta_theta &
     &                   / ( log(z_t/z_0_t) &
     &                       - psi_h(zeta_t) &
! extra term?
     &                       + psi_h(z_0_t/monin) &
     &                     )

              q_star = karman * delta_q &
     &                   / ( log(z_t/z_0_t) &
     &                       - psi_h(zeta_t) &
! extra term?
     &                       + psi_h(z_0_t/monin) &
     &                     )

            else if (zeta_t .lt. 1.0) then              ! neutral/stable

              theta_star = karman * delta_theta &
     &                   / ( log(z_t/z_0_t) &
     &                       + 5.0*zeta_t &
! extra term?
     &                       - 5.0*z_0_t/monin &
     &                     )

              q_star = karman * delta_q &
     &                   / ( log(z_t/z_0_t) &
     &                       + 5.0*zeta_t &
! extra term?
     &                       - 5.0*z_0_t/monin &
     &                     )

            else                                        ! very stable

              theta_star = karman * delta_theta &
     &                   / ( log(monin/z_0_t) + 5.0 &
     &                       + 5.0*log(zeta_t) &
! extra term?
     &                       - 5.0*z_0_t/monin &
     &                       + zeta_t - 1.0 )

              q_star = karman * delta_q &
     &                   / ( log(monin/z_0_t) + 5.0 &
     &                       + 5.0*log(zeta_t) &
! extra term?
     &                       - 5.0*z_0_t/monin &
     &                       + zeta_t - 1.0 )

            endif

! calculate theta_v_star and monin
            theta_v_star = theta_star * (1.0 + 0.608 * mix_ratio) &
     &                   + 0.608 * theta_air * q_star
            monin = theta_v_air * u_star * u_star &
     &            / (karman * g * theta_v_star)

! depending on surface layer stability, calculate the effective
! near-surface wind speed
! (ie relative to the flowing water surface)
            if (delta_theta_v .ge. 0.0) then                  ! stable

              speed = &
     &          max( sqrt( &
     &                 (u_air(i_node) - uu2(i_node, kfp(i_node)))**2 + &
     &                 (v_air(i_node) - vv2(i_node, kfp(i_node)))**2 ), &
     &               0.1)

            else                                              ! unstable

! calculate the convective velocity scale
              w_star = (-g * theta_v_star * u_star * z_i / theta_v_air) &
     &                 ** one_third

              speed = &
     &          sqrt( (u_air(i_node) - uu2(i_node, kfp(i_node)))**2 + &
     &                (v_air(i_node) - vv2(i_node, kfp(i_node)))**2 + &
     &                (beta * w_star)**2 )

            endif
            if (mod(i_node-1,printit) .eq. 0) then
              write(38,*) 'iter, u_star, q_star, theta_star = ', &
     &                     iter, u_star, q_star, theta_star
              write(38,*) 'iter, theta_v_star, monin, speed = ', &
     &                     iter, theta_v_star, monin, speed
              write(38,*) 'iter, zeta_u, zeta_t = ', &
     &                     iter, zeta_u, zeta_t
            endif

! bottom of main iteration loop
          if (.not. converged .and. iter .lt. max_iter) goto 100

! calculate fluxes
          sen_flux(i_node) = - rho_air * c_p_air * u_star * theta_star
          lat_flux(i_node) = - rho_air * latent * u_star * q_star

! calculate wind stresses
          speed_res = &
     &          sqrt( (u_air(i_node) - uu2(i_node, kfp(i_node)))**2 + &
     &                (v_air(i_node) - vv2(i_node, kfp(i_node)))**2 )
          if (speed_res .gt. 0.0) then
            tau = rho_air * u_star * u_star * speed_res / speed
            tau_xz(i_node) = - tau &
     &                     * (u_air(i_node) - uu2(i_node, kfp(i_node))) &
     &                     / speed_res
            tau_yz(i_node) = - tau &
     &                     * (v_air(i_node) - vv2(i_node, kfp(i_node))) &
     &                     / speed_res
          else
            tau_xz(i_node) = 0.0
            tau_yz(i_node) = 0.0
          endif

          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*) 'sen_flux, lat_flux = ', &
     &                   sen_flux(i_node), lat_flux(i_node)
            write(38,*) 'tau_xz, tau_yz = ', &
     &                   tau_xz(i_node), tau_yz(i_node)
          endif
!         write(38,*) tnd(i_node, kfp(i_node)), sen_flux(i_node),
!    +                lat_flux(i_node)

! end of wet/dry block
        endif

! end of loop over points
        enddo

        write(38,*) 'exit turb_fluxes'

      return
      end
!-----------------------------------------------------------------------
!
! Calculate saturation vapor pressure using the eighth order relative
! error norm method of Flatau et al., J Applied Meteo, v31, p 1507,
! Dec 1992.
!
      real*4 function esat_flat_r(t)
        implicit none
        real*4 t
        real*4 c0, c1, c2, c3, c4, c5, c6, c7, c8, t_eff
        parameter ( &
     &        c0= 6.11583699e+02,  c1= 0.444606896e+02, &
     &        c2= 0.143177157e+01, c3= 0.264224321e-01, &
     &        c4= 0.299291081e-03, c5= 0.203154182e-05, &
     &        c6= 0.702620698e-08, c7= 0.379534310e-11, &
     &        c8=-0.321582393e-13)

! t     : temperature in K
! t_eff : effective temperature in C

        t_eff = max(-85.,t-273.16)

        esat_flat_r = c0+t_eff*(c1+t_eff*(c2+t_eff*(c3+t_eff*(c4+t_eff* &
     &                         (c5+t_eff*(c6+t_eff*(c7+t_eff*c8)))))))

      return
      end
!-----------------------------------------------------------------------
      subroutine copy_arr(array_1, ni_1, nj_1, array_2, ni_2, nj_2)
        implicit none
        integer ni_1, nj_1, ni_2, nj_2, i, j
        real*4 array_1(ni_1,nj_1)
        real*8 array_2(ni_2,nj_2)

        do j = 1, nj_1
          do i = 1, ni_1
            array_2(i,j) = array_1(i,j)
          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine file_exst (file_name, exst, fail)
        implicit none
        character file_name*50
        logical exst, fail

        inquire(file=file_name, exist=exst)
        if (exst .and. fail) then
          write(*,*)
          write(*,*) file_name, ' already exists!'
          write(11,*)
          write(11,*) file_name, ' already exists!'
         stop
        endif

      return
      end
!-----------------------------------------------------------------------
      subroutine get_bracket (time, input_times, i_time, bracket, &
     &                        num_times)
        implicit none
        integer num_times, i_time
        real*4 time, input_times(num_times)
        logical bracket

        i_time = 0
        bracket = .false.
10      continue
          i_time = i_time + 1
          bracket = ( input_times(i_time) .le. time .and. &
     &                time .le. input_times(i_time+1) )
        if (.not. bracket .and. i_time .lt. num_times-1) goto 10

      return
      end
!-----------------------------------------------------------------------
      subroutine interp_set(in_set, data, data_label, time, &
     &                      data_1, data_2, input_times, &
     &                      i_time, num_times, ni, nj, &
     &                      time_files, num_files, max_files)
        implicit none
        integer ni, nj, i, j, num_times, i_time, num_files, max_files
        real*4 data_1(ni,nj), data_2(ni,nj)
        real*4 data(ni,nj), time, input_times(num_times), ratio
        character in_set*50, data_label*20
        character time_files(num_times)*50

! read in the data at the first of the bracketing times
        call read_2d_arr (time_files(i_time), data_1, data_label, &
     &                    input_times(i_time), ni, nj)

! read in the data at the second of the bracketing times
        call read_2d_arr (time_files(i_time+1), data_2, data_label, &
     &                    input_times(i_time+1), ni, nj)

! now interpolate these to the desired time
        ratio = (time                  - input_times(i_time)) &
     &        / (input_times(i_time+1) - input_times(i_time))
        do j = 1, nj
          do i = 1, ni
            data(i,j) = data_1(i,j) &
     &                + (data_2(i,j) - data_1(i,j)) * ratio
          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine get_dims(in_file, data_label, t, rank, dim_sizes)
        implicit none
        real*4 t
        integer ret, read_only, sfstart, sfend, sd_id, sds_id
        integer sds_index, sfn2index, sfendacc, sfselect
        integer max_rank, rank, data_type, n_attrs, sfginfo
        parameter (max_rank = 3)
        parameter (read_only = 1)
        integer dim_sizes(max_rank)
        character data_label*20, in_file*50
        character dat_time_label*33, sds_name*33

! open in_file in read only mode
        sd_id = sfstart(in_file, read_only)

! create the data-time label, which we'll use as the name
        call get_label (dat_time_label, data_label, t)

! find index for this data set
        if (t .lt. 0.0) then
          sds_index = sfn2index(sd_id, data_label)
        else
          sds_index = sfn2index(sd_id, dat_time_label)
        endif

! get info on dataset
        sds_id = sfselect(sd_id, sds_index)
        ret = sfginfo(sds_id, sds_name, rank, dim_sizes, &
     &                data_type, n_attrs)
        call checkret(ret)

! find the id for this index
        sds_id = sfselect(sd_id, sds_index)

! close access to the dataset
        ret = sfendacc(sds_id)
        call checkret(ret)

! close access to the file
        ret = sfend(sd_id)
        call checkret(ret)

      return
      end
!-----------------------------------------------------------------------
      subroutine list_nodes (node_i, node_j, node_num, &
     &                       n_nodes_in, ni, nj)
        implicit none
        integer ni, nj, n_nodes_in, i, j, i_node
        integer node_i(n_nodes_in), node_j(n_nodes_in)
        integer node_num(ni,nj)

        i_node = 0
        do j = 1, nj
          do i = 1, ni
            i_node = i_node + 1
            node_i(i_node) = i
            node_j(i_node) = j
            node_num(i,j) = i_node
          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine list_elems (elem_nodes_in, node_num_in, &
     &                       ni, nj, n_elems_in)
        implicit none
        integer ni, nj, n_elems_in, i, j, i_elem
        integer node_num_in(ni,nj), elem_nodes_in(n_elems_in,3)

        i_elem = 0
        do j = 1, nj-1
          do i = 1, ni-1

! define the first element in this grid box
            i_elem = i_elem + 1
            elem_nodes_in(i_elem,1) = node_num_in(i,j)
            elem_nodes_in(i_elem,2) = node_num_in(i+1,j+1)
            elem_nodes_in(i_elem,3) = node_num_in(i,j+1)

! define the second element in this grid box
            i_elem = i_elem + 1
            elem_nodes_in(i_elem,1) = node_num_in(i,j)
            elem_nodes_in(i_elem,2) = node_num_in(i+1,j)
            elem_nodes_in(i_elem,3) = node_num_in(i+1,j+1)

          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------

! note: this subroutine has been modified from the original version
! to account for full-sized data arrays (and still reduced size integer
! arrays)

      subroutine interp_data (data_out, data_in, weight, &
     &                        i_elem_ae_min, elem_nodes_in, &
     &                        node_i, node_j, &
     &                        n_nodes_in, n_nodes_out, n_elems_in, &
     &                        max_ni, max_nj, max_nodes_out)
        implicit none
        integer n_nodes_in, n_nodes_out, n_elems_in
        integer max_ni, max_nj, max_nodes_out
        integer i_elem_ae_min(n_nodes_out)
        integer elem_nodes_in(n_elems_in,3)
        integer node_i(n_nodes_in), node_j(n_nodes_in)
        real*8 data_in(max_ni,max_nj), data_out(max_nodes_out)
        real*8 weight(max_nodes_out,3)
        integer i_node, i_elem, i1, j1, i2, j2, i3, j3

! loop over the output nodes
        do i_node = 1, n_nodes_out

! get the locations of the nodes for the surrounding element on the
! input grid
          i_elem = i_elem_ae_min(i_node)
          i1 = node_i(elem_nodes_in(i_elem,1))
          j1 = node_j(elem_nodes_in(i_elem,1))
          i2 = node_i(elem_nodes_in(i_elem,2))
          j2 = node_j(elem_nodes_in(i_elem,2))
          i3 = node_i(elem_nodes_in(i_elem,3))
          j3 = node_j(elem_nodes_in(i_elem,3))

! the data on the output grid is simply the weighted data from the input
! grid
          data_out(i_node) = data_in(i1,j1) * weight(i_node,1) &
     &                     + data_in(i2,j2) * weight(i_node,2) &
     &                     + data_in(i3,j3) * weight(i_node,3)

        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine get_albedo (albedo, num_nodes, max_nodes)
        implicit none
        integer max_nodes, num_nodes, i_node
        real*8 albedo(max_nodes)

        do i_node = 1, num_nodes
          albedo(i_node) = 0.06
        enddo

      return
      end
!-----------------------------------------------------------------------
      real*8 function psi_m(zeta)
        implicit none
        real*8 zeta, chi, half_pi

        half_pi = 2.0 * atan(1.0)
        chi = (1.0 - 16.0 * zeta)**0.25
        psi_m = 2.0 * log( 0.5 * (1.0 + chi) ) + &
     &          log( 0.5 * (1.0 + chi*chi) ) - &
     &          2.0 * atan(chi) + half_pi

      return
      end
!-----------------------------------------------------------------------
      real*8 function psi_h(zeta)
        implicit none
        real*8 zeta, chi

        chi = (1.0 - 16.0 * zeta)**0.25
        psi_h = 2.0 * log( 0.5 * (1.0 + chi*chi) )

      return
      end
!-----------------------------------------------------------------------
      subroutine read_2d_arr(in_file, data, data_label, t, ni, nj)
        implicit none
        integer ni, nj, edges(2)
        real*4 data(ni,nj), t
        integer ret, read_only, sfstart, sfend, sd_id, sds_id
        integer start(2), stride(2), sds_index, sfn2index
        integer sfrdata, sfendacc, sfselect
        parameter (read_only = 1)
        character data_label*20, in_file*50, dat_time_label*33

! open in_file in read only mode
        sd_id = sfstart(in_file, read_only)

! create the data-time label, which we'll use as the name
        call get_label (dat_time_label, data_label, t)

! find index for this data set
        if (t .lt. 0.0) then
          sds_index = sfn2index(sd_id, data_label)
        else
          sds_index = sfn2index(sd_id, dat_time_label)
        endif

! find the id for this index
        sds_id = sfselect(sd_id, sds_index)

! set up start, stride, and edges to read in the entire dataset
        start(1) = 0
        start(2) = 0
        stride(1) = 1
        stride(2) = 1
        edges(1) = ni
        edges(2) = nj

! read in the dataset
        ret = sfrdata(sds_id, start, stride, edges, data)
        call checkret(ret)

! close access to the dataset
        ret = sfendacc(sds_id)
        call checkret(ret)

! close access to the file
        ret = sfend(sd_id)
        call checkret(ret)

      return
      end
!-----------------------------------------------------------------------
      subroutine read_vec_int(in_file, data, data_label, t, ni)
        implicit none
        integer ni
        real*4 t
        integer data(ni)
        integer ret, read_only, sfstart, sfend, sd_id, sds_id
        integer sds_index, sfn2index, sfrdata, sfendacc, sfselect
        parameter (read_only = 1)
        character data_label*20, in_file*50, dat_time_label*33
        
! open in_file in read only mode
        sd_id = sfstart(in_file, read_only)

! create the data-time label, which we'll use as the name
        call get_label (dat_time_label, data_label, t)
        
! find index for this data set
        if (t .lt. 0.0) then
          sds_index = sfn2index(sd_id, data_label)
        else
          sds_index = sfn2index(sd_id, dat_time_label)
        endif

! find the id for this index
        sds_id = sfselect(sd_id, sds_index)

! read in the dataset
        ret = sfrdata(sds_id, 0, 1, ni, data)
        call checkret(ret)

! close access to the dataset
        ret = sfendacc(sds_id)
        call checkret(ret)

! close access to the file
        ret = sfend(sd_id)
        call checkret(ret)

      return
      end
!-----------------------------------------------------------------------
      subroutine read_scalar(in_file, data, data_label, t)
        implicit none
        real*4 data, t
        integer ret, read_only, sfstart, sfend, sd_id, sds_id
        integer sds_index, sfn2index, sfrdata, sfendacc, sfselect
        parameter (read_only = 1)
        character data_label*20, in_file*50, dat_time_label*33
        
! open in_file in read only mode
        sd_id = sfstart(in_file, read_only)

! create the data-time label, which we'll use as the name
        call get_label (dat_time_label, data_label, t)
        
! find index for this data set
        if (t .lt. 0.0) then
          sds_index = sfn2index(sd_id, data_label)
        else
          sds_index = sfn2index(sd_id, dat_time_label)
        endif

! find the id for this index
        sds_id = sfselect(sd_id, sds_index)

! read in the dataset
        ret = sfrdata(sds_id, 0, 1, 1, data)
        call checkret(ret)

! close access to the dataset
        ret = sfendacc(sds_id)
        call checkret(ret)

! close access to the file
        ret = sfend(sd_id)
        call checkret(ret)

      return
      end
!-----------------------------------------------------------------------
      subroutine get_label (dat_time_label, data_label, t)
        implicit none
        character data_label*20, dat_time_label*33, time_label*12
        character zero_short*12
        real*4 t
        integer i
        logical nonzero
        parameter (zero_short = '      0.    ')

! create the time_label
        write(time_label,10) t
10      format(g12.6)

! since different platforms use different formats for zero, we'll use a
! standard form
        nonzero = .false.
        do i = 1, 12
          if (time_label(i:i) .ge. '1' .and. &
     &        time_label(i:i) .le. '9') nonzero = .true.
        enddo
        
        if (.not. nonzero) time_label = zero_short

! create the data-time label
        write(dat_time_label,20) data_label, ' ', time_label
20      format(a20,a1,a12)

      return
      end
!-----------------------------------------------------------------------
      subroutine checkret (ret)
        implicit none
        integer ret
        
        if (ret .ne. 0) then
          write(*,*) 'nonzero HDF return code!'
          write(*,*) 'ret = ', ret
          write(*,*)
          write(11,*) 'nonzero HDF return code!'
          write(11,*) 'ret = ', ret
          write(11,*)
          stop
!       else
!         write(*,*) 'HDF file access ok. . .'
        endif
        
      return
      end
!-----------------------------------------------------------------------
      subroutine get_weight (x_in, y_in, x_out, y_out, &
     &                       elem_nodes_in, node_i, node_j, &
     &                       max_ni, max_nj, &
     &                       n_elems_in, n_nodes_in, &
     &                       n_nodes_out, &
     &                       max_nodes_out, &
     &                       i_elem_ae_min, &
     &                       area_in, weight)
        implicit none
        integer max_ni, max_nj, n_elems_in, n_nodes_in
        integer n_nodes_out, max_nodes_out
        integer node_i(n_nodes_in), node_j(n_nodes_in)
        integer elem_nodes_in(n_elems_in,3)
        real*8 x_in(max_ni,max_nj), y_in(max_ni,max_nj)
        real*8 x_out(n_nodes_out), y_out(n_nodes_out)
        real*8 area_in(n_elems_in)
        real*8 weight(max_nodes_out,3)
        integer i_elem, i_node, i_elem_ae_min(n_nodes_out)
        integer i1, j1, i2, j2, i3, j3
        real*8 x1, y1, x2, y2, x3, y3, x4, y4, a1, a2, a3, aa, ae
        real*8 ae_min

! calculate and store the areas of the input grid elements
        do i_elem = 1, n_elems_in

          i1 = node_i(elem_nodes_in(i_elem,1))
          j1 = node_j(elem_nodes_in(i_elem,1))
          x1 = x_in(i1,j1)
          y1 = y_in(i1,j1)

          i2 = node_i(elem_nodes_in(i_elem,2))
          j2 = node_j(elem_nodes_in(i_elem,2))
          x2 = x_in(i2,j2)
          y2 = y_in(i2,j2)

          i3 = node_i(elem_nodes_in(i_elem,3))
          j3 = node_j(elem_nodes_in(i_elem,3))
          x3 = x_in(i3,j3)
          y3 = y_in(i3,j3)

          area_in(i_elem) = 0.5 * &
     &                      ( (x1-x3)*(y2-y3) + (x3-x2)*(y1-y3) )

        enddo

! now loop over the nodes of the output grid, searching for the
! surrounding elements on the input grid
        do i_node = 1, n_nodes_out

          ae_min = 1.0e25
          i_elem_ae_min(i_node) = 0

          do i_elem = 1, n_elems_in

! get the locations of the nodes for this element on the input grid
            i1 = node_i(elem_nodes_in(i_elem,1))
            j1 = node_j(elem_nodes_in(i_elem,1))
            x1 = x_in(i1,j1)
            y1 = y_in(i1,j1)

            i2 = node_i(elem_nodes_in(i_elem,2))
            j2 = node_j(elem_nodes_in(i_elem,2))
            x2 = x_in(i2,j2)
            y2 = y_in(i2,j2)

            i3 = node_i(elem_nodes_in(i_elem,3))
            j3 = node_j(elem_nodes_in(i_elem,3))
            x3 = x_in(i3,j3)
            y3 = y_in(i3,j3)

! get the locations of this node for the output grid
            x4 = x_out(i_node)
            y4 = y_out(i_node)

! calculate some type of areas
            a1 = (x4-x3)*(y2-y3) + (x2-x3)*(y3-y4)
            a2 = (x4-x1)*(y3-y1) - (y4-y1)*(x3-x1)
            a3 = (y4-y1)*(x2-x1) - (x4-x1)*(y2-y1)
            aa = abs(a1) + abs(a2) + abs(a3)
            if (area_in(i_elem) .gt. 0.0) then
              ae = abs(aa - 2.0*area_in(i_elem)) &
     &           / (2.0*area_in(i_elem))
            else
              ae = 1.0e25
            endif


! find the smallest ae
! (element of input grid that contains node for output grid)
            if (ae .lt. ae_min) then
              ae_min = ae
              i_elem_ae_min(i_node) = i_elem
            endif

          enddo

! if we didnt find a good ae_min, then there are problems
          if (i_elem_ae_min(i_node) .eq. 0) then
            write(*,*)
            write(*,*) 'Couldn''t find suitable element in input grid'
            write(*,*) 'for output node #', i_node
            write(11,*)
            write(11,*) 'Couldn''t find suitable element in input grid'
            write(11,*) 'for output node #', i_node
            stop
          endif

! the proper element has been found, now calculate the averaging weights

! get the locations of the nodes for this element on the input grid
          i_elem = i_elem_ae_min(i_node)
          i1 = node_i(elem_nodes_in(i_elem,1))
          j1 = node_j(elem_nodes_in(i_elem,1))
          x1 = x_in(i1,j1)
          y1 = y_in(i1,j1)

          i2 = node_i(elem_nodes_in(i_elem,2))
          j2 = node_j(elem_nodes_in(i_elem,2))
          x2 = x_in(i2,j2)
          y2 = y_in(i2,j2)

          i3 = node_i(elem_nodes_in(i_elem,3))
          j3 = node_j(elem_nodes_in(i_elem,3))
          x3 = x_in(i3,j3)
          y3 = y_in(i3,j3)

! get the locations of this node for the output grid
          x4 = x_out(i_node)
          y4 = y_out(i_node)

! now calculate the weighting functions
          weight(i_node,1) = ( (x4-x3)*(y2-y3) + (x2-x3)*(y3-y4) ) &
     &                     / ( 2.0*area_in(i_elem) )
          weight(i_node,2) = ( (x4-x1)*(y3-y1) - (y4-y1)*(x3-x1) ) &
     &                     / ( 2.0*area_in(i_elem) )
          weight(i_node,3) = ( -(x4-x1)*(y2-y1) + (y4-y1)*(x2-x1) ) &
     &                     / ( 2.0*area_in(i_elem) )

!999       format(2i8,g14.5,4f8.4)
!         if (mod(i_node-1,1000) .eq. 0)
!    +      write(*,999) i_node, i_elem_ae_min(i_node), ae_min,
!    +           weight(i_node,1), weight(i_node,2), weight(i_node,3),
!    +           weight(i_node,1) + weight(i_node,2) + weight(i_node,3)

        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine get_set (time, data_label, &
     &                    in_set, time_exst, &
     &                    data_out, temp_arr_1, temp_arr_2, &
     &                    input_times, num_times, ni, nj, &
     &                    num_files, max_files, window, &
     &                    time_files)
        implicit none
        integer ni, nj, num_times, num_files, max_files, time_num
        real*4 time, data_out(ni,nj), input_times(num_times)
        real*4 temp_arr_1(ni,nj), temp_arr_2(ni,nj), window
        character data_label*20, in_set*50
        character time_files(num_times)*50
        logical time_exst, exact_time_exst
        integer first_bracket
        logical bracket

! check to see if data exists at the desired time (and the file it's in)
        call does_exact_time_exst (time, input_times, num_times, &
     &                             time_num, exact_time_exst)

! if data exists at the exact time, read it in directly
        if (exact_time_exst) then
          call read_2d_arr(time_files(time_num), data_out, data_label, &
     &                     time, ni, nj)
          time_exst = .true.
          window = 0.0

! or, if data doesn't exist at this exact time
        else

! check to see if/when times exist which bracket the desired time
          call get_bracket (time, input_times, first_bracket, &
     &                      bracket, num_times)

! if bracketing times exist, then read in the data at the bracketing
! times and interpolate to the desired time
          if (bracket) then

            call interp_set(in_set, data_out, data_label, time, &
     &                      temp_arr_1, temp_arr_2, input_times, &
     &                      first_bracket, num_times, ni, nj, &
     &                      time_files, num_files, max_files)

! calculate the width of the bracketing window (convert to hours)
            window = ( input_times(first_bracket+1) - &
     &                 input_times(first_bracket) ) * 24.0

            time_exst = .true.
          else
            time_exst = .false.
          endif

        endif

      return
      end
!-----------------------------------------------------------------------
      subroutine combine_set (time, data_label, &
     &                        in_set_1, in_set_2, &
     &                        data_1, data_2, &
     &                        temp_arr_1, temp_arr_2, temp_arr_3, &
     &                        input_times_1, num_times_1, ni_1, nj_1, &
     &                        input_times_2, num_times_2, ni_2, nj_2, &
     &                        num_files_1, num_files_2, max_files, &
     &                        data_1_node, data_2_node, data_combo_node, &
     &                        in_elem_to_out_node_1, weight_node_1, &
     &                        in_elem_to_out_node_2, weight_node_2, &
     &                        elem_nodes_in_1, node_i_1, node_j_1, &
     &                        elem_nodes_in_2, node_i_2, node_j_2, &
     &                        max_ni, max_nj, max_nodes_in, &
     &                        num_nodes_out, max_nodes_out, &
     &                        num_nodes_in_1, num_elems_in_1, &
     &                        num_nodes_in_2, num_elems_in_2, &
     &                        max_elems_in, have_wind_2, &
     &                        relative_weight_1, relative_weight_2, &
     &                        max_window_1, max_window_2, &
     &                        time_files_1, time_files_2)
        implicit none
        integer max_ni, max_nj, max_elems_in, i_node
        integer num_nodes_out, max_nodes_out, max_nodes_in
        integer num_nodes_in_1, num_elems_in_1
        integer num_nodes_in_2, num_elems_in_2
        integer num_files_1, num_files_2, max_files
        integer ni_1, nj_1, num_times_1
        integer ni_2, nj_2, num_times_2
        real*4 input_times_1(num_times_1), input_times_2(num_times_2)
        real*4 data_1(ni_1, nj_1), data_2(ni_2, nj_2), time
        real*4 temp_arr_1(max_ni, max_nj), temp_arr_2(max_ni, max_nj)
        real*4 window_1, window_2
        real*8 max_window_1, max_window_2
        real*8 temp_arr_3(max_ni, max_nj)
        real*8 data_1_node(max_nodes_out), data_2_node(max_nodes_out)
        real*8 data_combo_node(max_nodes_out)
        real*8 weight_node_1(max_nodes_out,3)
        real*8 weight_node_2(max_nodes_out,3)
        real*8 relative_weight_1, relative_weight_2
        real*8 local_weight_2, sum_weights
        integer in_elem_to_out_node_1(max_nodes_out)
        integer in_elem_to_out_node_2(max_nodes_out)
        integer elem_nodes_in_1(max_elems_in,3)
        integer elem_nodes_in_2(max_elems_in,3)
        integer node_i_1(max_nodes_in), node_j_1(max_nodes_in)
        integer node_i_2(max_nodes_in), node_j_2(max_nodes_in)
        character data_label*20, in_set_1*50, in_set_2*50
        character time_files_1(num_times_1)*50
        character time_files_2(num_times_2)*50
        logical have_wind_2, time_exst_1, time_exst_2, bad_node_2
        logical use_wind_2

! wind_set_1
        call get_set (time, data_label, &
     &                in_set_1, time_exst_1, &
     &                data_1, temp_arr_1, temp_arr_2, &
     &                input_times_1, num_times_1, ni_1, nj_1, &
     &                num_files_1, max_files, window_1, &
     &                time_files_1)

! if we don't have data at this time for set 1, then bomb out with error
! message (or if window width is greater than maximum)
        if ( (.not. time_exst_1) .or. &
     &       (window_1 .gt. max_window_1) ) then
          write(*,*)
          write(*,*) 'Wind forcing data 1: time not present in files!'
          write(11,*)
          write(11,*) 'Wind forcing data 1: time not present in files!'
          stop
        endif

! copy the data to a full size real*8 array (temp_arr_3)
        call copy_arr(data_1, ni_1, nj_1, temp_arr_3, max_ni, max_nj)

! use the weightings to interpolate from the input grid to the node
! points
        call interp_data (data_1_node, temp_arr_3, weight_node_1, &
     &                    in_elem_to_out_node_1, elem_nodes_in_1, &
     &                    node_i_1, node_j_1, &
     &                    num_nodes_in_1, num_nodes_out, num_elems_in_1, &
     &                    max_ni, max_nj, max_nodes_out)

! wind_set_2
        if (have_wind_2) then

          call get_set (time, data_label, &
     &                  in_set_2, time_exst_2, &
     &                  data_2, temp_arr_1, temp_arr_2, &
     &                  input_times_2, num_times_2, ni_2, nj_2, &
     &                  num_files_2, max_files, window_2, &
     &                  time_files_2)

        else
          time_exst_2 = .false.
        endif

! even if we have data for this time, only use it if the interpolation
! window isn't too large
        use_wind_2 = time_exst_2 .and. (window_2 .le. max_window_2)

! if we have usable data at this time for set 2, then copy to full size
! array and interpolate to node points (as above)
        if (use_wind_2) then
          call copy_arr(data_2, ni_2, nj_2, temp_arr_3, max_ni, max_nj)

! use the weightings to interpolate from the input grid to the node
! points
          call interp_data (data_2_node, temp_arr_3, weight_node_2, &
     &                    in_elem_to_out_node_2, elem_nodes_in_2, &
     &                    node_i_2, node_j_2, &
     &                    num_nodes_in_2, num_nodes_out, num_elems_in_2, &
     &                    max_ni, max_nj, max_nodes_out)

        endif

! loop over all the nodes, calculating the combined value
        do i_node = 1, num_nodes_out

! determine if this node is not enclosed within the grid of data set 2
          bad_node_2 = ( &
     &        weight_node_2(i_node,1) .lt. 0.0 .or. &
     &        weight_node_2(i_node,1) .gt. 1.0 .or. &
     &        weight_node_2(i_node,2) .lt. 0.0 .or. &
     &        weight_node_2(i_node,2) .gt. 1.0 .or. &
     &        weight_node_2(i_node,3) .lt. 0.0 .or. &
     &        weight_node_2(i_node,3) .gt. 1.0 &
     &                 )

! calculate a local weight for the second data set (depending on whether
! the data exists at this time and place)
          if (use_wind_2 .and. .not. bad_node_2) then
            local_weight_2 = relative_weight_2
          else
            local_weight_2 = 0.0
          endif

! calculate the sum of the weightings
          sum_weights = relative_weight_1 + local_weight_2

! calculate a combined data value
          data_combo_node(i_node) = &
     &      ( relative_weight_1 * data_1_node(i_node) + &
     &        local_weight_2    * data_2_node(i_node) ) / sum_weights

        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine get_times(times, set_name, seek_label, &
     &                     time_files, num_times, &
     &                     num_files, max_times, max_files)
        implicit none
        integer max_times, num_times, file_number, max_files, num_files
        real*4 times(max_times), times_tmp
        character set_name*50, in_file*50, set_number*3, suffix*3
        character sds_name*33, data_label*20, time_label*12
        character seek_label*20, last_label*12, time_files(max_times)*50
        character test_time_label*12, test_dat_time_label*33
        logical exst, repeat, at_end
        integer sd_id_in, sfstart, read_only, ret, sffinfo, rank
        integer n_datasets, n_file_attrs, sds_index, sfginfo, max_rank
        integer data_type, n_attrs, sds_id, sfendacc, sfend, sfselect
        parameter (max_rank = 3)
        parameter (read_only = 1)
        integer dim_sizes(max_rank), i_time

        suffix = 'hdf'
        num_times = 0
        num_files = 0

        write(*,*)
        write(*,*) 'building database of times/files for: ', set_name
! loop over the possible file names
        file_number = 0
100     continue
          file_number = file_number + 1

110       format(a2,i1)
120       format(a1,i2)
130       format(i3)
          if (file_number .lt. 10) then
            write(set_number,110) '00', file_number
          else if (file_number .lt. 100) then
            write(set_number,120) '0', file_number
          else if (file_number .lt. 1000) then
            write(set_number,130) file_number
          else if (file_number .gt. max_files) then
            write(*,*)
            write(*,*) 'too many input files!'
            write(11,*)
            write(11,*) 'too many input files!'
            stop
          endif

! construct the file name
200       format(a15,a1,a3,a1,a3)
          write(in_file,200) set_name(1:15), '.', set_number, '.', &
     &                       suffix
!         write(*,*) in_file

! check to see if this file exists
          call file_exst (in_file, exst, .false.)
          if (exst) then

            num_files = num_files + 1
!           write(*,*) in_file,  ' exists'

! open file in read-only mode
            sd_id_in = sfstart(in_file, read_only)

! get info on input file
            ret = sffinfo(sd_id_in, n_datasets, n_file_attrs)
            call checkret(ret)

! step through the datasets
            do sds_index = 0, n_datasets-1

! get info on dataset
              sds_id = sfselect(sd_id_in, sds_index)
              ret = sfginfo(sds_id, sds_name, rank, dim_sizes, &
     &                      data_type, n_attrs)
              call checkret(ret)
!             write(*,*) sds_index, ': ', sds_name

! if this is the field we seek, then check the time label
              data_label = sds_name(1:20)
              if (data_label .eq. seek_label) then
                time_label = sds_name(22:33)
                read(time_label,*) times_tmp
!               write(*,*) 'time = ', times_tmp

! if this is the first time read, then store it and file name
                if (num_times .eq. 0) then

                  num_times = num_times + 1
                  times(num_times) = times_tmp
                  time_files(num_times) = in_file

! if this is not the first time read, then we need to make sure it's
! not a repeat
                else

! loop over previous times, comparing present time label with their
! time labels
                  repeat = .false.
                  i_time = 0
400               continue
                    i_time = i_time + 1

                    call get_label (test_dat_time_label, seek_label, &
     &                              times(i_time))

                    test_time_label = test_dat_time_label(22:33)

                    repeat = repeat .or. &
     &                       (test_time_label .eq. time_label)

                  if ( (i_time .lt. num_times) .and. (.not. repeat) ) &
     &            goto 400

!                 write(*,*)
!                 write(*,*) 'testing time = ', time_label
!                 write(*,*) 'label is a repeat: ', repeat

! if time is _not_ a repeat, then check to see if it would be at the
! end of the list of times
                  if (.not. repeat) then
                    at_end = (times_tmp .gt. times(num_times))
                  endif

! if the time label is _not_ a repeat, _and_ it would be at end of list,
! then add it and file name to the list
                  if ( (.not. repeat) .and. at_end ) then
                    num_times = num_times + 1
                    if (num_times .gt. max_times) then
                      write(*,*)
                      write(*,*) 'maximum number of times exceeded!'
                      write(11,*)
                      write(11,*) 'maximum number of times exceeded!'
                      stop
                    endif
                    times(num_times) = times_tmp
                    time_files(num_times) = in_file
!                   write(*,*) 'num_times, time: ',
!    +                         num_times, times(num_times)
                  else if (repeat) then
! if it is a repeat, then use this file for this time
                    time_files(i_time) = in_file
! endif for the not repeat block
                  endif

! endif for the first time block
                endif

! endif for the seek_label block
              endif

! terminate access to this dataset
              ret = sfendacc(sds_id)
              call checkret(ret)

            enddo

! close HDF input file
            ret = sfend(sd_id_in)
            call checkret(ret)

! endif for the file exst block
          endif

! end of loop over possible file names
        if (exst) goto 100

!     write(40,*)
!     write(40,*) 'data available at:'
!25    format ( 6(g12.6, 1x) )
!     write(40,25) (times(i_time), i_time = 1, num_times)
!     write(40,*)

      return
      end
!-----------------------------------------------------------------------
      subroutine does_exact_time_exst (time, times, num_times, &
     &                                 time_num, exact_time_exst)
        implicit none
        integer num_times, time_num, i
        real*4 times(num_times), time
        logical exact_time_exst, nonzero
        character time_label*12, test_label*12
        character zero_short*12
        parameter (zero_short = '      0.    ')

! create the time_label
        write(time_label,10) time
10      format(g12.6)

! use a standard form for a zero label
        nonzero = .false.
        do i = 1, 12
          if (time_label(i:i) .ge. '1' .and. &
     &        time_label(i:i) .le. '9') nonzero = .true.
        enddo
        if (.not. nonzero) time_label = zero_short

!       write(*,*)
!       write(*,*) 'seeking for time_label: ', time_label

! loop over input times, creating their labels, and comparing
        time_num = 0
        exact_time_exst = .false.
100     continue
          time_num = time_num + 1

! create the time_label
          write(test_label,10) times(time_num)

! once again, use a standard form for a zero label
          nonzero = .false.
          do i = 1, 12
            if (time_label(i:i) .ge. '1' .and. &
     &          time_label(i:i) .le. '9') nonzero = .true.
          enddo
          if (.not. nonzero) time_label = zero_short

! compare labels
          exact_time_exst = (time_label .eq. test_label)
!         write(*,*) 'test_label = ', test_label
!         write(*,*) 'exact_time_exst = ', exact_time_exst

! bottom of loop
          if ( (.not. exact_time_exst) .and. (time_num .lt. num_times) ) &
     &    goto 100
        
      return
      end
!-----------------------------------------------------------------------
