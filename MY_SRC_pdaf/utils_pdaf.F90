!> Utility Routines
!!
!! This module contains several routines useful for common
!! model tasks. The initial routines included output configuration
!! information about the PDAF library, and configuration information
!! about the assimilation parameters.
!! 
module utils_pdaf

  use mod_kind_pdaf

  implicit none
  save

contains

  !> This routine performs a model-sided screen output about
  !! the configuration of the data assimilation system.
  !!
  !! **Calling Sequence**
  !!
  !! - Called from: `init_pdaf`
  !!
  subroutine init_info_pdaf()

    use assimilation_pdaf, & ! Variables for assimilation
         only: filtertype, subtype, dim_ens, delt_obs, model_error, &
         model_err_amp, forget, rank_analysis_enkf

    ! *****************************
    ! *** Initial Screen output ***
    ! *****************************

    if (filtertype == 1) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: SEIK'
       if (subtype == 2) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed error-space basis'
       else if (subtype == 3) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed state covariance matrix'
       else if (subtype == 4) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- use ensemble transformation'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 2) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: EnKF'
       if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
       if (rank_analysis_enkf > 0) then
          write (*, '(a,6x, a, i5)') 'NEMO-PDAF', &
               'analysis with pseudo-inverse of HPH, rank:', rank_analysis_enkf
       end if
    else if (filtertype == 3) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LSEIK'
       if (subtype == 2) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed error-space basis'
       else if (subtype == 3) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed state covariance matrix'
       else if (subtype == 4) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- use ensemble transformation'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 4) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: ETKF'
       if (subtype == 0) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant using T-matrix'
       else if (subtype == 1) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant following Hunt et al. (2007)'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 5) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LETKF'
       if (subtype == 0) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant using T-matrix'
       else if (subtype == 1) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant following Hunt et al. (2007)'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 6) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: ESTKF'
       if (subtype == 0) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Standard mode'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    else if (filtertype == 7) then
       write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LESTKF'
       if (subtype == 0) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Standard mode'
       else if (subtype == 5) then
          write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
       end if
       write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
       if (subtype /= 5) write (*, '(a,6x, a, i5)') 'NEMO-PDAF', 'Assimilation interval:', delt_obs
       write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
       if (model_error) then
          write (*, '(a,6x, a, f5.2)') 'NEMO-PDAF', 'model error amplitude:', model_err_amp
       end if
    end if

  end subroutine init_info_pdaf

!-------------------------------------------------------------------------------

  !> This routine reads the namelist file with parameters
  !! controlling data assimilation with PDAF and outputs to
  !! screen.
  !!
  !! **Calling Sequence**
  !!
  !! - Called from: `init_pdaf`
  !!
  subroutine read_config_pdaf()

    use parallel_pdaf, &
         only: mype_ens
    use assimilation_pdaf, &
         only: filtertype, subtype, type_trans, type_sqrt, &
         locweight, screen, dim_ens, ensscale, delt_obs, &
         type_forget, forget, type_ens_init, type_central_state, &
         ens_restart, type_hyb, hyb_gamma, hyb_kappa
    use io_pdaf, &
         only: verbose_io, path_inistate, path_ens, file_ens, file_covar, &
         sgldbl_io, coupling_nemo, save_var, save_state, save_ens_sngl, &
         ids_write, add_slash
    use statevector_pdaf, &
         only: update_ssh, update_temp, update_salt, update_vel
#if defined key_top
    use statevector_pdaf, &
         only: update_trc
#endif
    use asminc_pdaf, &
         only: do_asmiau, steps_asmiau, shape_asmiau, iter_divdmp
    use obs_ssh_mgrid_pdafomi, &
         only: assim_ssh_mgrid, rms_ssh_mgrid, file_ssh_mgrid, &
         lradius_ssh_mgrid, sradius_ssh_mgrid, varname_ssh_mgrid
    use obs_sst_cmems_pdafomi, &
         only: assim_sst_cmems, path_sst_cmems, file_sst_cmems, rms_obs_sst_cmems, &
         lradius_sst_cmems, sradius_sst_cmems, mode_sst_cmems, dist_sst_cmems, &
         varname_sst_cmems
#if defined key_top
    use asminc_pdaf, &
         only: do_bgciau, steps_bgciau, shape_asmiau, shape_bgciau
#endif

    !< Namelist file
    character(lc) :: nmlfile

    namelist /pdaf_nml/ &
         screen, filtertype, subtype, type_trans, type_sqrt, &
         type_forget, forget, locweight, delt_obs, save_var, &
         save_state, save_ens_sngl, ids_write, verbose_io, sgldbl_io, &
         type_hyb, hyb_gamma, hyb_kappa

    namelist /init_nml/ &
         type_ens_init, type_central_state, ensscale, ens_restart, &
         path_inistate, path_ens, file_ens, file_covar, coupling_nemo

#if defined key_top
    namelist /update_nml/ &
         update_ssh, update_temp, update_salt, update_vel, &
         do_asmiau, steps_asmiau, shape_asmiau, iter_divdmp, &
         update_trc, do_bgciau, steps_bgciau, shape_bgciau
#else
    namelist /update_nml/ &
         update_ssh, update_temp, update_salt, update_vel, &
         do_asmiau, steps_asmiau, shape_asmiau, iter_divdmp
#endif

    namelist /obs_ssh_mgrid_nml/ &
         assim_ssh_mgrid, rms_ssh_mgrid, file_ssh_mgrid, &
         lradius_ssh_mgrid,  sradius_ssh_mgrid, varname_ssh_mgrid

    namelist /obs_sst_cmems_nml/ &
         assim_sst_cmems, path_sst_cmems, file_sst_cmems, rms_obs_sst_cmems, &
         lradius_sst_cmems, sradius_sst_cmems, mode_sst_cmems, dist_sst_cmems, &
         varname_sst_cmems


    ! ****************************************************
    ! ***   Initialize PDAF parameters from namelist   ***
    ! ****************************************************

    ! Initialize array of single field IDs for which ensemble could be written
    ids_write = 0

    nmlfile = 'namelist_cfg.pdaf'

    open (20, file=nmlfile)
    read (20, NML=pdaf_nml)
    rewind(20)
    read (20, NML=init_nml)
    rewind(20)
    read (20, NML=update_nml)
    rewind(20)
    read (20, NML=obs_ssh_mgrid_nml)
    rewind(20)
    read (20, NML=obs_sst_cmems_nml)
    rewind(20)
    close (20)

    ! *** Add trailing slash to paths ***
    call add_slash(path_sst_cmems)
    call add_slash(path_inistate)
    call add_slash(path_ens)

    ! *** Set flags for ensemble restart ***
    if (ens_restart) then
       type_ens_init = 4
       type_central_state = 0
    end if

    ! Print PDAF parameters to screen
    showconf: if (mype_ens == 0) then

       write (*, '(/a,1x,a)') 'NEMO-PDAF','-- Overview of PDAF configuration --'
       write (*, '(a,3x,a)') 'NEMO-PDAF','[pdaf_nml]:'
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','screen       ', screen
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','filtertype   ', filtertype
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','subtype      ', subtype
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_trans   ', type_trans
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_sqrt    ', type_sqrt
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_forget  ', type_forget
       write (*, '(a,5x,a,f10.3)') 'NEMO-PDAF','forget       ', forget
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','locweight    ', locweight
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','delt_obs     ', delt_obs
       write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','save_var     ', trim(save_var)
       write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','save_state   ', trim(save_state)
       write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','sgldbl_io    ', sgldbl_io
       if (filtertype==11) then
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_hyb       ', type_hyb
       write (*, '(a,5x,a,f10.3)') 'NEMO-PDAF','hyb_gamma      ', hyb_gamma
       write (*, '(a,5x,a,f10.3)') 'NEMO-PDAF','hyb_kappa      ', hyb_kappa
       end if
       write (*, *) ''
       write (*, '(a,3x,a)') 'NEMO-PDAF','[init_nml]:'
       write (*, '(a,5x,a,l)') 'NEMO-PDAF','ens_restart ', ens_restart
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_ens_init      ', type_ens_init
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_central_state ', type_central_state
       write (*, '(a,5x,a,5x,a)') 'NEMO-PDAF','coupling_nemo        ', coupling_nemo
       write (*, '(a,5x,a,5x,f10.2)') 'NEMO-PDAF','ensscale     ', ensscale
       write (*, *) ''
       write (*, '(a,3x,a)') 'NEMO-PDAF','[update_nml]:'
       write (*, '(a,5x,a,l)') 'NEMO-PDAF','update_ssh      ', update_ssh
       write (*, '(a,5x,a,l)') 'NEMO-PDAF','update_temp     ', update_temp
       write (*, '(a,5x,a,l)') 'NEMO-PDAF','update_salt     ', update_salt
       write (*, '(a,5x,a,l)') 'NEMO-PDAF','update_vel      ', update_vel
       write (*, *) ''
       write (*, '(a,3x,a)') 'NEMO-PDAF','[obs_ssh_mgrid_nml]:'
       write (*, '(a,5x,a,5x,l)') 'NEMO-PDAF','assim_ssh_mgrid      ', assim_ssh_mgrid
       if (assim_ssh_mgrid) then
          write (*, '(a,5x,a,f12.4)') 'NEMO-PDAF','rms_ssh_mgrid      ', rms_ssh_mgrid
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','file_ssh_mgrid        ', trim(file_ssh_mgrid)
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','varname_ssh_mgrid     ', trim(varname_ssh_mgrid)
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','lradius_ssh_mgrid      ', lradius_ssh_mgrid
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','sradius_ssh_mgrid      ', sradius_ssh_mgrid
       end if
       write (*, *) ''
       write (*, '(a,3x,a)') 'NEMO-PDAF','[obs_sst_cmems_nml]:'
       write (*, '(a,5x,a,5x,l)') 'NEMO-PDAF','assim_sst_cmems      ', assim_sst_cmems
       if (assim_sst_cmems) then
          write (*, '(a,5x,a,f12.4)') 'NEMO-PDAF','rms_obs_sst_cmems  ', rms_obs_sst_cmems
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','file_sst_cmems        ', trim(file_sst_cmems)
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','path_sst_cmems        ', trim(path_sst_cmems)
          write (*, '(a,5x,a,a)') 'NEMO-PDAF','varname_sst_cmems     ', trim(varname_sst_cmems)
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','lradius_sst_cmems      ', lradius_sst_cmems
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','sradius_sst_cmems      ', sradius_sst_cmems
          write (*, '(a,5x,a,i10)') 'NEMO-PDAF','mode_sst_cmems      ', mode_sst_cmems
          write (*, '(a,5x,a,8x,a)') 'NEMO-PDAF','dist_sst_cmems      ', dist_sst_cmems
       end if
       write (*, '(a,1x,a/)') 'NEMO-PDAF','-- End of PDAF configuration overview --'

    end if showconf

  end subroutine read_config_pdaf

!-------------------------------------------------------------------------------

!> Timing and clean-up of PDAF
!!
!! The routine prints timing and memory information.
!! It further deallocates PDAF internal arrays
!! and the ASMINC increment array for BGC variables
!!
!! - Called from: `nemogcm`
!!
  subroutine finalize_pdaf()

    use parallel_pdaf, &
         only: mype_ens, comm_ensemble, mpierr, mype_model
    use asminc_pdaf, &
         only: asm_inc_deallocate_pdaf
    use timer, &
         only: timeit, time_tot


    ! Show allocated memory for PDAF
    if (mype_ens==0) call PDAF_print_info(10)
    call PDAF_print_info(11)

    ! Print PDAF timings onto screen
    if (mype_ens==0) call PDAF_print_info(3)

    ! Deallocate PDAF arrays
    call PDAF_deallocate()

    ! Deallocate ASMINC arrays
    call asm_inc_deallocate_pdaf()

    call timeit(4,'old')
    call timeit(5,'old')

    if (mype_ens==0) then
       write (*, '(/a,10x,a)') 'NEMO-PDAF', 'Model-sided timing overview'
       write (*, '(a,2x,a)') 'NEMO-PDAF', '-----------------------------------'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize MPI:  ', time_tot(1), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize model:', time_tot(2), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize PDAF :', time_tot(3), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'main part:       ', time_tot(4), 's'
       write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'total:         ', time_tot(5), 's'
    end if

    call mpi_barrier(comm_ensemble, mpierr)

  end subroutine finalize_pdaf

end module utils_pdaf
