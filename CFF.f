      program prob_invagination_w_affected_rate

      implicit none
      integer i,j,k
      integer maxnpart,npart,nspart
      parameter (maxnpart=8192)
      integer IC_npart,IC_nspart
      double precision IC_x(maxnpart),IC_y(maxnpart)
      double precision IC_vx(maxnpart),IC_vy(maxnpart)
      double precision IC_radius(maxnpart),IC_mass(maxnpart)
      double precision Rratio,radius1
      double precision dt,dt_sigma1
      double precision gamman,gamman_sigma1
      common /damping/ gamman,gamman_sigma1
      double precision boxx,boxy
      double precision x(maxnpart),y(maxnpart)
      double precision radius(maxnpart),ini_radius(maxnpart)
      double precision vx(maxnpart),vy(maxnpart),mass(maxnpart)
      character*60 str1,str2,str3,str8,str9,str10,str11
      double precision V
      double precision attract_fx(maxnpart),attract_fy(maxnpart)
      double precision attract_V
      double precision repel_fx(maxnpart),repel_fy(maxnpart),repel_V
      double precision attract_dampfx(maxnpart),attract_dampfy(maxnpart)
      double precision repel_dampfx(maxnpart),repel_dampfy(maxnpart)
      integer shrink_idx,step_count,tot_step_count
      integer cycle_idx
      double precision V_old,V_rel,stop_V_rel
      parameter(stop_V_rel=1.d-3)
      integer iniVlist_update
      double precision list_range
      common /iniVlist1/ iniVlist_update
      common /iniVlist2/ list_range
      double precision attract_tot_deform(maxnpart)
      double precision repel_tot_deform(maxnpart)
      double precision tot_deform(maxnpart)
      double precision max_tot_deform
      integer constrict_label(maxnpart),ini_constrict_label(maxnpart)
      integer ran2_seed
      double precision ran2
      integer seed(12)
      double precision rand
      integer num_skip_seed,IC_idx
      integer ini_num_active_part
      integer invagination_type,num_active_part
      double precision active_zone_range, constriction_rate
      integer pbc_y_switch
      common /pbc/ pbc_y_switch
      double precision rcut,rlist,ranrlist
      double precision alpha,beta,cut_tot_deform,cut,constriction_prob
      integer num_constricted_part
      double precision trial_rand
      integer num_negative_prob,num_trial
      integer normalize_switch
      double precision ave_single_deform
      integer bdry_force_switch
      double precision bdry_y_thick
      double precision top_bdry,bottom_bdry
      double precision bdryfx(maxnpart),bdryfy(maxnpart)
      integer bdry_label(maxnpart)
      common /boundary0/ bdry_force_switch
      common /boundary1/ radius1,bdry_y_thick,top_bdry,bottom_bdry
      common /boundary2/ bdry_label
      double precision spring_k,spring_k_sigma1
      common /spring/ spring_k,spring_k_sigma1
      !********************** for affected prob ************************
      ! here shrink and constrict have the same physical meaning
      double precision affected_left_bdry1,affected_right_bdry1
      double precision affected_left_bdry2,affected_right_bdry2
      double precision affected_left_bdry3,affected_right_bdry3
      double precision affected_left_bdry4,affected_right_bdry4
      double precision affected_top_bdry,affected_bottom_bdry
      double precision affected_rate1,affected_rate2
      double precision affected_rate3,affected_rate4
      integer initiation_cycle_idx
      integer affected_label(maxnpart)
      integer num_affected_part
      integer num_constricted_affected_part
      integer num_unconstricted_affected_part ! after constriction
      !********************** for affected prob ************************
      integer saved_step_count,saved_num_trial
      integer saved_num_unconstricted_affected_part
      double precision saved_V
      ! status label: ini = 0; ua = 1; ca = 2; uu = 3; cu = 4
      integer status_label(maxnpart)
      !************** for deformation (stress) data analysis *************
      ! ua = 'u'nconstricted 'a'ffected
      integer num_ua_part                     ! before constriction
      integer ua_part_idx(maxnpart)
      double precision ua_deform_data(maxnpart)
      double precision ave_ua_deform
      double precision max_ua_deform

      ! ca = 'c'onstricted 'a'ffected
      integer num_ca_part
      integer ca_part_idx(maxnpart)
      double precision ca_deform_data(maxnpart)
      double precision ave_ca_deform
      double precision max_ca_deform
      
      ! uu = 'u'nconstricted 'u'naffected
      integer num_uu_part
      integer uu_part_idx(maxnpart)
      double precision uu_deform_data(maxnpart)
      double precision ave_uu_deform
      double precision max_uu_deform

      ! cu = 'c'onstricted 'u'naffected
      integer num_cu_part
      integer cu_part_idx(maxnpart)
      double precision cu_deform_data(maxnpart)
      double precision ave_cu_deform
      double precision max_cu_deform
      !************** for deformation (stress) data analysis *************
      integer movie_switch !0: do not save; 1: save movie configs 

      boxx = 1.d0
      boxy = 1.d0

      write(*,*) 'Enter name of IC file.'
      read(*,*) str1
      write(*,*) 'Enter num of skipped ran2 seed.'
      read(*,*) num_skip_seed
      write(*,*) 'Enter IC_idx.'
      read(*,*) IC_idx
      write(*,*) 'Enter movie_switch (0: no movie; 1: save movie).'
      read(*,*) movie_switch
      if(movie_switch .eq. 1) then
        write(*,*) 'Enter name of file for configs for movies.'
        read(*,*) str2
      endif
      write(*,*) 'Enter name of file for info.'
      read(*,*) str3
      ! file 8-11 are for deformation (stress) analysis
      write(*,*) 'Enter file name for unconstricted affected deform.'
      read(*,*) str8
      write(*,*) 'Enter file name for constricted affected deform.'
      read(*,*) str9
      write(*,*) 'Enter file name for unconstricted unaffected deform.'
      read(*,*) str10
      write(*,*) 'Enter file name for constricted unaffected deform.'
      read(*,*) str11
      write(*,*) 'Enter value of time step of particle size 1.'
      read(*,*) dt_sigma1
      write(*,*) 'Enter gamman_sigma1.'
      read(*,*) gamman_sigma1
      write(*,*) 'Enter list_range (rij < sigmaij*list_range = contact)'
      write(*,*) 'for generating an initial neighbor list of IC.'
      read(*,*) list_range
      write(*,*) 'Enter invagination_type. Part w/in the active zone'
      write(*,*) 'will be constricted --- 0: randomly; or by'
      write(*,*) '1: expansion in the xx dir.;'
      write(*,*) '20: contraction only w/o expansion;'
      write(*,*) '21: overall contraction from expansion + contraction;'
      write(*,*) '3: abs(deformation).'
      read(*,*) invagination_type
      write(*,*) 'Enter 0 <= active_zone_range <= 0.5, Part init w/in'
      write(*,*) '[-active_zone_range*boxy : active_zone_range*boxy]'
      write(*,*) 'will be constricted.'
      read(*,*) active_zone_range
      write(*,*) 'Enter 0 <= constriction_rate <= 1.'
      read(*,*) constriction_rate
      write(*,*) 'Enter pbc_y_switch (0: off; 1: on).'
      read(*,*) pbc_y_switch
      write(*,*) 'Enter rlist/rcut.'
      read(*,*) ranrlist
      write(*,*) 'Enter Rratio of large over small particles.'
      read(*,*) Rratio
      write(*,*) 'Enter alpha for constriction_prob.'
      read(*,*) alpha
      write(*,*) 'Enter beta for constriction_prob.'
      read(*,*) beta
      write(*,*) 'Enter cut_tot_deform for normalized tot_deform.'
      write(*,*) '(guarantees tot_deform <= cut_tot_deform)'
      read(*,*) cut_tot_deform
      write(*,*) 'Enter cut for constriction_prob.'
      write(*,*) '(guarantees 1 + beta * tot_deform**3 >= cut)'
      read(*,*) cut
      write(*,*) 'Enter normalization switch for tot_deform'
      write(*,*) '(0: using max_tot_deform; 1: using a constant value).'
      read(*,*) normalize_switch
      if(normalize_switch .eq. 1) then
        write(*,*) 'Enter the const value for tot_deform normalization.'
        read(*,*) ave_single_deform
      endif
      write(*,*) 'Enter top_bdry (< boxy/2.d0 = 0.5d0).'
      read(*,*) top_bdry
      write(*,*) 'Enter bottom_bdry (> -boxy/2.d0 = -0.5d0).'
      read(*,*) bottom_bdry
      write(*,*) 'Enter bdry_y_thick.'
      read(*,*) bdry_y_thick
      write(*,*) 'Enter spring_k for box size 1 from stress analysis.'
      read(*,*) spring_k
      !*********************** for affected prob ***********************
      write(*,*) '--- below are for affected constriction process ---'
      write(*,*) 'part w/ini_x=[affected_left_bdry:affected_right_bdry]'
      write(*,*) '&    w/ini_y=[affected_bottom_bdry:affected_top_bdry]'
      write(*,*) 'its constriction_prob=constriction_prob*affected_rate'
      write(*,*) 'Enter affected_left_bdry1 (> -boxx/2.d0 = -0.5d0).'
      read(*,*) affected_left_bdry1
      write(*,*) 'Enter affected_right_bdry1 (< boxx/2.d0 = 0.5d0).'
      read(*,*) affected_right_bdry1
      write(*,*) 'Enter affected_left_bdry2 (> -boxx/2.d0 = -0.5d0).'
      read(*,*) affected_left_bdry2
      write(*,*) 'Enter affected_right_bdry2 (< boxx/2.d0 = 0.5d0).'
      read(*,*) affected_right_bdry2
      write(*,*) 'Enter affected_left_bdry3 (> -boxx/2.d0 = -0.5d0).'
      read(*,*) affected_left_bdry3
      write(*,*) 'Enter affected_right_bdry3 (< boxx/2.d0 = 0.5d0).'
      read(*,*) affected_right_bdry3
      write(*,*) 'Enter affected_left_bdry4 (> -boxx/2.d0 = -0.5d0).'
      read(*,*) affected_left_bdry4
      write(*,*) 'Enter affected_right_bdry4 (< boxx/2.d0 = 0.5d0).'
      read(*,*) affected_right_bdry4
      write(*,*) 'Enter affected_top_bdry'
      write(*,*) '(< active_zone_range * boxy = active_zone_range).'
      read(*,*) affected_top_bdry
      write(*,*) 'Enter affected_bottom_bdry'
      write(*,*) '(> -active_zone_range * boxy = -active_zone_range).'
      read(*,*) affected_bottom_bdry
      write(*,*) 'Enter affected_rate1 for affected region 1'
      read(*,*) affected_rate1
      write(*,*) 'Enter affected_rate2 for affected region 2'
      read(*,*) affected_rate2
      write(*,*) 'Enter affected_rate3 for affected region 3'
      read(*,*) affected_rate3
      write(*,*) 'Enter affected_rate4 for affected region 4'
      read(*,*) affected_rate4
      write(*,*) 'Enter initiation_cycle_idx for affected regions 1 & 2'
      write(*,*) '(if cycle_idx >= initiation_cycle_idx, affected part'
      write(*,*) 'in these two regions are forced to constrict.)'
      read(*,*) initiation_cycle_idx
      !*********************** for affected prob ***********************

      open(unit=1,file=str1,status='old')
      if(movie_switch .eq. 1) then
        open(unit=2,file=str2,status='unknown')
      endif
      open(unit=3,file=str3,status='unknown')
      ! file 8-11 are for deformation (stress) analysis
      open(unit=8,file=str8,status='unknown')
      open(unit=9,file=str9,status='unknown')
      open(unit=10,file=str10,status='unknown')
      open(unit=11,file=str11,status='unknown')

      if(normalize_switch .eq. 0) then
        write(3,*) '# normalization switch = 0, using max_tot_deform'
      elseif(normalize_switch .eq. 1) then
        write(3,*) '# normalization switch = 1, using a constant'
        write(3,*) '# tot_deform normalization =',ave_single_deform
      endif
      write(3,*) '# *** (all data are relaxed) ***'
      write(3,*) '# (1) cycle_idx (clock of the simulation)'
      write(3,*) '# (2) shrink_idx'
      write(3,*) '# (3) step_count (4) tot_V'
      write(3,*) '# (5) #_ran_trial (6) #_unconstricted_affected_part'
      write(3,*) '# (7) ave unconstricted affected (ua) deform'
      write(3,*) '# (8) max unconstricted affected (ua) deform'
      write(3,*) '# (9) ave constricted affected (ca) deform'
      write(3,*) '# (10) max constricted affected (ca) deform'
      write(3,*) '# (11) ave unconstricted unaffected (uu) deform'
      write(3,*) '# (12) max unconstricted unaffected (uu) deform'
      write(3,*) '# (13) ave constricted unaffected (cu) deform'
      write(3,*) '# (14) max constricted unaffected (cu) deform'

      read(1,*) IC_npart
      read(1,*) IC_nspart
      do i=1,IC_npart
        read(1,*) IC_x(i),IC_y(i),IC_vx(i),IC_vy(i),IC_radius(i)

        ! Put particle back in the main box
        IC_x(i) = IC_x(i) - boxx*dnint(IC_x(i)/boxx)
        IC_y(i) = IC_y(i) - boxy*dnint(IC_y(i)/boxy)
      enddo

      !_______________________ setting parameters _______________________
      radius1 = IC_radius(1)

      do i=1,IC_npart
        IC_mass(i) = (IC_radius(i) / radius1)**2
      enddo

      rcut = 2.d0 * radius1 * Rratio
      rlist = ranrlist * rcut

      ! the following lines translate length scale from particle size 1
      ! to box size 1. Choose smaller particle diameter as the ruler.
      dt = dt_sigma1 * 2.d0 * radius1
      gamman = gamman_sigma1 / (2.d0 * radius1)
      !__________________________________________________________________

      !________________________ trimming IC ________________________
      npart = 0
      nspart = 0

      do i=1,IC_npart
        if(IC_y(i) .le. top_bdry + bdry_y_thick * radius1 .and.
     &     IC_y(i) .ge. bottom_bdry - bdry_y_thick * radius1) then

          npart = npart + 1
          if(i .le. IC_npart/2) nspart = nspart + 1

          x(npart) = IC_x(i)
          y(npart) = IC_y(i)
          vx(npart) = IC_vx(i)
          vy(npart) = IC_vy(i)
          radius(npart) = IC_radius(i)
          mass(npart) = IC_mass(i)

        endif
      enddo

      write(*,*) 'npart after trimming =',npart
      write(*,*) 'nspart after trimming =',nspart
      !write(*,*) 'x(npart),y(npart) =',x(npart),y(npart)
      !____________________ end of trimming IC _____________________

      ! save particle info in a file
 20   format(5(e21.14,3X),2(I5,3X))

      if(movie_switch .eq. 1) then
        write(2,*) npart, -1, -1 ! shrink_idx = cycle_idx = -1 (trimmed IC)
        write(2,*) nspart
        do i=1,npart
          ! x,y,vx,vy,radius,bdry_label,affected_label
          write(2,20) x(i),y(i),vx(i),vy(i),radius(i),0,0
        enddo
        write(2,*) ''
      endif

      !_______________________ create Vlist of IC _______________________
      iniVlist_update = 1
      call iniVlist_attract_force(x,y,attract_fx,attract_fy,attract_V,
     &                            npart,boxx,boxy,radius)

      call iniVlist_attract_dampforce(x,y,vx,vy,attract_dampfx,
     &                            attract_dampfy,npart,boxx,boxy,radius)

      call iniVlist_attract_tot_deform(x,y,npart,boxx,boxy,radius,
     &                                 attract_tot_deform)
      ! turn off all iniVlist updates
      iniVlist_update = 0
      !___________________ end of create Vlist of IC ____________________

      tot_step_count = 0

      !__________________________ relaxation of IC ______________________
      V_old = 0.d0
      V_rel = 1.d0
      step_count = 0

      do while(V_rel .gt. stop_V_rel)
        call verlet(x,y,vx,vy,npart,boxx,boxy,radius,dt,mass,V,
     &              rcut,rlist)

        step_count = step_count + 1

        if(mod(step_count,1000) .eq. 0) then

          V_rel = dabs(V - V_old)/V
          write(*,*) 'V_rel =',V_rel
          V_old = V

        endif
      enddo ! do while(V_rel .gt. stop_V_rel)

      !***** assigning bdry particles after IC relaxation *****
      bdry_force_switch = 1
      call bdry_force(x,y,bdryfx,bdryfy,npart,boxx,boxy,radius)
      ! turn off assigning bdry particles
      bdry_force_switch = 0
      !********************************************************

      !if(movie_switch .eq. 1) then
      !  write(2,*) npart, 0
      !  write(2,*) nspart
      !  do i=1,npart
      !    ! x,y,vx,vy,radius
      !    write(2,20) x(i),y(i),vx(i),vy(i),radius(i),bdry_label(i)
      !  enddo
      !  write(2,*) ''
      !endif

      write(*,*) 'relaxed after',step_count,'steps.'
      write(*,*) ''
      write(3,*) '# 0 (IC relaxation)',step_count,V
      tot_step_count = tot_step_count + step_count
      !____________________ end of relaxation of IC _____________________



      !_____________________ invagination process _______________________
      ini_num_active_part = 0

      do i=1,npart ! initialization
        ini_radius(i) = radius(i) ! initial particle sizes
        status_label(i) = 0 ! ini = 0

        if(y(i) .ge. - active_zone_range * boxy .and.
     &     y(i) .le.   active_zone_range * boxy) then
          ini_num_active_part = ini_num_active_part + 1
          ini_constrict_label(i) = 0 ! 0: unconstricted
          constrict_label(i) = 0     ! 0: unconstricted
        else
          ini_constrict_label(i) = 1 ! 1: will not be constricted
          constrict_label(i) = 1     ! 1: will not be constricted
        endif
      enddo

      write(*,*) 'initial num_active_part =',ini_num_active_part
      write(3,*) '# initial num_active_part ='
      write(3,*) '#',ini_num_active_part

      !*********************** for affected prob ************************
      num_affected_part = 0

      do i=1,npart
        ! affected region 1
        if(x(i) .ge. affected_left_bdry1 .and.
     &     x(i) .le. affected_right_bdry1 .and.
     &     y(i) .ge. affected_bottom_bdry .and.
     &     y(i) .le. affected_top_bdry .and.
     &     y(i) .ge. - active_zone_range * boxy .and.
     &     y(i) .le.   active_zone_range * boxy) then
          num_affected_part = num_affected_part + 1
          affected_label(i) = 1 ! constriction_prob will be affected
        ! affected region 2
        elseif(x(i) .ge. affected_left_bdry2 .and.
     &         x(i) .le. affected_right_bdry2 .and.
     &         y(i) .ge. affected_bottom_bdry .and.
     &         y(i) .le. affected_top_bdry .and.
     &         y(i) .ge. - active_zone_range * boxy .and.
     &         y(i) .le.   active_zone_range * boxy) then
          num_affected_part = num_affected_part + 1
          affected_label(i) = 2 ! constriction_prob will be affected
        ! affected region 3
        elseif(x(i) .ge. affected_left_bdry3 .and.
     &         x(i) .le. affected_right_bdry3 .and.
     &         y(i) .ge. affected_bottom_bdry .and.
     &         y(i) .le. affected_top_bdry .and.
     &         y(i) .ge. - active_zone_range * boxy .and.
     &         y(i) .le.   active_zone_range * boxy) then
          num_affected_part = num_affected_part + 1
          affected_label(i) = 3 ! constriction_prob will be affected
        ! affected region 4
        elseif(x(i) .ge. affected_left_bdry4 .and.
     &         x(i) .le. affected_right_bdry4 .and.
     &         y(i) .ge. affected_bottom_bdry .and.
     &         y(i) .le. affected_top_bdry .and.
     &         y(i) .ge. - active_zone_range * boxy .and.
     &         y(i) .le.   active_zone_range * boxy) then
          num_affected_part = num_affected_part + 1
          affected_label(i) = 4 ! constriction_prob will be affected
        else
          affected_label(i) = 0 ! constriction_prob will not be affected
        endif
      enddo

      write(*,*) 'initial num_affected_part =',num_affected_part
      write(3,*) '# initial num_affected_part ='
      write(3,*) '#',num_affected_part
      !*********************** for affected prob ************************

      ! initialization
      shrink_idx = 0
      cycle_idx = 0
      saved_num_trial = 0

      ran2_seed = -(num_skip_seed + IC_idx) ! using large negative int
      write(3,*) '# ini ran2_seed =',ran2_seed

      ! The RANDOM_NUMBER accepts 12 integer seeds by gfortran,
      ! which is reduced to the range [1, 2147483562].
      ! This means that the generator effectively uses two 31-bit seeds.
      seed(1)  = int(ran2(ran2_seed) * 1.0d6)
      seed(2)  = int(ran2(ran2_seed) * 2.0d6)
      seed(3)  = int(ran2(ran2_seed) * 3.0d6)
      seed(4)  = int(ran2(ran2_seed) * 4.0d6)
      seed(5)  = int(ran2(ran2_seed) * 5.0d6)
      seed(6)  = int(ran2(ran2_seed) * 6.0d6)
      seed(7)  = int(ran2(ran2_seed) * 7.0d6)
      seed(8)  = int(ran2(ran2_seed) * 8.0d6)
      seed(9)  = int(ran2(ran2_seed) * 9.0d6)
      seed(10) = int(ran2(ran2_seed) * 10.0d6)
      seed(11) = int(ran2(ran2_seed) * 11.0d6)
      seed(12) = int(ran2(ran2_seed) * 12.0d6)
      write(3,*) '# random_number seed1, seed2   =',seed(1),seed(2)
      write(3,*) '# random_number seed3, seed4   =',seed(3),seed(4)
      write(3,*) '# random_number seed5, seed6   =',seed(5),seed(6)
      write(3,*) '# random_number seed7, seed8   =',seed(7),seed(8)
      write(3,*) '# random_number seed9, seed10  =',seed(9),seed(10)
      write(3,*) '# random_number seed11, seed12 =',seed(11),seed(12)
      call random_seed(put = seed)

 1    continue

        !_______________ selecting part for constriction ________________
        !________________________________________________________________
        !________________________________________________________________
        ! calculate total deformation w.r.t. Vlisted neighbors
        call iniVlist_attract_tot_deform(x,y,npart,boxx,boxy,radius,
     &                                   attract_tot_deform)
        call Vlist_repel_tot_deform(rcut,rlist,x,y,npart,
     &                              boxx,boxy,radius,repel_tot_deform)
        do i=1,npart
          if(invagination_type .eq. 1) then
            tot_deform(i) = attract_tot_deform(i)
          elseif(invagination_type .eq. 20) then
            tot_deform(i) = repel_tot_deform(i)
          elseif(invagination_type .eq. 21) then
            tot_deform(i) = attract_tot_deform(i) + repel_tot_deform(i)
          elseif(invagination_type .eq. 3) then
            tot_deform(i) = attract_tot_deform(i) + repel_tot_deform(i)
          endif
        enddo

        do i=1,npart
          if(invagination_type .eq. 0) then
            ! part will be constricted randomly
            call random_number(rand)
            tot_deform(i) = rand
          elseif(invagination_type .eq. 1) then
            ! part w/ larger expansion will be constricted first
            tot_deform(i) = tot_deform(i)
          elseif(invagination_type .eq. 20) then
            ! part w/ larger contraction will be constricted first
            tot_deform(i) = - tot_deform(i)
          elseif(invagination_type .eq. 21) then
            ! part w/ larger contraction will be constricted first
            tot_deform(i) = - tot_deform(i)
          elseif(invagination_type .eq. 3) then
            ! part w/ larger deformation will be constricted first
            tot_deform(i) = dabs(tot_deform(i))
          endif
        enddo



        !****** outputting deformation (stress) analysis info *******
        !****************** (all data are relaxed) ******************
        ! initialize
        num_ua_part = 0 ! ua = 'u'nconstricted 'a'ffected
        num_ca_part = 0 ! ca = 'c'onstricted 'a'ffected
        num_uu_part = 0 ! uu = 'u'nconstricted 'u'naffected
        num_cu_part = 0 ! cu = 'c'onstricted 'u'naffected

        do i=1,npart
          ! 0: initially unconstricted part w/in the active region
          if(ini_constrict_label(i) .eq. 0) then

            ! 1, 2, or more: affected part
             if(affected_label(i) .ne. 0) then

              if(constrict_label(i) .eq. 0) then ! 0: unconstricted part
                status_label(i) = 1 ! ua = 1
                num_ua_part = num_ua_part + 1
                ua_part_idx(num_ua_part) = i
                ua_deform_data(num_ua_part) = tot_deform(i)
              elseif(constrict_label(i) .eq. 1) then ! 1: constricted part
                status_label(i) = 2 ! ca = 2
                num_ca_part = num_ca_part + 1
                ca_part_idx(num_ca_part) = i
                ca_deform_data(num_ca_part) = tot_deform(i)
              endif

            ! 0: unaffected part
            elseif(affected_label(i) .eq. 0) then

              if(constrict_label(i) .eq. 0) then ! 0: unconstricted part
                status_label(i) = 3 ! uu = 3
                num_uu_part = num_uu_part + 1
                uu_part_idx(num_uu_part) = i
                uu_deform_data(num_uu_part) = tot_deform(i)
              elseif(constrict_label(i) .eq. 1) then ! 1: constricted part
                status_label(i) = 4 ! cu = 4
                num_cu_part = num_cu_part + 1
                cu_part_idx(num_cu_part) = i
                cu_deform_data(num_cu_part) = tot_deform(i)
              endif

            endif

          endif !if(ini_constrict_label(i) .eq. 0)
        enddo ! do i=1,npart

        ! save relaxed configs for movies
        if(movie_switch .eq. 1) then
          write(2,*) npart, shrink_idx, cycle_idx
          write(2,*) nspart
          do i=1,npart
            ! x,y,vx,vy,radius,bdry_label,affected_label
            write(2,20) x(i),y(i),vx(i),vy(i),radius(i),bdry_label(i),
     &                  affected_label(i) 
          enddo
          write(2,*) ''
        endif

        ! save unconstricted affected (ua) part
        if(num_ua_part .gt. 0) then
          call sortmore_0_version(ua_deform_data,ua_part_idx,
     &                            maxnpart,num_ua_part)

          ave_ua_deform = 0.d0

          write(8,*) '# --- unconstricted affected (ua) deform ---:'
          write(8,*) '# num_ua_part'
          write(8,*) num_ua_part
          do i=1,num_ua_part
            ave_ua_deform = ave_ua_deform + ua_deform_data(i)
            write(8,*) ua_part_idx(i),ua_deform_data(i)
          enddo

          ave_ua_deform = ave_ua_deform / dble(num_ua_part)
          max_ua_deform = ua_deform_data(num_ua_part)

          write(8,*) '# ave =',ave_ua_deform
          write(8,*) '# max =',max_ua_deform
        else
          write(8,*) '# --- unconstricted affected (ua) deform ---:'
          write(8,*) '# num_ua_part'
          write(8,*) num_ua_part

          ave_ua_deform = 0.d0
          max_ua_deform = 0.d0

          write(8,*) '# ave =',ave_ua_deform,'(N/A)'
          write(8,*) '# max =',max_ua_deform,'(N/A)'
        endif

        ! save constricted affected (ca) part
        if(num_ca_part .gt. 0) then
          call sortmore_0_version(ca_deform_data,ca_part_idx,
     &                            maxnpart,num_ca_part)

          ave_ca_deform = 0.d0

          write(9,*) '# --- constricted affected (ca) deform ---:'
          write(9,*) '# num_ca_part'
          write(9,*) num_ca_part
          do i=1,num_ca_part
            ave_ca_deform = ave_ca_deform + ca_deform_data(i)
            write(9,*) ca_part_idx(i),ca_deform_data(i)
          enddo

          ave_ca_deform = ave_ca_deform / dble(num_ca_part)
          max_ca_deform = ca_deform_data(num_ca_part)

          write(9,*) '# ave =',ave_ca_deform
          write(9,*) '# max =',max_ca_deform
        else
          write(9,*) '# --- constricted affected (ca) deform ---:'
          write(9,*) '# num_ca_part'
          write(9,*) num_ca_part

          ave_ca_deform = 0.d0
          max_ca_deform = 0.d0

          write(9,*) '# ave =',ave_ca_deform,'(N/A)'
          write(9,*) '# max =',max_ca_deform,'(N/A)'
        endif

        ! save unconstricted unaffected (uu) part
        if(num_uu_part .gt. 0) then
          call sortmore_0_version(uu_deform_data,uu_part_idx,
     &                            maxnpart,num_uu_part)

          ave_uu_deform = 0.d0

          write(10,*) '# --- unconstricted unaffected (uu) deform ---:'
          write(10,*) '# num_uu_part'
          write(10,*) num_uu_part
          do i=1,num_uu_part
            ave_uu_deform = ave_uu_deform + uu_deform_data(i)
            write(10,*) uu_part_idx(i),uu_deform_data(i)
          enddo

          ave_uu_deform = ave_uu_deform /dble(num_uu_part)
          max_uu_deform = uu_deform_data(num_uu_part)

          write(10,*) '# ave =',ave_uu_deform
          write(10,*) '# max =',max_uu_deform
        else
          write(10,*) '# --- unconstricted unaffected (uu) deform ---:'
          write(10,*) '# num_uu_part'
          write(10,*) num_uu_part

          ave_uu_deform = 0.d0
          max_uu_deform = 0.d0

          write(10,*) '# ave =',ave_uu_deform,'(N/A)'
          write(10,*) '# max =',max_uu_deform,'(N/A)'
        endif

        ! save constricted unaffected (cu) part
        if(num_cu_part .gt. 0) then
          call sortmore_0_version(cu_deform_data,cu_part_idx,
     &                            maxnpart,num_cu_part)

          ave_cu_deform = 0.d0

          write(11,*) '# --- constricted unaffected (cu) deform ---:'
          write(11,*) '# num_cu_part'
          write(11,*) num_cu_part
          do i=1,num_cu_part
            ave_cu_deform = ave_cu_deform + cu_deform_data(i)
            write(11,*) cu_part_idx(i),cu_deform_data(i)
          enddo

          ave_cu_deform = ave_cu_deform / dble(num_cu_part)
          max_cu_deform = cu_deform_data(num_cu_part)

          write(11,*) '# ave =',ave_cu_deform
          write(11,*) '# max =',max_cu_deform
        else
          write(11,*) '# --- constricted unaffected (cu) deform ---:'
          write(11,*) '# num_cu_part'
          write(11,*) num_cu_part

          ave_cu_deform = 0.d0
          max_cu_deform = 0.d0

          write(11,*) '# ave =',ave_cu_deform,'(N/A)'
          write(11,*) '# max =',max_cu_deform,'(N/A)'
        endif

        ! save ave and max deformations
30      format(2(I10,3X),I10,3X,e21.14,3X,I10,3X,I5,3X,8(e21.14,3X))

        if(shrink_idx .gt. 0) then
          write(3,30) cycle_idx,shrink_idx,saved_step_count,saved_V,
     &                saved_num_trial,
     &                saved_num_unconstricted_affected_part,
     &                ave_ua_deform,max_ua_deform,
     &                ave_ca_deform,max_ca_deform,
     &                ave_uu_deform,max_uu_deform,
     &                ave_cu_deform,max_cu_deform
        endif
        !****************** (all data are relaxed) ******************
        !****** outputting deformation (stress) analysis info *******



        !------------- counting number of active part and --------------
        !-------- determine if quit the invagination procedure ---------
        !---------------------------------------------------------------
        num_active_part = 0
        do i=1,npart
          if(constrict_label(i) .eq. 0) then
            num_active_part = num_active_part + 1
          endif
        enddo
        write(*,*) 'current num_active_part =',num_active_part

        if(num_active_part .eq. 0) then
          ! quitting the invagination procedure if no more suitable part exist.
          write(*,*) 'No suitable unconstricted part,quit invagination.'
          goto 100
        endif
        !------------- counting number of active part and --------------
        !-------- determine if quit the invagination procedure ---------
        !---------------------------------------------------------------

        !------------------- normalizing tot_deform -------------------
        !------------------- normalizing tot_deform -------------------
        !------------------- normalizing tot_deform -------------------
        if(normalize_switch .eq. 0) then

          !----- find max_tot_deform for normalizing tot_deform -----
          ! initialization, some very large negative number
          max_tot_deform = -1.d5

          do i=1,npart
            if( constrict_label(i) .eq. 0 .and.
     &          tot_deform(i) .gt. max_tot_deform) then
              max_tot_deform = tot_deform(i)
            endif
          enddo
          write(*,*) 'max_tot_deform =',max_tot_deform
          !----------------------------------------------------------

          ! normalizing tot_deform by max_tot_deform
          do i=1,npart
            tot_deform(i) = tot_deform(i) / max_tot_deform
          enddo

        elseif(normalize_switch .eq. 1) then

          ! normalizing tot_deform by deformation of only one particle
          do i=1,npart
            tot_deform(i) = tot_deform(i) / ave_single_deform
          enddo

        endif
        !------------------- normalizing tot_deform -------------------
        !------------------- normalizing tot_deform -------------------
        !------------------- normalizing tot_deform -------------------

        !------------- constricting part by constriction_prob -----------
        ! initialization
        num_trial = 0

 2      continue

        ! initialization
        num_constricted_part = 0
        num_negative_prob = 0

        do i=1,npart
          if( constrict_label(i) .eq. 0 ) then

            ! evaluating constriction probability
            tot_deform(i) = min(tot_deform(i),cut_tot_deform)
            constriction_prob = max(cut, 1.d0 + beta * tot_deform(i)**3)
            constriction_prob = alpha * constriction_prob / (1.d0+beta)
            constriction_prob = constriction_prob /dble(num_active_part)

            !***************** for affected prob ******************
            if(affected_label(i) .eq. 1) then
              constriction_prob = constriction_prob * affected_rate1
            endif

            if(affected_label(i) .eq. 2) then
              constriction_prob = constriction_prob * affected_rate2
            endif

            if(affected_label(i) .eq. 3) then
              constriction_prob = constriction_prob * affected_rate3
            endif

            if(affected_label(i) .eq. 4) then
              constriction_prob = constriction_prob * affected_rate4
            endif
            !***************** for affected prob ******************

            call random_number(rand)
            trial_rand = rand

            if(constriction_prob .gt. trial_rand) then
              write(*,*) constriction_prob,'>',trial_rand
              write(*,*) 'normalized tot_deform',tot_deform(i)

              ! writing info on the screen.
              shrink_idx = shrink_idx + 1
              write(*,*) shrink_idx,'th constriction, part =',i

              ! shrinking selected particle size.
              radius(i) = radius(i) * constriction_rate

              ! labelling selected particle as constricted.
              constrict_label(i) = 1 ! constricted

              ! counting num_constricted_part in this do-loop
              num_constricted_part = num_constricted_part + 1
            elseif(affected_label(i) .eq. 1 .and.
     &             cycle_idx .ge. initiation_cycle_idx) then
              ! ### force initiators to constrict
              ! ### if cycle_idx >= initiation_cycle_idx
              write(*,*) cycle_idx,'>=',initiation_cycle_idx
              write(*,*) 'initiator (affected_label = 1) constricted'

              ! writing info on the screen.
              shrink_idx = shrink_idx + 1
              write(*,*) shrink_idx,'th constriction, part =',i

              ! shrinking selected particle size.
              radius(i) = radius(i) * constriction_rate

              ! labelling selected particle as constricted.
              constrict_label(i) = 1 ! constricted

              ! counting num_constricted_part in this do-loop
              num_constricted_part = num_constricted_part + 1
            elseif(affected_label(i) .eq. 2 .and.
     &             cycle_idx .ge. initiation_cycle_idx) then
              ! ### force initiators to constrict
              ! ### if cycle_idx >= initiation_cycle_idx
              write(*,*) cycle_idx,'>=',initiation_cycle_idx
              write(*,*) 'initiator (affected_label = 2) constricted'

              ! writing info on the screen.
              shrink_idx = shrink_idx + 1
              write(*,*) shrink_idx,'th constriction, part =',i

              ! shrinking selected particle size.
              radius(i) = radius(i) * constriction_rate

              ! labelling selected particle as constricted.
              constrict_label(i) = 1 ! constricted

              ! counting num_constricted_part in this do-loop
              num_constricted_part = num_constricted_part + 1
            endif

            ! counting num of part having negative constriction prob
            if(constriction_prob .le. 0.d0) then
              num_negative_prob = num_negative_prob + 1
            endif
          
          endif
        enddo

        ! if there is no more particles whose constriction prob > 0,
        ! quit the process
        if(num_negative_prob .eq. num_active_part) then
          write(*,*) 'All active part have negative constriction prob;'
          write(*,*) 'quit invagination.'
          write(3,*) '# All active part have negative constriction prob'
          write(3,*) '# quit invagination.'
          write(3,*) '#',num_active_part,'part are left unconstricted.'
          goto 100
        endif

        ! count num of trials
        num_trial = num_trial + 1

        ! make sure at least 1 part is constricted after the above do-loop
        ! as long as their constriction prob > 0.
        if(num_constricted_part .eq. 0) then
          !write(*,*) 'try more'
          goto 2
        endif
        !--------- end of constricting part by constriction_prob --------
        !________________________________________________________________
        !________________________________________________________________
        !_____________ end of selecting part for constriction ___________

        !_________ relaxation config containing constricted parts _______
        V_old = 0.d0
        V_rel = 1.d0
        step_count = 0

        do while(V_rel .gt. stop_V_rel)
          call verlet(x,y,vx,vy,npart,boxx,boxy,radius,dt,mass,V,
     &                rcut,rlist)

          step_count = step_count + 1

          if(mod(step_count,1000) .eq. 0) then

            V_rel = dabs(V - V_old)/V
            write(*,*) 'V_rel =',V_rel
            V_old = V

          endif
        enddo ! do while(V_rel .gt. stop_V_rel)

        !if(movie_switch .eq. 1) then
        !  write(2,*) npart, shrink_idx
        !  write(2,*) nspart
        !  do i=1,npart
        !    ! x,y,vx,vy,radius
        !    write(2,20) x(i),y(i),vx(i),vy(i),radius(i),bdry_label(i)
        !  enddo
        !  write(2,*) ''
        !endif

        !********************* for affected prob **********************
        num_constricted_affected_part = 0
        num_unconstricted_affected_part = 0

        do i=1,npart
          if(affected_label(i) .ne. 0 .and.
     &       constrict_label(i) .eq. 0) then
            num_unconstricted_affected_part
     &      = num_unconstricted_affected_part + 1
          elseif(affected_label(i) .ne. 0 .and.
     &           constrict_label(i) .eq. 1) then
            num_constricted_affected_part
     &      = num_constricted_affected_part + 1
          endif
        enddo

        if(num_unconstricted_affected_part 
     &     + num_constricted_affected_part .eq. num_affected_part) then
          !write(*,*) 'info of #_affected_part matches, ok!'
        else
          write(*,*) 'info of #_affected_part dismatches, stop!'
          stop
        endif
        !********************* for affected prob **********************

        saved_step_count = step_count
        saved_V = V
        saved_num_trial = num_trial
        saved_num_unconstricted_affected_part 
     &  = num_unconstricted_affected_part

        tot_step_count = tot_step_count + step_count

        write(8,*) '# above sorted info is used for shrink_idx becomes'
        write(8,*) shrink_idx
        write(8,*) ''
        write(9,*) '# above sorted info is used for shrink_idx becomes'
        write(9,*) shrink_idx
        write(9,*) ''
        write(10,*) '# above sorted info is used for shrink_idx becomes'
        write(10,*) shrink_idx
        write(10,*) ''
        write(11,*) '# above sorted info is used for shrink_idx becomes'
        write(11,*) shrink_idx
        write(11,*) ''

        write(*,*) 'relaxed after',step_count,'steps.'
        write(*,*) ''
        !____ end of relaxation config containing constricted parts _____

        cycle_idx = cycle_idx + saved_num_trial

        write(*,*) '~~~~~ cycle_idx  =',cycle_idx
        write(*,*) '~~~~~ shrink_idx =',shrink_idx
        write(*,*) '~~~~ (goto 1) ~~~~'
        write(*,*) ''

      goto 1
 100  continue
      !_____________________ invagination process _______________________



      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine verlet(x,y,vx,vy,npart,boxx,boxy,radius,dt,mass,V,
     &                  rcut,rlist)

      ! See Allen & Tildesley p.81, Eq.(3.19)
      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter (maxnpart=8192)
      double precision x(maxnpart),y(maxnpart),vx(maxnpart),vy(maxnpart)
      double precision radius(maxnpart)
      double precision boxx,boxy
      double precision dt,V
      double precision attract_fx(maxnpart),attract_fy(maxnpart)
      double precision attract_V
      double precision repel_fx(maxnpart),repel_fy(maxnpart),repel_V
      double precision attract_dampfx(maxnpart),attract_dampfy(maxnpart)
      double precision repel_dampfx(maxnpart),repel_dampfy(maxnpart)
      double precision mass(maxnpart)
      double precision old_ax(maxnpart),old_ay(maxnpart)
      double precision new_ax(maxnpart),new_ay(maxnpart)
      double precision rcut,rlist
      double precision bdryfx(maxnpart),bdryfy(maxnpart)

      do i=1,npart
        old_ax(i) = 0.d0
        old_ay(i) = 0.d0
        new_ax(i) = 0.d0
        new_ay(i) = 0.d0
      enddo

      call Vlist_repel_force(rcut,rlist,x,y,repel_fx,repel_fy,repel_V,
     &                       npart,boxx,boxy,radius)
      call iniVlist_attract_force(x,y,attract_fx,attract_fy,attract_V,
     &                            npart,boxx,boxy,radius)

      call Vlist_repel_dampforce(rcut,rlist,x,y,vx,vy,repel_dampfx,
     &                           repel_dampfy,npart,boxx,boxy,radius)
      call iniVlist_attract_dampforce(x,y,vx,vy,attract_dampfx,
     &                            attract_dampfy,npart,boxx,boxy,radius)

      call bdry_force(x,y,bdryfx,bdryfy,npart,boxx,boxy,radius)

      do i=1,npart
        old_ax(i) = ( attract_fx(i) + repel_fx(i)
     &            + attract_dampfx(i) + repel_dampfx(i) 
     &            + bdryfx(i) ) / mass(i)
        old_ay(i) = ( attract_fy(i) + repel_fy(i)
     &            + attract_dampfy(i) + repel_dampfy(i) 
     &            + bdryfy(i) ) / mass(i)
      enddo

      do i=1,npart
        x(i) = x(i) + vx(i) * dt + (old_ax(i) * dt * dt) / 2.d0
        y(i) = y(i) + vy(i) * dt + (old_ay(i) * dt * dt) / 2.d0
      enddo

      call Vlist_repel_force(rcut,rlist,x,y,repel_fx,repel_fy,repel_V,
     &                       npart,boxx,boxy,radius)
      call iniVlist_attract_force(x,y,attract_fx,attract_fy,attract_V,
     &                            npart,boxx,boxy,radius)

      call Vlist_repel_dampforce(rcut,rlist,x,y,vx,vy,repel_dampfx,
     &                           repel_dampfy,npart,boxx,boxy,radius)
      call iniVlist_attract_dampforce(x,y,vx,vy,attract_dampfx,
     &                            attract_dampfy,npart,boxx,boxy,radius)

      call bdry_force(x,y,bdryfx,bdryfy,npart,boxx,boxy,radius)

      do i=1,npart
        new_ax(i) = ( attract_fx(i) + repel_fx(i)
     &            + attract_dampfx(i) + repel_dampfx(i) 
     &            + bdryfx(i) ) / mass(i)
        new_ay(i) = ( attract_fy(i) + repel_fy(i)
     &            + attract_dampfy(i) + repel_dampfy(i) 
     &            + bdryfy(i) ) / mass(i)
      enddo

      do i=1,npart
        vx(i) = vx(i) + (new_ax(i) + old_ax(i)) / 2.d0 * dt
        vy(i) = vy(i) + (new_ay(i) + old_ay(i)) / 2.d0 * dt
      enddo

      V = attract_V + repel_V

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine iniVlist_attract_force(x,y,fx,fy,V,
     &                                  npart,boxx,boxy,radius)

      implicit none
      integer i,j,k
      integer npart,maxnpart,maxlist
      parameter(maxnpart=8192)
      parameter(maxlist=80000) ! increase it with increasing npart
      double precision x(maxnpart),y(maxnpart),fx(maxnpart),fy(maxnpart)
      double precision V
      double precision boxx,boxy
      double precision radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2,ff
      double precision x0(maxnpart),y0(maxnpart)
      integer nlist,point(maxlist),list(maxlist),jbegin,jend,jn
      integer iniVlist_update
      double precision list_range
      common /iniVlist1/ iniVlist_update
      common /iniVlist2/ list_range
      integer pbc_y_switch
      common /pbc/ pbc_y_switch
      save

      V = 0.d0
      
      do i=1,npart
        fx(i) = 0.d0
        fy(i) = 0.d0
      enddo

      if (iniVlist_update .eq. 1) then ! list need be updated

        do i=1,npart
          x0(i) = x(i)
          y0(i) = y(i)
        enddo

        nlist = 0

        do i=1,npart-1
          point(i) = nlist + 1
          do j=i+1,npart

            rx = x(i) - x(j)
            ry = y(i) - y(j)
            rx = rx - boxx*dnint(rx/boxx)
            if(pbc_y_switch .eq. 1) then
              ry = ry - boxy*dnint(ry/boxy)
            endif
            r2 = (rx * rx + ry * ry)
          
            sigma = radius(i) + radius(j)
            sigma2 = sigma * sigma

            if(r2 .lt. sigma2*list_range*list_range) then

              ! save neighbor list
              nlist = nlist + 1
              list(nlist) = j

            endif

          enddo
        enddo
        point(npart) = nlist + 1

      else ! list need not be updated

        do i=1,npart-1
          jbegin = point(i)
          jend = point(i+1) - 1
          if(jbegin .le. jend) then
            do jn=jbegin,jend
              j = list(jn)

              rx = x(i) - x(j)
              ry = y(i) - y(j)
              rx = rx - boxx*dnint(rx/boxx)
              if(pbc_y_switch .eq. 1) then
                ry = ry - boxy*dnint(ry/boxy)
              endif
              r2 = (rx*rx + ry*ry)

              sigma = radius(i) + radius(j)
              sigma2 = sigma*sigma

              if(r2 .gt. sigma2) then
                r = dsqrt(r2)

                ff = (1.d0 - r/sigma)

                fx(i) = fx(i) + (rx/r)*ff/sigma
                fy(i) = fy(i) + (ry/r)*ff/sigma
                fx(j) = fx(j) - (rx/r)*ff/sigma
                fy(j) = fy(j) - (ry/r)*ff/sigma

                V = V + 0.5d0*ff*ff
              endif
            enddo
          endif
        enddo

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine iniVlist_attract_dampforce(x,y,vx,vy,dampfx,dampfy,
     &                                      npart,boxx,boxy,radius)

      implicit none
      integer i,j,k
      integer npart,maxnpart,maxlist
      parameter(maxnpart=8192)
      parameter(maxlist=80000) ! increase it with increasing npart
      double precision x(maxnpart),y(maxnpart)
      double precision dampfx(maxnpart),dampfy(maxnpart)
      double precision vx(maxnpart),vy(maxnpart)
      double precision boxx,boxy,radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2
      double precision vxij,vyij
      double precision x0(maxnpart),y0(maxnpart)
      integer nlist,point(maxlist),list(maxlist),jbegin,jend,jn
      integer iniVlist_update
      double precision list_range
      common /iniVlist1/ iniVlist_update
      common /iniVlist2/ list_range
      double precision vn
      double precision gamman,gamman_sigma1
      common /damping/ gamman,gamman_sigma1
      integer pbc_y_switch
      common /pbc/ pbc_y_switch
      save

      !----------------------   setting parameters   -------------------
      ! the following line translates length scale from particle size 1
      ! to box size 1. Choose smaller particle diameter as the ruler.
      ! gamman = gamman_sigma1 / (2.d0 * radius(1))
      !-----------------------------------------------------------------

      do i=1,npart
        dampfx(i) = 0.d0
        dampfy(i) = 0.d0
      enddo

      if (iniVlist_update .eq. 1) then ! list need be updated

        do i=1,npart
          x0(i) = x(i)
          y0(i) = y(i)
        enddo

        nlist = 0

        do i=1,npart-1
          point(i) = nlist + 1
          do j=i+1,npart

            rx = x(i) - x(j)
            ry = y(i) - y(j)
            rx = rx - boxx*dnint(rx/boxx)
            if(pbc_y_switch .eq. 1) then
              ry = ry - boxy*dnint(ry/boxy)
            endif
            r2 = (rx*rx + ry*ry)

            vxij = vx(i) - vx(j)
            vyij = vy(i) - vy(j)

            sigma = radius(i) + radius(j)
            sigma2 = sigma*sigma

            if(r2.lt.sigma2*list_range*list_range) then

              ! save neighbor list
              nlist = nlist + 1
              list(nlist) = j

            endif

          enddo
        enddo
        point(npart) = nlist + 1

      else ! list need not be updated

        do i=1,npart-1
          jbegin = point(i)
          jend = point(i+1) - 1
          if(jbegin .le. jend) then
            do jn=jbegin,jend
              j = list(jn)

              rx = x(i) - x(j)
              ry = y(i) - y(j)
              rx = rx - boxx*dnint(rx/boxx)
              if(pbc_y_switch .eq. 1) then
                ry = ry - boxy*dnint(ry/boxy)
              endif
              r2 = (rx*rx + ry*ry)

              vxij = vx(i) - vx(j)
              vyij = vy(i) - vy(j)

              sigma = radius(i) + radius(j)
              sigma2 = sigma*sigma

              if(r2 .gt. sigma2) then
                r = dsqrt(r2) ! square root is expensive operation

                vn = (vxij*rx + vyij*ry)/r

                ! damping force is minus sign
                dampfx(i) = dampfx(i) - vn*(rx/r)*gamman
                dampfy(i) = dampfy(i) - vn*(ry/r)*gamman
                dampfx(j) = dampfx(j) + vn*(rx/r)*gamman
                dampfy(j) = dampfy(j) + vn*(ry/r)*gamman
              endif
            enddo
          endif
        enddo
        
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine iniVlist_attract_tot_deform(x,y,npart,boxx,boxy,radius,
     &                                       tot_deform)

      implicit none
      integer i,j,k
      integer npart,maxnpart,maxlist
      parameter(maxnpart=8192)
      parameter(maxlist=80000) ! increase it with increasing npart
      double precision x(maxnpart),y(maxnpart)
      double precision boxx,boxy
      double precision radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2
      double precision x0(maxnpart),y0(maxnpart)
      integer nlist,point(maxlist),list(maxlist),jbegin,jend,jn
      integer iniVlist_update
      double precision list_range
      common /iniVlist1/ iniVlist_update
      common /iniVlist2/ list_range
      double precision ff,tot_deform(maxnpart)
      integer pbc_y_switch
      common /pbc/ pbc_y_switch
      save

      ! initialization
      do i=1,npart
        tot_deform(i) = 0.d0
      enddo

      if (iniVlist_update .eq. 1) then ! list need be updated

        do i=1,npart
          x0(i) = x(i)
          y0(i) = y(i)
        enddo

        nlist = 0

        do i=1,npart-1
          point(i) = nlist + 1
          do j=i+1,npart

            rx = x(i) - x(j)
            ry = y(i) - y(j)
            rx = rx - boxx*dnint(rx/boxx)
            if(pbc_y_switch .eq. 1) then
              ry = ry - boxy*dnint(ry/boxy)
            endif
            r2 = (rx * rx + ry * ry)
          
            sigma = radius(i) + radius(j)
            sigma2 = sigma * sigma

            if(r2 .lt. sigma2*list_range*list_range) then

              ! save neighbor list
              nlist = nlist + 1
              list(nlist) = j

            endif

          enddo
        enddo
        point(npart) = nlist + 1

      else ! list need not be updated

        do i=1,npart-1
          jbegin = point(i)
          jend = point(i+1) - 1
          if(jbegin .le. jend) then
            do jn=jbegin,jend
              j = list(jn)

              rx = x(i) - x(j)
              ry = y(i) - y(j)
              rx = rx - boxx*dnint(rx/boxx)
              if(pbc_y_switch .eq. 1) then
                ry = ry - boxy*dnint(ry/boxy)
              endif
              r2 = (rx*rx + ry*ry)

              sigma = radius(i) + radius(j)
              sigma2 = sigma*sigma

              if(r2 .gt. sigma2) then
                r = dsqrt(r2)
                ff = (r/sigma - 1.d0) ! ff < 0: contraction; ff > 0: expansion
                ff = ff * (dabs(rx)/r)**2 ! contribution of ff in the xx dir.

                tot_deform(i) = tot_deform(i) + ff
                tot_deform(j) = tot_deform(j) + ff
              endif
            enddo
          endif
        enddo

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sortmore_0_version(array,A_follow_array,max,num)

      integer max,num,i,j,locsm
      double precision array(max),small
      integer A_follow_array(max),dummy_A

      !!! Be careful, the formats (double precision or integer) have to 
      !!! be exactly the same as input arrays, otherwise it will go wrong

      do i=1,num-1

         ! 'locsm','small' are dummy variables; locsm stores index of local min,
         ! small stores its value and so on.
         ! Sort 'array', then all the other matrices are sorted according to 
         ! 'array''s order.

         ! initialize, if real minimun is not in the remaining elements; then 
         ! the data stored here will be used.

         ! major matrix to be sorted
         locsm = i
         small = array(locsm)

         ! 1st follower matrix
         dummy_A = A_follow_array(locsm)
         
         ! find real minimum from remianing elements
         do j=i+1,num

            if ( array(j) .lt. small ) then
               ! major matrix to be sorted
               locsm = j ! index
               small = array(j) ! dummy value of 'array'

               ! 1st follower matrix
               dummy_A = A_follow_array(j) ! dummy value of 'A_follow_array'
            endif
         enddo

         ! exchange their values

         ! major matrix to be sorted
         array(locsm)=array(i)
         array(i) = small

         ! 1st follower matrix
         A_follow_array(locsm) = A_follow_array(i)
         A_follow_array(i) = dummy_A

      enddo

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check(rcut,rlist,x,y,x0,y0,npart,update)

      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter(maxnpart=8192)
      double precision x(maxnpart),y(maxnpart)
      double precision x0(maxnpart),y0(maxnpart)
      double precision rcut,rlist
      integer update
      double precision dispmax
      save

      dispmax = 0.d0

      do i=1,npart

        dispmax = max(dabs(x(i) - x0(i)),dispmax)
        dispmax = max(dabs(y(i) - y0(i)),dispmax)

      enddo

      dispmax = 2.d0 * 1.414213562373095d0 * dispmax

      if (dispmax .gt. (rlist-rcut)) then
        update = 1 ! list need be updated
      else
        update = 0 ! list need not be updated
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Vlist_repel_force(rcut,rlist,x,y,fx,fy,
     &                             V,npart,boxx,boxy,radius)

      implicit none
      integer i,j,k
      integer npart,maxnpart,maxlist
      parameter(maxnpart=8192)
      parameter(maxlist=80000) ! increase it with increasing npart
      double precision x(maxnpart),y(maxnpart),fx(maxnpart),fy(maxnpart)
      double precision V
      double precision boxx,boxy
      double precision radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2,ff
      double precision rcut,rlist,rlist2
      double precision x0(maxnpart),y0(maxnpart)
      integer update,nlist,point(maxlist),list(maxlist),jbegin,jend,jn
      integer pbc_y_switch
      common /pbc/ pbc_y_switch
      save

      rlist2 = rlist * rlist

      V = 0.d0
      
      do i=1,npart
        fx(i) = 0.d0
        fy(i) = 0.d0
      enddo

      call check(rcut,rlist,x,y,x0,y0,npart,update)

      if (update .eq. 1) then ! list need be updated

        do i=1,npart
          x0(i) = x(i)
          y0(i) = y(i)
        enddo

        nlist = 0

        do i=1,npart-1
          point(i) = nlist + 1
          do j=i+1,npart

            rx = x(i) - x(j)
            ry = y(i) - y(j)
            rx = rx - boxx*dnint(rx/boxx)
            if(pbc_y_switch .eq. 1) then
              ry = ry - boxy*dnint(ry/boxy)
            endif
            r2 = (rx * rx + ry * ry)

            sigma = radius(i) + radius(j)
            sigma2 = sigma * sigma

            if (r2 .lt. rlist2) then

              nlist = nlist + 1
              list(nlist) = j

              if(r2 .lt. sigma2) then

                r = dsqrt(r2)
           
                ff = (1.d0 - r/sigma)

                fx(i) = fx(i) + (rx/r)*ff/sigma
                fy(i) = fy(i) + (ry/r)*ff/sigma
                fx(j) = fx(j) - (rx/r)*ff/sigma
                fy(j) = fy(j) - (ry/r)*ff/sigma

                V = V + 0.5d0 * ff * ff

              endif

            endif

          enddo
        enddo
        point(npart) = nlist + 1

      else ! list need not be updated

        do i=1,npart-1
          jbegin = point(i)
          jend = point(i+1) - 1
          if(jbegin .le. jend) then
            do jn=jbegin,jend
              j = list(jn)

              rx = x(i) - x(j)
              ry = y(i) - y(j)
              rx = rx - boxx*dnint(rx/boxx)
              if(pbc_y_switch .eq. 1) then
                ry = ry - boxy*dnint(ry/boxy)
              endif
              r2 = (rx*rx + ry*ry)

              sigma = radius(i) + radius(j)
              sigma2 = sigma*sigma

              if(r2 .lt. sigma2) then

                r = dsqrt(r2)

                ff = (1.d0 - r/sigma)

                fx(i) = fx(i) + (rx/r)*ff/sigma
                fy(i) = fy(i) + (ry/r)*ff/sigma
                fx(j) = fx(j) - (rx/r)*ff/sigma
                fy(j) = fy(j) - (ry/r)*ff/sigma

                V = V + 0.5d0*ff*ff

              endif
            enddo
          endif
        enddo

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Vlist_repel_dampforce(rcut,rlist,x,y,vx,vy,dampfx,
     &                                 dampfy,npart,boxx,boxy,radius)

      implicit none
      integer i,j,k
      integer npart,maxnpart,maxlist
      parameter(maxnpart=8192)
      parameter(maxlist=80000) ! increase it with increasing npart
      double precision x(maxnpart),y(maxnpart)
      double precision dampfx(maxnpart),dampfy(maxnpart)
      double precision vx(maxnpart),vy(maxnpart)
      double precision boxx,boxy,radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2
      double precision vxij,vyij
      double precision rcut,rlist,rlist2
      double precision x0(maxnpart),y0(maxnpart)
      integer update,nlist,point(maxlist),list(maxlist),jbegin,jend,jn
      double precision vn
      double precision gamman,gamman_sigma1
      common /damping/ gamman,gamman_sigma1
      integer pbc_y_switch
      common /pbc/ pbc_y_switch
      save

      !----------------------   setting parameters   -------------------
      ! the following line translates length scale from particle size 1
      ! to box size 1. Choose smaller particle diameter as the ruler.
      ! gamman = gamman_sigma1 / (2.d0 * radius(1))
      !-----------------------------------------------------------------

      rlist2 = rlist*rlist

      do i=1,npart
        dampfx(i) = 0.d0
        dampfy(i) = 0.d0
      enddo

      call check(rcut,rlist,x,y,x0,y0,npart,update)

      if (update .eq. 1) then ! list need be updated

        do i=1,npart
          x0(i) = x(i)
          y0(i) = y(i)
        enddo

        nlist = 0

        do i=1,npart-1
          point(i) = nlist + 1
          do j=i+1,npart

            rx = x(i) - x(j)
            ry = y(i) - y(j)
            rx = rx - boxx*dnint(rx/boxx)
            if(pbc_y_switch .eq. 1) then
              ry = ry - boxy*dnint(ry/boxy)
            endif
            r2 = (rx*rx + ry*ry)

            vxij = vx(i) - vx(j)
            vyij = vy(i) - vy(j)

            sigma = radius(i) + radius(j)
            sigma2 = sigma*sigma

            if (r2 .lt. rlist2) then

              nlist = nlist + 1
              list(nlist) = j

              ! only contact produces damping force
              if(r2.lt.sigma2) then

                r = dsqrt(r2) ! square root is expensive operation

                vn = (vxij*rx + vyij*ry)/r

                ! damping force is minus sign
                dampfx(i) = dampfx(i) - vn*(rx/r)*gamman
                dampfy(i) = dampfy(i) - vn*(ry/r)*gamman
                dampfx(j) = dampfx(j) + vn*(rx/r)*gamman
                dampfy(j) = dampfy(j) + vn*(ry/r)*gamman

              endif

            endif

          enddo
        enddo
        point(npart) = nlist + 1

      else ! list need not be updated

        do i=1,npart-1
          jbegin = point(i)
          jend = point(i+1) - 1
          if(jbegin .le. jend) then
            do jn=jbegin,jend
              j = list(jn)

              rx = x(i) - x(j)
              ry = y(i) - y(j)
              rx = rx - boxx*dnint(rx/boxx)
              if(pbc_y_switch .eq. 1) then
                ry = ry - boxy*dnint(ry/boxy)
              endif
              r2 = (rx*rx + ry*ry)

              vxij = vx(i) - vx(j)
              vyij = vy(i) - vy(j)

              sigma = radius(i) + radius(j)
              sigma2 = sigma*sigma

              ! only contact produces damping force
              if(r2.lt.sigma2) then

                r = dsqrt(r2) ! square root is expensive operation

                vn = (vxij*rx + vyij*ry)/r

                ! damping force is minus sign
                dampfx(i) = dampfx(i) - vn*(rx/r)*gamman
                dampfy(i) = dampfy(i) - vn*(ry/r)*gamman
                dampfx(j) = dampfx(j) + vn*(rx/r)*gamman
                dampfy(j) = dampfy(j) + vn*(ry/r)*gamman

              endif
            enddo
          endif
        enddo
        
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Vlist_repel_tot_deform(rcut,rlist,x,y,npart,
     &                                  boxx,boxy,radius,tot_deform)

      implicit none
      integer i,j,k
      integer npart,maxnpart,maxlist
      parameter(maxnpart=8192)
      parameter(maxlist=80000) ! increase it with increasing npart
      double precision x(maxnpart),y(maxnpart)
      double precision boxx,boxy
      double precision radius(maxnpart)
      double precision rx,ry,r,r2,sigma,sigma2
      double precision rcut,rlist,rlist2
      double precision x0(maxnpart),y0(maxnpart)
      integer update,nlist,point(maxlist),list(maxlist),jbegin,jend,jn
      double precision ff,tot_deform(maxnpart)
      integer pbc_y_switch
      common /pbc/ pbc_y_switch
      save

      rlist2 = rlist * rlist
      
      ! initialization
      do i=1,npart
        tot_deform(i) = 0.d0
      enddo

      call check(rcut,rlist,x,y,x0,y0,npart,update)

      if (update .eq. 1) then ! list need be updated

        do i=1,npart
          x0(i) = x(i)
          y0(i) = y(i)
        enddo

        nlist = 0

        do i=1,npart-1
          point(i) = nlist + 1
          do j=i+1,npart

            rx = x(i) - x(j)
            ry = y(i) - y(j)
            rx = rx - boxx*dnint(rx/boxx)
            if(pbc_y_switch .eq. 1) then
              ry = ry - boxy*dnint(ry/boxy)
            endif
            r2 = (rx * rx + ry * ry)

            sigma = radius(i) + radius(j)
            sigma2 = sigma * sigma

            if (r2 .lt. rlist2) then

              nlist = nlist + 1
              list(nlist) = j

              if(r2 .lt. sigma2) then

                r = dsqrt(r2)
                ff = (r/sigma - 1.d0) ! ff < 0: contraction; ff > 0: expansion

                tot_deform(i) = tot_deform(i) + ff
                tot_deform(j) = tot_deform(j) + ff

              endif

            endif

          enddo
        enddo
        point(npart) = nlist + 1

      else ! list need not be updated

        do i=1,npart-1
          jbegin = point(i)
          jend = point(i+1) - 1
          if(jbegin .le. jend) then
            do jn=jbegin,jend
              j = list(jn)

              rx = x(i) - x(j)
              ry = y(i) - y(j)
              rx = rx - boxx*dnint(rx/boxx)
              if(pbc_y_switch .eq. 1) then
                ry = ry - boxy*dnint(ry/boxy)
              endif
              r2 = (rx*rx + ry*ry)

              sigma = radius(i) + radius(j)
              sigma2 = sigma*sigma

              if(r2 .lt. sigma2) then

                r = dsqrt(r2)
                ff = (r/sigma - 1.d0) ! ff < 0: contraction; ff > 0: expansion

                tot_deform(i) = tot_deform(i) + ff
                tot_deform(j) = tot_deform(j) + ff

              endif
            enddo
          endif
        enddo

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine bdry_force(x,y,bdryfx,bdryfy,npart,boxx,boxy,radius)

      implicit none
      integer i,j,k
      integer npart,maxnpart
      parameter(maxnpart=8192)
      integer bdry_label(maxnpart)
      double precision saved_ini_x(maxnpart),saved_ini_y(maxnpart)
      double precision x(maxnpart),y(maxnpart)
      double precision bdryfx(maxnpart),bdryfy(maxnpart)
      double precision boxx,boxy
      double precision radius(maxnpart),radius1
      double precision rx,ry,r,r2,sigma,ff
      integer bdry_force_switch
      double precision bdry_y_thick
      double precision top_bdry,bottom_bdry
      integer tot_top_bdry_part,tot_bottom_bdry_part
      common /boundary0/ bdry_force_switch
      common /boundary1/ radius1,bdry_y_thick,top_bdry,bottom_bdry
      common /boundary2/ bdry_label
      double precision spring_k,spring_k_sigma1
      common /spring/ spring_k,spring_k_sigma1
      save

      !----------------------   setting parameters   -------------------
      ! the following line translates length scale from particle size 1
      ! to box size 1. Choose smaller particle diameter as the ruler.

      ! spring_k = spring_k_sigma1 / (2.d0 * radius1)**2
      !-----------------------------------------------------------------

      if(bdry_force_switch .eq. 1) then ! assigning bdry particles

        ! initialize
        !top_bdry = boxy/2.d0
        !bottom_bdry = -boxy/2.d0

        tot_top_bdry_part = 0
        tot_bottom_bdry_part = 0

        ! assign bdry part
        do i=1,npart
          bdry_label(i) = 0 ! 0: initialize. non-bdry_particles
          saved_ini_x(i) = 0.d0
          saved_ini_y(i) = 0.d0

          if(y(i) .ge. top_bdry .and.
     &       y(i) .lt. top_bdry + bdry_y_thick * radius1) then
            bdry_label(i) = 1 ! 1: top bdry_particles
            saved_ini_x(i) = x(i)
            saved_ini_y(i) = y(i)
            tot_top_bdry_part = tot_top_bdry_part + 1
          elseif(y(i) .le. bottom_bdry .and.
     &           y(i) .gt. bottom_bdry - bdry_y_thick * radius1) then
            bdry_label(i) = 2 ! 2: bottom bdry_particles
            saved_ini_x(i) = x(i)
            saved_ini_y(i) = y(i)
            tot_bottom_bdry_part = tot_bottom_bdry_part + 1
          endif
        enddo

        write(*,*) 'tot_top_bdry_part =',tot_top_bdry_part
        write(*,*) 'tot_bottom_bdry_part =',tot_bottom_bdry_part

      else ! turn off assigning bdry particles, calculating forces

        ! initialize
        do i=1,npart
          bdryfx(i) = 0.d0
          bdryfy(i) = 0.d0
        enddo

        ! calculate force
        do i=1,npart

          if(bdry_label(i) .ne. 0) then ! bdry particles
            rx = x(i) - saved_ini_x(i)
            ry = y(i) - saved_ini_y(i)

            bdryfx(i) = -spring_k * rx
            bdryfy(i) = -spring_k * ry
          endif

        enddo

      endif

      return
      end
