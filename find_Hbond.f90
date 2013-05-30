Program Hbonds
   !
   use gen_lib
   !
   implicit none        
   !
   integer, parameter   :: namax = 2            !!!!Has only been tested with water and Ions!!!!
   integer, parameter   :: shell_max = 20       !max number of Oxygen assumed around the defect
   integer, parameter   :: numHref = 4          !max number H's associated with a Oxygen
   !
   real(DP), parameter  :: pi    = 4.d0*atan(1.d0),& !Pi
                           ao    = 0.52917720859_DP
   !
   !
   real(DP)                :: aprim(3,3),    &  !lattice prim vectors
                              apinv(3,3),    &  !inverse of prim vectors, for sca
                              box,           &  !length of the region 
                              Hbox,          &  !length of half the region 
                              omega=0.0,     &  !volume of cell
                              rcut,          &  !Covalent radius
                              Hbond_cut,     &  !Hbond cut-off distance
                              Hbond_angle,   &  !Hbond Angle cut-off
                              ref_cos,       &  !Cosine of the above Hbond_angle
                              time,          &  !current time of the step
                              dummy 
   !
   integer                 :: Ostar,         &  !Current Water molecule (oxygen) where Hbonds will be calculated 
                              Ospecies,      &  !what species number is the Oxygen
                              Hspecies,      &  !what species number is the Hydrogen
                              nat,           &  !total number of atoms
                              ntsp,          &  !number of total atomic species, limit namax
                              nsp(namax),    &  !number of each atomic species     
                              mic=1,         &  !Scaled length of the box region, (for future use)
                              step,          &  !Desired step 
                              Oshell_count,  &  !Total Number of Oxygen mol inside the Hbond_cut radius
                                                !WARNING: This establishes a NEW index set for the 
                                                !Oxygens inside Oshell
                              readstep,      &  !Current readstep in *pos and *cel file
                              ierror,        &  !error check
                              na,ns,         &  !Atomic, species index
                              i,j,k
   !
   real(DP),allocatable    :: tau(:,:,:),    &  !atomic positions (dim, index, atomic-species)
                              stau(:,:,:),   &  !scaled atomic positions (dim,index, atomic-species)
                              sOshell(:,:),  &  !Oxygens inside the Hbond_cut length (dim, Oshell-index) (assume index < 20)
                              sHshell(:,:,:)    !Hydrogens of Oxygens inside the HBond_cut length (dim, Oshell-index, H#)
                                                !WARNING: The indexes used above are for the Oshell index within the Oshell radius

   integer,allocatable     :: Href(:,:),     &  !Hydrogens associations with Oxygens (O-index, ref1-ref2-ref3 )   
                              O_index(:)        !The Orginal indexes of those oxygens in the Hbond_cut length
   !
   character(len=256)       :: filePos, fileCel, dum
   !
   namelist /system/ filePos, fileCel, ntsp, nsp, step, rcut, Hbond_cut, Hbond_angle, Ospecies, Hspecies

   !
   !
   !initialize
   CALL init_hbonds()
   !
   !Open Files
   open(unit=1,   file=(trim(filePos)), status='old')
   open(unit=2,   file=(trim(fileCel)), status='old')
   open(unit=11, file='nearest-O.dat', status='unknown')
   open(unit=12, file='Href.dat', status='unknown')
   !
   write(11,'("#Step:   ", I5)') step
   write(11,'("#Oxygen: ", I5)') Ostar
   write(12,'("#Step:   ", I5)') step
   write(11,'("#---------------------------------------------------------------------------------------")')
   !
   !******************************************************************
   Main_loop: do 
      !
      !-----read *.pos file-----
      read(1,*, iostat=ierror) readstep, time
      !
      !-----Check for end of *.pos file-----
      if (ierror < 0) then 
         write(*,*) ''
         write(*,*) ' ERROR: Position File Not Cannot be read' 
         stop
      endif
      !
      !Read in current position values
      do ns=1,ntsp,1
         do na=1,nsp(ns),1
            read(1,*) (tau(j,na,ns),j=1,3)   
         enddo
      enddo
      !
      !-----read *cel file------
      read(2, *) dummy, dummy
      do i=1,3,1
         read(2, *) (aprim(i,j),j=1,3)
      enddo
      !
      !----Check for desired step------
      if (readstep < step) then
         cycle Main_loop
      endif
      !
      write(*,*) " "
      write(*, '(1X, "Current Water Molecule: ", 2X, I3)') Ostar
      write(*,'(1X, 3(F10.7, 2X))') (tau(j,Ostar,Ospecies),j=1,3)
      write(*,*) ''
      write(*,'(1X, "Program start... ")') 
      write(*,*) ''
      !
      write(11,'("#Current molecule position: ", 3(F14.8))') (tau(k,Ostar,Ospecies),k=1,3)
      write(11,'("#---------------------------------------------------------------------------------------")')
      !
      !calculate the inverse
      call invert(aprim, apinv, omega)
      !
      !region dimensions for THIS loop
      box  = mic*omega**(1./3.)
      Hbox = box/2.
      !
      !Convert to scaled positions 
      do ns=1,ntsp,1
         do na=1,nsp(ns),1
            call r_to_s( tau(1:3,na,ns), stau(1:3,na,ns), apinv )
         enddo
      enddo
      !
      !Find Href array
      CALL find_Href(stau, Href)
      !
      CALL Find_Hbonds(Ostar, stau, Href)
      ! 
      !-----Check for stop condition in *.pos file-------
      if(readstep == step) then
         write(*,*) ''
         write(*,'(1X, "...Program Complete")') 
         write(*,*) ''
         exit
      endif !
      !
   enddo Main_loop
   !******************************************************************
   !
   close(1)
   close(2)
   close(11)
   close(12)
   !
   !
   Contains
      !
      !
      Subroutine Find_Hbonds(Ostar, stau,Href)
         !
         use gen_lib
         !
         implicit none
         !
         real(DP), intent(in) :: stau(:,:,:)
         integer, intent(in)  :: Ostar, Href(:,:)
         !
         real(DP)             :: rdist(3),            &  !square of the components distance in real coordinates
                                 r2                      !squared distanced
         !
         real(DP)             :: r2AB, r2AC, r2BC,    &  !for the law of cosines variables
                                 cos_angle
         !
         integer              :: na,ns,no,            &  !index: atoms, species, oxygen-shell
                                 num_donate,          &  !local donated Hbonds, for each configuration
                                 num_accept,          &  !local accepted Hbonds, for each configuration
                                 i,j,k    
         !
         !
         num_donate = 0
         num_accept = 0
         !
         !Reset the Oshell_count
         Oshell_count = 0
         sOshell(:,:) = 0.0_DP
         sHshell(:,:,:) = 0.0_DP
         !
         !
         !Loop over all Oxygen atoms (cycle when Ostar) find
         !All atoms within Hbond length
         Oloop: do na=1,nsp(Ospecies),1
            !
            if (na == Ostar) cycle
            !
            !Oxygen-Oxygen Distance
            CALL get_rdist(stau(1:3,Ostar,1),stau(1:3,na,1),rdist,mic,aprim)
            r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
            !
            if (r2 < (Hbond_cut/ao)**2 ) then 
               ! 
               !Current, running, number of of atoms in first Oshell 
               Oshell_count = Oshell_count + 1
               !
               !Create O_index (a array to keep track of the previous indexes of
               !those Oxygens inside the Hbond_cut length
               O_index(Oshell_count) = na
               !
               !Create sOshell
               sOshell(:,Oshell_count) = stau(1:3,na,1)
               !
               !Create sHshell, using Href(O-index, H#)
               Hloop: do j=1,numHref,1
               !
               !cycle if Href is 0 (no Hydrogen association) 
               if (Href(na,j) == 0) cycle Hloop 
                  !
                  !Assign Hydrogen to the sHshell array
                  sHshell(:,Oshell_count,j) = stau(1:3,Href(na,j),2)
                  !
               enddo Hloop
               ! 
            endif
            !
         enddo Oloop
         !
         write(*,*) " "
         write(*,'(2X, "Oshell_count: ", I5)') Oshell_count
         write(*,*) " "
         !
         !------------------------------------------
         !Loop over all donated H bonds (for each Hprime), check angle
         !write to coord-num.dat-> readstep -- accepted-Hbonds donated-Hbonds
         !r2AB = OstarHprime
         !r2BC = HprimeOw
         !r2AC = OstarOw
         donate_H: do i=1,numHref,1
            !
            !cycle if Href is 0 (no Hydrogen association) 
            if (Href(Ostar,i) == 0) cycle donate_H
            !
            CALL  get_rdist(stau(1:3,Ostar,1),stau(1:3,Href(Ostar,i),2),rdist,mic,aprim)
            r2AB = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
            !
            donate_O: do j=1,Oshell_count,1
               !
               CALL  get_rdist(stau(1:3,Href(Ostar,i),2),sOshell(1:3,j),rdist,mic,aprim)
               r2BC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               CALL  get_rdist(stau(1:3,Ostar,1),sOshell(1:3,j),rdist,mic,aprim)
               r2AC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               cos_angle = ( (r2AB + r2BC) - r2AC )/(2.0*sqrt(r2AB*r2BC))
               ! 
               !confirm Hbond
               if (cos_angle <= ref_cos) then
                  !
                  num_donate = num_donate + 1
                  write(11,'( 2X, I5, 3X, 3(F20.15), "  D")') O_index(j), (tau(k,O_index(j),Ospecies), k=1,3)
                  !
               endif  
               !
               enddo donate_O
         enddo donate_H
         !------------------------------------------
         ! 
         ! 
         !------------------------------------------
         !Loop over all  accepted H bonds, check angle
         !write to coord-num.dat -> readstep -- accepted-Hbonds donated-Hbonds
         !r2AB = OstarHw
         !r2BC = HwOw
         !r2AC = OstarOw
         accept_O: do i=1,Oshell_count,1
            accept_H: do j=1,numHref  !loop over each Hydrogen of normal water molecules
               ! 
               CALL  get_rdist(stau(1:3,Ostar,1),sHshell(1:3,i,j),rdist,mic,aprim)
               r2AB = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               ! 
               CALL  get_rdist(sHshell(1:3,i,j),sOshell(1:3,i),rdist,mic,aprim)
               r2BC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               CALL  get_rdist(stau(1:3,Ostar,1),sOshell(1:3,i),rdist,mic,aprim)
               r2AC = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               cos_angle = ( (r2AB + r2BC) - r2AC )/(2.0*sqrt(r2AB*r2BC))
               ! 
               if (cos_angle <= ref_cos) then
                  !
                  num_accept = num_accept + 1
                  write(11,'( 2X, I5, 3X, 3(F20.15), "  A")') O_index(i), (tau(k,O_index(i),Ospecies), k=1,3)
                  !
               endif  
               !
            enddo accept_H
         enddo accept_O
         !------------------------------------------
         !         
         write(*,*) " "
         write(*,*) "------------------------------------------------"
         write(*,*) " Total number accepted bonds: ", num_accept
         write(*,*) " Total number donated bonds:  ", num_donate
         write(*,*) "------------------------------------------------"
         write(*,*) " "
         ! 
         !write(*,*) "
         !----------------------------------------------------------------------------- 
         !
         return

      End Subroutine Find_Hbonds
      !
      !
      Subroutine Find_Href(stau, Href)
         !
         implicit none
         !
         real(DP), intent(in)    :: stau(:,:,:)
         integer, intent(out)    :: Href(:,:)
         !
         real(DP)             :: rdist(3),      &  !square of the components distance in real coordinates
                                 r2                !total squared distance
         !
         integer              :: H1,H2,H3,H4,   &  !Counters for the Hydrogen around O
                                 i,j,k             !general index
         !
         !Loop over all Oxygens
         Oloop: do i=1,nsp(Ospecies),1
            !
            !write(145,*) "Oxygen: ", i
            !
            !Hydrogen counters, set to zero for each Oxygen
            H1 = 0
            H2 = 0
            H3 = 0
            H4 = 0
            !
            Hloop: do j=1,nsp(Hspecies),1
               !
               !Get the real distance in each coordinate direction
               CALL get_rdist(stau(1:3,i,1),stau(1:3,j,2),rdist,mic,aprim)
               r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
               !
               if (r2 < (rcut/ao)**2 ) then
                  !How many Hydrogen bonds are already associated with the O*
                  if ( (H1==0) .and. (H2==0) .and. (H3==0) .and. (H4==0) ) then
                     H1 = j 
                     !write(145,*) 'Hydrogen1 ', i, j, ' ', (sqrt(r2)*ao)
                  elseif ( (H2==0) .and. (H3==0) .and. (H4==0) ) then
                     H2 = j 
                     !write(145,*) 'Hydrogen2 ', i, j, ' ', (sqrt(r2)*ao)
                  elseif ( (H3==0) .and. (H4==0) ) then
                     H3 = j 
                     !write(145,*) 'Hydrogen3 ', i, j, ' ', (sqrt(r2)*ao)
                  elseif ( (H4==0) )then !Error
                     !A warning message will be generated below
                     H4 = j 
                     !write(145,*) 'Hydrogen4', i, j, ' ', (sqrt(r2)*ao)
                  endif
               endif
               !
            enddo Hloop
            !
            !---------------------------------------------------
            ! Check to see how many H atoms are associated 
            ! with this Oxygen
            !
            !  1 H (Error)    : H1  = 0 H2  =0 H3  =0 H4  =0 
            !  OH             : H1 /= 0 H2  =0 H3  =0 H4  =0
            !  Normal Water   : H1 /= 0 H2 /=0 H3  =0 H4  =0 
            !  H3O            : H1 /= 0 H2 /=0 H3 /=0 H4  =0
            !  4+H (Error)    : H1 /= 0 H2 /=0 H3 /=0 H4 /=0
            !
            !---------------------------------------------------
            !1-H (Error): Check for a Single O Atom (Error!!)
            if ( (H1 == 0) .and. (H2 == 0) .and. (H3 == 0) .and. (H4 == 0) ) then
               !
               write(*,'(2X, "NO Hydrogens within Rcut, Oxygen Number",I5)' ) i
               !
            !OH: Check for a Single OH molecule (Double PT)
            elseif ( (H1 /= 0) .and. (H2 == 0) .and. (H3 == 0) .and. (H4 == 0) ) then
               !
               write(*,'(2X, "Found: OH- molecule, Oxygen Number",I5)' ) i
               write(*,'(3X, 3(F10.7, 2X))') (tau(j,i,Ospecies),j=1,3)
               !
               !Set reference array
               Href(i,1:1) = (/ H1 /)
               !
            !H2O: Check for a Single H2O molecule 
            elseif ( (H1 /= 0) .and. (H2 /= 0) .and. (H3 == 0) .and. (H4 == 0) ) then
               !
               !Set reference array
               Href(i,1:2) = (/ H1, H2 /)
               !
            !H3Check for a singular H3O
            elseif ((H1 /= 0) .and. (H2 /= 0) .and. (H3 /= 0) .and. (H4 == 0)  ) then
               !
               write(*,'(2X, "Found: H3O+ molecule, Oxygen Number",I5)' ) i
               write(*,'(3X, 3(F10.7, 2X))') (tau(j,i,Ospecies),j=1,3)
               !
               !Set reference array
               Href(i,1:3) = (/ H1, H2, H3 /)
               !
            elseif ((H1 /= 0) .and. (H2 /= 0) .and. (H3 /= 0) .and. (H4 /= 0)) then 
               !
               write(*,*) ' ERROR: Four Hydrogens (or more) within Rcut :'
               !
            endif
            !
         enddo Oloop
         !======================================================================
         !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         !
         return 
         !
         End Subroutine find_Href
         !
      !Initialize the hbond variables
      Subroutine init_hbonds
         !
         implicit none
         !
         character(len=256)      :: buffer
         !
         !Main variabales
         filePos    = 'cp.pos'
         fileCel    = 'cp.cel'
         Ostar      = 0
         Ospecies   = 1
         Hspecies   = 2
         ntsp       = 0
         nsp(:)     = 0
         step       = 0
         rcut       = 1.1655
         !
         !hbonds variables
         Hbond_cut  = 3.5
         Hbond_angle= 150.0
         !
         !Read in the command line argument for the Ostar
         CALL getarg(1,buffer)

         if (buffer == '') then 
            write(*,*) "ERROR: Ostar is a command line arugment"
            stop
         else
            read(buffer,*) Ostar
         endif
         !
         !Read in parameters from standard input
         read(*, nml=system)
         !
         !Total number of atoms
         do ns=1,ntsp
         nat = nat + nsp(ns)
         enddo
         !
         !Set the cos_angle problem
         ref_cos = cos(Hbond_angle*(pi/180.0))
         !
         !Allocations/Initialization
         allocate( tau(3,(mic**3*nat),ntsp)       )  
         allocate( stau(3,(mic**3*nat),ntsp)      )
         allocate( Href(nsp(Ospecies),numHref)    )  
         allocate( sOshell(3, shell_max)          )
         allocate( O_index(shell_max)             )
         allocate( sHshell(3, shell_max,numHref)  )
         !
         Href(:,:) = 0
         sOshell(:,:) = 0.0
         sHshell(:,:,:) = 0.0
         !
         !
         return
      End Subroutine init_hbonds
      !
      !!
END Program Hbonds
