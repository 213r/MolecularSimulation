C     Compile keyword 
C     "f2py -c --fcompiler=gfortran -m Potential_MM_f2py Potential_MM_f2py.f"

C     ****************************************************************** 
C     ****************************************************************** 
C     Lennard Jones Module 
C
C     ****************************************************************** 
C     Potential Energy and Interactive Force 
C      Between Ar and Ar

      subroutine get_lj_potential_arar(r, potential)
         implicit none
         real(8) :: eps, sigma
         real(8) :: r, m
         real(8) :: potential
Cf2py   intent(in) :: r
Cf2py   intent(out) :: potential 
      
        eps = 0.000375861d0  
        sigma = 7.180958757d0 
        m = 0.88d0 
        potential=eps*((1.d0-exp(-m * (r - sigma)))**2-1.d0)
    
      end subroutine get_lj_potential_arar
    
      subroutine get_lj_force_arar(r, force) 
        implicit none 
        real(8) :: eps, sigma  
        real(8) :: r, m, ex
        real(8) :: force 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: force 
 
        eps = 0.000375861d0  
        sigma = 7.180958757d0 
        m = 0.88d0 
  
        ex = exp(-m * (r - sigma))
        force = -2.d0 * m * eps * ex * (1.d0 - ex)
        
      end subroutine get_lj_force_arar

C     ****************************************************************** 
C     Potential Energy and Interactive Force 
C      Between H and Ar
     
      subroutine get_lj_potential_har_adjusted(r, r_hc, potential) 
        implicit none 
        real(8) :: eps0, sigma0, eps1,sigma1,m0, m1
        real(8) :: pot0,pot1 
        real(8) :: r, r_hc, adjust  
        real(8) :: potential 
Cf2py   intent(in) :: r, r_hc
Cf2py   intent(out) :: potential 
      
C       f(x) = 1 / (1 + exp(-7.0*(x-2.0))) at the unit of ang
C       the value of f(x) changes drastically in the range of [1, 3] ang  
        adjust = 1.d0 / (1.d0 + exp(-7.0d0 * (r_hc*0.529d0 - 2.0d0)))  

C       the H1(-c2) of HCOOH 
        eps0 = 0.000845957d0 
        sigma0 = 5.619710956d0
        m0 = 0.85d0 
        pot0=eps0*((1.d0 - exp(-m0 * (r - sigma0)))**2 - 1.d0)

C       the sole H
        eps1 = 0.00015622d0  
        sigma1 = 6.823643849d0 
        m1 = 0.77d0 
        pot1=eps1*((1.d0 - exp(-m1 * (r - sigma1)))**2 - 1.d0)
      
        potential = pot0 + (pot1 - pot0) * adjust
    
      end subroutine get_lj_potential_har_adjusted

      subroutine get_lj_force_har_adjusted(r, r_hc, force) 
        implicit none 
        real(8) :: eps0, sigma0, eps1,sigma1,m0, m1
        real(8) :: force0, force1, ex0, ex1
        real(8) :: r, r_hc, adjust  
        real(8) :: force 
Cf2py   intent(in) :: r, r_hc
Cf2py   intent(out) :: force 
      
C       f(x) = 1 / (1 + exp(-7.0*(x-2.0))) at the unit of ang
C       the value of f(x) changes drastically in the range of [1, 3] ang  
        adjust = 1.d0 / (1.d0 + exp(-7.0d0 * (r_hc*0.529d0 - 2.0d0)))  

C       the H1(-c2) of HCOOH 
        eps0 = 0.000845957d0 
        sigma0 = 5.619710956d0
        m0 = 0.85d0 
        ex0 = exp(-m0 * (r - sigma0))
        force0 = -2.d0 * m0 * eps0 * ex0 * (1.d0 - ex0)

C       the sole H
        eps1 = 0.00015622d0  
        sigma1 = 6.823643849d0 
        m1 = 0.77d0 
        ex1 = exp(-m1 * (r - sigma1))
        force1 = -2.d0 * m1 * eps1 * ex1 * (1.d0 - ex1)
      
        force = force0 + (force1 - force0) * adjust
    
      end subroutine get_lj_force_har_adjusted

C     ****************************************************************** 
C     Potential Energy and Interactive Force 
C      Between C and Ar

      subroutine get_lj_potential_car(r, potential) 
        implicit none 
        real(8) :: eps, sigma
        real(8) :: r, m  
        real(8) :: potential 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: potential 
      
        eps = 0.0010312d0  
        sigma = 6.653570943d0 
        m = 0.78d0 

       potential = eps * ((1.d0 - exp(-m * (r - sigma)))**2 - 1.d0)
    
      end subroutine get_lj_potential_car
    
      subroutine get_lj_force_car(r, force) 
        implicit none 
        real(8) :: eps, sigma, ex  
        real(8) :: r, m
        real(8) :: force 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: force 
        eps = 0.0010312d0  
        sigma = 6.653570943d0 
        m = 0.78d0 
        
        ex = exp(-m * (r - sigma))
        force = -2.d0 * m * eps * ex * (1.d0 - ex)
        
      end subroutine get_lj_force_car

C     ****************************************************************** 
C     Potential Energy and Interactive Force 
C      Between O and Ar
 
      subroutine get_lj_potential_oar(r, potential) 
        implicit none 
        real(8) :: eps, sigma
        real(8) :: r, m  
        real(8) :: potential 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: potential 
      
C      Between O3 and Ar
        eps = 0.000751319d0  
        sigma = 6.482703137d0 
        m = 0.86d0 

       potential = eps * ((1.d0 - exp(-m * (r - sigma)))**2 - 1.d0)
    
      end subroutine get_lj_potential_oar
    
      subroutine get_lj_force_oar(r, force) 
        implicit none 
        real(8) :: eps, sigma, ex  
        real(8) :: r, m
        real(8) :: force 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: force 
 
C      Between O3 and Ar
        eps = 0.000751319d0  
        sigma = 6.482703137d0 
        m = 0.86d0 

        ex = exp(-m * (r - sigma))
        force = -2.d0 * m * eps * ex * (1.d0 - ex)
        
      end subroutine get_lj_force_oar
    
C     ****************************************************************** 
C     this subroutine give both energy and forces in MM domain of Ar 
      
      subroutine potential_mm_f2py_ar(energy, force, natom,
     &           ndim, positions, rlimit, 
     &           check_pbc, box, invbox)
       
        implicit none 
        integer ndim 
        integer natom
        logical check_pbc 
        real(8) positions(ndim, natom)
        real(8) box(ndim, ndim), invbox(ndim, ndim)
        real(8) s(ndim), n(ndim) 
        real(8) force(ndim, natom) 
        real(8) rlimit, energy
Cf2py intent(in) natom, atom_numbers, positions, rlimit
Cf2py intent(in) ndim, check_pbc, box, invbox
Cf2py intent(out) energy, force

        integer :: i, j, k 
        real(8) :: rji, pot_lj, force_lj 
        real(8), dimension(ndim) :: vec_ji, force_ji 
      
        if ((ndim .ne. 3) .and. check_pbc) then
           write(*,*) 'pbc is valid only for 3-dimension'
           stop  
        endif 
        
        energy = 0.d0
        force(:,:) = 0.d0
        
        do i = 1, natom-1  
          do j = i+1, natom  
            
            rji = 0.d0  
            pot_lj = 0.d0 
            force_lj = 0.d0 
            vec_ji(:) = positions(:,i) - positions(:,j)
            if (check_pbc) then 
                s(:) = matmul(invbox(:,:), vec_ji(:)) 
                n(:) = - nint(s(:))
                s(:) = s(:) + n(:)
                vec_ji(:) = matmul(box(:,:),s(:))
            end if  
            do k = 1, ndim  
              rji = rji + vec_ji(k) ** 2
            end do 
            rji = sqrt(rji) 
            if (rji > rlimit) cycle 
           
            call get_lj_potential_arar(rji, pot_lj) 
            call get_lj_force_arar(rji, force_lj) 
            
            energy = energy + pot_lj 
C            force_lj = - force_lj 
            
            force_ji(:) = vec_ji(:) * (force_lj / rji) 
            force(:,i) = force(:,i) + force_ji(:) 
            force(:,j) = force(:,j) - force_ji(:) 
          end do 
      
        end do 
      
      end subroutine potential_mm_f2py_ar

C     ****************************************************************** 
C     ****************************************************************** 
C     this subroutine give both energy and forces between each atom of 
C     HCOOH and Ar in MM domain 

      subroutine potential_qmmm_f2py_hcooh_ar(energy_qmmm, force_hcooh, 
     &           force_ar,
     &           ndim,
     &           natom_hcooh,
     &           natom_ar,
     &           positions_hcooh, 
     &           positions_ar, 
     &           rlimit, check_pbc, box, invbox)
       
        implicit none 
        integer ndim, natom_hcooh, natom_ar
        logical check_pbc
        real(8) positions_hcooh(ndim,natom_hcooh)
        real(8) positions_ar(ndim, natom_ar)
        real(8) force_hcooh(ndim, natom_hcooh), force_ar(ndim, natom_ar)
        real(8) box(ndim, ndim), invbox(ndim, ndim)
        real(8) s(ndim), n(ndim) 
        real(8) rlimit, energy_qmmm
Cf2py intent(in) natom_ar, positions_ar,rlimit
Cf2py intent(in) natom_hcooh, positions_hcooh, 
Cf2py intent(in) check_pbc, box, invbox
Cf2py intent(out) energy_qmmm, force_hcooh, force_ar

        integer :: i, j, k 
        real(8) :: rji, r_h2, r_h5, pot_lj, force_lj 
        real(8) :: r_h2_save, r_h5_save 
        real(8), dimension(ndim) :: vec_ji, force_ji 
        real(8), dimension(ndim) :: vec_h2, vec_h5
        
        if ((ndim .ne. 3) .and. check_pbc) then
           write(*,*) 'pbc is valid only for 3-dimension'
           stop  
        endif 
       
        energy_qmmm = 0.d0
        force_hcooh(:,:) = 0.d0
        force_ar(:,:) = 0.d0
        r_h2_save = 1d5 
        r_h5_save = 1d5 

        do i = 1, natom_hcooh 
           r_h2 = 0.d0  
           if (i == 2) cycle 
           vec_h2(:) = positions_hcooh(:,2) - positions_hcooh(:,i)
           do k = 1, ndim  
              r_h2 = r_h2 + vec_h2(k) ** 2
           end do 
           r_h2_save = min(r_h2_save, sqrt(r_h2)) 
        end do

        do i = 1, natom_hcooh 
           r_h5 = 0.d0  
           if (i == 5) cycle 
           vec_h5(:) = positions_hcooh(:,5) - positions_hcooh(:,i)
           do k = 1, ndim  
              r_h5 = r_h5 + vec_h5(k) ** 2
           end do 
           r_h5_save = min(r_h5_save, sqrt(r_h5)) 
        end do 
           
        do i = 1, natom_hcooh 
           do j = 1, natom_ar  
            rji = 0.d0  
            pot_lj = 0.d0  
            force_lj = 0.d0  
            vec_ji(:) = positions_ar(:,j) - positions_hcooh(:,i)
            if (check_pbc) then 
                s(:) = matmul(invbox(:,:), vec_ji(:)) 
                n(:) = - nint(s(:))
                s(:) = s(:) + n(:)
                vec_ji(:) = matmul(box(:,:),s(:))
            end if  
            do k = 1, ndim  
              rji = rji + vec_ji(k) ** 2
            end do 
            rji = sqrt(rji) 
C           write(*,*) i,j,rji 
            if (rji > rlimit) cycle 
            
            if (i == 1) then 
              call get_lj_potential_car(rji, pot_lj) 
              call get_lj_force_car(rji, force_lj) 
            elseif (i == 2) then 
              call get_lj_potential_har_adjusted(rji, r_h2_save, pot_lj)
              call get_lj_force_har_adjusted(rji, r_h2_save, force_lj) 
            elseif (i == 3) then 
              call get_lj_potential_oar(rji, pot_lj) 
              call get_lj_force_oar(rji, force_lj) 
            elseif (i == 4) then 
              call get_lj_potential_oar(rji, pot_lj) 
              call get_lj_force_oar(rji, force_lj) 
            elseif (i == 5) then 
              call get_lj_potential_har_adjusted(rji, r_h5_save, pot_lj)
              call get_lj_force_har_adjusted(rji, r_h5_save, force_lj) 
            else
              write(*,*) "Parameters of these atom pairs haven't 
     &                    been set yet" 
            end if 
            force_lj = -force_lj 
            energy_qmmm = energy_qmmm + pot_lj 
C            write(*,*) i,j, force_lj 
            force_ji(:) = vec_ji(:) * (force_lj / rji) 
            force_hcooh(:,i) = force_hcooh(:,i) + force_ji(:) 
            force_ar(:,j) = force_ar(:,j) - force_ji(:) 
          end do 
      
        end do 
      
      end subroutine potential_qmmm_f2py_hcooh_ar

C     ****************************************************************** 
C     ****************************************************************** 
C         These below methods are only called for MCMC calculation
C     ****************************************************************** 
C     ****************************************************************** 

      subroutine potential_mm_f2py_ar_mc(energy, ind, natom,
     &           ndim, positions, rlimit, 
     &           check_pbc, box, invbox)
       
        implicit none 
        integer ind 
        integer ndim 
        integer natom
        logical check_pbc 
        real(8) positions(ndim, natom)
        real(8) box(ndim, ndim), invbox(ndim, ndim)
        real(8) s(ndim), n(ndim) 
        real(8) energy(natom)
        real(8) rlimit
Cf2py intent(in) natom, atom_numbers, positions, rlimit
Cf2py intent(in) ind, ndim, check_pbc, box, invbox
Cf2py intent(out) energy

        integer :: i, j, k 
        real(8) :: rji, pot_lj  
        real(8), dimension(ndim) :: vec_ji 
      
        if ((ndim .ne. 3) .and. check_pbc) then
           write(*,*) 'pbc is valid only for 3-dimension'
           stop  
        endif 
        
        energy(:) = 0.d0
       
        j = ind + 1

        do i = 1, natom  
          if (i .eq. j) cycle   
            
            rji = 0.d0  
            pot_lj = 0.d0 
            vec_ji(:) = positions(:,i) - positions(:,j)
            if (check_pbc) then 
                s(:) = matmul(invbox(:,:), vec_ji(:)) 
                n(:) = - nint(s(:))
                s(:) = s(:) + n(:)
                vec_ji(:) = matmul(box(:,:),s(:))
            end if  
            do k = 1, ndim  
              rji = rji + vec_ji(k) ** 2
            end do 
            rji = sqrt(rji) 
            if (rji > rlimit) cycle 
           
            call get_lj_potential_arar(rji, pot_lj) 
            energy(i) = pot_lj 
      
        end do 
      
      end subroutine potential_mm_f2py_ar_mc

C     ****************************************************************** 
C     ****************************************************************** 

      subroutine potential_qmmm_f2py_hcooh_ar_mc(energy_qmmm, 
     &           ind,
     &           ndim,
     &           natom_hcooh,
     &           natom_ar,
     &           positions_hcooh, 
     &           positions_ar, 
     &           rlimit, check_pbc, box, invbox)
       
        implicit none 
        integer ind, ndim, natom_hcooh, natom_ar
        logical check_pbc
        real(8) positions_hcooh(ndim,natom_hcooh)
        real(8) positions_ar(ndim, natom_ar)
        real(8) box(ndim, ndim), invbox(ndim, ndim)
        real(8) s(ndim), n(ndim) 
        real(8) rlimit, energy_qmmm
Cf2py intent(in) natom_ar, positions_ar,rlimit
Cf2py intent(in) natom_hcooh, positions_hcooh, 
Cf2py intent(in) check_pbc, box, invbox
Cf2py intent(out) energy_qmmm

        integer :: i, j, k 
        real(8) :: rji, r_h2, r_h5, pot_lj 
        real(8) :: r_h2_save, r_h5_save 
        real(8), dimension(ndim) :: vec_ji 
        real(8), dimension(ndim) :: vec_h2, vec_h5
        
        if ((ndim .ne. 3) .and. check_pbc) then
           write(*,*) 'pbc is valid only for 3-dimension'
           stop  
        endif 
       
        energy_qmmm = 0.d0
        r_h2_save = 1d5 
        r_h5_save = 1d5 

        do i = 1, natom_hcooh 
           r_h2 = 0.d0  
           if (i == 2) cycle 
           vec_h2(:) = positions_hcooh(:,2) - positions_hcooh(:,i)
           do k = 1, ndim  
              r_h2 = r_h2 + vec_h2(k) ** 2
           end do 
           r_h2_save = min(r_h2_save, sqrt(r_h2)) 
        end do

        do i = 1, natom_hcooh 
           r_h5 = 0.d0  
           if (i == 5) cycle 
           vec_h5(:) = positions_hcooh(:,5) - positions_hcooh(:,i)
           do k = 1, ndim  
              r_h5 = r_h5 + vec_h5(k) ** 2
           end do 
           r_h5_save = min(r_h5_save, sqrt(r_h5)) 
        end do 
    
        j = ind + 1

        do i = 1, natom_hcooh  
            rji = 0.d0  
            pot_lj = 0.d0  
            vec_ji(:) = positions_ar(:,j) - positions_hcooh(:,i)
            if (check_pbc) then 
                s(:) = matmul(invbox(:,:), vec_ji(:)) 
                n(:) = - nint(s(:))
                s(:) = s(:) + n(:)
                vec_ji(:) = matmul(box(:,:),s(:))
            end if  
            do k = 1, ndim  
              rji = rji + vec_ji(k) ** 2
            end do 
            rji = sqrt(rji) 
C           write(*,*) i,j,rji 
            if (rji > rlimit) cycle 
            
            if (i == 1) then 
              call get_lj_potential_car(rji, pot_lj) 
            elseif (i == 2) then 
              call get_lj_potential_har_adjusted(rji, r_h2_save, pot_lj)
            elseif (i == 3) then 
              call get_lj_potential_oar(rji, pot_lj) 
            elseif (i == 4) then 
              call get_lj_potential_oar(rji, pot_lj) 
            elseif (i == 5) then 
              call get_lj_potential_har_adjusted(rji, r_h5_save, pot_lj)
            else
              write(*,*) "Parameters of these atom pairs haven't 
     &                    been set yet" 
            end if 
            energy_qmmm = energy_qmmm + pot_lj 
C            write(*,*) i,j, force_lj 

          end do 
      
      
      end subroutine potential_qmmm_f2py_hcooh_ar_mc


