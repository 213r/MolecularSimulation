C     Compile keyword 
C     "f2py -c -m Potential_MM_f2py Potential_MM_f2py.f"

C     ****************************************************************** 
C     ****************************************************************** 
C     Lennard Jones Module 
C
C     ****************************************************************** 
C     Potential Energy and Interactive Force 
C      Between Ar and Ar

      subroutine get_lj_potential_arar(r, potential) 
        implicit none 
        real(8) :: eps, sigma, num6 
        real(8) :: r 
        real(8) :: potential 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: potential 
      
        eps = 0.000379386396928d0
        sigma = 6.43640672d0 
      
        num6 = (sigma/ r)**6 
        potential = 4.d0 * eps * num6 * (- 1.d0 + num6)
    
      end subroutine get_lj_potential_arar
    
      subroutine get_lj_force_arar(r, force) 
        implicit none 
        real(8) :: eps, sigma, num6  
        real(8) :: r
        real(8) :: force 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: force 
        
        eps = 0.000379386396928d0
        sigma = 6.43640672d0 
     
        num6 = (sigma/ r)**6 
        force = -24.d0 * eps * num6 * 
     &         (1.d0 - 2.d0 * num6) / r
        
      end subroutine get_lj_force_arar

C     ****************************************************************** 
C     Potential Energy and Interactive Force 
C      Between H and Ar

      subroutine get_lj_potential_har(r, potential) 
        implicit none 
        real(8) :: eps, sigma, num6 
        real(8) :: r 
        real(8) :: potential 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: potential 
      
        eps =  0.000602290708285d0   
        sigma = 3.23303770901d0
      
        num6 = (sigma/ r)**6 
        potential = 4.d0 * eps * num6 * (- 1.d0 + num6)
    
      end subroutine get_lj_potential_har
    
      subroutine get_lj_force_har(r, force) 
        implicit none 
        real(8) :: eps, sigma, num6  
        real(8) :: r
        real(8) :: force 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: force 
    
        eps =  0.000602290708285d0   
        sigma = 3.23303770901d0
     
        num6 = (sigma/ r)**6 
        force = -24.d0 * eps * num6 * 
     &         (1.d0 - 2.d0 * num6) / r
        
      end subroutine get_lj_force_har

C     ****************************************************************** 
C     Potential Energy and Interactive Force 
C      Between C and Ar

      subroutine get_lj_potential_car(r, potential) 
        implicit none 
        real(8) :: eps, sigma, num6 
        real(8) :: r 
        real(8) :: potential 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: potential 
      
        eps = 0.00107403818804d0   
        sigma = 3.29946157751d0
      
        num6 = (sigma/ r)**6 
        potential = 4.d0 * eps * num6 * (- 1.d0 + num6)
    
      end subroutine get_lj_potential_car
    
      subroutine get_lj_force_car(r, force) 
        implicit none 
        real(8) :: eps, sigma, num6  
        real(8) :: r
        real(8) :: force 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: force 
        
        eps = 0.00107403818804d0   
        sigma = 3.29946157751d0
     
        num6 = (sigma/ r)**6 
        force = -24.d0 * eps * num6 * 
     &         (1.d0 - 2.d0 * num6) / r
        
      end subroutine get_lj_force_car

C     ****************************************************************** 
C     Potential Energy and Interactive Force 
C      Between O and Ar
 
      subroutine get_lj_potential_oar(r, potential) 
        implicit none 
        real(8) :: eps, sigma, num6 
        real(8) :: r 
        real(8) :: potential 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: potential 
      
        eps =  0.00100217023372d0
        sigma = 3.4166245888d0
      
        num6 = (sigma/ r)**6 
        potential = 4.d0 * eps * num6 * (- 1.d0 + num6)
    
      end subroutine get_lj_potential_oar
    
      subroutine get_lj_force_oar(r, force) 
        implicit none 
        real(8) :: eps, sigma, num6  
        real(8) :: r
        real(8) :: force 
Cf2py   intent(in) :: r
Cf2py   intent(out) :: force 
        
        eps =  0.00100217023372d0
        sigma = 3.4166245888d0
     
        num6 = (sigma/ r)**6 
        force = -24.d0 * eps * num6 * 
     &         (1.d0 - 2.d0 * num6) / r
        
      end subroutine get_lj_force_oar

C     ****************************************************************** 
C     ****************************************************************** 
      
      subroutine potential_mm_f2py(energy, force, natom,
     &           ndim, atom_numbers, positions, rlimit, 
     &           check_pbc, box, invbox)
       
        implicit none 
        integer natom, ndim
        integer atom_numbers(natom)
        logical check_pbc
        real(8) positions(ndim, natom)
        real(8) box(ndim, ndim), invbox(ndim, ndim)
        real(8) s(ndim), n(ndim) 
        real(8) force(ndim, natom) 
        real(8) rlimit, energy
Cf2py intent(in) natom, ndim, atom_numbers, positions, rlimit
Cf2py intent(in) check_pbc, box, invbox
Cf2py intent(out) energy, force

        integer :: i, j, k, numi, numj, tmp 
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
            numi = atom_numbers(i)
            numj = atom_numbers(j)
            if (numi > numj) then 
              tmp = numj
              numj = numi 
              numi = tmp
            end if 
            if (numi == 1 .and. numj == 18) then 
              call get_lj_potential_har(rji, pot_lj) 
              call get_lj_force_har(rji, force_lj) 
            elseif (numi == 6 .and. numj == 18) then 
              call get_lj_potential_car(rji, pot_lj) 
              call get_lj_force_car(rji, force_lj) 
            elseif (numi == 8 .and. numj == 18) then 
              call get_lj_potential_oar(rji, pot_lj) 
              call get_lj_force_oar(rji, force_lj) 
            elseif (numi == 18 .and. numj == 18) then 
              call get_lj_potential_arar(rji, pot_lj) 
              call get_lj_force_arar(rji, force_lj) 
            else
              write(*,*) "Parameters of these atom pairs haven't 
     &                    been set yet" 
            end if 
            energy = energy + pot_lj 
            force_ji(:) = vec_ji(:) * (force_lj / rji) 
            force(:,i) = force(:,i) + force_ji(:) 
            force(:,j) = force(:,j) - force_ji(:) 
          end do 
      
        end do 
      
      end subroutine potential_mm_f2py

C     ****************************************************************** 
C     ****************************************************************** 
      
      subroutine potential_mm_f2py_mc(energy, natom,
     &           ndim, ind, atom_numbers, positions, rlimit,
     &           check_pbc, box, invbox)
       
        implicit none 
        integer natom, ndim, ind
        integer atom_numbers(natom)
        logical check_pbc
        real(8) positions(ndim, natom)
        real(8) box(ndim, ndim), invbox(ndim, ndim)
        real(8) s(ndim), n(ndim) 
        real(8) energy(natom)
        real(8) rlimit
Cf2py intent(in) natom, ndim, ind, atom_numbers, positions, rlimit
Cf2py intent(in) check_pbc, box, invbox
Cf2py intent(out) energy

        integer :: i, k, numi, numj, tmp 
        real(8) :: rji, pot_lj  
        real(8), dimension(ndim) :: vec_ji 
      
        if ((ndim .ne. 3) .and. check_pbc) then
           write(*,*) 'pbc is valid only for 3-dimension'
           stop  
        endif 

        energy(:) = 0.d0
        ind = ind + 1
        do i = 1, natom  
            numj = atom_numbers(ind)
            if (i .eq. ind) cycle  
            rji = 0.d0  
            vec_ji(:) = positions(:,i) - positions(:,ind)
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
            numi = atom_numbers(i)
            if (numi > numj) then 
              tmp = numj
              numj = numi 
              numi = tmp
            end if 
            if (numi == 1 .and. numj == 18) then 
              call get_lj_potential_har(rji, pot_lj) 
            elseif (numi == 6 .and. numj == 18) then 
              call get_lj_potential_car(rji, pot_lj) 
            elseif (numi == 8 .and. numj == 18) then 
              call get_lj_potential_oar(rji, pot_lj) 
            elseif (numi == 18 .and. numj == 18) then 
              call get_lj_potential_arar(rji, pot_lj) 
            else
              pot_lj = 0.d0 
            end if 
            energy(i) = pot_lj 
      
        end do 
      
      end subroutine potential_mm_f2py_mc

C     ****************************************************************** 
C     ****************************************************************** 
C     ****************************************************************** 

      subroutine potential_qmmm_f2py(energy, force1, force2,
     &           natom1, natom2, ndim, atom_numbers1,
     &           atom_numbers2, positions1, positions2, rlimit,
     &           check_pbc, box, invbox)
       
        implicit none 
        integer natom1, natom2, ndim
        integer atom_numbers1(natom1), atom_numbers2(natom2) 
        logical check_pbc
        real(8) positions1(ndim, natom1), positions2(ndim, natom2)
        real(8) force1(ndim, natom1), force2(ndim, natom2) 
        real(8) box(ndim, ndim), invbox(ndim, ndim)
        real(8) s(ndim), n(ndim) 
        real(8) rlimit, energy
Cf2py intent(in) natom1, natom2 ndim, atom_numbers1, atom_numbers2, positions,rlimit
Cf2py intent(in) positions1, positions2, rlimit
Cf2py intent(in) check_pbc, box, invbox
Cf2py intent(out) energy, force1, force2

        integer :: i, j, k, numi, numj, tmp 
        real(8) :: rji, pot_lj, force_lj 
        real(8), dimension(ndim) :: vec_ji, force_ji 
        
        if ((ndim .ne. 3) .and. check_pbc) then
           write(*,*) 'pbc is valid only for 3-dimension'
           stop  
        endif 
      
        energy = 0.d0
        force1(:,:) = 0.d0
        force2(:,:) = 0.d0
        do i = 1, natom1  
          do j = 1, natom2  
            rji = 0.d0  
            vec_ji(:) = positions1(:,i) - positions2(:,j)
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
            numi = atom_numbers1(i)
            numj = atom_numbers2(j)
            if (numi > numj) then 
              tmp = numj
              numj = numi 
              numi = tmp
            end if 
            if (numi == 1 .and. numj == 18) then 
              call get_lj_potential_har(rji, pot_lj) 
              call get_lj_force_har(rji, force_lj) 
            elseif (numi == 6 .and. numj == 18) then 
              call get_lj_potential_car(rji, pot_lj) 
              call get_lj_force_car(rji, force_lj) 
            elseif (numi == 8 .and. numj == 18) then 
              call get_lj_potential_oar(rji, pot_lj) 
              call get_lj_force_oar(rji, force_lj) 
            elseif (numi == 18 .and. numj == 18) then 
              call get_lj_potential_arar(rji, pot_lj) 
              call get_lj_force_arar(rji, force_lj) 
            else
              write(*,*) "Parameters of these atom pairs haven't 
     &                    been set yet" 
            end if 
            energy = energy + pot_lj 
            force_ji(:) = vec_ji(:) * (force_lj / rji) 
            force1(:,i) = force1(:,i) + force_ji(:) 
            force2(:,j) = force2(:,j) - force_ji(:) 
          end do 
      
        end do 
      
      end subroutine potential_qmmm_f2py


