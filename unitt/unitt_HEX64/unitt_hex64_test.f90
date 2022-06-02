program unitt_hex64_test

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Tests how well the element type HEX27 is implemented through     !
   ! geometrical tests over the reference element.                    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use mod_constants     ! Local version of mod_constants
   use elem_hex          ! Module with all basic element definitions
   use mod_maths         ! Module with all mathematical functions
   use quadrature_rules  ! Module with all quadrature rules
   use jacobian_oper     ! Module with all jacobian operators

   implicit none

   integer(4) :: idime, inode, igaus
   integer(4) :: jdime, jnode, jgaus
   integer(4) :: atoIJK(ngaus), vtk_atoIJK(ngaus), listHEX08(27,8)
   real(8)    :: s, t, z, aux
   real(8)    :: elcod(nnode,ndime)
   real(8)    :: xgp(ngaus,ndime), wgp(ngaus)
   real(8)    :: Nl(nnode), dNl(ndime,nnode)
   real(8)    :: N(ngaus,nnode), dN(ndime,nnode,ngaus)
   logical    :: flag_pass

   !
   ! Initialize check
   !
   flag_pass = .true. ! Will flip to false if a test fails

   !
   ! Element geometry: comes from the chebyshev_hex quadrature rule
   !
   call hex64(1.0d0,1.0d0,1.0d0,atoIJK,vtk_atoIJK,listHEX08)
   call chebyshev_hex(atoIJK,xgp,wgp)

   !
   ! Compute N_a(xi_b); check whether it is 1 only when a == b
   !
   outer: do igaus = 1,ngaus
      s = xgp(igaus,1)
      t = xgp(igaus,2)
      z = xgp(igaus,3)
      call hex64(s,t,z,atoIJK,vtk_atoIJK,listHEX08,N(igaus,:),dN(:,:,igaus),Nl,dNl)
      do jgaus = 1,ngaus
         aux = N(igaus,jgaus)
         if (jgaus == igaus .and. aux .ne. 1.0d0) then
            flag_pass = .false.
            exit outer
         else if (jgaus .ne. igaus .and. aux .ne. 0.0d0) then
            flag_pass = .false.
            exit outer
         end if
      end do
   end do outer

   !
   ! Check if test failed
   !
   if (flag_pass) then
      write(*,*) 'HEX64: passed'
      STOP
   else
      write(*,*) 'HEX64: failed'
      ERROR STOP
   end if

end program unitt_hex64_test
