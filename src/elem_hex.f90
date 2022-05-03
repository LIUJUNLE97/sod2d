module elem_hex

   use mod_constants
   use mod_maths

   contains  

      subroutine hex08(xi,eta,zeta,N,dN) ! HEX08 element     

         implicit none
         
         real(8), intent(in)  :: xi, eta, zeta
         real(8), intent(out) :: N(8), dN(3,8)

         ! Regarding coordinates for the isopar space:
         ! xi is left to right (-1,1)
         ! eta is bottom to top (-1,1)
         ! zeta is back to front (-1,1)
         ! These are like the Gmsh (x,y,z) system.
         ! These will apply to all orders of HEX

         N(1) = (1.0d0-xi)*(1.0d0-eta)*(1.0d0-zeta)
         N(2) = (1.0d0-xi)*(1.0d0-eta)*(1.0d0+zeta)
         N(3) = (1.0d0+xi)*(1.0d0-eta)*(1.0d0+zeta)
         N(4) = (1.0d0+xi)*(1.0d0-eta)*(1.0d0-zeta)
         N(5) = (1.0d0-xi)*(1.0d0+eta)*(1.0d0-zeta)
         N(6) = (1.0d0-xi)*(1.0d0+eta)*(1.0d0+zeta)
         N(7) = (1.0d0+xi)*(1.0d0+eta)*(1.0d0+zeta)
         N(8) = (1.0d0+xi)*(1.0d0+eta)*(1.0d0-zeta)
         N = 0.125d0*N
         
         dN(1,1) = -(1.0d0-eta)*(1.0d0-zeta)
         dN(2,1) = -(1.0d0-xi )*(1.0d0-zeta)
         dN(3,1) = -(1.0d0-xi )*(1.0d0-eta )

         dN(1,2) = -(1.0d0-eta)*(1.0d0+zeta)
         dN(2,2) = -(1.0d0-xi )*(1.0d0+zeta)
         dN(3,2) =  (1.0d0-xi )*(1.0d0-eta )

         dN(1,3) =  (1.0d0-eta)*(1.0d0+zeta)
         dN(2,3) = -(1.0d0+xi )*(1.0d0+zeta)
         dN(3,3) =  (1.0d0+xi )*(1.0d0-eta )

         dN(1,4) =  (1.0d0-eta)*(1.0d0-zeta)
         dN(2,4) = -(1.0d0+xi )*(1.0d0-zeta)
         dN(3,4) = -(1.0d0+xi )*(1.0d0-eta )

         dN(1,5) = -(1.0d0+eta)*(1.0d0-zeta)
         dN(2,5) =  (1.0d0-xi )*(1.0d0-zeta)
         dN(3,5) = -(1.0d0-xi )*(1.0d0+eta )

         dN(1,6) = -(1.0d0+eta)*(1.0d0+zeta)
         dN(2,6) =  (1.0d0-xi )*(1.0d0+zeta)
         dN(3,6) =  (1.0d0-xi )*(1.0d0+eta )

         dN(1,7) =  (1.0d0+eta)*(1.0d0+zeta)
         dN(2,7) =  (1.0d0+xi )*(1.0d0+zeta)
         dN(3,7) =  (1.0d0+xi )*(1.0d0+eta )

         dN(1,8) =  (1.0d0+eta)*(1.0d0-zeta)
         dN(2,8) =  (1.0d0+xi )*(1.0d0-zeta)
         dN(3,8) = -(1.0d0+xi )*(1.0d0+eta )
         dN = 0.125d0*dN

      end subroutine hex08
      
      subroutine hex27(xi,eta,zeta,N,dN) ! HEX27 element

         implicit none
         
         real(8), intent(in)  :: xi, eta, zeta
         real(8), intent(out) :: N(27), dN(3,27)
         real(8)              :: n1xi, n2xi, n3xi
         real(8)              :: n1eta, n2eta, n3eta
         real(8)              :: n1zeta, n2zeta, n3zeta
         real(8)              :: dn1xi, dn2xi, dn3xi
         real(8)              :: dn1eta, dn2eta, dn3eta
         real(8)              :: dn1zeta, dn2zeta, dn3zeta

         !
         ! Triquadratic brick
         !        

         ! 1D quadratic functions for xi, eta, zeta
         n1xi   = 0.5d0*xi*(xi-1.0d0)
         n1eta  = 0.5d0*eta*(eta-1.0d0)
         n1zeta = 0.5d0*zeta*(zeta-1.0d0)
         n2xi   = 1.0d0-xi*xi
         n2eta  = 1.0d0-eta*eta
         n2zeta = 1.0d0-zeta*zeta
         n3xi   = 0.5d0*xi*(xi+1.0d0)
         n3eta  = 0.5d0*eta*(eta+1.0d0)
         n3zeta = 0.5d0*zeta*(zeta+1.0d0)

         ! 1D quadratic function derivatives for xi, eta, zeta
         dn1xi   = 0.5d0*(2.0d0*xi-1.0d0)
         dn1eta  = 0.5d0*(2.0d0*eta-1.0d0)
         dn1zeta = 0.5d0*(2.0d0*zeta-1.0d0)
         dn2xi   = -2.0d0*xi
         dn2eta  = -2.0d0*eta
         dn2zeta = -2.0d0*zeta
         dn3xi   = 0.5d0*(2.0d0*xi+1.0d0)
         dn3eta  = 0.5d0*(2.0d0*eta+1.0d0)
         dn3zeta = 0.5d0*(2.0d0*zeta+1.0d0)

         !
         ! Shape functions
         !

         ! Bottom, corner
         N( 1) = n1xi*n1eta*n1zeta
         N( 2) = n1xi*n1eta*n3zeta
         N( 3) = n3xi*n1eta*n3zeta
         N( 4) = n3xi*n1eta*n1zeta

         ! Top, corner
         N( 5) = n1xi*n3eta*n1zeta
         N( 6) = n1xi*n3eta*n3zeta
         N( 7) = n3xi*n3eta*n3zeta
         N( 8) = n3xi*n3eta*n1zeta

         ! Bottom, edges
         N( 9) = n1xi*n1eta*n2zeta
         N(10) = n2xi*n1eta*n3zeta
         N(11) = n3xi*n1eta*n2zeta
         N(12) = n2xi*n1eta*n1zeta

         ! Top, edges
         N(13) = n1xi*n3eta*n2zeta
         N(14) = n2xi*n3eta*n3zeta
         N(15) = n3xi*n3eta*n2zeta
         N(16) = n2xi*n3eta*n1zeta

         ! Middle, edges
         N(17) = n1xi*n2eta*n1zeta
         N(18) = n1xi*n2eta*n3zeta
         N(19) = n3xi*n2eta*n3zeta
         N(20) = n3xi*n2eta*n1zeta

         ! Faces
         N(21) = n2xi*n2eta*n1zeta
         N(22) = n2xi*n2eta*n3zeta
         N(23) = n1xi*n2eta*n2zeta
         N(24) = n3xi*n2eta*n2zeta
         N(25) = n2xi*n1eta*n2zeta
         N(26) = n2xi*n3eta*n2zeta

         ! Volume
         N(27) = n2xi*n2eta*n2zeta

         !
         ! Shape function derivatives
         !

         ! Direction xi

         ! Bottom, corner
         dN(1, 1) = dn1xi*n1eta*n1zeta
         dN(1, 2) = dn1xi*n1eta*n3zeta
         dN(1, 3) = dn3xi*n1eta*n3zeta
         dN(1, 4) = dn3xi*n1eta*n1zeta

         ! Top, corner
         dN(1, 5) = dn1xi*n3eta*n1zeta
         dN(1, 6) = dn1xi*n3eta*n3zeta
         dN(1, 7) = dn3xi*n3eta*n3zeta
         dN(1, 8) = dn3xi*n3eta*n1zeta

         ! Bottom, edges
         dN(1, 9) = dn1xi*n1eta*n2zeta
         dN(1,10) = dn2xi*n1eta*n3zeta
         dN(1,11) = dn3xi*n1eta*n2zeta
         dN(1,12) = dn2xi*n1eta*n1zeta

         ! Top, edges
         dN(1,13) = dn1xi*n3eta*n2zeta
         dN(1,14) = dn2xi*n3eta*n3zeta
         dN(1,15) = dn3xi*n3eta*n2zeta
         dN(1,16) = dn2xi*n3eta*n1zeta

         ! Middle, edges
         dN(1,17) = dn1xi*n2eta*n1zeta
         dN(1,18) = dn1xi*n2eta*n3zeta
         dN(1,19) = dn3xi*n2eta*n3zeta
         dN(1,20) = dn3xi*n2eta*n1zeta

         ! Faces
         dN(1,21) = dn2xi*n2eta*n1zeta
         dN(1,22) = dn2xi*n2eta*n3zeta
         dN(1,23) = dn1xi*n2eta*n2zeta
         dN(1,24) = dn3xi*n2eta*n2zeta
         dN(1,25) = dn2xi*n1eta*n2zeta
         dN(1,26) = dn2xi*n3eta*n2zeta

         ! Volume
         dN(1,27) = dn2xi*n2eta*n2zeta

         ! Direction eta

         ! Bottom, corner
         dN(2, 1) = n1xi*dn1eta*n1zeta
         dN(2, 2) = n1xi*dn1eta*n3zeta
         dN(2, 3) = n3xi*dn1eta*n3zeta
         dN(2, 4) = n3xi*dn1eta*n1zeta

         ! Top, corner
         dN(2, 5) = n1xi*dn3eta*n1zeta
         dN(2, 6) = n1xi*dn3eta*n3zeta
         dN(2, 7) = n3xi*dn3eta*n3zeta
         dN(2, 8) = n3xi*dn3eta*n1zeta

         ! Bottom, edges
         dN(2, 9) = n1xi*dn1eta*n2zeta
         dN(2,10) = n2xi*dn1eta*n3zeta
         dN(2,11) = n3xi*dn1eta*n2zeta
         dN(2,12) = n2xi*dn1eta*n1zeta

         ! Top, edges
         dN(2,13) = n1xi*dn3eta*n2zeta
         dN(2,14) = n2xi*dn3eta*n3zeta
         dN(2,15) = n3xi*dn3eta*n2zeta
         dN(2,16) = n2xi*dn3eta*n1zeta

         ! Middle, edges
         dN(2,17) = n1xi*dn2eta*n1zeta
         dN(2,18) = n1xi*dn2eta*n3zeta
         dN(2,19) = n3xi*dn2eta*n3zeta
         dN(2,20) = n3xi*dn2eta*n1zeta

         ! Faces
         dN(2,21) = n2xi*dn2eta*n1zeta
         dN(2,22) = n2xi*dn2eta*n3zeta
         dN(2,23) = n1xi*dn2eta*n2zeta
         dN(2,24) = n3xi*dn2eta*n2zeta
         dN(2,25) = n2xi*dn1eta*n2zeta
         dN(2,26) = n2xi*dn3eta*n2zeta

         ! Volume
         dN(2,27) = n2xi*dn2eta*n2zeta

         ! Direction zeta

         ! Bottom, corner
         dN(3, 1) = n1xi*n1eta*dn1zeta
         dN(3, 2) = n1xi*n1eta*dn3zeta
         dN(3, 3) = n3xi*n1eta*dn3zeta
         dN(3, 4) = n3xi*n1eta*dn1zeta

         ! Top, corner
         dN(3, 5) = n1xi*n3eta*dn1zeta
         dN(3, 6) = n1xi*n3eta*dn3zeta
         dN(3, 7) = n3xi*n3eta*dn3zeta
         dN(3, 8) = n3xi*n3eta*dn1zeta

         ! Bottom, edges
         dN(3, 9) = n1xi*n1eta*dn2zeta
         dN(3,10) = n2xi*n1eta*dn3zeta
         dN(3,11) = n3xi*n1eta*dn2zeta
         dN(3,12) = n2xi*n1eta*dn1zeta

         ! Top, edges
         dN(3,13) = n1xi*n3eta*dn2zeta
         dN(3,14) = n2xi*n3eta*dn3zeta
         dN(3,15) = n3xi*n3eta*dn2zeta
         dN(3,16) = n2xi*n3eta*dn1zeta

         ! Middle, edges
         dN(3,17) = n1xi*n2eta*dn1zeta
         dN(3,18) = n1xi*n2eta*dn3zeta
         dN(3,19) = n3xi*n2eta*dn3zeta
         dN(3,20) = n3xi*n2eta*dn1zeta

         ! Faces
         dN(3,21) = n2xi*n2eta*dn1zeta
         dN(3,22) = n2xi*n2eta*dn3zeta
         dN(3,23) = n1xi*n2eta*dn2zeta
         dN(3,24) = n3xi*n2eta*dn2zeta
         dN(3,25) = n2xi*n1eta*dn2zeta
         dN(3,26) = n2xi*n3eta*dn2zeta

         ! Volume
         dN(3,27) = n2xi*n2eta*dn2zeta


      end subroutine hex27

      subroutine hex64(xi,eta,zeta,atoIJK,listHEX08,N,dN,N_lagrange,dN_lagrange)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Lagrangian HEX64 element model. Built using    !
         ! equispaced nodes between [-1,1] on             !
         ! (xi,eta,zeta). Ordering follows that of GMSH.  !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         integer(4), optional, intent(out) :: listHEX08(27,8)
         integer(4),           intent(out) :: atoIJK(64)
         real(8),              intent(in)  :: xi, eta, zeta
         real(8),    optional, intent(out) :: N(nnode), dN(ndime,nnode)
         real(8),    optional, intent(out) :: N_lagrange(nnode), dN_lagrange(ndime,nnode)
         real(8)                           :: xi_grid(porder+1)

         atoIJK = [1,4,11,12,2,3,15,16,9,20,33,34,10,19,36,35, &
                   5,8,27,28,6,7,29,30,25,32,53,56,26,31,54,55, &
                   13,23,41,44,17,21,45,46,37,50,57,60,38,49,58,59, &
                   14,24,42,43,18,22,48,47,40,51,61,64,39,52,62,63]

         if (present(listHEX08)) then
            listHEX08( 1,1:8) = [1,9,33,11,13,37,57,41]
            listHEX08( 2,1:8) = [9,10,36,33,37,38,58,57]
            listHEX08( 3,1:8) = [10,2,15,36,38,17,45,58]
            listHEX08( 4,1:8) = [11,33,34,12,41,57,60,44]
            listHEX08( 5,1:8) = [33,36,35,34,57,58,59,60]
            listHEX08( 6,1:8) = [36,15,16,35,58,45,46,59]
            listHEX08( 7,1:8) = [12,34,20,4,44,60,50,23]
            listHEX08( 8,1:8) = [34,35,19,20,60,59,49,50]
            listHEX08( 9,1:8) = [35,16,3,19,59,46,21,49]
            listHEX08(10,1:8) = [13,37,57,41,14,40,61,42]
            listHEX08(11,1:8) = [37,38,58,57,40,39,62,61]
            listHEX08(12,1:8) = [38,17,45,58,39,18,48,62]
            listHEX08(13,1:8) = [41,57,60,44,42,61,64,43]
            listHEX08(14,1:8) = [57,58,59,60,61,62,63,64]
            listHEX08(15,1:8) = [58,45,46,59,62,48,47,63]
            listHEX08(16,1:8) = [44,60,50,23,43,64,51,24]
            listHEX08(17,1:8) = [60,59,49,50,64,63,52,51]
            listHEX08(18,1:8) = [59,46,21,49,63,47,22,52]
            listHEX08(19,1:8) = [14,40,61,42,5,25,53,27]
            listHEX08(20,1:8) = [40,39,62,61,25,26,54,53]
            listHEX08(21,1:8) = [39,18,48,62,26,6,29,54]
            listHEX08(22,1:8) = [42,61,64,43,27,53,56,28]
            listHEX08(23,1:8) = [61,62,63,64,53,54,55,56]
            listHEX08(24,1:8) = [62,48,47,63,54,29,30,55]
            listHEX08(25,1:8) = [43,64,51,24,28,56,32,8]
            listHEX08(26,1:8) = [64,63,52,51,56,55,31,32]
            listHEX08(27,1:8) = [63,47,22,52,55,30,7,31]
         end if

         if (present(N) .and. present(dN)) then
            call lagrange_roots(xi_grid)
            call tripleTensorProduct(xi_grid,xi,eta,zeta,atoIJK,N,dN)
            if (flag_spectralElem == 1) then
               N_lagrange(:) = N(:)
               dN_lagrange(:,:) = dN(:,:)
               call chebyshev_roots(xi_grid)
               call tripleTensorProduct(xi_grid,xi,eta,zeta,atoIJK,N,dN)
            end if
         end if

      end subroutine hex64

      subroutine hexa_edges(ielem,nelem,npoin,connec,coord,ncorner,nedge,dist)

         implicit none
         
         integer(4), intent(in)            :: ielem, nelem, npoin
         integer(4), intent(in)            :: connec(nelem,nnode)
         real(8),    intent(in)            :: coord(npoin,ndime)
         integer(4), intent(out)           :: ncorner, nedge
         real(8),    intent(out)           :: dist(12,ndime)
         integer(4)                        :: ind(nnode)
         real(8)                           :: xp(12,ndime)
         
         ind = connec(ielem,:)
         ncorner = 8
         nedge = 12
         
         xp(1:8,1:ndime) = coord(ind(1:8),1:ndime) ! Corner coordinates
         dist(1,:) = xp(2,:)-xp(1,:)
         dist(2,:) = xp(3,:)-xp(2,:)
         dist(3,:) = xp(4,:)-xp(3,:)
         dist(4,:) = xp(1,:)-xp(4,:)
         
         dist(5,:) = xp(6,:)-xp(5,:)
         dist(6,:) = xp(7,:)-xp(6,:)
         dist(7,:) = xp(8,:)-xp(7,:)
         dist(8,:) = xp(5,:)-xp(8,:)
         
         dist(9,:) = xp(5,:)-xp(1,:)
         dist(10,:) = xp(6,:)-xp(2,:)
         dist(11,:) = xp(7,:)-xp(3,:)
         dist(12,:) = xp(8,:)-xp(4,:)
      
      end subroutine hexa_edges

end module
