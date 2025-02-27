module elem_stab

   use mod_nvtx
   use mod_numerical_params
   
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use mod_comms

      contains
        subroutine full_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Pr,rho_n,rho,u,Tem,Ml,ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ,tau_stab,Rmass,Rmom,Rener,initialze,fact)
             implicit none

             integer(4), intent(in)  :: nelem, npoin
             integer(4), intent(in)  :: connec(nelem,nnode)
             real(rp),   intent(in)  :: Ngp(ngaus,nnode)
             real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
             real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
             integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
             real(rp),   intent(in)  :: Cp,Pr,rho_n(npoin),rho(npoin),u(npoin,ndime), Tem(npoin),Ml(npoin)
             real(rp),   intent(in)  :: tau_stab(nelem),ProjMass(npoin,ndime),ProjEner(npoin,ndime),ProjMX(npoin,ndime),ProjMY(npoin,ndime),ProjMZ(npoin,ndime)
             real(rp),   intent(inout) :: Rmass(npoin)
             real(rp),   intent(inout) :: Rmom(npoin,ndime)
             real(rp),   intent(inout) :: Rener(npoin)
             logical, optional, intent(in)    :: initialze
             real(rp), optional, intent(in)  :: fact
             integer(4)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
             integer(4)              :: ipoin(nnode)
             real(rp)                :: gradU(ndime,ndime), gradT(ndime),tmp1,vol,arho
             real(rp)                :: gradIsoRho(ndime),gradIsoT(ndime),gradIsoU(ndime,ndime)
             real(rp)                :: gradRho(ndime),divDm(ndime),divDr,divDe,taustabl
             real(rp)                :: ul(nnode,ndime), rhol(nnode),Teml(nnode),mufluidl(nnode),rhonl(nnode)
             real(rp)                :: tauXl(nnode,ndime), tauYl(nnode,ndime), tauZl(nnode,ndime),divU
             real(rp)                :: gradTl(nnode,ndime),gradRhol(nnode,ndime),tauUl(nnode,ndime),tau(ndime,ndime)
             real(rp)                :: projMassl(nnode,ndime),ProjEnerl(nnode,ndime),ProjMXl(nnode,ndime),ProjMYl(nnode,ndime),ProjMZl(nnode,ndime)
             real(rp)  :: aux_fact = 1.0_rp, twoThirds


             twoThirds = 2.0_rp/3.0_rp

             call nvtxStartRange("Full stab")
             if(present(initialze)) then
               if (initialze .eqv. .true.) then
               !$acc kernels
                  Rmom(:,:) = 0.0_rp
                  Rmass(:) = 0.0_rp
                  Rener(:) = 0.0_rp
                  !$acc end kernels
               end if
            else
               !$acc kernels
               Rmom(:,:) = 0.0_rp
               Rmass(:) = 0.0_rp
               Rener(:) = 0.0_rp
               !$acc end kernels
            end if
            if(present(fact)) then
               aux_fact = fact
            end if

             !$acc parallel loop gang  private(ipoin,ul,Teml,rhol,rhonl,gradRhol,gradTl,tauXl,tauYl,tauZl,taustabl,projMassl,ProjEnerl,ProjMXl,ProjMYl,ProjMZl)
             do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
                   rhol(inode) = rho(ipoin(inode))
                   rhonl(inode) = rho_n(ipoin(inode))
                   Teml(inode) = Tem(ipoin(inode))      
                   !$acc loop seq          
                   do idime = 1,ndime
                      ul(inode,idime) = u(ipoin(inode),idime)
                      projMassl(inode,idime) = ProjMass(ipoin(inode),idime)
                      projEnerl(inode,idime) = ProjEner(ipoin(inode),idime)
                      projMXl(inode,idime) = ProjMX(ipoin(inode),idime)
                      projMYl(inode,idime) = ProjMY(ipoin(inode),idime)
                      projMZl(inode,idime) = ProjMZ(ipoin(inode),idime)
                   end do
                end do
                tauXl(:,:) = 0.0_rp
                tauYl(:,:) = 0.0_rp
                tauZl(:,:) = 0.0_rp
                gradTl(:,:) = 0.0_rp
                gradRhol(:,:) = 0.0_rp

                taustabl = tau_stab(ielem)

                !$acc loop vector private(tau,gradU,gradT,gradIsoRho,gradIsoT,gradIsoU,gradRho)
                do igaus = 1,ngaus

                   isoI = gmshAtoI(igaus) 
                   isoJ = gmshAtoJ(igaus) 
                   isoK = gmshAtoK(igaus) 

                   gradIsoRho(:) = 0.0_rp
                   gradIsoT(:) = 0.0_rp
                   gradIsoU(:,:) = 0.0_rp
                   !$acc loop seq
                   do ii=1,porder+1
                      gradIsoRho(1) = gradIsoRho(1) + dlxigp_ip(igaus,1,ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                      gradIsoRho(2) = gradIsoRho(2) + dlxigp_ip(igaus,2,ii)*rhol(invAtoIJK(isoI,ii,isoK))
                      gradIsoRho(3) = gradIsoRho(3) + dlxigp_ip(igaus,3,ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                      gradIsoT(1) = gradIsoT(1) + dlxigp_ip(igaus,1,ii)*Teml(invAtoIJK(ii,isoJ,isoK))
                      gradIsoT(2) = gradIsoT(2) + dlxigp_ip(igaus,2,ii)*Teml(invAtoIJK(isoI,ii,isoK))
                      gradIsoT(3) = gradIsoT(3) + dlxigp_ip(igaus,3,ii)*Teml(invAtoIJK(isoI,isoJ,ii))

                      !$acc loop seq
                      do idime=1,ndime
                         gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                         gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                         gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
                      end do
                   end do

                   gradRho(:) = 0.0_rp
                   gradT(:) = 0.0_rp
                   gradU(:,:) = 0.0_rp
                   !$acc loop seq
                   do idime=1, ndime
                      !$acc loop seq
                      do jdime=1, ndime
                         gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                         gradT(idime)   = gradT(idime)   + He(idime,jdime,igaus,ielem) * gradIsoT(jdime)
                         !$acc loop seq
                         do kdime=1,ndime
                            gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                         end do
                      end do
                   end do

                   divU = gradU(1,1)+gradU(2,2)+gradU(3,3)

                   !$acc loop seq
                   do idime = 1,ndime
                      !$acc loop seq
                      do jdime = 1,ndime
                         tau(idime,jdime) = (gradU(idime,jdime)+gradU(jdime,idime))
                      end do
                      tau(idime,idime) = tau(idime,idime)-twoThirds*divU
                   end do


                   !$acc loop seq
                   do idime = 1,ndime
                      tauXl(igaus,idime)    =  0.1_rp*taustabl*(projMXl(igaus,idime) - tau(1,idime))*rhonl(igaus)
                      tauYl(igaus,idime)    =  0.1_rp*taustabl*(projMYl(igaus,idime) - tau(2,idime))*rhonl(igaus)
                      tauZl(igaus,idime)    =  0.1_rp*taustabl*(projMZl(igaus,idime) - tau(3,idime))*rhonl(igaus)
                      gradTl(igaus,idime)   =  0.1_rp*taustabl*(projEnerl(igaus,idime) - gradT(idime))*rhonl(igaus)*Cp/Pr
                      gradRhol(igaus,idime) =  0.1_rp*taustabl*(projMassl(igaus,idime) - gradRho(idime))
                   end do
                end do

                !$acc loop vector private(divDm,divDr,divDe) 
                do igaus = 1,ngaus
                   isoI = gmshAtoI(igaus) 
                   isoJ = gmshAtoJ(igaus) 
                   isoK = gmshAtoK(igaus) 

                   divDe = 0.0_rp
                   divDr = 0.0_rp
                   divDm(:) = 0.0_rp
                   
                   !$acc loop seq
                   do ii=1,porder+1
                      !$acc loop seq
                      do idime=1,ndime
                         divDr = divDr + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*gradRhol(invAtoIJK(ii,isoJ,isoK),idime)
                         divDr = divDr + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*gradRhol(invAtoIJK(isoI,ii,isoK),idime)
                         divDr = divDr + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*gradRhol(invAtoIJK(isoI,isoJ,ii),idime)

                         divDe = divDe + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*(gradTl(invAtoIJK(ii,isoJ,isoK),idime))
                         divDe = divDe + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*(gradTl(invAtoIJK(isoI,ii,isoK),idime))
                         divDe = divDe + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*(gradTl(invAtoIJK(isoI,isoJ,ii),idime))
                        
                         divDm(1) = divDm(1) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauXl(invAtoIJK(ii,isoJ,isoK),idime)
                         divDm(1) = divDm(1) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauXl(invAtoIJK(isoI,ii,isoK),idime)
                         divDm(1) = divDm(1) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauXl(invAtoIJK(isoI,isoJ,ii),idime)

                         divDm(2) = divDm(2) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauYl(invAtoIJK(ii,isoJ,isoK),idime)
                         divDm(2) = divDm(2) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauYl(invAtoIJK(isoI,ii,isoK),idime)
                         divDm(2) = divDm(2) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauYl(invAtoIJK(isoI,isoJ,ii),idime)

                         divDm(3) = divDm(3) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauZl(invAtoIJK(ii,isoJ,isoK),idime)
                         divDm(3) = divDm(3) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauZl(invAtoIJK(isoI,ii,isoK),idime)
                         divDm(3) = divDm(3) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauZl(invAtoIJK(isoI,isoJ,ii),idime)
                      end do
                   end do

                   !$acc atomic update
                   Rmass(ipoin(igaus)) = Rmass(ipoin(igaus))+aux_fact*divDr
                   !$acc end atomic
                   !$acc atomic update
                   Rener(ipoin(igaus)) = Rener(ipoin(igaus))+aux_fact*divDe
                   !$acc end atomic
                   do idime = 1,ndime
                      !$acc atomic update
                      Rmom(ipoin(igaus),idime) = Rmom(ipoin(igaus),idime)+aux_fact*divDm(idime)
                      !$acc end atomic
                   end do
                end do
             end do
             !$acc end parallel loop
            call nvtxEndRange
        end subroutine full_stab_ijk

        subroutine full_diff_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Pr,rho_n,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Ml,ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ,tau_stab,Rmass,Rmom,Rener,initialze,fact)
         implicit none

         integer(4), intent(in)  :: nelem, npoin
         integer(4), intent(in)  :: connec(nelem,nnode)
         real(rp),   intent(in)  :: Ngp(ngaus,nnode)
         real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
         real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
         integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
         real(rp),   intent(in)  :: Cp,Pr,rho_n(npoin),rho(npoin),u(npoin,ndime), Tem(npoin),Ml(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
         real(rp),   intent(in)  :: tau_stab(nelem),ProjMass(npoin,ndime),ProjEner(npoin,ndime),ProjMX(npoin,ndime),ProjMY(npoin,ndime),ProjMZ(npoin,ndime)
         real(rp),   intent(in)  :: mu_fluid(npoin)
         real(rp),   intent(inout) :: Rmass(npoin)
         real(rp),   intent(inout) :: Rmom(npoin,ndime)
         real(rp),   intent(inout) :: Rener(npoin)
         logical, optional, intent(in)    :: initialze
         real(rp), optional, intent(in)  :: fact
         integer(4)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
         integer(4)              :: ipoin(nnode)
         real(rp)                :: gradU(ndime,ndime), gradT(ndime),tmp1,vol,arho,tauU(ndime)
         real(rp)                :: gradIsoRho(ndime),gradIsoT(ndime),gradIsoU(ndime,ndime),kappa_e, mu_fgp, mu_egp
         real(rp)                :: gradRho(ndime),divDm(ndime),divDr,divDe,taustabl
         real(rp)                :: ul(nnode,ndime), rhol(nnode),Teml(nnode),mufluidl(nnode),rhonl(nnode)
         real(rp)                :: tauXl(nnode,ndime), tauYl(nnode,ndime), tauZl(nnode,ndime),divU
         real(rp)                :: gradTl(nnode,ndime),gradRhol(nnode,ndime),tau(ndime,ndime),nu_e
         real(rp)                :: projMassl(nnode,ndime),ProjEnerl(nnode,ndime),ProjMXl(nnode,ndime),ProjMYl(nnode,ndime),ProjMZl(nnode,ndime)
         real(rp)  :: aux_fact = 1.0_rp, twoThirds


         twoThirds = 2.0_rp/3.0_rp

         call nvtxStartRange("Full stab")
         if(present(initialze)) then
           if (initialze .eqv. .true.) then
           !$acc kernels
              Rmom(:,:) = 0.0_rp
              Rmass(:) = 0.0_rp
              Rener(:) = 0.0_rp
              !$acc end kernels
           end if
        else
           !$acc kernels
           Rmom(:,:) = 0.0_rp
           Rmass(:) = 0.0_rp
           Rener(:) = 0.0_rp
           !$acc end kernels
        end if
        if(present(fact)) then
           aux_fact = fact
        end if

         !$acc parallel loop gang  private(ipoin,ul,Teml,rhol,rhonl,mufluidl,gradRhol,gradTl,tauXl,tauYl,tauZl,taustabl,projMassl,ProjEnerl,ProjMXl,ProjMYl,ProjMZl)
         do ielem = 1,nelem
            !$acc loop vector
            do inode = 1,nnode
              ipoin(inode) = connec(ielem,inode)
               rhol(inode) = rho(ipoin(inode))
               rhonl(inode) = rho_n(ipoin(inode))
               Teml(inode) = Tem(ipoin(inode))     
               mufluidl(inode) = mu_fluid(ipoin(inode)) 
               !$acc loop seq          
               do idime = 1,ndime
                  ul(inode,idime) = u(ipoin(inode),idime)
                  projMassl(inode,idime) = ProjMass(ipoin(inode),idime)
                  projEnerl(inode,idime) = ProjEner(ipoin(inode),idime)
                  projMXl(inode,idime) = ProjMX(ipoin(inode),idime)
                  projMYl(inode,idime) = ProjMY(ipoin(inode),idime)
                  projMZl(inode,idime) = ProjMZ(ipoin(inode),idime)
               end do
            end do
            tauXl(:,:) = 0.0_rp
            tauYl(:,:) = 0.0_rp
            tauZl(:,:) = 0.0_rp
            gradTl(:,:) = 0.0_rp
            gradRhol(:,:) = 0.0_rp

            taustabl = tau_stab(ielem)

            !$acc loop vector private(tau,gradU,gradT,gradIsoRho,gradIsoT,gradIsoU,gradRho)
            do igaus = 1,ngaus
               nu_e = c_rho*mu_e(ielem,igaus)/rhonl(igaus)
               mu_fgp = mufluidl(igaus)+rhonl(igaus)*mu_sgs(ielem,igaus)
               mu_egp = mu_e(ielem,igaus)
               kappa_e =mufluidl(igaus)*Cp/Pr+c_ener*mu_e(ielem,igaus)*Cp/Pr + Cp*rhonl(igaus)*mu_sgs(ielem,igaus)/0.9_rp

               isoI = gmshAtoI(igaus) 
               isoJ = gmshAtoJ(igaus) 
               isoK = gmshAtoK(igaus) 

               gradIsoRho(:) = 0.0_rp
               gradIsoT(:) = 0.0_rp
               gradIsoU(:,:) = 0.0_rp
               !$acc loop seq
               do ii=1,porder+1
                  gradIsoRho(1) = gradIsoRho(1) + dlxigp_ip(igaus,1,ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                  gradIsoRho(2) = gradIsoRho(2) + dlxigp_ip(igaus,2,ii)*rhol(invAtoIJK(isoI,ii,isoK))
                  gradIsoRho(3) = gradIsoRho(3) + dlxigp_ip(igaus,3,ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                  gradIsoT(1) = gradIsoT(1) + dlxigp_ip(igaus,1,ii)*Teml(invAtoIJK(ii,isoJ,isoK))
                  gradIsoT(2) = gradIsoT(2) + dlxigp_ip(igaus,2,ii)*Teml(invAtoIJK(isoI,ii,isoK))
                  gradIsoT(3) = gradIsoT(3) + dlxigp_ip(igaus,3,ii)*Teml(invAtoIJK(isoI,isoJ,ii))

                  !$acc loop seq
                  do idime=1,ndime
                     gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                     gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                     gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
                  end do
               end do

               gradRho(:) = 0.0_rp
               gradT(:) = 0.0_rp
               gradU(:,:) = 0.0_rp
               !$acc loop seq
               do idime=1, ndime
                  !$acc loop seq
                  do jdime=1, ndime
                     gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                     gradT(idime)   = gradT(idime)   + He(idime,jdime,igaus,ielem) * gradIsoT(jdime)
                     !$acc loop seq
                     do kdime=1,ndime
                        gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                     end do
                  end do
               end do

               divU = gradU(1,1)+gradU(2,2)+gradU(3,3)

               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     tau(idime,jdime) = (gradU(idime,jdime)+gradU(jdime,idime))
                  end do
                  tau(idime,idime) = tau(idime,idime)-twoThirds*divU
               end do

               !$acc loop seq
               do idime = 1,ndime
                  tauXl(igaus,idime)    =  -0.1_rp*taustabl*(projMXl(igaus,idime) - tau(1,idime))*rhonl(igaus)
                  tauYl(igaus,idime)    =  -0.1_rp*taustabl*(projMYl(igaus,idime) - tau(2,idime))*rhonl(igaus)
                  tauZl(igaus,idime)    =  -0.1_rp*taustabl*(projMZl(igaus,idime) - tau(3,idime))*rhonl(igaus)
                  gradTl(igaus,idime)   =  -0.1_rp*taustabl*(projEnerl(igaus,idime) - gradT(idime))*rhonl(igaus)*Cp/Pr
                  gradRhol(igaus,idime) =  -0.1_rp*taustabl*(projMassl(igaus,idime) - gradRho(idime))
               end do

               divU = gradU(1,1)+gradU(2,2)+gradU(3,3)

               tauU(:) = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     tauU(idime) = tauU(idime) + &
                        (mu_fgp+mu_egp)*(gradU(idime,jdime)+ gradU(jdime,idime))*ul(igaus,jdime)
                     tau(idime,jdime) = (mu_fgp+mu_egp)*(gradU(idime,jdime)+gradU(jdime,idime))
                  end do
                  tauU(idime) = tauU(idime)-(mu_fgp)*twoThirds*divU*ul(igaus,idime)
                  tau(idime,idime) = tau(idime,idime)-(mu_fgp)*twoThirds*divU
               end do

               !$acc loop seq
               do idime = 1,ndime
                  tauXl(igaus,idime) =  tauXl(igaus,idime) + tau(1,idime)
                  tauYl(igaus,idime) =  tauYl(igaus,idime) + tau(2,idime)
                  tauZl(igaus,idime) =  tauZl(igaus,idime) + tau(3,idime)
                  gradTl(igaus,idime) = gradTl(igaus,idime)+ kappa_e*gradT(idime) + tauU(idime)
                  gradRhol(igaus,idime) = gradRhol(igaus,idime) + nu_e*gradRho(idime)                  
               end do
            end do

            !$acc loop vector private(divDm,divDr,divDe) 
            do igaus = 1,ngaus
               isoI = gmshAtoI(igaus) 
               isoJ = gmshAtoJ(igaus) 
               isoK = gmshAtoK(igaus) 

               divDe = 0.0_rp
               divDr = 0.0_rp
               divDm(:) = 0.0_rp
               
               !$acc loop seq
               do ii=1,porder+1
                  !$acc loop seq
                  do idime=1,ndime
                     divDr = divDr + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*gradRhol(invAtoIJK(ii,isoJ,isoK),idime)
                     divDr = divDr + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*gradRhol(invAtoIJK(isoI,ii,isoK),idime)
                     divDr = divDr + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*gradRhol(invAtoIJK(isoI,isoJ,ii),idime)

                     divDe = divDe + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*(gradTl(invAtoIJK(ii,isoJ,isoK),idime))
                     divDe = divDe + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*(gradTl(invAtoIJK(isoI,ii,isoK),idime))
                     divDe = divDe + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*(gradTl(invAtoIJK(isoI,isoJ,ii),idime))
                    
                     divDm(1) = divDm(1) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauXl(invAtoIJK(ii,isoJ,isoK),idime)
                     divDm(1) = divDm(1) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauXl(invAtoIJK(isoI,ii,isoK),idime)
                     divDm(1) = divDm(1) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauXl(invAtoIJK(isoI,isoJ,ii),idime)

                     divDm(2) = divDm(2) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauYl(invAtoIJK(ii,isoJ,isoK),idime)
                     divDm(2) = divDm(2) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauYl(invAtoIJK(isoI,ii,isoK),idime)
                     divDm(2) = divDm(2) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauYl(invAtoIJK(isoI,isoJ,ii),idime)

                     divDm(3) = divDm(3) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauZl(invAtoIJK(ii,isoJ,isoK),idime)
                     divDm(3) = divDm(3) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauZl(invAtoIJK(isoI,ii,isoK),idime)
                     divDm(3) = divDm(3) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauZl(invAtoIJK(isoI,isoJ,ii),idime)
                  end do
               end do

               !$acc atomic update
               Rmass(ipoin(igaus)) = Rmass(ipoin(igaus))+aux_fact*divDr
               !$acc end atomic
               !$acc atomic update
               Rener(ipoin(igaus)) = Rener(ipoin(igaus))+aux_fact*divDe
               !$acc end atomic
               do idime = 1,ndime
                  !$acc atomic update
                  Rmom(ipoin(igaus),idime) = Rmom(ipoin(igaus),idime)+aux_fact*divDm(idime)
                  !$acc end atomic
               end do
            end do
         end do
         !$acc end parallel loop
        call nvtxEndRange
    end subroutine full_diff_stab_ijk        

        subroutine full_proj_ijk(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Pr,rho,u,Tem,Ml,ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ)
         use mod_solver, only : lumped_solver_vect
         implicit none

         integer(4), intent(in)  :: nelem, npoin,npoin_w
         integer(4), intent(in)  :: connec(nelem,nnode),lpoin_w(npoin_w)
         real(rp),   intent(in)  :: Ngp(ngaus,nnode)
         real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
         real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
         integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
         real(rp),   intent(in)  :: Cp,Pr,rho(npoin),u(npoin,ndime),Tem(npoin),Ml(npoin)
         real(rp),   intent(inout)  :: ProjMass(npoin,ndime),ProjEner(npoin,ndime),ProjMX(npoin,ndime),ProjMY(npoin,ndime),ProjMZ(npoin,ndime)
         integer(4)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
         integer(4)              :: ipoin(nnode)
         real(rp)                :: gradU(ndime,ndime), gradT(ndime),tmp1,vol,arho
         real(rp)                :: gradIsoRho(ndime),gradIsoT(ndime),gradIsoU(ndime,ndime)
         real(rp)                :: gradRho(ndime),divDm(ndime),divDr,divDe,taustabl,divU
         real(rp)                :: ul(nnode,ndime), rhol(nnode),Teml(nnode),mufluidl(nnode),tau(ndime,ndime)
         real(rp)                :: projMassl(nnode,ndime),ProjEnerl(nnode,ndime),ProjMXl(nnode,ndime),ProjMYl(nnode,ndime),ProjMZl(nnode,ndime)
         real(rp)  :: aux_fact = 1.0_rp, twoThirds

         call nvtxStartRange("Full stab")
         
         twoThirds = 2.0_rp/3.0_rp

         !$acc kernels
         ProjMass(:,:) = 0.0_rp
         ProjEner(:,:) = 0.0_rp
         ProjMX(:,:) = 0.0_rp
         ProjMY(:,:) = 0.0_rp
         ProjMZ(:,:) = 0.0_rp
         !$acc end kernels

         !$acc parallel loop gang  private(ipoin,ul,Teml,rhol)
         do ielem = 1,nelem
            !$acc loop vector
            do inode = 1,nnode
               ipoin(inode) = connec(ielem,inode)
               rhol(inode) = rho(ipoin(inode))
               Teml(inode) = Tem(ipoin(inode))
               !$acc loop seq
               do idime = 1,ndime
                  ul(inode,idime) = u(ipoin(inode),idime)
               end do
            end do

            !$acc loop vector private(tau,gradU,gradT,gradIsoRho,gradIsoT,gradIsoU,gradRho)
            do igaus = 1,ngaus
               isoI = gmshAtoI(igaus) 
               isoJ = gmshAtoJ(igaus) 
               isoK = gmshAtoK(igaus) 

               gradIsoRho(:) = 0.0_rp
               gradIsoT(:) = 0.0_rp
               gradIsoU(:,:) = 0.0_rp
               !$acc loop seq
               do ii=1,porder+1
                  gradIsoRho(1) = gradIsoRho(1) + dlxigp_ip(igaus,1,ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                  gradIsoRho(2) = gradIsoRho(2) + dlxigp_ip(igaus,2,ii)*rhol(invAtoIJK(isoI,ii,isoK))
                  gradIsoRho(3) = gradIsoRho(3) + dlxigp_ip(igaus,3,ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                  gradIsoT(1) = gradIsoT(1) + dlxigp_ip(igaus,1,ii)*Teml(invAtoIJK(ii,isoJ,isoK))
                  gradIsoT(2) = gradIsoT(2) + dlxigp_ip(igaus,2,ii)*Teml(invAtoIJK(isoI,ii,isoK))
                  gradIsoT(3) = gradIsoT(3) + dlxigp_ip(igaus,3,ii)*Teml(invAtoIJK(isoI,isoJ,ii))

                  !$acc loop seq
                  do idime=1,ndime
                     gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                     gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                     gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
                  end do
               end do

               gradRho(:) = 0.0_rp
               gradT(:) = 0.0_rp
               gradU(:,:) = 0.0_rp
               !$acc loop seq
               do idime=1, ndime
                  !$acc loop seq
                  do jdime=1, ndime
                     gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                     gradT(idime)   = gradT(idime)   + He(idime,jdime,igaus,ielem) * gradIsoT(jdime)
                     !$acc loop seq
                     do kdime=1,ndime
                        gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                     end do
                  end do
               end do

               divU = gradU(1,1)+gradU(2,2)+gradU(3,3)

               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     tau(idime,jdime) = (gradU(idime,jdime)+gradU(jdime,idime))
                  end do
                  tau(idime,idime) = tau(idime,idime)-twoThirds*divU
               end do


               !$acc loop seq
               do idime = 1,ndime
                     !$acc atomic update
                     ProjMX(ipoin(igaus),idime) = ProjMX(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*tau(1,idime)
                     !$acc end atomic
                     !$acc atomic update
                     ProjMY(ipoin(igaus),idime) = ProjMY(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*tau(2,idime)
                     !$acc end atomic
                     !$acc atomic update
                     ProjMZ(ipoin(igaus),idime) = ProjMZ(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*tau(3,idime)
                     !$acc end atomic
                     !$acc atomic update
                     ProjMass(ipoin(igaus),idime) = ProjMass(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*gradRho(idime)
                     !$acc end atomic
                     !$acc atomic update
                     ProjEner(ipoin(igaus),idime) = ProjEner(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*gradT(idime)
                     !$acc end atomic
               end do
            end do
         end do
         !$acc end parallel loop
         if(mpi_size.ge.2) then
            call nvtxStartRange("MPI_comms_tI")
            call mpi_halo_atomic_update_real_arrays(ndime,ProjMX(:,:))
            call mpi_halo_atomic_update_real_arrays(ndime,ProjMY(:,:))
            call mpi_halo_atomic_update_real_arrays(ndime,ProjMZ(:,:))
            call mpi_halo_atomic_update_real_arrays(ndime,ProjMass(:,:))
            call mpi_halo_atomic_update_real_arrays(ndime,ProjEner(:,:))
            call nvtxEndRange
         end if
         
         !
         ! Call lumped mass matrix solver
         !
      
         call nvtxStartRange("Call solver")
         call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,ProjMX(:,:))
         call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,ProjMY(:,:))
         call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,ProjMZ(:,:))
         call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,ProjMass(:,:))
         call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,ProjEner(:,:))
         
         call nvtxEndRange
    end subroutine full_proj_ijk


    subroutine comp_tau(nelem,npoin,connec,csound,u,helem_k,dt,tau)

      ! TODO: Compute element size h

      implicit none

      integer(4), intent(in)   :: nelem, npoin,connec(nelem,nnode)
      real(rp),    intent(in)  :: u(npoin,ndime),helem_k(nelem),dt,csound(npoin)
      real(rp),    intent(out) :: tau(nelem)
      integer(4)               :: ielem, inode
      real(rp)                 :: taul
      real(rp)                 :: aux1, fact_low_mach =1.0_rp

      if(flag_drop_c_in_envit) fact_low_mach = 0.0_rp

      !$acc parallel loop gang 
      do ielem = 1,nelem
          taul = 0.0_rp
          !$acc loop vector reduction(max:taul)
          do inode = 1,nnode
              aux1 = fact_low_mach*csound(connec(ielem,inode))+sqrt(dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))) ! Velocity mag. at element node
              taul = max(taul,(helem_k(ielem))*(c_species_stab/real(porder,rp))*aux1)
          end do
          tau(ielem) = taul
      end do
      !$acc end parallel loop

  end subroutine comp_tau
end module elem_stab
