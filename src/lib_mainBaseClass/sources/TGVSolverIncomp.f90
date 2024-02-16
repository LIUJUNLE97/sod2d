module TGVSolverIncomp_mod
   use mod_arrays
   use mod_nvtx
#ifndef NOACC
   use cudafor
#endif   
   use json_module
   use elem_qua
   use elem_hex
   use jacobian_oper
   use quadrature_rules
   use mod_inicond_reader
   use mass_matrix
   use mod_geom
   use time_integ
   use mod_analysis
   use mod_numerical_params
   use mod_time_ops
   use mod_fluid_viscosity
   use mod_postpro
   use mod_aver
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use CFDSolverPeriodicIncomp_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicIncomp) :: TGVSolverIncomp

      real(rp) , public  ::  Re

   contains
      procedure, public :: initializeParameters  => TGVSolverIncomp_initializeParameters
      procedure, public :: evalInitialConditions => TGVSolverIncomp_evalInitialConditions
   end type TGVSolverIncomp
contains

   subroutine TGVSolverIncomp_initializeParameters(this)
      class(TGVSolverIncomp), intent(inout) :: this
      logical :: found

      call this%loadJSONMeshFileData()
      call this%loadJSONResultsFileData()
      call this%loadJSONIOParams()

      call json%get("doGlobalAnalysis",this%doGlobalAnalysis, found)
      call json%get("doTimerAnalysis",this%doTimerAnalysis, found)

      call json%get("final_istep",this%final_istep, found)
      call json%get("maxPhysTime",this%maxPhysTime, found)

      call json%get("cfl_conv",this%cfl_conv, found)
      call json%get("cfl_diff",this%cfl_diff, found)

      call json%get("maxIter",maxIter, found)
      call json%get("tol",tol, found)

      call json%get("Re",this%Re, found)

      ! fixed by the type of base class parameters
      incomp_viscosity = 1.0_rp/this%Re
      flag_fs_fix_pressure = .false.
      flag_mu_factor = 1.0_rp
      nscbc_p_inf = 0.0_rp

   end subroutine TGVSolverIncomp_initializeParameters

   subroutine TGVSolverIncomp_evalInitialConditions(this)
      class(TGVSolverIncomp), intent(inout) :: this
      real(4) :: V0,L
      real(4) :: x,y,z
      integer(4) :: iNodeL,iNodeGSrl

      if(mpi_rank.eq.0) write(*,*) "--| TGV - Setting Initial Conditions..."

      V0 = 1.0_rp
      L  = 1.0_rp

      !$acc parallel loop
      do iNodeL=1,numNodesRankPar
         x = coordPar(iNodeL,1)
         y = coordPar(iNodeL,2)
         z = coordPar(iNodeL,3)

         u(iNodeL,1,2) =  V0*sin(x/(L))*cos(y/(L))*cos(z/(L))
         u(iNodeL,2,2) = -V0*cos(x/(L))*sin(y/(L))*cos(z/(L))
         u(iNodeL,3,2) = 0.0

         pr(iNodeL,2)  = nscbc_p_inf+((1.0_rp*V0*V0)/(16.0_rp))*(cos(2.0_rp*x/L)+cos(2.0_rp*y/L))*(cos(2.0_rp*z/L)+2.0_rp)
         rho(iNodeL,2) = 1.0_rp
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         eta(iNodeL,2) =0.5_rp*dot_product(u(iNodeL,1:ndime,2),u(iNodeL,1:ndime,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         u(iNodeL,1:ndime,3) = u(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         eta(iNodeL,3) =  eta(iNodeL,2)
         pr(iNodeL,3) =  pr(iNodeL,2)  

         q(iNodeL,1:ndime,4) = q(iNodeL,1:ndime,2)
         u(iNodeL,1:ndime,4) = u(iNodeL,1:ndime,2)
         rho(iNodeL,4) = rho(iNodeL,2)
         eta(iNodeL,4) =  eta(iNodeL,2)
         pr(iNodeL,4) =  pr(iNodeL,2)  

      end do
      !$acc end parallel loop

      !$acc kernels
      mu_e(:,:) = 0.0_rp ! Element syabilization viscosity
      mu_sgs(:,:) = 0.0_rp
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      !$acc end kernels
      call nvtxEndRange

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop
   end subroutine TGVSolverIncomp_evalInitialConditions

end module TGVSolverIncomp_mod
