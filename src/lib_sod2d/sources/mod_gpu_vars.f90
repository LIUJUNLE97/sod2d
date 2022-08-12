module mod_gpu_vars

   implicit none

#ifdef GPU

   ! Mesh information

   integer(4), device, allocatable :: connec_d(:,:)
   integer(4), device, allocatable :: lbnodes_d(:)

   ! Primary variables

   real(rp), device, allocatable    :: mu_e_d(:)
   real(rp), device, allocatable    :: rho_d(:,:)
   real(rp), device, allocatable    :: pr_d(:,:)
   real(rp), device, allocatable    :: E_d(:,:)
   real(rp), device, allocatable    :: Tem_d(:,:)
   real(rp), device, allocatable    :: e_int_d(:,:)
   real(rp), device, allocatable    :: u_d(:,:,:)
   real(rp), device, allocatable    :: q_d(:,:,:)

   ! Mass matrices

   real(rp), device, allocatable    :: Ml_d(:)
   real(rp), device, allocatable    :: Mc_d(:)

   ! Elemental info

   real(rp), device, allocatable    :: gpvol_d(:,:,:)
   real(rp), device, allocatable    :: gpcar_d(:,:,:,:)
   real(rp), device, allocatable    :: Ngp_d(:,:)

#endif

end module mod_gpu_vars
