subroutine gradKE(u, v, KE, KEx, KEy, k, G, GV)
  !$omp declare target
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  real, intent(in) :: u(:,:,:), v(:,:,:)
    !< Zonal and meridional velocity [L T-1 ~> m s-1]
  real, intent(out) :: KE(:,:)
    !< Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  real, intent(out) :: KEx(:,:), KEy(:,:)
    !< Acceleration due to kinetic energy gradient [L T-2 ~> m s-2]
  integer, intent(in) :: k
    !< Layer number to calculate for

  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  !! Can also be declared here

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  do concurrent (j=Jsq:Jeq+1, i=Isq:Ieq+1)
    KE(i,j) = ( ( (G%areaCu( I ,j)*(u( I ,j,k)*u( I ,j,k))) + &
                  (G%areaCu(I-1,j)*(u(I-1,j,k)*u(I-1,j,k))) ) + &
                ( (G%areaCv(i, J )*(v(i, J ,k)*v(i, J ,k))) + &
                  (G%areaCv(i,J-1)*(v(i,J-1,k)*v(i,J-1,k))) ) )*0.25*G%IareaT(i,j)
  enddo

  do concurrent (j=js:je, I=Isq:Ieq)
    KEx(I,j) = (KE(i+1,j) - KE(i,j)) * G%IdxCu(I,j)
  enddo

  do concurrent (J=Jsq:Jeq, i=is:ie)
    KEy(i,J) = (KE(i,j+1) - KE(i,j)) * G%IdyCv(i,J)
  enddo
end subroutine gradKE
