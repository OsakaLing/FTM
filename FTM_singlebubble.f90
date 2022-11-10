!-----------------------------------------------!
!   Program of cavity flow using SMAC and SOR   !
!-----------------------------------------------!

!計算に用いるパラメータの設定（グローバル変数として）
module parameter
    


    !-------------パラメータの設定-------------!
    integer, parameter ::       step_max = 10000        !最大ステップ
    integer, parameter ::       interval = 100          !データ出力の間隔
    double precision, parameter ::    dt = 0.005        !時間刻み
    integer, parameter ::             nx = 40           !x方向の格子数
    integer, parameter ::             ny = 80           !y方向の格子数
    double precision, parameter ::    Lx = 2.0          !x方向代表長さ
    double precision, parameter ::    Ly = 4.0          !y方向代表長さ
    double precision, parameter ::    Re = 400.0        !レイノルズ数
    double precision, parameter :: omega = 1.5          !加速係数(SOR法)
    double precision, parameter ::  rho1 = 2.0          !密度　(流体１)
    double precision, parameter ::  rho2 = 1.0          !密度　(流体２)
    double precision, parameter ::  myu1 = 0.2          !粘度　(流体１)
    double precision, parameter ::  myu2 = 0.1          !粘度　(流体２)
    double precision, parameter ::  sigma = 10.0         ! surface tension coefficient
    double precision, parameter ::  gx = 0.0            ! gravity coefficient of x direction
    double precision, parameter ::  gy = 100.0           ! gravity coefficient of y direction

    !計算に用いる微小量の設定
    double precision, parameter ::   dx = Lx / nx      !x刻み幅
    double precision, parameter ::   dy = Ly / ny      !y刻み幅
    double precision, parameter ::  eps = 1.0e-8        !収束判定用

    !計算に用いる逆数(inverse)を設定
    double precision ::   dxinv = 1.0 / dx         !dxの逆数(inverse)
    double precision ::   dyinv = 1.0 / dy         !dyの逆数
    double precision :: dxdxinv = 1.0 / dx / dx    !dx^2の逆数
    double precision :: dydyinv = 1.0 / dy / dy    !dy^2の逆数
    double precision ::   dtinv = 1.0 / dt         !dtの逆数
    double precision ::   Reinv = 1.0 / Re         !Reの逆数

    !気泡の初期設定
    integer         , parameter :: nf_first=100      ! Number of front element
    double precision, parameter :: pi=4.0*ATAN(1.0)
    double precision, parameter :: xcf=1.0,ycf=0.4  ! Axis of center of bubble
    double precision, parameter :: r=0.3           ! radius of bubble



    !---------境界条件の設定---------!
    !正方領域
    !上　==　一定速度で動く壁
    double precision, parameter ::    U_top = 0.0     !x方向に一定速度で動く
    double precision, parameter ::    V_top = 0.0     !壁面方向に速度を持たない
    !下，左，右　==　固定壁，滑りなし
    double precision, parameter :: U_bottom = 0.0     !滑りなし
    double precision, parameter :: V_bottom = 0.0     !壁面方向に速度を持たない．
    double precision, parameter ::   U_left = 0.0     !壁面方向
    double precision, parameter ::   V_left = 0.0     !滑りなし
    double precision, parameter ::  U_right = 0.0     !壁面方向
    double precision, parameter ::  V_right = 0.0     !滑りなし



    !---------the value need to be changed because of the reconstruction---------!
    double precision, allocatable :: xf(   :    )    !界面の x 座標
    double precision, allocatable :: yf(   :    )    !界面の y 座標
    double precision, allocatable :: tan_x(   :    ) ! x axis of tangent vectors of front elements
    double precision, allocatable :: tan_y(   :    ) ! y axis of tangent vectors of front elements

    integer                       :: nf              ! number of front elements


    
    

    


end module parameter




!-----------------------------------------------!
!                 Main Program                  !
!-----------------------------------------------!
program cavity
    use parameter
    implicit none

    !-------------変数の設定-------------!
    !スタッガード格子に基づく配列の設定
    double precision   u(0:nx   , 0:ny+1)   !x方向速度
    double precision   v(0:nx+1 , 0:ny  )   !y方向速度
    double precision   p(0:nx+1 , 0:ny+1)   !圧力場
    !SMAC法,SOR法で用いる変数
    double precision  up(0:nx   , 0:ny+1)   !x方向予測速度
    double precision  vp(0:nx+1 , 0:ny  )   !y方向予測速度
    double precision phi(0:nx+1 , 0:ny+1)   !補正圧力
    double precision rho(0:nx+1 , 0:ny+1)   !密度
    double precision myu(0:nx+1 , 0:ny+1)   !粘度
    ! parameters of Front tracking method 
    double precision  x (nx           )    !固定格子の x 座標
    double precision  y (ny           )    !固定格子の y 座標
    double precision surface_tension_x(0:nx+1 , 0:ny+1)   !surface tension of x direction
    double precision surface_tension_y(0:nx+1 , 0:ny+1)   !surface tension of y direction
    integer :: t = 0 , i , j
     
    nf = nf_first

    allocate(xf(0:nf*3))
    allocate(yf(0:nf*3))
    allocate(tan_x(0:nf+1))
    allocate(tan_y(0:nf+1))

    !open(1,file="front.txt",status="replace")
    !初期化
    u   = 0.0
    v   = 0.0
    p   = 0.0
    up  = 0.0
    vp  = 0.0
    phi = 0.0
    rho = rho1
    myu = myu1
    xf  = 0.0
    yf  = 0.0
    tan_x = 0.0
    tan_y = 0.0
    surface_tension_x = 0.0
    surface_tension_y = 0.0
    

    !-------------メインの処理-------------!
    call set_BC(u,v,p,rho,myu)
    call set_BC(up,vp,phi,rho,myu)
    ! call out_data(u,v,p,t)
    call initial_front()                                            !ok
    call initial_rho_density(rho,myu,x,y)                              !ok
    call find_tan()                                                          !ok
    call surface_tension (surface_tension_x,surface_tension_y)     !OK except specfic point of very small value

    !    do t = 1, 35
     do t = 1, 17
        call set_BC(u,v,p,rho,myu)
        call find_tan()
        call surface_tension (surface_tension_x,surface_tension_y)
     
        call calc_pre_velo(u,v,p,up,vp,rho,myu,surface_tension_x,surface_tension_y)     !予測速度の計算(SMAC)
        call poisson_SOR(up,vp,phi,rho,myu)         !補正圧力の計算(SOR)
        call upt_uvp(u,v,p,up,vp,phi,rho,myu)       !速度と圧力の更新
        
        call front_advect(u,v)
        call distribute_density_field(rho)
        call distribute_viscosity_field(rho,myu)
        call reconstruct()
        write(*,*)nf,t




    !     データ出力
    !     if (mod(t,interval)==0) then
            !  call out_data(u,v,p,t,rho,myu,phi)
            !  call check_div(u,v)
    !     end if

       end do

            !        do i = 1 , nf

            !     write(1,*)xf(i),yf(i)

            ! end do

            ! write(1,*)"                          "
            ! write(1,*)"                          "
            ! write(1,*)"                          "
            ! write(1,*)"                          "
    call out_data(u,v,p,t,rho,myu,phi)
    ! call out_data_vsGhia(u,v)

end program cavity


!-----------------------------------------------!
!                  Subroutines                  !
!-----------------------------------------------!


!-----------------------------------------------------!
!         Fluid field solution(SMAC and SOR)          !
!-----------------------------------------------------!

!スタッガード格子に基づき，境界条件および境界外側のパラメータを設定する
subroutine set_BC(u,v,p,rho,myu)
    use parameter
    implicit none

    double precision u(0:nx,0:ny+1), v(0:nx+1,0:ny), p(0:nx+1,0:ny+1), rho(0:nx+1,0:ny+1), myu(0:nx+1,0:ny+1)

    !----------上壁の境界条件の設定----------!
    u(0:nx,ny+1) = 2*U_top-u(0:nx,ny)   !境界条件と内部速度から外挿
    v(0:nx,ny  ) = V_top                !境界条件
    p(1:nx,ny+1) = p(1:nx,ny)           !壁近傍で圧力勾配なし
    rho(1:nx,ny+1) = rho(1:nx,ny)
    myu(1:nx,ny+1) = myu(1:nx,ny)

    !----------下壁の境界条件の設定----------!
    u(0:nx,0) = 2*U_bottom-u(0:nx,1)
    v(0:nx,0) = V_bottom
    p(1:nx,0) = p(1:nx,1)
    rho(1:nx,0) = rho(1:nx,1)
    myu(1:nx,0) = myu(1:nx,1)

    !----------右壁の境界条件の設定----------!
    u(nx  ,0:ny) = U_right
    v(nx+1,0:ny) = 2*V_right-v(nx,0:ny)
    p(nx+1,1:ny) = p(nx,1:ny)
    rho(nx+1,1:ny) = rho(nx,1:ny)
    myu(nx+1,1:ny) = myu(nx,1:ny)

    !----------左壁の境界条件の設定----------!
    u(0,0:ny) = U_left
    v(0,0:ny) = 2*V_left-v(1,0:ny)
    p(0,1:ny) = p(1,1:ny)
    rho(0,1:ny) = rho(1,1:ny)
    myu(0,1:ny) = myu(1,1:ny)

    !------------- 四隅の圧力 ---------------!
    p(0   ,0   ) = p(1 ,1 )
    p(nx+1,0   ) = p(nx ,1 )
    p(0   ,ny+1) = p(1 ,ny)
    p(nx+1,ny+1) = p(nx,ny)

    !------------- 四隅の密度 ---------------!
    rho(0   ,0   ) = rho(1 ,1 )
    rho(nx+1,0   ) = rho(nx ,1 )
    rho(0   ,ny+1) = rho(1 ,ny)
    rho(nx+1,ny+1) = rho(nx,ny)

    !------------- 四隅の粘度 ---------------!
    myu(0   ,0   ) = myu(1 ,1 )
    myu(nx+1,0   ) = myu(nx ,1 )
    myu(0   ,ny+1) = myu(1 ,ny)
    myu(nx+1,ny+1) = myu(nx,ny)

end subroutine set_BC


 !------------- adv和visc没有问题---------------!
 !------------- SMAC法における予測速度を求める ---------------!
!SMAC法における予測速度を求める．
subroutine calc_pre_velo(u,v,p,up,vp,rho,myu,surface_tension_x,surface_tension_y)
    use parameter
    implicit none
    
    !引数
    double precision u(0:nx,0:ny+1), v(0:nx+1,0:ny), p(0:nx+1,0:ny+1)
    double precision up(0:nx,0:ny+1), vp(0:nx+1,0:ny)
    double precision rho(0:nx+1 , 0:ny+1)   
    double precision myu(0:nx+1 , 0:ny+1)   
    double precision surface_tension_x(0:nx+1 , 0:ny+1)   !surface tension of x direction
    double precision surface_tension_y(0:nx+1 , 0:ny+1)   !surface tension of y direction
    double precision  a
    double precision uad(0:nx,0:ny+1) , vad(0:nx+1,0:ny)

    !計算に用いる変数
    double precision :: adv  = 0.0    !移流項
    double precision :: pres = 0.0    !圧力項
    double precision :: visc = 0.0    !粘性項

    double precision :: u_ave = 0.0   !平均流速u
    double precision :: v_ave = 0.0   !平均流速v

    integer i,j
    
    up = 0.0
    vp = 0.0
    a = 0.0
    uad = 0.0
    vad = 0.0





  
! open(100,file="compare.txt",status="replace")

    !予測速度upの計算
    do j = 1, ny
        do i = 1, nx-1
       
            !移流項の計算
            ! v_ave = 0.25*(v(i,j)+v(i+1,j)+v(i+1,j-1)+v(i,j-1))
            ! adv   = -0.5*(u(i,j)*(u(i+1,j)-u(i-1,j))*dxinv &
            !               +v_ave*(u(i,j+1)-u(i,j-1))*dyinv)
            adv   = 0.25*(-1)*(((u(i+1,j)+u(i,j))**2.0-(u(i,j)+u(i-1,j))**2.0)*dxinv &
                        +(((u(i,j+1)+u(i,j))*(v(i+1,j) + v(i,j))) &
                        -((u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1))))*dyinv) &
                        +surface_tension_x(i,j) / (0.5*(rho(i+1,j)+rho(i,j))) &
                        -gx &
                        +rho1 * gx / (0.5*(rho(i+1,j)+rho(i,j)))
            !圧力項の計算
            pres  = -(p(i+1,j)-p(i,j))*dxinv  / (0.5*(rho(i+1,j)+rho(i,j)))
            !粘性項の計算
            ! visc  = Reinv*((u(i+1,j)-2.0*u(i,j)+u(i-1,j))*dxdxinv &
            !               +(u(i,j+1)-2.0*u(i,j)+u(i,j-1))*dydyinv)
            visc =  (dxinv * ( 2.0 * myu(i+1,j) * (u(i+1,j) - u(i,j)) * dxinv &
                             - 2.0 * myu(i,j) * (u(i,j) - u(i-1,j)) * dxinv ) &
                   + dyinv * ( 0.25*( myu(i,j) + myu(i+1,j) + myu(i,j+1) + myu(i+1,j+1) ) &
                           * ( ( u(i,j+1) - u(i,j) )  * dyinv +  ( v(i+1,j) - v(i,j) ) * dxinv)  &
                           -   0.25*( myu(i,j) + myu(i+1,j) + myu(i,j-1) + myu(i+1,j-1) ) &
                           * ( ( u(i,j) - u(i,j-1) )  * dyinv +  ( v(i+1,j-1) - v(i,j-1) ) * dxinv)   )) &
                           / (0.5*(rho(i+1,j)+rho(i,j)))
            
            ! a =   dyinv * ( 0.25*( myu(i,j) + myu(i+1,j) + myu(i,j-1) + myu(i+1,j-1) ) &
            !                * ( ( u(i,j) - u(i,j-1) )  * dyinv +  ( v(i+1,j-1) - v(i,j-1) ) ) * dxinv )
             a =  dxinv * ( 2.0 * myu(i+1,j) * (u(i+1,j) - u(i,j)) * dxinv &
                             - 2.0 * myu(i,j) * (u(i,j) - u(i-1,j)) * dxinv )
                            
   
    ! write(100,*)i*dx,j*dy,visc
                           
          
            !予測速度upの計算
            ! up(i,j) = u(i,j) + dt*(adv + pres + visc)
            up(i,j) = u(i,j) + dt*(adv + pres + visc)
              
          
        end do
    end do

        !    close(100) 

     open(371,file="temporary.txt",status="replace")
    !予測速度vpの計算
    do i = 1, nx
        do j = 1, ny-1
        
            !移流項の計算
            ! u_ave = 0.25*(u(i,j)+u(i,j+1)+u(i-1,j+1)+u(i-1,j))
            ! adv   = -0.5*(u_ave*(v(i+1,j)-v(i-1,j))*dxinv &
            !            +v(i,j)*(v(i,j+1)-v(i,j-1))*dyinv)
            adv   = 0.25*(-1)*((((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j))) &
                                -((u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j))))*dxinv &
                                +((v(i,j+1)+v(i,j))**2.0-(v(i,j)+v(i,j-1))**2.0)*dyinv) &
                                +surface_tension_y(i,j) / (0.5*(rho(i,j)+rho(i,j+1))) &
                                -gy &
                                +rho1 * gy / (0.5*(rho(i,j+1)+rho(i,j)))
            !圧力項の計算
            pres  = -(p(i,j+1)-p(i,j))*dyinv / (0.5*(rho(i,j+1)+rho(i,j)))  
            !粘性項の計算
            ! visc  = Reinv*((v(i+1,j)-2*v(i,j)+v(i-1,j))*dxdxinv &
            !               +(v(i,j+1)-2*v(i,j)+v(i,j-1))*dydyinv)
            visc  =  (dyinv * ( 2.0 * myu(i,j+1) * (v(i,j+1) - v(i,j)) * dyinv &
                             - 2.0 * myu(i,j) * (v(i,j) - v(i,j-1)) * dyinv ) &
                   + dxinv * ( 0.25*( myu(i,j) + myu(i+1,j) + myu(i,j+1) + myu(i+1,j+1) ) &
                           * ( ( u(i,j+1) - u(i,j) )  * dyinv +  ( v(i+1,j) - v(i,j) )* dxinv )  &
                           -   0.25*( myu(i,j) + myu(i-1,j) + myu(i,j+1) + myu(i-1,j+1) ) &
                           * ( ( u(i-1,j+1) - u(i-1,j) )  * dyinv +  ( v(i,j) - v(i-1,j) ) * dxinv )  ) ) &
                           / (0.5*(rho(i,j+1)+rho(i,j)))   
            !予測速度vpの計算
             vp(i,j) = v(i,j) + dt*(adv + pres + visc)
        
                    write(371,*)i*dx,j*dy,surface_tension_x(i,j)

        end do
    end do

    !  open(371,file="temporary.txt",status="replace")
! ! !      open(101,file="python",status="old")
! ! !      do i = 1, nx
! ! !             do j = 1, ny
! ! ! read(101,*,end = 1)uad(i,j)
! ! !             end do
! ! !             end do
! ! !  1   close(101)

    !  do i = 1, nx
    !         do j = 1, ny
          
    !       write(371,*)i*dx,j*dy,up(i,j),vp(i,j)

    !     end do
    ! end do
close(371)

   
   write(*,*)"calc_pre_velo==================OK"
       
end subroutine calc_pre_velo


!SOR法を用いて補正圧力phiを計算する
subroutine poisson_SOR(up,vp,phi,rho,myu)
    use parameter
    implicit none

    double precision  up(0:nx   , 0:ny+1)   !x方向予測速度
    double precision  vp(0:nx+1 , 0:ny  )   !y方向予測速度
    double precision phi(0:nx+1 , 0:ny+1)   !補正圧力
    double precision phi_old(0:nx+1 , 0:ny+1)   !補正圧力
    double precision rho(0:nx+1 , 0:ny+1)   !密度
    double precision rhot(0:nx+1 , 0:ny+1)  !Density for calculation poisson
    double precision myu(0:nx+1 , 0:ny+1)   !粘度

    

    double precision ::   div1 = 0.0    ! coefficient of poisson 
    double precision ::   div2 = 0.0    ! coefficient of poisson 
    double precision ::  dphi = 0.0    !補正圧力の補正項 delta phi
    double precision :: resid = 0.0    !残差 residual error

    double precision :: error = 0.0    !誤差 -> 収束判定

    integer i,j,m,nmax
    

    phi = 0.0
    nmax = 2*nx*ny

    rhot = rho
    rhot(0     ,1:ny) = 1000.0
    rhot(nx+1  ,1:ny) = 1000.0
    rhot(1:nx,0)      = 1000.0
    rhot(1:nx,ny+1)   = 1000.0

   


open (442 , file = "phi.txt",status = "replace" )

    !SOR法による補正圧力 phiの求解
    do m = 1, nmax

        call set_BC(up,vp,phi,rho,myu)  !計算領域外の設定
        
        error = 0.0             !誤差（収束判定用）のリセット
        phi_old = phi

        
        do j = 1, ny
            do i = 1, nx
                !   予測速度を用いた連続の式 Div(up)の計算
                !-> 残差 resid の計算　->　補正項 dphi の計算
                !-> 補正圧力の更新
                    !    calculate the coefficient for poisson
                     div1 = (0.5*dtinv)*((up(i,j)-up(i-1,j))*dxinv+(vp(i,j)-vp(i,j-1))*dyinv)


                     div2=1.0/( dxinv*( dxinv/(rhot(i+1,j)+rhot(i,j))+&
                           dxinv/(rhot(i-1,j)+rhot(i,j)) )+           &
                           dyinv*(dyinv/(rhot(i,j+1)+rhot(i,j))+      &
                           dyinv/(rhot(i,j-1)+rhot(i,j)) ) )
                    !    calculate the next step pressure  
                             
                      if(m == 1)then 
                     write(442,*)i*dx,j*dy,up(i,j)

              
                        end if
                  
                     phi(i,j) = ( 1.0 - omega )*phi(i,j) + omega*div2*( dxinv*( phi(i+1,j)*dxinv/(rhot(i+1,j)+rhot(i,j))  &
                                                                          +     phi(i-1,j)*dxinv/(rhot(i-1,j)+rhot(i,j)) ) & 
                                                                       + dyinv*(phi(i,j+1)*dyinv/(rhot(i,j+1)+rhot(i,j))  &
                                                                          +     phi(i,j-1)*dyinv/(rhot(i,j-1)+rhot(i,j)))  -div1)                              
                   
                !     resid =  ((phi(i+1,j)-phi(i,j))/(rho(i+1,j)+rho(i,j))+(phi(i-1,j)-phi(i,j))/(rho(i-1,j)+rho(i,j)))*dxdxinv*2.0 &
                !             +((phi(i,j+1)-phi(i,j))/(rho(i+1,j)+rho(i,j))+(phi(i,j-1)-phi(i,j))/(rho(i+1,j)+rho(i,j)))*dydyinv*2.0 &
                !             - div1
                !     ! dphi = (0.5*resid)/(dxdxinv+dydyinv)
                !     dphi = resid &
                !     /(dxdxinv*(1.0/(rho(i+1,j)+rho(i,j))+1.0/(rho(i-1,j)+rho(i,j))) &
                !      +dydyinv*(1.0/(rho(i,j+1)+rho(i,j))+1.0/(rho(i,j-1)+rho(i,j))))
                ! phi(i,j) = phi(i,j) + omega * dphi
   
                  
                !誤差：各補正項の相対誤差の最大値とする
                error = max(error,abs(phi_old(i,j)-phi(i,j))/abs(phi_old(i,j)))
            end do
        end do
       
        phi(nx+1,1:ny) = phi(nx,1:ny)
        phi(0,1:ny)    = phi(1,1:ny)
        phi(1:nx,ny+1) = phi(1:nx,ny)
        phi(1:nx,0)    = phi(1:nx,1)
         
       
            
           

        


        if(error<eps) exit  !収束判定
        
    end do

    write(*,*)'SOR:m=',m,', Error=',error

    close(442)

   write(*,*)"poisson_SOR==================OK"
end subroutine poisson_SOR

!予測速度と補正圧力を用いて，速度u,v圧力pの時間ステップを進める．
subroutine upt_uvp(u,v,p,up,vp,phi,rho,myu)
    use parameter
    implicit none
    
    double precision   u(0:nx   , 0:ny+1)   !x方向速度
    double precision   v(0:nx+1 , 0:ny  )   !y方向速度
    double precision   p(0:nx+1 , 0:ny+1)   !圧力場
    double precision  up(0:nx   , 0:ny+1)   !x方向予測速度
    double precision  vp(0:nx+1 , 0:ny  )   !y方向予測速度
    double precision phi(0:nx+1 , 0:ny+1)   !補正圧力
    double precision rho(0:nx+1 , 0:ny+1)   !密度
    double precision myu(0:nx+1 , 0:ny+1)   !粘度

    integer i,j

    call set_BC(up,vp,phi,rho,myu)

    !x方向速度 uの更新：計算領域に注意
    do j = 1, ny
        do i = 1, nx-1
            u(i,j) = up(i,j)-(phi(i+1,j)-phi(i,j))*dxinv*dt/(0.5*(rho(i+1,j)+rho(i,j)))
            !u(i,j) = up(i,j)-phi(i+1,j)*dxinv*dtinv+phi(i,j)*dxinv*dt
        end do
    end do

    !y方向速度 vの更新
    do j = 1, ny-1
        do i = 1, nx
            v(i,j) = vp(i,j)-(phi(i,j+1)-phi(i,j))*dyinv*dt/(0.5*(rho(i,j+1)+rho(i,j)))
            !v(i,j) = vp(i,j)-phi(i,j+1)*dyinv*dtinv+phi(i,j)*dyinv*dt
        end do
    end do

    !圧力の更新
    do j = 1, ny
        do i = 1, nx
            p(i,j)=p(i,j)+phi(i,j)
        end do
    end do





   write(*,*)"upt_uvp==================OK"

end subroutine upt_uvp

!計算したデータのアウトプット
!速度はスタッガード格子からレギュラー格子に直す
subroutine out_data(u,v,p,t,rho,myu,phi)
    use parameter
    implicit none

    !スタッガード格子
    double precision   u(0:nx   , 0:ny+1)   !x方向速度
    double precision   v(0:nx+1 , 0:ny  )   !y方向速度
    double precision   p(0:nx+1 , 0:ny+1)   !圧力場
    double precision rho(0:nx+1 , 0:ny+1)   !密度
    double precision myu(0:nx+1 , 0:ny+1)   !粘度
    double precision phi(0:nx+1 , 0:ny+1)   !補正圧力

    !出力(レギュラー)
    double precision   u_out(1:nx , 1:ny)   !x方向速度
    double precision   v_out(1:nx , 1:ny)   !y方向速度

    double precision x,y    !座標

    integer i,j,t
    character fname*30

    u_out(1:nx,1:ny) = (u(0:nx-1,1:ny)+u(1:nx,1:ny))*0.5
    v_out(1:nx,1:ny) = (v(1:nx,0:ny-1)+v(1:nx,1:ny))*0.5

    !ファイルを上書き(replace)モードで開く．
    ! write(fname,'("data/uvp"i7.7".csv")')t     ! ex) fname=uvp0000100
    ! open(10, file = fname, status = 'replace')
    open(10, file = "front.txt", status = 'replace')
    open(100,file="compare.txt",status="replace")
    write(100,*) '# x , y , u , v , rho , myu'

    ! do i = 1, nx
    !     x = (i-0.5)*dx
    !     do j = 1, ny
    !         y = (j-0.5)*dy
    !         write(10,*) x,',',y,',',u_out(i,j),',',v_out(i,j),',',p(i,j)
    !     end do
    !     write(10,*)     !空白行を入れる（プロット用）
    ! end do

            Do i = 1 , nx
              do j = 1 , ny             
                      write(100,*)(real(i)-0.5)*dx,(real(j)-0.5)*dy,u_out(i,j),v_out(i,j),rho(i,j),myu(i,j),phi(i,j)
                end do
            end do
    close(100)

    Do i = 1 , nf+1
                          
        write(10,*)xf(i),yf(i)
               
    end do

    close(10)

 

end subroutine out_data

!Ghiaとの比較のためのデータ抽出
subroutine out_data_vsGhia(u,v)
    use parameter
    implicit none

    double precision   u(0:nx   , 0:ny+1)   !x方向速度
    double precision   v(0:nx+1 , 0:ny  )   !y方向速度

    double precision x,y    !座標

    integer :: i = 0
    character fname*30

    !ファイルを上書き(replace)モードで開く．
    write(fname,'("Ghia_yu.csv")')
    open(20, file = fname, status = 'replace')
    write(20,*) '# x=0.5'
    write(20,*) '# y , u'

    write(20,*)'0, 0'
    do i = 1, ny
        write(20,*) (i-0.5)*dy,',',u(nx/2,i)
    end do
    write(20,*)'1.0, ',U_top

    close(20)

    write(fname,'("Ghia_xv.csv")')
    open(21, file = fname, status = 'replace')
    write(21,*) '# y=0.5'
    write(21,*) '# x , v'

    write(21,*)'0, 0'
    do i = 1, nx
        write(21,*) (i-0.5)*dx,',',v(i,ny/2)
    end do
    write(21,*)'1.0, 0'

    close(21)

end subroutine out_data_vsGhia


subroutine check_div(u,v)
    use parameter
    implicit none

    double precision   u(0:nx   , 0:ny+1)   !x方向速度
    double precision   v(0:nx+1 , 0:ny  )   !y方向速度

    double precision   divergence(1:nx,1:ny)

    integer i,j

    do j = 1, ny
        do i = 1, nx
            divergence(i,j) &
            = (u(i,j)-u(i-1,j))*dxinv &
             +(v(i,j)-v(i,j-1))*dyinv
        end do
    end do

     write(*,*)MAXVAL(divergence)

end subroutine check_div



!-----------------------------------------------!
!                  Front tracking               !
!-----------------------------------------------!

!---------OK---------!
!---------initialize the front---------!
subroutine initial_front()
    use parameter
    implicit none

    Integer::i,j

    !Set initial front axis 
    Do i=1,nf

        xf(i)=xcf-r*sin(2*pi/real(nf)*real(i-1))

    end do 
    Do j=1,nf

        yf(j)=ycf+r*cos(2*pi/real(nf)*real(j-1))

    end do
    
    xf(nf+1)=xf(1)
    yf(nf+1)=yf(1)
    xf(0)   =xf(nf)
    yf(0)   =yf(nf)
    
    ! open(100,file="compare.txt",status="new")
    ! Do i = 0, nf+1
    ! write(100,*)xf(i),yf(i)
    ! end do
    ! close(100)

 write(*,*)"initial_front==================OK"
end subroutine initial_front


!---------OK---------!
!---------find tangent vector---------!
subroutine find_tan()

    use parameter
    implicit none



    double precision dist
    integer i
    
    open(736,file="tangent.txt",status="replace")
    ! initilize
    tan_x = 0.0
    tan_y = 0.0
    dist  = 0.0

    Do i = 1 , nf

        dist     = sqrt( (xf(i+1)-xf(i))**2 + (yf(i+1)-yf(i))**2 )
        tan_x(i) = ( xf(i+1) - xf(i) ) / dist
        tan_y(i) = ( yf(i+1) - yf(i) ) / dist
        write(736,*)dist
    end do
    ! connect the end and beginning
    tan_x(0) = tan_x(nf)
    tan_y(0) = tan_y(nf)
    tan_x(nf+1) = tan_x(1)
    tan_y(nf+1) = tan_y(1)

    close(736)
  
 write(*,*)"find_tan==================OK"
end subroutine find_tan


!---------OK---------!
!---------set initial density and viscosity---------!
subroutine initial_rho_density (rho,myu,x,y)
    use parameter 
    implicit none

    double precision rho(0:nx+1 , 0:ny+1)   !密度
    double precision myu(0:nx+1 , 0:ny+1)   !粘度
    double precision  x(nx           )    !固定格子の x 座標
    double precision  y(ny           )    !固定格子の y 座標

    integer i ,j
    

    ! find the location of the center of the grid
    Do i = 1 , nx
    
      x(i) = (real(i)-0.5)*dx

    End do

    Do j = 1 , ny
    
      y(j) = (real(j)-0.5)*dy

    End do



    Do i = 1 , nx
      Do j = 1 , ny
      
        if ((x(i)-xcf)**2 + (y(j)-ycf)**2 < r**2 ) then
        
          rho(i,j) = rho2
          myu(i,j) = myu2



        end if


      End do
    End do


  !    write (*,*)myu


 write(*,*)"initial_rho_density==================OK"
end subroutine initial_rho_density


!---------OK except specfic point of very small value---------!
!---------calculate the surface tension---------!

subroutine surface_tension (surface_tension_x,surface_tension_y)

    use parameter 
    implicit none


    double precision surface_tension_x(0:nx+1 , 0:ny+1)   !surface tension of x direction
    double precision surface_tension_y(0:nx+1 , 0:ny+1)   !surface tension of y direction
    double precision force_x , force_y
    double precision coeff_x , coeff_y
    double precision xi(0:nx+1) , yi(0:nx+1)

    integer :: i ,j , fi , fj

    surface_tension_x = 0.0
    surface_tension_y = 0.0
    force_x = 0.0
    force_y = 0.0
    xi = 0.0
    yi = 0.0

    Do i = 0 , nx
      xi(i) = i * dx
    End do

    Do i = 0 , nx
      yi(i) = i * dy - 0.5 * dy 
    End do

    DO i = 1 , nf
      ! calculate force of x direction
      force_x = sigma * ( tan_x(i) - tan_x(i-1) )
      ! distribute to Euler
      fi = floor( xf(i)/dx )
      fj = floor (( yf(i) + 0.5*dy ) /dy)

      coeff_x = xf(i)/dx - fi
      coeff_y = ( yf(i)+0.5*dy )/dy -fj


      surface_tension_x( fi , fj ) = surface_tension_x( fi , fj ) + &
                                     (1.0 - coeff_x)*(1.0 - coeff_y) * &
                                      force_x /dx/dy
      surface_tension_x( fi+1 , fj ) = surface_tension_x( fi+1 , fj ) + &
                                     coeff_x*(1.0 - coeff_y) * &
                                      force_x /dx/dy
      surface_tension_x( fi , fj+1 ) = surface_tension_x( fi , fj+1 ) + &
                                     (1.0 - coeff_x)*coeff_y * &
                                      force_x /dx/dy
      surface_tension_x( fi+1 , fj+1 ) = surface_tension_x( fi+1 , fj+1 ) + &
                                      coeff_x*coeff_y * &
                                      force_x /dx/dy

      ! write (*,*) surface_tension_x( fi , fj ) ,xi(fi) , yi (fj)
      ! write (*,*) surface_tension_x( fi+1 , fj ) ,xi(fi+1) , yi (fj)
      ! write (*,*) surface_tension_x( fi , fj+1 ) ,xi(fi) , yi (fj+1)
      ! write (*,*) surface_tension_x( fi+1 , fj+1 ) ,xi(fi+1) , yi (fj+1)

      ! calculate force of y direction
       force_y = sigma * ( tan_y(i) - tan_y(i-1) )
      !  distribute to lagrange
      fi = floor ( (xf(i) + 0.5*dx)/dx )
      fj = floor (( yf(i)  ) /dy)


      coeff_x = (xf(i) + 0.5*dx)/dx  - fi
      coeff_y = ( yf(i)  ) /dy -fj


      surface_tension_y( fi , fj ) = surface_tension_y( fi , fj ) + &
                                     (1.0 - coeff_x)*(1.0 - coeff_y) * &
                                      force_y /dx/dy
      surface_tension_y( fi+1 , fj ) = surface_tension_y( fi+1 , fj ) + &
                                     coeff_x*(1.0 - coeff_y) * &
                                      force_y /dx/dy
      surface_tension_y( fi , fj+1 ) = surface_tension_y( fi , fj+1 ) + &
                                     (1.0 - coeff_x)*coeff_y * &
                                      force_y /dx/dy
      surface_tension_y( fi+1 , fj+1 ) = surface_tension_y( fi+1 , fj+1 ) + &
                                      coeff_x*coeff_y * &
                                      force_y /dx/dy
    !   write (*,*) surface_tension_y( fi , fj ) ,xi(fi) , yi (fj)
    !   write (*,*) surface_tension_y( fi+1 , fj ) ,xi(fi+1) , yi (fj)
    !   write (*,*) surface_tension_y( fi , fj+1 ) ,xi(fi) , yi (fj+1)
    !   write (*,*) surface_tension_y( fi+1 , fj+1 ) ,xi(fi+1) , yi (fj+1)


    End do




        open(902,file="surfacetension.txt",status="replace")
    Do i = 1, nx
      do j = 1 , ny
    write(902,*)i*dx,j*dy,surface_tension_x(i,j),surface_tension_y(i,j)
    end do
    end do
    close(902)



   write(*,*)"surface_tension==================OK"
end subroutine surface_tension


!---------advect the front---------!

subroutine front_advect ( u , v)

    use parameter
    implicit none


    double precision u(0:nx   , 0:ny+1)   !velocity of grid of x direction
    double precision v(0:nx+1 , 0:ny  )   !velocity of grid of y direction
    double precision uf(nf+1           )   !velocity of front of x direction
    double precision vf(nf+1           )   !velocity of front of y direction
    double precision coeff_x , coeff_y
  
    integer :: i , fi , fj 

    uf = 0.0
    vf = 0.0


    DO i = 1 , nf
       ! distribute to lagrange of u
       ! determine the nearest grid
        fi = floor( xf(i)/dx )
        fj = floor (( yf(i) + 0.5*dy ) /dy)

        !calculate the weight
        coeff_x = xf(i)/dx - fi
        coeff_y = ( yf(i)+0.5*dy )/dy -fj

       ! interpolate
        uf(i) = uf(i) + (1.0 - coeff_x)*(1.0 - coeff_y) * u( fi , fj ) 
        uf(i) = uf(i) + coeff_x*(1.0 - coeff_y) * u( fi+1 , fj ) 
        uf(i) = uf(i) + (1.0 - coeff_x)* coeff_y * u( fi , fj+1 ) 
        uf(i) = uf(i) + coeff_x* coeff_y * u( fi+1 , fj+1 ) 

        ! distribute to lagrange of v
        ! determine the nearest grid
        fi = floor ( (xf(i) + 0.5*dx)/dx )
        fj = floor (( yf(i)  ) /dy)

        !calculate the weight
        coeff_x = (xf(i) + 0.5*dx)/dx - fi
        coeff_y = ( yf(i)  ) /dy -fj

        ! interpolate
        vf(i) = vf(i) + (1.0 - coeff_x)*(1.0 - coeff_y) * v( fi , fj ) 
        vf(i) = vf(i) + coeff_x*(1.0 - coeff_y) * v( fi+1 , fj ) 
        vf(i) = vf(i) + (1.0 - coeff_x)* coeff_y * v( fi , fj+1 ) 
        vf(i) = vf(i) + coeff_x* coeff_y * v( fi+1 , fj+1 ) 
        
        !advect the front
        xf(i) = xf(i) + uf(i) * dt
        yf(i) = yf(i) + vf(i) * dt

        !connect the end and beginning
        xf(0)    = xf(nf)
        yf(0)    = yf(nf)
        xf(nf+1) = xf(1)
        yf(nf+1) = yf(1)
   

    End do

    open(987,file="frontvelocity.txt",status="replace")


     do i = 1 , nf

    write(987,*)xf(i),yf(i),uf(i),vf(i),sqrt(uf(i)**2+vf(i)**2)
  
    end do
    close(987)
                 

                
   write(*,*)"front_advect==================OK"
end subroutine front_advect

!---------reconstruct front---------!
subroutine reconstruct()

    use parameter 
    implicit none




    double precision , allocatable :: xf_old(:)
    double precision , allocatable :: yf_old(:)

    integer :: i , nfback , j ,nfr
    double precision :: dis
    
    ! the array to store old front information
    allocate(xf_old(0:(nf)*3))
    allocate(yf_old(0:(nf)*3))

    !     open(1095,file="front.txt",status="replace")
   
    !   do i = 0 , nf+1
    !  write(1095,*) xf(i),yf(i)

    !   end do

    ! close(1095)
    
    xf_old = xf
    yf_old = yf
    j = 0
    nfr = nf
  
   
    !   open(1095,file="front.txt",status="replace")
   
    !   do i = 0 , nf+1
    !  write(1095,*)  xf_old(i), yf_old(i)

    !   end do

    ! close(1095)
    

    


    do i = 2 , nfr+1
        
        dis = sqrt ( ((xf_old(i)-xf(j))/dx)**2 + ((yf_old(i)-yf(j))/dy)**2 )

        
        if (dis > 0.5)then
            
             
            j = j+1
            xf(j) = 0.5 * (xf_old(i)+xf(j-1))
            yf(j) = 0.5 * (yf_old(i)+yf(j-1))
            j = j+1
            xf(j) = xf_old(i)
            yf(j) = yf_old(i) 
            nf = nf +1
           
        

            else if (dis < 0.25) then 
                
                nf = nf-1
                continue

                else
                    j = j+1
                    xf(j) = xf_old(i)
                    yf(j) = yf_old(i)
                   
         
               
        end if
   
        !   write(*,*)j,i      
    end do 
    ! connect end and beginning of front
    write(*,*)nf
    xf(0)=xf(nf)
    yf(0)=yf(nf)
    xf(nf+1)= xf(1)
    yf(nf+1)= yf(1)
   
      

     
    ! reshape the array
    deallocate(xf_old)
    deallocate(yf_old)
    allocate(xf_old(0:nf+1))
    allocate(yf_old(0:nf+1))

xf_old=0.0
yf_old=0.0
    do i = 0 , nf+1

        xf_old(i) = xf(i) 
        yf_old(i) = yf(i) 

    end do
    
    !   open(1095,file="front.txt",status="replace")
   
    !   do i = 0 , nf+1
    !  write(1095,*)  xf_old(i), yf_old(i)

    !   end do

    ! close(1095)

    deallocate(xf)
    deallocate(yf)
    allocate(xf(0:nf*3))
    allocate(yf(0:nf*3))
    
    xf = 0.0
    yf = 0.0
    do i = 0 , nf+1

        xf(i) = xf_old(i) 
        yf(i) = yf_old(i) 

    end do
    
    xf(0)=xf(nf)
    yf(0)=yf(nf)
    xf(nf+1)= xf(1)
    yf(nf+1)= yf(1)

    open(1095,file="front.txt",status="replace")
   
      do i = 0 , nf+1
     write(1095,*) xf(i),yf(i)

      end do

    close(1095)
  
   write(*,*)"reconstruct==================OK"
end subroutine reconstruct



!---------distribute new density---------!
subroutine distribute_density_field(rho)

use parameter 
implicit none

double precision rho(0:nx+1 , 0:ny+1)   
double precision drhox(0:nx+1 , 0:ny+1)  ! drivation of rho of x direction
double precision drhoy(0:nx+1 , 0:ny+1)  ! drivation of rho of y direction
double precision force_x , force_y
double precision coeff_x , coeff_y
double precision xi(0:nx+1) , yi(0:nx+1)
double precision dphi  , div  ,  error  ,  resid
Integer          nmax  , m

    integer :: i ,j , fi , fj
    
    drhox   = 0.0
    drhoy   = 0.0
    force_x = 0.0
    force_y = 0.0
  
           

  

    DO i = 1 , nf
      ! calculate jump of x direction
      force_x = -0.5*(yf(i+1)-yf(i-1))* ( rho2 - rho1 )

      ! distribute to Euler
      fi = floor( xf(i)/dx )
      fj = floor (( yf(i) + 0.5*dy ) /dy)
      coeff_x = xf(i)/dx - fi
      coeff_y = ( yf(i)+0.5*dy )/dy -fj

      drhox( fi , fj )     = drhox( fi , fj ) + &
                             (1.0 - coeff_x)*(1.0 - coeff_y) * &
                             force_x /dx/dy
      drhox( fi+1 , fj )   = drhox( fi+1 , fj ) + &
                             coeff_x*(1.0 - coeff_y) * &
                             force_x /dx/dy
      drhox( fi , fj+1 )   = drhox( fi , fj+1 ) + &
                             (1.0 - coeff_x)*coeff_y * &
                             force_x /dx/dy
      drhox( fi+1 , fj+1 ) = drhox( fi+1 , fj+1 ) + &
                             coeff_x*coeff_y * &
                             force_x /dx/dy

    

      ! calculate jump of y direction
       force_y = 0.5*(xf(i+1)-xf(i-1))* ( rho2 - rho1 )

      !  distribute to lagrange
      fi = floor ( (xf(i) + 0.5*dx)/dx )
      fj = floor (( yf(i)  ) /dy)
      coeff_x = (xf(i) + 0.5*dx)/dx - fi
      coeff_y = ( yf(i)  ) /dy -fj


      drhoy( fi , fj )     = drhoy( fi , fj ) + &
                             (1.0 - coeff_x)*(1.0 - coeff_y) * &
                             force_y /dx/dy
      drhoy( fi+1 , fj )   = drhoy( fi+1 , fj ) + &
                             coeff_x*(1.0 - coeff_y) * &
                             force_y /dx/dy
      drhoy( fi , fj+1 )   = drhoy( fi , fj+1 ) + &
                             (1.0 - coeff_x)*coeff_y * &
                             force_y /dx/dy
      drhoy( fi+1 , fj+1 ) = drhoy( fi+1 , fj+1 ) + &
                             coeff_x*coeff_y * &
                             force_y /dx/dy

    End do



                 

!  construct the density field by SOR

    nmax = 2*nx*ny


    !SOR法による補正圧力 phiの求解
    do m = 1, nmax

            !----------boundary condition of density----------!
            !----------上壁の境界条件の設定----------!
            rho(1:nx,ny+1) = rho(1:nx,ny)

            !----------下壁の境界条件の設定----------!
            rho(1:nx,0) = rho(1:nx,1)

            !----------右壁の境界条件の設定----------!
            rho(nx+1,1:ny) = rho(nx,1:ny)

            !----------左壁の境界条件の設定----------!
            rho(0,1:ny) = rho(1,1:ny)
    
            !------------- 四隅の密度 ---------------!
            rho(0   ,0   ) = rho(1 ,1 )
            rho(nx+1,0   ) = rho(nx ,1 )
            rho(0   ,ny+1) = rho(1 ,ny)
            rho(nx+1,ny+1) = rho(nx,ny)
    !----------boundary condition of density----------!

        error = 0.0             !誤差（収束判定用）のリセット

        do i = 1, nx
            do j = 1, ny
                !-> 残差 resid の計算　->　補正項 dphi の計算
                !-> 密度ジャンプの更新

                     div = (drhox(i-1,j)-drhox(i,j))*dx+(drhoy(i,j-1)-drhoy(i,j))*dy

                    resid = rho(i,j)

 
                rho(i,j) = (1.0-omega)*rho(i,j) + 0.25*omega * &
                           (rho(i+1,j)+rho(i-1,j)+rho(i,j+1)+rho(i,j-1)+div)

                !誤差：各補正項の相対誤差の最大値とする
                error = max(error,abs(resid-rho(i,j))/abs(resid))
    
            end do
        end do
        
        ! write(*,*)error
        if(error<eps) exit  !収束判定
    
    !----------boundary condition of density----------!
    !----------上壁の境界条件の設定----------!
    rho(1:nx,ny+1) = rho(1:nx,ny)

    !----------下壁の境界条件の設定----------!
    rho(1:nx,0) = rho(1:nx,1)

    !----------右壁の境界条件の設定----------!
    rho(nx+1,1:ny) = rho(nx,1:ny)

    !----------左壁の境界条件の設定----------!
    rho(0,1:ny) = rho(1,1:ny)
    
    !------------- 四隅の密度 ---------------!
    rho(0   ,0   ) = rho(1 ,1 )
    rho(nx+1,0   ) = rho(nx ,1 )
    rho(0   ,ny+1) = rho(1 ,ny)
    rho(nx+1,ny+1) = rho(nx,ny)
     !----------boundary condition of density----------!


    end do


   write(*,*)"distribute_density_field==================OK"
end subroutine distribute_density_field

subroutine distribute_viscosity_field(rho,myu)

use parameter
implicit none

    double precision rho(0:nx+1 , 0:ny+1)   !密度
    double precision myu(0:nx+1 , 0:ny+1)   !粘度

    integer          i , j

    do j = 1 , ny
        do i = 1 , nx

            myu(i,j) = myu1 + ( myu2 - myu1 ) * ( (rho(i,j)-rho1) / (rho2-rho1) )

        end do
    end do

    
   write(*,*)"distribute_viscosity_field==================OK"
end subroutine distribute_viscosity_field
