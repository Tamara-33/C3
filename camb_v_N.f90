module tmatrix
implicit none 

INTEGER, PARAMETER :: I32=4 
INTEGER, PARAMETER :: dp=8
INTEGER, PARAMETER :: FAILED=1

REAL(dp), PARAMETER :: Zero=0._8

COMPLEX, PARAMETER ::  C0=(0._8,0._8)
COMPLEX, PARAMETER :: C1= (1._8,0._8)
COMPLEX, PARAMETER :: CI=(0._8,1._8)


REAL(dp), PARAMETER :: M_PI=3.14159265358979323846264338328_8     

REAL(dp), PARAMETER :: M_DPI=6.28318530717958647692528676656_8   

REAL(dp), PARAMETER :: M_PI_2=1.57079632679489661923132169164_8      

  
  real(dp),public :: ZT=1.065_8                               !>!< Initial and perturbation target charge
  real(dp),public :: ZTf=1._8                                 !>!< Final target charge
  real(dp),public :: ZPf=1._8                                 !>!< Final projectile charge 
  real(dp),public :: ZNf=1._8                                 !>!< Final projectile-target charge interaction
  real(dp),public :: MT=4*1836.152_8                          !>!< Target Mass   (proton mass=1836.152 au)
  real(dp),public :: MP=1*1836.152_8                          !>!< Projectile Mass  
  real(dp),public :: Me=1._8                                  !>!< Electron Mass
  real(dp),public :: m_P                                      !>!< projectile-electron reduced mass
  real(dp),public :: m_T                                      !>!< target-electron reduced mass
  real(dp),public :: m_N                                      !>!< projectile-target reduced mass
  integer(dp), parameter :: imax=28
  integer(dp), parameter :: ip=24
  real(dp), private :: q, theta

  ! Errores sugeridos: 0.005,0.002,0.001 o 0.002,0.001,0.0007
  real(dp), public :: er_q,er_theta,er_phi                   !>!< Precision en la integracion
  data er_q/0.005/
  data er_theta/0.002/
  data er_phi/0.001/

  real(dp), dimension(3) :: veck_P, veck_T, veck_N, vecQ
  real(dp) :: nu_P, nu_N, nu_T
  private
  public :: cdwb1

contains
  
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!!!!                                                                                    C3 approximation

!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  !> Evaluates the Transition Matrix in the CDW-B1 (C3) approximation
  !! @param veck electron momentum
  !! @param vecQQ projectile's transferred momentum \f$ P_{i} - P_{f} \f$
  !! @param vel incident projectile's velocity
  !! @return Tmatrix
  !! Tif= e^(i kj rj + Kj Rj)/(2PI)^3 \PROD Gamma()e^(-PI nuj/) F1()
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  function cdwb1(veck, vecQQ, vel) result(Tmatrix) 
    implicit none
    complex(dp) :: Tmatrix
    real(dp), dimension(3), intent(IN) :: veck
    real(dp), dimension(3), intent(IN) :: vecQQ
    real(dp), intent(IN) :: vel
    real(dp) :: k_P, k_T, k_N
    complex(dp) :: Coeff
    real(dp), parameter :: iDPI3=1/(M_DPI**3)

    !<!> Masas reducidas
 
    m_P=MP*Me/(MP + Me)        
    m_T= MT*Me/(MT + Me)         
    m_N= MP*MT/(MP + MT)   

    vecQ= vecQQ
    veck_T= veck - (m_T/MT)*vecQ
    veck_P= m_P*(veck + vecQ/MP)
    veck_N= m_N*veck/MT - vecQ
    veck_P(3)= veck_P(3) - m_P*vel
    veck_N(3)= veck_N(3) + m_N*vel
    k_P= SQRT(dot_product(veck_P, veck_P));  k_T= SQRT(dot_product(veck_T, veck_T))
    k_N= SQRT(dot_product(veck_N, veck_N))
    nu_T = -m_T*ZTf/k_T  ;  nu_P = -m_P*ZPf/k_P   ;  nu_N = m_N*ZPf*ZNf/k_N

    !coeff= coufa(nu_P)*coufa(nu_N)*coufa(nu_T) !Coulomb Factors  
    coeff=1.0; 
    Tmatrix= iDPI3 * coeff * H()      !PLANTEA LA FORMA GENERAL DEL ELEMNTO DE MATRIZ 
  end function cdwb1



  !> Performs the 3D convolution to get the T-Matrix in C3 approximation
  !! @return T

!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!!!!  Auxiliary functions
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  !> Coulomb factor \f$ e^{-\pi \, \nu/2} \Gamma(1 + i \nu) \f$
  !! @param nu Sommerfeld parameter 
  !! 
  !! @return A complex number
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  function coufa(nu) result(CF)
    implicit none
    real(8), intent(IN) :: nu
    complex(8) :: CF
    CF= -nu*M_PI_2 + clgama(C1 + CI*nu)
    CF= exp(CF)
  end function coufa

!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  !> Logarithm of the gamma function
  !! @param CZ argument
  !! 
  !! REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATICAL FUNCTIONS'.
  !! DOVER, NEW YORK (1974). PP 255-257.
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  function clgama(CZ) result(LG)
    implicit none
    complex(8), intent(IN) :: CZ
    complex(8) :: LG
    real(8) :: AR, AI
    complex(8) :: CZA,CZFAC,CZFL, CZS, CZI2
    integer(4) :: iconj
    CZA= CZ
    AR= CZA
    LG= 36.84136149D0   !¿porque usa el logaritmo de este numero de donde salio?????? 
    IF(abs(CZA) < 1.d-16) return
    !
    AI= -CZA*CI
    if(AI > Zero) then
      iconj=0
    else
      iconj=1
      CZA= conjg(CZA)
    endif

    CZFAC= C1
    CZFL=  C0
    do   ! Recursion
      CZFAC=CZFAC/CZA
      if(abs(CZFAC) > 1.0d8) then
        CZFL=CZFL + log(CZFAC)
        CZFAC= C1
      endif
      CZA= CZA + 1._8
      AR= real(CZA, kind=8)
      IF(abs(CZA) < 1.d-16) return   ! Polo
      IF((abs(CZA) > 15) .and. (AR > Zero)) exit
    enddo
    ! Stirling's expansion of log(GAMMA(CZA)). [Abram 6.1.41] (para aproximación de factoriales grandes)
    CZI2=1.0D0/(CZA*CZA)
    CZS=(43867.0D0/244188.0D0)*CZI2
    CZS=(CZS-3617.0D0/122400.0D0)*CZI2
    CZS=(CZS+1.0D0/156.0D0)*CZI2
    CZS=(CZS-691.0D0/360360.0D0)*CZI2
    CZS=(CZS+1.0D0/1188.0D0)*CZI2
    CZS=(CZS-1.0D0/1680.0D0)*CZI2
    CZS=(CZS+1.0D0/1260.0D0)*CZI2
    CZS=(CZS-1.0D0/360.0D0)*CZI2
    CZS=(CZS+1.0D0/12.0D0)/CZA
    LG=(CZA-0.5D0)*log(CZA)-CZA + 9.1893853320467274D-1+CZS +CZFL+log(CZFAC)
    IF(iconj == 1) LG= conjg(LG)
  end function clgama
  





  function H() result(T)
    implicit none
    complex(dp) :: T
    real(dp) :: QM, DQ, error
    integer(I32) :: n
    QM= SQRT(dot_product(vecQ,vecQ))
    Dq= 1.d0 + QM               ! Step used in the integration
    call adsi_ZI(Zero, Dq, F3d_phi_theta_q, er_q, T, error, n)
  end function H

  function F3d_phi_theta_q(qq) result(y)
    implicit none
    real(dp), intent(IN) :: qq
    complex(dp) :: y
    integer n
    real(dp) :: error
    q=qq
    !     integracion en theta:
    call adsi_Z(Zero,M_PI, F3d_phi_theta, er_theta, y,error, n)

    IF(error > er_theta) print *,'# warning: theta integ', error, '>', er_theta

  end function F3d_phi_theta_q

  function F3d_phi_theta(ttheta) result(y)
    implicit none
    real(dp),intent(IN) :: ttheta
    complex(dp) :: y
    real(dp) :: error
    integer(I32) :: n
    theta=ttheta
    !     integracion en phi
    call adsi_Z(Zero, M_DPI,F3d_phi,er_phi, y, error,n)

    IF(error > er_phi) print *,'# warning: phi integ', error, '>', er_phi

  end function F3d_phi_theta

!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  !> Integrand as a function of phi_q, for fixed q, theta_q.
  !! The expression is:
  !! \f{align*}{
  !!   t_{if}&= \frac{1}{(2 \pi)^{9/2}} \int d \bm{p}\, F_{if}(\bm{p})\,\Big[ Z_{N} J_{1}(0^{+}, \bm{k}_{1},-\nu_{N}, \bm{k}_{N}) J_{0}(0^{+}, \bm{k}_{2},-\nu_{P}, \bm{k}_{P}) - Z_{P}\,  J_{0}(0^{+}, \bm{k}_{1},-\nu_{N}, \bm{k}_{N}) J_{1}(0^{+}, \bm{k}_{2},-\nu_{P}, \bm{k}_{P}) \Big]
  !!   \intertext{with}
  !!   \bm{k}_{1}&= \frac{m_{T}}{m}\bm{Q} - \bm{p} & \bm{k}_{2}&= - \left( \frac{m_{T}}{M_{T}}\bm{Q} + \bm{p}  \right)
  !!   \f}
  !! @param phi 
  !! @return 
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  function F3d_phi(phi) result(CInteg)
    implicit none
    real(dp), intent(IN) :: phi
    complex(dp) :: CInteg

    real(dp) :: lambda= 1.e-4_8  !>!< Convergence factor
    real(dp) :: s
    real(dp), dimension(3) :: q_v
    real(dp), dimension(3) :: veck1, veck2 !!, veck3
    complex(dp) :: J0_P, J0_N, J1_P, J1_N
    complex(dp) :: Fif
    !
    s= sin(theta)
    q_v(1)=q*s*cos(phi)
    q_v(2)=q*s*sin(phi)
    q_v(3)=q*cos(theta)
    veck1= (m_T/Me)*vecQ - q_v
    veck2= -(m_T/MT)*vecQ - q_v
    call J_01_nordsieck(lambda,veck1, -nu_N, veck_N, J0_N, J1_N)
    call J_01_nordsieck(lambda,veck2, -nu_P, veck_P, J0_P, J1_P)

    Fif= J_0_nordsieck(ZT, q_v - veck_T, -nu_T, veck_T) ! Form-factor
    Fif= Fif*(2*ZT)**1.5_8/(M_DPI*M_DPI)

    Cinteg= C0
    CInteg= Cinteg + ZPf*ZNf*J0_P*J1_N*Fif ! T_N
    CInteg= Cinteg - ZPf*J1_P*J0_N*Fif     ! T_P
    CInteg= CInteg*s*q**2
  end function F3d_phi





!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  !> Routine of integration
  !! res= \f$ \int_{xl}^{\infty}  fun(x) dx \f$
  !!
  !! It is suspected that the main contribution comes from the region
  !! FUN(xl) to FUN(xl+DX), afterwards moves by adding ranges of DX
  !! @param xl inferior limit
  !! @param dx 
  !! @param fun a complex valued function to integrate
  !! @param er relative error, suggested values 1.e-3,1.e-2,
  !! @param res Value of integral
  !! @param overall_error 
  !! @param nout If positive: total number of evaluations
  !! @return nout, res
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  subroutine adsi_ZI(xl, dx, fun, er, res, overall_error, nout) 
    implicit none
    real(dp), intent(IN) :: xl
    real(dp), intent(IN) :: dx
    complex(dp) :: fun
    real(dp), intent(IN) :: er
    real(dp), intent(OUT) :: overall_error
    complex(dp), intent(INOUT) :: res
    integer(I32), intent(OUT) :: nout
    external fun

    integer(4) :: in
    real(8) :: xu, lambda
    real(8) :: xu_i,xl_i,ERI,RATIO, error
    complex(8) :: R_i
    lambda=1.d0
    xu= xl + lambda*DX   ! Primera integracion : [xl,xl+DX]
    Res= C0;     R_i= C0
    call adsi_Z(xl,xu, fun, ER, R_i, Overall_error, nout)
    IF(nout < 0) return 
    Res= Res + R_i

    xu_i= xu    ;   ERI=ER    ! Sucesivos pasos de integracion
    do nout= 1, imax
      R_i= C0 ;   xl_i= xu_i;    xu_i= xl_i+ lambda*DX
      call adsi_Z(xl_i,xu_i,FUN,ERI,R_i,error, in)
      IF(in < 0) return
      Res= Res + R_i 
      RATIO=abs(R_i)/abs(Res)
      Overall_error= Overall_error + error*RATIO
      IF(RATIO < ER) return
      if(RATIO < 0.1d0) then
        ERI=ER/RATIO ! Relaja la condicion
      else
        lambda=max(1.5*lambda, 10.d0) ! Extends the new region
        ERI= ERI*1.5
      end if
      ERI= min(ERI,0.1d0) 
    enddo
    nout= FAILED !  No convergence was obtained
  end subroutine adsi_ZI

!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  !> Integration by Adaptive Simpson method of a function
  !! \f[ \int_{a}^{b} f(t)\, dt \f]
  !! @param xL inferior limit
  !! @param xU superior limit
  !! @param FUN \f$f(t)\f$: Function to integrate
  !! @param EPS Bound to estimation error desired
  !! @param[out] Res result
  !! @param[out] Error Estimated error
  !! @return  El numero de puntos usados es 3+2*ncalls
  !! @note Usa metodo de integracion por Simpson con 3 puntos
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  subroutine adsi_Z(xL,xU,FUN,EPS,Res,Error, ncalls)
    implicit none
    real(8), intent(IN) :: xL
    real(8), intent(IN) :: xU
    real(8), intent(IN) :: EPS
    real(8), intent(OUT) :: Error
    complex(8), intent(INOUT) :: Res
    complex(8) :: FUN
    integer(4), intent(OUT) :: ncalls
    integer(4) :: k,j
    real(8) :: XM
    complex(8) :: CR,CRR, Fswp
    real(8), dimension (imax) :: x_L,x_U
    complex(8), dimension(imax) :: FxL,FXM,FxU
    complex(8), dimension(imax) :: CINT

    !  ----VALORES INICIALES
    ncalls=0 ;      j=1
    Res= C0
    x_L(j)=xL            ;     FxL(j)=FUN(xL)
    x_U(j)=xU            ;     FxU(j)=FUN(xU)
    XM=.5d0*(xU+xL)      ;     FXM(j)=FUN(XM)
    CINT(j)= ((x_U(j)-x_L(j))/6.d0)* (FxL(j) + FxU(j) +4*FXM(j)) !Simpson
    do ! ---- comienzo de la iteracion
      ncalls= ncalls + 1
      k=j ;      j=j+1
      CR=CINT(k)

      x_U(j)= x_U(k)                  ;     FxU(j)= FxU(k)
      x_U(k)= 0.5d0*(x_L(k) + x_U(k)) ;     FxU(k)= FXM(k)
      x_L(j)= x_U(k)                  ;     FxL(j)= FXM(k)
      XM=  0.5d0*( x_U(k)+x_L(k) )    ;     FXM(k)= FUN(XM)
      CINT(k)= ((x_U(k)-x_L(k))/6.d0)* (FxL(k) + FxU(k) +4*FXM(k)) !Simpson
      XM= 0.5d0*( x_U(j)+x_L(j) )     ;     FXM(j)= FUN(XM)
      CINT(j)= ((x_U(j)-x_L(j))/6.d0)* (FxL(j) + FxU(j) +4*FXM(j)) !Simpson

      ! Elegimos para dividir el intervalo que mas contribuye a la integral
      if( abs(CINT(k)) > abs(CINT(j)) ) then ! SWAP k <--> j
        XM= x_U(j);     x_U(j)= x_U(k);   x_U(k)= XM
        Fswp= FxU(j);   FxU(j)= FxU(k);   FxU(k)= Fswp
        XM= x_L(j);     x_L(j)= x_L(k);   x_L(k)= XM
        Fswp= FxL(j);   FxL(j)= FxL(k);   FxL(k)= Fswp
        Fswp= FXM(j);   FXM(j)= FXM(k);   FXM(k)= Fswp
        Fswp= CINT(j);  CINT(j)=CINT(k);  CINT(k)=Fswp
      end if

      CRR=CINT(k)+CINT(j)
      Error= abs(CR-CRR)
      if( Error < (EPS*abs(Res+CRR)) )  then
        !  Fin de la iteracion en k para un subintervalo
        Res= Res + CRR
        Error= Error/abs(Res)
        j= j-2
        IF(j <= 0) return
      else
        if(j+1 > imax) then
          ncalls= FAILED ! This error message should mean that the obtained value 
          return   ! is meaningless. We should set some flags with values of X,INT
        end if
      end if
    enddo
  end subroutine adsi_Z
  
  !> \copydoc J_0_nordsieck
  !! @param Z 
  !! @param Q 
  !! @param A1 
  !! @param P1 
  !! @param J0 
  !! @param J1 
  !! 
  !! @return 
  subroutine J_01_nordsieck(Z,Q,A1,P1,J0,J1)
    ! EVALUATION OF J_0 AND J_1, (SEE BELOW) J_0 AND J_1 SUBROUTINES
    implicit none
    real(8), intent(IN) :: Z, A1
    real(8), dimension(3), intent(IN) :: Q, P1
    complex(8), intent(OUT) :: J0,J1
    real(8) :: P1M,QM2,D
    complex(8) :: CIA1,CIZ,u1,u1m1,CC1

    CIZ=CI*Z
    CIA1=CI*A1
    P1M= SQRT(dot_product(P1,P1))
    QM2= dot_product(Q,Q)
    D=2/(Z*Z + QM2)
    u1=D*(dot_product(P1,Q) - CIZ*P1M)
    u1m1= C1 + u1
    CC1= P1M - CIZ*u1

    J1= M_DPI*D*exp(-CIA1*log(u1m1))
    J0= D*(Z+A1*CC1/u1m1)*J1
  end subroutine J_01_nordsieck

!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  !> Calculates the analytical form of the simple Nordsieck integral
  !! \f[
  !! J_{0}( Z, \bm{p}, a_{1}, \bm{k}_{1} ) =
  !!\int d \bm{r} \; e^{i \bm{p} \cdot \bm{r}} \; e^{- Z\, r} \, {_{1}F_{1}} \left( i
  !!a_{1} ; 1 ; i \left(k r + \bm{k}_{1}\cdot \bm{r} \right)\right)
  !! \f]
  !! @param Z  real coefficient in the exponential
  !! @param p  real vector in the imaginary exponential
  !! @param a1 Sommerfeld parameter appearing in the hypergeometric 
  !! @param k1 momentum, argumento in the hypergeometric
  !! @return An scalar complex with the value of the integral
  !! @details The result of the integral is
  !! \f[
  !! J_{0}( Z, \bm{p}, a_{1}, \bm{k}_{1} ) = \left( \frac{1}{D}
  !!\frac{d D}{d Z} + \frac{i a_{1}}{1 + u_{1}} \frac{d u_{1}}{d Z} \right) J_{1}( Z,
  !!\bm{p}, a_{1}, \bm{k}_{1} )
  !! \f]
  !! with
  !! \f{eqnarray*}{
  !! u_{1} &=& 2 \left( \bm{k}_{1}\cdot \bm{p} - i Z k_{1} \right)/D \cr
  !! D &=& p^{2} + Z^{2} 
  !! \f}
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  function J_0_nordsieck(Z,p,a1,k1) result(J0)
    implicit none
    real(8), intent(IN) :: Z, a1
    real(8), dimension(3), intent(IN) :: p
    real(8), dimension(3), intent(IN) :: k1
    complex(8) :: J0
    real(8) :: P1M,QM2,D
    complex(8) :: CIZ
    complex(8) :: u1,u1m1,CC1

    CIZ=CI*Z
    P1M= SQRT(dot_product(k1,k1))
    QM2= dot_product(p,p)
    D= 2/(Z*Z + QM2)
    u1= (dot_product(k1,p) - CIZ*P1M )*D
    u1m1= C1 + u1
    CC1= P1M - CIZ*u1
    J0= M_DPI*(D**2)*(Z+a1*CC1/u1m1)*exp(-CI*a1*log(u1m1))
  end function J_0_nordsieck
  
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  !> Analytical form of:
  !! \f[ J_{1}( Z, \bm{Q}, a_{1}, \bm{k}_{1} ) = \int d \bm{r} \; e^{i \bm{p} \cdot
  !! \bm{r}} \; \frac{e^{- Z\, r}}{r} {_{1}F_{1}} \left( i a_{1} ; 1 ; i \left(k r +
  !! \bm{k}_{1}\cdot \bm{r} \right)\right) = \frac{4 \pi}{D} (1+a_{1})^{- i a_{1}}
  !! \f]
  !! @param Z  real coefficient in the exponential
  !! @param p  real vector in the imaginary exponential
  !! @param a1 Sommerfeld parameter appearing in the hypergeometric 
  !! @param k1 momentum, argumento in the hypergeometric
  !! @return An scalar complex with the value of the integral
!<!>--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  function J_1_nordsieck(Z,p,a1,k1) result(J1)
    implicit none
    real(8), intent(IN) :: Z, a1
    real(8), dimension(3), intent(IN) :: p, k1
    complex(8) :: J1
    real(8) :: P1M,QM2,D
    complex(8) :: u1

    P1M= SQRT(dot_product(k1,k1))
    QM2= dot_product(p,p)
    D=2/(Z*Z + QM2)
    u1= (dot_product(k1,p) - CI*Z*P1M )*D
    J1=(M_DPI*D)*exp(-CI*a1*log(C1 + u1))
  end function J_1_nordsieck

end module tmatrix


!<!><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
program tdcs
use tmatrix
implicit none 


  integer, parameter :: out_tif=20                                                 !>!< maximum length of the file name 
  real(8) :: E_proj                                                               !>!< Projectile energy in a.u
  real(8) :: vel                                                                  !>!< initial projectile velocity
  real(8) :: thp                                                                  !>!< final projectile angle polar
  real(8) :: the                                                                  !>!< final ejection electron angle polar
  real(8) :: phip                                                                 !>!< Angle azimutal of proyectil
  real(8) :: phie                                                                 !>!< Angle azimutal of electron
   real(8) :: phi                                                                 !>!< Angle between projectil and electron
  real(8),dimension(3) :: veck       	                                           !>!> Final momentum of the electron
  real(8),dimension(3) :: vecP     	                                               !>!> Final momentum of the proton (proyectil)
  real(8),dimension(3) :: vecQQ		                                               !>!> Momentum transfer
  real(8) :: K                                                                    !>!< modulo del momento del proyectil
  real(8) :: ke                                                                   !>!< Modulo del momento del electron
  real(8) :: Elig
  real(8) :: pN2, pk, Coef
  real(8), dimension(3) :: K_N, p, K_v
  complex(8) :: T
  real(8) :: tr, tred, traux
  integer :: i,j
 real(8) ::  mu_N, mu_T

  !<------------------------------------------------------------------------------------------------------------------------------------
  !            Fijamos datos del sistema de colisión (proyectil, target, etc)
  !>------------------------------------------------------------------------------------------------------------------------------------
  E_proj = (75000_8)*0.036749_8                                                      !>!< Energía del proyectil es 75 keV.
  Elig = (24.49046233_8)*0.036749_8                                                  !>!< Energía de ligadura del átomo de He (CAMBIE LA ENERGIA DEL ENLACE DEL HELIO)
  phip = 0.0_8                                                                   !>!< Angle azimutal del proyectil (protón)
  phie=0.0_8                                                                     !>!< Angle azimutal del electron
  the= 0.0_8                                                                     !>!< Angle polar del proyectil 
  phi= 0.0005
  thp= 0.0005


!	open (unit=out_tif,file="100.dat",action="write",status="replace")

  mu_N = (MP + MT) * Me / (MP + MT + Me)                                           !< Masa reducida (P+T)-e
  mu_T = (MT + Me) * MP / (MP + MT + Me)                                           !< Masa reducida (T+e)-P
  m_T= MT*Me/(MT + Me)    
  m_N= MP*MT/(MP + MT) 

  vel = sqrt(2 * E_proj / MP)                                                         !< Velocidad inicial del proyectil

   open (unit=out_tif,file="salidita.dat",action="write",status="replace")


  
  !<-----------------------------------------------------------------------------------------------------------------------------------------
  !< Momento transferido
  !!>-----------------------------------------------------------------------------------------------------------------------------------------
    
	do i=0,20
    
    	do j=0,20

         vecQQ(3)=4.7D0+254.03D0*(dble(i))                                             !>!< z-component of momentum transfer
         vecQQ(1)=(dsqrt((2540.3D0**2)-(vecQQ(3)-2545D0)**2))*(dble(j)/20D0)           !>!< x-component of momentum transfer
         vecQQ(2)=0.0                                                                  !>!< y-component of momentum transfer

  !<-----------------------------------------------------------------------------------------------------------------------------------------
  !< Momento del proyectil K
  !!>-----------------------------------------------------------------------------------------------------------------------------------------
         

         vecP(1)= -vecQQ(1)                                                            !>!< x-component of proton momentum
         vecP(2)= -vecQQ(2)                                                            !>!< y-component of proton momentum
         vecP(3)= MP*vel - vecQQ(3)                                                    !>!< z-component of proton momentum
         
         K=sqrt(vecP(1)**2 + vecP(2)**2 + vecP(3)**2 )


!<-----------------------------------------------------------------------------------------------------------------------------------------
  !< Momento del electrón k
  !!>-----------------------------------------------------------------------------------------------------------------------------------------
           
      ke=  (MP*m_T*vel*cos(the))/MT - (K*m_T*cos(the)*cos(thp))/MT - (K*m_T*cos(phi)*cos(phip)*sin(the)*sin(thp))/MT +& 
       m_T*Sqrt(((MP*vel*cos(the))/MT - (K*cos(the)*cos(thp))/MT - (K*cos(phi)*cos(phip)*sin(the)*sin(thp))/MT)**2 +& 
          (2*(-Elig + (MP*vel**2)/2. - (MP**2*vel**2)/(2.*MT) + (K*MP*vel*cos(thp))/MT - (K**2*cos(thp)**2)/(2.*m_N) -&
         (K**2*cos(phip)**2*sin(thp)**2)/(2.*m_N) - (K**2*sin(phip)**2*sin(thp)**2)/(2.*m_N)))/m_T)

    
        veck(1)= ke*sin(the)*cos(phie)                                                   !>!< x-component of electron momentum
        veck(2)= ke*sin(the)*sin(phie)                                                   !>!< y-component of electron momentum
        veck(3)= ke*cos(the)                                                             !>!< z-component of electron momentum

       T = cdwb1(veck, vecQQ , vel) 
       !write (out_tif,*) log(abs(T)**2)
       print *, log(abs(T)**2)
		 end do  

	end do  


 	close (out_tif)

 end program tdcs

