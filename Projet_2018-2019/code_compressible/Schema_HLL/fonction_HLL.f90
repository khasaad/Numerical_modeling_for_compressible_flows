Module fonction_HLL

  IMPLICIT NONE

CONTAINS

  FUNCTION pression(U,gamma)
    real*8 :: pression,b
    real*8,dimension(:), intent(in)::U
    real*8,intent(in):: gamma
    !real*8 :: gamma,  b
    !gamma = 1.4!; b = 0.0
    b = 0.0
    pression = (U(3) - 0.5*U(2)**2/U(1))*(gamma-1.0)/(1-b*U(1))
  END FUNCTION pression

  FUNCTION pression_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
    real*8,intent(in)::rhoL,rhoR,uR,uL,gamma,PL,PR
    real*8:: aL,aR,a_bar,rho_bar
    real*8:: pression_star
    aL = sqrt((gamma*PL)/rhoL)
    aR = sqrt((gamma*PR)/rhoR)
    a_bar = 0.5*(aL + aR)
    rho_bar = 0.5*(rhoL + rhoR)
    pression_star = 0.5*(PL + PR) - 0.5*(uR - uL)*rho_bar*a_bar
  END FUNCTION pression_star

  !FUNCTION pression(U)
  !  real*8 :: pression
  !  real*8,dimension(:), intent(in)::U
  !  real*8 :: gamma,  b
  !  gamma = 1.4; b = 0.0
  !  pression = (U(3) - 0.5*U(2)**2/U(1))*(gamma-1.0)/(1-b*U(1))
  !END FUNCTION pression

!Flux
FUNCTION SLR(rhoL,rhoR,uL,uR,PL,PR,gamma)
  real*8,intent(in)::UL,UR, PL,PR,rhoL,rhoR, gamma
  real*8::aL,aR,pstar,hL,hR,qL,qR,ustar,astar
  real*8,dimension(2)::SLR
  !
  aL = sqrt(gamma*pL/rhoL)
  aR = sqrt(gamma*pR/rhoR)
  !
  pstar = pression_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
  hL = pstar/pL
  hR = pstar/pR

  if (hL < 1.0) then
     qL = 1.0
  else
     qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(hL-1.0))
  end if
  if (hR < 1.0) then
     qR = 1.0
  else
     qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(hR-1.0))
  end if
  !
  !  SLR(1) = UL-aL*qL ! SL
  !  SLR(2) = uR+aR*qR ! SR
  ustar = 0.5*(UL-UR)+(aL-aR)/(gamma-1.0)
  astar = 0.5*(aL+aR)+0.25*(gamma-1.0)*(UL-UR)
  SLR(1) = minval([uL-aL,ustar-astar])
  SLR(2) = maxval([uR+aR,ustar+astar])
end FUNCTION SLR

FUNCTION EL(rhoL,UL,PL,gamma)
  real*8,intent(in)::gamma,UL,PL,rhoL
  real*8::EL
  !b = 0.0
  EL = 0.5*rhoL*UL**2 + PL/(gamma-1)
END FUNCTION EL

FUNCTION ER(rhoR,UR,PR,gamma)
  real*8,intent(in)::gamma,UR,PR,rhoR
  real*8::ER
  ER = 0.5*rhoR*UR**2 + PR/(gamma-1)
END FUNCTION ER


  FUNCTION Flux(U,gamma)
    real*8,dimension(:), intent(in)::U
    real*8,dimension(3):: Flux
    real*8,intent(in)::gamma
    real*8::b
    b = 0.0
    Flux(1)=U(2)
    Flux(2) = 0.5*(3-gamma)*(U(2)**2/U(1)) + (gamma - 1)*U(3)
    Flux(3) = gamma*(U(2)/U(1))*U(3) - 0.5*(gamma - 1)*(U(2)**3/U(1)**2)
    !Flux(2)= U(2)**2/U(1) + pression(U)
    !Flux(3)=U(2)/U(1)*(U(3) + pression(U))

  end FUNCTION Flux

  FUNCTION S_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
    real*8,intent(in)::UL,UR,PL,PR,rhoL,rhoR,gamma
    real*8::S_star
    real*8,dimension(2)::SSLR
    SSLR = SLR(rhoL,rhoR,uL,uR,PL,PR,gamma)
    S_star = (rhoR*UR*(SSLR(2)-UR)-PR-rhoL*UL*(SSLR(1)-UL)+PL)/(rhoR*(SSLR(2)-UR)-rhoL*(SSLR(1)-UL))
  end FUNCTION S_star

! vecteur vect=[rho, rho*u, E]^T
!  FUNCTION vect( p, ue, rho, gamma, b)
!    real*8,dimension(3)::vect
!    real*8,intent(in):: p, ue, rho, gamma,  b
!    real*8::q, ee, E
!    q = rho*ue
!    ee = p*(1-b*p)/((gamma - 1)*rho)
!    E=rho*(0.5*ue**2 + ee)
!    vect(1)= rho
!    vect(2)= q
!    vect(3) = E
!  end FUNCTION  VECT

  FUNCTION VECT0(ng,nd)
    !real*8,dimension(:),intent(in)::x
    integer,intent(in) :: ng,nd
    real*8,dimension(3,ng+nd+1)::VECT0
    real*8 :: p_left, u_left, rho_left
    real*8 :: p_right, u_right, rho_right
    real*8::q_left, ee_left, E_left
    real*8::q_right, ee_right, E_right
    real*8 :: gamma, b
    gamma = 1.4; b = 0.0
    rho_left = 1.0; u_left = 0.0; p_left = 1.0
    rho_right = 0.125; u_right = 0.0; p_right = 0.1
    !
    q_left = rho_left*u_left
    ee_left = p_left*(1-b*p_left)/((gamma - 1)*rho_left)
    E_left=rho_left*(0.5*u_left**2 + ee_left)
    !
    q_right = rho_right*u_right
    ee_right = p_right*(1-b*p_right)/((gamma - 1)*rho_right)
    E_right=rho_right*(0.5*u_right**2 + ee_right)
    !
    VECT0(1,1:ng) = rho_left
    VECT0(1,ng+1:ng+nd+1) = rho_right
    !
    VECT0(2,1:ng) = q_left
    VECT0(2,ng+1:ng+nd+1) = q_right
    !
    VECT0(3,1:ng) = E_left
    VECT0(3,ng+1:ng+nd+1) = E_right
  END FUNCTION VECT0

! valeur propre
  !!FUNCTION valp(U)
  !  real*8,dimension(:),intent(in)::U
  !  real*8,dimension(2):: valp
  !  real*8::  gamma
  !  real*8 :: aa
  !  gamma = 1.4
  !  aa = sqrt((gamma*pression(U))/U(1))
    !print*, 'aa=', pression(U)
  !  valp(1) = U(2)/U(1) - aa
  !  valp(2) = U(2)/U(1) + aa
  !end FUNCTION valp

!approximation de solveur de Riemann HLL
!  FUNCTION UHLL(UL,UR, SL, SR, p, gamma)
!    real*8,dimension(3):: UHLL
!    real*8,intent(in):: SL, SR, p, gamma
!    real*8,dimension(:),intent(in):: UL, UR
!    real*8,dimension(2):: tpr1, tpr2
!    tpr1=valp(UL,p,gamma)
!    tpr2=valp(UR,p,gamma)
!    if(SL>tpr1(1))                     UHLL = UL
!    if(SL<=tpr1(1) .and. SR>= tpr2(2)) UHLL = (SR*UR - SL*UL + Flux(UL,p) - Flux(UR,p))/(SR - SL)
!    if(SR<tpr2(2))                     UHLL = UR
!
!  END FUNCTION UHLL

! Flux HLL
  FUNCTION FHLL(Uip,Ui,gamma)
    !real*8 :: SL, SR
    !real*8,dimension(:),intent(in) :: UL, UR, vp
    real*8,dimension(3):: FHLL
    real*8,dimension(:),intent(in) :: Uip, Ui
    real*8,intent(in)::gamma
    real*8::UL,UR, PL,PR,rhoL,rhoR
    !real*8:: sstar,pstar
    real*8,dimension(2)::SSLR
    rhoL = Ui(1); rhoR = Uip(1);
    UL = Ui(2)/Ui(1); UR = Uip(2)/Uip(1);
    PL = pression(Ui,gamma); PR = pression(Uip,gamma)
    !
    SSLR = SLR(rhoL,rhoR,uL,uR,PL,PR,gamma)
    !real*8,dimension(3):: Utpr
    !Utpr = UHLL(UL, UR, SL, SR, p, gamma )
    !SL = minval(vp)-0.3
    !SR = maxval(vp)-0.3
    !
    if(SSLR(1)>=0) FHLL = Flux(Ui,gamma)
    if(SSLR(1)<=0 .and. SSLR(2)>=0) FHLL = (SSLR(2)*Flux(Ui,gamma) - SSLR(1)*Flux(Uip,gamma) + &
         SSLR(1)*SSLR(2)*(Uip - Ui))/(SSLR(2) - SSLR(1))
    if(SSLR(2)<=0) FHLL = Flux(Uip,gamma)
    !if(SL>=0) FHLL = Flux(UL)
    !if(SL<=0 .and. SR>=0) FHLL = (SR*Flux(UL) - SL*Flux(UR) + SL*SR*(UL - UR))/(SR - SL)
    !if(SR<=0) FHLL = Flux(UR)
  end FUNCTION FHLL

END MODULE fonction_HLL
