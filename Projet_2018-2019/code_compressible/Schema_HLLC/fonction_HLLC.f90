Module fonction_HLLC

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

  FUNCTION SLR(rhoL,rhoR,uL,uR,PL,PR,gamma)
    real*8,intent(in)::UL,UR, PL,PR,rhoL,rhoR, gamma
    real*8::aL,aR,pstar,hL,hR,qL,qR
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
    SLR(1) = UL-aL*qL ! SL
    SLR(2) = uR+aR*qR ! SR
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

!Flux
  FUNCTION Flux(U,gamma)
    real*8,dimension(:),intent(in)::U
    real*8,dimension(3):: Flux
    real*8,intent(in)::gamma
    real*8::b
    b = 0.0
    Flux(1) = U(2)
    !Flux(2) = U(2)**2/U(1) + pression(U,gamma)
    !Flux(3) = U(2)/U(1)*(U(3) + pression(U,gamma))
    Flux(2) = 0.5*(3-gamma)*(U(2)**2/U(1)) + (gamma - 1)*U(3) !+ pression(U,gamma)
    Flux(3) = gamma*(U(2)/U(1))*U(3) - 0.5*(gamma - 1)*(U(2)**3/U(1)**2) !+ pression(U,gamma))
  end FUNCTION Flux


  FUNCTION S_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
    real*8,intent(in)::UL,UR,PL,PR,rhoL,rhoR,gamma
    real*8::S_star
    real*8,dimension(2)::SSLR
    SSLR = SLR(rhoL,rhoR,uL,uR,PL,PR,gamma)
    S_star = (rhoR*UR*(SSLR(2)-UR)-PR-rhoL*UL*(SSLR(1)-UL)+PL)/(rhoR*(SSLR(2)-UR)-rhoL*(SSLR(1)-UL))
  end FUNCTION S_star


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

    FUNCTION Flux_starL(rhoL,rhoR,uL,uR,PL,PR,gamma)
      real*8,intent(in)::rhoL,rhoR,uR,uL,gamma,PL,PR
      real*8,dimension(3):: Flux_starL
      real*8:: sstar,ELL
      real*8:: pstar
      real*8,dimension(2)::SSLR
      ELL = EL(rhoL,UL,PL,gamma)
      SSLR = SLR(rhoL,rhoR,uL,uR,PL,PR,gamma)
      pstar =pression_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
      sstar = S_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
      Flux_starL(1)=(rhoL*(SSLR(1)-ul)/(SSLR(1)-Sstar))*(sstar)
      Flux_starL(2)=(rhoL*(SSLR(1)-ul)/(SSLR(1)-Sstar))*sstar**2 + pstar
      Flux_starL(3)= sstar*( (ELL*(SSLR(1) -uL) - uL*PL + Sstar*pstar)/(SSLR(1) - Sstar) + PL)
    end FUNCTION Flux_starL

    FUNCTION Flux_starR(rhoL,rhoR,uL,uR,PL,PR,gamma)
      real*8,intent(in)::rhoL,rhoR,uR,uL,gamma,PL,PR
      real*8,dimension(3):: Flux_starR
      real*8:: sstar,ERR
      real*8::pstar
      real*8,dimension(2)::SSLR
      ERR = ER(rhoR,UR,PR,gamma)
      pstar = pression_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
      sstar = S_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
      SSLR = SLR(rhoL,rhoR,uL,uR,PL,PR,gamma)
      !
      Flux_starR(1)=(rhoR*(SSLR(2)-uR)/(SSLR(2)-Sstar))*(sstar)
      Flux_starR(2)=(rhoR*(SSLR(2)-uR)/(SSLR(2)-Sstar))*sstar**2 + pstar
      Flux_starR(3)=sstar*( (ERR*(SSLR(2) -uR) - uR*PR + Sstar*pstar)/(SSLR(2) - Sstar) + PL)
    end FUNCTION Flux_starR


      FUNCTION FHLLC(Uip,Ui,gamma)!rhoL,rhoR,uL,uR,PL,PR,gamma)
        real*8,dimension(:),intent(in) :: Uip, Ui
        real*8,intent(in)::gamma
        real*8::UL,UR, PL,PR,rhoL,rhoR
        real*8:: sstar,pstar
        real*8,dimension(3):: FHLLC
        real*8,dimension(2)::SSLR
        !
        rhoL = Ui(1); rhoR = Uip(1);
        UL = Ui(2)/Ui(1); UR = Uip(2)/Uip(1);
        PL = pression(Ui,gamma); PR = pression(Uip,gamma)
        !
        SSLR = SLR(rhoL,rhoR,uL,uR,PL,PR,gamma)
        sstar = S_star(rhoL,rhoR,uL,uR,PL,PR,gamma)
        !
        if(SSLR(1)>=0) FHLLC = Flux(Ui,gamma)
        if(SSLR(1)<=0 .and. sstar>=0) FHLLC = Flux_starL(rhoL,rhoR,uL,uR,PL,PR,gamma)
        if(sstar<=0 .and. SSLR(2)>=0) FHLLC = Flux_starR(rhoL,rhoR,uL,uR,PL,PR,gamma)
        if(SSLR(2)<=0) FHLLC = Flux(Uip,gamma)
      end FUNCTION FHLLC

    END MODULE fonction_HLLC
