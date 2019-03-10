PROGRAM CODEHLL
  use fonction_HLL
  IMPLICIT NONE
  !
  integer::n,i,ng,nd , iter
  real*8:: dt,dx, tinit, tfinale,gamma,prL,prR
  real*8,dimension(:),allocatable::x,SR
  real*8,dimension(:,:),allocatable::w1,w0
  gamma = 1.4
  !
  ng = 50
  nd = 50
  n = nd + nd
  allocate(x(n+1))
  allocate(w1(3,n+1),w0(3,n+1))
  x(1)=0
  x(n) = 1.0
  dx = (x(n) - x(1))/n
  do i=2,n+1
     x(i) = x(i -1) + dx
  end do
!  print*, x

  w0 = VECT0(ng,nd)
  !
  w1 = 0
  ! temps initiale
  tinit = 0
  ! temps finale
  tfinale = 0.25
  iter = 0
  do while(tinit<tfinale)
     iter = iter + 1
     w1(:,1) = w0(:,1)
     do i =2,ng+nd        
        prL = pression(w0(:,i),gamma)
        prR = pression(w0(:,i+1),gamma)
        SR = SLR(w0(1,i),w0(1,i+1),w0(2,i)/w0(1,i),w0(2,i+1)/w0(1,i+1),prL,prR,gamma)
        dt=dx/2.36!SR(2)*0.5
        w1(:,i)=w0(:,i)-(dt/dx)*( FHLL(w0(:,i+1),w0(:,i),gamma)-FHLL(w0(:,i),w0(:,i-1),gamma))
     end do
     w1(:,size(x))=w0(:,size(x))-(dt/dx)*(FHLL(w0(:,size(x)),w0(:,size(x)-1),gamma)-FHLL(w0(:,size(x)-1),w0(:,size(x)-2),gamma))
     w0 = w1
     tinit = tinit + dt
     print*, tinit
  end do

  open(1,file='resultat_HLL.txt')
  do i=1,nd+ng+1
     write(1,*) x(i),W0(1,i),W0(2,i)/W0(1,i),(W0(3,i)-0.5*W0(2,i)*W0(2,i)/W0(1,i))/W0(1,i), pression(W0(:,i),gamma)
  end do
  close(1)
  print*,'nombre ditération = ',n
  deallocate(x,w1,w0)
  
!  deltaT = deltax/SR !Davis(1988)
   
END  PROGRAM CODEHLL
