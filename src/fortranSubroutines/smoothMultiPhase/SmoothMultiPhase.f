c   ***************************
c   Multi-Phase (Smooth Model)
c   ***************************

c$$$      program SmoothMultiPhaseModel
c$$$
c$$$      implicit none
c$$$
c$$$      real*8 EyeTens(3,3),Dt,K,mu,a0,b0,a1,b1,m,kappas,Jn,kappan
c$$$      real*8 Beprn(3,3),Fr(3,3),J,Bepr(3,3),kappa,TGpa(3,3)
c$$$      integer i1,i2
c$$$
c$$$      call IdentityTens(EyeTens)
c$$$
c$$$      K=1.60d0*(10.0d0**9.0d0)
c$$$      mu=1.0d0*(10.0d0**9.0d0)
c$$$      Dt=0.10d0
c$$$
c$$$      a0=1.0d0
c$$$      a1=2.0d0
c$$$      b0=3.0d0
c$$$      b1=20.0d0
c$$$      m=0.10d0
c$$$
c$$$      Jn=1.0d0
c$$$      kappan=0.010d0
c$$$      kappas=0.020d0
c$$$
c$$$      do i1=1,3
c$$$        do i2=1,3
c$$$          Beprn(i1,i2)=EyeTens(i1,i2)
c$$$        end do
c$$$      end do
c$$$
c$$$      Fr(1,1)=1.1051709180756480d0
c$$$      Fr(1,2)=0.0d0
c$$$      Fr(1,3)=0.0d0
c$$$      Fr(2,1)=0.0d0
c$$$      Fr(2,2)=0.9512294245007140d0
c$$$      Fr(2,3)=0.0d0
c$$$      Fr(3,1)=0.0d0
c$$$      Fr(3,2)=0.0d0
c$$$      Fr(3,3)=0.9512294245007140d0
c$$$
c$$$      call SmoothMultiPhase(Dt,K,mu,a0,b0,a1,b1,m,kappas,Jn,Beprn,
c$$$     # kappan,Fr,J,Bepr,kappa,TGpa)
c$$$
c$$$      end program

c   ***************************************************************************
c   ***************************************************************************

c   ************
c   SUBROUTINES
c   ************

c   Calculate J=J(t_n)*Jr

      subroutine ElasticTrialDil(Jr,J1,Jes)

      implicit none

c   Declarations
      real*8 J1,Jr,Jes

      Jes=Jr*J1

      end subroutine

c   ================
c   Main Subroutine
c   ================

      subroutine SmoothMultiPhase(Dt,K,mu,a0,b0,a1,b1,m,kappas,Jn,
     # Beprn,kappan,Fr,J,Bepr,kappa,TGpa)

      implicit none

      real*8  EyeTens(3,3),Dt,K,mu,a0,b0,a1,b1,m,kappas,Jn,kappan
      real*8  Beprn(3,3),Fr(3,3),FrT(3,3),Frpr(3,3),Jr,J,Beprs(3,3)
      real*8  devBeprs(3,3),hydBeprs(3,3),devgeprs(3,3),games
      real*8  Dbar(3,3),Br(3,3),Deps,c0,c1,c2,gam
      real*8  GAMMA,devBepr(3,3),devgepr(3,3),alpha1,Bepr(3,3)
      real*8  kappa,p,devT(3,3),T(3,3),TGpa(3,3)
      integer i1,i2

c$$$      PRINT *, 'Dt ', Dt, ' K ', K, ' mu ', mu, ' a0 ', a0, ' b0 ', b0
c$$$      PRINT *, 'a1 ', a1, ' b1 ', b1, ' m ', m, ' ks ', kappas, ' '
c$$$      PRINT *, 'Jn, ', Jn
c$$$      PRINT *, 'Beprn ', Beprn(1,1), ' ', Beprn(1,2), ' ', Beprn(1,3)
c$$$      PRINT *, 'Beprn ', Beprn(2,1), ' ', Beprn(2,2), ' ', Beprn(2,3)
c$$$      PRINT *, 'Beprn ', Beprn(3,1), ' ', Beprn(3,2), ' ', Beprn(3,3)
c$$$      PRINT *, 'kn ', kappan
c$$$      PRINT *, 'F', Fr(1,1), ' ', Fr(1,2), ' ', Fr(1,3)
c$$$      PRINT *, 'F', Fr(2,1), ' ', Fr(2,2), ' ', Fr(2,3)
c$$$      PRINT *, 'F', Fr(3,1), ' ', Fr(3,2), ' ', Fr(3,3)
c$$$      PRINT *, 'J ', J
c$$$      PRINT *, 'Bepr ', Bepr(1,1), ' ', Bepr(1,2), ' ', Bepr(1,3)
c$$$      PRINT *, 'Bepr ', Bepr(2,1), ' ', Bepr(2,2), ' ', Bepr(2,3)
c$$$      PRINT *, 'Bepr ', Bepr(3,1), ' ', Bepr(3,2), ' ', Bepr(3,3)
c$$$      PRINT *, 'k ', kappa
c$$$      PRINT *, 'TGPa ', TGPa(1,1), ' ', TGPa(1,2), ' ', TGPa(1,3)
c$$$      PRINT *, 'TGPa ', TGPa(2,1), ' ', TGPa(2,2), ' ', TGPa(2,3)
c$$$      PRINT *, 'TGPa ', TGPa(3,1), ' ', TGPa(3,2), ' ', TGPa(3,3)

c     Philipc bug-fix: EyeTens needs to be set here
      call IdentityTens(EyeTens)

      call DetTens(Fr,Jr)
      J=Jr*Jn

      call Unimodular(Fr,Frpr)
      call ElasticTrialDist(Frpr,Beprn,Beprs)
      call DevTens(Beprs,devBeprs,hydBeprs)
      do i1=1,3
        do i2=1,3
          devgeprs(i1,i2)=0.50d0*devBeprs(i1,i2)
        end do
      end do

      call EquivStrain(devBeprs,games)

      call TransTens(Fr,FrT)
      call JuxtaTensTens(Fr,FrT,Br)
      do i1=1,3
        do i2=1,3
          Dbar(i1,i2)=0.50d0/Dt*(Br(i1,i2)-EyeTens(i1,i2))
        end do
      end do
      call EffDistStrain(Dbar,Dt,Deps)

      call DeltatGamma(a0,a1,b0,b1,m,kappas,Dt,kappan,games,Deps,gam,
     #                 GAMMA)

      call Beprpr(devgeprs,Dt,GAMMA,devBepr)
      do i1=1,3
        do i2=1,3
          devgepr(i1,i2)=0.50d0*devBepr(i1,i2)
        end do
      end do
      call OneThirdAlpha(devBepr,alpha1)
      call Beprime(devBepr,alpha1,Bepr)

      call CauchyStress(K,mu,J,devgepr,p,devT,T)
      do i1=1,3
        do i2=1,3
          TGpa(i1,i2)=T(i1,i2)*10.0d0**(-9.0d0)
        end do
      end do

      call hardening(kappan,kappas,m,gam,kappa)

      end subroutine

c   ===========================
c   DISTORTIONAL ELASTIC TRIAL
c   ===========================

c   Compute the elastic trial value Be'* of Be'

      subroutine ElasticTrialDist(Frpr,Beprn,Bepret)

      implicit none

c   Declarations
      real*8 Frpr(3,3),FrprT(3,3),Beprn(3,3),aux(3,3),Bepret(3,3)

      call TransTens(Frpr,FrprT)
      call JuxtaTensTens(Frpr,Beprn,aux)
      call JuxtaTensTens(aux,FrprT,Bepret)

      end subroutine


C   ==========================
C   EQUIVALENT STRAIN gamma_e
C   ==========================

      subroutine EquivStrain(devBepr,game)

      implicit none

C   Declarations
      real*8 devBepr(3,3),game

      call DotTensTens(devBepr,devBepr,game)
      game=dsqrt((3.0d0/8.0d0)*game)

      end subroutine


C   ===============================================
C   EFFECTIVE TOTAL DISTORTIONAL DEFORMATION D(eps)
C   ===============================================

      subroutine EffDistStrain(Dbar,Dt,Deps)

      implicit none

C   Declarations
      real*8 Dbar(3,3),Dt,Deps,devDbar(3,3),hydDbar(3,3)

      call DevTens(Dbar,devDbar,hydDbar)
      call DotTensTens(devDbar,devDbar,Deps)
      Deps=Dt*dsqrt(2.0d0/3.0d0*Deps)

      end subroutine


c   ============================
c   Deviatoric part Be'' of Be'
c   ============================

      subroutine Beprpr(devgeprs,Dt,GAMMA,devBepr)

      implicit none

C   Declarations
      real*8 Dt,GAMMA,devgeprs(3,3),devBepr(3,3),fac1
      integer i1,i2

      fac1=2.0d0/(1.0d0+(Dt*GAMMA))

      do i1=1,3
        do i2=1,3
          devBepr(i1,i2)=fac1*devgeprs(i1,i2)
        end do
      end do

      end subroutine


C   =======================================================
C   CALCULATION OF THE 1ST NON-TRIVIAL INVARIANT (alpha_1)
C   =======================================================

C   Calculate alpha=trace(Bepr) by solving a quadratic equation

      subroutine OneThirdAlpha(devBepr,alpha1)

      implicit none

C   Declarations
      real*8 devBepr(3,3),alpha1
      real*8 fac1,fac2,detdevBepr,dotprod

      call DetTens(devBepr,detdevBepr)
      call DotTensTens(devBepr,devBepr,dotprod)

      fac1=2.0d0*dotprod/3.0d0
      if (fac1.eq.0.0d0) then
        alpha1=3.0d0
      else
        fac2=(4.0d0*(1.0d0-detdevBepr))/(fac1**1.5d0)
        if (fac2.ge.1.0d0) then
          alpha1=3.0d0*dsqrt(fac1)*dcosh(dacosh(fac2)/3.0d0)
        else
          alpha1=3.0d0*dsqrt(fac1)*dcos(dacos(fac2)/3.0d0)
        end if
      end if

      end subroutine


C   =====================================================
C   CALCULATION OF THE ELASTIC DISTORTIONAL TENSOR (Be')
C   =====================================================

      subroutine Beprime(devBepr,alpha1,Bepr)

      implicit none

C   Declarations
      real*8 EyeTens(3,3),devBepr(3,3),alpha1,Bepr(3,3)
      integer i1,i2

      call IdentityTens(EyeTens)
      do i1=1,3
        do i2=1,3
          Bepr(i1,i2)=alpha1/3.0d0*EyeTens(i1,i2)+devBepr(i1,i2)
        end do
      end do

      end subroutine


C   ===========================================================
C   CACULAITION OF THE AUXILIARY VARIABLE (gamma) and GAMMA
C   ===========================================================

      subroutine DeltatGamma(a0,a1,b0,b1,m,kappas,Dt,kappan,games,
     #                       Deps,gam,GAMMA)

      implicit none

C   Declarations
      real*8 a0,a1,b0,b1,Dt,games,kappan,gam,GAMMA,DtGam0,DtGam1
      real*8 c0,c1,c2,Deps,m,kappas

      DtGam0=Dt*a0+b0*Deps
      DtGam1=Dt*a1+b1*Deps
      c0=DtGam1*(games-(1.0d0+DtGam0)*kappan)

      if (c0.le.0.0d0) then
        gam=0.0d0
      else
        c1=games+DtGam1*(kappan-m*(games-(1+DtGam0)*kappas))
        if (m.eq.0.0d0) then
          gam=c0/c1
        else
          c2=m*(games+DtGam1*kappas)
          gam=(-c1+dsqrt(c1**2.0d0+4.0d0*c0*c2))/2.0d0/c2
        end if
      end if

      GAMMA=DtGam0/Dt+gam/Dt

      end subroutine


C   ========================================
C   CALCULATION OF THE CAUCHY STRESS TENSOR
C   ========================================

      subroutine CauchyStress(K,mu,J,devgepr,p,devT,T)

      implicit none

C   Declarations
      real*8 K,mu,J,devgepr(3,3),p,devT(3,3),T(3,3)
      real*8 EyeTens(3,3)
      integer i1,i2

      call IdentityTens(EyeTens)

      p=-K*(J-1.0d0)
      do i1=1,3
        do i2=1,3
          devT(i1,i2)=2.0d0*mu*devgepr(i1,i2)/J
          T(i1,i2)=-p*EyeTens(i1,i2)+devT(i1,i2)
        end do
      end do

      end subroutine


C   ======================================
C   CALCULATION OF THE HARDENING VARIABLE
C   ======================================

      subroutine hardening(kappan,kappas,m,gam,kappa)

      implicit none

C   Declarations
      real*8 kappan,kappas,m,gam,kappa

      kappa=(kappan+m*kappas*gam)/(1.0d0+m*gam)

      end subroutine


c   %%%%%%%%%%%%%%%%%%
c   TENSOR OPERATIONS
c   %%%%%%%%%%%%%%%%%%

c   Juxtaposition of a 3x3 tensor and a 3x1 vector
c   Also holds for a Dot Product between a 2nd order tensor and a vector

      subroutine JuxtaTensVec(A,v,Av)

      implicit none

c   Declarations
      real*8 A(3,3),v(3),Av(3)
      integer i,j

      do i=1,3
        Av(i)=0.0d0
        do j=1,3
          Av(i)=Av(i)+A(i,j)*v(j)
        end do
      end do

      end subroutine

c   ====================================================================

c   Juxtaposition between two 2nd order tensors

      subroutine JuxtaTensTens(A,B,AB)

      implicit none

c   Declarations
      real*8 A(3,3),B(3,3),AB(3,3)
      integer i,j,k

      do i=1,3
        do j=1,3
          AB(i,j)=0.0d0
          do k=1,3
            AB(i,j)=AB(i,j)+A(i,k)*B(k,j)
          end do
        end do
      end do

      end subroutine

c   ====================================================================

c   Transpose of a 2nd order tensor

      subroutine TransTens(A,At)

      implicit none

c   Declarations
      real*8 A(3,3),At(3,3)
      integer i,j

      do i=1,3
        do j=1,3
          At(i,j)=A(j,i)
        end do
      end do

      end subroutine

c   ====================================================================

c   Determinant of a 2nd order tensor

      subroutine DetTens(A,detA)

      implicit none

c   ----------------
c   Standard method
c   ----------------

c   Declarations
c      real*8 A(3,3),detA,detA11,detA12,detA13

c      detA11=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
c      detA12=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
c      detA13=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

c      detA=detA11-detA12+detA13

c   ------------------------------------
c   Using Using Cayley-Hamilton theorem
c   ------------------------------------

c   Declarations
      real*8 A(3,3),A2(3,3),A3(3,3),detA,trA,trA2,trA3

      call JuxtaTensTens(A,A,A2)
      call JuxtaTensTens(A,A2,A3)

      call TraceTens(A,trA)
      call TraceTens(A2,trA2)
      call TraceTens(A3,trA3)

      detA=1.0d0/3.0d0*(trA3-trA*trA2+0.50d0*(trA**2.0d0-trA2)*trA)

      end subroutine

c   ====================================================================

c   Trace of a 2nd order tensor

      subroutine TraceTens(A,trA)

      implicit none

c   Declarations
      real*8 A(3,3)
      real*8 trA

      trA=A(1,1)+A(2,2)+A(3,3)

      end subroutine

c   ====================================================================

c   Inverse of a 2nd order tensor

      subroutine InvTens(A,invA)

      implicit none

c   ----------------
c   Standard Method
c   ----------------
c   Declarations
c      real*8 A(3,3),detA,invA(3,3)
c
c   Compute the determinant of the tensor
c      call DetTens(A,detA)
c
c   Compute the components of the inverse tensor (adjoint)
c      invA(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))/detA
c      invA(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))/detA
c      invA(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))/detA
c
c      invA(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))/detA
c      invA(2,2)=(A(1,1)*A(3,3)-A(1,3)*A(3,1))/detA
c      invA(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))/detA
c
c      invA(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/detA
c      invA(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))/detA
c      invA(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/detA

c   ------------------------------
c   Using Cayley-Hamilton theorem
c   ------------------------------
c   Declarations
      real*8 EyeTens(3,3),A(3,3),detA,I11,A2(3,3),trA2,I22,invA(3,3)
      integer i1,i2

      call IdentityTens(EyeTens)
      call DetTens(A,detA)
      I11=A(1,1)+A(2,2)+A(3,3)
      call JuxtaTensTens(A,A,A2)
      trA2=A2(1,1)+A2(2,2)+A2(3,3)
      I22=0.50d0*(I11**2.0d0-trA2)

      do i1=1,3
        do i2=1,3
          invA(i1,i2)=1.0d0/detA*(A2(i1,i2)-I11*A(i1,i2)
     #               +I22*EyeTens(i1,i2))
        end do
      end do

      end subroutine

c   ====================================================================

c   Dot Product between two 2nd order tensors

      subroutine DotTensTens(A,B,AdotB)

      implicit none

c   -------------
c   FIRST METHOD
c   -------------
c   Declarations
c       real*8 A(3,3),B(3,3)
c       real*8 AdotB
c       integer i,j
c
c       AdotB=0.0d0;
c       do i=1,3
c         do j=1,3
c           AdotB=AdotB+A(i,j)*B(i,j)
c         end do
c       end do


c   --------------
c   SECOND METHOD
c   --------------
c   Declarations
      real*8 A(3,3),B(3,3),Bt(3,3),ABt(3,3)
      real*8 trABt,AdotB

      call TransTens(B,Bt)
      call JuxtaTensTens(A,Bt,ABt)
      call TraceTens(ABt,trABt)

      AdotB=trABt

      end subroutine

c   ====================================================================

c   2nd order identity tensor

      subroutine IdentityTens(EyeTens)

      implicit none

c   Declarations
      real*8 EyeTens(3,3)

      EyeTens(1,1)=1.0d0
      EyeTens(1,2)=0.0d0
      EyeTens(1,3)=0.0d0
      EyeTens(2,1)=0.0d0
      EyeTens(2,2)=1.0d0
      EyeTens(2,3)=0.0d0
      EyeTens(3,1)=0.0d0
      EyeTens(3,2)=0.0d0
      EyeTens(3,3)=1.0d0

      end subroutine

c   ====================================================================

c   2nd order sparse tensor

      subroutine ZeroTens(zero3)

      implicit none

c   Declarations
      real*8 zero3(3,3)

      zero3(1,1)=0.0d0
      zero3(1,2)=0.0d0
      zero3(1,3)=0.0d0
      zero3(2,1)=0.0d0
      zero3(2,2)=0.0d0
      zero3(2,3)=0.0d0
      zero3(3,1)=0.0d0
      zero3(3,2)=0.0d0
      zero3(3,3)=0.0d0

      end subroutine

c   ====================================================================

c   Deviator of a 2nd order tensor

      subroutine DevTens(A,devA,hydA)

      implicit none

c   Declarations
      real*8 A(3,3),devA(3,3),EyeTens(3,3),hydA(3,3)
      real*8 p,trA
      integer i,j

      call TraceTens(A,trA)
      p=-trA/3.0d0
      call IdentityTens(EyeTens)

      do i=1,3
        do j=1,3
          hydA(i,j)=-p*EyeTens(i,j)
          devA(i,j)=A(i,j)-hydA(i,j)
        end do
      end do

      end subroutine

c   ====================================================================

c   Unimodular 2nd order tensor

      subroutine Unimodular(A,Apr)

      implicit none

c   Declarations
      real*8 A(3,3),Apr(3,3)
      real*8 detA,power,fac
      integer i,j
      power=-1.0d0/3.0d0

      call DetTens(A,detA)

      do i=1,3
        do j=1,3
          Apr(i,j)=sign(dabs(detA)**power, detA)*A(i,j)
        end do
      end do

c       do i=1,3
c         do j=1,3
c           Apr(i,j)=fac*A(i,j)
c         enddo
c       enddo

      end subroutine

c   ====================================================================

c   Cross Product between two vectors

      subroutine CrossVecVec(a,b,axb)

      implicit none

c   Declarations
      real*8 a(3),b(3),axb(3)

      axb(1)=a(2)*b(3)-a(3)*b(2)
      axb(2)=a(3)*b(1)-a(1)*b(3)
      axb(3)=a(1)*b(2)-a(2)*b(1)

      end subroutine

c   ====================================================================

c   Dot Product between two vectors

      subroutine DotVecVec(a,b,adotb)

      implicit none

      real*8 a(3),b(3)
      real*8 adotb

      adotb=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

      end subroutine

c   ====================================================================

c   Tensor Product between two vectors

      subroutine TensProd(a,b,aoxb)

      implicit none

      real*8 a(3),b(3),aoxb(3,3)
      integer i,j

      do i=1,3
        do j=1,3
          aoxb(i,j)=a(i)*b(j)
        end do
      end do

      end subroutine

c   ====================================================================

c   Tensor Product berween two 2nd order tensors

      subroutine TensProd33(A,B,AoxB)

      implicit none

      real*8 A(3,3),B(3,3),AoxB(3,3,3,3)
      integer i,j,m,n

      do i=1,3
        do j=1,3
          do m=1,3
            do n=1,3
              AoxB(i,j,m,n)=A(i,j)*B(m,n)
            end do
          end do
        end do
      end do

      end subroutine

c   ====================================================================

c   The Operation "oplus" between two 2nd order tensors

      subroutine oplus(A,B,AopB)

      implicit none

      real*8 A(3,3),B(3,3),AopB(3,3,3,3)
      integer i,j,m,n

      do i=1,3
        do j=1,3
          do m=1,3
            do n=1,3
              AopB(i,j,m,n)=A(i,n)*B(j,m)
            end do
          end do
        end do
      end do

      end subroutine

c   ====================================================================

c   The Operation "ominus" between two 2nd order tensors

      subroutine ominus(A,B,AomB)

      implicit none

      real*8 A(3,3),B(3,3),AomB(3,3,3,3)
      integer i,j,m,n

      do i=1,3
        do j=1,3
          do m=1,3
            do n=1,3
              AomB(i,j,m,n)=A(i,m)*B(j,n)
            end do
          end do
        end do
      end do

      end subroutine

c   ====================================================================

c   Dot Product between a 2nd order tensor and a 4th order tensor

      subroutine DotTens2Tens4(A2,B4,A2dotB4)

      implicit none

c   Declarations
      real*8 A2(3,3),B4(3,3,3,3),A2dotB4(3,3)
      integer i,j,m,n

      do i=1,3
        do j=1,3
          A2dotB4(i,j)=0.0d0
          do m=1,3
            do n=1,3
              A2dotB4(i,j)=A2dotB4(i,j)+A2(m,n)*B4(m,n,i,j)
            end do
          end do
        end do
      end do

      end subroutine

c   ====================================================================

c   Dot Product between a 4th order tensor and a 2nd order tensor

      subroutine DotTens4Tens2(A4,B2,A4dotB2)

      implicit none

c   Declarations
      real*8 A4(3,3,3,3),B2(3,3),A4dotB2(3,3)
      integer i,j,m,n

      do i=1,3
        do j=1,3
          A4dotB2(i,j)=0.0d0
          do m=1,3
            do n=1,3
              A4dotB2(i,j)=A4dotB2(i,j)+A4(i,j,m,n)*B2(m,n)
            end do
          end do
        end do
      end do

      end subroutine

c   ====================================================================
