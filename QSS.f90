!Albert Isanta
!Gener 2026
!UB-Física computacional M1

module var
    implicit none

    ! Definim dp amb un numero random en dp per poder declarar complexos de doble precisió (no he trobat altra manera)
    integer,parameter::dp=kind(0.d0)
    double precision,parameter::pi=acos(-1.d0)
    complex(dp),parameter::iunit=cmplx(0.d0,1.d0,dp) ! Unitat imaginària en doble precisió
end module

module subrutinesP11
    use var
    implicit none
    contains
    subroutine init_cn_sch(nx,dx,dt,V,a,b,c)
        !! Inicialitza els coeficients de la matriu de Crank-Nicolson per a l'equació de Schrödinger.
        !! @param nx int,in: nombre de punts de la malla.
        !! @param dx dp,in: espaiat entre punts de la malla.
        !! @param dt dp,in: pas de temps.
        !! @param V(0:nx) dp,in: vector del potencial.
        !! @param a(1:nx-1) complex(dp),out: diagonal inferior de la matriu.
        !! @param b(1:nx-1) complex(dp),out: diagonal principal de la matriu.
        !! @param c(1:nx-1) complex(dp),out: diagonal superior de la matriu.
        use var
        implicit none
        integer,intent(in)::nx
        double precision,intent(in)::dx,dt,V(0:nx)
        complex(dp),intent(out)::a(1:nx-1),b(1:nx-1),c(1:nx-1)
        integer::i
        double precision::alpha

        alpha=dt/(4.d0*dx*dx)

        ! Coeficients de la matriu tridiagonal (I + i*dt/2 * H)
        do i=1,nx-1
            a(i)=-iunit*alpha                                              ! Subdiagonal
            c(i)=-iunit*alpha                                              ! Superdiagonal
            b(i)=cmplx(1.d0,0.d0,dp)+iunit*(dt/2.d0)*((1.d0/(dx*dx))+V(i)) ! Diagonal
        end do

        ! Condicions de contorn (extrems fixos a zero)
        a(1)=cmplx(0.d0,0.d0,dp)
        c(nx-1)=cmplx(0.d0,0.d0,dp)
    end subroutine

    subroutine pas_cn_sch(nx,dx,dt,V,a,b,c,psi,psi_nou,rhs)
        !! Realitza un pas temporal utilitzant el mètode de Crank-Nicolson.
        !! @param nx int,in: nombre de punts de la malla.
        !! @param dx dp,in: espaiat entre punts de la malla.
        !! @param dt dp,in: pas de temps.
        !! @param V(0:nx) dp,in: vector del potencial.
        !! @param a(1:nx-1) complex(dp),in: diagonal inferior.
        !! @param b(1:nx-1) complex(dp),in: diagonal principal.
        !! @param c(1:nx-1) complex(dp),in: diagonal superior.
        !! @param psi(0:nx) complex(dp),in: funció d'ona al pas actual.
        !! @param psi_nou(0:nx) complex(dp),out: funció d'ona al pas següent.
        !! @param rhs(1:nx-1) complex(dp),out: vector del terme de la dreta.
        use var
        implicit none
        integer,intent(in)::nx
        double precision,intent(in)::dx,dt,V(0:nx)
        complex(dp),intent(in)::a(1:nx-1),b(1:nx-1),c(1:nx-1)
        complex(dp),intent(in)::psi(0:nx)
        complex(dp),intent(out)::psi_nou(0:nx),rhs(1:nx-1)
        integer::i
        double precision::alpha
        complex(dp)::coef_centre,coef_veins

        alpha=dt/(4.d0*dx*dx)

        ! Càlcul del terme de la dreta (RHS): (I - i*dt/2 * H) * psi^n
        do i=1,nx-1
            coef_centre=cmplx(1.d0,0.d0,dp)-iunit*(dt/2.d0)*((1.d0/(dx*dx))+V(i))
            coef_veins=iunit*alpha
            rhs(i)=coef_centre*psi(i)+coef_veins*(psi(i-1)+psi(i+1))
        end do

        ! Resolució del sistema d'equacions tridiagonal (Algorisme de Thomas)
        call tridiag(nx-1,a,b,c,rhs,psi_nou(1:nx-1))

        ! Assegurem condicions de contorn a zero
        psi_nou(0)=cmplx(0.d0,0.d0,dp)
        psi_nou(nx)=cmplx(0.d0,0.d0,dp)
    end subroutine

    subroutine calcula_norma(nx,dx,psi,norma,Pesq,Pdre)
        !! Calcula la norma de la funció d'ona i les probabilitats integrant amb trapezis.
        !! @param nx int,in: nombre de punts de la malla.
        !! @param dx dp,in: espaiat entre punts.
        !! @param psi(0:nx) complex(dp),in: funció d'ona.
        !! @param norma dp,out: norma total.
        !! @param Pesq dp,out: probabilitat a l'esquerra (x<0).
        !! @param Pdre dp,out: probabilitat a la dreta (x>0).
        use var
        implicit none
        integer,intent(in)::nx
        double precision,intent(in)::dx
        complex(dp),intent(in)::psi(0:nx)
        double precision,intent(out)::norma,Pesq,Pdre
        integer::i,idx0
        double precision::s

        idx0=nx/2

        ! Norma total: integral de |psi|^2 en tot el domini
        s=0.5d0*abs(psi(0))**2
        do i=1,nx-1
            s=s+abs(psi(i))**2
        end do
        s=s+0.5d0*abs(psi(nx))**2
        norma=dx*s

        ! Probabilitat a l'esquerra (x < 0)
        s=0.5d0*abs(psi(0))**2
        do i=1,idx0-1
            s=s+abs(psi(i))**2
        end do
        s=0.5d0*abs(psi(idx0))**2
        Pesq=dx*s

        ! Probabilitat a la dreta (x > 0)
        s=0.5d0*abs(psi(idx0))**2
        do i=idx0+1,nx-1
            s=s+abs(psi(i))**2
        end do
        s=0.5d0*abs(psi(nx))**2
        Pdre=dx*s
    end subroutine

    subroutine energies(nx,dx,V,psi,Etot,Ecin,Epot)
        !! Calcula les energies cinètica, potencial i total del paquet d'ones.
        !! @param nx int,in: nombre de punts.
        !! @param dx dp,in: espaiat entre punts.
        !! @param V(0:nx) dp,in: potencial.
        !! @param psi(0:nx) complex(dp),in: funció d'ona.
        !! @param Etot dp,out: energia total.
        !! @param Ecin dp,out: energia cinètica.
        !! @param Epot dp,out: energia potencial.
        use var
        implicit none
        integer,intent(in)::nx
        double precision,intent(in)::dx,V(0:nx)
        complex(dp),intent(in)::psi(0:nx)
        double precision,intent(out)::Etot,Ecin,Epot
        integer::i
        complex(dp)::sumk
        double precision::sump

        ! Energia cinètica: valor d'esperança de l'operador laplacià
        sumk=cmplx(0.d0,0.d0,dp)
        do i=1,nx-1
            sumk=sumk+conjg(psi(i))*(-0.5d0*(psi(i+1)-2.d0*psi(i)+psi(i-1))/dx**2)
        end do
        Ecin=dble(sumk)*dx

        ! Energia potencial: integral de V(x) * |psi(x)|^2
        sump=0.d0
        do i=1,nx-1
            sump=sump+V(i)*abs(psi(i))**2
        end do
        Epot=sump*dx

        Etot=Ecin+Epot
    end subroutine

    subroutine pas_adi_2d(nx,ny,dx,dy,dt,V,psi,psi_nou)
        !! Implementa un pas temporal de l'equació de Schrödinger 2D usant el mètode ADI.
        !! @param nx int,in: nombre de punts en l'eix x.
        !! @param ny int,in: nombre de punts en l'eix y.
        !! @param dx dp,in: pas de malla en x.
        !! @param dy dp,in: pas de malla en y.
        !! @param dt dp,in: pas de temps total per a un pas complet.
        !! @param V(0:nx,0:ny) dp,in: potencial en 2D.
        !! @param psi(0:nx,0:ny) complex(dp),in: funció d'ona al pas n.
        !! @param psi_nou(0:nx,0:ny) complex(dp),out: funció d'ona al pas n+1.
        use var
        implicit none
        integer,intent(in)::nx,ny
        double precision,intent(in)::dx,dy,dt,V(0:nx,0:ny)
        complex(dp),intent(in)::psi(0:nx,0:ny)
        complex(dp),intent(out)::psi_nou(0:nx,0:ny)
        complex(dp)::psi_inter(0:nx,0:ny)
        complex(dp)::ax(1:nx-1),bx(1:nx-1),cx(1:nx-1),rhsx(1:nx-1)
        complex(dp)::ay(1:ny-1),by(1:ny-1),cy(1:ny-1),rhsy(1:ny-1)
        complex(dp)::rx,ry
        integer::i,j

        ! Coeficients de dispersió (el factor 4 ve de dt/2 i la discretització de la segona derivada)
        rx=iunit*dt/(4.d0*dx*dx)
        ry=iunit*dt/(4.d0*dy*dy)

        ! 1) Primer sub-pas: Implícit en X, Explícit en Y
        do j=1,ny-1
            do i=1,nx-1
                ax(i)=-rx
                cx(i)=-rx
                bx(i)=cmplx(1.d0,0.d0,dp)+2.d0*rx+iunit*(dt/4.d0)*V(i,j)
                rhsx(i)=(1.d0-2.d0*ry-iunit*(dt/4.d0)*V(i,j))*psi(i,j)+ry*(psi(i,j-1)+psi(i,j+1))
            end do
            
            call tridiag(nx-1,ax,bx,cx,rhsx,psi_inter(1:nx-1,j))
        end do

        ! Condicions de contorn per a l'estat intermedi (parets infinites)
        psi_inter(0,:)=cmplx(0.d0,0.d0,dp)
        psi_inter(nx,:)=cmplx(0.d0,0.d0,dp)
        psi_inter(:,0)=cmplx(0.d0,0.d0,dp)
        psi_inter(:,ny)=cmplx(0.d0,0.d0,dp)

        ! 2) Segon sub-pas: Implícit en Y, Explícit en X
        do i=1,nx-1
            do j=1,ny-1
                ay(j)=-ry
                cy(j)=-ry
                by(j)=cmplx(1.d0,0.d0,dp)+2.d0*ry+iunit*(dt/4.d0)*V(i,j)
                rhsy(j)=(1.d0-2.d0*rx-iunit*(dt/4.d0)*V(i,j))*psi_inter(i,j)+rx*(psi_inter(i-1,j)+psi_inter(i+1,j))
            end do
            
            call tridiag(ny-1,ay,by,cy,rhsy,psi_nou(i,1:ny-1))
        end do

        ! Condicions de contorn finals
        psi_nou(0,:)=cmplx(0.d0,0.d0,dp)
        psi_nou(nx,:)=cmplx(0.d0,0.d0,dp)
        psi_nou(:,0)=cmplx(0.d0,0.d0,dp)
        psi_nou(:,ny)=cmplx(0.d0,0.d0,dp)
    end subroutine

    subroutine tridiag(imax,a,b,c,r,psi)
        !! Resolució d'un sistema d'equacions lineals tridiagonal (Algorisme de Thomas) per a nombres complexos.
        !! @param imax int,in: dimensió de la matriu.
        !! @param a(imax) complex(dp),in: diagonal inferior.
        !! @param b(imax) complex(dp),in: diagonal principal.
        !! @param c(imax) complex(dp),in: diagonal superior.
        !! @param r(imax) complex(dp),in: vector del terme independent.
        !! @param psi(imax) complex(dp),out: vector solució.
        use var
        implicit none
        integer,intent(in)::imax
        complex(dp),dimension(imax),intent(in)::a,b,c,r
        complex(dp),dimension(imax),intent(out)::psi
        complex(dp)::bet
        complex(dp),dimension(imax)::gam
        integer::j

        if(abs(b(1))==0.d0) then
            stop "error: b(1) es zero a tridiag"
        end if

        bet=b(1)
        psi(1)=r(1)/bet

        do j=2,imax
            gam(j)=c(j-1)/bet
            bet=b(j)-a(j)*gam(j)
            if(abs(bet)==0.d0) then
                stop "error: divisio per zero (bet) a tridiag"
            end if
            psi(j)=(r(j)-a(j)*psi(j-1))/bet
        end do

        do j=imax-1,1,-1
            psi(j)=psi(j)-gam(j+1)*psi(j+1)
        end do
    end subroutine
end module

program P11
    use var
    use subrutinesP11
    implicit none

    integer,parameter::nx=500
    double precision,parameter::Lx=17.d0
    double precision,parameter::tmax = 20.d0
    double precision,parameter::dt=0.002d0

    integer::i,j,it,nt,nout,idx0
    double precision::dx,temps
    double precision::x(0:nx),V(0:nx)

    complex(dp)::psi_ini(0:nx),psi(0:nx),psi_nou(0:nx)
    complex(dp)::a(1:nx-1),b(1:nx-1),c(1:nx-1),rhs(1:nx-1)
    complex(dp)::prod

    double precision::norma,Pesq,Pdre
    double precision::Etot,Ecin,Epot,Etot_ini
    double precision::x0,p0,sigma,amp
    double precision::difmax,solapament

    dx=2.d0*Lx/dble(nx)
    nt=nint(tmax/dt)
    idx0=nx/2
    nout=20

    ! Mallat espacial: x in [-Lx,Lx]
    do i=0,nx
        x(i)=-Lx+dble(i)*dx
    end do

    ! 1) Potencial barrera (gaussiana centrada a 0) --------------------------------------------

    do i=0,nx
        V(i)=42.55d0*exp(-(x(i)**2)/0.1d0)/sqrt(0.1d0*pi)
    end do

    ! 2) Condició inicial: paquet d'ones i normalització ---------------------------------------

    x0=7.d0
    p0=200.d0/Lx
    sigma=1.d0

    ! Factor de normalització analític per a la gaussiana
    amp=(2.d0*pi*sigma**2)**(-0.25d0)

    do i=0,nx
        psi_ini(i)=amp*exp(-((x(i)-x0)**2)/(4.d0*sigma**2))*exp(-iunit*p0*x(i))
    end do

    ! Condicions de contorn: parets infinites (psi=0 als extrems)
    psi_ini(0)=cmplx(0.d0,0.d0,dp)
    psi_ini(nx)=cmplx(0.d0,0.d0,dp)

    ! Normalització numèrica inicial per assegurar norma = 1
    call calcula_norma(nx,dx,psi_ini,norma,Pesq,Pdre)
    psi_ini=psi_ini/sqrt(norma)

    psi=psi_ini

    ! Càlcul de l'energia inicial (cinètica + potencial)
    call energies(nx,dx,V,psi,Etot_ini,Ecin,Epot)

    write(*,*) "dx=",dx," dt=",dt," nt=",nt
    write(*,*) "E inicial=",Etot_ini,"  (Ecin=",Ecin," Epot=",Epot,")"

    ! 3) Matriu Crank-Nicolson (constant per V fix) + fitxers ----------------------------------
    ! Es precalculen les diagonals de la matriu (I + i*dt/2 * H)

    call init_cn_sch(nx,dx,dt,V,a,b,c)

    open(10,file="figP11.dat")
    open(11,file="figP12.dat")

    write(11,*) "# t  Etot  Ecin  Epot  norma  Pesq(x<0)  Pdre(x>0)"

    temps=0.d0
    call calcula_norma(nx,dx,psi,norma,Pesq,Pdre)
    call energies(nx,dx,V,psi,Etot,Ecin,Epot)
    write(11,*) temps,Etot,Ecin,Epot,norma,Pesq,Pdre

    ! Guardem l'estat inicial per a l'animació
    write(10,*) "# it=0 t=0"
    do i=0,nx
        write(10,*) x(i),dble(psi(i)),aimag(psi(i)),abs(psi(i))**2,V(i)
    end do

    write(10,*)
    write(10,*)

    ! 4) Evolució temporal Crank-Nicolson + energia + reflexió/transmissió --------------------
    ! Es resol el sistema tridiagonal a cada pas de temps

    do it=1,nt
        ! Pas de Crank-Nicolson: (I + i*dt/2 * H) psi^{n+1} = (I - i*dt/2 * H) psi^n
        call pas_cn_sch(nx,dx,dt,V,a,b,c,psi,psi_nou,rhs)

        psi=psi_nou
        temps=dble(it)*dt

        ! Monitoritzem la norma i les energies durant l'evolució
        call calcula_norma(nx,dx,psi,norma,Pesq,Pdre)
        call energies(nx,dx,V,psi,Etot,Ecin,Epot)
        write(11,*) temps,Etot,Ecin,Epot,norma,Pesq,Pdre

        ! Guardem dades per a l'animació cada nout passos
        if(mod(it,nout)==0) then
            write(10,*) "# it=",it," t=",temps
            do i=0,nx
                ! dble(psi(i)) assegura que no perdem precisió en escriure al fitxer
                write(10,*) x(i),dble(psi(i)),aimag(psi(i)),abs(psi(i))**2,V(i)
            end do
            write(10,*)
            write(10,*)
        end if
    end do

    close(10)
    close(11)

    write(*,*) "Final: norma=",norma,"  Pesq=",Pesq,"  Pdre=",Pdre
    write(*,*) "E final=",Etot,"  deltaE=",Etot-Etot_ini

    ! 5) Comprovació de reversibilitat: anar fins tmax i tornar invertint dt ------------------
    ! Si el mètode és unitari, hauríem de recuperar l'estat inicial
    psi_nou=psi

    ! Re-inicialitzem la matriu amb pas de temps negatiu
    call init_cn_sch(nx,dx,-dt,V,a,b,c)

    psi=psi_nou
    do it=1,nt
        call pas_cn_sch(nx,dx,-dt,V,a,b,c,psi,psi_nou,rhs)
        psi=psi_nou
    end do

    ! Calculem la diferència màxima i el solapament amb l'estat inicial
    difmax=0.d0
    do i=0,nx
        difmax=max(difmax,abs(psi(i)-psi_ini(i)))
    end do

    prod=cmplx(0.d0,0.d0,dp)
    do i=1,nx-1
        prod=prod+conjg(psi_ini(i))*psi(i)
    end do
    prod=prod*dx
    solapament=abs(prod)

    ! Guardem la diferència final per a la figura 13
    open(12,file="figP13.dat")
    write(12,*) "# x  |psi_back-psi_ini|^2"
    do i=0,nx
        write(12,*) x(i),abs(psi(i)-psi_ini(i))**2
    end do
    close(12)

    write(*,*) "Tornar enrere: max|dif|=",difmax,"  solapament=",solapament
    write(*,*) "Fet: figP11.dat, figP12.dat, figP13.dat"

    ! 6) EXTRA: Simulació en 2D usant el mètode ADI --------------------------------------------
    write(*,*) "Iniciant simulacio 2D (extra)..."
    
    block
        integer,parameter::nx2=150,ny2=150
        integer::nt2,nout2
        double precision,parameter::Lx2=10.d0,Ly2=10.d0
        double precision,parameter::tmax2=5.d0,dt2=0.01d0
        double precision::dx2,dy2
        double precision,allocatable::x2(:),y2(:),V2(:,:)
        complex(dp),allocatable::p2(:,:),p2_n(:,:)
        double precision::x02,y02,p0x,p0y,s2,norm2,a2
        
        allocate(x2(0:nx2),y2(0:ny2),V2(0:nx2,0:ny2))
        allocate(p2(0:nx2,0:ny2),p2_n(0:nx2,0:ny2))
        
        dx2=2.d0*Lx2/dble(nx2)
        dy2=2.d0*Ly2/dble(ny2)
        nt2=nint(tmax2/dt2)
        nout2=10
        
        do i=0,nx2
            x2(i)=-Lx2+dble(i)*dx2
        end do
        do j=0,ny2
            y2(j)=-Ly2+dble(j)*dy2
        end do
        
        ! Potencial: Una barrera circular al centre
        do i=0,nx2
            do j=0,ny2
                if(sqrt(x2(i)**2+y2(j)**2)<1.5d0) then
                    V2(i,j)=50.d0
                else
                    V2(i,j)=0.d0
                end if
            end do
        end do
        
        ! Condició inicial: Paquet gaussià 2D
        x02=-5.d0; y02=0.d0
        p0x=10.d0; p0y=0.d0
        s2=0.8d0
        a2=(2.d0*pi*s2**2)**(-0.5d0)
        
        do i=0,nx2
            do j=0,ny2
                p2(i,j)=a2*exp(-((x2(i)-x02)**2+(y2(j)-y02)**2)/(4.d0*s2**2))*exp(iunit*(p0x*x2(i)+p0y*y2(j)))
            end do
        end do
        
        ! Normalització numèrica 2D (trapezis simplificat)
        norm2=0.d0
        do i=1,nx2-1
            do j=1,ny2-1
                norm2=norm2+abs(p2(i,j))**2
            end do
        end do
        p2=p2/sqrt(norm2*dx2*dy2)
        
        open(20,file="figP11_2D.dat")
        write(*,*) "Guardant dades 2D a figP11_2D.dat..."
        
        do it=0,nt2
            if(mod(it,nout2)==0) then
                write(20,*) "# t=",dble(it)*dt2
                do i=0,nx2
                    do j=0,ny2
                        write(20,*) x2(i),y2(j),abs(p2(i,j))**2
                    end do
                    write(20,*) ""
                end do
                write(20,*) ""
            end if
            
            call pas_adi_2d(nx2,ny2,dx2,dy2,dt2,V2,p2,p2_n)
            p2=p2_n
        end do
        close(20)
        
        deallocate(x2,y2,V2,p2,p2_n)
        write(*,*) "Simulacio 2D finalitzada."
    end block

end program

