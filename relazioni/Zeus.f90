MODULE DATA
    integer :: N, i
    real*8 :: pi,cmpc,cmkpc,yr,kbol,mu,mp,msol, cv, gam
    parameter(N=2000)
    parameter(pi=3.141592)
    parameter(cmpc=3.085d18)
    parameter(cmkpc=1000.*cmpc)
    parameter(yr=3.156d7)
    parameter(kbol=1.38d-16)
    parameter(mu=1.4)
    parameter(mp=1.67d-24)
    parameter(msol= 1.989d33)
    parameter(gam = 1.667)
    parameter(cv = 2.d8)
END MODULE DATA

program zeus
    use data
    implicit real*8 (a-h,o-z) !  i l m n

    real*8 :: lx1, lum_x, kinetic, m_lost, lum_bol
    real*8 :: v(N), e(N), p(N), d(N), t(N), xa(N), xb(N), dxa(N), dxb(N)
    real*8 :: g2a(N), g2b(N), g31a(N), g31b(N), dvl1a(N), dvl1b(N), deltav(N)
    real*8 :: divV(N), dstar(N), s(N), F1(N), M(N), e_dstar(N), F2(N), vstar(N), F3(N)
    real*8 :: q(N)
    real*8 :: v1(N),v2(N),v3(N),v4(N),v5(N),p1(N),p2(N),p3(N),p4(N),p5(N)
    real*8 :: d1(N),d2(N),d3(N),d4(N),d5(N)
    real*8 :: t1(N), t2(N), t3(N), t4(N), t5(N), dn(N)
    !real*8 , dimension(:) :: sim_temp
    integer :: sdr, Num, ncicli
    real*8, EXTERNAL :: Cool, wsum
    real, dimension(5) :: tmax
    logical :: no_cooling, orangotango, winds
    character*12 :: nome_file, nome_file2



!============== GRID CREATION ===============

    !questo è un test per vedere se github funziona
    
    
    xmin=0.
    xmax=200.*cmpc
    xa(1)= xmin+(xmax-xmin)*(-1.)/(N-1.)
    xa(2)=0.
    !GRIGLIA "a"
    do i=2,N
        xa(i)= xmin+(xmax-xmin)*(i-2.)/(N-1.)
    end do
    xa(1)=-xa(3)
    deltax=xa(3)-xa(2)

    !GRIGLIA "b"
    do i=2, N-1
        xb(i)=0.5*(xa(i)+xa(i+1))
    end do
    xb(N)=xb(N-1)+(xb(N-1)-xb(N-2))   !! add the last calculated Delta_xb to xb(N-1)

    do i=2, N-1
        dxa(i)=xa(i+1)-xa(i) !delta spaziale
        dxb(i)=xb(i)-xb(i-1)
    end do
    dxa(1)=dxa(2)
    dxa(N)=dxa(N-1)
    dxb(1)=dxb(2)
    dxb(N)=dxb(N-1)

!
!=============== DEFINIZIONE FATTORI DI SCALA METRICI ========================

    sdr=1    !! this parameter selects the type of coordinates: 0 = Cartesian

    if (sdr==0) then  !! Cartesian !!

        do i=1, N
        g2a(i)=1.
        g2b(i)=1.
        g31a(i)=1.
        g31b(i)=1.
        end do

        do i=1, N-1
        dvl1a(i)=xa(i+1)-xa(i)   !! Note that is centered in xb(i)
        end do
        dvl1a(N)=dvl1a(N-1)

        do i=2, N
        dvl1b(i)=xb(i)-xb(i-1)  !! Note that it is centered in xa(i)
        end do
        dvl1b(1)=dvl1b(2)



    else if (sdr==1) then   !! spherical !!
        do i=1, N
        g2a(i)=xa(i)
        g31a(i)=xa(i)
        g2b(i)=xb(i)
        g31b(i)=xb(i)
        end do

        do i=1, N-1
        dvl1a(i)=(xa(i+1)**3-xa(i)**3)/3.
        end do
        dvl1a(N)=dvl1a(N-1)
        do i=2, N
        dvl1b(i)=(xb(i)**3-xb(i-1)**3)/3.
        end do
        dvl1b(1)=dvl1b(2)

    end if
!
!=============== SOME OTHER PARAMETER THAT CAN BE CHANGED WHEN NECESSARY =========
    write(*,*)'DO YOU WANT TO SKIP THE COOLING PHASE?'
    write(*,*)'YES == 0'
    WRITE(*,*)'NO == 1'
    read(*,*) scimpanzee

    
    write(*,*)'DO YOU WANT TO SKIP THE STELLAR WINDS PHASE?'
    write(*,*)'YES == 0'
    WRITE(*,*)'NO == 1'
    read(*,*) gorilla

    write(*,*)'PLEASE SELECT THE ISM CONDITIONS'
    write(*,*)'standard ism == 0'
    write(*,*)'hot ionised medium == 1'
    write(*,*)'cold ism == 2'
    read(*,*) ism

    if(scimpanzee == 1 ) then
        no_cooling = .false.  !.true. = skip the temperature evolution
    else 
        no_cooling = .true.
    end if

    if(gorilla == 1 ) then
        winds = .true.          !true = DO NOT skip winds phase
    else 
        winds = .false.
    end if

    if(ism == 0)then
        d0 = 2.d-24      !density of ISM    
        t0 = 1.d4        !temperature of ISM
        t0_min = t0
    else if(ism == 1)then
        !value for HIM test
        d0 = 2.d-26
        t0 = 1.d6
        t0_min = 1.d4
    else
        !value for CIM test
        d0 = 2.d-22
        t0 = 100.
        t0_min = t0
    end if

    !_____________________________________

    e0 = 1.d51       !energy of sna
    C2 = 3           !for artificial viscosity

    cfl = 0.01   !cfl coefficent
    
    orangotango = .false.         ! initial condition for the time writing file
    nome_file = 'sedov000.dat' !output file for enegies and Sedov SOlution and lum
    nome_file2 = 'energ000.dat' !output file for energies



    n0 = d0/(mu*mp)  !numerical density
    vol = 1.333*pi*xa(4)**3
    

    m_lost = 1.d19
    v_winds = 1.d8

!=============== INITIAL CONDITIONS OF THE GRID ==================

    
    do i=1, N
        v(i) = 0
        d(i) = d0                !initial conditions
        t(i) = t0
        e(i) = cv*d(i)*t(i)
        p(i) = (gam-1.)*e(i)
    enddo

    !appropriate condition for t=t0 in a small volume

    1000 format(8(1pe12.4))

    22 continue 

    call file('initial_pressure.dat',xa/cmpc ,p1,p2,p3,p4,p5)
    call file('initial_velocity.dat',xa/cmpc, v1/1d5,v2/1d5,v3/1d5,v4/1d5,v5/1d5)
    call file('initial_density.dat',xa/cmpc, d1,d2,d3,d4,d5)
    call file('initial_temp.dat',xa/cmpc, t1,t2,t3,t4,t5)

    if(winds .eqv. .false.) then
        
        do i=1, 3
            e(i) = e0/vol
            p(i) = e(i)*(gam -1)
            t(i) = e(i)/(cv*d(i))
        enddo

        tmax = (/2.d4*yr, 4.d4*yr, 6.d4*yr, 8.d4*yr, 1.d5*yr/)
    end if


    !STELLAR WINDS PHASE INITIAL CONDITIONS

    if(winds .eqv. .true.)then
      
      print*, 'STELLAR WINDS PHASE'

      

      !call winds_injection(v,d,e,p,t,xb,yr/100) !injection for 1 month 
      !non modifico le condizioni iniziali eccetto la velocità
      do i=1,3
        v(i) = v_winds
        t(i) = t0
        d(i) = d0
        e(i) = cv*d(i)*t(i)
      enddo

      tmax = (/1.d4*yr, 5.d4*yr, 1.d5*yr, 5.d5*yr, 1.d6*yr/)
        
    end if

    !print*, v(3)/1.d5, d(3), e(3), p(3)


!============== METALS INITIAL DISTRIBUTION =================


!============= TIME EVOLUTION ================
    
        ncicli = 0
        time = 0
        
               
    open(28, file=nome_file)
    open(77, file=nome_file2)
   

        

 do j=1,5


        
       do while (time<tmax(j))      !!!! HERE STARTS THE TIME INTEGRATION !!!!!


            ncicli=ncicli+1
            
            bonobo = bonobo +1
            !this condition create a counter in order to non print at each time cycle
            !writing on terminal = bottleneck, make the code faster
                if(bonobo.ge.500)then   
                    orangotango = .true.
                    bonobo = 0
                else
                    orangotango = .false.
                end if
            
                


    !CALC DTMIN

        if(winds .eqv. .true.) then
            do i=1,3
                v(i) = v_winds
                t(i) = t0
                d(i) = d(i) + (m_lost*dtmin/vol) !dens = dens precedente + materia aggiunta
                e(i) =  e(i) + 0.5*m_lost*dtmin*v_winds**2/vol !energia = energia meccanica dei venti
            enddo    
              
            
        end if

        do i=1, N        
            P(i)=(gam-1.)*e(i) !perfect gas equation
        end do



        dtmin=1.d30   !! any very large value !!

        do i=2, N-1
            dtmin=min(dtmin,(xb(i)-xb(i-1))/(abs(v(i))+sqrt(gam*p(i)/d(i))))
        end do

        dtmin=cfl*dtmin !cfl is the C constant 

        time=time+dtmin

     
                if(ncicli == 2)then
                    print*, 'FIRST TIME-STEP == ', dtmin/yr
                    print*, 'initial energy ==', e(2)*vol
                    if(winds .eqv. .false.)then
                        print*, 'initial velocity ==', v(4)/1.d5
                        !print*,'expected temperature'
                        !write(*,1034)  3*mu*mp/(16*kbol)*v(4)**2

                        1034 format(s1pe12.4)
                    else
                        print*, 'initial velocity ==', v(3)/1.d5
                    end if
                end if
       
        



    
    !============= SOURCE STEP ====================

        !SUBSTEP I: AGGIORNAMENTO DELLA VELOCITÃ PER GRADIENTE DI P
        DO i=2, N-1
            v(i)=v(i)-dtmin*2.*(P(i)-P(i-1))/((d(i)+d(i-1))*dxb(i))
            IF (v(i)<1.d-20) THEN
                v(i)=0.
            END IF
        ENDDO
        CALL BCa(v) 


        !CALC Q: ARTIFICIAL VISCOSITY (dv/dx<0)
            do i=2, N-1
                if ((v(i+1)-v(i))<0.) then
                    q(i)=C2*d(i)*(v(i+1)-v(i))**2
                else
                    q(i)=0.
                end if
            end do
            CALL BCb(q)


        ! SUBSTEP II: AGGIORNAMENTO DELLA VELOCITÃ FOR ARTIFICIAL VISC. (as additional term to pression) AND ENERGY
        DO i=2, N-1
            v(i)=v(i)-dtmin*2.*(q(i)-q(i-1))/((d(i)+d(i-1))*dxb(i))
            IF (v(i)<1.d-20) THEN
                v(i)=0.
            END IF
        ENDDO
        CALL BCa(v)

        ! SUBSTEP II: ENERGY 
        DO i=2, N-1
            e(i)=e(i)-dtmin*q(i)*(v(i+1)-v(i))/dxa(i)
        ENDDO
        CALL BCb(e)


        !SUBSTEP III: COMPRESSION HEATING
        DO i=2,N-1
            divV(i)=(g2a(i+1)*g31a(i+1)*v(i+1)-g2a(i)*g31a(i)*v(i))/dvl1a(i)
        ENDDO
        CALL BCa(divV)

        DO i=2, N-1
            e(i)=e(i)*(1.-0.5*dtmin*(gam-1.)*divV(i))/(1.+0.5*dtmin*(gam-1.)*divV(i))
        ENDDO
        CALL BCb(e)


    !========================== TEMPERATURE ===================================
      
    
        if(no_cooling .eqv. .TRUE.)then !to skip the cooling 
            go to 44
        end if

                

                !add the cooling function
                do i=2, N-1
                    e(i)=e(i)-dtmin*(d(i)/2.17d-24)**2 *Cool(t(i)) !energy with en. loss
                  
                end do
                CALL BCb(e)

        44 continue
    
        do i=1, N
            t(i)=e(i)/(cv*d(i)) !temp. with en.loss
        end do
    
    
        do i=1, N
        if (t(i) .le. t0_min) then !condition t>t0_min
            t(i)=t0_min
        end if
        end do
    
    
        do i=2, N-1
            e(i)=cv*d(i)*t(i) !new energy associated to new temperature condition
        end do
        CALL BCb(e)
     
   

    

    !================ SECOND PART ======================




    !TRANSPORT STEP (Upwind FIRST ORDER)

        do i=2, N-1       !! here define the momentum densityP
            s(i)=0.5*(d(i)+d(i-1))*v(i)  !! this is at "i" !!
        end do

        CALL BCa(s)


    !AGGIORNAMENTO DENSITÃ€: SELEZIONO IL VALORE DELLA DENSITÃ€ NEI BORDI DELLE CELLETTE IN BASE ALLA VELOCITÃ€ (UPWIND)
        DO i=2, N-1       !! here select the value of the density at the interface "i"
            IF (v(i)>0.) THEN
                dstar(i)=d(i-1)     !! at i !!
            ELSE
                dstar(i)=d(i)
            END IF
        ENDDO
        dstar(N)=dstar(N-1)
        dstar(1)=dstar(3)

        DO i=2, N
            F1(i)=dstar(i)*v(i)*g2a(i)*g31a(i)    !DENSITY FLUX ATTRAVERSO SURFACES    !! at i !!
        ENDDO


    !AGGIORNAMENTO ENERGIA: CALCOLO MOMENTO SUPERFICIALE ED ENERGIA SUI BORDI DELLE CELLE
            DO i=2, N-1
                M(i)=dstar(i)*v(i) !momento sulla superficie delle cellette
            ENDDO
            CALL BCa(M)


            DO i=2, N-1
                IF (v(i)>0.) THEN
                    e_dstar(i)=e(i-1)/d(i-1)   !! at i !!
                ELSE
                    e_dstar(i)=e(i)/d(i)
                END IF
            ENDDO
            e_dstar(N)=e_dstar(N-1)
            e_dstar(1)=e_dstar(3)
 

    !ORA AGGIORNO LA DENSITÃ€ E L'ENERGIA
        do i=2, N-1
            d(i)=d(i)-dtmin*(F1(i+1)-F1(i))/dvl1a(i) !DENSITÃ€ DOPO L'ARRIVO DEL FLUSSO DI MATERIA
        end do
        CALL BCb(d)

        do i=2, N
            F2(i)=e_dstar(i)*M(i)*g2a(i)*g31a(i) !FLUSSO DI ENERGIA
        end do
        CALL BCa(F2)

        do i=2, N-1
            e(i)=e(i)-dtmin*(F2(i+1)-F2(i))/dvl1a(i) !ENERGIA DOPO L'ARRIVO DEL FLUSSO DI ENERGIA
        end do
        CALL BCb(e)


    !AGGIORNAMENTO MOMENTO
        do i=2, N-1
            if ((v(i-1)+v(i))*0.5>0) then !CALCOLO LE VELOCITa IN TUTTE LE CELLE 
                vstar(i)=v(i-1)       !! at i-1/2  !!
            else
                vstar(i)=v(i)
            end if
        end do
        CALL BCb (vstar)

        do i=1, N-1
            F3(i)=vstar(i+1)*0.5*(M(i)+M(i+1))*g2b(i)*g31b(i)   !! questo e' a i+1/2, occhio !!
        end do

        do i=2, N-1
            s(i)=s(i)-dtmin/dvl1b(i)*(F3(i)-F3(i-1)) !AGGIORNO LA DENSITÃ  DI MOMENTO DAL FLUSSO CALCOLATO SOPRA
        end do
        CALL BCa(s)

        do i=2, N-1
            v(i)=2.*s(i)/(d(i)+d(i-1)) !VELOCITÃ  FINALE DELLE CELLETTE
        end do
        call BCa(v)

        if(cfl<0.5) then
            cfl=cfl+(cfl*0.1)
        else
            cfl=0.5
        end if

     
    
    !======== FINAL PART ==========

        R_shock=(2.*e0/d0)**(1./5.)*(time**(2./5.)) ! sedov solution

        !temperatura teorica post-shock
        !theo_temp = 3.3d6*(e0/1.d51)**(2/5) *(n0)**(-2/5) *(time/1d4/yr)**(-6/5)
        theo_temp = 3*mu*mp*(maxval(v, dim=1))**2 /(16*kbol) 

        !da confrontare con simulatione
        sim_temp = wsum(t,d) !media pesata sulla densità

        !evoluzione teorica wind bubble
        !kk = 0.5*m_lost*dtmin*v_winds**2/vol/1.d36   
        !R_winds = 27*((time/1.d6)**3*kk/n0)**(1/5) !this is already in parsec
        
    !====================== X-RAY LUMINOSITY ===========================
        
         
              

        !only gas with t>10e6 emitt  X rays 
   

        do i=1, N
            dn(i)=d(i)/(mu*mp)            
        enddo

        lum_x = 0   !reset lum x for each time step
        energy_x=0
        energy_bol = 0
        lum_bol = 0
        len = 0 
        do i=2, N-1
            if(t(i)>= 1.d6)then                
                !lum_x = lum_x + 4*pi*(xa(i)**2)*(dn(i)**2) *cool(t(i))*(xa(i)-xa(i-1))
                lum_x = lum_x + 1.333*pi*(xa(i+1)**3-xa(i)**3)*cool(t(i))
                len = len+1
            else
               lum_x = lum_x               
            end if       
            energy_x = energy_x + lum_x*dtmin 
        
        
            if(t(i) >= 1.d4)then
                !lum_bol = lum_bol + 4*pi*(xa(i)**2)*(dn(i)**2) *cool(t(i))*(xa(i)-xa(i-1))
                lum_bol = lum_bol + 1.333*pi*(xa(i+1)**3-xa(i)**3)*cool(t(i))
            else
                lum_bol = lum_bol
            end if                    
                energy_bol = energy_bol + lum_bol*dtmin
        end do
            !la bolometrica ha senso solo se no cooling

        
        !====================== ENERGY FRACTIONS ==================
        
                !similar to x ray lum
            

        thermal = 0
        kinetic = 0
        total = 0

        do i=2, N
            thermal = thermal + 1.333*pi*(xa(i+1)**3-xa(i)**3)*e(i)  
            
            kinetic = kinetic + (0.5)*d(i)*(4./3.)*pi*(xa(i+1)**3-xa(i)**3)*(v(i)**2)

            total = kinetic + thermal 

        enddo

        

        efficency = total/e0    !efficency

        if(winds .eqv. .true.)then
            shock = xa(maxloc(q, dim=1))
        else if(winds .eqv. .false.) then
            shock = xa(maxloc(v, dim=1))
        end if
              

        !temperatura simulata == media pesata sulle densità
        !EVOLUZIONE HOT BUBBLE
        !temp teorica è temp dietro lo shock
        
        !sim_temp = sum(t, mask=t>1.d6) !seleziono hot bubble
        
        !sim_temp = sum(t(maxloc(v, dim=1)-1 : maxloc(v, dim=1)))/ &
        !(size(t(maxloc(v,dim=1)-1: maxloc(v,dim=1))))
        !tre punti dietro fronte di shock

        
        !sim_temp = t(maxloc(v, dim=1))
        
        !sim_temp = sum(sim_temp)/size(sim_temp)
        
        !find a way to plot  the R_shock over time
        !approximation: r_shock = xa corresponding to max vaule of v or p    
        !shock come maggiore differenza di velocità(?)



        if(orangotango .eqv. .true. )then

            write(28,1000)time/yr, shock/cmpc, R_shock/cmpc, lum_x, theo_temp, sim_temp

            write(77,1000)time/yr,log10(thermal/e0) , log10(kinetic/e0), log10(total/e0),log10(energy_bol/e0)

        end if
        !===================================================================================
        
     
       
    enddo !end time cycle
    
        print*, 'TIME PASSED == ', int(tmax(j)/yr), 'YEARS' , dtmin/yr

    !conditions for different tmax
    do i=1,N-1
    if(j==1) then
    v1(i)=v(i)
    d1(i)=d(i)
    t1(i)=t(i)
    p1(i)=p(i)
    end if
    if(j==2) then
    v2(i)=v(i)
    d2(i)=d(i)
    t2(i)=t(i)
    p2(i)=p(i)
    end if
    if(j==3) then
    v3(i)=v(i)
    d3(i)=d(i)
    t3(i)=t(i)
    p3(i)=p(i)
    end if
    if(j==4) then
    v4(i)=v(i)
    d4(i)=d(i)
    t4(i)=t(i)
    p4(i)=p(i)
    end if
    if(j==5) then
    v5(i)=v(i)
    d5(i)=d(i)
    t5(i)=t(i)
    p5(i)=p(i)
    end if
    enddo

enddo !end of j-cycle

    if(winds .eqv. .true.)then
        print*, 'STELLAR WINDS PHASE TERMINATED'
        write(*,*)
        print*, '---BOOOOOOM---'
        write(*,*)
        winds = .false.
        cfl = 0.01
        d0 = sum(d)/N !approximation for mean density in the sedov solution
        print*, 'MEAN DENSITY AFTER STELLAR WINDS ==', real(d0)
        call incname(nome_file)
        call incname(nome_file2)
        close(28)
        close(77)
        go to 22
    end if

    write(*,*)
    print*, 'SUPERNOVAE EFFICENCY == ', efficency*100 ,'%'
   

    close(28)
    close(77)

!========== MAIN RESULTS =================


    open(30,file='speed_winds.dat')
    write(30,*)'#--t=1e4------t=2e4--------t=4e4-------t=6e4------t=8e4-------t=1e5'
    do i=2,N-1
        write(30,1000)xa(i)/cmpc,v1(i)/1.d5,v2(i)/1.d5,v3(i)/1.d5,v4(i)/1.d5,v5(i)/1.d5
    enddo
    close(30)

    open(40,file='pressure_winds.dat')
    write(40,*)'#--t=1e4------t=2e4--------t=4e4-------t=6e4------t=8e4-------t=1e5'
    do i=2,N-1
        write(40,1000)xa(i)/cmpc,p1(i),p2(i),p3(i),p4(i),p5(i)
    enddo
    close(40)

    open(25,file='density_winds.dat')
    write(25,*)'#--t=1e4------t=2e4--------t=4e4-------t=6e4------t=8e4-------t=1e5'
    do i=2,N-1
        write(25,1000)xa(i)/cmpc,d1(i)/(1.4*mp),d2(i)/(1.4*mp),d3(i)/(1.4*mp),d4(i)/(1.4*mp),d5(i)/(1.4*mp)
    enddo
    close(25)

    open(26,file='temperature_winds.dat')
    write(26,*)'#--t=1e4------t=2e4--------t=4e4-------t=6e4------t=8e4-------t=1e5'
    do i=2,N-1
         write(26,1000)xa(i)/cmpc,t1(i), t2(i), t3(i), t4(i), t5(i)
    enddo
    close(26)



!

end program zeus

!========== subroutines =============

    SUBROUTINE BCa(z1) 
    USE DATA
    IMPLICIT NONE
    real*8, dimension (N) :: z1


  !  z1(2)=0.
  !  z1(1)=-z1(3)
  !  z1(N)=-z1(N-2)
  !  z1(N-1)=0
   z1(1)=z1(2)                  
   z1(N)=z1(N-1)

    END SUBROUTINE BCa

    SUBROUTINE BCb(z2) ! BC di outflow tradizionali
    USE DATA
    IMPLICIT NONE
    real*8, dimension (N) :: z2
    z2(1)=z2(2)
    z2(N)=z2(N-1)
    END SUBROUTINE BCb
 
    Real*8 FUNCTION Cool(Temp1)
    USE DATA
    IMPLICIT NONE
    Real*8:: Temp1, Temp_kev
    
    Temp_kev=Temp1/(1.16d7)
    
    if (Temp_kev>0.02) then
     Cool= ( 8.6d-3*(Temp_kev**(-1.7)) + 0.058*(Temp_kev**(0.5)) + 0.063 )*1.d-22
    
    else if (Temp_kev<=0.02 .and. Temp_kev>=0.0017235) then
     Cool= 6.72d-22*((Temp_kev/0.02)**(0.6))
    
    else if (Temp_kev<0.0017235) then
     Cool= 1.544d-22*((Temp_kev/0.0017235)**6)
    
    end if
    
    END FUNCTION Cool

    subroutine file(nome_file, grid, a,b,c,d,e)
        use data
    character nome_file*(*)
    real*8, dimension(N):: a,b,c,d,e, grid

    1001 format(6(1pe12.4))
  
    open(24, file=nome_file)
        do i=2, N-2
            write(24,1001)grid(i),a(i), b(i), c(i), d(i), e(i)
        enddo
    close(24)

    end subroutine file

    subroutine incname(name)

	character*7 name

	if(name(7:7) .eq. '9') go to 20
		name(7:7) = char(1 + ichar(name(7:7)))
		return

    20  if(name(6:6) .eq. '9') go to 30
		name(6:6) = char(1 + ichar(name(6:6)))
		name(7:7) = '0'
		return
    30	if(name(5:5) .eq. '9') stop
		name(5:5) = char(1 + ichar(name(5:5)))
		name(6:6) = '0'
		name(7:7) = '0'

		return
    end subroutine incname
    
    real*8 FUNCTION wsum(temp, dens) !funzione per la media pesata
    USE DATA
    IMPLICIT NONE

    real*8, dimension(N):: temp, dens
    real*8 :: num, den
    integer :: nn, nf

    nn = maxloc(dens, dim = 1) !dimensione hot bubble
    nf = maxloc(dens, dim=1)-10  !punti post-shock
    
    num=0
    den=0
    do i=nf, nn
        num = num + (temp(i)*dens(i))
        
        den = den + dens(i)
    end do
        wsum = num/den
   end function wsum