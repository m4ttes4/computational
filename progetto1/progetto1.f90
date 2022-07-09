
MODULE DATA
   
   real*8 :: pi,cmpc,cmkpc,yr,boltz,mu,mp, guniv,msol
   parameter (guniv=6.6720e-8)
   parameter(pi=3.141592)
   parameter(cmpc=3.085d18)
   parameter(cmkpc=1000.*cmpc)
   parameter(yr=3.156d7)
   parameter(boltz=1.38d-16)
   parameter(mu=0.61)
   parameter(mp=1.67d-24)
   parameter(msol= 1.989d33)
   parameter(jmax=2500) 
END MODULE DATA




program progetto1
   use DATA



implicit real*8 (a-h,o-z) !ogni variabile dalla a alla h e o alla z sono real*8 !! i l m n
real*8, dimension(jmax) :: r(jmax),rr(jmax),vol(jmax),massaDM_num(jmax),&
        dens_gas(jmax),densDM(jmax),massaDM_ana(jmax),&
        grvnfw(jmax),lnd(jmax), dens_gas_ana(jmax), M_star(jmax), &
        pot_gal(jmax), M_tot(jmax), f_b(jmax), M_gas(jmax), num1(jmax), den1(jmax),&
        den2(jmax), T(jmax), h(jmax), dens_gas_old(jmax),dens_gas_gal(jmax),grav2(jmax),&
        lnd2(jmax), lnd3(jmax), dens_gas_temp(jmax), f_b2(jmax), M_gas2(jmax), dens_gas_ana2(jmax),&
        M_gas_temp(jmax), f_b3(jmax)

real*8 :: rmin,rmax,mvir,rvir, Mgal, r_half, a, Mgas, num, b2,Tcentr, kappa, lturb, k

real*8, dimension(jmax) :: ne(jmax),dens_gas_paper(jmax),temp_paper(jmax),zfe(jmax), dens_fe_paper(jmax),&
         Zfe_paper(jmax),M_fe(jmax), dens_star(jmax), gradzfe(jmax), dens_fe_diff(jmax), zfe_diff(jmax), &
         source(jmax), Zfe_source(jmax), dens_fe_source(jmax),dens_fe_source1(jmax),Zfe_source1(jmax),Zfe_paper_sub(jmax),&
         M_fe_tot(jmax), rho_fe(jmax), zfe_sorg_diff(jmax),zfe_time(jmax), source_time(jmax), rho_fe_time(jmax), &
         rho_fe_5gyr(jmax), zfe_1gyr(jmax),rho_fe_1gyr(jmax), zfe_5gyr(jmax),alpha_1a(jmax), kappa2(jmax)
      
integer :: n, nmax
                          
                           
!========= some paramenter of the problem =================0
   densDM0=7.35d-26
   rs=435.7*cmkpc
   dens_gas0=2.882d-26
   temp=8.9e7
   rvir=2797.*cmkpc
   fc=1.138799
   mvir=1.3e15*msol
   Mgas=1.d13*msol
   Mgal=1.d12*msol
   r_half=12*cmkpc
   aher=r_half/(1.+sqrt(2.))

   Zfe_sol=1.8d-3
   Zfe_out=0.4*Zfe_sol
   zfesn=0.744/1.4

   1003 format(1f8.3, 1pe12.4) !format for bar frac on terminal

   1005 format(i5, i10) !format  for i-cycle on terminal

!
   write(*,*) 'per quanti Gyr si vuole far evolvere il sistema?'
   read(*,*) k

   tmax=k*3.154d16  

   1008 format(1pe12.4)



!########## CREAZIONE GRGLIA ###############

   rmin = 0.*cmkpc
   rmax = 3000.*cmkpc
   do j=1,jmax
      r(j)=rmin+(j-1)*rmax/(jmax-1) !ottiene una griglia 
   enddo
   do j=1,jmax-1
      rr(j)=r(j)+0.5*(r(j+1)-r(j)) !questi sono i punti sulla griglia a distanza mezza per rendera piu densa
   enddo
   rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2)) 

   !PROBLEMA A SIMMETRIA SFERICA
   vol(1)=4.1888*r(1)**3 !volume di sfera, 4.1888 è quattro terzi pigraco
   do j=2,jmax
      vol(j)=4.1888*(r(j)**3-r(j-1)**3)  !VOLUME TRA DUE SHELL  !! centrato a rr(j-1) !! definisce le shell sferiche su cui andra a calcolare masse e desnsita
   enddo


!############ DISTRIBUZIONE TEORICA DENSITA DM ##############
   do j=1,jmax       
      x=rr(j)/rs
      densDM(j)=densDM0/(x*(1.+x)**2) 
   enddo

!######### SOLUZIONE NUMERICA E ANALITICA MASSA DM #######################

   massaDM_num(1)=0.
   massaDM_ana(j)=0.
   do j=2,jmax
      x=r(j)/rs
      massaDM_num(j)=massaDM_num(j-1)+densDM(j-1)*vol(j) !massa in una shell sferica soluzione num
      massaDM_ana(j)=mvir*(log(1.+x)-x/(1.+x))/fc !massa DM da soluzione analitica
   enddo


!############ PROFILI DI TEMPERATURA ############

   do j=1, jmax
      h(j)=((rr(j)/(1.4e3*cmkpc))/0.045)
      num1(j)=h(j)**(1.9)+0.45
      den1(j)=h(j)**(1.9)+1
      den2(j)=(1+(((h(j)*0.045)/0.6)**2))**(0.45)

      T(j)=temp*1.35*(num1(j)/den1(j))*(1/den2(j)) !temp profile for dens_gas_temp

      !valori paper per diffusione fe per confronto
      temp_paper(j)=7*(1+((rr(j)/cmkpc)/71)**3)/(2.3+((rr(j)/cmkpc)/71)**3)*(1.16d7) !profilo temp paper rebusco
   enddo

   !profilo teorico massa stelle 
   do j=1, jmax
      M_star(j)=Mgal*r(j)**2/((r(j)+aher)**2)
      M_tot(j)=massaDM_num(j)+M_star(j)
   enddo

   !per plot di confronto
   open(30, file='massaDM.dat')
      do j=1,  jmax
         write(30,1006)rr(j)/cmkpc, massaDM_ana(j)/msol, massaDM_num(j)/msol, M_star(j)/msol
      enddo
   close(30)
   1006 format(4(1pe12.4))

!################ CALCOLO DEI PROFILI DI DENSITA #############

   grvnfw(1)=0.
   grav2(1)=0.
   do j=2,jmax
      grvnfw(j)=guniv*massaDM_num(j)/r(j)**2 
      grav2(j)=guniv*M_tot(j)/r(j)**2  !per quando metto la galassia     
   enddo

   dens_gas0=4.22d-26  !original
   dens_gas1=8.03d-26  !with bcg
   dens_gas2=1.489d-25  !with bcg and dt/dr

   11 continue
   lnd(1)=log(dens_gas0)    
   lnd2(1)=log(dens_gas1) 
   lnd3(1)=log(dens_gas2)    
   do j=2,jmax
      gg=grvnfw(j)   
      lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*temp) !original

      gg2=grav2(j)
      lnd2(j)=lnd2(j-1)-grav2(j)*(mu*mp)*(rr(j)-rr(j-1))/(boltz*temp)  !per quando metto la galassia

      Tcentr=0.5*(T(j)+T(j-1))                                          !temp centrata a r 
      lnd3(j)=lnd3(j-1)-grav2(j)*(mu*mp)*(rr(j)-rr(j-1))/(boltz*Tcentr)-(log(T(j))-log(T(j-1)))  !con nuovo profilo di temp
   enddo


   !b=(8*pi*guniv*mu*mp*dens_gas0*rs**2.)/(27.*boltz*temp)
   !A=exp(-27.*b/2.)
   !cost=dens_gas0*A

   !dens_gas_ana(1)=dens_gas0
!====================  CALCOLO DELLE DENSITA + PROFILO FERRO ============================
   do j=1,jmax

      !a=r(j)/rs    
      !dens_gas_ana(j)=cost*((1.+a)**(27.*b/(2.*a)) ) !densità gas analitica

      dens_gas(j)=exp(lnd(j))        !original
      dens_gas_gal(j)=exp(lnd2(j))   !with BCG
      dens_gas_temp(j)=exp(lnd3(j))  !with temp profile

      dens_star(j)=Mgal/(2.*pi)*(aher/rr(j))/(aher+rr(j))**3  !hernquish profile for fe diffusion

      ne(j)=4.6d-2/((1+((rr(j)/cmkpc)/57)**2)**1.8) +4.8d-3/(1+((rr(j)/cmkpc)/200)**2)**0.87  !densità elettronica paper rebusco
      dens_gas_paper(j)=1.937d-24*ne(j)     !densità gas paper rebusco

      xa=rr(j)/(80.*cmkpc)
      zfe_paper(j)=Zfe_sol*0.3*1.4*1.15*(2.2+xa**3)/(1+xa**3)/1.15  !Perseus!
      Zfe_paper_sub(j)=Zfe_paper(j) - Zfe_out   !! subtract z_Fe,out

      !Zfe_paper(j)=0.3*(2.2+((rr(j)/cmkpc)/80)**3)/(1+(((rr(j)/cmkpc)/80)**3))*Zfe_sol !abbondanza osservata perseo
      !Zfe_paper_sub(j)=Zfe_paper(j)-Zfe_out  !setto condizioni iniziali sottraendo eccesso
      dens_fe_paper(j)=dens_gas_temp(j)*Zfe_paper(j)/1.4  !densità fe rebusco !da sostituire defe fe papaer per plot di confronto
   enddo


   !MASSA GAS
   M_gas(1)=vol(1)*dens_gas0
   M_gas2(1)=vol(1)*dens_gas1
   M_gas_temp(1)=vol(1)*dens_gas2
   do j=2, jmax
      M_gas(j)=M_gas(j-1)+dens_gas(j-1)*vol(j)
      M_gas2(j)=M_gas2(j-1)+dens_gas_gal(j-1)*vol(j) !per la galassia
      M_gas_temp(j)=M_gas_temp(j-1)+dens_gas_temp(j-1)*vol(j) !with temp profile
   enddo


   !FRAZIONI BARIONICHE
   write(*,*)'valori delle frazioni barioniche per controllo'
   write(*,*)'--f_b------',        'density'
  
   write(*,1003) f_b(2490) , dens_gas0 !2490 is r_vir
   write(*,1003) f_b2(2490), dens_gas1
   write(*,1003) f_b3(2490), dens_gas2


   !ciclo di controllo sulle frazioni barioniche
   if(f_b3(2490).ge.0.17)then
      dens_gas2=dens_gas2-1.d-27
      go to 11
   else if(f_b3(2490).le.0.17)then
   end if

   if(f_b2(2490).ge.0.17)then
      dens_gas1=dens_gas1-1.d-27
      go to 11
   else if(f_b2(2490).le.0.17)then
   end if

   if(f_b(2490).ge.0.17)then
      dens_gas0=dens_gas0-1.d-27
      go to 11
   else if(f_b(2490).le.0.17)then
   end if


   write(*,*) !spazio sul term 


   !salvo copie vecchi dati
   !dati sulle densità
   open(10, file='densit.dat')
   write(10,*)'#'
   do j=1, jmax
    write(10,1000)rr(j)/cmkpc,densDM(j), dens_gas(j), dens_gas_temp(j), dens_gas_gal(j)
   enddo
   close(10)
   1000 format(5(1pe12.4))

   !per confronto con paper rebusco
   open(30, file='temp_paper.dat')
   do j=1, jmax
      write(30,1000)rr(j)/cmkpc, T(j), temp_paper(j), dens_gas_paper(j), dens_gas_temp(j)
   enddo
   close(30)

   open(30, file='paper.dat')
   do j=1, jmax
    write(30,1000)rr(j)/cmkpc, Zfe_paper(j)/Zfe_sol, Zfe_paper_sub(j)/Zfe_sol
   enddo
   close(30)

!############################################################
!############### PARTE DELLA DIFFUSIONE #####################

  ! tmax = (/2.d9*yr, 3.d9*yr, 5.d9*yr, 8.d9*yr, 1.d10*yr/)

   vturb=260.e5   
   lturb=15.*cmkpc  
   rscala=30.*cmkpc
   kappa=0.11*vturb*lturb   !diffusion term

   dt=((r(5)-r(4))**2/(2*kappa))

   !calcolo massa Fe iniziale
 
   write(*,*)'massa ferro tot prima di evoluzione='
   call massa(dens_fe_paper,vol)
  

!############### evoluzione temporale con sorgente ############

   !condizioni iniziali + termine sorgente
   do j=1, jmax
      Zfe_source(j)=0
      dens_fe_source(j)=0
      source(j)=dens_star(j)*(6.d-23+4.7d-22)
   enddo
   

   tempo=0
   do while(tempo.le.tmax)

      do j=2, jmax-2
         dens_fe_source(j)=dens_fe_source(j)+dt*source(j) 
         Zfe_source(j)=1.4*(dens_fe_source(j)/dens_gas_temp(j))
      enddo

      tempo=tempo+dt
      
      call BCb(dens_fe_source)
      call BCb (Zfe_source)
   enddo


   !scrittura risultati
   open(30, file='sorgente.dat')
      do j=1, jmax
         write(30, 1002) rr(j)/cmkpc, Zfe_paper_sub(j)/Zfe_sol, Zfe_paper(j)/Zfe_sol, Zfe_source(j)/Zfe_sol
      enddo
   close(30)
   1002 format (5(1pe12.4))



!################# DIFFUSIONE FERRO NO SORGENTE ##########################

   !definisco variabili su cui lavoro
   do j=1, jmax
      zfe_diff(j)=Zfe_paper(j)
      dens_fe_diff(j)=dens_fe_paper(j)
   enddo

   !resetto tempo

   tempo=0
   do while(tempo.le.tmax)

      tempo=tempo+dt
      
      do j=2, jmax-1
         gradzfe(j)=(zfe_diff(j)-zfe_diff(j-1))/(rr(j)-rr(j-1))  !calcolo gradiente
      enddo
      gradzfe(1)=0
      gradzfe(jmax)=0

      do j=2, jmax-2
         densp1=0.5*(dens_gas_temp(j+1)+dens_gas_temp(j))
         dens=0.5*(dens_gas_temp(j-1)+dens_gas_temp(j))

         dens_fe_diff(j)=dens_fe_diff(j)&
         +(3*dt/1.4)*((r(j+1)**2*kappa*densp1*gradzfe(j+1))&
         -(r(j)**2*kappa*dens*gradzfe(j)))/(r(j+1)**3-r(j)**3)

         zfe_diff(j)=1.4*(dens_fe_diff(j))/dens_gas_temp(j)    
      enddo

      call BCb(zfe_diff)
      call BCb(dens_fe_diff)

   enddo

   !sottraggo eccesso
   do j=1, jmax
   zfe_diff(j)=zfe_diff(j)-Zfe_out
   enddo

   !controllo onservazione ferro totale
   write(*,*)'massa ferro tot prima di evoluzione='
   call massa(dens_fe_diff,vol)
   


   open(33, file='diffusion.dat')
   do j=1, jmax
      write(33,1001)rr(j)/cmkpc, Zfe_paper_sub(j)/Zfe_sol, zfe_diff(j)/Zfe_sol, zfe_paper(j)/Zfe_sol
   enddo
   1001 format(5(1pe12.4))
   close(33)

!################### DIFFUSIONE PIU' SORGENTE ############################

   !nuove condizioni iniziali
   do j=1, jmax
      zfe(j)=0
      rho_fe(j)=0
   enddo


   tempo=0 

   do while(tempo.le.tmax)
      
      do j=2, jmax-1
      gradzfe(j)=(zfe(j)-zfe(j-1))/(rr(j)-rr(j-1))
      enddo
      gradzfe(1)=0
      gradzfe(jmax)=0

      do j=2, jmax-2
      rho_fe(j)=rho_fe(j)+dt*source(j)
      enddo

      do j=2, jmax-2
      densp1=0.5*(dens_gas_temp(j+1)+dens_gas_temp(j))
      dens=0.5*(dens_gas_temp(j-1)+dens_gas_temp(j))

      rho_fe(j)=rho_fe(j)&
      +(3*dt/1.4)*((r(j+1)**2*kappa*densp1*gradzfe(j+1))&
      -(r(j)**2*kappa*dens*gradzfe(j)))/(r(j+1)**3-r(j)**3)

      zfe(j)=1.4*(rho_fe(j)/dens_gas_temp(j))
      enddo

      call BCb(rho_fe)
      call BCb(zfe)
      
      tempo=tempo+dt

   enddo

   open(33, file='diffusion_source.dat')
   do j=1, jmax
      write(33,1001)rr(j)/cmkpc, Zfe_paper_sub(j)/Zfe_sol, zfe(j)/Zfe_sol
   enddo


   !controllo massa ferro per paragone, non più conservazione
   write(*,*)'massa con diff e sorg'
   call massa(rho_fe,vol)

!############################################
   !TUTTO OK, LA SORGENTE CON LA DIFFUSIONE NON RIPRODUCE LE OSSERVAZIONI, CAMBIO DATI INIZIALI
   !prima cosa, modifico termine di sorgente

   !go to 20
   snu=0.25  !supernovae rate a t=tnow
   

   vturb=260.e5   
   lturb=25.*cmkpc  
   rscala=30.*cmkpc
   kappa=vturb*lturb   !diffusion term
   dt=((r(5)-r(4))**2/(2*kappa))

   !trovo una formula per kappa variabile


   tnow=1e9*yr !10Gyr

   !ridefinisco variabili di lavoro
   do j=1, jmax
      zfe(j)=0
      rho_fe(j)=0
   enddo


!write(*,*)kappa

   tempo=0 

   do while(tempo.le.tmax)

      tempo=tempo+dt

      ncicli=ncicli+1
      if(ncicli.ge.100)then
      ncicli=0
       write(*,1008)snu_time
       end if
      
      do j=1, jmax
         snu_time=snu*(1.17*exp((((-tempo)/tnow-0.261)**1.1)/(2*0.276)**2)) !supernovae rate variable in time
      
         alphast=4.7e-20
         alphasn=4.436e-20*(snu_time/7.4)

         source_time(j)=dens_star(j)*(alphast*Zfe_sol/1.4+alphasn*zfesn)
      enddo
      
      do j=2, jmax-1
         gradzfe(j)=(zfe_time(j)-zfe_time(j-1))/(rr(j)-rr(j-1))
      enddo
      gradzfe(1)=0
      gradzfe(jmax)=0 

      do j=2, jmax-1
         rho_fe_time(j)=rho_fe_time(j)+dt*source_time(j)
      enddo

      do j=2, jmax-1
      densp1=0.5*(dens_gas_temp(j+1)+dens_gas_temp(j))
      dens=0.5*(dens_gas_temp(j-1)+dens_gas_temp(j))

      rho_fe_time(j)=rho_fe_time(j)&
      +(3*dt/1.4)*((r(j+1)**2*kappa*densp1*gradzfe(j+1))&
      -(r(j)**2*kappa*dens*gradzfe(j)))/(r(j+1)**3-r(j)**3)

      zfe_time(j)=1.4*(rho_fe_time(j)/dens_gas_temp(j))
      enddo
      call BCb(rho_fe_time)
      call BCb(zfe_time)
      
      
   enddo

   open(33, file='snu.dat')
   do j=1, jmax
      write(33,1001)rr(j)/cmkpc, Zfe_paper_sub(j)/Zfe_sol, zfe_time(j)/Zfe_sol
   enddo

   close(33)
   !20 continue

   write(*,*)'massa con metodo matlab'
 call massa(rho_fe_time,vol)


end program progetto1

SUBROUTINE BCb(z2) ! BC di outflow tradizionali
   USE DATA
   IMPLICIT NONE
   REAL*8, DIMENSION (jmax):: z2
   
   z2(1)=z2(2)
   z2(jmax)=z2(jmax-1)
END SUBROUTINE BCb

subroutine massa(a1, b1)
   use data
   implicit none 
   integer :: j
   real*8, dimension(jmax):: a1, b1, mass
   real*8 :: result
  
   do j=2, jmax
   mass(j)=a1(j-1)*b1(j)
   enddo

   result=sum(mass)

    write(*,89) real(result/msol)
    89 format(1pe12.4)
   end subroutine massa