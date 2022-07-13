!******************************************************************
!  this program solves the hydrostatic equilibrium equation
!  for an isothermal gas in a NFW halo
!*****************************************************************
PROGRAM progetto1
!rhost non lo uso
PARAMETER(jmax=1500)
IMPLICIT REAL*8 (a-h,o-z)
REAL*8, DIMENSION(jmax) :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rho(jmax),mhern(jmax),rhonfw(jmax),mdark(jmax),&
        grvnfw(jmax),lnd(jmax),entot(jmax),encin(jmax),enpot(jmax),temp(jmax),&
		lndver(jmax),rhover(jmax),eq(jmax),lntemp(jmax),massgall(jmax),dgran(jmax),&
		gradzfe(jmax),zfe(jmax),rhofe(jmax),zfeoss(jmax),zfe1(jmax),zfe2(jmax),&
		zfe5(jmax),zfeSN(jmax),rho_source(jmax),sorgente_1(jmax),sorgente_2(jmax),sorgente_3(jmax),&
		T_rebusco(jmax),rho_analitic(jmax),rho_rebusco(jmax),sorgente_snfe(jmax),rhofeoss(jmax)
		
REAL*8:: f_b(2),massgas(2),rho0(2),mfepost(3),m100post(3)

REAL*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbgc,ahern,vturb,lturb,rscala,dc,tempo,&
          tmax,dt,zfein,tmaxv,mfe,m100,sorg_sn1,sn_rate_1,sorg_sn2,sn_rate_2,sorg_sn3,sn_rate_3,&
		  sorg_star,rs,cost,cost1,cost2,mfe_fin,sn_rate_fe,sorg_snfe
		  
INTEGER:: k		  
		 		  
!  constants
pi=3.14159
msol = 1.989d33
cmkpc = 3.085e21
mu=0.61
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
mbgc=(1.d12)*msol
ahern=12*cmkpc/(1+sqrt(2.))
!    set the grid

rmin = 0.*cmkpc
rmax = 3000.*cmkpc
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
enddo

do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
enddo

rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))

open(10,file='grid.dat',status='unknown')

do j=1,jmax
   write(10,*)real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo

close(10)

vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centrato a rr(j-1) !!
enddo

!  parametri del problema

rho0nfw=7.35d-26
rs=435.7*cmkpc
!creo due valori di densità iniziale in modo da considerare il gas non isotermo(1) e gas isotermo (2)
rho0(1)=1.295d-25
rho0(2)=7.320d-26
ticm=8.9e7
rvir=2797.*cmkpc
fc=1.138799
mvir=1.3e15*msol

do j=1,jmax
   x=rr(j)/rs
   rhonfw(j)=rho0nfw/(x*(1.+x)**2) !centrato in rr
enddo

open(20,file='masse.dat')

mnfw(1)=0.

do j=2,jmax
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j) ! centrate in r, massa con NFW
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc ! massa dark matter analitica
   write(20,1001)log10(r(j)/cmkpc),log10(mnfw(j)/msol),log10(mdark(j)/msol)
enddo

1001 format(3(1pe12.4))
close(20)

do j=1, jmax
   massgall(j)= mbgc*(r(j)**2)/(r(j)+ahern)**2 ! centrata in r
end do

open(20,file='grv.dat')

grvnfw(1)=0.          !! ok per alone NFW, isotermo o beta-model
do j=2,jmax
   grvnfw(j)=guniv*(mnfw(j)+massgall(j))/r(j)**2 ! centrata in r
   write(20,1002)r(j)/cmkpc,grvnfw(j)/msol
enddo

1002 format(2(1pe12.4))
close(20)

! calcolo il profilo di temperatura centrato in rr come la densità
do j=1, jmax
  eq(j)= rr(j)/(1400.*cmkpc)
  temp(j)= ticm*1.35*((((eq(j)/0.045)**(1.9))+0.45)/(((eq(j)/0.045)**(1.9))+1))*(1/((1+(eq(j)/0.6)**(2))**(0.45)))
  lntemp(j)=log10(temp(j))
end do


open(30,file='temp.dat',status='unknown')
do j=1,jmax
   write(30,1003) log10(rr(j)/cmkpc),lntemp(j),temp(j)
end do
close(30)
1003 format(3(1pe12.4))


!   calculate the gas density, assuming ticm, and temp variable
lndver(1)=log(rho0(1)) !rho realistica
lnd(1)=log(rho0(2))    !rho nel caso isotermo

rho(1)=exp(lnd(1))
rhover(1)=exp(lndver(1)) 

       !! mette il gas in eq. con il potenziale, densità centrate in rr
do j=2,jmax
   gg=grvnfw(j)
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm) 
   lndver(j)=lndver(j-1)-((rr(j)-rr(j-1))*gg*(mu*mp))/(boltz*temp(j))-(log(temp(j))-log(temp(j-1)))

   rho(j)=exp(lnd(j))
   rhover(j)=exp(lndver(j))
enddo

open(40,file='density.dat',status='unknown')
do j=1,jmax
   write(40,1000) log10(rr(j)/cmkpc), log10(rho(j)), log10(rhover(j)), log10(rhonfw(j))
enddo
close(40)
1000 format(4(1pe12.4))



!ora mi calcolo il profilo della temperatura e della densità ANALITICO assunto nell'articolo Rebusco et al., per fare la comparazione

cost=1.9d-24
cost1=4.6d-2
cost2=4.8d-3

do j=1,jmax

   t1=rr(j)/(57.*cmkpc)

   t2=rr(j)/(200.*cmkpc)

   t3=rr(j)/(71.*cmkpc)

   rho_rebusco(j)=(cost*cost1)/(1+t1**2)**1.8 + (cost*cost2)/(1+t2**2)**0.87

   T_rebusco(j)=((7.+7.*t3**3)/(2.3+t3**3))*1.2d7
   
end do

open(50,file='Rebusco.dat',status='unknown')

do j=1,jmax

   write(50,*) log10(rr(j)/cmkpc), log10(rho_rebusco(j)),log10(T_rebusco(j))

enddo
close(50)


!ora calcolo il profilo di densità ANALITICO mostrato nelle slide

b=(8*3.14*guniv*mu*mp*rho0nfw*(rs**2))/(27*boltz*ticm)

open(60,file='d_analitic.dat',status='unknown') 
do j=1,jmax
   
   x=rr(j)/rs
   rho_analitic(j)=rho0(2)*(exp(-27*b*0.5))*(1+x)**(27*b*0.5/x)

   write(60,*) log10(rr(j)/cmkpc), log10(rho_analitic(j))

end do
close(60)


!calcolo la massa del gas dell'ammasso
  massgas=(0,0)
do j=1, jmax
  massgas(1) = massgas(1) + rhover(j)*vol(j)
  massgas(2) = massgas(2) + rho(j)*vol(j)
end do
!valori di f utili per correggere le densita iniziali del gas nel caso realistico (1) e isotermo (2)
f_b(1)=(massgas(1)+massgall(jmax))/(massgall(jmax)+massgas(1)+mdark(jmax))
f_b(2)=(massgas(2)+massgall(jmax))/(massgall(jmax)+massgas(2)+mdark(jmax))


! printo i valori di f_b per stabilire se ho scelto correttamente il valore iniziale di rho_0
print*, "Per i valori di densita' iniziale scelti ottengo i seguenti valori di frazione barionica"
print*,f_b(1),"per il gas con T(r) e", f_b(2), "per quello isotermo"
         
write(*,*)


!"Ora passo alla diffusione della metallicita"
vturb=260.d5
lturb=15.*cmkpc
rscala=30.*cmkpc
dc=0.11*vturb*lturb

do j=1,jmax  
	 dgran(j)=dc-(0.6*dc)*exp(-((r(j)/rscala)**2))	
end do

!distribuzione del ferro osservata
do j=1,jmax
    zfe(j)=0.3*((2.2+(rr(j)/(cmkpc*80))**3)/(1+(rr(j)/(cmkpc*80))**3))*0.0018*1.4
	rhofe(j)=zfe(j)*rhover(j)/1.4
	zfeoss(j)=zfe(j)
	rhofeoss(j)=rhofe(j)
end do

! calcolo della massa totale per controllarne la conservazione  
   mfe=0
do j=1,jmax
    
    mfe=mfe+vol(j)*(zfeoss(j))*rhover(j)/1.4
	
	if(r(j)<100*cmkpc) then
	   m100=mfe
	endif

enddo

!calcolo della diffusione senza termini di sorgente a partire dalla zfe osservata come abbondanza iniziale
tmax=5*3.16d16
time=0.
dt=0.5*(r(5)-r(4))**2/(2.*dc)

k=0
do while(time<tmax) ! le densità del ferro sono centrate in rr
   
    do j=2,jmax-1
       gradzfe(j)=(zfe(j)-zfe(j-1))/(rr(j)-rr(j-1))
    enddo
    gradzfe(1)=0.
    gradzfe(jmax)=0.

    do j=2,jmax-1
      rhoj1=0.5*(rhover(j+1)+rhover(j))
      rhoj=0.5*(rhover(j-1)+rhover(j))
	  
      fj1=r(j+1)**2*dgran(j+1)*rhoj1*gradzfe(j+1)
      fj=r(j)**2*dgran(j)*rhoj*gradzfe(j)	  
	  dem=(r(j+1)**3-r(j)**3)/3
	  
	  rhofe(j)=(rhofe(j)+(dt/1.4)*(fj1-fj)/dem)
	  zfe(j)=1.4*rhofe(j)/rhover(j)
    enddo 
	
	zfe(1)=zfe(2)
    zfe(jmax)=zfe(jmax-1)
    rhofe(1)=rhofe(2)
    rhofe(jmax)=rhofe(jmax-1)
	
    time=time+dt
	
	!salvo i dati relativi a tre diverso tempi (1-2-5 Gyr)
	if (time>5*3.16d16 .and. k==2 ) then
	     zfe5=zfe
	     do j=1,jmax
    
           mfepost(3)=mfepost(3)+vol(j)*(zfe5(j))*rhover(j)/1.4
	
	       if(r(j)<100*cmkpc) then
	         m100post(3)=mfepost(3)
	       endif

         enddo
		 
		 k=3
		 
    else
	    if(time>2*3.16d16 .and. k==1) then
		     zfe2=zfe
			 do j=1,jmax
    
                mfepost(2)=mfepost(2)+vol(j)*(zfe2(j))*rhover(j)/1.4
	
	            if(r(j)<100*cmkpc) then
	              m100post(2)=mfepost(2)
	            endif

             enddo
			 
			 k=2
			 
	    else
			if(time>3.16d16 .and. k==0)then
			     zfe1=zfe
			     do j=1,jmax
    
                    mfepost(1)=mfepost(1)+vol(j)*(zfe1(j))*rhover(j)/1.4
	
	                if(r(j)<100*cmkpc) then
	                    m100post(1)=mfepost(1)
	                endif

                 enddo
				 
			     k=1
				 
			endif		 
		endif
    endif	
	
enddo


mfe_fin=0
do j=1,jmax
    
    mfe_fin=mfe_fin+vol(j)*(zfe(j))*rhover(j)/1.4
	
enddo

! printo per verificare che la massa si conservi 
print*,"Controllo della conservazione della massa nell'utilizzo del codice diffusivo"
print*, "Massa totale iniziale",mfe/(msol*(10**8)),"e finale",mfe_fin/(msol*(10**8)), "in unita di 10^8 Msol"
write(*,*)
write(*,*)
!printo la massa entro 100 kpc per vedere se la massa venga diiffusa al di fuori di tale distanza
print* ,"Massa entro 100 kpc all'inizio",m100/(msol*(10**8)), "in unita di 10^8 Msol"
print*, "e dopo 1-2-5 Gyr", m100post/(msol*(10**8)), "in unita di 10^8 Msol"
write(*,*)
write(*,*)
open(70,file="ferrodiff.dat") !scrivo i dati relativi alla sola difusione
do j=1,jmax
   write(70,1500) log10(rr(j)/cmkpc),zfeoss(j)/0.0018,zfe1(j)/0.0018,zfe2(j)/0.0018,zfe5(j)/0.0018, zfe(j)/0.0018
   
enddo

1500 format(6(1pe12.4))
close(70)


!considero ora l'evoluzione del ferro solo trammite il contributo delle supernove con zfe iniziale = 0
zfeSN=0

!inserisco i tre rate di supernova
sn_rate_1=0.15
sn_rate_2=0.3
sn_rate_3=0.5

!sotto sono i tre termini da inserire nella equazione differenizale
sorg_sn1=(0.74/7.5)*(1/1.d10)*(1/3.156d9)*sn_rate_1
sorg_sn2=(0.74/7.5)*(1/1.d10)*(1/3.156d9)*sn_rate_2
sorg_sn3=(0.74/7.5)*(1/1.d10)*(1/3.156d9)*sn_rate_3


!ora la sorgente

sorg_star=6.d-23 !sempre costante

! calcolo i termini sorgente per diversi rate di supernova
do j=1, jmax

   rho_source(j)=(mbgc/(2.*pi))*(ahern/rr(j))*(1./((rr(j)+ahern)**3))
   
   sorgente_1(j)=rho_source(j)*(sorg_star+sorg_sn1)
   sorgente_2(j)=rho_source(j)*(sorg_star+sorg_sn2)
   sorgente_3(j)=rho_source(j)*(sorg_star+sorg_sn3)
   
end do

time=0.
rhofe=0

! ora vedo come si forma il picco del ferro a partire dalle sole sorgenti con sn rate= 0.15, per i plot ricalcolo cambiando il valore
do while(time<tmax)

   time=time+dt !è un contatore
 
   !ora calcolo lo step temporale della densità del ferro


    do j=2,jmax-1 !è in cgs tutto

      rhofe(j)=rhofe(j) + dt*sorgente_1(j)
      zfe(j)=(1.4*rhofe(j))/(rhover(j))  

    end do
	
   !condizione al contorno
   rhofe(1)=rhofe(2)
   rhofe(jmax)=rhofe(jmax-1)  
   zfe(1)=zfe(2)
   zfe(jmax)=zfe(jmax-1)


end do
!FINE DEL CICLO while temporale.

mfe_fin=0.

do j=1,jmax
	
    mfe_fin=mfe_fin+vol(j)*rhofe(j)
	
enddo
write(*,*)
write(*,*)
print*,"La massa di ferro prodotta con un SN_rate = ", sn_rate_1,"e'",mfe_fin/(msol*(10**8)), "in unita di 10^8 Msol"
print*, "in confronto quella osservata e' di",mfe/(msol*(10**8)), "in unita di 10^8 Msol"
write(*,*)
write(*,*)

open(80, file='zfe_sorgente1.dat', status='unknown')
do j=1,jmax
   write(80,*) log10(rr(j)/cmkpc), zfe(j)/0.0018
end do
close(80)




!mi occupo ora di studiare il modello con diffusione+sorgente in cui ridefinisco i parametri per riprodurre il picco osservato

!inizalizzo le variabile del problema
time=0.
k=0



zfe=0
zfe1=0
zfe2=0
zfe5=0
rhofe=0
dc=0.15*vturb*lturb
dt=0.5*(r(5)-r(4))**2/(2.*dc)

do j=1,jmax  
	 dgran(j)=dc-(0.6*dc)*exp(-((r(j)/rscala)**2))	
end do



do while(time<tmax)

   time=time+dt
   
   ! caso più relistico con CN rate variabile nel tempo secondo Rebusco
   sn_rate_fe=0.038*((time/(13.6*3.16d16))**(-1.1))!sn rate che produce la quantità di ferro osservata entro il tempo di evoluzione

   sorg_snfe=(0.74/7.5)*(1/1.d10)*(1/3.156d9)*SN_rate_fe

   do j=1, jmax
   
     sorgente_snfe(j)=rho_source(j)*(sorg_star+sorg_snfe)
   
   end do
   
   
   
   !calcolo il gradiente di zfe ad ogni ciclo
   do j=2,jmax-1
      gradzfe(j)=(zfe(j)-zfe(j-1))/(rr(j)-rr(j-1)) 
   enddo

   gradzfe(1)=0.
   gradzfe(jmax)=0.

   !riscivo lo stesso algoritmo per la diffusione aggiungendo il termine di sorgente1
   do j=2,jmax-1 

      rhoj1=0.5*(rhover(j+1)+rhover(j))
      rhoj=0.5*(rhover(j-1)+rhover(j))
	  
      fj1=r(j+1)**2*dgran(j+1)*rhoj1*gradzfe(j+1)
      fj=r(j)**2*dgran(j)*rhoj*gradzfe(j)	  
	  dem=(r(j+1)**3-r(j)**3)/3
	  
	  rhofe(j)=(rhofe(j)+(dt/1.4)*(fj1-fj)/dem)+(sorgente_snfe(j)*dt)
	  zfe(j)=1.4*rhofe(j)/rhover(j)
      	 
   end do
    
   !condizioni al contorno
   rhofe(1)=rhofe(2)
   rhofe(jmax)=rhofe(jmax-1)  
   zfe(1)=zfe(2)
   zfe(jmax)=zfe(jmax-1)
   
   
   
   !salvoi dati a tre tempi tipici
	 if (time>5*3.16d16 .and. k==2) then
	  
	    zfe5=zfe
		k=3
	   
     else
	    if(time>2*3.16d16 .and. k==1) then
		
		    zfe2=zfe
		    k=2
			
	    else
			if(time>3.16d16 .and. k==0)then
			
			   zfe1=zfe
			   k=1
			  
			endif
		endif
     endif	
     
end do

! calcolol la massa totale prodotta del picco per confrontarla con quella osservata senza il platò
mfe_fin=0
mfe=0

do j=1,jmax

    if(r(j)<100.*cmkpc) then
	
    mfe_fin=mfe_fin+vol(j)*rhofe(j)
	
	end if
enddo

!la osservata riguarda solo quella prodotta dalle SN 1a cioè senza il contributo del platò legato alle SN 2
do j=1,jmax
    
	if(r(j)<100.*cmkpc) then
	   
    mfe=mfe+vol(j)*(zfeoss(j)-(0.42*0.0018))*rhover(j)/1.4
	
	end if
	
enddo

!calcolo la massa prodotta dal mio codice entro la zona del picco e la confronto con quella osservata
!per avere un indice per calibrare i parametri
print*, "Controllo di riprodurrre la stessa quantita' di ferro nella parte diff.+sorg."
write(*,*)
print*, "Massa di ferro osservata entro 100 Kpc",mfe_fin/(msol*(10**8)), "in unita di 10^8 Msol"
write(*,*)
print*, "Massa di ferro prodotta entro 100 Kpc", mfe/(msol*(10**8)), "in unita di 10^8 Msol"


open(90, file='zfe_FIN.dat', status='unknown')

do j=1,jmax
write(90,1600) log10(rr(j)/cmkpc), (zfe1(j)/0.0018), (zfe2(j)/0.0018), (zfe5(j)/0.0018), &
                         (zfe(j)/0.0018), (zfeoss(j)/0.0018)-0.42
end do   
1600 format(6(1pe12.4))
close(90)

END PROGRAM progetto1