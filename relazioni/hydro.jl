
function main()
#definizione de parametri
N = 2000
cmpc = 3.085e18
mu = 1.4
msol = 1.989e33
mp =  1.67e-24
cv = 2e8
gam = 1.667
yr = 3.156e7



xa = zeros(N)
xb = copy(xa)


#grid creation
xmin = 0
xmax = 200*cmpc


xa[1] = xmin +(xmax-xmin)*(-1) /(N-1)
xa[2] =0
for i = 1:N
    xa[i] = xmin+(xmax-xmin)*(i-2)/(N-1)
end
xa[1] = -xa[3]

deltax = xa[5]-xa[4]


for i=2:N-1
    xb[i] = 1/2 * (xa[i]+xa[i+1])
end
xb[N] = xb[N-1]+(xb[N-1]-xb[N-2])




#fattori di scala
g2a = zeros(N)
g2b = rand(N)
g31a = rand(N)
g31b = rand(N)
dvl1a = rand(N)
dvl1b = copy(g2a)

sdr = 1
if sdr==1  
    for i=1:N
        g2a[i]=xa[i]
        g31a[i]=xa[i]
        g2b[i]=xb[i]
        g31b[i]=xb[i]
    end 

    for i=1:N-1
        dvl1a[i]=(xa[i+1]^3-xa[i]^3)/3.
    end 
    dvl1a[N]=dvl1a[N-1]
    for i=2: N
        dvl1b[i]=(xb[i]^3-xb[i-1]^3)/3.
    end 
    dvl1b[1]=dvl1b[2]
end


#----------- fattori supernova e condizioni iniziali --------
e0 = 1e51
c2 = 3
cfl = 0.01

d0 = 2e-24
t0 = 1e4
t0_min = t0
n0=d0/(mu*mp)
vol = 1.333*pi*xa[4]^3

v = zeros(N)
d = copy(v)
e = copy(v)
t = copy(v)
p = copy(v)
q = copy(v)
divV = copy(v)

for i=1:N
    v[i] = 0
    d[i] = d0
    t[i] = t0
    e[i] = cv*d[i]*t[i]
    p[i] = (gam -1)*e[i]
end


for i=1:3
    e[i] = e0/vol
    p[i] = e[i]*(gam-1)
    t[i] = e[i]/(cv*d[i])
end

tmax=1e5*yr
tempo = 0. 
ncicli = 0
bonobo = 0

display("Initial condition")
display(v[2:4])
display(d[2:4])
display(t[2:4])
display(e[2:4]/e0)


#--------------- CICLO TEMPORALE PRIMCIPALE ----------------
while tempo <tmax

    #global bonobo
    #global ncicli
    #global tempo
    #global gam, cfl, deltax, yr

    ncicli += 1
    bonobo += 1
    if bonobo > 500
        orangotango = true
        bonobo = 0
    else 
        orangotango = false
    end

    #equazione dei gas perfetti
    for i=1:N
        p[i] = (gam-1)*e[i]
    end

    #definisco dtmin
    dtmin = 1e30
    for i=2:N-1        
        dtmin = min(dtmin, (xb[i]-xb[i-1])/(abs(v[i])+sqrt(gam*p[i]/d[i])))
    end

    dtmin = cfl*dtmin
    tempo += dtmin

    if ncicli == 2
        println("first time step = ", dtmin/yr)
        println("initial velocity = ", v[4]/1e5)
        #println(gam,yr,cfl)
    elseif ncicli > 1000
        display("cocco bello")
        break
    end
    end

    #------------ SOURCE STEP --------------
    for i=2:N-1
        v[i] = v[i]-dtmin*2*(p[i]-p[i-1])/((d[i]+d[i-1]))*deltax
        if v[i]<0
            v[i] = 0
        end
    end
    v[1]=v[2]
    v[N]=v[N-1]

    #artificial viscosity
    for i=2:N-1
        if v[i+1]-v[i]<0 
            q[i] = c2*d[i]*(v[i+1]-v[i])^2
        else
            q[i] = 0
        end
    end
    q[1] = q[2]
    q[N] = q[N-1]


        # SUBSTEP II: AGGIORNAMENTO DELLA VELOCITÃ FOR ARTIFICIAL VISC. [as additional term to pression] AND ENERGY
        for i=2: N-1
            v[i]=v[i]-dtmin*2*(q[i]-q[i-1])/((d[i]+d[i-1])*deltax)
            if v[i]<1e-20
                v[i]=0.
            end
        end
        v[1] = v[2]
        v[N] = v[N-1]
       

        # SUBSTEP II: ENERGY 
        for i=2: N-1
            e[i]=e[i]-dtmin*q[i]*(v[i+1]-v[i])/deltax
        end
        e[1] = e[2]
        e[N] = e[N-1]

        #SUBSTEP III: COMPRESSION HEATING
        for i=2:N-1
            divV[i]=(g2a[i+1]*g31a[i+1]*v[i+1]-g2a[i]*g31a[i]*v[i])/dvl1a[i]
        end
        divV[1] = divV[2]
        divV[N] =divV[N-1]
      

        for i=2: N-1
            e[i]=e[i]*(1-0.5*dtmin*(gam-1.)*divV[i])/(1+0.5*dtmin*(gam-1.)*divV[i])
        end
        e[1] = e[2]
        e[N] = e[N-1]

        println("se stai leggendo questo sei uhn bomber")

end #fine della funzione
main()


