using GLMakie



function main()

    function Cool(Temp1)
        
        Temp_kev = Temp1/1.16e7
        
        if Temp_kev > 0.02
            return (8.6e-3*(Temp_kev^(-1.7)) + 0.058*(Temp_kev^(0.5)) + 0.063)*1e-22
        elseif Temp_kev <= 0.02 && Temp_kev >= 0.0017235
            return 6.72e-22*((Temp_kev/0.02)^0.6)
        elseif Temp_kev < 0.0017235
            return 1.544e-22*((Temp_kev/0.0017235)^6)
        end
    end
    
#definizione de parametri
N = 1000
cmpc = 3.085e18
mu = 1.4
msol = 1.989e33
mp =  1.67e-24
cv = 2e8
gam = 1.667
yr = 3.156e7




xa = zeros(N)
xb = copy(xa)
dxa = copy(xa)
dxb = copy(xa)


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

for i in 2:N-1
    dxa[i] = xa[i+1] - xa[i] # delta spaziale
    dxb[i] = xb[i] - xb[i-1]
end
dxa[1] = dxa[2]
dxa[N] = dxa[N-1]
dxb[1] = dxb[2]
dxb[N] = dxb[N-1]







#fattori di scala
g2a = zeros(N)
g2b = zeros(N)
g31a = zeros(N)
g31b = zeros(N)
dvl1a = zeros(N)
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
C2 = 3
cfl = 0.01

d0 = 2e-24
t0 = 1e4
t0_min = 1e4
#n0=d0/(mu*mp)
vol = 1.333*pi*xa[4]^3

v = zeros(N)
d = copy(v)
e = copy(v)
t = copy(v)
p = copy(v)
q = copy(v)
divV = copy(v)
F3 = copy(v)

q = zeros(N)
s = copy(v)
dstar = copy(v)
F1 = copy(v)
M = copy(v)
e_dstar = copy(v)
vstar = copy(v)
F2 = copy(v)

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
winds = false

display("Initial condition")
display(xb[1])
display(v[2:4])
display(d[2:4])
display(t[2:4])
display(e[2:4]/e0)


fig = Figure()
ax = Axis(fig[1,1],
    title= "density",
    xlabel=" distance",
    ylabel = "density",
    #xscale = log10

    )
    ax1 = Axis(fig[2,1],
    title="velocity",
    xlabel="distance",
    ylabel="velocity")

ax2 = Axis(fig[3,1],
title="temperature",
xlabel="distance",
ylabel = "T")

plot = lines!(ax, xb/cmpc,d)
xlims!(ax, 0, 75)
#ylims!(ax, 1e-28 ,1e-23)
plot2 = lines!(ax1, xb/cmpc, v/1e5, color=:red)


plot3 = lines!(ax2, xb/cmpc, t, color= :green)


display(fig)
println("Hit enter when plot is ready")
readline()


v_winds = 1e8
m_lost = 1e19

if winds == true
    println("STELLAR WINDS PHASE")
    # call winds_injection(v, d, e, p, t, xb, yr/100) !injection for 1 month
    # non modifico le condizioni iniziali eccetto la velocità
    for i = 1:3
        v[i] = v_winds
        t[i] = t0
        d[i] = d0
        e[i] = cv * d[i] * t[i]
    end

    tmax = 1.0e6 * yr
end


dtmin = 0

#--------------- CICLO TEMPORALE PRIMCIPALE ----------------
while tempo <tmax

    if winds == true
        for i = 1:3
            v[i] = v_winds
            t[i] = t0
            d[i] += m_lost * dtmin / vol # dens = dens precedente + materia aggiunta
            e[i] += 0.5 * m_lost * dtmin * v_winds^2 / vol # energia = energia meccanica dei venti
        end
    end
    

    ncicli += 1
    bonobo += 1
    if bonobo > 10
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
        dtmin=min(dtmin,(xb[i]-xb[i-1])/(abs(v[i])+sqrt(gam*p[i]/d[i])))
    end
    
    dtmin=cfl*dtmin #cfl is the C constant 
    
    tempo+=dtmin
    
    if ncicli == 2
        println("FIRST TIME-STEP == ", dtmin/yr)
        println("initial energy ==", e[2]*vol)
        if winds == false
            println("initial velocity ==", v[4]/1e5)
            #println("expected temperature")
            #write(*,1034)  3*mu*mp/(16*kbol)*v(4)**2
        else
            println("initial velocity ==", v[3]/1e5)
        end
    end
    

    #------------ SOURCE STEP --------------
    #1 AGGIORNAMENTO DELLE VELOITA
    for i=2:N-1
        v[i] = v[i] - dtmin*2*(p[i]-p[i-1])/((d[i]+d[i-1])*dxb[i])
        if v[i] < 1e-20
            v[i] = 0
        end
        
    end
    v[1] = v[2]
    v[N] = v[N-1]
    
    #2 ARTIFICIAL VISCOSITY
    for i=2:N-1
        if v[i+1] - v[i] < 0
            q[i] = C2*d[i]*(v[i+1] - v[i])^2
        else
            q[i] = 0
        end
    end
    q[1] = q[2]
    q[N] = q[N-1]
    
    #AGG VEL CON VISCOSITA
    for i=2:N-1
        v[i] = v[i] - dtmin*2*(q[i]-q[i-1])/((d[i]+d[i-1])*dxb[i])
        if v[i] < 1e-20
            v[i] = 0
        end
    end
    v[1]=v[2]
    v[N] = v[N-1]
    
    #ENERGIA
    for i=2:N-1
        e[i] = e[i] - dtmin*q[i]*(v[i+1]-v[i])/dxa[i]
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

        
                # add the cooling function
        for i in 2:N-1
            e[i] -= dtmin*(d[i]/2.17e-24)^2 * Cool(t[i]) # energy with en. loss
        end
        e[1] = e[2]
        e[N] = e[N-1]



        
        for i=2: N-1
            t[i] = e[i]/(cv*d[i])
            if t[i] < t0_min
                t[i] = t0_min
            end
        end

        for i=2: N-1
            e[i]=cv*d[i]*t[i]
        end
        e[1] = e[2]
        e[N] = e[N-1]

                # Define the momentum density
        for i = 2:N-1
            s[i] = 0.5*(d[i] + d[i-1])*v[i]   # at "i"
        end
        s[1] = s[2]
        s[N] = s[N-1]


        # Update the density based on the velocity (upwind)
        for i = 2:N-1
            if v[i] > 0
                dstar[i] = d[i-1]   # at "i"
            else
                dstar[i] = d[i]
            end
        end
        dstar[N] = dstar[N-1]
        dstar[1] = dstar[3]

        for i = 2:N
            F1[i] = dstar[i]*v[i]*g2a[i]*g31a[i]   # density flux across surfaces at "i"
        end

        # Update the energy: calculate the surface momentum and energy on cell edges
        for i = 2:N-1
            M[i] = dstar[i]*v[i]   # surface momentum
        end
        M[1] = M[2]
        M[N] = M[N-1]

        for i = 2:N-1
            if v[i] > 0
                e_dstar[i] = e[i-1]/d[i-1]   # at "i"
            else
                e_dstar[i] = e[i]/d[i]
            end
        end
        e_dstar[N] = e_dstar[N-1]
        e_dstar[1] = e_dstar[3]

        # Calculate the density flux across surfaces

        # Update the density after the arrival of the matter flow
        for i = 2:N-1
            d[i] = d[i] - dtmin*(F1[i+1] - F1[i])/dvl1a[i]   
            if d[i] < 1e-40 
                d[i]= 0
            end # density after the arrival of matter flow
        end
        d[1] = d[2]
        d[N] = d[N-1]

        # Calculate the energy flux
        for i = 2:N
            F2[i] = e_dstar[i]*M[i]*g2a[i]*g31a[i]   # energy flux
        end
        F2[1] = F2[2]
        F2[N] = F2[N-1]

        # Update the energy after the arrival of the energy flow
        for i = 2:N-1
            e[i] = e[i] - dtmin*(F2[i+1] - F2[i])/dvl1a[i]   # energy after the arrival of energy flow
        end
        e[1] = e[2]
        e[N] = e[N-1]

        #-------------------

        for i=2:N-1
            if ((v[i-1]+v[i])*0.5 > 0)
                vstar[i] = v[i-1]  # at i-1/2
            else
                vstar[i] = v[i]
            end
        end
        vstar[1] = vstar[2]
        vstar[N] = vstar[N-1]
        
        for i=1:N-1
            F3[i] = vstar[i+1] * 0.5 * (M[i] + M[i+1]) * g2b[i] * g31b[i]  # at i+1/2
        end
        
        for i=2:N-1
            s[i] = s[i] - dtmin / dvl1b[i] * (F3[i] - F3[i-1])  # update momentum density
        end
        s[1] = s[2]
        s[N] = s[N-1]
        
        for i=2:N-1
            v[i] = 2. * s[i] / (d[i] + d[i-1])  # final velocity of the cells
        end
        v[1] = v[2]
        v[N] = v[N-1]
        
        if (cfl < 0.5)
            cfl = cfl + (cfl * 0.1)
        else
            cfl = 0.5
        end

        function text()
            println("hello world")
            
        end
        if orangotango == true
           # text()
            println("TIME PASSED == $(tempo/yr), $ncicli")
            
            delete!(plot.parent, plot)
            plot = lines!(ax, xb[2:end]/cmpc, d[2:end], color=:blue)
        
            xlims!(ax, 0, 75)

            delete!(plot2.parent, plot2)            
            plot2=lines!(ax1, xb/cmpc, v/1e5, color=:red)
        
           ylims!(ax, 1e-28, 1.6e-23)
           #ylims!(ax, 1e-28 ,1e-23)
           
           xlims!(ax1, 0, 75)
           ylims!(ax1, -10, 600)

           delete!(plot3.parent, plot3)
           plot3 = lines!(ax2, xb/cmpc, t, color= :green)
           xlims!(ax2, 0, 75)
           # sleep(0.01)
        end

        #delete!(makie_plot_3d_contour.parent, makie_plot_3d_contour)




        #println("se stai leggendo questo sei uhn bomber")

end 
end
main()


