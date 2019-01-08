dt=.0001;
tvec = 0:dt:4;
preSynFiring = zeros(size(tvec));
for i = 1:length(tvec)
    if i >= 0 && i <= 1/dt
       preSynFiring(1,i) = 20;
    end
    if i > 1/dt && i <= 2/dt
       preSynFiring(1,i) = 100;
    end
    if i > 2/dt && i <= 3/dt
       preSynFiring(1,i) = 10;
    end
    if i > 3/dt && i <= 4/dt
       preSynFiring(1,i) = 50;  
    end
end

spikes = rand(size(preSynFiring)) < preSynFiring * dt;

syn_g1 = zeros(size(tvec));

syn_tau = .100;
delta_g = 1e-9;

for i = 2:length(tvec)
    if spikes(i-1) == 1
        syn_g1(i) = syn_g1(i-1) + delta_g;
    else
        dsyngdt = -syn_g1(i-1)/syn_tau;
        syn_g1(i) = syn_g1(i-1) + dt*dsyngdt;
    end 
end

figure
subplot(2, 1, 1)
plot(tvec, preSynFiring)
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate over Time')
subplot(2, 1, 2)
plot(tvec, syn_g1)
xlabel('Time')
ylabel('Conductance (S)')
title('Conductance over Time')

%%
dt=.0001;
tvec = 0:dt:4;
preSynFiring = zeros(size(tvec));
for i = 1:length(tvec)
    if i >= 0 && i <= 1/dt
       preSynFiring(1,i) = 20;
    end
    if i > 1/dt && i <= 2/dt
       preSynFiring(1,i) = 100;
    end
    if i > 2/dt && i <= 3/dt
       preSynFiring(1,i) = 10;
    end
    if i > 3/dt && i <= 4/dt
       preSynFiring(1,i) = 50;  
    end
end

spikes = rand(size(preSynFiring)) < preSynFiring * dt;

syn_g1 = zeros(size(tvec));

syn_tau = .100;
delta_g = 1e-9; 

depression1 = zeros(size(tvec));
depression1(1) = 1;
p0 = .5;
d_tau = .25;

syn_g2 = zeros(size(tvec));
gmax = 4e-9;

for i = 2:length(tvec)
    if spikes(i-1) == 1
        depression1(i) = depression1(i-1) - p0*depression1(i-1);
        
        syn_g2(i) = syn_g2(i-1) + gmax*p0*depression1(i-1);
    else
        dDdt = (1-depression1(i-1))/d_tau;
        depression1(i) = depression1(i-1) + dt*dDdt;
        
        dsyngdt = -syn_g2(i-1)/syn_tau;
        syn_g2(i) = syn_g2(i-1) + dt*dsyngdt;
    end 
end

figure
subplot(2, 1, 1)
plot(tvec, preSynFiring)
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate over Time')
subplot(2, 1, 2)
plot(tvec, syn_g2)
xlabel('Time')
ylabel('Conductance (S)')
title('Conductance over Time')

%%

dt=.0001;
tvec = 0:dt:4;
preSynFiring = zeros(size(tvec));
for i = 1:length(tvec)
    if i >= 0 && i <= 1/dt
       preSynFiring(1,i) = 20;
    end
    if i > 1/dt && i <= 2/dt
       preSynFiring(1,i) = 100;
    end
    if i > 2/dt && i <= 3/dt
       preSynFiring(1,i) = 10;
    end
    if i > 3/dt && i <= 4/dt
       preSynFiring(1,i) = 50;  
    end
end

spikes = rand(size(preSynFiring)) < preSynFiring * dt;

syn_tau = .100; 

p0 = .2;
d_tau = .25;

gmax = 4e-9;

facilitation = ones(size(tvec));
f_tau = .25;
f_fac = .25;
f_max = 1/p0;

depression2 = ones(size(tvec));

syn_g3 = zeros(size(tvec));

for i = 2:length(tvec)
    if spikes(i-1) == 1
        facilitation(i) = facilitation(i-1) + f_fac*(f_max-facilitation(i-1));
        
        depression2(i) = depression2(i-1) - p0*facilitation(i-1)*depression2(i-1);
        
        syn_g3(i) = syn_g3(i-1) + gmax*p0*facilitation(i-1)*depression2(i-1);
    else
        dFdt = (1-facilitation(i-1))/f_tau;
        facilitation(i) = facilitation(i-1) + dt*dFdt;
        
        dDdt = (1-depression2(i-1))/d_tau;
        depression2(i) = depression2(i-1) + dt*dDdt;
        
        dsyngdt = -syn_g3(i-1)/syn_tau;
        syn_g3(i) = syn_g3(i-1) + dt*dsyngdt;
    end 
end

figure
subplot(3, 1, 1)
plot(tvec, preSynFiring)
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate over Time')
subplot(3, 1, 2)
plot(tvec, syn_g3)
xlabel('Time')
ylabel('Conductance (S)')
title('Conductance over Time')
subplot(3, 1, 3)
plot(tvec, depression2)

%%

dt=.0001;
tvec = 0:dt:4;
preSynFiring = zeros(size(tvec));
for i = 1:length(tvec)
    if i >= 0 && i <= 1/dt
       preSynFiring(1,i) = 20;
    end
    if i > 1/dt && i <= 2/dt
       preSynFiring(1,i) = 100;
    end
    if i > 2/dt && i <= 3/dt
       preSynFiring(1,i) = 10;
    end
    if i > 3/dt && i <= 4/dt
       preSynFiring(1,i) = 50;  
    end
end

syn_g1 = zeros(size(tvec));

syn_tau = .002;
delta_g = 1e-9; 

depression1 = ones(size(tvec));
p0 = .5;
d_tau = .25;

syn_g2 = zeros(size(tvec));
gmax = 2e-9;

g1total = zeros(size(tvec));
g2total = zeros(size(tvec));

for j = 1:50
    spikes = rand(size(preSynFiring)) < preSynFiring * dt;

    for i = 2:length(tvec)
        if spikes(i-1) == 1
            syn_g1(i) = syn_g1(i-1) + delta_g;

            depression1(i) = depression1(i-1) - p0*depression1(i-1);

            syn_g2(i) = syn_g2(i-1) + gmax*p0*depression1(i-1);
        else
            dsyngdt = -syn_g1(i-1)/syn_tau;
            syn_g1(i) = syn_g1(i-1) + dt*dsyngdt;

            dDdt = (1-depression1(i-1))/d_tau;
            depression1(i) = depression1(i-1) + dt*dDdt;

            dsyngdt = -syn_g2(i-1)/syn_tau;
            syn_g2(i) = syn_g2(i-1) + dt*dsyngdt;
        end 
        g1total(i-1) = g1total(i-1) + syn_g1(i-1);
        g2total(i-1) = g2total(i-1) + syn_g2(i-1);
    end
end   

figure
subplot(2, 1, 1)
plot(tvec, preSynFiring)
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate over Time')
subplot(2, 1, 2)
plot(tvec, g1total)
xlabel('Time')
ylabel('Conductance (S)')
title('Conductance over Time')

figure
subplot(2, 1, 1)
plot(tvec, preSynFiring)
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate over Time')
subplot(2, 1, 2)
plot(tvec, g2total)
xlabel('Time')
ylabel('Conductance (S)')
title('Conductance over Time')

%%

dt=.0001;
tvec = 0:dt:4;
preSynFiring = zeros(size(tvec));
for i = 1:length(tvec)
    if i >= 0 && i <= 1/dt
       preSynFiring(1,i) = 20;
    end
    if i > 1/dt && i <= 2/dt
       preSynFiring(1,i) = 100;
    end
    if i > 2/dt && i <= 3/dt
       preSynFiring(1,i) = 10;
    end
    if i > 3/dt && i <= 4/dt
       preSynFiring(1,i) = 50;  
    end
end

spikes = rand(size(preSynFiring)) < preSynFiring * dt;

syn_g1 = zeros(size(tvec));

syn_tau = .100;
delta_g = 1e-9; 

depression1 = ones(size(tvec));
p0 = .2;
d_tau = .25;

syn_g2 = zeros(size(tvec));
gmax = 4e-9;

facilitation = ones(size(tvec));
f_tau = .25;
f_fac = .25;
f_max = 1/p0;

depression2 = ones(size(tvec));

syn_g3 = zeros(size(tvec));

Vm = zeros(size(tvec));
G_leak = 2e-9;
E_leak = -.065;
E_Syn = 0;
THRESHOLD = -.050;
RESET = -.080;
Cm = 20e-12;
Vm(1) = E_leak;

trial = 'B';

for i = 2:length(tvec)
    if spikes(i-1) == 1
        syn_g1(i) = syn_g1(i-1) + delta_g;
        
        depression1(i) = depression1(i-1) - p0*depression1(i-1);
        
        syn_g2(i) = syn_g2(i-1) + gmax*p0*depression1(i-1);
        
        facilitation(i) = facilitation(i-1) + f_fac*(f_max-facilitation(i-1));
        
        depression2(i) = depression2(i-1) - p0*facilitation(i-1)*depression2(i-1);
        
        syn_g3(i) = syn_g3(i-1) + gmax*p0*facilitation(i-1)*depression2(i-1);
    else
        dsyngdt = -syn_g1(i-1)/syn_tau;
        syn_g1(i) = syn_g1(i-1) + dt*dsyngdt;
        
        dDdt = (1-depression1(i-1))/d_tau;
        depression1(i) = depression1(i-1) + dt*dDdt;
        
        dsyngdt = -syn_g2(i-1)/syn_tau;
        syn_g2(i) = syn_g2(i-1) + dt*dsyngdt;
        
        dFdt = (1-facilitation(i-1))/f_tau;
        facilitation(i) = facilitation(i-1) + dt*dFdt;
        
        dDdt = (1-depression2(i-1))/d_tau;
        depression2(i) = depression2(i-1) + dt*dDdt;
        
        dsyngdt = -syn_g3(i-1)/syn_tau;
        syn_g3(i) = syn_g3(i-1) + dt*dsyngdt;
    end 
    
    switch trial
        case 'A'
            syn_g = syn_g1(i-1);
        case 'B'
            syn_g = syn_g2(i-1);
        case 'C'
            syn_g = syn_g3(i-1);
    end
    
    dVdt = (G_leak*(E_leak-Vm(i-1)) + syn_g*(E_Syn-Vm(i-1)))/Cm; 
    Vm(i) = Vm(i-1) + dt*dVdt;
    if Vm(i) > THRESHOLD
           Vm(i) = RESET;
    end
end

figure
subplot(2, 1, 1)
plot(tvec, preSynFiring)
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate over Time')
subplot(2, 1, 2)
plot(tvec, Vm)
xlabel('Time')
ylabel('Potential (V)')
title('Potential over Time')
