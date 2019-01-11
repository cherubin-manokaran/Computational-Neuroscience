G_leak = 30e-9;
G_na = 12e-6;
G_k = 4e-6;
E_na = .05;
E_k = -.075;
E_leak = -.06;
C_M = 100e-12;

dt = .00001;
tvec = 0:dt:.3; 

vvec = zeros(size(tvec));
vvec(1) = E_leak;

mvec = zeros(size(tvec));
mvec(1) = 0;

hvec = zeros(size(tvec));
hvec(1) = 0;

nvec = zeros(size(tvec));
nvec(1) = 0;

Iapp = 0;

for i = 2:length(tvec)
    
    Vm = vvec(i-1);          % membane potential for calculations
    m = mvec(i-1);
    h = hvec(i-1);
    n = nvec(i-1);
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else                    % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.075)/0.08);     % potassium deactivation rate
    
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    dVmdt = (G_leak*(E_leak-Vm) + G_na*(m*m*m)*h*(E_na-Vm) + G_k*(n*n*n*n)*(E_k-Vm) + Iapp)/C_M;
    vvec(i) = Vm + dt * dVmdt;
    
    dmdt = alpha_m*(1-m) - beta_m*m;
    mvec(i) = m + dt * dmdt;
    
    dhdt = alpha_h*(1-h) - beta_h*h;
    hvec(i) = h + dt * dhdt;
    
    dndt = alpha_n*(1-n) - beta_n*n;
    nvec(i) = n + dt * dndt;
end

figure;

plot(tvec, vvec)
xlabel('Time')
ylabel('Potential')

%

G_leak = 30e-9;
G_na = 12e-6;
G_k = 4e-6;
E_na = .05;
E_k = -.075;
E_leak = -.06;
C_M = 100e-12;

dt = .00001;
tvec = 0:dt:.3; 

vvec = zeros(size(tvec));
vvec(1) = E_leak;

mvec = zeros(size(tvec));
mvec(1) = 0;

hvec = zeros(size(tvec));
hvec(1) = 0;

nvec = zeros(size(tvec));
nvec(1) = 0;

cvec = zeros(size(tvec));

for i = 2:length(tvec)
    
    Vm = vvec(i-1);          % membane potential for calculations
    m = mvec(i-1);
    h = hvec(i-1);
    n = nvec(i-1);
    
    time = tvec(i-1);
    if (time >= .05 && time <= .150)
        cvec(i-1) = .18e-9;
    end
    
    Iapp = cvec(i-1);
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else                    % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.075)/0.08);     % potassium deactivation rate
    
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    dVmdt = (G_leak*(E_leak-Vm) + G_na*(m*m*m)*h*(E_na-Vm) + G_k*(n*n*n*n)*(E_k-Vm) + Iapp)/C_M;
    vvec(i) = Vm + dt * dVmdt;
    
    dmdt = alpha_m*(1-m) - beta_m*m;
    mvec(i) = m + dt * dmdt;
    
    dhdt = alpha_h*(1-h) - beta_h*h;
    hvec(i) = h + dt * dhdt;
    
    dndt = alpha_n*(1-n) - beta_n*n;
    nvec(i) = n + dt * dndt;
end

figure;

subplot(2, 1, 1)
plot(tvec, cvec)
xlabel('Time')
ylabel('Current')

subplot(2, 1, 2)
plot(tvec, vvec)
xlabel('Time')
ylabel('Potential')

%%

G_leak = 30e-9;
G_na = 12e-6;
G_k = 4e-6;
E_na = .05;
E_k = -.075;
E_leak = -.06;
C_M = 100e-12;

dt = .00001;
tvec = 0:dt:.3; 
steplength = .005;
nsteplength = round(steplength/dt);

vvec = zeros(size(tvec));
vvec(1) = E_leak;

mvec = zeros(size(tvec));
mvec(1) = 0;

hvec = zeros(size(tvec));
hvec(1) = 0;

nvec = zeros(size(tvec));
nvec(1) = 0;

cvec = zeros(size(tvec));

for step = 1:10
    istart = (step+9)*nsteplength+(step-1)*1000;
    istop = (step+10)*nsteplength+(step-1)*1000;
    cvec(istart:istop) = .18e-9;
end

for i = 2:length(tvec)
    
    Vm = vvec(i-1);          % membane potential for calculations
    m = mvec(i-1);
    h = hvec(i-1);
    n = nvec(i-1);
    
    time = tvec(i-1);
    
    Iapp = cvec(i-1);
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else                    % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.075)/0.08);     % potassium deactivation rate
    
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    dVmdt = (G_leak*(E_leak-Vm) + G_na*(m*m*m)*h*(E_na-Vm) + G_k*(n*n*n*n)*(E_k-Vm) + Iapp)/C_M;
    vvec(i) = Vm + dt * dVmdt;
    
    dmdt = alpha_m*(1-m) - beta_m*m;
    mvec(i) = m + dt * dmdt;
    
    dhdt = alpha_h*(1-h) - beta_h*h;
    hvec(i) = h + dt * dhdt;
    
    dndt = alpha_n*(1-n) - beta_n*n;
    nvec(i) = n + dt * dndt;
end

figure;

subplot(2, 1, 1)
plot(tvec, cvec)
xlabel('Time')
ylabel('Current')

subplot(2, 1, 2)
plot(tvec, vvec)
xlabel('Time')
ylabel('Potential')

%%

G_leak = 30e-9;
G_na = 12e-6;
G_k = 4e-6;
E_na = .05;
E_k = -.075;
E_leak = -.06;
C_M = 100e-12;

dt = .00001;
tvec = 0:dt:.3; 
steplength = .005;
nsteplength = round(steplength/dt);

vvec = zeros(size(tvec));
vvec(1) = -.065;

mvec = zeros(size(tvec));
mvec(1) = .05;

hvec = zeros(size(tvec));
hvec(1) = .5;

nvec = zeros(size(tvec));
nvec(1) = .35;

cvec = zeros(size(tvec)) + .65e-9;

for step = 1:10
    istart = (step)*nsteplength+(step-1)*2000;
    istop = (step+1)*nsteplength+(step-1)*2000;
    cvec(istart:istop) = 0;
end

for i = 2:length(tvec)
    
    Vm = vvec(i-1);          % membane potential for calculations
    m = mvec(i-1);
    h = hvec(i-1);
    n = nvec(i-1);
    
    time = tvec(i-1);
    
    Iapp = cvec(i-1);
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else                    % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.075)/0.08);     % potassium deactivation rate
    
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    dVmdt = (G_leak*(E_leak-Vm) + G_na*(m*m*m)*h*(E_na-Vm) + G_k*(n*n*n*n)*(E_k-Vm) + Iapp)/C_M;
    vvec(i) = Vm + dt * dVmdt;
    
    dmdt = alpha_m*(1-m) - beta_m*m;
    mvec(i) = m + dt * dmdt;
    
    dhdt = alpha_h*(1-h) - beta_h*h;
    hvec(i) = h + dt * dhdt;
    
    dndt = alpha_n*(1-n) - beta_n*n;
    nvec(i) = n + dt * dndt;
end

figure;

subplot(2, 1, 1)
plot(tvec, cvec)
xlabel('Time')
ylabel('Current')

subplot(2, 1, 2)
plot(tvec, vvec)
xlabel('Time')
ylabel('Potential')

%%

G_leak = 30e-9;
G_na = 12e-6;
G_k = 4e-6;
E_na = .05;
E_k = -.075;
E_leak = -.06;
C_M = 100e-12;

dt = .00001;
tvec = 0:dt:.3; 
steplength = .005;
nsteplength = round(steplength/dt);

vvec = zeros(size(tvec));
vvec(1) = -.065;

mvec = zeros(size(tvec));
mvec(1) = .05;

hvec = zeros(size(tvec));
hvec(1) = .5;

nvec = zeros(size(tvec));
nvec(1) = .35;

cvec = zeros(size(tvec)) + .7e-9;

for i = 2:length(tvec)
    
    Vm = vvec(i-1);          % membane potential for calculations
    m = mvec(i-1);
    h = hvec(i-1);
    n = nvec(i-1);
    
    time = tvec(i-1);
    if (time >= .05 && time <= .055)
        cvec(i-1) = 1e-9;
    end
    
    Iapp = cvec(i-1);
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else                    % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.075)/0.08);     % potassium deactivation rate
    
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    dVmdt = (G_leak*(E_leak-Vm) + G_na*(m*m*m)*h*(E_na-Vm) + G_k*(n*n*n*n)*(E_k-Vm) + Iapp)/C_M;
    vvec(i) = Vm + dt * dVmdt;
    
    dmdt = alpha_m*(1-m) - beta_m*m;
    mvec(i) = m + dt * dmdt;
    
    dhdt = alpha_h*(1-h) - beta_h*h;
    hvec(i) = h + dt * dhdt;
    
    dndt = alpha_n*(1-n) - beta_n*n;
    nvec(i) = n + dt * dndt;
end

figure;

subplot(2, 1, 1)
plot(tvec, cvec)
xlabel('Time')
ylabel('Current')

subplot(2, 1, 2)
plot(tvec, vvec)
xlabel('Time')
ylabel('Potential')

%%

G_leak = 30e-9;
G_na = 12e-6;
G_k = 4e-6;
E_na = .05;
E_k = -.075;
E_leak = -.06;
C_M = 100e-12;

dt = .00001;
tvec = 0:dt:.3; 
steplength = .005;
nsteplength = round(steplength/dt);

vvec = zeros(size(tvec));
vvec(1) = -.065;

mvec = zeros(size(tvec));
mvec(1) = 0;

hvec = zeros(size(tvec));
hvec(1) = 0;

nvec = zeros(size(tvec));
nvec(1) = 0;

cvec = zeros(size(tvec)) + .7e-9;

for i = 2:length(tvec)
    
    Vm = vvec(i-1);          % membane potential for calculations
    m = mvec(i-1);
    h = hvec(i-1);
    n = nvec(i-1);
    
    time = tvec(i-1);
    if (time >= .05 && time <= .055)
        cvec(i-1) = 1e-9;
    end
    
    Iapp = cvec(i-1);
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else                    % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.075)/0.08);     % potassium deactivation rate
    
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    dVmdt = (G_leak*(E_leak-Vm) + G_na*(m*m*m)*h*(E_na-Vm) + G_k*(n*n*n*n)*(E_k-Vm) + Iapp)/C_M;
    vvec(i) = Vm + dt * dVmdt;
    
    dmdt = alpha_m*(1-m) - beta_m*m;
    mvec(i) = m + dt * dmdt;
    
    dhdt = alpha_h*(1-h) - beta_h*h;
    hvec(i) = h + dt * dhdt;
    
    dndt = alpha_n*(1-n) - beta_n*n;
    nvec(i) = n + dt * dndt;
end

figure;

subplot(2, 1, 1)
plot(tvec, cvec)
xlabel('Time')
ylabel('Current')

subplot(2, 1, 2)
plot(tvec, vvec)
xlabel('Time')
ylabel('Potential')
