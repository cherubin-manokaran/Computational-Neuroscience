Vm = -.085:.001:.050;
step_Ca = 2e-3/(length(Vm)-1);
Ca = 0:step_Ca:2e-3;

[ alpha_mca, beta_mca, alpha_kca, beta_kca, alpha_kahp, beta_kahp ] =...
    PR_dend_gating( Vm, Ca );
figure
plot(Vm, alpha_mca)
hold on
plot(Vm, beta_mca)
xlabel('Potential (V)')
ylabel('Gating Variable Values')
title('Gating Variable Values against Potential (mca)')
legend('alpha', 'beta') 

figure
plot(Vm, alpha_kca)
hold on
plot(Vm, beta_kca)
xlabel('Potential (V)')
ylabel('Gating Variable Values')
title('Gating Variable Values against Potential (kca)')
legend('alpha', 'beta') 

figure
plot(Vm, alpha_kahp)
hold on
plot(Vm, beta_kahp)
xlabel('Potential (V)')
ylabel('Gating Variable Values')
title('Gating Variable Values against Potential (kahp)')
legend('alpha', 'beta') 

plot(alpha_kahp, beta_kahp)
[ alpha_m, beta_m, alpha_h, beta_h , alpha_n, beta_n] =...
    PR_soma_gating( Vm );

figure
plot(Vm, alpha_m)
hold on
plot(Vm, beta_m)
xlabel('Potential (V)')
ylabel('Gating Variable Values')
title('Gating Variable Values against Potential (m)')
legend('alpha', 'beta') 

figure
plot(Vm, alpha_h)
hold on
plot(Vm, beta_h)
xlabel('Potential (V)')
ylabel('Gating Variable Values')
title('Gating Variable Values against Potential (h)')
legend('alpha', 'beta') 

figure
plot(Vm, alpha_n)
hold on
plot(Vm, beta_n)
xlabel('Potential (V)')
ylabel('Gating Variable Values')
title('Gating Variable Values against Potential (n)')
legend('alpha', 'beta') 
%%

A_S = 1/3;
A_D = 1-A_S;
G_SL = A_S*5e-9;
G_DL = A_D*5e-9;
G_Na = A_S*3e-6;
G_K = A_S*2e-6;
G_Ca = A_D*2e-6;
G_KCa = A_D*2.5e-6;
G_KAHP = A_D*40e-9;
G_Link = 20e-9;
E_Na = .060;
E_Ca = .080;
E_K = -.075;
E_L = -.060;
C_S = A_S*100e-12;
C_D = A_D*100e-12;
I_Sapp = 0;
I_Dapp = 0;
tau_ca = 50e-3;
k = 5e6/A_D;

dt = .000002;
tvec = 0:dt:2;

part = 'D';
switch (part)
    case 'A'
        G_Link = 20e-9;
    case 'B'
        G_Link = 0;
    case 'C'
        G_Link = 10e-9;
    case 'D'
        G_Link = 100e-8;
    case 'E'
        I_Dapp = 50e-12;
    case 'F'
        I_Dapp = 100e-12;
    case 'G'
        I_Dapp = 200e-12;
    case 'H'
        I_Sapp = 50e-12;
    case 'I'
        I_Sapp = 100e-12;
    case 'J'
        I_Sapp = 200e-12;
end

cavec = zeros(size(tvec));
icavec = zeros(size(tvec));

vdvec = zeros(size(tvec));
vdvec(1) = -.060;

vsvec = zeros(size(tvec));
vsvec(1) = -.060;

mvec = zeros(size(tvec));

hvec = zeros(size(tvec));

nvec = zeros(size(tvec));

mcavec = zeros(size(tvec));

kcavec = zeros(size(tvec));

mkahpvec = zeros(size(tvec));

for i = 2:length(tvec)
    Ca = cavec(i-1);
    I_Ca = icavec(i-1);
    Vm_D = vdvec(i-1);
    Vm_S = vsvec(i-1);
    m = mvec(i-1);
    h = hvec(i-1);
    n = nvec(i-1);
    mca = mcavec(i-1);
    kca = kcavec(i-1);
    mkahp = mkahpvec(i-1);
    
    [ alpha_mca, beta_mca, alpha_kca, beta_kca, alpha_kahp, beta_kahp ] = ...
    PR_dend_gating( Vm_D, Ca );
    
    [ alpha_m, beta_m, alpha_h, beta_h , alpha_n, beta_n] =...
    PR_soma_gating( Vm_S );

    m_inf = alpha_m./(alpha_m+beta_m);      % sodium activation variable

    tau_h = 1./(alpha_h+beta_h);            % sodium inactivation time constant
    h_inf = alpha_h./(alpha_h+beta_h);      % sodium inactivation variable

    tau_n = 1./(alpha_n+beta_n);            % potassium activation time constant
    n_inf = alpha_n./(alpha_n+beta_n);      % potassium activation variable
    
    tau_mca = 1./(alpha_mca+beta_mca);           
    mca_inf = alpha_mca./(alpha_mca+beta_mca);
    
    tau_kca = 1./(alpha_kca+beta_kca);           
    kca_inf = alpha_kca./(alpha_kca+beta_kca); 
    
    tau_mkahp = 1./(alpha_kahp+beta_kahp);
    mkahp_inf = alpha_kahp./(alpha_kahp+beta_kahp);
    
    dCadt = -Ca/tau_ca + k*I_Ca;
    cavec(i) = Ca + dt * dCadt;
    I_Ca = G_Ca*mca*mca*(E_Ca-Vm_D);
    icavec(i) = I_Ca;
    
    X = min(1,4000*Ca);
    
    dVsmdt = (G_SL*(E_L-Vm_S) + G_Na*(m*m)*h*(E_Na-Vm_S) + ...
        G_K*(n*n)*(E_K-Vm_S) + G_Link*(Vm_D-Vm_S) + I_Sapp)/C_S;
    vsvec(i) = Vm_S + dt * dVsmdt;
    
    dVdmdt = (G_DL*(E_L-Vm_D) + G_Ca*(mca*mca)*(E_Ca-Vm_D) +...
        G_KCa*kca*X*(E_K-Vm_D)+ G_KAHP*mkahp*(E_K-Vm_D) -...
        G_Link*(Vm_D-Vm_S) + I_Dapp)/C_D;
    vdvec(i) = Vm_D + dt * dVdmdt;
    
    dmcadt = alpha_mca*(1-mca) - beta_mca*mca;
    mcavec(i) = mca + dt * dmcadt;
    
    dkcadt = alpha_kca*(1-kca) - beta_kca*kca;
    kcavec(i) = kca + dt * dkcadt;
    
    dmkahpdt = alpha_kahp*(1-mkahp) - beta_kahp*mkahp;
    mkahpvec(i) = mkahp + dt * dmkahpdt;
    
    dmdt = alpha_m*(1-m) - beta_m*m;
    mvec(i) = m + dt * dmdt;
    
    dhdt = alpha_h*(1-h) - beta_h*h;
    hvec(i) = h + dt * dhdt;
    
    dndt = alpha_n*(1-n) - beta_n*n;
    nvec(i) = n + dt * dndt;
end

% Now plot the functions after setting up graphics parameters
% set(0,'DefaultLineLineWidth',2,...
%     'DefaultLineMarkerSize',8, ...
%     'DefaultAxesLineWidth',2, ...
%     'DefaultAxesFontSize',14,...
%     'DefaultAxesFontWeight','Bold');

figure
plot(tvec, vsvec)
xlabel('Time')
ylabel('Potential')
title('Somatic Potential over Time')


figure
subplot(2, 1, 1)
plot(tvec, vdvec)
xlabel('Time')
ylabel('Potential')
title('Dendritic Potential over Time')
subplot(2, 1, 2)
plot(tvec, cavec)
xlabel('Time')
ylabel('Calcium Current (amperes)')
title('Calcium Current over Time')



