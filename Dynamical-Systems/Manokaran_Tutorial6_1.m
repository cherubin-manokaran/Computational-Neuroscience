r_max = 100;
theta_e = -5;
theta_i = 0;
alpha_e = .05;
alpha_i = 1;
W_ee = 2;
W_ei = 2.5;
W_ie = -2.5;
W_ii = -2;

dt = 0.0001;
tvec = 0:dt:3;

I_e = zeros(size(tvec));
I_i = zeros(size(tvec));

r_evec = zeros(size(tvec));
r_ivec = zeros(size(tvec));

part = 'B';

switch part 
    case 'A'
        I_stim = 20;
        I_ebase = 0;
        I_ibase = 0;
        tau_e = 5e-3;
        tau_i = 5e-3;
    case 'B'
        I_stim = 20;
        I_ebase = 25;
        I_ibase = 15;
        tau_e = 5e-3;
        tau_i = 5e-3;
    case 'C'
       I_stim = 20;
       I_ebase = 0;
       I_ibase = 0;
       tau_e = 2e-3;
       tau_i = 10e-3; 
    case 'D'
       I_stim = 20;
       I_ebase = 25;
       I_ibase = 15;
       tau_e = 2e-3;
       tau_i = 10e-3; 
end

I_eapp = ones(size(tvec))*I_ebase;
I_iapp = ones(size(tvec))*I_ibase;

for t = 2:length(tvec)
    time = tvec(t);
    if (time >= 1 && time <= 2)
        I_iapp(t) = I_ibase + I_stim;
    end
    
    r_e = r_evec(t-1);
    r_i = r_ivec(t-1);
    
    I_e(t-1) = W_ee*r_e + W_ie*r_i + I_eapp(t-1);
    I_i(t-1) = W_ei*r_e + W_ii*r_i + I_iapp(t-1);
    
    if (I_e(t-1)-theta_e > 0)
        sign1 = 1;
    elseif (I_e(t-1)-theta_e < 0)
        sign1 = -1;
    else
        sign1 = 0;
    end
    dredt = (-r_e + alpha_e*(I_e(t-1) - theta_e) * sign1*(I_e(t-1) - theta_e))/tau_e;
    r_evec(t) = r_e + dredt*dt;
    
    dridt = (-r_i + alpha_i*(I_i(t-1) - theta_i))/tau_i;
    r_ivec(t) = r_i + dridt*dt;
    
    if (r_ivec(t) > r_max)
        r_ivec(t) = r_max; 
    end
    if (r_evec(t) > r_max)
        r_evec(t) = r_max;
    end
    
    if (r_ivec(t) < 0)
        r_ivec(t) = 0;
    end
    if (r_evec(t) < 0)
        r_evec(t) = 0;
    end
end

figure
subplot(2,1,1)
plot(tvec, I_iapp)
hold on
plot(tvec, I_eapp)
xlabel('Time')
ylabel('Current')
legend('Inhibitory Current', 'Excitatory Current')
subplot(2,1,2)
plot(tvec, r_ivec)
hold on
plot(tvec, r_evec)
xlabel('Time')
ylabel('Firing Rate')
legend('Inhibitory Firing Rate', 'Excitatory Firing Rate')

figure
plot(r_ivec, r_evec)
xlabel('Inhibitory Firing Rate')
ylabel('Excitatory Firing Rate')

%%

r_max = 100;
theta_e = -5;
theta_i = 0;
alpha_e = .05;
alpha_i = 1;
W_ee = 2;
W_ei = 2.5;
W_ie = -2.5;
W_ii = -2;
W_eix = 1.75;

dt = 0.0001;
tvec = 0:dt:3;

I1_e = zeros(size(tvec));
I1_i = zeros(size(tvec));
I2_e = zeros(size(tvec));
I2_i = zeros(size(tvec));

r1_evec = zeros(size(tvec));
r1_ivec = zeros(size(tvec));
r2_evec = zeros(size(tvec));
r2_ivec = zeros(size(tvec));

I_stim = 5;
I_ebase = 25;
I_ibase = 20;
tau_e = 5e-3;
tau_i = 5e-3; 

I1_eapp = ones(size(tvec))*I_ebase;
I1_iapp = ones(size(tvec))*I_ibase;

I2_eapp = ones(size(tvec))*I_ebase;
I2_iapp = ones(size(tvec))*I_ibase;

for t = 2:length(tvec)
    time = tvec(t);
    if (time >= 1 && time <= 1.1)
        I1_eapp(t) = I_ebase + I_stim;
    end
    
    if (time >= 2 && time <= 2.1)
        I2_eapp(t) = I_ebase + I_stim;
    end
    
    r1_e = r1_evec(t-1);
    r1_i = r1_ivec(t-1);
    r2_e = r2_evec(t-1);
    r2_i = r2_ivec(t-1);
    
    I1_e(t-1) = W_ee*r1_e + W_ie*r1_i + I1_eapp(t-1);
    I1_i(t-1) = W_eix*r2_e + W_ei*r1_e + W_ii*r1_i + I1_iapp(t-1);
    
    I2_e(t-1) = W_ee*r2_e + W_ie*r2_i + I2_eapp(t-1);
    I2_i(t-1) = W_eix*r1_e + W_ei*r2_e + W_ii*r2_i + I2_iapp(t-1);
    
    if (I1_e(t-1)-theta_e > 0)
        sign1 = 1;
    elseif (I1_e(t-1)-theta_e < 0)
        sign1 = -1;
    else
        sign1 = 0;
    end
    
    if (I2_e(t-1)-theta_e > 0)
        sign2 = 1;
    elseif (I2_e(t-1)-theta_e < 0)
        sign2 = -1;
    else
        sign2 = 0;
    end
    
    dr1edt = (-r1_e + alpha_e*(I1_e(t-1) - theta_e) * sign1*(I1_e(t-1) - theta_e))/tau_e;
    r1_evec(t) = r1_e + dr1edt*dt;
    
    dr1idt = (-r1_i + alpha_i*(I1_i(t-1) - theta_i))/tau_i;
    r1_ivec(t) = r1_i + dr1idt*dt;
    
    dr2edt = (-r2_e + alpha_e*(I2_e(t-1) - theta_e) * sign2*(I2_e(t-1) - theta_e))/tau_e;
    r2_evec(t) = r2_e + dr2edt*dt;
    
    dr2idt = (-r2_i + alpha_i*(I2_i(t-1) - theta_i))/tau_i;
    r2_ivec(t) = r2_i + dr2idt*dt;
    
    if (r1_ivec(t) > r_max)
        r1_ivec(t) = r_max; 
    end
    if (r1_evec(t) > r_max)
        r1_evec(t) = r_max;
    end
    
    if (r1_ivec(t) < 0)
        r1_ivec(t) = 0;
    end
    if (r1_evec(t) < 0)
        r1_evec(t) = 0;
    end
    
    if (r2_ivec(t) > r_max)
        r2_ivec(t) = r_max; 
    end
    if (r2_evec(t) > r_max)
        r2_evec(t) = r_max;
    end
    
    if (r2_ivec(t) < 0)
        r2_ivec(t) = 0;
    end
    if (r2_evec(t) < 0)
        r2_evec(t) = 0;
    end
end

figure
subplot(2,1,1)
plot(tvec, I1_eapp)
hold on
plot(tvec, I2_eapp)
xlabel('Time')
ylabel('Current')
legend('Excitatory Current Unit 1', 'Excitatory Current Unit 2')
subplot(2,1,2)
plot(tvec, r1_evec)
hold on
plot(tvec, r1_ivec)
hold on
plot(tvec, r2_evec)
hold on
plot(tvec, r2_ivec)
xlabel('Time')
ylabel('Firing Rate')
legend('Excitatory Firing Rate Unit 1', 'Inhibitory Firing Rate Unit 1',...
    'Excitatory Firing Rate Unit 2', 'Inhibitory Firing Rate Unit 2')

figure
plot(r1_ivec, r1_evec)
hold on
plot(r2_ivec, r2_evec)
xlabel('Inhibitory Firing Rate')
ylabel('Excitatory Firing Rate')
legend('Unit 1', 'Unit 2')