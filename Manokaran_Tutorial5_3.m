%% Part 1

V_th = -.050;
V_reset = -.080;
sigmaV = .001;
tau = .003;
tau_e = .002;
tau_i = .005;
E_l = -.070;
E_i = -.065;
E_e = 0;
G_l = .05e-9;

G_ee = 25e-9;
G_ei = 4e-9;
G_ie = 800e-9;

dt = .0001;
tvec = 0:dt:2.5;

G1_in = 1e-9;
G2_in = 0;

G1_e = zeros(size(tvec));
G1_i = zeros(size(tvec));
G2_e = zeros(size(tvec));
G2_i = zeros(size(tvec));

S1_e = zeros(size(tvec));
S2_i = zeros(size(tvec));

r1_i = zeros(size(tvec));
r2_i = zeros(size(tvec));

r2_ss = zeros(size(tvec));
r1_ss = zeros(size(tvec));
V1_ssi = zeros(size(tvec));
V2_ssi = zeros(size(tvec));

osc1 = 0;
osc2 = 0;

for t = 1:length(tvec)-1
    alpha = 0.2;

    G1_e(t) = G_ee*S1_e(t) + G1_in;
    G1_i(t) = G_ie*S2_i(t);
    G2_e(t) = G_ei*S1_e(t) + G2_in;
    G2_i(t) = 0;
    
    V1_ssi(t) = (G_l*E_l + G1_i(t)*E_i + G1_e(t)*E_e)/(G_l + G1_i(t) +G1_e(t));
    V2_ssi(t) = (G_l*E_l + G2_i(t)*E_i + G2_e(t)*E_e)/(G_l + G2_i(t) +G2_e(t));

    if (V1_ssi(t) == V_th)
        r1_ss(t) = sigmaV*(tau*(V1_ssi(t)-Vreset));
    else
        r1_ss(t) = (V1_ssi(t) - V_th)/(tau*(V_th-V_reset)*(1-exp(-(V1_ssi(t)-V_th)/sigmaV))); 
    end
    
    if (V2_ssi(t) == V_th)
        r2_ss(t) = sigmaV*(tau*(V2_ssi(t)-V_reset));
    else
        r2_ss(t) = (V2_ssi(t) - V_th)/(tau*(V_th-V_reset)*(1-exp(-(V2_ssi(t)-V_th)/sigmaV)));
    end
    
    if (r1_i(t) == r1_ss(t))
        osc1 = osc1+1;
    end
    
    dr1idt = (-r1_i(t)+r1_ss(t))/tau;
    r1_i(t+1) = r1_i(t) + dt*dr1idt;
    
    if (r2_i(t) == r2_ss(t))
        osc2 = osc2+1;
    end
    
    dr2idt = (-r2_i(t)+r2_ss(t))/tau;
    r2_i(t+1) = r2_i(t) + dt*dr2idt;

    dsedt = -(S1_e(t)/tau_e) + alpha*r1_i(t)*(1-S1_e(t));
    S1_e(t+1) = S1_e(t) + dsedt*dt;
    dsidt = -(S2_i(t)/tau_i) + alpha*r2_i(t)*(1-S2_i(t));
    S2_i(t+1) = S2_i(t) + dsidt*dt;
end

figure
plot(tvec,r1_i,'--')
hold on
plot(tvec,r2_i)
title('Part 1: Firing Rate over Time')
xlabel('Time')
ylabel('Firing Rate')
legend('Unit 1','Unit 2')

%% Part 2

V_th = -.050;
V_reset = -.080;
sigmaV = .001;
tau = .003;
tau_e = .002;
tau_i = .005;
E_l = -.070;
E_i = -.065;
E_e = 0;
G_l = .05e-9;

G_ee = 25e-9;
G_ei = 4e-9;
G_ie = 800e-9;

dt = .0001;
tvec = 0:dt:2.5;

G1_in = 1e-9;
G2_in = 0;

G1_e = zeros(size(tvec));
G1_i = zeros(size(tvec));
G2_e = zeros(size(tvec));
G2_i = zeros(size(tvec));

S1_e = zeros(size(tvec));
S2_i = zeros(size(tvec));

r1_i = zeros(size(tvec));
r2_i = zeros(size(tvec));

r2_ss = zeros(size(tvec));
r1_ss = zeros(size(tvec));
V1_ssi = zeros(size(tvec));
V2_ssi = zeros(size(tvec));

r_part2 = zeros(size(tvec));

osc1 = 0;
osc2 = 0;

peak_times = zeros(size(tvec));

for t = 1:length(tvec)-1
    alpha = 0.2;

    G1_e(t) = G_ee*S1_e(t) + G1_in;
    G1_i(t) = G_ie*S2_i(t);
    G2_e(t) = G_ei*S1_e(t) + G2_in;
    G2_i(t) = 0;
    
    V1_ssi(t) = (G_l*E_l + G1_i(t)*E_i + G1_e(t)*E_e)/(G_l + G1_i(t) +G1_e(t));
    V2_ssi(t) = (G_l*E_l + G2_i(t)*E_i + G2_e(t)*E_e)/(G_l + G2_i(t) +G2_e(t));

    if (V1_ssi(t) == V_th)
        r1_ss(t) = sigmaV*(tau*(V1_ssi(t)-Vreset));
    else
        r1_ss(t) = (V1_ssi(t) - V_th)/(tau*(V_th-V_reset)*(1-exp(-(V1_ssi(t)-V_th)/sigmaV))); 
    end
    
    if (V2_ssi(t) == V_th)
        r2_ss(t) = sigmaV*(tau*(V2_ssi(t)-V_reset));
    else
        r2_ss(t) = (V2_ssi(t) - V_th)/(tau*(V_th-V_reset)*(1-exp(-(V2_ssi(t)-V_th)/sigmaV)));
    end
    
    if (r1_i(t) == r1_ss(t))
        osc1 = osc1+1;
    end
    
    dr1idt = (-r1_i(t)+r1_ss(t))/tau;
    r1_i(t+1) = r1_i(t) + dt*dr1idt;

    if (r2_i(t) == r2_ss(t))
        osc2 = osc2+1;
    end
    
    dr2idt = (-r2_i(t)+r2_ss(t))/tau;
    r2_i(t+1) = r2_i(t) + dt*dr2idt;

    dsedt = -(S1_e(t)/tau_e) + alpha*r1_i(t)*(1-S1_e(t));
    S1_e(t+1) = S1_e(t) + dsedt*dt;
    dsidt = -(S2_i(t)/tau_i) + alpha*r2_i(t)*(1-S2_i(t));
    S2_i(t+1) = S2_i(t) + dsidt*dt;
end

r_part2(5001:25001) = r1_i(5001:25001);
high = max(r_part2) - .2;
low = min(r_part2(5001:25001)) + .2;

indicator = 1;
time = zeros(size(r_part2));

for i = 5001:length(r_part2)
    if (indicator == 0 && r_part2(i) > high)
        time(i) = tvec(i);
        high = r_part2(i-2);
        indicator = 1;
    end 
    if (r_part2(i) < low)
        low = r_part2(i-2);
        indicator = 0;
    end
end

osc_amp = max(r_part2);
oscillation_times = find(time);
osc_period = ((oscillation_times(length(oscillation_times)) - oscillation_times(1))*dt)/length(oscillation_times);

figure
plot(tvec,r_part2)
title('Part 2: Firing Rate over Time')
xlabel('Time')
ylabel('Firing Rate')

display('Part 2:')
fprintf('Oscillation Frequency = %i\n', 1/osc_period)
fprintf('Oscillation Amplitude = %i\n', osc_amp)

%% Part 3
sigma_v = 0.001;
tmax = 2.5;
dt = 0.1e-3;
time_omit = 0.5;
tvec = 0:dt:tmax;
G_e1 = zeros(size(tvec));
G_i1 = zeros(size(tvec));
G_e2 = zeros(size(tvec));
G_i2 = 0;

G_ee = 25e-9;
G_ei = 4e-9;
G_ie = 800e-9;
G_in1 = 1e-9;
G_in2 = 0e-9;
r1 = zeros(size(tvec));
r2 = zeros(size(tvec));
r_part3 = zeros(size(tvec));
tau_i = 5e-3;
tau_e = 2e-3;
alpha = 0.2;

s_e1 = zeros(size(tvec));
s_i2 = zeros(size(tvec));

osc_count1 = 0;
osc_count2 = 0;
time_start1 = tmax;
time_start2 = tmax;

threshold = 1;

frequency = 0:0.5:100; 

time = zeros(size(tvec));
peak = max(r_part2);
trough = min(r_part2);
thresh_high = peak;
thresh_low = trough;


for i = 2:length(tvec)
    
    
    V_ss1 = (G_l * E_l + G_i1(i-1)*E_i + G_e1(i-1) * E_e)/(G_l + G_i1(i-1) + G_e1(i-1));
    
    if (V_ss1 == V_th)
        r_ss1 = sigma_v/(tau*(V_ss1 - V_reset));
        
    else
        r_ss1 = (V_ss1 - V_th) / (tau*(V_th - V_reset) * (1 - exp(-(V_ss1 - V_th)/sigma_v)));
    end
    
    
    V_ss2 = (G_l * E_l + G_i2*E_i + G_e2(i-1) * E_e)/(G_l + G_i2 + G_e2(i-1));
    if (V_ss2 == V_th)
        r_ss2 = sigma_v/(tau*(V_ss2 - V_reset));
    else
        r_ss2 = (V_ss2 - V_th) / (tau*(V_th - V_reset) * (1 - exp(-(V_ss2 - V_th)/sigma_v)));
    end
    
    
    drdt1 = (-r1(i-1) + r_ss1)/tau;
    r1(i) = r1(i-1) + dt*drdt1;
    drdt2 = (-r2(i-1) + r_ss2)/tau;
    r2(i) = r2(i-1) + dt*drdt2;
    
    dse1dt = -s_e1(i-1)/tau_e + alpha * r1(i-1) *(1-s_e1(i-1));
    s_e1(i) = s_e1(i-1) + dt*dse1dt;
    
    dsi2dt = -s_i2(i-1)/tau_i + alpha * r2(i-1) *(1-s_i2(i-1));
    s_i2(i) = s_i2(i-1) + dt*dsi2dt;
    
    G_e1(i) = G_ee * s_e1(i-1) + G_in1;
    
    G_i1(i) = G_ie * s_i2(i-1);
    
    G_e2(i) = G_ei * s_e1(i-1) + G_in2;


    if(tvec(i)>time_omit)
       r_part3(i) = r1(i);
       
    end
    
end

vec1 = zeros(length(frequency),length(tvec));
vec2 = zeros(length(frequency),length(tvec));
A = zeros(1,length(frequency));
B = zeros(1,length(frequency));


for f = 1:length(frequency)
    for t = 5000:25001
        vec1(f,t) = sin(2*pi*frequency(f)*tvec(t));
        vec2(f,t) = cos(2*pi*frequency(f)*tvec(t));
        
    end
    A(f) = mean(r_part3 * vec1(f,:)');
    B(f) = mean(r_part3 * vec2(f,:)');
end


power = A.^2 + B.^2;

value = max(power(2:201));

find( power == max(power(2:201)));

[M,I] = max(power(2:201));

index = I/2;

figure
plot(frequency,power)
xlabel('Frequency (Hz)')
ylabel('Power')
title('Part 3')

%% Part 4

G_in1 = 0e-9:1e-9:10e-9;


r1 = zeros(length(G_in1),length(tvec));
r2 = zeros(size(r1));
r3_part4 = zeros(size(r1));
r4_part4 = zeros(size(r1));
threshold1 = 1;
threshold2 = 1;
thresh_high1 = max(r1);
thresh_low1 = min(r1);
thresh_high2 = max(r2);
thresh_low2 = min(r2);
oscamp1 = zeros(size(G_in1));
oscamp2 = zeros(size(G_in1));


time4= zeros(size(tvec));
time5=zeros(size(tvec));

oscfreq1 = zeros(size(G_in1));
oscfreq2 = zeros(size(G_in1));

mean_freq1 = zeros(size(G_in1));
mean_freq2 = zeros(size(G_in1));

for j = 1:length(G_in1)
    for i = 2:length(tvec)
        V_ss1 = (G_l * E_l + G_i1(i-1)*E_i + G_e1(i-1) * E_e)/(G_l + G_i1(i-1) + G_e1(i-1));
        if (V_ss1 == V_th)
            r_ss1 = sigma_v/(tau*(V_ss1 - V_reset));
        else
            r_ss1 = (V_ss1 - V_th) / (tau*(V_th - V_reset) * (1 - exp(-(V_ss1 - V_th)/sigma_v)));
        end

        V_ss2 = (G_l * E_l + G_i2*E_i + G_e2(i-1) * E_e)/(G_l + G_i2 + G_e2(i-1));
        if (V_ss2 == V_th)
            r_ss2 = sigma_v/(tau*(V_ss2 - V_reset));
        else
            r_ss2 = (V_ss2 - V_th) / (tau*(V_th - V_reset) * (1 - exp(-(V_ss2 - V_th)/sigma_v)));
        end

        drdt1 = (-r1(j,i-1) + r_ss1)/tau;
        r1(j,i) = r1(j,i-1) + dt*drdt1;
        drdt2 = (-r2(j,i-1) + r_ss2)/tau;
        r2(j,i) = r2(j,i-1) + dt*drdt2;

        dse1dt = -s_e1(i-1)/tau_e + alpha * r1(j,i-1) *(1-s_e1(i-1));
        s_e1(i) = s_e1(i-1) + dt*dse1dt;

        dsi2dt = -s_i2(i-1)/tau_i + alpha * r2(j,i-1) *(1-s_i2(i-1));
        s_i2(i) = s_i2(i-1) + dt*dsi2dt;

        G_e1(i) = G_ee * s_e1(i-1) + G_in1(j);

        G_i1(i) = G_ie * s_i2(i-1);

        G_e2(i) = G_ei * s_e1(i-1) + G_in2;

        if(tvec(i)>time_omit)
           r3_part4(j,i) = r1(j,i);
           
           if(threshold1 == 0 && r3_part4(j,i) > thresh_high1)
              time4(i) = tvec(i);
              threshold1 = 1;
           end
           
           if(r3_part4(j,i) < thresh_low1)
              threshold1 = 0;
           end

           if(r3_part4(j,i) < r3_part4(j,i-1) && r3_part4(j,i-1) > r3_part4(j,i-2))
              thresh_high1 = r3_part4(j,i-2);
           end
           
           if(r3_part4(j,i) > r3_part4(j,i-1) && r3_part4(j,i-1) < r3_part4(j,i-2))
              thresh_low1 = r3_part4(j,i-2);
           end

           r4_part4(j,i) = r2(j,i);
           
           if(threshold1 == 0 && r4_part4(j,i) > thresh_high2)
              time5(i) = tvec(i);
              threshold2 = 1;
           end
           
           if(r4_part4(j,i) < thresh_low2)
              threshold2 = 0;
           end

           if(r4_part4(j,i) < r4_part4(j,i-1) && r4_part4(j,i-1) > r4_part4(j,i-2))
              thresh_high2 = r4_part4(j,i-2);
           end

           if(r4_part4(j,i) > r4_part4(j,i-1) && r4_part4(j,i-1) < r4_part4(j,i-2))
              thresh_low2 = r4_part4(j,i-2);
           end
        end
    end
    
    mean_freq1(j) = (mean(r3_part4(j,5000:25001))/tmax);
    mean_freq2(j) = (mean(r4_part4(j,5000:25001))/tmax);
    osctimes1 = find(time4);
    osctimes2 = find(time5);
    
    if (osctimes1 ~= 0)
        oscfreq1(j) = ((osctimes1(end) - osctimes1(1)) * dt)/length(osctimes1);
    end
    if (osctimes2 ~= 0)
        oscfreq2(j) = ((osctimes2(end) - osctimes2(1)) * dt)/length(osctimes2);
    end
    
    oscamp1(j) = max(r3_part4(j,:));
    oscamp2(j) = max(r4_part4(j,:));
    
end

figure
plot(G_in1,oscfreq1)
hold on
plot(G_in1,oscfreq2)
xlabel('G_ in1 (S)')
ylabel('Oscillation Frequency (Hz)')
title('Part 4: Oscillation Frequency against Excitatory Cell Conductance')
legend('Unit 1', 'Unit 2')

figure
plot(G_in1,oscamp1)
hold on
plot(G_in1,oscamp2)
xlabel('G_in1 (S)')
ylabel('Oscillation Amplitude')
title('Part 4: Amplitude against Excitatory Cell Conductance')
legend('Unit 1', 'Unit 2')

figure
plot(G_in1,mean_freq1)
hold on
plot(G_in1,mean_freq2)
xlabel('G_in1 (S)')
ylabel('Mean Firing Rate (Hz)')
title('Part 4: Firing Rate against Excitatory Cell Conductance')
legend('Unit 1', 'Unit 2')

%% Part 5

G_in1 = 2e-9;
G_in2 = 0e-9:1e-9:10e-9;

r1 = zeros(length(G_in2),length(tvec));
r2 = zeros(size(r1));
r3_part4 = zeros(size(r1));
r4_part4 = zeros(size(r1));
threshold1 = 1;
threshold2 = 1;
thresh_high1 = max(r1);
thresh_low1 = min(r1);
thresh_high2 = max(r2);
thresh_low2 = min(r2);
oscamp1 = zeros(size(G_in2));
oscamp2 = zeros(size(G_in2));

time4= zeros(size(tvec));
time5=zeros(size(tvec));

oscfreq1 = zeros(size(G_in2));

mean_freq1 = zeros(size(G_in2));
mean_freq2 = zeros(size(G_in2));

for j = 1:length(G_in2)
    for i = 2:length(tvec)
        V_ss1 = (G_l * E_l + G_i1(i-1)*E_i + G_e1(i-1) * E_e)/(G_l + G_i1(i-1) + G_e1(i-1));
        if (V_ss1 == V_th)
            r_ss1 = sigma_v/(tau*(V_ss1 - V_reset));
        else
            r_ss1 = (V_ss1 - V_th) / (tau*(V_th - V_reset) * (1 - exp(-(V_ss1 - V_th)/sigma_v)));
        end

        V_ss2 = (G_l * E_l + G_i2*E_i + G_e2(i-1) * E_e)/(G_l + G_i2 + G_e2(i-1));
        if (V_ss2 == V_th)
            r_ss2 = sigma_v/(tau*(V_ss2 - V_reset));
        else
            r_ss2 = (V_ss2 - V_th) / (tau*(V_th - V_reset) * (1 - exp(-(V_ss2 - V_th)/sigma_v)));
        end

        drdt1 = (-r1(j,i-1) + r_ss1)/tau;
        r1(j,i) = r1(j,i-1) + dt*drdt1;
        drdt2 = (-r2(j,i-1) + r_ss2)/tau;
        r2(j,i) = r2(j,i-1) + dt*drdt2;

        dse1dt = -s_e1(i-1)/tau_e + alpha * r1(j,i-1) *(1-s_e1(i-1));
        s_e1(i) = s_e1(i-1) + dt*dse1dt;

        dsi2dt = -s_i2(i-1)/tau_i + alpha * r2(j,i-1) *(1-s_i2(i-1));
        s_i2(i) = s_i2(i-1) + dt*dsi2dt;

        G_e1(i) = G_ee * s_e1(i-1) + G_in1;

        G_i1(i) = G_ie * s_i2(i-1);

        G_e2(i) = G_ei * s_e1(i-1) + G_in2(j);

        if(tvec(i)>time_omit)
           r3_part4(j,i) = r1(j,i);
           
           if(threshold1 == 0 && r3_part4(j,i) > thresh_high1)
              time4(i) = tvec(i);
              threshold1 = 1;
           end
           
           if(r3_part4(j,i) < thresh_low1)
              threshold1 = 0;
           end

           if(r3_part4(j,i) < r3_part4(j,i-1) && r3_part4(j,i-1) > r3_part4(j,i-2))
              thresh_high1 = r3_part4(j,i-2);
           end
           
           if(r3_part4(j,i) > r3_part4(j,i-1) && r3_part4(j,i-1) < r3_part4(j,i-2))
              thresh_low1 = r3_part4(j,i-2);
           end

           r4_part4(j,i) = r2(j,i);
           
           if(threshold1 == 0 && r4_part4(j,i) > thresh_high2)
              time5(i) = tvec(i);
              threshold2 = 1;
           end
           
           if(r4_part4(j,i) < thresh_low2)
              threshold2 = 0;
           end

           if(r4_part4(j,i) < r4_part4(j,i-1) && r4_part4(j,i-1) > r4_part4(j,i-2))
              thresh_high2 = r4_part4(j,i-2);
           end

           if(r4_part4(j,i) > r4_part4(j,i-1) && r4_part4(j,i-1) < r4_part4(j,i-2))
              thresh_low2 = r4_part4(j,i-2);
           end
        end
    end
    
    mean_freq1(j) = (mean(r3_part4(j,5000:25001))/tmax);
    mean_freq2(j) = (mean(r4_part4(j,5000:25001))/tmax);
    osctimes1 = find(time4);
    osctimes2 = find(time5);
    
    if (osctimes1 ~= 0)
        oscfreq1(j) = ((osctimes1(end) - osctimes1(1)) * dt)/length(osctimes1);
    end
    if (osctimes2 ~= 0)
        oscfreq2(j) = ((osctimes2(end) - osctimes2(1)) * dt)/length(osctimes2);
    end
    
    oscamp1(j) = max(r3_part4(j,:));
    oscamp2(j) = max(r4_part4(j,:));
    
end

figure
plot(G_in2,oscfreq1)
hold on
plot(G_in2,oscfreq2)
xlabel('G_ in1 (S)')
ylabel('Oscillation Frequency (Hz)')
title('Part 5: Oscillation Frequency against Inhibitory Cell Conductance')
legend('Unit 1', 'Unit 2')

figure
plot(G_in2,oscamp1)
hold on
plot(G_in2,oscamp2)
xlabel('G_in1 (S)')
ylabel('Oscillation Amplitude')
title('Part 5: Amplitude against Inhibitory Cell Conductance')
legend('Unit 1', 'Unit 2')

figure
plot(G_in2,mean_freq1)
hold on
plot(G_in2,mean_freq2)
xlabel('G_in1 (S)')
ylabel('Mean Firing Rate (Hz)')
title('Part 5: Firing Rate against Inhibitory Cell Conductance')
legend('Unit 1', 'Unit 2')




