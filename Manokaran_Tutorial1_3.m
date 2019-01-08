% Manokaran_Tutorial1_3

% At constant current of 600 picoamps
% Over 100 ms interval

LEAK_POTENTIAL = -.075;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = 100e-12;
PEAK_POTENTIAL = .050;
POTASSIUM_POTENTIAL = -.080;
DELTA_CONDUCTANCE = 1e-9;

RESET = -.080;
THRESHOLD = -.050;

dt = .0001;
tvec = 0:dt:1.5;

constant = 500e-12;                                                         
constant_current = zeros(size(tvec));                                                                   

potential = zeros(size(tvec));                                          
potential(1) = LEAK_POTENTIAL;
conductance = zeros(size(tvec));

for i = 2:length(tvec)
   time = tvec(i-1);
   if (time >= .5 && time <= 1)
        constant_current(i) = constant;
   end
   dGsradt = -conductance(i-1)/.2;
   conductance(i) = conductance(i-1) + dt * dGsradt;
   dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + (conductance(i) *(POTASSIUM_POTENTIAL - potential(i-1))) + constant_current(i))/MEMBRANE_CAPACITANCE;
   potential(i) = potential(i-1) + dt * dvdt;
   if potential(i) > THRESHOLD
       conductance(i) = conductance(i) + DELTA_CONDUCTANCE;
       potential(i-1) = PEAK_POTENTIAL;
       potential(i) = RESET;
   end
end

subplot(3,1,2)
plot(tvec, constant_current)
xlabel('Current')
ylabel('Voltage (V)')

subplot(3,1,1)
plot(tvec, potential)
xlabel('Time')
ylabel('Voltage (V)')

subplot(3,1,3)
plot(tvec, conductance)
xlabel('Time')
ylabel('Conductance')

%

LEAK_POTENTIAL = -.075;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = 100e-12;
PEAK_POTENTIAL = .050;
POTASSIUM_POTENTIAL = -.080;
DELTA_CONDUCTANCE = 1e-9;

RESET = -.080;
THRESHOLD = -.050;

dt = .0001;
tvec = 0:dt:5;

I_app = zeros(size(tvec));
app_current = 200e-12:20e-12:600e-12;
firing_rate = zeros(size(app_current));
ISI = zeros(size(app_current));

for j = 1:length(app_current)
    potential = zeros(size(tvec));                                          
    potential(1) = LEAK_POTENTIAL;
    conductance = zeros(size(tvec));
    for i = 2:length(tvec)
       time = tvec(i-1);
       if (time >= .5 && time <= 1)
            I_app(i) = app_current(j);
       end
       dGsradt = -conductance(i-1)/.2;
       conductance(i) = conductance(i-1) + dt * dGsradt;
       dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + (conductance(i) *(POTASSIUM_POTENTIAL - potential(i-1))) + I_app(i))/MEMBRANE_CAPACITANCE;
       potential(i) = potential(i-1) + dt * dvdt;
       if potential(i) > THRESHOLD
           conductance(i) = conductance(i) + DELTA_CONDUCTANCE;
           potential(i) = RESET;
           firing_rate(j) = firing_rate(j) + 1;
       end
    end
    ISI(j) = 1/(MEMBRANE_RESISTANCE * MEMBRANE_CAPACITANCE * log((LEAK_POTENTIAL + (app_current(j)/LEAK_CONDUCTANCE) - RESET)/(LEAK_POTENTIAL + (app_current(j)/LEAK_CONDUCTANCE) - THRESHOLD)));
end
    
firing_rate = firing_rate ./ 5;
figure;
plot(app_current, firing_rate);
hold on 
plot(app_current, ISI, 'o');
hold on
legend('Firing Rate', 'ISI')
title('LIF Model with Adaptation')
xlabel('Current')
ylabel('Firing Rate')

%

LEAK_POTENTIAL = -.075;
LEAK_CONDUCTANCE = 10e-9;
MEMBRANE_CAPACITANCE = 100e-12;

a = 2e-9;
b = .02e-9;
SRA_TAU = .2;


RESET = -.080;
THRESHOLD = -.050;
V_MAX = .050;
THRESHOLD_CHANGE = .002;


dt = .0001;
tvec = 0:dt:1.5;

constant = 500e-12;                                                         
constant_current = zeros(size(tvec));                                                                   

potential = zeros(size(tvec));                                          
potential(1) = LEAK_POTENTIAL;
sra = zeros(size(tvec));

for i = 2:length(tvec)
   time = tvec(i-1);
   if (time >= .5 && time <= 1)
        constant_current(i) = constant;
   end
   dSRAdt = ((a *(potential(i) - LEAK_POTENTIAL) - sra(i-1))/SRA_TAU);
   sra(i) = sra(i-1) + dt * dSRAdt;
   dvdt = (LEAK_CONDUCTANCE * (LEAK_POTENTIAL - potential(i-1) + THRESHOLD_CHANGE * exp((potential(i-1)- THRESHOLD)/THRESHOLD_CHANGE)) - sra(i) + constant_current(i))/MEMBRANE_CAPACITANCE;
   potential(i) = potential(i-1) + dt * dvdt;
   if potential(i) > V_MAX
       potential(i) = RESET;
       sra(i) = sra(i) + b;
   end
end

figure;

subplot(3,1,2)
plot(tvec, constant_current)
xlabel('Current')
ylabel('Voltage (V)')

subplot(3,1,1)
plot(tvec, potential)
xlabel('Time')
ylabel('Voltage (V)')

subplot(3,1,3)
plot(tvec, sra)
xlabel('Time')
ylabel('Voltage (V)')

%

LEAK_POTENTIAL = -.075;
LEAK_CONDUCTANCE = 10e-9;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = 100e-12;
PEAK_POTENTIAL = .050;
POTASSIUM_POTENTIAL = -.080;

DELTA_CONDUCTANCE = 1e-9;
a = 2e-9;
b = .02e-9;
SRA_TAU = .2;


RESET = -.080;
THRESHOLD = -.050;
V_MAX = .050;
THRESHOLD_CHANGE = .002;


dt = .0001;
tvec = 0:dt:5;

I_app = zeros(size(tvec));
app_current = 200e-12:20e-12:600e-12;
firing_rate = zeros(size(app_current));
ISI = zeros(size(app_current));

for j = 1:length(app_current)
    potential = zeros(size(tvec));                                          
    potential(1) = LEAK_POTENTIAL;
    sra = zeros(size(tvec));
    for i = 2:length(tvec)
       time = tvec(i-1);
       if (time >= .5 && time <= 1)
            I_app(i) = app_current(j);
       end
       dSRAdt = ((a *(potential(i) - LEAK_POTENTIAL) - sra(i-1))/SRA_TAU);
       sra(i) = sra(i-1) + dt * dSRAdt;
       dvdt = (LEAK_CONDUCTANCE * (LEAK_POTENTIAL - potential(i-1) + THRESHOLD_CHANGE * exp((potential(i-1)- THRESHOLD)/THRESHOLD_CHANGE)) - sra(i) + I_app(i))/MEMBRANE_CAPACITANCE;
       potential(i) = potential(i-1) + dt * dvdt;
       if potential(i) > V_MAX
           potential(i) = RESET;
           sra(i) = sra(i) + b;
           firing_rate(j) = firing_rate(j) + 1;
       end
    end
    ISI(j) = 1/(MEMBRANE_RESISTANCE * MEMBRANE_CAPACITANCE * log((LEAK_POTENTIAL + (app_current(j)/LEAK_CONDUCTANCE) - RESET)/(LEAK_POTENTIAL + (app_current(j)/LEAK_CONDUCTANCE) - THRESHOLD)));
end

firing_rate = firing_rate ./ 5;
figure;
plot(app_current, firing_rate);
hold on
plot(app_current, ISI, 'o');
title('AELIF Model')
legend('Firing Rate', 'ISI');
xlabel('Current')
ylabel('Firing rate')