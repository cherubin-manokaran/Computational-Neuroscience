% Manokaran_Tutorial1_2.m

% Forced Voltage Clamp

% At constant current of 220 picoamps
% Records potential over 100 ms interval
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;

THRESHOLD = -.050;
RESET = -.065;
REFRACTORY_PERIOD = 2.5e-3;

dt = .0001;
tvec = 0:dt:.1;
spiketime = -2.6e-3;                                                        % Initializes variable for recording the time of spikes

constant = 220e-12;                                                         % Sets a constant current value
constant_current = ones(size(tvec));                                        % Creates a constant current vector
constant_current = constant * constant_current;                             % Sets values in vector to constant current 

potential = zeros(size(tvec));                                              % Creates a vector of equal size to the time vector for membrane potential
potential(1) = LEAK_POTENTIAL;                                              % Sets the initial potential to the leak potential

for i = 2:length(tvec)
    val = tvec(i);
    if (val > spiketime + REFRACTORY_PERIOD)                                % Tests whether the repractory period has been exited
       dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + constant_current(i))/MEMBRANE_CAPACITANCE;
       potential(i) = potential(i-1) + dt * dvdt;
       if potential(i) > THRESHOLD                                          % Tests whether potential is greater than threshold
           spiketime = val;                                                 % Sets spike time to the current time
           potential(i) = RESET;                                            % Sets current potential to reset value
           potential(i-1) = PEAK_POTENTIAL;                                 % Sets previous potential value to peak potential
       end
    else
        potential(i) = RESET;                                               % Puts into refractory period
    end
end

figure;
plot(tvec, potential)
hold on

% At constant current of 600 picoamps
% Over 100 ms inteval
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;

THRESHOLD = -.050;
RESET = -.065;
REFRACTORY_PERIOD = 2.5e-3;

dt = .0001;
tvec = 0:dt:.1;
spiketime = -2.6e-3;

constant = 600e-12;                                                         % Sets a constant current value
constant_current = ones(size(tvec));                                        % Creates of a vonstant current vector
constant_current = constant * constant_current; 

potential = zeros(size(tvec));                                              % Creates a vector of equal size to the time vector for membrane potential
potential(1) = LEAK_POTENTIAL;  

for i = 2:length(tvec)
    val = tvec(i);
    if (val > spiketime + REFRACTORY_PERIOD)
       dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + constant_current(i))/MEMBRANE_CAPACITANCE;
       potential(i) = potential(i-1) + dt * dvdt;
       if potential(i) > THRESHOLD
           spiketime = val;
           potential(i) = RESET;
           potential(i-1) = PEAK_POTENTIAL;
       end
    else
        potential(i) = RESET;
    end
end

plot(tvec, potential)
legend('220 picoamps', '600 picoamps')
title('Voltage Clamp')
xlabel('Time')
ylabel('Voltage (V)')

% At variable currents
% For 2 s interval
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;

THRESHOLD = -.050;
RESET = -.065;
REFRACTORY_PERIOD = 2.5e-3;

dt = .0001;
tvec = 0:dt:2;

variable_current = 100e-12:50e-12:600e-12;                                  % Creates a variable current vector
firing_rate = zeros(size(variable_current));                                % Creates a vector to record the firing rate for each trial at a unique current
mean_potential = zeros(size(variable_current));                             % Creates a vector to record the mean potential for each trial

for j = 1:length(variable_current)
    potential = zeros(size(tvec));                                          
    potential(1) = LEAK_POTENTIAL;
    Iapp = variable_current(j);
    spiketime = -2.6e-3;
    total_potential = -.070;                                                % Initializes variable to to calculate total potential for eahc trial
    for i = 2:length(tvec)
        val = tvec(i);
        if (val > spiketime + REFRACTORY_PERIOD)
           dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + Iapp)/MEMBRANE_CAPACITANCE;
           potential(i) = potential(i-1) + dt * dvdt;
           if potential(i) > THRESHOLD
               spiketime = val;
               potential(i) = RESET;
               firing_rate(j) = firing_rate(j) + 1;
           end
        else
            potential(i) = RESET;
        end
        total_potential = total_potential + potential(i);                   % Calculates total potential by adding the current potential
    end
    mean_potential(j) = total_potential/length(potential);                  % Calculates mean potential by dividing the total potential by the total number of potentials recorded
end

% Threshold Increase

% At constant current of 220 picoamps
% Over 100 ms interval
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;

BASE_THRESHOLD = -.050;
MAX_THRESHOLD = .2;
RESET = -.065;

dt = .0001;
tvec = 0:dt:.1;

constant = 220e-12;                                                         % Sets a constant current value
constant_current = ones(size(tvec));                                        % Creates of a vonstant current vector
constant_current = constant * constant_current; 

potential = zeros(size(tvec));                                          
potential(1) = LEAK_POTENTIAL;
threshold = zeros(size(tvec));                                              % Creates a dynamic threshold vector
threshold(1) = BASE_THRESHOLD;                                              % Sets the intital threshold to base threshold

for i = 2:length(tvec)
   dVthdt = (BASE_THRESHOLD - threshold(i-1))/.001;                         % Calculate the change in threshold over the time step
   threshold(i) = threshold(i-1) + dt * dVthdt;                             % Calculates the new threshold
   dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + constant_current(i))/MEMBRANE_CAPACITANCE;
   potential(i) = potential(i-1) + dt * dvdt;
   if potential(i) > threshold(i)                                           % Tests if potential is greater tha threshold
       threshold(i) = MAX_THRESHOLD;                                        % Sets threshold to the maximum threshold
       potential(i) = RESET;                                                % Sets current potential to the reset value
       potential(i-1) = PEAK_POTENTIAL;                                     % Sets the previous potential to the maxium
   end
end

figure;
plot(tvec, potential, 'b')
hold on
% At constant current of 600 picoamps
% Over 100 ms interval
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;

BASE_THRESHOLD = -.050;
MAX_THRESHOLD = .2;
RESET = -.065;

dt = .0001;
tvec = 0:dt:.1;

constant = 600e-12;                                                         
constant_current = ones(size(tvec));                                        
constant_current = constant * constant_current;                             

potential = zeros(size(tvec));                                          
potential(1) = LEAK_POTENTIAL;
threshold = zeros(size(tvec));
threshold(1) = BASE_THRESHOLD;

for i = 2:length(tvec)
   dVthdt = (BASE_THRESHOLD - threshold(i-1))/.001;
   threshold(i) = threshold(i-1) + dt * dVthdt;
   dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + constant_current(i))/MEMBRANE_CAPACITANCE;
   potential(i) = potential(i-1) + dt * dvdt;
   if potential(i) > threshold(i)
       threshold(i) = MAX_THRESHOLD;
       potential(i) = RESET;
       potential(i-1) = PEAK_POTENTIAL;
   end
end

plot(tvec, potential, 'r')
legend('220 picoamps', '600 picoamps')
title('Increased Threshold')
xlabel('Time')
ylabel('Voltage (V)')

% At variable currents
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;

BASE_THRESHOLD = -.050;
MAX_THRESHOLD = .2;
RESET = -.065;

dt = .0001;
tvec = 0:dt:2;

variable_current2 = 100e-12:50e-12:600e-12; 
firing_rate2 = zeros(size(variable_current)); 
mean_potential2 = zeros(size(variable_current));

for j = 1:length(variable_current)
    potential = zeros(size(tvec));                                          
    potential(1) = LEAK_POTENTIAL;
    threshold = zeros(size(tvec));                                          
    threshold(1) = BASE_THRESHOLD;
    Iapp = variable_current2(j);
    total_potential = -.070;
    for i = 2:length(tvec)
       dVthdt = (BASE_THRESHOLD - threshold(i-1))/.001;
       threshold(i) = threshold(i-1) + dt * dVthdt;
       dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + variable_current2(j))/MEMBRANE_CAPACITANCE;
       potential(i) = potential(i-1) + dt * dvdt;
       if potential(i) > threshold(i)
           threshold(i) = MAX_THRESHOLD;
           potential(i) = RESET;
           potential(i-1) = PEAK_POTENTIAL;
           firing_rate2(j) = firing_rate2(j) + 1;
       end
       total_potential = total_potential + potential(i);
    end
    mean_potential2(j) = total_potential/length(potential);
end

% Refractory conductance with threshold increase

% At constant current of 220 picoamps
% Over 100 ms interval
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;
POTASSIUM_POTENTIAL = -.080;
DELTA_CONDUCTANCE = 2e-6;

BASE_THRESHOLD = -.050;
MAX_THRESHOLD = .2;
RESET = -.065;

dt = .0001;
tvec = 0:dt:.1;

constant = 220e-12;                                                         
constant_current = ones(size(tvec));                                        
constant_current = constant * constant_current;                             

potential = zeros(size(tvec));                                          
potential(1) = LEAK_POTENTIAL;
threshold = zeros(size(tvec));
threshold(1) = BASE_THRESHOLD;
kconductance = zeros(size(tvec));                                           % Creates a vector for potassium conductance overtime

for i = 2:length(tvec)
   dGdt = -kconductance(i-1)/.0002;                                         % Calculates the change in potassium conductance over the time step
   kconductance(i) = kconductance(i-1) + dt * dGdt;                         % Calculates the current conductance
   dVthdt = (BASE_THRESHOLD - threshold(i-1))/.001;                         
   threshold(i) = threshold(i-1) + dt * dVthdt;
   dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + (kconductance(i) *(POTASSIUM_POTENTIAL - potential(i-1))) + constant_current(i))/MEMBRANE_CAPACITANCE; % Calculates the change in potential and adds a term for the effect of conductance 
   potential(i) = potential(i-1) + dt * dvdt;
   if potential(i) > threshold(i)
       threshold(i) = MAX_THRESHOLD;
       kconductance(i) = kconductance(i) + DELTA_CONDUCTANCE;               % Conductance is zero until the potential becomes greater than threshold
       potential(i-1) = PEAK_POTENTIAL;
   end
end

figure;
plot(tvec, potential,'b')
hold on

% At constant current of 600 picoamps
% Over 100 ms interval
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;
POTASSIUM_POTENTIAL = -.080;
DELTA_CONDUCTANCE = 2e-6;

BASE_THRESHOLD = -.050;
MAX_THRESHOLD = .2;
RESET = -.065;

dt = .0001;
tvec = 0:dt:.1;

constant = 600e-12;                                                         
constant_current = ones(size(tvec));                                        
constant_current = constant * constant_current;                             

potential = zeros(size(tvec));                                          
potential(1) = LEAK_POTENTIAL;
threshold = zeros(size(tvec));
threshold(1) = BASE_THRESHOLD;
kconductance = zeros(size(tvec));

for i = 2:length(tvec)
   dGdt = -kconductance(i-1)/.0002;
   kconductance(i) = kconductance(i-1) + dt * dGdt;
   dVthdt = (BASE_THRESHOLD - threshold(i-1))/.001;
   threshold(i) = threshold(i-1) + dt * dVthdt;
   dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + (kconductance(i) *(POTASSIUM_POTENTIAL - potential(i-1))) + constant_current(i))/MEMBRANE_CAPACITANCE;
   potential(i) = potential(i-1) + dt * dvdt;
   if potential(i) > threshold(i)
       threshold(i) = MAX_THRESHOLD;
       kconductance(i) = kconductance(i) + DELTA_CONDUCTANCE;
       potential(i-1) = PEAK_POTENTIAL;
   end
end

plot(tvec, potential, 'r')
legend('220 picoamps', '600 picoamps')
title('Potassium Conductance')
xlabel('Time')
ylabel('Voltage (V)')

% At variable currents
% Over 2 s interval
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 100e6;
MEMBRANE_CAPACITANCE = .1e-9;
PEAK_POTENTIAL = .050;
POTASSIUM_POTENTIAL = -.080;
DELTA_CONDUCTANCE = 2e-6;

BASE_THRESHOLD = -.050;
MAX_THRESHOLD = .2;
RESET = -.065;

dt = .0001;
tvec = 0:dt:2;

variable_current3 = 100e-12:50e-12:600e-12; 
firing_rate3 = zeros(size(variable_current)); 
mean_potential3 = zeros(size(variable_current));                           

for j = 1:length(variable_current3) 
    potential = zeros(size(tvec));                                          
    potential(1) = LEAK_POTENTIAL;
    threshold = zeros(size(tvec));
    threshold(1) = BASE_THRESHOLD;
    kconductance = zeros(size(tvec));
    total_potential = -.070;
    for i = 2:length(tvec)
       dGdt = -kconductance(i-1)/.0002;
       kconductance(i) = kconductance(i-1) + dt * dGdt;
       dVthdt = (BASE_THRESHOLD - threshold(i-1))/.001;
       threshold(i) = threshold(i-1) + dt * dVthdt;
       dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + (kconductance(i) *(POTASSIUM_POTENTIAL - potential(i-1))) + variable_current3(j))/MEMBRANE_CAPACITANCE;
       potential(i) = potential(i-1) + dt * dvdt;
       if potential(i) > threshold(i)
           threshold(i) = MAX_THRESHOLD;
           kconductance(i) = kconductance(i) + DELTA_CONDUCTANCE;
           potential(i-1) = PEAK_POTENTIAL;
           firing_rate3(j) = firing_rate3(j) + 1;
       end
       total_potential = total_potential + potential(i);
    end
    mean_potential3(j) = total_potential/length(potential);
end

% figure;
% plot(tvec, potential)
% xlabel('Time')
% ylabel('Voltage (V)')
firing_rate = firing_rate/2;
firing_rate2 = firing_rate2 / 2;
firing_rate3 = firing_rate3 / 2;

figure;
plot(variable_current, firing_rate, 'b')
hold on
plot(variable_current2, firing_rate2, 'r')
hold on
plot(variable_current3, firing_rate3, 'y')
legend('Voltage Clamp', 'Increased Threshold', 'Potassium Conductance')
xlabel('Current (amp)')
ylabel('Firing Rate')

figure;
plot(variable_current, mean_potential, 'b')
hold on
plot(variable_current2, mean_potential2, 'r')
hold on
plot(variable_current3, mean_potential3, 'y')
legend('Voltage Clamp', 'Increased Threshold', 'Potassium Conductance')
xlabel('Current (amp)')
ylabel('Mean Potential')

figure;
plot(firing_rate, mean_potential, 'b')
hold on
plot(firing_rate2, mean_potential2, 'r')
hold on
plot(firing_rate3, mean_potential3, 'y')
legend('Voltage Clamp', 'Increased Threshold', 'Potassium Conductance')
xlabel('Firing Rate')
ylabel('Mean Potential')



