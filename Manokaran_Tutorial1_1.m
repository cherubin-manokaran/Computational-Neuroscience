%% part a
% Model a leaky integrate-and-fire neuron
LEAK_POTENTIAL = -.070;
MEMBRANE_RESISTANCE = 5e6;
MEMBRANE_CAPACITANCE = 2e-9;

THRESHOLD = -.050;
RESET = -.065;

dt = .00001;                                                            % Sets the time step to .1 ms

tvec = 0:dt:2;                                                          % Creates a 2 s time vector with steps of .1 ms
potential = zeros(size(tvec));                                          % Creates a vector of equal size to the time vector for membrane potential
potential(1) = LEAK_POTENTIAL;                                          % Sets the initial membrane potential to the leak potential
constant = 4.05e-9;                                                     % Sets a constant current value
constant_current = ones(size(tvec));                                    % Creates of a vonstant current vector
constant_current = constant * constant_current;                         % Sets each value in the vector to the constant current
variable_current = zeros(1,10);                                         % Creates a variable current vector

for i = 2:length(tvec)
   dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + constant_current(i))/MEMBRANE_CAPACITANCE;
   potential(i) = potential(i-1) + dt * dvdt;
   if potential(i) > THRESHOLD
       potential(i) = RESET;
   end
end
figure;
plot(tvec, potential)
xlabel('Time')
ylabel('Voltage (V)')
drawnow

%% part b
% Appplies constant current over a 200 ms interval
% Trials will include currents above and below threshold current
tvec = 0:dt:.2;                                                         % Creates a 200 ms time vector with a steps of .1 ms
potential = zeros(size(tvec));                                          % Creates a vector of equal size to the time vector for membrane potential
potential(1) = LEAK_POTENTIAL;                                          % Sets the initial membrane potential to the leak potential
constant = 4.01e-9;                                                     % Sets a constant current value
constant_current = ones(size(tvec)) * constant;                         % Creates of a constant current vector

for i = 2:length(tvec)                                                  % Integrates over two seconds in .1 ms intervals
   dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + constant_current(i))/MEMBRANE_CAPACITANCE; % Calculates membrane potential over the time step
   potential(i) = potential(i-1) + dt * dvdt;                           % Calculates the resulting membrane potential
   if potential(i) > THRESHOLD                                          % Tests whether the membrane potential is above the threshold
       potential(i) = RESET;                                            % If true, sets the potential to the reset potential
   end
end
figure;                                                                 % Creates new figure
plot(tvec, potential)
xlabel('Time')
ylabel('Voltage (V)')
drawnow

%% Part 1c, d & 2

% Runs 10 trials using variable currents, plots the resulting changes 
% in potential and counts the number of spikes to calculate the 
% firing rate. Also implesments the firing rate equation provided and
% plots both the counted and calculated firing rates. Adds an changes
% the noise term. Have also adjusted the time step for the last part.

dt = .0000001;                                                           % Previously set the time step to .1 ms but for last set of trials decreased the time step by a factor of ten
tvec = 0:dt:2;                                                          % Creates a 2 s time vector with steps of .1 ms
potential = zeros(size(tvec));                                          % Creates a vector of equal size to the time vector for membrane potential
potential(1) = LEAK_POTENTIAL;                                          % Sets the initial membrane potential to the leak potential
variable_current = zeros(1,10);                                         % Creates a variable current vector
firing_rate = zeros(1,10);                                              % Creates a firing rate vector
firing_rate_c = zeros(1,10);                                            % Creates another firing rate vector for those calculated using the given equation
variable_current(1) = 4e-9;                                             % Sets initial input to threshold current.
sigma_I = 1e-11;                                                        % Sets sigma_I for noise term

for j = 1:10                                                            % Runs ten trials with variable currents
    for i = 2:length(tvec)
        dvdt = (((LEAK_POTENTIAL - potential(i-1))/ MEMBRANE_RESISTANCE) + variable_current(j))/MEMBRANE_CAPACITANCE;    
        potential(i) = potential(i-1) + dt * dvdt + randn()*sigma_I*sqrt(dt);
        % potential(i) = potential(i-1) + dt * dvdt;
        if potential(i) > THRESHOLD
           potential(i) = RESET;
           firing_rate(j) = firing_rate(j) + 1;                         % If above threshold, also increase firing_rate for trial by 1
        end
    end
    figure;
    plot(tvec, potential)
    xlabel('Time')
    ylabel('Voltage (V)')
    title(['Applied Current = ',num2str(variable_current(j))])
    drawnow                                                             % Directs program to draw plot immediately
    firing_rate_c(j) =(MEMBRANE_RESISTANCE * MEMBRANE_CAPACITANCE * log((variable_current(j)* MEMBRANE_RESISTANCE) + LEAK_POTENTIAL - RESET)) - (MEMBRANE_RESISTANCE * MEMBRANE_CAPACITANCE * log((variable_current(j)* MEMBRANE_RESISTANCE) + LEAK_POTENTIAL - THRESHOLD));
    if j < 10                                                           % Enusres that eleventh trial is not run
        variable_current(j+1) = variable_current(j) + 5e-11;            % Increases current for next trial
    end
end

firing_rate = firing_rate ./ 2;
firing_rate_c = 1 ./ firing_rate_c;                                     % Calculates firing rate using provided equation by finding the inverse of the values caluclated above
figure;
plot(variable_current, firing_rate, 'O')
hold on
plot(variable_current, firing_rate_c, 'X')
legend('firing rate counted', 'firing rate calculated')
xlabel('Current (amp)')
ylabel('Firing Rate')