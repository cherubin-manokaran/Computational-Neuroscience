LEAK_POTENTIAL = -.06;
LEAK_CONDUCTANCE = 10e-9;
MEMBRANE_CAPACITANCE = 100e-12;

a = 10e-9; % nS
b = .5e-9; % nA
SRA_TAU = .05;

RESET = -.080;
THRESHOLD = -.050;
V_MAX = .050;
THRESHOLD_CHANGE = .002;

dt = .00002;
T_ZERO = 0 + dt;
T_MAX = 40000 * .005;
tvec = T_ZERO:dt:T_MAX;                                                                

Iapp = zeros(size(tvec));
current = randn(1, 40000) * .5e-9;          % Creates a vector with randomly generated currents ranging from -0.5e-9 to 0.5e-9 amperes

spikes = zeros(size(tvec));
potential = zeros(size(tvec));                                          
potential(1) = LEAK_POTENTIAL;
sra = zeros(size(tvec));

j = 1;                                      % Used to iterate through variable currents

for i = 1:length(tvec)-1
   if potential(i) > THRESHOLD              % Tests if current value potential value is greater than threshold
       potential(i) = RESET;                % If so, resets potential to -0.080 V
       sra(i) = sra(i) + b;                 % Increases spike-dependent adaptive current
       spikes(i) = 1;                       % And indicates spike has occured by adding 1 to spike vector at the corresponding time
   end
   dvdt = (LEAK_CONDUCTANCE * (LEAK_POTENTIAL - potential(i)...
       + THRESHOLD_CHANGE * exp((potential(i)- THRESHOLD)/THRESHOLD_CHANGE))... 
       - sra(i) + Iapp(i))/MEMBRANE_CAPACITANCE;
   potential(i+1) = potential(i) + dt * dvdt;                                   % Calculates potential at next time value
   dSRAdt = ((a *(potential(i) - LEAK_POTENTIAL) - sra(i))/SRA_TAU);
   sra(i+1) = sra(i) + dt * dSRAdt;                                             % Calculates adaptive current at next time value
   time = tvec(i+1);                                                            % Records next time value
   temp = mod(time, .005);                                                      % Tests if it is a multiple of 5 ms
   if (temp == 0)
       Iapp(i + 1) = current(j + 1);                                            % If so, sets next applied curret to the next current from random current vector 
       j = j + 1;                                                               % And increases j by one
   else
       Iapp(i + 1) = current(j);                                                % If not, sets the next applied current to the previously applied current from the random current vector
   end
end

figure;

plot(tvec, potential)
xlabel('Time')
ylabel('Potential')
title('Potential over Time')

figure;

plot(tvec, spikes)
xlabel('Time')
ylabel('Spikes')
title('Spikes over Time')
ylim([-3 3]);

%
new_dt = .001;                                   % Sets new time step 
Iapp = expandbin(Iapp, dt, new_dt);              % Down-samples current applied current vector based on new time step
spikes = expandbin(spikes, dt, new_dt);          % Down-samples spikes vector

%

[sta, tcorr] = STA(Iapp, spikes, new_dt);        % Produces spike-triggered average from the two down-sampled vectors over a particular time window

tcorr = -1 .* tcorr;                             % Negates values in the time window in order that the x-axis goes from negative to positive values
figure;

plot(tcorr, sta)
xlabel('Time Window')
ylabel('Spike-triggered Average')
title('STA against Time Window')
