LEAK_POTENTIAL = -.07;
LEAK_CONDUCTANCE = 10e-9;
POTASSIUM_POTENTIAL = -.08;
MEMBRANE_CAPACITANCE = 100e-12;

a = 2e-9; % nS
b = 0; % nA
SRA_TAU = .05;

RESET = -.080;
THRESHOLD = -.050;
V_MAX = .050;
THRESHOLD_CHANGE = .002;

dt = .00001;
T_ZERO = 0 + dt;
T_MAX = 100;
tvec = T_ZERO:dt:T_MAX; 

sigma = 50e-12;
noise = sigma / sqrt(dt);

% Iapp = zeros(size(tvec));
Iapp = randn(size(tvec)) * noise;          % Creates a vector with randomly generated currents ranging from -0.5e-9 to 0.5e-9 amperes

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
end

spiketimes = find(spikes) * dt;
ISIs = diff(spiketimes);
totaltime = sum(ISIs);
mean = totaltime / length(ISIs);
deviation = std(ISIs);
C_V = deviation / mean;

figure;
hist(ISIs)

figure;
hist(ISIs, 25)

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

LEAK_POTENTIAL = -.07;
LEAK_CONDUCTANCE = 10e-9;
POTASSIUM_POTENTIAL = -.08;
MEMBRANE_CAPACITANCE = 100e-12;

a = 2e-9; % nS
b = 0; % nA
SRA_TAU = .05;

RESET = -.080;
THRESHOLD = -.050;
V_MAX = .050;
THRESHOLD_CHANGE = .002;

dt = .00001;
T_ZERO = 0 + dt;
T_MAX = 100;
tvec = T_ZERO:dt:T_MAX; 

sigma = 50e-12;
noise = sigma / sqrt(dt);

% Iapp = zeros(size(tvec));
Iapp = randn(size(tvec)) * noise;          % Creates a vector with randomly generated currents ranging from -0.5e-9 to 0.5e-9 amperes

spikes = zeros(size(tvec));
potential = zeros(size(tvec));                                          
potential(1) = LEAK_POTENTIAL;
sra = zeros(size(tvec));
spikes_100 = zeros(1, 1000);
counter = 0;
j = 1;                                      % Used to iterate through variable currents

for i = 1:length(tvec)-1
   if potential(i) > THRESHOLD              % Tests if current value potential value is greater than threshold
       potential(i) = RESET;                % If so, resets potential to -0.080 V
       sra(i) = sra(i) + b;                 % Increases spike-dependent adaptive current
       spikes(i) = 1;                       % And indicates spike has occured by adding 1 to spike vector at the corresponding time
       counter = counter + 1;
   end
   dvdt = (LEAK_CONDUCTANCE * (LEAK_POTENTIAL - potential(i)...
       + THRESHOLD_CHANGE * exp((potential(i)- THRESHOLD)/THRESHOLD_CHANGE))... 
       - sra(i) + Iapp(i))/MEMBRANE_CAPACITANCE;
   potential(i+1) = potential(i) + dt * dvdt;                                   % Calculates potential at next time value
   dSRAdt = ((a *(potential(i) - LEAK_POTENTIAL) - sra(i))/SRA_TAU);
   sra(i+1) = sra(i) + dt * dSRAdt;                                             % Calculates adaptive current at next time value
   time = tvec(i+1);                                                            % Records next time value
   temp = mod(time, .1);
   if (time == .1 * j)
       spikes_100(j) = counter;
       counter = 0;
       j = j + 1;
   end
end

spiketimes = find(spikes_100) * dt * 10000;
ISIs = diff(spiketimes);
totaltime = sum(ISIs);
mean_100 = totaltime / length(ISIs);
display(mean_100)
deviation_100 = std(spikes_100);
display(deviation_100)
fano_factor = deviation_100 / mean_100;


% figure;
% hist(ISIs)
% 
% figure;
% hist(ISIs, 25)
% 
% figure;
% 
% plot(tvec, potential)
% xlabel('Time')
% ylabel('Potential')
% title('Potential over Time')
% 
% figure;
% 
% plot(tvec, spikes)
% xlabel('Time')
% ylabel('Spikes')
% title('Spikes over Time')
% ylim([-3 3]);

%

LEAK_POTENTIAL = -.07;
LEAK_CONDUCTANCE = 10e-9;
POTASSIUM_POTENTIAL = -.08;
MEMBRANE_CAPACITANCE = 100e-12;

a = 2e-9; % nS
b = 1e-9; % nA
SRA_TAU = .05;

RESET = -.080;
THRESHOLD = -.050;
V_MAX = .050;
THRESHOLD_CHANGE = .002;

dt = .00001;
T_ZERO = 0 + dt;
T_MAX = 100;
tvec = T_ZERO:dt:T_MAX; 

sigma = 50e-12;
noise = sigma / sqrt(dt);

% Iapp = zeros(size(tvec));
Iapp = randn(size(tvec)) * noise;          % Creates a vector with randomly generated currents ranging from -0.5e-9 to 0.5e-9 amperes

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
end

spiketimes = find(spikes) * dt;
ISIs = diff(spiketimes);
totaltime = sum(ISIs);
mean = totaltime / length(ISIs);
deviation = std(ISIs);
C_V = deviation / mean;
display(C_V)

figure;
hist(ISIs)

figure;
hist(ISIs, 25)

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

