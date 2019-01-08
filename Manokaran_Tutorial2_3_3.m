clear

LEAK_POTENTIAL = -.07;
LEAK_CONDUCTANCE = 10e-9;
MEMBRANE_CAPACITANCE = 100e-12;

a = 2e-9; % nS
b = 0; % nA
SRA_TAU = .150;

RESET = -.080;
THRESHOLD = -.050;
V_MAX = .050;
THRESHOLD_CHANGE = .002;

dt = .00001;
T_ZERO = 0 + dt;
T_MAX = .5;
tvec = T_ZERO:dt:T_MAX;                                                                
       
sigma = 50e-12;
spikes = zeros(1, 1000);

for j = 1:length(spikes)
    potential = zeros(size(tvec));                                          
    potential(1) = LEAK_POTENTIAL;
    sra = zeros(size(tvec));
    Iapp = randn(size(tvec)) * sigma/sqrt(dt);
    for i = 1:length(tvec)-1
       if potential(i) > THRESHOLD              
           potential(i) = RESET;                
           sra(i) = sra(i) + b;                 
           spikes(j) = spikes(j) + 1;                       
       end
       dvdt = (LEAK_CONDUCTANCE * (LEAK_POTENTIAL - potential(i)...
           + THRESHOLD_CHANGE * exp((potential(i)- THRESHOLD)/THRESHOLD_CHANGE))... 
           - sra(i) + Iapp(i))/MEMBRANE_CAPACITANCE;
       potential(i+1) = potential(i) + dt * dvdt;                                   
       dSRAdt = ((a *(potential(i) - LEAK_POTENTIAL) - sra(i))/SRA_TAU);
       sra(i+1) = sra(i) + dt * dSRAdt;                                             
    end
end

mean = mean(spikes);
display(mean/T_MAX)

range1 = 0:max(spikes);
[count1, centers1] = hist(spikes, range1);
fraction = cumsum(hist(spikes, range1)) ./ 1000;
probablity1 = 1 - fraction;

%

LEAK_POTENTIAL = -.07;
LEAK_CONDUCTANCE = 10e-9;
MEMBRANE_CAPACITANCE = 100e-12;

a = 2e-9; % nS
b = 0; % nA
SRA_TAU = .150;

RESET = -.080;
THRESHOLD = -.050;
V_MAX = .050;
THRESHOLD_CHANGE = .002;

dt = .00001;
T_ZERO = 0;
T_MAX = .5;
tvec = T_ZERO:dt:T_MAX;
current2 = (randn(size(tvec)) * .2e-9) + .1e-9; 

sigma = 5e-12;

spikes2 = zeros(1, 1000);
for j = 1:length(spikes2) 
    potential = zeros(size(tvec));                                          
    potential(1) = LEAK_POTENTIAL;
    sra = zeros(size(tvec));
    Iapp2 = (randn(size(tvec)) * sigma/sqrt(dt)) + .5e-9;
    for i = 1:length(tvec)-1
       if potential(i) > THRESHOLD              
           potential(i) = RESET;                
           sra(i) = sra(i) + b;                 
           spikes2(j) = spikes2(j) + 1;                      
       end
       dvdt = (LEAK_CONDUCTANCE * (LEAK_POTENTIAL - potential(i)...
           + THRESHOLD_CHANGE * exp((potential(i)- THRESHOLD)/THRESHOLD_CHANGE))... 
           - sra(i) + Iapp2(i))/MEMBRANE_CAPACITANCE;
       potential(i+1) = potential(i) + dt * dvdt;                                   
       dSRAdt = ((a *(potential(i) - LEAK_POTENTIAL) - sra(i))/SRA_TAU);
       sra(i+1) = sra(i) + dt * dSRAdt;                                             
       Iapp2(i+1) = Iapp2(i+1);
    end
end
clear mean

mean = mean(spikes2);
display(mean/T_MAX)

range2 = 0:max(spikes2);
[count2, centers2] = hist(spikes2, range2);
fraction2 = cumsum(hist(spikes2, range2)) ./ 1000;
probablity2 = 1 - fraction2;

figure;

stairs(centers1, count1)
hold on
stairs(centers2, count2)
title('Number of Trials against Spike Count')
xlabel('Spike Count')
ylabel('Number of Trials');
legend('- stimulus', '+ stimulus')

figure;

stairs(centers1, probablity1)
hold on
stairs(centers2, probablity2)
title('Probability against Spike Count')
xlabel('Spike Count')
ylabel('Probability')
legend('- stimulus', '+ stimulus')

if (length(probablity2) > length(probablity1))
    temp = probablity1;
    newlength = length(probablity2) - length(probablity1);
    zerosvec = zeros(1, newlength);
    probablity1 = horzcat(temp, zerosvec);
end

figure;

plot(probablity1, probablity2)
title('ROC Curve')
xlabel('False +ve')
ylabel('True +ve')