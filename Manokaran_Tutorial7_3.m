%% Rule 1

tau = 20e-3;
r_max = 100;
I_th = 50;
I_sigma = 5;

tmax = 0.5;
dt = 0.001;                                                                                                                                                                                                                                                                                                                                                                                                                  
tvec = 0:dt:tmax;
sigma = 1;

numTrials = 800;

P = [0.05 0.95; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.95 0.05];

r_d = zeros(2,length(tvec));
r_i = zeros(5,length(tvec));

I_app = zeros(5,numTrials);
I_i = zeros(size(r_d));

W_i = zeros(5,2,800);
W_i(:,:,1) = 0.2;
W_d = [0.5 -0.5; -0.5 0.5];
threshold = 40;

epsilon = 0.04;
decision = 0;
correct = zeros(1,800);

for i = 2:numTrials
    weather = round(rand(1)+ 1);
    I_rnd1 = randn(size(tvec));
    I_rnd2 = randn(size(tvec));
    cueUnit1 = (rand(1)<P(1,weather))*50;
    cueUnit2 = (rand(1)<P(2,weather))*50;
    cueUnit3 = (rand(1)<P(3,weather))*50;
    cueUnit4 = (rand(1)<P(4,weather))*50;
    cueUnit5 = (rand(1)<P(5,weather))*50;
    
    for t = 2:length(tvec)
        if (tvec(t) >= 0.1)
            I_app(1,i) = cueUnit1;
            I_app(2,i) = cueUnit2;
            I_app(3,i) = cueUnit3;
            I_app(4,i) = cueUnit4;
            I_app(5,i) = cueUnit5;
        end
            
        for k = 1:5
            dridt = (-r_i(k,t-1) + r_max / (1 + exp((I_th - I_app(k,i))/I_sigma)))/tau;
            r_i(k,t) = r_i(k,(t-1)) + dt*dridt;    
        end
        
        for k = 1:2
            drddt = (-r_d(k,t-1) + r_max / (1 + exp((I_th - I_i(k,t-1))/I_sigma)))/tau;
            r_d(k,t) = r_d(k,(t-1)) + dt*drddt;
        end
        
        I_i(1,t) = W_i(1,1,i-1)*r_i(1,t-1) + W_i(2,1,i-1)*r_i(2,t-1) + W_i(3,1,i-1)*r_i(3,t-1) + ...
            W_i(4,1,i-1)*r_i(4,t-1) + W_i(5,1,i-1)*r_i(5,t-1) + W_d(1,1) * r_d(1,t-1) + W_d(2,1) *...
            r_d(2,t-1) + (I_rnd1(t-1) * sigma/sqrt(dt));
        I_i(2,t) = W_i(1,2,i-1)*r_i(1,t-1) + W_i(2,2,i-1)*r_i(2,t-1) + W_i(3,2,i-1)*r_i(3,t-1) + ...
            W_i(4,2,i-1)*r_i(4,t-1) + W_i(5,2,i-1)*r_i(5,t-1) + W_d(1,2) * r_d(1,t-1) + W_d(2,2) *...
            r_d(2,t-1) + (I_rnd2(t-1) * sigma/sqrt(dt)); 
    end
    
    if (r_d(1,end) > threshold)
        decision = 1;
    elseif (r_d(2,end) > threshold)
       decision = 2; 
    end
    
    if (decision == weather)
        correct(i) = 1;
        E = 0.5;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    elseif (decision ~= weather && decision ~=0)
        E = -0.5;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    else
        for k = 1:5
            W_i(k,1,i) = W_i(k,1,i-1);
            W_i(k,2,i) = W_i(k,2,i-1);
        end
    end
    
    for k = 1:5
       if (W_i(k,1,i)<0)
          W_i(k,1,i)=0; 
       end
       if (W_i(k,2,i)<0)
          W_i(k,2,i)=0; 
       end
       
       if (W_i(k,1,i) > 2*W_i(k,1,1))
           W_i(k,1,i) = 2*W_i(k,1,1);
       end
       
       if (W_i(k,2,i) > 2*W_i(k,2,1))
           W_i(k,2,i) = 2*W_i(k,2,1);
       end
    end
    
end


total_correct = sum(correct);
total_100 = sum(correct(701:800));
W_i1 = W_i(:,1,:);
W_i1 = reshape(W_i1, [5 800]);
W_i2 = W_i(:,2,:);
W_i2 = reshape(W_i2, [5 800]);

figure 
plot(log(P(:,2)/P(:,1)),W_i2(:,end) - W_i1(:,end))
xlabel('Likelihood Ratio')
ylabel('Synaptic Strength Difference')
title('Synaptic Strength Difference vs Likelihood Ratio')

figure
plot(1:800, W_i1)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 1 - Weights Decision Unit 1')

figure
plot(1:800, W_i2)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 1 - Weights Decision Unit 2')

%% Rule 2

tau = 20e-3;
r_max = 100;
I_th = 50;
I_sigma = 5;

tmax = 0.5;
dt = 0.001;
tvec = 0:dt:tmax;
sigma = 1;

numTrials = 800;

P = [0.05 0.95; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.95 0.05];

r_d = zeros(2,length(tvec));
r_i = zeros(5,length(tvec));

I_app = zeros(5,numTrials);
I_i = zeros(size(r_d));

W_i = zeros(5,2,800);
W_i(:,:,1) = 0.2;
W_d = [0.5 -0.5; -0.5 0.5];
threshold = 40;

epsilon = 0.04;
decision = 0;
correct = zeros(1,800);

mean_R = 0.5;
R = zeros(1,800);

for i = 2:numTrials
    weather = round(rand(1)+ 1);
    I_rnd1 = randn(size(tvec));
    I_rnd2 = randn(size(tvec));
    
    cueUnit1 = (rand(1)<P(1,weather))*50;
    cueUnit2 = (rand(1)<P(2,weather))*50;
    cueUnit3 = (rand(1)<P(3,weather))*50;
    cueUnit4 = (rand(1)<P(4,weather))*50;
    cueUnit5 = (rand(1)<P(5,weather))*50;
    
    for t = 2:length(tvec)
        if (tvec(t) >= 0.1)
            I_app(1,i) = cueUnit1;
            I_app(2,i) = cueUnit2;
            I_app(3,i) = cueUnit3;
            I_app(4,i) = cueUnit4;
            I_app(5,i) = cueUnit5;
        end
        for k = 1:5
            dridt = (-r_i(k,t-1) + r_max / (1 + exp((I_th - I_app(k,i))/I_sigma)))/tau;
            r_i(k,t) = r_i(k,(t-1)) + dt*dridt;  
        end
        
        for k = 1:2
            drddt = (-r_d(k,t-1) + r_max / (1 + exp((I_th - I_i(k,t-1))/I_sigma)))/tau;
            r_d(k,t) = r_d(k,(t-1)) + dt*drddt; 
        end
        
        I_i(1,t) = W_i(1,1,i-1)*r_i(1,t-1) + W_i(2,1,i-1)*r_i(2,t-1) + W_i(3,1,i-1)*r_i(3,t-1) + ...
            W_i(4,1,i-1)*r_i(4,t-1) + W_i(5,1,i-1)*r_i(5,t-1) + W_d(1,1) * r_d(1,t-1) + W_d(2,1) *...
            r_d(2,t-1) + (I_rnd1(t-1) * sigma/sqrt(dt));
        
        I_i(2,t) = W_i(1,2,i-1)*r_i(1,t-1) + W_i(2,2,i-1)*r_i(2,t-1) + W_i(3,2,i-1)*r_i(3,t-1) + ...
            W_i(4,2,i-1)*r_i(4,t-1) + W_i(5,2,i-1)*r_i(5,t-1) + W_d(1,2) * r_d(1,t-1) + W_d(2,2) *...
            r_d(2,t-1) + (I_rnd2(t-1) * sigma/sqrt(dt)); 
    end
    
    if (i>11)
       mean_R = mean(R(i-10:i)); 
        
    end
    
    if (r_d(1,end) > threshold)
        decision = 1;
    elseif (r_d(2,end) > threshold)
       decision = 2; 
    end

    if (decision == weather) 
        correct(i) = 1;
        R(i) = 1;
        E = R(i) - mean_R;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    elseif (decision ~= weather && decision ~=0) 
        R(i) = 0;
        E = R(i) - mean_R;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    else
        for k = 1:5
            W_i(k,1,i) = W_i(k,1,i-1);
            W_i(k,2,i) = W_i(k,2,i-1);
        end
    end
    
    for k = 1:5
       if (W_i(k,1,i)<0)
          W_i(k,1,i)=0; 
       end
       if (W_i(k,2,i)<0)
          W_i(k,2,i)=0; 
       end
       
       if (W_i(k,1,i) > 2*W_i(k,1,1))
           W_i(k,1,i) = 2*W_i(k,1,1);
       end
       
       if (W_i(k,2,i) > 2*W_i(k,2,1))
           W_i(k,2,i) = 2*W_i(k,2,1);
       end
    end
    
end

total_correct = sum(correct);
total_100 = sum(correct(701:800));
W_i1 = W_i(:,1,:);
W_i1 = reshape(W_i1, [5 800]);
W_i2 = W_i(:,2,:);
W_i2 = reshape(W_i2, [5 800]);

figure 
plot(log(P(:,2)/P(:,1)),W_i2(:,end) - W_i1(:,end))
xlabel('Likelihood Ratio')
ylabel('Synaptic Strength Difference')
title('Synaptic Strength Difference vs Likelihood Ratio')

figure
plot(1:800, W_i1)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 2 - Weights Decision Unit 1')

figure
plot(1:800, W_i2)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 2 - Weights Decision Unit 2')

%% Rule 4

tau = 20e-3;
r_max = 100;
I_th = 50;
I_sigma = 5;

tmax = 0.5;
dt = 0.001;
tvec = 0:dt:tmax;
sigma = 1;

numTrials = 800;

P = [0.05 0.95; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.95 0.05];

r_d = zeros(2, length(tvec));
r_i = zeros(5,length(tvec));

I_app = zeros(5,numTrials);
I_i = zeros(size(r_d));

W_i = zeros(5,2,800);
W_i(:,:,1) = 0.2;
W_d = [0.5 -0.5; -0.5 0.5];
threshold = 40;

epsilon = 0.04;
decision = 0;
correct = zeros(1,800);

mean_R = 0.5;
R = zeros(1,800);

for i = 2:numTrials
    weather = round(rand(1)+ 1);
    I_rnd1 = randn(size(tvec));
    I_rnd2 = randn(size(tvec));
    cueUnit1 = (rand(1)<P(1,weather))*50;
    cueUnit2 = (rand(1)<P(2,weather))*50;
    cueUnit3 = (rand(1)<P(3,weather))*50;
    cueUnit4 = (rand(1)<P(4,weather))*50;
    cueUnit5 = (rand(1)<P(5,weather))*50;
    
    for t = 2:length(tvec)
        if (tvec(t) >= 0.1)
            I_app(1,i) = cueUnit1;
            I_app(2,i) = cueUnit2;
            I_app(3,i) = cueUnit3;
            I_app(4,i) = cueUnit4;
            I_app(5,i) = cueUnit5;
        end
        for k = 1:5
            dridt = (-r_i(k,t-1) + r_max / (1 + exp((I_th - I_app(k,i))/I_sigma)))/tau;
            r_i(k,t) = r_i(k,(t-1)) + dt*dridt; 
        end
        
        for k = 1:2
            drddt = (-r_d(k,t-1) + r_max / (1 + exp((I_th - I_i(k,t-1))/I_sigma)))/tau;
            r_d(k,t) = r_d(k,(t-1)) + dt*drddt;
        end
        
        I_i(1,t) = W_i(1,1,i-1)*r_i(1,t-1) + W_i(2,1,i-1)*r_i(2,t-1) + W_i(3,1,i-1)*r_i(3,t-1) + ...
            W_i(4,1,i-1)*r_i(4,t-1) + W_i(5,1,i-1)*r_i(5,t-1) + W_d(1,1) * r_d(1,t-1) + W_d(2,1) *...
            r_d(2,t-1) + (I_rnd1(t-1) * sigma/sqrt(dt));
        
        I_i(2,t) = W_i(1,2,i-1)*r_i(1,t-1) + W_i(2,2,i-1)*r_i(2,t-1) + W_i(3,2,i-1)*r_i(3,t-1) + ...
            W_i(4,2,i-1)*r_i(4,t-1) + W_i(5,2,i-1)*r_i(5,t-1) + W_d(1,2) * r_d(1,t-1) + W_d(2,2) *...
            r_d(2,t-1) + (I_rnd2(t-1) * sigma/sqrt(dt));  
    end
    
    if (i>11)
       mean_R = mean(R(i-10:i)); 
        
    end
    
    if (r_d(1,end) > threshold)
        decision = 1;
    elseif (r_d(2,end) > threshold) 
       decision = 2; 
    end
    
    if (decision == weather) % correct response to cues
        correct(i) = 1;
        R(i) = 1;
        E = R(i) - mean_R;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    else
        for k = 1:5
            W_i(k,1,i) = W_i(k,1,i-1);
            W_i(k,2,i) = W_i(k,2,i-1);
        end
    end
    
    for k = 1:5
       if (W_i(k,1,i)<0)
          W_i(k,1,i)=0; 
       end
       if (W_i(k,2,i)<0)
          W_i(k,2,i)=0; 
       end
       
       if (W_i(k,1,i) > 2*W_i(k,1,1))
           W_i(k,1,i) = 2*W_i(k,1,1);
       end
       
       if (W_i(k,2,i) > 2*W_i(k,2,1))
           W_i(k,2,i) = 2*W_i(k,2,1);
       end
    end
    
end

total_correct = sum(correct);
total_100 = sum(correct(701:800));
W_i1 = W_i(:,1,:);
W_i1 = reshape(W_i1, [5 800]);
W_i2 = W_i(:,2,:);
W_i2 = reshape(W_i2, [5 800]);

figure 
plot(log(P(:,2)/P(:,1)),W_i2(:,end) - W_i1(:,end))
xlabel('Likelihood Ratio')
ylabel('Synaptic Strength Difference')
title('Synaptic Strength Difference vs Likelihood Ratio')

figure
plot(1:800, W_i1)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 4 - Weights Decision Unit 1')

figure
plot(1:800, W_i2)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 4 - Weights Decision Unit 2')

%% Alternative Training Method Part 1

tau = 20e-3;
r_max = 100;
I_th = 50;
I_sigma = 5;

tmax = 0.5;
dt = 0.001;
tvec = 0:dt:tmax;
sigma = 1;

numTrials = 800;

P = [0.05 0.95; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.95 0.05];

r_d = zeros(2, length(tvec));
r_i = zeros(5,length(tvec));

I_app = zeros(5,numTrials);
I_i = zeros(size(r_d));

W_i = zeros(5,2,800);
W_i(:,:,1) = 0.2;
W_d = [0.5 -0.5; -0.5 0.5];
threshold = 40;

epsilon = 0.04;
decision = 0;
correct = zeros(1,800);

for i = 2:numTrials
    I_rnd1 = randn(size(tvec));
    I_rnd2 = randn(size(tvec));
    if (rand<0.2)
       weather = 2; 
    else
        weather = 1;
    end
    
    cueUnit1 = (rand(1)<P(1,weather))*50;
    cueUnit2 = (rand(1)<P(2,weather))*50;
    cueUnit3 = (rand(1)<P(3,weather))*50;
    cueUnit4 = (rand(1)<P(4,weather))*50;
    cueUnit5 = (rand(1)<P(5,weather))*50;
    
    for t = 2:length(tvec)
        if (tvec(t) >= 0.1)
            I_app(1,i) = cueUnit1;
            I_app(2,i) = cueUnit2;
            I_app(3,i) = cueUnit3;
            I_app(4,i) = cueUnit4;
            I_app(5,i) = cueUnit5;
        end
        for k = 1:5
            dridt = (-r_i(k,t-1) + r_max / (1 + exp((I_th - I_app(k,i))/I_sigma)))/tau;
            r_i(k,t) = r_i(k,(t-1)) + dt*dridt;     
        end
        
        for k = 1:2
            drddt = (-r_d(k,t-1) + r_max / (1 + exp((I_th - I_i(k,t-1))/I_sigma)))/tau;
            r_d(k,t) = r_d(k,(t-1)) + dt*drddt;
        end

        I_i(1,t) = W_i(1,1,i-1)*r_i(1,t-1) + W_i(2,1,i-1)*r_i(2,t-1) + W_i(3,1,i-1)*r_i(3,t-1) + ...
            W_i(4,1,i-1)*r_i(4,t-1) + W_i(5,1,i-1)*r_i(5,t-1) + W_d(1,1) * r_d(1,t-1) + W_d(2,1) *...
            r_d(2,t-1) + (I_rnd1(t-1) * sigma/sqrt(dt));
        
        I_i(2,t) = W_i(1,2,i-1)*r_i(1,t-1) + W_i(2,2,i-1)*r_i(2,t-1) + W_i(3,2,i-1)*r_i(3,t-1) + ...
            W_i(4,2,i-1)*r_i(4,t-1) + W_i(5,2,i-1)*r_i(5,t-1) + W_d(1,2) * r_d(1,t-1) + W_d(2,2) *...
            r_d(2,t-1) + (I_rnd2(t-1) * sigma/sqrt(dt)); 
        
        
    end
    
    if (r_d(1,end) > threshold)
        decision = 1;
    elseif (r_d(2,end) > threshold)
       decision = 2; 
    end

    if (decision == weather)
        correct(i) = 1;
        E = 0.5;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    elseif (decision ~= weather && decision ~=0)
        E = -0.5;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    else
        for k = 1:5
            W_i(k,1,i) = W_i(k,1,i-1);
            W_i(k,2,i) = W_i(k,2,i-1);
        end
    end
    
    for k = 1:5
       if (W_i(k,1,i)<0)
          W_i(k,1,i)=0; 
       end
       if (W_i(k,2,i)<0)
          W_i(k,2,i)=0; 
       end
       
       if (W_i(k,1,i) > 2*W_i(k,1,1))
           W_i(k,1,i) = 2*W_i(k,1,1);
       end
       
       if (W_i(k,2,i) > 2*W_i(k,2,1))
           W_i(k,2,i) = 2*W_i(k,2,1);
       end
    end
    
end

total_correct = sum(correct);
total_100 = sum(correct(701:800));
W_i1 = W_i(:,1,:);
W_i1 = reshape(W_i1, [5 800]);
W_i2 = W_i(:,2,:);
W_i2 = reshape(W_i2, [5 800]);

figure 
plot(log(P(:,2)/P(:,1)),W_i2(:,end) - W_i1(:,end))
xlabel('Likelihood Ratio')
ylabel('Synaptic Strength Difference')
title('Synaptic Strength Difference vs Likelihood Ratio')

figure
plot(1:800, W_i1)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 1 - Weights Decision Unit 1')

figure
plot(1:800, W_i2)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 1 - Weights Decision Unit 2')

%% Alternative Training Method Part 2

tau = 20e-3;
r_max = 100;
I_th = 50;
I_sigma = 5;

tmax = 0.5;
dt = 0.001;
tvec = 0:dt:tmax;
sigma = 1;

numTrials = 4000;

P = [0.05 0.95; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.95 0.05];

r_d = zeros(2, length(tvec));
r_i = zeros(5,length(tvec));

I_app = zeros(5,numTrials);
I_i = zeros(size(r_d));

W_i = zeros(5,2,4000);
W_i(:,:,1) = 0.2;
W_d = [0.5 -0.5; -0.5 0.5];
threshold = 40;

epsilon = 0.04;
decision = 0;
correct = zeros(1,4000);

for i = 2:numTrials
    weather = round(rand(1)+ 1);
    I_rnd1 = randn(size(tvec));
    I_rnd2 = randn(size(tvec));
    
    cueUnit1 = (rand(1)<P(1,weather))*50;
    cueUnit2 = (rand(1)<P(2,weather))*50;
    cueUnit3 = (rand(1)<P(3,weather))*50;
    cueUnit4 = (rand(1)<P(4,weather))*50;
    cueUnit5 = (rand(1)<P(5,weather))*50;
    
    for t = 2:length(tvec)
        if (tvec(t) >= 0.1)
            I_app(1,i) = cueUnit1;
            I_app(2,i) = cueUnit2;
            I_app(3,i) = cueUnit3;
            I_app(4,i) = cueUnit4;
            I_app(5,i) = cueUnit5;
        end
        
        for k = 1:5
            dridt = (-r_i(k,t-1) + r_max / (1 + exp((I_th - I_app(k,i))/I_sigma)))/tau;
            r_i(k,t) = r_i(k,(t-1)) + dt*dridt;    
        end
        
        for k = 1:2
            drddt = (-r_d(k,t-1) + r_max / (1 + exp((I_th - I_i(k,t-1))/I_sigma)))/tau;
            r_d(k,t) = r_d(k,(t-1)) + dt*drddt;
        end
        
        I_i(1,t) = W_i(1,1,i-1)*r_i(1,t-1) + W_i(2,1,i-1)*r_i(2,t-1) + W_i(3,1,i-1)*r_i(3,t-1) + ...
            W_i(4,1,i-1)*r_i(4,t-1) + W_i(5,1,i-1)*r_i(5,t-1) + W_d(1,1) * r_d(1,t-1) + W_d(2,1) *...
            r_d(2,t-1) + (I_rnd1(t-1) * sigma/sqrt(dt));
        
        I_i(2,t) = W_i(1,2,i-1)*r_i(1,t-1) + W_i(2,2,i-1)*r_i(2,t-1) + W_i(3,2,i-1)*r_i(3,t-1) + ...
            W_i(4,2,i-1)*r_i(4,t-1) + W_i(5,2,i-1)*r_i(5,t-1) + W_d(1,2) * r_d(1,t-1) + W_d(2,2) *...
            r_d(2,t-1) + (I_rnd2(t-1) * sigma/sqrt(dt)); 
    end
    
    if (r_d(1,end) > threshold)
        decision = 1;
    

    elseif (r_d(2,end) > threshold)
       decision = 2; 
    end

    if (decision == weather) 
        correct(i) = 1;
        E = 0.5;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    elseif (decision ~= weather && decision ~=0)
        E = -0.5;
        for k = 1:5
            if (I_app(k,i) > 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) + epsilon * E;
            elseif (I_app(k,i) == 0)
                W_i(k,decision,i) = W_i(k,decision,i-1) - epsilon * E;
            end
        end
    else
        for k = 1:5
            W_i(k,1,i) = W_i(k,1,i-1);
            W_i(k,2,i) = W_i(k,2,i-1);
        end
    end
    
    for k = 1:5
       if (W_i(k,1,i)<0)
          W_i(k,1,i)=0; 
       end
       if (W_i(k,2,i)<0)
          W_i(k,2,i)=0; 
       end
       
       if (W_i(k,1,i) > 2*W_i(k,1,1))
           W_i(k,1,i) = 2*W_i(k,1,1);
       end
       
       if (W_i(k,2,i) > 2*W_i(k,2,1))
           W_i(k,2,i) = 2*W_i(k,2,1);
       end
    end
    
end

total_correct = sum(correct);
total_100 = sum(correct(3901:4000));
W_i1 = W_i(:,1,:);
W_i1 = reshape(W_i1, [5 4000]);
W_i2 = W_i(:,2,:);
W_i2 = reshape(W_i2, [5 4000]);

figure 
plot(log(P(:,2)/P(:,1)),W_i2(:,end) - W_i1(:,end))
xlabel('Likelihood Ratio')
ylabel('Synaptic Strength Difference')
title('Synaptic Strength Difference vs Likelihood Ratio')

figure
plot(1:4000, W_i1)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 1 - Weights Decision Unit 1')

figure
plot(1:4000, W_i2)
xlabel('Trial Number')
ylabel('Synaptic Weight')
legend('Cue A', 'Cue B','Cue C','Cue D','Cue E')
title('Rule 1 - Weights Decision Unit 2')
