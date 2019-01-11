A = 20;
B = 10;
f = 0.5;

dt = 0.001;
tvec = 0:dt:10;

I_a = A * sin(2*pi*f*tvec);
I_b = B * cos(2*pi*f*tvec);

numNeurons = 50;

I_0 = 50;

sigma = 10;

firing_rate = zeros(numNeurons, length(tvec));

eta = randn(size(firing_rate));

W_0 = randn(numNeurons, 1);
W_A = randn(numNeurons, 1);
W_B = randn(numNeurons, 1);

for t = 1:length(tvec)
  I_A = I_a(t);
  I_B = I_b(t);
  
  firing_rate(:,t) = 100 + W_0 * I_0 + W_A * I_A + W_B * I_B + sigma * eta(:,t);
end

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(firing_rate');

figure
subplot(1, 2, 1)
plot(W_A, COEFF(:,1))
title('First PC against W_A')
subplot(1, 2, 2)
plot(W_B, COEFF(:,2))
title('Second PC against W_B')

fractionExplained = EXPLAINED(1) + EXPLAINED(2);
figure
plot(EXPLAINED)
title(['Fraction Explained = ', num2str(fractionExplained)]);
transpose_SCORE = SCORE';

denoised_firing_rate = COEFF(:,1:2) * transpose_SCORE(1:2,:);
for i = 1:size(denoised_firing_rate,2)
    denoised_firing_rate(:,i) = denoised_firing_rate(:,i) + MU(1,:)';
end

figure
plot(tvec,firing_rate(1:2,:))
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate')

figure
plot(tvec,denoised_firing_rate(1:2,:))
xlabel('Time')
ylabel('Firing Rate')
title('Denoised Firing Rate')

figure
plot(firing_rate(1,:), firing_rate(2,:))
xlabel('Neuron 1')
ylabel('Neuron 2')
title('Firing Rate')

figure
plot(denoised_firing_rate(1,:), denoised_firing_rate(2,:))
xlabel('Neuron 1')
ylabel('Neuron 2')
title('Firing Rate')

figure
plot(tvec,SCORE(:,1))
xlabel('Time')
ylabel('Score')

figure
plot(tvec,SCORE(:,2))
xlabel('Time')
ylabel('Score')

%%

A = 20;
B = 10;
f_a = 1;
f_b = 0.5;

dt = 0.001;
tvec = 0:dt:10;

I_a = A * sin(2*pi*f_a*tvec);
I_b = B * cos(2*pi*f_b*tvec);

numNeurons = 50;

I_0 = 50;

sigma = 10;

firing_rate = zeros(numNeurons, length(tvec));

eta = randn(size(firing_rate));

W_0 = randn(numNeurons, 1);
W_A = randn(numNeurons, 1);
W_B = randn(numNeurons, 1);

for t = 1:length(tvec)
  I_A = I_a(t);
  I_B = I_b(t);
  
  firing_rate(:,t) = 100 + W_0 * I_0 + W_A * I_A + W_B * I_B + sigma * eta(:,t);
end

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(firing_rate');

figure
subplot(1, 2, 1)
plot(W_A, COEFF(:,1))
title('First PC against W_A')
subplot(1, 2, 2)
plot(W_B, COEFF(:,2))
title('Second PC against W_B')

fractionExplained = EXPLAINED(1) + EXPLAINED(2);
figure
plot(EXPLAINED)
title(['Fraction Explained = ', num2str(fractionExplained)]);
transpose_SCORE = SCORE';

denoised_firing_rate = COEFF(:,1:2) * transpose_SCORE(1:2,:);
for i = 1:size(denoised_firing_rate,2)
    denoised_firing_rate(:,i) = denoised_firing_rate(:,i) + MU(1,:)';
end

figure
plot(tvec,firing_rate(1:2,:))
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate')

figure
plot(tvec,denoised_firing_rate(1:2,:))
xlabel('Time')
ylabel('Firing Rate')
title('Denoised Firing Rate')

figure
plot(firing_rate(1,:), firing_rate(2,:))
xlabel('Neuron 1')
ylabel('Neuron 2')
title('Firing Rate')

figure
plot(denoised_firing_rate(1,:), denoised_firing_rate(2,:))
xlabel('Neuron 1')
ylabel('Neuron 2')
title('Firing Rate')

figure
plot(tvec,SCORE(:,1))
xlabel('Time')
ylabel('Score')

figure
plot(tvec,SCORE(:,2))
xlabel('Time')
ylabel('Score')

%%

A = 20;
B = 10;
f = 0.5;

dt = 0.001;
tvec = 0:dt:10;

I_a = A * tvec;
I_b = zeros(size(tvec));
start = round(4/dt);
stop = round(5/dt);
I_b(start:stop) = B * sin(2*pi*f*tvec(start:stop));

numNeurons = 50;

I_0 = 50;

sigma = 10;

firing_rate = zeros(numNeurons, length(tvec));

eta = randn(size(firing_rate));

W_0 = randn(numNeurons, 1);
W_A = randn(numNeurons, 1);
W_B = randn(numNeurons, 1);

for t = 1:length(tvec)
  I_A = I_a(t);
  I_B = I_b(t);
  
  firing_rate(:,t) = 100 + W_0 * I_0 + W_A * I_A + W_B * I_B + sigma * eta(:,t);
end

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(firing_rate');

figure
subplot(1, 2, 1)
plot(W_A, COEFF(:,1))
title('First PC against W_A')
subplot(1, 2, 2)
plot(W_B, COEFF(:,2))
title('Second PC against W_B')

fractionExplained = EXPLAINED(1) + EXPLAINED(2);
figure
plot(EXPLAINED)
title(['Fraction Explained = ', num2str(fractionExplained)]);
transpose_SCORE = SCORE';

denoised_firing_rate = COEFF(:,1:2) * transpose_SCORE(1:2,:);
for i = 1:size(denoised_firing_rate,2)
    denoised_firing_rate(:,i) = denoised_firing_rate(:,i) + MU(1,:)';
end

figure
plot(tvec,firing_rate(1:2,:))
xlabel('Time')
ylabel('Firing Rate')
title('Firing Rate')

figure
plot(tvec,denoised_firing_rate(1:2,:))
xlabel('Time')
ylabel('Firing Rate')
title('Denoised Firing Rate')

figure
plot(firing_rate(1,:), firing_rate(2,:))
xlabel('Neuron 1')
ylabel('Neuron 2')
title('Firing Rate')

figure
plot(denoised_firing_rate(1,:), denoised_firing_rate(2,:))
xlabel('Neuron 1')
ylabel('Neuron 2')
title('Firing Rate')

figure
plot(tvec,SCORE(:,1))
xlabel('Time')
ylabel('Score')

figure
plot(tvec,SCORE(:,2))
xlabel('Time')
ylabel('Score')



