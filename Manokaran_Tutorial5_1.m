%% Part 1
rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.5;
pR = 1;
Wee = 8;
tau_s = .002;
r0 = 0;
x = 1.2;
D = 1;

rrange = 0:1:rmax;
srange = 0:.01:1;

figure;
total_input = Wee*srange;
r_s = r0 + rmax*total_input.^x./(total_input.^x + sigma^x);

subplot(1,2,1)
plot(srange, r_s)
title('Part 1')
hold on

alpha = alpha0*D;
s_r = (alpha*pR*tau_s*rrange)./(1+alpha*pR*tau_s*rrange);

plot(s_r, rrange,'--')
ylabel('Rate (Hz)')
xlabel('Synaptic Input')

%%

rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.5;
D = 1;
pR = 1;
Wee = 8;
tau_s = .002;
dt = .0001;
tvec = 0:dt:20;
S_in = 0;
start = 10;
duration = .050;
stop = 10+duration;
r0 = 0;
rvec = zeros(size(tvec));
svec = zeros(size(tvec));

s_in = zeros(size(tvec));
x = 1.2;

for i = 1:length(tvec)-1
    time = tvec(i);
    if (time >= start && time <= stop)
        s_in(i) = .05;
    end
    
    alpha = alpha0*D;
    total_input = Wee*svec(i) + s_in(i);
    r_ss = r0 + rmax*total_input^x/(total_input^x + sigma^x);
    rvec(i+1) = rvec(i) + dt*(r_ss - rvec(i) )/tau_r;
    if (rvec(i+1) < 0)
        rvec(i+1) = 0;
    end
    dsdt = -svec(i)/tau_s + alpha*pR*rvec(i)*(1-svec(i));
    svec(i+1) = svec(i) + dt*dsdt;
end

subplot(1,2,2)
plot(tvec, rvec)
title('Part 1')
xlabel('Time')
ylabel('Rate (Hz)')

%% Part 2
rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.5;
pR = 0.2;
Wee = 60;
tau_s = .002;
r0 = 0.1;
x = 1.2;
tau_d = .250;

rrange = 0:1:rmax;
srange = 0:.01:1;

D = zeros(size(rrange));

figure;
total_input = Wee*srange;
r_s = r0 + rmax*total_input.^x./(total_input.^x + sigma^x);

subplot(1,2,1)
plot(srange, r_s)
title('Part 2')

hold on
for i = 1:length(rrange)-1
    alpha = alpha0*D(i);
    s_r(i) = (alpha*pR*tau_s*rrange(i))/(1+alpha*pR*tau_s*rrange(i));
    dDdt = (1-D(i))/tau_d - pR*D(i)*rrange(i);
    D(i+1) = D(i) + dt*dDdt;
end  
 
plot(s_r, rrange,'--')
xlabel('Synaptic Input')
ylabel('Rate (Hz)')

%%

rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.5;
D = 1;
pR = .2;
Wee = 60;
tau_s = .002;
dt = .0001;
tvec = 0:dt:20;
S_in = 0;
start = 10;
duration = 2;
stop = 10+duration;
r0 = .1;
rvec = zeros(size(tvec));
svec = zeros(size(tvec));

s_in = zeros(size(tvec));
x = 1.2;
tau_d = .250;

D = zeros(size(tvec));

for i = 1:length(tvec)-1
    time = tvec(i);
    if (time >= start && time <= stop)
        s_in(i) = .002;
    end
    alpha = alpha0*D(i);
    
    dDdt = (1-D(i))/tau_d - pR*D(i)*rvec(i);
    D(i+1) = D(i) + dt*dDdt;
    
    total_input = Wee*svec(i) + s_in(i);
    r_ss = r0 + rmax*total_input^x/(total_input^x + sigma^x);
    rvec(i+1) = rvec(i) + dt*(r_ss - rvec(i) )/tau_r;
    if (rvec(i+1) < 0)
        rvec(i+1) = 0;
    end
    dsdt = -svec(i)/tau_s + alpha*pR*rvec(i)*(1-svec(i));
    svec(i+1) = svec(i) + dt*dsdt;
end

subplot(1,2,2)
plot(tvec, rvec)
title('Part 2')
xlabel('Time')
ylabel('Rate (Hz)')

%% Part 3
rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.5;
pR = 0.5;
Wee = 35;
tau_s = .002;
r0 = -0.1;
x = 1.2;
tau_d = .250;

rrange = 0:1:rmax;
srange = 0:.01:1;

D = zeros(size(rrange));

figure;
total_input = Wee*srange;
r_s = r0 + rmax*total_input.^x./(total_input.^x + sigma^x);

subplot(1,3,1)
plot(srange, r_s)
title('Part 3')

hold on
for i = 1:length(rrange)-1
    alpha = alpha0*D(i);
    s_r(i) = (alpha*pR*tau_s*rrange(i))/(1+alpha*pR*tau_s*rrange(i));
    dDdt = (1-D(i))/tau_d - pR*D(i)*rrange(i);
    D(i+1) = D(i) + dt*dDdt;
end  

plot(s_r, rrange, '--')
ylabel('Rate (Hz)')
xlabel('Synaptic Input')

%%

rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.5;
D = 1;
pR = .5;
Wee = 35;
tau_s = .002;
dt = .0001;
tvec = 0:dt:20;
S_in = 0;
start = 10;
duration = .05;
stop = 10+duration;
r0 = -0.1;
rvec = zeros(size(tvec));
svec = zeros(size(tvec));

s_in = zeros(size(tvec));
x = 1.2;
tau_d = .250;

D = zeros(size(tvec));

for i = 1:length(tvec)-1
    time = tvec(i);
    if (time >= start && time <= stop)
        s_in(i) = .05;
    end
    alpha = alpha0*D(i);
    
    dDdt = (1-D(i))/tau_d - pR*D(i)*rvec(i);
    D(i+1) = D(i) + dt*dDdt;
    
    total_input = Wee*svec(i) + s_in(i);
    r_ss = r0 + rmax*total_input^x/(total_input^x + sigma^x);
    rvec(i+1) = rvec(i) + dt*(r_ss - rvec(i) )/tau_r;
    if (rvec(i+1) < 0)
        rvec(i+1) = 0;
    end
    dsdt = -svec(i)/tau_s + alpha*pR*rvec(i)*(1-svec(i));
    svec(i+1) = svec(i) + dt*dsdt;
end

subplot(1,3,2)
plot(tvec, rvec)
title('Part 3')
xlabel('Time')
ylabel('Rate (Hz)')

%% Part 3c

rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.5;
D = 1;
pR = .5;
Wee = 35;
tau_s = .002;
dt = .0001;
tvec = 0:dt:20;
S_in = 0;
start = 10;
duration = .05;
stop = 10+duration;
r0 = -0.1;
rvec = zeros(size(tvec));
rvec(1) = 9;
svec = zeros(size(tvec));
svec(1) = (alpha0*pR*rvec(1)*tau_s)/(1+alpha0*pR*rvec(1)*tau_s);

s_in = zeros(size(tvec));
x = 1.2;
tau_d = .250;

D = zeros(size(tvec));
D(1) = 1/(1+pR*rvec(1)*tau_d);

for i = 1:length(tvec)-1
    time = tvec(i);
    if (time >= start && time <= stop)
        s_in(i) = .05;
    end
    alpha = alpha0*D(i);
    
    dDdt = (1-D(i))/tau_d - pR*D(i)*rvec(i);
    D(i+1) = D(i) + dt*dDdt;
    
    total_input = Wee*svec(i) + s_in(i);
    r_ss = r0 + rmax*total_input^x/(total_input^x + sigma^x);
    rvec(i+1) = rvec(i) + dt*(r_ss - rvec(i) )/tau_r;
    if (rvec(i+1) < 0)
        rvec(i+1) = 0;
    end
    dsdt = -svec(i)/tau_s + alpha*pR*rvec(i)*(1-svec(i));
    svec(i+1) = svec(i) + dt*dsdt;
end

subplot(1,3,3)
title('Part 3')
plot(tvec, rvec)
title('Part 3')
xlabel('Time')
ylabel('Rate (Hz)')

%% Part 4
rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.25;
pR = 1;
Wee = 35;
tau_s = .002;
r0 = -0.01;
x = 1.2;
tau_d = .125;

rrange = 0:1:rmax;
srange = 0:.01:1;

D = zeros(size(rrange));

figure;
total_input = Wee*srange;
r_s = r0 + rmax*total_input.^x./(total_input.^x + sigma^x);

subplot(1,2,1)
plot(srange, r_s)
title('Part 4')

hold on
for i = 1:length(rrange)-1
    alpha = alpha0*D(i);
    s_r(i) = (alpha*pR*tau_s*rrange(i))/(1+alpha*pR*tau_s*rrange(i));
    dDdt = (1-D(i))/tau_d - pR*D(i)*rrange(i);
    D(i+1) = D(i) + dt*dDdt;
end  
    
plot(s_r, rrange,'--')
ylabel('Rate (Hz)')
xlabel('Synaptic Input')

%%

rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 0.25;
D = 1;
pR = 1;
Wee = 35;
tau_s = .002;
dt = .0001;
tvec = 0:dt:20;
S_in = 0;
start = 10;
duration = .05;
stop = 10+duration;
r0 = -0.1;
rvec = zeros(size(tvec));
svec = zeros(size(tvec));

s_in = zeros(size(tvec));
x = 1.2;
tau_d = .125;

D = zeros(size(tvec));

for i = 1:length(tvec)-1
    time = tvec(i);
    if (time >= start && time <= stop)
        s_in(i) = .05;
    end
    alpha = alpha0*D(i);
    
    dDdt = (1-D(i))/tau_d - pR*D(i)*rvec(i);
    D(i+1) = D(i) + dt*dDdt;
    
    total_input = Wee*svec(i) + s_in(i);
    r_ss = r0 + rmax*total_input^x/(total_input^x + sigma^x);
    rvec(i+1) = rvec(i) + dt*(r_ss - rvec(i) )/tau_r;
    if (rvec(i+1) < 0)
        rvec(i+1) = 0;
    end
    dsdt = -svec(i)/tau_s + alpha*pR*rvec(i)*(1-svec(i));
    svec(i+1) = svec(i) + dt*dsdt;
end

subplot(1,2,2)
plot(tvec, rvec)
title('Part 4')
xlabel('Time')
ylabel('Rate (Hz)')

%% Part 5
rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 1;
pR = 1;
Wee = 1.6;
tau_s = .050;
r0 = 0.5;
x = 1.2;
D = 1;

rrange = 0:1:rmax;
srange = 0:.01:1;

figure;
total_input = Wee*srange;
r_s = r0 + rmax*total_input.^x./(total_input.^x + sigma^x);

subplot(1,3,1)
plot(srange, r_s)
title('Part 5')

hold on
alpha = alpha0*D;
s_r = (alpha*pR*tau_s*rrange)./(1+alpha*pR*tau_s*rrange);

plot(s_r, rrange,'--')
ylabel('Rate (Hz)')
xlabel('Synaptic Input')

%%

rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 1;
D = 1;
pR = 1;
Wee = 1.6;
tau_s = .050;
dt = .0001;
tvec = 0:dt:20;
S_in = 0;
start = 10;
duration = 6;
stop = 10+duration;
r0 = 0.5;
rvec = zeros(size(tvec));
svec = zeros(size(tvec));
s_in = zeros(size(tvec));
x = 1.2;

for i = 1:length(tvec)-1
    time = tvec(i);
    if (time >= start && time <= stop)
        s_in(i) = .001;
    end
    
    alpha = alpha0*D;
    total_input = Wee*svec(i) + s_in(i);
    r_ss = r0 + rmax*total_input^x/(total_input^x + sigma^x);
    rvec(i+1) = rvec(i) + dt*(r_ss - rvec(i) )/tau_r;
    if (rvec(i+1) < 0)
        rvec(i+1) = 0;
    end
    dsdt = -svec(i)/tau_s + alpha*pR*rvec(i)*(1-svec(i));
    svec(i+1) = svec(i) + dt*dsdt;
end

subplot(1,3,2)
plot(tvec, rvec)
title('Part 5')
xlabel('Time')
ylabel('Rate (Hz)')

%%

rmax = 100;
sigma = 0.5;
tau_r = .010;
alpha0 = 1;
D = 1;
pR = 1;
Wee = 1.6;
tau_s = .050;
dt = .0001;
tvec = 0:dt:20;
S_in = 0;
start = 10;
duration = 6;
stop = 10+duration;
r0 = 0.5;
rvec = zeros(size(tvec));
svec = zeros(size(tvec));
s_in = zeros(size(tvec));
x = 1.2;
index = 0;

for i = 1:length(tvec)-1
    time = tvec(i);
    start_temp = start + .5*index;
    stop_temp = start + (.5*index) + .1;
    if (time <= stop)
        if (time >= start_temp)
            if (time == stop_temp)
                index = index + 1;  
                s_in(i) = 0;
            end
            s_in(i) = .005;
        end
    end
    
    alpha = alpha0*D;
    total_input = Wee*svec(i) + s_in(i);
    r_ss = r0 + rmax*total_input^x/(total_input^x + sigma^x);
    rvec(i+1) = rvec(i) + dt*(r_ss - rvec(i) )/tau_r;
    if (rvec(i+1) < 0)
        rvec(i+1) = 0;
    end
    dsdt = -svec(i)/tau_s + alpha*pR*rvec(i)*(1-svec(i));
    svec(i+1) = svec(i) + dt*dsdt;
end

subplot(1,3,3)
plot(tvec, rvec)
title('Part 5')
xlabel('Time')
ylabel('Rate (Hz)')