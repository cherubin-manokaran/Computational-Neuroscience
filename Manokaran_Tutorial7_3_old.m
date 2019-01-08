trials = 800;
dt = 0.001;
tvec = 0:dt:.5;

rivec = zeros(5,length(tvec));
rdvec = zeros(2,length(tvec));

P = [.05 .95; .25 .75; .5 .5; .75 .25; .95 .05]; 

civec = zeros(size(rivec));
cdvec = zeros(size(rdvec));

I_th = 50;
I_sigma = 5;
rmax = 100;
tau = .02;

Wi = zeros(5,800);
Wi(:,1) = .2;

Wdec1 = zeros(size(rdvec));
Wdec1(:,1) = .5;  

Wdec2 = zeros(size(rdvec));
Wdec2(:,1) = -.5;

sigma = 1;

sun = 1;
rain = 2;

E_correct = 0.5;
E_incorrect = -0.5;
epsilon = 0.4;

for i = 1:trials 
    cdrand1 = randn(size(tvec));
    cdrand2 = randn(size(tvec));
    weather = round(rand(1))+1;
    for t = 1:length(tvec)-1
        if (tvec(t) >= .1)
           civec(weather,t) = (rand(1)<P(1,weather)) * 50;
           civec(weather,t) = (rand(1)<P(2,weather)) * 50;
           civec(weather,t) = (rand(1)<P(3,weather)) * 50;
           civec(weather,t) = (rand(1)<P(4,weather)) * 50;
           civec(weather,t) = (rand(1)<P(5,weather)) * 50;
        end
        
        Wa1 = Wi(1,i);
        Wb1 = Wi(2,i);
        Wc1 = Wi(3,i);
        Wd1 = Wi(4,i);
        We1 = Wi(5,i);
        
        W11 = Wdec1(1,t);
        W21 = Wdec1(2,t);
        
        W12 = Wdec2(1,t);
        W22 = Wdec2(2,t);
        
        r_i = rivec(:,t);
        r_d = rdvec(:,t);
        
        if (r_d(sun,end) > 40)
            decision = sun;
        elseif (r_d(rain,end) > 40)
            decision = rain;
        end

        if (decision == weather)
            E = E_correct;
            for j = 1:5
                if (civec(j,t) > 0)
                    Wi(j,i+1) = Wi(j,i) + E*epsilon;
                elseif (civec(j,t) == 0)
                    Wi(j,i+1) = Wi(j,i) - E*epsilon;
                end
            end
        elseif (decision ~= weather)
            E = E_incorrect;
            for j = 1:5
                if (civec(j,t) > 0)
                    Wi(j,i+1) = Wi(j,i) + E*epsilon;
                elseif (civec(j,t) == 0)
                    Wi(j,i+1) = Wi(j,i) - E*epsilon;
                end
            end
        end
        
        I_app = civec(:,t);
        I_i = cdvec(:,t);
        
        dridt = (-r_i' + (rmax/((1+exp(I_th-I_app)/I_sigma))))/tau;
        rivec(:,t+1) = r_i' + dt*dridt;
        drddt = (-r_d' + (rmax/(1+exp(I_th-I_i)/I_sigma)))/tau;
        rdvec(:,t+1) = r_d' + dt*drddt;
        
        cdvec(1,1,t+1) = Wa1*rivec(1,t+1) + Wb1*rivec(2,t+1) + Wc1*rivec(3,t+1) + ...
            Wd1*rivec(4,t+1) + We1*rivec(5,t+1) + W11*rdvec(1,t+1) + W21*rdvec(2,t+1) + ...
            cdrand1(1,t);
        cdvec(2,1,t+1) = Wa1*rivec(1,t+1) + Wb1*rivec(2,t+1) + Wc1*rivec(3,t+1) + ...
            Wd1*rivec(4,t+1) + We1*rivec(5,t+1) + W21*rdvec(1,t+1) + W22*rdvec(2,t+1) + ...
            cdrand2(1,t)*sigma/sqrt(dt);
    end
end

rand(
