pattern_array = zeros(17,17,4);

for i = 1:size(pattern_array,1)
    for j = 1:size(pattern_array,1)
        if (i == j)
            pattern_array(i,j,1) = 1;
        elseif (i + j == 17)
            pattern_array(i,j,1) = 1;
        else
            pattern_array(i,j,1) = 0;
        end
    end
end

for i = 2:15
    pattern_array(2,i,2) = 1;
    pattern_array(i,2,2) = 1;
    pattern_array(i,16,2) = 1;
    pattern_array(16,i,2) = 1;
end

for i = 3:15
    for j = 3:15
      if (mod(j,2) ~= 0)
         pattern_array(i,j,3) = 1; 
      end
    end
end

for i = 2:16
    for j = 2:16
      if (mod(i,2) == 0)
         pattern_array(i,j,4) = 1; 
      end
    end
end

input_pattern = pattern_array(:,:,4);

figure
imagesc(input_pattern)

indices = rand(1,289)<0.1;
input_pattern(indices) = 1-input_pattern(indices);

figure;
imagesc(input_pattern)

nunits = size(input_pattern,1) * size(input_pattern,2);

%%

dt = 0.0001;
tvec = 0:dt:1;

tau_r = .010;
r_max = 50;
I_th = 10;
deltaI = 1;

cvec = zeros(length(tvec),nunits);

rvec = zeros(length(tvec),nunits);

ivec = zeros(length(tvec),nunits);

wvec = zeros(nunits,nunits,400);
wvec(:,:,1) = -0.3/nunits;

epsilon_plus = 0.1/nunits;
r_t = 25;
epsilon_minus = epsilon_plus/4;

w_max = 8/nunits;
w_min = -8/nunits;

for j = 1:400
    input_pattern = pattern_array(:,:,randi(4));
    input_pattern = input_pattern(:);

    indices = rand(1,289)<0.1;
    input_pattern(indices) = 1-input_pattern(indices);

    nunits = size(input_pattern,1) * size(input_pattern,2);
    
    w_ij = wvec(:,:,j);
    
    indices = find(input_pattern);
    
    for t = 2:length(tvec)
        if (tvec(t-1) > 0 && tvec(t-1) <.5)
            for index = 1:length(indices)
                cvec(t-1, indices(index)) = 50;
            end
        end
        
        I_app = cvec(t-1,:);
        r_i = rvec(t-1,:);
        I_i = ivec(t-1,:);

        dridt = (-r_i + (r_max./(1+exp(-(I_i - I_th)/deltaI)))) / tau_r;
        rvec(t, :) = r_i + dridt*dt;

        ivec(t,:) = sum(r_i*w_ij) + I_app;
    end
    dwdt =(epsilon_plus*double(r_i' > r_t)*double(r_i > r_t)*dt)...
        - (epsilon_minus*double(r_i' > r_t)*double(r_i > r_t)*dt);
    wvec(:,:,j+1) = wvec(:,:,j) + dwdt;
    if (wvec(:,:,j+1) > w_max)
        wvec(:,:,j+1) = w_max;
    elseif (wvec(j+1) < w_min)
        wvec(:,:,j+1) = w_min;
    end
end
